import streamlit as st
from Bio import Entrez
import requests
import json
from tenacity import retry, stop_after_attempt, wait_exponential
import re
from collections import Counter
import ast
import math

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br"

MAPA_SINONIMOS = {
    "BLADDER": "Bladder OR Urothelium OR Detrusor OR Vesical OR Urethra OR Micturition OR LUTS OR Cystitis OR Overactive Bladder OR OAB OR Urinary Tract",
    "PAIN": "Pain OR Nociception OR Analgesia OR Neuropathic OR TRP Channels",
    "INFLAMMATION": "Inflammation OR Cytokine OR Macrophage OR Sepsis OR Inflammasome OR NF-kB"
}

MODELOS_ATIVOS = ["gemini-2.0-flash", "gemini-2.0-flash-exp", "gemini-flash-latest"]

def montar_url_limpa(modelo, chave):
    return f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={chave.strip()}"

# --- 1. MOTOR DE MÉTRICAS ---
def calcular_metricas_originais(freq, total_docs, n_alvo_total):
    lambda_score = (freq / total_docs) * 100 if total_docs > 0 else 0
    p_value = math.exp(-freq/5) if freq > 0 else 1.0
    blue_ocean = max(0, 100 - (n_alvo_total / 10)) if n_alvo_total > 0 else 100.0
    status = "Saturado" if blue_ocean < 25 else "Blue Ocean" if blue_ocean > 75 else "Competitivo"
    return lambda_score, p_value, blue_ocean, status

# --- 2. SNIPER IA: ALVO | FÁRMACO | AÇÃO ---
def analisar_abstract_com_ia(titulo, dados, api_key, lang='pt'):
    key = api_key.strip()
    if not key: return "Erro | Key Ausente"
    
    prompt = f"""
    ROLE: Senior PhD Pharmacologist.
    INPUT: Title: {titulo} | Data: {dados}
    TASK: Extract strictly: TARGET | DRUG | ACTION. 
    RULES: 
    - Output ONLY the 3 fields separated by '|'.
    - DELETE clinical terms (OAB, LUTS, MRI, UTI).
    - FOCUS on molecular targets (TRP, Piezo, ROCK).
    """
    
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.0}}
    
    for m in MODELOS_ATIVOS:
        try:
            url = montar_url_limpa(m, key)
            resp = requests.post(url, headers=headers, json=data, timeout=12)
            if resp.status_code == 200:
                res = resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
                return res if "|" in res else f"{res} | Endógeno | Mecanismo"
        except: continue
    return "N/A | N/A | N/A"

# --- 3. FAXINEIRO IA (REFRENA A LISTA DE TAGS) ---
def _faxina_ia(lista_suja):
    api_key = st.session_state.get('api_key_usuario', '').strip()
    if not api_key: return lista_suja[:30]
    
    prompt = f"""
    ROLE: Senior Scientist in Pharmacology.
    INPUT: {", ".join(lista_suja)}
    TASK: Keep ONLY specific pharmacological targets (Channels, Receptors, Enzymes) and specific experimental drugs.
    DELETE ALL CLINICAL/DIAGNOSTIC JARGON: (OAB, LUTS, MRI, ICS, PTNS, BPS, VI-RADS, UTI, BPH, TURBT, SUI, SNM, COVID-19, BCG, WHO).
    DELETE: Animals, anatomy, generic biology.
    OUTPUT: Strictly a Python list [].
    """
    
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.0}}
    
    for m in MODELOS_ATIVOS:
        try:
            url = montar_url_limpa(m, api_key)
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=20)
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                match = re.search(r'\[.*\]', texto, re.DOTALL)
                if match: return ast.literal_eval(match.group())
        except: continue
    return lista_suja[:30]

# --- 4. BUSCA E CONTAGEM (FILTRO TEMPORAL FIXO) ---
@st.cache_data(ttl=86400, show_spinner=False)
def consultar_pubmed_count(termo, contexto, email, ano_ini=2015, ano_fim=2026):
    if email: Entrez.email = email
    query = f"({termo}) AND (2015:2026[Date - Publication])"
    if contexto:
        q_contexto = MAPA_SINONIMOS.get(contexto.upper(), contexto)
        query += f" AND ({q_contexto})"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle); handle.close()
        return int(record["Count"])
    except: return 0

@st.cache_data(ttl=86400, show_spinner=False)
def buscar_resumos_detalhados(termo, orgao, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    q_orgao = MAPA_SINONIMOS.get(orgao.upper(), orgao)
    query = f"({termo}) AND ({q_orgao}) AND (2018:2026[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5, sort="relevance")
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        artigos = []
        for raw in dados.split("\n\nPMID-"):
            tit, pmid, keywords, abstract = "", "", "", ""
            for line in raw.split("\n"):
                if line.startswith("TI  - "): tit = line[6:].strip()
                if line.startswith("AB  - "): abstract = line[6:500].strip()
                if line.startswith("OT  - ") or line.startswith("KW  - "): keywords += line[6:].strip() + ", "
            if tit:
                artigos.append({"Title": tit, "Info_IA": f"{keywords} {abstract}", "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return artigos
    except: return []

# --- 5. MINERAÇÃO (CORE COM COUNTER E BLACKLIST REFORÇADA) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    termo_upper = termo_base.upper().strip()
    query_string = MAPA_SINONIMOS.get(termo_upper, f"{termo_base}[Title/Abstract]")
    final_query = f"({query_string}) AND (2020:2026[Date - Publication]) AND (NOT Review[pt])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=final_query, retmax=1000, sort="relevance")
        record = Entrez.read(handle); handle.close()
        
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        # BLACKLIST DE EXTERMÍNIO CLÍNICO (Aumentada conforme sua reclamação)
        blacklist = {
            'OAB', 'LUTS', 'MRI', 'ICS', 'PTNS', 'BPS', 'LUT', 'VI-RADS', 'UTI', 'ICI-RS', 'BPH', 
            'EMG', 'LUTD', 'PET', 'TURBT', 'SUI', 'COVID-19', 'SNM', 'WHO', 'BCG', 'NOTNLM',
            'DNA', 'RNA', 'CELL', 'MUSCLE', 'BLADDER', 'URINARY', 'EPITHELIUM', 'STUDY', 'PATIENT'
        }

        candidatos = []
        for artigo in full_data.split("\n\nPMID-"):
            texto_sniper = ""
            for line in artigo.split("\n"):
                if line.strip().startswith("TI") or line.strip().startswith("KW") or line.strip().startswith("OT"):
                    texto_sniper += line[6:].strip() + " "
            
            siglas = re.findall(r'\b[A-Z0-9-]{3,15}\b', texto_sniper)
            for s in siglas:
                if s.upper() not in blacklist and not s.isdigit():
                    candidatos.append(s.upper())

        contagem = Counter(candidatos)
        top_nomes = [t for t, f in contagem.most_common(120)]
        
        if usar_ia:
            return _faxina_ia(top_nomes)
        return top_nomes[:40]
    except: return []

# --- 6. RADAR DE NOTÍCIAS 2025 ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_todas_noticias(lang='pt'):
    try:
        query = "(bladder pharmacology signaling OR urothelium ion channels) AND (2024:2026[Date])"
        handle = Entrez.esearch(db="pubmed", term=query, retmax=6, sort="pub_date")
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        news = []
        for art in dados.split("\n\nPMID-"):
            tit, journal, pmid = "", "", ""
            for line in art.split("\n"):
                if line.startswith("TI  - "): tit = line[6:].strip()
                if line.startswith("JT  - "): journal = line[3:].strip()
                if line.strip().isdigit() and not pmid: pmid = line.strip()
            if tit:
                news.append({"titulo": tit, "fonte": journal[:35], "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return news
    except: return []
