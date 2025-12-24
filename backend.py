import streamlit as st
from Bio import Entrez
import requests
import json
from tenacity import retry, stop_after_attempt, wait_exponential
import re
from collections import Counter
import time
import ast

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br"

MAPA_SINONIMOS = {
    "BLADDER": "Bladder OR Urothelium OR Detrusor OR Vesical OR Urethra OR Micturition OR LUTS OR Cystitis OR Overactive Bladder OR OAB OR Urinary Tract",
    "PAIN": "Pain OR Nociception OR Analgesia OR Neuropathic OR TRP Channels",
    "INFLAMMATION": "Inflammation OR Cytokine OR Macrophage OR Sepsis OR Inflammasome OR NF-kB"
}

# Foco nos modelos estáveis do seu tier pago
MODELOS_ATIVOS = ["gemini-2.0-flash", "gemini-2.0-flash-exp"]

def montar_url_limpa(modelo, chave):
    return f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={chave.strip()}"

# --- 1. FAXINEIRO IA (COM TRATAMENTO DE ERRO 429) ---
def _faxina_ia(lista_suja):
    api_key = st.session_state.get('api_key_usuario', '').strip()
    if not api_key: return lista_suja[:30] 

    prompt = f"""
    ROLE: Senior Pharmacologist (Molecular focus).
    INPUT: {", ".join(lista_suja)}
    TASK: Keep ONLY molecular targets (Channels, Receptors) and experimental drugs.
    DELETE ALL CLINICAL: (OAB, LUTS, MRI, TURBT, BCG, Botox, Patient, Surgery, Treatment).
    OUTPUT: Return strictly a Python list [].
    """
    
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.0}}

    for m in MODELOS_ATIVOS:
        try:
            url_final = montar_url_limpa(m, api_key)
            resp = requests.post(url_final, headers=headers, json=data, timeout=20)
            
            if resp.status_code == 429:
                time.sleep(2) # Espera estratégica para evitar bloqueio
                continue
                
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                texto = texto.replace("```python", "").replace("```", "").strip()
                return ast.literal_eval(texto)
        except: continue
    return lista_suja[:30]

# --- 2. ANÁLISE SNIPER: ALVO | FÁRMACO | AÇÃO (CORRIGE 429) ---
def analisar_abstract_com_ia(titulo, dados_curtos, api_key, lang='pt'):
    key = api_key.strip()
    if not key: return "⚠️ Key Ausente"
    idioma = "Português" if lang == 'pt' else "Inglês"
    
    prompt = f"Extract strictly: TARGET | DRUG | ACTION. Based on: {titulo}. Language: {idioma}."
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.0}}

    for m in MODELOS_ATIVOS:
        try:
            url_final = montar_url_limpa(m, key)
            resp = requests.post(url_final, headers=headers, json=data, timeout=12)
            
            if resp.status_code == 429:
                time.sleep(1) # Aguarda cota respirar
                continue
                
            if resp.status_code == 200:
                return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
        except: continue
    return "⚠️ Cota de IA Excedida (Aguarde 30s)"

# --- 3. MINERAÇÃO (FILTRO REFORÇADO +90%) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    termo_upper = termo_base.upper().strip()
    query_string = MAPA_SINONIMOS.get(termo_upper, f"{termo_base}[Title/Abstract]")
    
    # TRAVA TEMPORAL INTERNA 2020+: Elimina 90% do ruído histórico
    final_query = f"({query_string}) AND (2020:2026[Date - Publication]) AND (NOT Review[pt])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=final_query, retmax=1500, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        artigos_raw = full_data.split("\n\nPMID-")

        # BLACKLIST DE ELITE (Foco em limpar a Bancada)
        blacklist = {
            "OAB", "LUTS", "MRI", "ICS", "PTNS", "BPS", "LUT", "VI-RADS", "UTI", "ICI-RS", "BPH", 
            "EMG", "LUTD", "PET", "TURBT", "SUI", "COVID-19", "SNM", "WHO", "BCG", "NOTNLM",
            "DNA", "RNA", "CELL", "MUSCLE", "BLADDER", "URINARY", "EPITHELIUM", "ROLE", "STUDY", "PATIENT",
            "THE", "AND", "FOR", "WITH", "CASE", "MANAGEMENT", "CLINICAL", "TRIAL"
        }

        candidatos = []
        for artigo in artigos_raw:
            texto_sniper = ""
            for line in artigo.split("\n"):
                if line.startswith("TI  - ") or line.startswith("KW  - ") or line.startswith("OT  - "):
                    texto_sniper += line[6:].strip() + " "
            
            # REGEX SNIPER: Só pega siglas técnicas (ex: TRPV4, P2X3, ROCK1)
            # Elimina palavras comuns que começam com maiúscula
            encontrados = re.findall(r'\b[A-Z0-9-]{3,15}\b', texto_sniper)
            for s in encontrados:
                if s.upper() not in blacklist and not s.isdigit():
                    candidatos.append(s.upper())

        contagem = Counter(candidatos)
        top_nomes = [t for t, f in contagem.most_common(120)]
        
        if usar_ia:
            return _faxina_ia(top_nomes)
        return top_nomes[:40]
    except: return []

# --- 4. FUNÇÕES DE APOIO ---
@st.cache_data(ttl=86400, show_spinner=False)
def consultar_pubmed_count(termo, contexto, email, ano_ini=2015, ano_fim=2026):
    if email: Entrez.email = email
    query = f"({termo}) AND (2015:2026[Date])" # Força ciência moderna
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
    query = f"({termo}) AND ({q_orgao}) AND (2018:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        artigos = []
        for raw in dados.split("\n\nPMID-"):
            tit, pmid, abstract = "", "", ""
            for line in raw.split("\n"):
                if line.startswith("TI  - "): tit = line[6:].strip()
                if line.startswith("AB  - "): abstract = line[6:500].strip()
                if line.startswith("PMID- "): pmid = line[6:].strip()
            if tit:
                artigos.append({"Title": tit, "Info_IA": abstract, "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return artigos
    except: return []

def buscar_todas_noticias(lang='pt'):
    try:
        query = "(molecular pharmacology OR ion channels) AND (2024:2026[Date])"
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
