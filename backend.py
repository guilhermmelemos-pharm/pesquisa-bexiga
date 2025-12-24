import streamlit as st
from Bio import Entrez
import requests
import json
import re
from collections import Counter
import ast
import math
from tenacity import retry, stop_after_attempt, wait_exponential

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br"

# MODELOS DO SEU PAID TIER
MODELOS_ATIVOS = ["gemini-2.0-flash", "gemini-2.0-flash-exp", "gemini-flash-latest"]

# --- 1. MOTOR DE MÉTRICAS ---
def calcular_metricas_originais(freq, total_docs, n_alvo_total):
    lambda_score = (freq / total_docs) * 100 if total_docs > 0 else 0
    p_value = math.exp(-freq/5) if freq > 0 else 1.0
    blue_ocean = max(0, 100 - (n_alvo_total / 10)) if n_alvo_total > 0 else 100.0
    status = "Saturado" if blue_ocean < 25 else "Blue Ocean" if blue_ocean > 75 else "Competitivo"
    return lambda_score, p_value, blue_ocean, status

# --- 2. FAXINEIRO IA (ORDEM DE EXTERMÍNIO E FORMATAÇÃO) ---
def _faxina_ia_elite(lista_bruta):
    api_key = st.session_state.get('api_key_usuario', '').strip()
    if not api_key: return []
    
    prompt = f"""
    ROLE: Senior PhD in Pharmacology (Molecular Biology focus).
    INPUT: {lista_bruta[:120]}
    
    TASK: Extract exactly 20 specific molecular targets and drugs.
    STRICT FORMAT: Return a Python list of dictionaries:
    [{"Alvo": "TRPV4", "Farmaco": "GSK1016790A", "Acao": "Agonista"}, ...]
    
    RULES:
    1. DELETE: Toads, Frogs, Dogs, Sodium, Water, Muscle, Cell, Urinary, Metabolism.
    2. DELETE: General anatomy (Bladder, Kidney) and clinical terms (Botox, Surgery).
    3. KEEP: Receptors, Channels, Enzymes, and specific experimental compounds.
    """
    
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.1}}
    
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, json=data, timeout=20)
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                match = re.search(r'\[.*\]', texto, re.DOTALL)
                if match: return ast.literal_eval(match.group())
        except: continue
    return []

# --- 3. MINERAÇÃO E BUSCA (ESTATÍSTICA + SNIPER) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    
    # BUSCA RESTRITA: 2020-2026 para matar estudos de 100 anos atrás
    query = f"({termo_base} AND (Pharmacology OR Molecular OR Signaling)) AND (2020:2026[Date])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1000, sort="relevance")
        record = Entrez.read(handle); handle.close()
        total_docs = int(record["Count"])
        
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        # Blacklist de Extermínio (Fundida com a sua versão anterior)
        blacklist = {
            'TOAD', 'FROG', 'DOG', 'RABBIT', 'SODIUM', 'WATER', 'TRANSPORT', 'CELL', 'MUSCLE', 
            'BLADDER', 'URINARY', 'EPITHELIUM', 'STUDY', 'ROLE', 'SYSTEM', 'ACTION', 'EFFECT',
            'DNA', 'RNA', 'PROTEIN', 'GENE', 'PATIENT', 'HUMAN', 'CLINICAL', 'TREATMENT'
        }
        
        candidatos_raw = []
        for artigo in full_data.split("\n\nPMID-"):
            texto_focado = ""
            for line in artigo.split("\n"):
                if line.startswith("TI  - ") or line.startswith("KW  - ") or line.startswith("OT  - "):
                    texto_focado += line[6:].strip() + " "
            
            siglas = re.findall(r'\b[A-Z0-9-]{3,15}\b', texto_focado)
            for s in siglas:
                if s.upper() not in blacklist and not s.isdigit():
                    candidatos_raw.append(s.upper())

        contagem = Counter(candidatos_raw)
        top_nomes = [t for t, f in contagem.most_common(150)]
        
        if usar_ia:
            entidades = _faxina_ia_elite(top_nomes)
            res_finais = []
            for ent in entidades:
                freq = contagem.get(ent['Alvo'].upper(), 1)
                l_score, p_val, b_ocean, status = calcular_metricas_originais(freq, total_docs, freq*10)
                res_finais.append({
                    "Alvo": ent['Alvo'], "Farmaco": ent['Farmaco'], "Acao": ent['Acao'],
                    "Lambda": round(l_score, 2), "Blue Ocean": round(b_ocean, 1), "Status": status
                })
            return res_finais
        return []
    except: return []

# --- 4. RADAR E APOIO ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_todas_noticias(lang='pt'):
    try:
        query = "(bladder pharmacology OR ion channels signaling) AND (2024:2026[Date])"
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
                news.append({"titulo": tit, "fonte": journal[:40], "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return news
    except: return []

@st.cache_data(ttl=86400, show_spinner=False)
def buscar_resumos_detalhados(termo, orgao, email, ano_ini, ano_fim):
    query = f"({termo}) AND ({orgao}) AND (2020:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5, sort="relevance")
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        artigos = []
        for raw in dados.split("\n\nPMID-"):
            tit, pmid, abstract = "", "", ""
            for line in raw.split("\n"):
                if line.startswith("TI  - "): tit = line[6:].strip()
                if line.startswith("AB  - "): abstract = line[6:400].strip()
                if line.startswith("PMID- "): pmid = line[6:].strip()
            if tit:
                artigos.append({"Title": tit, "Info_IA": abstract, "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return artigos
    except: return []

def analisar_abstract_com_ia(titulo, dados, api_key, lang='pt'):
    # Retorna o formato Alvo | Fármaco | Ação para os detalhes
    prompt = f"Based on {titulo} {dados}, extract: Target | Drug | Action. Concise."
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}]}
    url = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.0-flash:generateContent?key={api_key}"
    try:
        resp = requests.post(url, headers=headers, json=data, timeout=10)
        return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
    except: return "N/A | N/A | N/A"

def consultar_pubmed_count(termo, contexto, email, ano_ini=2020, ano_fim=2026):
    if email: Entrez.email = email
    query = f"({termo}) AND ({ano_ini}:{ano_fim}[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle); handle.close()
        return int(record["Count"])
    except: return 0
