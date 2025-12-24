import streamlit as st
from Bio import Entrez
import requests
import json
import re
from collections import Counter
import ast

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br"

# --- MODELOS DA SUA LISTA (FOCO EM 2.5 E 3.0) ---
MODELOS_ATIVOS = [
    "gemini-3-flash-preview",
    "gemini-2.5-flash",
    "gemini-2.0-flash"
]

# --- 1. BLACKLIST DE PRE-FILTRAGEM (Python mata o lixo clínico antes) ---
BLACKLIST_RADICAL = {
    "VI-RADS", "RADS", "RCT", "TURP", "TURBT", "BCG", "III", "3-D", "FISH",
    "GUERIN", "CALMETTE", "DISEASE", "SURGERY", "CLINICAL", "PATIENTS",
    "THE", "CASE", "UPDATE", "MANAGEMENT", "PRIMARY", "MRI-", "RNA-"
}

# --- 2. O PROMPT DE FILTRAGEM RIGOROSA ---
def _faxina_ia(lista_suja):
    api_key = st.session_state.get('api_key_usuario', '').strip()
    if not api_key: 
        return [t for t in lista_suja if t.upper() not in BLACKLIST_RADICAL][:30]

    # PROMPT COM FOCO EM FARMACOLOGIA EXPERIMENTAL
    prompt = f"""
    ROLE: PhD Senior Researcher in Pharmacophysiology.
    CONTEXT: You are analyzing potential targets for isolated organ bath experiments and molecular signaling.
    
    INPUT LIST: {", ".join(lista_suja)}
    
    TASK: Strictly filter this list and return ONLY relevant molecular entities.
    
    ✅ KEEP (STRICT):
    - Receptors & Channels: TRPV4, TRPA1, 5-HT, P2X, Muscarinic receptors.
    - Signaling Proteins/Enzymes: ERK, STING, NLRP3, TGF-beta, NF-kB, BDNF.
    - Molecules/Drugs: ATP, Cyclophosphamide, Vitamin (if specific), Resveratrol.
    - Genetic Markers/EMT: EMT markers (Vimentin, E-cadherin), Genes.
    
    ❌ DELETE (STRICT):
    - Clinical/Surgical terms: VI-RADS, TURP, RCT, Surgery, Patients.
    - Anatomy/Disease names: Bladder, Cancer, Cystitis, SCI (unless target-related).
    - Vague terms: Disease, Case, Update, 3-D, III, Guerin, FISH.
    
    OUTPUT: Return strictly a Python list of strings. If the list is empty, return [].
    """
    
    payload = {
        "contents": [{"parts": [{"text": prompt}]}],
        "generationConfig": {"temperature": 0.0} # Precisão máxima
    }
    
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, json=payload, timeout=15)
            
            if resp.status_code == 200:
                raw_text = resp.json()['candidates'][0]['content']['parts'][0]['text']
                # Limpeza de resíduos de markdown (```python ... ```)
                clean_text = re.sub(r'```[a-zA-Z]*', '', raw_text).replace('```', '').strip()
                res = ast.literal_eval(clean_text)
                if isinstance(res, list): return res
        except: continue
    
    return [t for t in lista_suja if t.upper() not in BLACKLIST_RADICAL][:30]

# --- 3. MINERAÇÃO PUBMED (AJUSTADA) ---
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    # Busca focada em mecanismos e farmacologia
    query = f"({termo_base}) AND (2022:2026[Date - Publication]) NOT (Review[pt] OR Case Reports[pt])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1500, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        lines = handle.read().splitlines(); handle.close()
        
        candidatos = []
        for line in lines:
            if line.startswith("TI  - ") or line.startswith("OT  - "):
                # Pega siglas e termos técnicos Capitalizados
                found = re.findall(r'\b(?:[A-Z0-9-]{3,12}|[A-Z][a-z]{3,}(?:in|one|ol|ide|ase|ant|receptor|channel))\b', line)
                for f in found:
                    f_up = f.upper()
                    if f_up not in BLACKLIST_RADICAL and len(f) > 2:
                        candidatos.append(f)
        
        top_terms = [t for t, count in Counter(candidatos).most_common(120)]
        
        if usar_ia:
            return _faxina_ia(top_terms)
        return top_terms[:30]
    except: return []

# --- 4. ANALISAR ABSTRACT ---
def analisar_abstract_com_ia(titulo, resumo_texto, api_key):
    prompt = f"""
    Resuma como Doutor em Farmacologia (máx 20 palavras):
    Foque no Alvo Molecular e na sinalização intracelular (ex: ERK, NF-kB).
    TÍTULO: {titulo}
    TEXTO: {resumo_texto[:800]}
    """
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, json={"contents": [{"parts": [{"text": prompt}]}]}, timeout=10)
            if resp.status_code == 200:
                return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
        except: continue
    return "Erro de resposta da IA."

# --- FUNÇÕES SUPORTE ---
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    query = f"({termo}) AND ({contexto}) AND ({ano_ini}:{ano_fim}[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        return int(Entrez.read(handle)["Count"])
    except: return 0

def buscar_resumos_detalhados(termo, orgao, email, ano_ini, ano_fim):
    query = f"({termo}) AND ({orgao}) AND ({ano_ini}:{ano_fim}[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5)
        ids = Entrez.read(handle)["IdList"]
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
        dados = handle.read().split("\n\nPMID-")
        artigos = []
        for art in dados:
            tit = re.search(r"TI  - (.*)", art)
            ab = re.search(r"AB  - (.*)", art)
            if tit: artigos.append({"Title": tit.group(1), "Info_IA": ab.group(1) if ab else "", "Link": "https://pubmed.ncbi.nlm.nih.gov/"})
        return artigos
    except: return []

def buscar_todas_noticias(lang='pt'): return []
