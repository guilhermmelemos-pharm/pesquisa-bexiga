import streamlit as st
from Bio import Entrez
import requests
import json
import re
from collections import Counter
import ast

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br"

# --- LISTA DE MODELOS ATUALIZADA (Baseada na sua permissão) ---
# Usei os mais potentes e rápidos da sua lista
MODELOS_ATIVOS = [
    "gemini-3-flash-preview",
    "gemini-2.5-flash",
    "gemini-2.0-flash",
    "gemini-flash-latest"
]

def montar_url(modelo, chave):
    # Removendo prefixo 'models/' se existir para a URL
    mod = modelo.replace("models/", "")
    return f"https://generativelanguage.googleapis.com/v1beta/models/{mod}:generateContent?key={chave}"

# --- 1. BLACKLIST RADICAL (Pré-IA) ---
BLACKLIST_RADICAL = {
    "URINARY", "BLADDER", "NEUROGENIC", "SWIM", "WITH", "DYSFUNCTION", "SPINAL", "LOWER", 
    "ULTRASOUND", "COMMENT", "RECURRENT", "EDITORIAL", "IMAGING", "POSTERIOR", "ACUTE",
    "RADIATION", "TISSUE", "EFFECT", "EFFECTS", "EXPRESSION", "INDUCED", "INFLAMMATORY"
}

# --- 2. FAXINEIRO IA (PROMPT EXPERIMENTAL) ---
def _faxina_ia(lista_suja):
    api_key = st.session_state.get('api_key_usuario', '').strip()
    if not api_key: return lista_suja[:40]

    prompt = f"""
    ROLE: PhD Senior Pharmacologist (Expert in Organ Bath & Molecular Biology).
    TASK: Review this list of terms.
    GOAL: Keep ONLY specific molecular targets (receptors, channels, enzymes) and experimental drugs.
    
    STRICT RULES:
    - KEEP: Receptors (TRPV1, M3), Channels (BKCa), Enzymes (ROCK, mTOR), Drugs (Cyclophosphamide, Mirabegron), Signaling (ATP, cAMP).
    - DELETE: Anatomy (bladder, spinal), Clinical concepts (dysfunction, imaging), Filler (with, effects, editorial).
    - FOCUS: If it's not a molecule I can pipette or quantify in a Western Blot, DELETE it.
    
    INPUT: {", ".join(lista_suja)}
    
    OUTPUT: Return strictly a Python list of strings.
    """
    
    payload = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.1}}
    
    for m in MODELOS_ATIVOS:
        try:
            url = montar_url(m, api_key)
            resp = requests.post(url, json=payload, timeout=15)
            if resp.status_code == 200:
                raw_text = resp.json()['candidates'][0]['content']['parts'][0]['text']
                clean_list = re.sub(r'```[a-z]*', '', raw_text).replace('```', '').strip()
                res = ast.literal_eval(clean_list)
                return res
        except: continue
    return lista_suja[:40]

# --- 3. ANÁLISE DE ABSTRACT (FOCO EM BANCADA) ---
def analisar_abstract_com_ia(titulo, resumo_texto, api_key):
    if not api_key: return "Chave API não configurada."
    
    prompt = f"""
    Resuma como Doutor em Farmacologia (máx 20 palavras).
    Foque no Alvo (Receptor/Via) e na Resposta Mecânica (Contração/Relaxamento).
    Ignore dados epidemiológicos ou clínicos.
    TÍTULO: {titulo}
    TEXTO: {resumo_texto[:800]}
    """
    
    for m in MODELOS_ATIVOS:
        try:
            url = montar_url(m, api_key)
            resp = requests.post(url, json={"contents": [{"parts": [{"text": prompt}]}]}, timeout=10)
            if resp.status_code == 200:
                return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
        except: continue
    return "Erro de conexão com modelos 2.5/3.0."

# --- 4. BUSCA PUBMED ---
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    query = f"({termo_base}) AND (2020:2026[Date - Publication]) NOT Review[pt]"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1000, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        lines = handle.read().splitlines(); handle.close()
        
        candidatos = []
        for line in lines:
            if line.startswith("TI  - ") or line.startswith("OT  - "):
                found = re.findall(r'\b(?:[A-Z]{2,}[A-Z0-9-]*|[A-Z][a-z]{3,})\b', line)
                for f in found:
                    if f.upper() not in BLACKLIST_RADICAL and len(f) > 2:
                        candidatos.append(f)
        
        top_terms = [t for t, count in Counter(candidatos).most_common(120)]
        if usar_ia: return _faxina_ia(top_terms)
        return top_terms[:60]
    except: return []

# --- 5. BUSCA RESUMOS PARA O FRONT ---
def buscar_resumos_detalhados(termo, orgao, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    query = f"({termo}) AND ({orgao}) AND ({ano_ini}:{ano_fim}[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        artigos = []
        for raw in dados.split("\n\nPMID-"):
            tit, abstract, pmid = "", "", ""
            for line in raw.split("\n"):
                if line.startswith("TI  - "): tit = line[6:].strip()
                if line.startswith("AB  - "): abstract = line[6:500].strip()
                if line.startswith("PMID- "): pmid = line[6:].strip()
            if tit: artigos.append({"Title": tit, "Info_IA": abstract, "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return artigos
    except: return []

# --- 6. CONSULTA COUNT (FISHER) ---
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    query = f"({termo})" + (f" AND ({contexto})" if contexto else "")
    query += f" AND ({ano_ini}:{ano_fim}[Date - Publication]) AND (NOT Review[pt])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle); handle.close()
        return int(record["Count"])
    except: return 0

def buscar_todas_noticias(lang='pt'):
    return [] # Radar simplificado para focar no erro
