# backend.py
import streamlit as st
from Bio import Entrez
import requests
import re
import time
import ast
import json
import math
from collections import Counter
from difflib import SequenceMatcher

# ================= CONFIG =================
Entrez.email = "pesquisador_guest@unifesp.br"
MODELOS_ATIVOS = ["gemini-2.0-flash", "gemini-2.0-flash-exp"]

MAPA_SINONIMOS_BASE = {
    "BLADDER": "(Bladder OR Vesical OR Detrusor OR Urothelium)",
    "PAIN": "(Pain OR Nociception OR Analgesia)",
    "INFLAMMATION": "(Inflammation OR Cytokines OR Inflammasome)"
}

# ================= GEMINI CORE =================

def call_gemini(prompt, api_key, temperature=0.0):
    if not api_key: return ""
    headers = {"Content-Type": "application/json"}
    payload = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": temperature}}
    
    for modelo in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={api_key.strip()}"
            resp = requests.post(url, headers=headers, json=payload, timeout=25)
            if resp.status_code == 200:
                return resp.json()["candidates"][0]["content"]["parts"][0]["text"]
            elif resp.status_code == 429:
                time.sleep(1)
        except: continue
    return ""

def parse_structure(text, expected=list):
    if not text: return expected()
    text = re.sub(r"```.*?```", "", text, flags=re.DOTALL).strip()
    try:
        match = re.search(r"\[.*\]" if expected == list else r"\{.*\}", text, re.DOTALL)
        if match:
            obj = ast.literal_eval(match.group(0))
            if isinstance(obj, expected): return obj
    except: pass
    return expected()

# ================= FILTRO DE ELITE (MATA BCG, MRI, LUTS) =================

def ner_extraction_batch(artigos, api_key):
    texto_input = "\n".join([f"- {a['texto']}" for a in artigos[:25]])
    prompt = f"""
    ROLE: Senior Molecular Pharmacologist.
    TASK: Extract Molecular Targets and Experimental Drugs.
    
    STRICT EXCLUSION (DELETE):
    - Clinical terms: BCG, MRI, NMIBC, FDG-PET, ICER, TURBT, Cystoscopy.
    - Symptoms/Scores: LUTS, OAB, IPSS, QoL.
    - General: Bladder, Patient, Study, Treatment.
    
    KEEP ONLY: Receptors, Channels, Enzymes, Signaling Proteins (e.g., HSP90, LRP5, NF-kB, TRPV4).
    
    TEXT:
    {texto_input}
    
    OUTPUT: Python list of strings only.
    """
    raw = call_gemini(prompt, api_key, temperature=0.1)
    return parse_structure(raw, list)

# ================= ANALISAR ABSTRACT (FIX ERRO LINHA 260) =================

def analisar_abstract_com_ia(titulo, dados_curtos, api_key, lang='pt'):
    """
    FIX: Agora aceita os 4 argumentos exigidos pelo app_doutorado.py
    """
    if not api_key: return "N/A | N/A | N/A"
    idioma = "Português" if lang == 'pt' else "Inglês"
    
    prompt = f"""
    Analyze: {titulo}. Context: {dados_curtos}. 
    Extract: TARGET | DRUG | MECHANISM. 
    Strictly one line. Language: {idioma}.
    """
    res = call_gemini(prompt, api_key)
    return res.strip() if res else "N/A | N/A | N/A"

# ================= MÉTRICAS E PIPELINE =================

def calcular_metricas_originais(freq, total_docs, n_alvo_total):
    lambda_score = (freq / total_docs) * 100 if total_docs > 0 else 0
    # Fix p-value: escala logarítmica para evitar zeros infinitos
    p_val = math.exp(-freq/3.5) if freq > 0 else 1.0
    blue_ocean = max(10, 100 - (n_alvo_total * 3))
    return lambda_score, round(p_val, 4), round(blue_ocean, 2), "Competitivo"

@st.cache_data(ttl=3600)
def minerar_pubmed(termo_base, email):
    api_key = st.session_state.get("api_key_usuario", "").strip()
    Entrez.email = email
    query = f"({termo_base}) AND (2020:2026[Date])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=150)
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return {}

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        raw_medline = handle.read(); handle.close()

        artigos = []
        for bloco in raw_medline.split("\n\nPMID-"):
            titulo, texto = "", ""
            for line in bloco.split("\n"):
                if line.startswith("TI  - "): titulo = line[6:].strip()
                if line.startswith(("TI  - ", "KW  - ", "OT  - ")): texto += line[6:].strip() + " "
            if texto: artigos.append({"titulo": titulo, "texto": texto})

        entidades = ner_extraction_batch(artigos, api_key) if api_key else []
        counts = Counter(entidades)
        
        # Filtro de recorrência inteligente
        recorrentes = [e for e, c in counts.items() if c >= 2] or [e for e, _ in counts.most_common(15)]
        
        def normalize(entities):
            canon = []
            for e in entities:
                if not any(SequenceMatcher(None, e.lower(), c.lower()).ratio() > 0.85 for c in canon):
                    canon.append(e)
            return canon

        return {
            "termos_indicados": normalize(recorrentes),
            "counts": counts,
            "total_docs": len(artigos),
            "artigos_originais": artigos
        }
    except: return {}

# ================= FUNÇÕES DE APOIO (FRONTEND) =================

def buscar_alvos_emergentes_pubmed(alvo, email, usar_ia=True):
    res = minerar_pubmed(alvo, email)
    return res.get("termos_indicados", [])

@st.cache_data(ttl=3600)
def buscar_resumos_detalhados(termo, orgao, email, ano_ini, ano_fim):
    Entrez.email = email
    query = f"({termo}) AND ({orgao}) AND ({ano_ini}:{ano_fim}[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5)
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        raw = handle.read(); handle.close()
        artigos = []
        for bloco in raw.split("\n\nPMID-"):
            tit, pmid, abstract = "", "", ""
            for line in bloco.split("\n"):
                if line.startswith("TI  - "): tit = line[6:].strip()
                if line.startswith("AB  - "): abstract = line[6:500].strip()
                if line.startswith("PMID- "): pmid = line[6:].strip()
            if tit: artigos.append({"Title": tit, "Info_IA": abstract, "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return artigos
    except: return []

def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    Entrez.email = email
    query = f"({termo}) AND ({contexto}) AND ({ano_ini}:{ano_fim}[Date])"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    record = Entrez.read(handle); handle.close()
    return int(record["Count"])

def buscar_todas_noticias(email):
    Entrez.email = email
    query = "(molecular pharmacology OR bladder) AND (2024:2026[Date])"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=5, sort="pub_date")
    record = Entrez.read(handle); handle.close()
    handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
    raw = handle.read(); handle.close()
    news = []
    for bloco in raw.split("\n\nPMID-"):
        tit, journal, pmid = "", "", ""
        for line in bloco.split("\n"):
            if line.startswith("TI  - "): tit = line[6:].strip()
            if line.startswith("JT  - "): journal = line[3:].strip()
            if line.strip().isdigit() and not pmid: pmid = line.strip()
        if tit: news.append({"titulo": tit, "fonte": journal[:30], "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
    return news
