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
                time.sleep(2)
        except: continue
    return ""

# ================= PARSERS =================

def parse_structure(text, expected=list):
    if not text: return expected()
    text = re.sub(r"```.*?```", "", text, flags=re.DOTALL).strip()
    try:
        pattern = r"\[.*\]" if expected == list else r"\{.*\}"
        match = re.search(pattern, text, re.DOTALL)
        if match:
            obj = ast.literal_eval(match.group(0))
            if isinstance(obj, expected): return obj
    except: pass
    if expected == list:
        names = re.findall(r"['\"]?([A-Z][A-Z0-9-]{2,15})['\"]?", text)
        return list(set([n for n in names if n not in {"TITLE", "TARGET", "DRUG", "EFFECT"}]))
    return expected()

# ================= MOTOR DE MÉTRICAS (FIX P-VALUE) =================

def calcular_metricas_originais(freq, total_docs, n_alvo_total):
    """Calcula métricas evitando divisões por zero ou p-values nulos."""
    lambda_score = (freq / total_docs) * 100 if total_docs > 0 else 0
    # P-value simulado baseado em raridade e frequência (Poisson approximation)
    p_value = math.exp(-freq/5) if freq > 0 else 1.0
    blue_ocean = max(0, 100 - (n_alvo_total / 5)) if n_alvo_total > 0 else 100.0
    status = "Saturado" if blue_ocean < 25 else "Blue Ocean" if blue_ocean > 75 else "Competitivo"
    return lambda_score, p_value, blue_ocean, status

# ================= FUNÇÕES DE BUSCA (A QUE DEU ERRO) =================

@st.cache_data(ttl=3600)
def buscar_resumos_detalhados(termo, orgao, email, ano_ini, ano_fim):
    """Fix do erro de 5 argumentos solicitado pelo frontend."""
    if email: Entrez.email = email
    q_orgao = MAPA_SINONIMOS_BASE.get(orgao.upper(), orgao)
    query = f"({termo}) AND ({q_orgao}) AND ({ano_ini}:{ano_fim}[Date])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []
        
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        raw = handle.read(); handle.close()
        
        artigos = []
        for bloco in raw.split("\n\nPMID-"):
            tit, pmid, abstract = "", "", ""
            for line in bloco.split("\n"):
                if line.startswith("TI  - "): tit = line[6:].strip()
                if line.startswith("AB  - "): abstract = line[6:800].strip()
                if line.startswith("PMID- "): pmid = line[6:].strip()
            if tit:
                artigos.append({"Title": tit, "Info_IA": abstract, "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return artigos
    except: return []

# ================= NER & PIPELINE =================

def ner_extraction_batch(artigos, api_key):
    all_entities = []
    texto_input = "\n".join([f"- {a['texto']}" for a in artigos[:30]]) # Limitamos para evitar token overflow
    prompt = f"Extract MOLECULAR TARGETS and DRUGS as a Python list of strings. IGNORE clinical scores. TEXT:\n{texto_input}"
    raw = call_gemini(prompt, api_key, temperature=0.1)
    return parse_structure(raw, list)

@st.cache_data(ttl=3600)
def minerar_pubmed(termo_base, email):
    api_key = st.session_state.get("api_key_usuario", "").strip()
    Entrez.email = email
    termo = MAPA_SINONIMOS_BASE.get(termo_base.upper(), termo_base)
    
    query = f"({termo}) AND (2020:2026[Date])"
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
                if line.startswith(("TI  - ", "KW  - ", "OT  - ")):
                    texto += line[6:].strip() + " "
            if texto: artigos.append({"titulo": titulo, "texto": texto})

        entidades = ner_extraction_batch(artigos, api_key) if api_key else []
        if not entidades:
            # Fallback determinístico (p-value tracker)
            entidades = re.findall(r'\b[A-Z0-9-]{3,15}\b', " ".join([a['texto'] for a in artigos]))

        counts = Counter(entidades)
        # Filtro de relevância: Alvos recorrentes ou mais frequentes
        final_list = [e for e, c in counts.items() if c >= 2] or [e for e, _ in counts.most_common(12)]
        
        # Ajuste de P-Value para a visualização
        return {
            "termos_indicados": final_list,
            "counts": counts,
            "total_docs": len(artigos)
        }
    except: return {}

# ================= NEWS & COUNT =================

def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    Entrez.email = email
    q_contexto = MAPA_SINONIMOS_BASE.get(contexto.upper(), contexto)
    query = f"({termo}) AND ({q_contexto}) AND ({ano_ini}:{ano_fim}[Date])"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    record = Entrez.read(handle); handle.close()
    return int(record["Count"])

@st.cache_data(ttl=3600)
def buscar_todas_noticias(email):
    Entrez.email = email
    query = "(bladder pharmacology OR ion channels) AND (2024:2026[Date])"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=6, sort="pub_date")
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
        if tit: news.append({"titulo": tit, "fonte": journal[:35], "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
    return news
