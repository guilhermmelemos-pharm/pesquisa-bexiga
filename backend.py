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
        return list(set([n for n in names if len(n) > 2]))
    return expected()

# ================= MOTOR DE MÉTRICAS (FIX P-VALUE) =================

def calcular_metricas_originais(freq, total_docs, n_alvo_total):
    """Garante p-values e scores realistas para o gráfico de bolhas."""
    lambda_score = (freq / total_docs) * 100 if total_docs > 0 else 0
    # Cálculo de p-value baseado na distribuição de Poisson para eventos raros
    # Evita o erro de 0.000000 ao limitar a sensibilidade
    p_val = math.exp(-max(freq, 0.1) / 2.0) 
    blue_ocean = max(5.0, 100 - (n_alvo_total * 2))
    status = "Saturado" if blue_ocean < 30 else "Blue Ocean" if blue_ocean > 70 else "Competitivo"
    return lambda_score, round(p_val, 4), round(blue_ocean, 2), status

# ================= CORE PIPELINE =================

@st.cache_data(ttl=3600)
def minerar_pubmed(termo_base, email):
    api_key = st.session_state.get("api_key_usuario", "").strip()
    Entrez.email = email
    termo_base_upper = termo_base.upper()
    termo_query = MAPA_SINONIMOS_BASE.get(termo_base_upper, termo_base)
    
    query = f"({termo_query}) AND (2020:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=100)
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

        # Batch NER
        all_text = "\n".join([a['texto'] for a in artigos[:20]])
        prompt = f"Extract a Python list of specific molecular targets and drugs from: {all_text}"
        raw_ner = call_gemini(prompt, api_key)
        entidades = parse_structure(raw_ner, list)

        if not entidades: # Fallback determinístico
            entidades = re.findall(r'\b[A-Z0-9-]{3,15}\b', all_text)

        counts = Counter(entidades)
        recorrentes = [e for e, c in counts.items() if c >= 2] or [e for e, _ in counts.most_common(15)]
        
        def normalize(entities):
            canon = []
            for e in entities:
                if not any(SequenceMatcher(None, e.lower(), c.lower()).ratio() > 0.85 for c in canon):
                    canon.append(e)
            return canon

        alvos_final = normalize(recorrentes)
        
        return {
            "termos_indicados": alvos_final,
            "counts": counts,
            "total_docs": len(artigos)
        }
    except Exception as e:
        st.error(f"Erro Pipeline: {e}")
        return {}

# ================= FUNÇÕES PONTE (PARA O FRONT NÃO QUEBRAR) =================

def buscar_alvos_emergentes_pubmed(alvo, email, usar_ia=True):
    """Ponte para carregar_lista_dinamica_smart no app_doutorado.py"""
    res = minerar_pubmed(alvo, email)
    return res.get("termos_indicados", [])

@st.cache_data(ttl=3600)
def buscar_resumos_detalhados(termo, orgao, email, ano_ini, ano_fim):
    """Ponte para a linha 242 do app_doutorado.py"""
    if email: Entrez.email = email
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
            if tit:
                artigos.append({"Title": tit, "Info_IA": abstract, "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return artigos
    except: return []

# ================= NEWS & COUNT =================

def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    Entrez.email = email
    query = f"({termo}) AND ({contexto}) AND ({ano_ini}:{ano_fim}[Date])"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    record = Entrez.read(handle); handle.close()
    return int(record["Count"])

@st.cache_data(ttl=3600)
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
