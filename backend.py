# backend.py
import streamlit as st
from Bio import Entrez
import requests
import re
import time
import ast
import json
from collections import Counter
from difflib import SequenceMatcher

# ================= CONFIG =================
Entrez.email = "pesquisador_guest@unifesp.br"

MODELOS_ATIVOS = [
    "gemini-2.0-flash",
    "gemini-2.0-flash-exp"
]

MAPA_SINONIMOS_BASE = {
    "BLADDER": "(Bladder OR Vesical OR Detrusor OR Urothelium)",
    "PAIN": "(Pain OR Nociception OR Analgesia)",
    "INFLAMMATION": "(Inflammation OR Cytokines OR Inflammasome)"
}

# ================= GEMINI CORE =================

def montar_url(modelo, chave):
    return f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={chave.strip()}"

def call_gemini(prompt, api_key, temperature=0.0):
    headers = {"Content-Type": "application/json"}
    payload = {
        "contents": [{"parts": [{"text": prompt}]}],
        "generationConfig": {"temperature": temperature}
    }

    for modelo in MODELOS_ATIVOS:
        try:
            resp = requests.post(
                montar_url(modelo, api_key),
                headers=headers,
                json=payload,
                timeout=25
            )
            if resp.status_code == 200:
                return resp.json()["candidates"][0]["content"]["parts"][0]["text"]
            elif resp.status_code == 429:
                time.sleep(2)
        except:
            pass
    return ""

# ================= PARSERS =================

def clean_text(text):
    if not text:
        return ""
    text = re.sub(r"```.*?```", "", text, flags=re.DOTALL)
    return text.strip()

def parse_structure(text, expected=list):
    text = clean_text(text)
    try:
        data = json.loads(text)
        if isinstance(data, expected):
            return data
    except:
        pass
    try:
        data = ast.literal_eval(text)
        if isinstance(data, expected):
            return data
    except:
        pass
    try:
        pattern = r"\[.*\]" if expected == list else r"\{.*\}"
        match = re.search(pattern, text, re.DOTALL)
        if match:
            data = ast.literal_eval(match.group())
            if isinstance(data, expected):
                return data
    except:
        pass
    return expected()

# ================= QUERY INTELLIGENCE =================

def expandir_termo_com_ia(termo, api_key):
    prompt = f"""
Expand this biomedical concept into a valid PubMed Boolean query.
Return ONLY the query string.

TERM: {termo}
"""
    r = call_gemini(prompt, api_key, temperature=0.0)
    if r and r.count("(") == r.count(")"):
        return r.replace('"', '').strip()
    return termo

# ================= NORMALIZAÇÃO =================

def normalize_locally(entities):
    canonical = []
    for e in entities:
        if not any(SequenceMatcher(None, e.lower(), c.lower()).ratio() > 0.85 for c in canonical):
            canonical.append(e)
    return canonical

# ================= BATCH NER =================

def ner_extraction_batch(artigos, api_key, batch_size=15):
    all_entities = []

    for i in range(0, len(artigos), batch_size):
        batch = artigos[i:i + batch_size]
        texto = "\n".join([f"- {a['texto']}" for a in batch])

        prompt = f"""
You are a PhD-level biocurator.
Extract molecular targets (receptors, channels, enzymes) and drugs.
Ignore organs, diseases, scores, and methods.

Return ONLY a Python list of strings.

TEXT:
{texto}
"""
        raw = call_gemini(prompt, api_key, temperature=0.1)
        ents = parse_structure(raw, list)
        all_entities.extend([str(e).strip() for e in ents if len(str(e)) > 2])

    return all_entities

# ================= BATCH EFFECT INFERENCE =================

def inferir_efeitos_em_lote(artigos, api_key):
    if not artigos:
        return []

    prompt = """
Infer pharmacological effect ONLY from titles.

Return a Python list of dicts:
{ "title": "...", "target": "...", "drug": "...", "effect": "..." }

Use "N/A" if unclear.

TITLES:
"""
    for a in artigos:
        prompt += f"- {a['titulo']}\n"

    raw = call_gemini(prompt, api_key, temperature=0.0)
    return parse_structure(raw, list)

# ================= CORE PIPELINE =================

@st.cache_data(ttl=3600)
def minerar_pubmed(termo_base, email, api_key):
    if not email or not api_key:
        return {}

    Entrez.email = email

    termo = MAPA_SINONIMOS_BASE.get(termo_base.upper(), termo_base)
    termo = expandir_termo_com_ia(termo, api_key)
    query = f"({termo}) AND (2020:2026[Date - Publication])"

    handle = Entrez.esearch(db="pubmed", term=query, retmax=200)
    record = Entrez.read(handle)
    handle.close()

    if not record["IdList"]:
        return {}

    handle = Entrez.efetch(
        db="pubmed",
        id=record["IdList"],
        rettype="medline",
        retmode="text"
    )
    raw = handle.read()
    handle.close()

    artigos = []
    for bloco in raw.split("\n\nPMID-"):
        titulo, texto = "", ""
        for line in bloco.split("\n"):
            if line.startswith("TI  - "):
                titulo = line[6:].strip()
            if line.startswith(("TI  - ", "KW  - ", "OT  - ")):
                texto += line[6:].strip() + " "
        if texto:
            artigos.append({"titulo": titulo, "texto": texto})

    # === NER ===
    entidades = ner_extraction_batch(artigos, api_key)
    counts = Counter(entidades)

    recorrentes = [e for e, c in counts.items() if c >= 2]
    if not recorrentes:
        recorrentes = [e for e, _ in counts.most_common(15)]

    alvos_final = normalize_locally(recorrentes)

    # === INFERENCE ===
    artigos_top = artigos[:12]
    inferencias = inferir_efeitos_em_lote(artigos_top, api_key)

    return {
        "termos_indicados": alvos_final,
        "inferencias": inferencias
    }

# ================= COMPATIBILITY LAYER (FRONTEND ANTIGO) =================

def buscar_alvos_emergentes_pubmed(alvo, email, usar_ia=True):
    api_key = st.session_state.get("api_key_usuario", "")
    res = minerar_pubmed(alvo, email, api_key)
    return res.get("termos_indicados", [])

# ================= NEWS RADAR =================

@st.cache_data(ttl=3600)
def buscar_todas_noticias(email):
    if not email:
        return []

    Entrez.email = email
    query = "(molecular pharmacology OR ion channels OR signaling pathway) AND (2024:2026[Date - Publication])"

    handle = Entrez.esearch(db="pubmed", term=query, retmax=6, sort="pub_date")
    record = Entrez.read(handle)
    handle.close()

    if not record["IdList"]:
        return []

    handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
    raw = handle.read()
    handle.close()

    news = []
    for bloco in raw.split("\n\nPMID-"):
        tit, journal, pmid = "", "", ""
        for line in bloco.split("\n"):
            if line.startswith("TI  - "):
                tit = line[6:].strip()
            if line.startswith("JT  - "):
                journal = line[3:].strip()
            if line.strip().isdigit() and not pmid:
                pmid = line.strip()
        if tit:
            news.append({
                "titulo": tit,
                "fonte": journal[:35],
                "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            })
    return news
