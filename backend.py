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

# ================= GEMINI CORE (ROBUSTO) =================

def montar_url(modelo, chave):
    return f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={chave.strip()}"

def call_gemini(prompt, api_key, temperature=0.0):
    """
    Chamada IA SEM retry automático.
    Falha silenciosa: nunca quebra o pipeline.
    """
    if not api_key:
        return ""

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
                timeout=20
            )

            if resp.status_code == 200:
                data = resp.json()
                return data["candidates"][0]["content"]["parts"][0]["text"]

            if resp.status_code == 429:
                time.sleep(1)
                continue

        except Exception:
            continue

    return ""

# ================= PARSERS ROBUSTOS =================

def clean_text(text):
    if not text:
        return ""
    text = re.sub(r"```.*?```", "", text, flags=re.DOTALL)
    return text.strip()

def parse_structure(text, expected=list):
    text = clean_text(text)
    if not text:
        return expected()

    try:
        obj = json.loads(text)
        if isinstance(obj, expected):
            return obj
    except:
        pass

    try:
        obj = ast.literal_eval(text)
        if isinstance(obj, expected):
            return obj
    except:
        pass

    try:
        pattern = r"\[.*\]" if expected == list else r"\{.*\}"
        match = re.search(pattern, text, re.DOTALL)
        if match:
            obj = ast.literal_eval(match.group(0))
            if isinstance(obj, expected):
                return obj
    except:
        pass

    return expected()

# ================= QUERY INTELLIGENCE =================

def expandir_termo_com_ia(termo, api_key):
    """
    IA só melhora a query.
    Se falhar → usa termo original.
    """
    prompt = f"""
Expand this biomedical concept into a valid PubMed Boolean query.
Use ONLY synonyms.
Return ONLY the query.
TERM: {termo}
"""
    out = call_gemini(prompt, api_key, temperature=0.0)
    out = out.replace('"', '').strip()

    if not out:
        return termo

    if out.count("(") != out.count(")"):
        return termo

    return out

# ================= NORMALIZAÇÃO =================

def normalize_entities(entities):
    canon = []
    for e in entities:
        if not any(SequenceMatcher(None, e.lower(), c.lower()).ratio() > 0.85 for c in canon):
            canon.append(e)
    return canon

# ================= NER EM LOTE (TÍTULO + KEYWORDS) =================

def ner_extraction_batch(artigos, api_key, batch_size=12):
    entidades = []

    for i in range(0, len(artigos), batch_size):
        bloco = artigos[i:i+batch_size]
        texto = "\n".join(a["texto"] for a in bloco)

        prompt = f"""
You are a senior biocurator.
Extract molecular targets and drugs.
Ignore diseases, organs, methods.
Return ONLY a JSON list of strings.

TEXT:
{texto}
"""
        raw = call_gemini(prompt, api_key, temperature=0.1)
        ents = parse_structure(raw, list)
        entidades.extend([str(e).strip() for e in ents if len(str(e)) > 2])

    return entidades

# ================= INFERÊNCIA DE EFEITO =================

def inferir_efeitos_em_lote(artigos, api_key):
    if not artigos:
        return []

    prompt = """
From the PubMed TITLES below, infer pharmacological relationship.
Return ONLY a JSON list of objects:
[{ "title": "...", "target": "...", "drug": "...", "effect": "..." }]
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
    if not email:
        raise ValueError("Email obrigatório")

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

    # --- NER ---
    entidades = ner_extraction_batch(artigos, api_key)
    contagem = Counter(entidades)

    recorrentes = [e for e, c in contagem.items() if c >= 2]
    if not recorrentes:
        recorrentes = [e for e, _ in contagem.most_common(15)]

    alvos = normalize_entities(recorrentes)

    # --- Inferência ---
    artigos_top = artigos[:12]
    inferencias = inferir_efeitos_em_lote(artigos_top, api_key)

    return {
        "termos_indicados": alvos,
        "inferencias": inferencias
    }

# ================= FUNÇÃO ESPERADA PELO FRONT =================

def buscar_alvos_emergentes_pubmed(alvo, email, usar_ia=True):
    api_key = st.session_state.get("api_key_usuario", "") if usar_ia else ""
    res = minerar_pubmed(alvo, email, api_key)
    return res.get("termos_indicados", [])

# ================= CONTADOR PUBMED =================

def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    Entrez.email = email
    query = f"({termo}) AND ({ano_ini}:{ano_fim}[Date - Publication])"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    record = Entrez.read(handle)
    handle.close()
    return int(record["Count"])

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

    handle = Entrez.efetch(
        db="pubmed",
        id=record["IdList"],
        rettype="medline",
        retmode="text"
    )
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
                "fonte": journal[:40],
                "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            })

    return news
