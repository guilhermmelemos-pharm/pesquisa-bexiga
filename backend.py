"""
Lemos Lambda Backend v2.1
Stable Hybrid Deterministic + LLM Pipeline
"""

import streamlit as st
from Bio import Entrez, Medline
import requests, json, re
from collections import Counter
from typing import List, Dict

# ================= CONFIG =================
Entrez.email = "pesquisador_guest@unifesp.br"

MODEL_JSON = "gemini-1.5-flash"
MODEL_TEXT = "gemini-1.5-flash"

GEMINI_URL = "https://generativelanguage.googleapis.com/v1beta/models/{model}:generateContent?key={key}"

# ================= BLACKLIST =================
BLACKLIST_TOTAL = {
    "geral": {
        "STUDY","ANALYSIS","REVIEW","DATA","RESULTS","METHODS",
        "CLINICAL","TRIAL","PATIENT","SIGNIFICANT","EFFECT"
    }
}
UNIFIED_BLACKLIST = set().union(*BLACKLIST_TOTAL.values())

# ================= REGEX =================
REGEX_PATTERNS = [
    r'\b[A-Z]{2,6}\d{1,4}\b',
    r'\b[A-Z]{3,6}R\b',
    r'\b[A-Z][a-z]{4,}(?:ine|mab|ib|one|ol|ate|ase)\b'
]

# ================= GEMINI CORE =================

def gemini_request(prompt: str, api_key: str, expect_json: bool) -> str:
    if not api_key:
        return ""

    payload = {
        "contents": [{"parts": [{"text": prompt}]}],
        "generationConfig": {
            "temperature": 0.1 if expect_json else 0.3,
            "maxOutputTokens": 512
        }
    }

    if expect_json:
        payload["generationConfig"]["response_mime_type"] = "application/json"

    url = GEMINI_URL.format(
        model=MODEL_JSON if expect_json else MODEL_TEXT,
        key=api_key.strip()
    )

    try:
        r = requests.post(url, json=payload, timeout=35)
        if r.status_code != 200:
            return ""

        data = r.json()
        return data["candidates"][0]["content"]["parts"][0]["text"]

    except:
        return ""

# ================= IA VALIDADORA =================

def validar_com_ia(candidatos: List[str], api_key: str, contexto: str) -> List[str]:
    if not candidatos:
        return []

    prompt = f"""
    You are a senior pharmacologist.
    Context: {contexto}

    From the list below, return ONLY valid:
    - Drugs
    - Chemical compounds
    - Molecular targets (receptors, enzymes, channels)

    Return STRICT JSON array. No comments.

    List: {", ".join(candidatos[:80])}
    """

    raw = gemini_request(prompt, api_key, expect_json=True)

    try:
        clean = re.sub(r"```json|```", "", raw).strip()
        return json.loads(clean)
    except:
        return candidatos[:40]

# ================= PUBMED MINER =================

@st.cache_data(ttl=3600)
def minerar_pubmed(termo: str, email: str, usar_ia=True) -> Dict:
    api_key = st.session_state.get("api_key_usuario", "").strip()
    Entrez.email = email

    query = f"({termo}) AND (drug OR inhibitor OR target) AND (2020:2026[Date])"

    try:
        h = Entrez.esearch(db="pubmed", term=query, retmax=200)
        ids = Entrez.read(h)["IdList"]
        h.close()

        if not ids:
            return {}

        h = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
        records = Medline.parse(h)

        corpus = ""
        total = 0
        for r in records:
            corpus += f" {r.get('TI','')} {r.get('AB','')}"
            total += 1

        found = []
        for p in REGEX_PATTERNS:
            found.extend(re.findall(p, corpus))

        cleaned = [
            f.upper().replace("-", "")
            for f in found
            if f.upper() not in UNIFIED_BLACKLIST and len(f) > 2
        ]

        ranked = [k for k,_ in Counter(cleaned).most_common(100)]

        if usar_ia and api_key:
            ranked = validar_com_ia(ranked, api_key, termo)

        return {
            "termos_indicados": ranked[:50],
            "metadata": {"total_articles": total}
        }

    except:
        return {}

# ================= ABSTRACT ANALYSIS =================

def analisar_abstract_com_ia(titulo: str, resumo: str, api_key: str) -> str:
    if not api_key:
        return "API Key não configurada."

    prompt = f"""
    Você é uma Doutora em Farmacologia.

    Analise e deduza no formato EXATO:
    Órgão/Tecido - Alvo Molecular - Ação Farmacológica

    TÍTULO: {titulo}
    RESUMO: {resumo}

    Responda APENAS uma linha.
    """

    resposta = gemini_request(prompt, api_key, expect_json=False)

    if not resposta:
        return "Análise indisponível (timeout ou instabilidade da API)."

    return resposta.split("\n")[0].replace("*","").strip()

# ================= WRAPPERS =================

def buscar_alvos_emergentes_pubmed(alvo, email, usar_ia=True):
    return minerar_pubmed(alvo, email, usar_ia).get("termos_indicados", [])

