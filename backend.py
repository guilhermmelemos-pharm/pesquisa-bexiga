"""
Lemos Lambda Backend v2.0
Hybrid Deterministic-LLM Pipeline for Pharmacological Discovery
SDK Oficial Gemini 2.5 (Hibridização Pro/Flash para Streamlit)
"""

import streamlit as st
from Bio import Entrez, Medline
import google.generativeai as genai
import json
import re
from collections import Counter
from typing import List, Dict

# ================= CONFIGURAÇÃO =================

Entrez.email = "pesquisador_guest@unifesp.br"

# Definições Canônicas conforme decisão final
MODELO_PRO = "gemini-2.5-pro"
MODELO_FLASH = "gemini-2.5-flash"

# ================= BLACKLIST CATEGORIZADA =================

BLACKLIST_TOTAL = {
    "metodologia": {"STUDY", "ANALYSIS", "REVIEW", "DATA", "RESULTS", "CONCLUSION", "METHODS", "TRIAL", "RCT", "COHORT"},
    "estatistica": {"P-VALUE", "ANOVA", "RATIO", "ODDS", "STATISTICS", "SIGNIFICANT", "DIFFERENCE", "BASELINE", "SCORE", "KAPLAN"},
    "clinico": {"SURGERY", "RESECTION", "DIAGNOSIS", "PROGNOSIS", "MANAGEMENT", "THERAPY", "TREATMENT", "SYNDROME", "PATIENT", "HOSPITAL", "BIOPSY"},
    "anatomia_ruido": {"KIDNEY", "PROSTATE", "LIVER", "LUNG", "HEART", "BLOOD", "URINE", "CELLS", "PATIENTS"},
    "geral": {"ROLE", "EFFECT", "IMPACT", "POTENTIAL", "NOVEL", "ASSOCIATION", "EVALUATION", "IDENTIFICATION", "ACTIVATION"}
}
UNIFIED_BLACKLIST = set().union(*BLACKLIST_TOTAL.values())

# ================= REGEX FARMACÊUTICO =================

REGEX_PATTERNS = [
    r'\b[A-Z]{2,5}\d{1,4}[A-Z]?\b',
    r'\b[A-Z]{3,6}R\b',
    r'\b[A-Z]{2,4}[- ]?\d{3,6}\b',
    r'\b[A-Z][a-z]{3,}(?:ine|mab|ib|ol|on|one|ide|ate|ase|an)\b'
]

# ================= MOTOR HÍBRIDO GEMINI SDK =================

def configurar_gemini(api_key: str):
    if api_key:
        genai.configure(api_key=api_key.strip())

def gerar_com_gemini(prompt: str, is_json: bool = False) -> str:
    """
    Roteamento Inteligente: 
    - JSON/Validação -> Flash Primeiro (Velocidade)
    - Dedução/Abstract -> Pro Primeiro (Raciocínio PhD)
    """
    safety_settings = [
        {"category": "HARM_CATEGORY_HARASSMENT", "threshold": "BLOCK_NONE"},
        {"category": "HARM_CATEGORY_HATE_SPEECH", "threshold": "BLOCK_NONE"},
        {"category": "HARM_CATEGORY_SEXUALLY_EXPLICIT", "threshold": "BLOCK_NONE"},
        {"category": "HARM_CATEGORY_DANGEROUS_CONTENT", "threshold": "BLOCK_NONE"},
    ]

    generation_config = {
        "temperature": 0.05 if is_json else 0.3,
        "max_output_tokens": 1024,
    }

    # Lógica de priorização H2O (Hybrid High-Output)
    modelos_prioridade = [MODELO_FLASH, MODELO_PRO] if is_json else [MODELO_PRO, MODELO_FLASH]

    for modelo_nome in modelos_prioridade:
        try:
            model = genai.GenerativeModel(
                model_name=modelo_nome,
                generation_config=generation_config,
                safety_settings=safety_settings
            )
            response = model.generate_content(prompt)
            if response and response.text:
                return response.text.strip()
        except Exception:
            continue # Fallback automático para o próximo modelo da lista

    return ""

# ================= VALIDAÇÃO COM IA (FLASH DRIVEN) =================

def validar_com_ia(candidatos: List[str], contexto: str) -> List[str]:
    if not candidatos:
        return []

    prompt = f"""
    Expert Pharmacologist Review. Context: {contexto}
    Task: Return ONLY a pure JSON list of valid Drugs or Molecular Targets (Receptors, Channels, Enzymes).
    LIST: {", ".join(candidatos[:120])}
    """

    resposta = gerar_com_gemini(prompt, is_json=True)

    try:
        clean = re.sub(r"```json|```", "", resposta).strip()
        return json.loads(clean)
    except Exception:
        return []

# ================= PIPELINE PUBMED =================

@st.cache_data(ttl=3600)
def minerar_pubmed(termo_base: str, email: str, usar_ia: bool = True) -> Dict:
    api_key = st.session_state.get("api_key_usuario", "").strip()
    configurar_gemini(api_key)

    Entrez.email = email
    query = f"({termo_base}) AND (drug OR inhibitor OR compound OR target) AND (2020:2026[Date])"

    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=200)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            return {}

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)

        full_text = ""
        total_docs = 0
        for r in records:
            full_text += f" {r.get('TI','')} {r.get('AB','')} {' '.join(r.get('OT',[]))}"
            total_docs += 1

        raw_found = []
        for p in REGEX_PATTERNS:
            raw_found.extend(re.findall(p, full_text))

        normalized = [f.upper().replace("-", "").replace(" ", "") for f in raw_found]
        cleaned = [n for n in normalized if n not in UNIFIED_BLACKLIST and len(n) > 2]

        counts = Counter(cleaned)
        candidatos = [item for item, _ in counts.most_common(120)]

        if api_key and usar_ia:
            validados = validar_com_ia(candidatos, termo_base)
            if not validados: validados = candidatos[:50]
        else:
            validados = candidatos[:50]

        return {"termos_indicados": validados[:50], "metadata": {"total_articles": total_docs}}

    except Exception as e:
        return {"termos_indicados": [], "error": str(e)}

# ================= ANÁLISE DE ABSTRACT (PRO DRIVEN) =================

def analisar_abstract_com_ia(titulo: str, dados_curtos: str, api_key: str) -> str:
    if not api_key: return "API Key pendente."
    configurar_gemini(api_key)

    prompt = f"""
    You are a PhD Pharmacologist. Deduce the pharmacological mechanism.
    Return EXACTLY one line: Organ/Tissue - Molecular Target - Pharmacological Action
    TITLE: "{titulo}"
    ABSTRACT: "{dados_curtos}"
    """

    resposta = gerar_com_gemini(prompt, is_json=False)
    if resposta:
        return resposta.split("\n")[0].replace("*", "").strip()

    return "Análise indisponível no momento."

# ================= WRAPPERS =================

def buscar_alvos_emergentes_pubmed(alvo: str, email: str, usar_ia: bool = True) -> List[str]:
    res = minerar_pubmed(alvo, email, usar_ia)
    return res.get("termos_indicados", [])

def consultar_pubmed_count(termo: str, contexto: str, email: str, ano_ini: int, ano_fim: int) -> int:
    Entrez.email = email
    query = f"({termo}) AND ({contexto}) AND ({ano_ini}:{ano_fim}[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle)
        handle.close()
        return int(record["Count"])
    except Exception: return 0

@st.cache_data(ttl=3600)
def buscar_resumos_detalhados(termo: str, orgao: str, email: str, ano_ini: int, ano_fim: int) -> List[Dict]:
    Entrez.email = email
    query = f"({termo}) AND ({orgao}) AND ({ano_ini}:{ano_fim}[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=6)
        record = Entrez.read(handle)
        handle.close()
        if not record["IdList"]: return []
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        return [{"Title": r.get("TI", ""), "Info_IA": r.get("AB", "")[:1000], "Link": f"https://pubmed.ncbi.nlm.nih.gov/{r.get('PMID', '')}/"} for r in records]
    except Exception: return []

def buscar_todas_noticias(lang_code: str) -> List[Dict]: return []
