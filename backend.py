"""
Lemos Lambda Backend v3.1
Hybrid Deterministic-LLM Pipeline for Pharmacological Discovery
"""
import streamlit as st
from Bio import Entrez, Medline
import requests
import json
import re
from collections import Counter
from typing import List, Dict, Any

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br"
MODELOS_ATIVOS = ["gemini-2.0-flash", "gemini-1.5-flash"]

# --- BLACKLIST CATEGORIZADA (Mínima e Eficiente) ---
BLACKLIST_TOTAL = {
    "metodologia": {"STUDY", "ANALYSIS", "REVIEW", "DATA", "RESULTS", "CONCLUSION", "METHODS", "TRIAL", "RCT", "COHORT"},
    "estatistica": {"P-VALUE", "ANOVA", "RATIO", "ODDS", "STATISTICS", "SIGNIFICANT", "DIFFERENCE", "BASELINE", "SCORE"},
    "clinico": {"SURGERY", "RESECTION", "DIAGNOSIS", "PROGNOSIS", "MANAGEMENT", "THERAPY", "TREATMENT", "SYNDROME", "PATIENT"},
    "geral": {"ROLE", "EFFECT", "IMPACT", "POTENTIAL", "NOVEL", "ASSOCIATION", "EVALUATION", "IDENTIFICATION", "ACTIVATION"}
}
UNIFIED_BLACKLIST = set().union(*BLACKLIST_TOTAL.values())

# --- REGEX FARMACÊUTICO (Acurácia Técnica) ---
REGEX_PATTERNS = [
    r'\b[A-Z]{2,5}\d{1,4}\b',        # Ex: P2X3, TRPV1, HDAC6
    r'\b[A-Z]{3,6}R\b',              # Ex: AT1R, GPCR
    r'\b[A-Z]{2,4}[- ]?\d{3,6}\b',   # Ex: GYY-4137, BAY 123
    r'\b[A-Z][a-z]{3,}(?:ine|mab|ib|ol|on|one|ide|ate|ase|an)\b' # Ex: Mirabegron
]

# ================= PIPELINE DE PROCESSAMENTO =================

def normalizar_entidade(e: str) -> str:
    """Padronização para evitar duplicatas (TRPV-1 == TRPV1)."""
    return e.strip().upper().replace("-", "").replace(" ", "")

def call_gemini_json(prompt: str, api_key: str) -> List[str]:
    if not api_key: return []
    headers = {"Content-Type": "application/json"}
    payload = {
        "contents": [{"parts": [{"text": prompt}]}], 
        "generationConfig": {"temperature": 0.1, "response_mime_type": "application/json"}
    }
    for modelo in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={api_key.strip()}"
            resp = requests.post(url, headers=headers, json=payload, timeout=20)
            if resp.status_code == 200:
                text = resp.json()["candidates"][0]["content"]["parts"][0]["text"]
                clean_text = re.sub(r"```json|```", "", text).strip()
                return json.loads(clean_text)
        except: continue
    return []

def validar_com_ia(candidatos: List[str], api_key: str, contexto: str) -> List[str]:
    """IA como Juíza: Filtra ruído e valida farmacologia."""
    if not candidatos: return []
    prompt = f"""
    Expert Pharmacologist Review. Context: {contexto}.
    Task: From the list below, return ONLY valid Drugs, Chemical Compounds, or Molecular Targets (Receptors/Channels/Enzymes).
    Format: Pure JSON string list.
    List: {", ".join(candidatos[:120])}
    """
    return call_gemini_json(prompt, api_key)

@st.cache_data(ttl=3600)
def minerar_pubmed(termo_base: str, email: str, usar_ia: bool = True) -> Dict:
    api_key = st.session_state.get("api_key_usuario", "").strip()
    Entrez.email = email
    
    # Query otimizada
    query = f"({termo_base}) AND (drug OR inhibitor OR compound OR target) AND (2020:2026[Date])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=200)
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return {}

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        
        full_text = ""
        total_docs = 0
        for r in records:
            full_text += f" {r.get('TI','')} {r.get('AB','')} {' '.join(r.get('OT',[]))}"
            total_docs += 1

        # STEP 1: Extração Determinística
        raw_found = []
        for p in REGEX_PATTERNS:
            raw_found.extend(re.findall(p, full_text))
        
        # STEP 2: Normalização e Limpeza
        normalized = [normalizar_entidade(f) for f in raw_found]
        cleaned = [n for n in normalized if n not in UNIFIED_BLACKLIST and len(n) > 2]
        
        # STEP 3: Ranking por Frequência Relativa (Valor Científico)
        counts = Counter(cleaned)
        # Pega os 100 mais frequentes para validação
        candidates_to_validate = [item for item, _ in counts.most_common(100)]
        
        # STEP 4: Validação (LLM vs Fallback)
        if api_key and usar_ia:
            validated = validar_com_ia(candidates_to_validate, api_key, termo_base)
        else:
            validated = candidates_to_validate[:50]

        return {
            "termos_indicados": validated[:50],
            "metadata": {"total_articles": total_docs, "raw_candidates": len(candidates_to_validate)}
        }
    except Exception: return {}

# --- WRAPPERS DE COMPATIBILIDADE ---
def buscar_alvos_emergentes_pubmed(alvo: str, email: str, usar_ia: bool = True) -> List[str]:
    res = minerar_pubmed(alvo, email, usar_ia=usar_ia)
    return res.get("termos_indicados", [])

def consultar_pubmed_count(termo: str, contexto: str, email: str, ano_ini: int, ano_fim: int) -> int:
    Entrez.email = email
    query = f"({termo}) AND ({contexto}) AND ({ano_ini}:{ano_fim}[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle); handle.close()
        return int(record["Count"])
    except: return 0

@st.cache_data(ttl=3600)
def buscar_resumos_detalhados(termo: str, orgao: str, email: str, ano_ini: int, ano_fim: int) -> List[Dict]:
    Entrez.email = email
    query = f"({termo}) AND ({orgao}) AND ({ano_ini}:{ano_fim}[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=6)
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        return [{"Title": r.get("TI", ""), "Info_IA": r.get("AB", "")[:1000], "Link": f"https://pubmed.ncbi.nlm.nih.gov/{r.get('PMID', '')}/"} for r in records]
    except: return []

def analisar_abstract_com_ia(titulo: str, dados_curtos: str, api_key: str, lang: str = 'pt') -> str:
    if not api_key: return "API Key pendente."
    prompt = f"Pharmacology PhD. Analyze: {titulo}. Context: {dados_curtos}. Output: Organ - Target - Action."
    try:
        url = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash:generateContent?key={api_key.strip()}"
        resp = requests.post(url, headers={"Content-Type": "application/json"}, json={"contents": [{"parts": [{"text": prompt}]}]}, timeout=10)
        return resp.json()["candidates"][0]["content"]["parts"][0]["text"].strip()
    except: return "Análise indisponível."

def buscar_todas_noticias(lang_code: str) -> List[Dict]: return []
