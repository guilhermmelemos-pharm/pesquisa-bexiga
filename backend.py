"""
Lemos Lambda Backend v2.2 - Ultra Stable Build
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
# Prioridade para o Flash que é mais estável para requisições rápidas
MODELOS_ATIVOS = ["gemini-1.5-flash", "gemini-2.0-flash"]

# --- BLACKLIST CATEGORIZADA ---
BLACKLIST_TOTAL = {
    "metodologia": {"STUDY", "ANALYSIS", "REVIEW", "DATA", "RESULTS", "CONCLUSION", "METHODS", "TRIAL", "RCT", "COHORT"},
    "estatistica": {"P-VALUE", "ANOVA", "RATIO", "ODDS", "STATISTICS", "SIGNIFICANT", "DIFFERENCE", "BASELINE", "SCORE", "KAPLAN"},
    "clinico": {"SURGERY", "RESECTION", "DIAGNOSIS", "PROGNOSIS", "MANAGEMENT", "THERAPY", "TREATMENT", "SYNDROME", "PATIENT", "HOSPITAL", "BIOPSY"},
    "anatomia_ruido": {"KIDNEY", "PROSTATE", "LIVER", "LUNG", "HEART", "BLOOD", "URINE", "CELLS", "PATIENTS"},
    "geral": {"ROLE", "EFFECT", "IMPACT", "POTENTIAL", "NOVEL", "ASSOCIATION", "EVALUATION", "IDENTIFICATION", "ACTIVATION", "POUR", "VOLUME"}
}
UNIFIED_BLACKLIST = set().union(*BLACKLIST_TOTAL.values())

# --- REGEX FARMACÊUTICO REFINADO ---
REGEX_PATTERNS = [
    r'\b[A-Z]{2,5}\d{1,4}[A-Z]?\b',  # TRPV1, P2X3, HDAC6
    r'\b[A-Z]{3,6}R\b',              # AT1R, GPCR
    r'\b[A-Z]{2,4}[- ]?\d{3,6}\b',   # GYY-4137, BAY-123
    r'\b[A-Z][a-z]{3,}(?:ine|mab|ib|ol|on|one|ide|ate|ase|an)\b' # Mirabegron
]

# ================= MOTORES DE IA (VERSÃO ULTRA STABLE) =================

def call_gemini_api(prompt: str, api_key: str, is_json: bool = False) -> str:
    """Motor universal com tratamento de erro detalhado."""
    if not api_key: return ""
    
    headers = {"Content-Type": "application/json"}
    payload = {
        "contents": [{"parts": [{"text": prompt}]}],
        "generationConfig": {"temperature": 0.1 if is_json else 0.3}
    }
    if is_json:
        payload["generationConfig"]["response_mime_type"] = "application/json"

    for modelo in MODELOS_ATIVOS:
        try:
            # URL padronizada com timeout estendido para abstracts longos
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={api_key.strip()}"
            resp = requests.post(url, headers=headers, json=payload, timeout=20)
            
            if resp.status_code == 200:
                data = resp.json()
                if "candidates" in data and len(data["candidates"]) > 0:
                    return data["candidates"][0]["content"]["parts"][0]["text"].strip()
            elif resp.status_code == 429:
                return "Erro: Limite de cota atingido (Rate Limit). Aguarde 60s."
        except Exception as e:
            continue # Tenta o próximo modelo disponível
    return ""

def validar_com_ia(candidatos: List[str], api_key: str, contexto: str) -> List[str]:
    """IA como Juíza: valida os termos minerados."""
    if not candidatos: return []
    prompt = f"Pharmacologist Review for {contexto}. Return ONLY a JSON list of valid Drugs or Molecular Targets from: {', '.join(candidatos[:120])}"
    
    res = call_gemini_api(prompt, api_key, is_json=True)
    try:
        if res:
            # Proteção contra blocos de código markdown que a IA às vezes envia
            clean_text = re.sub(r"```json|```", "", res).strip()
            return json.loads(clean_text)
    except:
        pass
    return []

# ================= PIPELINE PRINCIPAL =================

@st.cache_data(ttl=3600)
def minerar_pubmed(termo_base: str, email: str, usar_ia: bool = True) -> Dict:
    api_key = st.session_state.get("api_key_usuario", "").strip()
    Entrez.email = email
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

        raw_found = []
        for p in REGEX_PATTERNS:
            raw_found.extend(re.findall(p, full_text))
        
        normalized = [f.strip().upper().replace("-", "").replace(" ", "") for f in raw_found]
        cleaned = [n for n in normalized if n not in UNIFIED_BLACKLIST and len(n) > 2]
        
        counts = Counter(cleaned)
        candidates_to_validate = [item for item, _ in counts.most_common(100)]
        
        if api_key and usar_ia:
            validated = validar_com_ia(candidates_to_validate, api_key, termo_base)
            if not validated: validated = candidates_to_validate[:50]
        else:
            validated = candidates_to_validate[:50]

        return {"termos_indicados": validated[:50], "metadata": {"total_articles": total_docs}}
    except Exception as e:
        return {"termos_indicados": [], "error": str(e)}

# ================= ANÁLISE DE ABSTRACT (ÓRGÃO - ALVO - AÇÃO) =================

def analisar_abstract_com_ia(titulo: str, dados_curtos: str, api_key: str, lang: str = 'pt') -> str:
    """Análise dedutiva PhD: Órgão - Alvo - Ação."""
    if not api_key: return "API Key pendente."
    
    prompt = f"""
    Você é uma Doutora em Farmacologia (PhD).
    Deduza e extraia o mecanismo farmacológico seguindo ESTE MODELO de uma única linha:
    Órgão/Tecido - Alvo Molecular - Ação Farmacológica
    
    TÍTULO: "{titulo}"
    RESUMO: "{dados_curtos}"
    """
    
    # Tentativa de resposta com timeout e erro amigável
    try:
        resposta = call_gemini_api(prompt, api_key, is_json=False)
        if resposta:
            # Retorna apenas a primeira linha e remove asteriscos
            return resposta.split('\n')[0].replace("*", "").strip()
        else:
            return "A IA não conseguiu processar este abstract (Cota ou Conexão)."
    except Exception as e:
        return f"Análise indisponível: {str(e)}"

# ================= WRAPPERS COMPATIBILIDADE =================

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

def buscar_todas_noticias(lang_code: str) -> List[Dict]: return []
