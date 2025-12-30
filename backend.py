import streamlit as st
from Bio import Entrez, Medline
import requests
import json
import re
from collections import Counter
from typing import List, Dict, Any

# ================= CONFIGURAÇÃO =================
Entrez.email = "pesquisador_guest@unifesp.br"

MODELOS_ATIVOS = ["gemini-2.0-flash", "gemini-1.5-flash"]

# 1. CATEGORIZAÇÃO DA BLACKLIST (Manutenção Sênior)
BLACKLIST_CATEGORIAS = {
    "metodologia": {"STUDY", "ANALYSIS", "REVIEW", "META-ANALYSIS", "DATA", "RESULTS", "CONCLUSION", "METHODS", "AIM", "TRIAL", "RCT", "COHORT"},
    "estatistica": {"P-VALUE", "ANOVA", "RATIO", "ODDS", "CONFIDENCE", "INTERVAL", "STATISTICS", "SIGNIFICANT", "DIFFERENCE", "BASELINE", "SCORE"},
    "clinico": {"SURGERY", "RESECTION", "DIAGNOSIS", "PROGNOSIS", "MANAGEMENT", "THERAPY", "TREATMENT", "PROTOCOL", "SYNDROME", "PATIENT", "HOSPITAL"},
    "geral": {"ROLE", "EFFECT", "IMPACT", "POTENTIAL", "NOVEL", "ASSOCIATION", "EVALUATION", "IDENTIFICATION", "ACTIVATION", "EXPRESSION", "PATHWAY"}
}

# 2. REGEX DETERMINÍSTICO (Farmacologia Real)
# Pega padrões como TRPV1, P2X3, mTOR, GYY-4137, GSK12345
REGEX_FARMA_ALVOS = [
    r'\b[A-Z]{2,5}\d{1,4}\b',        # P2X3, TRPV1, HDAC6
    r'\b[A-Z]{3,6}R\b',              # AT1R, GPCR
    r'\b[A-Z]{2,4}[- ]?\d{3,6}\b',   # GYY-4137, BAY 123
    r'\b[A-Z][a-z]{3,}(?:ine|mab|ib|ol|on|one|ide|ate|ase|an)\b' # Fármacos (Mirabegron, etc)
]

MAPA_SINONIMOS_BASE = {
    "BLADDER": "(Bladder OR Urothelial OR Urothelium OR Detrusor)",
    "PAIN": "(Pain OR Nociception OR Analgesia)",
    "INFLAMMATION": "(Inflammation OR Cytokines OR NF-kappaB)",
    "BRAIN": "(Brain OR Cerebral OR CNS OR Neuron OR Glia)"
}

# ================= MOTOR DE IA =================

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
                parsed = json.loads(clean_text)
                return parsed if isinstance(parsed, list) else []
        except: continue
    return []

# ================= PIPELINE DE EXTRAÇÃO =================

def validar_candidatos_com_ia(candidatos: List[str], api_key: str, contexto: str) -> List[str]:
    """
    IA atua como VALIDADORA de candidatos pré-selecionados pelo Regex.
    Isso economiza tokens e aumenta a precisão.
    """
    if not candidatos: return []
    
    lista_texto = ", ".join(candidatos[:100])
    prompt = f"""
    You are a Senior Pharmacologist. 
    Context: {contexto}.
    Below is a list of potential molecular targets and drugs.
    Task: Return ONLY the items that are real Drugs, Chemical Compounds, or Molecular Targets (Receptors, Channels, Enzymes).
    Discard verbs, adjectives, methods, or anatomy.
    Return a pure JSON list of strings.
    List: {lista_texto}
    """
    return call_gemini_json(prompt, api_key)

@st.cache_data(ttl=3600)
def minerar_pubmed(termo_base: str, email: str, usar_ia: bool = True) -> Dict:
    api_key = st.session_state.get("api_key_usuario", "").strip()
    Entrez.email = email
    
    query = f"({MAPA_SINONIMOS_BASE.get(termo_base.upper(), termo_base)}) AND (drug OR inhibitor OR compound) AND (2020:2026[Date])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=200)
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return {}

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        
        texto_full = ""
        artigos_completos = []
        for r in records:
            txt = f"{r.get('TI', '')} {r.get('AB', '')} {' '.join(r.get('OT', []))}"
            texto_full += " " + txt
            artigos_completos.append({"titulo": r.get('TI', ''), "texto": txt})

        # 1. Extração Determinística (Regex)
        candidatos_brutos = []
        for pattern in REGEX_FARMA_ALVOS:
            candidatos_brutos.extend(re.findall(pattern, texto_full))
        
        # 2. Limpeza Básica (Blacklist Universal)
        blacklist_total = set().union(*BLACKLIST_CATEGORIAS.values())
        candidatos_limpos = [
            c for c in candidatos_brutos 
            if c.upper() not in blacklist_total and len(c) > 2 and not c.isdigit()
        ]
        
        candidatos_unicos = list(set(candidatos_limpos))
        
        # 3. Validação por IA (O Juiz Final)
        if api_key and usar_ia and candidatos_unicos:
            final_list = validar_candidatos_com_ia(candidatos_unicos, api_key, termo_base)
        else:
            # Fallback: Se sem IA, usa os mais frequentes do Regex
            counts = Counter(candidatos_limpos)
            final_list = [item for item, count in counts.most_common(50)]

        return {
            "termos_indicados": final_list[:50],
            "total_docs": len(artigos_completos),
            "artigos_originais": artigos_completos
        }
    except Exception as e:
        print(f"Erro: {e}")
        return {}

# ================= WRAPPERS (Mantendo compatibilidade com Frontend) =================

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
    prompt = f"Pharmacology PhD. Analyze abstract: {titulo}. Context: {dados_curtos}. Output: Organ - Target - Action."
    headers = {"Content-Type": "application/json"}
    payload = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.3}}
    try:
        url = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash:generateContent?key={api_key.strip()}"
        resp = requests.post(url, headers=headers, json=payload, timeout=10)
        return resp.json()["candidates"][0]["content"]["parts"][0]["text"].strip()
    except: return "Erro na análise."

def buscar_todas_noticias(lang_code: str) -> List[Dict]: return []
