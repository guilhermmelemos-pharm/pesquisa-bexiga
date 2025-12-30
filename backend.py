import streamlit as st
from Bio import Entrez, Medline
import requests
import json
import re
import time
from collections import Counter
from typing import List, Dict, Any

# ================= CONFIGURAÇÃO =================
Entrez.email = "pesquisador_guest@unifesp.br"

MODELOS_ATIVOS = [
    "gemini-2.5-flash", 
    "gemini-2.0-flash", 
    "gemini-1.5-flash"
]

# LISTA NEGRA UNIVERSAL: O Exterminador de Ruído Metodológico
BLACKLIST_MONSTRO = {
    # 1. Termos Metodológicos e de Publicação
    "STUDY", "ANALYSIS", "REVIEW", "META-ANALYSIS", "DATA", "RESULTS", "CONCLUSION",
    "BACKGROUND", "METHODS", "OBJECTIVE", "AIM", "HYPOTHESIS", "INTRODUCTION",
    "CASE", "REPORT", "SERIES", "SURVEY", "QUESTIONNAIRE", "TRIAL", "RCT", "COHORT",
    "RETROSPECTIVE", "PROSPECTIVE", "CROSS-SECTIONAL", "MULTICENTER", "PILOT",
    
    # 2. Estatística e Métricas (Ruído Matemático)
    "P-VALUE", "ANOVA", "RATIO", "ODDS", "CONFIDENCE", "INTERVAL", "STATISTICS",
    "SIGNIFICANT", "DIFFERENCE", "INCREASED", "DECREASED", "LEVELS", "SCORE", 
    "SCALE", "INDEX", "MENDELIAN", "RANDOMIZATION", "MULTIVARIATE", "UNIVARIATE", 
    "VALIDATION", "BASELINE", "PREDICTION", "SIMULATION", "PRECISION", "ACCURACY",
    "SENSITIVITY", "SPECIFICITY", "PERCENT", "TOTAL", "SAMPLE", "SIZE",
    
    # 3. Verbos e Termos de Processo Biológico Genérico
    "ROLE", "EFFECT", "IMPACT", "POTENTIAL", "NOVEL", "ASSOCIATION", "EVALUATION",
    "IDENTIFICATION", "ACTIVATION", "DIVISION", "REGULATION", "FUNCTION", "ACTION",
    "PATHOGENESIS", "DEVELOPMENT", "PROGRESSION", "CHARACTERIZATION", "INVESTIGATION",
    "MODULATION", "INTERACTION", "OBSERVATION", "DEMONSTRATION", "DISTRIBUTION",
    "EXPRESSION", "ACCUMULATION", "ORIGIN", "PRESERVATION", "CORRECTION", "SYNERGY",
    
    # 4. Contexto Clínico e Hospitalar (Não são alvos)
    "SURGERY", "RESECTION", "INCISION", "OPERATION", "TRANSPLANT", "GRAFT", "STENT",
    "CATHETER", "BIOPSY", "IMAGING", "MRI", "CT", "PET", "ULTRASOUND", "DIAGNOSIS",
    "PROGNOSIS", "MANAGEMENT", "THERAPY", "TREATMENT", "PROTOCOL", "GUIDELINE",
    "COMPLICATION", "INFECTION", "DYSFUNCTION", "SYNDROME", "DISORDER", "DISEASE",
    "PATIENT", "PARTICIPANT", "CHILDREN", "ADULT", "WOMEN", "MEN", "ELDERLY",
    "HOSPITAL", "CLINIC", "CENTER", "DEPARTMENT", "UNIVERSITY", "OUTCOME", "SAFETY",
    "EFFICACY", "MORTALITY", "SURVIVAL", "ADMISSION", "DISCHARGE",
    
    # 5. Termos Biológicos Genéricos (Apenas Ruído de fundo)
    "DNA", "RNA", "MRNA", "PROTEIN", "CELL", "TISSUE", "SERUM", "PLASMA", "URINE", "BLOOD",
    "HISTONE", "FACTOR", "COMPONENT", "SYSTEM", "MODEL", "PATHWAY", "MECHANISM",
    "GUERIN", "BACILLUS", "CALMETTE", "BCG", "PLACEBO", "CONTROL", "SHAM", "VEHICLE", 
    "SALINE", "BUFFER", "IN-VITRO", "IN-VIVO", "ANIMAL", "HUMAN", "MICE", "RAT",
    
    # 6. Outros (Lixo de extração comum)
    "COVID-19", "COVID", "SARS-COV-2", "PANDEMIC", "VIRUS", "YEAR", "MONTH", "DAY",
    "HIGH", "LOW", "NEW", "OLD", "FIRST", "SECOND", "THIRD", "AMONG", "BETWEEN"
}

# Mapa de Sinônimos focado em melhorar a busca inicial do PubMed
MAPA_SINONIMOS_BASE = {
    "BLADDER": "(Bladder OR Urothelial OR Urothelium OR Detrusor)",
    "KIDNEY": "(Kidney OR Renal OR Nephron OR Glomerulus)",
    "BRAIN": "(Brain OR Cerebral OR CNS OR Neuron OR Glia)",
    "HEART": "(Heart OR Cardiac OR Myocardium)",
    "LIVER": "(Liver OR Hepatic OR Hepatocyte)"
}

# ================= GEMINI CORE =================

def clean_model_name(model_name: str) -> str:
    return model_name.replace("models/", "")

def call_gemini_json(prompt: str, api_key: str) -> List[str]:
    if not api_key: return []
    headers = {"Content-Type": "application/json"}
    payload = {
        "contents": [{"parts": [{"text": prompt}]}], 
        "generationConfig": {"temperature": 0.1, "response_mime_type": "application/json"}
    }
    for modelo_raw in MODELOS_ATIVOS:
        modelo = clean_model_name(modelo_raw)
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={api_key.strip()}"
            resp = requests.post(url, headers=headers, json=payload, timeout=30)
            if resp.status_code == 200:
                try:
                    text = resp.json()["candidates"][0]["content"]["parts"][0]["text"]
                    clean_text = re.sub(r"```json|```", "", text).strip()
                    parsed = json.loads(clean_text)
                    if isinstance(parsed, list): return [str(x) for x in parsed]
                    if isinstance(parsed, dict): 
                        for v in parsed.values(): 
                            if isinstance(v, list): return [str(x) for x in v]
                except: continue
            elif resp.status_code == 429: break 
        except: continue
    return []

def simple_gemini_text(prompt: str, api_key: str) -> str:
    if not api_key: return None
    headers = {"Content-Type": "application/json"}
    payload = {
        "contents": [{"parts": [{"text": prompt}]}], 
        "generationConfig": {"temperature": 0.3}
    }
    for modelo_raw in MODELOS_ATIVOS:
        modelo = clean_model_name(modelo_raw)
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={api_key.strip()}"
            resp = requests.post(url, headers=headers, json=payload, timeout=15)
            if resp.status_code == 200:
                return resp.json()["candidates"][0]["content"]["parts"][0]["text"]
        except: continue
    return None

# ================= MINERAÇÃO DE ALTA PRECISÃO =================

def ner_extraction_batch(textos_completos: List[str], api_key: str, contexto_alvo: str) -> List[str]:
    if not textos_completos: return []
    
    texto_input = "\n---\n".join(textos_completos[:60]) 
    
    # Prompt focado em QUALQUER área, mas APENAS em agentes químicos/alvos
    prompt = f"""
    ROLE: Elite Pharmacological Data Curator.
    TARGET CONTEXT: {contexto_alvo.upper()}.
    
    TASK: Scan the abstracts. Extract ONLY specific chemical agents and molecular targets.
    
    EXTRACT:
    1. DRUG NAMES (Approved, Experimental, Toxins, or Drugs of Abuse).
    2. CHEMICAL COMPOUNDS (e.g., GYY4137, Nimbolide, Curcumin).
    3. MOLECULAR ALVOS (Receptors, Ion Channels, Enzymes like TRPV1, P2X3, mTOR).
    
    STRICTLY IGNORE:
    - NO Clinical procedures, NO study methodologies, NO general biological processes.
    - NO anatomical parts unless they are part of a drug name.
    
    OUTPUT FORMAT: JSON list of strings only.
    
    INPUT:
    {texto_input}
    """
    return call_gemini_json(prompt, api_key)

@st.cache_data(ttl=3600)
def minerar_pubmed(termo_base: str, email: str, usar_ia: bool = True) -> Dict:
    api_key = st.session_state.get("api_key_usuario", "").strip()
    Entrez.email = email
    
    termo_expandido = MAPA_SINONIMOS_BASE.get(termo_base.upper(), termo_base)
    query = f"({termo_expandido}) AND (drug OR inhibitor OR agonist OR antagonist OR compound OR agent) AND (2020:2026[Date])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=250)
        record = Entrez.read(handle); handle.close()
        
        if not record["IdList"]: return {}

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        
        raw_texts_for_ai = []
        artigos_completos = []
        
        for r in records:
            titulo = r.get('TI', '')
            abstract = r.get('AB', '')
            keywords = ' '.join(r.get('OT', []))
            intro = abstract[:300] if len(abstract) > 300 else abstract
            conclusao = abstract[-300:] if len(abstract) > 300 else ""
            texto_rico = f"TITLE: {titulo}\nCONTEXT: {intro} ... {conclusao}\nKEYWORDS: {keywords}"
            raw_texts_for_ai.append(texto_rico)
            artigos_completos.append({"titulo": titulo, "texto": texto_rico})
        
        entidades = []
        if api_key and usar_ia:
            entidades = ner_extraction_batch(raw_texts_for_ai, api_key, termo_base)
        
        # Regex de segurança para capturar o que a IA possa pular
        if True: 
            texto_full = " ".join([a['texto'] for a in artigos_completos])
            regex_codigos = r'\b[A-Z]{2,4}[- ]?[0-9]{3,6}\b'
            entidades.extend(re.findall(regex_codigos, texto_full))
            sufixos = r'\b[A-Z][a-z]{3,}(?:ine|in|mab|ib|ol|on|one|il|ide|ate|ase|an)\b'
            entidades.extend(re.findall(sufixos, texto_full))
            regex_alvos = r'\b[A-Z0-9-]{3,8}\b'
            candidatos = re.findall(regex_alvos, texto_full)
            candidatos = [c for c in candidatos if (re.search(r'\d', c) or c.endswith("R")) and len(c)>2]
            entidades.extend(candidatos)

        # FILTRAGEM COM A BLACKLIST UNIVERSAL
        entidades_limpas = []
        for e in entidades:
            e = e.strip(".,-;:()[] ")
            if len(e) < 3 or e.isdigit(): continue
            if e.upper() in BLACKLIST_MONSTRO: continue
            if e.lower() in ["with", "from", "after", "during", "high", "low", "using", "treated", "group", "sham"]: continue
            entidades_limpas.append(e)

        counts = Counter(entidades_limpas)
        recorrentes = sorted([e for e, c in counts.items() if c >= 1], key=lambda x: counts[x], reverse=True)
        
        final = []
        seen = set()
        for item in recorrentes:
            if item.lower() not in seen:
                final.append(item); seen.add(item.lower())

        return {
            "termos_indicados": final[:55],
            "counts": counts,
            "total_docs": len(artigos_completos),
            "artigos_originais": artigos_completos
        }
    except Exception: return {}

# ================= WRAPPERS =================
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
        artigos = []
        for r in records:
            artigos.append({
                "Title": r.get("TI", "Sem Título"), 
                "Info_IA": r.get("AB", "Resumo indisponível.")[:1000], 
                "Link": f"https://pubmed.ncbi.nlm.nih.gov/{r.get('PMID', '')}/"
            })
        return artigos
    except: return []

def analisar_abstract_com_ia(titulo: str, dados_curtos: str, api_key: str, lang: str = 'pt') -> str:
    if api_key:
        prompt = f"""
        Você é Doutora em Farmacologia (PhD).
        DEDUZA o mecanismo farmacológico principal.
        TÍTULO: "{titulo}"
        RESUMO: "{dados_curtos}"
        SAÍDA: Órgão - Alvo - Ação (Fármaco/Composto)
        """
        resposta_ia = simple_gemini_text(prompt, api_key)
        if resposta_ia: return resposta_ia.strip()
    return "Análise pendente."

def buscar_todas_noticias(lang_code: str) -> List[Dict]:
    return []
