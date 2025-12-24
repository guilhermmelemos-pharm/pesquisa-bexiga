# backend.py
import streamlit as st
from Bio import Entrez, Medline
import requests
import json
import re
import time
from collections import Counter
from difflib import SequenceMatcher
from typing import List, Dict, Any

# ================= CONFIGURAÇÃO =================
Entrez.email = "pesquisador_guest@unifesp.br"

# Modelos rápidos e inteligentes
MODELOS_ATIVOS = [
    "gemini-2.5-flash",
    "gemini-2.0-flash",
    "gemini-1.5-flash",
    "gemini-2.0-flash-lite"
]

# Blacklist agressiva para limpar termos médicos/metodológicos que não são alvos
BLACKLIST_MONSTRO = {
    # Metodologia e Estrutura do Artigo
    "OBJECTIVE", "INTRODUCTION", "METHODS", "RESULTS", "CONCLUSION", "DISCUSSION",
    "BACKGROUND", "ABSTRACT", "HYPOTHESIS", "PURPOSE", "MATERIAL", "MATERIALS",
    "PRESENTATION", "REPORT", "CASE", "SERIES", "REVIEW", "META-ANALYSIS",
    "SYSTEMATIC", "COHORT", "RETROSPECTIVE", "PROSPECTIVE", "RANDOMIZED", "TRIAL",
    "STUDY", "ANALYSIS", "DATA", "STATISTICS", "KEY", "FINDINGS", "LIMITATIONS",
    "REPLY", "LETTER", "EDITOR", "EDITORIAL", "COMMENT", "NOTE", "FIGURE", "TABLE",
    "FIG", "OPEN", "ACCESS", "SUPPLEMENTARY", "APPENDIX", "FILES", "PREPRINT",
    
    # Datas e Locais
    "JANUARY", "FEBRUARY", "MARCH", "APRIL", "MAY", "JUNE", "JULY", "AUGUST", 
    "SEPTEMBER", "OCTOBER", "NOVEMBER", "DECEMBER", "2019", "2020", "2021", "2022", 
    "2023", "2024", "2025", "HOSPITAL", "UNIVERSITY", "CENTER", "DEPARTMENT", 
    "CHINA", "USA", "UK", "JAPAN", "INDIAN", "CHINESE", "EUROPEAN", "AMERICAN",
    "ASIAN", "AFRICA", "AUSTRALIA", "ZEALAND", "NATIONAL", "INTERNATIONAL",
    
    # Termos Clínicos Gerais (Não são alvos moleculares)
    "PATIENT", "PATIENTS", "GROUP", "GROUPS", "ADULT", "CHILD", "MALE", "FEMALE",
    "AGE", "YEARS", "OLD", "QUALITY", "LIFE", "QOL", "SURVIVAL", "RISK", "FACTOR",
    "FACTORS", "SCORE", "INDEX", "GRADE", "STAGE", "DIAGNOSIS", "PROGNOSIS",
    "THERAPY", "TREATMENT", "MANAGEMENT", "SURGERY", "SURGICAL", "RESECTION",
    "OUTCOME", "OUTCOMES", "EFFICACY", "SAFETY", "ADVERSE", "EVENTS", "AES",
    "SYMPTOMS", "SYNDROME", "DISEASE", "DISORDER", "CANCER", "CARCINOMA", "TUMOR",
    "METASTASIS", "INVASIVE", "BLADDER", "UROTHELIAL", "URINARY", "TRACT",
    "CYSTOSCOPY", "ULTRASOUND", "IMAGING", "MRI", "CT", "PET", "BIOPSY",
    "CLINICAL", "PRECLINICAL", "PRIMARY", "SECONDARY", "TERTIARY", "NORMAL",
    "CONTROL", "PLACEBO", "BASELINE", "FOLLOW-UP", "TOTAL", "RATE", "RATIO",
    "MEAN", "MEDIAN", "RANGE", "HIGH", "LOW", "LEVEL", "INCREASED", "DECREASED",
    "POSITIVE", "NEGATIVE", "SIGNIFICANT", "ASSOCIATED", "BETWEEN", "AMONG",
    "WITH", "WITHOUT", "DURING", "AFTER", "BEFORE", "WHILE", "WHEN", "WHERE",
    "THIS", "THAT", "THESE", "THOSE", "THEIR", "FROM", "BASED", "USING", "USED",
    "HUMAN", "MOUSE", "RAT", "ANIMAL", "MODEL", "CELL", "CELLS", "LINE", "LINES",
    "TISSUE", "BLOOD", "URINE", "SERUM", "PLASMA", "EXPRESSION", "ACTIVITY",
    "FUNCTION", "MECHANISM", "PATHWAY", "TARGET", "TARGETS", "POTENTIAL", "ROLE",
    "NOVEL", "NEW", "RECENT", "CURRENT", "FUTURE", "IMPACT", "EFFECT", "EFFECTS",
    "CHANGES", "RESPONSE", "RESISTANCE", "SENSITIVITY", "INHIBITION", "ACTIVATION",
    "REGULATION", "SIGNALING", "MEDIATED", "DEPENDENT", "INDEPENDENT", "INDUCED",
    "RELATED", "SPECIFIC", "NON", "ANTI", "PRO", "PRE", "POST", "TRANS", "CIS",
    "TYPE", "TYPES", "CLASS", "FAMILY", "GENE", "PROTEIN", "MOLECULE", "RECEPTOR",
    "ENZYME", "KINASE", "LIGAND", "CHANNEL", "TRANSPORTER", "INHIBITOR", "AGONIST",
    "ANTAGONIST", "BLOCKER", "AGENT", "DRUG", "COMPOUND", "MEDICINE", "PHARMACY",
    "PHARMACOLOGY", "CHEMISTRY", "BIOLOGY", "SCIENCE", "MEDICAL", "HEALTH", "CARE",
    "SYSTEM", "NETWORK", "DATABASE", "WEB", "TOOL", "METHOD", "ALGORITHM", "AI",
    "ARTIFICIAL", "INTELLIGENCE", "MACHINE", "LEARNING", "DEEP", "NEURAL",
    
    # Acrônimos de Doenças/Procedimentos (Não são alvos)
    "BCG", "TURBT", "NMIBC", "MIBC", "UTI", "OAB", "BPH", "LUTS", "BMI", "ECOG",
    "ASA", "TNM", "WHO", "AUA", "EAU", "NCCN", "FDA", "EMA", "DNA", "RNA", "MRNA",
    "MIRNA", "LNCRNA", "CIRCRNA", "SIRNA", "SHRNA", "CRISPR", "CAS9", "PCR",
    "QPCR", "RT-PCR", "WESTERN", "BLOT", "ELISA", "IHC", "IF", "H&E", "FACS"
}

# Mapeamento para expandir a busca inicial
MAPA_SINONIMOS_BASE = {
    "BLADDER": "(Bladder OR Urothelial OR Urothelium)",
    "PAIN": "(Pain OR Nociception OR Analgesia)",
    "INFLAMMATION": "(Inflammation OR Cytokines OR NF-kappaB)"
}

# ================= GEMINI CORE =================

def clean_model_name(model_name: str) -> str:
    return model_name.replace("models/", "")

def call_gemini_json(prompt: str, api_key: str) -> List[str]:
    if not api_key: return []
    headers = {"Content-Type": "application/json"}
    payload = {
        "contents": [{"parts": [{"text": prompt}]}], 
        "generationConfig": {
            "temperature": 0.1,
            "response_mime_type": "application/json"
        }
    }
    
    for modelo_raw in MODELOS_ATIVOS:
        modelo = clean_model_name(modelo_raw)
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={api_key.strip()}"
            resp = requests.post(url, headers=headers, json=payload, timeout=20)
            
            if resp.status_code == 200:
                try:
                    candidates = resp.json().get("candidates", [])
                    if not candidates: continue
                    text = candidates[0]["content"]["parts"][0]["text"]
                    clean_text = re.sub(r"```json|```", "", text).strip()
                    parsed = json.loads(clean_text)
                    
                    if isinstance(parsed, list): return [str(x) for x in parsed]
                    elif isinstance(parsed, dict):
                        for val in parsed.values():
                            if isinstance(val, list): return [str(x) for x in val]
                    return []
                except: continue
        except: continue
    return []

def simple_gemini_text(prompt: str, api_key: str) -> str:
    if not api_key: return ""
    headers = {"Content-Type": "application/json"}
    payload = {
        "contents": [{"parts": [{"text": prompt}]}], 
        "generationConfig": {"temperature": 0.2}
    }
    for modelo_raw in MODELOS_ATIVOS:
        modelo = clean_model_name(modelo_raw)
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={api_key.strip()}"
            resp = requests.post(url, headers=headers, json=payload, timeout=20)
            if resp.status_code == 200:
                return resp.json()["candidates"][0]["content"]["parts"][0]["text"]
        except: continue
    return ""

# ================= MINERAÇÃO CIRÚRGICA =================

def ner_extraction_batch(titulos_keywords: List[str], api_key: str) -> List[str]:
    """Extrai entidades farmacológicas apenas de Títulos e Keywords."""
    if not titulos_keywords: return []
    
    # Formata lista simples para o prompt
    texto_input = "\n".join([f"- {t}" for t in titulos_keywords[:40]])
    
    prompt = f"""
    You are a Molecular Pharmacologist. Analyze these scientific paper TITLES and KEYWORDS.
    
    TASK: Extract strictly MOLECULAR TARGETS (Receptors, Channels, Enzymes, Transporters) and SPECIFIC DRUGS/COMPOUNDS.
    
    STRICT RULES:
    1. IGNORE all clinical terms (e.g., "Quality of Life", "Surgery", "Diagnosis", "Bladder Cancer", "Patients").
    2. IGNORE study types (e.g., "Review", "Meta-analysis", "Case Report").
    3. IGNORE dates, locations, and hospital names.
    4. IGNORE general biology terms (e.g., "Cell", "Tissue", "Expression", "Apoptosis") unless part of a specific pathway name.
    5. OUTPUT format: A simple JSON list of strings.
    
    EXAMPLES OF WHAT TO KEEP:
    - "P2X3", "TRPV1", "Muscarinic Receptors", "M3 receptor"
    - "Mirabegron", "Botulinum toxin", "Resiniferatoxin", "Gemcitabine"
    - "NF-kappaB", "mTOR", "VEGF", "COX-2"
    
    EXAMPLES OF WHAT TO THROW AWAY:
    - "Objective", "Conclusion", "January", "Hospital", "Study", "Safety"
    - "Bladder", "Urothelium", "Cancer", "Pain", "Infection"
    
    INPUT DATA:
    {texto_input}
    """
    
    return call_gemini_json(prompt, api_key)

@st.cache_data(ttl=3600)
def minerar_pubmed(termo_base: str, email: str, usar_ia: bool = True) -> Dict:
    api_key = st.session_state.get("api_key_usuario", "").strip()
    Entrez.email = email
    
    termo_expandido = MAPA_SINONIMOS_BASE.get(termo_base.upper(), termo_base)
    # Busca focada nos últimos anos
    query = f"({termo_expandido}) AND (2020:2026[Date])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=100)
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]: return {}

        # Fetch APENAS Title (TI) e Keywords (OT/MH) - SEM ABSTRACT
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        
        raw_texts = []
        artigos_completos = []
        
        for r in records:
            # Pega Título e Keywords (Other Terms)
            titulo = r.get('TI', '')
            keywords = ' '.join(r.get('OT', [])) # Keywords do autor costumam ser ricas
            texto_limpo = f"{titulo} . {keywords}"
            
            if len(texto_limpo) > 10:
                raw_texts.append(texto_limpo)
                artigos_completos.append({"titulo": titulo, "texto": texto_limpo}) # Mantendo 'texto' p/ compatibilidade
        
        entidades = []
        
        # 1. Tenta IA primeiro (Se tiver chave e usar_ia=True)
        if api_key and usar_ia:
            entidades = ner_extraction_batch(raw_texts, api_key)
        
        # 2. Fallback Regex SUPER ESTRITO (se IA falhar ou estiver desligada)
        # Só pega siglas com números (ex: P2X3, TRPV1) ou sufixos farmacológicos óbvios
        if not entidades:
            texto_full = " ".join(raw_texts)
            
            # Padrão 1: Siglas de genes/receptores (ex: P2X3, CD44, EGFR, IL-6)
            # Letras maiúsculas seguidas de números ou hifens
            regex_receptores = r'\b[A-Z]{2,}[0-9]+[A-Z0-9-]*\b'
            candidatos_receptores = re.findall(regex_receptores, texto_full)
            
            # Padrão 2: Fármacos comuns (terminados em in, ol, ib, mab - arriscado mas melhor que pegar "January")
            # regex_farmacos = r'\b[A-Z][a-z]{4,}(?:ine|in|ol|an|ib|mab)\b'
            # candidatos_farmacos = re.findall(regex_farmacos, texto_full)
            
            candidatos = candidatos_receptores # + candidatos_farmacos
            entidades = candidatos

        # 3. Limpeza Final (O Grande Filtro)
        entidades_limpas = []
        for e in entidades:
            e = e.strip()
            # Remove caracteres estranhos
            e = e.strip(".,-;:()[]")
            
            if len(e) < 3: continue # Remove siglas muito curtas (ex: "II", "A")
            
            # Verifica contra a Blacklist Monstro
            if e.upper() in BLACKLIST_MONSTRO: continue
            
            # Verifica se é apenas um ano
            if e.isdigit() or (len(e)==4 and e.startswith("20")): continue
            
            entidades_limpas.append(e)

        # 4. Estatística
        counts = Counter(entidades_limpas)
        
        # Se IA foi usada, confiamos mais. Se foi regex, cortamos coisas raras.
        corte = 1 if (api_key and usar_ia) else 2
        recorrentes = [e for e, c in counts.items() if c >= corte]
        
        # Ordena por frequência
        recorrentes = sorted(recorrentes, key=lambda x: counts[x], reverse=True)

        return {
            "termos_indicados": recorrentes[:30], # Top 30 para não poluir
            "counts": counts,
            "total_docs": len(artigos_completos),
            "artigos_originais": artigos_completos
        }
    except Exception:
        return {}

# ================= WRAPPERS P/ FRONTEND (MANTIDOS IGUAIS) =================

def buscar_alvos_emergentes_pubmed(alvo: str, email: str, usar_ia: bool = True) -> List[str]:
    res = minerar_pubmed(alvo, email, usar_ia=usar_ia)
    return res.get("termos_indicados", [])

def consultar_pubmed_count(termo: str, contexto: str, email: str, ano_ini: int, ano_fim: int) -> int:
    Entrez.email = email
    query = f"({termo}) AND ({contexto}) AND ({ano_ini}:{ano_fim}[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle)
        handle.close()
        return int(record["Count"])
    except: return 0

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
        artigos = []
        for r in records:
            artigos.append({
                "Title": r.get("TI", "Sem Título"), 
                "Info_IA": r.get("AB", "Resumo indisponível.")[:800], 
                "Link": f"https://pubmed.ncbi.nlm.nih.gov/{r.get('PMID', '')}/"
            })
        return artigos
    except: return []

def analisar_abstract_com_ia(titulo: str, dados_curtos: str, api_key: str, lang: str = 'pt') -> str:
    if not api_key: return "Chave API necessária."
    prompt = f"Analyze: {titulo}\nContext: {dados_curtos}\nOutput one line: TARGET | DRUG | EFFECT"
    return simple_gemini_text(prompt, api_key).replace("\n", " ")

def buscar_todas_noticias(lang_code: str) -> List[Dict]:
    Entrez.email = "pesquisador_guest@unifesp.br"
    query = "(molecular pharmacology) AND (bladder) AND (2024:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=4, sort="pub_date")
        record = Entrez.read(handle)
        handle.close()
        if not record["IdList"]: return []
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        news = []
        for r in records:
            news.append({
                "titulo": r.get("TI", "Novo Artigo"), 
                "fonte": r.get("JT", "Journal")[:20], 
                "link": f"https://pubmed.ncbi.nlm.nih.gov/{r.get('PMID', '')}/"
            })
        return news
    except: return []
