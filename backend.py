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

# Modelos Flash são essenciais para o plano gratuito (gastam menos quota)
MODELOS_ATIVOS = [
    "gemini-2.5-flash",          # Mais novo e rápido
    "gemini-2.0-flash",          # Muito estável
    "gemini-2.0-flash-lite",     # Ultra econômico
    "gemini-1.5-flash"           # Fallback clássico
]

# A Blacklist continua para garantir que, mesmo sem abstract, não passe lixo
BLACKLIST_MONSTRO = {
    # Metodologia
    "OBJECTIVE", "INTRODUCTION", "METHODS", "RESULTS", "CONCLUSION", "DISCUSSION",
    "BACKGROUND", "ABSTRACT", "HYPOTHESIS", "PURPOSE", "MATERIAL", "MATERIALS",
    "PRESENTATION", "REPORT", "CASE", "SERIES", "REVIEW", "META-ANALYSIS",
    "SYSTEMATIC", "COHORT", "RETROSPECTIVE", "PROSPECTIVE", "RANDOMIZED", "TRIAL",
    "STUDY", "ANALYSIS", "DATA", "STATISTICS", "KEY", "FINDINGS", "LIMITATIONS",
    "REPLY", "LETTER", "EDITOR", "EDITORIAL", "DOI", "PMID", "PMC",
    
    # Tempo e Lugar
    "JANUARY", "FEBRUARY", "MARCH", "APRIL", "MAY", "JUNE", "JULY", "AUGUST", 
    "SEPTEMBER", "OCTOBER", "NOVEMBER", "DECEMBER", "2020", "2021", "2022", "2023", 
    "2024", "2025", "HOSPITAL", "UNIVERSITY", "CENTER", "DEPARTMENT", "CHINA", "USA", 
    "UK", "JAPAN", "INDIAN", "CHINESE", "EUROPEAN", "AMERICAN", "NATIONAL", "INTERNATIONAL",
    
    # Termos Gerais (Não moleculares)
    "PATIENT", "PATIENTS", "ADULT", "CHILD", "QUALITY", "LIFE", "QOL", "SURVIVAL", 
    "RISK", "FACTOR", "SCORE", "DIAGNOSIS", "PROGNOSIS", "THERAPY", "TREATMENT", 
    "SURGERY", "OUTCOME", "EFFICACY", "SAFETY", "ADVERSE", "EVENTS", "SYMPTOMS", 
    "DISEASE", "CANCER", "TUMOR", "BLADDER", "UROTHELIAL", "URINARY", "CYSTOSCOPY",
    "CLINICAL", "PRIMARY", "NORMAL", "CONTROL", "BASELINE", "TOTAL", "HIGH", "LOW",
    "LEVEL", "INCREASED", "DECREASED", "POSITIVE", "NEGATIVE", "SIGNIFICANT", 
    "ASSOCIATED", "BETWEEN", "AMONG", "WITH", "WITHOUT", "DURING", "AFTER", "BEFORE",
    "HUMAN", "MOUSE", "RAT", "ANIMAL", "MODEL", "CELL", "CELLS", "TISSUE", "BLOOD",
    "URINE", "EXPRESSION", "ACTIVITY", "FUNCTION", "MECHANISM", "PATHWAY", "TARGET", 
    "POTENTIAL", "ROLE", "NOVEL", "NEW", "RECENT", "FUTURE", "IMPACT", "EFFECT",
    "SPECIFIC", "NON", "ANTI", "PRO", "TYPE", "GENE", "PROTEIN", "MOLECULE", 
    "RECEPTOR", "ENZYME", "DRUG", "PHARMACOLOGY", "SCIENCE", "MEDICAL", "HEALTH",
    
    # Conjunções/Preposições comuns em títulos
    "THE", "AND", "FOR", "NOT", "BUT", "WAS", "ARE", "HAS", "ITS", "THAT", "THIS",
    "FROM", "INTO", "ONTO", "VIA", "ALL", "ANY", "TWO", "ONE", "BOTH"
}

MAPA_SINONIMOS_BASE = {
    "BLADDER": "(Bladder OR Urothelial OR Urothelium)",
    "PAIN": "(Pain OR Nociception OR Analgesia)",
    "INFLAMMATION": "(Inflammation OR Cytokines OR NF-kappaB)"
}

# ================= GEMINI CORE =================

def clean_model_name(model_name: str) -> str:
    return model_name.replace("models/", "")

def call_gemini_json(prompt: str, api_key: str) -> List[str]:
    """Tenta extrair JSON. Se der erro 429 (Cota excedida), retorna lista vazia para usar o Regex."""
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
                    text = resp.json()["candidates"][0]["content"]["parts"][0]["text"]
                    clean_text = re.sub(r"```json|```", "", text).strip()
                    parsed = json.loads(clean_text)
                    if isinstance(parsed, list): return [str(x) for x in parsed]
                    if isinstance(parsed, dict): 
                        for v in parsed.values(): 
                            if isinstance(v, list): return [str(x) for x in v]
                    return []
                except: continue
            elif resp.status_code == 429:
                # Se der erro de cota, para de tentar chamar a IA e deixa o Regex assumir
                break 
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

# ================= MINERAÇÃO OTIMIZADA (TOKEN SAVER) =================

def ner_extraction_batch(titulos_keywords: List[str], api_key: str) -> List[str]:
    """Usa IA apenas em Títulos e Keywords (poucos tokens)."""
    if not titulos_keywords: return []
    
    # Envia lote menor para não estourar tokens
    texto_input = "\n".join([f"- {t}" for t in titulos_keywords[:50]])
    
    prompt = f"""
    Role: Molecular Pharmacologist.
    Task: Analyze these paper TITLES and KEYWORDS. Extract a JSON list of SPECIFIC MOLECULAR TARGETS and DRUGS.
    
    RULES:
    1. Extract Receptors (e.g., P2X3, TRPV1), Enzymes (e.g., mTOR, COX-2), Genes, and Drugs (e.g., Mirabegron, Cisplatin).
    2. IGNORE general terms (e.g., "Patients", "Study", "Cancer", "Cell").
    3. IGNORE dates and locations.
    
    INPUT:
    {texto_input}
    
    OUTPUT: JSON list of strings.
    """
    return call_gemini_json(prompt, api_key)

@st.cache_data(ttl=3600)
def minerar_pubmed(termo_base: str, email: str, usar_ia: bool = True) -> Dict:
    api_key = st.session_state.get("api_key_usuario", "").strip()
    Entrez.email = email
    
    termo_expandido = MAPA_SINONIMOS_BASE.get(termo_base.upper(), termo_base)
    # Busca 100 artigos (equilíbrio entre volume e velocidade)
    query = f"({termo_expandido}) AND (2020:2026[Date])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=100)
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]: return {}

        # Fetch APENAS Título (TI) e Keywords (OT) - SEM ABSTRACT (AB)
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        
        raw_texts = []
        artigos_completos = []
        
        for r in records:
            titulo = r.get('TI', '')
            keywords = ' '.join(r.get('OT', []))
            
            # Texto ultra leve para a IA processar rápido
            texto_limpo = f"{titulo} . {keywords}"
            
            if len(texto_limpo) > 5:
                raw_texts.append(texto_limpo)
                # Salvamos no objeto apenas o que temos (sem abstract longo)
                artigos_completos.append({"titulo": titulo, "texto": texto_limpo})
        
        entidades = []
        
        # 1. IA Leve (Se possível)
        if api_key and usar_ia:
            entidades = ner_extraction_batch(raw_texts, api_key)
        
        # 2. Regex Avançado (Funciona sem IA ou se a cota acabar)
        # O "Meio Termo": Pega códigos E nomes de drogas sem pegar lixo
        if True: # Executa sempre para complementar a IA ou agir como fallback
            texto_full = " ".join(raw_texts)
            
            # A: Siglas Técnicas (Letras + Números) -> P2X3, HSP90, IL-6
            regex_codigos = r'\b[A-Z]{2,}[0-9]+[A-Z0-9-]*\b'
            candidatos_codigos = re.findall(regex_codigos, texto_full)
            
            # B: Fármacos (PascalCase + Sufixos Comuns) -> Cisplatin, Mirabegron, Solifenacin
            # Sufixos: -in, -ine, -ib, -mab, -ol, -on, -one, -il, -ide, -ate
            regex_farmacos = r'\b[A-Z][a-z]{3,}(?:in|ine|ib|mab|ol|on|one|il|ide|ate)\b'
            candidatos_farmacos = re.findall(regex_farmacos, texto_full)
            
            # C: Acrônimos Puros (3+ letras maiúsculas) -> VEGF, BDNF, EGFR
            regex_acronimos = r'\b[A-Z]{3,}\b'
            candidatos_acronimos = re.findall(regex_acronimos, texto_full)

            # Junta tudo
            entidades.extend(candidatos_codigos)
            entidades.extend(candidatos_farmacos)
            entidades.extend(candidatos_acronimos)

        # 3. Filtragem (Limpeza)
        entidades_limpas = []
        for e in entidades:
            e = e.strip(".,-;:()[] ")
            
            if len(e) < 3: continue 
            if e.isdigit(): continue
            
            # Verifica Blacklist
            if e.upper() in BLACKLIST_MONSTRO: continue
            
            entidades_limpas.append(e)

        # 4. Contagem e Retorno
        counts = Counter(entidades_limpas)
        
        # Se IA rodou, aceita ocorrência única. Se só Regex rodou, idealmente 2+, mas no Free liberamos 1
        recorrentes = sorted([e for e in counts], key=lambda x: counts[x], reverse=True)
        
        # Deduplicação visual (case insensitive)
        final = []
        seen = set()
        for item in recorrentes:
            if item.lower() not in seen:
                final.append(item)
                seen.add(item.lower())

        return {
            "termos_indicados": final[:40],
            "counts": counts,
            "total_docs": len(artigos_completos),
            "artigos_originais": artigos_completos
        }
    except Exception:
        return {}

# ================= WRAPPERS =================

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
    """Busca abstracts apenas para exibição final (leitura humana), não para mineração."""
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
    """Função sob demanda: Só gasta token se o usuário clicar para analisar um paper específico."""
    if not api_key: return "Chave API necessária."
    prompt = f"Analyze: {titulo}\nOutput one line: TARGET | DRUG | EFFECT"
    return simple_gemini_text(prompt, api_key).replace("\n", " ")

def buscar_todas_noticias(lang_code: str) -> List[Dict]:
    Entrez.email = "pesquisador_guest@unifesp.br"
    query = "(molecular pharmacology) AND (bladder) AND (2025:2026[Date])"
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
