# backend.py
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
    "gemini-2.5-flash", "gemini-2.0-flash", "gemini-1.5-flash"
]

# 1. LISTA NEGRA (O que deletar)
BLACKLIST_MONSTRO = {
    "ASSOCIATION", "EVALUATION", "UNCOMMON", "INTRAUTERINE", "GERMLINE", 
    "COLLISION", "NOTABLY", "IMPACT", "FOLLOW-UP", "KEY", "FINDINGS",
    "LIMITATIONS", "BACKGROUND", "OBJECTIVE", "METHODS", "RESULTS", "CONCLUSION",
    "ABSTRACT", "INTRODUCTION", "STUDY", "ANALYSIS", "DATA", "STATISTICS",
    "SIGNIFICANT", "DIFFERENCE", "BETWEEN", "AMONG", "WITHIN", "DURING",
    "PREVALENCE", "INCIDENCE", "RISK", "FACTOR", "ROLE", "POTENTIAL", "LEVEL",
    "PATIENT", "PATIENTS", "CLINICAL", "CASE", "REPORT", "REVIEW", "META",
    "SYSTEMATIC", "TREATMENT", "THERAPY", "SURGERY", "DIAGNOSIS", "PROGNOSIS",
    "QUALITY", "LIFE", "QOL", "SURVIVAL", "MORTALITY", "MORBIDITY", "SAFETY",
    "EFFICACY", "OUTCOME", "OUTCOMES", "PLACEBO", "CONTROL", "GROUP", "TOTAL",
    "JANUARY", "FEBRUARY", "MARCH", "APRIL", "MAY", "JUNE", "JULY", "AUGUST",
    "SEPTEMBER", "OCTOBER", "NOVEMBER", "DECEMBER", "2020", "2021", "2022", 
    "2023", "2024", "2025", "HOSPITAL", "UNIVERSITY", "NATIONAL", "INTERNATIONAL",
    "DNA", "RNA", "ATP", "MRI", "CT", "PET", "BCG", "BMI", "WHO", "HIC", "AE", "AES"
}

# 2. STOP WORDS (Palavras comuns do inglês para filtrar o regex solto)
STOP_WORDS = {
    "THE", "AND", "FOR", "NOT", "BUT", "WAS", "ARE", "HAS", "HAVE", "HAD",
    "WITH", "WITHOUT", "FROM", "INTO", "ONTO", "VIA", "ALL", "ANY", "TWO", 
    "ONE", "BOTH", "EACH", "WHICH", "THAT", "THIS", "THESE", "THOSE", "BEEN",
    "WERE", "CAN", "MAY", "SHOULD", "COULD", "USING", "USED", "BASED", "HIGH",
    "LOW", "NEW", "OLD", "NON", "ANTI", "PRO", "PRE", "POST", "SUB", "SUPRA"
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
                except: continue
            elif resp.status_code == 429: break 
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

# ================= MINERAÇÃO ROBUSTA =================

def ner_extraction_batch(titulos_keywords: List[str], api_key: str) -> List[str]:
    if not titulos_keywords: return []
    texto_input = "\n".join([f"- {t}" for t in titulos_keywords[:50]])
    prompt = f"""
    Role: Molecular Pharmacologist.
    Task: Extract MOLECULAR TARGETS (Receptors, Enzymes, Genes) and DRUGS/COMPOUNDS.
    
    IGNORE: "Study", "Analysis", "Patients", "Diagnosis", "Results".
    KEEP: "HSP90", "SGLT2", "Alantolactone", "Cisplatin", "mTOR", "VEGF".
    
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
    # Aumentei para 120 para garantir volume
    query = f"({termo_expandido}) AND (2020:2026[Date])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=120)
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]: return {}

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        
        raw_texts = []
        artigos_completos = []
        
        for r in records:
            titulo = r.get('TI', '')
            keywords = ' '.join(r.get('OT', []))
            texto_limpo = f"{titulo} . {keywords}"
            if len(texto_limpo) > 5:
                raw_texts.append(texto_limpo)
                artigos_completos.append({"titulo": titulo, "texto": texto_limpo})
        
        entidades = []
        
        # 1. IA
        if api_key and usar_ia:
            entidades = ner_extraction_batch(raw_texts, api_key)
        
        # 2. REGEX "REDE DE ARRASTÃO" (Agora pega tudo e filtra depois)
        texto_full = " ".join(raw_texts)
        
        # A: Siglas Técnicas (Letras maiúsculas + números) -> P2X3, HSP90
        regex_codigos = r'\b[A-Z]{2,}[0-9]+[A-Z0-9-]*\b'
        candidatos_codigos = re.findall(regex_codigos, texto_full)
        
        # B: Palavras Capitalizadas (Candidates a nomes próprios/drogas)
        # Pega qualquer palavra que comece com Maiúscula e tenha 4+ letras
        regex_capitalizadas = r'\b[A-Z][a-z]{3,}\b' 
        candidatos_capitalizadas = re.findall(regex_capitalizadas, texto_full)
        
        # C: Acrônimos (3+ Maiúsculas)
        regex_acronimos = r'\b[A-Z]{3,}\b'
        candidatos_acronimos = re.findall(regex_acronimos, texto_full)

        entidades.extend(candidatos_codigos)
        entidades.extend(candidatos_capitalizadas)
        entidades.extend(candidatos_acronimos)

        # 3. FILTRAGEM (Onde a mágica acontece)
        entidades_limpas = []
        for e in entidades:
            e = e.strip(".,-;:()[] ")
            if len(e) < 3: continue 
            if e.isdigit(): continue
            
            e_upper = e.upper()
            
            # Se estiver na Blacklist -> LIXO
            if e_upper in BLACKLIST_MONSTRO: continue
            
            # Se for Stop Word (The, And, With) -> LIXO
            if e_upper in STOP_WORDS: continue
            
            entidades_limpas.append(e)

        # 4. CONTAGEM E RETORNO
        counts = Counter(entidades_limpas)
        
        # Se a IA não rodou, exigimos pelo menos 2 aparições para reduzir ruído
        # MAS, se a lista estiver muito curta, aceitamos 1 aparição (Modo Resgate)
        min_count = 2 if (not api_key) else 1
        if len(counts) < 10: min_count = 1
        
        recorrentes = [e for e, c in counts.items() if c >= min_count]
        recorrentes = sorted(recorrentes, key=lambda x: counts[x], reverse=True)
        
        # Deduplicação Final
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
    except Exception as e:
        print(f"Erro no Backend: {e}") # Para debug no terminal
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
