# backend.py
import streamlit as st
from Bio import Entrez, Medline
import requests
import json
import re
import time
import math
from collections import Counter
from difflib import SequenceMatcher
from typing import List, Dict, Any, Union, Optional

# ================= CONFIG =================
Entrez.email = "pesquisador_guest@unifesp.br"
# Recomendado: Mover chaves para st.secrets em produção
MODELOS_ATIVOS = ["gemini-2.0-flash", "gemini-2.0-flash-exp", "gemini-1.5-flash"]

MAPA_SINONIMOS_BASE = {
    "BLADDER": "(Bladder OR Vesical OR Detrusor OR Urothelium)",
    "PAIN": "(Pain OR Nociception OR Analgesia)",
    "INFLAMMATION": "(Inflammation OR Cytokines OR Inflammasome)"
}

# ================= GEMINI CORE =================

def call_gemini_json(prompt: str, api_key: str, temperature: float = 0.1) -> Any:
    """Chamada robusta à API esperando resposta JSON."""
    if not api_key:
        return []

    headers = {"Content-Type": "application/json"}
    # Instrução explícita para JSON no payload se o modelo suportar, 
    # ou embutido no prompt (garantido abaixo)
    payload = {
        "contents": [{"parts": [{"text": prompt}]}], 
        "generationConfig": {
            "temperature": temperature,
            "response_mime_type": "application/json" # Força JSON estruturado
        }
    }
    
    for modelo in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={api_key.strip()}"
            resp = requests.post(url, headers=headers, json=payload, timeout=30)
            
            if resp.status_code == 200:
                result = resp.json()
                try:
                    text_content = result["candidates"][0]["content"]["parts"][0]["text"]
                    return json.loads(text_content)
                except (KeyError, json.JSONDecodeError, IndexError):
                    # Fallback: tentar limpar blocos de código markdown se o JSON falhar
                    clean_text = re.sub(r"```json|```", "", text_content).strip()
                    return json.loads(clean_text)
            
            elif resp.status_code == 429:
                time.sleep(1.5) # Backoff simples
                continue
                
        except Exception as e:
            # Em produção, usar logging.error(e)
            continue
            
    return []

def simple_gemini_text(prompt: str, api_key: str) -> str:
    """Chamada simples para texto livre (resumos)."""
    if not api_key: return ""
    headers = {"Content-Type": "application/json"}
    payload = {
        "contents": [{"parts": [{"text": prompt}]}], 
        "generationConfig": {"temperature": 0.2}
    }
    
    for modelo in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={api_key.strip()}"
            resp = requests.post(url, headers=headers, json=payload, timeout=25)
            if resp.status_code == 200:
                return resp.json()["candidates"][0]["content"]["parts"][0]["text"]
        except: continue
    return ""

# ================= NER (EXTRAÇÃO DE ALVOS) =================

def ner_extraction_batch(artigos: List[Dict], api_key: str) -> List[str]:
    # Limita o contexto para evitar estouro de tokens
    texto_amostra = [f"- {a['texto']}" for a in artigos[:30]] 
    texto_input = "\n".join(texto_amostra)
    
    prompt = f"""
    Role: Senior Molecular Pharmacologist.
    Task: Extract a JSON list of specific MOLECULAR TARGETS (Receptors, Enzymes, Channels) and DRUGS mentioned in the text.
    
    STRICT EXCLUSIONS (Do NOT include):
    - Clinical terms (MRI, CT, PET, BCG, TURBT, GWAS, RBE)
    - Cell lines (T24, UM-UC-3)
    - General terms (DNA, RNA, ATP, Gene, Protein, Patient, Study)
    - Disease/Anatomy (BLCA, NMIBC, OAB, LUTS, Bladder)
    
    Format: Return ONLY a JSON list of strings. Example: ["P2X3", "Mirabegron", "TRPV1"]
    
    TEXT DATA:
    {texto_input}
    """
    
    resultado = call_gemini_json(prompt, api_key)
    if isinstance(resultado, list):
        return [str(x) for x in resultado]
    return []

# ================= ANALISAR ABSTRACT =================

def analisar_abstract_com_ia(titulo: str, dados_curtos: str, api_key: str) -> str:
    if not api_key: return "N/A | N/A | N/A"
    prompt = f"Analyze: {titulo}. Context: {dados_curtos}. Extract: TARGET | DRUG | EFFECT. Output strictly one line separated by pipes."
    res = simple_gemini_text(prompt, api_key)
    if not res or "|" not in res: return "Alvo | Fármaco | Análise em andamento"
    return res.strip()

# ================= MÉTRICAS E PIPELINE =================

def calcular_metricas_originais(freq: int, total_docs: int, n_alvo_total: int):
    if total_docs == 0: return 0, 1.0, 0, "N/A"
    
    lambda_score = (freq / total_docs) * 100
    p_val = math.exp(-freq/2.5) if freq > 0 else 1.0 
    blue_ocean = max(10, 100 - (n_alvo_total * 2))
    
    return lambda_score, round(p_val, 4), round(blue_ocean, 2), "Competitivo"

@st.cache_data(ttl=3600)
def minerar_pubmed(termo_base: str, email: str) -> Dict:
    api_key = st.session_state.get("api_key_usuario", "").strip()
    Entrez.email = email
    
    termo_expandido = MAPA_SINONIMOS_BASE.get(termo_base.upper(), termo_base)
    query = f"({termo_expandido}) AND (2019:2026[Date])"
    
    try:
        # Busca IDs
        handle = Entrez.esearch(db="pubmed", term=query, retmax=150)
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]: return {}

        # Fetch e Parsing robusto com Bio.Medline
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        artigos_parsed = []
        
        for r in records:
            # Combina Título (TI), Palavras-chave (OT/KW) e Abstract (AB) se disponível
            texto_combinado = f"{r.get('TI', '')} {r.get('AB', '')} {' '.join(r.get('OT', []))}"
            if texto_combinado.strip():
                artigos_parsed.append({"titulo": r.get('TI', ''), "texto": texto_combinado})
        
        # 1. Tentativa via IA
        entidades = []
        if api_key:
            entidades = ner_extraction_batch(artigos_parsed, api_key)
        
        # 2. Fallback Determinístico (Se IA falhar ou não houver chave)
        if not entidades:
            texto_full = " ".join([a['texto'] for a in artigos_parsed])
            # Regex melhorado para capturar siglas e nomes compostos simples
            candidatos = re.findall(r'\b[A-Z][a-zA-Z0-9-]{2,15}\b', texto_full)
            
            blacklist_farmaco = {
                "MRI", "CT", "PET", "BCG", "DNA", "RNA", "ATP", "GWAS", "TURP", "NHANES",
                "LUTS", "OAB", "BPH", "IPSS", "AUA", "EAU", "NMIBC", "MIBC", "TURBT",
                "SBRT", "IMPT", "RBE", "BLCA", "T24", "HLA", "CALR", "SGLT2", "VIII",
                "STUDY", "ROLE", "CELL", "TISSUE", "PATIENT", "AND", "WITH", "THE", "FOR",
                "ABSTRACT", "BACKGROUND", "METHODS", "RESULTS", "CONCLUSION"
            }
            entidades = [e for e in candidatos if e.upper() not in blacklist_farmaco]

        # Contagem e Normalização
        counts = Counter(entidades)
        # Filtra ruído de ocorrência única se houver muitos dados
        limit_cut = 1 if len(entidades) < 50 else 2
        recorrentes = [e for e, c in counts.items() if c >= limit_cut]
        
        # Se filtrou demais, pega os top 15
        if not recorrentes: 
            recorrentes = [e for e, _ in counts.most_common(15)]
        
        # Deduplicação Fuzzy simples
        canonicos = []
        for e in recorrentes:
            # Se não for similar a nenhum já adicionado
            if not any(SequenceMatcher(None, e.lower(), c.lower()).ratio() > 0.90 for c in canonicos):
                canonicos.append(e)

        return {
            "termos_indicados": canonicos,
            "counts": counts,
            "total_docs": len(artigos_parsed),
            "artigos_originais": artigos_parsed
        }
    except Exception as e:
        # st.error(f"Erro na mineração: {e}") # Descomentar para debug
        return {}

# ================= FUNÇÕES DE APOIO =================

def buscar_alvos_emergentes_pubmed(alvo: str, email: str) -> List[str]:
    res = minerar_pubmed(alvo, email)
    return res.get("termos_indicados", [])

@st.cache_data(ttl=3600)
def buscar_resumos_detalhados(termo: str, orgao: str, email: str, ano_ini: int, ano_fim: int) -> List[Dict]:
    Entrez.email = email
    query = f"({termo}) AND ({orgao}) AND ({ano_ini}:{ano_fim}[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=10)
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        
        artigos = []
        for r in records:
            tit = r.get("TI", "Sem título")
            abstract = r.get("AB", "Resumo indisponível.")[:800] # Limita tamanho
            pmid = r.get("PMID", "")
            if pmid:
                artigos.append({
                    "Title": tit, 
                    "Info_IA": abstract, 
                    "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                })
        return artigos
    except: return []

def consultar_pubmed_count(termo: str, contexto: str, email: str, ano_ini: int, ano_fim: int) -> int:
    Entrez.email = email
    query = f"({termo}) AND ({contexto}) AND ({ano_ini}:{ano_fim}[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle)
        handle.close()
        return int(record["Count"])
    except:
        return 0

def buscar_todas_noticias(email: str) -> List[Dict]:
    Entrez.email = email
    query = "(molecular bladder pharmacology) AND (2024:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5, sort="pub_date")
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        
        news = []
        for r in records:
            tit = r.get("TI", "Sem título")
            journal = r.get("JT", "Journal desconhecido")
            pmid = r.get("PMID", "")
            if pmid:
                news.append({
                    "titulo": tit, 
                    "fonte": journal[:30], 
                    "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                })
        return news
    except: return []
