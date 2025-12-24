import streamlit as st
from Bio import Entrez
import requests
import json
import re
from collections import Counter
import ast

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br"

# --- MODELOS DA SUA LISTA (SÓ 2.5 E 3.0) ---
MODELOS_ATIVOS = [
    "gemini-2.5-flash",
    "gemini-3-flash-preview",
    "gemini-2.0-flash"
]

# --- 1. BLACKLIST DE EXTERMÍNIO (Mata o lixo clínico no Python) ---
BLACKLIST_RADICAL = {
    "THE", "AND", "WITH", "FOR", "FROM", "CASE", "REPORT", "STUDY", 
    "TURBT", "BCG", "NMIBC", "MIBC", "LUTS", "OAB", "UTI", "BPH", 
    "RADS", "SEER", "MRI", "PET", "BBN", "NEOADJUVANT", "PRIMARY",
    "MALIGNANT", "RADICAL", "SURGERY", "CLINICAL", "EFFICACY", "DIAGNOSIS",
    "MANAGEMENT", "UPDATE", "CURRENT", "SYNDROME", "OUTCOMES", "WOMEN"
}

# --- 2. FAXINEIRO IA (USANDO 2.5/3.0) ---
def _faxina_ia(lista_suja):
    api_key = st.session_state.get('api_key_usuario', '').strip()
    # Se não tem chave, faz uma limpeza básica via Python
    if not api_key: 
        return [t for t in lista_suja if t.upper() not in BLACKLIST_RADICAL][:30]

    prompt = f"""
    ACT AS: Senior PhD Pharmacologist.
    TASK: Extract ONLY molecular targets (receptors, channels, enzymes) and experimental compounds.
    
    STRICT RULES:
    - KEEP: TRPV1, TRPA1, NLRP3, STING, BDNF, PDE, Solifenacin, Cyclophosphamide.
    - DELETE: Clinical acronyms (BCG, OAB, LUTS), Imaging (MRI, PET), Procedures (TURBT), and filler words (THE).
    - If it's not a molecule for an 'Organ Bath' or 'Western Blot' experiment, DELETE it.
    
    INPUT: {", ".join(lista_suja)}
    
    OUTPUT: Return strictly a Python list of strings.
    """
    
    payload = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.0}}
    
    for m in MODELOS_ATIVOS:
        try:
            # URL pura e direta para evitar erros de biblioteca
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, json=payload, timeout=12)
            
            if resp.status_code == 200:
                raw_text = resp.json()['candidates'][0]['content']['parts'][0]['text']
                clean_text = re.sub(r'```[a-zA-Z]*', '', raw_text).replace('```', '').strip()
                res = ast.literal_eval(clean_text)
                if isinstance(res, list): return res
        except: continue
    
    # FALLBACK: Se a IA falhar (429 ou 404), o Python limpa o grosso do lixo
    return [t for t in lista_suja if t.upper() not in BLACKLIST_RADICAL][:30]

# --- 3. MINERAÇÃO PUBMED ---
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    # Foco em artigos muito recentes de ciência básica
    query = f"({termo_base}) AND (2023:2026[Date - Publication]) NOT Review[pt]"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1200, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        lines = handle.read().splitlines(); handle.close()
        
        candidatos = []
        for line in lines:
            if line.startswith("TI  - ") or line.startswith("OT  - "):
                # Pega siglas (ex: TRPV4) e nomes químicos que terminam em in, ol, receptor, etc.
                found = re.findall(r'\b(?:[A-Z0-9-]{3,10}|[A-Z][a-z]{3,}(?:in|ol|ide|ase|ant|receptor|channel))\b', line)
                for f in found:
                    if f.upper() not in BLACKLIST_RADICAL:
                        candidatos.append(f)
        
        top_terms = [t for t, count in Counter(candidatos).most_common(120)]
        
        if usar_ia:
            return _faxina_ia(top_terms)
        return top_terms[:30]
    except: return []

# --- 4. ANALISAR ARTIGO ---
def analisar_abstract_com_ia(titulo, resumo_texto, api_key):
    prompt = f"Como Farmacologista, resuma em 15 palavras o Alvo Molecular e efeito tecidual: {titulo}. Contexto: {resumo_texto[:600]}"
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, json={"contents": [{"parts": [{"text": prompt}]}]}, timeout=8)
            if resp.status_code == 200:
                return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
        except: continue
    return "Erro nos modelos 2.5/3.0. Verifique sua cota."

# --- FUNÇÕES SUPORTE (COUNT/RESUMOS) ---
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    query = f"({termo}) AND ({contexto}) AND ({ano_ini}:{ano_fim}[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        return int(Entrez.read(handle)["Count"])
    except: return 0

def buscar_resumos_detalhados(termo, orgao, email, ano_ini, ano_fim):
    query = f"({termo}) AND ({orgao}) AND ({ano_ini}:{ano_fim}[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5)
        ids = Entrez.read(handle)["IdList"]
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
        dados = handle.read().split("\n\nPMID-")
        artigos = []
        for art in dados:
            tit = re.search(r"TI  - (.*)", art)
            ab = re.search(r"AB  - (.*)", art)
            if tit: artigos.append({"Title": tit.group(1), "Info_IA": ab.group(1) if ab else "", "Link": "https://pubmed.ncbi.nlm.nih.gov/"})
        return artigos
    except: return []

def buscar_todas_noticias(lang='pt'): return []
