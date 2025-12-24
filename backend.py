import streamlit as st
from Bio import Entrez
import requests
import json
import re
from collections import Counter
import ast

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br"

# --- MODELOS QUE VOCÊ CONFIRMOU POSSUIR ---
MODELOS_ATIVOS = [
    "gemini-3-flash-preview",
    "gemini-2.5-flash",
    "gemini-2.0-flash"
]

# --- 1. BLACKLIST DE EXTERMÍNIO (Termos que nem o Regex deixa passar) ---
BLACKLIST_RADICAL = {
    "WITH", "AFTER", "SINGLE", "REPORT", "STUDY", "CASE", "MANAGEMENT", "PRIMARY",
    "RADICAL", "POTENTIAL", "LONG", "AFTER", "UNDER", "BETWEEN", "TREATMENT",
    "CLINICAL", "SURGERY", "ROBOTIC", "DIAGNOSIS", "OUTCOMES", "EFFICACY",
    "BLADDER", "URINARY", "DETRUSOR", "UROTHELIUM", "KIDNEY", "PROSTATE"
}

# --- 2. FAXINEIRO IA (PROMPT DE DOUTORADO) ---
def _faxina_ia(lista_suja):
    api_key = st.session_state.get('api_key_usuario', '').strip()
    if not api_key: return lista_suja[:30]

    # Só envia o que parece minimamente promissor para economizar cota e tempo
    prompt = f"""
    ACT AS: Senior PhD Pharmacologist.
    TASK: Extract ONLY specific molecular targets and drugs from this list.
    STRICT RULES:
    - KEEP: Receptors (M3, P2X, TRPV1), Enzymes (ROCK, mTOR), Drugs (Mirabegron, Cyclophosphamide), Ions (Ca2+, K+).
    - DELETE: Anatomy, Clinical words, Verbs, Prepositions, General Bio terms.
    - If it's not a molecule for an 'Organ Bath' experiment, DELETE it.
    
    INPUT: {", ".join(lista_suja)}
    
    OUTPUT: Return strictly a Python list of strings.
    """
    
    payload = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.0}}
    
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, json=payload, timeout=12)
            if resp.status_code == 200:
                raw_text = resp.json()['candidates'][0]['content']['parts'][0]['text']
                clean_list = re.sub(r'```[a-z]*', '', raw_text).replace('```', '').strip()
                return ast.literal_eval(clean_list)
        except: continue
    return lista_suja[:30]

# --- 3. MINERAÇÃO PUBMED (FILTRO CIRÚRGICO) ---
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    query = f"({termo_base}) AND (2021:2026[Date - Publication]) NOT Review[pt]"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1200, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        lines = handle.read().splitlines(); handle.close()
        
        candidatos = []
        for line in lines:
            if line.startswith("TI  - ") or line.startswith("OT  - "):
                # REGEX MELHORADO: 
                # Pega siglas (3+ letras) OU palavras que terminam com sufixos farmacológicos
                found = re.findall(r'\b(?:[A-Z]{3,}[A-Z0-9-]*|[A-Z][a-z]{3,}(?:in|one|ol|ide|ase|ant|receptor|channel))\b', line)
                for f in found:
                    f_up = f.upper()
                    if f_up not in BLACKLIST_RADICAL and len(f) > 2:
                        candidatos.append(f)
        
        # Filtro de abundância: se o termo aparece pouco, pode ser sujeira
        top_terms = [t for t, count in Counter(candidatos).most_common(100)]
        
        if usar_ia:
            return _faxina_ia(top_terms)
        return top_terms[:40]
    except: return []

# --- 4. ANALISAR ARTIGO (FOCO EM BANCADA) ---
def analisar_abstract_com_ia(titulo, resumo_texto, api_key):
    prompt = f"Extract molecular target and signaling pathway (max 15 words) from: {titulo}. Context: {resumo_texto[:600]}"
    
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, json={"contents": [{"parts": [{"text": prompt}]}]}, timeout=8)
            if resp.status_code == 200:
                return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
        except: continue
    return "Erro nos modelos 2.5/3.0. Verifique sua cota."

# --- FUNÇÕES SUPORTE ---
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
