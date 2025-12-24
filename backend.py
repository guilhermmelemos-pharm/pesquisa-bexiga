import streamlit as st
from Bio import Entrez
import requests
import json
from tenacity import retry, stop_after_attempt, wait_exponential
import re
from collections import Counter
import ast

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br"

# --- 1. SEUS MODELOS (PAID TIER) ---
MODELOS_ATIVOS = [
    "gemini-2.5-flash", "gemini-2.0-flash", "gemini-2.0-flash-exp", 
    "gemini-flash-latest", "gemini-2.5-flash-lite"
]

# --- 3. A DOUTORA EM FARMACOLOGIA (CURADORIA DE ELITE) ---
def _faxina_ia(lista_suja):
    api_key = st.session_state.get('api_key_usuario', '')
    if not api_key: return lista_suja[:40] 

    lista_str = ", ".join(lista_suja)
    prompt = f"""
    Amiga, tu é uma doutora em farmacologia e fisiologia experiente. 
    Lembre-se: FOCO EM BANCADA E PROMISSOR.
    
    REGRAS:
    1. DELETE: Animais, termos médicos (Diabetes, Incontinence, BPH), e lixo metodológico (Cystometrography, Experimental).
    2. MANTENHA: Alvos (TRPV4, NLRP3, PIEZO), Fármacos (GSK1016790A, Ouabain) e Isoformas/miRNAs.
    3. Tu não é médica, vacilão.
    
    LISTA: {lista_str}
    OUTPUT: Retorne APENAS a lista Python limpa.
    """
    
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.1}}
    
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=12)
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                texto = texto.replace("```python", "").replace("```", "").strip()
                return ast.literal_eval(texto)
        except: continue
    return lista_suja[:35]

# --- 6. MINERAÇÃO MASSIVA COM PRÉ-FILTRO RÍGIDO ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    query = f"({termo_base} AND (Physiology OR Pharmacology OR Molecular)) AND (2018:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1500, sort="relevance")
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        # FILTRO DE EXTERMÍNIO IMEDIATO (Para não poluir os resultados se a IA falhar)
        blacklist_hard = {
            "THE", "AND", "ROLE", "EFFECT", "STUDIES", "STUDY", "RESULTS", "DURING", "AFTER",
            "PRELIMINARY", "EXPERIMENTAL", "BIOLOGICAL", "TECHNIQUE", "MEASUREMENT", "RESPONSES",
            "RESPONSE", "PROPERTIES", "MEANS", "DIABETES", "INCONTINENCE", "REFLUX", "BPH", 
            "CYSTOMETROGRAPHY", "EMG", "FEMALE", "SLEEP", "HEALING", "SITE", "CLOSURE"
        }
        
        candidatos = []
        for artigo in full_data.split("\n\nPMID-"):
            texto = ""
            for line in artigo.split("\n"):
                if line.startswith("TI  - ") or line.startswith("KW  - "): texto += line[6:].strip() + " "
            
            encontrados = re.findall(r'\b(?:[A-Z0-9-]{2,}|[a-z]{1,2}-[A-Z0-9]{2,}|[a-z]{1,2}[A-Z][a-zA-Z0-9-]*)\b', texto)
            for t in encontrados:
                t_clean = t.strip("-").upper()
                # Só passa se não estiver na blacklist e tiver relevância molecular
                if len(t_clean) >= 3 and t_clean not in blacklist_hard:
                    candidatos.append(t_clean)

        contagem = Counter(candidatos)
        top = [termo for termo,freq in contagem.most_common(150)]
        
        if usar_ia and st.session_state.get('api_key_usuario'):
            return _faxina_ia(top)
        else:
            return top[:30] # Se a IA não estiver ativa, entrega o top 30 já sem o lixo da blacklist
    except: return []
