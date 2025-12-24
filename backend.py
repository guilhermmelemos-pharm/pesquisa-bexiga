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

MODELOS_ATIVOS = [
    "gemini-2.5-flash", "gemini-2.0-flash", "gemini-2.0-flash-exp", 
    "gemini-flash-latest", "gemini-2.5-flash-lite"
]

# --- A DOUTORA INVESTIGADORA (PENSAMENTO PRÉVIO E CRUZAMENTO) ---
def _doutora_investigadora(termo_base, lista_pubmed=None, fase="brainstorming"):
    api_key = st.session_state.get('api_key_usuario', '')
    if not api_key: return []

    if fase == "brainstorming":
        prompt = f"""
        Tu é uma doutora em farmacologia molecular. 
        MISSÃO: Gerar uma lista de 50 alvos, fármacos, micro-RNAs e sinalizadores celulares que são a "fronteira" da ciência para {termo_base}.
        FOCO: Bancada, banho de órgãos e biologia molecular. 
        EXEMPLOS DE PADRÃO: TRPV4, GSK1016790A, SPHK1, NLRP3, miR-132, SNARE, PIEZO1.
        OUTPUT: Retorne APENAS uma lista Python de strings.
        """
    else:
        lista_str = ", ".join(lista_pubmed)
        prompt = f"""
        Tu é a doutora investigadora. 
        CRUZAMENTO DE DADOS:
        1. Meus alvos previstos (brainstorming).
        2. Evidência real do PubMed: {lista_str}.
        
        MISSÃO: Criar a lista final (60-100 itens). 
        Priorize siglas técnicas e fármacos específicos. 
        Delete lixo clínico (Diabetes, Incontinence, BPH), lixo acadêmico (Role, Effect, Study) e animais (Toad, Rabbit).
        OUTPUT: Retorne APENAS a lista Python final.
        """

    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.5 if fase=="brainstorming" else 0.2}}
    
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=20)
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                return ast.literal_eval(texto.replace("```python", "").replace("```", "").strip())
        except: continue
    return []

# --- MINERAÇÃO MASSIVA COM FLUXO TRIPLO ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    
    # PASSO 1: IA PENSA PRIMEIRO (Brainstorming)
    alvos_previstos = []
    if usar_ia:
        alvos_previstos = _doutora_investigadora(termo_base, fase="brainstorming")
    
    # PASSO 2: BUSCA REAL NO PUBMED
    query = f"({termo_base} AND (Pharmacology OR Molecular)) AND (2018:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1500)
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        blacklist_hard = {"THE", "AND", "WITH", "ROLE", "EFFECT", "STUDY", "RESULTS", "DURING", "TOAD", "RABBIT"}
        
        candidatos_pubmed = []
        for artigo in full_data.split("\n\nPMID-"):
            texto = ""
            for line in artigo.split("\n"):
                if line.startswith("TI  - ") or line.startswith("KW  - "): texto += line[6:].strip() + " "
            encontrados = re.findall(r'\b(?:[A-Z0-9-]{2,}|[a-z]{1,2}-[A-Z0-9]{2,}|[a-z]{1,2}[A-Z][a-zA-Z0-9-]*)\b', texto)
            candidatos_pubmed.extend([t.upper() for t in encontrados if t.upper() not in blacklist_hard])

        contagem = Counter(candidatos_pubmed)
        top_pubmed = [termo for termo, freq in contagem.most_common(200)]
        
        # PASSO 3: CRUZAMENTO FINAL
        if usar_ia:
            # Mistura o brainstorming com os 200 do PubMed para a IA validar
            lista_para_cruzamento = list(set(alvos_previstos + top_pubmed))
            return _doutora_investigadora(termo_base, lista_pubmed=lista_para_cruzamento, fase="cruzamento")
        else:
            return top_pubmed[:40]
    except: return []
