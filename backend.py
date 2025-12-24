import streamlit as st
from Bio import Entrez
import requests
import json
import re
from collections import Counter
import ast
import math

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br"

MODELOS_ATIVOS = [
    "gemini-2.0-flash", "gemini-2.0-flash-exp", 
    "gemini-1.5-flash", "gemini-flash-latest"
]

# --- 1. MOTOR DE MÉTRICAS ---
def calcular_metricas_originais(freq, total_docs, n_alvo_total):
    lambda_score = (freq / total_docs) * 100 if total_docs > 0 else 0
    p_value = math.exp(-freq/5) if freq > 0 else 1.0
    blue_ocean = max(0, 100 - (n_alvo_total / 10)) if n_alvo_total > 0 else 100.0
    status = "Saturado" if blue_ocean < 25 else "Blue Ocean" if blue_ocean > 75 else "Competitivo"
    return lambda_score, p_value, blue_ocean, status

# --- 2. HIGIENIZAÇÃO ---
def limpar_termo_para_pubmed(termo):
    return re.sub(r'^[^a-zA-Z0-9]+|[^a-zA-Z0-9]+$', '', termo.strip())

# --- 3. RADAR DE NOTÍCIAS (2025) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_todas_noticias(lang='pt'):
    try:
        query = "(bladder physiology OR pharmacology OR molecular biology) AND (2024:2026[Date - Publication])"
        handle = Entrez.esearch(db="pubmed", term=query, retmax=6, sort="pub_date")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []
        
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        news = []
        for art in dados.split("\n\nPMID-"):
            tit, journal, pmid = "", "", ""
            for line in art.split("\n"):
                if line.startswith("TI  - "): tit = line[6:].strip()
                if line.startswith("JT  - "): journal = line[6:].strip()
                if line.strip().isdigit() and not pmid: pmid = line.strip()
            if tit and pmid:
                news.append({"titulo": tit, "fonte": journal[:40], "link": f"https://pubmed.ncbi.gov/{pmid}/"})
        return news
    except: return []

# --- 4. DOUTORA INVESTIGADORA (STRICT FILTER) ---
def _doutora_investigadora(termo_base, lista_pubmed=None, fase="brainstorming"):
    api_key = st.session_state.get('api_key_usuario', '')
    if not api_key: return []

    if fase == "brainstorming":
        prompt = f"Como PhD em farmacologia molecular, liste 50 alvos moleculares, fármacos e miRNAs para {termo_base}. Ex: TRPV4, GSK1016790A. Apenas lista Python []."
    else:
        lista_str = ", ".join(lista_pubmed[:100])
        prompt = f"""
        Cruze seu conhecimento com estes termos PubMed: {lista_str}. 
        ATENÇÃO: Ignore termos de metadados (PMID, FAU, NLM, JID). 
        Extraia apenas ALVOS (canais, receptores, enzimas) e FÁRMACOS específicos. 
        Apenas lista Python [].
        """

    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.2}}
    
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=15)
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                match = re.search(r'\[.*\]', texto, re.DOTALL)
                if match: return ast.literal_eval(match.group())
        except: continue
    return []

# --- 5. MINERAÇÃO E CRUZAMENTO (O CORAÇÃO DO FIX) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    
    alvos_previstos = []
    if usar_ia: alvos_previstos = _doutora_investigadora(termo_base, fase="brainstorming")
    
    query = f"({termo_base} AND (Pharmacology OR Molecular)) AND (2018:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1000, sort="relevance")
        record = Entrez.read(handle); handle.close()
        total_docs = int(record["Count"])
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        # --- NOVO FILTRO: EXTRAI APENAS TEXTO ÚTIL ---
        blacklist_tags = {'FAU', 'PMID', 'NLM', 'DCOM', 'JID', 'EDAT', 'MHDA', 'CRDT', 'PST', 'HHS', 'NIH', 'NIDDK', 'LID', 'NID', 'PMC', 'PMCR', 'OTO', 'AUID', 'ORCID', 'PHST', 'AID'}
        
        candidatos_pubmed = []
        for artigo in full_data.split("\n\nPMID-"):
            texto_util = ""
            for line in artigo.split("\n"):
                # SÓ OLHA TÍTULOS (TI), ABSTRACT (AB) E KEYWORDS (KW/OT)
                if line.startswith("TI  - ") or line.startswith("AB  - ") or line.startswith("KW  - ") or line.startswith("OT  - "):
                    texto_util += line[6:].strip() + " "
            
            # Captura siglas técnicos e fármacos no texto útil
            encontrados = re.findall(r'\b[A-Z0-9-]{3,}\b', texto_util)
            for t in encontrados:
                if t.upper() not in blacklist_tags:
                    candidatos_pubmed.append(t.upper())

        contagem = Counter(candidatos_pubmed)
        top_pubmed = [termo for termo, freq in contagem.most_common(150)]
        
        nomes_finais = []
        if usar_ia:
            lista_para_cruzamento = list(set(alvos_previstos + top_pubmed))
            nomes_finais = _doutora_investigadora(termo_base, lista_pubmed=lista_para_cruzamento, fase="cruzamento")
        
        if not nomes_finais: nomes_finais = top_pubmed[:60]

        res_finais = []
        for nome in nomes_finais:
            n_limpo = limpar_termo_para_pubmed(nome)
            freq = contagem.get(n_limpo.upper(), 1)
            l_score, p_val, b_ocean, status = calcular_metricas_originais(freq, total_docs, freq * 10)
            res_finais.append({
                "Alvo": n_limpo, "Lambda": round(l_score, 2), "P-value": round(p_val, 4), 
                "Blue Ocean": round(b_ocean, 1), "Status": status
            })
        return res_finais
    except: return []

# --- 6. ANÁLISE IA ---
def analisar_abstract_com_ia(titulo, dados, api_key, lang='pt'):
    if not api_key: return "Chave necessária."
    prompt = f"Analise: {titulo}. Resuma Alvo, Fármaco e Efeito em 20 palavras."
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}]}
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=12)
            if resp.status_code == 200:
                return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
        except: continue
    return "Erro IA."
