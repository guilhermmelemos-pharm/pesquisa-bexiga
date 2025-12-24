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

# --- 3. RADAR DE NOTÍCIAS (RESTAURADO) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_todas_noticias(lang='pt'):
    try:
        query = "(bladder physiology OR pharmacology) AND (2024:2026[Date - Publication])"
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
                if line.startswith("JT  - "): journal = line[3:].strip()
                if line.strip().isdigit() and not pmid: pmid = line.strip()
            if tit and pmid:
                news.append({"titulo": tit, "fonte": journal[:40], "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return news
    except: return []

# --- 4. DOUTORA INVESTIGADORA (FILTRO SEVERO) ---
def _doutora_investigadora(termo_base, lista_pubmed=None, fase="brainstorming"):
    api_key = st.session_state.get('api_key_usuario', '')
    if not api_key: return []

    if fase == "brainstorming":
        prompt = f"Como PhD em farmacologia, liste 50 alvos moleculares (canais, receptores, enzimas) para {termo_base}. Retorne APENAS uma lista Python de strings."
    else:
        # Enviamos apenas os top 100 para não estourar a IA
        lista_str = ", ".join(lista_pubmed[:100])
        prompt = f"""
        Cruze seu conhecimento com estes termos: {lista_str}.
        MISSÃO: Isolar apenas alvos de bancada e fármacos.
        REGRAS: Remova termos genéricos (RNA, DNA, Study, Effect, Results, Body), conectivos e números isolados.
        Retorne APENAS a lista Python final (max 60 itens).
        """

    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.2}}
    
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=20)
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                # Limpeza extra para garantir que o ast.literal_eval funcione
                texto_limpo = re.search(r'\[.*\]', texto, re.DOTALL)
                if texto_limpo:
                    return ast.literal_eval(texto_limpo.group())
        except: continue
    return []

# --- 5. MINERAÇÃO E CRUZAMENTO ---
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
        
        # PENEIRA DE REGEX: Só aceita termos com 3+ letras, ignora números puros e conectivos comuns
        blacklist = {'THE', 'AND', 'FOR', 'WITH', 'FROM', 'ROLE', 'STUDY', 'DATA', 'CELL', 'LEVEL'}
        candidatos_raw = re.findall(r'\b[A-Z][A-Z0-9-]{2,}\b', full_data) # Siglas técnicas de 3+ letras
        
        contagem = Counter([c.upper() for c in candidatos_raw if c.upper() not in blacklist])
        top_pubmed = [termo for termo, freq in contagem.most_common(150)]
        
        nomes_finais = []
        if usar_ia:
            lista_cruzamento = list(set(alvos_previstos + top_pubmed))
            nomes_finais = _doutora_investigadora(termo_base, lista_pubmed=lista_cruzamento, fase="cruzamento")
        else:
            nomes_finais = top_pubmed[:40]

        if not nomes_finais: return [] # Se falhar, o app carrega o preset

        res_finais = []
        for nome in nomes_finais:
            n_limpo = limpar_termo_para_pubmed(nome)
            freq = contagem.get(n_limpo.upper(), 1)
            l_score, p_val, b_ocean, status = calcular_metricas_originais(freq, total_docs, freq * 10)
            res_finais.append({"Alvo": n_limpo, "Lambda": round(l_score, 2), "P-value": round(p_val, 4), "Blue Ocean": round(b_ocean, 1), "Status": status})
        return res_finais
    except: return []

# --- 6. ANÁLISE IA ---
def analisar_abstract_com_ia(titulo, dados, api_key, lang='pt'):
    if not api_key: return "Chave necessária."
    prompt = f"Analise: {titulo}. Resuma Alvo, Fármaco e Efeito Funcional em 20 palavras."
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}]}
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelante.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=12)
            if resp.status_code == 200:
                return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
        except: continue
    return "Erro IA."
