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

# --- 3. RADAR DE NOTÍCIAS (RESTAURADO PARA 2025) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_todas_noticias(lang='pt'):
    try:
        # Busca focada em fisiologia e farmacologia molecular recente
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
                if line.startswith("JT  - "): journal = line[3:].strip()
                if line.strip().isdigit() and not pmid: pmid = line.strip()
            if tit and pmid:
                news.append({
                    "titulo": tit, 
                    "fonte": journal[:40], 
                    "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                })
        return news
    except: return []

# --- 4. DOUTORA INVESTIGADORA (FILTRO SEVERO E CRUZAMENTO) ---
def _doutora_investigadora(termo_base, lista_pubmed=None, fase="brainstorming"):
    api_key = st.session_state.get('api_key_usuario', '')
    if not api_key: return []

    if fase == "brainstorming":
        # PASSO 1: Pensamento prévio baseado em conhecimento de PhD
        prompt = f"""
        Tu é uma doutora em farmacologia molecular. 
        MISSÃO: Listar 50 alvos moleculares (canais, receptores, enzimas), micro-RNAs e fármacos que são a "fronteira" para {termo_base}.
        FOCO: Bancada, banho de órgãos e biologia molecular. 
        EXEMPLOS: TRPV4, GSK1016790A, SPHK1, NLRP3, miR-132, SNARE, PIEZO1.
        OUTPUT: Retorne APENAS uma lista Python de strings.
        """
    else:
        # PASSO 3: Cruzamento final com dados reais
        lista_str = ", ".join(lista_pubmed[:100])
        prompt = f"""
        Tu é a doutora investigadora. 
        CRUZAMENTO DE DADOS:
        1. Seus alvos previstos (brainstorming).
        2. Evidência real do PubMed: {lista_str}.
        
        MISSÃO: Criar a lista final técnica (60-80 itens). 
        REGRAS: Remova termos genéricos (RNA, DNA, Study, Effect, Results, Body), conectivos, animais (Toad, Rabbit) e siglas de 1 ou 2 letras.
        Mantenha siglas técnicas e fármacos específicos.
        OUTPUT: Retorne APENAS a lista Python final.
        """

    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.4 if fase=="brainstorming" else 0.2}}
    
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=20)
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                texto_limpo = re.search(r'\[.*\]', texto, re.DOTALL)
                if texto_limpo:
                    return ast.literal_eval(texto_limpo.group())
        except: continue
    return []

# --- 5. MINERAÇÃO E CRUZAMENTO MASSIVO ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    
    # 1. IA Pensa primeiro (Brainstorming)
    alvos_previstos = []
    if usar_ia: alvos_previstos = _doutora_investigadora(termo_base, fase="brainstorming")
    
    # 2. Busca Real no PubMed
    query = f"({termo_base} AND (Pharmacology OR Molecular)) AND (2018:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1500, sort="relevance")
        record = Entrez.read(handle); handle.close()
        total_docs = int(record["Count"])
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        # PENEIRA DE REGEX: Só aceita siglas/termos com 3+ letras (mata ruído como 'O', 'S', 'J')
        blacklist = {'THE', 'AND', 'FOR', 'WITH', 'FROM', 'ROLE', 'STUDY', 'DATA', 'CELL', 'LEVEL', 'RESULTS'}
        candidatos_raw = re.findall(r'\b[A-Z][A-Z0-9-]{2,}\b', full_data) 
        
        contagem = Counter([c.upper() for c in candidatos_raw if c.upper() not in blacklist])
        top_pubmed = [termo for termo, freq in contagem.most_common(150)]
        
        # 3. Cruzamento final
        nomes_finais = []
        if usar_ia:
            lista_para_cruzamento = list(set(alvos_previstos + top_pubmed))
            nomes_finais = _doutora_investigadora(termo_base, lista_pubmed=lista_para_cruzamento, fase="cruzamento")
        else:
            nomes_finais = top_pubmed[:40]

        if not nomes_finais: return [] 

        res_finais = []
        for nome in nomes_finais:
            n_limpo = limpar_termo_para_pubmed(nome)
            freq = contagem.get(n_limpo.upper(), 1)
            # Cálculo de Blue Ocean baseado na frequência encontrada
            l_score, p_val, b_ocean, status = calcular_metricas_originais(freq, total_docs, freq * 10)
            res_finais.append({
                "Alvo": n_limpo, 
                "Lambda": round(l_score, 2), 
                "P-value": round(p_val, 4), 
                "Blue Ocean": round(b_ocean, 1), 
                "Status": status
            })
        return res_finais
    except: return []

# --- 6. ANÁLISE IA ---
def analisar_abstract_com_ia(titulo, dados, api_key, lang='pt'):
    if not api_key: return "Chave necessária."
    prompt = f"Analise como PhD em farmacologia: {titulo}. Contexto: {dados}. Resuma Alvo, Fármaco e Efeito Funcional em 20 palavras."
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
