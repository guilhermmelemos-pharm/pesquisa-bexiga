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
    "gemini-2.5-flash", 
    "gemini-2.0-flash", 
    "gemini-2.0-flash-exp", 
    "gemini-flash-latest"
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

# --- 3. DOUTORA INVESTIGADORA (COM TIMEOUT E LIMPEZA ROBUSTA) ---
def _doutora_investigadora(termo_base, lista_pubmed=None, fase="brainstorming"):
    api_key = st.session_state.get('api_key_usuario', '').strip()
    if not api_key: return []

    if fase == "brainstorming":
        prompt = f"As a PhD in Pharmacology, list 60 specific molecular targets (channels, receptors, enzymes) and drugs for {termo_base}. Focus: bench research. Output: strictly a Python list of strings []."
    else:
        lista_str = ", ".join(lista_pubmed[:100])
        prompt = f"Cross-reference your knowledge with this PubMed list: {lista_str}. Keep only pharmacological targets and drugs. No clinical noise. Output: strictly a Python list of strings []."

    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.2}}
    
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            # Timeout aumentado para 30s para evitar o erro de "Vazio"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=30)
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                # Tenta capturar a lista mesmo se a IA vier com conversa
                match = re.search(r'\[.*\]', texto, re.DOTALL)
                if match:
                    res = ast.literal_eval(match.group())
                    if isinstance(res, list): return res
        except: continue
    return []

# --- 4. MINERAÇÃO MASSIVA (COM REDUNDÂNCIA) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    
    # PASSO 1: IA Pensa
    alvos_previstos = []
    if usar_ia: alvos_previstos = _doutora_investigadora(termo_base, fase="brainstorming")
    
    # PASSO 2: Busca PubMed
    query = f"({termo_base} AND (Pharmacology OR Molecular)) AND (2018:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1000, sort="relevance")
        record = Entrez.read(handle); handle.close()
        total_docs = int(record["Count"])
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        # Filtro inicial de lixo
        blacklist = {'THE', 'AND', 'FOR', 'WITH', 'FROM', 'ROLE', 'STUDY', 'DATA', 'CELL', 'LEVEL', 'RESULTS', 'TOAD', 'FROG', 'DOG'}
        candidatos_raw = re.findall(r'\b[A-Z0-9-]{3,}\b', full_data) 
        contagem = Counter([c.upper() for c in candidatos_raw if c.upper() not in blacklist])
        top_pubmed = [termo for termo, freq in contagem.most_common(150)]
        
        # PASSO 3: Cruzamento final
        nomes_finais = []
        if usar_ia:
            lista_para_cruzamento = list(set(alvos_previstos + top_pubmed))
            nomes_finais = _doutora_investigadora(termo_base, lista_pubmed=lista_para_cruzamento, fase="cruzamento")
        
        # REDUNDÂNCIA: Se a IA falhar (timeout ou erro), usamos o Top PubMed para não dar "Vazio"
        if not nomes_finais:
            nomes_finais = top_pubmed[:60]

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

# --- 5. FUNÇÕES DE APOIO ---
@st.cache_data(ttl=86400, show_spinner=False)
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    termo_limpo = limpar_termo_para_pubmed(termo)
    query = f"({termo_limpo})"
    if contexto: query += f" AND ({contexto})"
    query += f" AND ({ano_ini}:{ano_fim}[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle); handle.close()
        return int(record["Count"])
    except: return 0

def analisar_abstract_com_ia(titulo, dados, api_key, lang='pt'):
    # Mesma lógica robusta para análise individual
    prompt = f"Summarize Target and Drug for: {titulo}. Context: {dados}. Limit: 20 words."
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}]}
    url = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.0-flash:generateContent?key={api_key}"
    try:
        resp = requests.post(url, headers=headers, json=data, timeout=15)
        return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
    except: return "Analysis unavailable."

# RADAR
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_todas_noticias(lang='pt'):
    try:
        query = "(bladder physiology OR pharmacology) AND (2024:2026[Date - Publication])"
        handle = Entrez.esearch(db="pubmed", term=query, retmax=6, sort="pub_date")
        record = Entrez.read(handle); handle.close()
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
                news.append({"titulo": tit, "fonte": journal[:40], "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return news
    except: return []
