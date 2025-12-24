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

# --- 3. RADAR DE NOTÍCIAS (FOCO 2025/2026) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_todas_noticias(lang='pt'):
    try:
        query = "(bladder pharmacology OR urinary molecular targets) AND (2024:2026[Date - Publication])"
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
                news.append({"titulo": tit, "fonte": journal[:40], "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return news
    except: return []

# --- 4. CONTAGEM (FIX ERRO 126) ---
@st.cache_data(ttl=86400, show_spinner=False)
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    termo_limpo = limpar_termo_para_pubmed(termo)
    # Garante que a contagem respeite o foco em coisas novas (mínimo ano 2000)
    ano_partida = max(int(ano_ini), 2000)
    query = f"({termo_limpo}) AND ({ano_partida}:{ano_fim}[Date - Publication])"
    if contexto: query += f" AND ({contexto})"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle); handle.close()
        return int(record["Count"])
    except: return 0

# --- 5. DOUTORA INVESTIGADORA (FILTRAGEM PRÉVIA DE ELITE) ---
def _doutora_investigadora(termo_base, lista_pubmed=None, fase="brainstorming"):
    api_key = st.session_state.get('api_key_usuario', '')
    if not api_key: return []

    if fase == "brainstorming":
        prompt = f"""
        Role: PhD Senior Scientist in Pharmacology (2025 standards).
        Task: Create a list of 60 MODERN molecular targets and innovative drugs for "{termo_base}".
        Focus: Ion channels (TRP, Piezo, P2X, Nav, Kv), miRNAs, SPHK, Inflammasomes, and G-protein receptors.
        Strict Exclusion: NO classic 20th-century physiology (No Toads, Dogs, Vasopressin, general Sodium/Water transport).
        Output: ONLY a Python list of strings.
        """
    else:
        lista_str = ", ".join(lista_pubmed[:150])
        prompt = f"""
        Role: Senior Investigating Scientist.
        Task: Cross-reference your expert knowledge with this PubMed 2018-2025 list: {lista_str}.
        Rules: 
        1. Keep only targets with molecular bench-top relevance.
        2. Delete anything that sounds like basic anatomy or clinical surgery (Surgery, Bladder, Urine, Catheter).
        3. Prioritize "Hot Targets" like Piezo1, TRPV4, SPHK1, GSK1016790A.
        Output: ONLY a Python list of strings [].
        """

    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.3}}
    
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

# --- 6. MINERAÇÃO PROATIVA (O FILTRO DE 2000+) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    
    # 1. IA define o que é "Novo e Bom" antes de olhar o PubMed
    alvos_elite = []
    if usar_ia: alvos_elite = _doutora_investigadora(termo_base, fase="brainstorming")
    
    # 2. Busca PubMed focada no período moderno (2018-2026)
    query = f"({termo_base} AND (Pharmacology OR Molecular OR Signaling)) AND (2018:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1000, sort="relevance")
        record = Entrez.read(handle); handle.close()
        total_docs = int(record["Count"])
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        # Peneira de Época: Deleta termos que remetem à fisiologia datada
        blacklist = {
            'TOAD', 'FROG', 'DOG', 'CAT', 'SODIUM', 'WATER', 'TRANSPORT', 'KIDNEY', 
            'VASOPRESSIN', 'ALDOSTERONE', 'URETHRA', 'SURGERY', 'STUDY', 'ROLE', 
            'EFFECT', 'BLOOD', 'BLADDER', 'URINARY', 'PURPOSE', 'BACKGROUND'
        }
        
        candidatos_pubmed = []
        for artigo in full_data.split("\n\nPMID-"):
            texto_util = ""
            for line in artigo.split("\n"):
                if line.startswith("TI  - ") or line.startswith("KW  - ") or line.startswith("OT  - "):
                    texto_util += line[6:].strip() + " "
            
            encontrados = re.findall(r'\b[A-Z0-9-]{3,}\b', texto_util)
            for t in encontrados:
                if t.upper() not in blacklist:
                    candidatos_pubmed.append(t.upper())

        contagem = Counter(candidatos_pubmed)
        top_pubmed = [termo for termo, freq in contagem.most_common(180)]
        
        # 3. Cruzamento: O que a IA sugeriu vs O que o PubMed confirmou
        nomes_finais = []
        if usar_ia:
            # Mistura os alvos de elite da IA com os achados reais do PubMed
            lista_para_cruzamento = list(set(alvos_elite + top_pubmed))
            nomes_finais = _doutora_investigadora(termo_base, lista_pubmed=lista_para_cruzamento, fase="cruzamento")
        
        if not nomes_finais: nomes_finais = top_pubmed[:60]

        res_finais = []
        for nome in nomes_finais:
            n_limpo = limpar_termo_para_pubmed(nome)
            if n_limpo.upper() in blacklist: continue 
            freq = contagem.get(n_limpo.upper(), 1)
            l_score, p_val, b_ocean, status = calcular_metricas_originais(freq, total_docs, freq * 10)
            res_finais.append({"Alvo": n_limpo, "Lambda": round(l_score, 2), "P-value": round(p_val, 4), "Blue Ocean": round(b_ocean, 1), "Status": status})
        return res_finais
    except: return []

# --- 7. ANÁLISE IA (RESUMO PHD) ---
def analisar_abstract_com_ia(titulo, dados, api_key, lang='pt'):
    if not api_key: return "API Key Required."
    # Forçamos a análise a focar em inovação molecular
    prompt = f"Advanced Molecular Analysis: {titulo}. Context: {dados}. Identify Target, Novel Compound, and Benchwork Result. Max 20 words."
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}]}
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=12)
            if resp.status_code == 200:
                return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
        except: continue
    return "AI Error."
