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
MODELOS_ATIVOS = ["gemini-2.0-flash", "gemini-2.0-flash-exp"]

# --- 1. MOTOR DE MÉTRICAS ---
def calcular_metricas_originais(freq, total_docs, n_alvo_total):
    lambda_score = (freq / total_docs) * 100 if total_docs > 0 else 0
    p_value = math.exp(-freq/5) if freq > 0 else 1.0
    blue_ocean = max(0, 100 - (n_alvo_total / 10)) if n_alvo_total > 0 else 100.0
    status = "Saturado" if blue_ocean < 25 else "Blue Ocean" if blue_ocean > 75 else "Competitivo"
    return lambda_score, p_value, blue_ocean, status

# --- 2. SNIPER DE TABELA (ALVO | FÁRMACO | AÇÃO) ---
def analisar_abstract_com_ia(titulo, dados, api_key, lang='pt'):
    """
    Extrator de alta precisão para o formato Alvo | Fármaco | Ação.
    """
    if not api_key: return "Erro | Key Ausente | Configure a API"
    
    prompt = f"""
    ANALISE MOLECULAR:
    Título: {titulo}
    Metadados: {dados}
    
    TASK: Retorne EXCLUSIVAMENTE no formato: ALVO | FÁRMACO | AÇÃO.
    REGRAS: 
    - Se não houver fármaco específico, use 'Endógeno'.
    - Seja ultra-conciso (máximo 3 palavras por campo).
    - FOCO: Alvos modernos (TRP, Piezo, P2X, ROCK).
    - DELETE: Termos genéricos, métodos ou animais.
    """
    
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.0}}
    url = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.0-flash:generateContent?key={api_key}"
    
    try:
        resp = requests.post(url, headers=headers, json=data, timeout=15)
        if resp.status_code == 200:
            res = resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
            # Garante que o separador existe no retorno
            return res if "|" in res else f"{res} | N/A | N/A"
        return "Erro | API | Status " + str(resp.status_code)
    except:
        return "N/A | N/A | Timeout da Extração"

# --- 3. MINERAÇÃO SNIPER (SÓ TÍTULO/KW PÓS-2020) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    # Bloqueio histórico: só 2020 para frente
    query = f"({termo_base} AND (Pharmacology OR Molecular)) AND (2020:2026[Date])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1000, sort="relevance")
        record = Entrez.read(handle); handle.close()
        total_docs = int(record["Count"])
        
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        # Blacklist agressiva para limpar a lista de tags
        lixo = {'TOAD', 'FROG', 'DOG', 'RABBIT', 'SODIUM', 'WATER', 'MUSCLE', 'CELL', 'DNA', 'RNA', 'THE', 'AND', 'FOR'}
        
        candidatos = []
        for artigo in full_data.split("\n\nPMID-"):
            texto_sniper = ""
            for line in artigo.split("\n"):
                if line.startswith("TI  - ") or line.startswith("KW  - ") or line.startswith("OT  - "):
                    texto_sniper += line[6:].strip() + " "
            
            siglas = re.findall(r'\b[A-Z0-9-]{3,15}\b', texto_sniper)
            for s in siglas:
                if s.upper() not in lixo and not s.isdigit():
                    candidatos.append(s.upper())

        contagem = Counter(candidatos)
        top_nomes = [t for t, f in contagem.most_common(20)] # Aumentado para 20 alvos

        res_finais = []
        for nome in top_nomes:
            l_score, p_val, b_ocean, status = calcular_metricas_originais(contagem[nome], total_docs, 10)
            res_finais.append({
                "Alvo": nome, "Lambda": round(l_score, 2), 
                "Blue Ocean": round(b_ocean, 1), "Status": status
            })
        return res_finais
    except: return []

# --- 4. RADAR DE NOTÍCIAS (2025) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_todas_noticias(lang='pt'):
    query = "(bladder pharmacology OR urothelium signaling) AND (2024:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=6, sort="pub_date")
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        news = []
        for art in dados.split("\n\nPMID-"):
            tit, journal, pmid = "", "", ""
            for line in art.split("\n"):
                if line.startswith("TI  - "): tit = line[6:].strip()
                if line.startswith("JT  - "): journal = line[3:].strip()
                if line.strip().isdigit() and not pmid: pmid = line.strip()
            if tit:
                news.append({"titulo": tit, "fonte": journal[:40], "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return news
    except: return []

# --- 5. DETALHES ---
@st.cache_data(ttl=86400, show_spinner=False)
def buscar_resumos_detalhados(termo, orgao, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    query = f"({termo}) AND ({orgao}) AND (2018:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=3, sort="relevance")
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        artigos = []
        for raw in dados.split("\n\nPMID-"):
            tit, pmid, abstract, keywords = "", "", "", ""
            for line in raw.split("\n"):
                if line.startswith("TI  - "): tit = line[6:].strip()
                if line.startswith("AB  - "): abstract = line[6:400].strip()
                if line.startswith("PMID- "): pmid = line[6:].strip()
                if line.startswith("KW  - ") or line.startswith("OT  - "): keywords += line[6:].strip()
            if tit:
                artigos.append({"Title": tit, "Info_IA": f"KW: {keywords} AB: {abstract}", "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return artigos
    except: return []
