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

# --- 2. EXTRATOR SNIPER (ALVO / FÁRMACO / AÇÃO) ---
def analisar_abstract_com_ia(titulo, dados, api_key, lang='pt'):
    """
    Substitui a função de resumo antiga. 
    Agora ela extrai estritamente a relação Alvo/Fármaco/Ação do artigo.
    """
    prompt = f"""
    Baseado no Título: {titulo} e Dados: {dados}.
    Extraia APENAS: Alvo Molecular | Fármaco | Ação.
    Exemplo: TRPV4 | GSK1016790A | Agonista Urotelial.
    Se não houver fármaco, coloque 'N/A'. Seja extremamente curto.
    """
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.1}}
    url = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.0-flash:generateContent?key={api_key}"
    try:
        resp = requests.post(url, headers=headers, json=data, timeout=10)
        return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
    except: return "N/A | N/A | Erro na extração"

# --- 3. MINERAÇÃO SNIPER (SÓ TÍTULO/KW PÓS-2020) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    # CORTE TEMPORAL: Bloqueia qualquer coisa antes de 2020
    query = f"({termo_base} AND (Pharmacology OR Molecular)) AND (2020:2026[Date])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1000, sort="relevance")
        record = Entrez.read(handle); handle.close()
        total_docs = int(record["Count"])
        
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        # Blacklist de lixo histórico e genérico
        lixo = {'TOAD', 'FROG', 'DOG', 'RABBIT', 'SODIUM', 'WATER', 'MUSCLE', 'URINARY', 'CELL', 'DNA', 'RNA'}
        
        candidatos = []
        for artigo in full_data.split("\n\nPMID-"):
            texto_sniper = ""
            for line in artigo.split("\n"):
                # IGNORA O RESUMO (AB), foca em Título (TI) e Keywords (KW/OT)
                if line.startswith("TI  - ") or line.startswith("KW  - ") or line.startswith("OT  - "):
                    texto_sniper += line[6:].strip() + " "
            
            siglas = re.findall(r'\b[A-Z0-9-]{3,15}\b', texto_sniper)
            for s in siglas:
                if s.upper() not in lixo and not s.isdigit():
                    candidatos.append(s.upper())

        contagem = Counter(candidatos)
        # Retorna os top 15 alvos moleculares reais
        top_nomes = [t for t, f in contagem.most_common(15)]

        res_finais = []
        for nome in top_nomes:
            l_score, p_val, b_ocean, status = calcular_metricas_originais(contagem[nome], total_docs, 10)
            res_finais.append({
                "Alvo": nome, "Lambda": round(l_score, 2), 
                "Blue Ocean": round(b_ocean, 1), "Status": status
            })
        return res_finais
    except: return []

# --- 4. RADAR E APOIO ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_todas_noticias(lang='pt'):
    query = "(bladder pharmacology OR ion channels) AND (2024:2026[Date])"
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
