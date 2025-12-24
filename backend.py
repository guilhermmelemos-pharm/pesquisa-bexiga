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

# --- 2. EXTRATOR DE ELITE (ALVO / FÁRMACO / AÇÃO) ---
def _extrair_entidades_farmacologicas(lista_bruta):
    api_key = st.session_state.get('api_key_usuario', '').strip()
    if not api_key: return []
    
    prompt = f"""
    Como PhD em Farmacologia, analise estes Títulos e Keywords: {lista_bruta[:150]}
    
    TASK: Retorne uma lista de dicionários Python estritamente no formato:
    {{"Alvo": "Nome do Alvo", "Farmaco": "Fármaco associado ou 'N/A'", "Acao": "Mecanismo curto (ex: Agonista, Bloqueador, Inibidor)"}}
    
    REGRAS CRÍTICAS:
    1. DELETE: Toads, Frogs, Dogs, Rabbits, Sodium, Water, Muscle, Cell, Urinary, Metabolism.
    2. FOCO: Apenas farmacologia molecular e sinalização urotelial/detrusora moderna (2020-2025).
    3. FORMATO: Apenas a lista Python, sem texto explicativo.
    """
    
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.1}}
    
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, json=data, timeout=20)
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                match = re.search(r'\[.*\]', texto, re.DOTALL)
                if match: return ast.literal_eval(match.group())
        except: continue
    return []

# --- 3. MINERAÇÃO SNIPER (SÓ TÍTULO E KEYWORDS 2020+) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    query = f"({termo_base} AND (Pharmacology OR Molecular)) AND (2020:2026[Date])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1000, sort="relevance")
        record = Entrez.read(handle); handle.close()
        total_docs = int(record["Count"])
        
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        texto_acumulado = []
        for artigo in full_data.split("\n\nPMID-"):
            # SÓ TÍTULO E KEYWORDS
            for line in artigo.split("\n"):
                if line.startswith("TI  - ") or line.startswith("KW  - ") or line.startswith("OT  - "):
                    texto_acumulado.append(line[6:].strip())

        # Envia para a IA formatar como Alvo/Fármaco/Ação
        entidades = _extrair_entidades_farmacologicas(texto_acumulado)
        
        res_finais = []
        for ent in entidades:
            # Reutiliza o motor de métricas para cada alvo extraído
            l_score, p_val, b_ocean, status = calcular_metricas_originais(1, total_docs, 10)
            res_finais.append({
                "Alvo": ent.get("Alvo"),
                "Farmaco": ent.get("Farmaco"),
                "Acao": ent.get("Acao"),
                "Lambda": round(l_score, 2),
                "Blue Ocean": round(b_ocean, 1),
                "Status": status
            })
        return res_finais
    except:
        return []

# --- 4. RADAR DE NOTÍCIAS (MANTIDO) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_todas_noticias(lang='pt'):
    try:
        query = "(bladder pharmacology OR ion channels) AND (2024:2026[Date])"
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
