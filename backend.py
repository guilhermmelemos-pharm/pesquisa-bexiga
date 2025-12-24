import streamlit as st  # PRECISA ESTAR AQUI NO TOPO
from Bio import Entrez
import requests
import json
import re
from collections import Counter
import ast
import math

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br"

# MODELOS DO SEU PAID TIER
MODELOS_ATIVOS = ["gemini-2.0-flash", "gemini-2.0-flash-exp", "gemini-flash-latest"]

# --- 1. MOTOR DE MÉTRICAS ---
def calcular_metricas_originais(freq, total_docs, n_alvo_total):
    lambda_score = (freq / total_docs) * 100 if total_docs > 0 else 0
    p_value = math.exp(-freq/5) if freq > 0 else 1.0
    blue_ocean = max(0, 100 - (n_alvo_total / 10)) if n_alvo_total > 0 else 100.0
    status = "Saturado" if blue_ocean < 25 else "Blue Ocean" if blue_ocean > 75 else "Competitivo"
    return lambda_score, p_value, blue_ocean, status

# --- 2. DOUTORA INVESTIGADORA (EXPANSÃO CONTROLADA) ---
def _faxina_ia_elite(lista_bruta):
    api_key = st.session_state.get('api_key_usuario', '').strip()
    # Lista robusta de segurança caso a API falhe
    fallback = ["TRPV4", "SPHK1", "PIEZO1", "P2X3", "GSK1016790A", "NLRP3", "ROCK1", "MIR-132", "P2X7", "CGRP"]
    if not api_key: return fallback
    
    prompt = f"""
    Role: Senior PhD Pharmacologist.
    Input Data: {lista_bruta[:150]}
    
    TASK: Provide a Python list of the 40 MOST RELEVANT molecular targets and drugs.
    STRICT RULES:
    1. INCLUDE: Specific ion channels (TRP, Piezo, Nav, Kv), Receptors (P2X, P2Y, EP, Muscarinic), Signaling (SPHK1, RhoA, ROCK, AKT), and specific compounds.
    2. DELETE: General anatomy (Bladder, Muscle), animals (Toad, Frog), and generic biology (DNA, RNA, Cell).
    3. Output: Return ONLY the Python list [].
    """
    
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.4}}
    
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=20)
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                match = re.search(r'\[.*\]', texto, re.DOTALL)
                if match: return ast.literal_eval(match.group())
        except: continue
    return fallback

# --- 3. MINERAÇÃO SNIPER (SÓ CIÊNCIA MODERNA 2020+) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    
    # Foco total no moderno: 2020-2026
    query = f"({termo_base} AND (Pharmacology OR Molecular OR Signaling)) AND (2020:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1500, sort="relevance")
        record = Entrez.read(handle); handle.close()
        total_docs = int(record["Count"])
        
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        blacklist_fatal = {'TOAD', 'FROG', 'DOG', 'RABBIT', 'SODIUM', 'WATER', 'TRANSPORT', 'CELL', 'MUSCLE', 'BLADDER'}
        
        candidatos_pubmed = []
        for artigo in full_data.split("\n\nPMID-"):
            texto_sniper = ""
            for line in artigo.split("\n"):
                if line.startswith("TI  - ") or line.startswith("KW  - ") or line.startswith("OT  - "):
                    texto_sniper += line[6:].strip() + " "
            
            siglas = re.findall(r'\b[A-Z0-9-]{3,15}\b', texto_sniper)
            for t in siglas:
                if t.upper() not in blacklist_fatal and not t.isdigit():
                    candidatos_pubmed.append(t.upper())

        contagem = Counter(candidatos_pubmed)
        top_bruto = [termo for termo, freq in contagem.most_common(200)]
        
        nomes_finais = _faxina_ia_elite(top_bruto) if usar_ia else top_bruto[:60]

        res_finais = []
        for nome in nomes_finais:
            if len(nome) < 3 or nome.upper() in blacklist_fatal: continue
            freq = contagem.get(nome.upper(), 1)
            l_score, p_val, b_ocean, status = calcular_metricas_originais(freq, total_docs, freq * 10)
            res_finais.append({
                "Alvo": nome, "Lambda": round(l_score, 2), "P-value": round(p_val, 4), 
                "Blue Ocean": round(b_ocean, 1), "Status": status
            })
        return res_finais
    except: return []

# --- 4. FUNÇÕES DE APOIO ---
def analisar_abstract_com_ia(titulo, dados, api_key):
    prompt = f"Summarize mechanism of {titulo} in 15 words. Context: {dados}."
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}]}
    url = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.0-flash:generateContent?key={api_key}"
    try:
        resp = requests.post(url, headers=headers, json=data, timeout=12)
        return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
    except: return "Analysis unavailable."

@st.cache_data(ttl=86400, show_spinner=False)
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    query = f"({termo}) AND ({contexto}) AND ({ano_ini}:{ano_fim}[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle); handle.close()
        return int(record["Count"])
    except: return 0
