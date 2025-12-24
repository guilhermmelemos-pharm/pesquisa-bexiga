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
    "gemini-2.5-flash", "gemini-2.0-flash", "gemini-2.0-flash-exp"
]

# --- 1. MOTOR DE MÉTRICAS ---
def calcular_metricas_originais(freq, total_docs, n_alvo_total):
    lambda_score = (freq / total_docs) * 100 if total_docs > 0 else 0
    p_value = math.exp(-freq/5) if freq > 0 else 1.0
    blue_ocean = max(0, 100 - (n_alvo_total / 10)) if n_alvo_total > 0 else 100.0
    status = "Saturado" if blue_ocean < 25 else "Blue Ocean" if blue_ocean > 75 else "Competitivo"
    return lambda_score, p_value, blue_ocean, status

# --- 2. RADAR DE NOTÍCIAS (RESTAURADO 2025) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_todas_noticias(lang='pt'):
    try:
        query = "(bladder pharmacology OR urothelium signaling) AND (2024:2026[Date - Publication])"
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
            if tit and pmid:
                news.append({"titulo": tit, "fonte": journal[:40], "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return news
    except: return []

# --- 3. DOUTORA INVESTIGADORA (FILTRO FARMACÊUTICO RÍGIDO) ---
def _faxina_ia_sniper(lista_bruta):
    api_key = st.session_state.get('api_key_usuario', '').strip()
    if not api_key: return lista_bruta[:40]
    
    prompt = f"""
    Role: Senior Molecular Pharmacologist.
    Input: {lista_bruta[:150]}
    
    TASK: Strictly separate PHARMACOLOGICAL TARGETS from GARBAGE.
    
    RULES:
    1. KEEP ONLY: Ion Channels (TRPV4, P2X3), Receptors (M3, Beta-3), Specific Drugs (GSK1016790A, Mirabegron), Signaling (cAMP, NLRP3).
    2. DELETE GARBAGE: (Urea, Across, Report, Muscle, Water, Lithium, Urea, Studies, Isolated, Biological, III, Guinea, Rabbit).
    3. DELETE: All general physiology and anatomy terms.
    
    OUTPUT: Return strictly a Python list [].
    """
    
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.0}}
    
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

# --- 4. MINERAÇÃO SNIPER (SÓ TÍTULO E KEYWORDS) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    
    # Busca focada em 2015+ para evitar a "fisiologia de sapo" antiga
    query = f"({termo_base} AND (Pharmacology OR Molecular)) AND (2015:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1000, sort="relevance")
        record = Entrez.read(handle); handle.close()
        total_docs = int(record["Count"])
        
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        candidatos_pubmed = []
        for artigo in full_data.split("\n\nPMID-"):
            # FOCO ABSOLUTO: Só TI (Title) e OT/KW (Keywords)
            texto_sniper = ""
            for line in artigo.split("\n"):
                if line.startswith("TI  - ") or line.startswith("OT  - ") or line.startswith("KW  - "):
                    texto_sniper += line[6:].strip() + " "
            
            # Captura siglas e nomes químicos (3+ letras)
            encontrados = re.findall(r'\b[A-Z0-9-]{3,12}\b', texto_sniper)
            for t in encontrados:
                # Remove lixo óbvio antes da IA
                if not t.isdigit() and len(t) > 2:
                    candidatos_pubmed.append(t.upper())

        contagem = Counter(candidatos_pubmed)
        top_bruto = [termo for termo, freq in contagem.most_common(180)]
        
        # Filtro da IA
        nomes_finais = []
        if usar_ia:
            nomes_finais = _faxina_ia_sniper(top_bruto)
        
        if not nomes_finais: nomes_finais = top_bruto[:60]

        res_finais = []
        for nome in nomes_finais:
            freq = contagem.get(nome.upper(), 1)
            l_score, p_val, b_ocean, status = calcular_metricas_originais(freq, total_docs, freq * 10)
            res_finais.append({"Alvo": nome, "Lambda": round(l_score, 2), "P-value": round(p_val, 4), "Blue Ocean": round(b_ocean, 1), "Status": status})
        return res_finais
    except: return []

# --- 5. APOIO (FIX ERRO 126) ---
@st.cache_data(ttl=86400, show_spinner=False)
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    query = f"({termo}) AND ({ano_ini}:{ano_fim}[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle); handle.close()
        return int(record["Count"])
    except: return 0
