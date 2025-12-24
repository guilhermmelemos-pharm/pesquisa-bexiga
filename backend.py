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

MODELOS_ATIVOS = ["gemini-2.5-flash", "gemini-2.0-flash", "gemini-2.0-flash-exp"]

# --- 1. MOTOR DE MÉTRICAS ---
def calcular_metricas_originais(freq, total_docs, n_alvo_total):
    lambda_score = (freq / total_docs) * 100 if total_docs > 0 else 0
    p_value = math.exp(-freq/5) if freq > 0 else 1.0
    blue_ocean = max(0, 100 - (n_alvo_total / 10)) if n_alvo_total > 0 else 100.0
    status = "Saturado" if blue_ocean < 25 else "Blue Ocean" if blue_ocean > 75 else "Competitivo"
    return lambda_score, p_value, blue_ocean, status

# --- 2. FAXINEIRO IA (ORDEM DE EXTERMÍNIO) ---
def _faxina_ia_sniper(lista_bruta):
    api_key = st.session_state.get('api_key_usuario', '').strip()
    if not api_key: return lista_bruta[:40]
    
    prompt = f"""
    Como PhD em Farmacologia Molecular, analise esta lista: {lista_bruta[:150]}
    Mantenha APENAS alvos moleculares específicos (ex: TRPV4, P2X3, NLRP3) e fármacos experimentais (GSK1016790A).
    DELETE TUDO que for: fisiologia clássica, anatomia genérica, animais e termos de 100 anos atrás.
    Exemplos de LIXO para DELETAR: Toads, Frogs, Dogs, Muscle, Water, Sodium, Isolated, Studies, Metabolism.
    Output: Apenas lista Python [].
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

# --- 3. MINERAÇÃO SNIPER (FOCO 2018-2026 E EXCLUSÃO DE ABSTRACT) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    
    # Corte temporal radical: 2018 para frente
    query = f"({termo_base} AND (Pharmacology OR Molecular)) AND (2018:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1000, sort="relevance")
        record = Entrez.read(handle); handle.close()
        total_docs = int(record["Count"])
        
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        # BLACKLIST DE EXTERMÍNIO (O que você não aguenta mais)
        lixo_fatal = {
            'TOAD', 'TOADS', 'FROG', 'FROGS', 'DOG', 'DOGS', 'RABBIT', 'SODIUM', 'WATER', 'TRANSPORT', 
            'ALDOSTERONE', 'VASOPRESSIN', 'KIDNEY', 'RNA', 'DNA', 'CELL', 'STUDY', 'STUDIES', 'ROLE', 
            'LAB', 'EXPERIMENTAL', 'PHYSIOLOGY', 'HORMONES', 'BLADDER', 'URINARY', 'EPITHELIUM', 
            'SYSTEM', 'ACTION', 'EFFECT', 'ACTIVE', 'POSTERIOR', 'PITUITARY', 'ION', 'OSMOSIS', 
            'PACAP', 'ABSORPTION', 'MUSCLE', 'NERVE', 'NERVOUS', 'STIMULATION', 'ACID', 'COMPOUNDS'
        }
        
        candidatos_pubmed = []
        for artigo in full_data.split("\n\nPMID-"):
            texto_sniper = ""
            for line in artigo.split("\n"):
                # SÓ OLHA TÍTULO (TI) E KEYWORDS (KW/OT)
                if line.startswith("TI  - ") or line.startswith("KW  - ") or line.startswith("OT  - "):
                    texto_sniper += line[6:].strip() + " "
            
            # Captura siglas e nomes químicos (3+ letras)
            encontrados = re.findall(r'\b[A-Z0-9-]{3,12}\b', texto_sniper)
            for t in encontrados:
                t_up = t.upper()
                if t_up not in lixo_fatal and not t_up.isdigit():
                    candidatos_pubmed.append(t_up)

        contagem = Counter(candidatos_pubmed)
        top_bruto = [termo for termo, freq in contagem.most_common(180)]
        
        nomes_finais = _faxina_ia_sniper(top_bruto) if usar_ia else top_bruto[:60]

        res_finais = []
        for nome in nomes_finais:
            if nome.upper() in lixo_fatal: continue
            freq = contagem.get(nome.upper(), 1)
            l_score, p_val, b_ocean, status = calcular_metricas_originais(freq, total_docs, freq * 10)
            res_finais.append({"Alvo": nome, "Lambda": round(l_score, 2), "P-value": round(p_val, 4), "Blue Ocean": round(b_ocean, 1), "Status": status})
        return res_finais
    except: return []

# --- 4. RADAR DE NOTÍCIAS (FOCO EM CANAIS IÔNICOS 2025) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_todas_noticias(lang='pt'):
    try:
        query = "(bladder pharmacology OR urothelium ion channels OR mechanotransduction) AND (2024:2026[Date])"
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
