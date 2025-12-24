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
    "gemini-2.5-flash", "gemini-2.0-flash", "gemini-2.0-flash-exp", "gemini-flash-latest"
]

# --- 1. MOTOR DE MÉTRICAS ---
def calcular_metricas_originais(freq, total_docs, n_alvo_total):
    lambda_score = (freq / total_docs) * 100 if total_docs > 0 else 0
    p_value = math.exp(-freq/5) if freq > 0 else 1.0
    blue_ocean = max(0, 100 - (n_alvo_total / 10)) if n_alvo_total > 0 else 100.0
    status = "Saturado" if blue_ocean < 25 else "Blue Ocean" if blue_ocean > 75 else "Competitivo"
    return lambda_score, p_value, blue_ocean, status

def limpar_termo_para_pubmed(termo):
    return re.sub(r'^[^a-zA-Z0-9]+|[^a-zA-Z0-9]+$', '', termo.strip())

# --- 2. DOUTORA INVESTIGADORA (PROMPT RECALIBRADO) ---
def _chamar_ia_estrita(prompt):
    api_key = st.session_state.get('api_key_usuario', '').strip()
    if not api_key: return []
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.1}}
    
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, json=data, timeout=25)
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                match = re.search(r'\[.*\]', texto, re.DOTALL)
                if match: return ast.literal_eval(match.group())
        except: continue
    return []

# --- 3. MINERAÇÃO MASSIVA (PENEIRA DE 2025) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    
    # Brainstorming focado apenas no que é BOLA DA VEZ
    alvos_previstos = []
    if usar_ia:
        prompt_brain = f"""
        As a Senior PhD in Molecular Pharmacology, provide a Python list of 50 specific molecular targets (channels, receptors, enzymes) for research on {termo_base}.
        STRICT EXCLUSION: No classic physiology (Sodium, Water, Toad, Frog, Vasopressin, Aldosterone).
        FOCUS: Bench research, Ion Channels (TRP, Piezo, P2X), miRNAs, Inflammasomes.
        Output: ONLY a Python list [].
        """
        alvos_previstos = _chamar_ia_estrita(prompt_brain)
    
    # Busca PubMed focada do ano 2000 para frente (mata os sapos de 1960)
    query = f"({termo_base} AND (Pharmacology OR Molecular OR Signaling)) AND (2015:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1000, sort="relevance")
        record = Entrez.read(handle); handle.close()
        total_docs = int(record["Count"])
        
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        # SUPER BLACKLIST DE EXTERMÍNIO
        blacklist_pesada = {
            'TOAD', 'FROG', 'DOG', 'RABBIT', 'SODIUM', 'WATER', 'TRANSPORT', 'PERMEABILITY',
            'VASOPRESSIN', 'ALDOSTERONE', 'KIDNEY', 'PITUITARY', 'POSTERIOR', 'SYSTEM',
            'ACTION', 'EFFECT', 'STUDY', 'THE', 'AND', 'MUSCLE', 'CELL', 'BIOLOGICAL',
            'ISOLATED', 'EXPERIMENTAL', 'PHYSIOLOGY', 'HORMONES', 'ISOTOPES', 'ACID'
        }
        
        candidatos_pubmed = []
        for artigo in full_data.split("\n\nPMID-"):
            texto_util = ""
            for linha in artigo.split("\n"):
                if linha.startswith("TI  - ") or linha.startswith("KW  - ") or linha.startswith("OT  - "):
                    texto_util += linha[6:].strip() + " "
            
            # Captura siglas técnicas (3+ letras)
            encontrados = re.findall(r'\b[A-Z0-9-]{3,}\b', texto_util)
            for t in encontrados:
                if t.upper() not in blacklist_pesada:
                    candidatos_pubmed.append(t.upper())

        contagem = Counter(candidatos_pubmed)
        top_pubmed = [termo for termo, freq in contagem.most_common(150)]
        
        # Cruzamento final: Hipótese (IA) vs Evidência (PubMed)
        nomes_finais = []
        if usar_ia:
            lista_cruzamento = list(set(alvos_previstos + top_pubmed))
            prompt_cross = f"PhD Cross-reference: From this list {lista_cruzamento[:120]}, keep ONLY specific molecular targets and drugs for {termo_base}. Delete all general physiology and system terms. Output: ONLY a Python list []."
            nomes_finais = _chamar_ia_estrita(prompt_cross)
        
        if not nomes_finais: nomes_finais = top_pubmed[:60]

        res_finais = []
        for nome in nomes_finais:
            n_limpo = limpar_termo_para_pubmed(nome)
            if n_limpo.upper() in blacklist_pesada: continue
            freq = contagem.get(n_limpo.upper(), 1)
            l_score, p_val, b_ocean, status = calcular_metricas_originais(freq, total_docs, freq * 10)
            res_finais.append({"Alvo": n_limpo, "Lambda": round(l_score, 2), "P-value": round(p_val, 4), "Blue Ocean": round(b_ocean, 1), "Status": status})
        return res_finais
    except: return []
