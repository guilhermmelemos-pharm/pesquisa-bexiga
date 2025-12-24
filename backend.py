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

# --- 2. DOUTORA INVESTIGADORA (FILTRO DE ALTO IMPACTO) ---
def _chamar_ia_seletiva(prompt):
    api_key = st.session_state.get('api_key_usuario', '').strip()
    if not api_key: return []
    headers = {'Content-Type': 'application/json'}
    # Temperatura zero para máxima precisão técnica
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.0}}
    
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

# --- 3. MINERAÇÃO COM FILTRO DE ÉPOCA E QUALIDADE ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    
    # Brainstorming: A IA define o que é o "Estado da Arte" antes de ver o lixo do PubMed
    alvos_elite = []
    if usar_ia:
        prompt_brain = f"""
        Role: PhD Molecular Pharmacologist.
        Task: List 60 specific molecular targets (channels, receptors, enzymes) and specific experimental drugs for "{termo_base}" based ONLY on 2020-2026 science.
        STRICT EXCLUSION: No classic physiology (Toads, Frogs, Dogs, Sodium/Water transport, Vasopressin).
        CORE TARGETS: TRPV4, Piezo1, P2X3, ROCK, NLRP3, SPHK1, Rho-kinase, miRNAs.
        Output: ONLY a Python list of strings [].
        """
        alvos_elite = _chamar_ia_seletiva(prompt_brain)
    
    # Busca PubMed: Só 2018 para frente. Isso mata a fisiologia clássica de 1970 na fonte.
    query = f"({termo_base} AND (Pharmacology OR Molecular OR Signaling)) AND (2018:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1200, sort="relevance")
        record = Entrez.read(handle); handle.close()
        total_docs = int(record["Count"])
        
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        # Super Blacklist (Termos que não são alvos de bancada modernos)
        lixo = {
            'TOAD', 'FROG', 'DOG', 'RABBIT', 'SODIUM', 'WATER', 'TRANSPORT', 'KIDNEY', 
            'STUDY', 'ROLE', 'LAB', 'EXPERIMENTAL', 'PHYSIOLOGY', 'HORMONES', 
            'BLADDER', 'URINARY', 'EPITHELIUM', 'SYSTEM', 'ACTION', 'EFFECT'
        }
        
        candidatos_pubmed = []
        for artigo in full_data.split("\n\nPMID-"):
            # FOCO TOTAL: Título e Keywords (onde está a identidade molecular do artigo)
            texto_util = ""
            for linha in artigo.split("\n"):
                if linha.startswith("TI  - ") or linha.startswith("KW  - ") or line.startswith("OT  - "):
                    texto_util += linha[6:].strip() + " "
            
            # Captura siglas e nomes químicos (3+ caracteres)
            siglas = re.findall(r'\b[A-Z0-9-]{3,12}\b', texto_util)
            for s in siglas:
                if s.upper() not in lixo:
                    candidatos_pubmed.append(s.upper())

        contagem = Counter(candidatos_pubmed)
        top_pubmed = [termo for termo, freq in contagem.most_common(150)]
        
        # O PULO DO GATO: Cruzamento Rígido com ordem de extermínio
        nomes_finais = []
        if usar_ia:
            lista_para_cruzamento = list(set(alvos_elite + top_pubmed))
            prompt_cross = f"""
            Analyze this list of terms from PubMed: {lista_para_cruzamento[:120]}.
            MISSION: Keep ONLY specific Ion Channels, Receptors, Enzymes, and Novel Drugs.
            DELETE 80% of the list: Remove general physiology, anatomy, and classic 20th-century terms.
            Output: ONLY a Python list of strings [].
            """
            nomes_finais = _chamar_ia_seletiva(prompt_cross)
        
        if not nomes_finais: nomes_finais = top_pubmed[:60]

        res_finais = []
        for nome in nomes_finais:
            if len(nome) < 3 or nome.upper() in lixo: continue
            freq = contagem.get(nome.upper(), 1)
            l_score, p_val, b_ocean, status = calcular_metricas_originais(freq, total_docs, freq * 10)
            res_finais.append({"Alvo": nome, "Lambda": round(l_score, 2), "P-value": round(p_val, 4), "Blue Ocean": round(b_ocean, 1), "Status": status})
        return res_finais
    except: return []
