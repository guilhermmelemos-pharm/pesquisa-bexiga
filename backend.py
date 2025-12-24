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
    # Temperatura baixa (0.1) para evitar que a IA "invente" ou seja boazinha com lixo
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

# --- 3. MINERAÇÃO COM FILTRO DE ÉPOCA E QUALIDADE ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    
    # Brainstorming: A IA define o que é "bom" antes de ver o lixo do PubMed
    alvos_elite = []
    if usar_ia:
        prompt_brain = f"""
        Como PhD em Farmacologia Molecular, liste 50 alvos moleculares (canais iônicos, receptores, miRNAs) e fármacos específicos de 2020-2026 para {termo_base}.
        REGRAS DE EXCLUSÃO: Delete fisiologia animal (sapos, cães), termos anatômicos (bexiga, músculo), e biologia genérica (DNA, RNA, Proteína).
        FOCO: Piezo1, TRPV4, P2X3, ROCK, NLRP3, SPHK1, etc.
        Output: Apenas lista Python [].
        """
        alvos_elite = _chamar_ia_seletiva(prompt_brain)
    
    # Busca PubMed: Só 2018 para frente. Mata o passado.
    query = f"({termo_base} AND (Pharmacology OR Molecular OR Signaling)) AND (2018:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1200, sort="relevance")
        record = Entrez.read(handle); handle.close()
        total_docs = int(record["Count"])
        
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        # Super Blacklist (Filtro de Segurança Adicional)
        lixo = {'TOAD', 'FROG', 'DOG', 'RABBIT', 'SODIUM', 'WATER', 'TRANSPORT', 'KIDNEY', 'RNA', 'DNA', 'CELL', 'STUDY', 'ROLE', 'LAB', 'EXPERIMENTAL', 'PHYSIOLOGY', 'HORMONES', 'BLADDER', 'URINARY', 'EPITHELIUM'}
        
        candidatos_pubmed = []
        for artigo in full_data.split("\n\nPMID-"):
            # Analisa apenas Título e Keywords (onde está a essência)
            texto_util = ""
            for linha in artigo.split("\n"):
                if linha.startswith("TI  - ") or linha.startswith("KW  - ") or linha.startswith("OT  - "):
                    texto_util += linha[6:].strip() + " "
            
            siglas = re.findall(r'\b[A-Z0-9-]{3,12}\b', texto_util)
            for s in siglas:
                if s.upper() not in lixo:
                    candidatos_pubmed.append(s.upper())

        contagem = Counter(candidatos_pubmed)
        top_pubmed = [termo for termo, freq in contagem.most_common(150)]
        
        # O PULO DO GATO: Cruzamento Rígido
        nomes_finais = []
        if usar_ia:
            lista_cruzamento = list(set(alvos_elite + top_pubmed))
            prompt_cross = f"""
            Analise esta lista de termos extraídos do PubMed: {lista_cruzamento[:120]}.
            Sua missão é DELETAR 80% da lista. Mantenha APENAS:
            1. Canais Iônicos (ex: TRPV4, Piezo1).
            2. Receptores e Enzimas (ex: P2X3, ROCK).
            3. Fármacos experimentais e novos (ex: GSK1016790A).
            Delete qualquer termo genérico, animal, ou metabólico básico.
            Output: Apenas lista Python [].
            """
            nomes_finais = _chamar_ia_seletiva(prompt_cross)
        
        if not nomes_finais: nomes_finais = top_pubmed[:60]

        res_finais = []
        for nome in nomes_finais:
            # Filtro final triplo
            if nome.upper() in lixo or len(nome) < 3: continue
            freq = contagem.get(nome.upper(), 1)
            l_score, p_val, b_ocean, status = calcular_metricas_originais(freq, total_docs, freq * 10)
            res_finais.append({"Alvo": nome, "Lambda": round(l_score, 2), "P-value": round(p_val, 4), "Blue Ocean": round(b_ocean, 1), "Status": status})
        return res_finais
    except: return []
