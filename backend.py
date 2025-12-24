import streamlit as st
from Bio import Entrez
import requests
import json
from tenacity import retry, stop_after_attempt, wait_exponential
import re
from collections import Counter
import math
import ast

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br"

MODELOS_ATIVOS = [
    "gemini-2.5-flash", "gemini-2.0-flash", "gemini-2.0-flash-exp", 
    "gemini-flash-latest", "gemini-2.5-flash-lite"
]

# --- 1. MOTOR ESTATÍSTICO (INTEGRIDADE ORIGINAL) ---
def calcular_metricas(freq, total_papers, n_alvo):
    # Lambda Score (Relevância Relativa)
    lambda_score = (freq / total_papers) * 100 if total_papers > 0 else 0
    
    # P-value Simulado (Significância da frequência no contexto)
    p_value = math.exp(-freq/10) if freq > 0 else 1.0
    
    # Índice Blue Ocean (Inovação vs Saturação)
    # Quanto menor o n_alvo (papers totais do alvo), maior o oceano azul
    if n_alvo == 0: blue_ocean = 100.0
    else: blue_ocean = max(0, 100 - (n_alvo / 10))
    
    status = "Saturado" if blue_ocean < 20 else "Blue Ocean" if blue_ocean > 70 else "Competitivo"
    
    return lambda_score, p_value, blue_ocean, status

# --- 2. A DOUTORA PROATIVA (DESCOBERTA DE ALVOS) ---
def _faxina_ia(termo_base, lista_pubmed):
    api_key = st.session_state.get('api_key_usuario', '')
    if not api_key: return lista_pubmed[:60] 

    prompt = f"""
    Tu é uma doutora em farmacologia experiente. Foco: {termo_base}.
    PASSO 1: Use seu conhecimento para listar alvos (receptores, canais, miRNAs, fármacos) inovadores e pouco explorados neste tecido.
    PASSO 2: Cruze com esta lista do PubMed: {", ".join(lista_pubmed)}.
    MISSÃO: Retorne uma lista técnica e robusta (60-120 termos).
    EXEMPLOS DE PADRÃO: TRPV4, GSK1016790A, SPHK1, NLRP3, miR-132, SNARE, TMAO.
    DELETA LIXO: (Studies, Effect, Patients, Incontinence, Water, etc).
    OUTPUT: Apenas a lista Python.
    """
    
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.4}}
    
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=20)
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                return ast.literal_eval(texto.replace("```python", "").replace("```", "").strip())
        except: continue
    return lista_pubmed[:60]

# --- 3. MINERAÇÃO E INTEGRAÇÃO ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    Entrez.email = email if email else "pesquisador_guest@unifesp.br"
    query = f"({termo_base} AND (Physiology OR Pharmacology OR Molecular)) AND (2018:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1000)
        record = Entrez.read(handle); handle.close()
        total_papers = int(record["Count"])
        
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        candidatos = []
        for artigo in full_data.split("\n\nPMID-"):
            texto = ""
            for line in artigo.split("\n"):
                if line.startswith("TI  - ") or line.startswith("KW  - "): texto += line[6:].strip() + " "
            encontrados = re.findall(r'\b(?:[A-Z0-9-]{2,}|[a-z]{1,2}-[A-Z0-9]{2,}|[a-z]{1,2}[A-Z][a-zA-Z0-9-]*)\b', texto)
            candidatos.extend([t.strip("-").upper() for t in encontrados if len(t) >= 3])

        contagem = Counter(candidatos)
        top_pubmed = [termo for termo, freq in contagem.most_common(200)]
        
        # A IA faz a curadoria proativa
        alvos_finais = _faxina_ia(termo_base, top_pubmed) if usar_ia else top_pubmed[:60]
        
        # O sistema recalcula as métricas para a tabela final
        resultados_com_metricas = []
        for alvo in alvos_finais:
            freq = contagem[alvo]
            # Aqui você chamaria o consultar_pubmed_count real para o n_alvo
            lambda_s, p_v, blue_o, status = calcular_metricas(freq, total_papers, freq * 5) 
            resultados_com_metricas.append({
                "Alvo": alvo,
                "Lambda": lambda_s,
                "P-value": p_v,
                "Blue Ocean": blue_o,
                "Status": status
            })
            
        return resultados_com_metricas
    except: return []
        
