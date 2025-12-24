import streamlit as st  # FIX: Adicionado para evitar o erro de importação
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

def limpar_termo_para_pubmed(termo):
    return re.sub(r'^[^a-zA-Z0-9]+|[^a-zA-Z0-9]+$', '', termo.strip())

# --- 2. DOUTORA INVESTIGADORA (FILTRO DE ELITE) ---
def _chamar_ia_seletiva(prompt):
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

# --- 3. MINERAÇÃO E CRUZAMENTO (PENEIRA DE 2025) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    
    # Brainstorming: IA projeta alvos modernos antes de ver o lixo do passado
    alvos_elite = []
    if usar_ia:
        prompt_brain = f"""
        Como PhD em Farmacologia Molecular, liste 60 alvos (canais, receptores, miRNAs) e fármacos de 2020-2026 para {termo_base}.
        EXCLUSÃO CRÍTICA: Não use anatomia genérica (músculo, nervo, bexiga) nem animais (sapos, cães).
        FOCO: TRPV4, Piezo1, P2X3, NLRP3, SPHK1, Rho-kinase, ROCK.
        Output: Apenas lista Python [].
        """
        alvos_elite = _chamar_ia_seletiva(prompt_brain)
    
    # Busca PubMed: Só 2018 para frente para matar a fisiologia clássica
    query = f"({termo_base} AND (Pharmacology OR Molecular OR Signaling)) AND (2018:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1200, sort="relevance")
        record = Entrez.read(handle); handle.close()
        total_docs = int(record["Count"])
        
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        # Super Blacklist (Filtro de Extermínio Semântico)
        lixo = {
            'TOAD', 'FROG', 'DOG', 'RABBIT', 'SODIUM', 'WATER', 'TRANSPORT', 'KIDNEY', 'CELL', 'MUSCLE', 
            'NERVE', 'STUDY', 'ROLE', 'SYSTEM', 'ACTION', 'EFFECT', 'URINARY', 'BLADDER', 'EPITHELIUM',
            'PHARMACOLOGY', 'THE', 'AND', 'ISOLATED', 'EXPERIMENTAL', 'PHYSIOLOGY', 'HORMONES'
        }
        
        candidatos_pubmed = []
        for artigo in full_data.split("\n\nPMID-"):
            texto_util = ""
            for linha in artigo.split("\n"):
                # SÓ OLHA TÍTULO E KEYWORDS (Limpeza radical)
                if linha.startswith("TI  - ") or linha.startswith("KW  - ") or linha.startswith("OT  - "):
                    texto_util += linha[6:].strip() + " "
            
            # Regex: Captura siglas técnicas (3 a 12 caracteres, evita lixo)
            encontrados = re.findall(r'\b[A-Z0-9-]{3,12}\b', texto_util)
            for t in encontrados:
                if t.upper() not in lixo:
                    candidatos_pubmed.append(t.upper())

        contagem = Counter(candidatos_pubmed)
        top_pubmed = [termo for termo, freq in contagem.most_common(150)]
        
        # Cruzamento: IA atua como o "Filtro Final"
        nomes_finais = []
        if usar_ia:
            lista_para_cruzamento = list(set(alvos_elite + top_pubmed))
            prompt_cross = f"PhD Cross-ref: From this list {lista_para_cruzamento[:120]}, keep ONLY specific molecular targets and drugs. DELETE 80% (general biology/anatomy). Output: ONLY Python list []."
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

# --- 4. FUNÇÃO DE CONTAGEM (FIX ERRO LINHA 126) ---
@st.cache_data(ttl=86400, show_spinner=False)
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    termo_limpo = limpar_termo_para_pubmed(termo)
    query = f"({termo_limpo}) AND ({ano_ini}:{ano_fim}[Date - Publication])"
    if contexto: query += f" AND ({contexto})"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle); handle.close()
        return int(record["Count"])
    except: return 0
