import streamlit as st
from Bio import Entrez
import requests
import json
import re
from collections import Counter
import time
import ast

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br"

# --- LISTA DE MODELOS (Ordem de tentativa) ---
MODELOS_ATIVOS = [
    "gemini-2.0-flash", 
    "gemini-1.5-flash", 
    "gemini-1.5-pro"
]

# --- 1. FILTRO DE EXTERMÍNIO CLÍNICO ---
# Termos que morrem antes de chegar na IA
BLACKLIST_RADICAL = {
    "CANCER", "CARCINOMA", "TUMOR", "DIAGNOSIS", "THERAPY", "SURGERY", "SYNDROME",
    "PATIENT", "WOMEN", "MEN", "CHILDREN", "CLINICAL", "OUTCOME", "SAFETY", "EFFICACY",
    "UROLOGY", "CYSTECTOMY", "PROSTATE", "KIDNEY", "RENAL", "INCONTINENCE", "OAB", "LUTS",
    "OVERACTIVE", "TREATMENT", "MANAGEMENT", "CASE", "REPORT", "UPDATE", "REVIEW"
}

# --- 2. MONTAGEM DE URL SEGURA ---
def montar_url(modelo, chave):
    # Quebrado para evitar hyperlink automático no chat
    base = "https://generativelanguage" + ".googleapis.com/v1beta/models"
    return f"{base}/{modelo}:generateContent?key={chave}"

# --- 3. FAXINEIRO MOLECULAR (IA) ---
def _faxina_ia(lista_suja):
    api_key = st.session_state.get('api_key_usuario', '').strip()
    if not api_key: return lista_suja[:40]

    prompt = f"""
    ACT AS: PhD Senior Pharmacologist. 
    TASK: Review this list of terms from PubMed.
    GOAL: Keep ONLY specific molecular targets (receptors, enzymes, ions channels) and drugs.
    
    STRICT RULES:
    - KEEP: TRPV1, P2X3, ROCK, mTOR, Trehalose, Mirabegron, Muscarinic M3, etc.
    - DELETE: Anatomy (bladder, detrusor), Diseases (OAB, Cystitis), Concepts (Pathophysiology, Stress).
    - If you can't use it in an 'Organ Bath' or 'Western Blot', DELETE it.
    
    INPUT LIST: {", ".join(lista_suja)}
    
    OUTPUT: Return strictly a Python list of strings.
    """
    
    payload = {
        "contents": [{"parts": [{"text": prompt}]}],
        "generationConfig": {"temperature": 0.0}
    }
    
    for m in MODELOS_ATIVOS:
        try:
            url = montar_url(m, api_key)
            resp = requests.post(url, json=payload, timeout=15)
            if resp.status_code == 200:
                raw_text = resp.json()['candidates'][0]['content']['parts'][0]['text']
                # Limpeza de markdown
                clean_list = re.sub(r'```[a-z]*', '', raw_text).replace('```', '').strip()
                res = ast.literal_eval(clean_list)
                if isinstance(res, list): return res
        except:
            continue
    return lista_suja[:40]

# --- 4. RESUMO PARA BANCADA (IA) ---
def analisar_abstract_com_ia(titulo, resumo_texto, api_key):
    if not api_key: return "Chave API não configurada."
    
    prompt = f"""
    Resuma para um farmacologista experimental (máximo 15 palavras).
    Foque no Alvo (Receptor/Enzima) e na Resposta Tecidual (Contração/Relaxamento/Sinalização).
    IGNORE dados clínicos.
    
    TÍTULO: {titulo}
    TEXTO: {resumo_texto[:800]}
    """
    
    for m in MODELOS_ATIVOS:
        try:
            url = montar_url(m, api_key)
            resp = requests.post(url, json={"contents": [{"parts": [{"text": prompt}]}]}, timeout=10)
            if resp.status_code == 200:
                return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
        except:
            continue
    return "IA ocupada ou erro de conexão."

# --- 5. MINERAÇÃO DE ALVOS (PUBMED) ---
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    
    # Busca focada em artigos recentes de ciência básica (evitando reviews)
    query = f"({termo_base}) AND (2020:2026[Date - Publication]) NOT Review[pt]"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=500, sort="relevance")
        record = Entrez.read(handle); handle.close()
        
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        lines = handle.read().splitlines(); handle.close()
        
        candidatos = []
        for line in lines:
            if line.startswith("TI  - ") or line.startswith("OT  - "):
                # Regex para pegar siglas e nomes próprios (ex: TRPV1, Trehalose)
                found = re.findall(r'\b(?:[A-Z]{2,}[A-Z0-9-]*|[A-Z][a-z]{3,})\b', line)
                for f in found:
                    if f.upper() not in BLACKLIST_RADICAL and len(f) > 2:
                        candidatos.append(f)
        
        top_terms = [t for t, count in Counter(candidatos).most_common(100)]
        
        if usar_ia:
            return _faxina_ia(top_terms)
        return top_terms[:40]
    except Exception as e:
        st.error(f"Erro PubMed: {e}")
        return []

# =========================
# UI STREAMLIT
# =========================
st.title("🔬 Deep Science Prospector (Lemos Lambda)")

with st.sidebar:
    st.header("Configurações")
    st.session_state.api_key_usuario = st.text_input("Gemini API Key", type="password")
    email_user = st.text_input("Seu Email (Pubmed)", value="pesquisador@unifesp.br")

termo_busca = st.text_input("Alvo de Prospecção (ex: Bladder contractility)", value="Bladder contractility")

if st.button("Minerar Alvos e Mecanismos"):
    with st.spinner("Varrendo literatura e limpando ruído clínico..."):
        resultados = buscar_alvos_emergentes_pubmed(termo_busca, email_user)
        
        if resultados:
            st.subheader("🎯 Alvos e Fármacos Identificados")
            st.write(", ".join(resultados))
            
            
            
            st.divider()
            st.subheader("📄 Análise de Artigos Recentes")
            # Busca os resumos para análise individual
            artigos = buscar_alvos_emergentes_pubmed(termo_busca, email_user, usar_ia=False) # simplificado para exemplo
            # (Aqui você pode adicionar a lógica de exibir os abstracts e o botão de analisar com IA)
        else:
            st.warning("Nenhum alvo molecular puro encontrado.")
