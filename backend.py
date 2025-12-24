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

# --- 4. MINERAÇÃO MASSIVA (FOCO RESTRITO A TI E KW) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    
    # Passo 1: IA Proativa (Brainstorming)
    alvos_previstos = []
    if usar_ia:
        # Prompt em inglês para maior estabilidade no Paid Tier
        prompt_brain = f"As a PhD in Pharmacology, provide a Python list of 50 specific molecular targets (channels, receptors, enzymes) for research on {termo_base}. Focus on bench research (TRPV4, SPHK1, etc). No conversation, only []."
        alvos_previstos = _chamar_ia_simples(prompt_brain)
    
    # Passo 2: Busca PubMed com Filtro de Linha Real
    query = f"({termo_base} AND (Pharmacology OR Molecular)) AND (2018:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1000, sort="relevance")
        record = Entrez.read(handle); handle.close()
        total_docs = int(record["Count"])
        
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        # --- FILTRAGEM CIRÚRGICA DE LINHAS ---
        candidatos_pubmed = []
        # Lista de metadados para ignorar mesmo se aparecerem no título por erro
        metadados_sistema = {'FAU', 'PMID', 'STAT', 'OWN', 'DCOM', 'JID', 'EDAT', 'MHDA', 'CRDT', 'PHST', 'AID', 'NLM', 'USA', 'HHS', 'NIH'}
        
        for artigo in full_data.split("\n\nPMID-"):
            texto_util_do_artigo = ""
            for linha in artigo.split("\n"):
                # SÓ CONCATENA SE A LINHA COMEÇAR COM AS TAGS DE CONTEÚDO
                # TI = Título, KW/OT = Keywords (onde estão os fármacos e alvos)
                if linha.startswith("TI  - ") or linha.startswith("KW  - ") or linha.startswith("OT  - "):
                    texto_util_do_artigo += linha[6:].strip() + " "
            
            # Agora aplicamos o Regex APENAS no texto que extraímos de TI e KW
            # Ignoramos qualquer coisa com menos de 3 letras (limpa anos e siglas do sistema)
            encontrados = re.findall(r'\b[A-Z-]{3,}\b', texto_util_do_artigo)
            for t in encontrados:
                t_up = t.upper()
                if t_up not in metadados_sistema:
                    candidatos_pubmed.append(t_up)

        contagem = Counter(candidatos_pubmed)
        top_pubmed = [termo for termo, freq in contagem.most_common(150)]
        
        # Passo 3: Cruzamento via IA
        nomes_finais = []
        if usar_ia:
            lista_cruzamento = list(set(alvos_previstos + top_pubmed))
            prompt_cross = f"Cross-reference your knowledge with this list: {lista_cruzamento[:120]}. Keep only specific molecular targets and drugs for {termo_base}. Output: ONLY a Python list []."
            nomes_finais = _chamar_ia_simples(prompt_cross)
        
        # Redundância (se a IA falhar, o app não fica vazio)
        if not nomes_finais: nomes_finais = top_pubmed[:60]

        res_finais = []
        for nome in nomes_finais:
            n_limpo = limpar_termo_para_pubmed(nome)
            freq = contagem.get(n_limpo.upper(), 1)
            l_score, p_val, b_ocean, status = calcular_metricas_originais(freq, total_docs, freq * 10)
            res_finais.append({"Alvo": n_limpo, "Lambda": round(l_score, 2), "P-value": round(p_val, 4), "Blue Ocean": round(b_ocean, 1), "Status": status})
        return res_finais
    except: return []

# --- FUNÇÃO INTERNA PARA CHAMADA DE MODELOS (STABILITY) ---
def _chamar_ia_simples(prompt):
    api_key = st.session_state.get('api_key_usuario', '').strip()
    if not api_key: return []
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.2}}
    
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

# RADAR E APOIO (Mantidos conforme sua versão funcional)
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_todas_noticias(lang='pt'):
    try:
        query = "(bladder physiology OR pharmacology OR molecular biology) AND (2024:2026[Date - Publication])"
        handle = Entrez.esearch(db="pubmed", term=query, retmax=6, sort="pub_date")
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        news = []
        for art in dados.split("\n\nPMID-"):
            tit, journal, pmid = "", "", ""
            for line in art.split("\n"):
                if line.startswith("TI  - "): tit = line[6:].strip()
                if line.startswith("JT  - "): journal = line[6:].strip()
                if line.strip().isdigit() and not pmid: pmid = line.strip()
            if tit and pmid:
                news.append({"titulo": tit, "fonte": journal[:40], "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return news
    except: return []

@st.cache_data(ttl=86400, show_spinner=False)
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    termo_limpo = limpar_termo_para_pubmed(termo)
    query = f"({termo_limpo}) AND ({contexto}) AND ({ano_ini}:{ano_fim}[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle); handle.close()
        return int(record["Count"])
    except: return 0
