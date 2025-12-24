import streamlit as st
from Bio import Entrez
import requests
import json
from tenacity import retry, stop_after_attempt, wait_exponential
import re
from collections import Counter
import ast
import math

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br"

MODELOS_ATIVOS = [
    "gemini-2.5-flash", "gemini-2.0-flash", "gemini-2.0-flash-exp", 
    "gemini-flash-latest", "gemini-2.5-flash-lite"
]

# --- 1. MOTOR DE MÉTRICAS ---
def calcular_metricas_originais(freq, total_docs, n_alvo_total):
    lambda_score = (freq / total_docs) * 100 if total_docs > 0 else 0
    p_value = math.exp(-freq/5) if freq > 0 else 1.0
    blue_ocean = max(0, 100 - (n_alvo_total / 10)) if n_alvo_total > 0 else 100.0
    status = "Saturado" if blue_ocean < 25 else "Blue Ocean" if blue_ocean > 75 else "Competitivo"
    return lambda_score, p_value, blue_ocean, status

# --- 2. HIGIENIZAÇÃO RÍGIDA ---
def limpar_termo_para_pubmed(termo):
    # Remove qualquer caractere que não seja letra ou número no início/fim
    return re.sub(r'^[^a-zA-Z0-9]+|[^a-zA-Z0-9]+$', '', termo.strip())

# --- 3. A DOUTORA INVESTIGADORA (COM PENEIRA DIAMANTE) ---
def _doutora_investigadora(termo_base, lista_pubmed=None, fase="brainstorming"):
    api_key = st.session_state.get('api_key_usuario', '')
    if not api_key: return []

    if fase == "brainstorming":
        prompt = f"""
        Tu é uma doutora em farmacologia molecular de elite. 
        MISSÃO: Listar 60 alvos, fármacos e sinalizadores que são a "fronteira" para {termo_base}.
        FOCO: Alvos de bancada (Canais, Receptores, Enzimas).
        EXEMPLOS: TRPV4, GSK1016790A, SPHK1, NLRP3, miR-132, SNARE, PIEZO1, P2X3.
        OUTPUT: Retorne APENAS a lista Python de strings.
        """
    else:
        lista_str = ", ".join(lista_pubmed)
        prompt = f"""
        Tu é a doutora investigadora. 
        CRUZAMENTO: Combine suas previsões com estes dados do PubMed: {lista_str}.
        
        REGRAS DE EXTERMÍNIO (STRICT):
        - Delete termos de 1 ou 2 letras (O, J, II, IV, S).
        - Delete termos genéricos: (Study, Clinical, Results, During, UK, 2024, 2025, New, OF, IN).
        - Delete termos administrativos/clínicos: (NSQIP, GA-CARES, GWAS, GENIE, TCGA).
        - Mantenha apenas o OURO farmacológico (ex: GSK1016790A, Mirabegron, TRPV4).
        
        OUTPUT: Retorne APENAS a lista Python final (60-80 itens).
        """

    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.3}}
    
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=20)
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                texto = texto.replace("```python", "").replace("```", "").strip()
                return ast.literal_eval(texto)
        except: continue
    return []

# --- 4. RADAR DE NOTÍCIAS (RESTAURADO) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_todas_noticias(lang='pt'):
    try:
        query = "(molecular biology OR pharmacology OR bladder physiology) AND (2024:2026[Date - Publication])"
        handle = Entrez.esearch(db="pubmed", term=query, retmax=6, sort="pub_date")
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        news = []
        for art in dados.split("\n\nPMID-"):
            tit, pmid, journal = "", "", ""
            for line in art.split("\n"):
                if line.startswith("TI  - "): tit = line[6:].strip()
                if line.startswith("JT  - "): journal = line[3:].strip()
                if line.strip().isdigit() and not pmid: pmid = line.strip()
            if tit and pmid:
                news.append({
                    "titulo": tit, "fonte": journal[:40], 
                    "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                })
        return news
    except: return []

# --- 5. MINERAÇÃO E CRUZAMENTO ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    
    # PASSO 1: Brainstorming
    alvos_previstos = []
    if usar_ia: alvos_previstos = _doutora_investigadora(termo_base, fase="brainstorming")
    
    # PASSO 2: Busca PubMed
    query = f"({termo_base} AND (Pharmacology OR Molecular)) AND (2018:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=2000, sort="relevance")
        record = Entrez.read(handle); handle.close()
        total_docs = int(record["Count"])
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        # Blacklist de extermínio estatístico
        blacklist_hard = {"THE", "AND", "WITH", "ROLE", "EFFECT", "STUDY", "RESULTS", "DURING", "TOAD", "RABBIT", "NEW", "III", "IV"}
        
        candidatos_pubmed = []
        for artigo in full_data.split("\n\nPMID-"):
            texto = ""
            for line in artigo.split("\n"):
                if line.startswith("TI  - ") or line.startswith("KW  - "): texto += line[6:].strip() + " "
            encontrados = re.findall(r'\b(?:[A-Z0-9-]{3,}|[a-z]{1,2}-[A-Z0-9]{2,}|[a-z]{1,2}[A-Z][a-zA-Z0-9-]*)\b', texto)
            candidatos_pubmed.extend([t.upper() for t in encontrados if t.upper() not in blacklist_hard])

        contagem = Counter(candidatos_pubmed)
        top_pubmed = [termo for termo, freq in contagem.most_common(250)]
        
        # PASSO 3: Cruzamento Final
        nomes_finais = []
        if usar_ia:
            lista_para_cruzamento = list(set(alvos_previstos + top_pubmed))
            nomes_finais = _doutora_investigadora(termo_base, lista_pubmed=lista_para_cruzamento, fase="cruzamento")
        else:
            nomes_finais = top_pubmed[:40]

        resultados_finais = []
        for nome in nomes_finais:
            nome_limpo = limpar_termo_para_pubmed(nome)
            freq = contagem.get(nome_limpo.upper(), 1)
            l_score, p_val, b_ocean, status = calcular_metricas_originais(freq, total_docs, freq * 10)
            resultados_finais.append({
                "Alvo": nome_limpo, "Lambda": round(l_score, 2), "P-value": round(p_val, 4),
                "Blue Ocean": round(b_ocean, 1), "Status": status
            })
        return resultados_finais
    except: return []

# --- 6. ANÁLISE DE ABSTRACT ---
def analisar_abstract_com_ia(titulo, dados_curtos, api_key, lang='pt'):
    if not api_key: return "⚠️ Chave necessária."
    prompt_text = f"Analise como PhD em farmacologia: {titulo}. {dados_curtos}. Formato: Alvo: [Sigla] | Fármaco: [O que usaram] | Efeito: [Resposta funcional]. Max 20 palavras."
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt_text}]}]}
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=12)
            if resp.status_code == 200:
                return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
        except: continue
    return "⚠️ Erro IA."
