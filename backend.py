import streamlit as st
from Bio import Entrez
import requests
import json
from tenacity import retry, stop_after_attempt, wait_exponential
import re
from collections import Counter
import time
import ast

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br"

# --- 1. SEUS MODELOS (PAID TIER) ---
MODELOS_ATIVOS = [
    "gemini-2.5-flash",          
    "gemini-2.0-flash",          
    "gemini-2.0-flash-exp",      
    "gemini-flash-latest",       
    "gemini-2.5-flash-lite"
]

# --- 2. MAPA DE EXPANSÃO (FISIOLOGIA PURA) ---
MAPA_SINONIMOS = {
    "HEART": "Heart AND (Physiology OR Pharmacology OR Molecular)",
    "BLADDER": "Bladder AND (Physiology OR Pharmacology OR Molecular OR Detrusor OR Urothelium)",
    "KIDNEY": "Kidney AND (Physiology OR Pharmacology OR Molecular OR Nephron)",
    "BRAIN": "Brain AND (Physiology OR Pharmacology OR Molecular OR Neuron)",
}

# --- 3. A DOUTORA EM FARMACOLOGIA (FILTRO DE BANCADA) ---
def _faxina_ia(lista_suja):
    api_key = st.session_state.get('api_key_usuario', '')
    if not api_key: return lista_suja[:40] 

    lista_str = ", ".join(lista_suja)
    
    prompt = f"""
    Amiga, tu é uma doutora em farmacologia e fisiologia experiente. 
    Tu sabe que a gente quer o que é PROMISSOR para o laboratório (banho de órgãos e biologia molecular).
    
    REGRAS DE EXTERMÍNIO DA DOUTORA:
    1. DELETE SEM DÓ: Animais (Toad, Turtle, Rabbit, Mammalian), Termos de preenchimento (Studies, Role, Effect, Action, Biological, Experimental, Technique, Interpretation, Measurement, Properties, Responses, Responses, Movement, Water, Function, Pressure, Resistance, Consumption, Nervous).
    2. DELETE: Metodologia e solventes (DMSO, Vitro, Cystometrography, EMG, Technique).
    3. MANTENHA APENAS: Alvos (TRPV4, PIEZO1, SPHK1, NLRP3, Receptores, Genes) e Fármacos/Moléculas de ação (GSK1016790A, ATP, Aldosterone, Ouabain, Inibidores).
    4. Tu não é médica e nem vai escrever revisão, vacilão. Foca na bancada.
    
    LISTA PARA FILTRAR: {lista_str}
    
    OUTPUT: Retorne APENAS uma lista Python limpa. Exemplo: ['TRPV4', 'GSK1016790A', 'PIEZO1', 'SPHK1']
    """
    
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.1}}
    
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=12)
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                texto = texto.replace("```python", "").replace("```", "").strip()
                return ast.literal_eval(texto)
        except: continue
    return lista_suja[:30]

# --- 4. ANÁLISE DE RESUMOS (ALVO | FÁRMACO | EFEITO) ---
def analisar_abstract_com_ia(titulo, dados_curtos, api_key, lang='pt'):
    if not api_key: return "⚠️ Chave API necessária."
    
    idioma = "Português" if lang == 'pt' else "Inglês"
    prompt_text = f"""
    Doutora, analise esse paper com foco em farmacologia e fisiologia de lab:
    ARTIGO: {titulo}. {dados_curtos}
    
    FORMATO OBRIGATÓRIO:
    Alvo: [Sigla do Alvo] | Fármaco: [O que usaram no lab] | Efeito: [Resposta funcional observada].
    
    Máximo 20 palavras. Seja técnica e direta.
    Idioma: {idioma}.
    """
    
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt_text}]}]}
    
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=10)
            if resp.status_code == 200:
                return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
        except: continue
    return "⚠️ IA indisponível."

# --- 5. BUSCA PUBMED ---
@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
def _fetch_pubmed_count(query):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    record = Entrez.read(handle); handle.close()
    return int(record["Count"])

@st.cache_data(ttl=86400, show_spinner=False)
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    q_contexto = MAPA_SINONIMOS.get(contexto.upper(), contexto) if contexto else ""
    query = f"({termo}) AND ({q_contexto}) AND ({ano_ini}:{ano_fim}[Date - Publication])"
    try: return _fetch_pubmed_count(query)
    except: return 0

@st.cache_data(ttl=86400, show_spinner=False)
def buscar_resumos_detalhados(termo, orgao, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    q_orgao = MAPA_SINONIMOS.get(orgao.upper(), orgao)
    query = f"({termo}) AND ({q_orgao}) AND ({ano_ini}:{ano_fim}[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        artigos = []
        for raw in dados.split("\n\nPMID-"):
            tit, pmid, keywords, abstract = "", "", "", ""
            for line in raw.split("\n"):
                if line.startswith("TI  - "): tit = line[6:].strip()
                if line.startswith("AB  - "): abstract = line[6:500].strip()
                if line.startswith("OT  - ") or line.startswith("KW  - "): keywords += line[6:].strip() + ", "
            if tit:
                artigos.append({"Title": tit, "Info_IA": f"{keywords} {abstract}", "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return artigos
    except: return []

# --- 6. MINERAÇÃO MASSIVA ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    termo_up = termo_base.upper().strip()
    query = f"({MAPA_SINONIMOS.get(termo_up, termo_base)}) AND (2018:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1500, sort="relevance")
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        # Blacklist Pesada de Termos Irrelevantes
        blacklist = {
            "THE", "AND", "ROLE", "EFFECT", "ACTION", "STUDIES", "STUDY", "DATA", "RESULTS",
            "WATER", "PRESSURE", "RESISTANCE", "MOVEMENT", "CONSUMPTION", "EXPERIMENTAL",
            "BIOLOGICAL", "PROPERTIES", "MEASUREMENT", "INTERPRETATION", "TECHNIQUE", "RESPONSE",
            "RESPONSES", "VITRO", "VIVO", "MODEL", "MODELS", "RABBIT", "TURTLE", "MAMMALIAN",
            "TOAD", "MOUSE", "MICE", "RAT", "RATS", "DMSO", "PBS", "SALINE"
        }
        
        candidatos = []
        for artigo in full_data.split("\n\nPMID-"):
            texto = ""
            for line in artigo.split("\n"):
                if line.startswith("TI  - ") or line.startswith("KW  - "): texto += line[6:].strip() + " "
            encontrados = re.findall(r'\b(?:[A-Z]{2,}[A-Z0-9-]*|[a-z]{1,2}[A-Z][a-zA-Z0-9-]*)\b', texto)
            for t in encontrados:
                t_clean = re.sub(r'[^A-Z0-9]', '', t.upper())
                if len(t_clean) >= 3 and t_clean not in blacklist: candidatos.append(t_clean)

        contagem = Counter(candidatos)
        top = [termo for termo,freq in contagem.most_common(180)]
        
        return _faxina_ia(top) if usar_ia and st.session_state.get('api_key_usuario') else top[:30]
    except: return []
