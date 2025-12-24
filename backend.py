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

# --- 3. A DOUTORA EM FARMACOLOGIA (FILTRO COM REFERÊNCIA) ---
def _faxina_ia(lista_suja):
    api_key = st.session_state.get('api_key_usuario', '')
    if not api_key: return lista_suja[:50] 

    lista_str = ", ".join(lista_suja)
    
    prompt = f"""
    Amiga, tu é uma doutora em farmacologia e fisiologia experiente. 
    Tua missão é filtrar essa lista e deixar apenas o que é ÚTIL para quem faz pesquisa de bancada (molécula, fármaco, alvo).
    
    COMO IDENTIFICAR UM ALVO (MEUS EXEMPLOS):
    - Canais/Receptores: TRPV4, PIEZO1, P2X3, TRPM8, TRPV1.
    - Fármacos/Inibidores: GSK1016790A, Ouabain, Aldosterone, O-1602.
    - Sinalização/Enzimas: SPHK1, NLRP3, SNARE, MAPK, PI3K.
    - Reguladores: miR-132, GATA3, PPAR.
    - Metabólitos: TMAO.

    REGRA DE OURO: 
    Usa esses exemplos para entender o PADRÃO (siglas técnicas, códigos alfanuméricos), mas caça TUDO que for parecido com isso na lista. Não fica presa só nos exemplos!
    
    O QUE DELETAR (LIXO):
    1. Termos de "prosa" acadêmica: (Role, Effect, Action, Studies, Response, During, After, High, Low).
    2. Termos médicos/clínicos: (Incontinence, Diabetes, Patient, Treatment, Symptoms, Survey).
    3. Metodologia: (Cystometry, EMG, Questionnaire, Measurement, Technique).
    4. Animais e termos comuns: (Rabbit, Toad, Water, Female, Male, Body).
    
    LISTA PARA FILTRAR: {lista_str}
    
    OUTPUT: Retorne APENAS a lista Python limpa. Exemplo: ['TRPV4', 'ALVO_NOVO', 'SIGLA_TECNICA']
    """
    
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.2}}
    
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=15)
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                texto = texto.replace("```python", "").replace("```", "").strip()
                lista_limpa = ast.literal_eval(texto)
                if isinstance(lista_limpa, list): return lista_limpa
        except: continue
    return lista_suja[:40]

# --- 4. ANÁLISE DE RESUMOS ---
def analisar_abstract_com_ia(titulo, dados_curtos, api_key, lang='pt'):
    if not api_key: return "⚠️ IA não ativada."
    idioma = "Português" if lang == 'pt' else "Inglês"
    prompt_text = f"Doutora, analise: {titulo}. {dados_curtos}. FORMATO: Alvo: [Sigla] | Fármaco: [Agentes] | Efeito: [Resposta]. Max 20 palavras. Idioma: {idioma}."
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt_text}]}]}
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=12)
            if resp.status_code == 200:
                return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
        except: continue
    return "⚠️ IA indisponível."

# --- 5. BUSCA PUBMED ---
@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=1, min=4, max=20))
def _fetch_pubmed_count(query):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    record = Entrez.read(handle); handle.close()
    return int(record["Count"])

@st.cache_data(ttl=86400, show_spinner=False)
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    Entrez.email = email if email else "pesquisador_guest@unifesp.br"
    q_contexto = MAPA_SINONIMOS.get(contexto.upper(), contexto) if contexto else ""
    query = f"({termo})"
    if q_contexto: query += f" AND ({q_contexto})"
    query += f" AND ({ano_ini}:{ano_fim}[Date - Publication])"
    try: return _fetch_pubmed_count(query)
    except: return 0

# --- 6. MINERAÇÃO MASSIVA ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    Entrez.email = email if email else "pesquisador_guest@unifesp.br"
    termo_up = termo_base.upper().strip()
    query = f"({MAPA_SINONIMOS.get(termo_up, termo_base)}) AND (2018:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1500, sort="relevance")
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        candidatos = []
        for artigo in full_data.split("\n\nPMID-"):
            texto = ""
            for line in artigo.split("\n"):
                if line.startswith("TI  - ") or line.startswith("KW  - "): texto += line[6:].strip() + " "
            encontrados = re.findall(r'\b(?:[A-Z0-9-]{2,}|[a-z]{1,2}-[A-Z0-9]{2,}|[a-z]{1,2}[A-Z][a-zA-Z0-9-]*)\b', texto)
            for t in encontrados:
                t_clean = t.strip("-").upper()
                if len(t_clean) >= 3: candidatos.append(t_clean)
        contagem = Counter(candidatos)
        top = [termo for termo,freq in contagem.most_common(200)]
        if usar_ia and st.session_state.get('api_key_usuario'):
            return _faxina_ia(top)
        else: return top[:40]
    except: return []
        
