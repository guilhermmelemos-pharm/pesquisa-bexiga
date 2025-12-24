import streamlit as st
from Bio import Entrez
import requests
import json
from tenacity import retry, stop_after_attempt, wait_exponential
import re
from collections import Counter
import time

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br"

# --- IA: CONEXÃO DIRETA (MODELOS 2025) ---
def analisar_abstract_com_ia(titulo, dados_curtos, api_key, lang='pt'):
    if not api_key:
        return "⚠️ IA não ativada"
    
    idioma = "Português" if lang == 'pt' else "Inglês"
    
    prompt_text = f"""PhD em Farmacologia, analise:
FONTE: {titulo}. {dados_curtos}
FORMATO: Alvo: [Sigla] | Fármaco: [Nome] | Efeito: [Ação funcional].
REGRAS: Máximo 12 palavras. Seja técnico. Idioma: {idioma}."""

    # Headers e Config JSON
    headers = {'Content-Type': 'application/json'}
    data = {
        "contents": [{"parts": [{"text": prompt_text}]}],
        "safetySettings": [
            {"category": "HARM_CATEGORY_HARASSMENT", "threshold": "BLOCK_NONE"},
            {"category": "HARM_CATEGORY_HATE_SPEECH", "threshold": "BLOCK_NONE"},
            {"category": "HARM_CATEGORY_SEXUALLY_EXPLICIT", "threshold": "BLOCK_NONE"},
            {"category": "HARM_CATEGORY_DANGEROUS_CONTENT", "threshold": "BLOCK_NONE"}
        ],
        "generationConfig": {"temperature": 0.1}
    }

    modelos_disponiveis = [
        "gemini-2.5-flash", "gemini-2.0-flash", "gemini-2.0-flash-exp", 
        "gemini-flash-latest", "gemini-2.5-flash-lite"
    ]

    for modelo in modelos_disponiveis:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={api_key}"
            response = requests.post(url, headers=headers, data=json.dumps(data), timeout=10)
            if response.status_code == 200:
                resultado = response.json()
                try:
                    texto = resultado['candidates'][0]['content']['parts'][0]['text']
                    return texto.strip()
                except: return "⚠️ IA respondeu vazio."
            elif response.status_code == 429: return "💡 IA Ocupada. Aguarde..."
            elif response.status_code in [400, 403]: return "❌ Chave Inválida."
            continue
        except: continue
    return f"❌ Falha técnica. Título: {titulo[:30]}..."

# --- BUSCA PUBMED ---
@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
def _fetch_pubmed_count(query):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    record = Entrez.read(handle)
    handle.close()
    return int(record["Count"])

@st.cache_data(ttl=86400, show_spinner=False)
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    query = f"({termo})"
    if contexto: query += f" AND ({contexto})"
    query += f" AND ({ano_ini}:{ano_fim}[Date - Publication]) AND (NOT Review[pt])"
    try: return _fetch_pubmed_count(query)
    except: return 0

@st.cache_data(ttl=86400, show_spinner=False)
def buscar_resumos_detalhados(termo, orgao, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    query = f"({termo}[Title/Abstract]) AND ({orgao}[Title/Abstract]) AND ({ano_ini}:{ano_fim}[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        artigos = []
        for raw in dados.split("\n\nPMID-"):
            tit, pmid, keywords, fallback_text = "", "", "", ""
            lines = raw.split("\n")
            for line in lines:
                if line.strip().isdigit() and not pmid: pmid = line.strip()
                if line.startswith("TI  - "): tit = line[6:].strip()
                if line.startswith("OT  - ") or line.startswith("KW  - "): keywords += line[6:].strip() + ", "
                if line.startswith("AB  - ") and not fallback_text: fallback_text = line[6:500].strip()
            if tit:
                info_final = keywords if len(keywords) > 5 else fallback_text
                artigos.append({"Title": tit, "Info_IA": info_final if info_final else "Sem resumo.", "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return artigos
    except: return []

# --- MOTOR DE MINERAÇÃO (CORRIGIDO: FILTRO DE CONTEÚDO) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email):
    if email: Entrez.email = email
    
    # Busca agressiva: 200 artigos recentes
    query = f"({termo_base}[Title/Abstract]) AND (2018:2030[Date - Publication]) AND (NOT Review[pt])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=200, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        artigos_raw = full_data.split("\n\nPMID-")

        # NOVA BLACKLIST DE METADADOS (Remove o lixo que você encontrou)
        blacklist = {
            "PII", "DOI", "ISSN", "PMID", "PPUBLISH", "ACCEPTED", "RECEIVED", "REVISED", "COPYRIGHT", 
            "AUTHOR", "AUTHORS", "EDITOR", "EMAIL", "CORRESPONDENCE", "DEPARTMENT", "UNIVERSITY", 
            "PUBMED", "MEDLINE", "PMC", "NIHMS", "ARTICLE", "JOURNAL", "UROLOGY", "NEPHROLOGY", 
            "PHARMACOLOGY", "SCIENCE", "REPORTS", "NATURE", "CELL", "PLOS", "ONE", "FRONTIERS",
            "AND", "THE", "FOR", "NOT", "BUT", "WITH", "FROM", "STUDY", "RESULTS", "CELLS", "WAS", "WERE",
            "USING", "USED", "DATA", "ANALYSIS", "GROUP", "CONTROL", "MODEL", "RATS", "MICE", "HUMAN",
            "PATIENTS", "CLINICAL", "TREATMENT", "EFFECTS", "LEVELS", "EXPRESSION", "INCREASED", "DECREASED",
            "SIGNIFICANTLY", "CONCLUSION", "BACKGROUND", "METHODS", "OBJECTIVE", "AIM", "PURPOSE",
            "AFTER", "BEFORE", "DURING", "WITHIN", "BETWEEN", "AMONG", "UNDER", "ABOVE", "BELOW",
            "ACUTE", "CHRONIC", "DISEASE", "SYNDROME", "DISORDER", "INJURY", "HEALTH", "CARE", "RISK",
            "TOTAL", "MEAN", "RATIO", "SCORE", "FACTOR", "TYPE", "CASE", "SERIES", "TRIAL", "REVIEW"
        } 
        
        candidatos_por_artigo = []
        for artigo in artigos_raw:
            # --- MUDANÇA CRÍTICA: LÊ APENAS TÍTULO (TI) E ABSTRACT (AB) ---
            # Isso evita ler datas, códigos de publicação e afiliações
            linhas_uteis = []
            lendo_abstract = False
            
            for line in artigo.split("\n"):
                if line.startswith("TI  - "):
                    linhas_uteis.append(line[6:].strip())
                    lendo_abstract = False
                elif line.startswith("AB  - "):
                    linhas_uteis.append(line[6:].strip())
                    lendo_abstract = True
                elif lendo_abstract and line.startswith("      "): # Continuação do abstract
                    linhas_uteis.append(line.strip())
                elif line.startswith(" ") is False: # Nova tag começou (ex: FAU -)
                    lendo_abstract = False

            texto_focado = " ".join(linhas_uteis).upper()
            
            # Extrai siglas apenas do texto focado (Titulo + Abstract)
            encontrados = re.findall(r'\b[A-Z][A-Z0-9-]{2,14}\b', texto_focado)
            
            candidatos_locais = set()
            for t in encontrados:
                t_clean = re.sub(r'[^A-Z0-9]', '', t)
                
                # Filtros Rigorosos
                if t_clean in blacklist: continue
                if len(t_clean) < 3: continue 
                if t_clean == termo_base.upper().replace(" ", ""): continue
                if t_clean.isdigit(): continue 
                
                candidatos_locais.add(t_clean)
            candidatos_por_artigo.extend(list(candidatos_locais))

        if not candidatos_por_artigo: return []
        
        contagem = Counter(candidatos_por_artigo)
        total_docs = max(1, len(artigos_raw))
        
        # Filtro final: Top 50 candidatos, frequência < 80%
        return [termo for termo,freq in contagem.most_common(50) if (freq/total_docs)<0.80][:10]
    except: return []

# --- RADAR ---
def buscar_todas_noticias(lang='pt'):
    try:
        query = "(molecular biology OR pharmacology) AND (2024/09/01:2025/12/31[Date - Publication])"
        handle = Entrez.esearch(db="pubmed", term=query, retmax=6, sort="pub_date")
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        news = []
        for art in dados.split("\n\nPMID-"):
            tit, pmid, journal = "", "", ""
            for line in art.split("\n"):
                if re.match(r'^TI\s+-', line): tit = re.sub(r'^TI\s+-\s+', '', line).strip()
                if line.startswith("JT  - "): journal = line.replace("JT  - ", "").strip()
                if line.strip().isdigit() and not pmid: pmid = line.strip()
            if tit and pmid:
                news.append({"titulo": tit, "fonte": journal[:30], "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/", "img":"https://images.unsplash.com/photo-1532187863486-abf9dbad1b69?w=400"})
        return news
    except: return []
            
