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

    # LISTA ATUALIZADA COM SEUS MODELOS (Prioridade: 2.5 > 2.0 > Latest)
    modelos_disponiveis = [
        "gemini-2.5-flash",          # Prioridade máxima (Rápido e Novo)
        "gemini-2.0-flash",          # Backup estável
        "gemini-2.0-flash-exp",      # Experimental
        "gemini-flash-latest",       # Genérico
        "gemini-2.5-flash-lite"      # Super leve se tudo falhar
    ]

    for modelo in modelos_disponiveis:
        try:
            # URL Direta v1beta (Padrão para modelos novos)
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={api_key}"
            
            response = requests.post(url, headers=headers, data=json.dumps(data), timeout=10)
            
            if response.status_code == 200:
                resultado = response.json()
                try:
                    texto = resultado['candidates'][0]['content']['parts'][0]['text']
                    return texto.strip()
                except:
                    return "⚠️ IA respondeu vazio (Erro no parse)."
            
            elif response.status_code == 429:
                return "💡 IA Ocupada (Cota excedida). Aguarde..."
            
            elif response.status_code == 400 or response.status_code == 403:
                # Se a chave for inválida, para de tentar e avisa logo
                return "❌ Chave API Inválida/Permissão Negada."
            
            # Se der 404, apenas continua o loop para o próximo modelo
            continue

        except Exception:
            continue

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
    try:
        return _fetch_pubmed_count(query)
    except:
        return 0

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
                if line.startswith("OT  - ") or line.startswith("KW  - "):
                    keywords += line[6:].strip() + ", "
                if line.startswith("AB  - ") and not fallback_text:
                    fallback_text = line[6:500].strip()
            
            if tit:
                info_final = keywords if len(keywords) > 5 else fallback_text
                artigos.append({
                    "Title": tit, 
                    "Info_IA": info_final if info_final else "Sem resumo disponível.",
                    "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                })
        return artigos
    except: return []

# --- MOTOR DE MINERAÇÃO (CALIBRADO V2.2 - Janela Estendida) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email):
    if email: Entrez.email = email
    
    # 1. DATA AMPLIADA PARA 2016 (Janela de ~10 anos)
    query = f"({termo_base}[Title/Abstract]) AND (2016:2030[Date - Publication]) AND (NOT Review[pt])"
    
    try:
        # Mantive retmax 100 para ter boa amostragem
        handle = Entrez.esearch(db="pubmed", term=query, retmax=100, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        artigos_raw = full_data.split("\n\nPMID-")

        keywords_funcionais = {"SIGNALING", "PATHWAY", "RECEPTOR", "CHANNEL", "ACTIVATION", "INHIBITION", "EXPRESSION",
                               "REGULATION", "MEDIATED", "MECHANISM", "FUNCTION", "ROLE", "TARGET", "MOLECULAR", "GENE",
                               "PROTEIN", "ENZYME", "KINASE", "AUTOPHAGY", "APOPTOSIS", "DIFFERENTIATION", "HOMEOSTASIS", 
                               "STRESS", "CONTRACTILITY", "RELAXATION"}
        
        # Blacklist Expandida Mantida
        blacklist = {"AND", "THE", "FOR", "NOT", "BUT", "WITH", "FROM", "STUDY", "RESULTS", "CELLS", "WAS", "WERE", 
                     "CBS", "FAU", "AID", "USE", "AIM", "DUE", "VIA", "MAY", "CAN", "HAS", "HAD", "LOW", "HIGH",
                     "NEW", "TWO", "ONE", "RAT", "MICE", "DOG", "PIG", "CAT", "MAN", "AGE", "OLD", "DRUG", "USED",
                     "AFTER", "BEFORE", "DURING", "WITHIN", "BETWEEN", "AMONG", "UNDER", "ABOVE", "BELOW", "THESE",
                     "THOSE", "THIS", "THAT", "WHICH", "WHAT", "WHEN", "WHERE", "WHY", "HOW", "ALL", "ANY", "EACH",
                     "SIGNIFICANTLY", "RESPECTIVELY", "CONCLUSION", "BACKGROUND", "METHODS", "OBJECTIVE", "PURPOSE",
                     "DEPARTMENT", "UNIVERSITY", "HOSPITAL", "CENTER", "CENTRE", "SCHOOL", "COLLEGE", "INSTITUTE",
                     "RESEARCH", "SCIENCE", "SCIENCES", "MEDICINE", "MEDICAL", "CLINICAL", "PATIENTS", "GROUP", "GROUPS",
                     "CONTROL", "TREATED", "UNTREATED", "SHAM", "VEHICLE", "MODEL", "MODELS", "DATA", "ANALYSIS",
                     "USING", "USED", "PERFORMED", "OBSERVED", "SHOWED", "FOUND", "INDICATED", "SUGGESTED", "DEMONSTRATED",
                     "INVESTIGATED", "EVALUATED", "COMPARED", "INCREASED", "DECREASED", "LEVELS", "EFFECTS"} 
        
        unidades = {"MMHG","KPA","MIN","SEC","HRS","ML","MG","KG","NM","UM","MM","NMOL", "MOL", "PH"}

        candidatos_por_artigo = []
        for artigo in artigos_raw:
            texto_upper = artigo.upper()
            if not any(kw in texto_upper for kw in keywords_funcionais): continue

            # Regex 15 chars mantido
            encontrados = re.findall(r'\b[A-Z][A-Z0-9-]{2,14}\b', texto_upper)
            
            candidatos_locais = set()
            for t in encontrados:
                t_clean = re.sub(r'[^A-Z0-9]', '', t)
                
                if t_clean in blacklist or t_clean in unidades: continue
                if len(t_clean) < 3: continue 
                if t_clean == termo_base.upper().replace(" ", ""): continue
                if t_clean.isdigit(): continue 
                
                candidatos_locais.add(t_clean)
            candidatos_por_artigo.extend(list(candidatos_locais))

        if not candidatos_por_artigo: return []
        contagem = Counter(candidatos_por_artigo)
        total_docs = max(1, len(artigos_raw))
        
        # Filtro de corte mantido
        return [termo for termo,freq in contagem.most_common(15) if (freq/total_docs)<0.50][:10]
    except: return []

# --- RADAR DE NOTÍCIAS ---
def buscar_todas_noticias(lang='pt'):
    try:
        query = "(molecular biology OR pharmacology OR drug discovery) AND (2024/09/01:2025/12/31[Date - Publication])"
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
                news.append({"titulo": tit, "fonte": journal[:30], 
                             "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                             "img":"https://images.unsplash.com/photo-1532187863486-abf9dbad1b69?w=400"})
        return news
    except: return []
