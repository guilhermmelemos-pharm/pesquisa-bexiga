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

# --- 1. MAPA DE EXPANSÃO SEMÂNTICA (O "CÉREBRO" ANATÔMICO) ---
# Isso expande a busca para termos técnicos da área automaticamente
MAPA_SINONIMOS = {
    "HEART": "Heart OR Cardiac OR Myocardium OR Cardiomyocyte OR Coronary OR Artery OR Blood Pressure OR Hypertension",
    "BLADDER": "Bladder OR Urothelium OR Detrusor OR Vesical OR Urethra OR Micturition OR LUTS OR Cystitis",
    "KIDNEY": "Kidney OR Renal OR Nephron OR Glomerulus OR Tubular OR Podocyte",
    "BRAIN": "Brain OR CNS OR Neuron OR Glia OR Cortex OR Hippocampus OR Synaptic OR Neurotransmitter",
    "LIVER": "Liver OR Hepatic OR Hepatocyte OR Steatosis OR Fibrosis",
    "LUNG": "Lung OR Pulmonary OR Alveolar OR Bronchial OR Respiratory OR Asthma",
    "INTESTINE": "Intestine OR Gut OR Colon OR Bowel OR Enteric OR Colitis",
    "PAIN": "Pain OR Nociception OR Analgesia OR Neuropathic OR Hyperalgesia OR Dorsal Root Ganglion",
    "INFLAMMATION": "Inflammation OR Cytokine OR Macrophage OR Neutrophil OR Immune OR Sepsis",
    "METABOLISM": "Metabolism OR Obesity OR Diabetes OR Insulin OR Glucose OR Adipose",
    "CANCER": "Cancer OR Tumor OR Oncology OR Carcinoma OR Metastasis OR Proliferation"
}

# --- IA: CONEXÃO DIRETA (MODELOS 2025) ---
def analisar_abstract_com_ia(titulo, dados_curtos, api_key, lang='pt'):
    if not api_key: return "⚠️ IA não ativada"
    idioma = "Português" if lang == 'pt' else "Inglês"
    prompt_text = f"""PhD em Farmacologia, analise:
FONTE: {titulo}. {dados_curtos}
FORMATO: Alvo: [Sigla] | Fármaco: [Nome] | Efeito: [Ação funcional].
REGRAS: Máximo 12 palavras. Seja técnico. Idioma: {idioma}."""
    
    headers = {'Content-Type': 'application/json'}
    data = {
        "contents": [{"parts": [{"text": prompt_text}]}],
        "safetySettings": [{"category": "HARM_CATEGORY_DANGEROUS_CONTENT", "threshold": "BLOCK_NONE"}],
        "generationConfig": {"temperature": 0.1}
    }
    modelos = ["gemini-2.5-flash", "gemini-2.0-flash", "gemini-flash-latest"]
    
    for m in modelos:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=8)
            if resp.status_code == 200:
                return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
        except: continue
    return f"❌ Falha técnica. Título: {titulo[:30]}..."

# --- BUSCA PUBMED COUNTS (Para Estatística) ---
@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
def _fetch_pubmed_count(query):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    record = Entrez.read(handle)
    handle.close()
    return int(record["Count"])

@st.cache_data(ttl=86400, show_spinner=False)
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    # Aqui usamos o termo expandido se houver contexto
    q_contexto = MAPA_SINONIMOS.get(contexto.upper(), contexto) if contexto else ""
    
    query = f"({termo})"
    if q_contexto: query += f" AND ({q_contexto})"
    query += f" AND ({ano_ini}:{ano_fim}[Date - Publication]) AND (NOT Review[pt])"
    try: return _fetch_pubmed_count(query)
    except: return 0

# --- BUSCA DE ARTIGOS PARA LEITURA ---
@st.cache_data(ttl=86400, show_spinner=False)
def buscar_resumos_detalhados(termo, orgao, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    # Expande o órgão também na busca de leitura
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
                if line.strip().isdigit() and not pmid: pmid = line.strip()
                if line.startswith("TI  - "): tit = line[6:].strip()
                if line.startswith("AB  - "): abstract = line[6:500].strip() # Aqui lemos abstract só para exibir pro usuário
                if line.startswith("OT  - ") or line.startswith("KW  - "): keywords += line[6:].strip() + ", "
            
            if tit:
                artigos.append({
                    "Title": tit, 
                    "Info_IA": f"{keywords} {abstract}", # Para IA ler, mandamos tudo
                    "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                })
        return artigos
    except: return []

# --- MOTOR DE MINERAÇÃO: ALTA PRECISÃO (TITULO + KEYWORDS) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email):
    if email: Entrez.email = email
    
    # 1. VERIFICA SE EXISTE EXPANSÃO NO DICIONÁRIO
    termo_upper = termo_base.upper().strip()
    query_string = MAPA_SINONIMOS.get(termo_upper, f"{termo_base}[Title/Abstract]")
    
    # 2. BUSCA NO PUBMED (Janela 2018-2030)
    final_query = f"({query_string}) AND (2018:2030[Date - Publication]) AND (NOT Review[pt])"
    
    try:
        # Buscamos 300 artigos para ter volume, já que vamos descartar o Abstract
        handle = Entrez.esearch(db="pubmed", term=final_query, retmax=300, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        artigos_raw = full_data.split("\n\nPMID-")

        # 3. BLACKLIST LEVE (Como olhamos só titulo, a sujeira é menor)
        blacklist = {
            "PII", "DOI", "ISSN", "PMID", "RATS", "MICE", "HUMAN", "PATIENTS", "CLINICAL", "TRIAL",
            "REVIEW", "META", "ANALYSIS", "REPORT", "CASE", "SERIES", "STUDY", "RESULTS", "CONCLUSION",
            "BACKGROUND", "METHODS", "OBJECTIVE", "AIM", "PURPOSE", "AND", "THE", "FOR", "WITH", "NOT",
            "BUT", "FROM", "USING", "USED", "BETWEEN", "AMONG", "DURING", "AFTER", "BEFORE", "WITHIN",
            "ACUTE", "CHRONIC", "DISEASE", "SYNDROME", "DISORDER", "INJURY", "FAILURE", "DAMAGE",
            "TREATMENT", "THERAPY", "MANAGEMENT", "DIAGNOSIS", "PROGNOSIS", "OUTCOME", "RISK", "FACTOR",
            "ROLE", "EFFECT", "EFFECTS", "ACTION", "FUNCTION", "EXPRESSION", "REGULATION", "LEVELS",
            "NEW", "NOVEL", "RECENT", "INSIGHTS", "UPDATE", "PERSPECTIVE", "FUTURE", "CHALLENGES",
            "PATHWAY", "MECHANISM", "SIGNALING", "TARGET", "BIOMARKER", "POTENTIAL", "CANDIDATE",
            "ASSOCIATED", "MEDIATED", "INDUCED", "DEPENDENT", "INDEPENDENT", "RELATED", "LINKED",
            "UNIVERSITY", "DEPARTMENT", "HOSPITAL", "CENTER", "GROUP", "SOCIETY", "ASSOCIATION"
        }
        
        candidatos_por_artigo = []
        
        for artigo in artigos_raw:
            # --- ESTRATÉGIA "SNIPER": SÓ TI (TÍTULO) E OT/KW (KEYWORDS) ---
            # Ignoramos AB (Abstract) completamente para evitar ruído
            texto_focado = ""
            for line in artigo.split("\n"):
                if line.startswith("TI  - "):
                    texto_focado += line[6:].strip() + " "
                elif line.startswith("OT  - ") or line.startswith("KW  - "):
                    texto_focado += line[6:].strip() + " "
            
            # Se não achou título nem keywords, pula
            if not texto_focado: continue

            # Regex: Siglas com pelo menos 2 letras maiúsculas (TRPV1, mTOR, BDNF)
            # ou mistura de letra+número (P2X3, Nav1.8)
            encontrados = re.findall(r'\b(?:[A-Z]{2,}[A-Z0-9-]*|[a-z]*[A-Z][a-zA-Z0-9-]*[0-9]+[a-zA-Z0-9-]*)\b', texto_focado)
            
            for t in encontrados:
                t_clean = re.sub(r'[^A-Z0-9]', '', t).upper()
                
                if len(t_clean) < 3: continue
                if t_clean in blacklist: continue
                if t_clean == termo_upper.replace(" ", ""): continue # Remove o termo de busca (ex: HEART)
                if t_clean.isdigit(): continue
                
                candidatos_por_artigo.append(t_clean)

        if not candidatos_por_artigo: return []
        
        contagem = Counter(candidatos_por_artigo)
        total_docs = max(1, len(artigos_raw))
        
        # Filtro: Retorna Top 10 que apareçam em menos de 90% dos títulos (muito permissivo, pois título é curto)
        return [termo for termo,freq in contagem.most_common(50) if (freq/total_docs)<0.90][:10]

    except: return []

# --- RADAR DE NOTÍCIAS ---
def buscar_todas_noticias(lang='pt'):
    try:
        # Usa o termo base 'pharmacology' mas podemos expandir no futuro
        query = "(molecular biology OR pharmacology) AND (2024/09/01:2025/12/31[Date - Publication])"
        handle = Entrez.esearch(db="pubmed", term=query, retmax=6, sort="pub_date")
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        news = []
        for art in dados.split("\n\nPMID-"):
            tit, pmid, journal = "", "", ""
            for line in art.split("\n"):
                if line.startswith("TI  - "): tit = line[6:].strip()
                if line.startswith("JT  - "): journal = line.replace("JT  - ", "").strip()
                if line.strip().isdigit() and not pmid: pmid = line.strip()
            if tit and pmid:
                news.append({"titulo": tit, "fonte": journal[:30], "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/", "img":"https://images.unsplash.com/photo-1532187863486-abf9dbad1b69?w=400"})
        return news
    except: return []
