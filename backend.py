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

# --- 1. MAPA DE EXPANSÃO SEMÂNTICA ---
MAPA_SINONIMOS = {
    "HEART": "Heart OR Cardiac OR Myocardium OR Cardiomyocyte OR Coronary OR Artery OR Ventricular OR Atrial OR Ischemia OR Heart Failure OR Cardiotoxicity",
    "BLADDER": "Bladder OR Urothelium OR Detrusor OR Vesical OR Urethra OR Micturition OR LUTS OR Cystitis OR Overactive Bladder OR OAB OR Urinary Tract",
    "KIDNEY": "Kidney OR Renal OR Nephron OR Glomerulus OR Tubular OR Podocyte OR AKI OR CKD OR Renal Fibrosis",
    "BRAIN": "Brain OR CNS OR Neuron OR Glia OR Cortex OR Hippocampus OR Synaptic OR Neurotransmitter OR Cognitive OR Neurodegeneration",
    "LIVER": "Liver OR Hepatic OR Hepatocyte OR Steatosis OR Fibrosis OR Cirrhosis OR NASH OR NAFLD OR Cytochrome P450",
    "LUNG": "Lung OR Pulmonary OR Alveolar OR Bronchial OR Respiratory OR Asthma OR COPD OR Pulmonary Fibrosis",
    "INTESTINE": "Intestine OR Gut OR Colon OR Bowel OR Enteric OR Colitis OR Microbiota OR IBD OR Epithelial Barrier",
    "PAIN": "Pain OR Nociception OR Analgesia OR Neuropathic OR Hyperalgesia OR Dorsal Root Ganglion OR Allodynia OR TRP Channels",
    "INFLAMMATION": "Inflammation OR Cytokine OR Macrophage OR Neutrophil OR Immune OR Sepsis OR Inflammasome OR T-cell OR NF-kB",
    "METABOLISM": "Metabolism OR Obesity OR Diabetes OR Insulin OR Glucose OR Adipose OR Lipid OR Metabolic Syndrome OR Mitochondria",
    "CANCER": "Cancer OR Tumor OR Oncology OR Carcinoma OR Metastasis OR Proliferation OR Angiogenesis OR Apoptosis OR Microenvironment"
}

# --- LISTA DE MODELOS (Seus modelos potentes) ---
MODELOS_ATIVOS = [
    "gemini-2.5-flash",          
    "gemini-2.0-flash",          
    "gemini-2.0-flash-exp",      
    "gemini-flash-latest",       
    "gemini-1.5-pro"
]

# --- 2. FAXINEIRO IA (O CERÉBRO DA OPERAÇÃO) ---
def _faxina_ia(lista_suja):
    """
    Envia uma lista bruta gigante para a IA classificar semanticamente.
    """
    api_key = st.session_state.get('api_key_usuario', '')
    
    # Se não tiver chave, infelizmente temos que retornar a lista suja (cortada)
    if not api_key: return lista_suja[:40] 

    # Transformamos a lista em string para o prompt
    lista_str = ", ".join(lista_suja)
    
    prompt = f"""
    ROLE: Expert Senior Pharmacologist & Data Scientist.
    
    INPUT: The following list of terms extracted from PubMed titles:
    [{lista_str}]
    
    TASK: Filter this list rigorously. We are looking ONLY for Molecular Targets and Mechanisms.
    
    CRITERIA FOR KEEPING (INCLUDE):
    - Specific Genes (e.g., mTOR, TFEB, GATA3)
    - Receptors & Channels (e.g., TRPV1, P2X3, Beta-3-AR)
    - Signaling Molecules & Enzymes (e.g., ATP, NO, PGE2, COX2)
    - Specific RNAs (e.g., miRNA-132)
    - Drugs/Compounds (e.g., Mirabegron, Resiniferatoxin)
    
    CRITERIA FOR DELETING (EXCLUDE):
    - Locations/Geography (e.g., China, USA, London)
    - Clinical Conditions/Diseases (e.g., LUTS, OAB, Cancer, Infection, COVID)
    - Medical Procedures (e.g., TURBT, MRI, Surgery, Injection)
    - Study Types (e.g., RCT, Review, Meta-analysis)
    - Organizations (e.g., AHA, WHO)
    - General Biological Terms (e.g., Protein, Gene, Cell, Study, Data)
    
    OUTPUT FORMAT: 
    Return strictly a Python List of strings containing ONLY the valid targets.
    Example: ['TRPV1', 'NGF', 'ATP']
    """
    
    headers = {'Content-Type': 'application/json'}
    data = {
        "contents": [{"parts": [{"text": prompt}]}],
        "safetySettings": [{"category": "HARM_CATEGORY_DANGEROUS_CONTENT", "threshold": "BLOCK_NONE"}],
        "generationConfig": {"temperature": 0.0}
    }
    
    base_url = "https://generativelanguage.googleapis.com/v1beta/models"

    for m in MODELOS_ATIVOS:
        try:
            url = f"{base_url}/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=15) # Mais tempo para processar lista grande
            
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                # Limpeza cirúrgica da resposta
                texto = texto.replace("```python", "").replace("```json", "").replace("```", "").strip()
                
                # Tenta parsear a lista
                if texto.startswith("[") and texto.endswith("]"):
                    lista_limpa = ast.literal_eval(texto)
                    if isinstance(lista_limpa, list) and len(lista_limpa) > 0:
                        # Sucesso! Retorna a lista curada pela IA
                        return lista_limpa
        except: continue
    
    # Fallback se a IA falhar
    return lista_suja[:40]

# --- 3. ANÁLISE DE RESUMOS ---
def analisar_abstract_com_ia(titulo, dados_curtos, api_key, lang='pt'):
    if not api_key: return "⚠️ Chave API não detectada."
    
    idioma = "Português" if lang == 'pt' else "Inglês"
    
    prompt_text = f"""Atue como Pesquisador em Farmacologia. Analise:
    TITULO: {titulo}
    CONTEXTO: {dados_curtos}
    TAREFA: Resuma o Alvo Molecular, o Fármaco e o Efeito. 
    REGRAS: Máximo 20 palavras. Idioma: {idioma}."""
    
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt_text}]}]}
    
    base_url = "[https://generativelanguage.googleapis.com/v1beta/models](https://generativelanguage.googleapis.com/v1beta/models)"
    
    for m in MODELOS_ATIVOS:
        try:
            url = f"{base_url}/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=10)
            if resp.status_code == 200:
                return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
        except: continue
            
    return "⚠️ IA indisponível."

# --- 4. FUNÇÕES DE BUSCA ---
@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
def _fetch_pubmed_count(query):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    record = Entrez.read(handle)
    handle.close()
    return int(record["Count"])

@st.cache_data(ttl=86400, show_spinner=False)
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    q_contexto = MAPA_SINONIMOS.get(contexto.upper(), contexto) if contexto else ""
    query = f"({termo})"
    if q_contexto: query += f" AND ({q_contexto})"
    query += f" AND ({ano_ini}:{ano_fim}[Date - Publication]) AND (NOT Review[pt])"
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
                artigos.append({"Title": tit, "Info_IA": f"{keywords} {abstract}", "Link": f"[https://pubmed.ncbi.nlm.nih.gov/](https://pubmed.ncbi.nlm.nih.gov/){pmid}/"})
        return artigos
    except: return []

# --- 5. MINERAÇÃO IA-CENTRIC (A Nova Lógica) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    
    termo_upper = termo_base.upper().strip()
    query_string = MAPA_SINONIMOS.get(termo_upper, f"{termo_base}[Title/Abstract]")
    
    # 1. Busca Massiva (1500 abstracts)
    final_query = f"({query_string}) AND (2018:2030[Date - Publication]) AND (NOT Review[pt])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=final_query, retmax=1500, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        artigos_raw = full_data.split("\n\nPMID-")

        # --- BLACKLIST MÍNIMA (SÓ LIXO GRAMATICAL) ---
        # Não filtramos mais doenças ou procedimentos aqui. Deixamos a IA decidir.
        # Filtramos apenas o que é OBVIAMENTE lixo linguístico para não gastar token.
        blacklist_tecnica = {
            "AND", "THE", "FOR", "NOT", "BUT", "WITH", "FROM", "THIS", "THAT", "THESE", "THOSE",
            "WHICH", "WHAT", "WHEN", "WHERE", "WHO", "WHY", "HOW", "ANY", "ALL", "EACH", "EVERY",
            "HAVE", "HAS", "HAD", "WAS", "WERE", "BEEN", "BEING", "ARE", "IS", "CAN", "COULD",
            "SHOULD", "WOULD", "MAY", "MIGHT", "MUST", "WILL", "SHALL", "DOES", "DID", "DOING",
            "VIA", "DUE", "BETWEEN", "AMONG", "WITHIN", "WITHOUT", "UNDER", "ABOVE", "BELOW", 
            "AFTER", "BEFORE", "DURING", "SINCE", "UNTIL", "WHILE", "ONCE", "UPON", "INTO", "ONTO",
            "TOTAL", "MEAN", "RATIO", "SD", "SEM", "YEAR", "MONTH", "DAY", "HOUR", "MIN", "SEC",
            "USING", "USED", "USE", "DATA", "ANALYSIS", "STUDY", "RESULTS", "CONCLUSION", "BACKGROUND",
            "METHODS", "OBJECTIVE", "AIM", "PURPOSE", "DEPARTMENT", "UNIVERSITY", "HOSPITAL",
            "PUBLISH", "ACCEPTED", "RECEIVED", "REVISED", "CORRESPONDENCE", "EMAIL", "AUTHOR", "EDITOR",
            "PII", "DOI", "ISSN", "PMID", "PMC", "ISBN", "COPYRIGHT", "VOLUME", "ISSUE", "PAGE",
            "FIG", "FIGURE", "TABLE", "SUPPL", "TEXT", "FULL", "ABSTRACT", "TITLE"
        }

        candidatos_por_artigo = []
        for artigo in artigos_raw:
            texto_focado = ""
            for line in artigo.split("\n"):
                # Foco em Título e Keywords
                if line.startswith("TI  - "): texto_focado += line[6:].strip() + " "
                elif line.startswith("OT  - ") or line.startswith("KW  - "): texto_focado += line[6:].strip() + " "
            
            if not texto_focado: continue

            # Regex Permissivo: Pega qualquer coisa que pareça sigla ou termo técnico
            encontrados = re.findall(r'\b(?:[A-Z]{2,}[A-Z0-9-]*|[a-z]{1,2}[A-Z][a-zA-Z0-9-]*)\b', texto_focado)
            
            for t in encontrados:
                t_clean = re.sub(r'[^a-zA-Z0-9]', '', t).upper()
                if len(t_clean) < 3: continue 
                if t_clean in blacklist_tecnica: continue # Só remove preposições
                if t_clean.isdigit(): continue
                if t_clean == termo_upper.replace(" ", ""): continue
                
                candidatos_por_artigo.append(t_clean)

        if not candidatos_por_artigo: return []
        
        # Estatística
        contagem = Counter(candidatos_por_artigo)
        total_docs = max(1, len(artigos_raw))
        
        # PEGAMOS UMA AMOSTRA GRANDE (Top 150)
        # Enviamos "lixo" junto (LUTS, OAB, TURBT) propositalmente para a IA filtrar
        top_candidatos = [termo for termo,freq in contagem.most_common(150) if (freq/total_docs)<0.90]
        
        # --- A MÁGICA ACONTECE AQUI ---
        if usar_ia and st.session_state.get('api_key_usuario'):
            # Envia o "Pacotão" de 150 termos para a IA limpar
            return _faxina_ia(top_candidatos)
        else:
            # Se não tiver IA, temos que ser conservadores e cortar em 30
            return top_candidatos[:30]

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
                if line.startswith("TI  - "): tit = line[6:].strip()
                if line.startswith("JT  - "): journal = line.replace("JT  - ", "").strip()
                if line.strip().isdigit() and not pmid: pmid = line.strip()
            if tit and pmid:
                news.append({"titulo": tit, "fonte": journal[:30], "link": f"[https://pubmed.ncbi.nlm.nih.gov/](https://pubmed.ncbi.nlm.nih.gov/){pmid}/", "img":"[https://images.unsplash.com/photo-1532187863486-abf9dbad1b69?w=400](https://images.unsplash.com/photo-1532187863486-abf9dbad1b69?w=400)"})
        return news
    except: return []
