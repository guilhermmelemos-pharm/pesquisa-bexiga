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

# --- LISTA DE MODELOS ---
MODELOS_ATIVOS = [
    "gemini-2.0-flash",          
    "gemini-2.0-flash-exp",      
    "gemini-1.5-flash",          
    "gemini-1.5-pro",
    "gemini-flash-latest"
]

# --- FUNÇÃO AUXILIAR PARA LIMPAR URL ---
def montar_url_limpa(modelo, chave):
    parte1 = "https://generativelanguage"
    parte2 = ".googleapis.com/v1beta/models"
    base = parte1 + parte2
    url = f"{base}/{modelo}:generateContent?key={chave}"
    url = url.replace("[", "").replace("]", "").replace("(", "").replace(")", "").strip()
    return url

# --- 2. DOUTORA PROATIVA (CURADORIA + DESCOBERTA DE ALVOS) ---
def _faxina_ia(lista_suja):
    termo_busca = st.session_state.get('input_alvo', 'Baixo Trato Urinário')
    api_key = st.session_state.get('api_key_usuario', '').strip()
    if not api_key: return lista_suja[:30] 

    lista_str = ", ".join(lista_suja)
    
    # PROMPT DE ELITE: BRAINSTORMING + FILTRO DE BANCADA
    prompt = f"""
    ROLE: Senior Scientist & PhD in Pharmacology (Lab Focus: Organ Bath & Molecular Biology).
    
    CONTEXT: Research on "{termo_busca}".
    
    TASK (PROACTIVE DOCTOR): 
    1. BRAINSTORM: Use your internal PhD knowledge to identify under-exploited or novel molecular targets, receptors, and drugs for "{termo_busca}" not yet saturated in the literature.
    2. CROSS-CHECK: Combine your knowledge with this PubMed list: {lista_str}.
    
    ✅ USE THESE EXAMPLES AS QUALITY STANDARDS:
    - Targets: TRPV4, PIEZO1, P2X3, TRPM8, TRPV1, Beta-3-AR, ROCK, mTOR, SPHK1.
    - Molecules: ATP, NO, PGE2, TMAO, Aldosterone.
    - Drugs: GSK1016790A, Ouabain, O-1602, Mirabegron, Trehalose, Resveratrol.
    - Signaling/miRNAs: miR-132, SNARE, NLRP3.
    
    ❌ STRICT DELETION RULES:
    - No medical talk: Incontinence, BPH, Diabetes, Treatment, Patient, Symptoms, Clinicial.
    - No general biology: Microbiome, Stem, Cells, Organoids, Toxin, Immune, Bacillus.
    - No academic noise: Role, Effect, Study, Comparison, Mechanism, Action, Function.
    
    OUTPUT: Return strictly a Python list of strings (60 to 100 items). 
    Focus on "Hand-on" lab targets. No medical talk, vacilão.
    """
    
    headers = {'Content-Type': 'application/json'}
    data = {
        "contents": [{"parts": [{"text": prompt}]}],
        "generationConfig": {"temperature": 0.4} # Maior temperatura para proatividade
    }

    for m in MODELOS_ATIVOS:
        try:
            url_final = montar_url_limpa(m, api_key)
            resp = requests.post(url_final, headers=headers, data=json.dumps(data), timeout=25)
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                texto = texto.replace("```python", "").replace("```json", "").replace("```", "").strip()
                lista_limpa = ast.literal_eval(texto)
                if isinstance(lista_limpa, list) and len(lista_limpa) > 0:
                    return lista_limpa
        except: continue
    return lista_suja[:40]

# --- 3. ANÁLISE DE RESUMOS ---
def analisar_abstract_com_ia(titulo, dados_curtos, api_key, lang='pt'):
    key = api_key.strip()
    if not key: return "⚠️ Chave API não detectada."
    idioma = "Português" if lang == 'pt' else "Inglês"
    prompt_text = f"Doutora, resuma Alvo/Fármaco e Mecanismo funcional (lab focus) em 15 palavras. {titulo}. {dados_curtos}. Idioma: {idioma}."
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt_text}]}]}
    for m in MODELOS_ATIVOS:
        try:
            url_final = montar_url_limpa(m, key)
            resp = requests.post(url_final, headers=headers, data=json.dumps(data), timeout=10)
            if resp.status_code == 200:
                return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
        except: continue
    return "⚠️ Erro na análise da IA."

# --- 4. FUNÇÕES DE BUSCA ---
@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
def _fetch_pubmed_count(query):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    record = Entrez.read(handle); handle.close()
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
                artigos.append({
                    "Title": tit, 
                    "Info_IA": f"{keywords} {abstract}", 
                    "Link": f"[https://pubmed.ncbi.nlm.nih.gov/](https://pubmed.ncbi.nlm.nih.gov/){pmid}/"
                })
        return artigos
    except: return []

# --- 5. MINERAÇÃO (CORE) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    termo_upper = termo_base.upper().strip()
    query_string = MAPA_SINONIMOS.get(termo_upper, f"{termo_base}[Title/Abstract]")
    
    # 2500 Abstracts para dar volume à Doutora Proativa
    final_query = f"({query_string}) AND (2018:2030[Date - Publication]) AND (NOT Review[pt])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=final_query, retmax=2500, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        artigos_raw = full_data.split("\n\nPMID-")

        blacklist_exterminio = {
            "MICROBIOME", "MICROBIOTA", "BACILLUS", "CALMETTE", "GUERIN", "BCG", "UPEC",
            "STEM", "CELLS", "CELL", "ORGANOIDS", "TOXIN", "BOTULINUM", "ONABOTULINUM", "BTX",
            "CHECKPOINT", "IMMUNE", "ANTIBODIES", "VACCINE", "PLACEBO", "CLINICAL", "TRIAL", "PATIENT",
            "DETRUSOR", "UROTHELIUM", "NERVE", "MUCOSA", "SMOOTH", "MUSCLE", "GANGLION",
            "EPITHELIUM", "VESICAL", "URETHRA", "URINARY", "TRACT", "PELVIC", "LOWER",
            "CURRENT", "ROLE", "POTENTIAL", "EFFECT", "STUDY", "ANALYSIS", "REVIEW", "DATA",
            "TREATMENT", "MANAGEMENT", "DIAGNOSIS", "SYNDROME", "BPH", "OAB", "LUTS", "TURBT",
            "DNA", "RNA", "ROS", "VEGF", "CALCIUM", "PROTEIN", "GENE"
        }

        candidatos_por_artigo = []
        for artigo in artigos_raw:
            texto_focado = ""
            for line in artigo.split("\n"):
                if line.startswith("TI  - "): texto_focado += line[6:].strip() + " "
                elif line.startswith("OT  - ") or line.startswith("KW  - "): texto_focado += line[6:].strip() + " "
            
            # Regex para Siglas Técnicas (ex: TRPV4, miR-132, GSK-1016790A)
            encontrados = re.findall(r'\b(?:[A-Z0-9-]{2,}|[a-z]{1,2}-[A-Z0-9]{2,}|[a-z]{1,2}[A-Z][a-zA-Z0-9-]*)\b', texto_focado)
            
            for t in encontrados:
                t_clean = re.sub(r'[^a-zA-Z0-9-]', '', t)
                t_upper = t_clean.upper()
                if len(t_clean) < 3 or t_upper in blacklist_exterminio: continue
                candidatos_por_artigo.append(t_clean)

        if not candidatos_por_artigo: return []
        
        contagem = Counter(candidatos_por_artigo)
        top_candidatos = [termo for termo,freq in contagem.most_common(250)]
        
        if usar_ia and st.session_state.get('api_key_usuario'):
            return _faxina_ia(top_candidatos)
        else:
            return top_candidatos[:40]
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
                news.append({
                    "titulo": tit, "fonte": journal[:30], 
                    "link": f"[https://pubmed.ncbi.nlm.nih.gov/](https://pubmed.ncbi.nlm.nih.gov/){pmid}/"
                })
        return news
    except: return []
