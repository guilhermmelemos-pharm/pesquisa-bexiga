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
MODELOS_ATIVOS = ["gemini-2.0-flash", "gemini-1.5-flash", "gemini-1.5-pro"]

# --- FUNÇÃO AUXILIAR PARA LIMPAR URL ---
def montar_url_limpa(modelo, chave):
    p1, p2 = "https://generativelanguage", ".googleapis.com/v1beta/models"
    url = f"{p1}{p2}/{modelo}:generateContent?key={chave}"
    return url.replace("[", "").replace("]", "").replace("(", "").replace(")", "").strip()

# --- 2. FAXINEIRO IA (PROMPT DOUTORADO UNIFESP) ---
def _faxina_ia(lista_suja):
    api_key = st.session_state.get('api_key_usuario', '').strip()
    if not api_key: return lista_suja[:60] 

    lista_str = ", ".join(lista_suja)
    
    prompt = f"""
    ROLE: Senior PhD Researcher in Pharmacology and Physiopathology.
    TASK: Curate a list of scientific terms for Organ Bath and Molecular Biology research.
    
    INPUT: {lista_str}
    
    STRICT FILTERING RULES:
    1. KEEP ONLY: Specific Pharmacological Targets (TRPV1, P2X3, Beta-3-AR), Enzymes (mTOR, ROCK, COX-2), Specific Small Molecules/Metabolites (Trehalose, TMAO, Resveratrol), and experimental drugs.
    2. DISCARD ALL: Clinical terms (Overactive, Incontinence, LUTS, Treatment), Anatomy (Detrusor, Urothelium, Bladder), and vague biological concepts (Pathophysiology, Stress, Chronic, Outcome).
    3. DISCARD: General procedures (Surgical, Robotic, Management, Diagnosis).
    4. FOCUS: If you cannot test it in an isolated organ bath or quantify it via Western Blot/qPCR, delete it.
    
    OUTPUT: Return strictly a Python list of strings. Include up to 60 relevant items if possible.
    """
    
    headers = {'Content-Type': 'application/json'}
    data = {
        "contents": [{"parts": [{"text": prompt}]}],
        "safetySettings": [{"category": "HARM_CATEGORY_DANGEROUS_CONTENT", "threshold": "BLOCK_NONE"}],
        "generationConfig": {"temperature": 0.1}
    }

    for m in MODELOS_ATIVOS:
        try:
            url = montar_url_limpa(m, api_key)
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=20)
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                texto = texto.replace("```python", "").replace("```json", "").replace("```", "").strip()
                if texto.startswith("[") and texto.endswith("]"):
                    lista_limpa = ast.literal_eval(texto)
                    return [str(x) for x in lista_limpa if isinstance(x, str)]
        except: continue
    return lista_suja[:60]

# --- 3. ANÁLISE DE RESUMOS (FOCO EM BANCADA) ---
def analisar_abstract_com_ia(titulo, dados_curtos, api_key, lang='pt'):
    key = api_key.strip()
    if not key: return "⚠️ Chave API não detectada."
    idioma = "Português" if lang == 'pt' else "Inglês"
    
    prompt_text = f"""
    Atue como Pesquisador PhD em Farmacologia. 
    Analise o texto e extraia Alvo Molecular, Fármaco e Mecanismo de Ação celular.
    PROIBIDO: Dar conselhos médicos ou descrever terapias comportamentais (como treinamento vesical).
    FOCO: Sinalização intracelular e contratilidade tecidual.
    FORMATO: "Alvo: [X] | Fármaco: [Y] | Mecanismo: [Z]".
    DADOS: {titulo}. {dados_curtos}.
    RESPOSTA EM: {idioma}. Máximo 20 palavras.
    """
    
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt_text}]}]}

    for m in MODELOS_ATIVOS:
        try:
            url = montar_url_limpa(m, key)
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=10)
            if resp.status_code == 200:
                return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
        except: continue
    return "⚠️ IA indisponível no momento."

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
    query = f"({termo})" + (f" AND ({q_contexto})" if q_contexto else "")
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
            tit, keywords, abstract, pmid = "", "", "", ""
            for line in raw.split("\n"):
                if line.startswith("TI  - "): tit = line[6:].strip()
                if line.startswith("AB  - "): abstract = line[6:500].strip()
                if line.startswith("OT  - ") or line.startswith("KW  - "): keywords += line[6:].strip() + ", "
                if line.startswith("PMID- "): pmid = line[6:].strip()
            if tit:
                artigos.append({"Title": tit, "Info_IA": f"{keywords} {abstract}", "Link": f"[https://pubmed.ncbi.nlm.nih.gov/](https://pubmed.ncbi.nlm.nih.gov/){pmid}/"})
        return artigos
    except: return []

# --- 5. MINERAÇÃO ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    termo_upper = termo_base.upper().strip()
    query_string = MAPA_SINONIMOS.get(termo_upper, f"{termo_base}[Title/Abstract]")
    final_query = f"({query_string}) AND (2018:2030[Date - Publication]) AND (NOT Review[pt])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=final_query, retmax=2000, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        artigos_raw = handle.read().split("\n\nPMID-"); handle.close()

        blacklist_exterminio = {
            "LUTS", "OAB", "BPH", "UTI", "IC", "BPS", "DETRUSOR", "UROTHELIUM", "NERVE", "IMMUNE", 
            "CURRENT", "ROLE", "WHAT", "CASE", "CLINICAL", "EFFICACY", "TREATMENT", "MANAGEMENT", 
            "STUDY", "ANALYSIS", "REVIEW", "DATA", "RESULTS", "DNA", "RNA", "FACTOR", "LEVEL"
        }

        candidatos = []
        for artigo in artigos_raw:
            texto = ""
            for line in artigo.split("\n"):
                if line.startswith("TI  - ") or line.startswith("OT  - "): texto += line[6:].strip() + " "
            encontrados = re.findall(r'\b(?:[A-Z]{2,}[A-Z0-9-]*|[A-Z][a-z]{3,}[a-z0-9-]*)\b', texto)
            for t in encontrados:
                t_clean = re.sub(r'[^a-zA-Z0-9]', '', t)
                if len(t_clean) < 4 and t_clean.upper() not in ["NO", "CO", "H2S", "ATP", "BK", "M3"]: continue 
                if t_clean.upper() in blacklist_exterminio: continue
                candidatos.append(t_clean)

        if not candidatos: return []
        top_candidatos = [t for t, f in Counter(candidatos).most_common(200)]
        
        if usar_ia and st.session_state.get('api_key_usuario'):
            return _faxina_ia(top_candidatos)
        return top_candidatos[:60]
    except: return []

def buscar_todas_noticias(lang='pt'):
    try:
        query = "(molecular biology OR pharmacology) AND (2024/09/01:2025/12/31[Date - Publication])"
        handle = Entrez.esearch(db="pubmed", term=query, retmax=6, sort="pub_date")
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read().split("\n\nPMID-"); handle.close()
        news = []
        for art in dados:
            tit, pmid, journal = "", "", ""
            for line in art.split("\n"):
                if line.startswith("TI  - "): tit = line[6:].strip()
                if line.startswith("JT  - "): journal = line[6:36].strip()
                if "PMID-" in line: pmid = line.split("-")[-1].strip()
            if tit: news.append({"titulo": tit, "fonte": journal, "link": f"[https://pubmed.ncbi.nlm.nih.gov/](https://pubmed.ncbi.nlm.nih.gov/){pmid}/"})
        return news
    except: return []
