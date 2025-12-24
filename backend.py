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
# Expande termos simples para queries complexas de Fisiologia e Farmacologia
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

# --- 2. FAXINEIRO IA (PERSONA: FARM ACOLOGISTA SENIOR) ---
def _faxina_ia(lista_suja):
    """
    Usa a IA com Persona de Pesquisador Sênior para filtrar estritamente alvos moleculares.
    """
    api_key = st.session_state.get('api_key_usuario', '')
    
    # Se não tiver chave, retorna a lista bruta cortada (Fallback seguro)
    if not api_key: return lista_suja[:25] 

    lista_str = ", ".join(lista_suja)
    
    # PROMPT DE "ALTA FISIOLOGIA"
    prompt = f"""
    ACT AS: Senior Research Pharmacologist & Physiologist.
    CONTEXT: We are mining PubMed for novel molecular targets in specific organs.
    
    INPUT LIST: {lista_str}
    
    TASK: Curate this list to keep ONLY basic science/molecular targets.
    
    STRICT FILTERING RULES:
    1. KEEP ONLY: 
       - Specific Receptors (e.g., TRPV1, P2X3, Beta-3-AR, GLP1R)
       - Ion Channels (e.g., BK, SK, Nav1.8)
       - Enzymes/Kinases (e.g., mTOR, PI3K, PDE5, COX2)
       - Signaling Molecules/Gasotransmitters (e.g., ATP, NO, H2S, PGE2)
       - Transcription Factors (e.g., NF-kB, NRF2, HIF-1a)
       - Specific Genes/RNAs (e.g., GATA3, miRNA-132)
       
    2. DELETE IMMEDIATELY (ZERO TOLERANCE for Clinical/Surgical terms):
       - Diseases/Phenotypes (e.g., LUTS, OAB, BPH, Cancer, Infection, COVID, SUI)
       - Diagnostic/Imaging (e.g., MRI, CT, VIRADS, PIRADS, Ultrasound)
       - Procedures/Surgeries (e.g., TURBT, SNM, Botox, Transplantation)
       - Study Types (e.g., RCT, Review, Cohort, Analysis)
       - General Anatomy (e.g., Bladder, Heart, Urothelium)
       - Acronyms of Associations (e.g., EAU, AUA, SUFU)
    
    OUTPUT FORMAT: Return ONLY a clean Python list of strings.
    Example: ['TRPV4', 'NGF', 'BDNF', 'P2X3', 'YAP1']
    """
    
    headers = {'Content-Type': 'application/json'}
    data = {
        "contents": [{"parts": [{"text": prompt}]}],
        "safetySettings": [{"category": "HARM_CATEGORY_DANGEROUS_CONTENT", "threshold": "BLOCK_NONE"}],
        "generationConfig": {"temperature": 0.0} # Temperatura Zero = Rigidez Máxima
    }
    
    # Tenta modelos rápidos primeiro, depois os mais robustos
    modelos = ["gemini-2.5-flash", "gemini-2.0-flash", "gemini-flash-latest"]
    
    for m in modelos:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=10) # Timeout maior para processar listas grandes
            
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                # Limpeza da resposta para garantir formato lista
                texto = texto.replace("```python", "").replace("```json", "").replace("```", "").strip()
                if texto.startswith("[") and texto.endswith("]"):
                    lista_limpa = ast.literal_eval(texto)
                    if isinstance(lista_limpa, list) and len(lista_limpa) > 0:
                        return lista_limpa
        except: continue
    
    # Se falhar, retorna a lista original (parcialmente limpa pelo Regex)
    return lista_suja[:30]

# --- 3. ANÁLISE DE RESUMOS ---
def analisar_abstract_com_ia(titulo, dados_curtos, api_key, lang='pt'):
    if not api_key: return "⚠️ IA não ativada"
    idioma = "Português" if lang == 'pt' else "Inglês"
    prompt_text = f"""Atue como Pesquisador em Farmacologia. Analise:
FONTE: {titulo}. {dados_curtos}
FORMATO: Alvo: [Sigla] | Fármaco: [Nome] | Efeito: [Ação funcional].
REGRAS: Máximo 15 palavras. Foco em mecanismo de ação. Idioma: {idioma}."""
    
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt_text}]}]}
    
    try:
        url = f"[https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash:generateContent?key=](https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash:generateContent?key=){api_key}"
        resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=8)
        if resp.status_code == 200:
            return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
    except: pass
    return "⚠️ IA indisponível."

# --- 4. BUSCA PUBMED (CORE) ---
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
                if line.strip().isdigit() and not pmid: pmid = line.strip()
                if line.startswith("TI  - "): tit = line[6:].strip()
                if line.startswith("AB  - "): abstract = line[6:500].strip()
                if line.startswith("OT  - ") or line.startswith("KW  - "): keywords += line[6:].strip() + ", "
            if tit:
                artigos.append({"Title": tit, "Info_IA": f"{keywords} {abstract}", "Link": f"[https://pubmed.ncbi.nlm.nih.gov/](https://pubmed.ncbi.nlm.nih.gov/){pmid}/"})
        return artigos
    except: return []

# --- 5. MOTOR DE MINERAÇÃO (ULTIMATE EDITION) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    
    termo_upper = termo_base.upper().strip()
    query_string = MAPA_SINONIMOS.get(termo_upper, f"{termo_base}[Title/Abstract]")
    
    # VOLUME MASSIVO: Buscamos até 1500 abstracts recentes.
    # OBS: "Milhões" em tempo real é impossível via API gratuita (timeout). 
    # 1500 é o "Sweet Spot" para uma demo rápida e impressionante.
    final_query = f"({query_string}) AND (2018:2030[Date - Publication]) AND (NOT Review[pt])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=final_query, retmax=1500, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        artigos_raw = full_data.split("\n\nPMID-")

        # --- BLACKLIST DE EXTERMÍNIO (CLÍNICA/CIRURGIA/DIAGNÓSTICO) ---
        blacklist = {
            # 1. Termos que você flagrou (LUTS, OAB, PROCEDIMENTOS)
            "LUTS", "SUI", "UUI", "MUI", "OAB", "BPH", "UTI", "IC", "BPS", "ICIRS", "LUT", "LUTSFORTA", "BOO",
            "SNM", "PTNS", "TURBT", "TURP", "BOTOX", "INJECTION", "STENT", "CATHETER", "SLING", "MESH", "GRAFT",
            "COVID", "COVID19", "SARS", "VIRUS", "INFECTION", "SEPSIS", "CANCER", "TUMOR", "CIS", "MIBC", "NMIBC",
            "VIRADS", "PIRADS", "CEUS", "MRI", "CT", "PET", "ULTRASOUND", "ACR", "AUC", "ROC", "CI", "OR", "HR", "P",
            "SNAP25", "SV2", "VVF", "UVF", "SUFU", "ALE", "FOR", "AND", "THE", "WITH", "FROM", "BUT", "NOT",
            "TRANSPLANTATION", "NEPHRECTOMY", "LITHOTRIPSY", "URS", "RIRS", "HOLEP", "THULEP", "GREENLIGHT",
            "RADIOTHERAPY", "RADIATION", "TOXICITY", "SAFETY", "EFFICACY", "CLINICAL", "PRECLINICAL", "TRIAL",
            "SACRAL", "NERVE", "STIMULATION", "MODULATION", "DHEA", "BAG3", "CK20", "PSA", "PROSTATE",
            "RRNA", "TRNA", "DNA", "RNA", "GENE", "PROTEIN", "CELL", "TISSUE", "PATIENT", "HUMAN",
            "PUBLISH", "ACCEPTED", "RECEIVED", "REVISED", "CORRESPONDENCE", "EMAIL", "AUTHOR", "EDITOR",
            "PII", "DOI", "ISSN", "PMID", "PMC", "ISBN", "COPYRIGHT", "VOLUME", "ISSUE", "PAGE",
            "AHA", "ACC", "ESC", "WHO", "NIH", "CDC", "FDA", "EMA", "ISO", "AUA", "EAU", "ICS", "IUGA",
            "UNIVERSITY", "DEPARTMENT", "HOSPITAL", "CENTER", "CENTRE", "SCHOOL", "COLLEGE", "INSTITUTE",
            "REVIEW", "META", "ANALYSIS", "STUDY", "DATA", "USING", "USED", "BETWEEN", "AMONG", "WITHIN",
            "TOTAL", "MEAN", "RATIO", "SD", "SEM", "YEAR", "MONTH", "DAY", "HOUR", "MIN", "SEC",
            "WESTERN", "BLOT", "PCR", "QPCR", "ELISA", "STAINING", "IMMUNO", "HISTOLOGY", "ASSAY"
        }

        candidatos_por_artigo = []
        for artigo in artigos_raw:
            # Estratégia Sniper: Só lê TÍTULO e KEYWORDS (onde está o ouro)
            texto_focado = ""
            for line in artigo.split("\n"):
                if line.startswith("TI  - "): texto_focado += line[6:].strip() + " "
                elif line.startswith("OT  - ") or line.startswith("KW  - "): texto_focado += line[6:].strip() + " "
            
            if not texto_focado: continue

            # Regex Aprimorado: Pega siglas (TRPV1) e moléculas (TMAO, BDNF)
            # Exige pelo menos uma letra maiúscula e tamanho >= 3
            encontrados = re.findall(r'\b(?:[A-Z]{2,}[A-Z0-9-]*|[a-z]{1,2}[A-Z][a-zA-Z0-9-]*)\b', texto_focado)
            
            for t in encontrados:
                t_clean = re.sub(r'[^a-zA-Z0-9]', '', t).upper()
                
                # Filtros Locais
                if len(t_clean) < 3: continue 
                if t_clean in blacklist: continue
                if t_clean.isdigit(): continue
                if t_clean == termo_upper.replace(" ", ""): continue # Remove o próprio termo de busca (ex: HEART)
                
                candidatos_por_artigo.append(t_clean)

        if not candidatos_por_artigo: return []
        
        # Estatística de Frequência
        contagem = Counter(candidatos_por_artigo)
        total_docs = max(1, len(artigos_raw))
        
        # Pega os 120 termos mais relevantes (mas que não apareçam em 90% dos abstracts, pra evitar vício)
        top_candidatos = [termo for termo,freq in contagem.most_common(120) if (freq/total_docs)<0.90]
        
        # Filtro de Segurança Final (Garante que nada da blacklist passou)
        lista_final = [t for t in top_candidatos if t not in blacklist]
        
        # --- DECISÃO FINAL ---
        if usar_ia and st.session_state.get('api_key_usuario'):
            # Envia para o "Faxineiro PhD" (Prompt Rigoroso)
            # Mand
                                  
