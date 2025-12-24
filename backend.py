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
MAPA_SINONIMOS = {
    "HEART": "Heart OR Cardiac OR Myocardium OR Cardiomyocyte OR Coronary OR Artery OR Ventricular OR Atrial OR Ischemia",
    "BLADDER": "Bladder OR Urothelium OR Detrusor OR Vesical OR Urethra OR Micturition OR LUTS OR Cystitis OR Overactive",
    "KIDNEY": "Kidney OR Renal OR Nephron OR Glomerulus OR Tubular OR Podocyte OR AKI OR CKD",
    "BRAIN": "Brain OR CNS OR Neuron OR Glia OR Cortex OR Hippocampus OR Synaptic OR Neurotransmitter OR Cognitive",
    "LIVER": "Liver OR Hepatic OR Hepatocyte OR Steatosis OR Fibrosis OR Cirrhosis OR NASH",
    "LUNG": "Lung OR Pulmonary OR Alveolar OR Bronchial OR Respiratory OR Asthma OR COPD",
    "INTESTINE": "Intestine OR Gut OR Colon OR Bowel OR Enteric OR Colitis OR Microbiota",
    "PAIN": "Pain OR Nociception OR Analgesia OR Neuropathic OR Hyperalgesia OR Dorsal Root Ganglion",
    "INFLAMMATION": "Inflammation OR Cytokine OR Macrophage OR Neutrophil OR Immune OR Sepsis OR Inflammasome",
    "METABOLISM": "Metabolism OR Obesity OR Diabetes OR Insulin OR Glucose OR Adipose OR Lipid",
    "CANCER": "Cancer OR Tumor OR Oncology OR Carcinoma OR Metastasis OR Proliferation OR Angiogenesis"
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
    
    # Lista de Modelos (Prioridade para os mais novos e rápidos)
    modelos = ["gemini-2.5-flash", "gemini-2.0-flash", "gemini-2.0-flash-exp", "gemini-flash-latest"]
    
    for m in modelos:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=8)
            if resp.status_code == 200:
                return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
        except: continue
    return f"❌ Falha técnica. Título: {titulo[:30]}..."

# --- BUSCA PUBMED COUNTS (Para Estatística Lambda) ---
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

# --- BUSCA DE ARTIGOS PARA LEITURA (Frontend) ---
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
                artigos.append({
                    "Title": tit, 
                    "Info_IA": f"{keywords} {abstract}", 
                    "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                })
        return artigos
    except: return []

# --- MOTOR DE MINERAÇÃO: "HEAVY DUTY" (V6.0) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email):
    if email: Entrez.email = email
    
    # 1. EXPANSÃO SEMÂNTICA
    termo_upper = termo_base.upper().strip()
    query_string = MAPA_SINONIMOS.get(termo_upper, f"{termo_base}[Title/Abstract]")
    
    # 2. BUSCA MASSIVA (1000 abstracts para garantir volume estatístico)
    # Aumentamos para 1000 porque vamos filtrar muito agressivamente
    final_query = f"({query_string}) AND (2018:2030[Date - Publication]) AND (NOT Review[pt])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=final_query, retmax=1000, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        artigos_raw = full_data.split("\n\nPMID-")

        # 3. BLACKLIST "NUCLEAR" (FILTRA TUDO QUE NÃO É ALVO MOLECULAR)
        blacklist = {
            # -- LIXO DE PUBLICAÇÃO --
            "PII", "DOI", "ISSN", "PMID", "PMC", "NIHMS", "ISBN", "COPYRIGHT", "AUTHOR", "AUTHORS", 
            "EDITOR", "PUBLISHER", "RECEIVED", "ACCEPTED", "PUBLISH", "PUBLISHED", "REVISED", "CORRESPONDENCE", 
            "EMAIL", "ARTICLE", "JOURNAL", "REPORT", "REPORTS", "VOLUME", "ISSUE", "PAGE", "PAGES", "SUPPL",
            "FIG", "FIGURE", "TABLE", "REF", "REFERENCE", "REFERENCES", "ABSTRACT", "TEXT", "FULL", "TITLE",
            "INTRODUCTION", "BACKGROUND", "METHODS", "RESULTS", "DISCUSSION", "CONCLUSION", "CONCLUSIONS",
            "UNIVERSITY", "DEPARTMENT", "HOSPITAL", "CENTER", "CENTRE", "SCHOOL", "COLLEGE", "INSTITUTE",
            "FOUNDATION", "SOCIETY", "ASSOCIATION", "FEDERATION", "COMMITTEE", "TASK", "FORCE", "GROUP",
            "SCIENCE", "NATURE", "CELL", "PLOS", "ONE", "FRONTIERS", "MDPI", "WILEY", "ELSEVIER",

            # -- SIGLAS ADMINISTRATIVAS & CONTEXTO (Cardio/Uro/Geral) --
            "AHA", "ACC", "ESC", "HFSA", "ACHD", "WHO", "NIH", "CDC", "FDA", "EMA", "ISO", "AUA", "EAU", "ICS",
            "GUIDELINE", "GUIDELINES", "STATEMENT", "POSITION", "CONSENSUS", "REGISTRY", "TRIAL", "TRIALS",
            "STUDY", "STUDIES", "COHORT", "PROSPECTIVE", "RETROSPECTIVE", "RANDOMIZED", "BLIND", "DOUBLE", "META",
            "SYSTEMATIC", "REVIEW", "CASE", "SERIES", "PILOT", "MULTICENTER", "ANALYSIS", "DATA", "USING",

            # -- TERMOS CLÍNICOS --
            "PATIENT", "PATIENTS", "HUMAN", "MEN", "WOMEN", "ADULT", "CHILD", "INFANT", "ELDERLY", "MALE", "FEMALE",
            "CONTROL", "PLACEBO", "SHAM", "TREATED", "UNTREATED", "VEHICLE", "BASELINE", "NORMAL", "ABNORMAL",
            "CLINICAL", "PRECLINICAL", "MEDICAL", "MEDICINE", "HEALTH", "CARE", "PUBLIC", "COMMUNITY", "GLOBAL",
            "DISEASE", "DISORDER", "SYNDROME", "CONDITION", "FAILURE", "INJURY", "DAMAGE", "STRESS", "TRAUMA",
            "ACUTE", "CHRONIC", "SEVERE", "MILD", "MODERATE", "STAGE", "GRADE", "CLASS", "RISK", "FACTOR", "SCORE",
            "DIAGNOSIS", "PROGNOSIS", "THERAPY", "TREATMENT", "MANAGEMENT", "STRATEGY", "OUTCOME", "OUTCOMES",
            "EVENT", "EVENTS", "DEATH", "MORTALITY", "MORBIDITY", "SURVIVAL", "SAFETY", "EFFICACY", "QUALITY", "LIFE",
            "SYSTOLIC", "DIASTOLIC", "CARDIAC", "HEART", "LEFT", "RIGHT", "VENTRICULAR", "ATRIAL", "MYOCARDIAL",
            "VASCULAR", "ARTERY", "VEIN", "AORTIC", "VALVE", "BLOOD", "PRESSURE", "FLOW", "RATE", "RHYTHM",
            "EJECTION", "FRACTION", "PRESERVED", "REDUCED", "HFPEF", "HFREF", "HFMREF", "COVID", "COVID19", "SARS",
            "UROLOGY", "NEPHROLOGY", "URINARY", "BLADDER", "RENAL", "KIDNEY", "LOWER", "TRACT", "SYMPTOMS", "LUTS",
            "PAIN", "SEPSIS", "CANCER", "TUMOR", "METASTASIS", "OBESITY", "DIABETES", "INSULIN", "GLUCOSE",
            
            # -- TERMOS GENÉRICOS DE BIOLOGIA (O QUE SOBRA É ALVO) --
            "DNA", "RNA", "MRNA", "MIRNA", "SIRNA", "GENE", "GENES", "PROTEIN", "PROTEINS", "CELL", "CELLS", 
            "TISSUE", "TISSUES", "RATS", "MICE", "RABBIT", "PIG", "DOG", "MODEL", "MODELS", "VIVO", "VITRO", 
            "EX", "SITU", "LEVEL", "LEVELS", "EXPRESSION", "REGULATION", "ACTIVITY", "ACTIVATION", "INHIBITION",
            "PATHWAY", "SIGNALING", "MECHANISM", "TARGET", "ROLE", "FUNCTION", "ACTION", "EFFECT", "EFFECTS",
            "MEDIATED", "INDUCED", "DEPENDENT", "INDEPENDENT", "ASSOCIATED", "RELATED", "LINKED", "POTENTIAL",
            "BIOMARKER", "CANDIDATE", "NOVEL", "NEW", "RECENT", "INSIGHTS", "UPDATE", "FUTURE", "PERSPECTIVE",
            "MOLECULAR", "GENETIC", "EPIGENETIC", "PHARMACOLOGY", "DRUG", "DRUGS", "RECEPTOR", "CHANNEL", "ENZYME",

            # -- EXAMES & UNIDADES --
            "MRI", "CMR", "CT", "PET", "SPECT", "ECG", "EKG", "ECHO", "ULTRASOUND", "XRAY", "IMAGING",
            "WESTERN", "BLOT", "PCR", "QPCR", "RT", "ELISA", "STAINING", "IMMUNO", "HISTOLOGY", "ASSAY",
            "TOTAL", "MEAN", "RATIO", "SD", "SEM", "CI", "HR", "RR", "OR", "VS", "VERSUS", "LOG",
            "YEAR", "YEARS", "MONTH", "MONTHS", "DAY", "DAYS", "HOUR", "HOURS", "MIN", "SEC", "TIME",
            "KG", "MG", "ML", "NM", "UM", "MM", "CM", "HZ", "MV", "MMHG",

            # -- INGLÊS COMUM --
            "AND", "THE", "FOR", "NOT", "BUT", "WITH", "FROM", "THIS", "THAT", "THESE", "THOSE",
            "WHICH", "WHAT", "WHEN", "WHERE", "WHO", "WHY", "HOW", "ANY", "ALL", "EACH", "EVERY",
            "HAVE", "HAS", "HAD", "WAS", "WERE", "BEEN", "BEING", "ARE", "IS", "CAN", "COULD",
            "SHOULD", "WOULD", "MAY", "MIGHT", "MUST", "WILL", "SHALL", "DOES", "DID", "DOING",
            "VIA", "DUE", "BETWEEN", "AMONG", "WITHIN", "WITHOUT", "UNDER", "ABOVE", "BELOW", 
            "AFTER", "BEFORE", "DURING", "SINCE", "UNTIL", "WHILE", "ONCE", "UPON", "INTO", "ONTO"
        }
        
        candidatos_por_artigo = []
        for artigo in artigos_raw:
            # 4. EXTRAÇÃO: SÓ TÍTULO (TI) E KEYWORDS (OT/KW)
            texto_focado = ""
            for line in artigo.split("\n"):
                if line.startswith("TI  - "): texto_focado += line[6:].strip() + " "
                elif line.startswith("OT  - ") or line.startswith("KW  - "): texto_focado += line[6:].strip() + " "
            
            if not texto_focado: continue

            # 5. REGEX FORTE (Case Sensitive)
            # Pega: "TRPV1", "mTOR", "NF-kB", "TGF-beta", "SGLT2", "CYP2J2"
            # Ignora: "Cardiac", "Heart", "Renal" (Só 1 maiúscula)
            encontrados = re.findall(r'\b(?:[A-Z]{2,}[A-Z0-9-]*|[a-z]*[A-Z][a-zA-Z0-9-]*[0-9]+[a-zA-Z0-9-]*)\b', texto_focado)
            
            for t in encontrados:
                t_clean = re.sub(r'[^A-Z0-9]', '', t).upper() # Normaliza para checar blacklist
                
                if len(t_clean) < 3: continue
                if t_clean in blacklist: continue
                if t_clean == termo_upper.replace(" ", ""): continue
                if t_clean.isdigit(): continue 
                
                # Adiciona o termo limpo
                candidatos_por_artigo.append(t_clean)

        if not candidatos_por_artigo: return []
        
        contagem = Counter(candidatos_por_artigo)
        total_docs = max(1, len(artigos_raw))
        
        # 6. RETORNO MASSIVO: Top 70 (Para garantir >60 na tela)
        # Filtro de corte: Aceita termos
