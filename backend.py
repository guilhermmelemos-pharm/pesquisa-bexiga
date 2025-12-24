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
# O "Cérebro" que traduz termos simples para buscas complexas
MAPA_SINONIMOS = {
    "HEART": "Heart OR Cardiac OR Myocardium OR Cardiomyocyte OR Coronary OR Artery OR Ventricular OR Atrial OR Ischemia OR Heart Failure",
    "BLADDER": "Bladder OR Urothelium OR Detrusor OR Vesical OR Urethra OR Micturition OR LUTS OR Cystitis OR Overactive Bladder OR OAB",
    "KIDNEY": "Kidney OR Renal OR Nephron OR Glomerulus OR Tubular OR Podocyte OR AKI OR CKD OR Fibrosis",
    "BRAIN": "Brain OR CNS OR Neuron OR Glia OR Cortex OR Hippocampus OR Synaptic OR Neurotransmitter OR Cognitive OR Alzheimer",
    "LIVER": "Liver OR Hepatic OR Hepatocyte OR Steatosis OR Fibrosis OR Cirrhosis OR NASH OR NAFLD",
    "LUNG": "Lung OR Pulmonary OR Alveolar OR Bronchial OR Respiratory OR Asthma OR COPD OR Fibrosis",
    "INTESTINE": "Intestine OR Gut OR Colon OR Bowel OR Enteric OR Colitis OR Microbiota OR IBD",
    "PAIN": "Pain OR Nociception OR Analgesia OR Neuropathic OR Hyperalgesia OR Dorsal Root Ganglion OR Allodynia",
    "INFLAMMATION": "Inflammation OR Cytokine OR Macrophage OR Neutrophil OR Immune OR Sepsis OR Inflammasome OR T-cell",
    "METABOLISM": "Metabolism OR Obesity OR Diabetes OR Insulin OR Glucose OR Adipose OR Lipid OR Metabolic Syndrome",
    "CANCER": "Cancer OR Tumor OR Oncology OR Carcinoma OR Metastasis OR Proliferation OR Angiogenesis OR Apoptosis"
}

# --- 2. FAXINEIRO IA (LIMPEZA INTELIGENTE) ---
def _faxina_ia(lista_suja):
    """
    Usa a IA para separar Alvos Moleculares Reais de Lixo Clínico/Administrativo.
    """
    api_key = st.session_state.get('api_key_usuario', '')
    
    # Se não tiver chave, retorna a lista suja cortada (Backup)
    if not api_key: 
        return lista_suja[:25] 

    lista_str = ", ".join(lista_suja)
    
    prompt = f"""
    Task: Pharmacological Data Cleaning.
    Input: {lista_str}
    
    Instruction: Filter this list extracted from PubMed titles. 
    REMOVE: Clinical terms (OAB, BPH, UTI), Procedures (TURBT, MRI), Organizations (AHA, EAU), Statistics (ANOVA), and general words.
    KEEP ONLY: Genes (e.g., TFEB, mTOR), Proteins, Receptors (e.g., TRPV1), Enzymes, Metabolites (e.g., TMAO), and Pharmacological Targets.
    
    Output format: A clean Python list string, e.g., ['TMAO', 'TFEB', 'TRPV1']
    """
    
    headers = {'Content-Type': 'application/json'}
    data = {
        "contents": [{"parts": [{"text": prompt}]}],
        "safetySettings": [{"category": "HARM_CATEGORY_DANGEROUS_CONTENT", "threshold": "BLOCK_NONE"}],
        "generationConfig": {"temperature": 0.0}
    }
    
    # Tenta modelos rápidos primeiro
    modelos = ["gemini-2.5-flash", "gemini-2.0-flash", "gemini-flash-latest"]
    
    for m in modelos:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=6)
            
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                # Limpeza da string JSON/Markdown
                texto = texto.replace("```python", "").replace("```json", "").replace("```", "").strip()
                
                if texto.startswith("[") and texto.endswith("]"):
                    lista_limpa = ast.literal_eval(texto)
                    if isinstance(lista_limpa, list) and len(lista_limpa) > 0:
                        return lista_limpa
        except: continue
            
    # Se a IA falhar em todos os modelos, retorna a lista original
    return lista_suja[:25]

# --- 3. ANÁLISE DE RESUMOS COM IA ---
def analisar_abstract_com_ia(titulo, dados_curtos, api_key, lang='pt'):
    if not api_key: return "⚠️ IA não ativada"
    idioma = "Português" if lang == 'pt' else "Inglês"
    prompt_text = f"""PhD em Farmacologia, analise:
FONTE: {titulo}. {dados_curtos}
FORMATO: Alvo: [Sigla] | Fármaco: [Nome] | Efeito: [Ação funcional].
REGRAS: Máximo 12 palavras. Seja técnico. Idioma: {idioma}."""
    
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt_text}]}]}
    
    try:
        url = f"[https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash:generateContent?key=](https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash:generateContent?key=){api_key}"
        resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=8)
        if resp.status_code == 200:
            return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
    except: pass
    return "⚠️ IA indisponível no momento."

# --- 4. FUNÇÕES DE BUSCA PUBMED (CORE) ---

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

# --- 5. MOTOR DE MINERAÇÃO HÍBRIDO (HEAVY DUTY V8) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email):
    if email: Entrez.email = email
    
    termo_upper = termo_base.upper().strip()
    query_string = MAPA_SINONIMOS.get(termo_upper, f"{termo_base}[Title/Abstract]")
    
    # BUSCA MASSIVA: 1000 abstracts para garantir diversidade
    final_query = f"({query_string}) AND (2018:2030[Date - Publication]) AND (NOT Review[pt])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=final_query, retmax=1000, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        artigos_raw = full_data.split("\n\nPMID-")

        # --- BLACKLIST SUPREMA (Proteção contra Lixo de Todas as Áreas) ---
        blacklist = set()
        
        # A. Metadados, Editoras & Publicação
        blacklist.update([
            "PII", "DOI", "ISSN", "PMID", "PMC", "NIHMS", "ISBN", "COPYRIGHT", "AUTHOR", "AUTHORS", 
            "EDITOR", "PUBLISHER", "RECEIVED", "ACCEPTED", "PUBLISH", "PUBLISHED", "REVISED", "CORRESPONDENCE", 
            "EMAIL", "ARTICLE", "JOURNAL", "REPORT", "REPORTS", "VOLUME", "ISSUE", "PAGE", "PAGES", "SUPPL",
            "FIG", "FIGURE", "TABLE", "REF", "REFERENCE", "REFERENCES", "ABSTRACT", "TEXT", "FULL", "TITLE",
            "INTRODUCTION", "BACKGROUND", "METHODS", "RESULTS", "DISCUSSION", "CONCLUSION", "CONCLUSIONS",
            "UNIVERSITY", "DEPARTMENT", "HOSPITAL", "CENTER", "CENTRE", "SCHOOL", "COLLEGE", "INSTITUTE",
            "FOUNDATION", "SOCIETY", "ASSOCIATION", "FEDERATION", "COMMITTEE", "TASK", "FORCE", "GROUP",
            "SCIENCE", "NATURE", "CELL", "PLOS", "ONE", "FRONTIERS", "MDPI", "WILEY", "ELSEVIER"
        ])

        # B. Urologia & Nefrologia (Procedimentos & Doenças)
        blacklist.update([
            "BPH", "LUTS", "OAB", "IC", "BPS", "ICIRS", "UTI", "UPEC", "SCI", "NB", "UI", "SUI", "MUI",
            "MIBC", "NMIBC", "CIS", "TURBT", "TURP", "BCG", "RCC", "TCC", "UC", "PCAR", "PSA", "CRPC",
            "URODYNAMICS", "CYSTOSCOPY", "VIRADS", "PIRADS", "PTNS", "SNM", "BOTOX", "ONABOTULINUMTOXINA",
            "HOLEP", "THULEP", "GREENLIGHT", "RES", "RESECTION", "STENT", "CATHETER", "SLING", "MESH",
            "TRANSPLANTATION", "GRAFT", "DIVERSION", "CONDUIT", "NEPHRECTOMY", "LITHOTRIPSY"
        ])

        # C. Cardiologia & Metabólico
        blacklist.update([
            "AHA", "ACC", "ESC", "HFSA", "ACHD", "STEMI", "NSTEMI", "HFPEF", "HFREF", "HFMREF", "CABG", "PCI",
            "LVAD", "HVAD", "ECMO", "TAVI", "TAVR", "PACEMAKER", "ICD", "CRT", "AFIB", "VT", "VF", "SCD",
            "SYSTOLIC", "DIASTOLIC", "EJECTION", "FRACTION", "PRESSURE", "HYPERTENSION", "CAD", "ACS",
            "BMI", "OBESITY", "DIABETES", "INSULIN", "GLUCOSE", "HBA1C", "LIPID", "CHOLESTEROL"
        ])

        # D. Termos Clínicos Gerais (Todas as áreas)
        blacklist.update([
            "PATIENT", "PATIENTS", "HUMAN", "MEN", "WOMEN", "ADULT", "CHILD", "INFANT", "ELDERLY", "MALE", "FEMALE",
            "CONTROL", "PLACEBO", "SHAM", "TREATED", "UNTREATED", "VEHICLE", "BASELINE", "NORMAL", "ABNORMAL",
            "CLINICAL", "PRECLINICAL", "MEDICAL", "MEDICINE", "HEALTH", "CARE", "PUBLIC", "COMMUNITY", "GLOBAL",
            "DISEASE", "DISORDER", "SYNDROME", "CONDITION", "FAILURE", "INJURY", "DAMAGE", "STRESS", "TRAUMA",
            "ACUTE", "CHRONIC", "SEVERE", "MILD", "MODERATE", "STAGE", "GRADE", "CLASS", "RISK", "FACTOR", "SCORE",
            "DIAGNOSIS", "PROGNOSIS", "THERAPY", "TREATMENT", "MANAGEMENT", "STRATEGY", "OUTCOME", "OUTCOMES",
            "EVENT", "EVENTS", "DEATH", "MORTALITY", "MORBIDITY", "SURVIVAL", "SAFETY", "EFFICACY", "QUALITY", "LIFE",
            "COVID", "COVID19", "SARS", "PANDEMIC", "VIRUS", "INFECTION", "SEPSIS", "CANCER", "TUMOR", "METASTASIS"
        ])

        # E. Exames, Estatística & Genéricos
        blacklist.update([
            "MRI", "CMR", "CT", "PET", "SPECT", "ECG", "EKG", "ECHO", "ULTRASOUND", "XRAY", "IMAGING", "CEUS",
            "WESTERN", "BLOT", "PCR", "QPCR", "RT", "ELISA", "STAINING", "IMMUNO", "HISTOLOGY", "ASSAY",
            "TOTAL", "MEAN", "RATIO", "SD", "SEM", "CI", "HR", "RR", "OR", "VS", "VERSUS", "LOG", "ANOVA",
            "YEAR", "YEARS", "MONTH", "MONTHS", "DAY", "DAYS", "HOUR", "HOURS", "MIN", "SEC", "TIME",
            "KG", "MG", "ML", "NM", "UM", "MM", "CM", "HZ", "MV", "MMHG", "AUC", "ROC", "P", "VALUE",
            "AND", "THE", "FOR", "NOT", "BUT", "WITH", "FROM", "THIS", "THAT", "THESE", "THOSE",
            "WHICH", "WHAT", "WHEN", "WHERE", "WHO", "WHY", "HOW", "ANY", "ALL", "EACH", "EVERY",
            "HAVE", "HAS", "HAD", "WAS", "WERE", "BEEN", "BEING", "ARE", "IS", "CAN", "COULD",
            "SHOULD", "WOULD", "MAY", "MIGHT", "MUST", "WILL", "SHALL", "DOES", "DID", "DOING",
            "VIA", "DUE", "BETWEEN", "AMONG", "WITHIN", "WITHOUT", "UNDER", "ABOVE", "BELOW", 
            "AFTER", "BEFORE", "DURING", "SINCE", "UNTIL", "WHILE", "ONCE", "UPON", "INTO", "ONTO",
            "USING", "USED", "USE", "DATA", "ANALYSIS", "STUDY", "TRIALS", "TRIAL", "META", "REVIEW"
        ])
        
        # F. Termos Biológicos Genéricos (Para deixar passar apenas os Específicos)
        blacklist.update([
            "DNA", "RNA", "MRNA", "MIRNA", "SIRNA", "GENE", "GENES", "PROTEIN", "PROTEINS", "CELL", "CELLS", 
            "TISSUE", "TISSUES", "RATS", "MICE", "RABBIT", "PIG", "DOG", "MODEL", "MODELS", "VIVO", "VITRO", 
            "EX", "SITU", "LEVEL", "LEVELS", "EXPRESSION", "REGULATION", "ACTIVITY", "ACTIVATION", "INHIBITION",
            "PATHWAY", "SIGNALING", "MECHANISM", "TARGET", "ROLE", "FUNCTION", "ACTION", "EFFECT", "EFFECTS",
            "MEDIATED", "INDUCED", "DEPENDENT", "INDEPENDENT", "ASSOCIATED", "RELATED", "LINKED", "POTENTIAL",
            "BIOMARKER", "CANDIDATE", "NOVEL", "NEW", "RECENT", "INSIGHTS", "UPDATE", "FUTURE", "PERSPECTIVE",
            "MOLECULAR", "GENETIC", "EPIGENETIC", "PHARMACOLOGY", "DRUG", "DRUGS", "RECEPTOR", "CHANNEL", "ENZYME"
        ])

        candidatos_por_artigo = []
        for artigo in artigos_raw:
            # SÓ LÊ TÍTULO E KEYWORDS
            texto_focado = ""
            for line in artigo.split("\n"):
                if line.startswith("TI  - "): texto_focado += line[6:].strip() + " "
                elif line.startswith("OT  - ") or line.startswith("KW  - "): texto_focado += line[6:].strip() + " "
            
            if not texto_focado: continue

            # REGEX OTIMIZADO:
            # 1. Siglas clássicas (TRPV1, mTOR, P2X3)
            # 2. Siglas de letras (TMAO, BDNF, TFEB) com min 3 letras
            # 3. Evita números isolados
            encontrados = re.findall(r'\b(?:[A-Z]{2,}[A-Z0-9-]*|[a-z]{1,2}[A-Z][a-zA-Z0-9-]*)\b', texto_focado)
            
            for t in encontrados:
                t_clean = re.sub(r'[^a-zA-Z0-9]', '', t).upper() # Normaliza tudo para UPPER na checagem
                
                if len(t_clean) < 3: continue
                if t_clean in blacklist: continue
                if t_clean == termo_upper.replace(" ", ""): continue
                if t_clean.isdigit(): continue 
                
                # Mantém a formatação original (ex: 'mTOR' em vez de 'MTOR') se possível, 
                # mas aqui vamos padronizar para evitar duplicatas.
                candidatos_por_artigo.append(t_clean)

        if not candidatos_por_artigo: return []
        
        # 1. Contagem Bruta
        contagem = Counter(candidatos_por_artigo)
        total_docs = max(1, len(artigos_raw))
        
        # 2. Filtro Estatístico Inicial (Lista Suja mas Rica)
        # Pegamos Top 80 termos que aparecem em <90% dos títulos
        top_candidatos = [termo for termo,freq in contagem.most_common(120) if (freq/total_docs)<0.90][:80]
        
        # 3. LIMPEZA FINAL VIA IA (Se o usuário tiver chave)
        if st.session_state.get('api_key_usuario'):
            # Chama o "Faxineiro IA"
            lista_limpa = _faxina_ia(top_candidatos)
            return lista_limpa # Retorna a lista curada pela IA
        else:
            # Sem chave, retorna a lista bruta filtrada (Top 30)
            return top_candidatos[:30]

    except: return []

# --- RADAR DE NOTÍCIAS ---
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
                                  
