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

# --- 2. FAXINEIRO IA (PROMPT INOVADOR) ---
def _faxina_ia(lista_suja):
    api_key = st.session_state.get('api_key_usuario', '').strip()
    if not api_key: return lista_suja[:30] 

    lista_str = ", ".join(lista_suja)
    
    # Prompt focado em achar "A Agulha no Palheiro" (TMAO, TFEB, Trehalose)
    prompt = f"""
    ACT AS: PhD Researcher in Pharmacophysiology searching for NOVEL targets.
    
    INPUT LIST: {lista_str}
    
    TASK: Curate this list. Remove generic biology terms and keep specific pharmacological interests.
    
    ✅ KEEP (PRIORITY):
    - Novel Metabolites & Small Molecules: TMAO, Short-chain fatty acids, Resveratrol, Trehalose.
    - Specific Regulators: TFEB (Autophagy), mTOR, HIF-1a, NRF2.
    - Specific Receptors/Channels: TRPV1, P2X3, Beta-3-AR, Piezo1.
    - Experimental Drugs/Compounds: Mirabegron, Solifenacin, C29, GYY4137.
    
    ❌ DELETE (NOISE):
    - Ubiquitous Molecules: ATP, DNA, RNA, cAMP, VEGF.
    - Clinical/Diseases: LUTS, OAB, Cancer, Infection, COVID, UTI, BPH.
    - Procedures/Diagnostics: MRI, TURBT, Botox, EMG, Urodynamics.
    - General Terms: Protein, Gene, Study, Analysis, Group, Rat, Mouse.
    
    OUTPUT: Return strictly a Python list of strings.
    """
    
    headers = {'Content-Type': 'application/json'}
    data = {
        "contents": [{"parts": [{"text": prompt}]}],
        "safetySettings": [{"category": "HARM_CATEGORY_DANGEROUS_CONTENT", "threshold": "BLOCK_NONE"}],
        "generationConfig": {"temperature": 0.0}
    }
    
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=20) # Timeout maior
            
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                texto = texto.replace("```python", "").replace("```json", "").replace("```", "").strip()
                if texto.startswith("[") and texto.endswith("]"):
                    lista_limpa = ast.literal_eval(texto)
                    if isinstance(lista_limpa, list) and len(lista_limpa) > 0:
                        return lista_limpa
        except: continue
    
    return lista_suja[:30]

# --- 3. ANÁLISE DE RESUMOS ---
def analisar_abstract_com_ia(titulo, dados_curtos, api_key, lang='pt'):
    key = api_key.strip()
    if not key: return "⚠️ Chave API não detectada."
    idioma = "Português" if lang == 'pt' else "Inglês"
    
    prompt_text = f"Resuma Alvo/Fármaco e Mecanismo em 15 palavras. {titulo}. {dados_curtos}. Idioma: {idioma}."
    
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt_text}]}]}
    
    ultimo_erro = ""

    for m in MODELOS_ATIVOS:
        try:
            url = f"[https://generativelanguage.googleapis.com/v1beta/models/](https://generativelanguage.googleapis.com/v1beta/models/){m}:generateContent?key={key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=10)
            if resp.status_code == 200:
                return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
            else:
                ultimo_erro = f"HTTP {resp.status_code}"
        except Exception as e:
            ultimo_erro = str(e)
            continue
            
    return f"⚠️ Erro IA: {ultimo_erro}"

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
    
    # Busca 2000 artigos para garantir que pegamos termos raros como Trehalose
    final_query = f"({query_string}) AND (2018:2030[Date - Publication]) AND (NOT Review[pt])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=final_query, retmax=2000, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        artigos_raw = full_data.split("\n\nPMID-")

        # --- LISTA NEGRA DE EXTERMÍNIO ---
        # Removendo termos genéricos que você não quer ver
        blacklist_exterminio = {
            # Genéricos Biológicos (Feijão com Arroz)
            "ATP", "DNA", "RNA", "NO", "ROS", "CO2", "H2O", "PGE", "PGE2", "VEGF", "PD1", "PDL1", "HER2",
            "UPEC", "HBP", "SGC", "PDE5", "ALK", "TP63", "DHEA", "CD44", "IFN", "MIRNAS", "PRP", "SVF",
            "RRNA", "UCA", "NGAL", "EVA", "BTX", "BOTULINUM", "TOXIN", "ANTIBODIES", "NEUTRALIZING", 
            "GPA", "PACAP", "PAC1", "EG70", "CG0070", "CD90", "CXCL13", "AVP", "AQP2", "HCG", "DMSO",
            
            # Clínicos & Doenças (Sua lista anterior)
            "LUTS", "OAB", "BPH", "UTI", "IC", "BPS", "ICIRS", "LUT", "LUTD", "BOO", "SUI", "UUI", "MUI",
            "COVID", "COVID19", "SARS", "VIRUS", "INFECTION", "SEPSIS", "CANCER", "TUMOR", "CIS", "MIBC", "NMIBC",
            
            # Procedimentos & Diagnóstico
            "BCG", "TURBT", "TURP", "SNM", "PTNS", "BOTOX", "INJECTION", "STENT", "CATHETER", "MRI", "CT", "PET",
            "VIRADS", "PIRADS", "CEUS", "ULTRASOUND", "EMG", "URODYNAMICS", "ACR", "AUC", "ROC", "CI", "OR",
            
            # Lixo Gramatical
            "THE", "AND", "WITH", "FOR", "BUT", "NOT", "FROM", "USED", "USING", "DATA", "STUDY", "RESULTS",
            "UNIVERSITY", "DEPARTMENT", "CENTER", "HOSPITAL", "PUBLISH", "ACCEPTED", "RECEIVED", "PMC", "PMID",
            "WESTERN", "BLOT", "PCR", "ELISA", "ANALYSIS", "REVIEW", "META", "TRIAL", "COHORT", "GROUP", "CONTROL",
            "TOTAL", "MEAN", "RATIO", "YEAR", "MONTH", "DAY", "HOUR", "MIN", "SEC", "HIGH", "LOW", "LEVEL"
        }

        candidatos_por_artigo = []
        for artigo in artigos_raw:
            texto_focado = ""
            for line in artigo.split("\n"):
                if line.startswith("TI  - "): texto_focado += line[6:].strip() + " "
                elif line.startswith("OT  - ") or line.startswith("KW  - "): texto_focado += line[6:].strip() + " "
            
            if not texto_focado: continue

            # --- REGEX NOVO: Captura Siglas E Palavras estilo "Chemical" (Trehalose, Resveratrol) ---
            # 1. [A-Z]{2,} -> Siglas (TRPV1)
            # 2. [A-Z][a-z]{3,} -> Title Case (Trehalose, Metformin)
            encontrados = re.findall(r'\b(?:[A-Z]{2,}[A-Z0-9-]*|[A-Z][a-z]{3,}[a-z0-9-]*)\b', texto_focado)
            
            for t in encontrados:
                # Removemos plural simples ('s) e pontuação
                t_clean = re.sub(r'[^a-zA-Z0-9]', '', t)
                t_upper = t_clean.upper()
                
                # Regra 1: Tamanho mínimo (TMAO passa, NO morre)
                if len(t_clean) < 3: continue 
                
                # Regra 2: Blacklist (Verifica em Maiúsculo para garantir)
                if t_upper in blacklist_exterminio: continue
                if t_upper == termo_upper.replace(" ", ""): continue
                
                # Regra 3: Se for Title Case (Ex: Trehalose), mantemos a formatação original na lista
                # para que a IA entenda que é nome próprio/químico.
                candidatos_por_artigo.append(t_clean)

        if not candidatos_por_artigo: return []
        
        contagem = Counter(candidatos_por_artigo)
        total_docs = max(1, len(artigos_raw))
        
        # Pega uma amostra MAIOR (Top 200) para incluir os termos raros (TFEB, Trehalose)
        top_candidatos = [termo for termo,freq in contagem.most_common(200) if (freq/total_docs)<0.90]
        
        # Filtra redundante (caso regex tenha pegado algo da blacklist com casing diferente)
        lista_para_ia = [t for t in top_candidatos if t.upper() not in blacklist_exterminio]
        
        if usar_ia and st.session_state.get('api_key_usuario'):
            # Envia para a IA filtrar
            return _faxina_ia(lista_para_ia)
        else:
            return lista_para_ia[:40]

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
                    "link": f"[https://pubmed.ncbi.nlm.nih.gov/](https://pubmed.ncbi.nlm.nih.gov/){pmid}/", 
                    "img": "[https://images.unsplash.com/photo-1532187863486-abf9dbad1b69?w=400](https://images.unsplash.com/photo-1532187863486-abf9dbad1b69?w=400)"
                })
        return news
    except: return []
