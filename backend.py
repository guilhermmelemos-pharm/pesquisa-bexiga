import streamlit as st
from Bio import Entrez
import requests
import json
import re
from collections import Counter
import ast
import math

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br"

# --- 1. MAPA DE EXPANSÃO SEMÂNTICA ---
MAPA_SINONIMOS = {
    "BLADDER": "Bladder OR Urothelium OR Detrusor OR Vesical OR Urethra OR Micturition OR LUTS OR Cystitis OR Overactive Bladder OR OAB OR Urinary Tract",
    "PAIN": "Pain OR Nociception OR Analgesia OR Neuropathic OR Hyperalgesia OR Dorsal Root Ganglion OR Allodynia OR TRP Channels",
    "INFLAMMATION": "Inflammation OR Cytokine OR Macrophage OR Neutrophil OR Immune OR Sepsis OR Inflammasome OR T-cell OR NF-kB",
    "METABOLISM": "Metabolism OR Obesity OR Diabetes OR Insulin OR Glucose OR Adipose OR Lipid OR Metabolic Syndrome OR Mitochondria"
}

# --- SEUS MODELOS (PAID TIER) ---
MODELOS_ATIVOS = [
    "gemini-2.5-flash", 
    "gemini-2.0-flash", 
    "gemini-2.0-flash-exp", 
    "gemini-flash-latest", 
    "gemini-2.5-flash-lite"
]

# --- 2. FAXINEIRO IA (CURADORIA PhD) ---
def _faxina_ia(lista_suja):
    api_key = st.session_state.get('api_key_usuario', '').strip()
    if not api_key: return lista_suja[:40] 

    lista_str = ", ".join(lista_suja)
    prompt = f"""
    ROLE: Senior Scientist in Pharmacology & Molecular Biology.
    INPUT: {lista_str}
    TASK: Strictly filter this list for BENCH RESEARCH relevance.
    
    RULES:
    1. KEEP: Specific receptors (TRPV4, P2X3), enzymes (SPHK1, ROCK), signaling molecules (cAMP, NO), and specific bench drugs (GSK1016790A, Ouabain).
    2. DELETE: Any clinical term, surgery name, disease symptom, or general anatomy.
    3. DELETE: Standard clinical treatments like Botox, BCG, or generic classes like 'Anticholinergics'.
    
    OUTPUT: Return strictly a Python list of strings.
    """
    
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": 0.0}}

    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=15)
            if resp.status_code == 200:
                texto = resp.json()['candidates'][0]['content']['parts'][0]['text']
                texto = texto.replace("```python", "").replace("```", "").strip()
                return ast.literal_eval(texto)
        except: continue
    return lista_suja[:35]

# --- 3. MINERAÇÃO MASSIVA COM FILTRO RÍGIDO ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    termo_upper = termo_base.upper().strip()
    query_string = MAPA_SINONIMOS.get(termo_upper, f"{termo_base}[Title/Abstract]")
    
    final_query = f"({query_string}) AND (2018:2030[Date - Publication]) AND (Pharmacology[Filter] OR Molecular[Filter])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=final_query, retmax=1500, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()

        # SUPER BLACKLIST (Foco em Bancada)
        blacklist = {
            "DIABETES", "CANCER", "TRIAL", "PATIENT", "PATIENTS", "CLINICAL", "SYNDROME", "DIAGNOSIS", 
            "TREATMENT", "MANAGEMENT", "SURGERY", "SURGICAL", "OUTCOME", "OUTCOMES", "EFFICACY", "SAFETY",
            "PLACEBO", "RANDOMIZED", "CYSTITIS", "INCONTINENCE", "LUTS", "OAB", "BPH", "UTI", "ICBPS",
            "STUDY", "REVIEW", "ANALYSIS", "DATA", "RESULTS", "CONCLUSION", "OBJECTIVE", "METHODS",
            "BACKGROUND", "PURPOSE", "HYPOTHESIS", "EVALUATION", "ASSESSMENT", "COMPARISON", "ROLE",
            "NOVEL", "NEW", "POTENTIAL", "EFFECT", "EFFECTS", "IMPACT", "DURING", "AFTER", "LEVEL",
            "EXPRESSION", "PATHWAY", "MECHANISM", "ACTIVITY", "FUNCTION", "ASSOCIATION", "THERAPY",
            "BLADDER", "URINARY", "TRACT", "URETHRA", "URETER", "KIDNEY", "RENAL", "DETRUSOR", 
            "UROTHELIUM", "UROTHELIAL", "MUCOSA", "MUSCLE", "SMOOTH", "TISSUE", "CELLS", "CELL", 
            "ORGANOIDS", "BODY", "HUMAN", "MALE", "FEMALE", "ANIMAL", "MODEL", "RAT", "MOUSE",
            "BCG", "BACILLUS", "CALMETTE", "GUERIN", "BOTOX", "ONABOTULINUM", "BTX", "TOXIN",
            "MICROBIOME", "MICROBIOTA", "BACTERIA", "UPEC", "INFECTION", "IMMUNE", "ANTIBODY",
            "ANTICHOLINERGIC", "ANTICHOLINERGICS", "ANTIMUSCARINIC", "ANTIMUSCARINICS",
            "DNA", "RNA", "mRNA", "VEGF", "ROS", "ATP", "CO2", "CALCIUM", "GENE", "PROTEIN"
        }

        candidatos = []
        for artigo in full_data.split("\n\nPMID-"):
            texto = ""
            for line in artigo.split("\n"):
                if line.startswith("TI  - ") or line.startswith("KW  - ") or line.startswith("OT  - "):
                    texto += line[6:].strip() + " "
            
            if not texto: continue
            encontrados = re.findall(r'\b(?:[A-Z]{2,}[A-Z0-9-]*|[A-Z][a-z]{4,}[a-z0-9-]*)\b', texto)
            
            for t in encontrados:
                t_clean = re.sub(r'[^a-zA-Z0-9]', '', t)
                t_up = t_clean.upper()
                if len(t_clean) < 3 or t_up in blacklist: continue
                candidatos.append(t_clean)

        if not candidatos: return []
        contagem = Counter(candidatos)
        top_candidatos = [termo for termo,freq in contagem.most_common(150)]
        
        if usar_ia and st.session_state.get('api_key_usuario'):
            return _faxina_ia(top_candidatos)
        else:
            return top_candidatos[:40]
    except: return []

# --- 4. ANÁLISE DETALHADA ---
def analisar_abstract_com_ia(titulo, dados_curtos, api_key, lang='pt'):
    api_key = api_key.strip()
    if not api_key: return "Chave API necessária."
    
    prompt = f"Explique o mecanismo farmacológico/fisiológico de {titulo} em 15 palavras. Foco em bancada e molecular."
    headers = {'Content-Type': 'application/json'}
    data = {"contents": [{"parts": [{"text": prompt}]}]}
    
    # Tentativa sequencial nos seus modelos Paid Tier
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, headers=headers, data=json.dumps(data), timeout=12)
            if resp.status_code == 200:
                return resp.json()['candidates'][0]['content']['parts'][0]['text'].strip()
        except: continue
        
    return "Erro na análise (IA indisponível)."

# Funções de contagem mantidas para o radar e métricas
@st.cache_data(ttl=86400, show_spinner=False)
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    query = f"({termo}) AND ({contexto}) AND ({ano_ini}:{ano_fim}[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle); handle.close()
        return int(record["Count"])
    except: return 0
