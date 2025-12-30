import streamlit as st
from Bio import Entrez, Medline
import requests
import json
import re
import time
from collections import Counter
from typing import List, Dict, Any

# ================= CONFIGURAÇÃO =================
Entrez.email = "pesquisador_guest@unifesp.br"

MODELOS_ATIVOS = [
    "gemini-2.5-flash", 
    "gemini-2.0-flash", 
    "gemini-1.5-flash"
]

# LISTA NEGRA: Bloqueio de Ruído Clínico e Metodológico
BLACKLIST_MONSTRO = {
    # 1. Termos Genéricos de Pesquisa
    "STUDY", "ANALYSIS", "REVIEW", "META-ANALYSIS", "DATA", "RESULTS", "CONCLUSION",
    "BACKGROUND", "METHODS", "OBJECTIVE", "AIM", "HYPOTHESIS", "INTRODUCTION",
    "SIGNIFICANT", "DIFFERENCE", "INCREASED", "DECREASED", "LEVELS", "EXPRESSION",
    "ROLE", "EFFECT", "IMPACT", "POTENTIAL", "NOVEL", "ASSOCIATION", "EVALUATION",
    
    # 2. Estatística e Métricas
    "P-VALUE", "ANOVA", "RATIO", "ODDS", "CONFIDENCE", "INTERVAL", "STATISTICS",
    "COHORT", "POPULATION", "SAMPLE", "SIZE", "BASELINE", "PREDICTION", "SIMULATION",
    "PRECISION", "ACCURACY", "SENSITIVITY", "SPECIFICITY", "SCORE", "SCALE", "INDEX",
    
    # 3. Procedimentos e Clínica (Onde costuma sujar)
    "SURGERY", "RESECTION", "INCISION", "OPERATION", "TRANSPLANT", "GRAFT", "STENT",
    "CATHETER", "BIOPSY", "IMAGING", "MRI", "CT", "PET", "ULTRASOUND", "DIAGNOSIS",
    "PROGNOSIS", "MANAGEMENT", "THERAPY", "TREATMENT", "PROTOCOL", "GUIDELINE",
    "COMPLICATION", "INFECTION", "DYSFUNCTION", "SYNDROME", "DISORDER", "DISEASE",
    "PATIENT", "PARTICIPANT", "CHILDREN", "ADULT", "WOMEN", "MEN", "ELDERLY",
    "HOSPITAL", "CLINIC", "CENTER", "DEPARTMENT", "UNIVERSITY",
    
    # 4. Termos Específicos que Vazaram Antes
    "PROSTATE", "KIDNEY", "LIVER", "HEART", "LUNG", "BRAIN" if "BLADDER" not in st.session_state.get('input_alvo', '').upper() else "", # Contexto dinâmico
    "COVID", "SARS-COV-2", "PANDEMIC", "VIRUS",
    "ENFORTUMAB", "PEMBROLIZUMAB", "NIVOLUMAB", "ATEROLIZUMAB", # Chemo clássica (a menos que seja o foco)
    "GEMCITABINE", "CISPLATIN", "CARBOPLATIN", "DOXORUBICIN",
    "PLACEBO", "CONTROL", "SHAM", "VEHICLE", "SALINE",
    "DNA", "RNA", "MRNA", "PROTEIN", "CELL", "TISSUE", "SERUM", "PLASMA", "URINE", "BLOOD"
}

MAPA_SINONIMOS_BASE = {
    "BLADDER": "(Bladder OR Urothelial OR Urothelium OR Detrusor)",
    "PAIN": "(Pain OR Nociception OR Analgesia)",
    "INFLAMMATION": "(Inflammation OR Cytokines OR NF-kappaB)"
}

# ================= GEMINI CORE =================

def clean_model_name(model_name: str) -> str:
    return model_name.replace("models/", "")

def call_gemini_json(prompt: str, api_key: str) -> List[str]:
    if not api_key: return []
    headers = {"Content-Type": "application/json"}
    payload = {
        "contents": [{"parts": [{"text": prompt}]}], 
        "generationConfig": {"temperature": 0.1, "response_mime_type": "application/json"}
    }
    for modelo_raw in MODELOS_ATIVOS:
        modelo = clean_model_name(modelo_raw)
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={api_key.strip()}"
            resp = requests.post(url, headers=headers, json=payload, timeout=30)
            if resp.status_code == 200:
                try:
                    text = resp.json()["candidates"][0]["content"]["parts"][0]["text"]
                    clean_text = re.sub(r"```json|```", "", text).strip()
                    parsed = json.loads(clean_text)
                    if isinstance(parsed, list): return [str(x) for x in parsed]
                    if isinstance(parsed, dict): 
                        for v in parsed.values(): 
                            if isinstance(v, list): return [str(x) for x in v]
                except: continue
            elif resp.status_code == 429: break 
        except: continue
    return []

def simple_gemini_text(prompt: str, api_key: str) -> str:
    if not api_key: return None
    headers = {"Content-Type": "application/json"}
    payload = {
        "contents": [{"parts": [{"text": prompt}]}], 
        "generationConfig": {"temperature": 0.3}
    }
    for modelo_raw in MODELOS_ATIVOS:
        modelo = clean_model_name(modelo_raw)
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={api_key.strip()}"
            resp = requests.post(url, headers=headers, json=payload, timeout=15)
            if resp.status_code == 200:
                return resp.json()["candidates"][0]["content"]["parts"][0]["text"]
        except: continue
    return None

# ================= MINERAÇÃO DE ALTA PRECISÃO =================

def ner_extraction_batch(textos_completos: List[str], api_key: str, contexto_alvo: str) -> List[str]:
    """
    Prompt 'Drug Hunter'. Recebe trechos ricos (Título + Abstract Parcial)
    e extrai apenas química pesada.
    """
    if not textos_completos: return []
    
    # Prepara o input para a IA (limitado para não estourar tokens)
    texto_input = "\n---\n".join(textos_completos[:50]) 
    
    # --- O PROMPT "DRUG HUNTER" ---
    prompt = f"""
    ROLE: Elite Pharmacologist & Chemical Data Curator.
    TARGET ORGAN/DISEASE: {contexto_alvo.upper()}.
    
    TASK: Scan the abstracts provided below. Extract ONLY:
    1. SPECIFIC DRUG NAMES (e.g., Mirabegron, Tadalafil, Resiniferatoxin).
    2. EXPERIMENTAL COMPOUNDS (codes like GYY4137, BAY-1234, AL-353).
    3. SPECIFIC MOLECULAR TARGETS (Receptors/Channels/Enzymes like P2X3, TRPV1, mTOR).
    
    CRITICAL EXCLUSION RULES (IGNORE THESE):
    - NO Clinical Procedures (Surgery, Resection, Injection).
    - NO Study Types (RCT, Review, Meta-analysis).
    - NO General Biological Terms (Gene, Protein, Cell, Pathway, Expression).
    - NO General Drug Classes without specific names (e.g., ignore "Antibiotics", "Chemotherapy").
    - NO "Standard of Care" drugs unless used in a NOVEL way (Ignore Cisplatin/Gemcitabine if used as standard control).
    
    YOUR OUTPUT FORMAT:
    A pure JSON list of strings. Example: ["Mirabegron", "P2X3", "GYY4137", "Trealose"]
    
    INPUT DATA:
    {texto_input}
    """
    return call_gemini_json(prompt, api_key)

@st.cache_data(ttl=3600)
def minerar_pubmed(termo_base: str, email: str, usar_ia: bool = True) -> Dict:
    api_key = st.session_state.get("api_key_usuario", "").strip()
    Entrez.email = email
    
    termo_expandido = MAPA_SINONIMOS_BASE.get(termo_base.upper(), termo_base)
    # Busca 200 artigos (menos que 300, mas com leitura mais profunda)
    query = f"({termo_expandido}) AND (drug OR inhibitor OR agonist OR antagonist OR compound) AND (2020:2026[Date])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=200)
        record = Entrez.read(handle); handle.close()
        
        if not record["IdList"]: return {}

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        
        raw_texts_for_ai = []
        artigos_completos = []
        
        for r in records:
            titulo = r.get('TI', '')
            abstract = r.get('AB', '')
            keywords = ' '.join(r.get('OT', []))
            
            # --- ESTRATÉGIA DE LEITURA INTELIGENTE ---
            # Pega o Título + Primeiros 250 caracteres (Objetivo) + Últimos 250 (Conclusão)
            # É aqui que os fármacos se escondem.
            intro = abstract[:250] if len(abstract) > 250 else abstract
            conclusao = abstract[-250:] if len(abstract) > 250 else ""
            
            texto_rico = f"TITLE: {titulo}\nABSTRACT_START: {intro}\nABSTRACT_END: {conclusao}\nKEYWORDS: {keywords}"
            
            raw_texts_for_ai.append(texto_rico)
            artigos_completos.append({"titulo": titulo, "texto": texto_rico})
        
        entidades = []
        
        # 1. IA DRUG HUNTER
        if api_key and usar_ia:
            entidades = ner_extraction_batch(raw_texts_for_ai, api_key, termo_base)
        
        # 2. REGEX DE APOIO (Focado em Química)
        if True: 
            texto_full = " ".join([a['texto'] for a in artigos_completos])
            
            # Códigos Experimentais (Letra+Numero) ex: GYY4137
            regex_codigos = r'\b[A-Z]{2,4}[- ]?[0-9]{3,5}\b'
            entidades.extend(re.findall(regex_codigos, texto_full))
            
            # Sufixos Farmacológicos
            sufixos = r'\b[A-Z][a-z]{3,}(?:ine|in|mab|ib|ol|on|one|il|ide|ate|ase|an)\b'
            entidades.extend(re.findall(sufixos, texto_full))
            
            # Alvos (Genes/Receptores em CAPS)
            regex_alvos = r'\b[A-Z0-9-]{3,8}\b'
            # Filtragem leve no regex bruto antes de adicionar
            candidatos = re.findall(regex_alvos, texto_full)
            candidatos = [c for c in candidatos if re.search(r'\d', c) or c.endswith("R")] # Ex: P2X3, TRPV1, EGFR
            entidades.extend(candidatos)

        # 3. FILTRAGEM FINAL
        entidades_limpas = []
        for e in entidades:
            e = e.strip(".,-;:()[] ")
            if len(e) < 3: continue 
            if e.isdigit(): continue
            if e.upper() in BLACKLIST_MONSTRO: continue
            
            # Remove palavras comuns do inglês que escapam
            if e.lower() in ["with", "from", "after", "during", "high", "low", "using", "treated"]: continue
            
            entidades_limpas.append(e)

        counts = Counter(entidades_limpas)
        limit = 1 # Como filtramos muito bem, até 1 ocorrência pode ser Blue Ocean
        recorrentes = sorted([e for e, c in counts.items() if c >= limit], key=lambda x: counts[x], reverse=True)
        
        final = []
        seen = set()
        for item in recorrentes:
            if item.lower() not in seen:
                final.append(item)
                seen.add(item.lower())

        return {
            "termos_indicados": final[:50],
            "counts": counts,
            "total_docs": len(artigos_completos),
            "artigos_originais": artigos_completos
        }
    except Exception as e: 
        print(f"Erro mineração: {e}")
        return {}

# ================= WRAPPERS =================
def buscar_alvos_emergentes_pubmed(alvo: str, email: str, usar_ia: bool = True) -> List[str]:
    res = minerar_pubmed(alvo, email, usar_ia=usar_ia)
    return res.get("termos_indicados", [])

def consultar_pubmed_count(termo: str, contexto: str, email: str, ano_ini: int, ano_fim: int) -> int:
    Entrez.email = email
    query = f"({termo}) AND ({contexto}) AND ({ano_ini}:{ano_fim}[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle); handle.close()
        return int(record["Count"])
    except: return 0

@st.cache_data(ttl=3600)
def buscar_resumos_detalhados(termo: str, orgao: str, email: str, ano_ini: int, ano_fim: int) -> List[Dict]:
    Entrez.email = email
    query = f"({termo}) AND ({orgao}) AND ({ano_ini}:{ano_fim}[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=6)
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        artigos = []
        for r in records:
            artigos.append({
                "Title": r.get("TI", "Sem Título"), 
                "Info_IA": r.get("AB", "Resumo indisponível.")[:1000], 
                "Link": f"https://pubmed.ncbi.nlm.nih.gov/{r.get('PMID', '')}/"
            })
        return artigos
    except: return []

# ================= ANÁLISE DETALHADA =================

def analisar_abstract_com_ia(titulo: str, dados_curtos: str, api_key: str, lang: str = 'pt') -> str:
    if api_key:
        prompt = f"""
        Você é Doutora em Farmacologia (PhD).
        Leia o texto e DEDUZA o mecanismo.
        TÍTULO: "{titulo}"
        RESUMO: "{dados_curtos}"
        
        SAÍDA (Uma linha): Órgão - Alvo - Ação (Fármaco)
        Se não houver fármaco explícito, deduza pelo mecanismo.
        Ex: Bexiga - Receptor Beta-3 - Agonista (Mirabegron)
        """
        resposta_ia = simple_gemini_text(prompt, api_key)
        if resposta_ia:
            return resposta_ia.replace("\n", " ").strip().replace("Output:", "").replace("*", "")
    return "Análise pendente."

def buscar_todas_noticias(lang_code: str) -> List[Dict]:
    return []
