# backend.py
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

# LISTA NEGRA: Anatomia e Termos Clínicos (O que estava sujando sua lista)
BLACKLIST_MONSTRO = {
    # Anatomia e Doença (O que você viu e não gostou)
    "BLADDER", "CANCER", "URINARY", "UROTHELIAL", "MUSCLE", "OVERACTIVE", "TUMOR",
    "CARCINOMA", "PELVIC", "URETHRAL", "UROLOGIC", "NEUROGENIC", "DIABETES",
    "SJOGREN", "INTRAVESICAL", "VOID", "VOIDING", "DETRUSOR", "KIDNEY", "RENAL",
    
    # Metodologia e Tipos de Estudo
    "NOVEL", "DIAGNOSTIC", "ARTIFICIAL", "MANAGEMENT", "LAPAROSCOPIC", "PROGNOSTIC",
    "FACTORS", "NEOADJUVANT", "RANDOMIZED", "TYPE", "COHORT", "SURGICAL", 
    "ASSOCIATIONS", "SYMPTOMS", "IMPLICATIONS", "ONCOLOGY", "TRIAL", "RETROSPECTIVE",
    "EVALUATION", "REVIEW", "STUDY", "ANALYSIS", "DATA", "RESULTS", "CONCLUSION",
    
    # Grupos e Demografia
    "CHILDREN", "PATIENT", "PATIENTS", "WOMEN", "MEN", "ADULT", "NHANES",
    
    # Nomes Próprios/Históricos (Que não são alvos moleculares)
    "GUERIN", "CALMETTE", "BACILLUS", "BCG", 
    
    # Biologia Genérica
    "DNA", "RNA", "ATP", "GENE", "PROTEIN", "CELL", "EXPRESSION"
}

MAPA_SINONIMOS_BASE = {
    "BLADDER": "(Bladder OR Urothelial OR Urothelium)",
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
            resp = requests.post(url, headers=headers, json=payload, timeout=20)
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
    if not api_key: return ""
    headers = {"Content-Type": "application/json"}
    payload = {
        "contents": [{"parts": [{"text": prompt}]}], 
        "generationConfig": {"temperature": 0.2}
    }
    for modelo_raw in MODELOS_ATIVOS:
        modelo = clean_model_name(modelo_raw)
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={api_key.strip()}"
            resp = requests.post(url, headers=headers, json=payload, timeout=20)
            if resp.status_code == 200:
                return resp.json()["candidates"][0]["content"]["parts"][0]["text"]
        except: continue
    return ""

# ================= MINERAÇÃO ESTRUTURADA =================

def ner_extraction_batch(titulos_keywords: List[str], api_key: str) -> List[str]:
    """Prompt específico para separar o joio do trigo."""
    if not titulos_keywords: return []
    
    # Enviamos lotes pequenos de títulos
    texto_input = "\n".join([f"- {t}" for t in titulos_keywords[:60]])
    
    prompt = f"""
    You are a Pharmacologist. Analyze these paper titles.
    
    TASK: Extract strictly MOLECULAR TARGETS (Receptors, Enzymes, Channels, Genes) and DRUGS/COMPOUNDS.
    
    STRICT RULES:
    1. IGNORE Anatomy (Bladder, Muscle, Urothelial).
    2. IGNORE Disease types (Cancer, Carcinoma, Tumor).
    3. IGNORE Study types (Cohort, Trial, Analysis).
    4. EXTRACT ONLY Specific Targets (e.g., "HSP90", "P2X3", "mTOR") and Drugs (e.g., "Cisplatin", "Mirabegron").
    
    INPUT:
    {texto_input}
    
    OUTPUT: JSON list of strings (Target/Drug names only).
    """
    return call_gemini_json(prompt, api_key)

@st.cache_data(ttl=3600)
def minerar_pubmed(termo_base: str, email: str, usar_ia: bool = True) -> Dict:
    api_key = st.session_state.get("api_key_usuario", "").strip()
    Entrez.email = email
    
    termo_expandido = MAPA_SINONIMOS_BASE.get(termo_base.upper(), termo_base)
    query = f"({termo_expandido}) AND (2020:2026[Date])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=100)
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]: return {}

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        
        raw_texts = []
        artigos_completos = []
        
        for r in records:
            titulo = r.get('TI', '')
            keywords = ' '.join(r.get('OT', []))
            texto_limpo = f"{titulo} . {keywords}"
            if len(texto_limpo) > 5:
                raw_texts.append(texto_limpo)
                artigos_completos.append({"titulo": titulo, "texto": texto_limpo})
        
        entidades = []
        
        # 1. IA (Primeira linha de defesa)
        if api_key and usar_ia:
            entidades = ner_extraction_batch(raw_texts, api_key)
        
        # 2. Regex Seguro (Só entra se a IA falhar ou para complementar)
        # REMOVIDO: O regex genérico que pegava "Bladder" e "Cancer"
        if True:
            texto_full = " ".join(raw_texts)
            
            # A: Códigos (Letras maiúsculas + Números) -> P2X3, HSP90
            regex_codigos = r'\b[A-Z]{2,}[0-9]+[A-Z0-9-]*\b'
            candidatos_codigos = re.findall(regex_codigos, texto_full)
            
            # B: Fármacos via Sufixos Químicos (Segurança total)
            # Pega palavras PascalCase que terminam com sufixos de drogas
            sufixos = r'(?:ine|in|mab|ib|ol|on|one|il|ide|ate|ase)\b'
            regex_farmacos = r'\b[A-Z][a-z]{3,}' + sufixos
            candidatos_farmacos = re.findall(regex_farmacos, texto_full)
            
            # C: Acrônimos de Vias (mTOR, VEGF) - 3+ Letras Maiúsculas
            regex_acronimos = r'\b[A-Z]{3,}\b'
            candidatos_acronimos = re.findall(regex_acronimos, texto_full)

            entidades.extend(candidatos_codigos)
            entidades.extend(candidatos_farmacos)
            entidades.extend(candidatos_acronimos)

        # 3. Filtragem Final
        entidades_limpas = []
        for e in entidades:
            e = e.strip(".,-;:()[] ")
            if len(e) < 3: continue 
            if e.isdigit(): continue
            
            # Blacklist reforçada
            if e.upper() in BLACKLIST_MONSTRO: continue
            
            # Filtro extra para preposições
            if e.upper() in ["THE", "AND", "FOR", "NOT", "BUT", "VIA", "ALL"]: continue
            
            entidades_limpas.append(e)

        # 4. Contagem
        counts = Counter(entidades_limpas)
        
        # Exigimos pelo menos 2 ocorrências para limpar ruído (exceto se a lista for muito curta)
        limit = 2 if len(counts) > 20 else 1
        recorrentes = [e for e, c in counts.items() if c >= limit]
        
        # Ordenar por frequência
        recorrentes = sorted(recorrentes, key=lambda x: counts[x], reverse=True)
        
        # Deduplicação
        final = []
        seen = set()
        for item in recorrentes:
            if item.lower() not in seen:
                final.append(item)
                seen.add(item.lower())

        return {
            "termos_indicados": final[:40],
            "counts": counts,
            "total_docs": len(artigos_completos),
            "artigos_originais": artigos_completos
        }
    except Exception:
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
        record = Entrez.read(handle)
        handle.close()
        return int(record["Count"])
    except: return 0

@st.cache_data(ttl=3600)
def buscar_resumos_detalhados(termo: str, orgao: str, email: str, ano_ini: int, ano_fim: int) -> List[Dict]:
    Entrez.email = email
    query = f"({termo}) AND ({orgao}) AND ({ano_ini}:{ano_fim}[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=6)
        record = Entrez.read(handle)
        handle.close()
        if not record["IdList"]: return []
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        artigos = []
        for r in records:
            artigos.append({
                "Title": r.get("TI", "Sem Título"), 
                "Info_IA": r.get("AB", "Resumo indisponível.")[:800], 
                "Link": f"https://pubmed.ncbi.nlm.nih.gov/{r.get('PMID', '')}/"
            })
        return artigos
    except: return []

def analisar_abstract_com_ia(titulo: str, dados_curtos: str, api_key: str, lang: str = 'pt') -> str:
    """Aqui é onde entra a 'Função' que você pediu, sob demanda."""
    if not api_key: return "Chave API necessária."
    # Prompt explícito para categorizar
    prompt = f"Analyze: {titulo}\nOutput format: [TARGET/DRUG] | [MECHANISM/FUNCTION]"
    return simple_gemini_text(prompt, api_key).replace("\n", " ")

def buscar_todas_noticias(lang_code: str) -> List[Dict]:
    Entrez.email = "pesquisador_guest@unifesp.br"
    query = "(molecular pharmacology) AND (bladder) AND (2025:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=4, sort="pub_date")
        record = Entrez.read(handle)
        handle.close()
        if not record["IdList"]: return []
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        news = []
        for r in records:
            news.append({
                "titulo": r.get("TI", "Novo Artigo"), 
                "fonte": r.get("JT", "Journal")[:20], 
                "link": f"https://pubmed.ncbi.nlm.nih.gov/{r.get('PMID', '')}/"
            })
        return news
    except: return []
