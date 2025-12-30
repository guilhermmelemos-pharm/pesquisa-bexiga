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

# LISTA NEGRA: A Curadora Implacável
BLACKLIST_MONSTRO = {
    # 1. Termos de Imagem, Física e Estatística
    "MRI", "CT", "PET", "RBE", "IMPT", "RBE", "VIII", "RATIO", "CI", "HR", "OR",
    "SEER", "QSAR", "BMI", "ROC", "AUC", "ANOVA", "P-VALUE", "STATISTICS",
    
    # 2. Metodologia e Tipos de Estudo
    "ASSOCIATION", "EVALUATION", "UNCOMMON", "DISEASE", "BACKGROUND", 
    "OBJECTIVE", "METHODS", "RESULTS", "CONCLUSION", "ABSTRACT", 
    "INTRODUCTION", "STUDY", "ANALYSIS", "DATA", 
    "SIGNIFICANT", "DIFFERENCE", "BETWEEN", "AMONG", "WITHIN", 
    "DURING", "PREVALENCE", "INCIDENCE", "RISK", "FACTOR", "ROLE", 
    "POTENTIAL", "NOVEL", "DIAGNOSTIC", "ARTIFICIAL", "MANAGEMENT",
    "PROGNOSTIC", "FACTORS", "NEUROMODULATION", "IMPLICATIONS",
    "CLINICAL", "REVIEW", "META-ANALYSIS", "SYSTEMATIC", "SURVEY",
    "RCT", "COHORT", "TRIAL", "RETROSPECTIVE", "PROSPECTIVE", "MULTICENTER",
    "COMBINATION", "THERAPY", "TREATMENT", "EFFICACY", "SAFETY",
    
    # 3. Contexto Clínico, Anatomia e Doenças Irrelevantes
    "INTRAUTERINE", "GERMLINE", "BLADDER", "CANCER", "URINARY", 
    "UROTHELIAL", "MUSCLE", "OVERACTIVE", "TUMOR", "CARCINOMA", 
    "PELVIC", "URETHRAL", "UROLOGIC", "NEUROGENIC", "DIABETES", 
    "SJOGREN", "INTRAVESICAL", "VOID", "VOIDING", "DETRUSOR", 
    "PATIENT", "PATIENTS", "CHILDREN", "ADULT", "NHANES", "WOMEN", "MEN",
    "PROSTATE", "PARKINSON", "NMIBC", "PSA", "OAB", "COVID", "PFAS",
    "SYNDROME", "INJURY", "OBSTRUCTION", "HYPERPLASIA", "LUTS", "CAKUT",
    "DEPRESSION", "SUICIDE", "ADHD", "ATTENTION", "JAACAP", "CHRT",
    
    # 4. Biologia Genérica
    "DNA", "RNA", "ATP", "GENE", "PROTEIN", "CELL", "EXPRESSION",
    "BCG", "GUERIN", "CALMETTE", "BACILLUS", "HLA", "CAR",
    "PATHWAY", "MECHANISM", "RECEPTOR", "SIGNALING", "ACTIVITY",
    "LEVELS", "TISSUE", "SERUM", "PLASMA", "BLOOD", "URINE",
    "ERO1A", "ALOX5", "HDAC7", "CD8", "THBS1", "SOCS3", "ACTIVIN",
    "UBIQUITIN", "TRANSLATION", "YDSRN", "FGF",
    
    # 5. Stop Words e Conectivos
    "THE", "AND", "FOR", "NOT", "BUT", "VIA", "ALL", "WITH", "FROM", "AFTER",
    "RADIATION", "CHEMOTHERAPY", "SURGERY", "SHAM", "CONTROL", "GROUP",
    "PLACEBO", "SALINE", "DOSE", "DAILY", "WEEKLY", "MONTHLY"
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
            resp = requests.post(url, headers=headers, json=payload, timeout=25)
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

# ================= MINERAÇÃO ESTRUTURADA =================

def ner_extraction_batch(titulos_keywords: List[str], api_key: str, contexto_alvo: str) -> List[str]:
    """Extração focada APENAS em Fármacos do Alvo com Dedução PhD."""
    if not titulos_keywords: return []
    
    texto_input = "\n".join([f"- {t}" for t in titulos_keywords[:100]])
    
    # --- PROMPT V3.2: DOUTORA EM FARMACOLOGIA ---
    prompt = f"""
    Você é uma Doutora em Farmacologia (PhD).
    O seu foco de pesquisa atual é EXCLUSIVAMENTE: {contexto_alvo.upper()}.
    
    Analise a lista abaixo de termos extraídos da literatura recente.
    
    SUA TAREFA:
    Filtrar e retornar APENAS os Fármacos, Compostos e Alvos Moleculares que fazem sentido farmacológico para {contexto_alvo}.
    
    CRITÉRIOS DE DOUTORADO (Raciocínio Clínico):
    1. Se um fármaco é clássico de outra área (ex: Antidepressivo, Estatina) e não tem link claro com {contexto_alvo} no título, DESCARTE.
    2. Se for um termo genérico ("Therapy", "Protein"), DESCARTE.
    3. PRIORIZE: Novos compostos, inibidores específicos, canais iônicos e receptores.
    
    INPUT:
    {texto_input}
    
    OUTPUT: Apenas a lista JSON de strings (ex: ["Mirabegron", "P2X3", "Trealose"]).
    """
    return call_gemini_json(prompt, api_key)

@st.cache_data(ttl=3600)
def minerar_pubmed(termo_base: str, email: str, usar_ia: bool = True) -> Dict:
    api_key = st.session_state.get("api_key_usuario", "").strip()
    Entrez.email = email
    
    termo_expandido = MAPA_SINONIMOS_BASE.get(termo_base.upper(), termo_base)
    # Busca 300 artigos para ter volume de dados
    query = f"({termo_expandido}) AND (2020:2026[Date])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=300)
        record = Entrez.read(handle); handle.close()
        
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
        
        # 1. IA PhD: Filtra pelo alvo pedido
        if api_key and usar_ia:
            entidades = ner_extraction_batch(raw_texts, api_key, termo_base)
        
        # 2. REGEX: Rede de Segurança
        if True: 
            texto_full = " ".join(raw_texts)
            regex_codigos = r'\b[A-Z]{2,}[0-9]+[A-Z0-9-]*\b'
            entidades.extend(re.findall(regex_codigos, texto_full))
            
            sufixos = r'(?:ine|in|mab|ib|ol|on|one|il|ide|ate|ase|an)\b'
            regex_farmacos = r'\b[A-Z][a-z]{3,}' + sufixos
            entidades.extend(re.findall(regex_farmacos, texto_full))
            
            regex_acronimos = r'\b[A-Z]{3,}\b'
            entidades.extend(re.findall(regex_acronimos, texto_full))

        # 3. FILTRAGEM (Blacklist)
        entidades_limpas = []
        for e in entidades:
            e = e.strip(".,-;:()[] ")
            if len(e) < 3: continue 
            if e.isdigit(): continue
            if e.upper() in BLACKLIST_MONSTRO: continue
            if e.lower() in ["with", "from", "after", "during", "high", "low"]: continue
            entidades_limpas.append(e)

        counts = Counter(entidades_limpas)
        limit = 2 if len(counts) > 30 else 1
        recorrentes = sorted([e for e, c in counts.items() if c >= limit], key=lambda x: counts[x], reverse=True)
        
        final = []
        seen = set()
        for item in recorrentes:
            if item.lower() not in seen:
                final.append(item)
                seen.add(item.lower())

        return {
            "termos_indicados": final[:45],
            "counts": counts,
            "total_docs": len(artigos_completos),
            "artigos_originais": artigos_completos
        }
    except Exception: return {}

# ================= WRAPPERS (AQUI ESTAVA O ERRO DE ASSINATURA) =================
# Esta função agora aceita usar_ia, que o app está enviando
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
                "Info_IA": r.get("AB", "Resumo indisponível.")[:800], 
                "Link": f"https://pubmed.ncbi.nlm.nih.gov/{r.get('PMID', '')}/"
            })
        return artigos
    except: return []

# ================= ANÁLISE DETALHADA (COM DEDUÇÃO) =================

def analisar_abstract_com_ia(titulo: str, dados_curtos: str, api_key: str, lang: str = 'pt') -> str:
    """Gera a Análise Lemos Lambda com Dedução Farmacológica."""
    
    if api_key:
        prompt = f"""
        Você é Doutora em Farmacologia (PhD).
        Leia o título e resumo abaixo e DEDUZA o mecanismo farmacológico principal.
        
        TÍTULO: "{titulo}"
        RESUMO: "{dados_curtos}"
        
        SUA TAREFA:
        Preencha o trio: Órgão - Alvo - Ação.
        
        RACIOCÍNIO DE DEDUÇÃO (Use seu conhecimento de PhD):
        1. Se o artigo cita "Mirabegron", você SABE que o alvo é "Receptor Beta-3" e ação é "Agonista", mesmo que não esteja escrito. Use esse conhecimento!
        2. Se cita "Trealose", o alvo provável é "Agregação Proteica" ou "Autofagia".
        3. Se não houver fármaco, identifique a fisiopatologia.
        
        SAÍDA OBRIGATÓRIA (Apenas uma linha):
        Órgão/Tecido - Alvo (Receptor/Enzima) - Ação (Agonista/Inibidor/Expressão)
        
        Responda em Português. Se impossível deduzir, use N/A.
        """
        resposta_ia = simple_gemini_text(prompt, api_key)
        if resposta_ia:
            return resposta_ia.replace("\n", " ").strip().replace("Output:", "").replace("*", "")
    
    # Fallback
    palavras = titulo.split()
    palavras_chave = [p for p in palavras if len(p) > 5 and p.upper() not in BLACKLIST_MONSTRO]
    resumo_fallback = " | ".join(palavras_chave[:3])
    if not resumo_fallback: return "N/A - N/A - N/A"
    return f"Foco Provável: {resumo_fallback} (Sem IA)"

def buscar_todas_noticias(lang_code: str) -> List[Dict]:
    return []
