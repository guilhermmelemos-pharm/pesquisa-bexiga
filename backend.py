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
    "SYNDROME", "INJURY", "OBSTRUCTION", "HYPERPLASIA", "LUTS",
    
    # 4. Biologia Genérica (Termos que não são fármacos específicos)
    "DNA", "RNA", "ATP", "GENE", "PROTEIN", "CELL", "EXPRESSION",
    "BCG", "GUERIN", "CALMETTE", "BACILLUS", "HLA", "CAR",
    "PATHWAY", "MECHANISM", "RECEPTOR", "SIGNALING", "ACTIVITY",
    "LEVELS", "TISSUE", "SERUM", "PLASMA", "BLOOD", "URINE",
    "ERO1A", "ALOX5", "HDAC7", "CD8", "THBS1", # Genes genéricos que vazaram antes
    
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
    """Retorna lista JSON. Se falhar, retorna lista vazia."""
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
            resp = requests.post(url, headers=headers, json=payload, timeout=25) # Timeout maior para lista longa
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

def ner_extraction_batch(titulos_keywords: List[str], api_key: str) -> List[str]:
    """Extração focada APENAS em Fármacos e Alvos Específicos."""
    if not titulos_keywords: return []
    
    # Aumentei o input para a IA ter mais contexto
    texto_input = "\n".join([f"- {t}" for t in titulos_keywords[:80]])
    
    # --- PROMPT ATUALIZADO (MODO PESQUISADORA) ---
    prompt = f"""
    Atue como uma Pesquisadora Sênior em Farmacologia.
    Você recebeu esta lista de títulos e keywords de artigos recentes.
    
    Sua missão: Identificar FÁRMACOS (Drogas, Compostos, Moléculas Pequenas) e ALVOS MOLECULARES ESPECÍFICOS.
    
    REGRAS RÍGIDAS:
    1. NÃO INVENTE. Se não houver fármaco no título, ignore.
    2. SELECIONE APENAS OS TERMOS QUÍMICOS/MOLECULARES.
    3. Retorne APENAS a lista de nomes (JSON Array).
    
    O QUE IGNORAR (LIXO):
    - Termos genéricos (Study, Effect, Patient, Bladder, Cancer, Therapy).
    - Metodologias (RCT, Review, Analysis).
    - Doenças (Cistite, Dor, Inflamação).
    
    O QUE BUSCAR (BLUE OCEAN):
    - Fármacos Específicos: Mirabegron, Tadalafil, Resiniferatoxin, Trealose, TMAO.
    - Alvos Moleculares: P2X3, TRPV1, mTOR, NGF, BDNF, Piezo1.
    - Compostos Experimentais (Siglas): GYY4137, AL-353, BAY-1234.
    
    INPUT:
    {texto_input}
    
    OUTPUT JSON:
    """
    return call_gemini_json(prompt, api_key)

@st.cache_data(ttl=3600)
def minerar_pubmed(termo_base: str, email: str, usar_ia: bool = True) -> Dict:
    api_key = st.session_state.get("api_key_usuario", "").strip()
    Entrez.email = email
    
    termo_expandido = MAPA_SINONIMOS_BASE.get(termo_base.upper(), termo_base)
    # Janela de tempo focada em novidades
    query = f"({termo_expandido}) AND (2020:2026[Date])"
    
    try:
        # --- AUMENTO DE ESCOPO: 300 ARTIGOS ---
        handle = Entrez.esearch(db="pubmed", term=query, retmax=300)
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]: return {}

        # Busca metadados
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        
        raw_texts = []
        artigos_completos = []
        
        for r in records:
            titulo = r.get('TI', '')
            keywords = ' '.join(r.get('OT', []))
            # Combinando título e keywords para dar contexto à IA
            texto_limpo = f"{titulo} . {keywords}"
            if len(texto_limpo) > 5:
                raw_texts.append(texto_limpo)
                artigos_completos.append({"titulo": titulo, "texto": texto_limpo})
        
        entidades = []
        
        # 1. IA: O Cérebro da Pesquisadora
        if api_key and usar_ia:
            # Enviamos em batches maiores ou o top relevante
            entidades = ner_extraction_batch(raw_texts, api_key)
        
        # 2. REGEX: A Rede de Segurança (Garante que códigos como P2X3 não escapem)
        if True: 
            texto_full = " ".join(raw_texts)
            
            # Códigos (P2X3, GYY4137)
            regex_codigos = r'\b[A-Z]{2,}[0-9]+[A-Z0-9-]*\b'
            entidades.extend(re.findall(regex_codigos, texto_full))
            
            # Fármacos (Sufixos Químicos: -ine, -mab, -ol, -an, -ib)
            sufixos = r'(?:ine|in|mab|ib|ol|on|one|il|ide|ate|ase|an)\b'
            regex_farmacos = r'\b[A-Z][a-z]{3,}' + sufixos
            entidades.extend(re.findall(regex_farmacos, texto_full))
            
            # Acrônimos Fortes (VEGF, mTOR, BDNF)
            regex_acronimos = r'\b[A-Z]{3,}\b'
            entidades.extend(re.findall(regex_acronimos, texto_full))

        # 3. FILTRAGEM (A Blacklist Curadora)
        entidades_limpas = []
        for e in entidades:
            e = e.strip(".,-;:()[] ")
            if len(e) < 3: continue 
            if e.isdigit(): continue
            
            # O Grande Filtro
            if e.upper() in BLACKLIST_MONSTRO: continue
            if e.lower() in ["with", "from", "after", "during", "high", "low"]: continue
            
            entidades_limpas.append(e)

        # 4. Contagem e Seleção
        counts = Counter(entidades_limpas)
        # Se temos muitos dados (300 artigos), exigimos que apareça pelo menos 2x para ser relevante
        # Se for muito raro (Blue Ocean extremo), a IA pode ter pego, então relaxamos se a lista for curta
        limit = 2 if len(counts) > 30 else 1
        
        recorrentes = [e for e, c in counts.items() if c >= limit]
        recorrentes = sorted(recorrentes, key=lambda x: counts[x], reverse=True)
        
        # Deduplicação Final
        final = []
        seen = set()
        for item in recorrentes:
            if item.lower() not in seen:
                final.append(item)
                seen.add(item.lower())

        return {
            "termos_indicados": final[:45], # Retorna Top 45
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
    """Gera a Análise Lemos Lambda (Órgão - Alvo - Ação)."""
    
    # 1. Tenta via IA (Gemini)
    if api_key:
        prompt = f"""
        Atue como Farmacologista PhD.
        Analise o título: "{titulo}" e o resumo: "{dados_curtos}"
        
        SAÍDA OBRIGATÓRIA (Apenas uma linha):
        Órgão/Tecido - Alvo (Receptor/Enzima) - Ação (Agonista/Inibidor/Expressão)
        
        Regras:
        - Se não houver info, use "N/A".
        - Responda em Português.
        """
        resposta_ia = simple_gemini_text(prompt, api_key)
        
        if resposta_ia:
            clean = resposta_ia.replace("\n", " ").strip().replace("Output:", "").replace("*", "")
            return clean

    # 2. PLANO B (Fallback Manual)
    palavras = titulo.split()
    palavras_chave = [p for p in palavras if len(p) > 5 and p.upper() not in BLACKLIST_MONSTRO]
    resumo_fallback = " | ".join(palavras_chave[:3])
    
    if not resumo_fallback: return "N/A - N/A - N/A"
        
    return f"Foco Provável: {resumo_fallback} (Sem IA)"

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
