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

# LISTA NEGRA: Bloqueia metodologias, anatomia genérica e ruído estatístico
BLACKLIST_MONSTRO = {
    # 1. Termos de Imagem e Física (Lixo Técnico)
    "MRI", "CT", "PET", "RBE", "IMPT", "RBE", "VIII", "RATIO",
    
    # 2. Termos Metodológicos Genéricos (Lixo de Procedimento)
    "ASSOCIATION", "EVALUATION", "UNCOMMON", "DISEASE", "BACKGROUND", 
    "OBJECTIVE", "METHODS", "RESULTS", "CONCLUSION", "ABSTRACT", 
    "INTRODUCTION", "STUDY", "ANALYSIS", "DATA", "STATISTICS", 
    "SIGNIFICANT", "DIFFERENCE", "BETWEEN", "AMONG", "WITHIN", 
    "DURING", "PREVALENCE", "INCIDENCE", "RISK", "FACTOR", "ROLE", 
    "POTENTIAL", "NOVEL", "DIAGNOSTIC", "ARTIFICIAL", "MANAGEMENT",
    "PROGNOSTIC", "FACTORS", "NEUROMODULATION", "IMPLICATIONS",
    "CLINICAL", "REVIEW", "META-ANALYSIS", "SYSTEMATIC", "SURVEY",
    
    # 3. Contexto Clínico/Anatomia (Lixo de Contexto)
    "INTRAUTERINE", "GERMLINE", "BLADDER", "CANCER", "URINARY", 
    "UROTHELIAL", "MUSCLE", "OVERACTIVE", "TUMOR", "CARCINOMA", 
    "PELVIC", "URETHRAL", "UROLOGIC", "NEUROGENIC", "DIABETES", 
    "SJOGREN", "INTRAVESICAL", "VOID", "VOIDING", "DETRUSOR", 
    "PATIENT", "PATIENTS", "CHILDREN", "ADULT", "NHANES", "WOMEN", "MEN",
    
    # 4. Biologia Genérica e Imunologia Ampla (Lixo Biológico)
    "DNA", "RNA", "ATP", "GENE", "PROTEIN", "CELL", "EXPRESSION",
    "BCG", "GUERIN", "CALMETTE", "BACILLUS", "HLA", "CAR",
    
    # 5. Stop Words e Conectivos
    "THE", "AND", "FOR", "NOT", "BUT", "VIA", "ALL", "WITH", "FROM", "AFTER"
}

MAPA_SINONIMOS_BASE = {
    "BLADDER": "(Bladder OR Urothelial OR Urothelium)",
    "PAIN": "(Pain OR Nociception OR Analgesia)",
    "INFLAMMATION": "(Inflammation OR Cytokines OR NF-kappaB)"
}

# ================= GEMINI CORE (COM TRATAMENTO DE ERRO BLINDADO) =================

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
    """
    Tenta obter texto com insistência (Retry Logic). 
    Se der erro de limite (429), espera e tenta de novo até 3 vezes por modelo.
    """
    if not api_key: return None
    
    headers = {"Content-Type": "application/json"}
    payload = {
        "contents": [{"parts": [{"text": prompt}]}], 
        "generationConfig": {"temperature": 0.3}
    }
    
    for modelo_raw in MODELOS_ATIVOS:
        modelo = clean_model_name(modelo_raw)
        
        # Tenta até 3 vezes por modelo
        for tentativa in range(3):
            try:
                url = f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={api_key.strip()}"
                resp = requests.post(url, headers=headers, json=payload, timeout=15)
                
                if resp.status_code == 200:
                    return resp.json()["candidates"][0]["content"]["parts"][0]["text"]
                
                elif resp.status_code == 429:
                    # Cota excedida? Espera um pouco e tenta de novo
                    time.sleep(2 + (tentativa * 1)) 
                    continue 
                
                else:
                    # Erro fatal no modelo, pula para o próximo
                    break 
                    
            except Exception: 
                time.sleep(1)
                continue
            
    return None

# ================= MINERAÇÃO ESTRUTURADA =================

def ner_extraction_batch(titulos_keywords: List[str], api_key: str) -> List[str]:
    """Extração focada em Alvos Moleculares e Drogas."""
    if not titulos_keywords: return []
    
    texto_input = "\n".join([f"- {t}" for t in titulos_keywords[:60]])
    
    prompt = f"""
    Role: Molecular Pharmacologist.
    Task: Extract strictly MOLECULAR TARGETS (Receptors, Enzymes, Genes, Ion Channels) and DRUGS/COMPOUNDS.
    
    STRICT EXCLUSIONS (Do NOT extract):
    - Methods: "Association", "Evaluation", "MRI", "Analysis", "Study".
    - Context: "Disease", "Intrauterine", "Germline", "Neuromodulation".
    - General: "Patient", "Bladder", "Cancer", "Cell", "Protein".
    
    EXTRACT ONLY SPECIFIC ENTITIES:
    - Targets: "HSP90", "P2X3", "SGLT2", "LRP5", "mTOR", "VEGF".
    - Drugs: "Alantolactone", "Mirabegron", "Cisplatin", "Botox".
    
    INPUT:
    {texto_input}
    
    OUTPUT: JSON list of strings (entity names only).
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
        
        # 1. TENTATIVA VIA IA (Principal)
        if api_key and usar_ia:
            entidades = ner_extraction_batch(raw_texts, api_key)
        
        # 2. TENTATIVA VIA REGEX (Complementar/Fallback)
        if True:
            texto_full = " ".join(raw_texts)
            
            # A: Códigos (Letras + Números) -> ex: P2X3, HSP90
            regex_codigos = r'\b[A-Z]{2,}[0-9]+[A-Z0-9-]*\b'
            candidatos_codigos = re.findall(regex_codigos, texto_full)
            
            # B: Fármacos (PascalCase + Sufixos Químicos)
            sufixos = r'(?:ine|in|mab|ib|ol|on|one|il|ide|ate|ase)\b'
            regex_farmacos = r'\b[A-Z][a-z]{3,}' + sufixos
            candidatos_farmacos = re.findall(regex_farmacos, texto_full)
            
            # C: Acrônimos de Vias (3+ Maiúsculas) -> ex: VEGF, mTOR
            regex_acronimos = r'\b[A-Z]{3,}\b'
            candidatos_acronimos = re.findall(regex_acronimos, texto_full)

            entidades.extend(candidatos_codigos)
            entidades.extend(candidatos_farmacos)
            entidades.extend(candidatos_acronimos)

        # 3. FILTRAGEM FINAL
        entidades_limpas = []
        for e in entidades:
            e = e.strip(".,-;:()[] ")
            if len(e) < 3: continue 
            if e.isdigit(): continue
            
            # Bloqueio explícito da Blacklist
            if e.upper() in BLACKLIST_MONSTRO: continue
            
            # Bloqueio extra para preposições
            if e.lower() in ["with", "from", "after", "during"]: continue
            
            entidades_limpas.append(e)

        # 4. Contagem e Seleção
        counts = Counter(entidades_limpas)
        
        # Limiar adaptativo
        limit = 2 if len(counts) > 20 else 1
        recorrentes = [e for e, c in counts.items() if c >= limit]
        
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

# ================= WRAPPERS E ANÁLISE =================

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
    """Gera a Análise Lemos Lambda. Se a IA falhar, gera um resumo tático (Plano B)."""
    
    # 1. Tenta via IA (Gemini) com Prompt OTIMIZADO (Menos tokens = Mais rápido)
    if api_key:
        prompt = f"""
        Paper: "{titulo}"
        Task: 1 short sentence summarizing mechanism.
        Format: [TARGET] affects [PATHWAY] to improve [CONDITION].
        """
        resposta_ia = simple_gemini_text(prompt, api_key)
        
        if resposta_ia:
            return resposta_ia.replace("\n", " ").strip()

    # 2. PLANO B (Fallback): Extração manual se a IA falhar
    palavras = titulo.split()
    palavras_chave = [p for p in palavras if len(p) > 5 and p.upper() not in BLACKLIST_MONSTRO]
    
    resumo_fallback = " | ".join(palavras_chave[:3])
    
    if not resumo_fallback:
        return "Análise pendente (Verificar artigo original)."
        
    return f"🧬 Foco provável: {resumo_fallback} (Gerado automaticamente sem IA)"

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
