"""
Lemos Lambda Backend v2.5 - Edição Platinum
Foco: Farmacologia Estrita | Órgão Específico | Zero Ruído Acadêmico
"""

from Bio import Entrez, Medline
from google import genai
from google.genai import types
import json, re
from collections import Counter
from typing import List, Dict

# ================= CONFIGURAÇÃO =================
Entrez.email = "pesquisador_guest@unifesp.br"
MODELO_PRO = "gemini-2.5-pro"
MODELO_FLASH = "gemini-2.5-flash"

# ================= SUPER BLACKLIST (A "VAXINA" DO SISTEMA) =================

BLACKLIST_TOTAL = {
    "metodologia": {"STUDY", "ANALYSIS", "REVIEW", "DATA", "RESULTS", "CONCLUSION", "METHODS", "TRIAL", "RCT", "COHORT", "VALIDATION", "DETECTION", "INVESTIGATION", "EVALUATE", "PROVIDE", "ADDITION"},
    "estatistica": {"PVALUE", "ANOVA", "RATIO", "ODDS", "STATISTICS", "SIGNIFICANT", "DIFFERENCE", "BASELINE", "SCORE", "KAPLAN", "MEDIAN", "REGRESSION", "CORRELATION", "FRACTION"},
    "clinico": {"SURGERY", "RESECTION", "DIAGNOSIS", "PROGNOSIS", "MANAGEMENT", "THERAPY", "TREATMENT", "SYNDROME", "PATIENT", "HOSPITAL", "BIOPSY", "POPULATION", "INFECTION", "RETENTION", "DYSFUNCTION", "CONTROL", "PRECISION"},
    "processos_bio": {"EXPRESSION", "PROGRESSION", "FUNCTION", "INVASION", "REGULATION", "OVEREXPRESSION", "INFLAMMATION", "METHYLATION", "FORMATION", "TRANSCRIPTION", "MUTATION", "INCREASE", "INTERACTION", "SENSITIVITY"},
    "anatomia_ruido": {"KIDNEY", "PROSTATE", "LIVER", "LUNG", "HEART", "BLOOD", "URINE", "CELLS", "PATIENTS", "ORGAN", "TISSUE", "HUMAN", "MICE", "RAT"},
    "geral": {"ROLE", "EFFECT", "IMPACT", "POTENTIAL", "NOVEL", "ASSOCIATION", "EVALUATION", "IDENTIFICATION", "ACTIVATION", "DISEASE", "POUR", "VOLUME", "COMMON", "RADIATION", "MEDICINE"}
}
UNIFIED_BLACKLIST = set().union(*BLACKLIST_TOTAL.values())

# ================= REGEX RESTRITIVO =================

REGEX_PATTERNS = [
    r"\b[A-Z]{2,5}\d{1,4}[A-Z]?\b", # Siglas Alvo (TRPV1, P2X3)
    r"\b[A-Z]{2,4}[- ]?\d{3,6}\b",  # Códigos (GYY-4137)
    r"\b[A-Za-z]{4,}(?:ine|mab|ib|ol|on|one|ide|ate|ase|an|tin|pril|afil|arin)\b" # Sufixos Farmacológicos
]

# ================= GEMINI CLIENT =================

def configurar_gemini(api_key: str):
    if api_key:
        try: return genai.Client(api_key=api_key.strip())
        except: return None
    return None

def gerar_com_gemini(prompt: str, client, is_json: bool = False) -> str:
    if not client: return ""
    config = types.GenerateContentConfig(
        temperature=0.01 if is_json else 0.2, # Quase determinístico para evitar alucinação
        max_output_tokens=1024,
        safety_settings=[types.SafetySetting(category="HARM_CATEGORY_DANGEROUS_CONTENT", threshold="BLOCK_NONE")]
    )
    modelos = [MODELO_FLASH, MODELO_PRO] if is_json else [MODELO_PRO, MODELO_FLASH]
    for modelo in modelos:
        try:
            response = client.models.generate_content(model=modelo, contents=prompt, config=config)
            if response and response.text: return response.text.strip()
        except: continue
    return ""

# ================= VALIDAÇÃO SEMÂNTICA (PROMPT DE EXCLUSÃO) =================

def validar_com_ia(candidatos: List[str], contexto: str, api_key: str) -> List[str]:
    client = configurar_gemini(api_key)
    if not candidatos or not client: return []

    prompt = f"""
    ROLE: Expert PhD Pharmacologist & Molecular Biologist.
    CONTEXT: {contexto} (Lower Urinary Tract).

    TASK:
    From the list below, extract ONLY specific Drugs or druggable Molecular Targets.

    STRICT EXCLUSIONS (CRITICAL):
    - NO generic biological processes (ex: progression, expression, invasion).
    - NO anatomical parts or other organs (ex: prostate, kidney, organ).
    - NO clinical or statistical terms (ex: median, regression, population).
    - NO verbs or investigative actions (ex: provide, evaluate, interaction).

    LIST: {", ".join(candidatos[:120])}

    OUTPUT: Pure JSON list of strings only.
    """
    resposta = gerar_com_gemini(prompt, client, is_json=True)
    try:
        clean = re.sub(r"```json|```", "", resposta).strip()
        return json.loads(clean)
    except: return []

# ================= PIPELINE PUBMED =================

def minerar_pubmed(termo_base: str, email: str, api_key: str, usar_ia: bool = True) -> Dict:
    Entrez.email = email
    # Filtro de exclusão direto na fonte (PubMed NOT)
    query = f"({termo_base}) NOT (prostate[TI] OR kidney[TI] OR cancer progression[TI]) AND (drug OR inhibitor OR agonist OR target) AND (2020:2026[Date])"

    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=200)
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return {"termos_indicados": []}

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)

        raw_found = []
        for r in records:
            text = f" {r.get('TI','')} {r.get('AB','')} {' '.join(r.get('OT', []))}"
            for pattern in REGEX_PATTERNS:
                raw_found.extend(re.findall(pattern, text))

        # Normalização e Limpeza contra a Super Blacklist
        cleaned = []
        for x in raw_found:
            x_norm = x.strip().upper().replace("-", "").replace(" ", "")
            if x_norm not in UNIFIED_BLACKLIST and len(x_norm) > 2:
                # Mantém a capitalização original se não for sigla (para nomes de fármacos)
                cleaned.append(x if not x.isupper() else x_norm)

        counts = Counter(cleaned)
        candidatos = [k for k, _ in counts.most_common(120)]

        if usar_ia and api_key:
            validados = validar_com_ia(candidatos, termo_base, api_key)
            return {"termos_indicados": validados[:50] if validados else candidatos[:50]}
        
        return {"termos_indicados": candidatos[:50]}

    except Exception as e:
        return {"termos_indicados": [], "error": str(e)}

# ================= ABSTRACT ANALYSIS (PhD DEDUCTION) =================

def analisar_abstract_com_ia(titulo: str, dados_curtos: str, api_key: str, lang: str = "pt") -> str:
    client = configurar_gemini(api_key)
    if not client: return "API Key pendente."
    
    prompt = f"""
    PhD Pharmacologist: Deduce mechanism. One line only in {lang}.
    Format: Organ - Target - Action
    Title: {titulo}
    Abstract: {dados_curtos}
    """
    res = gerar_com_gemini(prompt, client)
    return res.split("\n")[0].replace("*", "").strip() if res else "Indisponível"

# Wrappers de contagem e resumos permanecem iguais...
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    Entrez.email = email
    query = f"({termo}) AND ({contexto}) NOT (prostate[TI] OR kidney[TI]) AND ({ano_ini}:{ano_fim}[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle); handle.close()
        return int(record["Count"])
    except: return 0

def buscar_resumos_detalhados(termo, orgao, email, y_ini, y_fim):
    Entrez.email = email
    query = f"({termo}) AND ({orgao}) AND ({y_ini}:{y_fim}[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=6)
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        return [{"Title": r.get("TI", ""), "Info_IA": r.get("AB", "")[:1000], "Link": f"https://pubmed.ncbi.nlm.nih.gov/{r.get('PMID', '')}/"} for r in Medline.parse(handle)]
    except: return []

def buscar_todas_noticias(lang_code): return []
