"""
Lemos Lambda Backend v2.0 - Platinum Edition (Refined)
Foco: Farmacologia Estrita | Zero Ruído Acadêmico | Anti-Linguagem Natural
Autor: Guilherme Lemos
"""

from Bio import Entrez, Medline
from collections import Counter
from typing import List, Dict
import json, re

try:
    from google import genai
    from google.genai import types
except Exception:
    genai = None; types = None

# ================= CONFIGURAÇÃO GLOBAL =================
Entrez.email = "pesquisador_guest@unifesp.br"
MODELO_PRO = "gemini-2.5-pro"
MODELO_FLASH = "gemini-2.5-flash"

# ================= SUPER BLACKLIST (BLOQUEIO DE RUÍDO RECENTE) =================
BLACKLIST_TOTAL = {
    "ruido_academico": {
        "CONTROL", "MEDICINE", "ALONGSIDE", "MULTIVARIATE", "DETERMINE", 
        "INDICATE", "ADEQUATE", "DETERMINATION", "EVALUATE", "FURTHERMORE",
        "STRICTLY", "RELEVANT", "CONSIDER", "PRESENT", "ASSOCIATED", "INDICATED"
    },
    "substancias_genericas": {
        "ETHANOL", "CHOLESTEROL", "GLUTAMINE", "TYROSINE", "SULFATE", 
        "GLUCOSE", "OXYGEN", "SODIUM", "POTASSIUM", "ACID", "ALCOHOL"
    },
    "estatistica_e_metodo": {
        "PVALUE", "ANOVA", "RATIO", "ODDS", "STATISTICS", "SIGNIFICANT", 
        "MEDIAN", "REGRESSION", "CORRELATION", "FRACTION", "EQD2", "GSE13507",
        "DATABASE", "POPULATION", "PROTOCOL", "VALIDATE", "IDENTIFY"
    },
    "anatomia_e_geral": {
        "KIDNEY", "PROSTATE", "LIVER", "LUNG", "BLOOD", "URINE", "CELLS",
        "PATIENTS", "ORGAN", "TISSUE", "HUMAN", "MICE", "RAT", "DISEASE",
        "CANCER", "POUR", "VOLUME", "COMMON", "GERMLINE", "GENOMIC"
    }
}
UNIFIED_BLACKLIST = set().union(*BLACKLIST_TOTAL.values())

# ================= REGEX FARMACOLÓGICO (FILTRO DE PRECISÃO) =================
REGEX_PATTERNS = [
    r"\b[A-Z]{2,5}\d{1,4}[A-Z]?\b", # Siglas Alvo: STAT1, TRPV4, CDK1
    r"\b[A-Z]{2,4}[- ]?\d{3,6}\b",  # Códigos: CG0070
    # Exigência de 6 letras para evitar verbos curtos como "Indicate" ou "Adequate"
    r"\b[A-Za-z]{6,}(?:mab|ib|ol|one|ide|ate|ase|pril|afil|tin|ine|arin|oside)\b"
]

# ================= GEMINI CLIENT =================
def configurar_gemini(api_key: str):
    if not api_key or not genai: return None
    try: return genai.Client(api_key=api_key.strip())
    except: return None

def gerar_com_gemini(prompt: str, client, is_json: bool = False) -> str:
    if not client: return ""
    config = types.GenerateContentConfig(temperature=0.0, max_output_tokens=1024)
    modelos = [MODELO_FLASH, MODELO_PRO] if is_json else [MODELO_PRO, MODELO_FLASH]
    for modelo in modelos:
        try:
            resp = client.models.generate_content(model=modelo, contents=prompt, config=config)
            if resp.text: return resp.text.strip()
        except: continue
    return ""

# ================= VALIDAÇÃO SEMÂNTICA (FILTRO CIRÚRGICO) =================
def validar_com_ia(candidatos: List[str], contexto: str, api_key: str) -> List[str]:
    client = configurar_gemini(api_key)
    if not candidatos or not client: return []

    prompt = f"""
ROLE: Senior PhD Pharmacologist.
TASK: Return ONLY a JSON list of:
1. Specific Drugs (ex: Enfortumab, Bupivacaine, Pembrolizumab).
2. Molecular Targets (ex: TRPV4, FGFR3, BCL6, STAT1).

STRICTLY EXCLUDE:
- Adverbs/Verbs (Determine, Indicate, Alongside, Adequate, Multivariate).
- Generic substances (Ethanol, Glutamine, Cholesterol).
- Generic medical terms (Medicine, Disease, Population).

LIST TO FILTER:
{", ".join(candidatos[:120])}

OUTPUT: Pure JSON list only.
"""
    res = gerar_com_gemini(prompt, client, is_json=True)
    try:
        clean = re.sub(r"```json|```", "", res).strip()
        return json.loads(clean)
    except: return []

# ================= PIPELINE PUBMED =================
def minerar_pubmed(termo_base: str, email: str, api_key: str, usar_ia: bool = True) -> Dict:
    Entrez.email = email
    query = f"({termo_base}) AND (drug OR inhibitor OR target) NOT (prostate[TI] OR kidney[TI]) AND (2020:2026[Date])"

    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=200)
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return {"termos_indicados": []}

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)

        raw_found = []
        for r in records:
            text = f"{r.get('TI','')} {r.get('AB','')} {' '.join(r.get('OT', []))}"
            for pattern in REGEX_PATTERNS:
                raw_found.extend(re.findall(pattern, text))

        cleaned = []
        for x in raw_found:
            norm = x.upper().replace("-", "").replace(" ", "").strip()
            if norm not in UNIFIED_BLACKLIST and len(norm) > 2:
                cleaned.append(norm)

        counts = Counter(cleaned)
        candidatos = [k for k, _ in counts.most_common(120)]

        if usar_ia and api_key:
            validados = validar_com_ia(candidatos, termo_base, api_key)
            return {"termos_indicados": validados[:50] or candidatos[:50]}

        return {"termos_indicados": candidatos[:50]}
    except: return {"termos_indicados": []}

# Funções auxiliares (Analisar Abstract, Contagem, Resumos) permanecem as mesmas v2.0
def analisar_abstract_com_ia(titulo, abstract, api_key, lang="pt"):
    client = configurar_gemini(api_key); prompt = f"PhD Pharmacologist. One line in {lang}: Organ – Target – Action. Title: {titulo}. Abstract: {abstract}"
    res = gerar_com_gemini(prompt, client); return res.split("\n")[0].strip() if res else "Indisponível"

def consultar_pubmed_count(termo, contexto, email, y_ini, y_fim):
    Entrez.email = email; query = f"({termo}) AND ({contexto}) NOT (prostate[TI] OR kidney[TI]) AND ({y_ini}:{y_fim}[Date])"
    try: handle = Entrez.esearch(db="pubmed", term=query, retmax=0); rec = Entrez.read(handle); handle.close(); return int(rec["Count"])
    except: return 0

def buscar_resumos_detalhados(termo, orgao, email, y_ini, y_fim):
    Entrez.email = email; query = f"({termo}) AND ({orgao}) AND ({y_ini}:{y_fim}[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=6); rec = Entrez.read(handle); handle.close()
        if not rec["IdList"]: return []
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        return [{"Title": r.get("TI", ""), "Abstract": r.get("AB", "")[:1000], "Link": f"https://pubmed.ncbi.nlm.nih.gov/{r.get('PMID','')}/"} for r in Medline.parse(handle)]
    except: return []

def buscar_todas_noticias(l): return []
