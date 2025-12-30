"""
Lemos Lambda Backend v2.0 - Platinum Edition Final
Foco: Farmacologia Estrita | Órgão Específico | Ineditismo Molecular
Autor: Guilherme Lemos
"""

from Bio import Entrez, Medline
from collections import Counter
from typing import List, Dict
import json
import re

try:
    from google import genai
    from google.genai import types
except Exception:
    genai = None
    types = None

# ================= CONFIGURAÇÃO GLOBAL =================
Entrez.email = "pesquisador_guest@unifesp.br"
MODELO_PRO = "gemini-2.5-pro"
MODELO_FLASH = "gemini-2.5-flash"

# ================= SUPER BLACKLIST CONSOLIDADA =================
BLACKLIST_TOTAL = {
    "metodologia": {
        "STUDY", "ANALYSIS", "REVIEW", "DATA", "RESULTS", "CONCLUSION",
        "METHODS", "TRIAL", "RCT", "COHORT", "VALIDATION", "DETECTION",
        "INVESTIGATION", "EVALUATE", "PROVIDE", "ADDITION", "RELEASE",
        "APPROPRIATE", "DATABASE", "GUIDELINE", "MACHINE", "ROUTINE", "PROTOCOL"
    },
    "estatistica": {
        "PVALUE", "ANOVA", "RATIO", "ODDS", "STATISTICS", "SIGNIFICANT",
        "DIFFERENCE", "BASELINE", "SCORE", "KAPLAN", "MEDIAN",
        "REGRESSION", "CORRELATION", "FRACTION", "VALIDATE", "INVESTIGATE",
        "ELUCIDATE", "FACILITATE", "CANDIDATE"
    },
    "conhecidos_ou_ruido": {
        "CISPLATIN", "ZSTK474", "GEMCITABINE", "ESM1", "ATF6", "TRPV3", "HDAC2",
        "CD8", "SGLT2", "HDAC7", "RNF139", "TMEM208", "PI3K", "FLLL31", "TPI1",
        "KRT14", "STAT3", "HT1376", "CHLOROPROCAINE", "FANCD2", "PSMB5", "TSP4",
        "SMK002", "TLR4", "TP53", "FGL1", "NSUN2", "PSMD8", "DHRS2", "CDC20",
        "CCPG1", "KRT17", "PINCH2", "TISLELIZUMAB", "KDM4A", "KIF26B", "ATF4", 
        "CISD2", "NAT10", "PRDX1", "GARCINOL", "FOXN3", "NEAT1", "CETRELIMAB", 
        "PDE5", "NEUROENDOCRINE"
    },
    "anatomia_ruido": {
        "KIDNEY", "PROSTATE", "LIVER", "LUNG", "HEART", "BLOOD", "URINE",
        "CELLS", "PATIENTS", "ORGAN", "TISSUE", "HUMAN", "MICE", "RAT", 
        "UTERINE", "CHROMATIN", "HISTONE"
    },
    "processos_bio": {
        "EXPRESSION", "PROGRESSION", "FUNCTION", "INVASION", "REGULATION",
        "OVEREXPRESSION", "INFLAMMATION", "METHYLATION", "FORMATION",
        "TRANSCRIPTION", "MUTATION", "INCREASE", "INTERACTION",
        "SENSITIVITY", "HORMONE", "SULFATE", "GLUTATHIONE", "PEROXIDASE"
    }
}
UNIFIED_BLACKLIST = set().union(*BLACKLIST_TOTAL.values())

# ================= REGEX FARMACOLÓGICO =================
REGEX_PATTERNS = [
    r"\b[A-Z]{2,5}\d{1,4}[A-Z]?\b", 
    r"\b[A-Z]{2,4}[- ]?\d{3,6}\b",  
    r"\b[A-Za-z]{5,}(?:mab|ib|ol|one|ide|ate|ase|pril|afil|tin|ine|vir|arin)\b"
]

# ================= GEMINI CLIENT =================
def configurar_gemini(api_key: str):
    if not api_key or not genai: return None
    try: return genai.Client(api_key=api_key.strip())
    except Exception: return None

def gerar_com_gemini(prompt: str, client, is_json: bool = False) -> str:
    if not client or not types: return ""
    config = types.GenerateContentConfig(temperature=0.01 if is_json else 0.2, max_output_tokens=1024)
    modelos = [MODELO_FLASH, MODELO_PRO] if is_json else [MODELO_PRO, MODELO_FLASH]
    for modelo in modelos:
        try:
            response = client.models.generate_content(model=modelo, contents=prompt, config=config)
            if response and response.text: return response.text.strip()
        except Exception: continue
    return ""

# ================= VALIDAÇÃO SEMÂNTICA ULTRA-RESTRITIVA =================
def validar_com_ia(candidatos: List[str], contexto: str, api_key: str) -> List[str]:
    client = configurar_gemini(api_key)
    if not candidatos or not client: return []

    prompt = f"""
ROLE: Senior PhD Pharmacologist.
CONTEXT: {contexto} (Lower Urinary Tract Research).

TASK:
Filter the list below. Return ONLY a pure JSON list of valid Drugs or Molecular Targets.

STRICT EXCLUSIONS (MANDATORY):
1. NO Common Drugs: Cisplatin, Gemcitabine, etc.
2. NO Overstudied Targets: PI3K, STAT3, TP53, TLR4, PDE5.
3. NO anatomical terms, biological processes, or statistical filler.
4. If a term is vague or generic, DELETE it.

LIST TO FILTER:
{", ".join(candidatos[:120])}

OUTPUT:
Pure JSON list of strings only.
"""
    resposta = gerar_com_gemini(prompt, client, is_json=True)
    try:
        clean = re.sub(r"```json|```", "", resposta).strip()
        return json.loads(clean)
    except Exception: return []

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
            texto = f"{r.get('TI','')} {r.get('AB','')} {' '.join(r.get('OT', []))}"
            for pattern in REGEX_PATTERNS:
                raw_found.extend(re.findall(pattern, texto))

        # Normalização Consistente (Ajuste recomendado)
        cleaned = []
        for x in raw_found:
            norm = x.upper().replace("-", "").replace(" ", "").strip()
            if norm not in UNIFIED_BLACKLIST and len(norm) > 2:
                # Armazenamos sempre o termo normalizado para evitar duplicidade semântica
                cleaned.append(norm)

        counts = Counter(cleaned)
        candidatos = [k for k, _ in counts.most_common(120)]

        if usar_ia and api_key:
            validados = validar_com_ia(candidatos, termo_base, api_key)
            return {"termos_indicados": validados[:50] or candidatos[:50]}

        return {"termos_indicados": candidatos[:50]}
    except Exception as e:
        return {"termos_indicados": [], "error": str(e)}

# Funções auxiliares mantidas conforme v2.0 Platinum...
def analisar_abstract_com_ia(titulo, abstract, api_key, lang="pt"):
    client = configurar_gemini(api_key)
    if not client: return "API Key ausente."
    prompt = f"PhD Pharmacologist. One line in {lang}: Organ – Target – Action. Title: {titulo}. Abstract: {abstract}"
    res = gerar_com_gemini(prompt, client)
    return res.split("\n")[0].strip() if res else "Indisponível"

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
        return [{"Title": r.get("TI", ""), "Abstract": r.get("AB", "")[:1000], "Link": f"https://pubmed.ncbi.nlm.nih.gov/{r.get('PMID','')}/"} for r in Medline.parse(handle)]
    except: return []

def buscar_todas_noticias(lang_code): return []
