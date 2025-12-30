"""
Lemos Lambda Backend v2.0
Hybrid Deterministic-LLM Pipeline for Pharmacological Discovery
SDK Oficial Gemini 2.5 (Pro + Flash) — API google-genai
"""

from Bio import Entrez, Medline
from google import genai
from google.genai import types
import json
import re
from collections import Counter
from typing import List, Dict

# ================= CONFIGURAÇÃO =================

Entrez.email = "pesquisador_guest@unifesp.br"

MODELO_PRO = "gemini-2.5-pro"
MODELO_FLASH = "gemini-2.5-flash"

# ================= BLACKLIST CATEGORIZADA =================

BLACKLIST_TOTAL = {
    "metodologia": {"STUDY", "ANALYSIS", "REVIEW", "DATA", "RESULTS", "CONCLUSION", "METHODS", "TRIAL", "RCT", "COHORT"},
    "estatistica": {"PVALUE", "ANOVA", "RATIO", "ODDS", "STATISTICS", "SIGNIFICANT", "DIFFERENCE", "BASELINE", "SCORE", "KAPLAN"},
    "clinico": {"SURGERY", "RESECTION", "DIAGNOSIS", "PROGNOSIS", "MANAGEMENT", "THERAPY", "TREATMENT", "SYNDROME", "PATIENT", "HOSPITAL", "BIOPSY"},
    "anatomia_ruido": {"KIDNEY", "PROSTATE", "LIVER", "LUNG", "HEART", "BLOOD", "URINE", "CELLS", "PATIENTS"},
    "geral": {"ROLE", "EFFECT", "IMPACT", "POTENTIAL", "NOVEL", "ASSOCIATION", "EVALUATION", "IDENTIFICATION", "ACTIVATION"}
}
UNIFIED_BLACKLIST = set().union(*BLACKLIST_TOTAL.values())

# ================= REGEX FARMACÊUTICO =================

REGEX_PATTERNS = [
    r"\b[A-Z]{2,5}\d{1,4}[A-Z]?\b",
    r"\b[A-Z]{3,6}R\b",
    r"\b[A-Z]{2,4}[- ]?\d{3,6}\b",
    r"\b[A-Z][a-z]{3,}(?:ine|mab|ib|ol|on|one|ide|ate|ase|an)\b"
]

# ================= GEMINI CLIENT =================

def configurar_gemini(api_key: str):
    """Inicializa o cliente genai sem depender de st.session_state."""
    if api_key:
        try:
            return genai.Client(api_key=api_key.strip())
        except:
            return None
    return None

def gerar_com_gemini(prompt: str, client, is_json: bool = False) -> str:
    """Motor de geração com roteamento inteligente Pro/Flash."""
    if not client:
        return ""

    generation_config = types.GenerateContentConfig(
        temperature=0.05 if is_json else 0.3,
        max_output_tokens=1024,
        safety_settings=[
            types.SafetySetting(
                category="HARM_CATEGORY_DANGEROUS_CONTENT",
                threshold="BLOCK_NONE"
            )
        ]
    )

    # Roteamento: JSON prefere Flash, Texto prefere Pro
    modelos = [MODELO_FLASH, MODELO_PRO] if is_json else [MODELO_PRO, MODELO_FLASH]

    for modelo in modelos:
        try:
            response = client.models.generate_content(
                model=modelo,
                contents=prompt,
                config=generation_config
            )
            if response and response.text:
                return response.text.strip()
        except:
            continue
    return ""

# ================= VALIDAÇÃO COM IA =================

def validar_com_ia(candidatos: List[str], contexto: str, api_key: str) -> List[str]:
    """Valida termos usando o cliente configurado localmente."""
    client = configurar_gemini(api_key)
    if not candidatos or not client:
        return []

    prompt = f"""
    Expert Pharmacologist Review. Context: {contexto}
    Return ONLY a pure JSON list of valid Drugs or Molecular Targets (Receptors, Channels, Enzymes).
    LIST: {", ".join(candidatos[:120])}
    """

    resposta = gerar_com_gemini(prompt, client, is_json=True)

    try:
        clean = re.sub(r"```json|```", "", resposta).strip()
        return json.loads(clean)
    except:
        return []

# ================= PIPELINE PUBMED =================

def minerar_pubmed(termo_base: str, email: str, api_key: str, usar_ia: bool = True) -> Dict:
    """Pipeline puro: recebe todos os parâmetros necessários do frontend."""
    Entrez.email = email
    query = f"({termo_base}) AND (drug OR inhibitor OR compound OR target) AND (2020:2026[Date])"

    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=200)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            return {"termos_indicados": [], "metadata": {"total_articles": 0}}

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)

        full_text = ""
        total_docs = 0
        for r in records:
            full_text += f" {r.get('TI','')} {r.get('AB','')} {' '.join(r.get('OT', []))}"
            total_docs += 1

        raw_found = []
        for pattern in REGEX_PATTERNS:
            raw_found.extend(re.findall(pattern, full_text))

        normalized = [x.upper().replace("-", "").replace(" ", "") for x in raw_found]
        cleaned = [x for x in normalized if x not in UNIFIED_BLACKLIST and len(x) > 2]

        counts = Counter(cleaned)
        candidatos = [k for k, _ in counts.most_common(120)]

        if usar_ia and api_key:
            validados = validar_com_ia(candidatos, termo_base, api_key)
            final_list = validados if validados else candidatos[:50]
        else:
            final_list = candidatos[:50]

        return {
            "termos_indicados": final_list[:50],
            "metadata": {"total_articles": total_docs}
        }

    except Exception as e:
        return {"termos_indicados": [], "error": str(e)}

# ================= ABSTRACT ANALYSIS =================

def analisar_abstract_com_ia(titulo: str, dados_curtos: str, api_key: str, lang: str = "pt") -> str:
    """Análise PhD injetando o idioma via parâmetro."""
    client = configurar_gemini(api_key)
    if not client:
        return "API Key pendente."

    prompt = f"""
    You are a PhD Pharmacologist. Deduce the primary pharmacological mechanism.
    Return EXACTLY one line in {lang} in this format:
    Organ/Tissue - Molecular Target - Pharmacological Action

    TITLE: {titulo}
    ABSTRACT: {dados_curtos}
    """

    resposta = gerar_com_gemini(prompt, client, is_json=False)
    if resposta:
        return resposta.split("\n")[0].replace("*", "").strip()
    return "Análise indisponível."

# ================= WRAPPERS =================

def buscar_alvos_emergentes_pubmed(alvo: str, email: str, api_key: str, usar_ia: bool = True) -> List[str]:
    res = minerar_pubmed(alvo, email, api_key, usar_ia)
    return res.get("termos_indicados", [])

def consultar_pubmed_count(termo: str, contexto: str, email: str, ano_ini: int, ano_fim: int) -> int:
    Entrez.email = email
    query = f"({termo}) AND ({contexto}) AND ({ano_ini}:{ano_fim}[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle)
        handle.close()
        return int(record["Count"])
    except:
        return 0

def buscar_resumos_detalhados(termo: str, orgao: str, email: str, ano_ini: int, ano_fim: int) -> List[Dict]:
    Entrez.email = email
    query = f"({termo}) AND ({orgao}) AND ({ano_ini}:{ano_fim}[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=6)
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        return [{"Title": r.get("TI", ""), "Info_IA": r.get("AB", "")[:1000], "Link": f"https://pubmed.ncbi.nlm.nih.gov/{r.get('PMID', '')}/"} for r in records]
    except:
        return []

def buscar_todas_noticias(lang_code: str) -> List[Dict]:
    return []
