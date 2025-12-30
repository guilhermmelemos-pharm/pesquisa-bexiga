"""
Lemos Lambda Backend v2.5
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
    "geral": {"ROLE", "EFFECT", "IMPACT", "POTENTIAL", "NOVEL", "ASSOCIATION", "EVALUATION", "IDENTIFICATION", "ACTIVATION", "DISEASE", "POUR", "VOLUME"}
}
UNIFIED_BLACKLIST = set().union(*BLACKLIST_TOTAL.values())

# ================= REGEX FARMACÊUTICO REFINADO =================

REGEX_PATTERNS = [
    # 1. Siglas e Alvos com Números (ex: TRPV1, P2X3)
    r"\b[A-Z]{2,5}\d{1,4}[A-Z]?\b",
    # 2. Códigos de Compostos Experimentais (ex: GYY-4137, BAY-123)
    r"\b[A-Z]{2,4}[- ]?\d{3,6}\b",
    # 3. NOVO: Nomes de Fármacos por Sufixos IUPAC/WHO (ex: Sildenafil, Mirabegron)
    r"\b[A-Za-z]{3,}(?:ine|mab|ib|ol|on|one|ide|ate|ase|an|tin|pril|afil|arin|vir|setron|statin)\b"
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

    # Roteamento: JSON (Validação) prefere Flash, Texto (Dedução) prefere Pro
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

# ================= VALIDAÇÃO COM IA (FOCO ANATÔMICO E QUÍMICO) =================

def validar_com_ia(candidatos: List[str], contexto: str, api_key: str) -> List[str]:
    """Valida termos focando em fármacos e eliminando sinônimos de outros órgãos."""
    client = configurar_gemini(api_key)
    if not candidatos or not client:
        return []

    prompt = f"""
    Expert Pharmacologist and Anatomist Review.
    Context: {contexto} (Lower Urinary Tract Focus).

    TASK:
    1. From the list below, return ONLY valid Drugs (names or codes) or Molecular Targets.
    2. STRICTLY REMOVE anatomical terms or synonyms of other organs (ex: PROSTATE, SEMINAL, KIDNEY, NEPHRON, URETER).
    3. REMOVE cell lines (HT1376) and statistical/generic noise.
    4. Group chemical synonyms into the most standard pharmacological name.

    LIST:
    {", ".join(candidatos[:120])}

    Return a pure JSON list of strings. No comments.
    """

    resposta = gerar_com_gemini(prompt, client, is_json=True)

    try:
        clean = re.sub(r"```json|```", "", resposta).strip()
        return json.loads(clean)
    except:
        return []

# ================= PIPELINE PUBMED (QUERY NEGATIVA) =================

def minerar_pubmed(termo_base: str, email: str, api_key: str, usar_ia: bool = True) -> Dict:
    """Busca com filtro de exclusão para evitar 'contaminação' de outros órgãos."""
    Entrez.email = email
    
    # Query v2.5: Foca no alvo mas exclui explicitamente Próstata e Rim do título para evitar falsos sinônimos
    query = f"({termo_base}) NOT (prostate[TI] OR kidney[TI] OR renal[TI]) AND (drug OR inhibitor OR compound OR target) AND (2020:2026[Date])"

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

        # Normalização inteligente (apenas para siglas, nomes de fármacos mantêm capitalização)
        cleaned = []
        for x in raw_found:
            x_clean = x.strip()
            # Se for sigla (toda maiúscula), normaliza. Se for nome de fármaco, limpa espaços.
            if x_clean.isupper():
                x_norm = x_clean.replace("-", "").replace(" ", "")
            else:
                x_norm = x_clean
            
            if x_norm.upper() not in UNIFIED_BLACKLIST and len(x_norm) > 2:
                cleaned.append(x_norm)

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

# ================= ABSTRACT ANALYSIS (TRADUÇÃO INJETADA) =================

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
    # Filtro para contagem também respeitar a exclusão anatômica
    query = f"({termo}) AND ({contexto}) NOT (prostate[TI] OR kidney[TI]) AND ({ano_ini}:{ano_fim}[Date])"
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
