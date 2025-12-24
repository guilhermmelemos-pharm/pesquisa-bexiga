import streamlit as st
from Bio import Entrez
import google.generativeai as genai
from tenacity import retry, stop_after_attempt, wait_exponential
import re
import time
import ast
from collections import Counter

# ======================================================
# CONFIGURAÇÃO GERAL
# ======================================================
Entrez.email = "pesquisador_guest@unifesp.br"

# ======================================================
# MODELOS GEMINI (ORDEM DE FALLBACK)
# ======================================================
MODELOS_GEMINI = [
    "models/gemini-2.5-flash",
    "models/gemini-2.5-pro",
    "models/gemini-2.0-flash-exp",
    "models/gemini-flash-latest",
]

# ======================================================
# MAPA SEMÂNTICO DE ÓRGÃOS
# ======================================================
MAPA_SINONIMOS = {
    "BLADDER": "Bladder OR Urothelium OR Detrusor OR Vesical OR Urethra OR Micturition",
    "KIDNEY": "Kidney OR Renal OR Nephron OR Glomerulus",
    "HEART": "Heart OR Cardiac OR Myocardium",
    "BRAIN": "Brain OR CNS OR Neuron OR Glia",
    "LIVER": "Liver OR Hepatic",
    "LUNG": "Lung OR Pulmonary"
}

# ======================================================
# UTIL — CONFIGURA GEMINI
# ======================================================
def configurar_gemini():
    api_key = st.session_state.get("api_key_usuario", "").strip()
    if not api_key:
        return False
    genai.configure(api_key=api_key)
    return True

# ======================================================
# 🔬 EXTRAÇÃO DE SINAIS MOLECULARES (ANTI-LIMITE)
# ======================================================
def extrair_sinais_moleculares(titulo, keywords="", abstract=""):
    sinais = []

    if keywords:
        sinais.append(f"Keywords: {keywords[:250]}")

    if abstract:
        frases = re.split(r'(?<=[.!?])\s+', abstract)
        if frases:
            sinais.append(f"Intro: {frases[0][:300]}")
            if len(frases) > 1:
                sinais.append(f"Conclusion: {frases[-1][:300]}")

    return f"{titulo}. " + " ".join(sinais)

# ======================================================
# 🧹 CURADORIA MOLECULAR COM GEMINI
# ======================================================
def _faxina_ia(lista_suja):
    if not configurar_gemini():
        return lista_suja[:60]

    lista_str = ", ".join(lista_suja)

    prompt = f"""
You are a Senior PhD in Experimental Pharmacology.

TASK:
Curate a list of STRICTLY MOLECULAR entities.

KEEP ONLY:
Genes, receptors, ion channels, enzymes, transporters,
ligands, agonists, antagonists, inhibitors, metabolites,
intracellular signaling pathways.

REMOVE COMPLETELY:
Diseases (BPH, OAB, LUTS),
clinical terms,
anatomy,
symptoms,
patients,
outcomes,
reviews,
methods.

INPUT LIST:
{lista_str}

OUTPUT RULES:
- Return ONLY a valid Python list
- Strings only
- Max 60 items
"""

    for modelo in MODELOS_GEMINI:
        try:
            model = genai.GenerativeModel(modelo)
            resp = model.generate_content(prompt, generation_config={"temperature": 0.1})
            texto = resp.text.strip().replace("```python", "").replace("```", "")
            if texto.startswith("["):
                return ast.literal_eval(texto)
        except:
            time.sleep(1)
            continue

    return lista_suja[:60]

# ======================================================
# 🤖 ANÁLISE MOLECULAR (TÍTULO + SINAIS)
# ======================================================
def analisar_abstract_com_ia(titulo, keywords, abstract, lang="pt"):
    if not configurar_gemini():
        return "⚠️ IA desativada"

    sinais = extrair_sinais_moleculares(titulo, keywords, abstract)
    idioma = "Português" if lang == "pt" else "English"

    prompt = f"""
You are a Senior PhD in Molecular Pharmacology.

TASK:
Infer molecular mechanism ONLY from the signals below.

ALLOWED:
Targets, genes, receptors, ion channels,
enzymes, pathways, drugs, agonists, antagonists.

FORBIDDEN:
Diseases, clinical terms, anatomy, therapy, patients.

MOLECULAR SIGNALS:
{sinais}

FORMAT:
Target → Molecule → Mechanism

RULES:
- One single line
- Max 20 words
- No clinical language
- Language: {idioma}
"""

    for modelo in MODELOS_GEMINI:
        try:
            model = genai.GenerativeModel(modelo)
            resp = model.generate_content(prompt, generation_config={"temperature": 0.3})
            if resp.text:
                return resp.text.strip()
        except:
            time.sleep(1)
            continue

    return "⚠️ IA indisponível"

# ======================================================
# 🔎 PUBMED COUNT
# ======================================================
@retry(stop=stop_after_attempt(3), wait=wait_exponential())
def _fetch_pubmed_count(query):
    h = Entrez.esearch(db="pubmed", term=query, retmax=0)
    r = Entrez.read(h)
    h.close()
    return int(r["Count"])

def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    if email:
        Entrez.email = email
    ctx = MAPA_SINONIMOS.get(contexto.upper(), contexto) if contexto else ""
    q = f"({termo})"
    if ctx:
        q += f" AND ({ctx})"
    q += f" AND ({ano_ini}:{ano_fim}[Date - Publication]) AND NOT Review[pt]"
    try:
        return _fetch_pubmed_count(q)
    except:
        return 0

# ======================================================
# 📄 BUSCA RESUMOS (ENXUTO)
# ======================================================
def buscar_resumos_detalhados(termo, orgao, email, ano_ini, ano_fim):
    if email:
        Entrez.email = email

    org = MAPA_SINONIMOS.get(orgao.upper(), orgao)
    query = f"({termo}) AND ({org}) AND ({ano_ini}:{ano_fim}[Date - Publication])"

    try:
        h = Entrez.esearch(db="pubmed", term=query, retmax=5, sort="relevance")
        r = Entrez.read(h)
        h.close()

        if not r["IdList"]:
            return []

        h = Entrez.efetch(db="pubmed", id=r["IdList"], rettype="medline", retmode="text")
        raw = h.read()
        h.close()

        artigos = []
        for bloco in raw.split("\n\nPMID-"):
            tit, abs_, kw, pmid = "", "", "", ""
            for l in bloco.split("\n"):
                if l.startswith("TI  - "): tit = l[6:].strip()
                if l.startswith("AB  - "): abs_ = l[6:1200].strip()
                if l.startswith(("OT  - ", "KW  - ")): kw += l[6:].strip() + ", "
                if l.startswith("PMID- "): pmid = l[6:].strip()
            if tit:
                artigos.append({
                    "Title": tit,
                    "Keywords": kw,
                    "Abstract": abs_,
                    "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                })
        return artigos
    except:
        return []

# ======================================================
# 🧪 MINERAÇÃO DE ALVOS
# ======================================================
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email:
        Entrez.email = email

    termo = MAPA_SINONIMOS.get(termo_base.upper(), termo_base)
    query = f"({termo}) AND (2018:2030[Date - Publication]) AND NOT Review[pt]"

    try:
        h = Entrez.esearch(db="pubmed", term=query, retmax=2000)
        r = Entrez.read(h)
        h.close()

        if not r["IdList"]:
            return []

        h = Entrez.efetch(db="pubmed", id=r["IdList"], rettype="medline", retmode="text")
        raw = h.read().split("\n\nPMID-")
        h.close()

        blacklist = {
            "OAB","BPH","LUTS","PATIENT","CLINICAL","STUDY","REVIEW",
            "TREATMENT","DISEASE","SYMPTOM","MANAGEMENT","OUTCOME"
        }

        candidatos = []
        for art in raw:
            texto = ""
            for l in art.split("\n"):
                if l.startswith(("TI  - ", "OT  - ")):
                    texto += l[6:] + " "
            achados = re.findall(r'\b[A-Z][A-Za-z0-9-]{2,}\b', texto)
            for t in achados:
                if t.upper() not in blacklist and len(t) >= 4:
                    candidatos.append(t)

        top = [x for x, _ in Counter(candidatos).most_common(200)]

        if usar_ia and st.session_state.get("api_key_usuario"):
            return _faxina_ia(top)

        return top[:60]

    except:
        return []
