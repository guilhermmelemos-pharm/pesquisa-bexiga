# backend.py
import time
import re
from difflib import SequenceMatcher
from collections import Counter

import streamlit as st
from Bio import Entrez
import google.generativeai as genai
from tenacity import retry, stop_after_attempt, wait_exponential


# =========================
# CONFIGURAÇÕES GERAIS
# =========================
Entrez.email = "pesquisador_guest@unifesp.br"

MAX_ARTIGOS = 60          # limite duro
BATCH_SIZE = 10           # batch de IA
DELAY_ENTRE_BATCH = 1.2   # proteção rate limit


# =========================
# HELPERS
# =========================
def normalizar_termo(t):
    return re.sub(r"\s+", " ", t.strip().lower())


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()


# =========================
# PUBMED
# =========================
@st.cache_data(ttl=86400, show_spinner=False)
def consultar_pubmed_count(termo, contexto, email, ano_ini=1900, ano_fim=2030):
    Entrez.email = email
    query = f"({termo}) AND ({ano_ini}:{ano_fim}[Date - Publication])"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    record = Entrez.read(handle)
    handle.close()
    return int(record.get("Count", 0))


def buscar_pubmed_ids(termo, email, retmax=MAX_ARTIGOS):
    Entrez.email = email
    handle = Entrez.esearch(
        db="pubmed",
        term=termo,
        retmax=retmax,
        sort="relevance"
    )
    record = Entrez.read(handle)
    handle.close()
    return record.get("IdList", [])


def buscar_titulos_keywords(pmids, email):
    Entrez.email = email
    handle = Entrez.efetch(
        db="pubmed",
        id=",".join(pmids),
        rettype="medline",
        retmode="text"
    )
    text = handle.read()
    handle.close()

    artigos = []
    blocos = text.split("\n\n")
    for bloco in blocos:
        titulo = ""
        keywords = ""
        for linha in bloco.split("\n"):
            if linha.startswith("TI  -"):
                titulo = linha.replace("TI  -", "").strip()
            if linha.startswith("OT  -"):
                keywords += " " + linha.replace("OT  -", "").strip()
        if titulo:
            artigos.append({
                "titulo": titulo,
                "keywords": keywords.strip()
            })
    return artigos


# =========================
# IA — CONFIG
# =========================
def configurar_ia(api_key):
    genai.configure(api_key=api_key)
    return genai.GenerativeModel("gemini-1.5-flash")


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=8))
def chamada_ia(model, prompt):
    resp = model.generate_content(prompt)
    return resp.text if resp and resp.text else ""


# =========================
# IA — EXPANSÃO DE QUERY
# =========================
def expandir_termo_com_ia(termo, model):
    prompt = f"""
Dado o termo biomédico: "{termo}"
Liste até 8 sinônimos ou termos relacionados usados em artigos científicos.
Retorne apenas os termos separados por ponto e vírgula.
"""
    out = chamada_ia(model, prompt)
    termos = [t.strip() for t in out.split(";") if len(t.strip()) > 3]
    return termos[:8]


# =========================
# IA — NER EM BATCH (TÍTULO + KEYWORDS)
# =========================
def ner_em_lote(artigos, model):
    texto = ""
    for i, art in enumerate(artigos):
        texto += f"[ID:{i}] {art['titulo']} | {art['keywords']}\n"

    prompt = f"""
Extraia ALVOS BIOLÓGICOS, FÁRMACOS e PROCESSOS/FUNÇÕES dos textos abaixo.
Ignore abstracts. Use apenas o que está explícito.
Formato de saída:
ID:X | ENTIDADE | TIPO

Textos:
{texto}
"""
    out = chamada_ia(model, prompt)

    resultados = []
    for linha in out.splitlines():
        if "|" in linha:
            partes = [p.strip() for p in linha.split("|")]
            if len(partes) == 3:
                resultados.append(partes[1])
    return resultados


# =========================
# IA — INFERÊNCIA EM BATCH
# =========================
def inferir_efeitos_em_lote(titulos, model):
    texto = ""
    for i, t in enumerate(titulos):
        texto += f"[ID:{i}] {t}\n"

    prompt = f"""
Para cada título abaixo, infira:
TARGET | DRUG | EFFECT (uma linha por ID)

Títulos:
{texto}
"""
    out = chamada_ia(model, prompt)

    inferencias = []
    for linha in out.splitlines():
        if "|" in linha:
            inferencias.append(linha.strip())
    return inferencias


# =========================
# PIPELINE PRINCIPAL
# =========================
def minerar_pubmed(termo_base, email, api_key):
    if not api_key:
        return {"termos_indicados": [], "inferencias": []}

    model = configurar_ia(api_key)

    # 1️⃣ Expansão de busca
    termos_expandidos = expandir_termo_com_ia(termo_base, model)
    query = " OR ".join([termo_base] + termos_expandidos)

    # 2️⃣ Buscar PubMed
    pmids = buscar_pubmed_ids(query, email)
    if not pmids:
        return {"termos_indicados": [], "inferencias": []}

    artigos = buscar_titulos_keywords(pmids, email)

    # 3️⃣ NER em batches
    entidades = []
    for i in range(0, len(artigos), BATCH_SIZE):
        lote = artigos[i:i+BATCH_SIZE]
        try:
            entidades.extend(ner_em_lote(lote, model))
            time.sleep(DELAY_ENTRE_BATCH)
        except Exception:
            continue

    # normalização + ranking
    entidades_norm = [normalizar_termo(e) for e in entidades]
    contagem = Counter(entidades_norm)
    termos_indicados = [t for t, _ in contagem.most_common(15)]

    # 4️⃣ Inferência de efeito (somente títulos)
    titulos = [a["titulo"] for a in artigos[:20]]
    inferencias = []
    for i in range(0, len(titulos), BATCH_SIZE):
        lote = titulos[i:i+BATCH_SIZE]
        try:
            inferencias.extend(inferir_efeitos_em_lote(lote, model))
            time.sleep(DELAY_ENTRE_BATCH)
        except Exception:
            continue

    return {
        "termos_indicados": termos_indicados,
        "inferencias": inferencias
    }


# =========================
# FUNÇÕES DE COMPATIBILIDADE (FRONTEND ANTIGO)
# =========================
def buscar_alvos_emergentes_pubmed(alvo, email, usar_ia=True):
    api_key = st.session_state.get("api_key_usuario", "")
    res = minerar_pubmed(alvo, email, api_key)
    return res.get("termos_indicados", [])


def resumir_artigos_pos_busca(titulos, api_key):
    if not api_key or not titulos:
        return ""

    model = configurar_ia(api_key)

    texto = "\n".join(f"- {t}" for t in titulos[:20])
    prompt = f"""
Resuma os principais achados biológicos considerando apenas os títulos abaixo.
Não use abstracts.

Títulos:
{texto}
"""
    return chamada_ia(model, prompt)
