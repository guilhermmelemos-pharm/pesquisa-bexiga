import streamlit as st
from Bio import Entrez
from google import genai
import re
import time

# =========================
# CONFIGURAÇÕES GERAIS
# =========================
Entrez.email = "pesquisador_guest@unifesp.br"

MODELOS_GEMINI = [
    "models/gemini-2.5-flash",
    "models/gemini-2.0-flash",
    "models/gemini-flash-latest"
]

CLINICAL_BLACKLIST = [
    "cancer", "carcinoma", "tumor", "diagnosis", "therapy", "surgery",
    "syndrome", "benign", "malignant", "case", "risk", "outcomes",
    "prostate", "overactive", "urothelial", "renal", "kidney",
    "immunotherapy", "bacillus", "calmette", "guérin", "hbp", "bph"
]


# =========================
# GEMINI CLIENT
# =========================
def get_gemini_client(api_key):
    return genai.Client(api_key=api_key)


# =========================
# CURADORIA MOLECULAR (IA)
# =========================
def curadoria_molecular_com_ia(titulo, abstract, api_key, lang="pt"):
    if not api_key:
        return []

    # Abstract longo → truncamento inteligente
    texto_base = abstract.strip() if abstract else ""
    texto_base = texto_base[:1200] if len(texto_base) > 1200 else texto_base

    prompt = f"""
Você é um PhD em Farmacologia Molecular.

TAREFA:
Extraia APENAS entidades moleculares RELEVANTES para farmacologia experimental.

INCLUIR:
- Genes
- Receptores
- Canais iônicos
- Enzimas
- Fármacos
- Ligantes
- Vias de sinalização

EXCLUIR TOTALMENTE:
- Doenças
- Síndromes
- Diagnósticos
- Procedimentos clínicos
- Termos anatômicos
- Populações ou sexo
- Termos médicos gerais

TEXTO:
TÍTULO: {titulo}
RESUMO: {texto_base}

FORMATO DE SAÍDA:
Lista simples separada por vírgula.
Se não houver alvos moleculares, responda apenas: NONE
"""

    client = get_gemini_client(api_key)

    for modelo in MODELOS_GEMINI:
        try:
            response = client.models.generate_content(
                model=modelo,
                contents=prompt
            )

            texto = response.text.strip()

            if texto.upper() == "NONE":
                return []

            candidatos = [t.strip() for t in texto.split(",") if len(t.strip()) > 2]

            # Filtro extra de segurança
            filtrados = [
                c for c in candidatos
                if not any(b in c.lower() for b in CLINICAL_BLACKLIST)
            ]

            return list(dict.fromkeys(filtrados))

        except Exception:
            time.sleep(0.8)
            continue

    # FALHA TOTAL → retorna vazio (nunca lixo)
    return []


# =========================
# RESUMO TÉCNICO (IA)
# =========================
def resumo_tecnico_com_ia(titulo, api_key, lang="pt"):
    if not api_key:
        return "IA desativada."

    prompt = f"""
Resuma tecnicamente este artigo em até 25 palavras.
Foque em mecanismo molecular ou farmacológico.

TÍTULO: {titulo}
"""

    client = get_gemini_client(api_key)

    for modelo in MODELOS_GEMINI:
        try:
            r = client.models.generate_content(model=modelo, contents=prompt)
            return r.text.strip()
        except Exception:
            time.sleep(0.5)

    return "Resumo indisponível (IA ocupada)."


# =========================
# PUBMED – BUSCA RESUMOS
# =========================
def buscar_resumos_pubmed(termo, ano_ini=2015, ano_fim=2025):
    query = f"({termo}) AND ({ano_ini}:{ano_fim}[Date - Publication])"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=3)
    record = Entrez.read(handle)
    handle.close()

    artigos = []
    if not record["IdList"]:
        return artigos

    handle = Entrez.efetch(
        db="pubmed",
        id=record["IdList"],
        rettype="medline",
        retmode="text"
    )
    dados = handle.read()
    handle.close()

    for bloco in dados.split("\n\nPMID-"):
        titulo, resumo = "", ""
        for linha in bloco.split("\n"):
            if linha.startswith("TI  - "):
                titulo = linha.replace("TI  - ", "").strip()
            if linha.startswith("AB  - "):
                resumo = linha.replace("AB  - ", "").strip()
        if titulo:
            artigos.append({"titulo": titulo, "resumo": resumo})

    return artigos


# =========================
# STREAMLIT UI
# =========================
st.set_page_config(page_title="Curadoria Molecular – Bexiga", layout="wide")

st.title("🔬 Curadoria Molecular com Gemini (Farmacologia da Bexiga)")

api_key = st.text_input("🔑 Google AI API Key", type="password")
termo = st.text_input("🔍 Termo PubMed", value="bladder mechanotransduction")

if st.button("Analisar"):
    artigos = buscar_resumos_pubmed(termo)

    if not artigos:
        st.warning("Nenhum artigo encontrado.")
    else:
        for art in artigos:
            st.subheader(art["titulo"])

            resumo = resumo_tecnico_com_ia(art["titulo"], api_key)
            st.markdown(f"**Resumo técnico:** {resumo}")

            alvos = curadoria_molecular_com_ia(
                art["titulo"],
                art["resumo"],
                api_key
            )

            st.markdown("**Alvos moleculares identificados:**")
            st.write(alvos if alvos else "Nenhum alvo molecular relevante.")
