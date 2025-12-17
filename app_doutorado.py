"""
Lemos Lambda: Deep Science Prospector
Copyright (c) 2025 Guilherme Lemos
Licensed under the MIT License.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Author: Guilherme Lemos (Unifesp)
Creation Date: December 2025
"""

import streamlit as st
import pandas as pd
from Bio import Entrez
import time
import plotly.express as px
import re
from deep_translator import GoogleTranslator
from datetime import datetime
import io
import feedparser
import random

# ==========================================
# 0. FUNÇÃO DE AUTOMAÇÃO DE TENDÊNCIAS
# ==========================================
def buscar_alvos_emergentes_pubmed(email):
    Entrez.email = email
    # Query de fronteira 2024-2025: Orphan receptors, Piezos e Ferroptose
    query = '("orphan receptor" OR "GPR" OR "Piezo channel" OR "TAS2R" OR "ferroptosis" OR "SPM mediator") AND ("2024"[Date - Publication] : "2025"[Date - Publication])'
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=20)
        record = Entrez.read(handle)
        if not record["IdList"]: return []
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="abstract", retmode="text")
        texto = handle.read()
        encontrados = re.findall(r'\b[A-Z]{2,6}[0-9]{1,4}\b', texto)
        blacklist = ["DNA", "RNA", "USA", "NCBI", "NIH", "ATP", "AMP", "GDP", "COVID", "SARS"]
        return sorted(list(set([t for t in encontrados if t not in blacklist and len(t) > 2])))
    except:
        return []

# ==========================================
# 1. CONFIGURAÇÃO GLOBAL
# ==========================================
st.set_page_config(page_title="Lemos Lambda", page_icon="λ", layout="wide")

# Inicialização do Session State
if 'alvos_val' not in st.session_state: st.session_state.alvos_val = ""
if 'fonte_val' not in st.session_state: st.session_state.fonte_val = ""
if 'alvo_val' not in st.session_state: st.session_state.alvo_val = ""

# ==========================================
# 2. BANCO DE DADOS: LISTA FIXA GIGANTE
# ==========================================
LISTA_FIXA_LEMOS = [
    "Hydrogen Sulfide (H2S)", "CBS", "CSE", "GYY4137", "AP39", "Nitric Oxide", "Riociguat", "Vericiguat", "Carbon Monoxide (CO)", "HO-1",
    "P2X1 receptor", "P2X3", "P2X7", "P2Y6", "P2Y12", "Adenosine A2A", "FAAH", "MAGL", "Anandamide", "2-AG", "GPR55",
    "KATP channel", "Kir6.1", "Kir6.2", "Glibenclamide", "Cromakalim", "SK channels", "SK3", "Kv7.4", "Retigabine", "BKCa",
    "MALAT1", "HOTAIR", "MEG3", "H19", "GAS5", "miR-29b", "miR-132", "miR-199a", "miR-21", "miR-145", "siRNA therapy",
    "Exosomes", "CD63", "CD9", "CD81", "TSG101", "Alix", "Extracellular Vesicles", "Gap Junctions", "Connexin 43",
    "PD-1", "PD-L1", "CTLA-4", "LAG-3", "TIM-3", "Siglec-8", "Mast Cell Tryptase", "IL-33", "ST2 receptor",
    "Olfactory Receptors", "OR51E2", "OR1D2", "Taste Receptors", "TAS2R", "TAS1R3", "TRPM5",
    "Clock genes", "BMAL1", "CLOCK", "PER1", "PER2", "CRY1", "Rev-erb alpha", "MT1", "MT2",
    "YAP", "TAZ", "Hippo pathway", "Piezo1", "Piezo2", "Integrin beta-1", "FAK", "CTGF", "LOX", "Caveolin-1", "Pirfenidone",
    "HDAC inhibitors", "HDAC1", "Valproic acid", "Vorinostat", "DNMT1", "TET2", "EZH2",
    "Mitochondrial dynamics", "Drp1", "Mfn2", "PGC-1alpha", "Sirtuin-1", "Sirtuin-3", "NAMPT",
    "Ferroptosis", "GPX4", "SLC7A11", "Pyroptosis", "Gasdermin D", "Necroptosis", "RIPK1", "RIPK3",
    "Microplastics", "Nanoplastics", "Bisphenol S", "Phthalates", "Glyphosate", "Acrolein", "Cadmium",
    "TMEM16A", "HCN1", "HCN4", "Kv7.1", "TREK-1", "TRAAK", "TRPML1"
]
SUGESTOES_TEXTO = ", ".join(LISTA_FIXA_LEMOS)

# ==========================================
# 3. INTERFACE (UI) E DICIONÁRIO
# ==========================================
TEXTOS = {
    "pt": {
        "titulo": "λ Lemos Lambda: Deep Science",
        "email_label": "Seu E-mail (NCBI):",
        "btn_restaurar": "📥 Termos indicados",
        "btn_minerar": "⛏️ Minerar 'Blue Oceans'",
        "label_lista": "**Palavras-chave:**",
        "btn_trend": "🔍 Tendências (2025)",
        "btn_avanco": "🚀 Rumo ao Avanço",
        "col_ratio": "Ratio (Oceano Azul)",
        "raio_x": "🔎 Raio-X Literário"
    },
    "en": {
        "titulo": "λ Lemos Lambda: Deep Science",
        "email_label": "Your E-mail (NCBI):",
        "btn_restaurar": "📥 Termos indicados",
        "btn_minerar": "⛏️ Mine 'Blue Oceans'",
        "label_lista": "**Keywords:**",
        "btn_trend": "🔍 Trends (2025)",
        "btn_avanco": "🚀 Launch Analysis",
        "col_ratio": "Blue Ocean Ratio",
        "raio_x": "🔎 Literary X-Ray"
    }
}

lang_sel = st.sidebar.radio("Idioma:", ["🇧🇷 Português", "🇺🇸 English"])
lang_code = "pt" if "Português" in lang_sel else "en"
t = TEXTOS[lang_code]

st.title(t["titulo"])

# Sidebar: Configuração
st.sidebar.header("1. Configurações")
email_user = st.sidebar.text_input(t["email_label"])
anos_range = st.sidebar.slider("Período:", 1990, 2025, (2010, 2025))

st.sidebar.header("2. Órgãos")
st.session_state.fonte_val = st.sidebar.text_input("Fonte:", value=st.session_state.fonte_val)
st.session_state.alvo_val = st.sidebar.text_input("Alvo:", value=st.session_state.alvo_val)

st.sidebar.header("3. Palavras-chave")
# Lógica de sincronização para evitar erro de API
alvos_input = st.sidebar.text_area(t["label_lista"], value=st.session_state.alvos_val, height=200)
st.session_state.alvos_val = alvos_input

# BOTÃO DE TENDÊNCIAS (Abaixo da caixa de palavras-chave)
if st.sidebar.button(t["btn_trend"]):
    if email_user:
        with st.sidebar.status("Buscando..."):
            novos = buscar_alvos_emergentes_pubmed(email_user)
            if novos:
                txt = ", ".join(novos)
                st.session_state.alvos_val = (st.session_state.alvos_val.strip(", ") + ", " + txt).strip(", ")
                st.rerun()

# Botões de Termos Indicados e Mineração
col_b1, col_b2 = st.sidebar.columns(2)
if col_b1.button(t["btn_restaurar"]):
    st.session_state.alvos_val = SUGESTOES_TEXTO
    st.rerun()

if col_b2.button(t["btn_minerar"]):
    # Lógica de mineração (Filtra a lista fixa contra o órgão alvo)
    with st.sidebar.status("Minerando..."):
        encontrados = []
        for termo in LISTA_FIXA_LEMOS:
            Entrez.email = email_user
            q = f"({termo}) AND ({st.session_state.alvo_val}) AND {anos_range[0]}:{anos_range[1]}[DP]"
            try:
                c = int(Entrez.read(Entrez.esearch(db="pubmed", term=q, retmax=0))["Count"])
                if 0 <= c < 150: encontrados.append(termo)
            except: continue
        st.session_state.alvos_val = ", ".join(encontrados)
        st.rerun()

# ==========================================
# 4. PROCESSAMENTO E RESULTADOS
# ==========================================
if st.sidebar.button(t["btn_avanco"], type="primary"):
    lista_alvos = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
    res = []
    bar = st.progress(0)
    for i, alvo in enumerate(lista_alvos):
        # Consulta Pubmed
        Entrez.email = email_user
        qf = f"({alvo}) AND ({st.session_state.fonte_val})"
        qa = f"({alvo}) AND ({st.session_state.alvo_val})"
        nf = int(Entrez.read(Entrez.esearch(db="pubmed", term=qf, retmax=0))["Count"])
        na = int(Entrez.read(Entrez.esearch(db="pubmed", term=qa, retmax=0))["Count"])
        ratio = nf/na if na > 0 else nf
        res.append({"Alvo": alvo, "Ratio": ratio, "Status": "💎 DIAMANTE" if ratio > 10 else "🥇 Ouro"})
        bar.progress((i+1)/len(lista_alvos))
    st.session_state.dados = pd.DataFrame(res).sort_values(by="Ratio", ascending=False)

if 'dados' in st.session_state:
    df = st.session_state.dados
    st.plotly_chart(px.bar(df.head(15), x="Alvo", y="Ratio", color="Status", title=t["col_ratio"]), use_container_width=True)
    st.dataframe(df, use_container_width=True, hide_index=True)
    
    st.divider()
    st.header(t["raio_x"])
    sel = st.selectbox("Selecione para ver a conclusão da IA:", df["Alvo"].tolist())
    if st.button("Ver Conclusão"):
        with st.spinner("IA analisando..."):
            query = f"({sel}) AND ({st.session_state.alvo_val})"
            h = Entrez.esearch(db="pubmed", term=query, retmax=1, sort="relevance")
            pmid = Entrez.read(h)["IdList"][0]
            abstract = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="text").read()
            traducao = GoogleTranslator(source='auto', target=lang_code).translate(abstract[-500:])
            st.success(traducao)
            st.link_button("Ver artigo completo", f"https://pubmed.ncbi.nlm.nih.gov/{pmid}")
