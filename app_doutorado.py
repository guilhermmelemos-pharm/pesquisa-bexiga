"""
Lemos Lambda: Deep Science Prospector
Copyright (c) 2025 Guilherme Lemos
Licensed under the MIT License.
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
# 0. FUNÇÃO DE AUTOMAÇÃO DE TERMOS (NOVO)
# ==========================================
def buscar_alvos_emergentes_pubmed(email):
    Entrez.email = email
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

st.markdown("""
    <style>
    div[data-testid="stImage"] img { height: 150px !important; object-fit: cover !important; border-radius: 8px !important; }
    .stButton button { width: 100%; }
    div[data-testid="stVerticalBlock"] > div { gap: 0.5rem; }
    </style>
""", unsafe_allow_html=True)

if 'alvos_val' not in st.session_state: st.session_state.alvos_val = ""
if 'fonte_val' not in st.session_state: st.session_state.fonte_val = ""
if 'alvo_val' not in st.session_state: st.session_state.alvo_val = ""
if 'news_index' not in st.session_state: st.session_state.news_index = 0

# ==========================================
# 2. DICIONÁRIO DE TRADUÇÃO E BANCO DE DADOS
# ==========================================
CANDIDATOS_MINERACAO = ["GPR37", "GPR17", "GPR55", "GPR84", "GPR35", "GPR183", "GPR119", "Piezo1", "Piezo2", "TREK-1", "TASK-1", "TAS2R10", "TAS2R14", "Ferroptosis", "GPX4", "SLC7A11", "Resolvin D1", "Maresin 1", "P2X4", "P2X7"]
LISTA_ALVOS_PRONTA = ", ".join(CANDIDATOS_MINERACAO)

TEXTOS = {
    "pt": {
        "titulo_desk": "λ Lemos Lambda: Deep Science",
        "subtitulo": "**Ferramenta de Prospecção de Alto Impacto**",
        "titulo_mob": "📱 Lemos Pocket",
        "credenciais": "1. Credenciais",
        "email_label": "Seu E-mail:",
        "periodo": "📅 Período:",
        "config": "2. Configuração (Órgãos)",
        "label_fonte": "**Fonte (Órgão, tecido, célula similar):**",
        "label_alvo": "**Alvo (Órgão de interesse):**",
        "btn_setup": "🎓 Doutorado Guilherme Lemos",
        "sec_alvos": "3. Palavras-chave",
        "label_lista": "**Palavras-chave de Pesquisa:**",
        "btn_restaurar": "📥 Termos indicados",
        "btn_minerar": "⛏️ Minerar 'Blue Oceans'",
        "btn_trend": "🔍 Injetar Tendências (2025)",
        "btn_avanco": "🚀 Rumo ao Avanço",
        "col_artigos": "Artigos",
        "col_ratio": "Ratio (Fonte/Alvo)",
        "raio_x": "🔎 Raio-X",
        "baixar": "📥 Baixar Planilha"
    },
    "en": {
        "titulo_desk": "λ Lemos Lambda: Deep Science",
        "subtitulo": "**High Impact Prospecting Tool**",
        "credenciais": "1. Credentials",
        "email_label": "Your E-mail:",
        "periodo": "📅 Timeframe:",
        "config": "2. Configuration (Organs)",
        "label_fonte": "**Source (Organ, tissue, similar cell):**",
        "label_alvo": "**Target (Organ of interest):**",
        "sec_alvos": "3. Keywords",
        "label_lista": "**Research Keywords:**",
        "btn_restaurar": "📥 Termos indicados",
        "btn_minerar": "⛏️ Mine 'Blue Oceans'",
        "btn_trend": "🔍 Inject Trends (2025)",
        "btn_avanco": "🚀 Launch Analysis",
        "col_ratio": "Ratio (Source/Target)",
        "raio_x": "🔎 X-Ray",
        "baixar": "📥 Download CSV"
    }
}

# ==========================================
# 3. FUNÇÕES DE SUPORTE
# ==========================================
@st.cache_data(ttl=3600)
def buscar_todas_noticias(lang_code):
    feeds = [{"url": "https://www.sciencedaily.com/rss/health_medicine/pharmacology.xml", "lang": "🇺🇸"}]
    noticias = []
    for fonte in feeds:
        try:
            feed = feedparser.parse(fonte["url"])
            for entry in feed.entries[:3]:
                noticias.append({"titulo": entry.title, "link": entry.link, "fonte": "ScienceDaily", "img": "https://images.unsplash.com/photo-1532094349884-543bc11b234d?w=400", "bandeira": fonte["lang"]})
        except: continue
    return noticias

@st.fragment(run_every=60) 
def exibir_radar_cientifico(lang_code):
    news_list = buscar_todas_noticias(lang_code)
    if not news_list: return
    with st.container(border=True):
        st.caption(f"📡 **Radar Científico**")
        cols = st.columns(3)
        for i, n in enumerate(news_list[:3]):
            with cols[i]: st.markdown(f"**{n['titulo'][:60]}...**"); st.link_button("Ler", n['link'])

def carregar_setup_lemos(t):
    st.session_state.alvos_val = LISTA_ALVOS_PRONTA
    st.session_state.fonte_val = "Brain OR Kidney OR Liver"
    st.session_state.alvo_val = "Bladder OR Vesical OR Cystitis"

def carregar_alvos_apenas(t): st.session_state.alvos_val = LISTA_ALVOS_PRONTA

def consultar_pubmed_count(termo_farmaco, termo_orgao, email, y_start, y_end):
    if not email: return 0
    Entrez.email = email
    query = f"({termo_farmaco}) AND ({termo_orgao}) AND {y_start}:{y_end}[DP]" if termo_orgao else f"({termo_farmaco}) AND {y_start}:{y_end}[DP]"
    try: return int(Entrez.read(Entrez.esearch(db="pubmed", term=query, retmax=0))["Count"])
    except: return 0

# ==========================================
# 4. INTERFACE (UI)
# ==========================================
lang_opt = st.sidebar.radio("Language:", ["🇧🇷 Português", "🇺🇸 English"])
lang = "pt" if "Português" in lang_opt else "en"
t = TEXTOS[lang]
modo = st.sidebar.radio("📱 Mode:", ["Desktop", "Mobile"], index=0)

if modo == "Desktop":
    st.title(t["titulo_desk"])
    st.sidebar.header(t["credenciais"])
    email_user = st.sidebar.text_input(t["email_label"], key="email_desk")
    anos = st.sidebar.slider(t["periodo"], 1990, 2025, (2010, 2025), key="anos_desk")
    
    st.sidebar.header(t["config"])
    t_fonte = st.sidebar.text_input(t["label_fonte"], value=st.session_state.fonte_val, key="fv")
    t_alvo = st.sidebar.text_input(t["label_alvo"], value=st.session_state.alvo_val, key="av")
    st.sidebar.button(t["btn_setup"], on_click=carregar_setup_lemos, args=(t,))

    st.sidebar.header(t["sec_alvos"])
    alvos_in = st.sidebar.text_area(t["label_lista"], value=st.session_state.alvos_val, height=150)
    st.session_state.alvos_val = alvos_in

    if st.sidebar.button(t["btn_trend"]):
        if email_user:
            novos = buscar_alvos_emergentes_pubmed(email_user)
            st.session_state.alvos_val = (st.session_state.alvos_val.strip(", ") + ", " + ", ".join(novos)).strip(", ")
            st.rerun()

    b1, b2 = st.sidebar.columns(2)
    b1.button(t["btn_restaurar"], on_click=carregar_alvos_apenas, args=(t,))
    
    if st.sidebar.button(t["btn_avanco"], type="primary"):
        lst = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
        res = []
        bar = st.progress(0)
        for i, item in enumerate(lst):
            nf = consultar_pubmed_count(item, t_fonte, email_user, anos[0], anos[1])
            na = consultar_pubmed_count(item, t_alvo, email_user, anos[0], anos[1])
            ratio = nf/na if na > 0 else nf
            res.append({"Alvo": item, "Status": "💎" if ratio > 10 else "🥇", "Ratio": ratio, "Fonte": nf, "Alvo_Interest": na})
            bar.progress((i+1)/len(lst))
        st.session_state['dados_res'] = pd.DataFrame(res).sort_values(by="Ratio", ascending=False)
        st.rerun()

    if 'dados_res' in st.session_state:
        df = st.session_state['dados_res']
        st.plotly_chart(px.bar(df.head(20), x="Alvo", y="Ratio", color="Status", title=t["col_ratio"]), use_container_width=True)
        st.dataframe(df, use_container_width=True, hide_index=True)
