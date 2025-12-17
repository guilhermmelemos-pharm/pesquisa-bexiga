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
# 0. FUNÇÃO DE AUTOMAÇÃO FOCADA NO ÓRGÃO (NOVO)
# ==========================================
def buscar_alvos_emergentes_focados(orgao_alvo, email):
    """Minera o PubMed em tempo real buscando alvos farmacológicos para o órgão específico."""
    if not orgao_alvo or not email:
        return []
    Entrez.email = email
    # Busca focada: Órgão + termos de fronteira 2024-2025
    query = f"({orgao_alvo}) AND (receptor OR channel OR protein OR signaling) AND (\"2024\"[Date - Publication] : \"2025\"[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=30, sort="relevance")
        record = Entrez.read(handle)
        ids = record["IdList"]
        if not ids: return []
        
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
        texto = handle.read()
        
        # Regex para capturar siglas proteicas (Ex: GPR35, TRPV4, P2X7)
        encontrados = re.findall(r'\b[A-Z]{2,6}[0-9]{0,4}\b', texto)
        blacklist = ["DNA", "RNA", "USA", "NCBI", "NIH", "ATP", "AMP", "GDP", "COVID", "SARS", "PMID", "DOI", "FAPESP"]
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
# 2. BANCO DE DADOS E DICIONÁRIO
# ==========================================
CANDIDATOS_MINERACAO = [
    "Hydrogen Sulfide (H2S)", "CBS", "CSE", "GYY4137", "AP39", "Nitric Oxide", "Riociguat", "Vericiguat", "Carbon Monoxide (CO)", "HO-1",
    "P2X1 receptor", "P2X3", "P2X7", "P2Y6", "P2Y12", "Adenosine A2A", "FAAH", "MAGL", "Anandamide", "2-AG", "GPR55",
    "KATP channel", "Kir6.1", "Kir6.2", "Glibenclamide", "Cromakalim", "SK channels", "SK3", "Kv7.4", "Retigabine", "BKCa",
    "MALAT1", "HOTAIR", "MEG3", "H19", "GAS5", "miR-29b", "miR-132", "miR-199a", "miR-21", "miR-145", "siRNA therapy",
    "Exosomes", "CD63", "CD9", "CD81", "TSG101", "Alix", "Extracellular Vesicles", "Gap Junctions", "Connexin 43",
    "PD-1", "PD-L1", "CTLA-4", "LAG-3", "TIM-3", "Siglec-8", "Mast Cell Tryptase", "IL-33", "ST2 receptor",
    "Olfactory Receptors", "OR51E2", "OR1D2", "Taste Receptors", "TAS2R", "TAS1R3", "TRPM5",
    "Clock genes", "BMAL1", "CLOCK", "PER1", "PER2", "CRY1", "Rev-erb alpha", "MT1", "MT2",
    "YAP", "TAZ", "Hippo pathway", "Piezo1", "Piezo2", "Integrin beta-1", "FAAK", "CTGF", "LOX", "Caveolin-1", "Pirfenidone",
    "HDAC inhibitors", "HDAC1", "Valproic acid", "Vorinostat", "DNMT1", "TET2", "EZH2",
    "Mitochondrial dynamics", "Drp1", "Mfn2", "PGC-1alpha", "Sirtuin-1", "Sirtuin-3", "NAMPT",
    "Ferroptosis", "GPX4", "SLC7A11", "Pyroptosis", "Gasdermin D", "Necroptosis", "RIPK1", "RIPK3",
    "Microplastics", "Nanoplastics", "Bisphenol S", "Phthalates", "Glyphosate", "Acrolein", "Cadmium",
    "TMEM16A", "HCN1", "HCN4", "Kv7.1", "TREK-1", "TRAAK", "TRPML1"
]
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
        "label_fonte": "**Fonte (Orgão, tecido, célula similar):**",
        "holder_fonte": "Ex: Kidney...",
        "label_alvo": "**Alvo (Orgão de interesse):**",
        "holder_alvo": "Ex: Bladder...",
        "btn_setup": "🎓 Doutorado Guilherme Lemos",
        "toast_setup": "Setup 'Deep Science' Carregado!",
        "sec_alvos": "3. Palavras-chave",
        "expander_upload": "📂 Importar Biblioteca (.csv/.txt)",
        "toast_upload": "Biblioteca importada!",
        "label_lista": "**Palavras-chave de Pesquisa:**",
        "holder_lista": "Insira os alvos ou use a automação...",
        "btn_restaurar": "📥 Termos indicados",
        "toast_restaurar": "Lista Inovadora Restaurada!",
        "btn_minerar": "⛏️ Minerar 'Blue Oceans'",
        "btn_trend": "🔍 Injetar Tendências (2025)",
        "toast_aviso_minerar": "⚠️ Preencha o 'Alvo' e 'E-mail' para minerar!",
        "prog_minerar": "⛏️ Procurando termos chave...",
        "prog_testando": "⛏️ Analisando: {termo} ({count} artigos)",
        "toast_sucesso_minerar": "✅ {qtd} termos encontrados!",
        "btn_avanco": "🚀 Rumo ao Avanço",
        "erro_email": "E-mail obrigatório!",
        "aviso_lista": "Lista de Palavras-chave vazia!",
        "prog_investigando": "⏳ Investigando {atual}/{total}: {alvo}",
        "analise_pronta": "✅ Análise Pronta. Destaque: **{top}**.",
        "col_artigos": "Artigos",
        "col_global": "Global",
        "col_ratio": "Ratio",
        "filtro": "🔍 Filtro:",
        "grafico_qtd": "📊 Qtd. no Gráfico:",
        "raio_x": "🔎 Raio-X",
        "btn_ler": "Ler Artigos",
        "btn_scholar": "🎓 Google Scholar",
        "sem_artigos": "Zero artigos encontrados.",
        "lendo": "Buscando e Traduzindo...",
        "baixar": "📥 Baixar Planilha",
        "citar_titulo": "📄 Como Citar",
        "citar_texto": "Lemos, G. (2025). Lemos Lambda: Deep Science Prospector [Software]. Versão 1.0.0. DOI: 10.5281/zenodo.17958507",
        "link_doi": "🔗 Ver no Zenodo (DOI)"
    },
    "en": {
        "titulo_desk": "λ Lemos Lambda: Deep Science",
        "subtitulo": "**High Impact Prospecting Tool**",
        "titulo_mob": "📱 Lemos Pocket",
        "credenciais": "1. Credentials",
        "email_label": "Your E-mail:",
        "periodo": "📅 Timeframe:",
        "config": "2. Configuration (Organs)",
        "label_fonte": "**Source (Organ, tissue, similar cell):**",
        "holder_fonte": "Ex: Kidney...",
        "label_alvo": "**Target (Organ of interest):**",
        "holder_alvo": "Ex: Bladder...",
        "btn_setup": "🎓 Guilherme Lemos PhD Setup",
        "toast_setup": "'Deep Science' Setup Loaded!",
        "sec_alvos": "3. Keywords",
        "expander_upload": "📂 Import Library (.csv/.txt)",
        "toast_upload": "Library imported!",
        "label_lista": "**Research Keywords:**",
        "holder_lista": "Load keywords...",
        "btn_restaurar": "📥 Termos indicados",
        "toast_restaurar": "Innovative List Restored!",
        "btn_minerar": "⛏️ Mine 'Blue Oceans'",
        "btn_trend": "🔍 Inject Trends (2025)",
        "toast_aviso_minerar": "⚠️ Fill in 'Target' and 'E-mail' to mine!",
        "prog_minerar": "⛏️ Searching for key terms...",
        "prog_testando": "⛏️ Analyzing: {termo} ({count} papers)",
        "toast_sucesso_minerar": "✅ {qtd} terms found!",
        "btn_avanco": "🚀 Launch Analysis",
        "erro_email": "E-mail required!",
        "aviso_lista": "Keyword list is empty!",
        "prog_investigando": "⏳ Investigating {atual}/{total}: {alvo}",
        "analise_pronta": "✅ Analysis Ready. Highlight: **{top}**.",
        "col_artigos": "Papers",
        "col_global": "Global",
        "col_ratio": "Ratio",
        "filtro": "🔍 Filter:",
        "grafico_qtd": "📊 Chart Qty:",
        "raio_x": "🔎 X-Ray",
        "btn_ler": "Read Papers",
        "btn_scholar": "🎓 Google Scholar",
        "sem_artigos": "Zero papers found.",
        "lendo": "Searching and Translating...",
        "baixar": "📥 Download CSV",
        "citar_titulo": "📄 How to Cite",
        "citar_texto": "Lemos, G. (2025). Lemos Lambda: Deep Science Prospector [Software]. Version 1.0.0. DOI: 10.5281/zenodo.17958507",
        "link_doi": "🔗 View on Zenodo (DOI)"
    }
}

# ==========================================
# 3. FUNÇÕES DE SUPORTE
# ==========================================
@st.cache_data(ttl=3600)
def buscar_todas_noticias(lang_code):
    feeds = [
        {"url": "https://www.sciencedaily.com/rss/health_medicine/pharmacology.xml", "lang": "🇺🇸"},
        {"url": "https://www.nature.com/nbt.rss", "lang": "🇬🇧"},
        {"url": "https://agencia.fapesp.br/rss/", "lang": "🇧🇷"},
    ]
    noticias = []
    translator = GoogleTranslator(source='auto', target=lang_code)
    for fonte in feeds:
        try:
            feed = feedparser.parse(fonte["url"])
            for entry in feed.entries[:3]:
                img_url = "https://images.unsplash.com/photo-1532094349884-543bc11b234d?w=400"
                titulo = entry.title
                if lang_code == 'pt' and fonte["lang"] != "🇧🇷":
                    try: titulo = translator.translate(titulo)
                    except: pass
                noticias.append({"titulo": titulo, "link": entry.link, "fonte": feed.feed.title[:20], "img": img_url, "bandeira": fonte["lang"]})
        except: continue
    return noticias

@st.fragment(run_every=60) 
def exibir_radar_cientifico(lang_code):
    news_list = buscar_todas_noticias(lang_code)
    if not news_list: return
    batch = news_list[:3]
    with st.container(border=True):
        st.caption(f"📡 **Radar Científico**")
        cols = st.columns(3)
        for i, n in enumerate(batch):
            with cols[i]:
                st.markdown(f"**{n['titulo'][:60]}...**")
                st.link_button("Ler", n['link'])

def carregar_setup_lemos(t):
    st.session_state.alvos_val = LISTA_ALVOS_PRONTA
    st.session_state.fonte_val = "Brain OR Kidney OR Liver"
    st.session_state.alvo_val = "Bladder OR Vesical OR Cystitis"
    st.toast(t["toast_setup"], icon="🧬")

def carregar_alvos_apenas(t): 
    st.session_state.alvos_val = LISTA_ALVOS_PRONTA
    st.toast(t["toast_restaurar"], icon="✨")

def limpar_campo_fonte(): st.session_state.fonte_val = ""
def limpar_campo_alvo(): st.session_state.alvo_val = ""
def limpar_campo_alvos(): st.session_state.alvos_val = ""

def consultar_pubmed_count(termo_farmaco, termo_orgao, email, y_start, y_end):
    if not email: return -1
    Entrez.email = email
    query = f"({termo_farmaco}) AND ({termo_orgao}) AND {y_start}:{y_end}[DP]" if termo_orgao else f"({termo_farmaco}) AND {y_start}:{y_end}[DP]"
    try: return int(Entrez.read(Entrez.esearch(db="pubmed", term=query, retmax=0))["Count"])
    except: return -1

def buscar_resumos_detalhados(termo_farmaco, termo_orgao, email, y_start, y_end, lang_target, limit=3):
    query = f"({termo_farmaco}) AND ({termo_orgao}) AND {y_start}:{y_end}[DP]" if termo_orgao else f"({termo_farmaco}) AND {y_start}:{y_end}[DP]"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=limit, sort="relevance")
        ids = Entrez.read(handle)["IdList"]
        if not ids: return []
        records = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text").read().split("\n\n")
        artigos = []
        for art_text in records:
            art_data = {"PMID": "N/A", "Title": "S/T", "Abstract": ""}
            for line in art_text.split("\n"):
                tag = line[:4].strip()
                if tag=="PMID": art_data["PMID"]=line[6:]
                elif tag=="TI": art_data["Title"]=line[6:]
                elif tag=="AB": art_data["Abstract"]=line[6:]
            if art_data["PMID"]!="N/A":
                res = art_data["Abstract"][-400:] if art_data["Abstract"] else "..."
                art_data["Resumo_IA"] = GoogleTranslator(source='auto', target=lang_target).translate(res)
                artigos.append(art_data)
        return artigos
    except: return []

# ==========================================
# 4. INTERFACE (UI)
# ==========================================
lang_opt = st.sidebar.radio("Language:", ["🇧🇷 Português", "🇺🇸 English"])
lang = "pt" if "Português" in lang_opt else "en"
t = TEXTOS[lang]
modo = st.sidebar.radio("📱 Mode:", ["Desktop", "Mobile"], index=0)

if modo == "Desktop":
    st.title(t["titulo_desk"])
    st.markdown(t["subtitulo"])
    if 'dados_desk' not in st.session_state: exibir_radar_cientifico(lang)
    
    st.sidebar.header(t["credenciais"])
    email_user = st.sidebar.text_input(t["email_label"], placeholder="pesquisador@unifesp.br", key="email_desk")
    anos = st.sidebar.slider(t["periodo"], 1990, 2025, (2010, 2025), key="anos_desk")
    
    st.sidebar.markdown("---")
    st.sidebar.header(t["config"])
    st.sidebar.markdown(t["label_fonte"])
    c1, c2 = st.sidebar.columns([6, 1], vertical_alignment="bottom")
    with c1: t_fonte = st.text_input("Fonte", key="fonte_val", placeholder=t["holder_fonte"], label_visibility="collapsed")
    with c2: st.button("🗑️", key="del_f", on_click=limpar_campo_fonte)

    st.sidebar.markdown(t["label_alvo"])
    c3, c4 = st.sidebar.columns([6, 1], vertical_alignment="bottom")
    with c3: t_alvo = st.text_input("Alvo", key="alvo_val", placeholder=t["holder_alvo"], label_visibility="collapsed")
    with c4: st.button("🗑️", key="del_a", on_click=limpar_campo_alvo)
    
    st.sidebar.button(t["btn_setup"], on_click=carregar_setup_lemos, args=(t,))
    
    st.sidebar.markdown("---")
    st.sidebar.header(t["sec_alvos"])
    st.sidebar.markdown(t["label_lista"])
    c5, c6 = st.sidebar.columns([6, 1], vertical_alignment="bottom")
    with c5: 
        alvos_in = st.text_area(t["label_lista"], value=st.session_state.alvos_val, height=150, label_visibility="collapsed")
        st.session_state.alvos_val = alvos_in 
    with c6: st.button("🗑️", key="del_l", on_click=limpar_campo_alvos)

    # BOTÃO TENDÊNCIAS FOCADO NO ÓRGÃO
    if st.sidebar.button(t["btn_trend"], key="trend_desk"):
        if email_user and t_alvo:
            with st.sidebar.status(f"Minerando termos para {t_alvo}..."):
                novos = buscar_alvos_emergentes_focados(t_alvo, email_user)
                if novos:
                    txt_novos = ", ".join(novos)
                    if st.session_state.alvos_val:
                        st.session_state.alvos_val = (st.session_state.alvos_val.strip(", ") + ", " + txt_novos)
                    else:
                        st.session_state.alvos_val = txt_novos
                    st.rerun()
        else: st.sidebar.error("E-mail e Órgão Alvo necessários!")

    b1, b2 = st.sidebar.columns(2)
    b1.button(t["btn_restaurar"], on_click=carregar_alvos_apenas, args=(t,))
    b2.button(t["btn_minerar"], on_click=minerar_blue_oceans, args=(t_alvo, email_user, t))
    
    st.sidebar.markdown("---")

    if st.sidebar.button(t["btn_avanco"], type="primary"):
        if not email_user: st.error(t["erro_email"])
        else:
            lst = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
            res = []; pg = st.empty(); bar = st.progress(0)
            for i, item in enumerate(lst):
                pg.text(t["prog_investigando"].format(atual=i+1, total=len(lst), alvo=item))
                nf = consultar_pubmed_count(item, t_fonte, email_user, anos[0], anos[1])
                na = consultar_pubmed_count(item, t_alvo, email_user, anos[0], anos[1])
                pot = nf/na if na > 0 else nf
                res.append({"Alvo": item, "Status": "💎 DIAMANTE" if pot > 10 else "🥇 Ouro", "Ratio": pot, "Fonte": nf, "Alvo_Interest": na})
                bar.progress((i+1)/len(lst))
            pg.empty()
            st.session_state['dados_desk'] = pd.DataFrame(res).sort_values(by="Ratio", ascending=False)
            st.rerun()

    if 'dados_desk' in st.session_state:
        df = st.session_state['dados_desk']
        st.plotly_chart(px.bar(df.head(20), x="Alvo", y="Ratio", color="Status"), use_container_width=True)
        st.dataframe(df, use_container_width=True, hide_index=True)
        st.divider()
        sel = st.selectbox("Raio-X:", sorted(df['Alvo'].unique().tolist()))
        if st.button(t["btn_ler"]):
            arts = buscar_resumos_detalhados(sel, t_alvo, email_user, anos[0], anos[1], lang)
            for a in arts:
                with st.expander(a['Title']):
                    st.success(a['Resumo_IA'])

elif modo == "Mobile":
    st.title(t["titulo_mob"])
    email_mob = st.text_input(t["email_label"], key="emob")
    alvo_mob = st.text_input(t["label_alvo"], key="amob")
    alvos_m = st.text_area(t["label_lista"], value=st.session_state.alvos_val, key="alm")
    st.session_state.alvos_val = alvos_m

    if st.button(t["btn_trend"]):
        n = buscar_alvos_emergentes_focados(alvo_mob, email_mob)
        st.session_state.alvos_val += ", ".join(n); st.rerun()

    if st.button(t["btn_avanco"], type="primary"):
        st.info("Processando...")
