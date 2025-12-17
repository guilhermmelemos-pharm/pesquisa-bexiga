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
# 0. FUNÇÕES DE INTELIGÊNCIA (DESCOBERTA)
# ==========================================
def buscar_alvos_emergentes(email):
    Entrez.email = email
    query = '("orphan receptor" OR "neglected target" OR "novel GPCR" OR "rarely studied") AND ("2024"[Date - Publication] : "2025"[Date - Publication])'
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=30)
        record = Entrez.read(handle)
        if not record["IdList"]: return []
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="abstract", retmode="text")
        texto = handle.read()
        encontrados = re.findall(r'\b[A-Z]{2,6}[0-9]{1,4}\b', texto)
        blacklist = ["DNA", "RNA", "USA", "NCBI", "NIH", "ATP", "AMP", "GDP", "COVID", "SARS", "FAPESP", "PMC", "UI", "UK", "USA"]
        return sorted(list(set([t for t in encontrados if t not in blacklist and len(t) > 2])))
    except:
        return []

# ==========================================
# 1. CONFIGURAÇÃO GLOBAL
# ==========================================
st.set_page_config(page_title="Lemos Lambda", page_icon="λ", layout="wide")

if 'alvos_val' not in st.session_state: st.session_state.alvos_val = ""
if 'fonte_val' not in st.session_state: st.session_state.fonte_val = ""
if 'alvo_val' not in st.session_state: st.session_state.alvo_val = ""

# ==========================================
# 2. DICIONÁRIO DE TRADUÇÃO (I18N)
# ==========================================
TEXTOS = {
    "pt": {
        "titulo_desk": "λ Lemos Lambda: Deep Science",
        "subtitulo": "**Ferramenta de Prospecção de Alto Impacto**",
        "credenciais": "1. Credenciais",
        "email_label": "Seu E-mail:",
        "periodo": "📅 Período:",
        "config": "2. Configuração (Órgãos)",
        "label_fonte": "**Fonte (Órgão, tecido, célula similar):**",
        "label_alvo": "**Alvo (Órgão de interesse):**",
        "sec_alvos": "3. Palavras-chave",
        "label_lista": "**Palavras-chave de Pesquisa:**",
        "btn_restaurar": "📥 Restaurar Padrão",
        "btn_trend": "🔍 Injetar Tendências (2025)",
        "btn_avanco": "🚀 Rumo ao Avanço",
        "col_ratio": "Ratio",
        "raio_x": "🔎 Raio-X",
        "btn_ler": "Ler Artigos",
        "lendo": "Buscando e Traduzindo...",
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
        "btn_restaurar": "📥 Restore Default",
        "btn_trend": "🔍 Inject Trends (2025)",
        "btn_avanco": "🚀 Launch Analysis",
        "col_ratio": "Ratio",
        "raio_x": "🔎 X-Ray",
        "btn_ler": "Read Papers",
        "lendo": "Searching and Translating...",
        "baixar": "📥 Download CSV"
    }
}

# ==========================================
# 3. FUNÇÕES DE SUPORTE
# ==========================================
def carregar_setup_lemos(t):
    st.session_state.alvos_val = "Piezo1, Piezo2, TREK-1, TASK-1, GPR35, GPR55, TAS2R14, P2X4"
    st.session_state.fonte_val = "Brain OR Kidney OR Liver OR Lung"
    st.session_state.alvo_val = "Bladder OR Vesical OR Urothelium"
    st.toast("Setup Carregado!", icon="🧬")

def carregar_alvos_apenas(t): 
    st.session_state.alvos_val = "Piezo1, Piezo2, TREK-1, TASK-1, GPR35, GPR55, TAS2R14, P2X4"
    st.toast("Lista Restaurada!", icon="✨")

def consultar_pubmed_count(termo_farmaco, termo_orgao, email, y_start, y_end):
    Entrez.email = email
    query = f"({termo_farmaco}) AND ({termo_orgao}) AND {y_start}:{y_end}[DP]" if termo_orgao else f"({termo_farmaco}) AND {y_start}:{y_end}[DP]"
    try: return int(Entrez.read(Entrez.esearch(db="pubmed", term=query, retmax=0))["Count"])
    except: return 0

def buscar_resumos_detalhados(termo_farmaco, termo_orgao, email, y_start, y_end, lang_target, limit=3):
    Entrez.email = email
    query = f"({termo_farmaco}) AND ({termo_orgao}) AND {y_start}:{y_end}[DP]"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=limit, sort="relevance")
        ids = Entrez.read(handle)["IdList"]
        if not ids: return []
        records = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text").read().split("\n\n")
        artigos = []
        for art_text in records:
            art_data = {"PMID": "N/A", "Title": "S/T", "Source": "N/A", "Abstract": ""}
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

st.title(t["titulo_desk"])
st.markdown(t["subtitulo"])

st.sidebar.header(t["credenciais"])
email_user = st.sidebar.text_input(t["email_label"], placeholder="email@exemplo.com", key="email_desk")
anos = st.sidebar.slider(t["periodo"], 1990, 2025, (2010, 2025), key="anos_desk")

st.sidebar.header(t["config"])
st.sidebar.markdown(t["label_fonte"])
t_fonte = st.sidebar.text_input("F", key="fonte_val", label_visibility="collapsed")
st.sidebar.markdown(t["label_alvo"])
t_alvo = st.sidebar.text_input("A", key="alvo_val", label_visibility="collapsed")

st.sidebar.button("🎓 Setup Guilherme Lemos", on_click=carregar_setup_lemos, args=(t,))
st.sidebar.markdown("---")

st.sidebar.header(t["sec_alvos"])
alvos_in = st.sidebar.text_area(t["label_lista"], key="alvos_val", height=150)

# Botão de Tendências abaixo das Palavras-chave
if st.sidebar.button(t["btn_trend"]):
    if email_user and "@" in email_user:
        with st.sidebar.status("Minerando Tendências..."):
            novos = buscar_alvos_emergentes(email_user)
            if novos:
                txt = ", ".join(novos)
                st.session_state.alvos_val = (st.session_state.alvos_val + ", " + txt).strip(", ")
                st.sidebar.success("Novos alvos injetados!"); time.sleep(1); st.rerun()
    else: st.sidebar.error("E-mail necessário.")

st.sidebar.button(t["btn_restaurar"], on_click=carregar_alvos_apenas, args=(t,))

if st.sidebar.button(t["btn_avanco"], type="primary"):
    if not email_user: st.error("E-mail obrigatório!")
    else:
        lst = sorted(list(set([x.strip() for x in alvos_in.split(",") if x.strip()])))
        res = []; bar = st.progress(0)
        for i, item in enumerate(lst):
            nf = consultar_pubmed_count(item, t_fonte, email_user, anos[0], anos[1]) if t_fonte else 0
            na = consultar_pubmed_count(item, t_alvo, email_user, anos[0], anos[1]) if t_alvo else 0
            
            # Lógica do Ratio e Status
            ratio = nf/na if na > 0 else (nf if nf > 0 else 0)
            status = "💎 DIAMANTE" if ratio > 10 and nf > 50 else "🥇 Ouro" if ratio > 2 else "📉 Raro"
            
            res.append({"Alvo": item, "Status": status, "Ratio": ratio, "Qtd_Fonte": nf, "Qtd_Alvo": na})
            bar.progress((i+1)/len(lst))
        st.session_state['dados_res'] = pd.DataFrame(res).sort_values(by="Ratio", ascending=False); st.rerun()

if 'dados_res' in st.session_state:
    df = st.session_state['dados_res']
    st.plotly_chart(px.bar(df.head(20), x="Alvo", y="Ratio", color="Status", color_discrete_map={"💎 DIAMANTE": "#00CC96", "🥇 Ouro": "#636EFA", "📉 Raro": "#EF553B"}), use_container_width=True)
    st.dataframe(df, use_container_width=True, hide_index=True)
    st.download_button(t["baixar"], df.to_csv(index=False).encode('utf-8'), "lemos_lambda_results.csv")
    
    st.divider()
    st.header(t["raio_x"])
    sel = st.selectbox("Escolha um alvo para ler resumos:", df['Alvo'].unique())
    if st.button(t["btn_ler"]):
        with st.spinner(t["lendo"]):
            arts = buscar_resumos_detalhados(sel, t_alvo, email_user, anos[0], anos[1], lang)
            if not arts: st.info("Sem artigos encontrados.")
            for a in arts:
                with st.expander(a['Title']):
                    st.success(a['Resumo_IA'])
                    st.markdown(f"[PubMed Link](https://pubmed.ncbi.nlm.nih.gov/{a['PMID']})")
