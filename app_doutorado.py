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
# 1. CONFIGURAÇÃO GLOBAL E ESTADO
# ==========================================
st.set_page_config(page_title="Lemos Lambda", page_icon="λ", layout="wide")

# Inicialização do Session State (CRÍTICO para evitar erros de API)
if 'alvos_val' not in st.session_state: st.session_state.alvos_val = ""
if 'fonte_val' not in st.session_state: st.session_state.fonte_val = ""
if 'alvo_val' not in st.session_state: st.session_state.alvo_val = ""

# CSS para melhorar o visual
st.markdown("""
    <style>
    .stButton button { width: 100%; }
    .stTextArea textarea { font-family: monospace; }
    </style>
""", unsafe_allow_html=True)

# ==========================================
# 2. DICIONÁRIO DE TRADUÇÃO (I18N)
# ==========================================
TEXTOS = {
    "pt": {
        "titulo": "λ Lemos Lambda: Deep Science",
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
        "raio_x": "🔎 Raio-X",
        "btn_ler": "Ler Artigos",
        "lendo": "Buscando e Traduzindo...",
        "baixar": "📥 Baixar Planilha"
    },
    "en": {
        "titulo": "λ Lemos Lambda: Deep Science",
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
        "raio_x": "🔎 X-Ray",
        "btn_ler": "Read Papers",
        "lendo": "Searching and Translating...",
        "baixar": "📥 Download CSV"
    }
}

# ==========================================
# 3. FUNÇÕES DE SUPORTE
# ==========================================
def carregar_setup_lemos():
    st.session_state.alvos_val = "Piezo1, Piezo2, TREK-1, TASK-1, GPR35, GPR55, TAS2R14, P2X4"
    st.session_state.fonte_val = "Brain OR Kidney OR Liver OR Lung"
    st.session_state.alvo_val = "Bladder OR Vesical OR Urothelium"

def carregar_alvos_apenas(): 
    st.session_state.alvos_val = "Piezo1, Piezo2, TREK-1, TASK-1, GPR35, GPR55, TAS2R14, P2X4"

def consultar_pubmed_count(termo_farmaco, termo_orgao, email, y_start, y_end):
    if not email: return 0
    Entrez.email = email
    query = f"({termo_farmaco}) AND ({termo_orgao}) AND {y_start}:{y_end}[DP]" if termo_orgao else f"({termo_farmaco}) AND {y_start}:{y_end}[DP]"
    try: 
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        return int(Entrez.read(handle)["Count"])
    except: return 0

def buscar_resumos_detalhados(termo_farmaco, termo_orgao, email, y_start, y_end, lang_target):
    Entrez.email = email
    query = f"({termo_farmaco}) AND ({termo_orgao}) AND {y_start}:{y_end}[DP]"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=3, sort="relevance")
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
# 4. INTERFACE PRINCIPAL
# ==========================================
lang_opt = st.sidebar.radio("Language:", ["🇧🇷 Português", "🇺🇸 English"])
lang = "pt" if "Português" in lang_opt else "en"
t = TEXTOS[lang]

st.title(t["titulo"])
st.markdown(t["subtitulo"])

# --- SIDEBAR: PASSO 1 E 2 ---
st.sidebar.header(t["credenciais"])
email_user = st.sidebar.text_input(t["email_label"], value="", placeholder="seu@email.com")
anos = st.sidebar.slider(t["periodo"], 1990, 2025, (2010, 2025))

st.sidebar.header(t["config"])
st.sidebar.markdown(t["label_fonte"])
st.session_state.fonte_val = st.sidebar.text_input("F", value=st.session_state.fonte_val, label_visibility="collapsed")
st.sidebar.markdown(t["label_alvo"])
st.session_state.alvo_val = st.sidebar.text_input("A", value=st.session_state.alvo_val, label_visibility="collapsed")

if st.sidebar.button("🎓 Setup Guilherme Lemos"):
    carregar_setup_lemos()
    st.rerun()

st.sidebar.markdown("---")

# --- SIDEBAR: PASSO 3 (PALAVRAS-CHAVE) ---
st.sidebar.header(t["sec_alvos"])

# IMPORTANTE: Usamos 'value' em vez de 'key' para permitir injeção via código sem erro
alvos_input = st.sidebar.text_area(t["label_lista"], value=st.session_state.alvos_val, height=150)
st.session_state.alvos_val = alvos_input # Mantém o estado sincronizado

# Botão de Tendências abaixo da caixa
if st.sidebar.button(t["btn_trend"]):
    if email_user and "@" in email_user:
        with st.sidebar.status("Minerando PubMed 2025..."):
            novos = buscar_alvos_emergentes(email_user)
            if novos:
                txt_novos = ", ".join(novos)
                # Concatena com o que já existe
                if st.session_state.alvos_val:
                    st.session_state.alvos_val = (st.session_state.alvos_val.strip(", ") + ", " + txt_novos)
                else:
                    st.session_state.alvos_val = txt_novos
                st.sidebar.success("Novidades injetadas!")
                time.sleep(1)
                st.rerun()
    else:
        st.sidebar.error("E-mail válido necessário para minerar.")

if st.sidebar.button(t["btn_restaurar"]):
    carregar_alvos_apenas()
    st.rerun()

# --- EXECUÇÃO ---
if st.sidebar.button(t["btn_avanco"], type="primary"):
    if not email_user:
        st.error("Por favor, insira seu e-mail no Passo 1.")
    elif not st.session_state.alvos_val:
        st.warning("A lista de palavras-chave está vazia.")
    else:
        lista_alvos = sorted(list(set([x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()])))
        resultados = []
        prog_bar = st.progress(0)
        
        for i, alvo in enumerate(lista_alvos):
            st.toast(f"Analisando: {alvo}")
            count_fonte = consultar_pubmed_count(alvo, st.session_state.fonte_val, email_user, anos[0], anos[1])
            count_alvo = consultar_pubmed_count(alvo, st.session_state.alvo_val, email_user, anos[0], anos[1])
            
            # Cálculo do Ratio (Se alvo for 0, ratio é o total da fonte como 'Oceano Azul')
            ratio = count_fonte / count_alvo if count_alvo > 0 else count_fonte
            
            # Lógica de Status
            if ratio > 10 and count_fonte > 50: status = "💎 DIAMANTE"
            elif ratio > 2: status = "🥇 Ouro"
            else: status = "📉 Raro"
            
            resultados.append({
                "Alvo": alvo,
                "Status": status,
                "Ratio": round(ratio, 2),
                "Qtd_Fonte": count_fonte,
                "Qtd_Alvo": count_alvo
            })
            prog_bar.progress((i + 1) / len(lista_alvos))
        
        st.session_state.dados_final = pd.DataFrame(resultados).sort_values(by="Ratio", ascending=False)
        st.rerun()

# --- EXIBIÇÃO DE RESULTADOS ---
if 'dados_final' in st.session_state:
    df = st.session_state.dados_final
    
    st.subheader("Análise de Potencial Prospectivo")
    fig = px.bar(
        df.head(20), 
        x="Alvo", y="Ratio", 
        color="Status",
        color_discrete_map={"💎 DIAMANTE": "#00CC96", "🥇 Ouro": "#636EFA", "📉 Raro": "#EF553B"},
        text_auto=True
    )
    st.plotly_chart(fig, use_container_width=True)
    
    st.dataframe(df, use_container_width=True, hide_index=True)
    
    st.download_button(
        t["baixar"], 
        df.to_csv(index=False).encode('utf-8'), 
        "analise_lemos_lambda.csv", 
        "text/csv"
    )
    
    # --- RAIO-X (ARTIGOS) ---
    st.divider()
    st.header(t["raio_x"])
    escolha = st.selectbox("Selecione um alvo para ler resumos recentes:", df['Alvo'].unique())
    
    if st.button(t["btn_ler"]):
        with st.spinner(t["lendo"]):
            artigos = buscar_resumos_detalhados(escolha, st.session_state.alvo_val, email_user, anos[0], anos[1], lang)
            if not artigos:
                st.info("Nenhum artigo detalhado encontrado para este cruzamento.")
            else:
                for a in artigos:
                    with st.expander(f"📄 {a['Title']}"):
                        st.success(a['Resumo_IA'])
                        st.markdown(f"[Abrir no PubMed](https://pubmed.ncbi.nlm.nih.gov/{a['PMID']})")
