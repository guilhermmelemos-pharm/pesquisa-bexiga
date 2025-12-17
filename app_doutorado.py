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
# 1. CONFIGURAÇÃO GLOBAL
# ==========================================
st.set_page_config(page_title="Lemos Doutorado", page_icon="🎓", layout="wide")

# CSS: Padroniza imagens e garante largura total dos botões
st.markdown("""
    <style>
    div[data-testid="stImage"] img {
        height: 150px !important;
        object-fit: cover !important;
        border-radius: 8px !important;
    }
    .stButton button {
        width: 100%;
    }
    </style>
""", unsafe_allow_html=True)

# Inicialização do Session State
if 'alvos_val' not in st.session_state: st.session_state.alvos_val = ""
if 'fonte_val' not in st.session_state: st.session_state.fonte_val = ""
if 'alvo_val' not in st.session_state: st.session_state.alvo_val = ""
if 'news_index' not in st.session_state: st.session_state.news_index = 0

# ==========================================
# 2. FUNÇÕES DE SUPORTE
# ==========================================
@st.cache_data(ttl=3600)
def buscar_todas_noticias():
    feeds = [
        {"url": "https://www.sciencedaily.com/rss/health_medicine/pharmacology.xml", "lang": "🇺🇸"},
        {"url": "https://www.fiercebiotech.com/rss/biotech", "lang": "🇺🇸"},
        {"url": "https://www.nature.com/nbt.rss", "lang": "🇬🇧"},
        {"url": "https://agencia.fapesp.br/rss/", "lang": "🇧🇷"},
    ]
    noticias = []
    backups = ["https://images.unsplash.com/photo-1532094349884-543bc11b234d?w=400&h=250&fit=crop"]
    
    for fonte in feeds:
        try:
            feed = feedparser.parse(fonte["url"])
            for entry in feed.entries[:3]:
                img = random.choice(backups)
                if 'media_content' in entry: img = entry.media_content[0]['url']
                noticias.append({"titulo": entry.title, "link": entry.link, "img": img, "fonte": fonte["lang"]})
        except: continue
    random.shuffle(noticias)
    return noticias

@st.fragment(run_every=60) 
def exibir_radar_cientifico():
    news = buscar_todas_noticias()
    if not news: return
    idx = st.session_state.news_index % len(news)
    batch = news[idx:idx+3]
    st.session_state.news_index += 3
    
    with st.container(border=True):
        st.caption("📡 **Radar Científico**")
        cols = st.columns(3)
        for i, n in enumerate(batch):
            with cols[i]:
                st.image(n['img'], use_container_width=True)
                st.caption(f"{n['fonte']} {n['titulo'][:60]}...")
                st.link_button("Ler", n['link'], use_container_width=True)

# BANCO DE DADOS
CANDIDATOS_MINERACAO = ["GPR37", "GPR17", "GPR55", "Piezo1", "Piezo2", "TMEM16A", "Ferroptosis", "Pyroptosis", "Gasdermin D", "Sirtuin-1"]
SUGESTOES_ALVOS_RAW = """
-- ALVOS EXEMPLO --
Piezo1, Piezo2, Ferroptosis, GPX4, P2X3, P2X7, TRPV4, NLRP3, Gasdermin D, Clock genes
"""
LISTA_ALVOS_PRONTA = ", ".join([x.strip() for x in SUGESTOES_ALVOS_RAW.split(',') if not x.strip().startswith("--")])

PRESETS_ORGAOS = {"(Sugestão Lemos)": {"fonte": "Kidney OR Liver", "alvo": "Bladder OR Urothelium"}}

def carregar_setup_lemos():
    st.session_state.alvos_val = LISTA_ALVOS_PRONTA
    st.session_state.fonte_val = PRESETS_ORGAOS["(Sugestão Lemos)"]["fonte"]
    st.session_state.alvo_val = PRESETS_ORGAOS["(Sugestão Lemos)"]["alvo"]

def carregar_alvos_apenas(): st.session_state.alvos_val = LISTA_ALVOS_PRONTA
def limpar_campo_fonte(): st.session_state.fonte_val = ""
def limpar_campo_alvo(): st.session_state.alvo_val = ""
def limpar_campo_alvos(): st.session_state.alvos_val = ""

def minerar_blue_oceans(orgao, email):
    if not orgao or not email:
        st.toast("Preencha Alvo e E-mail!", icon="⚠️"); return
    
    Entrez.email = email
    st.toast("⛏️ Mineração iniciada... Clique em 'Rumo ao Avanço' ao terminar.", icon="⏳")
    encontrados = []
    my_bar = st.progress(0, text="Procurando termos chave...")
    
    for i, termo in enumerate(CANDIDATOS_MINERACAO):
        try:
            handle = Entrez.esearch(db="pubmed", term=f"({termo}) AND ({orgao})", retmax=0)
            count = int(Entrez.read(handle)["Count"])
            if 0 <= count < 150: encontrados.append(termo)
            my_bar.progress((i+1)/len(CANDIDATOS_MINERACAO), text=f"Testando: {termo} ({count})")
            time.sleep(0.1)
        except: continue
    
    my_bar.empty()
    if encontrados:
        st.session_state.alvos_val = ", ".join(encontrados)
        st.toast(f"✅ {len(encontrados)} termos raros encontrados! Clique em 'Rumo ao Avanço'.", icon="🚀")

def processar_upload():
    uploaded = st.session_state.get('uploader_key')
    if uploaded:
        try:
            content = uploaded.getvalue().decode("utf-8")
            st.session_state.alvos_val = " ".join(content.replace("\n", ",").split())
            st.toast("Arquivo importado!", icon="📂")
        except: st.error("Erro ao ler arquivo.")

def consultar_pubmed_count(termo, orgao, email, y1, y2):
    if not email: return -1
    Entrez.email = email
    q = f"({termo}) AND ({orgao}) AND {y1}:{y2}[DP]" if orgao else f"({termo}) AND {y1}:{y2}[DP]"
    try: return int(Entrez.read(Entrez.esearch(db="pubmed", term=q, retmax=0))["Count"])
    except: return -1

def buscar_resumos(termo, orgao, email, y1, y2):
    # Função simplificada para leitura
    return [] # (Mantido a estrutura, implementar se necessário a leitura completa)

# ==========================================
# 3. INTERFACE
# ==========================================
modo = st.sidebar.radio("📱 Modo:", ["Desktop (Completo)", "Mobile (Pocket)"], index=0)
st.sidebar.markdown("---")

if modo == "Desktop (Completo)":
    st.title("🎓 Lemos Doutorado: Deep Science")
    if 'dados_desk' not in st.session_state: exibir_radar_cientifico()
    
    st.sidebar.header("1. Credenciais")
    email_user = st.sidebar.text_input("Seu E-mail:", key="email_desk")
    anos = st.sidebar.slider("📅 Período:", 1990, 2025, (2010, 2025))
    
    st.sidebar.markdown("---")
    st.sidebar.header("2. Configuração")
    
    # --- ALINHAMENTO PERFEITO (vertical_alignment="bottom") ---
    st.sidebar.markdown("**Fonte (Orgão, tecido, célula similar):**")
    c1, c2 = st.sidebar.columns([6, 1], vertical_alignment="bottom")
    with c1: t_fonte = st.text_input("Fonte", key="fonte_val", label_visibility="collapsed")
    with c2: st.button("🗑️", key="del_f", on_click=limpar_campo_fonte)

    st.sidebar.markdown("**Alvo (Orgão de interesse):**")
    c3, c4 = st.sidebar.columns([6, 1], vertical_alignment="bottom")
    with c3: t_alvo = st.text_input("Alvo", key="alvo_val", label_visibility="collapsed")
    with c4: st.button("🗑️", key="del_a", on_click=limpar_campo_alvo)
    
    st.sidebar.button("🎓 Doutorado Guilherme Lemos", on_click=carregar_setup_lemos)
    
    st.sidebar.markdown("---")
    st.sidebar.header("3. Alvos")
    
    with st.sidebar.expander("📂 Importar Biblioteca"):
        st.file_uploader("Upload", key="uploader_key", on_change=processar_upload)
    
    st.sidebar.markdown("**Lista de Pesquisa:**")
    c5, c6 = st.sidebar.columns([6, 1], vertical_alignment="bottom")
    with c5: alvos_in = st.text_area("Lista", key="alvos_val", height=150, label_visibility="collapsed")
    with c6: st.button("🗑️", key="del_l", on_click=limpar_campo_alvos)

    b1, b2 = st.sidebar.columns(2)
    b1.button("📥 Restaurar Padrão", on_click=carregar_alvos_apenas)
    b2.button("⛏️ Minerar 'Blue Oceans'", on_click=minerar_blue_oceans, args=(t_alvo, email_user))
    
    st.sidebar.markdown("---")

    if st.sidebar.button("🚀 Rumo ao Avanço", type="primary"):
        if not email_user: st.error("E-mail obrigatório!")
        elif not alvos_in: st.warning("Lista vazia!")
        else:
            lst = [x.strip() for x in alvos_in.split(",") if x.strip()]
            res = []
            pg = st.progress(0)
            for i, item in enumerate(lst):
                nf = consultar_pubmed_count(item, t_fonte, email_user, anos[0], anos[1]) if t_fonte else 0
                na = consultar_pubmed_count(item, t_alvo, email_user, anos[0], anos[1]) if t_alvo else 0
                ng = consultar_pubmed_count(item, "", email_user, anos[0], anos[1]) if not t_fonte and not t_alvo else 0
                
                pot = 0
                stat = "N/A"
                
                if t_fonte and t_alvo:
                    pot = nf/na if na > 0 else nf
                    stat = "💎 DIAMANTE" if pot > 10 and nf > 50 else "🔴 Saturado" if na >= nf else "🥇 Ouro"
                elif t_alvo:
                    pot = na
                    stat = "🔥 Hot" if na > 200 else "📉 Raro"
                else:
                    pot = ng
                    stat = "Global"

                res.append({
                    "Alvo": item, 
                    "Status": stat, 
                    "Potencial": pot, 
                    "Qtd_Fonte": nf, 
                    "Qtd_Alvo": na if t_alvo else ng
                })
                pg.progress((i+1)/len(lst))
            
            pg.empty()
            st.session_state['dados_desk'] = pd.DataFrame(res).sort_values(by="Potencial", ascending=False)
            st.rerun()

    if 'dados_desk' in st.session_state:
        df = st.session_state['dados_desk']
        
        # --- TABELA LIMPA E RENOMEADA ---
        # Renomeia as colunas dinamicamente baseada no que o usuário digitou
        nome_fonte = f"Artigos em: {t_fonte}" if t_fonte else "Sem Fonte"
        nome_alvo = f"Artigos em: {t_alvo}" if t_alvo else "Artigos (Global)"
        
        df_show = df.rename(columns={
            "Qtd_Fonte": nome_fonte,
            "Qtd_Alvo": nome_alvo,
            "Potencial": "Potencial (Ratio/Qtd)"
        })
        
        # Configuração das Colunas (Formatação Numérica)
        st.dataframe(
            df_show.style.format({
                "Potencial (Ratio/Qtd)": "{:.1f}", # Apenas 1 casa decimal
                nome_fonte: "{:.0f}", # Inteiro
                nome_alvo: "{:.0f}"   # Inteiro
            }),
            use_container_width=True,
            height=500,
            column_config={
                "Status": st.column_config.TextColumn("Status", help="Classificação de oportunidade"),
            }
        )
        
        csv = df_show.to_csv(index=False).encode('utf-8')
        st.download_button("📥 Baixar Planilha", csv, "resultado_lemos.csv", "text/csv")

elif modo == "Mobile (Pocket)":
    st.title("📱 Lemos Pocket")
    # ... (Mesma lógica simplificada para Mobile, usando vertical_alignment="bottom" nas colunas)
