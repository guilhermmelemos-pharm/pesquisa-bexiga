"""
Lemos Lambda: Deep Science Prospector
Copyright (c) 2025 Guilherme Lemos
Licensed under the MIT License.
Version: 2.0 (Stable)
"""
import streamlit as st

# --- 1. CONFIGURAÇÃO DA PÁGINA ---
st.set_page_config(
    page_title="λ Lemos Lambda v2.0: Deep Science Prospector", 
    page_icon="λ", 
    layout="wide", 
    initial_sidebar_state="collapsed"
)

import pandas as pd
import plotly.express as px
from datetime import datetime
import time
import scipy.stats as stats
import constantes as c
import backend as bk 

# --- 2. ESTADO E INICIALIZAÇÃO ---
defaults = {
    'pagina': 'home', 'alvos_val': "", 'resultado_df': None, 'news_index': 0,
    'input_alvo': "", 'input_fonte': "", 'input_email': "", 'artigos_detalhe': None,
    'email_guardado': "", 'alvo_guardado': "", 'lang': 'pt',
    'api_key_usuario': "", 'usar_ia_faxina': True, 'ia_global_switch': True 
}

for k, v in defaults.items():
    if k not in st.session_state: 
        st.session_state[k] = v

def get_textos(): 
    return c.TEXTOS.get(st.session_state.lang, c.TEXTOS["pt"])

t = get_textos()

# --- 3. CSS ---
st.markdown("""
<style>
.stButton button { border-radius: 12px; height: 50px; font-weight: bold; }
div[data-testid="stMetricValue"] { font-size: 1.8rem !important; }
.big-button button { background-color: #FF4B4B !important; color: white !important; border: none; font-size: 1.1rem !important; }
.stTextArea textarea { font-family: monospace; }
.header-style { font-size: 2.5rem; font-weight: 700; color: #FAFAFA; margin-bottom: 0px; }
.sub-header-style { font-size: 1.2rem; font-weight: 400; color: #A0A0A0; margin-bottom: 20px; }
</style>
""", unsafe_allow_html=True)

# ================= LÓGICA =================

def mudar_idioma(novo_lang): 
    st.session_state.lang = novo_lang
    resetar_pesquisa()

def limpar_campo(k): 
    st.session_state[k] = ""

def limpar_lista_total(): 
    st.session_state.alvos_val = ""

def adicionar_termos_seguro(lista, textos):
    atuais = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
    atuais_up = [x.upper() for x in atuais]
    blacklist = [b.lower() for b in c.BLACKLIST_GERAL]
    add = []
    for item in lista:
        cl = item.strip()
        if not cl or any(b in cl.lower() for b in blacklist): 
            continue
        if cl.upper() not in atuais_up:
            atuais.append(cl)
            atuais_up.append(cl.upper())
            add.append(cl)
    st.session_state.alvos_val = ", ".join(atuais)
    return len(add)

# ---------- (demais funções permanecem idênticas) ----------

# ====== TRECHO CRÍTICO CORRIGIDO ======

if st.session_state.artigos_detalhe:
    st.info(f"{len(st.session_state.artigos_detalhe)} papers found for {sel}.")
    for i, art in enumerate(st.session_state.artigos_detalhe):
        with st.expander(f"📄 {art['Title']}", expanded=False):
            st.caption(f"**Keywords/Context:** {art.get('Info_IA', 'N/A')[:200]}...")
            c_ia, c_link = st.columns([1, 1])
            with c_ia:
                if not st.session_state.api_key_usuario:
                    nova_chave = st.text_input(
                        "Google API Key:", 
                        type="password", 
                        key=f"key_input_{i}"
                    )
                    if nova_chave:
                        st.session_state.api_key_usuario = nova_chave
                        st.rerun()

                if st.session_state.api_key_usuario and st.session_state.ia_global_switch:
                    if st.button(f"🤖 {t['btn_investigar']}", key=f"btn_ia_{i}"):
                        with st.spinner("Analyzing..."):
                            resumo = bk.analisar_abstract_com_ia(
                                art['Title'],
                                art.get('Info_IA', ''),
                                st.session_state.api_key_usuario
                            )
                            st.markdown(
                                f"""
                                <div style='background-color:#262730;
                                            color:#ffffff;
                                            padding:15px;
                                            border-radius:8px;
                                            border-left:5px solid #FF4B4B;
                                            margin-top:10px;'>
                                <small style='color:#FF4B4B;'>🧠 <b>Lemos Lambda AI:</b></small><br>
                                <span style='font-size:1.1em;'>{resumo}</span>
                                </div>
                                """,
                                unsafe_allow_html=True
                            )
            with c_link:
                st.link_button(t["btn_pubmed"], art['Link'], use_container_width=True)
