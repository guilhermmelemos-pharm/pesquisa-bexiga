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
import plotly.express as px
from datetime import datetime
import time

# ==========================================
# IMPORTAÇÃO MODULAR
# ==========================================
import constantes as c
import backend as bk

# ==========================================
# CONFIGURAÇÃO VISUAL
# ==========================================
st.set_page_config(page_title="Lemos Lambda", page_icon="λ", layout="wide")

st.markdown("""
    <style>
    .stButton button { border-radius: 12px; height: 50px; font-weight: bold; }
    .css-1d391kg { padding-top: 2rem; }
    div[data-testid="stMetricValue"] { font-size: 1.8rem !important; }
    div[data-testid="stImage"] img { height: 150px !important; object-fit: cover !important; border-radius: 8px !important; }
    </style>
""", unsafe_allow_html=True)

# --- ESTADO ---
if 'alvos_val' not in st.session_state: st.session_state.alvos_val = ""
if 'resultado_df' not in st.session_state: st.session_state.resultado_df = None
if 'news_index' not in st.session_state: st.session_state.news_index = 0

# --- IDIOMA ---
lang_opt = st.sidebar.radio("🌐", ["🇧🇷 PT", "🇺🇸 EN"], horizontal=True)
lang = "pt" if "PT" in lang_opt else "en"
t = c.TEXTOS[lang]

# ==========================================
# UI PRINCIPAL (WIZARD)
# ==========================================
st.title(t["titulo_desk"])
st.caption(t["subtitulo"])

# --- ÁREA DE INPUT (PASSO 1) ---
with st.container(border=True):
    col_input, col_btn = st.columns([3, 1])
    
    with col_input:
        st.subheader(t["step_1"])
        email_user = st.text_input("E-mail (Obrigatório para PubMed)", placeholder="ex: gl@unifesp.br")
        alvo = st.text_input(t["label_alvo"], placeholder=t["holder_alvo"])
        contexto = st.text_input(t["label_fonte"], placeholder=t["holder_fonte"])
    
    with col_btn:
        st.write(" ")
        st.write(" ")
        st.write(" ") # Espaçamento visual
        # BOTÃO INTELIGENTE DE MINERAÇÃO
        if st.button(t["btn_magic"], type="primary"):
            if not email_user or not alvo:
                st.error("⚠️ E-mail e Alvo necessários!")
            else:
                with st.status(t["prog_magic"], expanded=True) as status:
                    st.write(t["status_minerando"])
                    # Chama o backend para buscar termos novos (Blue Ocean)
                    novos_termos = bk.buscar_alvos_emergentes_pubmed(alvo, email_user)
                    
                    st.write(t["status_filtrando"])
                    termos_base = c.CANDIDATOS_MINERACAO
                    
                    # Junta lista base + descobertas novas + ácidos/moléculas complexas
                    lista_final = list(set(termos_base + novos_termos))
                    st.session_state.alvos_val = ", ".join(lista_final)
                    
                    status.update(label=t["status_pronto"], state="complete", expanded=False)
                    st.toast(f"✅ {len(novos_termos)} novos termos complexos encontrados!", icon="🧬")

# --- ÁREA DE ANÁLISE (PASSO 2) ---
st.divider()
st.subheader(t["step_2"])

# Só mostra a área de análise se tiver termos carregados
if st.session_state.alvos_val:
    with st.expander("📝 Ver/Editar Lista de Palavras-Chave (Clique para abrir)", expanded=False):
        st.session_state.alvos_val = st.text_area("Termos", value=st.session_state.alvos_val, height=150)

    # BOTÃO DE CÁLCULO (RATIO)
    if st.button(t["analise_btn"], use_container_width=True):
        if not email_user: st.error("E-mail necessário.")
        else:
            lista = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
            resultados = []
            
            progresso = st.progress(0)
            status_text = st.empty()
            
            for i, item in enumerate(lista):
                status_text.caption(f"🔍 Investigando: **{item}**...")
                time.sleep(0.05) # Pequeno delay para UI fluida
                
                # Buscas no Backend
                n_global = bk.consultar_pubmed_count(item, contexto, email_user, 2015, 2025)
                n_especifico = bk.consultar_pubmed_count(item, alvo, email_user, 2015, 2025)
                
                # Lógica de Classificação Blue Ocean
                ratio = n_global / n_especifico if n_especifico > 0 else n_global
                
                if n_especifico == 0 and n_global > 50: tag = "💎 Blue Ocean (Inexplorado)"
                elif ratio > 10: tag = "🥇 Ouro (Promissor)"
                elif ratio < 2: tag = "🔴 Saturado"
                else: tag = "⚖️ Neutro"
                
                resultados.append({
                    "Molécula/Alvo": item, 
                    "Status": tag, 
                    "Potencial (Ratio)": round(ratio, 1),
                    "Artigos no Alvo": n_especifico,
                    "Global": n_global
                })
                progresso.progress((i+1)/len(lista))
            
            progresso.empty()
            status_text.empty()
            st.session_state.resultado_df = pd.DataFrame(resultados).sort_values(by="Potencial (Ratio)", ascending=False)
            st.rerun()

# --- EXIBIÇÃO DE RESULTADOS ---
if st.session_state.resultado_df is not None:
    df = st.session_state.resultado_df
    
    st.markdown("### 🎯 Resultados")
    
    # KPIs Rápidos (Destaques)
    if not df.empty:
        top_term = df.iloc[0]
        c1, c2, c3 = st.columns(3)
        c1.metric("🏆 Maior Potencial", top_term['Molécula/Alvo'])
        c2.metric("📊 Score", top_term['Potencial (Ratio)'])
        c3.metric("📚 Artigos Existentes", top_term['Artigos no Alvo'])
        
        # Gráfico Interativo
        fig = px.bar(df.head(20), x="Molécula/Alvo", y="Potencial (Ratio)", color="Status",
                     color_discrete_map={"💎 Blue Ocean (Inexplorado)": "#00CC96", "🥇 Ouro (Promissor)": "#636EFA", "🔴 Saturado": "#EF553B"})
        st.plotly_chart(fig, use_container_width=True)
        
        # Tabela Detalhada
        st.dataframe(df, use_container_width=True, hide_index=True)
    else:
        st.warning(t["tabela_vazia"])

# --- RODAPÉ ---
st.markdown("---")
st.caption(f"© 2025 Guilherme Lemos | {t['footer_citar']}")