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
    div[data-testid="stMetricValue"] { font-size: 1.8rem !important; }
    div[data-testid="stImage"] img { height: 150px !important; object-fit: cover !important; border-radius: 8px !important; }
    </style>
""", unsafe_allow_html=True)

# --- ESTADO ---
if 'alvos_val' not in st.session_state: st.session_state.alvos_val = ""
if 'resultado_df' not in st.session_state: st.session_state.resultado_df = None
if 'news_index' not in st.session_state: st.session_state.news_index = 0
if 'fonte_val' not in st.session_state: st.session_state.fonte_val = ""

# --- IDIOMA ---
lang_opt = st.sidebar.radio("🌐", ["🇧🇷 PT", "🇺🇸 EN"], horizontal=True)
lang = "pt" if "PT" in lang_opt else "en"
t = c.TEXTOS[lang]

# ==========================================
# UI HELPERS
# ==========================================
@st.fragment(run_every=60) 
def exibir_radar_cientifico(lang_code):
    """Exibe notícias científicas em um container moderno"""
    news_list = bk.buscar_todas_noticias(lang_code)
    if not news_list: return
    total_news = len(news_list)
    idx = st.session_state.news_index % total_news
    batch = news_list[idx:idx+3]
    st.session_state.news_index += 3
    with st.container(border=True):
        st.caption(f"📡 **Radar Científico** (Updates ao vivo)")
        cols = st.columns(3)
        for i, n in enumerate(batch):
            with cols[i]:
                st.image(n['img'], use_container_width=True)
                st.markdown(f"**{n['titulo'][:60]}...**")
                st.caption(f"{n['bandeira']} {n['fonte']}")
                st.link_button("Ler" if lang_code=='pt' else "Read", n['link'], use_container_width=True)

def processar_upload():
    uploaded_file = st.session_state.get('uploader_key')
    if uploaded_file is not None:
        try:
            content = uploaded_file.getvalue().decode("utf-8")
            # Limpa e formata a lista importada
            termos_importados = [x.strip() for x in content.replace("\n", ",").split(",") if x.strip()]
            
            # Adiciona aos existentes ou substitui
            if st.session_state.alvos_val:
                st.session_state.alvos_val += ", " + ", ".join(termos_importados)
            else:
                st.session_state.alvos_val = ", ".join(termos_importados)
                
            st.toast(f"✅ {len(termos_importados)} termos importados com sucesso!", icon="📂")
        except: st.error("Erro ao ler arquivo.")

# ==========================================
# UI PRINCIPAL
# ==========================================
st.title(t["titulo_desk"])
st.caption(t["subtitulo"])

# 1. RADAR CIENTÍFICO (FEED) - Restaurado no topo
exibir_radar_cientifico(lang)

st.divider()

# 2. ÁREA DE INPUT (WIZARD)
with st.container(border=True):
    col_input, col_config = st.columns([2, 1])
    
    with col_input:
        st.subheader(t["step_1"])
        email_user = st.text_input("E-mail (Obrigatório para PubMed)", placeholder="ex: gl@unifesp.br")
        alvo = st.text_input(t["label_alvo"], placeholder=t["holder_alvo"])
        
        # Botão Mágico (Mineração Automática)
        if st.button(t["btn_magic"], type="primary"):
            if not email_user or not alvo:
                st.error("⚠️ E-mail e Alvo necessários!")
            else:
                with st.status(t["prog_magic"], expanded=True) as status:
                    st.write(t["status_minerando"])
                    novos_termos = bk.buscar_alvos_emergentes_pubmed(alvo, email_user)
                    
                    st.write(t["status_filtrando"])
                    termos_base = c.CANDIDATOS_MINERACAO
                    lista_final = list(set(termos_base + novos_termos))
                    
                    # Preserva termos já digitados se houver
                    if st.session_state.alvos_val:
                        existentes = [x.strip() for x in st.session_state.alvos_val.split(",")]
                        lista_final = list(set(lista_final + existentes))
                        
                    st.session_state.alvos_val = ", ".join(lista_final)
                    status.update(label=t["status_pronto"], state="complete", expanded=False)
                    st.toast(f"✅ {len(novos_termos)} novos termos encontrados!", icon="🧬")

    with col_config:
        st.subheader("⚙️ Comparação & Fonte")
        # Restaurado: Comparação com Órgão Fonte
        contexto = st.text_input(t["label_fonte"], value=st.session_state.fonte_val, placeholder=t["holder_fonte"], help="Define o 'universo' de comparação (ex: Rim, Cérebro) para calcular o Ratio.")
        st.session_state.fonte_val = contexto
        
        st.markdown("---")
        # Restaurado: Upload de Lista Própria
        st.write("📂 **Importar Lista Própria**")
        st.file_uploader("Upload (.csv/.txt)", type=["csv", "txt"], key="uploader_key", on_change=processar_upload, label_visibility="collapsed")
        st.caption("Formato: Termos separados por vírgula ou um por linha.")

# 3. ÁREA DE ANÁLISE (PASSO 2)
st.divider()
st.subheader(t["step_2"])

if st.session_state.alvos_val:
    with st.expander("📝 Ver/Editar Lista de Palavras-Chave (Clique para abrir)", expanded=False):
        c1, c2 = st.columns([5,1])
        with c1:
            st.session_state.alvos_val = st.text_area("Termos", value=st.session_state.alvos_val, height=100, label_visibility="collapsed")
        with c2:
            if st.button("🗑️ Limpar"): st.session_state.alvos_val = ""
            st.caption(f"Total: {len(st.session_state.alvos_val.split(',')) if st.session_state.alvos_val else 0}")

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
                time.sleep(0.05) 
                
                # Buscas no Backend
                # Se o usuário não preencheu contexto, usa busca global (apenas o termo)
                termo_contexto = contexto if contexto else None
                
                n_global = bk.consultar_pubmed_count(item, termo_contexto, email_user, 2015, 2025)
                n_especifico = bk.consultar_pubmed_count(item, alvo, email_user, 2015, 2025)
                
                # Lógica de Classificação Blue Ocean
                ratio = n_global / n_especifico if n_especifico > 0 else n_global
                
                if n_especifico == 0 and n_global > 50: tag = "💎 Blue Ocean"
                elif ratio > 10: tag = "🥇 Ouro"
                elif ratio < 2: tag = "🔴 Saturado"
                else: tag = "⚖️ Neutro"
                
                resultados.append({
                    "Molécula/Alvo": item, 
                    "Status": tag, 
                    "Potencial (Ratio)": round(ratio, 1),
                    "Artigos no Alvo": n_especifico,
                    "Global/Fonte": n_global
                })
                progresso.progress((i+1)/len(lista))
            
            progresso.empty()
            status_text.empty()
            st.session_state.resultado_df = pd.DataFrame(resultados).sort_values(by="Potencial (Ratio)", ascending=False)
            st.rerun()

# 4. EXIBIÇÃO DE RESULTADOS
if st.session_state.resultado_df is not None:
    df = st.session_state.resultado_df
    
    st.markdown("### 🎯 Resultados da Prospecção")
    
    if not df.empty:
        # KPIs
        top_term = df.iloc[0]
        c1, c2, c3 = st.columns(3)
        c1.metric("🏆 Maior Potencial", top_term['Molécula/Alvo'])
        c2.metric("📊 Score (Ratio)", top_term['Potencial (Ratio)'])
        c3.metric("📚 Artigos (Alvo)", top_term['Artigos no Alvo'])
        
        # Gráfico
        fig = px.bar(df.head(20), x="Molécula/Alvo", y="Potencial (Ratio)", color="Status",
                     color_discrete_map={"💎 Blue Ocean": "#00CC96", "🥇 Ouro": "#636EFA", "🔴 Saturado": "#EF553B"})
        st.plotly_chart(fig, use_container_width=True)
        
        # Tabela
        st.dataframe(df, use_container_width=True, hide_index=True)
        
        # Download
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button("📥 Baixar Relatório CSV", csv, "lemos_lambda_report.csv", "text/csv")
    else:
        st.warning(t["tabela_vazia"])

# --- RODAPÉ ---
st.markdown("---")
st.caption(f"© 2025 Guilherme Lemos | {t['footer_citar']}")