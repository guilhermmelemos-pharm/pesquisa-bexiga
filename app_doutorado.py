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
    div[data-testid="stImage"] img { height: 160px !important; object-fit: cover !important; border-radius: 10px !important; }
    .stAlert { padding: 0.5rem; margin-bottom: 1rem; border-radius: 8px; }
    </style>
""", unsafe_allow_html=True)

# --- ESTADO ---
if 'alvos_val' not in st.session_state: st.session_state.alvos_val = ""
if 'resultado_df' not in st.session_state: st.session_state.resultado_df = None
if 'news_index' not in st.session_state: st.session_state.news_index = 0
if 'input_alvo' not in st.session_state: st.session_state.input_alvo = ""
if 'input_fonte' not in st.session_state: st.session_state.input_fonte = ""
if 'input_manual' not in st.session_state: st.session_state.input_manual = ""

# --- IDIOMA ---
lang_opt = st.sidebar.radio("🌐 Language:", ["🇧🇷 PT", "🇺🇸 EN"], horizontal=True)
lang = "pt" if "PT" in lang_opt else "en"
t = c.TEXTOS[lang]

# ==========================================
# BARRA LATERAL (HOW TO CITE)
# ==========================================
st.sidebar.markdown("---")
with st.sidebar.expander(t["citar_titulo"], expanded=True):
    st.code(t["citar_texto"], language="text")
    st.link_button(t["link_doi"], "https://doi.org/10.5281/zenodo.17958507")
st.sidebar.markdown("---")

# ==========================================
# UI HELPERS
# ==========================================
def limpar_campo(chave_session):
    st.session_state[chave_session] = ""

def aplicar_preset_lemos(textos):
    """Carrega as configurações do doutorado do Lemos via Callback"""
    st.session_state.input_alvo = c.PRESET_LEMOS["alvo"]
    st.session_state.input_fonte = c.PRESET_LEMOS["fonte"]
    # Carrega a lista base de forma segura
    adicionar_termos_seguro(c.CANDIDATOS_MINERACAO, textos)
    st.toast(textos["toast_preset"], icon="🎓")

def adicionar_termos_seguro(novos_termos_lista, textos):
    """Adiciona termos evitando duplicatas"""
    atuais = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
    atuais_upper = [x.upper() for x in atuais]
    adicionados = []
    
    for termo in novos_termos_lista:
        t_limpo = termo.strip()
        if t_limpo and (t_limpo.upper() not in atuais_upper):
            atuais.append(t_limpo)
            atuais_upper.append(t_limpo.upper())
            adicionados.append(t_limpo)
    
    st.session_state.alvos_val = ", ".join(atuais)
    # Feedback visual (Toast) só se houver novos, mas sem bloquear callback
    return len(adicionados)

@st.fragment(run_every=60) 
def exibir_radar_cientifico(lang_code, textos):
    news_list = bk.buscar_todas_noticias(lang_code)
    if not news_list: return
    total_news = len(news_list)
    idx = st.session_state.news_index % total_news
    batch = news_list[idx:idx+3]
    st.session_state.news_index += 3
    with st.container(border=True):
        st.caption(textos["radar_titulo"])
        cols = st.columns(3)
        for i, n in enumerate(batch):
            with cols[i]:
                st.image(n['img'], use_container_width=True)
                st.markdown(f"**{n['titulo'][:75]}...**")
                st.caption(f"{n['bandeira']} {n['fonte']}")
                st.link_button(textos["btn_ler_feed"], n['link'], use_container_width=True)

def processar_upload(textos):
    uploaded_file = st.session_state.get('uploader_key')
    if uploaded_file is not None:
        try:
            content = uploaded_file.getvalue().decode("utf-8")
            termos_importados = [x.strip() for x in content.replace("\n", ",").split(",") if x.strip()]
            count = adicionar_termos_seguro(termos_importados, textos)
            st.toast(f"{textos['toast_import']} ({count})", icon="📂")
        except: st.error(textos["erro_ler"])

# ==========================================
# UI PRINCIPAL
# ==========================================
st.title(t["titulo_desk"])
st.caption(t["subtitulo"])

# 1. RADAR CIENTÍFICO
exibir_radar_cientifico(lang, t)

st.divider()

# 2. ÁREA DE INPUT (WIZARD)
with st.container(border=True):
    col_input, col_config = st.columns([2, 1])
    
    with col_input:
        st.subheader(t["step_1"])
        st.warning(t["aviso_pubmed"])
        
        email_user = st.text_input(t["label_email"], placeholder=t["holder_email"])
        
        # Alvo
        c_in, c_bin = st.columns([8, 1], vertical_alignment="bottom")
        with c_in: alvo = st.text_input(t["label_alvo"], key="input_alvo", placeholder=t["holder_alvo"])
        with c_bin: st.button("🗑️", key="lixo_alvo", on_click=limpar_campo, args=("input_alvo",), help=t["btn_limpar"])

        # AÇÕES (Magic + Manual)
        c_magic, c_manual = st.columns(2)
        with c_magic:
            if st.button(t["btn_magic"], type="primary", use_container_width=True):
                if not email_user or not alvo:
                    st.error(t["erro_campos"])
                else:
                    with st.status(t["prog_magic"], expanded=True) as status:
                        st.write(t["status_minerando"])
                        novos_termos = bk.buscar_alvos_emergentes_pubmed(alvo, email_user)
                        st.write(t["status_filtrando"])
                        count = adicionar_termos_seguro(c.CANDIDATOS_MINERACAO + novos_termos, t)
                        status.update(label=t["status_pronto"], state="complete", expanded=False)
                        st.toast(f"✅ {len(novos_termos)} novos termos!", icon="🧬")

        with c_manual:
            with st.popover(t["label_manual"], use_container_width=True):
                termo_man = st.text_input("Termo", key="input_manual", placeholder=t["holder_manual"])
                if st.button(t["btn_add_manual"], use_container_width=True):
                    if termo_man:
                        count = adicionar_termos_seguro([x.strip() for x in termo_man.split(",")], t)
                        st.session_state.input_manual = ""
                        st.rerun()

    with col_config:
        st.subheader("Config")
        
        # --- CORREÇÃO DO ERRO ---
        # Botão Preset usando on_click para evitar o StreamlitAPIException
        st.button(t["btn_preset"], on_click=aplicar_preset_lemos, args=(t,), use_container_width=True)
            
        st.subheader(t["label_periodo"])
        anos_range = st.slider("Anos", 2000, datetime.now().year, (2015, datetime.now().year), label_visibility="collapsed")
        st.caption(f"{anos_range[0]} - {anos_range[1]}")
        
        st.markdown("---")
        st.caption(t["desc_fonte"])
        c_src, c_sbin = st.columns([8, 1], vertical_alignment="bottom")
        with c_src: contexto = st.text_input("Context", key="input_fonte", placeholder=t["holder_fonte"], label_visibility="collapsed")
        with c_sbin: st.button("🗑️", key="lixo_fonte", on_click=limpar_campo, args=("input_fonte",), help=t["btn_limpar"])
        
        st.markdown("---")
        st.write(t["titulo_import"])
        st.file_uploader(t["desc_import"], type=["csv", "txt"], key="uploader_key", on_change=processar_upload, args=(t,), label_visibility="collapsed")

# 3. ÁREA DE ANÁLISE
st.divider()
st.subheader(t["step_2"])

if st.session_state.alvos_val:
    with st.expander(t["ver_editar"], expanded=False):
        c1, c2 = st.columns([5,1])
        with c1: st.session_state.alvos_val = st.text_area("Termos", value=st.session_state.alvos_val, height=100, label_visibility="collapsed")
        with c2:
            if st.button(t["btn_limpar_tudo"]): st.session_state.alvos_val = ""
            st.caption(f"{t['qtd_termos']} {len(st.session_state.alvos_val.split(',')) if st.session_state.alvos_val else 0}")

    if st.button(t["analise_btn"], use_container_width=True):
        if not email_user: st.error(t["erro_email"])
        else:
            lista = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
            resultados = []
            prog = st.progress(0); stt = st.empty()
            ano_ini, ano_fim = anos_range
            
            for i, item in enumerate(lista):
                stt.caption(f"🔍 {item}...")
                time.sleep(0.05)
                
                termo_contexto = contexto if contexto else None
                n_global = bk.consultar_pubmed_count(item, termo_contexto, email_user, ano_ini, ano_fim)
                n_especifico = bk.consultar_pubmed_count(item, alvo, email_user, ano_ini, ano_fim)
                
                ratio = n_global / n_especifico if n_especifico > 0 else n_global
                
                if n_especifico == 0 and n_global > 20: tag = "💎 Blue Ocean"
                elif ratio > 10: tag = "🥇 Ouro"
                elif ratio < 2: tag = "🔴 Saturado"
                else: tag = "⚖️ Neutro"
                
                resultados.append({t["col_mol"]: item, t["col_status"]: tag, t["col_ratio"]: round(ratio, 1), t["col_art_alvo"]: n_especifico, t["col_global"]: n_global})
                prog.progress((i+1)/len(lista))
            
            prog.empty(); stt.empty()
            st.session_state.resultado_df = pd.DataFrame(resultados).sort_values(by=t["col_ratio"], ascending=False)
            st.rerun()

# 4. EXIBIÇÃO
if st.session_state.resultado_df is not None:
    df = st.session_state.resultado_df
    st.markdown(f"### {t['resultados']}")
    if not df.empty:
        top = df.iloc[0]
        c1, c2, c3 = st.columns(3)
        c1.metric(t["metrica_potencial"], top[t["col_mol"]])
        c2.metric(t["metrica_score"], top[t["col_ratio"]])
        c3.metric(t["metrica_artigos"], top[t["col_art_alvo"]])
        
        fig = px.bar(df.head(20), x=t["col_mol"], y=t["col_ratio"], color=t["col_status"],
                     color_discrete_map={"💎 Blue Ocean": "#00CC96", "🥇 Ouro": "#636EFA", "🔴 Saturado": "#EF553B"})
        st.plotly_chart(fig, use_container_width=True)
        st.dataframe(df, use_container_width=True, hide_index=True)
        st.download_button(t["btn_baixar"], df.to_csv(index=False).encode('utf-8'), "lemos_lambda_report.csv", "text/csv")
    else: st.warning(t["tabela_vazia"])

st.markdown("---"); st.caption(f"© 2025 Guilherme Lemos | {t['footer_citar']}")