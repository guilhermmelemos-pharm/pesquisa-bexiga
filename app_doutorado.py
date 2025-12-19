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
Version: 1.6.1 (Clean UX + International)
"""
import streamlit as st

# --- 1. CONFIGURAÇÃO DO CABEÇALHO CORRIGIDA ---
st.set_page_config(
    page_title="λ Lemos Lambda: Deep Science Prospector",
    page_icon="λ",
    layout="wide",
    initial_sidebar_state="collapsed"
)

import pandas as pd
import plotly.express as px
from datetime import datetime
import time
import constantes as c
import backend as bk

# --- 2. ESTADO E IDIOMA ---
defaults = {
    'pagina': 'home', 'alvos_val': "", 'resultado_df': None, 'news_index': 0,
    'input_alvo': "", 'input_fonte': "", 'input_email': "", 'artigos_detalhe': None,
    'email_guardado': "", 'alvo_guardado': "", 'lang': 'pt'
}
for k, v in defaults.items():
    if k not in st.session_state: st.session_state[k] = v

def get_textos():
    return c.TEXTOS.get(st.session_state.lang, c.TEXTOS["pt"])

t = get_textos()

# --- 3. CSS ---
st.markdown("""
    <style>
    .stButton button { border-radius: 12px; height: 50px; font-weight: bold; }
    div[data-testid="stMetricValue"] { font-size: 1.8rem !important; }
    div[data-testid="stImage"] img { height: 160px !important; object-fit: cover !important; width: 100% !important; }
    
    .big-button button {
        background-color: #FF4B4B !important;
        color: white !important;
        border: none;
        font-size: 1.1rem !important;
    }
    .flag-btn button {
        background: none !important; border: none !important; font-size: 1.5rem !important;
    }
    .footer-text {
        text-align: center; font-size: 0.9rem; color: #666; margin-top: 20px;
    }
    </style>
""", unsafe_allow_html=True)

# ================= LÓGICA =================

def mudar_idioma(novo_lang): st.session_state.lang = novo_lang
def limpar_campo(k): st.session_state[k] = ""
def limpar_lista_total(): st.session_state.alvos_val = ""

def adicionar_termos_seguro(lista, textos):
    atuais = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
    atuais_up = [x.upper() for x in atuais]
    blacklist = [b.lower() for b in c.BLACKLIST_GERAL]
    add = []
    for item in lista:
        cl = item.strip()
        if not cl or any(b in cl.lower() for b in blacklist): continue
        if cl.upper() not in atuais_up:
            atuais.append(cl); atuais_up.append(cl.upper()); add.append(cl)
    st.session_state.alvos_val = ", ".join(atuais)
    return len(add)

def carregar_lista_dinamica_smart(textos):
    email, alvo = st.session_state.input_email, st.session_state.input_alvo
    if not alvo: st.error("⚠️ Define Target first!"); return
    
    # Carrega o que já está na caixa + Lista de Elite
    existentes = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
    lista_mestra = list(set(existentes + c.CANDIDATOS_MINERACAO)) 
    
    # CHAMA A MINERAÇÃO PESADA DO BACKEND (Busca no órgão específico)
    with st.spinner(f"{textos.get('status_minerando', 'Mining...')} {alvo}..."):
        novos = bk.buscar_alvos_emergentes_pubmed(alvo, email)
        if novos: lista_mestra.extend(novos)
    
    qtd = adicionar_termos_seguro(lista_mestra, textos)
    st.toast("Updated list!", icon="✅")
    if qtd > 0:
        st.success(f"✅ {qtd} loaded (Deep Mining + Presets). Click EXECUTE below.")

def ir_para_analise(email, contexto, alvo, y_ini, y_fim):
    st.session_state.email_guardado = email
    st.session_state.alvo_guardado = alvo
    lista = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
    res = []
    
    with st.spinner("Calculating Lambda Score (v1.6.1)..."):
        n_total_alvo = bk.consultar_pubmed_count(alvo, "", email, 1900, 2030)
        if n_total_alvo == 0: n_total_alvo = 1
    
    N_PUBMED = 36000000
    placeholder = st.empty()
    with placeholder.container():
        st.markdown("## 🧬 Processing...")
        prog = st.progress(0)
    
    for i, item in enumerate(lista):
        time.sleep(0.01) 
        termo_ctx = contexto if contexto else None
        
        n_global = bk.consultar_pubmed_count(item, termo_ctx, email, y_ini, y_fim)
        n_especifico = bk.consultar_pubmed_count(item, alvo, email, y_ini, y_fim)
        
        # --- CÁLCULO INVERTIDO (OPORTUNIDADE) ---
        raw_score = (n_global + 1) / (n_especifico + 1)
        
        # Bônus para Blue Ocean
        if n_especifico == 0 and n_global > 0:
            raw_score = raw_score * 2 

        tag, score_sort = "⚖️ Neutral", 0
        if n_especifico == 0:
            if n_global > 100: tag, score_sort = "💎 Blue Ocean", 2000
            else: tag, score_sort = "👻 Ghost", 0
        elif n_especifico <= 15:
            if n_global > 50: tag, score_sort = "🌱 Embryonic", 500
            else: tag, score_sort = "⚖️ Neutral", 50
        else:
            if raw_score > 50: tag, score_sort = "🥇 Gold", 200
            elif raw_score > 10: tag, score_sort = "🚀 Trending", 100
            else: tag, score_sort = "🔴 Saturated", 10

        res.append({
            t["col_mol"]: item, 
            t["col_status"]: tag, 
            t["col_ratio"]: round(float(raw_score), 1), 
            t["col_art_alvo"]: n_especifico, 
            t["col_global"]: n_global, 
            "_sort": score_sort
        })
        prog.progress((i+1)/len(lista))
    
    placeholder.empty()
    st.session_state.resultado_df = pd.DataFrame(res).sort_values(by=["_sort", t["col_ratio"]], ascending=[False, False])
    st.session_state.pagina = 'resultados'
    st.rerun()

def resetar_pesquisa():
    st.session_state.pagina = 'home'
    st.session_state.resultado_df = None

def processar_upload(textos):
    f = st.session_state.get('uploader_key')
    if f:
        try:
            c = f.getvalue().decode("utf-8")
            l = [x.strip() for x in c.replace("\n", ",").split(",") if x.strip()]
            n = adicionar_termos_seguro(l, textos)
            st.toast(f"Imported ({n})", icon="📂")
        except: st.error("Error reading file.")

@st.fragment(run_every=60) 
def exibir_radar_cientifico(lang_code, textos):
    try:
        news = bk.buscar_todas_noticias(lang_code)
        if not news or not isinstance(news, list): return
        if 'news_index' not in st.session_state: st.session_state.news_index = 0
        idx = st.session_state.news_index % len(news)
        batch = news[idx:idx+3]
        st.session_state.news_index += 3
        
        with st.container(border=True):
            st.caption(textos.get("radar_titulo", "Radar"))
            cols = st.columns(3)
            for i, n in enumerate(batch):
                with cols[i]:
                    if n.get('img'): st.image(n['img'], use_container_width=True)
                    st.markdown(f"**{n['titulo'][:75]}...**")
                    st.caption(f"{n.get('bandeira','')} {n.get('fonte','')}")
                    st.link_button("Read", n['link'], use_container_width=True)
    except: return

# --- HEADER (Bandeiras) ---
c_logo, c_lang = st.columns([10, 2])
with c_lang:
    c1, c2 = st.columns(2)
    with c1: 
        if st.button("🇧🇷", key="head_btn_pt"): mudar_idioma("pt"); st.rerun()
    with c2: 
        if st.button("🇺🇸", key="head_btn_en"): mudar_idioma("en"); st.rerun()

# --- UI PRINCIPAL ---
if st.session_state.pagina == 'resultados':
    c_back, c_tit = st.columns([1, 5])
    with c_back: st.button("⬅ Back", on_click=resetar_pesquisa, use_container_width=True)
    with c_tit: st.title(t["titulo_desk"])
    
    df = st.session_state.resultado_df
    if df is not None and not df.empty:
        top = df.iloc[0]
        c1, c2, c3 = st.columns(3)
        c1.metric("Top Opportunity", top[t["col_mol"]], delta=top[t["col_status"]])
        c2.metric("Lambda Score", top[t["col_ratio"]])
        c3.metric("Hits (Local)", top[t["col_art_alvo"]])
        
        st.subheader("Heatmap")
        st.plotly_chart(px.bar(df.drop(columns=["_sort"]).head(25), x=t["col_mol"], y=t["col_ratio"], color=t["col_status"]), use_container_width=True)
        st.dataframe(df.drop(columns=["_sort"]), use_container_width=True, hide_index=True)
        csv = df.drop(columns=["_sort"]).to_csv(index=False).encode('utf-8')
        st.download_button(t["btn_baixar"], csv, "lemos_lambda_report.csv", "text/csv")
        
        st.divider()
        sel = st.selectbox("Investigate Target:", sorted(df[t["col_mol"]].unique().tolist()))
        if st.button(f"🔎 Analyze {sel}", type="secondary"):
            with st.spinner(f"Searching..."):
                st.session_state.artigos_detalhe = bk.buscar_resumos_detalhados(sel, st.session_state.alvo_guardado, st.session_state.email_guardado, 2015, 2025)
        
        if st.session_state.artigos_detalhe:
            for a in st.session_state.artigos_detalhe:
                with st.expander(f"{a['Title']}"):
                    st.info(f"{a['Resumo_IA']}"); st.link_button("PubMed", a['Link'])

else:
    # HOME
    st.title(t["titulo_desk"]); st.caption(t.get("subtitulo", ""))
    exibir_radar_cientifico(st.session_state.lang, t)
    st.divider()

    col_main, col_config = st.columns([2, 1])

    with col_main:
        st.subheader("1. Scope Definition")
        st.text_input(t["label_email"], key="input_email", placeholder=t.get("holder_email", ""))
        
        c_in, c_trash = st.columns([9, 1], vertical_alignment="bottom")
        with c_in: st.text_input(t["label_alvo"], key="input_alvo", placeholder=t.get("holder_alvo", ""))
        with c_trash: st.button("🗑️", on_click=lambda: limpar_campo("input_alvo"))

        st.write(" ")
        st.markdown('<div class="big-button">', unsafe_allow_html=True)
        if st.button("🧠 AUTO-DETECT TARGETS & START", use_container_width=True):
             carregar_lista_dinamica_smart(t)
        st.markdown('</div>', unsafe_allow_html=True)

        st.write(" ")
        with st.expander("💎 Blue Ocean Presets", expanded=False):
            # --- MUDANÇA: O botão agora chama a MINERAÇÃO (Deep Mining) em vez de só adicionar texto ---
            if st.button("📥 Minerar Blue Oceans (PubMed + Presets)", use_container_width=True):
                # Usa a função inteligente que varre o PubMed pelo órgão (termo_base)
                carregar_lista_dinamica_smart(t)

            st.divider()
            with st.popover("✍️ Manual Input"):
                tm = st.text_input("Term", key="input_manual")
                if st.button("Add Manual"):
                    if tm: adicionar_termos_seguro(tm.split(","), t); st.rerun()

    with col_config:
        st.subheader("Config")
        anos = st.slider("Time Range", 2000, datetime.now().year, (2015, datetime.now().year))
        st.text_input(t.get("label_contexto", "Context"), key="input_fonte", placeholder=t.get("holder_fonte", ""))
        st.file_uploader("Import (CSV/TXT)", type=["csv", "txt"], key="uploader_key", on_change=processar_upload, args=(t,))

    st.divider()
    if st.session_state.alvos_val:
        st.success(f"✅ **{len(st.session_state.alvos_val.split(','))} targets loaded.**")
        
        with st.expander("View/Edit List", expanded=True):
             st.text_area("Edit terms (comma separated):", key="alvos_val", height=150, disabled=False)
             if st.button("Clear List"): limpar_lista_total(); st.rerun()

        if st.button("🚀 EXECUTE ANALYSIS", type="primary", use_container_width=True):
            if not st.session_state.input_email: st.error("E-mail required!")
            else: ir_para_analise(st.session_state.input_email, st.session_state.input_fonte, st.session_state.input_alvo, anos[0], anos[1])

# --- RODAPÉ (FOOTER) ---
st.markdown("---")
st.markdown(f"<div class='footer-text'>{t.get('footer_rights', '© 2025 Guilherme Lemos | Unifesp')}</div>", unsafe_allow_html=True)
st.write("")

st.caption(t.get("footer_citar", "Academic Use"))
with st.expander(t.get("citar_titulo", "Citation"), expanded=False):
    st.code(t.get("citar_texto", ""), language="text")
    )
