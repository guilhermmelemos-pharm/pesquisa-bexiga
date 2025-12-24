"""
Lemos Lambda: Deep Science Prospector
Copyright (c) 2025 Guilherme Lemos
Licensed under the MIT License.
"""

# ================= CONFIG STREAMLIT =================
import streamlit as st
st.set_page_config(
    page_title="λ Lemos Lambda: Deep Science Prospector",
    page_icon="λ",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# ================= IMPORTS =================
import pandas as pd
import plotly.express as px
from datetime import datetime
import time
import scipy.stats as stats

import constantes as c
import backend as bk

# ================= SESSION STATE =================
defaults = {
    "pagina": "home",
    "alvos_val": "",
    "resultado_df": None,
    "artigos_detalhe": None,
    "news_index": 0,
    "input_email": "",
    "input_alvo": "",
    "input_fonte": "",
    "api_key_usuario": "",
    "usar_ia_faxina": True,
    "lang": "pt",
    "email_guardado": "",
    "alvo_guardado": ""
}

for k, v in defaults.items():
    if k not in st.session_state:
        st.session_state[k] = v

def textos():
    return c.TEXTOS.get(st.session_state.lang, c.TEXTOS["pt"])

t = textos()

# ================= CSS =================
st.markdown("""
<style>
.stButton button { border-radius: 12px; height: 48px; font-weight: bold; }
div[data-testid="stMetricValue"] { font-size: 1.7rem !important; }
.big-button button { background-color: #FF4B4B !important; color: white !important; }
.stTextArea textarea { font-family: monospace; }
.header-style { font-size: 2.4rem; font-weight: 700; }
.sub-header-style { font-size: 1.2rem; color: #999; }
</style>
""", unsafe_allow_html=True)

# ================= FUNÇÕES AUX =================
def mudar_idioma(lang):
    st.session_state.lang = lang
    st.session_state.pagina = "home"
    st.rerun()

def limpar_lista():
    st.session_state.alvos_val = ""

def adicionar_termos(lista):
    atuais = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
    atuais_up = {x.upper() for x in atuais}

    for item in lista:
        if item.upper() not in atuais_up:
            atuais.append(item)

    st.session_state.alvos_val = ", ".join(atuais)

# ================= MINERAÇÃO =================
def minerar_alvos():
    if not st.session_state.input_alvo:
        st.error(t["erro_alvo"])
        return

    with st.spinner("Minerando literatura..."):
        novos = bk.buscar_alvos_emergentes_pubmed(
            termo_base=st.session_state.input_alvo,
            email=st.session_state.input_email,
            api_key=st.session_state.api_key_usuario,
            usar_ia=st.session_state.usar_ia_faxina
        )

    if novos:
        adicionar_termos(novos)
        st.success(f"✅ {len(novos)} novos alvos adicionados")
    else:
        st.warning("Nenhum alvo encontrado")

# ================= ANÁLISE =================
def executar_analise(ano_ini, ano_fim):
    lista = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
    if not lista:
        st.error("Lista vazia.")
        return

    res = []
    email = st.session_state.input_email
    alvo = st.session_state.input_alvo

    with st.spinner("Processando análise estatística..."):
        n_total_alvo = bk.consultar_pubmed_count(alvo, "", email, 1900, 2030)
        if n_total_alvo == 0:
            n_total_alvo = 1

    N_PUBMED = 36000000
    prog = st.progress(0)

    for i, termo in enumerate(lista):
        n_global = bk.consultar_pubmed_count(termo, "", email, ano_ini, ano_fim)
        n_especifico = bk.consultar_pubmed_count(termo, alvo, email, ano_ini, ano_fim)

        a = n_especifico
        b = max(0, n_global - n_especifico)
        c_val = max(0, n_total_alvo - n_especifico)
        d = max(0, N_PUBMED - (a + b + c_val))

        try:
            _, p = stats.fisher_exact([[a, b], [c_val, d]], alternative="greater")
        except:
            p = 1.0

        expected = (n_global * n_total_alvo) / N_PUBMED
        enrichment = (a + 0.1) / max(expected, 1e-6)

        res.append({
            t["col_mol"]: termo,
            t["col_ratio"]: round(enrichment, 2),
            "P-Value": f"{p:.4f}",
            t["col_art_alvo"]: a,
            t["col_global"]: n_global
        })

        prog.progress((i + 1) / len(lista))

    st.session_state.resultado_df = pd.DataFrame(res).sort_values(
        by=[t["col_ratio"]], ascending=False
    )
    st.session_state.pagina = "resultados"
    st.rerun()

# ================= HEADER =================
c1, c2 = st.columns([8, 2])
with c2:
    st.button("🇧🇷", on_click=mudar_idioma, args=("pt",))
    st.button("🇺🇸", on_click=mudar_idioma, args=("en",))

# ================= UI =================
if st.session_state.pagina == "home":
    st.markdown(f"<p class='header-style'>λ {t['titulo_desk']}</p>", unsafe_allow_html=True)
    st.markdown(f"<p class='sub-header-style'>{t['subtitulo']}</p>", unsafe_allow_html=True)

    col1, col2 = st.columns([2, 1])

    with col1:
        st.text_input(t["label_email"], key="input_email")
        st.text_input(t["label_alvo"], key="input_alvo")

        st.markdown("<div class='big-button'>", unsafe_allow_html=True)
        st.button(t["btn_auto"], on_click=minerar_alvos, use_container_width=True)
        st.markdown("</div>", unsafe_allow_html=True)

        st.text_area("Alvos selecionados", key="alvos_val", height=160)
        st.button("🧹 Limpar lista", on_click=limpar_lista)

    with col2:
        st.subheader("Configurações")
        st.text_input("Google API Key", type="password", key="api_key_usuario")
        st.toggle("✨ Curadoria por IA", key="usar_ia_faxina", value=True)

        anos = st.slider(
            t["slider_tempo"],
            2000,
            datetime.now().year,
            (2015, datetime.now().year)
        )

        if st.button(t["btn_executar"], type="primary", use_container_width=True):
            if not st.session_state.input_email:
                st.error(t["erro_email"])
            else:
                executar_analise(anos[0], anos[1])

elif st.session_state.pagina == "resultados":
    st.button(t["btn_voltar"], on_click=lambda: st.session_state.update({"pagina": "home"}))
    st.title(t["resultados"])

    df = st.session_state.resultado_df
    st.plotly_chart(
        px.bar(df.head(25), x=t["col_mol"], y=t["col_ratio"]),
        use_container_width=True
    )

    st.dataframe(df, use_container_width=True, hide_index=True)

    st.divider()
    st.subheader("Leitura dirigida")

    sel = st.selectbox("Selecionar alvo", df[t["col_mol"]].tolist())

    if st.button("🔎 Buscar artigos"):
        with st.spinner("Buscando abstracts..."):
            st.session_state.artigos_detalhe = bk.buscar_resumos_detalhados(
                sel,
                st.session_state.input_alvo,
                st.session_state.input_email,
                2018,
                datetime.now().year
            )

    if st.session_state.artigos_detalhe:
        for art in st.session_state.artigos_detalhe:
            with st.expander(f"📄 {art['Title']}"):
                st.caption(art.get("Info_IA", "")[:300])

                if st.button("🤖 Analisar com IA", key=art["Title"]):
                    with st.spinner("Analisando..."):
                        resumo = bk.analisar_abstract_com_ia(
                            art["Title"],
                            art.get("Info_IA", ""),
                            st.session_state.api_key_usuario,
                            st.session_state.lang
                        )
                        st.success(resumo)
