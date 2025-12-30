"""
Lemos Lambda: Deep Science Prospector
Copyright (c) 2025 Guilherme Lemos
Licensed under the MIT License.
"""
import streamlit as st

# Configuração da página DEVE ser o primeiro comando Streamlit
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
import scipy.stats as stats
import constantes as c
import backend as bk

# --- 2. ESTADO ---
defaults = {
    'pagina': 'home', 
    'alvos_val': "", 
    'resultado_df': None, 
    'news_index': 0,
    'input_alvo': "", 
    'input_fonte': "", 
    'input_email': "", 
    'artigos_detalhe': None,
    'email_guardado': "", 
    'alvo_guardado': "", 
    'lang': 'pt',
    'api_key_usuario': "" 
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
    div[data-testid="stImage"] { height: 160px !important; overflow: hidden !important; border-radius: 8px !important; }
    div[data-testid="stImage"] img { height: 160px !important; object-fit: cover !important; width: 100% !important; }
    .big-button button { background-color: #FF4B4B !important; color: white !important; border: none; font-size: 1.1rem !important; }
    .stTextArea textarea { font-family: monospace; }
    .news-card { 
        background-color: #262730; padding: 15px; border-radius: 10px; 
        border-left: 5px solid #FF4B4B; margin-bottom: 15px;
    }
    </style>
""", unsafe_allow_html=True)

# ================= LÓGICA =================

def mudar_idioma(novo_lang): 
    st.session_state.lang = novo_lang
    st.rerun()

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
    if not alvo: 
        st.error(textos["erro_alvo"]); return
    
    with st.spinner(f"{textos['status_minerando']} {alvo}..."):
        # Minerando do PubMed através do Backend
        novos = bk.minerar_pubmed(alvo, email, st.session_state.api_key_usuario)
        termos = novos.get("termos_indicados", [])
        
        # Se falhar ou vier vazio, injeta os presets de segurança
        if not termos:
            st.warning("Mineração vazia. Injetando presets de fronteira...")
            for lista in c.PRESETS_FRONTEIRA.values():
                termos.extend(lista)
    
    qtd = adicionar_termos_seguro(termos, textos)
    st.toast(textos["toast_atualizado"], icon="✅")
    if qtd > 0: st.success(f"✅ {qtd} {textos['sucesso_carregado']}")

def ir_para_analise(email, contexto, alvo, y_ini, y_fim, textos):
    st.session_state.email_guardado = email
    st.session_state.alvo_guardado = alvo
    lista = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
    res = []
    
    N_PUBMED = 36000000
    placeholder = st.empty()
    with placeholder.container():
        st.markdown(textos["titulo_processando"])
        prog = st.progress(0)
    
    for i, item in enumerate(lista):
        n_global = bk.consultar_pubmed_count(item, contexto if contexto else "", email, y_ini, y_fim)
        n_especifico = bk.consultar_pubmed_count(item, alvo, email, y_ini, y_fim)
        
        # Estatística Fisher
        a, b = n_especifico, max(0, n_global - n_especifico)
        c_v, d = max(0, n_global - n_especifico), max(0, N_PUBMED - n_global)
        
        try: _, p_value = stats.fisher_exact([[a, b], [c_v, d]], alternative='greater')
        except: p_value = 1.0
        
        expected = (max(0.1, n_global) * 500) / N_PUBMED # Normalização Lambda
        enrichment = (n_especifico + 0.1) / expected
        
        tag, score_sort = textos["tag_neutral"], 0
        if n_especifico <= 5:
            if n_global > 20: tag, score_sort = textos["tag_blue_ocean"], 1000
            else: tag, score_sort = textos["tag_ghost"], 0
        elif n_especifico <= 25: tag, score_sort = textos["tag_embryonic"], 500
        else:
            if enrichment > 5: tag, score_sort = textos["tag_gold"], 100
            elif enrichment > 1.5: tag, score_sort = textos["tag_trending"], 200
            else: tag, score_sort = textos["tag_saturated"], 10

        res.append({
            textos["col_mol"]: item, textos["col_status"]: tag, 
            textos["col_ratio"]: round(float(enrichment), 1), "P-Value": f"{p_value:.4f}", 
            textos["col_art_alvo"]: n_especifico, textos["col_global"]: n_global, "_sort": score_sort
        })
        prog.progress((i+1)/len(lista))
    
    st.session_state.resultado_df = pd.DataFrame(res).sort_values(by=["_sort", textos["col_ratio"]], ascending=[False, False])
    st.session_state.pagina = 'resultados'
    st.rerun()

def processar_upload(textos):
    f = st.session_state.get('uploader_key')
    if f:
        try:
            content = f.getvalue().decode("utf-8")
            termos = [x.strip() for x in content.replace("\n", ",").split(",") if x.strip()]
            n = adicionar_termos_seguro(termos, textos)
            st.toast(f"Importados {n} termos!", icon="📂")
        except: st.error("Erro ao ler arquivo.")

@st.fragment(run_every=60) 
def exibir_radar_cientifico(textos):
    noticias = [
        {"t": "Novas evidências de Piezo1 na bexiga", "d": "Canais mecanossensores regulam a contratilidade do urotélio em modelos de órgão isolado.", "l": "https://pubmed.ncbi.nlm.nih.gov/"},
        {"t": "Inibidores de SGLT2 e Urodinâmica", "d": "Estudos recentes associam gliflozinas a melhoras na hiperatividade detrusora em ratos.", "l": "https://pubmed.ncbi.nlm.nih.gov/"},
        {"t": "H2S: O novo gasotransmissor vesical", "d": "Doador de H2S (GYY4137) reduz a inflamação em cistites hemorrágicas experimentais.", "l": "https://pubmed.ncbi.nlm.nih.gov/"}
    ]
    with st.container(border=True):
        st.caption(textos["radar_titulo"])
        cols = st.columns(3)
        for i, n in enumerate(noticias):
            with cols[i]:
                st.markdown(f"<div class='news-card'><b>{n['t']}</b><br><small>{n['d']}</small></div>", unsafe_allow_html=True)
                st.link_button(textos["btn_ler"], n['l'], use_container_width=True)

# --- HEADER ---
c_logo, c_lang = st.columns([10, 2])
with c_lang:
    c1, c2 = st.columns(2)
    if c1.button("🇧🇷", key="pt"): mudar_idioma("pt")
    if c2.button("🇺🇸", key="en"): mudar_idioma("en")

# --- UI PRINCIPAL ---
if st.session_state.pagina == 'resultados':
    if st.button(t["btn_voltar"]):
        st.session_state.pagina = 'home'
        st.rerun()
    st.title(t["resultados"])
    df = st.session_state.resultado_df
    if df is not None:
        top = df.iloc[0]
        c1, c2, c3 = st.columns(3)
        c1.metric(t["metric_top"], top[t["col_mol"]], delta=top[t["col_status"]])
        c2.metric(t["metric_score"], top[t["col_ratio"]])
        c3.metric("P-Value", top["P-Value"])
        st.plotly_chart(px.bar(df.head(20), x=t["col_mol"], y=t["col_ratio"], color=t["col_status"]), use_container_width=True)
        st.dataframe(df.drop(columns=["_sort"]), use_container_width=True, hide_index=True)
        
        st.divider()
        st.subheader(t["header_leitura"])
        alvo_sel = st.selectbox(t["label_investigar"], df[t["col_mol"]].tolist())
        if st.button(t["btn_investigar"]):
            st.session_state.artigos_detalhe = bk.buscar_resumos_detalhados(alvo_sel, st.session_state.alvo_guardado, st.session_state.email_guardado, 2015, 2026)
        
        if st.session_state.artigos_detalhe:
            for art in st.session_state.artigos_detalhe:
                with st.expander(f"📄 {art['Title']}"):
                    st.write(art['Abstract'])
                    if st.button(f"🤖 Analyze PhD", key=art['Title']):
                        st.info(bk.analisar_abstract_com_ia(art['Title'], art['Abstract'], st.session_state.api_key_usuario, st.session_state.lang))

else:
    # HOME
    st.title(t["titulo_desk"])
    st.caption(t["subtitulo"])
    exibir_radar_cientifico(t)
    st.divider()

    col_main, col_config = st.columns([2, 1])

    with col_main:
        st.subheader("1. Setup")
        st.text_input(t["label_email"], key="input_email")
        st.text_input(t["label_alvo"], key="input_alvo")
        
        if st.button(t["btn_auto"], use_container_width=True, type="primary"):
            carregar_lista_dinamica_smart(t)

        st.write(" ")
        with st.expander(t["expander_presets"]):
            cat = st.selectbox(t["label_categoria"], list(c.PRESETS_FRONTEIRA.keys()))
            if st.button(t["btn_add_preset"]):
                adicionar_termos_seguro(c.PRESETS_FRONTEIRA[cat], t)
                st.rerun()

        st.divider()
        if st.session_state.alvos_val:
            st.text_area(t["expander_lista"], key="alvos_val", height=150)
            if st.button(t["btn_executar"], use_container_width=True, type="primary"):
                ir_para_analise(st.session_state.input_email, st.session_state.input_fonte, st.session_state.input_alvo, 2015, 2026, t)

    with col_config:
        st.subheader(t["header_config"])
        st.session_state.api_key_usuario = st.text_input("Gemini API Key", type="password")
        st.text_input(t["label_contexto"], key="input_fonte")
        st.file_uploader(t["uploader_label"], type=["csv", "txt"], key="uploader_key", on_change=processar_upload, args=(t,))

# RODAPÉ
st.markdown("---")
cf1, cf2 = st.columns([2, 1])
with cf1:
    st.caption(t["footer_citar"])
    with st.expander(t["citar_titulo"]):
        st.code(t["citar_texto"], language="text")
with cf2:
    st.caption(t["apoio_titulo"])
    st.text_input("Pix (E-mail):", value="960f3f16-06ce-4e71-9b5f-6915b2a10b5a", disabled=True)
