"""
Lemos Lambda: Deep Science Prospector
Copyright (c) 2025 Guilherme Lemos
Licensed under the MIT License.
"""
import streamlit as st

# Configuração da página DEVE ser o primeiro comando Streamlit
st.set_page_config(page_title="λ Lemos Lambda: Deep Science Prospector", page_icon="λ", layout="wide", initial_sidebar_state="collapsed")

import pandas as pd
import plotly.express as px
from datetime import datetime
import time
import scipy.stats as stats
import constantes as c
import backend as bk

# --- 2. ESTADO ---
defaults = {
    'pagina': 'home', 'alvos_val': "", 'resultado_df': None, 'news_index': 0,
    'input_alvo': "", 'input_fonte': "", 'input_email': "", 'artigos_detalhe': None,
    'email_guardado': "", 'alvo_guardado': "", 'lang': 'pt',
    'api_key_usuario': "", 'usar_ia_faxina': True # Novo estado padrão
}
for k, v in defaults.items():
    if k not in st.session_state: st.session_state[k] = v

# Garante que a chave não suma ao recarregar
if 'api_key_usuario' not in st.session_state:
    st.session_state.api_key_usuario = ""

def get_textos(): return c.TEXTOS.get(st.session_state.lang, c.TEXTOS["pt"])
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
    </style>
""", unsafe_allow_html=True)

# ================= LÓGICA =================

def mudar_idioma(novo_lang): 
    st.session_state.lang = novo_lang
    resetar_pesquisa()
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
    # Pega o estado do botão de IA
    usar_ia = st.session_state.usar_ia_faxina
    
    if not alvo: st.error(textos["erro_alvo"]); return
    
    existentes = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
    lista_mestra = list(set(existentes)) 
    
    with st.spinner(f"{textos['status_minerando']} {alvo}..."):
        # Passa a preferência do usuário para o backend
        novos = bk.buscar_alvos_emergentes_pubmed(alvo, email, usar_ia=usar_ia)
        
        if novos: 
            lista_mestra.extend(novos)
        else:
            st.warning("Mineração retornou vazio. Carregando sugestões padrão...")
            for lista in c.PRESETS_FRONTEIRA.values():
                lista_mestra.extend(lista)
    
    qtd = adicionar_termos_seguro(lista_mestra, textos)
    st.toast(textos["toast_atualizado"], icon="✅")
    if qtd > 0: st.success(f"✅ {qtd} {textos['sucesso_carregado']}")

def ir_para_analise(email, contexto, alvo, y_ini, y_fim, textos):
    st.session_state.email_guardado = email
    st.session_state.alvo_guardado = alvo
    lista = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
    res = []
    
    with st.spinner(textos["spinner_analise"]):
        n_total_alvo = bk.consultar_pubmed_count(alvo, "", email, 1900, 2030)
        if n_total_alvo == 0: n_total_alvo = 1
    
    N_PUBMED = 36000000
    placeholder = st.empty()
    with placeholder.container():
        st.markdown(textos["titulo_processando"])
        prog = st.progress(0)
    
    for i, item in enumerate(lista):
        time.sleep(0.01) 
        n_global = bk.consultar_pubmed_count(item, "", email, y_ini, y_fim)
        n_especifico = bk.consultar_pubmed_count(item, alvo, email, y_ini, y_fim)
        
        a, b, c_val = n_especifico, n_global - n_especifico, n_total_alvo - n_especifico
        d = max(0, N_PUBMED - (a + max(0,b) + max(0,c_val)))
        
        try: _, p_value = stats.fisher_exact([[a, max(0,b)], [max(0,c_val), d]], alternative='greater')
        except: p_value = 1.0
        
        expected = (n_global * n_total_alvo) / N_PUBMED
        enrichment = (n_especifico + 0.1) / (expected if expected > 0 else 0.00001)
        
        tag, score_sort = textos["tag_neutral"], 0
        if n_especifico == 0:
            if n_global > 100: tag, score_sort = textos["tag_blue_ocean"], 1000
            else: tag, score_sort = textos["tag_ghost"], 0
        elif n_especifico <= 15:
            if n_global > 50: tag, score_sort = textos["tag_embryonic"], 500
            else: tag, score_sort = textos["tag_neutral"], 20
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
    
    placeholder.empty()
    st.session_state.resultado_df = pd.DataFrame(res).sort_values(by=["_sort", textos["col_ratio"]], ascending=[False, False])
    st.session_state.pagina = 'resultados'
    st.rerun()

def resetar_pesquisa():
    st.session_state.pagina = 'home'
    st.session_state.resultado_df = None
    st.session_state.artigos_detalhe = None

def processar_upload(textos):
    f = st.session_state.get('uploader_key')
    if f:
        try:
            content = f.getvalue().decode("utf-8")
            l = [x.strip() for x in content.replace("\n", ",").split(",") if x.strip()]
            adicionar_termos_seguro(l, textos)
            st.toast(textos['toast_importado'], icon="📂")
        except: st.error(textos["erro_arquivo"])

@st.fragment(run_every=60) 
def exibir_radar_cientifico(lang_code, textos):
    try:
        news = bk.buscar_todas_noticias(lang_code)
        if not news: return
        idx = st.session_state.get('news_index', 0) % len(news)
        batch = news[idx:idx+3]
        st.session_state.news_index = idx + 3
        with st.container(border=True):
            st.caption(textos["radar_titulo"])
            cols = st.columns(3)
            for i, n in enumerate(batch):
                with cols[i]:
                    st.markdown(f"**{n['titulo'][:75]}...**")
                    st.caption(f"🔬 {n.get('fonte','')}")
                    st.link_button(textos["btn_ler"], n['link'], use_container_width=True)
    except: pass

# --- HEADER ---
c_logo, c_lang = st.columns([10, 2])
with c_lang:
    c1, c2 = st.columns(2)
    with c1: st.button("🇧🇷", key="pt_btn", on_click=mudar_idioma, args=("pt",))
    with c2: st.button("🇺🇸", key="en_btn", on_click=mudar_idioma, args=("en",))

# --- UI PRINCIPAL ---
if st.session_state.pagina == 'resultados':
    c_back, c_tit = st.columns([1, 5])
    with c_back: st.button(t["btn_voltar"], on_click=resetar_pesquisa, use_container_width=True)
    with c_tit: st.title(t["resultados"])
    
    df = st.session_state.resultado_df
    if df is not None and not df.empty:
        top = df.iloc[0]
        c1, c2, c3 = st.columns(3)
        c1.metric(t["metric_top"], top[t["col_mol"]], delta=top[t["col_status"]])
        c2.metric(t["metric_score"], top[t["col_ratio"]])
        c3.metric("P-Valor", top["P-Value"])
        
        st.subheader(t["header_heatmap"])
        st.plotly_chart(px.bar(df.drop(columns=["_sort"]).head(25), x=t["col_mol"], y=t["col_ratio"], color=t["col_status"]), use_container_width=True)
        st.dataframe(df.drop(columns=["_sort"]), use_container_width=True, hide_index=True)
        
        st.divider()
        st.subheader(t["header_leitura"])
        lista_alvos = sorted(df[t["col_mol"]].unique().tolist())
        
        c_sel, c_btn = st.columns([3, 1], vertical_alignment="bottom")
        with c_sel: sel = st.selectbox(t["label_investigar"], lista_alvos)
        with c_btn: 
            if st.button(f"🔎 Buscar Papers", use_container_width=True):
                with st.spinner(t["spinner_investigando"]):
                    st.session_state.artigos_detalhe = bk.buscar_resumos_detalhados(sel, st.session_state.alvo_guardado, st.session_state.email_guardado, 2015, 2025)
                    st.rerun()

        if st.session_state.artigos_detalhe:
            st.info(f"Foram encontrados {len(st.session_state.artigos_detalhe)} artigos recentes sobre {sel}.")
            for i, art in enumerate(st.session_state.artigos_detalhe):
                with st.expander(f"📄 {art['Title']}", expanded=False):
                    st.caption(f"**Keywords/Contexto:** {art.get('Info_IA', 'N/A')[:200]}...")
                    c_ia, c_link = st.columns([1, 1])
                    with c_ia:
                        if not st.session_state.api_key_usuario:
                            nova_chave = st.text_input("Cole sua Google API Key:", type="password", key=f"key_input_{i}", help="Necessário para análise.")
                            if nova_chave:
                                st.session_state.api_key_usuario = nova_chave
                                st.rerun()
                        if st.session_state.api_key_usuario:
                            if st.button(f"🤖 Analisar este artigo", key=f"btn_ia_{i}"):
                                with st.spinner("Analisando..."):
                                    resumo = bk.analisar_abstract_com_ia(art['Title'], art.get('Info_IA', ''), st.session_state.api_key_usuario, st.session_state.lang)
                                    st.markdown(f"<div style='background-color: #262730; color: #ffffff; padding: 15px; border-radius: 8px; border-left: 5px solid #FF4B4B; margin-top: 10px;'><small style='color: #FF4B4B;'>🧠 <b>Análise Lemos Lambda:</b></small><br><span style='font-size: 1.1em;'>{resumo}</span></div>", unsafe_allow_html=True)
                    with c_link: st.link_button("🔗 Abrir no PubMed", art['Link'], use_container_width=True)

else:
    st.title(t["titulo_desk"]); st.caption(t["subtitulo"])
    exibir_radar_cientifico(st.session_state.lang, t)
    st.divider()
    col_main, col_config = st.columns([2, 1])

    with col_main:
        st.subheader("1. Definição do Escopo")
        st.text_input(t["label_email"], key="input_email", placeholder=t["holder_email"])
        c_in, c_trash = st.columns([9, 1], vertical_alignment="bottom")
        with c_in: st.text_input(t["label_alvo"], key="input_alvo", placeholder=t["holder_alvo"])
        with c_trash: st.button("🗑️", on_click=limpar_campo, args=("input_alvo",))
        
        st.markdown('<div class="big-button">', unsafe_allow_html=True)
        if st.button(t["btn_auto"], use_container_width=True): carregar_lista_dinamica_smart(t)
        st.markdown('</div>', unsafe_allow_html=True)

        with st.expander(t["expander_presets"], expanded=False):
            if st.button(t["btn_add_all"], use_container_width=True):
                todos = []
                for lista in c.PRESETS_FRONTEIRA.values(): todos.extend(lista)
                adicionar_termos_seguro(todos, t); st.rerun()
            st.divider() 
            escolha = st.selectbox(t["label_categoria"], options=list(c.PRESETS_FRONTEIRA.keys()))
            if st.button(f"{t['btn_add_preset']} {escolha}", use_container_width=True):
                adicionar_termos_seguro(c.PRESETS_FRONTEIRA[escolha], t); st.rerun()

    with col_config:
        st.subheader(t["header_config"])
        with st.expander(t["expander_ia"], expanded=True):
            st.session_state.api_key_usuario = st.text_input("Google API Key", type="password", value=st.session_state.api_key_usuario)
            
            # --- NOVO BOTÃO DE TOGGLE DA IA ---
            st.toggle("✨ Ativar Curadoria por IA", key="usar_ia_faxina", help="Usa a IA para limpar a lista de mineração, removendo termos clínicos irrelevantes.")
            
            st.markdown(f"[{t['link_key']}](https://aistudio.google.com/app/apikey)")
            
        st.divider()
        anos = st.slider(t["slider_tempo"], 2000, datetime.now().year, (2015, datetime.now().year))
        st.text_input(t["label_contexto"], key="input_fonte")
        st.file_uploader(t["uploader_label"], type=["csv", "txt"], key="uploader_key", on_change=processar_upload, args=(t,))

    st.divider()
    if st.session_state.alvos_val:
        st.success(f"✅ {len(st.session_state.alvos_val.split(','))} alvos prontos.")
        with st.expander(t["expander_lista"]):
            st.text_area("", key="alvos_val", height=150)
            st.button("termos indicados", on_click=limpar_lista_total)

        if st.button(t["btn_executar"], type="primary", use_container_width=True):
            if not st.session_state.input_email: st.error(t["erro_email"])
            else: ir_para_analise(st.session_state.input_email, st.session_state.input_fonte, st.session_state.input_alvo, anos[0], anos[1], t)

st.markdown("---")
cf1, cf2 = st.columns([2, 1])
with cf1:
    st.caption(t["footer_citar"])
    with st.expander(t["citar_titulo"]): st.code(t["citar_texto"], language="text")
with cf2:
    st.caption(t["apoio_titulo"])
    st.text_input("Chave Pix:", value="960f3f16-06ce-4e71-9b5f-6915b2a10b5a", disabled=False)
