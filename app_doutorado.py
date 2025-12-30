"""
Lemos Lambda: Deep Science Prospector
Copyright (c) 2025 Guilherme Lemos
Licensed under the MIT License.
Version: 2.0 (Stable)
"""
import streamlit as st
import pandas as pd
import plotly.express as px
from datetime import datetime
import time
import scipy.stats as stats

# Importação dos módulos locais
import constantes as c
import backend as bk 

# --- 1. CONFIGURAÇÃO DA PÁGINA ---
st.set_page_config(
    page_title="λ Lemos Lambda v2.0: Deep Science Prospector",
    page_icon="λ",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# --- 2. ESTADO E INICIALIZAÇÃO ---
DEFAULTS = {
    "pagina": "home",
    "alvos_val": "",
    "resultado_df": None,
    "news_index": 0,
    "input_alvo": "",
    "input_fonte": "",
    "input_email": "",
    "artigos_detalhe": None,
    "email_guardado": "",
    "alvo_guardado": "",
    "lang": "pt",
    "api_key_usuario": "",
    "usar_ia_faxina": True,
    "ia_global_switch": True,
}

for k, v in DEFAULTS.items():
    if k not in st.session_state:
        st.session_state[k] = v

def get_textos():
    return c.TEXTOS.get(st.session_state.lang, c.TEXTOS["pt"])

t = get_textos()

# --- 3. CSS CUSTOMIZADO ---
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

# --- 4. FUNÇÕES DE SUPORTE (FRONTEND) ---

def mudar_idioma(novo_lang):
    st.session_state.lang = novo_lang
    resetar_pesquisa()

def resetar_pesquisa():
    st.session_state.pagina = "home"
    st.session_state.resultado_df = None
    st.session_state.artigos_detalhe = None

def limpar_campo(k):
    st.session_state[k] = ""

def adicionar_termos_seguro(lista, textos):
    atuais = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
    atuais_up = [x.upper() for x in atuais]
    blacklist = [b.lower() for b in c.BLACKLIST_GERAL]
    adicionados = []
    for item in lista:
        termo = item.strip()
        if not termo or any(b in termo.lower() for b in blacklist):
            continue
        if termo.upper() not in atuais_up:
            atuais.append(termo)
            atuais_up.append(termo.upper())
            adicionados.append(termo)
    st.session_state.alvos_val = ", ".join(atuais)
    return len(adicionados)

@st.cache_data(ttl=3600)
def minerar_cached(alvo, email, api_key, usar_ia):
    # Passa a API Key explicitamente para o backend estável
    return bk.minerar_pubmed(alvo, email, api_key, usar_ia)

def carregar_lista_dinamica_smart(textos):
    if not st.session_state.input_alvo:
        st.error(textos["erro_alvo"])
        return
    
    with st.spinner(f"{textos['status_minerando']} {st.session_state.input_alvo}..."):
        resultado = minerar_cached(
            st.session_state.input_alvo, 
            st.session_state.input_email, 
            st.session_state.api_key_usuario, 
            st.session_state.usar_ia_faxina and st.session_state.ia_global_switch
        )
        novos = resultado.get("termos_indicados", [])
        
        if not novos:
            st.warning("Nenhum alvo novo detectado. Usando presets...")
            for lista in c.PRESETS_FRONTEIRA.values(): novos.extend(lista)
        
        qtd = adicionar_termos_seguro(novos, textos)
        st.toast(textos["toast_atualizado"], icon="✅")
        if qtd > 0: st.success(f"✅ {qtd} {textos['sucesso_carregado']}")

def ir_para_analise(email, contexto, alvo, y_ini, y_fim, textos):
    st.session_state.email_guardado = email
    st.session_state.alvo_guardado = alvo
    lista = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
    
    if not lista:
        st.error("Lista de termos vazia.")
        return

    res = []
    N_PUBMED = 36000000 
    placeholder = st.empty()
    
    with placeholder.container():
        st.markdown(textos["titulo_processando"])
        prog = st.progress(0)
        n_total_alvo = bk.consultar_pubmed_count(alvo, "", email, 1900, 2030)
        if n_total_alvo == 0: n_total_alvo = 1

        for i, item in enumerate(lista):
            n_base = bk.consultar_pubmed_count(item, contexto if contexto else "", email, y_ini, y_fim)
            n_especifico = bk.consultar_pubmed_count(item, alvo, email, y_ini, y_fim)
            
            # Estatística Fisher Exact
            a, b = n_especifico, max(0, n_base - n_especifico)
            c_val = max(0, n_total_alvo - n_especifico)
            d = max(0, N_PUBMED - (a + b + c_val))
            
            try: _, p_value = stats.fisher_exact([[a, b], [c_val, d]], alternative='greater')
            except: p_value = 1.0
            
            expected = (max(0.1, n_base) * n_total_alvo) / N_PUBMED
            enrichment = (n_especifico + 0.1) / expected
            
            # Classificação Lambda 2.0
            tag, score_sort = textos["tag_neutral"], 0
            if n_especifico <= 5:
                if n_base > 20: tag, score_sort = textos["tag_blue_ocean"], 1000
                else: tag, score_sort = textos["tag_ghost"], 0
            elif n_especifico <= 25: tag, score_sort = textos["tag_embryonic"], 500
            else:
                if enrichment > 5: tag, score_sort = textos["tag_gold"], 100
                elif enrichment > 1.5: tag, score_sort = textos["tag_trending"], 200 
                else: tag, score_sort = textos["tag_saturated"], 10

            res.append({
                textos["col_mol"]: item, textos["col_status"]: tag, 
                textos["col_ratio"]: round(float(enrichment), 1), "P-Value": f"{p_value:.4f}", 
                textos["col_art_alvo"]: n_especifico, textos["col_global"]: n_base, "_sort": score_sort
            })
            prog.progress((i+1)/len(lista))

    st.session_state.resultado_df = pd.DataFrame(res).sort_values(by=["_sort", textos["col_ratio"]], ascending=[False, False])
    st.session_state.pagina = 'resultados'
    st.rerun()

# --- 5. INTERFACE ---
c_logo, c_lang = st.columns([10, 2])
with c_lang:
    c1, c2 = st.columns(2)
    with c1: st.button("🇧🇷", on_click=mudar_idioma, args=("pt",), key="btn_pt")
    with c2: st.button("🇺🇸", on_click=mudar_idioma, args=("en",), key="btn_en")

if st.session_state.pagina == 'resultados':
    if st.button(t["btn_voltar"]):
        st.session_state.pagina = 'home'
        st.rerun()
    
    st.title(t["resultados"])
    df = st.session_state.resultado_df
    if df is not None:
        top = df.iloc[0]
        col1, col2, col3 = st.columns(3)
        col1.metric(t["metric_top"], top[t["col_mol"]], delta=top[t["col_status"]])
        col2.metric(t["metric_score"], top[t["col_ratio"]])
        col3.metric("P-Value", top["P-Value"])
        
        st.plotly_chart(px.bar(df.head(20), x=t["col_mol"], y=t["col_ratio"], color=t["col_status"]), use_container_width=True)
        st.dataframe(df.drop(columns=["_sort"]), use_container_width=True, hide_index=True)
        
        st.divider()
        st.subheader(t["header_leitura"])
        alvo_sel = st.selectbox(t["label_investigar"], df[t["col_mol"]].tolist())
        
        if st.button(t["btn_investigar"]):
            with st.spinner(t["spinner_investigando"]):
                st.session_state.artigos_detalhe = bk.buscar_resumos_detalhados(
                    alvo_sel, st.session_state.alvo_guardado, st.session_state.email_guardado, 2015, 2026
                )
        
        if st.session_state.artigos_detalhe:
            for i, art in enumerate(st.session_state.artigos_detalhe):
                with st.expander(f"📄 {art.get('Title', 'Untitled')}"):
                    # Ajustado para usar 'Abstract' conforme backend v2.6
                    st.write(art.get("Abstract", "Sem resumo disponível."))
                    if st.session_state.api_key_usuario and st.session_state.ia_global_switch:
                        if st.button(f"🤖 Analyze PhD", key=f"btn_ia_{i}"):
                            # Passagem corrigida para 4 argumentos
                            res = bk.analisar_abstract_com_ia(
                                art['Title'], art.get('Abstract',''), 
                                st.session_state.api_key_usuario, st.session_state.lang
                            )
                            st.info(res)
                    st.link_button("PubMed", art.get("Link", "#"))
else:
    st.markdown(f'<p class="header-style">{t["titulo_desk"]}</p>', unsafe_allow_html=True)
    st.markdown(f'<p class="sub-header-style">{t["subtitulo"]}</p>', unsafe_allow_html=True)
    
    col1, col2 = st.columns([2, 1])
    with col1:
        st.text_input(t["label_email"], key="input_email")
        st.text_input(t["label_alvo"], key="input_alvo")
        if st.button(t["btn_auto"], use_container_width=True, type="primary"):
            carregar_lista_dinamica_smart(t)
        
        if st.session_state.alvos_val:
            st.text_area(t["expander_lista"], key="alvos_val", height=150)
            if st.button(t["btn_executar"], use_container_width=True, type="primary"):
                ir_para_analise(st.session_state.input_email, st.session_state.input_fonte, st.session_state.input_alvo, 2015, 2026, t)
    
    with col2:
        st.session_state.api_key_usuario = st.text_input("Gemini API Key", type="password")
        st.toggle(t["label_ia_global"], key="ia_global_switch")
        st.text_input(t["label_contexto"], key="input_fonte")

# RODAPÉ v2.0
st.markdown("---")
c1, c2 = st.columns([2, 1])
with c1:
    st.caption(t["footer_citar"])
    with st.expander(t["citar_titulo"]): st.code(t["citar_texto"], language="text")
with c2:
    st.caption(t["apoio_titulo"]); st.text_input("Pix:", value="960f3f16-06ce-4e71-9b5f-6915b2a10b5a")
