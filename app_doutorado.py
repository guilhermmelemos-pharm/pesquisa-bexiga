"""
Lemos Lambda: Deep Science Prospector
Copyright (c) 2025 Guilherme Lemos
Licensed under the MIT License.
Version: 2.0 (Stable)

Citation:
Lemos, G. (2025). Lemos Lambda: Deep Science Prospector (v2.0). Zenodo. https://doi.org/10.5281/zenodo.18092141
"""
import streamlit as st

# --- 1. CONFIGURAÇÃO DA PÁGINA (Deve ser a primeira linha) ---
st.set_page_config(
    page_title="λ Lemos Lambda v2.0",
    page_icon="λ",
    layout="wide",
    initial_sidebar_state="collapsed"
)

import pandas as pd
import plotly.express as px
from datetime import datetime
import time
import scipy.stats as stats

# Blindagem de Dependências
try:
    import feedparser
except ImportError:
    feedparser = None

# Importação dos módulos locais
import constantes as c
import backend as bk

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
    .news-card { background-color: #262730; padding: 15px; border-radius: 10px; border-left: 5px solid #FF4B4B; margin-bottom: 10px; }
    
    /* Estilo para o Modal */
    div[data-testid="stDialog"] { border-radius: 15px; }
</style>
""", unsafe_allow_html=True)

# --- 4. FUNÇÕES DE SUPORTE ---

def mudar_idioma(novo_lang):
    st.session_state.lang = novo_lang
    st.rerun()

def resetar_pesquisa():
    st.session_state.pagina = "home"
    st.session_state.resultado_df = None
    st.session_state.artigos_detalhe = None

def limpar_campo(k):
    st.session_state[k] = ""

def limpar_lista_total():
    st.session_state.alvos_val = ""

def gerar_bibtex():
    return """@software{lemos_lambda_2025,
  author       = {Lemos, Guilherme},
  title        = {Lemos Lambda: Deep Science Prospector},
  version      = {2.0},
  year         = {2025},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.18092141},
  url          = {https://doi.org/10.5281/zenodo.18092141}
}"""

def adicionar_termos_seguro(lista):
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
            st.warning("Mineração vazia. Injetando Presets de Fronteira...")
            for lista in c.PRESETS_FRONTEIRA.values():
                novos.extend(lista)
        
        qtd = adicionar_termos_seguro(novos)
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
        n_total_alvo = bk.consultar_pubmed_count(alvo, "", email, 1900, 2030)
        if n_total_alvo == 0: n_total_alvo = 1

        for i, item in enumerate(lista):
            termo_ctx = contexto if contexto else ""
            n_base = bk.consultar_pubmed_count(item, termo_ctx, email, y_ini, y_fim)
            n_especifico = bk.consultar_pubmed_count(item, alvo, email, y_ini, y_fim)
            
            a = n_especifico
            b = max(0, n_base - n_especifico)
            c_val = max(0, n_total_alvo - n_especifico)
            d = max(0, N_PUBMED - (a + b + c_val))
            
            try: _, p_value = stats.fisher_exact([[a, b], [c_val, d]], alternative='greater')
            except: p_value = 1.0
            
            expected = (max(0.1, n_base) * n_total_alvo) / N_PUBMED
            enrichment = (n_especifico + 0.1) / expected
            
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

def processar_upload(textos):
    f = st.session_state.get('uploader_key')
    if f:
        try:
            content = f.getvalue().decode("utf-8")
            termos = [x.strip() for x in content.replace("\n", ",").split(",") if x.strip()]
            n = adicionar_termos_seguro(termos)
            st.toast(f"{textos.get('toast_importado','Importados')}: {n}", icon="📂")
        except: st.error("Erro ao ler arquivo.")

@st.fragment(run_every=3600)
def exibir_radar_cientifico(lang_code, textos):
    try:
        if not feedparser: raise ImportError
        news = bk.buscar_todas_noticias(lang_code)
        if not news: raise Exception
    except:
        news = [
            {"titulo": "TRPV4 modulation in Bladder Urothelium", "fonte": "Nature Urology", "link": "https://pubmed.ncbi.nlm.nih.gov"},
            {"titulo": "SGLT2 inhibitors and Voiding Dysfunction", "fonte": "PubMed Trending", "link": "https://pubmed.ncbi.nlm.nih.gov"}
        ]
    with st.container(border=True):
        st.caption(textos.get("radar_titulo", "📡 Science Radar"))
        cols = st.columns(2)
        for i, n in enumerate(news[:2]):
            with cols[i]:
                st.markdown(f"<div class='news-card'><a href='{n['link']}' target='_blank'>{n['titulo']}</a><br><small>{n['fonte']}</small></div>", unsafe_allow_html=True)

# === POP-UP (MODAL) PARA INVESTIGAÇÃO ===
@st.dialog("🔬 Leitura Guiada por IA", width="large")
def abrir_modal_investigacao(alvo, contexto, email, api_key, lang):
    st.markdown(f"### Investigando: **{alvo}**")
    st.caption(f"Contexto: {contexto}")
    
    with st.spinner("Buscando abstracts no PubMed..."):
        artigos = bk.buscar_resumos_detalhados(alvo, contexto, email, 2015, 2026)
    
    if not artigos:
        st.warning("Nenhum artigo relevante encontrado nos últimos 10 anos.")
        return

    for i, art in enumerate(artigos):
        with st.container(border=True):
            st.markdown(f"**📄 {art.get('Title')}**")
            with st.expander("Ler Abstract Original"):
                st.write(art.get('Abstract', 'Resumo não disponível.'))
            
            c1, c2 = st.columns([1, 1])
            with c1:
                st.link_button("🔗 Ver no PubMed", art.get('Link', '#'), use_container_width=True)
            with c2:
                # Ação da IA entra AQUI
                if st.button(f"🤖 Analisar (PhD Mode)", key=f"ai_{alvo}_{i}", use_container_width=True):
                    if not api_key:
                        st.error("Configure sua API Key no menu lateral.")
                    else:
                        with st.spinner("Interpretando mecanismo..."):
                            resumo = bk.analisar_abstract_com_ia(
                                art['Title'], art.get('Abstract', ''), api_key, lang
                            )
                            st.success(resumo)

# --- 5. INTERFACE ---
c_logo, c_lang = st.columns([10, 2])
with c_lang:
    c1, c2 = st.columns(2)
    if c1.button("🇧🇷", key="btn_pt"): mudar_idioma("pt")
    if c2.button("🇺🇸", key="btn_en"): mudar_idioma("en")

if st.session_state.pagina == 'resultados':
    # --- PÁGINA RESULTADOS ---
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
        col3.metric("P-Value (Fisher)", top["P-Value"])
        
        st.plotly_chart(px.bar(df.head(20), x=t["col_mol"], y=t["col_ratio"], color=t["col_status"]), use_container_width=True)
        
        # LEGENDA ESTATÍSTICA (Abaixo do Gráfico)
        with st.expander("ℹ️ Metodologia Estatística (Legenda)", expanded=False):
            st.markdown("""
            **1. Lambda Score (Enrichment):** Razão entre a frequência do termo no seu contexto específico versus a frequência global na ciência. Scores > 5 indicam altíssima especificidade.
            **2. P-Value (Fisher's Exact Test):** Probabilidade estatística da associação não ser aleatória. Valores < 0.05 são significativos.
            **3. Blue Ocean:** Termos com alta relevância global mas NENHUMA ou POUCAS publicações no seu contexto específico (Oportunidade de Ineditismo).
            """)

        # TABELA COM LUPA NA ÚLTIMA COLUNA
        st.subheader("Mapa de Descobertas")
        df_view = df.drop(columns=["_sort"]).copy()
        
        # Move a coluna de investigação para o final para ficar mais intuitivo
        # Streamlit dataframe selection event
        selection = st.dataframe(
            df_view,
            use_container_width=True,
            hide_index=True,
            on_select="rerun",
            selection_mode="single-row"
        )
        
        # GATILHO DO POP-UP
        if selection.selection.rows:
            idx = selection.selection.rows[0]
            alvo_sel = df.iloc[idx][t["col_mol"]]
            # Abre o Modal
            abrir_modal_investigacao(
                alvo_sel, 
                st.session_state.alvo_guardado, 
                st.session_state.email_guardado, 
                st.session_state.api_key_usuario, 
                st.session_state.lang
            )

        # Download CSV
        csv = df.drop(columns=["_sort"]).to_csv(index=False).encode('utf-8')
        st.download_button("📥 Download CSV", csv, "lemos_lambda_results.csv", "text/csv")

else:
    # --- PÁGINA HOME (Igual à anterior) ---
    st.markdown(f'<p class="header-style">{t["titulo_desk"]}</p>', unsafe_allow_html=True)
    st.markdown(f'<p class="sub-header-style">{t["subtitulo"]}</p>', unsafe_allow_html=True)
    
    exibir_radar_cientifico(st.session_state.lang, t)
    st.divider()
    
    col1, col2 = st.columns([2, 1])
    with col1:
        tab1, tab2, tab3 = st.tabs(["🔍 Auto-Miner", "📥 Import List", "🔥 Presets"])
        
        with tab1:
            st.text_input(t["label_email"], key="input_email")
            st.text_input(t["label_alvo"], key="input_alvo")
            if st.button(t["btn_auto"], use_container_width=True, type="primary"):
                carregar_lista_dinamica_smart(t)
        
        with tab2:
            st.write("Importe sua lista (.txt ou .csv)")
            st.file_uploader("Upload", type=['txt', 'csv'], key="uploader_key", on_change=processar_upload, args=(t,))

        with tab3:
            st.write("Presets de Fronteira")
            preset_select = st.selectbox("Selecione Categoria:", list(c.PRESETS_FRONTEIRA.keys()))
            if st.button("Adicionar Preset"):
                qtd = adicionar_termos_seguro(c.PRESETS_FRONTEIRA[preset_select])
                st.success(f"Adicionados {qtd} termos!")

        st.divider()
        if st.session_state.alvos_val:
            st.text_area(t["expander_lista"], key="alvos_val", height=150)
            c_run, c_clear = st.columns([4, 1])
            with c_run:
                if st.button(t["btn_executar"], use_container_width=True, type="primary"):
                    ir_para_analise(st.session_state.input_email, st.session_state.input_fonte, st.session_state.input_alvo, 2015, 2026, t)
            with c_clear:
                st.button("🗑️", on_click=limpar_lista_total)
    
    with col2:
        st.subheader("⚙️ Setup")
        st.session_state.api_key_usuario = st.text_input("Gemini API Key", type="password")
        st.markdown("[🔑 **Gerar API Key Grátis**](https://aistudio.google.com/app/apikey)")
        st.toggle(t["label_ia_global"], key="ia_global_switch")
        st.text_input(t["label_contexto"], key="input_fonte")

# --- RODAPÉ ---
st.markdown("---")
c1, c2 = st.columns([2, 1])
with c1:
    st.caption(t["footer_citar"])
    with st.expander("📚 Citação (Zenodo/Mendeley)"):
        citacao_texto = "Lemos, G. (2025). Lemos Lambda: Deep Science Prospector (v2.0). Zenodo. https://doi.org/10.5281/zenodo.18092141"
        st.code(citacao_texto, language="text")
        st.download_button("📥 Baixar Citação (.bib)", gerar_bibtex(), "lemos_lambda.bib", "application/x-bibtex")

with c2:
    st.caption(t["apoio_titulo"])
    st.text_input("Pix:", value="960f3f16-06ce-4e71-9b5f-6915b2a10b5a")
