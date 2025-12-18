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

# --- ESTADO GERAL ---
if 'pagina' not in st.session_state: st.session_state.pagina = 'home'
if 'alvos_val' not in st.session_state: st.session_state.alvos_val = ""
if 'resultado_df' not in st.session_state: st.session_state.resultado_df = None
if 'news_index' not in st.session_state: st.session_state.news_index = 0
if 'input_alvo' not in st.session_state: st.session_state.input_alvo = ""
if 'input_fonte' not in st.session_state: st.session_state.input_fonte = ""
if 'artigos_detalhe' not in st.session_state: st.session_state.artigos_detalhe = None
# Guarda o email na sessão para não pedir de novo na tela de detalhes
if 'email_guardado' not in st.session_state: st.session_state.email_guardado = ""

# --- IDIOMA ---
lang_opt = st.sidebar.radio("🌐 Language:", ["🇧🇷 PT", "🇺🇸 EN"], horizontal=True)
lang = "pt" if "PT" in lang_opt else "en"
t = c.TEXTOS[lang]

# ==========================================
# HELPERS
# ==========================================
def limpar_campo(chave_session):
    st.session_state[chave_session] = ""

def aplicar_preset_lemos(textos):
    st.session_state.input_alvo = c.PRESET_LEMOS["alvo"]
    st.session_state.input_fonte = c.PRESET_LEMOS["fonte"]
    adicionar_termos_seguro(c.CANDIDATOS_MINERACAO, textos)
    st.toast(textos["toast_preset"], icon="🎓")

def carregar_apenas_biblioteca(textos):
    adicionar_termos_seguro(c.CANDIDATOS_MINERACAO, textos)
    st.toast(textos["toast_lib"], icon="📚")

def adicionar_termos_seguro(novos_termos_lista, textos):
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
    return len(adicionados)

def resetar_pesquisa():
    st.session_state.pagina = 'home'
    st.session_state.resultado_df = None
    st.session_state.artigos_detalhe = None
    st.rerun()

def ir_para_analise(email_user, contexto, alvo, ano_ini, ano_fim):
    """
    Executa a análise com o 'Lemos Score' Otimizado.
    """
    # Guarda email para uso futuro
    st.session_state.email_guardado = email_user
    
    lista = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
    resultados = []
    
    # Placeholder de carregamento (Tela inteira limpa)
    placeholder = st.empty()
    with placeholder.container():
        st.markdown("## 🧬 Lemos Lambda Deep Processing...")
        st.markdown(f"Analisando **{len(lista)} moléculas** contra o alvo **'{alvo}'**...")
        st.info("Cruze os dedos! Procurando 'Blue Oceans' e Pesquisas Embrionárias...")
        prog = st.progress(0)
    
    for i, item in enumerate(lista):
        # Pequeno delay para não bloquear API
        time.sleep(0.05)
        
        termo_contexto = contexto if contexto else None
        
        # Consultas
        n_global = bk.consultar_pubmed_count(item, termo_contexto, email_user, ano_ini, ano_fim)
        n_especifico = bk.consultar_pubmed_count(item, alvo, email_user, ano_ini, ano_fim)
        
        # Cálculo do Ratio Simples
        ratio = n_global / n_especifico if n_especifico > 0 else n_global
        
        # --- LÓGICA OTIMIZADA (LEMOS SCORE v2.0) ---
        # Prioridade: Encontrar o "Sweet Spot" (Embrionário)
        
        if n_especifico == 0:
            if n_global > 50:
                tag = "💎 Blue Ocean (Inexplorado)"
                score_sort = 1000 # Força topo da lista
            else:
                tag = "👻 Fantasma (Sem relevância global)"
                score_sort = 0
        
        elif 1 <= n_especifico <= 15:
            # AQUI ESTÁ A "PESQUISA EMBRIONÁRIA"
            tag = "🌱 Embrionário (Nascendo agora)"
            score_sort = 500 # Segundo no topo
            
        elif ratio > 20:
            tag = "🚀 Tendência Global (Translação Alta)"
            score_sort = 100
            
        elif ratio > 5:
            tag = "🥇 Ouro (Promissor)"
            score_sort = 50
            
        elif ratio < 2:
            tag = "🔴 Saturado (Difícil publicar)"
            score_sort = 10
            
        else:
            tag = "⚖️ Neutro"
            score_sort = 20
        
        resultados.append({
            t["col_mol"]: item, 
            t["col_status"]: tag, 
            t["col_ratio"]: round(ratio, 1),
            t["col_art_alvo"]: n_especifico,
            t["col_global"]: n_global,
            "_sort": score_sort # Coluna oculta para ordenação inteligente
        })
        prog.progress((i+1)/len(lista))
    
    placeholder.empty()
    
    # Ordena pelo Score Lemos (Blue Ocean -> Embrionário -> Tendência -> Resto)
    df_final = pd.DataFrame(resultados).sort_values(by=["_sort", t["col_ratio"]], ascending=[False, False])
    st.session_state.resultado_df = df_final
    
    st.session_state.pagina = 'resultados' # Troca a tela
    st.rerun()

@st.fragment(run_every=60) 
def exibir_radar_cientifico(lang_code, textos):
    news_list = bk.buscar_todas_noticias(lang_code)
    if not news_list: return
    idx = st.session_state.news_index % len(news_list)
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
# ROTEAMENTO DE TELA
# ==========================================

# --- TELA 1: HOME (WIZARD) ---
if st.session_state.pagina == 'home':
    st.title(t["titulo_desk"])
    st.caption(t["subtitulo"])
    
    exibir_radar_cientifico(lang, t)
    st.divider()

    with st.container(border=True):
        col_input, col_config = st.columns([2, 1])
        
        with col_input:
            st.subheader(t["step_1"])
            st.warning(t["aviso_pubmed"])
            # Usa o email guardado se existir
            email_val = st.session_state.email_guardado if st.session_state.email_guardado else ""
            email_user = st.text_input(t["label_email"], value=email_val, placeholder=t["holder_email"])
            
            c_in, c_bin = st.columns([8, 1], vertical_alignment="bottom")
            with c_in: alvo = st.text_input(t["label_alvo"], key="input_alvo", placeholder=t["holder_alvo"])
            with c_bin: st.button("🗑️", key="lixo_alvo", on_click=limpar_campo, args=("input_alvo",))

            c_magic, c_manual = st.columns(2)
            with c_magic:
                if st.button(t["btn_magic"], type="primary", use_container_width=True):
                    if not email_user or not alvo: st.error(t["erro_campos"])
                    else:
                        with st.status(t["prog_magic"], expanded=True) as status:
                            st.write(t["status_minerando"])
                            novos = bk.buscar_alvos_emergentes_pubmed(alvo, email_user)
                            st.write(t["status_filtrando"])
                            count = adicionar_termos_seguro(c.CANDIDATOS_MINERACAO + novos, t)
                            status.update(label=t["status_pronto"], state="complete", expanded=False)
                            st.toast(f"✅ {len(novos)} novos!", icon="🧬")

            with c_manual:
                with st.popover(t["label_manual"], use_container_width=True):
                    termo_man = st.text_input("Termo", key="input_manual", placeholder=t["holder_manual"])
                    if st.button(t["btn_add_manual"], use_container_width=True):
                        if termo_man:
                            adicionar_termos_seguro([x.strip() for x in termo_man.split(",")], t)
                            st.session_state.input_manual = ""; st.rerun()

        with col_config:
            st.subheader("Config")
            c_btn1, c_btn2 = st.columns(2)
            with c_btn1: st.button(t["btn_lib"], on_click=carregar_apenas_biblioteca, args=(t,), use_container_width=True)
            with c_btn2: st.button(t["btn_preset"], on_click=aplicar_preset_lemos, args=(t,), use_container_width=True)
                
            st.subheader(t["label_periodo"])
            anos_range = st.slider("Anos", 2000, datetime.now().year, (2015, datetime.now().year), label_visibility="collapsed")
            
            st.markdown("---")
            c_src, c_sbin = st.columns([8, 1], vertical_alignment="bottom")
            with c_src: contexto = st.text_input("Context", key="input_fonte", placeholder=t["holder_fonte"], label_visibility="collapsed")
            with c_sbin: st.button("🗑️", key="lixo_fonte", on_click=limpar_campo, args=("input_fonte",))
            
            st.markdown("---")
            st.file_uploader(t["desc_import"], type=["csv", "txt"], key="uploader_key", on_change=processar_upload, args=(t,), label_visibility="collapsed")

    st.divider()
    if st.session_state.alvos_val:
        with st.expander(t["ver_editar"], expanded=False):
            c1, c2 = st.columns([5,1])
            with c1: st.session_state.alvos_val = st.text_area("Termos", value=st.session_state.alvos_val, height=100)
            with c2: 
                if st.button(t["btn_limpar_tudo"]): st.session_state.alvos_val = ""
        
        # BOTÃO PRINCIPAL DE AÇÃO
        if st.button(t["analise_btn"], type="primary", use_container_width=True):
            if not email_user: st.error(t["erro_email"])
            else:
                ir_para_analise(email_user, contexto, alvo, anos_range[0], anos_range[1])

# --- TELA 2: RESULTADOS (DASHBOARD) ---
elif st.session_state.pagina == 'resultados':
    
    # Header de Navegação
    c_back, c_tit = st.columns([1, 5])
    with c_back:
        st.button("⬅️ Nova Pesquisa", on_click=resetar_pesquisa, use_container_width=True)
    with c_tit:
        st.title(t["resultados"])
    
    df = st.session_state.resultado_df
    if df is not None and not df.empty:
        # Pega o melhor termo (ignorando colunas ocultas)
        top = df.iloc[0]
        
        # KPIS
        c1, c2, c3 = st.columns(3)
        c1.metric(t["metrica_potencial"], top[t["col_mol"]], delta=top[t["col_status"]])
        c2.metric(t["metrica_score"], top[t["col_ratio"]])
        c3.metric(t["metrica_artigos"], top[t["col_art_alvo"]])
        
        # CHART
        st.subheader("Mapa de Oportunidades")
        # Filtra o dataframe para exibição (remove coluna de sort interna)
        df_show = df.drop(columns=["_sort"])
        
        fig = px.bar(df_show.head(25), x=t["col_mol"], y=t["col_ratio"], color=t["col_status"],
                     color_discrete_map={
                         "💎 Blue Ocean (Inexplorado)": "#00CC96", 
                         "🌱 Embrionário (Nascendo agora)": "#00FF00", # Verde brilhante para destaque
                         "🚀 Tendência Global (Translação Alta)": "#AB63FA",
                         "🥇 Ouro (Promissor)": "#636EFA", 
                         "🔴 Saturado (Difícil publicar)": "#EF553B",
                         "👻 Fantasma (Sem relevância global)": "#808080"
                     })
        st.plotly_chart(fig, use_container_width=True)
        
        # TABELA
        st.dataframe(df_show, use_container_width=True, hide_index=True)
        st.download_button(t["btn_baixar"], df_show.to_csv(index=False).encode('utf-8'), "lemos_lambda_report.csv", "text/csv")
        
        # DETALHES DOS ARTIGOS
        st.divider()
        st.subheader("📄 Leitura Profunda: Investigação de Papers")
        st.info("Selecione um alvo acima para ver os artigos reais no PubMed com tradução das conclusões.")
        
        termos_disp = df[t["col_mol"]].unique().tolist()
        # Seleciona o primeiro da lista (o melhor) por padrão
        sel_mol = st.selectbox("Selecione o alvo:", termos_disp, index=0)
        
        if st.button(f"🔍 Carregar Artigos sobre: {sel_mol}", type="secondary"):
            with st.spinner(f"Buscando literatura sobre {sel_mol}..."):
                # Recupera o alvo original do input que está guardado na session
                alvo_analise = st.session_state.input_alvo
                email_analise = st.session_state.email_guardado
                
                if not alvo_analise or not email_analise:
                    st.warning("Dados da sessão perdidos. Por favor reinicie a pesquisa.")
                else:
                    arts = bk.buscar_resumos_detalhados(sel_mol, alvo_analise, email_analise, 2015, 2025, lang)
                    st.session_state.artigos_detalhe = arts
        
        if st.session_state.artigos_detalhe:
            st.markdown(f"### Artigos encontrados: {len(st.session_state.artigos_detalhe)}")
            if not st.session_state.artigos_detalhe:
                st.warning("Nenhum artigo encontrado com resumo disponível neste período.")
            
            for art in st.session_state.artigos_detalhe:
                with st.expander(f"{art['Title']}"):
                    st.write(f"**Fonte:** {art['Source']} | **PMID:** {art['PMID']}")
                    st.info(f"**Conclusão/Resumo (Traduzido):**\n\n{art['Resumo_IA']}")
                    st.link_button("Ler no PubMed 🔗", art['Link'])

# RODAPÉ COMUM
st.markdown("---"); st.caption(f"© 2025 Guilherme Lemos | {t['footer_citar']}")

# SIDEBAR COMUM
st.sidebar.markdown("---")
with st.sidebar.expander(t["citar_titulo"], expanded=True):
    st.code(t["citar_texto"], language="text")
    st.link_button(t["link_doi"], "https://doi.org/10.5281/zenodo.17958507")