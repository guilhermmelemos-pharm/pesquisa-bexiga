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
import constantes as c
import backend as bk

st.set_page_config(page_title="Lemos Lambda", page_icon="λ", layout="wide")

# --- CSS INJECTION PARA FAZER O BLUE OCEAN BRILHAR ---
st.markdown("""
    <style>
    /* Estilo Geral dos Botões */
    .stButton button { 
        border-radius: 12px; 
        height: 50px; 
        font-weight: bold; 
    }
    
    /* Animação de Pulso para o Blue Ocean */
    @keyframes pulse-blue {
        0% { box-shadow: 0 0 0 0 rgba(0, 204, 150, 0.7); transform: scale(1); }
        70% { box-shadow: 0 0 0 10px rgba(0, 204, 150, 0); transform: scale(1.02); }
        100% { box-shadow: 0 0 0 0 rgba(0, 204, 150, 0); transform: scale(1); }
    }
    
    /* Classe específica para o container do Blue Ocean */
    .blue-ocean-btn button {
        animation: pulse-blue 2s infinite;
        background: linear-gradient(90deg, #00C9FF 0%, #92FE9D 100%) !important;
        color: #004d40 !important;
        border: none !important;
        font-size: 1.1rem !important;
        text-transform: uppercase;
        letter-spacing: 1px;
    }
    
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
if 'input_email' not in st.session_state: st.session_state.input_email = ""
if 'artigos_detalhe' not in st.session_state: st.session_state.artigos_detalhe = None
if 'email_guardado' not in st.session_state: st.session_state.email_guardado = ""
if 'alvo_guardado' not in st.session_state: st.session_state.alvo_guardado = ""

lang_opt = st.sidebar.radio("🌐 Language:", ["🇧🇷 PT", "🇺🇸 EN"], horizontal=True)
lang = "pt" if "PT" in lang_opt else "en"
t = c.TEXTOS[lang]

# ==========================================
# LÓGICA DE NEGÓCIO
# ==========================================
def limpar_campo(chave_session):
    st.session_state[chave_session] = ""

def limpar_lista_total():
    st.session_state.alvos_val = ""

def carregar_lista_dinamica_smart(textos):
    email = st.session_state.input_email
    alvo = st.session_state.input_alvo
    
    if not alvo:
        st.error("⚠️ Preencha o campo 'Alvo Principal' (ex: Liver, Kidney) antes de buscar!")
        return

    existentes = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
    lista_mestra = list(set(existentes + c.CANDIDATOS_MINERACAO))
    
    msg_final = textos["msg_sucesso_base"]
    novos_encontrados = 0
    
    if alvo and email:
        with st.spinner(f"{textos['status_minerando']} {alvo}..."):
            novos = bk.buscar_alvos_emergentes_pubmed(alvo, email)
            if novos:
                lista_mestra.extend(novos)
                novos_encontrados = len(novos)
                msg_final = textos["msg_sucesso_dinamico"].format(qtd=novos_encontrados)
    
    adicionar_termos_seguro(lista_mestra, textos)
    st.toast(msg_final, icon="🧬")

def explorar_blue_ocean(textos):
    email = st.session_state.input_email
    alvo = st.session_state.input_alvo

    if not email or not alvo:
        st.error(textos["erro_campos"])
        return

    with st.spinner(textos["status_blue_ocean"]):
        novos = bk.buscar_alvos_emergentes_pubmed(alvo, email)
        if novos:
            count = adicionar_termos_seguro(novos, textos)
            st.success(textos["msg_sucesso_blue"].format(qtd=count))
        else:
            st.warning("Nenhum termo novo encontrado.")

def minerar_novidades_fonte(textos):
    fonte = st.session_state.input_fonte
    email = st.session_state.input_email
    if not fonte:
        st.error(textos["erro_fonte_vazia"])
        return
    if not email:
        st.error(textos["erro_email"])
        return

    with st.spinner(f"Minerando: {fonte}..."):
        novos_termos = bk.buscar_alvos_emergentes_pubmed(fonte, email)
        if novos_termos:
            count = adicionar_termos_seguro(novos_termos, textos)
            st.success(f"✅ {count} novos termos em '{fonte}'.")

def aplicar_preset_lemos(textos):
    st.session_state.input_alvo = c.PRESET_LEMOS["alvo"]
    st.session_state.input_fonte = c.PRESET_LEMOS["fonte"]
    carregar_lista_dinamica_smart(textos)

def adicionar_termos_seguro(novos_termos_lista, textos):
    atuais = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
    atuais_upper = [x.upper() for x in atuais]
    adicionados = []
    
    blacklist_lower = [x.lower() for x in c.BLACKLIST_GERAL]
    
    for termo in novos_termos_lista:
        t_limpo = termo.strip()
        t_lower = t_limpo.lower()
        if any(bad in t_lower for bad in blacklist_lower): continue
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

def ir_para_analise(email_user, contexto, alvo, ano_ini, ano_fim):
    st.session_state.email_guardado = email_user
    st.session_state.alvo_guardado = alvo
    
    lista = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
    resultados = []
    
    placeholder = st.empty()
    with placeholder.container():
        st.markdown("## 🧬 Lemos Lambda Deep Processing...")
        st.markdown(f"Analisando **{len(lista)} moléculas** contra o alvo **'{alvo}'**...")
        prog = st.progress(0)
    
    for i, item in enumerate(lista):
        time.sleep(0.05)
        termo_contexto = contexto if contexto else None
        n_global = bk.consultar_pubmed_count(item, termo_contexto, email_user, ano_ini, ano_fim)
        n_especifico = bk.consultar_pubmed_count(item, alvo, email_user, ano_ini, ano_fim)
        ratio = n_global / n_especifico if n_especifico > 0 else n_global
        
        if n_especifico == 0:
            if n_global > 50: tag, score_sort = "💎 Blue Ocean (Inexplorado)", 1000
            else: tag, score_sort = "👻 Fantasma (Sem relevância)", 0
        elif 1 <= n_especifico <= 15: tag, score_sort = "🌱 Embrionário (Nascendo agora)", 500
        elif ratio > 20: tag, score_sort = "🚀 Tendência (Translação)", 100
        elif ratio > 5: tag, score_sort = "🥇 Ouro", 50
        elif ratio < 2: tag, score_sort = "🔴 Saturado", 10
        else: tag, score_sort = "⚖️ Neutro", 20
        
        # FIX: Usamos chaves internas FIXAS (termo, status, ratio) para evitar KeyError
        resultados.append({
            "termo": item,
            "status": tag,
            "ratio": round(ratio, 1),
            "alvo_count": n_especifico,
            "fonte_count": n_global,
            "_sort": score_sort
        })
        prog.progress((i+1)/len(lista))
    
    placeholder.empty()
    df_final = pd.DataFrame(resultados).sort_values(by=["_sort", "ratio"], ascending=[False, False])
    st.session_state.resultado_df = df_final
    st.session_state.pagina = 'resultados'
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

# --- UI ---
if st.session_state.pagina == 'home':
    st.title(t["titulo_desk"]); st.caption(t["subtitulo"])
    exibir_radar_cientifico(lang, t)
    st.divider()
    with st.container(border=True):
        c1, c2 = st.columns([2, 1])
        with c1:
            st.subheader(t["step_1"]); st.warning(t["aviso_pubmed"])
            st.text_input(t["label_email"], key="input_email", placeholder=t["holder_email"])
            ca, cb = st.columns([8, 1], vertical_alignment="bottom")
            with ca: st.text_input(t["label_alvo"], key="input_alvo", placeholder=t["holder_alvo"])
            with cb: st.button(t["btn_limpar"], key="lixo_alvo", on_click=limpar_campo, args=("input_alvo",))

            st.write(" ")
            b_smart, b_preset = st.columns(2)
            with b_smart:
                st.button(t["btn_smart_load"], type="primary", on_click=carregar_lista_dinamica_smart, args=(t,), use_container_width=True)
            with b_preset:
                st.button(t["btn_preset"], type="secondary", on_click=aplicar_preset_lemos, args=(t,), use_container_width=True)
            
            st.write(" ")
            st.markdown('<div class="blue-ocean-btn">', unsafe_allow_html=True)
            if st.button(t["btn_blue_ocean"], on_click=explorar_blue_ocean, args=(t,), use_container_width=True):
                pass 
            st.markdown('</div>', unsafe_allow_html=True)
            
            st.write(" ")
            st.button(t["btn_lib"], on_click=minerar_novidades_fonte, args=(t,), use_container_width=True)

            with st.popover(t["label_manual"], use_container_width=True):
                termo_man = st.text_input("Termo", key="input_manual", placeholder=t["holder_manual"])
                if st.button(t["btn_add_manual"], use_container_width=True):
                    if termo_man:
                        adicionar_termos_seguro([x.strip() for x in termo_man.split(",")], t)
                        st.session_state.input_manual = ""; st.rerun()

        with c2:
            st.subheader("Config"); st.subheader(t["label_periodo"])
            anos_range = st.slider("Anos", 2000, datetime.now().year, (2015, datetime.now().year), label_visibility="collapsed")
            st.markdown("---"); st.caption(t["label_fonte"])
            ca, cb = st.columns([8, 1], vertical_alignment="bottom")
            with ca: st.text_input("Context", key="input_fonte", placeholder=t["holder_fonte"], label_visibility="collapsed")
            with cb: st.button(t["btn_limpar"], key="lixo_fonte", on_click=limpar_campo, args=("input_fonte",))
            st.markdown("---")
            st.file_uploader(t["desc_import"], type=["csv", "txt"], key="uploader_key", on_change=processar_upload, args=(t,), label_visibility="collapsed")

    st.divider()
    if st.session_state.alvos_val:
        with st.expander(t["ver_editar"], expanded=False):
            c1, c2 = st.columns([5,1])
            with c1: st.session_state.alvos_val = st.text_area("Termos", value=st.session_state.alvos_val, height=100)
            with c2: 
                st.button(t["btn_limpar_tudo"], on_click=limpar_lista_total)
                lista_txt = st.session_state.alvos_val.replace(", ", "\n").replace(",", "\n")
                st.download_button(t["btn_export_lista"], lista_txt, "lemos_lambda_list.csv", "text/csv")
        if st.button(t["analise_btn"], type="primary", use_container_width=True):
            email = st.session_state.input_email
            if not email: st.error(t["erro_email"])
            else: ir_para_analise(email, st.session_state.input_fonte, st.session_state.input_alvo, anos_range[0], anos_range[1])

elif st.session_state.pagina == 'resultados':
    c_back, c_tit = st.columns([1, 5])
    with c_back: st.button(t["btn_nova_pesquisa"], on_click=resetar_pesquisa, use_container_width=True)
    with c_tit: st.title(t["resultados"])
    
    df = st.session_state.resultado_df
    if df is not None and not df.empty:
        # FIX: Acessa com chaves FIXAS ("termo", "status"), o que resolve o KeyError
        top = df.iloc[0]
        c1, c2, c3 = st.columns(3)
        c1.metric(t["metrica_potencial"], top["termo"], delta=top["status"])
        c2.metric(t["metrica_score"], top["ratio"])
        c3.metric(t["metrica_artigos"], top["alvo_count"])
        
        st.subheader(t["titulo_mapa"])
        
        # Renomear colunas apenas para exibição
        df_show = df.rename(columns={
            "termo": t["col_mol"],
            "status": t["col_status"],
            "ratio": t["col_ratio"],
            "alvo_count": t["col_art_alvo"],
            "fonte_count": t["col_global"]
        }).drop(columns=["_sort"])
        
        # Gráfico usa as colunas renomeadas
        fig = px.bar(df_show.head(25), x=t["col_mol"], y=t["col_ratio"], color=t["col_status"], 
                     color_discrete_map={"💎 Blue Ocean (Inexplorado)": "#00CC96", "🌱 Embrionário (Nascendo agora)": "#00FF00", "🚀 Tendência (Translação)": "#AB63FA", "🥇 Ouro": "#636EFA", "🔴 Saturado": "#EF553B", "👻 Fantasma (Sem relevância)": "#808080"})
        st.plotly_chart(fig, use_container_width=True)
        st.dataframe(df_show, use_container_width=True, hide_index=True)
        st.download_button(t["btn_baixar"], df_show.to_csv(index=False).encode('utf-8'), "lemos_lambda_report.csv", "text/csv")
        
        st.divider(); st.subheader(t["titulo_leitura"]); st.info(t["info_leitura"])
        # Selectbox usa chave fixa "termo"
        termos_disp = sorted(df["termo"].unique().tolist())
        sel_mol = st.selectbox(t["sel_leitura"], termos_disp, index=0)
        
        if st.button(f"{t['btn_buscar_artigos']} {sel_mol}", type="secondary"):
            with st.spinner(f"{t['msg_buscando_lit']} {sel_mol}..."):
                alvo_analise = st.session_state.get('alvo_guardado', '')
                email_analise = st.session_state.get('email_guardado', '')
                if not alvo_analise or not email_analise: st.warning(t["erro_sessao"])
                else: st.session_state.artigos_detalhe = bk.buscar_resumos_detalhados(sel_mol, alvo_analise, email_analise, 2015, 2025, lang)
        
        if st.session_state.artigos_detalhe:
            st.markdown(f"### {t['header_artigos_enc']} {len(st.session_state.artigos_detalhe)}")
            if not st.session_state.artigos_detalhe: st.warning(t["aviso_sem_artigos"])
            for art in st.session_state.artigos_detalhe:
                with st.expander(f"{art['Title']}"):
                    st.info(f"**Abstract/Conclusão:**\n\n{art['Resumo_IA']}")
                    st.link_button("PubMed 🔗", art['Link'])

st.markdown("---"); st.caption(f"© 2025 Guilherme Lemos | {t['footer_citar']}")
st.sidebar.markdown("---")
with st.sidebar.expander(t["citar_titulo"], expanded=True):
    st.code(t["citar_texto"], language="text"); st.link_button(t["link_doi"], "https://doi.org/10.5281/zenodo.17958507")