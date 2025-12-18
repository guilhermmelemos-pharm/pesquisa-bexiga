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
Version: 1.7.2 (Statistical Update)
"""
import streamlit as st

# --- 1. CONFIGURAÇÃO (PRIMEIRA LINHA OBRIGATÓRIA) ---
st.set_page_config(page_title="Lemos Lambda", page_icon="λ", layout="wide")

import pandas as pd
import plotly.express as px
from datetime import datetime
import time
import math
import constantes as c
import backend as bk

# --- 2. INICIALIZAÇÃO SEGURA DE ESTADO (ANTI-TRAVAMENTO) ---
# Define valores padrão para evitar erros de "KeyError" ou telas vazias
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
    'alvo_guardado': ""
}

for key, val in defaults.items():
    if key not in st.session_state:
        st.session_state[key] = val

# --- 3. CSS INJECTION ---
st.markdown("""
    <style>
    .stButton button { border-radius: 12px; height: 50px; font-weight: bold; }
    
    @keyframes pulse-blue {
        0% { box-shadow: 0 0 0 0 rgba(0, 204, 150, 0.7); transform: scale(1); }
        70% { box-shadow: 0 0 0 10px rgba(0, 204, 150, 0); transform: scale(1.02); }
        100% { box-shadow: 0 0 0 0 rgba(0, 204, 150, 0); transform: scale(1); }
    }
    
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
    
    /* CSS controlando tamanho da imagem para evitar erro do Python */
    div[data-testid="stImage"] img { 
        height: 160px !important; 
        object-fit: cover !important; 
        border-radius: 10px !important; 
        width: 100% !important;
    }
    .stAlert { padding: 0.5rem; margin-bottom: 1rem; border-radius: 8px; }
    </style>
""", unsafe_allow_html=True)

lang_opt = st.sidebar.radio("🌐 Language:", ["🇧🇷 PT", "🇺🇸 EN"], horizontal=True)
lang = "pt" if "PT" in lang_opt else "en"
t = c.TEXTOS[lang]

# ================= LÓGICA DE NEGÓCIO =================

def limpar_campo(k): st.session_state[k] = ""
def limpar_lista_total(): st.session_state.alvos_val = ""

def carregar_lista_dinamica_smart(textos):
    email, alvo = st.session_state.input_email, st.session_state.input_alvo
    if not alvo:
        st.error("⚠️ Defina o Alvo Principal (ex: Liver, Kidney)!")
        return
    existentes = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
    lista_mestra = list(set(existentes + c.CANDIDATOS_MINERACAO))
    
    msg = textos["msg_sucesso_base"]
    if alvo and email:
        with st.spinner(f"{textos['status_minerando']} {alvo}..."):
            novos = bk.buscar_alvos_emergentes_pubmed(alvo, email)
            if novos:
                lista_mestra.extend(novos)
                msg = textos["msg_sucesso_dinamico"].format(qtd=len(novos))
    
    adicionar_termos_seguro(lista_mestra, textos)
    st.toast(msg, icon="🧬")

def explorar_blue_ocean(textos):
    email, alvo = st.session_state.input_email, st.session_state.input_alvo
    if not email or not alvo:
        st.error(textos["erro_campos"])
        return
    with st.spinner(textos["status_blue_ocean"]):
        novos = bk.buscar_alvos_emergentes_pubmed(alvo, email)
        if novos:
            c = adicionar_termos_seguro(novos, textos)
            st.success(textos["msg_sucesso_blue"].format(qtd=c))
        else:
            st.warning("Nada novo encontrado.")

def minerar_novidades_fonte(textos):
    fonte, email = st.session_state.input_fonte, st.session_state.input_email
    if not fonte or not email:
        st.error("Preencha e-mail e contexto.")
        return
    with st.spinner(f"Minerando {fonte}..."):
        novos = bk.buscar_alvos_emergentes_pubmed(fonte, email)
        if novos:
            c = adicionar_termos_seguro(novos, textos)
            st.success(f"✅ {c} novos termos.")

def aplicar_preset_lemos(textos):
    st.session_state.input_alvo = c.PRESET_LEMOS["alvo"]
    st.session_state.input_fonte = c.PRESET_LEMOS["fonte"]
    carregar_lista_dinamica_smart(textos)

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

def resetar_pesquisa():
    st.session_state.pagina = 'home'
    st.session_state.resultado_df = None
    st.session_state.artigos_detalhe = None

def ir_para_analise(email, contexto, alvo, y_ini, y_fim):
    st.session_state.email_guardado = email
    st.session_state.alvo_guardado = alvo
    lista = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
    res = []
    
    # --- CALIBRAÇÃO ESTATÍSTICA (Z-SCORE BASELINE) ---
    with st.spinner("Calibrando estatística..."):
        try:
            n_total_alvo = bk.consultar_pubmed_count(alvo, "", email, 1900, 2030)
            if n_total_alvo == 0: n_total_alvo = 1
        except:
            n_total_alvo = 100000 
            
    N_PUBMED = 36000000 
    
    placeholder = st.empty()
    with placeholder.container():
        st.markdown("## 🧬 Lemos Lambda Processing...")
        st.markdown(f"Analisando **{len(lista)} alvos**...")
        prog = st.progress(0)
    
    for i, item in enumerate(lista):
        time.sleep(0.05)
        termo_ctx = contexto if contexto else None
        
        n_global = bk.consultar_pubmed_count(item, termo_ctx, email, y_ini, y_fim)
        n_especifico = bk.consultar_pubmed_count(item, alvo, email, y_ini, y_fim)
        
        # --- CÁLCULO ESTATÍSTICO (Enrichment Score) ---
        expected = (n_global * n_total_alvo) / N_PUBMED
        if expected == 0: expected = 0.00001
        enrichment = (n_especifico + 0.1) / expected
        
        tag = "⚖️ Neutro"
        score_sort = 0
        
        # LÓGICA DE CLASSIFICAÇÃO
        if n_especifico == 0:
            if n_global > 100: 
                tag = "💎 Blue Ocean (Inexplorado)"
                score_sort = 1000
            else:
                tag = "👻 Fantasma"
                score_sort = 0
        
        elif n_especifico <= 15:
            # Resgate de emergentes
            if n_global > 50:
                tag = "🌱 Embrionário (Emergente)"
                score_sort = 500
            else:
                tag = "⚖️ Neutro"
                score_sort = 20
        
        else:
            # Estatística robusta
            if enrichment > 5: 
                tag = "🥇 Ouro (Alta Conexão)"
                score_sort = 100
            elif enrichment > 1.5:
                tag = "🚀 Tendência"
                score_sort = 200
            else:
                tag = "🔴 Saturado"
                score_sort = 10

        ratio = float(enrichment)
        
        res.append({
            t["col_mol"]: item, t["col_status"]: tag, t["col_ratio"]: round(ratio, 1), 
            t["col_art_alvo"]: n_especifico, t["col_global"]: n_global, "_sort": score_sort
        })
        prog.progress((i+1)/len(lista))
    
    placeholder.empty()
    df_final = pd.DataFrame(res).sort_values(by=["_sort", t["col_ratio"]], ascending=[False, False])
    st.session_state.resultado_df = df_final
    st.session_state.pagina = 'resultados'
    st.rerun()

@st.fragment(run_every=60) 
def exibir_radar_cientifico(lang_code, textos):
    news = bk.buscar_todas_noticias(lang_code)
    if not news: return
    idx = st.session_state.news_index % len(news)
    batch = news[idx:idx+3]
    st.session_state.news_index += 3
    with st.container(border=True):
        st.caption(textos["radar_titulo"])
        cols = st.columns(3)
        for i, n in enumerate(batch):
            with cols[i]:
                # IMAGEM SEGURA (CSS trata o tamanho)
                if n.get('img'): st.image(n['img'])
                st.markdown(f"**{n['titulo'][:75]}...**")
                st.caption(f"{n.get('bandeira','')} {n.get('fonte','')}")
                st.link_button(textos["btn_ler_feed"], n['link'], use_container_width=True)

def processar_upload(textos):
    f = st.session_state.get('uploader_key')
    if f:
        try:
            c = f.getvalue().decode("utf-8")
            l = [x.strip() for x in c.replace("\n", ",").split(",") if x.strip()]
            n = adicionar_termos_seguro(l, textos)
            st.toast(f"{textos['toast_import']} ({n})", icon="📂")
        except: st.error(textos["erro_ler"])

# --- SISTEMA DE NAVEGAÇÃO À PROVA DE FALHAS ---
# Verifica Explicitamente qual página carregar
if st.session_state.pagina == 'resultados':
    # --- PÁGINA DE RESULTADOS ---
    c_back, c_tit = st.columns([1, 5])
    with c_back: st.button(t["btn_nova_pesquisa"], on_click=resetar_pesquisa, use_container_width=True)
    with c_tit: st.title(t["resultados"])
    
    df = st.session_state.resultado_df
    if df is not None and not df.empty:
        top = df.iloc[0]
        c1, c2, c3 = st.columns(3)
        c1.metric("Top Potencial", top[t["col_mol"]], delta=top[t["col_status"]])
        c2.metric("Enrichment Score", top[t["col_ratio"]])
        c3.metric("Artigos Locais", top[t["col_art_alvo"]])
        
        st.subheader(t["titulo_mapa"])
        df_s = df.drop(columns=["_sort"])
        fig = px.bar(df_s.head(25), x=t["col_mol"], y=t["col_ratio"], color=t["col_status"], 
                     color_discrete_map={
                         "💎 Blue Ocean (Inexplorado)": "#00CC96", 
                         "🌱 Embrionário (Emergente)": "#00FF00", 
                         "🚀 Tendência": "#AB63FA", 
                         "🥇 Ouro (Alta Conexão)": "#636EFA", 
                         "🔴 Saturado": "#EF553B", 
                         "👻 Fantasma": "#808080",
                         "⚖️ Neutro": "#D3D3D3"
                     })
        st.plotly_chart(fig, use_container_width=True)
        st.dataframe(df_s, use_container_width=True, hide_index=True)
        st.download_button(t["btn_baixar"], df_s.to_csv(index=False).encode('utf-8'), "lemos_lambda_report.csv", "text/csv")
        
        st.divider(); st.subheader(t["titulo_leitura"]); st.info(t["info_leitura"])
        t_list = sorted(df[t["col_mol"]].unique().tolist())
        sel = st.selectbox(t["sel_leitura"], t_list, index=0)
        
        if st.button(f"{t['btn_buscar_artigos']} {sel}", type="secondary"):
            with st.spinner(f"{t['msg_buscando_lit']} {sel}..."):
                aa, ee = st.session_state.alvo_guardado, st.session_state.email_guardado
                if not aa or not ee: st.warning(t["erro_sessao"])
                else: st.session_state.artigos_detalhe = bk.buscar_resumos_detalhados(sel, aa, ee, 2015, 2025, lang)
        
        if st.session_state.artigos_detalhe:
            st.markdown(f"### {t['header_artigos_enc']} {len(st.session_state.artigos_detalhe)}")
            if not st.session_state.artigos_detalhe: st.warning(t["aviso_sem_artigos"])
            for a in st.session_state.artigos_detalhe:
                with st.expander(f"{a['Title']}"):
                    st.info(f"**Abstract:**\n\n{a['Resumo_IA']}")
                    st.link_button("PubMed 🔗", a['Link'])

else:
    # --- PÁGINA HOME (DEFAULT FALLBACK) ---
    # Qualquer estado inválido cai aqui, impedindo a tela preta.
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
            with cb: st.button(t["btn_limpar"], key="lixo_alvo", on_click=lambda: limpar_campo("input_alvo"))
            
            st.write(" ")
            b_smart, b_preset = st.columns(2)
            with b_smart: st.button(t["btn_smart_load"], type="primary", on_click=carregar_lista_dinamica_smart, args=(t,), use_container_width=True)
            with b_preset: st.button(t["btn_preset"], type="secondary", on_click=aplicar_preset_lemos, args=(t,), use_container_width=True)
            
            st.write(" ")
            st.markdown('<div class="blue-ocean-btn">', unsafe_allow_html=True)
            if st.button(t["btn_blue_ocean"], on_click=explorar_blue_ocean, args=(t,), use_container_width=True): pass 
            st.markdown('</div>', unsafe_allow_html=True)
            
            st.write(" ")
            st.button(t["btn_lib"], on_click=minerar_novidades_fonte, args=(t,), use_container_width=True)
            with st.popover(t["label_manual"], use_container_width=True):
                tm = st.text_input("Termo", key="input_manual", placeholder=t["holder_manual"])
                if st.button(t["btn_add_manual"], use_container_width=True):
                    if tm: 
                        adicionar_termos_seguro([x.strip() for x in tm.split(",")], t)
                        st.session_state.input_manual = ""; st.rerun()

        with c2:
            st.subheader("Config"); st.subheader(t["label_periodo"])
            anos = st.slider("Anos", 2000, datetime.now().year, (2015, datetime.now().year), label_visibility="collapsed")
            st.markdown("---"); st.caption(t["label_fonte"])
            ca, cb = st.columns([8, 1], vertical_alignment="bottom")
            with ca: st.text_input("Context", key="input_fonte", placeholder=t["holder_fonte"], label_visibility="collapsed")
            with cb: st.button(t["btn_limpar"], key="lixo_fonte", on_click=lambda: limpar_campo("input_fonte"))
            st.markdown("---"); st.file_uploader(t["desc_import"], type=["csv", "txt"], key="uploader_key", on_change=processar_upload, args=(t,), label_visibility="collapsed")

    st.divider()
    if st.session_state.alvos_val:
        with st.expander(t["ver_editar"], expanded=False):
            c1, c2 = st.columns([5,1])
            with c1: st.session_state.alvos_val = st.text_area("Termos", value=st.session_state.alvos_val, height=100)
            with c2: 
                st.button(t["btn_limpar_tudo"], on_click=limpar_lista_total)
                l_txt = st.session_state.alvos_val.replace(", ", "\n").replace(",", "\n")
                st.download_button(t["btn_export_lista"], l_txt, "lemos_lambda_list.csv", "text/csv")
        if st.button(t["analise_btn"], type="primary", use_container_width=True):
            if not st.session_state.input_email: st.error(t["erro_email"])
            else: ir_para_analise(st.session_state.input_email, st.session_state.input_fonte, st.session_state.input_alvo, anos[0], anos[1])

st.markdown("---"); st.caption(f"© 2025 Guilherme Lemos | {t['footer_citar']}")
st.sidebar.markdown("---")
with st.sidebar.expander(t["citar_titulo"], expanded=True):
    st.code(t["citar_texto"], language="text"); st.link_button(t["link_doi"], "https://doi.org/10.5281/zenodo.17958507")
