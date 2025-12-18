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
import time
import plotly.express as px
from datetime import datetime
from Bio import Entrez 

# ==========================================
# IMPORTAÇÃO MODULAR
# ==========================================
# Certifique-se que constantes.py e backend.py estão na mesma pasta
import constantes as c
import backend as bk

# ==========================================
# CONFIGURAÇÃO GLOBAL
# ==========================================
st.set_page_config(page_title="Lemos Lambda", page_icon="λ", layout="wide")

# CSS: Estilo
st.markdown("""
    <style>
    div[data-testid="stImage"] img { height: 150px !important; object-fit: cover !important; border-radius: 8px !important; }
    .stButton button { width: 100%; }
    div[data-testid="stVerticalBlock"] > div { gap: 0.5rem; }
    </style>
""", unsafe_allow_html=True)

# Inicialização do Session State
if 'alvos_val' not in st.session_state: st.session_state.alvos_val = ""
if 'fonte_val' not in st.session_state: st.session_state.fonte_val = ""
if 'alvo_val' not in st.session_state: st.session_state.alvo_val = ""
if 'news_index' not in st.session_state: st.session_state.news_index = 0

# ==========================================
# UI HELPERS (Funções de Interface)
# ==========================================
@st.fragment(run_every=60) 
def exibir_radar_cientifico(lang_code):
    news_list = bk.buscar_todas_noticias(lang_code)
    if not news_list: return
    total_news = len(news_list)
    idx = st.session_state.news_index % total_news
    batch = news_list[idx:idx+3]
    st.session_state.news_index += 3
    with st.container(border=True):
        st.caption(f"📡 **Radar Científico**")
        cols = st.columns(3)
        for i, n in enumerate(batch):
            with cols[i]:
                st.image(n['img'], use_container_width=True)
                st.markdown(f"**{n['titulo'][:60]}...**")
                st.caption(f"{n['bandeira']} {n['fonte']}")
                st.link_button("Ler" if lang_code=='pt' else "Read", n['link'], use_container_width=True)

def carregar_setup_lemos(t):
    st.session_state.alvos_val = ", ".join(c.CANDIDATOS_MINERACAO)
    st.session_state.fonte_val = "Brain OR Kidney OR Liver OR Intestine OR Lung OR Vascular OR Immune System"
    st.session_state.alvo_val = "Bladder OR Vesical OR Urothelium OR Detrusor OR Cystitis OR Overactive Bladder"
    st.toast(t["toast_setup"], icon="🧬")

def carregar_termos_indicados_orgao(orgao, email, t):
    if not orgao or not email:
        st.error("E-mail e Alvo necessários!")
        return
    with st.spinner("Minerando PubMed com precisão farmacológica..."):
        termos = bk.buscar_alvos_emergentes_pubmed(orgao, email)
        if termos:
            atuais = set([x.strip().upper() for x in st.session_state.alvos_val.split(",") if x.strip()])
            filtrados = [term for term in termos if term.upper() not in atuais]
            if filtrados:
                novo_texto = (st.session_state.alvos_val.strip(", ") + ", " + ", ".join(filtrados)).strip(", ")
                st.session_state.alvos_val = novo_texto
                st.toast(t["toast_restaurar"], icon="✨")
            else:
                st.toast("Nenhum termo novo para adicionar.", icon="ℹ️")

def minerar_blue_oceans_ui(orgao, email, t):
    """Lógica de UI para mineração, chama APIs diretas para controle de progresso"""
    if not orgao or not email:
        st.toast(t["toast_aviso_minerar"], icon="⚠️"); return
    
    encontrados = []
    Entrez.email = email
    my_bar = st.progress(0, text=t["prog_minerar"])
    amostra = c.CANDIDATOS_MINERACAO 
    total = len(amostra)
    
    for i, termo in enumerate(amostra):
        try:
            time.sleep(0.3) 
            query = f"({termo}) AND ({orgao}) AND 2010:2025[DP]"
            handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
            record = Entrez.read(handle)
            count = int(record["Count"])
            if 0 <= count < 150: encontrados.append(f"{termo}")
            my_bar.progress((i + 1) / total, text=t["prog_testando"].format(termo=termo, count=count))
        except: continue
    
    my_bar.empty()
    if encontrados:
        st.session_state.alvos_val = ", ".join(encontrados)
        st.toast(t["toast_sucesso_minerar"].format(qtd=len(encontrados)), icon="💡")

def processar_upload(t):
    uploaded_file = st.session_state.get('uploader_key')
    if uploaded_file is not None:
        try:
            content = uploaded_file.getvalue().decode("utf-8")
            st.session_state.alvos_val = " ".join(content.replace("\n", ",").split())
            st.toast(t["toast_upload"], icon="📂")
        except: st.error("Erro upload")

def limpar_campo(key): st.session_state[key] = ""

# ==========================================
# INTERFACE PRINCIPAL (UI)
# ==========================================
lang_opt = st.sidebar.radio("Language / Idioma:", ["🇧🇷 Português", "🇺🇸 English"])
lang = "pt" if "Português" in lang_opt else "en"
t = c.TEXTOS[lang]
modo = st.sidebar.radio("📱 Mode:", ["Desktop", "Mobile (Pocket)"], index=0)

st.sidebar.markdown("---")
with st.sidebar.expander(t["citar_titulo"]):
    st.code(t["citar_texto"], language="text")
    st.link_button(t["link_doi"], "https://doi.org/10.5281/zenodo.17958507")
st.sidebar.markdown("---")

# --- MODO DESKTOP ---
if modo == "Desktop":
    st.title(t["titulo_desk"])
    st.markdown(t["subtitulo"])
    if 'dados_desk' not in st.session_state: exibir_radar_cientifico(lang)
    
    st.sidebar.header(t["credenciais"])
    email_user = st.sidebar.text_input(t["email_label"], placeholder="pesquisador@unifesp.br", key="email_desk")
    anos = st.sidebar.slider(t["periodo"], 1900, datetime.now().year, (2010, datetime.now().year), key="anos_desk")
    
    st.sidebar.markdown("---")
    st.sidebar.header(t["config"])
    st.sidebar.markdown(t["label_fonte"])
    c1, c2 = st.sidebar.columns([6, 1], vertical_alignment="bottom")
    with c1: t_fonte = st.text_input("Fonte", key="fonte_val", value=st.session_state.fonte_val, label_visibility="collapsed")
    with c2: st.button("🗑️", key="del_f", on_click=limpar_campo, args=("fonte_val",))
    
    st.sidebar.markdown(t["label_alvo"])
    c3, c4 = st.sidebar.columns([6, 1], vertical_alignment="bottom")
    with c3: t_alvo = st.text_input("Alvo", key="alvo_val", value=st.session_state.alvo_val, label_visibility="collapsed")
    with c4: st.button("🗑️", key="del_a", on_click=limpar_campo, args=("alvo_val",))
    st.sidebar.button(t["btn_setup"], on_click=carregar_setup_lemos, args=(t,))
    
    st.sidebar.markdown("---")
    st.sidebar.header(t["sec_alvos"])
    with st.sidebar.expander(t["expander_upload"]):
        st.file_uploader("Upload", type=["csv", "txt"], key="uploader_key", on_change=processar_upload, args=(t,))
    st.sidebar.markdown(t["label_lista"])
    c5, c6 = st.sidebar.columns([6, 1], vertical_alignment="bottom")
    with c5: 
        alvos_in = st.text_area(t["label_lista"], value=st.session_state.alvos_val, height=150, label_visibility="collapsed")
        st.session_state.alvos_val = alvos_in 
    with c6: st.button("🗑️", key="del_l", on_click=limpar_campo, args=("alvos_val",))

    if st.sidebar.button(t["btn_trend"], key="trend_desk"):
        carregar_termos_indicados_orgao(t_alvo, email_user, t)

    b1, b2 = st.sidebar.columns(2)
    b1.button(t["btn_restaurar"], on_click=carregar_termos_indicados_orgao, args=(t_alvo, email_user, t))
    b2.button(t["btn_minerar"], on_click=minerar_blue_oceans_ui, args=(t_alvo, email_user, t))
    
    st.sidebar.markdown("---")
    if st.sidebar.button(t["btn_avanco"], type="primary"):
        if not email_user: st.error(t["erro_email"])
        elif not st.session_state.alvos_val: st.warning(t["aviso_lista"])
        else:
            lst = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
            res = []; pg = st.empty(); bar = st.progress(0)
            for i, item in enumerate(lst):
                pg.text(t["prog_investigando"].format(atual=i+1, total=len(lst), alvo=item))
                nf = bk.consultar_pubmed_count(item, t_fonte, email_user, anos[0], anos[1])
                na = bk.consultar_pubmed_count(item, t_alvo, email_user, anos[0], anos[1])
                pot = nf/na if na > 0 else nf
                stat = "💎 DIAMANTE" if pot > 10 else "🥇 Ouro" if pot > 2 else "🔴 Saturado"
                res.append({"Alvo": item, "Status": stat, "Ratio": pot, "Fonte": nf, "Alvo_Interest": na})
                bar.progress((i+1)/len(lst))
            st.session_state['dados_desk'] = pd.DataFrame(res).sort_values(by="Ratio", ascending=False)
            st.rerun()

    if 'dados_desk' in st.session_state:
        df = st.session_state['dados_desk']
        st.success(t["analise_pronta"].format(top=df.iloc[0]['Alvo']))
        
        st.subheader("📊 Visualização de Potencial")
        filt_ops = df['Status'].unique().tolist()
        escolha = st.multiselect(t["filtro"], filt_ops, default=filt_ops)
        df_f = df[df['Status'].isin(escolha)].head(20)
        
        st.plotly_chart(px.bar(df_f, x="Alvo", y="Ratio", color="Status", 
                               color_discrete_map={"💎 DIAMANTE": "#00CC96", "🥇 Ouro": "#636EFA", "🔴 Saturado": "#EF553B"}), use_container_width=True)
        st.dataframe(df, use_container_width=True, hide_index=True)
        st.download_button(t["baixar"], df.to_csv(index=False).encode('utf-8'), "prospeccao_lemos.csv", "text/csv")
        
        st.divider()
        st.header(t["raio_x"])
        sel = st.selectbox("Selecione para ler artigos focados no alvo:", sorted(df['Alvo'].unique().tolist()))
        if st.button(t["btn_ler"]):
            with st.spinner(t["lendo"]):
                arts = bk.buscar_resumos_detalhados(sel, t_alvo, email_user, anos[0], anos[1], lang)
                if not arts: st.warning(t["sem_artigos"])
                for a in arts:
                    with st.expander(f"📄 {a['Title']}"): st.success(a['Resumo_IA'])
        if st.button(t["btn_scholar"]):
             st.markdown(f"👉 [Google Scholar](https://scholar.google.com.br/scholar?q={sel}+{t_alvo})")

# --- MODO MOBILE ---
elif modo == "Mobile (Pocket)":
    st.title(t["titulo_mob"])
    if 'dados_mob' not in st.session_state: exibir_radar_cientifico(lang)
    email_mob = st.text_input(t["email_label"], key="email_mob")
    with st.expander("⚙️ Config"):
        anos_mob = st.slider(t["periodo"], 1900, datetime.now().year, (2010, datetime.now().year))
        st.markdown(t["label_fonte"]); c1,c2=st.columns([6,1], vertical_alignment="bottom")
        with c1: t_fonte_m=st.text_input("F",key="fm", value=st.session_state.fonte_val, label_visibility="collapsed")
        with c2: st.button("🗑️",key="xf",on_click=limpar_campo, args=("fonte_val",))
        st.markdown(t["label_alvo"]); c3,c4=st.columns([6,1], vertical_alignment="bottom")
        with c3: t_alvo_m=st.text_input("A",key="am", value=st.session_state.alvo_val, label_visibility="collapsed")
        with c4: st.button("🗑️",key="xa",on_click=limpar_campo, args=("alvo_val",))
        st.button(t["btn_setup"], on_click=carregar_setup_lemos, args=(t,))
        st.markdown(t["label_lista"]); c5,c6=st.columns([6,1], vertical_alignment="bottom")
        with c5: 
            alvos_m = st.text_area(t["label_lista"], value=st.session_state.alvos_val, height=100, label_visibility="collapsed")
            st.session_state.alvos_val = alvos_m
        with c6: st.button("🗑️",key="xl",on_click=limpar_campo, args=("alvos_val",))
        if st.button(t["btn_trend"], key="trend_mob"):
            carregar_termos_indicados_orgao(t_alvo_m, email_mob, t)
        b1,b2=st.columns(2)
        b1.button(t["btn_restaurar"],on_click=carregar_termos_indicados_orgao,args=(t_alvo_m,email_mob,t))
        b2.button(t["btn_minerar"],on_click=minerar_blue_oceans_ui,args=(t_alvo_m,email_mob,t))

    if st.button(t["btn_avanco"], type="primary"):
        if not email_mob: st.error(t["erro_email"])
        else:
            l = [x.strip() for x in st.session_state.alvos_val.split(",") if x.strip()]
            r=[]; p=st.progress(0)
            for i, x in enumerate(l):
                time.sleep(0.1)
                nf = bk.consultar_pubmed_count(x, t_fonte_m, email_mob, anos_mob[0], anos_mob[1])
                na = bk.consultar_pubmed_count(x, t_alvo_m, email_mob, anos_mob[0], anos_mob[1])
                ratio = nf/na if na > 0 else nf
                r.append({"Alvo":x, "P":ratio})
                p.progress((i+1)/len(l))
            st.session_state['dados_mob'] = pd.DataFrame(r).sort_values(by="P", ascending=False)
            st.rerun()
    if 'dados_mob' in st.session_state:
        d=st.session_state['dados_mob']
        st.metric("🏆 Top 1", d.iloc[0]['Alvo'], f"{d.iloc[0]['P']:.1f}")
        st.dataframe(d, use_container_width=True, hide_index=True)
        st.download_button("📥 CSV", d.to_csv(index=False).encode('utf-8'), "lemos_lambda_mob.csv", "text/csv")