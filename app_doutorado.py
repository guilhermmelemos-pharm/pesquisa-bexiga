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
Version: 2.3 (CSS-Only Image Sizing)
"""

import streamlit as st
import pandas as pd
import plotly.express as px
from datetime import datetime
import time
import constantes as c
import backend as bk

st.set_page_config(page_title="Lemos Lambda", page_icon="λ", layout="wide")

# --- CSS ---
st.markdown("""
<style>
.stButton button { border-radius:12px; height:50px; font-weight:bold; }
@keyframes pulse-blue {
  0% { box-shadow:0 0 0 0 rgba(0,204,150,.7); }
  70% { box-shadow:0 0 0 10px rgba(0,204,150,0); }
  100% { box-shadow:0 0 0 0 rgba(0,204,150,0); }
}
.blue-ocean-btn button {
  animation:pulse-blue 2s infinite;
  background:linear-gradient(90deg,#00C9FF,#92FE9D)!important;
  color:#004d40!important;
  border:none!important;
}
div[data-testid="stImage"] img{
  width:100%!important; height:150px!important; object-fit:cover!important;
}
</style>
""", unsafe_allow_html=True)

# --- SESSION STATE INIT ---
_defaults = {
    'pagina':'home','alvos_val':'','resultado_df':None,'news_index':0,
    'input_alvo':'','input_fonte':'','input_email':'',
    'artigos_detalhe':None,'email_guardado':'','alvo_guardado':''
}
for k,v in _defaults.items(): st.session_state.setdefault(k,v)

lang_opt = st.sidebar.radio("🌐 Language", ["🇧🇷 PT","🇺🇸 EN"], horizontal=True)
lang = "pt" if "PT" in lang_opt else "en"
t = c.TEXTOS[lang]

# ================= FUNÇÕES =================
def resetar_pesquisa():
    st.session_state.pagina='home'
    st.session_state.resultado_df=None
    st.session_state.artigos_detalhe=None

def limpar_campo(k): st.session_state[k]=''
def limpar_lista_total(): st.session_state.alvos_val=''

def adicionar_termos_seguro(lista):
    atuais=[x.strip() for x in st.session_state.alvos_val.split(',') if x.strip()]
    atuais_u=[x.upper() for x in atuais]
    adicionados=0
    for t_ in lista:
        if t_.upper() not in atuais_u:
            atuais.append(t_); atuais_u.append(t_.upper()); adicionados+=1
    st.session_state.alvos_val=', '.join(atuais)
    return adicionados

# ================= BLUE OCEAN =================
def ir_para_analise(email, contexto, alvo, ano_i, ano_f):
    st.session_state.email_guardado=email
    st.session_state.alvo_guardado=alvo
    lista=[x.strip() for x in st.session_state.alvos_val.split(',') if x.strip()]
    resultados=[]
    prog=st.progress(0)
    for i,item in enumerate(lista):
        time.sleep(0.03)
        n_global=bk.consultar_pubmed_count(item, contexto or None, email, ano_i, ano_f)
        n_especifico=bk.consultar_pubmed_count(item, alvo, email, ano_i, ano_f)
        status,score=bk.classificar_oportunidade(n_especifico,n_global)
        ratio=n_global/max(n_especifico,1)
        resultados.append({
            'term':item,'status':status,'ratio':round(ratio,1),
            'alvo_count':n_especifico,'fonte_count':n_global,'_sort':score
        })
        prog.progress((i+1)/len(lista))
    df=pd.DataFrame(resultados).sort_values(by=['_sort','ratio'],ascending=False)
    st.session_state.resultado_df=df
    st.session_state.pagina='resultados'
    st.rerun()

# ================= UI =================
if st.session_state.pagina=='home':
    st.title(t['titulo_desk']); st.caption(t['subtitulo'])
    st.text_input(t['label_email'], key='input_email')
    st.text_input(t['label_alvo'], key='input_alvo')
    st.text_area('Palavras-chave', key='alvos_val', height=120)
    anos=st.slider('Anos',2000,datetime.now().year,(2015,datetime.now().year))
    if st.button(t['analise_btn'], type='primary', use_container_width=True):
        ir_para_analise(st.session_state.input_email, st.session_state.input_fonte,
                        st.session_state.input_alvo, anos[0], anos[1])

elif st.session_state.pagina=='resultados':
    st.button(t['btn_nova_pesquisa'], on_click=resetar_pesquisa)
    df=st.session_state.resultado_df
    top=df.iloc[0]
    c1,c2,c3=st.columns(3)
    c1.metric(t['metrica_potencial'], top['term'], delta=top['status'])
    c2.metric(t['metrica_score'], top['ratio'])
    c3.metric(t['metrica_artigos'], top['alvo_count'])
    df_show=df.rename(columns={
        'term':t['col_mol'],'status':t['col_status'],'ratio':t['col_ratio'],
        'alvo_count':t['col_art_alvo'],'fonte_count':t['col_global']
    }).drop(columns=['_sort'])
    fig=px.bar(df_show.head(25), x=t['col_mol'], y=t['col_ratio'], color=t['col_status'])
    st.plotly_chart(fig, use_container_width=True)
    st.dataframe(df_show, use_container_width=True, hide_index=True)

st.markdown('---'); st.caption(f"© 2025 Guilherme Lemos")
