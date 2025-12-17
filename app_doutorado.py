"""
Lemos Lambda: Deep Science Prospector
Copyright (c) 2025 Guilherme Lemos
Licensed under the MIT License.
"""

import streamlit as st
import pandas as pd
from Bio import Entrez
import time
import plotly.express as px
import re
from deep_translator import GoogleTranslator
from datetime import datetime
import io
import feedparser
import random

# ==========================================
# 0. FUNÇÕES DE INTELIGÊNCIA (DESCOBERTA)
# ==========================================
def buscar_alvos_emergentes(email):
    Entrez.email = email
    # Busca por termos de inovação no último ano (2024-2025)
    query = '("orphan receptor" OR "neglected target" OR "novel GPCR" OR "emerging target") AND ("2024"[Date - Publication] : "2025"[Date - Publication])'
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=15)
        record = Entrez.read(handle)
        if not record["IdList"]: return []
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="abstract", retmode="text")
        texto = handle.read()
        # Captura padrões como GPR123, P2X4, TRPV1
        encontrados = re.findall(r'\b[A-Z]{2,6}[0-9]{1,4}\b', texto)
        blacklist = ["DNA", "RNA", "USA", "NCBI", "NIH", "ATP", "AMP", "GDP", "COVID", "SARS", "FAPESP"]
        return sorted(list(set([t for t in encontrados if t not in blacklist and len(t) > 2])))
    except:
        return []

# ==========================================
# 1. CONFIGURAÇÃO GLOBAL
# ==========================================
st.set_page_config(page_title="Lemos Lambda", page_icon="λ", layout="wide")

st.markdown("""
    <style>
    div[data-testid="stImage"] img { height: 150px !important; object-fit: cover !important; border-radius: 8px !important; }
    .stButton button { width: 100%; }
    div[data-testid="stVerticalBlock"] > div { gap: 0.5rem; }
    </style>
""", unsafe_allow_html=True)

if 'alvos_val' not in st.session_state: st.session_state.alvos_val = ""
if 'fonte_val' not in st.session_state: st.session_state.fonte_val = ""
if 'alvo_val' not in st.session_state: st.session_state.alvo_val = ""
if 'news_index' not in st.session_state: st.session_state.news_index = 0

# ==========================================
# 2. DICIONÁRIO DE TRADUÇÃO (I18N)
# ==========================================
TEXTOS = {
    "pt": {
        "titulo_desk": "λ Lemos Lambda: Deep Science",
        "subtitulo": "**Ferramenta de Prospecção de Alto Impacto**",
        "titulo_mob": "📱 Lemos Pocket",
        "credenciais": "1. Credenciais",
        "email_label": "Seu E-mail:",
        "periodo": "📅 Período:",
        "config": "2. Configuração (Órgãos)",
        "label_fonte": "**Fonte (Órgão, tecido, célula similar):**",
        "holder_fonte": "Ex: Kidney...",
        "label_alvo": "**Alvo (Órgão de interesse):**",
        "holder_alvo": "Ex: Bladder...",
        "btn_setup": "🎓 Doutorado Guilherme Lemos",
        "toast_setup": "Setup 'Deep Science' Carregado!",
        "sec_alvos": "3. Alvos",
        "expander_upload": "📂 Importar Biblioteca (.csv/.txt)",
        "toast_upload": "Biblioteca importada!",
        "label_lista": "**Lista de Pesquisa:**",
        "holder_lista": "Carregue a lista...",
        "btn_restaurar": "📥 Restaurar Padrão",
        "toast_restaurar": "Lista Inovadora Restaurada!",
        "btn_minerar": "⛏️ Minerar 'Blue Oceans'",
        "toast_aviso_minerar": "⚠️ Preencha o 'Alvo' e 'E-mail' para minerar!",
        "prog_minerar": "⛏️ Procurando termos chave, após isso clique em 'Rumo ao Avanço'...",
        "prog_testando": "⛏️ Analisando: {termo} ({count} artigos)",
        "toast_sucesso_minerar": "✅ {qtd} termos encontrados! Agora clique em 'Rumo ao Avanço' 🚀",
        "toast_fail_minerar": "Nenhum alvo raro encontrado.",
        "btn_avanco": "🚀 Rumo ao Avanço",
        "erro_email": "E-mail obrigatório!",
        "aviso_lista": "Lista de Alvos vazia!",
        "prog_investigando": "⏳ Investigando {atual}/{total}: {alvo}",
        "analise_pronta": "✅ Análise Pronta. Destaque: **{top}**.",
        "col_artigos": "Artigos",
        "col_global": "Global",
        "col_ratio": "Ratio",
        "filtro": "🔍 Filtro:",
        "grafico_qtd": "📊 Qtd. no Gráfico:",
        "raio_x": "🔎 Raio-X",
        "btn_ler": "Ler Artigos",
        "btn_scholar": "🎓 Google Scholar",
        "sem_artigos": "Zero artigos encontrados.",
        "lendo": "Buscando e Traduzindo...",
        "baixar": "📥 Baixar Planilha",
        "citar_titulo": "📄 Como Citar",
        "citar_texto": "Lemos, G. (2025). Lemos Lambda: Deep Science Prospector [Software]. Versão 1.0.0. DOI: 10.5281/zenodo.17958507",
        "link_doi": "🔗 Ver no Zenodo (DOI)"
    },
    "en": {
        "titulo_desk": "λ Lemos Lambda: Deep Science",
        "subtitulo": "**High Impact Prospecting Tool**",
        "titulo_mob": "📱 Lemos Pocket",
        "credenciais": "1. Credentials",
        "email_label": "Your E-mail:",
        "periodo": "📅 Timeframe:",
        "config": "2. Configuration (Organs)",
        "label_fonte": "**Source (Organ, tissue, similar cell):**",
        "holder_fonte": "Ex: Kidney...",
        "label_alvo": "**Target (Organ of interest):**",
        "holder_alvo": "Ex: Bladder...",
        "btn_setup": "🎓 Guilherme Lemos PhD Setup",
        "toast_setup": "'Deep Science' Setup Loaded!",
        "sec_alvos": "3. Targets",
        "expander_upload": "📂 Import Library (.csv/.txt)",
        "toast_upload": "Library imported!",
        "label_lista": "**Search List:**",
        "holder_lista": "Load the list...",
        "btn_restaurar": "📥 Restore Default",
        "toast_restaurar": "Innovative List Restored!",
        "btn_minerar": "⛏️ Mine 'Blue Oceans'",
        "toast_aviso_minerar": "⚠️ Fill in 'Target' and 'E-mail' to mine!",
        "prog_minerar": "⛏️ Searching for key terms, then click 'Launch'...",
        "prog_testando": "⛏️ Analyzing: {termo} ({count} papers)",
        "toast_sucesso_minerar": "✅ {qtd} terms found! Now click 'Launch' 🚀",
        "toast_fail_minerar": "No rare targets found.",
        "btn_avanco": "🚀 Launch Analysis",
        "erro_email": "E-mail required!",
        "aviso_lista": "Target list is empty!",
        "prog_investigando": "⏳ Investigating {atual}/{total}: {alvo}",
        "analise_pronta": "✅ Analysis Ready. Highlight: **{top}**.",
        "col_artigos": "Papers",
        "col_global": "Global",
        "col_ratio": "Ratio",
        "filtro": "🔍 Filter:",
        "grafico_qtd": "📊 Chart Qty:",
        "raio_x": "🔎 X-Ray",
        "btn_ler": "Read Papers",
        "btn_scholar": "🎓 Google Scholar",
        "sem_artigos": "Zero papers found.",
        "lendo": "Searching and Translating...",
        "baixar": "📥 Download CSV",
        "citar_titulo": "📄 How to Cite",
        "citar_texto": "Lemos, G. (2025). Lemos Lambda: Deep Science Prospector [Software]. Version 1.0.0. DOI: 10.5281/zenodo.17958507",
        "link_doi": "🔗 View on Zenodo (DOI)"
    }
}

# ==========================================
# 3. FUNÇÕES DE SUPORTE
# ==========================================
@st.cache_data(ttl=3600)
def buscar_todas_noticias(lang_code):
    feeds = [
        {"url": "https://www.sciencedaily.com/rss/health_medicine/pharmacology.xml", "lang": "🇺🇸"},
        {"url": "https://www.fiercebiotech.com/rss/biotech", "lang": "🇺🇸"},
        {"url": "https://www.nature.com/nbt.rss", "lang": "🇬🇧"},
        {"url": "https://agencia.fapesp.br/rss/", "lang": "🇧🇷"},
    ]
    noticias = []
    backups = ["https://images.unsplash.com/photo-1532094349884-543bc11b234d?w=400&h=250&fit=crop"]
    translator = GoogleTranslator(source='auto', target=lang_code)

    for fonte in feeds:
        try:
            feed = feedparser.parse(fonte["url"])
            for entry in feed.entries[:3]:
                img_url = random.choice(backups)
                if 'media_content' in entry: img_url = entry.media_content[0]['url']
                elif 'links' in entry:
                    for link in entry.links:
                        if link['type'].startswith('image'): img_url = link['href']; break
                
                titulo = entry.title
                if lang_code == 'pt' and fonte["lang"] != "🇧🇷":
                    try: titulo = translator.translate(titulo)
                    except: pass
                
                noticias.append({
                    "titulo": titulo, "link": entry.link,
                    "fonte": feed.feed.title.split("-")[0].strip()[:20], 
                    "img": img_url, "bandeira": fonte["lang"]
                })
        except: continue
    random.shuffle(noticias)
    return noticias

@st.fragment(run_every=60) 
def exibir_radar_cientifico(lang_code):
    news_list = buscar_todas_noticias(lang_code)
    if not news_list: st.caption("Loading feed..."); return
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

CANDIDATOS_MINERACAO = ["GPR37", "GPR17", "GPR55", "GPR84", "GPR35", "TAAR1", "P2Y14", "FFAR2", "FFAR3", "SUCNR1", "OXGR1", "HCAR1", "LGR4", "LGR5", "MRGPRX2", "MRGPRD", "ASIC1a", "ASIC2", "TRPML1", "TRPML2", "TRPML3", "TPC1", "TPC2", "TMEM16A", "TMEM16B", "Piezo1", "Piezo2", "TREK-1", "TREK-2", "TRAAK", "TASK-1", "TASK-3", "TWIK-1", "THIK-1", "KCa3.1", "Kv7.5", "ClC-2", "Bestrophin-1", "Pannexin-1", "S1P1", "S1P2", "S1P3", "LPA1", "LPA2", "CysLT1", "CysLT2", "FPR2", "ChemR23", "BLT1", "GYY4137", "AP39", "Drp1", "Mfn2", "Sirtuin-1", "NAMPT", "Ferroptosis", "Pyroptosis", "Necroptosis", "Gasdermin D", "TAS2R", "TAS2R10", "TAS2R14", "Olfactory Receptors", "OR51E2", "SLC7A11", "NLRP3", "Resolvin D1", "Maresin 1", "Lipoxin A4", "Itaconate", "P2X4"]

SUGESTOES_ALVOS_RAW = "Piezo1, Piezo2, TREK-1, TRAAK, TASK-1, GPR35, GPR55, GPR84, GPR183, TAS2R, TAS2R10, TAS2R14, Olfactory Receptors, OR51E2, Ferroptosis, GPX4, SLC7A11, Pyroptosis, Gasdermin D, NLRP3, H2S, GYY4137, CSE enzyme, CBS enzyme, Resolvin D1, Maresin 1, Lipoxin A4, Itaconate, Pannexin-1, P2X4, TMEM16A"
LISTA_ALVOS_PRONTA = ", ".join([x.strip() for x in SUGESTOES_ALVOS_RAW.split(',') if x.strip()])

PRESETS_ORGAOS = {"(Sugestão Lemos)": {"fonte": "Brain OR Kidney OR Liver OR Intestine OR Lung OR Vascular OR Immune System", "alvo": "Bladder OR Vesical OR Urothelium OR Detrusor OR Cystitis OR Overactive Bladder"}}

def carregar_setup_lemos(t):
    st.session_state.alvos_val = LISTA_ALVOS_PRONTA
    st.session_state.fonte_val = PRESETS_ORGAOS["(Sugestão Lemos)"]["fonte"]
    st.session_state.alvo_val = PRESETS_ORGAOS["(Sugestão Lemos)"]["alvo"]
    st.toast(t["toast_setup"], icon="🧬")

def carregar_alvos_apenas(t): 
    st.session_state.alvos_val = LISTA_ALVOS_PRONTA
    st.toast(t["toast_restaurar"], icon="✨")

def limpar_campo_fonte(): st.session_state.fonte_val = ""
def limpar_campo_alvo(): st.session_state.alvo_val = ""
def limpar_campo_alvos(): st.session_state.alvos_val = ""

def minerar_blue_oceans(orgao, email, t):
    if not orgao or not email: st.toast(t["toast_aviso_minerar"], icon="⚠️"); return
    encontrados = []
    Entrez.email = email
    my_bar = st.progress(0, text=t["prog_minerar"])
    amostra = CANDIDATOS_MINERACAO 
    total = len(amostra)
    for i, termo in enumerate(amostra):
        try:
            query = f"({termo}) AND ({orgao}) AND 2010:2025[DP]"
            record = Entrez.read(Entrez.esearch(db="pubmed", term=query, retmax=0))
            count = int(record["Count"])
            if 0 <= count < 150: encontrados.append(f"{termo}")
            my_bar.progress((i + 1) / total, text=t["prog_testando"].format(termo=termo, count=count))
            time.sleep(0.05)
        except: continue
    my_bar.empty()
    if encontrados:
        st.session_state.alvos_val = ", ".join(encontrados)
        st.toast(t["toast_sucesso_minerar"].format(qtd=len(encontrados)), icon="💡")
    else: st.toast(t["toast_fail_minerar"], icon="🤷")

def processar_upload(t):
    uploaded_file = st.session_state.get('uploader_key')
    if uploaded_file:
        try:
            content = uploaded_file.getvalue().decode("utf-8")
            st.session_state.alvos_val = " ".join(content.replace("\n", ",").split())
            st.toast(t["toast_upload"], icon="📂")
        except: st.error("Erro upload")

def consultar_pubmed_count(termo_farmaco, termo_orgao, email, y_start, y_end):
    if not email: return -1
    Entrez.email = email
    query = f"({termo_farmaco}) AND ({termo_orgao}) AND {y_start}:{y_end}[DP]" if termo_orgao else f"({termo_farmaco}) AND {y_start}:{y_end}[DP]"
    try: return int(Entrez.read(Entrez.esearch(db="pubmed", term=query, retmax=0))["Count"])
    except: return -1

def extrair_conclusao(abstract_text, lang_target):
    if not abstract_text: return "Resumo não disponível." if lang_target == 'pt' else "Abstract not available."
    match = re.search(r'(Conclusion|Conclusions|In conclusion|Summary|Results suggest that)(.*)', abstract_text, re.IGNORECASE | re.DOTALL)
    texto_final = match.group(2).strip()[:400] if match else abstract_text[-400:]
    return ("🇧🇷 " if lang_target=='pt' else "🇺🇸 ") + GoogleTranslator(source='auto', target=lang_target).translate(texto_final) + "..."

def buscar_resumos_detalhados(termo_farmaco, termo_orgao, email, y_start, y_end, lang_target, limit=5):
    if not email: return []
    query = f"({termo_farmaco}) AND ({termo_orgao}) AND {y_start}:{y_end}[DP]" if termo_orgao else f"({termo_farmaco}) AND {y_start}:{y_end}[DP]"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=limit, sort="relevance")
        ids = Entrez.read(handle)["IdList"]
        if not ids: return []
        records = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text").read().split("\n\n")
        artigos = []
        for art_text in records:
            art_data = {"PMID": "N/A", "Title": "S/T", "Source": "N/A", "Abstract": ""}
            for line in art_text.split("\n"):
                if len(line)<4: continue
                tag, content = line[:4].strip(), line[6:]
                if tag=="PMID": art_data["PMID"]=content
                elif tag=="TI": art_data["Title"]=content
                elif tag=="TA": art_data["Source"]=content
                elif tag=="AB": art_data["Abstract"]=content
            if art_data["PMID"]!="N/A":
                art_data["Resumo_IA"] = extrair_conclusao(art_data["Abstract"], lang_target)
                artigos.append(art_data)
        return artigos
    except: return []

# ==========================================
# 4. INTERFACE (UI)
# ==========================================
lang_opt = st.sidebar.radio("Language / Idioma:", ["🇧🇷 Português", "🇺🇸 English"])
lang = "pt" if "Português" in lang_opt else "en"
t = TEXTOS[lang]

modo = st.sidebar.radio("📱 Mode:", ["Desktop", "Mobile (Pocket)"], index=0)

st.sidebar.markdown("---")
with st.sidebar.expander(t["citar_titulo"]):
    st.code(t["citar_texto"], language="text")
    st.link_button(t["link_doi"], "https://doi.org/10.5281/zenodo.17958507")
st.sidebar.markdown("---")

if modo == "Desktop":
    st.title(t["titulo_desk"])
    st.markdown(t["subtitulo"])
    if 'dados_desk' not in st.session_state: exibir_radar_cientifico(lang)
    
    st.sidebar.header(t["credenciais"])
    email_user = st.sidebar.text_input(t["email_label"], placeholder="pesquisador@unifesp.br", key="email_desk")

    # Módulo de Inteligência (Desktop)
    st.sidebar.markdown("---")
    st.sidebar.subheader("💡 Inteligência Lemos Lambda")
    st.sidebar.caption("Minere alvos promissores/órfãos da literatura recente (2024-2025).")
    if st.sidebar.button("🔍 Descobrir Tendências", key="btn_trend_desk"):
        if email_user and "@" in email_user:
            with st.sidebar.status("Varrendo o PubMed..."):
                novos_termos = buscar_alvos_emergentes(email_user)
                if novos_termos:
                    st.sidebar.success(f"Encontrei {len(novos_termos)} termos!")
                    st.sidebar.text_area("Copie e cole na lista de alvos:", ", ".join(novos_termos), height=100)
                else: st.sidebar.warning("Nada novo encontrado.")
        else: st.sidebar.error("E-mail válido obrigatório.")
    st.sidebar.markdown("---")

    anos = st.sidebar.slider(t["periodo"], 1990, 2025, (2010, 2025), key="anos_desk")
    st.sidebar.header(t["config"])
    
    st.sidebar.markdown(t["label_fonte"])
    c1, c2 = st.sidebar.columns([6, 1], vertical_alignment="bottom")
    with c1: t_fonte = st.text_input("Fonte", key="fonte_val", placeholder=t["holder_fonte"], label_visibility="collapsed")
    with c2: st.button("🗑️", key="del_f", on_click=limpar_campo_fonte)

    st.sidebar.markdown(t["label_alvo"])
    c3, c4 = st.sidebar.columns([6, 1], vertical_alignment="bottom")
    with c3: t_alvo = st.text_input("Alvo", key="alvo_val", placeholder=t["holder_alvo"], label_visibility="collapsed")
    with c4: st.button("🗑️", key="del_a", on_click=limpar_campo_alvo)
    
    st.sidebar.button(t["btn_setup"], on_click=carregar_setup_lemos, args=(t,))
    st.sidebar.markdown("---")
    st.sidebar.header(t["sec_alvos"])
    
    with st.sidebar.expander(t["expander_upload"]):
        st.file_uploader("Upload", type=["csv", "txt"], key="uploader_key", on_change=processar_upload, args=(t,))
    
    c5, c6 = st.sidebar.columns([6, 1], vertical_alignment="bottom")
    with c5: alvos_in = st.text_area(t["label_lista"], key="alvos_val", height=150, placeholder=t["holder_lista"])
    with c6: st.button("🗑️", key="del_l", on_click=limpar_campo_alvos)

    b1, b2 = st.sidebar.columns(2)
    b1.button(t["btn_restaurar"], on_click=carregar_alvos_apenas, args=(t,))
    b2.button(t["btn_minerar"], on_click=minerar_blue_oceans, args=(t_alvo, email_user, t))

    if st.sidebar.button(t["btn_avanco"], type="primary"):
        if not email_user: st.error(t["erro_email"])
        elif not alvos_in: st.warning(t["aviso_lista"])
        else:
            lst = [x.strip() for x in alvos_in.split(",") if x.strip()]
            res = []; pg = st.empty(); bar = st.progress(0)
            for i, item in enumerate(lst):
                pg.text(t["prog_investigando"].format(atual=i+1, total=len(lst), alvo=item))
                nf = consultar_pubmed_count(item, t_fonte, email_user, anos[0], anos[1]) if t_fonte else 0
                na = consultar_pubmed_count(item, t_alvo, email_user, anos[0], anos[1]) if t_alvo else 0
                ng = consultar_pubmed_count(item, "", email_user, anos[0], anos[1]) if not t_fonte and not t_alvo else 0
                pot = nf/na if na > 0 and t_fonte else (na if t_alvo else ng)
                stat = "💎 DIAMANTE" if t_fonte and t_alvo and pot > 10 and nf > 50 else "🥇 Ouro"
                res.append({"Alvo": item, "Status": stat, "Potencial": pot, "Qtd_Fonte": nf, "Qtd_Alvo": na if t_alvo else ng})
                bar.progress((i+1)/len(lst))
            st.session_state['dados_desk'] = pd.DataFrame(res).sort_values(by="Potencial", ascending=False); st.rerun()

    if 'dados_desk' in st.session_state:
        df = st.session_state['dados_desk']
        st.success(t["analise_pronta"].format(top=df.iloc[0]['Alvo']))
        c_g1, c_g2 = st.columns(2)
        with c_g1: qtd_graf = st.slider(t["grafico_qtd"], 10, 100, 20)
        with c_g2: 
            ops = df['Status'].unique().tolist()
            filt = st.multiselect(t["filtro"], ops, default=ops)
        df_f = df[df['Status'].isin(filt)].head(qtd_graf)
        col1, col2 = st.columns([2, 1])
        with col1: st.plotly_chart(px.bar(df_f, x="Alvo", y="Potencial", color="Status"), use_container_width=True)
        with col2: 
            st.dataframe(df_f, use_container_width=True)
            st.download_button(t["baixar"], df.to_csv(index=False).encode('utf-8'), "lemos_analise.csv")

elif modo == "Mobile (Pocket)":
    st.title(t["titulo_mob"])
    if 'dados_mob' not in st.session_state: exibir_radar_cientifico(lang)
    email_mob = st.text_input(t["email_label"], key="email_mob")

    # Módulo de Inteligência (Mobile)
    if st.button("🔍 Descobrir Tendências", key="btn_trend_mob"):
        if email_mob and "@" in email_mob:
            with st.spinner("Buscando..."):
                novos = buscar_alvos_emergentes(email_mob)
                if novos:
                    st.success(f"{len(novos)} encontrados!")
                    st.text_area("Copie:", ", ".join(novos))
                else: st.warning("Nada novo.")
        else: st.error("E-mail primeiro.")

    with st.expander("⚙️ Config"):
        anos_mob = st.slider(t["periodo"], 1990, 2025, (2010, 2025))
        t_fonte_m = st.text_input("Fonte", key="fm", placeholder=t["holder_fonte"])
        t_alvo_m = st.text_input("Alvo", key="am", placeholder=t["holder_alvo"])
        st.button(t["btn_setup"], on_click=carregar_setup_lemos, args=(t,))
        alvos_m = st.text_area("Lista", key="alm")
        b1, b2 = st.columns(2)
        b1.button(t["btn_restaurar"], on_click=carregar_alvos_apenas, args=(t,))
        b2.button(t["btn_minerar"], on_click=minerar_blue_oceans, args=(t_alvo_m, email_mob, t))

    if st.button(t["btn_avanco"], type="primary"):
        if not email_mob: st.error(t["erro_email"])
        else:
            l = [x.strip() for x in alvos_m.split(",") if x.strip()]
            r = []; p = st.progress(0)
            for i, x in enumerate(l):
                nf = consultar_pubmed_count(x, t_fonte_m, email_mob, anos_mob[0], anos_mob[1]) if t_fonte_m else 0
                na = consultar_pubmed_count(x, t_alvo_m, email_mob, anos_mob[0], anos_mob[1]) if t_alvo_m else 0
                pot = nf/na if na > 0 and t_fonte_m else na
                r.append({"Alvo": x, "P": pot})
                p.progress((i+1)/len(l))
            st.session_state['dados_mob'] = pd.DataFrame(r).sort_values(by="P", ascending=False); st.rerun()

    if 'dados_mob' in st.session_state:
        st.dataframe(st.session_state['dados_mob'], use_container_width=True)
