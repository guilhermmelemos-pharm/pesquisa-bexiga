import streamlit as st
import pandas as pd
from Bio import Entrez
import time
import plotly.express as px
import re
from deep_translator import GoogleTranslator
from datetime import datetime
import io
import feedparser # BIBLIOTECA NOVA PARA NOTÍCIAS

# ==========================================
# 1. CONFIGURAÇÃO GLOBAL
# ==========================================
st.set_page_config(page_title="Lemos Doutorado", page_icon="🎓", layout="wide")

# Inicialização do Session State
if 'alvos_val' not in st.session_state: st.session_state.alvos_val = ""
if 'fonte_val' not in st.session_state: st.session_state.fonte_val = ""
if 'alvo_val' not in st.session_state: st.session_state.alvo_val = ""

# ==========================================
# 2. FUNÇÃO DE NOTÍCIAS (LIVE SCIENCE)
# ==========================================
@st.cache_data(ttl=3600) # Atualiza a cada 1 hora para não pesar
def buscar_noticias_ciencia():
    feeds = [
        # ScienceDaily - Pharmacology
        "https://www.sciencedaily.com/rss/health_medicine/pharmacology.xml",
        # FierceBiotech - Biotech Industry
        "https://www.fiercebiotech.com/rss/biotech",
    ]
    noticias = []
    for url in feeds:
        try:
            feed = feedparser.parse(url)
            for entry in feed.entries[:2]: # Pega as 2 mais recentes de cada
                noticias.append({
                    "titulo": entry.title,
                    "link": entry.link,
                    "fonte": feed.feed.title.split(" - ")[0] if " - " in feed.feed.title else feed.feed.title
                })
        except:
            continue
    return noticias

# ==========================================
# 3. BANCO DE DADOS (FRONTEIRA EXTREMA)
# ==========================================
SUGESTOES_ALVOS_RAW = """
-- GENÉTICA REGULATÓRIA (lncRNAs & microRNAs) --
MALAT1, HOTAIR, MEG3, H19, GAS5, miR-29b, miR-132, miR-199a, miR-21, miR-145, Antagomirs

-- COMUNICAÇÃO CELULAR (Exossomos & Vesículas) --
Exosomes, CD63, CD9, CD81, TSG101, Alix, Extracellular Vesicles, Microvesicles, MVBs

-- IMUNOLOGIA AVANÇADA (Checkpoints) --
PD-1 (Programmed cell death protein 1), PD-L1, CTLA-4, LAG-3, TIM-3, Siglec-8, Mast Cell Tryptase

-- SENSORS "EXÓTICOS" (Olfato/Sabor na Bexiga) --
Olfactory Receptors, OR51E2, OR1D2, Taste Receptors, TAS2R, TAS1R3, TRPM5

-- CRONOBIOLOGIA (Relógio Biológico) --
Clock genes, BMAL1, CLOCK, PER1, PER2, CRY1, Rev-erb alpha, Melatonin Receptor MT1

-- MECANO-BIOLOGIA & FIBROSE --
YAP, TAZ, Hippo pathway, Piezo1, Piezo2, Integrin beta-1, FAK, CTGF, LOX

-- EPIGENÉTICA --
HDAC inhibitors, HDAC1, Valproic acid, Vorinostat, DNMT1, TET2, EZH2

-- METABOLISMO MITOCONDRIAL --
Mitochondrial dynamics, Drp1, Mfn2, PGC-1alpha, Sirtuin-1, NAMPT, NAD+

-- NOVAS VIAS DE MORTE --
Ferroptosis, GPX4, SLC7A11, Pyroptosis, Gasdermin D, Necroptosis, RIPK1, MLKL

-- TOXICOLOGIA AMBIENTAL --
Microplastics, Nanoplastics, Bisphenol S, Phthalates, Glyphosate, Acrolein

-- CANAIS IÔNICOS RAROS --
TMEM16A, HCN1, HCN4, Kv7.1, TREK-1, TRAAK
"""

def limpar_lista_alvos(texto_bruto):
    linhas = texto_bruto.split('\n')
    alvos_limpos = []
    for linha in linhas:
        if linha.strip() and not linha.strip().startswith("--"):
            itens = linha.split(',')
            for item in itens:
                item_limpo = item.split('(')[0].strip()
                if item_limpo:
                    alvos_limpos.append(item_limpo)
    return ", ".join(alvos_limpos)

LISTA_ALVOS_PRONTA = limpar_lista_alvos(SUGESTOES_ALVOS_RAW)

PRESETS_ORGAOS = {
    "(Sugestão Lemos)": {
        "fonte": "Kidney OR Renal OR Vascular OR Intestine OR Lung OR Brain OR Liver OR Adipose Tissue OR Immune System",
        "alvo": "Bladder OR Vesical OR Urothelium OR Detrusor OR Cystitis OR Overactive Bladder"
    }
}

# --- CALLBACKS & HELPERS ---
def carregar_setup_lemos():
    st.session_state.alvos_val = LISTA_ALVOS_PRONTA
    st.session_state.fonte_val = PRESETS_ORGAOS["(Sugestão Lemos)"]["fonte"]
    st.session_state.alvo_val = PRESETS_ORGAOS["(Sugestão Lemos)"]["alvo"]
    st.toast("Setup 'Deep Science' Carregado!", icon="🧬")

def carregar_alvos_apenas(): 
    st.session_state.alvos_val = LISTA_ALVOS_PRONTA
    st.toast("Lista Inovadora Restaurada!", icon="✨")

def limpar_campo_fonte(): st.session_state.fonte_val = ""
def limpar_campo_alvo(): st.session_state.alvo_val = ""
def limpar_campo_alvos(): st.session_state.alvos_val = ""

def processar_upload():
    uploaded_file = st.session_state.get('uploader_key')
    if uploaded_file is not None:
        try:
            string_final = ""
            if uploaded_file.name.endswith('.csv'):
                df = pd.read_csv(uploaded_file, header=None)
                lista_itens = df.stack().dropna().astype(str).tolist()
                string_final = ", ".join(lista_itens)
            elif uploaded_file.name.endswith('.txt'):
                stringio = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
                conteudo = stringio.read()
                string_final = " ".join(conteudo.replace("\n", ",").split())
            
            if string_final:
                st.session_state.alvos_val = string_final
                st.toast(f"Biblioteca '{uploaded_file.name}' importada!", icon="📂")
            else:
                st.error("Arquivo vazio ou inválido.")
        except Exception as e:
            st.error(f"Erro: {e}")

# ==========================================
# 4. FUNÇÕES PUBMED & TRADUÇÃO
# ==========================================
def consultar_pubmed_count(termo_farmaco, termo_orgao, email, y_start, y_end):
    if not email: return -1
    Entrez.email = email
    query = f"({termo_farmaco}) AND ({termo_orgao}) AND {y_start}:{y_end}[DP]"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle)
        return int(record["Count"])
    except: return -1

def traduzir_para_pt(texto):
    try: return GoogleTranslator(source='auto', target='pt').translate(texto)
    except: return texto

def extrair_conclusao(abstract_text):
    if not abstract_text: return "Resumo não disponível."
    match = re.search(r'(Conclusion|Conclusions|In conclusion|Summary|Results suggest that)(.*)', abstract_text, re.IGNORECASE | re.DOTALL)
    texto_final = match.group(2).strip()[:400] if match else abstract_text[-400:]
    return "🇧🇷 " + traduzir_para_pt(texto_final) + "..."

def buscar_resumos_detalhados(termo_farmaco, termo_orgao, email, y_start, y_end, limit=5):
    if not email: return []
    query = f"({termo_farmaco}) AND ({termo_orgao}) AND {y_start}:{y_end}[DP]"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=limit, sort="relevance")
        record = Entrez.read(handle)
        id_list = record["IdList"]
        if not id_list: return []
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        records = handle.read()
        artigos = []
        raw_articles = records.split("\n\n")
        for art_text in raw_articles:
            lines = art_text.split("\n")
            art_data = {"PMID": "N/A", "Title": "S/T", "Source": "N/A", "Abstract": ""}
            for line in lines:
                if len(line)<4: continue
                tag, content = line[:4].strip(), line[6:]
                if tag=="PMID": art_data["PMID"]=content
                elif tag=="TI": art_data["Title"]=content
                elif tag=="TA": art_data["Source"]=content
                elif tag=="AB": art_data["Abstract"]=content
            if art_data["PMID"]!="N/A":
                art_data["Resumo_IA"] = extrair_conclusao(art_data["Abstract"])
                artigos.append(art_data)
        return artigos
    except: return []

# ==========================================
# 5. INTERFACE (FRONT-END)
# ==========================================
modo = st.sidebar.radio("📱 Modo:", ["Desktop (Completo)", "Mobile (Pocket)"], index=0)
st.sidebar.markdown("---")

if modo == "Desktop (Completo)":
    # --- CABEÇALHO COM NEWS TICKER ---
    st.title("🎓 Lemos Doutorado: Deep Science")
    st.markdown("**Ferramenta de Prospecção de Alto Impacto**")
    
    # RADAR CIENTÍFICO (NEWS)
    with st.expander("📡 Radar Científico (Últimas Notícias - ScienceDaily & FierceBiotech)", expanded=False):
        news = buscar_noticias_ciencia()
        if news:
            for n in news:
                st.markdown(f"**[{n['fonte']}]** [{n['titulo']}]({n['link']})")
        else:
            st.caption("Carregando notícias ou feed indisponível...")
    
    # --- SIDEBAR DESKTOP ---
    st.sidebar.header("1. Credenciais")
    email_user = st.sidebar.text_input("Seu E-mail:", placeholder="pesquisador@unifesp.br", key="email_desk")
    anos = st.sidebar.slider("📅 Período:", 1990, 2025, (2010, 2025), key="anos_desk")
    min_year, max_year = anos
    
    st.sidebar.markdown("---")
    st.sidebar.header("2. Configuração (Órgãos)")
    
    # Layout colado (Lixeiras perto) - Ratio [5, 1]
    col_fonte, col_limp_f = st.sidebar.columns([5, 1])
    with col_fonte: termo_fonte = st.text_input("Fonte (Comparação):", key="fonte_val", placeholder="Sistemas Consolidados...")
    with col_limp_f: 
        st.write(""); st.write("")
        st.button("🗑️", key="btn_cls_f_dk", on_click=limpar_campo_fonte, help="Limpar")

    col_alvo, col_limp_a = st.sidebar.columns([5, 1])
    with col_alvo: termo_alvo = st.text_input("Alvo (Seu Foco):", key="alvo_val", placeholder="Bexiga/Urotélio...")
    with col_limp_a: 
        st.write(""); st.write("")
        st.button("🗑️", key="btn_cls_a_dk", on_click=limpar_campo_alvo, help="Limpar")
    
    st.sidebar.caption("👇 Setup Automático:")
    # BOTÃO CINZA (SECUNDÁRIO)
    st.sidebar.button("🎓 Doutorado Guilherme Lemos", on_click=carregar_setup_lemos)
    
    st.sidebar.markdown("---")
    st.sidebar.header("3. Alvos")
    
    with st.sidebar.expander("📂 Importar Biblioteca (.csv/.txt)"):
        st.file_uploader("Upload:", type=["csv", "txt"], key="uploader_key", on_change=processar_upload)
    
    col_lista, col_limp_l = st.sidebar.columns([5, 1])
    with col_lista: alvos_input = st.text_area("Lista de Pesquisa:", key="alvos_val", height=150, placeholder="Carregue a lista...")
    with col_limp_l: 
        st.write(""); st.write("")
        st.button("🗑️", key="btn_cls_l_dk", on_click=limpar_campo_alvos, help="Limpar")

    st.sidebar.button("📥 Restaurar Lista Padrão", on_click=carregar_alvos_apenas)
    st.sidebar.markdown("---")

    # BOTÃO VERMELHO (PRIMÁRIO) - AÇÃO PRINCIPAL
    if st.sidebar.button("🚀 Rumo ao Avanço", type="primary"):
        if not email_user or "@" not in email_user: st.error("E-mail obrigatório!")
        elif not termo_fonte or not termo_alvo: st.warning("Configure os órgãos!")
        elif not alvos_input: st.warning("Lista vazia!")
        else:
            alvos_lista = [x.strip() for x in alvos_input.split(",") if x.strip()]
            resultados = []
            progresso_texto = st.empty()
            bar = st.progress(0)
            
            for i, alvo in enumerate(alvos_lista):
                progresso_texto.text(f"⏳ Investigando {i+1}/{len(alvos_lista)}: {alvo}")
                n_fonte = consultar_pubmed_count(alvo, termo_fonte, email_user, min_year, max_year)
                n_bexiga = consultar_pubmed_count(alvo, termo_alvo, email_user, min_year, max_year)
                
                if n_fonte != -1:
                    ratio = n_fonte / n_bexiga if n_bexiga > 0 else n_fonte
                    status = "N/A"
                    if n_bexiga >= n_fonte and n_bexiga > 10: status = "🔴 Saturado"
                    elif ratio > 10 and n_fonte > 100: status = "💎 DIAMANTE (Inédito)"
                    elif ratio > 5 and n_fonte > 50: status = "🥇 Ouro (Promissor)"
                    elif ratio > 2: status = "🥈 Prata"
                    else: status = "🥚 Embrionário (Risco)"
                    
                    resultados.append({
                        "Alvo": alvo, "Status": status, "Potencial (x)": round(ratio, 1),
                        "Fonte Total": n_fonte, "Bexiga Total": n_bexiga
                    })
                bar.progress((i+1)/len(alvos_lista))
            
            progresso_texto.empty()
            st.session_state['dados_desk'] = pd.DataFrame(resultados).sort_values(by="Potencial (x)", ascending=False)

    # --- RESULTADOS ---
    if 'dados_desk' in st.session_state:
        df = st.session_state['dados_desk']
        top = df.iloc[0]
        st.success(f"✅ Análise Pronta. Maior Potencial: **{top['Alvo']}**.")
        
        col1, col2 = st.columns([2, 1])
        with col1:
            fig = px.bar(df.head(20), x="Alvo", y="Potencial (x)", color="Status", 
                         title="Top 20 Alvos Inovadores", 
                         color_discrete_map={"💎 DIAMANTE (Inédito)": "#00CC96", "🥇 Ouro (Promissor)": "#636EFA", "🔴 Saturado": "#EF553B"})
            st.plotly_chart(fig, use_container_width=True)
        with col2:
            st.dataframe(df[["Alvo", "Status", "Potencial (x)", "Fonte Total", "Bexiga Total"]].hide(axis="index"), use_container_width=True, height=500)
            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button("📥 Baixar Planilha", csv, f'lemos_innov_{datetime.now().strftime("%Y%m%d")}.csv', 'text/csv', use_container_width=True)
            
        st.divider()
        st.header("🔎 Raio-X")
        sel = st.selectbox("Investigar Alvo:", sorted(df['Alvo'].unique().tolist()))
        col_btn1, col_btn2 = st.columns([1,4])
        if col_btn1.button("Ler Artigos"):
            with st.spinner("Buscando..."):
                arts = buscar_resumos_detalhados(sel, termo_alvo, email_user, min_year, max_year)
                if not arts: st.info(f"Zero artigos sobre {sel} na bexiga! Você pode ser o primeiro.")
                else:
                    for a in arts:
                        with st.expander(f"📄 {a['Title']}"):
                            st.write(f"**Revista:** {a['Source']}")
                            st.success(a['Resumo_IA'])
                            st.markdown(f"[Link PubMed](https://pubmed.ncbi.nlm.nih.gov/{a['PMID']})")
        if col_btn2.button("🎓 Google Scholar"):
             st.markdown(f"👉 [Abrir Scholar: **{sel} + Bexiga**](https://scholar.google.com.br/scholar?q={sel}+AND+bladder)", unsafe_allow_html=True)

elif modo == "Mobile (Pocket)":
    st.title("📱 Lemos Pocket")
    # News Ticker Mobile
    with st.expander("📡 News", expanded=False):
        news = buscar_noticias_ciencia()
        if news:
            for n in news: st.markdown(f"- [{n['titulo']}]({n['link']})")
            
    email_mob = st.text_input("📧 E-mail:", placeholder="pesquisador@unifesp.br", key="email_mob")
    with st.expander("⚙️ Configurar"):
        anos_mob = st.slider("📅 Anos:", 1990, 2025, (2010, 2025))
        
        c1, c2 = st.columns([5,1])
        with c1: t_fonte_mob = st.text_input("Fonte:", key="fonte_val", placeholder="Fonte...")
        with c2: st.button("🗑️", key="cls_f_mob", on_click=limpar_campo_fonte)
        
        c3, c4 = st.columns([5,1])
        with c3: t_alvo_mob = st.text_input("Alvo:", key="alvo_val", placeholder="Alvo...")
        with c4: st.button("🗑️", key="cls_a_mob", on_click=limpar_campo_alvo)
        
        st.button("🎓 Doutorado Guilherme Lemos", key="mob_lemos", on_click=carregar_setup_lemos)
        st.markdown("---")
        
        st.file_uploader("📂 Upload", type=["csv", "txt"], key="uploader_key_mob", on_change=processar_upload)
        
        c5, c6 = st.columns([5,1])
        with c5: alvos_mob = st.text_area("Alvos:", key="alvos_val", height=150)
        with c6: st.button("🗑️", key="cls_l_mob", on_click=limpar_campo_alvos)
        
        st.button("📥 Restaurar", key="mob_alvos", on_click=carregar_alvos_apenas)
        
    if st.button("🚀 Rumo ao Avanço", type="primary", use_container_width=True):
        if not email_mob: st.error("E-mail necessário")
        else:
            lst = [x.strip() for x in alvos_mob.split(",") if x.strip()]
            res = []
            pg = st.progress(0)
            for i, al in enumerate(lst):
                nf = consultar_pubmed_count(al, t_fonte_mob, email_mob, anos_mob[0], anos_mob[1])
                nb = consultar_pubmed_count(al, t_alvo_mob, email_mob, anos_mob[0], anos_mob[1])
                if nf!=-1:
                    rat = nf/nb if nb>0 else nf
                    stat = "N/A"
                    if nb >= nf: stat = "🔴"
                    elif rat > 10: stat = "💎"
                    else: stat = "👍"
                    res.append({"Alvo": al, "Status": stat, "Potencial": round(rat, 1)})
                pg.progress((i+1)/len(lst))
            st.session_state['dados_mob'] = pd.DataFrame(res).sort_values(by="Potencial", ascending=False)
            
    if 'dados_mob' in st.session_state:
        d = st.session_state['dados_mob']
        t = d.iloc[0]
        st.divider()
        st.metric("🏆 Top 1", t['Alvo'], f"{t['Potencial']}x ({t['Status']})")
        csv_mob = d.to_csv(index=False).encode('utf-8')
        st.download_button("📥 Baixar CSV", csv_mob, "mobile.csv", "text/csv", use_container_width=True)
        with st.expander("Ver Lista"): st.dataframe(d, use_container_width=True, hide_index=True)
        st.divider()
        sl = st.selectbox("Ler:", sorted(d['Alvo'].unique().tolist()))
        if st.button("Ler", use_container_width=True):
            with st.spinner("Traduzindo..."):
                as_mob = buscar_resumos_detalhados(sl, t_alvo_mob, email_mob, anos_mob[0], anos_mob[1], limit=3)
                if not as_mob: st.info("Sem artigos!")
                else:
                    for am in as_mob:
                        st.success(f"**{am['Title']}**\n\n{am['Resumo_IA']}")
                        st.write("---")
