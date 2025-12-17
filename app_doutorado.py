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
# 1. CONFIGURAÇÃO GLOBAL
# ==========================================
st.set_page_config(page_title="Lemos Doutorado", page_icon="🎓", layout="wide")

# CSS para padronizar imagens (Cards)
st.markdown("""
    <style>
    div[data-testid="stImage"] img {
        height: 150px !important;
        object-fit: cover !important;
        border-radius: 8px !important;
    }
    </style>
""", unsafe_allow_html=True)

# Inicialização do Session State
if 'alvos_val' not in st.session_state: st.session_state.alvos_val = ""
if 'fonte_val' not in st.session_state: st.session_state.fonte_val = ""
if 'alvo_val' not in st.session_state: st.session_state.alvo_val = ""
if 'news_index' not in st.session_state: st.session_state.news_index = 0

# ==========================================
# 2. FUNÇÃO DE NOTÍCIAS
# ==========================================
@st.cache_data(ttl=3600)
def buscar_todas_noticias():
    feeds = [
        {"url": "https://www.sciencedaily.com/rss/health_medicine/pharmacology.xml", "lang": "🇺🇸"},
        {"url": "https://www.fiercebiotech.com/rss/biotech", "lang": "🇺🇸"},
        {"url": "https://www.nature.com/nbt.rss", "lang": "🇬🇧"},
        {"url": "https://agencia.fapesp.br/rss/", "lang": "🇧🇷"},
    ]
    
    noticias = []
    backups = [
        "https://images.unsplash.com/photo-1532094349884-543bc11b234d?w=400&h=250&fit=crop", 
        "https://images.unsplash.com/photo-1576086213369-97a306d36557?w=400&h=250&fit=crop", 
        "https://images.unsplash.com/photo-1581093458791-9f3c3900df4b?w=400&h=250&fit=crop", 
        "https://images.unsplash.com/photo-1530026405186-ed1f139313f8?w=400&h=250&fit=crop"
    ]
    
    translator = GoogleTranslator(source='auto', target='pt')

    for fonte in feeds:
        try:
            feed = feedparser.parse(fonte["url"])
            for entry in feed.entries[:3]:
                img_url = None
                if 'media_content' in entry: img_url = entry.media_content[0]['url']
                elif 'links' in entry:
                    for link in entry.links:
                        if link['type'].startswith('image'): img_url = link['href']; break
                elif 'summary' in entry:
                    match = re.search(r'src="(http.*?)"', entry.summary)
                    if match: img_url = match.group(1)
                
                if not img_url: img_url = random.choice(backups)

                titulo_orig = entry.title
                try:
                    if fonte["lang"] != "🇧🇷": titulo_pt = translator.translate(titulo_orig)
                    else: titulo_pt = titulo_orig
                except: titulo_pt = titulo_orig

                noticias.append({
                    "titulo_pt": titulo_pt,
                    "titulo_orig": titulo_orig,
                    "link": entry.link,
                    "fonte": feed.feed.title.split("-")[0].strip()[:20],
                    "img": img_url,
                    "bandeira": fonte["lang"]
                })
        except: continue
            
    random.shuffle(noticias)
    return noticias

@st.fragment(run_every=30) 
def exibir_radar_cientifico():
    news_list = buscar_todas_noticias()
    if not news_list:
        st.caption("Carregando feed...")
        return

    total_news = len(news_list)
    qtd = 3
    start = st.session_state.news_index % total_news
    end = start + qtd
    
    batch = []
    for i in range(start, end):
        batch.append(news_list[i % total_news])
    
    st.session_state.news_index += qtd
    
    with st.container(border=True):
        st.caption(f"📡 **Radar Científico**")
        cols = st.columns(3)
        for i, n in enumerate(batch):
            with cols[i]:
                st.image(n['img'], use_container_width=True)
                st.markdown(f"**{n['titulo_pt'][:75]}...**")
                if n['titulo_pt'] != n['titulo_orig']:
                    st.caption(f"🇬🇧 *{n['titulo_orig'][:75]}...*")
                st.caption(f"{n['bandeira']} {n['fonte']}")
                st.link_button("Ler", n['link'], use_container_width=True)

# ==========================================
# 3. BANCO DE DADOS (MEGA EXPANDIDO)
# ==========================================
SUGESTOES_ALVOS_RAW = """
-- GASOTRANSMISSORES & SINALIZAÇÃO GASOSA (O Novo Hype) --
Hydrogen Sulfide (H2S), CBS (Cystathionine beta-synthase), CSE (Cystathionine gamma-lyase), GYY4137, AP39, Nitric Oxide, sGC stimulators, Riociguat, Vericiguat, Carbon Monoxide (CO), Heme Oxygenase-1 (HO-1)

-- PURINÉRGICOS & ENDOCANABINÓIDES --
P2X1 receptor, P2X3, P2X7, P2Y6, P2Y12, Adenosine A2A, FAAH (Fatty Acid Amide Hydrolase), MAGL (Monoacylglycerol lipase), Anandamide, 2-AG, GPR55

-- CANAIS DE POTÁSSIO ESPECÍFICOS --
KATP channel, Kir6.1, Kir6.2, Glibenclamide, Cromakalim, SK channels, SK3, KCa2.3, Kv7.2, Kv7.3, Kv7.4, Retigabine, BKCa (Slo1)

-- GENÉTICA REGULATÓRIA (lncRNAs & microRNAs) --
MALAT1, HOTAIR, MEG3, H19, GAS5, miR-29b, miR-132, miR-199a, miR-21, miR-145, Antagomirs, siRNA therapy

-- COMUNICAÇÃO CELULAR (Exossomos & Vesículas) --
Exosomes, CD63, CD9, CD81, TSG101, Alix, Extracellular Vesicles, Microvesicles, MVBs, Gap Junctions, Connexin 43

-- IMUNOLOGIA AVANÇADA (Checkpoints em Inflamação) --
PD-1, PD-L1, CTLA-4, LAG-3, TIM-3, Siglec-8, Mast Cell Tryptase, Eosinophil Cationic Protein, IL-33, ST2 receptor

-- SENSORS "EXÓTICOS" (Olfato & Sabor na Bexiga) --
Olfactory Receptors, OR51E2, OR1D2, Taste Receptors, TAS2R (Bitter), TAS1R3 (Sweet), TRPM5, VN1R1

-- CRONOBIOLOGIA (Relógio da Bexiga) --
Clock genes, BMAL1, CLOCK, PER1, PER2, CRY1, Rev-erb alpha, Melatonin Receptor MT1, MT2

-- MECANO-BIOLOGIA & FIBROSE --
YAP, TAZ, Hippo pathway, Piezo1, Piezo2, Integrin beta-1, FAK, CTGF, LOX, Caveolin-1, Pirfenidone, Nintedanib

-- EPIGENÉTICA --
HDAC inhibitors, HDAC1, Valproic acid, Vorinostat, DNMT1, TET2, EZH2, Bromodomain

-- METABOLISMO MITOCONDRIAL --
Mitochondrial dynamics, Drp1, Mfn2, PGC-1alpha, Sirtuin-1, Sirtuin-3, NAMPT, NAD+

-- NOVAS VIAS DE MORTE --
Ferroptosis, GPX4, SLC7A11, Pyroptosis, Gasdermin D, Necroptosis, RIPK1, RIPK3, MLKL

-- TOXICOLOGIA AMBIENTAL --
Microplastics, Nanoplastics, Bisphenol S, Phthalates, Glyphosate, Acrolein, Cadmium

-- CANAIS IÔNICOS RAROS --
TMEM16A, HCN1, HCN4, Kv7.1, TREK-1, TRAAK, TRPML1
"""

def limpar_lista_alvos(texto_bruto):
    linhas = texto_bruto.split('\n')
    alvos_limpos = []
    for linha in linhas:
        if linha.strip() and not linha.strip().startswith("--"):
            itens = linha.split(',')
            for item in itens:
                item_limpo = item.split('(')[0].strip()
                if item_limpo: alvos_limpos.append(item_limpo)
    return ", ".join(alvos_limpos)

LISTA_ALVOS_PRONTA = limpar_lista_alvos(SUGESTOES_ALVOS_RAW)

PRESETS_ORGAOS = {
    "(Sugestão Lemos)": {
        "fonte": "Kidney OR Renal OR Vascular OR Intestine OR Lung OR Brain OR Liver OR Adipose Tissue OR Immune System",
        "alvo": "Bladder OR Vesical OR Urothelium OR Detrusor OR Cystitis OR Overactive Bladder"
    }
}

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
            else: st.error("Arquivo vazio ou inválido.")
        except Exception as e: st.error(f"Erro: {e}")

# ==========================================
# 4. FUNÇÕES PUBMED
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
# 5. INTERFACE (UI)
# ==========================================
modo = st.sidebar.radio("📱 Modo:", ["Desktop (Completo)", "Mobile (Pocket)"], index=0)
st.sidebar.markdown("---")

if modo == "Desktop (Completo)":
    st.title("🎓 Lemos Doutorado: Deep Science")
    st.markdown("**Ferramenta de Prospecção de Alto Impacto**")
    
    # Feed condicional
    if 'dados_desk' not in st.session_state:
        exibir_radar_cientifico()
    
    st.sidebar.header("1. Credenciais")
    email_user = st.sidebar.text_input("Seu E-mail:", placeholder="pesquisador@unifesp.br", key="email_desk")
    anos = st.sidebar.slider("📅 Período:", 1990, 2025, (2010, 2025), key="anos_desk")
    min_year, max_year = anos
    
    st.sidebar.markdown("---")
    st.sidebar.header("2. Configuração (Órgãos)")
    
    # ALINHAMENTO NATIVO (vertical_alignment="bottom") - Correção da Lixeira
    col_fonte, col_limp_f = st.sidebar.columns([0.85, 0.15], vertical_alignment="bottom")
    with col_fonte: termo_fonte = st.text_input("Fonte (Comparação):", key="fonte_val", placeholder="Sistemas Consolidados...")
    with col_limp_f: st.button("🗑️", key="btn_cls_f_dk", on_click=limpar_campo_fonte, help="Limpar")

    col_alvo, col_limp_a = st.sidebar.columns([0.85, 0.15], vertical_alignment="bottom")
    with col_alvo: termo_alvo = st.text_input("Alvo (Seu Foco):", key="alvo_val", placeholder="Bexiga/Urotélio...")
    with col_limp_a: st.button("🗑️", key="btn_cls_a_dk", on_click=limpar_campo_alvo, help="Limpar")
    
    st.sidebar.caption("👇 Setup Automático:")
    st.sidebar.button("🎓 Doutorado Guilherme Lemos", on_click=carregar_setup_lemos)
    
    st.sidebar.markdown("---")
    st.sidebar.header("3. Alvos")
    
    with st.sidebar.expander("📂 Importar Biblioteca (.csv/.txt)"):
        st.file_uploader("Upload", type=["csv", "txt"], key="uploader_key", on_change=processar_upload)
    
    col_lista, col_limp_l = st.sidebar.columns([0.85, 0.15], vertical_alignment="bottom")
    with col_lista: alvos_input = st.text_area("Lista de Pesquisa:", key="alvos_val", height=150, placeholder="Carregue a lista...")
    with col_limp_l: st.button("🗑️", key="btn_cls_l_dk", on_click=limpar_campo_alvos, help="Limpar")

    st.sidebar.button("📥 Restaurar Lista Padrão", on_click=carregar_alvos_apenas)
    st.sidebar.markdown("---")

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
                    elif ratio > 10 and n_fonte > 100: status = "💎 DIAMANTE"
                    elif ratio > 5 and n_fonte > 50: status = "🥇 Ouro"
                    elif ratio > 2: status = "🥈 Prata"
                    else: status = "🥚 Embrionário"
                    
                    resultados.append({
                        "Alvo": alvo, "Status": status, "Potencial (x)": round(ratio, 1),
                        "Fonte Total": n_fonte, "Bexiga Total": n_bexiga
                    })
                bar.progress((i+1)/len(alvos_lista))
            
            progresso_texto.empty()
            st.session_state['dados_desk'] = pd.DataFrame(resultados).sort_values(by="Potencial (x)", ascending=False)
            st.rerun()

    if 'dados_desk' in st.session_state:
        df = st.session_state['dados_desk']
        top = df.iloc[0]
        st.success(f"✅ Análise Pronta. Maior Potencial: **{top['Alvo']}**.")
        
        col1, col2 = st.columns([2, 1])
        with col1:
            fig = px.bar(df.head(20), x="Alvo", y="Potencial (x)", color="Status", 
                         title="Top 20 Alvos Inovadores", 
                         color_discrete_map={"💎 DIAMANTE": "#00CC96", "🥇 Ouro": "#636EFA", "🔴 Saturado": "#EF553B"})
            st.plotly_chart(fig, use_container_width=True)
        with col2:
            st.dataframe(df[["Alvo", "Status", "Potencial (x)", "Fonte Total", "Bexiga Total"]].style.hide(axis="index"), use_container_width=True, height=500)
            
            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button("📥 Baixar Planilha", csv, f'lemos_innov_{datetime.now().strftime("%Y%m%d")}.csv', 'text/csv', use_container_width=True)
            
        st.divider()
        st.header("🔎 Raio-X")
        sel = st.selectbox("Investigar Alvo:", sorted(df['Alvo'].unique().tolist()))
        col_btn1, col_btn2 = st.columns([1,4])
        if col_btn1.button("Ler Artigos"):
            with st.spinner("Buscando..."):
                arts = buscar_resumos_detalhados(sel, termo_alvo, email_user, min_year, max_year)
                if not arts: st.info(f"Zero artigos encontrados.")
                else:
                    for a in arts:
                        with st.expander(f"📄 {a['Title']}"):
                            st.write(f"**Revista:** {a['Source']}")
                            st.success(a['Resumo_IA'])
                            st.markdown(f"[Link PubMed](https://pubmed.ncbi.nlm.nih.gov/{a['PMID']})")
        if col_btn2.button("🎓 Google Scholar"):
             st.markdown(f"👉 [Abrir Scholar](https://scholar.google.com.br/scholar?q={sel}+AND+bladder)", unsafe_allow_html=True)

elif modo == "Mobile (Pocket)":
    st.title("📱 Lemos Pocket")
    
    if 'dados_mob' not in st.session_state:
        exibir_radar_cientifico() 
            
    email_mob = st.text_input("📧 E-mail:", placeholder="pesquisador@unifesp.br", key="email_mob")
    with st.expander("⚙️ Configurar"):
        anos_mob = st.slider("📅 Anos:", 1990, 2025, (2010, 2025))
        
        c1, c2 = st.columns([0.85, 0.15], vertical_alignment="bottom")
        with c1: t_fonte_mob = st.text_input("Fonte:", key="fonte_val", placeholder="Fonte...")
        with c2: st.button("🗑️", key="cls_f_mob", on_click=limpar_campo_fonte)
        
        c3, c4 = st.columns([0.85, 0.15], vertical_alignment="bottom")
        with c3: t_alvo_mob = st.text_input("Alvo:", key="alvo_val", placeholder="Alvo...")
        with c4: st.button("🗑️", key="cls_a_mob", on_click=limpar_campo_alvo)
        
        st.button("🎓 Doutorado Guilherme Lemos", key="mob_lemos", on_click=carregar_setup_lemos)
        st.markdown("---")
        
        st.file_uploader("📂 Upload", type=["csv", "txt"], key="uploader_key_mob", on_change=processar_upload)
        
        c5, c6 = st.columns([0.85, 0.15], vertical_alignment="bottom")
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
            st.rerun()
            
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
