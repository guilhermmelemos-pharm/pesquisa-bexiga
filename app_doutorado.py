import streamlit as st
import pandas as pd
from Bio import Entrez
import time
import plotly.express as px
import re
from deep_translator import GoogleTranslator
from datetime import datetime

# ==========================================
# 1. CONFIGURAÇÃO GLOBAL
# ==========================================
st.set_page_config(page_title="Lemos Buscador UX", page_icon="🧬", layout="wide")

# Inicialização do Session State
if 'alvos_val' not in st.session_state: st.session_state.alvos_val = ""
if 'fonte_val' not in st.session_state: st.session_state.fonte_val = ""
if 'alvo_val' not in st.session_state: st.session_state.alvo_val = ""

# ==========================================
# 2. BANCO DE DADOS (SUA INTELIGÊNCIA)
# ==========================================
SUGESTOES_ALVOS = """
-- ALVOS MAIS PROMISSORES (PRIORIDADE) --
Autophagy, LC3B, Beclin-1, p62, mTOR, AMPK, VEGF, VEGFR2, TGF-beta1, CTGF, Galectin-3, P2X3, P2X7, TRPV1, TRPV4, TRPM8, Beta-3 Adrenergic, Muscarinic M3, Cannabinoid CB2

-- FÁRMACOS & TOXINAS --
Mirabegron, Solifenacin, Oxybutynin, Botulinum toxin A (BoNT/A), Resiniferatoxin (RTX), Tadalafil, Sildenafil, Rapamycin, Metformin, Silodosin, Tamsulosin

-- CANAIS IÔNICOS (MUSCULO LISO & NERVO) --
BK channel (KCa1.1), SK3 channel, Kv7.4 (KCNQ4), Kv7.5, KATP channel (Kir6.2), L-type Calcium Channel (Cav1.2), T-type Calcium Channel, Piezo1, Piezo2, ASIC1, ASIC3, TRPA1, TRPC6

-- RECEPTORES GPCRS --
Alpha-1A Adrenergic, Alpha-1D Adrenergic, Beta-2 Adrenergic, Muscarinic M2, Dopamine D2, Serotonin 5-HT, Adenosine A1, Adenosine A2A, P2Y receptors, Angiotensin II receptor (AT1R), Mas receptor

-- INFLAMAÇÃO, DOR & NEUROPEPTÍDEOS --
NLRP3, IL-1beta, IL-6, IL-17, IL-33, TNF-alpha, COX-2, PGE2, NGF, BDNF, CGRP, Substance P, VIP, PACAP, Bradykinin B1, Bradykinin B2

-- METABOLISMO, HORMONIOS & ENZIMAS --
Estrogen Receptor Alpha (ESR1), Estrogen Receptor Beta (ESR2), Androgen Receptor, ROCK (Rho-kinase), RhoA, PDE4, PDE5, nNOS, eNOS, iNOS, SGLT2, ACE2, Nrf2, HO-1, Sirtuin-1
"""
LISTA_ALVOS_LIMPA = " ".join(SUGESTOES_ALVOS.replace("\n", " ").split())

# PRESETS DE ÓRGÃOS
PRESETS_ORGAOS = {
    "(Sugestão Lemos)": {
        "fonte": "Kidney OR Renal OR Blood Vessels OR Vascular OR Intestine OR Gut OR Lung OR Airway OR Uterus OR Prostate OR Heart OR Cardiac OR Smooth Muscle",
        "alvo": "Bladder OR Vesical OR Urothelium OR Detrusor OR Cystitis OR Overactive Bladder OR Painful Bladder"
    },
    "Rim/Vaso -> Bexiga": {
        "fonte": "Kidney OR Renal OR Blood Vessels OR Vascular OR Hypertension",
        "alvo": "Bladder OR Vesical OR Urothelium OR Detrusor OR Cystitis"
    },
    "Próstata/Uretra -> Bexiga": {
        "fonte": "Prostate OR Prostatic OR Urethra OR Urethral OR BPH OR LUTS",
        "alvo": "Bladder OR Detrusor OR Overactive Bladder"
    },
    "Pulmão/Asma -> Bexiga": {
        "fonte": "Lung OR Pulmonary OR Airway OR Bronchial OR Asthma OR COPD",
        "alvo": "Bladder OR Detrusor OR Smooth Muscle"
    },
    "Útero -> Bexiga": {
        "fonte": "Uterus OR Uterine OR Myometrium OR Endometrium",
        "alvo": "Bladder OR Detrusor OR Smooth Muscle"
    },
    "Coração -> Bexiga": {
        "fonte": "Heart OR Cardiac OR Myocardium OR Cardiomyocyte",
        "alvo": "Bladder OR Detrusor OR Smooth Muscle"
    },
    "Intestino -> Bexiga": {
        "fonte": "Gut OR Intestine OR Colon OR Bowel OR Enteric",
        "alvo": "Bladder OR Cystitis"
    }
}

def carregar_setup_lemos():
    st.session_state.alvos_val = LISTA_ALVOS_LIMPA
    st.session_state.fonte_val = PRESETS_ORGAOS["(Sugestão Lemos)"]["fonte"]
    st.session_state.alvo_val = PRESETS_ORGAOS["(Sugestão Lemos)"]["alvo"]

def carregar_alvos_apenas(): st.session_state.alvos_val = LISTA_ALVOS_LIMPA
def carregar_orgaos(preset_name):
    dados = PRESETS_ORGAOS[preset_name]
    st.session_state.fonte_val = dados["fonte"]
    st.session_state.alvo_val = dados["alvo"]

# ==========================================
# 3. FUNÇÕES TÉCNICAS
# ==========================================
def consultar_pubmed_count(termo_farmaco, termo_orgao, email, y_start, y_end):
    if not email: return -1
    Entrez.email = email
    termo_farmaco = termo_farmaco.replace(",", "").strip()
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
    Entrez.email = email
    termo_farmaco = termo_farmaco.replace(",", "").strip()
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
            current_tag = ""
            for line in lines:
                if len(line)<4: continue
                tag, content = line[:4].strip(), line[6:]
                if tag=="PMID": art_data["PMID"]=content
                elif tag=="TI": art_data["Title"]=content; current_tag="TI"
                elif tag=="TA": art_data["Source"]=content
                elif tag=="AB": art_data["Abstract"]=content; current_tag="AB"
                elif tag=="" and current_tag=="AB": art_data["Abstract"]+=" "+line.strip()
                elif tag=="" and current_tag=="TI": art_data["Title"]+=" "+line.strip()
            if art_data["PMID"]!="N/A":
                art_data["Resumo_IA"] = extrair_conclusao(art_data["Abstract"])
                artigos.append(art_data)
        return artigos
    except: return []

# ==========================================
# 4. INTERFACE
# ==========================================
modo = st.sidebar.radio("📱 Modo:", ["Desktop (Completo)", "Mobile (Pocket)"], index=0)
st.sidebar.markdown("---")

if modo == "Desktop (Completo)":
    st.title("🔬 Lemos Ultimate Edition")
    st.markdown("**Ferramenta Bibliométrica Personalizada**")

    st.sidebar.header("1. Identificação")
    email_user = st.sidebar.text_input("Seu E-mail:", placeholder="pesquisador@unifesp.br", key="email_desk")
    anos = st.sidebar.slider("📅 Período:", 1990, 2025, (2010, 2025), key="anos_desk")
    min_year, max_year = anos
    
    st.sidebar.markdown("---")
    
    # --- SETUP MASTER ---
    st.sidebar.header("2. Configuração de Órgãos")
    
    # 1º: INPUTS VISÍVEIS (O Usuário vê onde o texto entra)
    termo_fonte = st.sidebar.text_input("Fonte:", key="fonte_val", placeholder="Orgãos de Comparação...")
    termo_alvo = st.sidebar.text_input("Alvo:", key="alvo_val", placeholder="Orgão Alvo (Bexiga)...")
    
    # 2º: BOTÕES EMBAIXO DOS INPUTS (Para preencher)
    st.sidebar.caption("👇 Ou preencha com um clique:")
    
    # Botão Master
    if st.sidebar.button("🧪 (Sugestão Lemos) - Comparar Tudo", type="primary"): 
        carregar_setup_lemos()
        st.toast("Setup Carregado!", icon="✅")
    
    # Botões Específicos
    c1, c2 = st.sidebar.columns(2)
    if c1.button("Rim ➡️ Bexiga"): carregar_orgaos("Rim/Vaso -> Bexiga")
    if c2.button("Próstata ➡️ Bexiga"): carregar_orgaos("Próstata/Uretra -> Bexiga")
    
    c3, c4 = st.sidebar.columns(2)
    if c3.button("Pulmão ➡️ Bexiga"): carregar_orgaos("Pulmão/Asma -> Bexiga")
    if c4.button("Intestino ➡️ Bexiga"): carregar_orgaos("Intestino -> Bexiga")
    
    st.sidebar.markdown("---")
    
    # --- CONFIGURAÇÃO DE ALVOS ---
    st.sidebar.header("3. Lista de Alvos")
    
    # 1º: TEXT AREA VISÍVEL
    alvos_input = st.sidebar.text_area("Lista de Pesquisa:", key="alvos_val", height=150, placeholder="Digite ou carregue a lista...")
    
    # 2º: BOTÃO EMBAIXO
    if st.sidebar.button("📥 Restaurar Lista Completa (+90 Alvos)"): 
        carregar_alvos_apenas()

    st.sidebar.markdown("---")

    if st.sidebar.button("🚀 INICIAR VARREDURA", type="secondary"):
        if not email_user or "@" not in email_user: st.error("E-mail obrigatório!")
        elif not termo_fonte or not termo_alvo: st.warning("Configure os órgãos!")
        elif not alvos_input: st.warning("Lista vazia!")
        else:
            alvos_lista = [x.strip() for x in alvos_input.split(",") if x.strip()]
            resultados = []
            
            progresso_texto = st.empty()
            bar = st.progress(0)
            
            for i, alvo in enumerate(alvos_lista):
                progresso_texto.text(f"⏳ Processando {i+1}/{len(alvos_lista)}: {alvo}")
                n_fonte = consultar_pubmed_count(alvo, termo_fonte, email_user, min_year, max_year)
                n_bexiga = consultar_pubmed_count(alvo, termo_alvo, email_user, min_year, max_year)
                if n_fonte != -1:
                    ratio = n_fonte / n_bexiga if n_bexiga > 0 else n_fonte
                    resultados.append({"Alvo": alvo, "Fonte Total": n_fonte, "Alvo Total": n_bexiga, "Potencial": round(ratio, 1)})
                bar.progress((i+1)/len(alvos_lista))
            
            progresso_texto.empty()
            st.session_state['dados_desk'] = pd.DataFrame(resultados).sort_values(by="Potencial", ascending=False)

    if 'dados_desk' in st.session_state:
        df = st.session_state['dados_desk']
        top = df.iloc[0]
        st.success(f"✅ Processado! Destaque: **{top['Alvo']}**.")
        
        col1, col2 = st.columns([2, 1])
        with col1:
            fig = px.bar(df.head(20), x="Alvo", y="Potencial", color="Potencial", title="Top 20 Oportunidades", color_continuous_scale="Bluered")
            st.plotly_chart(fig, use_container_width=True)
        with col2:
            st.dataframe(df[["Alvo", "Fonte Total", "Alvo Total", "Potencial"]].style.background_gradient(subset=['Potencial'], cmap="Greens").hide(axis="index"), use_container_width=True, height=500)
            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button("📥 Baixar Planilha Completa", csv, f'analise_lemos_{datetime.now().strftime("%Y%m%d")}.csv', 'text/csv', use_container_width=True)
            
        st.divider()
        st.header("🔎 Raio-X & Tradução")
        sel = st.selectbox("Investigar Alvo:", df['Alvo'].tolist())
        
        col_btn1, col_btn2 = st.columns([1,4])
        if col_btn1.button("Ler Artigos"):
            with st.spinner("Buscando..."):
                arts = buscar_resumos_detalhados(sel, termo_alvo, email_user, min_year, max_year)
                if not arts: st.info("Zero artigos encontrados.")
                else:
                    for a in arts:
                        with st.expander(f"📄 {a['Title']}"):
                            st.write(f"**Revista:** {a['Source']}")
                            st.success(a['Resumo_IA'])
                            st.markdown(f"[Link PubMed](https://pubmed.ncbi.nlm.nih.gov/{a['PMID']})")
        
        if col_btn2.button("🎓 Buscar no Google Scholar"):
             st.markdown(f"👉 [Clique aqui para abrir **{sel} + Bexiga** no Google Scholar](https://scholar.google.com.br/scholar?q={sel}+AND+bladder)", unsafe_allow_html=True)

elif modo == "Mobile (Pocket)":
    st.title("📱 Lemos Pocket")
    email_mob = st.text_input("📧 E-mail:", placeholder="pesquisador@unifesp.br", key="email_mob")
    
    with st.expander("⚙️ Configurar Busca"):
        anos_mob = st.slider("📅 Anos:", 1990, 2025, (2010, 2025))
        
        # UI Mobile corrigida também
        t_fonte_mob = st.text_input("Fonte:", key="fonte_val", placeholder="Fonte...")
        t_alvo_mob = st.text_input("Alvo:", key="alvo_val", placeholder="Alvo...")
        
        st.caption("Preenchimento Rápido:")
        if st.button("🧪 SETUP COMPLETO", key="mob_lemos", type="primary"): carregar_setup_lemos()
        
        st.markdown("---")
        alvos_mob = st.text_area("Alvos:", key="alvos_val", height=150)
        if st.button("📥 Restaurar Lista", key="mob_alvos"): carregar_alvos_apenas()
        
    if st.button("🚀 INICIAR", type="secondary", use_container_width=True):
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
                    res.append({"Alvo": al, "Potencial": round(rat, 1)})
                pg.progress((i+1)/len(lst))
            st.session_state['dados_mob'] = pd.DataFrame(res).sort_values(by="Potencial", ascending=False)
            
    if 'dados_mob' in st.session_state:
        d = st.session_state['dados_mob']
        t = d.iloc[0]
        st.divider()
        st.metric("🏆 Top 1", t['Alvo'], f"{t['Potencial']}x")
        csv_mob = d.to_csv(index=False).encode('utf-8')
        st.download_button("📥 Baixar CSV", csv_mob, "mobile.csv", "text/csv", use_container_width=True)
        with st.expander("Ver Lista"): st.dataframe(d, use_container_width=True, hide_index=True)
        st.divider()
        sl = st.selectbox("Ler:", d['Alvo'].tolist())
        if st.button("Ler", use_container_width=True):
            with st.spinner("Traduzindo..."):
                as_mob = buscar_resumos_detalhados(sl, t_alvo_mob, email_mob, anos_mob[0], anos_mob[1], limit=3)
                if not as_mob: st.info("Sem artigos!")
                else:
                    for am in as_mob:
                        st.success(f"**{am['Title']}**\n\n{am['Resumo_IA']}")
                        st.write("---")
