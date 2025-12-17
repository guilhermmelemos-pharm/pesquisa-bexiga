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
st.set_page_config(page_title="Lemos Doutorado", page_icon="🎓", layout="wide")

# Inicialização do Session State
if 'alvos_val' not in st.session_state: st.session_state.alvos_val = ""
if 'fonte_val' not in st.session_state: st.session_state.fonte_val = ""
if 'alvo_val' not in st.session_state: st.session_state.alvo_val = ""

# ==========================================
# 2. BANCO DE DADOS (LISTA INTERNA EXPANDIDA)
# ==========================================
# Aqui está a lista MESTRA. Se quiser adicionar algo novo, basta escrever aqui dentro.
SUGESTOES_ALVOS_RAW = """
-- ALVOS PRIORITÁRIOS (Mecanismos Celulares) --
Autophagy, LC3B, Beclin-1, p62, mTOR, AMPK, ULK1, VEGF, VEGFR2, TGF-beta1, CTGF, Galectin-3, P2X3, P2X7, TRPV1, TRPV4, TRPM8, Beta-3 Adrenergic, Muscarinic M3, Cannabinoid CB2

-- FÁRMACOS, TOXINAS & AGONISTAS --
Mirabegron, Solifenacin, Oxybutynin, Botulinum toxin A, Resiniferatoxin, Tadalafil, Sildenafil, Rapamycin, Metformin, Silodosin, Tamsulosin, Capsaicin, Beta-sitosterol

-- CANAIS IÔNICOS (O Nicho Eletrofisiológico) --
BK channel (KCa1.1), SK3 channel, Kv7.4, Kv7.5, KATP channel, Kv2.1, L-type Calcium Channel, T-type Calcium Channel, Piezo1, Piezo2, ASIC1, ASIC3, TRPA1, TRPC6, TMEM16A (Anoctamin-1), HCN1, HCN4

-- RECEPTORES (GPCRs & Nucleares) --
Alpha-1A Adrenergic, Alpha-1D Adrenergic, Beta-2 Adrenergic, Muscarinic M2, Dopamine D2, Serotonin 5-HT, Adenosine A1, Adenosine A2A, P2Y receptors, Angiotensin II receptor, Mas receptor, Vitamin D Receptor (VDR), GPER (Estrogen Receptor), EP1 receptor, EP2 receptor, EP3 receptor, EP4 receptor

-- INFLAMAÇÃO, BARREIRA & SINALIZAÇÃO --
NLRP3, IL-1beta, IL-6, IL-17, IL-33, TNF-alpha, COX-2, PGE2, NGF, BDNF, CGRP, Substance P, VIP, PACAP, Bradykinin B1, Bradykinin B2, ROCK, RhoA, PDE4, PDE5, nNOS, eNOS, iNOS, SGLT2, ACE2, Nrf2, HO-1, Sirtuin-1, Claudin-1, Occludin, Uroplakin
"""

# FUNÇÃO DE LIMPEZA (Transforma o texto acima em lista para busca)
def limpar_lista_alvos(texto_bruto):
    linhas = texto_bruto.split('\n')
    alvos_limpos = []
    for linha in linhas:
        if linha.strip() and not linha.strip().startswith("--"):
            # Separa por vírgula e remove parênteses explicativos se houver
            itens = linha.split(',')
            for item in itens:
                item_limpo = item.split('(')[0].strip() # Remove descrições entre parênteses para a busca ser exata
                if item_limpo:
                    alvos_limpos.append(item_limpo)
    return ", ".join(alvos_limpos)

LISTA_ALVOS_PRONTA = limpar_lista_alvos(SUGESTOES_ALVOS_RAW)

PRESETS_ORGAOS = {
    "(Sugestão Lemos)": {
        "fonte": "Kidney OR Renal OR Blood Vessels OR Vascular OR Intestine OR Gut OR Lung OR Airway OR Uterus OR Prostate OR Heart OR Cardiac OR Smooth Muscle",
        "alvo": "Bladder OR Vesical OR Urothelium OR Detrusor OR Cystitis OR Overactive Bladder OR Painful Bladder"
    }
}

# --- FUNÇÕES DE CALLBACK ---
def carregar_setup_lemos():
    st.session_state.alvos_val = LISTA_ALVOS_PRONTA
    st.session_state.fonte_val = PRESETS_ORGAOS["(Sugestão Lemos)"]["fonte"]
    st.session_state.alvo_val = PRESETS_ORGAOS["(Sugestão Lemos)"]["alvo"]
    st.toast("Setup Completo Carregado!", icon="🎓")

def carregar_alvos_apenas(): 
    st.session_state.alvos_val = LISTA_ALVOS_PRONTA
    st.toast("Lista de Alvos Restaurada!", icon="📋")

# ==========================================
# 3. FUNÇÕES TÉCNICAS
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
    st.title("🎓 Lemos Doutorado")
    st.markdown("**Ferramenta de Inteligência Bibliométrica**")

    st.sidebar.header("1. Credenciais")
    email_user = st.sidebar.text_input("Seu E-mail:", placeholder="pesquisador@unifesp.br", key="email_desk")
    anos = st.sidebar.slider("📅 Período:", 1990, 2025, (2010, 2025), key="anos_desk")
    min_year, max_year = anos
    
    st.sidebar.markdown("---")
    st.sidebar.header("2. Configuração")
    
    termo_fonte = st.sidebar.text_input("Fonte (Comparação):", key="fonte_val", placeholder="Sistemas Consolidados...")
    termo_alvo = st.sidebar.text_input("Alvo (Seu Foco):", key="alvo_val", placeholder="Bexiga/Urotélio...")
    
    st.sidebar.caption("👇 Configuração com um clique:")
    
    # BOTÃO PRINCIPAL
    st.sidebar.button("🎓 Doutorado Guilherme Lemos", type="primary", on_click=carregar_setup_lemos)
    
    st.sidebar.markdown("---")
    st.sidebar.header("3. Alvos")
    
    alvos_input = st.sidebar.text_area("Lista de Pesquisa:", key="alvos_val", height=150, placeholder="Carregue a lista...")
    st.sidebar.button("📥 Restaurar Lista Interna", on_click=carregar_alvos_apenas)

    st.sidebar.markdown("---")

    if st.sidebar.button("🚀 INICIAR VARREDURA", type="secondary"):
        if not email_user or "@" not in email_user: 
            st.error("E-mail obrigatório!")
        elif not termo_fonte or not termo_alvo: 
            st.warning("Configure os órgãos!")
        elif not alvos_input: 
            st.warning("Lista vazia!")
        else:
            # LIMPEZA EXTRA ANTES DA BUSCA
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
                    
                    status = "N/A"
                    if n_bexiga >= n_fonte: status = "🔴 Saturado"
                    elif ratio > 10 and n_fonte > 200: status = "💎 Chance de OURO"
                    elif ratio > 5 and n_fonte > 100: status = "🥇 Chance Alta"
                    elif ratio > 3: status = "🥇 Chance Média"
                    elif ratio > 1.5: status = "🥈 Chance Baixa"
                    else: status = "🔴 Saturado"
                    
                    resultados.append({
                        "Alvo": alvo, "Status": status, "Potencial (x)": round(ratio, 1),
                        "Fonte Total": n_fonte, "Bexiga Total": n_bexiga
                    })
                bar.progress((i+1)/len(alvos_lista))
            
            progresso_texto.empty()
            st.session_state['dados_desk'] = pd.DataFrame(resultados).sort_values(by="Potencial (x)", ascending=False)

    if 'dados_desk' in st.session_state:
        df = st.session_state['dados_desk']
        top = df.iloc[0]
        
        total_ouro = len(df[df['Status'].str.contains("OURO")])
        st.success(f"✅ Varredura Concluída! Encontramos **{total_ouro} Chances de Ouro**. O maior destaque é **{top['Alvo']}**.")
        
        col1, col2 = st.columns([2, 1])
        with col1:
            fig = px.bar(df.head(20), x="Alvo", y="Potencial (x)", color="Status", 
                         title="Top 20 Oportunidades (Por Status)", 
                         color_discrete_map={"💎 Chance de OURO": "#00CC96", "🥇 Chance Alta": "#636EFA", "🥇 Chance Média": "#AB63FA", "🥈 Chance Baixa": "#FFA15A", "🔴 Saturado": "#EF553B"})
            st.plotly_chart(fig, use_container_width=True)
            
        with col2:
            st.dataframe(
                df[["Alvo", "Status", "Potencial (x)", "Fonte Total", "Bexiga Total"]]
                .style.applymap(lambda v: 'color: red;' if 'Saturado' in str(v) else ('color: green; font-weight: bold;' if 'OURO' in str(v) else ''), subset=['Status'])
                .hide(axis="index"), use_container_width=True, height=500
            )
            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button("📥 Baixar Planilha Completa", csv, f'analise_lemos_{datetime.now().strftime("%Y%m%d")}.csv', 'text/csv', use_container_width=True)
            
        st.divider()
        st.header("🔎 Raio-X & Tradução")
        
        # LISTA ORDENADA A-Z
        lista_ordenada = sorted(df['Alvo'].unique().tolist())
        sel = st.selectbox("Investigar Alvo (A-Z):", lista_ordenada)
        
        col_btn1, col_btn2 = st.columns([1,4])
        if col_btn1.button("Ler Artigos"):
            with st.spinner("Buscando..."):
                arts = buscar_resumos_detalhados(sel, termo_alvo, email_user, min_year, max_year)
                if not arts: st.info("Zero artigos encontrados na bexiga. (Isso confirma o Status!)")
                else:
                    for a in arts:
                        with st.expander(f"📄 {a['Title']}"):
                            st.write(f"**Revista:** {a['Source']}")
                            st.success(a['Resumo_IA'])
                            st.markdown(f"[Link PubMed](https://pubmed.ncbi.nlm.nih.gov/{a['PMID']})")
        
        if col_btn2.button("🎓 Google Scholar"):
             st.markdown(f"👉 [Ver **{sel} + Bexiga** no Scholar](https://scholar.google.com.br/scholar?q={sel}+AND+bladder)", unsafe_allow_html=True)

elif modo == "Mobile (Pocket)":
    st.title("📱 Lemos Pocket")
    email_mob = st.text_input("📧 E-mail:", placeholder="pesquisador@unifesp.br", key="email_mob")
    
    with st.expander("⚙️ Configurar Busca"):
        anos_mob = st.slider("📅 Anos:", 1990, 2025, (2010, 2025))
        t_fonte_mob = st.text_input("Fonte:", key="fonte_val", placeholder="Fonte...")
        t_alvo_mob = st.text_input("Alvo:", key="alvo_val", placeholder="Alvo...")
        st.button("🎓 Doutorado Guilherme Lemos", key="mob_lemos", type="primary", on_click=carregar_setup_lemos)
        st.markdown("---")
        alvos_mob = st.text_area("Alvos:", key="alvos_val", height=150)
        st.button("📥 Restaurar Lista", key="mob_alvos", on_click=carregar_alvos_apenas)
        
    if st.button("🚀 INICIAR", type="secondary", use_container_width=True):
        if not email_mob: st.error("E-mail necessário")
        else:
            # LIMPEZA EXTRA MOBILE
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
                    elif rat > 10 and nf > 200: stat = "💎 OURO"
                    elif rat > 5: stat = "🥇 ALTA"
                    else: stat = "🥈 MÉDIA"
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
        # LISTA ORDENADA MOBILE
        lista_ordenada_mob = sorted(d['Alvo'].unique().tolist())
        sl = st.selectbox("Ler:", lista_ordenada_mob)
        if st.button("Ler", use_container_width=True):
            with st.spinner("Traduzindo..."):
                as_mob = buscar_resumos_detalhados(sl, t_alvo_mob, email_mob, anos_mob[0], anos_mob[1], limit=3)
                if not as_mob: st.info("Sem artigos!")
                else:
                    for am in as_mob:
                        st.success(f"**{am['Title']}**\n\n{am['Resumo_IA']}")
                        st.write("---")
