import streamlit as st
import pandas as pd
from Bio import Entrez
import time
import plotly.express as px
import re

# ==========================================
# 1. CONFIGURAÇÃO GLOBAL
# ==========================================
st.set_page_config(
    page_title="Lemos Buscador", 
    page_icon="🧬", 
    layout="wide" # Layout wide funciona bem para ambos (no mobile ele empilha)
)

# ==========================================
# 2. FUNÇÕES COMUNS (O CÉREBRO DO APP)
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
    except:
        return -1

def extrair_conclusao(abstract_text):
    if not abstract_text: return "Resumo não disponível."
    match = re.search(r'(Conclusion|Conclusions|In conclusion|Summary|Results suggest that)(.*)', abstract_text, re.IGNORECASE | re.DOTALL)
    if match: return "➡️ " + match.group(2).strip()[:300] + "..." 
    return "➡️ " + abstract_text[:200] + "..."

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
    except Exception as e: return []

# ==========================================
# 3. SELETOR DE MODO (A MÁGICA)
# ==========================================
# Menu lateral fixo para troca de versão
modo = st.sidebar.radio("📱 Escolha a Versão:", ["Desktop (Completo)", "Mobile (Pocket)"], index=0)
st.sidebar.markdown("---")

# Listas Padrão
lista_padrao = """Autophagy, LC3B (MAP1LC3B), Beclin-1 (BECN1), p62 (SQSTM1), ATG5, mTOR, VEGF, VEGFR1, VEGFR2, TGF-beta1, CTGF, Galectin-3, P2X3, TRPV1, TRPV4, Beta-3 Adrenergic, SGLT2, ROCK (Rho-kinase), NLRP3, IL-17"""
lista_limpa = " ".join(lista_padrao.replace("\n", " ").split())

# ==========================================
# 4. VERSÃO DESKTOP PRO (Lógica V6.0)
# ==========================================
if modo == "Desktop (Completo)":
    st.title("🧬 Lemos Buscador: Desktop Pro")
    st.markdown("**Ferramenta de Inteligência Bibliométrica Avançada**")

    # --- Sidebar Desktop ---
    st.sidebar.header("⚙️ Parâmetros")
    email_user = st.sidebar.text_input("Seu E-mail:", placeholder="pesquisador@unifesp.br", key="email_desk")
    anos = st.sidebar.slider("📅 Período:", 1990, 2025, (2010, 2025), key="anos_desk")
    min_year, max_year = anos
    
    alvos_input = st.sidebar.text_area("Lista de Alvos:", value=lista_limpa, height=250, key="alvos_desk")
    st.sidebar.info("Use vírgulas para separar os alvos.")
    
    termo_fonte = st.sidebar.text_input("Fonte (Comparação):", value="Kidney OR Renal OR Blood Vessels OR Vascular OR Lung OR Gut", key="fonte_desk")
    termo_alvo = st.sidebar.text_input("Alvo (Seu Foco):", value="Bladder OR Vesical OR Urothelium OR Detrusor OR Cystitis", key="alvo_desk")
    
    if st.sidebar.button("🚀 Iniciar Análise Completa", type="primary"):
        if not email_user or "@" not in email_user:
            st.error("E-mail obrigatório!")
        else:
            alvos_lista = [x.strip() for x in alvos_input.split(",") if x.strip()]
            resultados = []
            bar = st.progress(0)
            
            for i, alvo in enumerate(alvos_lista):
                n_fonte = consultar_pubmed_count(alvo, termo_fonte, email_user, min_year, max_year)
                n_bexiga = consultar_pubmed_count(alvo, termo_alvo, email_user, min_year, max_year)
                if n_fonte != -1:
                    ratio = n_fonte / n_bexiga if n_bexiga > 0 else n_fonte
                    resultados.append({"Alvo": alvo, "Fonte Total": n_fonte, "Bexiga Total": n_bexiga, "Potencial": round(ratio, 1)})
                bar.progress((i+1)/len(alvos_lista))
            
            st.session_state['dados_desk'] = pd.DataFrame(resultados).sort_values(by="Potencial", ascending=False)
            st.success("Análise finalizada!")

    # --- Resultados Desktop ---
    if 'dados_desk' in st.session_state:
        df = st.session_state['dados_desk']
        
        # Resumo Inteligente
        top = df.iloc[0]
        st.info(f"💡 **Insight Rápido:** A maior oportunidade encontrada foi **{top['Alvo']}**, que é **{top['Potencial']}x** mais estudado fora da bexiga. Existem {len(df[df['Potencial']>10])} alvos com alto potencial de ineditismo.")
        
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("📊 Gráfico de Nichos")
            fig = px.bar(df.head(15), x="Alvo", y="Potencial", color="Potencial", color_continuous_scale="Bluered")
            st.plotly_chart(fig, use_container_width=True)
        with col2:
            st.subheader("📋 Tabela de Dados")
            st.dataframe(df.style.background_gradient(subset=['Potencial'], cmap="Greens").hide(axis="index"), use_container_width=True, height=400)
            
        st.divider()
        st.header("🔎 Raio-X Detalhado")
        sel_alvo = st.selectbox("Selecione para ver os artigos:", df['Alvo'].tolist())
        
        if st.button("Buscar Artigos (Desktop)"):
            with st.spinner("Lendo conclusões..."):
                artigos = buscar_resumos_detalhados(sel_alvo, termo_alvo, email_user, min_year, max_year)
                if not artigos: st.balloons(); st.success("Nicho Confirmado! Zero artigos encontrados.")
                else:
                    for art in artigos:
                        with st.expander(f"📄 {art['Title']}"):
                            st.write(f"**Fonte:** {art['Source']}")
                            st.info(art['Resumo_IA'])
                            st.caption(art['Abstract'][:300] + "...")
                            st.markdown(f"[Link PubMed](https://pubmed.ncbi.nlm.nih.gov/{art['PMID']})")

# ==========================================
# 5. VERSÃO MOBILE POCKET (Lógica V7.0)
# ==========================================
elif modo == "Mobile (Pocket)":
    st.title("📱 Lemos Pocket")
    st.caption("Interface simplificada para uso em celular.")

    # --- Inputs Mobile (No centro, escondidos em Expander) ---
    email_mobile = st.text_input("📧 E-mail (Obrigatório):", placeholder="pesquisador@unifesp.br", key="email_mob")
    
    with st.expander("⚙️ Configurar Busca (Toque para abrir)"):
        anos_mob = st.slider("📅 Anos:", 1990, 2025, (2010, 2025), key="anos_mob")
        alvos_mob = st.text_area("Alvos:", value=lista_limpa, height=150, key="alvos_mob")
        t_fonte_mob = st.text_input("Fonte:", value="Kidney OR Vascular OR Lung", key="f_mob")
        t_alvo_mob = st.text_input("Alvo:", value="Bladder OR Cystitis", key="a_mob")
    
    if st.button("🚀 INICIAR (Modo Rápido)", type="primary", use_container_width=True):
        if not email_mobile or "@" not in email_mobile:
            st.error("Preencha o e-mail!")
        else:
            alvos_lista = [x.strip() for x in alvos_mob.split(",") if x.strip()]
            resultados = []
            progresso = st.progress(0)
            for i, alvo in enumerate(alvos_lista):
                n_fonte = consultar_pubmed_count(alvo, t_fonte_mob, email_mobile, anos_mob[0], anos_mob[1])
                n_bexiga = consultar_pubmed_count(alvo, t_alvo_mob, email_mobile, anos_mob[0], anos_mob[1])
                if n_fonte != -1:
                    ratio = n_fonte / n_bexiga if n_bexiga > 0 else n_fonte
                    resultados.append({"Alvo": alvo, "Potencial": round(ratio, 1)})
                progresso.progress((i+1)/len(alvos_lista))
            
            st.session_state['dados_mob'] = pd.DataFrame(resultados).sort_values(by="Potencial", ascending=False)
            st.toast("Busca Concluída!", icon="✅")

    # --- Resultados Mobile ---
    if 'dados_mob' in st.session_state:
        df_mob = st.session_state['dados_mob']
        top_mob = df_mob.iloc[0]
        
        st.divider()
        col_a, col_b = st.columns(2)
        col_a.metric("🏆 Top Alvo", top_mob['Alvo'])
        col_b.metric("Potencial", f"{top_mob['Potencial']}x")
        
        with st.expander("📋 Ver Lista Completa"):
            st.dataframe(df_mob, use_container_width=True, hide_index=True)
            
        st.divider()
        st.subheader("🔎 Raio-X Rápido")
        sel_mob = st.selectbox("Escolha o alvo:", df_mob['Alvo'].tolist(), key="sel_mob")
        
        if st.button(f"Ler sobre {sel_mob}", use_container_width=True):
            with st.spinner("Lendo..."):
                arts_mob = buscar_resumos_detalhados(sel_mob, t_alvo_mob, email_mobile, anos_mob[0], anos_mob[1], limit=3)
                if not arts_mob: st.info("Nenhum artigo encontrado!")
                else:
                    for art in arts_mob:
                        st.success(f"**{art['Title']}**\n\n{art['Resumo_IA']}")
                        st.link_button("Abrir PubMed", f"https://pubmed.ncbi.nlm.nih.gov/{art['PMID']}", use_container_width=True)
                        st.write("---")
