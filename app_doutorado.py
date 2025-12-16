import streamlit as st
import pandas as pd
from Bio import Entrez
import time
import plotly.express as px
import re
from deep_translator import GoogleTranslator

# ==========================================
# 1. CONFIGURAÇÃO GLOBAL
# ==========================================
st.set_page_config(
    page_title="Lemos Buscador", 
    page_icon="🇧🇷", 
    layout="wide" 
)

# ==========================================
# 2. FUNÇÕES DO SISTEMA
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

def traduzir_para_pt(texto):
    """ Traduz do Inglês para Português usando Deep Translator """
    try:
        tradutor = GoogleTranslator(source='auto', target='pt')
        return tradutor.translate(texto)
    except:
        return texto # Se der erro, devolve em inglês mesmo

def extrair_conclusao(abstract_text):
    if not abstract_text: return "Resumo não disponível."
    
    # Tenta achar a conclusão em inglês
    match = re.search(r'(Conclusion|Conclusions|In conclusion|Summary|Results suggest that)(.*)', abstract_text, re.IGNORECASE | re.DOTALL)
    
    texto_final = ""
    if match: 
        texto_final = match.group(2).strip()[:400] # Pega até 400 caracteres
    else:
        # Se não achar a palavra mágica, pega as últimas 3 frases
        frases = abstract_text.split(". ")
        if len(frases) > 3:
            texto_final = ". ".join(frases[-3:])
        else:
            texto_final = abstract_text[:300]
            
    # TRADUZ A CONCLUSÃO ENCONTRADA
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
                # Aqui chamamos a tradução
                art_data["Resumo_IA"] = extrair_conclusao(art_data["Abstract"])
                artigos.append(art_data)
        return artigos
    except Exception as e: return []

# ==========================================
# 3. INTERFACE E MODOS
# ==========================================
modo = st.sidebar.radio("📱 Modo de Visualização:", ["Desktop (Completo)", "Mobile (Pocket)"], index=0)
st.sidebar.markdown("---")

# LISTA COMPLETA (SUGESTÃO)
lista_sugestao_completa = """
Autophagy, LC3B (MAP1LC3B), Beclin-1 (BECN1), p62 (SQSTM1), ATG5, mTOR, 
VEGF, VEGFR1, VEGFR2, NRP1 (Neuropilin), VEGF-B, 
TGF-beta1, CTGF, Galectin-3, MMP-9, 
P2X3, TRPV1, TRPV4, Beta-3 Adrenergic, Muscarinic M3, 
SGLT2, ROCK (Rho-kinase), NLRP3, IL-17, Nrf2
"""
lista_limpa = " ".join(lista_sugestao_completa.replace("\n", " ").split())

# ==========================================
# 4. VERSÃO DESKTOP
# ==========================================
if modo == "Desktop (Completo)":
    st.title("🧬 Lemos Buscador: Versão Tradutor")
    st.markdown("**Ferramenta de Inteligência Bibliométrica com IA de Tradução**")

    # Sidebar
    st.sidebar.header("⚙️ Parâmetros")
    email_user = st.sidebar.text_input("Seu E-mail:", placeholder="pesquisador@unifesp.br", key="email_desk")
    anos = st.sidebar.slider("📅 Período:", 1990, 2025, (2010, 2025), key="anos_desk")
    min_year, max_year = anos
    
    # LÓGICA DO BOTÃO DE SUGESTÃO (DESKTOP)
    if "alvos_desk_val" not in st.session_state: st.session_state.alvos_desk_val = ""
    
    def carregar_sugestao_desk():
        st.session_state.alvos_desk_val = lista_limpa

    st.sidebar.markdown("### 🎯 Alvos Moleculares")
    # Botão que carrega a lista
    st.sidebar.button("📥 Carregar Sugestões do Doutorado", on_click=carregar_sugestao_desk)
    
    # Caixa de texto vinculada ao Session State (começa vazia, enche se clicar no botão)
    alvos_input = st.sidebar.text_area("Lista para Pesquisa:", key="alvos_desk_val", height=250)
    
    termo_fonte = st.sidebar.text_input("Fonte (Comparação):", value="Kidney OR Renal OR Blood Vessels OR Vascular OR Lung OR Gut", key="fonte_desk")
    termo_alvo = st.sidebar.text_input("Alvo (Seu Foco):", value="Bladder OR Vesical OR Urothelium OR Detrusor OR Cystitis", key="alvo_desk")
    
    if st.sidebar.button("🚀 Iniciar Análise", type="primary"):
        if not email_user or "@" not in email_user:
            st.error("E-mail obrigatório!")
        elif not alvos_input:
            st.warning("A lista de alvos está vazia! Digite algo ou clique em 'Carregar Sugestões'.")
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

    # Resultados Desktop
    if 'dados_desk' in st.session_state:
        df = st.session_state['dados_desk']
        top = df.iloc[0]
        st.info(f"💡 **Insight:** O maior nicho é **{top['Alvo']}** ({top['Potencial']}x).")
        
        col1, col2 = st.columns(2)
        with col1:
            fig = px.bar(df.head(15), x="Alvo", y="Potencial", color="Potencial", color_continuous_scale="Bluered")
            st.plotly_chart(fig, use_container_width=True)
        with col2:
            st.dataframe(df.style.background_gradient(subset=['Potencial'], cmap="Greens").hide(axis="index"), use_container_width=True, height=400)
            
        st.divider()
        st.header("🔎 Raio-X (Com Tradução PT-BR)")
        st.caption("O sistema vai ler o abstract em inglês e traduzir a conclusão para você.")
        
        sel_alvo = st.selectbox("Selecione o alvo:", df['Alvo'].tolist())
        
        if st.button("Buscar e Traduzir Artigos"):
            with st.spinner(f"Lendo e traduzindo artigos sobre {sel_alvo}... (Isso pode levar alguns segundos)"):
                artigos = buscar_resumos_detalhados(sel_alvo, termo_alvo, email_user, min_year, max_year)
                if not artigos: st.balloons(); st.success("Nicho Confirmado! Zero artigos encontrados.")
                else:
                    for art in artigos:
                        with st.expander(f"📄 {art['Title']}"):
                            st.write(f"**Fonte:** {art['Source']}")
                            # Exibe a tradução em destaque
                            st.success(art['Resumo_IA'])
                            st.caption(f"Original (En): {art['Abstract'][:200]}...")
                            st.markdown(f"[Link PubMed](https://pubmed.ncbi.nlm.nih.gov/{art['PMID']})")

# ==========================================
# 5. VERSÃO MOBILE
# ==========================================
elif modo == "Mobile (Pocket)":
    st.title("📱 Lemos Pocket Tradutor")
    
    email_mobile = st.text_input("📧 E-mail:", placeholder="pesquisador@unifesp.br", key="email_mob")
    
    # LÓGICA DO BOTÃO DE SUGESTÃO (MOBILE)
    if "alvos_mob_val" not in st.session_state: st.session_state.alvos_mob_val = ""
    def carregar_sugestao_mob(): st.session_state.alvos_mob_val = lista_limpa

    with st.expander("⚙️ Configurar Busca"):
        anos_mob = st.slider("📅 Anos:", 1990, 2025, (2010, 2025), key="anos_mob")
        st.button("📥 Usar Sugestões Padrão", on_click=carregar_sugestao_mob, key="btn_mob_sug")
        alvos_mob = st.text_area("Alvos:", key="alvos_mob_val", height=150)
        t_fonte_mob = st.text_input("Fonte:", value="Kidney OR Vascular OR Lung", key="f_mob")
        t_alvo_mob = st.text_input("Alvo:", value="Bladder OR Cystitis", key="a_mob")
    
    if st.button("🚀 INICIAR", type="primary", use_container_width=True):
        if not email_mobile or "@" not in email_mobile: st.error("Preencha o e-mail!")
        elif not alvos_mob: st.warning("Lista vazia!")
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

    if 'dados_mob' in st.session_state:
        df_mob = st.session_state['dados_mob']
        top_mob = df_mob.iloc[0]
        
        st.divider()
        col_a, col_b = st.columns(2)
        col_a.metric("🏆 Top 1", top_mob['Alvo'])
        col_b.metric("Potencial", f"{top_mob['Potencial']}x")
        
        with st.expander("📋 Ver Lista Completa"):
            st.dataframe(df_mob, use_container_width=True, hide_index=True)
            
        st.divider()
        st.subheader("🔎 Raio-X Traduzido")
        sel_mob = st.selectbox("Escolha:", df_mob['Alvo'].tolist(), key="sel_mob")
        
        if st.button(f"Ler sobre {sel_mob}", use_container_width=True):
            with st.spinner("Traduzindo..."):
                arts_mob = buscar_resumos_detalhados(sel_mob, t_alvo_mob, email_mobile, anos_mob[0], anos_mob[1], limit=3)
                if not arts_mob: st.info("Nenhum artigo!")
                else:
                    for art in arts_mob:
                        st.success(f"**{art['Title']}**\n\n{art['Resumo_IA']}")
                        st.link_button("Ver Original", f"https://pubmed.ncbi.nlm.nih.gov/{art['PMID']}", use_container_width=True)
                        st.write("---")
