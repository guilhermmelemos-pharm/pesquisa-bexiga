import streamlit as st
import pandas as pd
from Bio import Entrez
import time
import plotly.express as px

# ==========================================
# 1. CONFIGURAÇÃO INICIAL DA PÁGINA
# ==========================================
st.set_page_config(
    page_title="Lemos Buscador", 
    page_icon="🧬", 
    layout="wide"
)

st.title("🧬 Lemos Buscador")
st.markdown("""
**Ferramenta de Inteligência Bibliométrica**
1. Insira seu e-mail (obrigatório).
2. O sistema fará a varredura e gerará um **Relatório de Inteligência** automático.
3. No final, use o **Raio-X** para ver os artigos reais.
""")

# ==========================================
# 2. BARRA LATERAL (INPUTS)
# ==========================================
st.sidebar.header("⚙️ Parâmetros de Pesquisa")

email_user = st.sidebar.text_input("Seu E-mail (Obrigatório):", 
                                  value="", placeholder="ex: pesquisador@unifesp.br")

# LISTA MESTRA (Mantivemos a mesma lista completa)
lista_sugestao = """
-- AUTOFAGIA --
Autophagy, LC3B (MAP1LC3B), Beclin-1 (BECN1), p62 (SQSTM1), 
ATG5, ATG7, ULK1, LAMP2, TFEB, AMPK, mTOR,

-- FATORES DE CRESCIMENTO & FIBROSE --
VEGF, VEGFR1, VEGFR2, NRP1 (Neuropilin), VEGF-B,
TGF-beta1, CTGF, Galectin-3, MMP-9, NGF, BDNF,

-- CANAIS IÔNICOS & RECEPTORES --
P2X3, P2X7, TRPV1, TRPV4, BK channel, Kv7.4, SK3, 
Piezo1, Piezo2, Beta-3 Adrenergic, Muscarinic M3,
Cannabinoid CB1, Cannabinoid CB2,

-- ENZIMAS, INFLAMAÇÃO & OUTROS --
SGLT2, PDE5, ROCK (Rho-kinase), ACE2, Angiotensin II,
COX-2, NLRP3, IL-17, TLR4, Nrf2, PPAR-gamma
"""

lista_limpa = lista_sugestao.replace("\n", " ").replace("-- AUTOFAGIA --", "").replace("-- FATORES DE CRESCIMENTO & FIBROSE --", "").replace("-- CANAIS IÔNICOS & RECEPTORES --", "").replace("-- ENZIMAS, INFLAMAÇÃO & OUTROS --", "")
lista_limpa = " ".join(lista_limpa.split())

alvos_input = st.sidebar.text_area("Lista de Alvos:", value=lista_limpa, height=300)

st.sidebar.markdown("---")
st.sidebar.info("Buscas em INGLÊS com operadores booleanos.")

termo_fonte = st.sidebar.text_input("Termos Fonte (Modelos):", 
                                    value="Kidney OR Renal OR Blood Vessels OR Vascular OR Lung OR Airway OR Intestine OR Gut OR Diabetic Nephropathy")

termo_alvo = st.sidebar.text_input("Termos Alvo (Bexiga):", 
                                   value="Bladder OR Vesical OR Urothelium OR Detrusor OR Cystitis OR Painful Bladder OR Overactive Bladder")

botao_buscar = st.sidebar.button("🚀 Iniciar Lemos Buscador", type="primary")

# ==========================================
# 3. FUNÇÕES (PUBMED E ANÁLISE)
# ==========================================
def consultar_pubmed_count(termo_farmaco, termo_orgao, email):
    Entrez.email = email
    termo_farmaco = termo_farmaco.replace(",", "").strip()
    query = f"({termo_farmaco}) AND ({termo_orgao})"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle)
        return int(record["Count"])
    except:
        return -1

def buscar_resumos_bexiga(termo_farmaco, termo_orgao, email, max_results=5):
    Entrez.email = email
    termo_farmaco = termo_farmaco.replace(",", "").strip()
    query = f"({termo_farmaco}) AND ({termo_orgao})"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
        record = Entrez.read(handle)
        id_list = record["IdList"]
        
        if not id_list: return []
            
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        records = handle.read()
        
        artigos = []
        raw_articles = records.split("\n\n")
        
        for art in raw_articles:
            lines = art.split("\n")
            title = "Sem Título"
            source = "Fonte desconhecida"
            pmid = "N/A"
            for line in lines:
                if line.startswith("TI  - "): title = line[6:]
                if line.startswith("TA  - "): source = line[6:]
                if line.startswith("PMID- "): pmid = line[6:]
            if pmid != "N/A":
                artigos.append({"PMID": pmid, "Título": title, "Revista": source})
        return artigos
    except Exception as e:
        return [{"Erro": str(e)}]

def gerar_analise_textual(df):
    """ Gera um resumo inteligente baseado nos dados """
    top_nicho = df.iloc[0]
    total_nichos = len(df[df['Potencial'] > 10])
    
    texto = f"""
    ### 🧠 Análise Automática
    O algoritmo varreu **{len(df)} alvos farmacológicos**.
    
    **1. O Grande Destaque:**
    A maior oportunidade detectada foi para **{top_nicho['Alvo']}**. 
    Este alvo é **{top_nicho['Potencial']} vezes mais estudado** nos modelos comparativos (Rim/Vaso/Pulmão) do que na Bexiga.
    Isso indica uma maturidade científica alta em outras áreas, mas um terreno quase virgem no seu campo.
    
    **2. Volume de Oportunidades:**
    Encontramos **{total_nichos} alvos classificados como 'Nichos de Ouro'** (Ratio > 10x). 
    Esses são os candidatos ideais para reposicionamento imediato.
    
    **3. Sugestão de Próximo Passo:**
    Recomendamos focar a leitura nos resumos de **{top_nicho['Alvo']}** (usando a ferramenta abaixo) para verificar se os poucos artigos existentes ({top_nicho['Bexiga Total']}) já cobriram o mecanismo que você deseja propor.
    """
    return texto

# ==========================================
# 4. LÓGICA PRINCIPAL
# ==========================================
if botao_buscar:
    # --- BLOQUEIO DE E-MAIL ---
    if not email_user or "@" not in email_user or len(email_user) < 5:
        st.error("⛔ PARE! O preenchimento do E-mail é obrigatório para acessar o PubMed.")
        st.stop() # Para a execução aqui se não tiver e-mail
    
    else:
        alvos_lista = [x.strip() for x in alvos_input.split(",") if x.strip()]
        resultados = []
        progresso = st.progress(0)
        total = len(alvos_lista)
        
        for i, alvo in enumerate(alvos_lista):
            n_fonte = consultar_pubmed_count(alvo, termo_fonte, email_user)
            n_bexiga = consultar_pubmed_count(alvo, termo_alvo, email_user)
            
            if n_fonte != -1:
                ratio = n_fonte / n_bexiga if n_bexiga > 0 else n_fonte
                resultados.append({
                    "Alvo": alvo,
                    "Fonte Total": n_fonte,
                    "Bexiga Total": n_bexiga,
                    "Potencial": round(ratio, 1)
                })
            progresso.progress((i + 1) / total)
            time.sleep(0.1)

        st.session_state['dados'] = pd.DataFrame(resultados).sort_values(by="Potencial", ascending=False)
        st.success("Varredura concluída!")

# --- PARTE 2: EXIBIÇÃO ---
if 'dados' in st.session_state:
    df = st.session_state['dados']
    
    # 1. RESUMO DE INTELIGÊNCIA (NOVIDADE)
    st.divider()
    with st.container():
        st.markdown("## 📝 Resumo de Inteligência")
        st.info(gerar_analise_textual(df))
    
    st.divider()

    # 2. Gráfico e Tabela
    col_chart, col_table = st.columns([1, 1])
    
    with col_chart:
        st.subheader("📊 Ranking de Oportunidade")
        fig = px.bar(df.head(15), x="Alvo", y="Potencial", color="Potencial", 
                     title="Top 15 Nichos (Ratio Fonte/Bexiga)", color_continuous_scale="Bluered")
        st.plotly_chart(fig, use_container_width=True)
        
    with col_table:
        st.subheader("📋 Dados Brutos")
        st.dataframe(df.style.background_gradient(subset=['Potencial'], cmap="Greens").hide(axis="index"), 
                     use_container_width=True, height=400)

    st.divider()

    # 3. RAIO-X
    st.header("🔎 Raio-X do Nicho: Validação")
    st.markdown("Selecione um alvo para ler os resumos:")
    
    lista_alvos = df['Alvo'].tolist()
    alvo_selecionado = st.selectbox("Alvo:", lista_alvos)
    
    if st.button(f"Buscar Artigos sobre {alvo_selecionado}"):
        with st.spinner(f"Lemos Buscador investigando {alvo_selecionado}..."):
            artigos = buscar_resumos_bexiga(alvo_selecionado, termo_alvo, email_user)
            
            if not artigos:
                st.balloons()
                st.success(f"💎 CONFIRMADO! Zero artigos encontrados para '{alvo_selecionado}' na bexiga.")
            else:
                st.warning(f"Atenção: Já existem {len(artigos)} artigos principais. Verifique se não saturaram o tema.")
                for art in artigos:
                    with st.expander(f"📄 {art.get('Título', 'Sem Título')}"):
                        st.write(f"**Revista:** {art.get('Revista', 'N/A')}")
                        st.markdown(f"[Ler no PubMed](https://pubmed.ncbi.nlm.nih.gov/{art.get('PMID', '')})")
