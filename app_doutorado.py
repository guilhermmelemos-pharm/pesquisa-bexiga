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
**Ferramenta de Inteligência Bibliométrica Personalizada**
1. Identifique o nicho numérico na tabela.
2. **Selecione o alvo no final da página** para ler os resumos do que já foi publicado na Bexiga.
""")

# ==========================================
# 2. BARRA LATERAL (INPUTS)
# ==========================================
st.sidebar.header("⚙️ Parâmetros de Pesquisa")

email_user = st.sidebar.text_input("Seu E-mail (Obrigatório):", 
                                  value="pesquisador@unifesp.br")

# LISTA MESTRA
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
# 3. FUNÇÕES PUBMED (CONTAGEM E DETALHES)
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
    """ Busca os detalhes (Título/Abstract) dos artigos encontrados na Bexiga """
    Entrez.email = email
    termo_farmaco = termo_farmaco.replace(",", "").strip()
    query = f"({termo_farmaco}) AND ({termo_orgao})"
    
    try:
        # 1. Pega os IDs dos artigos
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
        record = Entrez.read(handle)
        id_list = record["IdList"]
        
        if not id_list:
            return []
            
        # 2. Baixa os detalhes desses IDs
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        records = handle.read()
        
        # Processamento simples do texto retornado
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

# ==========================================
# 4. LÓGICA PRINCIPAL
# ==========================================
if botao_buscar:
    if not email_user or "@" not in email_user:
        st.error("⚠️ Insira um e-mail válido.")
    else:
        # --- PARTE 1: A VARREDURA NUMÉRICA ---
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

        # Salva no Session State
        st.session_state['dados'] = pd.DataFrame(resultados).sort_values(by="Potencial", ascending=False)
        st.success("Varredura concluída!")

# --- PARTE 2: EXIBIÇÃO E RAIO-X ---
if 'dados' in st.session_state:
    df = st.session_state['dados']
    
    # 1. Gráfico e Tabela
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

    # --- RAIO-X DE ARTIGOS ---
    st.header("🔎 Raio-X do Nicho: O que já existe?")
    st.markdown("Selecione um alvo abaixo para ver os **5 artigos mais relevantes** publicados sobre ele na Bexiga.")
    
    lista_alvos = df['Alvo'].tolist()
    alvo_selecionado = st.selectbox("Selecione o Alvo para investigar:", lista_alvos)
    
    if st.button(f"Buscar Artigos sobre {alvo_selecionado} na Bexiga"):
        with st.spinner(f"O Lemos Buscador está lendo o PubMed sobre {alvo_selecionado}..."):
            # Busca os detalhes
            artigos = buscar_resumos_bexiga(alvo_selecionado, termo_alvo, email_user)
            
            if not artigos:
                st.balloons()
                st.success(f"💎 NICHO CONFIRMADO! Nenhum artigo encontrado sobre '{alvo_selecionado}' com os termos atuais.")
            else:
                st.write(f"Foram encontrados {len(artigos)} artigos principais. Confira se o tema já está saturado:")
                for art in artigos:
                    with st.expander(f"📄 {art.get('Título', 'Sem Título')}"):
                        st.write(f"**Revista:** {art.get('Revista', 'N/A')}")
                        st.write(f"**PMID:** {art.get('PMID', 'N/A')}")
                        st.markdown(f"[Ler no PubMed](https://pubmed.ncbi.nlm.nih.gov/{art.get('PMID', '')})")
