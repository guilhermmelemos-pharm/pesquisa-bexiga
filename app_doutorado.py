import streamlit as st
import pandas as pd
from Bio import Entrez
import time
import plotly.express as px

# ==========================================
# 1. CONFIGURAÇÃO INICIAL DA PÁGINA
# ==========================================
st.set_page_config(
    page_title="Lemos Buscador 1.0", 
    page_icon="🧬", 
    layout="wide"
)

st.title("🧬 Lemos Buscador")
st.markdown("""
**Ferramenta de Inteligência Bibliométrica**
Análise unificada de oportunidades farmacológicas, cruzando dados de 
**Rins, Vasos, Pulmão e Intestino** contra o **Baixo Trato Urinário**.
""")

# ==========================================
# 2. BARRA LATERAL (INPUTS)
# ==========================================
st.sidebar.header("⚙️ Parâmetros de Pesquisa")

# E-mail é obrigatório para a API do NCBI
email_user = st.sidebar.text_input("Seu E-mail (Obrigatório pelo PubMed):", 
                                  value="pesquisador@unifesp.br",
                                  help="O NCBI exige um e-mail para monitorar o uso da API.")

# --- LISTA MESTRA (TUDO INCLUÍDO) ---
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

# Limpeza da string para o formato de busca
lista_limpa = lista_sugestao.replace("\n", " ").replace("-- AUTOFAGIA --", "").replace("-- FATORES DE CRESCIMENTO & FIBROSE --", "").replace("-- CANAIS IÔNICOS & RECEPTORES --", "").replace("-- ENZIMAS, INFLAMAÇÃO & OUTROS --", "")
lista_limpa = " ".join(lista_limpa.split()) # Remove espaços duplos

alvos_input = st.sidebar.text_area("Lista de Alvos (Editável):", 
                                   value=lista_limpa, height=300)

st.sidebar.markdown("---")
st.sidebar.subheader("🔠 Termos de Busca (Query)")
st.sidebar.info("Buscas em INGLÊS com operadores booleanos.")

# --- TERMOS FONTE (INCLUI RIM, VASO, PULMÃO, INTESTINO) ---
termo_fonte = st.sidebar.text_input("Termos Fonte (Modelos Comparativos):", 
                                    value="Kidney OR Renal OR Blood Vessels OR Vascular OR Lung OR Airway OR Intestine OR Gut OR Diabetic Nephropathy OR Hypertension")

# --- TERMOS ALVO (BEXIGA E DOENÇAS) ---
termo_alvo = st.sidebar.text_input("Termos Alvo (Seu Foco):", 
                                   value="Bladder OR Vesical OR Urothelium OR Detrusor OR Cystitis OR Painful Bladder OR Overactive Bladder OR Lower Urinary Tract")

botao_buscar = st.sidebar.button("🚀 Iniciar Varredura Completa", type="primary")

# ==========================================
# 3. FUNÇÃO DE CONEXÃO COM PUBMED
# ==========================================
def consultar_pubmed(termo_farmaco, termo_orgao, email):
    Entrez.email = email
    # Remove caracteres especiais que podem quebrar a busca
    termo_farmaco_limpo = termo_farmaco.replace(",", "").strip()
    
    # Monta a query
    query_final = f"({termo_farmaco_limpo}) AND ({termo_orgao})"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query_final, retmax=0)
        record = Entrez.read(handle)
        handle.close()
        return int(record["Count"])
    except Exception as e:
        return -1 

# ==========================================
# 4. LÓGICA DE PROCESSAMENTO
# ==========================================
if botao_buscar:
    if not email_user or "@" not in email_user:
        st.error("⚠️ Por favor, insira um e-mail válido.")
    else:
        # Separa a lista por vírgulas e remove vazios
        alvos_lista = [x.strip() for x in alvos_input.split(",") if x.strip()]
        
        resultados = []
        progresso = st.progress(0)
        status = st.empty()
        total = len(alvos_lista)
        
        for i, alvo in enumerate(alvos_lista):
            status.markdown(f"🔍 Analisando: **{alvo}**...")
            
            # Buscas
            n_fonte = consultar_pubmed(alvo, termo_fonte, email_user)
            n_bexiga = consultar_pubmed(alvo, termo_alvo, email_user)
            
            if n_fonte != -1 and n_bexiga != -1:
                gap = n_fonte - n_bexiga
                # Ratio: Quantas vezes é mais estudado nos outros órgãos?
                ratio = n_fonte / n_bexiga if n_bexiga > 0 else n_fonte
                
                # Classificação de Nicho
                classificacao = "Neutro"
                if ratio > 15 and n_fonte > 300: classificacao = "💎 NICHO DE OURO"
                elif ratio > 5 and n_fonte > 100: classificacao = "🥇 Oportunidade Alta"
                elif ratio > 2 and n_fonte > 50: classificacao = "🥈 Oportunidade Média"
                elif n_bexiga >= n_fonte and n_bexiga > 50: classificacao = "🔴 Saturado na Bexiga"

                resultados.append({
                    "Alvo Molecular": alvo,
                    "Hits (Fonte Total)": n_fonte,
                    "Hits (Bexiga)": n_bexiga,
                    "Gap Absoluto": gap,
                    "Potencial (Ratio)": round(ratio, 1),
                    "Status": classificacao
                })
            
            progresso.progress((i + 1) / total)
            time.sleep(0.35) # Delay importante para lista grande

        status.success("✅ Varredura Completa Finalizada!")
        
        df = pd.DataFrame(resultados)
        
        if not df.empty:
            df = df.sort_values(by="Potencial (Ratio)", ascending=False)

            st.divider()
            
            # KPIs
            top_nicho = df.iloc[0]
            col1, col2, col3 = st.columns(3)
            col1.metric("Maior Nicho Encontrado", top_nicho['Alvo Molecular'])
            col2.metric("Potencial (x vezes)", f"{top_nicho['Potencial (Ratio)']}x")
            col3.metric("Total Artigos (Fonte)", top_nicho['Hits (Fonte Total)'])
            
            # Gráfico
            st.subheader("📊 Comparativo de Oportunidades")
            
            df_long = pd.melt(
                df, 
                id_vars=['Alvo Molecular'], 
                value_vars=['Hits (Fonte Total)', 'Hits (Bexiga)'],
                var_name='Origem', 
                value_name='Artigos'
            )
            
            fig = px.bar(
                df_long, 
                x="Alvo Molecular", 
                y="Artigos", 
                color="Origem", 
                barmode='group',
                color_discrete_map={
                    "Hits (Fonte Total)": "#FF4B4B", 
                    "Hits (Bexiga)": "#1E90FF"
                },
                height=600
            )
            
            st.plotly_chart(fig, use_container_width=True)

            # Tabela Limpa (Sem índice numérico)
            st.subheader("📋 Dados Detalhados")
            st.dataframe(
                df.style.background_gradient(subset=['Potencial (Ratio)'], cmap="Greens")
                        .format({"Potencial (Ratio)": "{:.1f}", "Hits (Fonte Total)": "{:,.0f}"})
                        .hide(axis="index"),
                use_container_width=True
            )
            
            # Download
            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button(
                label="📥 Baixar Tabela Completa (CSV)", 
                data=csv, 
                file_name='analise_completa_bexiga.csv', 
                mime='text/csv'
            )
        else:
            st.warning("Sem dados retornados. Verifique a conexão.")
