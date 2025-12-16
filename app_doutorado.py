import streamlit as st
import pandas as pd
from Bio import Entrez
import time
import plotly.express as px

# ==========================================
# 1. CONFIGURAÇÃO INICIAL DA PÁGINA
# ==========================================
st.set_page_config(
    page_title="Lemos Buscador 1.2", 
    page_icon="🧬", 
    layout="wide"
)

st.title("🧬 Lemos Buscador 1.2")
st.markdown("""
**Ferramenta de Inteligência Bibliométrica**
Focada em identificar oportunidades de reposicionamento farmacológico, comparando 
sistemas ricos em dados (Renal, Vascular, Respiratório) com a Bexiga.
""")

# ==========================================
# 2. BARRA LATERAL (INPUTS)
# ==========================================
st.sidebar.header("⚙️ Parâmetros de Pesquisa")

# E-mail é obrigatório para a API do NCBI
email_user = st.sidebar.text_input("Seu E-mail (Obrigatório pelo PubMed):", 
                                  value="pesquisador@unifesp.br",
                                  help="O NCBI exige um e-mail para monitorar o uso da API.")

# --- ATUALIZAÇÃO: AUTOFAGIA + ALVOS CLÁSSICOS ---
lista_sugestao = """
Autophagy, LC3B (MAP1LC3B), Beclin-1 (BECN1), p62 (SQSTM1), 
ATG5, ATG7, AMPK, mTOR, ULK1, LAMP2, TFEB,
VEGF, VEGFR2, TGF-beta1, CTGF, Galectin-3, 
P2X3, TRPV1, TRPV4, Beta-3 Adrenergic, Muscarinic M3, 
SGLT2, PDE5, ROCK (Rho-kinase), Nrf2, NLRP3
"""
# Limpeza da string
lista_sugestao = lista_sugestao.replace("\n", " ").strip()

alvos_input = st.sidebar.text_area("Lista de Alvos/Proteínas (separados por vírgula):", 
                                   value=lista_sugestao, height=250)

st.sidebar.markdown("---")
st.sidebar.subheader("🔠 Termos de Busca (Query)")
st.sidebar.info("Buscas em INGLÊS com operadores booleanos.")

# --- ATUALIZAÇÃO: INCLUSÃO DE PULMÃO E INTESTINO (MÚSCULO LISO) ---
# Adicionei Lung/Airway (Asma) e Gut (Intestino) pois são modelos fortes de músculo liso
termo_fonte = st.sidebar.text_input("Termos Fonte (Modelos Comparativos):", 
                                    value="Kidney OR Renal OR Blood Vessels OR Vascular OR Lung OR Airway OR Intestine OR Gut OR Diabetic Nephropathy")

# Termos Alvo (Bexiga)
termo_alvo = st.sidebar.text_input("Termos Alvo (Seu Foco):", 
                                   value="Bladder OR Vesical OR Urothelium OR Detrusor OR Cystitis OR Painful Bladder OR Overactive Bladder")

botao_buscar = st.sidebar.button("🚀 Iniciar Varredura de Autofagia", type="primary")

# ==========================================
# 3. FUNÇÃO DE CONEXÃO COM PUBMED
# ==========================================
def consultar_pubmed(termo_farmaco, termo_orgao, email):
    Entrez.email = email
    # Limpeza de parênteses (ex: "LC3B (MAP1LC3B)" vira "LC3B MAP1LC3B") para a busca funcionar melhor
    termo_farmaco_limpo = termo_farmaco.replace(",", "")
    
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
                if ratio > 15 and n_fonte > 300: classificacao = "💎 NICHO DE OURO (Autofagia?)"
                elif ratio > 5 and n_fonte > 100: classificacao = "🥇 Oportunidade Alta"
                elif ratio > 2 and n_fonte > 50: classificacao = "🥈 Oportunidade Média"
                elif n_bexiga >= n_fonte: classificacao = "🔴 Saturado na Bexiga"

                resultados.append({
                    "Alvo/Proteína": alvo,
                    "Hits (Fonte: Rim/Vaso/Pulmão)": n_fonte,
                    "Hits (Alvo: Bexiga)": n_bexiga,
                    "Gap Absoluto": gap,
                    "Potencial (Ratio)": round(ratio, 1),
                    "Status": classificacao
                })
            
            progresso.progress((i + 1) / total)
            time.sleep(0.35) 

        status.success("✅ Análise Finalizada!")
        
        df = pd.DataFrame(resultados)
        
        if not df.empty:
            df = df.sort_values(by="Potencial (Ratio)", ascending=False)

            st.divider()
            
            # KPIs
            top_nicho = df.iloc[0]
            col1, col2, col3 = st.columns(3)
            col1.metric("Maior Nicho", top_nicho['Alvo/Proteína'])
            col2.metric("Potencial (x vezes mais estudado fora)", f"{top_nicho['Potencial (Ratio)']}x")
            col3.metric("Total Artigos (Fonte)", top_nicho['Hits (Fonte: Rim/Vaso/Pulmão)'])
            
            # Gráfico
            st.subheader("📊 Comparativo de Maturidade (Autofagia & Outros)")
            
            df_long = pd.melt(
                df, 
                id_vars=['Alvo/Proteína'], 
                value_vars=['Hits (Fonte: Rim/Vaso/Pulmão)', 'Hits (Alvo: Bexiga)'],
                var_name='Origem', 
                value_name='Artigos'
            )
            
            fig = px.bar(
                df_long, 
                x="Alvo/Proteína", 
                y="Artigos", 
                color="Origem", 
                barmode='group',
                color_discrete_map={
                    "Hits (Fonte: Rim/Vaso/Pulmão)": "#FF4B4B", 
                    "Hits (Alvo: Bexiga)": "#1E90FF"
                },
                height=550
            )
            
            st.plotly_chart(fig, use_container_width=True)

            # Tabela Limpa
            st.subheader("📋 Dados Detalhados")
            st.dataframe(
                df.style.background_gradient(subset=['Potencial (Ratio)'], cmap="Greens")
                        .format({"Potencial (Ratio)": "{:.1f}", "Hits (Fonte: Rim/Vaso/Pulmão)": "{:,.0f}"})
                        .hide(axis="index"),
                use_container_width=True
            )
            
            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button(
                label="📥 Baixar Dados (CSV)", 
                data=csv, 
                file_name='analise_autofagia_bexiga.csv', 
                mime='text/csv'
            )
        else:
            st.warning("Sem dados. Verifique a conexão.")
