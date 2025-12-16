import streamlit as st
import pandas as pd
from Bio import Entrez
import time
import plotly.express as px

# ==========================================
# 1. CONFIGURAÇÃO INICIAL DA PÁGINA
# ==========================================
st.set_page_config(page_title="Pharma-Gap Hunter", page_icon="💊", layout="wide")

st.title("💊 Pharma-Gap Hunter: Detector de Nichos Farmacológicos")
st.markdown("""
**Ferramenta de Inteligência Bibliométrica para Farmacologia**
Compara o volume de publicações de alvos/fármacos em órgãos "Modelo" (Rins/Vasos) 
versus o órgão "Alvo" (Bexiga), identificando oportunidades de reposicionamento.
""")

# ==========================================
# 2. BARRA LATERAL (INPUTS)
# ==========================================
st.sidebar.header("⚙️ Parâmetros de Busca")

# E-mail é obrigatório para a API do NCBI
email_user = st.sidebar.text_input("Seu E-mail (Obrigatório pelo PubMed):", 
                                  value="pesquisador@unifesp.br",
                                  help="O NCBI exige um e-mail para monitorar o uso da API.")

# Lista padrão focada em farmacologia (Canais, Receptores, Enzimas)
lista_sugestao = """
P2X3, TRPV1, TRPV4, Beta-3 Adrenergic, Muscarinic M3, 
SGLT2, mTOR, PDE5, VEGF, NGF, BDNF, COX-2, 
ROCK (Rho-kinase), Cannabinoid CB1, Cannabinoid CB2
"""
lista_sugestao = lista_sugestao.replace("\n", "").strip()

alvos_input = st.sidebar.text_area("Lista de Alvos/Fármacos (separados por vírgula):", 
                                   value=lista_sugestao, height=150)

# Termos de Busca (Editáveis para flexibilidade)
st.sidebar.markdown("---")
st.sidebar.subheader("🔠 Definição dos Termos (Query)")
st.sidebar.info("As buscas são feitas em INGLÊS para máxima abrangência.")

termo_fonte = st.sidebar.text_input("Termos Fonte (Órgãos Análogos):", 
                                    value="Kidney OR Renal OR Nephron OR Blood Vessels OR Vascular OR Endothelial")

termo_alvo = st.sidebar.text_input("Termos Alvo (Seu Foco):", 
                                   value="Bladder OR Vesical OR Urothelium OR Detrusor OR Cystitis")

botao_buscar = st.sidebar.button("🚀 Iniciar Mineração de Dados", type="primary")

# ==========================================
# 3. FUNÇÃO DE CONEXÃO COM PUBMED
# ==========================================
def consultar_pubmed(termo_farmaco, termo_orgao, email):
    Entrez.email = email
    # Monta a query: (Fármaco) AND (Termos do Órgão)
    query_final = f"({termo_farmaco}) AND ({termo_orgao})"
    try:
        # Busca apenas a contagem (retmax=0)
        handle = Entrez.esearch(db="pubmed", term=query_final, retmax=0)
        record = Entrez.read(handle)
        handle.close()
        return int(record["Count"])
    except Exception as e:
        return -1 # Código de erro

# ==========================================
# 4. LÓGICA DE PROCESSAMENTO
# ==========================================
if botao_buscar:
    if not email_user or "@" not in email_user:
        st.error("⚠️ Por favor, insira um e-mail válido antes de continuar.")
    else:
        # Limpar e preparar lista de alvos
        alvos_lista = [x.strip() for x in alvos_input.split(",") if x.strip()]
        
        resultados = []
        progresso = st.progress(0)
        status = st.empty()
        
        total = len(alvos_lista)
        
        for i, alvo in enumerate(alvos_lista):
            status.text(f"🔍 Pesquisando no PubMed: {alvo}...")
            
            # 1. Busca nos Análogos (Fonte)
            n_fonte = consultar_pubmed(alvo, termo_fonte, email_user)
            
            # 2. Busca na Bexiga (Alvo)
            n_bexiga = consultar_pubmed(alvo, termo_alvo, email_user)
            
            if n_fonte != -1 and n_bexiga != -1:
                gap = n_fonte - n_bexiga
                # Ratio: Quantas vezes é mais estudado fora? (Evita divisão por zero)
                ratio = n_fonte / n_bexiga if n_bexiga > 0 else n_fonte
                
                # Classificação automática
                classificacao = "Neutro"
                if ratio > 5 and n_fonte > 100: classificacao = "💎 Oportunidade Alta"
                elif ratio > 2 and n_fonte > 50: classificacao = "🥇 Oportunidade Média"
                elif n_bexiga > n_fonte: classificacao = "🔴 Já Saturado na Bexiga"

                resultados.append({
                    "Alvo Molecular": alvo,
                    "Hits (Rins/Vasos)": n_fonte,
                    "Hits (Bexiga)": n_bexiga,
                    "Gap Absoluto": gap,
                    "Potencial (Ratio)": round(ratio, 1),
                    "Status": classificacao
                })
            
            # Atualiza barra de progresso
            progresso.progress((i + 1) / total)
            time.sleep(0.34) # Delay para respeitar limite de 3 req/s do NCBI sem API Key

        status.text("✅ Análise Bibliométrica Concluída!")
        
        # Cria DataFrame
        df = pd.DataFrame(resultados)
        
        if not df.empty:
            # Ordenar por Potencial
            df = df.sort_values(by="Potencial (Ratio)", ascending=False)

            # ==========================================
            # 5. DASHBOARD DE RESULTADOS
            # ==========================================
            
            # Métricas Principais
            top_nicho = df.iloc[0]
            col1, col2, col3 = st.columns(3)
            col1.metric("Maior Nicho Encontrado", top_nicho['Alvo Molecular'])
            col2.metric("Dominância Externa", f"{top_nicho['Potencial (Ratio)']}x maior")
            col3.metric("Publicações Totais Analisadas", df['Hits (Rins/Vasos)'].sum() + df['Hits (Bexiga)'].sum())
            
            st.divider()

            # Gráfico
            st.subheader("📊 Comparativo de Maturidade da Pesquisa")
            df_long = pd.melt(df, id_vars=['Alvo Molecular'], 
                              value_vars=['Hits (Rins/Vasos)', 'Hits (Bexiga)'],
                              var_name='Origem', value_name='Artigos')
            
            fig = px.bar(df_long, x="Alvo Molecular", y="Artigos", color="Origem", barmode='group',
                         color_discrete_map={"Hits (Rins/Vasos)": "#FF4B4B", "Hits (Bexiga)": "#1E90FF"})
            st.plotly_chart(fig, use_container_width=True)

            # Tabela Interativa
            st.subheader("📋 Dados Detalhados")
            st.dataframe(df.style.background_gradient(subset=['Potencial (Ratio)'], cmap="Greens"), use_container_width=True)
            
            # Botão de Download
            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button("📥 Baixar Relatório (CSV)", data=csv, file_name='analise_nichos_doutorado.csv', mime='text/csv')
        else:
            st.warning("Nenhum dado retornado. Verifique sua conexão.")