import streamlit as st
import pandas as pd
from Bio import Entrez
import time
import plotly.express as px

# ==========================================
# 1. CONFIGURAÇÃO INICIAL DA PÁGINA
# ==========================================
st.set_page_config(
    page_title="Pharma-Gap Hunter Pro", 
    page_icon="🧬", 
    layout="wide"
)

st.title("🧬 Pharma-Gap Hunter: Refinado")
st.markdown("""
**Ferramenta de Inteligência Bibliométrica Avançada**
Identifica nichos de pesquisa comparando a maturidade de alvos farmacológicos 
entre sistemas consolidados (Renal/Vascular) e o Baixo Trato Urinário.
""")

# ==========================================
# 2. BARRA LATERAL (INPUTS)
# ==========================================
st.sidebar.header("⚙️ Parâmetros de Pesquisa")

# E-mail é obrigatório para a API do NCBI
email_user = st.sidebar.text_input("Seu E-mail (Obrigatório pelo PubMed):", 
                                  value="pesquisador@unifesp.br",
                                  help="O NCBI exige um e-mail para monitorar o uso da API.")

# Lista Expandida com sugestões da Família VEGF e outros
lista_sugestao = """
VEGF, VEGFR1, VEGFR2, NRP1 (Neuropilin), VEGF-B,
P2X3, P2X7, TRPV1, TRPV4, Beta-3 Adrenergic, Muscarinic M3, 
SGLT2, mTOR, PDE5, NGF, BDNF, COX-2, 
ROCK (Rho-kinase), Cannabinoid CB1, Cannabinoid CB2,
ACE2, Bradykinin B1, Bradykinin B2
"""
lista_sugestao = lista_sugestao.replace("\n", " ").strip()

alvos_input = st.sidebar.text_area("Lista de Alvos/Fármacos (separados por vírgula):", 
                                   value=lista_sugestao, height=200)

st.sidebar.markdown("---")
st.sidebar.subheader("🔠 Termos de Busca (Query)")
st.sidebar.info("As buscas utilizam operadores booleanos (OR/AND) em INGLÊS.")

# Termos Fonte (Órgãos Análogos)
termo_fonte = st.sidebar.text_input("Termos Fonte (Órgãos Análogos):", 
                                    value="Kidney OR Renal OR Nephron OR Blood Vessels OR Vascular OR Endothelial OR Diabetic Nephropathy OR Hypertension")

# Termos Alvo (Seu Foco - Bexiga e Patologias)
termo_alvo = st.sidebar.text_input("Termos Alvo (Seu Foco):", 
                                   value="Bladder OR Vesical OR Urothelium OR Detrusor OR Cystitis OR Painful Bladder OR Overactive Bladder OR Lower Urinary Tract")

botao_buscar = st.sidebar.button("🚀 Iniciar Mineração Refinada", type="primary")

# ==========================================
# 3. FUNÇÃO DE CONEXÃO COM PUBMED
# ==========================================
def consultar_pubmed(termo_farmaco, termo_orgao, email):
    Entrez.email = email
    # Remove parênteses extras do termo do fármaco para evitar erro de sintaxe
    termo_farmaco_limpo = termo_farmaco.replace("(", "").replace(")", "")
    
    # Monta a query: (Fármaco) AND (Termos do Órgão)
    query_final = f"({termo_farmaco_limpo}) AND ({termo_orgao})"
    
    try:
        # Busca apenas a contagem (retmax=0)
        handle = Entrez.esearch(db="pubmed", term=query_final, retmax=0)
        record = Entrez.read(handle)
        handle.close()
        return int(record["Count"])
    except Exception as e:
        return -1 # Retorna -1 se der erro

# ==========================================
# 4. LÓGICA DE PROCESSAMENTO
# ==========================================
if botao_buscar:
    if not email_user or "@" not in email_user:
        st.error("⚠️ Por favor, insira um e-mail válido antes de continuar.")
    else:
        # Limpar e preparar lista de alvos (remove espaços vazios)
        alvos_lista = [x.strip() for x in alvos_input.split(",") if x.strip()]
        
        resultados = []
        progresso = st.progress(0)
        status = st.empty()
        
        total = len(alvos_lista)
        
        for i, alvo in enumerate(alvos_lista):
            status.markdown(f"🔍 Investigando: **{alvo}**...")
            
            # 1. Busca nos Análogos (Fonte)
            n_fonte = consultar_pubmed(alvo, termo_fonte, email_user)
            
            # 2. Busca na Bexiga (Alvo)
            n_bexiga = consultar_pubmed(alvo, termo_alvo, email_user)
            
            if n_fonte != -1 and n_bexiga != -1:
                gap = n_fonte - n_bexiga
                # Ratio: Quantas vezes é mais estudado fora? (Evita divisão por zero)
                ratio = n_fonte / n_bexiga if n_bexiga > 0 else n_fonte
                
                # Classificação automática de oportunidade
                classificacao = "Neutro"
                if ratio > 10 and n_fonte > 200: classificacao = "💎 NICHO DE OURO (Inexplorado)"
                elif ratio > 3 and n_fonte > 100: classificacao = "🥇 Oportunidade Alta"
                elif ratio > 1 and n_fonte > 50: classificacao = "🥈 Oportunidade Média"
                elif n_bexiga >= n_fonte and n_bexiga > 50: classificacao = "🔴 Saturado / Bem Estudado"

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
            time.sleep(0.35) # Delay crucial para a API do PubMed não bloquear

        status.success("✅ Varredura Bibliométrica Concluída com Sucesso!")
        
        # Cria DataFrame
        df = pd.DataFrame(resultados)
        
        if not df.empty:
            # Ordenar por Potencial (Ratio) do maior para o menor
            df = df.sort_values(by="Potencial (Ratio)", ascending=False)

            # ==========================================
            # 5. DASHBOARD DE RESULTADOS
            # ==========================================
            
            st.divider()
            
            # KPIs Principais
            top_nicho = df.iloc[0]
            col1, col2, col3, col4 = st.columns(4)
            col1.metric("Maior Nicho", top_nicho['Alvo Molecular'])
            col2.metric("Potencial de Transferência", f"{top_nicho['Potencial (Ratio)']}x")
            col3.metric("Total Artigos (Fonte)", top_nicho['Hits (Rins/Vasos)'])
            col4.metric("Total Artigos (Alvo)", top_nicho['Hits (Bexiga)'])
            
            # Gráfico Comparativo
            st.subheader("📊 Disparidade de Pesquisa (Gap Analysis)")
            df_long = pd.melt(df, id_vars=['Alvo Molecular'], 
                              value_vars=['Hits (Rins/Vasos)', 'Hits (Bexiga)'],
                              var_name='Origem', value_name='Artigos')
            
            fig = px.bar(df_long, x="Alvo Molecular", y="Artigos", color="Orig
