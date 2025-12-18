import streamlit as st
import pandas as pd
import plotly.express as px
import backend as bk
import constantes as c

# --- CONFIGURAÇÃO DA PÁGINA ---
st.set_page_config(
    page_title="Lemos Lambda v1.6",
    page_icon="λ",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- ESTILOS CSS (MODO ESCURO/CIENTÍFICO) ---
st.markdown("""
<style>
    .stApp {
        background-color: #0E1117;
        color: #FAFAFA;
    }
    .stMetric {
        background-color: #262730;
        border: 1px solid #4F4F4F;
        padding: 10px;
        border-radius: 5px;
    }
    h1, h2, h3 {
        color: #00B4D8; /* Azul Lemos */
    }
    .stButton>button {
        width: 100%;
        border-radius: 5px;
        font-weight: bold;
    }
    /* Destaque para Blue Ocean */
    .blue-ocean {
        color: #48CAE4;
        font-weight: bold;
    }
</style>
""", unsafe_allow_html=True)

# --- INICIALIZAÇÃO DE SESSÃO ---
if 'idioma' not in st.session_state: st.session_state.idioma = 'pt'
if 'termos_extras' not in st.session_state: st.session_state.termos_extras = []
t = c.TEXTOS[st.session_state.idioma]

# --- SIDEBAR (CONTROLES) ---
with st.sidebar:
    st.header(f"⚙️ {t['step_1']}")
    
    # 1. Dados de Entrada
    email = st.text_input(t['label_email'], placeholder=t['holder_email'], key="email_input")
    alvo = st.text_input(t['label_alvo'], placeholder=t['holder_alvo'], key="alvo_input")
    fonte = st.text_input(t['label_fonte'], value="Immune System", placeholder=t['holder_fonte'])
    
    st.info(t['aviso_pubmed'])
    
    st.markdown("---")
    
    # 2. Seleção de Modo
    modo_busca = st.radio("Modo de Mineração:", 
                         [t['btn_preset'], t['btn_blue_ocean'], t['btn_smart_load']],
                         index=0)
    
    st.markdown("---")
    
    # 3. Filtros Temporais
    ano_atual = 2026
    periodo = st.slider(t['label_periodo'], 2015, ano_atual, (ano_atual-5, ano_atual))
    
    # 4. Botão Principal
    executar = st.button(t['analise_btn'], type="primary")
    
    # Rodapé
    st.markdown("---")
    st.caption(t['footer_citar'])

# --- CORPO PRINCIPAL ---
st.title(t['titulo_desk'])
st.subheader(f"🧬 {t['subtitulo']}")

if executar:
    if not email or not alvo:
        st.error(t['erro_campos'])
    else:
        # 1. Definir a Lista de Candidatos
        candidatos = []
        
        with st.status(t['status_minerando'] + f" {alvo}...") as status:
            if modo_busca == t['btn_preset']:
                candidatos = c.CANDIDATOS_MINERACAO
                status.write("✅ Preset Lemos Carregado.")
                
            elif modo_busca == t['btn_smart_load']:
                status.write("🧠 Minerando o PubMed por termos associados...")
                novos = bk.buscar_alvos_emergentes_pubmed(alvo, email)
                candidatos = list(set(c.CANDIDATOS_MINERACAO + novos))
                status.write(f"✅ {len(novos)} novos termos encontrados!")
                
            elif modo_busca == t['btn_blue_ocean']:
                status.write(t['status_blue_ocean'])
                # No Blue Ocean, focamos 100% na Enciclopédia Sci-Fi
                candidatos = c.CANDIDATOS_MINERACAO
                status.write("🌊 Modo Descoberta Ativado.")
            
            # Adiciona manuais se houver
            if st.session_state.termos_extras:
                candidatos += st.session_state.termos_extras
            
            # Remove duplicatas e limpa
            candidatos = sorted(list(set(candidatos)))
            
            # 2. Processamento (Loop de Análise)
            dados = []
            progresso = st.progress(0)
            total = len(candidatos)
            
            for i, termo in enumerate(candidatos):
                # Busca contagens
                n_alvo = bk.consultar_pubmed_count(termo, alvo, email, periodo[0], periodo[1])
                n_fonte = bk.consultar_pubmed_count(termo, fonte, email, periodo[0], periodo[1]) # Contexto Global
                
                # Cálculo do Ratio (Potencial)
                # Se tem muito na fonte (global) e pouco no alvo = OPORTUNIDADE (Blue Ocean)
                if n_alvo == 0 and n_fonte > 50:
                    ratio = 100.0 # Blue Ocean Puro
                    status_termo = "💎 Blue Ocean"
                elif n_alvo > 0:
                    ratio = (n_fonte / n_alvo) if n_alvo > 0 else 0
                    if n_alvo < 10: status_termo = "🌱 Embrionário"
                    elif ratio < 2: status_termo = "🔴 Saturado"
                    else: status_termo = "🚀 Tendência"
                else:
                    ratio = 0
                    status_termo = "⚪ Sem Dados"
                
                if n_fonte > 0 or n_alvo > 0: # Só adiciona se existir na ciência
                    dados.append({
                        "Termo": termo,
                        "Alvo_Count": n_alvo,
                        "Fonte_Count": n_fonte,
                        "Ratio": round(ratio, 2),
                        "Status": status_termo
                    })
                
                progresso.progress((i + 1) / total)
            
            status.update(label="Análise Concluída!", state="complete", expanded=False)

        # 3. Visualização dos Resultados
        if dados:
            df_results = pd.DataFrame(dados)
            
            # --- CORREÇÃO DO KEYERROR: RENOMEAR COLUNAS PARA O IDIOMA ---
            # Isso garante que o DataFrame tenha as chaves que o 'constantes.py' espera
            mapa_colunas = {
                "Termo": t["col_mol"],         # "Molécula/Alvo"
                "Status": t["col_status"],     # "Status"
                "Ratio": t["col_ratio"],       # "Potencial (Ratio)"
                "Alvo_Count": t["col_art_alvo"], # "Artigos no Alvo"
                "Fonte_Count": t["col_global"]   # "Global/Fonte"
            }
            # Renomeia e reordena
            df_results = df_results.rename(columns=mapa_colunas)
            cols_order = [t["col_mol"], t["col_status"], t["col_ratio"], t["col_art_alvo"], t["col_global"]]
            df_results = df_results[cols_order].sort_values(by=t["col_ratio"], ascending=False)

            # Dashboard Top
            c1, c2, c3 = st.columns(3)
            
            # Pega o Top 1
            if not df_results.empty:
                top = df_results.iloc[0]
                # Agora as chaves batem porque renomeamos o DF acima!
                c1.metric(t["metrica_potencial"], top[t["col_mol"]], delta=top[t["col_status"]])
                c2.metric(t["metrica_score"], top[t["col_ratio"]])
                c3.metric(t["metrica_artigos"], top[t["col_art_alvo"]])

            st.markdown("---")
            
            # Gráfico de Bolhas (Scatter Plot)
            # Adaptado para usar os novos nomes de colunas
            fig = px.scatter(
                df_results, 
                x=t["col_art_alvo"], 
                y=t["col_ratio"], 
                size=t["col_global"], 
                color=t["col_status"],
                hover_name=t["col_mol"],
                log_x=True, 
                size_max=60,
                color_discrete_map={
                    "💎 Blue Ocean": "#00B4D8",
                    "🌱 Embrionário": "#90E0EF",
                    "🚀 Tendência": "#48CAE4",
                    "🔴 Saturado": "#F72585",
                    "⚪ Sem Dados": "#6C757D"
                },
                title=t["radar_titulo"]
            )
            fig.update_layout(paper_bgcolor="#0E1117", plot_bgcolor="#0E1117", font_color="white")
            st.plotly_chart(fig, use_container_width=True)
            
            # Tabela Interativa
            st.dataframe(
                df_results.style.background_gradient(subset=[t["col_ratio"]], cmap="Blues"),
                use_container_width=True,
                hide_index=True
            )
            
            # Download
            csv = df_results.to_csv(index=False).encode('utf-8')
            st.download_button(
                label=t["btn_baixar"],
                data=csv,
                file_name=f'lemos_lambda_{alvo}.csv',
                mime='text/csv',
            )
            
            # --- ÁREA DE LEITURA (PAPER FETCH) ---
            st.markdown("---")
            st.subheader(t["titulo_leitura"])
            st.info(t["info_leitura"])
            
            # Selectbox usa a coluna renomeada
            termo_leitura = st.selectbox(t["sel_leitura"], df_results[t["col_mol"]].unique())
            
            if st.button(f"{t['btn_buscar_artigos']} {termo_leitura}"):
                with st.spinner(f"{t['msg_buscando_lit']} {termo_leitura}..."):
                    artigos = bk.buscar_resumos_detalhados(termo_leitura, alvo, email, periodo[0], periodo[1], st.session_state.idioma)
                    
                    if artigos:
                        st.success(f"{t['header_artigos_enc']} {len(artigos)}")
                        for art in artigos:
                            with st.expander(f"📄 {art['Title']}"):
                                st.markdown(f"**PMID:** [{art['PMID']}]({art['Link']}) | **Journal:** {art['Source']}")
                                st.markdown(f"**🤖 Resumo IA:** {art['Resumo_IA']}")
                                st.caption("Original Abstract: " + art['Abstract'][:300] + "...")
                    else:
                        st.warning(t['aviso_sem_artigos'])

        else:
            st.warning("Nenhum dado encontrado com os filtros atuais.")