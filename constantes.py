# constantes.py
# Versão: 4.0 (Tradução Completa + Chaves Adicionais)

TEXTOS = {
    "pt": {
        # --- GERAL ---
        "titulo_desk": "λ Lemos Lambda: Deep Science Prospector",
        "subtitulo": "Ferramenta de Prospecção Farmacológica Dinâmica",
        "radar_titulo": "📡 Radar Científico (Últimas 24h)",
        "btn_ler": "Ler Artigo",
        
        # --- INPUTS ---
        "label_email": "E-mail do Pesquisador (Obrigatório NCBI)",
        "holder_email": "ex: guilherme@unifesp.br",
        "erro_email": "⚠️ E-mail é obrigatório para usar a API!",
        "label_alvo": "Alvo Principal (Órgão/Doença)",
        "holder_alvo": "ex: Bladder, Urothelium, Detrusor",
        "erro_alvo": "⚠️ Defina o Alvo Principal primeiro!",
        "label_contexto": "Órgão para Comparar (Contexto)",
        "uploader_label": "Importar Lista (CSV/TXT)",
        
        # --- BOTÕES ---
        "btn_auto": "🧠 AUTO-DETECTAR ALVOS & INICIAR",
        "btn_add_preset": "📥 Adicionar Categoria",
        "btn_add_all": "💎 ADICIONAR TODOS OS PRESETS",
        "btn_limpar": "Apagar Termos",
        "btn_executar": "🚀 EXECUTAR ANÁLISE ESTATÍSTICA",
        "btn_voltar": "⬅ Voltar",
        "btn_baixar": "📥 Baixar Relatório (CSV)",
        "btn_investigar": "🔎 Investigar com IA",
        "btn_pubmed": "Ver no PubMed",
        
        # --- CONFIGURAÇÃO ---
        "header_config": "Configuração",
        "expander_presets": "💎 Adicionar Fronteiras (Blue Ocean) & Presets",
        "label_categoria": "Categoria:",
        "popover_manual": "✍️ Adição Manual",
        "holder_manual": "Ex: TMAO, Trehalose",
        "btn_add_manual": "Adicionar Manual",
        "slider_tempo": "Janela de Tempo",
        
        # --- IA ---
        "label_ia_global": "🤖 Habilitar IA Generativa (Global)",
        "expander_ia": "🔑 Configurar Chave API",
        "caption_ia": "Insira sua chave gratuita do Google para resumos inteligentes.",
        "placeholder_key": "Cole sua chave AIza... aqui",
        "link_key": "Obter chave grátis",
        "toggle_ia_curadoria": "✨ Ativar Curadoria por IA (Limpeza)",
        
        # --- FEEDBACK & STATUS ---
        "status_minerando": "Minerando literatura para",
        "toast_atualizado": "Lista atualizada!",
        "sucesso_carregado": "alvos carregados. Clique em EXECUTAR abaixo.",
        "msg_alvos_ok": "alvos prontos para análise.",
        "expander_lista": "👀 Ver/Editar Lista Completa",
        "toast_importado": "Importação concluída",
        "erro_arquivo": "Erro ao ler arquivo.",
        "spinner_analise": "Calibrando estatística (Lemos Lambda v4.0)...",
        "titulo_processando": "## 🧬 Processando Estatística...",
        "spinner_investigando": "Lendo artigos com Inteligência Artificial...",

        # --- RESULTADOS ---
        "resultados": "Relatório de Inteligência",
        "metric_top": "Top Oportunidade",
        "metric_score": "Lambda Score",
        "metric_hits": "Artigos Locais",
        "header_heatmap": "Mapa de Calor",
        "header_leitura": "Leitura Guiada (IA)",
        "label_investigar": "Escolha o alvo:",
        "col_mol": "Alvo Molecular",
        "col_status": "Classificação",
        "col_ratio": "Lambda Score",
        "col_art_alvo": "Hits (Alvo)",
        "col_global": "Hits (Comparação)",

        # --- TAGS (Classificação) ---
        "tag_blue_ocean": "💎 Blue Ocean",
        "tag_ghost": "👻 Fantasma",
        "tag_embryonic": "🌱 Embrionário",
        "tag_neutral": "⚖️ Neutro",
        "tag_gold": "🥇 Ouro",
        "tag_trending": "🚀 Tendência",
        "tag_saturated": "🔴 Saturado",

        # --- RODAPÉ ---
        "footer_citar": "Uso Acadêmico",
        "citar_titulo": "📄 Como Citar",
        "citar_texto": "Lemos, G. (2025). Lemos Lambda v4.0. Zenodo. doi.org/10.5281/zenodo.18092141",
        "apoio_titulo": "Apoie o Projeto (Pix):",
        "apoio_desc": "Ajude a manter o código atualizado ☕"
    },
    "en": {
        "titulo_desk": "λ Lemos Lambda: Deep Science Prospector",
        "subtitulo": "Dynamic Pharmacological Prospecting Tool",
        "radar_titulo": "📡 Scientific Radar (Last 24h)",
        "btn_ler": "Read Article",
        "label_email": "Researcher E-mail (NCBI Required)",
        "holder_email": "ex: researcher@university.edu",
        "erro_email": "⚠️ E-mail is required for API usage!",
        "label_alvo": "Main Target (Organ/Disease)",
        "holder_alvo": "ex: Bladder, Urothelium, Detrusor",
        "erro_alvo": "⚠️ Define Target first!",
        "label_contexto": "Comparison Organ (Context)",
        "uploader_label": "Import List (CSV/TXT)",
        "btn_auto": "🧠 AUTO-DETECT TARGETS",
        "btn_add_preset": "📥 Add Category",
        "btn_add_all": "💎 ADD ALL PRESETS",
        "btn_limpar": "Clear Terms",
        "btn_executar": "🚀 EXECUTE ANALYSIS",
        "btn_voltar": "⬅ Back",
        "btn_baixar": "📥 Download Report (CSV)",
        "btn_investigar": "🔎 Investigate with AI",
        "btn_pubmed": "View on PubMed",
        "header_config": "Configuration",
        "expander_presets": "💎 Blue Ocean Frontiers & Presets",
        "label_categoria": "Category:",
        "popover_manual": "✍️ Manual Input",
        "holder_manual": "Ex: ATP, P2X3",
        "btn_add_manual": "Add Manual",
        "slider_tempo": "Time Range",
        "label_ia_global": "🤖 Enable Generative AI (Global)",
        "expander_ia": "🔑 API Key Config",
        "caption_ia": "Insert your free Google API Key for smart summaries.",
        "placeholder_key": "Paste your AIza... key here",
        "link_key": "Get free key",
        "toggle_ia_curadoria": "✨ Enable AI Curation (Cleaning)",
        "status_minerando": "Mining literature for",
        "toast_atualizado": "List updated!",
        "sucesso_carregado": "targets loaded. Click EXECUTE below.",
        "msg_alvos_ok": "targets ready for analysis.",
        "expander_lista": "👀 View/Edit Full List",
        "toast_importado": "Import successful",
        "erro_arquivo": "Error reading file.",
        "spinner_analise": "Calibrating statistics (Lemos Lambda v4.0)...",
        "titulo_processando": "## 🧬 Processing Statistics...",
        "spinner_investigando": "Reading articles with AI...",
        "resultados": "Intelligence Report",
        "metric_top": "Top Opportunity",
        "metric_score": "Lambda Score",
        "metric_hits": "Local Hits",
        "header_heatmap": "Heatmap",
        "header_leitura": "Guided Reading (AI)",
        "label_investigar": "Choose target:",
        "col_mol": "Molecular Target",
        "col_status": "Classification",
        "col_ratio": "Lambda Score",
        "col_art_alvo": "Hits (Target)",
        "col_global": "Hits (Comparison)",
        "tag_blue_ocean": "💎 Blue Ocean",
        "tag_ghost": "👻 Ghost",
        "tag_embryonic": "🌱 Embryonic",
        "tag_neutral": "⚖️ Neutral",
        "tag_gold": "🥇 Gold",
        "tag_trending": "🚀 Trending",
        "tag_saturated": "🔴 Saturated",
        "footer_citar": "Academic Use",
        "citar_titulo": "📄 How to Cite",
        "citar_texto": "Lemos, G. (2025). Lemos Lambda v4.0. Zenodo.",
        "apoio_titulo": "Support the Project (Pix):",
        "apoio_desc": "Help keep the code updated ☕"
    }
}

PRESETS_FRONTEIRA = {
    "🧪 Chaperonas & Estresse": ["Trehalose", "TMAO", "4-PBA", "Taurine", "Betaine", "Chemical Chaperones", "HSP70", "HSP90", "ER Stress"],
    "🔥 Gasotransmissores": ["Carbon Monoxide (CO)", "HO-1", "CORM-2", "CORM-401", "Hydrogen Sulfide (H2S)", "CBS", "CSE", "GYY4137", "AP39", "Hydrogen Gas (H2)"],
    "⚡ Mecanossensores & Canais": ["Piezo1", "Piezo2", "OSCA1", "TMEM63", "TRPA1", "TRPM8", "TRPV4", "TREK-1", "TREK-2", "TRAAK", "TMEM16A", "P2X3", "P2X7"],
    "🧬 Epigenética & RNAs": ["TET2", "HDAC6 inhibitor", "DNMT1", "MALAT1", "HOTAIR", "miR-29b", "miR-132", "Exosomes"],
    "🥑 Resolução & Lipídios": ["Resolvin D1", "Resolvin E1", "Maresin-1", "Lipoxin A4", "Anandamide (AEA)", "2-AG", "GPR55", "GPR18"],
    "👻 Receptores Órfãos": ["GPR35", "GPR84", "GPR68", "OR51E2", "OR1D2", "TAS2R14", "TAS2R38"]
}

CANDIDATOS_MINERACAO = []
for lista in PRESETS_FRONTEIRA.values():
    CANDIDATOS_MINERACAO.extend(lista)

BLACKLIST_GERAL = ["review", "meta-analysis", "rat", "mice", "human", "study"]
