# constantes.py
# Armazena textos, listas estáticas, traduções e PRESETS.

PRESET_LEMOS = {
    "alvo": "Bladder OR Vesical OR Urothelium OR Detrusor OR Cystitis OR Overactive Bladder",
    "fonte": "Brain OR Kidney OR Liver OR Intestine OR Lung OR Vascular OR Immune System"
}

CANDIDATOS_MINERACAO = [
    # --- FRONTEIRAS MODERNAS (Sci-Fi to Reality) ---
    "Microplastics", "Nanoplastics", "Bisphenol A", "Phthalates", "PFAS", # Ambiental
    "Trimethylamine N-oxide (TMAO)", "Gut microbiota metabolites", "Indole-3-propionic acid", # Microbioma
    "Cannabidiol (CBD)", "Cannabigerol (CBG)", "GPR55", "GPR18", "FAAH enzyme", "MAGL enzyme", "Anandamide", # Cannabis Expandido
    "Senolytics", "Klotho", "Nicotinamide mononucleotide (NMN)", "Sirtuins", # Longevidade
    "Ferroptosis", "Pyroptosis", "Necroptosis", "Exosomes", "Extracellular vesicles", # Morte Celular e Comunicação
    "GLP-1 receptor", "SGLT2 inhibitors", "Probiotics", "Prebiotics", # Metabólico Moderno
    
    # --- Ácidos e Metabólitos Clássicos ---
    "Alpha-lipoic acid", "Butyrate", "Short-chain fatty acids", "Sphingosine-1-phosphate",
    "Hyaluronic acid", "Succinate", "Lactate", "Kynurenic acid", "Prostaglandin E2",
    "Resolvin D1", "Lipoxin A4", "Melatonin", "Taurine", "Uric Acid",
    
    # --- Receptores e Canais ---
    "P2X3 receptor", "P2X7 receptor", "TRPV1 channel", "TRPM8", "Piezo1", "Piezo2",
    "Cannabinoid receptor 2", "GPR119", "GPR40", "GPR84", "GPR120",
    "Beta-3 adrenergic receptor", "Muscarinic M3", "Nicotinic alpha-7",
    "TMEM16A", "HCN channels", "Kv7 channels",
    
    # --- Vias e Enzimas ---
    "Rho-kinase (ROCK)", "mTOR pathway", "AMPK signaling", "NLRP3 inflammasome",
    "Nitric Oxide Synthase", "Heme Oxygenase-1", "Cyclooxygenase-2", "Phosphodiesterase-5",
    "YAP/TAZ pathway", "Hippo pathway",
    
    # --- Genética e RNA ---
    "MALAT1 lncRNA", "miR-21", "miR-145", "SIRT1", "NRF2 pathway", "NF-kappaB"
]

# --- FILTRO DE RUÍDO (BLACKLIST) ---
BLACKLIST_GERAL = [
    # Termos identificados como ruído
    "cross-sectional", "postmarketing", "survey", "adrs", "aes", "ars", "aims",
    "pmid", "authors", "areas covered", "accordingly", "induced",
    "pharma", "pharmaceuticals", "solutions", "abbvie", "aveo", "accord",
    "abu dhabi", "abe", "abbosh", "atlab", "abyost", "aiq solutions", "aikido", "akus-11",
    "apobec3-induced", "april", "asb3",
    
    # Termos de erro/metodologia
    "adverse event", "adverse effect", "ae rate", "safety", "efficacy", "placebo",
    "control group", "study design", "double-blind", "randomized", "clinical trial",
    "p-value", "confidence interval", "odds ratio", "hazard ratio", "standard deviation",
    "anova", "regression", "analysis", "data", "result", "conclusion", "method",
    "significant", "statistically", "increased", "decreased", "compared to",
    "associated with", "observed in", "related to", "due to",
    
    # Termos genéricos demais
    "signaling pathway", "signal transduction", "gene expression", "protein level",
    "messenger rna", "receptor agonist", "receptor antagonist", "inhibitor",
    "mechanism of action", "therapeutic target", "potential target", "biomarker",
    "pathophysiology", "metabolism", "oxidative stress", "inflammation",
    "cell culture", "in vivo", "in vitro", "western blot", "pcr", "elisa",
    "stem cell", "progenitor cell", "expression of", "activation of", "levels of",
    
    # Palavras de ligação/Lixo
    "the", "and", "with", "for", "that", "this", "were", "was", "have", "has",
    "between", "among", "during", "after", "before", "however", "therefore",
    "furthermore", "moreover", "additionally", "notably", "interestingly",
    
    # Institucional
    "department", "university", "hospital", "institute", "center", "usa", "china",
    "brazil", "europe", "funding", "grant", "review", "article", "copyright"
]

TEXTOS = {
    "pt": {
        "titulo_desk": "λ Lemos Lambda: Deep Science",
        "subtitulo": "Ferramenta de Prospecção Farmacológica Dinâmica",
        "step_1": "1️⃣ Definição de Contexto",
        "step_2": "2️⃣ Mineração & Análise",
        "label_email": "E-mail (Obrigatório para PubMed):",
        "holder_email": "ex: pesquisador@unifesp.br",
        "label_alvo": "Alvo Principal (Doença/Órgão):",
        "holder_alvo": "ex: Overactive Bladder, Fibrosis...",
        "aviso_pubmed": "⚠️ **Atenção:** Escreva os termos em **INGLÊS** para garantir a mineração correta.",
        
        "btn_smart_load": "🔄 Carregar Minha Lista (+ Dinâmica)",
        "btn_blue_ocean": "🌊 Explore o Blue Ocean (Apenas Novidades)",
        "btn_lib": "📚 Minerar no Contexto (Fonte)",
        "btn_preset": "🎓 Carregar Preset Doutorado (Lemos)",
        
        "status_blue_ocean": "🌊 Mergulhando no PubMed em busca de alvos inexplorados...",
        "msg_sucesso_blue": "🌊 {qtd} novos alvos do Blue Ocean adicionados!",
        "status_minerando": "🔍 Minerando novidades para:",
        "msg_sucesso_dinamico": "✅ Lista Base + {qtd} novidades específicas preservadas!",
        "msg_sucesso_base": "✅ Apenas Lista Base carregada (Preencha 'Alvo' para torná-la dinâmica).",
        "erro_fonte_vazia": "⚠️ Preencha o campo 'Contexto/Fonte' para buscar novidades nele.",
        "toast_preset": "🧬 Preset Lemos carregado com sucesso!",
        
        "analise_btn": "🚀 Executar Análise de Potencial",
        "resultados": "🎯 Dashboard de Prospecção",
        "label_periodo": "📅 Período de Análise",
        "label_manual": "🔎 Investigar Termo Específico",
        "holder_manual": "ex: Curcumina, Gene X...",
        "btn_add_manual": "➕ Adicionar",
        "toast_add": "✅ termo(s) adicionado(s)!",
        "toast_dup": "⚠️ Duplicatas ignoradas.",
        "label_fonte": "Contexto Comparativo (Opcional):",
        "holder_fonte": "ex: Brain, Kidney, Liver...",
        "desc_fonte": "Define o 'universo' de comparação. Também usado para mineração contextual.",
        "titulo_import": "📂 Importar Lista Extra",
        "desc_import": "Upload (.csv/.txt)",
        "toast_import": "✅ termos importados!",
        "erro_ler": "Erro ao ler arquivo.",
        "btn_limpar": "🗑️",
        "btn_limpar_tudo": "🗑️ Limpar Lista",
        "btn_export_lista": "💾 Salvar Lista (CSV)",
        "ver_editar": "📝 Ver/Editar Lista de Palavras-Chave",
        "qtd_termos": "Qtd:",
        "radar_titulo": "📡 Radar Científico",
        "btn_ler_feed": "Ler Completo",
        "metrica_potencial": "🏆 Maior Potencial",
        "metrica_score": "📊 Score (Ratio)",
        "metrica_artigos": "📚 Artigos (Alvo)",
        "col_mol": "Molécula/Alvo",
        "col_status": "Status",
        "col_ratio": "Potencial (Ratio)",
        "col_art_alvo": "Artigos no Alvo",
        "col_global": "Global/Fonte",
        "btn_baixar": "📥 Baixar Relatório CSV",
        "erro_email": "E-mail necessário para conectar ao NCBI.",
        "erro_campos": "⚠️ Preencha E-mail e Alvo (em Inglês)!",
        "erro_sessao": "⚠️ Dados da sessão expiraram. Por favor, faça a pesquisa novamente.",
        
        "titulo_mapa": "Mapa de Oportunidades",
        "titulo_leitura": "📄 Leitura Profunda: Investigação de Papers",
        "btn_nova_pesquisa": "⬅️ Nova Pesquisa",
        "info_leitura": "Selecione um alvo abaixo para buscar os artigos reais no PubMed com tradução.",
        "sel_leitura": "Selecione o alvo:",
        "btn_buscar_artigos": "🔍 Carregar Artigos sobre:",
        "msg_buscando_lit": "Buscando literatura sobre",
        "header_artigos_enc": "Artigos encontrados:",
        "aviso_sem_artigos": "Nenhum artigo encontrado com resumo disponível neste período.",
        
        "footer_citar": "Lemos Lambda v1.1.0 - Uso Acadêmico",
        "citar_titulo": "📄 Como Citar",
        "citar_texto": "Lemos, G. (2025). Lemos Lambda: Deep Science Prospector [Software]. Versão 1.1.0. DOI: 10.5281/zenodo.17958507",
        "link_doi": "🔗 Ver no Zenodo (DOI)"
    },
    "en": {
        "titulo_desk": "λ Lemos Lambda: Deep Science",
        "subtitulo": "Dynamic Pharmacological Prospecting Tool",
        "step_1": "1️⃣ Context Definition",
        "step_2": "2️⃣ Mining & Analysis",
        "label_email": "E-mail (Required for PubMed):",
        "holder_email": "ex: researcher@university.edu",
        "label_alvo": "Main Target (Disease/Organ):",
        "holder_alvo": "ex: Overactive Bladder, Fibrosis...",
        "aviso_pubmed": "⚠️ **Warning:** Please input terms in **ENGLISH**.",
        
        "btn_smart_load": "🔄 Load My List (+ Dynamic)",
        "btn_blue_ocean": "🌊 Explore Blue Ocean (Novelties Only)",
        "btn_lib": "📚 Context Mining (Source)",
        "btn_preset": "🎓 Load Lemos PhD Preset",
        
        "status_blue_ocean": "🌊 Diving into PubMed for unexplored targets...",
        "msg_sucesso_blue": "🌊 {qtd} new Blue Ocean targets added!",
        "status_minerando": "🔍 Mining novelties for:",
        "msg_sucesso_dinamico": "✅ Base List + {qtd} specific novelties preserved!",
        "msg_sucesso_base": "✅ Base List loaded only (Fill 'Target' to make it dynamic).",
        "erro_fonte_vazia": "⚠️ Fill in 'Context/Source' field to search novelties within it.",
        "toast_preset": "🧬 Lemos Preset loaded successfully!",
        
        "analise_btn": "🚀 Run Potential Analysis",
        "resultados": "🎯 Prospecting Dashboard",
        "label_periodo": "📅 Analysis Period",
        "label_manual": "🔎 Investigate Specific Term",
        "holder_manual": "ex: Curcumin, Gene X...",
        "btn_add_manual": "➕ Add",
        "toast_add": "✅ term(s) added!",
        "toast_dup": "⚠️ Duplicates ignored.",
        "label_fonte": "Comparative Context (Optional):",
        "holder_fonte": "ex: Brain, Kidney, Liver...",
        "desc_fonte": "Defines comparison universe. Also used for contextual mining.",
        "titulo_import": "📂 Import Extra List",
        "desc_import": "Upload (.csv/.txt)",
        "toast_import": "✅ terms imported!",
        "erro_ler": "Error reading file.",
        "btn_limpar": "🗑️",
        "btn_limpar_tudo": "🗑️ Clear List",
        "btn_export_lista": "💾 Save List (CSV)",
        "ver_editar": "📝 View/Edit Keywords List",
        "qtd_termos": "Qty:",
        "radar_titulo": "📡 Science Radar",
        "btn_ler_feed": "Read Full",
        "metrica_potencial": "🏆 Top Potential",
        "metrica_score": "📊 Score (Ratio)",
        "metrica_artigos": "📚 Papers (Target)",
        "col_mol": "Molecule/Target",
        "col_status": "Status",
        "col_ratio": "Potential (Ratio)",
        "col_art_alvo": "Papers on Target",
        "col_global": "Global/Source",
        "btn_baixar": "📥 Download CSV Report",
        "erro_email": "E-mail required.",
        "erro_campos": "⚠️ Fill in E-mail and Target (in English)!",
        "erro_sessao": "⚠️ Session data expired. Please run search again.",
        
        "titulo_mapa": "Opportunity Map",
        "titulo_leitura": "📄 Deep Reading: Paper Investigation",
        "btn_nova_pesquisa": "⬅️ New Search",
        "info_leitura": "Select a target below to fetch real articles from PubMed with translation.",
        "sel_leitura": "Select target:",
        "btn_buscar_artigos": "🔍 Load Articles on:",
        "msg_buscando_lit": "Searching literature for",
        "header_artigos_enc": "Articles found:",
        "aviso_sem_artigos": "No articles found with abstract available in this period.",
        
        "footer_citar": "Lemos Lambda v1.1.0 - Academic Use",
        "citar_titulo": "📄 How to Cite",
        "citar_texto": "Lemos, G. (2025). Lemos Lambda: Deep Science Prospector [Software]. Version 1.1.0. DOI: 10.5281/zenodo.17958507",
        "link_doi": "🔗 View on Zenodo (DOI)"
    }
}