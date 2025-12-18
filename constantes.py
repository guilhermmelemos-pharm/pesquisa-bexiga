# constantes.py
# Armazena textos, listas estáticas, traduções e PRESETS.

PRESET_LEMOS = {
    "alvo": "Bladder OR Vesical OR Urothelium OR Detrusor OR Cystitis OR Overactive Bladder",
    "fonte": "Brain OR Kidney OR Liver OR Intestine OR Lung OR Vascular OR Immune System"
}

# --- A ENCICLOPÉDIA DE INOVAÇÃO (LIMPA v1.6) ---
# Aqui ficam apenas os termos que você QUER ver.
CANDIDATOS_MINERACAO = [
    # --- 1. ESPECÍFICOS DA TESE (Trehalose & Musculatura Lisa) ---
    "Trehalose", "Autophagy flux", "TFEB (Transcription Factor EB)", 
    "Smooth muscle contraction", "Calcium sensitization", "Myosin phosphatase",
    "Rho-kinase (ROCK)", "mTOR pathway", "AMPK signaling",
    
    # --- 2. SCI-FI & BIOFÍSICA ---
    "Liquid-liquid phase separation", "Biomolecular condensates", "Stress granules", "P-bodies",
    "Cryptochromes", "Magnetoreception", "Quantum biology", "Radical pair mechanism",
    "Bioelectricity", "Resting membrane potential (Vmem)", "Gap junctional communication",
    "Optogenetics", "Sonogenetics", "Magnetogenetics", "Thermogenetics",
    "Tensegrity", "Nuclear mechanotransduction", "Lamin A/C", "Focal adhesions",
    
    # --- 3. RECEPTORES ECTÓPICOS & SENSORIAIS ---
    "Opsins", "Melanopsin (OPN4)", "Encephalopsin (OPN3)", "Neuropsin (OPN5)",
    "Olfactory receptor 51E2 (OR51E2)", "Olfactory receptor 2AT4 (OR2AT4)", "OR10J5", "OR1D2",
    "Bitter taste receptors (TAS2Rs)", "TAS2R14", "Sweet taste receptors (TAS1R2/TAS1R3)", "Umami receptor (TAS1R1)",
    "Mrgprs", "Mas-related G-protein coupled receptors", "MrgprX2", "MrgprD",
    "Piezo1", "Piezo2", "TMEM16A", "TMEM63", "OSCA/TMEM63 family",
    
    # --- 4. EXPOSSOMA & MICROPLÁSTICOS ---
    "Microplastics", "Nanoplastics", "Polystyrene nanoparticles", "PS-MPs", "Polyethylene microbeads",
    "Bisphenol A", "Phthalates", "PFAS", "Perfluorooctanoic acid", "Teflon breakdown products",
    "Endocrine disruptors", "Obesogens", "Glyphosate", "Airborne particulate matter (PM2.5)",
    
    # --- 5. DARK GENOME & EPIGENÉTICA ---
    "Retrotransposons", "LINE-1 elements", "Alu repeats", "Endogenous retroviruses",
    "G-quadruplexes", "R-loops", "Extrachromosomal circular DNA",
    "Circular RNA (circRNA)", "Piwi-interacting RNA (piRNA)", "Enhancer RNAs (eRNA)", "Super-enhancers",
    "MALAT1 lncRNA", "HOTAIR lncRNA", "H19 lncRNA", "NEAT1 lncRNA",
    "Histone lactylation", "Histone crotonylation", "Histone succinylation",
    
    # --- 6. MICROBIOMA & METABÓLITOS ---
    "Trimethylamine N-oxide (TMAO)", "Indole-3-propionic acid", "Urolithin A", "Equol",
    "Short-chain fatty acids", "Butyrate", "Propionate", "Acetate", "Valerate",
    "Secondary bile acids", "Lithocholic acid", "Deoxycholic acid",
    "Akkermansia muciniphila", "Faecalibacterium prausnitzii", "Outer membrane vesicles",
    "Peptidoglycan fragments", "Lipopolysaccharide (LPS)",
    
    # --- 7. LONGEVIDADE & MITOCÔNDRIA ---
    "Senolytics", "Senomorphics", "SASPy", "Klotho", "GDF11", "GDF15",
    "Nicotinamide mononucleotide (NMN)", "Nicotinamide riboside", "NAD+ metabolism",
    "Mitokines", "Humanin", "MOTS-c", "FGF21", 
    "Mitophagy", "Pink1/Parkin pathway", "Mitochondrial unfolded protein response",
    
    # --- 8. MORTE CELULAR & IMUNIDADE NOVA ---
    "Ferroptosis", "GPX4", "Lipid peroxidation",
    "Pyroptosis", "Gasdermin D", "NLRP3 inflammasome",
    "Necroptosis", "RIPK1", "RIPK3", "MLKL",
    "Cuproptosis", "Copper metabolism", "Parthanatos",
    "ILC3s (Innate Lymphoid Cells)", "ZG16 protein", "Galectin-1", "TREM1",
    
    # --- 9. CANNABIS & SINALIZADORES ---
    "Endocannabinoidome", "Anandamide", "2-Arachidonoylglycerol (2-AG)",
    "Oleoylethanolamide (OEA)", "Palmitoylethanolamide (PEA)",
    "Cannabidiol (CBD)", "Cannabigerol (CBG)", "Cannabinol (CBN)", "Tetrahydrocannabivarin (THCV)",
    "GPR55", "GPR18", "GPR119", "GPR110", "GPR120 (FFAR4)",
    "Resolvins", "Protectins", "Maresins", "Lipoxins", 
    
    # --- 10. CANAIS & GPCRs CLÁSSICOS ---
    "P2X3 receptor", "P2X7 receptor", "Purinergic signaling", "ATP release",
    "TRPV1", "TRPV4", "TRPM8", "TRPA1", "ASIC channels",
    "TREK-1 channel", "TRAAK channel", "HCN channels",
    "Muscarinic M3", "Beta-3 adrenergic receptor", "Nicotinic alpha-7", "PAC1 receptor",
    "Nitric Oxide Synthase", "Heme Oxygenase-1", "Hydrogen sulfide (H2S)",
    "YAP/TAZ pathway", "Hippo pathway", "WNT4 signaling"
]

# --- BLACKLIST GERAL (O Lixo morre aqui) ---
BLACKLIST_GERAL = [
    # Siglas Quebradas e Incompletas (Erro de Extração)
    "tgf-", "nf-", "il-", "rxr-", "ppar-", "tnf-", "tlr4", "ca2", "cd320", 
    "lps-induced", "cyp-induced", "sv-huc-1", "rna-seq", "pd-l1", "oab",
    
    # Locais Específicos que Vazaram
    "lublin", "berlin", "sakyo-ku", "hwasun-gun", "jeonnam-do", "sun yat-sen",
    "gustave roussy", "institut curie", "chongqing", "jiangsu", "heidelberg",
    "china", "usa", "japan", "germany", "uk", "france", "italy", "canada",
    "beijing", "shanghai", "guangzhou", "wuhan", "london", "boston", "new york",
    
    # Empresas e Editoras
    "bristol myers squibb", "elsevier", "roche", "mdpi", "springer", "wiley",
    
    # Revistas e Termos Genéricos
    "j mol sci", "medical science", "pediatrics", "genetics", "covid-19",
    "medline", "indexed", "electronic", "epub", "print", "pmid", "doi",
    "background", "methods", "results", "discussion", "conclusion", "abstract",
    "introduction", "references", "acknowledgements", "declaration", "conflict",
    "interest", "funding", "availability", "contributed", "author", "editor",
    
    # Stopwords Acadêmicas
    "however", "moreover", "furthermore", "additionally", "interestingly",
    "significantly", "respectively", "associated", "observed", "indicated"
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
        
        "btn_preset": "🎓 Guilherme Lemos Preset",
        "btn_smart_load": "🔍 Buscar com base no seu Alvo",
        "btn_blue_ocean": "🌊 EXPLORAR BLUE OCEAN (DESCOBERTA)",
        
        "btn_lib": "📚 Minerar no Contexto (Fonte)",
        
        "status_blue_ocean": "🌊 Mergulhando no PubMed em busca de Sci-Fi, Genes e Alvos Ocultos...",
        "msg_sucesso_blue": "🌊 {qtd} tesouros do Blue Ocean adicionados!",
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
        
        "footer_citar": "Lemos Lambda v1.6.0 - Uso Acadêmico",
        "citar_titulo": "📄 Como Citar",
        "citar_texto": "Lemos, G. (2025). Lemos Lambda: Deep Science Prospector [Software]. Versão 1.6.0. DOI: 10.5281/zenodo.17958507",
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
        
        "btn_preset": "🎓 Guilherme Lemos Preset",
        "btn_smart_load": "🔍 Search based on your Target",
        "btn_blue_ocean": "🌊 EXPLORE BLUE OCEAN (DISCOVERY)",
        
        "btn_lib": "📚 Context Mining (Source)",
        "btn_preset": "🎓 Load Lemos PhD Preset",
        
        "status_blue_ocean": "🌊 Diving into PubMed for unexplored targets and Sci-Fi...",
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
        
        "footer_citar": "Lemos Lambda v1.6.0 - Academic Use",
        "citar_titulo": "📄 How to Cite",
        "citar_texto": "Lemos, G. (2025). Lemos Lambda: Deep Science Prospector [Software]. Version 1.6.0. DOI: 10.5281/zenodo.17958507",
        "link_doi": "🔗 View on Zenodo (DOI)"
    }
}