# constantes.py
# Armazena textos, listas estáticas, traduções e PRESETS.

PRESET_LEMOS = {
    "alvo": "Bladder OR Vesical OR Urothelium OR Detrusor OR Cystitis OR Overactive Bladder",
    "fonte": "Brain OR Kidney OR Liver OR Intestine OR Lung OR Vascular OR Immune System"
}

# --- SUA LISTA ÚNICA (ENCICLOPÉDIA LEMOS) ---
CANDIDATOS_MINERACAO = [
    # --- 1. FRONTEIRAS AMBIENTAIS & EXPOSSOMA ---
    "Microplastics", "Nanoplastics", "Bisphenol A", "Phthalates", "PFAS", "Polystyrene nanoparticles",
    "Endocrine disruptors", "Heavy metals", "Glyphosate", "Airborne particulate matter",
    
    # --- 2. MICROBIOMA & METABÓLITOS ---
    "Trimethylamine N-oxide (TMAO)", "Indole-3-propionic acid", "Short-chain fatty acids", 
    "Butyrate", "Propionate", "Acetate", "Secondary bile acids", "Urolithin A",
    "Gut microbiota metabolites", "Lipopolysaccharide (LPS)", "Peptidoglycan",
    "Probiotics", "Prebiotics", "Postbiotics", "Akkermansia muciniphila",
    
    # --- 3. SISTEMA ENDOCANABINOIDE EXPANDIDO ---
    "Cannabidiol (CBD)", "Cannabigerol (CBG)", "Cannabinol (CBN)", "Anandamide", "2-Arachidonoylglycerol (2-AG)",
    "Cannabinoid receptor 1", "Cannabinoid receptor 2", "GPR55", "GPR18", "GPR119", 
    "FAAH enzyme", "MAGL enzyme", "N-acyl ethanolamines", "PPAR-alpha", "PPAR-gamma", "TRPV1 channel",
    
    # --- 4. LONGEVIDADE, SENESCÊNCIA & MORTE CELULAR ---
    "Senolytics", "Senomorphics", "Klotho", "Sirtuins", "SIRT1", "SIRT3", "SIRT6",
    "Nicotinamide mononucleotide (NMN)", "Nicotinamide riboside", "NAD+ metabolism",
    "Ferroptosis", "Pyroptosis", "Necroptosis", "Cuproptosis", "Autophagy", "Mitophagy",
    "mTOR pathway", "AMPK signaling", "Telomerase", "p16INK4a", "p21CIP1",
    
    # --- 5. RECEPTORES OLFATIVOS & GUSTATIVOS ECTÓPICOS ---
    "Olfactory receptor 51E2 (OR51E2)", "Olfactory receptor 2AT4 (OR2AT4)", "OR10J5", "OR1D2",
    "Bitter taste receptors (TAS2Rs)", "Sweet taste receptors (TAS1R2/TAS1R3)", 
    "Umami receptor", "Chemosensory receptors", "Ectopic olfactory receptors",
    
    # --- 6. MECANOTRANSDUÇÃO & CANAIS IÔNICOS ---
    "Piezo1", "Piezo2", "TMEM16A", "TMEM63", "TREK-1 channel", "TRAAK channel",
    "TRPV4", "TRPM8", "TRPA1", "P2X3 receptor", "P2X7 receptor", "ASIC channels",
    "Mechanosensitive channels", "Stretch-activated channels", "YAP/TAZ pathway", "Hippo pathway",
    "HCN channels", "Kv7 channels", "BK channels", "SK channels",
    
    # --- 7. COMUNICAÇÃO CELULAR & EXOSSOMAS ---
    "Exosomes", "Extracellular vesicles", "Microvesicles", "Exosomal miRNA", "Gap junctions", 
    "Connexin 43", "Pannexin 1", "Tunneling nanotubes",
    
    # --- 8. RECEPTORES ÓRFÃOS & GPCRs ---
    "GPR40 (FFAR1)", "GPR41 (FFAR3)", "GPR43 (FFAR2)", "GPR84", "GPR120 (FFAR4)", 
    "GPR35", "GPR183", "GPR17", "GPR30 (GPER)", "LGR5", "LGR6",
    "Muscarinic M3", "Beta-3 adrenergic receptor", "Nicotinic alpha-7", "Purinergic signaling",
    
    # --- 9. INFLAMAÇÃO & RESOLUÇÃO ---
    "NLRP3 inflammasome", "cGAS-STING pathway", "NF-kappaB", "HMGB1", 
    "Resolvin D1", "Resolvin D2", "Resolvin E1", "Lipoxin A4", "Maresins", "Protectins",
    "Specialized pro-resolving mediators (SPMs)", "Prostaglandin E2", "Cyclooxygenase-2",
    "Nitric Oxide Synthase", "Heme Oxygenase-1", "NRF2 pathway",
    
    # --- 10. METABÓLITOS CLÁSSICOS & SINALIZADORES ---
    "Kynurenic acid", "Succinate", "Lactate", "Fumarate", "Itaconate", "Alpha-lipoic acid",
    "Sphingosine-1-phosphate", "Ceramides", "Hyaluronic acid", "Taurine", "Uric Acid",
    "Melatonin", "Adenosine", "ATP", "Glutamate",
    
    # --- 11. GENÉTICA NÃO-CODIFICANTE ---
    "MALAT1 lncRNA", "HOTAIR lncRNA", "H19 lncRNA", "NEAT1 lncRNA",
    "miR-21", "miR-145", "miR-29", "miR-126", "Circular RNA (circRNA)", "Piwi-interacting RNA"
]

# --- FILTRO DE RUÍDO AUTOMATIZADO (BLACKLIST CATEGÓRICA) ---
# Se o termo contiver QUALQUER uma dessas palavras, ele é eliminado.
BLACKLIST_GERAL = [
    # 1. TEMPO & CALENDÁRIO (Remove "Apr", "Aug", "Summer 2024")
    "january", "february", "march", "april", "may", "june", "july", "august", "september", "october", "november", "december",
    "jan", "feb", "mar", "apr", "jun", "jul", "aug", "sep", "oct", "nov", "dec",
    "monday", "tuesday", "wednesday", "thursday", "friday", "saturday", "sunday",
    "spring", "summer", "autumn", "winter", "fiscal year", "annual",
    
    # 2. CORPORATIVO & INSTITUCIONAL (Remove "Astellas", "Merck", "University of...")
    "inc.", "ltd", "corp", "corporation", "company", "plc", "gmbh", "s.a.", "co.",
    "pharma", "pharmaceuticals", "biosciences", "biotechnologies", "therapeutics", "solutions",
    "biopharma", "medtech", "healthcare", "labs", "laboratories", "systems", "group",
    "university", "universidad", "universite", "college", "school", "academy", "institute",
    "department", "division", "center", "centre", "hospital", "clinic", "foundation",
    "association", "society", "board", "committee", "council", "federation", "assembly",
    "ministry", "agency", "authority", "commission", "consortium", "workshop",
    # Empresas específicas detectadas no seu log
    "astellas", "astrazeneca", "abbvie", "merck", "amgen", "allergan", "aveo", "accord", 
    "anchiano", "apogepha", "axonics", "atlab", "abyost", "aiq solutions", "astrin",
    
    # 3. PUBLICAÇÃO & METADADOS (Remove "Author", "PMID", "Review")
    "author", "editor", "publisher", "copyright", "rights reserved", "license",
    "pmid", "doi", "issn", "isbn", "vol.", "no.", "pp.", "page", "issue",
    "abstract", "introduction", "results", "discussion", "conclusion", "references",
    "acknowledgements", "supplementary", "figure", "table", "fig.", "tab.",
    "review", "article", "report", "letter", "editorial", "commentary", "survey",
    "received", "accepted", "published", "available online", "corresponding",
    "affiliated", "addressing", "approval no", "areas covered",
    
    # 4. GEOGRAFIA & DEMOGRAFIA (Remove "Americans", "Abu Dhabi", "Adults")
    "usa", "uk", "united states", "kingdom", "china", "japan", "germany", "france", "italy",
    "spain", "brazil", "russia", "india", "korea", "australia", "canada", "europe", "asia", "africa",
    "american", "european", "asian", "african", "japanese", "chinese", "brazilian",
    "abu dhabi", "los angeles", "new york", "london", "tokyo", "beijing", "seoul",
    "adults", "children", "adolescent", "infant", "elderly", "aged", "aging", "men", "women",
    "male", "female", "pediatric", "geriatric", "demographic", "population",
    
    # 5. METODOLOGIA GENÉRICA (Remove "Animal", "Study", "Adverse")
    "study", "trial", "experiment", "investigation", "analysis", "assessment", "evaluation",
    "method", "technique", "approach", "procedure", "protocol", "design",
    "animal", "human", "mouse", "mice", "rat", "murine", "rabbit", "dog", "monkey",
    "in vivo", "in vitro", "ex vivo", "cell culture", "clinical", "preclinical",
    "control", "placebo", "double-blind", "randomized", "cohort", "cross-sectional",
    "data", "statistics", "p-value", "significance", "confidence interval", "odds ratio",
    "adverse", "safety", "efficacy", "toxicity", "side effect", "tolerability",
    "treatment", "therapy", "administration", "dosage", "dose", "regimen",
    "pathophysiology", "mechanism", "pathway", "expression", "activation", "inhibition",
    "regulation", "modulation", "production", "synthesis", "secretion", "release",
    "level", "concentration", "activity", "function", "role", "effect", "impact",
    "associated with", "related to", "observed in", "induced by", "mediated by",
    "signaling", "transduction", "biomarker", "target", "potential", "candidate",
    
    # 6. PALAVRAS DE LIGAÇÃO & LIXO (Remove "The", "With", "Although")
    "the", "and", "of", "in", "to", "a", "an", "for", "with", "on", "at", "by", "from",
    "is", "are", "was", "were", "be", "been", "being", "have", "has", "had",
    "this", "that", "these", "those", "it", "its", "which", "who", "whom", "whose",
    "between", "among", "within", "during", "after", "before", "while", "when",
    "however", "although", "though", "nevertheless", "furthermore", "moreover",
    "additionally", "interestingly", "notably", "importantly", "specifically",
    "accordingly", "consequently", "therefore", "thus", "hence", "so",
    "advancing", "application", "applied", "aggregates", "anatomy", "angiogenesis"
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
        
        "footer_citar": "Lemos Lambda v1.3.0 - Uso Acadêmico",
        "citar_titulo": "📄 Como Citar",
        "citar_texto": "Lemos, G. (2025). Lemos Lambda: Deep Science Prospector [Software]. Versão 1.3.0. DOI: 10.5281/zenodo.17958507",
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
        
        "footer_citar": "Lemos Lambda v1.3.0 - Academic Use",
        "citar_titulo": "📄 How to Cite",
        "citar_texto": "Lemos, G. (2025). Lemos Lambda: Deep Science Prospector [Software]. Version 1.3.0. DOI: 10.5281/zenodo.17958507",
        "link_doi": "🔗 View on Zenodo (DOI)"
    }
}