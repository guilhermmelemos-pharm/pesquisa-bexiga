# constantes.py
# Armazena textos, listas estáticas, traduções e PRESETS.

# --- PRESET LEMOS (DOUTORADO) ---
# Seus parâmetros fixos de comparação
PRESET_LEMOS = {
    "alvo": "Bladder OR Vesical OR Urothelium OR Detrusor OR Cystitis OR Overactive Bladder",
    "fonte": "Brain OR Kidney OR Liver OR Intestine OR Lung OR Vascular OR Immune System"
}

# --- SUA LISTA ÚNICA (BASE) ---
CANDIDATOS_MINERACAO = [
    # --- Ácidos e Metabólitos ---
    "Alpha-lipoic acid", "Butyrate", "Short-chain fatty acids", "Sphingosine-1-phosphate",
    "Hyaluronic acid", "Succinate", "Lactate", "Kynurenic acid", "Prostaglandin E2",
    "Resolvin D1", "Lipoxin A4", "Melatonin", "Taurine", "Uric Acid",
    # --- Receptores e Canais ---
    "P2X3 receptor", "P2X7 receptor", "TRPV1 channel", "TRPM8", "Piezo1", "Piezo2",
    "Cannabinoid receptor 2", "GPR119", "GPR40", "GPR84", "GPR120",
    "Beta-3 adrenergic receptor", "Muscarinic M3", "Nicotinic alpha-7",
    # --- Vias e Enzimas ---
    "Rho-kinase (ROCK)", "mTOR pathway", "AMPK signaling", "NLRP3 inflammasome",
    "Nitric Oxide Synthase", "Heme Oxygenase-1", "Cyclooxygenase-2", "Phosphodiesterase-5",
    # --- Genética e RNA ---
    "MALAT1 lncRNA", "miR-21", "miR-145", "SIRT1", "NRF2 pathway", "NF-kappaB"
]

# --- FILTRO DE RUÍDO ---
BLACKLIST_GERAL = [
    "Adverse events", "Adverse effect", "AE rates", "Safety", "Efficacy", "Placebo",
    "Control group", "Study design", "Double-blind", "Randomized", "Clinical trial",
    "P-value", "Confidence interval", "Odds ratio", "Hazard ratio", "Standard deviation",
    "ANOVA", "Regression", "Analysis", "Data", "Results", "Conclusion", "Methods",
    "Significant", "Statistically", "Increased", "Decreased", "Compared to",
    "Signaling pathway", "Signal transduction", "Gene expression", "Protein levels",
    "Messenger RNA", "Receptor agonist", "Receptor antagonist", "Inhibitor",
    "Mechanism of action", "Therapeutic target", "Potential target", "Biomarker",
    "Pathophysiology", "Metabolism", "Oxidative stress", "Inflammation",
    "Cell culture", "In vivo", "In vitro", "Western blot", "PCR", "ELISA",
    "Stem cell", "Progenitor cell", "The", "And", "With", "For", "That", "This", 
    "Were", "Was", "Have", "Has", "Between", "Among", "During", "After", "Before", 
    "However", "Therefore", "Furthermore", "Moreover", "Additionally", "Notably", 
    "Interestingly", "Department", "University", "Hospital", "Institute", "Center", 
    "USA", "China", "Brazil", "Europe", "Funding", "Grant", "Review", "Article", "Copyright"
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
        
        # NOVOS BOTÕES INTELIGENTES
        "btn_smart_load": "🔄 Carregar Minha Lista (+ Dinâmica)",
        "desc_smart_load": "Carrega sua lista base e busca novidades automaticamente para o Alvo inserido.",
        "btn_preset": "🎓 Carregar Preset Doutorado (Bexiga)",
        
        "status_minerando": "🔍 Minerando novidades para:",
        "msg_sucesso_dinamico": "✅ Lista Base + {qtd} novidades específicas carregadas!",
        "msg_sucesso_base": "✅ Apenas Lista Base carregada (Preencha 'Alvo' para torná-la dinâmica).",
        
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
        "titulo_import": "📂 Importar Lista Extra",
        "desc_import": "Upload (.csv/.txt)",
        "toast_import": "✅ termos importados!",
        "erro_ler": "Erro ao ler arquivo.",
        "btn_limpar": "🗑️",
        "btn_limpar_tudo": "🗑️ Limpar Lista",
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
        "footer_citar": "Lemos Lambda v2.2 - Uso Acadêmico",
        "citar_titulo": "📄 Como Citar",
        "citar_texto": "Lemos, G. (2025). Lemos Lambda: Deep Science Prospector [Software]. Versão 2.2.0. DOI: 10.5281/zenodo.17958507",
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
        
        # NEW SMART BUTTONS
        "btn_smart_load": "🔄 Load My List (+ Dynamic)",
        "desc_smart_load": "Loads your base list and automatically mines novelties for the input Target.",
        "btn_preset": "🎓 Load Lemos PhD Preset",
        
        "status_minerando": "🔍 Mining novelties for:",
        "msg_sucesso_dinamico": "✅ Base List + {qtd} specific novelties loaded!",
        "msg_sucesso_base": "✅ Base List loaded only (Fill 'Target' to make it dynamic).",
        
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
        "titulo_import": "📂 Import Extra List",
        "desc_import": "Upload (.csv/.txt)",
        "toast_import": "✅ terms imported!",
        "erro_ler": "Error reading file.",
        "btn_limpar": "🗑️",
        "btn_limpar_tudo": "🗑️ Clear List",
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
        "footer_citar": "Lemos Lambda v2.2 - Academic Use",
        "citar_titulo": "📄 How to Cite",
        "citar_texto": "Lemos, G. (2025). Lemos Lambda: Deep Science Prospector [Software]. Version 2.2.0. DOI: 10.5281/zenodo.17958507",
        "link_doi": "🔗 View on Zenodo (DOI)"
    }
}