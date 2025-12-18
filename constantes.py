# constantes.py
# Armazena textos, listas estáticas, traduções e PRESETS.

# --- PRESET LEMOS (DOUTORADO) ---
PRESET_LEMOS = {
    "alvo": "Bladder OR Vesical OR Urothelium OR Detrusor OR Cystitis OR Overactive Bladder",
    "fonte": "Brain OR Kidney OR Liver OR Intestine OR Lung OR Vascular OR Immune System"
}

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

TEXTOS = {
    "pt": {
        "titulo_desk": "λ Lemos Lambda: Deep Science",
        "subtitulo": "Ferramenta de Prospecção Farmacológica",
        "step_1": "1️⃣ Defina seu Alvo",
        "step_2": "2️⃣ Mineração Profunda",
        "label_email": "E-mail (Obrigatório para PubMed):",
        "holder_email": "ex: pesquisador@unifesp.br",
        "label_alvo": "Qual órgão ou doença?",
        "holder_alvo": "ex: Overactive Bladder, Fibrosis...",
        "aviso_pubmed": "⚠️ **Atenção:** Para o PubMed funcionar, escreva os termos em **INGLÊS** (ex: *Kidney* em vez de Rim).",
        "btn_magic": "✨ Descobrir 'Blue Oceans' (Automático)",
        "prog_magic": "A IA está varrendo a literatura recente...",
        "status_minerando": "Lendo abstracts recentes...",
        "status_filtrando": "Identificando moléculas complexas...",
        "status_pronto": "Bibliotecas atualizadas!",
        "analise_btn": "🚀 Calcular Potencial (Ratio)",
        "resultados": "🎯 Resultados da Prospecção",
        "tabela_vazia": "Adicione termos ou use a descoberta automática acima.",
        "footer_citar": "Lemos Lambda v1.4 - Uso Acadêmico",
        
        # Novos campos v1.4
        "btn_preset": "🎓 Carregar Preset Doutorado (Lemos)",
        "toast_preset": "🧬 Preset Lemos carregado com sucesso!",
        "label_periodo": "📅 Período de Análise",
        "label_manual": "🔎 Investigar Termo Específico",
        "holder_manual": "ex: Curcumina, Gene X...",
        "btn_add_manual": "➕ Adicionar à Lista",
        "toast_add": "✅ termo(s) adicionado(s)!",
        "toast_dup": "⚠️ Duplicatas ignoradas.",
        
        "label_fonte": "Filtro de Fonte (Tecido/Célula)",
        "holder_fonte": "ex: Urothelium, Smooth Muscle...",
        "desc_fonte": "Opcional: Restringir comparação a um tecido específico.",
        "titulo_import": "📂 Importar Lista",
        "desc_import": "Upload (.csv/.txt)",
        "toast_import": "✅ termos importados!",
        "erro_ler": "Erro ao ler arquivo.",
        "btn_limpar": "🗑️ Limpar",
        "btn_limpar_tudo": "🗑️ Limpar Tudo",
        "ver_editar": "📝 Ver/Editar Lista de Palavras-Chave",
        "qtd_termos": "Qtd:",
        "radar_titulo": "📡 Radar Científico (Updates via RSS)",
        "btn_ler_feed": "Ler Completo",
        
        # Colunas e Métricas
        "metrica_potencial": "🏆 Maior Potencial",
        "metrica_score": "📊 Score (Ratio)",
        "metrica_artigos": "📚 Artigos (Alvo)",
        "col_mol": "Molécula/Alvo",
        "col_status": "Status",
        "col_ratio": "Potencial (Ratio)",
        "col_art_alvo": "Artigos no Alvo",
        "col_global": "Global/Fonte",
        "btn_baixar": "📥 Baixar Relatório CSV",
        "erro_email": "E-mail necessário.",
        "erro_campos": "⚠️ Preencha E-mail e Alvo (em Inglês)!",
        
        "citar_titulo": "📄 Como Citar",
        "citar_texto": "Lemos, G. (2025). Lemos Lambda: Deep Science Prospector [Software]. Versão 1.4.0. DOI: 10.5281/zenodo.17958507",
        "link_doi": "🔗 Ver no Zenodo (DOI)"
    },
    "en": {
        "titulo_desk": "λ Lemos Lambda: Deep Science",
        "subtitulo": "Pharmacological Prospecting Tool",
        "step_1": "1️⃣ Define Target",
        "step_2": "2️⃣ Deep Mining",
        "label_email": "E-mail (Required for PubMed):",
        "holder_email": "ex: researcher@university.edu",
        "label_alvo": "Target Organ or Disease?",
        "holder_alvo": "ex: Overactive Bladder, Fibrosis...",
        "aviso_pubmed": "⚠️ **Warning:** Please input terms in **ENGLISH** for PubMed accuracy (e.g., *Kidney* instead of Rim).",
        "btn_magic": "✨ Discover 'Blue Oceans' (Auto)",
        "prog_magic": "AI is scanning recent literature...",
        "status_minerando": "Reading recent abstracts...",
        "status_filtrando": "Identifying complex molecules...",
        "status_pronto": "Libraries updated!",
        "analise_btn": "🚀 Calculate Potential (Ratio)",
        "resultados": "🎯 Prospecting Results",
        "tabela_vazia": "Add terms or use automatic discovery above.",
        "footer_citar": "Lemos Lambda v1.4 - Academic Use",
        
        # New fields v1.4
        "btn_preset": "🎓 Load Lemos PhD Preset",
        "toast_preset": "🧬 Lemos Preset loaded successfully!",
        "label_periodo": "📅 Analysis Period",
        "label_manual": "🔎 Investigate Specific Term",
        "holder_manual": "ex: Curcumin, Gene X...",
        "btn_add_manual": "➕ Add to List",
        "toast_add": "✅ term(s) added!",
        "toast_dup": "⚠️ Duplicates ignored.",
        
        "label_fonte": "Source Filter (Tissue/Cell)",
        "holder_fonte": "ex: Urothelium, Smooth Muscle...",
        "desc_fonte": "Optional: Restrict comparison to specific tissue.",
        "titulo_import": "📂 Import List",
        "desc_import": "Upload (.csv/.txt)",
        "toast_import": "✅ terms imported!",
        "erro_ler": "Error reading file.",
        "btn_limpar": "🗑️ Clear",
        "btn_limpar_tudo": "🗑️ Clear All",
        "ver_editar": "📝 View/Edit Keywords List",
        "qtd_termos": "Qty:",
        "radar_titulo": "📡 Science Radar (RSS Updates)",
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
        
        "citar_titulo": "📄 How to Cite",
        "citar_texto": "Lemos, G. (2025). Lemos Lambda: Deep Science Prospector [Software]. Version 1.4.0. DOI: 10.5281/zenodo.17958507",
        "link_doi": "🔗 View on Zenodo (DOI)"
    }
}