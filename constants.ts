import { Language } from './types';

export const TEXTS = {
  [Language.PT]: {
    // --- GERAL ---
    titulo_desk: "λ Lemos Lambda: Deep Science Prospector",
    subtitulo: "Ferramenta de Prospecção Farmacológica Dinâmica (v2.1)",
    radar_titulo: "📡 Radar Científico (Últimas 24h)",
    btn_ler: "Ler Artigo",
    
    // --- INPUTS & CONFIG ---
    header_scope: "1. Definição do Escopo",
    label_email: "E-mail do Pesquisador (Obrigatório NCBI)",
    holder_email: "ex: guilherme@unifesp.br",
    erro_email: "⚠️ E-mail é obrigatório para usar a API!",
    label_alvo: "Alvo Principal (Órgão/Doença)",
    holder_alvo: "ex: Bladder, Urothelium, Detrusor, Heart",
    helper_alvo: "* Dica: Você pode inserir múltiplos órgãos ou doenças separados por vírgula (ex: Pulmão, Fibrose).",
    erro_alvo: "⚠️ Defina o Alvo Principal primeiro!",
    btn_auto: "🧠 AUTO-DETECTAR ALVOS & INICIAR",
    expander_presets: "💎 Adicionar Fronteiras (Blue Ocean) & Presets",
    label_categoria: "Categoria:",
    btn_add_preset: "📥 Adicionar Categoria",
    btn_add_all: "💎 ADICIONAR TODOS OS PRESETS (Varredura Total)",
    popover_manual: "✍️ Adição Manual",
    holder_manual: "Ex: ATP, P2X3",
    btn_add_manual: "Adicionar Manual",
    header_config: "Configuração",
    slider_tempo: "Janela de Tempo",
    label_contexto: "Contexto (Opcional)",
    upload_label: "Importar Lista (CSV/TXT)",
    
    // --- ESTRATÉGIAS IA ---
    lbl_ai_strategy: "Estratégia de Mineração IA",
    helper_ai_strategy: "Selecione e clique em Minerar para ADICIONAR à lista",
    strat_clean: "Faxina",
    strat_repurpose: "Reposicionar",
    strat_mechanism: "Mecanismo",
    strat_blue_ocean: "Blue Ocean",
    
    // TOOLTIPS ESTRATÉGIAS
    desc_clean: "Remove termos genéricos (ex: 'rato', 'estudo') e padroniza nomes (ex: 'ATP receptor' -> 'P2X3'). Não adiciona novos alvos.",
    desc_repurpose: "Busca fármacos já aprovados ou em fase clínica que atuam neste alvo/doença (Drug Repurposing).",
    desc_mechanism: "Foca em vias de sinalização (Quinases, Fatores de Transcrição) e mecanismos moleculares upstream/downstream.",
    desc_blue_ocean: "Modo Discovery: Busca alvos inéditos, pouco citados ou descobertos nos últimos 5 anos. Ignora o óbvio.",

    btn_run_prefix: "🤖 RODAR:",
    
    // --- CONTEXTO ---
    lbl_contexto_tit: "Comparador Biológico (Contexto)",
    lbl_contexto_help: "* A IA usa isso para buscar reposicionamento ou vias homólogas.",
    
    // --- AVISOS ---
    warn_no_key_tit: "Modo Limitado (Sem IA)",
    warn_no_key_desc: "Sem API Key, a mineração é apenas por texto. Para Blue Ocean, configure a chave.",
    upload_helper: "* Aceita termos exatos (vírgula/linha).",
    
    // --- LISTA ---
    lbl_targets_preview: "LISTA DE ALVOS",
    lbl_targets_count: "ALVOS",
    
    // --- CONFIG IA ---
    expander_ia: "🧠 Ativar Inteligência Artificial (Gemini)",
    caption_ia: "Insira sua chave gratuita do Google para resumos inteligentes.",
    placeholder_key: "Cole sua chave AIza... aqui",
    link_key: "Obter chave grátis",
    
    // --- FEEDBACK ---
    status_minerando: "Minerando literatura para",
    toast_atualizado: "Lista atualizada!",
    sucesso_carregado: "alvos carregados. Clique em EXECUTAR abaixo.",
    msg_alvos_ok: "alvos prontos para análise.",
    expander_lista: "👀 Ver/Editar Lista Completa",
    btn_limpar: "Limpar Lista",
    btn_executar: "🚀 EXECUTAR ANÁLISE ESTATÍSTICA",
    toast_importado: "Importação concluída",
    erro_arquivo: "Erro ao ler arquivo.",
    spinner_analise: "Calibrando estatística (Lemos Lambda v2.1)...",
    titulo_processando: "## 🧬 Processando Estatística...",

    // --- RESULTADOS ---
    resultados: "Relatório de Inteligência",
    btn_voltar: "⬅ Voltar",
    metric_top: "Top Oportunidade",
    metric_score: "Lambda Score",
    metric_hits: "Artigos Locais",
    header_heatmap: "Mapa de Calor",
    btn_baixar: "📥 Baixar Relatório (CSV)",
    header_leitura: "Leitura Guiada",
    label_investigar: "Escolha o alvo:",
    btn_ver_papers: "📄 Ver Abstracts",
    btn_investigar: "🤖 Analisar com IA",
    spinner_investigando: "Lendo artigos com Inteligência Artificial...",
    btn_pubmed: "Ver no PubMed",

    // --- COLUNAS ---
    col_mol: "Alvo Molecular",
    col_status: "Classificação",
    col_ratio: "Lambda Score",
    col_pvalue: "P-Valor",
    col_art_alvo: "Hits (Alvo)",
    col_global: "Hits (Global)",

    // --- TAGS ---
    tag_blue_ocean: "💎 Blue Ocean",
    tag_ghost: "👻 Fantasma",
    tag_embryonic: "🌱 Embrionário",
    tag_neutral: "⚖️ Neutro",
    tag_gold: "🥇 Ouro",
    tag_trending: "🚀 Tendência",
    tag_saturated: "🔴 Saturado",

    // --- RODAPÉ ---
    footer_citar: "Uso Acadêmico",
    citar_titulo: "📄 Como Citar",
    citar_texto: "Lemos, G. (2025). Lemos Lambda: Deep Science Prospector (v2.1.0). Zenodo. https://doi.org/10.5281/zenodo.18090934",
    apoio_titulo: "Apoie o Pesquisador (Pix):",
    apoio_desc: "Ajude a manter o código atualizado e os servidores ativos ☕",
    pix_key: "960f3f16-06ce-4e71-9b5f-6915b2a10b5a",
    btn_copy: "Copiar"
  },
  [Language.EN]: {
    // --- GENERAL ---
    titulo_desk: "λ Lemos Lambda: Deep Science Prospector",
    subtitulo: "Dynamic Pharmacological Prospecting Tool (v2.1)",
    radar_titulo: "📡 Scientific Radar (Last 24h)",
    btn_ler: "Read Article",
    
    // --- INPUTS ---
    header_scope: "1. Scope Definition",
    label_email: "Researcher E-mail (NCBI Required)",
    holder_email: "ex: researcher@university.edu",
    erro_email: "⚠️ E-mail is required for API usage!",
    label_alvo: "Main Target (Organ/Disease)",
    holder_alvo: "ex: Bladder, Urothelium, Detrusor",
    helper_alvo: "* Tip: You can enter multiple organs or diseases separated by commas (e.g., Lung, Asthma).",
    erro_alvo: "⚠️ Define Target first!",
    btn_auto: "🧠 AUTO-DETECT TARGETS & START",
    expander_presets: "💎 Blue Ocean Frontiers & Presets",
    label_categoria: "Category:",
    btn_add_preset: "📥 Add Category",
    btn_add_all: "💎 ADD ALL PRESETS (Full Scan)",
    popover_manual: "✍️ Manual Input",
    holder_manual: "Ex: ATP, P2X3",
    btn_add_manual: "Add Manual",
    header_config: "Configuration",
    slider_tempo: "Time Range",
    label_contexto: "Context (Optional)",
    upload_label: "Import List (CSV/TXT)",
    
    // --- AI STRATEGIES ---
    lbl_ai_strategy: "AI Mining Strategy",
    helper_ai_strategy: "Select tool & Click Mine to ADD to results",
    strat_clean: "Clean",
    strat_repurpose: "Repurpose",
    strat_mechanism: "Mechanism",
    strat_blue_ocean: "Blue Ocean",

    // TOOLTIPS STRATEGIES
    desc_clean: "Removes generic terms (e.g., 'study', 'rat') and standardizes names (e.g., 'ATP receptor' -> 'P2X3'). No new targets.",
    desc_repurpose: "Finds existing approved drugs or clinical candidates for this target/disease.",
    desc_mechanism: "Focuses on signaling pathways (Kinases, Transcription Factors) and upstream/downstream mechanisms.",
    desc_blue_ocean: "Discovery Mode: Finds novel, understudied targets discovered in the last 5 years. Ignores the obvious.",

    btn_run_prefix: "🤖 RUN:",

    // --- CONTEXT ---
    lbl_contexto_tit: "Biological Comparator (Context)",
    lbl_contexto_help: "* AI uses this to find repurposing or homologous pathways.",
    
    // --- WARNINGS ---
    warn_no_key_tit: "Limited Mode (No AI)",
    warn_no_key_desc: "Without API Key, mining is text-only. For Blue Ocean, set the key.",
    upload_helper: "* Accepts exact terms (comma/newline).",

    // --- LIST ---
    lbl_targets_preview: "TARGET LIST",
    lbl_targets_count: "TARGETS",

    // --- AI CONFIG ---
    expander_ia: "🧠 Activate Artificial Intelligence (Gemini)",
    caption_ia: "Insert your free Google API Key for smart summaries.",
    placeholder_key: "Paste your AIza... key here",
    link_key: "Get free key",
    
    // --- FEEDBACK ---
    status_minerando: "Mining literature for",
    toast_atualizado: "List updated!",
    sucesso_carregado: "targets loaded. Click EXECUTE below.",
    msg_alvos_ok: "targets ready for analysis.",
    expander_lista: "👀 View/Edit Full List",
    btn_limpar: "Clear List",
    btn_executar: "🚀 EXECUTE STATISTICAL ANALYSIS",
    toast_importado: "Import successful",
    erro_arquivo: "Error reading file.",
    spinner_analise: "Calibrating statistics (Lemos Lambda v2.1)...",
    titulo_processando: "## 🧬 Processing Statistics...",

    // --- RESULTS ---
    resultados: "Intelligence Report",
    btn_voltar: "⬅ Back",
    metric_top: "Top Opportunity",
    metric_score: "Lambda Score",
    metric_hits: "Local Hits",
    header_heatmap: "Heatmap",
    btn_baixar: "📥 Download Report (CSV)",
    header_leitura: "Guided Reading",
    label_investigar: "Choose target:",
    btn_ver_papers: "📄 View Abstracts",
    btn_investigar: "🤖 Analyze with AI",
    spinner_investigando: "Reading articles with AI...",
    btn_pubmed: "View on PubMed",

    // --- COLUMNS ---
    col_mol: "Molecular Target",
    col_status: "Classification",
    col_ratio: "Lambda Score",
    col_pvalue: "P-Value",
    col_art_alvo: "Hits (Target)",
    col_global: "Hits (Global)",

    // --- TAGS ---
    tag_blue_ocean: "💎 Blue Ocean",
    tag_ghost: "👻 Ghost",
    tag_embryonic: "🌱 Embryonic",
    tag_neutral: "⚖️ Neutral",
    tag_gold: "🥇 Gold",
    tag_trending: "🚀 Trending",
    tag_saturated: "🔴 Saturated",

    // --- FOOTER ---
    footer_citar: "Academic Use",
    citar_titulo: "📄 How to Cite",
    citar_texto: "Lemos, G. (2025). Lemos Lambda: Deep Science Prospector (v2.1.0). Zenodo. https://doi.org/10.5281/zenodo.18090934",
    apoio_titulo: "Support the Researcher (Pix):",
    apoio_desc: "Help keep the code updated and servers running ☕",
    pix_key: "960f3f16-06ce-4e71-9b5f-6915b2a10b5a",
    btn_copy: "Copy"
  }
};

export const PRESETS_FRONTEIRA: Record<string, string[]> = {
  "🔥 Blockbusters & Inovação (2025)": [
    "Semaglutide", "Tirzepatide", "Retatrutide", "Empagliflozin", "Dapagliflozin",
    "Iptacopan", "Lecanemab", "Donanemab", "Zuranolone", "Resmetirom",
    "TFEB", "TMAO", "PROTACs", "Molecular Glues", "Senolytics", "Finerenone", "Vericiguat"
  ],
  "❤️ Cardio-Renal & Metabólico": [
    "SGLT2", "SGLT1", "GLP1R", "GIPR", "GCGR", "PCSK9", "NLRP3", 
    "Neprilysin", "sGC", "HCN4", "RyR2", "SERCA2a", "Aldosterone Synthase", 
    "Lp(a)", "ANGPTL3", "Factor XI", "CETP", "Adiponectin", "Leptin"
  ],
  "⚡ Mecanossensores & Canais (Core)": [
    "Piezo1", "Piezo2", "TRPV1", "TRPV4", "TRPA1", "TRPM8", "P2X3", "P2X7", 
    "Nav1.7", "Nav1.8", "Kv7.2/7.3", "BK Channel", "SK Channel", "TMEM16A", "OSCA1", "TMEM63", 
    "HCN1", "TREK-1", "TREK-2", "TRAAK"
  ],
  "🌿 Canabinoides & Lipídios": [
    "Cannabidiol (CBD)", "Tetrahydrocannabinol (THC)", "CB1 Receptor", "CB2 Receptor", 
    "GPR55", "GPR18", "GPR119", "FAAH", "MAGL", "Anandamide (AEA)", "2-AG", 
    "Resolvin D1", "Resolvin E1", "Maresin-1", "Lipoxin A4", "Sphingosine-1-Phosphate"
  ],
  "💨 Gasotransmissores & Doadores": [
    "Carbon Monoxide (CO)", "HO-1", "CORM-2", "CORM-401", 
    "Hydrogen Sulfide (H2S)", "CBS", "CSE", "GYY4137", "AP39", 
    "Hydrogen Gas (H2)", "Sulfur Dioxide (SO2)", "Nitric Oxide (NO)", "sGC"
  ],
  "🧬 Epigenética, RNAs & Vias": [
    "TET2", "HDAC6 inhibitor", "DNMT1", "MALAT1", "HOTAIR", "miR-29b", "miR-132", "miR-145", "Exosomes",
    "YAP", "TAZ", "Hippo Pathway", "mTORC1", "AMPK", "STING", "cGAS", "Ferroptosis", "Autophagy"
  ],
  "👻 Receptores Órfãos": [
    "GPR35", "GPR84", "GPR68", "GPR120", "OR51E2", "OR1D2", "TAS2R14", "TAS2R38"
  ]
};

export const MOCK_NEWS = [
  { id: 1, title: "New CRISPR technique allows editing of RNA in living cells", source: "Nature Biotechnology", link: "#", image: "https://picsum.photos/400/200?random=1" },
  { id: 2, title: "Senolytics show promise in reversing age-related macular degeneration", source: "Cell Metabolism", link: "#", image: "https://picsum.photos/400/200?random=2" },
  { id: 3, title: "Gut microbiome diversity linked to longevity in centenarians", source: "Science", link: "#", image: "https://picsum.photos/400/200?random=3" },
];

export const BLACKLIST_MONSTRO = new Set([
  // Images & Physics
  "MRI", "CT", "PET", "RBE", "IMPT", "VIII", "RATIO",
  // Methods
  "ASSOCIATION", "EVALUATION", "UNCOMMON", "DISEASE", "BACKGROUND", 
  "OBJECTIVE", "METHODS", "RESULTS", "CONCLUSION", "ABSTRACT", 
  "INTRODUCTION", "STUDY", "ANALYSIS", "DATA", "STATISTICS", 
  "SIGNIFICANT", "DIFFERENCE", "BETWEEN", "AMONG", "WITHIN", 
  "DURING", "PREVALENCE", "INCIDENCE", "RISK", "FACTOR", "ROLE", 
  "POTENTIAL", "NOVEL", "DIAGNOSTIC", "ARTIFICIAL", "MANAGEMENT",
  "PROGNOSTIC", "FACTORS", "NEUROMODULATION", "IMPLICATIONS",
  "CLINICAL", "REVIEW", "META-ANALYSIS", "SYSTEMATIC", "SURVEY",
  // Context
  "INTRAUTERINE", "GERMLINE", "BLADDER", "CANCER", "URINARY", 
  "UROTHELIAL", "MUSCLE", "OVERACTIVE", "TUMOR", "CARCINOMA", 
  "PELVIC", "URETHRAL", "UROLOGIC", "NEUROGENIC", "DIABETES", 
  "SJOGREN", "INTRAVESICAL", "VOID", "VOIDING", "DETRUSOR", 
  "PATIENT", "PATIENTS", "CHILDREN", "ADULT", "NHANES", "WOMEN", "MEN",
  // Biology Generic
  "DNA", "RNA", "ATP", "GENE", "PROTEIN", "CELL", "EXPRESSION",
  "BCG", "GUERIN", "CALMETTE", "BACILLUS", "HLA", "CAR",
  // Stop Words
  "THE", "AND", "FOR", "NOT", "BUT", "VIA", "ALL", "WITH", "FROM", "AFTER"
]);