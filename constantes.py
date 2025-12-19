# constantes.py
# Configurações, Textos e Listas de Prospecção
# Versão: 1.6.1

TEXTOS = {
    "pt": {
        "titulo_desk": "λ Lemos Lambda: Deep Science Prospector",
        "subtitulo": "Ferramenta de Prospecção Farmacológica Dinâmica",
        "radar_titulo": "📡 Radar Científico (Últimas 24h)",
        "label_email": "E-mail do Pesquisador (Obrigatório PubMed)",
        "holder_email": "ex:guilherme@unifesp.br",
        "label_alvo": "Alvo Principal (Órgão/Doença)",
        "holder_alvo": "ex: Bladder, Urothelium, Detrusor",
        "label_contexto": "Contexto Comparativo (Opcional)",
        "holder_fonte": "ex: Heart, Kidney, Liver (Órgãos similares separados por vírgula)",
        # -----------------------
        "status_minerando": "Minerando literatura para",
        "resultados": "Relatório de Inteligência",
        "btn_baixar": "📥 Baixar Relatório (CSV)",
        "col_mol": "Alvo Molecular",
        "col_status": "Classificação",
        "col_ratio": "Lambda Score",
        "col_art_alvo": "Hits (Alvo)",
        "col_global": "Hits (Global)",
        "footer_rights": "© 2025 Guilherme Lemos | Uso Acadêmico - Unifesp",
        "footer_citar": "Uso Acadêmico - Unifesp",
        "citar_titulo": "📄 Como Citar",
        "citar_texto": "Lemos, G. (2025). Lemos Lambda: Deep Science Prospector (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.17958507",
        "link_doi": "🔗 Ver no Zenodo (DOI)",
        "apoio_titulo": "Gostou da ferramenta?",
        "apoio_btn": "Buy Me a Coffee"
    },
    "en": {
        "titulo_desk": "λ Lemos Lambda: Deep Science Prospector",
        "subtitulo": "Dynamic Pharmacological Prospecting Tool",
        "radar_titulo": "📡 Scientific Radar (Last 24h)",
        "label_email": "Researcher E-mail (NCBI Required)",
        "holder_email": "ex: researcher@university.edu",
        "label_alvo": "Main Target (Organ/Disease)",
        "holder_alvo": "ex: Bladder, Urothelium, Detrusor",
        "label_contexto": "Context (Optional)",
        "holder_fonte": "ex: Heart, Kidney (Similar organs comma separated)",
        "status_minerando": "Mining literature for",
        "resultados": "Intelligence Report",
        "btn_baixar": "📥 Download Report (CSV)",
        "col_mol": "Molecular Target",
        "col_status": "Classification",
        "col_ratio": "Lambda Score",
        "col_art_alvo": "Hits (Local)",
        "col_global": "Hits (Global)",
        "footer_rights": "© 2025 Guilherme Lemos | Academic Use - Unifesp",
        "footer_citar": "Academic Use",
        "citar_titulo": "📄 How to Cite",
        "citar_texto": "Lemos, G. (2025). Lemos Lambda: Deep Science Prospector (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.17958507",
        "link_doi": "🔗 View on Zenodo (DOI)",
        "apoio_titulo": "Enjoying the tool?",
        "apoio_btn": "Buy Me a Coffee"
    }
}

# --- LISTAS DE PROSPECÇÃO ---

PRESETS_FRONTEIRA = {
    "🔥 Gasotransmissores": [
        "Carbon Monoxide (CO)", "HO-1", "CORM-2", "CORM-401", "Hydrogen Sulfide (H2S)", 
        "CBS", "CSE", "GYY4137", "AP39", "Hydrogen Gas (H2)", "Sulfur Dioxide (SO2)"
    ],
    "⚡ Mecanossensores & Canais": [
        "Piezo1", "Piezo2", "OSCA1", "TMEM63", "TRPA1", "TRPM8", "TRPV4", 
        "TREK-1", "TREK-2", "TRAAK", "TMEM16A", "HCN1", "HCN4", "P2X3", "P2X7"
    ],
    "🧬 Epigenética & RNAs": [
        "TET2", "HDAC6 inhibitor", "DNMT1", "MALAT1", "HOTAIR", 
        "miR-29b", "miR-132", "miR-145", "Exosomes"
    ],
    "🥑 Resolução & Lipídios": [
        "Resolvin D1", "Resolvin E1", "Maresin-1", "Lipoxin A4", 
        "Anandamide (AEA)", "2-AG", "GPR55", "GPR18", "GPR120"
    ],
    "👻 Receptores Órfãos": [
        "GPR35", "GPR84", "GPR68", "OR51E2", "OR1D2", "TAS2R14", "TAS2R38"
    ]
}

CANDIDATOS_MINERACAO = (
    PRESETS_FRONTEIRA["🔥 Gasotransmissores"] + 
    PRESETS_FRONTEIRA["⚡ Mecanossensores & Canais"] +
    PRESETS_FRONTEIRA["👻 Receptores Órfãos"]
)

BLACKLIST_GERAL = ["review", "meta-analysis", "rat", "mice", "human", "study"]