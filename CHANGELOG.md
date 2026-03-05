# Changelog

All notable changes to the "Lemos Lambda" project will be documented in this file.

---

## [2.1.0] — 2026

### 🚀 New Features
- **Inverted Mining Logic:** System now searches 16 distant biological fields (proteostasis, microbiome, longevity, neuroprotection, gasotransmitters, resolution of inflammation, ferroptosis, mechanosensing, m6A epigenetics, etc.) without mentioning the target organ — inspired by the discovery of Trehalose and TMAO as candidates through cross-field literature reasoning.
- **Flash Bridge (mechanism → organ):** After collecting candidates from distant fields, Gemini Flash evaluates biological plausibility for the target organ. Reasoning example: *"Trehalose induces autophagy and reduces ER stress — urothelium under distension experiences ER stress — worth testing."*
- **16 Lateral Search Fields:** Defined in `CAMPOS_LATERAIS` in the backend, covering osmolytes, specialized pro-resolving mediators (SPMs), gasotransmitters, senolytics, m6A modulators, ferroptosis, mitochondrial metabolism, and approved drug repurposing.
- **Article Side Panel:** Clicking any row in the results table opens a lateral panel with related PubMed articles.
- **"View Abstract" Button:** Expands the article abstract without consuming any AI tokens.
- **"Investigate with AI" Button (lazy):** Flash reads the article and returns `[Tissue] → [Mechanism] → [Effect]`, called ONLY when clicked.
- **"Investigate Target" Button (Pro, lazy):** Gemini Pro performs deep analysis including mechanism, druggability, novelty score, suggested experiments, biological rationale, and red flags — called ONLY when clicked.
- **Real-time `st.status()` Progress:** Each mining step displays a live log (fields searched, candidates extracted, HPA check, virginity check).
- **Expanded "Type" Column:** Table now shows Receptor / Drug / Modulator / Gene / Ion Channel / Enzyme / Pathway / Biomarker.
- **Human Protein Atlas (HPA):** Integrated to verify target expression in the specific tissue of the target organ.
- **ClinicalTrials.gov Integration:** Fetches compounds in active trials for other indications (repurposing candidates).
- **OpenTargets Integration:** Fetches genetically validated targets.
- **Orbitron Font Logo:** New sci-fi logo with green glow effect.
- **PubMed Syntax Tooltip:** Help icon on the Main Target field explains OR, NOT, AND operators with examples.
- **Expanded Presets:** New categories — Autophagy & Metabolism, Vasculature & Fibrosis, Microbiome & Metabolites.
- **`biologicalRationale` Field:** Investigation JSON now includes an explanation of why the target makes sense in the organ even if never tested there.
- **Wide Funnel:** Returns up to 100 candidates instead of 16–20, letting the scientist decide.

### 🐛 Bug Fixes
- **Article card overflow:** Fixed article titles breaking outside the card container (`overflow:hidden`, `word-wrap:break-word`).
- **4-minute frozen screen:** Mining now shows real-time progress via `st.status()` — no more blank waiting.
- **Version label:** Corrected from v3.2 to v2.1.

### 🔧 Improvements
- **Flash prompt completely rewritten:** Perspective of a senior pharmacologist hunting for unexplored opportunities, with explicit instructions for cross-organ translational potential, animal-model-testable compounds, and "when in doubt — KEEP IT" criterion.
- **Zero tokens during statistical analysis:** No AI calls during Lambda Score calculation — Flash only called once for classification and on-demand per article.
- **Expanded blacklist:** Added academic text transition words (`multivariate`, `ameliorate`, `coordinate`, `alleviate`) and genomic noise patterns (lncRNAs, pseudogenes, antisense RNAs) via `LIXO_REGEX`.

---

## [2.0.0] — 2025-12-23

### 🚀 New Features
- **Statistical Engine:** Implemented Fisher's Exact Test to calculate the *Lambda Score*, providing statistical validation for target enrichment (p-value).
- **Generative Pharmacology (AI):** Integrated Google Gemini models to analyze abstracts. The AI now acts as a "Senior Researcher," extracting *Target → Drug → Physiological Effect* in natural language.
- **Universal API Compatibility:** Developed a "Cascade Fallback System" that automatically tests multiple Gemini versions (2.5, 2.0, 1.5, 1.0) to ensure compatibility with both legacy and new API keys.
- **Sherlock Holmes Mode (v2.10):** A new heuristic logic that infers pharmacological targets from the article Title when the Abstract is truncated or missing.
- **Scientific Radar:** New dashboard-style UI displaying the latest relevant news and findings with a modern card layout.
- **Editable Target List:** Users can now manually edit, copy, and paste the generated target list directly in the UI before execution.
- **3-File Architecture:** Codebase split into `app.py`, `backend.py`, and `constantes.py` for maintainability.
- **Real-time Progress Bar:** Statistical analysis updates a progress bar per term via Python generator.
- **Tag Classification:** Blue Ocean, Embryonic, Gold, Trending, Ghost, Saturated — with explicit thresholds.
- **CSV Export:** Download full results table with all metrics.
- **Temporal Trend:** Lambda Score now includes a recency ratio (2023–2026 articles).

### 🐛 Bug Fixes
- **Abstract Parsing:** Fixed a critical bug where the parser read only the first line of long PubMed abstracts; it now correctly handles multi-line indentation.
- **API Error Handling:** Resolved "Error 404" and "Error 400" by implementing better exception handling and user guidance for Google Cloud project creation.
- **Language Switching:** Fixed a crash that occurred when switching between English/Portuguese without clearing the previous results table.

### 🔧 Improvements
- **Prompt Engineering:** Refined the AI prompt (v2.9) to prevent "hallucinations" (e.g., writing emails) and focus strictly on data extraction.
- **Performance:** Optimized the Biopython search query to reduce latency during deep mining.

---

## [1.x] — Earlier versions (original `app_doutorado.py`)

- Single-file Python + Streamlit
- Basic PubMed mining with regex
- Simple Gemini classification
- Minimal UI without side panels
- Initial blacklist of generic terms

---

## Roadmap

- [ ] FAERS (FDA Adverse Events) — drugs with relevant side effects as mechanism signals
- [ ] Google Patents API — patented compounds not yet published for the target organ
- [ ] Automatic comparison with diseases of similar pathophysiology
- [ ] Cross-organ translatability score (kidney → bladder, intestine → bladder)
- [ ] Results caching to avoid re-querying the same terms
