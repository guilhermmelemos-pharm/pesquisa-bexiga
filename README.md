# λ Lemos Lambda v2.1

> **Deep Science Prospector** — Finds unexplored therapeutic opportunities by crossing distant literature fields with the biology of the target organ.

---

## What is it

Lemos Lambda automates the lateral reasoning an experienced researcher performs when crossing literature from different fields. Instead of searching "what is already known about the bladder", the system asks: **"what is being discovered in microbiome, autophagy, gasotransmitters, and longevity research that could work in the bladder and nobody has tested yet?"**

The result are candidates like Trehalose and TMAO — compounds with solid mechanisms in other contexts, plausible biology in the target organ, and zero competition in the literature.

---

## How it works

### Inverted logic: mechanism → organ

```
16 distant fields (PubMed without mentioning the target organ)
    ↓
Raw candidates (~1000–1500)
    ↓  Regex filter (lncRNAs, pseudogenes, academic text noise)
Clean candidates
    ↓  Human Protein Atlas (expression in the tissue)
Expressed candidates
    ↓  PubMed virginity check (≤5 articles in organ = Blue Ocean)
Ranked by novelty
    ↓  Gemini Flash — "does this make biological sense in the organ?"
Final list 80–100 candidates
    ↓  Lambda Score statistical analysis
Ranked report
```

### Lateral search fields

| Field | Example candidates |
|-------|--------------------|
| Proteostasis & Chaperones | HSP70, GRP78, 4-PBA, Trehalose |
| Autophagy | ULK1, ATG5, Beclin1, Spermidine |
| Osmolytes & Cytoprotection | Trehalose, TMAO, Betaine, Taurine |
| Microbiome & Metabolites | Butyrate, Urolithin A, Indole |
| Senolytics & Longevity | Quercetin, Rapamycin, NAD+ |
| Resolution of Inflammation | Resolvin D1, Lipoxin A4, Maresin-1 |
| Gasotransmitters | GYY4137, CORM-2, AP39, HO-1 |
| Mechanosensing | Piezo1, TRPV4, TRPM7, ANO6 |
| Ferroptosis | GPX4, FSP1, ACSL4, RSL3 |
| m6A Epigenetics | METTL3, ALKBH5, YTHDF2 |
| Mitochondrial Metabolism | PGC1α, CoQ10, NAD precursors |
| Neuroprotection | Sigma-1R, neurosteroids |
| Natural Anti-inflammatories | Resveratrol, Curcumin, Luteolin |
| Drug Repurposing | Compounds in trials for other indications |

---

## Installation

### Requirements
```
Python 3.9+
streamlit>=1.35.0
biopython>=1.83
requests>=2.31.0
pandas>=2.0.0
plotly>=5.18.0
```

### Local
```bash
pip install -r requirements.txt
streamlit run app.py
```

### Streamlit Cloud (recommended)
1. Upload the 4 files to GitHub: `app.py`, `backend.py`, `constantes.py`, `requirements.txt`
2. Go to [share.streamlit.io](https://share.streamlit.io) → New app → select your repository
3. Main file path: `app.py` (or `app_doutorado.py` if you rename it)
4. Under **Settings → Secrets** add:
```toml
GEMINI_API_KEY = "AIza..."
```
> If no Secret is set, the app asks for the key manually — it stays in the session only, never in code.

---

## Usage

### 1. Configure
- **NCBI Email** — required for the PubMed API
- **Main Target** — supports full PubMed syntax:
  - `Bladder` — simple
  - `Bladder OR Urethra OR Kidney` — multiple organs
  - `Bladder NOT Cancer` — exclude context (useful to focus on normal physiology)
  - `(Bladder OR Urethra) NOT (Cancer OR Tumor)` — combined
- Select the **System** to auto-fill the anatomical context

### 2. Auto-Detect Targets
Click **AUTO-DETECT TARGETS**. The system shows real-time progress:
```
🔭 Searching 16 distant fields...
   📚 proteostasis chaperone stress... ✅ 168 candidates
   📚 osmolyte cytoprotection... ✅ 166 candidates
🧹 Removing genomic noise...
🧬 Verifying tissue expression (HPA)...
🔬 Checking PubMed virginity...
🤖 Flash evaluating biological plausibility...
✅ 87 candidates ready
```

### 3. Statistical Analysis
Click **RUN STATISTICAL ANALYSIS**. Calculates for each candidate:
- **λ Score** — enrichment ratio in the context
- **P-Value** — statistical significance (approximated Fisher)
- **Target Hits / Global Hits** — articles in the organ vs. all contexts
- **Trend** — ratio of recent articles (2023–2026)

### 4. Explore Results
- Filter by tag: 💎 Blue Ocean | 🌱 Embryonic | 🥇 Gold | 🚀 Trending
- **Click any row** → side panel with PubMed articles
- **📖 View Abstract** — free, no AI tokens
- **🤖 Investigate with AI** — Flash: `[Tissue] → [Mechanism] → [Effect]`
- **🧠 Investigate Target (Pro)** — deep analysis: druggability, novelty score, suggested experiments, biological rationale

---

## Tag Classification

| Tag | Criterion | Meaning |
|-----|-----------|---------|
| 💎 Blue Ocean | ≤5 hits in target, >20 global | Virgin in this context |
| 🌱 Embryonic | 6–25 hits in target | Emerging |
| 🥇 Gold | Enrichment >5× | High specificity ratio |
| 🚀 Trending | Enrichment >1.5× | Growing above average |
| 👻 Ghost | ≤5 global hits | Insufficient evidence |
| 🔴 Saturated | Enrichment ≤1.5× | Already well explored |

---

## Data Sources

| Source | Use | Cost |
|--------|-----|------|
| PubMed/NCBI | Lateral mining + virginity check | Free |
| OpenTargets | Genetically validated targets | Free |
| Human Protein Atlas | Tissue expression | Free |
| ClinicalTrials.gov | Drug repurposing | Free |
| Gemini 2.5 Flash | Classification + article analysis (lazy) | Tokens |
| Gemini 2.5 Pro | Deep target investigation (lazy) | Tokens |

> Flash is called once for classification and on-demand per article. Pro only when you click. Statistical analysis = zero tokens.

---

## File Structure

```
app.py            → Streamlit interface
backend.py        → Mining pipeline, statistical analysis, AI calls
constantes.py     → Blacklist, presets, systems, tags, models
requirements.txt  → Python dependencies
```

> You can rename `app.py` to `app_doutorado.py` — just update the Main file path in Streamlit Cloud settings.

---

## Citation

```
Lemos, G. (2025). Lemos Lambda v2.1: Deep Science Prospector.
Pharmacological prospecting tool via cross-field lateral literature mining.
```

Support the project: **Pix** `960f3f16-06ce-4e71-9b5f-6915b2a10b5a`

---

## License

Academic use. Developed for translational pharmacology research.
