# λ Lemos Lambda: Deep Science Prospector

> **Version:** 1.6.1  
> **Author:** Guilherme Lemos (Unifesp)  
> **License:** MIT  

## 🧬 Overview
**Lemos Lambda** is a pharmacological prospecting tool designed to identify emerging molecular targets in biomedical literature. Unlike standard search engines, it uses a heuristic **"Semantic Wall"** architecture to filter out 99% of structural noise (administrative metadata, clinical acronyms, grant codes) and isolate high-value biological signals.

It calculates a **Lambda Score** to rank targets based on the "Blue Ocean" strategy: finding molecules with high global relevance but low saturation in the specific target organ/tissue.

## 🚀 Key Features

* **Deep Mining:** Real-time extraction from PubMed (Title/Abstract) using Biopython.
* **Functional Context Filter:** automatically discards articles that do not contain mechanistic vocabulary (e.g., *signaling, pathway, relaxation*), ensuring biological relevance.
* **Inverted Lambda Score:** A metric that highlights "Opportunity vs. Saturation".
* **Architecture of Silence:** A 7-layer filtration pipeline to ensure clean data.

## 🛠️ Filtration Pipeline (The Semantic Wall)

1.  **Query Restriction:** Limits search to `[Title/Abstract]` to avoid affiliation/grant noise.
2.  **Type Exclusion:** `NOT (Review OR Editorial OR Comment)`.
3.  **Functional Context:** Article text must contain functional keywords (e.g., *activation, inhibition*).
4.  **Structural Regex:** Captures only alphanumeric patterns typical of gene symbols (e.g., *TRPV1, HO-1*).
5.  **Semantic Blacklist:** Blocks clinical acronyms (OAB, UTI), statistics (ANOVA, CI), and geography (USA, NSW).
6.  **Biological Whitelist:** Rescues short chemical entities (NO, H2S, CO) often missed by regex.
7.  **Frequency Cutoff:** Removes terms appearing in >25% of the corpus (too generic).

## 📦 Installation

1.  Clone this repository or extract the files.
2.  Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```
3.  Run the application:
    ```bash
    streamlit run app_doutorado.py
    ```

## 📄 Citation

If you use this software in your research, please cite:

> Lemos, G. (2025). *Lemos Lambda: Deep Science Prospector* [Software]. Version 1.6.1. Unifesp.