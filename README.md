# λ Lemos Lambda: Deep Science Prospector

> **Version:** 2.0 (Stable)  
> **Author:** Guilherme Lemos (Unifesp)  
> **License:** MIT  

## 🧬 Overview
**Lemos Lambda** is a pharmacological prospecting tool designed to identify emerging molecular targets in biomedical literature. Integrating **Generative AI (Google Gemini)** with inferential statistics (**Fisher's Exact Test**), it goes beyond standard search engines.

It employs a heuristic **"Semantic Wall"** architecture to filter out 99% of structural noise (administrative metadata, clinical acronyms, grant codes) and isolate high-value biological signals, classifying them into strategic zones (e.g., "Blue Ocean", "Saturated").

## 🚀 Key Features

* **Deep Mining:** Real-time extraction from PubMed (Title/Abstract) using Biopython.
* **Generative Pharmacology (New in v2.0):** Uses LLMs (Gemini 2.5/1.5) to read abstracts and extract: *Target -> Drug -> Physiological Effect*.
* **Sherlock Holmes Mode:** Capable of inferring pharmacological targets even from incomplete or truncated abstracts.
* **Fisher's Enrichment Score:** A statistical metric that highlights "Opportunity vs. Saturation" with p-value validation.
* **Architecture of Silence:** A multi-layer filtration pipeline to ensure clean data.

## 🛠️ Filtration Pipeline (The Semantic Wall)

1.  **Query Restriction:** Limits search to `[Title/Abstract]` to avoid affiliation/grant noise.
2.  **Type Exclusion:** `NOT (Review OR Editorial OR Comment)`.
3.  **Functional Context Filter:** Articles must contain mechanistic vocabulary (e.g., *signaling, pathway, activation*).
4.  **Structural Regex:** Captures only alphanumeric patterns typical of gene symbols (e.g., *TRPV1, HO-1*).
5.  **Semantic Blacklist:** Blocks clinical acronyms (OAB, UTI), statistics (ANOVA, CI), and geography (USA, NSW).
6.  **Biological Whitelist:** Rescues short chemical entities (NO, H2S, CO) often missed by regex.

## 📦 Installation

1.  Clone this repository or extract the files.
2.  Install dependencies:
    ```bash
    pip install streamlit pandas plotly scipy biopython google-generativeai tenacity
    ```
3.  Run the application:
    ```bash
    streamlit run app_doutorado.py
    ```

## 🔑 AI Configuration (Essential)

To enable the **Generative AI** features (preventing Error 404), follow these steps:

1.  Go to [Google AI Studio](https://aistudio.google.com/app/apikey).
2.  Click **"Create API Key"**.
3.  **IMPORTANT:** Select **"Create API key in NEW project"**.
4.  Copy the key and paste it into the Lemos Lambda settings sidebar.

## 📄 Citation

If you use this software in your research, please cite:

> Lemos, G. (2025). *Lemos Lambda: Deep Science Prospector* [Software]. Version 2.0. Unifesp.