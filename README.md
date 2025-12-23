# λ Lemos Lambda: Deep Science Prospector

> **Version:** 2.0 (Stable - 2025 Edition)  
> **Author:** Guilherme Lemos (Unifesp)  
> **Field:** Pharmacology & Molecular Urology  
> **License:** MIT  

## 🧬 Overview
**Lemos Lambda** is a high-performance prospecting system designed to navigate the "infodemic" of biomedical literature. It identifies emerging molecular targets (receptors, ion channels, enzymes) by crossing **PubMed Big Data** with **Generative AI (Google Gemini 2.5/2.0)** and **Fisher's Exact Test** inferential statistics.

Unlike traditional keyword searches, Lemos Lambda quantifies the **Enrichment Ratio**, allowing researchers to distinguish between "Saturated Fields" and "Blue Ocean Opportunities" in pharmacology.



## 🚀 Key Features (v2.0 Updates)

* **RESTful Core:** Bypasses library dependencies by using direct HTTP/REST calls to Google’s servers, ensuring 100% uptime regardless of server environment.
* **Pharmacological On-Demand Analysis:** Surgically analyze high-potential papers: *Target → Drug → Functional Effect*.
* **Safety-Off Protocol:** Specifically tuned for Pharmacology. It bypasses standard AI medical filters to allow deep analysis of terms like "toxicity," "bladder dysfunction," and "cystitis."
* **Fisher’s Enrichment Score ($P < 0.05$):** Validates the significance of a target's appearance within a specific context (e.g., Piezo2 in Bladder) vs. the global literature.

## 🛠️ The Semantic Wall (Pipeline v2.0)

To ensure a clean signal, the system employs a 6-layer filtration architecture:

1.  **Query Precision:** Metadata isolation in `[Title/Abstract]`.
2.  **Mechanistic Bias:** Only articles with functional keywords (*pathway, signaling, activation*) are analyzed.
3.  **Denoising Regex:** Dynamic capture of biological entities (e.g., *TRPV1, P2X7, FAAH*).
4.  **Blacklist 2.0:** Blocks "False Positives" like ANOVA, CI, USA, and clinical acronyms (UTI, OAB).
5.  **Dynamic Whitelist:** Rescues essential gasotransmitters (NO, $H_{2}S$, CO).



## 📦 Installation & Setup

1.  **Environment:**
    ```bash
    pip install streamlit pandas plotly scipy biopython requests tenacity
    ```
2.  **Execution:**
    ```bash
    streamlit run app_doutorado.py
    ```

## 🔑 AI Configuration (The "2025 Protocol")

If you encounter **Error 404** or **Incompatibility**, ensure your API Key is generated correctly:

1.  Access [Google AI Studio](https://aistudio.google.com/app/apikey).
2.  Click **"Create API Key"** and select a **New Project**.
3.  Copy the key and paste it into the Lemos Lambda settings. The system will automatically detect the best available model (e.g., **Gemini 2.5 Flash**).

## 📊 Mathematical Logic
The system calculates the **Enrichment Ratio ($ER$)** as:

$$ER = \frac{\text{Observed specific hits}}{\text{Expected global hits}}$$

Where $p$-values are derived from a contingency table (Fisher’s Exact Test) comparing target hits within the user-defined context vs. the total PubMed database ($N \approx 36,000,000$).



## 📄 Citation

If you use this software in your research, please cite:

> Lemos, G. (2025). Lemos Lambda: Deep Science Prospector (v2.0). Zenodo. https://doi.org/10.5281/zenodo.18036690
