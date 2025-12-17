# λ Lemos Lambda: Deep Science Prospector

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://lemosgbuscador.streamlit.app/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17958507.svg)](https://doi.org/10.5281/zenodo.17958507)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Lemos Lambda** is a scientific competitive intelligence tool designed for pharmacology and physiology researchers. It automates PubMed data mining to identify "Blue Oceans"—promising therapeutic targets that are widely studied in general physiology but neglected in specific organs (e.g., Bladder), identifying high-impact research opportunities.

## 🚀 Key Features

* **Opportunity Ratio (Ratio):** Automatically calculates the publication ratio between a "Super Source" (e.g., Brain, Kidney) and your "Target Organ".
* **'Blue Ocean' Mining:** A built-in algorithm that screens a curated list of ~60 innovative targets (Orphan GPCRs, Piezo1/2, Ferroptosis) against your specific tissue to find rare opportunities.
* **Bibliographic X-Ray:** Reads, summarizes, and translates relevant abstracts automatically.
* **Smart Classification:** Tags terms as "Diamond" (Unexplored Opportunity), "Gold", or "Saturated" based on publication volume.

## 🛠️ How to Use

### Option 1: Web App (No installation required)
Access the tool directly via your browser:
👉 **[https://lemosgbuscador.streamlit.app/](https://lemosgbuscador.streamlit.app/)**

### Option 2: Run Locally (For developers)
1.  Clone this repository:
    ```bash
    git clone [https://github.com/guilhermelemos-pharm/pesquisa-bexiga.git](https://github.com/guilhermelemos-pharm/pesquisa-bexiga.git)
    ```
2.  Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```
3.  Run the application:
    ```bash
    streamlit run app_doutorado.py
    ```

## 📄 How to Cite

If **Lemos Lambda** has assisted your research, please cite it. This is crucial for supporting open-source academic software.

> Lemos, G. (2025). *Lemos Lambda: Deep Science Prospector* [Software]. Available at: https://github.com/guilhermelemos-pharm/pesquisa-bexiga

Or use the `CITATION.cff` file in this repository to export directly to BibTeX or EndNote.

## ⚖️ License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---
Developed by **Guilherme Lemos** (Federal University of São Paulo - Unifesp)
