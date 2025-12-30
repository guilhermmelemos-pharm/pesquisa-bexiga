# λ Lemos Lambda: Deep Science Prospector

> **Version:** 2.1 (React Edition - 2025)  
> **Author:** Guilherme Lemos (Unifesp)  
> **Stack:** React, TypeScript, TailwindCSS, Google GenAI SDK  
> **License:** MIT  

## 🧬 Overview
**Lemos Lambda v2.1** is a high-performance prospecting system designed to navigate the "infodemic" of biomedical literature. It identifies emerging molecular targets (receptors, ion channels, enzymes) by crossing **Biomedical Big Data** with **Generative AI (Google Gemini 2.0/3.0)** and inferential statistics.

Unlike the previous Python version, the **React Edition** runs entirely in the user's browser, providing a zero-latency experience for data mining and visualization.

## 🚀 Key Features (v2.1 Updates)

*   **Client-Side Engine:** No server processing required. All regular expressions and statistical calculations (Fisher's Exact Test simulation) happen locally.
*   **Gemini Web SDK:** Direct integration with `@google/genai` for analyzing abstracts without intermediate backends.
*   **Safety-Off Protocol:** Specifically tuned for Pharmacology. It bypasses standard AI medical filters to allow deep analysis of terms like "toxicity," "bladder dysfunction," and "cystitis."
*   **Interactive Heatmaps:** Dynamic visualization of "Blue Ocean" vs. "Saturated" targets using Recharts.

## 🛠️ The Semantic Wall (Pipeline v2.1)

To ensure a clean signal, the system employs a multi-layer filtration architecture:

1.  **AI Strategy Selector:** Choose between *Cleaning*, *Repurposing*, or *Mechanism Discovery* modes.
2.  **Denoising Regex:** Javascript-based extraction of biological entities (e.g., *TRPV1, P2X7, FAAH*).
3.  **Blacklist 2.1:** Blocks "False Positives" like ANOVA, CI, USA, and clinical acronyms.
4.  **Lambda Score:** Calculates the Enrichment Ratio ($ER$) vs. Global Hits.

## 📦 Installation & Setup

This is a **React/Vite** project. You need Node.js installed.

1.  **Install Dependencies:**
    ```bash
    npm install
    ```

2.  **Run Locally:**
    ```bash
    npm run dev
    ```

3.  **Build for Production (Vercel/Netlify):**
    ```bash
    npm run build
    ```

## 🔑 AI Configuration

To use the **Deep Mining** and **Abstract Analysis** features:

1.  Get your key at [Google AI Studio](https://aistudio.google.com/app/apikey).
2.  Paste it into the **Settings** panel within the app.
3.  The key is stored in your browser's memory for the session and is **never** sent to our servers (direct client-to-Google communication).

## 📊 Mathematical Logic
The system calculates the **Enrichment Ratio ($ER$)** as:

$$ER = \frac{\text{Observed specific hits}}{\text{Expected global hits}}$$

Where expected hits are derived from a contingency table comparing target hits within the user-defined context vs. the total PubMed database ($N \approx 36,000,000$).

## 📄 Citation

If you use this software in your research, please cite:

> Lemos, G. (2025). Lemos Lambda: Deep Science Prospector (v2.1). Zenodo. https://doi.org/10.5281/zenodo.18036690