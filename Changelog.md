# Changelog - Lemos Lambda

All notable changes to the **"Lemos Lambda"** project will be documented in this file.

## [2.0.0] - 2025-12-23 (Stable Edition)

### 🚀 New Features
- **Statistical Engine:** Implemented Fisher's Exact Test to calculate the **Lambda Score**, providing statistical validation for target enrichment (p-value).
- **Generative Pharmacology (AI):** Integrated Google Gemini models to analyze abstracts. The AI now acts as a "Senior Researcher," extracting **Target -> Drug -> Physiological Effect** in natural language.
- **RESTful Universal Compatibility:** Developed a "Cascade Fallback System" via direct HTTP/REST requests. This automatically tests multiple Gemini versions (2.5, 2.0, 1.5, Pro) to bypass legacy library issues and ensure 100% uptime.
- **Sherlock Holmes Mode (v2.10):** A new heuristic logic that infers pharmacological targets from the article Title when the Abstract is truncated or missing.
- **Scientific Radar:** New dashboard-style UI displaying the latest relevant news and findings with a modern card layout.
- **On-Demand Intelligence:** Analysis buttons integrated into individual article cards to optimize API quota usage and allow surgical investigation of targets.
- **Editable Target List:** Users can now manually edit, copy, and paste the generated target list directly in the UI before execution (Label: "termos indicados").

### 🐛 Bug Fixes
- **Visual Dark Mode:** Fixed the "white text on white background" bug in the AI analysis box; implemented high-contrast pharmacological labels.
- **Abstract Parsing:** Fixed a critical bug where the parser read only the first line of long PubMed abstracts; it now correctly handles multi-line indentation and captures keywords.
- **API Error Handling:** Resolved "Error 404", "400" and "IA Ocupada" by bypassing the `google-generativeai` package and implementing direct endpoint communication.
- **Language Switching:** Fixed a crash that occurred when switching between English/Portuguese without clearing the previous results table.

### 🔧 Improvements
- **Prompt Engineering (v2.9):** Refined the AI prompt to prevent "hallucinations" and strictly enforce the pharmacological extraction rules with **Safety-Off** settings for medical terms.
- **UX Refinement:** Moved the API Key configuration to the results screen, allowing "Just-in-Time" setup without losing the search context.
- **Performance:** Optimized the Biopython search query to reduce latency during deep mining and improved the "Semantic Wall" to block administrative metadata.

---

## [1.0.0] - 2025-11-05
- **Initial Release:** Foundation of the Semantic Wall and basic PubMed mining.
