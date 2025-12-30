# Changelog

All notable changes to the "Lemos Lambda" project will be documented in this file.

## [2.0.0] - 2025-12-23

### 🚀 New Features
- **Statistical Engine:** Implemented Fisher's Exact Test to calculate the *Lambda Score*, providing statistical validation for target enrichment (p-value).
- **Generative Pharmacology (AI):** Integrated Google Gemini models to analyze abstracts. The AI now acts as a "Senior Researcher," extracting *Target -> Drug -> Physiological Effect* in natural language.
- **Universal API Compatibility:** Developed a "Cascade Fallback System" that automatically tests multiple Gemini versions (2.5, 2.0, 1.5, 1.0) to ensure compatibility with both legacy and new API keys.
- **Sherlock Holmes Mode (v2.10):** A new heuristic logic that infers pharmacological targets from the article Title when the Abstract is truncated or missing.
- **Scientific Radar:** New dashboard-style UI displaying the latest relevant news and findings with a modern card layout.
- **Editable Target List:** Users can now manually edit, copy, and paste the generated target list directly in the UI before execution.

### 🐛 Bug Fixes
- **Abstract Parsing:** Fixed a critical bug where the parser read only the first line of long PubMed abstracts; it now correctly handles multi-line indentation.
- **API Error Handling:** Resolved "Error 404" and "Error 400" by implementing better exception handling and user guidance for Google Cloud project creation.
- **Language Switching:** Fixed a crash that occurred when switching between English/Portuguese without clearing the previous results table.

### 🔧 Improvements
- **Prompt Engineering:** Refined the AI prompt (v2.9) to prevent "hallucinations" (e.g., writing emails) and focus strictly on data extraction.
- **Performance:** Optimized the Biopython search query to reduce latency during deep mining.