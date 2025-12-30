# Changelog - Lemos Lambda

All notable changes to the **"Lemos Lambda"** project will be documented in this file.

## [2.1.0] - 2025-12-24 (React Edition)

### 🚀 Major Architecture Overhaul
- **React Migration:** Complete codebase rewrite from Python (Streamlit) to **React 19 + TypeScript + Vite**. The application now runs entirely on the client-side, offering instant feedback and 60fps animations.
- **Performance:** Removed server-side bottlenecks. The UI no longer reloads on every interaction.
- **Vercel Deployment:** Optimized for edge deployment on Vercel/Netlify.

### 🧠 AI & Logic
- **SDK Upgrade:** Migrated from Python `google-generativeai` to the new **`@google/genai` Web SDK (v1.34)**.
- **Client-Side Mining:** Implemented `services/scienceService.ts` to handle regex extraction and statistical processing directly in the browser.
- **Deep Mining Strategies:** Added "Mechanism", "Repurposing", and "Blue Ocean" specific prompt pipelines in `GeminiService`.

### ✨ UI/UX
- **TailwindCSS:** Modern, responsive design with dark mode natively supported.
- **Interactive Charts:** Replaced static Matplotlib images with interactive **Recharts** visualizations.
- **Drill-Down Modal:** Clicking "Analyze with AI" now opens a modal overlay instead of reloading the page.

---

## [2.0.0] - 2025-12-23 (Stable Python Edition)

### 🚀 New Features
- **Statistical Engine:** Implemented Fisher's Exact Test to calculate the **Lambda Score**, providing statistical validation for target enrichment (p-value).
- **Generative Pharmacology (AI):** Integrated Google Gemini models to analyze abstracts. The AI now acts as a "Senior Researcher," extracting **Target -> Drug -> Physiological Effect** in natural language.
- **RESTful Universal Compatibility:** Developed a "Cascade Fallback System" via direct HTTP/REST requests. This automatically tests multiple Gemini versions (2.5, 2.0, 1.5, Pro) to bypass legacy library issues and ensure 100% uptime.
- **Scientific Radar:** New dashboard-style UI displaying the latest relevant news and findings with a modern card layout.

### 🐛 Bug Fixes
- **Visual Dark Mode:** Fixed the "white text on white background" bug in the AI analysis box; implemented high-contrast pharmacological labels.
- **API Error Handling:** Resolved "Error 404", "400" and "IA Ocupada" by bypassing the `google-generativeai` package and implementing direct endpoint communication.

---

## [1.0.0] - 2025-11-05
- **Initial Release:** Foundation of the Semantic Wall and basic PubMed mining.