import { GoogleGenAI } from "@google/genai";
import { MiningStrategy } from "../types";

export class GeminiService {
  private ai: GoogleGenAI | null = null;

  constructor(apiKey: string) {
    if (apiKey) {
      this.ai = new GoogleGenAI({ apiKey });
    }
  }

  updateKey(apiKey: string) {
    this.ai = new GoogleGenAI({ apiKey });
  }

  /**
   * Deep Mining Strategy:
   * Branches based on the selected strategy.
   */
  async mineNovelTargets(target: string, context: string, currentList: string[], strategy: MiningStrategy): Promise<string[]> {
    if (!this.ai) return currentList;

    try {
      const isComparative = (context || "").length > 2;
      const baseContext = `- Primary Target Organ/Disease: "${target}"\n- Current Candidate List: ${currentList.join(", ")}`;
      let prompt = "";

      switch (strategy) {
        case 'conservative':
          prompt = `
            Role: Biomedical Data Curator.
            Task: Clean and standardize a list of pharmacological targets for: "${target}".
            ${baseContext}
            
            Rules:
            1. **REMOVE GARBAGE**: Delete generic terms like "Study", "Analysis", "Expression", "Rat", "Patient", "Method".
            2. **STANDARDIZE**: Convert "ATP receptor" to "P2X3", "Beta3-AR" to "Beta-3 adrenergic receptor".
            3. **NO NEW TARGETS**: Only clean the existing list.
            
            Output: Comma-separated list of cleaned targets.
          `;
          break;

        case 'repurposing':
          prompt = `
            Role: Drug Repurposing Specialist.
            Task: Identify EXISTING DRUGS or PHARMACOLOGICAL TOOLS that could be repurposed for "${target}".
            ${baseContext}
            ${isComparative ? `- Context: "${context}"` : ''}

            Rules:
            1. **Focus on DRUGS/MOLECULES**: Look for FDA-approved drugs, Phase 2/3 candidates, or well-known inhibitors.
            2. **Ignore Basic Biology**: Do not suggest "Inflammation" or "Apoptosis". Suggest specific drug targets (e.g., "GLP-1R", "SGLT2", "JAK1").
            3. **High Translatability**: Suggest targets with available ligands.
            
            Output: Comma-separated list of top 15 repurposing candidates.
          `;
          break;

        case 'mechanism':
          prompt = `
            Role: Molecular Biologist.
            Task: Identify UPSTREAM and DOWNSTREAM molecular pathways/mechanisms for "${target}".
            ${baseContext}

            Rules:
            1. **Focus on PATHWAYS**: Kinases (mTOR, PI3K), Transcription Factors (NF-kB, NRF2), Enzymes (PDE5, COX-2).
            2. **Deep Mechanism**: Look for specific subunits (e.g., instead of "G-protein", suggest "G alpha s").
            3. **Avoid Drugs**: Focus on the biological target, not the commercial drug name.
            
            Output: Comma-separated list of top 15 mechanistic targets.
          `;
          break;

        case 'blue_ocean':
          prompt = `
            Role: Elite Data Scientist & Explorer.
            Mission: Identify UNDERSTUDIED (Blue Ocean) targets for: "${target}".
            ${baseContext}
            
            Rules:
            1. **ANATOMICAL NICHE**: Focus on specific cell types (e.g., Podocytes, Interstitial Cells) within the target.
            2. **NOVELTY**: Prefer targets discovered in the last 5 years.
            3. **AVOID SATURATION**: Do NOT suggest "Insulin", "TNF", "IL-6".
            4. **BE SPECIFIC**: Good: "P2Y12", "Piezo2". Bad: "Purinergic receptor".
            
            Output: Comma-separated list of top 15 novel targets.
          `;
          break;
      }

      const response = await this.ai.models.generateContent({
        model: 'gemini-3-pro-preview', 
        contents: prompt,
      });

      const text = response.text || "";
      
      // FIX 1: Robust Parsing
      // AI often ignores "comma-separated" and uses bullets (* or -) or newlines.
      // We remove Markdown headers, brackets, bullets, and then split by newline OR comma.
      const cleanText = text.replace(/Output:|Here is the list:|\[|\]|\*|- /g, "");
      
      // Split by comma OR newline to handle "List format" responses
      const terms = cleanText.split(/,|\n/);
      
      return terms
        .map(t => t.trim())
        // Clean up numbering (e.g., "1. Target")
        .map(t => t.replace(/^\d+\.\s*/, ""))
        // Filter out if it contains a space AND is too long (hallucination sentence)
        .filter(t => t.length > 2 && !(t.includes(" ") && t.length > 30));
        
    } catch (error) {
      console.error("Gemini Mining Error:", error);
      // Fail gracefully returning current list
      return currentList;
    }
  }

  async analyzePaper(title: string, abstract: string): Promise<string> {
    if (!this.ai) return "API Key Required for Analysis.";
    
    // FIX 2: Crash Prevention for Null Abstracts
    if (!abstract || typeof abstract !== 'string') {
      return "Abstract not available for analysis.";
    }

    // SMART CONTEXT OPTIMIZATION:
    let optimizedAbstract = abstract;
    if (abstract.length > 600) {
      optimizedAbstract = abstract.substring(0, 300) + "\n...[middle section skipped]...\n" + abstract.substring(abstract.length - 300);
    }

    try {
      const response = await this.ai.models.generateContent({
        model: 'gemini-3-flash-preview',
        contents: `
        Act as a Pharmacologist. Analyze this truncated abstract.
        
        Title: "${title}"
        Text: "${optimizedAbstract}"
        
        Task:
        Identify the **Molecule/Target**, its **Action** (Agonist/Antagonist), and the biological **Outcome**.
        
        Output format (one sentence):
        "Mechanism: [Molecule] acts as [Mechanism] -> [Outcome]"
        `,
      });

      return response.text || "Analysis failed.";
    } catch (error) {
      return "Error analyzing paper (AI Limit or Network).";
    }
  }
}