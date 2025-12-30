
import { GoogleGenAI } from "@google/genai";
import { MiningStrategy } from "../types";

export class GeminiService {
  private ai: GoogleGenAI;

  constructor() {
    this.ai = new GoogleGenAI({ apiKey: process.env.API_KEY });
  }

  /**
   * Deep Mining Strategy:
   * Branches based on the selected strategy.
   */
  async mineNovelTargets(target: string, context: string, currentList: string[], strategy: MiningStrategy): Promise<string[]> {
    if (!this.ai) {
      console.warn("Gemini AI instance not initialized.");
      throw new Error("Gemini Client Error");
    }

    try {
      // If list is empty, we say "None".
      const currentListStr = currentList.length > 0 ? currentList.join(", ") : "None";
      
      const baseContext = `
      CONTEXT DATA:
      - Primary Disease/Organ: "${target}"
      - Additional Biological Context: "${context}"
      - Existing List (DO NOT REPEAT THESE): ${currentListStr}
      `;

      let prompt = "";
      let temperature = 0.5;
      let topK = 40;

      switch (strategy) {
        case 'conservative':
          temperature = 0.1; // Strict, deterministic
          prompt = `
            Role: Scientific Data Curator.
            Task: Clean and standardize the "Existing List" provided above.
            ${baseContext}
            
            Instructions:
            1. **REMOVE NOISE**: Delete terms that are NOT specific molecular targets (e.g., "Study", "Analysis", "Rats", "Patients", "Placebo", "Safety").
            2. **STANDARDIZE**: Convert names to official gene/protein symbols (e.g., "ATP receptor" -> "P2X3").
            3. **OUTPUT**: A strict comma-separated list of the cleaned molecules. Nothing else.
          `;
          break;

        case 'repurposing':
          temperature = 0.4; // Balanced for finding real drugs
          prompt = `
            Role: Senior Pharmacologist & Drug Hunter.
            Task: Suggest approved drugs or clinical candidates that could be REPURPOSED for "${target}".
            ${baseContext}

            Instructions:
            1. **IGNORE EXISTING**: Do not output items from the "Existing List". Find NEW candidates.
            2. **HIGH POTENTIAL**: Focus on drugs with known safety profiles (FDA approved) or Phase 2/3 data.
            3. **SPECIFICITY**: Return the DRUG NAME or SPECIFIC MOLECULE (e.g., "SGLT2 inhibitors", "Metformin", "Mirabegron").
            4. **OUTPUT**: A strict comma-separated list of 15 unique, high-value candidates. No numbering.
          `;
          break;

        case 'mechanism':
          temperature = 0.5; // Balanced for pathways
          prompt = `
            Role: Molecular Systems Biologist.
            Task: Identify upstream regulators and downstream effectors for "${target}".
            ${baseContext}

            Instructions:
            1. **PATHWAY MAPPING**: Suggest Kinases (e.g., mTOR, PI3K), Transcription Factors (NF-kB), or Ion Channels involved in the pathophysiology.
            2. **DEEP DIVE**: Look for specific subunits (e.g., "G alpha s", "Nav1.8", "TRPV4") rather than generic terms.
            3. **OUTPUT**: A strict comma-separated list of 15 mechanistic targets. No numbering.
          `;
          break;

        case 'blue_ocean':
          temperature = 0.9; // High creativity for discovery
          topK = 60;
          prompt = `
            Role: Elite Scientific Prospector.
            Task: Identify NOVEL, CONTROVERSIAL, or EMERGING targets for "${target}" (Blue Ocean Strategy).
            ${baseContext}
            
            Instructions:
            1. **NOVELTY**: Prioritize targets appearing in literature only in the last 3-5 years or those with low citation counts but high relevance.
            2. **IGNORE THE OBVIOUS**: Do NOT suggest well-known targets for this condition. Dig deeper.
            3. **ANATOMICAL NICHE**: Look for targets specific to cell types (e.g., "Interstitial Cells", "Urothelium", "Endothelium").
            4. **OUTPUT**: A strict comma-separated list of 15 novel targets. No bullet points.
          `;
          break;
      }

      // Using gemini-3-flash-preview for speed and reliability with lists
      const response = await this.ai.models.generateContent({
        model: 'gemini-3-flash-preview', 
        contents: prompt,
        config: {
          temperature: temperature,
          topK: topK,
        }
      });

      const text = response.text || "";
      console.log(`Gemini Raw Response (${strategy}):`, text); // Debugging
      
      // FIX: Robust Parsing
      // Remove common conversational fillers and Markdown
      const cleanText = text
        .replace(/^Here.*:/i, "")
        .replace(/Output:/i, "")
        .replace(/[\[\]*_`]/g, "") // Remove brackets, bold, italics, code blocks
        .trim();
      
      // Split by comma or newline
      const terms = cleanText.split(/,|\n/);
      
      return terms
        .map(t => t.trim())
        .map(t => t.replace(/^\d+[\).]\s*/, "")) // Remove "1. ", "2)"
        .filter(t => t.length > 2 && t.length < 50) // Length filter
        .filter(t => !t.toLowerCase().includes("context")) // Filter out hallucinations repeating "Context"
        .filter(t => !t.toLowerCase().includes("existing list")) // Filter out hallucinations
        .filter(t => !t.toLowerCase().includes("none")); // Filter out empty indicator
        
    } catch (error) {
      console.error("Gemini Mining Error:", error);
      throw error; // Throw so UI knows it failed
    }
  }

  async analyzePaper(title: string, abstract: string): Promise<string> {
    if (!this.ai) return "API Key Required for Analysis.";
    
    if (!abstract || typeof abstract !== 'string') {
      return "Abstract not available for analysis.";
    }

    let optimizedAbstract = abstract;
    if (abstract.length > 800) {
      optimizedAbstract = abstract.substring(0, 400) + "\n...\n" + abstract.substring(abstract.length - 400);
    }

    try {
      const response = await this.ai.models.generateContent({
        model: 'gemini-3-flash-preview',
        contents: `
        Analyze this scientific abstract.
        
        Title: "${title}"
        Text: "${optimizedAbstract}"
        
        Task: Identify the main MOLECULAR TARGET, the DRUG/MOLECULE used, and the EFFECT.
        Output exactly one short sentence: "M: [Molecule] | E: [Effect]"
        `,
      });

      return response.text || "Analysis failed.";
    } catch (error) {
      return "Error analyzing paper (AI Limit).";
    }
  }
}
