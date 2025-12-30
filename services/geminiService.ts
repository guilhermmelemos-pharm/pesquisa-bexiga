
import { GoogleGenAI } from "@google/genai";
import { MiningStrategy } from "../types";

export class GeminiService {
  
  private getClient(apiKey?: string): GoogleGenAI {
    // Prioritize key passed from UI, fallback to env, finally error if neither exists
    const key = apiKey || process.env.API_KEY;
    if (!key) {
      throw new Error("MISSING_KEY");
    }
    return new GoogleGenAI({ apiKey: key });
  }

  /**
   * Deep Mining Strategy:
   * Branches based on the selected strategy.
   */
  async mineNovelTargets(target: string, context: string, currentList: string[], strategy: MiningStrategy, apiKey?: string): Promise<string[]> {
    let ai: GoogleGenAI;
    try {
        ai = this.getClient(apiKey);
    } catch (e) {
        console.warn("API Key missing");
        throw new Error("API Key Missing");
    }

    try {
      // Create a set for case-insensitive local filtering of the output
      const lowerCaseCurrentSet = new Set(currentList.map(t => t.toLowerCase().trim()));
      
      const currentListStr = currentList.length > 0 ? currentList.join(", ") : "None";
      
      const baseContext = `
      CONTEXT:
      - Primary Disease/Organ: "${target}"
      - Bio Context: "${context}"
      
      EXCLUSION LIST (ALREADY KNOWN - DO NOT OUTPUT THESE):
      [${currentListStr}]
      `;

      let prompt = "";
      let temperature = 0.5;
      let topK = 40;

      switch (strategy) {
        case 'conservative':
          temperature = 0.1; // Strict
          prompt = `
            Role: Scientific Data Curator.
            Task: Clean, Standardize, and Convert the "EXCLUSION LIST" into a clean list of Molecular Targets.
            
            Input Data:
            ${baseContext}
            
            Instructions:
            1. **CONVERT DRUGS TO TARGETS**: If an item is a drug (e.g., "Semaglutide", "Viagra"), output its primary molecular target (e.g., "GLP1R", "PDE5").
            2. **REMOVE NOISE**: Delete terms that are NOT specific molecular targets (e.g., "Study", "Analysis", "Safety", "Placebo", "Rats").
            3. **STANDARDIZE**: Use official gene/protein symbols (e.g., "ATP receptor" -> "P2X3").
            4. **OUTPUT**: A single comma-separated list of the FINAL targets.
          `;
          break;

        case 'repurposing':
          temperature = 0.4;
          prompt = `
            Role: Senior Pharmacologist.
            Task: Suggest approved drugs or clinical candidates for "${target}".
            ${baseContext}

            Instructions:
            1. **STRICTLY NEW**: Do NOT output any item present in the "EXCLUSION LIST".
            2. **HIGH POTENTIAL**: Focus on drugs with FDA approval or Phase 2/3 data.
            3. **FORMAT**: Return the Drug Name or Specific Molecule.
            4. **OUTPUT**: A strict comma-separated list of 15 unique, high-value candidates.
          `;
          break;

        case 'mechanism':
          temperature = 0.5;
          prompt = `
            Role: Molecular Systems Biologist.
            Task: Identify upstream regulators and downstream effectors for "${target}".
            ${baseContext}

            Instructions:
            1. **STRICTLY NEW**: Do NOT output any item present in the "EXCLUSION LIST".
            2. **PATHWAY MAPPING**: Suggest Kinases, Transcription Factors, or Ion Channels.
            3. **DEEP DIVE**: Specific subunits (e.g., "G alpha s", "Nav1.8") rather than generic terms.
            4. **OUTPUT**: A strict comma-separated list of 15 mechanistic targets.
          `;
          break;

        case 'blue_ocean':
          temperature = 0.9;
          topK = 60;
          prompt = `
            Role: Elite Scientific Prospector.
            Task: Identify NOVEL, CONTROVERSIAL, or EMERGING targets for "${target}" (Blue Ocean Strategy).
            ${baseContext}
            
            Instructions:
            1. **STRICTLY NEW**: Do NOT output any item present in the "EXCLUSION LIST". This is critical.
            2. **NOVELTY**: Targets from literature in the last 3-5 years.
            3. **IGNORE THE OBVIOUS**: Dig for obscure receptors, non-coding RNAs, or orphan GPCRs.
            4. **OUTPUT**: A strict comma-separated list of 15 novel targets.
          `;
          break;
      }

      const response = await ai.models.generateContent({
        model: 'gemini-3-flash-preview', 
        contents: prompt,
        config: {
          temperature: temperature,
          topK: topK,
        }
      });

      const text = response.text || "";
      console.log(`Gemini Raw Response (${strategy}):`, text);
      
      const cleanText = text
        .replace(/^Here.*:/i, "")
        .replace(/Output:/i, "")
        .replace(/[\[\]*_`]/g, "")
        .trim();
      
      const terms = cleanText.split(/,|\n/);
      
      return terms
        .map(t => t.trim())
        .map(t => t.replace(/^\d+[\).]\s*/, ""))
        .filter(t => t.length > 2 && t.length < 50)
        .filter(t => !t.toLowerCase().includes("context"))
        .filter(t => !t.toLowerCase().includes("exclusion list"))
        .filter(t => !t.toLowerCase().includes("none"))
        // Strong client-side filter to prevent repetition
        .filter(t => {
            // If strategy is conservative, we ALLOW items from input (because we are cleaning/converting them)
            // If strategy is OTHERS, we STRICTLY BLOCK items from input
            if (strategy === 'conservative') return true;
            
            // For other strategies, check against the exclusion list
            return !lowerCaseCurrentSet.has(t.toLowerCase());
        });
        
    } catch (error) {
      console.error("Gemini Mining Error:", error);
      throw error;
    }
  }

  async analyzePaper(title: string, abstract: string, apiKey?: string): Promise<string> {
    let ai: GoogleGenAI;
    try {
        ai = this.getClient(apiKey);
    } catch (e) {
        return "API Key Required for Analysis.";
    }
    
    if (!abstract || typeof abstract !== 'string') {
      return "Abstract not available for analysis.";
    }

    let optimizedAbstract = abstract;
    if (abstract.length > 800) {
      optimizedAbstract = abstract.substring(0, 400) + "\n...\n" + abstract.substring(abstract.length - 400);
    }

    try {
      const response = await ai.models.generateContent({
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
