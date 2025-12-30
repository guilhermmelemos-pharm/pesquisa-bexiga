import { GoogleGenAI } from "@google/genai";
import { MiningStrategy } from "../types";

export class GeminiService {
  private ai: GoogleGenAI | null = null;

  // Mantivemos o construtor para pegar a chave do Input do usuário
  constructor(apiKey: string) {
    if (apiKey) {
      this.ai = new GoogleGenAI({ apiKey });
    }
  }

  updateKey(apiKey: string) {
    this.ai = new GoogleGenAI({ apiKey });
  }

  async mineNovelTargets(target: string, context: string, currentList: string[], strategy: MiningStrategy): Promise<string[]> {
    if (!this.ai) {
        console.warn("API Key missing");
        return currentList;
    }

    try {
      // 1. Prepara a lista de exclusão para a IA não repetir
      const currentListStr = currentList.length > 0 ? currentList.join(", ") : "None";
      
      const baseContext = `
      CONTEXT:
      - Primary Disease/Organ: "${target}"
      - Bio Context: "${context}"
      
      EXCLUSION LIST (DO NOT OUTPUT THESE):
      [${currentListStr}]
      `;

      let prompt = "";
      let temperature = 0.5;

      // Seus prompts excelentes (mantidos)
      switch (strategy) {
        case 'conservative':
          temperature = 0.1;
          prompt = `
            Role: Scientific Data Curator.
            Task: Clean and Standardize the candidate list.
            Input: ${baseContext}
            Instructions: 
            1. Convert Drugs to Targets (e.g. Semaglutide -> GLP1R).
            2. Remove noise. 
            3. Output ONLY the clean comma-separated list.
          `;
          break;

        case 'repurposing':
          temperature = 0.6;
          prompt = `
            Role: Pharmacologist.
            Task: Suggest approved drugs/candidates for "${target}".
            ${baseContext}
            Instructions:
            1. STRICTLY NEW items (Not in Exclusion List).
            2. Focus on FDA approved/Phase 2-3.
            3. Output: Comma-separated list of 30 items.
          `;
          break;

        case 'mechanism':
          temperature = 0.6;
          prompt = `
            Role: Systems Biologist.
            Task: Identify upstream/downstream pathways for "${target}".
            ${baseContext}
            Instructions:
            1. STRICTLY NEW items.
            2. Focus on Kinases, Transcription Factors, Ion Channels.
            3. Output: Comma-separated list of 30 items.
          `;
          break;

        case 'blue_ocean':
          temperature = 0.9; 
          prompt = `
            Role: Elite Scientific Prospector.
            Task: Identify NOVEL, CONTROVERSIAL, or EMERGING targets for "${target}".
            ${baseContext}
            Instructions:
            1. STRICTLY NEW items (Not in Exclusion List).
            2. Targets from last 2-5 years literature.
            3. Focus on Orphan GPCRs, lncRNAs, Ion Channels.
            4. Output: Comma-separated list of 30 items.
          `;
          break;
      }

      // Usando o modelo Flash Experimental (Rápido e Gratuito)
      const response = await this.ai.models.generateContent({
        model: 'gemini-2.0-flash-exp', 
        contents: prompt,
        config: {
          temperature: temperature,
        }
      });

      const text = response.text || "";
      
      // Limpeza robusta
      const cleanText = text
        .replace(/^Here.*:/i, "")
        .replace(/Output:/i, "")
        .replace(/[\[\]*_`]/g, "")
        .trim();
      
      const terms = cleanText.split(/,|\n/);
      
      // Filtragem final
      const lowerCaseCurrentSet = new Set(currentList.map(t => t.toLowerCase().trim()));
      
      return terms
        .map(t => t.trim())
        .map(t => t.replace(/^\d+[\).]\s*/, "")) // Remove "1. "
        .filter(t => t.length > 2 && t.length < 50)
        .filter(t => !t.toLowerCase().includes("context"))
        .filter(t => !lowerCaseCurrentSet.has(t.toLowerCase()));
        
    } catch (error) {
      console.error("Gemini Mining Error:", error);
      return currentList; // Retorna lista atual em caso de erro
    }
  }

  async analyzePaper(title: string, abstract: string): Promise<string> {
    if (!this.ai) return "API Key Required.";
    
    // ... (Mantive sua lógica de análise que estava boa)
    let optimizedAbstract = abstract;
    if (abstract.length > 800) {
      optimizedAbstract = abstract.substring(0, 400) + "\n...\n" + abstract.substring(abstract.length - 400);
    }

    try {
      const response = await this.ai.models.generateContent({
        model: 'gemini-2.0-flash-exp', // Atualizado para o modelo novo
        contents: `Analyze: Title: "${title}" Abstract: "${optimizedAbstract}". Task: Identify Target, Drug, and Effect. Output format: "M: [Molecule] | E: [Effect]"`,
      });
      return response.text || "Analysis failed.";
    } catch (error) {
      return "Error analyzing paper.";
    }
  }
}
