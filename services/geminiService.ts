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

  async mineNovelTargets(target: string, context: string, currentList: string[], strategy: MiningStrategy): Promise<string[]> {
    if (!this.ai) {
        console.warn("API Key missing");
        return currentList;
    }

    try {
      const currentListStr = currentList.length > 0 ? currentList.join(", ") : "None";
      
      const baseContext = `
      CONTEXT:
      - Primary Disease/Organ: "${target}"
      - Bio Context: "${context}"
      
      EXCLUSION LIST (DO NOT OUTPUT THESE):
      [${currentListStr}]
      `;

      let prompt = "";
      // Temperatura alta para o Gemini 3 ser criativo na busca
      let temperature = 0.8; 

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
          prompt = `
            Role: Senior Pharmacologist.
            Task: Suggest approved drugs or clinical candidates for "${target}".
            ${baseContext}
            Instructions:
            1. STRICTLY NEW items (Not in Exclusion List).
            2. Focus on FDA approved drugs and Phase 2/3 candidates.
            3. Output: Comma-separated list of 30 specific items.
          `;
          break;

        case 'mechanism':
          prompt = `
            Role: Systems Biologist.
            Task: Identify upstream/downstream signaling pathways for "${target}".
            ${baseContext}
            Instructions:
            1. STRICTLY NEW items.
            2. Focus on Kinases, Transcription Factors, Ion Channels, Enzymes.
            3. Output: Comma-separated list of 30 specific targets.
          `;
          break;

        case 'blue_ocean':
          temperature = 1.0; // Máxima criatividade
          prompt = `
            Role: Elite Scientific Prospector.
            Task: Identify NOVEL, CONTROVERSIAL, or EMERGING targets for "${target}".
            ${baseContext}
            Instructions:
            1. STRICTLY NEW items (Not in Exclusion List).
            2. Targets from literature in the last 2-5 years.
            3. Focus on Orphan GPCRs, lncRNAs, Ion Channels, and metabolic sensors.
            4. Output: Comma-separated list of 30 items.
          `;
          break;
      }

      // --- CONFIGURAÇÃO PARA O GEMINI 3 ---
      const response = await this.ai.models.generateContent({
        // Usando o modelo mais forte da sua lista
        model: 'gemini-3-pro-preview', 
        contents: prompt,
        config: {
          temperature: temperature,
          // CRÍTICO: Desativa o bloqueio "Médico/Perigoso" para permitir farmacologia
          safetySettings: [
            { category: 'HARM_CATEGORY_HARASSMENT', threshold: 'BLOCK_NONE' },
            { category: 'HARM_CATEGORY_HATE_SPEECH', threshold: 'BLOCK_NONE' },
            { category: 'HARM_CATEGORY_SEXUALLY_EXPLICIT', threshold: 'BLOCK_NONE' },
            { category: 'HARM_CATEGORY_DANGEROUS_CONTENT', threshold: 'BLOCK_NONE' },
            { category: 'HARM_CATEGORY_CIVIC_INTEGRITY', threshold: 'BLOCK_NONE' }
          ]
        }
      });

      const text = response.text || "";
      console.log(`[Gemini 3 Pro] Resposta (${strategy}):`, text);
      
      // Limpeza da resposta
      const cleanText = text
        .replace(/^Here.*:/i, "")
        .replace(/Output:/i, "")
        .replace(/[\[\]*_`]/g, "")
        .trim();
      
      const terms = cleanText.split(/,|\n/);
      
      const lowerCaseCurrentSet = new Set(currentList.map(t => t.toLowerCase().trim()));
      
      const finalTerms = terms
        .map(t => t.trim())
        .map(t => t.replace(/^\d+[\).]\s*/, "")) 
        .filter(t => t.length > 2 && t.length < 60)
        .filter(t => !t.toLowerCase().includes("context"))
        .filter(t => !lowerCaseCurrentSet.has(t.toLowerCase()));

      return finalTerms;
        
    } catch (error) {
      console.error("Gemini 3 Error:", error);
      // Fallback para o 2.5 Flash se o 3 der erro de cota ou instabilidade
      console.log("Tentando fallback para gemini-2.5-flash...");
      try {
          return await this.retryWithFallback(target, context, currentList, strategy);
      } catch (e) {
          alert("Erro na IA: Verifique a API Key (F12 para detalhes).");
          return currentList; 
      }
    }
  }

  // Fallback usando o modelo rápido da sua lista
  private async retryWithFallback(target: string, context: string, currentList: string[], strategy: MiningStrategy): Promise<string[]> {
      const response = await this.ai!.models.generateContent({
        model: 'gemini-2.5-flash', 
        contents: `List 10 novel pharmacological targets for ${target}. Comma separated.`,
      });
      const text = response.text || "";
      return text.split(',').map(t => t.trim());
  }

  async analyzePaper(title: string, abstract: string): Promise<string> {
    if (!this.ai) return "API Key Required.";
    
    let optimizedAbstract = abstract;
    if (abstract.length > 800) {
      optimizedAbstract = abstract.substring(0, 400) + "\n...\n" + abstract.substring(abstract.length - 400);
    }

    try {
      const response = await this.ai.models.generateContent({
        // Usando o Flash Preview para leitura rápida
        model: 'gemini-3-flash-preview', 
        contents: `Analyze: Title: "${title}" Abstract: "${optimizedAbstract}". Task: Identify Target, Drug, and Effect. Output format: "M: [Molecule] | E: [Effect]"`,
      });
      return response.text || "Analysis failed.";
    } catch (error) {
      return "Error analyzing paper.";
    }
  }
}
