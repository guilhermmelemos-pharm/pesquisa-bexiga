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
        alert("⚠️ API Key não configurada! Vá nas configurações.");
        return currentList;
    }

    try {
      const currentListStr = currentList.length > 0 ? currentList.join(", ") : "None";
      
      // --- O TRUQUE DO "CONTEXTO FANTASMA" ---
      // Damos exemplos para a IA entender o formato, mesmo com lista vazia.
      const baseContext = `
      TARGET CONTEXT:
      - Primary Disease/Organ: "${target}"
      - Biological Context: "${context}"
      
      KNOWN EXCLUSION LIST (DO NOT REPEAT THESE):
      [${currentListStr}]
      `;

      let prompt = "";
      // Temperatura alta para criatividade
      let temperature = 0.8; 

      switch (strategy) {
        case 'conservative':
          temperature = 0.2;
          prompt = `
            Role: Scientific Data Curator.
            Task: Clean and Standardize the candidate list.
            Input: ${baseContext}
            Instructions: 
            1. Convert Drugs to Targets (e.g. Semaglutide -> GLP1R).
            2. Remove noise terms. 
            3. Return a clean JSON Array of strings.
          `;
          break;

        case 'repurposing':
          prompt = `
            Role: Senior Pharmacologist.
            Task: Suggest approved drugs or clinical candidates for "${target}".
            ${baseContext}
            
            EXAMPLE OUTPUT FORMAT (Follow this style):
            ["SGLT2 inhibitors", "Metformin", "Atorvastatin"]

            Instructions:
            1. STRICTLY NEW items (Not in Exclusion List).
            2. Focus on FDA approved drugs and Phase 2/3 candidates.
            3. Return a JSON Array of strings (Max 30 items).
          `;
          break;

        case 'mechanism':
          prompt = `
            Role: Systems Biologist.
            Task: Identify upstream/downstream signaling pathways for "${target}".
            ${baseContext}

            EXAMPLE OUTPUT FORMAT (Follow this style):
            ["mTORC1", "NF-kB", "P2X3 receptor", "TRPV1"]

            Instructions:
            1. STRICTLY NEW items.
            2. Focus on Kinases, Transcription Factors, Ion Channels, Enzymes.
            3. Return a JSON Array of strings (Max 30 items).
          `;
          break;

        case 'blue_ocean':
          temperature = 1.0;
          prompt = `
            Role: Elite Scientific Prospector.
            Task: Identify NOVEL, CONTROVERSIAL, or EMERGING targets for "${target}".
            ${baseContext}

            EXAMPLE OUTPUT FORMAT (Follow this style):
            ["Piezo2", "LncRNA MALAT1", "GPR183", "OSCA1"]

            Instructions:
            1. STRICTLY NEW items (Not in Exclusion List).
            2. Targets from literature in the last 2-5 years.
            3. Focus on Orphan GPCRs, lncRNAs, Ion Channels.
            4. Return a JSON Array of strings (Max 30 items).
          `;
          break;
      }

      const response = await this.ai.models.generateContent({
        model: 'gemini-2.0-flash-exp', 
        contents: prompt + "\nOutput strictly a JSON Array of strings: [\"Item1\", \"Item2\"]",
        config: {
          temperature: temperature,
          responseMimeType: "application/json",
          safetySettings: [
            { category: 'HARM_CATEGORY_HARASSMENT', threshold: 'BLOCK_NONE' },
            { category: 'HARM_CATEGORY_HATE_SPEECH', threshold: 'BLOCK_NONE' },
            { category: 'HARM_CATEGORY_SEXUALLY_EXPLICIT', threshold: 'BLOCK_NONE' },
            { category: 'HARM_CATEGORY_DANGEROUS_CONTENT', threshold: 'BLOCK_NONE' },
            { category: 'HARM_CATEGORY_CIVIC_INTEGRITY', threshold: 'BLOCK_NONE' }
          ]
        }
      });

      const text = response.text || "[]";
      console.log(`[Gemini JSON] Resposta:`, text);
      
      let terms: string[] = [];
      try {
        terms = JSON.parse(text);
      } catch (jsonError) {
        console.error("Erro ao ler JSON da IA:", jsonError);
        terms = text.replace(/[\[\]"]/g, "").split(",");
      }
      
      // Limpeza final
      const lowerCaseCurrentSet = new Set(currentList.map(t => t.toLowerCase().trim()));
      
      return terms
        .map(t => t.trim())
        .filter(t => t.length > 2 && t.length < 80)
        .filter(t => !lowerCaseCurrentSet.has(t.toLowerCase()));
        
    } catch (error) {
      console.error("Gemini Critical Error:", error);
      alert("A IA falhou. Verifique se sua API Key tem acesso ao modelo 'gemini-2.0-flash-exp'. Detalhes no Console.");
      return currentList; 
    }
  }

  async analyzePaper(title: string, abstract: string): Promise<string> {
    if (!this.ai) return "API Key Required.";
    
    let optimizedAbstract = abstract;
    if (abstract.length > 800) {
      optimizedAbstract = abstract.substring(0, 400) + "\n...\n" + abstract.substring(abstract.length - 400);
    }

    try {
      const response = await this.ai.models.generateContent({
        model: 'gemini-2.0-flash-exp', 
        contents: `Analyze: Title: "${title}" Abstract: "${optimizedAbstract}". Task: Identify Target, Drug, and Effect. Output format: "M: [Molecule] | E: [Effect]"`,
      });
      return response.text || "Analysis failed.";
    } catch (error) {
      return "Error analyzing paper.";
    }
  }
}
