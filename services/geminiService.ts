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
      // --- LÓGICA DE PARTIDA (COLD START) ---
      const isListEmpty = currentList.length === 0;
      let contextInstruction = "";

      if (isListEmpty) {
        // Se a lista estiver vazia, pedimos para CRIAR a base
        contextInstruction = `
        STATE: COLD START (Zero knowledge).
        TASK: Generate the INITIAL foundation list of targets for this disease.
        NOTE: Do not ask for examples, just generate the most scientifically relevant targets.
        `;
      } else {
        // Se tiver itens, pedimos para EXPANDIR
        contextInstruction = `
        STATE: EXPANSION MODE.
        EXCLUSION LIST (IGNORE these, find NEW ones):
        [${currentList.join(", ")}]
        `;
      }
      
      const baseContext = `
      TARGET ORGAN/DISEASE: "${target}"
      BIOLOGICAL CONTEXT: "${context}"
      ${contextInstruction}
      `;

      let prompt = "";
      // Temperatura calibrada para o 2.0 Flash
      let temperature = 0.7; 

      switch (strategy) {
        case 'conservative':
          temperature = 0.1;
          prompt = `
            Role: Scientific Data Curator.
            Task: Standardize pharmacological targets.
            Input: ${baseContext}
            Instructions: 
            1. If Cold Start: List standard pharmacological targets.
            2. If Expansion: Clean and standardize the input.
            3. Return a JSON Array of strings.
          `;
          break;

        case 'repurposing':
          prompt = `
            Role: Senior Pharmacologist.
            Task: Suggest approved drugs or clinical candidates for "${target}".
            ${baseContext}
            
            Instructions:
            1. Focus on FDA approved drugs and Phase 2/3 candidates.
            2. High translatability candidates.
            3. Return a JSON Array of strings (Max 25 items).
          `;
          break;

        case 'mechanism':
          prompt = `
            Role: Systems Biologist.
            Task: Identify upstream/downstream signaling pathways for "${target}".
            ${baseContext}

            Instructions:
            1. Focus on Kinases, Transcription Factors, Ion Channels.
            2. Specific subunits (e.g. "G alpha s").
            3. Return a JSON Array of strings (Max 25 items).
          `;
          break;

        case 'blue_ocean':
          temperature = 0.9;
          prompt = `
            Role: Elite Scientific Prospector.
            Task: Identify NOVEL, CONTROVERSIAL, or EMERGING targets for "${target}".
            ${baseContext}

            Instructions:
            1. Targets from literature in the last 2-5 years.
            2. Focus on Orphan GPCRs, lncRNAs, Ion Channels.
            3. Return a JSON Array of strings (Max 25 items).
          `;
          break;
      }

      // --- USO DO GEMINI 2.0 FLASH EXPERIMENTAL ---
      // Este modelo está na sua lista "Check" ✅ e é muito poderoso.
      const response = await this.ai.models.generateContent({
        model: 'gemini-2.0-flash-exp', 
        contents: prompt + "\nOutput strictly a JSON Array of strings: [\"Item1\", \"Item2\"]",
        config: {
          temperature: temperature,
          maxOutputTokens: 2000,
          responseMimeType: "application/json",
          // Desativa bloqueios para permitir termos médicos
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
      console.log(`[Gemini 2.0] Resposta:`, text);
      
      let terms: string[] = [];
      try {
        terms = JSON.parse(text);
      } catch (jsonError) {
        console.error("Erro JSON:", jsonError);
        terms = text.replace(/[\[\]"]/g, "").split(",");
      }
      
      const lowerCaseCurrentSet = new Set(currentList.map(t => t.toLowerCase().trim()));
      
      return terms
        .map(t => t.trim())
        .filter(t => t.length > 2 && t.length < 80)
        .filter(t => !lowerCaseCurrentSet.has(t.toLowerCase()));
        
    } catch (error) {
      console.error("Gemini Critical Error:", error);
      alert("A IA falhou. Detalhes no console. Tente trocar para 'gemini-1.5-flash' se o '2.0-flash-exp' continuar instável.");
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
