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

  // --- FUNÇÃO DE DIAGNÓSTICO (VERSÃO SEGURA PARA BUILD) ---
  private async generateWithFallback(prompt: string): Promise<string> {
    const models = [
      'gemini-2.0-flash',       
      'gemini-2.0-flash-lite',  
      'gemini-flash-latest'     
    ];
    
    let fullErrorLog = ""; 

    for (const model of models) {
      try {
        console.log(`🔍 DIAGNÓSTICO: Tentando modelo ${model}...`);
        
        if (!this.ai) throw new Error("API Key is missing internally");

        const response = await this.ai.models.generateContent({
          model: model, 
          contents: prompt,
          config: {
            safetySettings: [
              { category: 'HARM_CATEGORY_HARASSMENT', threshold: 'BLOCK_NONE' },
              { category: 'HARM_CATEGORY_HATE_SPEECH', threshold: 'BLOCK_NONE' },
              { category: 'HARM_CATEGORY_SEXUALLY_EXPLICIT', threshold: 'BLOCK_NONE' },
              { category: 'HARM_CATEGORY_DANGEROUS_CONTENT', threshold: 'BLOCK_NONE' },
              { category: 'HARM_CATEGORY_CIVIC_INTEGRITY', threshold: 'BLOCK_NONE' }
            ]
          }
        });

        if (response.text) {
          console.log(`✅ SUCESSO com ${model}`);
          return response.text;
        }
        
      } catch (unknownError) {
        // TRUQUE PARA EVITAR ERRO DE BUILD: Cast manual
        const e = unknownError as any;
        const errorMsg = e.message || String(e);
        
        console.warn(`❌ Falha no ${model}: ${errorMsg}`);
        fullErrorLog += `\n[${model}]: ${errorMsg}`;
      }
    }
    
    throw new Error(fullErrorLog);
  }

  async mineNovelTargets(target: string, context: string, currentList: string[], strategy: MiningStrategy): Promise<string[]> {
    if (!this.ai) {
        alert("⚠️ API Key não configurada!");
        return currentList;
    }

    try {
      const listString = currentList.length > 0 ? currentList.join(", ") : "None";
      const baseContext = `- Target: "${target}"\n- List: ${listString}`;
      let prompt = "";

      switch (strategy) {
        case 'conservative':
          prompt = `Role: Data Curator. Clean this list: ${baseContext}. Output: Comma-separated strings.`;
          break;
        case 'repurposing':
          prompt = `Role: Pharmacologist. Suggest drugs for "${target}". ${baseContext}. Output: Comma-separated strings.`;
          break;
        case 'mechanism':
          prompt = `Role: Biologist. Identify pathways for "${target}". ${baseContext}. Output: Comma-separated strings.`;
          break;
        case 'blue_ocean':
          prompt = `Role: Scientist. Identify novel targets for "${target}". ${baseContext}. Output: Comma-separated strings.`;
          break;
      }

      const text = await this.generateWithFallback(prompt);
      
      const cleanText = text.replace(/Output:|Here is the list:|\[|\]|\*|- /g, "");
      const terms = cleanText.split(/,|\n/);
      
      return terms
        .map(t => t.trim())
        .map(t => t.replace(/^\d+\.\s*/, ""))
        .filter(t => t.length > 2 && !(t.includes(" ") && t.length > 50));
        
    } catch (unknownError) {
      // TRUQUE PARA EVITAR ERRO DE BUILD
      const error = unknownError as any;
      console.error("DIAGNÓSTICO FINAL:", error);
      
      alert(`🚨 RELATÓRIO DE ERRO 🚨\n\n${error.message || String(error)}`);
      
      return currentList;
    }
  }

  async analyzePaper(title: string, abstract: string): Promise<string> {
    if (!this.ai) return "API Key Required.";
    try {
      const response = await this.ai.models.generateContent({
        model: 'gemini-flash-latest', 
        contents: `Analyze: ${title}. ${abstract}`,
      });
      return response.text || "Analysis failed.";
    } catch (error) {
      return "Error analyzing paper.";
    }
  }
}
