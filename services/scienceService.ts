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

  // --- FUNÇÃO DE DIAGNÓSTICO ---
  private async generateWithFallback(prompt: string): Promise<string> {
    // Vamos testar os modelos principais. Se um falhar, guardamos o erro para te mostrar.
    const models = [
      'gemini-2.0-flash',       // Tentativa 1
      'gemini-2.0-flash-lite',  // Tentativa 2
      'gemini-flash-latest'     // Tentativa 3
    ];
    
    let fullErrorLog = ""; // Vai acumular os erros de cada tentativa

    for (const model of models) {
      try {
        console.log(`🔍 DIAGNÓSTICO: Tentando modelo ${model}...`);
        
        const response = await this.ai!.models.generateContent({
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
        
      } catch (e: any) {
        // Captura o erro bruto
        const errorMsg = e.message || JSON.stringify(e);
        const status = e.status || "Sem Status";
        
        console.warn(`❌ Falha no ${model}: ${errorMsg}`);
        
        // Adiciona ao relatório de erros
        fullErrorLog += `\n\n🔴 Modelo: ${model}\nStatus: ${status}\nErro: ${errorMsg}`;
      }
    }
    
    // Se chegou aqui, todos falharam. Lança o relatório completo.
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

      // Prompt simplificado para teste de conexão
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

      // Chama a função de diagnóstico
      const text = await this.generateWithFallback(prompt);
      
      const cleanText = text.replace(/Output:|Here is the list:|\[|\]|\*|- /g, "");
      const terms = cleanText.split(/,|\n/);
      
      return terms
        .map(t => t.trim())
        .map(t => t.replace(/^\d+\.\s*/, ""))
        .filter(t => t.length > 2 && !(t.includes(" ") && t.length > 50));
        
    } catch (error: any) {
      console.error("DIAGNÓSTICO FINAL:", error);
      
      // --- AQUI ESTÁ O QUE VOCÊ PEDIU ---
      // Um alerta gigante com o erro exato para você copiar/printar
      alert(`🚨 RELATÓRIO DE ERRO DA IA 🚨\nCopie isso e mande para o suporte:\n${error.message}`);
      
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
