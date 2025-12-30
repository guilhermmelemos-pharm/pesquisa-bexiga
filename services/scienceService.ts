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

  // --- SISTEMA DE ROTAÇÃO DE MODELOS (APENAS OS COMPATÍVEIS) ---
  private async generateWithFallback(prompt: string): Promise<string> {
    // LISTA SEGURA: Apenas modelos que apareceram como ✅ na sua lista
    const models = [
      'gemini-2.0-flash',       // O mais rápido e moderno disponível pra você
      'gemini-2.0-flash-lite',  // Versão leve (backup)
      'gemini-flash-latest'     // O ponteiro genérico do Google (último recurso)
    ];
    
    let lastError = "";

    for (const model of models) {
      try {
        console.log(`🤖 Tentando modelo: ${model}...`);
        
        const response = await this.ai!.models.generateContent({
          model: model, 
          contents: prompt,
          config: {
            // Travas desligadas para permitir Farmacologia
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
          console.log(`✅ Sucesso com ${model}!`);
          return response.text;
        }
        
      } catch (e: any) {
        lastError = e.message || String(e);
        console.warn(`⚠️ Falha no ${model}:`, lastError);
        // Continua para o próximo modelo da lista...
      }
    }
    
    // Se chegou aqui, nenhum funcionou
    throw new Error(`Falha em todos os modelos (2.0 e Latest). Erro: ${lastError}`);
  }

  async mineNovelTargets(target: string, context: string, currentList: string[], strategy: MiningStrategy): Promise<string[]> {
    if (!this.ai) return currentList;

    try {
      const isComparative = (context || "").length > 2;
      const listString = currentList.length > 0 ? currentList.join(", ") : "None";
      
      const baseContext = `- Primary Target Organ/Disease: "${target}"\n- Current Candidate List: ${listString}`;
      let prompt = "";

      switch (strategy) {
        case 'conservative':
          prompt = `
            Role: Scientific Data Curator.
            Task: Standardize pharmacological targets.
            Input: ${baseContext}
            Instructions: 
            1. If Cold Start: List standard pharmacological targets.
            2. If Expansion: Clean and standardize the input.
            3. Return a comma-separated list.
          `;
          break;

        case 'repurposing':
          prompt = `
            Role: Drug Repurposing Specialist.
            Task: Suggest approved drugs/candidates for "${target}".
            ${baseContext}
            ${isComparative ? `- Context: "${context}"` : ''}
            Instructions:
            1. Focus on FDA approved drugs and Phase 2/3 candidates.
            2. High translatability.
            3. Output: Comma-separated list of top 20 candidates.
          `;
          break;

        case 'mechanism':
          prompt = `
            Role: Molecular Biologist.
            Task: Identify signaling pathways for "${target}".
            ${baseContext}
            Instructions:
            1. Focus on Kinases, Transcription Factors, Ion Channels.
            2. Specific subunits.
            3. Output: Comma-separated list of top 20 targets.
          `;
          break;

        case 'blue_ocean':
          prompt = `
            Role: Elite Scientific Prospector.
            Mission: Identify UNDERSTUDIED (Blue Ocean) targets for: "${target}".
            ${baseContext}
            Instructions:
            1. Targets from literature (last 5 years).
            2. Focus on Orphan GPCRs, lncRNAs, Ion Channels.
            3. Output: Comma-separated list of top 20 novel targets.
          `;
          break;
      }

      // Chama a função robusta
      const text = await this.generateWithFallback(prompt);
      
      console.log("IA Respondeu:", text);

      // Limpeza Regex (Sua preferida)
      const cleanText = text.replace(/Output:|Here is the list:|\[|\]|\*|- /g, "");
      const terms = cleanText.split(/,|\n/);
      
      return terms
        .map(t => t.trim())
        .map(t => t.replace(/^\d+\.\s*/, ""))
        .filter(t => t.length > 2 && !(t.includes(" ") && t.length > 50));
        
    } catch (error: any) {
      console.error("Erro na Mineração:", error);
      alert(`Erro na IA: ${error.message || "Verifique o console"}`);
      return currentList;
    }
  }

  async analyzePaper(title: string, abstract: string): Promise<string> {
    if (!this.ai) return "API Key Required.";
    if (!abstract) return "No abstract.";

    let optimizedAbstract = abstract;
    if (abstract.length > 600) {
      optimizedAbstract = abstract.substring(0, 300) + "\n...\n" + abstract.substring(abstract.length - 300);
    }

    try {
      // Tenta analisar com o modelo mais leve disponível
      const response = await this.ai.models.generateContent({
        model: 'gemini-flash-latest', 
        contents: `Analyze: Title: "${title}" Abstract: "${optimizedAbstract}". Task: Identify Target, Drug, and Effect. Output format: "M: [Molecule] | E: [Effect]"`,
      });
      return response.text || "Analysis failed.";
    } catch (error) {
      return "Error analyzing paper.";
    }
  }
}
