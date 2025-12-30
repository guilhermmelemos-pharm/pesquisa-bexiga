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

  // --- FUNÇÃO DE BACKUP COM ALERTA DE ERRO ---
  private async generateWithBackup(prompt: string): Promise<string> {
    // Lista de Prioridade: Tenta o 3 (que você quer), se falhar, vai pro 2 (que funciona)
    const models = ['gemini-3-pro-preview', 'gemini-2.0-flash-exp']; 
    
    for (const model of models) {
      try {
        console.log(`🤖 Tentando modelo: ${model}...`);
        
        const response = await this.ai!.models.generateContent({
          model: model, 
          contents: prompt,
          // Travas desligadas para permitir Farmacologia
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
          console.log(`✅ Sucesso com ${model}!`);
          return response.text;
        }
        
      } catch (e: any) {
        const errorMsg = e.message || JSON.stringify(e);
        console.warn(`❌ Falha no ${model}:`, errorMsg);

        // SE O GEMINI 3 FALHAR, AVISAMOS VOCÊ NA TELA ANTES DE TROCAR
        if (model.includes('3-pro')) {
            alert(`⚠️ O Gemini 3 falhou.\nErro: ${errorMsg}\n\n🔄 Trocando automaticamente para o Gemini 2.0...`);
        }
      }
    }
    throw new Error("Todos os modelos falharam. Verifique sua API Key.");
  }

  /**
   * Deep Mining Strategy (Sua versão Clássica + Backup + Alerta)
   */
  async mineNovelTargets(target: string, context: string, currentList: string[], strategy: MiningStrategy): Promise<string[]> {
    if (!this.ai) {
        alert("⚠️ API Key não configurada!");
        return currentList;
    }

    try {
      const isComparative = (context || "").length > 2;
      // Se lista vazia, passa "None" para o prompt entender que é o início
      const listString = currentList.length > 0 ? currentList.join(", ") : "None";
      
      const baseContext = `- Primary Target Organ/Disease: "${target}"\n- Current Candidate List: ${listString}`;
      let prompt = "";

      // SEUS PROMPTS ORIGINAIS MANTIDOS
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

      // Chama a função blindada
      const text = await this.generateWithBackup(prompt);
      
      console.log("Texto Bruto da IA:", text);

      // SUA LÓGICA DE LIMPEZA ORIGINAL (REGEX)
      const cleanText = text.replace(/Output:|Here is the list:|\[|\]|\*|- /g, "");
      const terms = cleanText.split(/,|\n/);
      
      return terms
        .map(t => t.trim())
        .map(t => t.replace(/^\d+\.\s*/, ""))
        // Filtro levemente ajustado
        .filter(t => t.length > 2 && !(t.includes(" ") && t.length > 50));
        
    } catch (error) {
      console.error("Erro Geral:", error);
      alert("A IA não conseguiu responder com nenhum modelo. Verifique se sua conta tem cota disponível.");
      return currentList;
    }
  }

  async analyzePaper(title: string, abstract: string): Promise<string> {
    if (!this.ai) return "API Key Required.";
    
    if (!abstract || typeof abstract !== 'string') return "No abstract.";

    let optimizedAbstract = abstract;
    if (abstract.length > 600) {
      optimizedAbstract = abstract.substring(0, 300) + "\n...\n" + abstract.substring(abstract.length - 300);
    }

    try {
      // Tenta o flash para análise rápida
      const response = await this.ai.models.generateContent({
        model: 'gemini-1.5-flash', 
        contents: `Analyze: Title: "${title}" Abstract: "${optimizedAbstract}". Task: Identify Target, Drug, and Effect. Output format: "M: [Molecule] | E: [Effect]"`,
      });
      return response.text || "Analysis failed.";
    } catch (error) {
      return "Error analyzing paper.";
    }
  }
}
