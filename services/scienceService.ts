import { AnalysisResult, Article } from '../types';

// Lista de termos para ignorar (Blacklist interna de segurança)
const BLACKLIST_MONSTRO = new Set([
  "THE", "AND", "FOR", "NOT", "BUT", "VIA", "ALL", "WITH", "FROM", "AFTER",
  "STUDY", "ANALYSIS", "DATA", "RESULTS", "CONCLUSION", "ABSTRACT", "INTRODUCTION"
]);

// Constantes de Simulação
const N_PUBMED = 36000000;

// Helper de delay para não travar a UI
const delay = (ms: number) => new Promise(resolve => setTimeout(resolve, ms));

/**
 * Limpa e padroniza a lista de entrada
 */
export const parseSimpleUserList = (text: string): string[] => {
  if (!text) return [];
  return text.split(/[\n,;|]+/)
    .map(item => item.replace(/['"]/g, '').trim())
    .filter(item => item.length > 1) // Remove itens de 1 letra
    .filter(item => !BLACKLIST_MONSTRO.has(item.toUpperCase()));
};

/**
 * SIMULAÇÃO DO FISHER EXACT TEST (Blindada contra erros)
 * Gera estatísticas plausíveis baseadas no tamanho das palavras e aleatoriedade.
 */
export const performAnalysis = async (
  terms: string[],
  target: string,
  email: string,
  onProgress: (progress: number) => void
): Promise<AnalysisResult[]> => {
  
  const results: AnalysisResult[] = [];
  
  // Proteção contra lista vazia ou nula
  if (!terms || terms.length === 0) return [];

  // Simula "Hits" totais para a doença (baseado no nome para ser determinístico)
  let n_total_target = (target.length * 15432) % 500000; 
  if (n_total_target < 1000) n_total_target = 5000; // Mínimo de segurança

  for (let i = 0; i < terms.length; i++) {
    const term = terms[i];
    
    // Atualiza barra de progresso
    onProgress(Math.round(((i + 1) / terms.length) * 100));
    
    // Delay pequeno para a interface respirar (evita travamento visual)
    await delay(20); 

    try {
      // 1. Simula Global Count (Hits totais na ciência)
      const n_global = Math.floor((term.length * 98765) % 200000) + 50;
      
      // 2. Simula Specific Count (Interseção Alvo + Termo)
      // Lógica: 20% de chance de ser "Blue Ocean" (raro), senão proporcional
      let n_specific = 0;
      if (Math.random() > 0.85) {
         n_specific = Math.floor(Math.random() * 5); // 0 a 4 hits (Raro)
      } else {
         const baseRate = (n_global * n_total_target) / N_PUBMED;
         n_specific = Math.floor(baseRate * (Math.random() * 8)); // Simula enriquecimento
      }

      // 3. Calcula Enrichment (Ratio)
      const expected = (n_global * n_total_target) / N_PUBMED;
      const safeExpected = expected > 0 ? expected : 0.00001;
      const enrichment = (n_specific + 0.1) / safeExpected;

      // 4. Determina Status e Score
      let status: AnalysisResult['status'] = 'Neutral';
      let sortScore = 0;

      if (n_specific < 5) {
        if (n_global > 200) { 
            status = 'Blue Ocean'; 
            sortScore = 1000; // Prioridade máxima
        } else { 
            status = 'Ghost'; 
            sortScore = 0; 
        }
      } else if (n_specific <= 20) {
        if (n_global > 100) { status = 'Embryonic'; sortScore = 500; }
        else { status = 'Neutral'; sortScore = 20; }
      } else {
        if (enrichment > 5) { status = 'Gold'; sortScore = 100; }
        else if (enrichment > 1.5) { status = 'Trending'; sortScore = 200; }
        else { status = 'Saturated'; sortScore = 10; }
      }

      results.push({
        molecule: term,
        status,
        ratio: parseFloat(enrichment.toFixed(1)),
        pValue: (Math.random() * 0.05).toFixed(4), // P-valor simulado < 0.05
        targetArticles: n_specific,
        globalArticles: n_global,
        sortScore
      });

    } catch (innerError) {
      console.warn(`Erro ao calcular termo "${term}":`, innerError);
      // Não trava o loop, apenas ignora esse termo
    }
  }

  // Ordena por relevância (Score -> Ratio)
  return results.sort((a, b) => {
    if (a.sortScore !== b.sortScore) return b.sortScore - a.sortScore;
    return b.ratio - a.ratio;
  });
};

/**
 * Simula busca de artigos para leitura (Mock)
 */
export const fetchArticles = async (molecule: string, target: string): Promise<Article[]> => {
  await delay(500); // Simula rede
  
  const templates = [
    { type: 'agonist', verb: 'activates', outcome: 'increased activity' },
    { type: 'antagonist', verb: 'blocks', outcome: 'reduced inflammation' },
    { type: 'inhibitor', verb: 'inhibits', outcome: 'decreased expression' },
    { type: 'modulator', verb: 'modulates', outcome: 'altered signaling' },
    { type: 'expression', verb: 'is upregulated in', outcome: 'disease progression' }
  ];

  return Array.from({ length: 5 }).map((_, i) => {
    const tpl = templates[i % templates.length];
    return {
      id: `pmid-${Math.random().toString(36).substr(2, 9)}`,
      title: `Pharmacological evaluation of ${molecule} in ${target}: ${tpl.type} study ${i+1}`,
      journal: ['Nature', 'Cell', 'Science', 'PLOS One', 'BJP'][i % 5],
      year: 2024 - i,
      link: '#',
      abstract: `This study investigates the role of ${molecule}. We demonstrate that ${molecule} ${tpl.verb} the ${target} pathway, resulting in ${tpl.outcome}. These findings suggest ${molecule} as a potential therapeutic target.`
    };
  });
};

// Placeholder para compatibilidade caso o código tente chamar extractEntitiesRegex
export const extractEntitiesRegex = (text: string): string[] => {
    return parseSimpleUserList(text);
};
