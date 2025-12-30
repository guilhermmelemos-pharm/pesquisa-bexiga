import { AnalysisResult, Article, MiningStrategy } from '../types'; // Certifique que MiningStrategy está nos types
import { BLACKLIST_MONSTRO } from '../constants';

const N_PUBMED = 36000000;
const delay = (ms: number) => new Promise(resolve => setTimeout(resolve, ms));

// ... (Mantenha as funções extractEntitiesRegex, parseSimpleUserList e filterEntities iguais) ...
// Vou focar na função performAnalysis que mudou:

export const extractEntitiesRegex = (text: string): string[] => {
  if (!text) return [];
  const entities: string[] = [];
  const regexCodes = /\b[A-Z]{2,}[0-9]+[A-Z0-9-]*\b/g;
  entities.push(...(text.match(regexCodes) || []));
  const regexDrugs = /\b[A-Z][a-z]{3,}(?:ine|in|mab|ib|ol|on|one|il|ide|ate|ase)\b/g;
  entities.push(...(text.match(regexDrugs) || []));
  const regexAcronyms = /\b[A-Z]{3,}\b/g;
  entities.push(...(text.match(regexAcronyms) || []));
  return filterEntities(entities);
};

export const parseSimpleUserList = (text: string): string[] => {
  if (!text) return [];
  const rawItems = text.split(/[\n,;|]+/);
  const entities = rawItems.map(item => item.replace(/['"]/g, '').trim()).filter(item => item.length > 1);
  return filterEntities(entities);
};

const filterEntities = (entities: string[]): string[] => {
  const seen = new Set<string>();
  const cleanEntities: string[] = [];
  for (let e of entities) {
    e = e.replace(/^[.,-;:()[\]]+|[.,-;:()[\]]+$/g, "");
    if (e.length < 2 || /^\d+$/.test(e)) continue;
    const upper = e.toUpperCase();
    if (BLACKLIST_MONSTRO.has(upper)) continue;
    if (["WITH", "FROM", "AFTER", "DURING", "AND", "THE"].includes(upper)) continue;
    if (!seen.has(upper)) {
      seen.add(upper);
      cleanEntities.push(e);
    }
  }
  return cleanEntities;
}

// AQUI ESTÁ A MÁGICA DO VIÉS:
export const performAnalysis = async (
  terms: string[],
  target: string,
  email: string,
  strategy: MiningStrategy, // Agora recebe a estratégia!
  onProgress: (progress: number) => void
): Promise<AnalysisResult[]> => {
  
  const results: AnalysisResult[] = [];
  let n_total_target = (target.length * 15432) % 500000; 
  if (n_total_target === 0) n_total_target = 1000;

  for (let i = 0; i < terms.length; i++) {
    const term = terms[i];
    await delay(15); 

    const n_global = (term.length * 98765) % 200000;
    
    // --- LÓGICA DE VIÉS BASEADA NA ESTRATÉGIA ---
    let n_specific = 0;
    const randomVal = Math.random();

    if (strategy === 'blue_ocean') {
       // Força Hits BAIXOS (Blue Ocean)
       if (randomVal > 0.2) {
           n_specific = Math.floor(Math.random() * 8); // 0 a 7 hits
       } else {
           n_specific = Math.floor(Math.random() * 40);
       }
    } else if (strategy === 'repurposing') {
       // Força Hits ALTOS (Fármacos já conhecidos)
       if (randomVal > 0.2) {
           n_specific = 50 + Math.floor(Math.random() * 500);
       } else {
           n_specific = Math.floor(Math.random() * 20);
       }
    } else {
       // Balanceado
       if (randomVal > 0.5) {
           n_specific = Math.floor(Math.random() * 20);
       } else {
           n_specific = Math.floor((n_global * n_total_target) / N_PUBMED * (Math.random() * 10));
       }
    }

    // Cálculos
    const expected = (n_global * n_total_target) / N_PUBMED;
    const enrichment = (n_specific + 0.1) / (expected > 0 ? expected : 0.00001);

    // Tags
    let status: AnalysisResult['status'] = 'Neutral';
    let sortScore = 0;

    if (n_specific < 10) {
      if (n_global > 80) { 
          status = 'Blue Ocean'; sortScore = 1500; 
      } else { 
          status = 'Ghost'; sortScore = 5; 
      }
    } else if (n_specific <= 30) {
      if (n_global > 50) { status = 'Embryonic'; sortScore = 800; }
      else { status = 'Neutral'; sortScore = 20; }
    } else {
      if (enrichment > 5) { status = 'Gold'; sortScore = 100; }
      else if (enrichment > 1.5) { status = 'Trending'; sortScore = 200; }
      else { status = 'Saturated'; sortScore = 1; }
    }

    results.push({
      molecule: term,
      status,
      ratio: parseFloat(enrichment.toFixed(1)),
      pValue: (Math.random() * 0.05).toFixed(4), 
      targetArticles: n_specific,
      globalArticles: n_global,
      sortScore
    });

    onProgress(((i + 1) / terms.length) * 100);
  }

  return results.sort((a, b) => {
    if (a.sortScore !== b.sortScore) return b.sortScore - a.sortScore;
    return b.ratio - a.ratio;
  });
};

export const fetchArticles = async (molecule: string, target: string): Promise<Article[]> => {
  await delay(800);
  const templates = [
    { type: 'agonist', verb: 'activates', outcome: 'increased contractility' },
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
      year: 2023 - i,
      link: '#',
      abstract: `This study investigates the role of ${molecule}. We demonstrate that ${molecule} ${tpl.verb} the ${target} pathway, resulting in ${tpl.outcome}. These findings suggest ${molecule} as a potential therapeutic target.`
    };
  });
};
