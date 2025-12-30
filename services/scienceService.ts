
import { AnalysisResult, Article } from '../types';
import { BLACKLIST_MONSTRO } from '../constants';

// Constants for simulation
const N_PUBMED = 36000000;

// Helper to simulate network delay
const delay = (ms: number) => new Promise(resolve => setTimeout(resolve, ms));

/**
 * Extracts entities from text using the Regex logic from backend.py.
 * Identifies Codes (P2X3), Drugs (Mirabegron), and Acronyms (VEGF).
 * USE THIS FOR RAW TEXT BLOCKS (Abstracts).
 */
export const extractEntitiesRegex = (text: string): string[] => {
  if (!text) return [];

  const entities: string[] = [];
  
  // A: Codes (Letters + Numbers) -> ex: P2X3, HSP90
  const regexCodes = /\b[A-Z]{2,}[0-9]+[A-Z0-9-]*\b/g;
  const matchCodes = text.match(regexCodes) || [];
  entities.push(...matchCodes);

  // B: Drugs (PascalCase + Chemical Suffixes)
  const regexDrugs = /\b[A-Z][a-z]{3,}(?:ine|in|mab|ib|ol|on|one|il|ide|ate|ase)\b/g;
  const matchDrugs = text.match(regexDrugs) || [];
  entities.push(...matchDrugs);

  // C: Acronyms (3+ Uppercase) -> ex: VEGF, mTOR
  const regexAcronyms = /\b[A-Z]{3,}\b/g;
  const matchAcronyms = text.match(regexAcronyms) || [];
  entities.push(...matchAcronyms);

  return filterEntities(entities);
};

/**
 * Parses a simple user list (CSV, Newlines, etc).
 * USE THIS FOR FILE UPLOADS.
 * It is much more permissive than the Regex extractor.
 */
export const parseSimpleUserList = (text: string): string[] => {
  if (!text) return [];

  // Split by comma, newline, semicolon, or pipe
  const rawItems = text.split(/[\n,;|]+/);
  
  const entities = rawItems.map(item => {
    // Remove quotes and extra spaces
    return item.replace(/['"]/g, '').trim();
  }).filter(item => item.length > 1); // Allow small molecules like "NO", "CO"

  return filterEntities(entities);
};

/**
 * Shared filtering logic to remove blacklist items and duplicates
 */
const filterEntities = (entities: string[]): string[] => {
  const seen = new Set<string>();
  const cleanEntities: string[] = [];

  for (let e of entities) {
    // Clean punctuation
    e = e.replace(/^[.,-;:()[\]]+|[.,-;:()[\]]+$/g, "");
    
    if (e.length < 2) continue; // Allow length 2 (e.g. NO, CO)
    if (/^\d+$/.test(e)) continue; // Skip pure numbers

    const upper = e.toUpperCase();
    if (BLACKLIST_MONSTRO.has(upper)) continue;
    
    // Explicit preposition block
    if (["WITH", "FROM", "AFTER", "DURING", "AND", "THE"].includes(upper)) continue;

    if (!seen.has(upper)) {
      seen.add(upper);
      cleanEntities.push(e);
    }
  }
  return cleanEntities;
}

/**
 * Simulates the Fisher Exact Test and Enrichment Logic.
 */
export const performAnalysis = async (
  terms: string[],
  target: string,
  email: string,
  onProgress: (progress: number) => void
): Promise<AnalysisResult[]> => {
  
  const results: AnalysisResult[] = [];
  
  // Simulate fetching total articles for the target
  let n_total_target = (target.length * 15432) % 500000; 
  if (n_total_target === 0) n_total_target = 1000;

  for (let i = 0; i < terms.length; i++) {
    const term = terms[i];
    await delay(30); // Reduced delay slightly for speed

    // Simulate Global Count (Total papers for the molecule)
    const n_global = (term.length * 98765) % 200000;
    
    // Simulate Specific Count with a bias towards "Blue Ocean" (low hits)
    // 70% chance of being extremely low hit (Blue Ocean candidate) - HIGH PROBABILITY
    let n_specific = 0;
    if (Math.random() > 0.30) { // Threshold lowered to 0.30 -> 70% chance
        n_specific = Math.floor(Math.random() * 6); // 0 to 5 hits
    } else {
        n_specific = Math.floor((n_global * n_total_target) / N_PUBMED * (Math.random() * 10));
    }

    // Calculate Enrichment
    const expected = (n_global * n_total_target) / N_PUBMED;
    const enrichment = (n_specific + 0.1) / (expected > 0 ? expected : 0.00001);

    // Determine Status Tag
    let status: AnalysisResult['status'] = 'Neutral';
    let sortScore = 0;

    // Relaxed Blue Ocean Logic: Less than 10 hits
    if (n_specific < 10) {
      if (n_global > 80) { 
          status = 'Blue Ocean'; 
          sortScore = 1500; // Super Boost for Blue Ocean
      } else { 
          status = 'Ghost'; 
          sortScore = 5; 
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
      pValue: (Math.random() * 0.05).toFixed(4), // Simulated p-value
      targetArticles: n_specific,
      globalArticles: n_global,
      sortScore
    });

    onProgress(((i + 1) / terms.length) * 100);
  }

  // Sort by Importance (SortScore desc, Ratio desc)
  return results.sort((a, b) => {
    if (a.sortScore !== b.sortScore) return b.sortScore - a.sortScore;
    return b.ratio - a.ratio;
  });
};

export const fetchArticles = async (molecule: string, target: string): Promise<Article[]> => {
  await delay(1000);
  
  // Variations to prevent "same tags" AI output
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
