export interface AnalysisResult {
  molecule: string;
  status: 'Gold' | 'Trending' | 'Blue Ocean' | 'Embryonic' | 'Saturated' | 'Ghost' | 'Neutral';
  ratio: number;
  pValue: string;
  targetArticles: number;
  globalArticles: number;
  sortScore: number;
}

export interface Article {
  id: string;
  title: string;
  journal: string;
  year: number;
  link: string;
  abstract: string; // Simulated abstract/context
}

export interface NewsItem {
  id: number;
  title: string;
  source: string;
  link: string;
  image: string;
}

export type MiningStrategy = 'conservative' | 'repurposing' | 'mechanism' | 'blue_ocean';

export interface AppState {
  page: 'home' | 'results';
  email: string;
  target: string;
  context: string;
  yearStart: number;
  yearEnd: number;
  targetList: string; // Comma separated
  results: AnalysisResult[];
  apiKey: string;
  useAI: boolean;
  miningStrategy: MiningStrategy;
  selectedArticleDetails: Article[] | null;
}

export enum Language {
  PT = 'pt',
  EN = 'en'
}