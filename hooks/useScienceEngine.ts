import { useState, useEffect } from 'react';
import { AppState, AnalysisResult, Article } from '../types';
import { GeminiService } from '../services/geminiService';
import { performAnalysis, extractEntitiesRegex, fetchArticles, parseSimpleUserList } from '../services/scienceService';
import { PRESETS_FRONTEIRA } from '../constants';

export const useScienceEngine = () => {
  const [state, setState] = useState<AppState>({
    page: 'home',
    email: '',
    target: '',
    context: '',
    yearStart: 2015,
    yearEnd: 2025,
    targetList: '',
    results: [],
    apiKey: '',
    useAI: true,
    miningStrategy: 'blue_ocean', // Default
    selectedArticleDetails: null
  });

  const [isLoading, setIsLoading] = useState(false);
  const [progress, setProgress] = useState(0);
  const [gemini, setGemini] = useState<GeminiService>(new GeminiService(''));

  useEffect(() => {
    if (state.apiKey) {
      gemini.updateKey(state.apiKey);
    }
  }, [state.apiKey]);

  const updateState = (updates: Partial<AppState>) => {
    setState(prev => ({ ...prev, ...updates }));
  };

  const handlePresetAdd = (presetName: string) => {
    const items = PRESETS_FRONTEIRA[presetName];
    const current = state.targetList ? state.targetList.split(',').map(s => s.trim()) : [];
    const novel = items.filter(i => !current.includes(i));
    updateState({ targetList: [...current, ...novel].join(', ') });
  };

  const handleAddAllPresets = () => {
    const allItems = Object.values(PRESETS_FRONTEIRA).flat();
    const current = state.targetList ? state.targetList.split(',').map(s => s.trim()) : [];
    const novel = allItems.filter(i => !current.includes(i));
    updateState({ targetList: [...current, ...novel].join(', ') });
  };

  const handleFileUpload = (text: string) => {
    // Use the simple parser for uploads, respecting user's exact keywords
    const extracted = parseSimpleUserList(text);
    
    const current = state.targetList ? state.targetList.split(',').map(s => s.trim()) : [];
    const unique = [...new Set([...current, ...extracted])];
    
    updateState({ targetList: unique.join(', ') });
    return extracted.length;
  };

  const handleDeepMine = async () => {
    if (!state.target) throw new Error("Target is required");
    
    setIsLoading(true);
    let currentTerms = state.targetList 
      ? state.targetList.split(',').map(s => s.trim()).filter(x => x) 
      : [];

    // 1. Initial Seeding if empty
    if (currentTerms.length === 0) {
      const randomPreset = Object.values(PRESETS_FRONTEIRA)[0];
      currentTerms = [...randomPreset];
    }

    // 2. Cross-Tissue Intelligence
    if (state.useAI && state.apiKey) {
      try {
        const resultTerms = await gemini.mineNovelTargets(
          state.target,
          state.context,
          currentTerms,
          state.miningStrategy
        );
        
        if (state.miningStrategy === 'conservative') {
           // Conservative is a filter/cleaner, so we REPLACE the list with the cleaned version
           currentTerms = resultTerms;
        } else {
           // Repurposing, Mechanism, Blue Ocean: We ADD the new findings to the existing list
           // Ensure uniqueness case-insensitive
           const existingSet = new Set(currentTerms.map(t => t.toUpperCase()));
           const newItems = resultTerms.filter(t => !existingSet.has(t.toUpperCase()));
           currentTerms = [...currentTerms, ...newItems];
        }

      } catch (e) {
        console.warn("Deep Mining failed, falling back to basic list");
      }
    } else if (state.context) {
      // Fallback: simple regex on context if no AI
      const contextTerms = extractEntitiesRegex(state.context);
      const existingSet = new Set(currentTerms.map(t => t.toUpperCase()));
      const newItems = contextTerms.filter(t => !existingSet.has(t.toUpperCase()));
      currentTerms = [...currentTerms, ...newItems];
    }
    
    updateState({ targetList: currentTerms.join(', ') });
    setIsLoading(false);
  };

  const executeAnalysis = async () => {
    if (!state.target || !state.targetList) return;
    setIsLoading(true);
    setProgress(0);
    
    const terms = state.targetList.split(',').map(s => s.trim()).filter(s => s.length > 0);
    const results = await performAnalysis(terms, state.target, state.email, (p) => setProgress(p));
    
    updateState({ results, page: 'results' });
    setIsLoading(false);
  };

  const fetchDetails = async (molecule: string) => {
    const articles = await fetchArticles(molecule, state.target);
    updateState({ selectedArticleDetails: articles });
  };

  const analyzeArticle = async (title: string, abstract: string) => {
    if (!state.apiKey) return "API Key Required";
    return await gemini.analyzePaper(title, abstract);
  };

  return {
    state,
    updateState,
    isLoading,
    progress,
    handlePresetAdd,
    handleAddAllPresets,
    handleFileUpload,
    handleDeepMine,
    executeAnalysis,
    fetchDetails,
    analyzeArticle
  };
};