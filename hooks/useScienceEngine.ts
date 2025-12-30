
import { useState } from 'react';
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
    useAI: true,
    miningStrategy: 'blue_ocean', // Default
    selectedArticleDetails: null,
    miningStatus: ''
  });

  const [isLoading, setIsLoading] = useState(false);
  const [progress, setProgress] = useState(0);
  const [gemini] = useState<GeminiService>(() => new GeminiService());

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
    const extracted = parseSimpleUserList(text);
    const current = state.targetList ? state.targetList.split(',').map(s => s.trim()) : [];
    const unique = [...new Set([...current, ...extracted])];
    
    updateState({ targetList: unique.join(', ') });
    return extracted.length;
  };

  const handleDeepMine = async () => {
    if (!state.target) {
      alert("Por favor, defina o Alvo (Órgão ou Doença) primeiro.");
      return;
    }
    
    setIsLoading(true);
    updateState({ miningStatus: 'Conectando ao Gemini AI...' });
    
    let currentTerms = state.targetList 
      ? state.targetList.split(',').map(s => s.trim()).filter(x => x) 
      : [];

    if (currentTerms.length === 0 && state.miningStrategy === 'conservative') {
      alert("O modo 'Faxina' serve para limpar uma lista existente. Por favor, adicione alvos manualmente ou use 'Blue Ocean' para descobrir novos.");
      setIsLoading(false);
      updateState({ miningStatus: '' });
      return;
    }

    if (state.useAI) {
      try {
        updateState({ miningStatus: `Minerando (${state.miningStrategy})...` });
        
        const resultTerms = await gemini.mineNovelTargets(
          state.target,
          state.context || '', 
          currentTerms,
          state.miningStrategy
        );
        
        updateState({ miningStatus: 'Processando resultados...' });

        if (resultTerms.length === 0) {
          alert("A IA não encontrou novos alvos com os parâmetros atuais. Tente mudar a estratégia ou adicionar Contexto.");
        } else {
          if (state.miningStrategy === 'conservative') {
             currentTerms = resultTerms;
          } else {
             const existingSet = new Set(currentTerms.map(t => t.toUpperCase()));
             const newItems = resultTerms.filter(t => !existingSet.has(t.toUpperCase()));
             currentTerms = [...currentTerms, ...newItems];
          }
        }

      } catch (e) {
        console.warn("Deep Mining failed", e);
        alert("Erro na Mineração IA. Verifique se o ambiente possui a API Key configurada.");
      }
    } else if (state.context) {
      const contextTerms = extractEntitiesRegex(state.context);
      const existingSet = new Set(currentTerms.map(t => t.toUpperCase()));
      const newItems = contextTerms.filter(t => !existingSet.has(t.toUpperCase()));
      currentTerms = [...currentTerms, ...newItems];
    }
    
    updateState({ targetList: currentTerms.join(', '), miningStatus: '' });
    setIsLoading(false);
  };

  const executeAnalysis = async () => {
    if (!state.target) {
      alert("Defina um alvo antes de executar.");
      return;
    }
    if (!state.targetList) {
       alert("Sua lista de alvos está vazia. Use a IA para minerar ou adicione manualmente.");
       return;
    }

    setIsLoading(true);
    setProgress(0);
    updateState({ miningStatus: 'Calculando estatísticas (Lambda Score)...' });
    
    const terms = state.targetList.split(',').map(s => s.trim()).filter(s => s.length > 0);
    const results = await performAnalysis(terms, state.target, state.email, (p) => setProgress(p));
    
    updateState({ results, page: 'results', miningStatus: '' });
    setIsLoading(false);
  };

  const fetchDetails = async (molecule: string) => {
    const articles = await fetchArticles(molecule, state.target);
    updateState({ selectedArticleDetails: articles });
  };

  const analyzeArticle = async (title: string, abstract: string) => {
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
