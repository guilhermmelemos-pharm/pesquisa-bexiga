import { useState, useEffect } from 'react';
import { AppState, AnalysisResult, Article } from '../types';
import { GeminiService } from '../services/geminiService';
import { performAnalysis, fetchArticles, parseSimpleUserList } from '../services/scienceService';
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
    miningStrategy: 'blue_ocean',
    selectedArticleDetails: null
  });

  // SEPARAÇÃO DOS ESTADOS DE CARREGAMENTO
  const [isMining, setIsMining] = useState(false);     // Só para IA
  const [isAnalyzing, setIsAnalyzing] = useState(false); // Só para Estatística
  
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
    const unique = [...new Set([...current, ...allItems])];
    updateState({ targetList: unique.join(', ') });
  };

  const handleFileUpload = (text: string) => {
    const extracted = parseSimpleUserList(text);
    const current = state.targetList ? state.targetList.split(',').map(s => s.trim()) : [];
    const unique = [...new Set([...current, ...extracted])];
    updateState({ targetList: unique.join(', ') });
    return extracted.length;
  };

  // --- MINERAÇÃO (IA) ---
  const handleDeepMine = async () => {
    if (!state.target) {
      alert("⚠️ Defina o Alvo Principal (Doença/Órgão) antes de minerar.");
      return;
    }
    
    setIsMining(true); // Ativa apenas o loading da IA
    let currentTerms = state.targetList 
      ? state.targetList.split(',').map(s => s.trim()).filter(x => x) 
      : [];

    // Se estiver vazio e sem IA, aplica um preset básico
    if (currentTerms.length === 0 && (!state.useAI || !state.apiKey)) {
       const randomPreset = Object.values(PRESETS_FRONTEIRA)[0];
       currentTerms = [...randomPreset];
    }

    // Se tiver IA, chama o Gemini
    if (state.useAI && state.apiKey) {
      try {
        const resultTerms = await gemini.mineNovelTargets(
          state.target,
          state.context || '', 
          currentTerms,
          state.miningStrategy
        );
        
        if (state.miningStrategy === 'conservative') {
           currentTerms = resultTerms;
        } else {
           const existingSet = new Set(currentTerms.map(t => t.toUpperCase()));
           const newItems = resultTerms.filter(t => !existingSet.has(t.toUpperCase()));
           currentTerms = [...currentTerms, ...newItems];
        }

      } catch (e) {
        console.warn("Deep Mining falhou ou foi cancelado", e);
        // Não alertamos aqui pois o geminiService já alerta
      }
    }
    
    updateState({ targetList: currentTerms.join(', ') });
    setIsMining(false); // Destrava a interface
  };

  // --- ANÁLISE (ESTATÍSTICA LOCAL) ---
  const executeAnalysis = async () => {
    if (!state.target) {
        alert("Defina um alvo primeiro.");
        return;
    }
    if (!state.targetList) {
        alert("A lista de alvos está vazia. Adicione presets ou use a IA.");
        return;
    }

    setIsAnalyzing(true); // Ativa apenas o loading da Análise
    setProgress(0);
    
    try {
        const terms = state.targetList.split(',').map(s => s.trim()).filter(s => s.length > 0);
        // Adicionei um pequeno delay inicial para a UI atualizar
        await new Promise(r => setTimeout(r, 100));
        
        const results = await performAnalysis(terms, state.target, state.email, (p) => setProgress(p));
        
        updateState({ results, page: 'results' });
    } catch (error) {
        alert("Erro na análise estatística. Tente limpar a lista e recomeçar.");
        console.error(error);
    } finally {
        setIsAnalyzing(false); // Destrava a interface sempre
    }
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
    isMining,      // Exportando separado
    isAnalyzing,   // Exportando separado
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
