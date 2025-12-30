import React, { useState } from 'react';
import { Language } from './types';
import { TEXTS } from './constants';
import { useScienceEngine } from './hooks/useScienceEngine';
import HomeView from './components/HomeView';
import ResultsView from './components/ResultsView';

const App: React.FC = () => {
  const [lang, setLang] = useState<Language>(Language.PT);
  const t = TEXTS[lang];

  // Logic is now completely separated in the hook
  const engine = useScienceEngine();

  // --- Render Orchestration ---
  if (engine.state.page === 'home') {
    return (
      <HomeView 
        state={engine.state}
        t={t}
        updateState={engine.updateState}
        isLoading={engine.isLoading}
        progress={engine.progress}
        handleDeepMine={engine.handleDeepMine}
        executeAnalysis={engine.executeAnalysis}
        handlePresetAdd={engine.handlePresetAdd}
        handleAddAllPresets={engine.handleAddAllPresets}
        handleFileUpload={engine.handleFileUpload}
        setLang={setLang}
        lang={lang}
      />
    );
  }

  if (engine.state.page === 'results') {
    return (
      <ResultsView 
        state={engine.state}
        t={t}
        updateState={engine.updateState}
        fetchDetails={engine.fetchDetails}
        analyzeArticle={engine.analyzeArticle}
      />
    );
  }

  return null;
};

export default App;