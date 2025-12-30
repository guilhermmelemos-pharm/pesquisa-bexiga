import React from 'react';
import { 
  FlaskConical, Trash2, Settings, Key, Play, FileText, AlertTriangle, 
  Sparkles, ShieldCheck, Pill, Dna, Upload, Copy, Check, Loader2
} from 'lucide-react';
import { TEXTS, PRESETS_FRONTEIRA } from '../constants';
import Radar from './Radar';
import { AppState, Language, MiningStrategy } from '../types';

interface HomeViewProps {
  state: AppState;
  t: any;
  updateState: (updates: Partial<AppState>) => void;
  isLoading: boolean;
  progress: number;
  handleDeepMine: () => void;
  executeAnalysis: () => void;
  handlePresetAdd: (name: string) => void;
  handleAddAllPresets: () => void;
  handleFileUpload: (text: string) => number;
  setLang: (l: Language) => void;
  lang: Language;
}

const HomeView: React.FC<HomeViewProps> = ({
  state, t, updateState, isLoading, progress, 
  handleDeepMine, executeAnalysis, handlePresetAdd, 
  handleAddAllPresets, handleFileUpload, setLang, lang
}) => {

  const [copied, setCopied] = React.useState(false);

  const onFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (file) {
      const reader = new FileReader();
      reader.onload = (ev) => {
        const count = handleFileUpload(ev.target?.result as string);
        alert(`Success! ${count} items imported.`);
      };
      reader.readAsText(file);
    }
  };

  const copyPix = () => {
    navigator.clipboard.writeText(t.pix_key);
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  };

  const StrategyButton = ({ 
    id, icon: Icon, label, colorClass, desc 
  }: { 
    id: MiningStrategy, icon: any, label: string, colorClass: string, desc: string 
  }) => (
    <button 
      onClick={() => updateState({ miningStrategy: id })}
      title={desc}
      className={`relative p-2 rounded-lg border text-xs font-bold flex flex-col items-center justify-center gap-1 transition-all h-20 group
        ${state.miningStrategy === id 
          ? `bg-opacity-10 border-opacity-100 ${colorClass} bg-white ring-1 ring-offset-0 ring-white/20` 
          : 'bg-lemos-dark border-gray-700 text-gray-500 hover:border-gray-500 hover:bg-gray-800'
        }`}
    >
      <Icon size={22} className={`mb-1 transition-transform group-hover:scale-110 ${state.miningStrategy === id ? '' : 'opacity-60'}`} />
      <span className={`text-[10px] uppercase tracking-wide ${state.miningStrategy === id ? 'text-white' : ''}`}>{label}</span>
      
      {state.miningStrategy === id && (
        <span className="absolute -top-1 -right-1 flex h-3 w-3">
          <span className={`animate-ping absolute inline-flex h-full w-full rounded-full opacity-75 ${colorClass.split(' ')[0].replace('border', 'bg')}`}></span>
          <span className={`relative inline-flex rounded-full h-3 w-3 ${colorClass.split(' ')[0].replace('border', 'bg')}`}></span>
        </span>
      )}
    </button>
  );

  return (
    <div className="min-h-screen p-4 md:p-8 max-w-7xl mx-auto flex flex-col">
      {/* Header */}
      <header className="mb-8 flex justify-between items-center">
        <div>
          <div className="flex items-center gap-3">
            <span className="text-3xl text-lemos-red">λ</span>
            <h1 className="text-3xl md:text-4xl font-bold tracking-tight">{t.titulo_desk}</h1>
          </div>
          <p className="text-lemos-sub mt-1 text-lg">{t.subtitulo}</p>
        </div>
        <div className="flex gap-2">
          <button onClick={() => setLang(Language.PT)} className={`px-3 py-1 rounded ${lang === Language.PT ? 'bg-lemos-card border border-white/20' : 'opacity-50'}`}>🇧🇷</button>
          <button onClick={() => setLang(Language.EN)} className={`px-3 py-1 rounded ${lang === Language.EN ? 'bg-lemos-card border border-white/20' : 'opacity-50'}`}>🇺🇸</button>
        </div>
      </header>

      <Radar />

      <div className="grid grid-cols-1 lg:grid-cols-3 gap-8 flex-1">
        {/* Main Form */}
        <div className="lg:col-span-2 space-y-6">
          <div className="bg-lemos-card p-6 rounded-xl border border-gray-800">
            <h2 className="text-xl font-semibold mb-4 flex items-center gap-2">
              <FlaskConical size={20} className="text-lemos-red" />
              {t.header_scope}
            </h2>
            
            <div className="space-y-4">
              <div>
                <label className="block text-sm text-gray-400 mb-1">{t.label_email}</label>
                <input 
                  type="email" 
                  className="w-full bg-lemos-dark border border-gray-700 rounded-lg p-3 text-white focus:border-lemos-red focus:outline-none"
                  value={state.email}
                  onChange={(e) => updateState({ email: e.target.value })}
                  placeholder={t.holder_email}
                />
              </div>
              
              <div className="flex gap-2 items-end">
                <div className="flex-1">
                  <label className="block text-sm text-gray-400 mb-1">{t.label_alvo}</label>
                  <input 
                    type="text" 
                    className="w-full bg-lemos-dark border border-gray-700 rounded-lg p-3 text-white focus:border-lemos-red focus:outline-none"
                    value={state.target}
                    onChange={(e) => updateState({ target: e.target.value })}
                    placeholder={t.holder_alvo}
                  />
                  <p className="text-[10px] text-gray-500 mt-1">{t.helper_alvo}</p>
                </div>
                <button 
                  onClick={() => updateState({ target: '' })}
                  className="p-3 bg-lemos-dark border border-gray-700 rounded-lg hover:text-lemos-red h-[50px] mb-[22px]"
                >
                  <Trash2 size={20} />
                </button>
              </div>

              {/* Strategy Selector */}
              {state.useAI && state.apiKey && (
                <div className="space-y-2 pt-2">
                   <div className="flex justify-between items-end mb-1">
                     <label className="text-xs text-gray-400 font-bold uppercase tracking-wider">{t.lbl_ai_strategy}</label>
                     <span className="text-[10px] text-gray-500 italic">{t.helper_ai_strategy}</span>
                   </div>
                   <div className="grid grid-cols-4 gap-3">
                     <StrategyButton 
                       id="conservative" 
                       icon={ShieldCheck} 
                       label={t.strat_clean} 
                       desc={t.desc_clean} 
                       colorClass="border-green-500 text-green-400" 
                     />
                     <StrategyButton 
                       id="repurposing" 
                       icon={Pill} 
                       label={t.strat_repurpose} 
                       desc={t.desc_repurpose} 
                       colorClass="border-orange-500 text-orange-400" 
                     />
                     <StrategyButton 
                       id="mechanism" 
                       icon={Dna} 
                       label={t.strat_mechanism} 
                       desc={t.desc_mechanism} 
                       colorClass="border-purple-500 text-purple-400" 
                     />
                     <StrategyButton 
                       id="blue_ocean" 
                       icon={Sparkles} 
                       label={t.strat_blue_ocean} 
                       desc={t.desc_blue_ocean} 
                       colorClass="border-blue-500 text-blue-400" 
                     />
                   </div>
                </div>
              )}

              <div className="bg-lemos-dark/50 p-3 rounded-lg border border-gray-800">
                 <label className="block text-xs text-lemos-red font-bold mb-1 uppercase tracking-wide">
                   {t.lbl_contexto_tit}
                 </label>
                 <textarea 
                   value={state.context}
                   onChange={(e) => updateState({ context: e.target.value })}
                   className="w-full bg-lemos-dark border border-gray-700 rounded-lg p-2 text-sm h-16 focus:border-lemos-red focus:outline-none"
                   placeholder="Ex: 'Pulmão' / 'Hypoxia pathway'..."
                 />
                 <p className="text-[10px] text-gray-500 mt-1">
                   {t.lbl_contexto_help}
                 </p>
              </div>

              {/* Warning about Missing API Key */}
              {!state.apiKey && state.useAI && (
                <div className="bg-yellow-900/20 border border-yellow-600/50 p-3 rounded-lg flex gap-3 items-start animate-pulse">
                  <AlertTriangle className="text-yellow-500 shrink-0" size={18} />
                  <div>
                    <h4 className="text-yellow-400 text-xs font-bold uppercase">{t.warn_no_key_tit}</h4>
                    <p className="text-yellow-200/80 text-xs mt-1">
                      {t.warn_no_key_desc}
                    </p>
                  </div>
                </div>
              )}

              {/* MINING BUTTON / PROGRESS BAR AREA */}
              <div className="relative">
                {isLoading ? (
                  /* ESTADO DE CARREGAMENTO (Barra de Progresso) */
                  <div className="w-full bg-lemos-dark border border-blue-500/30 rounded-xl p-4 shadow-lg shadow-blue-900/10">
                    <div className="flex items-center justify-between mb-2">
                      <span className="text-xs font-bold text-blue-400 uppercase tracking-widest flex items-center gap-2">
                        <Loader2 className="animate-spin" size={14} />
                        Neuro-Mining
                      </span>
                      <span className="text-xs text-blue-300 font-mono animate-pulse">
                        Gemini 2.0 Thinking...
                      </span>
                    </div>
                    
                    {/* A Barra Animada */}
                    <div className="w-full h-2 bg-gray-800 rounded-full overflow-hidden relative">
                      <div className="absolute top-0 left-0 h-full bg-blue-500 w-1/3 animate-[slide_1.5s_ease-in-out_infinite]"></div>
                    </div>
                    
                    <p className="text-[10px] text-gray-500 mt-2 font-mono text-center">
                      Cruzando dados de "{state.target}" com literatura recente...
                    </p>
                  </div>
                ) : (
                  /* ESTADO NORMAL (Botão) */
                  <button 
                    onClick={handleDeepMine}
                    disabled={isLoading}
                    className={`w-full font-bold py-4 rounded-xl text-lg shadow-lg transition-all flex items-center justify-center gap-2 relative overflow-hidden group
                      ${!state.apiKey 
                        ? 'bg-gray-700 text-gray-300 hover:bg-gray-600' 
                        : 'bg-lemos-red hover:bg-red-600 text-white shadow-red-900/20'
                      }`}
                  >
                    <div className="absolute inset-0 bg-white/10 translate-y-full group-hover:translate-y-0 transition-transform duration-300"></div>
                    
                    {state.apiKey 
                      ? `${t.btn_run_prefix} ${t[`strat_${state.miningStrategy}`] || state.miningStrategy.toUpperCase()}`
                      : "🔍 Minerar (Modo Básico)"
                    }
                  </button>
                )}
              </div>
            </div>
          </div>

          {/* Presets */}
          <div className="bg-lemos-card p-6 rounded-xl border border-gray-800">
             <h3 className="text-lg font-semibold mb-4 text-gray-300">{t.expander_presets}</h3>
             <button 
              onClick={handleAddAllPresets}
              className="w-full mb-4 px-4 py-3 bg-blue-900/30 border border-blue-500/50 rounded-lg hover:bg-blue-900/50 text-blue-200 text-sm font-bold transition-colors"
             >
               {t.btn_add_all}
             </button>
             <div className="flex flex-wrap gap-2">
               {Object.keys(PRESETS_FRONTEIRA).map(key => (
                 <button 
                  key={key}
                  onClick={() => handlePresetAdd(key)}
                  className="px-4 py-2 bg-lemos-dark border border-gray-700 rounded-full hover:border-lemos-red hover:text-lemos-red transition-colors text-xs md:text-sm"
                 >
                   + {key}
                 </button>
               ))}
             </div>
          </div>
        </div>

        {/* Sidebar Config */}
        <div className="space-y-6">
          <div className="bg-lemos-card p-6 rounded-xl border border-gray-800">
            <h2 className="text-xl font-semibold mb-4 flex items-center gap-2">
              <Settings size={20} className="text-lemos-sub" />
              {t.header_config}
            </h2>

            <div className="space-y-4">
              <div>
                 <label className="block text-sm text-gray-400 mb-1 flex items-center gap-1">
                   <Key size={14} /> Google API Key
                 </label>
                 <input 
                  type="password" 
                  className={`w-full bg-lemos-dark border rounded-lg p-2 text-sm focus:outline-none ${!state.apiKey ? 'border-yellow-500/50' : 'border-gray-700'}`}
                  value={state.apiKey}
                  onChange={(e) => updateState({ apiKey: e.target.value })}
                  placeholder={t.placeholder_key}
                />
                <div className="mt-1 text-xs text-gray-500">
                  <a href="https://aistudio.google.com/app/apikey" target="_blank" rel="noreferrer" className="underline hover:text-lemos-red">{t.link_key}</a>
                </div>
                <div className="mt-3 flex items-center gap-2">
                  <input 
                    type="checkbox" 
                    id="aiToggle"
                    checked={state.useAI}
                    onChange={(e) => updateState({ useAI: e.target.checked })}
                    className="accent-lemos-red"
                  />
                  <label htmlFor="aiToggle" className="text-sm text-gray-300">{t.expander_ia}</label>
                </div>
              </div>
              
              <hr className="border-gray-700" />
              
              <div>
                 <label className="block text-sm text-gray-400 mb-1">{t.slider_tempo}</label>
                 <div className="flex items-center gap-2 text-sm">
                    <input 
                      type="number" 
                      value={state.yearStart} 
                      onChange={(e) => updateState({ yearStart: parseInt(e.target.value) })}
                      className="w-20 bg-lemos-dark border border-gray-700 rounded p-1"
                    />
                    <span>-</span>
                    <input 
                      type="number" 
                      value={state.yearEnd} 
                      onChange={(e) => updateState({ yearEnd: parseInt(e.target.value) })}
                      className="w-20 bg-lemos-dark border border-gray-700 rounded p-1"
                    />
                 </div>
              </div>

              <div>
                 <label className="block text-sm text-gray-400 mb-1 flex items-center gap-1">
                   <Upload size={14} /> {t.upload_label}
                 </label>
                 <input 
                  type="file"
                  accept=".txt,.csv" 
                  onChange={onFileChange}
                  className="w-full text-xs text-gray-400 file:mr-2 file:py-2 file:px-4 file:rounded-full file:border-0 file:text-xs file:font-semibold file:bg-lemos-dark file:text-white hover:file:bg-gray-700"
                 />
                 <p className="text-[10px] text-gray-500 mt-1">
                    {t.upload_helper}
                 </p>
              </div>
            </div>
          </div>

          {/* List Preview */}
          {state.targetList && (
            <div className="bg-lemos-dark p-4 rounded-xl border border-gray-800">
               <div className="flex justify-between items-center mb-2">
                 <span className="text-xs font-mono text-gray-400">{t.lbl_targets_preview}</span>
                 <span className="text-xs font-bold text-lemos-red">{state.targetList.split(',').length} {t.lbl_targets_count}</span>
               </div>
               <textarea 
                 className="w-full bg-transparent text-xs font-mono text-gray-300 h-24 focus:outline-none resize-none border border-transparent focus:border-gray-700 rounded p-1"
                 value={state.targetList}
                 onChange={(e) => updateState({ targetList: e.target.value })}
                 placeholder="Digite ou edite seus alvos aqui..."
               />
               <button 
                 onClick={executeAnalysis}
                 disabled={isLoading}
                 className="mt-2 w-full bg-blue-600 hover:bg-blue-500 text-white font-bold py-3 rounded-lg flex items-center justify-center gap-2"
               >
                 {isLoading ? (
                   <span>{t.titulo_processando} {Math.round(progress)}%</span>
                 ) : (
                   <>
                    <Play size={16} fill="white" />
                    {t.btn_executar}
                   </>
                 )}
               </button>
               <div className="mt-2 text-center">
                  <button onClick={() => updateState({ targetList: '' })} className="text-xs text-red-400 hover:text-red-300">
                     {t.btn_limpar}
                  </button>
               </div>
            </div>
          )}
        </div>
      </div>

      {/* FOOTER */}
      <footer className="mt-12 border-t border-gray-800 pt-8 pb-4">
        <div className="grid grid-cols-1 md:grid-cols-3 gap-8">
          <div className="md:col-span-2">
             <h4 className="text-sm font-bold text-gray-400 mb-2">{t.footer_citar}</h4>
             <details className="bg-lemos-dark border border-gray-700 rounded-lg p-3 group">
                <summary className="text-sm cursor-pointer font-medium hover:text-lemos-red flex items-center gap-2 text-gray-300">
                   <FileText size={16} />
                   {t.citar_titulo}
                </summary>
                <div className="mt-3 bg-black/30 p-2 rounded text-xs font-mono text-gray-300 break-all border border-gray-800 select-all">
                   {t.citar_texto}
                </div>
             </details>
          </div>
          <div>
             <h4 className="text-sm font-bold text-gray-400 mb-2">{t.apoio_titulo}</h4>
             <div className="flex gap-2">
                <input 
                  readOnly 
                  value={t.pix_key} 
                  className="bg-black/30 border border-gray-700 text-xs text-gray-300 p-2 rounded flex-1 focus:outline-none font-mono"
                />
                <button 
                  onClick={copyPix}
                  className="bg-lemos-card hover:bg-gray-700 border border-gray-700 text-white p-2 rounded transition-colors"
                  title={t.btn_copy}
                >
                  {copied ? <Check size={16} className="text-green-500" /> : <Copy size={16} />}
                </button>
             </div>
             <p className="text-xs text-gray-500 mt-2">{t.apoio_desc}</p>
          </div>
        </div>
      </footer>
    </div>
  );
};

export default HomeView;
