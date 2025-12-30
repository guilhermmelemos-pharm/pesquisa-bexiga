
import React, { useState, useMemo } from 'react';
import { AppState, AnalysisResult } from '../types';
import { 
  ArrowLeft, Search, ExternalLink, Bot, Download, Zap, Loader2, ArrowUpDown, ArrowUp, ArrowDown
} from 'lucide-react';
import { 
  BarChart, Bar, XAxis, YAxis, Tooltip, ResponsiveContainer, Cell 
} from 'recharts';

interface ResultsViewProps {
  state: AppState;
  t: any;
  updateState: (updates: Partial<AppState>) => void;
  fetchDetails: (mol: string) => void;
  analyzeArticle: (title: string, abs: string) => Promise<string>;
}

const ResultsView: React.FC<ResultsViewProps> = ({ 
  state, t, updateState, fetchDetails, analyzeArticle 
}) => {
  const topResult = state.results[0];
  
  // Local state to store AI analysis per article ID
  const [aiInsights, setAiInsights] = useState<Record<string, string>>({});
  const [analyzingIds, setAnalyzingIds] = useState<Record<string, boolean>>({});

  // Sorting State - Default by sortScore (Implied importance)
  const [sortConfig, setSortConfig] = useState<{ key: keyof AnalysisResult, direction: 'asc' | 'desc' }>({
    key: 'sortScore',
    direction: 'desc'
  });

  const sortedResults = useMemo(() => {
    let sortableItems = [...state.results];
    if (sortConfig !== null) {
      sortableItems.sort((a, b) => {
        let valA: any = a[sortConfig.key];
        let valB: any = b[sortConfig.key];

        // Ensure numbers are treated as numbers
        if (['ratio', 'targetArticles', 'globalArticles', 'sortScore'].includes(sortConfig.key)) {
            valA = Number(valA);
            valB = Number(valB);
        }
        // pValue is string but represents a number, better to sort numerically
        if (sortConfig.key === 'pValue') {
            valA = parseFloat(valA);
            valB = parseFloat(valB);
        }
        // Status string sort
        if (sortConfig.key === 'status') {
             valA = valA.toString();
             valB = valB.toString();
        }

        if (valA < valB) {
          return sortConfig.direction === 'asc' ? -1 : 1;
        }
        if (valA > valB) {
          return sortConfig.direction === 'asc' ? 1 : -1;
        }
        return 0;
      });
    }
    return sortableItems;
  }, [state.results, sortConfig]);

  const requestSort = (key: keyof AnalysisResult) => {
    let direction: 'asc' | 'desc' = 'desc'; // Default to desc for high numbers
    if (sortConfig && sortConfig.key === key && sortConfig.direction === 'desc') {
      direction = 'asc';
    }
    setSortConfig({ key, direction });
  };

  const getSortIcon = (key: keyof AnalysisResult) => {
    if (sortConfig?.key === key) {
        return sortConfig.direction === 'asc' 
          ? <ArrowUp size={14} className="inline ml-1 text-lemos-red" />
          : <ArrowDown size={14} className="inline ml-1 text-lemos-red" />;
    }
    return <ArrowUpDown size={14} className="inline ml-1 opacity-20 group-hover:opacity-100 transition-opacity" />;
  };

  const getStatusColor = (status: string) => {
    switch(status) {
      case 'Gold': return '#FFD700';
      case 'Trending': return '#FF4B4B';
      case 'Blue Ocean': return '#00BFFF';
      case 'Embryonic': return '#32CD32';
      case 'Saturated': return '#808080';
      default: return '#555';
    }
  };

  const getStatusLabel = (status: string) => {
    switch(status) {
      case 'Gold': return t.tag_gold;
      case 'Trending': return t.tag_trending;
      case 'Blue Ocean': return t.tag_blue_ocean;
      case 'Embryonic': return t.tag_embryonic;
      case 'Saturated': return t.tag_saturated;
      case 'Ghost': return t.tag_ghost;
      default: return t.tag_neutral;
    }
  };

  const handleAIAnalysis = async (art: any) => {
    // Set analyzing state
    setAnalyzingIds(prev => ({ ...prev, [art.id]: true }));
    
    try {
      const res = await analyzeArticle(art.title, art.abstract);
      setAiInsights(prev => ({ ...prev, [art.id]: res }));
    } catch (e) {
      console.error(e);
    } finally {
      setAnalyzingIds(prev => ({ ...prev, [art.id]: false }));
    }
  };

  const downloadCSV = () => {
    if (!state.results.length) return;
    
    const headers = [t.col_mol, "Status", "Ratio", "P-Value", "Hits Target", "Hits Global"];
    const rows = state.results.map(r => [
      r.molecule,
      r.status,
      r.ratio,
      r.pValue,
      r.targetArticles,
      r.globalArticles
    ]);

    const csvContent = [
      headers.join(","),
      ...rows.map(r => r.join(","))
    ].join("\n");

    const blob = new Blob([csvContent], { type: "text/csv;charset=utf-8;" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.setAttribute("href", url);
    link.setAttribute("download", `lemos_lambda_${state.target}_${new Date().toISOString().slice(0,10)}.csv`);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
  };

  return (
    <div className="min-h-screen p-4 md:p-8 max-w-7xl mx-auto">
      {/* Nav */}
      <div className="flex items-center justify-between mb-8">
        <div className="flex items-center gap-4">
          <button 
            onClick={() => updateState({ page: 'home' })}
            className="p-2 bg-lemos-card hover:bg-gray-700 rounded-lg transition-colors"
          >
            <ArrowLeft size={24} />
          </button>
          <h1 className="text-2xl font-bold">{t.resultados}: <span className="text-lemos-red">{state.target}</span></h1>
        </div>
        
        <button 
          onClick={downloadCSV}
          className="flex items-center gap-2 bg-green-700 hover:bg-green-600 text-white px-4 py-2 rounded-lg font-bold transition-colors"
        >
          <Download size={18} /> {t.btn_baixar}
        </button>
      </div>

      {/* Top Metric Cards */}
      {topResult && (
        <div className="grid grid-cols-1 md:grid-cols-3 gap-6 mb-8">
          <div className="bg-lemos-card p-6 rounded-xl border border-gray-700">
            <span className="text-gray-400 text-sm">{t.metric_top}</span>
            <div className="flex items-end justify-between mt-2">
               <span className="text-3xl font-bold truncate pr-2">{topResult.molecule}</span>
               <span className={`px-2 py-1 rounded text-xs font-bold bg-black/30 shrink-0`} style={{color: getStatusColor(topResult.status)}}>
                  {getStatusLabel(topResult.status)}
               </span>
            </div>
          </div>
          <div className="bg-lemos-card p-6 rounded-xl border border-gray-700">
            <span className="text-gray-400 text-sm">{t.metric_score}</span>
            <div className="mt-2 text-3xl font-bold text-blue-400">{topResult.ratio}x</div>
          </div>
          <div className="bg-lemos-card p-6 rounded-xl border border-gray-700">
            <span className="text-gray-400 text-sm">P-Value</span>
            <div className="mt-2 text-3xl font-bold text-emerald-400">{topResult.pValue}</div>
          </div>
        </div>
      )}

      {/* Charts */}
      <div className="bg-lemos-card p-6 rounded-xl border border-gray-800 mb-8 h-96">
         <h3 className="text-lg font-semibold mb-4 text-gray-300">{t.header_heatmap}</h3>
         <ResponsiveContainer width="100%" height="100%">
           <BarChart data={sortedResults.slice(0, 25)}>
             <XAxis dataKey="molecule" stroke="#666" fontSize={12} tick={{fill: '#999'}} interval={0} angle={-45} textAnchor="end" height={60} />
             <YAxis stroke="#666" fontSize={12} tick={{fill: '#999'}} />
             <Tooltip 
               contentStyle={{ backgroundColor: '#0E1117', border: '1px solid #333' }}
               itemStyle={{ color: '#fff' }}
             />
             <Bar dataKey="ratio">
               {sortedResults.slice(0, 25).map((entry, index) => (
                 <Cell key={`cell-${index}`} fill={getStatusColor(entry.status)} />
               ))}
             </Bar>
           </BarChart>
         </ResponsiveContainer>
      </div>

      {/* Data Table */}
      <div className="bg-lemos-card rounded-xl border border-gray-800 overflow-hidden mb-8">
        <div className="overflow-x-auto">
          <table className="w-full text-left">
            <thead className="bg-black/20 text-gray-400 text-sm uppercase cursor-pointer select-none">
              <tr>
                <th className="p-4 group hover:bg-white/5 whitespace-nowrap" onClick={() => requestSort('molecule')}>
                    {t.col_mol} {getSortIcon('molecule')}
                </th>
                <th className="p-4 group hover:bg-white/5 whitespace-nowrap" onClick={() => requestSort('status')}>
                    {t.col_status} {getSortIcon('status')}
                </th>
                <th className="p-4 group hover:bg-white/5 whitespace-nowrap" onClick={() => requestSort('ratio')}>
                    {t.col_ratio} {getSortIcon('ratio')}
                </th>
                <th className="p-4 group hover:bg-white/5 whitespace-nowrap" onClick={() => requestSort('pValue')}>
                    {t.col_pvalue} {getSortIcon('pValue')}
                </th>
                <th className="p-4 group hover:bg-white/5 whitespace-nowrap" onClick={() => requestSort('targetArticles')}>
                    {t.col_art_alvo} {getSortIcon('targetArticles')}
                </th>
                <th className="p-4 group hover:bg-white/5 whitespace-nowrap" onClick={() => requestSort('globalArticles')}>
                    {t.col_global} {getSortIcon('globalArticles')}
                </th>
                <th className="p-4">Action</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-gray-800 text-sm">
              {sortedResults.map((row, idx) => (
                <tr key={idx} className="hover:bg-white/5 transition-colors">
                  <td className="p-4 font-medium text-white">{row.molecule}</td>
                  <td className="p-4">
                    <span className="px-2 py-1 rounded text-xs font-bold bg-black/30 whitespace-nowrap" style={{color: getStatusColor(row.status)}}>
                      {getStatusLabel(row.status)}
                    </span>
                  </td>
                  <td className="p-4 text-white">{row.ratio}</td>
                  <td className="p-4 font-mono text-emerald-400">{row.pValue}</td>
                  <td className="p-4 text-gray-400">{row.targetArticles}</td>
                  <td className="p-4 text-gray-400">{row.globalArticles}</td>
                  <td className="p-4">
                     <button 
                      onClick={() => fetchDetails(row.molecule)}
                      className="text-blue-400 hover:text-blue-300 flex items-center gap-1 whitespace-nowrap"
                     >
                       <Search size={14} /> {t.btn_ver_papers}
                     </button>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>

      {/* Article Detail View (Drill Down) */}
      {state.selectedArticleDetails && (
        <div className="fixed inset-0 bg-black/80 backdrop-blur-sm flex items-center justify-center p-4 z-50">
          <div className="bg-lemos-card w-full max-w-3xl max-h-[80vh] overflow-y-auto rounded-xl border border-gray-700 p-6 shadow-2xl">
             <div className="flex justify-between items-center mb-6">
               <h2 className="text-xl font-bold text-white">{t.header_leitura}</h2>
               <button onClick={() => updateState({ selectedArticleDetails: null })} className="text-gray-400 hover:text-white bg-gray-800 p-2 rounded-full">✕</button>
             </div>
             
             <div className="space-y-4">
               {state.selectedArticleDetails.map((art) => (
                 <div key={art.id} className="bg-black/40 p-5 rounded-lg border border-gray-700/50">
                    <div className="flex justify-between items-start gap-4">
                      <h4 className="text-lg font-semibold text-blue-300 mb-1 leading-tight">{art.title}</h4>
                      <span className="text-xs text-gray-500 whitespace-nowrap bg-black/50 px-2 py-1 rounded">{art.year}</span>
                    </div>
                    <div className="text-xs text-gray-400 mb-3 italic">{art.journal}</div>
                    <p className="text-sm text-gray-300 mb-4 leading-relaxed">{art.abstract}</p>
                    
                    {/* Inline AI Insight */}
                    {aiInsights[art.id] && (
                       <div className="mb-4 bg-purple-900/20 border-l-4 border-purple-500 p-3 rounded-r-lg">
                         <div className="flex items-center gap-2 mb-1">
                           <Bot size={14} className="text-purple-400" />
                           <span className="text-xs font-bold text-purple-400 uppercase tracking-wider">AI Insight</span>
                         </div>
                         <p className="text-sm text-gray-200">{aiInsights[art.id]}</p>
                       </div>
                    )}

                    <div className="flex gap-3 pt-2 border-t border-gray-800/50">
                       <a href={art.link} className="flex items-center gap-1 text-xs text-white bg-gray-700 px-3 py-2 rounded hover:bg-gray-600 transition-colors font-medium">
                         <ExternalLink size={12} /> {t.btn_pubmed}
                       </a>
                       {!aiInsights[art.id] && (
                         <button 
                           onClick={() => handleAIAnalysis(art)}
                           disabled={analyzingIds[art.id]}
                           className="flex items-center gap-1 text-xs text-white bg-purple-900/50 border border-purple-500/50 px-3 py-2 rounded hover:bg-purple-900 transition-colors group disabled:opacity-50"
                           title="Consumes AI Token"
                         >
                           {analyzingIds[art.id] ? <Loader2 size={12} className="animate-spin" /> : <Bot size={12} className="group-hover:text-purple-300" />}
                           <span className="group-hover:text-purple-300">{analyzingIds[art.id] ? t.spinner_investigando : t.btn_investigar}</span>
                           <Zap size={10} className="text-yellow-500 ml-1" />
                         </button>
                       )}
                    </div>
                 </div>
               ))}
             </div>
          </div>
        </div>
      )}
    </div>
  );
};

export default ResultsView;
