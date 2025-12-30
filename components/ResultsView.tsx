import React from 'react';
import { AppState } from '../types';
import { 
  ArrowLeft, Search, ExternalLink, Bot, Download, Zap 
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

  // Maps the internal status ID to the translated text with Emoji
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
    if (!state.apiKey) {
      alert("Please add your API Key in the settings first.");
      return;
    }
    const confirm = window.confirm("This action uses 1 AI Credit (API Call). Continue?");
    if (confirm) {
      const res = await analyzeArticle(art.title, art.abstract);
      alert(`Lemos Lambda Insight:\n\n${res}`);
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
           <BarChart data={state.results.slice(0, 25)}>
             <XAxis dataKey="molecule" stroke="#666" fontSize={12} tick={{fill: '#999'}} interval={0} angle={-45} textAnchor="end" height={60} />
             <YAxis stroke="#666" fontSize={12} tick={{fill: '#999'}} />
             <Tooltip 
               contentStyle={{ backgroundColor: '#0E1117', border: '1px solid #333' }}
               itemStyle={{ color: '#fff' }}
             />
             <Bar dataKey="ratio">
               {state.results.slice(0, 25).map((entry, index) => (
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
            <thead className="bg-black/20 text-gray-400 text-sm uppercase">
              <tr>
                <th className="p-4">{t.col_mol}</th>
                <th className="p-4">{t.col_status}</th>
                <th className="p-4">{t.col_ratio}</th>
                <th className="p-4">{t.col_pvalue}</th>
                <th className="p-4">{t.col_art_alvo}</th>
                <th className="p-4">{t.col_global}</th>
                <th className="p-4">Action</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-gray-800 text-sm">
              {state.results.map((row, idx) => (
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
          <div className="bg-lemos-card w-full max-w-3xl max-h-[80vh] overflow-y-auto rounded-xl border border-gray-700 p-6">
             <div className="flex justify-between items-center mb-6">
               <h2 className="text-xl font-bold text-white">{t.header_leitura}</h2>
               <button onClick={() => updateState({ selectedArticleDetails: null })} className="text-gray-400 hover:text-white">Close</button>
             </div>
             
             <div className="space-y-4">
               {state.selectedArticleDetails.map((art) => (
                 <div key={art.id} className="bg-black/20 p-4 rounded-lg border border-gray-800">
                    <h4 className="text-lg font-semibold text-blue-300 mb-1">{art.title}</h4>
                    <div className="text-xs text-gray-500 mb-2">{art.journal} • {art.year}</div>
                    <p className="text-sm text-gray-300 mb-3 leading-relaxed">{art.abstract}</p>
                    
                    <div className="flex gap-3 pt-2">
                       <a href={art.link} className="flex items-center gap-1 text-xs text-white bg-gray-700 px-3 py-2 rounded hover:bg-gray-600 transition-colors">
                         <ExternalLink size={12} /> {t.btn_pubmed}
                       </a>
                       {state.apiKey && (
                         <button 
                           onClick={() => handleAIAnalysis(art)}
                           className="flex items-center gap-1 text-xs text-white bg-purple-900/50 border border-purple-500/50 px-3 py-2 rounded hover:bg-purple-900 transition-colors group"
                           title="Consumes AI Token"
                         >
                           <Bot size={12} className="group-hover:text-purple-300" /> 
                           <span className="group-hover:text-purple-300">{t.btn_investigar}</span>
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