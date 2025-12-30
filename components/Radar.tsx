import React from 'react';
import { MOCK_NEWS } from '../constants';
import { BookOpen } from 'lucide-react';

const Radar: React.FC = () => {
  return (
    <div className="w-full bg-lemos-card border border-gray-800 rounded-xl p-4 mb-6">
      <div className="flex items-center gap-2 mb-3 text-lemos-sub">
        <span className="text-sm uppercase tracking-widest font-bold">Scientific Radar</span>
      </div>
      <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
        {MOCK_NEWS.map((news) => (
          <div key={news.id} className="group cursor-pointer">
            <div className="h-32 w-full overflow-hidden rounded-lg mb-2 relative">
               <img 
                 src={news.image} 
                 alt={news.title}
                 className="w-full h-full object-cover transition-transform duration-500 group-hover:scale-110 opacity-80 group-hover:opacity-100"
               />
               <div className="absolute bottom-0 left-0 right-0 bg-gradient-to-t from-black/90 to-transparent p-2">
                 <span className="text-xs text-lemos-red font-mono">{news.source}</span>
               </div>
            </div>
            <h4 className="text-sm font-semibold text-gray-200 leading-snug group-hover:text-lemos-red transition-colors">
              {news.title}
            </h4>
            <a href={news.link} className="text-xs text-gray-500 flex items-center gap-1 mt-1 hover:text-white">
               <BookOpen size={12} /> Read Paper
            </a>
          </div>
        ))}
      </div>
    </div>
  );
};

export default Radar;