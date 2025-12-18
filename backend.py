import pandas as pd
from Bio import Entrez
import time
import re
from deep_translator import GoogleTranslator
from datetime import datetime
import feedparser
import random
import nltk
from nltk.corpus import words, stopwords
import constantes as c

# --- CONFIGURAÇÃO NLTK ---
try:
    nltk.data.find('corpora/words')
except LookupError:
    try:
        nltk.download('words', quiet=True)
        nltk.download('stopwords', quiet=True)
    except: pass # Fallback silencioso se não der para baixar

ENGLISH_WORDS = set(word.lower() for word in words.words()) if 'words' in locals() else set()
STOP_WORDS = set(stopwords.words('english')) if 'stopwords' in locals() else set()

# Regex ajustado
REGEX_COMPLEXO = r'\b([A-Z][a-zA-Z0-9\-]*(?:\s[A-Z][a-zA-Z0-9\-]+){0,3})\b'

IMAGENS_SCIENCE = [
    "https://images.unsplash.com/photo-1532094349884-543bc11b234d?auto=format&fit=crop&w=600&q=80", 
    "https://images.unsplash.com/photo-1576086213369-97a306d36557?auto=format&fit=crop&w=600&q=80", 
    "https://images.unsplash.com/photo-1507413245164-6160d8298b31?auto=format&fit=crop&w=600&q=80"
]

# --- LÓGICA CIENTÍFICA CENTRALIZADA (Reprodutibilidade) ---
def classificar_oportunidade(n_alvo, n_fonte):
    """
    Define matematicamente o status do alvo.
    Retorna: (Label_Visual, Score_Ordenacao)
    """
    if n_alvo == 0:
        if n_fonte >= 40: return "💎 Blue Ocean", 1000 # Alta relevância global, zero local
        return "👻 Fantasma", 0 # Irrelevante globalmente
    
    ratio = (n_fonte / n_alvo)
    
    if n_alvo <= 15: return "🌱 Embrionário", 500 # Poucos papers, campo nascendo
    if ratio >= 20: return "🚀 Tendência", 100 # Crescendo muito
    if ratio >= 5: return "🥇 Ouro", 50 # Sólido
    if ratio < 2: return "🔴 Saturado", 10 # Muita gente já estudou
    
    return "⚖️ Neutro", 20

def validar_termo_bioquimico(termo):
    # (Mantido igual à versão 1.6, apenas encapsulado para segurança)
    t_limpo = termo.strip()
    t_lower = t_limpo.lower()
    
    whitelist = ["tmao", "pfas", "nmn", "cbd", "cbg", "lps", "atp", "nad", "dna", "rna", "gpr", "rock", "mtor", "ampk", "yoda1"]
    if any(vip in t_lower for vip in whitelist): return True

    if len(t_limpo) < 3 or t_limpo.isdigit(): return False
    
    if t_lower in STOP_WORDS: return False
    
    # Filtro geográfico hardcoded para segurança extra
    if any(x in t_lower for x in ["university", "hospital", "china", "usa", "department", "fig", "table"]): return False

    return True

def buscar_alvos_emergentes_pubmed(orgao_alvo, email):
    if not orgao_alvo or not email: return []
    Entrez.email = email
    ano_fim = datetime.now().year
    ano_inicio = ano_fim - 3
    
    query = (
        f"({orgao_alvo}[Title/Abstract]) AND ("
        f"\"mechanism\" OR \"receptor\" OR \"pathway\" OR \"target\" OR \"expression\" OR "
        f"\"activation\" OR \"inhibition\" OR \"microplastics\" OR \"TMAO\" OR \"metabolite\")"
        f" AND (\"{ano_inicio}\"[Date - Publication] : \"{ano_fim}\"[Date - Publication])"
        f" NOT (Review[Publication Type])" 
    )
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=80, sort="relevance")
        ids = Entrez.read(handle)["IdList"]
        if not ids: return []
        
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
        texto = handle.read()
        
        candidatos = re.findall(REGEX_COMPLEXO, texto)
        frequencia = {t.strip(" .,();:"): 0 for t in candidatos}
        for t in candidatos: frequencia[t.strip(" .,();:")] += 1

        termos_limpos = []
        blacklist_lower = [x.lower() for x in c.BLACKLIST_GERAL]
        
        for t in frequencia:
            t_lower = t.lower()
            if any(bad in t_lower for bad in blacklist_lower): continue
            if not validar_termo_bioquimico(t): continue
            termos_limpos.append(t)
            
        return sorted(list(set(termos_limpos)), key=lambda x: frequencia[x], reverse=True)[:35]
    except Exception as e:
        print(f"Erro PubMed: {e}")
        return []

def consultar_pubmed_count(termo_farmaco, termo_orgao, email, y_start, y_end):
    if not email: return 0
    Entrez.email = email
    time.sleep(0.1) 
    query = f"(\"{termo_farmaco}\") AND ({termo_orgao}) AND {y_start}:{y_end}[DP]" if termo_orgao else f"(\"{termo_farmaco}\") AND {y_start}:{y_end}[DP]"
    try: return int(Entrez.read(Entrez.esearch(db="pubmed", term=query, retmax=0))["Count"])
    except: return 0

def extrair_conclusao_real(abstract_text, translator):
    if not abstract_text: return "Resumo não disponível."
    match = re.search(r'(Conclusion|Conclusions|In conclusion|Summary|Significance|Results suggest that)(.*)', abstract_text, re.IGNORECASE | re.DOTALL)
    texto_bruto = match.group(2).strip() if match else abstract_text[:500]
    try: return translator.translate(texto_bruto[:600]) + "..."
    except: return texto_bruto[:600] + "..."

def buscar_resumos_detalhados(termo_alvo, orgao, email, y_start, y_end, lang_target, limit=5):
    if not email: return []
    Entrez.email = email
    query = f"(\"{termo_alvo}\") AND ({orgao}) AND {y_start}:{y_end}[DP]"
    translator = GoogleTranslator(source='auto', target=lang_target)
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=limit, sort="relevance")
        ids = Entrez.read(handle)["IdList"]
        if not ids: return []
        records = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text").read().split("\n\n")
        artigos = []
        for art_text in records:
            art_data = {"PMID": "N/A", "Title": "S/T", "Source": "N/A", "Abstract": ""}
            for line in art_text.split("\n"):
                if len(line)<4: continue
                tag, content = line[:4].strip(), line[6:]
                if tag=="PMID": art_data["PMID"]=content
                elif tag=="TI": art_data["Title"]=content
                elif tag=="TA": art_data["Source"]=content
                elif tag=="AB": art_data["Abstract"]+=content + " "
            if art_data["PMID"]!="N/A":
                art_data["Resumo_IA"] = extrair_conclusao_real(art_data["Abstract"], translator)
                art_data["Link"] = f"https://pubmed.ncbi.nlm.nih.gov/{art_data['PMID']}/"
                artigos.append(art_data)
        return artigos
    except: return []

def buscar_todas_noticias(lang_code):
    # FALLBACK SILENCIOSO: Se der erro, retorna lista vazia e o app segue a vida.
    feeds = [{"url": "https://www.sciencedaily.com/rss/health_medicine/pharmacology.xml", "lang": "🇺🇸"},
             {"url": "https://agencia.fapesp.br/rss/", "lang": "🇧🇷"}]
    noticias = []
    try:
        for fonte in feeds:
            feed = feedparser.parse(fonte["url"])
            if feed.entries:
                for entry in feed.entries[:2]: 
                    noticias.append({"titulo": entry.title, "link": entry.link, "fonte": fonte["lang"], "img": random.choice(IMAGENS_SCIENCE)})
        random.shuffle(noticias)
        return noticias
    except:
        return []