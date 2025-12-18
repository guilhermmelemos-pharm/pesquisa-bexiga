# backend.py
import pandas as pd
from Bio import Entrez
import time
import re
from deep_translator import GoogleTranslator
from datetime import datetime
import feedparser
import random
import constantes as c

# Regex ajustado para pegar compostos químicos e genes, evitando texto comum
REGEX_COMPLEXO = r'\b([A-Z][a-zA-Z0-9\-]*(?:\s[A-Z][a-zA-Z0-9\-]+){0,3})\b'

IMAGENS_SCIENCE = [
    "https://images.unsplash.com/photo-1532094349884-543bc11b234d?auto=format&fit=crop&w=600&q=80", 
    "https://images.unsplash.com/photo-1576086213369-97a306d36557?auto=format&fit=crop&w=600&q=80", 
    "https://images.unsplash.com/photo-1507413245164-6160d8298b31?auto=format&fit=crop&w=600&q=80",
    "https://images.unsplash.com/photo-1581093458791-9f3c3900df4b?auto=format&fit=crop&w=600&q=80",
    "https://images.unsplash.com/photo-1530026405186-ed1f139313f8?auto=format&fit=crop&w=600&q=80",
    "https://images.unsplash.com/photo-1579154204601-01588f351e67?auto=format&fit=crop&w=600&q=80",
]

# --- LISTAS DE BLOQUEIO INTERNAS (HARDCODED PARA PERFORMANCE) ---
BLACKLIST_SOBRENOMES = [
    "wang", "zhang", "li", "liu", "yang", "huang", "wu", "zhou", "xu", "sun", "ma", "zhu", "hu", "guo", "he", "gao", "lin", "luo", "cheng",
    "kim", "lee", "park", "choi", "jeong", "chung", "song", "kang", "jang", "han", "lim",
    "singh", "kumar", "sharma", "patel", "gupta", "mishra", "yadav", "das",
    "smith", "johnson", "williams", "brown", "jones", "garcia", "miller", "davis", "rodriguez", "martinez",
    "chen", "zhao", "yu", "yin", "chang", "ng", "tran", "tan", "wong", "chan", "cai", "ao", "kr", "bm", "nd", "jr", "sr", "et", "al"
]

BLACKLIST_GEOGRAFIA = [
    "china", "usa", "japan", "germany", "uk", "france", "italy", "canada", "australia", "spain", "korea", "india", "brazil",
    "beijing", "shanghai", "guangzhou", "wuhan", "nanjing", "tianjin", "hangzhou",
    "london", "cambridge", "oxford", "boston", "new york", "california", "texas", "arizona", "massachusetts", "pennsylvania", "michigan",
    "tokyo", "osaka", "seoul", "busan", "daegu", "taipei", "hong kong", "singapore",
    "university", "hospital", "institute", "center", "school", "college", "academy", "department", "division", "ministry", "state", "provincial",
    "street", "road", "avenue", "district", "province", "city", "ltd", "inc", "corp", "co", "plc"
]

def validar_termo_inteligente(termo):
    """
    Filtro Heurístico Avançado: Mata nomes, cidades e lixo bibliográfico.
    """
    t_limpo = termo.strip()
    t_lower = t_limpo.lower()
    
    # 1. Tamanho e Estrutura
    if len(t_limpo) < 3: return False
    if t_limpo.isdigit(): return False
    if re.match(r'^CD\d+$', t_limpo): return False # Remove CD44, CD24, etc.

    # 2. Blacklist de Sobrenomes (Detecta "Chen", "Wang Y", "Li X")
    palavras = t_lower.split()
    if any(p in BLACKLIST_SOBRENOMES for p in palavras):
        # Se for só o sobrenome ou sobrenome + inicial, mata.
        # Exceção: Termos científicos que coincidem com nomes (raro em 2 palavras, mas possível)
        if len(palavras) <= 2: return False

    # 3. Blacklist Geográfica e Institucional
    if any(p in BLACKLIST_GEOGRAFIA for p in palavras):
        return False

    # 4. Regra de Siglas Curtas (ex: BPH, BJOG) - Exceções VIPs
    whitelist_siglas = ["PFAS", "TMAO", "NMN", "CBD", "CBG", "LPS", "ATP", "NAD", "RNA", "DNA", "ROS", "NO"]
    if t_limpo.isupper() and len(t_limpo) <= 5 and t_limpo not in whitelist_siglas: 
        return False
    
    # 5. Regra de Sufixos Corporativos
    if t_lower.endswith((" therapeutics", " bio", " pharma", " oncology", " inc", " co", " ltd", " group", " sci", " rep", " med", " biol", " chem")):
        return False

    return True

def buscar_alvos_emergentes_pubmed(orgao_alvo, email):
    """
    Busca profunda no PubMed exigindo contexto molecular para evitar lixo.
    """
    if not orgao_alvo or not email: return []
    Entrez.email = email
    ano_fim = datetime.now().year
    ano_inicio = ano_fim - 3 # Aumentei para 3 anos para pegar mais robustez
    
    # QUERY REFINADA: Exige que o termo alvo esteja perto de palavras de ação biológica
    # Isso evita pegar endereços, pois endereços não "ativam", "inibem" ou "regulam".
    query = (
        f"({orgao_alvo}[Title/Abstract]) AND ("
        f"\"molecular mechanism\" OR \"signaling pathway\" OR \"therapeutic target\" OR "
        f"\"gene expression\" OR \"upregulation\" OR \"downregulation\" OR \"activation\" OR "
        f"\"inhibition\" OR \"novel receptor\" OR \"metabolite\" OR \"microplastics\" OR \"TMAO\")"
        f" AND (\"{ano_inicio}\"[Date - Publication] : \"{ano_fim}\"[Date - Publication])"
        f" NOT (Review[Publication Type])" 
    )
    
    try:
        # Busca mais IDs para poder filtrar bastante
        handle = Entrez.esearch(db="pubmed", term=query, retmax=80, sort="relevance")
        ids = Entrez.read(handle)["IdList"]
        if not ids: return []
        
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
        texto = handle.read()
        
        candidatos = re.findall(REGEX_COMPLEXO, texto)
        
        termos_limpos = []
        blacklist_lower = [x.lower() for x in c.BLACKLIST_GERAL]
        
        # Contagem de frequência para priorizar termos que aparecem mais de uma vez (Sinal vs Ruído)
        frequencia = {}
        for t in candidatos:
            t = t.strip(" .,();:")
            if t in frequencia: frequencia[t] += 1
            else: frequencia[t] = 1

        for t, freq in frequencia.items():
            t_lower = t.lower()
            
            # FILTRAGEM EM CASCATA
            if any(bad in t_lower for bad in blacklist_lower): continue
            if not validar_termo_inteligente(t): continue
            if any(bad in t_lower.split() for bad in ["were", "was", "rates", "group", "study", "analysis", "doi", "vol", "no", "fig"]): continue
            
            # Regra de Ouro: Termos muito curtos que aparecem 1x só geralmente são ruído ou iniciais.
            # Termos bons tendem a repetir ou são compostos.
            if len(t) < 5 and freq < 2 and t not in ["TMAO", "PFAS", "NMN"]: continue

            if "-" in t or any(x.isupper() for x in t[1:]) or " " in t or t[0].isupper():
                termos_limpos.append(t)

        # Retorna os top 40 ordenados por frequência (os mais citados primeiro)
        termos_ordenados = sorted(termos_limpos, key=lambda x: frequencia[x], reverse=True)
        return termos_ordenados[:40] 
    except:
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
    texto_bruto = match.group(2).strip() if match else ". ".join(abstract_text.split(". ")[-3:])
    texto_final = texto_bruto[:600] # Aumentei um pouco
    try:
        return translator.translate(texto_final) + "..."
    except:
        return texto_final + "..."

def buscar_resumos_detalhados(termo_alvo_especifico, termo_orgao_interesse, email, y_start, y_end, lang_target, limit=5):
    if not email: return []
    Entrez.email = email
    query = f"(\"{termo_alvo_especifico}\") AND ({termo_orgao_interesse}) AND {y_start}:{y_end}[DP]"
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
    feeds = [{"url": "https://www.sciencedaily.com/rss/health_medicine/pharmacology.xml", "lang": "🇺🇸"},
             {"url": "https://agencia.fapesp.br/rss/", "lang": "🇧🇷"},
             {"url": "https://www.nature.com/nature.rss", "lang": "🌍"}]
    noticias = []
    translator = GoogleTranslator(source='auto', target=lang_code)
    for fonte in feeds:
        try:
            feed = feedparser.parse(fonte["url"])
            for entry in feed.entries[:2]: 
                img_url = random.choice(IMAGENS_SCIENCE)
                titulo = entry.title
                if lang_code == 'pt' and fonte["lang"] != "🇧🇷":
                    try: titulo = translator.translate(titulo)
                    except: pass
                noticias.append({"titulo": titulo, "link": entry.link, "fonte": feed.feed.title.split("-")[0].strip()[:20], "img": img_url, "bandeira": fonte["lang"]})
        except: continue
    random.shuffle(noticias)
    return noticias