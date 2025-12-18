# backend.py
import pandas as pd
from Bio import Entrez
import time
import re
from deep_translator import GoogleTranslator
from datetime import datetime
import feedparser
import random

# Regex melhorado para capturar termos compostos (ex: "Alpha-lipoic acid")
REGEX_COMPLEXO = r'\b([A-Z][a-zA-Z0-9\-]*(?:\s[a-zA-Z0-9\-]+){0,2})\b'

def buscar_alvos_emergentes_pubmed(orgao_alvo, email):
    if not orgao_alvo or not email: return []
    Entrez.email = email
    ano_fim = datetime.now().year
    ano_inicio = ano_fim - 2 # Aumentei para 2 anos para pegar mais contexto
    
    # Query focada em NOVIDADE e COMPLEXIDADE
    query = (
        f"({orgao_alvo}) AND ("
        f"'novel therapeutic target' OR 'emerging metabolite' OR 'signaling pathway' OR "
        f"'fatty acid' OR 'receptor agonist' OR 'orphan receptor')"
        f" AND (\"{ano_inicio}\"[Date - Publication] : \"{ano_fim}\"[Date - Publication])"
    )
    
    try:
        # Busca mais ampla
        handle = Entrez.esearch(db="pubmed", term=query, retmax=50, sort="relevance")
        ids = Entrez.read(handle)["IdList"]
        if not ids: return []
        
        # Leitura dos abstracts
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
        texto = handle.read()
        
        # Mineração de termos compostos
        candidatos = re.findall(REGEX_COMPLEXO, texto)
        
        # Filtro de qualidade (Blacklist)
        blacklist = ["The", "A", "An", "In", "Of", "For", "And", "With", "By", "To", "At", "On",
                     "PubMed", "Abstract", "Introduction", "Results", "Conclusion", "Department",
                     "University", "Study", "Group", "Data", "Analysis", "Significant", "Control",
                     "Figure", "Table", "Method", "Patient", "Rat", "Mouse", "Mice", "Cells",
                     "Level", "Effect", "Role", "High", "Low", "Expression", "Treatment"]
                     
        termos_limpos = []
        for t in candidatos:
            t = t.strip(" .,()")
            # Só aceita se tiver mais de 3 letras, não for blacklist e parecer termo científico
            if (len(t) > 3) and (t.split()[0] not in blacklist) and (t not in blacklist):
                # Prefere termos com hífens ou letras maiúsculas no meio (ex: mTOR, NF-kB) ou compostos
                if "-" in t or any(x.isupper() for x in t[1:]) or " " in t:
                    termos_limpos.append(t)

        return sorted(list(set(termos_limpos)))[:20] # Retorna top 20 para não poluir
    except:
        return []

def consultar_pubmed_count(termo_farmaco, termo_orgao, email, y_start, y_end):
    if not email: return 0
    Entrez.email = email
    time.sleep(0.1) 
    # Query precisa
    query = f"(\"{termo_farmaco}\") AND ({termo_orgao}) AND {y_start}:{y_end}[DP]"
    try: return int(Entrez.read(Entrez.esearch(db="pubmed", term=query, retmax=0))["Count"])
    except: return 0

def buscar_todas_noticias(lang_code):
    # (Mantido igual ao anterior, apenas para garantir integridade)
    feeds = [{"url": "https://www.sciencedaily.com/rss/health_medicine/pharmacology.xml", "lang": "🇺🇸"}]
    noticias = []
    try:
        feed = feedparser.parse(feeds[0]["url"])
        for entry in feed.entries[:3]:
            noticias.append({"titulo": entry.title, "link": entry.link, "img": "https://source.unsplash.com/random/400x200/?science", "fonte": "ScienceDaily", "bandeira": "🇺🇸"})
    except: pass
    return noticias