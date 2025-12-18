# backend.py
# Lógica de negócio, chamadas à API e processamento.

import pandas as pd
from Bio import Entrez
import time
import re
from deep_translator import GoogleTranslator
from datetime import datetime
import feedparser
import random

# Regex para capturar termos compostos
REGEX_COMPLEXO = r'\b([A-Z][a-zA-Z0-9\-]*(?:\s[a-zA-Z0-9\-]+){0,2})\b'

# Banco de imagens científicas confiáveis (Unsplash URLs diretas)
IMAGENS_SCIENCE = [
    "https://images.unsplash.com/photo-1532094349884-543bc11b234d?auto=format&fit=crop&w=600&q=80", # Lab
    "https://images.unsplash.com/photo-1576086213369-97a306d36557?auto=format&fit=crop&w=600&q=80", # Células
    "https://images.unsplash.com/photo-1507413245164-6160d8298b31?auto=format&fit=crop&w=600&q=80", # DNA
    "https://images.unsplash.com/photo-1581093458791-9f3c3900df4b?auto=format&fit=crop&w=600&q=80", # Microscópio
    "https://images.unsplash.com/photo-1530026405186-ed1f139313f8?auto=format&fit=crop&w=600&q=80", # Moléculas
    "https://images.unsplash.com/photo-1579154204601-01588f351e67?auto=format&fit=crop&w=600&q=80", # Pipeta
]

def buscar_alvos_emergentes_pubmed(orgao_alvo, email):
    if not orgao_alvo or not email: return []
    Entrez.email = email
    ano_fim = datetime.now().year
    ano_inicio = ano_fim - 2 
    
    query = (
        f"({orgao_alvo}) AND ("
        f"'novel therapeutic target' OR 'emerging metabolite' OR 'signaling pathway' OR "
        f"'fatty acid' OR 'receptor agonist' OR 'orphan receptor')"
        f" AND (\"{ano_inicio}\"[Date - Publication] : \"{ano_fim}\"[Date - Publication])"
    )
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=50, sort="relevance")
        ids = Entrez.read(handle)["IdList"]
        if not ids: return []
        
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
        texto = handle.read()
        
        candidatos = re.findall(REGEX_COMPLEXO, texto)
        blacklist = ["The", "A", "An", "In", "Of", "For", "And", "With", "By", "To", "At", "On",
                     "PubMed", "Abstract", "Introduction", "Results", "Conclusion", "Department",
                     "University", "Study", "Group", "Data", "Analysis", "Significant", "Control",
                     "Figure", "Table", "Method", "Patient", "Rat", "Mouse", "Mice", "Cells",
                     "Level", "Effect", "Role", "High", "Low", "Expression", "Treatment"]
                     
        termos_limpos = []
        for t in candidatos:
            t = t.strip(" .,()")
            if (len(t) > 3) and (t.split()[0] not in blacklist) and (t not in blacklist):
                if "-" in t or any(x.isupper() for x in t[1:]) or " " in t:
                    termos_limpos.append(t)

        return sorted(list(set(termos_limpos)))[:20] 
    except:
        return []

def consultar_pubmed_count(termo_farmaco, termo_orgao, email, y_start, y_end):
    if not email: return 0
    Entrez.email = email
    time.sleep(0.1) 
    query = f"(\"{termo_farmaco}\") AND ({termo_orgao}) AND {y_start}:{y_end}[DP]" if termo_orgao else f"(\"{termo_farmaco}\") AND {y_start}:{y_end}[DP]"
    try: return int(Entrez.read(Entrez.esearch(db="pubmed", term=query, retmax=0))["Count"])
    except: return 0

def buscar_todas_noticias(lang_code):
    feeds = [
        {"url": "https://www.sciencedaily.com/rss/health_medicine/pharmacology.xml", "lang": "🇺🇸"},
        {"url": "https://agencia.fapesp.br/rss/", "lang": "🇧🇷"},
        {"url": "https://www.nature.com/nature.rss", "lang": "🌍"}
    ]
    noticias = []
    translator = GoogleTranslator(source='auto', target=lang_code)
    
    for fonte in feeds:
        try:
            feed = feedparser.parse(fonte["url"])
            for entry in feed.entries[:2]: # Pega 2 de cada fonte
                # Tenta achar imagem no feed, senão usa aleatória científica
                img_url = random.choice(IMAGENS_SCIENCE)
                
                # Tradução do título se necessário
                titulo = entry.title
                if lang_code == 'pt' and fonte["lang"] != "🇧🇷":
                    try: titulo = translator.translate(titulo)
                    except: pass
                
                noticias.append({
                    "titulo": titulo, 
                    "link": entry.link,
                    "fonte": feed.feed.title.split("-")[0].strip()[:20], 
                    "img": img_url, 
                    "bandeira": fonte["lang"]
                })
        except: continue
    random.shuffle(noticias)
    return noticias

def buscar_resumos_detalhados(termo_alvo_especifico, termo_orgao_interesse, email, y_start, y_end, lang_target, limit=5):
    if not email: return []
    query = f"({termo_alvo_especifico}) AND ({termo_orgao_interesse}) AND {y_start}:{y_end}[DP]"
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
                elif tag=="AB": art_data["Abstract"]=content
            if art_data["PMID"]!="N/A":
                # Tradução simplificada da conclusão
                art_data["Resumo_IA"] = "Tradução indisponível no momento." # Placeholder para evitar erros de tradução em massa
                artigos.append(art_data)
        return artigos
    except: return []