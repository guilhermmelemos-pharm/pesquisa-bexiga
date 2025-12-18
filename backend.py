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

REGEX_COMPLEXO = r'\b([A-Z][a-zA-Z0-9\-]*(?:\s[A-Z][a-zA-Z0-9\-]+){0,3})\b'

IMAGENS_SCIENCE = [
    "https://images.unsplash.com/photo-1532094349884-543bc11b234d?auto=format&fit=crop&w=600&q=80", 
    "https://images.unsplash.com/photo-1576086213369-97a306d36557?auto=format&fit=crop&w=600&q=80", 
    "https://images.unsplash.com/photo-1507413245164-6160d8298b31?auto=format&fit=crop&w=600&q=80",
    "https://images.unsplash.com/photo-1581093458791-9f3c3900df4b?auto=format&fit=crop&w=600&q=80",
    "https://images.unsplash.com/photo-1530026405186-ed1f139313f8?auto=format&fit=crop&w=600&q=80",
    "https://images.unsplash.com/photo-1579154204601-01588f351e67?auto=format&fit=crop&w=600&q=80",
]

def buscar_alvos_emergentes_pubmed(orgao_alvo, email):
    if not orgao_alvo or not email: return []
    Entrez.email = email
    ano_fim = datetime.now().year
    ano_inicio = ano_fim - 2 
    
    query = (
        f"({orgao_alvo}) AND ("
        f"'novel target' OR 'emerging' OR 'receptor' OR 'pathway' OR 'metabolite' OR "
        f"'progenitor' OR 'stem cell' OR 'embryonic' OR 'organoid' OR 'developmental' OR "
        f"'mechanism')"
        f" AND (\"{ano_inicio}\"[Date - Publication] : \"{ano_fim}\"[Date - Publication])"
    )
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=40, sort="relevance")
        ids = Entrez.read(handle)["IdList"]
        if not ids: return []
        
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
        texto = handle.read()
        
        candidatos = re.findall(REGEX_COMPLEXO, texto)
        
        termos_limpos = []
        blacklist_lower = [x.lower() for x in c.BLACKLIST_GERAL]
        
        for t in candidatos:
            t = t.strip(" .,()")
            t_lower = t.lower()
            
            if len(t) < 3 or len(t) > 50: continue
            
            # FILTRO PARCIAL DE BLACKLIST (AQUI ESTÁ A CORREÇÃO)
            # Se QUALQUER termo da blacklist estiver contido no termo candidato, descarta.
            if any(bad in t_lower for bad in blacklist_lower): continue
            
            # Filtro de palavras proibidas soltas
            if any(bad in t_lower.split() for bad in ["were", "was", "rates", "group", "study", "analysis"]): continue
            
            if "-" in t or any(x.isupper() for x in t[1:]) or " " in t or t[0].isupper():
                termos_limpos.append(t)

        return sorted(list(set(termos_limpos)))[:25] 
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
    texto_final = texto_bruto[:500] 
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