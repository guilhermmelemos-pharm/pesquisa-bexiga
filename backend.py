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

def buscar_alvos_emergentes_pubmed(orgao_alvo, email):
    if not orgao_alvo or not email:
        return []
    Entrez.email = email
    ano_fim = datetime.now().year
    ano_inicio = ano_fim - 1
    
    query = (
        f"({orgao_alvo}) AND ('GPCR' OR 'ion channel' OR 'transporter' OR 'orphan receptor' OR "
        f"'gasotransmitter' OR 'signaling pathway') AND "
        f"(\"{ano_inicio}\"[Date - Publication] : \"{ano_fim}\"[Date - Publication])"
    )
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=35, sort="relevance")
        record = Entrez.read(handle)
        ids = record["IdList"]
        if not ids: return []
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
        texto = handle.read()
        encontrados = re.findall(r'\b[A-Z]{2,6}[0-9]{0,4}\b', texto)
        blacklist = ["DNA", "RNA", "USA", "NCBI", "NIH", "ATP", "AMP", "GDP", "COVID", "SARS", "PMID", "DOI", 
                     "FAPESP", "UNIFESP", "HPLC", "PCR", "ANOVA", "SD", "SEM", "GROUP", "MEAN"]
        return sorted(list(set([t for t in encontrados if t.upper() not in blacklist and len(t) > 2])))
    except:
        return []

def buscar_todas_noticias(lang_code):
    feeds = [
        {"url": "https://www.sciencedaily.com/rss/health_medicine/pharmacology.xml", "lang": "🇺🇸"},
        {"url": "https://agencia.fapesp.br/rss/", "lang": "🇧🇷"},
    ]
    noticias = []
    backups = ["https://images.unsplash.com/photo-1532094349884-543bc11b234d?w=400&h=250&fit=crop"]
    translator = GoogleTranslator(source='auto', target=lang_code)
    for fonte in feeds:
        try:
            feed = feedparser.parse(fonte["url"])
            for entry in feed.entries[:3]:
                img_url = random.choice(backups)
                titulo = entry.title
                if lang_code == 'pt' and fonte["lang"] != "🇧🇷":
                    try: titulo = translator.translate(titulo)
                    except: pass
                noticias.append({
                    "titulo": titulo, "link": entry.link,
                    "fonte": feed.feed.title.split("-")[0].strip()[:20], 
                    "img": img_url, "bandeira": fonte["lang"]
                })
        except: continue
    random.shuffle(noticias)
    return noticias

def consultar_pubmed_count(termo_farmaco, termo_orgao, email, y_start, y_end):
    if not email: return -1
    Entrez.email = email
    time.sleep(0.1) # Segurança NCBI
    query = f"({termo_farmaco}) AND ({termo_orgao}) AND {y_start}:{y_end}[DP]" if termo_orgao else f"({termo_farmaco}) AND {y_start}:{y_end}[DP]"
    try: return int(Entrez.read(Entrez.esearch(db="pubmed", term=query, retmax=0))["Count"])
    except: return -1

def extrair_conclusao(abstract_text, lang_target):
    if not abstract_text: return "Abstract not available."
    match = re.search(r'(Conclusion|Conclusions|In conclusion|Summary|Results suggest that)(.*)', abstract_text, re.IGNORECASE | re.DOTALL)
    texto_final = match.group(2).strip()[:400] if match else abstract_text[-400:]
    return ("🇧🇷 " if lang_target=='pt' else "🇺🇸 ") + GoogleTranslator(source='auto', target=lang_target).translate(texto_final) + "..."

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
                art_data["Resumo_IA"] = extrair_conclusao(art_data["Abstract"], lang_target)
                artigos.append(art_data)
        return artigos
    except: return []