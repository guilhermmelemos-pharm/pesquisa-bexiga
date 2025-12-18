# backend.py
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

# --- CONFIGURAÇÃO NLTK (BAIXA DICIONÁRIOS APENAS UMA VEZ) ---
try:
    nltk.data.find('corpora/words')
except LookupError:
    nltk.download('words')
    nltk.download('stopwords')

# Conjuntos para verificação rápida
ENGLISH_WORDS = set(word.lower() for word in words.words())
STOP_WORDS = set(stopwords.words('english'))

# Regex ajustado
REGEX_COMPLEXO = r'\b([A-Z][a-zA-Z0-9\-]*(?:\s[A-Z][a-zA-Z0-9\-]+){0,3})\b'

IMAGENS_SCIENCE = [
    "https://images.unsplash.com/photo-1532094349884-543bc11b234d?auto=format&fit=crop&w=600&q=80", 
    "https://images.unsplash.com/photo-1576086213369-97a306d36557?auto=format&fit=crop&w=600&q=80", 
    "https://images.unsplash.com/photo-1507413245164-6160d8298b31?auto=format&fit=crop&w=600&q=80",
    "https://images.unsplash.com/photo-1581093458791-9f3c3900df4b?auto=format&fit=crop&w=600&q=80",
    "https://images.unsplash.com/photo-1530026405186-ed1f139313f8?auto=format&fit=crop&w=600&q=80",
    "https://images.unsplash.com/photo-1579154204601-01588f351e67?auto=format&fit=crop&w=600&q=80",
]

# --- LISTA NEGRA GEOGRÁFICA & INSTITUCIONAL (Onde o lixo se esconde) ---
BLACKLIST_LOCAIS = [
    # Países e Regiões
    "china", "usa", "japan", "germany", "uk", "france", "italy", "canada", "australia", "spain", "korea", "india", "brazil",
    "united states", "america", "republic", "taiwan", "people", "bca",
    # Cidades Chinesas/Asiáticas (MUITO COMUM EM META-DADOS)
    "beijing", "shanghai", "guangzhou", "wuhan", "nanjing", "tianjin", "hangzhou", 
    "chongqing", "guangdong", "jiangsu", "zhengzhou", "shijiazhuang", "lanzhou", 
    "nanchang", "changsha", "zheng", "shenzhen", "kaohsiung", "hong kong", "singapore",
    # Cidades Ocidentais e Estados
    "london", "cambridge", "oxford", "boston", "new york", "california", "texas", "arizona", 
    "massachusetts", "pennsylvania", "michigan", "pittsburgh", "heidelberg", "paris", "san francisco",
    # Instituições e Departamentos
    "university", "hospital", "institute", "center", "school", "college", "academy", 
    "department", "division", "ministry", "state", "provincial", "laboratory", "faculty",
    "key laboratory", "national institutes", "foundation", "association", "society", "board",
    "education", "health", "engineering", "sciences", "technology", "agriculture"
]

BLACKLIST_ACADEMICA = [
    # Seções do Artigo
    "abstract", "introduction", "background", "methods", "results", "discussion", "conclusion",
    "declarations", "declaration", "conflict", "interest", "funding", "acknowledgements",
    "availability", "contributed", "corresponding", "author", "editor", "received", "accepted",
    # Termos de Indexação e Publicação
    "medline", "indexed", "electronic", "epub", "print", "pmid", "doi", "issn", "isbn",
    "vol", "issue", "suppl", "fig", "table", "copyright", "license", "open access",
    # Áreas Genéricas
    "oncology", "pathology", "urology", "gastroenterology", "hepatology", "ophthalmology", 
    "cardiology", "neurology", "cancer", "biology", "chemistry", "medicine"
]

def validar_termo_bioquimico(termo):
    """
    Algoritmo v1.6: Bio-Validação com NLTK e Filtro Geográfico Agressivo.
    """
    t_limpo = termo.strip()
    t_lower = t_limpo.lower()
    palavras = t_lower.split()
    
    # 1. WHITELIST VIP (Esses passam sempre)
    whitelist = ["tmao", "pfas", "nmn", "cbd", "cbg", "lps", "atp", "nad", "dna", "rna", "gpr", "rock", "mtor", "ampk"]
    if t_lower in whitelist: return True # Exato
    if any(vip in t_lower for vip in whitelist) and len(t_limpo) > 3: return True # Parcial

    # 2. FILTRO DE LIXO BÁSICO
    if len(t_limpo) < 3: return False
    if t_limpo.isdigit(): return False
    
    # 3. FILTRO DE DICIONÁRIO NLTK (Mata "However", "Moreover", "The")
    # Se é palavra comum e não tem sufixo bio, morre.
    sufixos_bio = ("ase", "in", "or", "ol", "ide", "ate", "one", "ics", "ome", "rna", "dna", "ab")
    
    # Verifica se é uma stopword (The, A, An, Our, These) - Mata direto
    if t_lower in STOP_WORDS: return False
    
    # Verifica se é palavra comum (University, Electronic, Indexed)
    if t_lower in ENGLISH_WORDS and not t_lower.endswith(sufixos_bio):
        return False 

    # 4. FILTRO GEOGRÁFICO E INSTITUCIONAL (NOVO)
    if any(loc in t_lower for loc in BLACKLIST_LOCAIS):
        return False

    # 5. FILTRO ACADÊMICO (NOVO)
    if any(acad in t_lower for acad in BLACKLIST_ACADEMICA):
        return False

    # 6. FILTRO DE PADRÃO DE AUTOR (Sobrenome + Iniciais)
    if re.search(r' [A-Z]{1,2}$', t_limpo): return False

    # 7. FILTRO DE MESES (Jan, Feb...)
    meses = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]
    if t_lower[:3] in meses and len(t_limpo) < 5: return False

    # 8. CRITÉRIO DE ACEITAÇÃO FINAL (BIO-SCORE)
    # Se passou por tudo isso, precisa parecer ciência.
    tem_numero = any(char.isdigit() for char in t_limpo)
    tem_hifen = "-" in t_limpo
    tem_sufixo_bio = t_lower.endswith(sufixos_bio)
    tem_grego = any(g in t_lower for g in ["alpha", "beta", "gamma", "delta", "kappa"])
    eh_composto = len(palavras) > 1
    
    # Se for palavra única simples (sem numero, hifen, sufixo), provavelmente é lixo que o dicionário não pegou
    if not eh_composto and not tem_numero and not tem_hifen and not tem_sufixo_bio and not tem_grego:
        # Última chance: está na whitelist interna de compostos?
        if t_lower not in ["exosomes", "ferroptosis", "pyroptosis", "autophagy", "mitophagy", "senolytics"]:
            return False

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
        handle = Entrez.esearch(db="pubmed", term=query, retmax=100, sort="relevance")
        ids = Entrez.read(handle)["IdList"]
        if not ids: return []
        
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
        texto = handle.read()
        
        candidatos = re.findall(REGEX_COMPLEXO, texto)
        
        termos_limpos = []
        blacklist_lower = [x.lower() for x in c.BLACKLIST_GERAL]
        
        frequencia = {}
        for t in candidatos:
            t = t.strip(" .,();:")
            if t in frequencia: frequencia[t] += 1
            else: frequencia[t] = 1

        for t, freq in frequencia.items():
            t_lower = t.lower()
            
            if any(bad in t_lower for bad in blacklist_lower): continue
            
            # CHAMA O NOVO VALIDADOR V1.6
            if not validar_termo_bioquimico(t): continue
            
            termos_limpos.append(t)

        termos_ordenados = sorted(termos_limpos, key=lambda x: frequencia[x], reverse=True)
        
        seen = set()
        resultado_final = []
        for item in termos_ordenados:
            if item.lower() not in seen:
                seen.add(item.lower())
                resultado_final.append(item)
                
        return resultado_final[:30] 
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
    texto_final = texto_bruto[:600]
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