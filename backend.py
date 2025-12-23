import streamlit as st
from Bio import Entrez
import google.generativeai as genai
from tenacity import retry, stop_after_attempt, wait_exponential
import re
from collections import Counter
import time

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br" 

def analisar_abstract_com_ia(titulo, abstract, api_key, lang='pt'):
    if not api_key: return "⚠️ IA não ativada"
    try:
        genai.configure(api_key=api_key)
        idioma = "Português" if lang == 'pt' else "Inglês"
        abs_input = abstract[:3000] if (abstract and len(abstract) > 30) else "Use apenas o título."
        
        prompt = f"""Como PhD em Farmacologia, analise:
        TÍTULO: {titulo} | RESUMO: {abs_input}
        FORMATO: Alvo → Fármaco → Efeito (Contextualizado ao Título).
        REGRAS: Máx 25 palavras. Resposta única e técnica. Idioma: {idioma}."""

        for mod in ['gemini-1.5-flash-8b', 'gemini-1.5-flash', 'gemini-1.5-pro']:
            try:
                model = genai.GenerativeModel(mod)
                response = model.generate_content(prompt, generation_config={"temperature": 0.7})
                return response.text.strip()
            except:
                time.sleep(1)
                continue 
        return f"💡 IA Ocupada: {titulo[:40]}..."
    except Exception as e:
        return f"❌ Erro: {str(e)[:30]}"

@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
def _fetch_pubmed_count(query):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    record = Entrez.read(handle); handle.close()
    return int(record["Count"])

@st.cache_data(ttl=86400, show_spinner=False)
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    query = f"({termo}) AND ({contexto}) AND ({ano_ini}:{ano_fim}[Date - Publication]) AND (NOT Review[pt])"
    try: return _fetch_pubmed_count(query)
    except: return 0

@st.cache_data(ttl=86400, show_spinner=False)
def buscar_resumos_detalhados(termo, orgao, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    query = f"({termo}) AND ({orgao}) AND ({ano_ini}:{ano_fim}[Date - Publication]) AND (NOT Review[pt])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5, sort="relevance")
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        artigos = []
        for raw in dados.split("\n\nPMID-"):
            tit, abs_txt, pmid, last_tag = "", "", "", ""
            for line in raw.split("\n"):
                if line.strip().isdigit() and not pmid: pmid = line.strip()
                if re.match(r'^TI\s+-', line): tit = re.sub(r'^TI\s+-\s+', '', line).strip(); last_tag = "TI"
                elif line.startswith("      ") and last_tag == "TI": tit += " " + line.strip()
                elif re.match(r'^AB\s+-', line): abs_txt = re.sub(r'^AB\s+-\s+', '', line).strip(); last_tag = "AB"
                elif line.startswith("      ") and last_tag == "AB": abs_txt += " " + line.strip()
            if tit: artigos.append({"Title": tit, "Resumo_Original": abs_txt, "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return artigos
    except: return []

@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email):
    if email: Entrez.email = email
    query = f"({termo_base}) AND (receptor OR pathway OR channel) AND (2024:2026[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=50)
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read().upper(); handle.close()
        
        # BLACKLIST AGRESSIVA (Filtra resíduos metodológicos e geográficos)
        blacklist = {
            "THE", "AND", "FOR", "NOT", "BUT", "WITH", "FROM", "AFTER", "BEFORE", "DURING", "BETWEEN",
            "STUDY", "GROUP", "DATA", "ANALYSIS", "RESULTS", "METHODS", "CONCLUSION", "AIMS", "SUMMARY",
            "HIGH", "LOW", "INCREASED", "DECREASED", "LEVELS", "EXPRESSION", "PATHWAY", "TARGET",
            "ROLE", "EFFECT", "IMPACT", "ACTION", "SYSTEM", "FUNCTION", "ACTIVITY", "POTENTIAL",
            "TREATED", "INDUCED", "MEDIATED", "ASSOCIATED", "OBSERVED", "COMPARED", "PERFORMED",
            "URINARY", "BLADDER", "URETHRA", "KIDNEY", "LIVER", "HEART", "BRAIN", "LUNG", "MUSCLE",
            "CANCER", "TUMOR", "DISEASE", "SYNDROME", "PAIN", "INFLAMMATION", "INFECTION", "INJURY",
            "PATIENT", "CLINICAL", "TRIAL", "THERAPY", "TREATMENT", "DRUG", "DOSE", "CONTROL",
            "MALE", "FEMALE", "ADULT", "CHILD", "AGE", "YEAR", "MONTH", "DAY", "TIME",
            "UNITED", "STATES", "CHINA", "JAPAN", "BRAZIL", "EUROPE", "ASIAN", "AMERICAN",
            "FIGURE", "TABLE", "PMID", "PMC", "DOI", "ISSN", "URL", "HTTP", "WWW", "TYPE", "CLASS", 
            "RECEPTOR", "CHANNEL", "PROTEIN", "GENE", "FACTOR", "RESPONSE", "CELLS", "TISSUE", "MODEL",
            "USING", "REVIEW", "MECHANISM", "SIGNALING", "NOVEL", "NEW", "ACUTE", "CHRONIC", "HUMAN",
            "MOUSE", "RAT", "SIGNIFICANT", "STATISTICAL", "VALUE", "MEAN", "RATE", "RATIO", "NORMAL",
            "POSITIVE", "NEGATIVE", "ACTIVE", "INACTIVE", "STABLE", "UNSTABLE", "TOTAL", "CASE", "REPORT"
        }
        
        # Filtro de Unidades e Termos Curtos Irrelevantes
        unidades = {"MMHG", "KPA", "MIN", "SEC", "HRS", "ML", "MG", "KG", "NM", "UM", "MM"}
        
        encontrados = re.findall(r'\b[A-Z][A-Z0-9-]{2,8}\b', dados)
        candidatos = []
        for t in encontrados:
            # Só aceita se: não estiver na blacklist, não for unidade, for longo o suficiente e não for o termo base
            if t not in blacklist and t not in unidades and len(t) >= 3 and t not in termo_base.upper():
                # Evita siglas de países e termos muito comuns de 3 letras que sobraram
                if t not in {"USA", "PRC", "UK", "EU", "APP", "ALL", "WAS", "ARE", "HAS"}:
                    candidatos.append(t)
        
        return [t for t, count in Counter(candidatos).most_common(7)]
    except: return []

def buscar_todas_noticias(lang='pt'):
    try:
        # Radar de Ciência Global: Nature, Science, Cell e Farmacologia Geral
        handle = Entrez.esearch(db="pubmed", term="(pharmacology[Filter]) AND (2025[Date - Publication])", retmax=3, sort="pub_date")
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        news = []
        for art in dados.split("\n\nPMID-"):
            tit, pmid = "", ""
            for line in art.split("\n"):
                if re.match(r'^TI\s+-', line): tit = re.sub(r'^TI\s+-\s+', '', line).strip()
                if line.strip().isdigit() and not pmid: pmid = line.strip()
            if tit: news.append({"titulo": tit, "fonte": "Science Frontier", "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return news
    except: return []
