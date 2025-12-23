import streamlit as st
from Bio import Entrez
import google.generativeai as genai
from tenacity import retry, stop_after_attempt, wait_exponential
import re
from collections import Counter
import time

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br" 

# --- IA: MODO INVESTIGAÇÃO SÊNIOR ---
def analisar_abstract_com_ia(titulo, abstract, api_key, lang='pt'):
    if not api_key: return "⚠️ IA não ativada"
    try:
        genai.configure(api_key=api_key)
        idioma = "Português" if lang == 'pt' else "Inglês"
        abs_input = abstract[:3000] if (abstract and len(abstract) > 30) else "Resumo incompleto. Use o título."
        
        prompt = f"""Como PhD em Farmacologia, analise:
        TÍTULO: {titulo} | RESUMO: {abs_input}
        TAREFA: Identifique o Alvo e o Fármaco. Descreva o efeito funcional único.
        FORMATO: Alvo → Fármaco → Efeito (Contextualizado ao Título).
        REGRAS: Máx 25 palavras. Proibido respostas genéricas. Idioma: {idioma}."""

        for mod in ['gemini-1.5-flash-8b', 'gemini-1.5-flash', 'gemini-1.5-pro']:
            try:
                model = genai.GenerativeModel(mod)
                response = model.generate_content(prompt, generation_config={"temperature": 0.7})
                return response.text.strip()
            except Exception as e:
                if "429" in str(e): time.sleep(2)
                continue 
        return f"💡 IA Ocupada: {titulo[:50]}..."
    except Exception as e:
        return f"❌ Erro Crítico: {str(e)[:40]}"

# --- BUSCA PUBMED ---
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

# --- MINERAÇÃO (BLACKLIST v3.0) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email):
    if email: Entrez.email = email
    query = f"({termo_base}) AND (receptor OR pathway OR channel) AND (2024:2026[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=50)
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read().upper(); handle.close()
        
        # --- BLACKLIST v4.0 (Agressiva: Limpeza Total) ---
# --- BLACKLIST v5.0 (Filtro Farmacológico Puro) ---
        blacklist = {
            # Ruídos que você identificou (CBS, FAU, BRISTOL, etc.)
            "CBS", "FAU", "AID", "BRISTOL", "COMPANY", "INC", "CORP", "LTD", "PHST", "AUID", "ORCID",
            "CANCER", "UROLOGY", "ONCOLOGY", "CLINICAL", "RESEARCH", "MEDICAL", "MEDICINE", "HOSPITAL",
            "UNIVERSITY", "INSTITUTE", "LABORATORY", "CENTER", "DEPT", "SCHOOL", "FOUNDATION",
            
            # Geográfico e Institucional
            "USA", "PRC", "UK", "EU", "CHINA", "NANJING", "FREIBURG", "UNIFESP", "BRAZIL", "SHANGHAI",
            
            # Conectivos e Verbos que a IA "pesca"
            "WERE", "THAT", "THIS", "THESE", "THOSE", "WHICH", "WHEN", "WHERE", "ALSO", "THAN", "BOTH",
            "UPON", "ONLY", "BEEN", "SOME", "COULD", "WELL", "VERY", "FROM", "INTO", "WITH", "NOT",
            "BUT", "WAS", "ARE", "HAS", "HAD", "ALL", "BEING", "FOR", "AND", "THE", "ABOUT",
            
            # Resíduos Metodológicos e Acadêmicos
            "STUDY", "GROUP", "DATA", "ANALYSIS", "RESULTS", "METHODS", "CONCLUSION", "AIMS", "SUMMARY",
            "PMID", "PMC", "DOI", "ISSN", "URL", "HTTP", "WWW", "PUBLISHED", "COPYRIGHT", "LICENSE",
            "FIGURE", "TABLE", "TYPE", "CLASS", "NORMAL", "TOTAL", "CASE", "REPORT", "ONLINE", "PRINT",
            "SIGNIFICANT", "STATISTICAL", "VALUE", "MEAN", "RATE", "RATIO", "STABLE", "UNSTABLE",
            
            # Termos Biológicos Genéricos (Para sobrar apenas os alvos específicos)
            "RECEPTOR", "CHANNEL", "PROTEIN", "GENE", "FACTOR", "RESPONSE", "CELLS", "TISSUE", "MODEL",
            "USING", "MECHANISM", "SIGNALING", "NOVEL", "NEW", "ACUTE", "CHRONIC", "HUMAN", "MOUSE",
            "EXPRESSION", "PATHWAY", "TARGET", "URINARY", "BLADDER", "URETHRA", "KIDNEY"
        }
        
        unidades = {"MMHG", "KPA", "MIN", "SEC", "HRS", "ML", "MG", "KG", "NM", "UM", "MM", "NMOL"}
        
        # Regex: Captura siglas, mas remove pontuações grudadas (como 'CBS,')
        encontrados = re.findall(r'\b[A-Z]{3,8}\b|\b[A-Z]{1,3}[0-9]{1,2}\b', dados)
        
        candidatos = []
        for t in encontrados:
            # Limpa qualquer resíduo de caractere não alfabético
            t_clean = re.sub(r'[^A-Z0-9]', '', t)
            
            if t_clean in blacklist or t_clean in unidades: continue
            if len(t_clean) < 3 or t_clean == termo_base.upper(): continue
            if t_clean.isdigit() or len(re.findall(r'[0-9]', t_clean)) > 2: continue
            
            # Filtro final para palavras comuns de 3 letras que a blacklist pode ter pulado
            if t_clean in {"APP", "THE", "AND", "FOR", "NOT", "BUT", "WAS", "ARE", "HAS", "HAD"}: continue
            
            candidatos.append(t_clean)
        
        return [t for t, count in Counter(candidatos).most_common(7)]
    except: return []

# --- RADAR ---
def buscar_todas_noticias(lang='pt'):
    try:
        query = "(molecular pharmacology OR drug targets) AND (2024/08/01:2025/12/31[Date - Publication])"
        handle = Entrez.esearch(db="pubmed", term=query, retmax=6, sort="pub_date")
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        
        news = []
        for art in dados.split("\n\nPMID-"):
            tit, pmid, journal = "", "", ""
            for line in art.split("\n"):
                if re.match(r'^TI\s+-', line): tit = re.sub(r'^TI\s+-\s+', '', line).strip()
                if line.startswith("JT  - "): journal = line.replace("JT  - ", "").strip()
                if line.strip().isdigit() and not pmid: pmid = line.strip()
            if tit and pmid:
                news.append({
                    "titulo": tit, "fonte": journal[:30], 
                    "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                    "img": "https://images.unsplash.com/photo-1532187863486-abf9dbad1b69?w=400"
                })
        return news
    except: return []


