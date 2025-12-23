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
        
        prompt = f"""Como PhD em Farmacologia, analise o seguinte paper de forma individualizada:
        TÍTULO: {titulo} | RESUMO: {abs_input}
        TAREFA: Identifique o Alvo Molecular e o Fármaco/Substância. Descreva o efeito funcional no sistema biológico citado.
        FORMATO: Alvo → Fármaco → Efeito (Contextualizado ao Título).
        REGRAS: Máx 25 palavras. Resposta técnica e única. Idioma: {idioma}."""

        for mod in ['gemini-1.5-flash-8b', 'gemini-1.5-flash', 'gemini-1.5-pro']:
            try:
                model = genai.GenerativeModel(mod)
                response = model.generate_content(prompt, generation_config={"temperature": 0.7})
                return response.text.strip()
            except:
                time.sleep(1.5); continue 
        return f"💡 IA Ocupada: {titulo[:40]}..."
    except Exception as e:
        return f"❌ Erro: {str(e)[:30]}"

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
    query = f"({termo}[Title/Abstract]) AND ({orgao}[Title/Abstract]) AND ({ano_ini}:{ano_fim}[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []
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

# --- MOTOR DE MINERAÇÃO UNIVERSAL ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email):
    if email: Entrez.email = email
    query = f"({termo_base}[Title/Abstract]) AND (2024:2030[Date - Publication]) AND (NOT Review[pt])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=50, sort="relevance")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        artigos_raw = full_data.split("\n\nPMID-")

        # FILTRO DE CONTEXTO UNIVERSAL (Mecanismos moleculares transversais)
        keywords_funcionais = {"SIGNALING", "PATHWAY", "RECEPTOR", "CHANNEL", "ACTIVATION", "INHIBITION", "EXPRESSION", "REGULATION", "MEDIATED", "MECHANISM", "FUNCTION", "ROLE", "TARGET", "MOLECULAR", "GENE", "PROTEIN", "ENZYME", "KINASE", "AUTOPHAGY", "APOPTOSIS", "DIFFERENTIATION", "HOMEOSTASIS", "STRESS"}

        # BLACKLIST v6.0 (Limpeza total de ruído linguístico e administrativo)
        blacklist = {"CBS", "FAU", "AID", "BRISTOL", "COMPANY", "INC", "CORP", "LTD", "PHST", "AUID", "ORCID", "CLINICAL", "RESEARCH", "MEDICAL", "MEDICINE", "HOSPITAL", "UNIVERSITY", "INSTITUTE", "LABORATORY", "CENTER", "DEPT", "SCHOOL", "FOUNDATION", "SCIENCES", "MEDLINE", "GERMANY", "MERCK", "DEC", "LID", "CELL", "MYERS", "FACULTY", "REPORTS", "PATIENTS", "ARTICLE", "FEES", "PERSONAL", "USA", "PRC", "UK", "EU", "CHINA", "NANJING", "FREIBURG", "UNIFESP", "BRAZIL", "SHANGHAI", "SCIENCE", "JOURNAL", "PFIZER", "PUBMED", "SQUIBB", "HEALTH", "WORK", "WERE", "THAT", "THIS", "THESE", "THOSE", "WHICH", "WHEN", "WHERE", "ALSO", "THAN", "BOTH", "UPON", "ONLY", "BEEN", "SOME", "COULD", "WELL", "VERY", "FROM", "INTO", "WITH", "NOT", "BUT", "WAS", "ARE", "HAS", "HAD", "ALL", "BEING", "FOR", "AND", "THE", "ABOUT", "STUDY", "DATA", "ANALYSIS", "RESULTS", "METHODS", "CONCLUSION", "AIMS", "SUMMARY", "PMID", "PMC", "DOI", "ISSN", "URL", "HTTP", "WWW", "PUBLISHED", "COPYRIGHT", "LICENSE", "FIGURE", "TABLE", "TYPE", "CLASS", "NORMAL", "TOTAL", "CASE", "REPORT", "ONLINE", "PRINT", "SIGNIFICANT", "STATISTICAL", "VALUE", "MEAN", "RATE", "RATIO", "STABLE", "UNSTABLE", "EXPRESSION", "PATHWAY"}
        unidades = {"MMHG", "KPA", "MIN", "SEC", "HRS", "ML", "MG", "KG", "NM", "UM", "MM", "NMOL"}

        candidatos_por_artigo = []
        for artigo in artigos_raw:
            texto_upper = artigo.upper()
            if not any(kw in texto_upper for kw in keywords_funcionais): continue 

            encontrados = re.findall(r'\b[A-Z][A-Z0-9-]{2,8}\b', texto_upper)
            candidatos_locais = set()
            for t in encontrados:
                t_clean = re.sub(r'[^A-Z0-9]', '', t)
                if t_clean in blacklist or t_clean in unidades or len(t_clean) < 3 or t_clean == termo_base.upper(): continue
                if t_clean.isdigit() or len(re.findall(r'[0-9]', t_clean)) > 2: continue
                candidatos_locais.add(t_clean)
            candidatos_por_artigo.extend(list(candidatos_locais))

        if not candidatos_por_artigo: return []
        contagem = Counter(candidatos_por_artigo)
        total_docs = max(1, len(artigos_raw))
        
        # Filtro de Frequência Relativa (0.35) para manter apenas alvos específicos do tema
        return [termo for termo, freq in contagem.most_common(12) if (freq/total_docs) < 0.35][:7]
    except: return []

# --- RADAR ---
def buscar_todas_noticias(lang='pt'):
    try:
        # Foco em Biologia Molecular e Farmacologia Geral 2025
        query = "(molecular biology OR pharmacology OR drug discovery) AND (2024/09/01:2025/12/31[Date - Publication])"
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
                news.append({"titulo": tit, "fonte": journal[:30], "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/", "img": "https://images.unsplash.com/photo-1532187863486-abf9dbad1b69?w=400"})
        return news
    except: return []
