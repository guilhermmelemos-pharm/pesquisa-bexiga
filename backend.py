import streamlit as st
from Bio import Entrez
import google.generativeai as genai
from tenacity import retry, stop_after_attempt, wait_exponential
import re
from collections import Counter
import time

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br"

# --- IA: EXTRAÇÃO VIA DADOS CURTOS (ANTI-COOLDOWN) ---
def analisar_abstract_com_ia(titulo, dados_curtos, api_key, lang='pt'):
    if not api_key:
        return "⚠️ IA não ativada"
    
    try:
        genai.configure(api_key=api_key)
        idioma = "Português" if lang == 'pt' else "Inglês"
        
        # O prompt agora processa apenas o "filé mignon" do paper (Keywords ou Início do Abstract)
        prompt = f"""Analise como PhD em Farmacologia:
FONTE: {titulo}. {dados_curtos}
FORMATO OBRIGATÓRIO: Alvo: [Sigla] | Fármaco: [Nome] | Efeito: [Ação funcional].
REGRAS: Máximo 12 palavras. Seja técnico e direto. Idioma: {idioma}."""

        modelos = ['gemini-2.0-flash-exp', 'gemini-1.5-flash-8b', 'gemini-1.5-flash']

        for mod in modelos:
            try:
                model = genai.GenerativeModel(mod)
                response = model.generate_content(prompt, generation_config={"temperature": 0.1})
                if response and response.text:
                    return response.text.strip()
            except:
                time.sleep(1.5)
                continue

        return f"💡 IA Ocupada. Alvo provável: {titulo[:30]}..."
    
    except Exception as e:
        return f"❌ Erro: {str(e)[:30]}"

# --- BUSCA PUBMED ---
@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
def _fetch_pubmed_count(query):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    record = Entrez.read(handle)
    handle.close()
    return int(record["Count"])

@st.cache_data(ttl=86400, show_spinner=False)
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    # Busca global real ignorando o contexto se ele for vazio
    query = f"({termo})"
    if contexto: query += f" AND ({contexto})"
    query += f" AND ({ano_ini}:{ano_fim}[Date - Publication]) AND (NOT Review[pt])"
    try:
        return _fetch_pubmed_count(query)
    except:
        return 0

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
            tit, pmid, keywords, fallback_text = "", "", "", ""
            lines = raw.split("\n")
            for line in lines:
                if line.strip().isdigit() and not pmid: pmid = line.strip()
                if line.startswith("TI  - "): tit = line[6:].strip()
                # 1. Captura Keywords (OT = Author Keywords, KW = Mesh Keywords)
                if line.startswith("OT  - ") or line.startswith("KW  - "):
                    keywords += line[6:].strip() + ", "
                # 2. Captura apenas o início do abstract como backup (primeiros 500 chars)
                if line.startswith("AB  - ") and not fallback_text:
                    fallback_text = line[6:500].strip()
            
            if tit:
                # Prioriza Keywords (mais limpo); se não houver, usa o início do abstract
                info_final = keywords if len(keywords) > 5 else fallback_text
                artigos.append({
                    "Title": tit, 
                    "Info_IA": info_final if info_final else "Sem resumo disponível.", 
                    "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                })
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

        keywords_funcionais = {"SIGNALING", "PATHWAY", "RECEPTOR", "CHANNEL", "ACTIVATION", "INHIBITION", "EXPRESSION",
                               "REGULATION", "MEDIATED", "MECHANISM", "FUNCTION", "ROLE", "TARGET", "MOLECULAR", "GENE",
                               "PROTEIN", "ENZYME", "KINASE", "AUTOPHAGY", "APOPTOSIS", "DIFFERENTIATION", "HOMEOSTASIS", "STRESS"}
        
        # Blacklist v6.1
        blacklist = {"AND", "THE", "FOR", "NOT", "BUT", "WITH", "FROM", "STUDY", "RESULTS", "CELLS", "WAS", "WERE", "CBS", "FAU", "AID"} 
        unidades = {"MMHG","KPA","MIN","SEC","HRS","ML","MG","KG","NM","UM","MM","NMOL"}

        candidatos_por_artigo = []
        for artigo in artigos_raw:
            texto_upper = artigo.upper()
            if not any(kw in texto_upper for kw in keywords_funcionais): continue

            encontrados = re.findall(r'\b[A-Z][A-Z0-9-]{2,8}\b', texto_upper)
            candidatos_locais = set()
            for t in encontrados:
                t_clean = re.sub(r'[^A-Z0-9]', '', t)
                if t_clean in blacklist or t_clean in unidades or len(t_clean)<3 or t_clean==termo_base.upper(): continue
                if t_clean.isdigit() or len(re.findall(r'[0-9]', t_clean))>2: continue
                candidatos_locais.add(t_clean)
            candidatos_por_artigo.extend(list(candidatos_locais))

        if not candidatos_por_artigo: return []
        contagem = Counter(candidatos_por_artigo)
        total_docs = max(1, len(artigos_raw))
        return [termo for termo,freq in contagem.most_common(12) if (freq/total_docs)<0.35][:7]
    except: return []

# --- RADAR DE NOTÍCIAS ---
def buscar_todas_noticias(lang='pt'):
    try:
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
                news.append({"titulo": tit, "fonte": journal[:30], 
                             "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                             "img":"https://images.unsplash.com/photo-1532187863486-abf9dbad1b69?w=400"})
        return news
    except: return []
        
