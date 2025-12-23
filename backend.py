import streamlit as st
from Bio import Entrez
import google.generativeai as genai
from tenacity import retry, stop_after_attempt, wait_exponential
import re
from collections import Counter
import time

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br" 

# --- IA: CÉREBRO DIGITAL (MODO DETETIVE V2.12) ---
def analisar_abstract_com_ia(titulo, abstract, api_key, lang='pt'):
    if not api_key:
        return "⚠️ IA não ativada (Insira a Chave na Configuração)"
    
    try:
        genai.configure(api_key=api_key)
        idioma_resp = "Português" if lang == 'pt' else "Inglês"
        
        # Modo Sherlock: Se o resumo for curto/ausente, foca no título
        abstract_limpo = abstract if (abstract and len(abstract) > 20) else "Abstract not available. Infer from title."
        
        # PROMPT ULTRA-FOCADO V2.12 (Consome menos cota e evita cortes)
        prompt = f"""
        Objective: Summary for a PhD researcher. 
        Respond in {idioma_resp}.
        
        DATA: {titulo} | {abstract_limpo[:4000]}
        
        TASK:
        Summarize in exactly 3 short bullet points:
        - ALVO PRINCIPAL (Ex: Canal Piezo1)
        - MECANISMO (Ação na bexiga/detrusor)
        - RELEVÂNCIA (Impacto clínico/funcional)
        
        Be extremely concise. Use technical terms.
        """

        # Modelos com nomes oficiais para evitar erro 404
        modelos_para_testar = ['gemini-1.5-flash', 'gemini-1.5-pro']
        
        for nome_modelo in modelos_para_testar:
            try:
                model = genai.GenerativeModel(nome_modelo)
                response = model.generate_content(prompt)
                return response.text.strip()
            except Exception as e:
                erro_msg = str(e)
                if "429" in erro_msg:
                    time.sleep(2) # Pausa maior se for cota
                continue 
        
        # FALLBACK: Se a IA falhar por cota, mostra o início do abstract original
        return f"📖 {abstract[:300]}... (IA ocupada, tente em 1 min)"
        
    except Exception as e:
        return f"❌ Erro Crítico: {str(e)[:50]}"

# --- FUNÇÕES DE BUSCA ---
@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
def _fetch_pubmed_count(query):
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle)
        handle.close()
        return int(record["Count"])
    except:
        return 0

@st.cache_data(ttl=86400, show_spinner=False)
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    if email and "@" in email: Entrez.email = email
    query_parts = [f"({termo})"]
    if contexto: query_parts.append(f"({contexto})") 
    query_parts.append(f"({ano_ini}:{ano_fim}[Date - Publication])")
    query_parts.append("(NOT Review[pt])") 
    query_final = " AND ".join(query_parts)
    return _fetch_pubmed_count(query_final)

# --- LEITOR DE RESUMOS ---
@st.cache_data(ttl=86400, show_spinner=False)
def buscar_resumos_detalhados(termo, orgao, email, ano_ini, ano_fim):
    if email and "@" in email: Entrez.email = email
    query = f"({termo}) AND ({orgao}) AND ({ano_ini}:{ano_fim}[Date - Publication]) AND (NOT Review[pt])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        if not id_list: return []
        
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        artigos_txt = handle.read()
        handle.close()
        
        artigos_finais = []
        raw_list = artigos_txt.split("\n\nPMID-")
        
        for raw in raw_list[:5]:
            tit, abs_txt, link_id, last_tag = "", "", "", ""
            lines = raw.split("\n")
            
            for line in lines:
                if line.strip().isdigit() and not link_id: link_id = line.strip()
                if line.startswith("TI  - "): 
                    tit = line.replace("TI  - ", "").strip()
                    last_tag = "TI"
                elif line.startswith("      ") and last_tag == "TI":
                    tit += " " + line.strip()
                elif line.startswith("AB  - "): 
                    abs_txt = line.replace("AB  - ", "").strip()
                    last_tag = "AB"
                elif line.startswith("      ") and last_tag == "AB":
                    abs_txt += " " + line.strip()
                elif len(line) > 4 and line[4] == "-":
                    last_tag = ""
            
            if tit:
                artigos_finais.append({
                    "Title": tit, "Resumo_Original": abs_txt, 
                    "Link": f"https://pubmed.ncbi.nlm.nih.gov/{link_id}/"
                })
        return artigos_finais
    except: return []

# --- MINERAÇÃO DE ALVOS (SMART MINER) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email):
    if email and "@" in email: Entrez.email = email
    # Busca focada em mecanismos e receptores recentes
    query = f"({termo_base}) AND (receptor OR pathway OR channel OR signaling) AND (2023:2030[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=50)
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]: return []
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read()
        handle.close()
        
        candidatos = re.findall(r'\b[A-Z][A-Z0-9]{2,8}\b', dados)
        blacklist = {"THE", "AND", "FOR", "PMID", "PMC", "DOI", "USA", "TYPE", "CELL", "ROLE", "DATA"}
        validos = [c for c in candidatos if c not in blacklist and c not in termo_base.upper()]
        
        return [item for item, count in Counter(validos).most_common(10)]
    except: return []

# --- RADAR CIENTÍFICO DINÂMICO ---
@st.cache_data(ttl=3600)
def buscar_todas_noticias(lang='pt'):
    # Radar focado no seu Doutorado: Bexiga, ROS e Farmacologia
    termos_radar = "((bladder OR detrusor) AND (oxidative stress OR reactive oxygen species OR pharmacology))"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=termos_radar, retmax=10, sort="pub_date")
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        if not id_list: return []
        
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        dados = handle.read()
        handle.close()
        
        artigos_raw = dados.split("\n\nPMID-")
        news = []
        
        for artigo in artigos_raw:
            tit = ""
            pmid = ""
            for line in artigo.split("\n"):
                if line.startswith("TI  - "): tit = line.replace("TI  - ", "").strip()
                if line.strip().isdigit() and not pmid: pmid = line.strip()
            
            if tit and len(news) < 6:
                news.append({
                    "titulo": tit,
                    "fonte": "PubMed Recent",
                    "img": "https://images.unsplash.com/photo-1532187863486-abf9dbad1b69?w=400",
                    "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                    "bandeira": "🧬"
                })
        return news
    except:
        return [{"titulo": "Aguardando novas publicações...", "fonte": "PubMed", "img": "", "link": "#", "bandeira": "⏳"}]
