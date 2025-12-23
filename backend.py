import streamlit as st
from Bio import Entrez
import google.generativeai as genai
from tenacity import retry, stop_after_attempt, wait_exponential
import re
from collections import Counter
import time

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br" 

# --- IA: CÉREBRO DIGITAL (MODO SHERLOCK HOLMES) ---
def analisar_abstract_com_ia(titulo, abstract, api_key, lang='pt'):
    if not api_key:
        return "⚠️ IA não ativada (Insira a Chave na Configuração)"
    
    try:
        genai.configure(api_key=api_key)
        idioma_resp = "Português" if lang == 'pt' else "Inglês"
        
        # Modo Sherlock: Se o resumo for curto/ausente, foca no título
        abstract_seguro = abstract if (abstract and len(abstract) > 20) else "Abstract not available. Infer from title."
        
        prompt = f"""
        Atue como um Pesquisador Sênior Investigativo em Farmacologia.
        Analise os dados abaixo e responda em {idioma_resp}.
        
        TÍTULO: {titulo}
        RESUMO/TEXTO: {abstract_seguro[:7000]}
        
        MISSÃO:
        1. Identifique o Alvo Molecular, o Fármaco/Substância e o Efeito Fisiológico.
        2. Se o resumo estiver incompleto, use o TÍTULO para inferir o mecanismo provável.
        3. Formato da resposta: Alvo → Fármaco → Efeito.
        """

        # Modelos com nomes simplificados para evitar erro 404/429
        modelos_para_testar = ['gemini-1.5-flash', 'gemini-1.5-pro', 'gemini-pro']
        
        for nome_modelo in modelos_para_testar:
            try:
                model = genai.GenerativeModel(nome_modelo)
                response = model.generate_content(prompt)
                return response.text.strip()
            except Exception as e:
                if "429" in str(e):
                    time.sleep(2)
                continue 
        
        return f"📖 {abstract[:250]}... (IA ocupada, tente em 1 min)"
        
    except Exception as e:
        return f"❌ Erro Crítico: {str(e)[:50]}"

# --- FUNÇÕES DE BUSCA ---
@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
def _fetch_pubmed_count(query):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    record = Entrez.read(handle)
    handle.close()
    return int(record["Count"])

@st.cache_data(ttl=86400, show_spinner=False)
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    if email and "@" in email: Entrez.email = email
    query_parts = [f"({termo})"]
    if contexto: query_parts.append(f"({contexto})") 
    query_parts.append(f"({ano_ini}:{ano_fim}[Date - Publication])")
    query_parts.append("(NOT Review[pt])") 
    query_final = " AND ".join(query_parts)
    try:
        return _fetch_pubmed_count(query_final)
    except:
        return 0

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
                    "Title": tit, 
                    "Resumo_Original": abs_txt, 
                    "Link": f"https://pubmed.ncbi.nlm.nih.gov/{link_id}/"
                })
        return artigos_finais
    except: return []

# --- MINERAÇÃO INTELIGENTE COM BLACKLIST CORRIGIDA ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email):
    if email and "@" in email: Entrez.email = email
    query = f"({termo_base}) AND (receptor OR pathway OR channel) AND (2023:2030[Date - Publication])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=50)
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]: return []
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read()
        handle.close()
        
        artigos_raw = dados.split("\n\nPMID-")
        candidatos_validos = []

        blacklist_vocabulario = {"THE", "AND", "FOR", "NOT", "BUT", "WITH", "FROM", "STUDY", "DATA", "RESULTS", "CELLS", "MODEL", "ROLE", "IMPACT"}
        whitelist_biologica = {"STING", "PIEZO", "YODA", "ROCK", "ASIC", "TRP", "GPR", "ROS", "NO", "ATP", "MTOR", "NFKB", "P2X", "P2Y"}

        for artigo in artigos_raw:
            # Extração simples de título e abstract para o minerador
            texto_artigo = artigo.upper()
            
            # Regex para encontrar potenciais alvos (Siglas em caixa alta com números ou hifens)
            encontrados = re.findall(r'\b[A-Z][A-Z0-9-]{2,8}\b', texto_artigo)

            for termo in encontrados:
                if termo in whitelist_biologica:
                    candidatos_validos.append(termo)
                elif termo not in blacklist_vocabulario and len(termo) >= 3 and termo not in termo_base.upper():
                    candidatos_validos.append(termo)
        
        contagem = Counter(candidatos_validos)
        return [termo for termo, qtd in contagem.most_common(7)]

    except:
        return []

@st.cache_data(ttl=3600)
def buscar_todas_noticias(lang='pt'):
    return [
        {"titulo": "Piezo1 channels: The future of mechanotransduction", "fonte": "Cell", "img": "https://images.unsplash.com/photo-1530026405186-ed1f139313f8?w=400", "link": "#", "bandeira": "⚡"},
        {"titulo": "New bladder targets identified in 2024", "fonte": "Nature Urology", "img": "https://images.unsplash.com/photo-1576086213369-97a306d36557?w=400", "link": "#", "bandeira": "🔬"},
        {"titulo": "H2S donors show promise in detrusor relaxation", "fonte": "ScienceDirect", "img": "https://images.unsplash.com/photo-1532187863486-abf9dbad1b69?w=400", "link": "#", "bandeira": "💊"}
    ]
