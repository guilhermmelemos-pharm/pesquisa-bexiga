import streamlit as st
from Bio import Entrez
import google.generativeai as genai
from tenacity import retry, stop_after_attempt, wait_exponential
import re
from collections import Counter
import time

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br" 

# --- IA: CÉREBRO DIGITAL (MODO SHERLOCK HOLMES + FAILOVER) ---
def analisar_abstract_com_ia(titulo, abstract, api_key, lang='pt'):
    if not api_key:
        return "⚠️ IA não ativada (Insira a Chave na Configuração)"
    
    genai.configure(api_key=api_key)
    idioma_resp = "Português" if lang == 'pt' else "Inglês"
    abstract_seguro = abstract[:3000] if abstract else "Abstract indisponível."
    
    prompt = f"""
    Atue como um Pesquisador Sênior Investigativo em Farmacologia.
    TÍTULO: {titulo}
    RESUMO: {abstract_seguro}
    
    MISSÃO:
    1. Se o resumo estiver cortado/ausente: Use o TÍTULO para inferir o mecanismo provável.
    Você é um Pesquisador Sênior em Farmacologia. 
        
        DADOS:
        TÍTULO: {titulo}
        RESUMO: {abstract_txt}
        
        TAREFA:
        1. Identifique o Alvo e o Fármaco.
        2. Descreva o efeito funcional.
        3. OBRIGATÓRIO: Sua resposta deve mencionar o CONTEXTO específico do título (ex: se o título fala de 'exercício' ou 'câncer', isso deve estar no efeito).
        
        FORMATO: Alvo → Fármaco → Efeito (Contextualizado ao Título).
        
        REGRAS: Máximo 25 palavras. Proibido repetir respostas de outros artigos. Idioma: {idioma}.
        """

    # LISTA DE MODELOS PARA TENTATIVAS (Evita 404 e 429)
    modelos_disponiveis = [
        'gemini-1.5-flash',
        'gemini-1.5-flash-8b', 
        'gemini-1.5-pro',
        'gemini-2.0-flash-exp'
    ]
    
    for nome_modelo in modelos_disponiveis:
        try:
            model = genai.GenerativeModel(nome_modelo)
            response = model.generate_content(prompt)
            return response.text.strip()
        except Exception as e:
            if "429" in str(e):
                time.sleep(1.5) # Pausa maior para limpar cota
            continue 
    
    # Fallback final baseado no conhecimento do Guilherme
    t_up = titulo.upper()
    if "PIEZO" in t_up: return "Piezo1 → Yoda1 / GsMTx4 → Mecanotransdução e sinalização de estiramento urotelial."
    if "ROS" in t_up or "OXIDATIVE" in t_up: return "ROS/NOX → Antioxidantes/SOD → Modulação do estresse oxidativo e contratilidade."
    
    return f"💡 Inferência: {titulo} → (Modelos ocupados, tente em 1 min)"

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
    query = f"({termo}) AND ({contexto}) AND ({ano_ini}:{ano_fim}[Date - Publication]) AND (NOT Review[pt])"
    try: return _fetch_pubmed_count(query)
    except: return 0

@st.cache_data(ttl=86400, show_spinner=False)
def buscar_resumos_detalhados(termo, orgao, email, ano_ini, ano_fim):
    if email and "@" in email: Entrez.email = email
    query = f"({termo}) AND ({orgao}) AND ({ano_ini}:{ano_fim}[Date - Publication]) AND (NOT Review[pt])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5, sort="relevance")
        record = Entrez.read(handle); handle.close()
        ids = record["IdList"]
        if not ids: return []
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        artigos = []
        for raw in dados.split("\n\nPMID-"):
            tit, abs_txt, pmid, last_tag = "", "", "", ""
            for line in raw.split("\n"):
                if line.strip().isdigit() and not pmid: pmid = line.strip()
                if line.startswith("TI  - "): tit = line.replace("TI  - ", "").strip(); last_tag = "TI"
                elif line.startswith("      ") and last_tag == "TI": tit += " " + line.strip()
                elif line.startswith("AB  - "): abs_txt = line.replace("AB  - ", "").strip(); last_tag = "AB"
                elif line.startswith("      ") and last_tag == "AB": abs_txt += " " + line.strip()
                elif len(line) > 4 and line[4] == "-": last_tag = ""
            if tit: artigos.append({"Title": tit, "Resumo_Original": abs_txt, "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return artigos
    except: return []

# --- MINERAÇÃO INTELIGENTE COM BLACKLIST INTEGRADA ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email):
    if email and "@" in email: Entrez.email = email
    query = f"({termo_base}) AND (receptor OR pathway OR channel) AND (2023:2030[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=50)
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        
        artigos_raw = dados.split("\n\nPMID-")
        candidatos_validos = []

        blacklist_vocabulario = {
            "THE", "AND", "FOR", "NOT", "BUT", "WITH", "FROM", "AFTER", "BEFORE", "DURING", "BETWEEN",
            "NOVEL", "NEW", "ACUTE", "CHRONIC", "HUMAN", "MOUSE", "RAT", "CELLS", "TISSUE", "MODEL",
            "STUDY", "GROUP", "DATA", "ANALYSIS", "RESULTS", "METHODS", "CONCLUSION", "AIMS", "SUMMARY",
            "ROLE", "EFFECT", "IMPACT", "ACTION", "SYSTEM", "FUNCTION", "ACTIVITY", "POTENTIAL",
            "URINARY", "BLADDER", "URETHRA", "CANCER", "TUMOR", "DISEASE", "PAIN", "MALE", "FEMALE",
            "PMID", "PMC", "DOI", "TYPE", "USING", "REVIEW", "MECHANISM", "SIGNALING", "PROTEIN"
        }

        whitelist_biologica = {"STING", "PIEZO", "YODA", "ROCK", "ASIC", "TRP", "GPR", "ROS", "NO", "ATP", "MTOR", "NFKB", "P2X", "P2Y"}

        for artigo in artigos_raw:
            # Pegamos o texto do artigo em Upper Case para facilitar a comparação
            texto_upper = artigo.upper()
            
            # Regex para siglas e alvos moleculares
            encontrados = re.findall(r'\b[A-Z][A-Z0-9-]{2,8}\b', texto_upper)

            for termo in encontrados:
                if termo in whitelist_biologica:
                    candidatos_validos.append(termo)
                elif termo not in blacklist_vocabulario and len(termo) >= 3 and termo not in termo_base.upper():
                    candidatos_validos.append(termo)
        
        return [t for t, count in Counter(candidatos_validos).most_common(7)]
    except: return []

@st.cache_data(ttl=3600)
def buscar_todas_noticias(lang='pt'):
    # Radar Dinâmico Real
    try:
        handle = Entrez.esearch(db="pubmed", term="(bladder) AND (oxidative stress OR mechanotransduction)", retmax=3, sort="pub_date")
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        news = []
        for art in dados.split("\n\nPMID-"):
            tit, pmid = "", ""
            for line in art.split("\n"):
                if line.startswith("TI  - "): tit = line.replace("TI  - ", "").strip()
                if line.strip().isdigit() and not pmid: pmid = line.strip()
            if tit: news.append({"titulo": tit, "fonte": "PubMed", "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return news
    except: return [{"titulo": "Pesquisando novidades...", "fonte": "PubMed", "link": "#"}]


