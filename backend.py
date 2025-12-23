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
    
    try:
        genai.configure(api_key=api_key)
        idioma = "Português" if lang == 'pt' else "Inglês"
        
        # Garante que temos conteúdo para enviar
        abstract_input = abstract[:3000] if (abstract and len(abstract) > 30) else "Abstract incompleto/ausente."
        
        # PROMPT REFINADO: Individualidade, Contexto e Farmacologia Sênior
        prompt = f"""
        Você é um Pesquisador Sênior em Farmacologia e Fisiologia.
        
        DADOS DO ARTIGO:
        TÍTULO: {titulo}
        RESUMO: {abstract_input}
        
        TAREFA:
        1. Identifique o Alvo Molecular e o Fármaco/Substância associada.
        2. Descreva o efeito funcional.
        3. OBRIGATÓRIO: Mencione o CONTEXTO específico do título (ex: se fala de 'exercício', 'cardiovascular', 'diabetes' ou 'estresse oxidativo').
        
        FORMATO OBRIGATÓRIO:
        Alvo → Fármaco → Efeito (Contextualizado ao Título).
        
        REGRAS:
        - Máximo 25 palavras. 
        - Respostas individuais e únicas para cada paper.
        - Idioma: {idioma}.
        - Nunca responda 'dados insuficientes'.
        """

        modelos_disponiveis = [
            'gemini-1.5-flash',
            'gemini-1.5-flash-8b', 
            'gemini-1.5-pro',
            'gemini-2.0-flash-exp'
        ]
        
        for nome_modelo in modelos_disponiveis:
            try:
                model = genai.GenerativeModel(nome_modelo)
                # Temperatura 0.7 para evitar respostas idênticas/clichês
                response = model.generate_content(prompt, generation_config={"temperature": 0.7})
                return response.text.strip()
            except Exception as e:
                if "429" in str(e):
                    time.sleep(1.5)
                continue 
        
        # Fallback inteligente baseado nos seus temas de estudo
        t_up = titulo.upper()
        if "PIEZO" in t_up:
            ctx = "cardiovascular" if "CARDIO" in t_up else "vesical"
            return f"Piezo1 → Yoda1/GsMTx4 → Modulação da mecanotransdução no sistema {ctx}."
        if "ROS" in t_up or "OXIDATIVE" in t_up:
            return f"ROS/NOX → Antioxidantes → Controle do estresse oxidativo no contexto de: {titulo[:30]}..."
            
        return f"💡 Inferência: {titulo[:50]}... → (Modelos ocupados, tente em 1 min)"
        
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
                # Captura flexível de tags para evitar falhas de espaços
                if re.match(r'^TI\s+-', line): 
                    tit = re.sub(r'^TI\s+-\s+', '', line).strip()
                    last_tag = "TI"
                elif line.startswith("      ") and last_tag == "TI": 
                    tit += " " + line.strip()
                elif re.match(r'^AB\s+-', line): 
                    abs_txt = re.sub(r'^AB\s+-\s+', '', line).strip()
                    last_tag = "AB"
                elif line.startswith("      ") and last_tag == "AB": 
                    abs_txt += " " + line.strip()
                elif len(line) > 4 and line[4] == "-": 
                    last_tag = ""
            if tit:
                artigos.append({"Title": tit, "Resumo_Original": abs_txt, "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
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
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read().upper(); handle.close()
        
        candidatos_validos = []
        blacklist = {"THE", "AND", "FOR", "NOT", "BUT", "WITH", "FROM", "AFTER", "BEFORE", "STUDY", "DATA", "CELLS", "MODEL", "ROLE", "URINARY", "BLADDER", "PMID", "PMC", "DOI", "TYPE", "REVIEW"}
        whitelist = {"STING", "PIEZO", "YODA", "ROCK", "ASIC", "TRP", "GPR", "ROS", "NO", "ATP", "MTOR", "NFKB", "P2X", "P2Y"}

        encontrados = re.findall(r'\b[A-Z][A-Z0-9-]{2,8}\b', dados)
        for termo in encontrados:
            if termo in whitelist or (termo not in blacklist and len(termo) >= 3 and termo not in termo_base.upper()):
                candidatos_validos.append(termo)
        
        return [t for t, count in Counter(candidatos_validos).most_common(7)]
    except: return []

@st.cache_data(ttl=3600)
def buscar_todas_noticias(lang='pt'):
    try:
        handle = Entrez.esearch(db="pubmed", term="(bladder) AND (pharmacology OR oxidative stress)", retmax=3, sort="pub_date")
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        news = []
        for art in dados.split("\n\nPMID-"):
            tit, pmid = "", ""
            for line in art.split("\n"):
                if re.match(r'^TI\s+-', line): tit = re.sub(r'^TI\s+-\s+', '', line).strip()
                if line.strip().isdigit() and not pmid: pmid = line.strip()
            if tit: news.append({"titulo": tit, "fonte": "PubMed", "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return news
    except: return [{"titulo": "Aguardando novas publicações...", "fonte": "PubMed", "link": "#"}]
