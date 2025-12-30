import streamlit as st
from Bio import Entrez
import google.generativeai as genai
from tenacity import retry, stop_after_attempt, wait_exponential
import re
from collections import Counter

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br" 

# --- IA: ANÁLISE INDIVIDUAL (FORMATO NOVO: ALVO - FÁRMACO - AÇÃO) ---
def analisar_abstract_com_ia(titulo, abstract, api_key, lang='pt'):
    if not api_key: return "⚠️ IA OFF"
    
    try:
        genai.configure(api_key=api_key)
        # Otimização: Corta abstracts gigantes
        abstract_seguro = abstract[:4000] if abstract else "" 
        
        prompt = f"""
        Atue como Farmacologista PhD focado em Drug Discovery.
        
        DADOS DO ARTIGO:
        TÍTULO: {titulo}
        TEXTO: {abstract_seguro}
        
        SUA MISSÃO: Identificar o trio farmacológico.
        
        FORMATO DE SAÍDA OBRIGATÓRIO (Apenas uma linha):
        Alvo Molecular - Fármaco/Composto - Ação
        
        REGRAS:
        1. Alvo: Receptor, Enzima ou Canal (ex: P2X3, mTOR, COX-2).
        2. Fármaco: O nome da droga, composto ou "Composto Experimental" (ex: Mirabegron, Trealose, AF-353).
        3. Ação: O que a droga faz no alvo (Agonista, Inibidor, Bloqueador).
        4. Se não houver droga específica, coloque "N/A".
        5. Responda em Português.
        
        Exemplo: "Receptor Beta-3 - Mirabegron - Agonista"
        """

        modelos = ['models/gemini-1.5-flash', 'models/gemini-2.0-flash-exp']
        for m in modelos:
            try:
                model = genai.GenerativeModel(m)
                resp = model.generate_content(prompt)
                return resp.text.strip().replace("Output:", "").replace("*", "")
            except: continue 
        return "Erro IA"
    except Exception as e: return f"Erro: {str(e)[:20]}"

# --- IA: CURADORIA MASSIVA DE LISTA ---
def filtrar_candidatos_com_ia(lista_suja, api_key):
    if not api_key or not lista_suja: return lista_suja[:15]

    try:
        genai.configure(api_key=api_key)
        
        # Prompt ajustado para priorizar DRUGS e TARGETS VALIDÁVEIS
        prompt_curadoria = f"""
        Sou um caçador de fármacos. Tenho uma lista suja de termos extraídos de 500 títulos do PubMed.
        
        LISTA SUJA: {', '.join(lista_suja)}

        TAREFA: 
        Filtre e retorne APENAS:
        1. Nomes de Fármacos/Compostos (ex: Trealose, Sildenafil, Resiniferatoxin).
        2. Alvos Moleculares acionáveis (ex: P2X3, TRPV1, Rho-Kinase).
        
        REMOVA IMEDIATAMENTE:
        - Termos anatômicos (Bexiga, Mucosa, Detrusor).
        - Termos genéricos (Estudo, Análise, Efeito, Ratos, Humanos).
        - Doenças (Cistite, Inflamação).
        
        Saída: Lista limpa separada por vírgula.
        """

        model = genai.GenerativeModel('models/gemini-1.5-flash')
        response = model.generate_content(prompt_curadoria)
        
        limpos = [x.strip() for x in response.text.split(",") if x.strip()]
        return limpos if limpos else lista_suja[:10]

    except: return lista_suja[:10]

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
    # Query focada em intervenção
    q_base = f"({termo}) AND ({contexto})" if contexto else f"({termo})"
    q_farmaco = " AND (drug OR inhibitor OR agonist OR antagonist OR treatment OR compound)"
    query = f"{q_base}{q_farmaco} AND ({ano_ini}:{ano_fim}[Date - Publication])"
    try: return _fetch_pubmed_count(query)
    except: return 0

# --- LEITOR DETALHADO (Mantém leitura do abstract para a tabela final) ---
@st.cache_data(ttl=86400, show_spinner=False)
def buscar_resumos_detalhados(termo, orgao, email, ano_ini, ano_fim):
    if email and "@" in email: Entrez.email = email
    query = f"({termo}) AND ({orgao}) AND ({ano_ini}:{ano_fim}[Date - Publication]) AND (NOT Review[pt])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5, sort="relevance") # Aqui mantemos 5 para leitura profunda
        record = Entrez.read(handle); handle.close()
        id_list = record["IdList"]
        if not id_list: return []
        
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        artigos_txt = handle.read(); handle.close()
        
        artigos_finais = []
        for raw in artigos_txt.split("\n\nPMID-")[:5]:
            tit, abs_txt, link_id = "", "", ""
            lines = raw.split("\n")
            last_tag = "" 
            for line in lines:
                if line.strip().isdigit(): link_id = line.strip()
                if line.startswith("TI  - "): tit = line[6:].strip(); last_tag="TI"
                elif line.startswith("AB  - "): abs_txt = line[6:].strip(); last_tag="AB"
                elif line.startswith("      "): 
                    if last_tag=="TI": tit+=" "+line.strip()
                    elif last_tag=="AB": abs_txt+=" "+line.strip()
                elif len(line)>4 and line[4]=="-": last_tag=""
            
            if not link_id and id_list: link_id=id_list[0]
            if tit:
                artigos_finais.append({
                    "Title": tit, "Resumo_Original": abs_txt, 
                    "Resumo_IA": "...", 
                    "Link": f"https://pubmed.ncbi.nlm.nih.gov/{link_id}/"
                })
        return artigos_finais
    except: return []

# --- MINERAÇÃO TURBO (500 ARTIGOS - SÓ TÍTULO/KEYWORDS) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email):
    api_key = st.session_state.get('api_key_usuario', '')
    if email and "@" in email: Entrez.email = email
    
    # 1. QUERY DE ALTA PRESSÃO (Focada 80% em Fármacos/Ação)
    # Exige termos de ação farmacológica para o artigo aparecer
    query = (
        f"({termo_base}) AND "
        "(inhibitor OR agonist OR antagonist OR blocker OR activator OR analogue OR derivative OR treatment OR drug OR compound) "
        "AND (2020:2030[Date - Publication])"
    )
    
    try:
        # AUMENTO PARA 500 ARTIGOS
        handle = Entrez.esearch(db="pubmed", term=query, retmax=500, sort="relevance")
        record = Entrez.read(handle); handle.close()
        id_list = record["IdList"]
        if not id_list: return []
        
        # OTIMIZAÇÃO: Buscamos apenas Medline, mas vamos focar em Títulos e Keywords
        # processar 500 abstracts inteiros demoraria muito. Vamos focar nos metadados.
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        
        artigos_raw = dados.split("\n\nPMID-")
        candidatos_brutos = []
        
        # Regex ajustada para pegar nomes químicos (hífen, números) e siglas
        regex_quimico = r'\b[A-Za-z0-9-]{3,20}\b' 
        
        ignore = set(["THE", "AND", "WITH", "FOR", "NOT", "BUT", "FROM", "AFTER", "DURING", "HIGH", "LOW", 
                  "STUDY", "DATA", "GROUP", "URINARY", "BLADDER", "CELLS", "RAT", "MICE", "HUMAN", 
                  "EXPRESSION", "ACTIVITY", "ROLE", "EFFECT", "LEVELS", "TREATED", "USING", "CONTROL", 
                  "SHAM", "WEEK", "DAY", "DOSE", "MGDL", "KG", "MG", "ML", "MODEL", "INDUCED"])

        # Loop de Extração Rápida
        for art in artigos_raw:
            txt_interesse = ""
            lines = art.split("\n")
            for line in lines:
                # SÓ LÊ TÍTULO (TI) E PALAVRAS-CHAVE (OT) - Ignora Abstract (AB) para ganhar velocidade
                if line.startswith("TI  - "): txt_interesse += line[6:].strip() + " "
                elif line.startswith("OT  - "): txt_interesse += line[6:].strip() + " "
            
            palavras = re.findall(regex_quimico, txt_interesse)
            for p in palavras:
                pu = p.upper()
                if pu not in ignore and len(pu) > 2:
                    # Filtro básico: se parece gene/droga (tem numero ou hifen) ou é uma palavra muito especifica
                    if re.search(r'\d', p) or "-" in p or pu in ["MIRABEGRON", "SOLIFENACIN", "TAMSULOSIN", "RESINIFERATOXIN", "BOTOX"]:
                        candidatos_brutos.append(p)
                    # Adiciona palavras em caixa alta (siglas)
                    elif p.isupper() and len(p) > 2:
                        candidatos_brutos.append(p)

        # Pega Top 50 termos brutos para a IA limpar
        top_brutos = [t for t, q in Counter(candidatos_brutos).most_common(50)]

        if api_key:
            return filtrar_candidatos_com_ia(top_brutos, api_key)
        else:
            return top_brutos[:15]

    except Exception as e:
        print(e)
        return []

@st.cache_data(ttl=3600)
def buscar_todas_noticias(lang='pt'):
    return [
        {"titulo": "Novos inibidores de P2X3 em teste", "fonte": "Nature Urology", "img": "https://images.unsplash.com/photo-1576086213369-97a306d36557?w=400", "link": "#", "bandeira": "🧪"},
        {"titulo": "Eficácia do Vibegron vs Mirabegron", "fonte": "ScienceDirect", "img": "https://images.unsplash.com/photo-1532187863486-abf9dbad1b69?w=400", "link": "#", "bandeira": "💊"},
        {"titulo": "Toxina Botulínica: Novos mecanismos", "fonte": "Cell", "img": "https://images.unsplash.com/photo-1530026405186-ed1f139313f8?w=400", "link": "#", "bandeira": "⚡"}
    ]
