import streamlit as st
from Bio import Entrez
import google.generativeai as genai
from tenacity import retry, stop_after_attempt, wait_exponential
import re
from collections import Counter

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br" 

# --- IA: ANÁLISE INDIVIDUAL ---
def analisar_abstract_com_ia(titulo, abstract, api_key, lang='pt'):
    if not api_key: return "⚠️ IA OFF"
    try:
        genai.configure(api_key=api_key)
        # Corta para economizar tokens, mas mantém o suficiente para contexto
        abstract_seguro = abstract[:3500] if abstract else "" 
        
        prompt = f"""
        Papel: Farmacologista Drug Discovery.
        Dado o abstract abaixo, identifique o trio principal.
        
        TÍTULO: {titulo}
        TEXTO: {abstract_seguro}
        
        SAÍDA (Uma linha): Alvo - Fármaco - Ação
        
        Regras:
        - Alvo: Receptor/Enzima (ex: P2X3, mTOR).
        - Fármaco: Nome da droga ou composto (ex: Gefapixant, Trealose). Se não houver, "N/A".
        - Ação: Inibidor, Agonista, Antagonista, Ativador.
        - Idioma: Português.
        """
        
        # Tenta modelos do mais barato/rápido para o mais parrudo
        modelos = ['models/gemini-1.5-flash', 'models/gemini-2.0-flash-exp']
        for m in modelos:
            try:
                model = genai.GenerativeModel(m)
                resp = model.generate_content(prompt)
                clean = resp.text.strip().replace("Output:", "").replace("*", "")
                return clean if len(clean) > 5 else "N/A - N/A - N/A"
            except: continue
        return "Erro IA"
    except: return "Erro Conexão"

# --- IA: FILTRO DE LISTA ---
def filtrar_candidatos_com_ia(lista_suja, api_key):
    if not api_key or not lista_suja: return lista_suja[:10]
    try:
        genai.configure(api_key=api_key)
        prompt = f"""
        Sou pesquisador. Tenho uma lista suja de termos extraídos do PubMed.
        LISTA SUJA: {', '.join(lista_suja)}
        
        TAREFA: Retorne APENAS os Fármacos, Compostos Químicos e Alvos Moleculares Reais.
        
        REGRAS DE EXCLUSÃO (CRÍTICO):
        - REMOVA: Termos anatômicos (Bexiga, Mucosa), genéricos (Estudo, Efeito, Ratos), doenças (Dor, Cistite).
        
        SAÍDA: Lista limpa separada por vírgula.
        """
        model = genai.GenerativeModel('models/gemini-1.5-flash')
        resp = model.generate_content(prompt)
        limpos = [x.strip() for x in resp.text.split(",") if x.strip()]
        return limpos if limpos else lista_suja[:10]
    except: return lista_suja[:10]

# --- BUSCA AUXILIAR ---
@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
def _fetch_pubmed_count(query):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    record = Entrez.read(handle); handle.close()
    return int(record["Count"])

@st.cache_data(ttl=86400, show_spinner=False)
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    if email and "@" in email: Entrez.email = email
    q_base = f"({termo}) AND ({contexto})" if contexto else f"({termo})"
    # Garante que estamos contando artigos farmacológicos/mecanísticos
    q_farmaco = " AND (drug OR treatment OR mechanism OR receptor OR inhibition)"
    query = f"{q_base}{q_farmaco} AND ({ano_ini}:{ano_fim}[Date - Publication])"
    try: return _fetch_pubmed_count(query)
    except: return 0

# --- LEITOR DETALHADO (Mantido igual) ---
@st.cache_data(ttl=86400, show_spinner=False)
def buscar_resumos_detalhados(termo, orgao, email, ano_ini, ano_fim):
    if email and "@" in email: Entrez.email = email
    query = f"({termo}) AND ({orgao}) AND ({ano_ini}:{ano_fim}[Date - Publication]) AND (NOT Review[pt])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5, sort="relevance")
        record = Entrez.read(handle); handle.close()
        id_list = record["IdList"]
        if not id_list: return []
        
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        txt = handle.read(); handle.close()
        
        finais = []
        for raw in txt.split("\n\nPMID-")[:5]:
            tit, abs_t, lid = "", "", ""
            lines = raw.split("\n")
            lt = ""
            for l in lines:
                if l.strip().isdigit(): lid = l.strip()
                if l.startswith("TI  - "): tit=l[6:].strip(); lt="TI"
                elif l.startswith("AB  - "): abs_t=l[6:].strip(); lt="AB"
                elif l.startswith("      "):
                    if lt=="TI": tit+=" "+l.strip()
                    elif lt=="AB": abs_t+=" "+l.strip()
                elif len(l)>4 and l[4]=="-": lt=""
            
            if not lid and id_list: lid=id_list[0]
            if tit:
                finais.append({"Title": tit, "Resumo_Original": abs_t, "Resumo_IA": "...", "Link": f"https://pubmed.ncbi.nlm.nih.gov/{lid}/"})
        return finais
    except: return []

# --- MINERAÇÃO CORRIGIDA (V2.4) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email):
    api_key = st.session_state.get('api_key_usuario', '')
    if email and "@" in email: Entrez.email = email
    
    # 1. Query: Força artigos que tenham "Drug/Inhibitor" E "Receptor/Pathway"
    query = (
        f"({termo_base}) AND (receptor OR channel OR enzyme) "
        f"AND (inhibitor OR agonist OR antagonist OR drug OR treatment) "
        f"AND (2022:2030[Date - Publication])"
    )
    
    try:
        # Baixamos 100 artigos (equilíbrio entre velocidade e variedade)
        handle = Entrez.esearch(db="pubmed", term=query, retmax=100, sort="relevance")
        record = Entrez.read(handle); handle.close()
        id_list = record["IdList"]
        if not id_list: return []
        
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        
        artigos_raw = dados.split("\n\nPMID-")
        candidatos_brutos = []
        
        ignore = set(["THE", "AND", "FOR", "NOT", "BUT", "WITH", "FROM", "AFTER", "HIGH", "LOW", "STUDY", "DATA", "GROUP", "BLADDER", "CELLS", "RATS", "MICE", "HUMAN", "EXPRESSION", "ACTIVITY", "EFFECT", "LEVELS", "TREATED", "USING", "CONTROL", "DOSE", "MODEL", "INDUCED", "PATHWAY", "RECEPTOR", "MECHANISM", "SIGNALING"])

        # Regex para palavras que parecem drogas ou alvos
        # 1. Siglas com números/hífens (P2X3, HO-1, NF-kB)
        # 2. Palavras terminadas em sufixos farmacológicos comuns (-in, -il, -on, -ol, -mab, -ax)
        # 3. Palavras em CAPS LOCK com pelo menos 3 letras (BDNF, NGF)
        
        regex_complexo = r'\b[A-Za-z0-9-]{3,20}\b'

        for art in artigos_raw:
            txt_interesse = ""
            lines = art.split("\n")
            for line in lines:
                # Agora lemos Abstract (AB) também, pois drogas aparecem lá
                if line.startswith("TI  - ") or line.startswith("AB  - ") or line.startswith("OT  - "):
                    txt_interesse += line[6:].strip() + " "
            
            palavras = re.findall(regex_complexo, txt_interesse)
            for p in palavras:
                pu = p.upper()
                if len(pu) < 3 or pu in ignore: continue
                
                # Critérios de Aceite (Heurística)
                is_gene = bool(re.search(r'\d', p) or "-" in p) # Tem numero ou hifen? (P2X3, IL-6)
                is_drug_suffix = any(pu.endswith(s) for s in ["IN", "IL", "ON", "OL", "MAB", "ATE", "IDE"]) # Sufixos (Mirabegron, Sildenafil)
                is_caps = p.isupper() # Siglas (BDNF, NGF)
                
                if is_gene or is_drug_suffix or is_caps:
                    candidatos_brutos.append(p)

        # Manda o Top 40 para a IA filtrar
        top_brutos = [t for t, q in Counter(candidatos_brutos).most_common(40)]

        if api_key:
            return filtrar_candidatos_com_ia(top_brutos, api_key)
        else:
            return top_brutos[:15]

    except Exception as e:
        print(e)
        return []

# Notícias (Placeholder mantido)
@st.cache_data(ttl=3600)
def buscar_todas_noticias(lang='pt'):
    return [{"titulo": "Carregando notícias...", "fonte": "PubMed", "img": "", "link": "#", "bandeira": "📡"}]
