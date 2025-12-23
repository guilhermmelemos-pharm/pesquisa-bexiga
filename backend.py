import streamlit as st
from Bio import Entrez
import google.generativeai as genai
from tenacity import retry, stop_after_attempt, wait_exponential
import re
from collections import Counter

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br" 

# --- IA: CÉREBRO DIGITAL (MODO DETETIVE) ---
def analisar_abstract_com_ia(titulo, abstract, api_key, lang='pt'):
    if not api_key:
        return "⚠️ IA não ativada (Insira a Chave na Configuração)"
    
    try:
        genai.configure(api_key=api_key)
        idioma_resp = "Português" if lang == 'pt' else "Inglês"
        abstract_seguro = abstract[:8000] if abstract else "" # Aumentado
        
        # --- PROMPT V2.10: SHERLOCK HOLMES ---
        prompt = f"""
        Atue como um Pesquisador Sênior Investigativo.
        
        Você recebeu os seguintes dados de um paper (pode ser apenas o Título ou um Resumo parcial):
        
        TÍTULO: {titulo}
        TEXTO DISPONÍVEL: {abstract_seguro}
        
        SUA MISSÃO:
        1. Se o resumo estiver completo: Faça a análise farmacológica padrão (Alvo -> Fármaco -> Efeito).
        2. Se o resumo estiver CORTADO ou AUSENTE: NÃO DESISTA. Use o TÍTULO para inferir o tema.
           - Exemplo: Se o título é "Effect of Mirabegron on Bladder", e o texto sumiu, diga: "O estudo foca no Mirabegron (agonista Beta-3), sugerindo investigação sobre relaxamento do detrusor, baseando-se no título."
        
        REGRAS:
        - Nunca responda apenas "Dados insuficientes" a menos que o título também seja vazio.
        - Seja direto e técnico.
        - Idioma: {idioma_resp}.
        """

        modelos_para_testar = [
            'models/gemini-2.5-flash',
            'models/gemini-2.0-flash',
            'models/gemini-1.5-flash',
            'models/gemini-1.5-pro',
            'models/gemini-pro'
        ]
        
        erros_coletados = []
        for nome_modelo in modelos_para_testar:
            try:
                model = genai.GenerativeModel(nome_modelo)
                response = model.generate_content(prompt)
                return response.text.strip()
            except Exception as e:
                erros_coletados.append(f"{nome_modelo}: {str(e)}")
                continue 
        
# Isso vai imprimir o erro técnico na tela (Ex: 400 Bad Request, 403 Forbidden)
    return f"❌ ERRO TÉCNICO: {str(erros_coletados)}"
        
    except Exception as e:
        return f"❌ Erro Crítico: {str(e)[:50]}..."

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

# --- LEITOR CORRIGIDO (BUG FIX CRÍTICO) ---
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
            tit = ""
            abs_txt = ""
            link_id = ""
            lines = raw.split("\n")
            
            # LÓGICA DE LEITURA REFINADA (Lê linhas quebradas)
            last_tag = ""
            for line in lines:
                # Captura ID
                if line.strip().isdigit(): 
                    link_id = line.strip()
                
                # Captura Título
                if line.startswith("TI  - "): 
                    tit = line.replace("TI  - ", "").strip()
                    last_tag = "TI"
                elif line.startswith("      ") and last_tag == "TI": # Continuação do título
                    tit += " " + line.strip()

                # Captura Abstract (CORREÇÃO AQUI)
                elif line.startswith("AB  - "): 
                    abs_txt = line.replace("AB  - ", "").strip()
                    last_tag = "AB"
                elif line.startswith("      ") and last_tag == "AB": # Continuação do abstract
                    abs_txt += " " + line.strip()
                
                # Se começar outra tag (ex: AU -), para de concatenar
                elif len(line) > 4 and line[4] == "-":
                    last_tag = ""
            
            if not link_id and id_list: link_id = id_list[0]
            
            if tit:
                artigos_finais.append({
                    "Title": tit, 
                    "Resumo_Original": abs_txt, 
                    "Resumo_IA": abs_txt[:300] + "...", 
                    "Link": f"https://pubmed.ncbi.nlm.nih.gov/{link_id}/"
                })
        return artigos_finais
    except: return []

# --- MINERAÇÃO INTELIGENTE (MANTIDA IGUAL) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email):
    if email and "@" in email: Entrez.email = email
    
    query = f"({termo_base}) AND (receptor OR pathway OR channel OR expression) AND (2023:2030[Date - Publication])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=60, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        if not id_list: return []
        
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        dados = handle.read()
        handle.close()
        
        artigos_raw = dados.split("\n\nPMID-")
        candidatos_validos = []

        keywords_funcionais = set([
            "RECEPTOR", "PATHWAY", "CHANNEL", "EXPRESSION", "ACTIVATION", "INHIBITION", 
            "MEDIATED", "SIGNALING", "PROTEIN", "GENE", "KINASE", "ENZYME", "BINDING",
            "REGULATION", "MODULATION", "TARGET", "AGONIST", "ANTAGONIST"
        ])

        blacklist_vocabulario = set([
            "THE", "AND", "FOR", "NOT", "BUT", "WITH", "FROM", "AFTER", "BEFORE", "DURING", "BETWEEN",
            "NOVEL", "NEW", "ACUTE", "CHRONIC", "HUMAN", "MOUSE", "RAT", "CELLS", "TISSUE", "MODEL",
            "STUDY", "GROUP", "DATA", "ANALYSIS", "RESULTS", "METHODS", "CONCLUSION", "AIMS", "SUMMARY",
            "HIGH", "LOW", "INCREASED", "DECREASED", "LEVELS", "EXPRESSION", "PATHWAY", "TARGET",
            "ROLE", "EFFECT", "IMPACT", "ACTION", "SYSTEM", "FUNCTION", "ACTIVITY", "POTENTIAL",
            "TREATED", "INDUCED", "MEDIATED", "ASSOCIATED", "OBSERVED", "COMPARED", "PERFORMED",
            "URINARY", "BLADDER", "URETHRA", "KIDNEY", "LIVER", "HEART", "BRAIN", "LUNG", "MUSCLE",
            "CANCER", "TUMOR", "DISEASE", "SYNDROME", "PAIN", "INFLAMMATION", "INFECTION", "INJURY",
            "PATIENT", "CLINICAL", "TRIAL", "THERAPY", "TREATMENT", "DRUG", "DOSE", "CONTROL",
            "MALE", "FEMALE", "ADULT", "CHILD", "AGE", "YEAR", "MONTH", "DAY", "TIME",
            "UNITED", "STATES", "CHINA", "JAPAN", "BRAZIL", "EUROPE", "ASIAN", "AMERICAN",
            "FIGURE", "TABLE", "PMID", "PMC", "DOI", "ISSN", "URL", "HTTP", "WWW",
            "RUPTURE", "COMPLETE", "PARTIAL", "SINGLE", "DOUBLE", "TRIPLE", "MULTIPLE",
            "SIGNIFICANT", "STATISTICAL", "DIFFERENCE", "VALUE", "MEAN", "SCORE", "RATE", "RATIO",
            "NORMAL", "ABNORMAL", "POSITIVE", "NEGATIVE", "ACTIVE", "INACTIVE", "STABLE", "UNSTABLE",
            "TYPE", "CLASS", "GROUP", "SUBGROUP", "CATEGORY", "VERSION", "EDITION", "VOLUME", "ISSUE",
            "USING", "REVIEW", "MECHANISM", "SIGNALING", "PROTEIN", "GENE", "FACTOR", "RESPONSE"
        ])

        whitelist_biologica = set([
            "STING", "PIEZO", "YODA", "TMAO", "ROCK", "ASIC", "TRP", "GPR", "TAS", 
            "BDNF", "NGF", "VEGF", "EGFR", "TNF", "IFN", "TGF", "IL", "ROS", "NO", "CO", "H2S",
            "ATP", "ADP", "AMP", "CAMP", "CGMP", "DNA", "RNA", "MIRNA", "SIRNA", "CRISPR",
            "MTOR", "AMPK", "MAPK", "ERK", "JNK", "NFKB", "STAT", "JAK", "PI3K", "AKT",
            "P53", "P2X", "P2Y"
        ])

        regex_num = r'\b[A-Z][a-zA-Z]*[0-9]+[a-zA-Z0-9]*\b'      
        regex_hifen = r'\b[A-Z][a-zA-Z0-9]{0,4}-[A-Z0-9]{1,4}\b' 
        regex_pure = r'\b[A-Z]{3,8}\b'                           

        for artigo in artigos_raw:
            texto_limpo = ""
            lines = artigo.split("\n")
            for line in lines:
                if line.startswith("TI  - "): tit = line.replace("TI  - ", "").strip()
                if line.startswith("AB  - "): abs_txt += line.replace("AB  - ", "").strip() + " "
                if line.strip().isdigit(): link_id = line.strip()
            
            if not texto_limpo: continue
            
            texto_upper = texto_limpo.upper()
            if not any(kw in texto_upper for kw in keywords_funcionais):
                continue

            encontrados = []
            encontrados.extend(re.findall(regex_num, texto_limpo))
            encontrados.extend(re.findall(regex_hifen, texto_limpo))
            encontrados.extend(re.findall(regex_pure, texto_limpo))

            for termo in encontrados:
                t_upper = termo.upper()
                
                if t_upper in whitelist_biologica:
                    candidatos_validos.append(termo)
                    continue
                if t_upper in blacklist_vocabulario:
                    continue
                if t_upper in termo_base.upper() or len(t_upper) < 3:
                    continue
                
                candidatos_validos.append(termo)
        
        contagem = Counter(candidatos_validos)
        return [termo for termo, qtd in contagem.most_common(12)][:7]

    except:
        return []

@st.cache_data(ttl=3600)
def buscar_todas_noticias(lang='pt'):
    # Notícias estáticas (mantidas)
    return [
        {"titulo": "New bladder targets identified in 2024 review", "fonte": "Nature Urology", "img": "https://images.unsplash.com/photo-1576086213369-97a306d36557?w=400", "link": "#", "bandeira": "🔬"},
        {"titulo": "H2S donors show promise in detrusor relaxation", "fonte": "ScienceDirect", "img": "https://images.unsplash.com/photo-1532187863486-abf9dbad1b69?w=400", "link": "#", "bandeira": "💊"},
        {"titulo": "Piezo1 channels: The future of mechanotransduction", "fonte": "Cell", "img": "https://images.unsplash.com/photo-1530026405186-ed1f139313f8?w=400", "link": "#", "bandeira": "⚡"}

    ]
