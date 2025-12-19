import streamlit as st
from Bio import Entrez
import time
import re
from collections import Counter

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador@unifesp.br"

# --- CACHE ---
@st.cache_data(ttl=86400, show_spinner=False)
def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    Entrez.email = email
    query_parts = [f"({termo}[Title/Abstract])"]
    if contexto:
        query_parts.append(f"({contexto}[Title/Abstract])") 
    query_parts.append(f"({ano_ini}:{ano_fim}[Date - Publication])")
    query_final = " AND ".join(query_parts)
    
    max_retries = 3
    for attempt in range(max_retries):
        try:
            handle = Entrez.esearch(db="pubmed", term=query_final, retmax=0)
            record = Entrez.read(handle)
            handle.close()
            return int(record["Count"])
        except Exception:
            time.sleep(1)
    return 0

@st.cache_data(ttl=86400, show_spinner=False)
def buscar_resumos_detalhados(termo, orgao, email, ano_ini, ano_fim):
    Entrez.email = email
    query = f"({termo}[Title/Abstract]) AND ({orgao}[Title/Abstract]) AND ({ano_ini}:{ano_fim}[Date - Publication])"
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
            for line in lines:
                if line.startswith("TI  - "): tit = line.replace("TI  - ", "").strip()
                if line.startswith("AB  - "): abs_txt += line.replace("AB  - ", "").strip() + " "
                if line.strip().isdigit(): link_id = line.strip()
            if not link_id and id_list: link_id = id_list[0]
            if tit:
                resumo_display = abs_txt[:500] + "..." if len(abs_txt) > 5 else "Resumo indisponível."
                artigos_finais.append({"Title": tit, "Resumo_IA": resumo_display, "Link": f"https://pubmed.ncbi.nlm.nih.gov/{link_id}/"})
        return artigos_finais
    except: return []

# --- MOTOR DE MINERAÇÃO: SEMANTIC WALL + FUNCTIONAL CONTEXT (v1.6.1) ---
def buscar_alvos_emergentes_pubmed(termo_base, email):
    Entrez.email = email
    # Busca 2024-2030 (Apenas novidades)
    query = f"({termo_base}[Title/Abstract]) AND (2024:2030[Date - Publication]) NOT (Review[Publication Type])"
    
    candidatos_por_artigo = []
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=60, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        if not id_list: return []

        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        full_data = handle.read()
        handle.close()
        
        artigos_raw = full_data.split("\n\nPMID-")
        
        # --- FILTRO 0: CONTEXTO FUNCIONAL (ITEM 6) ---
        # Só aceita artigos que contenham pelo menos UMA destas palavras de mecanismo.
        # Isso elimina artigos puramente administrativos, estatísticos ou clínicos genéricos.
        keywords_funcionais = set([
            "SIGNALING", "PATHWAY", "RECEPTOR", "CHANNEL", "ACTIVATION", "INHIBITION", 
            "EXPRESSION", "REGULATION", "MEDIATED", "MECHANISM", "FUNCTION", "ROLE", 
            "TARGET", "MOLECULAR", "GENE", "PROTEIN", "ENZYME", "KINASE", 
            "PHOSPHORYLATION", "APOPTOSIS", "PROLIFERATION", "MIGRATION",
            # Específicos de fisiologia (pedido no prompt)
            "RELAXATION", "CONTRACTION", "CONTRACTILITY", "MECHANOTRANSDUCTION", 
            "DETRUSOR", "UROTHELIUM", "UROTHELIAL", "SMOOTH MUSCLE", "ACTIVITY"
        ])

        # --- FILTRO 1: BLACKLIST ---
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
            "TYPE", "CLASS", "GROUP", "SUBGROUP", "CATEGORY", "VERSION", "EDITION", "VOLUME", "ISSUE"
        ])

        # --- FILTRO 2: WHITELIST ---
        whitelist_biologica = set([
            "STING", "PIEZO", "YODA", "TMAO", "ROCK", "ASIC", "TRP", "GPR", "TAS", 
            "BDNF", "NGF", "VEGF", "EGFR", "TNF", "IFN", "TGF", "IL", "ROS", "NO", "CO", "H2S",
            "ATP", "ADP", "AMP", "CAMP", "CGMP", "DNA", "RNA", "MIRNA", "SIRNA", "CRISPR",
            "MTOR", "AMPK", "MAPK", "ERK", "JNK", "NFKB", "STAT", "JAK", "PI3K", "AKT"
        ])

        # Regex Estrutural
        regex_num = r'\b[A-Z][a-zA-Z]*[0-9]+[a-zA-Z0-9]*\b'
        regex_hifen = r'\b[A-Z][a-zA-Z0-9]{0,4}-[A-Z0-9]{1,4}\b'
        regex_pure = r'\b[A-Z]{3,8}\b'

        for artigo in artigos_raw:
            texto_limpo = ""
            lines = artigo.split("\n")
            for line in lines:
                if line.startswith("TI  - ") or line.startswith("AB  - "):
                    texto_limpo += line[6:] + " "
            
            if not texto_limpo: continue
            
            # --- APLICAÇÃO DO ITEM 6 (Contexto Funcional) ---
            # Se o texto NÃO tiver nenhuma palavra funcional, pula o artigo inteiro.
            texto_upper = texto_limpo.upper()
            tem_contexto = any(kw in texto_upper for kw in keywords_funcionais)
            
            if not tem_contexto:
                continue # Pula artigos sem contexto mecanístico

            candidatos_locais = set()
            
            palavras = re.findall(regex_num, texto_limpo) + \
                       re.findall(regex_hifen, texto_limpo) + \
                       re.findall(regex_pure, texto_limpo)

            for p in palavras:
                clean_p = p.strip(".,;:()[]").upper()
                
                if clean_p in blacklist_vocabulario: continue
                if clean_p in termo_base.upper() or termo_base.upper() in clean_p: continue
                
                tem_numero = any(char.isdigit() for char in clean_p)
                tem_hifen = "-" in clean_p
                
                if not tem_numero and not tem_hifen:
                    if clean_p not in whitelist_biologica:
                        continue 
                
                if len(clean_p) < 2 or len(clean_p) > 8: continue
                if "PMC" in clean_p or "NCT" in clean_p: continue
                
                candidatos_locais.add(clean_p)
            
            if candidatos_locais:
                candidatos_por_artigo.extend(list(candidatos_locais))

        # --- FILTRO 3: FREQUÊNCIA ---
        if not candidatos_por_artigo: return []
        
        contagem = Counter(candidatos_por_artigo)
        total_docs = max(1, len(artigos_raw))
        
        lista_final = []
        for termo, freq in contagem.items():
            freq_relativa = freq / total_docs
            if freq_relativa < 0.25: 
                lista_final.append(termo)
                
        return list(set(lista_final))
        
    except Exception:
        return []

@st.cache_data(ttl=3600)
def buscar_todas_noticias(lang='pt'):
    return [
        {"titulo": "New bladder targets identified in 2024 review", "fonte": "Nature Urology", "img": "https://images.unsplash.com/photo-1576086213369-97a306d36557?w=400", "link": "#", "bandeira": "🔬"},
        {"titulo": "H2S donors show promise in detrusor relaxation", "fonte": "ScienceDirect", "img": "https://images.unsplash.com/photo-1532187863486-abf9dbad1b69?w=400", "link": "#", "bandeira": "💊"},
        {"titulo": "Piezo1 channels: The future of mechanotransduction", "fonte": "Cell", "img": "https://images.unsplash.com/photo-1530026405186-ed1f139313f8?w=400", "link": "#", "bandeira": "⚡"}
    ]