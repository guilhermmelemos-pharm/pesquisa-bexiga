import streamlit as st
from Bio import Entrez
import google.generativeai as genai
import re
from collections import Counter
import time

# --- CONFIGURAÇÃO ---
Entrez.email = "pesquisador_guest@unifesp.br" 

def analisar_abstract_com_ia(titulo, abstract, api_key, lang='pt'):
    if not api_key: return "⚠️ IA não ativada"
    
    try:
        genai.configure(api_key=api_key)
        idioma = "Português" if lang == 'pt' else "Inglês"
        
        # LÓGICA DE INFERÊNCIA SÊNIOR: Se o texto for curto ou nulo, foca no título
        abstract_txt = abstract if (abstract and len(abstract) > 30) else "Abstract indisponível. Infira pelo título."
        
        prompt = f"""Como PhD em Farmacologia, analise tecnicamente:
        TITULO: {titulo}
        RESUMO: {abstract_txt[:3000]}
        
        AÇÃO: Mesmo que o resumo esteja incompleto ou ausente, use o TITULO para inferir o mecanismo farmacológico provável.
        FORMATO OBRIGATÓRIO: Alvo → Fármaco/Substância → Efeito funcional na bexiga.
        REGRAS: Máximo 30 palavras. Idioma: {idioma}."""

        # LISTA DE MODELOS PARA TENTATIVAS (Evita Erro 404 e 429)
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
                erro_str = str(e)
                if "429" in erro_str:
                    time.sleep(1)
                continue

        # Fallback de segurança se todos os modelos falharem (Hardcoded para seus temas)
        t_upper = titulo.upper()
        if "ENDOCANNABINOID" in t_upper or "CB1" in t_upper:
            return "CB1/CB2 → Agonistas/Inibidores FAAH → Antinocicepção e redução de hiperatividade vesical."
        if "ROS" in t_upper or "OXIDATIVE" in t_upper:
            return "Enzimas Antioxidantes (SOD/CAT) → Inibidores NOX/Antioxidantes → Redução de fibrose e estresse oxidativo."
            
        return f"💡 Inferência: {titulo} → (Modelos ocupados, tente em 1 min)."

    except Exception as e:
        return f"❌ Erro Crítico: {str(e)[:30]}"

def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    query = f"({termo}) AND ({contexto}) AND ({ano_ini}:{ano_fim}[Date - Publication]) AND (NOT Review[pt])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle)
        handle.close()
        return int(record["Count"])
    except: return 0

def buscar_resumos_detalhados(termo, orgao, email, ano_ini, ano_fim):
    if email: Entrez.email = email
    query = f"({termo}) AND ({orgao}) AND ({ano_ini}:{ano_fim}[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=3)
        record = Entrez.read(handle)
        handle.close()
        ids = record["IdList"]
        if not ids: return []
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        artigos = []
        for raw in dados.split("\n\nPMID-"):
            tit, pmid, abs_txt, last_tag = "", "", "", ""
            for line in raw.split("\n"):
                if line.strip().isdigit() and not pmid: pmid = line.strip()
                if line.startswith("TI  - "): tit = line.replace("TI  - ", "").strip(); last_tag = "TI"
                elif line.startswith("      ") and last_tag == "TI": tit += " " + line.strip()
                if line.startswith("AB  - "): abs_txt = line.replace("AB  - ", "").strip(); last_tag = "AB"
                elif line.startswith("      ") and last_tag == "AB": abs_txt += " " + line.strip()
            if tit: artigos.append({"Title": tit, "Resumo_Original": abs_txt, "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return artigos
    except: return []

def buscar_alvos_emergentes_pubmed(termo_base, email):
    if email: Entrez.email = email
    query = f"({termo_base}) AND (receptor OR signaling OR channel) AND (2023:2026[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=30)
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        candidatos = re.findall(r'\b[A-Z][A-Z0-9]{2,8}\b', dados)
        blacklist = {"THE", "AND", "FOR", "PMID", "PMC", "DOI", "USA", "TYPE", "CELL", "ROLE", "DATA"}
        validos = [c for c in candidatos if c not in blacklist and c not in termo_base.upper()]
        return [item for item, count in Counter(validos).most_common(7)]
    except: return []

def buscar_todas_noticias(lang='pt'):
    try:
        query = "(bladder OR detrusor) AND (pharmacology OR oxidative stress OR mechanotransduction)"
        handle = Entrez.esearch(db="pubmed", term=query, retmax=3, sort="pub_date")
        record = Entrez.read(handle); handle.close()
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        dados = handle.read(); handle.close()
        news = []
        for art in dados.split("\n\nPMID-"):
            tit = ""; pmid = ""
            for line in art.split("\n"):
                if line.startswith("TI  - "): tit = line.replace("TI  - ", "").strip()
                if line.strip().isdigit() and not pmid: pmid = line.strip()
            if tit: news.append({"titulo": tit, "fonte": "PubMed", "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return news
    except: return []
