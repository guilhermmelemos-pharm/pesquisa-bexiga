# backend.py
import streamlit as st
from Bio import Entrez
import requests
import re
import time
import ast
import json
from collections import Counter
from difflib import SequenceMatcher

# ================= CONFIG =================
Entrez.email = "pesquisador_guest@unifesp.br"

MODELOS_ATIVOS = [
    "gemini-2.0-flash",
    "gemini-2.0-flash-exp"
]

MAPA_SINONIMOS_BASE = {
    "BLADDER": "(Bladder OR Vesical OR Detrusor OR Urothelium)",
    "PAIN": "(Pain OR Nociception OR Analgesia)",
    "INFLAMMATION": "(Inflammation OR Cytokines OR Inflammasome)"
}

# ================= GEMINI CORE =================

def montar_url(modelo, chave):
    return f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={chave.strip()}"

def call_gemini(prompt, api_key, temperature=0.0):
    headers = {"Content-Type": "application/json"}
    payload = {
        "contents": [{"parts": [{"text": prompt}]}],
        "generationConfig": {"temperature": temperature}
    }

    last_error = ""
    for modelo in MODELOS_ATIVOS:
        try:
            resp = requests.post(
                montar_url(modelo, api_key),
                headers=headers,
                json=payload,
                timeout=25
            )
            if resp.status_code == 200:
                return resp.json()["candidates"][0]["content"]["parts"][0]["text"]
            elif resp.status_code == 429:
                time.sleep(2)
                continue
            else:
                last_error = f"API {resp.status_code}: {resp.text}"
        except Exception as e:
            last_error = str(e)
            continue
    
    if last_error:
        st.error(f"Falha na comunicação com Gemini: {last_error}")
    return ""

# ================= ROBUST PARSERS =================

def clean_json_text(text):
    """Remove markdown code blocks and common noise before parsing."""
    if not text: return ""
    text = re.sub(r'```(?:json|python|list)?\n?', '', text)
    text = text.replace('```', '').strip()
    return text

def parse_structure(text, expected_type=list):
    """Generic robust parser for JSON and Python literals."""
    cleaned = clean_json_text(text)
    if not cleaned: return expected_type()
    
    # Tentativa 1: JSON puro
    try:
        data = json.loads(cleaned)
        if isinstance(data, expected_type): return data
    except: pass

    # Tentativa 2: ast.literal_eval (para listas/dicts Python)
    try:
        data = ast.literal_eval(cleaned)
        if isinstance(data, expected_type): return data
    except: pass

    # Tentativa 3: Regex para encontrar o bloco estruturado
    try:
        pattern = r'\[.*\]' if expected_type == list else r'\{.*\}'
        match = re.search(pattern, cleaned, re.DOTALL)
        if match:
            data = ast.literal_eval(match.group(0))
            if isinstance(data, expected_type): return data
    except: pass

    return expected_type()

# ================= QUERY INTELLIGENCE =================

def expandir_termo_com_ia(termo, api_key):
    prompt = f"""
    Expand this biomedical concept into a concise PubMed Boolean query with synonyms.
    Ensure parentheses are balanced and OR/AND operators are correct.
    Return ONLY the query string, no quotes.
    TERM: {termo}
    """
    r = call_gemini(prompt, api_key, temperature=0.0)
    query = r.replace('"', '').strip() if r else termo
    
    # Validação mínima de parênteses
    if query.count('(') != query.count(')'):
        return termo
    return query

# ================= NORMALIZAÇÃO =================

def normalize_locally(entities):
    canonical = []
    for e in entities:
        if not any(SequenceMatcher(None, e.lower(), c.lower()).ratio() > 0.85 for c in canonical):
            canonical.append(e)
    return canonical

# ================= BATCH NER =================

def ner_extraction_batch(artigos, api_key, batch_size=15):
    all_entities = []
    for i in range(0, len(artigos), batch_size):
        batch = artigos[i:i + batch_size]
        texto_input = ""
        for art in batch:
            texto_input += f"- {art['texto']}\n"

        prompt = f"""
        You are a PhD-level biocurator. Extract molecular targets (receptors, channels) and experimental drugs.
        Return ONLY a JSON list of strings. No explanations.
        TEXT:
        {texto_input}
        """
        raw = call_gemini(prompt, api_key, temperature=0.1)
        ents = parse_structure(raw, list)
        all_entities.extend([str(e).strip() for e in ents if len(str(e)) > 2])
    return all_entities

# ================= BATCH EFFECT INFERENCE =================

def inferir_efeitos_em_lote(artigos_selecionados, api_key):
    if not artigos_selecionados: return []
    
    prompt = """
    Analyze the following PubMed titles. For each, identify the Target, the Drug, and the specific Pharmacological Effect.
    Return strictly a JSON list of objects: [{"title": "...", "target": "...", "drug": "...", "effect": "..."}]
    If unclear, use "N/A". 
    TITLES:
    """
    for a in artigos_selecionados:
        prompt += f"- {a['titulo']}\n"

    raw = call_gemini(prompt, api_key, temperature=0.0)
    return parse_structure(raw, list)

# ================= CORE PIPELINE =================

@st.cache_data(ttl=3600)
def minerar_pubmed(termo_base, email, api_key):
    if not email: raise ValueError("Email obrigatório.")
    if not api_key: return {}

    Entrez.email = email
    termo = MAPA_SINONIMOS_BASE.get(termo_base.upper(), termo_base)
    query_ia = expandir_termo_com_ia(termo, api_key)
    query = f"({query_ia}) AND (2020:2026[Date - Publication])"

    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=200)
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return {}

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        raw = handle.read(); handle.close()

        artigos = []
        for bloco in raw.split("\n\nPMID-"):
            titulo, texto = "", ""
            for line in bloco.split("\n"):
                if line.startswith("TI  - "): titulo = line[6:].strip()
                if line.startswith(("TI  - ", "KW  - ", "OT  - ")):
                    texto += line[6:].strip() + " "
            if texto: artigos.append({"titulo": titulo, "texto": texto})

        # NER Batch
        all_entities = ner_extraction_batch(artigos, api_key)
        counts = Counter(all_entities)
        
        # Seleção de alvos para inferência baseada em recorrência (Top Overlap)
        alvos_recorrentes = [e for e, c in counts.most_common(20)]
        
        # Filtro de títulos mais informativos: aqueles que contêm os alvos mais recorrentes
        titulos_informativos = []
        for art in artigos:
            score = sum(1 for alvo in alvos_recorrentes if alvo.lower() in art['texto'].lower())
            titulos_informativos.append((score, art))
        
        # Ordena por score e pega os 12 melhores para amostragem
        artigos_amostragem = [art for score, art in sorted(titulos_informativos, key=lambda x: x[0], reverse=True)[:12]]

        # Processamento final
        recorrentes_final = [e for e, c in counts.items() if c >= 2]
        if not recorrentes_final: recorrentes_final = [e for e, c in counts.most_common(15)]
        
        alvos_final = normalize_locally(recorrentes_final)
        inferencias = inferir_efeitos_em_lote(artigos_amostragem, api_key)

        return {
            "termos_indicados": alvos_final,
            "inferencias": inferencias
        }
    except Exception as e:
        st.error(f"Erro no pipeline: {e}")
        return {}

# ================= NEWS RADAR =================

@st.cache_data(ttl=3600)
def buscar_todas_noticias(email):
    if not email: return []
    Entrez.email = email
    query = "(molecular pharmacology OR ion channels OR signaling pathway) AND (2024:2026[Date - Publication])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=6, sort="pub_date")
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return []
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        raw = handle.read(); handle.close()
        news = []
        for bloco in raw.split("\n\nPMID-"):
            tit, journal, pmid = "", "", ""
            for line in bloco.split("\n"):
                if line.startswith("TI  - "): tit = line[6:].strip()
                if line.startswith("JT  - "): journal = line[3:].strip()
                if line.strip().isdigit() and not pmid: pmid = line.strip()
            if tit:
                news.append({"titulo": tit, "fonte": journal[:35], "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        return news
    except: return []
