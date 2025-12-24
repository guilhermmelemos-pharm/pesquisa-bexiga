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
    if not api_key:
        return ""

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
                last_error = f"Status {resp.status_code}"
        except Exception as e:
            last_error = str(e)
            continue
    
    if last_error:
        st.sidebar.warning(f"IA instável ({last_error}). Tentando fallback...")
    return ""

# ================= PARSERS "ANTI-LIXO" =================

def clean_text(text):
    if not text: return ""
    # Remove blocos de código markdown que a IA adora colocar
    text = re.sub(r"```.*?```", "", text, flags=re.DOTALL)
    return text.strip()

def parse_structure(text, expected=list):
    """Parser robusto: tenta JSON, AST e, por fim, Regex Sniper."""
    text = clean_text(text)
    if not text: return expected()

    # Tentativa 1: JSON/AST puro
    try:
        # Tenta extrair o que está entre colchetes ou chaves
        pattern = r"\[.*\]" if expected == list else r"\{.*\}"
        match = re.search(pattern, text, re.DOTALL)
        if match:
            obj = ast.literal_eval(match.group(0))
            if isinstance(obj, expected): return obj
    except:
        pass

    # Tentativa 2: Regex Sniper (Se a IA falhar na sintaxe, pegamos os nomes na força bruta)
    if expected == list:
        # Busca siglas e nomes técnicos entre aspas ou soltos
        names = re.findall(r"['\"]?([A-Z][A-Z0-9-]{2,15})['\"]?", text)
        if names:
            # Remove lixo óbvio que o regex pode pegar
            blacklist = {"TITLE", "TARGET", "DRUG", "EFFECT", "JSON", "LIST", "TRUE", "FALSE", "PMID"}
            return list(set([n for n in names if n not in blacklist]))

    return expected()

# ================= INTELLIGENCE =================

def expandir_termo_com_ia(termo, api_key):
    prompt = f"Expand this biomedical concept into a PubMed Boolean query (synonyms only). Return ONLY the query string. TERM: {termo}"
    out = call_gemini(prompt, api_key, temperature=0.0)
    out = out.replace('"', '').strip()
    # Valida parênteses: se estiver quebrado, usa o original para não dar erro no PubMed
    if not out or out.count("(") != out.count(")"):
        return termo
    return out

def normalize_entities(entities):
    canon = []
    for e in entities:
        if not any(SequenceMatcher(None, e.lower(), c.lower()).ratio() > 0.85 for c in canon):
            canon.append(e)
    return canon

# ================= BATCH ENGINES =================

def ner_extraction_batch(artigos, api_key, batch_size=15):
    all_entities = []
    for i in range(0, len(artigos), batch_size):
        bloco = artigos[i:i+batch_size]
        texto_input = "\n".join([f"- {a['texto']}" for a in bloco])

        prompt = f"""
        Extract MOLECULAR TARGETS and EXPERIMENTAL DRUGS from the text below.
        Ignore diseases, organs, and clinical scores.
        Return ONLY a Python list of strings.
        TEXT:
        {texto_input}
        """
        raw = call_gemini(prompt, api_key, temperature=0.1)
        ents = parse_structure(raw, list)
        all_entities.extend(ents)
    return all_entities

def inferir_efeitos_em_lote(artigos, api_key):
    if not artigos: return []
    prompt = """
    Infer pharmacological relationship from TITLES.
    Return ONLY a JSON list of objects: [{"title": "...", "target": "...", "drug": "...", "effect": "..."}]
    TITLES:
    """
    for a in artigos: prompt += f"- {a['titulo']}\n"
    
    raw = call_gemini(prompt, api_key, temperature=0.0)
    return parse_structure(raw, list)

# ================= CORE PIPELINE =================

@st.cache_data(ttl=3600)
def minerar_pubmed(termo_base, email):
    if not email: raise ValueError("Email obrigatório.")
    api_key = st.session_state.get("api_key_usuario", "").strip()
    
    Entrez.email = email
    termo = MAPA_SINONIMOS_BASE.get(termo_base.upper(), termo_base)
    
    # Se houver API Key, tenta expandir, senão usa base
    if api_key:
        termo = expandir_termo_com_ia(termo, api_key)

    query = f"({termo}) AND (2020:2026[Date - Publication])"

    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=150)
        record = Entrez.read(handle); handle.close()
        if not record["IdList"]: return {}

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        raw_medline = handle.read(); handle.close()

        artigos = []
        for bloco in raw_medline.split("\n\nPMID-"):
            titulo, texto = "", ""
            for line in bloco.split("\n"):
                if line.startswith("TI  - "): titulo = line[6:].strip()
                if line.startswith(("TI  - ", "KW  - ", "OT  - ")):
                    texto += line[6:].strip() + " "
            if texto: artigos.append({"titulo": titulo, "texto": texto})

        # --- EXTRAÇÃO (NER) ---
        entidades = []
        if api_key:
            entidades = ner_extraction_batch(artigos, api_key)
        
        if not entidades:
            # Fallback determinístico caso a IA falhe completamente
            texto_bruto = " ".join([a['texto'] for a in artigos])
            entidades = re.findall(r'\b[A-Z0-9-]{3,15}\b', texto_bruto)

        contagem = Counter(entidades)

        # 🔥 TRAVA ANTI-VAZIO: Prioriza recorrência, mas aceita frequência bruta
        recorrentes = [e for e, c in contagem.items() if c >= 2]
        if not recorrentes:
            recorrentes = [e for e, _ in contagem.most_common(15)]

        alvos_final = normalize_entities(recorrentes)

        # --- INFERÊNCIA (AMOSTRAGEM INTELIGENTE) ---
        # Ordena artigos que possuem os alvos detectados (melhor overlap)
        artigos_rankeados = sorted(
            artigos, 
            key=lambda x: sum(1 for alvo in alvos_final if alvo.lower() in x['texto'].lower()), 
            reverse=True
        )
        
        inferencias = []
        if api_key:
            inferencias = inferir_efeitos_em_lote(artigos_rankeados[:12], api_key)

        return {
            "termos_indicados": alvos_final,
            "inferencias": inferencias
        }
    except Exception as e:
        st.error(f"Erro crítico: {e}")
        return {}

# ================= NEWS & UTILS =================

def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    Entrez.email = email
    query = f"({termo}) AND ({ano_ini}:{ano_fim}[Date - Publication])"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
    record = Entrez.read(handle); handle.close()
    return int(record["Count"])

@st.cache_data(ttl=3600)
def buscar_todas_noticias(email):
    if not email: return []
    Entrez.email = email
    query = "(molecular pharmacology OR ion channels OR signaling) AND (2024:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=6, sort="pub_date")
        record = Entrez.read(handle); handle.close()
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

# Função ponte para o front antigo não quebrar
def buscar_alvos_emergentes_pubmed(alvo, email, usar_ia=True):
    res = minerar_pubmed(alvo, email)
    return res.get("termos_indicados", [])
