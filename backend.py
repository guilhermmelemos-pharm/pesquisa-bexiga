# backend.py
import streamlit as st
from Bio import Entrez, Medline
import requests
import json
import re
import time
from collections import Counter
from difflib import SequenceMatcher
from typing import List, Dict, Any, Optional

# ================= CONFIGURAÇÃO =================
Entrez.email = "pesquisador_guest@unifesp.br"

# Lista atualizada com base NOS SEUS TESTES (priorizando velocidade e inteligência de texto)
MODELOS_ATIVOS = [
    "gemini-2.5-flash",          # Prioridade 1: Novo, rápido e inteligente
    "gemini-2.0-flash",          # Prioridade 2: Estável e muito capaz
    "gemini-2.0-flash-lite",     # Prioridade 3: Ultra rápido para volumes grandes
    "gemini-exp-1206",           # Experimental de alta performance
    "gemini-flash-latest"        # Fallback para a versão estável mais recente
]

MAPA_SINONIMOS_BASE = {
    "BLADDER": "(Bladder OR Vesical OR Detrusor OR Urothelium)",
    "PAIN": "(Pain OR Nociception OR Analgesia)",
    "INFLAMMATION": "(Inflammation OR Cytokines OR Inflammasome)"
}

# ================= GEMINI CORE (INTELIGÊNCIA ARTIFICIAL) =================

def clean_model_name(model_name: str) -> str:
    """Remove o prefixo 'models/' se existir, para a URL da API ficar correta."""
    return model_name.replace("models/", "")

def call_gemini_json(prompt: str, api_key: str, temperature: float = 0.1) -> List[str]:
    """
    Chamada à API Gemini forçando resposta JSON para extração de dados.
    """
    if not api_key:
        return []

    headers = {"Content-Type": "application/json"}
    
    # Payload configurado para JSON
    payload = {
        "contents": [{"parts": [{"text": prompt}]}], 
        "generationConfig": {
            "temperature": temperature,
            "response_mime_type": "application/json"
        }
    }
    
    for modelo_raw in MODELOS_ATIVOS:
        modelo = clean_model_name(modelo_raw)
        try:
            # URL Padrão do Google AI Studio / Vertex
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={api_key.strip()}"
            
            resp = requests.post(url, headers=headers, json=payload, timeout=25)
            
            if resp.status_code == 200:
                try:
                    candidates = resp.json().get("candidates", [])
                    if not candidates: continue
                    
                    content_text = candidates[0]["content"]["parts"][0]["text"]
                    
                    # Limpeza extra caso o modelo devolva markdown ```json ... ```
                    clean_text = re.sub(r"```json|```", "", content_text).strip()
                    parsed = json.loads(clean_text)
                    
                    # Normaliza a saída para garantir que seja uma lista de strings
                    if isinstance(parsed, list):
                        return [str(x) for x in parsed]
                    elif isinstance(parsed, dict):
                        # Se devolver {"alvos": [...]}, tenta pegar a primeira lista que achar
                        for val in parsed.values():
                            if isinstance(val, list): return [str(x) for x in val]
                    
                    return [] # Se não for lista nem dict com lista
                except:
                    continue # Falha no parse JSON, tenta próximo modelo
            
            elif resp.status_code == 429:
                time.sleep(1.5) # Rate limit, espera um pouco
                continue
                
        except Exception:
            continue # Erro de conexão, tenta próximo
            
    return []

def simple_gemini_text(prompt: str, api_key: str) -> str:
    """Chamada simples para análise de texto livre (resumos/chat)."""
    if not api_key: return ""
    headers = {"Content-Type": "application/json"}
    payload = {
        "contents": [{"parts": [{"text": prompt}]}], 
        "generationConfig": {"temperature": 0.2}
    }
    
    for modelo_raw in MODELOS_ATIVOS:
        modelo = clean_model_name(modelo_raw)
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{modelo}:generateContent?key={api_key.strip()}"
            resp = requests.post(url, headers=headers, json=payload, timeout=20)
            if resp.status_code == 200:
                candidates = resp.json().get("candidates", [])
                if candidates:
                    return candidates[0]["content"]["parts"][0]["text"]
        except: continue
    return ""

# ================= EXTRAÇÃO DE CONHECIMENTO (MINERAÇÃO) =================

def ner_extraction_batch(artigos: List[Dict], api_key: str) -> List[str]:
    """
    Usa IA para ler os abstracts e extrair entidades farmacológicas.
    """
    if not artigos: return []
    
    # Cria um texto único com bullet points para economizar chamadas
    texto_amostra = [f"- {a['texto']}" for a in artigos[:35]] 
    texto_input = "\n".join(texto_amostra)
    
    prompt = f"""
    Role: Senior Molecular Pharmacologist.
    Task: Extract a JSON list of MOLECULAR TARGETS (Receptors, Enzymes) and DRUGS mentioned.
    
    STRICT EXCLUSIONS (Ignore these):
    - Clinical: MRI, CT, PET, BCG, TURBT, GWAS, RBE, UTI.
    - Biology: DNA, RNA, ATP, Gene, Protein, Cell, Tissue, Mouse, Rat.
    - Acronyms: BLCA, NMIBC, OAB, LUTS, BPS.
    
    Output Format: JSON list of strings. Example: ["P2X3", "Mirabegron", "TRPV1"]
    
    TEXT:
    {texto_input}
    """
    
    return call_gemini_json(prompt, api_key)

@st.cache_data(ttl=3600)
def minerar_pubmed(termo_base: str, email: str, usar_ia: bool = True) -> Dict:
    """
    Core do sistema: Busca no PubMed, parser robusto e extração (IA ou Regex).
    """
    api_key = st.session_state.get("api_key_usuario", "").strip()
    Entrez.email = email
    
    termo_expandido = MAPA_SINONIMOS_BASE.get(termo_base.upper(), termo_base)
    query = f"({termo_expandido}) AND (2019:2026[Date])"
    
    try:
        # 1. Busca IDs
        handle = Entrez.esearch(db="pubmed", term=query, retmax=100)
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]: return {}

        # 2. Download e Parsing com Bio.Medline (Mais seguro)
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        
        artigos_parsed = []
        for r in records:
            # Junta título, abstract e keywords
            texto = f"{r.get('TI', '')} {r.get('AB', '')} {' '.join(r.get('OT', []))}"
            if len(texto) > 50:
                artigos_parsed.append({"titulo": r.get('TI', ''), "texto": texto})
        
        entidades = []
        
        # 3. Extração: Via IA (Se permitido e tiver chave)
        if api_key and usar_ia:
            entidades = ner_extraction_batch(artigos_parsed, api_key)
        
        # 4. Fallback (Regex) - Se IA falhar ou estiver desligada
        if not entidades:
            texto_full = " ".join([a['texto'] for a in artigos_parsed])
            # Regex para PascalCase e siglas maiúsculas
            candidatos = re.findall(r'\b[A-Z][a-zA-Z0-9-]{2,15}\b', texto_full)
            
            blacklist = {
                "MRI", "CT", "PET", "BCG", "DNA", "RNA", "ATP", "GWAS", "TURP", "NHANES",
                "LUTS", "OAB", "BPH", "IPSS", "AUA", "EAU", "NMIBC", "MIBC", "TURBT",
                "SBRT", "IMPT", "RBE", "BLCA", "T24", "HLA", "CALR", "SGLT2", "VIII",
                "STUDY", "ROLE", "CELL", "TISSUE", "PATIENT", "AND", "WITH", "THE", "FOR",
                "ABSTRACT", "BACKGROUND", "METHODS", "RESULTS", "CONCLUSION", "DATA"
            }
            entidades = [e for e in candidatos if e.upper() not in blacklist]

        # 5. Estatística e Limpeza
        counts = Counter(entidades)
        # Filtra termos muito raros (ruído)
        min_freq = 2 if len(entidades) > 60 else 1
        recorrentes = [e for e, c in counts.items() if c >= min_freq]
        
        if not recorrentes: 
            recorrentes = [e for e, _ in counts.most_common(20)]
        
        # Deduplicação (ex: remover duplicatas com case diferente)
        canonicos = []
        for e in recorrentes:
            if not any(SequenceMatcher(None, e.lower(), c.lower()).ratio() > 0.90 for c in canonicos):
                canonicos.append(e)

        return {
            "termos_indicados": canonicos,
            "counts": counts,
            "total_docs": len(artigos_parsed),
            "artigos_originais": artigos_parsed
        }
    except Exception:
        return {}

# ================= FUNÇÕES DE APOIO AO FRONTEND =================

def buscar_alvos_emergentes_pubmed(alvo: str, email: str, usar_ia: bool = True) -> List[str]:
    """Wrapper para o botão 'Mágica'."""
    res = minerar_pubmed(alvo, email, usar_ia=usar_ia)
    return res.get("termos_indicados", [])

def consultar_pubmed_count(termo: str, contexto: str, email: str, ano_ini: int, ano_fim: int) -> int:
    """Conta artigos para calcular o Score de Inovação."""
    Entrez.email = email
    query = f"({termo}) AND ({contexto}) AND ({ano_ini}:{ano_fim}[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle)
        handle.close()
        return int(record["Count"])
    except:
        return 0

@st.cache_data(ttl=3600)
def buscar_resumos_detalhados(termo: str, orgao: str, email: str, ano_ini: int, ano_fim: int) -> List[Dict]:
    """Busca abstracts para leitura humana."""
    Entrez.email = email
    query = f"({termo}) AND ({orgao}) AND ({ano_ini}:{ano_fim}[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=8)
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        
        artigos = []
        for r in records:
            artigos.append({
                "Title": r.get("TI", "Sem Título"), 
                "Info_IA": r.get("AB", "Resumo indisponível.")[:900], 
                "Link": f"https://pubmed.ncbi.nlm.nih.gov/{r.get('PMID', '')}/"
            })
        return artigos
    except: return []

def analisar_abstract_com_ia(titulo: str, dados_curtos: str, api_key: str, lang: str = 'pt') -> str:
    """Analisa um único paper quando o usuário pede."""
    if not api_key: return "Insira sua API Key."
    
    prompt = f"""
    Analyze this abstract.
    Title: {titulo}
    Context: {dados_curtos}
    
    Output strictly one line: "TARGET | DRUG | EFFECT"
    """
    res = simple_gemini_text(prompt, api_key)
    if not res: return "Análise indisponível."
    return res.replace("\n", " ").strip()

def buscar_todas_noticias(lang_code: str) -> List[Dict]:
    """Radar de novidades na Home."""
    Entrez.email = "pesquisador_guest@unifesp.br"
    # Query genérica para popular o radar
    query = "(molecular pharmacology) AND (bladder OR urothelium) AND (2024:2026[Date])"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5, sort="pub_date")
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]: return []

        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        records = Medline.parse(handle)
        
        news = []
        for r in records:
            news.append({
                "titulo": r.get("TI", "Novo Artigo"), 
                "fonte": r.get("JT", "Journal")[:30], 
                "link": f"https://pubmed.ncbi.nlm.nih.gov/{r.get('PMID', '')}/"
            })
        return news
    except: return []
