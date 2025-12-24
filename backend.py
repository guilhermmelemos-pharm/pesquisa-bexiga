import requests
import json
import re
import ast
from Bio import Entrez
from collections import Counter
from tenacity import retry, stop_after_attempt, wait_exponential

# ================= CONFIG =================
Entrez.email = "pesquisador_guest@unifesp.br"

MODELOS_ATIVOS = [
    "models/gemini-2.5-flash",
    "models/gemini-2.0-flash",
    "models/gemini-2.0-flash-exp",
    "models/gemini-flash-latest",
    "models/gemini-2.5-pro"
]

# ================= GEMINI =================
def montar_url(modelo, chave):
    return f"https://generativelanguage.googleapis.com/v1beta/{modelo}:generateContent?key={chave}"

def _chamar_gemini(prompt, api_key, temperatura=0.0, timeout=20):
    headers = {"Content-Type": "application/json"}
    data = {
        "contents": [{"parts": [{"text": prompt}]}],
        "generationConfig": {"temperature": temperatura}
    }

    ultimo_erro = None

    for modelo in MODELOS_ATIVOS:
        try:
            r = requests.post(
                montar_url(modelo, api_key),
                headers=headers,
                data=json.dumps(data),
                timeout=timeout
            )

            if r.status_code != 200:
                ultimo_erro = f"{modelo} → HTTP {r.status_code}"
                continue

            js = r.json()
            texto = js["candidates"][0]["content"]["parts"][0]["text"]
            if texto:
                return texto.strip()

        except Exception as e:
            ultimo_erro = f"{modelo} → {e}"
            continue

    raise RuntimeError(f"Falha Gemini: {ultimo_erro}")

# ================= IA – FAXINA =================
def faxina_ia(lista_suja, api_key):
    if not api_key or not lista_suja:
        return lista_suja[:30]

    prompt = f"""
    ROLE: Senior Scientist in Pharmacology & Physiopathology.
    TASK: From the list below, keep ONLY:
    - Pharmacological targets (receptors, channels, enzymes)
    - Specific drugs or bioactive molecules
    REMOVE:
    - Clinical terms
    - Anatomy
    - General biology
    - Study words

    RETURN: a pure Python list of strings.

    INPUT LIST:
    {", ".join(lista_suja)}
    """

    try:
        resposta = _chamar_gemini(prompt, api_key, temperatura=0.0)
        resposta = resposta.replace("```python", "").replace("```", "").strip()
        lista_limpa = ast.literal_eval(resposta)

        if isinstance(lista_limpa, list) and lista_limpa:
            return lista_limpa

    except Exception:
        pass

    return lista_suja[:30]

# ================= IA – ABSTRACT =================
def analisar_abstract_com_ia(titulo, dados, api_key, lang="pt"):
    idioma = "Portuguese" if lang == "pt" else "English"

    prompt = f"""
    Summarize the main target/drug and mechanism in 15 words.
    Title: {titulo}
    Context: {dados}
    Language: {idioma}
    """

    try:
        return _chamar_gemini(prompt, api_key, temperatura=0.2, timeout=10)
    except Exception as e:
        return f"⚠️ Erro IA: {e}"

# ================= PUBMED =================
@retry(stop=stop_after_attempt(3), wait=wait_exponential())
def _pubmed_count(query):
    h = Entrez.esearch(db="pubmed", term=query, retmax=0)
    r = Entrez.read(h)
    h.close()
    return int(r["Count"])

def consultar_pubmed_count(termo, contexto, email, ano_ini, ano_fim):
    if email:
        Entrez.email = email

    query = f"({termo}) AND ({ano_ini}:{ano_fim}[DP]) AND NOT Review[pt]"
    try:
        return _pubmed_count(query)
    except:
        return 0

# ================= MINERAÇÃO =================
def buscar_alvos_emergentes_pubmed(termo_base, email, api_key=None, usar_ia=True):
    if email:
        Entrez.email = email

    query = f"{termo_base}[Title/Abstract] AND (2018:2030[DP]) AND NOT Review[pt]"

    h = Entrez.esearch(db="pubmed", term=query, retmax=2000, sort="relevance")
    r = Entrez.read(h)
    h.close()

    if not r["IdList"]:
        return []

    h = Entrez.efetch(db="pubmed", id=r["IdList"], rettype="medline", retmode="text")
    raw = h.read()
    h.close()

    candidatos = []

    for art in raw.split("\n\nPMID-"):
        texto = " ".join(
            l[6:] for l in art.split("\n")
            if l.startswith(("TI  - ", "KW  - ", "OT  - "))
        )

        encontrados = re.findall(r"\b[A-Z]{2,}[A-Z0-9-]*\b", texto)
        candidatos.extend(encontrados)

    if not candidatos:
        return []

    freq = Counter(candidatos)
    total = max(1, len(raw))

    filtrados = [
        k for k, v in freq.items()
        if 2 <= len(k) <= 20 and (v / total) < 0.9
    ]

    if usar_ia and api_key:
        return faxina_ia(filtrados, api_key)

    return filtrados[:40]
