"""
Lemos Lambda v2.1 — Backend
Lógica INVERTIDA: busca em campos distantes → Flash faz ponte com o órgão
Inspirado em: Trehalose/TMAO descobertos por cruzamento lateral de literatura
"""

import re, json, time, requests
from typing import List, Dict, Generator, Callable
from Bio import Entrez, Medline
import constantes as c

# ─── Campos de busca lateral (distantes do órgão alvo) ───────────────────────
# Esses campos geram candidatos como Trehalose, TMAO, Spermidine, Urolithin A
CAMPOS_LATERAIS = [
    # Proteostase & Chaperonas
    ("proteostasis chaperone stress",           "chaperone OR proteostasis OR unfolded protein response OR ER stress"),
    ("autophagy inducer compound",              "autophagy OR mitophagy OR LC3 OR beclin OR mTOR inhibitor"),
    ("osmolyte cytoprotection",                 "osmolyte OR trehalose OR TMAO OR betaine OR taurine OR cytoprotection"),
    # Microbioma & Metabólitos
    ("microbiome metabolite systemic effect",   "microbiome metabolite OR butyrate OR propionate OR urolithin OR indole OR SCFA"),
    ("gut derived compound organ protection",   "gut-derived OR postbiotic OR microbial metabolite AND organ protection"),
    # Longevidade & Senescência
    ("senolytic senomorphic compound",          "senolytic OR senomorphic OR senostatic OR cellular senescence AND drug"),
    ("longevity compound aging pathway",        "longevity OR lifespan extension OR rapamycin OR spermidine OR NAD+ OR sirtuin activator"),
    # Neuroproteção (compostos com efeito sistêmico)
    ("neuroprotective compound translational",  "neuroprotective compound OR sigma-1 receptor OR neurosteroid OR BDNF modulator"),
    # Resolução de inflamação
    ("specialized pro-resolving mediator",      "resolvin OR maresin OR lipoxin OR protectin OR pro-resolving OR SPM"),
    ("anti-inflammatory natural compound",      "natural compound anti-inflammatory OR polyphenol OR flavonoid OR terpenoid AND mechanism"),
    # Gasotransmissores
    ("gasotransmitter H2S CO NO donor",         "hydrogen sulfide donor OR CO-releasing molecule OR nitric oxide donor OR gasotransmitter"),
    # Mecanossensação
    ("mechanosensing ion channel epithelium",   "Piezo1 OR Piezo2 OR TRPV4 OR TRPM7 OR mechanosensitive channel AND epithelium"),
    # Metabolismo mitocondrial
    ("mitochondrial uncoupler metabolic drug",  "mitochondrial OR uncoupler OR NAD precursor OR coenzyme Q OR PGC1alpha activator"),
    # Ferroptose & Estresse oxidativo
    ("ferroptosis GPX4 lipid peroxidation",     "ferroptosis OR GPX4 OR lipid peroxidation OR RSL3 OR ferrostatin"),
    # Epigenética emergente
    ("epigenetic modulator m6A RNA",            "m6A methylation OR METTL3 OR YTHDF OR RNA modification AND drug"),
    # Reposicionamento de drogas aprovadas
    ("approved drug unexpected mechanism",      "drug repurposing OR unexpected mechanism OR off-target effect AND organ"),
]

# ─── Gemini ───────────────────────────────────────────────────────────────────

def _gemini_call(api_key: str, model: str, prompt: str, max_tokens: int = 1500) -> str:
    url = f"https://generativelanguage.googleapis.com/v1beta/models/{model}:generateContent?key={api_key}"
    body = {
        "contents": [{"parts": [{"text": prompt}]}],
        "generationConfig": {"temperature": 0.2, "maxOutputTokens": max_tokens, "candidateCount": 1}
    }
    try:
        r = requests.post(url, json=body, timeout=40)
        r.raise_for_status()
        return r.json()["candidates"][0]["content"]["parts"][0]["text"].strip()
    except Exception as e:
        print(f"❌ Gemini [{model}]: {e}")
        return ""

# ─── Filtro de lixo genômico (regex, zero tokens) ────────────────────────────

def _e_lixo(termo: str) -> bool:
    for pat in c.LIXO_REGEX:
        if re.search(pat, termo, re.IGNORECASE):
            return True
    if re.search(r"[A-Z0-9]+-AS\d*$", termo): return True
    if termo.upper() in c.BLACKLIST: return True
    if len(termo) < 3 or termo.isdigit(): return True
    # Palavras comuns de texto acadêmico
    palavras_texto = {
        "multivariate","univariate","coordinate","alleviate","polymerase",
        "dehydrogenase","oxygenase","ubiquitin","carbachol","scopolamine",
        "adenylate","glutathione","bromadiolone","multivariate","cytokine",
        "chemokine","phosphatase","transferase","reductase","hydroxylase",
    }
    if termo.lower() in palavras_texto: return True
    return False

# ─── PubMed ───────────────────────────────────────────────────────────────────

def _pubmed_count(query: str, email: str) -> int:
    Entrez.email = email
    try:
        h = Entrez.esearch(db="pubmed", term=query, retmax=0)
        r = Entrez.read(h); h.close()
        return int(r["Count"])
    except: return 0

def _pubmed_search(query: str, email: str, retmax: int = 100) -> List[str]:
    Entrez.email = email
    try:
        h = Entrez.esearch(db="pubmed", term=query, retmax=retmax)
        r = Entrez.read(h); h.close()
        return r["IdList"]
    except: return []

def _pubmed_fetch_abstracts(ids: List[str], email: str) -> List[Dict]:
    if not ids: return []
    Entrez.email = email
    try:
        h = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
        records = list(Medline.parse(h))
        return [{"title": r.get("TI",""), "abstract": r.get("AB",""),
                 "pmid": r.get("PMID",""),
                 "link": f"https://pubmed.ncbi.nlm.nih.gov/{r.get('PMID','')}/"}
                for r in records if r.get("TI") or r.get("AB")]
    except: return []

def _extrair_candidatos(articles: List[Dict]) -> Dict[str, int]:
    raw: Dict[str, int] = {}
    for a in articles:
        text = f"{a.get('title','')} {a.get('abstract','')}"
        for pat in c.REGEX_PATTERNS:
            for m in re.findall(pat, text):
                if not _e_lixo(m):
                    raw[m] = raw.get(m, 0) + 1
    return raw

# ─── OpenTargets ──────────────────────────────────────────────────────────────

def _opentargets(disease: str) -> List[str]:
    query = """query($d:String!){search(queryString:$d,entityNames:["target"]){
    hits{object{...on Target{approvedSymbol}}}}}"""
    try:
        r = requests.post("https://api.platform.opentargets.org/api/v4/graphql",
                          json={"query": query, "variables": {"disease": disease}}, timeout=12)
        hits = r.json().get("data",{}).get("search",{}).get("hits",[])
        return [h["object"]["approvedSymbol"] for h in hits[:30]
                if h.get("object") and h["object"].get("approvedSymbol")]
    except: return []

# ─── ClinicalTrials reposicionamento ─────────────────────────────────────────

def _trials_reposicionamento(alvo: str) -> List[str]:
    termos_excluir = alvo.lower().split()
    candidatos = []
    try:
        r = requests.get(
            "https://clinicaltrials.gov/api/v2/studies?query.term=phase+2+drug&pageSize=25&format=json",
            headers={"Accept": "application/json"}, timeout=10
        )
        for s in r.json().get("studies", []):
            titulo = s["protocolSection"]["identificationModule"].get("briefTitle","").lower()
            if not any(t in titulo for t in termos_excluir):
                for iv in s["protocolSection"].get("armsInterventionsModule",{}).get("interventions",[]):
                    if iv.get("type") == "DRUG":
                        nome = iv.get("name","").strip()
                        if nome and len(nome) > 3 and not _e_lixo(nome):
                            candidatos.append(nome)
    except: pass
    return candidatos[:20]

# ─── HPA expressão no tecido ──────────────────────────────────────────────────

def _detectar_tecido(alvo: str) -> str:
    alvo_l = alvo.lower()
    for key, tissue in c.HPA_TISSUE_MAP.items():
        if key in alvo_l:
            return tissue
    return ""

def _hpa_expresso(gene: str, tecido: str) -> bool:
    if not tecido: return True
    try:
        r = requests.get(
            f"https://www.proteinatlas.org/api/search_download.php?search={gene}&format=json&columns=g,t&compress=no",
            timeout=6)
        if r.status_code != 200: return True
        data = r.json()
        if not data: return True
        for key, val in data[0].items():
            if tecido.lower() in key.lower() and "RNA" in key:
                try: return float(str(val).replace("<","").strip()) > 0.5
                except: return True
        return True
    except: return True

# ─── Flash: ponte mecanismo → órgão (coração do sistema) ─────────────────────

def flash_bridge(candidatos: List[str], alvo: str, api_key: str) -> Dict[str, str]:
    """
    Flash recebe candidatos de campos DISTANTES e decide:
    1. Esse composto/alvo tem mecanismo biológico plausível no órgão alvo?
    2. Já foi estudado nesse órgão? (se sim → menos interessante)
    3. Classifica o tipo

    Retorna {termo: tipo} — Junk se não faz sentido no órgão.
    Inspirado em: Trehalose → autofagia → urotélio sob distensão.
    """
    if not api_key or not candidatos:
        return {}

    prompt = (
        f'You are a senior pharmacologist and human physiologist with 15 years of bench experience. '
        f'You are hunting for completely unexplored therapeutic opportunities for "{alvo}". '
        f'You just ran a literature mining pipeline searching in DISTANT fields '
        f'(microbiome, proteostasis, longevity, neuroprotection, gasotransmitters, resolution of inflammation) '
        f'and got a list of compounds and targets — most have NEVER been tested in "{alvo}". '
        f'\n\nYour job: for each term, ask yourself — '
        f'"Does this compound or target have a plausible biological mechanism that could work in {alvo}?" '
        f'Think like a scientist who crosses literature from different fields. '
        f'Example reasoning: "Trehalose induces autophagy and reduces ER stress — '
        f'urothelium under distension experiences ER stress — this is worth testing." '
        f'Compounds testable in animal models (rodents, in vivo) are especially welcome. '
        f'Consider cross-organ translational potential: if it worked in kidney epithelium, '
        f'it might work in bladder urothelium. '
        f'Be generous: if there is ANY plausible mechanism — KEEP IT and classify it. '
        f'Be ruthless only with: statistical terms, bacterial strains, cell lines, '
        f'lncRNAs, pseudogenes, retroviral sequences, cosmetic compounds with zero pharmacology, '
        f'ambiguous acronyms with no known human target. '
        f'When in doubt — KEEP IT. '
        f'\n\nINPUT: {", ".join(candidatos[:80])}'
        f'\n\nReturn JSON only. Keys=exact term, Values=one of: '
        f'"Receptor"|"Drug"|"Modulator"|"Gene"|"Ion Channel"|"Enzyme"|"Pathway"|"Biomarker"|"Junk" '
        f'\nJSON ONLY:'
    )

    res = _gemini_call(api_key, c.MODELOS["flash"], prompt, 1500)
    try:
        clean = re.sub(r"```json|```", "", res).strip()
        s, e = clean.find("{"), clean.rfind("}")
        if s != -1 and e != -1:
            return json.loads(clean[s:e+1])
    except: pass
    return {}

# ─── Pipeline Principal ────────────────────────────────────────────────────────

def minerar_alvos(
    alvo: str, email: str, api_key: str,
    usar_ia: bool = True,
    status_cb: Callable[[str], None] = None
) -> Dict:
    """
    Lógica INVERTIDA: mecanismo → órgão
    Busca em campos distantes, Flash faz a ponte.
    """
    def log(msg):
        if status_cb: status_cb(msg)

    tecido_hpa = _detectar_tecido(alvo)
    todos: Dict[str, int] = {}
    typemap: Dict[str, str] = {}

    # ── Fase 1: Busca em campos distantes (SEM mencionar o órgão alvo) ────────
    log(f"🔭 Buscando em {len(CAMPOS_LATERAIS)} campos distantes do '{alvo}'...")
    for nome_campo, query_campo in CAMPOS_LATERAIS:
        log(f"   📚 {nome_campo}...")
        # Busca deliberadamente SEM o alvo principal
        query = f'({query_campo}) AND (2019:2026[Date]) NOT review[PT]'
        ids = _pubmed_search(query, email, 80)
        if ids:
            arts = _pubmed_fetch_abstracts(ids[:60], email)
            raw = _extrair_candidatos(arts)
            for k, v in raw.items():
                todos[k] = todos.get(k, 0) + v
            log(f"      ✅ {len(raw)} candidatos extraídos")
        time.sleep(0.35)

    log(f"📊 Total bruto: {len(todos)} candidatos únicos")

    # ── Fase 2: OpenTargets (alvos validados geneticamente) ───────────────────
    log("🎯 Consultando OpenTargets...")
    ot = _opentargets(alvo)
    for sym in ot:
        if sym and not _e_lixo(sym):
            todos[sym] = todos.get(sym, 0) + 3
    log(f"   ✅ {len(ot)} alvos do OpenTargets")

    # ── Fase 3: ClinicalTrials reposicionamento ───────────────────────────────
    log("🏥 Buscando reposicionamento (trials em outras indicações)...")
    drugs = _trials_reposicionamento(alvo)
    for d in drugs:
        todos[d] = todos.get(d, 0) + 2
    log(f"   ✅ {len(drugs)} compostos de trials")

    # ── Fase 4: Filtra lixo genômico (regex, zero tokens) ────────────────────
    log("🧹 Removendo lixo genômico...")
    antes = len(todos)
    todos_limpos = {k: v for k, v in todos.items() if not _e_lixo(k)}
    log(f"   ✅ {antes - len(todos_limpos)} removidos, {len(todos_limpos)} sobreviveram")

    # Ordena por frequência de aparição
    ranked = sorted(todos_limpos.keys(), key=lambda x: todos_limpos.get(x, 0), reverse=True)[:150]

    # ── Fase 5: HPA — expressão no tecido ────────────────────────────────────
    if tecido_hpa and ranked:
        log(f"🧬 Verificando expressão em '{tecido_hpa}' (HPA)...")
        confirmados, sem_dado = [], []
        for i, term in enumerate(ranked[:80]):
            if _hpa_expresso(term, tecido_hpa):
                confirmados.append(term)
            else:
                sem_dado.append(term)
            if (i+1) % 20 == 0:
                log(f"   HPA: {i+1}/80...")
            time.sleep(0.08)
        ranked = confirmados + sem_dado + ranked[80:]
        log(f"   ✅ {len(confirmados)} com expressão confirmada")

    # ── Fase 6: Virgindade no órgão alvo ─────────────────────────────────────
    log(f"🔬 Verificando virgindade em '{alvo}'...")
    virgens, parciais, conhecidos = [], [], []
    for term in ranked[:100]:
        n = _pubmed_count(f'("{term}") AND ({alvo})', email)
        time.sleep(0.25)
        if n <= 5:       virgens.append(term)
        elif n <= 25:    parciais.append(term)
        else:            conhecidos.append(term)

    log(f"   💎 {len(virgens)} virgens | 🌱 {len(parciais)} emergentes | 🔴 {len(conhecidos)} conhecidos")

    # Monta lista priorizando virgens
    final = virgens + parciais + conhecidos
    final = final[:100]

    # ── Fase 7: Flash faz a ponte mecanismo → órgão ───────────────────────────
    if usar_ia and api_key and final:
        log(f"🤖 Flash avaliando plausibilidade para '{alvo}'...")
        # Processa em lotes de 80
        for i in range(0, min(len(final), 160), 80):
            lote = final[i:i+80]
            mapa = flash_bridge(lote, alvo, api_key)
            typemap.update(mapa)

        # Remove apenas Junk confirmado
        final = [t for t in final if typemap.get(t, "?") != "Junk"]
        log(f"   ✅ {len(final)} candidatos após avaliação de plausibilidade")

    log(f"✅ Mineração concluída: {len(final)} candidatos prontos para análise")
    return {"terms": final[:100], "typemap": typemap}

# ─── Análise Estatística — Gerador (progresso em tempo real) ─────────────────

def analisar_lista_stream(
    terms: List[str], alvo: str, contexto: str,
    email: str, y_ini: int = 2015, y_fim: int = 2026
) -> Generator[Dict, None, None]:
    N_PUBMED = 36_000_000
    n_total  = _pubmed_count(f"({alvo})", email) or 1

    for term in terms:
        q_base     = (f'("{term}") AND ({contexto}) AND ({y_ini}:{y_fim}[Date])'
                      if contexto else f'("{term}") AND ({y_ini}:{y_fim}[Date])')
        q_specific = f'("{term}") AND ({alvo}) AND ({y_ini}:{y_fim}[Date])'
        q_recent   = f'("{term}") AND ({alvo}) AND (2023:2026[Date])'

        n_base     = _pubmed_count(q_base,     email); time.sleep(0.28)
        n_specific = _pubmed_count(q_specific, email); time.sleep(0.28)
        n_recent   = _pubmed_count(q_recent,   email); time.sleep(0.28)

        enrich = (n_specific + 0.1) / max(0.1, (max(0.1, n_base) * n_total) / N_PUBMED)
        a=n_specific; b=max(0,n_base-n_specific)
        cc=max(0,n_total-n_specific); d=N_PUBMED-a-b-cc
        p_val = min(1.0,(b*cc)/max(1,a*d)) if (a*d>b*cc) else 0.999

        if n_specific<=5:    tag="blue_ocean" if n_base>20 else "ghost"
        elif n_specific<=25: tag="embryonic"
        elif enrich>5:       tag="gold"
        elif enrich>1.5:     tag="trending"
        else:                tag="saturated"

        yield {
            "term":term, "tag":tag,
            "enrichment":round(enrich,1), "pValue":round(p_val,4),
            "hitsTarget":n_specific, "hitsGlobal":n_base,
            "nRecent":n_recent, "trendRatio":round(n_recent/max(n_base,1),2),
            "_score":{"blue_ocean":1000,"embryonic":500,"gold":300,
                      "trending":200,"ghost":0,"saturated":10}.get(tag,50)
        }

# ─── Busca artigos para painel lateral ───────────────────────────────────────

def buscar_artigos(termo: str, alvo: str, email: str) -> List[Dict]:
    ids = _pubmed_search(f'("{termo}") AND ({alvo}) AND (2015:2026[Date])', email, 8)
    return _pubmed_fetch_abstracts(ids, email)

# ─── Flash analisa artigo individual (lazy) ───────────────────────────────────

def analisar_artigo_lazy(titulo: str, abstract: str, alvo: str, api_key: str) -> str:
    if not api_key: return "⚠️ Chave API necessária."
    prompt = (
        f'Bench pharmacologist. Research interest: "{alvo}"\n'
        f'TITLE: {titulo}\nABSTRACT: {abstract[:600]}\n'
        f'One line: [Tissue/Cell] → [Target/Drug/Mechanism] → [Functional Effect]\n'
        f'ONE LINE ONLY:'
    )
    res = _gemini_call(api_key, c.MODELOS["flash"], prompt, 150)
    return res.strip() if res else "Não foi possível analisar."

# ─── Pro: investigação profunda (lazy) ───────────────────────────────────────

def investigar_alvo_profundo(termo: str, alvo: str, articles: List[Dict], api_key: str) -> Dict:
    if not api_key: return {}
    textos = "\n\n".join([
        f"[{i+1}] TITLE: {a.get('title','')}\nABSTRACT: {a.get('abstract','')[:400]}"
        for i,a in enumerate(articles[:4])
    ])
    prompt = (
        f'PhD pharmacologist. Target/Drug: {termo}. Organ/Disease: {alvo}.\n'
        f'ARTICLES:\n{textos}\n\n'
        f'Return JSON:\n'
        f'{{"mechanism":"[Cell]→[Target]→[Action]→[Effect]",'
        f'"noveltyScore":<1-10>,'
        f'"druggability":"High|Medium|Low — reason",'
        f'"biologicalRationale":"why this makes sense in {alvo} even if never tested there",'
        f'"keyFindings":["f1","f2"],'
        f'"suggestedExperiments":["e1","e2"],'
        f'"redFlags":["r1 or null"],'
        f'"overallAssessment":"2 sentences"}}\n'
        f'JSON ONLY:'
    )
    res = _gemini_call(api_key, c.MODELOS["pro"], prompt, 1200)
    if not res:
        res = _gemini_call(api_key, c.MODELOS["flash"], prompt, 800)
    try:
        clean = re.sub(r"```json|```","",res).strip()
        s,e = clean.find("{"),clean.rfind("}")
        if s!=-1 and e!=-1: return json.loads(clean[s:e+1])
    except: pass
    return {}

# ─── ClinicalTrials painel lateral ───────────────────────────────────────────

def buscar_trials(termo: str, condicao: str) -> List[Dict]:
    try:
        r = requests.get(
            f"https://clinicaltrials.gov/api/v2/studies?query.term={termo}+{condicao}&pageSize=5&format=json",
            headers={"Accept":"application/json"}, timeout=10)
        return [{"nctId": s["protocolSection"]["identificationModule"].get("nctId",""),
                 "title": s["protocolSection"]["identificationModule"].get("briefTitle",""),
                 "phase": s["protocolSection"].get("designModule",{}).get("phases",["N/A"])[0],
                 "status":s["protocolSection"]["statusModule"].get("overallStatus","")}
                for s in r.json().get("studies",[])]
    except: return []
