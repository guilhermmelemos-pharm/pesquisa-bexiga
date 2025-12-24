# --- 5. MINERAÇÃO E CRUZAMENTO (VERSÃO ELITE MOLECULAR) ---
@st.cache_data(ttl=3600, show_spinner=False)
def buscar_alvos_emergentes_pubmed(termo_base, email, usar_ia=True):
    if email: Entrez.email = email
    
    # Brainstorming: A IA projeta o que é "ponta" antes de ver o lixo do PubMed
    alvos_previstos = []
    if usar_ia:
        prompt_brain = f"""
        As a Senior PhD in Pharmacology, list 60 specific molecular targets (channels, receptors, enzymes) and new drugs for {termo_base}. 
        STRICT EXCLUSION: No general physiology (Sodium, Water, Toad, Frog, Muscle, Nerve, Cell).
        FOCUS: Bench research (TRPV4, SPHK1, GSK1016790A, P2X3, NLRP3, Piezo1).
        Output: ONLY a Python list [].
        """
        alvos_previstos = _chamar_ia_simples(prompt_brain)
    
    # Busca PubMed focada apenas em 2015-2026 (mata os artigos de sapo de 1970)
    query = f"({termo_base} AND (Pharmacology OR Molecular)) AND (2015:2026[Date])"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1000, sort="relevance")
        record = Entrez.read(handle); handle.close()
        total_docs = int(record["Count"])
        
        handle = Entrez.efetch(db="pubmed", id=record["IdList"], rettype="medline", retmode="text")
        full_data = handle.read(); handle.close()
        
        # --- BLACKLIST DE "TERMOS BURROS" (GENERALISTAS) ---
        blacklist_ruido = {
            'STUDIES', 'ACTIVE', 'POSTERIOR', 'PITUITARY', 'ION', 'OSMOSIS', 'MUSCLE', 'NERVE', 
            'CELL', 'STIMULATION', 'ACID', 'COMPOUNDS', 'NERVOUS', 'MUCOUS', 'CALCIUM', 'ISO',
            'ISOTOPES', 'SYSTEM', 'ISOLATED', 'EFFECT', 'ACTION', 'MEMBRANE', 'METABOLISM',
            'BIOLOGICAL', 'PERMEABILITY', 'THE', 'AND', 'FOR', 'DNA', 'RNA', 'FUNCTION'
        }
        # Blacklist de animais clássicos
        blacklist_animais = {'TOAD', 'TOADS', 'FROG', 'FROGS', 'DOGS', 'CATS', 'RABBITS'}
        
        candidatos_pubmed = []
        for artigo in full_data.split("\n\nPMID-"):
            texto_util = ""
            for linha in artigo.split("\n"):
                # SÓ OLHA TÍTULO E KEYWORDS
                if linha.startswith("TI  - ") or linha.startswith("KW  - ") or linha.startswith("OT  - "):
                    texto_util += linha[6:].strip() + " "
            
            # REGEX NINJA: Captura siglas (ex: TRPV4) ou Compostos (ex: GSK101)
            # Ignora palavras comuns minúsculas
            encontrados = re.findall(r'\b[A-Z0-9-]{3,}\b', texto_util)
            for t in encontrados:
                t_up = t.upper()
                if t_up not in blacklist_ruido and t_up not in blacklist_animais and not t_up.isdigit():
                    candidatos_pubmed.append(t_up)

        contagem = Counter(candidatos_pubmed)
        top_pubmed = [termo for termo, freq in contagem.most_common(180)]
        
        # Passo 3: Cruzamento (A IA atua como o filtro final "chato")
        nomes_finais = []
        if usar_ia:
            lista_cruzamento = list(set(alvos_previstos + top_pubmed))
            prompt_cross = f"PhD Review: From this list {lista_cruzamento[:120]}, keep ONLY specific pharmacological targets and drugs. DELETE all general biology and system terms. Output: ONLY a Python list []."
            nomes_finais = _chamar_ia_simples(prompt_cross)
        
        if not nomes_finais: nomes_finais = top_pubmed[:60]

        res_finais = []
        for nome in nomes_finais:
            n_limpo = limpar_termo_para_pubmed(nome)
            if n_limpo.upper() in blacklist_ruido: continue
            freq = contagem.get(n_limpo.upper(), 1)
            l_score, p_val, b_ocean, status = calcular_metricas_originais(freq, total_docs, freq * 10)
            res_finais.append({"Alvo": n_limpo, "Lambda": round(l_score, 2), "P-value": round(p_val, 4), "Blue Ocean": round(b_ocean, 1), "Status": status})
        return res_finais
    except: return []
