def _faxina_ia(lista_suja):
    api_key = st.session_state.get('api_key_usuario', '').strip()
    if not api_key: 
        return [t for t in lista_suja if t.upper() not in BLACKLIST_RADICAL][:30]

    # PROMPT DE "FILTRO PASSIVO": Proibido adicionar termos novos
    prompt = f"""
    ACT AS: A strict laboratory filter for a PhD Researcher.
    
    INPUT LIST: [{", ".join(lista_suja)}]
    
    TASK: 
    1. Look at the INPUT LIST provided.
    2. DELETE any term that is clinical, surgical, anatomical, or vague filler.
    3. KEEP only specific molecular targets, signaling proteins, or drugs present IN THE ORIGINAL LIST.
    
    CRITICAL RULE: 
    - DO NOT add any new terms. 
    - DO NOT suggest "trending" drugs if they are not in the input. 
    - If the input list is "Bladder, With, TRPV4", your output MUST be ["TRPV4"].
    
    OUTPUT: Return strictly a Python list of strings containing ONLY items from the input.
    """
    
    payload = {
        "contents": [{"parts": [{"text": prompt}]}],
        "generationConfig": {"temperature": 0.0} 
    }
    
    for m in MODELOS_ATIVOS:
        try:
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{m}:generateContent?key={api_key}"
            resp = requests.post(url, json=payload, timeout=15)
            
            if resp.status_code == 200:
                raw_text = resp.json()['candidates'][0]['content']['parts'][0]['text']
                clean_text = re.sub(r'```[a-zA-Z]*', '', raw_text).replace('```', '').strip()
                res = ast.literal_eval(clean_text)
                # Verifica se a IA tentou "vender" termos novos
                if isinstance(res, list):
                    return [t for t in res if t in lista_suja]
        except: continue
    
    return [t for t in lista_suja if t.upper() not in BLACKLIST_RADICAL][:30]
