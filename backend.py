def analisar_abstract_com_ia(titulo, abstract, api_key, lang='pt'):
    if not api_key: return "⚠️ IA não ativada"
    try:
        genai.configure(api_key=api_key)
        idioma = "Português" if lang == 'pt' else "Inglês"
        
        # Se o abstract for nulo ou cortado, focamos na Inferência Sênior
        abstract_txt = abstract if (abstract and len(abstract) > 20) else "Abstract incompleto. Infira pelo título."
        
        prompt = f"""Como PhD em Farmacologia, analise:
        TITULO: {titulo}
        RESUMO: {abstract_txt[:3000]}
        
        AÇÃO: Use o TITULO para inferir o mecanismo farmacológico caso o resumo esteja ausente.
        FORMATO OBRIGATÓRIO: Alvo → Fármaco/Substância → Efeito funcional na bexiga.
        REGRAS: Máximo 30 palavras. Idioma: {idioma}."""

        # --- MIGRAÇÃO PARA NOMES ESTÁVEIS (Evita 404) ---
        # Tentamos o Gemini 1.5 Flash estável (sem o prefixo models/ que às vezes buga)
        # Ou o novo Gemini 2.0 Flash se disponível na sua região/SDK
        try:
            model = genai.GenerativeModel('gemini-1.5-flash')
            response = model.generate_content(prompt)
        except:
            # Fallback para a versão estável específica
            model = genai.GenerativeModel('gemini-1.5-pro')
            response = model.generate_content(prompt)
            
        return response.text.strip()
    
    except Exception as e:
        erro_str = str(e)
        if "429" in erro_str: 
            return f"💡 Inferência: {titulo} → (Cota excedida. Tente em 1 min)"
        return f"❌ Erro na IA: {erro_str[:40]}"
