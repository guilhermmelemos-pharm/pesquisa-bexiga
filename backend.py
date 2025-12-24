def analisar_abstract_com_ia(titulo, abstract, api_key, lang="pt"):
    if not api_key:
        return "⚠️ IA desativada"

    try:
        genai.configure(api_key=api_key)

        idioma = "Português" if lang == "pt" else "Inglês"

        texto_base = (
            abstract if abstract and len(abstract) > 40
            else "Resumo indisponível. Inferir mecanismo com base no título."
        )

        prompt = f"""
Você é um PhD em Farmacologia Molecular e Biologia Translacional.

OBJETIVO:
Extrair APENAS entidades mecanísticas e farmacológicas do artigo.

ACEITO (priorizar):
- Alvos moleculares
- Receptores
- Canais iônicos
- Genes
- Proteínas
- Fármacos, agonistas, antagonistas, moduladores
- Vias de sinalização

PROIBIDO (NÃO mencionar):
- Diagnósticos clínicos (ex: HBP, OAB, LUTS, IC)
- Doenças
- Termos médicos sindrômicos
- Desfechos clínicos genéricos
- Epidemiologia

ARTIGO:
TÍTULO: {titulo}

RESUMO:
{texto_base[:3000]}

FORMATO OBRIGATÓRIO (se houver evidência clara):
Alvo/Gene/Receptor → Fármaco/Substância → Efeito molecular ou funcional na bexiga

REGRAS:
- Não inferir se não houver evidência
- Não repetir alvos clássicos sem base no texto
- Se não houver entidade molecular clara, responder:
  "Sem alvo molecular ou farmacológico identificável"
- Máximo 30 palavras
- Idioma: {idioma}
"""

        model = genai.GenerativeModel("models/gemini-2.5-flash")
        response = model.generate_content(prompt)

        return response.text.strip()

    except Exception:
        return "⚠️ IA indisponível no momento"
