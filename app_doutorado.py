import streamlit as st

# Configuração da página para parecer profissional
st.set_page_config(page_title="Lemos Lambda v2.1", page_icon="λ", layout="centered")

# CSS para o botão brilhar e centralizar tudo
st.markdown("""
<style>
.stApp {
    background-color: #0E1117;
    color: white;
    text-align: center;
    padding-top: 50px;
}
h1 {
    font-size: 3rem !important;
    font-weight: 800;
}
.big-btn {
    display: inline-block;
    background-color: #FF4B4B;
    color: white !important;
    padding: 20px 40px;
    border-radius: 12px;
    font-size: 24px;
    font-weight: bold;
    text-decoration: none;
    margin-top: 40px;
    box-shadow: 0 4px 20px rgba(255, 75, 75, 0.5);
    transition: transform 0.2s, box-shadow 0.2s;
}
.big-btn:hover {
    transform: scale(1.05);
    box-shadow: 0 6px 30px rgba(255, 75, 75, 0.7);
}
.info-box {
    background-color: #262730;
    padding: 20px;
    border-radius: 10px;
    border: 1px solid #444;
    margin-top: 30px;
    display: inline-block;
}
</style>
""", unsafe_allow_html=True)

# Conteúdo Principal
st.title("λ Lemos Lambda evoluiu.")
st.header("A versão 2.0 (Python) foi descontinuada.")

st.markdown("""
<div class="info-box">
    <h3>🚀 Migramos para v2.1 (React Engine)</h3>
    <p>A nova versão é <strong>10x mais rápida</strong>, possui interface moderna e inteligência artificial aprimorada.</p>
</div>
""", unsafe_allow_html=True)

# --- COLOQUE SEU LINK DA VERCEL AQUI EMBAIXO ---
LINK_NOVO = "https://seu-projeto-na-vercel.app" 

st.markdown(f'<a href="{LINK_NOVO}" class="big-btn" target="_blank">ACESSAR VERSÃO 2.1 AGORA</a>', unsafe_allow_html=True)

st.write("###")
st.write("###")
st.caption("Se você chegou aqui via DOI/Zenodo, o link acima é o oficial e atualizado.")
