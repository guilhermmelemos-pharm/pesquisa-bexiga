"""Lemos Lambda v2.1"""

import streamlit as st
import pandas as pd
import plotly.express as px
import constantes as c
import backend as bk

st.set_page_config(page_title="λ Lemos Lambda", page_icon="λ", layout="wide",
                   initial_sidebar_state="collapsed")

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Orbitron:wght@700;900&family=Space+Mono:wght@400;700&family=Syne:wght@600;800&display=swap');

html,body,[class*="css"]{font-family:'Syne',sans-serif;}
.stApp{background:#080c10;color:#e2e8f0;}
h1,h2,h3{font-family:'Syne',sans-serif!important;font-weight:800!important;}

.logo-main{font-family:'Orbitron',monospace;font-weight:900;font-size:1.9rem;
    color:#00ff88;letter-spacing:0.1em;text-shadow:0 0 28px rgba(0,255,136,0.45);line-height:1;}
.logo-sub{font-family:'Orbitron',monospace;font-weight:700;font-size:0.6rem;
    color:#64748b;letter-spacing:0.25em;text-transform:uppercase;margin-top:3px;}

.stButton>button{border-radius:10px;font-weight:700;height:44px;
    font-family:'Syne',sans-serif;letter-spacing:0.03em;transition:all 0.15s;}
.stButton>button[kind="primary"]{
    background:linear-gradient(135deg,#00ff88,#00cc6a)!important;
    color:#020a05!important;border:none!important;
    box-shadow:0 4px 16px rgba(0,255,136,0.25)!important;}
.stButton>button[kind="primary"]:hover{
    box-shadow:0 6px 24px rgba(0,255,136,0.4)!important;transform:translateY(-1px);}
div[data-testid="stMetricValue"]{font-family:'Space Mono',monospace!important;font-size:1.5rem!important;}

.mechanism-box{background:rgba(0,255,136,0.05);border:1px solid rgba(0,255,136,0.25);
    border-radius:10px;padding:12px 16px;font-family:'Space Mono',monospace;
    font-size:12px;color:#00ff88;margin:10px 0;}
.analysis-item{background:#161b22;border-left:3px solid #00ff88;
    padding:8px 12px;border-radius:0 8px 8px 0;font-size:12px;margin-bottom:5px;color:#94a3b8;}
.analysis-item.danger{border-color:#ff4444;}
.article-card{background:#0d1117;border:1px solid #1f2937;border-radius:10px;padding:14px;margin-bottom:10px;overflow:hidden;}
.article-title{font-size:13px;font-weight:700;color:#e2e8f0;line-height:1.4;margin-bottom:8px;word-wrap:break-word;overflow-wrap:break-word;white-space:normal;}
.ai-result{background:rgba(0,255,136,0.05);border:1px solid rgba(0,255,136,0.2);
    border-radius:8px;padding:10px 14px;font-family:'Space Mono',monospace;
    font-size:11px;color:#00ff88;margin-top:8px;}
.trial-row{background:#161b22;border-radius:8px;padding:8px 12px;margin-bottom:5px;
    font-size:11px;color:#94a3b8;}
.prog-line{font-family:'Space Mono',monospace;font-size:11px;color:#00ff88;}
.key-ok{color:#00ff88;font-size:12px;font-family:'Space Mono',monospace;}
.key-no{color:#64748b;font-size:12px;font-family:'Space Mono',monospace;}
</style>
""", unsafe_allow_html=True)

# ─── Estado ───────────────────────────────────────────────────────────────────
DEFAULTS = {
    "pagina":"home","alvos":[],"typemap":{},"resultados":None,
    "alvo_guardado":"","email_guardado":"",
    "input_alvo":"","input_contexto":"","input_email":"",
    "sistema_sel":"Selecionar Sistema...","ia_ativa":True,
    "filtro_ativo":None,"_api_key_temp":"",
    "painel_mol":None,"painel_artigos":[],"painel_analise":{},
    "ai_results":{},"show_abstract":{},
}
for k,v in DEFAULTS.items():
    if k not in st.session_state: st.session_state[k]=v

def get_api_key():
    try: return st.secrets["GEMINI_API_KEY"]
    except: return st.session_state.get("_api_key_temp","")

def has_key(): return bool(get_api_key())

def cb_sistema():
    sel=st.session_state.sistema_sel
    if sel in c.SISTEMAS: st.session_state.input_contexto=c.SISTEMAS[sel]

def adicionar_termos(lista):
    up={t.upper() for t in st.session_state.alvos}
    novos=[t for t in lista if t and t.upper() not in up]
    st.session_state.alvos.extend(novos)
    return len(novos)

def tipo_label(tipo): return c.TIPOS_LABEL.get(tipo,f"🔬 {tipo}")

# ─── Header ───────────────────────────────────────────────────────────────────
c1,c2=st.columns([3,1])
with c1:
    st.markdown('<div class="logo-main">λ LEMOS LAMBDA</div>'
                '<div class="logo-sub">Deep Science Prospector · v2.1</div>',
                unsafe_allow_html=True)
    st.write("")
with c2:
    if has_key():
        st.markdown('<p class="key-ok">🔒 Gemini 2.5 ativo</p>',unsafe_allow_html=True)
        st.caption(f"`{c.MODELOS['flash']}` flash · `{c.MODELOS['pro']}` pro")
    else:
        st.markdown('<p class="key-no">🔑 Sem chave</p>',unsafe_allow_html=True)
st.divider()

# ─── Painel lateral de artigos ────────────────────────────────────────────────
def renderizar_painel():
    mol=st.session_state.painel_mol
    if not mol: return
    st.markdown(f"### 🔬 `{mol}`")
    st.caption(f"Contexto: {st.session_state.alvo_guardado}")
    tipo_raw=st.session_state.typemap.get(mol,"Molecule")
    st.markdown(f"**Tipo:** {tipo_label(tipo_raw)}")

    if st.button("🧠 Investigar Alvo (Pro)",type="primary",use_container_width=True,key="btn_inv"):
        if not has_key(): st.warning("⚠️ Chave API necessária.")
        else:
            with st.spinner("Gemini Pro analisando..."):
                analise=bk.investigar_alvo_profundo(
                    mol,st.session_state.alvo_guardado,
                    st.session_state.painel_artigos,get_api_key())
            st.session_state.painel_analise=analise; st.rerun()

    an=st.session_state.painel_analise
    if an:
        if an.get("mechanism"):
            st.markdown(f'<div class="mechanism-box">⚡ {an["mechanism"]}</div>',unsafe_allow_html=True)
        cc1,cc2=st.columns(2)
        drugg=an.get("druggability","")
        color="#00ff88" if "High" in drugg else "#ff8c00" if "Medium" in drugg else "#94a3b8"
        cc1.markdown(f"**Druggability:** <span style='color:{color}'>{drugg.split('—')[0].strip()}</span>",
                     unsafe_allow_html=True)
        cc2.metric("Novelty",f"{an.get('noveltyScore','—')}/10")
        if an.get("biologicalRationale"): st.info(f"**Por que aqui:** {an['biologicalRationale']}")
        if an.get("keyFindings"):
            st.markdown("**Achados:**")
            for f in an["keyFindings"]:
                st.markdown(f'<div class="analysis-item">• {f}</div>',unsafe_allow_html=True)
        if an.get("suggestedExperiments"):
            with st.expander("🧪 Experimentos"):
                for e in an["suggestedExperiments"]:
                    st.markdown(f'<div class="analysis-item">🧪 {e}</div>',unsafe_allow_html=True)
        flags=[f for f in an.get("redFlags",[]) if f and f!="null"]
        if flags:
            with st.expander("⚠️ Red Flags"):
                for f in flags:
                    st.markdown(f'<div class="analysis-item danger">⚠ {f}</div>',unsafe_allow_html=True)
        if an.get("overallAssessment"): st.info(an["overallAssessment"])
        st.divider()

    trials=bk.buscar_trials(mol,st.session_state.alvo_guardado)
    if trials:
        st.markdown(f"**ClinicalTrials ({len(trials)}):**")
        for t in trials:
            st.markdown(f'<div class="trial-row"><b>{t["phase"]}</b> · {t["title"]}'
                        f' <span style="color:#64748b">· {t["status"]}</span></div>',
                        unsafe_allow_html=True)
        st.divider()

    artigos=st.session_state.painel_artigos
    if not artigos: st.info("Nenhum artigo encontrado."); return
    st.markdown(f"**Artigos ({len(artigos)}):**")
    for a in artigos:
        pmid=a.get("pmid",""); titulo=a.get("title",""); abstract=a.get("abstract","")
        st.markdown('<div class="article-card">',unsafe_allow_html=True)
        st.markdown(f'<div class="article-title">📄 {titulo}</div>',unsafe_allow_html=True)
        ca,cb=st.columns(2)
        with ca:
            if st.button("📖 Abstract",key=f"abs_{pmid}",use_container_width=True):
                st.session_state.show_abstract[pmid]=not st.session_state.show_abstract.get(pmid,False)
                st.rerun()
        with cb:
            ai_done=pmid in st.session_state.ai_results
            if st.button("✅ Ver IA" if ai_done else "🤖 Investigar com IA",
                         key=f"ai_{pmid}",use_container_width=True):
                if not has_key(): st.warning("⚠️ Chave necessária.")
                elif not ai_done:
                    with st.spinner("Flash lendo..."):
                        r=bk.analisar_artigo_lazy(titulo,abstract,st.session_state.alvo_guardado,get_api_key())
                    st.session_state.ai_results[pmid]=r; st.rerun()
        if st.session_state.show_abstract.get(pmid,False):
            st.markdown(f"> {abstract[:600]}..." if abstract else "> Sem abstract.")
        if pmid in st.session_state.ai_results:
            st.markdown(f'<div class="ai-result">⚡ {st.session_state.ai_results[pmid]}</div>',
                        unsafe_allow_html=True)
        st.markdown(f"[→ PubMed]({a.get('link','')})")
        st.markdown('</div>',unsafe_allow_html=True)

# ═══════════════════════════════════════════════════════════════════════════════
# RESULTADOS
# ═══════════════════════════════════════════════════════════════════════════════
if st.session_state.pagina=="resultados":
    if st.button("⬅ Voltar"):
        st.session_state.pagina="home"; st.session_state.painel_mol=None
        st.session_state.filtro_ativo=None; st.rerun()

    if st.session_state.painel_mol:
        col_t,col_p=st.columns([3,2])
    else:
        col_t=st.container(); col_p=None

    with col_t:
        st.markdown("## 📊 Relatório de Inteligência")
        df=st.session_state.resultados
        if df is not None and not df.empty:
            top=df.iloc[0]
            m1,m2,m3,m4=st.columns(4)
            m1.metric("🏆 Top",top["term"],c.TAGS.get(top["tag"],{}).get("label",""))
            m2.metric("λ Score",f"{top['enrichment']}×")
            m3.metric("💎 Blue Oceans",len(df[df["tag"]=="blue_ocean"]))
            m4.metric("🌱 Embrionários",len(df[df["tag"]=="embryonic"]))

            df_chart=df.head(20).copy()
            df_chart["classe"]=df_chart["tag"].map(lambda t:c.TAGS.get(t,{}).get("label",t))
            fig=px.bar(df_chart,x="term",y="enrichment",color="classe",
                color_discrete_map={v["label"]:v["color"] for v in c.TAGS.values()},
                hover_data=["hitsTarget","hitsGlobal","tipo"],
                labels={"term":"Alvo","enrichment":"λ Score","classe":"Classificação"},height=280)
            fig.update_layout(plot_bgcolor="#0d1117",paper_bgcolor="#0d1117",
                font_color="#94a3b8",xaxis_tickangle=-35,margin=dict(t=10,b=50))
            fig.update_traces(marker_line_width=0)
            st.plotly_chart(fig,use_container_width=True)

            fcols=st.columns(len(c.TAGS)+1)
            if fcols[0].button("Todos",key="f_all"):
                st.session_state.filtro_ativo=None; st.rerun()
            for i,(tk,tv) in enumerate(c.TAGS.items()):
                cnt=len(df[df["tag"]==tk])
                if cnt and fcols[i+1].button(f"{tv['label']} ({cnt})",key=f"f_{tk}"):
                    st.session_state.filtro_ativo=None if st.session_state.filtro_ativo==tk else tk
                    st.rerun()

            df_view=df if not st.session_state.filtro_ativo else df[df["tag"]==st.session_state.filtro_ativo]
            st.info(f"💡 {len(df_view)} alvos · Clique em uma linha para ver artigos e investigar")

            event=st.dataframe(
                df_view.drop(columns=["_score"],errors="ignore"),
                use_container_width=True,hide_index=True,
                on_select="rerun",selection_mode="single-row",
                column_config={
                    "term":       st.column_config.TextColumn("🔍 Alvo",width="medium"),
                    "tipo":       st.column_config.TextColumn("Tipo",width="medium"),
                    "tag":        st.column_config.TextColumn("Classificação"),
                    "enrichment": st.column_config.NumberColumn("λ Score",format="%.1f×"),
                    "pValue":     st.column_config.NumberColumn("P-Value",format="%.4f"),
                    "hitsTarget": st.column_config.NumberColumn("Hits Alvo"),
                    "hitsGlobal": st.column_config.NumberColumn("Hits Global"),
                    "nRecent":    st.column_config.NumberColumn("Rec.23-26"),
                    "trendRatio": st.column_config.NumberColumn("Trend",format="%.2f"),
                }
            )

            if event.selection.rows:
                mol=df_view.iloc[event.selection.rows[0]]["term"]
                if mol!=st.session_state.painel_mol:
                    st.session_state.painel_mol=mol
                    st.session_state.painel_analise={}
                    st.session_state.show_abstract={}
                    with st.spinner(f"Buscando artigos de {mol}..."):
                        st.session_state.painel_artigos=bk.buscar_artigos(
                            mol,st.session_state.alvo_guardado,st.session_state.email_guardado)
                    st.rerun()

            st.download_button("📥 Exportar CSV",
                df.to_csv(index=False).encode("utf-8"),
                file_name=f"lemos_{st.session_state.alvo_guardado.replace(' ','_')}.csv",
                mime="text/csv")

    if st.session_state.painel_mol and col_p:
        with col_p:
            st.divider(); renderizar_painel()

# ═══════════════════════════════════════════════════════════════════════════════
# HOME
# ═══════════════════════════════════════════════════════════════════════════════
else:
    cm,cs=st.columns([2,1])
    with cm:
        t1,t2,t3=st.tabs(["🧠 Auto-Miner","📥 Importar","💎 Presets"])

        with t1:
            st.text_input("E-mail NCBI (obrigatório)",key="input_email",placeholder="seu@email.edu")
            st.text_input("Alvo Principal (doença / órgão)",key="input_alvo",
                          placeholder="ex: Overactive Bladder, Heart Failure",
                          help=(
                              "Suporta sintaxe PubMed completa:\n\n"
                              "• Simples: `Bladder`\n"
                              "• Múltiplos órgãos: `Bladder OR Urethra OR Kidney`\n"
                              "• Excluir contexto: `Bladder NOT Cancer`\n"
                              "• Combinado: `(Bladder OR Urethra) NOT (Cancer OR Tumor)`\n"
                              "• Restringir: `Bladder AND inflammation`\n\n"
                              "Dica: use NOT Cancer se quiser focar em fisiologia normal."
                          ))
            cc1,cc2=st.columns(2)
            with cc1:
                st.selectbox("Sistema",list(c.SISTEMAS.keys()),key="sistema_sel",on_change=cb_sistema)
            with cc2:
                st.text_input("Contexto / Tecido",key="input_contexto",
                              placeholder="ex: Urothelium, Detrusor")

            if st.button("🧠 AUTO-DETECTAR ALVOS",type="primary",use_container_width=True):
                if not st.session_state.input_alvo:
                    st.error("⚠️ Defina o alvo principal.")
                elif not st.session_state.input_email:
                    st.error("⚠️ E-mail NCBI é obrigatório.")
                else:
                    logs=[]
                    def status_cb(msg):
                        logs.append(msg)

                    with st.status("🔬 Minerando candidatos...",expanded=True) as status_box:
                        def status_cb(msg):
                            st.write(msg)
                        res=bk.minerar_alvos(
                            st.session_state.input_alvo,
                            st.session_state.input_email,
                            get_api_key(),
                            st.session_state.ia_ativa,
                            status_cb=status_cb
                        )
                        status_box.update(label="✅ Mineração concluída!",state="complete")

                    qtd=adicionar_termos(res["terms"])
                    st.session_state.typemap.update(res.get("typemap",{}))
                    st.success(f"✅ {qtd} novos candidatos · Total: {len(st.session_state.alvos)}")

        with t2:
            up=st.file_uploader("Upload .txt ou .csv",type=["txt","csv"])
            if up:
                try:
                    content=up.read().decode("utf-8")
                    items=[x.strip() for line in content.splitlines()
                           for x in line.split(",") if x.strip()]
                    qtd=adicionar_termos(items)
                    st.success(f"✅ {qtd} termos importados")
                except: st.error("Erro ao ler arquivo.")

        with t3:
            for cat,items in c.PRESETS.items():
                with st.expander(cat):
                    cols=st.columns(4)
                    for i,item in enumerate(items):
                        with cols[i%4]:
                            sel=item in st.session_state.alvos
                            if st.button(f"{'✓ ' if sel else ''}{item}",
                                         key=f"p_{cat}_{item}",use_container_width=True):
                                if sel: st.session_state.alvos.remove(item)
                                else: adicionar_termos([item])
                                st.rerun()
            if st.button("💎 ADICIONAR TODOS",use_container_width=True):
                qtd=adicionar_termos([t for items in c.PRESETS.values() for t in items])
                st.success(f"✅ {qtd} presets adicionados")

        if st.session_state.alvos:
            st.divider()
            with st.expander(f"👀 Lista ({len(st.session_state.alvos)} alvos)",expanded=False):
                nova=st.text_area("Editar:",value=", ".join(st.session_state.alvos),height=100)
                if st.button("Salvar"): 
                    st.session_state.alvos=[x.strip() for x in nova.split(",") if x.strip()]
                    st.rerun()

            ce,cl=st.columns([3,1])
            with ce:
                if st.button("🚀 EXECUTAR ANÁLISE ESTATÍSTICA",type="primary",use_container_width=True):
                    if not st.session_state.input_email:
                        st.error("⚠️ E-mail NCBI obrigatório.")
                    else:
                        st.session_state.alvo_guardado=st.session_state.input_alvo
                        st.session_state.email_guardado=st.session_state.input_email
                        st.session_state.painel_mol=None; st.session_state.ai_results={}

                        termos=st.session_state.alvos; n=len(termos)
                        st.markdown("### 🧬 Analisando alvos...")
                        prog=st.progress(0); ptxt=st.empty(); tph=st.empty()
                        parciais=[]

                        for i,res in enumerate(bk.analisar_lista_stream(
                            termos,st.session_state.input_alvo,
                            st.session_state.input_contexto,st.session_state.input_email
                        )):
                            res["tipo"]=tipo_label(st.session_state.typemap.get(res["term"],"Molecule"))
                            res["tag_label"]=c.TAGS.get(res["tag"],{}).get("label",res["tag"])
                            parciais.append(res)
                            prog.progress((i+1)/n)
                            tag_info=c.TAGS.get(res["tag"],{})
                            ptxt.markdown(
                                f'<div class="prog-line">[{i+1}/{n}] {res["term"]} → '
                                f'{tag_info.get("label","?")} | λ {res["enrichment"]}× | '
                                f'Hits: {res["hitsTarget"]}</div>',unsafe_allow_html=True)
                            if (i+1)%5==0 or (i+1)==n:
                                df_p=pd.DataFrame(parciais).sort_values(
                                    by=["_score","enrichment"],ascending=[False,False])
                                tph.dataframe(
                                    df_p[["term","tipo","tag_label","enrichment","hitsTarget","hitsGlobal"]],
                                    use_container_width=True,hide_index=True)

                        df_final=pd.DataFrame(parciais).sort_values(
                            by=["_score","enrichment"],ascending=[False,False])
                        st.session_state.resultados=df_final
                        st.session_state.pagina="resultados"; st.rerun()

            with cl:
                if st.button("🗑 Limpar",use_container_width=True):
                    st.session_state.alvos=[]; st.session_state.typemap={}; st.rerun()

    with cs:
        st.markdown("### ⚙️ Config")
        try:
            st.secrets["GEMINI_API_KEY"]; st.success("🔒 Chave via Secrets")
        except:
            with st.expander("🔑 API Key Gemini",expanded=not has_key()):
                st.caption("Fica só nesta sessão.")
                kin=st.text_input("Chave:",type="password",key="key_sb",placeholder="AIza...")
                if kin: st.session_state._api_key_temp=kin; st.success("✅ Pronta")
                st.markdown("[🔑 Gerar grátis](https://aistudio.google.com/app/apikey)")

        st.toggle("🤖 IA na classificação",key="ia_ativa",value=True)
        st.divider()
        st.markdown("**Uso de IA:**")
        st.caption("⚡ Flash → classifica (1x) + artigo (lazy)\n🧠 Pro → investigação (lazy)\n✅ Análise = zero tokens")
        st.divider()
        st.markdown("**Fontes:**")
        st.caption("📚 PubMed · 🎯 OpenTargets\n🧬 HPA · 🏥 ClinicalTrials")
        st.divider()
        st.caption("Lemos Lambda v2.1")
        st.caption("Pix: `960f3f16-06ce-4e71-9b5f-6915b2a10b5a`")
