import streamlit as st
import pandas as pd
import numpy as np

st.set_page_config(page_title="Upload Data", layout="wide")

# ======================================================
# CSS
# ======================================================
st.markdown("""
<style>

.stApp{
background: linear-gradient(180deg,#030712 0%, #06111f 100%);
color:white;
}

.block-container{
max-width:1200px;
padding-top:2rem;
padding-bottom:2rem;
}

.title{
font-size:42px;
font-weight:800;
color:#2dd4bf;
}

.sub{
font-size:18px;
color:#cbd5e1;
margin-bottom:20px;
}

.card{
background:rgba(255,255,255,0.03);
padding:25px;
border-radius:22px;
border:1px solid rgba(45,212,191,0.15);
}

.footer{
text-align:center;
color:#94a3b8;
font-size:14px;
margin-top:40px;
}

</style>
""", unsafe_allow_html=True)

# ======================================================
# HEADER
# ======================================================
st.markdown("<div class='title'>Upload DEG Dataset</div>", unsafe_allow_html=True)
st.markdown(
"<div class='sub'>Upload differential expression results in the same format as Example Workflow.</div>",
unsafe_allow_html=True
)

# ======================================================
# TEMPLATE
# ======================================================
with st.expander("**Required Input Format**", expanded=True):

    temp = pd.DataFrame({
        "Gene":["MYC","RELA","MYD88"],
        "logFC":[0.20,0.22,0.24],
        "adj.P.Val":[1e-20,2e-12,4e-10]
    })

    st.dataframe(temp, use_container_width=True)

# ======================================================
# UPLOAD
# ======================================================
file = st.file_uploader(
    "Upload CSV or XLSX file",
    type=["csv","xlsx"]
)

# ======================================================
# COLUMN DETECTION
# ======================================================
def normalize(x):
    return x.lower().replace(" ","").replace("_","").replace(".","")

gene_opts = ["gene","symbol","genesymbol","hgnc"]
fc_opts   = ["logfc","log2fc","log2foldchange","fc","foldchange"]
p_opts    = ["adjpval","padj","fdr","qvalue","adjp"]

if file:

    try:
        df = pd.read_csv(file) if file.name.endswith(".csv") else pd.read_excel(file)

        cols = {normalize(c):c for c in df.columns}

        def find_col(options):
            for o in options:
                if o in cols:
                    return cols[o]
            return None

        gcol = find_col(gene_opts)
        fcol = find_col(fc_opts)
        pcol = find_col(p_opts)

        if gcol is None:
            st.error("No gene column detected.")
            st.stop()

        out = pd.DataFrame()
        out["Gene"] = df[gcol].astype(str).str.upper().str.strip()
        out["logFC"] = pd.to_numeric(df[fcol], errors="coerce") if fcol else 0
        out["adj.P.Val"] = pd.to_numeric(df[pcol], errors="coerce") if pcol else 1

        out = out.dropna(subset=["Gene"]).drop_duplicates("Gene")

        # store globally
        st.session_state["user_df"] = out
        st.session_state["gene_list"] = out["Gene"].tolist()

        st.success(f"Loaded {len(out)} genes successfully.")

        c1,c2,c3 = st.columns(3)

        c1.metric("Genes", len(out))
        c2.metric("Mean logFC", round(out["logFC"].mean(),3))
        c3.metric("Significant Hits", int((out["adj.P.Val"] < 0.05).sum()))

        st.markdown("### Preview")
        st.dataframe(out.head(20), use_container_width=True)

        csv = out.to_csv(index=False).encode()
        st.download_button(
            "⬇ Download Standardized File",
            csv,
            file_name="ubixplorer_cleaned.csv",
            mime="text/csv"
        )

        st.info("Next Step: Open Spatial Explorer from sidebar.")

    except Exception as e:
        st.error(f"Upload failed: {e}")

# ======================================================
# FOOTER
# ======================================================
st.markdown("""
<div class='footer'>
UbiXplorer: A Proof-Of-Concept Prototype created @ Molecular Neuroscience and Functional Genomics Lab   |     
                Department of Biotechnology |   Delhi Technological University
</div>
""", unsafe_allow_html=True)
