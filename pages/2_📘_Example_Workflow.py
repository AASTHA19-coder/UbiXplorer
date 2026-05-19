import streamlit as st
import pandas as pd
import math



st.markdown("""
<style>
.footer{
text-align:center;
color:#94a3b8;
font-size:14px;
margin-top:45px;
}
</style>
""", unsafe_allow_html=True)

st.set_page_config(layout="wide")

st.title("Example Workflow")

# =========================================================
# INTRO
# =========================================================
st.info(
"""
This example demonstrates **UbiXplorer** using our study dataset comparing
**Alzheimer’s disease vs Non-demented controls**. Differentially expressed ubiquitin-linked hub genes were prioritized and
projected onto HEALTHY Human and Mouse brain spatial 10X Visium Tissue Datasets.
"""
)

# =========================================================
# WORKFLOW
# =========================================================
st.markdown("""
### Workflow

**Upload DEG Data** → **Ub Gene Detection** → **Hub Prioritization** →  **Human / Mouse Spatial Maps** → **Hotspot Detection** → **STRING / KEGG**
""")

st.divider()

# =========================================================
# REAL DATA
# =========================================================
data = [
["TLR2",   0.4615600, 3.550780e-56, 27, 73],
["NFKBIA", 0.3603699, 3.325875e-72, 392,108],
["MYD88",  0.2470767, 6.546462e-47, 74,308],
["RELA",   0.2228550, 2.958587e-54, 640,85],
["FOXO1",  0.2459081, 7.728418e-60, 360,107],
["MYC",    0.2045997, 1.297806e-29, 56,170],
]

df = pd.DataFrame(
    data,
    columns=[
        "Gene","logFC","adj.P.Val",
        "Human Hotspots","Mouse Hotspots"
    ]
)

# Priority score
df["Priority Score"] = (
    abs(df["logFC"])*10 +
    (-df["adj.P.Val"].apply(lambda x: math.log10(x))) +
    ((df["Human Hotspots"]+df["Mouse Hotspots"])/50)
)

df = df.sort_values("Priority Score", ascending=False)

# =========================================================
# TABLES
# =========================================================
c1, c2 = st.columns(2)

with c1:
    st.subheader("**Study DEG Input**")
    st.dataframe(
        df[["Gene","logFC","adj.P.Val"]],
        use_container_width=True
    )

with c2:
    st.subheader("**Prioritized Ub Hub Genes**")
    st.dataframe(
        df[[
            "Gene",
            "Priority Score",
            "Human Hotspots",
            "Mouse Hotspots"
        ]],
        use_container_width=True
    )

st.caption(
"""
Priority Score integrates differential expression significance, fold-change magnitude, and cross-species hotspot evidence.
"""
)

st.divider()

# =========================================================
# GENE SELECTOR
# =========================================================
gene = st.selectbox(
    "**Select Gene**",
    df["Gene"].tolist()
)

mouse_gene = {
    "MYC":"Myc",
    "RELA":"Rela",
    "MYD88":"Myd88",
    "NFKBIA":"Nfkbia",
    "TLR2":"Tlr2",
    "FOXO1":"Foxo1"
}[gene]

row = df[df["Gene"] == gene].iloc[0]

# =========================================================
# TABS
# =========================================================
tab1, tab2, tab3, tab4 = st.tabs(
[
"**Spatial Maps**",
"**Hotspots**",
"**STRING / KEGG**",
"**Summary**"
]
)

# =========================================================
# MAPS
# =========================================================
with tab1:

    a,b = st.columns(2)

    with a:
        st.subheader("Human Brain")
        st.image(
            f"assets/human/maps/{gene}.png",
            use_container_width=True
        )

    with b:
        st.subheader("Mouse Brain")
        st.image(
            f"assets/mouse/maps/{mouse_gene}.png",
            use_container_width=True
        )

# =========================================================
# HOTSPOTS
# =========================================================
with tab2:

    a,b = st.columns(2)

    with a:
        st.subheader(
            f"Human Hotspots ({int(row['Human Hotspots'])})"
        )
        st.image(
            f"assets/human/hotspots/{gene}_hotspot.png",
            use_container_width=True
        )

    with b:
        st.subheader(
            f"Mouse Hotspots ({int(row['Mouse Hotspots'])})"
        )
        st.image(
            f"assets/mouse/hotspots/{mouse_gene}_hotspot.png",
            use_container_width=True
        )

# =========================================================
# PATHWAYS
# =========================================================
with tab3:

    st.subheader(f"{gene} Network & Pathway Links")

    a,b = st.columns(2)

    with a:
        st.markdown(
            f"""
            **STRING Database**

            https://string-db.org/cgi/network?identifier={gene}
            """
        )

    with b:
        st.markdown(
            f"""
            **KEGG Search**

            https://www.genome.jp/dbget-bin/www_bget?{gene}
            """
        )

# =========================================================
# SUMMARY
# =========================================================
with tab4:

    c1,c2,c3,c4 = st.columns(4)

    c1.metric("Gene", gene)
    c2.metric("logFC", round(row["logFC"],3))
    c3.metric(
        "Human Hotspots",
        int(row["Human Hotspots"])
    )
    c4.metric(
        "Mouse Hotspots",
        int(row["Mouse Hotspots"])
    )

    diff = abs(
        row["Human Hotspots"] -
        row["Mouse Hotspots"]
    )

    if diff < 80:
        msg = "Highly conserved hotspot architecture across species."
    elif row["Human Hotspots"] > row["Mouse Hotspots"]:
        msg = "Human-enriched hotspot burden detected."
    else:
        msg = "Mouse-enriched hotspot burden detected."

    st.success(msg)

    st.info(
        f"{gene} shows statistically significant differential "
        f"expression in AD vs Control with spatial localization "
        f"support across atlases."
    )
    # ======================================================
# FOOTER
# ======================================================
st.markdown("""
<div class='footer'>
UbiXplorer: A Proof-Of-Concept Prototype created @ Molecular Neuroscience and Functional Genomics Lab   |     
                Department of Biotechnology |   Delhi Technological University
</div>
""", unsafe_allow_html=True)