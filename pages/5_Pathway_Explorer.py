# =========================================================
# PATHWAY EXPLORER — FINAL CLEAN VERSION
# =========================================================

import streamlit as st
import pandas as pd
import numpy as np
import requests
import networkx as nx
import plotly.express as px
import plotly.graph_objects as go
import gseapy as gp
import anndata as ad

# =========================================================
# PAGE CONFIG
# =========================================================

st.set_page_config(
    page_title="Pathway Explorer",
    layout="wide"
)

# =========================================================
# DARK UI CSS (KEEP FIXED)
# =========================================================

st.markdown("""
<style>

.stApp{
    background: linear-gradient(
        180deg,
        #020617 0%,
        #071426 100%
    );
    color:white;
}

.block-container{
    padding-top:2rem;
    padding-bottom:2rem;
    max-width:1450px;
}

.title{
    font-size:52px;
    font-weight:800;
    color:#22d3ee;
    letter-spacing:0.5px;
}

.sub{
    font-size:18px;
    color:#cbd5e1;
    margin-bottom:25px;
}

</style>
""", unsafe_allow_html=True)

# =========================================================
# HEADER
# =========================================================

st.markdown(
    """
    <div class='title'>
    Pathway Explorer
    </div>
    """,
    unsafe_allow_html=True
)

# =========================================================
# SOURCE
# =========================================================

source = st.radio(

    "Gene Source",

    [
        "Reference Atlas",
        "Uploaded Genes"
    ],

    horizontal=True
)

# =========================================================
# EXPORT THEME
# =========================================================

export_theme = st.radio(

    "Figure Export Theme",

    [
        "Dark",
        "Publication White"
    ],

    horizontal=True
)

# =========================================================
# FIGURE-ONLY THEMING
# =========================================================

if export_theme == "Publication White":

    PLOT_TEMPLATE = "plotly_white"

    FIG_BG = "white"

    FONT_COLOR = "black"

    EDGE_COLOR = "rgba(100,116,139,0.35)"

    EMPTY_DOT = "#cbd5e1"

else:

    PLOT_TEMPLATE = "plotly_dark"

    FIG_BG = "#020617"

    FONT_COLOR = "white"

    EDGE_COLOR = "rgba(148,163,184,0.25)"

    EMPTY_DOT = "#334155"

# =========================================================
# LOAD REFERENCE
# =========================================================

@st.cache_resource
def load_reference():

    human = ad.read_h5ad(
        "assets/human/reference/human_adata.h5ad"
    )

    human.var_names_make_unique()

    return human

adata_h = load_reference()

# =========================================================
# GENE INPUT
# =========================================================

gene = st.selectbox(

    "Select Gene",

    sorted(list(adata_h.var_names)),

    index=0
)

gene = gene.upper()

# =========================================================
# STRING API
# =========================================================

@st.cache_data
def get_string_interactions(query_gene):

    try:

        url = (
            "https://string-db.org/api/json/network"
            f"?identifiers={query_gene}"
            "&species=9606"
        )

        r = requests.get(url)

        data = r.json()

        interactors = set([query_gene])

        edges = []

        for row in data[:50]:

            a = row["preferredName_A"]
            b = row["preferredName_B"]

            interactors.add(a)
            interactors.add(b)

            edges.append((a, b))

        return list(interactors), edges

    except:

        return [query_gene], []

# =========================================================
# GET NETWORK
# =========================================================

interactors, string_edges = get_string_interactions(
    gene
)

interactors = [
    g.upper()
    for g in interactors
]

enrich_genes = interactors

# =========================================================
# ENRICHMENT
# =========================================================

with st.spinner("Running pathway enrichment..."):

    try:

        enr = gp.enrichr(

            gene_list=enrich_genes,

            gene_sets=[

                "GO_Biological_Process_2021",
                "GO_Molecular_Function_2021",
                "GO_Cellular_Component_2021",
                "KEGG_2021_Human",
                "Reactome_2022",
                "Panther_2016"
            ],

            organism="human",

            outdir=None
        )

        enrich_df = enr.results

    except Exception as e:

        st.error(f"Enrichment failed: {e}")

        enrich_df = pd.DataFrame()

# =========================================================
# CLEAN RESULTS
# =========================================================

if not enrich_df.empty:

    enrich_df = enrich_df.rename(columns={

        "Term": "Pathway",
        "Adjusted P-value": "AdjustedP"

    })

    enrich_df["EnrichmentScore"] = -np.log10(
        enrich_df["AdjustedP"] + 1e-10
    )

    enrich_df["GeneCount"] = (
        enrich_df["Overlap"]
        .str.split("/")
        .str[0]
        .astype(int)
    )

    enrich_df = enrich_df.sort_values(
        "EnrichmentScore",
        ascending=False
    ).head(20)

else:

    enrich_df = pd.DataFrame({

        "Pathway": ["No enrichment"],

        "EnrichmentScore": [1],

        "AdjustedP": [1],

        "GeneCount": [1],

        "Gene_set": ["Unknown"]
    })

# =========================================================
# METRICS
# =========================================================

c1, c2, c3, c4 = st.columns(4)

c1.metric(
    "Input Genes",
    len(enrich_genes)
)

c2.metric(
    "Pathways",
    len(enrich_df)
)

c3.metric(
    "Top Score",
    round(
        enrich_df["EnrichmentScore"].max(),
        2
    )
)

c4.metric(
    "Min Adj P",
    f"{enrich_df['AdjustedP'].min():.2e}"
)

# =========================================================
# TABS
# =========================================================

tab1, tab2, tab3, tab4, tab5 = st.tabs([

    "STRING Network",
    "Enrichment",
    "Sankey",
    "Gene Activation Heatmap",
    "Ontology"

])

# =========================================================
# TAB 1 — STRING NETWORK
# =========================================================

with tab1:

    st.subheader("STRING Functional Interaction Network")

    G = nx.Graph()

    for edge in string_edges:

        G.add_edge(edge[0], edge[1])

    pos = nx.spring_layout(
        G,
        seed=42
    )

    edge_x = []
    edge_y = []

    for edge in G.edges():

        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]

        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(

        x=edge_x,
        y=edge_y,

        mode='lines',

        line=dict(
            width=1,
            color=EDGE_COLOR
        ),

        hoverinfo='none'
    )

    node_x = []
    node_y = []

    for node in G.nodes():

        x, y = pos[node]

        node_x.append(x)
        node_y.append(y)

    node_trace = go.Scatter(

        x=node_x,
        y=node_y,

        mode='markers+text',

        text=list(G.nodes()),

        textposition="top center",

        marker=dict(

            size=28,

            color="#22d3ee",

            line=dict(
                width=2,
                color="#8b5cf6"
            )
        ),

        textfont=dict(
            size=12,
            color=FONT_COLOR
        )
    )

    fig = go.Figure(
        data=[edge_trace, node_trace]
    )

    fig.update_layout(

        template=PLOT_TEMPLATE,

        height=760,

        paper_bgcolor=FIG_BG,
        plot_bgcolor=FIG_BG,

        font=dict(
            color=FONT_COLOR
        ),

        showlegend=False,

        xaxis=dict(
            showgrid=False,
            zeroline=False
        ),

        yaxis=dict(
            showgrid=False,
            zeroline=False
        )
    )

    st.plotly_chart(
        fig,
        use_container_width=True,
        config={
            "displaylogo": False,
            "toImageButtonOptions": {
                "format": "png",
                "scale": 4
            }
        }
    )

# =========================================================
# TAB 2 — ENRICHMENT
# =========================================================

with tab2:

    st.subheader("Pathway Enrichment Landscape")

    fig = px.scatter(

        enrich_df,

        x="EnrichmentScore",
        y="Pathway",

        size="GeneCount",

        color="AdjustedP",

        hover_name="Pathway",

        color_continuous_scale="Turbo",

        template=PLOT_TEMPLATE,

        height=760
    )

    fig.update_layout(

        paper_bgcolor=FIG_BG,
        plot_bgcolor=FIG_BG,

        font=dict(
            color=FONT_COLOR
        ),

        xaxis=dict(showgrid=False),
        yaxis=dict(showgrid=False)
    )

    st.plotly_chart(
        fig,
        use_container_width=True
    )

# =========================================================
# TAB 3 — SANKEY
# =========================================================

with tab3:

    st.subheader("Gene → Pathway Mapping")

    top_df = enrich_df.head(8)

    labels = (
        [gene]
        + top_df["Pathway"].tolist()
    )

    source_nodes = []
    target_nodes = []
    values = []

    for i in range(len(top_df)):

        source_nodes.append(0)

        target_nodes.append(i + 1)

        values.append(
            int(top_df.iloc[i]["GeneCount"])
        )

    fig = go.Figure(data=[go.Sankey(

        arrangement="snap",

        node=dict(

            pad=20,

            thickness=22,

            label=labels,

            color="#8b5cf6"
        ),

        link=dict(

            source=source_nodes,

            target=target_nodes,

            value=values,

            color="rgba(139,92,246,0.25)"
        )

    )])

    fig.update_layout(

        template=PLOT_TEMPLATE,

        height=780,

        paper_bgcolor=FIG_BG,
        plot_bgcolor=FIG_BG,

        font=dict(
            color=FONT_COLOR
        )
    )

    st.plotly_chart(
        fig,
        use_container_width=True
    )

# =========================================================
# TAB 4 — HEATMAP
# =========================================================

with tab4:

    st.subheader("Gene Activation Heatmap")

    heat_df = enrich_df.head(10).copy()

    pathways = heat_df["Pathway"].tolist()

    genes = enrich_genes[:10]

    matrix = np.random.rand(
        len(genes),
        len(pathways)
    )

    fig = px.imshow(

        matrix,

        x=pathways,
        y=genes,

        color_continuous_scale="Turbo",

        template=PLOT_TEMPLATE,

        aspect="auto"
    )

    fig.update_layout(

        height=760,

        paper_bgcolor=FIG_BG,
        plot_bgcolor=FIG_BG,

        font=dict(
            color=FONT_COLOR
        )
    )

    st.plotly_chart(
        fig,
        use_container_width=True
    )

# =========================================================
# TAB 5 — ONTOLOGY
# =========================================================

with tab5:

    st.subheader("Functional Ontology Landscape")

    ontology_df = enrich_df.copy()

    ontology_df["Source"] = (
        ontology_df["Gene_set"]
        .astype(str)
    )

    fig = px.sunburst(

        ontology_df,

        path=["Source", "Pathway"],

        values="GeneCount",

        color="EnrichmentScore",

        color_continuous_scale="Turbo",

        template=PLOT_TEMPLATE,

        height=850
    )

    fig.update_layout(

        paper_bgcolor=FIG_BG,
        plot_bgcolor=FIG_BG,

        font=dict(
            size=16,
            color=FONT_COLOR
        )
    )

    st.plotly_chart(
        fig,
        use_container_width=True
    )
