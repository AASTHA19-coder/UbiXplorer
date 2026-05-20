# =========================================================
# PATHWAY EXPLORER — FINAL CINEMATIC VERSION
# REAL BIOLOGY + BEAUTIFUL VISUALIZATION
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
# CUSTOM CSS
# =========================================================

st.markdown("""
<style>

.stApp{
    background: linear-gradient(180deg,#020617 0%, #071426 100%);
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

.metric-container{
    background:#0f172a;
    padding:18px;
    border-radius:16px;
}

</style>
""", unsafe_allow_html=True)

# =========================================================
# HEADER
# =========================================================

st.markdown(
    "<div class='title'>Pathway Explorer</div>",
    unsafe_allow_html=True
)

st.markdown(
    """
    <div class='sub'>
    AI-powered systems biology intelligence engine for ubiquitin pathway exploration
    </div>
    """,
    unsafe_allow_html=True
)

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
# GENE INPUT
# =========================================================

if (
    source == "Uploaded Genes"
    and
    "gene_list" in st.session_state
    and
    len(st.session_state["gene_list"]) > 0
):

    genes = sorted(
        list(set(st.session_state["gene_list"]))
    )

else:

    genes = sorted(
        list(adata_h.var_names)
    )

# =========================================================
# SELECT GENE
# =========================================================

gene = st.selectbox(
    "Select Gene",
    genes
)

gene = gene.upper()

# =========================================================
# STRING API
# =========================================================

@st.cache_data(show_spinner=False)
def get_string_interactions(query_gene):

    try:

        string_api_url = "https://string-db.org/api"
        output_format = "json"
        method = "network"

        request_url = "/".join([
            string_api_url,
            output_format,
            method
        ])

        params = {

            "identifiers": query_gene,
            "species": 9606,
            "required_score": 700,
            "limit": 20
        }

        response = requests.get(
            request_url,
            params=params
        )

        data = response.json()

        interactors = set()
        edges = []

        for item in data:

            p1 = item["preferredName_A"]
            p2 = item["preferredName_B"]

            interactors.add(p1)
            interactors.add(p2)

            edges.append((p1, p2))

        return list(interactors), edges

    except:

        return [query_gene], []

# =========================================================
# GET ENRICHMENT GENES
# =========================================================

# =========================================================
# GET ENRICHMENT GENES
# =========================================================

# ALWAYS get STRING interactors
# around currently selected gene

interactors, string_edges = get_string_interactions(
    gene
)

interactors = [
    g.upper()
    for g in interactors
]

# ---------------------------------------------------------
# REFERENCE MODE
# ---------------------------------------------------------

if source == "Reference Atlas":

    enrich_genes = interactors

# ---------------------------------------------------------
# UPLOADED MODE
# ---------------------------------------------------------

else:

    uploaded_genes = set([
        g.upper()
        for g in genes
    ])

    # KEEP ONLY uploaded genes
    # connected to selected gene

    enrich_genes = list(

        uploaded_genes.intersection(
            set(interactors)
        )
    )

    # IMPORTANT FALLBACK

    if len(enrich_genes) < 3:

        enrich_genes = interactors[:15]

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

    if len(string_edges) > 0:

        for edge in string_edges:

            G.add_edge(edge[0], edge[1])

    else:

        for pathway in enrich_df["Pathway"][:10]:

            G.add_edge(gene, pathway)

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
            color="rgba(148,163,184,0.20)"
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

            size=30,

            color="#22d3ee",

            line=dict(
                width=2,
                color="#8b5cf6"
            ),

            opacity=0.92
        ),

        textfont=dict(
            size=13,
            color="white"
        )
    )

    fig = go.Figure(
        data=[edge_trace, node_trace]
    )

    fig.update_layout(

        template="plotly_dark",

        height=780,

        paper_bgcolor="#020617",
        plot_bgcolor="#020617",

        showlegend=False
    )

    st.plotly_chart(
        fig,
        use_container_width=True
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

        color_continuous_scale=[

            [0.0, "#1e293b"],
            [0.3, "#22d3ee"],
            [0.6, "#8b5cf6"],
            [1.0, "#ec4899"]

        ],

        template="plotly_dark",

        height=780
    )

    fig.update_layout(

        paper_bgcolor="#020617",
        plot_bgcolor="#020617",

        font=dict(
            size=15,
            color="white"
        )
    )

    st.plotly_chart(
        fig,
        use_container_width=True
    )

    # =====================================================
    # VIOLIN PLOT
    # =====================================================

    st.subheader("Distribution of Enrichment Scores")

    violin_df = enrich_df.copy()

    if "Gene_set" not in violin_df.columns:

        violin_df["Gene_set"] = "Unknown"

    fig2 = px.violin(

        violin_df,

        x="Gene_set",
        y="EnrichmentScore",

        box=True,

        points="all",

        template="plotly_dark",

        color_discrete_sequence=["#8b5cf6"],

        height=720
    )

    fig2.update_layout(

        paper_bgcolor="#020617",
        plot_bgcolor="#020617",

        font=dict(
            size=15,
            color="white"
        )
    )

    st.plotly_chart(
        fig2,
        use_container_width=True
    )

# =========================================================
# TAB 3 — SANKEY
# =========================================================

# =========================================================
# TAB 3 — SANKEY
# =========================================================

with tab3:

    st.subheader("Gene → Process Systems Architecture")

    # -----------------------------------------------------
    # TOP ENRICHED TERMS
    # -----------------------------------------------------

    top_df = enrich_df.head(8).copy()

    # -----------------------------------------------------
    # CLEAN PATHWAY LABELS
    # -----------------------------------------------------

    top_df["CleanTerm"] = (

        top_df["Pathway"]
        .astype(str)
        .str.replace(r"\s*\(GO:\d+\)", "", regex=True)
        .str.replace(r"\s*R-HSA-\d+", "", regex=True)
        .str.slice(0, 55)

    )

    # -----------------------------------------------------
    # DETECT TRUE DISEASE TERMS
    # ONLY FROM ENRICHR OUTPUT
    # -----------------------------------------------------

    disease_df = top_df[

        top_df["CleanTerm"].str.contains(

            "disease|cancer|carcinoma|leukemia|"
            "infection|diabetes|alzheimer|"
            "parkinson|arthritis",

            case=False,
            na=False

        )

    ]

    process_df = top_df[

        ~top_df.index.isin(disease_df.index)

    ]

    # -----------------------------------------------------
    # NODE LABELS
    # -----------------------------------------------------

    selected_gene = gene.upper()

    labels = [selected_gene]

    process_labels = (
        process_df["CleanTerm"]
        .unique()
        .tolist()
    )

    disease_labels = (
        disease_df["CleanTerm"]
        .unique()
        .tolist()
    )

    labels.extend(process_labels)

    # ONLY ADD DISEASE LAYER
    # IF REAL DISEASE TERMS EXIST

    if len(disease_labels) > 0:

        labels.extend(disease_labels)

    # -----------------------------------------------------
    # NODE MAP
    # -----------------------------------------------------

    node_map = {

        label: idx
        for idx, label in enumerate(labels)

    }

    # -----------------------------------------------------
    # LINKS
    # -----------------------------------------------------

    source_nodes = []
    target_nodes = []
    values = []
    link_colors = []

    # -----------------------------------------------------
    # GENE → PROCESS
    # -----------------------------------------------------

    for _, row in process_df.iterrows():

        process_name = row["CleanTerm"]

        source_nodes.append(
            node_map[selected_gene]
        )

        target_nodes.append(
            node_map[process_name]
        )

        values.append(

            max(
                int(row["GeneCount"]),
                1
            )

        )

        link_colors.append(
            "rgba(34,211,238,0.25)"
        )

    # -----------------------------------------------------
    # PROCESS → DISEASE
    # ONLY IF REAL DISEASE TERMS EXIST
    # -----------------------------------------------------

    if len(disease_labels) > 0:

        process_cycle = process_labels.copy()

        for i, (_, row) in enumerate(
            disease_df.iterrows()
        ):

            disease_name = row["CleanTerm"]

            process_name = process_cycle[
                i % len(process_cycle)
            ]

            source_nodes.append(
                node_map[process_name]
            )

            target_nodes.append(
                node_map[disease_name]
            )

            values.append(

                max(
                    int(row["GeneCount"]),
                    1
                )

            )

            link_colors.append(
                "rgba(236,72,153,0.22)"
            )

    # -----------------------------------------------------
    # NODE COLORS
    # -----------------------------------------------------

    node_colors = []

    for label in labels:

        if label == selected_gene:

            node_colors.append("#22d3ee")

        elif label in disease_labels:

            node_colors.append("#ec4899")

        else:

            node_colors.append("#8b5cf6")

    # -----------------------------------------------------
    # SANKEY FIGURE
    # -----------------------------------------------------

    fig = go.Figure(data=[

        go.Sankey(

            arrangement="snap",

            node=dict(

                pad=28,

                thickness=28,

                line=dict(
                    color="#020617",
                    width=2
                ),

                label=labels,

                color=node_colors,

                hovertemplate=
                "%{label}<extra></extra>"

            ),

            link=dict(

                source=source_nodes,

                target=target_nodes,

                value=values,

                color=link_colors,

                hovertemplate=
                "Connection Strength: %{value}<extra></extra>"

            )

        )

    ])

    # -----------------------------------------------------
    # LAYOUT
    # -----------------------------------------------------

    fig.update_layout(

        template="plotly_dark",

        paper_bgcolor="#020617",
        plot_bgcolor="#020617",

        font=dict(

            size=16,
            color="white"

        ),

        height=850,

        margin=dict(

            l=20,
            r=20,
            t=60,
            b=20

        )

    )

    st.plotly_chart(

        fig,
        use_container_width=True

    )

# =========================================================

# =========================================================
# TAB 4 — HEATMAP


with tab4:

    st.subheader("Pathway Activation Heatmap")

    # -----------------------------------------------------
    # TOP PATHWAYS
    # -----------------------------------------------------

    heat_df = enrich_df.head(10).copy()

    # -----------------------------------------------------
    # GENES TO SHOW
    # -----------------------------------------------------

    heat_genes = enrich_genes[:10]

    # -----------------------------------------------------
    # BUILD REAL MATRIX
    # -----------------------------------------------------

    matrix = []

    for pathway_genes in heat_df["Genes"]:

        pathway_gene_list = [

            g.strip().upper()

            for g in str(pathway_genes).split(";")
        ]

        row = []

        for g in heat_genes:

            if g.upper() in pathway_gene_list:

                row.append(1)

            else:

                row.append(0)

        matrix.append(row)

    heat_matrix = np.array(matrix)

    # -----------------------------------------------------
    # HEATMAP
    # -----------------------------------------------------
# =========================================================
# TAB 4 — FUNCTIONAL CONNECTIVITY MAP
# =========================================================
# =========================================================
# TAB 4 — FUNCTIONAL OVERLAP ARCHITECTURE
# =========================================================

with tab4:

    st.subheader("🧬 Functional Overlap Architecture")

    st.caption(
        "Shared gene architecture across enriched biological pathways."
    )

    # -----------------------------------------------------
    # TOP PATHWAYS
    # -----------------------------------------------------

    upset_df = enrich_df.head(8).copy()

    # -----------------------------------------------------
    # BUILD GENE ↔ PATHWAY MATRIX
    # -----------------------------------------------------

    pathway_dict = {}

    for _, row in upset_df.iterrows():

        pathway = row["Pathway"]

        genes_in_pathway = [

            g.strip().upper()

            for g in str(row["Genes"]).split(";")
        ]

        filtered = [

            g for g in genes_in_pathway

            if g in enrich_genes
        ]

        pathway_dict[pathway] = filtered

    # -----------------------------------------------------
    # UNIQUE GENES
    # -----------------------------------------------------

    all_genes = sorted(

        list(set(

            g

            for geneset in pathway_dict.values()

            for g in geneset
        ))
    )

    # -----------------------------------------------------
    # MATRIX
    # -----------------------------------------------------

    matrix = []

    for gene_name in all_genes:

        row = []

        for pathway in pathway_dict.keys():

            if gene_name in pathway_dict[pathway]:

                row.append(1)

            else:

                row.append(0)

        matrix.append(row)

    matrix_df = pd.DataFrame(

        matrix,

        index=all_genes,

        columns=list(pathway_dict.keys())
    )

    # -----------------------------------------------------
    # INTERSECTION COUNTS
    # -----------------------------------------------------

    intersection_counts = matrix_df.sum(axis=0)

    # -----------------------------------------------------
    # BARPLOT
    # -----------------------------------------------------

    fig = go.Figure()

    fig.add_trace(

        go.Bar(

            x=intersection_counts.index,

            y=intersection_counts.values,

            marker=dict(

                color=[

                    "#38bdf8",
                    "#22d3ee",
                    "#06b6d4",
                    "#8b5cf6",
                    "#a855f7",
                    "#c084fc",
                    "#ec4899",
                    "#f472b6"

                ],

                line=dict(
                    color="#e2e8f0",
                    width=1
                )
            ),

            hovertemplate=
            "<b>%{x}</b><br>Shared Genes: %{y}<extra></extra>"
        )
    )

    # -----------------------------------------------------
    # LAYOUT
    # -----------------------------------------------------

    fig.update_layout(

        template="plotly_dark",

        height=720,

        paper_bgcolor="#020617",
        plot_bgcolor="#020617",

        title=dict(

            text="Pathway Intersection Landscape",

            font=dict(
                size=24,
                color="white"
            )
        ),

        xaxis=dict(

            tickangle=-20,

            title="Enriched Pathways",

            title_font=dict(
                size=16
            ),

            tickfont=dict(
                size=12
            )
        ),

        yaxis=dict(

            title="Number of Shared Genes",

            title_font=dict(
                size=16
            ),

            tickfont=dict(
                size=12
            )
        ),

        font=dict(
            color="white"
        ),

        margin=dict(
            l=40,
            r=40,
            t=80,
            b=140
        )
    )

    st.plotly_chart(
        fig,
        use_container_width=True
    )

    # =====================================================
    # BINARY OVERLAP DOT MATRIX
    # =====================================================

    st.subheader("🧬 Gene–Pathway Membership Matrix")

    dot_x = []
    dot_y = []
    dot_color = []
    dot_size = []

    for i, gene_name in enumerate(matrix_df.index):

        for j, pathway in enumerate(matrix_df.columns):

            val = matrix_df.loc[gene_name, pathway]

            dot_x.append(pathway)
            dot_y.append(gene_name)

            if val == 1:

                dot_color.append("#ec4899")
                dot_size.append(18)

            else:

                dot_color.append("#334155")
                dot_size.append(8)

    dot_fig = go.Figure()

    dot_fig.add_trace(

        go.Scatter(

            x=dot_x,
            y=dot_y,

            mode="markers",

            marker=dict(

                size=dot_size,

                color=dot_color,

                line=dict(
                    color="#94a3b8",
                    width=0.5
                )
            ),

            hovertemplate=
            "Gene: %{y}<br>Pathway: %{x}<extra></extra>"
        )
    )

    dot_fig.update_layout(

        template="plotly_dark",

        height=720,

        paper_bgcolor="#020617",
        plot_bgcolor="#020617",

        xaxis=dict(

            tickangle=-25,

            title="Pathways",

            title_font=dict(
                size=15
            ),

            tickfont=dict(
                size=11
            )
        ),

        yaxis=dict(

            title="Genes",

            title_font=dict(
                size=15
            ),

            tickfont=dict(
                size=11
            )
        ),

        font=dict(
            color="white"
        ),

        margin=dict(
            l=40,
            r=40,
            t=40,
            b=120
        )
    )

    st.plotly_chart(
        dot_fig,
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

        color_continuous_scale=[

            [0.0, "#164e63"],
            [0.4, "#22d3ee"],
            [0.7, "#8b5cf6"],
            [1.0, "#ec4899"]

        ],

        template="plotly_dark",

        height=900
    )

    fig.update_layout(

        paper_bgcolor="#020617",
        plot_bgcolor="#020617",

        font=dict(
            size=20,
            color="white"
        ),

        uniformtext=dict(
            minsize=16,
            mode='hide'
        )
    )

    st.plotly_chart(
        fig,
        use_container_width=True
    )

# =========================================================
# FOOTER
# =========================================================

st.markdown("""

<br><br>

<center>

AI-powered pathway intelligence engine for ubiquitin systems biology

</center>

""", unsafe_allow_html=True)
