import streamlit as st
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ======================================================
# PAGE CONFIG
# ======================================================
st.set_page_config(
    page_title="Spatial Explorer",
    #page_icon="🧠",
    layout="wide"
)

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
max-width:1350px;
padding-top:2rem;
padding-bottom:2rem;
}

.title{
font-size:42px;
font-weight:800;
color:#2dd4bf;
margin-bottom:8px;
}

.sub{
font-size:18px;
color:#cbd5e1;
margin-bottom:25px;
}

.footer{
text-align:center;
color:#94a3b8;
font-size:14px;
margin-top:50px;
}

</style>
""", unsafe_allow_html=True)

# ======================================================
# LOAD REFERENCE DATA
# ======================================================
@st.cache_resource
def load_data():

    human = ad.read_h5ad(
        "assets/human/reference/human_adata.h5ad"
    )

    mouse = ad.read_h5ad(
       "assets/mouse/reference/mouse_adata.h5ad"
    )

    human.var_names_make_unique()
    mouse.var_names_make_unique()

    return human, mouse

adata_h, adata_m = load_data()

# ======================================================
# HEADER
# ======================================================
st.markdown(
    "<div class='title'> Spatial Explorer</div>",
    unsafe_allow_html=True
)

st.markdown("""
<div class='sub'>
Disease-aware cross-species spatial transcriptomics explorer.
</div>
""", unsafe_allow_html=True)

# ======================================================
# SOURCE MODE
# ======================================================
source = st.radio(
    "Gene Source",
    ["Reference Atlas", "Uploaded Genes"],
    horizontal=True
)

# ======================================================
# GENE LIST
# ======================================================
if (
    source == "Uploaded Genes"
    and "gene_list" in st.session_state
):
    genes = sorted(st.session_state["gene_list"])

else:
    genes = sorted(list(adata_h.var_names))

# ======================================================
# GENE SELECTION
# ======================================================
gene = st.selectbox(
    "Select Gene",
    genes
)

# ======================================================
# INFO
# ======================================================
if source == "Uploaded Genes":

    st.info(
        "Uploaded DEG expression is projected onto the reference spatial atlas."
    )

# ======================================================
# HUMAN → MOUSE GENE MAP
# ======================================================
mouse_map = {
    "MYC": "Myc",
    "RELA": "Rela",
    "MYD88": "Myd88",
    "NFKBIA": "Nfkbia",
    "TLR2": "Tlr2",
    "FOXO1": "Foxo1"
}

mouse_gene = mouse_map.get(
    gene.upper(),
    gene.title()
)

# ======================================================
# GET XY COORDINATES
# ======================================================
# ======================================================
# GET XY COORDINATES
# ======================================================
def get_xy(adata):

    coords = adata.obsm["spatial"]

    x = coords[:, 0].astype(float)
    y = coords[:, 1].astype(float)

    

    return x, y
###############################################################
# =========================================================
# APPROXIMATE BRAIN REGION ANNOTATION
# =========================================================
# =========================================================
# APPROXIMATE BRAIN REGION ANNOTATION
# =========================================================

def assign_brain_region(x, y, species="human"):

    # =====================================================
    # HUMAN
    # =====================================================

    if species == "human":

        # LEFT CORTEX
        if x < 9000 and y < 9000:
            return "Prefrontal Cortex"

        # LOWER LEFT
        elif x < 9000 and y >= 9000:
            return "Temporal Cortex"

        # CENTRAL UPPER
        elif 9000 <= x < 15000 and y < 11000:
            return "Hippocampus"

        # CENTRAL LOWER
        elif 9000 <= x < 15000 and y >= 11000:
            return "Cingulate Cortex"

        # RIGHT SIDE
        elif x >= 15000:
            return "White Matter"

        else:
            return "Unknown"

    # =====================================================
    # MOUSE
    # =====================================================

    else:

        # LEFT UPPER
        if x < 3500 and y < 4500:
            return "Striatum"

        # LEFT LOWER
        elif x < 3500 and y >= 4500:
            return "Cortex"

        # CENTRAL
        elif 3500 <= x < 6000 and y < 6000:
            return "Hippocampus"

        # CENTRAL LOWER
        elif 3500 <= x < 6000 and y >= 6000:
            return "Thalamus"

        # RIGHT
        elif x >= 6000:
            return "White Matter"

        else:
            return "Unknown"
# ======================================================
def get_vals(adata, selected_gene):

    # ==================================================
    # GENE ABSENCE SAFETY
    # ==================================================
    if selected_gene not in adata.var_names:

        return np.zeros(adata.shape[0])

    # ==================================================
    # LOAD REFERENCE EXPRESSION
    # ==================================================
    vals = adata[:, selected_gene].X

    if hasattr(vals, "toarray"):

        vals = vals.toarray().flatten()

    else:

        vals = np.array(vals).flatten()

    vals = vals.astype(float)

    # ==================================================
    # REFERENCE MODE
    # ==================================================
    if source == "Reference Atlas":

        return vals

    # ==================================================
    # UPLOADED MODE
    # ==================================================
    if "user_df" not in st.session_state:

        return vals

    df = st.session_state["user_df"]

    # ==================================================
    # FIND MATCHING USER GENE
    # ==================================================
    row_df = df[
        df["Gene"].astype(str).str.upper()
        ==
        gene.upper()
    ]

    if row_df.empty:

        return vals

    # ==================================================
    # LOGFC
    # ==================================================
    fc = abs(
        float(
            row_df["logFC"].values[0]
        )
    )

    # ==================================================
    # NORMALIZATION
    # ==================================================
    vals = (
        vals - vals.mean()
    ) / (
        vals.std() + 1e-6
    )

    # ==================================================
    # DISEASE PERTURBATION
    # ==================================================
    weighted = np.power(
        np.maximum(vals, 0) + 1,
        1 + (fc / 2)
    )

    # ==================================================
    # LOCAL HOTSPOT BOOST
    # ==================================================
    hotspot_mask = weighted > np.percentile(
        weighted,
        80
    )

    weighted[hotspot_mask] *= (
        1 + (fc / 3)
    )

    # ==================================================
    # FINAL NORMALIZATION
    # ==================================================
    weighted = (
        weighted - weighted.mean()
    ) / (
        weighted.std() + 1e-6
    )

    weighted[
        weighted < 0
    ] = 0

    return weighted

# ======================================================
# CLEAN AXIS
# ======================================================
def clean_axis(ax, x, y):

    xt = np.linspace(
        np.min(x),
        np.max(x),
        5
    )

    yt = np.linspace(
        np.min(y),
        np.max(y),
        5
    )

    ax.set_xticks(xt)
    ax.set_yticks(yt)

    ax.set_xticklabels(
        [f"{int(i)}" for i in xt],
        color="white",
        fontsize=9
    )

    ax.set_yticklabels(
        [f"{int(i)}" for i in yt],
        color="white",
        fontsize=9
    )
  
# ======================================================
# SPATIAL MAP
# ======================================================
def spatial_map(ax, adata, vals, ttl):

    x, y = get_xy(adata)

    sc = ax.scatter(
        x,
        y,
        c=vals,
        cmap="magma",
        s=20,
        edgecolors="none"
    )

    ax.set_title(
        ttl,
        color="white",
        fontsize=18,
        weight="bold"
    )

    ax.set_facecolor("#06111f")

    ax.invert_yaxis()

    ax.set_xlabel(
        "Spatial X",
        color="white"
    )

    ax.set_ylabel(
        "Spatial Y",
        color="white"
    )

    clean_axis(ax, x, y)

    cbar = plt.colorbar(
        sc,
        ax=ax,
        fraction=0.046,
        pad=0.04
    )

    cbar.set_label(
        "Spatial Burden",
        color="white"
    )

    cbar.ax.yaxis.set_tick_params(
        color="white"
    )

    plt.setp(
        cbar.ax.get_yticklabels(),
        color="white"
    )

# ======================================================
# HOTSPOT MAP
# ======================================================
# ======================================================
# HOTSPOT MAP
# ======================================================
def hotspot_map(ax, adata, vals, ttl, species="human"):

    x, y = get_xy(adata)

    z = (
        vals - vals.mean()
    ) / (
        vals.std() + 1e-6
    )

    hot = z > 2

    # ==================================================
    # BACKGROUND
    # ==================================================

       # ==================================================
    # REGION-BASED BACKGROUND COLORING
    # ==================================================

  

    # ==================================================
    # HOTSPOTS
    # ==================================================

    ax.scatter(
        x[hot],
        y[hot],
        c="#ef4444",
        s=30,
        edgecolors="white",
        linewidth=0.2
    )

    # ==================================================
    # TITLES
    # ==================================================

    ax.set_title(
        ttl,
        color="white",
        fontsize=18,
        weight="bold"
    )

    ax.set_facecolor("#06111f")

    ax.invert_yaxis()

    ax.set_xlabel(
        "Spatial X",
        color="white"
    )

    ax.set_ylabel(
        "Spatial Y",
        color="white"
    )

    clean_axis(ax, x, y)
        # ==================================================
    # SUBTLE ANATOMICAL OVERLAYS
    # ==================================================

    overlays = []

    if species == "human":

        overlays = [

            ("Prefrontal Cortex", 7000, 5000, "#8b5cf6"),
            ("Temporal Cortex", 7000, 15000, "#2563eb"),
            ("Hippocampus", 11500, 8000, "#14b8a6"),
            ("Cingulate Cortex", 11500, 15000, "#f59e0b"),
            ("White Matter", 18000, 11000, "#94a3b8")
        ]

    else:

        overlays = [

            ("Striatum", 2600, 3000, "#8b5cf6"),
            ("Cortex", 2600, 7600, "#2563eb"),
            ("Hippocampus", 4700, 4300, "#14b8a6"),
            ("Thalamus", 4700, 7600, "#ef4444"),
            ("White Matter", 7200, 5000, "#94a3b8")
        ]

    # --------------------------------------------------
    # DRAW SOFT OVERLAYS
    # --------------------------------------------------

    for region, lx, ly, clr in overlays:

        ax.scatter(

            [lx],
            [ly],

            s=22000,

            c=clr,

            alpha=0.06,

            edgecolors="none"
        )

        ax.text(

            lx,
            ly,

            region,

            fontsize=8,

            color="white",

            ha="center",

            alpha=0.75,

            weight="bold"
        )
# ==================================================
# REGION LABEL OVERLAY
# ==================================================
    # ==================================================
    # REGION LABEL OVERLAY
    # ==================================================

    if species == "human":

        label_positions = {

            "Prefrontal Cortex": (6500, 5000),
            "Temporal Cortex": (6500, 15000),
            "Hippocampus": (11500, 7000),
            "Cingulate Cortex": (11500, 15000),
            "White Matter": (18000, 11000)
        }

    else:

        label_positions = {

            "Striatum": (2500, 3000),
            "Cortex": (2500, 7500),
            "Hippocampus": (4500, 4500),
            "Thalamus": (4500, 8000),
            "White Matter": (7200, 5000)
        }

    # --------------------------------------------------
    # DRAW LABELS
    # --------------------------------------------------

    for region, (lx, ly) in label_positions.items():

        ax.text(

            lx,
            ly,

            region,

            fontsize=9,

            color="white",

            weight="bold",

            ha="center",

            alpha=0.9,

            bbox=dict(

                facecolor="black",

                alpha=0.35,

                edgecolor="none",

                pad=2
            )
        )
    # ==================================================
    # REGION INTERPRETATION
    # ==================================================

    hot_x = x[hot]
    hot_y = y[hot]

    region_labels = []

    for xi, yi in zip(hot_x, hot_y):

        region = assign_brain_region(
            xi,
            yi,
            species=species
        )

        region_labels.append(region)

    region_df = pd.DataFrame({

        "Region": region_labels
    })

    region_counts = (

        region_df["Region"]
        .value_counts()
        .reset_index()
    )

    region_counts.columns = [

        "Brain Region",
        "Hotspot Count"
    ]

    return (
        int(hot.sum()),
        round(float(z.max()), 2),
        region_counts
    )
# ======================================================
# TABS
# ======================================================
tab1, tab2, tab3 = st.tabs([
    "Spatial Maps",
    "Hotspots",
    "Summary"
])

# ======================================================
# TAB 1
# ======================================================
with tab1:

    col1, col2 = st.columns(2)

    # HUMAN
    with col1:

        vals_h = get_vals(
            adata_h,
            gene
        )

        fig, ax = plt.subplots(
            figsize=(7, 7),
            facecolor="#06111f"
        )

        spatial_map(
            ax,
            adata_h,
            vals_h,
            "Human Disease Projection"
            if source == "Uploaded Genes"
            else f"Human - {gene}"
        )

        st.pyplot(fig)

    # MOUSE
    with col2:

        vals_m = get_vals(
            adata_m,
            mouse_gene
        )

        fig, ax = plt.subplots(
            figsize=(7, 7),
            facecolor="#06111f"
        )

        spatial_map(
            ax,
            adata_m,
            vals_m,
            "Mouse Disease Projection"
            if source == "Uploaded Genes"
            else f"Mouse - {mouse_gene}"
        )

        st.pyplot(fig)

# ======================================================
# TAB 2
# ======================================================
with tab2:

    col1, col2 = st.columns(2)

    # HUMAN
    with col1:

        vals_h = get_vals(
            adata_h,
            gene
        )

        fig, ax = plt.subplots(
            figsize=(7, 7),
            facecolor="#06111f"
        )

        human_hot, hz, human_regions = hotspot_map(
            ax,
            adata_h,
            vals_h,
            "Human Disease Hotspots"
            if source == "Uploaded Genes"
            else f"Human Hotspots - {gene}", 
            species="human"
        )

        st.pyplot(fig)

    # MOUSE
    with col2:

        vals_m = get_vals(
            adata_m,
            mouse_gene
        )

        fig, ax = plt.subplots(
            figsize=(7, 7),
            facecolor="#06111f"
        )

        mouse_hot, mz, mouse_regions = hotspot_map(
            ax,
            adata_m,
            vals_m,
            "Mouse Disease Hotspots"
            if source == "Uploaded Genes"
            else f"Mouse Hotspots - {mouse_gene}", 
            species="mouse"
        )

        st.pyplot(fig)

    # SUMMARY
    st.markdown(
        "### Hotspot Detection Summary"
    )

    st.info(f"""
Spatial vulnerability is computed by projecting uploaded disease-associated DEGs onto a reference spatial transcriptomics atlas.

Human hotspots: {human_hot}

Mouse hotspots: {mouse_hot}

Human Max Z-score: {hz}

Mouse Max Z-score: {mz}
""")

# ======================================================
# REGIONAL INTERPRETATION
# ======================================================

st.markdown(
    "### Human Regional Hotspots"
)

st.dataframe(
    human_regions,
    use_container_width=True
)

st.markdown(
    "### Mouse Regional Hotspots"
)

st.dataframe(
    mouse_regions,
    use_container_width=True
)

# ======================================================
# TAB 3
# ======================================================
with tab3:

    a, b, c, d = st.columns(4)

    a.metric(
        "Mode",
        source
    )

    b.metric(
        "Human Hotspots",
        human_hot
    )

    c.metric(
        "Mouse Hotspots",
        mouse_hot
    )

    d.metric(
        "Max Z-score",
        max(hz, mz)
    )

    st.markdown(
        "### Cross-Species Interpretation"
    )

    if human_hot > mouse_hot:

        st.info(
            "Human tissue shows stronger projected vulnerability."
        )

    elif mouse_hot > human_hot:

        st.info(
            "Mouse tissue shows stronger projected vulnerability."
        )

    else:

        st.success(
            "Comparable projected burden across species."
        )


# ======================================================
# FOOTER
# ======================================================
st.markdown("""
<div class='footer'>
UbiXplorer: A Proof-Of-Concept Prototype created @ Molecular Neuroscience and Functional Genomics Lab |
Department of Biotechnology | Delhi Technological University
</div>
""", unsafe_allow_html=True)
