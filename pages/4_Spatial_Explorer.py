import streamlit as st
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ======================================================
# PAGE CONFIG
# ======================================================f
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
# EXPORT THEME
# ======================================================

export_theme = st.radio(

    "Export Theme",

    ["Dark", "Publication White"],

    horizontal=True
)
# ======================================================
# THEME COLORS
# ======================================================

if export_theme == "Publication White":

    BG = "white"
    FONT = "black"
    GRID = "#d1d5db"

else:

    BG = "#06111f"
    FONT = "white"
    GRID = "#334155"
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
        color="FONT",
        fontsize=18,
        weight="bold"
    )

    ax.set_facecolor(BG)

    ax.invert_yaxis()

    ax.set_xlabel(
        "Spatial X",
        color="FONT"
    )

    ax.set_ylabel(
        "Spatial Y",
        color="FONT"
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
        color="FONT"
    )

    cbar.ax.yaxis.set_tick_params(
        color="FONT"
    )

    plt.setp(
        cbar.ax.get_yticklabels(),
        color="FONT"
    )

# ======================================================
# HOTSPOT MAP
# ======================================================
# ======================================================
# HOTSPOT MAP
# ======================================================
# ======================================================
# HOTSPOT MAP
# ======================================================
def hotspot_map(ax, adata, vals, ttl, species="human"):

    import seaborn as sns
    from pathlib import Path

    # ==================================================
    # COORDINATES
    # ==================================================

    x, y = get_xy(adata)

    z = (
        vals - vals.mean()
    ) / (
        vals.std() + 1e-6
    )

    hot = z > 2

    BASE_DIR = Path(__file__).resolve().parent.parent

    # ==================================================
    # LOAD TISSUE IMAGE
    # ==================================================

    tissue_img = None

    try:

        if species == "human":

            tissue_path = (

                BASE_DIR /
                "assets" /
                "human" /
                "reference" /
                "tissue_hires_image.png"
            )

        else:

            tissue_path = (

                BASE_DIR /
                "assets" /
                "mouse" /
                "reference" /
                "tissue_hires_image.png"
            )

        tissue_img = plt.imread(str(tissue_path))

    except Exception as e:

        st.warning(f"Tissue image not loaded: {e}")

    # ==================================================
    # BACKGROUND
    # ==================================================

    ax.set_facecolor(BG)

    # ==================================================
    # DRAW TISSUE IMAGE
    # ==================================================

    if tissue_img is not None:

        ax.imshow(

            tissue_img,

            extent=[
                x.min(),
                x.max(),
                y.max(),
                y.min()
            ],

            alpha=0.03,

            aspect="auto"
        )

    # ==================================================
    # TISSUE SPOTS
    # ==================================================

    ax.scatter(

        x,
        y,

        c="#94a3b8",

        s=12,

        alpha=0.03,

        edgecolors="none"
    )

    # ==================================================
    # KDE HOTSPOT MAP
    # ==================================================

    sns.kdeplot(

    x=x[hot],
    y=y[hot],

    fill=True,

    cmap="inferno",

    alpha=0.50,

    levels=45,

    thresh=0.12,

    bw_adjust=0.38,

    linewidths=0,

    ax=ax
)

    # ==================================================
    # HOTSPOT CORE POINTS
    # ==================================================

    ax.scatter(

        x[hot],
        y[hot],

        c="#ff8fa3",

        s=10,

        alpha=0.9,

        edgecolors="black" if export_theme == "Publication White" else "white",

        linewidth=0.2,

        zorder=5
    )

    # ==================================================
    # TITLE
    # ==================================================

    ax.set_title(

        ttl,

        color="white",

        fontsize=20,

        weight="bold"
    )

# ==================================================
# CLEAN MINIMAL AXES
# ==================================================
    # ==================================================
    # CLEAN MINIMAL AXES
    # ==================================================

    ax.set_xlabel(

        "Spatial X",

        color="white",

        fontsize=10,

        alpha=0.75
    )

    ax.set_ylabel(

        "Spatial Y",

        color="white",

        fontsize=10,

        alpha=0.75
    )

    xt = np.linspace(

        np.min(x),
        np.max(x),
        4
    )

    yt = np.linspace(

        np.min(y),
        np.max(y),
        4
    )

    ax.set_xticks(xt)
    ax.set_yticks(yt)

    ax.tick_params(

        colors="white",

        labelsize=8,

        length=0
    )

    for spine in ax.spines.values():

        spine.set_visible(False)

    ax.invert_yaxis()

    ax.set_aspect("auto")

    # ==================================================
    # REGION LABELS
    # ==================================================

    if species == "human":

        labels = [

            ("Prefrontal Cortex", 7000, 5000),
            ("Temporal Cortex", 7000, 15000),
            ("Hippocampus", 11500, 7000),
            ("Cingulate Cortex", 11500, 15000),
            ("White Matter", 18000, 11000)
        ]

    else:

        labels = [

            ("Striatum", 2600, 3500),
            ("Cortex", 2600, 7600),
            ("Hippocampus", 4700, 5000),
            ("Thalamus", 4700, 8200),
            ("White Matter", 7200, 5200)
        ]

    # ==================================================
    # DRAW LABELS
    # ==================================================

        # ==================================================
    # DRAW LABELS
    # ==================================================

    for txt, lx, ly in labels:

        ax.text(

            lx,
            ly,

            txt,

            fontsize=10,

            color="FONT",

            weight="bold",

            ha="center",

            va="center",

            alpha=0.92,

            zorder=10,

            bbox=dict(

                facecolor="white" if export_theme == "Publication White" else "black",

                alpha=0.18,

                edgecolor="none",

                boxstyle="round,pad=0.15"
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
            facecolor=BG
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

        fig, ax = plt.subplots(F
            figsize=(7, 7),
            facecolor=BG
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
    # ==================================================
    # REGIONAL HOTSPOT BARPLOTS
    # ==================================================

    import matplotlib.pyplot as plt
    import seaborn as sns

    human_plot = human_regions.sort_values(

        by="Hotspot Count",
        ascending=True
    )

    mouse_plot = mouse_regions.sort_values(

        by="Hotspot Count",
        ascending=True
    )

    fig, axes = plt.subplots(

        1,
        2,

        figsize=(14, 5),

        facecolor=BG
    )

    # ----------------------------------------------
    # HUMAN
    # ----------------------------------------------

    sns.barplot(

        data=human_plot,

        x="Hotspot Count",
        y="Brain Region",

        palette=sns.color_palette(

            "rocket",
            len(human_plot)
        ),

        ax=axes[0]
    )

    axes[0].set_title(

        "Human Regional Hotspots",

        color="FONT",

        fontsize=16,

        weight="bold"
    )

    axes[0].set_facecolor("black")

    axes[0].tick_params(

        colors="white",
        labelsize=10
    )

    axes[0].set_xlabel(

        "Hotspot Count",

        color="white"
    )

    axes[0].set_ylabel("")

    for spine in axes[0].spines.values():

        spine.set_visible(False)

    for i, v in enumerate(

        human_plot["Hotspot Count"]
    ):

        axes[0].text(

            v + 1,
            i,

            str(v),

            color="white",

            va="center",

            fontsize=9
        )

    # ----------------------------------------------
    # MOUSE
    # ----------------------------------------------

    sns.barplot(

        data=mouse_plot,

        x="Hotspot Count",
        y="Brain Region",

        palette=sns.color_palette(

            "mako",
            len(mouse_plot)
        ),

        ax=axes[1]
    )

    axes[1].set_title(

        "Mouse Regional Hotspots",

        color="white",

        fontsize=16,

        weight="bold"
    )

    axes[1].set_facecolor("black")

    axes[1].tick_params(

        colors="white",
        labelsize=10
    )

    axes[1].set_xlabel(

        "Hotspot Count",

        color="white"
    )

    axes[1].set_ylabel("")

    for spine in axes[1].spines.values():

        spine.set_visible(False)

    for i, v in enumerate(

        mouse_plot["Hotspot Count"]
    ):

        axes[1].text(

            v + 1,
            i,

            str(v),

            color="white",

            va="center",

            fontsize=9
        )

    plt.tight_layout()

    st.pyplot(fig)

# ======================================================
# FOOTER
# ======================================================
st.markdown("""
<div class='footer'>
UbiXplorer: A Proof-Of-Concept Prototype created @ Molecular Neuroscience and Functional Genomics Lab |
Department of Biotechnology | Delhi Technological University
</div>
""", unsafe_allow_html=True)
