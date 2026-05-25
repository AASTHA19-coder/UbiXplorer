import streamlit as st
from PIL import Image

st.set_page_config(
    page_title="UbiXplorer",
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
max-width:1250px;
padding-top:2rem;
padding-bottom:2rem;
}

section[data-testid="stSidebar"]{
background:#0b1120;
border-right:1px solid rgba(255,255,255,0.05);
}

.hero{
padding:35px;
border-radius:24px;
background:rgba(255,255,255,0.03);
border:1px solid rgba(45,212,191,0.18);
box-shadow:0 0 25px rgba(45,212,191,0.08);
}

.title{
font-size:58px;
font-weight:900;
color:#2dd4bf;
margin-bottom:0px;
}

.subtitle{
font-size:22px;
color:#cbd5e1;
margin-top:-8px;
}

.tag{
display:inline-block;
padding:8px 14px;
border-radius:20px;
background:rgba(45,212,191,0.12);
color:#5eead4;
font-size:14px;
font-weight:600;
margin-top:12px;
}

.section{
font-size:34px;
font-weight:800;
margin-top:35px;
margin-bottom:18px;
color:white;
}

.card{
background:rgba(255,255,255,0.03);
padding:24px;
border-radius:20px;
border:1px solid rgba(255,255,255,0.06);
height:190px;
transition:0.3s;
}

.card:hover{
border:1px solid rgba(45,212,191,0.5);
box-shadow:0 0 20px rgba(45,212,191,0.12);
}

.card h3{
color:#2dd4bf;
}

.small{
color:#cbd5e1;
font-size:16px;
}

.metric{
background:rgba(45,212,191,0.08);
padding:20px;
border-radius:18px;
text-align:center;
border:1px solid rgba(45,212,191,0.15);
}

.metricnum{
font-size:30px;
font-weight:800;
color:#5eead4;
}

.metrictxt{
font-size:14px;
color:#cbd5e1;
}

.ribbon{
padding:18px;
border-radius:18px;
background:rgba(255,255,255,0.03);
border:1px solid rgba(255,255,255,0.06);
font-size:18px;
color:#e2e8f0;
text-align:center;
}

.footer{
text-align:center;
color:#94a3b8;
font-size:14px;
margin-top:45px;
}

</style>
""", unsafe_allow_html=True)

# ======================================================
# HEADER
# ======================================================
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
ASSETS = BASE_DIR / "assets"

logo_path = ASSETS / "logo.svg"
#mnfg_path = ASSETS / "MNFG.svg"

a, b, c = st.columns([1.2,5,1.2], vertical_alignment="center")

with a:
    st.image(str(logo_path), width=240)

with b:
    st.markdown("<div class='title'>UbiXplorer</div>", unsafe_allow_html=True)
    st.markdown("<div class='subtitle'>Ubiquitin-Centric Spatial Analysis Platform</div>", unsafe_allow_html=True)
    st.markdown("<div class='tag'>Ubiquitin | Neurobiology | Cross-Species Analysis</div>", unsafe_allow_html=True)

#with c:
#    st.markdown("<div style='margin-top:18px;'></div>", unsafe_allow_html=True)
#    st.image(str(mnfg_path), width=170)
# ======================================================
# HERO
# ======================================================
st.write("")

st.markdown("""
<div class='hero'>

<h2 style='font-size:38px;color:white;margin-bottom:10px;'>
Decode Hidden Ubiquitin Biology from Gene Expression Data
</h2>

<p style='font-size:24px;color:#cbd5e1;'>

Upload DEG datasets, prioritize Ub-linked hub genes, visualize Human & Mouse hotspot maps, and interpret pathway intelligence instantly.

</p>

</div>
""", unsafe_allow_html=True)

# ======================================================
# CARDS
# ======================================================
st.markdown("<div class='section'>Conceptual Overview</div>", unsafe_allow_html=True)

c1,c2,c3,c4 = st.columns(4)

with c1:
    st.markdown("""
    <div class='card'>
    <h3>Upload Data</h3>
    <p class='small'>CSV / XLSX differential expression datasets</p>
    </div>
    """, unsafe_allow_html=True)

with c2:
    st.markdown("""
    <div class='card'>
    <h3>Hub Ranking</h3>
    <p class='small'>Ub-specific prioritization algorithm</p>
    </div>
    """, unsafe_allow_html=True)

with c3:
    st.markdown("""
    <div class='card'>
    <h3>Spatial Maps</h3>
    <p class='small'>Human + Mouse transcriptomic atlases</p>
    </div>
    """, unsafe_allow_html=True)

with c4:
    st.markdown("""
    <div class='card'>
    <h3>Pathways</h3>
    <p class='small'>STRING / KEGG intelligence layer</p>
    </div>
    """, unsafe_allow_html=True)

# ======================================================
# METRICS
# ======================================================
st.markdown("<div class='section'>Platform Summary </div>", unsafe_allow_html=True)

m1,m2,m3,m4 = st.columns(4)

with m1:
    st.markdown("""
    <div class='metric'>
    <div class='metricnum'>2</div>
    <div class='metrictxt'>Species Supported</div>
    </div>
    """, unsafe_allow_html=True)

with m2:
    st.markdown("""
    <div class='metric'>
    <div class='metricnum'>6+</div>
    <div class='metrictxt'>Hub Genes Validated</div>
    </div>
    """, unsafe_allow_html=True)

with m3:
    st.markdown("""
    <div class='metric'>
    <div class='metricnum'>1000+</div>
    <div class='metrictxt'>Spatial Spots Analysed</div>
    </div>
    """, unsafe_allow_html=True)

with m4:
    st.markdown("""
    <div class='metric'>
    <div class='metricnum'>Ubiquitin</div>
    <div class='metrictxt'>Priority Engine</div>
    </div>
    """, unsafe_allow_html=True)

# ======================================================


# ======================================================
# CTA
# ======================================================
st.markdown("<div class='section'>Next Step</div>", unsafe_allow_html=True)

st.info("Open 'Upload Data' from the side pane to continue further.")

# ======================================================
# FOOTER
# ======================================================
st.markdown("""
<div class='footer'>
UbiXplorer: A Proof-Of-Concept Prototype created @ Molecular Neuroscience and Functional Genomics Lab   |     
                Department of Biotechnology |   Delhi Technological University
</div>
""", unsafe_allow_html=True)
