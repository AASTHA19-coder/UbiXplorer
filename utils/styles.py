import streamlit as st

def load_css():

    st.markdown("""
    <style>

    .stApp{
        background: linear-gradient(180deg,#f9fffc 0%, #eef8f4 100%);
        font-family:Segoe UI;
    }

    .brand-title{
        font-size:42px;
        font-weight:800;
        color:#0f766e;
        margin-top:10px;
    }

    .brand-sub{
        font-size:17px;
        color:#475569;
        margin-top:-8px;
    }

    .hero-box{
        background:white;
        padding:45px;
        border-radius:24px;
        box-shadow:0 12px 28px rgba(0,0,0,0.06);
        text-align:center;
    }

    .hero-box h1{
        color:#134e4a;
        font-size:44px;
        margin-bottom:10px;
    }

    .hero-box p{
        font-size:20px;
        color:#475569;
    }

    .mini-card{
        background:white;
        padding:22px;
        border-radius:18px;
        text-align:center;
        box-shadow:0 8px 20px rgba(0,0,0,0.05);
        font-size:18px;
    }

    </style>
    """, unsafe_allow_html=True)