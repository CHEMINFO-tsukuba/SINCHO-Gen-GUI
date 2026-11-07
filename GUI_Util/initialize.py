import streamlit as st
import subprocess
import os
import tempfile
import io
from Bio.PDB import PDBParser
from fractions import Fraction
import py3Dmol
from streamlit.components.v1 import html

#Streamlitのセッション状態を初期化する関数
def init():
    # Streamlit のセッション状態を初期化
    if "uploaded_pdb_file" not in st.session_state:
        st.session_state.uploaded_pdb_file = None
    if "residues_list" not in st.session_state:
        st.session_state.residues_list = []
    if "hit_residue" not in st.session_state:
        st.session_state.hit_residue = None
    if "edited_pdb_content" not in st.session_state:
        st.session_state.edited_pdb_content = None
    if "md_settings" not in st.session_state:
        st.session_state.md_settings = {}
    if "p2c_sincho_settings" not in st.session_state:
        st.session_state.p2c_sincho_settings={}
    if "chemts_settings" not in st.session_state:
        st.session_state.chemts_settings = {}
    if "aascore_settings" not in st.session_state:
        st.session_state.aascore_settings = {}
    if "general_settings" not in st.session_state:
        st.session_state.general_settings = {}
    if "output_settings" not in st.session_state:
        st.session_state.output_settings = {}
    if "yaml_content" not in st.session_state:
        st.session_state.yaml_content = None
    st.session_state.general_settings.setdefault("dir_checker", False)
    st.session_state.general_settings.setdefault("dir_overwrite_confirm", False)

    st.session_state.input_state = False  # 初期化フラグ

    st.set_page_config(page_title="SINCHO-Gen-GUI", page_icon=os.path.join(os.path.dirname(__file__), "GEN_icon.png"), layout="wide")

    #if "main_tab" not in st.session_state:
    #    st.session_state["main_tab"] = "Home"


    return st
