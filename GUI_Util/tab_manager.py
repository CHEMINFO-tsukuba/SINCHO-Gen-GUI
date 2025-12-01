import streamlit as st
import subprocess
import os
import tempfile
import io
from Bio.PDB import PDBParser
from fractions import Fraction
import py3Dmol
from streamlit.components.v1 import html

"""
#Streamlitのセッション状態を初期化する関数
def init_tabs():

    st.sidebar.title("TABS")
    
    main_tab = st.sidebar.radio("Main Panel",["Input", "Output"])
    st.sidebar.image(os.path.join(os.path.dirname(__file__), main_tab+".png"), use_column_width=True)

    # 各タブ内で階層2を定義
    if main_tab == "Input":
        sub_tab = st.sidebar.radio("Details", ["General", "Upload Complex", "Select Hit Ligand", "MD Settings", "SINCHO Settings", "ChemTS Settings", "AAScore Settings", "Summary"])

            

    elif main_tab == "Output":
        sub_tab = st.sidebar.radio("Details", ["General", "MD", "SINCHO", "ChemTS", "AAScore"])


    return main_tab, sub_tab
"""

# tab_manager.py

import streamlit as st
import os
import base64

def init_tabs():

    # 初期化
    if "main_tab" not in st.session_state:
        st.session_state["main_tab"] = "Home"

    main_tab = st.session_state["main_tab"]

    # 階層2：Details
    if main_tab == "Home":
        sub_tab = None

        root_col1, root_col2 = st.columns([1, 1])
        with root_col1:
            try:
                st.image(os.path.join(os.path.dirname(__file__),"Input.png"), use_container_width=True)
            except Exception as e:
                st.image(os.path.join(os.path.dirname(__file__),"Input.png"), use_column_width=True)
            if st.button("Go to Input Generator"):
                st.session_state["main_tab"] = "Input"
                try:
                    st.rerun()
                except Exception as e:
                    st.experimental_rerun()
        with root_col2:
            try:
                st.image(os.path.join(os.path.dirname(__file__),"Output.png"), use_container_width=True)
            except Exception as e:
                st.image(os.path.join(os.path.dirname(__file__),"Output.png"), use_column_width=True)
            if st.button("Go to Output Analyzer"):
                st.session_state["main_tab"] = "Output"
                try:
                    st.rerun()
                except Exception as e:
                    st.experimental_rerun()
        



    elif main_tab == "Input":
        st.sidebar.title("")
        # 表示画像
        try:
            st.sidebar.image(os.path.join(os.path.dirname(__file__), main_tab + ".png"), use_container_width=True)
        except Exception as e:
            st.sidebar.image(os.path.join(os.path.dirname(__file__), main_tab + ".png"), use_column_width=True)
        sub_tab = st.sidebar.radio("Details", [
            "General", "Initial Upload", "Hit Residue Selection", "MD Settings",
            "SINCHO Settings", "ChemTS Settings", "AAScore Settings", "Summary"
        ])
    elif main_tab == "Output":
        st.sidebar.title("")
        # 表示画像
        try:
            st.sidebar.image(os.path.join(os.path.dirname(__file__), main_tab + ".png"), use_container_width=True)
        except Exception as e:
            st.sidebar.image(os.path.join(os.path.dirname(__file__), main_tab + ".png"), use_column_width=True)
        sub_tab = st.sidebar.radio("Details", [
            "General", "MD", "SINCHO", "ChemTS", "AAScore"
        ])

    # 切り替えボタン
    if main_tab == "Input":
        if st.sidebar.button("Change to Output Analyzer"):
            st.session_state["main_tab"] = "Output"
            try:
                st.rerun()
            except Exception as e:
                st.experimental_rerun()
    elif main_tab == "Output":
        if st.sidebar.button("Change to Input Generator"):
            st.session_state["main_tab"] = "Input"
            try:
                st.rerun()
            except Exception as e:
                st.experimental_rerun()



    return main_tab, sub_tab


