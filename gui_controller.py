import streamlit as st
import os
import tempfile
import io
import base64
from Bio.PDB import PDBParser
from fractions import Fraction
import py3Dmol
from streamlit.components.v1 import html

from GUI_Util import initialize, tab_manager
from GUI_Util.input_controller import InputController
from GUI_Util.output_controller import OutputController



#Initialize
initialize.init()

st.markdown(
    """
    <style>
    .stButton>button {
        background-color: #6B8E23;
        color: white;
        border-radius: 8px;
        border: none;
    }
    .stButton>button:hover {
        background-color: #A9BA9D;
        transform: scale(1.05);
    }
    </style>
    """,
    unsafe_allow_html=True
)

try:
    st.image(os.path.join(os.path.dirname(__file__),"GUI_Util/sincho-gen_logo_transparent.png"), use_container_width=True)
except Exception as e:
    st.image(os.path.join(os.path.dirname(__file__),"GUI_Util/sincho-gen_logo_transparent.png"), use_column_width=True)



main_tab, sub_tab = tab_manager.init_tabs()
if main_tab == "Input":
    InputController().process(sub_tab)

if main_tab == "Output":
    OutputController().process(sub_tab)
        
