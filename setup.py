from setuptools import setup
from pathlib import Path

setup(name="sincho-gen-gui",
      version=Path(__file__).with_name("VERSION").read_text().strip(),
      packages=["GUI_Util", "env"],
      py_modules=["main", "gui_controller"],
      entry_points={
          "console_scripts":[
              "sincho-gen-gui = main:main"
          ]
      },
      package_data={
          "GUI_Util": ["*.png", "*.yaml"],
          "env": [".streamlit/*"],
      },
      include_package_data=True,
      install_requires=[
          "streamlit",
            "biopython",
            "py3Dmol",
          "pandas<2",
            "numpy",
            "rdkit",
            "pyyaml",
            "matplotlib",
      ],)
