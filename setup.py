from setuptools import setup

setup(name="sincho-gen-gui",
      version="0.1.0",
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
            "pandas",
            "numpy",
            "rdkit",
            "pyyaml",
            "matplotlib",
      ],)
