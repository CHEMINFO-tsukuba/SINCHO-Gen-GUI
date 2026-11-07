import subprocess
import os
import shutil
from importlib import resources





def main():
    try:
        pkg_streamlit_dir = resources.files("env")/".streamlit"
    except Exception:
        pkg_streamlit_dir = None
    dst_dir = os.path.join(os.path.expanduser("~"), ".streamlit")
    if (
        pkg_streamlit_dir is not None
        and pkg_streamlit_dir.is_dir()
        and not os.path.exists(dst_dir)
    ):
        shutil.copytree(pkg_streamlit_dir, dst_dir)

    # ---- GUI 本体を起動 ----
    script = os.path.join(os.path.dirname(__file__), "gui_controller.py")
    port = os.environ.get("SINCHO_GUI_PORT", "8501")
    subprocess.run(
        [
            "streamlit",
            "run",
            script,
            "--server.address",
            "0.0.0.0",
            "--server.port",
            port,
        ],
        check=True,
    )


if __name__ == "__main__":
    main()
