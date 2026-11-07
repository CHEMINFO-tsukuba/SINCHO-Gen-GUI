# SINCHO-Gen-GUI

A Streamlit GUI add-in for SINCHO-Gen 0.1.0. It supports inspecting input PDB files, generating configuration YAML files, and visualizing MD, SINCHO, ChemTS, and AA Score results.

The GUI is maintained as an input/output add-in with the same version as SINCHO-Gen. Both releases are version `0.1.0`. The GUI image is built on the SINCHO-Gen image with the matching version tag.

## Supported SINCHO-Gen version

- Public model image: `sincho-gen:0.1.0`
- Computation driver: `/SINCHO-Gen/extend_driver.py`

The GUI uses the SINCHO-Gen computation driver installed in the same container. After generating the YAML configuration, run the driver as described below.

## Build the GUI image

```bash
cd /path/to/SINCHO-Gen-GUI
./build-image.sh
```

This creates the following image:

```text
sincho-gen-gui:0.1.0
```

To specify the base and output images explicitly:

```bash
SINCHO_GUI_BASE_IMAGE=sincho-gen:0.1.0 \
SINCHO_GUI_IMAGE=sincho-gen-gui:0.1.0 \
./build-image.sh
```

## Launch the GUI

Launch the GUI from the directory where you want to store calculation files.

A GPU is not required to use the GUI alone. By default, no GPU is assigned to the container.

```bash
cd /path/to/job
/path/to/SINCHO-Gen-GUI/run-gui.sh
```

Open the following URL in your browser:

```text
http://localhost:8501
```

You can change the port and image using environment variables. To run GPU computations in the same container, set `SINCHO_GEN_GPU` to the GPU index when launching it. This requires an NVIDIA driver and Docker GPU support on the host.

```bash
SINCHO_GUI_HOST_PORT=8502 \
SINCHO_GEN_GPU=1 \
SINCHO_GUI_IMAGE=sincho-gen-gui:0.1.0-private \
/path/to/SINCHO-Gen-GUI/run-gui.sh
```

## Start a calculation

After saving the YAML file on the GUI's Summary tab, open a shell in the container.
For calculations that require a GPU, launch the container with `SINCHO_GEN_GPU` set as shown above.

Use the filename you saved in the GUI in place of `conditions.yaml` below.

```bash
docker exec -it sincho-gen-gui /usr/sbin/gosu dockeruser bash
python /SINCHO-Gen/extend_driver.py conditions.yaml
```

To resume a calculation while skipping completed stages:

```bash
python /SINCHO-Gen/extend_driver.py conditions.yaml --skip-md
python /SINCHO-Gen/extend_driver.py conditions.yaml --skip-sincho
python /SINCHO-Gen/extend_driver.py conditions.yaml --skip-chemts
```

## Stop and remove the container

```bash
docker stop sincho-gen-gui
docker rm sincho-gen-gui
```
