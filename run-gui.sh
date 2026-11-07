#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
version="$(tr -d '[:space:]' < "${script_dir}/VERSION")"
name="${SINCHO_GUI_NAME:-sincho-gen-gui}"
image="${SINCHO_GUI_IMAGE:-sincho-gen-gui:${version}}"
port="${SINCHO_GUI_HOST_PORT:-8501}"

gpu_args=()
if [ -n "${SINCHO_GEN_GPU:-}" ]; then
  gpu_args=(--gpus "device=${SINCHO_GEN_GPU}")
fi

if docker ps -a --format '{{.Names}}' | grep -q "^${name}$"; then
  echo "[ERROR] A container with this name already exists: ${name}" >&2
  echo "Run docker rm -f ${name}, then try again." >&2
  exit 1
fi

docker run -d \
  --name "${name}" \
  "${gpu_args[@]}" \
  -e USER_ID="$(id -u)" \
  -e GROUP_ID="$(id -g)" \
  -p "${port}:8501" \
  -v "$(pwd):/work" \
  "${image}"

echo "SINCHO-Gen-GUI: http://localhost:${port}"
echo "View logs: docker logs -f ${name}"
