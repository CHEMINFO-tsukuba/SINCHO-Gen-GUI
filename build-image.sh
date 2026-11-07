#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
version="$(tr -d '[:space:]' < "${script_dir}/VERSION")"

base_image="${SINCHO_GUI_BASE_IMAGE:-sincho-gen:${version}}"
output_image="${SINCHO_GUI_IMAGE:-sincho-gen-gui:${version}}"

if ! docker image inspect "${base_image}" >/dev/null 2>&1; then
  echo "[ERROR] SINCHO-Gen ${version}の基盤イメージがありません: ${base_image}" >&2
  exit 1
fi

docker build \
  --build-arg "SINCHO_BASE_IMAGE=${base_image}" \
  --build-arg "SINCHO_VERSION=${version}" \
  -t "${output_image}" \
  "${script_dir}"
