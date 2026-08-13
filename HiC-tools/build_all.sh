#!/bin/bash
# ==============================================================================
# build_all.sh
# Build every tool's Docker image, then (optionally) convert each to a
# Singularity .sif so config/containers.env can point at them directly.
#
# Usage:
#   ./build_all.sh                 # build Docker images only
#   ./build_all.sh --to-singularity /opt/containers   # also convert + place .sif files
#
# Requires: docker (with buildkit), and singularity/apptainer + the
# docker-daemon:// build source if using --to-singularity (must run as a
# user who can both `docker build` and `singularity build`, e.g. on the
# same machine, or via a CI runner with Docker-in-Docker).
# ==============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

declare -A IMAGES=(
    [hicpro]="Dockerfile.hicpro:hicpro:3.1.0"
    [hicup]="Dockerfile.hicup:hicup:0.9.2"
    [juicer]="Dockerfile.juicer:juicer:1.9"
    [hicpipe]="Dockerfile.hicpipe:hicpipe:0.93"
    [hickit]="Dockerfile.hickit:hickit:latest"
    [bin3c]="Dockerfile.bin3c:bin3c:latest"
    [fanc]="Dockerfile.fanc:fanc:0.9.25"
)

SIF_OUTDIR=""
if [[ "${1:-}" == "--to-singularity" ]]; then
    SIF_OUTDIR="${2:?Provide an output directory for .sif files, e.g. ./build_all.sh --to-singularity /opt/containers}"
    mkdir -p "${SIF_OUTDIR}"
fi

for name in "${!IMAGES[@]}"; do
    IFS=":" read -r dockerfile repo tag <<< "${IMAGES[$name]}"
    image="${repo}:${tag}"

    echo "=============================================================="
    echo " Building ${name}  ->  ${image}  (from ${dockerfile})"
    echo "=============================================================="
    docker build -t "${image}" -f "${dockerfile}" .

    if [[ -n "${SIF_OUTDIR}" ]]; then
        sif_name="${SIF_OUTDIR}/${name}_$(echo "${tag}" | tr '/' '_').sif"
        echo "  Converting to Singularity: ${sif_name}"
        singularity build --force "${sif_name}" "docker-daemon://${image}"
    fi
done

echo "All images built."
if [[ -n "${SIF_OUTDIR}" ]]; then
    echo "Singularity images written to: ${SIF_OUTDIR}"
    echo "Point config/containers.env at these paths, e.g.:"
    echo "  export HICPRO_SIF=${SIF_OUTDIR}/hicpro_3.1.0.sif"
    echo "  export JUICER_SIF=${SIF_OUTDIR}/juicer_1.9.sif"
fi
