#!/bin/bash
# init.sh – environment setup for FrameD (Linux and macOS)
# Equivalent of init.csh but written in bash for broader platform support.
set -e

# Resolve the directory containing this script so it can be run from anywhere.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ---------------------------------------------------------------------------
# Environment variables (bash equivalent of dnastorage.env)
# ---------------------------------------------------------------------------
export DNASTORAGE_HOME="$PWD"
export DNASTORAGE_TOOLS="$PWD/tools"
export DNASTORAGE_LSF="$PWD/tools/lsf"
export DNArSimPath="$PWD/dnastorage/fi/DNArSim"
export FRAMED_CONFIGS="$PWD/examples"
export FRAMED_CONDA="$PWD/framed_conda"
export DNArSimData="$PWD/probEdit"
export FRAMED_IMAGE_FILES="$PWD/test_files"

# ---------------------------------------------------------------------------
# Clone DNArSim, apply patches, then remove the clone
# ---------------------------------------------------------------------------
git clone https://github.com/BHam-1/DNArSim/

rm -rf probEdit
mv DNArSim/simulator/probEdit "$PWD"

for i in channel functions loadProb; do
    patch_file="${DNArSimPath}/${i}.patch"
    original_file="DNArSim/simulator/${i}.jl"
    output_file="${DNArSimPath}/${i}.jl"
    patch -o "$output_file" "$original_file" "$patch_file"
done

rm -rf DNArSim

# ---------------------------------------------------------------------------
# Initialise git submodules (or clone schwimmbad if not inside a git repo)
# ---------------------------------------------------------------------------
if git rev-parse --git-dir > /dev/null 2>&1; then
    git submodule update --init --recursive
else
    git clone https://github.com/kvolkel/schwimmbad
fi

if [ "${1:-}" = "-no-env" ]; then
    exit 0
fi

# ---------------------------------------------------------------------------
# Create the conda environment
# ---------------------------------------------------------------------------
conda env create --prefix "$FRAMED_CONDA" --file "$PWD/dnastorage.yml"

# ---------------------------------------------------------------------------
# Install pip-only packages and the schwimmbad submodule inside the env.
# conda run is used so that this script works without interactive conda init.
# ---------------------------------------------------------------------------
conda run --live-stream --prefix "$FRAMED_CONDA" pip install -r requirements.txt
conda run --live-stream --prefix "$FRAMED_CONDA" pip install ./schwimmbad
conda run --live-stream --prefix "$FRAMED_CONDA" python -c "import julia; julia.install()"
