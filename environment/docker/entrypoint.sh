#!/bin/bash
set -e


# -------------------------
# Welcome message
# -------------------------
echo " ~~~~~~~  Welcome to the MURAVES Containter 0.0 ~~~~~~~~  "
echo "                      .-----. "
echo "             .----. .'       ' "
echo "            '      V           '  "
echo "          '                      ' "
echo "        '                          '   "
echo "      '                              ' "
echo "       _  _        _   _        _  _  "
echo "      |  V | |  | |_| |_| \\  / |_ |_  "
echo "      |    | |__| | \\ | |  \\/  |_  _| "
echo
echo " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  "


#-------------------------
# Optional: sanity check
# -------------------------
echo "Sanity check"
echo "Using python: $(type -p python)"
echo "Conda env: $CONDA_DEFAULT_ENV"

# -------------------------
# Install muraves_lib (editable mode): il container docker è configurato in modo tale che python/pip/jupiter puntano già a all'ambiente muraves.
#                               pertanto non è necessario attivarlo esplicitamente.
# -------------------------
if [ -d "/workspace/Software" ]; then
    cd /workspace/Software
    python -m pip install -e . --quiet
fi
pip install -U scikit-learn --quiet
#------------------------------------------------
# Check if /workspace is empty or not a git repo
#------------------------------------------------
if [ ! -d "/workspace/.git" ]; then
    echo "[WARNING:] No Git repo found in /workspace. "
    # git clone --depth 1 https://github.com/muraves/Software.git /workspace
    # chmod -R a-w /workspace   # make read-only
else
    echo "Detected existing Git repo in /workspace."
fi


# ---------------------------------------------------------
# Open interactive shell with muraves conda environment: 
# ---------------------------------------------------------
echo "[INFO] Container ready. Starting shell..."
if [ $# -gt 0 ]; then
    source /opt/conda/etc/profile.d/conda.sh && conda activate muraves
    exec "$@"
else
    exec bash --init-file <(echo "source /opt/conda/etc/profile.d/conda.sh && conda activate muraves")
fi


