#!/bin/bash
set -e
# -----------------------------
# Wrapper script to run muraves-env.sif on T2B
# -----------------------------
# This file isn't needed to build the image.
# It is needed to
# - run the container
# - automatically activate the mamba environment 
# - Check existance of Git repository in the binded worspace

echo "=== Muraves environment setup ==="

# Ask user for workspace location
read -rp "Enter workspace location [default: \$HOME]: " WORKSPACE_INPUT

# Use HOME if nothing is provided
if [ -z "$WORKSPACE_INPUT" ]; then
    WORKSPACE_HOST="$HOME"
else
    # Expand ~ if used
    WORKSPACE_HOST="$(eval echo "$WORKSPACE_INPUT")"
fi

export WORKSPACE_HOST
echo "Using workspace: $WORKSPACE_HOST"
# Create workspace if it doesn't exist
mkdir -p "$WORKSPACE_HOST"

# Clone repository if not present
REPO_URL="https://github.com/muraves/Software.git"
WORKSPACE_REPO="$WORKSPACE_HOST/Software"
if [ ! -d "$WORKSPACE_REPO" ]; then
    echo "Cloning Muraves repository..."
    git clone "$REPO_URL" "$WORKSPACE_REPO"
else
    echo "Repository already exists: $WORKSPACE_REPO"
fi

# Move into repo (important if script is run from elsewhere)
cd "$WORKSPACE_REPO"


# Ensure scripts are executable
chmod +x environment/docker/entrypoint.sh
#chmod +x environment/docker/run_muraves_env.sh

# Path to your Singularity image
SIF_PATH="/group/Muography/MURAVES/container/muraves-env.sif"

# Bind your workspace (host repo directory)
#WORKSPACE_HOST="path/to/your/git/repo"
WORKSPACE_CONTAINER="/workspace"

DIR_HOST="/pnfs/iihe/muraves/muraves_DATA"
DIR_CONTAINER="/data"

OUTPUT_HOST="/group/Muography/MURAVES/"
OUTPUT_CONTAINER="/outputs"



# If no command is passed, default to bash
if [ $# -eq 0 ]; then
    CMD="bash"
else
    CMD="$@"
fi


# Execute entrypoint inside the container
singularity exec --bind ${WORKSPACE_HOST}:${WORKSPACE_CONTAINER},${DIR_HOST}:${DIR_CONTAINER},${OUTPUT_HOST}:${OUTPUT_CONTAINER},  \
    --pwd ${WORKSPACE_CONTAINER} \
    ${SIF_PATH} /workspace/Software/environment/docker/entrypoint.sh ${CMD}
    
    # this would be the correct way to call it when entrypoint.sh is definitive.
    #${SIF_PATH} /usr/local/bin/entrypoint.sh ${CMD}


