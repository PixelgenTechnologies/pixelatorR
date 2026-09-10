#!/usr/bin/env bash
# Idempotent Cloud Agent setup for pixelatorR.
#
# Installs micromamba + task, creates the conda environment described in
# environment.yml (the repository's recommended Linux dev setup), adds the
# optional Suggests packages required to run the full test suite, and installs
# pixelatorR from the checked-out source. Safe to run repeatedly.
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOCAL_BIN="$HOME/.local/bin"
ENV_NAME="r-pixelator"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
mkdir -p "$LOCAL_BIN" "$MAMBA_ROOT_PREFIX"
export PATH="$LOCAL_BIN:$PATH"

# 1. micromamba (user-level, no root required)
if [ ! -x "$LOCAL_BIN/micromamba" ]; then
  echo "Installing micromamba..."
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -C "$HOME" -xj bin/micromamba
  mv "$HOME/bin/micromamba" "$LOCAL_BIN/micromamba"
  rmdir "$HOME/bin" 2>/dev/null || true
fi

# 2. task (go-task) so the documented DEVELOPERS.md tasks are available
if [ ! -x "$LOCAL_BIN/task" ]; then
  echo "Installing task..."
  sh -c "$(curl -fsSL https://taskfile.dev/install.sh)" -- -d -b "$LOCAL_BIN"
fi

# 3. Create/update the conda environment from environment.yml
if micromamba env list | grep -qE "\b${ENV_NAME}\b"; then
  echo "Updating existing '${ENV_NAME}' environment from environment.yml..."
  micromamba install -y -n "$ENV_NAME" -f "$REPO_DIR/environment.yml"
else
  echo "Creating '${ENV_NAME}' environment from environment.yml..."
  micromamba env create -y -n "$ENV_NAME" -f "$REPO_DIR/environment.yml"
fi

# 4. Optional (Suggests) packages needed to run the complete test suite.
#    environment.yml only pins the core dependencies; these extras let
#    devtools::test() and R CMD check exercise every code path.
echo "Installing optional Suggests packages..."
micromamba install -y -n "$ENV_NAME" -c conda-forge -c bioconda \
  bioconductor-complexheatmap bioconductor-pcamethods bioconductor-sparsematrixstats \
  r-reticulate r-dtplyr r-data.table r-rcppannoy r-gifski r-av r-magick r-ragg r-png \
  r-rcolorbrewer r-rcppml r-ggrepel r-rspectra r-irlba r-pls r-matrixstats r-fnn

# 5. Install pixelatorR itself from the checked-out source (deps already present)
echo "Installing pixelatorR from source..."
micromamba run -n "$ENV_NAME" R CMD INSTALL --no-multiarch --with-keep.source "$REPO_DIR"

# 6. Make the environment active in interactive shells
BASHRC="$HOME/.bashrc"
MARKER="# >>> pixelatorR dev environment >>>"
if ! grep -qF "$MARKER" "$BASHRC" 2>/dev/null; then
  {
    echo ""
    echo "$MARKER"
    echo "export PATH=\"$LOCAL_BIN:\$PATH\""
    echo "export MAMBA_ROOT_PREFIX=\"$MAMBA_ROOT_PREFIX\""
    echo "eval \"\$($LOCAL_BIN/micromamba shell hook --shell bash)\""
    echo "micromamba activate $ENV_NAME 2>/dev/null || true"
    echo "# <<< pixelatorR dev environment <<<"
  } >> "$BASHRC"
fi

echo "pixelatorR development environment is ready (conda env: ${ENV_NAME})."
