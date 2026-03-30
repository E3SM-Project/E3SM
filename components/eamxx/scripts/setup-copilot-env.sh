#!/bin/bash
#
# Setup script for the copilot-testing machine.
# Installs dependencies needed to configure/build/test EAMxx standalone.
#
# Usage:
#   source components/eamxx/scripts/setup-copilot-env.sh
#
# This script:
#   1. Installs system packages (apt-get on Debian/Ubuntu, with spack fallback)
#   2. Fixes library paths for Debian/Ubuntu multiarch layout
#   3. Installs Python packages
#   4. Initializes required git submodules (uses HTTPS if SSH fails)
#   5. Exports environment variables needed by copilot-testing.cmake
#
# After running this script, you can configure with:
#   cd components/eamxx
#   ./scripts/test-all-eamxx -m copilot-testing -t dbg --config-only

# Determine whether we are sourced or executed directly.
_setup_copilot_is_sourced=false
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    _setup_copilot_is_sourced=true
fi

# Helper: abort with a message.  Uses 'return 1' when sourced so the
# caller's shell survives, 'exit 1' otherwise.
_setup_copilot_fail() {
    echo "SETUP ERROR: $1" >&2
    cd "$_setup_copilot_orig_dir" 2>/dev/null
    if $_setup_copilot_is_sourced; then return 1; else exit 1; fi
}

# Only enable 'errexit' when running as a script, not when sourced.
if ! $_setup_copilot_is_sourced; then
    set -e
fi

# Save the caller's working directory so we can restore it at the end.
_setup_copilot_orig_dir="$(pwd)"

# Ensure we always restore the working directory, even on failure.
# Use RETURN trap when sourced (scoped to this script), EXIT when executed.
if $_setup_copilot_is_sourced; then
    trap 'cd "$_setup_copilot_orig_dir" 2>/dev/null' RETURN
else
    trap 'cd "$_setup_copilot_orig_dir" 2>/dev/null' EXIT
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EAMXX_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
E3SM_ROOT="$(cd "$EAMXX_DIR/../.." && pwd)"

echo "=== EAMxx Copilot Testing Environment Setup ==="
echo "E3SM root: $E3SM_ROOT"

###############################################################################
# 1. Install system packages
###############################################################################

install_with_apt() {
    echo "--- Installing dependencies via apt-get ---"
    # Use sudo -E to preserve proxy env vars (needed in some CI environments)
    sudo -E apt-get update -qq || sudo apt-get update -qq
    sudo -E apt-get install -y -qq \
        gfortran gcc g++ cmake make git \
        python3 python3-pip \
        libopenmpi-dev openmpi-bin \
        libnetcdf-dev libnetcdff-dev libpnetcdf-dev \
        libboost-dev libyaml-cpp-dev \
        libblas-dev liblapack-dev \
        libcurl4-openssl-dev zlib1g-dev \
        libhdf5-openmpi-dev \
        perl \
    || sudo apt-get install -y \
        gfortran gcc g++ cmake make git \
        python3 python3-pip \
        libopenmpi-dev openmpi-bin \
        libnetcdf-dev libnetcdff-dev libpnetcdf-dev \
        libboost-dev libyaml-cpp-dev \
        libblas-dev liblapack-dev \
        libcurl4-openssl-dev zlib1g-dev \
        libhdf5-openmpi-dev \
        perl
}

fix_multiarch_paths() {
    # Debian/Ubuntu multiarch puts libraries in lib/<triplet>/ instead of lib/.
    # Scorpio's FindNetCDF/FindPnetCDF only searches PATH_SUFFIXES "lib" with
    # NO_DEFAULT_PATH, so cmake can't find them. Create symlinks in /usr/lib/ to fix.
    echo "--- Fixing multiarch library paths ---"
    local arch_dir
    arch_dir="/usr/lib/$(dpkg-architecture -qDEB_HOST_MULTIARCH 2>/dev/null || echo x86_64-linux-gnu)"
    local target_dir="/usr/lib"
    for lib in libnetcdf.so libnetcdff.so libpnetcdf.so libhdf5.so libhdf5_hl.so; do
        if [ -f "$arch_dir/$lib" ] && [ ! -f "$target_dir/$lib" ]; then
            sudo ln -sf "$arch_dir/$lib" "$target_dir/$lib"
            echo "  Linked $lib"
        fi
    done
}

install_with_spack() {
    echo "--- Installing dependencies via spack ---"

    # Determine Spack root from existing installation when possible.
    local spack_root=""
    if command -v spack &> /dev/null; then
        spack_root="${SPACK_ROOT:-$(spack location -r 2>/dev/null || true)}"
    fi

    # Fall back to $HOME/spack, potentially cloning a fresh install.
    if [ -z "$spack_root" ]; then
        if [ -d "$HOME/spack" ]; then
            echo "Found existing spack installation at $HOME/spack, reusing it"
        else
            echo "Installing spack..."
            git clone --depth=2 https://github.com/spack/spack.git "$HOME/spack"
        fi
        spack_root="$HOME/spack"
    fi

    # Always source the setup script so 'spack load' updates the environment.
    local spack_setup="${spack_root}/share/spack/setup-env.sh"
    if [ ! -f "$spack_setup" ]; then
        _setup_copilot_fail "Spack setup script not found at: $spack_setup"
    fi
    . "$spack_setup"

    spack install --no-checksum \
        cmake \
        openmpi \
        netcdf-c +mpi \
        netcdf-fortran \
        parallel-netcdf \
        boost \
        yaml-cpp \
        openblas

    spack load cmake openmpi netcdf-c netcdf-fortran parallel-netcdf boost yaml-cpp openblas
}

if command -v apt-get &> /dev/null; then
    install_with_apt || _setup_copilot_fail "apt-get package installation failed"
    fix_multiarch_paths
elif command -v spack &> /dev/null || [ -d "$HOME/spack" ]; then
    install_with_spack || _setup_copilot_fail "spack package installation failed"
else
    echo "WARNING: Neither apt-get nor spack found."
    echo "Please install the following manually:"
    echo "  GNU compilers (gcc, g++, gfortran), cmake, make, MPI (openmpi),"
    echo "  netcdf-c, netcdf-fortran, pnetcdf, boost, yaml-cpp, blas/lapack"
    echo "Then re-run this script or set up the environment manually."
fi

###############################################################################
# 2. Install Python packages
###############################################################################

echo "--- Installing Python packages ---"
python3 -m pip install --quiet psutil pyyaml netCDF4 packaging \
    || _setup_copilot_fail "pip package installation failed"

###############################################################################
# 3. Initialize git submodules
###############################################################################

echo "--- Initializing required git submodules ---"
cd "$E3SM_ROOT"

# Use HTTPS instead of SSH if SSH is unavailable (common in CI containers)
if ! ssh -o BatchMode=yes -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
       -o ConnectTimeout=5 -T git@github.com 2>&1 | grep -q "successfully authenticated"; then
    git config url."https://github.com/".insteadOf "git@github.com:"
    echo "  Using HTTPS for GitHub (SSH unavailable)"
fi

git submodule update --init --recursive \
    externals/ekat externals/scorpio externals/mam4xx externals/haero \
    components/eam/src/physics/cosp2/external \
    components/eam/src/physics/rrtmgp/external \
    || _setup_copilot_fail "git submodule init failed"
git submodule update --init cime \
    || _setup_copilot_fail "cime submodule init failed"
cd cime && git submodule update --init CIME/non_py/cprnc \
    || _setup_copilot_fail "cprnc submodule init failed"

###############################################################################
# 4. Set environment variables
###############################################################################

# Restore the caller's working directory before exporting env vars.
cd "$_setup_copilot_orig_dir"

echo "--- Setting environment variables ---"

# Input data directory
export SCREAM_INPUT_ROOT="${SCREAM_INPUT_ROOT:-$HOME/e3sm-inputdata}"
mkdir -p "$SCREAM_INPUT_ROOT"

echo "=== Setup complete ==="
echo ""
echo "Environment variables set:"
echo "  SCREAM_INPUT_ROOT=$SCREAM_INPUT_ROOT"
echo ""
echo "Next steps:"
echo "  cd components/eamxx"
echo "  ./scripts/test-all-eamxx -m copilot-testing -t dbg --config-only"
echo ""
echo "Then for incremental build/test:"
echo "  cd ctest-build/copilot-testing/full_debug"
echo "  make -j\$(nproc) && ctest -j\$(nproc)"
