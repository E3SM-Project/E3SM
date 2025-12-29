#!/bin/bash
#==============================================================================
# Common Functions for Emulator Test Scripts
#==============================================================================

#------------------------------------------------------------------------------
# Colors (auto-detect terminal support)
# Use $'...' syntax for actual escape sequences
#------------------------------------------------------------------------------
if [[ -t 1 ]] && [[ "${TERM:-dumb}" != "dumb" ]]; then
    RED=$'\033[0;31m'
    GREEN=$'\033[0;32m'
    YELLOW=$'\033[1;33m'
    BLUE=$'\033[0;34m'
    BOLD=$'\033[1m'
    NC=$'\033[0m'
else
    # No color support
    RED=''
    GREEN=''
    YELLOW=''
    BLUE=''
    BOLD=''
    NC=''
fi

#------------------------------------------------------------------------------
# Logging Functions
#------------------------------------------------------------------------------
log_info() {
    if [[ "${QUIET}" != "true" ]]; then
        echo "${BLUE}[INFO]${NC} $1"
    fi
}

log_success() {
    echo "${GREEN}[OK]${NC} $1"
}

log_warning() {
    echo "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo "${RED}[ERROR]${NC} $1"
}

log_step() {
    echo "${BOLD}[$1]${NC} $2"
}

#------------------------------------------------------------------------------
# Usage
#------------------------------------------------------------------------------
show_usage() {
    cat << EOF
${BOLD}Emulator Components Test Tool${NC}

${BOLD}Usage:${NC} ./test [action] [options]

${BOLD}Actions:${NC}
    build       Configure and build (default)
    rebuild     Clean, then configure and build
    test        Run tests only (must be built first)
    configure   Run CMake configuration only
    clean       Remove build directory

${BOLD}Options:${NC}
    -h, --help              Show this help message
    -b, --backend BACKEND   Inference backend: stub (default), libtorch
    -m, --machine NAME      Specify machine (auto-detected)
    -j, --jobs N            Parallel jobs (default: 16)
    -t, --test-level N      Test level: 1=unit, 2=component, 3=all (default: 3)
    -q, --quiet             Minimal output (logs to files)
    -v, --verbose           Full output
    --mpi                   Enable MPI tests (requires compute node)
    --no-tests              Skip tests
    --debug                 Debug build (default: Release)
    --torch-dir PATH        LibTorch cmake directory (auto-detected on Perlmutter)
    --no-color              Disable colored output

${BOLD}Backends:${NC}
    stub        No-op backend for testing (no dependencies)
    libtorch    LibTorch C++ inference (auto-detects from pytorch module)

${BOLD}Examples:${NC}
    ./test                          # Build with stub backend
    ./test rebuild -b libtorch      # Build with LibTorch
    ./test rebuild --mpi            # Build with MPI tests enabled
    ./test test -t1                 # Run unit tests only

${BOLD}Logs:${NC} build/logs/
EOF
}

#------------------------------------------------------------------------------
# Machine Detection
#------------------------------------------------------------------------------
detect_machine() {
    if [[ -n "${NERSC_HOST}" ]]; then
        echo "${NERSC_HOST}"
    elif [[ -f /etc/ncar_system ]]; then
        echo "cheyenne"
    elif [[ -d /lcrc ]]; then
        echo "lcrc"
    else
        echo "generic"
    fi
}
