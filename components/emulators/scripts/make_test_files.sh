#!/bin/bash
#
# Populate a (build-tree) fixture directory for the emulator tests:
#   * copies the hand-written fixtures (configs/*.yaml, models/*.py)
#   * generates the coupling configs/plans
#   * with --libtorch, writes the TorchScript (.pt) fixtures
#
# Everything lands flat in OUTDIR; nothing is written to the source tree.
#
# Usage: make_test_files.sh --outdir DIR [--python PY] [--libtorch]

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
src_root="$(cd "${script_dir}/.." && pwd)"
fixture_src="${src_root}/tests/fixtures"

source "${script_dir}/shell/colors.sh"

libtorch=false
outdir=""
python="python3"

while [ $# -gt 0 ]; do
    case "$1" in
        --libtorch) libtorch=true; shift ;;
        --outdir)   outdir="$2"; shift 2 ;;
        --python)   python="$2"; shift 2 ;;
        *) print_error "Unknown option: $1"; exit 1 ;;
    esac
done

if [ -z "${outdir}" ]; then
    print_error "--outdir is required"
    exit 1
fi

mkdir -p "${outdir}"
outdir="$(cd "${outdir}" && pwd)"

# Hand-written fixtures
cp "${fixture_src}"/configs/*.yaml "${outdir}/"
cp "${fixture_src}"/models/*.py "${outdir}/"

# Generated configs / coupling plans
print_status "Writing emulator test configs to ${outdir}"
"${python}" "${script_dir}/config/create_coupling_test_files.py" "${outdir}"

if [ "${libtorch}" = true ]; then
    print_status "Writing LibTorch test fixtures to ${outdir}"
    "${python}" "${fixture_src}/models/make_libtorch_fixture.py" "${outdir}"
    "${python}" "${fixture_src}/models/make_coupling_fixtures.py" "${outdir}"
fi
