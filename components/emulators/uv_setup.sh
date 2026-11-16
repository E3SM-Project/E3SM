#!/bin/bash
#
# Create a Python virtual environment with uv for building and testing the
# emulator package (torch, numpy, PyYAML).
#
# Optional: CMake can also use an existing torch install. See README.md.
#

set -e

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/scripts/shell/colors.sh"

venv_dir="${script_dir}/.venv"
python_version="3.12"
torch_backend="auto"
# Pinned to the tested version. 2.7 is the first release whose Linux wheels all
# use the new C++ ABI (_GLIBCXX_USE_CXX11_ABI=1); older wheels make the package
# fail to compile.
torch_spec="torch==2.7.1"
force=false

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --backend <name>   Torch build: auto (default), cpu, cu126, cu128, ...
                     'auto' lets uv pick from the installed GPU driver.
  --torch <spec>     Torch requirement (default: '${torch_spec}')
  --python <ver>     Python version (default: ${python_version})
  --venv <dir>       Environment directory (default: ${venv_dir})
  --force            Recreate the environment if it already exists
  -h, --help         Show this help
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --backend) torch_backend="$2"; shift 2 ;;
        --torch)   torch_spec="$2"; shift 2 ;;
        --python)  python_version="$2"; shift 2 ;;
        --venv)    venv_dir="$2"; shift 2 ;;
        --force)   force=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) print_error "Unknown option: $1"; usage; exit 1 ;;
    esac
done

if ! command -v uv >/dev/null 2>&1; then
    print_error "uv is not installed. See https://docs.astral.sh/uv/getting-started/installation/"
    exit 1
fi

if [ -d "${venv_dir}" ]; then
    if [ "$force" = true ]; then
        print_status "Removing ${venv_dir}"
        rm -rf "${venv_dir}"
    else
        print_info "${venv_dir} exists; updating it (use --force to recreate)"
    fi
fi

if [ ! -d "${venv_dir}" ]; then
    print_status "Creating ${venv_dir} with Python ${python_version}"
    uv venv --prompt emu --python "${python_version}" "${venv_dir}"
fi

# --torch-backend picks the matching PyTorch index (cpu, cuXXX, rocmX.Y) for
# torch only, so numpy and PyYAML still come from PyPI.
print_status "Installing ${torch_spec} (backend: ${torch_backend}), numpy, PyYAML"
VIRTUAL_ENV="${venv_dir}" uv pip install \
    --torch-backend "${torch_backend}" \
    "${torch_spec}" numpy pyyaml

"${venv_dir}/bin/python" - <<'EOF'
import torch
print(f"torch {torch.__version__}, CUDA available: {torch.cuda.is_available()}, "
      f"C++11 ABI: {torch._C._GLIBCXX_USE_CXX11_ABI}")
import numpy, yaml
EOF

print_status "Done. To build and test with this environment:"
echo "    source ${venv_dir}/bin/activate"
echo "    ./test"
