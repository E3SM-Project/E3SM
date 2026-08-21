#! /bin/bash

set -euo pipefail

# Configure if needed; this is a no-op once build/ has a CMakeCache.txt.
cmake -S . -B build

cmake --build build --parallel
ctest --test-dir build --output-on-failure --parallel
