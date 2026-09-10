#! /bin/bash

set -euo pipefail

# Configure if needed; this is a no-op once build/ has a CMakeCache.txt.
cmake -S . -B build -DDEXPR_ENABLE_TESTS=ON -DDEXPR_ENABLE_TOOL=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

cmake --build build --parallel
ctest --test-dir build --output-on-failure --parallel
