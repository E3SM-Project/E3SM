#!/bin/bash
#
# Reusable argument parsing for test scripts
#
# Usage: source this file, then call parse_test_args "$@"
# After calling, the following variables will be set:
#   CLEAN_ONLY, BUILD_ONLY, VERBOSE, CMAKE_ARGS, TEST_REGEX
#

# Default values
CLEAN_ONLY=false
BUILD_ONLY=false
VERBOSE=false
CMAKE_ARGS=()
TEST_REGEX=""

parse_test_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --clean-only)
                CLEAN_ONLY=true
                shift
                ;;
            --build-only)
                BUILD_ONLY=true
                shift
                ;;
            --verbose|-v)
                VERBOSE=true
                shift
                ;;
            -D*)
                CMAKE_ARGS+=("$1")
                shift
                ;;
            --help|-h)
                echo "Usage: $0 [OPTIONS] [TEST_REGEX]"
                echo ""
                echo "Runs all tests, or those with names matching TEST_REGEX."
                echo ""
                echo "Options:"
                echo "  --clean-only  Remove build directory and exit"
                echo "  --build-only  Build without running tests"
                echo "  --verbose     Show verbose test output"
                echo "  -D<var>=<val> Pass an option to CMake, e.g."
                echo "                -DEMULATOR_ENABLE_LIBTORCH=ON"
                echo "  --help        Show this help message"
                exit 0
                ;;
            -*)
                echo "ERROR: Unknown option: $1" >&2
                exit 1
                ;;
            *)
                if [[ -n "$TEST_REGEX" ]]; then
                    echo "ERROR: Only one test regex allowed: $1" >&2
                    exit 1
                fi
                TEST_REGEX="$1"
                shift
                ;;
        esac
    done
}
