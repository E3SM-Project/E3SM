#!/bin/bash
#
# Reusable argument parsing for test scripts
#
# Usage: source this file, then call parse_test_args "$@"
# After calling, the following variables will be set:
#   CLEAN_ONLY, BUILD_ONLY, VERBOSE, ENABLE_PYTHON
#

# Default values
CLEAN_ONLY=false
BUILD_ONLY=false
VERBOSE=false
ENABLE_PYTHON=false

parse_test_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --clean-only)
                CLEAN_ONLY=true
                shift
                ;;
            --python)
                ENABLE_PYTHON=true
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
            --help|-h)
                echo "Usage: $0 [OPTIONS]"
                echo ""
                echo "Options:"
                echo "  --clean-only  Remove build directory and exit"
                echo "  --build-only  Build without running tests"
                echo "  --python      Build the embedded-Python inference backend"
                echo "  --verbose     Show verbose test output"
                echo "  --help        Show this help message"
                exit 0
                ;;
            *)
                echo "ERROR: Unknown option: $1" >&2
                exit 1
                ;;
        esac
    done
}
