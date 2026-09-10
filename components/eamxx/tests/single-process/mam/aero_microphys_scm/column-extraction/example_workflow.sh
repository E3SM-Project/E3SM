#!/bin/bash
#
# Example workflow for extracting columns and creating single-column tests
#
# This script demonstrates how to:
# 1. Extract a specific column from EAMxx output
# 2. Use it as initial conditions for a single-column test
#

set -e  # Exit on error

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXTRACT_SCRIPT="${SCRIPT_DIR}/extract_single_column.py"

# Default values
COLUMN_INDEX=0
TIME_INDEX=-1  # -1 means last time step
VERBOSE=false

# Usage function
usage() {
    cat << EOF
Usage: $0 [OPTIONS] INPUT_FILE OUTPUT_FILE

Extract a single column from EAMxx output for use in single-column tests.

Options:
  -c, --col-idx INDEX    Column index to extract (default: 0)
  -t, --time-idx INDEX   Time index to extract (default: last, -1)
  -v, --verbose          Enable verbose output
  -h, --help             Show this help message

Examples:
  # Extract column 0 from last time step
  $0 model_output.nc single_col_ic.nc

  # Extract column 50 from time step 5
  $0 -c 50 -t 5 model_output.nc single_col_ic.nc

  # Extract with verbose output
  $0 -v -c 100 model_output.nc single_col_ic.nc

EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--col-idx)
            COLUMN_INDEX="$2"
            shift 2
            ;;
        -t|--time-idx)
            TIME_INDEX="$2"
            shift 2
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -*)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
        *)
            if [ -z "$INPUT_FILE" ]; then
                INPUT_FILE="$1"
            elif [ -z "$OUTPUT_FILE" ]; then
                OUTPUT_FILE="$1"
            else
                echo "Too many arguments"
                usage
                exit 1
            fi
            shift
            ;;
    esac
done

# Check required arguments
if [ -z "$INPUT_FILE" ] || [ -z "$OUTPUT_FILE" ]; then
    echo "Error: INPUT_FILE and OUTPUT_FILE are required"
    usage
    exit 1
fi

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file not found: $INPUT_FILE"
    exit 1
fi

# Check if Python script exists
if [ ! -f "$EXTRACT_SCRIPT" ]; then
    echo "Error: Extraction script not found: $EXTRACT_SCRIPT"
    exit 1
fi

# Check for Python and required packages
if ! command -v python3 &> /dev/null; then
    echo "Error: python3 not found"
    exit 1
fi

# Build the command
CMD="python3 $EXTRACT_SCRIPT $INPUT_FILE $OUTPUT_FILE --col-idx $COLUMN_INDEX"

if [ "$TIME_INDEX" != "-1" ]; then
    CMD="$CMD --time-idx $TIME_INDEX"
fi

if [ "$VERBOSE" = true ]; then
    CMD="$CMD --verbose"
fi

# Show what we're doing
echo "Extracting column from EAMxx output..."
echo "  Input:  $INPUT_FILE"
echo "  Output: $OUTPUT_FILE"
echo "  Column: $COLUMN_INDEX"
if [ "$TIME_INDEX" = "-1" ]; then
    echo "  Time:   last time step"
else
    echo "  Time:   index $TIME_INDEX"
fi
echo ""

# Run the extraction
$CMD

# Check if successful
if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Successfully extracted column to: $OUTPUT_FILE"
    echo ""
    echo "Next steps:"
    echo "  1. Copy the file to your test data directory"
    echo "  2. Update your test's input.yaml to reference this file"
    echo "  3. Set number_of_global_columns: 1 in the grid configuration"
    echo ""
    echo "Example for MAM4 aerosol test:"
    echo "  cp $OUTPUT_FILE \$SCREAM_DATA_DIR/init/"
    echo "  # Then update input.yaml initial_conditions:filename"
else
    echo "Error: Extraction failed"
    exit 1
fi
