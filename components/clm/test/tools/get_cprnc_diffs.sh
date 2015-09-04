#!/bin/bash

# This script extracts lines from the output of cprnc that tell us
# which variables differ between two files
#
# Usage: get_cprnc_diffs filename

# ----------------------------------------------------------------------
# SET PARAMETERS HERE
# ----------------------------------------------------------------------

# maximum number of differences to extract from the cprnc output
maxdiffs=200

# ----------------------------------------------------------------------
# LOCAL FUNCTIONS DEFINED HERE
# ----------------------------------------------------------------------

# This function gets differences for one prefix (e.g., "RMS")
# Usage: get_diffs prefix
# (also uses $infile and $maxdiffs from the parent script)
function get_diffs {
    prefix=$1
    outfile=${infile}.${prefix}.$$
    grep "$prefix" $infile > $outfile
    numlines=`wc -l $outfile | awk '{print $1}'`
    if [ $numlines -gt $maxdiffs ]; then
	echo "WARNING: Too many instances of $prefix - only printing last $maxdiffs"
	tail -$maxdiffs $outfile
    else
	cat $outfile
    fi
    rm $outfile
}

# ----------------------------------------------------------------------
# BEGIN MAIN SCRIPT
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# Handle command-line arguments
# ----------------------------------------------------------------------

if [[ $# -ne 1 ]]; then
    echo "Usage: get_cprnc_diffs filename"
    exit 1
fi

infile=$1

# ----------------------------------------------------------------------
# Do the processing
# ----------------------------------------------------------------------

get_diffs RMS
get_diffs FILLDIFF
