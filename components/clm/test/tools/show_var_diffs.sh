#!/bin/bash

# This script processes a log file that was output by test_driver,
# giving lists of all variables with differences in values (those with
# RMS errors), and all variables with differences in fill patterns.
#
# This assumes that the log file contains output like:
#  RMS foo
#  RMS bar
#  FILLDIFF foo
#  FILLDIFF bar
# Some characteristics of these output lines are:
# - they begin with a leading space, followed by RMS or FILLDIFF
# - the variable name is in the second column of the line
#
# Note that (as of 4-5-12) the log file only contains output from the
# last file that didn't match, so this could potentially miss
# something -- especially if there are both h0 and h1 files in the
# comparison.

# Usage: show_var_diffs logfile

# ----------------------------------------------------------------------
# LOCAL FUNCTIONS DEFINED HERE
# ----------------------------------------------------------------------

# This function shows the differences for one prefix (e.g., "RMS")
# Usage: show_diffs prefix
# (also uses $logfile from the parent script)
# 
# Matches lines that start with the regular expression "^ ${prefix}"
# (note that one leading space is expected before the prefix)
#
# Assumes that the variable name is in the second column of matching lines
function show_diffs {
    prefix=$1

    # first determine if there were warnings relating to this prefix
    grep "WARNING: Too many instances of ${prefix}" $logfile > /dev/null
    if [ $? -eq 0 ]; then  # found a warning
	echo "WARNING: Some output was truncated; this may not be a complete list"
    fi
    
    # now make a list of all variables matching this prefix
    grep "^ ${prefix}" $logfile > $logfile.tmp.$$
    if [ $? -eq 0 ]; then
	awk '{print $2}' $logfile.tmp.$$ | sort | uniq
    else
	echo "(no differences)"
    fi

    rm $logfile.tmp.$$
}

# ----------------------------------------------------------------------
# BEGIN MAIN SCRIPT
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# Handle command-line arguments
# ----------------------------------------------------------------------

if [[ $# -ne 1 ]]; then
    echo "Usage: show_var_diffs logfile"
    exit 1
fi

logfile=$1

# ----------------------------------------------------------------------
# Do the processing
# ----------------------------------------------------------------------

echo "Variables with differences in values:"
show_diffs "RMS"

echo ""
echo "Variables with differences in fill patterns:"
show_diffs "FILLDIFF"