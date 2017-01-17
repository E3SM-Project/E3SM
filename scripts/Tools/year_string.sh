#!/bin/bash

# Convert an integer representing a year into a year string in the
# CESM format: yyyy (e.g., 0001 rather than 1)
#
# It is permitted for the provided integer to be greater than 9999; in
# this case, a string with greater than 4 digits is output

# Written by Bill Sacks, 8-17-12

# ----------------------------------------------------------------------
# LOCAL FUNCTIONS DEFINED HERE
# ----------------------------------------------------------------------

function Usage {
    progname=`basename $0`
    echo "Usage: $progname year"
    echo ""
}

# ----------------------------------------------------------------------
# BEGIN MAIN SCRIPT
# ----------------------------------------------------------------------

if [ $# -ne 1 ]; then
    Usage
    exit 1
fi

year=$1

if [ $year -lt 10 ]; then
    year_string="000$year"
elif [ $year -lt 100 ]; then
    year_string="00$year"
elif [ $year -lt 1000 ]; then 
    year_string="0$year"
else
    year_string="$year"
fi

echo $year_string