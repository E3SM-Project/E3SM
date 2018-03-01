#!/bin/bash

# This file contains functions used by other bash scripts in this directory.

# This function echoes the command given by $1 (cmd), then executes it.
# However, if $2 (dryrun) is non-zero, then it only does the echo, not the execution.
# Usage: do_cmd cmd dryrun
# Returns 0 on success, non-zero on failure; if there is an error, the error string is echoed.
function do_cmd {
    if [[ $# -ne 2 ]]; then
	echo "ERROR in do_cmd: wrong number of arguments: expected 2, received $#"
	exit 1
    fi

    local cmd=$1
    local dryrun=$2

    echo $cmd
    if [ $dryrun -eq 0 ]; then
	# We use 'eval $cmd' rather than just '$cmd', because the
	# latter doesn't work right if the command contains any quoted
	# strings (e.g., svn ci -m "this is my message")
	eval $cmd
	if [ $? -ne 0 ]; then
	    echo "ERROR in do_cmd: error executing command"
	    exit 2
	fi
    fi

    return 0
}

# make sure that the given file name argument was provided, and that
# the file exists; exit the script with a usage message if either of
# these is not true
# 
# Usage: check_file_arg filename_arg description
# (description is echoed if there is an error)
# Example: check_file_arg "$input_file" "input"
# (note that $input_file must be in quotes)
function check_file_arg {
    local filename=$1
    local description=$2
    
    if [ -z "$filename" ]; then
	echo "ERROR: Must specify $description file"
	Usage
	exit 1
    fi

    if [ ! -f $filename ]; then
	echo "ERROR: Can't find $description file: $filename"
	Usage
	exit 1
    fi
}

