#!/bin/bash

# Mariana Vertenstein
# 11-07-2014


#======================================================================
# Local functions
#======================================================================

function Usage {
    echo "SYNOPSIS"
    echo "     $progname [options]"
    echo ""
    echo "     Does a baseline comparison of history files and/or a baseline generation of history files"
    echo "     If there are multiple history file types for the component (e.g., h0 & h1),"
    echo "     it uses one of each type (e.g., the latest h0 file and the latest h1 file)."
    echo "     It then outputs the test status (PASS/BFAIL/GFAIL/etc.) for each comparison / generation."
    echo ""
    echo "OPTIONS"
    echo "     -testname <casename>  case name (identical to $CASE)"
    echo ""
    echo "     -test_dir <path>      Path to directory containing test run directories (optional)"
    echo "                           A given test's run directory can be found in: \$test_dir/\$CASE/run"
    echo ""
    echo "     -help                 Print this help message and exit"
    echo ""
    echo ""
}

# Given a status_and_info string, return just the status portion
# Inputs:
#   - status_and_info: colon-delimited string of the format
#     "STATUS:Extra info about failure". The string before the colon
#     is the test status. The string after the colon is extra
#     information abaout a failure (this may be blank, but the colon
#     must be present). (It is okay for the extra info to itself
#     contain one or more colons.)
function get_status {
    local status_and_info="$1"
    echo $status_and_info | cut -d ':' -f 1
}

# Given a status_and_info string, return just the info portion
# Inputs:
#   - status_and_info: colon-delimited string of the format
#     "STATUS:Extra info about failure". The string before the colon
#     is the test status. The string after the colon is extra
#     information abaout a failure (this may be blank, but the colon
#     must be present). (It is okay for the extra info to itself
#     contain one or more colons.)
function get_info {
    local status_and_info="$1"
    echo $status_and_info | cut -d ':' -f 2-
}

# Prints the status of the test in a standardized format for test results
# Inputs:
#   - status
#   - info
#   - testcase: name of test case
function print_status {
    local status="$1"
    local info="$2"
    local testcase="$3"

    # Enclose info in parentheses
    if [ -n "$info" ]; then
	info_str="($info)"
    else
	info_str=""
    fi    

    # Print formatted test result
    printf '%-3s %s %s\n' "$status" "${testcase}" "$info_str"
}

#======================================================================
# Begin main script
#======================================================================

tools_dir=`dirname $0`

#----------------------------------------------------------------------
# Define default values for command-line arguments
#----------------------------------------------------------------------
baseline_dir=''
test_dir=''
generate=''
compare_tag=''
msg='' 
 
#----------------------------------------------------------------------
# Process command-line arguments
#----------------------------------------------------------------------
while [ $# -gt 0 ]; do
    case $1 in
	-msg )
	    msg=$2
	    shift
	    ;;
	-baseline_dir )
	    baseline_dir=$2
	    shift
	    ;;
	-test_dir )
	    test_dir=$2
	    shift
	    ;;
	-testcase_base )
	    testcase_base=$2
	    shift
	    ;;
	-generate_tag )
	    generate_tag=$2
	    shift
	    ;;
	-compare_tag )
	    compare_tag=$2
	    shift
	    ;;

	-help )
	    Usage
	    exit 0
	    ;;
	* )
	    echo "Unknown argument: $1" >&2
	    echo "Run $progname -help for usage" >&2
	    exit 1
	    ;;
    esac
    shift
done

#----------------------------------------------------------------------
# Exit if required command-line arguments weren't provided 
#----------------------------------------------------------------------
error=0  # no errors yet
if [ -z "$baseline_dir" ]; then
    echo "$progname: baseline_dir must be provided" >&2
    error=1
fi
if [ -z "$test_dir" ]; then
    echo "$progname: test_dir must be provided" >&2
    error=1
fi
if [ $error -gt 0 ]; then
    echo "" >&2
    echo "Run $progname -help for usage" >&2
    exit 1
fi

#------------------------------------------------------------------
# Determine location of cprnc
#------------------------------------------------------------------
if [ -z $CASEROOT ]; then
    echo " environment variable CASEROOT is not defined "
    exit 1
fi
cd $CASEROOT
cprnc_exe=`./xmlquery CCSM_CPRNC -value`

#------------------------------------------------------------------
# Loop over models
#------------------------------------------------------------------

overall_compare_status='PASS'
overall_generate_status='PASS'

models=( cam cice cism clm2 cpl pop )
for model in ${models[*]}; do
    if [ "$model" = "cism" ]; then
	extensions=(h)
    elif [ "$model" = "clm2" ]; then
	extensions=(h0 h1 h2 h3 h4 h5)
    elif [ "$model" = "cpl" ]; then
	extensions=(hi)
    elif [ "$model" = "cice" ]; then
	extensions=(h)
    elif [ "$model" = "cam" ]; then
	extensions=(h0 h1 h2 h3 h4 h5 h6 h7)
    elif [ "$model" = "pop" ]; then
	extensions=(h)
    fi

    #------------------------------------------------------------------
    # Loop over history file extensions
    #------------------------------------------------------------------
    for extension in ${extensions[*]}; do
    
        # Set baseline history name 
	# Note that this name drops (1) the timestamp, and (2) the
        # instance number for multi-instance runs
	baseline_hist=${model}.${extension}.nc

        #--------------------------------------------------------------
        # Find last component hist files in this run directory, and
        # determine corresponding name of the baseline file (used for
        # either generation or comparison)
        #--------------------------------------------------------------
        # Note that we find the last alphabetically rather than by time
        # stamp, because the last by time stamp can be non-deterministic.
        # Note that we need a * after ${model} to capture multi-instance
        # output

	test_hist=`cd $test_dir; ls -1 *.${model}*.${extension}.*.nc.base 2>/dev/null | tail -1`

	if [ -n "$test_hist" ]; then

            #-------------------------------------------------------------
            # Do comparison, if desired
            #-------------------------------------------------------------
	    if [ -n "$compare_tag" ]; then
		compare_result=`${tools_dir}/component_compare.sh  -baseline_dir "$baseline_dir" -baseline_hist "$baseline_hist" -test_dir "$test_dir" -test_hist "$test_hist" -cprnc_exe "$cprnc_exe"`
		compare_status=`get_status "$compare_result"`
		compare_info=`get_info "$compare_result"`
		if [ "$compare_status" != "BFAIL_NA" ]; then
		    if [ -z "$msg" ]; then 
			print_status "$compare_status" "$compare_info" "${testcase_base}.$model.${extension}.nc : baseline compare ${model}.${extension} base with ${compare_tag}"
		    else
			print_status "$compare_status" "$compare_info" "${testcase_base}.$model.${extension}.nc : baseline compare ${model}.${extension} ($msg)"
		    fi
		    # if ANY teststatus is FAIL then overall status is FAIL
		    if [ "$compare_status" = "FAIL" ]; then
			overall_compare_status="FAIL"
		    fi
		    # if ANY teststatus is BFAIL then overall status is BFAIL
		    if [ "$overall_compare_status" != "FAIL" ]; then
			if [ "$compare_status" = "BFAIL" ]; then
			    overall_compare_status="BFAIL"
			fi
		    fi
		fi
	    fi

            #--------------------------------------------------------------
            # Do baseline generation, if desired
            #--------------------------------------------------------------
	    if [ -n "$generate_tag" ]; then
		generate_result=`${tools_dir}/component_generate.sh -baseline_dir "$baseline_dir" -baseline_hist "$baseline_hist" -test_dir "$test_dir" -test_hist "$test_hist"`
		generate_status=`get_status "$generate_result"`
		generate_info=`get_info "$generate_result"`
		if [ "$generate_status" != "GFAIL_NA" ]; then
		    print_status "$generate_status" "$generate_info" "${testcase_base}.$model.${extension}.nc : baseline generate ${model}.${extension} in baseline dir ${generate_tag}"
		fi
		# if ANY teststatus is GFAIL then overall status is GFAIL
		if [ "$generate_status" = "GFAIL" ]; then
		    overall_compare_status="GFAIL"
		fi
	    fi
	fi

    done  # loop over history file extensions

done  # loop over models

if [ -n "$compare_tag" ]; then
    if [ -z "$msg" ]; then 
	print_status "$overall_compare_status" "$compare_info" "${testcase_base} : baseline compare summary"
    else
	print_status "$overall_compare_status" "$compare_info" "${testcase_base} : baseline compare summary ($msg)"
    fi
fi
if [ -n "$generate_tag" ]; then
    print_status "$overall_generate_status" "$generate_info" "${testcase_base} : baseline generate summary"
fi


