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
    echo "     Does a baseline comparison of history files or a baseline generation of history files for one test"
    echo "     If there are multiple history file types for the component (e.g., h0 & h1),"
    echo "     it uses one of each type (e.g., the latest h0 file and the latest h1 file)."
    echo "     It then outputs the test status (PASS/BFAIL/GFAIL/etc.) for each comparison / generation."
    echo ""
    echo "     Exactly one of -generate_tag or -compare_tag should be provided."
    echo "     (They cannot both be provided, because only one -baseline_dir is supplied.)"
    echo ""
    echo "OPTIONS"
    echo "     -baseline_dir <path>  Path to directory containing baselines for this particular test (required)"
    echo "                           If -compare_tag is given, this gives the path to the"
    echo "                           baselines to compare against."
    echo "                           If -generate_tag is given, this gives the path in which"
    echo "                           the new baselines should be placed."
    echo ""
    echo "     -test_dir <path>      Path to the given test's run directory (required)"
    echo ""
    echo "     -testcase <name>      Full name of case including testid (required)"
    echo ""
    echo "     -testcase_base <name> Name of test case without testid, used for printing results (required)"
    echo ""
    echo "     -generate_tag <tag>   Tag to use for baseline generation (optional)"
    echo ""
    echo "     -compare_tag <tag>    Tag to use for baseline comparison (optional)"
    echo ""
    echo "     -cprnc_exe <name>     Full pathname to cprnc executable (optional)"
    echo "                           This is optional, but if not provided, then the environment"
    echo "                           variable CASEROOT must be defined (in which case, CCSM_CPRNC"
    echo "                           will be found from the xml files in CASEROOT)."
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
    local action="$4"

    # Enclose info in parentheses
    if [ -n "$info" ]; then
	info_str="($info)"
    else
	info_str=""
    fi    

    # Print formatted test result
    printf '%-3s %s %s\n' "$status" "${testcase}" "${action}"

    echo "COMMENT $info"
}

function print_comment {
    local status="$1"
    local info="$2"
    local testcase="$3"

    # Enclose info in parentheses
    if [ -n "$info" ]; then
	info_str="($info)"
    else
	info_str=""
    fi

    echo "COMMENT for $testcase $info"
}

#======================================================================
# Begin main script
#======================================================================

progname=`basename $0`
tools_dir=`dirname $0`

#----------------------------------------------------------------------
# Define default values for command-line arguments
#----------------------------------------------------------------------
cprnc_exe=''
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
	-testcase )
	    testcase=$2
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
	-cprnc_exe )
	    cprnc_exe=$2
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
if [ -z "$testcase" ]; then
    echo "$progname: testcase must be provided" >&2
    error=1
fi
if [ -z "$testcase_base" ]; then
    echo "$progname: testcase_base must be provided" >&2
    error=1
fi
if [ -z "$compare_tag" -a -z "$generate_tag" ]; then
    echo "$progname: Must provide either -compare_tag or -generate_tag" >&2
    error=1
fi
if [ -n "$compare_tag" -a -n "$generate_tag" ]; then
    echo "$progname: Can only provide -compare_tag OR -generate_tag, not both" >&2
    echo "(this is because only one -baseline_dir is provided to this script)" >&2
    error=1
fi

if [ $error -gt 0 ]; then
    echo "" >&2
    echo "Run $progname -help for usage" >&2
    exit 1
fi

#------------------------------------------------------------------
# Determine location of cprnc, if not provided
#------------------------------------------------------------------
if [ -z "$cprnc_exe" ]; then
    if [ -z $CASEROOT ]; then
	echo " environment variable CASEROOT is not defined "
	exit 1
    fi
    cd $CASEROOT
    cprnc_exe=`./xmlquery CCSM_CPRNC -value`
fi

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

	test_hist=`cd $test_dir; ls -1 ${testcase}.${model}*.${extension}.*.nc.base 2>/dev/null | tail -1`

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
			print_comment "$compare_status" "$compare_info : $model.${extension}.nc : baseline compare ${model}.${extension} base with ${compare_tag}" "${testcase_base}"
		    else
			print_comment "$compare_status" "$compare_info : $model.${extension}.nc : baseline compare ${model}.${extension} ($msg)" "${testcase_base}"
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
		    print_comment "$generate_status" "$generate_info : $model.${extension}.nc : baseline generate ${model}.${extension} in baseline dir ${generate_tag}" "${testcase_base}"
		fi
		# if ANY teststatus is GFAIL then overall status is GFAIL
		if [ "$generate_status" = "GFAIL" ]; then
		    overall_generate_status="GFAIL"
		fi
	    fi
	fi

    done  # loop over history file extensions

done  # loop over models

if [ -n "$compare_tag" ]; then
    if [ -z "$msg" ]; then 
	print_status "$overall_compare_status" "$compare_info : baseline compare summary" "${testcase_base}" compare
    else
	print_status "$overall_compare_status" "$compare_info : baseline compare summary ($msg)" "${testcase_base}" compare
    fi

    if [[ $overall_compare_status == PASS ]]; then
        exit 0
    else
        exit 2
    fi
fi
if [ -n "$generate_tag" ]; then
    print_status "$overall_generate_status" "$generate_info : baseline generate summary" "${testcase_base}" generate

    if [[ $overall_generate_status == PASS ]]; then
        exit 0
    else
        exit 2
    fi
fi


