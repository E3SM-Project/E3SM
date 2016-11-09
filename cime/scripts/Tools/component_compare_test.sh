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
    echo "     Compares two component history files in the testcase directory"
    echo "     If there are multiple history file types for the component (e.g., h0 & h1),"
    echo "     it uses one of each type (e.g., the latest h0 file and the latest h1 file)."
    echo "     It then outputs the test status (PASS/FAIL/etc.) for each comparison / generation."
    echo ""
    echo "OPTIONS"
    echo "     -rundir <path>        Path to directory containing test run directories (optional)"
    echo "                           A given test's run directory can be found in: \$rundir/\$CASE/run"
    echo ""
    echo "     -testcase <name>      Full name of testcase including testid"
    echo ""
    echo "     -testcase_base <name> Full name of testcase without testid"
    echo ""
    echo "     -suffix1 <name>       Suffix to attach to first file to compare"
    echo ""
    echo "     -suffix1 <name>       Suffix to attach to second file to compare"
    echo ""
    echo "     -help                 Print this help message and exit"
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
    printf '%-3s %s compare\n' "$status" "${testcase}"

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

tools_dir=`dirname $0`

#----------------------------------------------------------------------
# Define default values for command-line arguments
#----------------------------------------------------------------------
rundir=''
testcase=''
testcase_base=''
suffix1=''
suffix2=''
add_iop=''
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
	-rundir )
	    rundir=$2
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
	-suffix1 )
	    suffix1=$2
	    shift
	    ;;
	-suffix2 )
	    suffix2=$2
	    shift
	    ;;
	-add_iop )
	    add_iop=$2
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
if [ -z "$rundir" ]; then
    echo "$progname: rundir must be provided" >&2
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
if [ -z "$suffix1" ]; then
    echo "$progname: suffix1 must be provided" >&2
    error=1
fi
if [ -z "$suffix2" ]; then
    echo "$progname: suffix2 must be provided" >&2
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

cd $rundir
overall_status='PASS'

models=(cam cice clm2 pop cism cpl)
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
	extensions=(h0 h1 h2 h3 h4 h5 h6 h7 h8 h9 hs)
    elif [ "$model" = "pop" ]; then
	extensions=(h)
        if [ "$suffix2" = "rest" ]; then
            # Skip restart checks for pop! Temporary hack until MPAS goes in!
            continue
        fi
    fi

    #------------------------------------------------------------------
    # Loop over history file extensions
    #------------------------------------------------------------------
    for extension in ${extensions[*]}; do
    
        #--------------------------------------------------------------
        # Find last component hist files in this run directory, and
        # determine corresponding name of the baseline file (used for
        # either generation or comparison)
        #--------------------------------------------------------------
        # Note that we find the last alphabetically rather than by time
        # stamp, because the last by time stamp can be non-deterministic.
        # Note that we need a * after ${model} to capture multi-instance
        # output
        #-------------------------------------------------------------
        # Do comparison
        #-------------------------------------------------------------

	if [[ "$testcase" =~ ^NC[KRO]_.* ]]; then

	    hist1=`cd $rundir; ls -1 ${testcase}.${model}.${extension}.*.nc.${suffix1} 2>/dev/null | tail -1`
	    hist2_0001=`cd $rundir; ls -1 ${testcase}.${model}_0001.${extension}.*.nc.${suffix2} 2>/dev/null | tail -1`
	    hist2_0002=`cd $rundir; ls -1 ${testcase}.${model}_0002.${extension}.*.nc.${suffix2} 2>/dev/null | tail -1`

	    if  [[ -f ${hist1} ]] && [[ -f ${hist2_0001} ]] && [[ -f ${hist2_0002} ]] ; then

		if [ "$model" != "cpl" ]; then
		    # do all model comparisons except for cpl, since cpl history does not write out all instances - but just instance 1

		    compare_result=`${tools_dir}/component_compare.sh -baseline_dir "$rundir" -test_dir "$rundir" -baseline_hist "$hist1" -test_hist "$hist2_0001" -cprnc_exe "$cprnc_exe"`
		    compare_status=`get_status "$compare_result"`
		    compare_info=`get_info "$compare_result"`
		    print_comment "$compare_status" "$compare_info: ${model}.${extension}.nc : test compare ${model}.${extension} (.${suffix1} and .${suffix2} for _0001)" "${testcase_base}"
		    if [ "$compare_status" != "PASS" ]; then
			overall_status="FAIL"
		    fi

		    compare_result=`${tools_dir}/component_compare.sh -baseline_dir "$rundir" -test_dir "$rundir" -baseline_hist "$hist1" -test_hist "$hist2_0002" -cprnc_exe "$cprnc_exe"`
		    compare_status=`get_status "$compare_result"`
		    compare_info=`get_info "$compare_result"`
		    print_comment "$compare_status" "$compare_info: ${model}.${extension}.nc : test compare ${model}.${extension} (.${suffix1} and .${suffix2} for _0002)" "${testcase_base}"
		    if [ "$compare_status" != "PASS" ]; then
			overall_status="FAIL"
		    fi

		fi
	    fi		

	elif [[ "$testcase" =~ .*_N2.* ]]; then

	    hist1_0001=`cd $rundir; ls -1 ${testcase}.${model}_0001.${extension}.*.nc.${suffix1} 2>/dev/null | tail -1`
	    hist2_0001=`cd $rundir; ls -1 ${testcase}.${model}_0001.${extension}.*.nc.${suffix2} 2>/dev/null | tail -1`
	    hist1_0002=`cd $rundir; ls -1 ${testcase}.${model}_0002.${extension}.*.nc.${suffix1} 2>/dev/null | tail -1`
	    hist2_0002=`cd $rundir; ls -1 ${testcase}.${model}_0002.${extension}.*.nc.${suffix2} 2>/dev/null | tail -1`

	    if  [[ -f ${hist1_0001} ]] && [[ -f ${hist1_0002} ]] && [[ -f ${hist2_0001} ]] && [[ -f ${hist2_0002} ]] ; then

		compare_result=`${tools_dir}/component_compare.sh -baseline_dir "$rundir" -test_dir "$rundir" -baseline_hist "$hist1_0001" -test_hist "$hist1_0002" -cprnc_exe "$cprnc_exe"`
		compare_status=`get_status "$compare_result"`
		compare_info=`get_info "$compare_result"`
		print_comment "$compare_status" "$compare_info: ${model}.${extension}.nc : test compare ${model}.${extension} (.${suffix1} for _0001 and .${suffix2} for _0001)" "${testcase_base}"
		if [ "$compare_status" != "PASS" ]; then
		    overall_status="FAIL"
		fi

		compare_result=`${tools_dir}/component_compare.sh -baseline_dir "$rundir" -test_dir "$rundir" -baseline_hist "$hist2_0002" -test_hist "$hist2_0002" -cprnc_exe "$cprnc_exe"`
		compare_status=`get_status "$compare_result"`
		compare_info=`get_info "$compare_result"`
		print_comment "$compare_status" "$compare_info: ${model}.${extension}.nc : test compare ${model}.${extension} (.${suffix1} for _0002 and .${suffix2} for _0002)" "${testcase_base}"
		if [ "$compare_status" != "PASS" ]; then
		    overall_status="FAIL"
		fi

	    fi


	elif [[ "$testcase" =~ .*_IOP.* ]]; then

	    if [ -z "$add_iop" ]; then

		# determine status of non-IOP test

		hist1=`cd $rundir; ls -1 ${testcase}.${model}.${extension}.*.nc.${suffix1} 2>/dev/null | tail -1`
		hist2=`cd $rundir; ls -1 ${testcase}.${model}.${extension}.*.nc.${suffix2} 2>/dev/null | tail -1`
		
		if  [[ -f ${hist1} ]] && [[ -f ${hist2} ]] ; then
		    compare_result=`${tools_dir}/component_compare.sh -baseline_dir "$rundir" -test_dir "$rundir" -baseline_hist "$hist1" -test_hist "$hist2" -cprnc_exe "$cprnc_exe"`
		    compare_status=`get_status "$compare_result"`
		    compare_info=`get_info "$compare_result"`
		    print_comment "$compare_status" "$compare_info: ${model}.${extension}.nc : test compare ${model}.${extension} (.${suffix1} and .${suffix2} files)" "${testcase_base}"
		    if [ "$compare_status" != "PASS" ]; then
			overall_status="FAIL"
		    fi
		fi
	    fi

	    if [ -n "$add_iop" ]; then

		# determine status of iop test

		hist1=`cd $rundir; ls -1 ${testcase}.${model}*.${extension}.*.nc.${suffix1}            2>/dev/null | tail -1`
		hist2=`cd $rundir; ls -1 ${testcase}.${model}*.${extension}.*.nc.${suffix1}_${add_iop} 2>/dev/null | tail -1`

		if  [[ -f ${hist1} ]] && [[ -f ${hist2} ]] ; then
		    compare_result=`${tools_dir}/component_compare.sh  -baseline_dir "$rundir" -test_dir "$rundir" -baseline_hist "$hist1" -test_hist "$hist2" -cprnc_exe "$cprnc_exe"`
		    compare_status=`get_status "$compare_result"`
		    if [ "$compare_status" != "PASS" ]; then
			overall_status="FAIL"
		    fi
		    compare_info=`get_info "$compare_result"`
		    print_comment "$compare_status" "$compare_info: ${model}.${extension}.nc : test compare ${model}.${extension} (.${suffix1} and .${suffix1}_${add_iop} files)" "${testcase_base}"
		fi
		
		hist1=`cd $rundir; ls -1 ${testcase}.${model}*.${extension}.*.nc.${suffix2}            2>/dev/null | tail -1`
		hist2=`cd $rundir; ls -1 ${testcase}.${model}*.${extension}.*.nc.${suffix2}_${add_iop} 2>/dev/null | tail -1`

		if  [[ -f ${hist1} ]] && [[ -f ${hist2} ]] ; then
		    compare_result=`${tools_dir}/component_compare.sh -baseline_dir "$rundir" -test_dir "$rundir" -baseline_hist "$hist1" -test_hist "$hist2" -cprnc_exe "$cprnc_exe"`
		    compare_status=`get_status "$compare_result"`
		    if [ "$compare_status" != "PASS" ]; then
			overall_status="FAIL"
		    fi
		    compare_info=`get_info "$compare_result"`
		    print_comment "$compare_status" "$compare_info: ${model}.${extension}.nc : test compare ${model}.${extension} (.${suffix2} and .${suffix2}_${add_iop} files)" "${testcase_base}"
		fi

	    fi

	else

	    # No _IOP_ or _N2_ attributes or multi-instance NCK_ or NCR_ tests
	    hist1=`cd $rundir; ls -1 ${testcase}.${model}.${extension}.*.nc.${suffix1} 2>/dev/null | tail -1`
	    hist2=`cd $rundir; ls -1 ${testcase}.${model}.${extension}.*.nc.${suffix2} 2>/dev/null | tail -1`
	    if  [[ -f ${hist1} ]] && [[ -f ${hist2} ]] ; then
		compare_result=`${tools_dir}/component_compare.sh -baseline_dir "$rundir" -test_dir "$rundir" -baseline_hist "$hist1" -test_hist "$hist2" -cprnc_exe "$cprnc_exe"`
		compare_status=`get_status "$compare_result"`
		compare_info=`get_info "$compare_result"`
		print_comment "$compare_status" "$compare_info: ${model}.${extension}.nc : test compare ${model}.${extension} (.${suffix1} and .${suffix2} files)" "${testcase_base}"
		if [ "$compare_status" != "PASS" ]; then
		    overall_status="FAIL"
		fi
	    fi
	    
	fi

    done # loop over history file extensions

done  # loop over models

if [ -z "$msg" ]; then 
    print_status "$overall_status" "$compare_info: test functionality summary" "${testcase_base}"
else
    print_status "$overall_status" "$compare_info: test functionality summary ($msg)" "${testcase_base}"
fi


