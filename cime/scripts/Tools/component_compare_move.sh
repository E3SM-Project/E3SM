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
    echo "     Moves component history file to new name with suffix in the testcase directory"
    echo "     If there are multiple history file types for the component (e.g., h0 & h1),"
    echo "     it uses one of each type (e.g., the latest h0 file and the latest h1 file)."
    echo ""
    echo "OPTIONS"
    echo "     -rundir <path>        Path to directory containing test run directories (optional)"
    echo "                           A given test's run directory can be found in: \$rundir/\$CASE/run"
    echo ""
    echo "     -testcase <casename>  case name (identical to $CASE)"
    echo ""
    echo "     -suffix <name>        Suffix name to attache to testcase"
    echo ""
    echo "     -help                 Print this help message and exit"
    echo ""
    echo "EXAMPLES"
    echo "     $progname"
    echo "       -rundir /glade/scratch/\$USER -testcase ERS.f19_g16.FC5.yellowstone_intel -suffix init"
    echo ""
}

#======================================================================
# Begin main script
#======================================================================

tools_dir=`dirname $0`

#----------------------------------------------------------------------
# Define default values for command-line arguments
#----------------------------------------------------------------------
testcase=''
rundir=''
suffix=''
add_iop=''
 
#----------------------------------------------------------------------
# Process command-line arguments
#----------------------------------------------------------------------
while [ $# -gt 0 ]; do
    case $1 in
	-rundir )
	    rundir=$2
	    shift
	    ;;
	-testcase )
	    testcase=$2
	    shift
	    ;;
	-suffix )
	    suffix=$2
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
if [ -z "$suffix" ]; then
    echo "$progname: suffix must be provided" >&2
    error=1
fi
if [ $error -gt 0 ]; then
    echo "" >&2
    echo "Run $progname -help for usage" >&2
    exit 1
fi

#------------------------------------------------------------------
# Loop over models
#------------------------------------------------------------------

if [ -n "$add_iop" ]; then
    echo " add_iop is set"
    suffix="${suffix}_${add_iop}"
fi

cd $rundir
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
	extensions=(h0 h1 h2 h3 h4 h5 h6 h7)
    elif [ "$model" = "pop" ]; then
	extensions=(h)
    fi

    #------------------------------------------------------------------
    # Loop over history file extensions
    #------------------------------------------------------------------

    for extension in ${extensions[*]}; do
    
        #--------------------------------------------------------------
        # Find last component hist file in this run directory, and
        # determine corresponding name of the baseline file (used for
        # either generation or comparison)
        #--------------------------------------------------------------
        # Note that we find the last alphabetically rather than by time
        # stamp, because the last by time stamp can be non-deterministic.
        # Note that we need a * after ${model} to capture multi-instance
        # output

	if [ "$suffix" = "multiinst" ]; then
	    instances=(_0001 _0002)
	else
	    instances=(none _0001 _0002)
	fi
	for inst in ${instances[*]}; do 
	    if [ "$inst" = "none" ]; then
		test_hist=`ls -1 ${testcase}.${model}.${extension}.*.nc 2>/dev/null | tail -1`
	    else
		test_hist=`ls -1 ${testcase}.${model}${inst}.${extension}.*.nc 2>/dev/null | tail -1`
	    fi

	    if [[ -f ${test_hist} ]]; then 
		test_hist_suffix=${test_hist}.${suffix}
		cd "$rundir";

		# remove suffix file if already present
		if [[ -f ${test_hist_suffix} ]]; then
		    rm "$test_hist_suffix"
		fi

		# FIXME(bja, 2014-11) temp change to test if cp fixes the
		# problem with mfilt in restart tests. If it does, then
		# rename this script to component_compare_copy.sh !
		/bin/cp "$test_hist" "$test_hist_suffix"
	    fi
	done

    done

done  # loop over history file extensions



