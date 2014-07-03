#!/bin/bash

# Bill Sacks
# 8-22-12

#======================================================================
# Overview
#======================================================================

# This script does the '-generate' step for a single test case and a
# single history file.
#
# It is intended to be used by another script, such as
# component_gen_comp.
#
# The definition of the result is similar to that in testcase_end, but
# differs slightly. For example: In testcase_end, it is a FAIL for
# baselineroot to not exist; here we allow that, and happily generate
# that directory for you if it doesn't already exist.
# 
# Exit status will generally be 0 (even for test failure), but will be
# non-zero for incorrect usage.
#
# See usage message for details on inputs & return value.

#======================================================================
# Testing
#======================================================================

# There currently is no unit test script for this. Here are the unit
# tests that should be done on this script:
#
# - Missing arguments
#   - return status should be "UNDEF"
#
# - test_hist=''
#   - return status should be "BFAIL"
#   - example: tst=`component_generate.sh -baseline_dir /glade/scratch/sacks/cesm_baselines/test_script/ERI44y.f09_g16.TGRCP85.bluefire_ibm -baseline_hist cism.h.nc -test_dir /ptmp/sacks/ERI44y.f09_g16.TGRCP85.bluefire_ibm.C.114029/run -test_hist ''`
# 
# - success
#   - return status should be "PASS"
#   - example: tst=`component_generate.sh -baseline_dir /glade/scratch/sacks/cesm_baselines/test_script/ERI44y.f09_g16.TGRCP85.bluefire_ibm -baseline_hist cism.h.nc -test_dir /ptmp/sacks/ERI44y.f09_g16.TGRCP85.bluefire_ibm.C.114029/run -test_hist ERI44y.f09_g16.TGRCP85.bluefire_ibm.C.114029.cism.h.2046-01-01-00000.nc`

#======================================================================
# Local functions
#======================================================================

function Usage {
    echo "SYNOPSIS"
    echo "     $progname [options]"
    echo ""
    echo "RETURN VALUE"
    echo "     String of the format 'STATUS:Extra info about failure'"
    echo "     The string before the colon is the test status."
    echo "     The string after the colon is extra information about a failure"
    echo "     (this may be blank, but the colon will still be present)"
    echo ""
    echo "     Possible values for STATUS are:"
    echo "     UNDEF : undefined result; this includes incorrect usage"
    echo "     BFAIL1: no history file in test case"
    echo "       (note that BFAIL1 is only a problem if it occurs for all possible"
    echo "       history extensions for a given component)"
    echo "     BFAIL2: error creating baseline directory or copying baseline file into place"
    echo "     PASS  : success"
    echo ""
    echo "OPTIONS"
    echo "     -baseline_dir <name>   Full path to the baseline directory for this test (required)"
    echo "     -baseline_hist <name>  Name used for history file in the baseline directory (required)"
    echo "     -test_dir <name>       Full path to the directory containing history files for this test (required)"
    echo "     -test_hist <name>      Name of history file in the test directory (required)"
    echo "                            (it is NOT an error for this to be an empty string, but this will generate a BFAIL1;"
    echo "                            this will be the case when there are no history files for a component in the test directory)"
    echo "     -help                  Print this help message and exit"
    echo ""
    echo "EXAMPLES"
    echo "     TEST=SMS.T31_g37.IG4804.bluefire_ibm"
    echo "     CASE=\$TEST.GC.134426"
    echo ""
    echo "     $progname -baseline_dir /glade/scratch/\$USER/cesm_baselines/cesm1_1_beta16/\$TEST"
    echo "       -baseline_hist clm.h.nc"
    echo "       -test_dir /glade/scratch/\$USER/\$CASE/run"
    echo "       -test_hist \$CASE.clm.h0.0001-12.nc"
    echo "     This will copy /glade/scratch/\$USER/\$CASE/run/\$CASE.clm.h0.0001-12.nc"
    echo "     to /glade/scratch/\$USER/cesm_baselines/cesm1_1_beta16/\$TEST/clm.h.nc"
    echo ""
    echo "     $progname -baseline_dir /glade/scratch/\$USER/cesm_baselines/cesm1_1_beta16/\$TEST"
    echo "       -baseline_hist clm.h.nc"
    echo "       -test_dir /glade/scratch/\$USER/\$CASE/run"
    echo "       -test_hist ''"
    echo "     This will generate a BFAIL1 because there is no history file to copy to the baseline directory"
    echo ""
}

# Prints the test status and an info string, separated by ':'
# Inputs:
#   status: test status
#   info: optional auxiliary info about test failure
function print_result {
    status="$1"
    info="$2"

    echo "${status}:${info}"
}


#======================================================================
# Begin main script
#======================================================================

progname=`basename $0`

#----------------------------------------------------------------------
# Set default return values
#----------------------------------------------------------------------
status='UNDEF'
info=''

#----------------------------------------------------------------------
# Define default values for command-line arguments
#----------------------------------------------------------------------
baseline_dir=''
baseline_hist=''
test_dir=''

# test_hist will often be the empty string, so we initialize it to
# something else so we can check whether it has been provided
test_hist_init='XXX_NOTGIVEN_XXX'
test_hist=$test_hist_init

#----------------------------------------------------------------------
# Process command-line arguments
#----------------------------------------------------------------------
while [ $# -gt 0 ]; do
    case $1 in
	-baseline_dir )
	    baseline_dir=$2
	    shift
	    ;;
	-baseline_hist )
	    baseline_hist=$2
	    shift
	    ;;
	-test_dir )
	    test_dir=$2
	    shift
	    ;;
	-test_hist )
	    test_hist=$2
	    shift
	    ;;
	-help )
	    Usage
	    exit 0
	    ;;
	* )
	    echo "$progname: Unknown argument: $1" >&2
	    echo "Run $progname -help for usage" >&2
            # return default values for status & info
	    print_result $status "$info"
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
if [ -z "$baseline_hist" ]; then
    echo "$progname: baseline_hist must be provided" >&2
    error=1
fi
if [ -z "$test_dir" ]; then
    echo "$progname: test_dir must be provided" >&2
    error=1
fi
if [ "$test_hist" = "$test_hist_init" ]; then
    echo "$progname: test_hist must be provided" >&2
    error=1
fi

if [ $error -gt 0 ]; then
    echo "" >&2
    echo "Run $progname -help for usage" >&2
    # return default values for status & info
    print_result $status "$info"
    exit 1
fi

#----------------------------------------------------------------------
# Make sure there is a history file in the test case
#----------------------------------------------------------------------
if [ -z "$test_hist" ]; then
    status="BFAIL1"
    info="no history file in test case"
    print_result $status "$info"
    exit 0
fi

#----------------------------------------------------------------------
# Create baseline directory, if necessary
#----------------------------------------------------------------------
mkdir -p $baseline_dir
if [ $? -ne 0 ]; then
    status="BFAIL2"
    info="error creating baseline directory"
    print_result $status "$info"
    exit 0
fi
chmod ug+w,a+rx $baseline_dir
chmod ug+w,a+rx $baseline_dir/..

#----------------------------------------------------------------------
# Copy history file to baseline directory
#----------------------------------------------------------------------
cp $test_dir/$test_hist $baseline_dir/$baseline_hist
if [ $? -ne 0 ]; then
    status="BFAIL2"
    info="error copying history file to baseline directory"
    print_result $status "$info"
    exit 0
fi
chmod ug+w,a+r $baseline_dir/$baseline_hist

#----------------------------------------------------------------------
# Return final test status
#----------------------------------------------------------------------
status="PASS"
print_result $status "$info"
exit 0

