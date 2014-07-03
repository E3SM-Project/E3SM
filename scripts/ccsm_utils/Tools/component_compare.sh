#!/bin/bash

# Bill Sacks
# 8-23-12

#======================================================================
# Overview
#======================================================================

# This script does the '-compare' step for a single test case and a
# single history file.
#
# It is intended to be used by another script, such as
# component_gen_comp.
#
# The definition of the result is similar to the compare_hist status
# in testcase_end, but may differ slightly. Other than possible minor
# differences, though, one use of this script can be: You start
# baseline runs and test runs at the same time (both using
# 'create_test_suite ... -clean off'); test runs finish before
# baseline generation finished (leading to BFAILs); then you can run
# this script with model=cpl to do the comparisons after the baseline
# generation completes (but first you will need to copy the cpl.hi.nc
# files in each baseline directory to cpl.h.nc).
# 
# Exit status will generally be 0 (even for test failure), but will be
# non-zero for incorrect usage.
#
# cprnc should be in your path, but if it isn't, the script tries to
# find cprnc in its yellowstone location
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
# - test_hist='', no baseline
#   - return status should be "BFAIL_NA"
#   - example: tst=`component_compare.sh -baseline_dir /glade/scratch/sacks/cesm_baselines/test_script/ERI44y.f09_g16.TGRCP85.bluefire_ibm -baseline_hist DOES_NOT_EXIST -test_dir /ptmp/sacks/ERI44y.f09_g16.TGRCP85.bluefire_ibm.C.114029/run -test_hist ''`
# 
# - Given test history file doesn't exist, no baseline
#   - return status should be "BFAIL_NA"
#   - example: tst=`component_compare.sh -baseline_dir /glade/scratch/sacks/cesm_baselines/test_script/ERI44y.f09_g16.TGRCP85.bluefire_ibm -baseline_hist DOES_NOT_EXIST -test_dir /ptmp/sacks/ERI44y.f09_g16.TGRCP85.bluefire_ibm.C.114029/run -test_hist DOES_NOT_EXIST`
#
# - test_hist='' (but baseline file exists)
#   - return status should be "FAIL"
#   - example: tst=`component_compare.sh -baseline_dir /glade/scratch/sacks/cesm_baselines/test_script/ERI44y.f09_g16.TGRCP85.bluefire_ibm -baseline_hist cism.h.nc -test_dir /ptmp/sacks/ERI44y.f09_g16.TGRCP85.bluefire_ibm.C.114029/run -test_hist ''`
#
# - Given test history file doesn't exist (but baseline file exists)
#   - return status should be "FAIL"
#   - example: tst=`component_compare.sh -baseline_dir /glade/scratch/sacks/cesm_baselines/test_script/ERI44y.f09_g16.TGRCP85.bluefire_ibm -baseline_hist cism.h.nc -test_dir /ptmp/sacks/ERI44y.f09_g16.TGRCP85.bluefire_ibm.C.114029/run -test_hist DOES_NOT_EXIST`
#
# - Given test history file exists, but no baseline
#   - return status should be "BFAIL"
#   - example: component_compare.sh -baseline_dir $baseline_dir -baseline_hist DOES_NOT_EXIST -test_dir $test_dir -test_hist ERS_Ld211_P192x1.f19_g16.I_2000_CNDV_CROP.yellowstone_intel_crop.G.061555.clm2.h0.1998-07-30-00000.nc
#
# - History comparison fails: difference in values (e.g., set up by
# manually modifying values in a file)
#   - return status should be "FAIL"
#   - example (after running the PASS example from component_generate.sh, then running: ncap2 -s 'vvel[$time,$level,$y0,$x0]=1.7' ERI44y.f09_g16.TGRCP85.bluefire_ibm.C.114029.cism.h.2046-01-01-00000.nc test.nc): tst=`component_compare.sh -baseline_dir /glade/scratch/sacks/cesm_baselines/test_script/ERI44y.f09_g16.TGRCP85.bluefire_ibm -baseline_hist cism.h.nc -test_dir /ptmp/sacks/ERI44y.f09_g16.TGRCP85.bluefire_ibm.C.114029/run -test_hist test.nc`
#
# - History comparison fails: different time stamp on file (e.g., from
# a different case from the one used to generate the baseline)
# (this is really testing the behavior of cprnc more than this script
# per se)
#   - return status should be "FAIL"
#   - example (after running the PASS example from component_generate.sh): tst=`component_compare.sh -baseline_dir /glade/scratch/sacks/cesm_baselines/test_script/ERI44y.f09_g16.TGRCP85.bluefire_ibm -baseline_hist cism.h.nc -test_dir /ptmp/sacks/SMS_D.f09_g16.TG.bluefire_ibm.C.114029/run -test_hist SMS_D.f09_g16.TG.bluefire_ibm.C.114029.cism.h.0011-01-01-00000.nc`
#
# - History comparison passes
#   - return status should be "PASS"
#   - example (after running the PASS example from component_generate.sh): tst=`component_compare.sh -baseline_dir /glade/scratch/sacks/cesm_baselines/test_script/ERI44y.f09_g16.TGRCP85.bluefire_ibm -baseline_hist cism.h.nc -test_dir /ptmp/sacks/ERI44y.f09_g16.TGRCP85.bluefire_ibm.C.114029/run -test_hist ERI44y.f09_g16.TGRCP85.bluefire_ibm.C.114029.cism.h.2046-01-01-00000.nc`

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
    echo "     UNDEF   : undefined result; this includes incorrect usage"
    echo "     BFAIL_NA: no baseline file and no history file for the test case"
    echo "       (note that BFAIL_NA is only a problem if it occurs for all possible"
    echo "       history extensions for a given component)"
    echo "     BFAIL   : no baseline file, but there IS a history file for the test case"
    echo "     FAIL    : comparison fails (including: baseline exists but no history file for the test case)"
    echo "     PASS    : success"
    echo ""
    echo "OPTIONS"
    echo "     -baseline_dir <name>   Full path to the baseline directory for this test (required)"
    echo "     -baseline_hist <name>  Name used for history file in the baseline directory (required)"
    echo "     -test_dir <name>       Full path to the directory containing history files for this test (required)"
    echo "                            NOTE: You must have write permission in this directory (cprnc.out file is put there)"
    echo "     -test_hist <name>      Name of history file in the test directory (required)"
    echo "                            (it is NOT an error for this to be an empty string, but this may generate a FAIL;"
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
    echo "     This will compare /glade/scratch/\$USER/\$CASE/run/\$CASE.clm.h0.0001-12.nc"
    echo "     with /glade/scratch/\$USER/cesm_baselines/cesm1_1_beta16/\$TEST/clm.h.nc"
    echo ""
    echo "     $progname -baseline_dir /glade/scratch/\$USER/cesm_baselines/cesm1_1_beta16/\$TEST"
    echo "       -baseline_hist clm.h.nc"
    echo "       -test_dir /glade/scratch/\$USER/\$CASE/run"
    echo "       -test_hist ''"
    echo "     This will generate a FAIL because there is no history file to compare with the baseline"
    echo "     (or a BFAIL_NA if /glade/scratch/\$USER/cesm_baselines/cesm1_1_beta16/\$TEST/clm.h.nc also doesn't exist)"
    echo ""
}

# Given a relative path, convert to an absolute path
# (from http://www.lancejian.com/2011/04/13/get-absolute-path-of-the-running-bash-script.html)
function absolute_path {
    relative_path="$1"
    absolute_path=`cd "$relative_path"; pwd`
    echo "$absolute_path"
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
# Set parameters
#======================================================================
# If cprnc isn't in your path, this is where to find it
cprnc_default=/glade/p/cesm/cseg/tools/cprnc/cprnc

#======================================================================
# Begin main script
#======================================================================

progname=`basename $0`

# need absolute path (rather than relative path) because we use this
# path after we have cd'ed to another location
tools_dir=$(absolute_path `dirname $0`)  

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
# Determine path to cprnc
#----------------------------------------------------------------------
cprnc_path=`command -v cprnc`
if [ $? -gt 0 ]; then
    # cprnc not found in path; use default
    cprnc_path=$cprnc_default
fi

#----------------------------------------------------------------------
# Determine whether baseline file exists
#----------------------------------------------------------------------
if [ -f $baseline_dir/$baseline_hist ]; then
    baseline_exists=1
else
    baseline_exists=0
fi

#----------------------------------------------------------------------
# Determine whether test history file exists
#----------------------------------------------------------------------
# Test history file may be the empty string if there is known to be no
# history file for this component. Another possibility is that the
# given test history file doesn't exist; this is more likely to be an
# error in the usage of the script rather than a real test failure, so
# we may want to rethink the status value in that case.
if [ -z "$test_hist" ]; then
    test_exists=0
elif [ -f $test_dir/$test_hist ]; then
    test_exists=1
else
    test_exists=0
fi

#----------------------------------------------------------------------
# If baseline and/or test history file don't exist, return appropriate
# error status
#----------------------------------------------------------------------

if [[ $baseline_exists -eq 0 && $test_exists -eq 0 ]]; then
    status="BFAIL_NA"
    info="neither baseline nor test history file exists"
    print_result $status "$info"
    exit 0
elif [[ $baseline_exists -eq 0 && $test_exists -eq 1 ]]; then
    status="BFAIL"
    info="baseline history file does not exist"
    print_result $status "$info"
    exit 0
elif [[ $baseline_exists -eq 1 && $test_exists -eq 0 ]]; then
    status="FAIL"
    info="no history file in test case"
    print_result $status "$info"
    exit 0
fi

# Note: at this point, we know that there is both a baseline history
# file and a test history file

#----------------------------------------------------------------------
# Compare history files, get test status
# Put output in a file named ${test_dir}/${test_hist}.cprnc.out
#----------------------------------------------------------------------
# We cd to test_dir so that cprnc.out is put there (note that this
# assumes that the user has write permission in test_dir)
curdir=`pwd`
cd $test_dir

status=`CCSM_CPRNC=$cprnc_path $tools_dir/hist_compare.csh $test_dir/$test_hist $baseline_dir/$baseline_hist | tail -1`
mv cprnc.out ${test_hist}.cprnc.out

cd $curdir

print_result $status "$info"
exit 0
