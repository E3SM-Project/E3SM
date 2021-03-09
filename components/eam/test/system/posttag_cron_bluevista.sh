#!/bin/sh 
# script to test cam on bluevista - triggered from bluesky nightly cron job

if [ $# -ne 2 ]; then
    echo "td_nightly_bluevista.sh: incorrect number of input args" 
    exit 1
fi

tag_to_grab=$1
tag_for_baseline=$2

collections=/fis/cgd/cseg/csm/models/atm/cam

#prepare results directory for new tag
results_dir=${collections}/test_results/${tag_to_grab}_bluevista
if [ -d $results_dir ]; then
    rm -rf $results_dir
fi
mkdir $results_dir
cd $results_dir

env BL_ROOT=${collections}/${tag_for_baseline} \
    CAM_ROOT=${collections}/${tag_to_grab} \
    CAM_INPUT_TESTS=${collections}/${tag_to_grab}/models/atm/cam/test/system/tests_posttag_bluevista \
    ${collections}/${tag_to_grab}/models/atm/cam/test/system/test_driver.sh -f

exit 0
