#!/bin/sh

echo `date`
githash=`git rev-parse --short HEAD`
gitdate=`git log -1 --date=format:"%Y_%m_%d" --format="%ad"`
testname="$gitdate-$githash"
echo "Running test for hash: $githash, date: $gitdate"
./create_test \
  SMS_P48x1.ne4_ne4.FSCREAM-SA \
  --machine quartz \
  --baseline-root /p/lustre2/donahue5/E3SM_simulations/SCREAM/developer_tests/test_root/ \
  -r /p/lustre2/donahue5/E3SM_simulations/SCREAM/developer_tests/test_root/ \
  --output-root /p/lustre2/donahue5/E3SM_simulations/SCREAM/developer_tests/test_root/ \
  --project cbronze \
  -q pdebug \
  -t $testname
#  ERS_P48x1.ne4_ne4.FSCREAM-LR.quartz_intel.cam-force_smp \
#  ERS_P96x1.ne4_ne4.FSCREAM-LR.quartz_intel.cam-force_smp \
#  ERS_P24x2.ne4_ne4.FSCREAM-LR.quartz_intel.cam-force_smp \
#  ERS_P48x2.ne4_ne4.FSCREAM-LR.quartz_intel.cam-force_smp \

#  ERS_P36x2.ne4_ne4.FC5AV1C-L \
#  ERS_P48x1.ne4_ne4.FC5AV1C-L \
#  ERS_P96x1.ne4_ne4.FC5AV1C-L \
#  ERS_P24x2.ne4_ne4.FSCREAM-LR.quartz_intel.cam-force_smp \
#  ERS_P48x2.ne4_ne4.FSCREAM-LR.quartz_intel.cam-force_smp \
