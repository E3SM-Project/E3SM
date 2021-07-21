#!/bin/bash -xe

# Check if the env var PR_LABELS is defined, and contains something meaningful
IFS=';' read -r -a labels <<< "$PR_LABELS";

# default values
skip_testing=0
test_scripts=0
test_cime=0
if [ ${#labels[@]} -gt 0 ]; then
  # We do have some labels. Look for some supported ones.
  for label in "${labels[@]}"
  do
    if [ "$label" == "AT: Integrate Without Testing" ]; then
      skip_testing=1
    elif [ "$label" == "scripts" ]; then
      test_scripts=1
    elif [ "$label" == "CIME" ]; then
      test_cime=1
    fi
  done
fi

if [ $skip_testing -eq 0 ]; then
  # User did not request to skip tests. Proceed with testing.

  cd $JENKINS_SCRIPT_DIR/../..

  source scripts/jenkins/${NODE_NAME}_setup

  if [[ "$(whoami)" == "e3sm-jenkins" ]]; then
      git config --local user.email "jenkins@ignore.com"
      git config --local user.name "Jenkins Jenkins"
  fi

  SUBMIT="--submit"
  AUTOTESTER_CMAKE=""
  # If this is a nightly run, we do NOT want this script to stop on the first failed command
  # since that will prevent subsequent dashboard submissions
  set +e
  if [ -n "$PULLREQUESTNUM" ]; then
      SUBMIT="" # We don't submit AT runs
      AUTOTESTER_CMAKE="-c SCREAM_AUTOTESTER=ON"
      set -e # This is an AT run, not nightly, so it's OK to stop on first fail
  fi

  # The special string "AUTO" makes test-all-scream look for a baseline dir in the machine_specs.py file.
  # IF such dir is not found, then the default (ctest-build/baselines) is used
  BASELINES_DIR=AUTO

  ./scripts/gather-all-data "./scripts/test-all-scream --baseline-dir $BASELINES_DIR \$compiler -c EKAT_DISABLE_TPL_WARNINGS=ON $AUTOTESTER_CMAKE -p -i -m \$machine $SUBMIT" -l -m $SCREAM_MACHINE

  # Add a valgrind and coverage tests for mappy for nightlies
  if [[ -n "$SUBMIT" && "$SCREAM_MACHINE" == "mappy" ]]; then
    ./scripts/gather-all-data "./scripts/test-all-scream -t valg -t cov --baseline-dir $BASELINES_DIR \$compiler -c EKAT_DISABLE_TPL_WARNINGS=ON -p -i -m \$machine $SUBMIT" -l -m $SCREAM_MACHINE
  fi

  # scripts-tests is pretty expensive, so we limit this testing to mappy
  if [[ $test_scripts == 1 && "$SCREAM_MACHINE" == "mappy" ]]; then
    # JGF: I'm not sure there's much value in these dry-run comparisons
    # since we aren't changing HEADs
    ./scripts/scripts-tests -g
    ./scripts/scripts-tests -c

    ./scripts/scripts-tests -f -m $SCREAM_MACHINE
  fi

  # Run SCREAM CIME suite
  if [[ $test_cime == 1 && "$SCREAM_MACHINE" == "mappy" ]]; then
    ../../cime/scripts/create_test e3sm_scream -c -b master
  fi

else
  echo "Tests were skipped, since the Github label 'AT: Integrate Without Testing' was found.\n"
fi
