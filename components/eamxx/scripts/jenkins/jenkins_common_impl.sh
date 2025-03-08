#!/bin/bash -xe

# Check if the env var PR_LABELS is defined, and contains something meaningful
IFS=';' read -r -a labels <<< "$PR_LABELS";

# default values. By default, test only Stand-Alone (SA) on all machines
skip_testing=0
test_scripts=0
test_v0=0
test_v1=1
test_SA=1
skip_mappy=0
skip_weaver=0
skip_blake=0
is_at_run=0
if [ ${#labels[@]} -gt 0 ]; then
  # We do have some labels. Look for some supported ones.
  for label in "${labels[@]}"
  do
    if [ "$label" == "AT: skip eamxx-all" ]; then
      skip_testing=1
    elif [ "$label" == "AT: skip eamxx-sa" ]; then
      test_SA=0
    elif [ "$label" == "AT: skip eamxx-v1" ]; then
      test_v1=0
    elif [ "$label" == "scripts" ]; then
      test_scripts=1
    elif [ "$label" == "SCREAMv0" ]; then
      test_v0=1
    elif [ "$label" == "AT: Skip mappy" ]; then
      skip_mappy=1
    elif [ "$label" == "AT: Skip weaver" ]; then
      skip_weaver=1
    elif [ "$label" == "AT: Skip blake" ]; then
      skip_blake=1
    fi
  done
fi

if [ $skip_testing -eq 0 ]; then
  # User did not request to skip tests. Proceed with testing.

  v0_fail=0
  v1_fail=0
  sa_fail=0
  cov_fail=0
  scripts_fail=0
  memcheck_fail=0
  cd $JENKINS_SCRIPT_DIR/../..

  source scripts/jenkins/${NODE_NAME}_setup

  # Check machine-specific skipping
  if [[ $skip_mappy == 1 && "$SCREAM_MACHINE" == "mappy" ]]; then
    echo "Tests were skipped, since the Github label 'AT: Skip mappy' was found.\n"
    exit 0
  elif [[ $skip_weaver == 1 && "$SCREAM_MACHINE" == "weaver" ]]; then
    echo "Tests were skipped, since the Github label 'AT: Skip weaver' was found.\n"
    exit 0
  elif [[ $skip_blake == 1 && "$SCREAM_MACHINE" == "blake" ]]; then
    echo "Tests were skipped, since the Github label 'AT: Skip blake' was found.\n"
    exit 0
  fi

  if [[ "$(whoami)" == "e3sm-jenkins" ]]; then
      git config --local user.email "jenkins@ignore.com"
      git config --local user.name "Jenkins Jenkins"
  fi

  declare -i fails=0
  # The special string "AUTO" makes test-all-eamxx look for a baseline dir in the machine_specs.py file.
  # IF such dir is not found, then the default (ctest-build/baselines) is used
  BASELINES_DIR=AUTO

  TAS_ARGS="--baseline-dir $BASELINES_DIR \$compiler -p -c EKAT_DISABLE_TPL_WARNINGS=ON -m \$machine"
  # pm-gpu needs to do work in scratch area in order not to fill home quota
  if [[ "$SCREAM_MACHINE" == "pm-gpu" ]]; then
      TAS_ARGS="${TAS_ARGS} -w /pscratch/sd/e/e3smtest/e3sm_scratch/pm-gpu/ctest-build"
  fi

  # Now that we are starting to run things that we expect could fail, we
  # do not want the script to exit on any fail since this will prevent
  # later tests from running.
  set +e
  if [ -n "$PULLREQUESTNUM" ]; then
      is_at_run=1
  else
      TAS_ARGS="${TAS_ARGS} --test-level nightly"
      # We never want to submit a fake run to the dashboard
      if [ -z "$SCREAM_FAKE_ONLY" ]; then
          # Run EKAT tests for real nightly runs
          TAS_ARGS="${TAS_ARGS} --submit -c EKAT_ENABLE_TESTS=ON"
      fi
  fi

  SA_FAILURES_DETAILS=""
  # Run scream stand-alone tests (SA)
  if [ $test_SA -eq 1 ]; then
    this_output=$(./scripts/gather-all-data "./scripts/test-all-eamxx ${TAS_ARGS}" -l -m $SCREAM_MACHINE)
    if [[ $? != 0 ]]; then
      fails=$fails+1;
      sa_fail=1
      if [[ $is_at_run == 1 ]]; then
        errors=$(echo "$this_output" | grep -m1 -A 100000 'Build type ')
        SA_FAILURES_DETAILS+="$errors"
      fi
    fi

    # Add memcheck and coverage tests for nightlies on specific machines
    if [[ $is_at_run == 0 ]]; then
      if [[ "$SCREAM_MACHINE" == "mappy" ]]; then
        ./scripts/gather-all-data "./scripts/test-all-eamxx -t cov ${TAS_ARGS}" -l -m $SCREAM_MACHINE
        if [[ $? != 0 ]]; then
          fails=$fails+1;
          cov_fail=1
        fi
      fi

      # Add a memcheck test for mappy for nightlies
      if [[ "$SCREAM_MACHINE" == "mappy" ]]; then
        ./scripts/gather-all-data "./scripts/test-all-eamxx -t mem ${TAS_ARGS}" -l -m $SCREAM_MACHINE
        if [[ $? != 0 ]]; then
          fails=$fails+1;
          memcheck_fail=1
        fi
      fi

      if [[ -z "$SCREAM_FAKE_ONLY" && "$SCREAM_MACHINE" == "weaver" ]]; then
        # The fake-only tests don't launch any kernels which will cause all
        # the compute-sanitizer runs to fail.
        ./scripts/gather-all-data "./scripts/test-all-eamxx -t csm -t csr -t csi -t css ${TAS_ARGS}" -l -m $SCREAM_MACHINE
        if [[ $? != 0 ]]; then
          fails=$fails+1;
          memcheck_fail=1
        fi
      fi

    fi
  else
    echo "SCREAM Stand-Alone tests were skipped, since the Github label 'AT: skip eamxx-sa' was found.\n"
  fi

  # Run expensive tests on mappy only
  if [[ "$SCREAM_MACHINE" == "mappy" ]]; then

    # Run scripts-tests
    if [[ $test_scripts == 1 ]]; then
      ./scripts/scripts-tests -f -m $SCREAM_MACHINE
      if [[ $? != 0 ]]; then
        fails=$fails+1;
        scripts_fail=1
      fi

      ./scripts/cime-nml-tests
      if [[ $? != 0 ]]; then
        fails=$fails+1;
        scripts_fail=1
      fi
    fi

    # We do NOT want to do any of the items below if we are running
    # scripts-tests because the merge could change the state of the
    # repo and you never want to do that during testing. Also, CIME
    # V1 cases do not build with SCREAM_FAKE_ONLY on.
    #
    # Also, for the nightlies, we use a separate job to run CIME tests
    if [[ -z "$SCREAM_FAKE_ONLY" && $is_at_run == 1 ]]; then

      if [[ $test_v0 == 1 ]]; then
        ../../cime/scripts/create_test e3sm_scream_v0 -c -b master --wait
        if [[ $? != 0 ]]; then
          fails=$fails+1;
          v0_fail=1
        fi
      fi

      if [[ $test_v1 == 1 ]]; then
        # AT runs should be fast. => run only low resolution
        this_output=$(../../cime/scripts/create_test e3sm_scream_v1_at -c -b master --wait)
        if [[ $? != 0 ]]; then
          fails=$fails+1;
          v1_fail=1
          if [[ $is_at_run == 1 ]]; then
            errors=$(echo "$this_output" | grep -m1 -A 100000 'Waiting for tests to finish')
            V1_FAILURES_DETAILS+="$errors"
          fi
        fi
      else
          echo "SCREAM v1 tests were skipped, since the Github label 'AT: Skip v1 Testing' was found.\n"
      fi
    fi
  fi

  # Disable bash tracing to make the FAILS message stand out more
  set +x

  if [[ $fails > 0 ]]; then
    echo "######################################################"
    echo "FAILS DETECTED:"
    if [[ $sa_fail == 1 ]]; then
      echo "  SCREAM STANDALONE TESTING FAILED!"
      echo "$SA_FAILURES_DETAILS"
    fi
    if [[ $v1_fail == 1 ]]; then
      echo "  SCREAM V1 TESTING FAILED!"
      echo "$V1_FAILURES_DETAILS"
    fi
    if [[ $v0_fail == 1 ]]; then
      echo "  SCREAM V0 TESTING FAILED!"
    fi
    if [[ $memcheck_fail == 1 ]]; then
      echo "  SCREAM MEM CHECK TESTING FAILED!"
    fi
    if [[ $cov_fail == 1 ]]; then
      echo "  SCREAM COVERAGE BUILD FAILED!"
    fi
    if [[ $scripts_fail == 1 ]]; then
      echo "  SCREAM SCRIPTS TESTING FAILED!"
    fi
    echo "######################################################"
    exit 1
  fi

else
  echo "Tests were skipped, since the Github label 'AT: Integrate Without Testing' was found.\n"
fi
