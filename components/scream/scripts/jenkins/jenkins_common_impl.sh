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
    if [ "$label" == "AT: Integrate Without Testing" ]; then
      skip_testing=1
    elif [ "$label" == "AT: Skip Stand-Alone Testing" ]; then
      test_SA=0
    elif [ "$label" == "AT: Skip v1 Testing" ]; then
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
  # The special string "AUTO" makes test-all-scream look for a baseline dir in the machine_specs.py file.
  # IF such dir is not found, then the default (ctest-build/baselines) is used
  BASELINES_DIR=AUTO

  TAS_ARGS="--baseline-dir $BASELINES_DIR \$compiler -c EKAT_DISABLE_TPL_WARNINGS=ON ${TAS_ARGS} -p -i -m \$machine"
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
          TAS_ARGS="${TAS_ARGS} --submit"
      fi
  fi

  # Run scream stand-alone tests (SA)
  if [ $test_SA -eq 1 ]; then
    ./scripts/gather-all-data "./scripts/test-all-scream ${TAS_ARGS}" -l -m $SCREAM_MACHINE
    if [[ $? != 0 ]]; then fails=$fails+1; fi

    # Add a valgrind and coverage tests for mappy for nightlies
    if [[ $is_at_run == 0 && "$SCREAM_MACHINE" == "mappy" ]]; then
      ./scripts/gather-all-data "./scripts/test-all-scream -t dbg --mem-check ${TAS_ARGS}" -l -m $SCREAM_MACHINE
      if [[ $? != 0 ]]; then fails=$fails+1; fi

      ./scripts/gather-all-data "./scripts/test-all-scream -t cov ${TAS_ARGS}" -l -m $SCREAM_MACHINE
      if [[ $? != 0 ]]; then fails=$fails+1; fi
    fi

    # Add a cuda-memcheck test for weaver for nightlies
    if [[ $is_at_run == 0 && "$SCREAM_MACHINE" == "weaver" ]]; then
      ./scripts/gather-all-data "./scripts/test-all-scream -t dbg --mem-check ${TAS_ARGS}" -l -m $SCREAM_MACHINE
      if [[ $? != 0 ]]; then fails=$fails+1; fi
    fi
  else
    echo "SCREAM Stand-Alone tests were skipped, since the Github label 'AT: Skip Stand-Alone Testing' was found.\n"
  fi

  # Run expensive tests on mappy only
  if [[ "$SCREAM_MACHINE" == "mappy" ]]; then

    # Run scripts-tests
    if [[ $test_scripts == 1 ]]; then
      # JGF: I'm not sure there's much value in these dry-run comparisons
      # since we aren't changing HEADs
      ./scripts/scripts-tests -g -m $SCREAM_MACHINE
      if [[ $? != 0 ]]; then fails=$fails+1; fi
      ./scripts/scripts-tests -c -m $SCREAM_MACHINE
      if [[ $? != 0 ]]; then fails=$fails+1; fi

      ./scripts/scripts-tests -f -m $SCREAM_MACHINE
      if [[ $? != 0 ]]; then fails=$fails+1; fi
    fi

    # We do NOT want to do any of the items below if we are running
    # scripts-tests because the merge could change the state of the
    # repo and you never want to do that during testing. Also, CIME
    # V1 cases do not build with SCREAM_FAKE_ONLY on.
    if [ -z "$SCREAM_FAKE_ONLY" ]; then

      if [[ $test_v0 == 1 || $test_v1 == 1 ]]; then
        # AT CIME runs may need an upstream merge in order to ensure that any DIFFs
        # are caused by this PR and not simply because the PR is too far behind master.
        if [ -n "$PULLREQUESTNUM" ]; then
          ./scripts/git-merge-ref origin/master
          if [[ $? != 0 ]]; then
              echo "MERGE FAILED! Please resolve conflicts"
              exit 1
          fi
        fi
      fi

      if [[ $test_v0 == 1 ]]; then
        ../../cime/scripts/create_test e3sm_scream_v0 -c -b master
        if [[ $? != 0 ]]; then fails=$fails+1; fi
      fi

      if [[ $test_v1 == 1 ]]; then
        if [[ $is_at_run == 1 ]]; then
          # AT runs should be fast. => run only low resolution
          ../../cime/scripts/create_test e3sm_scream_v1_lowres --compiler=gnu9 -c -b master
        else
          # For nightlies, run both lowres and midres
          ../../cime/scripts/create_test e3sm_scream_v1 --compiler=gnu9 -c -b master
        fi
        if [[ $? != 0 ]]; then fails=$fails+1; fi
      else
          echo "SCREAM v1 tests were skipped, since the Github label 'AT: Skip v1 Testing' was found.\n"
      fi
    fi
  fi

  if [[ $fails > 0 ]]; then
      exit 1
  fi

else
  echo "Tests were skipped, since the Github label 'AT: Integrate Without Testing' was found.\n"
fi
