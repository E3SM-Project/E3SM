#!/bin/bash

# Strip the end of the stdout file which contains lsf dump
stripAppendage() {
  sed -i -e '/^ exiting/,/^ number of MPI/{/^ exiting/!{/^ number of MPI/!d}}' -e '/^process/d' $1
}

readTestList() {

  source ./test-list.in

  # Define some arrays
  TEST_NAMES=()
  EXEC_NAMES=()
  OMP_TEST=()

  # Then read the file from the test input files
  for testFile in "${TESTS[@]}"
  do
    source $testFile
    TEST_NAMES+=( "${test_name}" )
    EXEC_NAMES+=( "${exec_name}" )
    if [ -n "$OMP_SUB_TESTS" ] ; then
      OMP_SUB_TEST+=( "$OMP_SUB_TESTS" )
    else 
      OMP_SUB_TEST+=( "false" )
    fi
  done

}

setTestDirs() {

  # Determine some locations
  DIR_NAME=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd -P)

  # Determine the location of the tests results (for yellowstone)
  RESULT_DIR=$(cd "${DIR_NAME}/../../results/yellowstone" && pwd -P)

  # Set the location of the "build" base directory
  if [ -n "$1" -a -d "$1" ]; then
    # Use this location as the base of the file structure for the tests
    BUILD_DIR=$1
  else
    # Set the build directory from the set file structure
    BUILD_DIR=$(cd `dirname $DIR_NAME/../../../../..` && pwd -P)/build
  fi

}

submitTestsToLSF() {

  # Set up an array to catch the job names and bsub ids
  SUBMIT_TEST=()
  SUBMIT_TEST_DIR=()
  SUBMIT_JOB_ID=()

  # Loop through all of the tests
  for testNum in $(seq 0 $(( ${#TEST_NAMES[@]} - 1)))
  do

    test_name=${TEST_NAMES[$testNum]}
    exec_name=${EXEC_NAMES[$testNum]}
    test_dir="${BUILD_DIR}/${exec_name}/${test_name}"
    script_name="${test_dir}/${test_name}-run.sh"

    # the omp tests use the same executable as the non-omp tests check those executables
    if [ "${test_name:(-4)}" == "-omp" ]; then
      stripOmp=${test_name%"-omp"}
      executable="${BUILD_DIR}/${exec_name}/${stripOmp}/${exec_name}"
    else
      executable="${BUILD_DIR}/${exec_name}/${test_name}/${exec_name}"
    fi

    #ensure that the script and the executable exist
    if [ ! -f "${script_name}" ]; then      
      echo "Error: script `basename ${script_name}` does not exist"
      echo "  File should be ${script_name}"
      exit -1
    fi
    if [ ! -f "${executable}" ]; then
      echo "Error: executable `basename ${executable}` does not exist"
      echo "  File should be ${executable}"
      exit -2
    fi

    # setup file for stdout and stderr redirection
    THIS_STDOUT=$DIR_NAME/${test_name}.out
    THIS_STDERR=$DIR_NAME/${test_name}.err

    # Run the command
    #echo "Submitting bsub < $script_name"
    # For some reason bsub must not be part of a string
    echo -n "Submitting test ${exec_name}/${test_name} to the queue... "
    bsub < ${script_name} > $THIS_STDOUT 2> $THIS_STDERR
    BSUB_STAT=$?

    # Do some error checking
    if [ $BSUB_STAT == 0 ]; then
      # the command was succesful
      BSUB_ID=`cat $THIS_STDOUT | awk '{print $2}' | sed  's/<//' | sed -e 's/>//'`
      echo "successful job id = $BSUB_ID"
      SUBMIT_TEST+=( "${test_name}" )
      SUBMIT_JOB_ID+=( "$BSUB_ID" )
      SUBMIT_TEST_DIR+=( "$test_dir" )
    else 
      echo "failed with message:"
      cat $THIS_STDERR
      exit -1
    fi
    rm $THIS_STDOUT
    rm $THIS_STDERR
  done

}

printSubmissionSummary() {

  # Output a summary of the test name along with the bsub job id for easy reference
  echo "" # newline
  echo "############################################################################"
  echo "Summary of submissions"
  echo "############################################################################"
  for i in $(seq 0 $(( ${#SUBMIT_TEST[@]} - 1)))
  do
    echo "Test ${SUBMIT_TEST[$i]} has ID ${SUBMIT_JOB_ID[i]}"
  done
  echo "############################################################################"

}

queueWait() {

  for i in $(seq 0 $(( ${#SUBMIT_TEST[@]} - 1)))
  do
    echo -n "Examining status of ${SUBMIT_TEST[$i]}..."
    jobID=${SUBMIT_JOB_ID[i]}
    jobFinished=false

    while ! $jobFinished;
    do
      # Test if the job exists
      jobStat=`bjobs -a $jobID | tail -n 1 | awk '{print $3}'`

      # Print the status of the job
      echo -n "$jobStat..."

      # if the job is registered in the queue and the status is PEND or RUN then wait
      if [ -n "$jobStat" -a "$jobStat" == "PEND" -o "$jobStat" == "RUN" ]; then
        # Job still in queue or running
        sleep 60 # sleep for 60s
      else # if jobStat=DONE, EXIT or it is finished and no longer registered in the queue
        jobFinished=true
        echo "FINISHED..."
      fi
    done
  done
  echo "############################################################################"

}

diffStdOut() {

  # Should be a unique file
  diffFile="diff.${SUBMIT_JOB_ID[0]}"
  echo "Concatenating all diff output into $diffFile"

  # Then diff with the stored results (yellowstone only)
  for i in $(seq 0 $(( ${#SUBMIT_TEST[@]} - 1)))
  do
    THIS_TEST=${SUBMIT_TEST[$i]}
    # The following is not very clean
    NEW_RESULT=${DIR_NAME}/${THIS_TEST}.stdout.${SUBMIT_JOB_ID[$i]}
    SAVED_RESULT=${RESULT_DIR}/${THIS_TEST}/${THIS_TEST}.stdout
    # TODO: Ensure that teh files exist



    stripAppendage $NEW_RESULT
    echo "diff $NEW_RESULT $SAVED_RESULT" >> $diffFile
    # append the output to 
    diff $NEW_RESULT $SAVED_RESULT >> $diffFile
  done
}


