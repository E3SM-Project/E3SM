createLSFHeader() {

  RUN_SCRIPT=$1

  #delete the file if it exists
  rm -f $RUN_SCRIPT

  # Set up some yellowstone boiler plate
  echo "#!/bin/bash" >> $RUN_SCRIPT
  echo ""  >> $RUN_SCRIPT # newlines

  echo "#BSUB -a poe" >> $RUN_SCRIPT

  if [ -n "$HOMME_ACCOUNT" ]; then
    echo "#BSUB -P $HOMME_ACCOUNT" >> $RUN_SCRIPT
  else
    echo "PROJECT CHARGE ID (HOMME_PROJID) not set"
    exit -1
  fi 

  echo "" >> $RUN_SCRIPT # newline

  echo "#BSUB -q small" >> $RUN_SCRIPT
  echo "#BSUB -W 0:40" >> $RUN_SCRIPT

  echo "" >> $RUN_SCRIPT # newline

  echo "#BSUB -x" >> $RUN_SCRIPT

  echo "" >> $RUN_SCRIPT # newline

  # Set the job name
  echo "#BSUB -J $testName" >> $RUN_SCRIPT
  echo "" >> $RUN_SCRIPT

  # Set the output and error filenames
  echo "#BSUB -o $testName.stdout.%J" >> $RUN_SCRIPT
  echo "#BSUB -e $testName.stderr.%J" >> $RUN_SCRIPT
  echo "" >> $RUN_SCRIPT

  # Set the ncpus and ranks per MPI
  echo "#BSUB -n $num_cpus" >> $RUN_SCRIPT
  echo '#BSUB -R "span[ptile='$num_cpus']" ' >> $RUN_SCRIPT

  echo "" >> $RUN_SCRIPT

}

createPBSHeader() {

  RUN_SCRIPT=$1

  #delete the file if it exists
  rm -f $RUN_SCRIPT

  # Set up some yellowstone boiler plate
  echo "#!/bin/bash -l" >> $RUN_SCRIPT
  echo ""  >> $RUN_SCRIPT # newlines
  
  if [ -n "$HOMME_ACCOUNT" ]; then
    echo "#PBS -A $HOMME_ACCOUNT" >> $RUN_SCRIPT
  else
    echo "PROJECT CHARGE ID (HOMME_PROJID) not set"
    exit -2
  fi 

  # Set the output and error filenames
  echo "#PBS -o $testName.stdout.\${PBS_JOBID}" >> $RUN_SCRIPT
  echo "#PBS -e $testName.stderr.\${PBS_JOBID}" >> $RUN_SCRIPT

  echo "#PBS -N $testName" >> $RUN_SCRIPT

  echo "#PBS -l nodes=1" >> $RUN_SCRIPT
  echo "#PBS -l walltime=0:40:00" >> $RUN_SCRIPT
  # Not sure how to make the following portable
  #echo "#PBS -l gres=widow1" >> $RUN_SCRIPT

  echo "" >> $RUN_SCRIPT

}

createStdHeader() {

  RUN_SCRIPT=$1

  echo "#!/bin/bash " >> $RUN_SCRIPT

  echo "" >> $RUN_SCRIPT

}

printSubmissionSummary() {

  # Output a summary of the test name along with the bsub job id for easy reference
  echo "" # newline

  if [ "${SUBMIT_ALL_AT_ONCE}" == true ] ; then
    echo "############################################################################"
    echo "Summary of submissions"
  fi
  echo "############################################################################"
  for i in $(seq 0 $(( ${#SUBMIT_TEST[@]} - 1)))
  do
    echo "Test ${SUBMIT_TEST[$i]} has ID ${SUBMIT_JOB_ID[i]}"
  done
  echo "############################################################################"

}

examineJobStat() {
  # submit the job to the queue
  if [ "$HOMME_Submission_Type" = lsf ]; then
    jobStat=`bjobs -a $jobID 2>&1 | tail -n 1 | awk '{print $3}'`
    if [ -n "$jobStat" ] ; then
      if [ "$jobStat" == "PEND" ] ; then
        jobStat="pending" 
      elif [ "$jobStat" == "RUN" ]; then
        jobStat="running" 
      fi
    else
      jobStat="completed" 
    fi
  elif [ "$HOMME_Submission_Type" = pbs ]; then
    jobStat=`qstat $jobID 2>&1 | tail -n 1 | awk '{print $5}'`
    if [ -n "$jobStat" ] ; then
      if [ "$jobStat" == "Q" ] ; then
        jobStat="pending" 
      elif [ "$jobStat" == "R" ]; then
        jobStat="running" 
      fi
    else
      jobStat="completed" 
    fi
  elif [ "$HOMME_Submission_Type" = moab ] || [ "$HOMME_Submission_Type" = slurm ] ; then
    jobStat=`squeue -j $jobID 2>&1 | tail -n 1 | awk '{print $5}'`
    if [ -n "$jobStat" ] ; then
      if [ "$jobStat" == "PD" ] ; then
        jobStat="pending" 
      elif [ "$jobStat" == "R" ]; then
        jobStat="running" 
      fi
    else
      jobStat="completed" 
    fi
  else
    echo "Error: queue type not recognized"
    exit -3
  fi
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
      examineJobStat
      #jobStat=`bjobs -a $jobID | tail -n 1 | awk '{print $3}'`

      # if the job is registered in the queue and the status is PEND or RUN then wait
      if [ "$jobStat" == "pending" -o "$jobStat" == "running" ]; then
        # Print the status of the job
        echo -n "$jobStat..."

        # Job still in queue or running
        sleep 20 # sleep for 20s
      else # if jobStat=DONE, EXIT or it is finished and no longer registered in the queue
        jobFinished=true
        echo "finished."
      fi
    done
  done
  echo "############################################################################"

}

submitToQueue() {
  # submit the job to the queue
  if [ "$HOMME_Submission_Type" = lsf ]; then
    bsub < ${subFile} > $THIS_STDOUT 2> $THIS_STDERR
  elif [ "$HOMME_Submission_Type" = pbs ]; then
    qsub ${subFile} > $THIS_STDOUT 2> $THIS_STDERR
  elif [ "$HOMME_Submission_Type" = moab ]; then
    msub ${subFile} > $THIS_STDOUT 2> $THIS_STDERR
  elif [ "$HOMME_Submission_Type" = slurm ]; then
    sbatch ${subFile} > $THIS_STDOUT 2> $THIS_STDERR
  else
    echo "Error: queue type not recognized"
    exit -4
  fi
}

getJobID() {
  # submit the job to the queue
  if [ "$HOMME_Submission_Type" = lsf ]; then
    SUB_ID=`cat $THIS_STDOUT | awk '{print $2}' | sed  's/<//' | sed -e 's/>//'`
  elif [ "$HOMME_Submission_Type" = pbs ]; then
    SUB_ID=`cat $THIS_STDOUT | awk '{print $1}'`
  elif [ "$HOMME_Submission_Type" = moab ]; then
    SUB_ID=`cat $THIS_STDOUT | awk '{print $1}'`
  elif [ "$HOMME_Submission_Type" = slurm ]; then
    SUB_ID=`cat $THIS_STDOUT | tail -n 1 | awk '{print $4}'`
  else
    echo "Error: queue type not recognized"
    exit -5
  fi
}

submitTestsToQueue() {

  if [ "${CREATE_BASELINE}" == true ] ; then
    cd ${HOMME_DEFAULT_BASELINE_DIR}
  else
    cd ${HOMME_TESTING_DIR}
  fi

  if [ "${SUBMIT_ALL_AT_ONCE}" == true ] ; then
    echo "Submitting ${num_submissions} jobs to queue"
  fi

  SUBMIT_TEST=()
  SUBMIT_JOB_ID=()
  SUBMIT_TEST=()

  # Loop through all of the tests
  for subNum in $(seq 1 ${num_submissions})
  do

    subFile=subFile${subNum}
    subFile=${!subFile}
    subJobName=`basename ${subFile} .sh`

    # setup file for stdout and stderr redirection
    THIS_STDOUT=${subJobName}.out
    THIS_STDERR=${subJobName}.err

    # Run the command
    echo -n "Submitting test ${subJobName} to the queue... "
    submitToQueue ${subFile} $THIS_STDOUT $THIS_STDERR

    SUB_STAT=$?

    # Do some error checking
    if [ $SUB_STAT == 0 ]; then # the command was succesful
      # Get the queue job ID
      getJobID 

      echo "successful job id = $SUB_ID"
      SUBMIT_TEST+=( "${subJobName}" )
      SUBMIT_JOB_ID+=( "$SUB_ID" )
    else 
      echo "failed with message:"
      cat $THIS_STDERR
      exit -6
    fi
    rm $THIS_STDOUT
    rm $THIS_STDERR
  done
}

runTestsStd() {

  echo "Submitting ${num_submissions} jobs"

  SUBMIT_TEST=()
  SUBMIT_JOB_ID=()
  SUBMIT_TEST=()

  # Loop through all of the tests
  for subNum in $(seq 1 ${num_submissions})
  do

    subFile=subFile${subNum}
    subFile=${!subFile}
    #echo "subFile=${subFile}"
    subJobName=`basename ${subFile} .sh`
    #echo "subJobName=$subJobName"

    # setup file for stdout and stderr redirection
    THIS_STDOUT=${subJobName}.out
    THIS_STDERR=${subJobName}.err

    # Run the command
    # For some reason bsub must not be part of a string
    echo -n "Running test ${subJobName} ... "
    #echo "${subFile} > $THIS_STDOUT 2> $THIS_STDERR"
    chmod u+x ${subFile}
    cmd="${subFile} > $THIS_STDOUT 2> $THIS_STDERR"
    echo "$cmd"
    $cmd
    # Get the status of the run
    RUN_STAT=$?
    # Do some error checking
    if [ $RUN_STAT = 0 ]; then
      # the command was succesful
      echo "test ${subJobName} was run successfully"
      SUBMIT_TEST+=( "${subJobName}" )
      SUBMIT_JOB_ID+=( "${RUN_PID}" )
    else 
      echo "failed with message:"
      cat $THIS_STDERR
      exit -7
    fi
    rm $THIS_STDOUT
    rm $THIS_STDERR
  done
}

createAllRunScripts() {
  touch $lsfListFile
  echo "num_submissions=$num_test_files" > $lsfListFile

  for testFileNum in $(seq 1 $num_test_files)
  do

    testFile=test_file${testFileNum}
    source ${!testFile}

    testName=`basename ${!testFile} .sh`

    echo "Test $testName has $num_tests pure MPI tests"
    if [ -n "$omp_num_tests" ]; then
      echo "  and $omp_num_tests Hybrid MPI + OpenMP tests"
    fi

    if [ "${CREATE_BASELINE}" == true ] ; then
      # Create the run script
      thisRunScript="${HOMME_DEFAULT_BASELINE_DIR}/$testName/$testName-run.sh"
      outputDir="${HOMME_DEFAULT_BASELINE_DIR}/$testName"
    else
      # Create the run script
      thisRunScript=`dirname ${!testFile}`/$testName-run.sh
      outputDir=`dirname ${!testFile}`
    fi

    # delete the run script file if it exists
    rm -f ${thisRunScript}

    # Set up header
    submissionHeader $thisRunScript

    # cd into the correct dir
    echo "cd $outputDir" >> $RUN_SCRIPT
    echo "" >> $RUN_SCRIPT # new line

    for testNum in $(seq 1 $num_tests)
    do
      testExec=test${testNum}
      echo "# Pure MPI test ${testNum}" >> $thisRunScript
      #echo "mpiexec -n $num_cpus ${!testExec} > $testName.out 2> $testName.err" >> $thisRunScript
      execLine $thisRunScript "${!testExec}" $num_cpus
      echo "" >> $thisRunScript # new line
    done

    if [ -n "$omp_num_tests" -a "${RUN_OPENMP}" = TRUE ]; then
      echo "export OMP_NUM_THREADS=$omp_number_threads" >> $thisRunScript
      echo "export OMP_STACKSIZE=128M" >> $thisRunScript
      echo "" >> $thisRunScript # new line
      for testNum in $(seq 1 $omp_num_tests)
      do
         testExec=omp_test${testNum}
         echo "# Hybrid test ${testNum}" >> $thisRunScript
         #echo "mpiexec -n $omp_num_mpi ${!testExec} > $testName.out 2> $testName.err" >> $thisRunScript
         execLine $thisRunScript "${!testExec}" $omp_num_mpi
         echo "" >> $thisRunScript # new line
      done
    fi

    ############################################################
    # At this point, the submission files have been created.
    # For the baseline, we are finished. 
    # For the cprnc diffing we need to add that
    ############################################################
    if [ "${CREATE_BASELINE}" == false ] ; then 
      mkdir -p ${HOMME_DEFAULT_BASELINE_DIR}/${testName}
      echo "cp $thisRunScript ${HOMME_BASELINE_DIR}/${testName}/"
      cp $thisRunScript ${HOMME_BASELINE_DIR}/${testName}/

      ############################################################
      # Now set up the cprnc diffing
      ############################################################
      # load the cprnc files for this run
      FILES="${nc_output_files}"

      if [ -z "${FILES}" ] ; then
          echo "Test ${TEST_NAME} doesn't have Netcdf output files"
      fi

      # for files in movies
      for file in $FILES 
      do
        echo "file = ${file}"
        baseFilename=`basename $file`

        # new result
        newFile=${HOMME_TESTING_DIR}/${testName}/movies/$file

        #if [ ! -f "${newFile}" ] ; then
        #  echo "Error: The result file ${newFile} does not exist. Exiting" 
        #  exit -1
        #fi

        # result in the repo
        repoFile=${HOMME_BASELINE_DIR}/${testName}/movies/${baseFilename}

        #if [ ! -f "${newFile}" ] ; then
        #  echo "Warning: The repo file ${repoFile} does not exist in the baseline dir"
        #  exit 0
        #fi


        diffStdout=${testName}.${baseFilename}.out
        diffStderr=${testName}.${baseFilename}.err

        echo "# Running cprnc to difference ${baseFilename} against baseline " >> $thisRunScript
        #echo "$cmd > $diffStdout 2> $diffStderr" >> $thisRunScript
        cmd="${CPRNC_BINARY} ${newFile} ${repoFile} > $diffStdout 2> $diffStderr"
        #echo "  $cmd"
        serExecLine $thisRunScript "$cmd"
        echo "" >> $thisRunScript # blank line
      done

    fi
    ############################################################
    # Finished setting up cprnc
    ############################################################

    echo "subFile$testFileNum=$thisRunScript" >>  $lsfListFile

    # Make the script executable
    chmod u+x ${thisRunScript}

    # Reset the variables (in case they are not redefined in the next iteration)
    unset omp_num_tests
    unset num_tests

  done

}

createRunScript() {

  source ${THIS_TEST_FILE}

  testName=`basename ${THIS_TEST_FILE} .sh`

  echo "Test $testName has $num_tests pure MPI tests"
  if [ -n "$omp_num_tests" ]; then
    echo "  and $omp_num_tests Hybrid MPI + OpenMP tests"
  fi

  # Create the run script
  thisRunScript=`dirname ${THIS_TEST_FILE}`/$testName-run.sh
  rm -f ${thisRunScript}

  outputDir=`dirname ${THIS_TEST_FILE}`

  # Set up header
  submissionHeader $thisRunScript

  for testNum in $(seq 1 $num_tests)
  do
    testExec=test${testNum}
    echo "# Pure MPI test ${testNum}" >> $thisRunScript
    execLine $thisRunScript "${!testExec}" $num_cpus
    echo "" >> $thisRunScript # new line
  done

  if [ -n "$omp_num_tests" -a "${RUN_OPENMP}" == true ]; then
    echo "export OMP_NUM_THREADS=$omp_number_threads" >> $thisRunScript
    echo "" >> $thisRunScript # new line
    for testNum in $(seq 1 $omp_num_tests)
    do
       testExec=omp_test${testNum}
       echo "# Hybrid test ${testNum}" >> $thisRunScript
       execLine $thisRunScript "${!testExec}" $omp_num_mpi
       echo "" >> $thisRunScript # new line
    done
  fi

  # make it executable
  chmod u+x ${thisRunScript}

  # Reset the variables (in case they are not redefined in the next iteration)
  #unset omp_num_tests
  #unset num_tests

}

parseHeader() {
  cat $1 | \
  sed -e "s/\${testName}/${testName}/g" \
      -e "s;\${testDir};${outputDir};g" \
      -e "s;\${HOMME_PROJID};${HOMME_ACCOUNT};g" > $2
  echo ""  >> $2 # newline
  echo "########################################"  >> $2
  echo "# Above is header, below is run data"  >> $2
  echo "########################################"  >> $2
  echo ""  >> $2 # newline
}

submissionHeader() {
  RUN_SCRIPT=$1

  if [ "$HOMME_Submission_Type" != none ]; then 
    # Use a user defined header
    if [ -n "$HOMME_Submission_Header" ]; then
      parseHeader ${HOMME_Submission_Header} ${RUN_SCRIPT}
    else
      # Build a header for lsf or pbs
      if [ "$HOMME_Submission_Type" = lsf ]; then
        createLSFHeader $RUN_SCRIPT
      elif [ "$HOMME_Submission_Type" = pbs ]; then
        createPBSHeader $RUN_SCRIPT
      fi
    fi
  else
    createStdHeader $RUN_SCRIPT
  fi

}

execLine() {
  RUN_SCRIPT=$1
  EXEC=$2
  NUM_CPUS=$3

  if [ -n "${MPI_EXEC}" ]; then
    # mpirun.lsf is a special case
    if [ "${MPI_EXEC}" = "mpirun.lsf" ] ; then
      echo "mpirun.lsf -pam \"-n ${NUM_CPUS}\" ${MPI_OPTIONS} $EXEC" >> $RUN_SCRIPT
    else
      echo "${MPI_EXEC} -n ${NUM_CPUS} ${MPI_OPTIONS} $EXEC" >> $RUN_SCRIPT
    fi
  else
    if [ "$HOMME_Submission_Type" = lsf ]; then
      echo "mpirun.lsf -pam \"-n ${NUM_CPUS}\" ${MPI_OPTIONS} $EXEC" >> $RUN_SCRIPT

    elif [ "$HOMME_Submission_Type" = pbs ]; then
      echo "aprun -n ${NUM_CPUS} ${MPI_OPTIONS} $EXEC" >> $RUN_SCRIPT

    else
      echo "mpiexec -n ${NUM_CPUS} ${MPI_OPTIONS} $EXEC" >> $RUN_SCRIPT

    fi
  fi
}

serExecLine() {
  RUN_SCRIPT=$1
  EXEC=$2

  if [ -n "${MPI_EXEC}" ]; then
    # mpirun.lsf is a special case
    if [ "${MPI_EXEC}" = "mpirun.lsf" ] ; then
      echo "$EXEC" >> $RUN_SCRIPT
    else
      echo "${MPI_EXEC} -n 1 ${MPI_OPTIONS} $EXEC" >> $RUN_SCRIPT
    fi

  else
    if [ "$HOMME_Submission_Type" = lsf ]; then
      echo "$EXEC" >> $RUN_SCRIPT
    elif [ "$HOMME_Submission_Type" = pbs ]; then
      echo "aprun -n 1 ${MPI_OPTIONS} $EXEC" >> $RUN_SCRIPT
    else
      echo "mpiexec -n 1 ${MPI_OPTIONS} $EXEC" >> $RUN_SCRIPT
    fi
  fi
}



diffCprnc() {

  if [ ! -f "${CPRNC_BINARY}" ] ; then
    echo "Netcdf differencing tool cprnc not found"
    exit -8
  fi

  # source the test.sh file to get the names of the nc_output_files
  source ${HOMME_TESTING_DIR}/${TEST_NAME}/${TEST_NAME}.sh

  # nc_output_files is defined in the .sh file
  FILES="${nc_output_files}"

  if [ -z "${FILES}" ] ; then
      echo "Test ${TEST_NAME} doesn't have Netcdf output files"
  fi

  # for files in movies
  for file in $FILES 
  do
    echo "file = ${file}"
    baseFilename=`basename $file`

    # new result
    newFile=${HOMME_TESTING_DIR}/${TEST_NAME}/movies/$file
    if [ ! -f "${newFile}" ] ; then
      echo "ERROR: The result file ${newFile} does not exist exiting" 
      exit -9
    fi

    # result in the repo
    #repoFile=${HOMME_NC_RESULTS_DIR}/${TEST_NAME}/${baseFilename}
    repoFile=${HOMME_BASELINE_DIR}/${TEST_NAME}/movies/${baseFilename}

    if [ ! -f "${newFile}" ] ; then
      echo "ERROR: The repo file ${repoFile} does not exist exiting" 
      exit -10
    fi

    cmd="${CPRNC_BINARY} ${newFile} ${repoFile}"

    diffStdout=${TEST_NAME}.${baseFilename}.out
    diffStderr=${TEST_NAME}.${baseFilename}.err

    echo "Running cprnc:"
    echo "  $cmd"
    $cmd > $diffStdout 2> $diffStderr

    # Parse the output file to determine if they were identical
    DIFF_RESULT=`grep -e 'diff_test' $diffStdout | awk '{ print $8 }'`

    if [ "${DIFF_RESULT}" == IDENTICAL ] ; then
      echo "The files are identical: DIFF_RESULT=${DIFF_RESULT}"
      # Delete the output file to remove clutter
      rm $diffStdout
      rm $diffStderr
    else
      echo "The files are different: DIFF_RESULT=${DIFF_RESULT}"
      #echo "  The diff output is available in $diffStdout"
      echo "############################################################################"
      echo "CPRNC returned the following RMS differences"
      grep RMS ${diffStdout}
      echo "############################################################################"
      exit -11
    fi

    
  done
}

diffCprncOutput() {

  # source the test.sh file to get the names of the nc_output_files
  source ${HOMME_TESTING_DIR}/${TEST_NAME}/${TEST_NAME}.sh

  # nc_output_files is defined in the .sh file
  FILES="${nc_output_files}"

  if [ -z "${FILES}" ] ; then
      echo "Test ${TEST_NAME} doesn't have Netcdf output files"
  fi

  # for files in movies
  for file in $FILES 
  do
    echo "file = ${file}"
    cprncOutputFile="${HOMME_TESTING_DIR}/${TEST_NAME}/${TEST_NAME}.`basename $file`.out"
    cprncErrorFile="${HOMME_TESTING_DIR}/${TEST_NAME}/${TEST_NAME}.`basename $file`.err"

    # ensure that cprncOutputFile exists
    if [ ! -f "${cprncOutputFile}" ]; then
      echo "Error: cprnc output file ${cprncOutputFile} not found. Exiting."
      exit -12
    fi

    # Parse the output file to determine if they were identical
    DIFF_RESULT=`grep -e 'diff_test' ${cprncOutputFile} | awk '{ print $8 }'`

    if [ "${DIFF_RESULT}" == IDENTICAL ] ; then
      echo "The files are identical: DIFF_RESULT=${DIFF_RESULT}"
      # Delete the output file to remove clutter
    else
      echo "The files are different: DIFF_RESULT=${DIFF_RESULT}"
      echo "############################################################################"
      echo "CPRNC returned the following RMS differences"
      grep RMS ${cprncOutputFile}
      echo "############################################################################"
      exit -13
    fi
    
  done
}

moveBaseline() {

  source ${HOMME_TESTING_DIR}/test_list.sh

  for subNum in $(seq 1 ${num_test_files})
  do

    subFile=test_file${subNum}
    subFile=${!subFile}

    subDirName=`dirname ${subFile}`
    subBaseName=`basename ${subDirName}`

    baselineDir=${HOMME_BASELINE_DIR}/$subBaseName
    mkdir -p $baselineDir

    # source the test.sh file to get the name of the nc_output_files
    source ${subFile}

    # nc_output_files is defined in the .sh file
    FILES="${nc_output_files}"

    if [ -z "${FILES}" ] ; then
      : # pass for now
      #echo "Test ${subBaseName} doesn't have Netcdf output files"
    fi

    # for files in movies
    for file in $FILES 
    do
      #echo "file = ${file}"
      baseFilename=`basename $file`

      # new result
      newFile=${subDirName}/movies/$file
      newFileBasename=`basename $newFile`

      if [ ! -f "${newFile}" ] ; then
        echo "ERROR: The result file ${newFile} does not exist exiting" 
        exit -17
      fi

      # Make the directory to store the netcdf files
      cmd="mkdir -p $baselineDir/movies"
      #echo $cmd
      $cmd
      # Move the netcdf files
      cmd="mv $newFile $baselineDir/movies/$newFileBasename"
      #echo "$cmd"
      $cmd
      
    done
  done
}

checkBaselineResults() {

  filesNotFound=false

  for testFileNum in $(seq 1 $num_test_files)
  do

    testFile=test_file${testFileNum}
    source ${!testFile}

    testName=`basename ${!testFile} .sh`

    FILES="${nc_output_files}"


    if [ -z "${FILES}" ] ; then
      echo "Test ${TEST_NAME} doesn't have Netcdf output files"
    fi

    # for files in movies
    for file in $FILES 
    do
      baseFilename=`basename $file`

      # result in the repo
      repoFile=${HOMME_BASELINE_DIR}/${testName}/movies/${baseFilename}

      if [ ! -f "${repoFile}" ] ; then
        echo "Error: The Netcdf file ${repoFile} does not exist in the baseline dir"
        filesNotFound=true
      fi

    done

  done

  # Throw an error if any of the files to be compared were not found
  if $filesNotFound ; then
    exit -18
  fi

}


