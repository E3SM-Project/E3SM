parseHeader() {
  cat $1 | \
  sed -e "s/\${TEST_NAME}/${TEST_NAME}/g" \
      -e "s;\${TEST_DIR};${outputDir};g" \
      -e "s;\${NUM_CPUS};${NUM_CPUS};g" \
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
    # Use a provided or user defined header
    if [ -n "$HOMME_Submission_Header" ]; then
      parseHeader ${HOMME_Submission_Header} ${RUN_SCRIPT}
    else
      echo "Error: No queue subimission header supplied"
      exit -19
    fi
  else
    createStdHeader $RUN_SCRIPT
  fi

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
    chmod u+x ${subFile}
    #cmd="${subFile} > $THIS_STDOUT 2> $THIS_STDERR"
    cmd="${subFile}"
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
      echo "ERROR: run failed. check out/err files in:"
      echo `dirname ${subFile}`
      exit -7
    fi
  done
}


createAllRunScripts() {
  touch $subListFile
  echo "num_submissions=${NUM_TEST_FILES}" > $subListFile

  for testFileNum in $(seq 1 ${NUM_TEST_FILES})
  do
    #add summit-p9 here later
    if [ "${HOMME_MACHINE}" == "summit-gpu" ] ; then
       createRunScriptSummit
    else
       createAllRunScriptsGeneric
    fi

    echo "subFile$testFileNum=$thisRunScript" >>  $subListFile

    # Make the script executable
    chmod u+x ${thisRunScript}

    # Reset the variables (in case they are not redefined in the next iteration)
    unset OMP_NUM_TESTS
    unset NUM_TESTS

  done
}

emptyLine(){
echo "" >> $1
}


createRunScriptSummit() {
    testFile=TEST_FILE_${testFileNum}
    source ${!testFile}

    TEST_NAME=`basename ${!testFile} .sh`

    if [ "${CREATE_BASELINE}" == true ] ; then
      # Create the run script
      thisRunScript="${HOMME_DEFAULT_BASELINE_DIR}/${TEST_NAME}/${TEST_NAME}-run.sh"
      outputDir="${HOMME_DEFAULT_BASELINE_DIR}/${TEST_NAME}"
    else
      # Create the run script
      thisRunScript=`dirname ${!testFile}`/${TEST_NAME}-run.sh
      outputDir=`dirname ${!testFile}`
    fi

    # delete the run script file if it exists
    rm -f ${thisRunScript}

    # Set up header
    #add more for sbatch
    echo "#!/bin/bash" > $thisRunScript
    emptyLine $thisRunScript

    # cd into the correct dir
    echo "cd $outputDir " >> $thisRunScript
    emptyLine $thisRunScript

    # Remove all existing netcdf files
    echo "rm -f movies/* " >> $thisRunScript
    emptyLine $thisRunScript

    #check if this is a gpu exec
    if [[ "$EXEC" == *"kokkos"* ]] ; then
      echo "source ${HOMME_DIR}/${SUMMIT_MODULES_GPU}">> $thisRunScript
    else
      echo "source ${HOMME_DIR}/${SUMMIT_MODULES_P9}">> $thisRunScript
    fi
    emptyLine $thisRunScript

    #kokkos needs omp_num_threads set
    if [ -n "${OMP_NUMBER_THREADS_KOKKOS}" ]; then
      echo "export OMP_NUM_THREADS=${OMP_NUMBER_THREADS_KOKKOS}" >> $thisRunScript
      emptyLine $thisRunScript
    fi

    testExec=$TEST_1

    #check if this is a gpu exec
    if [[ "$EXEC" == *"kokkos"* ]] ; then
      echo "${SUMMIT_JSRUN_GPU} \\" >> $thisRunScript
      echo "${testExec} \\" >> $thisRunScript
      echo "${SUMMIT_JSRUN_TAIL} > ${TEST_NAME}_1.out 2> ${TEST_NAME}_1.err" >>$thisRunScript
      emptyLine $thisRunScript
    fi   

#cprnc
if true; then
    if [ "${CREATE_BASELINE}" == false ] ; then
      mkdir -p ${HOMME_DEFAULT_BASELINE_DIR}/${TEST_NAME}

      ############################################################
      # Now set up the cprnc diffing against baslines
      ############################################################
      # load the cprnc files for this run
      FILES="${NC_OUTPUT_FILES}"

      if [ -z "${FILES}" ] ; then
          echo "Test ${TEST_NAME} doesn't specify any baseline comparison tests"
      fi

      # for files in movies
      for file in $FILES
      do
        echo "file = ${file}"
        baseFilename=`basename $file`

        # new result
        newFile=${HOMME_TESTING_DIR}/${TEST_NAME}/movies/$file

        # result in the repo
        repoFile=${HOMME_BASELINE_DIR}/${TEST_NAME}/movies/${baseFilename}

        diffStdout=${TEST_NAME}.${baseFilename}.out
        diffStderr=${TEST_NAME}.${baseFilename}.err

        echo "# Running cprnc to difference ${baseFilename} against baseline " >> $thisRunScript
        #echo "$cmd > $diffStdout 2> $diffStderr" >> $thisRunScript
        cmd="${SUMMIT_JSRUN_SERIAL} ${CPRNC_BINARY} -m ${repoFile} ${newFile} > $diffStdout 2> $diffStderr"
        #echo "  $cmd"
        echo ${cmd} >> $thisRunScript
        emptyLine $thisRunScript
      done
      ############################################################
      # Now set up the cprnc diffing against REF solutions
      # these could be internally generated, or included in the repo
      ############################################################
      # load the cprnc files for this run
      REFFILES=( ${NC_OUTPUT_REF} )
      FILES="${NC_OUTPUT_CHECKREF}"

      if [ -z "${FILES}" ] ; then
          echo "Test ${TEST_NAME} doesn't specify any reference file tests"
      fi

      # for files in movies
      COUNT=0
      for file in $FILES
      do
        refname=${REFFILES[$COUNT]}
        echo "ref = ${refname}"
        echo "file = ${file}"
        baseFilename=`basename $file`

        # new result
        newFile=${HOMME_TESTING_DIR}/${TEST_NAME}/movies/${file}
        refFile=${HOMME_TESTING_DIR}/${TEST_NAME}/movies/${refname}

        diffStdout=${TEST_NAME}.ref.${baseFilename}.out
        diffStderr=${TEST_NAME}.ref.${baseFilename}.err

        echo "# Running cprnc to difference ${baseFilename} against reference " >> $thisRunScript
        cmd="${SUMMIT_JSRUN_SERIAL} ${CPRNC_BINARY} -m ${refFile} ${newFile} > $diffStdout 2> $diffStderr"
        echo ${cmd} >> $thisRunScript
        emptyLine $thisRunScript
        echo "" >> $thisRunScript # blank line
        let COUNT+=1
      done

    fi
    ############################################################
    # Finished setting up cprnc
    ############################################################

fi



}




createAllRunScriptsGeneric() {

    testFile=TEST_FILE_${testFileNum}
    source ${!testFile}

    TEST_NAME=`basename ${!testFile} .sh`

    echo "Test ${TEST_NAME} has ${NUM_TESTS} pure MPI tests"
    if [ -n "${OMP_NUM_TESTS}" ]; then
      echo "  and ${OMP_NUM_TESTS} Hybrid MPI + OpenMP tests"
    fi

    if [ "${CREATE_BASELINE}" == true ] ; then
      # Create the run script
      thisRunScript="${HOMME_DEFAULT_BASELINE_DIR}/${TEST_NAME}/${TEST_NAME}-run.sh"
      outputDir="${HOMME_DEFAULT_BASELINE_DIR}/${TEST_NAME}"
    else
      # Create the run script
      thisRunScript=`dirname ${!testFile}`/${TEST_NAME}-run.sh
      outputDir=`dirname ${!testFile}`
    fi

    # delete the run script file if it exists
    rm -f ${thisRunScript}

    # Set up header
    submissionHeader $thisRunScript

    # cd into the correct dir
    echo "cd $outputDir" >> $RUN_SCRIPT
    echo "" >> $RUN_SCRIPT # new line

    # Remove all existing netcdf files
    echo "rm -f movies/*" >> $RUN_SCRIPT
    echo "" >> $RUN_SCRIPT # new line

    #kokkos needs omp_num_threads set
    if [ -n "${OMP_NUMBER_THREADS_KOKKOS}" ]; then
      echo "export OMP_NUM_THREADS=${OMP_NUMBER_THREADS_KOKKOS}" >> $thisRunScript
      #do we need this?
      echo "export OMP_STACKSIZE=128M" >> $thisRunScript
      echo "" >> $thisRunScript # new line
    fi

    for testNum in $(seq 1 ${NUM_TESTS})
    do
      testExec=TEST_${testNum}
      echo "# Pure MPI test ${testNum}" >> $thisRunScript
      #echo "mpiexec -n ${NUM_CPUS }${!testExec} > ${TEST_NAME}.out 2> ${TEST_NAME}.err" >> $thisRunScript
      execLine $thisRunScript "${!testExec}" ${NUM_CPUS} ${TEST_NAME}_${testNum}
      echo "" >> $thisRunScript # new line
    done

    if [ -n "${OMP_NUM_TESTS}" -a "${RUN_OPENMP}" = ON ]; then
      echo "export OMP_NUM_THREADS=${OMP_NUMBER_THREADS}" >> $thisRunScript
      echo "export OMP_STACKSIZE=128M" >> $thisRunScript
      echo "" >> $thisRunScript # new line
      for testNum in $(seq 1 ${OMP_NUM_TESTS})
      do
         testExec=OMP_TEST_${testNum}
         echo "# Hybrid test ${testNum}" >> $thisRunScript
         #echo "mpiexec -n $omp_num_mpi ${!testExec} > ${TEST_NAME}.out 2> ${TEST_NAME}.err" >> $thisRunScript
         execLine $thisRunScript "${!testExec}" ${OMP_NUM_MPI} ${TEST_NAME}_OMP${testNum}
         echo "" >> $thisRunScript # new line
      done
    fi

    ############################################################
    # At this point, the submission files have been created.
    # For the baseline, we are finished. 
    # For the cprnc diffing we need to add that
    ############################################################
    if [ "${CREATE_BASELINE}" == false ] ; then 
      mkdir -p ${HOMME_DEFAULT_BASELINE_DIR}/${TEST_NAME}

#      MT: why do this? it will trash the run scripts that were created by the baseline code?
#      echo "cp $thisRunScript ${HOMME_BASELINE_DIR}/${TEST_NAME}/"
#      cp $thisRunScript ${HOMME_BASELINE_DIR}/${TEST_NAME}/

      ############################################################
      # Now set up the cprnc diffing against baslines
      ############################################################
      # load the cprnc files for this run
      FILES="${NC_OUTPUT_FILES}"

      if [ -z "${FILES}" ] ; then
          echo "Test ${TEST_NAME} doesn't specify any baseline comparison tests"
      fi

      # for files in movies
      for file in $FILES 
      do
        echo "file = ${file}"
        baseFilename=`basename $file`

        # new result
        newFile=${HOMME_TESTING_DIR}/${TEST_NAME}/movies/$file

        # result in the repo
        repoFile=${HOMME_BASELINE_DIR}/${TEST_NAME}/movies/${baseFilename}

        diffStdout=${TEST_NAME}.${baseFilename}.out
        diffStderr=${TEST_NAME}.${baseFilename}.err

        echo "# Running cprnc to difference ${baseFilename} against baseline " >> $thisRunScript
        #echo "$cmd > $diffStdout 2> $diffStderr" >> $thisRunScript
        cmd="${CPRNC_BINARY} -m ${repoFile} ${newFile} > $diffStdout 2> $diffStderr"
        #echo "  $cmd"
        serExecLine $thisRunScript "$cmd"
        echo "" >> $thisRunScript # blank line
      done



      ############################################################
      # Now set up the cprnc diffing against REF solutions
      # these could be internally generated, or included in the repo
      ############################################################
      # load the cprnc files for this run
      REFFILES=( ${NC_OUTPUT_REF} )
      FILES="${NC_OUTPUT_CHECKREF}"

      if [ -z "${FILES}" ] ; then
          echo "Test ${TEST_NAME} doesn't specify any reference file tests"
      fi

      # for files in movies
      COUNT=0
      for file in $FILES 
      do
        refname=${REFFILES[$COUNT]}
        echo "ref = ${refname}"
        echo "file = ${file}"
        baseFilename=`basename $file`

        # new result
        newFile=${HOMME_TESTING_DIR}/${TEST_NAME}/movies/${file}
        refFile=${HOMME_TESTING_DIR}/${TEST_NAME}/movies/${refname}

        diffStdout=${TEST_NAME}.ref.${baseFilename}.out
        diffStderr=${TEST_NAME}.ref.${baseFilename}.err

        echo "# Running cprnc to difference ${baseFilename} against reference " >> $thisRunScript
        cmd="${CPRNC_BINARY} -m ${refFile} ${newFile} > $diffStdout 2> $diffStderr"
        serExecLine $thisRunScript "$cmd"
        echo "" >> $thisRunScript # blank line
        let COUNT+=1
      done

    fi
    ############################################################
    # Finished setting up cprnc
    ############################################################

}

createRunScript() {

  source ${THIS_TEST_FILE}

  TEST_NAME=`basename ${THIS_TEST_FILE} .sh`

  echo "Test ${TEST_NAME} has ${NUM_TESTS} pure MPI tests"
  if [ -n "${OMP_NUM_TESTS}" ]; then
    echo "  and ${OMP_NUM_TESTS} Hybrid MPI + OpenMP tests"
  fi

  # Create the run script
  thisRunScript=`dirname ${THIS_TEST_FILE}`/${TEST_NAME}-run.sh
  rm -f ${thisRunScript}

  outputDir=`dirname ${THIS_TEST_FILE}`

  # Set up header
  submissionHeader $thisRunScript

  for testNum in $(seq 1 ${NUM_TESTS})
  do
    testExec=TEST_${testNum}
    echo "# Pure MPI test ${testNum}" >> $thisRunScript
    execLine $thisRunScript "${!testExec}" ${NUM_CPUS} ${TEST_NAME}_${testNum}
    echo "" >> $thisRunScript # new line
  done

  if [ -n "${OMP_NUM_TESTS}" -a "${RUN_OPENMP}" = ON ]; then
    echo "export OMP_NUM_THREADS=${OMP_NUMBER_THREADS}" >> $thisRunScript
    echo "" >> $thisRunScript # new line
    for testNum in $(seq 1 ${OMP_NUM_TESTS})
    do
       testExec=${OMP_TEST_${testNum}}
       echo "# Hybrid test ${testNum}" >> $thisRunScript
       execLine $thisRunScript "${!testExec}" ${OMP_NUM_MPI} ${TEST_NAME}_OMP${testNum}
       echo "" >> $thisRunScript # new line
    done
  fi

  # make it executable
  chmod u+x ${thisRunScript}

  # Reset the variables (in case they are not redefined in the next iteration)
  #unset omp_num_tests
  #unset num_tests

}

execLine() {
  RUN_SCRIPT=$1
  EXEC=$2
  NUM_MPI_PROCS=$3
  OPT="> $4.out 2> $4.err"


  if [ -n "${MPI_EXEC}" ]; then
    # mpirun.lsf is a special case
    if [ "${MPI_EXEC}" = "mpirun.lsf" ] ; then
      echo "mpirun.lsf -pam \"-n ${NUM_MPI_PROCS}\" ${MPI_OPTIONS} $EXEC $OPT" >> $RUN_SCRIPT
    elif [ "${MPI_EXEC}" = "runjob" ]; then
      echo "runjob -n ${NUM_MPI_PROCS} ${MPI_OPTIONS} --block \$COBALT_PARTNAME --verbose=INFO : $EXEC $OPT" >> $RUN_SCRIPT
    elif [ "${MPI_EXEC}" = "aprun" ] ; then
      if [[ $4 == *"_OMP"* ]]; then
        echo "aprun -n ${NUM_MPI_PROCS} -d ${OMP_NUMBER_THREADS} ${MPI_OPTIONS} $EXEC $OPT" >> $RUN_SCRIPT
      else
        echo "aprun -n ${NUM_MPI_PROCS} ${MPI_OPTIONS} $EXEC $OPT" >> $RUN_SCRIPT
      fi
    else
      echo "${MPI_EXEC} -n ${NUM_MPI_PROCS} ${MPI_OPTIONS} $EXEC $OPT" >> $RUN_SCRIPT
    fi
  else
    if [ "$HOMME_Submission_Type" = lsf ]; then
      echo "mpirun.lsf -pam \"-n ${NUM_MPI_PROCS}\" ${MPI_OPTIONS} $EXEC $OPT" >> $RUN_SCRIPT
    elif [ "$HOMME_Submission_Type" = pbs ]; then
        echo "aprun -n ${NUM_MPI_PROCS} ${MPI_OPTIONS} $EXEC $OPT" >> $RUN_SCRIPT
    else
      echo "mpiexec -n ${NUM_MPI_PROCS} ${MPI_OPTIONS} $EXEC $OPT" >> $RUN_SCRIPT

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
    elif [ "${MPI_EXEC}" = "runjob" ]; then
      echo "runjob -n 1 ${MPI_OPTIONS} --block \$COBALT_PARTNAME --verbose=INFO : $EXEC" >> $RUN_SCRIPT
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



diffCprncOutput() {

  # source the test.sh file to get the names of the NC_OUTPUT_FILES
  source ${HOMME_TESTING_DIR}/${TEST_NAME}/${TEST_NAME}.sh

  # NC_OUTPUT_FILES is defined in the .sh file
  FILES="${NC_OUTPUT_FILES}"

  if [ -z "${FILES}" ] ; then
      echo "Test ${TEST_NAME} doesn't have Netcdf output files"
  fi

  # for files in movies
  exitcode=0 
  for file in $FILES 
  do
    echo "file = ${file}"
    cprncOutputFile="${HOMME_TESTING_DIR}/${TEST_NAME}/${TEST_NAME}.`basename $file`.out"
    cprncErrorFile="${HOMME_TESTING_DIR}/${TEST_NAME}/${TEST_NAME}.`basename $file`.err"

    # ensure that cprncOutputFile exists
    if [ ! -f "${cprncOutputFile}" ]; then
      echo "Error: cprnc output file ${cprncOutputFile} not found. Exiting."
      ((exitcode=exitcode-1))
    fi

    # Parse the output file to determine if they were identical
    DIFF_RESULT=`grep -ae 'diff_test' ${cprncOutputFile} | awk '{ print $8 }'`

    if [ "${DIFF_RESULT}" == IDENTICAL ] ; then
      # check for missing variables
      NUMVARS_RESULT=`grep -ae 'were not found on' ${cprncOutputFile} | awk '{ print $5 }'`
      missing=0  
      for num in $NUMVARS_RESULT
      do
          if [ "${num}" -ne 0 ] ; then
              ((missing++))
          fi
      done
      if [ "${missing}" -eq 0 ] ; then
         echo "The files are identical: DIFF_RESULT=${DIFF_RESULT} missing vars=${missing}"
      else
         echo "The files are identical: DIFF_RESULT=${DIFF_RESULT}"
         echo "But there were missing variables: =${missing}"
         ((exitcode=exitcode-10))
      fi
    else
      echo "The files are different or missing: DIFF_RESULT=${DIFF_RESULT}"
      echo "############################################################################"
      echo "CPRNC returned the following RMS differences"
      grep -a RMS ${cprncOutputFile}
      echo "############################################################################"
      ((exitcode=exitcode-10))
    fi
  done
  exit ${exitcode}
}




diffCprncRef() {
  # source the test.sh file to get the names of the NC_OUTPUT_FILES
  source ${HOMME_TESTING_DIR}/${TEST_NAME}/${TEST_NAME}.sh
  if [ -z "${TESTCASE_REF_TOL}" ]; then
    tol=0.0E0
  else
     tol=${TESTCASE_REF_TOL}
  fi

  # NC_OUTPUT_FILES is defined in the .sh file
  FILES="${NC_OUTPUT_CHECKREF}"

  if [ -z "${FILES}" ] ; then
      echo Test ${TEST_NAME}: no netcdf CHECKREF output files, skipping CHECKREF test.
  fi

  # for files in movies
  for file in $FILES 
  do
    echo "file = ${file}"
    cprncOutputFile="${HOMME_TESTING_DIR}/${TEST_NAME}/${TEST_NAME}.ref.`basename $file`.out"
    cprncErrorFile="${HOMME_TESTING_DIR}/${TEST_NAME}/${TEST_NAME}.ref.`basename $file`.err"

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
      echo  Checking RMS differences with tol = $tol 
      echo "CPRNC returned the following RMS differences"
      if grep RMS ${cprncOutputFile} ; then
        RMSnormalized=`grep RMS ${cprncOutputFile} | awk '{print $5}' -`
        for rmserror in $RMSnormalized 
        do
           check=`awk "BEGIN{ print ( $rmserror <= $tol ) }"`
           if [ $check == "1" ]; then
              echo $rmserror '<='  $tol  OK
           else
              echo $rmserror '> '  $tol  ERROR: TOL EXCEEDED
              exit -14
           fi  
        done     
      else
         echo "Error: No RMS differences found"
         exit -13
      fi
    fi
echo "############################################################################"
echo "  The diff using CPRNC has passed"
echo "############################################################################"
  done
}

checkBaselineResults() {

  filesNotFound=false

  for testFileNum in $(seq 1 ${NUM_TEST_FILES})
  do

    testFile=TEST_FILE_${testFileNum}
    source ${!testFile}

    TEST_NAME=`basename ${!testFile} .sh`

    FILES="${NC_OUTPUT_FILES}"


    if [ -z "${FILES}" ] ; then
      echo "Test ${TEST_NAME} doesn't have Netcdf output files"
    fi

    # for files in movies
    for file in $FILES 
    do
      baseFilename=`basename $file`

      # result in the repo
      repoFile=${HOMME_BASELINE_DIR}/${TEST_NAME}/movies/${baseFilename}

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


