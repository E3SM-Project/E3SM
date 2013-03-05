#!/bin/bash

# This is required to properly scope the variables
exportVars() {
  export HOMME_ROOT

  export namelist_dir

  export config_dir
  export exec_dir
  export this_dir

  export exec_name
  export test_name
}

initBuild () {

  printBuildInit 

  exportVars

  namelist_dir=namelists/little_endian

  if [ `uname` == AIX ]; then
    namelist_dir=namelists/big_endian
  fi

  
  HOMME_ROOT=`(cd "$PWD/../../../../" && pwd -P)`
  echo HOMME_ROOT=$HOMME_ROOT

  # Initialize some variables 
  export HOMME_REBUILD=false
  export CONFIG_ONLY=false

  USAGE='Usage:\n  "-r" Force reconfigure/rebuild of all\n  "-c" Only configuring and copying files, no building (fast)\n  "--help" Print this message'

  while test -n "$1"; do
    case "$1" in
     "-r") 
        echo "Forcing reconfigure/rebuild of all"
        HOMME_REBUILD=true
        shift
        ;;
     "-c")
        echo "Only configuring and copying files, no building"
        CONFIG_ONLY=true
        shift
        ;;
     "--help")
        echo -e $USAGE
        exit 0
        ;;
     *)
        echo "Option $1 not supported"
        echo -e $USAGE
        exit -1
    esac
  done
}

setHommeDirs () 
{

  printSetupMsg

  config_dir=$HOMME_ROOT/build/${exec_name}
  #exec_dir=$HOMME_ROOT/build/${exec_name}/${exec_name}_${NP}_${PLEV}_${USE_PIO}
  exec_dir=$HOMME_ROOT/build/${exec_name}/${test_name}
  this_dir=$HOMME_ROOT/test/reg_test/individual_tests/ysScripts
  mkdir -p ${exec_dir}

  ## out_dir = where output will be saved
  if [ -n "$HOMME_OUT" ]; then
    out_dir=$HOMME_OUT/${exec_name}/${test_name}
    echo Found \$HOMME_OUT, output will be stored in $out_dir
  else
    out_dir=${exec_dir}
    echo Did not find \$HOMME_OUT, output will be stored in
    echo $out_dir instead
  fi

  echo "Test executable and run script being produced in"
  echo "$exec_dir"

}

addSuffixToHommeDirs()
{
  if [ -n "$1" ]; then
    exec_dir=${exec_dir}$1
    out_dir=${out_dir}$1
  else
    echo "Suffix not provided"
    exit -1
  fi
}

printPoundLine() {
  echo "#######################################################################"
}

printMessage() {
  echo "###### $1"
}

printSetupMsg() {
  
  printPoundLine
  printMessage "Starting build/rebuild of ${exec_name}/${test_name}" 
  printPoundLine

}

printTestConfigMsg() {
  
  printPoundLine
  printMessage "Querying the previous build of ${exec_name}/${test_name}" 
  printPoundLine

}

printReconfigMsg() {
  
  printPoundLine
  printMessage "Reconfiguring ${exec_name}/${test_name}"
  printPoundLine

}

printRebuildMsg() {
  
  printPoundLine
  printMessage "Building ${exec_name}/${test_name}" 
  printPoundLine

}

printBuildInit() {
  
  printPoundLine
  printMessage "Initializing build data"
  printPoundLine

}

# Not yet working
cmakeConfigure() {

  if [ "$USE_PIO" == 1 ]; then
    PIO_FLAG="-DUSE_PIO:BOOL=TRUE"
  fi

  if [ "$WITH_ENERGY" == 1 ]; then
    ENERGY_FLAG="-DENABLE_ENERGY_DIAGNOSTICS:BOOL=TRUE"
  fi

  CONFIG_COMMAND="$config_dir/cmake -DCMAKE_INSTALL_PREFIX:PATH=${exec_dir} \n\
  NP=$NP \n\
  PLEV=$PLEV \n\
  $PIO_FLAG \n\
  $ENERGY_FLAG \n\
  HOMME_ROOT"
}

autoConfCommand() {

  if [ "$USE_PIO" == 1 ]; then
    PIO_FLAG="--enable-pio"
  else
    PIO_FLAG=""
  fi

  if [ "$WITH_ENERGY" == 1 ]; then
    ENERGY_FLAG="--enable-energy-diagnostics"
  else
    ENERGY_FLAG=""
  fi

  CONFIG_COMMAND="$config_dir/configure --prefix=${exec_dir} --enable-blas --enable-lapack --with-netcdf=$NETCDF_PATH --with-pnetcdf=$PNETCDF_PATH NP=$NP PLEV=$PLEV $PIO_FLAG $ENERGY_FLAG"
}

setPIOFlag() {
  if [ "$USE_PIO" == 1 ]; then
    if "useCmake" ; then
      PIO_FLAG="-DUSE_PIO:BOOL=TRUE"
    else
      PIO_FLAG="--enable-pio"
    fi 
  fi
}

setEnergyFlag() {
  if [ "$WITH_ENERGY" == 1 ]; then
    if "useCmake" ; then
      ENERGY_FLAG="-DENABLE_ENERGY_DIAGNOSTICS:BOOL=TRUE"
    else
      ENERGY_FLAG="--enable-energy-diagnostics"
    fi
  fi
}

createConfigCommand() {

  # Set CONFIG_COMMAND
  autoConfCommand
  #cmakeConfCommand
}

resetConfOptions() {
  USE_PIO=0
  WITH_ENERGY=0
}

testConfigAndBuild ()
{

  setHommeDirs 

  # Sets CONFIG_COMMAND
  createConfigCommand

  configFile=${exec_dir}/config.h

  if [ -f $configFile -a -f ${exec_dir}/Makefile ] && ! $HOMME_REBUILD ; then
  
    printTestConfigMsg

    echo config.h exists, parsing for configuration variables
    PREV_NP=`sed -n 's/#define NP \([0-9]\)/\1/p' $configFile`
    PREV_PLEV=`sed -n 's/#define PLEV \([0-9]\)/\1/p' $configFile`
    PREV_PIO=`sed -n 's/#define PIO \([0-9]\)/\1/p' $configFile`
    PREV_ENERGY=`sed -n 's/#define ENERGY_DIAGNOSTICS \([0-9]\)/\1/p' $configFile`

    if [ "$PREV_PIO" == 1 ]; then
      PREV_USE_PIO=1
    else
      PREV_USE_PIO=0
    fi

    if [ "$PREV_ENERGY" == 1 ]; then
      PREV_WITH_ENERGY=1
    else
      PREV_WITH_ENERGY=0
    fi

    echo "Required values (NP, PLEV, USE_PIO, WITH_ENERGY) = ($NP,$PLEV,$USE_PIO,$WITH_ENERGY)"
    echo "Previous values (NP, PLEV, USE_PIO, WITH_ENERGY) = ($PREV_NP,$PREV_PLEV,$PREV_USE_PIO,$PREV_WITH_ENERGY)"

    if [ "$PREV_NP" == "$NP" -a "$PREV_PLEV" == "$PLEV" -a "$PREV_USE_PIO" == "$USE_PIO" -a \
        "$PREV_WITH_ENERGY" == "$WITH_ENERGY" -o "$PREV_WITH_ENERGY" == 0 -a -z "$WITH_ENERGY" ]; then
      echo "Previous configure matches, Skipping configure stage"
      rebuild
    else
      echo Need to reconfigure / rebuild to change NP, PLEV, PIO, or WITH_ENERGY
      reconfigRebuild
    fi
  else
    echo The file config.h or the Makefile doesn\'t exist or force rebuild, need to configure / build from start
    reconfigRebuild
  fi
  resetConfOptions
}

reconfigRebuild() {
  cd ${exec_dir}
  printReconfigMsg
  echo $CONFIG_COMMAND
  $CONFIG_COMMAND | tee config.out
  make clean
  rebuild
  cd ${this_dir}
}

rebuild() {
  cd ${exec_dir}
  # Run configure command
  printRebuildMsg
  if $CONFIG_ONLY ; then
    echo "Not building since configure only option (-c) present"
  else
    echo "Running \"make depends\"..."
    make depends > make.deps.out
    make -j 4
    #make install
  fi
  cd ${this_dir}
}

setupMPI () 
{
  if [ -f "${PBS_NODEFILE}" ]; then
     NCPU=`wc $PBS_NODEFILE | awk '{print $1}' - `
  fi
  if [ -n "${SLURM_NNODES}" ]; then
     NCPU=${SLURM_JOB_NUM_NODES}
     cores=8
     MPI_RUN="mpiexec --npernode $cores numa_wrapper --ppn $cores"
  else
    if [[ `uname -n` =~ yslogin. ]]; then
      MPI_RUN="mpirun.lsf"
    else
      MPI_RUN="mpirun -np $NCPU"
    fi
  fi
}


setupFileStructure () {

  # copy the namelist file and set a soft link
  #   a softlink will be set in the MPI run generation
  if [ -n "$nameListFiles" ]; then
    # Loop through files (there may be more than one)
    for nlFile in $nameListFiles
    do
      cp $nlFile ${out_dir}/
    done
  else
    echo "No namelist files given. Please set nameListFiles"
    exit -1
  fi

  # Make a movie directory
  mkdir -p $out_dir/movies
  rm -rf $out_dir/movies/*.nc

  mkdir -p $out_dir/restart
  rm -rf $out_dir/restart/R0000*

  if [ -n "$vcoordFiles" ]; then 
    # Create directories for output
    mkdir -p $out_dir/vcoord
    cp $vcoordFiles ${out_dir}/vcoord
  fi

  # Set soft links for reference solutions
  if [ -n "$refSolnFiles" ]; then
    for ref in $refSolnFiles
    do
      ln -sf $ref ${out_dir}/
    done
  fi

  if [ -n "$nclFiles" ]; then 
    for file in $nclFiles
    do 
      ln -sf $file ${out_dir}/`basename $file`
    done
  fi
 

  if [ -n "$OMP_SUB_TESTS" -a "$OMP_SUB_TESTS" == true ] ; then
    for nlFile in $ompNameListFiles
    do
      cp $nlFile ${out_dir}/
    done
  fi
}

yellowstoneSetupLSF () {

  RUN_SCRIPT=$out_dir/$test_name-run.sh

  #delete the file if it exists
  rm -f $RUN_SCRIPT

  # Set up some yellowstone boiler plate
  echo "#!/bin/bash" >> $RUN_SCRIPT
  echo ""  >> $RUN_SCRIPT # newlines

  echo "#BSUB -a poe" >> $RUN_SCRIPT

  # To do: move this check up and properly handle the error status
  if [ -n "$HOMME_PROJID" ]; then
    echo "#BSUB -P $HOMME_PROJID" >> $RUN_SCRIPT
  else
    echo "PROJECT CHARGE ID (HOMME_PROJID) not set"
    exit -1
  fi 

  echo "" >> $RUN_SCRIPT # newline

  echo "#BSUB -q small" >> $RUN_SCRIPT
  echo "#BSUB -W 0:20" >> $RUN_SCRIPT

  echo "" >> $RUN_SCRIPT # newline

  echo "#BSUB -x" >> $RUN_SCRIPT

  echo "" >> $RUN_SCRIPT # newline

  echo "#BSUB -R \"select[scratch_ok > 0 ]\"" >> $RUN_SCRIPT

  echo "" >> $RUN_SCRIPT # newline

  # This is the place to copy the lsf file if OMP
  if [ -n "$OMP_SUB_TESTS" -a "$OMP_SUB_TESTS" == true ] ; then
    # Copy the runscript
    OMP_RUN_SCRIPT=$out_dir/${test_name}-omp-run.sh
    cp $RUN_SCRIPT $OMP_RUN_SCRIPT
  fi

  # Set the job name
  echo "#BSUB -J $exec_name.$test_name" >> $RUN_SCRIPT
  echo "" >> $RUN_SCRIPT

  # Set the output and error filenames
  echo "#BSUB -o $test_name.stdout.%J" >> $RUN_SCRIPT
  echo "#BSUB -e $test_name.stderr.%J" >> $RUN_SCRIPT
  echo "" >> $RUN_SCRIPT

  # Set the ncpus and ranks per MPI
  echo "#BSUB -n $NCPU" >> $RUN_SCRIPT
  echo '#BSUB -R "span[ptile='$NCPU']" ' >> $RUN_SCRIPT

  echo "" >> $RUN_SCRIPT

  echo "cd $out_dir" >> $RUN_SCRIPT

  # Loop through input files and set soft links
  for nlFile in $nameListFiles
  do
    fullNlFile=${out_dir}/`basename $nlFile`
    echo "" >> $RUN_SCRIPT # blank line
    #echo "ln -sf ${out_dir}/`basename $nlFile` ${out_dir}/input.nl" >> $RUN_SCRIPT
    #echo "" >> $RUN_SCRIPT # blank line
    # Output gets automatically redirected to no redirection here
    echo "$MPI_RUN ${exec_dir}/${exec_name} < $fullNlFile" >> $RUN_SCRIPT
  done

  # Now set up the OMP script
  if [ -n "$OMP_SUB_TESTS" -a "$OMP_SUB_TESTS" == true ] ; then

    # set the job name
    echo "#BSUB -J $exec_name.$test_name-omp" >> $OMP_RUN_SCRIPT
    echo "" >> $OMP_RUN_SCRIPT

    # set the output and error filenames
    echo "#BSUB -o $test_name-omp.stdout.%J" >> $OMP_RUN_SCRIPT
    echo "#BSUB -e $test_name-omp.stderr.%J" >> $OMP_RUN_SCRIPT

    echo "" >> $OMP_RUN_SCRIPT

    calculateOMPVariables

    echo "#BSUB -n $OMP_NUM_MPI_JOBS" >> $OMP_RUN_SCRIPT
    echo '#BSUB -R "span[ptile='$OMP_NUM_MPI_JOBS']" ' >> $OMP_RUN_SCRIPT
    echo "" >> $OMP_RUN_SCRIPT # blank line
    echo "export OMP_NUM_THREADS=$OMP_NUM_THREADS" >> $OMP_RUN_SCRIPT

    for nlFile in $ompNameListFiles
    do
      echo "" >> $OMP_RUN_SCRIPT # blank line
      fullNlFile=${out_dir}/`basename $nlFile`
      #echo "ln -sf ${out_dir}/`basename $nlFile` ${out_dir}/input.nl" >> $OMP_RUN_SCRIPT
      #echo "" >> $OMP_RUN_SCRIPT # blank line
      # Output gets automatically redirected to no redirection here
      echo "$MPI_RUN ${exec_dir}/${exec_name} < $fullNlFile" >> $OMP_RUN_SCRIPT
    done
  fi

}

# Useful for debugging
printDirs () {
  echo namelist_dir = $namelist_dir

  echo Building from $HOMME_ROOT
  echo exec_dir = $exec_dir

}


calculateOMPVariables() {
  # Set the ncpus and ranks per MPI
  if  [ `expr $NCPU % $OMP_NUM_THREADS` != 0 ]; then
    echo NCPU not divisible by OMP_NUM_THREADS
    echo "NCPU = $NCPU"
    echo "OMP_NUM_THREADS = $OMP_NUM_THREADS"
    exit -1
  fi

  OMP_NUM_MPI_JOBS=`expr $NCPU / $OMP_NUM_THREADS`
}

resetVariables() {
  nameListFiles=""
  vcoordFiles=""
  nclFiles=""
  refSolnFiles=""
  ompNameListFiles=""
}

