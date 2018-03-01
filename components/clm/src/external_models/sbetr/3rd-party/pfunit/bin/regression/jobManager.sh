#!/bin/bash
#SBATCH --account=k3002
#SBATCH --job-name=pFUnit
#SBATCH --time=02:00:00
#SBATCH --ntasks=16
#SBATCH --constraint=hasw

# This script manages the jobs in a batch environment.
# It gets called from mainRegress.sh.
# It is tailored to work on NCCS's DISCOVER machine.

OK=0
ERR=1
umask 022

function abortNotify {
   if [ $# -eq 1 ]; then
      if [ -e $EmailLog ]; then
         cat $EmailLog >> $DebugLog
      fi
      echo -e "$1" >> $DebugLog
   fi 
   exit -1
}


function setModule {
   if [ $# -lt 3 ]; then
      abortNotify "Invalid number of arguments ($#) in function 'setModule'"
   fi

   local fortranCompiler=$1
   local version=$2
   local parallel=$3

   . "$MODULEINIT"
   if [ $? -ne 0 ]; then
      abortNotify "Problem with starting up the module environment"
   fi

   moduleList=''
   moduleMPI=''
   
   export PATH=/usr/local/other/SLES11/SIVO-PyD/1.9.0/bin:$PATH
   if [ "$fortranCompiler" == "INTEL" ]; then
      if [ "$version" == "13.1" ]; then
         moduleFortran='comp/intel-13.1.3.192'
      elif [ "$version" == "13.0" ]; then
         moduleFortran='comp/intel-13.0.1.117'
      elif [ "$version" == "14.0" ]; then
         moduleFortran='comp/intel-14.0.2.144'
      elif [ "$version" == "15.0" ]; then
         moduleFortran='comp/intel-15.0.3.187'
      else
         msg="$fortranCompiler version $version is not supported yet"
         echo -e "$msg\n\n" >> $DebugLog
      fi
      if [[ "$parallel" == "mpi" || "$parallel" == "hybrid" ]]; then
         moduleMPI=' mpi/impi-5.0.3.048'
      fi
   elif [ "$fortranCompiler" == "PGI" ]; then
      moduleFortran='comp/pgi-15.1.0'
      if [[ "$parallel" == "mpi"  || "$parallel" == "hybrid" ]]; then
         moduleMPI=' other/mpi/openmpi/1.8.1-pgi-15.1.0'
      fi
   elif [ "$fortranCompiler" == "NAG" ]; then
      moduleFortran='comp/nag-6.0'
      if [[ "$parallel" == "mpi"  || "$parallel" == "hybrid" ]]; then
         moduleMPI=' other/mpi/openmpi/1.8.1-nag-6.0'
      fi
   elif [ "$fortranCompiler" == "GNU" ]; then
      if [ "$version" == "4.9.1" ]; then
         moduleFortran='other/comp/gcc-4.9.1'
         if [[ "$parallel" == "mpi"  || "$parallel" == "hybrid" ]]; then
           moduleMPI=' other/mpi/openmpi/1.8.1-gcc-4.9.1'
         fi
      elif [ "$version" == "4.8.1" ]; then
         moduleFortran='other/comp/gcc-4.8.1'
         if [[ "$parallel" == "mpi"  || "$parallel" == "hybrid" ]]; then
           moduleMPI=' other/mpi/openmpi/1.7.2-gcc-4.8.1-shared'
         fi
      else
         msg="$fortranCompiler version $version is not supported yet"
         echo -e "$msg\n\n" >> $DebugLog
      fi
   else
      abortNotify "$fortranCompiler is not supported yet"
   fi
   moduleList="other/cmake-2.8.11.2 $moduleFortran $moduleMPI"

   module purge
   if [ -n "$moduleList" ]; then  
      module load $moduleList
      if [ $? -ne 0 ]; then
         abortNotify "Failure loading "$moduleList
      fi
   fi
   echo " -- module list: "$moduleList
}

function doMake {
   if [ $# -lt 4 ]; then
      abortNotify "Invalid number of arguments ($#) in function 'doCmake'"
   fi
   local COM=$1
   local VER=$2
   local PAR=$3
   local MAK=$4
   local BUILDOK
   local RUNOK

   # Avoid repeating some MPI combinations
   if [[ "$PAR" == "mpi"  || "$PAR" == "hybrid" ]]; then
      if [[ "$VER" == "13.0" && "$COM" == "INTEL" ]]; then
         return 0
      fi
      if [[ "$VER" == "13.1" && "$COM" == "INTEL" ]]; then
         return 0
      fi
   fi
   if [ "$BRANCH" == "pfunit_2.1.0" ]; then
      if [ "$VER" == "4.8.1" ]; then
         return 0
      fi
   fi

   # Default
   USEMPI="NO"
   USEOPENMP="NO"
   if [ "$PAR" == "omp" ]; then
      USEOPENMP="YES"
   fi
   if [ "$PAR" == "mpi" ]; then
      USEMPI="YES"
   fi
   if [ "$PAR" == "hybrid" ]; then
      USEMPI="YES"
      USEOPENMP="YES"
   fi

   makeLog=$LOG_DIR/${MAK}_${COM}_${VER}_PAR-${PAR}.log
   makeExLog=$LOG_DIR/${MAK}_${COM}_${VER}_PAR-${PAR}_Ex.log

   MAKE=/usr/bin/make

   if [ "$MAK" == "cmake" ]; then
     mkdir -p ${COM}_${VER}_PAR-${PAR}; cd ${COM}_${VER}_PAR-${PAR}
     echo " -- cmake -DMPI=$USEMPI -DOPENMP=$USEOPENMP ../"
     cmake -DMPI=$USEMPI -DOPENMP=$USEOPENMP ../ 1> $makeLog 2>&1
     echo " -- $MAKE all"
     $MAKE -j all 1>> $makeLog 2>&1
     echo " - Parse for build errors..." 
     buildErrorFile "$makeLog"
     BUILDOK=$?
     echo " --- bld RC = $BUILDOK"
     if [ "$BUILDOK" == "$OK" ]; then
       echo " -- $MAKE tests"
       $MAKE -j tests 1>> $makeLog 2>&1
       echo " - Parse for runtime errors..." 
       runError "$makeLog"
       RUNOK=$?
       echo " --- runError RC = $RUNOK"
       if [ "$RUNOK" == "-1" ]; then
         results[1]="Run_err"
       elif [ "$RUNOK" == "-2" ]; then
         results[1]="Tst_err"
       elif [ "$RUNOK" == "255" ]; then
         results[1]="Run_err"
       fi
     else
       results[1]="Bld_err"
       results[2]="na"
     fi
   else
     $MAKE --quiet distclean F90_VENDOR=$COM 1> $makeLog 2>&1
     echo " -- $MAKE all F90_VENDOR=$COM MPI=$USEMPI OPENMP=$USEOPENMP"
     $MAKE -j all F90_VENDOR=$COM MPI=$USEMPI OPENMP=$USEOPENMP 1> $makeLog 2>&1
     echo " - Parse for build errors..." 
     buildErrorFile "$makeLog"
     BUILDOK=$?
     echo " --- bld RC = $BUILDOK"
     if [ "$BUILDOK" == "$OK" ]; then
       echo " -- $MAKE tests F90_VENDOR=$COM MPI=$USEMPI OPENMP=$USEOPENMP"
       $MAKE -j tests F90_VENDOR=$COM MPI=$USEMPI OPENMP=$USEOPENMP 1> $makeLog 2>&1
       echo " - Parse for runtime errors..." 
       runError "$makeLog"
       RUNOK=$?
       echo " --- runError RC = $RUNOK"
       if [ "$RUNOK" == "-1" ]; then
         results[1]="Run_err"
       elif [ "$RUNOK" == "-2" ]; then
         results[1]="Tst_err"
       elif [ "$RUNOK" == "255" ]; then
         results[1]="Run_err"
       fi
     else
       results[1]="Bld_err"
       results[2]="na"
     fi
     if [ "$BUILDOK" == "$OK" ]; then
       # Test examples
       echo " - Test examples..."
       export PFUNIT=$SCR_DIR/pFUnit_${COM}_${VER}_PAR-${PAR}
       make install INSTALL_DIR=$PFUNIT 1>> $makeLog 2>&1
       cd $SCR_DIR/$BRANCH/Examples
       $MAKE --quiet clean
       echo " -- $MAKE all MPI=$USEMPI SKIP_INTENTIONALLY_BROKEN=YES"
       $MAKE all MPI=$USEMPI OPENMP=$USEOPENMP SKIP_INTENTIONALLY_BROKEN=YES 1>> $makeExLog 2>&1
       echo " - Parse for build errors..." 
       buildErrorFile "$makeExLog"
       BUILDOK=$?
       echo " --- bld RC = $BUILDOK"
       if [ $BUILDOK == "$ERR" ]; then
         results[2]="Bld_err"
       else
         echo " - Parse for runtime errors..." 
         runError "$makeExLog"
         RUNOK=$?
         echo " --- runError RC = $RUNOK"
         if [ $RUNOK != "$OK" ]; then
           results[2]="Run_err"
         fi
       fi
       cd -
     fi
   fi
   cd $SCR_DIR/$BRANCH

}

buildErrorFile() {
   local makeLog=$1
   local noExe=`cat $makeLog | grep tests.x | grep 'Command not found'`
   local intError=`grep '**Internal compiler error' $makeLog`
   if [[ "$rc" == "" && "$intError" == "" ]]; then
     echo " --- Build phase OK" 
     return $OK
   else
     echo " --- tests.x not found" 
     touch $SCR_DIR/.fail
     cd $SCR_DIR/$BRANCH
     return $ERR
   fi
}
runError() {
   local makeLog=$1
   local anyError=`grep 'make: ***' $makeLog`
   local failMsg=`grep Failures $makeLog`
   local segMsg=`grep SIGSEG $makeLog`
   local errorMsg=`grep Errors $makeLog`
   if [[ "$anyError" == "" && "$failMsg" == "" && "$errorMsg" == "" && "$segMsg" == "" ]]; then
     echo " --- Run phase OK" 
     return $OK
   else
     touch $SCR_DIR/.fail
     cd $SCR_DIR/$BRANCH
     if [[ "$anyError" != "" || "$segMsg" != "" ]]; then
       echo " --- runError: runtime error"
       return -1
     elif [[ "$failMsg" != "" || "$errorMsg" != "" ]]; then
       echo " --- runError: some tests failed" 
       return -2
     else
       echo " --- runError: unrecognized error"
       return 1
     fi
   fi
}

# -------------------------------------------------------------------
createEmailReport()
# -------------------------------------------------------------------
{
  local numLines=${#report[*]}
  local i=0
  while [ $i -lt $numLines ]; do
     echo "${report[$i]}" 
     echo "${report[$i]}" >> $LOG_DIR/.report
     let i++
  done

echo "pFUnit test results, branch=$BRANCH" >> $LOG_DIR/.foo
echo "------------------------------------------------------------" >> $LOG_DIR/.foo
echo "Make    ""Compiler  ""Version   ""Parallel  " "|" " Code    ""Examples" >> $LOG_DIR/.foo
echo "------------------------------------------------------------" >> $LOG_DIR/.foo
awk '{ printf "%-8s%-10s%-10s%-10s %s  %-8s%-8s \n", $1, $2, $3, $4, $5, $6, $7}' $LOG_DIR/.report  >> $LOG_DIR/.foo

   mv $LOG_DIR/.foo $EmailLog
}

# -------------------------------------------------------------------
# MAIN PROGRAM
# -------------------------------------------------------------------


DebugLog=$LOG_DIR/debug.log
EmailLog=$LOG_DIR/email.log

declare -a F90_VERSIONS
declare -a GNU_VERSIONS
declare -a INTEL_VERSIONS

# We support two build types
MAKE_TYPE=(cmake gmake)

# Three compilers
COMPILERS=(INTEL GNU NAG PGI)

# Compiler versions are separated into those that work
# with pfunit_2.1.0 and those that work with the master,
# development and release-3.0 branches.
INTEL_VERSIONS_master=(13.1 14.0 15.0)
INTEL_VERSIONS_2_1_0=(13.0 13.1 14.0)
GNU_VERSIONS_master=(4.9.1)
GNU_VERSIONS_2_1_0=(4.8.1 4.9.1)
NAG_VERSIONS=(6.0)
PGI_VERSIONS=(15.1.0)

# serial, mpi, omp, omp+mpi
PARALLEL=(off mpi omp hybrid)

curDate=`date +"%Y_%m_%d"`

cd $SCR_DIR/$BRANCH

# BRANCH is an input argument to mainRegress.sh
if [[ "$BRANCH" == "pfunit_2.1.0" ]]; then
  GNU_VERSIONS=( "${GNU_VERSIONS_2_1_0[@]}" )
  INTEL_VERSIONS=( "${INTEL_VERSIONS_2_1_0[@]}" )
else
  GNU_VERSIONS=( "${GNU_VERSIONS_master[@]}" )
  INTEL_VERSIONS=( "${INTEL_VERSIONS_master[@]}" )
fi

declare -a results
declare -a report

for eachMake in "${MAKE_TYPE[@]}"; do

  for eachFC in "${COMPILERS[@]}"; do
    if [ "$eachFC" == "INTEL" ]; then
      F90_VERSIONS=( "${INTEL_VERSIONS[@]}" )
    elif [ "$eachFC" == "GNU" ]; then 
      F90_VERSIONS=( "${GNU_VERSIONS[@]}" )
    elif [ "$eachFC" == "NAG" ]; then 
      F90_VERSIONS=( "${NAG_VERSIONS[@]}" )
    elif [ "$eachFC" == "PGI" ]; then 
      F90_VERSIONS=( "${PGI_VERSIONS[@]}" )
    fi

    for version in "${F90_VERSIONS[@]}"; do

      for eachPara in "${PARALLEL[@]}"; do
        logFile=$LOG_DIR/pFUnit_${eachFC}_${version}_MPI-${eachPara}.log
        echo " - TEST: $eachMake with Compiler=$eachFC, Version=$version, Parallel=$eachPara"
        # Initialize results
        unitTests=OK
        examples=OK
        [ $eachMake == "cmake" ] && examples=na || examples=OK
        results=("|" $unitTests $examples)
        setModule $eachFC $version $eachPara
        doMake $eachFC $version $eachPara $eachMake
        # Update report
        resultString="$eachMake $eachFC $version $eachPara ${results[@]}"
        report=( "${report[@]}" "$resultString" )

      done # eachPara

    done # eachVersion
    F90_VERSIONS=( )

  done #eachCompiler

done #eachMake

createEmailReport

exit 0
