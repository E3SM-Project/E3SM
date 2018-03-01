#!/bin/bash

# This script runs the pFUnit unit tests. It gets executed by a cron
# job and it requires one argument: the pFUnit branch to test.
# It is tailored to work on NCCS's DISCOVER machine.

# -------------------------------------------------------------------
gatherForEmail()
# -------------------------------------------------------------------
{
    echo $1 >> $EmailLog 2>&1
}

# -------------------------------------------------------------------
function notify 
# -------------------------------------------------------------------
{
   # List of recipients (need mailing list!):
   email="pfunit-regression@lists.sourceforge.net"
   if [ -e $SCR_DIR/.fail ];then
     SUBJECT="$HEADING failed"
   else
     SUBJECT="$HEADING successful!"
   fi
   if [ -e $DebugLog ]; then
      cat $DebugLog >> $EmailLog
      mail -s "$SUBJECT" $email < $DebugLog
   else
      mail -s "$SUBJECT" $email < $EmailLog
   fi
   cleanUp
   exit 0
}

# -------------------------------------------------------------------
function gitClone 
# -------------------------------------------------------------------
{
   echo "Clone sourceforge repository..."
   cd $SCR_DIR
   if [[ "$NODE" =~ discover ]]; then
     GIT="/usr/local/other/SLES11.1/git/1.8.5.2/libexec/git-core/git"   
   else
     GIT=git
   fi
   echo $GIT clone --quiet -b $BRANCH $GITREPO $BRANCH
   $GIT clone --quiet -b $BRANCH $GITREPO $BRANCH
   if [ ! -d $SCR_DIR/$BRANCH ]; then
      gatherForEmail "Not able to checkout pFUnit"
      touch $SCR_DIR/.fail
      notify 1
   fi
}


# -------------------------------------------------------------------
function handlePBS 
# -------------------------------------------------------------------
{
   echo "Submit job..."
   QSUB=sbatch

   jobScript=$HOME/bin/jobManager.sh
   if [ $USEBATCH -eq 1 ];then
     jobID=`$QSUB $jobScript | awk '{print $4}'`
   else
     $jobScript
     wait
     return 0
   fi
   if [ -z "$jobID" ]; then
      gatherForEmail "ERROR: not able to get job ID from qsub"
      touch $SCR_DIR/.fail
      notify 1
   fi

   # Max 12hr wait
   MAX_SECONDS=43200
   INCREMENT_SEC=30
   seconds=0
   done=0
   echo "Monitoring $jobID"
   while [ $seconds -lt $MAX_SECONDS ]; do
       qStatus=`qstat -a | grep $jobID | awk '{print $10}'`
       echo " - status = ($qStatus), time elapsed (seconds) = $seconds"
       if [ -z "$qStatus" ]; then
          done=1
          break
       fi
       sleep $INCREMENT_SEC
       let seconds=$seconds+$INCREMENT_SEC  
   done

   jobFile=`find . -type f -name \*$jobID\*`
   msg=''
   if [ -z "$jobFile" ]; then
      echo "FYI.  No output file from PBS"
   else
      mv $jobFile $LOG_DIR/.
      if [ "$?" -ne "0" ]; then
         gatherForEmail "ERROR: Not able to move $jobFile to $LOG_DIR"
         touch $SCR_DIR/.fail
         notify 1
      fi
   fi

   if [ $done -eq 0 ]; then
      gatherForEmail "ERROR: running out... not able to get the result from PBS script$msg"
      touch $SCR_DIR/.fail
      notify 1
   fi
}


# -------------------------------------------------------------------
function cleanUp 
# -------------------------------------------------------------------
{
   # If there were no errors then clean up scratch space
   if [ ! -e $SCR_DIR/.fail ];then
      rm -rf $SCR_DIR
   fi
}

# -------------------------------------------------------------------
# MAIN PROGRAM
# -------------------------------------------------------------------
if [ $# -lt 1 ]; then
   echo "First argument should be the branch name"
   exit 1
fi
BRANCH=$1
# If a second (optional) argument is supplied then it is used. 
# Second argument is the URL of the pFUnit git repository. Default:
GITREPO="git://git.code.sf.net/p/pfunit/code"
if [ -z "$2" ]; then
  echo "No repository supplied. Will use default"
  HEADING="pFUnit regression tests (branch: $BRANCH)"
else
  GITREPO="$2"
  HEADING="pFUnit regression tests (repository: $GITREPO)"
fi

ARCH=`uname -s`
NODE=`uname -n`
curDate=`date +"%Y_%m_%d_%m_%s"`

# Regression scripts rely on computational environment determined by
# by modules (modules.sourceforge.net). We support two machines:
#if [[ "$NODE" =~ discover || "$NODE" =~ dali ]]; then  # discover nodes
echo $NODE
MODULEINIT=/usr/share/modules/init/bash
SCRATCH=$NOBACKUP
if [[ "$NODE" =~ "ip-10" ]]; then # AWS nodes
  MODULEINIT=/usr/local/Modules/default/init/bash
  SCRATCH=/data
fi
export MODULEINIT

# Set to 0 to run interactively
USEBATCH=1

SCR_DIR=$SCRATCH/$curDate
LOG_DIR=$SCR_DIR/log
mkdir -p $SCR_DIR $LOG_DIR
export SCR_DIR=$SCR_DIR
export LOG_DIR=$LOG_DIR
export BRANCH=$BRANCH
DebugLog=$LOG_DIR/debug.log
EmailLog=$LOG_DIR/email.log

rm -f $DebugLog $EmailLog
rm -f $SCR_DIR/.fail $LOG_DIR/*

gitClone
handlePBS
notify
cleanUp

