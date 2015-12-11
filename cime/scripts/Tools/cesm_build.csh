#! /bin/csh -f

# This wrapper script assures that modules in env_mach_specific are properly 
# loaded before compiling

set xmlquery_data=`./xmlquery -s JGFSEP MPILIB COMPILER CASEROOT`
set MPILIB=`echo $xmlquery_data | awk -F'JGFSEP' '{print $1}'`
set COMPILER=`echo $xmlquery_data | awk -F'JGFSEP' '{print $2}'`
set CASEROOT=`echo $xmlquery_data | awk -F'JGFSEP' '{print $3}'`

setenv MPILIB $MPILIB
setenv COMPILER	$COMPILER
setenv CASEROOT $CASEROOT

echo " .... determining environment variables from env_mach_specific "

source env_mach_specific || exit -1

echo " .... building model executable (calling ./Buildconf/cesm_build.pl) "

./Buildconf/cesm_build.pl $CASEROOT || exit -1

echo " .... successfully built model executable"

exit 0;

