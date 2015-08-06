#! /bin/csh -f

# This wrapper script assures that modules in env_mach_specific are properly 
# loaded before compiling

setenv MPILIB     `./xmlquery  MPILIB   -value `
setenv COMPILER	  `./xmlquery  COMPILER	-value `
setenv CASEROOT   `./xmlquery  CASEROOT	-value `

echo " .... determining environment variables from env_mach_specific "

source env_mach_specific || exit -1

echo " .... building model executable (calling ./Buildconf/cesm_build.pl) "

./Buildconf/cesm_build.pl $CASEROOT || exit -1

echo " .... successfully built model executable"

exit 0;

