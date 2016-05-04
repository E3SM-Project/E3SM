#!/bin/csh -f

# Make copies of the namelist files for each part of the test. Enclose the
# copies in conditionals so that we only do this namelist setup the first time
# the build script is invoked - otherwise, if the build is rerun, the namelist
# files would build up repeated instances of the setting of force_init_intep.
#
# Note the use of shell wildcards to make sure we apply these mods to
# multi-instance versions

if ( ! -e user_nl_nointerp ) then
    mkdir user_nl_nointerp
    cp user_nl_clm* user_nl_nointerp
    foreach file (user_nl_nointerp/user_nl_clm*)
	# False is the default, but set it explicitly to be sure
    	echo "use_init_interp = .false." >> $file
    end
endif

if ( ! -e user_nl_interp ) then
    mkdir user_nl_interp
    cp user_nl_clm* user_nl_interp
    foreach file (user_nl_interp/user_nl_clm*)
    	echo "use_init_interp = .true." >> $file
    end
endif


# The actual build for this test is the same as for the standard test, so the
# contents of test_build.csh are copied below:
setenv CIMEROOT `./xmlquery CIMEROOT    -value`

./Tools/check_lockedfiles || exit -1

# NOTE - Are assumming that are already in $CASEROOT here
set CASE     = `./xmlquery CASE    -value`

./case.build -testmode
if ($status != 0) then
   echo "Error: build failed" >! ./TestStatus
   echo "CFAIL $CASE" > ./TestStatus
   exit -1    
endif 

