#
# Environment variables CIME_COMPILER and CIME_MPILIB
# can be used to send --compiler and --mpilib settings to scripts_regression_tests.py
#
set(CTEST_SITE "$CIME_MACHINE $CIME_COMPILER")
set(CTEST_PROJECT_NAME "CIME")
string(TIMESTAMP CURRTIME "%H:%M:%S" UTC)
set(CTEST_NIGHTLY_START_TIME "${CURRTIME} UTC")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "my.cdash.org")
set(CTEST_DROP_LOCATION "/submit.php?project=CIME")
set(CTEST_DROP_SITE_CDASH TRUE)

set(CTEST_TEST_TIMEOUT 82800 CACHE STRING "")
set(DART_TESTING_TIMEOUT 82800 CACHE STRING "")

set(shell $ENV{SHELL})

if (DEFINED ENV{CIME_COMPILER})
   set(compiler $ENV{CIME_COMPILER})
else()
   set(compiler "default")
endif()
if (DEFINED ENV{CIME_MPILIB})
   set(mpilib $ENV{CIME_MPILIB})
else()
   set(mpilib "default")
endif()

set(BUILDNAME "scripts_regression_${shell}_${compiler}_${mpilib}")
