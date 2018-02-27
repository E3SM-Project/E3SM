## This file should be placed in the root directory of your cmake project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
##
## # The following are required to submit to the CDash dashboard:
##   ENABLE_TESTING()
##   INCLUDE(CTest)
#
# No modification should be needed for specific machines

set(CTEST_PROJECT_NAME "ACME_Climate")
set(CTEST_NIGHTLY_START_TIME "00:00:00 EST")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "my.cdash.org")
set(CTEST_DROP_LOCATION "/submit.php?project=ACME_Climate")
set(CTEST_DROP_SITE_CDASH TRUE)
