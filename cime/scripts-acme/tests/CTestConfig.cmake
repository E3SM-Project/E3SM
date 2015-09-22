set(CTEST_PROJECT_NAME "ACME_Climate")
string(TIMESTAMP CURRTIME "%H:%M:%S" UTC)
set(CTEST_NIGHTLY_START_TIME "${CURRTIME} UTC")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "my.cdash.org")
set(CTEST_DROP_LOCATION "/submit.php?project=ACME_Climate")
set(CTEST_DROP_SITE_CDASH TRUE)

set(BUILDNAME "acme_scripts_regression")
