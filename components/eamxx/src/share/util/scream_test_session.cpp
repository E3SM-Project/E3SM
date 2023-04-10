#include "share/scream_session.hpp"

/*
 * This small file contains the default implementation of a test session
 * initialization/finalization. These implementation simply call the
 * correspondinf ekat session intialization/finalization.
 *
 * If your application needs to perform additional initialization or finalization
 * work, you MUST define these functions in a cpp file. Your implementation should
 * likely call EKAT's session initialization/finalization, but are allowed to do
 * additional work.
 * When calling EkatCreateUnitTest, you must specify the option EXCLUDE_TEST_SESSION,
 * which will prevent Ekat from linking in this file (avoiding multiple definition
 * of these routines).
 */

void ekat_initialize_test_session (int argc, char** argv, const bool print_config) {
  scream::initialize_scream_session (argc,argv,print_config);
}

void ekat_finalize_test_session () {
  scream::finalize_scream_session ();
}
