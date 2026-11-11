# System tests: tests that run a Homme executable (from test/execs) with some namelists.
# Each folder contains the tests of a model (or of a group of test cases), and their namelists.
# The tests run in <build-dir>/tests/<test-name>. See README for a description of the tests.

# The testing profiles to use for the tests created for each profile (see homme_set_profile_params)
homme_get_profiles_up_to (${HOMME_TESTING_PROFILE} HOMME_TESTING_PROFILES)

add_subdirectory (sweqx)
add_subdirectory (preqx)
add_subdirectory (theta-l/baroclinic)
add_subdirectory (theta-l/held-suarez)
add_subdirectory (theta-l/dcmip2012)
add_subdirectory (theta-l/nhgw)
add_subdirectory (theta-l/bubble)
add_subdirectory (theta-l/transport)
add_subdirectory (tool)
