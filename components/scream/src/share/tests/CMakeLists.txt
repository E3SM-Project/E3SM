include(ScreamUtils)

# Test fields
CreateUnitTest(field "field_tests.cpp" scream_share)

# Test atmosphere processes
CreateUnitTest(atm_proc "atm_process_tests.cpp" scream_share)
