# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)
  # Test atmosphere processes
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/atm_process_tests_named_procs.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/atm_process_tests_named_procs.yaml COPYONLY)

  CreateUnitTest(atm_proc "atm_process_tests.cpp" LIBS eamxx_atm_process)
endif()
