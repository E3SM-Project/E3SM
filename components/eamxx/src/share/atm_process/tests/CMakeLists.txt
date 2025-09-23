# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  # Test atmosphere processes
  CreateUnitTest(atm_proc "atm_proc_tests.cpp" LIBS eamxx_atm_process)

  # Test atm proc DAG
  CreateUnitTest(atm_proc_dag "atm_proc_dag_tests.cpp" LIBS eamxx_atm_process)

  # Test diag interfaces
  CreateUnitTest(atm_diag "atm_diag_tests.cpp" LIBS eamxx_atm_process)
endif()
