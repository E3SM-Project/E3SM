# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  if (EAMXX_ENABLE_PYTHON)
    configure_file (${CMAKE_CURRENT_SOURCE_DIR}/my_module.py
                    ${CMAKE_CURRENT_BINARY_DIR}/my_module.py)
    CreateUnitTest(pysession "pysession.cpp")
  endif()
endif()
