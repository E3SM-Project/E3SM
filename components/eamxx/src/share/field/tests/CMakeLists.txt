# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  # Test fields
  CreateUnitTest(field "field_tests.cpp" LIBS eamxx_field)

  # Test fields
  CreateUnitTest(subfield "subfield_tests.cpp" LIBS eamxx_field)

  # Test field utils
  CreateUnitTest(field_utils "field_utils.cpp" LIBS eamxx_field
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS})

  if (EAMXX_ENABLE_PYTHON)
    CreateUnitTest(pyfield "pyfield.cpp" LIBS eamxx_field)
  endif()
endif()
