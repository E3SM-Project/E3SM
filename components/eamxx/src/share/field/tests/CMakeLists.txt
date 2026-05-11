# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  # Test field header metadata
  CreateUnitTest(field_header
    SOURCES field_header_tests.cpp
    LIBS eamxx_field
  )

  # Test field class
  CreateUnitTest(field
    SOURCES field_tests.cpp
    LIBS eamxx_field
  )

  # Test field update methods
  CreateUnitTest(field_update
    SOURCES field_update_tests.cpp
    LIBS eamxx_field
  )

  # Test masked field operations
  CreateUnitTest(field_masked_ops
    SOURCES field_masked_ops_tests.cpp
    LIBS eamxx_field
  )

  # Test field groups
  CreateUnitTest(field_group
    SOURCES field_group_tests.cpp
    LIBS eamxx_field
  )

  # Test subfields
  CreateUnitTest(subfield
    SOURCES subfield_tests.cpp
    LIBS eamxx_field
  )

  # Test compute mask
  CreateUnitTest(compute_mask
    SOURCES compute_mask.cpp
    LIBS eamxx_field
  )

  # Test field reductions
  CreateUnitTest(reductions
    SOURCES reductions.cpp
    LIBS eamxx_field
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  )

  # Test field contractions
  CreateUnitTest(contractions
    SOURCES contractions.cpp
    LIBS eamxx_field
  )

  # Test field utils
  CreateUnitTest(field_utils
    SOURCES field_utils.cpp
    LIBS eamxx_field
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  )

  if (EAMXX_ENABLE_PYTHON)
    CreateUnitTest(pyfield
      SOURCES pyfield.cpp
      LIBS eamxx_field
    )
  endif()
endif()
