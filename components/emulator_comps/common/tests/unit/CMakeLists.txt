#==============================================================================
# Common Unit Tests
#==============================================================================

message(STATUS "  Configuring common unit tests...")

# Configuration parsing tests
EmulatorUnitTest(test_emulator_config
  SOURCES test_emulator_config.cpp
  LIBS emulator_common
  LABELS config
)

# Coupling fields tests
EmulatorUnitTest(test_coupling_fields
  SOURCES test_coupling_fields.cpp
  LIBS emulator_common
  LABELS coupling
)

# Diagnostics tests
EmulatorUnitTest(test_diagnostics
  SOURCES test_diagnostics.cpp
  LIBS emulator_common
  LABELS diagnostics
)

# Output averaging tests
EmulatorUnitTest(test_output_averaging
  SOURCES test_output_averaging.cpp
  LIBS emulator_common
  LABELS diagnostics
)

# Derived diagnostics tests (pattern matching, factory)
EmulatorUnitTest(test_derived_diagnostics
  SOURCES test_derived_diagnostics.cpp
  LIBS emulator_common
  LABELS diagnostics
)

# MPI edge cases tests (single rank, uneven distribution, validation)
EmulatorUnitTest(test_mpi_edge_cases
  SOURCES test_mpi_edge_cases.cpp
  LIBS emulator_common
  LABELS mpi validation
)

# Inference backend tests
add_subdirectory(inference)
