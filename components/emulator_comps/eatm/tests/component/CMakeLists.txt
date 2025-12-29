#==============================================================================
# EATM Component Tests
#==============================================================================

message(STATUS "  Configuring EATM component tests...")

# EmulatorAtm component tests
EmulatorComponentTest(test_emulator_atm
  SOURCES test_emulator_atm.cpp
  LIBS eatm
  LABELS eatm
  MPI_RANKS 1 2
)
