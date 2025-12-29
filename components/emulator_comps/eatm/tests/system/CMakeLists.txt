#==============================================================================
# EATM System Tests
#==============================================================================

message(STATUS "  Configuring EATM system tests...")

# Full EATM standalone run
EmulatorSystemTest(test_eatm_standalone
  SOURCES test_eatm_standalone.cpp
  LIBS eatm
  LABELS eatm system
  MPI_RANKS 1 4
)

# Multi-timestep evolution test
EmulatorSystemTest(test_multi_timestep
  SOURCES test_multi_timestep.cpp
  LIBS eatm
  LABELS eatm system
  MPI_RANKS 1 2
)
