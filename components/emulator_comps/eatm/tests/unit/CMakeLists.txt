#==============================================================================
# EATM Unit Tests
#==============================================================================

message(STATUS "  Configuring EATM unit tests...")

# AtmFieldManager tests
EmulatorUnitTest(test_atm_field_manager
  SOURCES test_atm_field_manager.cpp
  LIBS eatm
  LABELS eatm fields
)

# AtmCoupling tests
EmulatorUnitTest(test_atm_coupling
  SOURCES test_atm_coupling.cpp
  LIBS eatm
  LABELS eatm coupling
)
