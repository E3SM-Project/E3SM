#==============================================================================
# Inference Backend Unit Tests
#==============================================================================

message(STATUS "    Configuring inference backend tests...")

# Stub backend tests
EmulatorUnitTest(test_stub_backend
  SOURCES test_stub_backend.cpp
  LIBS emulator_common
  LABELS inference
)

# Inference factory tests
EmulatorUnitTest(test_inference_factory
  SOURCES test_inference_factory.cpp
  LIBS emulator_common
  LABELS inference
)

# LibTorch backend tests (optional, requires LibTorch)
if (EMULATOR_HAS_LIBTORCH)
  EmulatorUnitTest(test_libtorch_backend
    SOURCES test_libtorch_backend.cpp
    LIBS emulator_common
    LABELS inference libtorch
  )
  message(STATUS "      Added LibTorch backend test")

endif()
