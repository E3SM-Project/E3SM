# This macro configures an emulator-capable Python interpreter and uses it to generate test fixtures
# so we don't store binary files in the repo.
#
# The macro defines the following CMake variables:
# * Python3_EXECUTABLE: a Python interpreter in the appropriate virtual environment
# * EMULATOR_PYTHON_HOME: the base prefix for the Python interpreter
# * EMULATOR_TEST_FIXTURE_DIR: the directory where test fixtures are installed
# * EMULATOR_COUPLING_FIXTURES: a list of generated test fixtures for emulator coupling tests
#
# It also defines an `emulator_coupling_test_fixtures` target that can be used to express that
# an emulator test or library requires the generated fixtures.

macro(configure_emulators)

# Our emulators all require Python at the moment. And specifically, they need Python running in an
# appropriate virtual environment.

set(Python3_FIND_VIRTUALENV ONLY)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy"
  RESULT_VARIABLE EMULATOR_NUMPY_MISSING
  OUTPUT_QUIET
  ERROR_QUIET
)
execute_process(COMMAND ${Python3_EXECUTABLE} -c "import site; print(':'.join(site.getsitepackages()))"
  OUTPUT_VARIABLE EMULATOR_PYTHON_SITE_PACKAGES
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(COMMAND ${Python3_EXECUTABLE} -c "import torch"
  RESULT_VARIABLE EMULATOR_PYTHON_TORCH_MISSING
  OUTPUT_QUIET
  ERROR_QUIET
)
if(EMULATOR_PYTHON_TORCH_MISSING)
  message(FATAL_ERROR "Python interpreter (${Python3_EXECUTABLE}) does not have Torch installed.")
endif()
# A relocatable Python does not find its standard library when embedded
execute_process(COMMAND ${Python3_EXECUTABLE} -c "import sys; print(sys.base_prefix)"
  OUTPUT_VARIABLE EMULATOR_PYTHON_HOME
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Copy test fixtures into a build directory counterpart that also holds generated files.
set(EMULATOR_TEST_FIXTURE_DIR ${PROJECT_BINARY_DIR}/tests/fixtures)
file(COPY ${PROJECT_SOURCE_DIR}/tests/fixtures DESTINATION ${PROJECT_BINARY_DIR}/tests)
set(EMULATOR_COUPLING_FIXTURES
  ${EMULATOR_TEST_FIXTURE_DIR}/forcing_emulator.pt
  ${EMULATOR_TEST_FIXTURE_DIR}/land_emulator.pt
)
add_custom_command(
  OUTPUT ${EMULATOR_COUPLING_FIXTURES}
  COMMAND python ${PROJECT_SOURCE_DIR}/tests/fixtures/make_coupling_fixtures.py ${EMULATOR_TEST_FIXTURE_DIR}
  COMMENT "Generating emulator coupling test fixtures..."
)
add_custom_target(emulator_coupling_test_fixtures DEPENDS ${EMULATOR_COUPLING_FIXTURES})

endmacro()
