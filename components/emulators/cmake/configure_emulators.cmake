# This macro configures an emulator-capable Python interpreter and uses it to generate test fixtures
# so we don't store binary files in the repo.
#
# The macro defines the following CMake variables:
# * Python3_EXECUTABLE: a Python interpreter in the appropriate virtual environment
# * EMULATOR_PYTHON_HOME: the base prefix for the Python interpreter
# * EMULATOR_TEST_FIXTURE_DIR: the directory where test fixtures are installed
#
# Must be called after EMULATOR_ENABLE_LIBTORCH is defined: it decides whether the .pt fixtures
# are generated (scripts/make_test_files.sh --libtorch).

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

# Populate the build-tree fixture directory (configs, Python models and,
# with LibTorch, the .pt files). Nothing generated is written to the source tree.
set(EMULATOR_TEST_FIXTURE_DIR ${PROJECT_BINARY_DIR}/tests/fixtures)
set(_emu_fixture_args --outdir ${EMULATOR_TEST_FIXTURE_DIR} --python ${Python3_EXECUTABLE})
if(EMULATOR_ENABLE_LIBTORCH)
  list(APPEND _emu_fixture_args --libtorch)
endif()
message(STATUS "Generating emulator test fixtures in ${EMULATOR_TEST_FIXTURE_DIR}")
execute_process(
  COMMAND bash ${PROJECT_SOURCE_DIR}/scripts/make_test_files.sh ${_emu_fixture_args}
  RESULT_VARIABLE _emu_fixture_result
)
if(NOT _emu_fixture_result EQUAL 0)
  message(FATAL_ERROR "scripts/make_test_files.sh failed (${_emu_fixture_result})")
endif()

# Re-run CMake (and so the script) when any fixture source changes.
file(GLOB _emu_fixture_sources CONFIGURE_DEPENDS
  ${PROJECT_SOURCE_DIR}/scripts/make_test_files.sh
  ${PROJECT_SOURCE_DIR}/scripts/config/*.py
  ${PROJECT_SOURCE_DIR}/tests/fixtures/configs/*
  ${PROJECT_SOURCE_DIR}/tests/fixtures/models/*.py
)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${_emu_fixture_sources})

endmacro()
