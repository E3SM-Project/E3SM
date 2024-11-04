include(ScreamUtils)

if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.11.0")
    message(STATUS "Downloading Pybind11")
    include(FetchContent)

    FetchContent_Declare(pybind11
      GIT_REPOSITORY https://github.com/pybind/pybind11.git GIT_TAG v2.10.4)

    FetchContent_MakeAvailable(pybind11)
else()
    message(FATAL_ERROR "pybind11 is missing. Use CMake >= 3.11 or download it")
endif()

find_package(Python REQUIRED COMPONENTS Interpreter Development)

CreateUnitTest(ml_correction_standalone "ml_correction_standalone.cpp"
  LIBS pybind11::pybind11 Python::Python ml_correction scream_control scream_share
  LABELS ml_correction physics driver)

target_compile_definitions(ml_correction_standalone PRIVATE -DCUSTOM_SYS_PATH="${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(ml_correction_standalone SYSTEM PRIVATE ${PYTHON_INCLUDE_DIRS})

# Set AD configurable options
set(NUM_STEPS 1)
set(ATM_TIME_STEP 1800)
set (RUN_T0 2021-10-12-45000)

# Configure yaml input file to run directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
