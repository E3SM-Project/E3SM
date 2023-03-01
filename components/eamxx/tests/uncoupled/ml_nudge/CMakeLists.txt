include(ScreamUtils)
set(PYBIND11_PYTHON_VERSION 3.9 CACHE STRING "")
add_subdirectory(pybind11)

find_package(Python REQUIRED COMPONENTS Interpreter Development)

set(NEED_LIBS pybind11::pybind11 Python::Python ml_nudging scream_control scream_share)
CreateUnitTest(ml_nudge_standalone "ml_nudge_standalone.cpp" "${NEED_LIBS}" LABELS "ml_nudge;physics;driver")

target_compile_definitions(ml_nudge_standalone PRIVATE -DCUSTOM_SYS_PATH="${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(ml_nudge_standalone SYSTEM PRIVATE ${PYTHON_INCLUDE_DIRS})

# Set AD configurable options
set(NUM_STEPS 1)
set(ATM_TIME_STEP 1800)

# Configure yaml input file to run directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
    ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)