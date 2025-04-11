set(MLCORRECTION_SRCS
  eamxx_ml_correction_process_interface.cpp
)

set(MLCORRECTION_HEADERS
  eamxx_ml_correction_process_interface.hpp
)
include(ScreamUtils)
    if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.11.0")
    message(STATUS "Downloading Pybind11")
    include(FetchContent)

    FetchContent_Declare(pybind11 GIT_REPOSITORY https://github.com/pybind/pybind11.git GIT_TAG v2.10.4)
    FetchContent_MakeAvailable(pybind11)
else()
    message(FATAL_ERROR "pybind11 is missing. Use CMake >= 3.11 or download it")
endif()
find_package(Python REQUIRED COMPONENTS Interpreter Development)

add_library(ml_correction ${MLCORRECTION_SRCS})
target_compile_definitions(ml_correction PUBLIC EAMXX_HAS_ML_CORRECTION)
target_compile_definitions(ml_correction PRIVATE -DML_CORRECTION_CUSTOM_PATH="${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(ml_correction SYSTEM PUBLIC ${PYTHON_INCLUDE_DIRS})
target_link_libraries(ml_correction physics_share scream_share pybind11::pybind11 Python::Python)

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE ml_correction)
endif()
