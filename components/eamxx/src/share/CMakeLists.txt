# Core types and configuration options
add_subdirectory(core)

# Low-level scorpio interfaces
add_subdirectory(scorpio_interface)

# Small/generic utilities
add_subdirectory(util)

# Basic physics constant/functions
add_subdirectory(physics)

# Field structures
add_subdirectory(field)

# Grid structures
add_subdirectory(grid)

# Diagnostic classes
add_subdirectory(diagnostics)

# Property checks structures
add_subdirectory(property_checks)

# Remapper structures
add_subdirectory(remap)

# Data managers structures
add_subdirectory(data_managers)

# The python inference backend of components/emulators, for processes that run ML
# models. It embeds the python of Python3_EXECUTABLE (or Python3_ROOT_DIR).
if (EAMXX_ENABLE_INFERENCE AND NOT TARGET emulator_inference)
  set (EMULATOR_ENABLE_PYTHON ON CACHE BOOL "Build the Python inference backend" FORCE)
  add_subdirectory(${SCREAM_BASE_DIR}/../emulators/common/src/inference
                   ${CMAKE_BINARY_DIR}/emulator_inference)
  # Where that python has its standard library and its packages: when embedded, a
  # relocatable python does not see the former by itself, nor a virtual environment
  # the latter. Processes give them to the backend.
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  execute_process(COMMAND ${Python3_EXECUTABLE} -c "import sys; print(sys.base_prefix)"
    OUTPUT_VARIABLE PYTHON_HOME OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${Python3_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_paths()['purelib'])"
    OUTPUT_VARIABLE PYTHON_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
  set (EAMXX_INFERENCE_PYTHON ${Python3_EXECUTABLE} CACHE INTERNAL "Python the inference backend embeds")
  set (EAMXX_INFERENCE_PYTHON_HOME ${PYTHON_HOME} CACHE INTERNAL "Its home")
  set (EAMXX_INFERENCE_PYTHON_PACKAGES ${PYTHON_PACKAGES} CACHE INTERNAL "Its packages")
endif()

# Atm process
add_subdirectory(atm_process)

# IO library
add_subdirectory(io)

# Algorithm library
add_subdirectory(algorithm)

# For bwd compatibility: a "catch" all interface lib of all the stuff in the share folder
add_library(scream_share INTERFACE)
target_link_libraries(scream_share INTERFACE
  eamxx_core
  eamxx_scorpio_interface
  eamxx_utils
  eamxx_field
  eamxx_grid
  eamxx_remap
  eamxx_property_checks
  eamxx_data_managers
  eamxx_atm_process
  eamxx_diagnostics
  eamxx_io
  eamxx_algorithm)

# The "build_cf_dictionary" target downloads an XML file containing valid field
# names and aliases based on the CF conventions, and transforms it into a YAML
# file that we use to validate our own fields against it. This target is only
# available if CMake finds wget on your system.
find_program(WGET wget)
if (NOT WGET STREQUAL "WGET-NOTFOUND")
  set(CF_CONVENTIONS_XML_FILE cf-standard-name-table.xml)
  set(CF_CONVENTIONS_URL "http://cfconventions.org/Data/cf-standard-names/76/src/${CF_CONVENTIONS_XML_FILE}")
  add_custom_target(build_cf_dictionary
    ${WGET} ${CF_CONVENTIONS_URL}
    COMMAND ${PROJECT_SOURCE_DIR}/scripts/cf-xml-to-yaml ${CF_CONVENTIONS_XML_FILE}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/data
    COMMENT "Downloading CF conventions XML file and building YAML dictionary...")
endif()
