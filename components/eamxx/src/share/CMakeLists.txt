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

# Property checks structures
add_subdirectory(property_checks)

# Remapper structures
add_subdirectory(remap)

# Data managers structures
add_subdirectory(data_managers)

# Atm process
add_subdirectory(atm_process)

# Diagnostic classes
add_subdirectory(diagnostics)

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
