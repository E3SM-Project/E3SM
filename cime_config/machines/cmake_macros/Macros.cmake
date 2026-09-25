#
# Use this file to include the relevant macros based on
# machine/compiler settings. This file gets copied to CASEROOT
# and that's the one that gets included by the build system. Feel free
# to modify this file in the CASEROOT.
#
cmake_policy(SET CMP0057 NEW)

set(MACROS_DIR ${CASEROOT}/cmake_macros)

set(UNIVERSAL_MACRO ${MACROS_DIR}/universal.cmake)
set(COMPILER_MACRO ${MACROS_DIR}/${COMPILER}.cmake)
set(MACHINE_MACRO ${MACROS_DIR}/${MACH}.cmake)
set(COMPILER_MACHINE_MACRO ${MACROS_DIR}/${MACH}_${COMPILER}.cmake)
set(POST_PROCESS_MACRO ${SRCROOT}/cime_config/machines/cmake_macros/post_process.cmake)

if (CONVERT_TO_MAKE)
  get_cmake_property(E3SM_CMAKE_INTERNAL_VARS_BEFORE_BUILD_INTERNAL_IGNORE VARIABLES)
  foreach (VAR_BEFORE IN LISTS E3SM_CMAKE_INTERNAL_VARS_BEFORE_BUILD_INTERNAL_IGNORE)
    set("E3SM_CMAKE_INTERNAL_${VAR_BEFORE}" "${${VAR_BEFORE}}")
  endforeach()
  list(APPEND E3SM_CMAKE_INTERNAL_VARS_BEFORE_BUILD_INTERNAL_IGNORE "VAR_BEFORE")
  list(APPEND E3SM_CMAKE_INTERNAL_VARS_BEFORE_BUILD_INTERNAL_IGNORE "MACRO_FILE")
endif()

# Include order defines precedence
foreach (MACRO_FILE ${UNIVERSAL_MACRO} ${COMPILER_MACRO} ${MACHINE_MACRO} ${COMPILER_MACHINE_MACRO} ${POST_PROCESS_MACRO})
  if (EXISTS ${MACRO_FILE})
    include(${MACRO_FILE})
  else()
    message("No macro file found: ${MACRO_FILE}")
  endif()
endforeach()

# When building the eamxx component, also source the eamxx-specific machine
# file (if one exists) so that eamxx picks up settings such as EKAT/Kokkos
# includes and MPI launcher configuration. These are appended AFTER the CIME
# macros so eamxx uses the E3SM-provided compiler and flag settings by default,
# with eamxx-specific additions layered on top. Prefer a machine+compiler file
# over a plain machine file, and fall back to common.cmake if neither exists.
if (COMP_NAME STREQUAL "eamxx")
  set(EAMXX_MACH_FILES_DIR ${SRCROOT}/components/eamxx/cmake/machine-files)
  set(EAMXX_MACHINE_COMPILER_MACRO ${EAMXX_MACH_FILES_DIR}/${MACH}-${COMPILER}.cmake)
  set(EAMXX_MACHINE_MACRO          ${EAMXX_MACH_FILES_DIR}/${MACH}.cmake)
  set(EAMXX_COMMON_MACRO           ${EAMXX_MACH_FILES_DIR}/common.cmake)
  if (EXISTS ${EAMXX_MACHINE_COMPILER_MACRO})
    include(${EAMXX_MACHINE_COMPILER_MACRO})
  elseif (EXISTS ${EAMXX_MACHINE_MACRO})
    include(${EAMXX_MACHINE_MACRO})
  elseif (EXISTS ${EAMXX_COMMON_MACRO})
    include(${EAMXX_COMMON_MACRO})
    common_setup()
  endif()
endif()

if (CONVERT_TO_MAKE)
  get_cmake_property(VARS_AFTER VARIABLES)

  foreach (VAR_AFTER IN LISTS VARS_AFTER)
    if (VAR_AFTER MATCHES "^E3SM_CMAKE_INTERNAL_")
      # skip
    else()
      if (NOT VAR_AFTER IN_LIST E3SM_CMAKE_INTERNAL_VARS_BEFORE_BUILD_INTERNAL_IGNORE)
        message("CIME_SET_MAKEFILE_VAR ${VAR_AFTER} := ${${VAR_AFTER}}")
        list(APPEND E3SM_CMAKE_INTERNAL_VARS_BEFORE_BUILD_INTERNAL_IGNORE "${VAR_AFTER}")
        set("E3SM_CMAKE_INTERNAL_${VAR_AFTER}" "${${VAR_AFTER}}")
      elseif (NOT "${${VAR_AFTER}}" STREQUAL "${E3SM_CMAKE_INTERNAL_${VAR_AFTER}}")
        message("CIME_SET_MAKEFILE_VAR ${VAR_AFTER} := ${${VAR_AFTER}}")
        set("E3SM_CMAKE_INTERNAL_${VAR_AFTER}" "${${VAR_AFTER}}")
      endif()
    endif()
  endforeach()
endif()
