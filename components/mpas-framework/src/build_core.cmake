function(build_core CORE)
  set(EXE_NAME ${CORE}_model)
  set(NAMELIST_SUFFIX ${CORE})

  # Map the ESM component corresponding to each MPAS core
  if (CORE STREQUAL "ocean")
    set(COMPONENT "ocn")
  elseif(CORE STREQUAL "landice")
    set(COMPONENT "glc")
  elseif(CORE STREQUAL "seaice")
    set(COMPONENT "ice")
  else()
    message(FATAL_ERROR "Unrecognized core: ${CORE}")
  endif()

  # Gather sources
  set(CORE_BLDDIR ${CMAKE_BINARY_DIR}/core_${CORE})
  if (NOT EXISTS ${CORE_BLDDIR})
    file(MAKE_DIRECTORY ${CORE_BLDDIR})
  endif()

  set(CORE_INPUT_DIR ${CORE_BLDDIR}/default_inputs)
  if (NOT EXISTS ${CORE_INPUT_DIR})
    file(MAKE_DIRECTORY ${CORE_INPUT_DIR})
  endif()

  # Provides us RAW_SOURCES, CPPDEFS, and INCLUDES
  include(${CMAKE_CURRENT_SOURCE_DIR}/core_${CORE}/${CORE}.cmake)

  add_library(${COMPONENT})
  target_compile_definitions(${COMPONENT} PRIVATE ${CPPDEFS})
  target_include_directories(${COMPONENT} PRIVATE ${INCLUDES})
  target_link_libraries(${COMPONENT} PUBLIC ${LIBRARIES} common)

  # Make .inc files
  add_custom_command (
    OUTPUT ${CORE_BLDDIR}/Registry_processed.xml
    COMMAND cpp -P -traditional ${CPPDEFS} -Uvector
    ${CMAKE_CURRENT_SOURCE_DIR}/core_${CORE}/Registry.xml > Registry_processed.xml
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/core_${CORE}/Registry.xml
    WORKING_DIRECTORY ${CORE_BLDDIR}
    )

  set(INC_DIR ${CORE_BLDDIR}/inc)
  if (NOT EXISTS ${INC_DIR})
    file(MAKE_DIRECTORY ${INC_DIR})
  endif()

  add_custom_command(
    OUTPUT ${INC_DIR}/core_variables.inc
    COMMAND ${CMAKE_BINARY_DIR}/mpas-framework/src/tools/parse < ${CORE_BLDDIR}/Registry_processed.xml
    DEPENDS parse ${CORE_BLDDIR}/Registry_processed.xml
    WORKING_DIRECTORY ${INC_DIR}
  )

  # Disable qsmp for some files
  if (FFLAGS MATCHES ".*-qsmp.*")
    foreach(DISABLE_QSMP_FILE IN LISTS DISABLE_QSMP)
      get_filename_component(SOURCE_EXT ${DISABLE_QSMP_FILE} EXT)
      string(REPLACE "${SOURCE_EXT}" ".f90" SOURCE_F90 ${DISABLE_QSMP_FILE})
      set_property(SOURCE ${CMAKE_BINARY_DIR}/${SOURCE_F90} APPEND_STRING PROPERTY COMPILE_FLAGS " -qnosmp")
    endforeach()
  endif()

  genf90_targets("${RAW_SOURCES}" "${INCLUDES}" "${CPPDEFS}" "${NO_PREPROCESS}" "${INC_DIR}")
  target_sources(${COMPONENT} PRIVATE ${SOURCES})

endfunction(build_core)
