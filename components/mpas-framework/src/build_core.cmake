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

  # The registry parser writes the core's .inc files here
  set(INC_DIR ${CORE_BLDDIR}/inc)
  if (NOT EXISTS ${INC_DIR})
    file(MAKE_DIRECTORY ${INC_DIR})
  endif()
  list(APPEND INCLUDES ${INC_DIR})

  # Provides us RAW_SOURCES, CPPDEFS, and INCLUDES
  include(${CMAKE_CURRENT_SOURCE_DIR}/core_${CORE}/${CORE}.cmake)

  add_library(${COMPONENT})
  target_compile_definitions(${COMPONENT} PRIVATE ${CPPDEFS})
  target_include_directories(${COMPONENT} PRIVATE ${INCLUDES})
  target_link_libraries(${COMPONENT} PUBLIC ${LIBRARIES} common)

  # Make .inc files
  if (MPAS_USE_DEPFILE)
    # The depfile makes the processed registry depend on every registry file it includes
    add_custom_command (
      OUTPUT ${CORE_BLDDIR}/Registry_processed.xml
      COMMAND cpp -P -traditional ${CPPDEFS} -Uvector
      -MD -MF ${CORE_BLDDIR}/Registry_processed.xml.d -MT ${CORE_BLDDIR}/Registry_processed.xml
      ${CMAKE_CURRENT_SOURCE_DIR}/core_${CORE}/Registry.xml -o ${CORE_BLDDIR}/Registry_processed.xml
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/core_${CORE}/Registry.xml
      DEPFILE ${CORE_BLDDIR}/Registry_processed.xml.d
      WORKING_DIRECTORY ${CORE_BLDDIR}
      )
  else()
    add_custom_command (
      OUTPUT ${CORE_BLDDIR}/Registry_processed.xml
      COMMAND cpp -P -traditional ${CPPDEFS} -Uvector
      ${CMAKE_CURRENT_SOURCE_DIR}/core_${CORE}/Registry.xml > Registry_processed.xml
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/core_${CORE}/Registry.xml
      WORKING_DIRECTORY ${CORE_BLDDIR}
      )
  endif()

  if (MPAS_BUILD_MODE STREQUAL "E3SM")
    # These reproduce what coupled builds have always generated: the namelist
    # suffix and executable name are empty, and the git version is the name of
    # the undefined macro. The MCT drivers set the namelist and streams names.
    set(PARSE_ARGS "--namelist-suffix=" "--exe-name=" "--git-version=MPAS_GIT_VERSION")
  endif()

  set(PARSE_OUTPUTS
    block_dimension_routines.inc
    config_declare.inc
    config_get.inc
    core_variables.inc
    define_packages.inc
    domain_variables.inc
    namelist_call.inc
    namelist_defines.inc
    setup_immutable_streams.inc
    structs_and_variables.inc
  )
  list(TRANSFORM PARSE_OUTPUTS PREPEND ${INC_DIR}/)

  add_custom_command(
    OUTPUT ${PARSE_OUTPUTS}
    COMMAND ${MPAS_PARSE} ${PARSE_ARGS} < ${CORE_BLDDIR}/Registry_processed.xml
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

  # Disable optimizations on some files that would take too long to compile, expect these to all be fortran files
  foreach (SOURCE_FILE IN LISTS NOOPT_FILES)
    get_filename_component(SOURCE_EXT ${SOURCE_FILE} EXT)
    string(REPLACE "${SOURCE_EXT}" ".f90" SOURCE_F90 ${SOURCE_FILE})
    e3sm_deoptimize_file(${CMAKE_BINARY_DIR}/${SOURCE_F90})
  endforeach()

  genf90_targets("${RAW_SOURCES}" "${INCLUDES}" "${CPPDEFS}" "${NO_PREPROCESS}" "${INC_DIR}")
  target_sources(${COMPONENT} PRIVATE ${SOURCES})

endfunction(build_core)
