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
  endif()

  # build_options.mk stuff handled here
  if (CORE STREQUAL "ocean")
    list(APPEND CPPDEFS "-DCORE_OCEAN")
    list(APPEND INCLUDES "${CMAKE_BINARY_DIR}/core_ocean/BGC" "${CMAKE_BINARY_DIR}/core_ocean/shared" "${CMAKE_BINARY_DIR}/core_ocean/analysis_members" "${CMAKE_BINARY_DIR}/core_ocean/cvmix" "${CMAKE_BINARY_DIR}/core_ocean/mode_forward" "${CMAKE_BINARY_DIR}/core_ocean/mode_analysis" "${CMAKE_BINARY_DIR}/core_ocean/mode_init")

  elseif (CORE STREQUAL "seaice")
    list(APPEND CPPDEFS "-DCORE_SEAICE" "-Dcoupled" "-DCCSMCOUPLED")
    list(APPEND INCLUDES "${CMAKE_BINARY_DIR}/core_seaice/column" "${CMAKE_BINARY_DIR}/core_seaice/shared" "${CMAKE_BINARY_DIR}/core_seaice/analysis_members" "${CMAKE_BINARY_DIR}/core_seaice/model_forward")

  elseif (CORE STREQUAL "landice")
    list(APPEND CPPDEFS "-DCORE_LANDICE")
    list(APPEND INCLUDES "${CMAKE_BINARY_DIR}/core_landice/shared" "${CMAKE_BINARY_DIR}/core_landice/analysis_members" "${CMAKE_BINARY_DIR}/core_landice/mode_forward")

    #
    # Check if building with LifeV, Albany, and/or PHG external libraries
    #

    if (LIFEV)
      # LifeV can solve L1L2 or FO
      list(APPEND CPPDEFS "-DLIFEV" "-DUSE_EXTERNAL_L1L2" "-DUSE_EXTERNAL_FIRSTORDER" "-DMPAS_LI_BUILD_INTERFACE")
    endif()

    # Albany can only solve FO at present
    if (ALBANY)
      list(APPEND CPPDEFS "-DUSE_EXTERNAL_FIRSTORDER" "-DMPAS_LI_BUILD_INTERFACE")
    endif()

    if (LIFEV AND ALBANY)
      message(FATAL "Compiling with both LifeV and Albany is not allowed at this time.")
    endif()

    # PHG currently requires LifeV
    if (PHG AND NOT LIFEV)
      message(FATAL "Compiling with PHG requires LifeV at this time.")
    endif()

    # PHG can only Stokes at present
    if (PHG)
      list(APPEND CPPDEFS "-DUSE_EXTERNAL_STOKES" "-DMPAS_LI_BUILD_INTERFACE")
    endif()
  endif()

  add_library(${COMPONENT})
  target_compile_definitions(${COMPONENT} PRIVATE ${CPPDEFS})
  target_include_directories(${COMPONENT} PRIVATE ${INCLUDES})

  # Gather sources
  set(CORE_BLDDIR ${CMAKE_BINARY_DIR}/core_${CORE})
  if (NOT EXISTS ${CORE_BLDDIR})
    file(MAKE_DIRECTORY ${CORE_BLDDIR})
  endif()

  set(CORE_INPUT_DIR ${CORE_BLDDIR}/default_inputs)
  if (NOT EXISTS ${CORE_INPUT_DIR})
    file(MAKE_DIRECTORY ${CORE_INPUT_DIR})
  endif()

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
    COMMAND ${CMAKE_BINARY_DIR}/mpas-source/src/tools/parse < ${CORE_BLDDIR}/Registry_processed.xml
    DEPENDS parse ${CORE_BLDDIR}/Registry_processed.xml
    WORKING_DIRECTORY ${INC_DIR}
  )

  include(${CMAKE_CURRENT_SOURCE_DIR}/core_${CORE}/${CORE}.cmake)

  # Disable qsmp for some files
  if (FFLAGS MATCHES ".*-qsmp.*")
    foreach(DISABLE_QSMP_FILE IN LISTS DISABLE_QSMP)
      get_filename_component(SOURCE_EXT ${DISABLE_QSMP_FILE} EXT)
      string(REPLACE "${SOURCE_EXT}" ".f90" SOURCE_F90 ${DISABLE_QSMP_FILE})
      set_property(SOURCE ${CMAKE_BINARY_DIR}/${SOURCE_F90} APPEND_STRING PROPERTY COMPILE_FLAGS " -nosmp")
    endforeach()
  endif()

  genf90_targets("${RAW_SOURCES}" "${INCLUDES}" "${CPPDEFS}" "${NO_PREPROCESS}" "${INC_DIR}")
  target_sources(${COMPONENT} PRIVATE ${SOURCES} $<TARGET_OBJECTS:common>)

endfunction(build_core)
