function(build_model COMP_CLASS COMP_NAME)

  # We support component-specific configuration of flags, etc, so this setup
  # need to be done here.
  include(${CMAKE_SOURCE_DIR}/cmake/common_setup.cmake)

  set(MODELCONF_DIR "${BUILDCONF}/${COMP_NAME}conf")

  # Load dependency search path.
  file(STRINGS ${MODELCONF_DIR}/Filepath FILEPATH_DIRS)
  set(CPP_DIRS ${FILEPATH_DIRS})
  list(APPEND CPP_DIRS ".")

  # Load cppdefs
  # Tries CIME_cppdefs first, and if the file does not exist then tries CCSM_cppdefs
  set(CIME_CPPDEFS_FILE "${MODELCONF_DIR}/CIME_cppdefs")
  set(CCSM_CPPDEFS_FILE "${MODELCONF_DIR}/CCSM_cppdefs")
  if (EXISTS ${CIME_CPPDEFS_FILE})
    file(READ ${CIME_CPPDEFS_FILE} CIME_CPPDEFS)
    string(STRIP "${CIME_CPPDEFS}" CIME_CPPDEFS)
    set(CPPDEFS "${CPPDEFS} ${CIME_CPPDEFS}")
  elseif (EXISTS ${CCSM_CPPDEFS_FILE})
    file(READ ${CCSM_CPPDEFS_FILE} CCSM_CPPDEFS)
    string(STRIP "${CCSM_CPPDEFS}" CCSM_CPPDEFS)
    set(CPPDEFS "${CPPDEFS} ${CCSM_CPPDEFS}")
  endif()

  if (COMP_NAME STREQUAL "cpl")
    list(APPEND INCLDIR "${EXEROOT}/cmake-bld/mpas-source/src")
    foreach(ITEM IN LISTS COMP_CLASSES)
      list(APPEND INCLDIR "${EXEROOT}/cmake-bld/cmake/${ITEM}")
    endforeach()
  endif()

  # Source files don't live in components/cmake/$COMP_CLASS. This path won't work for
  # generated files.
  set(SOURCE_PATH "../..")
  set(CIMESRC_PATH "../cime/src")

  #-------------------------------------------------------------------------------
  # Build & include dependency files
  #-------------------------------------------------------------------------------

  gather_sources("${FILEPATH_DIRS}" "${CIMEROOT}")
  set(SOURCES ${SOURCES_RESULT})
  set(GEN_F90_SOURCES ${GEN_F90_SOURCES_RESULT})

  foreach(ITEM IN LISTS CPP_DIRS)
    if (EXISTS ${ITEM})
      list(APPEND INCLDIR "${ITEM}")
    endif()
  endforeach()

  #-------------------------------------------------------------------------------
  # Cam needs some special handling for cosp and turning off opts for some files
  #-------------------------------------------------------------------------------

  if (COMP_NAME STREQUAL "eam")
    # These RRTMG files take an extraordinarily long time to compile with optimization.
    # Until mods are made to read the data from files, just remove optimization from
    # their compilation.
    set(NOOPT_FILES "eam/src/physics/rrtmg/ext/rrtmg_lw/rrtmg_lw_k_g.f90;eam/src/physics/rrtmg/ext/rrtmg_sw/rrtmg_sw_k_g.f90")

    if (USE_COSP)
      include(${PROJECT_SOURCE_DIR}/eam/src/physics/cosp2/Cosp.cmake)
    endif()

    # If YAKL is needed, then set YAKL CMake vars
    if (USE_YAKL)
      # ARCH can be CUDA, HIP, or unset
      if (USE_CUDA)
        set(ARCH "CUDA")
        # CUDA_FLAGS is set through Macros.cmake / config_compilers.xml
        # We can't have duplicate flags with nvcc, so we only specify CPPDEFS,
        # and the rest is up to CUDAFLAGS
        set(YAKL_CXX_FLAGS "${CPPDEFS}")
      else()
        # For normal C++ compilers duplicate flags are fine, the last ones win typically
        set(YAKL_CXX_FLAGS "${CPPDEFS} ${CXXFLAGS}")
        set(ARCH "")
      endif()
      message(STATUS "Building YAKL")
      # Build YAKL as a static library
      # YAKL_HOME is YAKL's source directlry
      set(YAKL_HOME ${CMAKE_CURRENT_SOURCE_DIR}/../../../externals/YAKL)
      # YAKL_BIN is where we're placing the YAKL library
      set(YAKL_BIN  ${CMAKE_CURRENT_BINARY_DIR}/yakl)
      # YAKL_CUB_HOME is where Nvidia's cub repo lives (submodule)
      if (USE_CUDA)
        set(YAKL_CUB_HOME ${CMAKE_CURRENT_SOURCE_DIR}/../../../externals/cub)
      endif()
      # Build the YAKL static library
      add_subdirectory(${YAKL_HOME} ${YAKL_BIN})
      # Add both YAKL source and YAKL binary directories to the include paths
      include_directories(${YAKL_HOME} ${YAKL_BIN})
    endif()

    # if samxx is needed, build samxx as a static library
    if (USE_SAMXX)
      message(STATUS "Building SAMXX")
      # SAMXX_HOME is where the samxx source code lives
      set(SAMXX_HOME ${CMAKE_CURRENT_SOURCE_DIR}/../../eam/src/physics/crm/samxx)
      # SAMXX_BIN is where the samxx library will live
      set(SAMXX_BIN  ${CMAKE_CURRENT_BINARY_DIR}/samxx)
      # Build the static samxx library
      add_subdirectory(${SAMXX_HOME} ${SAMXX_BIN})
      # Add samxx F90 files to the main E3SM build
      set(SOURCES ${SOURCES} cmake/atm/../../eam/src/physics/crm/samxx/cpp_interface_mod.F90
                             cmake/atm/../../eam/src/physics/crm/samxx/params.F90
                             cmake/atm/../../eam/src/physics/crm/samxx/crm_ecpp_output_module.F90 )
    endif()

  endif()

  #-------------------------------------------------------------------------------
  # create list of component libraries - hard-wired for current e3sm components
  #-------------------------------------------------------------------------------

  if (NOT USE_SHARED_CLM)
    set(LNDOBJDIR "${EXEROOT}/lnd/obj")
    set(LNDLIBDIR "${LIBROOT}")
    if (COMP_LND STREQUAL "elm")
      set(LNDLIB "libelm.a")
    else()
      set(LNDLIB "liblnd.a")
    endif()
    list(APPEND INCLDIR "${LNDOBJDIR}")
  else()
    set(LNDLIB "libclm.a")
    set(LNDOBJDIR "${SHAREDLIBROOT}/${SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/clm/obj")
    set(LNDLIBDIR "${EXEROOT}/${SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/lib")
    list(APPEND INCLDIR "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/include")
    if (COMP_NAME STREQUAL "clm")
      set(INCLUDE_DIR "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/include")
    endif()
  endif()

  if (NOT ULIBDEP)
    if (LIBROOT)
      list(APPEND INCLDIR "${LNDOBJDIR}")
    endif()
  endif()

  if (COMP_GLC STREQUAL "cism")
    set(ULIBDEP "${ULIBDEP} ${CISM_LIBDIR}/libglimmercismfortran.a")
    if (CISM_USE_TRILINOS)
      set(ULIBDEP "${ULIBDEP} ${CISM_LIBDIR}/libglimmercismcpp.a")
    endif()
  endif()
  if (OCN_SUBMODEL STREQUAL "moby")
    set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libmoby.a")
  endif()

  set(CSMSHARELIB "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/lib/libcsm_share.a")
  set(ULIBDEP "${ULIBDEP} ${CSMSHARELIB}")

  # do the necessary genf90s
  foreach (SRC_FILE IN LISTS GEN_F90_SOURCES)
    get_filename_component(BASENAME ${SRC_FILE} NAME)
    add_custom_command (
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${BASENAME}
      COMMAND ${CIMEROOT}/src/externals/genf90/genf90.pl
      ${PROJECT_SOURCE_DIR}/${SRC_FILE}.in > ${CMAKE_CURRENT_BINARY_DIR}/${BASENAME}
      DEPENDS ${PROJECT_SOURCE_DIR}/${SRC_FILE}.in genf90)
    list(REMOVE_ITEM SOURCES ${SRC_FILE})
    list(APPEND SOURCES ${CMAKE_CURRENT_BINARY_DIR}/${BASENAME})
  endforeach ()

  # Set compiler flags for source files
  foreach (SOURCE_FILE IN LISTS SOURCES)
    get_filename_component(SOURCE_EXT ${SOURCE_FILE} EXT)

    # Set flags based on file extension. File extensions may change if the globs used by
    # gather_sources changes.
    if (SOURCE_EXT STREQUAL ".c")
      e3sm_add_flags("${SOURCE_FILE}" "${CFLAGS}")
    elseif (SOURCE_EXT STREQUAL ".cpp")
      e3sm_add_flags("${SOURCE_FILE}" "${CXXFLAGS}")
    else()
      # This is a fortran source
      e3sm_add_flags("${SOURCE_FILE}" "${FFLAGS}")

      # Cosp manages its own flags
      if (NOT SOURCE_FILE IN_LIST COSP_SOURCES)
        # Flags are slightly different for different fortran extensions
        if (SOURCE_EXT STREQUAL ".F" OR SOURCE_EXT STREQUAL ".f")
          e3sm_add_flags("${SOURCE_FILE}" "${FIXEDFLAGS}")
        elseif (SOURCE_EXT STREQUAL ".f90")
          e3sm_add_flags("${SOURCE_FILE}" "${FREEFLAGS}")
        elseif (SOURCE_EXT STREQUAL ".F90")
          e3sm_add_flags("${SOURCE_FILE}" "${FREEFLAGS} ${CONTIGUOUS_FLAG}")
        endif()
      endif()
    endif()
  endforeach()

  # Load machine/compiler specific settings
  set(COMPILER_SPECIFIC_DEPENDS ${CASEROOT}/Depends.${COMPILER}.cmake)
  set(MACHINE_SPECIFIC_DEPENDS ${CASEROOT}/Depends.${MACH}.cmake)
  set(PLATFORM_SPECIFIC_DEPENDS ${CASEROOT}/Depends.${MACH}.${COMPILER}.cmake)
  set(TRY_TO_LOAD ${COMPILER_SPECIFIC_DEPENDS} ${MACHINE_SPECIFIC_DEPENDS} ${PLATFORM_SPECIFIC_DEPENDS})
  foreach(ITEM IN LISTS TRY_TO_LOAD)
    if (EXISTS ${ITEM})
      include(${ITEM})
    endif()
  endforeach()

  # Disable optimizations on some files that would take too long to compile, expect these to all be fortran files
  foreach (SOURCE_FILE IN LISTS NOOPT_FILES)
    e3sm_deoptimize_file("${SOURCE_FILE}" "${FFLAGS_NOOPT}")
  endforeach()

  #-------------------------------------------------------------------------------
  # build rules:
  #-------------------------------------------------------------------------------

  if (MPILIB STREQUAL "mpi-serial")
    set(MPISERIAL "${INSTALL_SHAREDPATH}/lib/libmpi-serial.a")
    set(MLIBS "${MLIBS} ${MPISERIAL}")
  endif()

  foreach(ITEM IN LISTS SOURCES)
    if (ITEM MATCHES "${CMAKE_BINARY_DIR}/.*") # is generated
      list(APPEND REAL_SOURCES ${ITEM})
    else()
      list(APPEND REAL_SOURCES "${SOURCE_PATH}/${ITEM}")
    endif()
  endforeach()

  if (COMP_NAME STREQUAL "cpl")
    set(TARGET_NAME "${CIME_MODEL}.exe")
    add_executable(${TARGET_NAME})
    target_sources(${TARGET_NAME} PRIVATE ${REAL_SOURCES})
    set(ALL_LIBS "${ULIBDEP} ${MCTLIBS} ${PIOLIB} ${GPTLLIB} ${SLIBS} ${MLIBS}")
    separate_arguments(ALL_LIBS_LIST UNIX_COMMAND "${ALL_LIBS}")
    foreach(ITEM IN LISTS COMP_CLASSES)
      if (NOT ITEM STREQUAL "cpl")
        target_link_libraries(${TARGET_NAME} ${ITEM})
      endif()
    endforeach()
    foreach(ITEM IN LISTS ALL_LIBS_LIST)
      target_link_libraries(${TARGET_NAME} ${ITEM})
    endforeach()
    set_target_properties(${TARGET_NAME} PROPERTIES LINKER_LANGUAGE ${LD})
  else()
    set(TARGET_NAME ${COMP_CLASS})
    add_library(${TARGET_NAME})
    target_sources(${TARGET_NAME} PRIVATE ${REAL_SOURCES})
    if (COMP_NAME STREQUAL "eam")
      if (USE_YAKL)
        target_link_libraries(${TARGET_NAME} PRIVATE yakl)
      endif()
      if (USE_SAMXX)
        target_link_libraries(${TARGET_NAME} PRIVATE samxx)
      endif()
    endif()
    if (USE_KOKKOS)
      target_link_libraries (${TARGET_NAME} PRIVATE Kokkos::kokkos)
    endif ()
  endif()

  # Subtle: In order for fortran dependency scanning to work, our CPPFPP/DEFS must be registered
  # as COMPILE_DEFINITIONS, not simple added via CMAKE_Fortran_Flags. Also, CPPDEFS *must*
  # be provided as a list, not a whitespace-separated string; otherwise, things get wonky.
  separate_arguments(CPPDEFS_LIST UNIX_COMMAND "${CPPDEFS}")
  target_compile_definitions(${TARGET_NAME} PRIVATE ${CPPDEFS_LIST})
  add_dependencies(${TARGET_NAME} genf90)

  # Set flags for target
  target_include_directories(${TARGET_NAME} PRIVATE ${INCLDIR})

endfunction(build_model)
