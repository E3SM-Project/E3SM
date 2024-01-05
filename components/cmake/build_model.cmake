# This function must be a macro so that it does not get its own scope.
# This way, changes to CMAKE_* vars affect the directory which is what
# we want.
macro(build_model COMP_CLASS COMP_NAME)

  # We want real variables, not macro expansions
  set(COMP_CLASS ${COMP_CLASS})
  set(COMP_NAME ${COMP_NAME})

  # We support component-specific configuration of flags, etc, so this setup
  # need to be done here.
  include(${SRCROOT}/components/cmake/common_setup.cmake)

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
    list(APPEND INCLDIR "${EXEROOT}/cmake-bld/mpas-framework/src")
    foreach(ITEM IN LISTS COMP_CLASSES)
      list(APPEND INCLDIR "${EXEROOT}/cmake-bld/cmake/${ITEM}")
    endforeach()
  endif()

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
      # YAKL_ARCH can be CUDA, HIP, SYCL, OPENMP45, or empty
      # USE_CUDA or USE_HIP are set through Macros.cmake
      if (USE_CUDA)
        set(YAKL_ARCH "CUDA")
        # CUDA_FLAGS is set through Macros.cmake
        # For instance: cime_config/machines/cmake_macros/gnugpu_summit.cmake
        set(YAKL_CUDA_FLAGS "${CPPDEFS} ${CUDA_FLAGS}")
      elseif (USE_HIP)
        set(YAKL_ARCH "HIP")
        # HIP_FLAGS are set through Macros.cmake
        # For instance: cime_config/machines/cmake_macros/crayclanggpu_frontier.cmake
        set(YAKL_HIP_FLAGS "${CPPDEFS} ${HIP_FLAGS}")
      elseif (USE_SYCL)
        set(YAKL_ARCH "SYCL")
        # SYCL_FLAGS is set through Macros.cmake
        # For instance: cime_config/machines/cmake_macros/oneapi-ifxgpu_sunspot.cmake
        set(YAKL_SYCL_FLAGS "${CPPDEFS} ${SYCL_FLAGS}")
      else()
        # For CPU C++ compilers duplicate flags are fine, the last ones win typically
        set(YAKL_CXX_FLAGS "${CPPDEFS} ${CMAKE_CXX_FLAGS}")
        set(YAKL_ARCH "")
      endif()
      message(STATUS "Building YAKL")
      # Build YAKL as a static library
      # YAKL_HOME is YAKL's source directlry
      set(YAKL_HOME ${CMAKE_CURRENT_SOURCE_DIR}/../../../externals/YAKL)
      # YAKL_BIN is where we're placing the YAKL library
      set(YAKL_BIN  ${CMAKE_CURRENT_BINARY_DIR}/yakl)
      # Build the YAKL static library
      add_subdirectory(${YAKL_HOME} ${YAKL_BIN})
      # Add the YAKL bin directory, mainly due to a Fortran module if it's needed
      include_directories(${YAKL_BIN})
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
                             cmake/atm/../../eam/src/physics/crm/crm_ecpp_output_module.F90 )
    endif()

    if (USE_PAM)
      message(STATUS "Building PAM")
      # PAM_HOME is where the samxx source code lives
      set(PAM_HOME ${CMAKE_CURRENT_SOURCE_DIR}/../../eam/src/physics/crm/pam)
      # PAM_BIN is where the samxx library will live
      set(PAM_BIN  ${CMAKE_CURRENT_BINARY_DIR}/pam)
      # Build the static samxx library
      add_subdirectory(${PAM_HOME} ${PAM_BIN})
      # Add samxx F90 files to the main E3SM build
      set(SOURCES ${SOURCES} cmake/atm/../../eam/src/physics/crm/pam/params.F90
                             cmake/atm/../../eam/src/physics/crm/crm_ecpp_output_module.F90
                             cmake/atm/../../eam/src/physics/crm/pam/pam_driver.F90)
      # Pam interface need to include modules from pam
      include_directories(${PAM_BIN}/external/pam_core)
    endif()

    # Add rrtmgp++ source code if asked for
    if (USE_RRTMGPXX)
      message(STATUS "Building RRTMGPXX")
      # Build the static rrtmgpxx library
      set(RRTMGPXX_BIN ${CMAKE_CURRENT_BINARY_DIR}/rrtmgp)
      add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../eam/src/physics/rrtmgp/external/cpp ${RRTMGPXX_BIN})
      # Build the interface code
      set(RRTMGPXX_INTERFACE_BIN ${CMAKE_CURRENT_BINARY_DIR}/rrtmgp_interface)
      add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../eam/src/physics/rrtmgp/cpp ${RRTMGPXX_INTERFACE_BIN})
      # Interface code needs some additional headers
      target_include_directories(rrtmgp_interface PRIVATE
          ${CMAKE_CURRENT_SOURCE_DIR}/../../eam/src/physics/rrtmgp/external/cpp/extensions/fluxes_byband
          ${CMAKE_CURRENT_SOURCE_DIR}/../../eam/src/physics/rrtmgp/external/cpp/extensions/cloud_optics
          ${CMAKE_CURRENT_SOURCE_DIR}/../../eam/src/physics/rrtmgp/cpp
      )
      # The interface code needs to know about the NETCDF includes defined
      # above. The easiest way I know of to do this is to pass all of the
      # accumulated includes to the target.
      # TODO: this can go away if the above NETCDF section is refactored to
      # use find_library instead of appending to INCLDIR.
      target_include_directories(rrtmgp_interface PRIVATE ${INCLDIR})
      # Add the source files for the interface code to the main E3SM build
      set(RRTMGPXX_F90 cmake/atm/../../eam/src/physics/rrtmgp/cpp/rrtmgp_interface.F90)
      set(SOURCES ${SOURCES} ${RRTMGPXX_F90})
    endif()
  endif()

  #-------------------------------------------------------------------------------
  # create list of component libraries - hard-wired for current e3sm components
  #-------------------------------------------------------------------------------

  # do the necessary genf90s
  foreach (SRC_FILE IN LISTS GEN_F90_SOURCES)
    get_filename_component(BASENAME ${SRC_FILE} NAME)
    add_custom_command (
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${BASENAME}
      COMMAND ${CIMEROOT}/CIME/non_py/externals/genf90/genf90.pl
      ${PROJECT_SOURCE_DIR}/${SRC_FILE}.in > ${CMAKE_CURRENT_BINARY_DIR}/${BASENAME}
      DEPENDS ${PROJECT_SOURCE_DIR}/${SRC_FILE}.in genf90)
    list(REMOVE_ITEM SOURCES ${SRC_FILE})
    list(APPEND SOURCES ${CMAKE_CURRENT_BINARY_DIR}/${BASENAME})
  endforeach ()

  # Set special fixed/free compiler flags for non-cosp fortran files. This logic
  # can maybe be replaced/reworked once we have Fortran_FORMAT property
  # working
  foreach (SOURCE_FILE IN LISTS SOURCES)
    get_filename_component(SOURCE_EXT ${SOURCE_FILE} EXT)

    # Cosp manages its own flags
    if (NOT SOURCE_FILE IN_LIST COSP_SOURCES)
      # Flags are slightly different for different fortran extensions
      if (SOURCE_EXT STREQUAL ".F" OR SOURCE_EXT STREQUAL ".f")
        e3sm_set_source_property(${SOURCE_FILE} Fortran_FORMAT FIXED FALSE)
      elseif (SOURCE_EXT STREQUAL ".f90")
        e3sm_set_source_property(${SOURCE_FILE} Fortran_FORMAT FREE FALSE)
      elseif (SOURCE_EXT STREQUAL ".F90")
        e3sm_set_source_property(${SOURCE_FILE} Fortran_FORMAT FREE FALSE)
        e3sm_add_flags(${SOURCE_FILE} "${CONTIGUOUS_FLAG}")
      endif()
    endif()
  endforeach()

  # Load machine/compiler specific settings
  if (COMP_NAME STREQUAL "csm_share")
    # csm_share uses special Depends files, not customizable per case
    set(DEPENDS_LOC "${SRCROOT}/share")
  else()
    set(DEPENDS_LOC "${CASEROOT}")
  endif()

  set(COMPILER_SPECIFIC_DEPENDS ${DEPENDS_LOC}/Depends.${COMPILER}.cmake)
  set(MACHINE_SPECIFIC_DEPENDS ${DEPENDS_LOC}/Depends.${MACH}.cmake)
  set(PLATFORM_SPECIFIC_DEPENDS ${DEPENDS_LOC}/Depends.${MACH}.${COMPILER}.cmake)
  set(TRY_TO_LOAD ${COMPILER_SPECIFIC_DEPENDS} ${MACHINE_SPECIFIC_DEPENDS} ${PLATFORM_SPECIFIC_DEPENDS})
  foreach(ITEM IN LISTS TRY_TO_LOAD)
    if (EXISTS ${ITEM})
      include(${ITEM})
    endif()
  endforeach()

  # Disable optimizations on some files that would take too long to compile, expect these to all be fortran files
  foreach (SOURCE_FILE IN LISTS NOOPT_FILES)
    e3sm_deoptimize_file("${SOURCE_FILE}")
  endforeach()

  #-------------------------------------------------------------------------------
  # build rules:
  #-------------------------------------------------------------------------------

  foreach(ITEM IN LISTS SOURCES)
    if (ITEM MATCHES "${CMAKE_BINARY_DIR}/.*") # is generated
      list(APPEND REAL_SOURCES ${ITEM})
    else()
      list(APPEND REAL_SOURCES "${PROJECT_SOURCE_DIR}/${ITEM}")
    endif()
  endforeach()

  if (COMP_NAME STREQUAL "cpl")
    set(TARGET_NAME "${CIME_MODEL}.exe")
    add_executable(${TARGET_NAME})
    target_sources(${TARGET_NAME} PRIVATE ${REAL_SOURCES})

    separate_arguments(ALL_LIBS_LIST UNIX_COMMAND "${SLIBS}")

    foreach(ITEM IN LISTS COMP_CLASSES)
      if (NOT ITEM STREQUAL "cpl")
        target_link_libraries(${TARGET_NAME} ${ITEM})
      endif()
    endforeach()

    foreach(ITEM IN LISTS ALL_LIBS_LIST)
      target_link_libraries(${TARGET_NAME} ${ITEM})
    endforeach()

    # Make sure we link blas/lapack
    target_link_libraries(${TARGET_NAME} BLAS::BLAS LAPACK::LAPACK)

    if (E3SM_LINK_WITH_FORTRAN)
      set_target_properties(${TARGET_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

      # A bit hacky, some platforms need help with the fortran linker
      if (COMPILER STREQUAL "intel")
        string(APPEND CMAKE_EXE_LINKER_FLAGS " -cxxlib")
      endif()

    else()
      set_target_properties(${TARGET_NAME} PROPERTIES LINKER_LANGUAGE CXX)
    endif()

  else()
    set(TARGET_NAME ${COMP_CLASS})
    add_library(${TARGET_NAME})
    target_sources(${TARGET_NAME} PRIVATE ${REAL_SOURCES})
    if (COMP_NAME STREQUAL "csm_share")
      find_package(NETCDF REQUIRED)
      target_link_libraries(${TARGET_NAME} PRIVATE netcdf)
    else()
      target_link_libraries(${TARGET_NAME} PRIVATE csm_share)
      if (COMP_NAME STREQUAL "eam")
        if (USE_YAKL)
          target_link_libraries(${TARGET_NAME} PRIVATE yakl)
        endif()
        if (USE_SAMXX)
          target_link_libraries(${TARGET_NAME} PRIVATE samxx)
        endif()
        if (USE_PAM)
          target_link_libraries(${TARGET_NAME} PRIVATE pam_driver)
        endif()
        if (USE_RRTMGPXX)
          target_link_libraries(${TARGET_NAME} PRIVATE rrtmgp rrtmgp_interface)
        endif()
      endif()
      if (COMP_NAME STREQUAL "elm")
        if (USE_PETSC)
          target_link_libraries(${TARGET_NAME} PRIVATE "${PETSC_LIBRARIES}")
          target_include_directories(${TARGET_NAME} PRIVATE "${PETSC_INCLUDES}")
        endif()
      endif()
      if (USE_KOKKOS)
        target_link_libraries (${TARGET_NAME} PRIVATE Kokkos::kokkos)
      endif ()
    endif()
  endif()

  # Subtle: In order for fortran dependency scanning to work, our CPPFPP/DEFS must be registered
  # as COMPILE_DEFINITIONS, not simple added via CMAKE_Fortran_Flags. Also, CPPDEFS *must*
  # be provided as a list, not a whitespace-separated string; otherwise, things get wonky.
  separate_arguments(CPPDEFS_LIST UNIX_COMMAND "${CPPDEFS}")
  target_compile_definitions(${TARGET_NAME} PRIVATE ${CPPDEFS_LIST})
  add_dependencies(${TARGET_NAME} genf90)

  # Set flags for target
  target_include_directories(${TARGET_NAME} PRIVATE ${INCLDIR})

endmacro(build_model)
