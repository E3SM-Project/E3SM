function(build_model MODEL_ARG MODELCONF_ARG FILEPATH_ARG)
  # Load dependency search path.
  file(STRINGS ${FILEPATH_ARG} FILEPATH_DIRS)
  set(CPP_DIRS ${FILEPATH_DIRS})
  list(APPEND CPP_DIRS ".")

  # Load cppdefs
  set(CCSM_CPPDEFS_FILE "${MODELCONF_ARG}/CCSM_cppdefs")
  if (EXISTS ${CCSM_CPPDEFS_FILE})
    file(READ ${CCSM_CPPDEFS_FILE} CCSM_CPPDEFS)
    string(STRIP "${CCSM_CPPDEFS}" CCSM_CPPDEFS)
    set(CPPDEFS "${CPPDEFS} ${CCSM_CPPDEFS}")
  endif()

  if (MODEL_ARG STREQUAL "cpl")
    set(INCLDIR "${INCLDIR} -I${EXEROOT}/atm/obj -I${EXEROOT}/ice/obj -I${EXEROOT}/ocn/obj -I${EXEROOT}/glc/obj -I${EXEROOT}/rof/obj -I${EXEROOT}/wav/obj -I${EXEROOT}/esp/obj -I${EXEROOT}/iac/obj")
  endif()

  if (MODEL_ARG STREQUAL "cam")
    # These RRTMG files take an extraordinarily long time to compile with optimization.
    # Until mods are made to read the data from files, just remove optimization from
    # their compilation.
    set(NOOPT_FILES "cam/src/physics/rrtmg/ext/rrtmg_lw/rrtmg_lw_k_g.f90;cam/src/physics/rrtmg/ext/rrtmg_sw/rrtmg_sw_k_g.f90")

    if (COSP_LIBDIR)
      include(${CMAKE_CURRENT_SOURCE_DIR}/cam/src/physics/cosp/Cosp.cmake)
    endif()

  endif()

  #-------------------------------------------------------------------------------
  # Build & include dependency files
  #-------------------------------------------------------------------------------

  # Get src files
  # JGF: Cmake does not have a VPATH concept, so we need relative/absolute paths to source files
  # Note: Using absolute paths seems to wreck CMake's ability to do a dep analysis on fortran
  # sources and compile them in the right order.
  #
  # One additional subtley is that, when mkSrcfiles found multiple files with the same basename,
  # only the one found first gets compiled.

  set(SRCROOT_REL "${CIMEROOT}/..")
  set(BASENAME_SET)
  file(TO_CMAKE_PATH ${SRCROOT_REL} SRCROOT_ABS)
  foreach(DIRSEARCH ${FILEPATH_DIRS})
    file(GLOB MATCHES RELATIVE "${SRCROOT_ABS}/components" "${DIRSEARCH}/*.[Ffc]" "${DIRSEARCH}/*.[Ff]90" "${DIRSEARCH}/*.cpp" "${DIRSEARCH}/*.F90.in")
    if (MATCHES)
      foreach (MATCH IN LISTS MATCHES)
        get_filename_component(BASENAME ${MATCH} NAME)
        list(FIND BASENAME_SET ${BASENAME} BASENAME_WAS_FOUND)
        if (BASENAME_WAS_FOUND EQUAL -1)
          list(APPEND SOURCES ${MATCH})
          list(APPEND BASENAME_SET ${BASENAME})
        else()
          message(WARNING "Skipping repeated base filename ${BASENAME} for ${MATCH}")
        endif()
      endforeach()
    endif()
  endforeach()

  foreach(SOURCE_FILE IN LISTS SOURCES)
i    get_filename_component(SOURCE_EXT ${SOURCE_FILE} EXT)
    if (SOURCE_EXT STREQUAL ".F90.in")
      string(REPLACE ".in" "" SOURCE_NO_IN ${SOURCE_FILE})
      list(APPEND GEN_F90_SOURCES ${SOURCE_NO_IN})
      list(APPEND SOURCES ${SOURCE_NO_IN})
      list(REMOVE_ITEM SOURCES ${SOURCE_FILE})
    endif()
  endforeach()

  foreach(ITEM IN LISTS CPP_DIRS)
    if (EXISTS ${ITEM})
      set(INCLDIR "${INCLDIR} -I${ITEM}")
    endif()
  endforeach()

  #-------------------------------------------------------------------------------
  # create list of component libraries - hard-wired for current ccsm components
  #-------------------------------------------------------------------------------

  if (CIME_MODEL STREQUAL "cesm")
    if (COMP_LND STREQUAL "clm")
      set(USE_SHARED_CLM TRUE)
    else()
      set(USE_SHARED_CLM FALSE)
    endif()
  else()
    set(USE_SHARED_CLM FALSE)
  endif()

  if (NOT USE_SHARED_CLM)
    set(LNDOBJDIR "${EXEROOT}/lnd/obj")
    set(LNDLIBDIR "${LIBROOT}")
    if (COMP_LND STREQUAL "clm")
      set(LNDLIB "libclm.a")
    else()
      set(LNDLIB "liblnd.a")
    endif()
    set(INCLDIR "${INCLDIR} -I${LNDOBJDIR}")
  else()
    set(LNDLIB "libclm.a")
    set(LNDOBJDIR "${SHAREDLIBROOT}/${SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/clm/obj")
    set(LNDLIBDIR "${EXEROOT}/${SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/lib")
    set(INCLDIR "${INCLDIR} -I${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/include")
    if (MODEL_ARG STREQUAL "clm")
      set(INCLUDE_DIR "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/include")
    endif()
  endif()

  if (NOT ULIBDEP)
    if (LIBROOT)
      set(ULIBDEP "${LIBROOT}/libatm.a")
      set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libice.a")
      set(ULIBDEP "${ULIBDEP} ${LNDLIBDIR}/${LNDLIB}")
      set(INCLDIR "${INCLDIR} -I${LNDOBJDIR}")
      set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libocn.a")
      set(ULIBDEP "${ULIBDEP} ${LIBROOT}/librof.a")
      set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libglc.a")
      set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libwav.a")
      if (NOT COMP_INTERFACE STREQUAL "nuopc")
        set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libiac.a")
      endif()
      set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libesp.a")
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

  #------------------------------------------------------------------------------
  # Set key cmake vars
  #------------------------------------------------------------------------------

  set(CMAKE_Fortran_FLAGS "${FFLAGS} ${INCLDIR}")
  set(CMAKE_C_FLAGS "${CFLAGS} ${INCLDIR}")
  set(CMAKE_CXX_FLAGS "${CXXFLAGS} ${INCLDIR}")
  set(CMAKE_VERBOSE_MAKEFILE TRUE)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../../lib)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../../lib)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../..)

  set(CMAKE_C_COMPILER ${CC})
  set(CMAKE_CXX_COMPILER ${CXX})
  set(CMAKE_Fortran_COMPILER ${FC})
  set(CMAKE_EXE_LINKER_FLAGS "${LDFLAGS}")

  # do the necessary genf90s
  add_custom_target(genf90
    DEPENDS ${CIMEROOT}/src/externals/genf90/genf90.pl)
  foreach (SRC_FILE IN LISTS GEN_F90_SOURCES)
    get_filename_component(BASENAME ${SRC_FILE} NAME)
    add_custom_command (
      OUTPUT ${CMAKE_BINARY_DIR}/${BASENAME}
      COMMAND ${CIMEROOT}/src/externals/genf90/genf90.pl
      ${CMAKE_CURRENT_SOURCE_DIR}/${SRC_FILE}.in > ${CMAKE_BINARY_DIR}/${BASENAME}
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${SRC_FILE}.in genf90)
    list(REMOVE_ITEM SOURCES ${SRC_FILE})
    list(APPEND SOURCES ${CMAKE_BINARY_DIR}/${BASENAME})
  endforeach ()

  # Flags are slightly different for different fortran extensions
  foreach (SOURCE_FILE IN LISTS SOURCES)
    # Cosp manages its own flags
    if (NOT SOURCE_FILE IN_LIST COSP_SOURCES)
      get_filename_component(SOURCE_EXT ${SOURCE_FILE} EXT)
      if (SOURCE_EXT STREQUAL ".F" OR SOURCE_EXT STREQUAL ".f")
        set_property(SOURCE ${SOURCE_FILE} APPEND_STRING PROPERTY COMPILE_FLAGS " ${FIXEDFLAGS}")
      elseif(SOURCE_EXT STREQUAL ".f90")
        set_property(SOURCE ${SOURCE_FILE} APPEND_STRING PROPERTY COMPILE_FLAGS " ${FREEFLAGS}")
      elseif(SOURCE_EXT STREQUAL ".F90")
        set_property(SOURCE ${SOURCE_FILE} APPEND_STRING PROPERTY COMPILE_FLAGS " ${FREEFLAGS} ${CONTIGUOUS_FLAG}")
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
    set_property(SOURCE ${SOURCE_FILE} APPEND_STRING PROPERTY COMPILE_FLAGS " ${FFLAGS_NOOPT}")
  endforeach()

  #-------------------------------------------------------------------------------
  # build rules:
  #-------------------------------------------------------------------------------

  if (MPILIB STREQUAL "mpi-serial")
    set(MPISERIAL "${INSTALL_SHAREDPATH}/lib/lib-mpi-serial.a")
    set(MLIBS "${MLIBS} ${MPISERIAL}")
  endif()

  add_library(${MODEL_ARG})
  target_sources(${MODEL_ARG} PRIVATE ${SOURCES})

  # Subtle: In order for fortran dependency scanning to work, our CPPFPP/DEFS must be registered
  # as COMPILE_DEFINITIONS, not simple added via CMAKE_Fortran_Flags. Also, CPPDEFS *must*
  # be provided as a list, not a whitespace-separated string; otherwise, things get wonky.
  separate_arguments(CPPDEFS_LIST UNIX_COMMAND "${CPPDEFS}")
  target_compile_definitions(${MODEL_ARG} PRIVATE ${CPPDEFS_LIST})
  add_dependencies(${MODEL_ARG} genf90)

endfunction(build_model)
