# build Omega

# create library flag target for LibOmega
add_library(OmegaLibFlags INTERFACE)

target_include_directories(
    OmegaLibFlags
    INTERFACE
    ${OMEGA_SOURCE_DIR}/src/analysis
    ${OMEGA_SOURCE_DIR}/src/base
    ${OMEGA_SOURCE_DIR}/src/infra
    ${OMEGA_SOURCE_DIR}/src/ocn
    ${OMEGA_SOURCE_DIR}/src/timeStepping
    ${E3SM_ROOT}/share/timing
    ${E3SM_ROOT}/share/pacer
    ${E3SM_ROOT}/share/util_cxx
    ${PIOC_SOURCE_DIR}
    ${PIOC_BINARY_DIR}
    ${Parmetis_INCLUDE_DIRS}
    ${OMEGA_BUILD_DIR}/external/GSW-C
)

target_compile_definitions(
    OmegaLibFlags
    INTERFACE
    OMEGA_ARCH=${OMEGA_ARCH}
    OMEGA_VECTOR_LENGTH=${OMEGA_VECTOR_LENGTH}
)

if(OMEGA_SINGLE_PRECISION)
  target_compile_definitions(
      OmegaLibFlags
      INTERFACE
      OMEGA_SINGLE_PRECISION=1
  )
endif()


target_link_options(
    OmegaLibFlags
    INTERFACE
    ${OMEGA_LINK_OPTIONS}
)

target_link_directories(
	OmegaLibFlags
	INTERFACE
	${OMEGA_BUILD_DIR}/external/GSW-C
)

target_link_libraries(
    OmegaLibFlags
    INTERFACE
    Kokkos::kokkos
    spdlog
    pioc
    yaml-cpp
    parmetis
    metis
    gptl
    gswteos-10
    pacer
    stdc++fs
    cpptrace::cpptrace
)

if(GKlib_FOUND)
    target_link_libraries(
        OmegaLibFlags
        INTERFACE
        gklib
    )
endif()

# Add source files for the library
file(GLOB_RECURSE _LIBSRC_FILES analysis/*.cpp infra/*.cpp base/*.cpp ocn/*.cpp timeStepping/*.cpp)

add_library(${OMEGA_LIB_NAME} ${_LIBSRC_FILES})

if(OMEGA_BUILD_TEST)
  add_dependencies(${OMEGA_LIB_NAME} gswteos)
endif()

target_link_libraries(
    ${OMEGA_LIB_NAME}
    PRIVATE
    OmegaLibFlags
)

# some tests require single precision version of the library
# if OmegaLib is configured to use double precision we need to
# build another version
if(OMEGA_BUILD_TEST AND NOT OMEGA_SINGLE_PRECISION)
  add_library(${OMEGA_LIB_NAME}_SP ${_LIBSRC_FILES})

  target_link_libraries(
      ${OMEGA_LIB_NAME}_SP
      PRIVATE
      OmegaLibFlags
  )

  target_compile_definitions(
      ${OMEGA_LIB_NAME}_SP
      PUBLIC
      OMEGA_SINGLE_PRECISION=1
  )
endif()


# build Omega as an E3SM component
if("${OMEGA_BUILD_MODE}" STREQUAL "E3SM")

    set(OCN_SRC
      drivers/coupled/ocn_comp_mct.F90
      drivers/coupled/omega_cxx2f_interface.cpp
      drivers/coupled/omega_f2cxx_interface.F90
      drivers/coupled/omega_cpl_indices.F90
    )

    # MOAB is optional: only present when COMP_INTERFACE=moab (buildlib_cmake
    # passes OMEGA_MOAB_ROOT in that case). The non-MOAB build must be
    # unaffected when it isn't found.
    if(OMEGA_MOAB_ROOT)
      find_package(MOAB CONFIG PATHS ${OMEGA_MOAB_ROOT} NO_DEFAULT_PATH)
    endif()
    if(MOAB_FOUND)
      list(APPEND OCN_SRC drivers/coupled/MoabInterface.cpp)
    endif()

    # create ocn lib
    add_library(ocn ${OCN_SRC})

    # Set f90 compiled modules dir
    set_target_properties(
      ocn
      PROPERTIES
      Fortran_MODULE_DIRECTORY
      ${CMAKE_BINARY_DIR}/cmake/ocn
    )
    target_include_directories(ocn PUBLIC ${CMAKE_BINARY_DIR}/cmake/ocn)

    # Link libraries
    target_link_libraries(
      ocn
      PRIVATE
      ${OMEGA_LIB_NAME}
      OmegaLibFlags
      csm_share
    )

    if(MOAB_FOUND)
      target_compile_definitions(ocn PRIVATE HAVE_MOAB)
      target_include_directories(ocn PRIVATE ${MOAB_INCLUDE_DIRS})
      target_link_libraries(ocn PRIVATE ${MOAB_LIBRARIES})
    endif()

endif()

# build Omega executable
if(OMEGA_BUILD_EXECUTABLE)

  set(EXESRC_FILES
    drivers/standalone/OceanDriver.cpp
  )

  # Create the executable target
  add_executable(${OMEGA_EXE_NAME} ${EXESRC_FILES})

  target_link_libraries(
    ${OMEGA_EXE_NAME}
    PRIVATE
    ${OMEGA_LIB_NAME}
    OmegaLibFlags
  )

  set_target_properties(${OMEGA_EXE_NAME} PROPERTIES LINKER_LANGUAGE CXX)

endif()
