# build Omega

# create library flag target for LibOmega
add_library(OmegaLibFlags INTERFACE)

target_include_directories(
    OmegaLibFlags
    INTERFACE
    ${OMEGA_SOURCE_DIR}/src/base
    ${OMEGA_SOURCE_DIR}/src/infra
    ${OMEGA_SOURCE_DIR}/src/ocn
    ${OMEGA_SOURCE_DIR}/src/timeStepping
    ${PIOC_SOURCE_DIR}
    ${PIOC_BINARY_DIR}
    ${Parmetis_INCLUDE_DIRS}
)

target_compile_definitions(
    OmegaLibFlags
    INTERFACE
    OMEGA_ARCH=${OMEGA_ARCH}
)

target_link_options(
    OmegaLibFlags
    INTERFACE
    ${OMEGA_LINK_OPTIONS}
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
)

if(GKlib_FOUND)
    target_link_libraries(
        OmegaLibFlags
        PUBLIC
        gklib
    )
endif()

# Add source files for the library
file(GLOB_RECURSE _LIBSRC_FILES infra/*.cpp base/*.cpp ocn/*.cpp timeStepping/*.cpp)

add_library(${OMEGA_LIB_NAME} ${_LIBSRC_FILES})

target_link_libraries(
    ${OMEGA_LIB_NAME}
    PRIVATE
    OmegaLibFlags
)

# build Omega executable
if(OMEGA_BUILD_EXECUTABLE)

  set(EXESRC_FILES
    drivers/DrvDummy.cpp
  )

  # Create the executable target
  add_executable(${OMEGA_EXE_NAME} ${EXESRC_FILES})

  target_link_libraries(
    ${OMEGA_EXE_NAME}
    PRIVATE
    ${OMEGA_LIB_NAME}
    OmegaLibFlags
  )

  set_target_properties(${OMEGA_EXE_NAME} PROPERTIES LINKER_LANGUAGE C)

endif()
