# build Omega

# Add source files for the library
file(GLOB _LIBSRC_FILES infra/*.cpp base/*.cpp ocn/*.cpp)

# Create the library target
add_library(${OMEGA_LIB_NAME} ${_LIBSRC_FILES})

# add compiler options
target_compile_options(
    ${OMEGA_LIB_NAME}
    PRIVATE
    "-I${OMEGA_SOURCE_DIR}/src/base"
    "-I${OMEGA_SOURCE_DIR}/src/infra"
     ${OMEGA_CXX_FLAGS}
)

# add linker options
target_link_options(
    ${OMEGA_LIB_NAME}
    PRIVATE
    ${OMEGA_LINK_OPTIONS}
)

target_link_libraries(${OMEGA_LIB_NAME} spdlog yakl pioc)

# build Omega executable
if(OMEGA_BUILD_EXECUTABLE)

  set(EXESRC_FILES
    drivers/DrvDummy.cpp
  )

  # Create the executable target
  add_executable(${OMEGA_EXE_NAME} ${EXESRC_FILES})

  target_compile_options(
		${OMEGA_EXE_NAME}
		PRIVATE
		"-L${OMEGA_SOURCE_DIR}/src/base"
		"-L${CMAKE_CURRENT_SOURCE_DIR}/infra"
  )

  target_link_libraries(${OMEGA_EXE_NAME} ${OMEGA_LIB_NAME})

endif()
