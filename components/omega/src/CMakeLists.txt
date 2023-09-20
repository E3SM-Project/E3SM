# build Omega

# Add source files for the library
file(GLOB _LIBSRC_FILES base/*.cpp infra/*.cpp ocn/*.cpp)

# Create the library target
add_library(${OMEGA_LIB_NAME} ${_LIBSRC_FILES})
target_link_libraries(${OMEGA_LIB_NAME} spdlog "${OMEGA_LINK_FLAGS}")

# add include directories
target_include_directories(
    ${OMEGA_LIB_NAME}
    PRIVATE
    "${OMEGA_SOURCE_DIR}/src/base"
    "${OMEGA_SOURCE_DIR}/src/infra"
     ${OMEGA_INCLUDE_DIRECTORIES}
)

target_link_directories(
    ${OMEGA_LIB_NAME}
    PRIVATE
    ${OMEGA_LINK_DIRECTORIES}
)

target_link_options(
    ${OMEGA_LIB_NAME}
    PRIVATE
    ${OMEGA_LINK_OPTIONS}
)

# include yakl cmake utility
include(${E3SM_EXTERNALS_ROOT}/YAKL/yakl_utils.cmake)

# help the library target to be built with yakl
yakl_process_target(${OMEGA_LIB_NAME})

# handles cuda cases
if (YAKL_ARCH STREQUAL "CUDA")
  set_target_properties(${OMEGA_LIB_NAME} PROPERTIES LINKER_LANGUAGE CXX)
  if (CMAKE_VERSION VERSION_GREATER "3.18.0")
    set_target_properties(${OMEGA_LIB_NAME} PROPERTIES CUDA_ARCHITECTURES OFF)
  endif()
endif()

# build Omega executable
if(OMEGA_BUILD_EXECUTABLE)
  set(EXESRC_FILES
    drivers/DrvDummy.cpp
  )
  # Create the executable target
  add_executable(${OMEGA_EXE_NAME} ${EXESRC_FILES})
  target_include_directories(
		${OMEGA_EXE_NAME}
		PRIVATE
		"${OMEGA_SOURCE_DIR}/src/base"
		"${CMAKE_CURRENT_SOURCE_DIR}/infra"
  )
  target_link_libraries(${OMEGA_EXE_NAME} ${OMEGA_LIB_NAME})

  # help the executable target to be built with yakl
  yakl_process_target(${OMEGA_EXE_NAME})

  # handles cuda cases
  if (YAKL_ARCH STREQUAL "CUDA")
    set_target_properties(${OMEGA_EXE_NAME} PROPERTIES LINKER_LANGUAGE CXX)
    if (CMAKE_VERSION VERSION_GREATER "3.18.0")
      set_target_properties(${OMEGA_EXE_NAME} PROPERTIES CUDA_ARCHITECTURES OFF)
    endif()
  endif()

endif()
