# build Omega

# Add source files for the library
file(GLOB _LIBSRC_FILES ocn/*.cpp)

list(JOIN OMEGA_CXX_FLAGS " " _CXX_FLAGS_EXPANDED)

# Create the library target
add_library(${OMEGA_LIB_NAME} ${_LIBSRC_FILES})
set_target_properties(${OMEGA_LIB_NAME} PROPERTIES COMPILE_FLAGS "${_CXX_FLAGS_EXPANDED}")

# add include directories
target_include_directories(${OMEGA_LIB_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/infra")

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
    drivers/drvdummy.cpp
  )
  # Create the executable target
  add_executable(${OMEGA_EXE_NAME} ${EXESRC_FILES})
  set_target_properties(${OMEGA_EXE_NAME} PROPERTIES COMPILE_FLAGS "${_CXX_FLAGS_EXPANDED}")
  target_include_directories(${OMEGA_EXE_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/infra")
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
