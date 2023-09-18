set(ZM_SRCS
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam/physics_utils.F90
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam/scream_abortutils.F90
  zm_conv.F90
  eamxx_zm_process_interface.cpp
  scream_zm_interface.F90
)

set(ZM_HEADERS
  zm.hpp
  eamxx_zm_process_interface.hpp
  scream_zm_interface.hpp
)

# Add ETI source files if not on CUDA/HIP
if (NOT EAMXX_ENABLE_GPU)
  list(APPEND ZM_SRCS
      )
endif()

add_library(zm ${ZM_SRCS})
set_target_properties(zm PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)
target_include_directories(zm PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
)
target_link_libraries(zm physics_share scream_share)
