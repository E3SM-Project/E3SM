set(CIMEROOT ${SCREAM_BASE_DIR}/../../cime)
list(APPEND CMAKE_MODULE_PATH ${CIMEROOT}/src/CMake)

set(GENF90 ${CIMEROOT}/src/externals/genf90/genf90.pl)
set(ENABLE_GENF90 True)
include(genf90_utils)
include(Sourcelist_utils)

set(ZM_SRCS
  ${SCREAM_BASE_DIR}/../cam/src/physics/cam/physics_utils.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/cam/scream_abortutils.F90
  zm_conv.F90
  atmosphere_deep_convection.cpp
  scream_zm_interface.F90
  zm_inputs_initializer.cpp
)

set(ZM_HEADERS
  zm.hpp
  atmosphere_deep_convection.hpp
  scream_zm_interface.hpp
)

# Add ETI source files if not on CUDA
if (NOT CUDA_BUILD)
  list(APPEND ZM_SRCS
      )
endif()

add_library(zm ${ZM_SRCS})
target_include_directories(zm PUBLIC ${SCREAM_INCLUDE_DIRS})
target_include_directories(zm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../common)
target_include_directories(zm SYSTEM PUBLIC ${SCREAM_TPL_INCLUDE_DIRS} ${CIMEROOT}/src/share/include)
set_target_properties(zm PROPERTIES
  Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
target_link_libraries(zm physics_share scream_share ${SCREAM_TPL_LIBRARIES})

