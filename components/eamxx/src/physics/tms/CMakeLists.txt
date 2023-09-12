set(TMS_SRCS
  tms_iso_c.f90
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam/trb_mtn_stress.F90
  eamxx_tms_process_interface.cpp
)

if (NOT SCREAM_LIB_ONLY)
  list(APPEND TMS_SRCS
    tms_functions_f90.cpp
    )
endif()

set(TMS_HEADERS
  eamxx_tms_process_interface.hpp
)

# Add ETI source files if not on CUDA/HIP
if (NOT EAMXX_ENABLE_GPU)
  list(APPEND TMS_SRCS
    eti/compute_tms.cpp
  ) # TMS ETI SRCS
endif()

set(TMS_LIBS "tms")
add_library(tms ${TMS_SRCS})
target_compile_definitions(tms PUBLIC EAMXX_HAS_TMS)

foreach (TMS_LIB IN LISTS TMS_LIBS)
  set_target_properties(${TMS_LIB} PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${TMS_LIB}_modules
  )
  target_include_directories(${TMS_LIB} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../share
    ${CMAKE_CURRENT_BINARY_DIR}/${TMS_LIB}_modules
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/impl
  )
  target_link_libraries(${TMS_LIB} physics_share scream_share)
endforeach()

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
