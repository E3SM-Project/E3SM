# Create tms library
add_library(tms eamxx_tms_process_interface.cpp)
target_link_libraries(tms physics_share scream_share)
target_compile_definitions(tms PUBLIC EAMXX_HAS_TMS)
target_include_directories(tms PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/impl
)

# Add ETI source files if not on CUDA/HIP
if (NOT EAMXX_ENABLE_GPU)
  target_sources(tms PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/eti/compute_tms.cpp
  )
endif()

if (NOT SCREAM_LIB_ONLY)
  # For testing, add some more sources and modules directory
  target_sources (tms PRIVATE
    ${SCREAM_BASE_DIR}/../eam/src/physics/cam/trb_mtn_stress.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/tms_iso_c.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/tms_functions_f90.cpp
  )
  set_target_properties(tms PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tms_modules
  )
  target_include_directories (tms PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/tms_modules)

  add_subdirectory(tests)
endif()

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE tms)
endif()
