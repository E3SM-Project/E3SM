set(SHOC_SRCS
  eamxx_shoc_process_interface.cpp
)

set(SHOC_HEADERS
  shoc.hpp
  eamxx_shoc_process_interface.hpp
  shoc_constants.hpp
)

# Add ETI source files if not on CUDA/HIP
if (NOT EAMXX_ENABLE_GPU)
  list(APPEND SHOC_SRCS
    eti/shoc_adv_sgs_tke.cpp
    eti/shoc_assumed_pdf.cpp
    eti/shoc_calc_shoc_varorcovar.cpp
    eti/shoc_calc_shoc_vertflux.cpp
    eti/shoc_check_length_scale_shoc_length.cpp
    eti/shoc_check_tke.cpp
    eti/shoc_clipping_diag_third_shoc_moments.cpp
    eti/shoc_compute_brunt_shoc_length.cpp
    eti/shoc_compute_diag_third_shoc_moment.cpp
    eti/shoc_compute_l_inf_shoc_length.cpp
    eti/shoc_compute_shoc_mix_shoc_length.cpp
    eti/shoc_compute_shoc_vapor.cpp
    eti/shoc_compute_shoc_temperature.cpp
    eti/shoc_compute_shr_prod.cpp
    eti/shoc_compute_tmpi.cpp
    eti/shoc_diag_obklen.cpp
    eti/shoc_diag_second_moments.cpp
    eti/shoc_diag_second_moments_lbycond.cpp
    eti/shoc_diag_second_moments_srf.cpp
    eti/shoc_diag_second_moments_ubycond.cpp
    eti/shoc_diag_second_shoc_moments.cpp
    eti/shoc_diag_third_shoc_moments.cpp
    eti/shoc_dp_inverse.cpp
    eti/shoc_eddy_diffusivities.cpp
    eti/shoc_energy_fixer.cpp
    eti/shoc_energy_integrals.cpp
    eti/shoc_grid.cpp
    eti/shoc_integ_column_stability.cpp
    eti/shoc_isotropic_ts.cpp
    eti/shoc_length.cpp
    eti/shoc_linear_interp.cpp
    eti/shoc_main.cpp
    eti/shoc_pblintd.cpp
    eti/shoc_pblintd_check_pblh.cpp
    eti/shoc_pblintd_cldcheck.cpp
    eti/shoc_pblintd_height.cpp
    eti/shoc_pblintd_init_pot.cpp
    eti/shoc_pblintd_surf_temp.cpp
    eti/shoc_tke.cpp
    eti/shoc_tridiag_solver.cpp
    eti/shoc_update_host_dse.cpp
    eti/shoc_update_prognostics_implicit.cpp
  ) # SHOC ETI SRCS
endif()

# List of dispatch source files if monolithic kernels are off
set(SHOC_SK_SRCS
    disp/shoc_energy_integrals_disp.cpp
    disp/shoc_energy_fixer_disp.cpp
    disp/shoc_check_tke_disp.cpp
    disp/shoc_grid_disp.cpp
    disp/shoc_compute_shoc_vapor_disp.cpp
    disp/shoc_compute_shoc_temperature_disp.cpp
    disp/shoc_diag_obklen_disp.cpp
    disp/shoc_pblintd_disp.cpp
    disp/shoc_length_disp.cpp
    disp/shoc_tke_disp.cpp
    disp/shoc_update_prognostics_implicit_disp.cpp
    disp/shoc_diag_second_shoc_moments_disp.cpp
    disp/shoc_diag_third_shoc_moments_disp.cpp
    disp/shoc_assumed_pdf_disp.cpp
    disp/shoc_update_host_dse_disp.cpp
    )

if (NOT SCREAM_DEBUG)
  # Mods to compiler flags for individual files to resolve internal compiler errors.
  if ("${SCREAM_MACHINE}" STREQUAL "cori-knl")
    set_source_files_properties (
      eti/shoc_assumed_pdf.cpp
      eti/shoc_calc_shoc_varorcovar.cpp
      PROPERTIES COMPILE_FLAGS "-O1"
    )
  endif()
endif()

if (HIP_BUILD)
  #this is needed for crusher even with small kernels
  set_source_files_properties(shoc_diag_second_shoc_moments_disp.cpp  PROPERTIES COMPILE_FLAGS -O1)
endif()

set(SHOC_LIBS "shoc")
if (SCREAM_SHOC_SMALL_KERNELS)
  add_library(shoc ${SHOC_SRCS} ${SHOC_SK_SRCS})
else()
  add_library(shoc ${SHOC_SRCS})
  if (NOT SCREAM_LIB_ONLY AND NOT SCREAM_ONLY_GENERATE_BASELINES)
    add_library(shoc_sk ${SHOC_SRCS} ${SHOC_SK_SRCS})
    # Always build shoc_sk with SCREAM_SHOC_SMALL_KERNELS on
    target_compile_definitions(shoc_sk PUBLIC "SCREAM_SHOC_SMALL_KERNELS")
    list(APPEND SHOC_LIBS "shoc_sk")
  endif()
endif()
target_compile_definitions(shoc PUBLIC EAMXX_HAS_SHOC)

foreach (SHOC_LIB IN LISTS SHOC_LIBS)
  set_target_properties(${SHOC_LIB} PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${SHOC_LIB}_modules
  )
  target_include_directories(${SHOC_LIB} PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/${SHOC_LIB}_modules
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/impl
  )
  target_link_libraries(${SHOC_LIB} physics_share scream_share)
endforeach()

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()

if (TARGET eamxx_physics)
  # Add this library to eamxx_physics
  target_link_libraries(eamxx_physics INTERFACE shoc)
endif()
