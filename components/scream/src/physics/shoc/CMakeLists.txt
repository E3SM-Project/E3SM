set(SHOC_SRCS
  shoc_f90.cpp
  shoc_ic_cases.cpp
  shoc_iso_c.f90
  shoc_iso_f.f90
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam/shoc.F90
  atmosphere_macrophysics.cpp
)

if (NOT SCREAM_LIB_ONLY)
  list(APPEND SHOC_SRCS
    shoc_functions_f90.cpp
    shoc_main_wrap.cpp
    ) # Add f90 bridges needed for testing
endif()

set(SHOC_HEADERS
  shoc.hpp
  atmosphere_macrophysics.hpp
  shoc_constants.hpp
)

# Add ETI source files if not on CUDA/HIP
if (NOT EAMXX_ENABLE_GPU)
  list(APPEND SHOC_SRCS
    shoc_adv_sgs_tke.cpp
    shoc_assumed_pdf.cpp
    shoc_calc_shoc_varorcovar.cpp
    shoc_calc_shoc_vertflux.cpp
    shoc_check_length_scale_shoc_length.cpp
    shoc_check_tke.cpp
    shoc_clipping_diag_third_shoc_moments.cpp
    shoc_compute_brunt_shoc_length.cpp
    shoc_compute_diag_third_shoc_moment.cpp
    shoc_compute_l_inf_shoc_length.cpp
    shoc_compute_shoc_mix_shoc_length.cpp
    shoc_compute_shoc_vapor.cpp
    shoc_compute_shr_prod.cpp
    shoc_compute_tmpi.cpp
    shoc_diag_obklen.cpp
    shoc_diag_second_moments.cpp
    shoc_diag_second_moments_lbycond.cpp
    shoc_diag_second_moments_lbycond.cpp
    shoc_diag_second_moments_srf.cpp
    shoc_diag_second_moments_ubycond.cpp
    shoc_diag_second_shoc_moments.cpp
    shoc_diag_third_shoc_moments.cpp
    shoc_dp_inverse.cpp
    shoc_eddy_diffusivities.cpp
    shoc_energy_fixer.cpp
    shoc_energy_integrals.cpp
    shoc_grid.cpp
    shoc_integ_column_stability.cpp
    shoc_isotropic_ts.cpp
    shoc_length.cpp
    shoc_linear_interp.cpp
    shoc_main.cpp
    shoc_pblintd.cpp
    shoc_pblintd_check_pblh.cpp
    shoc_pblintd_cldcheck.cpp
    shoc_pblintd_height.cpp
    shoc_pblintd_init_pot.cpp
    shoc_pblintd_surf_temp.cpp
    shoc_tke.cpp
    shoc_tridiag_solver.cpp
    shoc_update_host_dse.cpp
    shoc_update_prognostics_implicit.cpp
  ) # SHOC ETI SRCS
endif()

# List of dispatch source files if monolithic kernels are off
set(SHOC_SK_SRCS
    shoc_energy_integrals_disp.cpp
    shoc_energy_fixer_disp.cpp
    shoc_check_tke_disp.cpp
    shoc_grid_disp.cpp
    shoc_compute_shoc_vapor_disp.cpp
    shoc_diag_obklen_disp.cpp
    shoc_pblintd_disp.cpp
    shoc_length_disp.cpp
    shoc_tke_disp.cpp
    shoc_update_prognostics_implicit_disp.cpp
    shoc_diag_second_shoc_moments_disp.cpp
    shoc_diag_third_shoc_moments_disp.cpp
    shoc_assumed_pdf_disp.cpp
    shoc_update_host_dse_disp.cpp
    )

if (NOT SCREAM_DEBUG)
  # Mods to compiler flags for individual files to resolve internal compiler errors.
  if ("${SCREAM_MACHINE}" STREQUAL "cori-knl")
    set_source_files_properties (
      shoc_assumed_pdf.cpp shoc_calc_shoc_varorcovar.cpp
      PROPERTIES COMPILE_FLAGS "-O1"
    )
  endif()
endif()

if (HIP_BUILD)
  #this is needed for crusher even with small kernels
  set_source_files_properties(shoc_diag_second_shoc_moments_disp.cpp  PROPERTIES COMPILE_FLAGS -O1)
endif()

set(SHOC_LIBS "shoc")
if (SCREAM_SMALL_KERNELS)
  add_library(shoc ${SHOC_SRCS} ${SHOC_SK_SRCS})
else()
  add_library(shoc ${SHOC_SRCS})
  if (NOT SCREAM_LIBS_ONLY AND NOT SCREAM_BASELINES_ONLY)
    add_library(shoc_sk ${SHOC_SRCS} ${SHOC_SK_SRCS})
    # Always build shoc_sk with SCREAM_SMALL_KERNELS on
    target_compile_definitions(shoc_sk PUBLIC "SCREAM_SMALL_KERNELS")
    list(APPEND SHOC_LIBS "shoc_sk")
  endif()
endif()

foreach (SHOC_LIB IN LISTS SHOC_LIBS)
  set_target_properties(${SHOC_LIB} PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${SHOC_LIB}_modules
  )
  target_include_directories(${SHOC_LIB} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../share
    ${CMAKE_CURRENT_BINARY_DIR}/${SHOC_LIB}_modules
  )
  target_link_libraries(${SHOC_LIB} physics_share scream_share)
endforeach()

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
