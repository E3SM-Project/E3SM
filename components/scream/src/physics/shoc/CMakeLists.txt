set(CIMEROOT ${SCREAM_BASE_DIR}/../../cime)
list(APPEND CMAKE_MODULE_PATH ${CIMEROOT}/src/CMake)

set(GENF90 ${CIMEROOT}/src/externals/genf90/genf90.pl)
set(ENABLE_GENF90 True)
include(genf90_utils)
include(Sourcelist_utils)

set(SHOC_SRCS
  shoc_f90.cpp
  shoc_functions_f90.cpp
  shoc_ic_cases.cpp
  shoc_iso_c.f90
  shoc_iso_f.f90
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam/shoc.F90
  atmosphere_macrophysics.cpp
  scream_shoc_interface.F90
  shoc_diag_second_moments_srf.cpp
  shoc_diag_second_moments_ubycond.cpp
  shoc_diag_second_moments_lbycond.cpp
  shoc_pblintd_init_pot.cpp
  shoc_pblintd_cldcheck.cpp
  shoc_check_tke.cpp
)

set(SHOC_HEADERS
  shoc.hpp
  atmosphere_macrophysics.hpp
  scream_shoc_interface.hpp
  shoc_constants.hpp
)

# Add ETI source files if not on CUDA
if (NOT CUDA_BUILD)
  list(APPEND SHOC_SRCS
    shoc_update_host_dse.cpp
    shoc_calc_shoc_varorcovar.cpp
    shoc_calc_shoc_vertflux.cpp
    shoc_linear_interp.cpp
    shoc_compute_shoc_mix_shoc_length.cpp
    shoc_compute_diag_third_shoc_moment.cpp
    shoc_clipping_diag_third_shoc_moments.cpp
    shoc_energy_integrals.cpp
    shoc_diag_second_moments_lbycond.cpp
    shoc_diag_second_moments.cpp
    shoc_diag_second_shoc_moments.cpp
    shoc_compute_brunt_shoc_length.cpp
    shoc_check_length_scale_shoc_length.cpp
    shoc_compute_conv_vel_shoc_length.cpp
    shoc_compute_l_inf_shoc_length.cpp
    shoc_diag_obklen.cpp
    shoc_compute_shr_prod.cpp
    shoc_compute_conv_time_shoc_length.cpp
    shoc_length.cpp
    shoc_energy_fixer.cpp
    shoc_compute_shoc_vapor.cpp
    shoc_update_prognostics_implicit.cpp
    shoc_diag_third_shoc_moments.cpp
    shoc_assumed_pdf.cpp
    shoc_adv_sgs_tke.cpp
    shoc_compute_tmpi.cpp
    shoc_integ_column_stability.cpp
    shoc_isotropic_ts.cpp
    shoc_dp_inverse.cpp
    shoc_main.cpp
    shoc_pblintd_height.cpp
    shoc_tridiag_solver.cpp
    shoc_pblintd_surf_temp.cpp
    shoc_pblintd_check_pblh.cpp
    shoc_pblintd.cpp
    shoc_grid.cpp
    shoc_eddy_diffusivities.cpp
    shoc_tke.cpp
  ) # SHOC ETI SRCS
endif()

add_library(shoc ${SHOC_SRCS})
target_include_directories(shoc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../share)
target_include_directories(shoc SYSTEM PUBLIC ${CIMEROOT}/src/share/include)
set_target_properties(shoc PROPERTIES Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
target_link_libraries(shoc physics_share scream_share)
target_compile_options(shoc PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${SCREAM_Fortran_FLAGS}>)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
