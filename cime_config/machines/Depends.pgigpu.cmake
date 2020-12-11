list(APPEND NOOPT_FILES
  eam/src/dynamics/eul/dyn_comp.F90
  eam/src/dynamics/fv/dyn_comp.F90
  eam/src/dynamics/se/dyn_comp.F90
  eam/src/dynamics/sld/dyn_comp.F90
  eam/src/physics/cam/microp_aero.F90)

set(FILES_NEED_OPENACC_FLAGS
  eam/src/physics/crm/sam/ADV_MPDATA/advect_scalar.F90
  eam/src/physics/crm/sam/ADV_MPDATA/advect_scalar2D.F90
  eam/src/physics/crm/sam/ADV_MPDATA/advect_scalar3D.F90
  eam/src/physics/crm/sam/ADV_MPDATA/advection.F90
  eam/src/physics/crm/sam/accelerate_crm.F90
  eam/src/physics/crm/sam/adams.F90
  eam/src/physics/crm/sam/MICRO_SAM1MOM/cloud.F90
  eam/src/physics/crm/sam/MICRO_SAM1MOM/micro_params.F90
  eam/src/physics/crm/sam/MICRO_SAM1MOM/microphysics.F90
  eam/src/physics/crm/sam/MICRO_SAM1MOM/precip_init.F90
  eam/src/physics/crm/sam/MICRO_SAM1MOM/precip_proc.F90
  eam/src/physics/crm/sam/advect2_mom_xy.F90
  eam/src/physics/crm/sam/SGS_TKE/diffuse_mom.F90
  eam/src/physics/crm/sam/SGS_TKE/diffuse_mom2D.F90
  eam/src/physics/crm/sam/SGS_TKE/diffuse_mom3D.F90
  eam/src/physics/crm/sam/SGS_TKE/diffuse_scalar.F90
  eam/src/physics/crm/sam/SGS_TKE/diffuse_scalar2D.F90
  eam/src/physics/crm/sam/SGS_TKE/diffuse_scalar3D.F90
  eam/src/physics/crm/sam/SGS_TKE/sgs.F90
  eam/src/physics/crm/sam/SGS_TKE/shear_prod2D.F90
  eam/src/physics/crm/sam/SGS_TKE/shear_prod3D.F90
  eam/src/physics/crm/sam/SGS_TKE/tke_full.F90
  eam/src/physics/crm/sam/abcoefs.F90
  eam/src/physics/crm/sam/advect2_mom_z.F90
  eam/src/physics/crm/sam/advect_all_scalars.F90
  eam/src/physics/crm/sam/buoyancy.F90
  eam/src/physics/crm/sam/crm_module.F90
  eam/src/physics/crm/sam/advect_mom.F90
  eam/src/physics/crm/sam/atmosphere.F90
  eam/src/physics/crm/sam/bound_duvdt.F90
  eam/src/physics/crm/sam/bound_exchange.F90
  eam/src/physics/crm/sam/boundaries.F90
  eam/src/physics/crm/sam/coriolis.F90
  eam/src/physics/crm/sam/crmtracers.F90
  eam/src/physics/crm/sam/crm_ecpp_output_module.F90
  eam/src/physics/crm/crm_input_module.F90
  eam/src/physics/crm/sam/crm_surface.F90
  eam/src/physics/crm/crm_output_module.F90
  eam/src/physics/crm/crm_rad_module.F90
  eam/src/physics/crm/crm_state_module.F90
  eam/src/physics/crm/sam/damping.F90
  eam/src/physics/crm/sam/grid.F90
  eam/src/physics/crm/sam/diagnose.F90
  eam/src/physics/crm/sam/params.F90
  eam/src/physics/crm/sam/dmdf.F90
  eam/src/physics/crm/sam/domain.F90
  eam/src/physics/crm/sam/ecppvars.F90
  eam/src/physics/crm/sam/fft.F90
  eam/src/physics/crm/sam/fftpack5.F90
  eam/src/physics/crm/sam/fftpack5_1d.F90
  eam/src/physics/crm/sam/forcing.F90
  eam/src/physics/crm/sam/ice_fall.F90
  eam/src/physics/crm/sam/kurant.F90
  eam/src/physics/crm/sam/press_grad.F90
  eam/src/physics/crm/sam/module_ecpp_stats.F90
  eam/src/physics/crm/sam/setparm.F90
  eam/src/physics/crm/sam/module_ecpp_crm_driver.F90
  eam/src/physics/crm/sam/press_rhs.F90
  eam/src/physics/crm/sam/pressure.F90
  eam/src/physics/crm/sam/periodic.F90
  eam/src/physics/crm/sam/scalar_momentum.F90
  eam/src/physics/crm/sam/random.F90
  eam/src/physics/crm/sam/setperturb.F90
  eam/src/physics/crm/sam/task_init.F90
  eam/src/physics/crm/sam/task_util_NOMPI.F90
  eam/src/physics/crm/sam/utils.F90
  eam/src/physics/crm/sam/uvw.F90
  eam/src/physics/crm/sam/vars.F90
  eam/src/physics/crm/sam/zero.F90
  eam/src/physics/crm/openacc_utils.F90
  eam/src/physics/crm/sam/sat.F90 )

# add accelerator/gpu flags for MPAS files
set(CPPDEFS "${CPPDEFS} -DMPAS_OPENACC")
list(APPEND MPAS_ADD_ACC_FLAGS
  ${CMAKE_BINARY_DIR}/core_ocean/shared/mpas_ocn_equation_of_state_jm.f90
  ${CMAKE_BINARY_DIR}/core_ocean/shared/mpas_ocn_mesh.f90
  ${CMAKE_BINARY_DIR}/core_ocean/shared/mpas_ocn_surface_bulk_forcing.f90
  ${CMAKE_BINARY_DIR}/core_ocean/shared/mpas_ocn_surface_land_ice_fluxes.f90
  ${CMAKE_BINARY_DIR}/core_ocean/shared/mpas_ocn_tendency.f90
  ${CMAKE_BINARY_DIR}/core_ocean/shared/mpas_ocn_vel_forcing_explicit_bottom_drag.f90
  ${CMAKE_BINARY_DIR}/core_ocean/shared/mpas_ocn_vel_forcing_surface_stress.f90
  ${CMAKE_BINARY_DIR}/core_ocean/shared/mpas_ocn_vel_hadv_coriolis.f90
  ${CMAKE_BINARY_DIR}/core_ocean/shared/mpas_ocn_vel_hmix_del2.f90
  ${CMAKE_BINARY_DIR}/core_ocean/shared/mpas_ocn_vel_hmix_del4.f90
  ${CMAKE_BINARY_DIR}/core_ocean/shared/mpas_ocn_vel_hmix_leith.f90
  ${CMAKE_BINARY_DIR}/core_ocean/shared/mpas_ocn_vel_pressure_grad.f90
  ${CMAKE_BINARY_DIR}/core_ocean/shared/mpas_ocn_vel_vadv.f90
  # seaice
  ${CMAKE_BINARY_DIR}/core_seaice/shared/mpas_seaice_mesh_pool.f90
  ${CMAKE_BINARY_DIR}/core_seaice/shared/mpas_seaice_velocity_solver_variational.f90
  ${CMAKE_BINARY_DIR}/core_seaice/shared/mpas_seaice_velocity_solver.f90
)

foreach(ITEM IN LISTS FILES_NEED_OPENACC_FLAGS)
  e3sm_add_flags("${ITEM}" "-Minline -acc -ta=tesla:ccall,fastmath,loadcache:L1,unroll,fma,managed,ptxinfo -Mcuda -Minfo=accel")
endforeach()

foreach(ITEM IN LISTS MPAS_ADD_ACC_FLAGS)
  e3sm_add_flags("${ITEM}" "-Minline -acc -ta=tesla:ccall,fastmath,loadcache:L1,unroll,fma,deepcopy,nonvvm -Minfo=accel")
endforeach()

