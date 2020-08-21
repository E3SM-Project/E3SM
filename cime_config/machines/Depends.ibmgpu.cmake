
set(CPPDEFS "${CPPDEFS} -DMPAS_OPENMP_OFFLOAD")
foreach(ITEM IN LISTS MPAS_ADD_ACC_FLAGS)
  e3sm_add_flags("${ITEM}" "-qsmp -qoffload")
endforeach()

set(MMF_FILES_NEED_OPENMP_FLAGS
  cam/src/physics/crm/sam/ADV_MPDATA/advect_scalar.F90
  cam/src/physics/crm/sam/ADV_MPDATA/advect_scalar2D.F90
  cam/src/physics/crm/sam/ADV_MPDATA/advect_scalar3D.F90
  cam/src/physics/crm/sam/ADV_MPDATA/advection.F90
  cam/src/physics/crm/sam/accelerate_crm.F90
  cam/src/physics/crm/sam/adams.F90
  cam/src/physics/crm/sam/MICRO_SAM1MOM/cloud.F90
  cam/src/physics/crm/sam/MICRO_SAM1MOM/micro_params.F90
  cam/src/physics/crm/sam/MICRO_SAM1MOM/microphysics.F90
  cam/src/physics/crm/sam/MICRO_SAM1MOM/precip_init.F90
  cam/src/physics/crm/sam/MICRO_SAM1MOM/precip_proc.F90
  cam/src/physics/crm/sam/advect2_mom_xy.F90
  cam/src/physics/crm/sam/SGS_TKE/diffuse_mom.F90
  cam/src/physics/crm/sam/SGS_TKE/diffuse_mom2D.F90
  cam/src/physics/crm/sam/SGS_TKE/diffuse_mom3D.F90
  cam/src/physics/crm/sam/SGS_TKE/diffuse_scalar.F90
  cam/src/physics/crm/sam/SGS_TKE/diffuse_scalar2D.F90
  cam/src/physics/crm/sam/SGS_TKE/diffuse_scalar3D.F90
  cam/src/physics/crm/sam/SGS_TKE/sgs.F90
  cam/src/physics/crm/sam/SGS_TKE/shear_prod2D.F90
  cam/src/physics/crm/sam/SGS_TKE/shear_prod3D.F90
  cam/src/physics/crm/sam/SGS_TKE/tke_full.F90
  cam/src/physics/crm/sam/abcoefs.F90
  cam/src/physics/crm/sam/advect2_mom_z.F90
  cam/src/physics/crm/sam/advect_all_scalars.F90
  cam/src/physics/crm/sam/buoyancy.F90
  cam/src/physics/crm/sam/crm_module.F90
  cam/src/physics/crm/sam/advect_mom.F90
  cam/src/physics/crm/sam/atmosphere.F90
  cam/src/physics/crm/sam/bound_duvdt.F90
  cam/src/physics/crm/sam/bound_exchange.F90
  cam/src/physics/crm/sam/boundaries.F90
  cam/src/physics/crm/sam/coriolis.F90
  cam/src/physics/crm/sam/crmtracers.F90
  cam/src/physics/crm/sam/crm_ecpp_output_module.F90
  cam/src/physics/crm/sam/crm_input_module.F90
  cam/src/physics/crm/sam/crm_surface.F90
  cam/src/physics/crm/sam/crm_output_module.F90
  cam/src/physics/crm/sam/crm_rad_module.F90
  cam/src/physics/crm/sam/crm_state_module.F90
  cam/src/physics/crm/sam/damping.F90
  cam/src/physics/crm/sam/grid.F90
  cam/src/physics/crm/sam/diagnose.F90
  cam/src/physics/crm/sam/params.F90
  cam/src/physics/crm/sam/dmdf.F90
  cam/src/physics/crm/sam/domain.F90
  cam/src/physics/crm/sam/ecppvars.F90
  cam/src/physics/crm/sam/fft.F90
  cam/src/physics/crm/sam/fftpack5.F90
  cam/src/physics/crm/sam/fftpack5_1d.F90
  cam/src/physics/crm/sam/forcing.F90
  cam/src/physics/crm/sam/ice_fall.F90
  cam/src/physics/crm/sam/kurant.F90
  cam/src/physics/crm/sam/press_grad.F90
  cam/src/physics/crm/sam/module_ecpp_stats.F90
  cam/src/physics/crm/sam/setparm.F90
  cam/src/physics/crm/sam/module_ecpp_crm_driver.F90
  cam/src/physics/crm/sam/press_rhs.F90
  cam/src/physics/crm/sam/pressure.F90
  cam/src/physics/crm/sam/periodic.F90
  cam/src/physics/crm/sam/scalar_momentum.F90
  cam/src/physics/crm/sam/random.F90
  cam/src/physics/crm/sam/setperturb.F90
  cam/src/physics/crm/sam/task_init.F90
  cam/src/physics/crm/sam/task_util_NOMPI.F90
  cam/src/physics/crm/sam/utils.F90
  cam/src/physics/crm/sam/uvw.F90
  cam/src/physics/crm/sam/vars.F90
  cam/src/physics/crm/sam/zero.F90
  cam/src/physics/crm/sam/openacc_utils.F90
  cam/src/physics/crm/sam/crm_data_module.F90
  cam/src/physics/crm/sam/sat.F90 )

foreach(ITEM IN LISTS MMF_FILES_NEED_OPENMP_FLAGS)
  e3sm_add_flags("${ITEM}" "-qsmp=omp -qoffload")
endforeach()

