# These routines have problems with stacksize when omp is invoked add -qsmallstack to resolve
set(SSOBJS
  eam/src/chemistry/mozart/mo_sethet.F90
  eam/src/chemistry/mozart/mo_drydep.F90)

foreach(ITEM IN LISTS SSOBJS)
  e3sm_add_flags("${ITEM}" "-qsmallstack")
endforeach()

# These routines benefit from -qnostrict without violating the bfb test
set(PERFOBJS
  homme/src/share/prim_advection_base.F90
  homme/src/share/vertremap_base.F90
  homme/src/share/edge_mod_base.F90
  homme/src/share/derivative_mod_base.F90
  homme/src/share/bndry_mod_base.F90
  homme/src/theta-l/share/prim_advance_mod.F90
  homme/src/preqx/share/prim_advance_mod.F90
  eam/src/physics/cam/uwshcu.F90
  eam/src/chemistry/aerosol/wetdep.F90)

#Model crashes if these files are compiled with O3(default) optimizations
set(REDUCEOPT
  eam/src/chemistry/bulk_aero/seasalt_model.F90
  eam/src/chemistry/modal_aero/seasalt_model.F90
  eam/src/chemistry/mozart/linoz_data.F90)

set(NOINLINE
  eam/src/physics/clubb/advance_xm_wpxp_module.F90
  eam/src/physics/clubb/advance_wp2_wp3_module.F90)

foreach(ITEM IN LISTS NOINLINE)
  e3sm_add_flags("${ITEM}" "-Q!")
endforeach()

if (NOT DEBUG)
  foreach(ITEM IN LISTS PERFOBJS)
    e3sm_add_flags("${ITEM}" "-qnostrict")
  endforeach()

  foreach(ITEM IN LISTS REDUCEOPT)
    e3sm_add_flags("${ITEM}" "-O2")
  endforeach()

endif()

# These files take long time to compile with default optimization flags.
# Reducing optimization gives <1min build-times and little impact on model run time.
list(APPEND NOOPT_FILES ${CMAKE_CURRENT_BINARY_DIR}/buffer.F90)

# xlf2008_r: qsmp and O0 are incompatible. Option O0 is ignored.
set(NOQSMP
  ${CMAKE_CURRENT_BINARY_DIR}/buffer.F90
  eam/src/dynamics/eul/dyn_comp.F90
  eam/src/dynamics/fv/dyn_comp.F90
  eam/src/dynamics/se/dyn_comp.F90
  eam/src/dynamics/sld/dyn_comp.F90
  eam/src/physics/cam/microp_aero.F90

  # These take >28min with -qsmp
  eam/src/physics/rrtmg/ext/rrtmg_lw/rrtmg_lw_k_g.f90
  eam/src/physics/rrtmg/ext/rrtmg_sw/rrtmg_sw_k_g.f90
)
if (compile_threaded)
  foreach(ITEM IN LISTS NOQSMP)
    e3sm_add_flags("${ITEM}" "-qnosmp")
  endforeach()
endif()

set(CPPDEFS "${CPPDEFS} -DMPAS_OPENMP_OFFLOAD")
list(APPEND MPAS_ADD_ACC_FLAGS
  ${CMAKE_BINARY_DIR}/core_seaice/shared/mpas_seaice_mesh_pool.f90
  ${CMAKE_BINARY_DIR}/core_seaice/shared/mpas_seaice_velocity_solver_variational.f90
  ${CMAKE_BINARY_DIR}/core_seaice/shared/mpas_seaice_velocity_solver.f90
)

set(FILES_NEED_OPENMP_FLAGS
  eam/src/physics/crm/samomp/ADV_MPDATA/advect_scalar.F90
  eam/src/physics/crm/samomp/ADV_MPDATA/advect_scalar2D.F90
  eam/src/physics/crm/samomp/ADV_MPDATA/advect_scalar3D.F90
  eam/src/physics/crm/samomp/ADV_MPDATA/advection.F90
  eam/src/physics/crm/samomp/accelerate_crm.F90
  eam/src/physics/crm/samomp/adams.F90
  eam/src/physics/crm/samomp/MICRO_SAM1MOM/cloud.F90
  eam/src/physics/crm/samomp/MICRO_SAM1MOM/micro_params.F90
  eam/src/physics/crm/samomp/MICRO_SAM1MOM/microphysics.F90
  eam/src/physics/crm/samomp/MICRO_SAM1MOM/precip_init.F90
  eam/src/physics/crm/samomp/MICRO_SAM1MOM/precip_proc.F90
  eam/src/physics/crm/samomp/advect2_mom_xy.F90
  eam/src/physics/crm/samomp/SGS_TKE/diffuse_mom.F90
  eam/src/physics/crm/samomp/SGS_TKE/diffuse_mom2D.F90
  eam/src/physics/crm/samomp/SGS_TKE/diffuse_mom3D.F90
  eam/src/physics/crm/samomp/SGS_TKE/diffuse_scalar.F90
  eam/src/physics/crm/samomp/SGS_TKE/diffuse_scalar2D.F90
  eam/src/physics/crm/samomp/SGS_TKE/diffuse_scalar3D.F90
  eam/src/physics/crm/samomp/SGS_TKE/sgs.F90
  eam/src/physics/crm/samomp/SGS_TKE/shear_prod2D.F90
  eam/src/physics/crm/samomp/SGS_TKE/shear_prod3D.F90
  eam/src/physics/crm/samomp/SGS_TKE/tke_full.F90
  eam/src/physics/crm/samomp/abcoefs.F90
  eam/src/physics/crm/samomp/advect2_mom_z.F90
  eam/src/physics/crm/samomp/advect_all_scalars.F90
  eam/src/physics/crm/samomp/buoyancy.F90
  eam/src/physics/crm/samomp/crm_module.F90
  eam/src/physics/crm/samomp/advect_mom.F90
  eam/src/physics/crm/samomp/atmosphere.F90
  eam/src/physics/crm/samomp/bound_duvdt.F90
  eam/src/physics/crm/samomp/bound_exchange.F90
  eam/src/physics/crm/samomp/boundaries.F90
  eam/src/physics/crm/samomp/coriolis.F90
  eam/src/physics/crm/samomp/crmtracers.F90
  eam/src/physics/crm/crm_ecpp_output_module.F90
  eam/src/physics/crm/crm_input_module.F90
  eam/src/physics/crm/samomp/crm_surface.F90
  eam/src/physics/crm/crm_output_module.F90
  eam/src/physics/crm/crm_rad_module.F90
  eam/src/physics/crm/crm_state_module.F90
  eam/src/physics/crm/samomp/damping.F90
  eam/src/physics/crm/samomp/grid.F90
  eam/src/physics/crm/samomp/diagnose.F90
  eam/src/physics/crm/samomp/params.F90
  eam/src/physics/crm/samomp/dmdf.F90
  eam/src/physics/crm/samomp/domain.F90
  eam/src/physics/crm/ecppvars.F90
  eam/src/physics/crm/samomp/fft.F90
  eam/src/physics/crm/samomp/fftpack5.F90
  eam/src/physics/crm/samomp/fftpack5_1d.F90
  eam/src/physics/crm/samomp/forcing.F90
  eam/src/physics/crm/samomp/ice_fall.F90
  eam/src/physics/crm/samomp/kurant.F90
  eam/src/physics/crm/samomp/press_grad.F90
  eam/src/physics/crm/samomp/module_ecpp_stats.F90
  eam/src/physics/crm/samomp/setparm.F90
  eam/src/physics/crm/samomp/module_ecpp_crm_driver.F90
  eam/src/physics/crm/samomp/press_rhs.F90
  eam/src/physics/crm/samomp/pressure.F90
  eam/src/physics/crm/samomp/periodic.F90
  eam/src/physics/crm/samomp/scalar_momentum.F90
  eam/src/physics/crm/samomp/random.F90
  eam/src/physics/crm/samomp/setperturb.F90
  eam/src/physics/crm/samomp/task_init.F90
  eam/src/physics/crm/samomp/task_util_NOMPI.F90
  eam/src/physics/crm/samomp/utils.F90
  eam/src/physics/crm/samomp/uvw.F90
  eam/src/physics/crm/samomp/vars.F90
  eam/src/physics/crm/samomp/zero.F90
  eam/src/physics/crm/openacc_utils.F90
  eam/src/physics/crm/samomp/sat.F90 )

foreach(ITEM IN LISTS MPAS_ADD_ACC_FLAGS)
  e3sm_add_flags("${ITEM}" "-qsmp -qoffload")
endforeach()

foreach(ITEM IN LISTS FILES_NEED_OPENMP_FLAGS)
  e3sm_add_flags("${ITEM}" "-qsmp=omp -qoffload")
endforeach()

