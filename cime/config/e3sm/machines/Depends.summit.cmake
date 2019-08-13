list(APPEND NOOPT_FILES
  cam/src/dynamics/eul/dyn_comp.F90
  cam/src/dynamics/fv/dyn_comp.F90
  cam/src/dynamics/se/dyn_comp.F90
  cam/src/dynamics/sld/dyn_comp.F90
  cam/src/physics/cam/microp_aero.F90)

set(FILES_NEED_CUDA_FLAGS
  homme/src/preqx_acc/bndry_mod.F90
  homme/src/preqx_acc/derivative_mod.F90
  homme/src/preqx_acc/edge_mod.F90
  homme/src/share/element_mod.F90
  homme/src/preqx_acc/element_state.F90
  homme/src/preqx_acc/openacc_utils_mod.F90
  homme/src/preqx_acc/prim_advection_mod.F90
  homme/src/share/prim_si_mod.F90
  homme/src/preqx_acc/model_init_mod.F90
  homme/src/preqx_acc/viscosity_mod.F90
  homme/src/preqx_acc/prim_driver_mod.F90
  homme/src/share/prim_driver_base.F90
  homme/src/share/physics_mod.F90
  cam/src/control/physconst.F90)
