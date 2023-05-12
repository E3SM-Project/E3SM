list(APPEND REDUCE_OPT_LIST
  homme/src/share/derivative_mod_base.F90
)

# Can use this flag to avoid internal compiler error for this file (with nvidia/21.11)
# Still needed with nvidia/22.5
if (NOT DEBUG)
  foreach(ITEM IN LISTS REDUCE_OPT_LIST)
    e3sm_add_flags("${ITEM}" " -Mnovect")
  endforeach()
endif()

# Use -O2 for a few files already found to benefit from increased optimization in Intel Depends file
set(PERFOBJS
  homme/src/share/prim_advection_base.F90
  homme/src/share/vertremap_base.F90
  homme/src/share/edge_mod_base.F90
  homme/src/share/bndry_mod_base.F90
  homme/src/theta-l/share/prim_advance_mod.F90
  homme/src/preqx/share/prim_advance_mod.F90
  homme/src/preqx/share/viscosity_preqx_base.F90
  homme/src/share/viscosity_base.F90
  homme/src/theta-l/share/viscosity_theta.F90
  homme/src/theta-l/share/eos.F90
  eam/src/physics/cam/uwshcu.F90)

if (NOT DEBUG)
  foreach(ITEM IN LISTS PERFOBJS)
    e3sm_remove_flags("${ITEM}" "-O1")
    e3sm_add_flags("${ITEM}" "-O2")
  endforeach()
endif()
