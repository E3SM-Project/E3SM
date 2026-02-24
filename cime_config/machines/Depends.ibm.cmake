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
