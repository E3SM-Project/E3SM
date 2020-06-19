# These routines have problems with stacksize when omp is invoked add -qsmallstack to resolve
set(SSOBJS
  cam/src/chemistry/mozart/mo_sethet.F90
  cam/src/chemistry/mozart/mo_drydep.F90)

foreach(ITEM IN LISTS SSOBJS)
  e3sm_add_flags("${ITEM}" "-qsmallstack")
endforeach()

if (compile_threaded)
  e3sm_add_flags("${CIMESRC_PATH}/share/util/shr_reprosum_mod.F90" "-qsmp=noauto:noomp")
endif()

# These routines benefit from -qnostrict without violating the bfb test
set(PERFOBJS
  homme/src/share/prim_advection_base.F90
  homme/src/share/vertremap_base.F90
  homme/src/share/edge_mod_base.F90
  homme/src/share/derivative_mod_base.F90
  homme/src/share/bndry_mod_base.F90
  homme/src/theta-l/prim_advance_mod.F90
  homme/src/preqx/share/prim_advance_mod.F90
  cam/src/physics/cam/uwshcu.F90
  cam/src/chemistry/aerosol/wetdep.F90)

#Model crashes if these files are compiled with O3(default) optimizations
set(REDUCEOPT
  cam/src/chemistry/bulk_aero/seasalt_model.F90
  cam/src/chemistry/modal_aero/seasalt_model.F90
  cam/src/chemistry/mozart/linoz_data.F90)

set(NOINLINE
  cam/src/physics/clubb/advance_xm_wpxp_module.F90
  cam/src/physics/clubb/advance_wp2_wp3_module.F90)

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
# begin
list(APPEND NOOPT_FILES ${CMAKE_CURRENT_BINARY_DIR}/buffer.F90)

foreach(ITEM IN LISTS NOINLINE)
  e3sm_add_flags("${ITEM}" "-Q!")
endforeach()
