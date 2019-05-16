# These routines have problems with stacksize when omp is invoked add -qsmallstack to resolve
set(SSOBJS
  cam/src/chemistry/mozart/mo_sethet.F90
  cam/src/chemistry/mozart/mo_drydep.F90)

foreach(ITEM IN LISTS SSOBJS)
  set_property(SOURCE ${ITEM} APPEND_STRING PROPERTY COMPILE_FLAGS " -qsmallstack ")
endforeach()

if (compile_threaded)
  set_property(SOURCE share/util/shr_reprosum_mod.F90 APPEND_STRING PROPERTY COMPILE_FLAGS " -qsmp=noauto:noomp ")
endif()

# These routines benefit from -qnostrict without violating the bfb test
set(PERFOBJS
  homme/src/share/prim_advection_base.F90
  homme/src/share/vertremap_base.F90
  homme/src/share/edge_mod_base.F90
  homme/src/share/derivative_mod_base.F90
  homme/src/share/bndry_mod_base.F90
  homme/src/theta-l/prim_advance_mod.F90
  homme/src/pese/prim_advance_mod.F90
  homme/src/theta/prim_advance_mod.F90
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
    set_property(SOURCE ${ITEM} APPEND_STRING PROPERTY COMPILE_FLAGS " -qnostrict ")
  endforeach()

  foreach(ITEM IN LISTS REDUCEOPT)
    set_property(SOURCE ${ITEM} APPEND_STRING PROPERTY COMPILE_FLAGS " -O2 ")
  endforeach()

endif()

# These files take long time to compile with default optimization flags.
# Reducing optimization gives <1min build-times and little impact on model run time.
# begin
list(APPEND NOOPT_FILES cam/src/utils/buffer.F90)

foreach(ITEN IN LISTS NOINLINE)
  set_property(SOURCE ${ITEM} APPEND_STRING PROPERTY COMPILE_FLAGS " -Q! ")
endforeach()
