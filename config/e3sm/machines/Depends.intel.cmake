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
  homme/src/preqx/share/viscosity_preqx_base.F90
  homme/src/share/viscosity_base.F90
  homme/src/theta-l/viscosity_theta.F90
  homme/src/theta/viscosity_theta.F90
  homme/src/theta-l/eos.F90
  homme/src/theta/eos.F90
  homme/src/theta/hevi_mod.F90
  cam/src/physics/cam/uwshcu.F90)

set(REDUCED_PRECISION_OBJS share/util/shr_wv_sat_mod.F90)

set(SHR_RANDNUM_FORT_OBJS
  share/RandNum/src/kissvec/kissvec_mod.F90
  share/RandNum/src/mt19937/mersennetwister_mod.F90
  share/RandNum/src/dsfmt_f03/dSFMT_interface.F90
  share/RandNum/src/shr_RandNum_mod.F90)

set(SHR_RANDNUM_C_OBJS
  share/RandNum/src/dsfmt_f03/dSFMT.c
  share/RandNum/src/dsfmt_f03/dSFMT_utils.c
  share/RandNum/src/kissvec/kissvec.c)

if (NOT DEBUG)
  foreach(ITEM IN LISTS PERFOBJS)
    set_property(SOURCE ${ITEM} APPEND_STRING PROPERTY COMPILE_FLAGS " -O3 -fp-model fast -no-prec-div ")
  endforeach()

  foreach(ITEM IN LISTS REDUCED_PRECISION_OBJS)
    set_property(SOURCE ${ITEM} APPEND_STRING PROPERTY COMPILE_FLAGS " -fimf-precision=low -fp-model fast ")
  endforeach()

  foreach(ITEM IN LISTS SHR_RANDNUM_FORT_OBJS)
    set_property(SOURCE ${ITEM} APPEND_STRING PROPERTY COMPILE_FLAGS " -O3 -fp-model fast -no-prec-div -no-prec-sqrt -qoverride-limits ")
  endforeach()

  foreach(ITEM IN LISTS SHR_RANDNUM_C_OBJS)
    set_property(SOURCE ${ITEM} APPEND_STRING PROPERTY COMPILE_FLAGS " -O3 -fp-model fast ")
  endforeach()

endif()