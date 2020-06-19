set(PERFOBJS
  homme/src/share/prim_advection_base.F90
  homme/src/share/vertremap_base.F90
  homme/src/share/edge_mod_base.F90
  homme/src/share/derivative_mod_base.F90
  homme/src/share/bndry_mod_base.F90
  homme/src/theta-l/prim_advance_mod.F90
  homme/src/preqx/share/prim_advance_mod.F90
  homme/src/preqx/share/viscosity_preqx_base.F90
  homme/src/share/viscosity_base.F90
  homme/src/theta-l/viscosity_theta.F90
  homme/src/theta-l/eos.F90
  cam/src/physics/cam/uwshcu.F90)

set(REDUCED_PRECISION_OBJS ${CIMESRC_PATH}/share/util/shr_wv_sat_mod.F90)

set(SHR_RANDNUM_FORT_OBJS
  ${CIMESRC_PATH}/share/RandNum/src/kissvec/kissvec_mod.F90
  ${CIMESRC_PATH}/share/RandNum/src/mt19937/mersennetwister_mod.F90
  ${CIMESRC_PATH}/share/RandNum/src/dsfmt_f03/dSFMT_interface.F90
  ${CIMESRC_PATH}/share/RandNum/src/shr_RandNum_mod.F90)

set(SHR_RANDNUM_C_OBJS
  ${CIMESRC_PATH}/share/RandNum/src/dsfmt_f03/dSFMT.c
  ${CIMESRC_PATH}/share/RandNum/src/dsfmt_f03/dSFMT_utils.c
  ${CIMESRC_PATH}/share/RandNum/src/kissvec/kissvec.c)

if (NOT DEBUG)
  foreach(ITEM IN LISTS PERFOBJS)
    e3sm_add_flags("${ITEM}" "-O3 -fp-model fast -no-prec-div")
  endforeach()

  foreach(ITEM IN LISTS REDUCED_PRECISION_OBJS)
    e3sm_add_flags("${ITEM}" "-fimf-precision=low -fp-model fast")
  endforeach()

  foreach(ITEM IN LISTS SHR_RANDNUM_FORT_OBJS)
    e3sm_add_flags("${ITEM}" "-O3 -fp-model fast -no-prec-div -no-prec-sqrt -qoverride-limits")
  endforeach()

  foreach(ITEM IN LISTS SHR_RANDNUM_C_OBJS)
    e3sm_add_flags("${ITEM}" "-O3 -fp-model fast")
  endforeach()

endif()
