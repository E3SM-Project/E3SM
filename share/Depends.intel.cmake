if (NOT DEBUG)
  # Note: FFLAGS contains flags such as -fp-model consistent (and -fimf-use-svml for intel version 18)
  # The -fp-model fast flags below will effectively override other -fp-model settings.

  # shr_wv_sat_mod does not need to have better than ~0.1% precision, and benefits
  # enormously from a lower precision in the vector functions.
  e3sm_add_flags("util/shr_wv_sat_mod.F90" "-fimf-precision=low -fp-model fast")

  set(SHR_RANDNUM_FORT_SRCS
    "RandNum/src/kissvec/kissvec_mod.F90"
    "RandNum/src/mt19937/mersennetwister_mod.F90"
    "RandNum/src/dsfmt_f03/dSFMT_interface.F90"
    "RandNum/src/shr_RandNum_mod.F90")

  foreach(SHR_RANDNUM_FORT_SRC IN LISTS SHR_RANDNUM_FORT_SRCS)
    e3sm_add_flags("${SHR_RANDNUM_FORT_SRC}" "-fp-model fast -no-prec-div -no-prec-sqrt -qoverride-limits")
  endforeach()

  set(SHR_RANDNUM_C_SRCS
    "RandNum/src/dsfmt_f03/dSFMT.c"
    "RandNum/src/dsfmt_f03/dSFMT_utils.c"
    "RandNum/src/kissvec/kissvec.c")

  foreach(SHR_RANDNUM_C_SRC IN LISTS SHR_RANDNUM_C_SRCS)
    e3sm_add_flags("${SHR_RANDNUM_C_SRC}" "-fp-model fast")
  endforeach()
endif()
