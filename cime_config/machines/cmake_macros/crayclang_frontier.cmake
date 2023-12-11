if (compile_threaded)
  #string(APPEND CFLAGS " -fopenmp")
  string(APPEND CMAKE_Fortran_FLAGS " -fopenmp")
  string(APPEND CMAKE_CXX_FLAGS " -fopenmp")
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -fopenmp")
endif()

if (COMP_NAME STREQUAL elm)
  # See Land NaNs in conditionals: https://github.com/E3SM-Project/E3SM/issues/4996
  string(APPEND CMAKE_Fortran_FLAGS " -hfp0")
endif()
# Disable ipa and zero initialization are for other NaN isues:
# https://github.com/E3SM-Project/E3SM/pull/5208
string(APPEND CMAKE_Fortran_FLAGS " -hipa0 -hzero")
# -em -ef generates modulename.mod (lowercase files) to support
# Scorpio installs
string(APPEND CMAKE_Fortran_FLAGS " -em -ef")
