if (COMP_NAME STREQUAL elm)
  # See Land NaNs in conditionals: https://github.com/E3SM-Project/E3SM/issues/4996
  string(APPEND FFLAGS " -hfp0")
endif()
# Disable ipa and zero initialization are for other NaN isues:
# https://github.com/E3SM-Project/E3SM/pull/5208
string(APPEND FFLAGS " -hipa0 -hzero")
# -em -ef generates modulename.mod (lowercase files) to support
# Scorpio installs
string(APPEND FFLAGS " -em -ef")

string(APPEND SLIBS " -L$ENV{PNETCDF_PATH}/lib -lpnetcdf")
set(NETCDF_PATH "$ENV{NETCDF_DIR}")
set(PNETCDF_PATH "$ENV{PNETCDF_DIR}")
set(PIO_FILESYSTEM_HINTS "gpfs")
string(APPEND CXX_LIBS " -lstdc++")
