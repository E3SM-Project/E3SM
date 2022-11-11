string(APPEND LDFLAGS " -framework Accelerate")
set(NETCDF_PATH "$ENV{NETCDF_PATH}")
string(APPEND SLIBS " -L${NETCDF_PATH}/lib -lnetcdff -lnetcdf")
