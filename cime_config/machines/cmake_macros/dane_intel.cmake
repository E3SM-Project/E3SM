list(APPEND CMAKE_BUILD_RPATH "$ENV{HDF5_ROOT}/lib" )
list(APPEND CMAKE_BUILD_RPATH "$ENV{NETCDF_C_PATH}/lib" )
list(APPEND CMAKE_BUILD_RPATH "$ENV{NETCDF_FORTRAN_PATH}/lib" )
list(APPEND CMAKE_BUILD_RPATH "$ENV{PNETCDF_PATH}/lib" )

# get_cmake_property(_variableNames VARIABLES)
# foreach (_variableName ${_variableNames})
#     message("${_variableName}=${${_variableName}}")
# endforeach()
