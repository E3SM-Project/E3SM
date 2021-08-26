set (F90_SRC rrtmgp_interface.F90)
set (CXX_SRC
    mo_garand_atmos_io.cpp
    mo_load_cloud_coefficients.cpp
    mo_load_coefficients.cpp  
    rrtmgp_interface.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/cpp/extensions/fluxes_byband/mo_fluxes_byband_kernels.cpp
)
set (RRTMGPXX_HEADERS    
    simple_netcdf.hpp
    mo_garand_atmos_io.h
    mo_load_cloud_coefficients.h
    mo_load_coefficients.h
)

# Set compile flags for cxx source
if ("${YAKL_ARCH}" STREQUAL "CUDA")
    message(STATUS "rrtmgp_interface flags: ${YAKL_CUDA_FLAGS}")
    set_source_files_properties(${CXX_SRC} PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(${CXX_SRC} PROPERTIES COMPILE_FLAGS "-DYAKL_ARCH_CUDA --expt-extended-lambda --expt-relaxed-constexpr  ${YAKL_CUDA_FLAGS}")
else ()
    message(STATUS "rrtmgp_interface flags: ${YAKL_CXX_FLAGS}")
    set_source_files_properties(${CXX_SRC} PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS}")
endif ()

# Add library for interface code
add_library(rrtmgp_interface STATIC ${CXX_SRC})

# Libraries to link
#find_library(NETCDF_C netcdf HINTS ${NetCDF_C_PATHS}/lib)
find_library(
    NETCDF_C_LIBRARY NAMES netcdf libnetcdf 
    HINTS ${LIB_NETCDF} ${LIB_NETCDF_C}
)
target_link_libraries(rrtmgp_interface ${NETCDF_C_LIBRARY} rrtmgp yakl)

# Where to find includes
target_include_directories(rrtmgp_interface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../external/cpp)
target_include_directories(rrtmgp_interface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../external/cpp/rte)
target_include_directories(rrtmgp_interface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../external/cpp/rrtmgp)
target_include_directories(rrtmgp_interface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../external/cpp/extensions/cloud_optics)
target_include_directories(rrtmgp_interface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../external/cpp/extensions/fluxes_byband)

# Set fortran compiler flags
set_source_files_properties(${F90_SRC} PROPERTIES COMPILE_FLAGS "${CPPDEFS} ${FFLAGS}")
