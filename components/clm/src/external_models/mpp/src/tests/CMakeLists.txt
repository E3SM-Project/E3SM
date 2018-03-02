
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)

link_libraries(${MPP_LIBRARIES})
link_libraries(${PETSC_LIBRARIES})
link_libraries(${NETCDF_LIBRARY})
link_libraries(${NETCDF_F90_LIBRARIES})

include(add_mpp_unit_test)
add_mpp_unit_test(test_eos_constant_density test_eos_constant_density.F90)
add_mpp_unit_test(test_eos_tgdp01_density test_eos_tgdp01_density.F90)
add_mpp_unit_test(test_eos_ifc67_density test_eos_ifc67_density.F90)

