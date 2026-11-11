# Homme system tests: theta-l semi-Lagrangian tracer transport tests.
# See test/README.md for how tests are declared, and cmake/HommeTesting.cmake for the functions' docs.

if (NOT HOMME_ENABLE_COMPOSE)
  return()
endif()

homme_add_system_test (thetah-sl-test11conv-r1t2-cdr20
  EXEC          theta-l-nlev30
  RANKS         16
  NAMELISTS     thetah-sl-test11conv-r1t2-cdr20.nl
  OUTPUTS       dcmip2012_test1_3a_conv1.nc)

# The C++ version is only tested in BFB builds, to compare against F90
homme_add_f90_cxx_system_tests (thetah-sl-test11conv-r0t1-cdr30-rrm CXX_IF_BFB CXX_VS_F90
  EXEC          theta-l-nlev30
  RANKS         16
  MAX_RANKS     16
  NAMELISTS     thetah-sl-test11conv-r0t1-cdr30-rrm.nl
  MESH_FILES    ${HOMME_TEST_DATA_DIR}/mesh_files/mountain_10_x2.g
  OUTPUTS       dcmip2012_test1_3a_conv1.nc)

# QLT vertical_levels not yet supported on GPU (compose_cedr_qlt.cpp:182)
if (Kokkos_ENABLE_CUDA OR Kokkos_ENABLE_HIP)
  set (no_cxx_on_gpu NO_CXX)
endif()
homme_add_f90_cxx_system_tests (thetah-sl-dcmip16_test1pg2 CXX_IF_BFB CXX_VS_F90 ${no_cxx_on_gpu}
  EXEC          theta-l-nlev30
  RANKS         16
  NAMELISTS     thetah-sl-dcmip16_test1pg2.nl
  VCOORD_FILES  ${HOMME_TEST_DATA_DIR}/vcoord/*30*
  OUTPUTS       dcmip2016_test1_pg21.nc)

# These also check the error norms printed at the end of the run
homme_add_f90_cxx_system_tests (thetah-sl-testconv-3e
  EXEC          theta-l-nlev30
  RANKS         16
  MAX_RANKS     800
  NAMELISTS     thetah-sl-testconv-3e.nl
  OUTPUTS       dcmip2012_test1_3e_conv1.nc)
homme_add_error_norms_test (thetah-sl-testconv-3e test1_conv thetah-sl-testconv-3e.log "9.1e-2,1.1e-1,4.0e-2,8.1e-2")
homme_add_error_norms_test (thetah-sl-testconv-3e-kokkos test1_conv thetah-sl-testconv-3e-kokkos.log "9.1e-2,1.1e-1,4.0e-2,8.1e-2")

homme_add_f90_cxx_system_tests (planar-transport-a-etm
  EXEC          theta-l-nlev128-native
  RANKS         16
  MAX_RANKS     160
  NAMELISTS     planar-transport-a-etm.nl
  OUTPUTS       planar_transport_a1.nc)
homme_add_error_norms_test (planar-transport-a-etm planar_conv planar-transport-a-etm.log "3.9e-3,3.9e-3,4.1e-3")
homme_add_error_norms_test (planar-transport-a-etm-kokkos planar_conv planar-transport-a-etm-kokkos.log "3.9e-3,3.9e-3,4.1e-3")
