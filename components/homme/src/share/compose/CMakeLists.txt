include_directories (
  ${CMAKE_CURRENT_BINARY_DIR}
  ${HOMME_SOURCE_DIR}/utils/cime/CIME/non_py/src/timing)
if (ENABLE_HORIZ_OPENMP)
  add_definitions(-DHORIZ_OPENMP)
endif ()
if (ENABLE_COLUMN_OPENMP)
  add_definitions(-DCOLUMN_OPENMP)
endif ()
if (COMPOSE_PORT)
  add_definitions(-DCOMPOSE_PORT)
endif ()
if (NOT HOMMEXX_MPI_ON_DEVICE)
  message("compose> MPI on host")
  add_definitions(-DCOMPOSE_MPI_ON_HOST)
endif ()
add_library (${COMPOSE_LIBRARY}
  compose_test.cpp
  compose_homme.cpp
  compose_hommexx.cpp
  compose_cedr.cpp
  compose_slmm.cpp
  compose_slmm_util.cpp
  compose_slmm_siqk.cpp
  compose_slmm_departure_point.cpp
  compose_slmm_advecter.cpp
  compose_slmm_islmpi.cpp
  compose_slmm_islmpi_adp.cpp
  compose_slmm_islmpi_comm.cpp
  compose_slmm_islmpi_pack.cpp
  compose_slmm_islmpi_q.cpp
  compose_slmm_islmpi_qextrema.cpp
  compose_slmm_islmpi_interpolate.cpp
  compose_slmm_islmpi_step.cpp
  compose_cedr_sl_run_global.cpp
  compose_cedr_sl_run_local.cpp
  compose_cedr_sl_run_check.cpp
  compose_cedr_qlt.cpp
  compose_cedr_caas.cpp
  compose_slmm_islmpi_calc_trajectory.cpp
  cedr_util.cpp
  cedr_mpi.cpp
  cedr_local.cpp
  cedr_bfb_tree_allreduce.cpp
  cedr_tree.cpp
  cedr_caas.cpp
  cedr_qlt.cpp
  cedr_test_randomized.cpp)

target_link_libraries(${COMPOSE_LIBRARY} Kokkos::kokkos)
