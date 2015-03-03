create_pFUnit_test(seq_map seq_map_exe
  "test_seq_map.pf" "")

target_link_libraries(seq_map_exe drv csm_share esmf_wrf_timemgr mct mpeu mpi-serial)
if (CESM_NETCDF_LIBS)
  # Not all machines define CESM_NETCDF_LIBS - e.g., on yellowstone, the correct
  # libraries are implicitly included on the link line via the necessary module
  # load.
  target_link_libraries(seq_map_exe ${CESM_NETCDF_LIBS})
endif()