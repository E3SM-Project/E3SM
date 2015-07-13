create_pFUnit_test(seq_map seq_map_exe
  "test_seq_map.pf" "")

target_link_libraries(seq_map_exe drv csm_share esmf_wrf_timemgr mct mpeu)
target_link_libraries(seq_map_exe ${NETCDF_LIBRARIES})
