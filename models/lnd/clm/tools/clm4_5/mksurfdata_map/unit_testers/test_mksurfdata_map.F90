! Run unit tests for mksurfdata_map
program mksurfdata_map_unit_tester
   use test_mkdomainMod
   use test_mkutilsMod
   use test_mkgridmapMod
   use test_mkindexmapMod
   use test_mkchecksMod
   use test_mkurbanparMod
   use test_mkncdio
   use test_mod, only : test_init, test_final

   call test_init

   ! Test mkdomainMod
   call test_domain_read_dims

   ! Test mkutilsMod
   call test_slightly_below
   call test_slightly_above

   ! Test mkgridmapMod
   call test_gridmap_areaave_default
   call test_gridmap_areaave_srcmask
   call test_gridmap_areaave_srcmask2
   call test_gridmap_areastddev

   ! Test mkindexmapMod
   call test_get_dominant_indices
   call test_filter_same
   call test_lookup_2d
   call test_lookup_2d_netcdf
   call test_which_max

   ! Test mkchecksMod
   call test_min_bad
   call test_max_bad
   
   ! Test mkurbanparMod
   call test_normalize_urbn_by_tot

   ! Test mkncdio
   call test_get_dim_lengths

   call test_final

end program mksurfdata_map_unit_tester
