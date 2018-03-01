! Run unit tests for mksurfdata_map
program mksurfdata_map_unit_tester
   use test_mkutilsMod
   use test_mkindexmapMod
   use test_mkurbanparDomMod
   use test_mkncdio
   use test_mod, only : test_init, test_final

   call test_init

   call test_slightly_below
   call test_slightly_above

   call test_get_dominant_indices
   call test_filter_same
   call test_lookup_2d
   call test_lookup_2d_netcdf
   call test_which_max

   call test_mkurban_dominant_density

   call test_get_dim_lengths

   call test_final

end program mksurfdata_map_unit_tester
