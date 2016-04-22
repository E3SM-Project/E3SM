! 1 -> 2

  last = regions % num_regions_per(regions % group_index)
  do b = 1, last
    ! get the output array, which should be the first block
    call mpas_pool_get_subpool(start_block % structs, &
      REGIONAL_STATS_POOL, amPool)
    call mpas_pool_get_array(amPool, variable % output_names(b), &
      out_array, 1)

    m = regions % groups(b, regions % group_index)

    call mpas_pool_get_array(amPool, regions % count_zerod_names(b), &
      count_array, 1)
    if (regions % function_oned == MUL_FUNC) then
      call mpas_pool_get_array(amPool, regions % weight_zerod_names(b), &
        weight_total, 1)
    end if

! 2 -> 3
