! 1 -> 2

  block => start_block
  do while (associated(block))
    call mpas_pool_get_subpool(block % structs, TIME_SERIES_STATS_POOL, amPool) 
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(amPool, variable % output_names(b), &
        out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else

! 2 -> 3
