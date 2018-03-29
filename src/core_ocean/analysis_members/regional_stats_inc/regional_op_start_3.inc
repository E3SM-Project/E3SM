! 2 -> 3

    ! iterate over blocks
    block => start_block
    do while (associated(block))
      ! get the dimensions for this block
      if (regions % region_element == CELL_REGION) then
        call mpas_pool_get_dimension(block % dimensions, CELL_SOLVE, solve)
      else 
        call mpas_pool_get_dimension(block % dimensions, VERTEX_SOLVE, solve)
      end if

      ! get the subpools
      call mpas_pool_get_subpool(block % structs, MASK_POOL_NAME, maskPool)

      ! get the arrays
      call mpas_pool_get_array(block % allFields, &
        variable % input_name, in_array, 1)
      call mpas_pool_get_array(maskPool, regions % masking_field, mask, 1)  

