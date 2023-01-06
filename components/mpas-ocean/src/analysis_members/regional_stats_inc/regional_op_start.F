  type (dm_info), pointer, intent(in) :: dminfo
  type (block_type), pointer, intent(in) :: start_block
  type (regional_type), intent(in) :: regions
  type (regional_variable_type), intent(in) :: variable

  real (kind=RKIND), dimension(:), allocatable :: flattened, reduced
  integer :: b, m, i, last
  integer, pointer :: solve
  type (block_type), pointer :: block
  type (mpas_pool_type), pointer :: amPool, maskPool
  real (kind=RKIND), dimension(:), pointer :: weights
  real (kind=RKIND), pointer :: weight_total
  integer, dimension(:,:), pointer :: mask
  integer :: flat_integer
  real (kind=RKIND) :: flat_real
  integer, pointer :: count_array

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

