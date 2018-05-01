  character (len=StrKIND), intent(in) :: inname, outname, elem_name
  type (mpas_pool_type), pointer, intent(inout) :: inpool, outpool, all_fields
  integer, intent(out), optional :: has_vertical
  character (len=StrKIND), intent(in), optional :: vertical_dim
  integer :: i

! 1 -> 2

  if (mpas_threading_get_thread_num() == 0) then
    ! get the field for this pool
    call mpas_pool_get_field(inpool, inname, src, 1)

    ! check that elem_name is in the dim list
    if (.not. check_element_dim(src % dimNames, elem_name)) then
      call mpas_dmpar_global_abort( &
        trim(CURRENT_CORE_NAME) // ' ERROR: element dimension "' // &
        elem_name // '" was not found in field "' // src % fieldName // &
        ', which is required for the current regional stats AM' // &
        'configuration.')
    end if

    ! allocate the linked list for the field blocks
    allocate(dst) 
    nullify(dst % prev)  
    nullify(dst % next)
  
    ! copy field info 
    dst % fieldName = outname ! name is different than src
    dst % isDecomposed = .false. ! NOT decomposed! reducing across space

    dst % block => src % block
    dst % isVarArray = src % isVarArray
    dst % defaultValue = src % defaultValue
    dst % isActive = src % isActive
    dst % hasTimeDimension = src % hasTimeDimension
    dst % sendList => src % sendList
    dst % recvList => src % recvList
    dst % copyList => src % copyList

    allocate(dst % attLists(1))

    ! no constituent names because dst is a scalar

    ! no isPersistent, dimNames, dimSizes because dst is scalar

    ! no need to allocate memory or loop over blocks, 
    ! as 0D has it already in dst % scalar

    ! add field to this block's pool
    call mpas_pool_add_field(outpool, dst % fieldName, dst)  
    call mpas_pool_add_field(all_fields, dst % fieldName, dst)  

    if (present(has_vertical)) then
      has_vertical = 0
    end if
  end if
