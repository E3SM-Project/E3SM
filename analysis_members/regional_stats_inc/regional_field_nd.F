  character (len=StrKIND), intent(in) :: inname, outname, elem_name
  type (mpas_pool_type), pointer, intent(inout) :: inpool, outpool, all_fields
  integer, intent(out), optional :: has_vertical
  character (len=StrKIND), intent(in), optional :: vertical_dim
  integer :: i

! 1 -> 2

  if (mpas_threading_get_thread_num() == 0 ) then
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
    dst % fieldName = outname ! has a different name
    dst % isDecomposed = .false. ! NOT decomposed! reducing across space

    dst % block => src % block
    dst % isVarArray = src % isVarArray
    dst % defaultValue = src % defaultValue
    dst % isActive = src % isActive
    dst % hasTimeDimension = src % hasTimeDimension
    dst % sendList => src % sendList
    dst % recvList => src % recvList
    dst % copyList => src % copyList

    ! copy constitutent names
    if (associated(src % constituentNames)) then
      allocate(dst % constituentNames( &
        size(src % constituentNames, dim=1)))
      allocate(dst % attLists(size(src % constituentNames, dim=1)))

      do i = 1, size(dst % constituentNames, dim=1)
        dst % constituentNames(i) = trim(outname) // '_' // &
          trim(src % constituentNames(i))
      end do
    else
      nullify(dst % constituentNames)
      allocate(dst % attLists(1))
    end if

    ! field values unique to non-scalars (rank-1+ arrays)
    dst % isPersistent = src % isPersistent

    i = size(src % dimNames) - 1
    dst % dimNames(1:i) = src % dimNames(1:i)
    dst % dimSizes(1:i) = src % dimSizes(1:i)

    ! allocate memory
    if (src % isActive) then
      ! create dimensions skipping the element dimension
      src_dims = shape(src % array)

! 2 -> 3

    else
      nullify(dst % array)
    end if

    ! add field to this block's pool
    call mpas_pool_add_field(outpool, dst % fieldName, dst)
    call mpas_pool_add_field(all_fields, dst % fieldName, dst)

    if (present(has_vertical)) then
      has_vertical = check_vertical_dim(src % dimNames, vertical_dim)
    end if
  end if
