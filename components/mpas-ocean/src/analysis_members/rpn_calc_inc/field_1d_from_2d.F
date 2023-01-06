  type (field2DReal), pointer, intent(in) :: top_head
  type (field1DReal), pointer, intent(out) :: head

  type (field1DReal), pointer :: dst, prev
  type (field2DReal), pointer :: top

  if (mpas_threading_get_thread_num() == 0 ) then
  nullify(head)
  nullify(prev)

  top => top_head
  do while (associated(top))

    ! allocate the linked list for the field blocks
    allocate(dst) 
    if (.not. associated(head)) then
      head => dst
    end if

    if (.not. associated(prev)) then
      nullify(dst % prev)
    else
      prev % next => dst
      dst % prev => prev
    end if
    nullify(dst % next)

    ! copy field info
    dst % fieldName = '_' // trim(top % fieldName)
    dst % isDecomposed = top % isDecomposed

    dst % block => top % block
    dst % isVarArray = .false.
    dst % defaultValue = top % defaultValue
    dst % isActive = top % isActive 
    dst % hasTimeDimension = top % hasTimeDimension
    dst % sendList => top % sendList
    dst % recvList => top % recvList
    dst % copyList => top % copyList
    dst % isPersistent = top % isPersistent

    allocate(dst % attLists(1))

    nullify(dst % constituentNames)

    dst % dimNames(1) = top % dimNames(2)
    dst % dimSizes(1) = top % dimSizes(2)

    ! allocate memory
    if (top % isActive) then
      allocate(dst % array(size(top % array, 2)))
    else
      nullify(dst % array)
    end if

    top => top % next
    prev => dst
  end do
  end if
