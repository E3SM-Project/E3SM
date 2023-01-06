  type (field1DReal), pointer, intent(in) :: second, top_head
  type (field2DReal), pointer, intent(out) :: head

  type (field2DReal), pointer :: dst, prev
  type (field1DReal), pointer :: top

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
    dst % fieldName = trim(second % fieldName) // '_' // trim(top % fieldName)
    dst % isDecomposed = top % isDecomposed

    dst % block => top % block
    dst % isVarArray = second % isVarArray
    dst % defaultValue = second % defaultValue
    dst % isActive = top % isActive .and. second % isActive
    dst % hasTimeDimension = &
      top % hasTimeDimension .or. second % hasTimeDimension
    dst % sendList => top % sendList
    dst % recvList => top % recvList
    dst % copyList => top % copyList
    dst % isPersistent = top % isPersistent .or. second % isPersistent


    ! copy constitutent names if second has them
    if (associated(second % constituentNames)) then
      allocate(dst % constituentNames(size(second % constituentNames, dim=1)))
      allocate(dst % attLists(size(second % constituentNames, dim=1)))

      dst % constituentNames(:) = second % constituentNames(:)
    else
      nullify(dst % constituentNames)
      allocate(dst % attLists(1))
    end if

    dst % dimNames(1) = second % dimNames(1)
    dst % dimNames(2) = top % dimNames(1)
    dst % dimSizes(1) = second % dimSizes(1)
    dst % dimSizes(2) = top % dimSizes(1)

    ! allocate memory
    if (top % isActive .and. second % isActive) then
      allocate(dst % array(size(second % array), size(top % array)))
    else
      nullify(dst % array)
    end if

    top => top % next
    prev => dst
  end do
  end if
