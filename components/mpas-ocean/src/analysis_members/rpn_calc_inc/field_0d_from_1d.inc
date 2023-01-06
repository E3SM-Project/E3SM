  type (field0DReal), pointer, intent(out) :: dst
  type (field1DReal), pointer, intent(inout) :: top

  if (mpas_threading_get_thread_num() == 0 ) then
    ! allocate the linked list for the field blocks
    allocate(dst) 
    nullify(dst % prev)
    nullify(dst % next)

    ! copy field info
    dst % fieldName = '_' // trim(top % fieldName)
    dst % isDecomposed = .false.

    dst % block => top % block
    dst % isVarArray = .false.
    dst % defaultValue = top % defaultValue
    dst % isActive = top % isActive 
    dst % hasTimeDimension = top % hasTimeDimension
    dst % sendList => top % sendList
    dst % recvList => top % recvList
    dst % copyList => top % copyList

    allocate(dst % attLists(1))

    nullify(dst % constituentNames)
  end if
