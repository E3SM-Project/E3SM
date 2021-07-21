  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field1DReal), pointer :: temp, temp_iter
  type (field0DReal), pointer :: top_iter
  type (field1DReal), pointer :: second_iter
  real (kind=RKIND), pointer :: top
  real (kind=RKIND), dimension(:), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer - 1) % d1, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d0
  second_iter => stack(stack_pointer - 1) % d1

  top => top_iter % scalar 

  do while (associated(temp_iter))
    ! do operation
    second => second_iter % array 
    temp_iter % array = &

! 1-2 break

    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d0)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d1)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d1 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 1
