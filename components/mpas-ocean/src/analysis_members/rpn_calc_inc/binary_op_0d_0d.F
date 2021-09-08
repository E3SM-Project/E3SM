  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field0DReal), pointer :: temp, temp_iter
  type (field0DReal), pointer :: top_iter
  type (field0DReal), pointer :: second_iter
  real (kind=RKIND), pointer :: top
  real (kind=RKIND), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d0, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d0
  second_iter => stack(stack_pointer - 1) % d0

  ! do operation
  top => top_iter % scalar
  second => second_iter % scalar
  temp_iter % scalar = &

! 1-2 break

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d0)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d0)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d0 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 0
