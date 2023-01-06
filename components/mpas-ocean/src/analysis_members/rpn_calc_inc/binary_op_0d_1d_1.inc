  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field1DReal), pointer :: temp, temp_iter
  type (field1DReal), pointer :: top_iter
  type (field0DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), pointer :: second

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d1, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1
  second_iter => stack(stack_pointer - 1) % d0

  second => second_iter % scalar

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    temp_iter % array = &

! 1-2 break
