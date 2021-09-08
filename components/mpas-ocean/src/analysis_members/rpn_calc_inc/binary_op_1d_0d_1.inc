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
