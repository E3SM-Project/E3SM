  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field1DReal), pointer :: top_iter
  type (field1DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), dimension(:), pointer :: second
  integer :: i, j, iend, jend

  ! allocate a temp for result
  call create_2d_field_from_1ds( &
    stack(stack_pointer - 1) % d1, stack(stack_pointer) % d1, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1
  second_iter => stack(stack_pointer - 1) % d1

  second => second_iter % array
  iend = size(second)

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    jend = size(top)

    do j = 1, jend
      do i = 1, iend
        temp_iter % array(i,j) = &

! 1-2 break
