  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field2DReal), pointer :: temp, temp_iter
  type (field2DReal), pointer :: top_iter
  type (field1DReal), pointer :: second_iter
  real (kind=RKIND), dimension(:,:), pointer :: top
  real (kind=RKIND), dimension(:), pointer :: second
  integer :: j, jend

  ! allocate a temp for result
  call mpas_duplicate_field(stack(stack_pointer) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d2
  second_iter => stack(stack_pointer - 1) % d1

  second => second_iter % array

  do while (associated(temp_iter))
    ! do operation
    top => top_iter % array
    jend = size(top, 2)

    do j = 1, jend
      temp_iter % array(:,j) = &

! 1-2 break
