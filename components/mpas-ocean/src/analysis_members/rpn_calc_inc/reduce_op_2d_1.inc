  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field1DReal), pointer :: temp, temp_iter
  type (field2DReal), pointer :: top_iter
  real (kind=RKIND), dimension(:,:), pointer :: top
  real (kind=RKIND), dimension(:), pointer :: reduced
  integer :: j, jend

  ! allocate a temp for result
  call create_1d_field_from_2d(stack(stack_pointer) % d2, temp)
  temp_iter => temp

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d2

  ! initial value
  reduced => temp_iter % array
  temp_iter % array = &

! 1-2 break
