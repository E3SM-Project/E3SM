  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer
  type (field0DReal), pointer :: temp
  type (field1DReal), pointer :: top_iter
  real (kind=RKIND), dimension(:), pointer :: top
  real (kind=RKIND), pointer :: reduced

  ! allocate a temp for result
  call create_0d_field_from_1d(stack(stack_pointer) % d1, temp)

  ! get pointers for computation
  top_iter => stack(stack_pointer) % d1

  ! initial value
  reduced => temp % scalar
  temp % scalar = &

! 1-2 break
