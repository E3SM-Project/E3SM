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

  do while (associated(top_iter))
    ! do operation
    top => top_iter % array

    temp % scalar = &

! 2-3 break

    top_iter => top_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d1)
  end if

  ! set stack
  stack(stack_pointer) % d0 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 0
