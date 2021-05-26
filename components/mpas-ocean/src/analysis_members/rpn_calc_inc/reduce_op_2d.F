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

  do while (associated(top_iter))
    ! do operation
    top => top_iter % array
    reduced => temp_iter % array

    jend = size(top, 2)

    do j = 1, jend
      temp_iter % array(j) = &

! 2-3 break

    end do

    top_iter => top_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d2)
  end if

  ! set stack
  stack(stack_pointer) % d1 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 1
