! 0d 0d

  if (type_checking) then
    ! check size of stack
    if (stack_pointer < 2) then
      call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
        'expression #' // trim(expression_names(exp_number)) // &
        ' tried to ' // trim(op_name) // ' when there ' // &
        'were less than two operands on the stack, in the RPN calculator AM')
    end if
  end if

  a_nd = stack(stack_pointer - 1) % number_of_dims
  b_nd = stack(stack_pointer) % number_of_dims

  ! call the right one
  if (a_nd == 0) then
    if (b_nd == 0) then

! 0d 1d
