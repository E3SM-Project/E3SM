! start -> 1d

  if (type_checking) then
    ! check size of stack
    if (stack_pointer < 1) then
      call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
        'expression #' // trim(expression_names(exp_number)) // &
        ' tried to ' // trim(op_name) // ' when there ' // &
        'were no operands on the stack, in the RPN calculator AM')
    end if
  end if

  a_nd = stack(stack_pointer) % number_of_dims

  if (type_checking) then
    if (a_nd < 1) then
      call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
        'expression #' // trim(expression_names(exp_number)) // &
        ' tried to ' // trim(op_name) // ' when the ' // &
        'operand on the stack is 0d, in the RPN calculator AM')
    end if
  end if

  ! call the right one
  if (a_nd == 1) then

! 1d -> 2d
