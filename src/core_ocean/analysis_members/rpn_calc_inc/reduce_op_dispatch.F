  integer, intent(in) :: exp_number
  logical, intent(in) :: type_checking

  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer

  integer :: a_nd
  character (len=StrKIND) ::  op_name

! start -> 1d

  if (type_checking) then
    ! check size of stack
    if (stack_pointer < 1) then
      call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
        'expression #' // trim(expression_names(exp_number)) // &
        ' tried to ' // trim(op_name) // ' when there ' // &
        'were no operands on the stack, in the RPN calculator AM')
    end if
  end if

  a_nd = stack(stack_pointer) % number_of_dims

  if (type_checking) then
    if (a_nd < 1) then
      call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
        'expression #' // trim(expression_names(exp_number)) // &
        ' tried to ' // trim(op_name) // ' when the ' // &
        'operand on the stack is 0d, in the RPN calculator AM')
    end if
  end if

  ! call the right one
  if (a_nd == 1) then

! 1d -> 2d

  else

! 2d -> end

  end if
