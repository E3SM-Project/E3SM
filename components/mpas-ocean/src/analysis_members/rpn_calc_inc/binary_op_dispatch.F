  integer, intent(in) :: exp_number
  logical, intent(in) :: type_checking

  type (rpn_stack_value_type), dimension(MAX_STACK_SIZE), intent(inout) :: stack
  integer, intent(inout) :: stack_pointer

  integer :: a_nd, b_nd
  character (len=StrKIND) ::  op_name

! 0d 0d

  if (type_checking) then
    ! check size of stack
    if (stack_pointer < 2) then
      call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
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

    else if (b_nd == 1) then

! 0d 2d

    else

! 1d 0d

    end if
  else if (a_nd == 1) then
    if (b_nd == 0) then

! 1d 1d same

    else if (b_nd == 1) then
      if (trim(stack(stack_pointer - 1) % d1 % dimNames(1)) == &
          trim(stack(stack_pointer) % d1 % dimNames(1))) then

! 1d 1d diff

      else
        if (type_checking) then
          if (stack(stack_pointer - 1) % d1 % isDecomposed) then
            call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
            trim(op_name) // ' in expression #' // &
            trim(expression_names(exp_number)) // ' tried to operate ' // &
            'on two 1d arrays, with different dimensions, where the first ' // &
            'operand (1d array) is decomposed -- only the ' // &
            'second operand (the top of the stack) can be decomposed')
          end if
        end if

! 1d 2d first

      end if
    else
      if (trim(stack(stack_pointer - 1) % d1 % dimNames(1)) == &
          trim(stack(stack_pointer) % d2 % dimNames(1))) then

! 1d 2d second

      else
        if (type_checking) then
          if (trim(stack(stack_pointer - 1) % d1 % dimNames(1)) /= &
            trim(stack(stack_pointer) % d2 % dimNames(2))) then
            call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
            trim(op_name) // ' in expression #' // &
            trim(expression_names(exp_number)) // ' tried to operate ' // &
            'with a 1d array on a 2d array when none of the dimensions ' // &
            'match between the two arrays')
          end if
        end if

! 2d 0d

      end if
    end if
  else 
    if (b_nd == 0) then

! 2d 1d first

    else if (b_nd == 1) then
      if (trim(stack(stack_pointer - 1) % d2 % dimNames(1)) == &
          trim(stack(stack_pointer) % d1 % dimNames(1))) then

! 2d 1d second

      else
        if (type_checking) then
          if (trim(stack(stack_pointer - 1) % d2 % dimNames(2)) /= &
            trim(stack(stack_pointer) % d1 % dimNames(1))) then
            call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
            trim(op_name) // ' in expression #' // &
            trim(expression_names(exp_number)) // ' tried to operate ' // &
            'with a 1d array on a 2d array when none of the dimensions ' // &
            'match between the two arrays')
          end if
        end if

! 2d 2d

      end if
    else
      if (type_checking) then
        if ((trim(stack(stack_pointer - 1) % d2 % dimNames(1)) /= &
             trim(stack(stack_pointer) % d2 % dimNames(1))) .or. &
            (trim(stack(stack_pointer - 1) % d2 % dimNames(2)) /= &
             trim(stack(stack_pointer) % d2 % dimNames(2)))) then
           call mpas_dmpar_global_abort(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
             trim(op_name) // ' in expression #' // &
             trim(expression_names(exp_number)) // ' tried to operate ' // &
             'on two 2d arrays when their dimension names do not match')
        end if
      end if

! end

    end if
  end if
