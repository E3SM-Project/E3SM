! 1d 1d diff

      else
        if (type_checking) then
          if (stack(stack_pointer - 1) % d1 % isDecomposed) then
            call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
            trim(op_name) // ' in expression #' // &
            trim(expression_names(exp_number)) // ' tried to operate ' // &
            'on two 1d arrays, with different dimensions, where the first ' // &
            'operand (1d array) is decomposed -- only the ' // &
            'second operand (the top of the stack) can be decomposed')
          end if
        end if

! 1d 2d first
