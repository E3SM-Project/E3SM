! 1d 2d second

      else
        if (type_checking) then
          if (trim(stack(stack_pointer - 1) % d1 % dimNames(1)) /= &
            trim(stack(stack_pointer) % d2 % dimNames(2))) then
            call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
            trim(op_name) // ' in expression #' // &
            trim(expression_names(exp_number)) // ' tried to operate ' // &
            'with a 1d array on a 2d array when none of the dimensions ' // &
            'match between the two arrays')
          end if
        end if

! 2d 0d
