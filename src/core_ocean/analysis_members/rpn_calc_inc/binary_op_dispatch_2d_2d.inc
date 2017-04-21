! 2d 2d

      end if
    else
      if (type_checking) then
        if ((trim(stack(stack_pointer - 1) % d2 % dimNames(1)) /= &
             trim(stack(stack_pointer) % d2 % dimNames(1))) .or. &
            (trim(stack(stack_pointer - 1) % d2 % dimNames(2)) /= &
             trim(stack(stack_pointer) % d2 % dimNames(2)))) then
           call mpas_log_write(trim(MPAS_CORE_NAME) // ' ERROR: ' // &
             trim(op_name) // ' in expression #' // &
             trim(expression_names(exp_number)) // ' tried to operate ' // &
             'on two 2d arrays when their dimension names do not match')
        end if
      end if

! end
