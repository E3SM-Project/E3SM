! 1-2 break

    end do

    second_iter => second_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d1)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d2)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
