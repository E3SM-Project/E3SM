! 1-2 break

    top_iter => top_iter % next
    temp_iter => temp_iter % next
  end do

  ! clean up old
  if (stack(stack_pointer) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer) % d2)
  end if

  if (stack(stack_pointer - 1) % symbol_type == IS_TEMPORARY) then
    call mpas_deallocate_field(stack(stack_pointer - 1) % d0)
  end if

  ! set stack
  stack_pointer = stack_pointer - 1
  stack(stack_pointer) % d2 => temp
  stack(stack_pointer) % symbol_type = IS_TEMPORARY
  stack(stack_pointer) % number_of_dims = 2
