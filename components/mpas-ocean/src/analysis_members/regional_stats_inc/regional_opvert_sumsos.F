      ! no weights
      if (regions % function_twod == ID_FUNC) then
        do i = 1, solve
          do v = 1, levels

! 1 -> 2

          end do
        end do
      ! real weights
      else 
        call mpas_pool_get_array(block % allFields, &
          regions % weights_twod, weights, 1)

        do i = 1, solve
          do v = 1, levels

! 2 -> 3

          end do
        end do
      end if

      block => block % next
    end do

    allocate(flattened(size(out_array)))
    allocate(reduced(size(out_array)))

    ! sum across processors
    flattened = reshape(out_array, shape(flattened))
    call mpas_dmpar_sum_real_array(dminfo, size(flattened), flattened, reduced)
    out_array = reshape(reduced, shape(out_array))

    deallocate(reduced)
    deallocate(flattened)
  end do
