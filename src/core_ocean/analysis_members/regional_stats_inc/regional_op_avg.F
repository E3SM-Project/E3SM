      ! no weights
      if (regions % function_oned == ID_FUNC) then
        do i = 1, solve

! 1 -> 2

        end do
      ! real weights
      else 
        call mpas_pool_get_array(block % allFields, &
          regions % weights_oned, weights, 1)

        do i = 1, solve

! 2 -> 3

        end do
      end if

      block => block % next
    end do

    allocate(flattened(size(out_array)))
    allocate(reduced(size(out_array)))

    ! sum across processors and divide by total
    flattened = reshape(out_array, shape(flattened))
    call mpas_dmpar_sum_real_array(dminfo, size(flattened), flattened, reduced)
    out_array = reshape(reduced, shape(out_array))
    if (count_array > 0) then
      if (regions % function_oned == ID_FUNC) then
        out_array = out_array / count_array
      else
        out_array = out_array / weight_total
      end if
    end if

    deallocate(reduced)
    deallocate(flattened)

  end do
