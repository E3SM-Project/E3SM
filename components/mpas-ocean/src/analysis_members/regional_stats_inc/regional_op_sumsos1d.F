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

    ! sum across processors
    call mpas_dmpar_sum_real(dminfo, out_array, flat_real)
    out_array = flat_real
  end do
