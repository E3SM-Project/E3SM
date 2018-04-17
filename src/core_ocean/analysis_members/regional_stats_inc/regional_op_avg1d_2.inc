! 1 -> 2

        end do
      ! real weights
      else 
        call mpas_pool_get_array(block % allFields, &
          regions % weights_oned, weights, 1)

        do i = 1, solve

! 2 -> 3
