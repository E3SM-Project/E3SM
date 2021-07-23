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
