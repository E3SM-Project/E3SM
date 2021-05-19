! 2 -> 3

        end do
      end if

      block => block % next
    end do

    ! sum across processors
    call mpas_dmpar_sum_real(dminfo, out_array, flat_real)
    out_array = flat_real
  end do
