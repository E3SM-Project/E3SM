! 1 -> 2

      end do

      block => block % next
    end do

    ! min across processors 
    call mpas_dmpar_min_real(dminfo, out_array, flat_real)
    out_array = flat_real
  end do
