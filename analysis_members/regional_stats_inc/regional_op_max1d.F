      do i = 1, solve

! 1 -> 2

      end do

      block => block % next
    end do

    ! max across processors 
    call mpas_dmpar_max_real(dminfo, out_array, flat_real)
    out_array = flat_real
  end do
