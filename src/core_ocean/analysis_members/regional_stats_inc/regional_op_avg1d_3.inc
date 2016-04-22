! 2 -> 3

        end do
      end if

      block => block % next
    end do

    ! sum across processors and divide by total
    call mpas_dmpar_sum_real(dminfo, out_array, flat_real)
    out_array = flat_real
    if (count_array > 0) then
      if (regions % function_oned == ID_FUNC) then
        out_array = out_array / count_array
      else
        out_array = out_array / weight_total
      end if
    end if
  end do
