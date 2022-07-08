! 2 -> 3

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
