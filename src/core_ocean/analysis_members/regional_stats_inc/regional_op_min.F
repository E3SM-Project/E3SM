      do i = 1, solve

! 1 -> 2

      end do

      block => block % next
    end do

    allocate(flattened(size(out_array)))
    allocate(reduced(size(out_array)))

    ! min across processors
    flattened = reshape(out_array, shape(flattened))
    call mpas_dmpar_min_real_array(dminfo, size(flattened), flattened, reduced)
    out_array = reshape(reduced, shape(out_array))

    deallocate(reduced)
    deallocate(flattened)
  end do
