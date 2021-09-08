      ! no weights
      do i = 1, solve
        do v = 1, levels

! 1 -> 2

        end do
      end do
      block => block % next
    end do

    allocate(flattened(size(out_array)))
    allocate(reduced(size(out_array)))

    ! max across processors
    flattened = reshape(out_array, shape(flattened))
    call mpas_dmpar_max_real_array(dminfo, size(flattened), flattened, reduced)
    out_array = reshape(reduced, shape(out_array))

    deallocate(reduced)
    deallocate(flattened)
  end do
