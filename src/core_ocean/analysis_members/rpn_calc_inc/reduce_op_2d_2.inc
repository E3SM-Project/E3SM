! 1-2 break

  do while (associated(top_iter))
    ! do operation
    top => top_iter % array
    reduced => temp_iter % array

    jend = size(top, 2)

    do j = 1, jend
      temp_iter % array(j) = &

! 2-3 break
