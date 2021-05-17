! 2 -> 3

    else
      nullify(dst % array)
    end if

    ! add field to this block's pool
    call mpas_pool_add_field(outpool, dst % fieldName, dst)
    call mpas_pool_add_field(all_fields, dst % fieldName, dst)

    if (present(has_vertical)) then
      has_vertical = check_vertical_dim(src % dimNames, vertical_dim)
    end if
  end if
