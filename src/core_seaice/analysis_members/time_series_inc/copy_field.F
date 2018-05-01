  character (len=StrKIND), intent(in) :: inname, outname
  type (mpas_pool_type), intent(inout) :: all_fields, amPool
  integer :: i

! 1 -> 2

  call mpas_pool_get_field(all_fields, inname, src, 1)
  call mpas_duplicate_field(src, dst)

  dst % fieldName = outname

  if (associated(dst % constituentNames)) then
    do i = 1, size(dst % constituentNames, dim=1)
      dst % constituentNames(i) = trim(outname) // '_' // &
        trim(dst % constituentNames(i))
    end do
  end if

  call mpas_pool_add_field(all_fields, dst % fieldName, dst)
  call mpas_pool_add_field(amPool, dst % fieldName, dst)
