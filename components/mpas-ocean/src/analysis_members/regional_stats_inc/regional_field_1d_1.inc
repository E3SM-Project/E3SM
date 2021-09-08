  character (len=StrKIND), intent(in) :: inname, outname, elem_name
  type (mpas_pool_type), pointer, intent(inout) :: inpool, outpool, all_fields
  integer, intent(out), optional :: has_vertical
  character (len=StrKIND), intent(in), optional :: vertical_dim
  integer :: i

! 1 -> 2
