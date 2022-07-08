  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  integer :: b
  type (block_type), pointer :: block
  type (mpas_pool_type), pointer :: amPool

! 1 -> 2
