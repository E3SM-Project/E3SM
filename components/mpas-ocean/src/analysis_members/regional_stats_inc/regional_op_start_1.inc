  type (dm_info), pointer, intent(in) :: dminfo
  type (block_type), pointer, intent(in) :: start_block
  type (regional_type), intent(in) :: regions
  type (regional_variable_type), intent(in) :: variable

  real (kind=RKIND), dimension(:), allocatable :: flattened, reduced
  integer :: b, m, i, last
  integer, pointer :: solve
  type (block_type), pointer :: block
  type (mpas_pool_type), pointer :: amPool, maskPool
  real (kind=RKIND), dimension(:), pointer :: weights
  real (kind=RKIND), pointer :: weight_total
  integer, dimension(:,:), pointer :: mask
  integer :: flat_integer
  real (kind=RKIND) :: flat_real
  integer, pointer :: count_array

! 1 -> 2
