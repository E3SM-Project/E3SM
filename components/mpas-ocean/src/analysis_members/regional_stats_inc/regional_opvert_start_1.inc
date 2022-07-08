  type (dm_info), pointer, intent(in) :: dminfo
  type (block_type), pointer, intent(in) :: start_block
  type (regional_type), intent(in) :: regions
  type (regional_variable_type), intent(in) :: variable
  integer, intent(in) :: levels

  real (kind=RKIND), dimension(:), allocatable :: flattened, reduced
  integer :: b, m, i, last, v
  integer, pointer :: solve
  type (block_type), pointer :: block
  type (mpas_pool_type), pointer :: amPool, maskPool
  real (kind=RKIND), dimension(:,:), pointer :: weights
  real (kind=RKIND), dimension(:), pointer :: weight_total
  integer, dimension(:,:), pointer :: mask, vertical_mask
  integer, dimension(:), pointer :: count_array
  integer, dimension(:), allocatable :: copy_count

! 1 -> 2
