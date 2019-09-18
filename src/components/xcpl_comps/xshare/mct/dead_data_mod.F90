MODULE dead_data_mod

  ! !DESCRIPTION:
  !  Declares data that's shared across the init, run, and finalize methods

  ! !USES:
  implicit none

  ! !PUBLIC TYPES:
  ! no public types

  ! !PUBLIC MEMBER FUNCTIONS:
  ! no public functions

  ! !PUBLIC DATA MEMBERS:
  integer :: dead_grid_lat    = 1   ! lat from component
  integer :: dead_grid_lon    = 2   ! lon from component
  integer :: dead_grid_area   = 3   ! area from component
  integer :: dead_grid_mask   = 4   ! mask, 0 = inactive cell
  integer :: dead_grid_frac   = 5   ! fractional area coverage
  integer :: dead_grid_aream  = 6   ! area from mapping file
  integer :: dead_grid_index  = 7   ! global index
  integer :: dead_grid_pid    = 8   ! proc id number
  integer :: dead_grid_total  = 8

END MODULE dead_data_mod
