!===============================================================================
!BOP ===========================================================================
!
! !MODULE: dead_data_mod -- data declaration
!
! !DESCRIPTION:
!  Declares data that's shared across the init, run, and finalize methods
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

MODULE dead_data_mod

! !USES:

   use shr_kind_mod, only : IN => SHR_KIND_IN

   implicit none

! !PUBLIC TYPES:

   ! no public types

! !PUBLIC MEMBER FUNCTIONS:

   ! no public functions

! !PUBLIC DATA MEMBERS:

  integer(IN) :: dead_grid_lat    = 1   ! lat from component
  integer(IN) :: dead_grid_lon    = 2   ! lon from component
  integer(IN) :: dead_grid_area   = 3   ! area from component
  integer(IN) :: dead_grid_mask   = 4   ! mask, 0 = inactive cell
  integer(IN) :: dead_grid_frac   = 5   ! fractional area coverage
  integer(IN) :: dead_grid_aream  = 6   ! area from mapping file
  integer(IN) :: dead_grid_index  = 7   ! global index
  integer(IN) :: dead_grid_pid    = 8   ! proc id number
  integer(IN) :: dead_grid_total  = 8    

!EOP

END MODULE dead_data_mod
!===============================================================================
