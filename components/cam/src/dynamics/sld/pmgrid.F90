
module pmgrid

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize grid point resolution parameters
! 
! Author: 
! 
!-----------------------------------------------------------------------
  implicit none

  private

  public plev, plnlv, plevp, plat, plon
  public begirow, endirow, beglat, endlat
  public numlats, dyndecomp_set

! Grid point resolution parameters

  integer, parameter :: plon       = PLON  ! number of longitudes
  integer, parameter :: plev       = PLEV  ! number of vertical levels
  integer, parameter :: plat       = PLAT  ! number of latitudes
  integer, parameter :: plevp      = plev + 1 ! plev + 1
  integer, parameter :: plnlv      = plon*plev     ! Length of multilevel field slice


  integer begirow    ! beg. index for lat pairs owned by a proc
  integer endirow    ! end. index for lat pairs owned by a proc
  integer beglat     ! beg. index for lats owned by a given proc
  integer endlat     ! end. index for lats owned by a given proc
  integer numlats    ! number of lats owned by a given proc

  logical :: dyndecomp_set = .false. ! flag indicates dynamics grid has been set
#if ( ! defined SPMD )
  parameter (begirow    = 1)
  parameter (endirow    = plat/2)
  parameter (beglat     = 1)
  parameter (endlat     = plat)
  parameter (numlats    = plat)
#endif
end module pmgrid


