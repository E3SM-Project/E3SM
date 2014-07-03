
module pmgrid

!----------------------------------------------------------------------- 
! 
! Purpose: Parameters and variables related to the dynamics grid
! 
! Author: 
! 
!-----------------------------------------------------------------------

   implicit none

   public 

   integer, parameter :: plon   = PLON                     ! number of longitudes
   integer, parameter :: plev   = PLEV                     ! number of vertical levels
   integer, parameter :: plat   = PLAT                     ! number of latitudes
   integer, parameter :: plevp  = plev + 1                 ! plev + 1
   integer, parameter :: plnlv  = plon*plev                ! Length of multilevel field slice
!
   integer :: beglat     ! beg. index for latitudes owned by a given proc
   integer :: endlat     ! end. index for latitudes owned by a given proc
   integer :: begirow    ! beg. index for latitude pairs owned by a given proc
   integer :: endirow    ! end. index for latitude pairs owned by a given proc
   integer :: numlats    ! number of latitudes owned by a given proc
   logical :: dyndecomp_set = .false. ! flag indicates dynamics grid has been set for history
!
#if ( ! defined SPMD )
   parameter (beglat   = 1)
   parameter (endlat   = plat)
   parameter (begirow  = 1)
   parameter (endirow  = plat/2)
   parameter (numlats  = plat)
#endif
end module pmgrid

