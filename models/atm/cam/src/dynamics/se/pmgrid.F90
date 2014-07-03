module pmgrid
  use edge_mod, only : EdgeBuffer_t
  use dimensions_mod, only : nnodes,npart,nmpi_per_node
  use reduction_mod, only : ReductionBuffer_ordered_1d_t
  
!----------------------------------------------------------------------- 
! 
! Purpose: Parameters and variables related to the dynamics grid
! 
! Author: Jim Edwards
! 
! PLON and PLAT do not correspond to the number of latitudes and longitudes in
! this version of dynamics. 
! 
! 
!-----------------------------------------------------------------------
   save

   integer, parameter :: plon   = 1                     ! number of longitudes
   integer, parameter :: plev   = PLEV                     ! number of vertical levels
   integer, parameter :: plat   = 1                     ! number of latitudes


   integer, parameter :: plevp  = plev + 1                 ! plev + 1


   type (EdgeBuffer_t) :: edge_3    ! 3 layer          edge buffer  (shared)
   type (EdgeBuffer_t) :: edge_1lp1 ! nlev+1 layer     edge buffer  (shared)
   type (EdgeBuffer_t) :: edge_3lp1 ! 3*(nlev)+1 layer edge buffer  (shared)

   type (ReductionBuffer_ordered_1d_t)  :: red    ! reduction buffer for cg
! End additions
!
!JPE: i1 and j1, nxpt, jintmx are defined as placeholders - these values are needed
!     to satisfy code in advection/slt 
!
#ifndef STAGGERED
   integer, parameter :: nxpt=0
   integer, parameter :: jintmx=0
   integer, parameter :: i1     = 1 
   integer, parameter :: j1     = 1 
#endif
   integer, parameter :: numbnd = 0          ! no.of latitudes passed N and S of forecast lat

!
   integer :: beglat     ! beg. index for latitudes owned by a given proc
   integer :: endlat     ! end. index for latitudes owned by a given proc
   integer :: beglatex   ! extended grid beglat
   integer :: endlatex   ! extended grid endlat
   integer :: numlats    ! number of latitudes owned by a given proc
   logical :: dyndecomp_set = .false. ! flag indicates dynamics grid has been set
!
#if ( ! defined SPMD )
   parameter (beglat   = 1)
   parameter (endlat   = plat)
   parameter (numlats  = plat)
   parameter (iam      = 0)
   parameter (masterproc = .true.)
#endif
end module pmgrid

