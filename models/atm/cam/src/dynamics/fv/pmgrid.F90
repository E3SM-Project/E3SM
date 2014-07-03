
module pmgrid
!BOP
!
! !MODULE: pmgrid --- Initialize grid point resolution parameters
!
! !USES:
   use decompmodule, only : decomptype
   use shr_kind_mod, only : r8=>shr_kind_r8
   implicit none

!
! !DESCRIPTION:  Initialize grid point resolution parameters
! 
! !REVISION HISTORY:
!
!   ??.??.??   ??????     Creation
!   01.03.26   Sawyer     Added ProTeX documentation
!   01.06.27   Mirin      Added 2D decomposition material
!   01.10.16   Sawyer     Added Y-at-Z subdomain decompositions
!   02.05.02   Sawyer     Added XY variable definitions for non-SPMD
!   03.05.07   Sawyer     Removed strip??q3 (?? = 3d, 3k, 4d, 4k)
!   03.07.22   Sawyer     Removed strip3zaty*j?, and t4 (outdated)
!
! !PUBLIC DATA MEMBERS:
   integer, parameter :: plon   = PLON       ! number of longitudes
   integer, parameter :: plev   = PLEV       ! number of vertical levels
   integer, parameter :: plat   = PLAT       ! number of latitudes

   integer, parameter :: plevp  = plev+1     ! plev + 1
   integer, parameter :: numbnd = 0          ! no.of latitudes passed N and S of forecast lat
   integer, parameter :: plnlv = plon*plev   ! Length of multilevel field slice

   integer beglat     ! beg. index for latitudes owned by a given proc
   integer endlat     ! end. index for latitudes owned by a given proc
   integer numlats    ! number of latitudes owned by a given proc

   integer beglev     ! beg. index for levels owned by a given task
   integer endlev     ! end. index for levels owned by a given task
   integer endlevp1   ! end. index + 1 for levels owned by a given task
   integer endlevp    ! equals endlev, except in last subdomain where equals endlevp1

   integer myid_y     ! subdomain index (0-based) in latitude (y)
   integer myid_z     ! subdomain index (0 based) in level (z)
   integer npr_y      ! number of subdomains in y
   integer npr_z      ! number of subdomains in z
   integer :: twod_decomp = 0  ! 1 for multi-2D decomposition with transposes, 0 otherwise
   integer :: mod_transpose = 0  ! method for computing mod_comm transposes
   integer :: mod_geopk = 0  ! method for computing mod_comm geopk
   integer :: mod_gatscat = 0  ! method for computing mod_comm gather/scatters

!  secondary xy decomposition used for remapping

   integer myidxy_x     ! subdomain index (0-based) in longitude (x) (second. decomp.)
   integer myidxy_y     ! subdomain index (0 based) in latitude (y) (second. decomp.)
   integer nprxy_x      ! number of subdomains in x (second. decomp.)
   integer nprxy_y      ! number of subdomains in y (second. decomp.)
   integer beglonxy     ! beg. index for longitudes (second. decomp.)
   integer endlonxy     ! end. index for longitudes (second. decomp.)
   integer beglatxy     ! beg. index for latitudes (second. decomp.)
   integer endlatxy     ! end. index for latitudes (second. decomp.)

   logical :: dyndecomp_set = .false. ! flag indicates dynamics grid has been set
   integer :: spmd_on = 0 ! 1 for Spmd, 0 for non-Spmd

!
#if ( ! defined SPMD )
   parameter (beglat   = 1)
   parameter (endlat   = plat)
   parameter (numlats  = plat)
   parameter (beglev = 1)
   parameter (endlev = plev)
   parameter (endlevp1 = plev+1)
   parameter (endlevp = plev+1)
   parameter (myid_y = 0)
   parameter (myid_z = 0)
   parameter (npr_y = 1)
   parameter (npr_z = 1)
!
! These are needed to pass strict run-time error checking
!
   parameter (myidxy_x=0)     ! subdomain index (0-based) in longitude (x) (second. decomp.)
   parameter (myidxy_y=0)     ! subdomain index (0 based) in latitude (y) (second. decomp.)
   parameter (nprxy_x=1)      ! number of subdomains in x (second. decomp.)
   parameter (nprxy_y=1)      ! number of subdomains in y (second. decomp.)
   parameter (beglonxy=1)     ! beg. index for longitudes (second. decomp.)
   parameter (endlonxy=plon)     ! end. index for longitudes (second. decomp.)
   parameter (beglatxy=1)     ! beg. index for latitudes (second. decomp.)
   parameter (endlatxy=plat)     ! end. index for latitudes (second. decomp.)
#endif

! Staggered grid parameters
! splon and splat may eventually need to become new parameters
! in params.h to define the size of the staggered grid arrays. - gg

   integer, parameter :: splon = plon     ! Number of longitudes on the staggered grid
   integer, parameter :: splat = plat     ! Number of latitudes on the staggered grid

! Note: In reality, the staggered latitude array for Lin-Rood dynamics only
! uses PLAT-1 latitudes, the first one being ignored.  So ideally the line
! above should read:
!   parameter (splat = plat-1)
! to define the staggered latitude winds with the correct dimension.
! However, the assumption that the staggered latitude grid has one extra
! latitude (making it the same dimension as the non-staggered grid) is
! pervasive throughout the Lin-Rood dynamical core, necessitating the
! extra latitude.
!
! A temporary parameter for I/O purposes is defined below. When the Lin-Rood
! dynamical core has been changed to use a staggered-U wind arrays with the
! correct dimension (plat-1), replace all occurances of platm1 with splat,
! set splat to plat-1, and delete the following definition for platm1:

!  integer platm1
!  parameter (platm1 = plat-1)

!EOP
end module pmgrid
