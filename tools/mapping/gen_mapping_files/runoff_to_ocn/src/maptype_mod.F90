MODULE maptype_mod

   use kind_mod

   implicit none

   real(r8), parameter :: pi =  3.14159265358979323846_r8
   real(r8), parameter :: rEarth     =  6.37122e+6         ! radius of earth (m)
   real(r8), parameter :: DEGtoRAD = pi/180.0_r8           ! degrees to radians
   real(r8), parameter :: RADtoDEG = 180.0_r8/pi           ! radians to degrees

   integer,  parameter :: strLen = 240
   integer,  parameter :: nibx = 360          ! size of sort bin boxes
   integer,  parameter :: njbx = nibx/2
   integer,  parameter :: vcells_req = 20     ! number of valid cells required per bin


   !----------------------------------------------------------------------------
   ! sparse matrix data type
   !----------------------------------------------------------------------------
   TYPE sMatrix

     !--- global text attributes ---
     character(strLen) :: title
     character(strLen) :: normal
     character(strLen) :: method
     character(strLen) :: history
     character(strLen) :: convention
     character(strLen) :: domain_a
     character(strLen) :: domain_b

     !--- domain a ---
     integer         ::    n_a      ! number of non-zero matrix elements
     integer         ::   ni_a      ! number of 2d array i indicies
     integer         ::   nj_a      ! number of 2d array j indicies
     integer         ::   nv_a      ! number of vertices per cell on a grid
     real(r8)   ,pointer ::   xc_a(:)   ! x-coords of centers   ~ deg east
     real(r8)   ,pointer ::   yc_a(:)   ! y-coords of centers   ~ deg north
     real(r8)   ,pointer ::   xv_a(:,:) ! x-coords of verticies ~ deg east, (nv_a,n)
     real(r8)   ,pointer ::   yv_a(:,:) ! y-coords of verticies ~ deg north (nv_a,n)
     integer,pointer :: mask_a(:)   ! mask: 0 <=> out-of-domain (invalid data)
     real(r8)   ,pointer :: area_a(:)   ! area of grid cell ~ radians squared
     integer         :: dims_a(2)       ! hardwire to 2 for now

     !--- domain b ---
     integer         ::    n_b      ! number of non-zero matrix elements
     integer         ::   ni_b      ! number of 2d array i indicies
     integer         ::   nj_b      ! number of 2d array j indicies
     integer         ::   nv_b      ! number of vertices per cell on b grid
     real(r8)   ,pointer ::   xc_b(:)   ! x-coords of centers   ~ deg east
     real(r8)   ,pointer ::   yc_b(:)   ! y-coords of centers   ~ deg north
     real(r8)   ,pointer ::   xv_b(:,:) ! x-coords of verticies ~ deg east, (nv_b,n)
     real(r8)   ,pointer ::   yv_b(:,:) ! y-coords of verticies ~ deg north (nv_b,n)
     integer,pointer :: mask_b(:)   ! mask: 0 <=> out-of-domain (invalid data)
     real(r8)   ,pointer :: area_b(:)   ! area of grid cell ~ radians squared
     integer         :: dims_b(2)       ! hardwire to 2 for now

     !--- fraction of cell mapped to domain b or from domain a ---
     real(r8)   ,pointer :: frac_a(:)   ! area of grid cell ~ radians squared
     real(r8)   ,pointer :: frac_b(:)   ! area of grid cell ~ radians squared

     !--- map: a->b ---
     integer         :: n_s         ! number of non-zero matrix elements
     real(kind=r8)   ,pointer :: s  (:)      ! the non-zero matrix elements
     integer,pointer :: row(:)      ! matrix row corresponding to each element
     integer,pointer :: col(:)      ! matrix col corresponding to each element

     !--- used for OMP/threading in mat-mult ---
     integer,pointer :: sn1(:)      ! # links in a given row
     integer,pointer :: sn2(:)      ! # links previous to a given row

     !--- required for computing NN maps ---
     integer :: imin_b(nibx,njbx)   ! xc_b least index for lat 0:360
     integer :: imax_b(nibx,njbx)   ! xc_b max index for lat 0:360
     integer :: jmin_b(nibx,njbx)   ! yc_b least index for lat 0:90
     integer :: jmax_b(nibx,njbx)   ! yc_b max index for lat 0:90

   END TYPE sMatrix

!===============================================================================
!===============================================================================
END MODULE maptype_mod
!===============================================================================
