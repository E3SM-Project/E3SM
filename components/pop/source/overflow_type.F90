!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 MODULE overflow_type

!BOP
! !MODULE: overflow_type
! !DESCRIPTION:
!  This module contains the overflow data types 
!   (previously located in module overflows).
!
! !REVISION HISTORY:
!  SVN:
!  

! !USES:

   use POP_KindsMod

   use blocks, only: nx_block, ny_block
   use domain_size, only: max_blocks_clinic, nt, km 
   use kinds_mod, only: r4, r8, int_kind, log_kind
   use mpi2s_gshalo, only: Schedule_t, Groups_t

   implicit none
   private
   save


! !PUBLIC DATA MEMBERS:

!-----------------------------------------------------------------------
!     list of nomenclature definitions
!-----------------------------------------------------------------------
! 
! ovf    = overflow
! inf    = inflow (refering to inflow region)
! src    = source (either region or grid box)
! ent    = entrainment (either region or grid box)
! prd    = product (either region or grid box)
! reg    = region (for averaging density and tracers over region)
! adj    = adjacent (for averaging density and tracers over adjacent boxes)
! num    = number (usually refers to actual number used based on input)
! no.    = number (usually refers to actual number used based on input)
! locs   = locations (i.e. grid boxes)
! orient = orientation (1,2,3 or 4; refers to grid box sidewall)
! params = parameters
! ssb    = shelf-slope break- shelf/slope transition to abyssal depth
! 
!-----------------------------------------------------------------------
!     define overflow types and parameters
!-----------------------------------------------------------------------

   logical (log_kind),   public  :: &
      overflows_on,     &         ! true=on, false=off
      overflows_interactive       ! true=interactive ovf

   character (POP_charLength), public  :: &
      overflows_infile,           &! overflow info file
      overflows_diag_outfile,     &! current filename for overflow output diagnostics file
      outfile_tmp                  ! temp for appending to outfile name

   character (POP_charLength), public  :: &
      overflows_restart_type,    &! restart type (ccsm_startup, ccsm_continue, ccsm_hybrid, ccsm_branch)
      overflows_restfile          ! overflow restart file name

   integer (int_kind), parameter, public :: &
      max_ovf      =    10,&  ! max no. ocean overflows
      max_kmt      =   200,&  ! max no. overflow kmt changes
      max_src      =    50,&  ! max no. overflow src locations
      max_ent      =    50,&  ! max no. overflow ent locations
      max_prd_sets =    20,&  ! max no. overflow prd sets
      max_prd      =    50    ! max no. overflow prd locs each set

   integer (int_kind), public :: &
      num_ovf                        ! no. of overflows from ovf info file

   type, public :: ovf_params        ! parameters for each overflow
      real      (r8)    :: & 
        lat               ,&  ! latitude (degrees)
        width             ,&  ! strait width (cm)
        source_thick      ,&  ! source water thickness (cm)
        distnc_str_ssb    ,&  ! distance strait to ssb (cm)
        bottom_slope      ,&  ! bottom slope beyond ssb 
        bottom_drag           ! bottom drag coefficient
   end type ovf_params

   type, public :: ovf_kmtbox        ! overflow grid-box for kmt changes
      integer   (int_kind)  :: &
        i                 ,&  ! x index
        j                 ,&  ! y index
        korg              ,&  ! original kmt value
        knew                  ! new kmt value
   end type ovf_kmtbox


   type, public :: ovf_region        ! overflow regional boundaries
      integer   (int_kind)  :: & 
        imin              ,&  ! x index min
        imax              ,&  ! x index max
        jmin              ,&  ! y index min
        jmax              ,&  ! y index ma
        kmin              ,&  ! z index min
        kmax                  ! z index max
   end type ovf_region

   type, public :: ovf_mask_reg      ! overflow regional mask
      real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
        inf               ,&  ! inflow region mask
        src               ,&  ! source region mask
        ent                   ! entrainment region mask
   end type ovf_mask_reg

   type, public :: ovf_mask_adj      ! overflow adjacent mask
      real (r8) :: src(nx_block,ny_block,max_blocks_clinic), &           ! src adj mask
                   ent(nx_block,ny_block,max_blocks_clinic)              ! ent adj mask
      real (r8) :: prd(nx_block,ny_block,max_blocks_clinic,max_prd_sets) ! prd adj mask(s)
   end type ovf_mask_adj

   type, public :: ovf_mask_reg_wght ! overflow regional mask weight
      real (r8) :: &
        inf               ,&  ! inflow region mask weight
        src               ,&  ! source region mask weight
        ent                   ! entrainment region mask weight
   end type ovf_mask_reg_wght

   type, public :: ovf_mask_adj_wght ! overflow adjacent mask weight
      real (r8) :: &
        src               ,&  ! source adj mask weight
        ent               ,&  ! entrainment adj mask weight
        prd(max_prd_sets)     ! product adj mask weight(s)
   end type ovf_mask_adj_wght

   type, public :: ovf_trcr_reg      ! overflow regional tracers
      real (r8), dimension(nt) :: &
        inf               ,&  ! inflow region tracers
        src               ,&  ! source region tracers
        ent               ,&  ! entrainment region tracers
        prd                   ! product region tracers
   end type ovf_trcr_reg

   type, public :: ovf_trcr_adj      ! overflow adjacent tracers
      real (r8), dimension(nt) :: &
        src               ,&  ! source adj tracers
        ent               ,&  ! entrainment adj tracers
        prd                   ! product adj tracers
   end type ovf_trcr_adj

   type, public :: ovf_rho_reg       ! overflow regional density
      real (r8) :: &
        inf               ,&  ! inflow region density
        src               ,&  ! source region density
        ent                   ! entrainment region density
   end type ovf_rho_reg

   type, public :: ovf_rho_adj       ! overflow adj density
      real (r8) :: &
        prd(max_prd_sets)     ! product region density(s)
   end type ovf_rho_adj

   type, public :: ovf_gridbox    ! overflow grid-box info
      integer   (int_kind)  :: &
        i                 ,&  ! x index for t grid
        j                 ,&  ! y index for t grid
        k                 ,&  ! z index for t grid
        orient            ,&  ! sidewall orientation of t grid box
        i_adv             ,&  ! x index for t grid advection
        j_adv             ,&  ! y index for t grid advection
        i_u               ,&  ! x index for u grid
        j_u               ,&  ! y index for u grid
        task_u                ! task number for (i_u,j_u)
      real (r8) :: &
        Utlda(km)         ,&  ! UVEL "tilda" at (n+1) column speed on u grid
        Vtlda(km)         ,&  ! VVEL "tilda" at (n+1) column speed on u grid
        Uovf_nm1          ,&  ! U at (n-1) speed on u grid
        Uovf_n            ,&  ! U at (n)   speed on u grid
        Uovf              ,&  ! U at (n+1) speed on u grid
        Wovf                  ! W at (n+1) vert speed on t grid
   end type ovf_gridbox

!-------------------------------------------------------------------------------
!     type overflow that follows contains all diagnostic and prognostic
!     data for each overflow; complete list for all overflows is contained 
!     in the array ovf. each overflow is specified by regions, grid locations,
!     and adjacent locations. inf (inflow) and src (source) regions are
!     geographically specified volumes from which density differences
!     determine source transport Ms. this transport is assumed to flow into
!     sidewall locations (possibly modified from original topography by any
!     kmt changes) given by src locations, transporting mean tracers from
!     adjacent grid boxes to the sidewall specified by adjacent boundaries.
!     this transport moves unimpeded to the ent (entrainment) locations, where
!     an entrainment region density along with source density (adjusted for
!     depth changes) determines mixing. entrainment tracers from means along
!     adjacent entrainment grid boxes are mixed with source tracers resulting 
!     in a total transport Mp and mixed product tracers. this product is then 
!     injected from a product sidewall for which the product density is neutral 
!     with the adjacent mean density. each product set is a group of points, 
!     and the collection of sets represents a product path of increasing depth.
!
!     the reader is to be commended for taking in this tedious explanation.
!     it is unfortunately necessary to explain the complexity of the overflow
!     parameterization.
!-------------------------------------------------------------------------------

   type, public :: overflow_t          ! individual overflow info  
  ! logicals and name
      logical   (log_kind)          :: interactive           ! T=ovf active with ocn
      character (32)                :: name                  ! name of ovf
  ! parameters
      type      (ovf_params)        :: ovf_params            ! ovf specific params
  ! kmt mods
   integer   (int_kind)          :: num_kmt               ! no. of kmt changes
      type      (ovf_kmtbox)        :: loc_kmt(max_kmt)      ! kmt locs
  ! source locations
      integer   (int_kind)          :: num_src               ! no. of src locs
      type      (ovf_gridbox)       :: loc_src(max_src)      ! src locs
  ! entrainment locations
      integer   (int_kind)          :: num_ent               ! no. of ent locs
      type      (ovf_gridbox)       :: loc_ent(max_ent)      ! ent locs
  ! product sets (various injection locations) and point for each
      integer   (int_kind)          :: num_prd_sets          ! no. prd sets of pnts
      integer   (int_kind)          :: num_prd(max_prd_sets) ! no. prd locs each set
      type      (ovf_gridbox)       :: loc_prd(max_prd_sets,max_prd) ! prd locs
  ! region locations, masks, tracer and density means for inf, src and ent
      type      (ovf_region)        :: reg_inf               ! inf reg boundaries
      type      (ovf_region)        :: reg_src               ! src reg boundaries
      type      (ovf_region)        :: reg_ent               ! ent reg boundaries
      type      (ovf_mask_reg)      :: mask_reg              ! regional inf, src, ent masks
      type      (ovf_mask_reg_wght) :: wght_reg              ! regional mask weights
      type      (ovf_trcr_reg)      :: trcr_reg              ! regional tracers
      type      (ovf_rho_reg)       :: rho_reg               ! regional densities
  ! adjacent locations, masks, tracer and density means for src, ent and prd
      type      (ovf_region)        :: adj_src               ! src adj boundaries
      type      (ovf_region)        :: adj_ent               ! ent adj boundaries
      type      (ovf_region)        :: adj_prd(max_prd_sets) ! prd adj boundaries
      type      (ovf_mask_adj)      :: mask_adj              ! adj mask
      type      (ovf_mask_adj_wght) :: wght_adj              ! adj mask weights
      type      (ovf_trcr_adj)      :: trcr_adj              ! adjacent tracers
      type      (ovf_rho_adj)       :: rho_adj               ! regional densities
  ! overflow transports and state
      real      (r8)                :: Ms                    ! src mass flux (Sv)
      real      (r8)                :: Ms_n                  ! src mass flux (Sv) at n
      real      (r8)                :: Ms_nm1                ! src mass flux (Sv) at n-1
      real      (r8)                :: Me                    ! ent mass flux (Sv)
      real      (r8)                :: Me_n                  ! ent mass flux (Sv) at n
      real      (r8)                :: Me_nm1                ! ent mass flux (Sv) at n-1
      real      (r8)                :: phi                   ! ent parameter (Me/Mp)
      real      (r8)                :: Mp                    ! prd mass flux (Sv)
      real      (r8)                :: Mp_n                  ! prd mass flux (Sv) at n
      real      (r8)                :: Mp_nm1                ! prd mass flux (Sv) at n-1
      real      (r8)                :: Tp                    ! prd temperature (C)
      real      (r8)                :: Sp                    ! prd salinity (ppt)
      integer   (int_kind)          :: prd_set_n             ! prd set index previous time step
      integer   (int_kind)          :: prd_set               ! prd set index
   end type overflow_t

   type (overflow_t), dimension(max_ovf), public, target :: ovf 
   ! contains all overflow info

   ! group information for the overflows
   type (Groups_t), public:: ovf_groups
 



!***********************************************************************

 end module overflow_type

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||










