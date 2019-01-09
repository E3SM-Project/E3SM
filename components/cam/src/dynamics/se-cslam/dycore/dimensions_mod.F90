module dimensions_mod
  use shr_kind_mod, only: r8=>shr_kind_r8
#ifdef FVM_TRACERS
  use constituents, only: ntrac_d=>pcnst ! _EXTERNAL
#else
  use constituents, only: qsize_d=>pcnst ! _EXTERNAL
#endif

  implicit none
  private

! set MAX number of tracers.  actual number of tracers is a run time argument  
#ifdef FVM_TRACERS
  integer, parameter         :: qsize_d = 6 ! SE tracers (currently SE supports 6 condensate loading tracers)
#else
  integer, parameter         :: ntrac_d = 0 ! No fvm tracers if CSLAM is off
#endif

  integer,               public :: qsize_condensate_loading = 1 !how many water variables to include in full density
  !
  ! The variables below hold indices of water vapor and condensate loading tracers as well as
  ! associated heat capacities (initialized in dyn_init):
  !
  !   qsize_condensate_loading_idx     = index of water tracers included in condensate loading according to CAM physics
  !   qsize_condensate_loading_idx_gll = index of water tracers included in condensate loading terms for SE tracers
  !
  ! Note that when running without CSLAM then
  !
  !   qsize_condensate_loading_idx_gll = qsize_condensate_loading_idx
  !
  ! but when running with CSLAM then SE tracers are only the water tracers included in the condensate loading
  !
  integer,            allocatable, public :: qsize_condensate_loading_idx(:)    
  integer,            allocatable, public :: qsize_condensate_loading_idx_gll(:)
  real(r8),           allocatable, public :: qsize_condensate_loading_cp(:)
  character(len=16),  allocatable, public :: cnst_name_gll(:)     ! constituent names for SE tracers
  character(len=128), allocatable, public :: cnst_longname_gll(:) ! long name of SE tracers
  !
  !moist cp in energy conversion term
  !
  ! .false.: force dycore to use cpd (cp dry) instead of moist cp
  ! .true. : use moist cp in dycore
  !
  logical           , public :: lcp_moist = .true. 
 
  integer, parameter, public :: np = NP
  integer, parameter, public :: nc = 3       !cslam resolution
  integer           , public :: fv_nphys !physics-grid resolution - the "MAX" is so that the code compiles with NC=0

  integer         :: ntrac = 0 !ntrac is set in dyn_comp
  integer         :: qsize = 0 !qsize is set in dyn_comp
  !
  ! hyperviscosity is applied on approximate pressure levels
  ! Similar to CAM-EUL; see CAM5 scietific documentation (Note TN-486), equation (3.09), page 58.
  ! 
  logical,            public :: hypervis_on_plevs = .true.  
  ! fvm dimensions:
  logical, public :: lprint!for debugging
  integer, parameter, public :: ngpc=3          !number of Gausspoints for the fvm integral approximation   !phl change from 4
  integer, parameter, public :: irecons_tracer=6!=1 is PCoM, =3 is PLM, =6 is PPM for tracer reconstruction
  integer, parameter, public :: nhe=1           !Max. Courant number
  integer, parameter, public :: nhr=2           !halo width needed for reconstruction - phl
  integer, parameter, public :: nht=nhe+nhr     !total halo width where reconstruction is needed (nht<=nc) - phl
  integer, parameter, public :: ns=3!quadratic halo interpolation - recommended setting for nc=3
  !nhc determines width of halo exchanged with neighboring elements
  integer, parameter, public :: nhc = nhr+(nhe-1)+(ns-MOD(ns,2))/2
                                                !(different from halo needed for elements on edges and corners
  integer, parameter, public :: lbc = 1-nhc
  integer, parameter, public :: ubc = nc+nhc


  integer, parameter, public :: kmin_jet=1,kmax_jet=PLEV !min and max level index for the jet

  integer,  public :: nhc_phys 
  integer,  public :: nhe_phys 
  integer,  public :: nhr_phys 
  integer,  public :: ns_phys  

  integer, public :: npdg = 0  ! dg degree for hybrid cg/dg element  0=disabled 

  integer, parameter, public :: npsq = np*np
  integer, parameter, public :: nlev=PLEV
  integer, parameter, public :: nlevp=nlev+1


!  params for a mesh 
!  integer, public, parameter :: max_elements_attached_to_node = 7
!  integer, public, parameter :: s_nv = 2*max_elements_attached_to_node 

  !default for non-refined mesh (note that these are *not* parameters now)
  integer, public  :: max_elements_attached_to_node = 4
  integer, public  :: s_nv = 6
  integer, public  :: max_corner_elem               = 1 !max_elements_attached_to_node-3
  integer, public  :: max_neigh_edges               = 8 !4 + 4*max_corner_elem

  public :: qsize,qsize_d,ntrac_d,ntrac

  integer, public :: ne
  integer, public :: nelem       ! total number of elements
  integer, public :: nelemd      ! number of elements per MPI task
  integer, public :: nelemdmax   ! max number of elements on any MPI task
  integer, public :: nPhysProc                          ! This is the number of physics processors/ per dynamics processor
  integer, public :: nnodes,npart,nmpi_per_node
  integer, public :: GlobalUniqueCols



  public :: set_mesh_dimensions

contains

  subroutine set_mesh_dimensions()

    ! new "params"
    max_elements_attached_to_node = 7  ! variable resolution
    s_nv = 2*max_elements_attached_to_node 

    !recalculate these
    max_corner_elem               = max_elements_attached_to_node-3
    max_neigh_edges               = 4 + 4*max_corner_elem


  end subroutine set_mesh_dimensions


end module dimensions_mod

