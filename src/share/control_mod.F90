#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

! This module contains constants and namelist variables used through out the model
! to avoid circular dependancies please do not 'use' any further modules here.
!
module control_mod
  use kinds, only : real_kind

  integer, public, parameter :: MAX_STRING_LEN=80
  integer, public, parameter :: MAX_FILE_LEN=240
  character(len=MAX_STRING_LEN)    , public :: integration    ! time integration (explicit or semi_implicit)
  integer, public, parameter :: TRACERADV_UGRADQ=0            !  u grad(Q) formulation
  integer, public, parameter :: TRACERADV_TOTAL_DIVERGENCE=1   ! div(u dp/dn Q ) formulation                                                    
  integer, public  :: tracer_advection_formulation  = TRACERADV_TOTAL_DIVERGENCE

!OG
!kmass points to a level with density so that not (\rho q) is limited
!but (\rho q)/\pho= q is limited. really, would only be used in swirl test case for now
!because other test cases have to be reviewed in order to add a const field
  integer, public  :: kmass  = -1


  integer, public  :: tstep_type= 0                           ! 0 = leapfrog
                                                              ! 1 = RK (foward-in-time)
  integer, public  :: rk_stage_user  = 0                      ! number of RK stages to use  
  integer, public  :: ftype = 0                                ! Forcing Type
                                                               ! ftype = 0  HOMME ApplyColumn() type forcing process split
                                                               ! ftype = -1   ignore forcing  (used for testing energy balance)
  integer, public  :: use_cpstar=0                             ! use cp or cp* in T equation                               
  integer, public  :: energy_fixer = 0    ! options appropriate for leapfrog:
                                          !   0:  no fixer, compute energy staggered in time
                                          !   1:  Energy with cp_star, staggered (for debug only)
                                          !   2:  Energy with cp, staggered
                                          ! options appropreate for forward-in-time
                                          ! -1,-2:  no fixer, compute energy non-staggered in time
                                          !   3:  Energy with cp_star, non-staggered (for debug only)
                                          !   4:  Energy with cp, non-staggered


                                          ! 1 = CAM style fixer  Tnew = T + beta
                                              
  integer, public :: qsplit = 1           ! ratio of dynamics tsteps to tracer tsteps        
  integer, public :: LFTfreq=0            ! leapfrog-trapazoidal frequency
                                          ! interspace a lf-trapazoidal step every LFTfreq leapfrogs    
                                          ! 0 = disabled

! compute_mean_flux:  -1  default chosen in namelist_mod.F90
!                      0  use instantaneous velocities for tracer advection
!                      1  use mean flux from dynamics for tracer advection
  integer, public :: compute_mean_flux=-1

! vert_remap_q_alg:    0  default value, Zerroukat monotonic splines
!                      1  PPM vertical remap with mirroring at the boundaries (solid wall bc's, high-order throughout)
!                      2  PPM vertical remap without mirroring at the boundaries (no bc's enforced, first-order at two cells bordering top and bottom boundaries)
  integer, public :: vert_remap_q_alg = 0

!tolerance to define smth small, was introduced for lim 8 in 2d and 3d
  real (kind=real_kind), public, parameter :: tol_limiter=1e-13

  integer              , public :: limiter_option = 0
  character(len=8)     , public :: filter_type
  character(len=8)     , public :: transfer_type
  integer              , public :: filter_freq
  integer              , public :: filter_freq_advection
  integer              , public :: filter_counter
  real (kind=real_kind), public :: filter_mu
  real (kind=real_kind), public :: filter_mu_advection
  character(len=MAX_STRING_LEN)    , public :: precon_method  ! if semi_implicit, type of preconditioner:
                                                  ! choices block_jacobi or identity

  integer              , public :: partmethod     ! partition methods
  character(len=MAX_STRING_LEN)    , public :: topology       ! options: "cube" is supported
  character(len=MAX_STRING_LEN)    , public :: test_case      ! options: if cube: "swtc1","swtc2",or "swtc6"  
  integer              , public :: sub_case                   ! generic test case param 
  integer              , public :: tasknum
  integer              , public :: remapfreq      ! remap frequency of synopsis of system state (steps)
  character(len=MAX_STRING_LEN) :: remap_type     ! selected remapping option
  integer              , public :: statefreq      ! output frequency of synopsis of system state (steps)
  integer              , public :: accumfreq      ! frequency in steps of field accumulation
  integer              , public :: accumstart     ! model day to start accumulation
  integer              , public :: accumstop      ! model day to stop  accumulation
  integer              , public :: restartfreq
  integer              , public :: runtype 
  integer              , public :: timerdetail 
  integer              , public :: numnodes 
  integer              , public :: multilevel
  logical              , public :: uselapi
  character(len=MAX_STRING_LEN)    , public :: restartfile 
  character(len=MAX_STRING_LEN)    , public :: restartdir

  character(len=MAX_STRING_LEN)    , public :: columnpackage
  character(len=MAX_STRING_LEN)    , public :: moisture
  
  integer              , public :: maxits         ! max iterations of solver
  real (kind=real_kind), public :: tol            ! solver tolerance (convergence criteria)
  integer              , public :: debug_level    ! debug level of CG solver


  ! Boyd Vandeven filter Transfer fn parameters

  real (kind=real_kind), public :: p_bv
  real (kind=real_kind), public :: s_bv

  ! Fischer-Mullen filter Transfer fn parameters

  real (kind=real_kind), public :: wght_fm
  integer              , public :: kcut_fm

  character(len=MAX_STRING_LEN)    ,public  :: vfile_int=""  ! vertical formulation (ecmwf,ccm1)
  character(len=MAX_STRING_LEN)    ,public  :: vfile_mid=""  ! vertical grid spacing (equal,unequal)
  character(len=MAX_STRING_LEN)    ,public  :: vform = ""    ! vertical coordinate system (sigma,hybrid)

  integer              , public :: while_iter
  integer              , public :: fine_ne = -1              ! set for refined exodus meshes (variable viscosity)
  real (kind=real_kind), public :: max_hypervis_courant = 1d99 ! upper bound for Courant number (only used for variable viscosity, recommend 1.9 in namelist)
  real (kind=real_kind), public :: nu      = 7.0D5           ! viscosity (momentum equ)
  real (kind=real_kind), public :: nu_div  = -1              ! viscsoity (momentum equ, div component)
  real (kind=real_kind), public :: nu_s    = -1              ! default = nu   T equ. viscosity
  real (kind=real_kind), public :: nu_q    = -1              ! default = nu   tracer viscosity
  real (kind=real_kind), public :: nu_p    = 0.0D5           ! default = 0    ps equ. viscosity
  real (kind=real_kind), public :: nu_top  = 0.0D5           ! top-of-the-model viscosity
  integer, public :: hypervis_subcycle=1    ! number of subcycles for hyper viscsosity timestep
  integer, public :: hypervis_subcycle_q=1  ! number of subcycles for hyper viscsosity timestep on TRACERS
  integer, public :: hypervis_order=0     ! laplace**hypervis_order.  0=not used  1=regular viscosity, 2=grad**4
  integer, public :: psurf_vis = 0        ! 0 = use laplace on eta surfaces
  real (kind=real_kind), public :: hypervis_power=0     ! if not 0, use variable hyperviscosity based on element area          
                                          ! 1 = use (approx.) laplace on p surfaces

  ! hyperviscosity parameters used for smoothing topography
  integer, public :: smooth_phis_numcycle = 0   ! 0 = disable
  integer, public :: smooth_sgh_numcycle = 0   ! 0 = disabled
  real (kind=real_kind), public :: smooth_phis_nudt = 0


  integer, public :: prescribed_wind=0    ! fix the velocities?

  real (kind=real_kind), public :: initial_total_mass = 0    ! initial perturbation in JW test case
  real (kind=real_kind), public :: u_perturb   = 0         ! initial perturbation in JW test case

  integer, public, parameter :: west  = 1
  integer, public, parameter :: east  = 2
  integer, public, parameter :: south = 3
  integer, public, parameter :: north = 4

  integer, public, parameter :: swest = 5
  integer, public, parameter :: seast = 6
  integer, public, parameter :: nwest = 7
  integer, public, parameter :: neast = 8
  
  logical, public            :: test_cfldep=.FALSE.

end module control_mod
