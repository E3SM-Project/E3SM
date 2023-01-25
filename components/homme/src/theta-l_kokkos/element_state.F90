#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module element_state

  use kinds,                  only: real_kind, long_kind, int_kind
  use dimensions_mod,         only: np, npsq, nlev, nlevp, qsize_d

  implicit none
  private

  public :: allocate_element_arrays
  public :: deallocate_element_arrays
  public :: setup_element_pointers_ie

#ifdef ARKODE
  integer, public, parameter :: timelevels = 50
#else
  integer, public, parameter :: timelevels = 3
#endif
  integer, public, parameter :: diagtimes = 6

  ! maximum number of Newton iterations taken for an IMEX-RK stage per time-step
  integer, public               :: max_itercnt=0
  real (kind=real_kind), public :: max_deltaerr=0
  real (kind=real_kind), public :: max_reserr=0

  ! pressure based TOM sponge layer
  real (kind=real_kind),public :: nu_scale_top(nlev)
  integer, public              :: nlev_tom

  ! flattened arrays for all state, derived, and accum quantities that need to be passed back and forth to CXX

  real (kind=real_kind), allocatable, target, public :: elem_state_v    (:,:,:,:,:,:)           ! horizontal velocity 
  real (kind=real_kind), allocatable, target, public :: elem_state_w_i  (:,:,:,:,:)             ! vertical velocity at interfaces
  real (kind=real_kind), allocatable, target, public :: elem_state_vtheta_dp  (:,:,:,:,:)       ! virtual potential temperature (mass)
  real (kind=real_kind), allocatable, target, public :: elem_state_phinh_i  (:,:,:,:,:)         ! geopotential used by NH model at interfaces
  real (kind=real_kind), allocatable, target, public :: elem_state_dp3d (:,:,:,:,:)             ! delta p on levels                  
  real (kind=real_kind), allocatable, target, public :: elem_state_ps_v (:,:,:,:)               ! surface pressure                   
  real (kind=real_kind), allocatable, target, public :: elem_state_phis (:,:,:)                 ! surface geopotential (prescribed)  
  real (kind=real_kind), allocatable, target, public :: elem_state_Q    (:,:,:,:,:)             ! Tracer concentration               
  real (kind=real_kind), allocatable, target, public :: elem_state_Qdp  (:,:,:,:,:,:)           ! Tracer mass                        

  real (kind=real_kind), allocatable, target, public :: elem_derived_omega_p (:,:,:,:)          ! vertical tendency (derived)
  real (kind=real_kind), allocatable, target, public :: elem_derived_eta_dot_dpdn (:,:,:,:)     ! used with prescribed_wind
  real (kind=real_kind), allocatable, target, public :: elem_derived_vn0 (:,:,:,:,:)            ! used with prescribed_wind

  real (kind=real_kind), allocatable, target, public :: elem_accum_KEner     (:,:,:,:)
  real (kind=real_kind), allocatable, target, public :: elem_accum_PEner     (:,:,:,:)
  real (kind=real_kind), allocatable, target, public :: elem_accum_IEner     (:,:,:,:)
  real (kind=real_kind), allocatable, target, public :: elem_accum_Qvar      (:,:,:,:,:)        ! Q variance at half time levels
  real (kind=real_kind), allocatable, target, public :: elem_accum_Qmass     (:,:,:,:,:)        ! Q mass at half time levels
  real (kind=real_kind), allocatable, target, public :: elem_accum_Q1mass    (:,:,:,:)          ! Q mass at full time levels

  real (kind=real_kind), allocatable, target, public :: elem_derived_FQ(:,:,:,:,:)              ! tracer forcing
  real (kind=real_kind), allocatable, target, public :: elem_derived_FM(:,:,:,:,:)              ! momentum forcing
  real (kind=real_kind), allocatable, target, public :: elem_derived_FT(:,:,:,:)                ! temperature forcing
  real (kind=real_kind), allocatable, target, public :: elem_derived_FVTheta(:,:,:,:)           ! potential temperature forcing
  real (kind=real_kind), allocatable, target, public :: elem_derived_FPHI(:,:,:,:)              ! PHI (NH) forcing
  real (kind=real_kind), allocatable, target, public :: elem_derived_FQps(:,:,:)                ! forcing of FQ on ps_v
  real (kind=real_kind), allocatable, target, public :: elem_derived_lap_p_wk(:,:,:,:)

  ! Reference states (computed from phis)
  real (kind=real_kind), allocatable, target, public :: elem_theta_ref(:,:,:,:)
  real (kind=real_kind), allocatable, target, public :: elem_dp_ref(:,:,:,:)
  real (kind=real_kind), allocatable, target, public :: elem_phi_ref(:,:,:,:)

! =========== PRIMITIVE-EQUATION DATA-STRUCTURES =====================

  type, public :: elem_state_t

    ! prognostic variables for preqx solver

    ! prognostics must match those in prim_restart_mod.F90
    ! vertically-lagrangian code advects dp3d instead of ps_v
    ! tracers Q, Qdp always use 2 level time scheme

    real (kind=real_kind), pointer :: v   (:,:,:,:,:)       ! horizontal velocity 
    real (kind=real_kind), pointer :: w_i (:,:,:,:)         ! vertical velocity at interfaces
    real (kind=real_kind), pointer :: vtheta_dp(:,:,:,:)    ! virtual potential temperature (mass)
    real (kind=real_kind), pointer :: phinh_i(:,:,:,:)      ! geopotential used by NH model at interfaces
    real (kind=real_kind), pointer :: dp3d(:,:,:,:)         ! delta p on levels                  
    real (kind=real_kind), pointer :: ps_v(:,:,:)           ! surface pressure                   
    real (kind=real_kind), pointer :: phis(:,:)             ! surface geopotential (prescribed)  
    real (kind=real_kind), pointer :: Q   (:,:,:,:)         ! Tracer concentration               
    real (kind=real_kind), pointer :: Qdp (:,:,:,:,:)       ! Tracer mass                        
  end type elem_state_t

  !___________________________________________________________________
  type, public :: derived_state_t

    ! diagnostic variables for preqx solver

    ! storage for subcycling tracers/dynamics

    real (kind=real_kind), pointer :: vn0(:,:,:,:)                    ! velocity for SE tracer advection
    real (kind=real_kind) :: vstar(np,np,2,nlev)                      ! velocity on Lagrangian surfaces
    real (kind=real_kind) :: dpdiss_biharmonic(np,np,nlev)            ! mean dp dissipation tendency, if nu_p>0
    real (kind=real_kind) :: dpdiss_ave(np,np,nlev)                   ! mean dp used to compute psdiss_tens

    ! diagnostics
    real (kind=real_kind), pointer :: omega_p(:,:,:)                  ! vertical tendency (derived)
    real (kind=real_kind), pointer :: eta_dot_dpdn(:,:,:)             ! mean vertical flux from dynamics
    real (kind=real_kind) :: eta_dot_dpdn_prescribed(np,np,nlevp)     ! prescribed wind test cases

    ! tracer advection fields used for consistency and limiters
    real (kind=real_kind) :: dp(np,np,nlev)                           ! for dp_tracers at physics timestep
    real (kind=real_kind) :: divdp(np,np,nlev)                        ! divergence of dp
    real (kind=real_kind) :: divdp_proj(np,np,nlev)                   ! DSSed divdp

    real (kind=real_kind), pointer :: FQ(:,:,:,:)                     ! tracer forcing
    real (kind=real_kind), pointer :: FM(:,:,:,:)                     ! momentum forcing
    real (kind=real_kind), pointer :: FT(:,:,:)                       ! temperature forcing
    real (kind=real_kind), pointer :: FVTheta(:,:,:)                  ! potential temperature forcing
    real (kind=real_kind), pointer :: FPHI(:,:,:)                     ! PHI (NH) forcing
    real (kind=real_kind), pointer :: FQps(:,:)                       ! forcing of FQ on ps_v
    real (kind=real_kind), pointer :: lap_p_wk(:,:,:)                 ! 

    real (kind=real_kind) :: gradphis(np,np,2)   ! grad phi at the surface, computed once in model initialization
    real (kind=real_kind), pointer :: dp_ref(:,:,:)         ! ref states based on PHIS
    real (kind=real_kind), pointer :: theta_ref(:,:,:)
    real (kind=real_kind), pointer :: phi_ref(:,:,:)

  end type derived_state_t
  

  !___________________________________________________________________
  type, public :: elem_accum_t

#ifdef ENERGY_DIAGNOSTICS
    ! Energy equation:
    real (kind=real_kind) :: KEu_horiz1(np,np)
    real (kind=real_kind) :: KEu_horiz2(np,np)
    real (kind=real_kind) :: KEu_vert1(np,np)
    real (kind=real_kind) :: KEu_vert2(np,np)
    real (kind=real_kind) :: KEw_horiz1(np,np)  ! nonhydro only
    real (kind=real_kind) :: KEw_horiz2(np,np)  ! nonhydro only
    real (kind=real_kind) :: KEw_horiz3(np,np)  ! nonhydro only
    real (kind=real_kind) :: KEw_vert1(np,np)   ! nonhydro only
    real (kind=real_kind) :: KEw_vert2(np,np)   ! nonhydro only

    real (kind=real_kind) :: IEvert1(np,np)
    real (kind=real_kind) :: IEvert2(np,np)     ! nonhydro only
    real (kind=real_kind) :: PEvert1(np,np)
    real (kind=real_kind) :: PEvert2(np,np)
    real (kind=real_kind) :: PEhoriz1(np,np)
    real (kind=real_kind) :: PEhoriz2(np,np)

    real (kind=real_kind) :: T01(np,np)
    real (kind=real_kind) :: T2(np,np)
    real (kind=real_kind) :: S1(np,np)
    real (kind=real_kind) :: S2(np,np)
    real (kind=real_kind) :: P1(np,np)
    real (kind=real_kind) :: P2(np,np)
    real (kind=real_kind) :: T2_nlevp_term(np,np)

    real (kind=real_kind) :: CONV(np,np,2,nlev)                       ! dpdn u dot CONV = T1 + T2
#endif

    ! the last dimension is "4" (timelevels) represents data computed at:
    !  1  t-.5
    !  2  t+.5   after dynamics
    !  3  t+.5   after forcing
    !  4  t+.5   after Robert
    ! after calling TimeLevelUpdate, all times above decrease by 1.0

    real (kind=real_kind), pointer :: KEner     (:,:,:)
    real (kind=real_kind), pointer :: PEner     (:,:,:)
    real (kind=real_kind), pointer :: IEner     (:,:,:)
    real (kind=real_kind), pointer :: Qvar      (:,:,:,:)             ! Q variance at half time levels
    real (kind=real_kind), pointer :: Qmass     (:,:,:,:)             ! Q mass at half time levels
    real (kind=real_kind), pointer :: Q1mass    (:,:,:)               ! Q mass at full time levels

  end type elem_accum_t

contains

  subroutine allocate_element_arrays (nelemd)
    !
    ! Inputs
    !
    integer, intent(in) :: nelemd

    ! State
    allocate(elem_state_v         (np,np,2,nlev, timelevels,nelemd) )
    allocate(elem_state_w_i       (np,np,  nlevp,timelevels,nelemd) )
    allocate(elem_state_vtheta_dp (np,np,  nlev, timelevels,nelemd) )
    allocate(elem_state_phinh_i   (np,np,  nlevp,timelevels,nelemd) )
    allocate(elem_state_dp3d      (np,np,  nlev, timelevels,nelemd) )
    allocate(elem_state_ps_v      (np,np,        timelevels,nelemd) )
    allocate(elem_state_phis      (np,np,                   nelemd) )
    allocate(elem_state_Q         (np,np,  nlev, qsize_d,   nelemd) )
    allocate(elem_state_Qdp       (np,np,  nlev, qsize_d,2, nelemd) )

    ! Derived
    allocate(elem_derived_omega_p (np,np,nlev,nelemd)       )
    allocate(elem_derived_vn0     (np,np,2,nlev,nelemd)     )
    allocate(elem_derived_eta_dot_dpdn (np,np,nlevp,nelemd) )

    ! Accum
    allocate(elem_accum_kener     (np,np,        diagtimes,nelemd) )
    allocate(elem_accum_pener     (np,np,        diagtimes,nelemd) )
    allocate(elem_accum_iener     (np,np,        diagtimes,nelemd) )
    allocate(elem_accum_qvar      (np,np,qsize_d,diagtimes,nelemd) )
    allocate(elem_accum_qmass     (np,np,qsize_d,diagtimes,nelemd) )
    allocate(elem_accum_Q1mass    (np,np,qsize_d,  nelemd) )

    ! Forcing
    allocate(elem_derived_FM      (np,np,3,nlev,nelemd) )
    allocate(elem_derived_FT      (np,np,nlev,nelemd) )
    allocate(elem_derived_FVTheta (np,np,nlev,nelemd) )
    allocate(elem_derived_FPHI    (np,np,nlevp,nelemd) )
    allocate(elem_derived_FQ      (np,np,nlev,qsize_d,nelemd) )
    allocate(elem_derived_FQps    (np,np,nelemd) )
    allocate(elem_derived_lap_p_wk(np,np,nlevp,nelemd) )

    ! Reference states
    allocate(elem_theta_ref (np,np,nlev,nelemd) )
    allocate(elem_dp_ref    (np,np,nlev,nelemd) )
    allocate(elem_phi_ref   (np,np,nlevp,nelemd) )

  end subroutine allocate_element_arrays

  subroutine deallocate_element_arrays ()
    ! This subroutine is pointless in a normal run, since all
    ! allocatable variables are automatically deallocated upon
    ! program termination. However, it is useful in cxx unit
    ! tests that compare against f90, due to how catch2 testing
    ! mechanism works. There, the same block of code is executed multiple
    ! times, which means the init functions are called multiple
    ! times. To make the 'allocate' calls above succeed, the
    ! variable must not be already allocated. Therefore, at the
    ! end of the unit test code, we clean up the f90 stuff,
    ! so that any potential subsequent call to initialization
    ! subroutines will not generate an error

    ! State
    deallocate(elem_state_v         )
    deallocate(elem_state_w_i       )
    deallocate(elem_state_vtheta_dp )
    deallocate(elem_state_phinh_i   )
    deallocate(elem_state_dp3d      )
    deallocate(elem_state_ps_v      )
    deallocate(elem_state_phis      )
    deallocate(elem_state_Q         )
    deallocate(elem_state_Qdp       )

    ! Derived
    deallocate(elem_derived_omega_p )
    deallocate(elem_derived_vn0     )
    deallocate(elem_derived_eta_dot_dpdn )

    ! Accum
    deallocate(elem_accum_kener     )
    deallocate(elem_accum_pener     )
    deallocate(elem_accum_iener     )
    deallocate(elem_accum_qvar      )
    deallocate(elem_accum_qmass     )
    deallocate(elem_accum_Q1mass    )

    ! Forcing
    deallocate(elem_derived_FM      )
    deallocate(elem_derived_FT      )
    deallocate(elem_derived_FVTheta )
    deallocate(elem_derived_FPHI    )
    deallocate(elem_derived_FQ      )
    deallocate(elem_derived_FQps    )
    deallocate(elem_derived_lap_p_wk)

    ! Reference states
    deallocate(elem_theta_ref       )
    deallocate(elem_dp_ref          )
    deallocate(elem_phi_ref         )
  end subroutine deallocate_element_arrays

  subroutine setup_element_pointers_ie (ie, state, derived, accum)
    !
    ! Inputs
    !
    integer, intent(in) :: ie
    type (elem_state_t),    intent(inout) :: state
    type (derived_state_t), intent(inout) :: derived
    type (elem_accum_t),    intent(inout) :: accum

    ! State
    state%v         => elem_state_v(:,:,:,:,:,ie)
    state%w_i       => elem_state_w_i(:,:,:,:,ie)
    state%vtheta_dp => elem_state_vtheta_dp(:,:,:,:,ie)
    state%phinh_i   => elem_state_phinh_i(:,:,:,:,ie)
    state%dp3d      => elem_state_dp3d(:,:,:,:,ie)
    state%ps_v      => elem_state_ps_v(:,:,:,ie)
    state%Q         => elem_state_Q(:,:,:,:,ie)
    state%Qdp       => elem_state_Qdp(:,:,:,:,:,ie)
    state%phis      => elem_state_phis(:,:,ie)

    ! Derived
    derived%omega_p => elem_derived_omega_p(:,:,:,ie)
    derived%vn0 => elem_derived_vn0(:,:,:,:,ie)
    derived%eta_dot_dpdn => elem_derived_eta_dot_dpdn(:,:,:,ie)

    ! Accum
    accum%KEner     => elem_accum_KEner    (:,:,:,ie)
    accum%PEner     => elem_accum_PEner    (:,:,:,ie)
    accum%IEner     => elem_accum_IEner    (:,:,:,ie)
    accum%Qvar      => elem_accum_Qvar     (:,:,:,:,ie)
    accum%Qmass     => elem_accum_Qmass    (:,:,:,:,ie)
    accum%Q1mass    => elem_accum_Q1mass   (:,:,:,ie)

    ! Forcing
    derived%FM      => elem_derived_FM(:,:,:,:,ie)
    derived%FT      => elem_derived_FT(:,:,:,ie)
    derived%FVTheta => elem_derived_FVTheta(:,:,:,ie)
    derived%FPHI    => elem_derived_FPHI(:,:,:,ie)
    derived%FQ      => elem_derived_FQ(:,:,:,:,ie)
    derived%FQps    => elem_derived_FQps(:,:,ie)
    derived%lap_p_wk=> elem_derived_lap_p_wk(:,:,:,ie)

    ! Reference states
    derived%theta_ref => elem_theta_ref(:,:,:,ie)
    derived%dp_ref    => elem_dp_ref(:,:,:,ie)
    derived%phi_ref   => elem_phi_ref(:,:,:,ie)
  end subroutine setup_element_pointers_ie
end module 
