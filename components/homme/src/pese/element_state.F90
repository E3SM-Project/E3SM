#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module element_state

  use kinds,                  only: real_kind, long_kind, int_kind
  use dimensions_mod,         only: np, npsq, nlev, nlevp, qsize_d

  implicit none
  private
  integer, public, parameter :: timelevels = 3

! =========== PRIMITIVE-EQUATION DATA-STRUCTURES =====================

  type, public :: elem_state_t

    ! prognostic variables for pese solver

    ! prognostics must match those in prim_restart_mod.F90
    ! vertically-lagrangian code advects dp3d instead of ps_v
    ! tracers Q, Qdp always use 2 level time scheme

    real (kind=real_kind) :: v   (np,np,2,nlev,timelevels)            ! velocity                           1
    real (kind=real_kind) :: T   (np,np,nlev,timelevels)              ! temperature                        2
    real (kind=real_kind) :: dp3d(np,np,nlev,timelevels)              ! delta p on levels                  8
    real (kind=real_kind) :: ps_v(np,np,timelevels)                   ! surface pressure                   4
    real (kind=real_kind) :: phis(np,np)                              ! surface geopotential (prescribed)  5
    real (kind=real_kind) :: Q   (np,np,nlev,qsize_d)                 ! Tracer concentration               6
    real (kind=real_kind) :: Qdp (np,np,nlev,qsize_d,2)               ! Tracer mass                        7

  end type elem_state_t

  !___________________________________________________________________
  type, public :: derived_state_t

    ! diagnostic variables for pese solver

    real (kind=real_kind) :: omega(np,np,nlev)                        ! dp/dt

    ! storage for subcycling tracers/dynamics

    real (kind=real_kind) :: vn0  (np,np,2,nlev)                      ! velocity for SE tracer advection
    real (kind=real_kind) :: vstar(np,np,2,nlev)                      ! velocity on Lagrangian surfaces
    real (kind=real_kind) :: dpdiss_biharmonic(np,np,nlev)            ! mean dp dissipation tendency, if nu_p>0
    real (kind=real_kind) :: dpdiss_ave(np,np,nlev)                   ! mean dp used to compute psdiss_tens

    ! diagnostics for explicit timestep
    real (kind=real_kind) :: phi(np,np,nlev)                          ! geopotential
    real (kind=real_kind) :: omega_p(np,np,nlev)                      ! vertical tendency (derived)
    real (kind=real_kind) :: eta_dot_dpdn(np,np,nlevp)                ! mean vertical flux from dynamics
    real (kind=real_kind) :: eta_dot_dpdn_prescribed(np,np,nlevp)     ! prescribed wind test cases

    ! semi-implicit diagnostics: computed in explict-component, reused in Helmholtz-component.
    real (kind=real_kind) :: grad_lnps(np,np,2)                       ! gradient of log surface pressure
    real (kind=real_kind) :: zeta(np,np,nlev)                         ! relative vorticity
    real (kind=real_kind) :: div(np,np,nlev,timelevels)               ! divergence

    ! tracer advection fields used for consistency and limiters
    real (kind=real_kind) :: dp(np,np,nlev)                           ! for dp_tracers at physics timestep
    real (kind=real_kind) :: divdp(np,np,nlev)                        ! divergence of dp
    real (kind=real_kind) :: divdp_proj(np,np,nlev)                   ! DSSed divdp

    ! forcing terms for CAM
    real (kind=real_kind) :: FQ(np,np,nlev,qsize_d)                   ! tracer forcing
    real (kind=real_kind) :: FM(np,np,2,nlev)                         ! momentum forcing
    real (kind=real_kind) :: FT(np,np,nlev)                           ! temperature forcing
    real (kind=real_kind) :: pecnd(np,np,nlev)                        ! pressure perturbation from condensate
    real (kind=real_kind) :: FQps(np,np)                              ! forcing of FQ on ps_v

  end type derived_state_t
  

  !___________________________________________________________________
  type, public :: elem_accum_t

#ifdef ENERGY_DIAGNOSTICS

    ! Energy equation:
    ! KE_t  = T1 + T2  + D1   + Err   +  vertical & horizontal advection terms
    ! IE_t  = S1 + D2                 +  vertical & horizontal advection terms
    ! PE_t  = S2
    !
    ! KEvert*  =  KE net vertical advection    (should be zero)
    ! KEhoriz* =  KE net horizonatl advection  (should be zero)
    ! IEvert*  =  IE net vertical advection    (should be zero)
    ! IEhoriz* =  IE net horizonatl advection  (should be zero)
    !
    ! With leapfrog, energy equations are all exact except KE
    ! (has an Err term that goes to zero as dt**2)
    !
    ! Transfer terms:
    ! T1   = -< dp/dn u, RT_v/p grad_p >     KE<->IE:   T1 + T2-T2_s = S1
    ! T2   = -< dp/dn u, grad_phi >          KE<->PE:   T2_s         = S2
    ! T2_s = -< dp/dn u, grad_phis >
    ! S1   = < Cp_star dp/dn , RT omega_p/Cp_star >
    ! S2   = -< div (u dp/dn), phis >

    real (kind=real_kind) :: KEvert1(np,np)                           ! term from continuity equ
    real (kind=real_kind) :: KEvert2(np,np)                           ! term from momentum equ
    real (kind=real_kind) :: IEvert1(np,np)                           ! term from continuity equ
    real (kind=real_kind) :: IEvert2(np,np)                           ! term from T equ
    real (kind=real_kind) :: IEvert1_wet(np,np)                       ! wet term from continuity equ
    real (kind=real_kind) :: IEvert2_wet(np,np)                       ! wet term from T equ

    real (kind=real_kind) :: KEhorz1(np,np)                           ! at time t
    real (kind=real_kind) :: KEhorz2(np,np)                           ! after calling time_advance, these will be at time t-1
    real (kind=real_kind) :: IEhorz1(np,np)
    real (kind=real_kind) :: IEhorz2(np,np)
    real (kind=real_kind) :: IEhorz1_wet(np,np)
    real (kind=real_kind) :: IEhorz2_wet(np,np)

    real (kind=real_kind) :: T1(np,np)
    real (kind=real_kind) :: T2(np,np)
    real (kind=real_kind) :: T2_s(np,np)
    real (kind=real_kind) :: S1(np,np)
    real (kind=real_kind) :: S1_wet(np,np)
    real (kind=real_kind) :: S2(np,np)

    ! the KE conversion term and diffusion term
    real (kind=real_kind) :: DIFF(np,np,2,nlev)                       ! net hypervis term
    real (kind=real_kind) :: DIFFT(np,np,nlev)                        ! net hypervis term
    real (kind=real_kind) :: CONV(np,np,2,nlev)                       ! dpdn u dot CONV = T1 + T2
#endif

    ! the "4" timelevels represents data computed at:
    !  1  t-.5
    !  2  t+.5   after dynamics
    !  3  t+.5   after forcing
    !  4  t+.5   after Robert
    ! after calling TimeLevelUpdate, all times above decrease by 1.0

    real (kind=real_kind) :: KEner(np,np,4)
    real (kind=real_kind) :: PEner(np,np,4)
    real (kind=real_kind) :: IEner(np,np,4)
    real (kind=real_kind) :: IEner_wet(np,np,4)
    real (kind=real_kind) :: Qvar(np,np,qsize_d,4)                    ! Q variance at half time levels
    real (kind=real_kind) :: Qmass(np,np,qsize_d,4)                   ! Q mass at half time levels
    real (kind=real_kind) :: Q1mass(np,np,qsize_d)                    ! Q mass at full time levels

  end type elem_accum_t

contains
end module 
