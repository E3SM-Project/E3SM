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

    ! prognostic variables for preqx solver

    ! prognostics must match those in prim_restart_mod.F90
    ! vertically-lagrangian code advects dp3d instead of ps_v
    ! tracers Q, Qdp always use 2 level time scheme

    real (kind=real_kind) :: v   (np,np,2,nlev,timelevels)            ! horizontal velocity 
    real (kind=real_kind) :: w   (np,np,nlev,timelevels)              ! vertical velocity                  
    real (kind=real_kind) :: theta(np,np,nlev,timelevels)             ! potential temperature                       
    real (kind=real_kind) :: phi(np,np,nlev,timelevels)               ! geopotential 
    real (kind=real_kind) :: dp3d(np,np,nlev,timelevels)              ! delta p on levels                  
    real (kind=real_kind) :: ps_v(np,np,timelevels)                   ! surface pressure                   
    real (kind=real_kind) :: phis(np,np)                              ! surface geopotential (prescribed)  
    real (kind=real_kind) :: Q   (np,np,nlev,qsize_d)                 ! Tracer concentration               
    real (kind=real_kind) :: Qdp (np,np,nlev,qsize_d,2)               ! Tracer mass                        

  end type elem_state_t

  !___________________________________________________________________
  type, public :: derived_state_t

    ! diagnostic variables for preqx solver

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

    real (kind=real_kind) :: FQ(np,np,nlev,qsize_d)                ! tracer forcing
    real (kind=real_kind) :: FM(np,np,2,nlev)                      ! momentum forcing
    real (kind=real_kind) :: FT(np,np,nlev)                        ! temperature forcing
    real (kind=real_kind) :: pecnd(np,np,nlev)                        ! pressure perturbation from condensate
    real (kind=real_kind) :: FQps(np,np)                   ! forcing of FQ on ps_v

  end type derived_state_t
  

  !___________________________________________________________________
  type, public :: elem_accum_t

#ifdef ENERGY_DIAGNOSTICS

    ! Energy equation:
    ! KE_t  = KE1 + KE2 + KEhoriz2 +  KEvert1 + KEvert2 + P1 + T1 + T2 + D1 + Err
    ! IE_t  = S1 + S2 + IEvert1 + d(p dphi/dt)/deta + ptop dphitop/dt + D2 
    ! PE_t  = PEhoriz1 + PEvert1 + P1
    !
    ! d(p dphi/dt)/deta + ptop dphitop/dt should vertically integrate to zero
    !
    ! KEvert*  =  KE net vertical advection    (should be zero)
    ! KEhoriz* =  KE net horizontal advection  (* = 2 might be zero due to horizontal gradient of vertical velocity)
    ! IEvert*  =  IE net vertical advection    (should be zero)
    ! IEhoriz* =  IE net horizontal advection  (should be zero)
    ! PEhoriz* =  PE net horizontal advection  (should be zero)
    !
    ! With leapfrog, energy equations are all exact except KE
    ! (has an Err term that goes to zero as dt**2)
    !
    ! non-transfer terms
    ! KE1 = -0.5*(dpi/deta)*u*grad(u^2)
    ! KE2 = -0.5*u^2 * div( u dpi/deta)
    ! KEhoriz2 = dpi/deta w grad(w)^2 u -0.5 w^2 div(dpi/deta u)
    ! PEhoriz1 = -phi div(u dpi/deta) - grad(phi)^T u dpi/deta
    !
    ! KEvert1  = - etadot u du/deta dpi/deta - 0.5*u^2 d(etadot dpi/deta )/deta
    ! KEvert2  = -etadot w dw/deta dpi/deta - 0.5 w^2 d(etadot dpi/deta)/deta
    ! IEvert1  = -p^kappa d(theta etadot)/deta - theta etadot d p^kappa / deta
    ! PEvert1  = -phi d(edtadot dpi/deta)deta -etadot dpi/deta dphi/deta
    !
    ! Transfer terms:
    ! T1 = -< theta grad_exner,u >             (KE<->IE)_1: T1 + S1 = 0
    ! T2 = gw dp/ds - dp/ds < u,grad(phi)>     (KE<->IE)_2: T2 + S2 = 0
    ! S1 = - exner div(theta u)
    ! S2 = -T2 (the terms are exactly opposite without integration by parts)
    ! P1 = -gw dp/deta 
    ! P2 =  gw dp/deta

    real (kind=real_kind) :: KEvert1(np,np)
    real (kind=real_kind) :: KEvert2(np,np)
    real (kind=real_kind) :: IEvert1(np,np)
    real (kind=real_kind) :: PEvert1(np,np)
    real (kind=real_kind) :: PEvert2(np,np)

    real (kind=real_kind) :: KEhoriz1(np,np)
    real (kind=real_kind) :: KEhoriz2(np,np)
    real (kind=real_kind) :: KE1(np,np)
    real (kind=real_kind) :: KE2(np,np)
    real (kind=real_kind) :: PEhoriz1(np,np)

    real (kind=real_kind) :: T1(np,np)
    real (kind=real_kind) :: T2(np,np)
    real (kind=real_kind) :: S1(np,np)
    real (kind=real_kind) :: S2(np,np)
    real (kind=real_kind) :: P1(np,np)
    real (kind=real_kind) :: P2(np,np)

    ! the KE conversion term and diffusion term
    real (kind=real_kind) :: DIFF(np,np,2,nlev)                       ! net hypervis term
    real (kind=real_kind) :: DIFFTHETA(np,np,nlev)                    ! net hypervis term
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
    real (kind=real_kind) :: Qvar(np,np,qsize_d,4)                    ! Q variance at half time levels
    real (kind=real_kind) :: Qmass(np,np,qsize_d,4)                   ! Q mass at half time levels
    real (kind=real_kind) :: Q1mass(np,np,qsize_d)                    ! Q mass at full time levels

  end type elem_accum_t



contains
end module 
