#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!
! This module should 'use' only module 'kinds'
!
module physical_constants
  ! ------------------------------
  use kinds, only : real_kind, longdouble_kind
#ifdef CAM

  ! set physical constants from CAM variables
  use physconst, only:    &
    pi,                   &  ! _EXTERNAL
    g => gravit,          &
    rearth,               &
    omega,                &
    Rgas => rair,         &
    p0 => pstd,           &
    MWDAIR => mwdry,      &
    Rwater_vapor => rh2o, &
    Cpwater_vapor => cpwv,&
    kappa => cappa,       &
    Rd_on_Rv => epsilo,   &
    Cpd_on_Cpv,           &
    rrearth => ra
  use shr_const_mod, only: shr_const_cpdair
#else
  ! numerical values of the constants come from the Physical Constants Dictionary (share/util/pcd_const.F90)
  use pcd, only:                                                       &
    pcd_pi    => pi,                                                   &
    pcd_g     => standard_acceleration_of_gravity,                     &
    pcd_rearth=> mean_radius,                                          &
    pcd_omega => angular_velocity,                                     &
    pcd_rd    => dry_air_specific_gas_constant_reference,              &
    pcd_rv    => water_vapor_specific_gas_constant_reference,          &
    pcd_cp    => dry_air_specific_heat_capacity_reference,             &
    pcd_cpv   => water_vapor_specific_heat_capacity_reference,         &
    pcd_mwdry => dry_air_molar_mass_reference
#endif
  ! -----------------------------
  implicit none

  private

#ifdef CAM

  ! physical constants inherited from CAM
  real (kind=real_kind), public, parameter :: DD_PI = pi
  real (kind=longdouble_kind), public, parameter :: QQ_PI = 3.141592653589793238462643383279_longdouble_kind
  public                                   :: rearth                    ! m
  public                                   :: g                         ! m s^-2
  public                                   :: omega                     ! s^-1
  public                                   :: Rgas
  real (kind=real_kind), public, parameter :: Cp = shr_const_cpdair     ! cpair from the "physconst" module in CAM is not a constant
  public                                   :: p0                        ! Pa
  public                                   :: MWDAIR
  public                                   :: Rwater_vapor
  public                                   :: Cpwater_vapor
  public                                   :: kappa
  public                                   :: Rd_on_Rv
  public                                   :: Cpd_on_Cpv
  public                                   :: rrearth                   ! 1/m
#else

  ! physical constants used in HOMME stand-alone simulations
  real (kind=real_kind), public, parameter :: DD_PI        = real(pcd_pi, real_kind)
  real (kind=longdouble_kind), public, parameter :: QQ_PI  = 3.141592653589793238462643383279_longdouble_kind
  real (kind=real_kind), public, parameter :: rearth0      = real(pcd_rearth, real_kind) ! m
  real (kind=real_kind), public            :: rearth       = rearth0    ! m
  real (kind=real_kind), public, parameter :: g            = real(pcd_g, real_kind) ! m s^-2
  real (kind=real_kind), public            :: ginv         = 1.0_real_kind/g
  real (kind=real_kind), public, parameter :: omega0       = real(pcd_omega, real_kind) ! s^-1
  real (kind=real_kind), public            :: omega        = omega0
  real (kind=real_kind), public, parameter :: Rgas         = real(pcd_rd, real_kind)
  real (kind=real_kind), public, parameter :: Cp           = real(pcd_cp, real_kind)
  real (kind=real_kind), public, parameter :: p0           = 100000.0D0 ! [Pa] (no PCD equivalent: standard_atmosphere is 101325 Pa)
  real (kind=real_kind), public, parameter :: MWDAIR       = real(1000*pcd_mwdry, real_kind) ! g mol^-1
  real (kind=real_kind), public, parameter :: Rwater_vapor = real(pcd_rv, real_kind)
  real (kind=real_kind), public, parameter :: Cpwater_vapor= real(pcd_cpv, real_kind)
  real (kind=real_kind), public, parameter :: kappa        = Rgas/Cp
  real (kind=real_kind), public, parameter :: Rd_on_Rv     = Rgas/Rwater_vapor
  real (kind=real_kind), public, parameter :: Cpd_on_Cpv   = Cp/Cpwater_vapor
  real (kind=real_kind), public, parameter :: rrearth0     = 1.0_real_kind/rearth0
  real (kind=real_kind), public            :: rrearth      = rrearth0

#endif

real (kind=real_kind), public :: scale_factor, scale_factor_inv, domain_size, laplacian_rigid_factor
real (kind=real_kind), public :: Lx = 0.0
real (kind=real_kind), public :: Ly = 0.0 
real (kind=real_kind), public :: Sx = 0.0
real (kind=real_kind), public :: Sy = 0.0
real (kind=real_kind), public :: dx, dy, dx_ref, dy_ref

  real (kind=real_kind), public, parameter :: TREF = 288

end module physical_constants
