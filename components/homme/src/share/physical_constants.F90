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
    cpair,                &
    p0 => pstd,           &
    MWDAIR => mwdry,      &
    Rwater_vapor => rh2o, &
    Cpwater_vapor => cpwv,&
    kappa => cappa,       &
    Rd_on_Rv => epsilo,   &
    Cpd_on_Cpv,           &
    rrearth => ra
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
  real (kind=real_kind), public, parameter :: Cp           = cpair
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
  real (kind=real_kind), public, parameter :: DD_PI        = 3.141592653589793238462643383279_real_kind
  real (kind=longdouble_kind), public, parameter :: QQ_PI  = 3.141592653589793238462643383279_longdouble_kind
  real (kind=real_kind), public, parameter :: rearth0      = 6.376D6    ! m
  real (kind=real_kind), public            :: rearth       = rearth0    ! m
  real (kind=real_kind), public, parameter :: g            = 9.80616D0  ! m s^-2
  real (kind=real_kind), public            :: ginv         = 1.0_real_kind/g
  real (kind=real_kind), public, parameter :: omega0       = 7.292D-5   ! s^-1
  real (kind=real_kind), public            :: omega        = omega0
  real (kind=real_kind), public, parameter :: Rgas         = 287.04D0        
  real (kind=real_kind), public, parameter :: Cp           = 1005.0D0
  real (kind=real_kind), public, parameter :: p0           = 100000.0D0 ! mbar
  real (kind=real_kind), public, parameter :: MWDAIR       = 28.966D0
  real (kind=real_kind), public, parameter :: Rwater_vapor = 461.50D0
  real (kind=real_kind), public, parameter :: Cpwater_vapor= 1870.0D0
  real (kind=real_kind), public, parameter :: kappa        = Rgas/Cp
  real (kind=real_kind), public, parameter :: Rd_on_Rv     = Rgas/Rwater_vapor	
  real (kind=real_kind), public, parameter :: Cpd_on_Cpv   = Cp/Cpwater_vapor
  real (kind=real_kind), public, parameter :: rrearth0     = 1.0_real_kind/rearth0
  real (kind=real_kind), public            :: rrearth      = rrearth0
#endif

end module physical_constants
