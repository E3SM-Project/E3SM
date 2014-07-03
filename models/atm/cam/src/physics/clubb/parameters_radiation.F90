!-------------------------------------------------------------------------------
! $Id: parameters_radiation.F90 5623 2012-01-17 17:55:26Z connork@uwm.edu $
module parameters_radiation

! Description:
!   Parameters for radiation schemes

! References:
!   None
!-------------------------------------------------------------------------------

  use clubb_precision, only: &
    dp, & ! double precision
    core_rknd

  implicit none

  character(len=20), public :: & 
    rad_scheme  ! Either BUGSrad, simplified, or simplied_bomex

  real( kind = dp ), dimension(1), public :: &
    sol_const ! Solar constant

  real( kind = core_rknd ), public :: &
    radiation_top ! The top of the atmosphere fed into a radiation scheme.
    !               The computational grid should be extended to reach this
    !               altitude.

  ! Albedo values (alvdr is used in the simplifed schemes as well)
  real( kind = dp ), public :: &
    alvdr, &   !Visible direct surface albedo   [-]
    alndr, &   !Near-IR direct surface albedo   [-]
    alvdf, &   !Visible diffuse surface albedo  [-]
    alndf      !Near-IR diffuse surface albedo  [-]


  ! Long-wave constants (simplified radiation)
  real( kind = core_rknd ), public :: &
    kappa, & ! A constant (Duynkerke eqn. 5)                   [m^2/kg]
    F0,    & ! Coefficient for cloud top heating (see Stevens) [W/m^2] 
    F1       ! Coefficient for cloud base heating (see Stevens)[W/m^2]

  ! Short-wave constants
  real( kind = core_rknd ), public :: &
    eff_drop_radius, & ! Effective droplet radius [m] 
    gc, & ! Asymmetry parameter, "g" in Duynkerke           [-]
    omega ! Single-scattering albedo                        [-] 

  real( kind = dp ), public :: &
    slr     ! Fraction of daylight

  real( kind = core_rknd ), public, dimension(20) :: &
    Fs_values, &           ! List of Fs0 values for simplified radiation
    cos_solar_zen_times, & ! List of cosine of the solar zenith angle times
    cos_solar_zen_values   ! List of cosine of the solar zenith angle values

  logical, public :: &
    l_fix_cos_solar_zen, l_sw_radiation

  logical, public :: &
    l_rad_above_cloud ! Use DYCOMS II RF02 heaviside step function

  integer, public :: &
    nparam

  ! Flag to signal the use of the U.S. Standard Atmosphere Profile, 1976
  logical, public :: l_use_default_std_atmosphere

  private ! Default Scope

! OpenMP directives. The first column of these cannot be indented.
!$omp threadprivate(rad_scheme, sol_const, alvdr, alvdf, alndr, alndf, &
!$omp   kappa, F0, F1, eff_drop_radius, gc, omega, radiation_top, Fs_values, &
!$omp   l_rad_above_cloud, cos_solar_zen_times, cos_solar_zen_values, &
!$omp   l_fix_cos_solar_zen, nparam, &
!$omp   l_sw_radiation, l_use_default_std_atmosphere)

end module parameters_radiation
