! $Id: pdf_parameter_module.F90 5623 2012-01-17 17:55:26Z connork@uwm.edu $
module pdf_parameter_module
! Description:
!   This module defines the derived type pdf_parameter.
! References:
!   None
!-------------------------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd

  implicit none

  private ! Default scope

  public :: pdf_parameter

  type pdf_parameter
    real( kind = core_rknd ) :: &
      w1,          & ! Mean of w for 1st normal distribution                 [m/s]
      w2,          & ! Mean of w for 2nd normal distribution                 [m/s]
      varnce_w1,   & ! Variance of w for 1st normal distribution         [m^2/s^2]
      varnce_w2,   & ! Variance of w for 2nd normal distribution         [m^2/s^2]
      rt1,         & ! Mean of r_t for 1st normal distribution             [kg/kg]
      rt2,         & ! Mean of r_t for 2nd normal distribution             [kg/kg]
      varnce_rt1,  & ! Variance of r_t for 1st normal distribution     [kg^2/kg^2]
      varnce_rt2,  & ! Variance of r_t for 2nd normal distribution     [kg^2/kg^2]
      crt1,        & ! Coefficient for s'                                      [-]
      crt2,        & ! Coefficient for s'                                      [-]
      cthl1,       & ! Coefficient for s'                                    [1/K]
      cthl2,       & ! Coefficient for s'                                    [1/K]
      thl1,        & ! Mean of th_l for 1st normal distribution                [K]
      thl2,        & ! Mean of th_l for 2nd normal distribution                [K]
      varnce_thl1, & ! Variance of th_l for 1st normal distribution          [K^2]
      varnce_thl2, & ! Variance of th_l for 2nd normal distribution          [K^2]
      mixt_frac,   & ! Weight of 1st normal distribution (Sk_w dependent)      [-]
      rc1,         & ! Mean of r_c for 1st normal distribution             [kg/kg]
      rc2,         & ! Mean of r_c for 2nd normal distribution             [kg/kg]
      rsl1,        & ! Mean of r_sl for 1st normal distribution            [kg/kg]
      rsl2,        & ! Mean of r_sl for 2nd normal distribution            [kg/kg]
      cloud_frac1, & ! Cloud fraction for 1st normal distribution              [-]
      cloud_frac2, & ! Cloud fraction for 2nd normal distribution              [-]
      s1,          & ! Mean of s for 1st normal distribution               [kg/kg]
      s2,          & ! Mean of s for 2nd normal distribution               [kg/kg]
      stdev_s1,    & ! Standard deviation of s for 1st normal distribution [kg/kg]
      stdev_s2,    & ! Standard deviation of s for 2nd normal distribution [kg/kg]
      rrtthl,      & ! Within-a-normal correlation of r_t and th_l             [-]
      alpha_thl,   & ! Factor relating to normalized variance for th_l         [-]
      alpha_rt       ! Factor relating to normalized variance for r_t          [-]
  end type pdf_parameter

end module pdf_parameter_module
