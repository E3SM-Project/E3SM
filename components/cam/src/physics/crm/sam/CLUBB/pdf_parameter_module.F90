! $Id: pdf_parameter_module.F90 5668 2012-01-29 03:40:28Z bmg2@uwm.edu $
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
      w1,          & ! Mean of w (1st PDF component)                       [m/s]
      w2,          & ! Mean of w (2nd PDF component)                       [m/s]
      varnce_w1,   & ! Variance of w (1st PDF component)               [m^2/s^2]
      varnce_w2,   & ! Variance of w (2nd PDF component)               [m^2/s^2]
      rt1,         & ! Mean of r_t (1st PDF component)                   [kg/kg]
      rt2,         & ! Mean of r_t (2nd PDF component)                   [kg/kg]
      varnce_rt1,  & ! Variance of r_t (1st PDF component)           [kg^2/kg^2]
      varnce_rt2,  & ! Variance of r_t (2nd PDF component)           [kg^2/kg^2]
      thl1,        & ! Mean of th_l (1st PDF component)                      [K]
      thl2,        & ! Mean of th_l (2nd PDF component)                      [K]
      varnce_thl1, & ! Variance of th_l (1st PDF component)                [K^2]
      varnce_thl2, & ! Variance of th_l (2nd PDF component)                [K^2]
      rrtthl,      & ! Correlation between r_t and th_l (both components)    [-]
      alpha_thl,   & ! Factor relating to normalized variance for th_l       [-]
      alpha_rt,    & ! Factor relating to normalized variance for r_t        [-]
      crt1,        & ! Coef. on r_t in s/t eqns. (1st PDF comp.)             [-]
      crt2,        & ! Coef. on r_t in s/t eqns. (2nd PDF comp.)             [-]
      cthl1,       & ! Coef. on th_l in s/t eqns. (1st PDF comp.)    [(kg/kg)/K]
      cthl2,       & ! Coef. on th_l in s/t eqns. (2nd PDF comp.)    [(kg/kg)/K]
      s1,          & ! Mean of s (1st PDF component)                     [kg/kg]
      s2,          & ! Mean of s (2nd PDF component)                     [kg/kg]
      stdev_s1,    & ! Standard deviation of s (1st PDF component)       [kg/kg]
      stdev_s2,    & ! Standard deviation of s (2nd PDF component)       [kg/kg]
      stdev_t1,    & ! Standard deviation of t (1st PDF component)       [kg/kg]
      stdev_t2,    & ! Standard deviation of t (2nd PDF component)       [kg/kg]
      covar_st_1,  & ! Covariance of s and t (1st PDF component)     [kg^2/kg^2]
      covar_st_2,  & ! Covariance of s and t (2nd PDF component)     [kg^2/kg^2]
      corr_st_1,   & ! Correlation between s and t (1st PDF component)       [-]
      corr_st_2,   & ! Correlation between s and t (2nd PDF component)       [-]
      rsl1,        & ! Mean of r_sl (1st PDF component)                  [kg/kg]
      rsl2,        & ! Mean of r_sl (2nd PDF component)                  [kg/kg]
      rc1,         & ! Mean of r_c (1st PDF component)                   [kg/kg]
      rc2,         & ! Mean of r_c (2nd PDF component)                   [kg/kg]
      cloud_frac1, & ! Cloud fraction (1st PDF component)                    [-]
      cloud_frac2, & ! Cloud fraction (2nd PDF component)                    [-]
      mixt_frac      ! Weight of 1st PDF component (Sk_w dependent)          [-]
  end type pdf_parameter

end module pdf_parameter_module
