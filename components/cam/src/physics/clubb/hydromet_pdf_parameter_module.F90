!---------------------------------------------------------------------------
! $Id: hydromet_pdf_parameter_module.F90 7284 2014-09-11 02:52:58Z bmg2@uwm.edu $
!===============================================================================
module hydromet_pdf_parameter_module

  ! Description:
  ! This module defines the derived type hydromet_pdf_parameter.

  ! References:
  !   None
  !-------------------------------------------------------------------------

  use clubb_precision, only: &
      core_rknd  ! Variable(s)

  implicit none

  private ! Default scope

  public :: hydromet_pdf_parameter,   & ! Variable type
            init_hydromet_pdf_params    ! Procedure

  integer, parameter, private :: &
    max_hydromet_dim = 8

  type hydromet_pdf_parameter

    real( kind = core_rknd ), dimension(max_hydromet_dim) :: &
      hm1,           & ! Mean of hydrometeor, hm (1st PDF component)   [un vary]
      hm2,           & ! Mean of hydrometeor, hm (2nd PDF component)   [un vary]
      mu_hm_1,       & ! Mean of hm (1st PDF component) in-precip (ip) [un vary]
      mu_hm_2,       & ! Mean of hm (2nd PDF component) ip             [un vary]
      sigma_hm_1,    & ! Standard deviation of hm (1st PDF comp.) ip   [un vary]
      sigma_hm_2,    & ! Standard deviation of hm (2nd PDF comp.) ip   [un vary]
      corr_w_hm_1,   & ! Correlation of w and hm (1st PDF component) ip      [-]
      corr_w_hm_2,   & ! Correlation of w and hm (2nd PDF component) ip      [-]
      corr_chi_hm_1, & ! Correlation of chi and hm (1st PDF component) ip    [-]
      corr_chi_hm_2, & ! Correlation of chi and hm (2nd PDF component) ip    [-]
      corr_eta_hm_1, & ! Correlation of eta and hm (1st PDF component) ip    [-]
      corr_eta_hm_2    ! Correlation of eta and hm (2nd PDF component) ip    [-]

    real( kind = core_rknd ) :: &
      mu_Ncn_1,    & ! Mean of Ncn (1st PDF component)                  [num/kg]
      mu_Ncn_2,    & ! Mean of Ncn (2nd PDF component)                  [num/kg]
      sigma_Ncn_1, & ! Standard deviation of Ncn (1st PDF component)    [num/kg]
      sigma_Ncn_2    ! Standard deviation of Ncn (2nd PDF component)    [num/kg]

    real( kind = core_rknd ) :: &
      precip_frac,   & ! Precipitation fraction (overall)           [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component) [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component) [-]

  end type hydromet_pdf_parameter

contains

  !=============================================================================
  subroutine init_hydromet_pdf_params( hydromet_pdf_params )

    ! Description:
    ! Initialize the elements of hydromet_pdf_params.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        zero  ! Constant(s)

    implicit none

    ! Output Variable
    type(hydromet_pdf_parameter), intent(out) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters      [units vary]

    ! Initialize hydromet_pdf_params.
    hydromet_pdf_params%hm1 = zero
    hydromet_pdf_params%hm2 = zero
    hydromet_pdf_params%mu_hm_1 = zero
    hydromet_pdf_params%mu_hm_2 = zero
    hydromet_pdf_params%sigma_hm_1 = zero
    hydromet_pdf_params%sigma_hm_2 = zero
    hydromet_pdf_params%corr_w_hm_1 = zero
    hydromet_pdf_params%corr_w_hm_2 = zero
    hydromet_pdf_params%corr_chi_hm_1 = zero
    hydromet_pdf_params%corr_chi_hm_2 = zero
    hydromet_pdf_params%corr_eta_hm_1 = zero
    hydromet_pdf_params%corr_eta_hm_2 = zero

    hydromet_pdf_params%mu_Ncn_1 = zero
    hydromet_pdf_params%mu_Ncn_2 = zero
    hydromet_pdf_params%sigma_Ncn_1 = zero
    hydromet_pdf_params%sigma_Ncn_2 = zero

    hydromet_pdf_params%precip_frac = zero
    hydromet_pdf_params%precip_frac_1 = zero
    hydromet_pdf_params%precip_frac_2 = zero


    return

  end subroutine init_hydromet_pdf_params

!===============================================================================

end module hydromet_pdf_parameter_module
