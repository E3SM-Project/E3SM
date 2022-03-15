!---------------------------------------------------------------------------
! $Id$
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
            precipitation_fractions, &
            init_hydromet_pdf_params, & ! Procedure
            init_precip_fracs

  integer, parameter, private :: &
    max_hydromet_dim = 8

  type hydromet_pdf_parameter

    real( kind = core_rknd ), dimension(max_hydromet_dim) :: &
      hm_1,          & ! Mean of hydrometeor, hm (1st PDF component)   [un vary]
      hm_2,          & ! Mean of hydrometeor, hm (2nd PDF component)   [un vary]
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

    real( kind = core_rknd ), dimension(max_hydromet_dim,max_hydromet_dim) :: &
      corr_hmx_hmy_1, & ! Correlation of hmx and hmy (1st PDF component) ip  [-]
      corr_hmx_hmy_2    ! Correlation of hmx and hmy (2nd PDF component) ip  [-]

    real( kind = core_rknd ) :: &
      mu_Ncn_1,    & ! Mean of Ncn (1st PDF component)                  [num/kg]
      mu_Ncn_2,    & ! Mean of Ncn (2nd PDF component)                  [num/kg]
      sigma_Ncn_1, & ! Standard deviation of Ncn (1st PDF component)    [num/kg]
      sigma_Ncn_2    ! Standard deviation of Ncn (2nd PDF component)    [num/kg]

  end type hydromet_pdf_parameter
  
  type precipitation_fractions
    
    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      precip_frac,   & ! Precipitation fraction (overall)           [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component) [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component) [-]
      
  end type 
    

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
    hydromet_pdf_params%hm_1 = zero
    hydromet_pdf_params%hm_2 = zero
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

    hydromet_pdf_params%corr_hmx_hmy_1 = zero
    hydromet_pdf_params%corr_hmx_hmy_2 = zero

    hydromet_pdf_params%mu_Ncn_1 = zero
    hydromet_pdf_params%mu_Ncn_2 = zero
    hydromet_pdf_params%sigma_Ncn_1 = zero
    hydromet_pdf_params%sigma_Ncn_2 = zero

    return

  end subroutine init_hydromet_pdf_params
  
  !=============================================================================
  subroutine init_precip_fracs( nz, ngrdcol, &
                                precip_fracs )

    ! Description:
    ! Initialize the elements of precip_fracs.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        zero  ! Constant(s)

    implicit none
    
    ! Input Variable(s)
    integer, intent(in) :: &
      nz,     & ! Number of vertical grid levels    [-]
      ngrdcol   ! Number of grid columns            [-]

    ! Output Variable
    type(precipitation_fractions), intent(out) :: &
      precip_fracs    ! Hydrometeor PDF parameters      [units vary]

    ! Allocate precip frac arrays
    allocate( precip_fracs%precip_frac(ngrdcol,nz), &
              precip_fracs%precip_frac_1(ngrdcol,nz), &
              precip_fracs%precip_frac_2(ngrdcol,nz)  )

    ! Initialize precip_fracs.
    precip_fracs%precip_frac   = zero
    precip_fracs%precip_frac_1 = zero
    precip_fracs%precip_frac_2 = zero

    return

  end subroutine init_precip_fracs

!===============================================================================

end module hydromet_pdf_parameter_module
