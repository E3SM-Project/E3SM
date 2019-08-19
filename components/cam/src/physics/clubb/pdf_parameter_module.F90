!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module pdf_parameter_module

  ! Description:
  ! This module defines the derived type pdf_parameter.

  ! References:
  !   None
  !-----------------------------------------------------------------------

  use clubb_precision, only: &
      core_rknd

  implicit none

  private ! Default scope

  public :: pdf_parameter,        & ! Variable Type(s)
            implicit_coefs_terms, &
            init_pdf_params         ! Procedure(s)

  ! CLUBB's PDF parameters.
  type pdf_parameter

    real( kind = core_rknd ), dimension(:), allocatable :: &
      w_1,             & ! Mean of w (1st PDF component)                   [m/s]
      w_2,             & ! Mean of w (2nd PDF component)                   [m/s]
      varnce_w_1,      & ! Variance of w (1st PDF component)           [m^2/s^2]
      varnce_w_2,      & ! Variance of w (2nd PDF component)           [m^2/s^2]
      rt_1,            & ! Mean of r_t (1st PDF component)               [kg/kg]
      rt_2,            & ! Mean of r_t (2nd PDF component)               [kg/kg]
      varnce_rt_1,     & ! Variance of r_t (1st PDF component)       [kg^2/kg^2]
      varnce_rt_2,     & ! Variance of r_t (2nd PDF component)       [kg^2/kg^2]
      thl_1,           & ! Mean of th_l (1st PDF component)                  [K]
      thl_2,           & ! Mean of th_l (2nd PDF component)                  [K]
      varnce_thl_1,    & ! Variance of th_l (1st PDF component)            [K^2]
      varnce_thl_2,    & ! Variance of th_l (2nd PDF component)            [K^2]
      corr_w_rt_1,     & ! Correlation of w and r_t (1st PDF component)      [-]
      corr_w_rt_2,     & ! Correlation of w and r_t (2nd PDF component)      [-]
      corr_w_thl_1,    & ! Correlation of w and th_l (1st PDF component)     [-]
      corr_w_thl_2,    & ! Correlation of w and th_l (2nd PDF component)     [-]
      corr_rt_thl_1,   & ! Correlation of r_t and th_l (1st PDF component)   [-]
      corr_rt_thl_2,   & ! Correlation of r_t and th_l (2nd PDF component)   [-]
      alpha_thl,       & ! Factor relating to normalized variance for th_l   [-]
      alpha_rt,        & ! Factor relating to normalized variance for r_t    [-]
      crt_1,           & ! r_t coef. in chi/eta eqns. (1st PDF comp.)        [-]
      crt_2,           & ! r_t coef. in chi/eta eqns. (2nd PDF comp.)        [-]
      cthl_1,          & ! th_l coef.: chi/eta eqns. (1st PDF comp.) [(kg/kg)/K]
      cthl_2,          & ! th_l coef.: chi/eta eqns. (2nd PDF comp.) [(kg/kg)/K]
      chi_1,           & ! Mean of chi (old s) (1st PDF component)       [kg/kg]
      chi_2,           & ! Mean of chi (old s) (2nd PDF component)       [kg/kg]
      stdev_chi_1,     & ! Standard deviation of chi (1st PDF component) [kg/kg]
      stdev_chi_2,     & ! Standard deviation of chi (2nd PDF component) [kg/kg]
      stdev_eta_1,     & ! Standard dev. of eta (old t) (1st PDF comp.)  [kg/kg]
      stdev_eta_2,     & ! Standard dev. of eta (old t) (2nd PDF comp.)  [kg/kg]
      covar_chi_eta_1, & ! Covariance of chi and eta (1st PDF comp.) [kg^2/kg^2]
      covar_chi_eta_2, & ! Covariance of chi and eta (2nd PDF comp.) [kg^2/kg^2]
      corr_w_chi_1,    & ! Correlation of w and chi (1st PDF component)      [-]
      corr_w_chi_2,    & ! Correlation of w and chi (2nd PDF component)      [-]
      corr_w_eta_1,    & ! Correlation of w and eta (1st PDF component)      [-]
      corr_w_eta_2,    & ! Correlation of w and eta (2nd PDF component)      [-]
      corr_chi_eta_1,  & ! Correlation of chi and eta (1st PDF component)    [-]
      corr_chi_eta_2,  & ! Correlation of chi and eta (2nd PDF component)    [-]
      rsatl_1,         & ! Saturation mixing ratio r_sat(mu_Tl_1,p)      [kg/kg]
      rsatl_2,         & ! Saturation mixing ratio r_sat(mu_Tl_2,p)      [kg/kg]
      rc_1,            & ! Mean of r_c (1st PDF component)               [kg/kg]
      rc_2,            & ! Mean of r_c (2nd PDF component)               [kg/kg]
      cloud_frac_1,    & ! Cloud fraction (1st PDF component)                [-]
      cloud_frac_2,    & ! Cloud fraction (2nd PDF component)                [-]
      mixt_frac,       & ! Weight of 1st PDF component (Sk_w dependent)      [-]
      ice_supersat_frac_1, & ! Ice supersaturation fraction (1st PDF comp.)  [-]
      ice_supersat_frac_2    ! Ice supersaturation fraction (2nd PDF comp.)  [-]

  end type pdf_parameter
  
  ! The implicit coefficients, semi-implicit coefficients and terms, and
  ! explicit terms for turbulent advection of turbulent fields are calculated
  ! from the PDF and the resulting PDF parameters.
  type implicit_coefs_terms

    real ( kind = core_rknd ) :: &
      coef_wp4_implicit,     & ! <w'^4> = coef_wp4_implicit * <w'^2>^2       [-]
      coef_wprtp2_implicit,  & ! <w'rt'^2> = coef_wprtp2_implicit*<rt'^2>  [m/s]
      coef_wpthlp2_implicit    ! <w'thl'^2>=coef_wpthlp2_implicit*<thl'^2> [m/s]

    ! <w'^2 rt'> = coef_wp2rtp_implicit * <w'rt'> + term_wp2rtp_explicit
    real ( kind = core_rknd ) :: &
      coef_wp2rtp_implicit, & ! Coefficient that is multiplied by <w'rt'>  [m/s]
      term_wp2rtp_explicit    ! Term that is on the RHS          [m^2/s^2 kg/kg]

    ! <w'^2 thl'> = coef_wp2thlp_implicit * <w'thl'> + term_wp2thlp_explicit
    real ( kind = core_rknd ) :: &
      coef_wp2thlp_implicit, & ! Coef. that is multiplied by <w'thl'>      [m/s]
      term_wp2thlp_explicit    ! Term that is on the RHS             [m^2/s^2 K]

    ! <w'rt'thl'> = coef_wprtpthlp_implicit*<rt'thl'> + term_wprtpthlp_explicit
    real ( kind = core_rknd ) :: &
      coef_wprtpthlp_implicit, & ! Coef. that is multiplied by <rt'thl'>   [m/s]
      term_wprtpthlp_explicit    ! Term that is on the RHS         [m/s(kg/kg)K]

  end type implicit_coefs_terms

! The CLUBB_CAM preprocessor directives are being commented out because this
! code is now also used for WRF-CLUBB.
!#ifdef CLUBB_CAM /* Code for storing pdf_parameter structs in pbuf as array */

  public :: pack_pdf_params, unpack_pdf_params

  integer, public, parameter :: num_pdf_params = 47

!#endif /* CLUBB_CAM */

  contains
  
  !=============================================================================
  subroutine init_pdf_params( nz, pdf_params )

    ! Description:
    ! Initializes all PDF parameters in the variable type pdf_parameter.

    ! References:
    !--------------------------------------------------------------------

    use constants_clubb, only: &
        zero    ! Constant(s)

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: &
      nz    ! Number of vertical grid levels    [-]

    ! Output Variable(s)
    type(pdf_parameter), intent(out) :: &
      pdf_params    ! PDF parameters            [units vary]

    allocate( pdf_params%w_1(nz), &
              pdf_params%w_2(nz), &
              pdf_params%varnce_w_1(nz), &
              pdf_params%varnce_w_2(nz), &
              pdf_params%rt_1(nz), &
              pdf_params%rt_2(nz), &
              pdf_params%varnce_rt_1(nz), &
              pdf_params%varnce_rt_2(nz), &
              pdf_params%thl_1(nz), &
              pdf_params%thl_2(nz), &
              pdf_params%varnce_thl_1(nz), &
              pdf_params%varnce_thl_2(nz), &
              pdf_params%corr_w_rt_1(nz), &
              pdf_params%corr_w_rt_2(nz), &
              pdf_params%corr_w_thl_1(nz), &
              pdf_params%corr_w_thl_2(nz), &
              pdf_params%corr_rt_thl_1(nz), &
              pdf_params%corr_rt_thl_2(nz), &
              pdf_params%alpha_thl(nz), &
              pdf_params%alpha_rt(nz), &
              pdf_params%crt_1(nz), &
              pdf_params%crt_2(nz), &
              pdf_params%cthl_1(nz), &
              pdf_params%cthl_2(nz), &
              pdf_params%chi_1(nz), &
              pdf_params%chi_2(nz), &
              pdf_params%stdev_chi_1(nz), &
              pdf_params%stdev_chi_2(nz), &
              pdf_params%stdev_eta_1(nz), &
              pdf_params%stdev_eta_2(nz), &
              pdf_params%covar_chi_eta_1(nz), &
              pdf_params%covar_chi_eta_2(nz), &
              pdf_params%corr_w_chi_1(nz), & 
              pdf_params%corr_w_chi_2(nz), &
              pdf_params%corr_w_eta_1(nz), & 
              pdf_params%corr_w_eta_2(nz), & 
              pdf_params%corr_chi_eta_1(nz), & 
              pdf_params%corr_chi_eta_2(nz), &
              pdf_params%rsatl_1(nz), &
              pdf_params%rsatl_2(nz), &
              pdf_params%rc_1(nz), &
              pdf_params%rc_2(nz), &
              pdf_params%cloud_frac_1(nz), &
              pdf_params%cloud_frac_2(nz), &
              pdf_params%mixt_frac(nz), &
              pdf_params%ice_supersat_frac_1(nz), &
              pdf_params%ice_supersat_frac_2(nz) )

    pdf_params%w_1 = zero
    pdf_params%w_2 = zero
    pdf_params%varnce_w_1 = zero
    pdf_params%varnce_w_2 = zero
    pdf_params%rt_1 = zero
    pdf_params%rt_2 = zero
    pdf_params%varnce_rt_1 = zero
    pdf_params%varnce_rt_2 = zero
    pdf_params%thl_1 = zero
    pdf_params%thl_2 = zero
    pdf_params%varnce_thl_1 = zero
    pdf_params%varnce_thl_2 = zero
    pdf_params%corr_w_rt_1 = zero
    pdf_params%corr_w_rt_2 = zero
    pdf_params%corr_w_thl_1 = zero
    pdf_params%corr_w_thl_2 = zero
    pdf_params%corr_rt_thl_1 = zero
    pdf_params%corr_rt_thl_2 = zero
    pdf_params%alpha_thl = zero
    pdf_params%alpha_rt = zero
    pdf_params%crt_1 = zero
    pdf_params%crt_2 = zero
    pdf_params%cthl_1 = zero
    pdf_params%cthl_2 = zero
    pdf_params%chi_1 = zero
    pdf_params%chi_2 = zero
    pdf_params%stdev_chi_1 = zero
    pdf_params%stdev_chi_2 = zero
    pdf_params%stdev_eta_1 = zero
    pdf_params%stdev_eta_2 = zero
    pdf_params%covar_chi_eta_1 = zero
    pdf_params%covar_chi_eta_2 = zero
    pdf_params%corr_w_chi_1 = zero 
    pdf_params%corr_w_chi_2 = zero 
    pdf_params%corr_w_eta_1 = zero 
    pdf_params%corr_w_eta_2 = zero 
    pdf_params%corr_chi_eta_1 = zero 
    pdf_params%corr_chi_eta_2 = zero 
    pdf_params%rsatl_1 = zero
    pdf_params%rsatl_2 = zero
    pdf_params%rc_1 = zero
    pdf_params%rc_2 = zero
    pdf_params%cloud_frac_1 = zero
    pdf_params%cloud_frac_2 = zero
    pdf_params%mixt_frac = zero
    pdf_params%ice_supersat_frac_1 = zero
    pdf_params%ice_supersat_frac_2 = zero


    return

  end subroutine init_pdf_params
  
  !=============================================================================

! The CLUBB_CAM preprocessor directives are being commented out because this
! code is now also used for WRF-CLUBB.
!#ifdef CLUBB_CAM /* Code for storing pdf_parameter structs in pbuf as array */

  subroutine pack_pdf_params(pdf_params, nz, r_param_array, &
                             k_start_in, k_end_in )
    implicit none
    
    integer, intent(in) :: nz ! Num Vert Model Levs
    
    ! Input a pdf_parameter array with nz instances of pdf_parameter
    type (pdf_parameter), intent(in) :: pdf_params

    ! Output a two dimensional real array with all values
    real (kind = core_rknd), dimension(nz,num_pdf_params), intent(out) :: &
       r_param_array  

    integer, optional, intent(in) :: k_start_in, k_end_in
    
    integer :: k_start, k_end
    
    if( present( k_start_in ) .and. present( k_end_in ) ) then
        k_start = k_start_in
        k_end = k_end_in
    else
        k_start = 1
        k_end = nz
    end if

    r_param_array(:,1) = pdf_params%w_1(k_start:k_end)
    r_param_array(:,2) = pdf_params%w_2(k_start:k_end)
    r_param_array(:,3) = pdf_params%varnce_w_1(k_start:k_end)
    r_param_array(:,4) = pdf_params%varnce_w_2(k_start:k_end)
    r_param_array(:,5) = pdf_params%rt_1(k_start:k_end)
    r_param_array(:,6) = pdf_params%rt_2(k_start:k_end)
    r_param_array(:,7) = pdf_params%varnce_rt_1(k_start:k_end)
    r_param_array(:,8) = pdf_params%varnce_rt_2(k_start:k_end)
    r_param_array(:,9) = pdf_params%thl_1(k_start:k_end)
    r_param_array(:,10) = pdf_params%thl_2(k_start:k_end)
    r_param_array(:,11) = pdf_params%varnce_thl_1(k_start:k_end)
    r_param_array(:,12) = pdf_params%varnce_thl_2(k_start:k_end)
    r_param_array(:,13) = pdf_params%corr_w_rt_1(k_start:k_end)
    r_param_array(:,14) = pdf_params%corr_w_rt_2(k_start:k_end)
    r_param_array(:,15) = pdf_params%corr_w_thl_1(k_start:k_end)
    r_param_array(:,16) = pdf_params%corr_w_thl_2(k_start:k_end)
    r_param_array(:,17) = pdf_params%corr_rt_thl_1(k_start:k_end)
    r_param_array(:,18) = pdf_params%corr_rt_thl_2(k_start:k_end)
    r_param_array(:,19) = pdf_params%alpha_thl(k_start:k_end)
    r_param_array(:,20) = pdf_params%alpha_rt(k_start:k_end)
    r_param_array(:,21) = pdf_params%crt_1(k_start:k_end)
    r_param_array(:,22) = pdf_params%crt_2(k_start:k_end)
    r_param_array(:,23) = pdf_params%cthl_1(k_start:k_end)
    r_param_array(:,24) = pdf_params%cthl_2(k_start:k_end)
    r_param_array(:,25) = pdf_params%chi_1(k_start:k_end)
    r_param_array(:,26) = pdf_params%chi_2(k_start:k_end)
    r_param_array(:,27) = pdf_params%stdev_chi_1(k_start:k_end)
    r_param_array(:,28) = pdf_params%stdev_chi_2(k_start:k_end)
    r_param_array(:,29) = pdf_params%stdev_eta_1(k_start:k_end)
    r_param_array(:,30) = pdf_params%stdev_eta_2(k_start:k_end)
    r_param_array(:,31) = pdf_params%covar_chi_eta_1(k_start:k_end)
    r_param_array(:,32) = pdf_params%covar_chi_eta_2(k_start:k_end)
    r_param_array(:,33) = pdf_params%corr_w_chi_1(k_start:k_end) 
    r_param_array(:,34) = pdf_params%corr_w_chi_2(k_start:k_end) 
    r_param_array(:,35) = pdf_params%corr_w_eta_1(k_start:k_end) 
    r_param_array(:,36) = pdf_params%corr_w_eta_2(k_start:k_end) 
    r_param_array(:,37) = pdf_params%corr_chi_eta_1(k_start:k_end) 
    r_param_array(:,38) = pdf_params%corr_chi_eta_2(k_start:k_end) 
    r_param_array(:,39) = pdf_params%rsatl_1(k_start:k_end)
    r_param_array(:,40) = pdf_params%rsatl_2(k_start:k_end)
    r_param_array(:,41) = pdf_params%rc_1(k_start:k_end)
    r_param_array(:,42) = pdf_params%rc_2(k_start:k_end)
    r_param_array(:,43) = pdf_params%cloud_frac_1(k_start:k_end)
    r_param_array(:,44) = pdf_params%cloud_frac_2(k_start:k_end)
    r_param_array(:,45) = pdf_params%mixt_frac(k_start:k_end)
    r_param_array(:,46) = pdf_params%ice_supersat_frac_1(k_start:k_end)
    r_param_array(:,47) = pdf_params%ice_supersat_frac_2(k_start:k_end)

  end subroutine pack_pdf_params

  subroutine unpack_pdf_params(r_param_array, nz, pdf_params, &
                               k_start_in, k_end_in )
    implicit none
    
    integer, intent(in) :: nz ! Num Vert Model Levs
    
    ! Input a two dimensional real array with pdf values
    real (kind = core_rknd), dimension(nz,num_pdf_params), intent(in) :: &
       r_param_array 

    ! Output a pdf_parameter array with nz instances of pdf_parameter
    type (pdf_parameter), intent(inout) :: pdf_params
    
    integer, optional, intent(in) :: k_start_in, k_end_in
    
    integer :: k_start, k_end
    
    if( present( k_start_in ) .and. present( k_end_in ) ) then
        k_start = k_start_in
        k_end = k_end_in
    else
        k_start = 1
        k_end = nz
    end if
    
    pdf_params%w_1(k_start:k_end) = r_param_array(:,1)
    pdf_params%w_2(k_start:k_end) = r_param_array(:,2)
    pdf_params%varnce_w_1(k_start:k_end) = r_param_array(:,3)
    pdf_params%varnce_w_2(k_start:k_end) = r_param_array(:,4)
    pdf_params%rt_1(k_start:k_end) = r_param_array(:,5)
    pdf_params%rt_2(k_start:k_end) = r_param_array(:,6)
    pdf_params%varnce_rt_1(k_start:k_end) = r_param_array(:,7)
    pdf_params%varnce_rt_2(k_start:k_end)= r_param_array(:,8)
    pdf_params%thl_1(k_start:k_end) = r_param_array(:,9)
    pdf_params%thl_2(k_start:k_end) = r_param_array(:,10)
    pdf_params%varnce_thl_1(k_start:k_end) = r_param_array(:,11)
    pdf_params%varnce_thl_2(k_start:k_end) = r_param_array(:,12)
    pdf_params%corr_w_rt_1(k_start:k_end) = r_param_array(:,13)
    pdf_params%corr_w_rt_2(k_start:k_end) = r_param_array(:,14)
    pdf_params%corr_w_thl_1(k_start:k_end) = r_param_array(:,15)
    pdf_params%corr_w_thl_2(k_start:k_end) = r_param_array(:,16)
    pdf_params%corr_rt_thl_1(k_start:k_end) = r_param_array(:,17)
    pdf_params%corr_rt_thl_2(k_start:k_end) = r_param_array(:,18)
    pdf_params%alpha_thl(k_start:k_end) = r_param_array(:,19)
    pdf_params%alpha_rt(k_start:k_end) = r_param_array(:,20)
    pdf_params%crt_1(k_start:k_end) = r_param_array(:,21)
    pdf_params%crt_2(k_start:k_end) = r_param_array(:,22)
    pdf_params%cthl_1(k_start:k_end) = r_param_array(:,23)
    pdf_params%cthl_2(k_start:k_end) = r_param_array(:,24)
    pdf_params%chi_1(k_start:k_end) = r_param_array(:,25)
    pdf_params%chi_2(k_start:k_end) = r_param_array(:,26)
    pdf_params%stdev_chi_1(k_start:k_end) = r_param_array(:,27)
    pdf_params%stdev_chi_2(k_start:k_end) = r_param_array(:,28)
    pdf_params%stdev_eta_1(k_start:k_end) = r_param_array(:,29)
    pdf_params%stdev_eta_2(k_start:k_end) = r_param_array(:,30)
    pdf_params%covar_chi_eta_1(k_start:k_end) = r_param_array(:,31)
    pdf_params%covar_chi_eta_2(k_start:k_end) = r_param_array(:,32)
    pdf_params%corr_w_chi_1(k_start:k_end) = r_param_array(:,33)
    pdf_params%corr_w_chi_2(k_start:k_end) = r_param_array(:,34)
    pdf_params%corr_w_eta_1(k_start:k_end) = r_param_array(:,35)
    pdf_params%corr_w_eta_2(k_start:k_end) = r_param_array(:,36)
    pdf_params%corr_chi_eta_1(k_start:k_end) = r_param_array(:,37)
    pdf_params%corr_chi_eta_2(k_start:k_end) = r_param_array(:,38)
    pdf_params%rsatl_1(k_start:k_end) = r_param_array(:,39)
    pdf_params%rsatl_2(k_start:k_end) = r_param_array(:,40)
    pdf_params%rc_1(k_start:k_end) = r_param_array(:,41)
    pdf_params%rc_2(k_start:k_end) = r_param_array(:,42)
    pdf_params%cloud_frac_1(k_start:k_end) = r_param_array(:,43)
    pdf_params%cloud_frac_2(k_start:k_end) = r_param_array(:,44)
    pdf_params%mixt_frac(k_start:k_end) = r_param_array(:,45)
    pdf_params%ice_supersat_frac_1(k_start:k_end) = r_param_array(:,46)
    pdf_params%ice_supersat_frac_2(k_start:k_end) = r_param_array(:,47)

  end subroutine unpack_pdf_params

!#endif /* CLUBB_CAM */

end module pdf_parameter_module
