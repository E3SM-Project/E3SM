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

  public :: pdf_parameter,                 & ! Variable Type(s)
            implicit_coefs_terms,          &
            init_pdf_params,               & ! Procedure(s)
            print_pdf_params,              &
            copy_single_pdf_params_to_multi, &
            copy_multi_pdf_params_to_single, &
            init_pdf_implicit_coefs_terms

  ! CLUBB's PDF parameters.
  type pdf_parameter

    real( kind = core_rknd ), dimension(:,:), allocatable :: &
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

    real( kind = core_rknd ), dimension(:), allocatable :: &
      coef_wp4_implicit    ! <w'^4> = coef_wp4_implicit * <w'^2>^2       [-]

    ! <w'^2 rt'> = coef_wp2rtp_implicit * <w'rt'> + term_wp2rtp_explicit
    real( kind = core_rknd ), dimension(:), allocatable :: &
      coef_wp2rtp_implicit, & ! Coefficient that is multiplied by <w'rt'>  [m/s]
      term_wp2rtp_explicit    ! Term that is on the RHS          [m^2/s^2 kg/kg]

    ! <w'^2 thl'> = coef_wp2thlp_implicit * <w'thl'> + term_wp2thlp_explicit
    real( kind = core_rknd ), dimension(:), allocatable :: &
      coef_wp2thlp_implicit, & ! Coef. that is multiplied by <w'thl'>      [m/s]
      term_wp2thlp_explicit    ! Term that is on the RHS             [m^2/s^2 K]

    ! <w'^2 u'> = coef_wp2up_implicit * <u'w'> + term_wp2up_explicit
    real( kind = core_rknd ), dimension(:), allocatable :: &
      coef_wp2up_implicit, & ! Coefficient that is multiplied by <u'w'>    [m/s]
      term_wp2up_explicit    ! Term that is on the RHS                 [m^3/s^3]

    ! <w'^2 v'> = coef_wp2vp_implicit * <v'w'> + term_wp2vp_explicit
    real( kind = core_rknd ), dimension(:), allocatable :: &
      coef_wp2vp_implicit, & ! Coefficient that is multiplied by <v'w'>    [m/s]
      term_wp2vp_explicit    ! Term that is on the RHS                 [m^3/s^3]

    ! <w'rt'^2> = coef_wprtp2_implicit * <rt'^2> + term_wprtp2_explicit
    real( kind = core_rknd ), dimension(:), allocatable :: &
      coef_wprtp2_implicit, & ! Coefficient that is multiplied by <rt'^2>  [m/s]
      term_wprtp2_explicit    ! Term that is on the RHS          [m/s kg^2/kg^2]

    ! <w'thl'^2> = coef_wpthlp2_implicit * <thl'^2> + term_wpthlp2_explicit
    real( kind = core_rknd ), dimension(:), allocatable :: &
      coef_wpthlp2_implicit, & ! Coef. that is multiplied by <thl'^2>      [m/s]
      term_wpthlp2_explicit    ! Term that is on the RHS               [m/s K^2]

    ! <w'rt'thl'> = coef_wprtpthlp_implicit*<rt'thl'> + term_wprtpthlp_explicit
    real( kind = core_rknd ), dimension(:), allocatable :: &
      coef_wprtpthlp_implicit, & ! Coef. that is multiplied by <rt'thl'>   [m/s]
      term_wprtpthlp_explicit    ! Term that is on the RHS         [m/s(kg/kg)K]

    ! <w'u'^2> = coef_wpup2_implicit * <u'^2> + term_wpup2_explicit
    real( kind = core_rknd ), dimension(:), allocatable :: &
      coef_wpup2_implicit, & ! Coefficient that is multiplied by <u'^2>    [m/s]
      term_wpup2_explicit    ! Term that is on the RHS                 [m^3/s^3]

    ! <w'v'^2> = coef_wpvp2_implicit * <v'^2> + term_wpvp2_explicit
    real( kind = core_rknd ), dimension(:), allocatable :: &
      coef_wpvp2_implicit, & ! Coefficient that is multiplied by <v'^2>    [m/s]
      term_wpvp2_explicit    ! Term that is on the RHS                 [m^3/s^3]

    ! <w'^2 sclr'> = coef_wp2sclrp_implicit * <w'sclr'> + term_wp2sclrp_explicit
    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      coef_wp2sclrp_implicit, & ! Coef. that is multiplied by <w'sclr'>    [m/s]
      term_wp2sclrp_explicit    ! Term that is on the RHS   [m^2/s^2 (un. vary)]

    ! <w'sclr'^2> = coef_wpsclrp2_implicit * <sclr'^2> + term_wpsclrp2_explicit
    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      coef_wpsclrp2_implicit, & ! Coef. that is multiplied by <sclr'^2>    [m/s]
      term_wpsclrp2_explicit    ! Term that is on the RHS    [m/s(units vary)^2]

    ! <w'rt'sclr'> = coef_wprtpsclrp_implicit * <sclr'rt'>
    !                + term_wprtpsclrp_explicit
    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      coef_wprtpsclrp_implicit, & ! Coef. that is multiplied by <sclr'rt'> [m/s]
      term_wprtpsclrp_explicit    ! Term that is on the RHS [m/s(kg/kg)(un. v.)]

    ! <w'thl'sclr'> = coef_wpthlpsclrp_implicit * <sclr'thl'>
    !                 + term_wpthlpsclrp_explicit
    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      coef_wpthlpsclrp_implicit, & ! Coef. that is mult. by <sclr'thl'>    [m/s]
      term_wpthlpsclrp_explicit    ! Term that is on the RHS  [(m/s)K(un. vary)]

  end type implicit_coefs_terms

! The CLUBB_CAM preprocessor directives are being commented out because this
! code is now also used for WRF-CLUBB.
!#ifdef CLUBB_CAM /* Code for storing pdf_parameter structs in pbuf as array */

  public :: pack_pdf_params, unpack_pdf_params

  integer, public, parameter :: num_pdf_params = 47

!#endif /* CLUBB_CAM */

  contains
  
  !=============================================================================
  subroutine init_pdf_params( nz, ngrdcol, &
                              pdf_params )

    ! Description:
    ! Initializes all PDF parameters in the variable type pdf_parameter.

    ! References:
    !--------------------------------------------------------------------

    use constants_clubb, only: &
        zero    ! Constant(s)

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: &
      nz,   & ! Number of vertical grid levels    [-]
      ngrdcol ! Number of grid columns            [-]

    ! Output Variable(s)
    type(pdf_parameter), intent(out) :: &
      pdf_params    ! PDF parameters            [units vary]

    allocate( pdf_params%w_1(ngrdcol,nz), &
              pdf_params%w_2(ngrdcol,nz), &
              pdf_params%varnce_w_1(ngrdcol,nz), &
              pdf_params%varnce_w_2(ngrdcol,nz), &
              pdf_params%rt_1(ngrdcol,nz), &
              pdf_params%rt_2(ngrdcol,nz), &
              pdf_params%varnce_rt_1(ngrdcol,nz), &
              pdf_params%varnce_rt_2(ngrdcol,nz), &
              pdf_params%thl_1(ngrdcol,nz), &
              pdf_params%thl_2(ngrdcol,nz), &
              pdf_params%varnce_thl_1(ngrdcol,nz), &
              pdf_params%varnce_thl_2(ngrdcol,nz), &
              pdf_params%corr_w_rt_1(ngrdcol,nz), &
              pdf_params%corr_w_rt_2(ngrdcol,nz), &
              pdf_params%corr_w_thl_1(ngrdcol,nz), &
              pdf_params%corr_w_thl_2(ngrdcol,nz), &
              pdf_params%corr_rt_thl_1(ngrdcol,nz), &
              pdf_params%corr_rt_thl_2(ngrdcol,nz), &
              pdf_params%alpha_thl(ngrdcol,nz), &
              pdf_params%alpha_rt(ngrdcol,nz), &
              pdf_params%crt_1(ngrdcol,nz), &
              pdf_params%crt_2(ngrdcol,nz), &
              pdf_params%cthl_1(ngrdcol,nz), &
              pdf_params%cthl_2(ngrdcol,nz), &
              pdf_params%chi_1(ngrdcol,nz), &
              pdf_params%chi_2(ngrdcol,nz), &
              pdf_params%stdev_chi_1(ngrdcol,nz), &
              pdf_params%stdev_chi_2(ngrdcol,nz), &
              pdf_params%stdev_eta_1(ngrdcol,nz), &
              pdf_params%stdev_eta_2(ngrdcol,nz), &
              pdf_params%covar_chi_eta_1(ngrdcol,nz), &
              pdf_params%covar_chi_eta_2(ngrdcol,nz), &
              pdf_params%corr_w_chi_1(ngrdcol,nz), & 
              pdf_params%corr_w_chi_2(ngrdcol,nz), &
              pdf_params%corr_w_eta_1(ngrdcol,nz), & 
              pdf_params%corr_w_eta_2(ngrdcol,nz), & 
              pdf_params%corr_chi_eta_1(ngrdcol,nz), & 
              pdf_params%corr_chi_eta_2(ngrdcol,nz), &
              pdf_params%rsatl_1(ngrdcol,nz), &
              pdf_params%rsatl_2(ngrdcol,nz), &
              pdf_params%rc_1(ngrdcol,nz), &
              pdf_params%rc_2(ngrdcol,nz), &
              pdf_params%cloud_frac_1(ngrdcol,nz), &
              pdf_params%cloud_frac_2(ngrdcol,nz), &
              pdf_params%mixt_frac(ngrdcol,nz), &
              pdf_params%ice_supersat_frac_1(ngrdcol,nz), &
              pdf_params%ice_supersat_frac_2(ngrdcol,nz) )

    pdf_params%w_1(:,:) = zero
    pdf_params%w_2(:,:) = zero
    pdf_params%varnce_w_1(:,:) = zero
    pdf_params%varnce_w_2(:,:) = zero
    pdf_params%rt_1(:,:) = zero
    pdf_params%rt_2(:,:) = zero
    pdf_params%varnce_rt_1(:,:) = zero
    pdf_params%varnce_rt_2(:,:) = zero
    pdf_params%thl_1(:,:) = zero
    pdf_params%thl_2(:,:) = zero
    pdf_params%varnce_thl_1(:,:) = zero
    pdf_params%varnce_thl_2(:,:) = zero
    pdf_params%corr_w_rt_1(:,:) = zero
    pdf_params%corr_w_rt_2(:,:) = zero
    pdf_params%corr_w_thl_1(:,:) = zero
    pdf_params%corr_w_thl_2(:,:) = zero
    pdf_params%corr_rt_thl_1(:,:) = zero
    pdf_params%corr_rt_thl_2(:,:) = zero
    pdf_params%alpha_thl(:,:) = zero
    pdf_params%alpha_rt(:,:) = zero
    pdf_params%crt_1(:,:) = zero
    pdf_params%crt_2(:,:) = zero
    pdf_params%cthl_1(:,:) = zero
    pdf_params%cthl_2(:,:) = zero
    pdf_params%chi_1(:,:) = zero
    pdf_params%chi_2(:,:) = zero
    pdf_params%stdev_chi_1(:,:) = zero
    pdf_params%stdev_chi_2(:,:) = zero
    pdf_params%stdev_eta_1(:,:) = zero
    pdf_params%stdev_eta_2(:,:) = zero
    pdf_params%covar_chi_eta_1(:,:) = zero
    pdf_params%covar_chi_eta_2(:,:) = zero
    pdf_params%corr_w_chi_1(:,:) = zero 
    pdf_params%corr_w_chi_2(:,:) = zero 
    pdf_params%corr_w_eta_1(:,:) = zero 
    pdf_params%corr_w_eta_2(:,:) = zero 
    pdf_params%corr_chi_eta_1(:,:) = zero 
    pdf_params%corr_chi_eta_2(:,:) = zero 
    pdf_params%rsatl_1(:,:) = zero
    pdf_params%rsatl_2(:,:) = zero
    pdf_params%rc_1(:,:) = zero
    pdf_params%rc_2(:,:) = zero
    pdf_params%cloud_frac_1(:,:) = zero
    pdf_params%cloud_frac_2(:,:) = zero
    pdf_params%mixt_frac(:,:) = zero
    pdf_params%ice_supersat_frac_1(:,:) = zero
    pdf_params%ice_supersat_frac_2(:,:) = zero


    return

  end subroutine init_pdf_params
  
  !=============================================================================
  subroutine init_pdf_implicit_coefs_terms( nz, sclr_dim, &
                                            pdf_implicit_coefs_terms )

    ! Description:
    ! Initializes all PDF implicit coefficients and explicit terms in the
    ! variable type implicit_coefs_terms.

    ! References:
    !--------------------------------------------------------------------

    use constants_clubb, only: &
        zero    ! Constant(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,       & ! Number of vertical grid levels    [-]
      sclr_dim    ! Number of scalar variables        [-]

    ! Output Variable
    type(implicit_coefs_terms), intent(out) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]


    ! Allocate pdf_implicit_coefs_terms
    allocate( pdf_implicit_coefs_terms%coef_wp4_implicit(1:nz), &
              pdf_implicit_coefs_terms%coef_wp2rtp_implicit(1:nz), &
              pdf_implicit_coefs_terms%term_wp2rtp_explicit(1:nz), &
              pdf_implicit_coefs_terms%coef_wp2thlp_implicit(1:nz), &
              pdf_implicit_coefs_terms%term_wp2thlp_explicit(1:nz), &
              pdf_implicit_coefs_terms%coef_wp2up_implicit(1:nz), &
              pdf_implicit_coefs_terms%term_wp2up_explicit(1:nz), &
              pdf_implicit_coefs_terms%coef_wp2vp_implicit(1:nz), &
              pdf_implicit_coefs_terms%term_wp2vp_explicit(1:nz), &
              pdf_implicit_coefs_terms%coef_wprtp2_implicit(1:nz), &
              pdf_implicit_coefs_terms%term_wprtp2_explicit(1:nz), &
              pdf_implicit_coefs_terms%coef_wpthlp2_implicit(1:nz), &
              pdf_implicit_coefs_terms%term_wpthlp2_explicit(1:nz), &
              pdf_implicit_coefs_terms%coef_wprtpthlp_implicit(1:nz), &
              pdf_implicit_coefs_terms%term_wprtpthlp_explicit(1:nz), &
              pdf_implicit_coefs_terms%coef_wpup2_implicit(1:nz), &
              pdf_implicit_coefs_terms%term_wpup2_explicit(1:nz), &
              pdf_implicit_coefs_terms%coef_wpvp2_implicit(1:nz), &
              pdf_implicit_coefs_terms%term_wpvp2_explicit(1:nz) )

    if ( sclr_dim > 0 ) then
       allocate( &
          pdf_implicit_coefs_terms%coef_wp2sclrp_implicit(1:nz,1:sclr_dim), &
          pdf_implicit_coefs_terms%term_wp2sclrp_explicit(1:nz,1:sclr_dim), &
          pdf_implicit_coefs_terms%coef_wpsclrp2_implicit(1:nz,1:sclr_dim), &
          pdf_implicit_coefs_terms%term_wpsclrp2_explicit(1:nz,1:sclr_dim), &
          pdf_implicit_coefs_terms%coef_wprtpsclrp_implicit(1:nz,1:sclr_dim), &
          pdf_implicit_coefs_terms%term_wprtpsclrp_explicit(1:nz,1:sclr_dim), &
          pdf_implicit_coefs_terms%coef_wpthlpsclrp_implicit(1:nz,1:sclr_dim), &
          pdf_implicit_coefs_terms%term_wpthlpsclrp_explicit(1:nz,1:sclr_dim) )
    endif ! sclr_dim > 0

    ! Initialize pdf_implicit_coefs_terms
    pdf_implicit_coefs_terms%coef_wp4_implicit = zero
    pdf_implicit_coefs_terms%coef_wp2rtp_implicit = zero
    pdf_implicit_coefs_terms%term_wp2rtp_explicit = zero
    pdf_implicit_coefs_terms%coef_wp2thlp_implicit = zero
    pdf_implicit_coefs_terms%term_wp2thlp_explicit = zero
    pdf_implicit_coefs_terms%coef_wp2up_implicit = zero
    pdf_implicit_coefs_terms%term_wp2up_explicit = zero
    pdf_implicit_coefs_terms%coef_wp2vp_implicit = zero
    pdf_implicit_coefs_terms%term_wp2vp_explicit = zero
    pdf_implicit_coefs_terms%coef_wprtp2_implicit = zero
    pdf_implicit_coefs_terms%term_wprtp2_explicit = zero
    pdf_implicit_coefs_terms%coef_wpthlp2_implicit = zero
    pdf_implicit_coefs_terms%term_wpthlp2_explicit = zero
    pdf_implicit_coefs_terms%coef_wprtpthlp_implicit = zero
    pdf_implicit_coefs_terms%term_wprtpthlp_explicit = zero
    pdf_implicit_coefs_terms%coef_wpup2_implicit = zero
    pdf_implicit_coefs_terms%term_wpup2_explicit = zero
    pdf_implicit_coefs_terms%coef_wpvp2_implicit = zero
    pdf_implicit_coefs_terms%term_wpvp2_explicit = zero
    if ( sclr_dim > 0 ) then
       pdf_implicit_coefs_terms%coef_wp2sclrp_implicit = zero
       pdf_implicit_coefs_terms%term_wp2sclrp_explicit = zero
       pdf_implicit_coefs_terms%coef_wpsclrp2_implicit = zero
       pdf_implicit_coefs_terms%term_wpsclrp2_explicit = zero
       pdf_implicit_coefs_terms%coef_wprtpsclrp_implicit = zero
       pdf_implicit_coefs_terms%term_wprtpsclrp_explicit = zero
       pdf_implicit_coefs_terms%coef_wpthlpsclrp_implicit = zero
       pdf_implicit_coefs_terms%term_wpthlpsclrp_explicit = zero
    endif ! sclr_dim > 0


    return

  end subroutine init_pdf_implicit_coefs_terms

  !=============================================================================
  subroutine print_pdf_params( pdf_params, ngrdcol )

    ! Description:
    ! Prints every value found in type pdf_parameter

    !-------------------------------------------------------------------

    use constants_clubb, only: &
        fstderr

    implicit none

    ! Input Variable(s)
    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters            [units vary]

    integer, intent(in) :: &
      ngrdcol    ! Number of horizontal grid columns

    ! Local Variable(s)
    integer :: i    ! Loop index


    do i = 1, ngrdcol

       write(fstderr,*) "PDF parameters for column ", i
       write(fstderr,*) "---------------------------------------------------"
       write(fstderr,*) "w_1(i,:) = ", pdf_params%w_1(i,:)
       write(fstderr,*) "w_2(i,:) = ", pdf_params%w_2(i,:)
       write(fstderr,*) "varnce_w_1(i,:) = ", pdf_params%varnce_w_1(i,:)
       write(fstderr,*) "varnce_w_2(i,:) = ", pdf_params%varnce_w_2(i,:)
       write(fstderr,*) "rt_1(i,:) = ", pdf_params%rt_1(i,:)
       write(fstderr,*) "rt_2(i,:) = ", pdf_params%rt_2(i,:)
       write(fstderr,*) "varnce_rt_1(i,:) = ", pdf_params%varnce_rt_1(i,:)
       write(fstderr,*) "varnce_rt_2(i,:) = ", pdf_params%varnce_rt_2(i,:)
       write(fstderr,*) "thl_1(i,:) = ", pdf_params%thl_1(i,:)
       write(fstderr,*) "thl_2(i,:) = ", pdf_params%thl_2(i,:)
       write(fstderr,*) "varnce_thl_1(i,:) = ", pdf_params%varnce_thl_1(i,:)
       write(fstderr,*) "varnce_thl_2(i,:) = ", pdf_params%varnce_thl_2(i,:)
       write(fstderr,*) "corr_w_rt_1(i,:) = ", pdf_params%corr_w_rt_1(i,:)
       write(fstderr,*) "corr_w_rt_2(i,:) = ", pdf_params%corr_w_rt_2(i,:)
       write(fstderr,*) "corr_w_thl_1(i,:) = ", pdf_params%corr_w_thl_1(i,:)
       write(fstderr,*) "corr_w_thl_2(i,:) = ", pdf_params%corr_w_thl_2(i,:)
       write(fstderr,*) "corr_rt_thl_1(i,:) = ", pdf_params%corr_rt_thl_1(i,:)
       write(fstderr,*) "corr_rt_thl_2(i,:) = ", pdf_params%corr_rt_thl_2(i,:)
       write(fstderr,*) "alpha_thl(i,:) = ", pdf_params%alpha_thl(i,:)
       write(fstderr,*) "alpha_rt(i,:) = ", pdf_params%alpha_rt(i,:)
       write(fstderr,*) "crt_1(i,:) = ", pdf_params%crt_1(i,:)
       write(fstderr,*) "crt_2(i,:) = ", pdf_params%crt_2(i,:)
       write(fstderr,*) "cthl_1(i,:) = ", pdf_params%cthl_1(i,:)
       write(fstderr,*) "cthl_2(i,:) = ", pdf_params%cthl_2(i,:)
       write(fstderr,*) "chi_1(i,:) = ", pdf_params%chi_1(i,:)
       write(fstderr,*) "chi_2(i,:) = ", pdf_params%chi_2(i,:)
       write(fstderr,*) "stdev_chi_1(i,:) = ", pdf_params%stdev_chi_1(i,:)
       write(fstderr,*) "stdev_chi_2(i,:) = ", pdf_params%stdev_chi_2(i,:)
       write(fstderr,*) "stdev_eta_1(i,:) = ", pdf_params%stdev_eta_1(i,:)
       write(fstderr,*) "stdev_eta_2(i,:) = ", pdf_params%stdev_eta_2(i,:)
       write(fstderr,*) "covar_chi_eta_1(i,:) = ", &
                        pdf_params%covar_chi_eta_1(i,:)
       write(fstderr,*) "covar_chi_eta_2(i,:) = ", &
                        pdf_params%covar_chi_eta_2(i,:)
       write(fstderr,*) "corr_w_chi_1(i,:) = ", pdf_params%corr_w_chi_1(i,:)
       write(fstderr,*) "corr_w_chi_2(i,:) = ", pdf_params%corr_w_chi_2(i,:)
       write(fstderr,*) "corr_w_eta_1(i,:) = ", pdf_params%corr_w_eta_1(i,:)
       write(fstderr,*) "corr_w_eta_2(i,:) = ", pdf_params%corr_w_eta_2(i,:)
       write(fstderr,*) "corr_chi_eta_1(i,:) = ", &
                        pdf_params%corr_chi_eta_1(i,:)
       write(fstderr,*) "corr_chi_eta_2(i,:) = ", &
                        pdf_params%corr_chi_eta_2(i,:)
       write(fstderr,*) "rsatl_1(i,:) = ", pdf_params%rsatl_1(i,:)
       write(fstderr,*) "rsatl_2(i,:) = ", pdf_params%rsatl_2(i,:)
       write(fstderr,*) "rc_1(i,:) = ", pdf_params%rc_1(i,:)
       write(fstderr,*) "rc_2(i,:) = ", pdf_params%rc_2(i,:)
       write(fstderr,*) "cloud_frac_1(i,:) = ", pdf_params%cloud_frac_1(i,:)
       write(fstderr,*) "cloud_frac_2(i,:) = ", pdf_params%cloud_frac_2(i,:)
       write(fstderr,*) "mixt_frac(i,:) = ", pdf_params%mixt_frac(i,:)
       write(fstderr,*) "ice_supersat_frac_1(i,:) = ", &
                        pdf_params%ice_supersat_frac_1(i,:)
       write(fstderr,*) "ice_supersat_frac_2(i,:) = ", &
                        pdf_params%ice_supersat_frac_2(i,:)

    enddo ! i = 1, ngrdcol


    return

  end subroutine print_pdf_params

  !=============================================================================

! The CLUBB_CAM preprocessor directives are being commented out because this
! code is now also used for WRF-CLUBB.
!#ifdef CLUBB_CAM /* Code for storing pdf_parameter structs in pbuf as array */

  subroutine pack_pdf_params(pdf_params, nz, &
                             r_param_array, &
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

    r_param_array(:,1) = pdf_params%w_1(1,k_start:k_end)
    r_param_array(:,2) = pdf_params%w_2(1,k_start:k_end)
    r_param_array(:,3) = pdf_params%varnce_w_1(1,k_start:k_end)
    r_param_array(:,4) = pdf_params%varnce_w_2(1,k_start:k_end)
    r_param_array(:,5) = pdf_params%rt_1(1,k_start:k_end)
    r_param_array(:,6) = pdf_params%rt_2(1,k_start:k_end)
    r_param_array(:,7) = pdf_params%varnce_rt_1(1,k_start:k_end)
    r_param_array(:,8) = pdf_params%varnce_rt_2(1,k_start:k_end)
    r_param_array(:,9) = pdf_params%thl_1(1,k_start:k_end)
    r_param_array(:,10) = pdf_params%thl_2(1,k_start:k_end)
    r_param_array(:,11) = pdf_params%varnce_thl_1(1,k_start:k_end)
    r_param_array(:,12) = pdf_params%varnce_thl_2(1,k_start:k_end)
    r_param_array(:,13) = pdf_params%corr_w_rt_1(1,k_start:k_end)
    r_param_array(:,14) = pdf_params%corr_w_rt_2(1,k_start:k_end)
    r_param_array(:,15) = pdf_params%corr_w_thl_1(1,k_start:k_end)
    r_param_array(:,16) = pdf_params%corr_w_thl_2(1,k_start:k_end)
    r_param_array(:,17) = pdf_params%corr_rt_thl_1(1,k_start:k_end)
    r_param_array(:,18) = pdf_params%corr_rt_thl_2(1,k_start:k_end)
    r_param_array(:,19) = pdf_params%alpha_thl(1,k_start:k_end)
    r_param_array(:,20) = pdf_params%alpha_rt(1,k_start:k_end)
    r_param_array(:,21) = pdf_params%crt_1(1,k_start:k_end)
    r_param_array(:,22) = pdf_params%crt_2(1,k_start:k_end)
    r_param_array(:,23) = pdf_params%cthl_1(1,k_start:k_end)
    r_param_array(:,24) = pdf_params%cthl_2(1,k_start:k_end)
    r_param_array(:,25) = pdf_params%chi_1(1,k_start:k_end)
    r_param_array(:,26) = pdf_params%chi_2(1,k_start:k_end)
    r_param_array(:,27) = pdf_params%stdev_chi_1(1,k_start:k_end)
    r_param_array(:,28) = pdf_params%stdev_chi_2(1,k_start:k_end)
    r_param_array(:,29) = pdf_params%stdev_eta_1(1,k_start:k_end)
    r_param_array(:,30) = pdf_params%stdev_eta_2(1,k_start:k_end)
    r_param_array(:,31) = pdf_params%covar_chi_eta_1(1,k_start:k_end)
    r_param_array(:,32) = pdf_params%covar_chi_eta_2(1,k_start:k_end)
    r_param_array(:,33) = pdf_params%corr_w_chi_1(1,k_start:k_end) 
    r_param_array(:,34) = pdf_params%corr_w_chi_2(1,k_start:k_end) 
    r_param_array(:,35) = pdf_params%corr_w_eta_1(1,k_start:k_end) 
    r_param_array(:,36) = pdf_params%corr_w_eta_2(1,k_start:k_end) 
    r_param_array(:,37) = pdf_params%corr_chi_eta_1(1,k_start:k_end) 
    r_param_array(:,38) = pdf_params%corr_chi_eta_2(1,k_start:k_end) 
    r_param_array(:,39) = pdf_params%rsatl_1(1,k_start:k_end)
    r_param_array(:,40) = pdf_params%rsatl_2(1,k_start:k_end)
    r_param_array(:,41) = pdf_params%rc_1(1,k_start:k_end)
    r_param_array(:,42) = pdf_params%rc_2(1,k_start:k_end)
    r_param_array(:,43) = pdf_params%cloud_frac_1(1,k_start:k_end)
    r_param_array(:,44) = pdf_params%cloud_frac_2(1,k_start:k_end)
    r_param_array(:,45) = pdf_params%mixt_frac(1,k_start:k_end)
    r_param_array(:,46) = pdf_params%ice_supersat_frac_1(1,k_start:k_end)
    r_param_array(:,47) = pdf_params%ice_supersat_frac_2(1,k_start:k_end)

  end subroutine pack_pdf_params
!===============================================================!
  subroutine unpack_pdf_params(r_param_array, nz, &
                               pdf_params, &
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
    
    pdf_params%w_1(1,k_start:k_end) = r_param_array(:,1)
    pdf_params%w_2(1,k_start:k_end) = r_param_array(:,2)
    pdf_params%varnce_w_1(1,k_start:k_end) = r_param_array(:,3)
    pdf_params%varnce_w_2(1,k_start:k_end) = r_param_array(:,4)
    pdf_params%rt_1(1,k_start:k_end) = r_param_array(:,5)
    pdf_params%rt_2(1,k_start:k_end) = r_param_array(:,6)
    pdf_params%varnce_rt_1(1,k_start:k_end) = r_param_array(:,7)
    pdf_params%varnce_rt_2(1,k_start:k_end)= r_param_array(:,8)
    pdf_params%thl_1(1,k_start:k_end) = r_param_array(:,9)
    pdf_params%thl_2(1,k_start:k_end) = r_param_array(:,10)
    pdf_params%varnce_thl_1(1,k_start:k_end) = r_param_array(:,11)
    pdf_params%varnce_thl_2(1,k_start:k_end) = r_param_array(:,12)
    pdf_params%corr_w_rt_1(1,k_start:k_end) = r_param_array(:,13)
    pdf_params%corr_w_rt_2(1,k_start:k_end) = r_param_array(:,14)
    pdf_params%corr_w_thl_1(1,k_start:k_end) = r_param_array(:,15)
    pdf_params%corr_w_thl_2(1,k_start:k_end) = r_param_array(:,16)
    pdf_params%corr_rt_thl_1(1,k_start:k_end) = r_param_array(:,17)
    pdf_params%corr_rt_thl_2(1,k_start:k_end) = r_param_array(:,18)
    pdf_params%alpha_thl(1,k_start:k_end) = r_param_array(:,19)
    pdf_params%alpha_rt(1,k_start:k_end) = r_param_array(:,20)
    pdf_params%crt_1(1,k_start:k_end) = r_param_array(:,21)
    pdf_params%crt_2(1,k_start:k_end) = r_param_array(:,22)
    pdf_params%cthl_1(1,k_start:k_end) = r_param_array(:,23)
    pdf_params%cthl_2(1,k_start:k_end) = r_param_array(:,24)
    pdf_params%chi_1(1,k_start:k_end) = r_param_array(:,25)
    pdf_params%chi_2(1,k_start:k_end) = r_param_array(:,26)
    pdf_params%stdev_chi_1(1,k_start:k_end) = r_param_array(:,27)
    pdf_params%stdev_chi_2(1,k_start:k_end) = r_param_array(:,28)
    pdf_params%stdev_eta_1(1,k_start:k_end) = r_param_array(:,29)
    pdf_params%stdev_eta_2(1,k_start:k_end) = r_param_array(:,30)
    pdf_params%covar_chi_eta_1(1,k_start:k_end) = r_param_array(:,31)
    pdf_params%covar_chi_eta_2(1,k_start:k_end) = r_param_array(:,32)
    pdf_params%corr_w_chi_1(1,k_start:k_end) = r_param_array(:,33)
    pdf_params%corr_w_chi_2(1,k_start:k_end) = r_param_array(:,34)
    pdf_params%corr_w_eta_1(1,k_start:k_end) = r_param_array(:,35)
    pdf_params%corr_w_eta_2(1,k_start:k_end) = r_param_array(:,36)
    pdf_params%corr_chi_eta_1(1,k_start:k_end) = r_param_array(:,37)
    pdf_params%corr_chi_eta_2(1,k_start:k_end) = r_param_array(:,38)
    pdf_params%rsatl_1(1,k_start:k_end) = r_param_array(:,39)
    pdf_params%rsatl_2(1,k_start:k_end) = r_param_array(:,40)
    pdf_params%rc_1(1,k_start:k_end) = r_param_array(:,41)
    pdf_params%rc_2(1,k_start:k_end) = r_param_array(:,42)
    pdf_params%cloud_frac_1(1,k_start:k_end) = r_param_array(:,43)
    pdf_params%cloud_frac_2(1,k_start:k_end) = r_param_array(:,44)
    pdf_params%mixt_frac(1,k_start:k_end) = r_param_array(:,45)
    pdf_params%ice_supersat_frac_1(1,k_start:k_end) = r_param_array(:,46)
    pdf_params%ice_supersat_frac_2(1,k_start:k_end) = r_param_array(:,47)

  end subroutine unpack_pdf_params
  
  
  !================================================================================================
  ! copy_single_pdf_params_to_multi - copies values of a single column version of pdf_params
  !   to a multiple column version for a specified column.
  !
  ! NOTE: THIS SUBROUTINE IS INTENDED TO BE TEMPORARY AND SHOULD BECOME UNNECESSARY ONCE 
  !       CLUBB IS ABLE TO OPERATE OVER MULTIPLE COLUMNS.
  !       See https://github.com/larson-group/cam/issues/129#issuecomment-827944454
  !================================================================================================
  subroutine copy_single_pdf_params_to_multi( pdf_params_single, icol, &
                                              pdf_params_multi )
    
    implicit none

    ! Input Variable(s)
    integer, intent(in) :: &
      icol   ! Column number to copy to
      
    type(pdf_parameter), intent(in) :: &
      pdf_params_single  ! PDF parameters            [units vary]

    ! Output Variable(s)
    type(pdf_parameter), intent(inout) :: &
      pdf_params_multi  ! PDF parameters            [units vary]
      
    pdf_params_multi%w_1(icol,:)                  = pdf_params_single%w_1(1,:) 
    pdf_params_multi%w_2(icol,:)                  = pdf_params_single%w_2(1,:) 
    pdf_params_multi%varnce_w_1(icol,:)           = pdf_params_single%varnce_w_1(1,:) 
    pdf_params_multi%varnce_w_2(icol,:)           = pdf_params_single%varnce_w_2(1,:) 
    pdf_params_multi%rt_1(icol,:)                 = pdf_params_single%rt_1(1,:) 
    pdf_params_multi%rt_2(icol,:)                 = pdf_params_single%rt_2(1,:) 
    pdf_params_multi%varnce_rt_1(icol,:)          = pdf_params_single%varnce_rt_1(1,:) 
    pdf_params_multi%varnce_rt_2(icol,:)          = pdf_params_single%varnce_rt_2(1,:) 
    pdf_params_multi%thl_1(icol,:)                = pdf_params_single%thl_1(1,:) 
    pdf_params_multi%thl_2(icol,:)                = pdf_params_single%thl_2(1,:) 
    pdf_params_multi%varnce_thl_1(icol,:)         = pdf_params_single%varnce_thl_1(1,:) 
    pdf_params_multi%varnce_thl_2(icol,:)         = pdf_params_single%varnce_thl_2(1,:) 
    pdf_params_multi%corr_w_rt_1(icol,:)          = pdf_params_single%corr_w_rt_1(1,:) 
    pdf_params_multi%corr_w_rt_2(icol,:)          = pdf_params_single%corr_w_rt_2(1,:) 
    pdf_params_multi%corr_w_thl_1(icol,:)         = pdf_params_single%corr_w_thl_1(1,:) 
    pdf_params_multi%corr_w_thl_2(icol,:)         = pdf_params_single%corr_w_thl_2(1,:) 
    pdf_params_multi%corr_rt_thl_1(icol,:)        = pdf_params_single%corr_rt_thl_1(1,:) 
    pdf_params_multi%corr_rt_thl_2(icol,:)        = pdf_params_single%corr_rt_thl_2(1,:) 
    pdf_params_multi%alpha_thl(icol,:)            = pdf_params_single%alpha_thl(1,:) 
    pdf_params_multi%alpha_rt(icol,:)             = pdf_params_single%alpha_rt(1,:) 
    pdf_params_multi%crt_1(icol,:)                = pdf_params_single%crt_1(1,:) 
    pdf_params_multi%crt_2(icol,:)                = pdf_params_single%crt_2(1,:) 
    pdf_params_multi%cthl_1(icol,:)               = pdf_params_single%cthl_1(1,:) 
    pdf_params_multi%cthl_2(icol,:)               = pdf_params_single%cthl_2(1,:) 
    pdf_params_multi%chi_1(icol,:)                = pdf_params_single%chi_1(1,:) 
    pdf_params_multi%chi_2(icol,:)                = pdf_params_single%chi_2(1,:) 
    pdf_params_multi%stdev_chi_1(icol,:)          = pdf_params_single%stdev_chi_1(1,:) 
    pdf_params_multi%stdev_chi_2(icol,:)          = pdf_params_single%stdev_chi_2(1,:) 
    pdf_params_multi%stdev_eta_1(icol,:)          = pdf_params_single%stdev_eta_1(1,:) 
    pdf_params_multi%stdev_eta_2(icol,:)          = pdf_params_single%stdev_eta_2(1,:) 
    pdf_params_multi%covar_chi_eta_1(icol,:)      = pdf_params_single%covar_chi_eta_1(1,:) 
    pdf_params_multi%covar_chi_eta_2(icol,:)      = pdf_params_single%covar_chi_eta_2(1,:) 
    pdf_params_multi%corr_w_chi_1(icol,:)         = pdf_params_single%corr_w_chi_1(1,:) 
    pdf_params_multi%corr_w_chi_2(icol,:)         = pdf_params_single%corr_w_chi_2(1,:) 
    pdf_params_multi%corr_w_eta_1(icol,:)         = pdf_params_single%corr_w_eta_1(1,:) 
    pdf_params_multi%corr_w_eta_2(icol,:)         = pdf_params_single%corr_w_eta_2(1,:) 
    pdf_params_multi%corr_chi_eta_1(icol,:)       = pdf_params_single%corr_chi_eta_1(1,:) 
    pdf_params_multi%corr_chi_eta_2(icol,:)       = pdf_params_single%corr_chi_eta_2(1,:) 
    pdf_params_multi%rsatl_1(icol,:)              = pdf_params_single%rsatl_1(1,:) 
    pdf_params_multi%rsatl_2(icol,:)              = pdf_params_single%rsatl_2(1,:) 
    pdf_params_multi%rc_1(icol,:)                 = pdf_params_single%rc_1(1,:) 
    pdf_params_multi%rc_2(icol,:)                 = pdf_params_single%rc_2(1,:) 
    pdf_params_multi%cloud_frac_1(icol,:)         = pdf_params_single%cloud_frac_1(1,:) 
    pdf_params_multi%cloud_frac_2(icol,:)         = pdf_params_single%cloud_frac_2(1,:) 
    pdf_params_multi%mixt_frac(icol,:)            = pdf_params_single%mixt_frac(1,:) 
    pdf_params_multi%ice_supersat_frac_1(icol,:)  = pdf_params_single%ice_supersat_frac_1(1,:) 
    pdf_params_multi%ice_supersat_frac_2(icol,:)  = pdf_params_single%ice_supersat_frac_2(1,:) 
    
  end subroutine copy_single_pdf_params_to_multi
  
  !================================================================================================
  ! copy_multi_pdf_params_to_single - copies values of a multiple column version of pdf_params
  !   at a specified column to a single column version.
  !
  ! NOTE: THIS SUBROUTINE IS INTENDED TO BE TEMPORARY AND SHOULD BECOME UNNECESSARY ONCE 
  !       CLUBB IS ABLE TO OPERATE OVER MULTIPLE COLUMNS.
  !       See https://github.com/larson-group/cam/issues/129#issuecomment-827944454
  !================================================================================================
  subroutine copy_multi_pdf_params_to_single( pdf_params_multi, icol, &
                                              pdf_params_single )
    
    implicit none

    ! Input Variable(s)
    integer, intent(in) :: &
      icol   ! Column number to copy to
      
    type(pdf_parameter), intent(in) :: &
      pdf_params_multi  ! PDF parameters            [units vary]

    ! Output Variable(s)
    type(pdf_parameter), intent(inout) :: &
      pdf_params_single   ! PDF parameters            [units vary]
      
    pdf_params_single%w_1(1,:)                  = pdf_params_multi%w_1(icol,:) 
    pdf_params_single%w_2(1,:)                  = pdf_params_multi%w_2(icol,:) 
    pdf_params_single%varnce_w_1(1,:)           = pdf_params_multi%varnce_w_1(icol,:) 
    pdf_params_single%varnce_w_2(1,:)           = pdf_params_multi%varnce_w_2(icol,:) 
    pdf_params_single%rt_1(1,:)                 = pdf_params_multi%rt_1(icol,:) 
    pdf_params_single%rt_2(1,:)                 = pdf_params_multi%rt_2(icol,:) 
    pdf_params_single%varnce_rt_1(1,:)          = pdf_params_multi%varnce_rt_1(icol,:) 
    pdf_params_single%varnce_rt_2(1,:)          = pdf_params_multi%varnce_rt_2(icol,:) 
    pdf_params_single%thl_1(1,:)                = pdf_params_multi%thl_1(icol,:) 
    pdf_params_single%thl_2(1,:)                = pdf_params_multi%thl_2(icol,:) 
    pdf_params_single%varnce_thl_1(1,:)         = pdf_params_multi%varnce_thl_1(icol,:) 
    pdf_params_single%varnce_thl_2(1,:)         = pdf_params_multi%varnce_thl_2(icol,:) 
    pdf_params_single%corr_w_rt_1(1,:)          = pdf_params_multi%corr_w_rt_1(icol,:) 
    pdf_params_single%corr_w_rt_2(1,:)          = pdf_params_multi%corr_w_rt_2(icol,:) 
    pdf_params_single%corr_w_thl_1(1,:)         = pdf_params_multi%corr_w_thl_1(icol,:) 
    pdf_params_single%corr_w_thl_2(1,:)         = pdf_params_multi%corr_w_thl_2(icol,:) 
    pdf_params_single%corr_rt_thl_1(1,:)        = pdf_params_multi%corr_rt_thl_1(icol,:) 
    pdf_params_single%corr_rt_thl_2(1,:)        = pdf_params_multi%corr_rt_thl_2(icol,:) 
    pdf_params_single%alpha_thl(1,:)            = pdf_params_multi%alpha_thl(icol,:) 
    pdf_params_single%alpha_rt(1,:)             = pdf_params_multi%alpha_rt(icol,:) 
    pdf_params_single%crt_1(1,:)                = pdf_params_multi%crt_1(icol,:) 
    pdf_params_single%crt_2(1,:)                = pdf_params_multi%crt_2(icol,:) 
    pdf_params_single%cthl_1(1,:)               = pdf_params_multi%cthl_1(icol,:) 
    pdf_params_single%cthl_2(1,:)               = pdf_params_multi%cthl_2(icol,:) 
    pdf_params_single%chi_1(1,:)                = pdf_params_multi%chi_1(icol,:) 
    pdf_params_single%chi_2(1,:)                = pdf_params_multi%chi_2(icol,:) 
    pdf_params_single%stdev_chi_1(1,:)          = pdf_params_multi%stdev_chi_1(icol,:) 
    pdf_params_single%stdev_chi_2(1,:)          = pdf_params_multi%stdev_chi_2(icol,:) 
    pdf_params_single%stdev_eta_1(1,:)          = pdf_params_multi%stdev_eta_1(icol,:) 
    pdf_params_single%stdev_eta_2(1,:)          = pdf_params_multi%stdev_eta_2(icol,:) 
    pdf_params_single%covar_chi_eta_1(1,:)      = pdf_params_multi%covar_chi_eta_1(icol,:) 
    pdf_params_single%covar_chi_eta_2(1,:)      = pdf_params_multi%covar_chi_eta_2(icol,:) 
    pdf_params_single%corr_w_chi_1(1,:)         = pdf_params_multi%corr_w_chi_1(icol,:) 
    pdf_params_single%corr_w_chi_2(1,:)         = pdf_params_multi%corr_w_chi_2(icol,:) 
    pdf_params_single%corr_w_eta_1(1,:)         = pdf_params_multi%corr_w_eta_1(icol,:) 
    pdf_params_single%corr_w_eta_2(1,:)         = pdf_params_multi%corr_w_eta_2(icol,:) 
    pdf_params_single%corr_chi_eta_1(1,:)       = pdf_params_multi%corr_chi_eta_1(icol,:) 
    pdf_params_single%corr_chi_eta_2(1,:)       = pdf_params_multi%corr_chi_eta_2(icol,:) 
    pdf_params_single%rsatl_1(1,:)              = pdf_params_multi%rsatl_1(icol,:) 
    pdf_params_single%rsatl_2(1,:)              = pdf_params_multi%rsatl_2(icol,:) 
    pdf_params_single%rc_1(1,:)                 = pdf_params_multi%rc_1(icol,:) 
    pdf_params_single%rc_2(1,:)                 = pdf_params_multi%rc_2(icol,:) 
    pdf_params_single%cloud_frac_1(1,:)         = pdf_params_multi%cloud_frac_1(icol,:) 
    pdf_params_single%cloud_frac_2(1,:)         = pdf_params_multi%cloud_frac_2(icol,:) 
    pdf_params_single%mixt_frac(1,:)            = pdf_params_multi%mixt_frac(icol,:) 
    pdf_params_single%ice_supersat_frac_1(1,:)  = pdf_params_multi%ice_supersat_frac_1(icol,:) 
    pdf_params_single%ice_supersat_frac_2(1,:)  = pdf_params_multi%ice_supersat_frac_2(icol,:) 
    
  end subroutine copy_multi_pdf_params_to_single

!#endif /* CLUBB_CAM */

end module pdf_parameter_module
