!---------------------------------------------------------------------------
! $Id$
!===============================================================================
module pdf_closure_module

  ! Options for the two component normal (double Gaussian) PDF type to use for
  ! the w, rt, and theta-l (or w, chi, and eta) portion of CLUBB's multivariate,
  ! two-component PDF.
  use model_flags, only: &
      iiPDF_ADG1,       & ! ADG1 PDF
      iiPDF_ADG2,       & ! ADG2 PDF
      iiPDF_3D_Luhar,   & ! 3D Luhar PDF
      iiPDF_new,        & ! new PDF
      iiPDF_TSDADG,     & ! new TSDADG PDF
      iiPDF_LY93,       & ! Lewellen and Yoh (1993)
      iiPDF_new_hybrid    ! new hybrid PDF

  implicit none

  public :: pdf_closure, &
            calc_wp4_pdf, &
            calc_wp2xp_pdf, &
            calc_wpxp2_pdf, &
            calc_wpxpyp_pdf, &
            calc_vert_avg_cf_component, &
            calc_w_up_in_cloud

  private ! Set Default Scope

  contains
!------------------------------------------------------------------------

  !#######################################################################
  !#######################################################################
  ! If you change the argument list of pdf_closure you also have to
  ! change the calls to this function in the host models CAM, WRF, SAM
  ! and GFDL.
  !#######################################################################
  !#######################################################################
  subroutine pdf_closure( gr, nz, ngrdcol,                          &
                          hydromet_dim, p_in_Pa, exner, thv_ds,     &
                          wm, wp2, wp3, sigma_sqd_w,                &
                          Skw, Skthl_in, Skrt_in, Sku_in, Skv_in,   &
                          rtm, rtp2, wprtp,                         &
                          thlm, thlp2, wpthlp,                      &
                          um, up2, upwp,                            &
                          vm, vp2, vpwp,                            &
                          rtpthlp,                                  &
                          sclrm, wpsclrp, sclrp2,                   &
                          sclrprtp, sclrpthlp, Sksclr_in,           &
                          gamma_Skw_fnc,                            &
#ifdef GFDL
                          RH_crit, do_liquid_only_in_clubb,         & ! h1g, 2010-06-15
#endif
                          wphydrometp, wp2hmp,                      &
                          rtphmp, thlphmp,                          &
                          clubb_params,                             &
                          iiPDF_type,                               &
                          wpup2, wpvp2,                             &
                          wp2up2, wp2vp2, wp4,                      &
                          wprtp2, wp2rtp,                           &
                          wpthlp2, wp2thlp, wprtpthlp,              &
                          cloud_frac, ice_supersat_frac,            &
                          rcm, wpthvp, wp2thvp, rtpthvp,            &
                          thlpthvp, wprcp, wp2rcp, rtprcp,          &
                          thlprcp, rcp2,                            &
                          uprcp, vprcp, w_up_in_cloud,              &
                          pdf_params, pdf_implicit_coefs_terms,     &
                          F_w, F_rt, F_thl,                         &
                          min_F_w, max_F_w,                         &
                          min_F_rt, max_F_rt,                       &
                          min_F_thl, max_F_thl,                     &
                          wpsclrprtp, wpsclrp2, sclrpthvp,          &
                          wpsclrpthlp, sclrprcp, wp2sclrp,          &
                          rc_coef                                   )


    ! Description:
    ! Subroutine that computes pdf parameters analytically.
    !
    ! Based of the original formulation, but with some tweaks
    ! to remove some of the less realistic assumptions and
    ! improve transport terms.

    !   Corrected version that should remove inconsistency

    ! References:
    !   The shape of CLUBB's PDF is given by the expression in
    !   https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:clubb_pdf

    !   Eqn. 29, 30, 31, 32 & 33  on p. 3547 of
    !   ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    !   Method and Model Description'' Golaz, et al. (2002)
    !   JAS, Vol. 59, pp. 3540--3551.
    !----------------------------------------------------------------------

    use grid_class, only: &
        grid ! Type

    use constants_clubb, only: &  ! Constants
        three,          & ! 3
        one,            & ! 1
        one_half,       & ! 1/2
        zero,           & ! 0
        Cp,             & ! Dry air specific heat at constant p [J/kg/K]
        Lv,             & ! Latent heat of vaporization         [J/kg]
        ep1,            & ! (1.0-ep)/ep; ep1 = 0.61             [-]
        ep2,            & ! 1.0/ep;      ep2 = 1.61             [-]
        rt_tol,         & ! Tolerance for r_t                   [kg/kg]
        thl_tol,        & ! Tolerance for th_l                  [K]
        T_freeze_K,     & ! Freezing point of water             [K]
        fstderr,        &
        zero_threshold, &
        eps, &
        w_tol

    use parameters_model, only: &
        mixt_frac_max_mag, & ! Variable(s)
        sclr_dim             ! Number of passive scalar variables

    use parameter_indices, only: &
        nparams,                       & ! Variable(s)
        ibeta,                         &
        iSkw_denom_coef,               &
        islope_coef_spread_DG_means_w, &
        ipdf_component_stdev_factor_w, &
        icoef_spread_DG_means_rt,      &
        icoef_spread_DG_means_thl

    use pdf_parameter_module, only:  &
        pdf_parameter,        & ! Variable Type
        implicit_coefs_terms

    use new_pdf_main, only: &
        new_pdf_driver    ! Procedure(s)

    use new_hybrid_pdf_main, only: &
        new_hybrid_pdf_driver    ! Procedure(s)

    use adg1_adg2_3d_luhar_pdf, only: &
        ADG1_pdf_driver,     & ! Procedure(s)
        ADG2_pdf_driver,     &
        Luhar_3D_pdf_driver

    use new_tsdadg_pdf, only: &
        tsdadg_pdf_driver    ! Procedure(s)

    use LY93_pdf, only: &
        LY93_driver    ! Procedure(s)

    use pdf_utilities, only: &
        calc_comp_corrs_binormal, & ! Procedure(s)
        calc_corr_chi_x,          &
        calc_corr_eta_x

    use array_index, only: &
        l_mix_rat_hm  ! Variable(s)

    use model_flags, only: &
        l_explicit_turbulent_adv_xpyp ! Variable(s)

    use numerical_check, only:  & 
        pdf_closure_check ! Procedure(s)

    use saturation, only:  & 
        sat_mixrat_liq, & ! Procedure(s)
        sat_mixrat_ice

    use stats_variables, only: &
        ircp2,      & ! Variables
        iwprtp2,    &
        iwprtpthlp, &
        iwpthlp2,   &
        iw_up_in_cloud

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_fatal_error              ! Constant

    implicit none

    intrinsic :: sqrt, exp, min, max, abs, present

    ! Input Variables
    integer, intent(in) :: &
      hydromet_dim, & ! Number of hydrometeor species              [#]
      nz, &
      ngrdcol
      
    type (grid), target, dimension(ngrdcol), intent(in) :: &
      gr

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  & 
      p_in_Pa,     & ! Pressure                                   [Pa]
      exner,       & ! Exner function                             [-]
      thv_ds,      & ! Dry, base-state theta_v (ref. th_l here)   [K]
      wm,          & ! mean w-wind component (vertical velocity)  [m/s] 
      wp2,         & ! w'^2                                       [m^2/s^2] 
      wp3,         & ! w'^3                                       [m^3/s^3]
      Skw,         & ! Skewness of w                              [-]
      Skthl_in,    & ! Skewness of thl                            [-]
      Skrt_in,     & ! Skewness of rt                             [-]
      Sku_in,      & ! Skewness of u                              [-]
      Skv_in,      & ! Skewness of v                              [-]
      rtm,         & ! Mean total water mixing ratio              [kg/kg]
      rtp2,        & ! r_t'^2                                     [(kg/kg)^2]
      wprtp,       & ! w'r_t'                                     [(kg/kg)(m/s)]
      thlm,        & ! Mean liquid water potential temperature    [K]
      thlp2,       & ! th_l'^2                                    [K^2]
      wpthlp,      & ! w'th_l'                                    [K(m/s)]
      rtpthlp        ! r_t'th_l'                                  [K(kg/kg)]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
      um,          & ! Grid-mean eastward wind     [m/s]
      up2,         & ! u'^2                        [(m/s)^2]
      upwp,        & ! u'w'                        [(m/s)^2]
      vm,          & ! Grid-mean northward wind    [m/s]
      vp2,         & ! v'^2                        [(m/s)^2]
      vpwp           ! v'w'                        [(m/s)^2]

    real( kind = core_rknd ), dimension(ngrdcol,nz, sclr_dim), intent(in) ::  & 
      sclrm,       & ! Mean passive scalar        [units vary]
      wpsclrp,     & ! w' sclr'                   [units vary]
      sclrp2,      & ! sclr'^2                    [units vary]
      sclrprtp,    & ! sclr' r_t'                 [units vary]
      sclrpthlp,   & ! sclr' th_l'                [units vary]
      Sksclr_in      ! Skewness of sclr           [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
      gamma_Skw_fnc    ! Gamma as a function of skewness            [-]

#ifdef  GFDL
    ! critial relative humidity for nucleation
    real( kind = core_rknd ), dimension(ngrdcol, nz, min(1,sclr_dim), 2 ), intent(in) ::  & ! h1g, 2010-06-15
       RH_crit     ! critical relative humidity for droplet and ice nucleation
! ---> h1g, 2012-06-14
    logical, intent(in)                 ::  do_liquid_only_in_clubb
! <--- h1g, 2012-06-14
#endif

    real( kind = core_rknd ), dimension(ngrdcol,nz, hydromet_dim), intent(in) :: &
      wphydrometp, & ! Covariance of w and a hydrometeor    [(m/s) <hm units>]
      wp2hmp,      & ! Third-order moment:  < w'^2 hm' >    [(m/s)^2 <hm units>]
      rtphmp,      & ! Covariance of rt and a hydrometeor   [(kg/kg) <hm units>]
      thlphmp        ! Covariance of thl and a hydrometeor  [K <hm units>]

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    integer, intent(in) :: &
      iiPDF_type    ! Selected option for the two-component normal (double
                    ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                    ! w, chi, and eta) portion of CLUBB's multivariate,
                    ! two-component PDF.

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(inout) :: &
      ! If iiPDF_type == iiPDF_ADG2, this gets overwritten. Therefore,
      ! intent(inout). Otherwise it should be intent(in)
      sigma_sqd_w   ! Width of individual w plumes               [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) ::  & 
      wpup2,              & ! w'u'^2                [m^3/s^3]
      wpvp2,              & ! w'v'^2                [m^3/s^3]
      wp2up2,             & ! w'^2u'^2              [m^2/s^4]
      wp2vp2,             & ! w'^2v'^2              [m^2/s^4]
      wp4,                & ! w'^4                  [m^4/s^4]
      wprtp2,             & ! w' r_t'               [(m kg)/(s kg)]
      wp2rtp,             & ! w'^2 r_t'             [(m^2 kg)/(s^2 kg)]
      wpthlp2,            & ! w' th_l'^2            [(m K^2)/s]
      wp2thlp,            & ! w'^2 th_l'            [(m^2 K)/s^2]
      cloud_frac,         & ! Cloud fraction        [-]
      ice_supersat_frac,  & ! Ice cloud fracion     [-]
      rcm,                & ! Mean liquid water     [kg/kg]
      wpthvp,             & ! Buoyancy flux         [(K m)/s] 
      wp2thvp,            & ! w'^2 th_v'            [(m^2 K)/s^2]
      rtpthvp,            & ! r_t' th_v'            [(kg K)/kg]
      thlpthvp,           & ! th_l' th_v'           [K^2]
      wprcp,              & ! w' r_c'               [(m kg)/(s kg)]
      wp2rcp,             & ! w'^2 r_c'             [(m^2 kg)/(s^2 kg)]
      rtprcp,             & ! r_t' r_c'             [(kg^2)/(kg^2)]
      thlprcp,            & ! th_l' r_c'            [(K kg)/kg]
      rcp2,               & ! r_c'^2                [(kg^2)/(kg^2)]
      wprtpthlp,          & ! w' r_t' th_l'         [(m kg K)/(s kg)]
      w_up_in_cloud         ! wm over Clouds        [m/s] TODO Correct place?

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) ::  &
      uprcp,              & ! u' r_c'               [(m kg)/(s kg)]
      vprcp                 ! v' r_c'               [(m kg)/(s kg)]

    type(pdf_parameter), intent(inout) :: & 
      pdf_params     ! pdf paramters         [units vary]

    type(implicit_coefs_terms), dimension(ngrdcol), intent(inout) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    ! Parameters output only for recording statistics (new PDF).
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      F_w,   & ! Parameter for the spread of the PDF component means of w    [-]
      F_rt,  & ! Parameter for the spread of the PDF component means of rt   [-]
      F_thl    ! Parameter for the spread of the PDF component means of thl  [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      min_F_w,   & ! Minimum allowable value of parameter F_w      [-]
      max_F_w,   & ! Maximum allowable value of parameter F_w      [-]
      min_F_rt,  & ! Minimum allowable value of parameter F_rt     [-]
      max_F_rt,  & ! Maximum allowable value of parameter F_rt     [-]
      min_F_thl, & ! Minimum allowable value of parameter F_thl    [-]
      max_F_thl    ! Maximum allowable value of parameter F_thl    [-]

    ! Output (passive scalar variables)

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,nz, sclr_dim) ::  & 
      sclrpthvp, & 
      sclrprcp, & 
      wpsclrp2, & 
      wpsclrprtp, & 
      wpsclrpthlp, & 
      wp2sclrp

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      rc_coef    ! Coefficient on X'r_c' in X'th_v' equation    [K/(kg/kg)]

    ! Local Variables

    ! Variables that are stored in derived data type pdf_params.
    real( kind = core_rknd ), dimension(ngrdcol,nz) ::  &
      u_1,           & ! Mean of eastward wind (1st PDF component)         [m/s]
      u_2,           & ! Mean of eastward wind (2nd PDF component)         [m/s]
      varnce_u_1,    & ! Variance of u (1st PDF component)             [m^2/s^2]
      varnce_u_2,    & ! Variance of u (2nd PDF component)             [m^2/s^2]
      v_1,           & ! Mean of northward wind (1st PDF component)        [m/s]
      v_2,           & ! Mean of northward wind (2nd PDF component)        [m/s]
      varnce_v_1,    & ! Variance of v (1st PDF component)             [m^2/s^2]
      varnce_v_2,    & ! Variance of v (2nd PDF component)             [m^2/s^2]
      alpha_u,       & ! Factor relating to normalized variance for u        [-]
      alpha_v          ! Factor relating to normalized variance for v        [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      corr_u_w_1,      & ! Correlation of u and w   (1st PDF component)      [-]
      corr_u_w_2,      & ! Correlation of u and w   (2nd PDF component)      [-]
      corr_v_w_1,      & ! Correlation of v and w   (1st PDF component)      [-]
      corr_v_w_2         ! Correlation of v and w   (2nd PDF component)      [-]

    ! Note:  alpha coefficients = 0.5 * ( 1 - correlations^2 ).
    !        These are used to calculate the scalar widths
    !        varnce_thl_1, varnce_thl_2, varnce_rt_1, and varnce_rt_2 as in
    !        Eq. (34) of Larson and Golaz (2005)

    ! Passive scalar local variables

    real( kind = core_rknd ), dimension(ngrdcol,nz,sclr_dim) ::  & 
      sclr1, sclr2,  &
      varnce_sclr1, varnce_sclr2, & 
      alpha_sclr,  & 
      corr_sclr_thl_1, corr_sclr_thl_2, &
      corr_sclr_rt_1, corr_sclr_rt_2, &
      corr_w_sclr_1, corr_w_sclr_2

    logical :: &
      l_scalar_calc, &  ! True if sclr_dim > 0
      l_calc_ice_supersat_frac ! True if we should calculate ice_supersat_frac

    ! Quantities needed to predict higher order moments
    real( kind = core_rknd ), dimension(ngrdcol,nz) ::  & 
      tl1, tl2

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      sqrt_wp2, & ! Square root of wp2          [m/s]
      Skthl,    & ! Skewness of thl             [-]
      Skrt,     & ! Skewness of rt              [-]
      Sku,      & ! Skewness of u               [-]
      Skv         ! Skewness of v               [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz,sclr_dim) :: &
      Sksclr      ! Skewness of rt              [-]

    ! Thermodynamic quantity

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      wprcp_contrib_comp_1,   & ! <w'rc'> contrib. (1st PDF comp.)  [m/s(kg/kg)]
      wprcp_contrib_comp_2,   & ! <w'rc'> contrib. (2nd PDF comp.)  [m/s(kg/kg)]
      wp2rcp_contrib_comp_1,  & ! <w'^2rc'> contrib. (1st comp) [m^2/s^2(kg/kg)]
      wp2rcp_contrib_comp_2,  & ! <w'^2rc'> contrib. (2nd comp) [m^2/s^2(kg/kg)]
      rtprcp_contrib_comp_1,  & ! <rt'rc'> contrib. (1st PDF comp.)  [kg^2/kg^2]
      rtprcp_contrib_comp_2,  & ! <rt'rc'> contrib. (2nd PDF comp.)  [kg^2/kg^2]
      thlprcp_contrib_comp_1, & ! <thl'rc'> contrib. (1st PDF comp.)  [K(kg/kg)]
      thlprcp_contrib_comp_2, & ! <thl'rc'> contrib. (2nd PDF comp.)  [K(kg/kg)]
      uprcp_contrib_comp_1,   & ! <u'rc'> contrib. (1st PDF comp.)  [m/s(kg/kg)]
      uprcp_contrib_comp_2,   & ! <u'rc'> contrib. (2nd PDF comp.)  [m/s(kg/kg)]
      vprcp_contrib_comp_1,   & ! <v'rc'> contrib. (1st PDF comp.)  [m/s(kg/kg)]
      vprcp_contrib_comp_2      ! <v'rc'> contrib. (2nd PDF comp.)  [m/s(kg/kg)]

    ! variables for computing ice cloud fraction
    real( kind = core_rknd), dimension(ngrdcol,nz) :: &
      rc_1_ice, rc_2_ice
    
    ! To test pdf parameters
    real( kind = core_rknd ), dimension(ngrdcol,nz) ::  &
      wm_clubb_pdf,    &
      rtm_clubb_pdf,   &
      thlm_clubb_pdf,  &
      wp2_clubb_pdf,   &
      rtp2_clubb_pdf,  &
      thlp2_clubb_pdf, &
      wp3_clubb_pdf,   &
      rtp3_clubb_pdf,  &
      thlp3_clubb_pdf, &
      Skw_clubb_pdf,   &
      Skrt_clubb_pdf,  &
      Skthl_clubb_pdf

    real( kind = core_rknd ) :: &
      chi_at_ice_sat1, chi_at_ice_sat2

    real( kind = core_rknd ) :: &
      beta,                         & ! CLUBB tunable parameter beta
      Skw_denom_coef,               & ! CLUBB tunable parameter Skw_denom_coef
      slope_coef_spread_DG_means_w, & ! CLUBB tunable parameter
      pdf_component_stdev_factor_w, & ! CLUBB tunable parameter
      coef_spread_DG_means_rt,      & ! CLUBB tunable parameter
      coef_spread_DG_means_thl        ! CLUBB tunable parameter

    logical, parameter :: &
      l_liq_ice_loading_test = .false. ! Temp. flag liq./ice water loading test

    integer :: k, i, j, hm_idx   ! Indices

#ifdef GFDL
    real ( kind = core_rknd ), parameter :: t1_combined = 273.16, &
                                            t2_combined = 268.16, &
                                            t3_combined = 238.16 
#endif

  !------------------------ Code Begins ----------------------------------

    ! Check whether the passive scalars are present.
    if ( sclr_dim > 0 ) then
      l_scalar_calc = .true.
    else
      l_scalar_calc = .false.
    end if


    ! Initialize to default values to prevent a runtime error
    if ( ( iiPDF_type /= iiPDF_ADG1 ) .and. ( iiPDF_type /= iiPDF_ADG2 ) ) then

      do k = 1, nz
        do i = 1, ngrdcol
          pdf_params%alpha_thl(i,k) = one_half
          pdf_params%alpha_rt(i,k) = one_half
        end do
      end do
      
      ! This allows for skewness to be clipped locally without passing the updated
      ! value back out.
      do k = 1, nz
        do i = 1, ngrdcol
          Skrt(i,k) = Skrt_in(i,k)
          Skthl(i,k) = Skthl_in(i,k)
          Sku(i,k) = Sku_in(i,k)
          Skv(i,k) = Skv_in(i,k)
        end do
      end do
        
      do j = 1, sclr_dim
        do k = 1, nz
          do i = 1, ngrdcol
            
            Sksclr(i,k,j) = Sksclr_in(i,k,j)
            
            if ( l_scalar_calc ) then
                alpha_sclr(i,k,j) = one_half
            end if
            
          end do
        end do
      end do

    end if


    ! Initialize to 0 to prevent a runtime error
    if ( iiPDF_type /= iiPDF_new .and. iiPDF_type /= iiPDF_new_hybrid ) then

      do k = 1, nz
        do i = 1, ngrdcol
          F_w(i,k) = zero
          F_rt(i,k) = zero
          F_thl(i,k) = zero
          min_F_w(i,k) = zero
          max_F_w(i,k) = zero
          min_F_rt(i,k) = zero
          max_F_rt(i,k) = zero
          min_F_thl(i,k) = zero
          max_F_thl(i,k) = zero
        end do
      end do

    end if

    ! Unpack CLUBB's tunable parameters
    if ( ( iiPDF_type == iiPDF_ADG1 ) .or. ( iiPDF_type == iiPDF_ADG2 ) ) then
       beta = clubb_params(ibeta)
    elseif ( iiPDF_type == iiPDF_new ) then
       slope_coef_spread_DG_means_w = clubb_params(islope_coef_spread_DG_means_w)
       pdf_component_stdev_factor_w = clubb_params(ipdf_component_stdev_factor_w)
       coef_spread_DG_means_rt = clubb_params(icoef_spread_DG_means_rt)
       coef_spread_DG_means_thl = clubb_params(icoef_spread_DG_means_thl)
    elseif ( iiPDF_type == iiPDF_new_hybrid ) then
       slope_coef_spread_DG_means_w = clubb_params(islope_coef_spread_DG_means_w)
       pdf_component_stdev_factor_w = clubb_params(ipdf_component_stdev_factor_w)
    end if
      

    ! To avoid recomputing
    do k = 1, nz
      do i = 1, ngrdcol
        sqrt_wp2(i,k) = sqrt( wp2(i,k) )
      end do
    end do

    ! Select the PDF closure method for the two-component PDF used by CLUBB for
    ! w, rt, theta-l, and passive scalar variables.
    ! Calculate the mixture fraction for the multivariate PDF, as well as both
    ! PDF component means and both PDF component variances for each of w, rt,
    ! theta-l, and passive scalar variables.
    if ( iiPDF_type == iiPDF_ADG1 ) then ! use ADG1
      
      call ADG1_pdf_driver( nz, ngrdcol,                                        & ! In
                            wm, rtm, thlm, um, vm,                              & ! In
                            wp2, rtp2, thlp2, up2, vp2,                         & ! In
                            Skw, wprtp, wpthlp, upwp, vpwp, sqrt_wp2,           & ! In
                            sigma_sqd_w, beta, mixt_frac_max_mag,               & ! In
                            sclrm, sclrp2, wpsclrp, l_scalar_calc,              & ! In
                            pdf_params%w_1, pdf_params%w_2,                     & ! Out
                            pdf_params%rt_1, pdf_params%rt_2,                   & ! Out
                            pdf_params%thl_1, pdf_params%thl_2,                 & ! Out
                            u_1, u_2, v_1, v_2,                                 & ! Out
                            pdf_params%varnce_w_1, pdf_params%varnce_w_2,       & ! Out
                            pdf_params%varnce_rt_1, pdf_params%varnce_rt_2,     & ! Out
                            pdf_params%varnce_thl_1, pdf_params%varnce_thl_2,   & ! Out
                            varnce_u_1, varnce_u_2,                             & ! Out
                            varnce_v_1, varnce_v_2,                             & ! Out
                            pdf_params%mixt_frac,                               & ! Out
                            pdf_params%alpha_rt, pdf_params%alpha_thl,          & ! Out
                            alpha_u, alpha_v,                                   & ! Out
                            sclr1, sclr2, varnce_sclr1,                         & ! Out
                            varnce_sclr2, alpha_sclr )                            ! Out
                            
    elseif ( iiPDF_type == iiPDF_ADG2 ) then ! use ADG2
      
      call ADG2_pdf_driver( gr, nz, ngrdcol,                                  & ! In
                            wm, rtm, thlm, wp2, rtp2, thlp2,                  & ! In
                            Skw, wprtp, wpthlp, sqrt_wp2, beta,               & ! In
                            sclrm, sclrp2, wpsclrp, l_scalar_calc,            & ! In
                            pdf_params%w_1, pdf_params%w_2,                   & ! Out
                            pdf_params%rt_1, pdf_params%rt_2,                 & ! Out
                            pdf_params%thl_1, pdf_params%thl_2,               & ! Out
                            pdf_params%varnce_w_1, pdf_params%varnce_w_2,     & ! Out
                            pdf_params%varnce_rt_1, pdf_params%varnce_rt_2,   & ! Out
                            pdf_params%varnce_thl_1, pdf_params%varnce_thl_2, & ! Out
                            pdf_params%mixt_frac,                             & ! Out
                            pdf_params%alpha_rt, pdf_params%alpha_thl,        & ! Out
                            sigma_sqd_w, sclr1, sclr2,                        & ! Out
                            varnce_sclr1, varnce_sclr2, alpha_sclr )            ! Out
                            
    elseif ( iiPDF_type == iiPDF_3D_Luhar ) then ! use 3D Luhar
      do i = 1, ngrdcol
        call Luhar_3D_pdf_driver( gr(i), &
                           wm(i,:), rtm(i,:), thlm(i,:), wp2(i,:), rtp2(i,:), thlp2(i,:),                             & ! In
                           Skw(i,:), Skrt(i,:), Skthl(i,:), wprtp(i,:), wpthlp(i,:),                             & ! In
                           pdf_params%w_1(i,:), pdf_params%w_2(i,:),                    & ! Out
                           pdf_params%rt_1(i,:), pdf_params%rt_2(i,:),                  & ! Out
                           pdf_params%thl_1(i,:), pdf_params%thl_2(i,:),                & ! Out
                           pdf_params%varnce_w_1(i,:), pdf_params%varnce_w_2(i,:),      & ! Out
                           pdf_params%varnce_rt_1(i,:), pdf_params%varnce_rt_2(i,:),    & ! Out
                           pdf_params%varnce_thl_1(i,:), pdf_params%varnce_thl_2(i,:),  & ! Out
                           pdf_params%mixt_frac(i,:) )                                    ! Out
      end do
    elseif ( iiPDF_type == iiPDF_new ) then ! use new PDF
      do i = 1, ngrdcol
        call new_pdf_driver( gr(i), wm(i,:), rtm(i,:), thlm(i,:), wp2(i,:), rtp2(i,:), thlp2(i,:), Skw(i,:),                   & ! In
                            wprtp(i,:), wpthlp(i,:), rtpthlp(i,:),                                     & ! In
                            slope_coef_spread_DG_means_w,                               & ! In
                            pdf_component_stdev_factor_w,                               & ! In
                            coef_spread_DG_means_rt,                                    & ! In
                            coef_spread_DG_means_thl,                                   & ! In
                            Skrt(i,:), Skthl(i,:),                                                & ! In/Out
                            pdf_params%w_1(i,:), pdf_params%w_2(i,:),                   & ! Out
                            pdf_params%rt_1(i,:), pdf_params%rt_2(i,:),                 & ! Out
                            pdf_params%thl_1(i,:), pdf_params%thl_2(i,:),               & ! Out
                            pdf_params%varnce_w_1(i,:), pdf_params%varnce_w_2(i,:),     & ! Out
                            pdf_params%varnce_rt_1(i,:), pdf_params%varnce_rt_2(i,:),   & ! Out
                            pdf_params%varnce_thl_1(i,:), pdf_params%varnce_thl_2(i,:), & ! Out
                            pdf_params%mixt_frac(i,:),                                  & ! Out
                            pdf_implicit_coefs_terms(i),                                   & ! Out
                            F_w(i,:), F_rt(i,:), F_thl(i,:), min_F_w(i,:), max_F_w(i,:),                         & ! Out
                            min_F_rt(i,:), max_F_rt(i,:), min_F_thl(i,:), max_F_thl(i,:) )                    ! Out
      end do
    elseif ( iiPDF_type == iiPDF_TSDADG ) then
      do i = 1, ngrdcol
        call tsdadg_pdf_driver( gr(i), &
                          wm(i,:), rtm(i,:), thlm(i,:), wp2(i,:), rtp2(i,:), thlp2(i,:),                                & ! In
                          Skw(i,:), Skrt(i,:), Skthl(i,:), wprtp(i,:), wpthlp(i,:),                                & ! In
                          pdf_params%w_1(i,:), pdf_params%w_2(i,:),                       & ! Out
                          pdf_params%rt_1(i,:), pdf_params%rt_2(i,:),                     & ! Out
                          pdf_params%thl_1(i,:), pdf_params%thl_2(i,:),                   & ! Out
                          pdf_params%varnce_w_1(i,:), pdf_params%varnce_w_2(i,:),         & ! Out
                          pdf_params%varnce_rt_1(i,:), pdf_params%varnce_rt_2(i,:),       & ! Out
                          pdf_params%varnce_thl_1(i,:), pdf_params%varnce_thl_2(i,:),     & ! Out
                          pdf_params%mixt_frac(i,:) )                                       ! Out
      end do
    elseif ( iiPDF_type == iiPDF_LY93 ) then ! use LY93
      do i = 1, ngrdcol
        call LY93_driver( gr(i), wm(i,:), rtm(i,:), thlm(i,:), wp2(i,:), rtp2(i,:),                          & ! In
                          thlp2(i,:), Skw(i,:), Skrt(i,:), Skthl(i,:),                           & ! In
                          pdf_params%w_1(i,:), pdf_params%w_2(i,:),                    & ! Out
                          pdf_params%rt_1(i,:), pdf_params%rt_2(i,:),                  & ! Out
                          pdf_params%thl_1(i,:), pdf_params%thl_2(i,:),                & ! Out
                          pdf_params%varnce_w_1(i,:), pdf_params%varnce_w_2(i,:),      & ! Out
                          pdf_params%varnce_rt_1(i,:), pdf_params%varnce_rt_2(i,:),    & ! Out
                          pdf_params%varnce_thl_1(i,:), pdf_params%varnce_thl_2(i,:),  & ! Out
                          pdf_params%mixt_frac(i,:) )                               ! Out
      end do
    elseif ( iiPDF_type == iiPDF_new_hybrid ) then ! use new hybrid PDF
      do i = 1, ngrdcol
        call new_hybrid_pdf_driver( gr(i), wm(i,:), rtm(i,:), thlm(i,:), um(i,:), vm(i,:),          & ! In
                                    wp2(i,:), rtp2(i,:), thlp2(i,:), up2(i,:), vp2(i,:),         & ! In
                                    Skw(i,:), wprtp(i,:), wpthlp(i,:), upwp(i,:), vpwp(i,:),     & ! In
                                    sclrm(i,:,:), sclrp2(i,:,:), wpsclrp(i,:,:),             & ! In
                                    gamma_Skw_fnc(i,:),                      & ! In
                                    slope_coef_spread_DG_means_w,       & ! In
                                    pdf_component_stdev_factor_w,       & ! In
                                    Skrt(i,:), Skthl(i,:), Sku(i,:), Skv(i,:), Sksclr(i,:,:),      & ! I/O
                                    pdf_params%w_1(i,:), pdf_params%w_2(i,:),     & ! Out
                                    pdf_params%rt_1(i,:), pdf_params%rt_2(i,:),   & ! Out
                                    pdf_params%thl_1(i,:), pdf_params%thl_2(i,:), & ! Out
                                    u_1(i,:), u_2(i,:), v_1(i,:), v_2(i,:),                 & ! Out
                                    pdf_params%varnce_w_1(i,:),         & ! Out
                                    pdf_params%varnce_w_2(i,:),         & ! Out
                                    pdf_params%varnce_rt_1(i,:),        & ! Out
                                    pdf_params%varnce_rt_2(i,:),        & ! Out
                                    pdf_params%varnce_thl_1(i,:),       & ! Out
                                    pdf_params%varnce_thl_2(i,:),       & ! Out
                                    varnce_u_1(i,:), varnce_u_2(i,:),             & ! Out
                                    varnce_v_1(i,:), varnce_v_2(i,:),             & ! Out
                                    sclr1(i,:,:), sclr2(i,:,:),                       & ! Out
                                    varnce_sclr1(i,:,:), varnce_sclr2(i,:,:),         & ! Out
                                    pdf_params%mixt_frac(i,:),          & ! Out
                                    pdf_implicit_coefs_terms(i),           & ! Out
                                    F_w(i,:), min_F_w(i,:), max_F_w(i,:)               ) ! Out
      end do
      
      ! The calculation of skewness of rt, thl, u, v, and scalars is hard-wired
      ! for use with the ADG1 code, which contains the variable sigma_sqd_w.
      ! In order to use an equivalent expression for these skewnesses using the
      ! new hybrid PDF (without doing more recoding), set the value of
      ! sigma_sqd_w to 1 - F_w.
      do k = 1, nz
        do i = 1, ngrdcol
          sigma_sqd_w(i,k) = one - F_w(i,k)
        end do
      end do

    end if ! iiPDF_type
    
    ! Calculate the PDF component correlations of rt and thl.
    call calc_comp_corrs_binormal( &
                         rtpthlp, rtm, thlm,                                    & ! In
                         pdf_params%rt_1, pdf_params%rt_2,                      & ! In
                         pdf_params%thl_1, pdf_params%thl_2,                    & ! In
                         pdf_params%varnce_rt_1, pdf_params%varnce_rt_2,        & ! In
                         pdf_params%varnce_thl_1, pdf_params%varnce_thl_2,      & ! In
                         pdf_params%mixt_frac,                                  & ! In
                         pdf_params%corr_rt_thl_1, pdf_params%corr_rt_thl_2 )     ! Out

    if ( iiPDF_type == iiPDF_ADG1 .or. iiPDF_type == iiPDF_ADG2 &
         .or. iiPDF_type == iiPDF_new_hybrid ) then

      ! These PDF types define corr_w_rt_1, corr_w_rt_2, corr_w_thl_1, and
      ! corr_w_thl_2 to all have a value of 0, so skip the calculation.
      ! The values of corr_u_w_1, corr_u_w_2, corr_v_w_1, and corr_v_w_2 are
      ! all defined to be 0, as well.
      do k = 1, nz
        do i = 1, ngrdcol
          pdf_params%corr_w_rt_1(i,k)  = zero
          pdf_params%corr_w_rt_2(i,k)  = zero
          pdf_params%corr_w_thl_1(i,k) = zero
          pdf_params%corr_w_thl_2(i,k) = zero
          corr_u_w_1(i,k)   = zero
          corr_u_w_2(i,k)   = zero
          corr_v_w_1(i,k)   = zero
          corr_v_w_2(i,k)   = zero
        end do
      end do

    else

      ! Calculate the PDF component correlations of w and rt.
      call calc_comp_corrs_binormal( & 
                              wprtp, wm, rtm,                                   & ! In
                              pdf_params%w_1, pdf_params%w_2,                   & ! In
                              pdf_params%rt_1, pdf_params%rt_2,                 & ! In
                              pdf_params%varnce_w_1, pdf_params%varnce_w_2,     & ! In
                              pdf_params%varnce_rt_1, pdf_params%varnce_rt_2,   & ! In
                              pdf_params%mixt_frac,                             & ! In
                              pdf_params%corr_w_rt_1, pdf_params%corr_w_rt_2 )    ! Out
    
      ! Calculate the PDF component correlations of w and thl.
      call calc_comp_corrs_binormal( &
                        wpthlp, wm, thlm,                                     & ! In
                        pdf_params%w_1, pdf_params%w_2,                       & ! In
                        pdf_params%thl_1, pdf_params%thl_2,                   & ! In
                        pdf_params%varnce_w_1, pdf_params%varnce_w_2,         & ! In
                        pdf_params%varnce_thl_1, pdf_params%varnce_thl_2,     & ! In
                        pdf_params%mixt_frac,                                 & ! In
                        pdf_params%corr_w_thl_1, pdf_params%corr_w_thl_2 )      ! Out
    end if
      
    if ( l_scalar_calc ) then

      ! Calculate the PDF component correlations of a passive scalar and thl.
      do j = 1, sclr_dim
        call calc_comp_corrs_binormal( &
                            sclrpthlp(:,:,j), sclrm(:,:,j), thlm,               & ! In
                             sclr1(:,:,j), sclr2(:,:,j),                        & ! In
                             pdf_params%thl_1, pdf_params%thl_2,                & ! In
                             varnce_sclr1(:,:,j), varnce_sclr2(:,:,j),          & ! In
                             pdf_params%varnce_thl_1, pdf_params%varnce_thl_2,  & ! In
                             pdf_params%mixt_frac,                              & ! In
                             corr_sclr_thl_1(:,:,j), corr_sclr_thl_2(:,:,j) )     ! Out
      end do

      ! Calculate the PDF component correlations of a passive scalar and rt.
      do j = 1, sclr_dim
        call calc_comp_corrs_binormal( &
                             sclrprtp(:,:,j), sclrm(:,:,j), rtm,                & ! In
                             sclr1(:,:,j), sclr2(:,:,j),                        & ! In
                             pdf_params%rt_1, pdf_params%rt_2,                  & ! In
                             varnce_sclr1(:,:,j), varnce_sclr2(:,:,j),          & ! In
                             pdf_params%varnce_rt_1, pdf_params%varnce_rt_2,    & ! In
                             pdf_params%mixt_frac,                              & ! In
                             corr_sclr_rt_1(:,:,j), corr_sclr_rt_2(:,:,j) )       ! Out
      end do

      if ( iiPDF_type == iiPDF_ADG1 .or. iiPDF_type == iiPDF_ADG2 &
           .or. iiPDF_type == iiPDF_new_hybrid ) then

        ! These PDF types define all PDF component correlations involving w
        ! to have a value of 0, so skip the calculation.
        do j = 1, sclr_dim
          do k = 1, nz
            do i = 1, ngrdcol
              corr_w_sclr_1(i,k,j) = zero
              corr_w_sclr_2(i,k,j) = zero
            end do
          end do
        end do

      else

        ! Calculate the PDF component correlations of w and a passive scalar.
        do j = 1, sclr_dim
          call calc_comp_corrs_binormal( &
                                 wpsclrp(:,:,j), wm, sclrm(:,:,j),              & ! In
                                 pdf_params%w_1, pdf_params%w_2,                & ! In
                                 sclr1(:,:,j), sclr2(:,:,j),                    & ! In
                                 pdf_params%varnce_w_1, pdf_params%varnce_w_2,  & ! In
                                 varnce_sclr1(:,:,j), varnce_sclr2(:,:,j),      & ! In
                                 pdf_params%mixt_frac,                          & ! In
                                 corr_w_sclr_1(:,:,j), corr_w_sclr_2(:,:,j) )     ! Out

        end do
        
      end if

    end if


    ! Compute higher order moments (these are interactive)
    call calc_wp2xp_pdf( nz, ngrdcol,                                       &
                         wm, rtm, pdf_params%w_1, pdf_params%w_2,           &
                         pdf_params%rt_1, pdf_params%rt_2,                  &
                         pdf_params%varnce_w_1, pdf_params%varnce_w_2,      &
                         pdf_params%varnce_rt_1, pdf_params%varnce_rt_2,    &
                         pdf_params%corr_w_rt_1, pdf_params%corr_w_rt_2,    &
                         pdf_params%mixt_frac,                              &
                         wp2rtp )

    call calc_wp2xp_pdf( nz, ngrdcol,                                          &
                         wm, thlm, pdf_params%w_1, pdf_params%w_2,             &
                         pdf_params%thl_1, pdf_params%thl_2,                   &
                         pdf_params%varnce_w_1, pdf_params%varnce_w_2,         &
                         pdf_params%varnce_thl_1, pdf_params%varnce_thl_2,     &
                         pdf_params%corr_w_thl_1, pdf_params%corr_w_thl_2,     &
                         pdf_params%mixt_frac,                                 &
                         wp2thlp )
    
    ! Compute higher order moments (these may be interactive)
    call calc_wpxp2_pdf( nz, ngrdcol, &
                         wm, um, pdf_params%w_1, pdf_params%w_2, &
                         u_1, u_2, &
                         pdf_params%varnce_w_1, pdf_params%varnce_w_2,   &
                         varnce_u_1, varnce_u_2, &
                         corr_u_w_1, corr_u_w_2, &
                         pdf_params%mixt_frac, &
                         wpup2 )
                         
    call calc_wpxp2_pdf( nz, ngrdcol, &
                         wm, vm, pdf_params%w_1, pdf_params%w_2, &
                         v_1, v_2, &
                         pdf_params%varnce_w_1, pdf_params%varnce_w_2,   &
                         varnce_v_1, varnce_v_2, &
                         corr_v_w_1, corr_v_w_2, &
                         pdf_params%mixt_frac, &
                         wpvp2 )
    
    call calc_wp2xp2_pdf( nz, ngrdcol, &
                          wm, um, pdf_params%w_1, pdf_params%w_2, &
                          u_1, u_2, &
                          pdf_params%varnce_w_1, pdf_params%varnce_w_2, &
                          varnce_u_1, varnce_u_2, &
                          corr_u_w_1, corr_u_w_2, &
                          pdf_params%mixt_frac, &
                          wp2up2 )

    call calc_wp2xp2_pdf( nz, ngrdcol, &
                          wm, vm, pdf_params%w_1, pdf_params%w_2, &
                          v_1, v_2, &
                          pdf_params%varnce_w_1, pdf_params%varnce_w_2, &
                          varnce_v_1, varnce_v_2, &
                          corr_v_w_1, corr_v_w_2, &
                          pdf_params%mixt_frac, &
                          wp2vp2 )
                          
    call calc_wp4_pdf( nz, ngrdcol, &
                       wm, pdf_params%w_1, pdf_params%w_2, &
                       pdf_params%varnce_w_1, pdf_params%varnce_w_2,    &
                       pdf_params%mixt_frac, &
                       wp4 )

    if ( l_explicit_turbulent_adv_xpyp .or. iwprtp2 > 0 ) then
      call calc_wpxp2_pdf( nz, ngrdcol, &
                           wm, rtm, pdf_params%w_1, pdf_params%w_2,        &
                           pdf_params%rt_1, pdf_params%rt_2,               &
                           pdf_params%varnce_w_1, pdf_params%varnce_w_2,   &
                           pdf_params%varnce_rt_1, pdf_params%varnce_rt_2, &
                           pdf_params%corr_w_rt_1, pdf_params%corr_w_rt_2, &
                           pdf_params%mixt_frac, &
                           wprtp2 )
    end if

    if ( l_explicit_turbulent_adv_xpyp .or. iwpthlp2 > 0 ) then
      call calc_wpxp2_pdf( nz, ngrdcol, &
                           wm, thlm, pdf_params%w_1, pdf_params%w_2,          &
                           pdf_params%thl_1, pdf_params%thl_2,                &
                           pdf_params%varnce_w_1, pdf_params%varnce_w_2,      &
                           pdf_params%varnce_thl_1, pdf_params%varnce_thl_2,  &
                           pdf_params%corr_w_thl_1, pdf_params%corr_w_thl_2,  &
                           pdf_params%mixt_frac, &
                           wpthlp2 )
    end if

    if ( l_explicit_turbulent_adv_xpyp .or. iwprtpthlp > 0 ) then
      
      call calc_wpxpyp_pdf( nz, ngrdcol, &
                            wm, rtm, thlm, pdf_params%w_1, pdf_params%w_2,      &
                            pdf_params%rt_1, pdf_params%rt_2,                   &
                            pdf_params%thl_1, pdf_params%thl_2,                 &
                            pdf_params%varnce_w_1, pdf_params%varnce_w_2,       &
                            pdf_params%varnce_rt_1, pdf_params%varnce_rt_2,     &
                            pdf_params%varnce_thl_1, pdf_params%varnce_thl_2,   &
                            pdf_params%corr_w_rt_1, pdf_params%corr_w_rt_2,     &
                            pdf_params%corr_w_thl_1, pdf_params%corr_w_thl_2,   &
                            pdf_params%corr_rt_thl_1, pdf_params%corr_rt_thl_2, &
                            pdf_params%mixt_frac, &
                            wprtpthlp )
    end if


    ! Scalar Addition to higher order moments
    if ( l_scalar_calc ) then

      do j = 1, sclr_dim
        call calc_wp2xp_pdf( nz, ngrdcol,                                       &
                             wm, sclrm(:,:,j), pdf_params%w_1, pdf_params%w_2,  &
                             sclr1(:,:,j), sclr2(:,:,j),                        &
                             pdf_params%varnce_w_1, pdf_params%varnce_w_2,      &
                             varnce_sclr1(:,:,j), varnce_sclr2(:,:,j),          &
                             corr_w_sclr_1(:,:,j), corr_w_sclr_2(:,:,j),        &
                             pdf_params%mixt_frac,                              &
                             wp2sclrp(:,:,j) )
      end do
      
      do j = 1, sclr_dim 
        call calc_wpxp2_pdf( nz, ngrdcol, &
                             wm, sclrm(:,:,j), pdf_params%w_1, pdf_params%w_2, &
                             sclr1(:,:,j), sclr2(:,:,j),                         &
                             pdf_params%varnce_w_1, pdf_params%varnce_w_2,   &
                             varnce_sclr1(:,:,j), varnce_sclr2(:,:,j),           &
                             corr_w_sclr_1(:,:,j), corr_w_sclr_2(:,:,j),         &
                             pdf_params%mixt_frac, &
                             wpsclrp2(:,:,j) )
      end do
       
      do j = 1, sclr_dim
        call calc_wpxpyp_pdf( nz, ngrdcol, &
                              wm, sclrm(:,:,j), rtm, pdf_params%w_1, pdf_params%w_2,  &
                              sclr1(:,:,j), sclr2(:,:,j),                             &
                              pdf_params%rt_1, pdf_params%rt_2,                       &
                              pdf_params%varnce_w_1, pdf_params%varnce_w_2,           &
                              varnce_sclr1(:,:,j), varnce_sclr2(:,:,j),               &
                              pdf_params%varnce_rt_1, pdf_params%varnce_rt_2,         &
                              corr_w_sclr_1(:,:,j), corr_w_sclr_2(:,:,j),             &
                              pdf_params%corr_w_rt_1, pdf_params%corr_w_rt_2,         &
                              corr_sclr_rt_1(:,:,j), corr_sclr_rt_2(:,:,j),           &
                              pdf_params%mixt_frac, &
                              wpsclrprtp(:,:,j) )
      end do
        
      do j = 1, sclr_dim
        call calc_wpxpyp_pdf( nz, ngrdcol, &
                              wm, sclrm(:,:,j), thlm, pdf_params%w_1, pdf_params%w_2,   &
                              sclr1(:,:,j), sclr2(:,:,j),                               &
                              pdf_params%thl_1, pdf_params%thl_2,                       &
                              pdf_params%varnce_w_1, pdf_params%varnce_w_2,             &
                              varnce_sclr1(:,:,j), varnce_sclr2(:,:,j),                 &
                              pdf_params%varnce_thl_1, pdf_params%varnce_thl_2,         &
                              corr_w_sclr_1(:,:,j), corr_w_sclr_2(:,:,j),               &
                              pdf_params%corr_w_thl_1, pdf_params%corr_w_thl_2,         &
                              corr_sclr_thl_1(:,:,j), corr_sclr_thl_2(:,:,j),           &
                              pdf_params%mixt_frac, &
                              wpsclrpthlp(:,:,j) )
      end do
      
    end if

    ! Compute higher order moments that include theta_v.

    ! First compute some preliminary quantities.
    ! "1" denotes first Gaussian; "2" denotes 2nd Gaussian
    ! liq water temp (Sommeria & Deardorff 1977 (SD), eqn. 3)
    
    do k = 1, nz
      do i = 1, ngrdcol
        tl1(i,k)  = pdf_params%thl_1(i,k)*exner(i,k)
        tl2(i,k)  = pdf_params%thl_2(i,k)*exner(i,k)
      end do
    end do

#ifdef GFDL
    if ( sclr_dim > 0  .and.  (.not. do_liquid_only_in_clubb) ) then ! h1g, 2010-06-16 begin mod

      do i = 1, ngrdcol
        where ( tl1(i,:) > t1_combined )
          pdf_params%rsatl_1(i,:) = sat_mixrat_liq( p_in_Pa(i,:), tl1(i,:) )
        elsewhere ( tl1(i,:) > t2_combined )
          pdf_params%rsatl_1(i,:) = sat_mixrat_liq( p_in_Pa(i,:), tl1(i,:) ) &
                    * (tl1(i,:) - t2_combined)/(t1_combined - t2_combined) &
                    + sat_mixrat_ice( p_in_Pa(i,:), tl1(i,:) ) &
                      * (t1_combined - tl1(i,:))/(t1_combined - t2_combined)
        elsewhere ( tl1(i,:) > t3_combined )
          pdf_params%rsatl_1(i,:) = sat_mixrat_ice( p_in_Pa(i,:), tl1(i,:) ) &
                    + sat_mixrat_ice( p_in_Pa(i,:), tl1(i,:) ) * (RH_crit(i, :, 1, 1) -one ) &
                      * ( t2_combined -tl1(i,:))/(t2_combined - t3_combined)
        elsewhere
          pdf_params%rsatl_1(i,:) = sat_mixrat_ice( p_in_Pa(i,:), tl1(i,:) ) * RH_crit(i, :, 1, 1)
        endwhere

        where ( tl2(i,:) > t1_combined )
          pdf_params%rsatl_2(i,:) = sat_mixrat_liq( p_in_Pa(i,:), tl2(i,:) )
        elsewhere ( tl2(i,:) > t2_combined )
          pdf_params%rsatl_2(i,:) = sat_mixrat_liq( p_in_Pa(i,:), tl2(i,:) ) &
                    * (tl2(i,:) - t2_combined)/(t1_combined - t2_combined) &
                    + sat_mixrat_ice( p_in_Pa(i,:), tl2(i,:) ) &
                      * (t1_combined - tl2(i,:))/(t1_combined - t2_combined)
        elsewhere ( tl2(i,:) > t3_combined )
          pdf_params%rsatl_2(i,:) = sat_mixrat_ice( p_in_Pa(i,:), tl2(i,:) ) &
                    + sat_mixrat_ice( p_in_Pa(i,:), tl2(i,:) )* (RH_crit(i, :, 1, 2) -one) &
                      * ( t2_combined -tl2(i,:))/(t2_combined - t3_combined)
        elsewhere
          pdf_params%rsatl_2(i,:) = sat_mixrat_ice( p_in_Pa(i,:), tl2(i,:) ) * RH_crit(i, :, 1, 2)
        endwhere
        
      end do

    else ! sclr_dim <= 0  or  do_liquid_only_in_clubb = .T.

      pdf_params%rsatl_1(:,:) = sat_mixrat_liq( p_in_Pa(:,:), tl1(:,:) )
      pdf_params%rsatl_2(:,:) = sat_mixrat_liq( p_in_Pa(:,:), tl2(:,:) )

    end if !sclr_dim > 0
      
    ! Determine whether to compute ice_supersat_frac. We do not compute
    ! ice_supersat_frac for GFDL (unless do_liquid_only_in_clubb is true),
    ! because liquid and ice are both fed into rtm, ruining the calculation.
    if (do_liquid_only_in_clubb) then
      l_calc_ice_supersat_frac = .true.
    else
      l_calc_ice_supersat_frac = .false.
    end if

#else
    pdf_params%rsatl_1(:,:) = sat_mixrat_liq( p_in_Pa(:,:), tl1(:,:) )
    pdf_params%rsatl_2(:,:) = sat_mixrat_liq( p_in_Pa(:,:), tl2(:,:) ) ! h1g, 2010-06-16 end mod

    l_calc_ice_supersat_frac = .true.
#endif

    call transform_pdf_chi_eta_component( nz, ngrdcol, &
                                          tl1, pdf_params%rsatl_1, pdf_params%rt_1, exner,  & ! In
                                          pdf_params%varnce_thl_1, pdf_params%varnce_rt_1,  & ! In
                                          pdf_params%corr_rt_thl_1, pdf_params%chi_1,       & ! In
                                          pdf_params%crt_1, pdf_params%cthl_1,              & ! Out
                                          pdf_params%stdev_chi_1, pdf_params%stdev_eta_1,   & ! Out
                                          pdf_params%covar_chi_eta_1,                       & ! Out
                                          pdf_params%corr_chi_eta_1 )                         ! Out
    
      
    ! Calculate cloud fraction component for pdf 1
    call calc_liquid_cloud_frac_component( nz, ngrdcol, &
                                           pdf_params%chi_1, pdf_params%stdev_chi_1,    & ! In
                                           pdf_params%cloud_frac_1, pdf_params%rc_1 )     ! Out

    ! Calc ice_supersat_frac
    if ( l_calc_ice_supersat_frac ) then

      do k = 1, nz
        do i = 1, ngrdcol

          if ( tl1(i,k) <= T_freeze_K ) then

            ! Temperature is freezing, we must compute chi_at_ice_sat and
            ! calculate the new cloud_frac_component
            chi_at_ice_sat1 = ( sat_mixrat_ice( p_in_Pa(i,k), tl1(i,k) ) - pdf_params%rsatl_1(i,k) ) &
                              * pdf_params%crt_1(i,k)

            call calc_cloud_frac_component( pdf_params%chi_1(i,k), pdf_params%stdev_chi_1(i,k), &!in
                                            chi_at_ice_sat1, & ! intent(in)
                                            pdf_params%ice_supersat_frac_1(i,k), rc_1_ice(i,k) )! out
          else

              ! Temperature is warmer than freezing, the ice_supersat_frac calculation is
              ! the same as cloud_frac
              pdf_params%ice_supersat_frac_1(i,k) = pdf_params%cloud_frac_1(i,k)
              rc_1_ice(i,k) = pdf_params%rc_1(i,k)

          end if
          
        end do
      end do

    end if

    call transform_pdf_chi_eta_component( nz, ngrdcol, &
                                          tl2, pdf_params%rsatl_2, pdf_params%rt_2, exner,  & ! In
                                          pdf_params%varnce_thl_2, pdf_params%varnce_rt_2,  & ! In
                                          pdf_params%corr_rt_thl_2, pdf_params%chi_2,       & ! In
                                          pdf_params%crt_2, pdf_params%cthl_2,              & ! Out
                                          pdf_params%stdev_chi_2, pdf_params%stdev_eta_2,   & ! Out
                                          pdf_params%covar_chi_eta_2,                       & ! Out
                                          pdf_params%corr_chi_eta_2 )                         ! Out

      
    ! Calculate cloud fraction component for pdf 2
    call calc_liquid_cloud_frac_component( nz, ngrdcol, &
                                           pdf_params%chi_2, pdf_params%stdev_chi_2,    & ! In
                                           pdf_params%cloud_frac_2, pdf_params%rc_2 )     ! Out

    ! Calc ice_supersat_frac
    if ( l_calc_ice_supersat_frac ) then

        do k = 1, nz
          do i = 1, ngrdcol

            if ( tl2(i,k) <= T_freeze_K ) then

              ! Temperature is freezing, we must compute chi_at_ice_sat and 
              ! calculate the new cloud_frac_component
              chi_at_ice_sat2 = ( sat_mixrat_ice( p_in_Pa(i,k), tl2(i,k) ) - pdf_params%rsatl_2(i,k) ) &
                                * pdf_params%crt_2(i,k)

              call calc_cloud_frac_component( pdf_params%chi_2(i,k), & ! intent(in)
                                              pdf_params%stdev_chi_2(i,k), & ! intent(in)
                                              chi_at_ice_sat2, & ! intent(in)
                                              pdf_params%ice_supersat_frac_2(i,k), & ! intent(out)
                                              rc_2_ice(i,k) ) ! intent(out)
            else

                ! Temperature is warmer than freezing, the ice_supersat_frac calculation is 
                ! the same as cloud_frac
                pdf_params%ice_supersat_frac_2(i,k) = pdf_params%cloud_frac_2(i,k)
                rc_2_ice(i,k) = pdf_params%rc_2(i,k)

            end if
            
          end do
        end do

        ! Compute ice cloud fraction, ice_supersat_frac
        do k = 1, nz
          do i = 1, ngrdcol
            ice_supersat_frac(i,k) = pdf_params%mixt_frac(i,k) &
                                     * pdf_params%ice_supersat_frac_1(i,k) &
                                     + ( one - pdf_params%mixt_frac(i,k) ) &
                                       * pdf_params%ice_supersat_frac_2(i,k)
          end do
        end do

    else 

      ! ice_supersat_frac will be garbage if computed as above
      do k = 1, nz
        do i = 1, ngrdcol
          ice_supersat_frac(i,k) = 0.0_core_rknd
        end do
      end do

      if (clubb_at_least_debug_level( 1 )) then
          write(fstderr,*) "Warning: ice_supersat_frac has garbage values if &
                          & do_liquid_only_in_clubb = .false."
      end if

    end if ! l_calc_ice_supersat_frac


    ! Compute cloud fraction and mean cloud water mixing ratio.
    ! Reference:
    ! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:anl_int_cloud_terms
    do k = 1, nz
      do i = 1, ngrdcol
        cloud_frac(i,k) = pdf_params%mixt_frac(i,k) * pdf_params%cloud_frac_1(i,k) &
                     + ( one - pdf_params%mixt_frac(i,k) ) * pdf_params%cloud_frac_2(i,k)
        rcm(i,k) = pdf_params%mixt_frac(i,k) * pdf_params%rc_1(i,k) + ( one - pdf_params%mixt_frac(i,k) ) &
                                                                 * pdf_params%rc_2(i,k)
        rcm(i,k) = max( zero_threshold, rcm(i,k) )
      end do
    end do

    if ( iiPDF_type == iiPDF_ADG1 .or. iiPDF_type == iiPDF_ADG2 &
         .or. iiPDF_type == iiPDF_new_hybrid ) then

      ! corr_w_rt and corr_w_thl are zero for these pdf types so
      ! corr_w_chi and corr_w_eta are zero as well
      do k = 1, nz
        do i = 1, ngrdcol
          pdf_params%corr_w_chi_1(i,k) = zero
          pdf_params%corr_w_chi_2(i,k) = zero
          pdf_params%corr_w_eta_1(i,k) = zero
          pdf_params%corr_w_eta_2(i,k) = zero
        end do
      end do

    else 
        
      ! Correlation of w and chi for each component.
      pdf_params%corr_w_chi_1 &
      = calc_corr_chi_x( pdf_params%crt_1, pdf_params%cthl_1, &
                         sqrt(pdf_params%varnce_rt_1), sqrt(pdf_params%varnce_thl_1), &
                         pdf_params%stdev_chi_1, &
                         pdf_params%corr_w_rt_1, pdf_params%corr_w_thl_1 )

      pdf_params%corr_w_chi_2 &
      = calc_corr_chi_x( pdf_params%crt_2, pdf_params%cthl_2, &
                         sqrt(pdf_params%varnce_rt_2), sqrt(pdf_params%varnce_thl_2), &
                         pdf_params%stdev_chi_2, pdf_params%corr_w_rt_2, &
                         pdf_params%corr_w_thl_2 )

      ! Correlation of w and eta for each component.
      pdf_params%corr_w_eta_1 &
      = calc_corr_eta_x( pdf_params%crt_1, pdf_params%cthl_1, &
                         sqrt(pdf_params%varnce_rt_1), sqrt(pdf_params%varnce_thl_1), &
                         pdf_params%stdev_eta_1, pdf_params%corr_w_rt_1, &
                         pdf_params%corr_w_thl_1 )

      pdf_params%corr_w_eta_2 &
      = calc_corr_eta_x( pdf_params%crt_2, pdf_params%cthl_2, &
                         sqrt(pdf_params%varnce_rt_2), sqrt(pdf_params%varnce_thl_2), &
                         pdf_params%stdev_eta_2, pdf_params%corr_w_rt_2, &
                         pdf_params%corr_w_thl_2 )

    end if

    
    ! Compute moments that depend on theta_v
    ! 
    ! The moments that depend on th_v' are calculated based on an approximated
    ! and linearized form of the theta_v equation:
    ! 
    ! theta_v = theta_l + { (R_v/R_d) - 1 } * thv_ds * r_t
    !                   + [ {L_v/(C_p*exner)} - (R_v/R_d) * thv_ds ] * r_c;
    ! 
    ! and therefore:
    ! 
    ! th_v' = th_l' + { (R_v/R_d) - 1 } * thv_ds * r_t'
    !               + [ {L_v/(C_p*exner)} - (R_v/R_d) * thv_ds ] * r_c';
    ! 
    ! where thv_ds is used as a reference value to approximate theta_l.
    ! 
    ! Reference:
    ! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:anl_int_buoy_terms
    
    ! Calculate the contributions to <w'rc'>, <w'^2 rc'>, <rt'rc'>, and
    ! <thl'rc'> from the 1st PDF component.
    call calc_xprcp_component( nz, ngrdcol,                                          & ! In
                               wm, rtm, thlm, um, vm, rcm,                           & ! In
                               pdf_params%w_1, pdf_params%rt_1,                      & ! In
                               pdf_params%thl_1, u_1, v_1,                           & ! In
                               pdf_params%varnce_w_1, pdf_params%chi_1,              & ! In
                               pdf_params%stdev_chi_1, pdf_params%stdev_eta_1,       & ! In
                               pdf_params%corr_w_chi_1, pdf_params%corr_chi_eta_1,   & ! In
                               pdf_params%crt_1, pdf_params%cthl_1,                  & ! In
                               pdf_params%rc_1, pdf_params%cloud_frac_1, iiPDF_type, & ! In
                               wprcp_contrib_comp_1, wp2rcp_contrib_comp_1,          & ! Out
                               rtprcp_contrib_comp_1, thlprcp_contrib_comp_1,        & ! Out
                               uprcp_contrib_comp_1, vprcp_contrib_comp_1 )            ! Out

    call calc_xprcp_component( nz, ngrdcol,                                          & ! In 
                               wm, rtm, thlm, um, vm, rcm,                           & ! In
                               pdf_params%w_2, pdf_params%rt_2,                      & ! In
                               pdf_params%thl_2, u_2, v_2,                           & ! In
                               pdf_params%varnce_w_2, pdf_params%chi_2,              & ! In
                               pdf_params%stdev_chi_2, pdf_params%stdev_eta_2,       & ! In
                               pdf_params%corr_w_chi_2, pdf_params%corr_chi_eta_2,   & ! In
                               pdf_params%crt_2, pdf_params%cthl_2,                  & ! In
                               pdf_params%rc_2, pdf_params%cloud_frac_2, iiPDF_type, & ! In
                               wprcp_contrib_comp_2, wp2rcp_contrib_comp_2,          & ! Out
                               rtprcp_contrib_comp_2, thlprcp_contrib_comp_2,        & ! Out
                               uprcp_contrib_comp_2, vprcp_contrib_comp_2 )            ! Out

    
    ! Calculate rc_coef, which is the coefficient on <x'rc'> in the <x'thv'> equation.
    do k = 1, nz
      do i = 1, ngrdcol
        
        rc_coef(i,k) = Lv / ( exner(i,k) * Cp ) - ep2 * thv_ds(i,k)

        ! Calculate <w'rc'>, <w'^2 rc'>, <rt'rc'>, and <thl'rc'>.
        wprcp(i,k) = pdf_params%mixt_frac(i,k) * wprcp_contrib_comp_1(i,k) &
                     + ( one - pdf_params%mixt_frac(i,k) ) * wprcp_contrib_comp_2(i,k)

        wp2rcp(i,k) = pdf_params%mixt_frac(i,k) * wp2rcp_contrib_comp_1(i,k) &
                      + ( one - pdf_params%mixt_frac(i,k) ) * wp2rcp_contrib_comp_2(i,k)

        rtprcp(i,k) = pdf_params%mixt_frac(i,k) * rtprcp_contrib_comp_1(i,k) & 
                      + ( one - pdf_params%mixt_frac(i,k) ) * rtprcp_contrib_comp_2(i,k)

        thlprcp(i,k) = pdf_params%mixt_frac(i,k) * thlprcp_contrib_comp_1(i,k) &
                       + ( one - pdf_params%mixt_frac(i,k) ) * thlprcp_contrib_comp_2(i,k)

        uprcp(i,k) = pdf_params%mixt_frac(i,k) * uprcp_contrib_comp_1(i,k) &
                     + ( one - pdf_params%mixt_frac(i,k) ) * uprcp_contrib_comp_2(i,k)

        vprcp(i,k) = pdf_params%mixt_frac(i,k) * vprcp_contrib_comp_1(i,k) &
                     + ( one - pdf_params%mixt_frac(i,k) ) * vprcp_contrib_comp_2(i,k)
      end do
    end do

    ! Calculate <w'thv'>, <w'^2 thv'>, <rt'thv'>, and <thl'thv'>.
    do k = 1, nz
      do i = 1, ngrdcol
        wpthvp(i,k)  = wpthlp(i,k)  + ep1 * thv_ds(i,k) * wprtp(i,k)   + rc_coef(i,k) * wprcp(i,k)
        wp2thvp(i,k) = wp2thlp(i,k) + ep1 * thv_ds(i,k) * wp2rtp(i,k)  + rc_coef(i,k) * wp2rcp(i,k)
        rtpthvp(i,k) = rtpthlp(i,k) + ep1 * thv_ds(i,k) * rtp2(i,k)    + rc_coef(i,k) * rtprcp(i,k)
        thlpthvp(i,k)= thlp2(i,k)   + ep1 * thv_ds(i,k) * rtpthlp(i,k) + rc_coef(i,k) * thlprcp(i,k)
      end do
    end do


    ! Add the precipitation loading term in the <x'thv'> equation.
    if ( l_liq_ice_loading_test ) then

       do hm_idx = 1, hydromet_dim, 1

          if ( l_mix_rat_hm(hm_idx) ) then

            do k = 1, nz
              do i = 1, ngrdcol
                wp2thvp(i,k)  = wp2thvp(i,k)  - thv_ds(i,k) * wp2hmp(i,k,hm_idx)
                wpthvp(i,k)   = wpthvp(i,k)   - thv_ds(i,k) * wphydrometp(i,k,hm_idx)
                thlpthvp(i,k) = thlpthvp(i,k) - thv_ds(i,k) * thlphmp(i,k,hm_idx)
                rtpthvp(i,k)  = rtpthvp(i,k)  - thv_ds(i,k) * rtphmp(i,k,hm_idx)
              end do
            end do

          end if

       end do

    end if
    
    ! Account for subplume correlation of scalar, theta_v.
    ! See Eqs. A13, A8 from Larson et al. (2002) ``Small-scale...''
    !  where the ``scalar'' in this paper is w.
    if ( l_scalar_calc ) then
      
      do j = 1, sclr_dim
        do i = 1, ngrdcol
          
          sclrprcp(i,:,j) &
          = pdf_params%mixt_frac(i,:) * ( ( sclr1(i,:,j) - sclrm(i,:,j) ) * pdf_params%rc_1(i,:) ) &
            + ( one - pdf_params%mixt_frac(i,:) ) * ( ( sclr2(i,:,j) - sclrm(i,:,j) ) &
                                                      * pdf_params%rc_2(i,:) ) &
            + pdf_params%mixt_frac(i,:) * corr_sclr_rt_1(i,:,j) * pdf_params%crt_1(i,:) &
              * sqrt( varnce_sclr1(i,:,j) * pdf_params%varnce_rt_1(i,:) ) &
              * pdf_params%cloud_frac_1(i,:) & 
            + ( one - pdf_params%mixt_frac(i,:) ) * corr_sclr_rt_2(i,:,j) * pdf_params%crt_2(i,:) &
              * sqrt( varnce_sclr2(i,:,j) * pdf_params%varnce_rt_2(i,:) ) &
              * pdf_params%cloud_frac_2(i,:) & 
            - pdf_params%mixt_frac(i,:) * corr_sclr_thl_1(i,:,j) * pdf_params%cthl_1(i,:) &
              * sqrt( varnce_sclr1(i,:,j) * pdf_params%varnce_thl_1(i,:) ) &
              * pdf_params%cloud_frac_1(i,:) & 
            - ( one - pdf_params%mixt_frac(i,:) ) * corr_sclr_thl_2(i,:,j) * pdf_params%cthl_2(i,:) &
              * sqrt( varnce_sclr2(i,:,j) * pdf_params%varnce_thl_2(i,:) ) &
              * pdf_params%cloud_frac_2(i,:)

          sclrpthvp(i,:,j) = sclrpthlp(i,:,j) + ep1*thv_ds(i,:)*sclrprtp(i,:,j) &
                           + rc_coef(i,:)*sclrprcp(i,:,j)
                           
        end do
      end do ! i=1, sclr_dim
      
    end if ! l_scalar_calc
    
    

      
    ! Compute variance of liquid water mixing ratio.
    ! This is not needed for closure.  Statistical Analysis only.

#ifndef CLUBB_CAM
      !  if CLUBB is used in CAM we want this variable computed no matter what
      if ( ircp2 > 0 ) then
#endif

    do k = 1,nz
      do i = 1, ngrdcol
        rcp2(i,k) = pdf_params%mixt_frac(i,k) &
                    * ( pdf_params%chi_1(i,k)*pdf_params%rc_1(i,k) &
                        + pdf_params%cloud_frac_1(i,k)*pdf_params%stdev_chi_1(i,k)**2 ) &
                    + ( one-pdf_params%mixt_frac(i,k) ) &
                      * ( pdf_params%chi_2(i,k)*pdf_params%rc_2(i,k) &
                          + pdf_params%cloud_frac_2(i,k)*pdf_params%stdev_chi_2(i,k)**2 ) &
                    - rcm(i,k)**2
        rcp2(i,k) = max( zero_threshold, rcp2(i,k) )
        
      end do
    end do

#ifndef CLUBB_CAM
      !  if CLUBB is used in CAM we want this variable computed no matter what
      end if
#endif

    if ((iiPDF_type == iiPDF_ADG1 .or. iiPDF_type == iiPDF_ADG2 &
                                 .or. iiPDF_type == iiPDF_new_hybrid) &
                                 .and. iw_up_in_cloud > 0) then
                 
      call calc_w_up_in_cloud( nz, ngrdcol, &                                      ! In
                               pdf_params%mixt_frac, &                             ! In
                               pdf_params%cloud_frac_1, pdf_params%cloud_frac_2, & ! In
                               pdf_params%w_1, pdf_params%w_2, &                   ! In
                               pdf_params%varnce_w_1, pdf_params%varnce_w_2, &     ! In
                               w_up_in_cloud )                                     ! Out

    else
      w_up_in_cloud(:,:) = zero
    end if

    if ( clubb_at_least_debug_level( 2 ) ) then
      do i = 1, ngrdcol
          
        call pdf_closure_check( & 
               gr(i), wp4(i,:), wprtp2(i,:), wp2rtp(i,:), wpthlp2(i,:), & ! intent(in)
               wp2thlp(i,:), cloud_frac(i,:), rcm(i,:), wpthvp(i,:), wp2thvp(i,:), &  ! intent(in)
               rtpthvp(i,:), thlpthvp(i,:), wprcp(i,:), wp2rcp(i,:), & ! intent(in)
               rtprcp(i,:), thlprcp(i,:), rcp2(i,:), wprtpthlp(i,:), & ! intent(in)
               pdf_params%crt_1(i,:), pdf_params%crt_2(i,:), & ! intent(in)
               pdf_params%cthl_1(i,:), pdf_params%cthl_2(i,:), & ! intent(in)
               pdf_params, & ! intent(in)
               sclrpthvp(i,:,:), sclrprcp(i,:,:), wpsclrp2(i,:,:), &  ! intent(in)
               wpsclrprtp(i,:,:), wpsclrpthlp(i,:,:), wp2sclrp(i,:,:) ) ! intent(in)

        ! Error Reporting
        ! Joshua Fasching February 2008

        if ( err_code == clubb_fatal_error ) then

          write(fstderr,*) "Error in pdf_closure_new"

          write(fstderr,*) "Intent(in)"

          write(fstderr,*) "p_in_Pa = ", p_in_Pa
          write(fstderr,*) "exner = ", exner
          write(fstderr,*) "thv_ds = ", thv_ds
          write(fstderr,*) "wm = ", wm
          write(fstderr,*) "wp2 = ", wp2
          write(fstderr,*) "wp3 = ", wp3
          write(fstderr,*) "sigma_sqd_w = ", sigma_sqd_w
          write(fstderr,*) "rtm = ", rtm
          write(fstderr,*) "rtp2 = ", rtp2
          write(fstderr,*) "wprtp = ", wprtp
          write(fstderr,*) "thlm = ", thlm
          write(fstderr,*) "thlp2 = ", thlp2
          write(fstderr,*) "wpthlp = ", wpthlp
          write(fstderr,*) "rtpthlp = ", rtpthlp

          if ( sclr_dim > 0 ) then
            write(fstderr,*) "sclrm = ", sclrm
            write(fstderr,*) "wpsclrp = ", wpsclrp
            write(fstderr,*) "sclrp2 = ", sclrp2
            write(fstderr,*) "sclrprtp = ", sclrprtp
            write(fstderr,*) "sclrpthlp = ", sclrpthlp
          end if

          write(fstderr,*) "Intent(out)"

          write(fstderr,*) "wp4 = ", wp4
          write(fstderr,*) "wprtp2 = ", wprtp2
          write(fstderr,*) "wp2rtp = ", wp2rtp
          write(fstderr,*) "wpthlp2 = ", wpthlp2
          write(fstderr,*) "cloud_frac = ", cloud_frac
          write(fstderr,*) "ice_supersat_frac = ", ice_supersat_frac
          write(fstderr,*) "rcm = ", rcm
          write(fstderr,*) "wpthvp = ", wpthvp
          write(fstderr,*) "wp2thvp = ", wp2thvp
          write(fstderr,*) "rtpthvp = ", rtpthvp
          write(fstderr,*) "thlpthvp = ", thlpthvp
          write(fstderr,*) "wprcp = ", wprcp
          write(fstderr,*) "wp2rcp = ", wp2rcp
          write(fstderr,*) "rtprcp = ", rtprcp
          write(fstderr,*) "thlprcp = ", thlprcp
          write(fstderr,*) "rcp2 = ", rcp2
          write(fstderr,*) "wprtpthlp = ", wprtpthlp
          write(fstderr,*) "pdf_params%w_1(i,:) = ", pdf_params%w_1(i,:)
          write(fstderr,*) "pdf_params%w_2(i,:) = ", pdf_params%w_2(i,:)
          write(fstderr,*) "pdf_params%varnce_w_1(i,:) = ", pdf_params%varnce_w_1(i,:)
          write(fstderr,*) "pdf_params%varnce_w_2(i,:) = ", pdf_params%varnce_w_2(i,:)
          write(fstderr,*) "pdf_params%rt_1(i,:) = ", pdf_params%rt_1(i,:)
          write(fstderr,*) "pdf_params%rt_2(i,:) = ", pdf_params%rt_2(i,:)
          write(fstderr,*) "pdf_params%varnce_rt_1(i,:) = ", pdf_params%varnce_rt_1(i,:)
          write(fstderr,*) "pdf_params%varnce_rt_2(i,:) = ", pdf_params%varnce_rt_2(i,:)
          write(fstderr,*) "pdf_params%thl_1(i,:) = ", pdf_params%thl_1(i,:)
          write(fstderr,*) "pdf_params%thl_2(i,:) = ", pdf_params%thl_2(i,:)
          write(fstderr,*) "pdf_params%varnce_thl_1(i,:) = ", pdf_params%varnce_thl_1(i,:)
          write(fstderr,*) "pdf_params%varnce_thl_2(i,:) = ", pdf_params%varnce_thl_2(i,:)
          write(fstderr,*) "pdf_params%corr_w_rt_1(i,:) = ", pdf_params%corr_w_rt_1(i,:)
          write(fstderr,*) "pdf_params%corr_w_rt_2(i,:) = ", pdf_params%corr_w_rt_2(i,:)
          write(fstderr,*) "pdf_params%corr_w_thl_1(i,:) = ", pdf_params%corr_w_thl_1(i,:)
          write(fstderr,*) "pdf_params%corr_w_thl_2(i,:) = ", pdf_params%corr_w_thl_2(i,:)
          write(fstderr,*) "pdf_params%corr_rt_thl_1(i,:) = ", pdf_params%corr_rt_thl_1(i,:)
          write(fstderr,*) "pdf_params%corr_rt_thl_2(i,:) = ", pdf_params%corr_rt_thl_2(i,:)
          write(fstderr,*) "pdf_params%alpha_thl(i,:) = ", pdf_params%alpha_thl(i,:)
          write(fstderr,*) "pdf_params%alpha_rt(i,:) = ", pdf_params%alpha_rt(i,:)
          write(fstderr,*) "pdf_params%crt_1(i,:) = ", pdf_params%crt_1(i,:)
          write(fstderr,*) "pdf_params%crt_2(i,:) = ", pdf_params%crt_2(i,:)
          write(fstderr,*) "pdf_params%cthl_1(i,:) = ", pdf_params%cthl_1(i,:)
          write(fstderr,*) "pdf_params%cthl_2(i,:) = ", pdf_params%cthl_2(i,:)
          write(fstderr,*) "pdf_params%chi_1(i,:) = ", pdf_params%chi_1(i,:)
          write(fstderr,*) "pdf_params%chi_2(i,:) = ", pdf_params%chi_2(i,:)
          write(fstderr,*) "pdf_params%stdev_chi_1(i,:) = ", pdf_params%stdev_chi_1(i,:)
          write(fstderr,*) "pdf_params%stdev_chi_2(i,:) = ", pdf_params%stdev_chi_2(i,:)
          write(fstderr,*) "pdf_params%stdev_eta_1(i,:) = ", pdf_params%stdev_eta_1(i,:)
          write(fstderr,*) "pdf_params%stdev_eta_2(i,:) = ", pdf_params%stdev_eta_2(i,:)
          write(fstderr,*) "pdf_params%covar_chi_eta_1(i,:) = ", &
                           pdf_params%covar_chi_eta_1(i,:)
          write(fstderr,*) "pdf_params%covar_chi_eta_2(i,:) = ", &
                           pdf_params%covar_chi_eta_2(i,:)
          write(fstderr,*) "pdf_params%corr_w_chi_1(i,:) = ", pdf_params%corr_w_chi_1(i,:)
          write(fstderr,*) "pdf_params%corr_w_chi_2(i,:) = ", pdf_params%corr_w_chi_2(i,:)
          write(fstderr,*) "pdf_params%corr_w_eta_1(i,:) = ", pdf_params%corr_w_eta_1(i,:)
          write(fstderr,*) "pdf_params%corr_w_eta_2(i,:) = ", pdf_params%corr_w_eta_2(i,:)
          write(fstderr,*) "pdf_params%corr_chi_eta_1(i,:) = ", &
                           pdf_params%corr_chi_eta_1(i,:)
          write(fstderr,*) "pdf_params%corr_chi_eta_2(i,:) = ", &
                           pdf_params%corr_chi_eta_2(i,:)
          write(fstderr,*) "pdf_params%rsatl_1(i,:) = ", pdf_params%rsatl_1(i,:)
          write(fstderr,*) "pdf_params%rsatl_2(i,:) = ", pdf_params%rsatl_2(i,:)
          write(fstderr,*) "pdf_params%rc_1(i,:) = ", pdf_params%rc_1(i,:)
          write(fstderr,*) "pdf_params%rc_2(i,:) = ", pdf_params%rc_2(i,:)
          write(fstderr,*) "pdf_params%cloud_frac_1(i,:) = ", pdf_params%cloud_frac_1(i,:)
          write(fstderr,*) "pdf_params%cloud_frac_2(i,:) = ", pdf_params%cloud_frac_2(i,:)
          write(fstderr,*) "pdf_params%mixt_frac(i,:) = ", pdf_params%mixt_frac(i,:)
          write(fstderr,*) "pdf_params%ice_supersat_frac_1(i,:) = ", &
                           pdf_params%ice_supersat_frac_1(i,:)
          write(fstderr,*) "pdf_params%ice_supersat_frac_2(i,:) = ", &
                           pdf_params%ice_supersat_frac_2(i,:)

          if ( sclr_dim > 0 )then
            write(fstderr,*) "sclrpthvp = ", sclrpthvp
            write(fstderr,*) "sclrprcp = ", sclrprcp
            write(fstderr,*) "wpsclrp2 = ", wpsclrp2
            write(fstderr,*) "wpsclrprtp = ", wpsclrprtp
            write(fstderr,*) "wpsclrpthlp = ", wpsclrpthlp
            write(fstderr,*) "wp2sclrp = ", wp2sclrp
          end if

        end if ! Fatal error

        ! Error check pdf parameters and moments to ensure consistency
        if ( iiPDF_type == iiPDF_3D_Luhar ) then

          ! Means
          wm_clubb_pdf(i,:) = pdf_params%mixt_frac(i,:) * pdf_params%w_1(i,:) &
                         + ( one - pdf_params%mixt_frac(i,:) ) * pdf_params%w_2(i,:)

          do k = 1, nz, 1
             if ( abs( ( wm_clubb_pdf(i,k) - wm(i,k) ) &
                       / max( wm(i,k), eps ) ) > .05_core_rknd ) then
                write(fstderr,*) "wm error at thlm = ", thlm(i,k), &
                                 ( ( wm_clubb_pdf(i,k) - wm(i,k) ) &
                                   / max( wm(i,k), eps ) )
             end if
          end do ! k = 1, nz, 1

          rtm_clubb_pdf(i,:) = pdf_params%mixt_frac(i,:) * pdf_params%rt_1(i,:) &
                          + ( one - pdf_params%mixt_frac(i,:) ) * pdf_params%rt_2(i,:)

          do k = 1, nz, 1
             if ( abs( ( rtm_clubb_pdf(i,k) - rtm(i,k) ) &
                       / max( rtm(i,k), eps ) ) > .05_core_rknd ) then
                write(fstderr,*) "rtm error at thlm = ", thlm(i,k), &
                                 ( ( rtm_clubb_pdf(i,k) - rtm(i,k) ) &
                                   / max( rtm(i,k), eps ) )
             end if
          end do ! k = 1, nz, 1

          thlm_clubb_pdf(i,:) = pdf_params%mixt_frac(i,:) * pdf_params%thl_1(i,:) &
                           + ( one - pdf_params%mixt_frac(i,:) ) * pdf_params%thl_2(i,:)

          do k = 1, nz, 1
             if ( abs( ( thlm_clubb_pdf(i,k) - thlm(i,k) ) / thlm(i,k) ) &
                  > .05_core_rknd ) then
                write(fstderr,*) "thlm error at thlm = ", thlm(i,k), &
                                 ( ( thlm_clubb_pdf(i,k) - thlm(i,k) ) / thlm(i,k) )
             end if
          end do ! k = 1, nz, 1

          ! Variances
          wp2_clubb_pdf(i,:) = pdf_params%mixt_frac(i,:) &
                          * ( ( pdf_params%w_1(i,:) - wm(i,:) )**2 + pdf_params%varnce_w_1(i,:) ) &
                          + ( one - pdf_params%mixt_frac(i,:) ) &
                            * ( ( pdf_params%w_2(i,:) - wm(i,:) )**2 + pdf_params%varnce_w_2(i,:) )

          do k = 1, nz, 1
             if ( wp2(i,k) > w_tol**2 ) then
                if ( abs( ( wp2_clubb_pdf(i,k) - wp2(i,k) ) / wp2(i,k) ) &
                     > .05_core_rknd ) then
                   write(fstderr,*) "wp2 error at thlm = ", thlm(i,k), &
                                    ( ( wp2_clubb_pdf(i,k) - wp2(i,k) ) / wp2(i,k) )
                end if
             end if
          end do ! k = 1, nz, 1

          rtp2_clubb_pdf(i,:) &
          = pdf_params%mixt_frac(i,:) &
            * ( ( pdf_params%rt_1(i,:) - rtm(i,:) )**2 + pdf_params%varnce_rt_1(i,:) ) &
            + ( one - pdf_params%mixt_frac(i,:) ) &
              * ( ( pdf_params%rt_2(i,:) - rtm(i,:) )**2 + pdf_params%varnce_rt_2(i,:) )

          do k = 1, nz, 1
             if ( rtp2(i,k) > rt_tol**2 ) then
                if ( abs( ( rtp2_clubb_pdf(i,k) - rtp2(i,k) ) / rtp2(i,k) ) &
                     > .05_core_rknd ) then
                   write(fstderr,*) "rtp2 error at thlm = ", thlm(i,k), &
                   "Error = ", ( ( rtp2_clubb_pdf(i,k) - rtp2(i,k) ) / rtp2(i,k) )
                end if
             end if
          end do ! k = 1, nz, 1

          thlp2_clubb_pdf(i,:) &
          = pdf_params%mixt_frac(i,:) &
            * ( ( pdf_params%thl_1(i,:) - thlm(i,:) )**2 + pdf_params%varnce_thl_1(i,:) ) &
            + ( one - pdf_params%mixt_frac(i,:) ) &
              * ( ( pdf_params%thl_2(i,:) - thlm(i,:) )**2 + pdf_params%varnce_thl_2(i,:) )

          do k = 1, nz, 1
             if( thlp2(i,k) > thl_tol**2 ) then
                if ( abs( ( thlp2_clubb_pdf(i,k) - thlp2(i,k) ) / thlp2(i,k) ) &
                     > .05_core_rknd ) then
                   write(fstderr,*) "thlp2 error at thlm = ", thlm(i,k), &
                   "Error = ", ( ( thlp2_clubb_pdf(i,k) - thlp2(i,k) ) / thlp2(i,k) )
                end if
             end if
          end do ! k = 1, nz, 1

          ! Third order moments
          wp3_clubb_pdf(i,:) &
          = pdf_params%mixt_frac(i,:) * ( pdf_params%w_1(i,:) - wm(i,:) ) &
            * ( ( pdf_params%w_1(i,:) - wm(i,:) )**2 + three * pdf_params%varnce_w_1(i,:) ) &
            + ( one - pdf_params%mixt_frac(i,:) ) * ( pdf_params%w_2(i,:) - wm(i,:) ) &
              * ( ( pdf_params%w_2(i,:) - wm(i,:) )**2 + three * pdf_params%varnce_w_2(i,:) )

          rtp3_clubb_pdf(i,:) &
          = pdf_params%mixt_frac(i,:) * ( pdf_params%rt_1(i,:) - rtm(i,:) ) &
            * ( ( pdf_params%rt_1(i,:) - rtm(i,:) )**2 + three * pdf_params%varnce_rt_1(i,:) ) &
            + ( one - pdf_params%mixt_frac(i,:) ) * ( pdf_params%rt_2(i,:) - rtm(i,:) ) &
              * ( ( pdf_params%rt_2(i,:) - rtm(i,:) )**2 + three * pdf_params%varnce_rt_2(i,:) )

          thlp3_clubb_pdf(i,:) &
          = pdf_params%mixt_frac(i,:) * ( pdf_params%thl_1(i,:) - thlm(i,:) ) &
            * ( ( pdf_params%thl_1(i,:) - thlm(i,:) )**2 + three * pdf_params%varnce_thl_1(i,:) ) &
            + ( one - pdf_params%mixt_frac(i,:) ) * ( pdf_params%thl_2(i,:) - thlm(i,:) ) &
              * ( ( pdf_params%thl_2(i,:) - thlm(i,:) )**2 + three * pdf_params%varnce_thl_2(i,:) )

          ! Skewness
          Skw_denom_coef = clubb_params(iSkw_denom_coef)

          Skw_clubb_pdf(i,:) &
          = wp3_clubb_pdf(i,:) &
            / ( wp2_clubb_pdf(i,:) + Skw_denom_coef * w_tol**2 )**1.5_core_rknd

          do k = 1, nz, 1
             if ( Skw(i,k) > .05_core_rknd ) then
                if( abs( ( Skw_clubb_pdf(i,k) - Skw(i,k) ) / Skw(i,k) ) &
                    > .25_core_rknd ) then
                   write(fstderr,*) "Skw error at thlm = ", thlm(i,k), &
                   "Error = ", ( ( Skw_clubb_pdf(i,k) - Skw(i,k) ) / Skw(i,k) ), &
                   Skw_clubb_pdf(i,k), Skw(i,k)
                end if
             end if
          end do ! k = 1, nz, 1

          Skrt_clubb_pdf(i,:) &
          = rtp3_clubb_pdf(i,:) &
            / ( rtp2_clubb_pdf(i,:) + Skw_denom_coef * rt_tol**2 )**1.5_core_rknd

          do k = 1, nz, 1
             if ( Skrt(i,k) > .05_core_rknd ) then
                if( abs( ( Skrt_clubb_pdf(i,k) - Skrt(i,k) ) / Skrt(i,k) ) &
                    > .25_core_rknd ) then
                   write(fstderr,*) "Skrt error at thlm = ", thlm(i,k), &
                   "Error = ", ( ( Skrt_clubb_pdf(i,k) - Skrt(i,k) ) / Skrt(i,k) ), &
                   Skrt_clubb_pdf(i,k), Skrt(i,k)
                end if
             end if
          end do ! k = 1, nz, 1

          Skthl_clubb_pdf(i,:) &
          = thlp3_clubb_pdf(i,:) &
            / ( thlp2_clubb_pdf(i,:) + Skw_denom_coef * thl_tol**2 )**1.5_core_rknd

          do k = 1, nz, 1
             if ( Skthl(i,k) > .05_core_rknd ) then
                if ( abs( ( Skthl_clubb_pdf(i,k) - Skthl(i,k) ) / Skthl(i,k) ) &
                     > .25_core_rknd ) then
                   write(fstderr,*) "Skthl error at thlm = ", thlm(i,k), &
                   "Error = ", ( ( Skthl_clubb_pdf(i,k) - Skthl(i,k) ) / Skthl(i,k) ), &
                   Skthl_clubb_pdf(i,k), Skthl(i,k)
                end if
             end if
          end do ! k = 1, nz, 1

        end if ! iiPDF_type == iiPDF_3D_Luhar
        
      end do

    end if ! clubb_at_least_debug_level

    return

  end subroutine pdf_closure

  !===============================================================================================
  subroutine transform_pdf_chi_eta_component( nz, ngrdcol, &
                                              tl, rsatl, rt, exner,     & ! In
                                              varnce_thl, varnce_rt,    & ! In
                                              corr_rt_thl, chi,         & ! In
                                              crt, cthl,                & ! Out
                                              stdev_chi, stdev_eta,     & ! Out
                                              covar_chi_eta,            & ! Out
                                              corr_chi_eta )              ! Out

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use constants_clubb, only: &
        zero, one, two, &
        ep, Lv, Rd, Cp, &
        chi_tol, &
        eta_tol, &
        max_mag_correlation

    implicit none

    integer, intent(in) :: &
      ngrdcol,  & ! Number of grid columns
      nz          ! Number of vertical level

    ! ----------- Input Variables -----------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
      tl, &
      rsatl, &
      rt, &
      varnce_thl, &
      varnce_rt, &
      corr_rt_thl, &
      exner
    
    ! ----------- Output Variables -----------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) ::  &
      chi, &            ! s from Lewellen and Yoh 1993 (LY) eqn. 1
      crt, &            ! Coefficients for s'
      cthl, &           ! Coefficients for s'
      stdev_chi, &      ! Standard deviation of chi for each component.
      stdev_eta, &      ! Standard deviation of eta for each component.
      covar_chi_eta, &  ! Covariance of chi and eta for each component.
      corr_chi_eta      ! Correlation of chi and eta for each component.

    ! ----------- Local Variables -----------
    real( kind = core_rknd ) :: &
      varnce_rt_term, &
      corr_rt_thl_term, &
      varnce_thl_term, &
      varnce_chi, &
      varnce_eta, &
      beta, &
      invrs_beta_rsatl_p1

    real( kind = core_rknd ), parameter :: &
      chi_tol_sqd = chi_tol**2, &
      eta_tol_sqd = eta_tol**2, &
      Cp_on_Lv = Cp / Lv

    ! Loop variable
    integer :: k, i

    ! ----------- Begin Code -----------

    do k = 1, nz
      do i = 1, ngrdcol

        ! SD's beta (eqn. 8)
        beta = ep * Lv**2 / ( Rd * Cp * tl(i,k)**2 )

        invrs_beta_rsatl_p1 = one / ( one + beta * rsatl(i,k) )

        ! s from Lewellen and Yoh 1993 (LY) eqn. 1
        chi(i,k) = ( rt(i,k) - rsatl(i,k) ) * invrs_beta_rsatl_p1

        ! For each normal distribution in the sum of two normal distributions,
        ! s' = crt * rt'  +  cthl * thl';
        ! therefore, x's' = crt * x'rt'  +  cthl * x'thl'.
        ! Larson et al. May, 2001.
        crt(i,k)  = invrs_beta_rsatl_p1
        cthl(i,k) = ( one + beta * rt(i,k) ) * invrs_beta_rsatl_p1**2 &
                    * Cp_on_Lv * beta * rsatl(i,k) * exner(i,k)
                  
      end do
    end do

    ! Calculate covariance, correlation, and standard deviation of 
    ! chi and eta for each component
    ! Include subplume correlation of qt, thl
    do k = 1, nz
      do i = 1, ngrdcol
       
        varnce_rt_term = crt(i,k)**2 * varnce_rt(i,k)
        varnce_thl_term = cthl(i,k)**2 * varnce_thl(i,k)

        covar_chi_eta(i,k) = varnce_rt_term - varnce_thl_term

        corr_rt_thl_term = two * corr_rt_thl(i,k) * crt(i,k) * cthl(i,k) &
                           * sqrt( varnce_rt(i,k) * varnce_thl(i,k) )

        varnce_chi = varnce_rt_term - corr_rt_thl_term + varnce_thl_term
        varnce_eta = varnce_rt_term + corr_rt_thl_term + varnce_thl_term

        ! We need to introduce a threshold value for the variance of chi and eta
        if ( varnce_chi < chi_tol_sqd .or. varnce_eta < eta_tol_sqd ) then

            if ( varnce_chi < chi_tol_sqd ) then
                stdev_chi(i,k) = zero  ! Treat chi as a delta function
            else
                stdev_chi(i,k) = sqrt( varnce_chi )
            end if

            if ( varnce_eta < eta_tol_sqd ) then
                stdev_eta(i,k) = zero  ! Treat eta as a delta function
            else
                stdev_eta(i,k) = sqrt( varnce_eta )
            end if

            corr_chi_eta(i,k) = zero

        else

            stdev_chi(i,k) = sqrt( varnce_chi )
            stdev_eta(i,k) = sqrt( varnce_eta )

            corr_chi_eta(i,k) = covar_chi_eta(i,k) / ( stdev_chi(i,k) * stdev_eta(i,k) )
            corr_chi_eta(i,k) = min( max_mag_correlation, &
                                   max( -max_mag_correlation, corr_chi_eta(i,k) ) )

        end if

      end do
    end do

  end subroutine transform_pdf_chi_eta_component
  
  !=============================================================================
  subroutine calc_wp4_pdf( nz, ngrdcol, &
                           wm, w_1, w_2, &
                           varnce_w_1, varnce_w_2,    &
                           mixt_frac, &
                           wp4 )

    ! Description:
    ! Calculates <w'^4> by integrating over the PDF of w.  The integral is:
    !
    ! <w'^4> = INT(-inf:inf) ( w - <w> )^4 P(w) dw;
    !
    ! where <w> is the overall mean of w and P(w) is a two-component normal
    ! distribution of w.  The integrated equation is:
    !
    ! <w'^4> = mixt_frac * ( 3 * sigma_w_1^4
    !                        + 6 * ( mu_w_1 - <w> )^2 * sigma_w_1^2
    !                        + ( mu_w_1 - <w> )^4 )
    !          + ( 1 - mixt_frac ) * ( 3 * sigma_w_2^4
    !                                  + 6 * ( mu_w_2 - <w> )^2 * sigma_w_2^2
    !                                  + ( mu_w_2 - <w> )^4 );
    !
    ! where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
    ! of w in the 2nd PDF component, sigma_w_1 is the standard deviation of w in
    ! the 1st PDF component, sigma_w_2 is the standard deviation of w in the 2nd
    ! PDF component, and mixt_frac is the mixture fraction, which is the weight
    ! of the 1st PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        six,   & ! Variable(s)
        three, &
        one

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    integer, intent(in) :: &
      ngrdcol,  & ! Number of grid columns
      nz          ! Number of vertical level

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
      wm,         & ! Mean of w (overall)                           [m/s]
      w_1,        & ! Mean of w (1st PDF component)                 [m/s]
      w_2,        & ! Mean of w (2nd PDF component)                 [m/s]
      varnce_w_1, & ! Variance of w (1st PDF component)             [m^2/s^2]
      varnce_w_2, & ! Variance of w (2nd PDF component)             [m^2/s^2]
      mixt_frac     ! Mixture fraction                              [-]

    ! Output Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz) ::  & 
      wp4    ! <w'^4>                   [m^4/s^4]
      
    ! Local Variables
    integer :: i, k

    
    do k = 1, nz
      do i = 1, ngrdcol

        ! Calculate <w'^4> by integrating over the PDF.
        wp4(i,k) = mixt_frac(i,k) * ( three * varnce_w_1(i,k)**2 &
                            + six * ( ( w_1(i,k) - wm(i,k) )**2 ) * varnce_w_1(i,k) &
                            + ( w_1(i,k) - wm(i,k) )**4 ) & 
                   + ( one - mixt_frac(i,k) ) * ( three * varnce_w_2(i,k)**2 &
                                        + six * ( (w_2(i,k) - wm(i,k) )**2 )*varnce_w_2(i,k) &
                                        + ( w_2(i,k) - wm(i,k) )**4 )
      end do
    end do

    return

  end subroutine calc_wp4_pdf

  !=============================================================================
  subroutine calc_wp2xp2_pdf( nz, ngrdcol,             &
                              wm, xm, w_1,             &
                              w_2, x_1, x_2,           &
                              varnce_w_1, varnce_w_2,  &
                              varnce_x_1, varnce_x_2,  &
                              corr_w_x_1, corr_w_x_2,  &
                              mixt_frac, &
                              wp2xp2 )

    ! Description:
    ! Calculates <w'^2x'^2> by integrating over the PDF of w and x.  The
    ! integral
    ! is:
    !
    ! <w'^2x'^2>
    ! = INT(-inf:inf) INT(-inf:inf) ( w - <w> )^2 ( x - <x> )^2 P(w,x) dx dw;
    !
    ! where <w> is the overall mean of w, <x> is the overall mean of x, and
    ! P(w,x) is a two-component bivariate normal distribution of w and x.  The
    ! integrated equation is:
    !
    ! <w'^2x'^2>
    !   = mixt_frac
    !      * ( ( mu_w_1 - <w> )**2 * ( ( mu_x_1 - <x> )**2 + sigma_x_1^2 )
    !      + four * corr_w_x_1 * sigma_w_1 * sigma_x_1 * ( mu_x_1 - <x> ) * (
    !      mu_w_1 - <w> )
    !      + ( ( mu_x_1 - <x> )**2 + ( 1 + 2*corr_w_x_1**2 ) * sigma_x_1^2 ) *
    !      sigma_w_1^2 )
    !    + ( one - mixt_frac )
    !      * ( ( mu_w_2 - <w> )**2 * ( ( mu_x_2 - <x> )**2 + sigma_x_2^2 )
    !      + four * corr_w_x_2 * sigma_w_2 * sigma_x_2 * ( mu_x_2 - <x> ) * (
    !      mu_w_2 - <w> )
    !      + ( ( mu_x_2 - <x> )**2 + ( 1 + 2*corr_w_x_2**2 ) * sigma_x_2^2 ) *
    !      sigma_w_2^2 )
    !
    ! where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
    ! of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
    ! component, mu_x_2 is the mean of x in the 2nd PDF component, sigma_w_1 is
    ! the standard deviation of w in the 1st PDF component, sigma_w_2 is the
    ! standard deviation of w in the 2nd PDF component, sigma_x_1 is the
    ! standard deviation of x in the 1st PDF component, sigma_x_2 is the
    ! standard deviation of x in the 2nd PDF component, corr_w_x_1 is the
    ! correlation of w and x in the 1st PDF component, corr_w_x_2 is the
    ! correlation of w and x in the 2nd PDF component, and mixt_frac is the
    ! mixture fraction, which is the weight of the 1st PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,   & ! Variable(s)
        four

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    integer, intent(in) :: &
      ngrdcol,  & ! Number of grid columns
      nz          ! Number of vertical level

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
      wm,         & ! Mean of w (overall)                       [m/s]
      xm,         & ! Mean of x (overall)                       [units vary]
      w_1,        & ! Mean of w (1st PDF component)             [m/s]
      w_2,        & ! Mean of w (2nd PDF component)             [m/s]
      x_1,        & ! Mean of x (1st PDF component)             [units vary]
      x_2,        & ! Mean of x (2nd PDF component)             [units vary]
      varnce_w_1, & ! Variance of w (1st PDF component)         [m^2/s^2]
      varnce_w_2, & ! Variance of w (2nd PDF component)         [m^2/s^2]
      varnce_x_1, & ! Variance of x (1st PDF component)         [(units vary)^2]
      varnce_x_2, & ! Variance of x (2nd PDF component)         [(units vary)^2]
      corr_w_x_1, & ! Correlation of w and x (1st PDF comp.)    [-]
      corr_w_x_2, & ! Correlation of w and x (2nd PDF comp.)    [-]
      mixt_frac     ! Mixture fraction                          [-]

    ! Output Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) ::  &
      wp2xp2        ! <w'^2x'^2>                   [m^2/s^2 (units vary)^2]
    
    ! Local Variable
    integer :: i, k


    do k = 1, nz
      do i = 1, ngrdcol

        ! Calculate <w'x'^2> by integrating over the PDF.
        wp2xp2(i,k) = mixt_frac(i,k) &
                * ( ( w_1(i,k) - wm(i,k) )**2 * ( ( x_1(i,k) - xm(i,k) )**2 + varnce_x_1(i,k) ) &
                + four * corr_w_x_1(i,k) * sqrt( varnce_w_1(i,k) * varnce_x_1(i,k) ) &
                                         * ( x_1(i,k) - xm(i,k) ) * ( w_1(i,k) - wm(i,k) ) &
                + ( ( x_1(i,k) - xm(i,k) )**2 &
                + ( 1 + 2*corr_w_x_1(i,k)**2 ) * varnce_x_1(i,k) ) * varnce_w_1(i,k) ) &
                + ( one - mixt_frac(i,k) ) &
                * ( ( w_2(i,k) - wm(i,k) )**2 * ( ( x_2(i,k) - xm(i,k) )**2 + varnce_x_2(i,k) ) &
                + four * corr_w_x_2(i,k) * sqrt( varnce_w_2(i,k) * varnce_x_2(i,k) ) &
                                         * ( x_2(i,k) - xm(i,k) ) * ( w_2(i,k) - wm(i,k) ) &
                + ( ( x_2(i,k) - xm(i,k) )**2 &
                + ( 1 + 2*corr_w_x_2(i,k)**2 ) * varnce_x_2(i,k) ) * varnce_w_2(i,k) )
      end do
    end do

    return

  end subroutine calc_wp2xp2_pdf

  !=============================================================================
  subroutine calc_wp2xp_pdf( nz, ngrdcol,             &
                             wm, xm, w_1, w_2,        &
                             x_1, x_2,                &
                             varnce_w_1, varnce_w_2,  &
                             varnce_x_1, varnce_x_2,  &
                             corr_w_x_1, corr_w_x_2,  &
                             mixt_frac, &
                             wp2xp ) 

    ! Description:
    ! Calculates <w'^2 x'> by integrating over the PDF of w and x.  The integral
    ! is:
    !
    ! <w'^2 x'>
    ! = INT(-inf:inf) INT(-inf:inf) ( w - <w> )^2 ( x - <x> ) P(w,x) dx dw;
    !
    ! where <w> is the overall mean of w, <x> is the overall mean of x, and
    ! P(w,x) is a two-component bivariate normal distribution of w and x.  The
    ! integrated equation is:
    !
    ! <w'^2 x'>
    ! = mixt_frac * ( ( mu_x_1 - <x> ) * ( ( mu_w_1 - <w> )^2 + sigma_w_1^2 )
    !                 + 2 * corr_w_x_1 * sigma_w_1 * sigma_x_1
    !                   * ( mu_w_1 - <w> ) )
    !   + ( 1 - mixt_frac ) * ( ( mu_x_2 - <x> )
    !                           * ( ( mu_w_2 - <w> )^2 + sigma_w_2^2 )
    !                           + 2 * corr_w_x_2 * sigma_w_2 * sigma_x_2
    !                             * ( mu_w_2 - <w> ) );
    !
    ! where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
    ! of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
    ! component, mu_x_2 is the mean of x in the 2nd PDF component, sigma_w_1 is
    ! the standard deviation of w in the 1st PDF component, sigma_w_2 is the
    ! standard deviation of w in the 2nd PDF component, sigma_x_1 is the
    ! standard deviation of x in the 1st PDF component, sigma_x_2 is the
    ! standard deviation of x in the 2nd PDF component, corr_w_x_1 is the
    ! correlation of w and x in the 1st PDF component, corr_w_x_2 is the
    ! correlation of w and x in the 2nd PDF component, and mixt_frac is the
    ! mixture fraction, which is the weight of the 1st PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two,   & ! Variable(s)
        one

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    integer, intent(in) :: &
      ngrdcol,  & ! Number of grid columns
      nz          ! Number of vertical level

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
      wm,         & ! Mean of w (overall)                       [m/s]
      xm,         & ! Mean of x (overall)                       [units vary]
      w_1,        & ! Mean of w (1st PDF component)             [m/s]
      w_2,        & ! Mean of w (2nd PDF component)             [m/s]
      x_1,        & ! Mean of x (1st PDF component)             [units vary]
      x_2,        & ! Mean of x (2nd PDF component)             [units vary]
      varnce_w_1, & ! Variance of w (1st PDF component)         [m^2/s^2]
      varnce_w_2, & ! Variance of w (2nd PDF component)         [m^2/s^2]
      varnce_x_1, & ! Variance of x (1st PDF component)         [(units vary)^2]
      varnce_x_2, & ! Variance of x (2nd PDF component)         [(units vary)^2]
      corr_w_x_1, & ! Correlation of w and x (1st PDF comp.)    [-]
      corr_w_x_2, & ! Correlation of w and x (2nd PDF comp.)    [-]
      mixt_frac     ! Mixture fraction                          [-]

    ! Output Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) ::  & 
      wp2xp    ! <w'^2 x'>                   [m^2/s^2 (units vary)]

    ! Local Variables
    integer :: i, k


    ! Calculate <w'^2 x'> by integrating over the PDF.
    do k = 1, nz
      do i = 1, ngrdcol
        
        wp2xp(i,k)  = mixt_frac(i,k) &
                   * ( ( ( w_1(i,k) - wm(i,k) )**2 + varnce_w_1(i,k) ) * ( x_1(i,k) - xm(i,k) ) &
                       + two * corr_w_x_1(i,k) * sqrt( varnce_w_1(i,k) * varnce_x_1(i,k) ) &
                         * ( w_1(i,k) - wm(i,k) ) ) &
                   + ( one - mixt_frac(i,k) ) &
                     * ( ( ( w_2(i,k) - wm(i,k) )**2 + varnce_w_2(i,k) ) * ( x_2(i,k) - xm(i,k) ) &
                         + two * corr_w_x_2(i,k) * sqrt( varnce_w_2(i,k) * varnce_x_2(i,k) ) &
                           * ( w_2(i,k) - wm(i,k) ) )
      end do
    end do

    return

  end subroutine calc_wp2xp_pdf

  !=============================================================================
  subroutine calc_wpxp2_pdf( nz, ngrdcol,             &
                             wm, xm, w_1,             &
                             w_2, x_1, x_2,           &
                             varnce_w_1, varnce_w_2,  &
                             varnce_x_1, varnce_x_2,  &
                             corr_w_x_1, corr_w_x_2,  &
                             mixt_frac, &
                             wpxp2 )

    ! Description:
    ! Calculates <w'x'^2> by integrating over the PDF of w and x.  The integral
    ! is:
    !
    ! <w'x'^2>
    ! = INT(-inf:inf) INT(-inf:inf) ( w - <w> ) ( x - <x> )^2 P(w,x) dx dw;
    !
    ! where <w> is the overall mean of w, <x> is the overall mean of x, and
    ! P(w,x) is a two-component bivariate normal distribution of w and x.  The
    ! integrated equation is:
    !
    ! <w'x'^2>
    ! = mixt_frac * ( ( mu_w_1 - <w> ) * ( ( mu_x_1 - <x> )^2 + sigma_x_1^2 )
    !                 + 2 * corr_w_x_1 * sigma_w_1 * sigma_x_1
    !                   * ( mu_x_1 - <x> ) )
    !   + ( 1 - mixt_frac ) * ( ( mu_w_2 - <w> )
    !                           * ( ( mu_x_2 - <x> )^2 + sigma_x_2^2 )
    !                           + 2 * corr_w_x_2 * sigma_w_2 * sigma_x_2
    !                             * ( mu_x_2 - <x> ) );
    !
    ! where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
    ! of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
    ! component, mu_x_2 is the mean of x in the 2nd PDF component, sigma_w_1 is
    ! the standard deviation of w in the 1st PDF component, sigma_w_2 is the
    ! standard deviation of w in the 2nd PDF component, sigma_x_1 is the
    ! standard deviation of x in the 1st PDF component, sigma_x_2 is the
    ! standard deviation of x in the 2nd PDF component, corr_w_x_1 is the
    ! correlation of w and x in the 1st PDF component, corr_w_x_2 is the
    ! correlation of w and x in the 2nd PDF component, and mixt_frac is the
    ! mixture fraction, which is the weight of the 1st PDF component.

    ! References:
    !-----------------------------------------------------------------------
    
    use constants_clubb, only: &
        two,   & ! Variable(s)
        one

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    integer, intent(in) :: &
      ngrdcol,  & ! Number of grid columns
      nz          ! Number of vertical level

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
      wm,         & ! Mean of w (overall)                       [m/s]
      xm,         & ! Mean of x (overall)                       [units vary]
      w_1,        & ! Mean of w (1st PDF component)             [m/s]
      w_2,        & ! Mean of w (2nd PDF component)             [m/s]
      x_1,        & ! Mean of x (1st PDF component)             [units vary]
      x_2,        & ! Mean of x (2nd PDF component)             [units vary]
      varnce_w_1, & ! Variance of w (1st PDF component)         [m^2/s^2]
      varnce_w_2, & ! Variance of w (2nd PDF component)         [m^2/s^2]
      varnce_x_1, & ! Variance of x (1st PDF component)         [(units vary)^2]
      varnce_x_2, & ! Variance of x (2nd PDF component)         [(units vary)^2]
      corr_w_x_1, & ! Correlation of w and x (1st PDF comp.)    [-]
      corr_w_x_2, & ! Correlation of w and x (2nd PDF comp.)    [-]
      mixt_frac     ! Mixture fraction                          [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) ::  & 
      wpxp2    ! <w'x'^2>                   [m/s (units vary)^2]
      
    ! Local Variables
    integer :: i, k

    do k = 1, nz
      do i = 1, ngrdcol

        ! Calculate <w'x'^2> by integrating over the PDF.
        wpxp2(i,k) = mixt_frac(i,k) &
                * ( ( w_1(i,k) - wm(i,k) ) * ( ( x_1(i,k) - xm(i,k) )**2 + varnce_x_1(i,k) ) &
                    + two * corr_w_x_1(i,k) * sqrt( varnce_w_1(i,k) * varnce_x_1(i,k) ) &
                      * ( x_1(i,k) - xm(i,k) ) ) &
                + ( one - mixt_frac(i,k) ) &
                  * ( ( w_2(i,k) - wm(i,k) ) * ( ( x_2(i,k) - xm(i,k) )**2 + varnce_x_2(i,k) ) &
                      + two * corr_w_x_2(i,k) * sqrt( varnce_w_2(i,k) * varnce_x_2(i,k) ) &
                        * ( x_2(i,k) - xm(i,k) ) )
      end do
    end do

    return

  end subroutine calc_wpxp2_pdf

  !=============================================================================
  subroutine calc_wpxpyp_pdf( nz, ngrdcol, &
                              wm, xm, ym, w_1, w_2,   &
                              x_1, x_2,               &
                              y_1, y_2,               &
                              varnce_w_1, varnce_w_2, &
                              varnce_x_1, varnce_x_2, &
                              varnce_y_1, varnce_y_2, &
                              corr_w_x_1, corr_w_x_2, &
                              corr_w_y_1, corr_w_y_2, &
                              corr_x_y_1, corr_x_y_2, &
                              mixt_frac, &
                              wpxpyp )

    ! Description:
    ! Calculates <w'x'y'> by integrating over the PDF of w, x, and y.  The
    ! integral is:
    !
    ! <w'x'y'>
    ! = INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    !   ( w - <w> ) ( x - <x> ) ( y - <y> ) P(w,x,y) dy dx dw;
    !
    ! where <w> is the overall mean of w, <x> is the overall mean of x, <y> is
    ! the overall mean of y, and P(w,x,y) is a two-component trivariate normal
    ! distribution of w, x, and y.  The integrated equation is:
    !
    ! <w'x'y'>
    ! = mixt_frac 
    !   * ( ( mu_w_1 - <w> ) * ( mu_x_1 - <x> ) * ( mu_y_1 - <y> )
    !       + corr_x_y_1 * sigma_x_1 * sigma_y_1 * ( mu_w_1 - <w> )
    !       + corr_w_y_1 * sigma_w_1 * sigma_y_1 * ( mu_x_1 - <x> )
    !       + corr_w_x_1 * sigma_w_1 * sigma_x_1 * ( mu_y_1 - <y> ) )
    !   + ( 1 - mixt_frac )
    !     * ( ( mu_w_2 - <w> ) * ( mu_x_2 - <x> ) * ( mu_y_2 - <y> )
    !         + corr_x_y_2 * sigma_x_2 * sigma_y_2 * ( mu_w_2 - <w> )
    !         + corr_w_y_2 * sigma_w_2 * sigma_y_2 * ( mu_x_2 - <x> )
    !         + corr_w_x_2 * sigma_w_2 * sigma_x_2 * ( mu_y_2 - <y> ) );
    !
    ! where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
    ! of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
    ! component, mu_x_2 is the mean of x in the 2nd PDF component, mu_y_1 is the
    ! mean of y in the 1st PDF component, mu_y_2 is the mean of y in the 2nd PDF
    ! component, sigma_w_1 is the standard deviation of w in the 1st PDF
    ! component, sigma_w_2 is the standard deviation of w in the 2nd PDF
    ! component, sigma_x_1 is the standard deviation of x in the 1st PDF
    ! component, sigma_x_2 is the standard deviation of x in the 2nd PDF
    ! component, sigma_y_1 is the standard deviation of y in the 1st PDF
    ! component, sigma_y_2 is the standard deviation of y in the 2nd PDF
    ! component, corr_w_x_1 is the correlation of w and x in the 1st PDF
    ! component, corr_w_x_2 is the correlation of w and x in the 2nd PDF
    ! component, corr_w_y_1 is the correlation of w and y in the 1st PDF
    ! component, corr_w_y_2 is the correlation of w and y in the 2nd PDF
    ! component, corr_x_y_1 is the correlation of x and y in the 1st PDF
    ! component, corr_x_y_2 is the correlation of x and y in the 2nd PDF
    ! component, and mixt_frac is the mixture fraction, which is the weight of
    ! the 1st PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid ! Type

    use constants_clubb, only: &
        one    ! Variable(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    integer, intent(in) :: &
      ngrdcol,  & ! Number of grid columns
      nz          ! Number of vertical level

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
      wm,         & ! Mean of w (overall)                          [m/s]
      xm,         & ! Mean of x (overall)                          [x units]
      ym,         & ! Mean of y (overall)                          [y units]
      w_1,        & ! Mean of w (1st PDF component)                [m/s]
      w_2,        & ! Mean of w (2nd PDF component)                [m/s]
      x_1,        & ! Mean of x (1st PDF component)                [x units]
      x_2,        & ! Mean of x (2nd PDF component)                [x units]
      y_1,        & ! Mean of y (1st PDF component)                [y units]
      y_2,        & ! Mean of y (2nd PDF component)                [y units]
      varnce_w_1, & ! Variance of w (1st PDF component)            [m^2/s^2]
      varnce_w_2, & ! Variance of w (2nd PDF component)            [m^2/s^2]
      varnce_x_1, & ! Variance of x (1st PDF component)            [(x units)^2]
      varnce_x_2, & ! Variance of x (2nd PDF component)            [(x units)^2]
      varnce_y_1, & ! Variance of y (1st PDF component)            [(y units)^2]
      varnce_y_2, & ! Variance of y (2nd PDF component)            [(y units)^2]
      corr_w_x_1, & ! Correlation of w and x (1st PDF component)   [-]
      corr_w_x_2, & ! Correlation of w and x (2nd PDF component)   [-]
      corr_w_y_1, & ! Correlation of w and y (1st PDF component)   [-]
      corr_w_y_2, & ! Correlation of w and y (2nd PDF component)   [-]
      corr_x_y_1, & ! Correlation of x and y (1st PDF component)   [-]
      corr_x_y_2, & ! Correlation of x and y (2nd PDF component)   [-]
      mixt_frac     ! Mixture fraction                             [-]

    ! Output Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz) ::  & 
      wpxpyp    ! <w'x'y'>                   [m/s (units vary)]

    ! Local Variables
    integer :: i, k
    

    ! Calculate <w'x'y'> by integrating over the PDF.
    do k = 1, nz
      do i = 1, ngrdcol
        wpxpyp(i,k) &
        = mixt_frac(i,k) &
          * ( ( w_1(i,k) - wm(i,k) ) * ( x_1(i,k) - xm(i,k) ) * ( y_1(i,k) - ym(i,k) ) &
              + corr_x_y_1(i,k)*sqrt( varnce_x_1(i,k)*varnce_y_1(i,k) )*( w_1(i,k)-wm(i,k) ) &
              + corr_w_y_1(i,k)*sqrt( varnce_w_1(i,k)*varnce_y_1(i,k) )*( x_1(i,k)-xm(i,k) ) &
              + corr_w_x_1(i,k)*sqrt( varnce_w_1(i,k)*varnce_x_1(i,k) )*( y_1(i,k)-ym(i,k) ) ) &
          + ( one - mixt_frac(i,k) ) &
            * ( ( w_2(i,k) - wm(i,k) )*( x_2(i,k) - xm(i,k) ) * ( y_2(i,k) - ym(i,k) ) &
                + corr_x_y_2(i,k)*sqrt( varnce_x_2(i,k)*varnce_y_2(i,k) )*( w_2(i,k)-wm(i,k) ) &
                + corr_w_y_2(i,k)*sqrt( varnce_w_2(i,k)*varnce_y_2(i,k) )*( x_2(i,k)-xm(i,k) ) &
                + corr_w_x_2(i,k)*sqrt( varnce_w_2(i,k)*varnce_x_2(i,k) )*( y_2(i,k)-ym(i,k) ) )
      end do
    end do

    return

  end subroutine calc_wpxpyp_pdf

  !=============================================================================
  subroutine calc_liquid_cloud_frac_component( nz, ngrdcol, &
                                               mean_chi, stdev_chi, &
                                               cloud_frac, rc )
    ! Description:
    ! Calculates the PDF component cloud water mixing ratio, rc_i, and cloud
    ! fraction, cloud_frac_i, for the ith PDF component.
    !
    ! The equation for cloud water mixing ratio, rc, at any point is:
    !
    ! rc = chi * H(chi);
    !
    ! and the equation for cloud fraction at a point, fc, is:
    !
    ! fc = H(chi);
    !
    ! where where extended liquid water mixing ratio, chi, is equal to cloud
    ! water mixing ratio, rc, when positive.  When the atmosphere is saturated
    ! at this point, cloud water is found, and rc = chi, while fc = 1.
    ! Otherwise, clear air is found at this point, and rc = fc = 0.
    !
    ! The mean of rc and fc is calculated by integrating over the PDF, such
    ! that:
    !
    ! <rc> = INT(-inf:inf) chi * H(chi) * P(chi) dchi; and
    !
    ! cloud_frac = <fc> = INT(-inf:inf) H(chi) * P(chi) dchi.
    !
    ! This can be rewritten as:
    !
    ! <rc> = INT(0:inf) chi * P(chi) dchi; and
    !
    ! cloud_frac = <fc> = INT(0:inf) P(chi) dchi;
    !
    ! and further rewritten as:
    !
    ! <rc> = SUM(i=1,N) mixt_frac_i INT(0:inf) chi * P_i(chi) dchi; and
    !
    ! cloud_frac = SUM(i=1,N) mixt_frac_i INT(0:inf) P_i(chi) dchi;
    !
    ! where N is the number of PDF components.  The equation for mean rc in the
    ! ith PDF component is:
    !
    ! rc_i = INT(0:inf) chi * P_i(chi) dchi;
    !
    ! and the equation for cloud fraction in the ith PDF component is:
    ! 
    ! cloud_frac_i = INT(0:inf) P_i(chi) dchi.
    !
    ! The component values are related to the overall values by:
    !
    ! <rc> = SUM(i=1,N) mixt_frac_i * rc_i; and
    !
    ! cloud_frac = SUM(i=1,N) mixt_frac_i * cloud_frac_i.

    ! References:
    !----------------------------------------------------------------------

    use constants_clubb, only: &
        chi_tol,        & ! Tolerance for pdf parameter chi       [kg/kg]
        sqrt_2pi,       & ! sqrt(2*pi)
        sqrt_2,       & ! sqrt(2*pi)
        one,            & ! 1
        one_half,       & ! 1/2
        zero,           & ! 0
        max_num_stdevs, &
        eps

    use clubb_precision, only: &
        core_rknd     ! Precision

    implicit none

    integer, intent(in) :: &
      ngrdcol,  & ! Number of grid columns
      nz          ! Number of vertical level

    !----------- Input Variables -----------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      mean_chi,  & ! Mean of chi (old s) (ith PDF component)           [kg/kg]
      stdev_chi    ! Standard deviation of chi (ith PDF component)     [kg/kg]

    !----------- Output Variables -----------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      cloud_frac, & ! Cloud fraction (ith PDF component)               [-]
      rc            ! Mean cloud water mixing ratio (ith PDF comp.)    [kg/kg]

    !----------- Local Variables -----------
    real( kind = core_rknd), parameter :: &
      invrs_sqrt_2 = one / sqrt_2, &
      invrs_sqrt_2pi = one / sqrt_2pi

    real( kind = core_rknd ) :: &
      zeta

    integer :: k, i    ! Vertical loop index

    !----------- Begin Code -----------

    do k = 1, nz
      do i = 1, ngrdcol

        if ( ( abs( mean_chi(i,k) ) <= eps .and. stdev_chi(i,k) <= chi_tol ) &
               .or. ( mean_chi(i,k) < - max_num_stdevs * stdev_chi(i,k) ) ) then

            ! The mean of chi is at saturation and does not vary in the ith PDF component
            cloud_frac(i,k) = zero
            rc(i,k)         = zero

        elseif ( mean_chi(i,k) > max_num_stdevs * stdev_chi(i,k) ) then

            ! The mean of chi is multiple standard deviations above the saturation point.
            ! Thus, all cloud in the ith PDF component.
            cloud_frac(i,k) = one
            rc(i,k)         = mean_chi(i,k)

        else

            ! The mean of chi is within max_num_stdevs of the saturation point.
            ! Thus, layer is partly cloudy, requires calculation.

            zeta = mean_chi(i,k) / stdev_chi(i,k)

            cloud_frac(i,k) = one_half * ( one + erf( zeta * invrs_sqrt_2 )  )

            rc(i,k) = mean_chi(i,k) * cloud_frac(i,k) &
                      + stdev_chi(i,k) * exp( - one_half * zeta**2 ) * invrs_sqrt_2pi

        end if
        
      end do
    end do

    return

  end subroutine calc_liquid_cloud_frac_component

  !=============================================================================
  elemental subroutine calc_cloud_frac_component( mean_chi, stdev_chi, &
                                                  chi_at_sat, &
                                                  cloud_frac, rc )
    ! Description:
    !   An elemental version of the cloud fraction calculation. This
    !   procedure takes an extra argument, chi_at_sat, allowing it to
    !   be used for the ice_supersat_frac calculation as well. This is
    !   because the saturation point will be non-zero when calculating
    !   ice_supersat_frac if tl is below freezing on that grid level,
    !   unlike the calculation of the liquid cloud fraction, where the
    !   saturation point is always 0. Additionally, the ice_supersat_frac
    !   only needs to be calculated when tl is below freezing, otherwise
    !   it is equal to the liquid cloud fraction component, so being
    !   elemental allows this procedure to be called only for the grid
    !   levels where tl < T_freeze_K.
    !
    !   The description of the equations are located in the description
    !   of calc_liquid_cloud_frac_component.
    !----------------------------------------------------------------------

    use constants_clubb, only: &
        chi_tol,        & ! Tolerance for pdf parameter chi       [kg/kg]
        sqrt_2pi,       & ! sqrt(2*pi)
        sqrt_2,         & ! sqrt(2)
        one,            & ! 1
        one_half,       & ! 1/2
        zero,           & ! 0
        max_num_stdevs, &
        eps

    use clubb_precision, only: &
        core_rknd     ! Precision

    implicit none

    !----------- Input Variables -----------
    real( kind = core_rknd ), intent(in) :: &
      mean_chi,    & ! Mean of chi (old s) (ith PDF component)           [kg/kg]
      stdev_chi,   & ! Standard deviation of chi (ith PDF component)     [kg/kg]
      chi_at_sat     ! Value of chi at saturation (0--liquid; neg.--ice) [kg/kg]

    !----------- Output Variables -----------
    real( kind = core_rknd ), intent(out) :: &
      cloud_frac, & ! Cloud fraction (ith PDF component)               [-]
      rc            ! Mean cloud water mixing ratio (ith PDF comp.)    [kg/kg]

    !----------- Local Variables -----------
    real( kind = core_rknd), parameter :: &
      invrs_sqrt_2 = one / sqrt_2, &
      invrs_sqrt_2pi = one / sqrt_2pi

    real( kind = core_rknd) :: zeta

    !----------- Begin Code -----------

    if ( ( abs( mean_chi - chi_at_sat ) <= eps .and. stdev_chi <= chi_tol ) &
           .or. ( mean_chi - chi_at_sat < - max_num_stdevs * stdev_chi ) ) then

        ! The mean of chi is at saturation and does not vary in the ith PDF component
        cloud_frac = zero
        rc         = zero

    elseif ( mean_chi - chi_at_sat > max_num_stdevs * stdev_chi ) then

        ! The mean of chi is multiple standard deviations above the saturation point.
        ! Thus, all cloud in the ith PDF component.
        cloud_frac = one
        rc         = mean_chi - chi_at_sat
    
    else

        ! The mean of chi is within max_num_stdevs of the saturation point.
        ! Thus, layer is partly cloudy, requires calculation.
        zeta = ( mean_chi - chi_at_sat ) / stdev_chi

        cloud_frac = one_half * ( one + erf( zeta * invrs_sqrt_2 )  )

        rc = ( mean_chi - chi_at_sat ) * cloud_frac &
             + stdev_chi * exp( - one_half * zeta**2 ) * invrs_sqrt_2pi

    end if

    return
    
  end subroutine calc_cloud_frac_component

  !=============================================================================
  subroutine calc_xprcp_component( nz, ngrdcol,                                     & ! In
                                   wm, rtm, thlm, um, vm, rcm,                      & ! In
                                   w_i, rt_i,                                       & ! In
                                   thl_i, u_i, v_i,                                 & ! In
                                   varnce_w_i, chi_i,                               & ! In
                                   stdev_chi_i, stdev_eta_i,                        & ! In
                                   corr_w_chi_i, corr_chi_eta_i,                    & ! In
!                                  corr_u_w_i, corr_v_w_i,                          & ! In
                                   crt_i, cthl_i,                                   & ! In
                                   rc_i, cloud_frac_i, iiPDF_type,                  & ! In
                                   wprcp_contrib_comp_i, wp2rcp_contrib_comp_i,     & ! Out
                                   rtprcp_contrib_comp_i, thlprcp_contrib_comp_i,   & ! Out
                                   uprcp_contrib_comp_i, vprcp_contrib_comp_i )       ! Out

    ! Description:
    ! Calculates the contribution to <w'rc'>, <w'^2 rc'>, <rt'rc'>, and
    ! <thl'rc'> from the ith PDF component.
    !
    !
    ! <w'rc'>
    ! -------
    !
    ! The value of <w'rc'> is calculated by integrating over the PDF:
    !
    ! <w'rc'>
    ! = INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    !   ( w - <w> ) ( rc - <rc> ) P(w,rt,thl) dthl drt dw;
    !
    ! where <w> is the overall mean of w, <rc> is the overall mean of rc, and
    ! P(w,rt,thl) is a two-component trivariate normal distribution of w, rt,
    ! and thl.  This equation is rewritten as:
    !
    ! <w'rc'>
    ! = mixt_frac 
    !   * INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    !     ( w - <w> ) ( rc - <rc> ) P_1(w,rt,thl) dthl drt dw
    !   + ( 1 - mixt_frac )
    !     * INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    !       ( w - <w> ) ( rc - <rc> ) P_2(w,rt,thl) dthl drt dw;
    !
    ! where mixt_frac is the mixture fraction, which is the weight of the 1st
    ! PDF component, and where P_1(w,rt,thl) and P_2(w,rt,thl) are the equations
    ! for the trivariate normal PDF of w, rt, and thl in the 1st and 2nd PDF
    ! components, respectively.  The contribution from the ith PDF component is:
    !
    ! INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    ! ( w - <w> ) ( rc - <rc> ) P_i(w,rt,thl) dthl drt dw;
    !
    ! where P_i(w,rt,thl) is the trivariate normal PDF of w, rt, and thl in the
    ! ith PDF component.  The PDF undergoes a PDF transformation in each PDF
    ! component, which is a change of variables and a translation, stretching,
    ! and rotation of the axes.  The PDF becomes a trivariate normal PDF that is
    ! written in terms of w, chi, and eta coordinates.  Cloud water mixing
    ! ratio, rc, is written in terms of extended liquid water mixing ratio, chi,
    ! such that:
    !
    ! rc = chi H(chi);
    !
    ! where H(chi) is the Heaviside step function.  The contribution from the
    ! ith PDF component to <w'rc'> can be written as:
    !
    ! INT(-inf:inf) INT(-inf:inf)
    ! ( w - <w> ) ( chi H(chi) - <rc> ) P_i(w,chi) dchi dw;
    !
    ! where P_i(w,chi) is the bivariate normal PDF of w and chi in the ith PDF
    ! component.  The solved equation for the <w'rc'> contribution from the ith
    ! PDF component (wprcp_contrib_comp_i) is:
    !
    ! wprcp_contrib_comp_i
    ! = INT(-inf:inf) INT(-inf:inf)
    !   ( w - <w> ) ( chi H(chi) - <rc> ) P_i(w,chi) dchi dw
    ! = ( mu_w_i - <w> )
    !   * ( mu_chi_i * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
    !       + 1/sqrt(2*pi) * sigma_chi_i
    !         * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) } - <rc> )
    !   + corr_w_chi_i * sigma_w_i * sigma_chi_i
    !     * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) );
    !
    ! where mu_w_i is the mean of w in the ith PDF component, mu_chi_i is the
    ! mean of chi in the ith PDF component, sigma_w_i is the standard deviation
    ! of w in the ith PDF component, sigma_chi_i is the standard deviation of
    ! chi in the ith PDF component, and corr_w_chi_i is the correlation of w and
    ! chi in the ith PDF component.
    !
    ! Special case:  sigma_chi_i = 0.
    !
    ! In the special case that sigma_chi_i = 0, chi, as well as rc, are constant
    ! in the ith PDF component.  The equation becomes:
    !
    ! wprcp_contrib_comp_i
    ! = | ( mu_w_i - <w> ) * ( mu_chi_i - <rc> ); when mu_chi_i > 0;
    !   | ( mu_w_i - <w> ) * ( -<rc> ); when mu_chi_i <= 0.
    !
    !
    ! <w'^2 rc'>
    ! ----------
    !
    ! The value of <w'^2 rc'> is calculated by integrating over the PDF:
    !
    ! <w'^2 rc'>
    ! = INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    !   ( w - <w> )^2 ( rc - <rc> ) P(w,rt,thl) dthl drt dw.
    !
    ! This equation is rewritten as:
    !
    ! <w'^2 rc'>
    ! = mixt_frac 
    !   * INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    !     ( w - <w> )^2 ( rc - <rc> ) P_1(w,rt,thl) dthl drt dw
    !   + ( 1 - mixt_frac )
    !     * INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    !       ( w - <w> )^2 ( rc - <rc> ) P_2(w,rt,thl) dthl drt dw.
    !
    ! The contribution from the ith PDF component is:
    !
    ! INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    ! ( w - <w> )^2 ( rc - <rc> ) P_i(w,rt,thl) dthl drt dw.
    !
    ! The PDF undergoes a PDF transformation in each PDF component, and becomes
    ! a trivariate normal PDF that is written in terms of w, chi, and eta
    ! coordinates.  The contribution from the ith PDF component to <w'^2 rc'>
    ! can be written as:
    !
    ! INT(-inf:inf) INT(-inf:inf)
    ! ( w - <w> )^2 ( chi H(chi) - <rc> ) P_i(w,chi) dchi dw.
    !
    ! The solved equation for the <w'^2 rc'> contribution from the ith PDF
    ! component (wp2rcp_contrib_comp_i) is:
    !
    ! wp2rcp_contrib_comp_i
    ! = INT(-inf:inf) INT(-inf:inf)
    !   ( w - <w> )^2 ( chi H(chi) - <rc> ) P_i(w,chi) dchi dw
    ! = ( ( mu_w_i - <w> )^2 + sigma_w_i^2 )
    !   * ( mu_chi_i * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
    !       + 1/sqrt(2*pi) * sigma_chi_i
    !         * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) } - <rc> )
    !   + ( mu_w_i - <w> ) * corr_w_chi_i * sigma_w_i * sigma_chi_i
    !     * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
    !   + 1/sqrt(2*pi) * corr_w_chi_i^2 * sigma_w_i^2 * sigma_chi_i
    !     * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) }.
    !
    ! Special case:  sigma_chi_i = 0.
    !
    ! In the special case that sigma_chi_i = 0, chi, as well as rc, are constant
    ! in the ith PDF component.  The equation becomes:
    !
    ! wp2rcp_contrib_comp_i
    ! = | ( ( mu_w_i - <w> )^2 + sigma_w_i^2 ) * ( mu_chi_i - <rc> );
    !   |     when mu_chi_i > 0;
    !   | ( ( mu_w_i - <w> )^2 + sigma_w_i^2 ) * ( -<rc> );
    !   |     when mu_chi_i <= 0.
    !
    !
    ! <rt'rc'>
    ! --------
    !
    ! The value of <rt'rc'> is calculated by integrating over the PDF:
    !
    ! <rt'rc'>
    ! = INT(-inf:inf) INT(-inf:inf)
    !   ( rt - <rt> ) ( rc - <rc> ) P(rt,thl) dthl drt;
    !
    ! where <rt> is the overall mean of rt, and where P(rt,thl) is a
    ! two-component bivariate normal distribution of rt and thl.  This equation
    ! is rewritten as:
    !
    ! <rt'rc'>
    ! = mixt_frac 
    !   * INT(-inf:inf) INT(-inf:inf)
    !     ( rt - <rt> ) ( rc - <rc> ) P_1(rt,thl) dthl drt
    !   + ( 1 - mixt_frac )
    !     * INT(-inf:inf) INT(-inf:inf)
    !       ( rt - <rt> ) ( rc - <rc> ) P_2(rt,thl) dthl drt;
    !
    ! where P_1(rt,thl) and P_2(rt,thl) are the equations for the bivariate
    ! normal PDF of rt and thl in the 1st and 2nd PDF components, respectively.
    ! The contribution from the ith PDF component is:
    !
    ! INT(-inf:inf) INT(-inf:inf)
    ! ( rt - <rt> ) ( rc - <rc> ) P_i(rt,thl) dthl drt;
    !
    ! where P_i(rt,thl) is the bivariate normal PDF of rt and thl in the ith PDF
    ! component.  The PDF undergoes a PDF transformation in each PDF component,
    ! and becomes a bivariate normal PDF that is written in terms of chi and
    ! eta coordinates.  Total water mixing ratio, rt, is rewritten in terms of
    ! chi and eta by:
    !
    ! rt = mu_rt_i
    !      + ( ( eta - mu_eta_i ) + ( chi - mu_chi_i ) ) / ( 2 * crt_i );
    !
    ! where mu_rt_i is the mean of rt in the ith PDF component, mu_eta_i is the
    ! mean of eta in the ith PDF component, and crt_i is a coefficient on rt in
    ! the chi/eta transformation equations.  The contribution from the ith PDF
    ! component to <rt'rc'> can be written as:
    !
    ! INT(-inf:inf) INT(-inf:inf)
    ! ( mu_rt_i - <rt> + ( eta - mu_eta_i ) / ( 2 * crt_i )
    !   + ( chi - mu_chi_i ) / ( 2 * crt_i ) )
    ! * ( chi H(chi) - <rc> ) P_i(chi,eta) deta dchi;
    !
    ! where P_i(chi,eta) is the bivariate normal PDF of chi and eta in the ith
    ! PDF component.  The solved equation for the <rt'rc'> contribution from the
    ! ith PDF component (rtprcp_contrib_comp_i) is:
    !
    ! rtprcp_contrib_comp_i
    ! = INT(-inf:inf) INT(-inf:inf)
    !   ( mu_rt_i - <rt> + ( eta - mu_eta_i ) / ( 2 * crt_i )
    !     + ( chi - mu_chi_i ) / ( 2 * crt_i ) )
    !   * ( chi H(chi) - <rc> ) P_i(chi,eta) deta dchi
    ! = ( mu_rt_i - <rt> )
    !   * ( mu_chi_i * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
    !       + 1/sqrt(2*pi) * sigma_chi_i
    !         * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) } - <rc> )
    !   + ( corr_chi_eta_i * sigma_eta_i + sigma_chi_i ) / ( 2 * crt_i )
    !     * sigma_chi_i
    !     * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) );
    !
    ! where sigma_eta_i is the standard deviation of eta in the ith PDF
    ! component and corr_chi_eta_i is the correlation of chi and eta in the ith
    ! PDF component.
    !
    ! Special case:  sigma_chi_i = 0.
    !
    ! In the special case that sigma_chi_i = 0, chi, as well as rc, are constant
    ! in the ith PDF component.  The equation becomes:
    !
    ! rtprcp_contrib_comp_i
    ! = | ( mu_rt_i - <rt> ) * ( mu_chi_i - <rc> ); when mu_chi_i > 0;
    !   | ( mu_rt_i - <rt> ) * ( -<rc> ); when mu_chi_i <= 0.
    !
    !
    ! <thl'rc'>
    ! ---------
    !
    ! The value of <thl'rc'> is calculated by integrating over the PDF:
    !
    ! <thl'rc'>
    ! = INT(-inf:inf) INT(-inf:inf)
    !   ( thl - <thl> ) ( rc - <rc> ) P(rt,thl) dthl drt;
    !
    ! where <thl> is the overall mean of thl.  This equation is rewritten as:
    !
    ! <thl'rc'>
    ! = mixt_frac 
    !   * INT(-inf:inf) INT(-inf:inf)
    !     ( thl - <thl> ) ( rc - <rc> ) P_1(rt,thl) dthl drt
    !   + ( 1 - mixt_frac )
    !     * INT(-inf:inf) INT(-inf:inf)
    !       ( thl - <thl> ) ( rc - <rc> ) P_2(rt,thl) dthl drt.
    !
    ! The contribution from the ith PDF component is:
    !
    ! INT(-inf:inf) INT(-inf:inf)
    ! ( thl - <thl> ) ( rc - <rc> ) P_i(rt,thl) dthl drt.
    !
    ! The PDF undergoes a PDF transformation in each PDF component, and becomes
    ! a bivariate normal PDF that is written in terms of chi and eta
    ! coordinates.  Liquid water potential temperature, thl, is rewritten in
    ! terms of chi and eta by:
    !
    ! thl = mu_thl_i
    !       + ( ( eta - mu_eta_i ) - ( chi - mu_chi_i ) ) / ( 2 * cthl_i );
    !
    ! where mu_thl_i is the mean of thl in the ith PDF component and cthl_i is a
    ! coefficient on thl in the chi/eta transformation equations.  The
    ! contribution from the ith PDF component to <thl'rc'> can be written as:
    !
    ! INT(-inf:inf) INT(-inf:inf)
    ! ( mu_thl_i - <thl> + ( eta - mu_eta_i ) / ( 2 * cthl_i )
    !   - ( chi - mu_chi_i ) / ( 2 * cthl_i ) )
    ! * ( chi H(chi) - <rc> ) P_i(chi,eta) deta dchi.
    !
    ! The solved equation for the <thl'rc'> contribution from the ith PDF
    ! component (thlprcp_contrib_comp_i) is:
    !
    ! thlprcp_contrib_comp_i
    ! = INT(-inf:inf) INT(-inf:inf)
    !   ( mu_thl_i - <thl> + ( eta - mu_eta_i ) / ( 2 * cthl_i )
    !     - ( chi - mu_chi_i ) / ( 2 * cthl_i ) )
    !   * ( chi H(chi) - <rc> ) P_i(chi,eta) deta dchi
    ! = ( mu_thl_i - <thl> )
    !   * ( mu_chi_i * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
    !       + 1/sqrt(2*pi) * sigma_chi_i
    !         * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) } - <rc> )
    !   + ( corr_chi_eta_i * sigma_eta_i - sigma_chi_i ) / ( 2 * cthl_i )
    !     * sigma_chi_i
    !     * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) ).
    !
    ! Special case:  sigma_chi_i = 0.
    !
    ! In the special case that sigma_chi_i = 0, chi, as well as rc, are constant
    ! in the ith PDF component.  The equation becomes:
    !
    ! thlprcp_contrib_comp_i
    ! = | ( mu_thl_i - <thl> ) * ( mu_chi_i - <rc> ); when mu_chi_i > 0;
    !   | ( mu_thl_i - <thl> ) * ( -<rc> ); when mu_chi_i <= 0.
    !
    !
    ! Use equations for PDF component cloud fraction cloud water mixing ratio
    ! -----------------------------------------------------------------------
    !
    ! The equation for cloud fraction in the ith PDF component, fc_i, is:
    !
    ! fc_i = 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) ).
    !
    ! In the special case that sigma_chi_i = 0, the equation becomes:
    !
    ! fc_i = | 1; when mu_chi_i > 0;
    !        | 0; when mu_chi_i <= 0.
    !
    ! The equation for mean cloud water mixing ratio in the ith PDF component,
    ! rc_i, is:
    !
    ! rc_i
    ! = mu_chi_i * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
    !   + 1/sqrt(2*pi) * sigma_chi_i
    !     * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) }
    ! = mu_chi_i * fc_i
    !   + 1/sqrt(2*pi) * sigma_chi_i
    !     * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) }.
    !
    ! In the special case that sigma_chi_i = 0, the equation becomes:
    !
    ! rc_i = | mu_chi_i; when mu_chi_i > 0;
    !        | 0; when mu_chi_i <= 0.
    !
    ! The above equations can be substituted into the equations for
    ! wprcp_contrib_comp_i, wp2rcp_contrib_comp_i, rtprcp_contrib_comp_i, and
    ! thlprcp_contrib_comp_i.  The new equations are:
    !
    ! wprcp_contrib_comp_i
    ! = ( mu_w_i - <w> ) * ( rc_i - <rc> )
    !   + corr_w_chi_i * sigma_w_i * sigma_chi_i * fc_i;
    !
    ! wp2rcp_contrib_comp_i
    ! = ( ( mu_w_i - <w> )^2 + sigma_w_i^2 ) * ( rc_i - <rc> )
    !   + 2 * ( mu_w_i - <w> ) * corr_w_chi_i * sigma_w_i * sigma_chi_i * fc_i
    !   + 1/sqrt(2*pi) * corr_w_chi_i^2 * sigma_w_i^2 * sigma_chi_i
    !     * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) };
    !
    ! rtprcp_contrib_comp_i
    ! = ( mu_rt_i - <rt> ) * ( rc_i - <rc> )
    !   + ( corr_chi_eta_i * sigma_eta_i + sigma_chi_i ) / ( 2 * crt_i )
    !     * sigma_chi_i * fc_i; and
    !
    ! thlprcp_contrib_comp_i
    ! = ( mu_thl_i - <thl> ) * ( rc_i - <rc> )
    !   + ( corr_chi_eta_i * sigma_eta_i - sigma_chi_i ) / ( 2 * cthl_i )
    !     * sigma_chi_i * fc_i.
    !
    ! While the above equations reduce to their listed versions in the special
    ! case that sigma_chi_i = 0, those versions are faster to calculate.  When
    ! mu_chi_i > 0, they are:
    !
    ! wprcp_contrib_comp_i = ( mu_w_i - <w> ) * ( mu_chi_i - <rc> );
    ! wp2rcp_contrib_comp_i
    ! = ( ( mu_w_i - <w> )^2 + sigma_w_i^2 ) * ( mu_chi_i - <rc> );
    ! rtprcp_contrib_comp_i = ( mu_rt_i - <rt> ) * ( mu_chi_i - <rc> ); and
    ! thlprcp_contrib_comp_i = ( mu_thl_i - <thl> ) * ( mu_chi_i - <rc> );
    !
    ! and when mu_chi_i <= 0, they are:
    !
    ! wprcp_contrib_comp_i = - ( mu_w_i - <w> ) * <rc>;
    ! wp2rcp_contrib_comp_i = - ( ( mu_w_i - <w> )^2 + sigma_w_i^2 ) * <rc>;
    ! rtprcp_contrib_comp_i = - ( mu_rt_i - <rt> ) * <rc>; and
    ! thlprcp_contrib_comp_i = - ( mu_thl_i - <thl> ) * <rc>.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid ! Type

    use constants_clubb, only: &
        sqrt_2pi,       & ! Variable(s)
        two,            &
        zero,           &
        chi_tol

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    integer, intent(in) :: &
      ngrdcol,  & ! Number of grid columns
      nz          ! Number of vertical level

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      wm,             & ! Mean of w (overall)                          [m/s]
      rtm,            & ! Mean of rt (overall)                         [kg/kg]
      thlm,           & ! Mean of thl (overall)                        [K]
      um,             & ! Mean of eastward wind (overall)              [m/s]
      vm,             & ! Mean of northward wind (overall)             [m/s]
      rcm,            & ! Mean of rc (overall)                         [kg/kg]
      w_i,            & ! Mean of w (ith PDF component)                [m/s]
      rt_i,           & ! Mean of rt (ith PDF component)               [kg/kg]
      thl_i,          & ! Mean of thl (ith PDF component)              [K]
      u_i,            & ! Mean of eastward wind (ith PDF component)    [m/s]
      v_i,            & ! Mean of northward wind (ith PDF component)   [m/s]
      varnce_w_i,     & ! Variance of w (ith PDF component)            [m^2/s^2]
      chi_i,          & ! Mean of chi (ith PDF component)              [kg/kg]
      stdev_chi_i,    & ! Standard deviation of chi (ith PDF comp.)    [kg/kg]
      stdev_eta_i,    & ! Standard deviation of eta (ith PDF comp.)    [kg/kg]
      corr_w_chi_i,   & ! Correlation of w and chi (ith PDF component) [-]
      corr_chi_eta_i, & ! Correlation of chi and eta (ith PDF comp.)   [-]
!     corr_u_w_i,     & ! Correlation of u and w (ith PDF component)   [-]
!     corr_v_w_i,     & ! Correlation of v and w (ith PDF component)   [-]
      crt_i,          & ! Coef. on rt in chi/eta eqns. (ith PDF comp.) [-]
      cthl_i,         & ! Coef. on thl: chi/eta eqns. (ith PDF comp.)  [kg/kg/K]
      rc_i,           & ! Mean of rc (ith PDF component)               [kg/kg]
      cloud_frac_i      ! Cloud fraction (ith PDF component)           [-]

    integer, intent(in) :: &
      iiPDF_type    ! Selected option for the two-component normal (double
                    ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                    ! w, chi, and eta) portion of CLUBB's multivariate,
                    ! two-component PDF.

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      wprcp_contrib_comp_i,   & ! <w'rc'> contrib. (ith PDF comp.)  [m/s(kg/kg)]
      wp2rcp_contrib_comp_i,  & ! <w'^2rc'> contrib. (ith comp) [m^2/s^2(kg/kg)]
      rtprcp_contrib_comp_i,  & ! <rt'rc'> contrib. (ith PDF comp.)  [kg^2/kg^2]
      thlprcp_contrib_comp_i, & ! <thl'rc'> contrib. (ith PDF comp.)  [K(kg/kg)]
      uprcp_contrib_comp_i,   & ! <u'rc'> contrib. (ith PDF comp.)  [m/s(kg/kg)]
      vprcp_contrib_comp_i      ! <v'rc'> contrib. (ith PDF comp.)  [m/s(kg/kg)]
      
    ! Local Variables
    integer :: i, k

    ! ---------------------- Begin Code ------------------
    
    ! Changing these conditionals may result in inconsistencies with the conditional
    ! statements located in calc_cloud_frac_component
    
    do k = 1, nz
      do i = 1, ngrdcol

        wprcp_contrib_comp_i(i,k) = ( w_i(i,k) - wm(i,k) ) * ( rc_i(i,k) - rcm(i,k) )

        wp2rcp_contrib_comp_i(i,k) = ( ( w_i(i,k) - wm(i,k) )**2 + varnce_w_i(i,k) ) &
                                     * ( rc_i(i,k) - rcm(i,k) )

        rtprcp_contrib_comp_i(i,k) = ( rt_i(i,k) - rtm(i,k) ) * ( rc_i(i,k) - rcm(i,k) ) &
                                + ( corr_chi_eta_i(i,k) * stdev_eta_i(i,k) + stdev_chi_i(i,k) ) &
                                  / ( two * crt_i(i,k) ) * stdev_chi_i(i,k) * cloud_frac_i(i,k)

        thlprcp_contrib_comp_i(i,k) = ( thl_i(i,k) - thlm(i,k) ) * ( rc_i(i,k) - rcm(i,k) ) &
                                 + ( corr_chi_eta_i(i,k) * stdev_eta_i(i,k) - stdev_chi_i(i,k) ) &
                                   / ( two * cthl_i(i,k) ) * stdev_chi_i(i,k) * cloud_frac_i(i,k)

        uprcp_contrib_comp_i(i,k) = ( u_i(i,k) - um(i,k) ) * ( rc_i(i,k) - rcm(i,k) )

        vprcp_contrib_comp_i(i,k) = ( v_i(i,k) - vm(i,k) ) * ( rc_i(i,k) - rcm(i,k) )
        
      end do
    end do

    ! If iiPDF_type isn't iiPDF_ADG1, iiPDF_ADG2, or iiPDF_new_hybrid, so
    ! corr_w_chi_i /= 0 (and perhaps corr_u_w_i /= 0).
    if ( .not. ( iiPDF_type == iiPDF_ADG1 .or. iiPDF_type == iiPDF_ADG2 &
                 .or. iiPDF_type == iiPDF_new_hybrid ) ) then

        ! Chi varies significantly in the ith PDF component (stdev_chi > chi_tol)
        ! and there is some cloud (0 < cloud_frac <= 1)
        do k = 1, nz
          do i = 1, ngrdcol
            if ( stdev_chi_i(i,k) > chi_tol .and. cloud_frac_i(i,k) > zero ) then

              wprcp_contrib_comp_i(i,k) = wprcp_contrib_comp_i(i,k) &
                                          + corr_w_chi_i(i,k) * sqrt( varnce_w_i(i,k) ) &
                                            * stdev_chi_i(i,k) * cloud_frac_i(i,k)

              wp2rcp_contrib_comp_i(i,k) = wp2rcp_contrib_comp_i(i,k) &
                                           + two * ( w_i(i,k) - wm(i,k) ) * corr_w_chi_i(i,k) &
                                             * sqrt( varnce_w_i(i,k) ) * stdev_chi_i(i,k) &
                                             * cloud_frac_i(i,k) &
                                           + corr_w_chi_i(i,k)**2 * varnce_w_i(i,k) &
                                             * stdev_chi_i(i,k) &
                                             * exp( -chi_i(i,k)**2 / ( two*stdev_chi_i(i,k)**2 ) ) &
                                               / sqrt_2pi

            ! In principle, uprcp_contrib_comp_i might depend on corr_u_w_i here.
          end if
        end do
      end do

    end if 

    return

  end subroutine calc_xprcp_component

  !=============================================================================
  subroutine calc_vert_avg_cf_component &
                  ( nz, k, z_vals, chi, stdev_chi, chi_at_sat, &
                    cloud_frac_i, rc_i )
  ! Description:
  !   This subroutine is similar to calc_cloud_frac_component, but
  !   resolves cloud_frac and rc at an arbitrary number of vertical levels
  !   in the vicinity of the desired level. This may give a better
  !   parameterization of sub-grid atmospheric conditions.
  !
  ! References:
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    intrinsic :: sum

    ! Local Constants
    integer, parameter :: &
      n_points = 9       ! Number of vertical levels to use in averaging
                         ! (arbitrary, but must be odd)

    ! Input Variables
    integer, intent(in) :: &
      nz, &       ! Number of vertical levels                         [count]
      k           ! Level at which cloud_frac is to be computed       [count]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      z_vals,    & ! Height at each vertical level                   [m]
      chi,       & ! Value of chi (old s)                            [kg/kg]
      stdev_chi, & ! Standard deviation of chi                       [kg/kg]
      chi_at_sat   ! Value of chi at saturation with respect to ice  [kg/kg]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      cloud_frac_i, & ! Vertically averaged cloud fraction               [-]
      rc_i            ! Vertically averaged cloud water mixing ratio     [kg/kg]

    ! Local Variables
    real( kind = core_rknd ), dimension(n_points) :: &
      chi_ref,           &   ! chi (old s) evaluated on refined grid     [kg/kg]
      stdev_chi_ref,     &   ! stdev_chi evaluated on refined grid       [kg/kg]
      cloud_frac_ref,    &   ! cloud_frac evaluated on refined grid      [-]
      rc_ref                 ! r_c evaluated on refined grid             [kg/kg]
      
  !-----------------------------------------------------------------------

    !----- Begin Code -----
    chi_ref = interp_var_array( n_points, nz, k, z_vals, chi )
    stdev_chi_ref = interp_var_array( n_points, nz, k, z_vals, stdev_chi )
    ! We could optionally compute chi_at_sat in an analogous manner. For now,
    ! use chi_at_sat(k) as an approximation.

    ! Compute cloud_frac and r_c at each refined grid level
    call calc_cloud_frac_component( chi_ref(:), stdev_chi_ref(:), chi_at_sat(k), & ! Intent(in)
                                    cloud_frac_ref(:), rc_ref(:) )                  ! Intent(out)

    cloud_frac_i = sum( cloud_frac_ref(:) ) / real( n_points, kind=core_rknd )
    rc_i = sum( rc_ref(:) ) / real( n_points, kind=core_rknd )

    return
  end subroutine calc_vert_avg_cf_component

  !=============================================================================
  subroutine calc_w_up_in_cloud( nz, ngrdcol, &
                                 mixt_frac, cloud_frac_1, cloud_frac_2, &
                                 w_1, w_2, varnce_w_1, varnce_w_2, w_up_in_cloud)
    ! Description:
    ! Subroutine that computes the mean cloudy updraft.
    !
    ! In order to activate aerosol, we'd like to feed the activation scheme
    ! a vertical velocity that's representative of cloudy updrafts. For skewed
    ! layers, like cumulus layers, this might be an improvement over the square
    ! root of wp2 that's currently used. At the same time, it would be simpler
    ! and less expensive than feeding SILHS samples into the aerosol code
    ! (see larson-group/e3sm#19 and larson-group/e3sm#26).
    !
    ! The formulas are only valid for certain PDFs in CLUBB (ADG1, ADG2,
    ! new hybrid), hence we omit calculation if another PDF type is used.
    !
    ! References: https://www.overleaf.com/project/614a136d47846639af22ae34
    !----------------------------------------------------------------------

    use constants_clubb, only: &
        sqrt_2pi, & ! sqrt(2*pi)
        sqrt_2,   & ! sqrt(2)
        one,      & ! 1
        one_half, & ! 1/2
        eps

    use clubb_precision, only: &
        core_rknd     ! Precision

    implicit none

    integer, intent(in) :: &
      ngrdcol,  & ! Number of grid columns
      nz          ! Number of vertical level

    !----------- Input Variables -----------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      mixt_frac, &      ! mixture fraction                             [-]
      cloud_frac_1, &   ! cloud fraction (1st PDF component)           [-]
      cloud_frac_2, &   ! cloud fraction (2nd PDF component)           [-]
      w_1, &            ! upward velocity (1st PDF component)          [m/s]
      w_2, &            ! upward velocity (2nd PDF component)          [m/s]
      varnce_w_1, &     ! standard deviation of w (1st PDF component)  [m^2/s^2]
      varnce_w_2        ! standard deviation of w (2nd PDF component)  [m^2/s^2]

    !----------- Output Variables -----------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      w_up_in_cloud ! mean updraft over clouds                         [m/s]

    !----------- Local Variables -----------
    real( kind = core_rknd ) :: &
      w_up_1, w_up_2, &       ! product of w and Heaviside function
      stdev_w_1, stdev_w_2  ! Standard deviation of w
      
    integer :: i, k
      
    do k = 1, nz
      do i = 1, ngrdcol
        
        stdev_w_1 = sqrt(varnce_w_1(i,k))
        stdev_w_2 = sqrt(varnce_w_2(i,k))
        
        w_up_1 &
        = one_half * w_1(i,k) &
            * (one + erf(w_1(i,k) / (sqrt_2 * max(eps, stdev_w_1)))) &
          + stdev_w_1 / sqrt_2pi &
              * exp(-one_half * (w_1(i,k) / max(eps, stdev_w_1)) ** 2)
        w_up_2 &
        = one_half * w_2(i,k) &
            * (one + erf(w_2(i,k) / (sqrt_2 * max(eps, stdev_w_2)))) &
          + stdev_w_2 / sqrt_2pi &
              * exp(-one_half * (w_2(i,k) / max(eps, stdev_w_2)) ** 2)

        w_up_in_cloud &
        = (mixt_frac(i,k) * cloud_frac_1(i,k) * w_up_1 &
            + (one - mixt_frac(i,k)) * cloud_frac_2(i,k) * w_up_2) &
          / (mixt_frac(i,k) * max(eps, cloud_frac_1(i,k)) &
            + (one - mixt_frac(i,k)) * max(eps, cloud_frac_2(i,k)))
            
      end do
    end do

  end subroutine calc_w_up_in_cloud

  !=============================================================================
  function interp_var_array( n_points, nz, k, z_vals, var )

  ! Description:
  !   Interpolates a variable to an array of values about a given level

  ! References
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd           ! Constant

    implicit none

  ! Input Variables
    integer, intent(in) :: &
      n_points, & ! Number of points to interpolate to (must be odd and >= 3)
      nz,       & ! Total number of vertical levels
      k           ! Center of interpolation array

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      z_vals, &         ! Height at each vertical level           [m]
      var               ! Variable values on grid                 [units vary]

  ! Output Variables
    real( kind = core_rknd ), dimension(n_points) :: &
      interp_var_array  ! Interpolated values of variable         [units vary]

  ! Local Variables
    real( kind = core_rknd ) :: &
      dz    ! Distance between vertical levels

    real( kind = core_rknd ) :: &
      z_val             ! Height at some sub-grid level

    integer :: &
      i, &                      ! Loop iterator

      subgrid_lev_count         ! Number of refined grid points located between
                              ! two defined grid levels

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    ! Place a point at each of k-1, k, and k+1.
    interp_var_array(1) = var_value_integer_height( nz, k-1, z_vals, var )
    interp_var_array((n_points+1)/2) = var_value_integer_height( nz, k, z_vals, var )
    interp_var_array(n_points) = var_value_integer_height( nz, k+1, z_vals, var )

    subgrid_lev_count = (n_points - 3) / 2

    ! Lower half
    if ( k == 1 ) then
      dz = (z_vals(2) - z_vals(1)) / real( subgrid_lev_count+1, kind=core_rknd )
    else
      dz = (z_vals(k) - z_vals(k-1)) / real( subgrid_lev_count+1, kind=core_rknd )
    end if
    do i=1, subgrid_lev_count
      z_val = z_vals(k) - real( i, kind=core_rknd ) * dz
      interp_var_array(1+i) &
      = var_subgrid_interp( nz, k, z_vals, var, z_val, l_below=.true. )
    end do

    ! Upper half
    if ( k == nz ) then
      dz = ( z_vals(nz) - z_vals(nz-1) ) / real( subgrid_lev_count+1, kind=core_rknd )
    else
      dz = ( z_vals(k+1) - z_vals(k) ) / real( subgrid_lev_count+1, kind=core_rknd )
    end if
    do i=1, (n_points-3)/2
      z_val = z_vals(k) + real( i, kind=core_rknd ) * dz
      interp_var_array((n_points+1)/2+i) &
      = var_subgrid_interp( nz, k, z_vals, var, z_val, l_below=.false. )
    end do

    return
  end function interp_var_array

  !=============================================================================
  function var_value_integer_height( nz, k, z_vals, var_grid_value ) result( var_value )

  ! Description
  !   Returns the value of a variable at an integer height between 0 and
  !   nz+1 inclusive, using extrapolation when k==0 or k==nz+1

  ! References
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd       ! Constant

    use interpolation, only: &
        mono_cubic_interp  ! Procedure

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,    & ! Total number of vertical levels
      k        ! Level to resolve variable value

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      z_vals,            & ! Height at each vertical level                  [m]
      var_grid_value       ! Value of variable at each grid level           [units vary]

    ! Output Variables
    real( kind = core_rknd ) :: &
      var_value            ! Value of variable at height level              [units vary]

    ! Local Variables
    integer :: km1, k00, kp1, kp2
  !-----------------------------------------------------------------------

    !----- Begin Code -----

    if ( k >= 1 .and. k <= nz ) then
      ! This is the simple case. No extrapolation necessary.
      var_value = var_grid_value(k)
    else if ( k == 0 ) then
      ! Extrapolate below the lower boundary
      km1 = nz
      k00 = 1
      kp1 = 2
      kp2 = 3
      var_value = mono_cubic_interp( z_vals(1)-(z_vals(2)-z_vals(1)), &
                                     km1, k00, kp1, kp2, &
                                     z_vals(km1), z_vals(k00), z_vals(kp1), z_vals(kp2), &
                                     var_grid_value(km1), var_grid_value(k00), &
                                     var_grid_value(kp1), var_grid_value(kp2) )
    else if ( k == nz+1 ) then
      ! Extrapolate above the upper boundary
      km1 = nz
      k00 = nz-1
      kp1 = nz
      kp2 = nz
      var_value = mono_cubic_interp( z_vals(nz)+(z_vals(nz)-z_vals(nz-1)), &
                                     km1, k00, kp1, kp2, &
                                     z_vals(km1), z_vals(k00), z_vals(kp1), z_vals(kp2), &
                                     var_grid_value(km1), var_grid_value(k00), &
                                     var_grid_value(kp1), var_grid_value(kp2) )
    else
      ! Invalid height requested
      var_value = -999._core_rknd
    end if ! k > 1 .and. k < nz
    return
  end function var_value_integer_height

  !=============================================================================
  function var_subgrid_interp( nz, k, z_vals, var, z_interp, l_below ) result( var_value )

  ! Description
  !   Interpolates (or extrapolates) a variable to a value between grid
  !   levels

  ! References
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd       ! Constant

    use interpolation, only: &
        mono_cubic_interp   ! Procedure

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz, &         ! Number of vertical levels
      k             ! Grid level near interpolation target

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      z_vals, &     ! Height at each grid level          [m]
      var           ! Variable values at grid levels     [units vary]

    real( kind = core_rknd ), intent(in) :: &
      z_interp      ! Interpolation target height        [m]

    logical, intent(in) :: &
      l_below       ! True if z_interp < z_vals(k), false otherwise

    ! Output Variable
    real( kind = core_rknd ) :: &
      var_value     ! Interpolated value of variable     [units vary]

    ! Local Variables
    integer :: km1, k00, kp1, kp2 ! Parameters for call to mono_cubic_interp
  !----------------------------------------------------------------------

    !----- Begin Code -----
    if ( l_below ) then

      if ( k == 1 ) then ! Extrapolation
        km1 = nz
        k00 = 1
        kp1 = 2
        kp2 = 3
      else if ( k == 2 ) then
        km1 = 1
        k00 = 1
        kp1 = 2
        kp2 = 3
      else if ( k == nz ) then
        km1 = nz-2
        k00 = nz-1
        kp1 = nz
        kp2 = nz
      else
        km1 = k-2
        k00 = k-1
        kp1 = k
        kp2 = k+1
      end if ! k == 1

    else ! .not. l_below

      if ( k == 1 ) then
        km1 = 1
        k00 = 1
        kp1 = 2
        kp2 = 3
      else if ( k == nz-1 ) then
        km1 = nz-2
        k00 = nz-1
        kp1 = nz
        kp2 = nz
      else if ( k == nz ) then ! Extrapolation
        km1 = nz
        k00 = nz-1
        kp1 = nz
        kp2 = nz
      else
        km1 = k-1
        k00 = k
        kp1 = k+1
        kp2 = k+2
      end if ! k == 1

    end if ! l_below

    ! Now perform the interpolation
    var_value = mono_cubic_interp( z_interp, km1, k00, kp1, kp2, &
                                   z_vals(km1), z_vals(k00), z_vals(kp1), z_vals(kp2), &
                                   var(km1), var(k00), var(kp1), var(kp2) )

    return

  end function var_subgrid_interp

  !=============================================================================

end module pdf_closure_module
