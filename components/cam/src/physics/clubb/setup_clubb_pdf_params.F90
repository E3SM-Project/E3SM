!-------------------------------------------------------------------------
! $Id: setup_clubb_pdf_params.F90 8196 2016-07-13 18:18:23Z raut@uwm.edu $
!===============================================================================
module setup_clubb_pdf_params

  implicit none

  private

  public :: setup_pdf_parameters,      &
            compute_mean_stdev,        &
            calc_comp_mu_sigma_hm,     &
            norm_transform_mean_stdev, &
            comp_corr_norm,            &
            denorm_transform_corr

  private :: calc_mu_sigma_two_comps,     &
             component_corr_w_x,          &
             component_corr_chi_eta,      &
             component_corr_w_hm_n_ip,    &
             component_corr_x_hm_n_ip,    &
             component_corr_hmx_hmy_n_ip, &
             component_corr_eta_hm_n_ip,  &
             calc_corr_w_hm_n,            &
             pdf_param_hm_stats,          &
             pdf_param_ln_hm_stats,       &
             pack_pdf_params,             &
             compute_rtp2_from_chi

  ! Prescribed parameters are set to in-cloud or outside-cloud (below-cloud)
  ! values based on whether or not cloud water mixing ratio has a value of at
  ! least rc_tol.  However, this does not take into account the amount of
  ! cloudiness in a component, just whether or not there is any cloud in the
  ! component.  The option l_interp_prescribed_params allows for an interpolated
  ! value between the in-cloud and below-cloud parameter value based on the
  ! component cloud fraction.
  logical, parameter, private :: &
    l_interp_prescribed_params = .false.

  contains

  !=============================================================================
  subroutine setup_pdf_parameters( nz, d_variables, dt, &                      ! Intent(in)
                                   Nc_in_cloud, rcm, cloud_frac, &             ! Intent(in)
                                   ice_supersat_frac, hydromet, wphydrometp, & ! Intent(in)
                                   corr_array_n_cloud, corr_array_n_below, &   ! Intent(in)
                                   pdf_params, l_stats_samp, &                 ! Intent(in)
                                   hydrometp2, &                               ! Intent(out)
                                   mu_x_1_n, mu_x_2_n, &                       ! Intent(out)
                                   sigma_x_1_n, sigma_x_2_n, &                 ! Intent(out)
                                   corr_array_1_n, corr_array_2_n, &           ! Intent(out)
                                   corr_cholesky_mtx_1, corr_cholesky_mtx_2, & ! Intent(out)
                                   hydromet_pdf_params )                       ! Intent(out)

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr,    & ! Variable(s)
        zm2zt, & ! Procedure(s)
        zt2zm

    use constants_clubb, only: &
        one,            & ! Constant(s)
        zero,           &
        rc_tol,         &
        Ncn_tol,        &
        cloud_frac_min, &
        fstderr,        &
        zero_threshold

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s)

    use hydromet_pdf_parameter_module, only: &
        hydromet_pdf_parameter,   &  ! Type
        init_hydromet_pdf_params     ! Procedure

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use model_flags, only: &
        l_use_precip_frac,   & ! Flag(s)
        l_calc_w_corr

    use array_index, only: &
        hydromet_list, &       ! Variable(s)
        hydromet_tol

    use model_flags, only: &
        l_const_Nc_in_cloud    ! Flag(s)

    use precipitation_fraction, only: &
        precip_fraction

    use Nc_Ncn_eqns, only: &
        Nc_in_cloud_to_Ncnm  ! Procedure(s)

    use advance_windm_edsclrm_module, only: &
        xpwp_fnc

    use variables_diagnostic_module, only: &
        Kh_zm

    use parameters_tunable, only: &
        c_K_hm

    use pdf_utilities, only: &
        calc_xp2,                  &  ! Procedure(s)
        compute_mean_binormal,     &
        compute_variance_binormal, &
        stdev_L2N,                 &
        corr_NN2NL

    use clip_explicit, only: &
        clip_covar_level, & ! Procedure(s)
        clip_wphydrometp    ! Variables(s)

    use clubb_precision, only: &
        core_rknd,      & ! Variable(s)
        dp

    use matrix_operations, only: &
        Cholesky_factor, & ! Procedure(s)
        mirror_lower_triangular_matrix

    use stats_type_utilities, only: &
        stat_update_var,    & ! Procedure(s)
        stat_update_var_pt

    use stats_variables, only: &
        iprecip_frac,   & ! Variable(s)
        iprecip_frac_1, &
        iprecip_frac_2, &
        iNcnm,          &
        ihmp2_zt,       &
        irtp2_from_chi, &
        stats_zt,             &
        stats_zm

    use model_flags, only: &
        l_diagnose_correlations ! Variable(s)

    use diagnose_correlations_module, only: &
        diagnose_correlations, & ! Procedure(s)
        calc_cholesky_corr_mtx_approx

    use corr_varnce_module, only: &
        assert_corr_symmetric, & ! Procedure(s)
        iiPDF_Ncn,             & ! Variable(s)
        iiPDF_chi,             &
        iiPDF_eta,             &
        hmp2_ip_on_hmm2_ip,    &
        Ncnp2_on_Ncnm2

    use error_code, only : &
        clubb_at_least_debug_level   ! Procedure(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,          & ! Number of model vertical grid levels
      d_variables    ! Number of variables in the correlation array

    real( kind = core_rknd ), intent(in) ::  &
      dt    ! Model timestep                                           [s]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      Nc_in_cloud,       & ! Mean (in-cloud) cloud droplet conc.       [num/kg]
      rcm,               & ! Mean cloud water mixing ratio, < r_c >    [kg/kg]
      cloud_frac,        & ! Cloud fraction                            [-]
      ice_supersat_frac    ! Ice supersaturation fraction              [-]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet,    & ! Mean of hydrometeor, hm (overall) (t-levs.) [units]
      wphydrometp    ! Covariance < w'h_m' > (momentum levels)     [(m/s)units]

    real( kind = core_rknd ), dimension(d_variables,d_variables), &
    intent(in) :: &
      corr_array_n_cloud, & ! Prescribed normal space corr. array in cloud  [-]
      corr_array_n_below    ! Prescribed normal space corr. array below cl. [-]

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params    ! PDF parameters                               [units vary]

    logical, intent(in) :: &
      l_stats_samp    ! Flag to sample statistics

    ! Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      hydrometp2    ! Variance of a hydrometeor (overall) (m-levs.)   [units^2]

    ! Output Variables
    real( kind = core_rknd ), dimension(d_variables, nz), intent(out) :: &
      mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(d_variables,d_variables,nz), &
    intent(out) :: &
      corr_array_1_n, & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
      corr_array_2_n    ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]

    real( kind = core_rknd ), dimension(d_variables,d_variables,nz), &
    intent(out) :: &
      corr_cholesky_mtx_1, & ! Transposed corr. cholesky matrix, 1st comp. [-]
      corr_cholesky_mtx_2    ! Transposed corr. cholesky matrix, 2nd comp. [-]

    type(hydromet_pdf_parameter), dimension(nz), intent(out) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]

    ! Local Variables

    real( kind = core_rknd ), dimension(d_variables,d_variables) :: &
      corr_mtx_approx_1, & ! Approximated corr. matrix (C = LL'), 1st comp. [-]
      corr_mtx_approx_2    ! Approximated corr. matrix (C = LL'), 2nd comp. [-]

    real( kind = core_rknd ), dimension(nz) :: &
      mu_w_1,       & ! Mean of w (1st PDF component)                    [m/s]
      mu_w_2,       & ! Mean of w (2nd PDF component)                    [m/s]
      mu_chi_1,     & ! Mean of chi (old s) (1st PDF component)          [kg/kg]
      mu_chi_2,     & ! Mean of chi (old s) (2nd PDF component)          [kg/kg]
      sigma_w_1,    & ! Standard deviation of w (1st PDF component)      [m/s]
      sigma_w_2,    & ! Standard deviation of w (2nd PDF component)      [m/s]
      sigma_chi_1,  & ! Standard deviation of chi (1st PDF component)    [kg/kg]
      sigma_chi_2,  & ! Standard deviation of chi (2nd PDF component)    [kg/kg]
      rc_1,         & ! Mean of r_c (1st PDF component)                  [kg/kg]
      rc_2,         & ! Mean of r_c (2nd PDF component)                  [kg/kg]
      cloud_frac_1, & ! Cloud fraction (1st PDF component)               [-]
      cloud_frac_2, & ! Cloud fraction (2nd PDF component)               [-]
      mixt_frac       ! Mixture fraction                                 [-]

    real( kind = core_rknd ), dimension(nz) :: &
      ice_supersat_frac_1, & ! Ice supersaturation fraction (1st PDF comp.)  [-]
      ice_supersat_frac_2    ! Ice supersaturation fraction (2nd PDF comp.)  [-]

    real( kind = core_rknd ), dimension(nz) :: &
      Ncnm    ! Mean cloud nuclei concentration, < N_cn >        [num/kg]

    real( kind = core_rknd ), dimension(nz) ::  &
      wpchip_zm, & ! Covariance of chi and w (momentum levels)   [(m/s)(kg/kg)]
      wpNcnp_zm, & ! Covariance of N_cn and w (momentum levs.)   [(m/s)(num/kg)]
      wpchip_zt, & ! Covariance of chi and w on t-levs           [(m/s)(kg/kg)]
      wpNcnp_zt    ! Covariance of N_cn and w on t-levs          [(m/s)(num/kg)]

    real( kind = core_rknd ), dimension(nz,hydromet_dim) :: &
      hm_1, & ! Mean of a precip. hydrometeor (1st PDF component)   [units vary]
      hm_2    ! Mean of a precip. hydrometeor (2nd PDF component)   [units vary]

    real( kind = core_rknd ), dimension(nz,hydromet_dim) :: &
      hydrometp2_zt,  & ! Variance of a hydrometeor (overall); t-lev   [units^2]
      wphydrometp_zt    ! Covariance of w and hm interp. to t-levs. [(m/s)units]

    real( kind = core_rknd ), dimension(nz) :: &
      precip_frac,   & ! Precipitation fraction (overall)           [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component) [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component) [-]

    real( kind = core_rknd ), dimension(d_variables,d_variables) :: &
      corr_array_1, & ! Correlation array of PDF vars. (comp. 1)             [-]
      corr_array_2    ! Correlation array of PDF vars. (comp. 2)             [-]

    real( kind = core_rknd ), dimension(d_variables) :: &
      mu_x_1,    & ! Mean array of PDF vars. (1st PDF component)    [units vary]
      mu_x_2,    & ! Mean array of PDF vars. (2nd PDF component)    [units vary]
      sigma_x_1, & ! Standard deviation array of PDF vars (comp. 1) [units vary]
      sigma_x_2    ! Standard deviation array of PDF vars (comp. 2) [units vary]

    real( kind = core_rknd ), dimension(d_variables) :: &
      corr_array_scaling

    real( kind = core_rknd ), dimension(d_variables) :: &
      sigma2_on_mu2_ip_1, & ! Ratio array sigma_hm_1^2/mu_hm_1^2             [-]
      sigma2_on_mu2_ip_2    ! Ratio array sigma_hm_2^2/mu_hm_2^2             [-]

    real( kind = core_rknd ) :: &
      const_Ncnp2_on_Ncnm2, & ! Prescribed ratio of <Ncn'^2> to <Ncn>^2      [-]
      const_corr_chi_Ncn,   & ! Prescribed correlation of chi (old s) & Ncn  [-]
      precip_frac_tol         ! Min. precip. frac. when hydromet. present    [-]

    real( kind = core_rknd ), dimension(nz,hydromet_dim) :: &
      wphydrometp_chnge    ! Change in wphydrometp_zt: covar. clip. [(m/s)units]

    real( kind = core_rknd ), dimension(nz) :: &
      wm_zt,  & ! Mean vertical velocity, <w>, on thermo. levels  [m/s]
      wp2_zt    ! Variance of w, <w'^2> (interp. to t-levs.)      [m^2/s^2]

    real( kind = core_rknd ), dimension(nz) :: &
      rtp2_zt_from_chi

    logical :: l_corr_array_scaling

    ! Flags used for covariance clipping of <w'hm'>.
    logical, parameter :: &
      l_first_clip_ts = .true., & ! First instance of clipping in a timestep.
      l_last_clip_ts  = .true.    ! Last instance of clipping in a timestep.

    character(len=10) :: &
      hydromet_name    ! Name of a hydrometeor

    integer :: k, i  ! Loop indices

    ! ---- Begin Code ----

    ! Assertion check
    ! Check that all hydrometeors are positive otherwise exit the program
    if ( clubb_at_least_debug_level( 2 ) ) then
       do i = 1, hydromet_dim
          if ( any( hydromet(:,i) < zero_threshold ) ) then
             hydromet_name = hydromet_list(i)
             do k = 1, nz
                if ( hydromet(k,i) < zero_threshold ) then

                   ! Write error message
                   write(fstderr,*) trim( hydromet_name )//" = ", &
                                    hydromet(k,i), " < ", zero_threshold, &
                                    " at beginning of setup_pdf_parameters" &
                                    //" at k = ", k

                   ! Exit program
                   stop "Exiting..."

                endif ! hydromet(k,i) < 0
             enddo ! k = 1, nz
          endif ! hydromet(:,i) < 0
       enddo ! i = 1, hydromet_dim

    endif !clubb_at_least_debug_level( 2 )

    ! Setup some of the PDF parameters
    mu_w_1       = pdf_params%w_1
    mu_w_2       = pdf_params%w_2
    mu_chi_1     = pdf_params%chi_1
    mu_chi_2     = pdf_params%chi_2
    sigma_w_1    = sqrt( pdf_params%varnce_w_1 )
    sigma_w_2    = sqrt( pdf_params%varnce_w_2 )
    sigma_chi_1  = pdf_params%stdev_chi_1
    sigma_chi_2  = pdf_params%stdev_chi_2
    rc_1         = pdf_params%rc_1
    rc_2         = pdf_params%rc_2
    cloud_frac_1 = pdf_params%cloud_frac_1
    cloud_frac_2 = pdf_params%cloud_frac_2
    mixt_frac    = pdf_params%mixt_frac

    ice_supersat_frac_1 = pdf_params%ice_supersat_frac_1
    ice_supersat_frac_2 = pdf_params%ice_supersat_frac_2

    ! Recalculate wm_zt and wp2_zt.  Mean vertical velocity may not be easy to
    ! pass into this subroutine from a host model, and wp2_zt needs to have a
    ! value consistent with the value it had when the PDF parameters involving w
    ! were originally set in subroutine pdf_closure.  The variable wp2 has since
    ! been advanced, resulting a new wp2_zt.  However, the value of wp2 here
    ! needs to be consistent with wp2 at the time the PDF parameters were
    ! calculated.
    do k = 1, nz, 1

       ! Calculate the overall mean of vertical velocity, w, on thermodynamic
       ! levels.
       wm_zt(k) = compute_mean_binormal( mu_w_1(k), mu_w_2(k), mixt_frac(k) )

       ! Calculate the overall variance of vertical velocity on thermodynamic
       ! levels.
       wp2_zt(k) = compute_variance_binormal( wm_zt(k), mu_w_1(k), mu_w_2(k), &
                                              sigma_w_1(k), sigma_w_2(k), &
                                              mixt_frac(k) )

    enddo

    ! Note on hydrometeor PDF shape:
    ! To use a single lognormal over the entire grid level, turn off the
    ! l_use_precip_frac flag and set omicron to 1 and zeta_vrnce_rat to 0 in
    ! tunable_parameters.in.
    ! To use a single delta-lognormal (single lognormal in-precip.), enable the
    ! l_use_precip_frac flag and set omicron to 1 and zeta_vrnce_rat to 0 in
    ! tunable_parameters.in.
    ! Otherwise, with l_use_precip_frac enabled and omicron and zeta_vrnce_rat
    ! values that are not 1 and 0, respectively, the PDF shape is a double
    ! delta-lognormal (two independent lognormals in-precip.).

    ! Calculate precipitation fraction.
    if ( l_use_precip_frac ) then

       call precip_fraction( nz, hydromet, cloud_frac, cloud_frac_1, &    ! In
                             cloud_frac_2, ice_supersat_frac, &           ! In
                             ice_supersat_frac_1, ice_supersat_frac_2, &  ! In
                             mixt_frac, l_stats_samp, &                   ! In
                             precip_frac, precip_frac_1, precip_frac_2, & ! Out
                             precip_frac_tol )                            ! Out

    else

       precip_frac     = one
       precip_frac_1   = one
       precip_frac_2   = one
       precip_frac_tol = cloud_frac_min

    endif

    ! Calculate <N_cn> from Nc_in_cloud, whether Nc_in_cloud is predicted or
    ! based on a prescribed value, and whether the value is constant or varying
    ! over the grid level.
    if ( .not. l_const_Nc_in_cloud ) then
       ! Ncn varies at each vertical level.
       const_Ncnp2_on_Ncnm2 = Ncnp2_on_Ncnm2
    else  ! l_const_Nc_in_cloud
       ! Ncn is constant at each vertical level.
       const_Ncnp2_on_Ncnm2 = zero
    endif

    const_corr_chi_Ncn = corr_NN2NL( corr_array_n_cloud(iiPDF_Ncn, iiPDF_chi), &
                                     stdev_L2N( const_Ncnp2_on_Ncnm2 ), &
                                     const_Ncnp2_on_Ncnm2 )

    do k = 2, nz

       Ncnm(k) &
       = Nc_in_cloud_to_Ncnm( mu_chi_1(k), mu_chi_2(k), sigma_chi_1(k), &
                              sigma_chi_2(k), mixt_frac(k), Nc_in_cloud(k), &
                              cloud_frac_1(k), cloud_frac_2(k), &
                              const_Ncnp2_on_Ncnm2, const_corr_chi_Ncn )

    enddo ! k = 2, nz

    ! Boundary Condition.
    ! At thermodynamic level k = 1, which is below the model lower boundary, the
    ! value of Ncnm does not matter.
    Ncnm(1) = Nc_in_cloud(1)

    ! Calculate the overall variance of a precipitating hydrometeor (hm),
    !<hm'^2>.
    do i = 1, hydromet_dim, 1

       do k = 1, nz, 1
          if ( hydromet(k,i) >= hydromet_tol(i) ) then
             ! There is some of the hydrometeor species found at level k.
             ! Calculate the variance (overall) of the hydrometeor.
             hydrometp2_zt(k,i) &
             = ( ( hmp2_ip_on_hmm2_ip(i) + one ) / precip_frac(k) - one ) &
               * hydromet(k,i)**2
          else
             hydrometp2_zt(k,i) = zero
          endif
       enddo ! k = 1, nz, 1

       ! Statistics
       if ( l_stats_samp ) then
          if ( ihmp2_zt(i) > 0 ) then
             ! Variance (overall) of the hydrometeor, <hm'^2>.
             call stat_update_var( ihmp2_zt(i), hydrometp2_zt(:,i), stats_zt )
          endif
       endif ! l_stats_samp

       ! Interpolate the covariances (overall) of w and precipitating
       ! hydrometeors to thermodynamic grid levels.
       wphydrometp_zt(:,i) = zm2zt( wphydrometp(:,i) )

       ! When the mean value of a precipitating hydrometeor is below tolerance
       ! value, it is considered to have a value of 0, and the precipitating
       ! hydrometeor does not vary over the grid level.  Any covariances
       ! involving that precipitating hydrometeor also have values of 0 at that
       ! grid level.
       do k = 1, nz, 1

          if ( hydromet(k,i) < hydromet_tol(i) ) then
             wphydrometp_zt(k,i) = zero
          endif

          ! Clip the value of covariance <w'hm'> on thermodynamic levels.
          call clip_covar_level( clip_wphydrometp, k, l_first_clip_ts, &
                                 l_last_clip_ts, dt, wp2_zt(k), &
                                 hydrometp2_zt(k,i), &
                                 wphydrometp_zt(k,i), wphydrometp_chnge(k,i) )

       enddo ! k = 1, nz, 1

    enddo ! i = 1, hydromet_dim, 1

    ! Calculate correlations involving w by first calculating total covariances
    ! involving w (<w'r_r'>, etc.) using the down-gradient approximation.
    if ( l_calc_w_corr ) then

       ! Calculate the covariances of w with the hydrometeors
       do k = 1, nz
          wpchip_zm(k) = pdf_params(k)%mixt_frac &
                       * ( one - pdf_params(k)%mixt_frac ) &
                       * ( pdf_params(k)%chi_1 - pdf_params(k)%chi_2 ) &
                       * ( pdf_params(k)%w_1 - pdf_params(k)%w_2 )
       enddo

       wpNcnp_zm(1:nz-1) = xpwp_fnc( -c_K_hm * Kh_zm(1:nz-1), Ncnm(1:nz-1), &
                                     Ncnm(2:nz), gr%invrs_dzm(1:nz-1) )

       ! Boundary conditions; We are assuming zero flux at the top.
       wpNcnp_zm(nz) = zero

       ! Interpolate the covariances to thermodynamic grid levels.
       wpchip_zt = zm2zt( wpchip_zm )
       wpNcnp_zt = zm2zt( wpNcnp_zm )

       ! When the mean value of Ncn is below tolerance value, it is considered
       ! to have a value of 0, and Ncn does not vary over the grid level.  Any
       ! covariance involving Ncn also has a value of 0 at that grid level.
       do k = 1, nz, 1
          if ( Ncnm(k) <= Ncn_tol ) then
             wpNcnp_zt(k) = zero
          endif
       enddo ! k = 1, nz, 1

    endif ! l_calc_w_corr

    ! Statistics
    if ( l_stats_samp ) then

       if ( iprecip_frac > 0 ) then
          ! Overall precipitation fraction.
          call stat_update_var( iprecip_frac, precip_frac, stats_zt )
       endif

       if ( iprecip_frac_1 > 0 ) then
          ! Precipitation fraction in PDF component 1.
          call stat_update_var( iprecip_frac_1, precip_frac_1, stats_zt )
       endif

       if ( iprecip_frac_2 > 0 ) then
          ! Precipitation fraction in PDF component 2.
          call stat_update_var( iprecip_frac_2, precip_frac_2, stats_zt )
       endif

       if ( iNcnm > 0 ) then
          ! Mean simplified cloud nuclei concentration (overall).
          call stat_update_var( iNcnm, Ncnm, stats_zt )
       endif

    endif


    !!! Setup PDF parameters loop.
    ! Loop over all model thermodynamic level above the model lower boundary.
    ! Now also including "model lower boundary" -- Eric Raut Aug 2013
    ! Now not  including "model lower boundary" -- Eric Raut Aug 2014
    do k = 2, nz, 1

       !!! Calculate the means and standard deviations involving PDF variables
       !!! -- w, chi, eta, N_cn, and any precipitating hydrometeors (hm
       !!! in-precip) -- for each PDF component.
       call compute_mean_stdev( hydromet(k,:), hydrometp2_zt(k,:),     & ! In
                                Ncnm(k), mixt_frac(k), precip_frac(k), & ! In
                                precip_frac_1(k), precip_frac_2(k),    & ! In
                                precip_frac_tol,                       & ! In
                                pdf_params(k), d_variables,            & ! In
                                mu_x_1, mu_x_2,                        & ! Out
                                sigma_x_1, sigma_x_2,                  & ! Out
                                hm_1(k,:), hm_2(k,:),                  & ! Out
                                sigma2_on_mu2_ip_1,                    & ! Out
                                sigma2_on_mu2_ip_2                     ) ! Out

       !!! Transform the component means and standard deviations involving
       !!! precipitating hydrometeors (hm in-precip) and N_cn -- ln hm and
       !!! ln N_cn -- to normal space for each PDF component.
       call norm_transform_mean_stdev( hm_1(k,:), hm_2(k,:), &
                                       Ncnm(k), d_variables, &
                                       mu_x_1, mu_x_2, &
                                       sigma_x_1, sigma_x_2, &
                                       sigma2_on_mu2_ip_1, &
                                       sigma2_on_mu2_ip_2, &
                                       mu_x_1_n(:,k), mu_x_2_n(:,k), &
                                       sigma_x_1_n(:,k), sigma_x_2_n(:,k) )

       !!! Calculate the normal space correlations.
       !!! The normal space correlations are the same as the true correlations
       !!! except when at least one of the variables involved is a precipitating
       !!! hydrometeor or Ncn.  In these cases, the normal space correlation
       !!! involves the natural logarithm of the precipitating hydrometeors,
       !!! ln hm (for example, ln r_r and ln N_r), and ln N_cn for each PDF
       !!! component.
       if ( l_diagnose_correlations ) then

          if ( rcm(k) > rc_tol ) then

             call diagnose_correlations( d_variables, corr_array_n_cloud, & ! Intent(in)
                                         corr_array_1_n )                   ! Intent(out)

             call diagnose_correlations( d_variables, corr_array_n_cloud, & ! Intent(in)
                                         corr_array_2_n )                   ! Intent(out)

          else

             call diagnose_correlations( d_variables, corr_array_n_below, & ! Intent(in)
                                         corr_array_1_n )                   ! Intent(out)

             call diagnose_correlations( d_variables, corr_array_n_below, & ! Intent(in)
                                         corr_array_2_n )                   ! Intent(out)

          endif

       else ! if .not. l_diagnose_correlations

          call comp_corr_norm( wm_zt(k), rc_1(k), rc_2(k), cloud_frac_1(k), &
                               cloud_frac_2(k), wpchip_zt(k), wpNcnp_zt(k), &
                               sqrt(wp2_zt(k)), mixt_frac(k), precip_frac_1(k),&
                               precip_frac_2(k), wphydrometp_zt(k,:), &
                               mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &
                               sigma_x_1_n(:,k), sigma_x_2_n(:,k), &
                               corr_array_n_cloud, corr_array_n_below, &
                               pdf_params(k), d_variables, &
                               corr_array_1_n(:,:,k), corr_array_2_n(:,:,k) )

       endif ! l_diagnose_correlations

       !!! Calculate the true correlations for each PDF component.
       call denorm_transform_corr( d_variables, &
                                   sigma_x_1_n(:,k), sigma_x_2_n(:,k), &
                                   sigma2_on_mu2_ip_1, sigma2_on_mu2_ip_2, &
                                   corr_array_1_n(:,:,k), &
                                   corr_array_2_n(:,:,k), &
                                   corr_array_1, corr_array_2 )

       !!! Statistics for standard PDF parameters involving hydrometeors.
       call pdf_param_hm_stats( d_variables, k, hm_1(k,:), hm_2(k,:), &
                                mu_x_1, mu_x_2, &
                                sigma_x_1, sigma_x_2, &
                                corr_array_1, corr_array_2, &
                                l_stats_samp )

       !!! Statistics for normal space PDF parameters involving hydrometeors.
       call pdf_param_ln_hm_stats( d_variables, k, mu_x_1_n(:,k), &
                                   mu_x_2_n(:,k), sigma_x_1_n(:,k), &
                                   sigma_x_2_n(:,k), corr_array_1_n(:,:,k), &
                                   corr_array_2_n(:,:,k), l_stats_samp )

       !!! Pack the PDF parameters
       call pack_pdf_params( hm_1(k,:), hm_2(k,:), d_variables, &          ! In
                             mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &       ! In
                             corr_array_1, corr_array_2, precip_frac(k), & ! In
                             precip_frac_1(k), precip_frac_2(k), &         ! In
                             hydromet_pdf_params(k) )                      ! Out

       if ( l_diagnose_correlations ) then

          call calc_cholesky_corr_mtx_approx &
                         ( d_variables, corr_array_1_n(:,:,k), &           ! intent(in)
                           corr_cholesky_mtx_1(:,:,k), corr_mtx_approx_1 ) ! intent(out)

          call calc_cholesky_corr_mtx_approx &
                         ( d_variables, corr_array_2_n(:,:,k), &           ! intent(in)
                           corr_cholesky_mtx_2(:,:,k), corr_mtx_approx_2 ) ! intent(out)

          corr_array_1_n(:,:,k) = corr_mtx_approx_1
          corr_array_2_n(:,:,k) = corr_mtx_approx_2

       else

          ! Compute choleksy factorization for the correlation matrix (out of
          ! cloud)
          call Cholesky_factor( d_variables, corr_array_1_n(:,:,k), & ! In
                                corr_array_scaling, corr_cholesky_mtx_1(:,:,k), &  ! Out
                                l_corr_array_scaling ) ! Out

          call Cholesky_factor( d_variables, corr_array_2_n(:,:,k), & ! In
                                corr_array_scaling, corr_cholesky_mtx_2(:,:,k), &  ! Out
                                l_corr_array_scaling ) ! Out
       endif

       ! For ease of use later in the code, we make the correlation arrays
       ! symmetrical
       call mirror_lower_triangular_matrix( d_variables, corr_array_1_n(:,:,k) )
       call mirror_lower_triangular_matrix( d_variables, corr_array_2_n(:,:,k) )

    enddo  ! Setup PDF parameters loop: k = 2, nz, 1

    ! Interpolate the overall variance of a hydrometeor, <hm'^2>, to its home on
    ! momentum grid levels.
    do i = 1, hydromet_dim, 1
       hydrometp2(:,i)  = zt2zm( hydrometp2_zt(:,i) )
       hydrometp2(nz,i) = zero
    enddo

    if ( l_stats_samp ) then
       if ( irtp2_from_chi > 0 ) then
          rtp2_zt_from_chi &
          = compute_rtp2_from_chi( pdf_params(:), &
                                   corr_array_1_n(iiPDF_chi,iiPDF_eta,:), &
                                   corr_array_2_n(iiPDF_chi,iiPDF_eta,:) )
          call stat_update_var( irtp2_from_chi, zt2zm( rtp2_zt_from_chi ), &
                                stats_zm )
       endif
    endif


    ! Boundary conditions for the output variables at k=1.
    mu_x_1_n(:,1) = zero
    mu_x_2_n(:,1) = zero
    sigma_x_1_n(:,1) = zero
    sigma_x_2_n(:,1) = zero
    corr_array_1_n(:,:,1) = zero
    corr_array_2_n(:,:,1) = zero
    corr_cholesky_mtx_1(:,:,1) = zero
    corr_cholesky_mtx_2(:,:,1) = zero
    call init_hydromet_pdf_params( hydromet_pdf_params(1) )

    if (clubb_at_least_debug_level( 2 )) then
       do k = 2, nz
          call assert_corr_symmetric( corr_array_1_n(:,:,k), d_variables )
          call assert_corr_symmetric( corr_array_2_n(:,:,k), d_variables )
       enddo
    endif


    return

  end subroutine setup_pdf_parameters

  !=============================================================================
  subroutine compute_mean_stdev( hydromet, hydrometp2_zt,       & ! Intent(in)
                                 Ncnm, mixt_frac, precip_frac,  & ! Intent(in)
                                 precip_frac_1, precip_frac_2,  & ! Intent(in)
                                 precip_frac_tol,               & ! Intent(in)
                                 pdf_params, d_variables,       & ! Intent(in)
                                 mu_x_1, mu_x_2,                & ! Intent(out)
                                 sigma_x_1, sigma_x_2,          & ! Intent(out)
                                 hm_1, hm_2,                    & ! Intent(out)
                                 sigma_hm_1_sqd_on_mu_hm_1_sqd, & ! Intent(out)
                                 sigma_hm_2_sqd_on_mu_hm_2_sqd  ) ! Intent(out)

    ! Description:
    ! Calculates the means and standard deviations (for each PDF component) of
    ! chi, eta, w, Ncn, and the precipitating hydrometeors.  For the
    ! precipitating hydrometeors, the component means and standard deviations
    ! are in-precip.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one,  & ! Constant(s)
        zero

    use array_index, only: &
        hydromet_tol

    use model_flags, only: &
        l_const_Nc_in_cloud ! Variable(s)

    use index_mapping, only: &
        pdf2hydromet_idx  ! Procedure(s)

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s) type

    use parameters_tunable, only: &
        omicron,        & ! Variable(s)
        zeta_vrnce_rat

    use corr_varnce_module, only: &
        iiPDF_chi,          & ! Variable(s)
        iiPDF_eta,          &
        iiPDF_w,            &
        iiPDF_Ncn,          &
        hmp2_ip_on_hmm2_ip, &
        Ncnp2_on_Ncnm2

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(hydromet_dim), intent(in) :: &
      hydromet,       & ! Mean of a hydrometeor (overall)             [hm units]
      hydrometp2_zt     ! Variance of a hydrometeor (overall)     [(hm units)^2]

    real( kind = core_rknd ), intent(in) :: &
      Ncnm,            & ! Mean simplified cloud nuclei concentration   [num/kg]
      mixt_frac,       & ! Mixture fraction                                  [-]
      precip_frac,     & ! Precipitation fraction (overall)                  [-]
      precip_frac_1,   & ! Precipitation fraction (1st PDF component)        [-]
      precip_frac_2,   & ! Precipitation fraction (2nd PDF component)        [-]
      precip_frac_tol    ! Minimum precip. frac. when hydromet. are present  [-]

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                                [units vary]

    integer, intent(in) :: &
      d_variables    ! Number of PDF variables

    ! Output Variables
    ! Note:  This code assumes to be these arrays in the same order as the
    ! correlation arrays, etc., which is determined by the iiPDF indices.
    ! The order should be as follows:  chi, eta, w, Ncn, <precip. hydrometeors>
    ! (indices increasing from left to right).
    real( kind = core_rknd ), dimension(d_variables), intent(out) :: &
      mu_x_1,    & ! Mean array of PDF vars. (1st PDF component)    [units vary]
      mu_x_2,    & ! Mean array of PDF vars. (2nd PDF component)    [units vary]
      sigma_x_1, & ! Standard deviation array of PDF vars (comp. 1) [units vary]
      sigma_x_2    ! Standard deviation array of PDF vars (comp. 2) [units vary]

    real( kind = core_rknd ), dimension(hydromet_dim), intent(out) :: &
      hm_1, & ! Mean of a precip. hydrometeor (1st PDF component)   [units vary]
      hm_2    ! Mean of a precip. hydrometeor (2nd PDF component)   [units vary]

    real( kind = core_rknd ), dimension(d_variables), intent(out) :: &
      sigma_hm_1_sqd_on_mu_hm_1_sqd, & ! Ratio sigma_hm_1^2 / mu_hm_1^2      [-]
      sigma_hm_2_sqd_on_mu_hm_2_sqd    ! Ratio sigma_hm_2^2 / mu_hm_2^2      [-]

    ! Local Variables
    integer :: ivar ! Loop iterator

    integer :: hm_idx  ! Hydrometeor array index.


    !!! Initialize output variables.
    mu_x_1 = zero
    mu_x_2 = zero
    sigma_x_1 = zero
    sigma_x_2 = zero
    hm_1 = zero
    hm_2 = zero
    sigma_hm_1_sqd_on_mu_hm_1_sqd = zero
    sigma_hm_2_sqd_on_mu_hm_2_sqd = zero


    !!! Enter the PDF parameters.

    !!! Vertical velocity, w.

    ! Mean of vertical velocity, w, in PDF component 1.
    mu_x_1(iiPDF_w) = pdf_params%w_1

    ! Mean of vertical velocity, w, in PDF component 2.
    mu_x_2(iiPDF_w) = pdf_params%w_2

    ! Standard deviation of vertical velocity, w, in PDF component 1.
    sigma_x_1(iiPDF_w) = sqrt( pdf_params%varnce_w_1 )

    ! Standard deviation of vertical velocity, w, in PDF component 2.
    sigma_x_2(iiPDF_w) = sqrt( pdf_params%varnce_w_2 )


    !!! Extended liquid water mixing ratio, chi.

    ! Mean of extended liquid water mixing ratio, chi (old s),
    ! in PDF component 1.
    mu_x_1(iiPDF_chi) = pdf_params%chi_1

    ! Mean of extended liquid water mixing ratio, chi (old s),
    ! in PDF component 2.
    mu_x_2(iiPDF_chi) = pdf_params%chi_2

    ! Standard deviation of extended liquid water mixing ratio, chi (old s),
    ! in PDF component 1.
    sigma_x_1(iiPDF_chi) = pdf_params%stdev_chi_1

    ! Standard deviation of extended liquid water mixing ratio, chi (old s),
    ! in PDF component 2.
    sigma_x_2(iiPDF_chi) = pdf_params%stdev_chi_2


    !!! Coordinate orthogonal to chi, eta.

    ! Mean of eta (old t) in PDF component 1.
    ! Set the component mean values of eta to 0.
    ! The component mean values of eta are not important.  They can be set to
    ! anything.  They cancel out in the model code.  However, the best thing to
    ! do is to set them to 0 and avoid any kind of numerical error.
    mu_x_1(iiPDF_eta) = zero

    ! Mean of eta (old t) in PDF component 2.
    ! Set the component mean values of eta to 0.
    ! The component mean values of eta are not important.  They can be set to
    ! anything.  They cancel out in the model code.  However, the best thing to
    ! do is to set them to 0 and avoid any kind of numerical error.
    mu_x_2(iiPDF_eta) = zero

    ! Standard deviation of eta (old t) in PDF component 1.
    sigma_x_1(iiPDF_eta) = pdf_params%stdev_eta_1

    ! Standard deviation of eta (old t) in PDF component 2.
    sigma_x_2(iiPDF_eta) = pdf_params%stdev_eta_2


    !!! Simplified cloud nuclei concentration, Ncn.

    ! Mean of simplified cloud nuclei concentration, Ncn, in PDF component 1.
    mu_x_1(iiPDF_Ncn) = Ncnm

    ! Mean of simplified cloud nuclei concentration, Ncn, in PDF component 2.
    mu_x_2(iiPDF_Ncn) = Ncnm

    ! Standard deviation of simplified cloud nuclei concentration, Ncn,
    ! in PDF component 1.
    if ( .not. l_const_Nc_in_cloud ) then

       ! Ncn varies in both PDF components.
       sigma_x_1(iiPDF_Ncn) = sqrt( Ncnp2_on_Ncnm2 ) * Ncnm

       sigma_x_2(iiPDF_Ncn) = sqrt( Ncnp2_on_Ncnm2 ) * Ncnm

       ! Ncn is not an official hydrometeor.  However, both the
       ! sigma_hm_1_sqd_on_mu_hm_1_sqd and sigma_hm_2_sqd_on_mu_hm_2_sqd arrays
       ! have size d_variables, and both sigma_Ncn_1^2/mu_Ncn_1^2 and
       ! sigma_Ncn_2^2/mu_Ncn_2^2 need to be output as part of these arrays.
       sigma_hm_1_sqd_on_mu_hm_1_sqd(iiPDF_Ncn) = Ncnp2_on_Ncnm2
       sigma_hm_2_sqd_on_mu_hm_2_sqd(iiPDF_Ncn) = Ncnp2_on_Ncnm2

    else ! l_const_Nc_in_cloud

       ! Ncn is constant in both PDF components.
       sigma_x_1(iiPDF_Ncn) = zero

       sigma_x_2(iiPDF_Ncn) = zero

       ! Ncn is not an official hydrometeor.  However, both the
       ! sigma_hm_1_sqd_on_mu_hm_1_sqd and sigma_hm_2_sqd_on_mu_hm_2_sqd arrays
       ! have size d_variables, and both sigma_Ncn_1^2/mu_Ncn_1^2 and
       ! sigma_Ncn_2^2/mu_Ncn_2^2 need to be output as part of these arrays.
       sigma_hm_1_sqd_on_mu_hm_1_sqd(iiPDF_Ncn) = zero
       sigma_hm_2_sqd_on_mu_hm_2_sqd(iiPDF_Ncn) = zero

    endif ! .not. l_const_Nc_in_cloud


    !!! Precipitating hydrometeor species.
    do ivar = iiPDF_Ncn+1, d_variables, 1

       hm_idx = pdf2hydromet_idx(ivar)

       call calc_comp_mu_sigma_hm( hydromet(hm_idx), hydrometp2_zt(hm_idx), &
                                   hmp2_ip_on_hmm2_ip(hm_idx), &
                                   mixt_frac, precip_frac, &
                                   precip_frac_1, precip_frac_2, &
                                   hydromet_tol(hm_idx), precip_frac_tol, &
                                   omicron, zeta_vrnce_rat, &
                                   mu_x_1(ivar), mu_x_2(ivar), &
                                   sigma_x_1(ivar), sigma_x_2(ivar), &
                                   hm_1(hm_idx), hm_2(hm_idx), &
                                   sigma_hm_1_sqd_on_mu_hm_1_sqd(ivar), &
                                   sigma_hm_2_sqd_on_mu_hm_2_sqd(ivar) )

    enddo ! ivar = iiPDF_Ncn+1, d_variables, 1


    return

  end subroutine compute_mean_stdev

  !=============================================================================
  subroutine comp_corr_norm( wm_zt, rc_1, rc_2, cloud_frac_1, &
                             cloud_frac_2, wpchip, wpNcnp, &
                             stdev_w, mixt_frac, precip_frac_1, &
                             precip_frac_2, wphydrometp_zt, &
                             mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &
                             sigma_x_1_n, sigma_x_2_n, &
                             corr_array_n_cloud, corr_array_n_below, &
                             pdf_params, d_variables, &
                             corr_array_1_n, corr_array_2_n )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        Ncn_tol,      &
        w_tol,        & ! [m/s]
        chi_tol, & ! [kg/kg]
        one,          &
        zero

    use model_flags, only: &
        l_calc_w_corr

    use diagnose_correlations_module, only: &
        calc_mean,        & ! Procedure(s)
        calc_w_corr

    use index_mapping, only: &
        pdf2hydromet_idx  ! Procedure(s)

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s) type

    use corr_varnce_module, only: &
        iiPDF_chi, & ! Variable(s)
        iiPDF_eta, &
        iiPDF_w,        &
        iiPDF_Ncn

    use array_index, only: &
        hydromet_tol

    implicit none

    ! Input Variables
    integer, intent(in) :: d_variables ! Number of variables in the corr/mean/stdev arrays

    real( kind = core_rknd ), intent(in) :: &
      wm_zt,         & ! Mean vertical velocity, <w>, on thermo. levels    [m/s]
      rc_1,          & ! Mean of r_c (1st PDF component)                 [kg/kg]
      rc_2,          & ! Mean of r_c (2nd PDF component)                 [kg/kg]
      cloud_frac_1,  & ! Cloud fraction (1st PDF component)                  [-]
      cloud_frac_2,  & ! Cloud fraction (2nd PDF component)                  [-]
      wpchip,        & ! Covariance of w and chi (old s)            [(m/s)kg/kg]
      wpNcnp,        & ! Covariance of w and N_cn (overall)       [(m/s) num/kg]
      stdev_w,       & ! Standard deviation of w                           [m/s]
      mixt_frac,     & ! Mixture fraction                                    [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)          [-]

    real( kind = core_rknd ), dimension(hydromet_dim), intent(in) :: &
      wphydrometp_zt    ! Covariance of w and hm interp. to t-levs.  [(m/s)u.v.]

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      mu_x_1,      & ! Mean of x array (1st PDF component)          [units vary]
      mu_x_2,      & ! Mean of x array (2nd PDF component)          [units vary]
      sigma_x_1,   & ! Standard deviation of x array (1st PDF comp.)  [un. vary]
      sigma_x_2,   & ! Standard deviation of x array (2nd PDF comp.)  [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(d_variables, d_variables), &
    intent(in) :: &
      corr_array_n_cloud, & ! Prescribed correlation array in cloud        [-]
      corr_array_n_below    ! Prescribed correlation array below cloud     [-]

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                                [units vary]

    ! Output Variables
    real( kind = core_rknd ), dimension(d_variables, d_variables), &
    intent(out) :: &
      corr_array_1_n, & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
      corr_array_2_n    ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]

    ! Local Variables
    real( kind = core_rknd ), dimension(d_variables)  :: &
      corr_w_hm_1_n, & ! Correlation of w and ln hm (1st PDF component) ip   [-]
      corr_w_hm_2_n    ! Correlation of w and ln hm (2nd PDF component) ip   [-]

    real( kind = core_rknd ) :: &
      chi_m,          & ! Mean of chi (s_mellor)                         [kg/kg]
      stdev_chi,      & ! Standard deviation of chi (s_mellor)           [kg/kg]
      corr_w_chi,     & ! Correlation of w and chi (overall)                 [-]
      corr_w_Ncn_1_n, & ! Correlation of w and ln Ncn (1st PDF component)    [-]
      corr_w_Ncn_2_n    ! Correlation of w and ln Ncn (2nd PDF component)    [-]

    logical :: &
      l_limit_corr_chi_eta    ! Flag to limit the correlation of chi and eta [-]

    integer :: ivar, jvar ! Loop iterators

    ! ---- Begin Code ----

    !!! Normal space correlations

    ! Initialize corr_w_hm_1_n and corr_w_hm_2_n arrays to 0.
    corr_w_hm_1_n = zero
    corr_w_hm_2_n = zero

    ! Calculate normal space correlations involving w by first calculating total
    ! covariances involving w (<w'Ncn'>, etc.) using the down-gradient
    ! approximation.
    if ( l_calc_w_corr ) then

       ! Approximate the correlation between w and chi.
       chi_m &
       = calc_mean( pdf_params%mixt_frac, pdf_params%chi_1, pdf_params%chi_2 )

       stdev_chi &
       = sqrt( pdf_params%mixt_frac &
               * ( ( pdf_params%chi_1 - chi_m )**2 &
                   + pdf_params%stdev_chi_1**2 ) &
             + ( one - pdf_params%mixt_frac ) &
               * ( ( pdf_params%chi_2 - chi_m )**2 &
                   + pdf_params%stdev_chi_2**2 ) &
             )

       corr_w_chi &
       = calc_w_corr( wpchip, stdev_w, stdev_chi, w_tol, chi_tol )

       ! Calculate the correlation of w and ln Ncn in each PDF component.
       ! The subroutine calc_corr_w_hm_n can be used to do this as long as a
       ! value of 1 is sent in for precip_frac_1 and precip_frac_2.
       jvar = iiPDF_Ncn
       call calc_corr_w_hm_n( wm_zt, wpNcnp, &
                              mu_x_1(iiPDF_w), mu_x_2(iiPDF_w), &
                              mu_x_1(jvar), mu_x_2(jvar), &
                              sigma_x_1(iiPDF_w), sigma_x_2(iiPDF_w), &
                              sigma_x_1(jvar), sigma_x_2(jvar), &
                              sigma_x_1_n(jvar), sigma_x_2_n(jvar), &
                              mixt_frac, one, one, &
                              corr_w_Ncn_1_n, corr_w_Ncn_2_n, &
                              Ncn_tol )

       ! Calculate the correlation of w and the natural logarithm of the
       ! hydrometeor for each PDF component and each hydrometeor type.
       do jvar = iiPDF_Ncn+1, d_variables

          call calc_corr_w_hm_n( wm_zt, wphydrometp_zt(pdf2hydromet_idx(jvar)),&
                                 mu_x_1(iiPDF_w), mu_x_2(iiPDF_w), &
                                 mu_x_1(jvar), mu_x_2(jvar), &
                                 sigma_x_1(iiPDF_w), sigma_x_2(iiPDF_w), &
                                 sigma_x_1(jvar), sigma_x_2(jvar), &
                                 sigma_x_1_n(jvar), sigma_x_2_n(jvar), &
                                 mixt_frac, precip_frac_1, precip_frac_2, &
                                 corr_w_hm_1_n(jvar), corr_w_hm_2_n(jvar), &
                                 hydromet_tol(pdf2hydromet_idx(jvar)) )

       enddo ! jvar = iiPDF_Ncn+1, d_variables

    endif

    ! In order to decompose the normal space correlation matrix,
    ! we must not have a perfect correlation of chi and
    ! eta. Thus, we impose a limitation.
    l_limit_corr_chi_eta = .true.


    ! Initialize the normal space correlation arrays
    corr_array_1_n = zero
    corr_array_2_n = zero

    !!! The corr_arrays are assumed to be lower triangular matrices
    ! Set diagonal elements to 1
    do ivar=1, d_variables
      corr_array_1_n(ivar, ivar) = one
      corr_array_2_n(ivar, ivar) = one
    end do


    !!! This code assumes the following order in the prescribed correlation
    !!! arrays (iiPDF indices):
    !!! chi, eta, w, Ncn, <hydrometeors> (indices increasing from left to right)

    ! Correlation of chi (old s) and eta (old t)
    corr_array_1_n(iiPDF_eta, iiPDF_chi) &
    = component_corr_chi_eta( pdf_params%corr_chi_eta_1, rc_1, cloud_frac_1, &
                              corr_array_n_cloud(iiPDF_eta, iiPDF_chi), &
                              corr_array_n_below(iiPDF_eta, iiPDF_chi), &
                              l_limit_corr_chi_eta )

    corr_array_2_n(iiPDF_eta, iiPDF_chi) &
    = component_corr_chi_eta( pdf_params%corr_chi_eta_2, rc_2, cloud_frac_2, &
                              corr_array_n_cloud(iiPDF_eta, iiPDF_chi), &
                              corr_array_n_below(iiPDF_eta, iiPDF_chi), &
                              l_limit_corr_chi_eta )

    ! Correlation of chi (old s) and w
    corr_array_1_n(iiPDF_w, iiPDF_chi) &
    = component_corr_w_x( corr_w_chi, rc_1, cloud_frac_1, &
                          corr_array_n_cloud(iiPDF_w, iiPDF_chi), &
                          corr_array_n_below(iiPDF_w, iiPDF_chi) )

    corr_array_2_n(iiPDF_w, iiPDF_chi) &
    = component_corr_w_x( corr_w_chi, rc_2, cloud_frac_2, &
                          corr_array_n_cloud(iiPDF_w, iiPDF_chi), &
                          corr_array_n_below(iiPDF_w, iiPDF_chi) )


    ! Correlation of chi (old s) and ln Ncn
    corr_array_1_n(iiPDF_Ncn, iiPDF_chi) &
    = component_corr_x_hm_n_ip( rc_1, one, &
                                corr_array_n_cloud(iiPDF_Ncn, iiPDF_chi), &
                                corr_array_n_cloud(iiPDF_Ncn, iiPDF_chi) )

    corr_array_2_n(iiPDF_Ncn, iiPDF_chi) &
    = component_corr_x_hm_n_ip( rc_2, one, &
                                corr_array_n_cloud(iiPDF_Ncn, iiPDF_chi), &
                                corr_array_n_cloud(iiPDF_Ncn, iiPDF_chi) )

    ! Correlation of chi (old s) and the natural logarithm of the hydrometeors
    ivar = iiPDF_chi
    do jvar = iiPDF_Ncn+1, d_variables
       corr_array_1_n(jvar, ivar) &
       = component_corr_x_hm_n_ip( rc_1, cloud_frac_1,&
                                   corr_array_n_cloud(jvar, ivar), &
                                   corr_array_n_below(jvar, ivar) )

       corr_array_2_n(jvar, ivar) &
       = component_corr_x_hm_n_ip( rc_2, cloud_frac_2,&
                                   corr_array_n_cloud(jvar, ivar), &
                                   corr_array_n_below(jvar, ivar) )
    enddo

    ! Correlation of eta (old t) and w
    corr_array_1_n(iiPDF_w, iiPDF_eta) = zero
    corr_array_2_n(iiPDF_w, iiPDF_eta) = zero

    ! Correlation of eta (old t) and ln Ncn
    corr_array_1_n(iiPDF_Ncn, iiPDF_eta) &
    = component_corr_x_hm_n_ip( rc_1, one, &
                                corr_array_n_cloud(iiPDF_Ncn, iiPDF_eta), &
                                corr_array_n_cloud(iiPDF_Ncn, iiPDF_eta) )

    corr_array_2_n(iiPDF_Ncn, iiPDF_eta) &
    = component_corr_x_hm_n_ip( rc_2, one, &
                                corr_array_n_cloud(iiPDF_Ncn, iiPDF_eta), &
                                corr_array_n_cloud(iiPDF_Ncn, iiPDF_eta) )

    ! Correlation of eta (old t) and the natural logarithm of the hydrometeors
    ivar = iiPDF_eta
    do jvar = iiPDF_Ncn+1, d_variables
      corr_array_1_n(jvar, ivar) &
      = component_corr_eta_hm_n_ip( corr_array_1_n( iiPDF_eta, iiPDF_chi), &
                                    corr_array_1_n( jvar, iiPDF_chi) )

      corr_array_2_n(jvar, ivar) &
      = component_corr_eta_hm_n_ip( corr_array_2_n( iiPDF_eta, iiPDF_chi), &
                                    corr_array_2_n( jvar, iiPDF_chi) )
    enddo


    ! Correlation of w and ln Ncn
    corr_array_1_n(iiPDF_Ncn, iiPDF_w) &
    = component_corr_w_hm_n_ip( corr_w_Ncn_1_n, rc_1, one, &
                                corr_array_n_cloud(iiPDF_Ncn, iiPDF_w), &
                                corr_array_n_below(iiPDF_Ncn, iiPDF_w) )

    corr_array_2_n(iiPDF_Ncn, iiPDF_w) &
    = component_corr_w_hm_n_ip( corr_w_Ncn_2_n, rc_2, one, &
                                corr_array_n_cloud(iiPDF_Ncn, iiPDF_w), &
                                corr_array_n_below(iiPDF_Ncn, iiPDF_w) )

    ! Correlation of w and the natural logarithm of the hydrometeors
    ivar = iiPDF_w
    do jvar = iiPDF_Ncn+1, d_variables

       corr_array_1_n(jvar, ivar) &
       = component_corr_w_hm_n_ip( corr_w_hm_1_n(jvar), rc_1, cloud_frac_1, &
                                   corr_array_n_cloud(jvar, ivar), &
                                   corr_array_n_below(jvar, ivar) )

       corr_array_2_n(jvar, ivar) &
       = component_corr_w_hm_n_ip( corr_w_hm_2_n(jvar), rc_2, cloud_frac_2, &
                                   corr_array_n_cloud(jvar, ivar), &
                                   corr_array_n_below(jvar, ivar) )

    enddo

    ! Correlation of ln Ncn and the natural logarithm of the hydrometeors
    ivar = iiPDF_Ncn
    do jvar = iiPDF_Ncn+1, d_variables
       corr_array_1_n(jvar, ivar) &
       = component_corr_hmx_hmy_n_ip( rc_1, cloud_frac_1, &
                                      corr_array_n_cloud(jvar, ivar), &
                                      corr_array_n_below(jvar, ivar) )

       corr_array_2_n(jvar, ivar) &
       = component_corr_hmx_hmy_n_ip( rc_2, cloud_frac_2, &
                                      corr_array_n_cloud(jvar, ivar), &
                                      corr_array_n_below(jvar, ivar) )
    enddo

    ! Correlation of the natural logarithm of two hydrometeors
    do ivar = iiPDF_Ncn+1, d_variables-1
       do jvar = ivar+1, d_variables

          corr_array_1_n(jvar, ivar) &
          = component_corr_hmx_hmy_n_ip( rc_1, cloud_frac_1, &
                                         corr_array_n_cloud(jvar, ivar), &
                                         corr_array_n_below(jvar, ivar) )

          corr_array_2_n(jvar, ivar) &
          = component_corr_hmx_hmy_n_ip( rc_2, cloud_frac_2, &
                                         corr_array_n_cloud(jvar, ivar), &
                                         corr_array_n_below(jvar, ivar) )

       enddo ! jvar
    enddo ! ivar


    return

  end subroutine comp_corr_norm

  !=============================================================================
  subroutine calc_comp_mu_sigma_hm( hmm, hmp2, &                     ! In
                                    hmp2_ip_on_hmm2_ip, &            ! In
                                    mixt_frac, precip_frac, &        ! In
                                    precip_frac_1, precip_frac_2, &  ! In
                                    hm_tol, precip_frac_tol, &       ! In
                                    omicron, zeta_vrnce_rat, &       ! In
                                    mu_hm_1, mu_hm_2, &              ! Out
                                    sigma_hm_1, sigma_hm_2, &        ! Out
                                    hm_1, hm_2, &                    ! Out
                                    sigma_hm_1_sqd_on_mu_hm_1_sqd, & ! Out
                                    sigma_hm_2_sqd_on_mu_hm_2_sqd )  ! Out

    ! Description:
    ! When precipitation is found in both PDF components (precip_frac_1 > 0 and
    ! precip_frac_2 > 0), the method that solves for in-precip. mean and
    ! in-precip. standard deviation in each PDF component, preserving overall
    ! mean and overall variance, is used.  When precipitation fraction is found
    ! in one PDF component but not the other one (precip_frac_1 > 0 and
    ! precip_frac_2 = 0, or precip_frac_1 = 0 and precip_frac_2 > 0), the
    ! calculation of component in-precip. mean and in-precip. standard deviation
    ! is simple.  When precipitation is not found in either component
    ! (precip_frac_1 = 0 and precip_frac_2 = 0), there isn't any precipitation
    ! found overall (at that grid level).

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,  & ! Constant(s)
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      hmm,                & ! Hydrometeor mean (overall), <hm>        [hm units]
      hmp2,               & ! Hydrometeor variance (overall), <hm'^2> [hm un.^2]
      hmp2_ip_on_hmm2_ip, & ! Ratio <hm|_ip'^2> / <hm|_ip>^2                 [-]
      mixt_frac,          & ! Mixture fraction                               [-]
      precip_frac,        & ! Precipitation fraction (overall)               [-]
      precip_frac_1,      & ! Precipitation fraction (1st PDF component)     [-]
      precip_frac_2,      & ! Precipitation fraction (2nd PDF component)     [-]
      hm_tol,             & ! Tolerance value of hydrometeor          [hm units]
      precip_frac_tol       ! Min. precip. frac. when hydromet. are present  [-]

    real( kind = core_rknd ), intent(in) :: &
      omicron,        & ! Relative width parameter, omicron = R / Rmax       [-]
      zeta_vrnce_rat    ! Width parameter for sigma_hm_1^2 / mu_hm_1^2       [-]


    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      mu_hm_1,    & ! Mean of hm (1st PDF component) in-precip (ip)   [hm units]
      mu_hm_2,    & ! Mean of hm (2nd PDF component) ip               [hm units]
      sigma_hm_1, & ! Standard deviation of hm (1st PDF component) ip [hm units]
      sigma_hm_2, & ! Standard deviation of hm (2nd PDF component) ip [hm units]
      hm_1,       & ! Mean of hm (1st PDF component)                  [hm units]
      hm_2          ! Mean of hm (2nd PDF component)                  [hm units]

    real( kind = core_rknd ), intent(out) :: &
      sigma_hm_1_sqd_on_mu_hm_1_sqd, & ! Ratio sigma_hm_1**2 / mu_hm_1**2    [-]
      sigma_hm_2_sqd_on_mu_hm_2_sqd    ! Ratio sigma_hm_2**2 / mu_hm_2**2    [-]


    if ( hmm >= hm_tol &
         .and. precip_frac_1 >= precip_frac_tol &
         .and. precip_frac_2 >= precip_frac_tol ) then

       ! Precipitation is found in both PDF components.
       call calc_mu_sigma_two_comps( hmm, hmp2, hmp2_ip_on_hmm2_ip, &
                                     mixt_frac, precip_frac, precip_frac_1, &
                                     precip_frac_2, hm_tol, &
                                     omicron, zeta_vrnce_rat, &
                                     mu_hm_1, mu_hm_2, sigma_hm_1, &
                                     sigma_hm_2, hm_1, hm_2, &
                                     sigma_hm_1_sqd_on_mu_hm_1_sqd, &
                                     sigma_hm_2_sqd_on_mu_hm_2_sqd )


    elseif ( hmm >= hm_tol .and. precip_frac_1 >= precip_frac_tol ) then

       ! Precipitation is found in the 1st PDF component, but not in the 2nd
       ! PDF component (precip_frac_2 = 0).
       mu_hm_1 = hmm / ( mixt_frac * precip_frac_1 )
       mu_hm_2 = zero

       sigma_hm_1 = sqrt( max( ( hmp2 + hmm**2 &
                                 - mixt_frac * precip_frac_1 * mu_hm_1**2 ) &
                               / ( mixt_frac * precip_frac_1 ), &
                               zero ) )
       sigma_hm_2 = zero

       hm_1 = mu_hm_1 * precip_frac_1
       hm_2 = zero

       sigma_hm_1_sqd_on_mu_hm_1_sqd = sigma_hm_1**2 / mu_hm_1**2 
       ! The ratio sigma_hm_2^2 / mu_hm_2^2 is undefined.
       sigma_hm_2_sqd_on_mu_hm_2_sqd = zero


    elseif ( hmm >= hm_tol .and. precip_frac_2 >= precip_frac_tol ) then

       ! Precipitation is found in the 2nd PDF component, but not in the 1st
       ! PDF component (precip_frac_1 = 0).
       mu_hm_1 = zero
       mu_hm_2 = hmm / ( ( one - mixt_frac ) * precip_frac_2 )

       sigma_hm_1 = zero
       sigma_hm_2 &
       = sqrt( max( ( hmp2 + hmm**2 &
                      - ( one - mixt_frac ) * precip_frac_2 * mu_hm_2**2 ) &
                    / ( ( one - mixt_frac ) * precip_frac_2 ), &
                    zero ) )

       hm_1 = zero
       hm_2 = mu_hm_2 * precip_frac_2

       ! The ratio sigma_hm_1^2 / mu_hm_1^2 is undefined.
       sigma_hm_1_sqd_on_mu_hm_1_sqd = zero
       sigma_hm_2_sqd_on_mu_hm_2_sqd = sigma_hm_2**2 / mu_hm_2**2


    else ! hm < hm_tol or ( precip_frac_1 = 0 and precip_frac_2 = 0 ).

       ! Precipitation is not found in either PDF component.
       mu_hm_1 = zero
       mu_hm_2 = zero

       sigma_hm_1 = zero
       sigma_hm_2 = zero

       hm_1 = zero
       hm_2 = zero

       ! The ratio sigma_hm_1^2 / mu_hm_1^2 is undefined.
       sigma_hm_1_sqd_on_mu_hm_1_sqd = zero
       ! The ratio sigma_hm_2^2 / mu_hm_2^2 is undefined.
       sigma_hm_2_sqd_on_mu_hm_2_sqd = zero


    endif ! hmm >= hm_tol and precip_frac_1 >= precip_frac_tol
          ! and precip_frac_2 >= precip_frac_tol


    return

  end subroutine calc_comp_mu_sigma_hm

  !=============================================================================
  subroutine calc_mu_sigma_two_comps( hmm, hmp2, hmp2_ip_on_hmm2_ip, &
                                      mixt_frac, precip_frac, precip_frac_1, &
                                      precip_frac_2, hm_tol, &
                                      omicron, zeta_vrnce_rat, &
                                      mu_hm_1, mu_hm_2, sigma_hm_1, &
                                      sigma_hm_2, hm_1, hm_2, &
                                      sigma_hm_1_sqd_on_mu_hm_1_sqd, &
                                      sigma_hm_2_sqd_on_mu_hm_2_sqd )

    ! Description:
    !
    ! OVERVIEW
    !
    ! The goal is to calculate the in-precip. mean of the hydrometeor field in
    ! each PDF component (mu_hm_1 and mu_hm_2) in a scenario when there is
    ! precipitation found in both PDF components.  The fields provided are the
    ! overall mean of the hydrometeor, <hm>, the overall variance of the
    ! hydrometeor, <hm’^2>, the mixture fraction, a, the overall precipitation
    ! fraction, f_p, and the precipitation fraction in each PDF component
    ! (f_p_1 and f_p_2).
    !
    ! The PDF equation for <hm> is:
    !
    ! <hm> = a * f_p_1 * mu_hm_1 + ( 1- a ) * f_p_2 * mu_hm_2.
    !
    ! Likewise, the PDF equation for <hm’^2> is:
    !
    ! <hm’^2> = a * f_p_1 * ( mu_hm_1^2 + sigma_hm_1^2 )
    !           + ( 1 - a ) * f_p_2 * ( mu_hm_2^2 + sigma_hm_2^2 )
    !           - <hm>^2;
    !
    ! where sigma_hm_1 and sigma_hm_2 are the in-precip. standard deviations of
    ! the hydrometeor field in each PDF component.  This can be rewritten as:
    !
    ! <hm’^2>
    ! = a * f_p_1 * ( 1 + sigma_hm_1^2 / mu_hm_1^2 ) * mu_hm_1^2
    !   + ( 1 - a ) * f_p_2 * ( 1 + sigma_hm_2^2 / mu_hm_2^2 ) * mu_hm_2^2
    !   - <hm>^2.
    !
    ! The ratio of sigma_hm_2^2 to mu_hm_2^2 is denoted R:
    !
    ! R = sigma_hm_2^2 / mu_hm_2^2.
    !
    ! In order to allow sigma_hm_1^2 / mu_hm_1^2 to have a different ratio, the
    ! parameter zeta is introduced, such that:
    !
    ! R * ( 1 + zeta ) = sigma_hm_1^2 / mu_hm_1^2;
    !
    ! where zeta > -1.  When -1 < zeta < 0, the ratio sigma_hm_2^2 / mu_hm_2^2
    ! grows at the expense of sigma_hm_1^2 / mu_hm_1^2, which narrows.  When
    ! zeta = 0, the ratio sigma_hm_1^2 / mu_hm_1^2 is the same as
    ! sigma_hm_2^2 / mu_hm_2^2.  When zeta > 0, sigma_hm_1^2 / mu_hm_1^2 grows
    ! at the expense of sigma_hm_2^2 / mu_hm_2^2, which narrows.  The component
    ! variances are written as:
    !
    ! sigma_hm_1^2 = R * ( 1 + zeta ) * mu_hm_1^2; and
    ! sigma_hm_2^2 = R * mu_hm_2^2,
    !
    ! and the component standard deviations are simply:
    !
    ! sigma_hm_1 = sqrt( R * ( 1 + zeta ) ) * mu_hm_1; and
    ! sigma_hm_2 = sqrt( R ) * mu_hm_2.
    !
    ! The equation for <hm’^2> can be rewritten as:
    !
    ! <hm’^2> = a * f_p_1 * ( 1 + R * ( 1 + zeta ) ) * mu_hm_1^2
    !           + ( 1 - a ) * f_p_2 * ( 1 + R ) * mu_hm_2^2
    !           - <hm>^2.
    !
    !
    ! HYDROMETEOR IN-PRECIP. VARIANCE:
    ! THE SPREAD OF THE MEANS VS. THE STANDARD DEVIATIONS
    !
    ! Part I:  Minimum and Maximum Values for R
    !
    ! The in-precip. variance of the hydrometeor is accounted for through a
    ! combination of the variance of each PDF component and the spread between
    ! the means of each PDF component.  At one extreme, the standard deviation
    ! of each component could be set to 0 and the in-prccip. variance could be
    ! accounted for by spreading the PDF component (in-precip.) means far apart.
    ! The value of R in this scenario would be its minimum possible value, which
    ! is 0.  At the other extreme, the means of each component could be set
    ! equal to each other and the in-precip. variance could be accounted for
    ! entirely by the PDF component (in-precip.) standard deviations.  The value
    ! of R in this scenario would be its maximum possible value, which is Rmax.
    !
    ! In order to calculate the value of Rmax, use the equation set but set
    ! mu_hm_1 = mu_hm_2 and R = Rmax.  When this happens:
    !
    ! <hm> = ( a * f_p_1 + ( 1- a ) * f_p_2 ) * mu_hm_i;
    !
    ! and since f_p = a * f_p_1 + ( 1 - a ) * f_p_2:
    !
    ! mu_hm_i = <hm> / f_p = <hm|_ip>;
    !
    ! where <hm|_ip> is the in-precip. mean of the hydrometeor.  The equation
    ! for hydrometeor variance in this scenario becomes:
    !
    ! <hm’^2> = <hm|_ip>^2 * ( a * f_p_1 * ( 1 + Rmax * ( 1 + zeta ) )
    !                          + ( 1 - a ) * f_p_2 * ( 1 + Rmax ) )
    !           - <hm>^2.
    !
    ! The general equation for the in-precip. variance of a hydrometeor,
    ! <hm|_ip’^2>, is given by:
    !
    ! <hm|_ip’^2> = ( <hm’^2> + <hm>^2 - f_p * <hm|_ip>^2 ) / f_p;
    !
    ! which can be rewritten as:
    !
    ! <hm’^2> + <hm>^2 = f_p * ( <hm|_ip’^2> + <hm|_ip>^2 ).
    !
    ! When the above equation is substituted into the modified PDF equation for
    ! <hm’^2>, Rmax is solved for and the equation is:
    !
    ! Rmax = ( f_p / ( a * f_p_1 * ( 1 + zeta ) + ( 1 - a ) * f_p_2 ) )
    !        * ( <hm|_ip’^2> / <hm|_ip>^2 ).
    !
    ! Here, in the scenario that zeta = 0, both PDF components have the same
    ! mean and same variance, which reduces the in-precip. distribution to an
    ! assumed single lognormal, and the above equation reduces to:
    !
    ! Rmax = <hm|_ip’^2> / <hm|_ip>^2;
    !
    ! which is what is expected in that case.
    !
    !
    ! Part II:  Enter omicron
    !
    ! A parameter is used to prescribe the ratio of R to its maximum value,
    ! Rmax.  The prescribed parameter is called omicron, where:
    !
    ! R = omicron * Rmax;
    !
    ! where 0 <= omicron <= 1.  When omicron = 0, the standard deviation of each
    ! PDF component is 0, and mu_hm_1 is spread as far away from mu_hm_2 as it
    ! needs to be to account for the in-precip. variance.  When omicron = 1,
    ! mu_hm_1 is equal to mu_hm_2, and the standard deviations of the PDF
    ! components account for all of the in-precip. variance (and when zeta = 0,
    ! the PDF shape is a single lognormal in-precip.).  At intermediate values
    ! of omicron, the means of each PDF component are somewhat spread and each
    ! PDF component has some width.  The modified parameters are listed below.
    !
    ! The ratio of sigma_hm_2^2 to mu_hm_2^2 is:
    !
    ! sigma_hm_2^2 / mu_hm_2^2 = omicron * Rmax;
    !
    ! and the ratio of sigma_hm_1^2 / mu_hm_1^2 is:
    !
    ! sigma_hm_1^2 / mu_hm_1^2 = omicron * Rmax * ( 1 + zeta ). 
    ! 
    ! The component variances are written as:
    !
    ! sigma_hm_1^2 = omicron * Rmax * ( 1 + zeta ) * mu_hm_1^2; and
    ! sigma_hm_2^2 = omicron * Rmax * mu_hm_2^2,
    !
    ! and the component standard deviations are simply:
    !
    ! sigma_hm_1 = sqrt( omicron * Rmax * ( 1 + zeta ) ) * mu_hm_1; and
    ! sigma_hm_2 = sqrt( omicron * Rmax ) * mu_hm_2.
    !
    ! The equation set becomes:
    !
    ! [1] <hm> = a * f_p_1 * mu_hm_1 + ( 1- a ) * f_p_2 * mu_hm_2; and
    !
    ! [2] <hm’^2>
    !     = a * f_p_1 * ( 1 + omicron * Rmax * ( 1 + zeta ) ) * mu_hm_1^2
    !       + ( 1 - a ) * f_p_2 * ( 1 + omicron * Rmax ) * mu_hm_2^2
    !       - <hm>^2.
    !
    !
    ! SOLVING THE EQUATION SET FOR MU_HM_1 AND MU_HM_2.
    !
    ! The above system of two equations can be solved for mu_hm_1 and mu_hm_2.
    ! All other quantities in the equation set are known quantities.  The
    ! equation for <hm> is rewritten to isolate mu_hm_2:
    !
    ! mu_hm_2 = ( <hm> - a * f_p_1 * mu_hm_1 ) / ( ( 1 - a ) * f_p_2 ).
    !
    ! The above equation is substituted into the equation for <hm’^2>.  The
    ! equation for <hm’^2> is rewritten, resulting in:
    !
    ! [ a * f_p_1 * ( 1 + omicron * Rmax * ( 1 + zeta ) )
    !   + a^2 * f_p_1^2 * ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 ) ]
    ! * mu_hm_1^2
    ! + [ - 2 * <hm> * a * f_p_1 * ( 1 + omicron * Rmax )
    !     / ( ( 1 - a ) * f_p_2 ) ] * mu_hm_1
    ! + [ - ( <hm’^2>
    !         + ( 1 - ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 ) )
    !           * <hm>^2 ) ]
    ! = 0.
    !
    ! This equation is of the form:
    !
    ! A * mu_hm_1^2 + B * mu_hm_1 + C = 0;
    !
    ! so the solution for mu_hm_1 is:
    !
    ! mu_hm_1 = ( -B +/- sqrt( B^2 - 4*A*C ) ) / (2*A);
    !
    ! where:
    !
    ! A = a * f_p_1 * ( 1 + omicron * Rmax * ( 1 + zeta ) )
    !     + a^2 * f_p_1^2 * ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 );
    !
    ! B = - 2 * <hm> * a * f_p_1 * ( 1 + omicron * Rmax )
    !     / ( ( 1 - a ) * f_p_2 );
    !
    ! and
    !
    ! C = - ( <hm’^2>
    !         + ( 1 - ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 ) )
    !           * <hm>^2 ).
    !
    ! The signs of the coefficients:
    !
    ! 1) coefficient A is always positive,
    ! 2) coefficient B is always negative (this means that -B is always
    !    positive), and
    ! 3) coefficient C can be positive, negative, or zero.
    !
    ! Since ( 1 - ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 ) ) * <hm>^2 is
    ! always negative and <hm’^2> is always positive, the sign of coefficient C
    ! depends on which term is greater in magnitude.
    !
    ! When <hm’^2> is greater, the sign of coefficient C is negative.  This
    ! means that -4*A*C is positive, which in turn means that
    ! sqrt( B^2 - 4*A*C ) is greater in magnitude than -B.  If the subtraction
    ! option of the +/- were to be chosen, the value of mu_hm_1 would be
    ! negative in this scenerio.  So the natural thing to do would be to always
    ! choose the addition option.  However, this method requires that mu_hm_1
    ! equals mu_hm_2 when omicron = 1.  When zeta >= 0, this happens when the
    ! addition option is chosen, but not when the subtraction option is chosen.
    ! However, when zeta < 0, this happens when the subtraction option is
    ! chosen, but not when the addition option is chosen.  So, the equation for
    ! mu_hm_1 becomes:
    !
    ! mu_hm_1 = ( -B + sqrt( B^2 - 4*A*C ) ) / (2*A); when zeta >= 0; and
    ! mu_hm_1 = ( -B - sqrt( B^2 - 4*A*C ) ) / (2*A); when zeta < 0.
    !
    ! Once this is set, of course:
    !
    ! mu_hm_2 = ( <hm> - a * f_p_1 * mu_hm_1 ) / ( ( 1 - a ) * f_p_2 ).
    !
    ! The system has been solved and the in-precip. PDF component means have
    ! been found!
    !
    !
    ! NOTES
    !
    ! Note 1:
    !
    ! The term B^2 - 4*A*C has been analyzed, and mathematically:
    !
    ! B^2 - 4*A*C >= 0
    !
    ! always holds true.  Additionally, the minimum value:
    !
    ! B^2 - 4*A*C = 0,
    !
    ! can only occur when omicron = 1 and zeta = 0 (or alternatively to
    ! zeta = 0, Rmax = 0, but this only occurs when <hm|_ip'^2> / <hm|_ip>^2 has
    ! a value of 0).
    !
    ! Numerically, when omicron = 1 and zeta = 0, B^2 - 4*A*C can produce very
    ! small (on the order of epsilon) negative values.  This is due to numerical
    ! round off error.  When this happens, the erroneous small, negative value
    ! of B^2 - 4*A*C is simply reset to the value it's supposed to have, which
    ! is 0.
    !
    !
    ! Note 2:
    !
    ! As the value of <hm|_ip'^2> / <hm|_ip>^2 increases and as the value of
    ! omicron decreases (narrowing the in-precip standard deviations and
    ! increasing the spread between the in-precip means), a situtation arises
    ! where the value of one of the component means will become negative.  This
    ! is because there is a limit to the amount of in-precip variance that can
    ! be represented by this kind of distribution.  In order to prevent
    ! out-of-bounds values of mu_hm_1 or mu_hm_2, lower limits will be
    ! declared, called mu_hm_1_min and mu_hm_2_min.  The value of the
    ! hydrometeor in-precip. component mean will be limited from going any
    ! smaller (or negative) at this value.  From there, the value of the other
    ! hydrometeor in-precip. component mean is easy to calculate.  Then, both
    ! values will be entered into the calculation of hydrometeor variance, which
    ! will be rewritten to solve for R.  Then, both the hydrometeor mean and
    ! hydrometeor variance will be preserved with a valid distribution.
    !
    ! In this emergency scenario, the value of R is:
    !
    ! R = ( <hm'^2> + <hm>^2 - a * f_p_1 * mu_hm_1^2
    !       - ( 1 - a ) * f_p_2 * mu_hm_2^2 )
    !     / ( a * f_p_1 * ( 1 + zeta ) * mu_hm_1^2
    !         + ( 1 - a ) * f_p_2 * mu_hm_2^2 ).
    !
    ! The minimum values of the in-precip. component means are bounded by:
    !
    ! mu_hm_1_min >= hm_tol / f_p_1; and
    ! mu_hm_2_min >= hm_tol / f_p_2.
    !
    ! These are set this way because hm_1 ( = mu_hm_1 * f_p_1 ) and
    ! hm_2 ( = mu_hm_2 * f_p_2 ) need to have values of at least hm_tol when
    ! precipitation is found in both PDF components.
    !
    ! However, an in-precip. component mean value of hm_tol / f_p_1 or
    ! hm_tol / f_p_2 often produces a distribution where one component centers
    ! around values that are too small to be a good match with data taken from
    ! Large Eddy Simulations (LES).  It is desirable to increase the minimum
    ! threshold of mu_hm_1 and mu_hm_2.
    !
    ! As the minimum threshold increases, the value of the in-precip. component
    ! mean that is from the component that is not being set to the minimum
    ! threshold decreases.  If the minimum threshold were to be boosted as high
    ! as <hm> / f_p (in most cases, <hm> / f_p >> hm_tol / f_p_i), both
    ! components would have a value of <hm> / f_p.  The minimum threshold should
    ! not be set this high.
    !
    ! Additionally, the minimum threshold for one in-precip. component mean
    ! cannot be set so high as to drive the other in-precip. component mean
    ! below hm_tol / f_p_i.  (This doesn't come into play unless <hm> is close
    ! to hm_tol.)  The upper limit for the in-precip. mean values are:
    !
    ! mu_hm_1|_(upper. lim.) = ( <hm> - ( 1 - a ) * f_p_2 * ( hm_tol / f_p_2 ) )
    !                          / ( a * f_p_1 ); and
    !
    ! mu_hm_2|_(upper. lim.) = ( <hm> - a * f_p_1 * ( hm_tol / f_p_1 ) )
    !                          / ( ( 1 - a ) * f_p_2 );
    !
    ! which reduces to:
    !
    ! mu_hm_1|_(upper. lim.) = ( <hm> - ( 1 - a ) * hm_tol ) / ( a * f_p_1 );
    ! and
    ! mu_hm_2|_(upper. lim.) = ( <hm> - a * hm_tol ) / ( ( 1 - a ) * f_p_2 ).
    !
    ! An appropriate minimum value for mu_hm_1 can be set by:
    !
    ! mu_hm_1_min = | min( hm_tol / f_p_1
    !               |      + mu_hm_min_coef * ( <hm> / f_p - hm_tol / f_p_1 ),
    !               |      ( <hm> - ( 1 - a ) * hm_tol ) / ( a * f_p_1 ) );
    !               |    where <hm> / f_p > hm_tol / f_p_1;
    !               | hm_tol / f_p_1;
    !               |    where <hm> / f_p <= hm_tol / f_p_1;
    !
    ! and similarly for mu_hm_2:
    !
    ! mu_hm_2_min = | min( hm_tol / f_p_2
    !               |      + mu_hm_min_coef * ( <hm> / f_p - hm_tol / f_p_2 ),
    !               |      ( <hm> - a * hm_tol ) / ( ( 1 - a ) * f_p_2 ) );
    !               |    where <hm> / f_p > hm_tol / f_p_2;
    !               | hm_tol / f_p_2;
    !               |    where <hm> / f_p <= hm_tol / f_p_2;
    !
    ! where mu_hm_min_coef is a coefficient that has a value
    ! 0 <= mu_hm_min_coef < 1.  When the value of mu_hm_min_coef is 0,
    ! mu_hm_1_min reverts to hm_tol / f_p_1 and mu_hm_2_min reverts to
    ! hm_tol / f_p_2.  An appropriate value for mu_hm_min_coef should be small,
    ! such as 0.01 - 0.05.
    !
    !
    ! Note 3:
    !
    ! When the value of zeta >= 0, the value of mu_hm_1 tends to be larger than
    ! the value of mu_hm_2.  Likewise when the value of zeta < 0, the value of
    ! mu_hm_2 tends to be larger than the value of mu_hm_1.  Since most cloud
    ! water and cloud fraction tends to be found in PDF component 1, it is
    ! advantageous to have the larger in-precip. component mean of the
    ! hydrometeor also found in PDF component 1.  The recommended value of zeta
    ! is a value greater than or equal to 0.

    ! References:
    !----------------------------------------------------------------------- 

    use constants_clubb, only: &
        four,    & ! Constant(s)
        two,     &
        one,     &
        zero,    &
        fstderr

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      hmm,                & ! Hydrometeor mean (overall), <hm>           [hm un]
      hmp2,               & ! Hydrometeor variance (overall), <hm'^2>  [hm un^2]
      hmp2_ip_on_hmm2_ip, & ! Ratio <hm|_ip'^2> / <hm|_ip>^2                 [-]
      mixt_frac,          & ! Mixture fraction                               [-]
      precip_frac,        & ! Precipitation fraction (overall)               [-]
      precip_frac_1,      & ! Precipitation fraction (1st PDF component)     [-]
      precip_frac_2,      & ! Precipitation fraction (2nd PDF component)     [-]
      hm_tol                ! Tolerance value of hydrometeor             [hm un]

    real( kind = core_rknd ), intent(in) :: &
      omicron,        & ! Relative width parameter, omicron = R / Rmax       [-]
      zeta_vrnce_rat    ! Width parameter for sigma_hm_1^2 / mu_hm_1^2       [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      mu_hm_1,    & ! Mean of hm (1st PDF component) in-precip (ip)      [hm un]
      mu_hm_2,    & ! Mean of hm (2nd PDF component) ip                  [hm un]
      sigma_hm_1, & ! Standard deviation of hm (1st PDF component) ip    [hm un]
      sigma_hm_2, & ! Standard deviation of hm (2nd PDF component) ip    [hm un]
      hm_1,       & ! Mean of hm (1st PDF component)                     [hm un]
      hm_2          ! Mean of hm (2nd PDF component)                     [hm un]

    real( kind = core_rknd ), intent(out) :: &
      sigma_hm_1_sqd_on_mu_hm_1_sqd, & ! Ratio sigma_hm_1**2 / mu_hm_1**2    [-]
      sigma_hm_2_sqd_on_mu_hm_2_sqd    ! Ratio sigma_hm_2**2 / mu_hm_2**2    [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      Rmax,       & ! Maximum possible value of ratio R                      [-]
      coef_A,     & ! Coefficient A in A*mu_hm_1^2 + B*mu_hm_1 + C = 0       [-]
      coef_B,     & ! Coefficient B in A*mu_hm_1^2 + B*mu_hm_1 + C = 0   [hm un]
      coef_C,     & ! Coefficient C in A*mu_hm_1^2 + B*mu_hm_1 + C = 0 [hm un^2]
      Bsqd_m_4AC    ! Value B^2 - 4*A*C in quadratic eqn. for mu_hm_1  [hm un^2]

    real( kind = core_rknd ) :: &
      mu_hm_1_min, & ! Minimum value of mu_hm_1 (precip. in both comps.) [hm un]
      mu_hm_2_min    ! Minimum value of mu_hm_2 (precip. in both comps.) [hm un]

    real( kind = core_rknd ), parameter :: &
      mu_hm_min_coef = 0.01_core_rknd  ! Coef. for mu_hm_1_min and mu_hm_2_min


    ! Calculate the value of Rmax.
    ! Rmax = ( f_p / ( a * f_p_1 * ( 1 + zeta ) + ( 1 - a ) * f_p_2 ) )
    !        * ( <hm|_ip’^2> / <hm|_ip>^2 ).
    ! The parameter zeta is written in the code as zeta_vrnce_rat.
    Rmax = ( precip_frac &
             / ( mixt_frac * precip_frac_1 * ( one + zeta_vrnce_rat ) &
                 + ( one - mixt_frac ) * precip_frac_2 ) ) &
           * hmp2_ip_on_hmm2_ip

    ! Calculate the value of coefficient A.
    ! A = a * f_p_1 * ( 1 + omicron * Rmax * ( 1 + zeta ) )
    !     + a^2 * f_p_1^2 * ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 ).
    coef_A = mixt_frac * precip_frac_1 &
             * ( one + omicron * Rmax * ( one + zeta_vrnce_rat ) ) &
             + mixt_frac**2 * precip_frac_1**2 &
               * ( one + omicron * Rmax ) &
               / ( ( one - mixt_frac ) * precip_frac_2 )

    ! Calculate the value of coefficient B.
    ! B = - 2 * <hm> * a * f_p_1 * ( 1 + omicron * Rmax )
    !     / ( ( 1 - a ) * f_p_2 ).
    coef_B = -two * hmm * mixt_frac * precip_frac_1 &
              * ( one + omicron * Rmax ) &
              / ( ( one - mixt_frac ) * precip_frac_2 )

    ! Calculate the value of coefficient C.
    ! C = - ( <hm’^2>
    !         + ( 1 - ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 ) )
    !           * <hm>^2 ).
    coef_C = - ( hmp2 + ( one &
                          - ( one + omicron * Rmax ) &
                            / ( ( one - mixt_frac ) * precip_frac_2 ) &
                        ) * hmm**2 )

    ! Calculate value of B^2 - 4*A*C.
    Bsqd_m_4AC = coef_B**2 - four * coef_A * coef_C

    ! Mathematically, the value of B^2 - 4*A*C cannot be less than 0.
    ! Numerically, this can happen when numerical round off error causes an
    ! epsilon-sized negative value.  When this happens, reset the value of
    ! B^2 - 4*A*C to 0.
    if ( Bsqd_m_4AC < zero ) then
       Bsqd_m_4AC = zero
    endif

    ! Calculate the mean (in-precip.) of the hydrometeor in the 1st PDF
    ! component.
    if ( zeta_vrnce_rat >= zero ) then
       mu_hm_1 = ( -coef_B + sqrt( Bsqd_m_4AC ) ) / ( two * coef_A )
    else
       mu_hm_1 = ( -coef_B - sqrt( Bsqd_m_4AC ) ) / ( two * coef_A )
    endif

    ! Calculate the mean (in-precip.) of the hydrometeor in the 2nd PDF
    ! component.
    mu_hm_2 = ( hmm - mixt_frac * precip_frac_1 * mu_hm_1 ) &
              / ( ( one - mixt_frac ) * precip_frac_2 )

    ! Calculate the value of the ratio R (which is sigma_hm_2^2 / mu_hm_2^2),
    ! where R = omicron * Rmax.  The name of the variable used for R is
    ! sigma_hm_2_sqd_on_mu_hm_2_sqd.
    sigma_hm_2_sqd_on_mu_hm_2_sqd = omicron * Rmax

    ! Calculate minimum allowable values for mu_hm_1 and mu_hm_2.
    if ( hmm / precip_frac > hm_tol / precip_frac_1 ) then
       mu_hm_1_min &
       = min( hm_tol / precip_frac_1 &
              + mu_hm_min_coef * ( hmm / precip_frac &
                                   - hm_tol / precip_frac_1 ), &
              ( hmm - ( one - mixt_frac ) * hm_tol ) &
              / ( mixt_frac * precip_frac_1 ) )
    else ! hmm / precip_frac <= hm_tol / precip_frac_1
       mu_hm_1_min = hm_tol / precip_frac_1
    endif
    if ( hmm / precip_frac > hm_tol / precip_frac_2 ) then
       mu_hm_2_min &
       = min( hm_tol / precip_frac_2 &
              + mu_hm_min_coef * ( hmm / precip_frac &
                                   - hm_tol / precip_frac_2 ), &
              ( hmm - mixt_frac * hm_tol ) &
              / ( ( one - mixt_frac ) * precip_frac_2 ) )
    else ! hmm / precip_frac <= hm_tol / precip_frac_2
       mu_hm_2_min = hm_tol / precip_frac_2
    endif

    ! Handle the "emergency" situation when the specified value of omicron is
    ! too small for the value of <hm|_ip'^2> / <hm|_ip>^2, resulting in a
    ! component mean that is too small (below tolerance value) or negative.
    if ( mu_hm_1 < mu_hm_1_min ) then

       ! Set the value of mu_hm_1 to the threshold positive value.
       mu_hm_1 = mu_hm_1_min

       ! Recalculate the mean (in-precip.) of the hydrometeor in the 2nd PDF
       ! component.
       mu_hm_2 = ( hmm - mixt_frac * precip_frac_1 * mu_hm_1 ) &
                 / ( ( one - mixt_frac ) * precip_frac_2 )

       ! Recalculate the value of R ( sigma_hm_2^2 / mu_hm_2^2 ) in this
       ! scenario.
       ! R = ( <hm'^2> + <hm>^2 - a * f_p_1 * mu_hm_1^2
       !       - ( 1 - a ) * f_p_2 * mu_hm_2^2 )
       !     / ( a * f_p_1 * ( 1 + zeta ) * mu_hm_1^2
       !         + ( 1 - a ) * f_p_2 * mu_hm_2^2 ).
       sigma_hm_2_sqd_on_mu_hm_2_sqd &
       = ( hmp2 + hmm**2 - mixt_frac * precip_frac_1 * mu_hm_1**2 &
           - ( one - mixt_frac ) * precip_frac_2 * mu_hm_2**2 ) &
         / ( mixt_frac * precip_frac_1 * ( one + zeta_vrnce_rat ) * mu_hm_1**2 &
             + ( one - mixt_frac ) * precip_frac_2 * mu_hm_2**2 )

       ! Mathematically, this ratio can never be less than 0.  In case numerical
       ! round off error produces a negative value in extreme cases, reset the
       ! value of R to 0.
       if ( sigma_hm_2_sqd_on_mu_hm_2_sqd < zero ) then
          sigma_hm_2_sqd_on_mu_hm_2_sqd = zero
       endif

    elseif ( mu_hm_2 < mu_hm_2_min ) then

       ! Set the value of mu_hm_2 to the threshold positive value.
       mu_hm_2 = mu_hm_2_min

       ! Recalculate the mean (in-precip.) of the hydrometeor in the 1st PDF
       ! component.
       mu_hm_1 = ( hmm - ( one - mixt_frac ) * precip_frac_2 * mu_hm_2 ) &
                 / ( mixt_frac * precip_frac_1 )

       ! Recalculate the value of R ( sigma_hm_2^2 / mu_hm_2^2 ) in this
       ! scenario.
       ! R = ( <hm'^2> + <hm>^2 - a * f_p_1 * mu_hm_1^2
       !       - ( 1 - a ) * f_p_2 * mu_hm_2^2 )
       !     / ( a * f_p_1 * ( 1 + zeta ) * mu_hm_1^2
       !         + ( 1 - a ) * f_p_2 * mu_hm_2^2 ).
       sigma_hm_2_sqd_on_mu_hm_2_sqd &
       = ( hmp2 + hmm**2 - mixt_frac * precip_frac_1 * mu_hm_1**2 &
           - ( one - mixt_frac ) * precip_frac_2 * mu_hm_2**2 ) &
         / ( mixt_frac * precip_frac_1 * ( one + zeta_vrnce_rat ) * mu_hm_1**2 &
             + ( one - mixt_frac ) * precip_frac_2 * mu_hm_2**2 )

       ! Mathematically, this ratio can never be less than 0.  In case numerical
       ! round off error produces a negative value in extreme cases, reset the
       ! value of R to 0.
       if ( sigma_hm_2_sqd_on_mu_hm_2_sqd < zero ) then
          sigma_hm_2_sqd_on_mu_hm_2_sqd = zero
       endif

    endif
 
    ! Calculate the standard deviation (in-precip.) of the hydrometeor in the
    ! 1st PDF component.
    sigma_hm_1 = sqrt( sigma_hm_2_sqd_on_mu_hm_2_sqd &
                       * ( one + zeta_vrnce_rat ) ) &
                 * mu_hm_1

    ! Calculate the standard deviation (in-precip.) of the hydrometeor in the
    ! 2nd PDF component.
    sigma_hm_2 = sqrt( sigma_hm_2_sqd_on_mu_hm_2_sqd ) * mu_hm_2

    ! Calculate the mean of the hydrometeor in the 1st PDF component.
    hm_1 = max( mu_hm_1 * precip_frac_1, hm_tol )

    ! Calculate the mean of the hydrometeor in the 1st PDF component.
    hm_2 = max( mu_hm_2 * precip_frac_2, hm_tol )

    ! Calculate the ratio of sigma_hm_1^2 / mu_hm_1^2.
    sigma_hm_1_sqd_on_mu_hm_1_sqd = sigma_hm_1**2 / mu_hm_1**2

    ! The value of R, sigma_hm_2_sqd_on_mu_hm_2_sqd, has already been
    ! calculated.


    return

  end subroutine calc_mu_sigma_two_comps

  !=============================================================================
  function component_corr_w_x( corr_w_x, rc_i, cloud_frac_i, &
                               corr_w_x_NN_cloud, corr_w_x_NN_below ) &
  result( corr_w_x_i )

    ! Description:
    ! Calculates the correlation of w and x within the ith PDF component.
    ! Here, x is a variable with a normally distributed individual marginal PDF,
    ! such as chi or eta.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one,    & ! Constant(s)
        zero,   &
        rc_tol

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use model_flags, only: &
        l_calc_w_corr

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      corr_w_x,     & ! Correlation of w and x (overall)               [-]
      rc_i,         & ! Mean cloud water mixing ratio (ith PDF comp.)  [kg/kg]
      cloud_frac_i    ! Cloud fraction (ith PDF component)             [-]

    real( kind = core_rknd ), intent(in) :: &
      corr_w_x_NN_cloud, & ! Corr. of w and x (ith PDF comp.); cloudy levs [-]
      corr_w_x_NN_below    ! Corr. of w and x (ith PDF comp.); clear levs  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      corr_w_x_i    ! Correlation of w and x (ith PDF component)  [-]

    ! Local Variables

    ! The component correlations of w and r_t and the component correlations of
    ! w and theta_l are both set to be 0 within the CLUBB model code.  In other
    ! words, w and r_t (theta_l) have overall covariance w'r_t' (w'theta_l'),
    ! but the single component covariance and correlation are defined to be 0.
    ! Since the component covariances (or correlations) of w and chi (old s) and
    ! of w and eta (old t) are based on the covariances (or correlations) of w
    ! and r_t and of w and theta_l, the single component correlation and
    ! covariance of w and chi, as well of as w and eta, are defined to be 0.
    logical, parameter :: &
      l_follow_CLUBB_PDF_standards = .true.


    ! Correlation of w and x in the ith PDF component.
    if ( l_follow_CLUBB_PDF_standards ) then

       ! The component correlations of w and r_t and the component correlations
       ! of w and theta_l are both set to be 0 within the CLUBB model code.  In
       ! other words, w and r_t (theta_l) have overall covariance w'r_t'
       ! (w'theta_l'), but the single component covariance and correlation are
       ! defined to be 0.  Since the component covariances (or correlations)
       ! of w and chi (old s) and of w and eta (old t) are based on the
       ! covariances (or correlations) of w and r_t and of w and theta_l, the
       ! single component correlation and covariance of w and chi, as well as of
       ! w and eta, are defined to be 0.
       corr_w_x_i = zero

    else ! not following CLUBB PDF standards

       ! WARNING:  the standards used in the generation of the two-component
       !           CLUBB PDF are not being obeyed.  The use of this code is
       !           inconsistent with the rest of CLUBB's PDF.
       if ( l_calc_w_corr ) then
          corr_w_x_i = corr_w_x
       else ! use prescribed parameter values
          if ( l_interp_prescribed_params ) then
             corr_w_x_i = cloud_frac_i * corr_w_x_NN_cloud &
                          + ( one - cloud_frac_i ) * corr_w_x_NN_below
          else
             if ( rc_i > rc_tol ) then
                corr_w_x_i = corr_w_x_NN_cloud
             else
                corr_w_x_i = corr_w_x_NN_below
             endif
          endif ! l_interp_prescribed_params
       endif ! l_calc_w_corr

    endif ! l_follow_CLUBB_PDF_standards


    return

  end function component_corr_w_x

  !=============================================================================
  function component_corr_chi_eta( pdf_corr_chi_eta_i, rc_i, cloud_frac_i, &
                                   corr_chi_eta_NN_cloud, &
                                   corr_chi_eta_NN_below, &
                                   l_limit_corr_chi_eta ) &
  result( corr_chi_eta_i )

    ! Description:
    ! Calculates the correlation of chi (old s) and eta (old t) within the
    ! ith PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one,    & ! Constant(s)
        rc_tol, &
        max_mag_correlation

    use model_flags, only: &
        l_fix_chi_eta_correlations  ! Variable(s)

    use clubb_precision, only: &
        core_rknd  ! Constant

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      pdf_corr_chi_eta_i, & ! Correlation of chi and eta (ith PDF component) [-]
      rc_i,               & ! Mean cloud water mix. rat. (ith PDF comp.) [kg/kg]
      cloud_frac_i          ! Cloud fraction (ith PDF component)             [-]

    real( kind = core_rknd ), intent(in) :: &
      corr_chi_eta_NN_cloud, & ! Corr. of chi & eta (ith PDF comp.); cloudy  [-]
      corr_chi_eta_NN_below    ! Corr. of chi & eta (ith PDF comp.); clear   [-]

    logical, intent(in) :: &
      l_limit_corr_chi_eta    ! We must limit the correlation of chi and eta if
                              ! we are to take the Cholesky decomposition of the
                              ! resulting correlation matrix. This is because a
                              ! perfect correlation of chi and eta was found to
                              ! be unrealizable.

    ! Return Variable
    real( kind = core_rknd ) :: &
      corr_chi_eta_i    ! Correlation of chi and eta (ith PDF component)     [-]


    ! Correlation of chi (old s) and eta (old t) in the ith PDF component.

    ! The PDF variables chi and eta result from a transformation of the PDF
    ! involving r_t and theta_l.  The correlation of chi and eta depends on the
    ! correlation of r_t and theta_l, as well as the variances of r_t and
    ! theta_l, and other factors.  The correlation of chi and eta is subject to
    ! change at every vertical level and model time step, and is calculated as
    ! part of the CLUBB PDF parameters.
    if ( .not. l_fix_chi_eta_correlations ) then

       ! Preferred, more accurate version.
       corr_chi_eta_i = pdf_corr_chi_eta_i

    else ! fix the correlation of chi (old s) and eta (old t).

       ! WARNING:  this code is inconsistent with the rest of CLUBB's PDF.  This
       !           code is necessary because SILHS is lazy and wussy, and only
       !           wants to declare correlation arrays at the start of the model
       !           run, rather than updating them throughout the model run.
       if ( l_interp_prescribed_params ) then
          corr_chi_eta_i = cloud_frac_i * corr_chi_eta_NN_cloud &
                           + ( one - cloud_frac_i ) * corr_chi_eta_NN_below
       else
          if ( rc_i > rc_tol ) then
             corr_chi_eta_i = corr_chi_eta_NN_cloud
          else
             corr_chi_eta_i = corr_chi_eta_NN_below
          endif
       endif

    endif

    ! We cannot have a perfect correlation of chi (old s) and eta (old t) if we
    ! plan to decompose this matrix and we don't want the Cholesky_factor code
    ! to throw a fit.
    if ( l_limit_corr_chi_eta ) then

       corr_chi_eta_i = max( min( corr_chi_eta_i, max_mag_correlation ), &
                             -max_mag_correlation )

    endif


    return

  end function component_corr_chi_eta

  !=============================================================================
  function component_corr_w_hm_n_ip( corr_w_hm_i_n_in, rc_i, cloud_frac_i, &
                                     corr_w_hm_n_NL_cloud, &
                                     corr_w_hm_n_NL_below ) &
  result( corr_w_hm_i_n )

    ! Description:
    ! Calculates the in-precip correlation of w and the natural logarithm of a
    ! hydrometeor species within the ith PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one,    & ! Constant(s)
        rc_tol

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use model_flags, only: &
        l_calc_w_corr

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      corr_w_hm_i_n_in, & ! Correlation of w and ln hm (ith PDF comp.) ip    [-]
      rc_i,             & ! Mean cloud water mix. ratio (ith PDF comp.)  [kg/kg]
      cloud_frac_i        ! Cloud fraction (ith PDF component)               [-]

    real( kind = core_rknd ), intent(in) :: &
      corr_w_hm_n_NL_cloud, & ! Corr. of w & ln hm (ith PDF comp.) ip; cloud [-]
      corr_w_hm_n_NL_below    ! Corr. of w & ln hm (ith PDF comp.) ip; clear [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      corr_w_hm_i_n    ! Correlation of w and ln hm (ith PDF component) ip   [-]


    ! Correlation (in-precip) of w and the natural logarithm of the hydrometeor
    ! in the ith PDF component.
    if ( l_calc_w_corr ) then
       corr_w_hm_i_n = corr_w_hm_i_n_in
    else ! use prescribed parameter values
       if ( l_interp_prescribed_params ) then
          corr_w_hm_i_n = cloud_frac_i * corr_w_hm_n_NL_cloud &
                          + ( one - cloud_frac_i ) * corr_w_hm_n_NL_below
       else
          if ( rc_i > rc_tol ) then
             corr_w_hm_i_n = corr_w_hm_n_NL_cloud
          else
             corr_w_hm_i_n = corr_w_hm_n_NL_below
          endif
       endif ! l_interp_prescribed_params
    endif ! l_calc_w_corr

    return

  end function component_corr_w_hm_n_ip

  !=============================================================================
  function component_corr_x_hm_n_ip( rc_i, cloud_frac_i, &
                                     corr_x_hm_n_NL_cloud, &
                                     corr_x_hm_n_NL_below ) &
  result( corr_x_hm_i_n )

    ! Description:
    ! Calculates the in-precip correlation of x and a hydrometeor species
    ! within the ith PDF component.  Here, x is a variable with a normally
    ! distributed individual marginal PDF, such as chi or eta.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one,    & ! Constant(s)
        rc_tol

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      rc_i,         & ! Mean cloud water mixing ratio (ith PDF comp.)   [kg/kg]
      cloud_frac_i    ! Cloud fraction (ith PDF component)              [-]

    real( kind = core_rknd ), intent(in) :: &
      corr_x_hm_n_NL_cloud, & ! Corr. of x and ln hm (ith PDF comp.) ip     [-]
      corr_x_hm_n_NL_below    ! Corr. of x and ln hm (ith PDF comp.) ip     [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      corr_x_hm_i_n    ! Correlation of x and ln hm (ith PDF component) ip  [-]


    ! Correlation (in-precip) of x and the hydrometeor in the ith PDF component.
    if ( l_interp_prescribed_params ) then
       corr_x_hm_i_n = cloud_frac_i * corr_x_hm_n_NL_cloud &
                       + ( one - cloud_frac_i ) * corr_x_hm_n_NL_below
    else
       if ( rc_i > rc_tol ) then
          corr_x_hm_i_n = corr_x_hm_n_NL_cloud
       else
          corr_x_hm_i_n = corr_x_hm_n_NL_below
       endif
    endif


    return

  end function component_corr_x_hm_n_ip

  !=============================================================================
  function component_corr_hmx_hmy_n_ip( rc_i, cloud_frac_i, &
                                        corr_hmx_hmy_n_LL_cloud, &
                                        corr_hmx_hmy_n_LL_below ) &
  result( corr_hmx_hmy_i_n )

    ! Description:
    ! Calculates the in-precip correlation of the natural logarithms of
    ! hydrometeor x and hydrometeor y within the ith PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one,    & ! Constant(s)
        rc_tol

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      rc_i,         & ! Mean cloud water mixing ratio (ith PDF comp.) [kg/kg]
      cloud_frac_i    ! Cloud fraction (ith PDF component)            [-]

    real( kind = core_rknd ), intent(in) :: &
      corr_hmx_hmy_n_LL_cloud, & ! Corr.: ln hmx & ln hmy (ith PDF comp.) ip [-]
      corr_hmx_hmy_n_LL_below    ! Corr.: ln hmx & ln hmy (ith PDF comp.) ip [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      corr_hmx_hmy_i_n    ! Corr. of ln hmx & ln hmy (ith PDF comp.) ip      [-]


    ! Correlation (in-precip) of the natural logarithms of hydrometeor x and
    ! hydrometeor y in the ith PDF component.
    if ( l_interp_prescribed_params ) then
       corr_hmx_hmy_i_n = cloud_frac_i * corr_hmx_hmy_n_LL_cloud &
                        + ( one - cloud_frac_i ) * corr_hmx_hmy_n_LL_below
    else
       if ( rc_i > rc_tol ) then
          corr_hmx_hmy_i_n = corr_hmx_hmy_n_LL_cloud
       else
          corr_hmx_hmy_i_n = corr_hmx_hmy_n_LL_below
       endif
    endif


    return

  end function component_corr_hmx_hmy_n_ip

  !=============================================================================
  pure function component_corr_eta_hm_n_ip( corr_chi_eta_i, corr_chi_hm_n_i ) &
  result( corr_eta_hm_n_i )

    ! Description:
    ! Estimates the correlation of eta and the natural logarithm of a
    ! hydrometeor species using the correlation of chi and eta and the
    ! correlation of chi and the natural logarithm of the hydrometeor.  This
    ! facilities the Cholesky decomposability of the correlation array that will
    ! inevitably be decomposed for SILHS purposes. Without this estimation, we
    ! have found that the resulting correlation matrix cannot be decomposed.

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd       ! Constant

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      corr_chi_eta_i,  & ! Component correlation of chi and eta              [-]
      corr_chi_hm_n_i    ! Component correlation of chi and ln hm            [-]

    ! Output Variables
    real( kind = core_rknd ) :: &
      corr_eta_hm_n_i    ! Component correlation of eta and ln hm            [-]


    corr_eta_hm_n_i = corr_chi_eta_i * corr_chi_hm_n_i


    return

  end function component_corr_eta_hm_n_ip

  !=============================================================================
  subroutine norm_transform_mean_stdev( hm_1, hm_2, &
                                        Ncnm, d_variables, &
                                        mu_x_1, mu_x_2, &
                                        sigma_x_1, sigma_x_2, &
                                        sigma2_on_mu2_ip_1, &
                                        sigma2_on_mu2_ip_2, &
                                        mu_x_1_n, mu_x_2_n, &
                                        sigma_x_1_n, sigma_x_2_n )

    ! Description:
    ! Transforms the means and the standard deviations of PDF variables that
    ! have assumed lognormal distributions -- which are precipitating
    ! hydrometeors (in precipitation) and N_cn -- to normal space for each PDF
    ! component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        Ncn_tol, &  ! Constant(s)
        zero

    use pdf_utilities, only: &
        mean_L2N,  & ! Procedure(s)
        stdev_L2N

    use index_mapping, only: &
        pdf2hydromet_idx  ! Procedure(s)

    use corr_varnce_module, only: &
        iiPDF_Ncn  ! Variable(s)

    use array_index, only: &
        hydromet_tol  ! Variable(s)

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use model_flags, only: &
        l_const_Nc_in_cloud ! Variable

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(hydromet_dim), intent(in) :: &
      hm_1, & ! Mean of a precip. hydrometeor (1st PDF component)   [units vary]
      hm_2    ! Mean of a precip. hydrometeor (2nd PDF component)   [units vary]

    real( kind = core_rknd ), intent(in) :: &
      Ncnm    ! Mean cloud nuclei concentration, < N_cn >               [num/kg]

    integer, intent(in) :: &
      d_variables    ! Number of variables in CLUBB's PDF

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      mu_x_1,    & ! Mean array of PDF vars. (1st PDF component)    [units vary]
      mu_x_2,    & ! Mean array of PDF vars. (2nd PDF component)    [units vary]
      sigma_x_1, & ! Standard deviation array of PDF vars (comp. 1) [units vary]
      sigma_x_2    ! Standard deviation array of PDF vars (comp. 2) [units vary]

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      sigma2_on_mu2_ip_1, & ! Prescribed ratio array: sigma_hm_1^2/mu_hm_1^2 [-]
      sigma2_on_mu2_ip_2    ! Prescribed ratio array: sigma_hm_2^2/mu_hm_2^2 [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(d_variables), intent(out) :: &
      mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    ! Local Variable
    integer :: ivar  ! Loop index


    ! The means and standard deviations in each PDF component of w, chi (old s),
    ! and eta (old t) do not need to be transformed to normal space, since w,
    ! chi, and eta already follow assumed normal distributions in each PDF
    ! component.  The normal space means and standard deviations are the same as
    ! the actual means and standard deviations.    
    mu_x_1_n = mu_x_1
    mu_x_2_n = mu_x_2
    sigma_x_1_n = sigma_x_1
    sigma_x_2_n = sigma_x_2

    !!! Transform the mean and standard deviation to normal space in each PDF
    !!! component for variables that have an assumed lognormal distribution,
    !!! given the mean and standard deviation in each PDF component for those
    !!! variables.  A precipitating hydrometeor has an assumed lognormal
    !!! distribution in precipitation in each PDF component.  Simplified cloud
    !!! nuclei concentration, N_cn, has an assumed lognormal distribution in
    !!! each PDF component, and furthermore, mu_Ncn_1 = mu_Ncn_2 and
    !!! sigma_Ncn_1 = sigma_Ncn_2, so N_cn has an assumed single lognormal
    !!! distribution over the entire domain.

    ! Normal space mean of simplified cloud nuclei concentration, N_cn,
    ! in PDF component 1.
    if ( Ncnm >= Ncn_tol ) then

       mu_x_1_n(iiPDF_Ncn) = mean_L2N( mu_x_1(iiPDF_Ncn), &
                                       sigma2_on_mu2_ip_1(iiPDF_Ncn) )

    else

       ! Mean simplified cloud nuclei concentration in PDF component 1 is less
       ! than the tolerance amount.  It is considered to have a value of 0.
       ! There are not any cloud nuclei or cloud at this grid level.  The value
       ! of mu_Ncn_1_n should be -inf.  It will be set to -huge for purposes of
       ! assigning it a value.
       mu_x_1_n(iiPDF_Ncn) = -huge( mu_x_1(iiPDF_Ncn) )

    endif

    ! Normal space mean of simplified cloud nuclei concentration, N_cn,
    ! in PDF component 2.
    if ( Ncnm >= Ncn_tol ) then

       mu_x_2_n(iiPDF_Ncn) = mean_L2N( mu_x_2(iiPDF_Ncn), &
                                       sigma2_on_mu2_ip_1(iiPDF_Ncn) )

    else

       ! Mean simplified cloud nuclei concentration in PDF component 1 is less
       ! than the tolerance amount.  It is considered to have a value of 0.
       ! There are not any cloud nuclei or cloud at this grid level.  The value
       ! of mu_Ncn_1_n should be -inf.  It will be set to -huge for purposes of
       ! assigning it a value.
       mu_x_2_n(iiPDF_Ncn) = -huge( mu_x_2(iiPDF_Ncn) )

    endif

    ! Normal space standard deviation of simplified cloud nuclei concentration,
    ! N_cn, in PDF components 1 and 2.
    if ( l_const_Nc_in_cloud ) then
      ! Ncn does not vary in the grid box.
      sigma_x_1_n(iiPDF_Ncn) = zero
      sigma_x_2_n(iiPDF_Ncn) = zero
    else
      ! Ncn (perhaps) varies in the grid box.
      sigma_x_1_n(iiPDF_Ncn) = stdev_L2N( sigma2_on_mu2_ip_1(iiPDF_Ncn) )
      sigma_x_2_n(iiPDF_Ncn) = stdev_L2N( sigma2_on_mu2_ip_2(iiPDF_Ncn) )
    end if

    ! Normal space precipitating hydrometeor means and standard deviations.
    do ivar = iiPDF_Ncn+1, d_variables, 1

       ! Normal space mean of a precipitating hydrometeor, hm, in PDF
       ! component 1.
       if ( hm_1(pdf2hydromet_idx(ivar)) &
            >= hydromet_tol(pdf2hydromet_idx(ivar)) ) then

          mu_x_1_n(ivar) = mean_L2N( mu_x_1(ivar), sigma2_on_mu2_ip_1(ivar) )

       else

          ! The mean of a precipitating hydrometeor in PDF component 1 is less
          ! than its tolerance amount.  It is considered to have a value of 0.
          ! There is not any of this precipitating hydrometeor in the 1st PDF
          ! component at this grid level.  The in-precip mean of this
          ! precipitating hydrometeor (1st PDF component) is also 0.  The value
          ! of mu_hm_1_n should be -inf.  It will be set to -huge for purposes
          ! of assigning it a value.
          mu_x_1_n(ivar) = -huge( mu_x_1(ivar) )

       endif

       ! Normal space standard deviation of a precipitating hydrometeor, hm, in
       ! PDF component 1.
       sigma_x_1_n(ivar) = stdev_L2N( sigma2_on_mu2_ip_1(ivar) )

       ! Normal space mean of a precipitating hydrometeor, hm, in PDF
       ! component 2.
       if ( hm_2(pdf2hydromet_idx(ivar)) &
            >= hydromet_tol(pdf2hydromet_idx(ivar)) ) then

          mu_x_2_n(ivar) = mean_L2N( mu_x_2(ivar), sigma2_on_mu2_ip_2(ivar) )

       else

          ! The mean of a precipitating hydrometeor in PDF component 2 is less
          ! than its tolerance amount.  It is considered to have a value of 0.
          ! There is not any of this precipitating hydrometeor in the 2nd PDF
          ! component at this grid level.  The in-precip mean of this
          ! precipitating hydrometeor (2nd PDF component) is also 0.  The value
          ! of mu_hm_2_n should be -inf.  It will be set to -huge for purposes
          ! of assigning it a value.
          mu_x_2_n(ivar) = -huge( mu_x_2(ivar) )

       endif

       ! Normal space standard deviation of a precipitating hydrometeor, hm, in
       ! PDF component 2.
       sigma_x_2_n(ivar) = stdev_L2N( sigma2_on_mu2_ip_2(ivar) )

    enddo ! ivar = iiPDF_Ncn+1, d_variables, 1


    return

  end subroutine norm_transform_mean_stdev

  !=============================================================================
  subroutine denorm_transform_corr( d_variables, &
                                    sigma_x_1_n, sigma_x_2_n, &
                                    sigma2_on_mu2_ip_1, sigma2_on_mu2_ip_2, &
                                    corr_array_1_n, &
                                    corr_array_2_n, &
                                    corr_array_1, corr_array_2 )

    ! Description:
    ! Calculates the true or "real-space" correlations between PDF variables,
    ! where at least one of the variables that is part of a correlation has an
    ! assumed lognormal distribution -- which are the precipitating hydrometeors
    ! (in precipitation) and N_cn.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        zero          ! Constant

    use pdf_utilities, only: &
        corr_NN2NL, & ! Procedure(s)
        corr_NN2LL

    use corr_varnce_module, only: &
        iiPDF_chi, & ! Variable(s)
        iiPDF_eta, &
        iiPDF_w,   &
        iiPDF_Ncn

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      d_variables ! Number of PDF variables

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real ( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      sigma2_on_mu2_ip_1, & ! Ratio array sigma_hm_1^2/mu_hm_1^2             [-]
      sigma2_on_mu2_ip_2    ! Ratio array sigma_hm_2^2/mu_hm_2^2             [-]

    real( kind = core_rknd ), dimension(d_variables, d_variables), &
    intent(in) :: &
      corr_array_1_n, & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
      corr_array_2_n    ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(d_variables, d_variables), &
    intent(out) :: &
      corr_array_1, & ! Correlation array of PDF vars. (comp. 1)             [-]
      corr_array_2    ! Correlation array of PDF vars. (comp. 2)             [-]

    ! Local Variables
    integer :: ivar, jvar ! Loop indices


    ! The correlations in each PDF component between two of w, chi (old s), and
    ! eta (old t) do not need to be transformed to standard space, since w, chi,
    ! and eta follow assumed normal distributions in each PDF component.  The
    ! normal space correlations between any two of these variables are the same
    ! as the actual correlations.    
    corr_array_1 = corr_array_1_n
    corr_array_2 = corr_array_2_n

    !!! Calculate the true correlation of variables that have an assumed normal
    !!! distribution and variables that have an assumed lognormal distribution
    !!! for the ith PDF component, given their normal space correlation and the
    !!! normal space standard deviation of the variable with the assumed
    !!! lognormal distribution.

    ! Transform the correlations between chi/eta/w and N_cn to standard space.

    ! Transform the correlation of w and N_cn to standard space in PDF
    ! component 1.
    corr_array_1(iiPDF_Ncn, iiPDF_w) &
    = corr_NN2NL( corr_array_1_n(iiPDF_Ncn, iiPDF_w), &
                  sigma_x_1_n(iiPDF_Ncn), sigma2_on_mu2_ip_1(iiPDF_Ncn) )

    ! Transform the correlation of w and N_cn to standard space in PDF
    ! component 2.
    corr_array_2(iiPDF_Ncn, iiPDF_w) &
    = corr_NN2NL( corr_array_2_n(iiPDF_Ncn, iiPDF_w), &
                  sigma_x_2_n(iiPDF_Ncn), sigma2_on_mu2_ip_1(iiPDF_Ncn) )

    ! Transform the correlation of chi (old s) and N_cn to standard space in
    ! PDF component 1.
    corr_array_1(iiPDF_Ncn, iiPDF_chi) &
    = corr_NN2NL( corr_array_1_n(iiPDF_Ncn, iiPDF_chi), &
                  sigma_x_1_n(iiPDF_Ncn), sigma2_on_mu2_ip_1(iiPDF_Ncn) )

    ! Transform the correlation of chi (old s) and N_cn to standard space in
    ! PDF component 2.
    corr_array_2(iiPDF_Ncn, iiPDF_chi) &
    = corr_NN2NL( corr_array_2_n(iiPDF_Ncn, iiPDF_chi), &
                  sigma_x_2_n(iiPDF_Ncn), sigma2_on_mu2_ip_1(iiPDF_Ncn) )

    ! Transform the correlation of eta (old t) and N_cn to standard space in
    ! PDF component 1.
    corr_array_1(iiPDF_Ncn, iiPDF_eta) &
    = corr_NN2NL( corr_array_1_n(iiPDF_Ncn, iiPDF_eta), &
                  sigma_x_1_n(iiPDF_Ncn), sigma2_on_mu2_ip_1(iiPDF_Ncn) )

    ! Transform the correlation of eta (old t) and N_cn to standard space in
    ! PDF component 2.
    corr_array_2(iiPDF_Ncn, iiPDF_eta) &
    = corr_NN2NL( corr_array_2_n(iiPDF_Ncn, iiPDF_eta), &
                  sigma_x_2_n(iiPDF_Ncn), sigma2_on_mu2_ip_1(iiPDF_Ncn) )

    ! Transform the correlations (in-precip) between chi/eta/w and the
    ! precipitating hydrometeors to standard space.
    do ivar = iiPDF_chi, iiPDF_w
       do jvar = iiPDF_Ncn+1, d_variables

          ! Transform the correlation (in-precip) between w, chi, or eta and a
          ! precipitating hydrometeor, hm, to standard space in PDF component 1.
          corr_array_1(jvar, ivar) &
          = corr_NN2NL( corr_array_1_n(jvar, ivar), sigma_x_1_n(jvar), &
                        sigma2_on_mu2_ip_1(jvar) )

          ! Transform the correlation (in-precip) between w, chi, or eta and a
          ! precipitating hydrometeor, hm, to standard space in PDF component 2.
          corr_array_2(jvar, ivar) &
          = corr_NN2NL( corr_array_2_n(jvar, ivar), sigma_x_2_n(jvar), &
                        sigma2_on_mu2_ip_2(jvar) )

       enddo ! jvar = iiPDF_Ncn+1, d_variables
    enddo ! ivar = iiPDF_chi, iiPDF_w


    !!! Calculate the true correlation of two variables that both have an
    !!! assumed lognormal distribution for the ith PDF component, given their
    !!! normal space correlation and both of their normal space standard
    !!! deviations.

    ! Transform the correlations (in-precip) between N_cn and the precipitating
    ! hydrometeors to standard space.
    ivar = iiPDF_Ncn
    do jvar = ivar+1, d_variables

       ! Transform the correlation (in-precip) between N_cn and a precipitating
       ! hydrometeor, hm, to standard space in PDF component 1.
       corr_array_1(jvar, ivar) &
       = corr_NN2LL( corr_array_1_n(jvar, ivar), &
                     sigma_x_1_n(ivar), sigma_x_1_n(jvar), &
                     sigma2_on_mu2_ip_1(iiPDF_Ncn), sigma2_on_mu2_ip_1(jvar) )

       ! Transform the correlation (in-precip) between N_cn and a precipitating
       ! hydrometeor, hm, to standard space in PDF component 2.
       corr_array_2(jvar, ivar) &
       = corr_NN2LL( corr_array_2_n(jvar, ivar), &
                     sigma_x_2_n(ivar), sigma_x_2_n(jvar), &
                     sigma2_on_mu2_ip_1(iiPDF_Ncn), sigma2_on_mu2_ip_2(jvar) )

    enddo ! jvar = ivar+1, d_variables

    ! Transform the correlations (in-precip) between two precipitating
    ! hydrometeors to standard space.
    do ivar = iiPDF_Ncn+1, d_variables-1
       do jvar = ivar+1, d_variables

          ! Transform the correlation (in-precip) between two precipitating
          ! hydrometeors (for example, r_r and N_r) to standard space in PDF
          ! component 1.
          corr_array_1(jvar, ivar) &
          = corr_NN2LL( corr_array_1_n(jvar, ivar), &
                        sigma_x_1_n(ivar), sigma_x_1_n(jvar), &
                        sigma2_on_mu2_ip_1(ivar), sigma2_on_mu2_ip_1(jvar) )

          ! Transform the correlation (in-precip) between two precipitating
          ! hydrometeors (for example, r_r and N_r) to standard space in PDF
          ! component 2.
          corr_array_2(jvar, ivar) &
          = corr_NN2LL( corr_array_2_n(jvar, ivar), &
                        sigma_x_2_n(ivar), sigma_x_2_n(jvar), &
                        sigma2_on_mu2_ip_2(ivar), sigma2_on_mu2_ip_2(jvar) )

       enddo ! jvar = ivar+1, d_variables
    enddo ! ivar = iiPDF_Ncn+1, d_variables-1


    return

  end subroutine denorm_transform_corr

  !=============================================================================
  subroutine calc_corr_w_hm_n( wm, wphydrometp, &
                               mu_w_1, mu_w_2, &
                               mu_hm_1, mu_hm_2, &
                               sigma_w_1, sigma_w_2, &
                               sigma_hm_1, sigma_hm_2, &
                               sigma_hm_1_n, sigma_hm_2_n, &
                               mixt_frac, precip_frac_1, precip_frac_2, &
                               corr_w_hm_1_n, corr_w_hm_2_n, &
                               hm_tol )

    ! Description:
    ! Calculates the PDF component correlation (in-precip) between vertical
    ! velocity, w, and the natural logarithm of a hydrometeor, ln hm.  The
    ! overall covariance of w and hm, <w'hm'> can be written in terms of the PDF
    ! parameters.  When both w and hm vary in both PDF components, the equation
    ! is written as:
    !
    ! <w'hm'> = mixt_frac * precip_frac_1
    !           * ( mu_w_1 - <w>
    !               + corr_w_hm_1_n * sigma_w_1 * sigma_hm_1_n ) * mu_hm_1
    !           + ( 1 - mixt_frac ) * precip_frac_2
    !             * ( mu_w_2 - <w>
    !                 + corr_w_hm_2_n * sigma_w_2 * sigma_hm_2_n ) * mu_hm_2.
    !
    ! The overall covariance is provided, so the component correlation is solved
    ! by setting corr_w_hm_1_n = corr_w_hm_2_n ( = corr_w_hm_n ).  The equation
    ! is:
    !
    ! corr_w_hm_n
    ! = ( <w'hm'>
    !     - mixt_frac * precip_frac_1 * ( mu_w_1 - <w> ) * mu_hm_1
    !     - ( 1 - mixt_frac ) * precip_frac_2 * ( mu_w_2 - <w> ) * mu_hm_2 )
    !   / ( mixt_frac * precip_frac_1 * sigma_w_1 * sigma_hm_1_n * mu_hm_1
    !       + ( 1 - mixt_frac ) * precip_frac_2
    !         * sigma_w_2 * sigma_hm_2_n * mu_hm_2 );
    !
    ! again, where corr_w_hm_1_n = corr_w_hm_2_n = corr_w_hm_n.  When either w
    ! or hm is constant in one PDF component, but both w and hm vary in the
    ! other PDF component, the equation for <w'hm'> is written as:
    !
    ! <w'hm'> = mixt_frac * precip_frac_1
    !           * ( mu_w_1 - <w>
    !               + corr_w_hm_1_n * sigma_w_1 * sigma_hm_1_n ) * mu_hm_1
    !           + ( 1 - mixt_frac ) * precip_frac_2
    !             * ( mu_w_2 - <w> ) * mu_hm_2.
    !
    ! In the above equation, either w or hm (or both) is (are) constant in PDF
    ! component 2, but both w and hm vary in PDF component 1.  When both w and
    ! hm vary in PDF component 2, but at least one of w or hm is constant in PDF
    ! component 1, the equation is similar.  The above equation can be rewritten
    ! to solve for corr_w_hm_1_n, such that:
    !
    ! corr_w_hm_1_n
    ! = ( <w'hm'>
    !     - mixt_frac * precip_frac_1 * ( mu_w_1 - <w> ) * mu_hm_1
    !     - ( 1 - mixt_frac ) * precip_frac_2 * ( mu_w_2 - <w> ) * mu_hm_2 )
    !   / ( mixt_frac * precip_frac_1 * sigma_w_1 * sigma_hm_1_n * mu_hm_1 ).
    !
    ! Since either w or hm is constant in PDF component 2, corr_w_hm_2_n is
    ! undefined.  When both w and hm vary in PDF component 2, but at least one
    ! of w or hm is constant in PDF component 1, the equation is similar, but
    ! is in terms of corr_w_hm_2_n, while corr_w_hm_1_n is undefined.  When
    ! either w or hm is constant in both PDF components, the equation for
    ! <w'hm'> is:
    !
    ! <w'hm'> = mixt_frac * precip_frac_1
    !           * ( mu_w_1 - <w> ) * mu_hm_1
    !           + ( 1 - mixt_frac ) * precip_frac_2
    !             * ( mu_w_2 - <w> ) * mu_hm_2.
    !
    ! When this is the case, both corr_w_hm_1_n and corr_w_hm_2_n are undefined.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one,                 & ! Constant(s)
        zero,                &
        max_mag_correlation, &
        w_tol

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      wm,            & ! Mean vertical velocity (overall), <w>             [m/s]
      wphydrometp,   & ! Covariance of w and hm (overall), <w'hm'>  [m/s(hm un)]
      mu_w_1,        & ! Mean of w (1st PDF component)                     [m/s]
      mu_w_2,        & ! Mean of w (2nd PDF component)                     [m/s]
      mu_hm_1,       & ! Mean of hm (1st PDF component) in-precip (ip)   [hm un]
      mu_hm_2,       & ! Mean of hm (2nd PDF component) ip               [hm un]
      sigma_w_1,     & ! Standard deviation of w (1st PDF component)       [m/s]
      sigma_w_2,     & ! Standard deviation of w (2nd PDF component)       [m/s]
      sigma_hm_1,    & ! Standard deviation of hm (1st PDF component) ip [hm un]
      sigma_hm_2,    & ! Standard deviation of hm (2nd PDF component) ip [hm un]
      sigma_hm_1_n,  & ! Standard deviation of ln hm (1st PDF component) ip  [-]
      sigma_hm_2_n,  & ! Standard deviation of ln hm (2nd PDF component) ip  [-]
      mixt_frac,     & ! Mixture fraction                                    [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2, & ! Precipitation fraction (2nd PDF component)          [-]
      hm_tol           ! Hydrometeor tolerance value                     [hm un]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      corr_w_hm_1_n, & ! Correlation of w and ln hm (1st PDF component) ip   [-]
      corr_w_hm_2_n    ! Correlation of w and ln hm (2nd PDF component) ip   [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      corr_w_hm_n    ! Correlation of w and ln hm (both PDF components) ip   [-]


    ! Calculate the PDF component correlation of vertical velocity, w, and the
    ! natural logarithm of a hydrometeor, ln hm, in precipitation.
    if ( sigma_w_1 > w_tol .and. sigma_hm_1 > hm_tol .and. &
         sigma_w_2 > w_tol .and. sigma_hm_2 > hm_tol ) then

       ! Both w and hm vary in both PDF components.
       ! Calculate corr_w_hm_n (where corr_w_hm_1_n = corr_w_hm_2_n
       ! = corr_w_hm_n).
       corr_w_hm_n &
       = ( wphydrometp &
           - mixt_frac * precip_frac_1 * ( mu_w_1 - wm ) * mu_hm_1 &
           - ( one - mixt_frac ) * precip_frac_2 * ( mu_w_2 - wm ) * mu_hm_2 ) &
         / ( mixt_frac * precip_frac_1 * sigma_w_1 * sigma_hm_1_n * mu_hm_1 &
             + ( one - mixt_frac ) * precip_frac_2 &
               * sigma_w_2 * sigma_hm_2_n * mu_hm_2 )

       ! Check that the PDF component correlations have reasonable values.
       if ( corr_w_hm_n > max_mag_correlation ) then
          corr_w_hm_n = max_mag_correlation
       elseif ( corr_w_hm_n < -max_mag_correlation ) then
          corr_w_hm_n = -max_mag_correlation
       endif

       ! The PDF component correlations between w and ln hm (in-precip) are
       ! equal.
       corr_w_hm_1_n = corr_w_hm_n
       corr_w_hm_2_n = corr_w_hm_n


    elseif ( sigma_w_1 > w_tol .and. sigma_hm_1 > hm_tol ) then

       ! Both w and hm vary in PDF component 1, but at least one of w and hm is
       ! constant in PDF component 2.
       ! Calculate the PDF component 1 correlation of w and ln hm (in-precip).
       corr_w_hm_1_n &
       = ( wphydrometp &
           - mixt_frac * precip_frac_1 * ( mu_w_1 - wm ) * mu_hm_1 &
           - ( one - mixt_frac ) * precip_frac_2 * ( mu_w_2 - wm ) * mu_hm_2 ) &
         / ( mixt_frac * precip_frac_1 * sigma_w_1 * sigma_hm_1_n * mu_hm_1 )

       ! Check that the PDF component 1 correlation has a reasonable value.
       if ( corr_w_hm_1_n > max_mag_correlation ) then
          corr_w_hm_1_n = max_mag_correlation
       elseif ( corr_w_hm_1_n < -max_mag_correlation ) then
          corr_w_hm_1_n = -max_mag_correlation
       endif

       ! The PDF component 2 correlation is undefined.
       corr_w_hm_2_n = zero
       

    elseif ( sigma_w_2 > w_tol .and. sigma_hm_2 > hm_tol ) then

       ! Both w and hm vary in PDF component 2, but at least one of w and hm is
       ! constant in PDF component 1.
       ! Calculate the PDF component 2 correlation of w and ln hm (in-precip).
       corr_w_hm_2_n &
       = ( wphydrometp &
           - mixt_frac * precip_frac_1 * ( mu_w_1 - wm ) * mu_hm_1 &
           - ( one - mixt_frac ) * precip_frac_2 * ( mu_w_2 - wm ) * mu_hm_2 ) &
         / ( ( one - mixt_frac ) * precip_frac_2 &
             * sigma_w_2 * sigma_hm_2_n * mu_hm_2 )

       ! Check that the PDF component 2 correlation has a reasonable value.
       if ( corr_w_hm_2_n > max_mag_correlation ) then
          corr_w_hm_2_n = max_mag_correlation
       elseif ( corr_w_hm_2_n < -max_mag_correlation ) then
          corr_w_hm_2_n = -max_mag_correlation
       endif

       ! The PDF component 1 correlation is undefined.
       corr_w_hm_1_n = zero
       

    else    ! sigma_w_1 * sigma_hm_1 = 0 .and. sigma_w_2 * sigma_hm_2 = 0.

       ! At least one of w and hm is constant in both PDF components.

       ! The PDF component 1 and component 2 correlations are both undefined.
       corr_w_hm_1_n = zero
       corr_w_hm_2_n = zero


    endif


    return

  end subroutine calc_corr_w_hm_n

  !=============================================================================
  subroutine pdf_param_hm_stats( d_variables, level, hm_1, hm_2, &
                                 mu_x_1, mu_x_2, &
                                 sigma_x_1, sigma_x_2, &
                                 corr_array_1, corr_array_2, &
                                 l_stats_samp )

    ! Description:
    ! Record statistics for standard PDF parameters involving hydrometeors.

    ! References:
    !-----------------------------------------------------------------------

    use index_mapping, only: &
        pdf2hydromet_idx  ! Procedure(s)

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use corr_varnce_module, only: &
        iiPDF_w,   & ! Variable(s)
        iiPDF_chi, &
        iiPDF_eta, &
        iiPDF_Ncn

    use clubb_precision, only: &
        core_rknd   ! Variable(s)

    use stats_type_utilities, only: &
        stat_update_var_pt  ! Procedure(s)

    use stats_variables, only : &
        ihm_1,        & ! Variable(s)
        ihm_2,        &
        imu_hm_1,     &
        imu_hm_2,     &
        imu_Ncn_1,    &
        imu_Ncn_2,    &
        isigma_hm_1,  &
        isigma_hm_2,  &
        isigma_Ncn_1, &
        isigma_Ncn_2

    use stats_variables, only : &
        icorr_w_chi_1,      & ! Variable(s)
        icorr_w_chi_2,      &
        icorr_w_eta_1,      &
        icorr_w_eta_2,      &
        icorr_w_hm_1,       &
        icorr_w_hm_2,       &
        icorr_w_Ncn_1,      &
        icorr_w_Ncn_2,      &
        icorr_chi_eta_1_ca, &
        icorr_chi_eta_2_ca, &
        icorr_chi_hm_1,     &
        icorr_chi_hm_2,     &
        icorr_chi_Ncn_1,    &
        icorr_chi_Ncn_2,    &
        icorr_eta_hm_1,     &
        icorr_eta_hm_2,     &
        icorr_eta_Ncn_1,    &
        icorr_eta_Ncn_2,    &
        icorr_Ncn_hm_1,     &
        icorr_Ncn_hm_2,     &
        icorr_hmx_hmy_1,    &
        icorr_hmx_hmy_2,    &
        stats_zt

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      d_variables, & ! Number of variables in the correlation array
      level          ! Vertical level index 

    real( kind = core_rknd ), dimension(hydromet_dim), intent(in) :: &
      hm_1, & ! Mean of a precip. hydrometeor (1st PDF component)   [units vary]
      hm_2    ! Mean of a precip. hydrometeor (2nd PDF component)   [units vary]

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      mu_x_1,    & ! Mean array of PDF vars. (1st PDF component)    [units vary]
      mu_x_2,    & ! Mean array of PDF vars. (2nd PDF component)    [units vary]
      sigma_x_1, & ! Standard deviation array of PDF vars (comp. 1) [units vary]
      sigma_x_2    ! Standard deviation array of PDF vars (comp. 2) [units vary]

    real( kind = core_rknd ), dimension(d_variables, d_variables), &
    intent(in) :: &
      corr_array_1, & ! Correlation array of PDF vars. (comp. 1)             [-]
      corr_array_2    ! Correlation array of PDF vars. (comp. 2)             [-]

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.

    ! Local Variable
    integer :: ivar, jvar  ! Loop indices


    !!! Output the statistics for hydrometeor PDF parameters.

    ! Statistics
    if ( l_stats_samp ) then

       do ivar = 1, hydromet_dim, 1

          if ( ihm_1(ivar) > 0 ) then
             ! Mean of the precipitating hydrometeor in PDF component 1.
             call stat_update_var_pt( ihm_1(ivar), level, hm_1(ivar), stats_zt )
          endif

          if ( ihm_2(ivar) > 0 ) then
             ! Mean of the precipitating hydrometeor in PDF component 2.
             call stat_update_var_pt( ihm_2(ivar), level, hm_2(ivar), stats_zt )
          endif

       enddo ! ivar = 1, hydromet_dim, 1

       do ivar = iiPDF_Ncn+1, d_variables, 1

          ! Mean of the precipitating hydrometeor (in-precip)
          ! in PDF component 1.
          if ( imu_hm_1(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt( imu_hm_1(pdf2hydromet_idx(ivar)), &
                                      level, mu_x_1(ivar), stats_zt )
          endif

          ! Mean of the precipitating hydrometeor (in-precip)
          ! in PDF component 2.
          if ( imu_hm_2(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt( imu_hm_2(pdf2hydromet_idx(ivar)), &
                                      level, mu_x_2(ivar), stats_zt )
          endif

       enddo ! ivar = iiPDF_Ncn+1, d_variables, 1

       ! Mean of cloud nuclei concentration in PDF component 1.
       if ( imu_Ncn_1 > 0 ) then
          call stat_update_var_pt( imu_Ncn_1, level, mu_x_1(iiPDF_Ncn), &
                                   stats_zt )
       endif

       ! Mean of cloud nuclei concentration in PDF component 2.
       if ( imu_Ncn_2 > 0 ) then
          call stat_update_var_pt( imu_Ncn_2, level, mu_x_2(iiPDF_Ncn), &
                                   stats_zt )
       endif

       do ivar = iiPDF_Ncn+1, d_variables, 1

          ! Standard deviation of the precipitating hydrometeor (in-precip)
          ! in PDF component 1.
          if ( isigma_hm_1(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt( isigma_hm_1(pdf2hydromet_idx(ivar)), &
                                      level, sigma_x_1(ivar), stats_zt )
          endif

          ! Standard deviation of the precipitating hydrometeor (in-precip)
          ! in PDF component 2.
          if ( isigma_hm_2(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt( isigma_hm_2(pdf2hydromet_idx(ivar)), &
                                      level, sigma_x_2(ivar), stats_zt )
          endif

       enddo ! ivar = iiPDF_Ncn+1, d_variables, 1

       ! Standard deviation of cloud nuclei concentration in PDF component 1.
       if ( isigma_Ncn_1 > 0 ) then
          call stat_update_var_pt( isigma_Ncn_1, level, &
                                   sigma_x_1(iiPDF_Ncn), stats_zt )
       endif

       ! Standard deviation of cloud nuclei concentration in PDF component 2.
       if ( isigma_Ncn_2 > 0 ) then
          call stat_update_var_pt( isigma_Ncn_2, level, &
                                   sigma_x_2(iiPDF_Ncn), stats_zt )
       endif

       ! Correlation of w and chi (old s) in PDF component 1.
       ! This correlation should always be 0 because both the correlation
       ! between w and rt and the correlation of w and theta-l within each
       ! PDF component are defined to be 0 by CLUBB standards.
       if ( icorr_w_chi_1 > 0 ) then
          call stat_update_var_pt( icorr_w_chi_1, level, &
                                   corr_array_1(iiPDF_w,iiPDF_chi), stats_zt )
       endif

       ! Correlation of w and chi (old s) in PDF component 2.
       ! This correlation should always be 0 because both the correlation
       ! between w and rt and the correlation of w and theta-l within each
       ! PDF component are defined to be 0 by CLUBB standards.
       if ( icorr_w_chi_2 > 0 ) then
          call stat_update_var_pt( icorr_w_chi_2, level, &
                                   corr_array_2(iiPDF_w,iiPDF_chi), stats_zt )
       endif

       ! Correlation of w and eta (old t) in PDF component 1.
       ! This correlation should always be 0 because both the correlation
       ! between w and rt and the correlation of w and theta-l within each
       ! PDF component are defined to be 0 by CLUBB standards.
       if ( icorr_w_eta_1 > 0 ) then
          call stat_update_var_pt( icorr_w_eta_1, level, &
                                   corr_array_1(iiPDF_w,iiPDF_eta), stats_zt )
       endif

       ! Correlation of w and eta (old t) in PDF component 2.
       ! This correlation should always be 0 because both the correlation
       ! between w and rt and the correlation of w and theta-l within each
       ! PDF component are defined to be 0 by CLUBB standards.
       if ( icorr_w_eta_2 > 0 ) then
          call stat_update_var_pt( icorr_w_eta_2, level, &
                                   corr_array_2(iiPDF_w,iiPDF_eta), stats_zt )
       endif

       do ivar = iiPDF_Ncn+1, d_variables, 1

          ! Correlation (in-precip) of w and the precipitating hydrometeor
          ! in PDF component 1.
          if ( icorr_w_hm_1(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt( icorr_w_hm_1(pdf2hydromet_idx(ivar)), &
                                      level, corr_array_1(ivar,iiPDF_w), &
                                      stats_zt )
          endif

          ! Correlation (in-precip) of w and the precipitating hydrometeor
          ! in PDF component 2.
          if ( icorr_w_hm_2(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt( icorr_w_hm_2(pdf2hydromet_idx(ivar)), &
                                      level, corr_array_2(ivar,iiPDF_w), &
                                      stats_zt )
          endif

       enddo ! ivar = iiPDF_Ncn+1, d_variables, 1

       ! Correlation of w and N_cn in PDF component 1.
       if ( icorr_w_Ncn_1 > 0 ) then
          call stat_update_var_pt( icorr_w_Ncn_1, level, &
                                   corr_array_1(iiPDF_Ncn,iiPDF_w), stats_zt )
       endif

       ! Correlation of w and N_cn in PDF component 2.
       if ( icorr_w_Ncn_2 > 0 ) then
          call stat_update_var_pt( icorr_w_Ncn_2, level, &
                                   corr_array_2(iiPDF_Ncn,iiPDF_w), stats_zt )
       endif

       ! Correlation of chi (old s) and eta (old t) in PDF component 1 found in
       ! the correlation array.
       ! The true correlation of chi and eta in each PDF component is solved for
       ! by an equation and is part of CLUBB's PDF parameters.  However, there
       ! is an option in CLUBB, l_fix_chi_eta_correlations, that sets the
       ! component correlation of chi and eta to a constant, prescribed value
       ! because of SILHS.  The correlation of chi and eta in PDF component 1
       ! that is calculated by an equation is stored in stats as
       ! "corr_chi_eta_1".  Here, "corr_chi_eta_1_ca" outputs whatever value is
       ! found in the correlation array, whether or not it matches
       ! "corr_chi_eta_1".
       if ( icorr_chi_eta_1_ca > 0 ) then
          call stat_update_var_pt( icorr_chi_eta_1_ca, level, &
                                   corr_array_1(iiPDF_eta,iiPDF_chi), stats_zt )
       endif

       ! Correlation of chi (old s) and eta (old t) in PDF component 2 found in
       ! the correlation array.
       ! The true correlation of chi and eta in each PDF component is solved for
       ! by an equation and is part of CLUBB's PDF parameters.  However, there
       ! is an option in CLUBB, l_fix_chi_eta_correlations, that sets the
       ! component correlation of chi and eta to a constant, prescribed value
       ! because of SILHS.  The correlation of chi and eta in PDF component 2
       ! that is calculated by an equation is stored in stats as
       ! "corr_chi_eta_2".  Here, "corr_chi_eta_2_ca" outputs whatever value is
       ! found in the correlation array, whether or not it matches
       ! "corr_chi_eta_2".
       if ( icorr_chi_eta_2_ca > 0 ) then
          call stat_update_var_pt( icorr_chi_eta_2_ca, level, &
                                   corr_array_2(iiPDF_eta,iiPDF_chi), stats_zt )
       endif

       do ivar = iiPDF_Ncn+1, d_variables, 1

          ! Correlation (in-precip) of chi (old s) and the precipitating
          ! hydrometeor in PDF component 1.
          if ( icorr_chi_hm_1(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt( icorr_chi_hm_1(pdf2hydromet_idx(ivar)), &
                                      level, corr_array_1(ivar,iiPDF_chi), &
                                      stats_zt )
          endif

          ! Correlation (in-precip) of chi (old s) and the precipitating
          ! hydrometeor in PDF component 2.
          if ( icorr_chi_hm_2(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt( icorr_chi_hm_2(pdf2hydromet_idx(ivar)), &
                                      level, corr_array_2(ivar,iiPDF_chi), &
                                      stats_zt )
          endif

       enddo ! ivar = iiPDF_Ncn+1, d_variables, 1

       ! Correlation of chi (old s) and N_cn in PDF component 1.
       if ( icorr_chi_Ncn_1 > 0 ) then
          call stat_update_var_pt( icorr_chi_Ncn_1, level, &
                                   corr_array_1(iiPDF_Ncn,iiPDF_chi), stats_zt )
       endif

       ! Correlation of chi (old s) and N_cn in PDF component 2.
       if ( icorr_chi_Ncn_2 > 0 ) then
          call stat_update_var_pt( icorr_chi_Ncn_2, level, &
                                   corr_array_2(iiPDF_Ncn,iiPDF_chi), stats_zt )
       endif

       do ivar = iiPDF_Ncn+1, d_variables, 1

          ! Correlation (in-precip) of eta (old t) and the precipitating
          ! hydrometeor in PDF component 1.
          if ( icorr_eta_hm_1(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt( icorr_eta_hm_1(pdf2hydromet_idx(ivar)), &
                                      level, corr_array_1(ivar,iiPDF_eta), &
                                      stats_zt )
          endif

          ! Correlation (in-precip) of eta (old t) and the precipitating
          ! hydrometeor in PDF component 2.
          if ( icorr_eta_hm_2(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt( icorr_eta_hm_2(pdf2hydromet_idx(ivar)), &
                                      level, corr_array_2(ivar,iiPDF_eta), &
                                      stats_zt )
          endif

       enddo ! ivar = iiPDF_Ncn+1, d_variables, 1

       ! Correlation of eta (old t) and N_cn in PDF component 1.
       if ( icorr_eta_Ncn_1 > 0 ) then
          call stat_update_var_pt( icorr_eta_Ncn_1, level, &
                                   corr_array_1(iiPDF_Ncn,iiPDF_eta), stats_zt )
       endif

       ! Correlation of eta (old t) and N_cn in PDF component 2.
       if ( icorr_eta_Ncn_2 > 0 ) then
          call stat_update_var_pt( icorr_eta_Ncn_2, level, &
                                   corr_array_2(iiPDF_Ncn,iiPDF_eta), stats_zt )
       endif

       do ivar = iiPDF_Ncn+1, d_variables, 1

          ! Correlation (in-precip) of N_cn and the precipitating
          ! hydrometeor in PDF component 1.
          if ( icorr_Ncn_hm_1(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt( icorr_Ncn_hm_1(pdf2hydromet_idx(ivar)), &
                                      level, corr_array_1(ivar,iiPDF_Ncn), &
                                      stats_zt )
          endif

          ! Correlation (in-precip) of N_cn and the precipitating
          ! hydrometeor in PDF component 2.
          if ( icorr_Ncn_hm_2(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt( icorr_Ncn_hm_2(pdf2hydromet_idx(ivar)), &
                                      level, corr_array_2(ivar,iiPDF_Ncn), &
                                      stats_zt )
          endif

       enddo ! ivar = iiPDF_Ncn+1, d_variables, 1

       do ivar = iiPDF_Ncn+1, d_variables, 1
         do jvar = ivar+1, d_variables, 1

           ! Correlation (in-precip) of two different hydrometeors (hmx and
           ! hmy) in PDF component 1.
           if ( icorr_hmx_hmy_1(pdf2hydromet_idx(jvar),pdf2hydromet_idx(ivar)) &
                > 0 ) then
             call stat_update_var_pt( &
               icorr_hmx_hmy_1(pdf2hydromet_idx(jvar),pdf2hydromet_idx(ivar)), &
               level, corr_array_1(jvar,ivar), stats_zt )
           endif

           ! Correlation (in-precip) of two different hydrometeors (hmx and
           ! hmy) in PDF component 2.
           if ( icorr_hmx_hmy_2(pdf2hydromet_idx(jvar),pdf2hydromet_idx(ivar)) &
                > 0 ) then
             call stat_update_var_pt( &
               icorr_hmx_hmy_2(pdf2hydromet_idx(jvar),pdf2hydromet_idx(ivar)), &
               level, corr_array_2(jvar,ivar), stats_zt )
           endif

         enddo ! jvar = ivar+1, d_variables, 1
       enddo ! ivar = iiPDF_Ncn+1, d_variables, 1

    endif ! l_stats_samp


    return

  end subroutine pdf_param_hm_stats

  !=============================================================================
  subroutine pdf_param_ln_hm_stats( d_variables, level, mu_x_1_n, &
                                    mu_x_2_n, sigma_x_1_n, &
                                    sigma_x_2_n, corr_array_1_n, &
                                    corr_array_2_n, l_stats_samp )

    ! Description:
    ! Record statistics for normal space PDF parameters involving hydrometeors.

    ! References:
    !-----------------------------------------------------------------------

    use index_mapping, only: &
        pdf2hydromet_idx  ! Procedure(s)

    use corr_varnce_module, only: &
        iiPDF_w,   & ! Variable(s)
        iiPDF_chi, &
        iiPDF_eta, &
        iiPDF_Ncn

    use clubb_precision, only: &
        core_rknd   ! Variable(s)

    use stats_type_utilities, only: &
        stat_update_var_pt  ! Procedure(s)

    use stats_variables, only : &
        imu_hm_1_n,     & ! Variable(s)
        imu_hm_2_n,     &
        imu_Ncn_1_n,    &
        imu_Ncn_2_n,    &
        isigma_hm_1_n,  &
        isigma_hm_2_n,  &
        isigma_Ncn_1_n, &
        isigma_Ncn_2_n

    use stats_variables, only : &
        icorr_w_hm_1_n,    & ! Variables
        icorr_w_hm_2_n,    &
        icorr_w_Ncn_1_n,   &
        icorr_w_Ncn_2_n,   &
        icorr_chi_hm_1_n,  &
        icorr_chi_hm_2_n,  &
        icorr_chi_Ncn_1_n, &
        icorr_chi_Ncn_2_n, &
        icorr_eta_hm_1_n,  &
        icorr_eta_hm_2_n,  &
        icorr_eta_Ncn_1_n, &
        icorr_eta_Ncn_2_n, &
        icorr_Ncn_hm_1_n,  &
        icorr_Ncn_hm_2_n,  &
        icorr_hmx_hmy_1_n, &
        icorr_hmx_hmy_2_n, &
        stats_zt

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      d_variables, & ! Number of variables in the correlation array
      level          ! Vertical level index 

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(d_variables, d_variables), &
    intent(in) :: &
      corr_array_1_n, & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
      corr_array_2_n    ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.

    ! Local Variable
    integer :: ivar, jvar  ! Loop indices


    !!! Output the statistics for normal space hydrometeor PDF parameters.

    ! Statistics
    if ( l_stats_samp ) then

       do ivar = iiPDF_Ncn+1, d_variables, 1

          ! Mean (in-precip) of ln hm in PDF component 1.
          if ( imu_hm_1_n(pdf2hydromet_idx(ivar)) > 0 ) then
             if ( mu_x_1_n(ivar) > real( -huge( 0.0 ), kind = core_rknd ) ) then
                call stat_update_var_pt( imu_hm_1_n(pdf2hydromet_idx(ivar)), &
                                         level, mu_x_1_n(ivar), stats_zt )
             else
                ! When hm_1 is 0 (or below tolerance value), mu_hm_1_n is -inf,
                ! and is set to -huge for the default CLUBB kind.  Some
                ! compilers have issues outputting to stats files (in single
                ! precision) when the default CLUBB kind is in double precision.
                ! Set to -huge for single precision.
                call stat_update_var_pt( imu_hm_1_n(pdf2hydromet_idx(ivar)), &
                                         level, real( -huge( 0.0 ), &
                                                      kind = core_rknd ), &
                                         stats_zt )
             endif
          endif

          ! Mean (in-precip) of ln hm in PDF component 2.
          if ( imu_hm_2_n(pdf2hydromet_idx(ivar)) > 0 ) then
             if ( mu_x_2_n(ivar) > real( -huge( 0.0 ), kind = core_rknd ) ) then
                call stat_update_var_pt( imu_hm_2_n(pdf2hydromet_idx(ivar)), &
                                         level, mu_x_2_n(ivar), stats_zt )
             else
                ! When hm_2 is 0 (or below tolerance value), mu_hm_2_n is -inf,
                ! and is set to -huge for the default CLUBB kind.  Some
                ! compilers have issues outputting to stats files (in single
                ! precision) when the default CLUBB kind is in double precision.
                ! Set to -huge for single precision.
                call stat_update_var_pt( imu_hm_2_n(pdf2hydromet_idx(ivar)), &
                                         level, real( -huge( 0.0 ), &
                                                      kind = core_rknd ), &
                                         stats_zt )
             endif
          endif

       enddo ! ivar = iiPDF_Ncn+1, d_variables, 1

       ! Mean of ln N_cn in PDF component 1.
       if ( imu_Ncn_1_n > 0 ) then
          if ( mu_x_1_n(iiPDF_Ncn) &
                  > real( -huge( 0.0 ), kind = core_rknd ) ) then
             call stat_update_var_pt( imu_Ncn_1_n, level, &
                                      mu_x_1_n(iiPDF_Ncn), stats_zt )
          else
             ! When Ncnm is 0 (or below tolerance value), mu_Ncn_1_n is -inf,
             ! and is set to -huge for the default CLUBB kind.  Some compilers
             ! have issues outputting to stats files (in single precision) when
             ! the default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             call stat_update_var_pt( imu_Ncn_1_n, level, &
                                      real( -huge( 0.0 ), kind = core_rknd ), &
                                      stats_zt )
          endif
       endif

       ! Mean of ln N_cn in PDF component 2.
       if ( imu_Ncn_2_n > 0 ) then
          if ( mu_x_2_n(iiPDF_Ncn) &
                  > real( -huge( 0.0 ), kind = core_rknd ) ) then
             call stat_update_var_pt( imu_Ncn_2_n, level, &
                                      mu_x_2_n(iiPDF_Ncn), stats_zt )
          else
             ! When Ncnm is 0 (or below tolerance value), mu_Ncn_2_n is -inf,
             ! and is set to -huge for the default CLUBB kind.  Some compilers
             ! have issues outputting to stats files (in single precision) when
             ! the default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             call stat_update_var_pt( imu_Ncn_2_n, level, &
                                      real( -huge( 0.0 ), kind = core_rknd ), &
                                      stats_zt )
          endif
       endif

       do ivar = iiPDF_Ncn+1, d_variables, 1

          ! Standard deviation (in-precip) of ln hm in PDF component 1.
          if ( isigma_hm_1_n(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt( isigma_hm_1_n(pdf2hydromet_idx(ivar)), &
                                      level, sigma_x_1_n(ivar), stats_zt )
          endif

          ! Standard deviation (in-precip) of ln hm in PDF component 2.
          if ( isigma_hm_2_n(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt( isigma_hm_2_n(pdf2hydromet_idx(ivar)), &
                                      level, sigma_x_2_n(ivar), stats_zt )
          endif

       enddo ! ivar = iiPDF_Ncn+1, d_variables, 1

       ! Standard deviation of ln N_cn in PDF component 1.
       if ( isigma_Ncn_1_n > 0 ) then
          call stat_update_var_pt( isigma_Ncn_1_n, level, &
                                   sigma_x_1_n(iiPDF_Ncn), stats_zt )
       endif

       ! Standard deviation of ln N_cn in PDF component 2.
       if ( isigma_Ncn_2_n > 0 ) then
          call stat_update_var_pt( isigma_Ncn_2_n, level, &
                                   sigma_x_2_n(iiPDF_Ncn), stats_zt )
       endif

       do ivar = iiPDF_Ncn+1, d_variables, 1

          ! Correlation (in-precip) of w and ln hm in PDF component 1.
          if ( icorr_w_hm_1_n(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt( icorr_w_hm_1_n(pdf2hydromet_idx(ivar)), &
                                      level, corr_array_1_n(ivar,iiPDF_w), &
                                      stats_zt )
          endif

          ! Correlation (in-precip) of w and ln hm in PDF component 2.
          if ( icorr_w_hm_2_n(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt( icorr_w_hm_2_n(pdf2hydromet_idx(ivar)), &
                                      level, corr_array_2_n(ivar,iiPDF_w), &
                                      stats_zt )
          endif

       enddo ! ivar = iiPDF_Ncn+1, d_variables, 1

       ! Correlation of w and ln N_cn in PDF component 1.
       if ( icorr_w_Ncn_1_n > 0 ) then
          call stat_update_var_pt( icorr_w_Ncn_1_n, level, &
                                   corr_array_1_n(iiPDF_Ncn,iiPDF_w), stats_zt )
       endif

       ! Correlation of w and ln N_cn in PDF component 2.
       if ( icorr_w_Ncn_2_n > 0 ) then
          call stat_update_var_pt( icorr_w_Ncn_2_n, level, &
                                   corr_array_2_n(iiPDF_Ncn,iiPDF_w), stats_zt )
       endif

       do ivar = iiPDF_Ncn+1, d_variables, 1

          ! Correlation (in-precip) of chi (old s) and ln hm in PDF component 1.
          if ( icorr_chi_hm_1_n(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt(icorr_chi_hm_1_n(pdf2hydromet_idx(ivar)), &
                                     level, corr_array_1_n(ivar,iiPDF_chi), &
                                     stats_zt )
          endif

          ! Correlation (in-precip) of chi( old s) and ln hm in PDF component 2.
          if ( icorr_chi_hm_2_n(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt(icorr_chi_hm_2_n(pdf2hydromet_idx(ivar)), &
                                     level, corr_array_2_n(ivar,iiPDF_chi), &
                                     stats_zt )
          endif

       enddo ! ivar = iiPDF_Ncn+1, d_variables, 1

       ! Correlation of chi (old s) and ln N_cn in PDF component 1.
       if ( icorr_chi_Ncn_1_n > 0 ) then
          call stat_update_var_pt( icorr_chi_Ncn_1_n, level, &
                                   corr_array_1_n(iiPDF_Ncn,iiPDF_chi), &
                                   stats_zt )
       endif

       ! Correlation of chi(old s) and ln N_cn in PDF component 2.
       if ( icorr_chi_Ncn_2_n > 0 ) then
          call stat_update_var_pt( icorr_chi_Ncn_2_n, level, &
                                   corr_array_2_n(iiPDF_Ncn,iiPDF_chi), &
                                   stats_zt )
       endif

       do ivar = iiPDF_Ncn+1, d_variables, 1

          ! Correlation (in-precip) of eta (old t) and ln hm in PDF component 1.
          if ( icorr_eta_hm_1_n(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt(icorr_eta_hm_1_n(pdf2hydromet_idx(ivar)), &
                                     level, corr_array_1_n(ivar,iiPDF_eta), &
                                     stats_zt )
          endif

          ! Correlation (in-precip) of eta (old t) and ln hm in PDF component 2.
          if ( icorr_eta_hm_2_n(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt(icorr_eta_hm_2_n(pdf2hydromet_idx(ivar)), &
                                     level, corr_array_2_n(ivar,iiPDF_eta), &
                                     stats_zt )
          endif

       enddo ! ivar = iiPDF_Ncn+1, d_variables, 1

       ! Correlation of eta (old t) and ln N_cn in PDF component 1.
       if ( icorr_eta_Ncn_1_n > 0 ) then
          call stat_update_var_pt( icorr_eta_Ncn_1_n, level, &
                                   corr_array_1_n(iiPDF_Ncn,iiPDF_eta), &
                                   stats_zt )
       endif

       ! Correlation of eta (old t) and ln N_cn in PDF component 2.
       if ( icorr_eta_Ncn_2_n > 0 ) then
          call stat_update_var_pt( icorr_eta_Ncn_2_n, level, &
                                   corr_array_2_n(iiPDF_Ncn,iiPDF_eta), &
                                   stats_zt )
       endif

       do ivar = iiPDF_Ncn+1, d_variables, 1

          ! Correlation (in-precip) of ln N_cn and ln hm in PDF
          ! component 1.
          if ( icorr_Ncn_hm_1_n(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt(icorr_Ncn_hm_1_n(pdf2hydromet_idx(ivar)), &
                                     level, corr_array_1_n(ivar,iiPDF_Ncn), &
                                     stats_zt )
          endif

          ! Correlation (in-precip) of ln N_cn and ln hm in PDF
          ! component 2.
          if ( icorr_Ncn_hm_2_n(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt(icorr_Ncn_hm_2_n(pdf2hydromet_idx(ivar)), &
                                     level, corr_array_2_n(ivar,iiPDF_Ncn), &
                                     stats_zt )
          endif

       enddo ! ivar = iiPDF_Ncn+1, d_variables, 1

       do ivar = iiPDF_Ncn+1, d_variables, 1
         do jvar = ivar+1, d_variables, 1

           ! Correlation (in-precip) of ln hmx and ln hmy (two different
           ! hydrometeors) in PDF component 1.
           if (icorr_hmx_hmy_1_n(pdf2hydromet_idx(jvar),pdf2hydromet_idx(ivar))&
                > 0 ) then
             call stat_update_var_pt( &
             icorr_hmx_hmy_1_n(pdf2hydromet_idx(jvar),pdf2hydromet_idx(ivar)), &
             level, corr_array_1_n(jvar,ivar), stats_zt )
           endif

           ! Correlation (in-precip) of ln hmx and ln hmy (two different
           ! hydrometeors) in PDF component 2.
           if (icorr_hmx_hmy_2_n(pdf2hydromet_idx(jvar),pdf2hydromet_idx(ivar))&
                > 0 ) then
             call stat_update_var_pt( &
             icorr_hmx_hmy_2_n(pdf2hydromet_idx(jvar),pdf2hydromet_idx(ivar)), &
             level, corr_array_2_n(jvar,ivar), stats_zt )
           endif

         enddo ! jvar = ivar+1, d_variables, 1
       enddo ! ivar = iiPDF_Ncn+1, d_variables, 1

    endif ! l_stats_samp


    return

  end subroutine pdf_param_ln_hm_stats

  !=============================================================================
  subroutine pack_pdf_params( hm_1, hm_2, d_variables, &                 ! In
                              mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &    ! In
                              corr_array_1, corr_array_2, precip_frac, & ! In
                              precip_frac_1, precip_frac_2, &            ! In
                              hydromet_pdf_params )                      ! Out

    ! Description:
    ! Pack the standard means and variances involving hydrometeors, as well as a
    ! few other variables, into the structure hydromet_pdf_params.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one  ! Constant(s)

    use hydromet_pdf_parameter_module, only: &
        hydromet_pdf_parameter  ! Variable(s)

    use index_mapping, only: &
        hydromet2pdf_idx  ! Procedure(s)

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use corr_varnce_module, only: &
        iiPDF_w,   & ! Variable(s)
        iiPDF_chi, &
        iiPDF_eta, &
        iiPDF_Ncn

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(hydromet_dim), intent(in) :: &
      hm_1, & ! Mean of a precip. hydrometeor (1st PDF component)  [units vary]
      hm_2    ! Mean of a precip. hydrometeor (2nd PDF component)  [units vary]

    integer, intent(in) :: &
      d_variables    ! Number of variables in the mean/stdev arrays

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      mu_x_1,    & ! Mean array of PDF vars. (1st PDF component)    [units vary]
      mu_x_2,    & ! Mean array of PDF vars. (2nd PDF component)    [units vary]
      sigma_x_1, & ! Standard deviation array of PDF vars (comp. 1) [units vary]
      sigma_x_2    ! Standard deviation array of PDF vars (comp. 2) [units vary]

    real( kind = core_rknd ), dimension(d_variables,d_variables), &
    intent(in) :: &
      corr_array_1, & ! Correlation array of PDF vars. (comp. 1)    [-]
      corr_array_2    ! Correlation array of PDF vars. (comp. 2)    [-]

    real( kind = core_rknd ), intent(in) :: &
      precip_frac,   & ! Precipitation fraction (overall)           [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component) [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component) [-]

    ! Output Variable
    type(hydromet_pdf_parameter), intent(out) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]

    ! Local Variables
    integer :: ivar, jvar  ! Loop indices


    ! Pack remaining means and standard deviations into hydromet_pdf_params.
    do ivar = 1, hydromet_dim, 1

       ! Mean of a hydrometeor (overall) in the 1st PDF component.
       hydromet_pdf_params%hm_1(ivar) = hm_1(ivar)
       ! Mean of a hydrometeor (overall) in the 2nd PDF component.
       hydromet_pdf_params%hm_2(ivar) = hm_2(ivar)

       ! Mean of a hydrometeor (in-precip) in the 1st PDF component.
       hydromet_pdf_params%mu_hm_1(ivar) = mu_x_1(hydromet2pdf_idx(ivar))
       ! Mean of a hydrometeor (in-precip) in the 2nd PDF component.
       hydromet_pdf_params%mu_hm_2(ivar) = mu_x_2(hydromet2pdf_idx(ivar))

       ! Standard deviation of a hydrometeor (in-precip) in the
       ! 1st PDF component.
       hydromet_pdf_params%sigma_hm_1(ivar) = sigma_x_1(hydromet2pdf_idx(ivar))
       ! Standard deviation of a hydrometeor (in-precip) in the
       ! 2nd PDF component.
       hydromet_pdf_params%sigma_hm_2(ivar) = sigma_x_2(hydromet2pdf_idx(ivar))

       ! Correlation (in-precip) of w and a hydrometeor in the 1st PDF
       ! component.
       hydromet_pdf_params%corr_w_hm_1(ivar) &
       = corr_array_1( hydromet2pdf_idx(ivar), iiPDF_w )

       ! Correlation (in-precip) of w and a hydrometeor in the 2nd PDF
       ! component.
       hydromet_pdf_params%corr_w_hm_2(ivar) &
       = corr_array_2( hydromet2pdf_idx(ivar), iiPDF_w )

       ! Correlation (in-precip) of chi and a hydrometeor in the 1st PDF
       ! component.
       hydromet_pdf_params%corr_chi_hm_1(ivar) &
       = corr_array_1( hydromet2pdf_idx(ivar), iiPDF_chi )

       ! Correlation (in-precip) of chi and a hydrometeor in the 2nd PDF
       ! component.
       hydromet_pdf_params%corr_chi_hm_2(ivar) &
       = corr_array_2( hydromet2pdf_idx(ivar), iiPDF_chi )

       ! Correlation (in-precip) of eta and a hydrometeor in the 1st PDF
       ! component.
       hydromet_pdf_params%corr_eta_hm_1(ivar) &
       = corr_array_1( hydromet2pdf_idx(ivar), iiPDF_eta )

       ! Correlation (in-precip) of eta and a hydrometeor in the 2nd PDF
       ! component.
       hydromet_pdf_params%corr_eta_hm_2(ivar) &
       = corr_array_2( hydromet2pdf_idx(ivar), iiPDF_eta )

       ! Correlation (in-precip) of two hydrometeors, hmx and hmy, in the 1st
       ! PDF component.
       hydromet_pdf_params%corr_hmx_hmy_1(ivar,ivar) = one

       do jvar = ivar+1, hydromet_dim, 1

          hydromet_pdf_params%corr_hmx_hmy_1(jvar,ivar) &
          = corr_array_1( hydromet2pdf_idx(jvar), hydromet2pdf_idx(ivar) )

          hydromet_pdf_params%corr_hmx_hmy_1(ivar,jvar) &
          = hydromet_pdf_params%corr_hmx_hmy_1(jvar,ivar)

       enddo ! jvar = ivar+1, hydromet_dim, 1

       ! Correlation (in-precip) of two hydrometeors, hmx and hmy, in the 2nd
       ! PDF component.
       hydromet_pdf_params%corr_hmx_hmy_2(ivar,ivar) = one

       do jvar = ivar+1, hydromet_dim, 1

          hydromet_pdf_params%corr_hmx_hmy_2(jvar,ivar) &
          = corr_array_2( hydromet2pdf_idx(jvar), hydromet2pdf_idx(ivar) )

          hydromet_pdf_params%corr_hmx_hmy_2(ivar,jvar) &
          = hydromet_pdf_params%corr_hmx_hmy_2(jvar,ivar)

       enddo ! jvar = ivar+1, hydromet_dim, 1

    enddo ! ivar = 1, hydromet_dim, 1

    ! Mean of Ncn (overall) in the 1st PDF component.
    hydromet_pdf_params%mu_Ncn_1 = mu_x_1(iiPDF_Ncn)
    ! Mean of Ncn (overall) in the 2nd PDF component.
    hydromet_pdf_params%mu_Ncn_2 = mu_x_2(iiPDF_Ncn)

    ! Standard deviation of Ncn (overall) in the 1st PDF component.
    hydromet_pdf_params%sigma_Ncn_1 = sigma_x_1(iiPDF_Ncn)
    ! Standard deviation of Ncn (overall) in the 2nd PDF component.
    hydromet_pdf_params%sigma_Ncn_2 = sigma_x_2(iiPDF_Ncn)

    ! Precipitation fraction (overall).
    hydromet_pdf_params%precip_frac   = precip_frac
    ! Precipitation fraction (1st PDF component).
    hydromet_pdf_params%precip_frac_1 = precip_frac_1
    ! Precipitation fraction (2nd PDF component).
    hydromet_pdf_params%precip_frac_2 = precip_frac_2


    return

  end subroutine pack_pdf_params

  !=============================================================================
  elemental function compute_rtp2_from_chi( pdf_params, corr_chi_eta_1, &
                                            corr_chi_eta_2 ) &
  result( rtp2_zt_from_chi )

    ! Description:
    ! Compute the variance of rt from the distribution of chi and eta. The
    ! resulting variance will be consistent with CLUBB's extended PDF
    ! involving chi and eta, including if l_fix_chi_eta_correlations = .true. .

    ! References:
    !   None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd    ! Constant

    use pdf_utilities, only: &
        compute_variance_binormal   ! Procedure

    use constants_clubb, only: &
        one_half, & ! Constant(s)
        one,      &
        two

    use pdf_parameter_module, only: &
        pdf_parameter    ! Type

    implicit none

    ! Input Variables
    type(pdf_parameter), intent(in) :: &
      pdf_params

    real( kind = core_rknd ), intent(in) :: &
      corr_chi_eta_1, & ! Correlation of chi and eta in 1st PDF component [-]
      corr_chi_eta_2    ! Correlation of chi and eta in 2nd PDF component [-]

    ! Output Variable
    real( kind = core_rknd ) :: &
      rtp2_zt_from_chi    ! Grid-box variance of rtp2 on thermo. levels  [kg/kg]

    ! Local Variables
    real( kind = core_rknd ) :: &
      varnce_rt_1_zt_from_chi, varnce_rt_2_zt_from_chi

    real( kind = core_rknd ) :: &
      sigma_chi_1,         & ! Standard deviation of chi (1st PDF comp.) [kg/kg]
      sigma_chi_2,         & ! Standard deviation of chi (2nd PDF comp.) [kg/kg]
      sigma_eta_1,         & ! Standard deviation of eta (1st PDF comp.) [kg/kg]
      sigma_eta_2,         & ! Standard deviation of eta (2nd PDF comp.) [kg/kg]
      crt_1,               & ! Coef. of r_t in chi/eta eqns. (1st comp.) [-]
      crt_2,               & ! Coef. of r_t in chi/eta eqns. (2nd comp.) [-]
      rt_1,                & ! Mean of rt (1st PDF component)            [kg/kg]
      rt_2,                & ! Mean of rt (2nd PDF component)            [kg/kg]
      rtm,                 & ! Mean of rt (overall)                      [kg/kg]
      sigma_rt_1_from_chi, & ! Standard deviation of rt (1st PDF comp.)  [kg/kg]
      sigma_rt_2_from_chi, & ! Standard deviation of rt (2nd PDF comp.)  [kg/kg]
      mixt_frac              ! Weight of 1st gaussian PDF component      [-]

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    ! Enter some PDF parameters
    sigma_chi_1 = pdf_params%stdev_chi_1
    sigma_chi_2 = pdf_params%stdev_chi_2
    sigma_eta_1 = pdf_params%stdev_eta_1
    sigma_eta_2 = pdf_params%stdev_eta_2
    rt_1        = pdf_params%rt_1
    rt_2        = pdf_params%rt_2
    crt_1       = pdf_params%crt_1
    crt_2       = pdf_params%crt_2
    mixt_frac   = pdf_params%mixt_frac

    varnce_rt_1_zt_from_chi &
    = ( corr_chi_eta_1 * sigma_chi_1 * sigma_eta_1 &
        + one_half * sigma_chi_1**2 + one_half * sigma_eta_1**2 ) &
        / ( two * crt_1**2 )

    varnce_rt_2_zt_from_chi &
    = ( corr_chi_eta_2 * sigma_chi_2 * sigma_eta_2 &
        + one_half * sigma_chi_2**2 + one_half * sigma_eta_2**2 ) &
        / ( two * crt_2**2 )

    rtm = mixt_frac*rt_1 + (one-mixt_frac)*rt_2

    sigma_rt_1_from_chi = sqrt( varnce_rt_1_zt_from_chi )
    sigma_rt_2_from_chi = sqrt( varnce_rt_2_zt_from_chi )

    rtp2_zt_from_chi &
    = compute_variance_binormal( rtm, rt_1, rt_2, sigma_rt_1_from_chi, &
                                 sigma_rt_2_from_chi, mixt_frac )


    return

  end function compute_rtp2_from_chi

!===============================================================================

end module setup_clubb_pdf_params
