!-------------------------------------------------------------------------
! $Id: setup_clubb_pdf_params.F90 7379 2014-11-11 05:32:53Z bmg2@uwm.edu $
!===============================================================================
module setup_clubb_pdf_params

  implicit none

  private

  public :: setup_pdf_parameters, &
            compute_mean_stdev,   &
            normalize_mean_stdev, &
            compute_corr,         &
            normalize_corr

  private :: component_means_hydromet_orig, &
             component_means_hydromet_corr, &
             precip_fraction,               &
             component_mean_hm_ip,          &
             component_stdev_hm_ip,         &
             component_corr_w_x,            &
             component_corr_chi_eta,        &
             component_corr_w_hm_ip,        &
             component_corr_x_hm_ip,        &
             component_corr_hmx_hmy_ip,     &
             calc_corr_w_hm,                &
             pdf_param_hm_stats,            &
             pdf_param_ln_hm_stats,         &
             pack_pdf_params

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
  subroutine setup_pdf_parameters( nz, d_variables, dt, rho, &                 ! Intent(in)
                                   Nc_in_cloud, rcm, cloud_frac, &             ! Intent(in)
                                   ice_supersat_frac, hydromet, wphydrometp, & ! Intent(in)
                                   corr_array_cloud, corr_array_below, &       ! Intent(in)
                                   pdf_params, l_stats_samp, &                 ! Intent(in)
                                   hydrometp2, &                               ! Intent(inout)
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
        fstderr,        &
        zero_threshold

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s)

    use hydromet_pdf_parameter_module, only: &
        hydromet_pdf_parameter, &  ! Type
        init_hydromet_pdf_params   ! Procedure

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
        compute_variance_binormal

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
        ihm1,           & ! Variable(s)
        ihm2,           &
        iprecip_frac,   &
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
        assert_corr_symmetric,        & ! Procedure(s)
        sigma2_on_mu2_ip_array_cloud, & ! Variable(s)
        sigma2_on_mu2_ip_array_below, &
        iiPDF_Ncn,                    &
        iiPDF_chi,                    &
        iiPDF_eta

    use index_mapping, only: &
        hydromet2pdf_idx    ! Procedure(s)

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
      rho,         & ! Density                                         [kg/m^3]
      Nc_in_cloud    ! Mean (in-cloud) cloud droplet concentration     [num/kg]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rcm,               & ! Mean cloud water mixing ratio, < r_c >    [kg/kg]
      cloud_frac,        & ! Cloud fraction                            [-]
      ice_supersat_frac    ! Ice supersaturation fraction              [-]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet,    & ! Mean of hydrometeor, hm (overall) (t-levs.) [units]
      wphydrometp    ! Covariance < w'h_m' > (momentum levels)     [(m/s)units]

    real( kind = core_rknd ), dimension(d_variables,d_variables), &
    intent(in) :: &
      corr_array_cloud, & ! Prescribed correlation array in cloud      [-]
      corr_array_below    ! Prescribed correlation array below cloud   [-]

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params    ! PDF parameters                               [units vary]

    logical, intent(in) :: &
      l_stats_samp    ! Flag to sample statistics

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(inout) :: &
      hydrometp2    ! Variance of a hydrometeor (overall) (m-levs.)   [units^2]

    ! Output Variables
    real( kind = core_rknd ), dimension(d_variables,d_variables,nz), &
    intent(out) :: &
      corr_array_1_n, & ! Corr. array (normalized) of PDF vars. (comp. 1)    [-]
      corr_array_2_n    ! Corr. array (normalized) of PDF vars. (comp. 2)    [-]

    real( kind = core_rknd ), dimension(d_variables, nz), intent(out) :: &
      mu_x_1_n,    & ! Mean array (normalized) of PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normalized) of PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normalized) of PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normalized) of PDF vars (comp. 2) [u.v.]

    type(hydromet_pdf_parameter), dimension(nz), intent(out) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]

    real( kind = core_rknd ), dimension(d_variables,d_variables,nz), &
    intent(out) :: &
      corr_cholesky_mtx_1, & ! Transposed corr. cholesky matrix, 1st comp. [-]
      corr_cholesky_mtx_2    ! Transposed corr. cholesky matrix, 2nd comp. [-]

    ! Local Variables
    real( kind = dp ), dimension(d_variables,d_variables,nz) :: &
      corr_cholesky_mtx_1_dp, & ! Used for call to Cholesky_factor, requires dp
      corr_cholesky_mtx_2_dp

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
      Ncnm    ! Mean cloud nuclei concentration, < N_cn >        [num/kg]

    real( kind = core_rknd ), dimension(nz) ::  &
      wpchip_zm, & ! Covariance of chi and w (momentum levels)   [(m/s)(kg/kg)]
      wpNcnp_zm, & ! Covariance of N_cn and w (momentum levs.)   [(m/s)(num/kg)]
      wpchip_zt, & ! Covariance of chi and w on t-levs           [(m/s)(kg/kg)]
      wpNcnp_zt    ! Covariance of N_cn and w on t-levs          [(m/s)(num/kg)]

    real( kind = core_rknd ), dimension(nz,hydromet_dim) :: &
      hm1, & ! Mean of a precip. hydrometeor (1st PDF component)    [units vary]
      hm2    ! Mean of a precip. hydrometeor (2nd PDF component)    [units vary]

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

    real( kind = dp ), dimension(d_variables) :: &
      corr_array_scaling

    real( kind = core_rknd ), dimension(d_variables) :: &
      sigma2_on_mu2_ip_1, & ! Prescribed ratio array: sigma_hm_1^2/mu_hm_1^2 [-]
      sigma2_on_mu2_ip_2    ! Prescribed ratio array: sigma_hm_2^2/mu_hm_2^2 [-]

    real( kind = core_rknd ) :: &
      const_Ncnp2_on_Ncnm2, & ! Prescribed ratio of <Ncn'^2> to <Ncn>^2      [-]
      const_corr_chi_Ncn      ! Prescribed correlation of chi (old s) & Ncn  [-]

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

    integer, parameter :: &
      comp_means_hm_type = 1  ! Option used to calculated hm1 and hm2.

    integer :: pdf_idx  ! Index of precipitating hydrometeor in PDF array.

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

    ! Interpolate the variances (overall) of precipitating hydrometeors and the
    ! covariances (overall) of w and precipitating hydrometeors to thermodynamic
    ! grid levels.
    do i = 1, hydromet_dim, 1

       hydrometp2_zt(:,i)  = max( zm2zt( hydrometp2(:,i) ), zero_threshold )
       wphydrometp_zt(:,i) = zm2zt( wphydrometp(:,i) )

       ! When the mean value of a precipitating hydrometeor is below tolerance
       ! value, it is considered to have a value of 0, and the precipitating
       ! hydrometeor does not vary over the grid level.  The variance of that
       ! precipitating hydrometeor and any covariance involving that
       ! precipitating hydrometeor also have values of 0 at that grid level.
       do k = 1, nz, 1
          if ( hydromet(k,i) <= hydromet_tol(i) ) then
             hydrometp2_zt(k,i)  = zero
             wphydrometp_zt(k,i) = zero
          endif
       enddo ! k = 1, nz, 1

    enddo ! i = 1, hydromet_dim, 1

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

    ! Component mean values for r_r and N_r, and precipitation fraction.
    if ( l_use_precip_frac ) then

       if ( comp_means_hm_type == 1 ) then

          ! Original formulation to calculate hm1 and hm2 based on
          ! liquid water path at or above a given grid level.
          call component_means_hydromet_orig( nz, hydromet, rho, rc_1, rc_2, &
                                              mixt_frac, l_stats_samp, &
                                              hm1, hm2 )

       elseif ( comp_means_hm_type == 2 ) then

          ! New formulations to calculate hm1 and hm2 based on the overall
          ! correlation between w and the hydrometeor species.
          call component_means_hydromet_corr( nz, hydromet, wphydrometp_zt, &
                                              hydrometp2_zt, wp2_zt, &
                                              mixt_frac, l_stats_samp, &
                                              hm1, hm2 )

       else

          write(fstderr,*) "Invalid option to calculate hm1 and hm2."
          stop

       endif ! comp_means_hm_type

       call precip_fraction( nz, hydromet, hm1, hm2, &
                             cloud_frac, cloud_frac_1, mixt_frac, &
                             ice_supersat_frac, &
                             precip_frac, precip_frac_1, precip_frac_2 )

    else

       hm1 = hydromet
       hm2 = hydromet

       precip_frac   = one
       precip_frac_1 = one
       precip_frac_2 = one

    endif

    ! Calculate <N_cn> from Nc_in_cloud, whether Nc_in_cloud is predicted or
    ! based on a prescribed value, and whether the value is constant or varying
    ! over the grid level.
    if ( .not. l_const_Nc_in_cloud ) then
       ! Ncn varies at each vertical level.
       const_Ncnp2_on_Ncnm2 = sigma2_on_mu2_ip_array_cloud(iiPDF_Ncn)
    else  ! l_const_Nc_in_cloud
       ! Ncn is constant at each vertical level.
       const_Ncnp2_on_Ncnm2 = zero
    endif

    const_corr_chi_Ncn = corr_array_cloud(iiPDF_Ncn, iiPDF_chi)

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

       do i = 1, hydromet_dim, 1

          if ( ihm1(i) > 0 ) then
             ! Mean of the precipitating hydrometeor in PDF component 1.
             call stat_update_var( ihm1(i), hm1(:,i), stats_zt )
          endif

          if ( ihm2(i) > 0 ) then
             ! Mean of the precipitating hydrometeor in PDF component 2.
             call stat_update_var( ihm2(i), hm2(:,i), stats_zt )
          endif

       enddo ! i = 1, hydromet_dim, 1

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

       if ( rc_1(k) > rc_tol ) then
          sigma2_on_mu2_ip_1 = sigma2_on_mu2_ip_array_cloud
       else
          sigma2_on_mu2_ip_1 = sigma2_on_mu2_ip_array_below
       endif

       if ( rc_2(k) > rc_tol ) then
          sigma2_on_mu2_ip_2 = sigma2_on_mu2_ip_array_cloud
       else
          sigma2_on_mu2_ip_2 = sigma2_on_mu2_ip_array_below
       endif

       !!! Calculate the means and standard deviations involving PDF variables
       !!! -- w, chi, eta, N_cn, and any precipitating hydrometeors (hm in-precip)
       !!! -- for each PDF component.
       call compute_mean_stdev( Ncnm(k), rc_1(k), rc_2(k), &          ! Intent(in)
                                cloud_frac_1(k), cloud_frac_2(k), &   ! Intent(in)
                                hm1(k,:), hm2(k,:), &                 ! Intent(in)
                                precip_frac_1(k), precip_frac_2(k), & ! Intent(in)
                                sigma2_on_mu2_ip_array_cloud, &       ! Intent(in)
                                sigma2_on_mu2_ip_array_below, &       ! Intent(in)
                                pdf_params(k), d_variables, &         ! Intent(in)
                                mu_x_1, mu_x_2, &                     ! Intent(out)
                                sigma_x_1, sigma_x_2 )                ! Intent(out)


       !!! Calculate the normalized means and normalized standard deviations
       !!! involving precipitating hydrometeors (hm in-precip) and N_cn --
       !!! ln hm and ln N_cn -- for each PDF component.
       call normalize_mean_stdev( hm1(k,:), hm2(k,:), Ncnm(k), d_variables, &
                                  mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &
                                  sigma2_on_mu2_ip_1, sigma2_on_mu2_ip_2, &
                                  mu_x_1_n(:,k), mu_x_2_n(:,k), &
                                  sigma_x_1_n(:,k), sigma_x_2_n(:,k) )

       ! Calculate the overall variance of a precipitating hydrometeor (hm),
       ! <hm'^2>.
       do i = 1, hydromet_dim, 1

          if ( hydromet(k,i) > hydromet_tol(i) ) then

             ! There is some of the hydrometeor species found at level k.
             ! Calculate the variance (overall) of the hydrometeor.

             pdf_idx = hydromet2pdf_idx(i)

             hydrometp2_zt(k,i) &
             = calc_xp2( mu_x_1(pdf_idx), mu_x_2(pdf_idx), &
                         mu_x_1_n(pdf_idx,k), mu_x_2_n(pdf_idx,k), &
                         sigma_x_1(pdf_idx), sigma_x_2(pdf_idx), &
                         sigma_x_1_n(pdf_idx,k), sigma_x_2_n(pdf_idx,k), &
                         mixt_frac(k), precip_frac_1(k), precip_frac_2(k), &
                         hydromet(k,i) )

          else ! hydromet(k,i) = 0.

             hydrometp2_zt(k,i) = zero

          endif

          ! Statistics
          if ( l_stats_samp ) then

             if ( ihmp2_zt(i) > 0 ) then
                ! Variance (overall) of the hydrometeor, <hm'^2>.
                call stat_update_var_pt( ihmp2_zt(i), k, &
                                         hydrometp2_zt(k,i), stats_zt )
             endif

          endif ! l_stats_samp

          ! Clip the value of covariance <w'hm'> on thermodynamic levels.
          call clip_covar_level( clip_wphydrometp, k, l_first_clip_ts, &
                                 l_last_clip_ts, dt, wp2_zt(k), &
                                 hydrometp2_zt(k,i), &
                                 wphydrometp_zt(k,i), wphydrometp_chnge(k,i) )

       enddo ! i = 1, hydromet_dim, 1

       if ( l_diagnose_correlations ) then

          if ( rcm(k) > rc_tol ) then

             call diagnose_correlations( d_variables, corr_array_cloud, & ! Intent(in)
                                         corr_array_1 )                   ! Intent(out)

             call diagnose_correlations( d_variables, corr_array_cloud, & ! Intent(in)
                                         corr_array_2 )                   ! Intent(out)

          else

             call diagnose_correlations( d_variables, corr_array_below, & ! Intent(in)
                                         corr_array_1 )                   ! Intent(out)

             call diagnose_correlations( d_variables, corr_array_below, & ! Intent(in)
                                         corr_array_2 )                   ! Intent(out)

          endif

       else ! if .not. l_diagnose_correlations

          call compute_corr( wm_zt(k), rc_1(k), rc_2(k), cloud_frac_1(k), &
                             cloud_frac_2(k), wpchip_zt(k), wpNcnp_zt(k), &
                             sqrt(wp2_zt(k)), mixt_frac(k), precip_frac_1(k), &
                             precip_frac_2(k), wphydrometp_zt(k,:), &
                             mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &
                             corr_array_cloud, corr_array_below, &
                             pdf_params(k), d_variables, &
                             corr_array_1, corr_array_2 )

       endif ! l_diagnose_correlations

       !!! Statistics for standard PDF parameters involving hydrometeors.
       call pdf_param_hm_stats( d_variables, k, mu_x_1, mu_x_2, &
                                sigma_x_1, sigma_x_2, &
                                corr_array_1, corr_array_2, &
                                l_stats_samp )

       !!! Calculate the correlations involving the natural logarithm of
       !!! precipitating hydrometeors, ln hm (for example, ln r_r and ln N_r),
       !!! and ln N_cn for each PDF component.
       call normalize_corr( d_variables, sigma_x_1_n(:,k), sigma_x_2_n(:,k), &
                            sigma2_on_mu2_ip_1, sigma2_on_mu2_ip_2, &
                            corr_array_1, corr_array_2, &
                            corr_array_1_n(:,:,k), corr_array_2_n(:,:,k) )


       !!! Statistics for normalized PDF parameters involving hydrometeors.
       call pdf_param_ln_hm_stats( d_variables, k, mu_x_1_n(:,k), &
                                   mu_x_2_n(:,k), sigma_x_1_n(:,k), &
                                   sigma_x_2_n(:,k), corr_array_1_n(:,:,k), &
                                   corr_array_2_n(:,:,k), l_stats_samp )

       !!! Pack the PDF parameters
       call pack_pdf_params( hm1(k,:), hm2(k,:), d_variables, &            ! In
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
          call Cholesky_factor( d_variables, real(corr_array_1_n(:,:,k), kind = dp), & ! In
                                corr_array_scaling, corr_cholesky_mtx_1_dp(:,:,k), &  ! Out
                                l_corr_array_scaling ) ! Out

          call Cholesky_factor( d_variables, real(corr_array_2_n(:,:,k), kind = dp), & ! In
                                corr_array_scaling, corr_cholesky_mtx_2_dp(:,:,k), &  ! Out
                                l_corr_array_scaling ) ! Out
          corr_cholesky_mtx_1(:,:,k) = real( corr_cholesky_mtx_1_dp(:,:,k), kind = core_rknd )
          corr_cholesky_mtx_2(:,:,k) = real( corr_cholesky_mtx_2_dp(:,:,k), kind = core_rknd )
       endif

       ! For ease of use later in the code, we make the correlation arrays
       ! symmetrical
       call mirror_lower_triangular_matrix( d_variables, corr_array_1_n(:,:,k) )
       call mirror_lower_triangular_matrix( d_variables, corr_array_2_n(:,:,k) )

    enddo  ! Setup PDF parameters loop: k = 2, nz, 1

    ! Boundary condition for the variance (overall) of a hydrometeor, <hm'^2>,
    ! on thermodynamic grid levels at the lowest thermodynamic grid level, k = 1
    ! (which is below the model lower boundary).
    hydrometp2_zt(1,:) = hydrometp2_zt(2,:)

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
          call stat_update_var( irtp2_from_chi, zt2zm( rtp2_zt_from_chi ), stats_zm )
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
  subroutine component_means_hydromet_orig( nz, hydromet, rho, rc_1, rc_2, &
                                            mixt_frac, l_stats_samp, &
                                            hm1, hm2 )

    ! Description:
    ! The values of grid-level mean hydrometeor fields, <hm>, (for example,
    ! grid-level mean rain water mixing ratio, <r_r>, and grid-level mean rain
    ! drop concentration, <N_r>) are solved as part of the predictive equation
    ! set, based on the microphysics scheme.  However, CLUBB has a two component
    ! PDF.  The grid-level means of all hydrometeors must be subdivided into
    ! component means for each PDF component.  The equation relating the overall
    ! mean to the component means (for any hydrometeor, hm) is:
    !
    ! <hm> = a * hm1 + (1-a) * hm2;
    !
    ! where "a" is the mixture fraction (weight of the 1st PDF component), hm1
    ! is the mean of the hydrometeor in PDF component 1, and hm2 is the mean of
    ! the hydrometeor in PDF component 2.  This equation can be rewritten as:
    !
    ! <hm> = hm1 * ( a + (1-a) * hm2/hm1 ).
    !
    ! One way to solve for a component mean is to relate the ratio hm2/hm1 to
    ! other factors.  For now, this ratio based on other factors will be called
    ! hm2_hm1_ratio.  This ratio is entered into the above equation, allowing
    ! the equation to be solved for hm1:
    !
    ! hm1 = <hm> / ( a + (1-a) * hm2_hm1_ratio ).
    !
    ! Once hm1 has been solved for, hm2 can be solved by:
    !
    ! hm2 = ( <hm> - a * hm1 ) / (1-a).
    !
    ! At a grid level that is at least mostly cloudy, the simplest way to handle
    ! the ratio hm2/hm1 is to set it equal to the ratio rc_2/rc_1, where rc_1 is
    ! the mean cloud water mixing ratio in PDF component 1 and rc_2 is the mean
    ! cloud water mixing ratio in PDF component 2.  However, a precipitating
    ! hydrometeor sediments, falling from higher altitudes downwards.  The
    ! values of cloud water mixing ratio at a given grid level are not
    ! necessarily indicative of the amount of cloud water at higher levels.  A
    ! precipitating hydrometeor may have been already produced from cloud water
    ! at a higher altitude (vertical level) and fallen downwards to the given
    ! grid level.  Additionally, using grid-level cloud water mixing ratio
    ! especially does not work for a precipitating hydrometeor below cloud base
    ! (near the ground).
    !
    ! However, an alternative to component cloud water mixing ratio is component
    ! liquid water path.  Liquid water path accounts for the cloud water mixing
    ! ratio at the given grid level and at all grid levels higher in altitude.
    !
    ! In a stratocumulus case, the cloud water is spread out over all or almost
    ! all of the horizontal domain over a group of vertical levels.  At a given
    ! vertical level, the component mean cloud water mixing ratios should be
    ! almost equal, although usually slightly larger in the component with the
    ! larger component mean extended liquid water mixing ratio, s.  Likewise,
    ! the component liquid water paths should be nearly equal, with one
    ! component having a slightly larger liquid water path than the other
    ! component.
    !
    ! In a case of cumulus rising into stratocumulus, the upper portion of the
    ! cloudy domain will be very similar to the stratocumulus case described
    ! above, with similar cloud water mixing ratio and liquid water path
    ! results.  However, below the base of the stratocumulus clouds, where the
    ! cumulus clouds are found, the horizontal domain at each vertical level is
    ! only partially cloudy.  At these levels, any precipitating hydrometeor
    ! that was produced in the stratocumulus clouds above and fallen downwards
    ! is evaporating in the clear-air portions, while not evaporating in the
    ! cloudy portions.  Additionally, new amounts of a hydrometeor are being
    ! produced in the cloudy portions.  The amount of a hydrometeor in the
    ! cloudy portions becomes significantly larger than the amount of a
    ! hydrometeor in the clear portions.  The partially cloudy levels usually
    ! have a PDF where one component is significantly more saturated than the
    ! other component.  By the time the cloud base of the cumulus clouds is
    ! reached, the liquid water path for one PDF component should be
    ! significantly greater than the liquid water path for the other PDF
    ! component.
    !
    ! In a cumulus case, the horizontal domain at each level is usually partly
    ! cloudy.  Throughout the entire vertical domain, at every vertical level,
    ! one component usually is much more saturated than the other component.
    ! The liquid water path for one component is much greater than the liquid
    ! water path in the other component.  Likewise, a precipitating hydrometeor
    ! that is formed in cloud and falls preferentially through cloud will have
    ! large values in a portion of the horizontal domain and very small or 0
    ! values over the rest of the horizontal domain.
    !
    ! In order to estimate the amount of a hydrometeor in each PDF component,
    ! the ratio hm2/hm1 is going to be set equal to the ratio LWP2/LWP1, where
    ! LWP1 is the liquid water path in PDF component 1 and LWP2 is the liquid
    ! water path in PDF component 2.  LWP1 will be computed by taking the
    ! vertical integral of cloud water (see equation below) through the 1st PDF
    ! component from the given vertical level all the way to the top of the
    ! model.  LWP2 will be computed in the same manner.   It should be noted
    ! that this method makes the poor assumption that PDF component 1 always
    ! overlaps PDF component 1 between vertical levels, and likewise for PDF
    ! component 2.
    !
    ! Total liquid water path, LWP, is given by the following equation:
    !
    ! LWP(z) = INT(z:z_top) rho_a <r_c> dz';
    !
    ! where z is the altitude of the vertical level for which LWP is desired,
    ! z_top is the altitude at the top of the model domain, and z' is the
    ! dummy variable of integration.  Mean cloud water mixing ratio can be
    ! written as:
    !
    ! <r_c> = a * rc_1 + (1-a) * rc_2.
    !
    ! The equation for liquid water path is rewritten as:
    !
    ! LWP(z) = INT(z:z_top) rho_a ( a rc_1 + (1-a) rc_2 ) dz'; or
    !
    ! LWP(z) = INT(z:z_top) a rho_a rc_1 dz'
    !          + INT(z:z_top) (1-a) rho_a rc_2 dz'.
    !
    ! This can be rewritten as:
    !
    ! LWP(z) = LWP1(z) + LWP2(z);
    !
    ! where:
    !
    ! LWP1(z) = INT(z:z_top) a rho_a rc_1 dz'; and
    ! LWP2(z) = INT(z:z_top) (1-a) rho_a rc_2 dz'.
    !
    ! The trapezoidal rule will be used to numerically integrate for LWP1
    ! and LWP2.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable(s)

    use constants_clubb, only: &
        one,      & ! Constant(s)
        one_half, &
        zero

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use array_index, only: &
        hydromet_tol  ! Variable(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    use stats_type_utilities, only: &
        stat_update_var  ! Procedure(s)

    use stats_variables, only : &
        iLWP1, & ! Variable(s)
        iLWP2, &
        stats_zt

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz    ! Number of model vertical grid levels

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet    ! Mean of hydrometeor, hm (overall)           [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rho,       & ! Air density                                        [kg/m^3]
      rc_1,      & ! Mean cloud water mixing ratio (1st PDF component)  [kg/kg]
      rc_2,      & ! Mean cloud water mixing ratio (2nd PDF component)  [kg/kg]
      mixt_frac    ! Mixture fraction                                   [-]

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.

    ! Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      hm1, & ! Mean of hydrometeor (1st PDF component)          [units vary]
      hm2    ! Mean of hydrometeor (2nd PDF component)          [units vary]

    ! Local Variable
    real( kind = core_rknd ), dimension(nz) :: &
      LWP1, & ! Liquid water path (1st PDF component) on thermo. levs.  [kg/m^2]
      LWP2    ! Liquid water path (2nd PDF component) on thermo. levs.  [kg/m^2]

    integer :: k, i  ! Array index

    real( kind = core_rknd ), parameter :: &
      LWP_tol = 5.0e-7_core_rknd  ! Tolerance value for component LWP


    !!! Compute component liquid water paths using trapezoidal rule for
    !!! numerical integration.

    ! At the uppermost thermodynamic level (k = nz), use the trapezoidal rule:
    !
    ! 0.5 * (integrand_a + integrand_b) * delta_z,
    !
    ! where integrand_a is the integrand at thermodynamic level k = nz,
    ! integrand_b is the integrand at momentum level k = nz (model upper
    ! boundary), and delta_z = zm(nz) - zt(nz).  At the upper boundary, r_c is
    ! set to 0, and the form of the trapezoidal rule is simply:
    !
    ! 0.5 * integrand_a * delta_z.

    ! Liquid water path in PDF component 1.
    LWP1(nz) &
    = one_half * mixt_frac(nz) * rho(nz) * rc_1(nz) * ( gr%zm(nz) - gr%zt(nz) )

    ! Liquid water path in PDF component 2.
    LWP2(nz) &
    = one_half * ( one - mixt_frac(nz) ) * rho(nz) * rc_2(nz) &
      * ( gr%zm(nz) - gr%zt(nz) )

    ! At all other thermodynamic levels, compute liquid water path using the
    ! trapezoidal rule:
    !
    ! 0.5 * (integrand_a + integrand_b) * delta_z,
    !
    ! where integrand_a is the integrand at thermodynamic level k, integrand_b
    ! is the integrand at thermodynamic level k+1, and
    ! delta_z = zt(k+1) - zt(k), or 1/invrs_dzm(k).  The total for the segment
    ! is added to the sum total of all higher vertical segments to compute the
    ! total vertical integral.
    do k = nz-1, 1, -1

       ! Liquid water path in PDF component 1.
       LWP1(k) &
       = LWP1(k+1) &
         + one_half * ( mixt_frac(k+1) * rho(k+1) * rc_1(k+1) &
                        + mixt_frac(k) * rho(k) * rc_1(k) ) / gr%invrs_dzm(k)

       ! Liquid water path in PDF component 2.
       LWP2(k) &
       = LWP2(k+1) &
         + one_half * ( ( one - mixt_frac(k+1) ) * rho(k+1) * rc_2(k+1) &
                        + ( one - mixt_frac(k) ) * rho(k) * rc_2(k) ) &
           / gr%invrs_dzm(k)

    enddo ! k = nz-1, 1, -1


    !!! Find hm1 and hm2 based on the ratio of LWP2/LWP1, such that:
    !!! hm2/hm1 ( = rr2/rr1 = Nr2/Nr1, etc. ) = LWP2/LWP1.
    do i = 1, hydromet_dim, 1

       do k = 1, nz, 1

          !!! Calculate the component means for the hydrometeor.
          if ( hydromet(k,i) > hydromet_tol(i) ) then

             if ( LWP1(k) <= LWP_tol .and. LWP2(k) <= LWP_tol ) then

                ! Both LWP1 and LWP2 are 0 (or an insignificant amount).
                !
                ! The hydrometeor is found at this level, yet there is no cloud
                ! at or above the current level.  This is usually due to a
                ! numerical artifact.  For example, the hydrometeor is diffused
                ! above cloud top.  Simply set each component mean equal to the
                ! overall mean.
                hm1(k,i) = hydromet(k,i)
                hm2(k,i) = hydromet(k,i)

             elseif ( LWP1(k) > LWP_tol .and. LWP2(k) <= LWP_tol ) then

                ! LWP1 is (significantly) greater than 0, while LWP2 is 0 (or an
                ! insignificant amount).
                !
                ! The hydrometeor is found at this level, and all cloud water at
                ! or above this level is found in the 1st PDF component.  All of
                ! the hydrometeor is found in the 1st PDF component.
                hm1(k,i) = hydromet(k,i) / mixt_frac(k)
                hm2(k,i) = zero

             elseif ( LWP2(k) > LWP_tol .and. LWP1(k) <= LWP_tol ) then

                ! LWP2 is (significantly) greater than 0, while LWP1 is 0 (or an
                ! insignificant amount).
                !
                ! The hydrometeor is found at this level, and all cloud water at
                ! or above this level is found in the 2nd PDF component.  All of
                ! the hydrometeor is found in the 2nd PDF component.
                hm1(k,i) = zero
                hm2(k,i) = hydromet(k,i) / ( one - mixt_frac(k) )

             else ! LWP1(k) > LWP_tol and LWP2(k) > LWP_tol

                ! Both LWP1 and LWP2 are (significantly) greater than 0.
                !
                ! The hydrometeor is found at this level, and there is
                ! sufficient cloud water at or above this level in both PDF
                ! components to find the hydrometeor in both PDF components.
                ! Delegate the hydrometeor between the 1st and 2nd PDF
                ! components according to the above equations.
                hm1(k,i) &
                = hydromet(k,i) &
                  / ( mixt_frac(k) + ( one - mixt_frac(k) ) * LWP2(k)/LWP1(k) )

                hm2(k,i) &
                = ( hydromet(k,i) - mixt_frac(k) * hm1(k,i) ) &
                  / ( one - mixt_frac(k) )

                if ( hm1(k,i) <= hydromet_tol(i) ) then

                   ! The mean value of the hydrometeor within the 1st PDF
                   ! component is below the tolerance value for the hydrometeor.
                   ! It is considered to have a value of 0.  All the the
                   ! hydrometeor is found within the 2nd PDF component.
                   hm1(k,i) = zero
                   hm2(k,i) = hydromet(k,i) / ( one - mixt_frac(k) )

                elseif ( hm2(k,i) <= hydromet_tol(i) ) then

                   ! The mean value of the hydrometeor within the 2nd PDF
                   ! component is below the tolerance value for the hydrometeor.
                   ! It is considered to have a value of 0.  All the the
                   ! hydrometeor is found within the 1st PDF component.
                   hm1(k,i) = hydromet(k,i) / mixt_frac(k)
                   hm2(k,i) = zero

                endif

             endif


          else ! hydromet(k,i) <= hydromet_tol(i)

             ! The overall hydrometeor is either 0 or below tolerance value (any
             ! postive value is considered to be a numerical artifact).  Simply
             ! set each pdf component mean equal to 0.
             hm1(k,i) = zero
             hm2(k,i) = zero

          endif

       enddo ! k = 1, nz, 1

    enddo ! i = 1, hydromet_dim, 1


    ! Statistics
    if ( l_stats_samp ) then

       if ( iLWP1 > 0 ) then
          ! Liquid water path in PDF component 1.
          call stat_update_var( iLWP1, LWP1, stats_zt )
       endif

       if ( iLWP2 > 0 ) then
          ! Liquid water path in PDF component 2.
          call stat_update_var( iLWP2, LWP2, stats_zt )
       endif
       
    endif


    return

  end subroutine component_means_hydromet_orig

  !=============================================================================
  subroutine component_means_hydromet_corr( nz, hydromet, wphydrometp_zt, &
                                            hydrometp2_zt, wp2_zt, &
                                            mixt_frac, l_stats_samp, &
                                            hm_1, hm_2 )

    ! Description:
    ! The values of grid-level mean hydrometeor fields, <hm>, (for example,
    ! grid-level mean rain water mixing ratio, <r_r>, and grid-level mean rain
    ! drop concentration, <N_r>) are solved as part of the predictive equation
    ! set, based on the microphysics scheme.  However, CLUBB has a two component
    ! PDF.  The grid-level means of all hydrometeors must be subdivided into
    ! component means for each PDF component.  The equation relating the overall
    ! mean to the component means (for any hydrometeor, hm) is:
    !
    ! <hm> = a * hm_1 + (1-a) * hm_2;
    !
    ! where "a" is the mixture fraction (weight of the 1st PDF component), hm_1
    ! is the mean of the hydrometeor in PDF component 1, and hm_2 is the mean of
    ! the hydrometeor in PDF component 2.  Both of these component means include
    ! any precipitationless regions in each PDF component (when component
    ! precipitation fraction < 1).
    !
    ! The challenge is to divide <hm> into hm_1 and hm_2.  One way to do this is
    ! to base this on the overall correlation of vertical velocity, w, and the
    ! hydrometeor, hm.  When the overall correlation of w and hm is positive,
    ! hm_1 > hm_2.  Likewise, when the overall correlation of w and hm is
    ! negative, hm_1 < hm_2.  When the overall correlation of w and hm is 0,
    ! hm_1 = hm_2.  This method has the following advantages.
    !
    ! 1) The main advantage is that this method aids the realizability of the
    !    multivariate PDF in each PDF component, when the PDF is considered in
    !    conjunction with the value of <w'hm'> produced by the microphysics.
    !    The ith PDF component, within-precip. correlation of w and hm
    !    (corr_w_hm_i) can be calculated based on the overall covariance of w
    !    and hm (<w'hm'>) and the other PDF parameters involving w and hm.  The
    !    value of <w'hm'> is produced when <hm> is advanced one model timestep
    !    in CLUBB's microphysics.  In the past, the calculated value of
    !    corr_w_hm_i has been unrealizable at some grid levels.  This was
    !    primarily due to the following issue.  The code that calculated hm_1
    !    and hm_2 based on integrated rc in each component always placed a great
    !    majority or all of the hydrometeor in PDF component 1 in cumulus cases,
    !    regardless of grid level.  Even in stratocumulus cases, hm_1 > hm_2.
    !    In CLUBB, the 1st PDF component mean of w (w_1) is defined around the
    !    updraft, while the 2nd PDF component mean of w (w_2) is defined around
    !    the downdraft, which means that w_1 is always greater than or equal to
    !    w_2.  Since hm_1 > hm_2 and w_1 > w_2, the means of the components are
    !    naturally associated with a positive value of covariance <w'hm'>.  In
    !    the scenario where <w'hm'> is negative at a grid level, the
    !    within-component correlation corr_w_hm_i needed to be so negative to
    !    produce <w'hm'>, because of the hm_1/hm_2 and w1/w2 values, that it had
    !    to be less than -1, which produces an unrealizable PDF.
    !
    !    In this method, when microphysics produces a positive value of <w'hm'>,
    !    hm_1 > hm_2, and when microphysics produces a negtive value of <w'hm'>,
    !    hm_1 < hm_2.  When microphysics produces a <w'hm'> of 0, hm_1 = hm_2.
    !    This will help keep the calculated values of corr_w_hm_i at realizable
    !    values.
    !
    ! 2) I have proposed a method to determine hm_1 amd hm_2 based on wind shear
    !    (the change in speed and/or direction of horizontal winds with
    !    altitude), which causes separation of updrafts and downdrafts in
    !    nature.  I am convinced that the profiles of hm_1 and hm_2 produced by
    !    this overall-correlation-based method would be roughly similar to those
    !    produced by a wind-shear-based method.
    !
    ! 3) This method involves minimal calculations and is conceptually simple.
    !    Any shear-based method would be conceptually complicated, and most
    !    likely involve more calculation.  This method is also a bit less
    !    numerically expensive than the integrated rc method, which involved
    !    extra vertical looping.
    !
    ! The value of hm_1 and hm_2 will be calculated by the following method.
    !
    ! When the overall correlation of w and hm (based on <w'hm'> provided by the
    ! microphysics) is exactly 1, all the hydrometeor will be found in the 1st
    ! PDF component.  In this scenario, hm_1 = <hm>/a and hm_2 = 0.  Likewise,
    ! when the overall correlation of w and hm is exactly -1, all the
    ! hydrometeor will be found in the 2nd PDF component.  In this scenario,
    ! hm_1 = 0 and hm_2 = <hm>/(1-a).  When the overall correlation of w and hm
    ! is exactly 0, hm_1 = hm_2 = <hm>.
    !
    ! What happens when the overall correlation of w and hm is at some
    ! intermediate value?  A function, based on the value of corr_w_hm_overall,
    ! is used to connect the three points listed above.  The function, when
    ! written to calculate hm_1, must be MONOTONICALLY INCREASING over the
    ! domain -1 <= corr_w_hm_overall <= 1.  A quadratic polynomial used to
    ! connect the three points for hm_1 (those points are (-1,0), (0,<hm>), and
    ! (1,<hm>/a)) is only monotonically increasing over the domain when
    ! 0.25 <= a <= 0.75.  Since "a" is often outside that range in highly skewed
    ! cases, a quadratic polynomial cannot be used.  Other options include a
    ! power-law fit and a piecewise linear fit.  I have opted for the power-law
    ! fit.
    !
    ! A power law is given by the equation:
    !
    ! hm_1 = A * x^kappa.
    !
    ! Since hm_1 is based on corr_w_hm_overall, and hm_1 must be positive and
    ! monotonically increasing over the domain -1 <= corr_w_hm_overall <= 1,
    ! the coefficient A, the value of x, and the exponent kappa must be
    ! positive.  The equation for hm_1 is given by:
    !
    ! hm_1 = A * ( 1 + corr_w_hm_overall )^kappa.
    !
    ! The three points listed above result in:
    !
    ! <hm>/a = A * ( 1 + 1 )^kappa;
    ! <hm>   = A * ( 1 + 0 )^kappa;
    ! 0      = A * ( 1 + -1 )^kappa;
    !
    ! and since 1^kappa = 1:
    !
    ! <hm>/a = A * 2^kappa;
    ! <hm>   = A;
    ! 0      = A * 0^kappa;
    !
    ! which further simplifies to (since A = <hm>):
    !
    ! <hm>/a = <hm> * 2^kappa;
    ! 0      = <hm> * 0^kappa.
    !
    ! As long as kappa > 0, the equation will work out for 0 = <hm> * 0^kappa.
    ! This leaves solving for kappa to <hm>/a = <hm> * 2^kappa.  Dividing both
    ! sides by <hm>, the equation reduces to:
    !
    ! 1/a = 2^kappa;
    ! ln( 1/a ) = ln( 2^kappa );
    ! ln( 1/a ) = kappa * ln( 2 ); and
    ! kappa = ln( 1/a ) / ln( 2 ).
    !
    ! The equation for hm_1 becomes:
    !
    ! hm_1 = <hm> * ( 1 + corr_w_hm_overall )^( ln( 1/a ) / ln( 2 ) ).
    !
    ! Since 1/a > 1, ln( 1/a ) > 0, and the exponent is always positive.
    !
    ! However, there have been issues in the past (in both the accuracy of the
    ! PDF shape and in problems resulting from extreme SILHS sample points) when
    ! all the hydrometeor is found in one PDF component.  In this method, that
    ! will occur when w and hm are perfectly correlated or perfectly
    ! anti-correlated.  In order to allow for a limit to be placed on the hm_1
    ! and hm_2 distribution, a new tunable parameter is introduced.  The method
    ! becomes the following.
    !
    ! First, calculate the overall correlation of w and hm:
    !
    ! corr_w_hm_overall = <w'hm'> / ( sqrt( <w'^2> ) * sqrt( <hm'^2> ) ).
    !
    ! Then, find the adjusted overall correlation of w and hm:
    !
    ! corr_w_hm_overall_adj = coef_hm_1_hm_2_corr_adj * corr_w_hm_overall;
    !
    ! where 0 <= coef_hm_1_hm_2_corr_adj <= 1.  Here, coef_hm_1_hm_2_corr_adj
    ! is a tunable parameter.  When it is equal to 1, the adjusted overall
    ! correlation is equal to the overall correlation.  When it is equal to 0,
    ! the adjusted correlation is always 0, resulting in hm_1 = hm_2.
    !
    ! Next, calculated hm_1 based on the adjusted overall correlation:
    !
    ! hm_1 = <hm> * ( 1 + corr_w_hm_overall_adj )^( ln( 1/a ) / ln( 2 ) ).
    !
    ! Once hm_1 has been solved for, hm_2 can be solved by:
    !
    ! hm_2 = ( <hm> - a * hm_1 ) / (1-a).

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two,   & ! Constant(s)
        one,   &
        zero,  &
        w_tol

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use array_index, only: &
        hydromet_tol  ! Variable(s)

    use parameters_tunable, only: &
        coef_hm_1_hm_2_corr_adj  ! Variable(s)

    use stats_type_utilities, only: &
        stat_update_var_pt  ! Procedure(s)

    use stats_variables, only: &
        icorr_w_hm_ov_adj, & ! Variable(s)
        stats_zt

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz    ! Number of model vertical grid levels

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet,       & ! Mean of hydrometeor, hm (overall) (t-levs.)    [units]
      wphydrometp_zt, & ! Covariance of w and hm interp. to t-levs. [(m/s)units]
      hydrometp2_zt     ! Variance of hm (overall) interp. to t-levs.  [units^2]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      wp2_zt,    & ! Variance of w, <w'^2> (interp. to t-levs.)  [m^2/s^2]
      mixt_frac    ! Mixture fraction                            [-]

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.

    ! Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      hm_1, & ! Mean of hydrometeor (1st PDF component)          [units vary]
      hm_2    ! Mean of hydrometeor (2nd PDF component)          [units vary]

    ! Local Variables
    real( kind = core_rknd ) :: &
      corr_w_hm_overall, & ! Overall correlation of w and hm                [-]
      kappa_exp            ! Exponent kappa = ln( 1/mixt_frac ) / ln( 2 )   [-]

    real( kind = core_rknd ), dimension(nz,hydromet_dim) :: &
      corr_w_hm_overall_adj    ! Adjusted overall correlation of w and hm   [-]

    real( kind = core_rknd ) :: &
      ln_2    ! Natural logarithm of 2                                      [-]

    integer :: k, i  ! Loop indices


    ! Initialize the adjusted overall correlation of w and the hydrometeor to 0.
    corr_w_hm_overall_adj = zero

    ! Calculate the Natural logarithm of 2.
    ln_2 = log( two )

    !!! Find hm_1 and hm_2 based on the overall correlation of w and the
    !!! hydrometeor, corr_w_hm_overall.
    do k = 1, nz, 1

       ! Calculate the value of the exponent kappa, where
       ! kappa = ln( 1/mixt_frac ) / ln( 2 ).
       ! This exponent is the same regardless of the hydrometeor type.
       kappa_exp = log( one / mixt_frac(k) ) / ln_2

       do i = 1, hydromet_dim, 1

          !!! Calculate the component means for the hydrometeor.
          if ( hydromet(k,i) > hydromet_tol(i) ) then

             ! Calculate the overall calculation of w and hm.
             if ( sqrt( wp2_zt(k) ) > w_tol .and. &
                  sqrt( hydrometp2_zt(k,i) ) > hydromet_tol(i) ) then

                ! Both w and the hydrometeor vary at this grid level.  The
                ! overall correlation between them is defined.
                ! Calculate the overall correlation of w and hm.
                corr_w_hm_overall &
                = wphydrometp_zt(k,i) &
                  / ( sqrt( wp2_zt(k) ) * sqrt( hydrometp2_zt(k,i) ) )

                ! Keep values realizable.
                if ( corr_w_hm_overall > one ) then
                   corr_w_hm_overall = one
                elseif ( corr_w_hm_overall < -one ) then
                   corr_w_hm_overall = -one
                endif

             else ! sqrt(wp2_zt) <= w_tol or sqrt(hydrometp2_zt) <= hydromet_tol

                ! Either w or the hydrometeor is constant at this grid level.
                ! This means that <w'hm'> must also have a value of 0, making
                ! the correlation undefined.  In the scenario that <hm'^2> = 0,
                ! the hydrometeor is constant at this grid level, which means
                ! that hm_1 = hm_2 = <hm>.  This is also the result when the
                ! correlation has a value of 0.  So, set the correlation to 0 in
                ! order to achieve the result hm_1 = hm_2 = <hm>.  In the
                ! scenario that <w'^2> = 0, w is constant at is grid level.
                ! To simplify matters, the undefined correlation will be set to
                ! 0 in order to produce hm_1 = hm_2 = <hm>.
                corr_w_hm_overall = zero

             endif ! sqrt(wp2_zt) > w_tol and sqrt(hydrometp2_zt) > hydromet_tol

             ! Calculate the adjusted overall correlation of w and hm.
             corr_w_hm_overall_adj(k,i) &
             = coef_hm_1_hm_2_corr_adj * corr_w_hm_overall

             ! Calculate the mean of the hydrometeor in the 1st PDF component.
             hm_1(k,i) &
             = hydromet(k,i) * ( one + corr_w_hm_overall_adj(k,i) )**kappa_exp

             ! Calculate the mean of the hydrometeor in the 2nd PDF component.
             hm_2(k,i) &
             = ( hydromet(k,i) - mixt_frac(k) * hm_1(k,i) ) &
               / ( one - mixt_frac(k) )

             if ( hm_1(k,i) < zero ) then

                ! The mean value of the hydrometeor within the 1st PDF component
                ! is below 0 due to numerical roundoff error.  Reset its value
                ! to 0.  All the the hydrometeor is found within the 2nd PDF
                ! component.
                hm_1(k,i) = zero
                hm_2(k,i) = hydromet(k,i) / ( one - mixt_frac(k) )

             elseif ( hm_2(k,i) < zero ) then

                ! The mean value of the hydrometeor within the 2nd PDF component
                ! is below 0 due to numerical roundoff error.  Reset its value
                ! to 0.  All the the hydrometeor is found within the 1st PDF
                ! component.
                hm_1(k,i) = hydromet(k,i) / mixt_frac(k)
                hm_2(k,i) = zero

             endif


          else ! hydromet(k,i) <= hydromet_tol(i)

             ! The overall hydrometeor is either 0 or below tolerance value (any
             ! postive value is considered to be a numerical artifact).  Simply
             ! set each PDF component mean equal to 0.  These values will not
             ! play into any further calculations.
             hm_1(k,i) = zero
             hm_2(k,i) = zero

          endif  ! hydromet(k,i) > hydromet_tol(i)

          ! Statistics
          if ( l_stats_samp ) then

             if ( icorr_w_hm_ov_adj(i) > 0 ) then
                ! Adjusted overall correlation of w and hm.
                call stat_update_var_pt( icorr_w_hm_ov_adj(i), k, &
                                         corr_w_hm_overall_adj(k,i), stats_zt )
             endif
       
          endif ! l_stats_samp

       enddo ! i = 1, hydromet_dim, 1

    enddo ! k = 1, nz, 1


    return

  end subroutine component_means_hydromet_corr

  !=============================================================================
  subroutine precip_fraction( nz, hydromet, hm1, hm2, &
                              cloud_frac, cloud_frac_1, mixt_frac, &
                              ice_supersat_frac, &
                              precip_frac, precip_frac_1, precip_frac_2 )

    ! Description:
    ! Determines (overall) precipitation fraction over the horizontal domain, as
    ! well as the precipitation fraction within each PDF component, at every
    ! vertical grid level.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,            & ! Constant(s)
        zero,           &
        cloud_frac_min, &
        fstderr

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use array_index, only: &
        l_mix_rat_hm, &  ! Variable(s)
        hydromet_tol  ! Variable(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz          ! Number of model vertical grid levels

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet, & ! Mean of hydrometeor, hm (overall)           [units vary]
      hm1,      & ! Mean of hydrometeor (1st PDF component)     [units vary]
      hm2         ! Mean of hydrometeor (2nd PDF component)     [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      cloud_frac,        &  ! Cloud fraction (overall)                     [-] 
      cloud_frac_1,      &  ! Cloud fraction (1st PDF component)           [-]
      mixt_frac,         &  ! Mixture fraction                             [-]
      ice_supersat_frac     ! Ice cloud fraction                           [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      precip_frac,   & ! Precipitation fraction (overall)               [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)     [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)     [-]

    ! Local Variables
    real( kind = core_rknd ), dimension(nz) :: &
      weighted_pfrac1    ! Product of mixt_frac and precip_frac_1       [-]

    real( kind = core_rknd ) :: &
      r_tot_hm_1, & ! Mean total hydromet mixing ratio (1st PDF comp.)  [kg/kg]
      r_tot_hm_2, & ! Mean total hydromet mixing ratio (2nd PDF comp.)  [kg/kg]
      N_tot_hm_1, & ! Mean total hydromet concentration (1st PDF comp.) [num/kg]
      N_tot_hm_2    ! Mean total hydromet concentration (2nd PDF comp.) [num/kg]

    real( kind = core_rknd ), parameter :: &
      precip_frac_tol = cloud_frac_min  ! Minimum precip. frac.         [-]
    
    ! "Maximum allowable" hydrometeor mixing ratio in-precip component mean.
    real( kind = core_rknd ), parameter :: &
      max_hm_ip_comp_mean = 0.0025_core_rknd  ! [kg/kg]

    integer, parameter :: &
      precip_frac_calc_type = 2  ! Option used to calc. component precip_frac

    integer :: &
      k, i   ! Loop indices


    ! Initialize the precipitation fraction variables (precip_frac,
    ! precip_frac_1, and precip_frac_2) to 0.
    precip_frac   = zero
    precip_frac_1 = zero
    precip_frac_2 = zero

    !!! Find overall precipitation fraction.
    do k = nz, 1, -1

       ! The precipitation fraction is the greatest cloud fraction at or above a
       ! vertical level.
       if ( k < nz ) then
          precip_frac(k) = max( precip_frac(k+1), cloud_frac(k) )
       else  ! k = nz
          precip_frac(k) = cloud_frac(k)
       endif

       if ( any( hydromet(k,:) > hydromet_tol(:) ) &
            .and. precip_frac(k) < precip_frac_tol ) then

          ! In a scenario where we find any hydrometeor at this grid level, but
          ! no cloud at or above this grid level, set precipitation fraction to
          ! a minimum threshold value.
          precip_frac(k) = precip_frac_tol

       elseif ( all( hydromet(k,:) <= hydromet_tol(:) ) &
                .and. precip_frac(k) < precip_frac_tol ) then

          ! The means (overall) of every precipitating hydrometeor are all less
          ! than their respective tolerance amounts.  They are all considered to
          ! have values of 0.  There are not any hydrometeor species found at
          ! this grid level.  There is also no cloud at or above this grid
          ! level, so set precipitation fraction to 0.
          precip_frac(k) = zero

       endif

    enddo ! Overall precipitation fraction loop: k = nz, 1, -1.

    !!! Account for ice cloud fraction
    do k = nz, 1, -1
      precip_frac(k) = max( precip_frac(k), ice_supersat_frac(k) )
    enddo


    !!! Find precipitation fraction within each PDF component.
    !
    ! The overall precipitation fraction, f_p, is given by the equation:
    !
    ! f_p = a * f_p(1) + ( 1 - a ) * f_p(2);
    !
    ! where a is the mixture fraction (weight of PDF component 1), f_p(1) is
    ! the precipitation fraction within PDF component 1, and f_p(2) is the
    ! precipitation fraction within PDF component 2.  Overall precipitation
    ! fraction is found according the method above, and mixture fraction is
    ! already determined, leaving f_p(1) and f_p(2) to be solved for.  The
    ! values for f_p(1) and f_p(2) must satisfy the above equation.
    if ( precip_frac_calc_type == 1 ) then

       ! This method needs some improvements -- Brian; 11/10/2014.

       !!! Find precipitation fraction within PDF component 1.
       ! The method used to find overall precipitation fraction will also be to
       ! find precipitation fraction within PDF component 1.  In order to do so,
       ! it is assumed (poorly) that PDF component 1 overlaps PDF component 1 at
       ! every vertical level in the vertical profile.
       do k = nz, 1, -1

          ! The weighted precipitation fraction (PDF component 1) is the
          ! greatest value of the product of mixture fraction and cloud fraction
          ! (PDF component 1) at or above a vertical level.
          if ( k < nz ) then
             weighted_pfrac1(k) = max( weighted_pfrac1(k+1), &
                                       mixt_frac(k) * cloud_frac_1(k) )
          else  ! k = nz
             weighted_pfrac1(k) = mixt_frac(k) * cloud_frac_1(k)
          endif

          precip_frac_1(k) = weighted_pfrac1(k) / mixt_frac(k)

          ! Special cases for precip_frac_1.
          if ( precip_frac_1(k) > one ) then

             ! Using the above method, it is possible for precip_frac_1 to be
             ! greater than 1.  For example, the mixture fraction at level k+1
             ! is 0.10 and the cloud_frac_1 at level k+1 is 1, resulting in a
             ! weighted_pfrac1 of 0.10.  This product is greater than the
             ! product of mixt_frac and cloud_frac_1 at level k.  The mixture
             ! fraction at level k is 0.05, resulting in a precip_frac_1 of 2.
             ! The value of precip_frac_1 is limited at 1.  The leftover
             ! precipitation fraction (a result of the decreasing weight of PDF
             ! component 1 between the levels) is applied to PDF component 2.
             precip_frac_1(k) = one

          elseif ( any( hm1(k,:) > hydromet_tol(:) ) &
                   .and. precip_frac_1(k) <= precip_frac_tol ) then

             ! In a scenario where we find any hydrometeor in the 1st PDF
             ! component at this grid level, but no cloud in the 1st PDF
             ! component at or above this grid level, set precipitation fraction
             ! (in the 1st PDF component) to a minimum threshold value.
             precip_frac_1(k) = precip_frac_tol

          elseif ( all( hm1(k,:) <= hydromet_tol(:) ) &
                   .and. precip_frac_1(k) <= precip_frac_tol ) then

             ! The means of every precipitating hydrometeor in the 1st PDF
             ! component are all less than their respective tolerance amounts.
             ! They are all considered to have values of 0.  There are not any
             ! hydrometeor species found in the 1st PDF component at this grid
             ! level.  There is also no cloud at or above this grid level, so
             ! set precipitation fraction (in the 1st PDF component) to 0.
             precip_frac_1(k) = zero

          endif

       enddo ! Precipitation fraction (1st PDF component) loop: k = nz, 1, -1.


       !!! Find precipitation fraction within PDF component 2.
       ! The equation for precipitation fraction within PDF component 2 is:
       !
       ! f_p(2) = ( f_p - a * f_p(1) ) / ( 1 - a );
       !
       ! given the overall precipitation fraction, f_p (calculated above), the
       ! precipitation fraction within PDF component 1, f_p(1) (calculated
       ! above), and mixture fraction, a.  Any leftover precipitation fraction
       ! from precip_frac_1 will be included in this calculation of
       ! precip_frac_2.
       do k = 1, nz, 1

          precip_frac_2(k) &
          = ( precip_frac(k) - mixt_frac(k) * precip_frac_1(k) ) &
            / ( one - mixt_frac(k) )

          ! Special cases for precip_frac_2.
          if ( precip_frac_2(k) > one ) then

             ! Again, it is possible for precip_frac_2 to be greater than 1.
             ! For example, the mixture fraction at level k+1 is 0.10 and the
             ! cloud_frac_1 at level k+1 is 1, resulting in a weighted_pfrac1 of
             ! 0.10.  This product is greater than the product of mixt_frac and
             ! cloud_frac_1 at level k.  Additionally, precip_frac (overall) is 1
             ! for level k.  The mixture fraction at level k is 0.5, resulting
             ! in a precip_frac_1 of 0.2.  Using the above equation,
             ! precip_frac_2 is calculated to be 1.8.  The value of
             ! precip_frac_2 is limited at 1.  The leftover precipitation
             ! fraction (as a result of the increasing weight of component 1
             ! between the levels) is applied to PDF component 1.
             precip_frac_2(k) = one

             ! Recalculate the precipitation fraction in PDF component 1.
             precip_frac_1(k) &
             = ( precip_frac(k) - ( one - mixt_frac(k) ) * precip_frac_2(k) ) &
               / mixt_frac(k)

             ! Double check for errors in PDF component 1.
             if ( precip_frac_1(k) > one ) then
                precip_frac_1(k) = one
             elseif ( any( hm1(k,:) > hydromet_tol(:) ) &
                      .and. precip_frac_1(k) <= precip_frac_tol ) then
                precip_frac_1(k) = precip_frac_tol
             elseif ( all( hm1(k,:) <= hydromet_tol(:) ) &
                      .and. precip_frac_1(k) <= precip_frac_tol ) then
                precip_frac_1(k) = zero
             endif

          elseif ( any( hm2(k,:) > hydromet_tol(:) ) &
                   .and. precip_frac_2(k) <= precip_frac_tol ) then

             ! In a scenario where we find any hydrometeor in the 2nd PDF
             ! component at this grid level, but no cloud in the 2nd PDF
             ! component at or above this grid level, set precipitation fraction
             ! (in the 2nd PDF component) to a minimum threshold value.
             precip_frac_2(k) = precip_frac_tol

          elseif ( all( hm2(k,:) <= hydromet_tol(:) ) &
                   .and. precip_frac_2(k) <= precip_frac_tol ) then

             ! The means of every precipitating hydrometeor in the 2nd PDF
             ! component are all less than their respective tolerance amounts.
             ! They are all considered to have values of 0.  There are not any
             ! hydrometeor species found in the 2nd PDF component at this grid
             ! level.  There is also no cloud at or above this grid level, so
             ! set precipitation fraction (in the 2nd PDF component) to 0.
             precip_frac_2(k) = zero

          endif

       enddo ! Precipitation fraction (2nd PDF component) loop: k = 1, nz, 1.


    elseif ( precip_frac_calc_type == 2 ) then

       ! This method needs to be eliminated.  I will keep it in the code as a
       ! temporary stopgap, since it is currently enabled by default.
       ! Brian; 11/10/2014.

       ! Precipitation fraction in each PDF component is based on the mean total
       ! hydrometeor mixing ratio in each PDF component, where total hydrometeor
       ! mixing ratio, r_Thm, is the sum of all precipitating hydrometeor
       ! species mixing ratios (which doesn't include cloud water), such that:
       !
       ! r_Thm = r_r + r_i + r_s + r_g;
       !
       ! where r_r is rain water mixing ratio, r_i is ice mixing ratio, r_s is
       ! snow mixing ratio, and r_g is graupel mixing ratio.
       !
       ! Precipitation fraction in each PDF component is based on the ratio:
       !
       ! r_Thm_1/f_p(1) = r_Thm_2/f_p(2);
       !
       ! where r_Thm_1 is mean total hydrometeor mixing ratio is the 1st PDF
       ! component and r_Thm_2 is mean total hydrometeor mixing ratio in the 2nd
       ! PDF component.  The equation can be rewritten as:
       !
       ! f_p(2)/f_p(1) = r_Thm_2/r_Thm_1.
       !
       ! Since overall precipitation fraction is given by the equation:
       !
       ! f_p = a f_p(1) + (1-a) f_p(2);
       !
       ! it can be rewritten as:
       !
       ! f_p = f_p(1) ( a + (1-a) f_p(2)/f_p(1) ).
       !
       ! Substituting the ratio r_Thm_2/r_Thm_1 for the ratio f_p(2)/f_p(1), the
       ! above equation can be solved for f_p(1):
       !
       ! f_p(1) = f_p / ( a + (1-a) r_Thm_2/r_Thm_1 ).
       !
       ! Then, f_p(2) can be solved for according to the equation:
       !
       ! f_p(2) = ( f_p - a f_p(1) ) / (1-a).
       !
       ! In the event where hydrometeor concentrations are found at a given
       ! vertical level, but not hydrometeor mixing ratios (due to numerical
       ! artifacts), the mean total hydrometeor concentrations in each PDF
       ! component will be used in place of mean total hydrometeor mixing ratios
       ! in the above equations to solve for component precipitation fractions.
       do k = 1, nz, 1

          if ( all( hm1(k,:) <= hydromet_tol(:) ) &
               .and. all( hm2(k,:) <= hydromet_tol(:) ) ) then

             ! There are no hydrometeors found in each PDF component.
             ! Precipitation fraction within each component is set to 0.
             precip_frac_1(k) = zero
             precip_frac_2(k) = zero

          elseif ( any( hm1(k,:) > hydromet_tol(:) ) &
                   .and. all( hm2(k,:) <= hydromet_tol(:) ) ) then

             ! All the hydrometeors are found within the 1st PDF component.
             precip_frac_1(k) = precip_frac(k) / mixt_frac(k)
             precip_frac_2(k) = zero

             ! Using the above method, it is possible for precip_frac_1 to be
             ! greater than 1. The value of precip_frac_1 is limited at 1.
             if ( precip_frac_1(k) > one ) then
                precip_frac_1(k) = one
                precip_frac(k) = mixt_frac(k)
             endif

          elseif ( any( hm2(k,:) > hydromet_tol(:) ) &
                   .and. all( hm1(k,:) <= hydromet_tol(:) ) ) then

             ! All the hydrometeors are found within the 2nd PDF component.
             precip_frac_1(k) = zero
             precip_frac_2(k) = precip_frac(k) / ( one - mixt_frac(k) )

             ! Using the above method, it is possible for precip_frac_2 to be
             ! greater than 1. The value of precip_frac_2 is limited at 1.
             if ( precip_frac_2(k) > one ) then
                precip_frac_2(k) = one
                precip_frac(k) = one - mixt_frac(k)
             endif

          else

             ! any( hm1(k,:) > hydromet_tol(:) )
             ! AND any( hm2(k,:) > hydromet_tol(:) )

             ! Hydrometeors are found within both PDF components.
             r_tot_hm_1 = zero
             r_tot_hm_2 = zero
             N_tot_hm_1 = zero
             N_tot_hm_2 = zero
             do i = 1, hydromet_dim, 1

                if ( l_mix_rat_hm(i) ) then

                   ! The hydrometeor is a mixing ratio.
                   ! Find total hydrometeor mixing ratio in each PDF component.
                   if ( hm1(k,i) > hydromet_tol(i) ) then
                      r_tot_hm_1 = r_tot_hm_1 + hm1(k,i)
                   endif
                   if ( hm2(k,i) > hydromet_tol(i) ) then
                      r_tot_hm_2 = r_tot_hm_2 + hm2(k,i)
                   endif

                else ! l_mix_rat_hm(i) is false

                   ! The hydrometeor is a concentration.
                   ! Find total hydrometeor concentration in each PDF component.
                   if ( hm1(k,i) > hydromet_tol(i) ) then
                      N_tot_hm_1 = N_tot_hm_1 + hm1(k,i)
                   endif
                   if ( hm2(k,i) > hydromet_tol(i) ) then
                      N_tot_hm_2 = N_tot_hm_2 + hm2(k,i)
                   endif

                endif ! l_mix_rat_hm(i)

             enddo ! i = 1, hydromet_dim, 1

             !!! Find precipitation fraction within PDF component 1.
             if ( r_tot_hm_1 > zero ) then
                precip_frac_1(k) &
                = precip_frac(k) &
                  / ( mixt_frac(k) &
                      + ( one - mixt_frac(k) ) * r_tot_hm_2/r_tot_hm_1 )
             else ! N_tot_hm_1 > zero 
                precip_frac_1(k) &
                = precip_frac(k) &
                  / ( mixt_frac(k) &
                      + ( one - mixt_frac(k) ) * N_tot_hm_2/N_tot_hm_1 )
             endif

             ! Using the above method, it is possible for precip_frac_1 to be
             ! greater than 1.  The value of precip_frac_1 is limited at 1.
             if ( precip_frac_1(k) > one ) then
                precip_frac_1(k) = one
             endif

             !!! Find precipitation fraction within PDF component 2.
             precip_frac_2(k) &
             = ( precip_frac(k) - mixt_frac(k) *  precip_frac_1(k) ) &
               / ( one - mixt_frac(k) )

             ! Using the above method, it is possible for precip_frac_2 to be
             ! greater than 1.  The value of precip_frac_2 is limited at 1.
             if ( precip_frac_2(k) > one ) then

                precip_frac_2(k) = one

                ! Recalculate the precipitation fraction in PDF component 1.
                precip_frac_1(k) &
                = ( precip_frac(k) &
                    - ( one - mixt_frac(k) ) * precip_frac_2(k) ) &
                  / mixt_frac(k)

             endif

          endif


          ! Special cases for PDF component 1.
          if ( any( hm1(k,:) > hydromet_tol(:) ) &
               .and. precip_frac_1(k) <= precip_frac_tol ) then

             ! In a scenario where we find any hydrometeor in the 1st PDF
             ! component at this grid level, but no cloud in the 1st PDF
             ! component at or above this grid level, set precipitation fraction
             ! (in the 1st PDF component) to a minimum threshold value.
             precip_frac_1(k) = precip_frac_tol

          elseif ( all( hm1(k,:) <= hydromet_tol(:) ) &
                   .and. precip_frac_1(k) <= precip_frac_tol ) then

             ! The means of every precipitating hydrometeor in the 1st PDF
             ! component are all less than their respective tolerance amounts.
             ! They are all considered to have values of 0.  There is not any
             ! hydrometeor species found in the 1st PDF component at this grid
             ! level.  There is also no cloud at or above this grid level, so
             ! set precipitation fraction (in the 1st PDF component) to 0.
             precip_frac_1(k) = zero

          endif


          ! Special cases for PDF component 2.
          if ( any( hm2(k,:) > hydromet_tol(:) ) &
               .and. precip_frac_2(k) <= precip_frac_tol ) then

             ! In a scenario where we find any hydrometeor in the 2nd PDF
             ! component at this grid level, but no cloud in the 2nd PDF
             ! component at or above this grid level, set precipitation fraction
             ! (in the 2nd PDF component) to a minimum threshold value.
             precip_frac_2(k) = precip_frac_tol

          elseif ( all( hm2(k,:) <= hydromet_tol(:) ) &
                   .and. precip_frac_2(k) <= precip_frac_tol ) then

             ! The means of every precipitating hydrometeor in the 2nd PDF
             ! component are all less than their respective tolerance amounts.
             ! They are all considered to have values of 0.  There is not any
             ! hydrometeor species found in the 2nd PDF component at this grid
             ! level.  There is also no cloud at or above this grid level, so
             ! set precipitation fraction (in the 2nd PDF component) to 0.
             precip_frac_2(k) = zero

          endif


       enddo ! Component precipitation fraction loop: k = 1, nz, 1.


    elseif ( precip_frac_calc_type == 3 ) then

      ! Temporary third option to test setting precip_frac_1 = precip_frac_2
      ! ( = precip_frac ).  Brian; 11/10/2014.
      precip_frac_1 = precip_frac
      precip_frac_2 = precip_frac


    else ! Invalid option selected.

       write(fstderr,*) "Invalid option to calculate precip_frac_1 " &
                        // "and precip_frac_2."
       stop


    endif ! precip_frac_calc_type


    ! Increase Precipiation Fraction under special conditions.
    !
    ! There are scenarios that sometimes occur that require precipitation
    ! fraction to be boosted.  Precipitation fraction is calculated from cloud
    ! fraction and ice supersaturation fraction.  For numerical reasons, CLUBB's
    ! PDF may become entirely subsaturated with respect to liquid and ice,
    ! resulting in both a cloud fraction of 0 and an ice supersaturation
    ! fraction of 0.  When this happens, precipitation fraction drops to 0 when
    ! there aren't any hydrometeors present at that grid level, or to
    ! precip_frac_tol when there is at least one hydrometeor present at that
    ! grid level.  However, sometimes there are large values of hydrometeors
    ! found at that grid level.  When this occurs, the PDF component in-precip
    ! mean of a hydrometeor can become ridiculously large.  This is because the
    ! ith PDF component in-precip mean of a hydrometeor, mu_hm_i,  is given by
    ! the equation:
    !
    ! mu_hm_i = hmi / precip_frac_i;
    !
    ! where hmi is the overall ith PDF component mean of the hydrometeor, and
    ! precip_frac_i is the ith PDF component precipitation fraction.  When
    ! precip_frac_i has a value of precip_frac_tol and hmi is large, mu_hm_i can
    ! be huge.  This can cause enormous microphysical process rates and result
    ! in numerical instability.  It is also very inaccurate.
    !
    ! In order to limit this problem, the ith PDF component precipitation
    ! fraction is increased in order to decrease mu_hm_i.  First, an "upper
    ! limit" is set for mu_hm_i when the hydrometeor is a mixing ratio.  This is
    ! called max_hm_ip_comp_mean.  At every vertical level and for every
    ! hydrometeor mixing ratio, a check is made to try to prevent mu_hm_i from
    ! exceeding the "upper limit".  The check is:
    ! hmi / precip_frac_i ( which = mu_hm_i ) > max_hm_ip_comp_mean, which can
    ! be rewritten:  hmi > precip_frac_i * max_hm_ip_comp_mean.  When this
    ! occurs, precip_frac_i is increased to hmi/max_hm_ip_comp_mean.  Of course,
    ! precip_frac_i is not allowed to exceed 1, so when hmi is already greater
    ! than max_hm_ip_comp_mean, mu_hm_i will also have to be greater than
    ! max_hm_ip_comp_mean.  However, the value of mu_hm_i is still reduced when
    ! compared to what it would have been using precip_frac_tol.  In the event
    ! that multiple hydrometeor mixing ratios violate the check, the code is set
    ! up so that precip_frac_i is increased based on the highest hmi.
    do k = 1, nz, 1

       do i = 1, hydromet_dim, 1

          if ( l_mix_rat_hm(i) ) then

             ! The hydrometeor is a mixing ratio.

             if ( hm1(k,i) > precip_frac_1(k) * max_hm_ip_comp_mean ) then

                ! Increase precipitation fraction in the 1st PDF component.
                precip_frac_1(k) = min( hm1(k,i)/max_hm_ip_comp_mean, one )

                ! Recalculate overall precipitation fraction.
                precip_frac(k) = mixt_frac(k) * precip_frac_1(k) &
                                 + ( one - mixt_frac(k) ) * precip_frac_2(k)

             endif ! mu_hm_1 = hm1/precip_frac_1 > max_hm_ip_comp_mean

             if ( hm2(k,i) > precip_frac_2(k) * max_hm_ip_comp_mean ) then

                ! Increase precipitation fraction in the 2nd PDF component.
                precip_frac_2(k) = min( hm2(k,i)/max_hm_ip_comp_mean, one )

                ! Recalculate overall precipitation fraction.
                precip_frac(k) = mixt_frac(k) * precip_frac_1(k) &
                                 + ( one - mixt_frac(k) ) * precip_frac_2(k)

             endif ! mu_hm_2 = hm2/precip_frac_2 > max_hm_ip_comp_mean

          endif ! l_mix_rat_hm(i)

       enddo ! i = 1, hydromet_dim, 1

    enddo ! k = 1, nz, 1


    return

  end subroutine precip_fraction

  !=============================================================================
  subroutine compute_mean_stdev( Ncnm, rc_1, rc_2, &                      ! Intent(in)
                                 cloud_frac_1, cloud_frac_2, &          ! Intent(in)
                                 hm1, hm2, &                            ! Intent(in)
                                 precip_frac_1, precip_frac_2, &        ! Intent(in)
                                 sigma2_on_mu2_ip_array_cloud, &        ! Intent(in)
                                 sigma2_on_mu2_ip_array_below, &        ! Intent(in)
                                 pdf_params, d_variables, &             ! Intent(in)
                                 mu_x_1, mu_x_2, sigma_x_1, sigma_x_2 ) ! Intent(out)
       
    ! Description:
    ! Calculates the means and standard deviations (for each PDF component) of
    ! chi, eta, w, Ncn, and the precipitating hydrometeors.  For the precipitating
    ! hydrometeors, the component means and standard deviations are in-precip. 

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

    use corr_varnce_module, only: &
        iiPDF_chi, & ! Variable(s)
        iiPDF_eta, &
        iiPDF_w,        &
        iiPDF_Ncn

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: d_variables ! Number of PDF variables

    real( kind = core_rknd ), intent(in) :: &
      Ncnm,          & ! Mean cloud nuclei concentration                [num/kg]
      rc_1,          & ! Mean of r_c (1st PDF component)                 [kg/kg]
      rc_2,          & ! Mean of r_c (2nd PDF component)                 [kg/kg]
      cloud_frac_1,  & ! Cloud fraction (1st PDF component)                  [-]
      cloud_frac_2,  & ! Cloud fraction (2nd PDF component)                  [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)          [-]

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      sigma2_on_mu2_ip_array_cloud, & ! Prescribed ratio array: cloudy levs. [-]
      sigma2_on_mu2_ip_array_below    ! Prescribed ratio array: clear levs.  [-]

    real( kind = core_rknd ), dimension(hydromet_dim), intent(in) :: &
      hm1, & ! Mean of a precip. hydrometeor (1st PDF component)    [units vary]
      hm2    ! Mean of a precip. hydrometeor (2nd PDF component)    [units vary]

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                                [units vary]

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

    ! Local Variables
    integer :: ivar ! Loop iterator


    !!! Enter the PDF parameters.

    !!! Means.

    ! Mean of vertical velocity, w, in PDF component 1.
    mu_x_1(iiPDF_w) = pdf_params%w_1

    ! Mean of vertical velocity, w, in PDF component 2.
    mu_x_2(iiPDF_w) = pdf_params%w_2

    ! Mean of extended liquid water mixing ratio, chi (old s),
    ! in PDF component 1.
    mu_x_1(iiPDF_chi) = pdf_params%chi_1

    ! Mean of extended liquid water mixing ratio, chi (old s),
    ! in PDF component 2.
    mu_x_2(iiPDF_chi) = pdf_params%chi_2

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

    ! Mean of simplified cloud nuclei concentration, Ncn, in PDF component 1.
    mu_x_1(iiPDF_Ncn) = Ncnm

    ! Mean of simplified cloud nuclei concentration, Ncn, in PDF component 2.
    mu_x_2(iiPDF_Ncn) = Ncnm

    ! Mean of the hydrometeor species
    do ivar = iiPDF_Ncn+1, d_variables

       ! Mean of hydrometeor, hm, in PDF component 1.
       mu_x_1(ivar) &
       = component_mean_hm_ip( hm1(pdf2hydromet_idx(ivar)), precip_frac_1, &
                               hydromet_tol(pdf2hydromet_idx(ivar)) )

       ! Mean of hydrometeor, hm, in PDF component 2.
       mu_x_2(ivar) &
       = component_mean_hm_ip( hm2(pdf2hydromet_idx(ivar)), precip_frac_2, &
                               hydromet_tol(pdf2hydromet_idx(ivar)) )

    enddo


    !!! Standard deviations.

    ! Standard deviation of vertical velocity, w, in PDF component 1.
    sigma_x_1(iiPDF_w) = sqrt( pdf_params%varnce_w_1 )

    ! Standard deviation of vertical velocity, w, in PDF component 2.
    sigma_x_2(iiPDF_w) = sqrt( pdf_params%varnce_w_2 )

    ! Standard deviation of extended liquid water mixing ratio, chi (old s),
    ! in PDF component 1.
    sigma_x_1(iiPDF_chi) = pdf_params%stdev_chi_1

    ! Standard deviation of extended liquid water mixing ratio, chi (old s),
    ! in PDF component 2.
    sigma_x_2(iiPDF_chi) = pdf_params%stdev_chi_2

    ! Standard deviation of eta (old t) in PDF component 1.
    sigma_x_1(iiPDF_eta) = pdf_params%stdev_eta_1

    ! Standard deviation of eta (old t) in PDF component 2.
    sigma_x_2(iiPDF_eta) = pdf_params%stdev_eta_2

    ! Standard deviation of simplified cloud nuclei concentration, Ncn,
    ! in PDF component 1.
    if ( .not. l_const_Nc_in_cloud ) then

       ! Ncn varies in both PDF components.
       sigma_x_1(iiPDF_Ncn) &
       = component_stdev_hm_ip( mu_x_1(iiPDF_Ncn), rc_1, one, &
                                sigma2_on_mu2_ip_array_cloud(iiPDF_Ncn), &
                                sigma2_on_mu2_ip_array_below(iiPDF_Ncn) )

    else ! l_const_Nc_in_cloud

       ! Ncn is constant in both PDF components.
       sigma_x_1(iiPDF_Ncn) = zero

    endif ! .not. l_const_Nc_in_cloud

    ! Standard deviation of simplified cloud nuclei concentration, Ncn,
    ! in PDF component 2.
    if ( .not. l_const_Nc_in_cloud ) then

       ! Ncn varies in both PDF components.
       sigma_x_2(iiPDF_Ncn) &
       = component_stdev_hm_ip( mu_x_2(iiPDF_Ncn), rc_2, one, &
                                sigma2_on_mu2_ip_array_cloud(iiPDF_Ncn), &
                                sigma2_on_mu2_ip_array_cloud(iiPDF_Ncn) )

    else ! l_const_Nc_in_cloud

       ! Ncn is constant in both PDF components.
       sigma_x_2(iiPDF_Ncn) = zero

    endif ! .not. l_const_Nc_in_cloud

    ! Set up the values of the statistical correlations and variances.  Since we
    ! currently do not have enough variables to compute the correlations and
    ! variances directly, we have obtained these values by analyzing LES runs of
    ! certain cases.  We have divided those results into an inside-cloud average
    ! and an outside-cloud (or below-cloud) average.  This coding leaves the
    ! software architecture in place in case we ever have the variables in place
    ! to compute these values directly.  It also allows us to use separate
    ! inside-cloud and outside-cloud parameter values.
    ! Brian Griffin; February 3, 2007.

    do ivar = iiPDF_Ncn+1, d_variables

       ! Standard deviation of hydrometeor, hm, in PDF component 1.
       sigma_x_1(ivar) &
       =  component_stdev_hm_ip( mu_x_1(ivar), &
                                 rc_1, cloud_frac_1, &
                                 sigma2_on_mu2_ip_array_cloud(ivar), &
                                 sigma2_on_mu2_ip_array_below(ivar) )

       ! Standard deviation of hydrometeor, hm, in PDF component 2.
       sigma_x_2(ivar) &
       =  component_stdev_hm_ip( mu_x_2(ivar), &
                                 rc_2, cloud_frac_2, &
                                 sigma2_on_mu2_ip_array_cloud(ivar), &
                                 sigma2_on_mu2_ip_array_below(ivar) )

    enddo


    return

  end subroutine compute_mean_stdev

  !=============================================================================
  subroutine compute_corr( wm_zt, rc_1, rc_2, cloud_frac_1, &
                           cloud_frac_2, wpchip, wpNcnp, &
                           stdev_w, mixt_frac, precip_frac_1, &
                           precip_frac_2, wphydrometp_zt, &
                           mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &
                           corr_array_cloud, corr_array_below, &
                           pdf_params, d_variables, &
                           corr_array_1, corr_array_2 )

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
      mu_x_1,    & ! Mean of x array (1st PDF component)            [units vary]
      mu_x_2,    & ! Mean of x array (2nd PDF component)            [units vary]
      sigma_x_1, & ! Standard deviation of x array (1st PDF comp.)  [units vary]
      sigma_x_2    ! Standard deviation of x array (2nd PDF comp.)  [units vary]

    real( kind = core_rknd ), dimension(d_variables, d_variables), &
    intent(in) :: &
      corr_array_cloud, & ! Prescribed correlation array in cloud        [-]
      corr_array_below    ! Prescribed correlation array below cloud     [-]

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                                [units vary]

    ! Output Variables
    real( kind = core_rknd ), dimension(d_variables, d_variables), &
    intent(out) :: &
      corr_array_1, & ! Correlation array (1st PDF component) [-]
      corr_array_2    ! Correlation array (2nd PDF component) [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      sigma_Ncn_1

    real( kind = core_rknd ), dimension(d_variables)  :: &
      corr_w_hm_1, & ! Correlation of w and hm (1st PDF component) ip    [-]
      corr_w_hm_2    ! Correlation of w and hm (2nd PDF component) ip    [-]

    real( kind = core_rknd ) :: &
      chi_m,      & ! Mean of chi (s_mellor)                    [kg/kg]
      stdev_chi,  & ! Standard deviation of chi (s_mellor)      [kg/kg]
      corr_w_chi, & ! Correlation of w and chi (overall)        [-]
      corr_w_Ncn    ! Correlation of w and Ncn (overall)        [-]

    logical :: &
      l_limit_corr_chi_eta    ! Flag to limit the correlation of chi and eta [-]

    integer :: ivar, jvar ! Loop iterators

    ! ---- Begin Code ----

    !!! Enter the PDF parameters.
    sigma_Ncn_1 = sigma_x_1(iiPDF_Ncn)

    !!! Correlations

    ! Initialize corr_w_hm_1 and corr_w_hm_2 arrays to 0.
    corr_w_hm_1 = zero
    corr_w_hm_2 = zero

    ! Calculate correlations involving w by first calculating total covariances
    ! involving w (<w'r_r'>, etc.) using the down-gradient approximation.
    if ( l_calc_w_corr ) then

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

       corr_w_Ncn = calc_w_corr( wpNcnp, stdev_w, sigma_Ncn_1, w_tol, Ncn_tol )

       do jvar = iiPDF_Ncn+1, d_variables

          call calc_corr_w_hm( wm_zt, wphydrometp_zt(pdf2hydromet_idx(jvar)), &
                               mu_x_1(iiPDF_w), mu_x_2(iiPDF_w), &
                               mu_x_1(jvar), mu_x_2(jvar), &
                               sigma_x_1(iiPDF_w), sigma_x_2(iiPDF_w), &
                               sigma_x_1(jvar), sigma_x_2(jvar), &
                               mixt_frac, precip_frac_1, precip_frac_2, &
                               corr_w_hm_1(jvar), corr_w_hm_2(jvar), &
                               hydromet_tol(pdf2hydromet_idx(jvar)) )

       enddo ! jvar = iiPDF_Ncn+1, d_variables


    endif

    ! In order to decompose the correlation matrix,
    ! we must not have a perfect correlation of chi and
    ! eta. Thus, we impose a limitation.
    l_limit_corr_chi_eta = .true.


    ! Initialize the correlation arrays
    corr_array_1 = zero
    corr_array_2 = zero

    !!! The corr_arrays are assumed to be lower triangular matrices
    ! Set diagonal elements to 1
    do ivar=1, d_variables
      corr_array_1(ivar, ivar) = one
      corr_array_2(ivar, ivar) = one
    end do


    !!! This code assumes the following order in the prescribed correlation
    !!! arrays (iiPDF indices):
    !!! chi, eta, w, Ncn, <hydrometeors> (indices increasing from left to right)

    ! Correlation of chi (old s) and eta (old t)
    corr_array_1(iiPDF_eta, iiPDF_chi) &
    = component_corr_chi_eta( pdf_params%corr_chi_eta_1, rc_1, cloud_frac_1, &
                              corr_array_cloud(iiPDF_eta, iiPDF_chi), &
                              corr_array_below(iiPDF_eta, iiPDF_chi), &
                              l_limit_corr_chi_eta )

    corr_array_2(iiPDF_eta, iiPDF_chi) &
    = component_corr_chi_eta( pdf_params%corr_chi_eta_2, rc_2, cloud_frac_2, &
                              corr_array_cloud(iiPDF_eta, iiPDF_chi), &
                              corr_array_below(iiPDF_eta, iiPDF_chi), &
                              l_limit_corr_chi_eta )

    ! Correlation of chi (old s) and w
    corr_array_1(iiPDF_w, iiPDF_chi) &
    = component_corr_w_x( corr_w_chi, rc_1, cloud_frac_1, &
                          corr_array_cloud(iiPDF_w, iiPDF_chi), &
                          corr_array_below(iiPDF_w, iiPDF_chi) )

    corr_array_2(iiPDF_w, iiPDF_chi) &
    = component_corr_w_x( corr_w_chi, rc_2, cloud_frac_2, &
                          corr_array_cloud(iiPDF_w, iiPDF_chi), &
                          corr_array_below(iiPDF_w, iiPDF_chi) )


    ! Correlation of chi (old s) and Ncn
    corr_array_1(iiPDF_Ncn, iiPDF_chi) &
    = component_corr_x_hm_ip( rc_1, one, &
                              corr_array_cloud(iiPDF_Ncn, iiPDF_chi), &
                              corr_array_cloud(iiPDF_Ncn, iiPDF_chi) )

    corr_array_2(iiPDF_Ncn, iiPDF_chi) &
    = component_corr_x_hm_ip( rc_2, one, &
                              corr_array_cloud(iiPDF_Ncn, iiPDF_chi), &
                              corr_array_cloud(iiPDF_Ncn, iiPDF_chi) )

    ! Correlation of chi (old s) and the hydrometeors
    ivar = iiPDF_chi
    do jvar = iiPDF_Ncn+1, d_variables
       corr_array_1(jvar, ivar) &
       = component_corr_x_hm_ip( rc_1, cloud_frac_1,&
                                 corr_array_cloud(jvar, ivar), &
                                 corr_array_below(jvar, ivar) )

       corr_array_2(jvar, ivar) &
       = component_corr_x_hm_ip( rc_2, cloud_frac_2,&
                                 corr_array_cloud(jvar, ivar), &
                                 corr_array_below(jvar, ivar) )
    enddo

    ! Correlation of eta (old t) and w
    corr_array_1(iiPDF_w, iiPDF_eta) = zero
    corr_array_2(iiPDF_w, iiPDF_eta) = zero

    ! Correlation of eta (old t) and Ncn
    corr_array_1(iiPDF_Ncn, iiPDF_eta) &
    = component_corr_x_hm_ip( rc_1, one, &
                              corr_array_cloud(iiPDF_Ncn, iiPDF_eta), &
                              corr_array_cloud(iiPDF_Ncn, iiPDF_eta) )

    corr_array_2(iiPDF_Ncn, iiPDF_eta) &
    = component_corr_x_hm_ip( rc_2, one, &
                              corr_array_cloud(iiPDF_Ncn, iiPDF_eta), &
                              corr_array_cloud(iiPDF_Ncn, iiPDF_eta) )

    ! Correlation of eta (old t) and the hydrometeors
    ivar = iiPDF_eta
    do jvar = iiPDF_Ncn+1, d_variables
      corr_array_1(jvar, ivar) &
      = component_corr_eta_hm_ip( corr_array_1( iiPDF_eta, iiPDF_chi), &
                                  corr_array_1( jvar, iiPDF_chi) )

      corr_array_2(jvar, ivar) &
      = component_corr_eta_hm_ip( corr_array_2( iiPDF_eta, iiPDF_chi), &
                                  corr_array_2( jvar, iiPDF_chi) )
    enddo


    ! Correlation of w and Ncn
    corr_array_1(iiPDF_Ncn, iiPDF_w) &
    = component_corr_w_hm_ip( corr_w_Ncn, rc_1, one, &
                              corr_array_cloud(iiPDF_Ncn, iiPDF_w), &
                              corr_array_below(iiPDF_Ncn, iiPDF_w) )

    corr_array_2(iiPDF_Ncn, iiPDF_w) &
    = component_corr_w_hm_ip( corr_w_Ncn, rc_2, one, &
                              corr_array_cloud(iiPDF_Ncn, iiPDF_w), &
                              corr_array_below(iiPDF_Ncn, iiPDF_w) )

    ! Correlation of w and the hydrometeors
    ivar = iiPDF_w
    do jvar = iiPDF_Ncn+1, d_variables

       corr_array_1(jvar, ivar) &
       = component_corr_w_hm_ip( corr_w_hm_1(jvar), rc_1, cloud_frac_1, &
                                 corr_array_cloud(jvar, ivar), &
                                 corr_array_below(jvar, ivar) )

       corr_array_2(jvar, ivar) &
       = component_corr_w_hm_ip( corr_w_hm_2(jvar), rc_2, cloud_frac_2, &
                                 corr_array_cloud(jvar, ivar), &
                                 corr_array_below(jvar, ivar) )

    enddo

    ! Correlation of Ncn and the hydrometeors
    ivar = iiPDF_Ncn
    do jvar = iiPDF_Ncn+1, d_variables
       corr_array_1(jvar, ivar) &
       = component_corr_hmx_hmy_ip( rc_1, cloud_frac_1, &
                                    corr_array_cloud(jvar, ivar), &
                                    corr_array_below(jvar, ivar) )

       corr_array_2(jvar, ivar) &
       = component_corr_hmx_hmy_ip( rc_2, cloud_frac_2, &
                                    corr_array_cloud(jvar, ivar), &
                                    corr_array_below(jvar, ivar) )
    enddo

    ! Correlation of two hydrometeors
    do ivar = iiPDF_Ncn+1, d_variables-1
       do jvar = ivar+1, d_variables

          corr_array_1(jvar, ivar) &
          = component_corr_hmx_hmy_ip( rc_1, cloud_frac_1, &
                                       corr_array_cloud(jvar, ivar), &
                                       corr_array_below(jvar, ivar) )

          corr_array_2(jvar, ivar) &
          = component_corr_hmx_hmy_ip( rc_2, cloud_frac_2, &
                                       corr_array_cloud(jvar, ivar), &
                                       corr_array_below(jvar, ivar) )

       enddo ! jvar
    enddo ! ivar


    return

  end subroutine compute_corr

  !=============================================================================
  function component_mean_hm_ip( hmi, precip_frac_i, hydromet_tol )  &
  result( mu_hm_i )

    ! Description:
    ! Calculates the in-precip mean of a hydrometeor species within the ith
    ! PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        zero  ! Constant(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      hmi,           & ! Mean of hydrometeor, hm (ith PDF component) [hm units]
      precip_frac_i, & ! Precipitation fraction (ith PDF component)  [-]
      hydromet_tol     ! Tolerance value for the hydrometeor         [hm units]

    ! Return Variable
    real( kind = core_rknd ) :: &
      mu_hm_i    ! Mean of hm (ith PDF component) in-precip (ip)     [hm units]


    ! Mean of the hydrometeor (in-precip) in the ith PDF component.
    if ( hmi > hydromet_tol ) then
       mu_hm_i = hmi / precip_frac_i
    else
       ! The mean of the hydrometeor in the ith PDF component is less than the
       ! tolerance amount for the particular hydrometeor.  It is considered to
       ! have a value of 0.  There is not any of this hydrometeor species in the
       ! ith PDF component at this grid level.
       mu_hm_i = zero
    endif


    return

  end function component_mean_hm_ip

  !=============================================================================
  function component_stdev_hm_ip( mu_hm_i, rci, cloud_fraci, &
                                  hm_sigma2_on_mu2_cloud, &
                                  hm_sigma2_on_mu2_below )  &
  result( sigma_hm_i )

    ! Description:
    ! Calculates the in-precip standard deviation of a hydrometeor species
    ! within the ith PDF component.

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
      mu_hm_i,     & ! Mean of hm (ith PDF component) in-precip (ip) [hm units]
      rci,         & ! Mean cloud water mixing ratio (ith PDF comp.) [kg/kg]
      cloud_fraci    ! Cloud fraction (ith PDF component)            [-]

    real( kind = core_rknd ), intent(in) :: &
      hm_sigma2_on_mu2_cloud, & ! Ratio sigma_hm_1^2/mu_hm_1^2; cloudy levs. [-]
      hm_sigma2_on_mu2_below    ! Ratio sigma_hm_2^2/mu_hm_2^2; clear levs.  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      sigma_hm_i    ! Standard deviation of hm (ith PDF component) ip [hm units]


    ! Standard deviation of the hydrometeor (in-precip) in the
    ! ith PDF component.
    if ( l_interp_prescribed_params ) then
       sigma_hm_i = sqrt( cloud_fraci * hm_sigma2_on_mu2_cloud &
                        + ( one - cloud_fraci ) * hm_sigma2_on_mu2_below ) &
                    * mu_hm_i
    else
       if ( rci > rc_tol ) then
          sigma_hm_i = sqrt( hm_sigma2_on_mu2_cloud ) * mu_hm_i
       else
          sigma_hm_i = sqrt( hm_sigma2_on_mu2_below ) * mu_hm_i
       endif
    endif

    return

  end function component_stdev_hm_ip

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
  function component_corr_w_hm_ip( corr_w_hm_i_in, rc_i, cloud_frac_i, &
                                   corr_w_hm_NL_cloud, corr_w_hm_NL_below ) &
  result( corr_w_hm_i )

    ! Description:
    ! Calculates the in-precip correlation of w and a hydrometeor species
    ! within the ith PDF component.

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
      corr_w_hm_i_in, & ! Correlation of w and hm (ith PDF comp.) ip     [-]
      rc_i,           & ! Mean cloud water mixing ratio (ith PDF comp.)  [kg/kg]
      cloud_frac_i      ! Cloud fraction (ith PDF component)             [-]

    real( kind = core_rknd ), intent(in) :: &
      corr_w_hm_NL_cloud, & ! Corr. of w and hm (ith PDF comp.) ip; cloudy [-]
      corr_w_hm_NL_below    ! Corr. of w and hm (ith PDF comp.) ip; clear  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      corr_w_hm_i    ! Correlation of w and hm (ith PDF component) ip  [-]


    ! Correlation (in-precip) of w and the hydrometeor in the ith PDF component.
    if ( l_calc_w_corr ) then
       corr_w_hm_i = corr_w_hm_i_in
    else ! use prescribed parameter values
       if ( l_interp_prescribed_params ) then
          corr_w_hm_i = cloud_frac_i * corr_w_hm_NL_cloud &
                        + ( one - cloud_frac_i ) * corr_w_hm_NL_below
       else
          if ( rc_i > rc_tol ) then
             corr_w_hm_i = corr_w_hm_NL_cloud
          else
             corr_w_hm_i = corr_w_hm_NL_below
          endif
       endif ! l_interp_prescribed_params
    endif ! l_calc_w_corr

    return

  end function component_corr_w_hm_ip

  !=============================================================================
  function component_corr_x_hm_ip( rc_i, cloud_frac_i, &
                                   corr_x_hm_NL_cloud, corr_x_hm_NL_below ) &
  result( corr_x_hm_i )

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
      rc_i,         & ! Mean cloud water mixing ratio (ith PDF comp.) [kg/kg]
      cloud_frac_i    ! Cloud fraction (ith PDF component)            [-]

    real( kind = core_rknd ), intent(in) :: &
      corr_x_hm_NL_cloud, & ! Corr. of x and hm (ith PDF comp.) ip; cloudy [-]
      corr_x_hm_NL_below    ! Corr. of x and hm (ith PDF comp.) ip; clear  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      corr_x_hm_i    ! Correlation of x and hm (ith PDF component) ip  [-]


    ! Correlation (in-precip) of x and the hydrometeor in the ith PDF component.
    if ( l_interp_prescribed_params ) then
       corr_x_hm_i = cloud_frac_i * corr_x_hm_NL_cloud &
                     + ( one - cloud_frac_i ) * corr_x_hm_NL_below
    else
       if ( rc_i > rc_tol ) then
          corr_x_hm_i = corr_x_hm_NL_cloud
       else
          corr_x_hm_i = corr_x_hm_NL_below
       endif
    endif

    return

  end function component_corr_x_hm_ip

  !=============================================================================
  function component_corr_hmx_hmy_ip( rc_i, cloud_frac_i, &
                                      corr_hmx_hmy_LL_cloud, &
                                      corr_hmx_hmy_LL_below ) &
  result( corr_hmx_hmy_i )

    ! Description:
    ! Calculates the in-precip correlation of hydrometeor x and
    ! hydrometeor y within the ith PDF component.

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
      corr_hmx_hmy_LL_cloud, & ! Corr.: hmx & hmy (ith PDF comp.) ip; cloudy [-]
      corr_hmx_hmy_LL_below    ! Corr.: hmx & hmy (ith PDF comp.) ip; clear  [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      corr_hmx_hmy_i   ! Correlation of hmx & hmy (ith PDF component) ip [-]


    ! Correlation (in-precip) of hydrometeor x and hydrometeor y in the
    ! ith PDF component.
    if ( l_interp_prescribed_params ) then
       corr_hmx_hmy_i = cloud_frac_i * corr_hmx_hmy_LL_cloud &
                        + ( one - cloud_frac_i ) * corr_hmx_hmy_LL_below
    else
       if ( rc_i > rc_tol ) then
          corr_hmx_hmy_i = corr_hmx_hmy_LL_cloud
       else
          corr_hmx_hmy_i = corr_hmx_hmy_LL_below
       endif
    endif

    return

  end function component_corr_hmx_hmy_ip

  !=============================================================================
  pure function component_corr_eta_hm_ip( corr_chi_eta_i, corr_chi_hm_i ) &
  result( corr_eta_hm_i )

    ! Description:
    ! Estimates the correlation of eta and a hydrometeor species using the
    ! correlation of chi and eta and the correlation of chi and the hydrometeor.
    ! This facilities the Cholesky decomposability of the correlation array that
    ! will inevitably be decomposed for SILHS purposes. Without this estimation,
    ! we have found that the resulting correlation matrix cannot be decomposed.

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd       ! Constant

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      corr_chi_eta_i, & ! Component correlation of chi and eta              [-]
      corr_chi_hm_i     ! Component correlation of chi and the hydrometeor  [-]

    ! Output Variables
    real( kind = core_rknd ) :: &
      corr_eta_hm_i     ! Component correlation of eta and the hydrometeor  [-]


    corr_eta_hm_i = corr_chi_eta_i * corr_chi_hm_i


    return

  end function component_corr_eta_hm_ip

  !=============================================================================
  subroutine normalize_mean_stdev( hm1, hm2, Ncnm, d_variables, &
                                   mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &
                                   sigma2_on_mu2_ip_1, sigma2_on_mu2_ip_2, &
                                   mu_x_1_n, mu_x_2_n, &
                                   sigma_x_1_n, sigma_x_2_n )

    ! Description:
    ! Calculates the normalized means and the normalized standard deviations
    ! of PDF variables that have assumed lognormal distributions -- which are
    ! precipitating hydrometeors (in precipitation) and N_cn.

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
      hm1, & ! Mean of a precip. hydrometeor (1st PDF component)    [units vary]
      hm2    ! Mean of a precip. hydrometeor (2nd PDF component)    [units vary]

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
      mu_x_1_n,    & ! Mean array (normalized) of PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normalized) of PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normalized) of PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normalized) of PDF vars (comp. 2) [u.v.]

    ! Local Variable
    integer :: ivar  ! Loop index


    ! The means and standard deviations in each PDF component of w, chi (old s),
    ! and eta (old t) do not need to be normalized, since w, chi, and eta
    ! already follow assumed normal distributions in each PDF component.  The
    ! normalized means and standard deviations are the same as the actual means
    ! and standard deviations.    
    mu_x_1_n = mu_x_1
    mu_x_2_n = mu_x_2
    sigma_x_1_n = sigma_x_1
    sigma_x_2_n = sigma_x_2

    !!! Calculate the normalized mean and standard deviation in each PDF
    !!! component for variables that have an assumed lognormal distribution,
    !!! given the mean and standard deviation in each PDF component for those
    !!! variables.  A precipitating hydrometeor has an assumed lognormal
    !!! distribution in precipitation in each PDF component.  Simplified cloud
    !!! nuclei concentration, N_cn, has an assumed lognormal distribution in
    !!! each PDF component, and furthermore, mu_Ncn_1 = mu_Ncn_2 and
    !!! sigma_Ncn_1 = sigma_Ncn_2, so N_cn has an assumed single lognormal
    !!! distribution over the entire domain.

    ! Normalized mean of simplified cloud nuclei concentration, N_cn,
    ! in PDF component 1.
    if ( Ncnm > Ncn_tol ) then

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

    ! Normalized mean of simplified cloud nuclei concentration, N_cn,
    ! in PDF component 2.
    if ( Ncnm > Ncn_tol ) then

       mu_x_2_n(iiPDF_Ncn) = mean_L2N( mu_x_2(iiPDF_Ncn), sigma2_on_mu2_ip_1(iiPDF_Ncn) )

    else

       ! Mean simplified cloud nuclei concentration in PDF component 1 is less
       ! than the tolerance amount.  It is considered to have a value of 0.
       ! There are not any cloud nuclei or cloud at this grid level.  The value
       ! of mu_Ncn_1_n should be -inf.  It will be set to -huge for purposes of
       ! assigning it a value.
       mu_x_2_n(iiPDF_Ncn) = -huge( mu_x_2(iiPDF_Ncn) )

    endif

    ! Normalized standard deviation of simplified cloud nuclei concentration,
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

    ! Normalize precipitating hydrometeor means and standard deviations.
    do ivar = iiPDF_Ncn+1, d_variables, 1

       ! Normalized mean of a precipitating hydrometeor, hm, in PDF component 1.
       if ( hm1(pdf2hydromet_idx(ivar)) &
            > hydromet_tol(pdf2hydromet_idx(ivar)) ) then

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

       ! Normalized standard deviation of a precipitating hydrometeor, hm, in
       ! PDF component 1.
       sigma_x_1_n(ivar) = stdev_L2N( sigma2_on_mu2_ip_1(ivar) )

       ! Normalized mean of a precipitating hydrometeor, hm, in PDF component 2.
       if ( hm2(pdf2hydromet_idx(ivar)) &
            > hydromet_tol(pdf2hydromet_idx(ivar)) ) then

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

       ! Normalized standard deviation of a precipitating hydrometeor, hm, in
       ! PDF component 2.
       sigma_x_2_n(ivar) = stdev_L2N( sigma2_on_mu2_ip_2(ivar) )

    enddo ! ivar = iiPDF_Ncn+1, d_variables, 1


    return

  end subroutine normalize_mean_stdev

  !=============================================================================
  subroutine normalize_corr( d_variables, sigma_x_1_n, sigma_x_2_n, &
                             sigma2_on_mu2_ip_1, sigma2_on_mu2_ip_2, &
                             corr_array_1, corr_array_2, &
                             corr_array_1_n, corr_array_2_n )

    ! Description:
    ! Calculates the normalized correlations between PDF variables, where at
    ! least one of the variables that is part of a correlation has an assumed
    ! lognormal distribution -- which are the precipitating hydrometeors (in
    ! precipitation) and N_cn.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        zero          ! Constant

    use pdf_utilities, only: &
        corr_NL2NN, & ! Procedure(s)
        corr_LL2NN

    use corr_varnce_module, only: &
        iiPDF_chi, & ! Variable(s)
        iiPDF_eta, &
        iiPDF_w,  &
        iiPDF_Ncn

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use model_flags, only: &
        l_const_Nc_in_cloud  ! Variable!!

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      d_variables ! Number of PDF variables

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      sigma_x_1_n, & ! Std. dev. array (normalized) of PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normalized) of PDF vars (comp. 2) [u.v.]

    real ( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      sigma2_on_mu2_ip_1, & ! Prescribed ratio array: sigma_hm_1^2/mu_hm_1^2 [-]
      sigma2_on_mu2_ip_2    ! Prescribed ratio array: sigma_hm_2^2/mu_hm_2^2 [-]

    real( kind = core_rknd ), dimension(d_variables, d_variables), &
    intent(in) :: &
      corr_array_1, & ! Correlation array of PDF vars. (comp. 1)             [-]
      corr_array_2    ! Correlation array of PDF vars. (comp. 2)             [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(d_variables, d_variables), &
    intent(out) :: &
      corr_array_1_n, & ! Corr. array (normalized) of PDF vars. (comp. 1)    [-]
      corr_array_2_n    ! Corr. array (normalized) of PDF vars. (comp. 2)    [-]

    ! Local Variables
    integer :: ivar, jvar ! Loop indices


    ! The correlations in each PDF component between two of w, chi (old s), and
    ! eta (old t) do not need to be normalized, since w, chi, and eta already
    ! follow assumed normal distributions in each PDF component.  The normalized
    ! correlations between any two of these variables are the same as the actual
    ! correlations.    
    corr_array_1_n = corr_array_1
    corr_array_2_n = corr_array_2

    !!! Calculate the normalized correlation of variables that have
    !!! an assumed normal distribution and variables that have an assumed
    !!! lognormal distribution for the ith PDF component, given their
    !!! correlation and the normalized standard deviation of the variable with
    !!! the assumed lognormal distribution.

    if ( l_const_Nc_in_cloud ) then

      ! Ncn does not vary in the grid box. Consequently, the correlation between
      ! Ncn and any other variate is not defined. Here, we set the correlations
      ! between Ncn and chi/eta/w to zero.
      corr_array_1_n(iiPDF_Ncn, iiPDF_w) = zero
      corr_array_2_n(iiPDF_Ncn, iiPDF_w) = zero
      corr_array_1_n(iiPDF_Ncn, iiPDF_chi) = zero
      corr_array_2_n(iiPDF_Ncn, iiPDF_chi) = zero
      corr_array_1_n(iiPDF_Ncn, iiPDF_eta) = zero
      corr_array_2_n(iiPDF_Ncn, iiPDF_eta) = zero

    else ! .not. l_const_Nc_in_cloud

      ! Normalize the correlations between chi/eta/w and N_cn.

      ! Normalize the correlation of w and N_cn in PDF component 1.
      corr_array_1_n(iiPDF_Ncn, iiPDF_w) &
      = corr_NL2NN( corr_array_1(iiPDF_Ncn, iiPDF_w), sigma_x_1_n(iiPDF_Ncn), &
                    sigma2_on_mu2_ip_1(iiPDF_Ncn) )

      ! Normalize the correlation of w and N_cn in PDF component 2.
      corr_array_2_n(iiPDF_Ncn, iiPDF_w) &
      = corr_NL2NN( corr_array_2(iiPDF_Ncn, iiPDF_w), sigma_x_2_n(iiPDF_Ncn), &
                    sigma2_on_mu2_ip_1(iiPDF_Ncn) )

      ! Normalize the correlation of chi (old s) and N_cn in PDF component 1.
      corr_array_1_n(iiPDF_Ncn, iiPDF_chi) &
      = corr_NL2NN( corr_array_1(iiPDF_Ncn, iiPDF_chi), &
                    sigma_x_1_n(iiPDF_Ncn), sigma2_on_mu2_ip_1(iiPDF_Ncn) )

      ! Normalize the correlation of chi (old s) and N_cn in PDF component 2.
      corr_array_2_n(iiPDF_Ncn, iiPDF_chi) &
      = corr_NL2NN( corr_array_2(iiPDF_Ncn, iiPDF_chi), &
                    sigma_x_2_n(iiPDF_Ncn), sigma2_on_mu2_ip_1(iiPDF_Ncn) )

      ! Normalize the correlation of eta (old t) and N_cn in PDF component 1.
      corr_array_1_n(iiPDF_Ncn, iiPDF_eta) &
      = corr_NL2NN( corr_array_1(iiPDF_Ncn, iiPDF_eta), &
                    sigma_x_1_n(iiPDF_Ncn), sigma2_on_mu2_ip_1(iiPDF_Ncn) )

      ! Normalize the correlation of eta (old t) and N_cn in PDF component 2.
      corr_array_2_n(iiPDF_Ncn, iiPDF_eta) &
      = corr_NL2NN( corr_array_2(iiPDF_Ncn, iiPDF_eta), &
                    sigma_x_2_n(iiPDF_Ncn), sigma2_on_mu2_ip_1(iiPDF_Ncn) )

    end if ! l_const_Nc_in_cloud

    ! Normalize the correlations (in-precip) between chi/eta/w and the
    ! precipitating hydrometeors.
    do ivar = iiPDF_chi, iiPDF_w
       do jvar = iiPDF_Ncn+1, d_variables

          ! Normalize the correlation (in-precip) between w, chi, or eta and a
          ! precipitating hydrometeor, hm, in PDF component 1.
          corr_array_1_n(jvar, ivar) &
          = corr_NL2NN( corr_array_1(jvar, ivar), sigma_x_1_n(jvar), &
                        sigma2_on_mu2_ip_1(jvar) )

          ! Normalize the correlation (in-precip) between w, chi, or eta and a
          ! precipitating hydrometeor, hm, in PDF component 2.
          corr_array_2_n(jvar, ivar) &
          = corr_NL2NN( corr_array_2(jvar, ivar), sigma_x_2_n(jvar), &
                        sigma2_on_mu2_ip_2(jvar) )

       enddo ! jvar = iiPDF_Ncn+1, d_variables
    enddo ! ivar = iiPDF_chi, iiPDF_w


    !!! Calculate the normalized correlation of two variables that both
    !!! have an assumed lognormal distribution for the ith PDF component, given
    !!! their correlation and both of their normalized standard deviations.

    ! Normalize the correlations (in-precip) between N_cn and the precipitating
    ! hydrometeors.
    ivar = iiPDF_Ncn
    do jvar = ivar+1, d_variables

       if ( l_const_Nc_in_cloud ) then

         ! Ncn does not vary, so these correlations are undefined. Set them to
         ! zero.
         corr_array_1_n(jvar,ivar) = zero
         corr_array_2_n(jvar,ivar) = zero

       else ! .not. l_const_Nc_in_cloud

         ! Normalize the correlation (in-precip) between N_cn and a precipitating
         ! hydrometeor, hm, in PDF component 1.
         corr_array_1_n(jvar, ivar) &
         = corr_LL2NN( corr_array_1(jvar, ivar), &
                       sigma_x_1_n(ivar), sigma_x_1_n(jvar), &
                       sigma2_on_mu2_ip_1(iiPDF_Ncn), sigma2_on_mu2_ip_1(jvar) )

         ! Normalize the correlation (in-precip) between N_cn and a precipitating
         ! hydrometeor, hm, in PDF component 2.
         corr_array_2_n(jvar, ivar) &
         = corr_LL2NN( corr_array_2(jvar, ivar), &
                       sigma_x_2_n(ivar), sigma_x_2_n(jvar), &
                       sigma2_on_mu2_ip_1(iiPDF_Ncn), sigma2_on_mu2_ip_2(jvar) )

       end if ! l_const_Nc_in_cloud

    enddo ! jvar = ivar+1, d_variables

    ! Normalize the correlations (in-precip) between two precipitating
    ! hydrometeors.
    do ivar = iiPDF_Ncn+1, d_variables-1
       do jvar = ivar+1, d_variables

          ! Normalize the correlation (in-precip) between two precipitating
          ! hydrometeors (for example, r_r and N_r) in PDF component 1.
          corr_array_1_n(jvar, ivar) &
          = corr_LL2NN( corr_array_1(jvar, ivar), &
                        sigma_x_1_n(ivar), sigma_x_1_n(jvar), &
                        sigma2_on_mu2_ip_1(ivar), sigma2_on_mu2_ip_1(jvar) )

          ! Normalize the correlation (in-precip) between two precipitating
          ! hydrometeors (for example, r_r and N_r) in PDF component 2.
          corr_array_2_n(jvar, ivar) &
          = corr_LL2NN( corr_array_2(jvar, ivar), &
                        sigma_x_2_n(ivar), sigma_x_2_n(jvar), &
                        sigma2_on_mu2_ip_2(ivar), sigma2_on_mu2_ip_2(jvar) )

       enddo ! jvar = ivar+1, d_variables
    enddo ! ivar = iiPDF_Ncn+1, d_variables-1


    return

  end subroutine normalize_corr

  !=============================================================================
  subroutine calc_corr_w_hm( wm, wphydrometp, &
                             mu_w_1, mu_w_2, &
                             mu_hm_1, mu_hm_2, &
                             sigma_w_1, sigma_w_2, &
                             sigma_hm_1, sigma_hm_2, &
                             mixt_frac, precip_frac_1, precip_frac_2, &
                             corr_w_hm_1, corr_w_hm_2, &
                             hm_tol )

    ! Description:
    ! Calculates the PDF component correlation (in-precip) between vertical
    ! velocity, w, and a hydrometeor, hm.  The overall covariance of w and hm,
    ! <w'hm'> can be written in terms of the PDF parameters.  When both w and hm
    ! vary in both PDF components, the equation is written as:
    !
    ! <w'hm'> = mixt_frac * precip_frac_1
    !           * ( ( mu_w_1 - <w> ) * mu_hm_1
    !               + corr_w_rr_1 * sigma_w_1 * sigma_rr_1 )
    !           + ( 1 - mixt_frac ) * precip_frac_2
    !             * ( ( mu_w_2 - <w> ) * mu_hm_2
    !                 + corr_w_rr_2 * sigma_w_2 * sigma_rr_2 ).
    !
    ! The overall covariance is provided, so the component correlation is solved
    ! by setting corr_w_rr_1 = corr_w_rr_2 ( = corr_w_rr ).  The equation is:
    !
    ! corr_w_rr
    ! = ( <w'hm'>
    !     - mixt_frac * precip_frac_1 * ( mu_w_1 - <w> ) * mu_hm_1
    !     - ( 1 - mixt_frac ) * precip_frac_2 * ( mu_w_2 - <w> ) * mu_hm_2 )
    !   / ( mixt_frac * precip_frac_1 * sigma_w_1 * sigma_hm_1
    !       + ( 1 - mixt_frac ) * precip_frac_2 * sigma_w_2 * sigma_hm_2 );
    !
    ! again, where corr_w_rr_1 = corr_w_rr_2 = corr_w_rr.  When either w or hm
    ! isbconstant in one PDF component, but both w and hm vary in the other PDF
    ! component, the equation for <w'hm'> is written as:
    !
    ! <w'hm'> = mixt_frac * precip_frac_1
    !           * ( ( mu_w_1 - <w> ) * mu_hm_1
    !               + corr_w_rr_1 * sigma_w_1 * sigma_rr_1 )
    !           + ( 1 - mixt_frac ) * precip_frac_2
    !             * ( mu_w_2 - <w> ) * mu_hm_2.
    !
    ! In the above equation, either w or hm (or both) is (are) constant in PDF
    ! component 2, but both w and hm vary in PDF component 1.  When both w and
    ! hm vary in PDF component 2, but at least one of w or hm is constant in PDF
    ! component 1, the equation is similar.  The above equation can be rewritten
    ! to solve for corr_w_rr_1, such that:
    !
    ! corr_w_rr_1
    ! = ( <w'hm'>
    !     - mixt_frac * precip_frac_1 * ( mu_w_1 - <w> ) * mu_hm_1
    !     - ( 1 - mixt_frac ) * precip_frac_2 * ( mu_w_2 - <w> ) * mu_hm_2 )
    !   / ( mixt_frac * precip_frac_1 * sigma_w_1 * sigma_hm_1 ).
    !
    ! Since either w or hm is constant in PDF component 2, corr_w_rr_2 is
    ! undefined.  When both w and hm vary in PDF component 2, but at least one
    ! of w or hm is constant in PDF component 1, the equation is similar, but
    ! is in terms of corr_w_rr_2, while corr_w_rr_1 is undefined.  When either w
    ! or hm is constant in both PDF components, the equation for <w'hm'> is:
    !
    ! <w'hm'> = mixt_frac * precip_frac_1
    !           * ( mu_w_1 - <w> ) * mu_hm_1
    !           + ( 1 - mixt_frac ) * precip_frac_2
    !             * ( mu_w_2 - <w> ) * mu_hm_2.
    !
    ! When this is the case, both corr_w_rr_1 and corr_w_rr_2 are undefined.

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
      mixt_frac,     & ! Mixture fraction                                    [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2, & ! Precipitation fraction (2nd PDF component)          [-]
      hm_tol           ! Hydrometeor tolerance value                     [hm un]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      corr_w_hm_1, & ! Correlation of w and hm (1st PDF component) ip    [-]
      corr_w_hm_2    ! Correlation of w and hm (2nd PDF component) ip    [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      corr_w_hm    ! Correlation of w and hm (both PDF components) ip    [-]


    ! Calculate the PDF component correlation of vertical velocity, w, and
    ! a hydrometeor, hm, in precipitation.
    if ( sigma_w_1 > w_tol .and. sigma_hm_1 > hm_tol .and. &
         sigma_w_2 > w_tol .and. sigma_hm_2 > hm_tol ) then

       ! Both w and hm vary in both PDF components.
       ! Calculate corr_w_hm (where corr_w_hm_1 = corr_w_hm_2 = corr_w_hm).
       corr_w_hm &
       = ( wphydrometp &
           - mixt_frac * precip_frac_1 * ( mu_w_1 - wm ) * mu_hm_1 &
           - ( one - mixt_frac ) * precip_frac_2 * ( mu_w_2 - wm ) * mu_hm_2 ) &
         / ( mixt_frac * precip_frac_1 * sigma_w_1 * sigma_hm_1 &
             + ( one - mixt_frac ) * precip_frac_2 * sigma_w_2 * sigma_hm_2 )

       ! Check that the PDF component correlations have reasonable values.
       if ( corr_w_hm > max_mag_correlation ) then
          corr_w_hm = max_mag_correlation
       elseif ( corr_w_hm < -max_mag_correlation ) then
          corr_w_hm = -max_mag_correlation
       endif

       ! The PDF component correlations between w and hm (in-precip) are equal.
       corr_w_hm_1 = corr_w_hm
       corr_w_hm_2 = corr_w_hm


    elseif ( sigma_w_1 > w_tol .and. sigma_hm_1 > hm_tol ) then

       ! Both w and hm vary in PDF component 1, but at least one of w and hm is
       ! constant in PDF component 2.
       ! Calculate the PDF component 1 correlation of w and hm (in-precip).
       corr_w_hm_1 &
       = ( wphydrometp &
           - mixt_frac * precip_frac_1 * ( mu_w_1 - wm ) * mu_hm_1 &
           - ( one - mixt_frac ) * precip_frac_2 * ( mu_w_2 - wm ) * mu_hm_2 ) &
         / ( mixt_frac * precip_frac_1 * sigma_w_1 * sigma_hm_1 )

       ! Check that the PDF component 1 correlation has a reasonable value.
       if ( corr_w_hm_1 > max_mag_correlation ) then
          corr_w_hm_1 = max_mag_correlation
       elseif ( corr_w_hm_1 < -max_mag_correlation ) then
          corr_w_hm_1 = -max_mag_correlation
       endif

       ! The PDF component 2 correlation is undefined.
       corr_w_hm_2 = zero
       

    elseif ( sigma_w_2 > w_tol .and. sigma_hm_2 > hm_tol ) then

       ! Both w and hm vary in PDF component 2, but at least one of w and hm is
       ! constant in PDF component 1.
       ! Calculate the PDF component 2 correlation of w and hm (in-precip).
       corr_w_hm_2 &
       = ( wphydrometp &
           - mixt_frac * precip_frac_1 * ( mu_w_1 - wm ) * mu_hm_1 &
           - ( one - mixt_frac ) * precip_frac_2 * ( mu_w_2 - wm ) * mu_hm_2 ) &
         / ( ( one - mixt_frac ) * precip_frac_2 * sigma_w_2 * sigma_hm_2 )

       ! Check that the PDF component 2 correlation has a reasonable value.
       if ( corr_w_hm_2 > max_mag_correlation ) then
          corr_w_hm_2 = max_mag_correlation
       elseif ( corr_w_hm_2 < -max_mag_correlation ) then
          corr_w_hm_2 = -max_mag_correlation
       endif

       ! The PDF component 1 correlation is undefined.
       corr_w_hm_1 = zero
       

    else    ! sigma_w_1 * sigma_hm_1 = 0 .and. sigma_w_2 * sigma_hm_2 = 0.

       ! At least one of w and hm is constant in both PDF components.

       ! The PDF component 1 and component 2 correlations are both undefined.
       corr_w_hm_1 = zero
       corr_w_hm_2 = zero


    endif


    return

  end subroutine calc_corr_w_hm

  !=============================================================================
  subroutine pdf_param_hm_stats( d_variables, level, mu_x_1, mu_x_2, &
                                 sigma_x_1, sigma_x_2, &
                                 corr_array_1, corr_array_2, &
                                 l_stats_samp )

    ! Description:
    ! Record statistics for standard PDF parameters involving hydrometeors.

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
        imu_hm_1,     & ! Variable(s)
        imu_hm_2,     &
        imu_Ncn_1,    &
        imu_Ncn_2,    &
        isigma_hm_1,  &
        isigma_hm_2,  &
        isigma_Ncn_1, &
        isigma_Ncn_2

    use stats_variables, only : &
        icorr_w_chi_1,       & ! Variable(s)
        icorr_w_chi_2,       &
        icorr_w_eta_1,       &
        icorr_w_eta_2,       &
        icorr_w_hm_1,        &
        icorr_w_hm_2,        &
        icorr_w_Ncn_1,       &
        icorr_w_Ncn_2,       &
        icorr_chi_eta_1_ca,  &
        icorr_chi_eta_2_ca,  &
        icorr_chi_hm_1,      &
        icorr_chi_hm_2,      &
        icorr_chi_Ncn_1,     &
        icorr_chi_Ncn_2,     &
        icorr_eta_hm_1,      &
        icorr_eta_hm_2,      &
        icorr_eta_Ncn_1,     &
        icorr_eta_Ncn_2,     &
        icorr_Ncn_hm_1,       &
        icorr_Ncn_hm_2,       &
        icorr_hmx_hmy_1,     &
        icorr_hmx_hmy_2,     &
        stats_zt

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      d_variables, & ! Number of variables in the correlation array
      level          ! Vertical level index 

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
          call stat_update_var_pt( imu_Ncn_1, level, mu_x_1(iiPDF_Ncn), stats_zt )
       endif

       ! Mean of cloud nuclei concentration in PDF component 2.
       if ( imu_Ncn_2 > 0 ) then
          call stat_update_var_pt( imu_Ncn_2, level, mu_x_2(iiPDF_Ncn), stats_zt )
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
                                      level, corr_array_1(ivar,iiPDF_w), stats_zt )
          endif

          ! Correlation (in-precip) of w and the precipitating hydrometeor
          ! in PDF component 2.
          if ( icorr_w_hm_2(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt( icorr_w_hm_2(pdf2hydromet_idx(ivar)), &
                                      level, corr_array_2(ivar,iiPDF_w), stats_zt )
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
                                      level, corr_array_1(ivar,iiPDF_chi), stats_zt )
          endif

          ! Correlation (in-precip) of chi (old s) and the precipitating
          ! hydrometeor in PDF component 2.
          if ( icorr_chi_hm_2(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt( icorr_chi_hm_2(pdf2hydromet_idx(ivar)), &
                                      level, corr_array_2(ivar,iiPDF_chi), stats_zt )
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
                                      level, corr_array_1(ivar,iiPDF_eta), stats_zt )
          endif

          ! Correlation (in-precip) of eta (old t) and the precipitating
          ! hydrometeor in PDF component 2.
          if ( icorr_eta_hm_2(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt( icorr_eta_hm_2(pdf2hydromet_idx(ivar)), &
                                      level, corr_array_2(ivar,iiPDF_eta), stats_zt )
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
                                      level, corr_array_1(ivar,iiPDF_Ncn), stats_zt )
          endif

          ! Correlation (in-precip) of N_cn and the precipitating
          ! hydrometeor in PDF component 2.
          if ( icorr_Ncn_hm_2(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt( icorr_Ncn_hm_2(pdf2hydromet_idx(ivar)), &
                                      level, corr_array_2(ivar,iiPDF_Ncn), stats_zt )
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
    ! Record statistics for normalized PDF parameters involving hydrometeors.

    ! References:
    !-----------------------------------------------------------------------

    use index_mapping, only: &
        pdf2hydromet_idx  ! Procedure(s)

    use corr_varnce_module, only: &
        iiPDF_w,                   & ! Variable(s)
        iiPDF_chi,            &
        iiPDF_eta,            &
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
        icorr_chi_hm_1_n,    &
        icorr_chi_hm_2_n,    &
        icorr_chi_Ncn_1_n,   &
        icorr_chi_Ncn_2_n,   &
        icorr_eta_hm_1_n,    &
        icorr_eta_hm_2_n,    &
        icorr_eta_Ncn_1_n,   &
        icorr_eta_Ncn_2_n,   &
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
      mu_x_1_n,    & ! Mean array (normalized) of PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normalized) of PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normalized) of PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normalized) of PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(d_variables, d_variables), &
    intent(in) :: &
      corr_array_1_n, & ! Corr. array (normalized) of PDF vars. (comp. 1)    [-]
      corr_array_2_n    ! Corr. array (normalized) of PDF vars. (comp. 2)    [-]

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.

    ! Local Variable
    integer :: ivar, jvar  ! Loop indices


    !!! Output the statistics for normalized hydrometeor PDF parameters.

    ! Statistics
    if ( l_stats_samp ) then

       do ivar = iiPDF_Ncn+1, d_variables, 1

          ! Mean (in-precip) of ln hm in PDF component 1.
          if ( imu_hm_1_n(pdf2hydromet_idx(ivar)) > 0 ) then
             if ( mu_x_1_n(ivar) > real( -huge( 0.0 ), kind = core_rknd ) ) then
                call stat_update_var_pt( imu_hm_1_n(pdf2hydromet_idx(ivar)), &
                                         level, mu_x_1_n(ivar), stats_zt )
             else
                ! When hm1 is 0 (or below tolerance value), mu_hm_1_n is -inf,
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
                ! When hm2 is 0 (or below tolerance value), mu_hm_2_n is -inf,
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
                                      level, corr_array_1_n(ivar,iiPDF_w), stats_zt )
          endif

          ! Correlation (in-precip) of w and ln hm in PDF component 2.
          if ( icorr_w_hm_2_n(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt( icorr_w_hm_2_n(pdf2hydromet_idx(ivar)), &
                                      level, corr_array_2_n(ivar,iiPDF_w), stats_zt )
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
                                     level, corr_array_1_n(ivar,iiPDF_chi), stats_zt )
          endif

          ! Correlation (in-precip) of chi( old s) and ln hm in PDF component 2.
          if ( icorr_chi_hm_2_n(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt(icorr_chi_hm_2_n(pdf2hydromet_idx(ivar)), &
                                     level, corr_array_2_n(ivar,iiPDF_chi), stats_zt )
          endif

       enddo ! ivar = iiPDF_Ncn+1, d_variables, 1

       ! Correlation of chi (old s) and ln N_cn in PDF component 1.
       if ( icorr_chi_Ncn_1_n > 0 ) then
          call stat_update_var_pt( icorr_chi_Ncn_1_n, level, &
                                   corr_array_1_n(iiPDF_Ncn,iiPDF_chi), stats_zt )
       endif

       ! Correlation of chi(old s) and ln N_cn in PDF component 2.
       if ( icorr_chi_Ncn_2_n > 0 ) then
          call stat_update_var_pt( icorr_chi_Ncn_2_n, level, &
                                   corr_array_2_n(iiPDF_Ncn,iiPDF_chi), stats_zt )
       endif

       do ivar = iiPDF_Ncn+1, d_variables, 1

          ! Correlation (in-precip) of eta (old t) and ln hm in PDF component 1.
          if ( icorr_eta_hm_1_n(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt(icorr_eta_hm_1_n(pdf2hydromet_idx(ivar)), &
                                     level, corr_array_1_n(ivar,iiPDF_eta), stats_zt )
          endif

          ! Correlation (in-precip) of eta (old t) and ln hm in PDF component 2.
          if ( icorr_eta_hm_2_n(pdf2hydromet_idx(ivar)) > 0 ) then
             call stat_update_var_pt(icorr_eta_hm_2_n(pdf2hydromet_idx(ivar)), &
                                     level, corr_array_2_n(ivar,iiPDF_eta), stats_zt )
          endif

       enddo ! ivar = iiPDF_Ncn+1, d_variables, 1

       ! Correlation of eta (old t) and ln N_cn in PDF component 1.
       if ( icorr_eta_Ncn_1_n > 0 ) then
          call stat_update_var_pt( icorr_eta_Ncn_1_n, level, &
                                   corr_array_1_n(iiPDF_Ncn,iiPDF_eta), stats_zt )
       endif

       ! Correlation of eta (old t) and ln N_cn in PDF component 2.
       if ( icorr_eta_Ncn_2_n > 0 ) then
          call stat_update_var_pt( icorr_eta_Ncn_2_n, level, &
                                   corr_array_2_n(iiPDF_Ncn,iiPDF_eta), stats_zt )
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
  subroutine pack_pdf_params( hm1, hm2, d_variables, &                   ! In
                              mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &    ! In
                              corr_array_1, corr_array_2, precip_frac, & ! In
                              precip_frac_1, precip_frac_2, &            ! In
                              hydromet_pdf_params )                      ! Out

    ! Description:
    ! Pack the standard means and variances involving hydrometeors, as well as a
    ! few other variables, into the structure hydromet_pdf_params.

    ! References:
    !-----------------------------------------------------------------------

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
      hm1, & ! Mean of a precip. hydrometeor (1st PDF component)  [units vary]
      hm2    ! Mean of a precip. hydrometeor (2nd PDF component)  [units vary]

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
    integer :: ivar  ! Loop index


    ! Pack remaining means and standard deviations into hydromet_pdf_params.
    do ivar = 1, hydromet_dim, 1

       ! Mean of a hydrometeor (overall) in the 1st PDF component.
       hydromet_pdf_params%hm1(ivar) = hm1(ivar)
       ! Mean of a hydrometeor (overall) in the 2nd PDF component.
       hydromet_pdf_params%hm2(ivar) = hm2(ivar)

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
      sigma_chi_1,        & ! Standard deviation of chi (1st PDF comp.)  [kg/kg]
      sigma_chi_2,        & ! Standard deviation of chi (2nd PDF comp.)  [kg/kg]
      sigma_eta_1,        & ! Standard deviation of eta (1st PDF comp.)  [kg/kg]
      sigma_eta_2,        & ! Standard deviation of eta (2nd PDF comp.)  [kg/kg]
      crt_1,               & ! Coef. of r_t in chi/eta eqns. (1st comp.)  [-]
      crt_2,               & ! Coef. of r_t in chi/eta eqns. (2nd comp.)  [-]
      rt_1,                & ! Mean of rt (1st PDF component)             [kg/kg]
      rt_2,                & ! Mean of rt (2nd PDF component)             [kg/kg]
      rtm,                & ! Mean of rt (overall)                       [kg/kg]
      sigma_rt_1_from_chi, & ! Standard deviation of rt (1st PDF comp.)   [kg/kg]
      sigma_rt_2_from_chi, & ! Standard deviation of rt (2nd PDF comp.)   [kg/kg]
      mixt_frac             ! Weight of 1st gaussian PDF component       [-]

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    ! Enter some PDF parameters
    sigma_chi_1 = pdf_params%stdev_chi_1
    sigma_chi_2 = pdf_params%stdev_chi_2
    sigma_eta_1 = pdf_params%stdev_eta_1
    sigma_eta_2 = pdf_params%stdev_eta_2
    rt_1         = pdf_params%rt_1
    rt_2         = pdf_params%rt_2
    crt_1        = pdf_params%crt_1
    crt_2        = pdf_params%crt_2
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
