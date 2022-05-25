!-------------------------------------------------------------------------
! $Id$
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

  private :: component_corr_w_x,          &
             component_corr_chi_eta,      &
             component_corr_w_hm_n_ip,    &
             component_corr_x_hm_n_ip,    &
             component_corr_hmx_hmy_n_ip, &
             component_corr_eta_hm_n_ip,  &
             calc_corr_w_hm_n,            &
             pdf_param_hm_stats,          &
             pdf_param_ln_hm_stats,       &
             pack_hydromet_pdf_params,    &
             compute_rtp2_from_chi
    
  ! The component correlations of w and r_t and the component correlations of
  ! w and theta_l must both be 0 when using the ADG1 PDF.  In other words, w
  ! and r_t (theta_l) have overall covariance w'r_t' (w'theta_l'), but the
  ! individual component covariance and correlation are defined to be 0.
  ! Since the component covariances (or correlations) of w and chi (old s) and
  ! of w and eta (old t) are based on the covariances (or correlations) of w
  ! and r_t and of w and theta_l, the individual component correlation and
  ! covariance of w and chi, as well of as w and eta, are defined to be 0.
  !
  ! The PDF component correlation of w and x, calculated as part of the PDF
  ! parameters, comes out to be 0 (within numerical roundoff) when the ADG1
  ! PDF is used.  However, when l_fix_w_chi_eta_correlations is enabled, the
  ! component correlations of PDF variables are specified, and the values
  ! that were calculated as part of the PDF parameters are ignored in the
  ! correlation array.  When this happens, and when the ADG1 PDF is selected,
  ! the l_follow_ADG1_PDF_standards flag can be used to force the component
  ! correlations of w and x to have a value of 0, following the ADG1 standard,
  ! rather than whatever value might be specified in the correlation array.
  !
  ! Note:  the ADG2 PDF follows the same standards as the ADG1 PDF.
  logical, parameter :: &
    l_follow_ADG1_PDF_standards = .true.

  contains

  !=============================================================================
  subroutine setup_pdf_parameters( gr, nz, ngrdcol, pdf_dim, dt, &             ! Intent(in)
                                   Nc_in_cloud, rcm, cloud_frac, Kh_zm, &      ! Intent(in)
                                   ice_supersat_frac, hydromet, wphydrometp, & ! Intent(in)
                                   corr_array_n_cloud, corr_array_n_below, &   ! Intent(in)
                                   pdf_params, l_stats_samp, &                 ! Intent(in)
                                   clubb_params, &                             ! Intent(in)
                                   iiPDF_type, &                               ! Intent(in)
                                   l_use_precip_frac, &                        ! Intent(in)
                                   l_predict_upwp_vpwp, &                      ! Intent(in)
                                   l_diagnose_correlations, &                  ! Intent(in)
                                   l_calc_w_corr, &                            ! Intent(in)
                                   l_const_Nc_in_cloud, &                      ! Intent(in)
                                   l_fix_w_chi_eta_correlations, &             ! Intent(in)
                                   stats_zt, stats_zm, stats_sfc, &            ! intent(inout)
                                   hydrometp2, &                               ! Intent(out)
                                   mu_x_1_n, mu_x_2_n, &                       ! Intent(out)
                                   sigma_x_1_n, sigma_x_2_n, &                 ! Intent(out)
                                   corr_array_1_n, corr_array_2_n, &           ! Intent(out)
                                   corr_cholesky_mtx_1, corr_cholesky_mtx_2, & ! Intent(out)
                                   precip_fracs, &                             ! Intent(out)
                                   hydromet_pdf_params )                       ! Intent(out)

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid, & ! Type
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
        precipitation_fractions,  &
        init_hydromet_pdf_params     ! Procedure

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use array_index, only: &
        hydromet_list, &       ! Variable(s)
        hydromet_tol,  &
        iiPDF_Ncn,     & 
        iiPDF_chi,     &
        iiPDF_eta

    use precipitation_fraction, only: &
        precip_fraction

    use Nc_Ncn_eqns, only: &
        Nc_in_cloud_to_Ncnm  ! Procedure(s)

    use parameter_indices, only: &
        nparams,         & ! Variable(s)
        ic_K_hm,         &
        iomicron,        &
        izeta_vrnce_rat

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
        core_rknd        ! Variable(s)

    use matrix_operations, only: &
        Cholesky_factor  ! Procedure(s)

    use stats_type_utilities, only: &
!        stat_update_var,    & ! Procedure(s)
        stat_update_var_pt

    use stats_variables, only: &
        iprecip_frac,   & ! Variable(s)
        iprecip_frac_1, &
        iprecip_frac_2, &
        iNcnm,          &
        ihmp2_zt,       &
        irtp2_from_chi

    use diagnose_correlations_module, only: &
        diagnose_correlations, & ! Procedure(s)
        calc_cholesky_corr_mtx_approx

    use corr_varnce_module, only: &
        assert_corr_symmetric, & ! Procedure(s)
        hmp2_ip_on_hmm2_ip,    & ! Variable(s)
        Ncnp2_on_Ncnm2

    use error_code, only: &
        clubb_at_least_debug_level, &   ! Procedure
        err_code, &                     ! Error Indicator
        clubb_fatal_error               ! Constant

    use stats_type, only: stats ! Type
    
    use advance_helper_module, only : &
        calc_xpwp  ! Procedure(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,          & ! Number of model vertical grid levels
      pdf_dim,     & ! Number of variables in the correlation array
      ngrdcol        ! Number of grid columns
      
    type (grid), target, dimension(ngrdcol), intent(in) :: gr

    real( kind = core_rknd ), intent(in) ::  &
      dt    ! Model timestep                                           [s]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      Nc_in_cloud,       & ! Mean (in-cloud) cloud droplet conc.       [num/kg]
      rcm,               & ! Mean cloud water mixing ratio, < r_c >    [kg/kg]
      cloud_frac,        & ! Cloud fraction                            [-]
      Kh_zm,             & ! Eddy diffusivity coef. on momentum levels [m^2/s]
      ice_supersat_frac    ! Ice supersaturation fraction              [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz,hydromet_dim), intent(in) :: &
      hydromet,    & ! Mean of hydrometeor, hm (overall) (t-levs.) [units]
      wphydrometp    ! Covariance < w'h_m' > (momentum levels)     [(m/s)units]

    real( kind = core_rknd ), dimension(pdf_dim,pdf_dim), &
    intent(in) :: &
      corr_array_n_cloud, & ! Prescribed normal space corr. array in cloud  [-]
      corr_array_n_below    ! Prescribed normal space corr. array below cl. [-]

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                               [units vary]

    logical, intent(in) :: &
      l_stats_samp    ! Flag to sample statistics

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    integer, intent(in) :: &
      iiPDF_type    ! Selected option for the two-component normal (double
                    ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                    ! w, chi, and eta) portion of CLUBB's multivariate,
                    ! two-component PDF.

    logical, intent(in) :: &
      l_use_precip_frac,            & ! Flag to use precipitation fraction in KK microphysics. The
                                      ! precipitation fraction is automatically set to 1 when this
                                      ! flag is turned off.
      l_predict_upwp_vpwp,          & ! Flag to predict <u'w'> and <v'w'> along with <u> and <v>
                                      ! alongside the advancement of <rt>, <w'rt'>, <thl>,
                                      ! <wpthlp>, <sclr>, and <w'sclr'> in subroutine
                                      ! advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still
                                      ! approximated by eddy diffusivity when <u> and <v> are
                                      ! advanced in subroutine advance_windm_edsclrm.
      l_diagnose_correlations,      & ! Diagnose correlations instead of using fixed ones
      l_calc_w_corr,                & ! Calculate the correlations between w and the hydrometeors
      l_const_Nc_in_cloud,          & ! Use a constant cloud droplet conc. within cloud (K&K)
      l_fix_w_chi_eta_correlations    ! Use a fixed correlation for s and t Mellor(chi/eta)

    ! Input/Output Variables
    type (stats), target, dimension(ngrdcol), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc
      
    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz,hydromet_dim), intent(out) :: &
      hydrometp2    ! Variance of a hydrometeor (overall) (m-levs.)   [units^2]

    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim), intent(out) :: &
      mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim,pdf_dim), &
    intent(out) :: &
      corr_array_1_n, & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
      corr_array_2_n    ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim,pdf_dim), &
    intent(out) :: &
      corr_cholesky_mtx_1, & ! Transposed corr. cholesky matrix, 1st comp. [-]
      corr_cholesky_mtx_2    ! Transposed corr. cholesky matrix, 2nd comp. [-]
      
    ! This is only an output, but it contains allocated arrays, so we need to treat it as inout
    type(precipitation_fractions), intent(inout) :: &
      precip_fracs           ! Precipitation fractions      [-]

    type(hydromet_pdf_parameter), dimension(ngrdcol,nz), optional, intent(out) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]

    ! Local Variables
    
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      Kh_zm_c_K_hm    ! Eddy diffusivity coef. on momentum levels [m^2/s]

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      sigma_w_1,    & ! Standard deviation of w (1st PDF component)      [m/s]
      sigma_w_2       ! Standard deviation of w (2nd PDF component)      [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      ice_supersat_frac_1, & ! Ice supersaturation fraction (1st PDF comp.)  [-]
      ice_supersat_frac_2    ! Ice supersaturation fraction (2nd PDF comp.)  [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      Ncnm    ! Mean cloud nuclei concentration, < N_cn >        [num/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nz) ::  &
      wpNcnp_zm, & ! Covariance of N_cn and w (momentum levs.)   [(m/s)(num/kg)]
      wpNcnp_zt    ! Covariance of N_cn and w on thermo. levels  [(m/s)(num/kg)]

    real( kind = core_rknd ), dimension(ngrdcol,nz,hydromet_dim) :: &
      hm_1, & ! Mean of a precip. hydrometeor (1st PDF component)   [units vary]
      hm_2    ! Mean of a precip. hydrometeor (2nd PDF component)   [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nz,hydromet_dim) :: &
      hydrometp2_zt,  & ! Variance of a hydrometeor (overall); t-lev   [units^2]
      wphydrometp_zt    ! Covariance of w and hm interp. to t-levs. [(m/s)units]

    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim,pdf_dim) :: &
      corr_array_1, & ! Correlation array of PDF vars. (comp. 1)             [-]
      corr_array_2    ! Correlation array of PDF vars. (comp. 2)             [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim) :: &
      mu_x_1,    & ! Mean array of PDF vars. (1st PDF component)    [units vary]
      mu_x_2,    & ! Mean array of PDF vars. (2nd PDF component)    [units vary]
      sigma_x_1, & ! Standard deviation array of PDF vars (comp. 1) [units vary]
      sigma_x_2    ! Standard deviation array of PDF vars (comp. 2) [units vary]

    real( kind = core_rknd ), dimension(pdf_dim) :: &
      corr_array_scaling

    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim) :: &
      sigma2_on_mu2_ip_1, & ! Ratio array sigma_hm_1^2/mu_hm_1^2             [-]
      sigma2_on_mu2_ip_2    ! Ratio array sigma_hm_2^2/mu_hm_2^2             [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      const_Ncnp2_on_Ncnm2,       & ! Prescribed ratio: <Ncn'^2> to <Ncn>^2  [-]
      stdev_const_Ncnp2_on_Ncnm2, & ! Standard deviation of const_Ncnp2_on_Ncnm2
      const_corr_chi_Ncn,         & ! Prescribed corr. of chi (old s) & Ncn  [-]
      const_corr_chi_Ncn_n_cloud    ! Prescribed corr. of chi & Ncn in cloud [-]

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      precip_frac_tol         ! Min. precip. frac. when hydromet. present    [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz,hydromet_dim) :: &
      wphydrometp_chnge    ! Change in wphydrometp_zt: covar. clip. [(m/s)units]

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      wm_zt,  & ! Mean vertical velocity, <w>, on thermo. levels  [m/s]
      wp2_zt    ! Variance of w, <w'^2> (interp. to t-levs.)      [m^2/s^2]

    real( kind = core_rknd ), dimension(nz) :: &
      rtp2_zt_from_chi
 
    real( kind = core_rknd ) :: &
      omicron,        & ! Relative width parameter, omicron = R / Rmax    [-]
      zeta_vrnce_rat    ! Width parameter for sigma_hm_1^2 / mu_hm_1^2    [-]

    logical :: l_corr_array_scaling

    ! Flags used for covariance clipping of <w'hm'>.
    logical, parameter :: &
      l_first_clip_ts = .true., & ! First instance of clipping in a timestep.
      l_last_clip_ts  = .true.    ! Last instance of clipping in a timestep.

    character(len=10) :: &
      hydromet_name    ! Name of a hydrometeor

    integer :: k, i, j  ! Loop indices

    ! ---- Begin Code ----
    
    ! Assertion check
    ! Check that all hydrometeors are positive otherwise exit the program
    if ( clubb_at_least_debug_level( 0 ) ) then
      if ( any(hydromet < zero_threshold) ) then
        do j = 1, hydromet_dim
          hydromet_name = hydromet_list(j)
          do k = 1, nz
            do i = 1, ngrdcol
              if ( hydromet(i,k,j) < zero_threshold ) then
                ! Write error message
                write(fstderr,*) " at beginning of setup_pdf_parameters: ", &
                                  trim( hydromet_name )//" = ", &
                                  hydromet(i,k,j), " < ", zero_threshold, &
                                  " at k = ", k
                err_code = clubb_fatal_error
                return
              end if ! hydromet(k,i) < 0
            end do ! k = 1, nz
          end do ! i = 1, hydromet_dim
        end do
      end if
    end if !clubb_at_least_debug_level( 0 )

    ! Setup some of the PDF parameters
    do k = 1, nz
      do i = 1, ngrdcol
        sigma_w_1(i,k)    = sqrt( pdf_params%varnce_w_1(i,k) )
        sigma_w_2(i,k)    = sqrt( pdf_params%varnce_w_2(i,k) )
      end do
    end do

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
      
      call precip_fraction( nz, ngrdcol,                                                  & ! In
                hydromet(:,:,:), cloud_frac(:,:), pdf_params%cloud_frac_1(:,:),           & ! In
                pdf_params%cloud_frac_2(:,:), ice_supersat_frac(:,:),                     & ! In
                pdf_params%ice_supersat_frac_1(:,:), pdf_params%ice_supersat_frac_2(:,:), & ! In
                pdf_params%mixt_frac(:,:), clubb_params, l_stats_samp,                    & ! In
                stats_sfc(:),                                                             & ! Inout
                precip_fracs%precip_frac(:,:),                                            & ! Out
                precip_fracs%precip_frac_1(:,:),                                          & ! Out
                precip_fracs%precip_frac_2(:,:),                                          & ! Out
                precip_frac_tol(:) )                                                        ! Out
      
      if ( err_code == clubb_fatal_error ) then
        write(fstderr,*) " in setup_pdf_parameters after calling precip_fraction"
        return
      end if

    else

      precip_fracs%precip_frac(:,:)   = one
      precip_fracs%precip_frac_1(:,:) = one
      precip_fracs%precip_frac_2(:,:) = one
      precip_frac_tol(:) = cloud_frac_min

    end if
    
    ! Calculate <N_cn> from Nc_in_cloud, whether Nc_in_cloud is predicted or
    ! based on a prescribed value, and whether the value is constant or varying
    ! over the grid level.
    if ( .not. l_const_Nc_in_cloud ) then
      ! Ncn varies at each vertical level.
      do k = 1, nz
        do i = 1, ngrdcol
         const_Ncnp2_on_Ncnm2(i,k) = Ncnp2_on_Ncnm2
        end do
      end do
       
      stdev_const_Ncnp2_on_Ncnm2(:,:) = stdev_L2N( const_Ncnp2_on_Ncnm2(:,:) )
       
    else  ! l_const_Nc_in_cloud
      ! Ncn is constant at each vertical level.
      const_Ncnp2_on_Ncnm2(:,:)       = zero
      stdev_const_Ncnp2_on_Ncnm2(:,:) = zero
    end if

    const_corr_chi_Ncn_n_cloud(:,:) = corr_array_n_cloud(iiPDF_Ncn,iiPDF_chi)
    
    call corr_NN2NL( nz, ngrdcol, & ! intent(in)
                     const_corr_chi_Ncn_n_cloud(:,:), & ! intent(in)
                     stdev_const_Ncnp2_on_Ncnm2(:,:), const_Ncnp2_on_Ncnm2(:,:), & ! intent(in)
                     const_corr_chi_Ncn(:,:)  ) ! intent(out)

    Ncnm(:,:) = Nc_in_cloud_to_Ncnm( &
                    pdf_params%chi_1(:,:), pdf_params%chi_2(:,:), pdf_params%stdev_chi_1(:,:), &
                    pdf_params%stdev_chi_2(:,:), pdf_params%mixt_frac(:,:), Nc_in_cloud(:,:), &
                    pdf_params%cloud_frac_1(:,:), pdf_params%cloud_frac_2(:,:), &
                    const_Ncnp2_on_Ncnm2(:,:), const_corr_chi_Ncn(:,:) )

    ! Boundary Condition.
    ! At thermodynamic level k = 1, which is below the model lower boundary, the
    ! value of Ncnm does not matter.
    Ncnm(:,1) = Nc_in_cloud(:,1)

    
    ! Calculate the overall variance of a precipitating hydrometeor (hm),
    !<hm'^2>.
    do j = 1, hydromet_dim
       do k = 1, nz, 1
         do i = 1, ngrdcol
             
            if ( hydromet(i,k,j) >= hydromet_tol(j) ) then
              ! There is some of the hydrometeor species found at level k.
              ! Calculate the variance (overall) of the hydrometeor.
              hydrometp2_zt(i,k,j) &
              = ( ( hmp2_ip_on_hmm2_ip(j) + one ) / precip_fracs%precip_frac(i,k) - one ) &
                * hydromet(i,k,j)**2
            else
              hydrometp2_zt(i,k,j) = zero
            end if
            
         end do ! k = 1, nz, 1
      end do
    end do
    
    ! Interpolate the overall variance of a hydrometeor, <hm'^2>, to its home on
    ! momentum grid levels.
    do i = 1, hydromet_dim, 1
      hydrometp2(:,:,i)  = zt2zm( nz, ngrdcol, gr, hydrometp2_zt(:,:,i) )
      hydrometp2(:,nz,i) = zero
    end do


    ! Calculate correlations involving w and Ncn by first calculating the
    ! overall covariance of w and Ncn using the down-gradient approximation.
    if ( l_calc_w_corr ) then
      
      ! Recalculate wm_zt and wp2_zt.  Mean vertical velocity may not be easy to
      ! pass into this subroutine from a host model, and wp2_zt needs to have a
      ! value consistent with the value it had when the PDF parameters involving w
      ! were originally set in subroutine pdf_closure.  The variable wp2 has since
      ! been advanced, resulting a new wp2_zt.  However, the value of wp2 here
      ! needs to be consistent with wp2 at the time the PDF parameters were
      ! calculated.

      ! Calculate the overall mean of vertical velocity, w, on thermodynamic
      ! levels.
      wm_zt(:,:) = compute_mean_binormal( pdf_params%w_1(:,:), pdf_params%w_2(:,:), &
                                          pdf_params%mixt_frac(:,:) )

      ! Calculate the overall variance of vertical velocity on thermodynamic
      ! levels.
      wp2_zt(:,:) = compute_variance_binormal( wm_zt(:,:), &
                                               pdf_params%w_1(:,:), pdf_params%w_2(:,:), &
                                               sigma_w_1(:,:), sigma_w_2(:,:), &
                                               pdf_params%mixt_frac(:,:) )
      
      ! Interpolate the covariances (overall) of w and precipitating
      ! hydrometeors to thermodynamic grid levels.
      do i = 1, hydromet_dim
        wphydrometp_zt(:,:,i) = zm2zt( nz, ngrdcol, gr, wphydrometp(:,:,i) )
      end do
          
      do j = 1, hydromet_dim
        do k = 1, nz
          do i = 1, ngrdcol
            if ( hydromet(i,k,j) < hydromet_tol(j) ) then
              wphydrometp_zt(i,k,j) = zero
            end if
          end do
        end do
      end do
      
      ! When the mean value of a precipitating hydrometeor is below tolerance
      ! value, it is considered to have a value of 0, and the precipitating
      ! hydrometeor does not vary over the grid level.  Any covariances
      ! involving that precipitating hydrometeor also have values of 0 at that
      ! grid level.
      do j = 1, hydromet_dim
        do k = 1, nz, 1
          do i = 1, ngrdcol

            ! Clip the value of covariance <w'hm'> on thermodynamic levels.
            call clip_covar_level( clip_wphydrometp, k, l_first_clip_ts,  & ! In
                                   l_last_clip_ts, dt, wp2_zt(i,k),       & ! In
                                   hydrometp2_zt(i,k,j),                  & ! In
                                   l_predict_upwp_vpwp,                   & ! In
                                   stats_zm(i),                           & ! Inout
                                   wphydrometp_zt(i,k,j),                 & ! Inout
                                   wphydrometp_chnge(i,k,j) )               ! Out

           end do ! k = 1, nz, 1
        end do ! i = 1, hydromet_dim, 1
      end do
      
      Kh_zm_c_K_hm(:,:) = -clubb_params(ic_K_hm) * Kh_zm(:,:)
      
      call calc_xpwp( nz, ngrdcol, gr, &
                      Kh_zm_c_K_hm, Ncnm, &
                      wpNcnp_zm )

      ! Boundary conditions; We are assuming zero flux at the top.
      wpNcnp_zm(:,nz) = zero

      ! Interpolate the covariances to thermodynamic grid levels.
      wpNcnp_zt(:,:) = zm2zt( nz, ngrdcol, gr, wpNcnp_zm(:,:) )

      ! When the mean value of Ncn is below tolerance value, it is considered
      ! to have a value of 0, and Ncn does not vary over the grid level.  Any
      ! covariance involving Ncn also has a value of 0 at that grid level.
      do k = 1, nz, 1
        do i = 1, ngrdcol
          if ( Ncnm(i,k) <= Ncn_tol ) then
            wpNcnp_zt(i,k) = zero
          end if
        end do ! k = 1, nz, 1
      end do
      
    end if ! l_calc_w_corr

    ! Unpack CLUBB parameters
    omicron = clubb_params(iomicron)
    zeta_vrnce_rat = clubb_params(izeta_vrnce_rat)

    !!! Calculate the means and standard deviations involving PDF variables
    !!! -- w, chi, eta, N_cn, and any precipitating hydrometeors (hm in-precip)
    !!! -- for each PDF component.
    call compute_mean_stdev( nz, ngrdcol, &
                             hydromet(:,:,:), hydrometp2_zt(:,:,:),         & ! Intent(in)
                             Ncnm(:,:), pdf_params%mixt_frac(:,:),          & ! Intent(in)
                             precip_fracs%precip_frac(:,:),                 & ! Intent(in)
                             precip_fracs%precip_frac_1(:,:),               & ! Intent(in)
                             precip_fracs%precip_frac_2(:,:),               & ! Intent(in)
                             precip_frac_tol(:),                            & ! Intent(in)
                             pdf_params%w_1(:,:),                           & ! Intent(in)
                             pdf_params%w_2(:,:),                           & ! Intent(in)
                             sigma_w_1(:,:),                                & ! Intent(in)
                             sigma_w_1(:,:),                                & ! Intent(in)
                             pdf_params%chi_1(:,:),                         & ! Intent(in)
                             pdf_params%chi_2(:,:),                         & ! Intent(in)
                             pdf_params%stdev_chi_1(:,:),                   & ! Intent(in)
                             pdf_params%stdev_chi_2(:,:),                   & ! Intent(in)
                             pdf_params%stdev_eta_1(:,:),                   & ! Intent(in)
                             pdf_params%stdev_eta_2(:,:),                   & ! Intent(in)
                             pdf_params%thl_1(:,:),                         & ! Intent(in)
                             pdf_params%thl_2(:,:),                         & ! Intent(in)
                             pdf_dim,                                       & ! Intent(in)
                             omicron, zeta_vrnce_rat,                       & ! Intent(in)
                             l_const_Nc_in_cloud,                           & ! Intent(in)
                             mu_x_1(:,:,:), mu_x_2(:,:,:),                  & ! Intent(out)
                             sigma_x_1(:,:,:), sigma_x_2(:,:,:),            & ! Intent(out)
                             hm_1(:,:,:), hm_2(:,:,:),                      & ! Intent(out)
                             sigma2_on_mu2_ip_1(:,:,:),                     & ! Intent(out)
                             sigma2_on_mu2_ip_2(:,:,:) )                      ! Intent(out)

    !!! Transform the component means and standard deviations involving
    !!! precipitating hydrometeors (hm in-precip) and N_cn -- ln hm and
    !!! ln N_cn -- to normal space for each PDF component.
    call norm_transform_mean_stdev( nz, ngrdcol,                            & ! Intent(in)
                                    hm_1(:,:,:), hm_2(:,:,:),               & ! Intent(in)
                                    Ncnm(:,:), pdf_dim,                     & ! Intent(in)
                                    mu_x_1(:,:,:), mu_x_2(:,:,:),           & ! Intent(in)
                                    sigma_x_1(:,:,:), sigma_x_2(:,:,:),     & ! Intent(in)
                                    sigma2_on_mu2_ip_1(:,:,:),              & ! Intent(in)
                                    sigma2_on_mu2_ip_2(:,:,:),              & ! Intent(in)
                                    l_const_Nc_in_cloud,                    & ! Intent(in)
                                    mu_x_1_n(:,:,:), mu_x_2_n(:,:,:),       & ! Intent(out)
                                    sigma_x_1_n(:,:,:), sigma_x_2_n(:,:,:) )  ! Intent(out)

    
    !!! Calculate the normal space correlations.
    !!! The normal space correlations are the same as the true correlations
    !!! except when at least one of the variables involved is a precipitating
    !!! hydrometeor or Ncn.  In these cases, the normal space correlation
    !!! involves the natural logarithm of the precipitating hydrometeors, ln hm
    !!! (for example, ln r_r and ln N_r), and ln N_cn for each PDF component.
    if ( l_diagnose_correlations ) then

      do k = 2, nz, 1
        do i = 1, ngrdcol

          if ( rcm(i,k) > rc_tol ) then

            call diagnose_correlations( pdf_dim, corr_array_n_cloud, & ! In
                                        l_calc_w_corr, &               ! In
                                        corr_array_1_n(i,k,:,:) )      ! Out

            call diagnose_correlations( pdf_dim, corr_array_n_cloud, & ! In
                                        l_calc_w_corr, &               ! In
                                        corr_array_2_n(i,k,:,:) )      ! Out

         else

            call diagnose_correlations( pdf_dim, corr_array_n_below, & ! In
                                        l_calc_w_corr, &               ! In
                                        corr_array_1_n(i,k,:,:) )      ! Out

            call diagnose_correlations( pdf_dim, corr_array_n_below, & ! In
                                        l_calc_w_corr, &               ! In
                                        corr_array_2_n(i,k,:,:) )      ! Out

          end if

        end do
      end do
      
      do k = 2, nz, 1
        do i = 1, ngrdcol
          call calc_cholesky_corr_mtx_approx &
                         ( pdf_dim, corr_array_1_n(i,k,:,:), &                      ! intent(in)
                           corr_cholesky_mtx_1(i,k,:,:), corr_array_1_n(i,k,:,:) )  ! intent(out)
        end do
      end do

      do k = 2, nz, 1
        do i = 1, ngrdcol
          call calc_cholesky_corr_mtx_approx &
                         ( pdf_dim, corr_array_2_n(i,k,:,:), &                      ! intent(in)
                           corr_cholesky_mtx_2(i,k,:,:), corr_array_2_n(i,k,:,:) )  ! intent(out)
        end do
      end do

    else ! if .not. l_diagnose_correlations
      
      if ( l_fix_w_chi_eta_correlations &
           .and. .not. l_calc_w_corr ) then
        
        ! When the flags are set this way, the correlation matrices do not vary with any vertical
        ! values, and instead are determined entirely by prescribed values. This results in there
        ! being only two unique correlation matrices, one for when the grid box is in cloud
        ! and one for when it is not. So instead of setting up correlation matrices for all
        ! grid boxes then calculating their Cholesky decomps, we can simply set up two correlation
        ! matrices, one for in cloud and one for out cloud, calculate the corresponding 
        ! Cholesky decompositions, then use the value of rc at each grid box to determine whether
        ! we assign the in cloud or out of cloud matrices to that grid box. 
        call calc_corr_norm_and_cholesky_factor( nz, ngrdcol, pdf_dim, iiPDF_type, & ! intent(in)
                                                 pdf_params%rc_1, pdf_params%rc_2, & ! intent(in)
                                                 corr_array_n_cloud, corr_array_n_below, &
                                                 corr_array_1_n, corr_array_2_n, & ! intent(out)
                                                 corr_cholesky_mtx_1, corr_cholesky_mtx_2 )!out
        
      else
        
        ! The correlation matrices can vary with vertical values. So we need to set the 
        ! correlation matrices up for each grid box, then find the Cholesky decomp for each
        ! grid box individually. This is very computationally expensive.

        call comp_corr_norm( nz, pdf_dim, ngrdcol, wm_zt(:,:),                        & ! In  
                             pdf_params%rc_1(:,:), pdf_params%rc_2(:,:),              & ! In  
                             pdf_params%mixt_frac(:,:),                               & ! In
                             precip_fracs%precip_frac_1(:,:),                         & ! In
                             precip_fracs%precip_frac_2(:,:),                         & ! In
                             wpNcnp_zt(:,:), wphydrometp_zt(:,:,:),                   & ! In
                             mu_x_1(:,:,:), mu_x_2(:,:,:),                            & ! In
                             sigma_x_1(:,:,:), sigma_x_2(:,:,:),                      & ! In
                             sigma_x_1_n(:,:,:), sigma_x_2_n(:,:,:),                  & ! In
                             corr_array_n_cloud, corr_array_n_below,                  & ! In
                             pdf_params,                                              & ! In
                             iiPDF_type,                                              & ! In
                             l_calc_w_corr,                                           & ! In
                             l_fix_w_chi_eta_correlations,                            & ! In
                             corr_array_1_n(:,:,:,:), corr_array_2_n(:,:,:,:) )         ! Out
                             
        ! Compute choleksy factorization for the correlation matrix of 1st PDF component
        do k = 2, nz, 1
          do i = 1, ngrdcol
            call Cholesky_factor( pdf_dim, corr_array_1_n(i,k,:,:),                 & ! In
                                  corr_array_scaling, corr_cholesky_mtx_1(i,k,:,:), & ! Out
                                  l_corr_array_scaling )                              ! Out
          end do
        end do
        
        ! Compute choleksy factorization for the correlation matrix of 2nd PDF component
        do k = 2, nz, 1
          do i = 1, ngrdcol
            call Cholesky_factor( pdf_dim, corr_array_2_n(i,k,:,:),                 & ! In
                                  corr_array_scaling, corr_cholesky_mtx_2(i,k,:,:), & ! Out
                                  l_corr_array_scaling )                              ! Out
          end do
        end do
        
      end if
      
    end if ! l_diagnose_correlations

    ! hydromet_pdf_params is optional, so if it is not present we simply skip the steps 
    ! to compute it.
    if ( present(hydromet_pdf_params) ) then

      !!! Calculate the true correlations for each PDF component.
      call denorm_transform_corr( nz, ngrdcol, pdf_dim,                                   & ! In
                                  sigma_x_1_n(:,:,:), sigma_x_2_n(:,:,:),                 & ! In
                                  sigma2_on_mu2_ip_1(:,:,:), sigma2_on_mu2_ip_2(:,:,:),   & ! In
                                  corr_array_1_n(:,:,:,:),                                & ! In
                                  corr_array_2_n(:,:,:,:),                                & ! In
                                  corr_array_1(:,:,:,:), corr_array_2(:,:,:,:) )            ! Out

      !!! Pack the PDF parameters
      call pack_hydromet_pdf_params( nz, ngrdcol, hm_1(:,:,:), hm_2(:,:,:), pdf_dim,  & ! In
                                     mu_x_1(:,:,:), mu_x_2(:,:,:),                    & ! In
                                     sigma_x_1(:,:,:), sigma_x_2(:,:,:),              & ! In
                                     corr_array_1(:,:,:,:), corr_array_2(:,:,:,:),    & ! In
                                     hydromet_pdf_params(:,:) )                         ! Out
                                     
      do i = 1, ngrdcol
        call init_hydromet_pdf_params( hydromet_pdf_params(i,1) ) ! intent(out)
      end do

    end if

    if ( l_stats_samp ) then
      
      do j = 1, hydromet_dim
        if ( ihmp2_zt(j) > 0 ) then
          ! Variance (overall) of the hydrometeor, <hm'^2>.
          ! call stat_update_var( ihmp2_zt(i), hydrometp2_zt(:,i), stats_zt )
          ! Switch back to using stat_update_var once the code is generalized
          ! to pass in the number of vertical levels.
            do k = 1, nz, 1
              do i = 1, ngrdcol
                call stat_update_var_pt( ihmp2_zt(j), k, hydrometp2_zt(i,k,j), & ! intent(in)
                                         stats_zt(i) ) ! intent(inout)
             end do ! k = 1, nz, 1
           end do
        end if
      end do
      
      do i = 1, ngrdcol
        call pdf_param_hm_stats( nz, pdf_dim, hm_1(i,:,:), hm_2(i,:,:), & ! intent(in)
                                 mu_x_1(i,:,:), mu_x_2(i,:,:), & ! intent(in)
                                 sigma_x_1(i,:,:), sigma_x_2(i,:,:), & ! intent(in)
                                 corr_array_1(i,:,:,:), corr_array_2(i,:,:,:), & ! intent(in)
                                 l_stats_samp, & ! intent(in)
                                 stats_zt(i) ) ! intent(inout)
      end do

      !!! Statistics for normal space PDF parameters involving hydrometeors.
      do i = 1, ngrdcol
        call pdf_param_ln_hm_stats( nz, pdf_dim, mu_x_1_n(i,:,:), & ! intent(in)
                                    mu_x_2_n(i,:,:), sigma_x_1_n(i,:,:), & ! intent(in)
                                    sigma_x_2_n(i,:,:), corr_array_1_n(i,:,:,:), & ! intent(in)
                                    corr_array_2_n(i,:,:,:), l_stats_samp, & ! intent(in)
                                    stats_zt(i) ) ! intent(inout)
      end do
      
      if ( irtp2_from_chi > 0 ) then

        do i = 1, ngrdcol
          rtp2_zt_from_chi &
            = compute_rtp2_from_chi( pdf_params%stdev_chi_1(i,:), pdf_params%stdev_chi_2(i,:), &
                                     pdf_params%stdev_eta_1(i,:), pdf_params%stdev_eta_2(i,:), &
                                     pdf_params%rt_1(i,:), pdf_params%rt_2(j,:),               &
                                     pdf_params%crt_1(i,:), pdf_params%crt_2(i,:),             &
                                     pdf_params%mixt_frac(i,:),                                &   
                                     corr_array_1_n(i,:,iiPDF_chi,iiPDF_eta),                  &
                                     corr_array_2_n(i,:,iiPDF_chi,iiPDF_eta) )

          ! Switch back to using stat_update_var once the code is generalized
          ! to pass in the number of vertical levels.
          ! call stat_update_var( irtp2_from_chi, zt2zm( rtp2_zt_from_chi ), &
          ! stats_zm )
          do k = 1, nz, 1
            call stat_update_var_pt( irtp2_from_chi, k, zt2zm( gr(i), rtp2_zt_from_chi, k ), & !in
                                     stats_zm(i) ) ! intent(inout)
          end do ! k = 1, nz, 1
        end do
        
      end if

      do i = 1, ngrdcol
        ! Switch back to using stat_update_var once the code is generalized
        ! to pass in the number of vertical levels.
        if ( iprecip_frac > 0 ) then
          ! Overall precipitation fraction.
          ! call stat_update_var( iprecip_frac, precip_frac, stats_zt )
          do k = 1, nz, 1
            call stat_update_var_pt( iprecip_frac, k, precip_fracs%precip_frac(i,k), & ! intent(in)
                                     stats_zt(i) ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

        if ( iprecip_frac_1 > 0 ) then
          ! Precipitation fraction in PDF component 1.
          ! call stat_update_var( iprecip_frac_1, precip_frac_1, stats_zt )
          do k = 1, nz, 1
            call stat_update_var_pt( iprecip_frac_1, k, precip_fracs%precip_frac_1(i,k), & ! In
                                     stats_zt(i) ) ! intent(inout)
          end do ! k = 1, nz, 1
        end if

        if ( iprecip_frac_2 > 0 ) then
          ! Precipitation fraction in PDF component 2.
          ! call stat_update_var( iprecip_frac_2, precip_frac_2, stats_zt )
          do k = 1, nz, 1
            call stat_update_var_pt( iprecip_frac_2, k, precip_fracs%precip_frac_2(i,k), & ! In
                                     stats_zt(i) ) ! intent(inout)
          end do ! k = 1, nz, 1
        end if

        if ( iNcnm > 0 ) then
          ! Mean simplified cloud nuclei concentration (overall).
          ! call stat_update_var( iNcnm, Ncnm, stats_zt )
          do k = 1, nz, 1
            call stat_update_var_pt( iNcnm, k, Ncnm(i,k), & ! intent(in)
                                     stats_zt(i) ) ! intent(inout)
          end do ! k = 1, nz, 1
        end if
      end do
      
    end if

    ! Boundary conditions for the output variables at k=1.
    mu_x_1_n(:,1,:) = zero
    mu_x_2_n(:,1,:) = zero
    sigma_x_1_n(:,1,:) = zero
    sigma_x_2_n(:,1,:) = zero
    corr_array_1_n(:,1,:,:) = zero
    corr_array_2_n(:,1,:,:) = zero
    corr_cholesky_mtx_1(:,1,:,:) = zero
    corr_cholesky_mtx_2(:,1,:,:) = zero
    
    precip_fracs%precip_frac(:,1)   = zero
    precip_fracs%precip_frac_1(:,1) = zero
    precip_fracs%precip_frac_2(:,1) = zero

    if ( clubb_at_least_debug_level( 2 ) ) then
      do k = 2, nz
        do i = 1, ngrdcol
          call assert_corr_symmetric( corr_array_1_n(i,k,:,:), pdf_dim) ! intent(in)
          call assert_corr_symmetric( corr_array_2_n(i,k,:,:), pdf_dim) ! intent(in)
        end do
      end do
    end if
    
    return

  end subroutine setup_pdf_parameters

  !=============================================================================
  subroutine compute_mean_stdev( nz, ngrdcol,                   & ! Intent(in)
                                 hydromet, hydrometp2_zt,       & ! Intent(in)
                                 Ncnm, mixt_frac,               & ! Intent(in)
                                 precip_frac,                   & ! Intent(in)
                                 precip_frac_1,                 & ! Intent(in)
                                 precip_frac_2,                 & ! Intent(in)
                                 precip_frac_tol,               & ! Intent(in)
                                 w_1, w_2,                      & ! Intent(in)
                                 stdev_w_1, stdev_w_2,          & ! Intent(in)
                                 chi_1, chi_2,                  & ! Intent(in)
                                 stdev_chi_1, stdev_chi_2,      & ! Intent(in)
                                 stdev_eta_1, stdev_eta_2,      & ! Intent(in)
                                 thl_1, thl_2,                  & ! Intent(in)
                                 pdf_dim,                       & ! Intent(in)
                                 omicron, zeta_vrnce_rat,       & ! Intent(in)
                                 l_const_Nc_in_cloud,           & ! Intent(in)
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
        zero     ! Constant(s)

    use array_index, only: &
        hydromet_tol, &
        iiPDF_chi,    & 
        iiPDF_eta,    &
        iiPDF_w,      &
        iiPDF_Ncn

    use index_mapping, only: &
        pdf2hydromet_idx  ! Procedure(s)

    use corr_varnce_module, only: &
        hmp2_ip_on_hmm2_ip, & ! Variable(s)
        Ncnp2_on_Ncnm2

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,      & ! Number of model vertical grid levels
      pdf_dim, & ! Number of PDF variables
      ngrdcol    ! Number of grid columns

    real( kind = core_rknd ), dimension(ngrdcol,nz,hydromet_dim), intent(in) :: &
      hydromet,       & ! Mean of a hydrometeor (overall)         [hm units]
      hydrometp2_zt     ! Variance of a hydrometeor (overall)     [(hm units)^2]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      Ncnm,          & ! Mean simplified cloud nuclei concentration   [num/kg]
      mixt_frac,     & ! Mixture fraction                             [-]
      precip_frac,   & ! Precipitation fraction (overall)             [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)   [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)   [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &  
      w_1,         & ! Mean of w (1st PDF component)                  [m/s]
      w_2,         & ! Mean of w (2nd PDF component)                  [m/s]
      stdev_w_1,  & ! Variance of w (1st PDF component)              [m^2/s^2]
      stdev_w_2,  & ! Variance of w (2nd PDF component)              [m^2/s^2]
      chi_1,       & ! Mean of chi (1st PDF component)                [kg/kg]
      chi_2,       & ! Mean of chi (2nd PDF component)                [kg/kg]
      stdev_chi_1, & ! Standard deviation of chi (1st PDF component)  [kg/kg]
      stdev_chi_2, & ! Standard deviation of chi (2nd PDF component)  [kg/kg]
      stdev_eta_1, & ! Standard deviation of eta (1st PDF component)  [kg/kg]
      stdev_eta_2, & ! Standard deviation of eta (2nd PDF component)  [kg/kg]
      thl_1,       & ! Mean of thl (1st PDF component)                [K]
      thl_2          ! Mean of thl (2nd PDF component)                [K]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      precip_frac_tol    ! Minimum precip. frac. when hydromet. are present  [-]

    real( kind = core_rknd ), intent(in) :: &
      omicron,        & ! Relative width parameter, omicron = R / Rmax    [-]
      zeta_vrnce_rat    ! Width parameter for sigma_hm_1^2 / mu_hm_1^2    [-]

    logical, intent(in) :: &
      l_const_Nc_in_cloud ! Use a constant cloud droplet conc. within cloud (K&K)

    ! Output Variables
    ! Note:  This code assumes to be these arrays in the same order as the
    ! correlation arrays, etc., which is determined by the iiPDF indices.
    ! The order should be as follows:  chi, eta, w, Ncn, <precip. hydrometeors>
    ! (indices increasing from left to right).
    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim), intent(out) :: &
      mu_x_1,    & ! Mean array of PDF vars. (1st PDF component)    [units vary]
      mu_x_2,    & ! Mean array of PDF vars. (2nd PDF component)    [units vary]
      sigma_x_1, & ! Standard deviation array of PDF vars (comp. 1) [units vary]
      sigma_x_2    ! Standard deviation array of PDF vars (comp. 2) [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nz,hydromet_dim), intent(out) :: &
      hm_1, & ! Mean of a precip. hydrometeor (1st PDF component)   [units vary]
      hm_2    ! Mean of a precip. hydrometeor (2nd PDF component)   [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim), intent(out) :: &
      sigma_hm_1_sqd_on_mu_hm_1_sqd, & ! Ratio sigma_hm_1^2 / mu_hm_1^2      [-]
      sigma_hm_2_sqd_on_mu_hm_2_sqd    ! Ratio sigma_hm_2^2 / mu_hm_2^2      [-]

    integer :: ivar ! Loop iterator

    integer :: hm_idx  ! Hydrometeor array index.


    !!! Initialize output variables.
    hm_1(:,:,:)       = zero
    hm_2(:,:,:)       = zero
    sigma_hm_1_sqd_on_mu_hm_1_sqd(:,:,:) = zero
    sigma_hm_2_sqd_on_mu_hm_2_sqd(:,:,:) = zero


    !!! Enter the PDF parameters.

    !!! Vertical velocity, w.

    ! Mean of vertical velocity, w, in PDF component 1.
    mu_x_1(:,:,iiPDF_w) = w_1(:,:)

    ! Mean of vertical velocity, w, in PDF component 2.
    mu_x_2(:,:,iiPDF_w) = w_2(:,:)

    ! Standard deviation of vertical velocity, w, in PDF component 1.
    sigma_x_1(:,:,iiPDF_w) = stdev_w_1(:,:)

    ! Standard deviation of vertical velocity, w, in PDF component 2.
    sigma_x_2(:,:,iiPDF_w) = stdev_w_2(:,:)


    !!! Extended liquid water mixing ratio, chi.

    ! Mean of extended liquid water mixing ratio, chi (old s),
    ! in PDF component 1.
    mu_x_1(:,:,iiPDF_chi) = chi_1(:,:)

    ! Mean of extended liquid water mixing ratio, chi (old s),
    ! in PDF component 2.
    mu_x_2(:,:,iiPDF_chi) = chi_2(:,:)

    ! Standard deviation of extended liquid water mixing ratio, chi (old s),
    ! in PDF component 1.
    sigma_x_1(:,:,iiPDF_chi) = stdev_chi_1(:,:)

    ! Standard deviation of extended liquid water mixing ratio, chi (old s),
    ! in PDF component 2.
    sigma_x_2(:,:,iiPDF_chi) = stdev_chi_2(:,:)


    !!! Coordinate orthogonal to chi, eta.

    ! Mean of eta (old t) in PDF component 1.
    ! Set the component mean values of eta to 0.
    ! The component mean values of eta are not important.  They can be set to
    ! anything.  They cancel out in the model code.  However, the best thing to
    ! do is to set them to 0 and avoid any kind of numerical error.
    mu_x_1(:,:,iiPDF_eta) = zero

    ! Mean of eta (old t) in PDF component 2.
    ! Set the component mean values of eta to 0.
    ! The component mean values of eta are not important.  They can be set to
    ! anything.  They cancel out in the model code.  However, the best thing to
    ! do is to set them to 0 and avoid any kind of numerical error.
    mu_x_2(:,:,iiPDF_eta) = zero

    ! Standard deviation of eta (old t) in PDF component 1.
    sigma_x_1(:,:,iiPDF_eta) = stdev_eta_1(:,:)

    ! Standard deviation of eta (old t) in PDF component 2.
    sigma_x_2(:,:,iiPDF_eta) = stdev_eta_2(:,:)


    !!! Simplified cloud nuclei concentration, Ncn.

    ! Mean of simplified cloud nuclei concentration, Ncn, in PDF component 1.
    mu_x_1(:,:,iiPDF_Ncn) = Ncnm(:,:)

    ! Mean of simplified cloud nuclei concentration, Ncn, in PDF component 2.
    mu_x_2(:,:,iiPDF_Ncn) = Ncnm(:,:)

    ! Standard deviation of simplified cloud nuclei concentration, Ncn,
    ! in PDF component 1.
    if ( .not. l_const_Nc_in_cloud ) then

      ! Ncn varies in both PDF components.
      sigma_x_1(:,:,iiPDF_Ncn) = sqrt( Ncnp2_on_Ncnm2 ) * Ncnm(:,:)

      sigma_x_2(:,:,iiPDF_Ncn) = sqrt( Ncnp2_on_Ncnm2 ) * Ncnm(:,:)

      ! Ncn is not an official hydrometeor.  However, both the
      ! sigma_hm_1_sqd_on_mu_hm_1_sqd and sigma_hm_2_sqd_on_mu_hm_2_sqd arrays
      ! have size pdf_dim, and both sigma_Ncn_1^2/mu_Ncn_1^2 and
      ! sigma_Ncn_2^2/mu_Ncn_2^2 need to be output as part of these arrays.
      sigma_hm_1_sqd_on_mu_hm_1_sqd(:,:,iiPDF_Ncn) = Ncnp2_on_Ncnm2
      sigma_hm_2_sqd_on_mu_hm_2_sqd(:,:,iiPDF_Ncn) = Ncnp2_on_Ncnm2

    else ! l_const_Nc_in_cloud

      ! Ncn is constant in both PDF components.
      sigma_x_1(:,:,iiPDF_Ncn) = zero

      sigma_x_2(:,:,iiPDF_Ncn) = zero

      ! Ncn is not an official hydrometeor.  However, both the
      ! sigma_hm_1_sqd_on_mu_hm_1_sqd and sigma_hm_2_sqd_on_mu_hm_2_sqd arrays
      ! have size pdf_dim, and both sigma_Ncn_1^2/mu_Ncn_1^2 and
      ! sigma_Ncn_2^2/mu_Ncn_2^2 need to be output as part of these arrays.
      sigma_hm_1_sqd_on_mu_hm_1_sqd(:,:,iiPDF_Ncn) = zero
      sigma_hm_2_sqd_on_mu_hm_2_sqd(:,:,iiPDF_Ncn) = zero

    end if ! .not. l_const_Nc_in_cloud


    !!! Precipitating hydrometeor species.
    do ivar = iiPDF_Ncn+1, pdf_dim, 1

      hm_idx = pdf2hydromet_idx(ivar)
       
      call calc_comp_mu_sigma_hm( nz, ngrdcol,                                       & ! In
                                  hydromet(:,:,hm_idx), hydrometp2_zt(:,:,hm_idx),   & ! In
                                  hmp2_ip_on_hmm2_ip(hm_idx),                        & ! In
                                  mixt_frac(:,:), precip_frac(:,:),                  & ! In
                                  precip_frac_1(:,:), precip_frac_2(:,:),            & ! In
                                  hydromet_tol(hm_idx), precip_frac_tol(:),          & ! In
                                  thl_1(:,:), thl_2(:,:),                            & ! In
                                  omicron, zeta_vrnce_rat,                           & ! In
                                  mu_x_1(:,:,ivar), mu_x_2(:,:,ivar),                & ! Out
                                  sigma_x_1(:,:,ivar), sigma_x_2(:,:,ivar),          & ! Out
                                  hm_1(:,:,hm_idx), hm_2(:,:,hm_idx),                & ! Out
                                  sigma_hm_1_sqd_on_mu_hm_1_sqd(:,:,ivar),           & ! Out
                                  sigma_hm_2_sqd_on_mu_hm_2_sqd(:,:,ivar) )            ! Out

    end do ! ivar = iiPDF_Ncn+1, pdf_dim, 1


    return

  end subroutine compute_mean_stdev
  
  !=============================================================================
  subroutine calc_corr_norm_and_cholesky_factor( nz, ngrdcol, pdf_dim, iiPDF_type, &
                                                 rc_1, rc_2, &
                                                 corr_array_n_cloud, corr_array_n_below, &
                                                 corr_array_1_n, corr_array_2_n, &
                                                 corr_cholesky_mtx_1, corr_cholesky_mtx_2 )

    ! Description: This subroutine computes the correlation arrays and correlation
    !   Cholesky matrices of PDF vars for both components. Here, we assume that
    !   there are only two unique correlation arrays, which allows us to compute 
    !   these two unique arrays and their corresponding Cholesky decompositions,
    !   then use rc to determine which one to assign to each grid column and 
    !   vertical level. If the correlation arrays vary based on vertically varying
    !   values, then this subroutine is not appropriate.
    !   
    ! References:
    !   https://github.com/larson-group/cam/issues/129#issuecomment-816205563
    !-----------------------------------------------------------------------
    
    use constants_clubb, only:  &
        rc_tol,      &
        zero

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use array_index, only: &
        iiPDF_chi, & ! Variable(s)
        iiPDF_eta, &
        iiPDF_w,   &
        iiPDF_Ncn
        
    use model_flags, only: &
        iiPDF_ADG1,       & ! Variable(s)
        iiPDF_ADG2,       &
        iiPDF_new_hybrid

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s)    
        
    use matrix_operations, only: &
        Cholesky_factor ! Procedure(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,      & ! Number of vertical levels
      pdf_dim, & ! Number of variables in the corr/mean/stdev arrays
      ngrdcol    ! Number of grid columns

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      rc_1,   & ! Mean of r_c (1st PDF component)                 [kg/kg]
      rc_2      ! Mean of r_c (2nd PDF component)                 [kg/kg]

    real( kind = core_rknd ), dimension(pdf_dim,pdf_dim), intent(in) :: &
      corr_array_n_cloud, & ! Prescribed correlation array in cloud        [-]
      corr_array_n_below    ! Prescribed correlation array below cloud     [-]

    integer, intent(in) :: &
      iiPDF_type    ! Selected option for the two-component normal (double
                    ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                    ! w, chi, and eta) portion of CLUBB's multivariate,
                    ! two-component PDF.
                    
    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim,pdf_dim), intent(out) :: &
      corr_array_1_n,       & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
      corr_array_2_n,       & ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]
      corr_cholesky_mtx_1,  & ! Transposed corr. cholesky matrix, 1st comp.        [-]
      corr_cholesky_mtx_2     ! Transposed corr. cholesky matrix, 2nd comp.        [-]

    ! Local Variables
    real( kind = core_rknd ), dimension(pdf_dim,pdf_dim) :: &
      corr_array_cloud,         & ! General in cloud corr. matrix
      corr_array_below,         & ! General out of cloud corr. matrix
      corr_cholesky_mtx_cloud,  & ! General in cloud Cholesky matrix
      corr_cholesky_mtx_below     ! General out of cloud Cholesky matrix
      
    logical :: &
      l_corr_array_scaling  ! Dummy variable that we need for calling Cholesky_factor
      
    real( kind = core_rknd ), dimension(pdf_dim) :: &
      corr_array_scaling    ! Dummy variable that we need for calling Cholesky_factor

    integer :: ivar, jvar, i, k ! Indices
                                    
    !-------------------- Begin Code --------------------
    
    ! Initialize correlation arrays with prescribed values
    corr_array_cloud(:,:) = corr_array_n_cloud(:,:)
    corr_array_below(:,:) = corr_array_n_below(:,:)
    
    ! The ADG1 PDF fixes the correlation of w and rt and the correlation of
    ! w and theta_l to be 0, which means the correlation of w and chi and the
    ! correlation of w and eta must also be 0.
    if ( ( iiPDF_type == iiPDF_ADG1 .or. iiPDF_type == iiPDF_ADG2 &
            .or. iiPDF_type == iiPDF_new_hybrid ) &
          .and. l_follow_ADG1_PDF_standards ) then
          
      corr_array_cloud(iiPDF_w,iiPDF_chi) = zero
      corr_array_below(iiPDF_w,iiPDF_chi) = zero
      
      corr_array_cloud(iiPDF_w,iiPDF_eta) = zero
      corr_array_below(iiPDF_w,iiPDF_eta) = zero
      
    end if
    
    ! Ncn is an inherently in-cloud property, so replace out of cloud correlation values
    ! with in cloud ones.
    corr_array_below(iiPDF_Ncn,iiPDF_chi) = corr_array_cloud(iiPDF_Ncn,iiPDF_chi)
    corr_array_below(iiPDF_Ncn,iiPDF_eta) = corr_array_cloud(iiPDF_Ncn,iiPDF_eta)
    
    
    ! Estimates the correlation of the natural logarithm of a
    ! hydrometeor species and eta using the correlation of chi and eta and the
    ! correlation of chi and the natural logarithm of the hydrometeor.  This
    ! facilitates the Cholesky decomposability of the correlation array that will
    ! inevitably be decomposed for SILHS purposes. Without this estimation, we
    ! have found that the resulting correlation matrix cannot be decomposed.
    do jvar = iiPDF_Ncn+1, pdf_dim
      
      corr_array_cloud(jvar,iiPDF_eta) = corr_array_cloud(iiPDF_eta,iiPDF_chi) &
                                         * corr_array_cloud(jvar,iiPDF_chi)
      
      corr_array_below(jvar,iiPDF_eta) = corr_array_below(iiPDF_eta,iiPDF_chi) &
                                         * corr_array_below(jvar,iiPDF_chi)
      
    end do
    
    ! Calc in cloud Cholesky 
    call Cholesky_factor( pdf_dim, corr_array_cloud(:,:), & ! In
                          corr_array_scaling, corr_cholesky_mtx_cloud(:,:), &  ! Out
                          l_corr_array_scaling ) ! Out
                
    ! Calc out of cloud Cholesky           
    call Cholesky_factor( pdf_dim, corr_array_below(:,:), & ! In
                          corr_array_scaling, corr_cholesky_mtx_below(:,:), &  ! Out
                          l_corr_array_scaling ) ! Out
                          
    ! Use rc_1 to determine which correlation and Cholesky matrices to assign to 1st PDF 
    ! Correlation matrices are symmetric, so we copy ij values to ji indices as well
    ! Cholesky matrices are lower triangular, so we only copy those values 
    do ivar = 1, pdf_dim
      do jvar = ivar, pdf_dim
        do k = 1, nz
          do i = 1, ngrdcol
            
            if ( rc_1(i,k) > rc_tol ) then
              ! Assign in cloud matrices to 1st PDF component
              corr_array_1_n(i,k,jvar,ivar)      = corr_array_cloud(jvar,ivar)
              corr_array_1_n(i,k,ivar,jvar)      = corr_array_cloud(jvar,ivar)
              corr_cholesky_mtx_1(i,k,jvar,ivar) = corr_cholesky_mtx_cloud(jvar,ivar)
            else
              ! Assign out of cloud matrices to 1st PDF component
              corr_array_1_n(i,k,jvar,ivar)      = corr_array_below(jvar,ivar)
              corr_array_1_n(i,k,ivar,jvar)      = corr_array_below(jvar,ivar)
              corr_cholesky_mtx_1(i,k,jvar,ivar) = corr_cholesky_mtx_below(jvar,ivar)
            end if
            
          end do
        end do
      end do
    end do
            
    ! Use rc_1 to determine which correlation and Cholesky matrices to assign to 2nd PDF 
    ! Correlation matrices are symmetric, so we copy ij values to ji indices as well
    ! Cholesky matrices are lower triangular, so we only copy those values 
    do ivar = 1, pdf_dim
      do jvar = ivar, pdf_dim
        do k = 1, nz
          do i = 1, ngrdcol
            
            if ( rc_2(i,k) > rc_tol ) then
              ! Assign in cloud matrices to 2nd PDF component
              corr_array_2_n(i,k,jvar,ivar)      = corr_array_cloud(jvar,ivar)
              corr_array_2_n(i,k,ivar,jvar)      = corr_array_cloud(jvar,ivar)
              corr_cholesky_mtx_2(i,k,jvar,ivar) = corr_cholesky_mtx_cloud(jvar,ivar)
            else
              ! Assign out of cloud matrices to 2nd PDF component
              corr_array_2_n(i,k,jvar,ivar)      = corr_array_below(jvar,ivar)
              corr_array_2_n(i,k,ivar,jvar)      = corr_array_below(jvar,ivar)
              corr_cholesky_mtx_2(i,k,jvar,ivar) = corr_cholesky_mtx_below(jvar,ivar)
            end if
            
          end do
        end do
      end do
    end do
    
    ! Set upper triangular parts of Cholesky matrices to zero
    do ivar = 1, pdf_dim
      do jvar = 1, ivar-1
        corr_cholesky_mtx_1(:,:,jvar,ivar) = zero
        corr_cholesky_mtx_2(:,:,jvar,ivar) = zero
      end do
    end do
  
  end subroutine calc_corr_norm_and_cholesky_factor
  
  !=============================================================================
  subroutine comp_corr_norm( nz, pdf_dim, ngrdcol, wm_zt, rc_1, rc_2, &
                             mixt_frac, &
                             precip_frac_1, &
                             precip_frac_2, &
                             wpNcnp_zt, wphydrometp_zt, &
                             mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &
                             sigma_x_1_n, sigma_x_2_n, &
                             corr_array_n_cloud, corr_array_n_below, &
                             pdf_params, &
                             iiPDF_type, &
                             l_calc_w_corr, &
                             l_fix_w_chi_eta_correlations, &
                             corr_array_1_n, corr_array_2_n )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        Ncn_tol,      &
        one,          &
        zero

    use diagnose_correlations_module, only: &
        calc_mean,        & ! Procedure(s)
        calc_w_corr

    use index_mapping, only: &
        pdf2hydromet_idx  ! Procedure(s)

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use array_index, only: &
        iiPDF_chi, & ! Variable(s)
        iiPDF_eta, &
        iiPDF_w,   &
        iiPDF_Ncn, &
        hydromet_tol

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s)    
        
    use matrix_operations, only: &
        mirror_lower_triangular_matrix

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,      & ! Number of vertical levels
      pdf_dim, & ! Number of variables in the corr/mean/stdev arrays
      ngrdcol    ! Number of grid columns

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      wm_zt,         & ! Mean vertical velocity, <w>, on thermo. levels    [m/s]
      rc_1,          & ! Mean of r_c (1st PDF component)                 [kg/kg]
      rc_2,          & ! Mean of r_c (2nd PDF component)                 [kg/kg]
      mixt_frac,     & ! Mixture fraction                                    [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2, & ! Precipitation fraction (2nd PDF component)          [-]
      wpNcnp_zt        ! Covariance of w and N_cn on t-levs.      [(m/s) num/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nz,hydromet_dim), intent(in) :: &
      wphydrometp_zt    ! Covariance of w and hm interp. to t-levs.  [(m/s)u.v.]

    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim), intent(in) :: &
      mu_x_1,      & ! Mean of x array (1st PDF component)          [units vary]
      mu_x_2,      & ! Mean of x array (2nd PDF component)          [units vary]
      sigma_x_1,   & ! Standard deviation of x array (1st PDF comp.)  [un. vary]
      sigma_x_2,   & ! Standard deviation of x array (2nd PDF comp.)  [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(pdf_dim,pdf_dim), &
    intent(in) :: &
      corr_array_n_cloud, & ! Prescribed correlation array in cloud        [-]
      corr_array_n_below    ! Prescribed correlation array below cloud     [-]
      
    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                               [units vary]

    integer, intent(in) :: &
      iiPDF_type    ! Selected option for the two-component normal (double
                    ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                    ! w, chi, and eta) portion of CLUBB's multivariate,
                    ! two-component PDF.

    logical, intent(in) :: &
      l_calc_w_corr, &                ! Calculate the correlations between w and the hydrometeors
      l_fix_w_chi_eta_correlations    ! Use a fixed correlation for s and t Mellor(chi/eta)

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim,pdf_dim), &
    intent(out) :: &
      corr_array_1_n, & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
      corr_array_2_n    ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]

    ! Local Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim)  :: &
      corr_w_hm_1_n, & ! Correlation of w and ln hm (1st PDF component) ip   [-]
      corr_w_hm_2_n    ! Correlation of w and ln hm (2nd PDF component) ip   [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      corr_w_Ncn_1_n, & ! Correlation of w and ln Ncn (1st PDF component)    [-]
      corr_w_Ncn_2_n    ! Correlation of w and ln Ncn (2nd PDF component)    [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
       ones_vector,     & ! Vector of 1s
       Ncn_tol_in,      & ! Tolerance value for Ncn
       hydromet_tol_in    ! Tolerance value for hydromet

    logical :: &
      l_limit_corr_chi_eta    ! Flag to limit the correlation of chi and eta [-]

    integer :: ivar, jvar, hm_idx, i, k ! Indices

    ! ---- Begin Code ----

    
    !!! Normal space correlations

    ! Initialize corr_w_hm_1_n and corr_w_hm_2_n arrays to 0.
    corr_w_hm_1_n(:,:,:) = zero
    corr_w_hm_2_n(:,:,:) = zero


    ! Set ones_vector to a vector of 1s.
    ones_vector = one

    ! Calculate normal space correlations involving w by first calculating total
    ! covariances involving w (<w'Ncn'>, etc.) using the down-gradient
    ! approximation.
    if ( l_calc_w_corr ) then

      Ncn_tol_in = Ncn_tol

      ! Calculate the correlation of w and ln Ncn in each PDF component.
      ! The subroutine calc_corr_w_hm_n can be used to do this as long as a
      ! value of 1 is sent in for precip_frac_1 and precip_frac_2.
      jvar = iiPDF_Ncn
      do i = 1, ngrdcol
        call calc_corr_w_hm_n( wm_zt(i,:), wpNcnp_zt(i,:), &
                               mu_x_1(i,:,iiPDF_w), mu_x_2(i,:,iiPDF_w), &
                               mu_x_1(i,:,jvar), mu_x_2(i,:,jvar), &
                               sigma_x_1(i,:,iiPDF_w), sigma_x_2(i,:,iiPDF_w), &
                               sigma_x_1(i,:,jvar), sigma_x_2(i,:,jvar), &
                               sigma_x_1_n(i,:,jvar), sigma_x_2_n(i,:,jvar), &
                               mixt_frac(i,:), ones_vector(i,:), &
                               ones_vector(i,:), Ncn_tol_in(i,:), &
                               corr_w_Ncn_1_n(i,:), corr_w_Ncn_2_n(i,:) )
      end do

      ! Calculate the correlation of w and the natural logarithm of the
      ! hydrometeor for each PDF component and each hydrometeor type.
      do jvar = iiPDF_Ncn+1, pdf_dim

        hm_idx = pdf2hydromet_idx(jvar)

        hydromet_tol_in(:,:) = hydromet_tol(hm_idx)

        do i = 1, ngrdcol
          call calc_corr_w_hm_n( wm_zt(i,:), wphydrometp_zt(i,:,hm_idx), & ! intent(in)
                                 mu_x_1(i,:,iiPDF_w), mu_x_2(i,:,iiPDF_w), & ! intent(in)
                                 mu_x_1(i,:,jvar), mu_x_2(i,:,jvar), & ! intent(in)
                                 sigma_x_1(i,:,iiPDF_w), sigma_x_2(i,:,iiPDF_w), & ! intent(in)
                                 sigma_x_1(i,:,jvar), sigma_x_2(i,:,jvar), & ! intent(in)
                                 sigma_x_1_n(i,:,jvar), sigma_x_2_n(i,:,jvar), & ! intent(in)
                                 mixt_frac(i,:), precip_frac_1(i,:), & ! intent(in)
                                 precip_frac_2(i,:), hydromet_tol_in(i,:), & ! intent(in)
                                 corr_w_hm_1_n(i,:,jvar), corr_w_hm_2_n(i,:,jvar) ) ! intent(out)
        end do

      end do ! jvar = iiPDF_Ncn+1, pdf_dim

    end if ! l_calc_w_corr

    ! In order to decompose the normal space correlation matrix,
    ! we must not have a perfect correlation of chi and
    ! eta. Thus, we impose a limitation.
    l_limit_corr_chi_eta = .true.


    ! Initialize the normal space correlation arrays
    corr_array_1_n(:,:,:,:) = zero
    corr_array_2_n(:,:,:,:) = zero

    !!! The corr_arrays are assumed to be lower triangular matrices
    ! Set diagonal elements to 1
    do ivar=1, pdf_dim
      do k = 1, nz
        do i = 1, ngrdcol
          corr_array_1_n(i,k,ivar,ivar) = one
          corr_array_2_n(i,k,ivar,ivar) = one
        end do
      end do
    end do


    !!! This code assumes the following order in the prescribed correlation
    !!! arrays (iiPDF indices):
    !!! chi, eta, w, Ncn, <hydrometeors> (indices increasing from left to right)
    if ( l_fix_w_chi_eta_correlations ) then
      ! Correlation of chi (old s) and eta (old t)
      call component_corr_chi_eta( nz, ngrdcol, & ! intent(in)
                                  rc_1(:,:), & ! intent(in)
                                  rc_2(:,:), & ! intent(in)
                                  corr_array_n_cloud(iiPDF_eta,iiPDF_chi), & ! intent(in)
                                  corr_array_n_below(iiPDF_eta,iiPDF_chi), & ! intent(in)
                                  l_limit_corr_chi_eta, & ! intent(in)
                                  corr_array_1_n(:,:,iiPDF_eta,iiPDF_chi), & ! intent(out)
                                  corr_array_2_n(:,:,iiPDF_eta,iiPDF_chi) ) ! intent(out)
    else
      
      ! Preferred, more accurate version.
      do i = 1, ngrdcol
        corr_array_1_n(i,:,iiPDF_eta,iiPDF_chi) = pdf_params%corr_chi_eta_1(i,:)
        corr_array_2_n(i,:,iiPDF_eta,iiPDF_chi) = pdf_params%corr_chi_eta_2(i,:)
      end do
      
    end if

    if ( l_fix_w_chi_eta_correlations ) then
      ! Correlation of chi (old s) and w
      call component_corr_w_x( nz, ngrdcol, & ! intent(in)
                               rc_1(:,:), & ! intent(in)
                               rc_2(:,:), & ! intent(in)
                               corr_array_n_cloud(iiPDF_w,iiPDF_chi), & ! intent(in)
                               corr_array_n_below(iiPDF_w,iiPDF_chi), & ! intent(in)
                               iiPDF_type, & ! intent(in)
                               corr_array_1_n(:,:,iiPDF_w,iiPDF_chi), & ! intent(out)
                               corr_array_2_n(:,:,iiPDF_w,iiPDF_chi) ) ! intent(out)
    else
      
      ! Preferred, more accurate version.
      do i = 1, ngrdcol
        corr_array_1_n(i,:,iiPDF_w,iiPDF_chi) = pdf_params%corr_w_chi_1(i,:)
        corr_array_2_n(i,:,iiPDF_w,iiPDF_chi) = pdf_params%corr_w_chi_2(i,:)
      end do
      
    end if 
                             

    ! Correlation of chi (old s) and ln Ncn, corr_array_n_cloud used twice because 
    ! Ncn is an inherently in-cloud property.
    call component_corr_x_hm_n_ip( nz, ngrdcol, & ! intent(in)
                                   rc_1(:,:), & ! intent(in)
                                   rc_2(:,:), & ! intent(in)
                                   corr_array_n_cloud(iiPDF_Ncn,iiPDF_chi), & ! intent(in)
                                   corr_array_n_cloud(iiPDF_Ncn,iiPDF_chi), & ! intent(in)
                                   corr_array_1_n(:,:,iiPDF_Ncn,iiPDF_chi), & ! intent(out)
                                   corr_array_2_n(:,:,iiPDF_Ncn,iiPDF_chi) ) ! intent(out)

    ! Correlation of chi (old s) and the natural logarithm of the hydrometeors
    do jvar = iiPDF_Ncn+1, pdf_dim
      call component_corr_x_hm_n_ip( nz, ngrdcol, & ! intent(in)
                                     rc_1(:,:), & ! intent(in)
                                     rc_2(:,:), & ! intent(in)
                                     corr_array_n_cloud(jvar,iiPDF_chi), & ! intent(in)
                                     corr_array_n_below(jvar,iiPDF_chi), & ! intent(in)
                                     corr_array_1_n(:,:,jvar,iiPDF_chi), & ! intent(out)
                                     corr_array_2_n(:,:,jvar,iiPDF_chi) ) ! intent(out)
    end do

    ! Correlation of eta (old t) and w
    if ( l_fix_w_chi_eta_correlations ) then
      ! Correlation of chi (old s) and w
      call component_corr_w_x( nz, ngrdcol, & ! intent(in)
                               rc_1(:,:), & ! intent(in)
                               rc_2(:,:), & ! intent(in)
                               corr_array_n_cloud(iiPDF_w,iiPDF_eta), & ! intent(in)
                               corr_array_n_below(iiPDF_w,iiPDF_eta), & ! intent(in)
                               iiPDF_type, & ! intent(in)
                               corr_array_1_n(:,:,iiPDF_w,iiPDF_eta), & ! intent(out)
                               corr_array_2_n(:,:,iiPDF_w,iiPDF_eta) ) ! intent(out)
    else
      
      ! Preferred, more accurate version.
      do i = 1, ngrdcol
        corr_array_1_n(i,:,iiPDF_w,iiPDF_chi) = pdf_params%corr_w_chi_1(i,:)
        corr_array_2_n(i,:,iiPDF_w,iiPDF_chi) = pdf_params%corr_w_chi_2(i,:)
      end do
      
    end if 
    
    ! Correlation of eta (old t) and ln Ncn, corr_array_n_cloud used twice because 
    ! Ncn is an inherently in-cloud property.
    call component_corr_x_hm_n_ip( nz, ngrdcol, & ! intent(in)
                                   rc_1(:,:), & ! intent(in)
                                   rc_2(:,:), & ! intent(in)
                                   corr_array_n_cloud(iiPDF_Ncn,iiPDF_eta), & ! intent(in)
                                   corr_array_n_cloud(iiPDF_Ncn,iiPDF_eta), & ! intent(in)
                                   corr_array_1_n(:,:,iiPDF_Ncn,iiPDF_eta), & ! intent(out)
                                   corr_array_2_n(:,:,iiPDF_Ncn,iiPDF_eta) ) ! intent(out)
                                   

    ! Correlation of eta (old t) and the natural logarithm of the hydrometeors
    do jvar = iiPDF_Ncn+1, pdf_dim
      call component_corr_eta_hm_n_ip( nz, ngrdcol, & ! intent(in)
                                       corr_array_1_n(:,:,iiPDF_eta,iiPDF_chi), & ! intent(in)
                                       corr_array_1_n(:,:,jvar,iiPDF_chi), & ! intent(in)
                                       corr_array_2_n(:,:,iiPDF_eta,iiPDF_chi), & ! intent(in)
                                       corr_array_2_n(:,:,jvar,iiPDF_chi), & ! intent(in)
                                       corr_array_1_n(:,:,jvar,iiPDF_eta), & ! intent(out)
                                       corr_array_2_n(:,:,jvar,iiPDF_eta) ) ! intent(out)
    end do
    

    ! Correlation of w and ln Ncn
    call component_corr_w_hm_n_ip( nz, ngrdcol, & ! intent(in)
                                   corr_w_Ncn_1_n(:,:), rc_1(:,:), & ! intent(in)
                                   corr_w_Ncn_2_n(:,:), rc_2(:,:), & ! intent(in)
                                   corr_array_n_cloud(iiPDF_Ncn,iiPDF_w), & ! intent(in)
                                   corr_array_n_below(iiPDF_Ncn,iiPDF_w), & ! intent(in)
                                   l_calc_w_corr, & ! intent(in)
                                   corr_array_1_n(:,:,iiPDF_Ncn,iiPDF_w), & ! intent(out)
                                   corr_array_2_n(:,:,iiPDF_Ncn,iiPDF_w) ) ! intent(out)

    ! Correlation of w and the natural logarithm of the hydrometeors
    do jvar = iiPDF_Ncn+1, pdf_dim
      call component_corr_w_hm_n_ip( nz, ngrdcol, & ! intent(in)
                                     corr_w_hm_1_n(:,:,jvar), rc_1(:,:), & ! intent(in)
                                     corr_w_hm_2_n(:,:,jvar), rc_2(:,:), & ! intent(in)
                                     corr_array_n_cloud(jvar,iiPDF_w), & ! intent(in)
                                     corr_array_n_below(jvar,iiPDF_w), & ! intent(in)
                                     l_calc_w_corr, & ! intent(in)
                                     corr_array_1_n(:,:,jvar,iiPDF_w), & ! intent(out)
                                     corr_array_2_n(:,:,jvar,iiPDF_w) ) ! intent(out)
    end do

    ! Correlation of ln Ncn and the natural logarithm of the hydrometeors
    do jvar = iiPDF_Ncn+1, pdf_dim
      
      call component_corr_hmx_hmy_n_ip( nz, ngrdcol, & ! intent(in)
                                        rc_1(:,:), & ! intent(in)
                                        rc_2(:,:), & ! intent(in)
                                        corr_array_n_cloud(jvar,iiPDF_Ncn), & ! intent(in)
                                        corr_array_n_below(jvar,iiPDF_Ncn), & ! intent(in)
                                        corr_array_1_n(:,:,jvar,iiPDF_Ncn), & ! intent(out)
                                        corr_array_2_n(:,:,jvar,iiPDF_Ncn) ) ! intent(out)
    end do

    ! Correlation of the natural logarithm of two hydrometeors
    do ivar = iiPDF_Ncn+1, pdf_dim-1
      do jvar = ivar+1, pdf_dim


        call component_corr_hmx_hmy_n_ip( nz, ngrdcol, & ! intent(in)
                                          rc_1(:,:), & ! intent(in)
                                          rc_2(:,:), & ! intent(in)
                                          corr_array_n_cloud(jvar,ivar), & ! intent(in)
                                          corr_array_n_below(jvar,ivar), & ! intent(in)
                                          corr_array_1_n(:,:,jvar,ivar), & ! intent(out)
                                          corr_array_2_n(:,:,jvar,ivar) ) ! intent(out)

       end do ! jvar
    end do ! ivar
    
    ! For ease of use later in the code, we make the correlation arrays
    ! symmetrical
    do k = 2, nz, 1
      do i = 1, ngrdcol
        call mirror_lower_triangular_matrix( pdf_dim, corr_array_1_n(i,k,:,:) )
        call mirror_lower_triangular_matrix( pdf_dim, corr_array_2_n(i,k,:,:) )
      end do
    end do

    return

  end subroutine comp_corr_norm
  
  !=============================================================================
  subroutine calc_comp_mu_sigma_hm( nz, ngrdcol, &                   ! In
                                    hmm, hmp2, &                     ! In
                                    hmp2_ip_on_hmm2_ip, &            ! In
                                    mixt_frac, precip_frac, &        ! In
                                    precip_frac_1, precip_frac_2, &  ! In
                                    hm_tol, precip_frac_tol, &       ! In
                                    mu_thl_1, mu_thl_2, &            ! In
                                    omicron, zeta_vrnce_rat_in, &    ! In
                                    mu_hm_1, mu_hm_2, &              ! Out
                                    sigma_hm_1, sigma_hm_2, &        ! Out
                                    hm_1, hm_2, &                    ! Out
                                    sigma_hm_1_sqd_on_mu_hm_1_sqd, & ! Out
                                    sigma_hm_2_sqd_on_mu_hm_2_sqd )  ! Out


    ! Description:
    ! Calculates the in-precipitation mean and in-precipitation standard
    ! deviation in both PDF components for a precipitating hydrometeor.
    !
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
    !
    !
    ! DESCRIPTION OF THE METHOD THAT SOLVES FOR TWO IN-PRECIPITATION COMPONENTS
    ! =========================================================================
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
    !
    !
    ! Update:
    !
    ! In order to better represent the increase in <th_l'^2> near the ground
    ! from the evaporation of rain, the code is modified to tend towards a
    ! negative correlation of th_l and hm.  When mu_thl_1 <= mu_thl_2,
    ! mu_hm_1 >= mu_hm_2; otherwise, when mu_thl_1 > mu_thl_2,
    ! mu_hm_1 <= mu_hm_2.
    !
    ! In the original derivation, mu_hm_1 >= mu_hm_2 when zeta >= 0, and
    ! mu_hm_1 < mu_hm_2 when zeta < 0, where zeta is a tunable or adjustable
    ! parameter.  In order to allow the relationship of mu_hm_1 to mu_hm_2 to
    ! depend on the relationship of mu_thl_1 to mu_thl_2, the value of zeta must
    ! also depend on the relationship of mu_thl_1 to mu_thl_2.
    !
    ! When mu_thl_1 <= mu_thl_2:
    !
    ! The relationship of mu_hm_1 to mu_hm_2 is mu_hm_1 >= mu_hm_2, so
    ! zeta >= 0.  The tunable value of zeta is referred to as zeta_in.  When
    ! zeta_in is already greater than 0 (meaning sigma_hm_1^2 / mu_hm_1^2 is
    ! greater than sigma_hm_2^2 / mu_hm_2^2), zeta is simply set to zeta_in.  In
    ! other words, the component with the greater mean value of the hydrometeor
    ! also has the greater value of component variance to the square of the
    ! component mean.  However, when zeta_in is less than 0, zeta must be
    ! adjusted to be greater than 0.  The following equation is used to set the
    ! value of zeta:
    !
    !        | zeta_in; when zeta_in >= 0
    ! zeta = |
    !        | ( 1 / ( 1 + zeta_in ) ) - 1; when zeta_in < 0.
    !
    ! Previously, when zeta_in < 0, mu_hm_1 < mu_hm_2, and
    ! sigma_hm_1^2 / mu_hm_1^2 < sigma_hm_2^2 / mu_hm_2^2.  Now that zeta has to
    ! be greater than 0, sigma_hm_1^2 / mu_hm_1^2 > sigma_hm_2^2 / mu_hm_2^2.
    ! The ratio of the greater variance-over-mean-squared to the smaller
    ! variance-over-mean-squared remains the same when using the equation for
    ! zeta listed above.
    !
    ! When mu_thl_1 > mu_thl_2:
    !
    ! The relationship of mu_hm_1 to mu_hm_2 is mu_hm_1 <= mu_hm_2, so
    ! zeta <= 0.  When zeta_in is already less than 0 (meaning
    ! sigma_hm_1^2 / mu_hm_1^2 is less than sigma_hm_2^2 / mu_hm_2^2), zeta is
    ! simply set to zeta_in.  In other words, the component with the greater
    ! mean value of the hydrometeor also has the greater value of component
    ! variance to the square of the component mean.  However, when zeta_in is
    ! greater than 0, zeta must be adjusted to be less than 0.  The following
    ! equation is used to set the value of zeta:
    !
    !        | zeta_in; when zeta_in <= 0
    ! zeta = |
    !        | ( 1 / ( 1 + zeta_in ) ) - 1; when zeta_in > 0.
    !
    ! Previously, when zeta_in > 0, mu_hm_1 > mu_hm_2, and
    ! sigma_hm_1^2 / mu_hm_1^2 > sigma_hm_2^2 / mu_hm_2^2.  Now that zeta has to
    ! be less than 0, sigma_hm_1^2 / mu_hm_1^2 < sigma_hm_2^2 / mu_hm_2^2.
    ! The ratio of the greater variance-over-mean-squared to the smaller
    ! variance-over-mean-squared remains the same when using the equation for
    ! zeta listed above.
    !
    !
    ! Brian Griffin; February 2015.

    ! References:
    !----------------------------------------------------------------------- 

    use constants_clubb, only: &
        four, & ! Constant(s)
        two,  &
        one,  &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none
    
    ! Input Variables
    integer, intent(in) :: &
      nz,      & ! Number of model vertical grid levels
      ngrdcol    ! Number of grid columns

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      hmm,           & ! Hydrometeor mean (overall), <hm>           [hm units]
      hmp2,          & ! Hydrometeor variance (overall), <hm'^2>    [hm units^2]
      mixt_frac,     & ! Mixture fraction                           [-]
      precip_frac,   & ! Precipitation fraction (overall)           [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component) [-]
      precip_frac_2, & ! Precipitation fraction (2nd PDF component) [-]
      mu_thl_1,      & ! Mean of th_l (1st PDF component)           [K]
      mu_thl_2         ! Mean of th_l (2nd PDF component)           [K]

    real( kind = core_rknd ), intent(in) :: &
      hmp2_ip_on_hmm2_ip, & ! Ratio <hm|_ip'^2> / <hm|_ip>^2                 [-]
      hm_tol                ! Tolerance value of hydrometeor          [hm units]
      
    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      precip_frac_tol       ! Min. precip. frac. when hydromet. are present  [-]

    real( kind = core_rknd ), intent(in) :: &
      omicron,           & ! Relative width parameter, omicron = R / Rmax    [-]
      zeta_vrnce_rat_in    ! Width parameter for sigma_hm_1^2 / mu_hm_1^2    [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      mu_hm_1,    & ! Mean of hm (1st PDF component) in-precip (ip)      [hm un]
      mu_hm_2,    & ! Mean of hm (2nd PDF component) ip                  [hm un]
      sigma_hm_1, & ! Standard deviation of hm (1st PDF component) ip    [hm un]
      sigma_hm_2, & ! Standard deviation of hm (2nd PDF component) ip    [hm un]
      hm_1,       & ! Mean of hm (1st PDF component)                     [hm un]
      hm_2          ! Mean of hm (2nd PDF component)                     [hm un]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      sigma_hm_1_sqd_on_mu_hm_1_sqd, & ! Ratio sigma_hm_1**2 / mu_hm_1**2    [-]
      sigma_hm_2_sqd_on_mu_hm_2_sqd    ! Ratio sigma_hm_2**2 / mu_hm_2**2    [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      Rmax,       & ! Maximum possible value of ratio R                [-]
      coef_A,     & ! Coefficient A in A*mu_hm_1^2 + B*mu_hm_1 + C = 0 [-]
      coef_B,     & ! Coefficient B in A*mu_hm_1^2 + B*mu_hm_1 + C = 0 [hm un]
      coef_C,     & ! Coefficient C in A*mu_hm_1^2 + B*mu_hm_1 + C = 0 [hm un^2]
      Bsqd_m_4AC    ! Value B^2 - 4*A*C in quadratic eqn. for mu_hm_1  [hm un^2]
      
    real( kind = core_rknd ) :: &
      zeta_vrnce_rat    ! Width parameter for sigma_hm_1^2 / mu_hm_1^2       [-]

    real( kind = core_rknd ) :: &
      mu_hm_1_min, & ! Minimum value of mu_hm_1 (precip. in both comps.) [hm un]
      mu_hm_2_min    ! Minimum value of mu_hm_2 (precip. in both comps.) [hm un]

    real( kind = core_rknd ), parameter :: &
      mu_hm_min_coef = 0.01_core_rknd  ! Coef. for mu_hm_1_min and mu_hm_2_min

    integer :: i, k ! Loop iterators
      
    do k = 1, nz
      do i = 1, ngrdcol

        if ( hmm(i,k) >= hm_tol &
             .and. precip_frac_1(i,k) >= precip_frac_tol(i) &
             .and. precip_frac_2(i,k) >= precip_frac_tol(i) ) then

           ! Adjust the value of zeta based on the relationship of mu_thl_1 to
           ! mu_thl_2.
           if ( mu_thl_1(i,k) <= mu_thl_2(i,k) ) then
              if ( zeta_vrnce_rat_in >= zero ) then
                 zeta_vrnce_rat = zeta_vrnce_rat_in
              else ! zeta_vrnce_rat_in < 0
                 zeta_vrnce_rat = ( one / ( one + zeta_vrnce_rat_in ) ) - one
              end if ! zeta_vrnce_rat_in >= 0
           else ! mu_thl_1 > mu_thl_2
              if ( zeta_vrnce_rat_in <= zero ) then
                 zeta_vrnce_rat = zeta_vrnce_rat_in
              else ! zeta_vrnce_rat_in > 0
                 zeta_vrnce_rat = ( one / ( one + zeta_vrnce_rat_in ) ) - one
              end if ! zeta_vrnce_rat_in <= 0
           end if ! mu_thl_1 <= mu_thl_2

           ! Calculate the value of Rmax.
           ! Rmax = ( f_p / ( a * f_p_1 * ( 1 + zeta ) + ( 1 - a ) * f_p_2 ) )
           !        * ( <hm|_ip’^2> / <hm|_ip>^2 ).
           ! The parameter zeta is written in the code as zeta_vrnce_rat.
           Rmax = ( precip_frac(i,k) &
                    / ( mixt_frac(i,k) * precip_frac_1(i,k) * ( one + zeta_vrnce_rat ) &
                        + ( one - mixt_frac(i,k) ) * precip_frac_2(i,k) ) ) &
                  * hmp2_ip_on_hmm2_ip

           ! Calculate the value of coefficient A.
           ! A = a * f_p_1 * ( 1 + omicron * Rmax * ( 1 + zeta ) )
           !     + a^2 * f_p_1^2 * ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 ).
           coef_A = mixt_frac(i,k) * precip_frac_1(i,k) &
                    * ( one + omicron * Rmax * ( one + zeta_vrnce_rat ) ) &
                    + mixt_frac(i,k)**2 * precip_frac_1(i,k)**2 &
                      * ( one + omicron * Rmax ) &
                      / ( ( one - mixt_frac(i,k) ) * precip_frac_2(i,k) )

           ! Calculate the value of coefficient B.
           ! B = - 2 * <hm> * a * f_p_1 * ( 1 + omicron * Rmax )
           !     / ( ( 1 - a ) * f_p_2 ).
           coef_B = -two * hmm(i,k) * mixt_frac(i,k) * precip_frac_1(i,k) &
                     * ( one + omicron * Rmax ) &
                     / ( ( one - mixt_frac(i,k) ) * precip_frac_2(i,k) )

           ! Calculate the value of coefficient C.
           ! C = - ( <hm’^2>
           !         + ( 1 - ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 ) )
           !           * <hm>^2 ).
           coef_C = - ( hmp2(i,k) + ( one &
                                 - ( one + omicron * Rmax ) &
                                   / ( ( one - mixt_frac(i,k) ) * precip_frac_2(i,k) ) &
                               ) * hmm(i,k)**2 )

           ! Calculate value of B^2 - 4*A*C.
           Bsqd_m_4AC = coef_B**2 - four * coef_A * coef_C

           ! Mathematically, the value of B^2 - 4*A*C cannot be less than 0.
           ! Numerically, this can happen when numerical round off error causes an
           ! epsilon-sized negative value.  When this happens, reset the value of
           ! B^2 - 4*A*C to 0.
           if ( Bsqd_m_4AC < zero ) then
              Bsqd_m_4AC = zero
           end if

           ! Calculate the mean (in-precip.) of the hydrometeor in the 1st PDF
           ! component.
           if ( mu_thl_1(i,k) <= mu_thl_2(i,k) ) then
              mu_hm_1(i,k) = ( -coef_B + sqrt( Bsqd_m_4AC ) ) / ( two * coef_A )
           else ! mu_thl_1 > mu_thl_2
              mu_hm_1(i,k) = ( -coef_B - sqrt( Bsqd_m_4AC ) ) / ( two * coef_A )
           end if ! mu_thl_1 <= mu_thl_2

           ! Calculate the mean (in-precip.) of the hydrometeor in the 2nd PDF
           ! component.
           mu_hm_2(i,k) = ( hmm(i,k) - mixt_frac(i,k) * precip_frac_1(i,k) * mu_hm_1(i,k) ) &
                     / ( ( one - mixt_frac(i,k) ) * precip_frac_2(i,k) )

           ! Calculate the value of the ratio R (which is sigma_hm_2^2 / mu_hm_2^2),
           ! where R = omicron * Rmax.  The name of the variable used for R is
           ! sigma_hm_2_sqd_on_mu_hm_2_sqd.
           sigma_hm_2_sqd_on_mu_hm_2_sqd(i,k) = omicron * Rmax

           ! Calculate minimum allowable values for mu_hm_1 and mu_hm_2.
           if ( hmm(i,k) / precip_frac(i,k) > hm_tol / precip_frac_1(i,k) ) then
              mu_hm_1_min &
              = min( hm_tol / precip_frac_1(i,k) &
                     + mu_hm_min_coef * ( hmm(i,k) / precip_frac(i,k) &
                                          - hm_tol / precip_frac_1(i,k) ), &
                     ( hmm(i,k) - ( one - mixt_frac(i,k) ) * hm_tol ) &
                     / ( mixt_frac(i,k) * precip_frac_1(i,k) ) )
           else ! hmm / precip_frac <= hm_tol / precip_frac_1
              mu_hm_1_min = hm_tol / precip_frac_1(i,k)
           end if
           if ( hmm(i,k) / precip_frac(i,k) > hm_tol / precip_frac_2(i,k) ) then
              mu_hm_2_min &
              = min( hm_tol / precip_frac_2(i,k) &
                     + mu_hm_min_coef * ( hmm(i,k) / precip_frac(i,k) &
                                          - hm_tol / precip_frac_2(i,k) ), &
                     ( hmm(i,k) - mixt_frac(i,k) * hm_tol ) &
                     / ( ( one - mixt_frac(i,k) ) * precip_frac_2(i,k) ) )
           else ! hmm / precip_frac <= hm_tol / precip_frac_2
              mu_hm_2_min = hm_tol / precip_frac_2(i,k)
           end if

           ! Handle the "emergency" situation when the specified value of omicron is
           ! too small for the value of <hm|_ip'^2> / <hm|_ip>^2, resulting in a
           ! component mean that is too small (below tolerance value) or negative.
           if ( mu_hm_1(i,k) < mu_hm_1_min ) then

              ! Set the value of mu_hm_1 to the threshold positive value.
              mu_hm_1(i,k) = mu_hm_1_min

              ! Recalculate the mean (in-precip.) of the hydrometeor in the 2nd PDF
              ! component.
              mu_hm_2(i,k) = ( hmm(i,k) - mixt_frac(i,k) * precip_frac_1(i,k) * mu_hm_1(i,k) ) &
                        / ( ( one - mixt_frac(i,k) ) * precip_frac_2(i,k) )

              ! Recalculate the value of R ( sigma_hm_2^2 / mu_hm_2^2 ) in this
              ! scenario.
              ! R = ( <hm'^2> + <hm>^2 - a * f_p_1 * mu_hm_1^2
              !       - ( 1 - a ) * f_p_2 * mu_hm_2^2 )
              !     / ( a * f_p_1 * ( 1 + zeta ) * mu_hm_1^2
              !         + ( 1 - a ) * f_p_2 * mu_hm_2^2 ).
              sigma_hm_2_sqd_on_mu_hm_2_sqd(i,k) &
              = ( hmp2(i,k) + hmm(i,k)**2 - mixt_frac(i,k) * precip_frac_1(i,k) * mu_hm_1(i,k)**2 &
                  - ( one - mixt_frac(i,k) ) * precip_frac_2(i,k) * mu_hm_2(i,k)**2 ) &
                / ( mixt_frac(i,k) * precip_frac_1(i,k) &
                    * ( one + zeta_vrnce_rat ) * mu_hm_1(i,k)**2 &
                    + ( one - mixt_frac(i,k) ) * precip_frac_2(i,k) * mu_hm_2(i,k)**2 )

              ! Mathematically, this ratio can never be less than 0.  In case
              ! numerical round off error produces a negative value in extreme
              ! cases, reset the value of R to 0.
              if ( sigma_hm_2_sqd_on_mu_hm_2_sqd(i,k) < zero ) then
                 sigma_hm_2_sqd_on_mu_hm_2_sqd(i,k) = zero
              end if

           elseif ( mu_hm_2(i,k) < mu_hm_2_min ) then

              ! Set the value of mu_hm_2 to the threshold positive value.
              mu_hm_2(i,k) = mu_hm_2_min

              ! Recalculate the mean (in-precip.) of the hydrometeor in the 1st PDF
              ! component.
              mu_hm_1(i,k) = ( hmm(i,k) - ( one - mixt_frac(i,k) ) &
                                          * precip_frac_2(i,k) * mu_hm_2(i,k) ) &
                             / ( mixt_frac(i,k) * precip_frac_1(i,k) )

              ! Recalculate the value of R ( sigma_hm_2^2 / mu_hm_2^2 ) in this
              ! scenario.
              ! R = ( <hm'^2> + <hm>^2 - a * f_p_1 * mu_hm_1^2
              !       - ( 1 - a ) * f_p_2 * mu_hm_2^2 )
              !     / ( a * f_p_1 * ( 1 + zeta ) * mu_hm_1^2
              !         + ( 1 - a ) * f_p_2 * mu_hm_2^2 ).
              sigma_hm_2_sqd_on_mu_hm_2_sqd(i,k) &
              = ( hmp2(i,k) + hmm(i,k)**2 - mixt_frac(i,k) * precip_frac_1(i,k) * mu_hm_1(i,k)**2 &
                  - ( one - mixt_frac(i,k) ) * precip_frac_2(i,k) * mu_hm_2(i,k)**2 ) &
                / ( mixt_frac(i,k) * precip_frac_1(i,k) &
                    * ( one + zeta_vrnce_rat ) * mu_hm_1(i,k)**2 &
                    + ( one - mixt_frac(i,k) ) * precip_frac_2(i,k) * mu_hm_2(i,k)**2 )

              ! Mathematically, this ratio can never be less than 0.  In case
              ! numerical round off error produces a negative value in extreme
              ! cases, reset the value of R to 0.
              if ( sigma_hm_2_sqd_on_mu_hm_2_sqd(i,k) < zero ) then
                 sigma_hm_2_sqd_on_mu_hm_2_sqd(i,k) = zero
              end if

           end if
     
           ! Calculate the standard deviation (in-precip.) of the hydrometeor in the
           ! 1st PDF component.
           sigma_hm_1(i,k) = sqrt( sigma_hm_2_sqd_on_mu_hm_2_sqd(i,k) &
                              * ( one + zeta_vrnce_rat ) ) &
                        * mu_hm_1(i,k)

           ! Calculate the standard deviation (in-precip.) of the hydrometeor in the
           ! 2nd PDF component.
           sigma_hm_2(i,k) = sqrt( sigma_hm_2_sqd_on_mu_hm_2_sqd(i,k) ) * mu_hm_2(i,k)

           ! Calculate the mean of the hydrometeor in the 1st PDF component.
           hm_1(i,k) = max( mu_hm_1(i,k) * precip_frac_1(i,k), hm_tol )

           ! Calculate the mean of the hydrometeor in the 1st PDF component.
           hm_2(i,k) = max( mu_hm_2(i,k) * precip_frac_2(i,k), hm_tol )

           ! Calculate the ratio of sigma_hm_1^2 / mu_hm_1^2.
           sigma_hm_1_sqd_on_mu_hm_1_sqd(i,k) = sigma_hm_1(i,k)**2 / mu_hm_1(i,k)**2

           ! The value of R, sigma_hm_2_sqd_on_mu_hm_2_sqd, has already been
           ! calculated.


        elseif ( hmm(i,k) >= hm_tol .and. precip_frac_1(i,k) >= precip_frac_tol(i) ) then

           ! Precipitation is found in the 1st PDF component, but not in the 2nd
           ! PDF component (precip_frac_2 = 0).
           mu_hm_1(i,k) = hmm(i,k) / ( mixt_frac(i,k) * precip_frac_1(i,k) )
           mu_hm_2(i,k) = zero

           sigma_hm_1(i,k) = sqrt( max( ( hmp2(i,k) + hmm(i,k)**2 &
                                     - mixt_frac(i,k) * precip_frac_1(i,k) * mu_hm_1(i,k)**2 ) &
                                   / ( mixt_frac(i,k) * precip_frac_1(i,k) ), &
                                   zero ) )
           sigma_hm_2(i,k) = zero

           hm_1(i,k) = mu_hm_1(i,k) * precip_frac_1(i,k)
           hm_2(i,k) = zero

           sigma_hm_1_sqd_on_mu_hm_1_sqd(i,k) = sigma_hm_1(i,k)**2 / mu_hm_1(i,k)**2
           ! The ratio sigma_hm_2^2 / mu_hm_2^2 is undefined.
           sigma_hm_2_sqd_on_mu_hm_2_sqd(i,k) = zero


        elseif ( hmm(i,k) >= hm_tol .and. precip_frac_2(i,k) >= precip_frac_tol(i) ) then

           ! Precipitation is found in the 2nd PDF component, but not in the 1st
           ! PDF component (precip_frac_1 = 0).
           mu_hm_1(i,k) = zero
           mu_hm_2(i,k) = hmm(i,k) / ( ( one - mixt_frac(i,k) ) * precip_frac_2(i,k) )

           sigma_hm_1(i,k) = zero
           sigma_hm_2(i,k) &
           = sqrt( max( ( hmp2(i,k) + hmm(i,k)**2 &
                          - ( one - mixt_frac(i,k) ) * precip_frac_2(i,k) * mu_hm_2(i,k)**2 ) &
                        / ( ( one - mixt_frac(i,k) ) * precip_frac_2(i,k) ), &
                        zero ) )

           hm_1(i,k) = zero
           hm_2(i,k) = mu_hm_2(i,k) * precip_frac_2(i,k)

           ! The ratio sigma_hm_1^2 / mu_hm_1^2 is undefined.
           sigma_hm_1_sqd_on_mu_hm_1_sqd(i,k) = zero
           sigma_hm_2_sqd_on_mu_hm_2_sqd(i,k) = sigma_hm_2(i,k)**2 / mu_hm_2(i,k)**2


        else ! hm < hm_tol or ( precip_frac_1 = 0 and precip_frac_2 = 0 ).

           ! Precipitation is not found in either PDF component.
           mu_hm_1(i,k) = zero
           mu_hm_2(i,k) = zero

           sigma_hm_1(i,k) = zero
           sigma_hm_2(i,k) = zero

           hm_1(i,k) = zero
           hm_2(i,k) = zero

           ! The ratio sigma_hm_1^2 / mu_hm_1^2 is undefined.
           sigma_hm_1_sqd_on_mu_hm_1_sqd(i,k) = zero
           ! The ratio sigma_hm_2^2 / mu_hm_2^2 is undefined.
           sigma_hm_2_sqd_on_mu_hm_2_sqd(i,k) = zero


        end if ! hmm >= hm_tol and precip_frac_1 >= precip_frac_tol
              ! and precip_frac_2 >= precip_frac_tol
      end do
    end do

    return

  end subroutine calc_comp_mu_sigma_hm

  !=============================================================================
  subroutine component_corr_w_x( nz, ngrdcol, &
                                 rc_1, &
                                 rc_2, &
                                 corr_w_x_NN_cloud, corr_w_x_NN_below, &
                                 iiPDF_type, &
                                 corr_w_x_1, &
                                 corr_w_x_2 )

    ! Description:
    ! Calculates the correlation of w and x within the ith PDF component.
    ! Here, x is a variable with a normally distributed individual marginal PDF,
    ! such as chi or eta.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        zero,   &
        rc_tol

    use model_flags, only: &
        iiPDF_ADG1,       & ! Variable(s)
        iiPDF_ADG2,       &
        iiPDF_new_hybrid

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,     & ! Number of model vertical grid levels
      ngrdcol   ! Number of grid columns

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      rc_1,           & ! Mean cloud water mixing ratio (1st PDF comp.)  [kg/kg]
      rc_2              ! Mean cloud water mixing ratio (2nd PDF comp.)  [kg/kg]

    real( kind = core_rknd ), intent(in) :: &
      corr_w_x_NN_cloud, & ! Corr. of w and x (ith PDF comp.); cloudy levs [-]
      corr_w_x_NN_below    ! Corr. of w and x (ith PDF comp.); clear levs  [-]

    integer, intent(in) :: &
      iiPDF_type    ! Selected option for the two-component normal (double
                    ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                    ! w, chi, and eta) portion of CLUBB's multivariate,
                    ! two-component PDF.

    ! Output Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      corr_w_x_1, & ! Correlation of w and x (1st PDF component)  [-]
      corr_w_x_2    ! Correlation of w and x (2nd PDF component)  [-]
      
    ! Local Variables
    integer :: i, k 
    
    ! Correlation of w and x in the ith PDF component.

    ! The PDF variables chi and eta result from a transformation of the PDF
    ! involving r_t and theta_l.  The correlation of w and x (whether x is chi
    ! or eta) depends on the correlation of w and r_t, the correlation of w and
    ! theta_l, as well as the variances of r_t and theta_l, and other factors.
    ! The correlation of w and x is subject to change at every vertical level
    ! and model time step, and is calculated as part of the CLUBB PDF
    ! parameters.

    ! The ADG1 PDF fixes the correlation of w and rt and the correlation of
    ! w and theta_l to be 0, which means the correlation of w and chi and the
    ! correlation of w and eta must also be 0.
    if ( ( iiPDF_type == iiPDF_ADG1 .or. iiPDF_type == iiPDF_ADG2 &
            .or. iiPDF_type == iiPDF_new_hybrid ) &
          .and. l_follow_ADG1_PDF_standards ) then
          
      corr_w_x_1(:,:) = zero
      corr_w_x_2(:,:) = zero

    else ! use prescribed paramter values
        
      do k = 1, nz
        do i = 1, ngrdcol 
          
          if ( rc_1(i,k) > rc_tol ) then
            corr_w_x_1(i,k) = corr_w_x_NN_cloud
          else
            corr_w_x_1(i,k) = corr_w_x_NN_below
          end if
          
          if ( rc_2(i,k) > rc_tol ) then
            corr_w_x_2(i,k) = corr_w_x_NN_cloud
          else
            corr_w_x_2(i,k) = corr_w_x_NN_below
          end if
          
        end do
      end do
       
    end if ! iiPDF_type

    return

  end subroutine component_corr_w_x

  !=============================================================================
  subroutine component_corr_chi_eta( nz, ngrdcol, &
                                     rc_1, &
                                     rc_2, &
                                     corr_chi_eta_NN_cloud, &
                                     corr_chi_eta_NN_below, &
                                     l_limit_corr_chi_eta, &
                                     corr_chi_eta_1,  &
                                     corr_chi_eta_2 )

    ! Description:
    ! Calculates the correlation of chi (old s) and eta (old t) within the
    ! ith PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        rc_tol, &
        max_mag_correlation

    use clubb_precision, only: &
        core_rknd  ! Constant

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,     & ! Number of model vertical grid levels
      ngrdcol   ! Number of grid columns

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      rc_1,               & ! Mean cloud water mix. rat. (1st PDF comp.) [kg/kg]&
      rc_2                  ! Mean cloud water mix. rat. (2nd PDF comp.) [kg/kg]

    real( kind = core_rknd ), intent(in) :: &
      corr_chi_eta_NN_cloud, & ! Corr. of chi & eta (ith PDF comp.); cloudy  [-]
      corr_chi_eta_NN_below    ! Corr. of chi & eta (ith PDF comp.); clear   [-]

    logical, intent(in) :: &
      l_limit_corr_chi_eta    ! We must limit the correlation of chi and eta if
                              ! we are to take the Cholesky decomposition of the
                              ! resulting correlation matrix. This is because a
                              ! perfect correlation of chi and eta was found to
                              ! be unrealizable.

    ! Output Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      corr_chi_eta_1, & ! Correlation of chi and eta (1st PDF component)     [-]
      corr_chi_eta_2    ! Correlation of chi and eta (2nd PDF component)     [-]

    ! Local Variables
    integer :: i, k
    
    ! Correlation of chi (old s) and eta (old t) in the ith PDF component.

    ! The PDF variables chi and eta result from a transformation of the PDF
    ! involving r_t and theta_l.  The correlation of chi and eta depends on the
    ! correlation of r_t and theta_l, as well as the variances of r_t and
    ! theta_l, and other factors.  The correlation of chi and eta is subject to
    ! change at every vertical level and model time step, and is calculated as
    ! part of the CLUBB PDF parameters.

    ! WARNING:  this code is inconsistent with the rest of CLUBB's PDF.  This
    !           code is necessary because SILHS is lazy and wussy, and only
    !           wants to declare correlation arrays at the start of the model
    !           run, rather than updating them throughout the model run.


    do k = 1, nz
      do i = 1, ngrdcol
           
        if ( rc_1(i,k) > rc_tol ) then
          corr_chi_eta_1(i,k) = corr_chi_eta_NN_cloud
        else
          corr_chi_eta_1(i,k) = corr_chi_eta_NN_below
        end if
      end do
    end do
      
    do k = 1, nz
      do i = 1, ngrdcol
        if ( rc_2(i,k) > rc_tol ) then
          corr_chi_eta_2(i,k) = corr_chi_eta_NN_cloud
        else
          corr_chi_eta_2(i,k) = corr_chi_eta_NN_below
        end if
      end do
    end do
        
   

    ! We cannot have a perfect correlation of chi (old s) and eta (old t) if we
    ! plan to decompose this matrix and we don't want the Cholesky_factor code
    ! to throw a fit.
    if ( l_limit_corr_chi_eta ) then

      do k = 1, nz
        do i = 1, ngrdcol
          corr_chi_eta_1(i,k) = max( min( corr_chi_eta_1(i,k), max_mag_correlation ), &
                                     -max_mag_correlation )
        end do
      end do
      
      do k = 1, nz
        do i = 1, ngrdcol
          corr_chi_eta_2(i,k) = max( min( corr_chi_eta_2(i,k), max_mag_correlation ), &
                                     -max_mag_correlation )
        end do
      end do

    end if
    
    return

  end subroutine component_corr_chi_eta

  !=============================================================================
  subroutine component_corr_w_hm_n_ip( nz, ngrdcol, &
                                       corr_w_hm_1_n_in, rc_1, &
                                       corr_w_hm_2_n_in, rc_2, &
                                       corr_w_hm_n_NL_cloud, &
                                       corr_w_hm_n_NL_below, &
                                       l_calc_w_corr, &
                                       corr_w_hm_1_n, &
                                       corr_w_hm_2_n )
                                       
    ! Description:
    ! Calculates the in-precip correlation of w and the natural logarithm of a
    ! hydrometeor species within the ith PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        rc_tol

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,     & ! Number of model vertical grid levels
      ngrdcol   ! Number of grid columns

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      corr_w_hm_1_n_in, & ! Correlation of w and ln hm (1st PDF comp.) ip    [-]
      corr_w_hm_2_n_in, & ! Correlation of w and ln hm (2nd PDF comp.) ip    [-]
      rc_1,             & ! Mean cloud water mix. ratio (1st PDF comp.)  [kg/kg]
      rc_2                ! Mean cloud water mix. ratio (2nd PDF comp.)  [kg/kg]

    real( kind = core_rknd ), intent(in) :: &
      corr_w_hm_n_NL_cloud, & ! Corr. of w & ln hm (ith PDF comp.) ip; cloud [-]
      corr_w_hm_n_NL_below    ! Corr. of w & ln hm (ith PDF comp.) ip; clear [-]

    logical, intent(in) :: &
      l_calc_w_corr ! Calculate the correlations between w and the hydrometeors

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      corr_w_hm_1_n, & ! Correlation of w and ln hm (1st PDF component) ip   [-]
      corr_w_hm_2_n    ! Correlation of w and ln hm (2nd PDF component) ip   [-]
      
    ! Local Variables
    integer :: i, k


    ! Correlation (in-precip) of w and the natural logarithm of the hydrometeor
    ! in the ith PDF component.
    if ( l_calc_w_corr ) then
      
      corr_w_hm_1_n(:,:) = corr_w_hm_1_n_in(:,:)
      corr_w_hm_2_n(:,:) = corr_w_hm_2_n_in(:,:)


    else

        do k = 1, nz
          do i = 1, ngrdcol
            
            if ( rc_1(i,k) > rc_tol ) then
               corr_w_hm_1_n(i,k) = corr_w_hm_n_NL_cloud
            else
               corr_w_hm_1_n(i,k) = corr_w_hm_n_NL_below
            end if
            
            if ( rc_2(i,k) > rc_tol ) then
               corr_w_hm_2_n(i,k) = corr_w_hm_n_NL_cloud
            else
               corr_w_hm_2_n(i,k) = corr_w_hm_n_NL_below
            end if
          end do
        end do
       
    end if ! l_calc_w_corr

    return

  end subroutine component_corr_w_hm_n_ip

  !=============================================================================
  subroutine component_corr_x_hm_n_ip( nz, ngrdcol, &
                                       rc_1, &
                                       rc_2, &
                                       corr_x_hm_n_NL_cloud, &
                                       corr_x_hm_n_NL_below, &
                                       corr_x_hm_1_n, &
                                       corr_x_hm_2_n ) 

    ! Description:
    ! Calculates the in-precip correlation of x and a hydrometeor species
    ! within the ith PDF component.  Here, x is a variable with a normally
    ! distributed individual marginal PDF, such as chi or eta.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        rc_tol

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,     & ! Number of model vertical grid levels
      ngrdcol   ! Number of grid columns

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      rc_1,         & ! Mean cloud water mixing ratio (1st PDF comp.)   [kg/kg]
      rc_2            ! Mean cloud water mixing ratio (2nd PDF comp.)   [kg/kg]


    real( kind = core_rknd ), intent(in) :: &
      corr_x_hm_n_NL_cloud, & ! Corr. of x and ln hm (ith PDF comp.) ip     [-]
      corr_x_hm_n_NL_below    ! Corr. of x and ln hm (ith PDF comp.) ip     [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      corr_x_hm_1_n, & ! Correlation of x and ln hm (1st PDF component) ip  [-]
      corr_x_hm_2_n    ! Correlation of x and ln hm (2nd PDF component) ip  [-]

    ! Local Variables
    integer :: i, k

    do k = 1, nz
      do i = 1, ngrdcol
      
        if ( rc_1(i,k) > rc_tol ) then
          corr_x_hm_1_n(i,k) = corr_x_hm_n_NL_cloud
        else
          corr_x_hm_1_n(i,k) = corr_x_hm_n_NL_below
        end if
        
        if ( rc_2(i,k) > rc_tol ) then
          corr_x_hm_2_n(i,k) = corr_x_hm_n_NL_cloud
        else
          corr_x_hm_2_n(i,k) = corr_x_hm_n_NL_below
        end if
        
      end do
    end do

    return

  end subroutine component_corr_x_hm_n_ip

  !=============================================================================
  subroutine component_corr_hmx_hmy_n_ip( nz, ngrdcol, &  
                                          rc_1, &
                                          rc_2, &
                                          corr_hmx_hmy_n_LL_cloud, &
                                          corr_hmx_hmy_n_LL_below, &
                                          corr_hmx_hmy_1_n, &
                                          corr_hmx_hmy_2_n ) 

    ! Description:
    ! Calculates the in-precip correlation of the natural logarithms of
    ! hydrometeor x and hydrometeor y within the ith PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        rc_tol

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,     & ! Number of model vertical grid levels
      ngrdcol   ! Number of grid columns

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      rc_1,         & ! Mean cloud water mixing ratio (ith PDF comp.) [kg/kg]
      rc_2            ! Mean cloud water mixing ratio (ith PDF comp.) [kg/kg]

    real( kind = core_rknd ), intent(in) :: &
      corr_hmx_hmy_n_LL_cloud, & ! Corr.: ln hmx & ln hmy (ith PDF comp.) ip [-]
      corr_hmx_hmy_n_LL_below    ! Corr.: ln hmx & ln hmy (ith PDF comp.) ip [-]

    ! Output Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      corr_hmx_hmy_1_n, & ! Corr. of ln hmx & ln hmy (ith PDF comp.) ip      [-]
      corr_hmx_hmy_2_n    ! Corr. of ln hmx & ln hmy (ith PDF comp.) ip      [-]

    ! Local Variable
    integer :: i, k
    

    ! Correlation (in-precip) of the natural logarithms of hydrometeor x and
    ! hydrometeor y in the ith PDF component.
    do k = 1, nz
      do i = 1, ngrdcol
        
        if ( rc_1(i,k) > rc_tol ) then
          corr_hmx_hmy_1_n(i,k) = corr_hmx_hmy_n_LL_cloud
        else
          corr_hmx_hmy_1_n(i,k) = corr_hmx_hmy_n_LL_below
        end if
        
        if ( rc_2(i,k) > rc_tol ) then
          corr_hmx_hmy_2_n(i,k) = corr_hmx_hmy_n_LL_cloud
        else
          corr_hmx_hmy_2_n(i,k) = corr_hmx_hmy_n_LL_below
        end if
        
      end do
    end do

    return

  end subroutine component_corr_hmx_hmy_n_ip

  !=============================================================================
  subroutine component_corr_eta_hm_n_ip( nz, ngrdcol, &
                                         corr_chi_eta_1, &
                                         corr_chi_hm_n_1, &
                                         corr_chi_eta_2, &
                                         corr_chi_hm_n_2, &
                                         corr_eta_hm_n_1, &
                                         corr_eta_hm_n_2) 
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
    integer, intent(in) :: &
      nz,     & ! Number of model vertical grid levels
      ngrdcol   ! Number of grid columns


    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      corr_chi_eta_1,  & ! Component correlation of chi and eta              [-]
      corr_chi_eta_2,  & ! Component correlation of chi and eta              [-]
      corr_chi_hm_n_1, & ! Component correlation of chi and ln hm            [-]
      corr_chi_hm_n_2    ! Component correlation of chi and ln hm            [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      corr_eta_hm_n_1, & ! Component correlation of eta and ln hm            [-]
      corr_eta_hm_n_2    ! Component correlation of eta and ln hm            [-]


    corr_eta_hm_n_1 = corr_chi_eta_1 * corr_chi_hm_n_1
    corr_eta_hm_n_2 = corr_chi_eta_2 * corr_chi_hm_n_2


    return

  end subroutine component_corr_eta_hm_n_ip

  !=============================================================================
  subroutine norm_transform_mean_stdev( nz, ngrdcol, &
                                        hm_1, hm_2, &
                                        Ncnm, pdf_dim, &
                                        mu_x_1, mu_x_2, &
                                        sigma_x_1, sigma_x_2, &
                                        sigma2_on_mu2_ip_1, &
                                        sigma2_on_mu2_ip_2, &
                                        l_const_Nc_in_cloud, &
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

    use array_index, only: &
        iiPDF_w,   & 
        iiPDF_chi, &
        iiPDF_eta, &
        iiPDF_Ncn, &  ! Variable(s)
        hydromet_tol  ! Variable(s)

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,      & ! Number of model vertical grid levels
      pdf_dim, & ! Number of variables in CLUBB's PDF
      ngrdcol    ! Number of grid columns

    real( kind = core_rknd ), dimension(ngrdcol,nz,hydromet_dim), intent(in) :: &
      hm_1, & ! Mean of a precip. hydrometeor (1st PDF component)   [units vary]
      hm_2    ! Mean of a precip. hydrometeor (2nd PDF component)   [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      Ncnm    ! Mean cloud nuclei concentration, < N_cn >               [num/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim), intent(in) :: &
      mu_x_1,    & ! Mean array of PDF vars. (1st PDF component)    [units vary]
      mu_x_2,    & ! Mean array of PDF vars. (2nd PDF component)    [units vary]
      sigma_x_1, & ! Standard deviation array of PDF vars (comp. 1) [units vary]
      sigma_x_2    ! Standard deviation array of PDF vars (comp. 2) [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim), intent(in) :: &
      sigma2_on_mu2_ip_1, & ! Prescribed ratio array: sigma_hm_1^2/mu_hm_1^2 [-]
      sigma2_on_mu2_ip_2    ! Prescribed ratio array: sigma_hm_2^2/mu_hm_2^2 [-]

    logical, intent(in) :: &
      l_const_Nc_in_cloud ! Use a constant cloud droplet conc. within cloud (K&K)

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim), intent(out) :: &
      mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    ! Local Variable
    integer :: ivar, hm_idx, k, i  ! Indices


    ! The means and standard deviations in each PDF component of w, chi (old s),
    ! and eta (old t) do not need to be transformed to normal space, since w,
    ! chi, and eta already follow assumed normal distributions in each PDF
    ! component.  The normal space means and standard deviations are the same as
    ! the actual means and standard deviations.    
    do k = 1, nz
      mu_x_1_n(:,k,iiPDF_chi) = mu_x_1(:,k,iiPDF_chi)
      mu_x_1_n(:,k,iiPDF_eta) = mu_x_1(:,k,iiPDF_eta)
      mu_x_1_n(:,k,iiPDF_w)   = mu_x_1(:,k,iiPDF_w)
    end do
      
    do k = 1, nz  
      mu_x_2_n(:,k,iiPDF_chi) = mu_x_2(:,k,iiPDF_chi)
      mu_x_2_n(:,k,iiPDF_eta) = mu_x_2(:,k,iiPDF_eta)
      mu_x_2_n(:,k,iiPDF_w)   = mu_x_2(:,k,iiPDF_w)
    end do
    
    do k = 1, nz
      sigma_x_1_n(:,k,iiPDF_chi) = sigma_x_1(:,k,iiPDF_chi)
      sigma_x_1_n(:,k,iiPDF_eta) = sigma_x_1(:,k,iiPDF_eta)
      sigma_x_1_n(:,k,iiPDF_w)   = sigma_x_1(:,k,iiPDF_w)
    end do
      
    do k = 1, nz
      sigma_x_2_n(:,k,iiPDF_chi) = sigma_x_2(:,k,iiPDF_chi)
      sigma_x_2_n(:,k,iiPDF_eta) = sigma_x_2(:,k,iiPDF_eta)
      sigma_x_2_n(:,k,iiPDF_w)   = sigma_x_2(:,k,iiPDF_w)
    end do
    
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
    do k = 1, nz
      do i = 1, ngrdcol
        
        if ( Ncnm(i,k) >= Ncn_tol ) then

           mu_x_1_n(i,k,iiPDF_Ncn) = mean_L2N( mu_x_1(i,k,iiPDF_Ncn), &
                                               sigma2_on_mu2_ip_1(i,k,iiPDF_Ncn) )

        else

           ! Mean simplified cloud nuclei concentration in PDF component 1 is less
           ! than the tolerance amount.  It is considered to have a value of 0.
           ! There are not any cloud nuclei or cloud at this grid level.  The value
           ! of mu_Ncn_1_n should be -inf.  It will be set to -huge for purposes of
           ! assigning it a value.
           mu_x_1_n(i,k,iiPDF_Ncn) = -huge( mu_x_1(i,k,iiPDF_Ncn) )

        end if
        
      end do
    end do

    ! Normal space mean of simplified cloud nuclei concentration, N_cn,
    ! in PDF component 2.
    do k = 1, nz
      do i = 1, ngrdcol
        
        if ( Ncnm(i,k) >= Ncn_tol ) then

           mu_x_2_n(i,k,iiPDF_Ncn) = mean_L2N( mu_x_2(i,k,iiPDF_Ncn), &
                                               sigma2_on_mu2_ip_1(i,k,iiPDF_Ncn) )

        else

           ! Mean simplified cloud nuclei concentration in PDF component 1 is less
           ! than the tolerance amount.  It is considered to have a value of 0.
           ! There are not any cloud nuclei or cloud at this grid level.  The value
           ! of mu_Ncn_1_n should be -inf.  It will be set to -huge for purposes of
           ! assigning it a value.
           mu_x_2_n(i,k,iiPDF_Ncn) = -huge( mu_x_2(i,k,iiPDF_Ncn) )

        end if
        
      end do
    end do

    ! Normal space standard deviation of simplified cloud nuclei concentration,
    ! N_cn, in PDF components 1 and 2.
    if ( l_const_Nc_in_cloud ) then
      ! Ncn does not vary in the grid box.
      sigma_x_1_n(:,:,iiPDF_Ncn) = zero
      sigma_x_2_n(:,:,iiPDF_Ncn) = zero
    else
       ! Ncn (perhaps) varies in the grid box.
       sigma_x_1_n(:,:,iiPDF_Ncn) = stdev_L2N( sigma2_on_mu2_ip_1(:,:,iiPDF_Ncn) )
       sigma_x_2_n(:,:,iiPDF_Ncn) = stdev_L2N( sigma2_on_mu2_ip_2(:,:,iiPDF_Ncn) )
    end if

    ! Normal space precipitating hydrometeor means and standard deviations.
    do ivar = iiPDF_Ncn+1, pdf_dim, 1

      hm_idx = pdf2hydromet_idx(ivar)

      ! Normal space mean of a precipitating hydrometeor, hm, in PDF
      ! component 1.
      do k = 1, nz
        do i = 1, ngrdcol
          
          if ( hm_1(i,k,hm_idx) >= hydromet_tol(hm_idx) ) then

            mu_x_1_n(i,k,ivar) = mean_L2N( mu_x_1(i,k,ivar), &
                                           sigma2_on_mu2_ip_1(i,k,ivar) )

          else

            ! The mean of a precipitating hydrometeor in PDF component 1 is less
            ! than its tolerance amount.  It is considered to have a value of 0.
            ! There is not any of this precipitating hydrometeor in the 1st PDF
            ! component at this grid level.  The in-precip mean of this
            ! precipitating hydrometeor (1st PDF component) is also 0.  The value
            ! of mu_hm_1_n should be -inf.  It will be set to -huge for purposes
            ! of assigning it a value.
            mu_x_1_n(i,k,ivar) = -huge( mu_x_1(i,k,ivar) )

          end if
          
        end do
      end do

      ! Normal space standard deviation of a precipitating hydrometeor, hm, in
      ! PDF component 1.
      sigma_x_1_n(:,:,ivar) = stdev_L2N( sigma2_on_mu2_ip_1(:,:,ivar) )

      ! Normal space mean of a precipitating hydrometeor, hm, in PDF
      ! component 2.
      do k = 1, nz
        do i = 1, ngrdcol
          
          if ( hm_2(i,k,hm_idx) >= hydromet_tol(hm_idx) ) then

            mu_x_2_n(i,k,ivar) = mean_L2N( mu_x_2(i,k,ivar), &
                                         sigma2_on_mu2_ip_2(i,k,ivar) )

          else

            ! The mean of a precipitating hydrometeor in PDF component 2 is less
            ! than its tolerance amount.  It is considered to have a value of 0.
            ! There is not any of this precipitating hydrometeor in the 2nd PDF
            ! component at this grid level.  The in-precip mean of this
            ! precipitating hydrometeor (2nd PDF component) is also 0.  The value
            ! of mu_hm_2_n should be -inf.  It will be set to -huge for purposes
            ! of assigning it a value.
            mu_x_2_n(i,k,ivar) = -huge( mu_x_2(i,k,ivar) )

          end if
            
        end do
      end do

      ! Normal space standard deviation of a precipitating hydrometeor, hm, in
      ! PDF component 2.
      sigma_x_2_n(:,:,ivar) = stdev_L2N( sigma2_on_mu2_ip_2(:,:,ivar) )

    end do ! ivar = iiPDF_Ncn+1, pdf_dim, 1


    return

  end subroutine norm_transform_mean_stdev
  
  !=============================================================================
  subroutine denorm_transform_corr( nz, ngrdcol, pdf_dim, &
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

    use constants_clubb, only : &
        one

    use pdf_utilities, only: &
        corr_NN2NL, & ! Procedure(s)
        corr_NN2LL

    use array_index, only: &
        iiPDF_chi, & ! Variable(s)
        iiPDF_eta, &
        iiPDF_w,   &
        iiPDF_Ncn

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,      & ! Number of vertical levels
      pdf_dim, & ! Number of variables
      ngrdcol    ! Number of grid columns

    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim), intent(in) :: &
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim), intent(in) :: &
      sigma2_on_mu2_ip_1, & ! Ratio array sigma_hm_1^2/mu_hm_1^2             [-]
      sigma2_on_mu2_ip_2    ! Ratio array sigma_hm_2^2/mu_hm_2^2             [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim,pdf_dim), &
    intent(in) :: &
      corr_array_1_n, & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
      corr_array_2_n    ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim,pdf_dim), &
    intent(out) :: &
      corr_array_1, & ! Correlation array of PDF vars. (comp. 1)             [-]
      corr_array_2    ! Correlation array of PDF vars. (comp. 2)             [-]

    ! Local Variables
    integer :: ivar, jvar ! Loop indices

    ! ----------- Begin Code -----------

    ! Initialize diagonal elements to one
    do ivar = 1, pdf_dim
      corr_array_1(:,:,ivar,ivar) = one
      corr_array_2(:,:,ivar,ivar) = one
    end do

    ! The correlations in each PDF component between two of w, chi (old s), and
    ! eta (old t) do not need to be transformed to standard space, since w, chi,
    ! and eta follow assumed normal distributions in each PDF component.  The
    ! normal space correlations between any two of these variables are the same
    ! as the actual correlations.    
    corr_array_1(:,:,iiPDF_eta,iiPDF_chi) = corr_array_1_n(:,:,iiPDF_eta,iiPDF_chi)
    corr_array_1(:,:,iiPDF_w,iiPDF_chi)   = corr_array_1_n(:,:,iiPDF_w,iiPDF_chi)
    corr_array_1(:,:,iiPDF_w,iiPDF_eta)   = corr_array_1_n(:,:,iiPDF_w,iiPDF_eta)
    
    corr_array_2(:,:,iiPDF_eta,iiPDF_chi) = corr_array_2_n(:,:,iiPDF_eta,iiPDF_chi)
    corr_array_2(:,:,iiPDF_w,iiPDF_chi)   = corr_array_2_n(:,:,iiPDF_w,iiPDF_chi)
    corr_array_2(:,:,iiPDF_w,iiPDF_eta)   = corr_array_2_n(:,:,iiPDF_w,iiPDF_eta)

    !!! Calculate the true correlation of variables that have an assumed normal
    !!! distribution and variables that have an assumed lognormal distribution
    !!! for the ith PDF component, given their normal space correlation and the
    !!! normal space standard deviation of the variable with the assumed
    !!! lognormal distribution.

    ! Transform the correlations between chi/eta/w and N_cn to standard space.

    ! Transform the correlation of chi (old s) and N_cn to standard space in PDF component 1.
    call corr_NN2NL( nz, ngrdcol, & ! intent(in)
                     corr_array_1_n(:,:,iiPDF_Ncn,iiPDF_chi), & ! intent(in)
                     sigma_x_1_n(:,:,iiPDF_Ncn), sigma2_on_mu2_ip_1(:,:,iiPDF_Ncn), & ! intent(in)
                     corr_array_1(:,:,iiPDF_Ncn,iiPDF_chi)  ) ! intent(out)
                        
    ! Transform the correlation of eta (old t) and N_cn to standard space in PDF component 1.
    call corr_NN2NL( nz, ngrdcol, & ! intent(in)
                     corr_array_1_n(:,:,iiPDF_Ncn,iiPDF_eta), & ! intent(in)
                     sigma_x_1_n(:,:,iiPDF_Ncn), sigma2_on_mu2_ip_1(:,:,iiPDF_Ncn), & ! intent(in)
                     corr_array_1(:,:,iiPDF_Ncn,iiPDF_eta) )! intent(out)
 
    ! Transform the correlation of w and N_cn to standard space in PDF component 1.
    call corr_NN2NL( nz, ngrdcol, & ! intent(in)
                     corr_array_1_n(:,:,iiPDF_Ncn,iiPDF_w), & ! intent(in)
                     sigma_x_1_n(:,:,iiPDF_Ncn), sigma2_on_mu2_ip_1(:,:,iiPDF_Ncn), & ! intent(in)
                     corr_array_1(:,:,iiPDF_Ncn,iiPDF_w)  ) ! intent(out)
 
    ! Transform the correlation of chi (old s) and N_cn to standard space in PDF component 2.
    call corr_NN2NL( nz, ngrdcol, & ! intent(in)
                     corr_array_2_n(:,:,iiPDF_Ncn,iiPDF_chi), & ! intent(in)
                     sigma_x_2_n(:,:,iiPDF_Ncn), sigma2_on_mu2_ip_1(:,:,iiPDF_Ncn), & ! intent(in)
                     corr_array_2(:,:,iiPDF_Ncn,iiPDF_chi)  ) ! intent(out)

    ! Transform the correlation of eta (old t) and N_cn to standard space in PDF component 2.
    call corr_NN2NL( nz, ngrdcol, & ! intent(in)
                     corr_array_2_n(:,:,iiPDF_Ncn,iiPDF_eta), & ! intent(in)
                     sigma_x_2_n(:,:,iiPDF_Ncn), sigma2_on_mu2_ip_1(:,:,iiPDF_Ncn), & ! intent(in)
                     corr_array_2(:,:,iiPDF_Ncn,iiPDF_eta)  ) ! intent(out)
                        
    ! Transform the correlation of w and N_cn to standard space in PDF component 2.
    call corr_NN2NL( nz, ngrdcol, & ! intent(in)
                     corr_array_2_n(:,:,iiPDF_Ncn,iiPDF_w), & ! intent(in)
                     sigma_x_2_n(:,:,iiPDF_Ncn), sigma2_on_mu2_ip_1(:,:,iiPDF_Ncn), & ! intent(in)
                     corr_array_2(:,:,iiPDF_Ncn,iiPDF_w)  ) ! intent(out)

    ! Transform the correlations (in-precip) between chi/eta/w and the
    ! precipitating hydrometeors to standard space.
    do ivar = iiPDF_chi, iiPDF_w
      do jvar = iiPDF_Ncn+1, pdf_dim

        ! Transform the correlation (in-precip) between w, chi, or eta and a
        ! precipitating hydrometeor, hm, to standard space in PDF component 1.
        call corr_NN2NL( nz, ngrdcol, & ! intent(in)
                         corr_array_1_n(:,:,jvar,ivar), & ! intent(in)
                         sigma_x_1_n(:,:,jvar), sigma2_on_mu2_ip_1(:,:,jvar), & ! intent(in)
                         corr_array_1(:,:,jvar,ivar)  ) ! intent(out)

        ! Transform the correlation (in-precip) between w, chi, or eta and a
        ! precipitating hydrometeor, hm, to standard space in PDF component 2.
        call corr_NN2NL( nz, ngrdcol, & ! intent(in)
                         corr_array_2_n(:,:,jvar,ivar), & ! intent(in)
                         sigma_x_2_n(:,:,jvar), sigma2_on_mu2_ip_2(:,:,jvar), & ! intent(in)
                         corr_array_2(:,:,jvar,ivar)  ) ! intent(out)

      end do ! jvar = iiPDF_Ncn+1, pdf_dim
    end do ! ivar = iiPDF_chi, iiPDF_w


    !!! Calculate the true correlation of two variables that both have an
    !!! assumed lognormal distribution for the ith PDF component, given their
    !!! normal space correlation and both of their normal space standard
    !!! deviations.

    ! Transform the correlations (in-precip) between N_cn and the precipitating
    ! hydrometeors to standard space.
    do jvar = iiPDF_Ncn+1, pdf_dim

      ! Transform the correlation (in-precip) between N_cn and a precipitating
      ! hydrometeor, hm, to standard space in PDF component 1.
      call corr_NN2LL( nz, ngrdcol, & ! intent(in)
                       corr_array_1_n(:,:,jvar,iiPDF_Ncn), & ! intent(in)
                       sigma_x_1_n(:,:,iiPDF_Ncn), sigma_x_1_n(:,:,jvar), & ! intent(in)
                       sigma2_on_mu2_ip_1(:,:,iiPDF_Ncn), sigma2_on_mu2_ip_1(:,:,jvar), & ! in
                       corr_array_1(:,:,jvar,iiPDF_Ncn)) ! intent(out)

      ! Transform the correlation (in-precip) between N_cn and a precipitating
      ! hydrometeor, hm, to standard space in PDF component 2.
      call corr_NN2LL( nz, ngrdcol, & ! intent(in)
                       corr_array_2_n(:,:,jvar,iiPDF_Ncn), & ! intent(in)
                       sigma_x_2_n(:,:,iiPDF_Ncn), sigma_x_2_n(:,:,jvar), & ! intent(in)
                       sigma2_on_mu2_ip_1(:,:,iiPDF_Ncn), sigma2_on_mu2_ip_2(:,:,jvar), & ! in
                       corr_array_2(:,:,jvar,iiPDF_Ncn)) ! intent(out)

    end do ! jvar = ivar+1, pdf_dim

    ! Transform the correlations (in-precip) between two precipitating
    ! hydrometeors to standard space.
    do ivar = iiPDF_Ncn+1, pdf_dim-1
      do jvar = ivar+1, pdf_dim

        ! Transform the correlation (in-precip) between two precipitating
        ! hydrometeors (for example, r_r and N_r) to standard space in PDF
        ! component 1.
        call corr_NN2LL( nz, ngrdcol, & ! intent(in)
                         corr_array_1_n(:,:,jvar,ivar), & ! intent(in)
                         sigma_x_1_n(:,:,ivar), sigma_x_1_n(:,:,jvar), & ! intent(in)
                         sigma2_on_mu2_ip_1(:,:,ivar), sigma2_on_mu2_ip_1(:,:,jvar), & ! in
                         corr_array_1(:,:,jvar,ivar)) ! intent(out)

        ! Transform the correlation (in-precip) between two precipitating
        ! hydrometeors (for example, r_r and N_r) to standard space in PDF
        ! component 2.
        call corr_NN2LL( nz, ngrdcol, & ! intent(in)
                         corr_array_2_n(:,:,jvar,ivar), & ! intent(in)
                         sigma_x_2_n(:,:,ivar), sigma_x_2_n(:,:,jvar), & ! intent(in)
                         sigma2_on_mu2_ip_2(:,:,ivar), sigma2_on_mu2_ip_2(:,:,jvar), & ! in
                         corr_array_2(:,:,jvar,ivar)) ! intent(out)

      end do ! jvar = ivar+1, pdf_dim
    end do ! ivar = iiPDF_Ncn+1, pdf_dim-1

    return

  end subroutine denorm_transform_corr

  !=============================================================================
  elemental subroutine calc_corr_w_hm_n( wm, wphydrometp, &
                                         mu_w_1, mu_w_2, &
                                         mu_hm_1, mu_hm_2, &
                                         sigma_w_1, sigma_w_2, &
                                         sigma_hm_1, sigma_hm_2, &
                                         sigma_hm_1_n, sigma_hm_2_n, &
                                         mixt_frac, precip_frac_1, &
                                         precip_frac_2, hm_tol, &
                                         corr_w_hm_1_n, corr_w_hm_2_n )

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
      precip_frac_2    ! Precipitation fraction (2nd PDF component)          [-]

    real( kind = core_rknd ), intent(in) :: &
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
       end if

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
       end if

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
       end if

       ! The PDF component 1 correlation is undefined.
       corr_w_hm_1_n = zero
       

    else    ! sigma_w_1 * sigma_hm_1 = 0 .and. sigma_w_2 * sigma_hm_2 = 0.

       ! At least one of w and hm is constant in both PDF components.

       ! The PDF component 1 and component 2 correlations are both undefined.
       corr_w_hm_1_n = zero
       corr_w_hm_2_n = zero


    end if


    return

  end subroutine calc_corr_w_hm_n

  !=============================================================================
  subroutine pdf_param_hm_stats( nz, pdf_dim, hm_1, hm_2, &
                                 mu_x_1, mu_x_2, &
                                 sigma_x_1, sigma_x_2, &
                                 corr_array_1, corr_array_2, &
                                 l_stats_samp, &
                                 stats_zt )

    ! Description:
    ! Record statistics for standard PDF parameters involving hydrometeors.

    ! References:
    !-----------------------------------------------------------------------

    use index_mapping, only: &
        pdf2hydromet_idx  ! Procedure(s)

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use array_index, only: &
        iiPDF_w,   & ! Variable(s)
        iiPDF_chi, &
        iiPDF_eta, &
        iiPDF_Ncn

    use clubb_precision, only: &
        core_rknd   ! Variable(s)

    use stats_type_utilities, only: &
!        stat_update_var  ! Procedure(s)
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
        icorr_w_chi_1_ca,   & ! Variable(s)
        icorr_w_chi_2_ca,   &
        icorr_w_eta_1_ca,   &
        icorr_w_eta_2_ca,   &
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
        icorr_hmx_hmy_2

    use stats_type, only: stats ! Type

    implicit none

    type (stats), target, intent(inout) :: &
      stats_zt

    ! Input Variables
    integer, intent(in) :: &
      nz,      & ! Number of vertical levels
      pdf_dim    ! Number of variables in the correlation array

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hm_1, & ! Mean of a precip. hydrometeor (1st PDF component)   [units vary]
      hm_2    ! Mean of a precip. hydrometeor (2nd PDF component)   [units vary]

    real( kind = core_rknd ), dimension(nz,pdf_dim), intent(in) :: &
      mu_x_1,    & ! Mean array of PDF vars. (1st PDF component)    [units vary]
      mu_x_2,    & ! Mean array of PDF vars. (2nd PDF component)    [units vary]
      sigma_x_1, & ! Standard deviation array of PDF vars (comp. 1) [units vary]
      sigma_x_2    ! Standard deviation array of PDF vars (comp. 2) [units vary]

    real( kind = core_rknd ), dimension(nz,pdf_dim,pdf_dim), &
    intent(in) :: &
      corr_array_1, & ! Correlation array of PDF vars. (comp. 1)             [-]
      corr_array_2    ! Correlation array of PDF vars. (comp. 2)             [-]

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.

    ! Local Variable
    integer :: ivar, jvar, hm_idx, hm_idx_ivar, hm_idx_jvar, k  ! Indices


    !!! Output the statistics for hydrometeor PDF parameters.

    ! Statistics
    if ( l_stats_samp ) then

       ! Switch back to using stat_update_var once the code is generalized
       ! to pass in the number of vertical levels.
       do ivar = 1, hydromet_dim, 1

          if ( ihm_1(ivar) > 0 ) then
             ! Mean of the precipitating hydrometeor in PDF component 1.
!             call stat_update_var( ihm_1(ivar), hm_1(:,ivar), stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( ihm_1(ivar), k, hm_1(k,ivar), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

          if ( ihm_2(ivar) > 0 ) then
             ! Mean of the precipitating hydrometeor in PDF component 2.
!             call stat_update_var( ihm_2(ivar), hm_2(:,ivar), stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( ihm_2(ivar), k, hm_2(k,ivar), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

       end do ! ivar = 1, hydromet_dim, 1

       do ivar = iiPDF_Ncn+1, pdf_dim, 1

          hm_idx = pdf2hydromet_idx(ivar)

          ! Mean of the precipitating hydrometeor (in-precip)
          ! in PDF component 1.
          if ( imu_hm_1(hm_idx) > 0 ) then
!             call stat_update_var( imu_hm_1(hm_idx), mu_x_1(ivar,:), stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( imu_hm_1(hm_idx), k, mu_x_1(k,ivar), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

          ! Mean of the precipitating hydrometeor (in-precip)
          ! in PDF component 2.
          if ( imu_hm_2(hm_idx) > 0 ) then
!             call stat_update_var( imu_hm_2(hm_idx), mu_x_2(ivar,:), stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( imu_hm_2(hm_idx), k, mu_x_2(k,ivar), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

       end do ! ivar = iiPDF_Ncn+1, pdf_dim, 1

       ! Mean of cloud nuclei concentration in PDF component 1.
       if ( imu_Ncn_1 > 0 ) then
!          call stat_update_var( imu_Ncn_1, mu_x_1(iiPDF_Ncn,:), stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( imu_Ncn_1, k, mu_x_1(k,iiPDF_Ncn), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       ! Mean of cloud nuclei concentration in PDF component 2.
       if ( imu_Ncn_2 > 0 ) then
!          call stat_update_var( imu_Ncn_2, mu_x_2(iiPDF_Ncn,:), stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( imu_Ncn_2, k, mu_x_2(k,iiPDF_Ncn), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       do ivar = iiPDF_Ncn+1, pdf_dim, 1

          hm_idx = pdf2hydromet_idx(ivar)

          ! Standard deviation of the precipitating hydrometeor (in-precip)
          ! in PDF component 1.
          if ( isigma_hm_1(hm_idx) > 0 ) then
!             call stat_update_var( isigma_hm_1(hm_idx), sigma_x_1(ivar,:), &
!                                   stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( isigma_hm_1(hm_idx), k, sigma_x_1(k,ivar), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

          ! Standard deviation of the precipitating hydrometeor (in-precip)
          ! in PDF component 2.
          if ( isigma_hm_2(hm_idx) > 0 ) then
!             call stat_update_var( isigma_hm_2(hm_idx), sigma_x_2(ivar,:), &
!                                   stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( isigma_hm_2(hm_idx), k, sigma_x_2(k,ivar), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

       end do ! ivar = iiPDF_Ncn+1, pdf_dim, 1

       ! Standard deviation of cloud nuclei concentration in PDF component 1.
       if ( isigma_Ncn_1 > 0 ) then
!          call stat_update_var( isigma_Ncn_1, sigma_x_1(iiPDF_Ncn,:), stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( isigma_Ncn_1, k, sigma_x_1(k,iiPDF_Ncn), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       ! Standard deviation of cloud nuclei concentration in PDF component 2.
       if ( isigma_Ncn_2 > 0 ) then
!          call stat_update_var( isigma_Ncn_2, sigma_x_2(iiPDF_Ncn,:), stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( isigma_Ncn_2, k, sigma_x_2(k,iiPDF_Ncn), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       ! Correlation of w and chi (old s) in PDF component 1 found in the
       ! correlation array.
       ! The true correlation of w and chi in each PDF component is solved for
       ! by an equation and is part of CLUBB's PDF parameters.  However, there
       ! is an option in CLUBB, l_fix_w_chi_eta_correlations, that sets the
       ! component correlation of w and chi to a constant, prescribed value
       ! because of SILHS.  The correlation of w and chi in PDF component 1
       ! that is calculated by an equation is stored in stats as "corr_w_chi_1".
       ! Here, "corr_w_chi_1_ca" outputs whatever value is found in the
       ! correlation array, whether or not it matches "corr_w_chi_1".
       if ( icorr_w_chi_1_ca > 0 ) then
!          call stat_update_var( icorr_w_chi_1_ca, &
!                                corr_array_1(iiPDF_w,iiPDF_chi,:), stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( icorr_w_chi_1_ca, k, &
                                      corr_array_1(k,iiPDF_w,iiPDF_chi), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       ! Correlation of w and chi (old s) in PDF component 2 found in the
       ! correlation array.
       ! The true correlation of w and chi in each PDF component is solved for
       ! by an equation and is part of CLUBB's PDF parameters.  However, there
       ! is an option in CLUBB, l_fix_w_chi_eta_correlations, that sets the
       ! component correlation of w and chi to a constant, prescribed value
       ! because of SILHS.  The correlation of w and chi in PDF component 2
       ! that is calculated by an equation is stored in stats as "corr_w_chi_2".
       ! Here, "corr_w_chi_2_ca" outputs whatever value is found in the
       ! correlation array, whether or not it matches "corr_w_chi_2".
       if ( icorr_w_chi_2_ca > 0 ) then
!          call stat_update_var( icorr_w_chi_2_ca, &
!                                corr_array_2(iiPDF_w,iiPDF_chi,:), stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( icorr_w_chi_2_ca, k, & ! intent(in)
                                      corr_array_2(k,iiPDF_w,iiPDF_chi), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       ! Correlation of w and eta (old t) in PDF component 1 found in the
       ! correlation array.
       ! The true correlation of w and eta in each PDF component is solved for
       ! by an equation and is part of CLUBB's PDF parameters.  However, there
       ! is an option in CLUBB, l_fix_w_chi_eta_correlations, that sets the
       ! component correlation of w and eta to a constant, prescribed value
       ! because of SILHS.  The correlation of w and eta in PDF component 1
       ! that is calculated by an equation is stored in stats as "corr_w_eta_1".
       ! Here, "corr_w_eta_1_ca" outputs whatever value is found in the
       ! correlation array, whether or not it matches "corr_w_eta_1".
       if ( icorr_w_eta_1_ca > 0 ) then
!          call stat_update_var( icorr_w_eta_1_ca, &
!                                corr_array_1(iiPDF_w,iiPDF_eta,:), stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( icorr_w_eta_1_ca, k, & ! intent(in)
                                      corr_array_1(k,iiPDF_w,iiPDF_eta), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       ! Correlation of w and eta (old t) in PDF component 2 found in the
       ! correlation array.
       ! The true correlation of w and eta in each PDF component is solved for
       ! by an equation and is part of CLUBB's PDF parameters.  However, there
       ! is an option in CLUBB, l_fix_w_chi_eta_correlations, that sets the
       ! component correlation of w and eta to a constant, prescribed value
       ! because of SILHS.  The correlation of w and eta in PDF component 2
       ! that is calculated by an equation is stored in stats as "corr_w_eta_2".
       ! Here, "corr_w_eta_2_ca" outputs whatever value is found in the
       ! correlation array, whether or not it matches "corr_w_eta_2".
       if ( icorr_w_eta_2_ca > 0 ) then
!          call stat_update_var( icorr_w_eta_2_ca, &
!                                corr_array_2(iiPDF_w,iiPDF_eta,:), stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( icorr_w_eta_2_ca, k, & ! intent(in)
                                      corr_array_2(k,iiPDF_w,iiPDF_eta), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       do ivar = iiPDF_Ncn+1, pdf_dim, 1

          hm_idx = pdf2hydromet_idx(ivar)

          ! Correlation (in-precip) of w and the precipitating hydrometeor
          ! in PDF component 1.
          if ( icorr_w_hm_1(hm_idx) > 0 ) then
!             call stat_update_var( icorr_w_hm_1(hm_idx), &
!                                   corr_array_1(ivar,iiPDF_w,:), stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( icorr_w_hm_1(hm_idx), k, & ! intent(in)
                                         corr_array_1(k,ivar,iiPDF_w), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

          ! Correlation (in-precip) of w and the precipitating hydrometeor
          ! in PDF component 2.
          if ( icorr_w_hm_2(hm_idx) > 0 ) then
!             call stat_update_var( icorr_w_hm_2(hm_idx), &
!                                   corr_array_2(ivar,iiPDF_w,:), stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( icorr_w_hm_2(hm_idx), k, & ! intent(in)
                                         corr_array_2(k,ivar,iiPDF_w), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

       end do ! ivar = iiPDF_Ncn+1, pdf_dim, 1

       ! Correlation of w and N_cn in PDF component 1.
       if ( icorr_w_Ncn_1 > 0 ) then
!          call stat_update_var( icorr_w_Ncn_1, &
!                                corr_array_1(iiPDF_Ncn,iiPDF_w,:), stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( icorr_w_Ncn_1, k, & ! intent(in)
                                      corr_array_1(k,iiPDF_Ncn,iiPDF_w), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       ! Correlation of w and N_cn in PDF component 2.
       if ( icorr_w_Ncn_2 > 0 ) then
!          call stat_update_var( icorr_w_Ncn_2, &
!                                corr_array_2(iiPDF_Ncn,iiPDF_w,:), stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( icorr_w_Ncn_2, k, & ! intent(in)
                                      corr_array_2(k,iiPDF_Ncn,iiPDF_w), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       ! Correlation of chi (old s) and eta (old t) in PDF component 1 found in
       ! the correlation array.
       ! The true correlation of chi and eta in each PDF component is solved for
       ! by an equation and is part of CLUBB's PDF parameters.  However, there
       ! is an option in CLUBB, l_fix_w_chi_eta_correlations, that sets the
       ! component correlation of chi and eta to a constant, prescribed value
       ! because of SILHS.  The correlation of chi and eta in PDF component 1
       ! that is calculated by an equation is stored in stats as
       ! "corr_chi_eta_1".  Here, "corr_chi_eta_1_ca" outputs whatever value is
       ! found in the correlation array, whether or not it matches
       ! "corr_chi_eta_1".
       if ( icorr_chi_eta_1_ca > 0 ) then
!          call stat_update_var( icorr_chi_eta_1_ca, &
!                                corr_array_1(iiPDF_eta,iiPDF_chi,:), stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( icorr_chi_eta_1_ca, k, & ! intent(in)
                                      corr_array_1(k,iiPDF_eta,iiPDF_chi), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       ! Correlation of chi (old s) and eta (old t) in PDF component 2 found in
       ! the correlation array.
       ! The true correlation of chi and eta in each PDF component is solved for
       ! by an equation and is part of CLUBB's PDF parameters.  However, there
       ! is an option in CLUBB, l_fix_w_chi_eta_correlations, that sets the
       ! component correlation of chi and eta to a constant, prescribed value
       ! because of SILHS.  The correlation of chi and eta in PDF component 2
       ! that is calculated by an equation is stored in stats as
       ! "corr_chi_eta_2".  Here, "corr_chi_eta_2_ca" outputs whatever value is
       ! found in the correlation array, whether or not it matches
       ! "corr_chi_eta_2".
       if ( icorr_chi_eta_2_ca > 0 ) then
!          call stat_update_var( icorr_chi_eta_2_ca, &
!                                corr_array_2(iiPDF_eta,iiPDF_chi,:), stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( icorr_chi_eta_2_ca, k, & ! intent(in)
                                      corr_array_2(k,iiPDF_eta,iiPDF_chi), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       do ivar = iiPDF_Ncn+1, pdf_dim, 1

          hm_idx = pdf2hydromet_idx(ivar)

          ! Correlation (in-precip) of chi (old s) and the precipitating
          ! hydrometeor in PDF component 1.
          if ( icorr_chi_hm_1(hm_idx) > 0 ) then
!             call stat_update_var( icorr_chi_hm_1(hm_idx), &
!                                   corr_array_1(ivar,iiPDF_chi,:), stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( icorr_chi_hm_1(hm_idx), k, & ! intent(in)
                                         corr_array_1(k,ivar,iiPDF_chi), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

          ! Correlation (in-precip) of chi (old s) and the precipitating
          ! hydrometeor in PDF component 2.
          if ( icorr_chi_hm_2(hm_idx) > 0 ) then
!             call stat_update_var( icorr_chi_hm_2(hm_idx), &
!                                   corr_array_2(ivar,iiPDF_chi,:), stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( icorr_chi_hm_2(hm_idx), k, & ! intent(in)
                                         corr_array_2(k,ivar,iiPDF_chi), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

       end do ! ivar = iiPDF_Ncn+1, pdf_dim, 1

       ! Correlation of chi (old s) and N_cn in PDF component 1.
       if ( icorr_chi_Ncn_1 > 0 ) then
!          call stat_update_var( icorr_chi_Ncn_1, &
!                                corr_array_1(iiPDF_Ncn,iiPDF_chi,:), stats_zt )
          do k = 1, nz, 1 
             call stat_update_var_pt( icorr_chi_Ncn_1, k, & ! intent(in)
                                      corr_array_1(k,iiPDF_Ncn,iiPDF_chi), & ! intent(inout)
                                      stats_zt ) ! intent(in)
          end do ! k = 1, nz, 1
       end if

       ! Correlation of chi (old s) and N_cn in PDF component 2.
       if ( icorr_chi_Ncn_2 > 0 ) then
!          call stat_update_var( icorr_chi_Ncn_2, &
!                                corr_array_2(iiPDF_Ncn,iiPDF_chi,:), stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( icorr_chi_Ncn_2, k, & ! intent(in)
                                      corr_array_2(k,iiPDF_Ncn,iiPDF_chi), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       do ivar = iiPDF_Ncn+1, pdf_dim, 1

          hm_idx = pdf2hydromet_idx(ivar)

          ! Correlation (in-precip) of eta (old t) and the precipitating
          ! hydrometeor in PDF component 1.
          if ( icorr_eta_hm_1(hm_idx) > 0 ) then
!             call stat_update_var( icorr_eta_hm_1(hm_idx), &
!                                   corr_array_1(ivar,iiPDF_eta,:), stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( icorr_eta_hm_1(hm_idx), k, & ! intent(in)
                                         corr_array_1(k,ivar,iiPDF_eta), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

          ! Correlation (in-precip) of eta (old t) and the precipitating
          ! hydrometeor in PDF component 2.
          if ( icorr_eta_hm_2(hm_idx) > 0 ) then
!             call stat_update_var( icorr_eta_hm_2(hm_idx), &
!                                   corr_array_2(ivar,iiPDF_eta,:), stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( icorr_eta_hm_2(hm_idx), k, & ! intent(in)
                                         corr_array_2(k,ivar,iiPDF_eta), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

       end do ! ivar = iiPDF_Ncn+1, pdf_dim, 1

       ! Correlation of eta (old t) and N_cn in PDF component 1.
       if ( icorr_eta_Ncn_1 > 0 ) then
!          call stat_update_var( icorr_eta_Ncn_1, &
!                                corr_array_1(iiPDF_Ncn,iiPDF_eta,:), stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( icorr_eta_Ncn_1, k, & ! intent(in)
                                      corr_array_1(k,iiPDF_Ncn,iiPDF_eta), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       ! Correlation of eta (old t) and N_cn in PDF component 2.
       if ( icorr_eta_Ncn_2 > 0 ) then
!          call stat_update_var( icorr_eta_Ncn_2, &
!                                corr_array_2(iiPDF_Ncn,iiPDF_eta,:), stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( icorr_eta_Ncn_2, k, & ! intent(in)
                                      corr_array_2(k,iiPDF_Ncn,iiPDF_eta), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       do ivar = iiPDF_Ncn+1, pdf_dim, 1

          hm_idx = pdf2hydromet_idx(ivar)

          ! Correlation (in-precip) of N_cn and the precipitating
          ! hydrometeor in PDF component 1.
          if ( icorr_Ncn_hm_1(hm_idx) > 0 ) then
!             call stat_update_var( icorr_Ncn_hm_1(hm_idx), &
!                                   corr_array_1(ivar,iiPDF_Ncn,:), stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( icorr_Ncn_hm_1(hm_idx), k, & ! intent(in)
                                         corr_array_1(k,ivar,iiPDF_Ncn), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

          ! Correlation (in-precip) of N_cn and the precipitating
          ! hydrometeor in PDF component 2.
          if ( icorr_Ncn_hm_2(hm_idx) > 0 ) then
!             call stat_update_var( icorr_Ncn_hm_2(hm_idx), &
!                                   corr_array_2(ivar,iiPDF_Ncn,:), stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( icorr_Ncn_hm_2(hm_idx), k, & ! intent(in)
                                         corr_array_2(k,ivar,iiPDF_Ncn), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

       end do ! ivar = iiPDF_Ncn+1, pdf_dim, 1

       do ivar = iiPDF_Ncn+1, pdf_dim, 1

          hm_idx_ivar = pdf2hydromet_idx(ivar)

          do jvar = ivar+1, pdf_dim, 1

             hm_idx_jvar = pdf2hydromet_idx(jvar)

             ! Correlation (in-precip) of two different hydrometeors (hmx and
             ! hmy) in PDF component 1.
             if ( icorr_hmx_hmy_1(hm_idx_jvar,hm_idx_ivar) > 0 ) then
!               call stat_update_var( icorr_hmx_hmy_1(hm_idx_jvar,hm_idx_ivar), &
!                                     corr_array_1(jvar,ivar,:), stats_zt )
               do k = 1, nz, 1
                  call stat_update_var_pt( icorr_hmx_hmy_1(hm_idx_jvar,hm_idx_ivar), k, & ! in
                                           corr_array_1(k,jvar,ivar), & ! intent(in)
                                           stats_zt ) ! intent(inout)
               end do ! k = 1, nz, 1
             end if

             ! Correlation (in-precip) of two different hydrometeors (hmx and
             ! hmy) in PDF component 2.
             if ( icorr_hmx_hmy_2(hm_idx_jvar,hm_idx_ivar) > 0 ) then
!               call stat_update_var( icorr_hmx_hmy_2(hm_idx_jvar,hm_idx_ivar), &
!                                     corr_array_2(jvar,ivar,:), stats_zt )
               do k = 1, nz, 1
                  call stat_update_var_pt( icorr_hmx_hmy_2(hm_idx_jvar,hm_idx_ivar), k, & ! in
                                           corr_array_2(k,jvar,ivar), & ! intent(in)
                                           stats_zt ) ! intent(inout)
               end do ! k = 1, nz, 1
             end if

          end do ! jvar = ivar+1, pdf_dim, 1
       end do ! ivar = iiPDF_Ncn+1, pdf_dim, 1

    end if ! l_stats_samp


    return

  end subroutine pdf_param_hm_stats

  !=============================================================================
  subroutine pdf_param_ln_hm_stats( nz, pdf_dim, mu_x_1_n, &
                                    mu_x_2_n, sigma_x_1_n, &
                                    sigma_x_2_n, corr_array_1_n, &
                                    corr_array_2_n, l_stats_samp, & 
                                    stats_zt )

    ! Description:
    ! Record statistics for normal space PDF parameters involving hydrometeors.

    ! References:
    !-----------------------------------------------------------------------

    use index_mapping, only: &
        pdf2hydromet_idx  ! Procedure(s)

    use array_index, only: &
        iiPDF_w,   & ! Variable(s)
        iiPDF_chi, &
        iiPDF_eta, &
        iiPDF_Ncn

    use clubb_precision, only: &
        core_rknd   ! Variable(s)

    use stats_type_utilities, only: &
!        stat_update_var  ! Procedure(s)
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
        icorr_hmx_hmy_2_n

    use stats_type, only: stats ! Type

    implicit none

    type (stats), target, intent(inout) :: &
      stats_zt

    ! Input Variables
    integer, intent(in) :: &
      nz,      & ! Number of vertical levels
      pdf_dim    ! Number of variables in the correlation array

    real( kind = core_rknd ), dimension(nz,pdf_dim), intent(in) :: &
      mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(nz,pdf_dim,pdf_dim), &
    intent(in) :: &
      corr_array_1_n, & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
      corr_array_2_n    ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.

    ! Local Variable
    real( kind = core_rknd ), dimension(nz) :: &
      mu_hm_1_n,  & ! Mean of ln hm (1st PDF component)    [units vary]
      mu_hm_2_n,  & ! Mean of ln hm (2nd PDF component)    [units vary]
      mu_Ncn_1_n, & ! Mean of ln Ncn (1st PDF component)   [ln(num/kg)]
      mu_Ncn_2_n    ! Mean of ln Ncn (2nd PDF component)   [ln(num/kg)]

    integer :: ivar, jvar, hm_idx, hm_idx_ivar, hm_idx_jvar, k  ! Indices


    !!! Output the statistics for normal space hydrometeor PDF parameters.

    ! Statistics
    if ( l_stats_samp ) then

       ! Switch back to using stat_update_var once the code is generalized
       ! to pass in the number of vertical levels.
       do ivar = iiPDF_Ncn+1, pdf_dim, 1

          hm_idx = pdf2hydromet_idx(ivar)

          ! Mean (in-precip) of ln hm in PDF component 1.
          if ( imu_hm_1_n(hm_idx) > 0 ) then
             where ( mu_x_1_n(:,ivar) > real( -huge( 0.0 ), kind = core_rknd ) )
                mu_hm_1_n = mu_x_1_n(:,ivar)
             elsewhere
                ! When hm_1 is 0 (or below tolerance value), mu_hm_1_n is -inf,
                ! and is set to -huge for the default CLUBB kind.  Some
                ! compilers have issues outputting to stats files (in single
                ! precision) when the default CLUBB kind is in double precision.
                ! Set to -huge for single precision.
                mu_hm_1_n = real( -huge( 0.0 ), kind = core_rknd )
             endwhere
!             call stat_update_var( imu_hm_1_n(hm_idx), mu_hm_1_n, stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( imu_hm_1_n(hm_idx), k, mu_hm_1_n(k), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

          ! Mean (in-precip) of ln hm in PDF component 2.
          if ( imu_hm_2_n(hm_idx) > 0 ) then
             where ( mu_x_2_n(:,ivar) > real( -huge( 0.0 ), kind = core_rknd ) )
                mu_hm_2_n = mu_x_2_n(:,ivar)
             elsewhere
                ! When hm_2 is 0 (or below tolerance value), mu_hm_2_n is -inf,
                ! and is set to -huge for the default CLUBB kind.  Some
                ! compilers have issues outputting to stats files (in single
                ! precision) when the default CLUBB kind is in double precision.
                ! Set to -huge for single precision.
                mu_hm_2_n = real( -huge( 0.0 ), kind = core_rknd )
             endwhere
!             call stat_update_var( imu_hm_2_n(hm_idx), mu_hm_2_n, stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( imu_hm_2_n(hm_idx), k, mu_hm_2_n(k), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

       end do ! ivar = iiPDF_Ncn+1, pdf_dim, 1

       ! Mean of ln N_cn in PDF component 1.
       if ( imu_Ncn_1_n > 0 ) then
          where ( mu_x_1_n(:,iiPDF_Ncn) &
                  > real( -huge( 0.0 ), kind = core_rknd ) )
             mu_Ncn_1_n = mu_x_1_n(:,iiPDF_Ncn)
          elsewhere
             ! When Ncnm is 0 (or below tolerance value), mu_Ncn_1_n is -inf,
             ! and is set to -huge for the default CLUBB kind.  Some compilers
             ! have issues outputting to stats files (in single precision) when
             ! the default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             mu_Ncn_1_n = real( -huge( 0.0 ), kind = core_rknd )
          endwhere
!          call stat_update_var( imu_Ncn_1_n, mu_Ncn_1_n, stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( imu_Ncn_1_n, k, mu_Ncn_1_n(k), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       ! Mean of ln N_cn in PDF component 2.
       if ( imu_Ncn_2_n > 0 ) then
          where ( mu_x_2_n(:,iiPDF_Ncn) &
                  > real( -huge( 0.0 ), kind = core_rknd ) )
             mu_Ncn_2_n = mu_x_2_n(:,iiPDF_Ncn)
          elsewhere
             ! When Ncnm is 0 (or below tolerance value), mu_Ncn_2_n is -inf,
             ! and is set to -huge for the default CLUBB kind.  Some compilers
             ! have issues outputting to stats files (in single precision) when
             ! the default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             mu_Ncn_2_n = real( -huge( 0.0 ), kind = core_rknd )
          endwhere
!          call stat_update_var( imu_Ncn_2_n, mu_Ncn_2_n, stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( imu_Ncn_2_n, k, mu_Ncn_2_n(k), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       do ivar = iiPDF_Ncn+1, pdf_dim, 1

          hm_idx = pdf2hydromet_idx(ivar)

          ! Standard deviation (in-precip) of ln hm in PDF component 1.
          if ( isigma_hm_1_n(hm_idx) > 0 ) then
!             call stat_update_var( isigma_hm_1_n(hm_idx), sigma_x_1_n(ivar,:), &
!                                   stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( isigma_hm_1_n(hm_idx), k, & ! intent(in)
                                         sigma_x_1_n(k,ivar), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

          ! Standard deviation (in-precip) of ln hm in PDF component 2.
          if ( isigma_hm_2_n(hm_idx) > 0 ) then
!             call stat_update_var( isigma_hm_2_n(hm_idx), sigma_x_2_n(ivar,:), &
!                                   stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( isigma_hm_2_n(hm_idx), k, & ! intent(in)
                                         sigma_x_2_n(k,ivar), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

       end do ! ivar = iiPDF_Ncn+1, pdf_dim, 1

       ! Standard deviation of ln N_cn in PDF component 1.
       if ( isigma_Ncn_1_n > 0 ) then
!          call stat_update_var( isigma_Ncn_1_n, sigma_x_1_n(iiPDF_Ncn,:), &
!                                stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( isigma_Ncn_1_n, k, & ! intent(in)
                                      sigma_x_1_n(k,iiPDF_Ncn), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       ! Standard deviation of ln N_cn in PDF component 2.
       if ( isigma_Ncn_2_n > 0 ) then
!          call stat_update_var( isigma_Ncn_2_n, sigma_x_2_n(iiPDF_Ncn,:), &
!                                stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( isigma_Ncn_2_n, k, & ! intent(in)
                                      sigma_x_2_n(k,iiPDF_Ncn), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       do ivar = iiPDF_Ncn+1, pdf_dim, 1

          hm_idx = pdf2hydromet_idx(ivar)

          ! Correlation (in-precip) of w and ln hm in PDF component 1.
          if ( icorr_w_hm_1_n(hm_idx) > 0 ) then
!             call stat_update_var( icorr_w_hm_1_n(hm_idx), &
!                                   corr_array_1_n(ivar,iiPDF_w,:), stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( icorr_w_hm_1_n(hm_idx), k, & ! intent(in)
                                         corr_array_1_n(k,ivar,iiPDF_w), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

          ! Correlation (in-precip) of w and ln hm in PDF component 2.
          if ( icorr_w_hm_2_n(hm_idx) > 0 ) then
!             call stat_update_var( icorr_w_hm_2_n(hm_idx), &
!                                   corr_array_2_n(ivar,iiPDF_w,:), stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( icorr_w_hm_2_n(hm_idx), k, & ! intent(in)
                                         corr_array_2_n(k,ivar,iiPDF_w), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

       end do ! ivar = iiPDF_Ncn+1, pdf_dim, 1

       ! Correlation of w and ln N_cn in PDF component 1.
       if ( icorr_w_Ncn_1_n > 0 ) then
!          call stat_update_var( icorr_w_Ncn_1_n, &
!                                corr_array_1_n(iiPDF_Ncn,iiPDF_w,:), stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( icorr_w_Ncn_1_n, k, & ! intent(in)
                                      corr_array_1_n(k,iiPDF_Ncn,iiPDF_w), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       ! Correlation of w and ln N_cn in PDF component 2.
       if ( icorr_w_Ncn_2_n > 0 ) then
!          call stat_update_var( icorr_w_Ncn_2_n, &
!                                corr_array_2_n(iiPDF_Ncn,iiPDF_w,:), stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( icorr_w_Ncn_2_n, k, & ! intent(in)
                                      corr_array_2_n(k,iiPDF_Ncn,iiPDF_w), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       do ivar = iiPDF_Ncn+1, pdf_dim, 1

          hm_idx = pdf2hydromet_idx(ivar)

          ! Correlation (in-precip) of chi (old s) and ln hm in PDF component 1.
          if ( icorr_chi_hm_1_n(hm_idx) > 0 ) then
!             call stat_update_var( icorr_chi_hm_1_n(hm_idx), &
!                                   corr_array_1_n(ivar,iiPDF_chi,:), stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( icorr_chi_hm_1_n(hm_idx), k, & ! intent(in)
                                         corr_array_1_n(k,ivar,iiPDF_chi), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

          ! Correlation (in-precip) of chi( old s) and ln hm in PDF component 2.
          if ( icorr_chi_hm_2_n(hm_idx) > 0 ) then
!             call stat_update_var( icorr_chi_hm_2_n(hm_idx), &
!                                   corr_array_2_n(ivar,iiPDF_chi,:), stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( icorr_chi_hm_2_n(hm_idx), k, & ! intent(in)
                                         corr_array_2_n(k,ivar,iiPDF_chi), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

       end do ! ivar = iiPDF_Ncn+1, pdf_dim, 1

       ! Correlation of chi (old s) and ln N_cn in PDF component 1.
       if ( icorr_chi_Ncn_1_n > 0 ) then
!          call stat_update_var( icorr_chi_Ncn_1_n, &
!                                corr_array_1_n(iiPDF_Ncn,iiPDF_chi,:), &
!                                stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( icorr_chi_Ncn_1_n, k, & ! intent(in)
                                      corr_array_1_n(k,iiPDF_Ncn,iiPDF_chi), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1 
       end if

       ! Correlation of chi(old s) and ln N_cn in PDF component 2.
       if ( icorr_chi_Ncn_2_n > 0 ) then
!          call stat_update_var( icorr_chi_Ncn_2_n, &
!                                corr_array_2_n(iiPDF_Ncn,iiPDF_chi,:), &
!                                stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( icorr_chi_Ncn_2_n, k, & ! intent(in)
                                      corr_array_2_n(k,iiPDF_Ncn,iiPDF_chi), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       do ivar = iiPDF_Ncn+1, pdf_dim, 1

          hm_idx = pdf2hydromet_idx(ivar)

          ! Correlation (in-precip) of eta (old t) and ln hm in PDF component 1.
          if ( icorr_eta_hm_1_n(hm_idx) > 0 ) then
!             call stat_update_var( icorr_eta_hm_1_n(hm_idx), &
!                                   corr_array_1_n(ivar,iiPDF_eta,:), stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( icorr_eta_hm_1_n(hm_idx), k, & ! intent(in)
                                         corr_array_1_n(k,ivar,iiPDF_eta), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

          ! Correlation (in-precip) of eta (old t) and ln hm in PDF component 2.
          if ( icorr_eta_hm_2_n(hm_idx) > 0 ) then
!             call stat_update_var( icorr_eta_hm_2_n(hm_idx), &
!                                   corr_array_2_n(ivar,iiPDF_eta,:), stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( icorr_eta_hm_2_n(hm_idx), k, & ! intent(in)
                                         corr_array_2_n(k,ivar,iiPDF_eta), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

       end do ! ivar = iiPDF_Ncn+1, pdf_dim, 1

       ! Correlation of eta (old t) and ln N_cn in PDF component 1.
       if ( icorr_eta_Ncn_1_n > 0 ) then
!          call stat_update_var( icorr_eta_Ncn_1_n, &
!                                corr_array_1_n(iiPDF_Ncn,iiPDF_eta,:), &
!                                stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( icorr_eta_Ncn_1_n, k, & ! intent(in)
                                      corr_array_1_n(k,iiPDF_Ncn,iiPDF_eta), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1 
       end if

       ! Correlation of eta (old t) and ln N_cn in PDF component 2.
       if ( icorr_eta_Ncn_2_n > 0 ) then
!          call stat_update_var( icorr_eta_Ncn_2_n, &
!                                corr_array_2_n(iiPDF_Ncn,iiPDF_eta,:), &
!                                stats_zt )
          do k = 1, nz, 1
             call stat_update_var_pt( icorr_eta_Ncn_2_n, k, & ! intent(in) 
                                      corr_array_2_n(k,iiPDF_Ncn,iiPDF_eta), & ! intent(in)
                                      stats_zt ) ! intent(inout)
          end do ! k = 1, nz, 1
       end if

       do ivar = iiPDF_Ncn+1, pdf_dim, 1

          hm_idx = pdf2hydromet_idx(ivar)

          ! Correlation (in-precip) of ln N_cn and ln hm in PDF
          ! component 1.
          if ( icorr_Ncn_hm_1_n(hm_idx) > 0 ) then
!             call stat_update_var( icorr_Ncn_hm_1_n(hm_idx), &
!                                   corr_array_1_n(ivar,iiPDF_Ncn,:), stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( icorr_Ncn_hm_1_n(hm_idx), k, & ! intent(in)
                                         corr_array_1_n(k,ivar,iiPDF_Ncn), & ! intent(in)
                                         stats_zt ) ! intent(inout)
             end do ! k = 1, nz, 1
          end if

          ! Correlation (in-precip) of ln N_cn and ln hm in PDF
          ! component 2.
          if ( icorr_Ncn_hm_2_n(hm_idx) > 0 ) then
!             call stat_update_var( icorr_Ncn_hm_2_n(hm_idx), &
!                                   corr_array_2_n(ivar,iiPDF_Ncn,:), stats_zt )
             do k = 1, nz, 1
                call stat_update_var_pt( icorr_Ncn_hm_2_n(hm_idx), k, & ! intent(in)
                                         corr_array_2_n(k,ivar,iiPDF_Ncn), & ! intent(in)
                                         stats_zt )! intent(inout)
             end do ! k = 1, nz, 1 
          end if

       end do ! ivar = iiPDF_Ncn+1, pdf_dim, 1

       do ivar = iiPDF_Ncn+1, pdf_dim, 1

          hm_idx_ivar = pdf2hydromet_idx(ivar)

          do jvar = ivar+1, pdf_dim, 1

             hm_idx_jvar= pdf2hydromet_idx(jvar)

             ! Correlation (in-precip) of ln hmx and ln hmy (two different
             ! hydrometeors) in PDF component 1.
             if ( icorr_hmx_hmy_1_n(hm_idx_jvar,hm_idx_ivar) > 0 ) then
!                call stat_update_var( &
!                        icorr_hmx_hmy_1_n(hm_idx_jvar,hm_idx_ivar), &
!                        corr_array_1_n(jvar,ivar,:), stats_zt )
                do k = 1, nz, 1
                   call stat_update_var_pt( & 
                        icorr_hmx_hmy_1_n(hm_idx_jvar,hm_idx_ivar), k, & ! intent(in)
                        corr_array_1_n(k,jvar,ivar), & ! intent(in)
                        stats_zt ) ! intent(inout)
                end do ! k = 1, nz, 1
             end if

             ! Correlation (in-precip) of ln hmx and ln hmy (two different
             ! hydrometeors) in PDF component 2.
             if ( icorr_hmx_hmy_2_n(hm_idx_jvar,hm_idx_ivar) > 0 ) then
!                call stat_update_var( &
!                        icorr_hmx_hmy_2_n(hm_idx_jvar,hm_idx_ivar), &
!                        corr_array_2_n(jvar,ivar,:), stats_zt )
                do k = 1, nz, 1
                   call stat_update_var_pt( &
                        icorr_hmx_hmy_2_n(hm_idx_jvar,hm_idx_ivar), k, & ! intent(in)
                        corr_array_2_n(k,jvar,ivar), & ! intent(in)
                        stats_zt ) ! intent(inout)
                end do ! k = 1, nz, 1
             end if

          end do ! jvar = ivar+1, pdf_dim, 1
       end do ! ivar = iiPDF_Ncn+1, pdf_dim, 1

    end if ! l_stats_samp


    return

  end subroutine pdf_param_ln_hm_stats

  !=============================================================================
  subroutine pack_hydromet_pdf_params( nz, ngrdcol, hm_1, hm_2, pdf_dim, mu_x_1, & ! In
                                       mu_x_2, sigma_x_1, sigma_x_2, &             ! In
                                       corr_array_1, corr_array_2, &               ! In
                                       hydromet_pdf_params )                       ! Out

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

    use array_index, only: &
        iiPDF_w,   & ! Variable(s)
        iiPDF_chi, &
        iiPDF_eta, &
        iiPDF_Ncn

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz, &     ! Number of vertical grid levels
      ngrdcol   ! Number of grid columns

    real( kind = core_rknd ), dimension(ngrdcol,nz,hydromet_dim), intent(in) :: &
      hm_1, & ! Mean of a precip. hydrometeor (1st PDF component)  [units vary]
      hm_2    ! Mean of a precip. hydrometeor (2nd PDF component)  [units vary]

    integer, intent(in) :: &
      pdf_dim   ! Number of variables in the mean/stdev arrays

    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim), intent(in) :: &
      mu_x_1,    & ! Mean array of PDF vars. (1st PDF component)    [units vary]
      mu_x_2,    & ! Mean array of PDF vars. (2nd PDF component)    [units vary]
      sigma_x_1, & ! Standard deviation array of PDF vars (comp. 1) [units vary]
      sigma_x_2    ! Standard deviation array of PDF vars (comp. 2) [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim,pdf_dim), intent(in) :: &
      corr_array_1, & ! Correlation array of PDF vars. (comp. 1)    [-]
      corr_array_2    ! Correlation array of PDF vars. (comp. 2)    [-]

    ! Output Variable
    type(hydromet_pdf_parameter), dimension(ngrdcol,nz), intent(out) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]

    ! Local Variables
    integer :: ivar, jvar, pdf_idx  ! Indices


    ! Pack remaining means and standard deviations into hydromet_pdf_params.
    do ivar = 1, hydromet_dim, 1

      pdf_idx = hydromet2pdf_idx(ivar)

      ! Mean of a hydrometeor (overall) in the 1st PDF component.
      hydromet_pdf_params(:,:)%hm_1(ivar) = hm_1(:,:,ivar)
      ! Mean of a hydrometeor (overall) in the 2nd PDF component.
      hydromet_pdf_params(:,:)%hm_2(ivar) = hm_2(:,:,ivar)

      ! Mean of a hydrometeor (in-precip) in the 1st PDF component.
      hydromet_pdf_params(:,:)%mu_hm_1(ivar) = mu_x_1(:,:,pdf_idx)
      ! Mean of a hydrometeor (in-precip) in the 2nd PDF component.
      hydromet_pdf_params(:,:)%mu_hm_2(ivar) = mu_x_2(:,:,pdf_idx)

      ! Standard deviation of a hydrometeor (in-precip) in the
      ! 1st PDF component.
      hydromet_pdf_params(:,:)%sigma_hm_1(ivar) = sigma_x_1(:,:,pdf_idx)
      ! Standard deviation of a hydrometeor (in-precip) in the
      ! 2nd PDF component.
      hydromet_pdf_params(:,:)%sigma_hm_2(ivar) = sigma_x_2(:,:,pdf_idx)

      ! Correlation (in-precip) of w and a hydrometeor in the 1st PDF
      ! component.
      hydromet_pdf_params(:,:)%corr_w_hm_1(ivar) = corr_array_1(:,:,pdf_idx,iiPDF_w)

      ! Correlation (in-precip) of w and a hydrometeor in the 2nd PDF
      ! component.
      hydromet_pdf_params(:,:)%corr_w_hm_2(ivar) = corr_array_2(:,:,pdf_idx,iiPDF_w)

      ! Correlation (in-precip) of chi and a hydrometeor in the 1st PDF
      ! component.
      hydromet_pdf_params(:,:)%corr_chi_hm_1(ivar) &
      = corr_array_1(:,:,pdf_idx,iiPDF_chi)

      ! Correlation (in-precip) of chi and a hydrometeor in the 2nd PDF
      ! component.
      hydromet_pdf_params(:,:)%corr_chi_hm_2(ivar) &
      = corr_array_2(:,:,pdf_idx,iiPDF_chi)

      ! Correlation (in-precip) of eta and a hydrometeor in the 1st PDF
      ! component.
      hydromet_pdf_params(:,:)%corr_eta_hm_1(ivar) &
      = corr_array_1(:,:,pdf_idx,iiPDF_eta)

      ! Correlation (in-precip) of eta and a hydrometeor in the 2nd PDF
      ! component.
      hydromet_pdf_params(:,:)%corr_eta_hm_2(ivar) &
      = corr_array_2(:,:,pdf_idx,iiPDF_eta)

      ! Correlation (in-precip) of two hydrometeors, hmx and hmy, in the 1st
      ! PDF component.
      hydromet_pdf_params(:,:)%corr_hmx_hmy_1(ivar,ivar) = one

      do jvar = ivar+1, hydromet_dim, 1

        hydromet_pdf_params(:,:)%corr_hmx_hmy_1(jvar,ivar) &
        = corr_array_1(:,:,hydromet2pdf_idx(jvar),pdf_idx)

        hydromet_pdf_params(:,:)%corr_hmx_hmy_1(ivar,jvar) &
        = hydromet_pdf_params(:,:)%corr_hmx_hmy_1(jvar,ivar)

      end do ! jvar = ivar+1, hydromet_dim, 1

      ! Correlation (in-precip) of two hydrometeors, hmx and hmy, in the 2nd
      ! PDF component.
      hydromet_pdf_params(:,:)%corr_hmx_hmy_2(ivar,ivar) = one

      do jvar = ivar+1, hydromet_dim, 1

        hydromet_pdf_params(:,:)%corr_hmx_hmy_2(jvar,ivar) &
        = corr_array_2(:,:,hydromet2pdf_idx(jvar),pdf_idx)

        hydromet_pdf_params(:,:)%corr_hmx_hmy_2(ivar,jvar) &
        = hydromet_pdf_params(:,:)%corr_hmx_hmy_2(jvar,ivar)

      end do ! jvar = ivar+1, hydromet_dim, 1

    end do ! ivar = 1, hydromet_dim, 1

    ! Mean of Ncn (overall) in the 1st PDF component.
    hydromet_pdf_params(:,:)%mu_Ncn_1 = mu_x_1(:,:,iiPDF_Ncn)
    ! Mean of Ncn (overall) in the 2nd PDF component.
    hydromet_pdf_params(:,:)%mu_Ncn_2 = mu_x_2(:,:,iiPDF_Ncn)

    ! Standard deviation of Ncn (overall) in the 1st PDF component.
    hydromet_pdf_params(:,:)%sigma_Ncn_1 = sigma_x_1(:,:,iiPDF_Ncn)
    ! Standard deviation of Ncn (overall) in the 2nd PDF component.
    hydromet_pdf_params(:,:)%sigma_Ncn_2 = sigma_x_2(:,:,iiPDF_Ncn)


    return

  end subroutine pack_hydromet_pdf_params

  !=============================================================================
  elemental function compute_rtp2_from_chi( sigma_chi_1, sigma_chi_2,        &
                                            sigma_eta_1, sigma_eta_2,        &
                                            rt_1, rt_2,                      &
                                            crt_1, crt_2,                    &
                                            mixt_frac,                       &   
                                            corr_chi_eta_1, corr_chi_eta_2 ) &
  result( rtp2_zt_from_chi )

    ! Description:
    ! Compute the variance of rt from the distribution of chi and eta. The
    ! resulting variance will be consistent with CLUBB's extended PDF
    ! involving chi and eta, including if l_fix_w_chi_eta_correlations = .true..

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

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      sigma_chi_1,         & ! Standard deviation of chi (1st PDF comp.) [kg/kg]
      sigma_chi_2,         & ! Standard deviation of chi (2nd PDF comp.) [kg/kg]
      sigma_eta_1,         & ! Standard deviation of eta (1st PDF comp.) [kg/kg]
      sigma_eta_2,         & ! Standard deviation of eta (2nd PDF comp.) [kg/kg]
      crt_1,               & ! Coef. of r_t in chi/eta eqns. (1st comp.) [-]
      crt_2,               & ! Coef. of r_t in chi/eta eqns. (2nd comp.) [-]
      rt_1,                & ! Mean of rt (1st PDF component)            [kg/kg]
      rt_2,                & ! Mean of rt (2nd PDF component)            [kg/kg]
      mixt_frac              ! Weight of 1st gaussian PDF component      [-]

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
      rtm, &                 ! Mean of rt (overall) [kg/kg]
      sigma_rt_1_from_chi, & ! Standard deviation of rt (1st PDF comp.)  [kg/kg]
      sigma_rt_2_from_chi    ! Standard deviation of rt (2nd PDF comp.)  [kg/kg]


  !-----------------------------------------------------------------------

    !----- Begin Code -----

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
