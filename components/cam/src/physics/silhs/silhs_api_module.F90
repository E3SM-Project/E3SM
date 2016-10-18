!--------------------------------------------------------------------------------------------------
! $Id: silhs_api_module.F90 8019 2016-03-20 22:17:45Z raut@uwm.edu $
!==================================================================================================
!
!       ######## ########### ###        ###    ###  ########            ###     ######### #########
!     ###    ###    ###     ###        ###    ### ###    ###         ### ###   ###    ###   ###
!    ###           ###     ###        ###    ### ###               ###   ###  ###    ###   ###
!   ##########    ###     ###        ########## ##########       ########### #########    ###
!         ###    ###     ###        ###    ###        ###       ###     ### ###          ###
! ###    ###    ###     ###        ###    ### ###    ###       ###     ### ###          ###
! ######## ########### ########## ###    ###  ########        ###     ### ###       #########
!
! The SILHS API serves as the doorway through which external models can interact with SILHS.
!
!               PLEASE REMEMBER, IF ANY CODE IS CHANGED IN THIS DOCUMENT,
!                   THE CHANGES MUST BE PROPOGATED TO ALL HOST MODELS.
!
!
! Cloud Layers Unified By Binormals (CLUBB) user license 
! agreement.
!
! Thank you for your interest in CLUBB. We work hard to create a
! code that implements the best software engineering practices,
! is supported to the extent allowed by our limited resources,
! and is available without cost to non-commercial users. You may
! use CLUBB if, in return, you abide by these conditions:
!
! 1. Please cite CLUBB in presentations and publications that
!  contain results obtained using CLUBB.
!
! 2. You may not use any part of CLUBB to create or modify
!  another single-column (1D) model that is not called CLUBB.
!  However, you may modify or augment CLUBB or parts of CLUBB if
!  you include "CLUBB" in the name of the resulting single-column
!  model. For example, a user at MIT might modify CLUBB and call
!  the modified version "CLUBB-MIT." Or, for example, a user of
!  the CLM land-surface model might interface CLM to CLUBB and
!  call it "CLM-CLUBB." This naming convention recognizes the
!  contributions of both sets of developers.
!
! 3. You may implement CLUBB as a parameterization in a large-
!  scale host model that has 2 or 3 spatial dimensions without 
!  including "CLUBB" in the combined model name, but please 
!  acknowledge in presentations and publications that CLUBB has 
!  been included as a parameterization.
!
! 4. You may not provide all or part of CLUBB to anyone without 
!  prior permission from Vincent Larson (vlarson@uwm.edu). If 
!  you wish to share CLUBB with your collaborators without 
!  seeking permission, please ask your collaborators to register 
!  as CLUBB users at http://clubb.larson-group.com and to 
!  download CLUBB from there.
!
! 5. You may not use CLUBB for commercial purposes unless you 
!  receive permission from Vincent Larson.
!
! 6. You may not re-license all or any part of CLUBB.
!
! 7. CLUBB is provided "as is" and without warranty.
!
! We hope that CLUBB will develop into a community resource. We 
! encourage users to contribute their CLUBB modifications or 
! extensions to the CLUBB development group. We will then 
! consider them for inclusion in CLUBB. Such contributions will 
! benefit all CLUBB users. We would be pleased to acknowledge 
! contributors and list their CLUBB-related papers on our "About 
! CLUBB" webpage (http://clubb.larson-group.com/about.php) for 
! those contributors who so desire.
!
! Thanks so much and best wishes for your research!
!
! The CLUBB Development Group
! (Present and past contributors to the source code include 
! Vincent Larson, Chris Golaz, David Schanen, Brian Griffin, 
! Joshua Fasching, Adam Smith, and Michael Falk).
!------------------------------------------------------------------

module silhs_api_module

#ifdef SILHS
  use parameters_silhs, only: &
    l_lh_importance_sampling, & ! Variable(s)
    l_Lscale_vert_avg

  use latin_hypercube_driver_module, only: &
    lh_clipped_variables_type ! Type

#endif

  implicit none

  private

#ifdef SILHS

  public  &
    lh_subcolumn_generator_api, &
    stats_accumulate_lh_api, &
    est_kessler_microphys_api, &
    l_lh_importance_sampling, &
    l_Lscale_vert_avg, &
    clip_transform_silhs_output_api, &
    lh_clipped_variables_type, &
    lh_microphys_var_covar_driver_api

contains

  !================================================================================================
  ! lh_subcolumn_generator - Generates sample points of moisture, temperature, et cetera.
  !================================================================================================

  subroutine lh_subcolumn_generator_api( &
    iter, d_variables, num_samples, sequence_length, nz, & ! In
    pdf_params, delta_zm, rcm, Lscale, & ! In
    rho_ds_zt, mu1, mu2, sigma1, sigma2, & ! In
    corr_cholesky_mtx_1, corr_cholesky_mtx_2, & ! In
    hydromet_pdf_params, & ! In
    X_nl_all_levs, X_mixt_comp_all_levs, & ! Out
    lh_sample_point_weights ) ! Out

    use latin_hypercube_driver_module, only : lh_subcolumn_generator

    use pdf_parameter_module, only: &
      pdf_parameter  ! Type

    use hydromet_pdf_parameter_module, only: &
      hydromet_pdf_parameter ! Type

    use clubb_precision, only: &
      core_rknd

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      iter,            & ! Model iteration number
      d_variables,     & ! Number of variables to sample
      num_samples,     & ! Number of samples per variable
      sequence_length, & ! nt_repeat/num_samples; number of timesteps before sequence repeats.
      nz                 ! Number of vertical model levels

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params ! PDF parameters       [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      delta_zm, &  ! Difference in moment. altitudes    [m]
      rcm          ! Liquid water mixing ratio          [kg/kg]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rho_ds_zt    ! Dry, static density on thermo. levels    [kg/m^3]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      Lscale       ! Turbulent Mixing Length  [m]

    ! Output Variables
    real( kind = core_rknd ), intent(out), dimension(nz,num_samples,d_variables) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    integer, intent(out), dimension(nz,num_samples) :: &
      X_mixt_comp_all_levs ! Which mixture component we're in

    real( kind = core_rknd ), intent(out), dimension(num_samples) :: &
      lh_sample_point_weights

    ! More Input Variables!
    real( kind = core_rknd ), dimension(d_variables,d_variables,nz), intent(in) :: &
      corr_cholesky_mtx_1, & ! Correlations Cholesky matrix (1st comp.)  [-]
      corr_cholesky_mtx_2    ! Correlations Cholesky matrix (2nd comp.)  [-]

    real( kind = core_rknd ), dimension(d_variables,nz), intent(in) :: &
      mu1,    & ! Means of the hydrometeors, 1st comp. (chi, eta, w, <hydrometeors>)  [units vary]
      mu2,    & ! Means of the hydrometeors, 2nd comp. (chi, eta, w, <hydrometeors>)  [units vary]
      sigma1, & ! Stdevs of the hydrometeors, 1st comp. (chi, eta, w, <hydrometeors>) [units vary]
      sigma2    ! Stdevs of the hydrometeors, 2nd comp. (chi, eta, w, <hydrometeors>) [units vary]

    type(hydromet_pdf_parameter), dimension(nz), intent(in) :: &
      hydromet_pdf_params

    call lh_subcolumn_generator( &
      iter, d_variables, num_samples, sequence_length, nz, & ! In
      pdf_params, delta_zm, rcm, Lscale, & ! In
      rho_ds_zt, mu1, mu2, sigma1, sigma2, & ! In
      corr_cholesky_mtx_1, corr_cholesky_mtx_2, & ! In
      hydromet_pdf_params, & ! In
      X_nl_all_levs, X_mixt_comp_all_levs, & ! Out
      lh_sample_point_weights ) ! Out

  end subroutine lh_subcolumn_generator_api

  !================================================================================================
  ! stats_accumulate_lh - Clips subcolumns from latin hypercube and creates stats.
  !================================================================================================

  subroutine stats_accumulate_lh_api( &
    nz, num_samples, d_variables, rho_ds_zt, &
    lh_sample_point_weights, X_nl_all_levs, &
    lh_clipped_vars )

    use latin_hypercube_driver_module, only : stats_accumulate_lh, lh_clipped_variables_type

    use clubb_precision, only: &
      core_rknd    ! Constant

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      d_variables,     & ! Number of variables to sample
      num_samples,   & ! Number of calls to microphysics per timestep (normally=2)
      nz                 ! Number of vertical model levels

    real( kind = core_rknd ), intent(in), dimension(nz) :: &
      rho_ds_zt  ! Dry, static density (thermo. levs.) [kg/m^3]

    real( kind = core_rknd ), intent(in), dimension(num_samples) :: &
      lh_sample_point_weights

    real( kind = core_rknd ), intent(in), dimension(nz,num_samples,d_variables) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    type(lh_clipped_variables_type), intent(in), dimension(nz,num_samples) :: &
      lh_clipped_vars

    call stats_accumulate_lh( &
      nz, num_samples, d_variables, rho_ds_zt, &
      lh_sample_point_weights, X_nl_all_levs, &
      lh_clipped_vars )

  end subroutine stats_accumulate_lh_api

  !================================================================================================
  ! est_kessler_microphys - Computes microphysical grid box means of Kessler autoconversion scheme.
  !================================================================================================

  subroutine est_kessler_microphys_api( &
    nz, num_samples, d_variables, &
    X_nl_all_levs, pdf_params, rcm, cloud_frac, &
    X_mixt_comp_all_levs, lh_sample_point_weights, &
    lh_AKm, AKm, AKstd, AKstd_cld, &
    AKm_rcm, AKm_rcc, lh_rcm_avg )

    use est_kessler_microphys_module, only : est_kessler_microphys

    use pdf_parameter_module, only:  &
      pdf_parameter  ! Type

    use clubb_precision, only: &
      core_rknd

    implicit none

    ! Input Variables

    integer, intent(in) :: &
      nz, &          ! Number of vertical levels
      num_samples, & ! Number of sample points
      d_variables    ! Number of variates

    real( kind = core_rknd ), dimension(nz,num_samples,d_variables), intent(in) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      cloud_frac    ! Cloud fraction           [-]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rcm          ! Liquid water mixing ratio                [kg/kg]

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params ! PDF parameters       [units vary]

    integer, dimension(nz,num_samples), intent(in) :: &
      X_mixt_comp_all_levs ! Whether we're in mixture component 1 or 2

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      lh_sample_point_weights ! Weight for cloud weighted sampling

    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      lh_AKm,    & ! Monte Carlo estimate of Kessler autoconversion [kg/kg/s]
      AKm,       & ! Exact Kessler autoconversion, AKm,             [kg/kg/s]
      AKstd,     & ! Exact standard deviation of gba Kessler        [kg/kg/s]
      AKstd_cld, & ! Exact w/in cloud std of gba Kessler            [kg/kg/s]
      AKm_rcm,   & ! Exact local gba Kessler auto based on rcm      [kg/kg/s]
      AKm_rcc      ! Exact local gba Kessler based on w/in cloud rc [kg/kg/s]

    ! For comparison, estimate kth liquid water using Monte Carlo
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      lh_rcm_avg ! lh estimate of grid box avg liquid water [kg/kg]

    call est_kessler_microphys( &
      nz, num_samples, d_variables, &
      X_nl_all_levs, pdf_params, rcm, cloud_frac, &
      X_mixt_comp_all_levs, lh_sample_point_weights, &
      lh_AKm, AKm, AKstd, AKstd_cld, &
      AKm_rcm, AKm_rcc, lh_rcm_avg )

  end subroutine est_kessler_microphys_api

  !================================================================================================
  ! clip_transform_silhs_output - Computes extra SILHS sample variables, such as rt and thl.
  !================================================================================================

  subroutine clip_transform_silhs_output_api( nz, num_samples, d_variables, &        ! In
                                              X_mixt_comp_all_levs, X_nl_all_levs, & ! In
                                              pdf_params, l_use_Ncn_to_Nc, &         ! In
                                              lh_clipped_vars )                      ! Out

    use latin_hypercube_driver_module, only : clip_transform_silhs_output, lh_clipped_variables_type

    use clubb_precision, only: &
      core_rknd       ! Our awesome generalized precision (constant)

    use pdf_parameter_module, only: &
      pdf_parameter

    implicit none

    ! Input Variables
    logical, intent(in) :: l_use_Ncn_to_Nc

    integer, intent(in) :: &
      nz,          &         ! Number of vertical levels
      num_samples, &         ! Number of SILHS sample points
      d_variables            ! Number of variates in X_nl_one_lev

    integer, dimension(nz,num_samples), intent(in) :: &
      X_mixt_comp_all_levs   ! Which component this sample is in (1 or 2)

    real( kind = core_rknd ), dimension(nz,num_samples,d_variables) :: &
      X_nl_all_levs          ! A SILHS sample

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params             ! **The** PDF parameters!

    ! Output Variables
    type(lh_clipped_variables_type), dimension(nz,num_samples), intent(out) :: &
      lh_clipped_vars        ! SILHS clipped and transformed variables

    call clip_transform_silhs_output( nz, num_samples, d_variables, &         ! In
                                      X_mixt_comp_all_levs, X_nl_all_levs, &  ! In
                                      pdf_params, l_use_Ncn_to_Nc, &          ! In
                                      lh_clipped_vars )                       ! Out

  end subroutine clip_transform_silhs_output_api

  !-----------------------------------------------------------------
  ! lh_microphys_var_covar_driver: Computes the effect of microphysics on gridbox covariances
  !-----------------------------------------------------------------

  subroutine lh_microphys_var_covar_driver_api &
             ( nz, num_samples, dt, lh_sample_point_weights, &
               lh_rt_all, lh_thl_all, lh_w_all, &
               lh_rcm_mc_all, lh_rvm_mc_all, lh_thlm_mc_all, &
               lh_rtp2_mc_zt, lh_thlp2_mc_zt, lh_wprtp_mc_zt, &
               lh_wpthlp_mc_zt, lh_rtpthlp_mc_zt )

    use lh_microphys_var_covar_module, only: &
      lh_microphys_var_covar_driver  ! Procedure

    use clubb_precision, only: &
      core_rknd   ! Constant

    implicit none

    ! Input Variables!
    integer, intent(in) :: &
      nz,           &                  ! Number of vertical levels
      num_samples                      ! Number of SILHS sample points

    real( kind = core_rknd ), intent(in) :: &
      dt                               ! Model time step                             [s]

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      lh_sample_point_weights          ! Weight of SILHS sample points

    real( kind = core_rknd ), dimension(nz,num_samples), intent(in) :: &
      lh_rt_all, &                     ! SILHS samples of total water                [kg/kg]
      lh_thl_all, &                    ! SILHS samples of potential temperature      [K]
      lh_w_all, &                      ! SILHS samples of vertical velocity          [m/s]
      lh_rcm_mc_all, &                 ! SILHS microphys. tendency of rcm            [kg/kg/s]
      lh_rvm_mc_all, &                 ! SILHS microphys. tendency of rvm            [kg/kg/s]
      lh_thlm_mc_all                   ! SILHS microphys. tendency of thlm           [K/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      lh_rtp2_mc_zt,   &               ! SILHS microphys. est. tendency of <rt'^2>   [(kg/kg)^2/s]
      lh_thlp2_mc_zt,  &               ! SILHS microphys. est. tendency of <thl'^2>  [K^2/s]
      lh_wprtp_mc_zt,  &               ! SILHS microphys. est. tendency of <w'rt'>   [m*(kg/kg)/s^2]
      lh_wpthlp_mc_zt, &               ! SILHS microphys. est. tendency of <w'thl'>  [m*K/s^2]
      lh_rtpthlp_mc_zt                 ! SILHS microphys. est. tendency of <rt'thl'> [K*(kg/kg)/s]


    call lh_microphys_var_covar_driver &
         ( nz, num_samples, dt, lh_sample_point_weights, &
           lh_rt_all, lh_thl_all, lh_w_all, &
           lh_rcm_mc_all, lh_rvm_mc_all, lh_thlm_mc_all, &
           lh_rtp2_mc_zt, lh_thlp2_mc_zt, lh_wprtp_mc_zt, &
           lh_wpthlp_mc_zt, lh_rtpthlp_mc_zt )

  end subroutine lh_microphys_var_covar_driver_api

#endif

end module silhs_api_module
