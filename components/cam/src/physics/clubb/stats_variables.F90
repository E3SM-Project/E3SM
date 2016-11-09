!-------------------------------------------------------------------------------
! $Id: stats_variables.F90 7383 2014-11-13 17:43:38Z schemena@uwm.edu $
!-------------------------------------------------------------------------------

! Description:
!   Holds pointers and other variables for statistics to be written to 
!   GrADS files and netCDF files.
!-------------------------------------------------------------------------------
module stats_variables


  use stats_type, only:  &
      stats ! Type

  use clubb_precision, only:  & 
      core_rknd ! Variable(s)

  implicit none

  private ! Set Default Scope

  ! Sampling and output frequencies
  real( kind = core_rknd ), public :: &
    stats_tsamp = 0._core_rknd, & ! Sampling interval   [s]
    stats_tout  = 0._core_rknd ! Output interval     [s]

!$omp   threadprivate(stats_tsamp, stats_tout)

  logical, public ::  &
    l_stats            = .false., & ! Main flag to turn statistics on/off
    l_output_rad_files = .false., & ! Flag to turn off radiation statistics output
    l_netcdf           = .false., & ! Output to NetCDF format
    l_grads            = .false., & ! Output to GrADS format
    l_silhs_out        = .false., & ! Output SILHS files (stats_lh_zt and stats_lh_sfc)
    l_allow_small_stats_tout = .false. ! Do not stop if output timestep is too low for
                      ! requested format, e.g. l_grads = .true. and
                      ! stats_tout < 60.0

!$omp   threadprivate(l_stats, l_output_rad_files, l_netcdf, l_grads, l_silhs_out, &
!$omp     l_allow_small_stats_tout)

  logical, public :: & 
    l_stats_samp = .false., & ! Sample flag for current time step
    l_stats_last = .false.    ! Last time step of output period

!$omp   threadprivate(l_stats_samp, l_stats_last)

  character(len=200), public ::  & 
    fname_zt     = '', & ! Name of the stats file for thermodynamic grid fields
    fname_lh_zt  = '', & ! Name of the stats file for LH variables on the stats_zt grid
    fname_lh_sfc = '', & ! Name of the stats file for LH variables on the stats_zt grid
    fname_zm     = '', & ! Name of the stats file for momentum grid fields
    fname_rad_zt = '', & ! Name of the stats file for the stats_zt radiation grid fields
    fname_rad_zm = '', & ! Name of the stats file for the stats_zm radiation grid fields
    fname_sfc    = ''    ! Name of the stats file for surface only fields

!$omp   threadprivate(fname_zt, fname_lh_zt, fname_lh_sfc, fname_zm, fname_rad_zt, &
!$omp     fname_rad_zm, fname_sfc)

!       Indices for statistics in stats_zt file

  integer, public :: & 
     ithlm = 0, & 
     ithvm = 0, & 
     irtm = 0, & 
     ircm = 0, &
     irvm = 0, & 
     ium = 0, & 
     ivm = 0, & 
     iwm_zt = 0, &
     iwm_zm = 0, &
     ium_ref = 0,&
     ivm_ref = 0, & 
     iug = 0, & 
     ivg = 0, & 
     icloud_frac = 0, &
     iice_supersat_frac = 0, &
     ircm_in_layer = 0, &
     ircm_in_cloud = 0, &
     icloud_cover = 0, &
     ip_in_Pa = 0, & 
     iexner = 0, & 
     irho_ds_zt = 0, &
     ithv_ds_zt = 0, &
     iLscale = 0, & 
     iwp3 = 0, & 
     iwpthlp2 = 0, & 
     iwp2thlp = 0, & 
     iwprtp2 = 0, & 
     iwp2rtp = 0, &
     iSkw_zt = 0
!$omp threadprivate(ithlm, ithvm, irtm, ircm, irvm, ium, ivm, ium_ref, ivm_ref, &
!$omp   iwm_zt, iwm_zm, iug, ivg, icloud_frac, iice_supersat_frac, ircm_in_layer, &
!$omp   ircm_in_cloud, icloud_cover, &
!$omp   ip_in_Pa, iexner, irho_ds_zt, ithv_ds_zt, iLscale, iwp3, &
!$omp   iwpthlp2, iwp2thlp, iwprtp2, iwp2rtp, iSkw_zt)

  integer, public :: & 
     iLscale_up = 0, & 
     iLscale_down = 0, & 
     iLscale_pert_1 = 0, & 
     iLscale_pert_2 = 0, & 
     itau_zt = 0, & 
     iKh_zt = 0, & 
     iwp2thvp = 0, & 
     iwp2rcp = 0, & 
     iwprtpthlp = 0, & 
     isigma_sqd_w_zt = 0, & 
     irho = 0
!$omp threadprivate( iLscale_up, iLscale_down, &
!$omp   iLscale_pert_1, iLscale_pert_2, &
!$omp   itau_zt, iKh_zt, iwp2thvp, iwp2rcp, iwprtpthlp, isigma_sqd_w_zt, irho )

  integer, dimension(:), allocatable, public :: & 
     icorr_w_hm_ov_adj, &
     ihm1, &
     ihm2
!$omp threadprivate( icorr_w_hm_ov_adj, ihm1, ihm2 )

  integer, public :: & 
     iLWP1 = 0, &
     iLWP2 = 0, &
     iprecip_frac = 0, &
     iprecip_frac_1 = 0, &
     iprecip_frac_2 = 0, &
     iNcnm = 0 
!$omp threadprivate( iLWP1, iLWP2, iprecip_frac, &
!$omp   iprecip_frac_1, iprecip_frac_2, iNcnm )

  integer, dimension(:), allocatable, public :: &
     imu_hm_1,         &
     imu_hm_2,         &
     imu_hm_1_n,       &
     imu_hm_2_n,       &
     isigma_hm_1,      &
     isigma_hm_2,      &
     isigma_hm_1_n,    &
     isigma_hm_2_n,    &
     icorr_w_hm_1,     &
     icorr_w_hm_2,     &
     icorr_chi_hm_1,   &
     icorr_chi_hm_2,   &
     icorr_eta_hm_1,   &
     icorr_eta_hm_2,   &
     icorr_Ncn_hm_1,   &
     icorr_Ncn_hm_2,   &
     icorr_w_hm_1_n,   &
     icorr_w_hm_2_n,   &
     icorr_chi_hm_1_n, &
     icorr_chi_hm_2_n, &
     icorr_eta_hm_1_n, &
     icorr_eta_hm_2_n, &
     icorr_Ncn_hm_1_n, &
     icorr_Ncn_hm_2_n
!$omp threadprivate( imu_hm_1, imu_hm_2, imu_hm_1_n, imu_hm_2_n, &
!$omp   isigma_hm_1, isigma_hm_2, isigma_hm_1_n, isigma_hm_2_n, &
!$omp   icorr_w_hm_1, icorr_w_hm_2, icorr_chi_hm_1, icorr_chi_hm_2, &
!$omp   icorr_eta_hm_1, icorr_eta_hm_2, icorr_Ncn_hm_1, icorr_Ncn_hm_2, &
!$omp   icorr_w_hm_1_n, icorr_w_hm_2_n, icorr_chi_hm_1_n, icorr_chi_hm_2_n, &
!$omp   icorr_eta_hm_1_n, icorr_eta_hm_2_n, icorr_Ncn_hm_1_n, icorr_Ncn_hm_2_n )

  integer, dimension(:,:), allocatable, public :: &
     icorr_hmx_hmy_1,   &
     icorr_hmx_hmy_2,   &
     icorr_hmx_hmy_1_n, &
     icorr_hmx_hmy_2_n
!$omp threadprivate( icorr_hmx_hmy_1, icorr_hmx_hmy_2, &
!$omp   icorr_hmx_hmy_1_n, icorr_hmx_hmy_2_n )

  integer, public :: &
     imu_Ncn_1 = 0,      &
     imu_Ncn_2 = 0,      &
     imu_Ncn_1_n = 0,    &
     imu_Ncn_2_n = 0,    &
     isigma_Ncn_1 = 0,   &
     isigma_Ncn_2 = 0,   &
     isigma_Ncn_1_n = 0, &
     isigma_Ncn_2_n = 0
!$omp threadprivate( imu_Ncn_1, imu_Ncn_2, imu_Ncn_1_n, imu_Ncn_2_n, &
!$omp   isigma_Ncn_1, isigma_Ncn_2, isigma_Ncn_1_n, isigma_Ncn_2_n )

  integer, public :: &
     icorr_w_chi_1 = 0,    &
     icorr_w_chi_2 = 0,    &
     icorr_w_eta_1 = 0,    &
     icorr_w_eta_2 = 0,    &
     icorr_w_Ncn_1 = 0,  &
     icorr_w_Ncn_2 = 0,  &
     icorr_chi_eta_1_ca = 0, &
     icorr_chi_eta_2_ca = 0, &
     icorr_chi_Ncn_1 = 0,  &
     icorr_chi_Ncn_2 = 0,  &
     icorr_eta_Ncn_1 = 0,  &
     icorr_eta_Ncn_2 = 0
!$omp threadprivate( icorr_w_chi_1, icorr_w_chi_2, icorr_w_eta_1, &
!$omp   icorr_w_eta_2, icorr_w_Ncn_1, icorr_w_Ncn_2, icorr_chi_eta_1_ca, &
!$omp   icorr_chi_eta_2_ca, icorr_chi_Ncn_1, icorr_chi_Ncn_2, icorr_eta_Ncn_1, &
!$omp   icorr_eta_Ncn_2 )

  integer, public :: &
     icorr_w_Ncn_1_n = 0, &
     icorr_w_Ncn_2_n = 0, &
     icorr_chi_Ncn_1_n = 0, &
     icorr_chi_Ncn_2_n = 0, &
     icorr_eta_Ncn_1_n = 0, &
     icorr_eta_Ncn_2_n = 0
!$omp threadprivate( icorr_w_Ncn_1_n, icorr_w_Ncn_2_n, icorr_chi_Ncn_1_n, &
!$omp   icorr_chi_Ncn_2_n, icorr_eta_Ncn_1_n, icorr_eta_Ncn_2_n )

  integer, public :: & 
     iNcm = 0,             & ! Brian
     iNccnm = 0,           & 
     iNc_in_cloud = 0,     &
     iNc_activated = 0,    &
     isnowslope = 0,       & ! Adam Smith, 22 April 2008
     ised_rcm = 0,         & ! Brian
     irsat = 0,            & ! Brian
     irsati = 0,           & 
     irrm = 0,          & ! Brian
     im_vol_rad_rain = 0,  & ! Brian
     im_vol_rad_cloud = 0, & ! COAMPS only. dschanen 6 Dec 2006
     iprecip_rate_zt = 0,    & ! Brian
     iAKm = 0,             & ! analytic Kessler.  Vince Larson 22 May 2005 
     ilh_AKm = 0,          & ! LH Kessler.  Vince Larson  22 May 2005
     iradht = 0,           & ! Radiative heating.
     iradht_LW = 0,        & !   "           "   Long-wave component
     iradht_SW = 0,        & !   "           "   Short-wave component
     irel_humidity = 0
!$omp  threadprivate( iNcm, iNccnm, iNc_in_cloud, iNc_activated, isnowslope, &
!$omp    ised_rcm, irsat, irsati, irrm, &
!$omp    im_vol_rad_rain, im_vol_rad_cloud, &
!$omp    iprecip_rate_zt, iAKm, ilh_AKm, &
!$omp    iradht, iradht_LW, iradht_SW, &
!$omp    irel_humidity )

  integer, public :: & 
     iAKstd = 0,     &
     iAKstd_cld = 0, & 
     iAKm_rcm = 0, & 
     iAKm_rcc = 0
!$omp threadprivate( iAKstd, iAKstd_cld, iAKm_rcm, iAKm_rcc )


  integer, public :: & 
   irfrzm = 0
!$omp threadprivate(irfrzm)

  ! Skewness functions on stats_zt grid
  integer, public :: &
    iC11_Skw_fnc = 0

!$omp threadprivate(iC11_Skw_fnc)

  integer, public :: &
    icloud_frac_zm = 0, &
    iice_supersat_frac_zm = 0, &
    ircm_zm = 0, &
    irtm_zm = 0, &
    ithlm_zm = 0

!$omp threadprivate(icloud_frac_zm, iice_supersat_frac_zm, ircm_zm, irtm_zm, ithlm_zm)

  integer, public :: &
    ilh_rcm_avg = 0, &
    ik_lh_start = 0

!$omp threadprivate(ilh_rcm_avg, ik_lh_start)

  integer, public :: & 
     iNrm = 0,       & ! Rain droplet number concentration
     iNim = 0,       & ! Ice number concentration
     iNsm = 0,    & ! Snow number concentration
     iNgm = 0    ! Graupel number concentration
!$omp   threadprivate(iNrm, iNim, iNsm, iNgm)

  integer, public :: & 
     iT_in_K      ! Absolute temperature
!$omp   threadprivate(iT_in_K)

  integer, public :: &
    ieff_rad_cloud = 0, &
    ieff_rad_ice = 0, &
    ieff_rad_snow = 0, &
    ieff_rad_rain = 0, &
    ieff_rad_graupel = 0

!$omp   threadprivate(ieff_rad_cloud, ieff_rad_ice, ieff_rad_snow) 
!$omp   threadprivate(ieff_rad_rain, ieff_rad_graupel)

  integer, public :: & 
    irsm = 0, &
    irgm = 0, & 
    irim = 0, & 
    idiam = 0,           & ! Diameter of ice crystal           [m]
    imass_ice_cryst = 0, & ! Mass of a single ice crystal      [kg]
    ircm_icedfs = 0,     & ! Change in liquid water due to ice [kg/kg/s]
    iu_T_cm = 0         ! Fallspeed of ice crystal in cm/s  [cm s^{-1}]

!$omp threadprivate(irsm, irgm, irim, idiam, &
!$omp   imass_ice_cryst, ircm_icedfs, iu_T_cm)


  ! thlm/rtm budget terms
  integer, public :: & 
    irtm_bt = 0,         & ! rtm total time tendency
    irtm_ma = 0,         & ! rtm mean advect. term
    irtm_ta = 0,         & ! rtm turb. advect. term
    irtm_forcing = 0,    & ! rtm large scale forcing term
    irtm_mc = 0,         & ! rtm change from microphysics
    irtm_sdmp = 0,       & ! rtm change from sponge damping
    irvm_mc = 0,         & ! rvm change from microphysics
    ircm_mc = 0,         & ! rcm change from microphysics
    ircm_sd_mg_morr = 0, & ! rcm sedimentation tendency
    irtm_mfl = 0,        & ! rtm change due to monotonic flux limiter
    irtm_tacl = 0,       & ! rtm correction from turbulent advection (wprtp) clipping
    irtm_cl = 0,         & ! rtm clipping term
    irtm_pd = 0,         & ! thlm postive definite adj term
    ithlm_bt = 0,        & ! thlm total time tendency
    ithlm_ma = 0,        & ! thlm mean advect. term
    ithlm_ta = 0,        & ! thlm turb. advect. term
    ithlm_forcing = 0,   & ! thlm large scale forcing term
    ithlm_sdmp = 0,      & ! thlm change from sponge damping
    ithlm_mc = 0,        & ! thlm change from microphysics
    ithlm_mfl = 0,       & ! thlm change due to monotonic flux limiter
    ithlm_tacl = 0,      & ! thlm correction from turbulent advection (wpthlp) clipping
    ithlm_cl = 0           ! thlm clipping term

!$omp   threadprivate(irtm_bt, irtm_ma, irtm_ta, irtm_forcing, &
!$omp     irtm_mc, irtm_sdmp, irtm_mfl, irtm_tacl, irtm_cl, irtm_pd, &
!$omp     irvm_mc, ircm_mc, ircm_sd_mg_morr, &
!$omp     ithlm_bt, ithlm_ma, ithlm_ta, ithlm_forcing, &
!$omp     ithlm_mc, ithlm_sdmp, ithlm_mfl, ithlm_tacl, ithlm_cl)

  !monatonic flux limiter diagnostic terms
  integer, public :: &
    ithlm_mfl_min = 0, &
    ithlm_mfl_max = 0, &
    iwpthlp_entermfl = 0, &
    iwpthlp_exit_mfl = 0, &
    iwpthlp_mfl_min = 0, &
    iwpthlp_mfl_max = 0, &
    irtm_mfl_min = 0, &
    irtm_mfl_max = 0, &
    iwprtp_enter_mfl = 0, &
    iwprtp_exit_mfl = 0, &
    iwprtp_mfl_min = 0, &
    iwprtp_mfl_max = 0, &
    ithlm_enter_mfl = 0, &
    ithlm_exit_mfl = 0, &
    ithlm_old = 0, &
    ithlm_without_ta = 0, &
    irtm_enter_mfl = 0, &
    irtm_exit_mfl = 0, &
    irtm_old = 0, &
    irtm_without_ta = 0

!$omp   threadprivate(ithlm_mfl_min, ithlm_mfl_max, iwpthlp_entermfl)
!$omp   threadprivate(iwpthlp_exit_mfl, iwpthlp_mfl_min, iwpthlp_mfl_max)
!$omp   threadprivate(irtm_mfl_min, irtm_mfl_max, iwprtp_enter_mfl)
!$omp   threadprivate(iwprtp_exit_mfl, iwprtp_mfl_min, iwprtp_mfl_max)
!$omp   threadprivate(ithlm_enter_mfl, ithlm_exit_mfl, ithlm_old, ithlm_without_ta)
!$omp   threadprivate(irtm_enter_mfl, irtm_exit_mfl, irtm_old, irtm_without_ta)

  integer, public :: & 
     iwp3_bt  = 0, & 
     iwp3_ma  = 0, & 
     iwp3_ta  = 0, & 
     iwp3_tp  = 0, & 
     iwp3_ac  = 0, & 
     iwp3_bp1 = 0, & 
     iwp3_bp2 = 0, & 
     iwp3_pr1 = 0, & 
     iwp3_pr2 = 0, & 
     iwp3_dp1 = 0, &
     iwp3_cl  = 0

!$omp   threadprivate(iwp3_bt, iwp3_ma, iwp3_ta, iwp3_tp, iwp3_ac, iwp3_bp1)
!$omp   threadprivate(iwp3_bp2, iwp3_pr1, iwp3_pr2, iwp3_dp1, iwp3_cl)

  ! Rain mixing ratio budgets
  integer, public :: & 
     irrm_bt = 0, &
     irrm_ma = 0, &
     irrm_ta = 0, &
     irrm_sd = 0, &
     irrm_ts = 0, &
     irrm_sd_morr = 0, &
     irrm_cond = 0, &
     irrm_auto = 0, &
     irrm_accr = 0, &
     irrm_cond_adj = 0, &
     irrm_src_adj = 0, &
     irrm_mc = 0, &
     irrm_hf = 0, &
     irrm_wvhf = 0, &
     irrm_cl = 0

!$omp   threadprivate(irrm_bt, irrm_ma, irrm_ta, irrm_sd)
!$omp   threadprivate(irrm_ts, irrm_sd_morr)
!$omp   threadprivate(irrm_cond, irrm_auto, irrm_accr)
!$omp   threadprivate(irrm_cond_adj, irrm_src_adj )
!$omp   threadprivate(irrm_mc, irrm_hf, irrm_wvhf, irrm_cl)

  integer, public :: &
     iNrm_bt = 0, &
     iNrm_ma = 0, &
     iNrm_ta = 0, &
     iNrm_sd = 0, &
     iNrm_ts = 0, &
     iNrm_cond = 0, &
     iNrm_auto = 0, &
     iNrm_cond_adj = 0, &
     iNrm_src_adj = 0, &
     iNrm_mc = 0, &
     iNrm_cl = 0

!$omp   threadprivate(iNrm_bt, iNrm_ma, iNrm_ta, iNrm_sd, iNrm_ts, iNrm_cond)
!$omp   threadprivate(iNrm_auto, iNrm_cond_adj, iNrm_src_adj )
!$omp   threadprivate(iNrm_mc, iNrm_cl)


  ! Snow/Ice/Graupel mixing ratio budgets
  integer, public :: & 
     irsm_bt = 0, & 
     irsm_ma = 0, & 
     irsm_sd = 0, & 
     irsm_sd_morr = 0, &
     irsm_ta = 0, & 
     irsm_mc = 0, & 
     irsm_hf = 0, &
     irsm_wvhf = 0, &
     irsm_cl = 0, &
     irsm_sd_morr_int = 0

!$omp   threadprivate(irsm_bt, irsm_ma, irsm_sd, irsm_sd_morr, irsm_ta)
!$omp   threadprivate(irsm_mc, irsm_hf, irsm_wvhf, irsm_cl)
!$omp   threadprivate(irsm_sd_morr_int)

  integer, public :: & 
     irgm_bt = 0, & 
     irgm_ma = 0, & 
     irgm_sd = 0, & 
     irgm_sd_morr = 0, &
     irgm_ta = 0, & 
     irgm_mc = 0, & 
     irgm_hf = 0, &
     irgm_wvhf = 0, &
     irgm_cl = 0

!$omp   threadprivate(irgm_bt, irgm_ma, irgm_sd, irgm_sd_morr)
!$omp   threadprivate(irgm_ta, irgm_mc)
!$omp   threadprivate(irgm_hf, irgm_wvhf, irgm_cl)

  integer, public :: & 
     irim_bt = 0, & 
     irim_ma = 0, & 
     irim_sd = 0, & 
     irim_sd_mg_morr = 0, &
     irim_ta = 0, & 
     irim_mc = 0, & 
     irim_hf = 0, &
     irim_wvhf = 0, &
     irim_cl = 0

!$omp   threadprivate(irim_bt, irim_ma, irim_sd, irim_sd_mg_morr, irim_ta)
!$omp   threadprivate(irim_mc, irim_hf, irim_wvhf, irim_cl)

  integer, public :: &
    iNsm_bt = 0,  &
    iNsm_ma = 0,  &
    iNsm_sd = 0,  &
    iNsm_ta = 0,  &
    iNsm_mc = 0,  &
    iNsm_cl = 0

!$omp threadprivate(iNsm_bt, iNsm_ma, iNsm_sd, iNsm_ta, &
!$omp   iNsm_mc, iNsm_cl)

  integer, public :: &
    iNgm_bt = 0, &
    iNgm_ma = 0, &
    iNgm_sd = 0, &
    iNgm_ta = 0, &
    iNgm_mc = 0, &
    iNgm_cl = 0

!$omp threadprivate(iNgm_bt, iNgm_ma, iNgm_sd, &
!$omp   iNgm_ta, iNgm_mc, iNgm_cl)

  integer, public :: &
    iNim_bt = 0, &
    iNim_ma = 0, &
    iNim_sd = 0, &
    iNim_ta = 0, &
    iNim_mc = 0, &
    iNim_cl = 0

!$omp threadprivate(iNim_bt, iNim_ma, iNim_sd, iNim_ta, &
!$omp   iNim_mc, iNim_cl)

  integer, public :: &
    iNcm_bt = 0, &
    iNcm_ma = 0, &
    iNcm_ta = 0, &
    iNcm_mc = 0, &
    iNcm_cl = 0, &
    iNcm_act = 0

!$omp threadprivate(iNcm_bt, iNcm_ma, iNcm_ta, &
!$omp   iNcm_mc, iNcm_cl, iNcm_act)

  ! Covariances between w, r_t, theta_l and KK microphysics tendencies.
  ! Additionally, covariances between r_r and N_r and KK rain drop mean
  ! volume radius.  These are all calculated on thermodynamic grid levels.
  integer, public :: &
    iw_KK_evap_covar_zt = 0,   & ! Covariance of w and KK evaporation tendency.
    irt_KK_evap_covar_zt = 0,  & ! Covariance of r_t and KK evaporation tendency.
    ithl_KK_evap_covar_zt = 0, & ! Covariance of theta_l and KK evap. tendency.
    iw_KK_auto_covar_zt = 0,   & ! Covariance of w and KK autoconversion tendency.
    irt_KK_auto_covar_zt = 0,  & ! Covariance of r_t and KK autoconversion tendency.
    ithl_KK_auto_covar_zt = 0, & ! Covariance of theta_l and KK autoconv. tendency.
    iw_KK_accr_covar_zt = 0,   & ! Covariance of w and KK accretion tendency.
    irt_KK_accr_covar_zt = 0,  & ! Covariance of r_t and KK accretion tendency.
    ithl_KK_accr_covar_zt = 0, & ! Covariance of theta_l and KK accretion tendency.
    irr_KK_mvr_covar_zt = 0,   & ! Covariance of r_r and KK mean volume radius.
    iNr_KK_mvr_covar_zt = 0,   & ! Covariance of N_r and KK mean volume radius.
    iKK_mvr_variance_zt = 0      ! Variance of KK rain drop mean volume radius.

!$omp threadprivate( iw_KK_evap_covar_zt, irt_KK_evap_covar_zt, &
!$omp   ithl_KK_evap_covar_zt, iw_KK_auto_covar_zt, irt_KK_auto_covar_zt, &
!$omp   ithl_KK_auto_covar_zt, iw_KK_accr_covar_zt, irt_KK_accr_covar_zt, &
!$omp   ithl_KK_accr_covar_zt, irr_KK_mvr_covar_zt, iNr_KK_mvr_covar_zt, &
!$omp   iKK_mvr_variance_zt )

  ! Wind budgets
  integer, public :: & 
     ivm_bt = 0, & 
     ivm_ma = 0, & 
     ivm_ta = 0, & 
     ivm_gf = 0, & 
     ivm_cf = 0, &
     ivm_f = 0, &
     ivm_sdmp = 0, &
     ivm_ndg = 0

!$omp   threadprivate(ivm_bt, ivm_ma, ivm_ta, ivm_gf, ivm_cf, ivm_f, ivm_sdmp, ivm_ndg)

  integer, public :: & 
     ium_bt = 0, & 
     ium_ma = 0, & 
     ium_ta = 0, & 
     ium_gf = 0, & 
     ium_cf = 0, & 
     ium_f = 0, &
     ium_sdmp = 0, &
     ium_ndg = 0

!$omp   threadprivate(ium_bt, ium_ma, ium_ta, ium_gf, ium_cf, ium_f, ium_sdmp, ium_ndg)


  ! PDF parameters
  integer, public :: & 
     imixt_frac = 0, & 
     iw_1 = 0, & 
     iw_2 = 0, & 
     ivarnce_w_1 = 0, & 
     ivarnce_w_2 = 0, & 
     ithl_1 = 0, & 
     ithl_2 = 0, & 
     ivarnce_thl_1 = 0, & 
     ivarnce_thl_2 = 0, & 
     irt_1 = 0, & 
     irt_2 = 0, & 
     ivarnce_rt_1 = 0, & 
     ivarnce_rt_2 = 0, & 
     irc_1 = 0, & 
     irc_2 = 0, & 
     irsatl_1 = 0, & 
     irsatl_2 = 0, & 
     icloud_frac_1 = 0, & 
     icloud_frac_2 = 0
!$omp  threadprivate(imixt_frac, iw_1, iw_2, ivarnce_w_1, ivarnce_w_2, ithl_1, ithl_2, &
!$omp  ivarnce_thl_1, ivarnce_thl_2, irt_1, irt_2, ivarnce_rt_1, ivarnce_rt_2, irc_1, irc_2, &
!$omp  irsatl_1, irsatl_2, icloud_frac_1, icloud_frac_2 )

  integer, public :: & 
     ichi_1 = 0, &
     ichi_2 = 0, &
     istdev_chi_1 = 0, & 
     istdev_chi_2 = 0, &
     ichip2 = 0, &
     istdev_eta_1 = 0, &
     istdev_eta_2 = 0, &
     icovar_chi_eta_1 = 0, &
     icovar_chi_eta_2 = 0, &
     icorr_chi_eta_1 = 0, &
     icorr_chi_eta_2 = 0, &
     irrtthl = 0, &
     icrt_1 = 0, &
     icrt_2 = 0, &
     icthl_1 = 0, &
     icthl_2 = 0
!$omp  threadprivate( ichi_1, ichi_2, istdev_chi_1, istdev_chi_2, ichip2, &
!$omp    istdev_eta_1, istdev_eta_2, icovar_chi_eta_1, icovar_chi_eta_2, &
!$omp    icorr_chi_eta_1, icorr_chi_eta_2, irrtthl, icrt_1, icrt_2, icthl_1, &
!$omp    icthl_2 )

  integer, public :: & 
    iwp2_zt = 0, & 
    ithlp2_zt = 0, & 
    iwpthlp_zt = 0, & 
    iwprtp_zt = 0, & 
    irtp2_zt = 0, & 
    irtpthlp_zt = 0, &
    iup2_zt = 0, &
    ivp2_zt = 0, &
    iupwp_zt = 0, &
    ivpwp_zt = 0

!$omp   threadprivate( iwp2_zt, ithlp2_zt, iwpthlp_zt, iwprtp_zt, irtp2_zt, &
!$omp                  irtpthlp_zt, iup2_zt, ivp2_zt, iupwp_zt, ivpwp_zt )

  integer, dimension(:), allocatable, public :: &
    iwp2hmp

!$omp   threadprivate( iwp2hmp )

  integer, dimension(:), allocatable, public :: &
    ihydrometp2, &
    iwphydrometp, &
    irtphmp,      &
    ithlphmp

!$omp   threadprivate( ihydrometp2, iwphydrometp, irtphmp, ithlphmp )

  integer, dimension(:), allocatable, public :: &
    ihmp2_zt

!$omp  threadprivate( ihmp2_zt )

  integer, public :: &
    ichi = 0
!$omp threadprivate(ichi)

  integer, target, allocatable, dimension(:), public :: & 
    isclrm,   & ! Passive scalar mean (1)
    isclrm_f    ! Passive scalar forcing (1)
!$omp   threadprivate(isclrm, isclrm_f)

! Used to calculate clear-sky radiative fluxes.
  integer, public :: &
    ifulwcl = 0, ifdlwcl = 0, ifdswcl = 0, ifuswcl = 0

!$omp   threadprivate(ifulwcl, ifdlwcl, ifdswcl, ifuswcl)

  integer, target, allocatable, dimension(:), public :: & 
    iedsclrm,   & ! Eddy-diff. scalar term (1)
    iedsclrm_f    ! Eddy-diffusivity scalar forcing (1)

!$omp   threadprivate(iedsclrm, iedsclrm_f)

  integer, public :: &
    ilh_thlm_mc = 0,      & ! Latin hypercube estimate of thlm_mc
    ilh_rvm_mc = 0,       & ! Latin hypercube estimate of rvm_mc
    ilh_rcm_mc = 0,       & ! Latin hypercube estimate of rcm_mc
    ilh_Ncm_mc = 0,       & ! Latin hypercube estimate of Ncm_mc
    ilh_rrm_mc = 0,    & ! Latin hypercube estimate of rrm_mc
    ilh_Nrm_mc = 0,       & ! Latin hypercube estimate of Nrm_mc
    ilh_rsm_mc = 0,    & ! Latin hypercube estimate of rsm_mc
    ilh_Nsm_mc = 0,    & ! Latin hypercube estimate of Nsm_mc
    ilh_rgm_mc = 0, & ! Latin hypercube estimate of rgm_mc
    ilh_Ngm_mc = 0, & ! Latin hypercube estimate of Ngm_mc
    ilh_rim_mc = 0,     & ! Latin hypercube estimate of rim_mc
    ilh_Nim_mc = 0          ! Latin hypercube estimate of Nim_mc
!$omp   threadprivate( ilh_thlm_mc, ilh_rvm_mc, ilh_rcm_mc, ilh_Ncm_mc, &
!$omp     ilh_rrm_mc,  ilh_Nrm_mc, ilh_rsm_mc, ilh_Nsm_mc, &
!$omp     ilh_rgm_mc, ilh_Ngm_mc, ilh_rim_mc, ilh_Nim_mc )

  integer, public :: &
    ilh_rrm_auto = 0, & ! Latin hypercube estimate of autoconversion
    ilh_rrm_accr = 0, & ! Latin hypercube estimate of accretion
    ilh_rrm_evap = 0, & ! Latin hypercube estimate of evaporation
    ilh_Nrm_auto    = 0, & ! Latin hypercube estimate of Nrm autoconversion
    ilh_Nrm_cond    = 0, & ! Latin hypercube estimate of Nrm evaporation
    ilh_m_vol_rad_rain = 0

!$omp   threadprivate( ilh_rrm_auto, ilh_rrm_accr, ilh_rrm_evap, &
!$omp                  ilh_Nrm_auto, ilh_Nrm_cond, & 
!$omp                  ilh_m_vol_rad_rain )

  integer, public :: &
    ilh_rrm_src_adj  = 0, & ! Latin hypercube estimate of source adjustment (KK only!)
    ilh_rrm_cond_adj = 0, & ! Latin hypercube estimate of evap adjustment (KK only!)
    ilh_Nrm_src_adj     = 0, & ! Latin hypercube estimate of Nrm source adjustmet (KK only!)
    ilh_Nrm_cond_adj    = 0    ! Latin hypercube estimate of Nrm evap adjustment (KK only!)
!$omp   threadprivate( ilh_rrm_src_adj, ilh_rrm_cond_adj, ilh_Nrm_src_adj, &
!$omp                  ilh_Nrm_cond_adj     )

  integer, public :: &
    ilh_Vrr = 0, & ! Latin hypercube estimate of rrm sedimentation velocity
    ilh_VNr = 0    ! Latin hypercube estimate of Nrm sedimentation velocity
!$omp   threadprivate(ilh_Vrr,  ilh_VNr)

  integer, public :: &
    ilh_rrm = 0, &
    ilh_Nrm = 0, &
    ilh_rim = 0, &
    ilh_Nim = 0, &
    ilh_rsm = 0, &
    ilh_Nsm = 0, &
    ilh_rgm = 0, &
    ilh_Ngm = 0, &
    ilh_thlm = 0, &
    ilh_rcm = 0, &
    ilh_Ncm = 0, &
    ilh_Ncnm = 0, &
    ilh_rvm = 0, &
    ilh_wm = 0, &
    ilh_cloud_frac = 0, &
    ilh_chi = 0, &
    ilh_eta = 0, &
    ilh_precip_frac = 0, &
    ilh_mixt_frac = 0

!$omp threadprivate(ilh_rrm, ilh_Nrm, ilh_rim, ilh_Nim, ilh_rsm, ilh_Nsm, &
!$omp   ilh_rgm, ilh_Ngm, &
!$omp   ilh_thlm, ilh_rcm, ilh_Ncm, ilh_Ncnm, ilh_rvm, ilh_wm, ilh_cloud_frac, &
!$omp   ilh_chi, ilh_eta, ilh_precip_frac, ilh_mixt_frac )

  integer, public :: &
    ilh_cloud_frac_unweighted  = 0,  &
    ilh_precip_frac_unweighted = 0,  &
    ilh_mixt_frac_unweighted   = 0

!$omp threadprivate( ilh_cloud_frac_unweighted, ilh_precip_frac_unweighted, &
!$omp                ilh_mixt_frac_unweighted )

  integer, public :: &
    ilh_wp2_zt = 0, &
    ilh_Nrp2_zt = 0, &
    ilh_Ncnp2_zt = 0, &
    ilh_Ncp2_zt = 0, &
    ilh_rcp2_zt = 0, &
    ilh_rtp2_zt = 0, &
    ilh_thlp2_zt = 0, &
    ilh_rrp2_zt = 0, &
    ilh_chip2 = 0 ! Eric Raut
!$omp threadprivate( ilh_wp2_zt, ilh_Nrp2_zt, ilh_Ncnp2_zt, ilh_Ncp2_zt, &
!$omp                ilh_rcp2_zt, ilh_rtp2_zt, ilh_thlp2_zt, ilh_rrp2_zt, ilh_chip2 )


  ! Indices for Morrison budgets
  integer, public :: &
    iPSMLT = 0, &
    iEVPMS = 0, &
    iPRACS = 0, &
    iEVPMG = 0, &
    iPRACG = 0, &
    iPGMLT = 0, &
    iMNUCCC = 0, &
    iPSACWS = 0, &
    iPSACWI = 0, &
    iQMULTS = 0, &
    iQMULTG = 0, &
    iPSACWG = 0, &
    iPGSACW = 0, &
    iPRD = 0, &
    iPRCI = 0, &
    iPRAI = 0, &
    iQMULTR = 0, &
    iQMULTRG = 0, &
    iMNUCCD = 0, &
    iPRACI = 0, &
    iPRACIS = 0, &
    iEPRD = 0, &
    iMNUCCR = 0, &
    iPIACR = 0, &
    iPIACRS = 0, &
    iPGRACS = 0, &
    iPRDS = 0, &
    iEPRDS = 0, &
    iPSACR = 0, &
    iPRDG = 0, &
    iEPRDG = 0

!$omp threadprivate( iPSMLT, iEVPMS, iPRACS, iEVPMG, iPRACG, iPGMLT, iMNUCCC, iPSACWS, iPSACWI, &
!$omp   iQMULTS, iQMULTG, iPSACWG, iPGSACW, iPRD, iPRCI, iPRAI, iQMULTR, &
!$omp   iQMULTRG, iMNUCCD, iPRACI, iPRACIS, iEPRD, iMNUCCR, iPIACR, iPIACRS, &
!$omp   iPGRACS, iPRDS, iEPRDS, iPSACR, iPRDG, iEPRDG  )

  ! More indices for Morrison budgets!!
  integer, public :: &
    iNGSTEN = 0, &
    iNRSTEN = 0, &
    iNISTEN = 0, &
    iNSSTEN = 0, &
    iNCSTEN = 0, &
    iNPRC1 = 0,  &
    iNRAGG = 0,  &
    iNPRACG = 0, &
    iNSUBR = 0,  &
    iNSMLTR = 0, &
    iNGMLTR = 0, &
    iNPRACS = 0, &
    iNNUCCR = 0, &
    iNIACR = 0,  &
    iNIACRS = 0, &
    iNGRACS = 0, &
    iNSMLTS = 0, &
    iNSAGG = 0,  &
    iNPRCI = 0, &
    iNSCNG = 0, &
    iNSUBS = 0, &
    iPRC = 0, &
    iPRA = 0, &
    iPRE = 0

!$omp threadprivate( iNGSTEN, iNRSTEN, iNISTEN, iNSSTEN, iNCSTEN, iNPRC1, iNRAGG, &
!$omp   iNPRACG, iNSUBR,  iNSMLTR, iNGMLTR, iNPRACS, iNNUCCR, iNIACR, &
!$omp   iNIACRS, iNGRACS, iNSMLTS, iNSAGG, iNPRCI, iNSCNG, iNSUBS, iPRC, iPRA, iPRE )

  ! More indices for Morrison budgets!!
  integer, public :: &
    iPCC = 0, & 
    iNNUCCC = 0, & 
    iNPSACWS = 0, &
    iNPRA = 0, & 
    iNPRC = 0, & 
    iNPSACWI = 0, &
    iNPSACWG = 0, &
    iNPRAI = 0, &
    iNMULTS = 0, & 
    iNMULTG = 0, & 
    iNMULTR = 0, & 
    iNMULTRG = 0, & 
    iNNUCCD = 0, & 
    iNSUBI = 0, & 
    iNGMLTG = 0, &
    iNSUBG = 0, &
    iNACT = 0

  integer, public :: &
    iSIZEFIX_NR = 0, &
    iSIZEFIX_NC = 0, &
    iSIZEFIX_NI = 0, &
    iSIZEFIX_NS = 0, &
    iSIZEFIX_NG = 0, &
    iNEGFIX_NR = 0, &
    iNEGFIX_NC = 0, &
    iNEGFIX_NI = 0, &
    iNEGFIX_NS = 0, &
    iNEGFIX_NG = 0, &
    iNIM_MORR_CL = 0, &
    iQC_INST = 0, &
    iQR_INST = 0, &
    iQI_INST = 0, &
    iQS_INST = 0, &
    iQG_INST = 0, &
    iNC_INST = 0, &
    iNR_INST = 0, &
    iNI_INST = 0, &
    iNS_INST = 0, &
    iNG_INST = 0, &
    iT_in_K_mc = 0, &
    ihl_on_Cp_residual = 0, &
    iqto_residual = 0

!$omp threadprivate(iPCC, iNNUCCC, iNPSACWS, iNPRA, iNPRC, iNPSACWI, iNPSACWG, iNPRAI, &
!$omp   iNMULTS, iNMULTG, iNMULTR, iNMULTRG, iNNUCCD, iNSUBI, iNGMLTG, iNSUBG, iNACT, &
!$omp   iSIZEFIX_NR, iSIZEFIX_NC, iSIZEFIX_NI, iSIZEFIX_NS, iSIZEFIX_NG, iNEGFIX_NR, &
!$omp   iNEGFIX_NC, iNEGFIX_NI, iNEGFIX_NS, iNEGFIX_NG, iNIM_MORR_CL, iQC_INST, iQR_INST, &
!$omp   iQI_INST, iQS_INST, iQG_INST, iNC_INST, iNR_INST, iNI_INST, iNS_INST, &
!$omp   iNG_INST, iT_in_K_mc, ihl_on_Cp_residual, iqto_residual  )

  ! Indices for statistics in stats_zm file
  integer, public :: & 
     iwp2 = 0, & 
     irtp2 = 0, & 
     ithlp2 = 0, & 
     irtpthlp = 0, & 
     iwprtp = 0, & 
     iwpthlp = 0, & 
     iwp4 = 0, & 
     iwpthvp = 0, & 
     irtpthvp = 0, & 
     ithlpthvp = 0, & 
     itau_zm = 0, & 
     iKh_zm = 0, & 
     iwprcp = 0, & 
     irc_coef = 0, &
     ithlprcp = 0, & 
     irtprcp = 0, & 
     ircp2 = 0, & 
     iupwp = 0, & 
     ivpwp = 0, &
     iSkw_zm = 0
!$omp threadprivate(iSkw_zm)

  integer, public :: &
     irho_zm = 0, & 
     isigma_sqd_w = 0, &
     irho_ds_zm = 0, &
     ithv_ds_zm = 0, &
     iem = 0, & 
     ishear = 0, & ! Brian
     imean_w_up = 0, &
     imean_w_down = 0, &
     iFrad = 0, & 
     iFrad_LW = 0,   & ! Brian
     iFrad_SW = 0,   & ! Brian
     iFrad_LW_up = 0,   & 
     iFrad_SW_up = 0,   & 
     iFrad_LW_down = 0,   & 
     iFrad_SW_down = 0,   & 
     iFprec = 0,          & ! Brian
     iFcsed = 0             ! Brian

   ! Stability correction applied to Kh_N2_zm (diffusion on rtm and thlm)
   integer, public :: &
     istability_correction = 0 ! schemena

!$omp   threadprivate(istability_correction)
!$omp   threadprivate(iwp2, irtp2, ithlp2, irtpthlp, iwprtp, iwpthlp)
!$omp   threadprivate(iwp4, iwpthvp, irtpthvp, ithlpthvp, itau_zm, iKh_zm)
!$omp   threadprivate(iwprcp, irc_coef, ithlprcp, irtprcp, ircp2, iupwp, ivpwp)
!$omp   threadprivate(irho_zm, isigma_sqd_w, irho_ds_zm, ithv_ds_zm, iem, ishear)
!$omp   threadprivate(imean_w_up, imean_w_down)
!$omp   threadprivate(iFrad, iFrad_LW, iFrad_SW, iFrad_SW_up, iFrad_SW_down)
!$omp   threadprivate(iFrad_LW_up, iFrad_LW_down, iFprec, iFcsed)

  integer, dimension(:), allocatable, public :: &
    iK_hm
!$omp   threadprivate(iK_hm)

  ! Skewness Functions on stats_zm grid
  integer, public :: &
    igamma_Skw_fnc = 0,  &
    iC6rt_Skw_fnc = 0,   &
    iC6thl_Skw_fnc = 0,  &
    iC7_Skw_fnc = 0,     &
    iC1_Skw_fnc = 0

!$omp   threadprivate(igamma_Skw_fnc, iC6rt_Skw_fnc, iC6thl_Skw_fnc)
!$omp   threadprivate(iC7_Skw_fnc, iC1_Skw_fnc)

  ! Covariance of w and cloud droplet concentration, < w'N_c' >
  integer, public :: &
    iwpNcp = 0

!$omp   threadprivate( iwpNcp )

  ! Sedimentation velocities
  integer, public :: & 
    iVNr = 0,    &
    iVrr = 0,    &
    iVNc = 0,    &
    iVrc = 0,    &
    iVNs = 0, &
    iVrs = 0, &
    iVNi = 0,  &
    iVri = 0,  &
    iVrg = 0

!$omp   threadprivate(iVNr, iVrr, iVNc, iVrc, iVNs, iVrs, iVNi, iVri, iVrg)

  ! Covariance of sedimentation velocity and hydrometeor, <V_xx'x_x'>.
  integer, public :: &
    iVrrprrp = 0,         &
    iVNrpNrp = 0,         &
    iVrrprrp_expcalc = 0, &
    iVNrpNrp_expcalc = 0

!$omp   threadprivate(iVrrprrp, iVNrpNrp, iVrrprrp_expcalc, iVNrpNrp_expcalc)

  integer, public :: & 
     iwp2_bt = 0, & 
     iwp2_ma = 0, & 
     iwp2_ta = 0, & 
     iwp2_ac = 0, & 
     iwp2_bp = 0, & 
     iwp2_pr1 = 0, & 
     iwp2_pr2 = 0, & 
     iwp2_pr3 = 0, & 
     iwp2_dp1 = 0, & 
     iwp2_dp2 = 0, &
     iwp2_pd = 0, & 
     iwp2_cl = 0, &
     iwp2_sf = 0

!$omp   threadprivate(iwp2_bt, iwp2_ma, iwp2_ta, iwp2_ac, iwp2_bp)
!$omp   threadprivate(iwp2_pr1, iwp2_pr2, iwp2_pr3)
!$omp   threadprivate(iwp2_dp1, iwp2_dp2)
!$omp   threadprivate(iwp2_pd, iwp2_cl, iwp2_sf)

  integer, public :: & 
     iwprtp_bt = 0,      & 
     iwprtp_ma = 0,      & 
     iwprtp_ta = 0,      & 
     iwprtp_tp = 0,      & 
     iwprtp_ac = 0,      & 
     iwprtp_bp = 0,      & 
     iwprtp_pr1 = 0,     & 
     iwprtp_pr2 = 0,     & 
     iwprtp_pr3 = 0,     & 
     iwprtp_dp1 = 0,     &
     iwprtp_mfl = 0,     &
     iwprtp_cl = 0,      & 
     iwprtp_sicl = 0,    & 
     iwprtp_pd = 0,      &
     iwprtp_forcing = 0, &
     iwprtp_mc = 0

!$omp   threadprivate(iwprtp_bt, iwprtp_ma, iwprtp_ta, iwprtp_tp)
!$omp   threadprivate(iwprtp_ac, iwprtp_bp, iwprtp_pr1, iwprtp_pr2)
!$omp   threadprivate(iwprtp_pr3, iwprtp_dp1, iwprtp_mfl, iwprtp_cl)
!$omp   threadprivate(iwprtp_sicl, iwprtp_pd, iwprtp_forcing, iwprtp_mc)

  integer, public :: & 
     iwpthlp_bt = 0,      & 
     iwpthlp_ma = 0,      & 
     iwpthlp_ta = 0,      & 
     iwpthlp_tp = 0,      & 
     iwpthlp_ac = 0,      & 
     iwpthlp_bp = 0,      & 
     iwpthlp_pr1 = 0,     & 
     iwpthlp_pr2 = 0,     & 
     iwpthlp_pr3 = 0,     & 
     iwpthlp_dp1 = 0,     &
     iwpthlp_mfl = 0,     &
     iwpthlp_cl = 0,      & 
     iwpthlp_sicl = 0,    &
     iwpthlp_forcing = 0, &
     iwpthlp_mc = 0

!$omp   threadprivate(iwpthlp_bt, iwpthlp_ma, iwpthlp_ta, iwpthlp_tp)
!$omp   threadprivate(iwpthlp_ac, iwpthlp_bp, iwpthlp_pr1, iwpthlp_pr2)
!$omp   threadprivate(iwpthlp_pr3, iwpthlp_dp1, iwpthlp_mfl, iwpthlp_cl)
!$omp   threadprivate(iwpthlp_sicl, iwpthlp_forcing, iwpthlp_mc)

!    Dr. Golaz's new variance budget terms
!    qt was changed to rt to avoid confusion

  integer, public :: & 
     irtp2_bt = 0,      & 
     irtp2_ma = 0,      & 
     irtp2_ta = 0,      & 
     irtp2_tp = 0,      & 
     irtp2_dp1 = 0,     & 
     irtp2_dp2 = 0,     & 
     irtp2_pd = 0,      & 
     irtp2_cl = 0,      &
     irtp2_sf = 0,      &
     irtp2_forcing = 0, &
     irtp2_mc = 0
     
!$omp   threadprivate(irtp2_bt, irtp2_ma, irtp2_ta, irtp2_tp, irtp2_dp1)
!$omp   threadprivate(irtp2_dp2, irtp2_pd, irtp2_cl, irtp2_sf, irtp2_forcing)
!$omp   threadprivate(irtp2_mc)

  integer, public :: & 
     ithlp2_bt = 0,      & 
     ithlp2_ma = 0,      & 
     ithlp2_ta = 0,      & 
     ithlp2_tp = 0,      & 
     ithlp2_dp1 = 0,     & 
     ithlp2_dp2 = 0,     & 
     ithlp2_pd = 0,      & 
     ithlp2_cl = 0,      &
     ithlp2_sf = 0,      &
     ithlp2_forcing = 0, &
     ithlp2_mc = 0

!$omp   threadprivate(ithlp2_bt, ithlp2_ma, ithlp2_ta, ithlp2_tp, ithlp2_dp1)
!$omp   threadprivate(ithlp2_dp2, ithlp2_pd, ithlp2_cl, ithlp2_sf)
!$omp   threadprivate(ithlp2_forcing, ithlp2_mc)

  integer, public :: & 
    irtpthlp_bt = 0,      & 
    irtpthlp_ma = 0,      & 
    irtpthlp_ta = 0,      & 
    irtpthlp_tp1 = 0,     & 
    irtpthlp_tp2 = 0,     & 
    irtpthlp_dp1 = 0,     & 
    irtpthlp_dp2 = 0,     & 
    irtpthlp_cl = 0,      &
    irtpthlp_sf = 0,      &
    irtpthlp_forcing = 0, &
    irtpthlp_mc = 0

!$omp   threadprivate(irtpthlp_bt, irtpthlp_ma, irtpthlp_ta)
!$omp   threadprivate(irtpthlp_tp1, irtpthlp_tp2, irtpthlp_dp1)
!$omp   threadprivate(irtpthlp_dp2, irtpthlp_cl, irtpthlp_sf, irtpthlp_forcing)
!$omp   threadprivate(irtpthlp_mc)

  integer, public :: & 
    iup2 = 0, & 
    ivp2 = 0

!$omp   threadprivate(iup2, ivp2)

  integer, public :: & 
    iup2_bt = 0, & 
    iup2_ta = 0, & 
    iup2_tp = 0, & 
    iup2_ma = 0, & 
    iup2_dp1 = 0, & 
    iup2_dp2 = 0, & 
    iup2_pr1 = 0, & 
    iup2_pr2 = 0, & 
    iup2_pd = 0, & 
    iup2_cl = 0, &
    iup2_sf = 0, &
    ivp2_bt = 0, & 
    ivp2_ta = 0, & 
    ivp2_tp = 0, & 
    ivp2_ma = 0, & 
    ivp2_dp1 = 0, & 
    ivp2_dp2 = 0, & 
    ivp2_pr1 = 0, & 
    ivp2_pr2 = 0, & 
    ivp2_pd = 0, & 
    ivp2_cl = 0, &
    ivp2_sf = 0

!$omp   threadprivate(iup2_bt, iup2_ta, iup2_tp, iup2_ma, iup2_dp1)
!$omp   threadprivate(iup2_dp2, iup2_pr1, iup2_pr2, iup2_cl, iup2_sf)
!$omp   threadprivate(ivp2_bt, ivp2_ta, ivp2_tp, ivp2_ma, ivp2_dp1)
!$omp   threadprivate(ivp2_dp2, ivp2_pr1, ivp2_pr2, ivp2_cl)
!$omp   threadprivate(iup2_pd, ivp2_pd, ivp2_sf)

!       Passive scalars.  Note that floating point roundoff may make
!       mathematically equivalent variables different values.
  integer,target, allocatable, dimension(:), public :: & 
    isclrprtp,           & ! sclr'(1)rt'     / rt'^2
    isclrp2,             & ! sclr'(1)^2      / rt'^2
    isclrpthvp,          & ! sclr'(1)th_v'   / rt'th_v' 
    isclrpthlp,          & ! sclr'(1)th_l'   / rt'th_l' 
    isclrprcp,           & ! sclr'(1)rc'     / rt'rc'
    iwpsclrp,            & ! w'slcr'(1)      / w'rt'
    iwp2sclrp,           & ! w'^2 sclr'(1)   / w'^2 rt'
    iwpsclrp2,           & ! w'sclr'(1)^2    / w'rt'^2
    iwpsclrprtp,         & ! w'sclr'(1)rt'   / w'rt'^2
    iwpsclrpthlp           ! w'sclr'(1)th_l' / w'rt'th_l'

!$omp   threadprivate(isclrprtp, isclrp2, isclrpthvp, isclrpthlp) 
!$omp   threadprivate(isclrprcp, iwpsclrp, iwp2sclrp, iwpsclrp2)
!$omp   threadprivate(iwpsclrprtp, iwpsclrpthlp)

  integer, target, allocatable, dimension(:), public :: & 
     iwpedsclrp ! eddy sclr'(1)w'

!$omp threadprivate(iwpedsclrp)

  ! Indices for statistics in stats_rad_zt file
  integer, public :: &
    iT_in_K_rad = 0, &
    ircil_rad = 0, &
    io3l_rad = 0, &
    irsm_rad = 0, &
    ircm_in_cloud_rad = 0, &
    icloud_frac_rad = 0, & 
    iice_supersat_frac_rad = 0, &
    iradht_rad = 0, &
    iradht_LW_rad = 0, &
    iradht_SW_rad = 0, &
    ip_in_mb_rad = 0, &
    isp_humidity_rad = 0

!$omp threadprivate( iT_in_K_rad, ircil_rad, io3l_rad, &
!$omp   irsm_rad, ircm_in_cloud_rad, icloud_frac_rad, &
!$omp   iice_supersat_frac_rad, &
!$omp   iradht_rad, iradht_LW_rad, iradht_SW_rad, &
!$omp   ip_in_mb_rad, isp_humidity_rad )

  ! Indices for statistics in stats_rad_zm file
  integer, public :: &
    iFrad_LW_rad = 0, &
    iFrad_SW_rad = 0, &
    iFrad_SW_up_rad = 0, &
    iFrad_LW_up_rad = 0, &
    iFrad_SW_down_rad = 0, &
    iFrad_LW_down_rad = 0

!$omp threadprivate(iFrad_LW_rad, iFrad_SW_rad, iFrad_SW_up_rad)
!$omp threadprivate(iFrad_LW_up_rad, iFrad_SW_down_rad, iFrad_LW_down_rad)

  ! Indices for statistics in stats_sfc file

  integer, public :: & 
    iustar = 0, &
    isoil_heat_flux = 0,&
    iveg_T_in_K = 0,&
    isfc_soil_T_in_K = 0, &
    ideep_soil_T_in_K = 0,& 
    ilh = 0, & 
    ish = 0, & 
    icc = 0, & 
    ilwp = 0, &
    ivwp = 0, &        ! nielsenb
    iiwp = 0, &        ! nielsenb
    iswp = 0, &        ! nielsenb
    irwp = 0, &
    iz_cloud_base = 0, & 
    iz_inversion = 0, & 
    iprecip_rate_sfc = 0,    &    ! Brian
    irain_flux_sfc = 0,   &    ! Brian
    irrm_sfc = 0, & ! Brian
    iwpthlp_sfc = 0
!$omp threadprivate(iustar, isoil_heat_flux, iveg_T_in_K, isfc_soil_T_in_K, ideep_soil_T_in_K, &
!$omp   ilh, ish, icc, ilwp, ivwp, iiwp, iswp, irwp, iz_cloud_base, iz_inversion, &
!$omp   iprecip_rate_sfc, irain_flux_sfc, irrm_sfc, &
!$omp   iwpthlp_sfc )

  integer, public :: &
    iwprtp_sfc = 0, &
    iupwp_sfc = 0, &
    ivpwp_sfc = 0, &
    ithlm_vert_avg = 0, &
    irtm_vert_avg = 0, &
    ium_vert_avg = 0, &
    ivm_vert_avg = 0, &
    iwp2_vert_avg = 0, & ! nielsenb
    iup2_vert_avg = 0, &
    ivp2_vert_avg = 0, &
    irtp2_vert_avg = 0, &
    ithlp2_vert_avg = 0, &
    iT_sfc         ! kcwhite
!$omp threadprivate(iwprtp_sfc, iupwp_sfc, ivpwp_sfc, &
!$omp   ithlm_vert_avg, irtm_vert_avg, ium_vert_avg, ivm_vert_avg, &
!$omp   iwp2_vert_avg, iup2_vert_avg, ivp2_vert_avg, irtp2_vert_avg, ithlp2_vert_avg, iT_sfc)

  integer, public :: & 
    iwp23_matrix_condt_num = 0, & 
    irtm_matrix_condt_num = 0, & 
    ithlm_matrix_condt_num = 0, & 
    irtp2_matrix_condt_num = 0, & 
    ithlp2_matrix_condt_num = 0, & 
    irtpthlp_matrix_condt_num = 0, & 
    iup2_vp2_matrix_condt_num = 0, & 
    iwindm_matrix_condt_num = 0
!$omp threadprivate(iwp23_matrix_condt_num, irtm_matrix_condt_num, ithlm_matrix_condt_num, &
!$omp   irtp2_matrix_condt_num, ithlp2_matrix_condt_num, irtpthlp_matrix_condt_num, &
!$omp   iup2_vp2_matrix_condt_num, iwindm_matrix_condt_num)

  integer, public :: & 
    imorr_snow_rate = 0

!$omp threadprivate( imorr_snow_rate)

  integer, public :: &
    irtm_spur_src = 0,    &
    ithlm_spur_src = 0

!$omp threadprivate(irtm_spur_src, ithlm_spur_src)

  integer, public :: &
    iSkw_velocity = 0, & ! Skewness velocity
    iwp3_zm = 0, &
    ia3_coef = 0, &
    ia3_coef_zt = 0
!$omp threadprivate(iSkw_velocity, iwp3_zm, ia3_coef, ia3_coef_zt)

  integer, public :: &
    iwp3_on_wp2 = 0, &  ! w'^3 / w'^2 [m/s]
    iwp3_on_wp2_zt = 0  ! w'^3 / w'^2 [m/s]
!$omp threadprivate(iwp3_on_wp2, iwp3_on_wp2_zt)

  integer, public :: & 
    ilh_morr_snow_rate = 0
!$omp threadprivate( ilh_morr_snow_rate )

  integer, public :: & 
    ilh_vwp = 0, &
    ilh_lwp = 0
!$omp threadprivate( ilh_vwp, ilh_lwp )


  integer, public :: &
    icloud_frac_refined = 0, &
    ircm_refined = 0
!$omp threadprivate( icloud_frac_refined, ircm_refined )

  integer, public :: &
    irtp2_from_chi = 0

!$omp threadprivate( irtp2_from_chi )

  ! Variables that contains all the statistics

  type (stats), target, public :: stats_zt,   &    ! stats_zt grid
                                  stats_zm,   &    ! stats_zm grid
                                  stats_lh_zt,  &  ! stats_lh_zt grid
                                  stats_lh_sfc,  & ! stats_lh_sfc grid
                                  stats_rad_zt,  & ! stats_rad_zt grid
                                  stats_rad_zm,  & ! stats_rad_zm grid
                                  stats_sfc        ! stats_sfc

!$omp threadprivate(stats_zt, stats_zm, stats_lh_zt, stats_lh_sfc)
!$omp threadprivate(stats_rad_zt, stats_rad_zm, stats_sfc)

  ! Scratch space

  real( kind = core_rknd ), dimension(:), allocatable, public :: &
    ztscr01, ztscr02, ztscr03, & 
    ztscr04, ztscr05, ztscr06, & 
    ztscr07, ztscr08, ztscr09, & 
    ztscr10, ztscr11, ztscr12, & 
    ztscr13, ztscr14, ztscr15, & 
    ztscr16, ztscr17, ztscr18, &
    ztscr19, ztscr20, ztscr21

!$omp threadprivate(ztscr01, ztscr02, ztscr03, ztscr04, ztscr05)
!$omp threadprivate(ztscr06, ztscr07, ztscr08, ztscr09, ztscr10)
!$omp threadprivate(ztscr11, ztscr12, ztscr13, ztscr14, ztscr15)
!$omp threadprivate(ztscr16, ztscr17, ztscr18, ztscr19, ztscr20)
!$omp threadprivate(ztscr21)

  real( kind = core_rknd ), dimension(:), allocatable, public :: &
    zmscr01, zmscr02, zmscr03, &
    zmscr04, zmscr05, zmscr06, & 
    zmscr07, zmscr08, zmscr09, & 
    zmscr10, zmscr11, zmscr12, & 
    zmscr13, zmscr14, zmscr15, &
    zmscr16, zmscr17

!$omp   threadprivate(zmscr01, zmscr02, zmscr03, zmscr04, zmscr05)
!$omp   threadprivate(zmscr06, zmscr07, zmscr08, zmscr09, zmscr10)
!$omp   threadprivate(zmscr11, zmscr12, zmscr13, zmscr14, zmscr15)
!$omp   threadprivate(zmscr16, zmscr17)

end module stats_variables
