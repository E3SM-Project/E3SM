!-------------------------------------------------------------------------------
! $Id: stats_variables.F90 5633 2012-01-19 01:00:48Z vondeyle@uwm.edu $
!-------------------------------------------------------------------------------

! Description:
!   Holds pointers and other variables for statistics to be written to 
!   GrADS files and netCDF files.
!-------------------------------------------------------------------------------
module stats_variables


  use stats_type, only:  & 
      stats ! Type

  use clubb_precision, only:  & 
      time_precision, & ! Variable
      core_rknd

  implicit none

  private ! Set Default Scope

  ! Sampling and output frequencies
  real(kind=time_precision), public :: &
    stats_tsamp, & ! Sampling interval   [s]
    stats_tout     ! Output interval     [s]

!$omp   threadprivate(stats_tsamp, stats_tout)

  logical, public ::  & 
    l_stats,  & ! Main flag to turn statistics on/off
    l_output_rad_files, & ! Flag to turn off radiation statistics output
    l_netcdf, & ! Output to NetCDF format
    l_grads     ! Output to GrADS format

!$omp   threadprivate(l_stats, l_netcdf, l_grads)

  logical, public :: & 
    l_stats_samp,   & ! Sample flag for current time step
    l_stats_last      ! Last time step of output period

!$omp   threadprivate(l_stats_samp, l_stats_last)

  character(len=200), public ::  & 
    fname_zt,  &    ! Name of the stats file for thermodynamic grid fields
    fname_zm,  &    ! Name of the stats file for momentum grid fields
    fname_rad_zt, & ! Name of the stats file for the zt radiation grid fields
    fname_rad_zm, & ! Name of the stats file for the zm radiation grid fields
    fname_sfc       ! Name of the stats file for surface only fields

!$omp   threadprivate(fname_zt, fname_zm, fname_rad_zt, fname_rad_zm, fname_sfc)

!       Indices for statistics in zt file

  integer, public :: & 
     ithlm, & 
     ithvm, & 
     irtm, & 
     ircm, &
     irvm, & 
     ium, & 
     ivm, & 
     iwm_zt, &
     iwm_zm, &
     ium_ref,&
     ivm_ref, & 
     iug, & 
     ivg, & 
     icloud_frac, &
     ircm_in_layer, &
     ircm_in_cloud, &
     icloud_cover, &
     ip_in_Pa, & 
     iexner, & 
     irho_ds_zt, &
     ithv_ds_zt, &
     iLscale, & 
     iwp3, & 
     iwpthlp2, & 
     iwp2thlp, & 
     iwprtp2, & 
     iwp2rtp, & 
     iLscale_up, & 
     iLscale_down, & 
     itau_zt, & 
     iKh_zt, & 
     iwp2thvp, & 
     iwp2rcp, & 
     iwprtpthlp, & 
     isigma_sqd_w_zt, & 
     irho, &
     irel_humidity

  integer, public :: & 
     iNcm,            & ! Brian
     iNcnm,           & 
     iNcm_in_cloud,   &
     iNc_activated,   &
     isnowslope,      & ! Adam Smith, 22 April 2008
     ised_rcm,        & ! Brian
     irsat,           & ! Brian
     irsati, & 
     irrainm,         & ! Brian
     im_vol_rad_rain,   & ! Brian
     im_vol_rad_cloud,  & ! COAMPS only. dschanen 6 Dec 2006
     irain_rate_zt,      & ! Brian
     iAKm,            & ! analytic Kessler.  Vince Larson 22 May 2005 
     iLH_AKm,         & ! LH Kessler.  Vince Larson  22 May 2005
     iradht,          & ! Radiative heating.
     iradht_LW,       & !   "           "   Long-wave component
     iradht_SW          !   "           "   Short-wave component

  integer, public :: & 
     iAKstd,     & 
     iAKstd_cld, & 
     iAKm_rcm, & 
     iAKm_rcc

!$omp threadprivate(ithlm, ithvm, irtm, ircm, irvm, ium, ivm, ium_ref, ivm_ref, &
!$omp   iwm_zt, iwm_zm, iug, ivg, icloud_frac, ircm_in_layer, ircm_in_cloud, icloud_cover, &
!$omp   ip_in_Pa, iexner, irho_ds_zt, ithv_ds_zt, iLscale, iwp3, &
!$omp   iwpthlp2, iwp2thlp, iwprtp2, iwp2rtp, iLscale_up, iLscale_down, &
!$omp   itau_zt, iKh_zt, iwp2thvp, iwp2rcp, iwprtpthlp, isigma_sqd_w_zt, &
!$omp   irho, irel_humidity, iNcm, iNcnm, isnowslope, &
!$omp   ised_rcm, irsat, irsati, irrainm, &
!$omp   im_vol_rad_rain, im_vol_rad_cloud, &
!$omp   irain_rate_zt, iAKm, iLH_AKm, &
!$omp   iradht, iradht_LW, iradht_SW, &
!$omp   iAKstd, iAKstd_cld, iAKm_rcm, iAKm_rcc)

  ! Skewness functions on zt grid
  integer, public :: &
    iC11_Skw_fnc

!$omp threadprivate(iC11_Skw_fnc)

  integer, public :: &
    icloud_frac_zm, &
    ircm_zm, &
    irtm_zm, &
    ithlm_zm

!$omp threadprivate(icloud_frac_zm, ircm_zm, irtm_zm, ithlm_zm)

  integer, public :: &
    iLH_rcm_avg

!$omp threadprivate(iLH_rcm_avg)

  integer, public :: & 
     iNrm,       & ! Rain droplet number concentration
     iNim,       & ! Ice number concentration
     iNsnowm,    & ! Snow number concentration
     iNgraupelm    ! Graupel number concentration
!$omp   threadprivate(iNrm, iNim, iNsnowm, iNgraupelm)

  integer, public :: & 
     iT_in_K      ! Absolute temperature
!$omp   threadprivate(iT_in_K)

  integer, public :: &
    ieff_rad_cloud, &
    ieff_rad_ice, &
    ieff_rad_snow, &
    ieff_rad_rain, &
    ieff_rad_graupel

!$omp   threadprivate(ieff_rad_cloud, ieff_rad_ice, ieff_rad_snow) 
!$omp   threadprivate(ieff_rad_rain, ieff_rad_graupel)

  integer, public :: & 
    irsnowm, & 
    irgraupelm, & 
    iricem, & 
    idiam,           & ! Diameter of ice crystal           [m]
    imass_ice_cryst, & ! Mass of a single ice crystal      [kg]
    ircm_icedfs,     & ! Change in liquid water due to ice [kg/kg/s]
    iu_T_cm         ! Fallspeed of ice crystal in cm/s  [cm s^{-1}]

!$omp   threadprivate(irsnowm, irgraupelm, iricem, idiam)
!$omp   threadprivate(imass_ice_cryst, ircm_icedfs, iu_T_cm)


  ! thlm/rtm budget terms
  integer, public :: & 
    irtm_bt,         & ! rtm total time tendency
    irtm_ma,         & ! rtm mean advect. term
    irtm_ta,         & ! rtm turb. advect. term
    irtm_forcing,    & ! rtm large scale forcing term
    irtm_mc,         & ! rtm change from microphysics
    irtm_sdmp,       & ! rtm change from sponge damping
    irvm_mc,         & ! rvm change from microphysics
    ircm_mc,         & ! rcm change from microphysics
    ircm_sd_mg_morr, & ! rcm sedimentation tendency
    irtm_mfl,        & ! rtm change due to monotonic flux limiter
    irtm_tacl,       & ! rtm correction from turbulent advection (wprtp) clipping
    irtm_cl,         & ! rtm clipping term
    irtm_pd,         & ! thlm postive definite adj term
    ithlm_bt,        & ! thlm total time tendency
    ithlm_ma,        & ! thlm mean advect. term
    ithlm_ta,        & ! thlm turb. advect. term
    ithlm_forcing,   & ! thlm large scale forcing term
    ithlm_sdmp,      & ! thlm change from sponge damping
    ithlm_mc,        & ! thlm change from microphysics
    ithlm_mfl,       & ! thlm change due to monotonic flux limiter
    ithlm_tacl,      & ! thlm correction from turbulent advection (wpthlp) clipping
    ithlm_cl           ! thlm clipping term

!$omp   threadprivate(irtm_bt, irtm_ma, irtm_ta, irtm_forcing, &
!$omp     irtm_mc, irtm_sdmp, irtm_mfl, irtm_tacl, irtm_cl, irtm_pd, &
!$omp     irvm_mc, ircm_mc, ircm_sd_mg_morr, &
!$omp     ithlm_bt, ithlm_ma, ithlm_ta, ithlm_forcing, &
!$omp     ithlm_mc, ithlm_sdmp, ithlm_mfl, ithlm_tacl, ithlm_cl)

  !monatonic flux limiter diagnostic terms
  integer, public :: &
    ithlm_mfl_min, &
    ithlm_mfl_max, &
    iwpthlp_entermfl, &
    iwpthlp_exit_mfl, &
    iwpthlp_mfl_min, &
    iwpthlp_mfl_max, &
    irtm_mfl_min, &
    irtm_mfl_max, &
    iwprtp_enter_mfl, &
    iwprtp_exit_mfl, &
    iwprtp_mfl_min, &
    iwprtp_mfl_max, &
    ithlm_enter_mfl, &
    ithlm_exit_mfl, &
    ithlm_old, &
    ithlm_without_ta, &
    irtm_enter_mfl, &
    irtm_exit_mfl, &
    irtm_old, &
    irtm_without_ta

!$omp   threadprivate(ithlm_mfl_min, ithlm_mfl_max, iwpthlp_entermfl)
!$omp   threadprivate(iwpthlp_exit_mfl, iwpthlp_mfl_min, iwpthlp_mfl_max)
!$omp   threadprivate(irtm_mfl_min, irtm_mfl_max, iwprtp_enter_mfl)
!$omp   threadprivate(iwprtp_exit_mfl, iwprtp_mfl_min, iwprtp_mfl_max)
!$omp   threadprivate(ithlm_enter_mfl, ithlm_exit_mfl, ithlm_old, ithlm_without_ta)
!$omp   threadprivate(irtm_enter_mfl, irtm_exit_mfl, irtm_old, irtm_without_ta)

  integer, public :: & 
     iwp3_bt, & 
     iwp3_ma, & 
     iwp3_ta, & 
     iwp3_tp, & 
     iwp3_ac, & 
     iwp3_bp1, & 
     iwp3_bp2, & 
     iwp3_pr1, & 
     iwp3_pr2, & 
     iwp3_dp1, &
     iwp3_4hd, & 
     iwp3_cl

!$omp   threadprivate(iwp3_bt, iwp3_ma, iwp3_ta, iwp3_tp, iwp3_ac, iwp3_bp1)
!$omp   threadprivate(iwp3_bp2, iwp3_pr1, iwp3_pr2, iwp3_dp1, iwp3_4hd, iwp3_cl)

  ! Rain mixing ratio budgets
  integer, public :: & 
     irrainm_bt, &
     irrainm_ma, &
     irrainm_sd, &
     irrainm_sd_morr, &
     irrainm_dff, &
     irrainm_cond, &
     irrainm_auto, &
     irrainm_accr, &
     irrainm_cond_adj, &
     irrainm_src_adj, &
     irrainm_mc, &
     irrainm_cl

!$omp   threadprivate(irrainm_bt, irrainm_ma, irrainm_sd, irrainm_sd_morr, irrainm_dff)
!$omp   threadprivate(irrainm_cond, irrainm_auto, irrainm_accr)
!$omp   threadprivate(irrainm_cond_adj, irrainm_src_adj, irrainm_mc, irrainm_cl)

  integer, public :: &
     iNrm_bt, &
     iNrm_ma, &
     iNrm_sd, &
     iNrm_dff, &
     iNrm_cond, &
     iNrm_auto, &
     iNrm_cond_adj, &
     iNrm_src_adj, &
     iNrm_mc, &
     iNrm_cl

!$omp   threadprivate(iNrm_bt, iNrm_ma, iNrm_sd, iNrm_dff, iNrm_cond)
!$omp   threadprivate(iNrm_auto, iNrm_cond_adj, iNrm_src_adj, iNrm_mc, iNrm_cl)


  ! Snow/Ice/Graupel mixing ratio budgets
  integer, public :: & 
     irsnowm_bt, & 
     irsnowm_ma, & 
     irsnowm_sd, & 
     irsnowm_sd_morr, &
     irsnowm_dff, & 
     irsnowm_mc, & 
     irsnowm_cl

!$omp   threadprivate(irsnowm_bt, irsnowm_ma, irsnowm_sd, irsnowm_sd_morr, irsnowm_dff)
!$omp   threadprivate(irsnowm_mc, irsnowm_cl)

  integer, public :: & 
     irgraupelm_bt, & 
     irgraupelm_ma, & 
     irgraupelm_sd, & 
     irgraupelm_sd_morr, &
     irgraupelm_dff, & 
     irgraupelm_mc, & 
     irgraupelm_cl

!$omp   threadprivate(irgraupelm_bt, irgraupelm_ma, irgraupelm_sd, irgraupelm_sd_morr)
!$omp   threadprivate(irgraupelm_dff, irgraupelm_mc, irgraupelm_cl)

  integer, public :: & 
     iricem_bt, & 
     iricem_ma, & 
     iricem_sd, & 
     iricem_sd_mg_morr, &
     iricem_dff, & 
     iricem_mc, & 
     iricem_cl

!$omp   threadprivate(iricem_bt, iricem_ma, iricem_sd, iricem_sd_mg_morr, iricem_dff)
!$omp   threadprivate(iricem_mc, iricem_cl)

  integer, public :: &
    iNsnowm_bt,  &
    iNsnowm_ma,  &
    iNsnowm_sd,  &
    iNsnowm_dff, &
    iNsnowm_mc,  &
    iNsnowm_cl

!$omp threadprivate(iNsnowm_bt, iNsnowm_ma, iNsnowm_sd, iNsnowm_dff, &
!$omp   iNsnowm_mc, iNsnowm_cl)

  integer, public :: &
    iNgraupelm_bt,  &
    iNgraupelm_ma,  &
    iNgraupelm_sd,  &
    iNgraupelm_dff, &
    iNgraupelm_mc,  &
    iNgraupelm_cl

!$omp threadprivate(iNgraupelm_bt, iNgraupelm_ma, iNgraupelm_sd, &
!$omp   iNgraupelm_dff, iNgraupelm_mc, iNgraupelm_cl)

  integer, public :: &
    iNim_bt,  &
    iNim_ma,  &
    iNim_sd,  &
    iNim_dff, &
    iNim_mc,  &
    iNim_cl

!$omp threadprivate(iNim_bt, iNim_ma, iNim_sd, iNim_dff, &
!$omp   iNim_mc, iNim_cl)

  integer, public :: &
    iNcm_bt,  &
    iNcm_ma,  &
    iNcm_dff, &
    iNcm_mc,  &
    iNcm_cl,  &
    iNcm_act

!$omp threadprivate(iNcm_bt, iNcm_ma, iNcm_dff, &
!$omp   iNcm_mc, iNcm_cl)


  ! Wind budgets
  integer, public :: & 
     ivm_bt, & 
     ivm_ma, & 
     ivm_ta, & 
     ivm_gf, & 
     ivm_cf, &
     ivm_f, &
     ivm_sdmp, &
     ivm_ndg

!$omp   threadprivate(ivm_bt, ivm_ma, ivm_ta, ivm_gf, ivm_cf, ivm_f, ivm_sdmp, ivm_ndg)

  integer, public :: & 
     ium_bt, & 
     ium_ma, & 
     ium_ta, & 
     ium_gf, & 
     ium_cf, & 
     ium_f, &
     ium_sdmp, &
     ium_ndg

!$omp   threadprivate(ium_bt, ium_ma, ium_ta, ium_gf, ium_cf, ium_f, ium_sdmp, ium_ndg)


  ! PDF parameters
  integer, public :: & 
     imixt_frac, & 
     iw1, & 
     iw2, & 
     ivarnce_w1, & 
     ivarnce_w2, & 
     ithl1, & 
     ithl2, & 
     ivarnce_thl1, & 
     ivarnce_thl2, & 
     irt1, & 
     irt2, & 
     ivarnce_rt1, & 
     ivarnce_rt2, & 
     irc1, & 
     irc2, & 
     irsl1, & 
     irsl2, & 
     icloud_frac1, & 
     icloud_frac2, & 
     is1, & 
     is2, & 
     istdev_s1, & 
     istdev_s2, & 
     irrtthl

!$omp   threadprivate(imixt_frac, iw1, iw2, ivarnce_w1, ivarnce_w2, ithl1, ithl2, ivarnce_thl1, &
!$omp     ivarnce_thl2, irt1, irt2, ivarnce_rt1, ivarnce_rt2, irc1, irc2, &
!$omp     irsl1, irsl2, icloud_frac1, icloud_frac2, is1, is2, istdev_s1, istdev_s2, &
!$omp     irrtthl)

  integer, public :: & 
     iwp2_zt, & 
     ithlp2_zt, & 
     iwpthlp_zt, & 
     iwprtp_zt, & 
     irtp2_zt, & 
     irtpthlp_zt, &
     iup2_zt, &
     ivp2_zt, &
     iupwp_zt, &
     ivpwp_zt

!$omp   threadprivate(iwp2_zt, ithlp2_zt, iwpthlp_zt, iwprtp_zt, irtp2_zt, irtpthlp_zt, &
!$omp     iup2_zt, ivp2_zt, iupwp_zt, ivpwp_zt)

  integer, public :: &
    is_mellor
!$omp threadprivate(is_mellor)

  integer, target, allocatable, dimension(:), public :: & 
    isclrm,    & ! Passive scalar mean (1)
    isclrm_f     ! Passive scalar forcing (1)

! Used to calculate clear-sky radiative fluxes.
  integer, public :: &
    ifulwcl, ifdlwcl, ifdswcl, ifuswcl

!$omp   threadprivate(isclrm, isclrm_f)

  integer, target, allocatable, dimension(:), public :: & 
    iedsclrm,   & ! Eddy-diff. scalar term (1)
    iedsclrm_f    ! Eddy-diffusivity scalar forcing (1)

!$omp   threadprivate(iedsclrm, iedsclrm_f)

  integer, public :: &
    iLH_thlm_mc,      & ! Latin hypercube estimate of thlm_mc
    iLH_rvm_mc,       & ! Latin hypercube estimate of rvm_mc
    iLH_rcm_mc,       & ! Latin hypercube estimate of rcm_mc
    iLH_Ncm_mc,       & ! Latin hypercube estimate of Ncm_mc
    iLH_rrainm_mc,    & ! Latin hypercube estimate of rrainm_mc
    iLH_Nrm_mc,       & ! Latin hypercube estimate of Nrm_mc
    iLH_rsnowm_mc,    & ! Latin hypercube estimate of rsnowm_mc
    iLH_Nsnowm_mc,    & ! Latin hypercube estimate of Nsnowm_mc
    iLH_rgraupelm_mc, & ! Latin hypercube estimate of rgraupelm_mc
    iLH_Ngraupelm_mc, & ! Latin hypercube estimate of Ngraupelm_mc
    iLH_ricem_mc,     & ! Latin hypercube estimate of ricem_mc
    iLH_Nim_mc          ! Latin hypercube estimate of Nim_mc
!$omp   threadprivate( iLH_thlm_mc, iLH_rvm_mc, iLH_rcm_mc, iLH_Ncm_mc, &
!$omp     iLH_rrainm_mc,  iLH_Nrm_mc, iLH_rsnowm_mc, iLH_Nsnowm_mc, &
!$omp     iLH_rgraupelm_mc, iLH_Ngraupelm_mc, iLH_ricem_mc, iLH_Nim_mc )

  integer, public :: &
    iLH_Vrr, & ! Latin hypercube estimate of rrainm sedimentation velocity
    iLH_VNr    ! Latin hypercube estimate of Nrm sedimentation velocity  
!$omp   threadprivate(iLH_Vrr,  iLH_VNr)

  integer, public :: &
    iLH_rrainm, &
    iLH_Nrm, &
    iLH_ricem, &
    iLH_Nim, &
    iLH_rsnowm, &
    iLH_Nsnowm, &
    iLH_rgraupelm, &
    iLH_Ngraupelm, &
    iLH_thlm, &
    iLH_rcm, &
    iLH_Ncm, &
    iLH_rvm, &
    iLH_wm, &
    iLH_cloud_frac

!$omp threadprivate(iLH_rrainm, iLH_Nrm, iLH_ricem, iLH_Nim, iLH_rsnowm, iLH_Nsnowm, &
!$omp   iLH_rgraupelm, iLH_Ngraupelm, &
!$omp   iLH_thlm, iLH_rcm, iLH_Ncm, iLH_rvm, iLH_wm, iLH_cloud_frac )

  integer, public :: &
    iLH_wp2_zt, &
    iLH_Nrp2_zt, &
    iLH_Ncp2_zt, &
    iLH_rcp2_zt, &
    iLH_rtp2_zt, &
    iLH_thlp2_zt, &
    iLH_rrainp2_zt

!$omp threadprivate(iLH_wp2_zt, iLH_Nrp2_zt, iLH_Ncp2_zt, iLH_rcp2_zt, iLH_rtp2_zt, &
!$omp               iLH_thlp2_zt, iLH_rrainp2_zt)

  integer, public :: &
    itp2_mellor_1, &
    itp2_mellor_2, &
    isptp_mellor_1, &
    isptp_mellor_2, &
    icorr_st_mellor1, &
    icorr_st_mellor2

!$omp threadprivate(itp2_mellor_1, itp2_mellor_2, isptp_mellor_1, &
!$omp   isptp_mellor_2, icorr_st_mellor1, icorr_st_mellor2)

  ! Indices for statistics in zm file
  integer, public :: & 
     iwp2, & 
     irtp2, & 
     ithlp2, & 
     irtpthlp, & 
     iwprtp, & 
     iwpthlp, & 
     iwp4, & 
     iwpthvp, & 
     irtpthvp, & 
     ithlpthvp, & 
     itau_zm, & 
     iKh_zm, & 
     iwprcp, & 
     ithlprcp, & 
     irtprcp, & 
     ircp2, & 
     iupwp, & 
     ivpwp, & 
     irho_zm, & 
     isigma_sqd_w, &
     irho_ds_zm, &
     ithv_ds_zm, &
     iem, & 
     ishear,     & ! Brian
     imean_w_up, &
     imean_w_down, &
     iFrad, & 
     iFrad_LW,   & ! Brian
     iFrad_SW,   & ! Brian
     iFrad_LW_up,   & 
     iFrad_SW_up,   & 
     iFrad_LW_down,   & 
     iFrad_SW_down,   & 
     iFprec,          & ! Brian
     iFcsed             ! Brian

!$omp   threadprivate(iwp2, irtp2, ithlp2, irtpthlp, iwprtp, iwpthlp)
!$omp   threadprivate(iwp4, iwpthvp, irtpthvp, ithlpthvp, itau_zm, iKh_zm)
!$omp   threadprivate(iwprcp, ithlprcp, irtprcp, ircp2, iupwp, ivpwp)
!$omp   threadprivate(irho_zm, isigma_sqd_w, irho_ds_zm, ithv_ds_zm, iem, ishear)
!$omp   threadprivate(iFrad, iFrad_LW, iFrad_SW, iFrad_SW_up, iFrad_SW_down)
!$omp   threadprivate(iFrad_LW_up, iFrad_LW_down, iFprec, iFcsed)

  ! Skewness Functions on zm grid
  integer, public :: &
    igamma_Skw_fnc,  &
    iC6rt_Skw_fnc,   &
    iC6thl_Skw_fnc,  &
    iC7_Skw_fnc,     &
    iC1_Skw_fnc

!$omp   threadprivate(igamma_Skw_fnc, iC6rt_Skw_fnc, iC6thl_Skw_fnc)
!$omp   threadprivate(iC7_Skw_fnc, iC1_Skw_fnc)

  ! Sedimentation velocities
  integer, public :: & 
    iVNr,    &
    iVrr,    &
    iVNc,    &
    iVrc,    &
    iVNsnow, &
    iVrsnow, &
    iVNice,  &
    iVrice,  &
    iVrgraupel

!$omp   threadprivate(iVNr, iVrr, iVNc, iVrc, iVNsnow, iVrsnow, iVNice, iVrice, iVrgraupel)

  integer, public :: & 
     iwp2_bt, & 
     iwp2_ma, & 
     iwp2_ta, & 
     iwp2_ac, & 
     iwp2_bp, & 
     iwp2_pr1, & 
     iwp2_pr2, & 
     iwp2_pr3, & 
     iwp2_dp1, & 
     iwp2_dp2, &
     iwp2_4hd, &
     iwp2_pd, & 
     iwp2_cl, &
     iwp2_sf

!$omp   threadprivate(iwp2_bt, iwp2_ma, iwp2_ta, iwp2_ac, iwp2_bp)
!$omp   threadprivate(iwp2_pr1, iwp2_pr2, iwp2_pr3)
!$omp   threadprivate(iwp2_dp1, iwp2_dp2, iwp2_4hd)
!$omp   threadprivate(iwp2_pd, iwp2_cl)

  integer, public :: & 
     iwprtp_bt, & 
     iwprtp_ma, & 
     iwprtp_ta, & 
     iwprtp_tp, & 
     iwprtp_ac, & 
     iwprtp_bp, & 
     iwprtp_pr1, & 
     iwprtp_pr2, & 
     iwprtp_pr3, & 
     iwprtp_dp1, &
     iwprtp_mfl, &
     iwprtp_cl, & 
     iwprtp_sicl, & 
     iwprtp_pd

!$omp   threadprivate(iwprtp_bt, iwprtp_ma, iwprtp_ta, iwprtp_tp)
!$omp   threadprivate(iwprtp_ac, iwprtp_bp, iwprtp_pr1, iwprtp_pr2)
!$omp   threadprivate(iwprtp_pr3, iwprtp_dp1, iwprtp_mfl, iwprtp_cl)
!$omp   threadprivate(iwprtp_sicl, iwprtp_pd)

  integer, public :: & 
     iwpthlp_bt, & 
     iwpthlp_ma, & 
     iwpthlp_ta, & 
     iwpthlp_tp, & 
     iwpthlp_ac, & 
     iwpthlp_bp, & 
     iwpthlp_pr1, & 
     iwpthlp_pr2, & 
     iwpthlp_pr3, & 
     iwpthlp_dp1, &
     iwpthlp_mfl, &
     iwpthlp_cl, & 
     iwpthlp_sicl

!$omp   threadprivate(iwpthlp_bt, iwpthlp_ma, iwpthlp_ta, iwpthlp_tp)
!$omp   threadprivate(iwpthlp_ac, iwpthlp_bp, iwpthlp_pr1, iwpthlp_pr2)
!$omp   threadprivate(iwpthlp_pr3, iwpthlp_dp1, iwpthlp_mfl, iwpthlp_cl)
!$omp   threadprivate(iwpthlp_sicl)

!    Dr. Golaz's new variance budget terms
!    qt was changed to rt to avoid confusion

  integer, public :: & 
     irtp2_bt, & 
     irtp2_ma, & 
     irtp2_ta, & 
     irtp2_tp, & 
     irtp2_dp1, & 
     irtp2_dp2, & 
     irtp2_pd, & 
     irtp2_cl, &
     irtp2_sf
     
!$omp   threadprivate(irtp2_bt, irtp2_ma, irtp2_ta, irtp2_tp)
!$omp   threadprivate(irtp2_dp1, irtp2_dp2, irtp2_pd, irtp2_cl)

  integer, public :: & 
     ithlp2_bt, & 
     ithlp2_ma, & 
     ithlp2_ta, & 
     ithlp2_tp, & 
     ithlp2_dp1, & 
     ithlp2_dp2, & 
     ithlp2_pd, & 
     ithlp2_cl, &
     ithlp2_sf

!$omp   threadprivate(ithlp2_bt, ithlp2_ma, ithlp2_ta, ithlp2_tp)
!$omp   threadprivate(ithlp2_dp1, ithlp2_dp2, ithlp2_pd, ithlp2_cl)

  integer, public :: & 
    irtpthlp_bt, & 
    irtpthlp_ma, & 
    irtpthlp_ta, & 
    irtpthlp_tp1, & 
    irtpthlp_tp2, & 
    irtpthlp_dp1, & 
    irtpthlp_dp2, & 
    irtpthlp_cl, &
    irtpthlp_sf

!$omp   threadprivate(irtpthlp_bt, irtpthlp_ma, irtpthlp_ta)
!$omp   threadprivate(irtpthlp_tp1, irtpthlp_tp2, irtpthlp_dp1)
!$omp   threadprivate(irtpthlp_dp2, irtpthlp_cl)

  integer, public :: & 
    iup2, & 
    ivp2

!$omp   threadprivate(iup2, ivp2)

  integer, public :: & 
    iup2_bt, & 
    iup2_ta, & 
    iup2_tp, & 
    iup2_ma, & 
    iup2_dp1, & 
    iup2_dp2, & 
    iup2_pr1, & 
    iup2_pr2, & 
    iup2_pd, & 
    iup2_cl, &
    iup2_sf, &
    ivp2_bt, & 
    ivp2_ta, & 
    ivp2_tp, & 
    ivp2_ma, & 
    ivp2_dp1, & 
    ivp2_dp2, & 
    ivp2_pr1, & 
    ivp2_pr2, & 
    ivp2_pd, & 
    ivp2_cl, &
    ivp2_sf

!$omp   threadprivate(iup2_bt, iup2_ta, iup2_tp, iup2_ma, iup2_dp1)
!$omp   threadprivate(iup2_dp2, iup2_pr1, iup2_pr2, iup2_cl)
!$omp   threadprivate(ivp2_bt, ivp2_ta, ivp2_tp, ivp2_ma, ivp2_dp1)
!$omp   threadprivate(ivp2_dp2, ivp2_pr1, ivp2_pr2, ivp2_cl)
!$omp   threadprivate(iup2_pd, ivp2_pd)

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

!$omp   threadprivate(iwpedsclrp)
  ! Indices for statistics in rad_zt file
  integer, public :: &
    iT_in_K_rad, &
    ircil_rad, &
    io3l_rad, &
    irsnowm_rad, &
    ircm_in_cloud_rad, &
    icloud_frac_rad, & 
    iradht_rad, &
    iradht_LW_rad, &
    iradht_SW_rad

!$omp threadprivate(iT_in_K_rad, ircil_rad, io3l_rad)
!$omp threadprivate(irsnowm_rad, ircm_in_cloud_rad, icloud_frac_rad)
!$omp threadprivate(iradht_rad, iradht_LW_rad, iradht_SW_rad)

  ! Indices for statistics in rad_zm file
  integer, public :: &
    iFrad_LW_rad, &
    iFrad_SW_rad, &
    iFrad_SW_up_rad, &
    iFrad_LW_up_rad, &
    iFrad_SW_down_rad, &
    iFrad_LW_down_rad

!$omp threadprivate(iFrad_LW_rad, iFrad_SW_rad, iFrad_SW_up_rad)
!$omp threadprivate(iFrad_LW_up_rad, iFrad_SW_down_rad, iFrad_LW_down_rad)

  ! Indices for statistics in sfc file

  integer, public :: & 
    iustar, &
    isoil_heat_flux,&
    iveg_T_in_K,&
    isfc_soil_T_in_K, &
    ideep_soil_T_in_K,& 
    ilh, & 
    ish, & 
    icc, & 
    ilwp, &
    ivwp, &        ! nielsenb
    iiwp, &        ! nielsenb
    iswp, &        ! nielsenb
    irwp, &
    iz_cloud_base, & 
    iz_inversion, & 
    irain_rate_sfc,    &    ! Brian
    irain_flux_sfc,   &    ! Brian
    irrainm_sfc, & ! Brian
    iwpthlp_sfc, &
    iwprtp_sfc, &
    iupwp_sfc, &
    ivpwp_sfc, &
    ithlm_vert_avg, &
    irtm_vert_avg, &
    ium_vert_avg, &
    ivm_vert_avg, &
    iwp2_vert_avg, & ! nielsenb
    iup2_vert_avg, &
    ivp2_vert_avg, &
    irtp2_vert_avg, &
    ithlp2_vert_avg, &
    iT_sfc         ! kcwhite

  integer, public :: & 
    iwp23_matrix_condt_num, & 
    irtm_matrix_condt_num, & 
    ithlm_matrix_condt_num, & 
    irtp2_matrix_condt_num, & 
    ithlp2_matrix_condt_num, & 
    irtpthlp_matrix_condt_num, & 
    iup2_vp2_matrix_condt_num, & 
    iwindm_matrix_condt_num

  integer, public :: & 
    imorr_rain_rate, &
    imorr_snow_rate

  integer, public :: &
    irtm_spur_src,    &
    ithlm_spur_src
!$omp threadprivate(iustar, isoil_heat_flux, iveg_T_in_K, isfc_soil_T_in_K, ideep_soil_T_in_K, &
!$omp   ilh, ish, icc, ilwp, ivwp, iiwp, iswp, irwp, iz_cloud_base, iz_inversion, &
!$omp   irain_rate_sfc, irain_flux_sfc, irrainm_sfc, &
!$omp   iwpthlp_sfc, iwprtp_sfc, iupwp_sfc, ivpwp_sfc, &
!$omp   ithlm_vert_avg, irtm_vert_avg, ium_vert_avg, ivm_vert_avg, &
!$omp   iwp2_vert_avg, iup2_vert_avg, ivp2_vert_avg, irtp2_vert_avg, ithlp2_vert_avg, iT_sfc, &
!$omp   iwp23_matrix_condt_num, irtm_matrix_condt_num, ithlm_matrix_condt_num, &
!$omp   irtp2_matrix_condt_num, ithlp2_matrix_condt_num, irtpthlp_matrix_condt_num, &
!$omp   iup2_vp2_matrix_condt_num, iwindm_matrix_condt_num, &
!$omp   imorr_rain_rate, imorr_snow_rate)

  integer, public :: &
    iSkw_velocity, & ! Skewness velocity
    iwp3_zm, &
    ia3_coef, &
    ia3_coef_zt
!$omp threadprivate(iSkw_velocity, iwp3_zm, ia3_coef, ia3_coef_zt)

  integer, public :: &
    iwp3_on_wp2, &  ! w'^3 / w'^2 [m/s]
    iwp3_on_wp2_zt  ! w'^3 / w'^2 [m/s]
!$omp threadprivate(iwp3_on_wp2, iwp3_on_wp2_zt)

  ! Variables that contains all the statistics

  type (stats), target, public :: zt,   &    ! zt grid
                                  zm,   &    ! zm grid
                                  rad_zt,  & ! rad_zt grid
                                  rad_zm,  & ! rad_zm grid
                                  sfc        ! sfc

!$omp   threadprivate(zt, zm, rad_zt, rad_zm, sfc)

  ! Scratch space

  real( kind = core_rknd ), dimension(:), allocatable, public :: &
    ztscr01, ztscr02, ztscr03, & 
    ztscr04, ztscr05, ztscr06, & 
    ztscr07, ztscr08, ztscr09, & 
    ztscr10, ztscr11, ztscr12, & 
    ztscr13, ztscr14, ztscr15, & 
    ztscr16, ztscr17, ztscr18, &
    ztscr19, ztscr20, ztscr21

!$omp   threadprivate(ztscr01, ztscr02, ztscr03, ztscr04, ztscr05)
!$omp   threadprivate(ztscr06, ztscr07, ztscr08, ztscr09, ztscr10)
!$omp   threadprivate(ztscr11, ztscr12, ztscr13, ztscr14, ztscr15)
!$omp   threadprivate(ztscr16, ztscr17, ztscr18, ztscr19, ztscr20)
!$omp   threadprivate(ztscr21)

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

  real( kind = core_rknd ), dimension(:), allocatable, public :: &
    radscr01, radscr02, radscr03, &
    radscr04, radscr05, radscr06, & 
    radscr07, radscr08, radscr09, & 
    radscr10, radscr11, radscr12, & 
    radscr13, radscr14, radscr15, &
    radscr16, radscr17

!$omp   threadprivate(radscr01, radscr02, radscr03, radscr04, radscr05)
!$omp   threadprivate(radscr06, radscr07, radscr08, radscr09, radscr10)
!$omp   threadprivate(radscr11, radscr12, radscr13, radscr14, radscr15)
!$omp   threadprivate(radscr16, radscr17)

end module stats_variables
