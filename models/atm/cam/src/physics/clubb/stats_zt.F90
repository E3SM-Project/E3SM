!-----------------------------------------------------------------------
! $Id: stats_zt.F90 5633 2012-01-19 01:00:48Z vondeyle@uwm.edu $

module stats_zt

  implicit none

  private ! Default Scope

  public :: stats_init_zt

! Constant parameters
  integer, parameter, public :: nvarmax_zt = 300 ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_zt( vars_zt, l_error )

! Description:
!   Initializes array indices for zt

! Note:
!   All code that is within subroutine stats_init_zt, including variable
!   allocation code, is not called if l_stats is false.  This subroutine is
!   called only when l_stats is true.

!-----------------------------------------------------------------------

    use constants_clubb, only:  &
        fstderr ! Constant(s)

    use stats_variables, only: & 
        ithlm,  & ! Variable(s)
        iT_in_K, & 
        ithvm, & 
        irtm, & 
        ircm, &
        irvm, & 
        ium, & 
        ivm, & 
        iwm_zt, & 
        ium_ref, &
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
        isigma_sqd_w_zt

    use stats_variables, only: & 
        irel_humidity, &
        irho, & 
        iNcm, &
        iNcm_in_cloud, &
        iNc_activated, &
        iNcnm, & 
        isnowslope, & 
        ised_rcm, & 
        irsat, & 
        irsati, & 
        irrainm, & 
        iNrm, & 
        irain_rate_zt, & 
        iAKm, & 
        iLH_AKm, & 
        iAKstd, & 
        iAKstd_cld, & 
        iAKm_rcm, & 
        iAKm_rcc, & 
        iradht, & 
        iradht_LW, & 
        iradht_SW, & 
        idiam, & 
        imass_ice_cryst, & 
        ircm_icedfs, & 
        iu_T_cm, & 
        im_vol_rad_rain, & 
        im_vol_rad_cloud, & 
        irsnowm, & 
        irgraupelm, & 
        iricem

    use stats_variables, only: & 
      ieff_rad_cloud, &
      ieff_rad_ice, &
      ieff_rad_snow, &
      ieff_rad_rain, &
      ieff_rad_graupel

    use stats_variables, only: & 
        irtm_bt, & 
        irtm_ma, & 
        irtm_ta, & 
        irtm_forcing, & 
        irtm_mc, &
        irtm_sdmp, &
        ircm_mc, &
        ircm_sd_mg_morr, &
        irvm_mc, &
        irtm_mfl, &
        irtm_tacl, & 
        irtm_cl, & 
        irtm_pd, & 
        ithlm_bt, & 
        ithlm_ma, & 
        ithlm_ta, & 
        ithlm_forcing, & 
        ithlm_mc, &
        ithlm_sdmp, &
        ithlm_mfl, &
        ithlm_tacl, &
        ithlm_cl, &
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

    ! Monotonic flux limiter diagnostic variables
    use stats_variables, only: &
        ithlm_mfl_min, &
        ithlm_mfl_max, &
        irtm_mfl_min, &
        irtm_mfl_max, &
        ithlm_enter_mfl, &
        ithlm_exit_mfl, &
        ithlm_old, &
        ithlm_without_ta, &
        irtm_enter_mfl, &
        irtm_exit_mfl, &
        irtm_old, &
        irtm_without_ta

    use stats_variables, only: & 
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
        irrainm_cl, & 
        iNrm_bt, & 
        iNrm_ma, & 
        iNrm_sd, & 
        iNrm_dff, & 
        iNrm_cond, & 
        iNrm_auto, & 
        iNrm_cond_adj, & 
        iNrm_src_adj, & 
        iNrm_mc, & 
        iNrm_cl, & 
        irsnowm_bt, & 
        irsnowm_ma, & 
        irsnowm_sd, &
        irsnowm_sd_morr, &
        irsnowm_dff

    use stats_variables, only: & 
        irsnowm_mc, & 
        irsnowm_cl, & 
        irgraupelm_bt, & 
        irgraupelm_ma, & 
        irgraupelm_sd, &
        irgraupelm_sd_morr, &
        irgraupelm_dff, & 
        irgraupelm_mc, & 
        irgraupelm_cl, & 
        iricem_bt, & 
        iricem_ma, & 
        iricem_sd, &
        iricem_sd_mg_morr, &
        iricem_dff, & 
        iricem_mc, & 
        iricem_cl, & 
        ivm_bt, & 
        ivm_ma, & 
        ivm_gf, & 
        ivm_cf, & 
        ivm_ta, &
        ivm_f, & 
        ivm_sdmp, &
        ivm_ndg, &
        ium_bt, & 
        ium_ma, & 
        ium_gf, & 
        ium_cf, & 
        ium_ta, &
        ium_f, &
        ium_sdmp, &
        ium_ndg

    use stats_variables, only: & 
        imixt_frac, & ! Variable(s) 
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
        irrtthl, & 
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

    use stats_variables, only: & 
        zt, & 
        isclrm, & 
        isclrm_f, & 
        iedsclrm, & 
        iedsclrm_f

    use stats_variables, only: & 
      iNsnowm, & ! Variable(s)
      iNrm, &
      iNgraupelm, &
      iNim, & 
      iNsnowm_bt, &
      iNsnowm_mc, &
      iNsnowm_ma, &
      iNsnowm_dff, &
      iNsnowm_sd, &
      iNsnowm_cl, &
      iNgraupelm_bt, &
      iNgraupelm_mc, &
      iNgraupelm_ma, &
      iNgraupelm_dff, &
      iNgraupelm_sd, &
      iNgraupelm_cl, &
      iNim_bt, &
      iNim_mc, &
      iNim_ma, &
      iNim_dff, &
      iNim_sd, &
      iNim_cl

    use stats_variables, only: & 
      iNcm_bt, &
      iNcm_mc, &
      iNcm_ma, &
      iNcm_dff, &
      iNcm_cl, &
      iNcm_act

    use stats_variables, only: & 
      ieff_rad_cloud, &
      ieff_rad_ice, &
      ieff_rad_snow, &
      ieff_rad_rain, &
      ieff_rad_graupel

    use stats_variables, only: &
      iLH_thlm_mc, &  ! Variable(s)
      iLH_rvm_mc, & 
      iLH_rcm_mc, & 
      iLH_Ncm_mc, & 
      iLH_rrainm_mc, & 
      iLH_Nrm_mc, & 
      iLH_rsnowm_mc, & 
      iLH_Nsnowm_mc, & 
      iLH_rgraupelm_mc, & 
      iLH_Ngraupelm_mc, & 
      iLH_ricem_mc, & 
      iLH_Nim_mc, & 
      iLH_Vrr, &
      iLH_VNr

    use stats_variables, only: &
      iLH_rcm_avg

    use stats_variables, only: &
      iLH_rrainm, & ! Variable(s)
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
      iLH_wp2_zt, &
      iLH_rcp2_zt, &
      iLH_rtp2_zt, &
      iLH_thlp2_zt, &
      iLH_rrainp2_zt, &
      iLH_Nrp2_zt, &
      iLH_Ncp2_zt, &
      iLH_cloud_frac

    use stats_variables, only: &
      iC11_Skw_fnc, & ! Variable(s)
      is_mellor, &
      iwp3_on_wp2_zt, &
      ia3_coef_zt
      

    use stats_type, only: & 
        stat_assign ! Procedure

    use parameters_model, only: &
        sclr_dim,&
        edsclr_dim
!use error_code, only: &
!    clubb_at_least_debug_level ! Function


    implicit none

    ! Input Variable
    character(len= * ), dimension(nvarmax_zt), intent(in) :: vars_zt

    ! Output Variable        
    logical, intent(inout) :: l_error

! Local Varables
    integer :: i, j, k

    logical :: l_found

    character(len=50) :: sclr_idx

! Default initialization for array indices for zt

    ithlm           = 0
    iT_in_K         = 0
    ithvm           = 0
    irtm            = 0
    ircm            = 0
    irvm            = 0
    ium             = 0
    ivm             = 0
    iwm_zt          = 0
    ium_ref         = 0
    ivm_ref         = 0
    iug             = 0
    ivg             = 0
    icloud_frac     = 0
    ircm_in_layer   = 0
    ircm_in_cloud   = 0
    icloud_cover    = 0
    ip_in_Pa        = 0
    iexner          = 0
    irho_ds_zt      = 0
    ithv_ds_zt      = 0
    iLscale         = 0
    iwp3            = 0
    iwpthlp2        = 0
    iwp2thlp        = 0
    iwprtp2         = 0
    iwp2rtp         = 0
    iLscale_up      = 0
    iLscale_down    = 0
    itau_zt         = 0
    iKh_zt          = 0
    iwp2thvp        = 0
    iwp2rcp         = 0
    iwprtpthlp      = 0
    isigma_sqd_w_zt = 0
    irho            = 0
    irel_humidity   = 0
    iNcm            = 0  ! Brian
    iNcm_in_cloud   = 0
    iNc_activated   = 0
    iNcnm           = 0
    iNim            = 0
    isnowslope      = 0  ! Adam Smith, 22 April 2008
    ised_rcm        = 0  ! Brian
    irsat           = 0  ! Brian
    irrainm      = 0  ! Brian
    irain_rate_zt   = 0  ! Brian
    iAKm            = 0  ! analytic Kessler.  Vince Larson 22 May 2005
    iLH_AKm         = 0  ! LH Kessler.  Vince Larson 22 May 2005
    iAKstd          = 0
    iAKstd_cld      = 0
    iAKm_rcm        = 0
    iAKm_rcc        = 0
    iradht          = 0
    iradht_LW       = 0
    iradht_SW       = 0

! Number concentrations
    iNsnowm    = 0  ! Adam Smith, 22 April 2008
    iNrm       = 0  ! Brian
    iNgraupelm = 0
    iNim       = 0

    idiam           = 0
    imass_ice_cryst = 0
    ircm_icedfs     = 0
    iu_T_cm         = 0

! From K&K microphysics
    im_vol_rad_rain  = 0  ! Brian
    im_vol_rad_cloud = 0

! From Morrison microphysics
    ieff_rad_cloud   = 0
    ieff_rad_ice     = 0
    ieff_rad_snow    = 0
    ieff_rad_rain    = 0
    ieff_rad_graupel = 0

    irsnowm         = 0
    irgraupelm      = 0
    iricem          = 0

    irtm_bt         = 0
    irtm_ma         = 0
    irtm_ta         = 0
    irtm_forcing    = 0
    irtm_sdmp       = 0
    irtm_mc         = 0
    ircm_mc         = 0 ! For the change due to COAMPS/Morrison microphysics
    ircm_sd_mg_morr = 0
    irvm_mc         = 0 ! For the change due to COAMPS/Morrison microphysics
    irtm_mfl        = 0
    irtm_tacl       = 0
    irtm_cl         = 0 ! Josh
    irtm_pd         = 0
    ithlm_bt        = 0
    ithlm_ma        = 0
    ithlm_ta        = 0
    ithlm_forcing   = 0
    ithlm_mc        = 0
    ithlm_sdmp      = 0
    ithlm_mfl       = 0
    ithlm_tacl      = 0
    ithlm_cl        = 0 ! Josh

    ithlm_mfl_min = 0
    ithlm_mfl_max = 0
    irtm_mfl_min = 0
    irtm_mfl_max = 0
    ithlm_enter_mfl = 0
    ithlm_exit_mfl = 0
    ithlm_old = 0
    ithlm_without_ta = 0
    irtm_enter_mfl = 0
    irtm_exit_mfl = 0
    irtm_old = 0
    irtm_without_ta = 0

    iwp3_bt       = 0
    iwp3_ma       = 0
    iwp3_ta       = 0
    iwp3_tp       = 0
    iwp3_ac       = 0
    iwp3_bp1      = 0
    iwp3_bp2      = 0
    iwp3_pr1      = 0
    iwp3_pr2      = 0
    iwp3_dp1      = 0
    iwp3_4hd      = 0
    iwp3_cl       = 0

    irrainm_bt       = 0
    irrainm_ma       = 0
    irrainm_sd       = 0
    irrainm_sd_morr  = 0
    irrainm_dff      = 0
    irrainm_cond     = 0
    irrainm_auto     = 0
    irrainm_accr     = 0
    irrainm_cond_adj = 0
    irrainm_src_adj  = 0
    irrainm_mc       = 0
    irrainm_cl       = 0

    iNrm_bt       = 0
    iNrm_ma       = 0
    iNrm_sd       = 0
    iNrm_dff      = 0
    iNrm_cond     = 0
    iNrm_auto     = 0
    iNrm_cond_adj = 0
    iNrm_src_adj  = 0
    iNrm_mc       = 0
    iNrm_cl       = 0

    iNsnowm_bt    = 0
    iNsnowm_ma    = 0
    iNsnowm_sd    = 0
    iNsnowm_dff   = 0
    iNsnowm_mc    = 0
    iNsnowm_cl    = 0

    iNim_bt    = 0
    iNim_ma    = 0
    iNim_sd    = 0
    iNim_dff   = 0
    iNim_mc    = 0
    iNim_cl    = 0

    iNcm_bt    = 0
    iNcm_ma    = 0
    iNcm_dff   = 0
    iNcm_mc    = 0
    iNcm_cl    = 0
    iNcm_act   = 0

    irsnowm_bt      = 0
    irsnowm_ma      = 0
    irsnowm_sd      = 0
    irsnowm_sd_morr = 0
    irsnowm_dff     = 0
    irsnowm_mc      = 0
    irsnowm_cl      = 0

    irgraupelm_bt      = 0
    irgraupelm_ma      = 0
    irgraupelm_sd      = 0
    irgraupelm_sd_morr = 0
    irgraupelm_dff     = 0
    irgraupelm_mc      = 0
    irgraupelm_cl      = 0

    iricem_bt         = 0
    iricem_ma         = 0
    iricem_sd         = 0
    iricem_sd_mg_morr = 0
    iricem_dff        = 0
    iricem_mc         = 0
    iricem_cl         = 0

    ivm_bt   = 0
    ivm_ma   = 0
    ivm_gf   = 0
    ivm_cf   = 0
    ivm_ta   = 0
    ivm_f    = 0
    ivm_sdmp = 0
    ivm_ndg   = 0

    ium_bt   = 0
    ium_ma   = 0
    ium_gf   = 0
    ium_cf   = 0
    ium_ta   = 0
    ium_f    = 0
    ium_sdmp = 0
    ium_ndg  = 0

    imixt_frac    = 0
    iw1           = 0
    iw2           = 0
    ivarnce_w1    = 0
    ivarnce_w2    = 0
    ithl1         = 0
    ithl2         = 0
    ivarnce_thl1  = 0
    ivarnce_thl2  = 0
    irt1          = 0
    irt2          = 0
    ivarnce_rt1   = 0
    ivarnce_rt2   = 0
    irc1          = 0
    irc2          = 0
    irsl1         = 0
    irsl2         = 0
    icloud_frac1  = 0
    icloud_frac2  = 0
    is1           = 0
    is2           = 0
    istdev_s1     = 0
    istdev_s2     = 0
    irrtthl       = 0

    is_mellor = 0

    iwp2_zt     = 0
    ithlp2_zt   = 0
    iwpthlp_zt  = 0
    iwprtp_zt   = 0
    irtp2_zt    = 0
    irtpthlp_zt = 0
    iup2_zt     = 0
    ivp2_zt     = 0
    iupwp_zt    = 0
    ivpwp_zt    = 0

    iLH_thlm_mc      = 0
    iLH_rvm_mc       = 0
    iLH_rcm_mc       = 0
    iLH_Ncm_mc       = 0
    iLH_rrainm_mc    = 0
    iLH_Nrm_mc       = 0
    iLH_rsnowm_mc    = 0
    iLH_Nsnowm_mc    = 0   
    iLH_rgraupelm_mc = 0 
    iLH_Ngraupelm_mc = 0 
    iLH_ricem_mc     = 0
    iLH_Nim_mc       = 0

    iLH_rcm_avg = 0

    iLH_Vrr = 0
    iLH_VNr = 0

    iLH_rrainm = 0
    iLH_ricem = 0
    iLH_rsnowm = 0
    iLH_rgraupelm = 0

    iLH_Nrm = 0
    iLH_Nim = 0
    iLH_Nsnowm = 0
    iLH_Ngraupelm = 0

    iLH_thlm = 0
    iLH_rcm = 0
    iLH_rvm = 0
    iLH_wm = 0
    iLH_cloud_frac = 0

    iLH_wp2_zt = 0
    iLH_rcp2_zt = 0
    iLH_rtp2_zt = 0
    iLH_thlp2_zt = 0
    iLH_rrainp2_zt = 0
    iLH_Nrp2_zt = 0
    iLH_Ncp2_zt = 0

    iC11_Skw_fnc = 0
    ia3_coef_zt = 0
    iwp3_on_wp2_zt = 0

    allocate( isclrm(1:sclr_dim) )
    allocate( isclrm_f(1:sclr_dim) )

    isclrm     = 0
    isclrm_f   = 0

    allocate( iedsclrm(1:edsclr_dim) )
    allocate( iedsclrm_f(1:edsclr_dim) )

    iedsclrm   = 0

    iedsclrm_f = 0

!     Assign pointers for statistics variables zt

    k = 1
    do i=1,zt%nn

      select case ( trim(vars_zt(i)) )
      case ('thlm')
        ithlm = k
        call stat_assign( ithlm, "thlm",  & 
              "Liquid water potential temperature (theta_l) [K]", "K", zt)
        k = k + 1

      case ('T_in_K')
        iT_in_K = k
        call stat_assign( iT_in_K, "T_in_K",  & 
              "Absolute temperature [K]", "K", zt )
        k = k + 1

      case ('thvm')
        ithvm = k
        call stat_assign( ithvm, "thvm", & 
              "Virtual potential temperature [K]", "K", zt )
        k = k + 1

      case ('rtm')
        irtm = k

        call stat_assign( irtm, "rtm", & 
              "Total (vapor+liquid) water mixing ratio [kg/kg]", "kg/kg", zt )

        !zt%f%var(irtm)%ptr => zt%x(:,k)
        !zt%f%var(irtm)%name = "rtm"
        !zt%f%var(irtm)%description
        != "total water mixing ratio (kg/kg)"
        !zt%f%var(irtm)%units = "kg/kg"

        k = k + 1

      case ('rcm')
        ircm = k
        call stat_assign( ircm, "rcm", & 
              "Cloud water mixing ratio [kg/kg]", "kg/kg", zt )
        k = k + 1
      case ('rvm')
        irvm = k
        call stat_assign( irvm, "rvm", & 
              "Vapor water mixing ratio [kg/kg]", "kg/kg", zt )
        k = k + 1
      case ('rel_humidity')
        irel_humidity = k
        call stat_assign( irel_humidity, "rel_humidity", & 
              "Relative humidity w.r.t. liquid (range [0,1]) [-]", "[-]", zt )
        k = k + 1
      case ('um')
        ium = k
        call stat_assign( ium, "um", & 
              "East-west (u) wind [m/s]", "m/s", zt )
        k = k + 1
      case ('vm')
        ivm = k
        call stat_assign( ivm,  "vm", & 
              "North-south (v) wind [m/s]", "m/s", zt )
        k = k + 1
      case ('wm_zt')
        iwm_zt = k
        call stat_assign( iwm_zt, "wm", & 
              "Vertical (w) wind [m/s]", "m/s", zt )
        k = k + 1
      case ('um_ref')
        ium_ref = k
        call stat_assign( ium_ref, "um_ref", & 
             "reference u wind (m/s) [m/s]", "m/s", zt)
        k = k + 1
      case ('vm_ref')
        ivm_ref = k
        call stat_assign( ivm_ref, "vm_ref", & 
             "reference v wind (m/s) [m/s]", "m/s", zt)
        k = k + 1
      case ('ug')
        iug = k
        call stat_assign( iug, "ug", & 
             "u geostrophic wind [m/s]", "m/s", zt)
        k = k + 1
      case ('vg')
        ivg = k
        call stat_assign( ivg, "vg", & 
             "v geostrophic wind [m/s]", "m/s", zt )
        k = k + 1
      case ('cloud_frac')
        icloud_frac = k
        call stat_assign( icloud_frac, "cloud_frac", & 
             "Cloud fraction (between 0 and 1) [-]", "count", zt )
        k = k + 1

      case ('rcm_in_layer')
        ircm_in_layer = k
        call stat_assign( ircm_in_layer, "rcm_in_layer", &
             "rcm in cloud layer [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('rcm_in_cloud')
        ircm_in_cloud = k
        call stat_assign( ircm_in_cloud, "rcm_in_cloud", &
             "in-cloud value of rcm (for microphysics) [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('cloud_cover')
        icloud_cover = k
        call stat_assign( icloud_cover, "cloud_cover", &
             "Cloud cover (between 0 and 1) [-]", "count", zt )
        k = k + 1
      case ('p_in_Pa')
        ip_in_Pa = k
        call stat_assign( ip_in_Pa, "p_in_Pa", & 
             "Pressure [Pa]", "Pa", zt )
        k = k + 1
      case ('exner')
        iexner = k
        call stat_assign( iexner, "exner", & 
             "Exner function = (p/p0)**(rd/cp) [-]", "count", zt )
        k = k + 1
      case ('rho_ds_zt')
        irho_ds_zt = k
        call stat_assign( irho_ds_zt, "rho_ds_zt", &
             "Dry, static, base-state density [kg/m^3]", "kg m^{-3}", zt )
        k = k + 1
      case ('thv_ds_zt')
        ithv_ds_zt = k
        call stat_assign( ithv_ds_zt, "thv_ds_zt", &
             "Dry, base-state theta_v [K]", "K", zt )
        k = k + 1
      case ('Lscale')
        iLscale = k
        call stat_assign( iLscale, "Lscale", & 
             "Mixing length [m]", "m", zt )
        k = k + 1
      case ('thlm_forcing')
        ithlm_forcing = k
        call stat_assign( ithlm_forcing, "thlm_forcing", & 
             "thlm budget: thetal forcing [K s^{-1}]", "K s^{-1}", zt )
        k = k + 1
      case ('thlm_mc')
        ithlm_mc = k
        call stat_assign( ithlm_mc, "thlm_mc", & 
             "Change in thlm due to microphysics (not in budget) [K s^{-1}]", "K s^{-1}", zt )
        k = k + 1
      case ('rtm_forcing')
        irtm_forcing = k
        call stat_assign( irtm_forcing, "rtm_forcing", & 
             "rtm budget: rt forcing [kg kg^{-1} s^{-1}]", "kg kg^{-1} s^{-1}", zt )
        k = k + 1

      case ('rtm_mc')
        irtm_mc = k
        call stat_assign( irtm_mc, "rtm_mc", & 
             "Change in rt due to microphysics (not in budget) [kg kg^{-1} s^{-1}]", &
             "kg kg^{-1} s^{-1}", zt )
        k = k + 1

      case ('rvm_mc')
        irvm_mc = k
        call stat_assign( irvm_mc, "rvm_mc", & 
             "Time tendency of vapor mixing ratio due to microphysics [kg/kg/s]", "kg/(kg s)", zt )
        k = k + 1

      case ('rcm_mc')
        ircm_mc = k
        call stat_assign( ircm_mc, "rcm_mc", & 
             "Time tendency of liquid water mixing ratio due microphysics [kg/kg/s]", &
             "kg/kg/s", zt )
        k = k + 1

      case ('rcm_sd_mg_morr')
        ircm_sd_mg_morr = k
        call stat_assign( ircm_sd_mg_morr, "rcm_sd_mg_morr", & 
             "rcm sedimentation when using morrision or MG microphysics (not in budget," &
             // " included in rcm_mc) [kg kg^{-1} s^{-1}]", "kg kg^{-1} s^{-1}", zt )
        k = k + 1

      case ('thlm_mfl_min')
        ithlm_mfl_min = k
        call stat_assign( ithlm_mfl_min, "thlm_mfl_min", & 
             "Minimum allowable thlm [K]", "K", zt )
        k = k + 1

      case ('thlm_mfl_max')
        ithlm_mfl_max = k
        call stat_assign( ithlm_mfl_max, "thlm_mfl_max", & 
             "Maximum allowable thlm [K]", "K", zt )
        k = k + 1

      case ('thlm_enter_mfl')
        ithlm_enter_mfl = k
        call stat_assign( ithlm_enter_mfl, "thlm_enter_mfl", & 
             "Thlm before flux-limiter [K]", "K", zt )
        k = k + 1

      case ('thlm_exit_mfl')
        ithlm_exit_mfl = k
        call stat_assign( ithlm_exit_mfl, "thlm_exit_mfl", & 
             "Thlm exiting flux-limiter [K]", "K", zt )
        k = k + 1

      case ('thlm_old')
        ithlm_old = k
        call stat_assign( ithlm_old, "thlm_old", & 
             "Thlm at previous timestep [K]", "K", zt )
        k = k + 1

      case ('thlm_without_ta')
        ithlm_without_ta = k
        call stat_assign( ithlm_without_ta, "thlm_without_ta", & 
             "Thlm without turbulent advection contribution [K]", "K", zt )
        k = k + 1

      case ('rtm_mfl_min')
        irtm_mfl_min = k
        call stat_assign( irtm_mfl_min, "rtm_mfl_min", & 
             "Minimum allowable rtm  [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('rtm_mfl_max')
        irtm_mfl_max = k
        call stat_assign( irtm_mfl_max, "rtm_mfl_max", & 
             "Maximum allowable rtm  [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('rtm_enter_mfl')
        irtm_enter_mfl = k
        call stat_assign( irtm_enter_mfl, "rtm_enter_mfl", & 
             "Rtm before flux-limiter  [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('rtm_exit_mfl')
        irtm_exit_mfl = k
        call stat_assign( irtm_exit_mfl, "rtm_exit_mfl", & 
             "Rtm exiting flux-limiter  [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('rtm_old')
        irtm_old = k
        call stat_assign( irtm_old, "rtm_old", & 
             "Rtm at previous timestep  [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('rtm_without_ta')
        irtm_without_ta = k
        call stat_assign( irtm_without_ta, "rtm_without_ta", & 
             "Rtm without turbulent advection contribution  [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('wp3')
        iwp3 = k
        call stat_assign( iwp3, "wp3", & 
             "w third order moment [m^3/s^3]", "m^3/s^3", zt )
        k = k + 1

      case ('wpthlp2')
        iwpthlp2 = k
        call stat_assign( iwpthlp2, "wpthlp2", & 
             "w'thl'^2 [(m K^2)/s]", "(m K^2)/s", zt )
        k = k + 1

      case ('wp2thlp')
        iwp2thlp = k
        call stat_assign( iwp2thlp, "wp2thlp", & 
             "w'^2thl' [(m^2 K)/s^2]", "(m^2 K)/s^2", zt )
        k = k + 1

      case ('wprtp2')
        iwprtp2 = k
        call stat_assign( iwprtp2, "wprtp2", & 
             "w'rt'^2 [(m kg)/(s kg)]", "(m kg)/(s kg)", zt )
        k = k + 1

      case ('wp2rtp')
        iwp2rtp = k
        call stat_assign( iwp2rtp, "wp2rtp", & 
             "w'^2rt' [(m^2 kg)/(s^2 kg)]", "(m^2 kg)/(s^2 kg)", zt )
        k = k + 1

      case ('Lscale_up')
        iLscale_up = k
        call stat_assign( iLscale_up, "Lscale_up", & 
             "Upward mixing length [m]", "m", zt )
        k = k + 1

      case ('Lscale_down')
        iLscale_down = k
        call stat_assign( iLscale_down, "Lscale_down", & 
             "Downward mixing length [m]", "m", zt )
        k = k + 1

      case ('tau_zt')
        itau_zt = k
        call stat_assign( itau_zt, "tau_zt", & 
             "Dissipation time [s]", "s", zt )
        k = k + 1

      case ('Kh_zt')
        iKh_zt = k
        call stat_assign( iKh_zt, "Kh_zt", & 
             "Eddy diffusivity [m^2/s]", "m^2/s", zt )
        k = k + 1

      case ('wp2thvp')
        iwp2thvp = k
        call stat_assign( iwp2thvp, "wp2thvp", & 
             "w'^2thv' [K m^2/s^2]", "K m^2/s^2", zt )
        k = k + 1

      case ('wp2rcp')
        iwp2rcp = k
        call stat_assign( iwp2rcp, "wp2rcp", & 
             "w'^2rc' [(m^2 kg)/(s^2 kg)]", "(m^2 kg)/(s^2 kg)", zt )
        k = k + 1

      case ('wprtpthlp')
        iwprtpthlp = k
        call stat_assign( iwprtpthlp, "wprtpthlp", & 
             "w'rt'thl' [(m kg K)/(s kg)]", "(m kg K)/(s kg)", zt )
        k = k + 1

      case ('sigma_sqd_w_zt')
        isigma_sqd_w_zt = k
        call stat_assign( isigma_sqd_w_zt, "sigma_sqd_w_zt", & 
             "Nondimensionalized w variance of Gaussian component [-]", "-", zt )
        k = k + 1

      case ('rho')
        irho = k
        call stat_assign( irho, "rho", & 
             "Air density [kg/m^3]", "kg m^{-3}", zt )
        k = k + 1

      case ('Ncm')           ! Brian
        iNcm = k
        call stat_assign( iNcm, "Ncm", & 
             "Cloud droplet number concentration [num/kg]", & 
             "num/kg", zt )
        k = k + 1

      case ('Ncm_in_cloud')
        iNcm_in_cloud = k

        call stat_assign( iNcm_in_cloud, "Ncm_in_cloud", &
             "In cloud droplet concentration [num/kg]", "num/kg", zt )

        k = k + 1

      case ('Nc_activated')
        iNc_activated = k

        call stat_assign( iNc_activated, "Nc_activated", &
             "Droplets activated by GFDL activation [num/kg]", "num/kg", zt )

        k = k + 1

      case ('Ncnm')
        iNcnm = k
        call stat_assign( iNcnm, "Ncnm", & 
             "Cloud nuclei number concentration [num/kg]", & 
             "num/kg", zt )
        k = k + 1

      case ('Nim')           ! Brian
        iNim = k
        call stat_assign( iNim, "Nim", & 
             "Ice crystal number concentration [num/kg]", & 
             "num/kg", zt )
        k = k + 1

      case ('snowslope')     ! Adam Smith, 22 April 2008
        isnowslope = k
        call stat_assign( isnowslope, "snowslope", & 
             "COAMPS microphysics snow slope parameter [1/m]", & 
             "1/m", zt )
        k = k + 1

      case ('Nsnowm')        ! Adam Smith, 22 April 2008
        iNsnowm = k
        call stat_assign( iNsnowm, "Nsnowm", & 
             "Snow particle number concentration [num/kg]", & 
             "num/kg", zt )
        k = k + 1

      case ('Ngraupelm')
        iNgraupelm = k
        call stat_assign( iNgraupelm, "Ngraupelm", & 
             "Graupel number concentration  [num/kg]", & 
             "num/kg", zt )
        k = k + 1

      case ('sed_rcm')       ! Brian
        ised_rcm = k
        call stat_assign( ised_rcm, "sed_rcm", & 
             "d(rcm)/dt due to cloud sedimentation [kg / (m^2 s)]", & 
             "kg / [m^2 s]", zt )
        k = k + 1

      case ('rsat')           ! Brian
        irsat = k
        call stat_assign( irsat, "rsat", & 
             "Saturation mixing ratio over liquid [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('rsati')
        irsati = k
        call stat_assign( irsati, "rsati", & 
             "Saturation mixing ratio over ice [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('rrainm')           ! Brian
        irrainm = k
        call stat_assign( irrainm, "rrainm", & 
             "Rain water mixing ratio [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('rsnowm')
        irsnowm = k
        call stat_assign( irsnowm, "rsnowm", & 
             "Snow water mixing ratio [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('ricem')
        iricem = k
        call stat_assign( iricem, "ricem", & 
             "Pristine ice water mixing ratio [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('rgraupelm')
        irgraupelm = k
        call stat_assign( irgraupelm, "rgraupelm", & 
             "Graupel water mixing ratio [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('Nrm')           ! Brian
        iNrm = k
        call stat_assign( iNrm, "Nrm", & 
             "Rain drop number concentration [num/kg]", & 
             "num/kg", zt )
        k = k + 1

      case ('m_vol_rad_rain')  ! Brian
        im_vol_rad_rain = k
        call stat_assign( im_vol_rad_rain, "mvrr", & 
             "Rain drop mean volume radius [m]", "m", zt )
        k = k + 1

      case ('m_vol_rad_cloud')
        im_vol_rad_cloud = k
        call stat_assign( im_vol_rad_cloud, "m_vol_rad_cloud", & 
             "Cloud drop mean volume radius [m]", "m", zt )
        k = k + 1

      case ('eff_rad_cloud')
        ieff_rad_cloud = k
        call stat_assign( ieff_rad_cloud, "eff_rad_cloud", & 
             "Cloud drop effective volume radius [microns]", "microns", zt )
        k = k + 1

      case ('eff_rad_ice')
        ieff_rad_ice = k

        call stat_assign( ieff_rad_ice, "eff_rad_ice", & 
             "Ice effective volume radius [microns]", "microns", zt )
        k = k + 1

      case ('eff_rad_snow')
        ieff_rad_snow = k
        call stat_assign( ieff_rad_snow, "eff_rad_snow", & 
             "Snow effective volume radius [microns]", "microns", zt )
        k = k + 1

      case ('eff_rad_rain')
        ieff_rad_rain = k
        call stat_assign( ieff_rad_rain, "eff_rad_rain", & 
             "Rain drop effective volume radius [microns]", "microns", zt )
        k = k + 1

      case ('eff_rad_graupel')
        ieff_rad_graupel = k
        call stat_assign( ieff_rad_graupel, "eff_rad_graupel", & 
             "Graupel effective volume radius [microns]", "microns", zt )
        k = k + 1

      case ('rain_rate_zt')     ! Brian
        irain_rate_zt = k

        call stat_assign( irain_rate_zt, "rain_rate_zt", & 
             "Rain rate [mm/day]", "mm/day", zt )
        k = k + 1

      case ( 'AKm' )           ! Vince Larson 22 May 2005
        iAKm = k
        call stat_assign( iAKm, "AKm", & 
             "Analytic Kessler ac [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ( 'LH_AKm' )       ! Vince Larson 22 May 2005
        iLH_AKm = k

        call stat_assign( iLH_AKm, "LH_AKm", & 
             "LH Kessler estimate  [kg/kg/s]", "kg/kg/s", zt )
        k = k + 1

      case ( 'AKstd' )
        iAKstd = k

        call stat_assign( iAKstd, "AKstd", & 
             "Exact standard deviation of gba Kessler [kg/kg/s]", "kg/kg/s", zt )
        k = k + 1

      case ( 'AKstd_cld' )
        iAKstd_cld = k

        call stat_assign( iAKstd_cld, "AKstd_cld", & 
             "Exact w/in cloud std of gba Kessler [kg/kg/s]", "kg/kg/s", zt )
        k = k + 1

      case ( 'AKm_rcm' )
        iAKm_rcm = k

        call stat_assign( iAKm_rcm, "AKm_rcm", & 
             "Exact local gba auto based on rcm [kg/kg/s]", "kg/kg/s", zt )
        k = k + 1

      case ( 'AKm_rcc' )
        iAKm_rcc = k

        call stat_assign( iAKm_rcc, "AKm_rcc", & 
             "Exact local gba based on w/in cloud rc [kg/kg/s]", "kg/kg/s", zt )
        k = k + 1

      case ('radht')
        iradht = k

        call stat_assign( iradht, "radht", & 
             "Total (sw+lw) radiative heating rate [K/s]", "K/s", zt )
        k = k + 1

      case ('radht_LW')
        iradht_LW = k

        call stat_assign( iradht_LW, "radht_LW", & 
             "Long-wave radiative heating rate [K/s]", "K/s", zt )

        k = k + 1

      case ('radht_SW')
        iradht_SW = k
        call stat_assign( iradht_SW, "radht_SW", & 
             "Short-wave radiative heating rate [K/s]", "K/s", zt )
        k = k + 1

      case ('diam')
        idiam = k

        call stat_assign( idiam, "diam", & 
             "Ice crystal diameter [m]", "m", zt )
        k = k + 1

      case ('mass_ice_cryst')
        imass_ice_cryst = k
        call stat_assign( imass_ice_cryst, "mass_ice_cryst", & 
             "Mass of a single ice crystal [kg]", "kg", zt )
        k = k + 1

      case ('rcm_icedfs')

        ircm_icedfs = k
        call stat_assign( ircm_icedfs, "rcm_icedfs", & 
             "Change in liquid due to ice [kg/kg/s]", "kg/kg/s", zt )
        k = k + 1

      case ('u_T_cm')
        iu_T_cm = k
        call stat_assign( iu_T_cm, "u_T_cm", & 
             "Ice crystal fallspeed [cm s^{-1}]", "cm s^{-1}", zt )
        k = k + 1

      case ('rtm_bt')
        irtm_bt = k

        call stat_assign( irtm_bt, "rtm_bt", & 
             "rtm budget: rtm time tendency [kg kg^{-1} s^{-1}]", "kg kg^{-1} s^{-1}", zt)
        k = k + 1

      case ('rtm_ma')
        irtm_ma = k

        call stat_assign( irtm_ma, "rtm_ma", & 
             "rtm budget: rtm vertical mean advection [kg kg^{-1} s^{-1}]", &
             "kg kg^{-1} s^{-1}", zt)
        k = k + 1

      case ('rtm_ta')
        irtm_ta = k

        call stat_assign( irtm_ta, "rtm_ta", & 
             "rtm budget: rtm turbulent advection [kg kg^{-1} s^{-1}]", "kg kg^{-1} s^{-1}", zt)
        k = k + 1

      case ('rtm_mfl')
        irtm_mfl = k

        call stat_assign( irtm_mfl, "rtm_mfl", &
              "rtm budget: rtm correction due to monotonic flux limiter [kg kg^{-1} s^{-1}]", &
              "kg kg^{-1} s^{-1}", zt)
        k = k + 1

      case ('rtm_tacl')
        irtm_tacl = k

        call stat_assign( irtm_tacl, "rtm_tacl", & 
             "rtm budget: rtm correction due to ta term (wprtp) clipping [kg kg^{-1} s^{-1}]", &
             "kg kg^{-1} s^{-1}", zt)

        k = k + 1

      case ('rtm_cl')
        irtm_cl = k

        call stat_assign( irtm_cl, "rtm_cl", & 
             "rtm budget: rtm clipping [kg kg^{-1} s^{-1}]", "kg kg^{-1} s^{-1}", zt)

        k = k + 1
      case ('rtm_sdmp')
        irtm_sdmp = k

        call stat_assign( irtm_sdmp, "rtm_sdmp", & 
             "rtm budget: rtm correction due to sponge damping [kg kg^{-1} s^{-1}]", &
             "kg kg^{-1} s^{-1}", zt)
        k = k + 1


      case ('rtm_pd')
        irtm_pd = k

        call stat_assign( irtm_pd, "rtm_pd", & 
             "rtm budget: rtm positive definite adjustment [kg kg^{-1} s^{-1}]", &
             "kg kg^{-1} s^{-1}", zt)

        k = k + 1

      case ('thlm_bt')
        ithlm_bt = k

        call stat_assign( ithlm_bt, "thlm_bt", & 
             "thlm budget: thlm time tendency [K s^{-1}]", "K s^{-1}", zt)
        k = k + 1

      case ('thlm_ma')
        ithlm_ma = k

        call stat_assign( ithlm_ma, "thlm_ma", & 
             "thlm budget: thlm vertical mean advection [K s^{-1}]", "K s^{-1}", zt)
        k = k + 1

      case ('thlm_sdmp')
        ithlm_sdmp = k

        call stat_assign( ithlm_sdmp, "thlm_sdmp", & 
             "thlm budget: thlm correction due to sponge damping [K s^{-1}]", "K s^{-1}", zt)
        k = k + 1


      case ('thlm_ta')
        ithlm_ta = k

        call stat_assign( ithlm_ta, "thlm_ta", & 
             "thlm budget: thlm turbulent advection [K s^{-1}]", "K s^{-1}", zt)
        k = k + 1

      case ('thlm_mfl')
        ithlm_mfl = k

        call stat_assign( ithlm_mfl, "thlm_mfl", &
              "thlm budget: thlm correction due to monotonic flux limiter [K s^{-1}]", &
              "K s^{-1}", zt)
        k = k + 1

      case ('thlm_tacl')
        ithlm_tacl = k

        call stat_assign( ithlm_tacl, "thlm_tacl", & 
              "thlm budget: thlm correction due to ta term (wpthlp) clipping [K s^{-1}]", &
              "K s^{-1}", zt)
        k = k + 1

      case ('thlm_cl')
        ithlm_cl = k

        call stat_assign( ithlm_cl, "thlm_cl", & 
              "thlm budget: thlm_cl [K s^{-1}]", "K s^{-1}", zt)
        k = k + 1

      case ('wp3_bt')
        iwp3_bt = k

        call stat_assign( iwp3_bt, "wp3_bt", & 
             "wp3 budget: wp3 time tendency [m^{3} s^{-4}]", "m^{3} s^{-4}", zt )
        k = k + 1

      case ('wp3_ma')
        iwp3_ma = k

        call stat_assign( iwp3_ma, "wp3_ma", & 
             "wp3 budget: wp3 vertical mean advection [m^{3} s^{-4}]", "m^{3} s^{-4}", zt )
        k = k + 1

      case ('wp3_ta')
        iwp3_ta = k

        call stat_assign( iwp3_ta, "wp3_ta", & 
             "wp3 budget: wp3 turbulent advection [m^{3} s^{-4}]", "m^{3} s^{-4}", zt )

        k = k + 1

      case ('wp3_tp')
        iwp3_tp = k
        call stat_assign( iwp3_tp, "wp3_tp", & 
             "wp3 budget: wp3 turbulent transport [m^{3} s^{-4}]", "m^{3} s^{-4}", zt )
        k = k + 1

      case ('wp3_ac')
        iwp3_ac = k
        call stat_assign( iwp3_ac, "wp3_ac", & 
             "wp3 budget: wp3 accumulation term [m^{3} s^{-4}]", "m^{3} s^{-4}", zt )
        k = k + 1

      case ('wp3_bp1')
        iwp3_bp1 = k
        call stat_assign( iwp3_bp1, "wp3_bp1", & 
             "wp3 budget: wp3 buoyancy production [m^{3} s^{-4}]", "m^{3} s^{-4}", zt )
        k = k + 1

      case ('wp3_bp2')
        iwp3_bp2 = k
        call stat_assign( iwp3_bp2, "wp3_bp2", & 
             "wp3 budget: wp3 2nd buoyancy production term [m^{3} s^{-4}]", "m^{3} s^{-4}", zt )
        k = k + 1

      case ('wp3_pr1')
        iwp3_pr1 = k
        call stat_assign( iwp3_pr1, "wp3_pr1", & 
             "wp3 budget: wp3 pressure term 1 [m^{3} s^{-4}]", "m^{3} s^{-4}", zt )
        k = k + 1

      case ('wp3_pr2')
        iwp3_pr2 = k
        call stat_assign( iwp3_pr2, "wp3_pr2", & 
             "wp3 budget: wp3 pressure term 2 [m^{3} s^{-4}]", "m^{3} s^{-4}", zt )

        k = k + 1

      case ('wp3_dp1')
        iwp3_dp1 = k
        call stat_assign( iwp3_dp1, "wp3_dp1", & 
             "wp3 budget: wp3 dissipation term 1 [m^{3} s^{-4}]", "m^{3} s^{-4}", zt )
        k = k + 1

      case ('wp3_4hd')
        iwp3_4hd = k
        call stat_assign( iwp3_4hd, "wp3_4hd", & 
             "wp3 budget: wp3 4th-order hyper-diffusion [m^{3} s^{-4}]", "m^{3} s^{-4}", zt )
        k = k + 1

      case ('wp3_cl')
        iwp3_cl = k
        call stat_assign( iwp3_cl, "wp3_cl", & 
             "wp3 budget: wp3 clipping term [m^{3} s^{-4}]", "m^{3} s^{-4}", zt )
        k = k + 1

      case ('rrainm_bt')
        irrainm_bt = k
        call stat_assign( irrainm_bt, "rrainm_bt", & 
             "rrainm budget: rrainm time tendency [kg kg^{-1} s^{-1}]", "kg kg^{-1} s^{-1}", zt )
        k = k + 1

      case ('rrainm_ma')
        irrainm_ma = k

        call stat_assign( irrainm_ma, "rrainm_ma", & 
             "rrainm budget: rrainm vertical mean advection [kg kg^{-1} s^{-1}]", &
             "kg kg^{-1} s^{-1}", zt )
        k = k + 1

      case ('rrainm_sd')
        irrainm_sd = k

        call stat_assign( irrainm_sd, "rrainm_sd", & 
             "rrainm budget: rrainm sedimentation [kg kg^{-1} s^{-1}]", "kg kg^{-1} s^{-1}", zt )
        k = k + 1

      case ('rrainm_sd_morr')
        irrainm_sd_morr = k

        call stat_assign( irrainm_sd_morr, "rrainm_sd_morr", & 
             "rrainm sedimentation when using morrision microphysics (not in budget, included" &
             // " in rrainm_mc) [kg kg^{-1} s^{-1}]", "kg kg^{-1} s^{-1}", zt )
        k = k + 1

      case ('rrainm_dff')
        irrainm_dff = k

        call stat_assign( irrainm_dff, "rrainm_dff", & 
             "rrainm budget: rrainm diffusion [kg kg^{-1} s^{-1}]", "kg kg^{-1} s^{-1}", zt )
        k = k + 1

      case ('rrainm_cond')
        irrainm_cond = k

        call stat_assign( irrainm_cond, "rrainm_cond", & 
             "rrainm condensation/evaporation [kg kg^{-1} s^{-1}]", &
             "kg kg^{-1} s^{-1}", zt )
        k = k + 1

      case ('rrainm_auto')
        irrainm_auto = k

        call stat_assign( irrainm_auto, "rrainm_auto", & 
             "rrainm autoconversion [kg kg^{-1} s^{-1}]", "kg kg^{-1} s^{-1}", zt )
        k = k + 1

      case ('rrainm_accr')
        irrainm_accr = k
        call stat_assign( irrainm_accr, "rrainm_accr", & 
             "rrainm accretion [kg kg^{-1} s^{-1}]", "kg kg^{-1} s^{-1}", zt )
        k = k + 1

      case ('rrainm_cond_adj')
        irrainm_cond_adj = k

        call stat_assign( irrainm_cond_adj, "rrainm_cond_adj", & 
             "rrainm budget: rrainm cond/evap adjustment due to over-evaporation " // &
             "[kg kg^{-1} s^{-1}]", "kg kg^{-1} s^{-1}", zt )
        k = k + 1

      case ('rrainm_src_adj')
        irrainm_src_adj = k

        call stat_assign( irrainm_src_adj, "rrainm_src_adj", & 
             "rrainm budget: rrainm source term adjustment due to over-depletion " // &
             "[kg kg^{-1} s^{-1}]", "kg kg^{-1} s^{-1}", zt )
        k = k + 1

      case ('rrainm_cl')
        irrainm_cl = k
        call stat_assign( irrainm_cl, "rrainm_cl", & 
             "rrainm budget: rrainm clipping term [kg kg^{-1} s^{-1}]", "kg kg^{-1} s^{-1}", zt )
        k = k + 1

      case ('rrainm_mc')
        irrainm_mc = k

        call stat_assign( irrainm_mc, "rrainm_mc", & 
             "rrainm budget: Change in rrainm due to microphysics [kg kg^{-1} s^{-1}]", &
             "kg kg^{-1} s^{-1}", zt )
        k = k + 1

      case ('Nrm_bt')
        iNrm_bt = k
        call stat_assign( iNrm_bt, "Nrm_bt", & 
             "Nrm budget: Nrm time tendency [(num/kg)/s]", "(num/kg)/s", zt )

        k = k + 1

      case ('Nrm_ma')
        iNrm_ma = k

        call stat_assign( iNrm_ma, "Nrm_ma", & 
             "Nrm budget: Nrm vertical mean advection [(num/kg)/s]", "(num/kg)/s", zt )
        k = k + 1

      case ('Nrm_sd')
        iNrm_sd = k

        call stat_assign( iNrm_sd, "Nrm_sd", & 
             "Nrm budget: Nrm sedimentation [(num/kg)/s]", "(num/kg)/s", zt )

        k = k + 1

      case ('Nrm_dff')
        iNrm_dff = k
        call stat_assign( iNrm_dff, "Nrm_dff", & 
             "Nrm budget: Nrm diffusion [(num/kg)/s]", "(num/kg)/s", zt )

        k = k + 1

      case ('Nrm_cond')
        iNrm_cond = k

        call stat_assign( iNrm_cond, "Nrm_cond", & 
             "Change in Nrm due to condensation [(num/kg)/s]", "(num/kg)/s", zt )
        k = k + 1

      case ('Nrm_auto')
        iNrm_auto = k

        call stat_assign( iNrm_auto, "Nrm_auto", & 
             "Change in Nrm due to autoconversion [(num/kg)/s]", "(num/kg)/s", zt )

        k = k + 1

      case ('Nrm_cond_adj')
        iNrm_cond_adj = k

        call stat_assign( iNrm_cond_adj, "Nrm_cond_adj", & 
             "Nrm budget: Nrm cond/evap adjustment due to over-evaporation [(num/kg)/s]", & 
             "(num/kg)/s", zt )
        k = k + 1

      case ('Nrm_src_adj')
        iNrm_src_adj = k

        call stat_assign( iNrm_src_adj, "Nrm_src_adj", & 
             "Nrm budget: Nrm source term adjustment due to over-depletion [(num/kg)/s]", &
             "(num/kg)/s", zt )
        k = k + 1

      case ('Nrm_cl')
        iNrm_cl = k
        call stat_assign( iNrm_cl, "Nrm_cl", & 
             "Nrm budget: Nrm clipping term [(num/kg)/s]", "(num/kg)/s", zt )
        k = k + 1

      case ('Nrm_mc')
        iNrm_mc = k
        call stat_assign( iNrm_mc, "Nrm_mc", & 
             "Nrm budget: Change in Nrm due to microphysics (Not in budget) [(num/kg)/s]", &
             "(num/kg)/s", zt )
        k = k + 1

      case ('rsnowm_bt')
        irsnowm_bt = k
        call stat_assign( irsnowm_bt, "rsnowm_bt", & 
             "rsnowm budget: rsnowm time tendency [(kg/kg)/s]", "(kg/kg)/s", zt )

        k = k + 1

      case ('rsnowm_ma')
        irsnowm_ma = k

        call stat_assign( irsnowm_ma, "rsnowm_ma", & 
             "rsnowm budget: rsnowm vertical mean advection [(kg/kg)/s]", "(kg/kg)/s", zt )
        k = k + 1

      case ('rsnowm_sd')
        irsnowm_sd = k
        call stat_assign( irsnowm_sd, "rsnowm_sd", & 
             "rsnowm budget: rsnowm sedimentation [(kg/kg)/s]", "(kg/kg)/s", zt )
        k = k + 1

      case ('rsnowm_sd_morr')
        irsnowm_sd_morr = k
        call stat_assign( irsnowm_sd_morr, "rsnowm_sd_morr", & 
             "rsnowm sedimentation when using morrison microphysics (Not in budget, included in" &
             // " rsnowm_mc) [(kg/kg)/s]", "(kg/kg)/s", zt )
        k = k + 1

      case ('rsnowm_dff')
        irsnowm_dff = k

        call stat_assign( irsnowm_dff, "rsnowm_dff", & 
             "rsnowm budget: rsnowm diffusion [(kg/kg)/s]", "(kg/kg)/s", zt )
        k = k + 1

      case ('rsnowm_mc')
        irsnowm_mc = k

        call stat_assign( irsnowm_mc, "rsnowm_mc", & 
             "rsnowm budget: Change in rsnowm due to microphysics [(kg/kg)/s]", "(kg/kg)/s", zt )
        k = k + 1

      case ('rsnowm_cl')
        irsnowm_cl = k

        call stat_assign( irsnowm_cl, "rsnowm_cl", & 
             "rsnowm budget: rsnowm clipping term [(kg/kg)/s]", "(kg/kg)/s", zt )
        k = k + 1

      case ('Nsnowm_bt')
        iNsnowm_bt = k
        call stat_assign( iNsnowm_bt, "Nsnowm_bt", & 
             "Nsnowm budget: [(num/kg)/s]", "(num/kg)/s", zt )

        k = k + 1

      case ('Nsnowm_ma')
        iNsnowm_ma = k

        call stat_assign( iNsnowm_ma, "Nsnowm_ma", & 
             "Nsnowm budget: Nsnowm mean advection [(num/kg)/s]", "(num/kg)/s", zt )
        k = k + 1

      case ('Nsnowm_sd')
        iNsnowm_sd = k

        call stat_assign( iNsnowm_sd, "Nsnowm_sd", & 
             "Nsnowm budget: Nsnowm sedimentation [(num/kg)/s]", "(num/kg)/s", zt )

        k = k + 1

      case ('Nsnowm_dff')
        iNsnowm_dff = k
        call stat_assign( iNsnowm_dff, "Nsnowm_dff", & 
             "Nsnowm budget: Nsnowm diffusion [(num/kg)/s]", "(num/kg)/s", zt )

        k = k + 1

      case ('Nsnowm_mc')
        iNsnowm_mc = k
        call stat_assign( iNsnowm_mc, "Nsnowm_mc", & 
             "Nsnowm budget: Nsnowm microphysics [(num/kg)/s]", "(num/kg)/s", zt )

        k = k + 1

      case ('Nsnowm_cl')
        iNsnowm_cl = k

        call stat_assign( iNsnowm_cl, "Nsnowm_cl", & 
             "Nsnowm budget: Nsnowm clipping term [(num/kg)/s]", "(num/kg)/s", zt )
        k = k + 1

      case ('ricem_bt')
        iricem_bt = k

        call stat_assign( iricem_bt, "ricem_bt", & 
             "ricem budget: ricem time tendency [(kg/kg)/s]", "(kg/kg)/s", zt )

        k = k + 1

      case ('ricem_ma')
        iricem_ma = k

        call stat_assign( iricem_ma, "ricem_ma", & 
             "ricem budget: ricem vertical mean advection [(kg/kg)/s]", "(kg/kg)/s", zt )
        k = k + 1

      case ('ricem_sd')
        iricem_sd = k

        call stat_assign( iricem_sd, "ricem_sd", & 
             "ricem budget: ricem sedimentation [(kg/kg)/s]", "(kg/kg)/s", zt )
        k = k + 1

      case ('ricem_sd_mg_morr')
        iricem_sd_mg_morr = k

        call stat_assign( iricem_sd_mg_morr, "ricem_sd_mg_morr", & 
             "ricem sedimentation when using morrison or MG microphysics (not in budget," &
             // " included in ricem_mc) [(kg/kg)/s]", "(kg/kg)/s", zt )
        k = k + 1

      case ('ricem_dff')
        iricem_dff = k

        call stat_assign( iricem_dff, "ricem_dff", & 
             "ricem budget: ricem diffusion [(kg/kg)/s]", "(kg/kg)/s", zt )
        k = k + 1

      case ('ricem_mc')
        iricem_mc = k

        call stat_assign( iricem_mc, "ricem_mc", & 
             "ricem budget: Change in ricem due to microphysics [(kg/kg)/s]", "(kg/kg)/s", zt )
        k = k + 1

      case ('ricem_cl')
        iricem_cl = k

        call stat_assign( iricem_cl, "ricem_cl", & 
             "ricem budget: ricem clipping term [(kg/kg)/s]", "(kg/kg)/s", zt )
        k = k + 1

      case ('rgraupelm_bt')
        irgraupelm_bt = k

        call stat_assign( irgraupelm_bt, "rgraupelm_bt", & 
             "rgraupelm budget: rgraupelm time tendency [(kg/kg)/s]", "(kg/kg)/s", zt )
        k = k + 1

      case ('rgraupelm_ma')
        irgraupelm_ma = k

        call stat_assign( irgraupelm_ma, "rgraupelm_ma", & 
             "rgraupelm budget: rgraupelm vertical mean advection [(kg/kg)/s]", "(kg/kg)/s", zt )
        k = k + 1

      case ('rgraupelm_sd')
        irgraupelm_sd = k

        call stat_assign( irgraupelm_sd, "rgraupelm_sd", & 
             "rgraupelm budget: rgraupelm sedimentation [(kg/kg)/s]", "(kg/kg)/s", zt )
        k = k + 1

      case ('rgraupelm_sd_morr')
        irgraupelm_sd_morr = k

        call stat_assign( irgraupelm_sd_morr, "rgraupelm_sd_morr", & 
             "rgraupelm sedimentation when using morrison microphysics (not in budget, included" &
             // " in rgraupelm_mc) [(kg/kg)/s]", "(kg/kg)/s", zt )
        k = k + 1

      case ('rgraupelm_dff')
        irgraupelm_dff = k

        call stat_assign( irgraupelm_dff, "rgraupelm_dff", & 
             "rgraupelm budget: rgraupelm diffusion [(kg/kg)/s]", "(kg/kg)/s", zt )
        k = k + 1

      case ('rgraupelm_mc')
        irgraupelm_mc = k

        call stat_assign( irgraupelm_mc, "rgraupelm_mc", & 
             "rgraupelm budget: Change in rgraupelm due to microphysics [(kg/kg)/s]", &
             "(kg/kg)/s", zt )
        k = k + 1

      case ('rgraupelm_cl')
        irgraupelm_cl = k

        call stat_assign( irgraupelm_cl, "rgraupelm_cl", & 
             "rgraupelm budget: rgraupelm clipping term [(kg/kg)/s]", "(kg/kg)/s", zt )
        k = k + 1

      case ('Ngraupelm_bt')
        iNgraupelm_bt = k
        call stat_assign( iNgraupelm_bt, "Ngraupelm_bt", & 
             "Ngraupelm budget: [(num/kg)/s]", "(num/kg)/s", zt )

        k = k + 1

      case ('Ngraupelm_ma')
        iNgraupelm_ma = k

        call stat_assign( iNgraupelm_ma, "Ngraupelm_ma", & 
             "Ngraupelm budget: Ngraupelm mean advection [(num/kg)/s]", "(num/kg)/s", zt )
        k = k + 1

      case ('Ngraupelm_sd')
        iNgraupelm_sd = k

        call stat_assign( iNgraupelm_sd, "Ngraupelm_sd", & 
             "Ngraupelm budget: Ngraupelm sedimentation [(num/kg)/s]", "(num/kg)/s", zt )

        k = k + 1

      case ('Ngraupelm_dff')
        iNgraupelm_dff = k
        call stat_assign( iNgraupelm_dff, "Ngraupelm_dff", & 
             "Ngraupelm budget: Ngraupelm diffusion [(num/kg)/s]", "(num/kg)/s", zt )

        k = k + 1

      case ('Ngraupelm_mc')
        iNgraupelm_mc = k

        call stat_assign( iNgraupelm_mc, "Ngraupelm_mc", & 
             "Ngraupelm budget: Ngraupelm microphysics term [(num/kg)/s]", "(num/kg)/s", zt )
        k = k + 1

      case ('Ngraupelm_cl')
        iNgraupelm_cl = k

        call stat_assign( iNgraupelm_cl, "Ngraupelm_cl", & 
             "Ngraupelm budget: Ngraupelm clipping term [(num/kg)/s]", "(num/kg)/s", zt )
        k = k + 1

      case ('Nim_bt')
        iNim_bt = k
        call stat_assign( iNim_bt, "Nim_bt", & 
             "Nim budget: [(num/kg)/s]", "(num/kg)/s", zt )

        k = k + 1

      case ('Nim_ma')
        iNim_ma = k

        call stat_assign( iNim_ma, "Nim_ma", & 
             "Nim budget: Nim mean advection [(num/kg)/s]", "(num/kg)/s", zt )
        k = k + 1

      case ('Nim_sd')
        iNim_sd = k

        call stat_assign( iNim_sd, "Nim_sd", & 
             "Nim budget: Nim sedimentation [(num/kg)/s]", "(num/kg)/s", zt )

        k = k + 1

      case ('Nim_dff')
        iNim_dff = k
        call stat_assign( iNim_dff, "Nim_dff", & 
             "Nim budget: Nim diffusion [(num/kg)/s]", "(num/kg)/s", zt )

        k = k + 1

      case ('Nim_mc')
        iNim_mc = k

        call stat_assign( iNim_mc, "Nim_mc", & 
             "Nim budget: Nim microphysics term [(num/kg)/s]", "(num/kg)/s", zt )
        k = k + 1

      case ('Nim_cl')
        iNim_cl = k

        call stat_assign( iNim_cl, "Nim_cl", & 
             "Nim budget: Nim clipping term [(num/kg)/s]", "(num/kg)/s", zt )
        k = k + 1

      case ('Ncm_bt')
        iNcm_bt = k
        call stat_assign( iNcm_bt, "Ncm_bt", & 
             "Ncm budget: Cloud droplet number concentration budget [(num/kg)/s]", &
             "(num/kg)/s", zt )

        k = k + 1

      case ('Ncm_ma')
        iNcm_ma = k

        call stat_assign( iNcm_ma, "Ncm_ma", & 
             "Ncm budget: Ncm vertical mean advection [(num/kg)/s]", "(num/kg)/s", zt )
        k = k + 1

      case ('Ncm_act')
        iNcm_act = k

        call stat_assign( iNcm_act, "Ncm_act", &
             "Ncm budget: Change in Ncm due to activation [(num/kg)/s]", "(num/kg)/s", zt )

        k = k + 1

      case ('Ncm_dff')
        iNcm_dff = k
        call stat_assign( iNcm_dff, "Ncm_dff", & 
             "Ncm budget: Ncm diffusion [(num/kg)/s]", "(num/kg)/s", zt )

        k = k + 1

      case ('Ncm_mc')
        iNcm_mc = k

        call stat_assign( iNcm_mc, "Ncm_mc", & 
             "Ncm budget: Change in Ncm due to microphysics [(num/kg)/s]", "(num/kg)/s", zt )
        k = k + 1

      case ('Ncm_cl')
        iNcm_cl = k

        call stat_assign( iNcm_cl, "Ncm_cl", & 
             "Ncm budget: Ncm clipping term [(num/kg)/s]", "(num/kg)/s", zt )
        k = k + 1

      case ('vm_bt')
        ivm_bt = k

        call stat_assign( ivm_bt, "vm_bt", & 
             "vm budget: vm time tendency [m s^{-2}]", "m s^{-2}", zt )
        k = k + 1

      case ('vm_ma')
        ivm_ma = k
        call stat_assign( ivm_ma, "vm_ma", & 
             "vm budget: vm vertical mean advection [m s^{-2}]", "m s^{-2}", zt )
        k = k + 1

      case ('vm_gf')
        ivm_gf = k

        call stat_assign( ivm_gf, "vm_gf", & 
             "vm budget: vm geostrophic forcing [m s^{-2}]", "m s^{-2}", zt )
        k = k + 1

      case ('vm_cf')
        ivm_cf = k

        call stat_assign( ivm_cf, "vm_cf", & 
             "vm budget: vm coriolis forcing [m s^{-2}]", "m s^{-2}", zt )
        k = k + 1

      case ('vm_ta')
        ivm_ta = k

        call stat_assign( ivm_ta, "vm_ta", & 
             "vm budget: vm turbulent transport [m s^{-2}]", "m s^{-2}", zt )
        k = k + 1

      case ('vm_f')
        ivm_f = k
        call stat_assign( ivm_f, "vm_f", & 
             "vm budget: vm forcing [m s^{-2}]", "m s^{-2}", zt )
        k = k + 1

      case ('vm_sdmp')
        ivm_sdmp = k
        call stat_assign( ivm_sdmp, "vm_sdmp", & 
             "vm budget: vm sponge damping [m s^{-2}]", "m s^{-2}", zt )
        k = k + 1

      case ('vm_ndg')
        ivm_ndg = k
        call stat_assign( ivm_ndg, "vm_ndg", & 
             "vm budget: vm nudging [m s^{-2}]", "m s^{-2}", zt )
        k = k + 1

      case ('um_bt')
        ium_bt = k

        call stat_assign( ium_bt, "um_bt", & 
             "um budget: um time tendency [m s^{-2}]", "m s^{-2}", zt )
        k = k + 1

      case ('um_ma')
        ium_ma = k

        call stat_assign( ium_ma, "um_ma", & 
             "um budget: um vertical mean advection [m s^{-2}]", "m s^{-2}", zt )
        k = k + 1

      case ('um_gf')
        ium_gf = k
        call stat_assign( ium_gf, "um_gf", & 
             "um budget: um geostrophic forcing [m s^{-2}]", "m s^{-2}", zt )
        k = k + 1

      case ('um_cf')
        ium_cf = k
        call stat_assign( ium_cf, "um_cf", & 
             "um budget: um coriolis forcing [m s^{-2}]", "m s^{-2}", zt )
        k = k + 1

      case ('um_ta')
        ium_ta = k
        call stat_assign( ium_ta, "um_ta", & 
             "um budget: um turbulent advection [m s^{-2}]", "m s^{-2}", zt )
        k = k + 1

      case ('um_f')
        ium_f = k
        call stat_assign( ium_f, "um_f", & 
             "um budget: um forcing [m s^{-2}]", "m s^{-2}", zt )
        k = k + 1

      case ('um_sdmp')
        ium_sdmp = k
        call stat_assign( ium_sdmp, "um_sdmp", & 
             "um budget: um sponge damping [m s^{-2}]", "m s^{-2}", zt )
        k = k + 1

      case ('um_ndg')
        ium_ndg = k
        call stat_assign( ium_ndg, "um_ndg", & 
             "um budget: um nudging [m s^{-2}]", "m s^{-2}", zt )
        k = k + 1

      case ('mixt_frac')
        imixt_frac = k
        call stat_assign( imixt_frac, "mixt_frac", & 
             "pdf parameter: mixture fraction [count]", "count", zt )
        k = k + 1

      case ('w1')
        iw1 = k
        call stat_assign( iw1, "w1", & 
             "pdf parameter: mean w of component 1 [m/s]", "m/s", zt )

        k = k + 1

      case ('w2')
        iw2 = k

        call stat_assign( iw2, "w2", & 
             "pdf paramete: mean w of component 2 [m/s]", "m/s", zt )
        k = k + 1

      case ('varnce_w1')
        ivarnce_w1 = k
        call stat_assign( ivarnce_w1, "varnce_w1", & 
             "pdf parameter: w variance of component 1 [m^2/s^2]", "m^2/s^2", zt )

        k = k + 1

      case ('varnce_w2')
        ivarnce_w2 = k

        call stat_assign( ivarnce_w2, "varnce_w2", & 
             "pdf parameter: w variance of component 2 [m^2/s^2]", "m^2/s^2", zt )
        k = k + 1

      case ('thl1')
        ithl1 = k

        call stat_assign( ithl1, "thl1", & 
             "pdf parameter: mean thl of component 1 [K]", "K", zt )

        k = k + 1

      case ('thl2')
        ithl2 = k

        call stat_assign( ithl2, "thl2", & 
             "pdf parameter: mean thl of component 2 [K]", "K", zt )
        k = k + 1

      case ('varnce_thl1')
        ivarnce_thl1 = k

        call stat_assign( ivarnce_thl1, "varnce_thl1", & 
             "pdf parameter: thl variance of component 1 [K^2]", "K^2", zt )

        k = k + 1

      case ('varnce_thl2')
        ivarnce_thl2 = k
        call stat_assign( ivarnce_thl2, "varnce_thl2", & 
             "pdf parameter: thl variance of component 2 [K^2]", "K^2", zt )

        k = k + 1

      case ('rt1')
        irt1 = k
        call stat_assign( irt1, "rt1", & 
             "pdf parameter: mean rt of component 1 [kg/kg]", "kg/kg", zt )

        k = k + 1

      case ('rt2')
        irt2 = k

        call stat_assign( irt2, "rt2", & 
             "pdf parameter: mean rt of component 2 [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('varnce_rt1')
        ivarnce_rt1 = k
        call stat_assign( ivarnce_rt1, "varnce_rt1", & 
             "pdf parameter: rt variance of component 1 [(kg^2)/(kg^2)]", "(kg^2)/(kg^2)", zt )
        k = k + 1

      case ('varnce_rt2')
        ivarnce_rt2 = k

        call stat_assign( ivarnce_rt2, "varnce_rt2", & 
             "pdf parameter: rt variance of component 2 [(kg^2)/(kg^2)]", "(kg^2)/(kg^2)", zt )
        k = k + 1

      case ('rc1')
        irc1 = k

        call stat_assign( irc1, "rc1", & 
             "pdf parameter: mean rc of component 1 [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('rc2')
        irc2 = k

        call stat_assign( irc2, "rc2", & 
             "pdf parameter: mean rc of component 2 [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('rsl1')
        irsl1 = k

        call stat_assign( irsl1, "rsl1", & 
             "pdf parameter: sat mix rat based on tl1 [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('rsl2')
        irsl2 = k

        call stat_assign( irsl2, "rsl2", & 
             "pdf parameter: sat mix rat based on tl2 [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('cloud_frac1')
        icloud_frac1 = k
        call stat_assign( icloud_frac1, "cloud_frac1", & 
             "pdf parameter cloud_frac1 [count]", "count", zt )
        k = k + 1

      case ('cloud_frac2')
        icloud_frac2 = k

        call stat_assign( icloud_frac2, "cloud_frac2", & 
             "pdf parameter cloud_frac2 [count]", "count", zt )
        k = k + 1

      case ('s1')
        is1 = k

        call stat_assign( is1, "s1", & 
             "pdf parameter: Mellor's s (extended liq) for component 1 [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('s2')
        is2 = k

        call stat_assign( is2, "s2", & 
             "pdf parameter: Mellor's s (extended liq) for component 2 [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('stdev_s1')
        istdev_s1 = k

        call stat_assign( istdev_s1, "stdev_s1", & 
             "pdf parameter: Std dev of s1 [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('stdev_s2')
        istdev_s2 = k

        call stat_assign( istdev_s2, "stdev_s2", & 
             "pdf parameter: Std dev of s2 [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ('rrtthl')
        irrtthl = k

        call stat_assign( irrtthl, "rrtthl", & 
             "pdf parameter: Within-component correlation of rt and thl [count]", "count", zt )
        k = k + 1

      case('wp2_zt')
        iwp2_zt = k

        call stat_assign( iwp2_zt, "wp2_zt", & 
             "w'^2 interpolated to thermodyamic levels [m^2/s^2]", "m^2/s^2", zt )
        k = k + 1

      case('thlp2_zt')
        ithlp2_zt = k

        call stat_assign( ithlp2_zt, "thlp2_zt", & 
             "thl'^2 interpolated to thermodynamic levels [K^2]", "K^2", zt )
        k = k + 1

      case('wpthlp_zt')
        iwpthlp_zt = k

        call stat_assign( iwpthlp_zt, "wpthlp_zt", & 
             "w'thl' interpolated to thermodynamic levels [(m K)/s]", "(m K)/s", zt )
        k = k + 1

      case('wprtp_zt')
        iwprtp_zt = k

        call stat_assign( iwprtp_zt, "wprtp_zt", & 
             "w'rt' interpolated to thermodynamic levels [(m kg)/(s kg)]", "(m kg)/(s kg)", zt )
        k = k + 1

      case('rtp2_zt')
        irtp2_zt = k

        call stat_assign( irtp2_zt, "rtp2_zt", & 
             "rt'^2 interpolated to thermodynamic levels [kg/kg]", "kg/kg", zt )
        k = k + 1

      case('rtpthlp_zt')
        irtpthlp_zt = k

        call stat_assign( irtpthlp_zt, "rtpthlp_zt", & 
             "rt'thl' interpolated to thermodynamic levels [(kg K)/kg]", "(kg K)/kg", zt )
        k = k + 1

      case ('up2_zt')
        iup2_zt = k
        call stat_assign( iup2_zt, "up2_zt", & 
             "u'^2 interpolated to thermodynamic levels [m^2/s^2]", "m^2/s^2", zt )
        k = k + 1

      case ('vp2_zt')
        ivp2_zt = k
        call stat_assign( ivp2_zt, "vp2_zt", & 
             "v'^2 interpolated to thermodynamic levels [m^2/s^2]", "m^2/s^2", zt )
        k = k + 1

      case ('upwp_zt')
        iupwp_zt = k
        call stat_assign( iupwp_zt, "upwp_zt", & 
             "u'w' interpolated to thermodynamic levels [m^2/s^2]", "m^2/s^2", zt )
        k = k + 1

      case ('vpwp_zt')
        ivpwp_zt = k
        call stat_assign( ivpwp_zt, "vpwp_zt", & 
             "v'w' interpolated to thermodynamic levels [m^2/s^2]", "m^2/s^2", zt )
        k = k + 1

      case('LH_rvm_mc')
        iLH_rvm_mc = k

        call stat_assign( iLH_rvm_mc, "LH_rvm_mc", & 
             "Latin hypercube estimate of rvm_mc [kg/kg/s]", "kg/kg/s", zt )
        k = k + 1

      case('LH_thlm_mc')
        iLH_thlm_mc = k

        call stat_assign( iLH_thlm_mc, "LH_thlm_mc", & 
             "Latin hypercube estimate of thlm_mc [kg/kg/s]", "kg/kg/s", zt )
        k = k + 1

      case('LH_rcm_mc')
        iLH_rcm_mc = k

        call stat_assign( iLH_rcm_mc, "LH_rcm_mc", & 
             "Latin hypercube estimate of rcm_mc [kg/kg/s]", "kg/kg/s", zt )
        k = k + 1

      case('LH_Ncm_mc')
        iLH_Ncm_mc = k

        call stat_assign( iLH_Ncm_mc, "LH_Ncm_mc", & 
             "Latin hypercube estimate of Ncm_mc [kg/kg/s]", "kg/kg/s", zt )
        k = k + 1

      case('LH_rrainm_mc')
        iLH_rrainm_mc = k

        call stat_assign( iLH_rrainm_mc, "LH_rrainm_mc", & 
             "Latin hypercube estimate of rrainm_mc [kg/kg/s]", "kg/kg/s", zt )
        k = k + 1

      case('LH_Nrm_mc')
        iLH_Nrm_mc = k

        call stat_assign( iLH_Nrm_mc, "LH_Nrm_mc", & 
             "Latin hypercube estimate of Nrm_mc [kg/kg/s]", "kg/kg/s", zt )
        k = k + 1

      case('LH_rsnowm_mc')
        iLH_rsnowm_mc = k

        call stat_assign( iLH_rsnowm_mc, "LH_rsnowm_mc", & 
             "Latin hypercube estimate of rsnowm_mc [kg/kg/s]", "kg/kg/s", zt )
        k = k + 1

      case('LH_Nsnowm_mc')
        iLH_Nsnowm_mc = k

        call stat_assign( iLH_Nsnowm_mc, "LH_Nsnowm_mc", & 
             "Latin hypercube estimate of Nsnowm_mc [kg/kg/s]", "kg/kg/s", zt )
        k = k + 1

      case('LH_rgraupelm_mc')
        iLH_rgraupelm_mc = k

        call stat_assign( iLH_rgraupelm_mc, "LH_rgraupelm_mc", & 
             "Latin hypercube estimate of rgraupelm_mc [kg/kg/s]", "kg/kg/s", zt )
        k = k + 1

      case('LH_Ngraupelm_mc')
        iLH_Ngraupelm_mc = k

        call stat_assign( iLH_Ngraupelm_mc, "LH_Ngraupelm_mc", & 
             "Latin hypercube estimate of Ngraupelm_mc [kg/kg/s]", "kg/kg/s", zt )
        k = k + 1

      case('LH_ricem_mc')
        iLH_ricem_mc = k

        call stat_assign( iLH_ricem_mc, "LH_ricem_mc", & 
             "Latin hypercube estimate of ricem_mc [kg/kg/s]", "kg/kg/s", zt )
        k = k + 1

      case('LH_Nim_mc')
        iLH_Nim_mc = k

        call stat_assign( iLH_Nim_mc, "LH_Nim_mc", & 
             "Latin hypercube estimate of Nim_mc [kg/kg/s]", "kg/kg/s", zt )
        k = k + 1

      case ( 'LH_Vrr' )
        iLH_Vrr = k

        call stat_assign( iLH_Vrr, "LH_Vrr", & 
             "Latin hypercube estimate of rrainm sedimentation velocity [m/s]", "m/s", zt )
        k = k + 1

      case ( 'LH_VNr' )
        iLH_VNr = k

        call stat_assign( iLH_VNr, "LH_VNr", & 
             "Latin hypercube estimate of Nrm sedimentation velocity [m/s]", "m/s", zt )
        k = k + 1

      case ( 'LH_rcm_avg' )
        iLH_rcm_avg = k

        call stat_assign( iLH_rcm_avg, "LH_rcm_avg", & 
             "Latin hypercube average estimate of rcm [kg/kg]", "kg/kg", zt )

        k = k + 1

      case ( 'LH_rrainm' )
        iLH_rrainm = k

        call stat_assign( iLH_rrainm, "LH_rrainm", & 
             "Latin hypercube estimate of rrainm [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ( 'LH_Nrm' )
        iLH_Nrm = k

        call stat_assign( iLH_Nrm, "LH_Nrm", & 
             "Latin hypercube estimate of Nrm [count/kg]", "count/kg", zt )
        k = k + 1

      case ( 'LH_ricem' )
        iLH_ricem = k

        call stat_assign( iLH_ricem, "LH_ricem", & 
             "Latin hypercube estimate of ricem [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ( 'LH_Nim' )
        iLH_Nim = k

        call stat_assign( iLH_Nim, "LH_Nim", & 
             "Latin hypercube estimate of Nim [count/kg]", "count/kg", zt )
        k = k + 1

      case ( 'LH_rsnowm' )
        iLH_rsnowm = k

        call stat_assign( iLH_rsnowm, "LH_rsnowm", & 
             "Latin hypercube estimate of rsnowm [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ( 'LH_Nsnowm' )
        iLH_Nsnowm = k

        call stat_assign( iLH_Nsnowm, "LH_Nsnowm", & 
             "Latin hypercube estimate of Nsnowm [count/kg]", "count/kg", zt )
        k = k + 1


      case ( 'LH_rgraupelm' )
        iLH_rgraupelm = k

        call stat_assign( iLH_rgraupelm, "LH_rgraupelm", & 
             "Latin hypercube estimate of rgraupelm [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ( 'LH_Ngraupelm' )
        iLH_Ngraupelm = k

        call stat_assign( iLH_Ngraupelm, "LH_Ngraupelm", & 
             "Latin hypercube estimate of Ngraupelm [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ( 'LH_thlm' )
        iLH_thlm = k

        call stat_assign( iLH_thlm, "LH_thlm", & 
             "Latin hypercube estimate of thlm [K]", "K", zt )
        k = k + 1

      case ( 'LH_rcm' )
        iLH_rcm = k

        call stat_assign( iLH_rcm, "LH_rcm", & 
             "Latin hypercube estimate of rcm [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ( 'LH_Ncm' )
        iLH_Ncm = k

        call stat_assign( iLH_Ncm, "LH_Ncm", & 
             "Latin hypercube estimate of Ncm [count/kg]", "count/kg", zt )
        k = k + 1


      case ( 'LH_rvm' )
        iLH_rvm = k

        call stat_assign( iLH_rvm, "LH_rvm", & 
             "Latin hypercube estimate of rvm [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ( 'LH_wm' )
        iLH_wm = k

        call stat_assign( iLH_wm, "LH_wm", & 
             "Latin hypercube estimate of vertical velocity [m/s]", "m/s", zt )
        k = k + 1

      case ( 'LH_cloud_frac' )
        iLH_cloud_frac = k

        ! Note: count is the udunits compatible unit
        call stat_assign( iLH_cloud_frac, "LH_cloud_frac", & 
             "Latin hypercube estimate of cloud fraction [count]", "count", zt )
        k = k + 1

      case ( 'LH_wp2_zt' )
        iLH_wp2_zt = k
        call stat_assign( iLH_wp2_zt, "LH_wp2_zt", & 
             "Variance of the latin hypercube estimate of w [m^2/s^2]", "m^2/s^2", zt )
        k = k + 1

      case ( 'LH_Ncp2_zt' )
        iLH_Ncp2_zt = k
        call stat_assign( iLH_Ncp2_zt, "LH_Ncp2_zt", & 
             "Variance of the latin hypercube estimate of Nc [count^2/kg^2]", "count^2/kg^2", zt )
        k = k + 1

      case ( 'LH_Nrp2_zt' )
        iLH_Nrp2_zt = k
        call stat_assign( iLH_Nrp2_zt, "LH_Nrp2_zt", & 
             "Variance of the latin hypercube estimate of Nr [count^2/kg^2]", "count^2/kg^2", zt )
        k = k + 1

      case ( 'LH_rcp2_zt' )
        iLH_rcp2_zt = k
        call stat_assign( iLH_rcp2_zt, "LH_rcp2_zt", & 
             "Variance of the latin hypercube estimate of rc [kg^2/kg^2]", "kg^2/kg^2", zt )
        k = k + 1

      case ( 'LH_rtp2_zt' )
        iLH_rtp2_zt = k
        call stat_assign( iLH_rtp2_zt, "LH_rtp2_zt", & 
             "Variance of the latin hypercube estimate of rt [kg^2/kg^2]", "kg^2/kg^2", zt )
        k = k + 1

      case ( 'LH_thlp2_zt' )
        iLH_thlp2_zt = k
        call stat_assign( iLH_thlp2_zt, "LH_thlp2_zt", & 
             "Variance of the latin hypercube estimate of thl [K^2]", "K^2", zt )
        k = k + 1

      case ( 'LH_rrainp2_zt' )
        iLH_rrainp2_zt = k
        call stat_assign( iLH_rrainp2_zt, "LH_rrainp2_zt", & 
             "Variance of the latin hypercube estimate of rrain [kg^2/kg^2]", "kg^2/kg^2", zt )
        k = k + 1

      case ('C11_Skw_fnc')
        iC11_Skw_fnc = k

        call stat_assign( iC11_Skw_fnc, "C11_Skw_fnc", & 
             "C_11 parameter with Sk_w applied [-]", "count", zt )
        k = k + 1

      case ('s_mellor')
        is_mellor = k

        call stat_assign( is_mellor, "s_mellor", & 
             "Mellor's s (extended liq) [kg/kg]", "kg/kg", zt )
        k = k + 1

      case ( 'a3_coef_zt' )
        ia3_coef_zt = k
        call stat_assign( ia3_coef_zt, "a3_coef_zt", & 
             "The a3 coefficient interpolated the the zt grid [-]", "count", zt )
        k = k + 1

      case ( 'wp3_on_wp2_zt' )
        iwp3_on_wp2_zt = k
        call stat_assign( iwp3_on_wp2_zt, "wp3_on_wp2_zt", & 
             "Smoothed version of wp3 / wp2 [m/s]", "m/s", zt )
        k = k + 1

      case default

        l_found =.false.

        j=1

        do while( j <= sclr_dim .and. .not. l_found)
          write(sclr_idx, * ) j

          sclr_idx = adjustl(sclr_idx)

          if(trim(vars_zt(i)) == "sclr"//trim(sclr_idx)//"m" .and. .not. l_found) then

            isclrm(j) = k

            call stat_assign( isclrm(j) , "sclr"//trim(sclr_idx)//"m",&
              "passive scalar "//trim(sclr_idx), "unknown", zt )

            k = k + 1

            l_found = .true.

          else if(trim(vars_zt(i)) == "sclr"//trim(sclr_idx)//"m_f" .and. .not. l_found) then

            isclrm_f(j) = k

            call stat_assign( isclrm_f(j) , "sclr"//trim(sclr_idx)//"m_f", &
              "passive scalar forcing "//trim(sclr_idx), "unknown", zt )

            k = k + 1

            l_found = .true.

          endif

          j = j + 1
        end do

        j = 1

        do while( j <= edsclr_dim .and. .not. l_found)

          write(sclr_idx, * ) j

          sclr_idx = adjustl(sclr_idx)

          if(trim(vars_zt(i)) == "edsclr"//trim(sclr_idx)//"m" .and. .not. l_found ) then

            iedsclrm(j) = k

            call stat_assign( iedsclrm(j) , "edsclr"//trim(sclr_idx)//"m", &
               "passive scalar "//trim(sclr_idx), "unknown", zt )

            k = k + 1

            l_found = .true.

          else if(trim(vars_zt(i)) == "edsclr"//trim(sclr_idx)//"m_f" .and. .not. l_found) then

            iedsclrm_f(j) = k

            call stat_assign( iedsclrm_f(j) , "edsclr"//trim(sclr_idx)//"m_f", & 
              "passive scalar forcing "//trim(sclr_idx), "unknown", zt )

            k = k + 1

            l_found = .true.

          endif

          j = j + 1

        end do

        if (.not. l_found ) then

          write(fstderr,*) 'Error:  unrecognized variable in vars_zt:  ', trim( vars_zt(i) )

          l_error = .true.  ! This will stop the run.

        end if

      end select

    end do

    return
  end subroutine stats_init_zt

end module stats_zt
