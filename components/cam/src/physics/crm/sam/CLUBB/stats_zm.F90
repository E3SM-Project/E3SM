!-----------------------------------------------------------------------
! $Id: stats_zm.F90 6146 2013-04-05 18:02:22Z raut@uwm.edu $
module stats_zm

  implicit none

  private ! Default Scope

  public :: stats_init_zm

  ! Constant parameters
  integer, parameter, public :: nvarmax_zm = 250  ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_zm( vars_zm, l_error )

! Description:
!   Initializes array indices for zm

! Note:
!   All code that is within subroutine stats_init_zm, including variable
!   allocation code, is not called if l_stats is false.  This subroutine is
!   called only when l_stats is true.

!-----------------------------------------------------------------------

    use constants_clubb, only: &
        fstderr ! Constant(s)

    use stats_variables, only: & 
        zm, & 
        iwp2, & 
        irtp2, & 
        ithlp2, & 
        irtpthlp, & 
        iwprtp, & 
        iwpthlp, & 
        iwp3_zm, & 
        iwp4, & 
        iwpthvp, & 
        irtpthvp, & 
        ithlpthvp, & 
        itau_zm, & 
        iKh_zm, & 
        iwprcp, & 
        irc_coef, &
        ithlprcp, & 
        irtprcp, & 
        ircp2

    use stats_variables, only: &
        iupwp, & 
        ivpwp, & 
        irho_zm, & 
        isigma_sqd_w, &
        irho_ds_zm, &
        ithv_ds_zm, &
        iem, & 
        ishear, &
        imean_w_up, &
        imean_w_down, & 
        iFrad, & 
        iFrad_LW, & 
        iFrad_SW, & 
        iFrad_LW_up, & 
        iFrad_SW_up, & 
        iFrad_LW_down, & 
        iFrad_SW_down, & 
        iFprec, & 
        iFcsed

    use stats_variables, only: & 
        iup2, & 
        ivp2, & 
        iup2_bt, & 
        iup2_ta, & 
        iup2_tp, & 
        iup2_ma, & 
        iup2_dp1, & 
        iup2_dp2, & 
        iup2_pr1, & 
        iup2_pr2, & 
        iup2_cl, & 
        iup2_pd, &
        iup2_sf, &
        ivp2_bt, & 
        ivp2_ta, & 
        ivp2_tp, & 
        ivp2_ma, & 
        ivp2_dp1, & 
        ivp2_dp2, & 
        ivp2_pr1, & 
        ivp2_pr2, & 
        ivp2_cl, & 
        ivp2_pd, & 
        ivp2_sf

    use stats_variables, only: & 
        iVNr,  & 
        iVrr, &
        iVNc, & 
        iVrc, &
        iVNice, & 
        iVrice, &
        iVNsnow, &
        iVrsnow, &
        iVrgraupel, &
        iVrrprrp, &
        iVNrpNrp, &
        iVrrprrp_net, &
        iVNrpNrp_net

    use stats_variables, only: & 
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
        iwp2_cl, & 
        iwp2_pd, &
        iwp2_sf

    use stats_variables, only: & 
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
        iwprtp_pd, &
        iwprtp_forcing, & 
        iwprtp_mc, & 
        iwpthlp_bt, & 
        iwpthlp_ma, & 
        iwpthlp_ta

    use stats_variables, only: &
        iwpthlp_tp, & 
        iwpthlp_ac, & 
        iwpthlp_bp, & 
        iwpthlp_pr1, & 
        iwpthlp_pr2, & 
        iwpthlp_pr3, & 
        iwpthlp_dp1, & 
        iwpthlp_mfl, & 
        iwpthlp_cl, & 
        iwpthlp_sicl, &
        iwpthlp_forcing, &
        iwpthlp_mc

    use stats_variables, only: & 
        irtp2_bt, & 
        irtp2_ma, & 
        irtp2_ta, & 
        irtp2_tp, & 
        irtp2_dp1, & 
        irtp2_dp2, & 
        irtp2_cl, & 
        irtp2_pd, &
        irtp2_sf, &
        irtp2_forcing, &
        irtp2_mc, &
        ithlp2_bt, & 
        ithlp2_ma, & 
        ithlp2_ta, & 
        ithlp2_tp, & 
        ithlp2_dp1, & 
        ithlp2_dp2, & 
        ithlp2_cl, & 
        ithlp2_pd

    use stats_variables, only: &
        ithlp2_sf, &
        ithlp2_forcing, &
        ithlp2_mc, &
        irtpthlp_bt, & 
        irtpthlp_ma, & 
        irtpthlp_ta, & 
        irtpthlp_tp1, & 
        irtpthlp_tp2, & 
        irtpthlp_dp1, & 
        irtpthlp_dp2, & 
        irtpthlp_cl, &
        irtpthlp_sf, &
        irtpthlp_forcing, &
        irtpthlp_mc
    
    use stats_variables, only: & 
        iwpthlp_entermfl, & ! Variable(s)
        iwpthlp_exit_mfl, &
        iwpthlp_mfl_min, &
        iwpthlp_mfl_max, &
        iwprtp_enter_mfl, &
        iwprtp_exit_mfl, &
        iwprtp_mfl_min, &
        iwprtp_mfl_max

    use stats_variables, only: & 
        iwm_zm,  &  ! Variable
        icloud_frac_zm, &
        iice_supersat_frac_zm, &
        ircm_zm, &
        irtm_zm, &
        ithlm_zm

    use stats_variables, only: & 
        isclrprtp, & 
        isclrp2, & 
        isclrpthvp, & 
        isclrpthlp, & 
        isclrprcp, & 
        iwpsclrp, & 
        iwp2sclrp, & 
        iwpsclrp2, & 
        iwpsclrprtp, & 
        iwpsclrpthlp, & 
        iwpedsclrp

    use stats_variables, only: &
        ia3_coef, &
        iwp3_on_wp2, &
        iSkw_velocity, &
        igamma_Skw_fnc, &
        iC6rt_Skw_fnc, &
        iC6thl_Skw_fnc, &
        iC7_Skw_fnc, &
        iC1_Skw_fnc

    use stats_type, only: & 
        stat_assign ! Procedure

    use parameters_model, only: &
        sclr_dim, &
        edsclr_dim

!   use error_code, only: &
!       clubb_at_least_debug_level ! Function

    implicit none

    ! Input Variable
    ! zm variable names

    character(len= * ), dimension(nvarmax_zm), intent(in) :: vars_zm

    ! Output Variable
    logical, intent(inout) :: l_error

    ! Local Varables
    integer :: i,j, k

    logical :: l_found

    character(len=50) :: sclr_idx

!     Default initialization for array indices for zm

    iwp2          = 0
    irtp2         = 0
    ithlp2        = 0
    irtpthlp      = 0
    iwprtp        = 0
    iwpthlp       = 0
    iwp3_zm       = 0
    iwp4          = 0
    iwpthvp       = 0
    irtpthvp      = 0
    ithlpthvp     = 0
    itau_zm       = 0
    iKh_zm        = 0
    iwprcp        = 0
    irc_coef      = 0
    ithlprcp      = 0
    irtprcp       = 0
    ircp2         = 0
    iupwp         = 0
    ivpwp         = 0
    irho_zm       = 0
    isigma_sqd_w  = 0
    irho_ds_zm    = 0
    ithv_ds_zm    = 0
    iem           = 0
    ishear        = 0  ! Brian
    imean_w_up    = 0
    imean_w_down  = 0
    iFrad         = 0
    iFrad_LW      = 0  ! Brian
    iFrad_SW      = 0  ! Brian
    iFrad_LW_up   = 0  ! Brian
    iFrad_SW_up   = 0  ! Brian
    iFrad_LW_down = 0  ! Brian
    iFrad_SW_down = 0  ! Brian
    iFprec        = 0  ! Brian
    iFcsed        = 0  ! Brian


    iup2 = 0
    ivp2 = 0

    iup2_bt  = 0
    iup2_ta  = 0
    iup2_tp  = 0
    iup2_ma  = 0
    iup2_dp1 = 0
    iup2_dp2 = 0
    iup2_pr1 = 0
    iup2_pr2 = 0
    iup2_cl  = 0
    iup2_sf  = 0

    ivp2_bt  = 0
    ivp2_ta  = 0
    ivp2_tp  = 0
    ivp2_ma  = 0
    ivp2_dp1 = 0
    ivp2_dp2 = 0
    ivp2_pr1 = 0
    ivp2_pr2 = 0
    ivp2_cl  = 0
    ivp2_sf  = 0

    ! Sedimentation velocities
    iVNr       = 0
    iVrr       = 0
    iVNc       = 0
    iVrc       = 0
    iVNice     = 0
    iVrice     = 0
    iVrgraupel = 0
    iVNsnow    = 0
    iVrsnow    = 0

    ! Covariance of sedimentation velocity and hydrometeor, <V_xx'x_x'>
    iVrrprrp     = 0
    iVNrpNrp     = 0
    iVrrprrp_net = 0
    iVNrpNrp_net = 0

    ! Vertical velocity budgets
    iwp2_bt   = 0
    iwp2_ma   = 0
    iwp2_ta   = 0
    iwp2_ac   = 0
    iwp2_bp   = 0
    iwp2_pr1  = 0
    iwp2_pr2  = 0
    iwp2_pr3  = 0
    iwp2_dp1  = 0
    iwp2_dp2  = 0
    iwp2_4hd  = 0
    iwp2_cl   = 0
    iwp2_pd   = 0
    iwp2_sf   = 0

    ! Flux budgets
    iwprtp_bt      = 0
    iwprtp_ma      = 0
    iwprtp_ta      = 0
    iwprtp_tp      = 0
    iwprtp_ac      = 0
    iwprtp_bp      = 0
    iwprtp_pr1     = 0
    iwprtp_pr2     = 0
    iwprtp_pr3     = 0
    iwprtp_dp1     = 0
    iwprtp_mfl     = 0
    iwprtp_cl      = 0
    iwprtp_sicl    = 0
    iwprtp_pd      = 0
    iwprtp_forcing = 0
    iwprtp_mc      = 0

    iwpthlp_bt      = 0
    iwpthlp_ma      = 0
    iwpthlp_ta      = 0
    iwpthlp_tp      = 0
    iwpthlp_ac      = 0
    iwpthlp_bp      = 0
    iwpthlp_pr1     = 0
    iwpthlp_pr2     = 0
    iwpthlp_pr3     = 0
    iwpthlp_dp1     = 0
    iwpthlp_mfl     = 0
    iwpthlp_cl      = 0
    iwpthlp_sicl    = 0
    iwpthlp_forcing = 0
    iwpthlp_mc      = 0

    ! Variance budgets
    irtp2_bt      = 0
    irtp2_ma      = 0
    irtp2_ta      = 0
    irtp2_tp      = 0
    irtp2_dp1     = 0
    irtp2_dp2     = 0
    irtp2_cl      = 0
    irtp2_pd      = 0
    irtp2_sf      = 0
    irtp2_forcing = 0
    irtp2_mc      = 0

    ithlp2_bt      = 0
    ithlp2_ma      = 0
    ithlp2_ta      = 0
    ithlp2_tp      = 0
    ithlp2_dp1     = 0
    ithlp2_dp2     = 0
    ithlp2_cl      = 0
    ithlp2_pd      = 0
    ithlp2_sf      = 0
    ithlp2_forcing = 0
    ithlp2_mc      = 0

    irtpthlp_bt      = 0
    irtpthlp_ma      = 0
    irtpthlp_ta      = 0
    irtpthlp_tp1     = 0
    irtpthlp_tp2     = 0
    irtpthlp_dp1     = 0
    irtpthlp_dp2     = 0
    irtpthlp_cl      = 0
    irtpthlp_sf      = 0
    irtpthlp_forcing = 0
    irtpthlp_mc      = 0

    !Monatonic flux limiter diagnostic output
    iwpthlp_mfl_min = 0
    iwpthlp_mfl_max = 0
    iwpthlp_entermfl = 0
    iwpthlp_exit_mfl = 0
    iwprtp_mfl_min = 0
    iwprtp_mfl_max = 0
    iwprtp_enter_mfl = 0
    iwprtp_exit_mfl = 0

    ! Skewness velocity
    iSkw_velocity = 0

    ! Skewness function
    igamma_Skw_fnc = 0
    iC6rt_Skw_fnc = 0
    iC6thl_Skw_fnc = 0
    iC7_Skw_fnc = 0
    iC1_Skw_fnc = 0

    ia3_coef = 0
    iwp3_on_wp2 = 0

    allocate(isclrprtp(1:sclr_dim))
    allocate(isclrp2(1:sclr_dim))
    allocate(isclrpthvp(1:sclr_dim))
    allocate(isclrpthlp(1:sclr_dim))
    allocate(isclrprcp(1:sclr_dim))
    allocate(iwpsclrp(1:sclr_dim))
    allocate(iwp2sclrp(1:sclr_dim))
    allocate(iwpsclrp2(1:sclr_dim))
    allocate(iwpsclrprtp(1:sclr_dim))
    allocate(iwpsclrpthlp(1:sclr_dim))

    allocate(iwpedsclrp(1:edsclr_dim))

!     Assign pointers for statistics variables zm

    isclrprtp    = 0
    isclrp2      = 0
    isclrpthvp   = 0
    isclrpthlp   = 0
    isclrprcp    = 0
    iwpsclrp     = 0
    iwp2sclrp    = 0
    iwpsclrp2    = 0
    iwpsclrprtp  = 0
    iwpsclrpthlp = 0

    iwpedsclrp   = 0

!     Assign pointers for statistics variables zm

    k = 1
    do i=1,zm%nn

      select case ( trim(vars_zm(i)) )

      case ('wp2')
        iwp2 = k
        call stat_assign(iwp2,"wp2", & 
             "w'^2, Variance of vertical air velocity [m^2/s^2]","m^2/s^2",zm)
        k = k + 1

      case ('rtp2')
        irtp2 = k
        call stat_assign(irtp2,"rtp2", & 
             "rt'^2, Variance of rt [(kg/kg)^2]","(kg/kg)^2",zm)
        k = k + 1

      case ('thlp2')
        ithlp2 = k
        call stat_assign(ithlp2,"thlp2", & 
             "thl'^2, Variance of thl [K^2]","K^2",zm)
        k = k + 1

      case ('rtpthlp')
        irtpthlp = k
        call stat_assign(irtpthlp,"rtpthlp", & 
             "rt'thl', Covariance of rt and thl [(kg K)/kg]","(kg K)/kg",zm)
        k = k + 1

      case ('wprtp')
        iwprtp = k

        call stat_assign(iwprtp,"wprtp", & 
             "w'rt', Vertical turbulent flux of rt [(kg/kg) m/s]","(m kg)/(s kg)",zm)
        k = k + 1

      case ('wpthlp')
        iwpthlp = k

        call stat_assign(iwpthlp,"wpthlp", & 
             "w'thl', Vertical turbulent flux of thl [K m/s]","(m K)/s",zm)
        k = k + 1

      case ('wp3_zm')
        iwp3_zm = k
        call stat_assign( iwp3_zm, "wp3_zm", & 
             "w'^3 interpolated to moment. levels [m^3/s^3]", "(m^3)/(s^3)", zm )
        k = k + 1

      case ('wp4')
        iwp4 = k
        call stat_assign(iwp4,"wp4", & 
             "w'^4 [m^4/s^4]","(m^4)/(s^4)",zm)
        k = k + 1

      case ('wpthvp')
        iwpthvp = k
        call stat_assign(iwpthvp,"wpthvp", & 
             "Buoyancy flux [K m/s]","K m/s",zm)
        k = k + 1

      case ('rtpthvp')
        irtpthvp = k
        call stat_assign(irtpthvp,"rtpthvp", & 
             "rt'thv' [(kg/kg) K]","(kg/kg) K",zm)
        k = k + 1

      case ('thlpthvp')
        ithlpthvp = k
        call stat_assign(ithlpthvp,"thlpthvp", & 
             "thl'thv' [K^2]","K^2",zm)
        k = k + 1

      case ('tau_zm')
        itau_zm = k

        call stat_assign(itau_zm,"tau_zm", & 
             "Time-scale tau on momentum levels [s]","s",zm)
        k = k + 1

      case ('Kh_zm')
        iKh_zm = k

        call stat_assign(iKh_zm,"Kh_zm", & 
             "Eddy diffusivity on momentum levels [m^2/s]","m^2/s",zm)
        k = k + 1

      case ('wprcp')
        iwprcp = k
        call stat_assign(iwprcp,"wprcp", & 
             "w' rc' [(m/s) (kg/kg)]","(m/s) (kg/kg)",zm)
        k = k + 1
      
      case ('rc_coef')
        irc_coef = k
        call stat_assign(irc_coef, "rc_coef", &
            "Coefficient of X' R_l' in Eq. (34)", "[-]", zm)
        k = k + 1

      case ('thlprcp')
        ithlprcp = k
        call stat_assign(ithlprcp,"thlprcp", & 
             "thl' rc' [K (kg/kg)]","K (kg/kg)",zm)
        k = k + 1

      case ('rtprcp')
        irtprcp = k

        call stat_assign(irtprcp,"rtprcp", & 
             "rt'rc' [(kg^2)/(kg^2)]","(kg^2)/(kg^2)",zm)
        k = k + 1

      case ('rcp2')
        ircp2 = k
        call stat_assign(ircp2,"rcp2", & 
             "rc'^2 [(kg^2)/(kg^2)]","(kg^2)/(kg^2)",zm)
        k = k + 1
      case ('upwp')
        iupwp = k
        call stat_assign(iupwp,"upwp", & 
             "u'w', Vertical east-west momentum flux [m^2/s^2]","m^2/s^2",zm)
        k = k + 1
      case ('vpwp')
        ivpwp = k
        call stat_assign(ivpwp,"vpwp", & 
             "v'w', Vertical north-south momentum flux [m^2/s^2]","m^2/s^2",zm)
        k = k + 1
      case ('rho_zm')
        irho_zm = k
        call stat_assign(irho_zm,"rho_zm", & 
             "Density on momentum levels [kg/m^3]","kg m^{-3}",zm)
        k = k + 1
      case ('sigma_sqd_w')
        isigma_sqd_w = k
        call stat_assign(isigma_sqd_w,"sigma_sqd_w", & 
             "Nondimensionalized w variance of Gaussian component [-]","-",zm)
        k = k + 1
      case ('rho_ds_zm')
        irho_ds_zm = k
        call stat_assign(irho_ds_zm,"rho_ds_zm", &
             "Dry, static, base-state density [kg/m^3]","kg m^{-3}",zm)
        k = k + 1
      case ('thv_ds_zm')
        ithv_ds_zm = k
        call stat_assign(ithv_ds_zm,"thv_ds_zm", &
             "Dry, base-state theta_v [K]","K",zm)
        k = k + 1
      case ('em')
        iem = k
        call stat_assign(iem,"em", & 
             "Turbulent kinetic energy, usu. 0.5*(u'^2+v'^2+w'^2) [m^2/s^2]","m^2/s^2",zm)
        k = k + 1
      case ('shear')      ! Brian
        ishear = k
        call stat_assign(ishear,"shear", & 
             "Wind shear production term [m^2/s^3]","m^2/s^3",zm)
        k = k + 1
      case ('mean_w_up')
        imean_w_up = k
        call stat_assign(imean_w_up, "mean_w_up", & 
             "Mean w >= w_ref [m/s]", "m/s", zm)
        k = k + 1
      case ('mean_w_down')
        imean_w_down = k
        call stat_assign(imean_w_down, "mean_w_down", & 
             "Mean w <= w_ref [m/s]", "m/s", zm)
        k = k + 1
      case ('Frad')
        iFrad = k
        call stat_assign(iFrad,"Frad", & 
             "Total (sw+lw) net (up+down) radiative flux [W/m^2]","W/m^2",zm)
        k = k + 1
      case ('Frad_LW')    ! Brian
        iFrad_LW = k
        call stat_assign(iFrad_LW,"Frad_LW", & 
             "Net long-wave radiative flux [W/m^2]","W/m^2",zm)
        k = k + 1
      case ('Frad_SW')    ! Brian
        iFrad_SW = k

        call stat_assign(iFrad_SW,"Frad_SW", & 
             "Net short-wave radiative flux [W/m^2]","W/m^2",zm)
        k = k + 1

      case ('Frad_LW_up')    ! Brian
        iFrad_LW_up = k
        call stat_assign(iFrad_LW_up,"Frad_LW_up", & 
             "Long-wave upwelling radiative flux [W/m^2]","W/m^2",zm)
        k = k + 1
      case ('Frad_SW_up')    ! Brian
        iFrad_SW_up = k

        call stat_assign(iFrad_SW_up,"Frad_SW_up", & 
             "Short-wave upwelling radiative flux [W/m^2]","W/m^2",zm)
        k = k + 1

      case ('Frad_LW_down')    ! Brian
        iFrad_LW_down = k
        call stat_assign(iFrad_LW_down,"Frad_LW_down", & 
        "Long-wave downwelling radiative flux [W/m^2]", "W/m^2", zm )
        k = k + 1
      case ('Frad_SW_down')    ! Brian
        iFrad_SW_down = k

        call stat_assign(iFrad_SW_down,"Frad_SW_down", & 
        "Short-wave downwelling radiative flux [W/m^2]", "W/m^2", zm )
        k = k + 1


      case ('Fprec')      ! Brian
        iFprec = k

        call stat_assign(iFprec,"Fprec", & 
             "Rain flux [W/m^2]","W/m^2",zm)
        k = k + 1

      case ('Fcsed')      ! Brian
        iFcsed = k

        call stat_assign(iFcsed,"Fcsed", & 
             "cloud water sedimentation flux [kg/(s*m^2)]", & 
             "kg/(s*m^2)",zm)
        k = k + 1

      case ('VNr')
        iVNr = k

        call stat_assign(iVNr,"VNr", & 
             "rrainm concentration fallspeed [m/s]","m/s",zm)
        k = k + 1

      case ('Vrr')
        iVrr = k

        call stat_assign(iVrr,"Vrr", & 
             "rrainm mixing ratio fallspeed [m/s]","m/s",zm)
        k = k + 1

      case ('VNc')
        iVNc = k

        call stat_assign(iVNc,"VNc", & 
             "Nrm concentration fallspeed [m/s]","m/s",zm)
        k = k + 1

      case ('Vrc')
        iVrc = k

        call stat_assign(iVrc,"Vrc", & 
             "Nrm mixing ratio fallspeed [m/s]","m/s",zm)
        k = k + 1

      case ('VNsnow')
        iVNsnow = k

        call stat_assign(iVNsnow,"VNsnow", & 
             "Snow concentration fallspeed [m/s]","m/s",zm)
        k = k + 1

      case ('Vrsnow')
        iVrsnow = k

        call stat_assign(iVrsnow,"Vrsnow", & 
             "Snow mixing ratio fallspeed [m/s]","m/s",zm)
        k = k + 1

      case ('Vrgraupel')
        iVrgraupel = k

        call stat_assign(iVrgraupel,"Vrgraupel", & 
             "Graupel sedimentation velocity [m/s]","m/s",zm)
        k = k + 1

      case ('VNice')
        iVNice = k

        call stat_assign(iVNice,"VNice", & 
             "Cloud ice concentration fallspeed [m/s]","m/s",zm)
        k = k + 1

      case ('Vrice')
        iVrice = k

        call stat_assign(iVrice,"Vrice", & 
             "Cloud ice mixing ratio fallspeed [m/s]","m/s",zm)
        k = k + 1

      case ('Vrrprrp')
        iVrrprrp = k

        call stat_assign( iVrrprrp, "Vrrprrp", & 
             "Covariance of V_rr (r_r sed. vel.) and r_r [(m/s)(kg/kg)]", &
             "(m/s)(kg/kg)", zm )
        k = k + 1

      case ('VNrpNrp')
        iVNrpNrp = k

        call stat_assign( iVNrpNrp, "VNrpNrp", & 
             "Covariance of V_Nr (N_r sed. vel.) and N_r [(m/s)(num/kg)]", &
             "(m/s)(num/kg)", zm )
        k = k + 1

      case ('Vrrprrp_net')
        iVrrprrp_net = k

        call stat_assign( iVrrprrp_net, "Vrrprrp_net", & 
             "Adjusted value of < V_rr'r_r' > (turb. sed. flux limiter)" &
             //" [(m/s)(kg/kg)]", "(m/s)(kg/kg)", zm )
        k = k + 1

      case ('VNrpNrp_net')
        iVNrpNrp_net = k

        call stat_assign( iVNrpNrp_net, "VNrpNrp_net", & 
             "Adjusted value of < V_Nr'N_r' > (turb. sed. flux limiter)" &
             //" [(m/s)(num/kg)]", "(m/s)(num/kg)", zm )
        k = k + 1

      case ('wp2_bt')
        iwp2_bt = k

        call stat_assign(iwp2_bt,"wp2_bt", & 
             "wp2 budget: wp2 time tendency [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('wp2_ma')
        iwp2_ma = k

        call stat_assign(iwp2_ma,"wp2_ma", & 
             "wp2 budget: wp2 vertical mean advection [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('wp2_ta')
        iwp2_ta = k

        call stat_assign(iwp2_ta,"wp2_ta", & 
             "wp2 budget: wp2 turbulent advection [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('wp2_ac')
        iwp2_ac = k

        call stat_assign(iwp2_ac,"wp2_ac", & 
             "wp2 budget: wp2 accumulation term [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('wp2_bp')
        iwp2_bp = k

        call stat_assign(iwp2_bp,"wp2_bp", & 
             "wp2 budget: wp2 buoyancy production [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('wp2_pr1')
        iwp2_pr1 = k

        call stat_assign(iwp2_pr1,"wp2_pr1", & 
             "wp2 budget: wp2 pressure term 1 [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('wp2_pr2')
        iwp2_pr2 = k
        call stat_assign(iwp2_pr2,"wp2_pr2", & 
             "wp2 budget: wp2 pressure term 2 [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('wp2_pr3')
        iwp2_pr3 = k
        call stat_assign(iwp2_pr3,"wp2_pr3", & 
             "wp2 budget: wp2 pressure term 3 [m^2/s^3]","m^2/s^3",zm)

        k = k + 1

      case ('wp2_dp1')
        iwp2_dp1 = k
        call stat_assign(iwp2_dp1,"wp2_dp1", & 
             "wp2 budget: wp2 dissipation term 1 [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('wp2_dp2')
        iwp2_dp2 = k
        call stat_assign(iwp2_dp2,"wp2_dp2", & 
             "wp2 budget: wp2 dissipation term 2 [m^2/s^3]","m^2/s^3",zm)

        k = k + 1

      case ('wp2_4hd')
        iwp2_4hd = k
        call stat_assign(iwp2_4hd,"wp2_4hd", & 
             "wp2 budget: wp2 4th-order hyper-diffusion [m^2/s^3]","m^2/s^3",zm)

        k = k + 1

      case ('wp2_cl')
        iwp2_cl = k

        call stat_assign(iwp2_cl,"wp2_cl", & 
             "wp2 budget: wp2 clipping term [m^2/s^3]","m^2/s^3",zm)

        k = k + 1

      case ('wp2_pd')
        iwp2_pd = k

        call stat_assign(iwp2_pd,"wp2_pd", & 
             "wp2 budget: wp2 positive definite adjustment [m^2/s^3]","m2/s3",zm)

        k = k + 1
        
      case ('wp2_sf')
        iwp2_sf = k
        
        call stat_assign( iwp2_sf, "wp2_sf", & 
             "wp2 budget: wp2 surface variance [m^2/s^3]","m2/s3",zm)
             
        k = k + 1

      case ('wprtp_bt')
        iwprtp_bt = k
        call stat_assign(iwprtp_bt,"wprtp_bt", & 
             "wprtp budget: wprtp time tendency [(m kg)/(s^2 kg)]","(m kg)/(s^2 kg)",zm)
        k = k + 1

      case ('wprtp_ma')
        iwprtp_ma = k

        call stat_assign(iwprtp_ma,"wprtp_ma", & 
             "wprtp budget: wprtp mean advection [(m kg)/(s^2 kg)]","(m kg)/(s^2 kg)",zm)
        k = k + 1

      case ('wprtp_ta')
        iwprtp_ta = k

        call stat_assign(iwprtp_ta,"wprtp_ta", & 
             "wprtp budget: wprtp turbulent advection [(m kg)/(s^2 kg)]","(m kg)/(s^2 kg)",zm)
        k = k + 1

      case ('wprtp_tp')
        iwprtp_tp = k

        call stat_assign(iwprtp_tp,"wprtp_tp", & 
             "wprtp budget: wprtp turbulent production [(m kg)/(s^2 kg)]","(m kg)/(s^2 kg)",zm)
        k = k + 1

      case ('wprtp_ac')
        iwprtp_ac = k

        call stat_assign(iwprtp_ac,"wprtp_ac", & 
             "wprtp budget: wprtp accumulation term [(m kg)/(s^2 kg)]","(m kg)/(s^2 kg)",zm)
        k = k + 1

      case ('wprtp_bp')
        iwprtp_bp = k

        call stat_assign(iwprtp_bp,"wprtp_bp", & 
             "wprtp budget: wprtp buoyancy production [(m kg)/(s^2 kg)]","(m kg)/(s^2 kg)",zm)
        k = k + 1

      case ('wprtp_pr1')
        iwprtp_pr1 = k

        call stat_assign(iwprtp_pr1,"wprtp_pr1", & 
             "wprtp budget: wprtp pressure term 1 [(m kg)/(s^2 kg)]","(m kg)/(s^2 kg)",zm)
        k = k + 1

      case ('wprtp_pr2')
        iwprtp_pr2 = k

        call stat_assign(iwprtp_pr2,"wprtp_pr2", & 
             "wprtp budget: wprtp pressure term 2 [(m kg)/(s^2 kg)]","(m kg)/(s^2 kg)",zm)
        k = k + 1

      case ('wprtp_pr3')
        iwprtp_pr3 = k

        call stat_assign(iwprtp_pr3,"wprtp_pr3", & 
             "wprtp budget: wprtp pressure term 3 [(m kg)/(s^2 kg)]","(m kg)/(s^2 kg)",zm)
        k = k + 1

      case ('wprtp_dp1')
        iwprtp_dp1 = k

        call stat_assign(iwprtp_dp1,"wprtp_dp1", & 
             "wprtp budget: wprtp dissipation term 1 [(m kg)/(s^2 kg)]","(m kg)/(s^2 kg)",zm)
        k = k + 1

      case ('wprtp_mfl')
        iwprtp_mfl = k

        call stat_assign(iwprtp_mfl,"wprtp_mfl", & 
             "wprtp budget: wprtp monotonic flux limiter [(m kg)/(s^2 kg)]","(m kg)/(s^2 kg)",zm)
        k = k + 1

      case ('wprtp_cl')
        iwprtp_cl = k

        call stat_assign(iwprtp_cl,"wprtp_cl", & 
             "wprtp budget: wprtp clipping term [(m kg)/(s^2 kg)]","(m kg)/(s^2 kg)",zm)
        k = k + 1

      case ('wprtp_sicl')
        iwprtp_sicl = k

        call stat_assign(iwprtp_sicl,"wprtp_sicl", & 
             "wprtp budget: wprtp semi-implicit clipping term [(m kg)/(s^2 kg)]", &
             "(m kg)/(s^2 kg)",zm)
        k = k + 1

      case ('wprtp_pd')
        iwprtp_pd = k

        call stat_assign(iwprtp_pd,"wprtp_pd", & 
             "wprtp budget: wprtp flux corrected trans. term [(m kg)/(s^2 kg)]", &
             "(m kg)/(s^2 kg)",zm)
        k = k + 1

      case ('wprtp_forcing')
        iwprtp_forcing = k

        call stat_assign( iwprtp_forcing, "wprtp_forcing", & 
             "wprtp budget: wprtp forcing (includes microphysics tendency) [(m kg/kg)/s^2]", &
             "(m kg/kg)/s^2", zm )
        k = k + 1

      case ('wprtp_mc')
        iwprtp_mc = k

        call stat_assign( iwprtp_mc, "wprtp_mc", & 
             "Microphysics tendency for wprtp (not in budget) [(m kg/kg)/s^2]", &
             "(m kg/kg)/s^2", zm )
        k = k + 1

      case ('wpthlp_bt')
        iwpthlp_bt = k

        call stat_assign(iwpthlp_bt,"wpthlp_bt", & 
             "wpthlp budget: [(m K)/s^2]","(m K)/s^2",zm)
        k = k + 1

      case ('wpthlp_ma')
        iwpthlp_ma = k
        call stat_assign(iwpthlp_ma,"wpthlp_ma", & 
             "wpthlp budget: wpthlp mean advection [(m K)/s^2]","(m K)/s^2",zm)

        k = k + 1

      case ('wpthlp_ta')
        iwpthlp_ta = k
        call stat_assign(iwpthlp_ta,"wpthlp_ta", & 
             "wpthlp budget: wpthlp turbulent advection [(m K)/s^2]","(m K)/s^2",zm)

        k = k + 1

      case ('wpthlp_tp')
        iwpthlp_tp = k
        call stat_assign(iwpthlp_tp,"wpthlp_tp", & 
             "wpthlp budget: wpthlp turbulent production [(m K)/s^2]","(m K)/s^2",zm)

        k = k + 1

      case ('wpthlp_ac')
        iwpthlp_ac = k
        call stat_assign(iwpthlp_ac,"wpthlp_ac", & 
             "wpthlp budget: wpthlp accumulation term [(m K)/s^2]","(m K)/s^2",zm)

        k = k + 1

      case ('wpthlp_bp')
        iwpthlp_bp = k
        call stat_assign(iwpthlp_bp,"wpthlp_bp", & 
             "wpthlp budget: wpthlp buoyancy production [(m K)/s^2]","(m K)/s^2",zm)
        k = k + 1

      case ('wpthlp_pr1')
        iwpthlp_pr1 = k

        call stat_assign(iwpthlp_pr1,"wpthlp_pr1", & 
             "wpthlp budget: wpthlp pressure term 1 [(m K)/s^2]","(m K)/s^2",zm)
        k = k + 1

      case ('wpthlp_pr2')
        iwpthlp_pr2 = k

        call stat_assign(iwpthlp_pr2,"wpthlp_pr2", & 
             "wpthlp budget: wpthlp pressure term 2 [(m K)/s^2]","(m K)/s^2",zm)
        k = k + 1

      case ('wpthlp_pr3')
        iwpthlp_pr3 = k
        call stat_assign(iwpthlp_pr3,"wpthlp_pr3", & 
             "wpthlp budget: wpthlp pressure term 3 [(m K)/s^2]","(m K)/s^2",zm)
        k = k + 1

      case ('wpthlp_dp1')
        iwpthlp_dp1 = k
        call stat_assign(iwpthlp_dp1,"wpthlp_dp1", & 
             "wpthlp budget: wpthlp dissipation term 1 [(m K)/s^2]","(m K)/s^2",zm)
        k = k + 1

      case ('wpthlp_mfl')
        iwpthlp_mfl = k
        call stat_assign(iwpthlp_mfl,"wpthlp_mfl", & 
             "wpthlp budget: wpthlp monotonic flux limiter [(m K)/s^2]","(m K)/s^2",zm)
        k = k + 1

      case ('wpthlp_cl')
        iwpthlp_cl = k
        call stat_assign(iwpthlp_cl,"wpthlp_cl", & 
             "wpthlp budget: wpthlp clipping term [(m K)/s^2]","(m K)/s^2",zm)
        k = k + 1

      case ('wpthlp_sicl')
        iwpthlp_sicl = k
        call stat_assign(iwpthlp_sicl,"wpthlp_sicl", & 
             "wpthlp budget: wpthlp semi-implicit clipping term [(m K)/s^2]","(m K)/s^2",zm)
        k = k + 1

      case ('wpthlp_forcing')
        iwpthlp_forcing = k

        call stat_assign( iwpthlp_forcing, "wpthlp_forcing", & 
             "wpthlp budget: wpthlp forcing (includes microphysics tendency) [(m K)/s^2]", &
             "(m K)/s^2", zm )
        k = k + 1

      case ('wpthlp_mc')
        iwpthlp_mc = k

        call stat_assign( iwpthlp_mc, "wpthlp_mc", & 
             "Microphysics tendency for wpthlp (not in budget) [(m K)/s^2]", &
             "(m K)/s^2", zm )
        k = k + 1

        ! Variance budgets
      case ('rtp2_bt')
        irtp2_bt = k
        call stat_assign(irtp2_bt,"rtp2_bt", & 
             "rtp2 budget: rtp2 time tendency [(kg^2)/(kg^2 s)]","(kg^2)/(kg^2 s)",zm)
        k = k + 1
      case ('rtp2_ma')
        irtp2_ma = k
        call stat_assign(irtp2_ma,"rtp2_ma", & 
             "rtp2 budget: rtp2 mean advection [(kg^2)/(kg^2 s)]","(kg^2)/(kg^2 s)",zm)
        k = k + 1
      case ('rtp2_ta')
        irtp2_ta = k
        call stat_assign(irtp2_ta,"rtp2_ta", & 
             "rtp2 budget: rtp2 turbulent advection [(kg^2)/(kg^2 s)]","(kg^2)/(kg^2 s)",zm)
        k = k + 1
      case ('rtp2_tp')
        irtp2_tp = k
        call stat_assign(irtp2_tp,"rtp2_tp", & 
             "rtp2 budget: rtp2 turbulent production [(kg^2)/(kg^2 s)]","(kg^2)/(kg^2 s)",zm)
        k = k + 1
      case ('rtp2_dp1')
        irtp2_dp1 = k
        call stat_assign(irtp2_dp1,"rtp2_dp1", & 
             "rtp2 budget: rtp2 dissipation term 1 [(kg^2)/(kg^2 s)]","(kg^2)/(kg^2 s)",zm)
        k = k + 1
      case ('rtp2_dp2')
        irtp2_dp2 = k
        call stat_assign(irtp2_dp2,"rtp2_dp2", & 
             "rtp2 budget: rtp2 dissipation term 2 [(kg^2)/(kg^2 s)]","(kg^2)/(kg^2 s)",zm)
        k = k + 1
      case ('rtp2_cl')
        irtp2_cl = k
        call stat_assign(irtp2_cl,"rtp2_cl", & 
             "rtp2 budget: rtp2 clipping term [(kg^2)/(kg^2 s)]","(kg^2)/(kg^2 s)",zm)
        k = k + 1

      case ('rtp2_pd')
        irtp2_pd = k
        call stat_assign( irtp2_pd, "rtp2_pd", & 
             "rtp2 budget: rtp2 positive definite adjustment [(kg^2)/(kg^2 s)]", &
             "(kg^2)/(kg^2 s)", zm )
        k = k + 1
        
      case ('rtp2_sf')
        irtp2_sf = k
        call stat_assign( irtp2_sf, "rtp2_sf", & 
             "rtp2 budget: rtp2 surface variance [(kg^2)/(kg^2 s)]", &
             "(kg^2)/(kg^2 s)", zm )
        k = k + 1

      case ('rtp2_forcing')
        irtp2_forcing = k

        call stat_assign( irtp2_forcing, "rtp2_forcing", & 
             "rtp2 budget: rtp2 forcing (includes microphysics tendency) [(kg/kg)^2/s]", &
             "(kg/kg)^2/s", zm )
        k = k + 1

      case ('rtp2_mc')
        irtp2_mc = k

        call stat_assign( irtp2_mc, "rtp2_mc", & 
             "Microphysics tendency for rtp2 (not in budget) [(kg/kg)^2/s]", &
             "(kg/kg)^2/s", zm )
        k = k + 1

      case ('thlp2_bt')
        ithlp2_bt = k
        call stat_assign(ithlp2_bt,"thlp2_bt", & 
             "thlp2 budget: thlp2 time tendency [(K^2)/s]","(K^2)/s",zm)
        k = k + 1
      case ('thlp2_ma')
        ithlp2_ma = k
        call stat_assign(ithlp2_ma,"thlp2_ma", & 
             "thlp2 budget: thlp2 mean advection [(K^2)/s]","(K^2)/s",zm)
        k = k + 1
      case ('thlp2_ta')
        ithlp2_ta = k
        call stat_assign(ithlp2_ta,"thlp2_ta", & 
             "thlp2 budget: thlp2 turbulent advection [(K^2)/s]","(K^2)/s",zm)
        k = k + 1
      case ('thlp2_tp')
        ithlp2_tp = k
        call stat_assign(ithlp2_tp,"thlp2_tp", & 
             "thlp2 budget: thlp2 turbulent production [(K^2)/s]","(K^2)/s",zm)
        k = k + 1
      case ('thlp2_dp1')
        ithlp2_dp1 = k
        call stat_assign(ithlp2_dp1,"thlp2_dp1", & 
             "thlp2 budget: thlp2 dissipation term 1 [(K^2)/s]","(K^2)/s",zm)
        k = k + 1
      case ('thlp2_dp2')
        ithlp2_dp2 = k
        call stat_assign(ithlp2_dp2,"thlp2_dp2", & 
             "thlp2 budget: thlp2 dissipation term 2 [(K^2)/s]","(K^2)/s",zm)
        k = k + 1
      case ('thlp2_cl')
        ithlp2_cl = k
        call stat_assign(ithlp2_cl,"thlp2_cl", & 
             "thlp2 budget: thlp2 clipping term [(K^2)/s]","(K^2)/s",zm)
        k = k + 1

      case ('thlp2_pd')
        ithlp2_pd = k
        call stat_assign( ithlp2_pd, "thlp2_pd", & 
             "thlp2 budget: thlp2 positive definite adjustment [(K^2)/s]", "K^2/s", zm )
        k = k + 1
        
      case ('thlp2_sf')
        ithlp2_sf = k
        call stat_assign( ithlp2_sf, "thlp2_sf", & 
             "thlp2 budget: thlp2 surface variance [(K^2)/s]", "K^2/s", zm )
        k = k + 1
      case ('thlp2_forcing')
        ithlp2_forcing = k
        call stat_assign( ithlp2_forcing, "thlp2_forcing", & 
             "thlp2 budget: thlp2 forcing (includes microphysics tendency) [K^2/s]", &
             "K^2/s", zm )
        k = k + 1
      case ('thlp2_mc')
        ithlp2_mc = k
        call stat_assign( ithlp2_mc, "thlp2_mc", & 
             "Microphysics tendency for thlp2 (not in budget) [K^2/s]", &
             "K^2/s", zm )
        k = k + 1

      case ('rtpthlp_bt')
        irtpthlp_bt = k
        call stat_assign(irtpthlp_bt,"rtpthlp_bt", & 
             "rtpthlp budget: rtpthlp time tendency [(kg K)/(kg s)]","(kg K)/(kg s)",zm)
        k = k + 1
      case ('rtpthlp_ma')
        irtpthlp_ma = k
        call stat_assign(irtpthlp_ma,"rtpthlp_ma", & 
             "rtpthlp budget: rtpthlp mean advection [(kg K)/(kg s)]","(kg K)/(kg s)",zm)
        k = k + 1
      case ('rtpthlp_ta')
        irtpthlp_ta = k
        call stat_assign(irtpthlp_ta,"rtpthlp_ta", & 
             "rtpthlp budget: rtpthlp turbulent advection [](kg K)/(kg s)","(kg K)/(kg s)",zm)
        k = k + 1
      case ('rtpthlp_tp1')
        irtpthlp_tp1 = k
        call stat_assign(irtpthlp_tp1,"rtpthlp_tp1", & 
             "rtpthlp budget: rtpthlp turbulent production 1 [(kg K)/(kg s)]","(kg K)/(kg s)",zm)
        k = k + 1
      case ('rtpthlp_tp2')
        irtpthlp_tp2 = k
        call stat_assign(irtpthlp_tp2,"rtpthlp_tp2", & 
             "rtpthlp budget: rtpthlp turbulent production 2 [(kg K)/(kg s)]","(kg K)/(kg s)",zm)
        k = k + 1
      case ('rtpthlp_dp1')
        irtpthlp_dp1 = k
        call stat_assign(irtpthlp_dp1,"rtpthlp_dp1", & 
             "rtpthlp budget: rtpthlp dissipation term 1 [(kg K)/(kg s)]","(kg K)/(kg s)",zm)
        k = k + 1
      case ('rtpthlp_dp2')
        irtpthlp_dp2 = k
        call stat_assign(irtpthlp_dp2,"rtpthlp_dp2", & 
             "rtpthlp budget: rtpthlp dissipation term 2 [(kg K)/(kg s)]","(kg K)/(kg s)",zm)
        k = k + 1
      case ('rtpthlp_cl')
        irtpthlp_cl = k
        call stat_assign(irtpthlp_cl,"rtpthlp_cl", & 
             "rtpthlp budget: rtpthlp clipping term [(kg K)/(kg s)]","(kg K)/(kg s)",zm)
        k = k + 1
      case ('rtpthlp_sf')
        irtpthlp_sf = k
        call stat_assign(irtpthlp_sf,"rtpthlp_sf", & 
             "rtpthlp budget: rtpthlp surface variance [(kg K)/(kg s)]","(kg K)/(kg s)",zm)
        k = k + 1
      case ('rtpthlp_forcing')
        irtpthlp_forcing = k
        call stat_assign( irtpthlp_forcing, "rtpthlp_forcing", & 
             "rtpthlp budget: rtpthlp forcing (includes microphysics tendency) [(K kg/kg)/s]", &
             "(K kg/kg)/s", zm )
        k = k + 1
      case ('rtpthlp_mc')
        irtpthlp_mc = k
        call stat_assign( irtpthlp_mc, "rtpthlp_mc", & 
             "Microphysics tendency for rtpthlp (not in budget) [(K kg/kg)/s]", &
             "(K kg/kg)/s", zm )
        k = k + 1

      case ('up2')
        iup2 = k
        call stat_assign(iup2,"up2", & 
             "u'^2 (momentum levels) [m^2/s^2]","m^2/s^2",zm)
        k = k + 1

      case ('vp2')
        ivp2 = k
        call stat_assign(ivp2,"vp2", & 
             "v'^2 (momentum levels) [m^2/s^2]","m^2/s^2",zm)
        k = k + 1

      case ('up2_bt')
        iup2_bt = k
        call stat_assign(iup2_bt,"up2_bt", & 
             "up2 budget: up2 time tendency [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('up2_ma')
        iup2_ma = k
        call stat_assign(iup2_ma,"up2_ma", & 
             "up2 budget: up2 mean advection [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('up2_ta')
        iup2_ta = k
        call stat_assign(iup2_ta,"up2_ta", & 
             "up2 budget: up2 turbulent advection [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('up2_tp')
        iup2_tp = k
        call stat_assign(iup2_tp,"up2_tp", & 
             "up2 budget: up2 turbulent production [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('up2_dp1')
        iup2_dp1 = k
        call stat_assign(iup2_dp1,"up2_dp1", & 
             "up2 budget: up2 dissipation term 1 [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('up2_dp2')
        iup2_dp2 = k
        call stat_assign(iup2_dp2,"up2_dp2", & 
             "up2 budget: up2 dissipation term 2 [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('up2_pr1')
        iup2_pr1 = k
        call stat_assign(iup2_pr1,"up2_pr1", & 
             "up2 budget: up2 pressure term 1 [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('up2_pr2')
        iup2_pr2 = k
        call stat_assign(iup2_pr2,"up2_pr2", & 
             "up2 budget: up2 pressure term 2 [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('up2_cl')
        iup2_cl = k
        call stat_assign(iup2_cl,"up2_cl", & 
             "up2 budget: up2 clipping [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('up2_pd')
        iup2_pd = k
        call stat_assign( iup2_pd, "up2_pd", & 
             "up2 budget: up2 positive definite adjustment [m^2/s^3]", "m^2/s^3", zm )
        k = k + 1
   
      case ('up2_sf')
        iup2_sf = k
        call stat_assign(iup2_sf,"up2_sf", & 
             "up2 budget: up2 surface variance [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('vp2_bt')
        ivp2_bt = k
        call stat_assign(ivp2_bt,"vp2_bt", & 
             "vp2 budget: vp2 time tendency [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('vp2_ma')
        ivp2_ma = k
        call stat_assign(ivp2_ma,"vp2_ma", & 
             "vp2 budget: vp2 mean advection [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('vp2_ta')
        ivp2_ta = k
        call stat_assign(ivp2_ta,"vp2_ta", & 
             "vp2 budget: vp2 turbulent advection [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('vp2_tp')
        ivp2_tp = k
        call stat_assign(ivp2_tp,"vp2_tp", & 
             "vp2 budget: vp2 turbulent production [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('vp2_dp1')
        ivp2_dp1 = k
        call stat_assign(ivp2_dp1,"vp2_dp1", & 
             "vp2 budget: vp2 dissipation term 1 [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('vp2_dp2')
        ivp2_dp2 = k
        call stat_assign(ivp2_dp2,"vp2_dp2", & 
             "vp2 budget: vp2 dissipation term 2 [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('vp2_pr1')
        ivp2_pr1 = k
        call stat_assign(ivp2_pr1,"vp2_pr1", & 
             "vp2 budget: vp2 pressure term 1 [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('vp2_pr2')
        ivp2_pr2 = k
        call stat_assign(ivp2_pr2,"vp2_pr2", & 
             "vp2 budget: vp2 pressure term 2 [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('vp2_cl')
        ivp2_cl = k
        call stat_assign(ivp2_cl,"vp2_cl", & 
             "vp2 budget: vp2 clipping [m^2/s^3]","m^2/s^3",zm)
        k = k + 1

      case ('vp2_pd')
        ivp2_pd = k
        call stat_assign( ivp2_pd, "vp2_pd", & 
             "vp2 budget: vp2 positive definite adjustment [m^2/s^3]", "m^2/s^3", zm )
        k = k + 1
        
      case ('vp2_sf')
        ivp2_sf = k
        call stat_assign( ivp2_sf, "vp2_sf", & 
             "vp2 budget: vp2 surface variance [m^2/s^3]", "m^2/s^3", zm )
        k = k + 1

      case ('wpthlp_entermfl')
        iwpthlp_entermfl = k
        call stat_assign( iwpthlp_entermfl, "wpthlp_entermfl", & 
             "Wpthlp entering flux limiter [(m K)/s]", "(m K)/s", zm )
        k = k + 1

      case ('wpthlp_exit_mfl')
        iwpthlp_exit_mfl = k
        call stat_assign( iwpthlp_exit_mfl, "wpthlp_exit_mfl", & 
             "Wpthlp exiting flux limiter [](m K)/s", "(m K)/s", zm )
        k = k + 1

      case ('wpthlp_mfl_min')
        iwpthlp_mfl_min = k
        call stat_assign( iwpthlp_mfl_min, "wpthlp_mfl_min", & 
             "Minimum allowable wpthlp [(m K)/s]", "(m K)/s", zm )
        k = k + 1

      case ('wpthlp_mfl_max')
        iwpthlp_mfl_max = k
        call stat_assign( iwpthlp_mfl_max, "wpthlp_mfl_max", & 
             "Maximum allowable wpthlp ((m K)/s) [(m K)/s]", "(m K)/s", zm )
        k = k + 1

      case ('wprtp_mfl_min')
        iwprtp_mfl_min = k
        call stat_assign( iwprtp_mfl_min, "wprtp_mfl_min", & 
             "Minimum allowable wprtp [(m kg)/(s kg)]", "(m kg)/(s kg)", zm )
        k = k + 1

      case ('wprtp_mfl_max')
        iwprtp_mfl_max = k
        call stat_assign( iwprtp_mfl_max, "wprtp_mfl_max", & 
             "Maximum allowable wprtp [(m kg)/(s kg)]", "(m kg)/(s kg)", zm )
        k = k + 1

      case ('wprtp_enter_mfl')
        iwprtp_enter_mfl = k
        call stat_assign( iwprtp_enter_mfl, "wprtp_enter_mfl", & 
             "Wprtp entering flux limiter [(m kg)/(s kg)]", "(m kg)/(s kg)", zm )
        k = k + 1

      case ('wprtp_exit_mfl')
        iwprtp_exit_mfl = k
        call stat_assign( iwprtp_exit_mfl, "wprtp_exit_mfl", & 
             "Wprtp exiting flux limiter [(m kg)/(s kg)]", "(m kg)/(s kg)", zm )
        k = k + 1        

      case ('wm_zm')
        iwm_zm = k
        call stat_assign( iwm_zm, "wm_zm", & 
             "Vertical (w) wind [m/s]", "m/s", zm )
        k = k + 1

      case ('cloud_frac_zm')
        icloud_frac_zm = k
        call stat_assign( icloud_frac_zm, "cloud_frac_zm", & 
                          "Cloud fraction", "count", zm )
        k = k + 1
      
      case ('ice_supersat_frac_zm')
        iice_supersat_frac_zm = k
        call stat_assign( iice_supersat_frac_zm, "ice_supersat_frac_zm", & 
                          "Ice cloud fraction", "count", zm )
        k = k + 1

      case ('rcm_zm')
        ircm_zm = k
        call stat_assign( ircm_zm, "rcm_zm", & 
             "Total water mixing ratio [kg/kg]", "kg/kg", zm )
        k = k + 1

      case ('rtm_zm')
        irtm_zm = k
        call stat_assign( irtm_zm, "rtm_zm", & 
             "Total water mixing ratio [kg/kg]", "kg/kg", zm )
        k = k + 1

      case ('thlm_zm')
        ithlm_zm = k
        call stat_assign( ithlm_zm, "thlm_zm", & 
             "Liquid potential temperature [K]", "K", zm )
        k = k + 1

      case ( 'Skw_velocity' )
        iSkw_velocity = k
        call stat_assign( iSkw_velocity, "Skw_velocity", & 
             "Skewness velocity [m/s]", "m/s", zm )
        k = k + 1

      case ( 'gamma_Skw_fnc' )
        igamma_Skw_fnc = k
        call stat_assign( igamma_Skw_fnc, "gamma_Skw_fnc", & 
             "Gamma as a function of skewness [-]", "count", zm )
        k = k + 1

      case ( 'C6rt_Skw_fnc' )
        iC6rt_Skw_fnc = k
        call stat_assign( iC6rt_Skw_fnc, "C6rt_Skw_fnc", & 
             "C_6rt parameter with Sk_w applied [-]", "count", zm )
        k = k + 1

      case ( 'C6thl_Skw_fnc' )
        iC6thl_Skw_fnc = k
        call stat_assign( iC6thl_Skw_fnc, "C6thl_Skw_fnc", & 
             "C_6thl parameter with Sk_w applied [-]", "count", zm )
        k = k + 1

      case ( 'C7_Skw_fnc' )
        iC7_Skw_fnc = k
        call stat_assign( iC7_Skw_fnc, "C7_Skw_fnc", & 
             "C_7 parameter with Sk_w applied [-]", "count", zm )
        k = k + 1

      case ( 'C1_Skw_fnc' )
        iC1_Skw_fnc = k
        call stat_assign( iC1_Skw_fnc, "C1_Skw_fnc", & 
             "C_1 parameter with Sk_w applied [-]", "count", zm )
        k = k + 1

      case ( 'a3_coef' )
        ia3_coef = k
        call stat_assign( ia3_coef, "a3_coef", & 
             "Quantity in formula 25 from Equations for CLUBB [-]", "count", zm )
        k = k + 1

      case ( 'wp3_on_wp2' )
        iwp3_on_wp2 = k
        call stat_assign( iwp3_on_wp2, "wp3_on_wp2", & 
             "Smoothed version of wp3 / wp2 [m/s]", "m/s", zm )
        k = k + 1

      case default
        l_found = .false.

        j = 1

        do while( j <= sclr_dim .and. .not. l_found )
          write( sclr_idx, * ) j
          sclr_idx = adjustl(sclr_idx)

          if( trim(vars_zm(i)) == 'sclr'//trim(sclr_idx)//'prtp'.and. .not. l_found ) then
            isclrprtp(j) = k

            call stat_assign(isclrprtp(j),"sclr"//trim(sclr_idx)//"prtp", & 
               "scalar("//trim(sclr_idx)//")'rt'","unknown",zm)
            k = k + 1
            l_found = .true.
          end if
          if( trim(vars_zm(i)) == 'sclr'//trim(sclr_idx)//'p2'.and. .not. l_found ) then
            isclrp2(j) = k
            call stat_assign(isclrp2(j) ,"sclr"//trim(sclr_idx)//"p2", & 
               "scalar("//trim(sclr_idx)//")'^2'","unknown",zm)
            k = k + 1
            l_found = .true.
          end if
          if( trim(vars_zm(i)) == 'sclr'//trim(sclr_idx)//'pthvp'.and. .not. l_found ) then
            isclrpthvp(j) = k
            call stat_assign(isclrpthvp(j),"sclr"//trim(sclr_idx)//"pthvp", & 
               "scalar("//trim(sclr_idx)//")'th_v'","unknown",zm)
            k = k + 1
            l_found = .true.
          end if
          if( trim(vars_zm(i)) == 'sclr'//trim(sclr_idx)//'pthlp'.and. .not. l_found ) then
            isclrpthlp(j) = k

            call stat_assign(isclrpthlp(j),"sclr"//trim(sclr_idx)//"pthlp", & 
               "scalar("//trim(sclr_idx)//")'th_l'","unknown",zm)
            k = k + 1
            l_found = .true.
          end if
          if( trim(vars_zm(i)) == 'sclr'//trim(sclr_idx)//'prcp'.and. .not. l_found ) then

            isclrprcp(j) = k

            call stat_assign(isclrprcp(j),"sclr"//trim(sclr_idx)//"prcp", & 
               "scalar("//trim(sclr_idx)//")'rc'","unknown",zm)
            k = k + 1
            l_found = .true.
          end if
          if( trim(vars_zm(i)) == 'wpsclr'//trim(sclr_idx)//'p'.and. .not. l_found ) then
            iwpsclrp(j) = k

            call stat_assign(iwpsclrp(j),"wpsclr"//trim(sclr_idx)//"p", & 
               "'w'scalar("//trim(sclr_idx)//")","unknown",zm)
            k = k + 1
            l_found = .true.
          end if
          if( trim(vars_zm(i)) == 'wpsclr'//trim(sclr_idx)//'p2'.and. .not. l_found ) then

            iwpsclrp2(j) = k

            call stat_assign(iwpsclrp2(j),"wpsclr"//trim(sclr_idx)//"p2", & 
               "'w'scalar("//trim(sclr_idx)//")'^2'","unknown",zm)
            k = k + 1
            l_found = .true.
          end if
          if( trim(vars_zm(i)) == 'wp2sclr'//trim(sclr_idx)//'p'.and. .not. l_found ) then

            iwp2sclrp(j) = k

            call stat_assign(iwp2sclrp(j) ,"wp2sclr"//trim(sclr_idx)//"p", & 
                 "'w'^2 scalar("//trim(sclr_idx)//")","unknown",zm)
            k = k + 1
            l_found = .true.
          end if
          if( trim(vars_zm(i)) == 'wpsclr'//trim(sclr_idx)//'prtp'.and. .not. l_found ) then
            iwpsclrprtp(j) = k

            call stat_assign( iwpsclrprtp(j),"wpsclr"//trim(sclr_idx)//"prtp", & 
               "'w' scalar("//trim(sclr_idx)//")'rt'","unknown",zm )
            k = k + 1
            l_found = .true.
          end if
          if( trim(vars_zm(i)) == 'wpsclr'//trim(sclr_idx)//'pthlp'.and. .not. l_found ) then
            iwpsclrpthlp(j) = k

            call stat_assign(iwpsclrpthlp(j),"wpsclr"//trim(sclr_idx)//"pthlp", & 
               "'w' scalar("//trim(sclr_idx)//")'th_l'","unknown",zm)
            k = k + 1
            l_found = .true.
          end if
          j = j + 1
        end do

        j = 1

        do while( j <= edsclr_dim .and. .not. l_found )

          write( sclr_idx, * ) j
          sclr_idx = adjustl(sclr_idx)

          if( trim(vars_zm(i)) == 'wpedsclr'//trim(sclr_idx)//'p'.and. .not. l_found ) then
            iwpedsclrp(j) = k

            call stat_assign(iwpedsclrp(j),"wpedsclr"//trim(sclr_idx)//"p", & 
               "eddy scalar("//trim(sclr_idx)//")'w'","unknown",zm)
            k = k + 1
            l_found = .true.
          end if

          j = j + 1

        end do

        if( .not. l_found ) then
          write(fstderr,*) 'Error:  unrecognized variable in vars_zm:  ',  trim(vars_zm(i))
          l_error = .true.  ! This will stop the run.
        end if
      end select

    end do

!   Non-interative diagnostics (zm)
!   iwp4, ircp2

!   if ( .not. clubb_at_least_debug_level( 1 ) ) then
!     if ( iwp4 + ircp2 + ishear > 0 ) then
!       write(fstderr,'(a)') &
!         "Warning: at debug level 0.  Non-interactive diagnostics will not be computed, "
!       write(fstderr,'(a)') "but some appear in the stats_zm namelist variable."
!     end if
!   end if

    return
  end subroutine stats_init_zm

end module stats_zm
