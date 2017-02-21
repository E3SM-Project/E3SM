!-----------------------------------------------------------------------
! $Id: stats_zm_module.F90 8116 2016-05-25 18:38:29Z raut@uwm.edu $
!===============================================================================
module stats_zm_module

  implicit none

  private ! Default Scope

  public :: stats_init_zm

  ! Constant parameters
  integer, parameter, public :: nvarmax_zm = 300  ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_zm( vars_zm, l_error )

! Description:
!   Initializes array indices for stats_zm

! Note:
!   All code that is within subroutine stats_init_zm, including variable
!   allocation code, is not called if l_stats is false.  This subroutine is
!   called only when l_stats is true.

!-----------------------------------------------------------------------

    use constants_clubb, only: &
        fstderr ! Constant(s)

    use stats_variables, only: & 
        stats_zm, & 
        iwp2, & 
        irtp2, & 
        ithlp2, & 
        irtpthlp, & 
        iwprtp, & 
        iwpthlp, & 
        iwp3_zm, &
        ithlp3_zm, &
        irtp3_zm, &
        iwp4, & 
        iwpthvp, & 
        irtpthvp, & 
        ithlpthvp, & 
        itau_zm, & 
        iKh_zm, & 
        iK_hm, & 
        iwprcp, & 
        irc_coef, &
        ithlprcp, & 
        irtprcp, & 
        ircp2,   &
        iSkw_zm, &
        iSkthl_zm, &
        iSkrt_zm

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
        iFcsed, &
        istability_correction

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
        iwpNcp

    use stats_variables, only: & 
        iVNr,  & 
        iVrr, &
        iVNc, & 
        iVrc, &
        iVNi, & 
        iVri, &
        iVNs, &
        iVrs, &
        iVrg, &
        iVrrprrp, &
        iVNrpNrp, &
        iVrrprrp_expcalc, &
        iVNrpNrp_expcalc

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
        iC1_Skw_fnc, &
        ibrunt_vaisala_freq_sqd, &
        iRichardson_num, &
        ishear_sqd, &
        ihydrometp2, &
        iwphydrometp, &
        irtphmp, &
        ithlphmp, &
        ihmxphmyp

    use stats_variables, only: &
        irtp2_from_chi

    use stats_variables, only: &
        ilh_rtp2_mc, &
        ilh_thlp2_mc, &
        ilh_wprtp_mc, &
        ilh_wpthlp_mc, &
        ilh_rtpthlp_mc

    use stats_type_utilities, only: & 
        stat_assign ! Procedure

    use parameters_model, only: &
        hydromet_dim, & ! Variable(s)
        sclr_dim,     &
        edsclr_dim

    use array_index, only: &
        hydromet_list, & ! Variable(s)
        l_mix_rat_hm

    implicit none

    ! External
    intrinsic :: trim

    ! Input Variable
    character(len= * ), dimension(nvarmax_zm), intent(in) :: vars_zm ! stats_zm variable names

    ! Input / Output Variable
    logical, intent(inout) :: l_error

    ! Local Varables
    integer :: tot_zm_loops

    integer :: hm_idx, hmx_idx, hmy_idx

    character(len=10) :: hm_type, hmx_type, hmy_type

    integer :: i, j, k

    character(len=50) :: sclr_idx

    ! The default initialization for array indices for stats_zm is zero (see module
    ! stats_variables)

    allocate( ihydrometp2(1:hydromet_dim) )
    allocate( iwphydrometp(1:hydromet_dim) )
    allocate( irtphmp(1:hydromet_dim) )
    allocate( ithlphmp(1:hydromet_dim) )
    allocate( ihmxphmyp(1:hydromet_dim,1:hydromet_dim) )
    allocate( iK_hm(1:hydromet_dim) )

    ihydrometp2(:) = 0
    iwphydrometp(:) = 0
    irtphmp(:) = 0
    ithlphmp(:) = 0
    ihmxphmyp(:,:) = 0
    iK_hm(:) = 0

    ! Allocate and then zero out passive scalar arrays on the stats_zm grid (fluxes,
    ! variances and other high-order moments)
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

    isclrprtp(:)    = 0
    isclrp2(:)      = 0
    isclrpthvp(:)   = 0
    isclrpthlp(:)   = 0
    isclrprcp(:)    = 0
    iwpsclrp(:)     = 0
    iwp2sclrp(:)    = 0
    iwpsclrp2(:)    = 0
    iwpsclrprtp(:)  = 0
    iwpsclrpthlp(:) = 0

    iwpedsclrp(:)   = 0

    ! Assign pointers for statistics variables stats_zm using stat_assign

    tot_zm_loops = stats_zm%num_output_fields

    if ( any( vars_zm == "hydrometp2" ) ) then
       ! Correct for number of variables found under "hydrometp2".
       ! Subtract 1 from the loop size for each hydrometeor.
       tot_zm_loops = tot_zm_loops - hydromet_dim
       ! Add 1 for "hydrometp2" to the loop size.
       tot_zm_loops = tot_zm_loops + 1
    endif

    if ( any( vars_zm == "wphydrometp" ) ) then
       ! Correct for number of variables found under "wphydrometp".
       ! Subtract 1 from the loop size for each hydrometeor.
       tot_zm_loops = tot_zm_loops - hydromet_dim
       ! Add 1 for "wphydrometp" to the loop size.
       tot_zm_loops = tot_zm_loops + 1
    endif

    if ( any( vars_zm == "rtphmp" ) ) then
       ! Correct for number of variables found under "rtphmp".
       ! Subtract 1 from the loop size for each hydrometeor.
       tot_zm_loops = tot_zm_loops - hydromet_dim
       ! Add 1 for "rtphmp" to the loop size.
       tot_zm_loops = tot_zm_loops + 1
    endif

    if ( any( vars_zm == "thlphmp" ) ) then
       ! Correct for number of variables found under "thlphmp".
       ! Subtract 1 from the loop size for each hydrometeor.
       tot_zm_loops = tot_zm_loops - hydromet_dim
       ! Add 1 for "thlphmp" to the loop size.
       tot_zm_loops = tot_zm_loops + 1
    endif

    if ( any( vars_zm == "hmxphmyp" ) ) then
       ! Correct for number of variables found under "hmxphmyp".
       ! Subtract the number of overall covariances of two hydrometeors, which
       ! is found by:  (1/2) * hydromet_dim * ( hydromet_dim - 1 );
       ! from the loop size.
       tot_zm_loops = tot_zm_loops - hydromet_dim * ( hydromet_dim - 1 ) / 2
       ! Add 1 for "hmxphmyp" to the loop size.
       tot_zm_loops = tot_zm_loops + 1
    endif

    if ( any( vars_zm == "K_hm" ) ) then
       ! Correct for number of variables found under "K_hm".
       ! Subtract 1 from the loop size for each hydrometeor.
       tot_zm_loops = tot_zm_loops - hydromet_dim
       ! Add 1 for "K_hm" to the loop size.
       tot_zm_loops = tot_zm_loops + 1
    endif

     if ( any( vars_zm == "sclrprtp" ) ) then
       ! Correct for number of variables found under "sclrprtp".
       ! Subtract 1 from the loop size for each scalar.
       tot_zm_loops = tot_zm_loops - sclr_dim
       ! Add 1 for "sclrprtp" to the loop size.
       tot_zm_loops = tot_zm_loops + 1
    endif   

    if ( any( vars_zm == "sclrp2" ) ) then
       ! Correct for number of variables found under "sclrp2".
       ! Subtract 1 from the loop size for each scalar.
       tot_zm_loops = tot_zm_loops - sclr_dim
       ! Add 1 for "sclrp2" to the loop size.
       tot_zm_loops = tot_zm_loops + 1
    endif

    if ( any( vars_zm == "sclrpthvp" ) ) then 
       ! Correct for number of variables found under "sclrpthvp".
       ! Subtract 1 from the loop size for each scalar.
       tot_zm_loops = tot_zm_loops - sclr_dim
       ! Add 1 for "sclrpthvp" to the loop size.
       tot_zm_loops = tot_zm_loops + 1
    endif

    if ( any( vars_zm == "sclrpthlp" ) ) then 
       ! Correct for number of variables found under "sclrpthlp".
       ! Subtract 1 from the loop size for each scalar.
       tot_zm_loops = tot_zm_loops - sclr_dim
       ! Add 1 for "sclrpthlp" to the loop size.
       tot_zm_loops = tot_zm_loops + 1
    endif

    if ( any( vars_zm == "sclrprcp" ) ) then 
       ! Correct for number of variables found under "sclrprcp".
       ! Subtract 1 from the loop size for each scalar.
       tot_zm_loops = tot_zm_loops - sclr_dim
       ! Add 1 for "sclrprcp" to the loop size.
       tot_zm_loops = tot_zm_loops + 1
    endif

    if ( any( vars_zm == "wpsclrp" ) ) then 
       ! Correct for number of variables found under "wpsclrp".
       ! Subtract 1 from the loop size for each scalar.
       tot_zm_loops = tot_zm_loops - sclr_dim
       ! Add 1 for "wpsclrp" to the loop size.
       tot_zm_loops = tot_zm_loops + 1
    endif

    if ( any( vars_zm == "wpsclrp2" ) ) then 
       ! Correct for number of variables found under "wpsclrp2".
       ! Subtract 1 from the loop size for each scalar.
       tot_zm_loops = tot_zm_loops - sclr_dim
       ! Add 1 for "wpsclrp2" to the loop size.
       tot_zm_loops = tot_zm_loops + 1
    endif

    if ( any( vars_zm == "wp2sclrp" ) ) then 
       ! Correct for number of variables found under "wp2sclrp".
       ! Subtract 1 from the loop size for each scalar.
       tot_zm_loops = tot_zm_loops - sclr_dim
       ! Add 1 for "wp2sclrp" to the loop size.
       tot_zm_loops = tot_zm_loops + 1
    endif

    if ( any( vars_zm == "wpsclrprtp" ) ) then 
       ! Correct for number of variables found under "wpsclrprtp".
       ! Subtract 1 from the loop size for each scalar.
       tot_zm_loops = tot_zm_loops - sclr_dim
       ! Add 1 for "wpsclrprtp" to the loop size.
       tot_zm_loops = tot_zm_loops + 1
    endif

    if ( any( vars_zm == "wpsclrpthlp" ) ) then 
       ! Correct for number of variables found under "wpsclrpthlp".
       ! Subtract 1 from the loop size for each scalar.
       tot_zm_loops = tot_zm_loops - sclr_dim
       ! Add 1 for "wpsclrpthlp" to the loop size.
       tot_zm_loops = tot_zm_loops + 1
    endif

    if ( any( vars_zm == "wpedsclrp" ) ) then 
       ! Correct for number of variables found under "wpedsclrp".
       ! Subtract 1 from the loop size for each scalar.
       tot_zm_loops = tot_zm_loops - edsclr_dim
       ! Add 1 for "wpedsclrp" to the loop size.
       tot_zm_loops = tot_zm_loops + 1
  endif



    k = 1

    do i = 1, tot_zm_loops

      select case ( trim( vars_zm(i) ) )

      case ('wp2')
        iwp2 = k
        call stat_assign( var_index=iwp2, var_name="wp2", &
             var_description="w'^2, Variance of vertical air velocity [m^2/s^2]", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('rtp2')
        irtp2 = k
        call stat_assign( var_index=irtp2, var_name="rtp2", &
             var_description="rt'^2, Variance of rt [(kg/kg)^2]", var_units="(kg/kg)^2", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('thlp2')
        ithlp2 = k
        call stat_assign( var_index=ithlp2, var_name="thlp2", &
             var_description="thl'^2, Variance of thl [K^2]", var_units="K^2", l_silhs=.false., &
             grid_kind=stats_zm )
        k = k + 1

      case ('rtpthlp')
        irtpthlp = k
        call stat_assign( var_index=irtpthlp, var_name="rtpthlp", &
             var_description="rt'thl', Covariance of rt and thl [(kg K)/kg]", &
             var_units="(kg K)/kg", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wprtp')
        iwprtp = k

        call stat_assign( var_index=iwprtp, var_name="wprtp", &
             var_description="w'rt', Vertical turbulent flux of rt [(kg/kg) m/s]", &
             var_units="(m kg)/(s kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wpthlp')
        iwpthlp = k

        call stat_assign( var_index=iwpthlp, var_name="wpthlp", &
             var_description="w'thl', Vertical turbulent flux of thl [K m/s]", &
             var_units="(m K)/s", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wp3_zm')
        iwp3_zm = k
        call stat_assign( var_index=iwp3_zm, var_name="wp3_zm", &
             var_description="w'^3 interpolated to moment. levels [m^3/s^3]", &
             var_units="(m^3)/(s^3)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('thlp3_zm')
        ithlp3_zm = k
        call stat_assign( var_index=ithlp3_zm, var_name="thlp3_zm", &
             var_description="thl'^3 interpolated to moment. levels [K^3]", &
             var_units="K^3", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('rtp3_zm')
        irtp3_zm = k
        call stat_assign( var_index=irtp3_zm, var_name="rtp3_zm", &
             var_description="rt'^3 interpolated to moment. levels [kg^3/kg^3]", &
             var_units="(kg^3)/(kg^3)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wp4')
        iwp4 = k
        call stat_assign( var_index=iwp4, var_name="wp4", var_description="w'^4 [m^4/s^4]", &
             var_units="(m^4)/(s^4)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wpthvp')
        iwpthvp = k
        call stat_assign( var_index=iwpthvp, var_name="wpthvp", &
             var_description="Buoyancy flux [K m/s]", var_units="K m/s", l_silhs=.false., &
             grid_kind=stats_zm )
        k = k + 1

      case ('rtpthvp')
        irtpthvp = k
        call stat_assign( var_index=irtpthvp, var_name="rtpthvp", &
             var_description="rt'thv' [(kg/kg) K]", var_units="(kg/kg) K", l_silhs=.false., &
             grid_kind=stats_zm )
        k = k + 1

      case ('thlpthvp')
        ithlpthvp = k
        call stat_assign( var_index=ithlpthvp, var_name="thlpthvp", &
          var_description="thl'thv' [K^2]", var_units="K^2", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('tau_zm')
        itau_zm = k

        call stat_assign( var_index=itau_zm, var_name="tau_zm", &
             var_description="Time-scale tau on momentum levels [s]", var_units="s", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('Kh_zm')
        iKh_zm = k

        call stat_assign( var_index=iKh_zm, var_name="Kh_zm", &
             var_description="Eddy diffusivity on momentum levels [m^2/s]", var_units="m^2/s", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('K_hm')

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            iK_hm(hm_idx) = k


            call stat_assign( var_index=iK_hm(hm_idx), &
                                 var_name="K_hm_"//trim( hm_type(1:2) ), &
                                 var_description="Eddy. diff. coef. of "  &
                                 // trim(hm_type(1:2)) &
                                 // " [m^2/s]", &
                                 var_units="[m^2/s]", &
                                 l_silhs=.false., grid_kind=stats_zm )

            k = k + 1

          end do


      case ('wprcp')
        iwprcp = k
        call stat_assign( var_index=iwprcp, var_name="wprcp", &
             var_description="w' rc' [(m/s) (kg/kg)]", var_units="(m/s) (kg/kg)", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      
      case ('rc_coef')
        irc_coef = k
        call stat_assign( var_index=irc_coef, var_name="rc_coef", &
             var_description="Coefficient of X' R_l' in Eq. (34)", var_units="[-]", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('thlprcp')
        ithlprcp = k
        call stat_assign( var_index=ithlprcp, var_name="thlprcp", &
             var_description="thl' rc' [K (kg/kg)]", var_units="K (kg/kg)", l_silhs=.false., &
             grid_kind=stats_zm )
        k = k + 1

      case ('rtprcp')
        irtprcp = k

        call stat_assign( var_index=irtprcp, var_name="rtprcp", &
             var_description="rt'rc' [(kg^2)/(kg^2)]", var_units="(kg^2)/(kg^2)", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('rcp2')
        ircp2 = k
        call stat_assign( var_index=ircp2, var_name="rcp2", &
             var_description="rc'^2 [(kg^2)/(kg^2)]", var_units="(kg^2)/(kg^2)", l_silhs=.false., &
             grid_kind=stats_zm )
        k = k + 1
      case ('upwp')
        iupwp = k
        call stat_assign( var_index=iupwp, var_name="upwp", &
             var_description="u'w', Vertical east-west momentum flux [m^2/s^2]", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('vpwp')
        ivpwp = k
        call stat_assign( var_index=ivpwp, var_name="vpwp", &
             var_description="v'w', Vertical north-south momentum flux [m^2/s^2]", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('rho_zm')
        irho_zm = k
        call stat_assign( var_index=irho_zm, var_name="rho_zm", &
             var_description="Density on momentum levels [kg/m^3]", var_units="kg m^{-3}", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('sigma_sqd_w')
        isigma_sqd_w = k
        call stat_assign( var_index=isigma_sqd_w, var_name="sigma_sqd_w", &
             var_description="Nondimensionalized w variance of Gaussian component [-]", &
             var_units="-", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('rho_ds_zm')
        irho_ds_zm = k
        call stat_assign( var_index=irho_ds_zm, var_name="rho_ds_zm", &
             var_description="Dry, static, base-state density [kg/m^3]", var_units="kg m^{-3}", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('thv_ds_zm')
        ithv_ds_zm = k
        call stat_assign( var_index=ithv_ds_zm, var_name="thv_ds_zm", &
             var_description="Dry, base-state theta_v [K]", var_units="K", l_silhs=.false., &
             grid_kind=stats_zm )
        k = k + 1
      case ('em')
        iem = k
        call stat_assign( var_index=iem, var_name="em", &
             var_description="Turbulent kinetic energy, usu. 0.5*(u'^2+v'^2+w'^2) [m^2/s^2]", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('shear')      ! Brian
        ishear = k
        call stat_assign( var_index=ishear, var_name="shear", &
             var_description="Wind shear production term [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('mean_w_up')
        imean_w_up = k
        call stat_assign( var_index=imean_w_up, var_name="mean_w_up", &
             var_description="Mean w >= w_ref [m/s]", var_units="m/s", l_silhs=.false., &
             grid_kind=stats_zm )
        k = k + 1
      case ('mean_w_down')
        imean_w_down = k
        call stat_assign( var_index=imean_w_down, var_name="mean_w_down", &
             var_description="Mean w <= w_ref [m/s]", var_units="m/s", l_silhs=.false., &
             grid_kind=stats_zm )
        k = k + 1
      case ('Frad')
        iFrad = k
        call stat_assign( var_index=iFrad, var_name="Frad", &
             var_description="Total (sw+lw) net (up+down) radiative flux [W/m^2]", &
             var_units="W/m^2", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('Frad_LW')    ! Brian
        iFrad_LW = k
        call stat_assign( var_index=iFrad_LW, var_name="Frad_LW", &
             var_description="Net long-wave radiative flux [W/m^2]", var_units="W/m^2", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('Frad_SW')    ! Brian
        iFrad_SW = k

        call stat_assign( var_index=iFrad_SW, var_name="Frad_SW", &
             var_description="Net short-wave radiative flux [W/m^2]", var_units="W/m^2", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('Frad_LW_up')    ! Brian
        iFrad_LW_up = k
        call stat_assign( var_index=iFrad_LW_up, var_name="Frad_LW_up", &
             var_description="Long-wave upwelling radiative flux [W/m^2]", var_units="W/m^2", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('Frad_SW_up')    ! Brian
        iFrad_SW_up = k

        call stat_assign( var_index=iFrad_SW_up, var_name="Frad_SW_up", &
             var_description="Short-wave upwelling radiative flux [W/m^2]", var_units="W/m^2", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('Frad_LW_down')    ! Brian
        iFrad_LW_down = k
        call stat_assign( var_index=iFrad_LW_down, var_name="Frad_LW_down", &
             var_description="Long-wave downwelling radiative flux [W/m^2]", var_units="W/m^2", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('Frad_SW_down')    ! Brian
        iFrad_SW_down = k

        call stat_assign( var_index=iFrad_SW_down, var_name="Frad_SW_down", &
             var_description="Short-wave downwelling radiative flux [W/m^2]", var_units="W/m^2", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1


      case ('Fprec')      ! Brian
        iFprec = k

        call stat_assign( var_index=iFprec, var_name="Fprec", &
             var_description="Rain flux [W/m^2]", var_units="W/m^2", l_silhs=.false., &
             grid_kind=stats_zm )
        k = k + 1

      case ('Fcsed')      ! Brian
        iFcsed = k

        call stat_assign( var_index=iFcsed, var_name="Fcsed", &
             var_description="cloud water sedimentation flux [kg/(s*m^2)]", &
             var_units="kg/(s*m^2)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case('hydrometp2')

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The overall variance of the hydrometeor.
            ihydrometp2(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=ihydrometp2(hm_idx), &
                                 var_name=trim( hm_type(1:2) )//"p2", &
                                 var_description="<" &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // "'^2> [(kg/kg)^2]", &
                                 var_units="(kg/kg)^2", &
                                 l_silhs=.false., grid_kind=stats_zm )

            else ! Concentration

               call stat_assign( var_index=ihydrometp2(hm_idx), &
                                 var_name=trim( hm_type(1:2) )//"p2", &
                                 var_description="<" &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // "'^2> [(num/kg)^2]", &
                                 var_units="(num/kg)^2", &
                                 l_silhs=.false., grid_kind=stats_zm )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ('wphydrometp')

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            iwphydrometp(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=iwphydrometp(hm_idx), &
                                 var_name="wp"//trim( hm_type(1:2) )//"p", &
                                 var_description="Covariance of w and " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " [(m/s) kg/kg]", &
                                 var_units="(m/s) kg/kg", &
                                 l_silhs=.false., grid_kind=stats_zm )

            else ! Concentration

               call stat_assign( var_index=iwphydrometp(hm_idx), &
                                 var_name="wp"//trim( hm_type(1:2) )//"p", &
                                 var_description="Covariance of w and " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " [(m/s) num/kg]", &
                                 var_units="(m/s) num/kg", &
                                 l_silhs=.false., grid_kind=stats_zm )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ('wpNcp')
        iwpNcp = k

        call stat_assign( var_index=iwpNcp, var_name="wpNcp", &
                          var_description="Covariance of w and " &
                                          // "N_c [(m/s) num/kg]", &
                          var_units="(m/s) num/kg", &
                          l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('rtphmp')

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            irtphmp(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=irtphmp(hm_idx), &
                                 var_name="rtp"//trim( hm_type(1:2) )//"p", &
                                 var_description="Covariance of r_t and " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " [kg^2/kg^2]", &
                                 var_units="kg^2/kg^2", &
                                 l_silhs=.false., grid_kind=stats_zm )

            else ! Concentration

               call stat_assign( var_index=irtphmp(hm_idx), &
                                 var_name="rtp"//trim( hm_type(1:2) )//"p", &
                                 var_description="Covariance of r_t and " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " [(kg/kg) num/kg]", &
                                 var_units="(kg/kg) num/kg", &
                                 l_silhs=.false., grid_kind=stats_zm )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ('thlphmp')

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ithlphmp(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=ithlphmp(hm_idx), &
                                 var_name="thlp"//trim( hm_type(1:2) )//"p", &
                                 var_description="Covariance of th_l and " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " [K kg/kg]", &
                                 var_units="K kg/kg", &
                                 l_silhs=.false., grid_kind=stats_zm )

            else ! Concentration

               call stat_assign( var_index=ithlphmp(hm_idx), &
                                 var_name="thlp"//trim( hm_type(1:2) )//"p", &
                                 var_description="Covariance of th_l and " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " [K num/kg]", &
                                 var_units="K num/kg", &
                                 l_silhs=.false., grid_kind=stats_zm )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ('hmxphmyp')

         do hmx_idx = 1, hydromet_dim, 1

            hmx_type = hydromet_list(hmx_idx)

            do hmy_idx = hmx_idx+1, hydromet_dim, 1

               hmy_type = hydromet_list(hmy_idx)

               ! The covariance (overall) of hmx and hmy.
               ihmxphmyp(hmy_idx,hmx_idx) = k

               if ( l_mix_rat_hm(hmx_idx) .and. l_mix_rat_hm(hmy_idx) ) then

                  ! Both hydrometeors are mixing ratios.
                  call stat_assign( var_index=ihmxphmyp(hmy_idx,hmx_idx), &
                                    var_name=trim( hmx_type(1:2) )//"p" &
                                    // trim( hmy_type(1:2) )//"p", &
                                    var_description="Covariance of " &
                                    // hmx_type(1:1)//"_"//trim(hmx_type(2:2)) &
                                    // " and " &
                                    // hmy_type(1:1)//"_"//trim(hmy_type(2:2)) &
                                    // " [(kg/kg)^2]", &
                                    var_units="(kg/kg)^2", l_silhs=.false., &
                                    grid_kind=stats_zm )

               elseif ( ( .not. l_mix_rat_hm(hmx_idx) ) &
                        .and. ( .not. l_mix_rat_hm(hmy_idx) ) ) then

                  ! Both hydrometeors are concentrations.
                  call stat_assign( var_index=ihmxphmyp(hmy_idx,hmx_idx), &
                                    var_name=trim( hmx_type(1:2) )//"p" &
                                    // trim( hmy_type(1:2) )//"p", &
                                    var_description="Covariance of " &
                                    // hmx_type(1:1)//"_"//trim(hmx_type(2:2)) &
                                    // " and " &
                                    // hmy_type(1:1)//"_"//trim(hmy_type(2:2)) &
                                    // " [(num/kg)^2]", &
                                    var_units="(num/kg)^2", l_silhs=.false., &
                                    grid_kind=stats_zm )

               else

                  ! One hydrometeor is a mixing ratio and the other hydrometeor
                  ! is a concentration.
                  call stat_assign( var_index=ihmxphmyp(hmy_idx,hmx_idx), &
                                    var_name=trim( hmx_type(1:2) )//"p" &
                                    // trim( hmy_type(1:2) )//"p", &
                                    var_description="Covariance of " &
                                    // hmx_type(1:1)//"_"//trim(hmx_type(2:2)) &
                                    // " and " &
                                    // hmy_type(1:1)//"_"//trim(hmy_type(2:2)) &
                                    // " [(kg/kg) num/kg]", &
                                    var_units="(kg/kg) num/kg", &
                                    l_silhs=.false., grid_kind=stats_zm )

               endif

               k = k + 1

            enddo ! hmy_idx = hmx_idx+1, hydromet_dim, 1

         enddo ! hmx_idx = 1, hydromet_dim, 1

      case ('VNr')
        iVNr = k

        call stat_assign( var_index=iVNr, var_name="VNr", &
             var_description="rrm concentration fallspeed [m/s]", var_units="m/s", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('Vrr')
        iVrr = k

        call stat_assign( var_index=iVrr, var_name="Vrr", &
             var_description="rrm mixing ratio fallspeed [m/s]", var_units="m/s", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('VNc')
        iVNc = k

        call stat_assign( var_index=iVNc, var_name="VNc", &
             var_description="Nrm concentration fallspeed [m/s]", var_units="m/s", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('Vrc')
        iVrc = k

        call stat_assign( var_index=iVrc, var_name="Vrc", &
             var_description="Nrm mixing ratio fallspeed [m/s]", var_units="m/s", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('VNs')
        iVNs = k

        call stat_assign( var_index=iVNs, var_name="VNs", &
             var_description="Snow concentration fallspeed [m/s]", var_units="m/s", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('Vrs')
        iVrs = k

        call stat_assign( var_index=iVrs, var_name="Vrs", &
             var_description="Snow mixing ratio fallspeed [m/s]", var_units="m/s", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('Vrg')
        iVrg = k

        call stat_assign( var_index=iVrg, var_name="Vrg", &
             var_description="Graupel sedimentation velocity [m/s]", var_units="m/s", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('VNi')
        iVNi = k

        call stat_assign( var_index=iVNi, var_name="VNi", &
             var_description="Cloud ice concentration fallspeed [m/s]", var_units="m/s", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('Vri')
        iVri = k

        call stat_assign( var_index=iVri, var_name="Vri", &
             var_description="Cloud ice mixing ratio fallspeed [m/s]", var_units="m/s", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('Vrrprrp')
        iVrrprrp = k

        call stat_assign( var_index=iVrrprrp, var_name="Vrrprrp", &
             var_description="Covariance of V_rr (r_r sed. vel.) and r_r [(m/s)(kg/kg)]", &
             var_units="(m/s)(kg/kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('VNrpNrp')
        iVNrpNrp = k

        call stat_assign( var_index=iVNrpNrp, var_name="VNrpNrp", &
             var_description="Covariance of V_Nr (N_r sed. vel.) and N_r [(m/s)(num/kg)]", &
             var_units="(m/s)(num/kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('Vrrprrp_expcalc')
        iVrrprrp_expcalc = k

        call stat_assign( var_index=iVrrprrp_expcalc, var_name="Vrrprrp_expcalc", &
             var_description="< V_rr'r_r' > (completely explicit calculation) [(m/s)(kg/kg)]", &
             var_units="(m/s)(kg/kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('VNrpNrp_expcalc')
        iVNrpNrp_expcalc = k

        call stat_assign( var_index=iVNrpNrp_expcalc, var_name="VNrpNrp_expcalc", &
             var_description="< V_Nr'N_r' > (completely explicit calculation) [(m/s)(num/kg)]", &
             var_units="(m/s)(num/kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wp2_bt')
        iwp2_bt = k

        call stat_assign( var_index=iwp2_bt, var_name="wp2_bt", &
             var_description="wp2 budget: wp2 time tendency [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wp2_ma')
        iwp2_ma = k

        call stat_assign( var_index=iwp2_ma, var_name="wp2_ma", &
             var_description="wp2 budget: wp2 vertical mean advection [m^2/s^3]", &
             var_units="m^2/s^3", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wp2_ta')
        iwp2_ta = k

        call stat_assign( var_index=iwp2_ta, var_name="wp2_ta", &
             var_description="wp2 budget: wp2 turbulent advection [m^2/s^3]", &
             var_units="m^2/s^3", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wp2_ac')
        iwp2_ac = k

        call stat_assign( var_index=iwp2_ac, var_name="wp2_ac", &
             var_description="wp2 budget: wp2 accumulation term [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wp2_bp')
        iwp2_bp = k

        call stat_assign( var_index=iwp2_bp, var_name="wp2_bp", &
             var_description="wp2 budget: wp2 buoyancy production [m^2/s^3]", &
             var_units="m^2/s^3", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wp2_pr1')
        iwp2_pr1 = k

        call stat_assign( var_index=iwp2_pr1, var_name="wp2_pr1", &
             var_description="wp2 budget: wp2 pressure term 1 [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wp2_pr2')
        iwp2_pr2 = k
        call stat_assign( var_index=iwp2_pr2, var_name="wp2_pr2", &
             var_description="wp2 budget: wp2 pressure term 2 [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wp2_pr3')
        iwp2_pr3 = k
        call stat_assign( var_index=iwp2_pr3, var_name="wp2_pr3", &
             var_description="wp2 budget: wp2 pressure term 3 [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )

        k = k + 1

      case ('wp2_dp1')
        iwp2_dp1 = k
        call stat_assign( var_index=iwp2_dp1, var_name="wp2_dp1", &
             var_description="wp2 budget: wp2 dissipation term 1 [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wp2_dp2')
        iwp2_dp2 = k
        call stat_assign( var_index=iwp2_dp2, var_name="wp2_dp2", &
             var_description="wp2 budget: wp2 dissipation term 2 [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )

        k = k + 1

      case ('wp2_cl')
        iwp2_cl = k

        call stat_assign( var_index=iwp2_cl, var_name="wp2_cl", &
             var_description="wp2 budget: wp2 clipping term [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )

        k = k + 1

      case ('wp2_pd')
        iwp2_pd = k

        call stat_assign( var_index=iwp2_pd, var_name="wp2_pd", &
             var_description="wp2 budget: wp2 positive definite adjustment [m^2/s^3]", &
             var_units="m2/s3", l_silhs=.false., grid_kind=stats_zm )

        k = k + 1
        
      case ('wp2_sf')
        iwp2_sf = k
        
        call stat_assign( var_index=iwp2_sf, var_name="wp2_sf", &
             var_description="wp2 budget: wp2 surface variance [m^2/s^3]", var_units="m2/s3", &
             l_silhs=.false., grid_kind=stats_zm )
             
        k = k + 1

      case ('wprtp_bt')
        iwprtp_bt = k
        call stat_assign( var_index=iwprtp_bt, var_name="wprtp_bt", &
             var_description="wprtp budget: wprtp time tendency [(m kg)/(s^2 kg)]", &
             var_units="(m kg)/(s^2 kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wprtp_ma')
        iwprtp_ma = k

        call stat_assign( var_index=iwprtp_ma, var_name="wprtp_ma", &
             var_description="wprtp budget: wprtp mean advection [(m kg)/(s^2 kg)]", &
             var_units="(m kg)/(s^2 kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wprtp_ta')
        iwprtp_ta = k

        call stat_assign( var_index=iwprtp_ta, var_name="wprtp_ta", &
             var_description="wprtp budget: wprtp turbulent advection [(m kg)/(s^2 kg)]", &
             var_units="(m kg)/(s^2 kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wprtp_tp')
        iwprtp_tp = k

        call stat_assign( var_index=iwprtp_tp, var_name="wprtp_tp", &
             var_description="wprtp budget: wprtp turbulent production [(m kg)/(s^2 kg)]", &
             var_units="(m kg)/(s^2 kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wprtp_ac')
        iwprtp_ac = k

        call stat_assign( var_index=iwprtp_ac, var_name="wprtp_ac", &
             var_description="wprtp budget: wprtp accumulation term [(m kg)/(s^2 kg)]", &
             var_units="(m kg)/(s^2 kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wprtp_bp')
        iwprtp_bp = k

        call stat_assign( var_index=iwprtp_bp, var_name="wprtp_bp", &
             var_description="wprtp budget: wprtp buoyancy production [(m kg)/(s^2 kg)]", &
             var_units="(m kg)/(s^2 kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wprtp_pr1')
        iwprtp_pr1 = k

        call stat_assign( var_index=iwprtp_pr1, var_name="wprtp_pr1", &
             var_description="wprtp budget: wprtp pressure term 1 [(m kg)/(s^2 kg)]", &
             var_units="(m kg)/(s^2 kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wprtp_pr2')
        iwprtp_pr2 = k

        call stat_assign( var_index=iwprtp_pr2, var_name="wprtp_pr2", &
             var_description="wprtp budget: wprtp pressure term 2 [(m kg)/(s^2 kg)]", &
             var_units="(m kg)/(s^2 kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wprtp_pr3')
        iwprtp_pr3 = k

        call stat_assign( var_index=iwprtp_pr3, var_name="wprtp_pr3", &
             var_description="wprtp budget: wprtp pressure term 3 [(m kg)/(s^2 kg)]", &
             var_units="(m kg)/(s^2 kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wprtp_dp1')
        iwprtp_dp1 = k

        call stat_assign( var_index=iwprtp_dp1, var_name="wprtp_dp1", &
             var_description="wprtp budget: wprtp dissipation term 1 [(m kg)/(s^2 kg)]", &
             var_units="(m kg)/(s^2 kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wprtp_mfl')
        iwprtp_mfl = k

        call stat_assign( var_index=iwprtp_mfl, var_name="wprtp_mfl", &
             var_description="wprtp budget: wprtp monotonic flux limiter [(m kg)/(s^2 kg)]", &
             var_units="(m kg)/(s^2 kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wprtp_cl')
        iwprtp_cl = k

        call stat_assign( var_index=iwprtp_cl, var_name="wprtp_cl", &
             var_description="wprtp budget: wprtp clipping term [(m kg)/(s^2 kg)]", &
             var_units="(m kg)/(s^2 kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wprtp_sicl')
        iwprtp_sicl = k

        call stat_assign( var_index=iwprtp_sicl, var_name="wprtp_sicl", &
             var_description="wprtp budget: wprtp semi-implicit clipping term [(m kg)/(s^2 kg)]", &
             var_units="(m kg)/(s^2 kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wprtp_pd')
        iwprtp_pd = k

        call stat_assign( var_index=iwprtp_pd, var_name="wprtp_pd", &
             var_description="wprtp budget: wprtp flux corrected trans. term [(m kg)/(s^2 kg)]", &
             var_units="(m kg)/(s^2 kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wprtp_forcing')
        iwprtp_forcing = k

        call stat_assign( var_index=iwprtp_forcing, var_name="wprtp_forcing", &
             var_description="wprtp budget: wprtp forcing (includes microphysics tendency) &
             &[(m kg/kg)/s^2]", &
             var_units="(m kg/kg)/s^2", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wprtp_mc')
        iwprtp_mc = k

        call stat_assign( var_index=iwprtp_mc, var_name="wprtp_mc", &
             var_description="Microphysics tendency for wprtp (not in budget) [(m kg/kg)/s^2]", &
             var_units="(m kg/kg)/s^2", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wpthlp_bt')
        iwpthlp_bt = k

        call stat_assign( var_index=iwpthlp_bt, var_name="wpthlp_bt", &
             var_description="wpthlp budget: [(m K)/s^2]", var_units="(m K)/s^2", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wpthlp_ma')
        iwpthlp_ma = k
        call stat_assign( var_index=iwpthlp_ma, var_name="wpthlp_ma", &
             var_description="wpthlp budget: wpthlp mean advection [(m K)/s^2]", &
             var_units="(m K)/s^2", l_silhs=.false., grid_kind=stats_zm )

        k = k + 1

      case ('wpthlp_ta')
        iwpthlp_ta = k
        call stat_assign( var_index=iwpthlp_ta, var_name="wpthlp_ta", &
             var_description="wpthlp budget: wpthlp turbulent advection [(m K)/s^2]", &
             var_units="(m K)/s^2", l_silhs=.false., grid_kind=stats_zm )

        k = k + 1

      case ('wpthlp_tp')
        iwpthlp_tp = k
        call stat_assign( var_index=iwpthlp_tp, var_name="wpthlp_tp", &
             var_description="wpthlp budget: wpthlp turbulent production [(m K)/s^2]", &
             var_units="(m K)/s^2", l_silhs=.false., grid_kind=stats_zm )

        k = k + 1

      case ('wpthlp_ac')
        iwpthlp_ac = k
        call stat_assign( var_index=iwpthlp_ac, var_name="wpthlp_ac", &
             var_description="wpthlp budget: wpthlp accumulation term [(m K)/s^2]", &
             var_units="(m K)/s^2", l_silhs=.false., grid_kind=stats_zm )

        k = k + 1

      case ('wpthlp_bp')
        iwpthlp_bp = k
        call stat_assign( var_index=iwpthlp_bp, var_name="wpthlp_bp", &
             var_description="wpthlp budget: wpthlp buoyancy production [(m K)/s^2]", &
             var_units="(m K)/s^2", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wpthlp_pr1')
        iwpthlp_pr1 = k

        call stat_assign( var_index=iwpthlp_pr1, var_name="wpthlp_pr1", &
             var_description="wpthlp budget: wpthlp pressure term 1 [(m K)/s^2]", &
             var_units="(m K)/s^2", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wpthlp_pr2')
        iwpthlp_pr2 = k

        call stat_assign( var_index=iwpthlp_pr2, var_name="wpthlp_pr2", &
             var_description="wpthlp budget: wpthlp pressure term 2 [(m K)/s^2]", &
             var_units="(m K)/s^2", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wpthlp_pr3')
        iwpthlp_pr3 = k
        call stat_assign( var_index=iwpthlp_pr3, var_name="wpthlp_pr3", &
             var_description="wpthlp budget: wpthlp pressure term 3 [(m K)/s^2]", &
             var_units="(m K)/s^2", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wpthlp_dp1')
        iwpthlp_dp1 = k
        call stat_assign( var_index=iwpthlp_dp1, var_name="wpthlp_dp1", &
             var_description="wpthlp budget: wpthlp dissipation term 1 [(m K)/s^2]", &
             var_units="(m K)/s^2", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wpthlp_mfl')
        iwpthlp_mfl = k
        call stat_assign( var_index=iwpthlp_mfl, var_name="wpthlp_mfl", &
             var_description="wpthlp budget: wpthlp monotonic flux limiter [(m K)/s^2]", &
             var_units="(m K)/s^2", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wpthlp_cl')
        iwpthlp_cl = k
        call stat_assign( var_index=iwpthlp_cl, var_name="wpthlp_cl", &
             var_description="wpthlp budget: wpthlp clipping term [(m K)/s^2]", &
             var_units="(m K)/s^2", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wpthlp_sicl')
        iwpthlp_sicl = k
        call stat_assign( var_index=iwpthlp_sicl, var_name="wpthlp_sicl", &
             var_description="wpthlp budget: wpthlp semi-implicit clipping term [(m K)/s^2]", &
             var_units="(m K)/s^2", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wpthlp_forcing')
        iwpthlp_forcing = k

        call stat_assign( var_index=iwpthlp_forcing, var_name="wpthlp_forcing", &
             var_description="wpthlp budget: wpthlp forcing (includes microphysics tendency) &
             &[(m K)/s^2]", &
             var_units="(m K)/s^2", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wpthlp_mc')
        iwpthlp_mc = k

        call stat_assign( var_index=iwpthlp_mc, var_name="wpthlp_mc", &
             var_description="Microphysics tendency for wpthlp (not in budget) [(m K)/s^2]", &
             var_units="(m K)/s^2", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

        ! Variance budgets
      case ('rtp2_bt')
        irtp2_bt = k
        call stat_assign( var_index=irtp2_bt, var_name="rtp2_bt", &
             var_description="rtp2 budget: rtp2 time tendency [(kg^2)/(kg^2 s)]", &
             var_units="(kg^2)/(kg^2 s)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('rtp2_ma')
        irtp2_ma = k
        call stat_assign( var_index=irtp2_ma, var_name="rtp2_ma", &
             var_description="rtp2 budget: rtp2 mean advection [(kg^2)/(kg^2 s)]", &
             var_units="(kg^2)/(kg^2 s)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('rtp2_ta')
        irtp2_ta = k
        call stat_assign( var_index=irtp2_ta, var_name="rtp2_ta", &
             var_description="rtp2 budget: rtp2 turbulent advection [(kg^2)/(kg^2 s)]", &
             var_units="(kg^2)/(kg^2 s)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('rtp2_tp')
        irtp2_tp = k
        call stat_assign( var_index=irtp2_tp, var_name="rtp2_tp", &
             var_description="rtp2 budget: rtp2 turbulent production [(kg^2)/(kg^2 s)]", &
             var_units="(kg^2)/(kg^2 s)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('rtp2_dp1')
        irtp2_dp1 = k
        call stat_assign( var_index=irtp2_dp1, var_name="rtp2_dp1", &
             var_description="rtp2 budget: rtp2 dissipation term 1 [(kg^2)/(kg^2 s)]", &
             var_units="(kg^2)/(kg^2 s)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('rtp2_dp2')
        irtp2_dp2 = k
        call stat_assign( var_index=irtp2_dp2, var_name="rtp2_dp2", &
             var_description="rtp2 budget: rtp2 dissipation term 2 [(kg^2)/(kg^2 s)]", &
             var_units="(kg^2)/(kg^2 s)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('rtp2_cl')
        irtp2_cl = k
        call stat_assign( var_index=irtp2_cl, var_name="rtp2_cl", &
             var_description="rtp2 budget: rtp2 clipping term [(kg^2)/(kg^2 s)]", &
             var_units="(kg^2)/(kg^2 s)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('rtp2_pd')
        irtp2_pd = k
        call stat_assign( var_index=irtp2_pd, var_name="rtp2_pd", &
             var_description="rtp2 budget: rtp2 positive definite adjustment [(kg^2)/(kg^2 s)]", &
             var_units="(kg^2)/(kg^2 s)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
        
      case ('rtp2_sf')
        irtp2_sf = k
        call stat_assign( var_index=irtp2_sf, var_name="rtp2_sf", &
             var_description="rtp2 budget: rtp2 surface variance [(kg^2)/(kg^2 s)]", &
             var_units="(kg^2)/(kg^2 s)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('rtp2_forcing')
        irtp2_forcing = k

        call stat_assign( var_index=irtp2_forcing, var_name="rtp2_forcing", &
             var_description="rtp2 budget: rtp2 forcing (includes microphysics tendency) &
             &[(kg/kg)^2/s]", &
             var_units="(kg/kg)^2/s", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('rtp2_mc')
        irtp2_mc = k

        call stat_assign( var_index=irtp2_mc, var_name="rtp2_mc", &
             var_description="Microphysics tendency for rtp2 (not in budget) [(kg/kg)^2/s]", &
             var_units="(kg/kg)^2/s", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('thlp2_bt')
        ithlp2_bt = k
        call stat_assign( var_index=ithlp2_bt, var_name="thlp2_bt", &
             var_description="thlp2 budget: thlp2 time tendency [(K^2)/s]", var_units="(K^2)/s", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('thlp2_ma')
        ithlp2_ma = k
        call stat_assign( var_index=ithlp2_ma, var_name="thlp2_ma", &
             var_description="thlp2 budget: thlp2 mean advection [(K^2)/s]", var_units="(K^2)/s", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('thlp2_ta')
        ithlp2_ta = k
        call stat_assign( var_index=ithlp2_ta, var_name="thlp2_ta", &
             var_description="thlp2 budget: thlp2 turbulent advection [(K^2)/s]", &
             var_units="(K^2)/s", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('thlp2_tp')
        ithlp2_tp = k
        call stat_assign( var_index=ithlp2_tp, var_name="thlp2_tp", &
             var_description="thlp2 budget: thlp2 turbulent production [(K^2)/s]", &
             var_units="(K^2)/s", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('thlp2_dp1')
        ithlp2_dp1 = k
        call stat_assign( var_index=ithlp2_dp1, var_name="thlp2_dp1", &
             var_description="thlp2 budget: thlp2 dissipation term 1 [(K^2)/s]", &
             var_units="(K^2)/s", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('thlp2_dp2')
        ithlp2_dp2 = k
        call stat_assign( var_index=ithlp2_dp2, var_name="thlp2_dp2", &
             var_description="thlp2 budget: thlp2 dissipation term 2 [(K^2)/s]", &
             var_units="(K^2)/s", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('thlp2_cl')
        ithlp2_cl = k
        call stat_assign( var_index=ithlp2_cl, var_name="thlp2_cl", &
             var_description="thlp2 budget: thlp2 clipping term [(K^2)/s]", var_units="(K^2)/s", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('thlp2_pd')
        ithlp2_pd = k
        call stat_assign( var_index=ithlp2_pd, var_name="thlp2_pd", &
             var_description="thlp2 budget: thlp2 positive definite adjustment [(K^2)/s]", &
             var_units="K^2/s", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
        
      case ('thlp2_sf')
        ithlp2_sf = k
        call stat_assign( var_index=ithlp2_sf, var_name="thlp2_sf", &
             var_description="thlp2 budget: thlp2 surface variance [(K^2)/s]", var_units="K^2/s", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('thlp2_forcing')
        ithlp2_forcing = k
        call stat_assign( var_index=ithlp2_forcing, var_name="thlp2_forcing", &
             var_description="thlp2 budget: thlp2 forcing (includes microphysics tendency) &
             &[K^2/s]", &
             var_units="K^2/s", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('thlp2_mc')
        ithlp2_mc = k
        call stat_assign( var_index=ithlp2_mc, var_name="thlp2_mc", &
             var_description="Microphysics tendency for thlp2 (not in budget) [K^2/s]", &
             var_units="K^2/s", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('rtpthlp_bt')
        irtpthlp_bt = k
        call stat_assign( var_index=irtpthlp_bt, var_name="rtpthlp_bt", &
             var_description="rtpthlp budget: rtpthlp time tendency [(kg K)/(kg s)]", &
             var_units="(kg K)/(kg s)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('rtpthlp_ma')
        irtpthlp_ma = k
        call stat_assign( var_index=irtpthlp_ma, var_name="rtpthlp_ma", &
             var_description="rtpthlp budget: rtpthlp mean advection [(kg K)/(kg s)]", &
             var_units="(kg K)/(kg s)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('rtpthlp_ta')
        irtpthlp_ta = k
        call stat_assign( var_index=irtpthlp_ta, var_name="rtpthlp_ta", &
             var_description="rtpthlp budget: rtpthlp turbulent advection [](kg K)/(kg s)", &
             var_units="(kg K)/(kg s)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('rtpthlp_tp1')
        irtpthlp_tp1 = k
        call stat_assign( var_index=irtpthlp_tp1, var_name="rtpthlp_tp1", &
             var_description="rtpthlp budget: rtpthlp turbulent production 1 [(kg K)/(kg s)]", &
             var_units="(kg K)/(kg s)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('rtpthlp_tp2')
        irtpthlp_tp2 = k
        call stat_assign( var_index=irtpthlp_tp2, var_name="rtpthlp_tp2", &
             var_description="rtpthlp budget: rtpthlp turbulent production 2 [(kg K)/(kg s)]", &
             var_units="(kg K)/(kg s)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('rtpthlp_dp1')
        irtpthlp_dp1 = k
        call stat_assign( var_index=irtpthlp_dp1, var_name="rtpthlp_dp1", &
             var_description="rtpthlp budget: rtpthlp dissipation term 1 [(kg K)/(kg s)]", &
             var_units="(kg K)/(kg s)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('rtpthlp_dp2')
        irtpthlp_dp2 = k
        call stat_assign( var_index=irtpthlp_dp2, var_name="rtpthlp_dp2", &
             var_description="rtpthlp budget: rtpthlp dissipation term 2 [(kg K)/(kg s)]", &
             var_units="(kg K)/(kg s)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('rtpthlp_cl')
        irtpthlp_cl = k
        call stat_assign( var_index=irtpthlp_cl, var_name="rtpthlp_cl", &
             var_description="rtpthlp budget: rtpthlp clipping term [(kg K)/(kg s)]", &
             var_units="(kg K)/(kg s)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('rtpthlp_sf')
        irtpthlp_sf = k
        call stat_assign( var_index=irtpthlp_sf, var_name="rtpthlp_sf", &
             var_description="rtpthlp budget: rtpthlp surface variance [(kg K)/(kg s)]", &
             var_units="(kg K)/(kg s)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('rtpthlp_forcing')
        irtpthlp_forcing = k
        call stat_assign( var_index=irtpthlp_forcing, var_name="rtpthlp_forcing", &
             var_description="rtpthlp budget: rtpthlp forcing (includes microphysics tendency) &
             &[(K kg/kg)/s]", &
             var_units="(K kg/kg)/s", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
      case ('rtpthlp_mc')
        irtpthlp_mc = k
        call stat_assign( var_index=irtpthlp_mc, var_name="rtpthlp_mc", &
             var_description="Microphysics tendency for rtpthlp (not in budget) [(K kg/kg)/s]", &
             var_units="(K kg/kg)/s", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('up2')
        iup2 = k
        call stat_assign( var_index=iup2, var_name="up2", &
             var_description="u'^2 (momentum levels) [m^2/s^2]", var_units="m^2/s^2", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('vp2')
        ivp2 = k
        call stat_assign( var_index=ivp2, var_name="vp2", &
             var_description="v'^2 (momentum levels) [m^2/s^2]", var_units="m^2/s^2", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('up2_bt')
        iup2_bt = k
        call stat_assign( var_index=iup2_bt, var_name="up2_bt", &
             var_description="up2 budget: up2 time tendency [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('up2_ma')
        iup2_ma = k
        call stat_assign( var_index=iup2_ma, var_name="up2_ma", &
             var_description="up2 budget: up2 mean advection [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('up2_ta')
        iup2_ta = k
        call stat_assign( var_index=iup2_ta, var_name="up2_ta", &
             var_description="up2 budget: up2 turbulent advection [m^2/s^3]", &
             var_units="m^2/s^3", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('up2_tp')
        iup2_tp = k
        call stat_assign( var_index=iup2_tp, var_name="up2_tp", &
             var_description="up2 budget: up2 turbulent production [m^2/s^3]", &
             var_units="m^2/s^3", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('up2_dp1')
        iup2_dp1 = k
        call stat_assign( var_index=iup2_dp1, var_name="up2_dp1", &
             var_description="up2 budget: up2 dissipation term 1 [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('up2_dp2')
        iup2_dp2 = k
        call stat_assign( var_index=iup2_dp2, var_name="up2_dp2", &
             var_description="up2 budget: up2 dissipation term 2 [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('up2_pr1')
        iup2_pr1 = k
        call stat_assign( var_index=iup2_pr1, var_name="up2_pr1", &
             var_description="up2 budget: up2 pressure term 1 [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('up2_pr2')
        iup2_pr2 = k
        call stat_assign( var_index=iup2_pr2, var_name="up2_pr2", &
             var_description="up2 budget: up2 pressure term 2 [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('up2_cl')
        iup2_cl = k
        call stat_assign( var_index=iup2_cl, var_name="up2_cl", &
             var_description="up2 budget: up2 clipping [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('up2_pd')
        iup2_pd = k
        call stat_assign( var_index=iup2_pd, var_name="up2_pd", &
             var_description="up2 budget: up2 positive definite adjustment [m^2/s^3]", &
             var_units="m^2/s^3", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
   
      case ('up2_sf')
        iup2_sf = k
        call stat_assign( var_index=iup2_sf, var_name="up2_sf", &
             var_description="up2 budget: up2 surface variance [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('vp2_bt')
        ivp2_bt = k
        call stat_assign( var_index=ivp2_bt, var_name="vp2_bt", &
             var_description="vp2 budget: vp2 time tendency [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('vp2_ma')
        ivp2_ma = k
        call stat_assign( var_index=ivp2_ma, var_name="vp2_ma", &
             var_description="vp2 budget: vp2 mean advection [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('vp2_ta')
        ivp2_ta = k
        call stat_assign( var_index=ivp2_ta, var_name="vp2_ta", &
             var_description="vp2 budget: vp2 turbulent advection [m^2/s^3]", &
             var_units="m^2/s^3", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('vp2_tp')
        ivp2_tp = k
        call stat_assign( var_index=ivp2_tp, var_name="vp2_tp", &
             var_description="vp2 budget: vp2 turbulent production [m^2/s^3]", &
             var_units="m^2/s^3", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('vp2_dp1')
        ivp2_dp1 = k
        call stat_assign( var_index=ivp2_dp1, var_name="vp2_dp1", &
             var_description="vp2 budget: vp2 dissipation term 1 [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('vp2_dp2')
        ivp2_dp2 = k
        call stat_assign( var_index=ivp2_dp2, var_name="vp2_dp2", &
             var_description="vp2 budget: vp2 dissipation term 2 [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('vp2_pr1')
        ivp2_pr1 = k
        call stat_assign( var_index=ivp2_pr1, var_name="vp2_pr1", &
             var_description="vp2 budget: vp2 pressure term 1 [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('vp2_pr2')
        ivp2_pr2 = k
        call stat_assign( var_index=ivp2_pr2, var_name="vp2_pr2", &
             var_description="vp2 budget: vp2 pressure term 2 [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('vp2_cl')
        ivp2_cl = k
        call stat_assign( var_index=ivp2_cl, var_name="vp2_cl", &
             var_description="vp2 budget: vp2 clipping [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('vp2_pd')
        ivp2_pd = k
        call stat_assign( var_index=ivp2_pd, var_name="vp2_pd", &
             var_description="vp2 budget: vp2 positive definite adjustment [m^2/s^3]", &
             var_units="m^2/s^3", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1
        
      case ('vp2_sf')
        ivp2_sf = k
        call stat_assign( var_index=ivp2_sf, var_name="vp2_sf", &
             var_description="vp2 budget: vp2 surface variance [m^2/s^3]", var_units="m^2/s^3", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wpthlp_entermfl')
        iwpthlp_entermfl = k
        call stat_assign( var_index=iwpthlp_entermfl, var_name="wpthlp_entermfl", &
             var_description="Wpthlp entering flux limiter [(m K)/s]", var_units="(m K)/s", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wpthlp_exit_mfl')
        iwpthlp_exit_mfl = k
        call stat_assign( var_index=iwpthlp_exit_mfl, var_name="wpthlp_exit_mfl", &
             var_description="Wpthlp exiting flux limiter [](m K)/s", var_units="(m K)/s", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wpthlp_mfl_min')
        iwpthlp_mfl_min = k
        call stat_assign( var_index=iwpthlp_mfl_min, var_name="wpthlp_mfl_min", &
             var_description="Minimum allowable wpthlp [(m K)/s]", var_units="(m K)/s", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wpthlp_mfl_max')
        iwpthlp_mfl_max = k
        call stat_assign( var_index=iwpthlp_mfl_max, var_name="wpthlp_mfl_max", &
             var_description="Maximum allowable wpthlp ((m K)/s) [(m K)/s]", var_units="(m K)/s", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wprtp_mfl_min')
        iwprtp_mfl_min = k
        call stat_assign( var_index=iwprtp_mfl_min, var_name="wprtp_mfl_min", &
             var_description="Minimum allowable wprtp [(m kg)/(s kg)]", &
             var_units="(m kg)/(s kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wprtp_mfl_max')
        iwprtp_mfl_max = k
        call stat_assign( var_index=iwprtp_mfl_max, var_name="wprtp_mfl_max", &
             var_description="Maximum allowable wprtp [(m kg)/(s kg)]", &
             var_units="(m kg)/(s kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wprtp_enter_mfl')
        iwprtp_enter_mfl = k
        call stat_assign( var_index=iwprtp_enter_mfl, var_name="wprtp_enter_mfl", &
             var_description="Wprtp entering flux limiter [(m kg)/(s kg)]", &
             var_units="(m kg)/(s kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wprtp_exit_mfl')
        iwprtp_exit_mfl = k
        call stat_assign( var_index=iwprtp_exit_mfl, var_name="wprtp_exit_mfl", &
             var_description="Wprtp exiting flux limiter [(m kg)/(s kg)]", &
             var_units="(m kg)/(s kg)", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('wm_zm')
        iwm_zm = k
        call stat_assign( var_index=iwm_zm, var_name="wm_zm", &
             var_description="Vertical (w) wind [m/s]", var_units="m/s", l_silhs=.false., &
             grid_kind=stats_zm )
        k = k + 1

      case ('cloud_frac_zm')
        icloud_frac_zm = k
        call stat_assign( var_index=icloud_frac_zm, var_name="cloud_frac_zm", &
          var_description="Cloud fraction [-]", var_units="-", l_silhs=.false., grid_kind=stats_zm)
        k = k + 1
      
      case ('ice_supersat_frac_zm')
        iice_supersat_frac_zm = k
        call stat_assign( var_index=iice_supersat_frac_zm, var_name="ice_supersat_frac_zm", &
             var_description="Ice cloud fraction", var_units="count", l_silhs=.false., &
             grid_kind=stats_zm )
        k = k + 1

      case ('rcm_zm')
        ircm_zm = k
        call stat_assign( var_index=ircm_zm, var_name="rcm_zm", &
             var_description="Total water mixing ratio [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('rtm_zm')
        irtm_zm = k
        call stat_assign( var_index=irtm_zm, var_name="rtm_zm", &
             var_description="Total water mixing ratio [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ('thlm_zm')
        ithlm_zm = k
        call stat_assign( var_index=ithlm_zm, var_name="thlm_zm", &
             var_description="Liquid potential temperature [K]", var_units="K", l_silhs=.false., &
             grid_kind=stats_zm )
        k = k + 1

      case ( 'Skw_velocity' )
        iSkw_velocity = k
        call stat_assign( var_index=iSkw_velocity, var_name="Skw_velocity", &
             var_description="Skewness velocity [m/s]", var_units="m/s", l_silhs=.false., &
             grid_kind=stats_zm )
        k = k + 1

      case ( 'gamma_Skw_fnc' )
        igamma_Skw_fnc = k
        call stat_assign( var_index=igamma_Skw_fnc, var_name="gamma_Skw_fnc", &
             var_description="Gamma as a function of skewness [-]", var_units="-", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ( 'C6rt_Skw_fnc' )
        iC6rt_Skw_fnc = k
        call stat_assign( var_index=iC6rt_Skw_fnc, var_name="C6rt_Skw_fnc", &
             var_description="C_6rt parameter with Sk_w applied [-]", var_units="-", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ( 'C6thl_Skw_fnc' )
        iC6thl_Skw_fnc = k
        call stat_assign( var_index=iC6thl_Skw_fnc, var_name="C6thl_Skw_fnc", &
             var_description="C_6thl parameter with Sk_w applied [-]", var_units="-", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ( 'C7_Skw_fnc' )
        iC7_Skw_fnc = k
        call stat_assign( var_index=iC7_Skw_fnc, var_name="C7_Skw_fnc", &
             var_description="C_7 parameter with Sk_w applied [-]", var_units="-", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ( 'C1_Skw_fnc' )
        iC1_Skw_fnc = k
        call stat_assign( var_index=iC1_Skw_fnc, var_name="C1_Skw_fnc", &
             var_description="C_1 parameter with Sk_w applied [-]", var_units="-", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ( 'brunt_vaisala_freq_sqd' )
        ibrunt_vaisala_freq_sqd = k
        call stat_assign( var_index=ibrunt_vaisala_freq_sqd, var_name="brunt_vaisala_freq_sqd", &
             var_description="Brunt-Vaisala freqency squared, N^2 [1/s^2]", var_units="1/s^2", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ( 'Richardson_num' )
        iRichardson_num = k
        call stat_assign( var_index=iRichardson_num, var_name="Richardson_num", &
             var_description="Richardson number [-]", var_units="-", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ( 'shear_sqd' )
        ishear_sqd = k
        call stat_assign( var_index=ishear_sqd, var_name="shear_sqd", &
             var_description="shear_sqd [-]", var_units="-", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ( 'a3_coef' )
        ia3_coef = k
        call stat_assign( var_index=ia3_coef, var_name="a3_coef", &
             var_description="Quantity in formula 25 from Equations for CLUBB [-]", &
             var_units="count", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ( 'wp3_on_wp2' )
        iwp3_on_wp2 = k
        call stat_assign( var_index=iwp3_on_wp2, var_name="wp3_on_wp2", &
             var_description="Smoothed version of wp3 / wp2 [m/s]", var_units="m/s", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ( 'Skw_zm' )
        iSkw_zm = k
        call stat_assign( var_index=iSkw_zm, var_name="Skw_zm", &
             var_description="Skewness of w on momentum levels [-]", var_units="-", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ( 'Skthl_zm' )
        iSkthl_zm = k
        call stat_assign( var_index=iSkthl_zm, var_name="Skthl_zm", &
             var_description="Skewness of thl on momentum levels [-]", var_units="-", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ( 'Skrt_zm' )
        iSkrt_zm = k
        call stat_assign( var_index=iSkrt_zm, var_name="Skrt_zm", &
             var_description="Skewness of rt on momentum levels [-]", var_units="-", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ( 'stability_correction' )
        istability_correction = k
        call stat_assign( var_index=istability_correction, var_name="stability_correction", &
             var_description="Stability applied to diffusion of rtm and thlm [-]", var_units="-", &
             l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ( 'rtp2_from_chi' )
        irtp2_from_chi = k
        call stat_assign( var_index=irtp2_from_chi, var_name="rtp2_from_chi", &
             var_description="Variance of rt, computed from the chi/eta distribution [(kg/kg)^2]", &
             var_units="(kg/kg)^2", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ( 'lh_rtp2_mc' )
        ilh_rtp2_mc = k
        call stat_assign( var_index=ilh_rtp2_mc, var_name="lh_rtp2_mc", &
             var_description="LH est. of rtp2_mc [(kg/kg)^2/s]", &
             var_units="(kg/kg)^2/s", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ( 'lh_thlp2_mc' )
        ilh_thlp2_mc = k
        call stat_assign( var_index=ilh_thlp2_mc, var_name="lh_thlp2_mc", &
             var_description="LH est. of thlp2_mc [K^2/s]", &
             var_units="K^2/s", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ( 'lh_wprtp_mc' )
        ilh_wprtp_mc = k
        call stat_assign( var_index=ilh_wprtp_mc, var_name="lh_wprtp_mc", &
             var_description="LH est. of wprtp_mc [(m kg/kg)/s^2]", &
             var_units="(m kg/kg)/s^2", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ( 'lh_wpthlp_mc' )
        ilh_wpthlp_mc = k
        call stat_assign( var_index=ilh_wpthlp_mc, var_name="lh_wpthlp_mc", &
             var_description="LH est. of wpthlp_mc [(m K)/s^2]", &
             var_units="(m K)/s^2", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ( 'lh_rtpthlp_mc' )
        ilh_rtpthlp_mc = k
        call stat_assign( var_index=ilh_rtpthlp_mc, var_name="lh_rtpthlp_mc", &
             var_description="LH est. of rtpthlp_mc [(K kg/kg)/s]", &
             var_units="(K kg/kg)/s", l_silhs=.false., grid_kind=stats_zm )
        k = k + 1

      case ( 'sclrprtp' )
        do j = 1, sclr_dim, 1
          write( sclr_idx, * ) j
          sclr_idx = adjustl(sclr_idx)
          isclrprtp(j) = k
          call stat_assign( var_index=isclrprtp(j), var_name="sclr"//trim(sclr_idx)//"prtp", &
            var_description="scalar("//trim(sclr_idx)//")'rt'", var_units="unknown", &
            l_silhs=.false., grid_kind=stats_zm )
          k = k + 1
        end do

      case ( 'sclrp2' )
        do j = 1, sclr_dim, 1
          write( sclr_idx, * ) j
          sclr_idx = adjustl(sclr_idx)
          isclrp2(j) = k
          call stat_assign( var_index=isclrp2(j), var_name="sclr"//trim(sclr_idx)//"p2", &
            var_description="scalar("//trim(sclr_idx)//")'^2'", var_units="unknown", &
            l_silhs=.false., grid_kind=stats_zm )
          k = k + 1
        end do

      case ( 'sclrpthvp' )
        do j = 1, sclr_dim, 1
          write( sclr_idx, * ) j
          sclr_idx = adjustl(sclr_idx)
          isclrpthvp(j) = k
          call stat_assign( var_index=isclrpthvp(j), var_name="sclr"//trim(sclr_idx)//"pthvp", &
            var_description="scalar("//trim(sclr_idx)//")'th_v'", var_units="unknown", &
            l_silhs=.false., grid_kind=stats_zm )
          k = k + 1
        end do

      case ( 'sclrpthlp' )
        do j = 1, sclr_dim, 1
          write( sclr_idx, * ) j
          sclr_idx = adjustl(sclr_idx)
          isclrpthlp(j) = k
          call stat_assign( var_index=isclrpthlp(j), var_name="sclr"//trim(sclr_idx)//"pthlp", &
            var_description="scalar("//trim(sclr_idx)//")'th_l'", var_units="unknown", &
            l_silhs=.false., grid_kind=stats_zm )
          k = k + 1
        end do

      case ( 'sclrprcp' )
        do j = 1, sclr_dim, 1
          write( sclr_idx, * ) j
          sclr_idx = adjustl(sclr_idx)
          isclrprcp(j) = k
          call stat_assign( var_index=isclrprcp(j), var_name="sclr"//trim(sclr_idx)//"prcp", &
            var_description="scalar("//trim(sclr_idx)//")'rc'", var_units="unknown", &
            l_silhs=.false., grid_kind=stats_zm )
          k = k + 1
        end do

      case ( 'wpsclrp' )
        do j = 1, sclr_dim, 1
          write( sclr_idx, * ) j
          sclr_idx = adjustl(sclr_idx)
          iwpsclrp(j) = k
          call stat_assign( var_index=iwpsclrp(j), var_name="wpsclr"//trim(sclr_idx)//"p", &
            var_description="'w'scalar("//trim(sclr_idx)//")", var_units="unknown", &
            l_silhs=.false., grid_kind=stats_zm )
          k = k + 1
        end do

      case ( 'wpsclrp2' )
        do j = 1, sclr_dim, 1
          write( sclr_idx, * ) j
          sclr_idx = adjustl(sclr_idx)
          iwpsclrp2(j) = k
          call stat_assign( var_index=iwpsclrp2(j), var_name="wpsclr"//trim(sclr_idx)//"p2", &
            var_description="'w'scalar("//trim(sclr_idx)//")'^2'", var_units="unknown", &
            l_silhs=.false., grid_kind=stats_zm )
          k = k + 1
        end do

      case ( 'wp2sclrp' )
        do j = 1, sclr_dim, 1
          write( sclr_idx, * ) j
          sclr_idx = adjustl(sclr_idx)
          iwp2sclrp(j) = k          
          call stat_assign( var_index=iwp2sclrp(j), var_name="wp2sclr"//trim(sclr_idx)//"p", &
            var_description="'w'^2 scalar("//trim(sclr_idx)//")", var_units="unknown", &
            l_silhs=.false., grid_kind=stats_zm )
          k = k + 1
        end do

      case ( 'wpsclrprtp' )
        do j = 1, sclr_dim, 1
          write( sclr_idx, * ) j
          sclr_idx = adjustl(sclr_idx)
          iwpsclrprtp(j) = k
          call stat_assign( var_index=iwpsclrprtp(j), var_name="wpsclr"//trim(sclr_idx)//"prtp", &
            var_description="'w' scalar("//trim(sclr_idx)//")'rt'", var_units="unknown", &
            l_silhs=.false., grid_kind=stats_zm )
          k = k + 1
        end do

      case ( 'wpsclrpthlp' )
        do j = 1, sclr_dim, 1
          write( sclr_idx, * ) j
          sclr_idx = adjustl(sclr_idx)
          iwpsclrpthlp(j) = k
          call stat_assign( var_index=iwpsclrpthlp(j), var_name="wpsclr"//trim(sclr_idx)//"pthlp", &
            var_description="'w' scalar("//trim(sclr_idx)//")'th_l'", var_units="unknown", &
            l_silhs=.false., grid_kind=stats_zm )
          k = k + 1
        end do

      case ( 'wpedsclrp' )
        do j = 1, edsclr_dim, 1
          write( sclr_idx, * ) j
          sclr_idx = adjustl(sclr_idx)
          iwpedsclrp(j) = k
          call stat_assign( var_index=iwpedsclrp(j), var_name="wpedsclr"//trim(sclr_idx)//"p", &
            var_description="eddy scalar("//trim(sclr_idx)//")'w'", var_units="unknown", &
            l_silhs=.false., grid_kind=stats_zm )
          k = k + 1
        end do
          
      case default
        write(fstderr,*) 'Error:  unrecognized variable in vars_zm:  ',  trim(vars_zm(i))
        l_error = .true.  ! This will stop the run.

      end select

    end do ! i = 1 .. stats_zm%num_output_fields

    return
  end subroutine stats_init_zm

end module stats_zm_module
