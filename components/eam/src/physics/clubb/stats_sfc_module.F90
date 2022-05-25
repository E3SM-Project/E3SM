!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module stats_sfc_module


  implicit none

  private ! Set Default Scope

  public :: stats_init_sfc

  ! Constant parameters
  integer, parameter, public :: nvarmax_sfc = 250  ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_sfc( vars_sfc, l_error, & !intent(in)
                             stats_sfc ) ! intent(inout)

! Description:
!   Initializes array indices for stats_sfc
! References:
!   None
!-----------------------------------------------------------------------

    use constants_clubb, only: &
        fstderr ! Constant(s)

    use stats_variables, only: &
        iustar, &
        isoil_heat_flux, &
        iveg_T_in_K, &
        isfc_soil_T_in_K,&
        ideep_soil_T_in_K, &
        ilh, &
        ish, &
        icc, &
        ilwp, &
        ivwp, &
        iiwp, &
        iswp, &
        irwp, &
        iz_cloud_base, &
        iz_inversion, &
        iprecip_rate_sfc, &
        irain_flux_sfc, &
        irrm_sfc, &
        iprecip_frac_tol

    use stats_variables, only: &
        iwpthlp_sfc, &
        iwprtp_sfc, &
        iupwp_sfc, &
        ivpwp_sfc, &
        ithlm_vert_avg, &
        irtm_vert_avg, &
        ium_vert_avg, &
        ivm_vert_avg, &
        iwp2_vert_avg, &
        iup2_vert_avg, &
        ivp2_vert_avg, &
        irtp2_vert_avg, &
        ithlp2_vert_avg, &
        iT_sfc

    use stats_variables, only: &
        iwp23_matrix_condt_num, &
        irtm_matrix_condt_num, &
        ithlm_matrix_condt_num, &
        irtp2_matrix_condt_num, &
        ithlp2_matrix_condt_num, &
        irtpthlp_matrix_condt_num, &
        iup2_vp2_matrix_condt_num, &
        iwindm_matrix_condt_num

    use stats_variables, only: &
        imorr_snow_rate ! Variable(s)

    use stats_variables, only: &
        irtm_spur_src,            &
        ithlm_spur_src, &
        irsm_sd_morr_int
        
    use stats_variables, only: &
        itot_vartn_normlzd_rtm, &
        itot_vartn_normlzd_thlm, &
        itot_vartn_normlzd_wprtp

    use stats_type_utilities, only: &
        stat_assign ! Procedure

    use stats_type, only: stats ! Type

    implicit none

    type (stats), target, intent(inout) :: &
      stats_sfc

    ! External
    intrinsic :: trim

    ! Input Variable
    character(len= * ), dimension(nvarmax_sfc), intent(in) :: vars_sfc

    ! Input / Output Variable
    logical, intent(inout) :: l_error

    ! Local Varables
    integer :: i, k

    ! ---- Begin Code ----

    ! Default initialization for array indices for stats_sfc is zero (see module
    ! stats_variables)

    ! Assign pointers for statistics variables stats_sfc using stat_assign

    k = 1
    do i = 1, stats_sfc%num_output_fields

      select case ( trim( vars_sfc(i) ) )
      case ('soil_heat_flux')
        isoil_heat_flux = k

        call stat_assign( var_index=isoil_heat_flux, var_name="soil_heat_flux", &
             var_description="soil_heat_flux, soil_heat_flux", &
             var_units="W/m^2", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1
      case ('ustar')
        iustar = k

        call stat_assign( var_index=iustar, var_name="ustar", &
             var_description="ustar, Friction velocity", var_units="m/s", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1
      case ('veg_T_in_K')
        iveg_T_in_K = k

        call stat_assign( var_index=iveg_T_in_K, var_name="veg_T_in_K", &
             var_description="veg_T_in_K, Surface Vegetation Temperature", var_units="K", &
             l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1
      case ('sfc_soil_T_in_K')
        isfc_soil_T_in_K = k

        call stat_assign( var_index=isfc_soil_T_in_K, var_name="sfc_soil_T_in_K", &
             var_description="sfc_soil_T_in_K, Surface soil temperature", &
             var_units="K", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1
      case ('deep_soil_T_in_K')
        ideep_soil_T_in_K = k

        call stat_assign( var_index=ideep_soil_T_in_K, var_name="deep_soil_T_in_K", &
             var_description="deep_soil_T_in_K, Deep soil Temperature", &
             var_units="K", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('lh')
        ilh = k
        call stat_assign( var_index=ilh, var_name="lh", &
             var_description="lh, Surface latent heating", var_units="W/m2", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('sh')
        ish = k
        call stat_assign( var_index=ish, var_name="sh", &
             var_description="sh, Surface sensible heating", var_units="W/m2", &
             l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('cc')
        icc = k
        call stat_assign( var_index=icc, var_name="cc", var_description="cc, Cloud cover", &
             var_units="count", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('lwp')
        ilwp = k
        call stat_assign( var_index=ilwp, var_name="lwp", &
             var_description="lwp, Liquid water path", var_units="kg/m2", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('vwp')
        ivwp = k
        call stat_assign( var_index=ivwp, var_name="vwp", &
             var_description="vwp, Vapor water path", var_units="kg/m2", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('iwp')
        iiwp = k
        call stat_assign( var_index=iiwp, var_name="iwp", &
             var_description="iwp, Ice water path", var_units="kg/m2", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('swp')
        iswp = k
        call stat_assign( var_index=iswp, var_name="swp", &
             var_description="swp, Snow water path", var_units="kg/m2", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('rwp')
        irwp = k
        call stat_assign( var_index=irwp, var_name="rwp", &
             var_description="rwp, Rain water path", var_units="kg/m2", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('z_cloud_base')
        iz_cloud_base = k
        call stat_assign( var_index=iz_cloud_base, var_name="z_cloud_base", &
             var_description="z_cloud_base, Cloud base altitude", &
             var_units="m", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('z_inversion')
        iz_inversion = k
        call stat_assign( var_index=iz_inversion, var_name="z_inversion", &
             var_description="z_inversion, Inversion altitude", var_units="m", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('precip_rate_sfc')          ! Brian
        iprecip_rate_sfc = k
        call stat_assign( var_index=iprecip_rate_sfc, var_name="precip_rate_sfc", &
             var_description="precip_rate_sfc, Surface rainfall rate", var_units="mm/day", &
             l_silhs=.true., grid_kind=stats_sfc )
        k = k + 1

      case ('rain_flux_sfc')         ! Brian
        irain_flux_sfc = k

        call stat_assign( var_index=irain_flux_sfc, var_name="rain_flux_sfc", &
             var_description="rain_flux_sfc, Surface rain flux", &
             var_units="W/m^2", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('rrm_sfc')       ! Brian
        irrm_sfc = k

        call stat_assign( var_index=irrm_sfc, var_name="rrm_sfc", &
             var_description="rrm_sfc, Surface rain water mixing ratio", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('precip_frac_tol')
        iprecip_frac_tol = k

        call stat_assign( var_index=iprecip_frac_tol, &
                          var_name="precip_frac_tol", &
                          var_description="Smallest allowable precipitation " &
                          // "fraction when hydrometeors are present [-]", &
                          var_units="-", &
                          l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ( 'morr_snow_rate' )
        imorr_snow_rate = k
        call stat_assign( var_index=imorr_snow_rate, var_name="morr_snow_rate", &
             var_description="morr_snow_rate, Snow+Ice+Graupel fallout rate " &
             // "from Morrison scheme", &
             var_units="mm/day", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('wpthlp_sfc')
        iwpthlp_sfc = k

        call stat_assign( var_index=iwpthlp_sfc, var_name="wpthlp_sfc", &
             var_description="w'thl'_sfc, wpthlp surface flux", &
             var_units="K m/s", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('wprtp_sfc')
        iwprtp_sfc = k

        call stat_assign( var_index=iwprtp_sfc, var_name="wprtp_sfc", &
             var_description="w'rt'_sfc, wprtp surface flux", var_units="(kg/kg) m/s", &
             l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('upwp_sfc')
        iupwp_sfc = k

        call stat_assign( var_index=iupwp_sfc, var_name="upwp_sfc", &
             var_description="u'w'_sfc, upwp surface flux", &
             var_units="m^2/s^2", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('vpwp_sfc')
        ivpwp_sfc = k

        call stat_assign( var_index=ivpwp_sfc, var_name="vpwp_sfc", &
             var_description="v'w'_sfc, vpwp surface flux", &
             var_units="m^2/s^2", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('thlm_vert_avg')
        ithlm_vert_avg = k

        call stat_assign( var_index=ithlm_vert_avg, var_name="thlm_vert_avg", &
             var_description="thlm_vert_avg, Vertical average (density-weighted) of thlm", &
             var_units="K", &
             l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('rtm_vert_avg')
        irtm_vert_avg = k

        call stat_assign( var_index=irtm_vert_avg, var_name="rtm_vert_avg", &
             var_description="rtm_vert_avg, Vertical average (density-weighted) of rtm", &
             var_units="kg/kg", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('um_vert_avg')
        ium_vert_avg = k

        call stat_assign( var_index=ium_vert_avg, var_name="um_vert_avg", &
             var_description="um_vert_avg, Vertical average (density-weighted) of um", &
             var_units="m/s", &
             l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('vm_vert_avg')
        ivm_vert_avg = k

        call stat_assign( var_index=ivm_vert_avg, var_name="vm_vert_avg", &
             var_description="vm_vert_avg, Vertical average (density-weighted) of vm", &
             var_units="m/s", &
             l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('wp2_vert_avg')
        iwp2_vert_avg = k

        call stat_assign( var_index=iwp2_vert_avg, var_name="wp2_vert_avg", &
             var_description="Density-weighted vertical average of w'^2", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('up2_vert_avg')
        iup2_vert_avg = k

        call stat_assign( var_index=iup2_vert_avg, var_name="up2_vert_avg", &
             var_description="u'^2_vert_avg, Vertical average (density-weighted) of up2", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('vp2_vert_avg')
        ivp2_vert_avg = k

        call stat_assign( var_index=ivp2_vert_avg, var_name="vp2_vert_avg", &
             var_description="v'^2_vert_avg, Vertical average (density-weighted) of vp2", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('rtp2_vert_avg')
        irtp2_vert_avg = k

        call stat_assign( var_index=irtp2_vert_avg, var_name="rtp2_vert_avg", &
             var_description="rt'^2_vert_avg, Vertical average (density-weighted) of rtp2", &
             var_units="kg^2/kg^2", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('thlp2_vert_avg')
        ithlp2_vert_avg = k

        call stat_assign( var_index=ithlp2_vert_avg, var_name="thlp2_vert_avg", &
             var_description="thl'^2_vert_avg, Vertical average (density-weighted) of thlp2", &
             var_units="K^2", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('T_sfc')
        iT_sfc = k

        call stat_assign( var_index=iT_sfc, var_name="T_sfc", &
             var_description="T_sfc, Surface Temperature", var_units="K", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('wp23_matrix_condt_num')
        iwp23_matrix_condt_num = k
        call stat_assign( var_index=iwp23_matrix_condt_num, var_name="wp23_matrix_condt_num", &
             var_description="w'^23_matrix_condt_num, Estimate of the condition number for wp2/3", &
             var_units="count", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('thlm_matrix_condt_num')
        ithlm_matrix_condt_num = k
        call stat_assign( var_index=ithlm_matrix_condt_num, var_name="thlm_matrix_condt_num", &
             var_description="thlm_matrix_condt_num, Estimate of the condition " &
             // "number for thlm/wpthlp", &
             var_units="count", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('rtm_matrix_condt_num')
        irtm_matrix_condt_num = k

        call stat_assign( var_index=irtm_matrix_condt_num, var_name="rtm_matrix_condt_num", &
             var_description="rtm_matrix_condt_num, Estimate of the condition number " &
             // "for rtm/wprtp", &
             var_units="count", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('thlp2_matrix_condt_num')
        ithlp2_matrix_condt_num = k

        call stat_assign( var_index=ithlp2_matrix_condt_num, var_name="thlp2_matrix_condt_num", &
             var_description="thl'^2_matrix_condt_num, Estimate of the condition " &
             // "number for thlp2", &
             var_units="count", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('rtp2_matrix_condt_num')
        irtp2_matrix_condt_num = k
        call stat_assign( var_index=irtp2_matrix_condt_num, var_name="rtp2_matrix_condt_num", &
             var_description="rt'^2_matrix_condt_num, Estimate of the condition " &
             // "number for rtp2", &
             var_units="count", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('rtpthlp_matrix_condt_num')
        irtpthlp_matrix_condt_num = k
        call stat_assign( var_index=irtpthlp_matrix_condt_num, &
             var_name="rtpthlp_matrix_condt_num", &
             var_description="rt'thl'_matrix_condt_num, Estimate of the condition " &
             // "number for rtpthlp", &
             var_units="count", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('up2_vp2_matrix_condt_num')
        iup2_vp2_matrix_condt_num = k
        call stat_assign( var_index=iup2_vp2_matrix_condt_num, &
             var_name="up2_vp2_matrix_condt_num", &
             var_description="u'^2_v'^2_matrix_condt_num, Estimate of the condition " &
             // "number for up2/vp2", &
             var_units="count", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('windm_matrix_condt_num')
        iwindm_matrix_condt_num = k
        call stat_assign( var_index=iwindm_matrix_condt_num, var_name="windm_matrix_condt_num", &
             var_description="windm_matrix_condt_num, Estimate of the condition " &
             // "number for the mean wind", &
             var_units="count", l_silhs=.false., grid_kind=stats_sfc )

        k = k + 1

      case ('rtm_spur_src')
        irtm_spur_src = k

        call stat_assign( var_index=irtm_spur_src, var_name="rtm_spur_src", &
             var_description="rtm_spur_src, rtm spurious source", var_units="kg/(m^2 s)", &
             l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('thlm_spur_src')
        ithlm_spur_src = k

        call stat_assign( var_index=ithlm_spur_src, var_name="thlm_spur_src", &
             var_description="thlm_spur_src, thlm spurious source", &
             var_units="(K kg) / (m^2 s)", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('rs_sd_morr_int')
        irsm_sd_morr_int = k

        call stat_assign( var_index=irsm_sd_morr_int, var_name="rs_sd_morr_int", &
             var_description="rs_sd_morr_int, rsm_sd_morr vertical integral", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_sfc )
        k = k + 1
        
      case ('tot_vartn_normlzd_rtm')
        itot_vartn_normlzd_rtm = k

        call stat_assign( var_index=itot_vartn_normlzd_rtm, var_name="tot_vartn_normlzd_rtm", &
             var_description="Total variation of rtm in the vertical normalized", &
             var_units="-", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1
        
      case ('tot_vartn_normlzd_thlm')
        itot_vartn_normlzd_thlm = k

        call stat_assign( var_index=itot_vartn_normlzd_thlm, var_name="tot_vartn_normlzd_thlm", &
             var_description="Total variation of thlm in the vertical normalized", &
             var_units="-", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1
        
      case ('tot_vartn_normlzd_wprtp')
        itot_vartn_normlzd_wprtp = k

        call stat_assign( var_index=itot_vartn_normlzd_wprtp, var_name="tot_vartn_normlzd_wprtp", &
             var_description="Total variation of wprtp in the vertical normalized", &
             var_units="-", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case default
        write(fstderr,*) 'Error:  unrecognized variable in vars_sfc:  ',  &
              trim( vars_sfc(i) )
        l_error = .true.  ! This will stop the run.

      end select

    end do ! 1 .. stats_sfc%num_output_fields

    return

  end subroutine stats_init_sfc


end module stats_sfc_module
