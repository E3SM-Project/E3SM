!-----------------------------------------------------------------------
! $Id: stats_sfc_module.F90 7315 2014-09-30 20:49:54Z schemena@uwm.edu $
!===============================================================================
module stats_sfc_module


  implicit none

  private ! Set Default Scope

  public :: stats_init_sfc

  ! Constant parameters
  integer, parameter, public :: nvarmax_sfc = 250  ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_sfc( vars_sfc, l_error )

! Description:
!   Initializes array indices for stats_sfc
! References:
!   None
!-----------------------------------------------------------------------

    use constants_clubb, only: &
        fstderr ! Constant(s)

    use stats_variables, only: & 
        stats_sfc,  & ! Variables
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
        irrm_sfc

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

    use stats_type_utilities, only: & 
        stat_assign ! Procedure

    implicit none

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
             var_description="soil_heat_flux[W/m^2]", var_units="W/m^2", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1
      case ('ustar')
        iustar = k

        call stat_assign( var_index=iustar, var_name="ustar", &
             var_description="Friction velocity [m/s]", var_units="m/s", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1
      case ('veg_T_in_K')
        iveg_T_in_K = k

        call stat_assign( var_index=iveg_T_in_K, var_name="veg_T_in_K", &
             var_description="Surface Vegetation Temperature [K]", var_units="K", &
             l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1
      case ('sfc_soil_T_in_K')
        isfc_soil_T_in_K = k

        call stat_assign( var_index=isfc_soil_T_in_K, var_name="sfc_soil_T_in_K", &
             var_description="Surface soil temperature [K]", var_units="K", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1
      case ('deep_soil_T_in_K')
        ideep_soil_T_in_K = k

        call stat_assign( var_index=ideep_soil_T_in_K, var_name="deep_soil_T_in_K", &
             var_description="Deep soil Temperature [K]", var_units="K", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('lh')
        ilh = k
        call stat_assign( var_index=ilh, var_name="lh", &
             var_description="Surface latent heating [W/m^2]", var_units="W/m2", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('sh')
        ish = k
        call stat_assign( var_index=ish, var_name="sh", &
             var_description="Surface sensible heating [W/m^2]", var_units="W/m2", &
             l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('cc')
        icc = k
        call stat_assign( var_index=icc, var_name="cc", var_description="Cloud cover [count]", &
             var_units="count", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('lwp')
        ilwp = k
        call stat_assign( var_index=ilwp, var_name="lwp", &
             var_description="Liquid water path [kg/m^2]", var_units="kg/m2", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('vwp')
        ivwp = k
        call stat_assign( var_index=ivwp, var_name="vwp", &
             var_description="Vapor water path [kg/m^2]", var_units="kg/m2", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('iwp')
        iiwp = k
        call stat_assign( var_index=iiwp, var_name="iwp", &
             var_description="Ice water path [kg/m^2]", var_units="kg/m2", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('swp')
        iswp = k
        call stat_assign( var_index=iswp, var_name="swp", &
             var_description="Snow water path [kg/m^2]", var_units="kg/m2", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('rwp')
        irwp = k
        call stat_assign( var_index=irwp, var_name="rwp", &
             var_description="Rain water path [kg/m^2]", var_units="kg/m2", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('z_cloud_base')
        iz_cloud_base = k
        call stat_assign( var_index=iz_cloud_base, var_name="z_cloud_base", &
             var_description="Cloud base altitude [m]", var_units="m", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('z_inversion')
        iz_inversion = k
        call stat_assign( var_index=iz_inversion, var_name="z_inversion", &
             var_description="Inversion altitude [m]", var_units="m", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('precip_rate_sfc')          ! Brian
        iprecip_rate_sfc = k
        call stat_assign( var_index=iprecip_rate_sfc, var_name="precip_rate_sfc", &
             var_description="Surface rainfall rate [mm/day]", var_units="mm/day", &
             l_silhs=.true., grid_kind=stats_sfc )
        k = k + 1

      case ('rain_flux_sfc')         ! Brian
        irain_flux_sfc = k

        call stat_assign( var_index=irain_flux_sfc, var_name="rain_flux_sfc", &
             var_description="Surface rain flux [W/m^2]", var_units="W/m^2", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('rrm_sfc')       ! Brian
        irrm_sfc = k

        call stat_assign( var_index=irrm_sfc, var_name="rrm_sfc", &
             var_description="Surface rain water mixing ratio [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ( 'morr_snow_rate' )
        imorr_snow_rate = k
        call stat_assign( var_index=imorr_snow_rate, var_name="morr_snow_rate", &
             var_description="Snow+Ice+Graupel fallout rate from Morrison scheme [mm/day]", &
             var_units="mm/day", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('wpthlp_sfc')
        iwpthlp_sfc = k

        call stat_assign( var_index=iwpthlp_sfc, var_name="wpthlp_sfc", &
             var_description="wpthlp surface flux [K m/s]", var_units="K m/s", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('wprtp_sfc')
        iwprtp_sfc = k

        call stat_assign( var_index=iwprtp_sfc, var_name="wprtp_sfc", &
             var_description="wprtp surface flux [kg/kg]", var_units="(kg/kg) m/s", &
             l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('upwp_sfc')
        iupwp_sfc = k

        call stat_assign( var_index=iupwp_sfc, var_name="upwp_sfc", &
             var_description="upwp surface flux [m^2/s^2]", var_units="m^2/s^2", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('vpwp_sfc')
        ivpwp_sfc = k

        call stat_assign( var_index=ivpwp_sfc, var_name="vpwp_sfc", &
             var_description="vpwp surface flux [m^2/s^2]", var_units="m^2/s^2", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1

      case ('thlm_vert_avg')
        ithlm_vert_avg = k

        call stat_assign( var_index=ithlm_vert_avg, var_name="thlm_vert_avg", &
             var_description="Vertical average (density-weighted) of thlm [K]", var_units="K", &
             l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('rtm_vert_avg')
        irtm_vert_avg = k

        call stat_assign( var_index=irtm_vert_avg, var_name="rtm_vert_avg", &
             var_description="Vertical average (density-weighted) of rtm [kg/kg]", &
             var_units="kg/kg", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('um_vert_avg')
        ium_vert_avg = k

        call stat_assign( var_index=ium_vert_avg, var_name="um_vert_avg", &
             var_description="Vertical average (density-weighted) of um [m/s]", var_units="m/s", &
             l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('vm_vert_avg')
        ivm_vert_avg = k

        call stat_assign( var_index=ivm_vert_avg, var_name="vm_vert_avg", &
             var_description="Vertical average (density-weighted) of vm [m/s]", var_units="m/s", &
             l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('wp2_vert_avg')
        iwp2_vert_avg = k

        call stat_assign( var_index=iwp2_vert_avg, var_name="wp2_vert_avg", &
             var_description="Vertical average (density-weighted) of wp2 [m^2/s^2]", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('up2_vert_avg')
        iup2_vert_avg = k

        call stat_assign( var_index=iup2_vert_avg, var_name="up2_vert_avg", &
             var_description="Vertical average (density-weighted) of up2 [m^2/s^2]", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('vp2_vert_avg')
        ivp2_vert_avg = k

        call stat_assign( var_index=ivp2_vert_avg, var_name="vp2_vert_avg", &
             var_description="Vertical average (density-weighted) of vp2 [m^2/s^2]", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('rtp2_vert_avg')
        irtp2_vert_avg = k

        call stat_assign( var_index=irtp2_vert_avg, var_name="rtp2_vert_avg", &
             var_description="Vertical average (density-weighted) of rtp2 [kg^2/kg^2]", &
             var_units="kg^2/kg^2", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('thlp2_vert_avg')
        ithlp2_vert_avg = k

        call stat_assign( var_index=ithlp2_vert_avg, var_name="thlp2_vert_avg", &
             var_description="Vertical average (density-weighted) of thlp2 [K^2]", &
             var_units="K^2", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('T_sfc')
        iT_sfc = k

        call stat_assign( var_index=iT_sfc, var_name="T_sfc", &
             var_description="Surface Temperature [K]", var_units="K", l_silhs=.false., &
             grid_kind=stats_sfc )
        k = k + 1 

      case ('wp23_matrix_condt_num')
        iwp23_matrix_condt_num = k
        call stat_assign( var_index=iwp23_matrix_condt_num, var_name="wp23_matrix_condt_num", &
             var_description="Estimate of the condition number for wp2/3 [count]", &
             var_units="count", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('thlm_matrix_condt_num')
        ithlm_matrix_condt_num = k
        call stat_assign( var_index=ithlm_matrix_condt_num, var_name="thlm_matrix_condt_num", &
             var_description="Estimate of the condition number for thlm/wpthlp [count]", &
             var_units="count", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('rtm_matrix_condt_num')
        irtm_matrix_condt_num = k

        call stat_assign( var_index=irtm_matrix_condt_num, var_name="rtm_matrix_condt_num", &
             var_description="Estimate of the condition number for rtm/wprtp [count]", &
             var_units="count", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('thlp2_matrix_condt_num')
        ithlp2_matrix_condt_num = k

        call stat_assign( var_index=ithlp2_matrix_condt_num, var_name="thlp2_matrix_condt_num", &
             var_description="Estimate of the condition number for thlp2 [count]", &
             var_units="count", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('rtp2_matrix_condt_num')
        irtp2_matrix_condt_num = k
        call stat_assign( var_index=irtp2_matrix_condt_num, var_name="rtp2_matrix_condt_num", &
             var_description="Estimate of the condition number for rtp2 [count]", &
             var_units="count", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('rtpthlp_matrix_condt_num')
        irtpthlp_matrix_condt_num = k
        call stat_assign( var_index=irtpthlp_matrix_condt_num, &
             var_name="rtpthlp_matrix_condt_num", &
             var_description="Estimate of the condition number for rtpthlp [count]", &
             var_units="count", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('up2_vp2_matrix_condt_num')
        iup2_vp2_matrix_condt_num = k
        call stat_assign( var_index=iup2_vp2_matrix_condt_num, &
             var_name="up2_vp2_matrix_condt_num", &
             var_description="Estimate of the condition number for up2/vp2 [count]", &
             var_units="count", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('windm_matrix_condt_num')
        iwindm_matrix_condt_num = k
        call stat_assign( var_index=iwindm_matrix_condt_num, var_name="windm_matrix_condt_num", &
             var_description="Estimate of the condition number for the mean wind [count]", &
             var_units="count", l_silhs=.false., grid_kind=stats_sfc )

        k = k + 1
        
      case ('rtm_spur_src')
        irtm_spur_src = k

        call stat_assign( var_index=irtm_spur_src, var_name="rtm_spur_src", &
             var_description="rtm spurious source [kg/(m^2 s)]", var_units="kg/(m^2 s)", &
             l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1
        
      case ('thlm_spur_src')
        ithlm_spur_src = k

        call stat_assign( var_index=ithlm_spur_src, var_name="thlm_spur_src", &
             var_description="thlm spurious source [(K kg) / (m^2 s)]", &
             var_units="(K kg) / (m^2 s)", l_silhs=.false., grid_kind=stats_sfc )
        k = k + 1

      case ('rs_sd_morr_int')
        irsm_sd_morr_int = k

        call stat_assign( var_index=irsm_sd_morr_int, var_name="rs_sd_morr_int", &
             var_description="rsm_sd_morr vertical integral [(kg/kg)/s]", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_sfc )
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

