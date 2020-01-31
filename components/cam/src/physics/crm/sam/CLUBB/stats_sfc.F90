!-----------------------------------------------------------------------
! $Id: stats_sfc.F90 6146 2013-04-05 18:02:22Z raut@uwm.edu $

module stats_sfc


  implicit none

  private ! Set Default Scope

  public :: stats_init_sfc

  ! Constant parameters
  integer, parameter, public :: nvarmax_sfc = 250  ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_sfc( vars_sfc, l_error )

! Description:
!   Initializes array indices for sfc
! References:
!   None
!-----------------------------------------------------------------------

    use constants_clubb, only: &
        fstderr ! Constant(s)

    use stats_variables, only: & 
        sfc,  & ! Variables
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
        irain_rate_sfc, & 
        irain_flux_sfc, & 
        irrainm_sfc

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
      imorr_rain_rate, &
      imorr_snow_rate
      
    use stats_variables, only: &
      irtm_spur_src,            &
      ithlm_spur_src

    use stats_type, only: & 
        stat_assign ! Procedure

    implicit none

    ! Input Variable
    character(len= * ), dimension(nvarmax_sfc), intent(in) :: vars_sfc

    ! Output Variable        
    logical, intent(inout) :: l_error

    ! Local Varables
    integer :: i, k

    ! ---- Begin Code ----

    ! Default initialization for array indices for sfc

    isoil_heat_flux   = 0
    iveg_T_in_K       = 0
    isfc_soil_T_in_K  = 0
    ideep_soil_T_in_K = 0

    iustar          = 0
    ilh             = 0
    ish             = 0
    icc             = 0
    ilwp            = 0
    irwp            = 0
    ivwp            = 0   ! nielsenb
    iiwp            = 0   ! nielsenb
    iswp            = 0   ! nielsenb
    iz_cloud_base   = 0
    iz_inversion             = 0
    irain_rate_sfc  = 0   ! Brian
    irain_flux_sfc  = 0   ! Brian
    irrainm_sfc     = 0   ! Brian
    iwpthlp_sfc     = 0
    iwprtp_sfc      = 0
    iupwp_sfc       = 0
    ivpwp_sfc       = 0
    ithlm_vert_avg  = 0
    irtm_vert_avg   = 0
    ium_vert_avg    = 0
    ivm_vert_avg    = 0
    iwp2_vert_avg   = 0   ! nielsenb
    iup2_vert_avg   = 0
    ivp2_vert_avg   = 0
    irtp2_vert_avg  = 0
    ithlp2_vert_avg = 0
    iT_sfc            = 0   ! kcwhite
    
    ! These are estimates of the condition number on each LHS
    ! matrix, and not located at the surface of the domain.
    iwp23_matrix_condt_num    = 0
    irtm_matrix_condt_num     = 0
    ithlm_matrix_condt_num    = 0
    irtp2_matrix_condt_num    = 0
    ithlp2_matrix_condt_num   = 0
    irtpthlp_matrix_condt_num = 0
    iup2_vp2_matrix_condt_num = 0
    iwindm_matrix_condt_num   = 0

    imorr_rain_rate = 0
    imorr_snow_rate = 0

    irtm_spur_src = 0
    ithlm_spur_src = 0

    ! Assign pointers for statistics variables sfc

    k = 1
    do i=1,sfc%nn

      select case ( trim(vars_sfc(i)) )
      case ('soil_heat_flux')
        isoil_heat_flux = k

        call stat_assign(isoil_heat_flux, "soil_heat_flux", & 
             "soil_heat_flux[W/m^2]","W/m^2",sfc )
        k = k + 1
      case ('ustar')
        iustar = k

        call stat_assign(iustar,"ustar", & 
             "Friction velocity [m/s]","m/s",sfc)
        k = k + 1
      case ('veg_T_in_K')
        iveg_T_in_K = k

        call stat_assign(iveg_T_in_K,"veg_T_in_K", & 
             "Surface Vegetation Temperature [K]","K",sfc)
        k = k + 1
      case ('sfc_soil_T_in_K')
        isfc_soil_T_in_K = k

        call stat_assign(isfc_soil_T_in_K,"sfc_soil_T_in_K", & 
             "Surface soil temperature [K]","K",sfc)
        k = k + 1
      case ('deep_soil_T_in_K')
        ideep_soil_T_in_K = k

        call stat_assign(ideep_soil_T_in_K,"deep_soil_T_in_K", & 
             "Deep soil Temperature [K]","K",sfc)
        k = k + 1

      case ('lh')
        ilh = k
        call stat_assign(ilh,"lh", & 
             "Surface latent heating [W/m^2]","W/m2",sfc)
        k = k + 1

      case ('sh')
        ish = k
        call stat_assign(ish,"sh", & 
             "Surface sensible heating [W/m^2]","W/m2",sfc)
        k = k + 1

      case ('cc')
        icc = k
        call stat_assign(icc,"cc", & 
             "Cloud cover [count]","count",sfc)
        k = k + 1

      case ('lwp')
        ilwp = k
        call stat_assign(ilwp,"lwp", & 
             "Liquid water path [kg/m^2]","kg/m2",sfc)
        k = k + 1

      case ('vwp')
        ivwp = k
        call stat_assign(ivwp,"vwp", & 
             "Vapor water path [kg/m^2]","kg/m2",sfc)
        k = k + 1

      case ('iwp')
        iiwp = k
        call stat_assign(iiwp,"iwp", & 
             "Ice water path [kg/m^2]","kg/m2",sfc)
        k = k + 1

      case ('swp')
        iswp = k
        call stat_assign(iswp,"swp", & 
             "Snow water path [kg/m^2]","kg/m2",sfc)
        k = k + 1

      case ('rwp')
        irwp = k
        call stat_assign(irwp,"rwp", & 
             "Rain water path [kg/m^2]","kg/m2",sfc)
        k = k + 1

      case ('z_cloud_base')
        iz_cloud_base = k
        call stat_assign(iz_cloud_base,"z_cloud_base", & 
             "Cloud base altitude [m]","m",sfc)
        k = k + 1

      case ('z_inversion')
        iz_inversion = k
        call stat_assign(iz_inversion,"z_inversion", & 
             "Inversion altitude [m]","m",sfc)
        k = k + 1

      case ('rain_rate_sfc')          ! Brian
        irain_rate_sfc = k
        call stat_assign(irain_rate_sfc,"rain_rate_sfc", & 
             "Surface rainfall rate [mm/day]","mm/day",sfc)
        k = k + 1

      case ('rain_flux_sfc')         ! Brian
        irain_flux_sfc = k

        call stat_assign( irain_flux_sfc,"rain_flux_sfc", & 
             "Surface rain flux [W/m^2]", "W/m^2", sfc )
        k = k + 1

      case ('rrainm_sfc')       ! Brian
        irrainm_sfc = k

        call stat_assign(irrainm_sfc,"rrainm_sfc", & 
             "Surface rain water mixing ratio [kg/kg]","kg/kg",sfc)
        k = k + 1

      case ( 'morr_rain_rate' )
        imorr_rain_rate = k
        call stat_assign( imorr_rain_rate, "morr_rain_rate", & 
             "Total precip fallout rate from Morrison scheme [mm/day]","mm/day", sfc )
        k = k + 1

      case ( 'morr_snow_rate' )
        imorr_snow_rate = k
        call stat_assign( imorr_snow_rate, "morr_snow_rate", & 
             "Snow+Ice+Graupel fallout rate from Morrison scheme [mm/day]","mm/day", sfc )
        k = k + 1

      case ('wpthlp_sfc')
        iwpthlp_sfc = k

        call stat_assign(iwpthlp_sfc,"wpthlp_sfc", &
             "wpthlp surface flux [K m/s]","K m/s",sfc)
        k = k + 1

      case ('wprtp_sfc')
        iwprtp_sfc = k

        call stat_assign(iwprtp_sfc,"wprtp_sfc", &
             "wprtp surface flux [kg/kg]","(kg/kg) m/s",sfc)
        k = k + 1

      case ('upwp_sfc')
        iupwp_sfc = k

        call stat_assign(iupwp_sfc,"upwp_sfc", &
             "upwp surface flux [m^2/s^2]","m^2/s^2",sfc)
        k = k + 1

      case ('vpwp_sfc')
        ivpwp_sfc = k

        call stat_assign(ivpwp_sfc,"vpwp_sfc", &
             "vpwp surface flux [m^2/s^2]","m^2/s^2",sfc)
        k = k + 1

      case ('thlm_vert_avg')
        ithlm_vert_avg = k

        call stat_assign( ithlm_vert_avg, "thlm_vert_avg", &
             "Vertical average (density-weighted) of thlm [K]", "K", sfc )
        k = k + 1

      case ('rtm_vert_avg')
        irtm_vert_avg = k

        call stat_assign( irtm_vert_avg, "rtm_vert_avg", &
             "Vertical average (density-weighted) of rtm [kg/kg]", "kg/kg", sfc )
        k = k + 1

      case ('um_vert_avg')
        ium_vert_avg = k

        call stat_assign( ium_vert_avg, "um_vert_avg", &
             "Vertical average (density-weighted) of um [m/s]", "m/s", sfc )
        k = k + 1

      case ('vm_vert_avg')
        ivm_vert_avg = k

        call stat_assign( ivm_vert_avg, "vm_vert_avg", &
             "Vertical average (density-weighted) of vm [m/s]", "m/s", sfc )
        k = k + 1

      case ('wp2_vert_avg')
        iwp2_vert_avg = k

        call stat_assign( iwp2_vert_avg, "wp2_vert_avg", &
             "Vertical average (density-weighted) of wp2 [m^2/s^2]", "m^2/s^2", &
                          sfc )
        k = k + 1

      case ('up2_vert_avg')
        iup2_vert_avg = k

        call stat_assign( iup2_vert_avg, "up2_vert_avg", &
             "Vertical average (density-weighted) of up2 [m^2/s^2]", "m^2/s^2", &
                          sfc )
        k = k + 1

      case ('vp2_vert_avg')
        ivp2_vert_avg = k

        call stat_assign( ivp2_vert_avg, "vp2_vert_avg", &
             "Vertical average (density-weighted) of vp2 [m^2/s^2]", "m^2/s^2", &
                          sfc )
        k = k + 1

      case ('rtp2_vert_avg')
        irtp2_vert_avg = k

        call stat_assign( irtp2_vert_avg, "rtp2_vert_avg", &
             "Vertical average (density-weighted) of rtp2 [kg^2/kg^2]", &
                          "kg^2/kg^2", sfc )
        k = k + 1

      case ('thlp2_vert_avg')
        ithlp2_vert_avg = k

        call stat_assign( ithlp2_vert_avg, "thlp2_vert_avg", &
             "Vertical average (density-weighted) of thlp2 [K^2]", "K^2", sfc )
        k = k + 1

      case ('T_sfc')
        iT_sfc = k

        call stat_assign( iT_sfc, "T_sfc", "Surface Temperature [K]", "K", sfc )
        k = k + 1 

      case ('wp23_matrix_condt_num')
        iwp23_matrix_condt_num = k
        call stat_assign(iwp23_matrix_condt_num,"wp23_matrix_condt_num", & 
             "Estimate of the condition number for wp2/3 [count]","count",sfc)
        k = k + 1

      case ('thlm_matrix_condt_num')
        ithlm_matrix_condt_num = k
        call stat_assign(ithlm_matrix_condt_num,"thlm_matrix_condt_num", & 
             "Estimate of the condition number for thlm/wpthlp [count]", & 
             "count",sfc)
        k = k + 1

      case ('rtm_matrix_condt_num')
        irtm_matrix_condt_num = k

        call stat_assign(irtm_matrix_condt_num,"rtm_matrix_condt_num", & 
             "Estimate of the condition number for rtm/wprtp [count]", & 
             "count",sfc)
        k = k + 1

      case ('thlp2_matrix_condt_num')
        ithlp2_matrix_condt_num = k

        call stat_assign(ithlp2_matrix_condt_num,"thlp2_matrix_condt_num", & 
             "Estimate of the condition number for thlp2 [count]", & 
             "count",sfc)
        k = k + 1

      case ('rtp2_matrix_condt_num')
        irtp2_matrix_condt_num = k
        call stat_assign(irtp2_matrix_condt_num,"rtp2_matrix_condt_num", & 
             "Estimate of the condition number for rtp2 [count]", & 
             "count",sfc)
        k = k + 1

      case ('rtpthlp_matrix_condt_num')
        irtpthlp_matrix_condt_num = k
        call stat_assign(irtpthlp_matrix_condt_num,"rtpthlp_matrix_condt_num", & 
            "Estimate of the condition number for rtpthlp [count]", & 
            "count",sfc)
        k = k + 1

      case ('up2_vp2_matrix_condt_num')
        iup2_vp2_matrix_condt_num = k
        call stat_assign(iup2_vp2_matrix_condt_num,"up2_vp2_matrix_condt_num", & 
             "Estimate of the condition number for up2/vp2 [count]","count",sfc)
        k = k + 1

      case ('windm_matrix_condt_num')
        iwindm_matrix_condt_num = k
        call stat_assign(iwindm_matrix_condt_num,"windm_matrix_condt_num", & 
             "Estimate of the condition number for the mean wind [count]","count",sfc)

        k = k + 1
        
      case ('rtm_spur_src')
        irtm_spur_src = k

        call stat_assign(irtm_spur_src, "rtm_spur_src", & 
             "rtm spurious source [kg/(m^2 s)]", "kg/(m^2 s)",sfc )
        k = k + 1
        
      case ('thlm_spur_src')
        ithlm_spur_src = k

        call stat_assign(ithlm_spur_src, "thlm_spur_src", & 
             "thlm spurious source [(K kg) / (m^2 s)]", "(K kg) / (m^2 s)",sfc )
        k = k + 1

      case default
        write(fstderr,*) 'Error:  unrecognized variable in vars_sfc:  ',  &
              trim( vars_sfc(i) )
        l_error = .true.  ! This will stop the run.

      end select

    end do

    return

  end subroutine stats_init_sfc


end module stats_sfc

