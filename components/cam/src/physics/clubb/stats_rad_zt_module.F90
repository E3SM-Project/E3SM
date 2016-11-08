!-----------------------------------------------------------------------
! $Id: stats_rad_zt_module.F90 7315 2014-09-30 20:49:54Z schemena@uwm.edu $
!===============================================================================

module stats_rad_zt_module

  implicit none

  private ! Default Scope

  public :: stats_init_rad_zt

  ! Constant parameters
  integer, parameter, public :: nvarmax_rad_zt = 250 ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_rad_zt( vars_rad_zt, l_error )

! Description:
!   Initializes array indices for stats_zt
!
! References:
!   None
!-----------------------------------------------------------------------

    use constants_clubb, only:  &
        fstderr ! Constant(s)

    use stats_variables, only: & 
        stats_rad_zt, &
        iT_in_K_rad, & ! Variable(s)
        ircil_rad, &
        io3l_rad, &
        irsm_rad, &
        ircm_in_cloud_rad, &
        icloud_frac_rad, & 
        iice_supersat_frac_rad, &
        iradht_rad, &
        iradht_LW_rad, &
        iradht_SW_rad, &
        ip_in_mb_rad, &
        isp_humidity_rad

    use stats_type_utilities, only: & 
        stat_assign ! Procedure

    implicit none

    ! Input Variable
    character(len= * ), dimension(nvarmax_rad_zt), intent(in) :: vars_rad_zt

    ! Input/Output Variable
    logical, intent(inout) :: l_error

    ! Local Varables
    integer :: i, k

    ! ---- Begin Code ----

    ! Default initialization for array indices for stats_rad_zt

    iT_in_K_rad = 0
    ircil_rad = 0
    io3l_rad = 0
    irsm_rad = 0
    ircm_in_cloud_rad = 0
    icloud_frac_rad = 0
    iice_supersat_frac_rad = 0
    iradht_rad = 0
    iradht_LW_rad = 0
    iradht_SW_rad = 0
    ip_in_mb_rad = 0
    isp_humidity_rad = 0


    ! Assign pointers for statistics variables stats_rad_zt

    k = 1
    do i=1,stats_rad_zt%num_output_fields

      select case ( trim(vars_rad_zt(i)) )

      case ('T_in_K_rad')
        iT_in_K_rad = k

        call stat_assign( var_index=iT_in_K_rad, var_name="T_in_K_rad", &
             var_description="Temperature [K]", var_units="K", l_silhs=.false., &
                grid_kind=stats_rad_zt )
        k = k + 1

      case ('rcil_rad')
        ircil_rad = k

        call stat_assign( var_index=ircil_rad, var_name="rcil_rad", &
             var_description="Ice mixing ratio [kg/kg]", var_units="kg/kg", l_silhs=.false., &
             grid_kind=stats_rad_zt )
        k = k + 1

      case ('o3l_rad')
        io3l_rad = k

        call stat_assign( var_index=io3l_rad, var_name="o3l_rad", &
             var_description="Ozone mixing ratio [kg/kg]", var_units="kg/kg", l_silhs=.false., &
             grid_kind=stats_rad_zt )
        k = k + 1

      case ('rsm_rad')
        irsm_rad = k

        call stat_assign( var_index=irsm_rad, var_name="rsm_rad", &
             var_description="Snow water mixing ratio [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_rad_zt )
        k = k + 1

      case ('rcm_in_cloud_rad')
        ircm_in_cloud_rad = k

        call stat_assign( var_index=ircm_in_cloud_rad, var_name="rcm_in_cloud_rad", &
             var_description="rcm in cloud layer [kg/kg]", var_units="kg/kg", l_silhs=.false., &
             grid_kind=stats_rad_zt )
        k = k + 1

      case ('cloud_frac_rad')
        icloud_frac_rad = k

        call stat_assign( var_index=icloud_frac_rad, var_name="cloud_frac_rad", &
             var_description="Cloud fraction (between 0 and 1) [-]", var_units="count", &
             l_silhs=.false., grid_kind=stats_rad_zt )
        k = k + 1
      
      case ('ice_supersat_frac_rad')
        iice_supersat_frac_rad = k

        call stat_assign( var_index=iice_supersat_frac_rad, var_name="ice_supersat_frac_rad", &
             var_description="Ice cloud fraction (between 0 and 1) [-]", var_units="count", &
             l_silhs=.false., grid_kind=stats_rad_zt )
        k = k + 1

      case ('radht_rad')
        iradht_rad = k

        call stat_assign( var_index=iradht_rad, var_name="radht_rad", &
             var_description="Total radiative heating rate [K/s]", var_units="K/s", &
             l_silhs=.false., grid_kind=stats_rad_zt )
        k = k + 1

      case ('radht_LW_rad')
        iradht_LW_rad = k

        call stat_assign( var_index=iradht_LW_rad, var_name="radht_LW_rad", &
             var_description="Long-wave radiative heating rate [K/s]", var_units="K/s", &
             l_silhs=.false., grid_kind=stats_rad_zt )
        k = k + 1

      case ('radht_SW_rad')
        iradht_SW_rad = k

        call stat_assign( var_index=iradht_SW_rad, var_name="radht_SW_rad", &
             var_description="Short-wave radiative heating rate [K/s]", var_units="K/s", &
             l_silhs=.false., grid_kind=stats_rad_zt )
        k = k + 1

      case ('p_in_mb_rad')
        ip_in_mb_rad = k

        call stat_assign( var_index=ip_in_mb_rad, var_name="p_in_mb_rad", &
             var_description="Pressure [hPa]", var_units="hPa", &
             l_silhs=.false., grid_kind=stats_rad_zt )
        k = k + 1

      case ('sp_humidity_rad')
        isp_humidity_rad = k

        call stat_assign( var_index=isp_humidity_rad, var_name="sp_humidity_rad", &
             var_description="Specific humidity [kg/kg]", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_rad_zt )
        k = k + 1

      case default

        write(fstderr,*) 'Error:  unrecognized variable in vars_rad_zt:  ', trim( vars_rad_zt(i) )

        l_error = .true.  ! This will stop the run.


      end select

    end do

    return
  end subroutine stats_init_rad_zt

end module stats_rad_zt_module
