!-----------------------------------------------------------------------
! $Id: stats_rad_zt.F90 4032 2009-08-17 21:45:29Z senkbeil@uwm.edu $

module stats_rad_zt

  implicit none

  private ! Default Scope

  public :: stats_init_rad_zt

  ! Constant parameters
  integer, parameter, public :: nvarmax_rad_zt = 250 ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_rad_zt( vars_rad_zt, l_error )

! Description:
!   Initializes array indices for zt
!
! References:
!   None
!-----------------------------------------------------------------------

    use constants_clubb, only:  &
        fstderr ! Constant(s)

    use stats_variables, only: & 
        rad_zt, &
        iT_in_K_rad, & ! Variable(s)
        ircil_rad, &
        io3l_rad, &
        irsnowm_rad, &
        ircm_in_cloud_rad, &
        icloud_frac_rad, & 
        iice_supersat_frac_rad, &
        iradht_rad, &
        iradht_LW_rad, &
        iradht_SW_rad

    use stats_type, only: & 
        stat_assign ! Procedure

    implicit none

    ! Input Variable
    character(len= * ), dimension(nvarmax_rad_zt), intent(in) :: vars_rad_zt

    ! Input/Output Variable
    logical, intent(inout) :: l_error

    ! Local Varables
    integer :: i, k

    ! ---- Begin Code ----

    ! Default initialization for array indices for rad_zt

    iT_in_K_rad = 0
    ircil_rad = 0
    io3l_rad = 0
    irsnowm_rad = 0
    ircm_in_cloud_rad = 0
    icloud_frac_rad = 0
    iice_supersat_frac_rad = 0
    iradht_rad = 0
    iradht_LW_rad = 0
    iradht_SW_rad = 0

    ! Assign pointers for statistics variables rad_zt

    k = 1
    do i=1,rad_zt%nn

      select case ( trim(vars_rad_zt(i)) )

      case ('T_in_K_rad')
        iT_in_K_rad = k

        call stat_assign( iT_in_K_rad, "T_in_K_rad", & 
             "Temperature [K]", "K", rad_zt )
        k = k + 1

      case ('rcil_rad')
        ircil_rad = k

        call stat_assign( ircil_rad, "rcil_rad", & 
             "Ice mixing ratio [kg/kg]", "kg/kg", rad_zt )
        k = k + 1

      case ('o3l_rad')
        io3l_rad = k

        call stat_assign( io3l_rad, "o3l_rad", & 
             "Ozone mixing ratio [kg/kg]", "kg/kg", rad_zt )
        k = k + 1

      case ('rsnowm_rad')
        irsnowm_rad = k

        call stat_assign( irsnowm_rad, "rsnowm_rad", & 
             "Snow water mixing ratio [kg/kg]", "kg/kg", rad_zt )
        k = k + 1

      case ('rcm_in_cloud_rad')
        ircm_in_cloud_rad = k

        call stat_assign( ircm_in_cloud_rad, "rcm_in_cloud_rad", & 
             "rcm in cloud layer [kg/kg]", "kg/kg", rad_zt )
        k = k + 1

      case ('cloud_frac_rad')
        icloud_frac_rad = k

        call stat_assign( icloud_frac_rad, "cloud_frac_rad", & 
             "Cloud fraction (between 0 and 1) [-]", "count", rad_zt )
        k = k + 1
      
      case ('ice_supersat_frac_rad')
        iice_supersat_frac_rad = k

        call stat_assign( iice_supersat_frac_rad, "ice_supersat_frac_rad", & 
             "Ice cloud fraction (between 0 and 1) [-]", "count", rad_zt )
        k = k + 1

      case ('radht_rad')
        iradht_rad = k

        call stat_assign( iradht_rad, "radht_rad", & 
             "Total radiative heating rate [K/s]", "K/s", rad_zt )
        k = k + 1

      case ('radht_LW_rad')
        iradht_LW_rad = k

        call stat_assign( iradht_LW_rad, "radht_LW_rad", & 
             "Long-wave radiative heating rate [K/s]", "K/s", rad_zt )
        k = k + 1

      case ('radht_SW_rad')
        iradht_SW_rad = k

        call stat_assign( iradht_SW_rad, "radht_SW_rad", & 
             "Short-wave radiative heating rate [K/s]", "K/s", rad_zt )
        k = k + 1

      case default

        write(fstderr,*) 'Error:  unrecognized variable in vars_rad_zt:  ', trim( vars_rad_zt(i) )

        l_error = .true.  ! This will stop the run.


      end select

    end do

    return
  end subroutine stats_init_rad_zt

end module stats_rad_zt
