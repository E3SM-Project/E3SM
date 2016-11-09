!-----------------------------------------------------------------------
! $Id: stats_rad_zm_module.F90 7315 2014-09-30 20:49:54Z schemena@uwm.edu $
!===============================================================================

module stats_rad_zm_module

  implicit none

  private ! Default Scope

  public :: stats_init_rad_zm

! Constant parameters
  integer, parameter, public :: nvarmax_rad_zm = 250 ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_rad_zm( vars_rad_zm, l_error )

!     Description:
!     Initializes array indices for stats_rad_zm variables
!-----------------------------------------------------------------------

    use constants_clubb, only:  &
        fstderr ! Constant(s)

    use stats_variables, only: & 
        stats_rad_zm, &
        iFrad_LW_rad, & ! Variable(s)
        iFrad_SW_rad, &
        iFrad_SW_up_rad, &
        iFrad_LW_up_rad, &
        iFrad_SW_down_rad, &
        iFrad_LW_down_rad

    use stats_variables, only: &
      ifulwcl, ifdlwcl, ifdswcl, ifuswcl ! Variable(s)

    use stats_type_utilities, only: & 
        stat_assign ! Procedure


    implicit none

    ! Input Variable
    character(len= * ), dimension(nvarmax_rad_zm), intent(in) :: vars_rad_zm

    ! Input/Output Variable        
    logical, intent(inout) :: l_error

    ! Local Varables
    integer :: i, k

    ! ---- Begin Code ----

    ! Default initialization for array indices for stats_rad_zm

    iFrad_LW_rad = 0
    iFrad_SW_rad = 0
    iFrad_SW_up_rad = 0
    iFrad_LW_up_rad = 0
    iFrad_SW_down_rad = 0
    iFrad_LW_down_rad = 0

    ifulwcl = 0
    ifdlwcl = 0
    ifdswcl = 0
    ifuswcl = 0

!     Assign pointers for statistics variables stats_rad_zm

    k = 1
    do i=1,stats_rad_zm%num_output_fields

      select case ( trim(vars_rad_zm(i)) )

      case('fulwcl')
        ifulwcl = k
        call stat_assign( var_index=ifulwcl, var_name="fulwcl", &
             var_description="Upward clear-sky LW flux [W/m^2]", var_units="W/m^2", &
             l_silhs=.false., grid_kind=stats_rad_zm )
        k = k + 1

      case( 'fdlwcl' )
        ifdlwcl = k
        call stat_assign( var_index=ifdlwcl, var_name="fdlwcl", &
             var_description="Downward clear-sky LW flux [W/m^2]", var_units="W/m^2", &
             l_silhs=.false., grid_kind=stats_rad_zm )
        k = k + 1

      case( 'fdswcl' )
        ifdswcl = k
        call stat_assign( var_index=ifdswcl, var_name="fdswcl", &
             var_description="Downward clear-sky SW flux [W/m^2]", var_units="W/m^2", &
             l_silhs=.false., grid_kind=stats_rad_zm )
        k = k + 1

      case( 'fuswcl' )
        ifuswcl = k
        call stat_assign( var_index=ifuswcl, var_name="fuswcl", &
             var_description="Upward clear-sky SW flux [W/m^2]", var_units="W/m^2", &
             l_silhs=.false., grid_kind=stats_rad_zm )
        k = k + 1

      case ('Frad_LW_rad')
        iFrad_LW_rad = k

        call stat_assign( var_index=iFrad_LW_rad, var_name="Frad_LW_rad", &
             var_description="Net long-wave radiative flux [W/m^2]", var_units="W/m^2", &
             l_silhs=.false., grid_kind=stats_rad_zm )
        k = k + 1

      case ('Frad_SW_rad')
        iFrad_SW_rad = k

        call stat_assign( var_index=iFrad_SW_rad, var_name="Frad_SW_rad", &
             var_description="Net short-wave radiative flux [W/m^2]", var_units="W/m^2", &
             l_silhs=.false., grid_kind=stats_rad_zm )
        k = k + 1

      case ('Frad_SW_up_rad')
        iFrad_SW_up_rad = k

        call stat_assign( var_index=iFrad_SW_up_rad, var_name="Frad_SW_up_rad", &
             var_description="Short-wave upwelling radiative flux [W/m^2]", var_units="W/m^2", &
             l_silhs=.false., grid_kind=stats_rad_zm )
        k = k + 1

      case ('Frad_LW_up_rad')
        iFrad_LW_up_rad = k

        call stat_assign( var_index=iFrad_LW_up_rad, var_name="Frad_LW_up_rad", &
             var_description="Long-wave upwelling radiative flux [W/m^2]", var_units="W/m^2", &
             l_silhs=.false., grid_kind=stats_rad_zm )
        k = k + 1

      case ('Frad_SW_down_rad')
        iFrad_SW_down_rad = k

        call stat_assign( var_index=iFrad_SW_down_rad, var_name="Frad_SW_down_rad", &
             var_description="Short-wave downwelling radiative flux [W/m^2]", var_units="W/m^2", &
             l_silhs=.false., grid_kind=stats_rad_zm )
        k = k + 1

      case ('Frad_LW_down_rad')
        iFrad_LW_down_rad = k

        call stat_assign( var_index=iFrad_LW_down_rad, var_name="Frad_LW_down_rad", &
             var_description="Long-wave downwelling radiative flux [W/m^2]", var_units="W/m^2", &
             l_silhs=.false., grid_kind=stats_rad_zm )
        k = k + 1

      case default

        write(fstderr,*) 'Error:  unrecognized variable in vars_rad_zm:  ', trim( vars_rad_zm(i) )

        l_error = .true.  ! This will stop the run.


      end select

    end do

    return
  end subroutine stats_init_rad_zm

end module stats_rad_zm_module
