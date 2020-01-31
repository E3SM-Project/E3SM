!-----------------------------------------------------------------------
! $Id: stats_LH_sfc.F90 6100 2013-03-08 17:53:44Z dschanen@uwm.edu $

module stats_LH_sfc


  implicit none

  private ! Set Default Scope

  public :: stats_init_LH_sfc

  ! Constant parameters
  integer, parameter, public :: nvarmax_LH_sfc = 10  ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_LH_sfc( vars_LH_sfc, l_error )

! Description:
!   Initializes array indices for LH_sfc
! References:
!   None
!-----------------------------------------------------------------------

    use constants_clubb, only: &
      fstderr ! Constant(s)

    use stats_variables, only: & 
      LH_sfc ! Variable(s)

    use stats_variables, only: & 
      iLH_morr_rain_rate, & ! Variable(s)
      iLH_morr_snow_rate, &
      iLH_vwp, &
      iLH_lwp
      
    use stats_type, only: & 
        stat_assign ! Procedure

    implicit none

    ! Input Variable
    character(len= * ), dimension(nvarmax_LH_sfc), intent(in) :: vars_LH_sfc

    ! Output Variable        
    logical, intent(inout) :: l_error

    ! Local Varables
    integer :: i, k

    ! ---- Begin Code ----

    ! Default initialization for array indices for sfc

    iLH_morr_rain_rate = 0
    iLH_morr_snow_rate = 0
    iLH_vwp = 0
    iLH_lwp = 0

    ! Assign pointers for statistics variables sfc

    k = 1
    do i=1,LH_sfc%nn

      select case ( trim( vars_LH_sfc(i) ) )

      case ( 'LH_morr_rain_rate' )
        iLH_morr_rain_rate = k
        call stat_assign( iLH_morr_rain_rate, "LH_morr_rain_rate", & 
             "Total precip fallout rate from Morrison scheme [mm/day]","mm/day", LH_sfc )
        k = k + 1

      case ( 'LH_morr_snow_rate' )
        iLH_morr_snow_rate = k
        call stat_assign( iLH_morr_snow_rate, "LH_morr_snow_rate", & 
             "Snow+Ice+Graupel fallout rate from Morrison scheme [mm/day]","mm/day", LH_sfc )
        k = k + 1

      case ( 'LH_vwp' )
        iLH_vwp = k
        call stat_assign( iLH_vwp, "LH_vwp", & 
             "Vapor water path [kg/m^2]","kg/m^2", LH_sfc )
        k = k + 1

      case ( 'LH_lwp' )
        iLH_lwp = k
        call stat_assign( iLH_lwp, "LH_lwp", & 
             "Liquid water path [kg/m^2]","kg/m^2", LH_sfc )
        k = k + 1

      case default
        write(fstderr,*) 'Error:  unrecognized variable in vars_LH_sfc:  ',  &
              trim( vars_LH_sfc(i) )
        l_error = .true.  ! This will stop the run.

      end select

    end do

    return
  end subroutine stats_init_LH_sfc

end module stats_LH_sfc

