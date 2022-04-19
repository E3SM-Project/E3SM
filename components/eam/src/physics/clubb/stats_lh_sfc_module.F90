!-----------------------------------------------------------------------
! $Id$
!===============================================================================

module stats_lh_sfc_module


  implicit none

  private ! Set Default Scope

  public :: stats_init_lh_sfc

  ! Constant parameters
  integer, parameter, public :: nvarmax_lh_sfc = 10  ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_lh_sfc( vars_lh_sfc, l_error, & !intent(in)
                                stats_lh_sfc ) ! intent(inout)

! Description:
!   Initializes array indices for stats_lh_sfc
! References:
!   None
!-----------------------------------------------------------------------

    use constants_clubb, only: &
        fstderr ! Constant(s)

    use stats_variables, only: & 
        ilh_morr_snow_rate, & ! Variable(s)
        ilh_vwp, &
        ilh_lwp, &
        ilh_sample_weights_sum, &
        ilh_sample_weights_avg, &
        ik_lh_start
      
    use stats_type_utilities, only: & 
        stat_assign ! Procedure

    use stats_type, only: stats ! Type

    implicit none

    type (stats), target, intent(inout) :: &
      stats_lh_sfc

    ! External
    intrinsic :: trim

    ! Input Variable
    character(len= * ), dimension(nvarmax_lh_sfc), intent(in) :: vars_lh_sfc

    ! Input / Output Variable        
    logical, intent(inout) :: l_error

    ! Local Varables
    integer :: i, k

    ! ---- Begin Code ----

    ! Default initialization for array indices for stats_sfc is zero (see module
    ! stats_variables)

    ! Assign pointers for statistics variables stats_sfc

    k = 1
    do i = 1, stats_lh_sfc%num_output_fields

      select case ( trim( vars_lh_sfc(i) ) )

      case ( 'lh_morr_snow_rate' )
        ilh_morr_snow_rate = k
        call stat_assign( var_index=ilh_morr_snow_rate, var_name="lh_morr_snow_rate", & !intent(in)
             var_description="Snow+Ice+Graupel fallout rate from Morrison scheme [mm/day]", & ! In
             var_units="mm/day", l_silhs=.true., & ! intent(in)
             grid_kind=stats_lh_sfc ) ! intent(inout)
        k = k + 1

      case ( 'lh_vwp' )
        ilh_vwp = k
        call stat_assign( var_index=ilh_vwp, var_name="lh_vwp", & ! intent(in)
             var_description="Vapor water path [kg/m^2]", var_units="kg/m^2", & ! intent(in)
             l_silhs=.true., & ! intent(in)
             grid_kind=stats_lh_sfc ) ! intent(inout)
        k = k + 1

      case ( 'lh_lwp' )
        ilh_lwp = k
        call stat_assign( var_index=ilh_lwp, var_name="lh_lwp", & ! intent(in)
             var_description="Liquid water path [kg/m^2]", var_units="kg/m^2", & ! intent(in)
             l_silhs=.true., & ! intent(in)
             grid_kind=stats_lh_sfc ) ! intent(inout)
        k = k + 1

      case ( 'k_lh_start' )
        ik_lh_start = k
        call stat_assign( var_index=ik_lh_start, var_name="k_lh_start", & ! intent(in)
             var_description="Index of height level for SILHS sampling preferentially within &
                             &cloud [integer]", var_units="integer", l_silhs=.true., & 
             grid_kind=stats_lh_sfc ) ! intent(inout)
        k = k + 1

      case ( 'lh_sample_weights_sum' )
        ilh_sample_weights_sum = k
        call stat_assign( var_index=ilh_sample_weights_sum, & ! intent(in)
             var_name="lh_sample_weights_sum", & ! intent(in)
             var_description="Sum of the sample point weights [-]", var_units="-", & ! intent(in)
             l_silhs=.true., & ! intent(in)
             grid_kind=stats_lh_sfc ) ! intent(inout)
        k = k + 1
        
      case ( 'lh_sample_weights_avg' )
        ilh_sample_weights_avg = k
        call stat_assign( var_index=ilh_sample_weights_avg, & ! intent(in)
             var_name="lh_sample_weights_avg", & ! intent(in)
             var_description="Average of the sample point weights [-]", &  ! intent(in)
             var_units="-", & ! intent(in)
             l_silhs=.true., &  ! intent(in)
             grid_kind=stats_lh_sfc ) ! intent(inout)
        k = k + 1

      case default
        write(fstderr,*) 'Error:  unrecognized variable in vars_lh_sfc:  ',  &
              trim( vars_lh_sfc(i) )
        l_error = .true.  ! This will stop the run.

      end select

    end do ! i = 1, stats_lh_sfc%num_output_fields

    return
  end subroutine stats_init_lh_sfc

end module stats_lh_sfc_module

