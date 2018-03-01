!-------------------------------------------------------------------------
! $Id: advance_helper_module.F90 7381 2014-11-11 23:59:39Z schemena@uwm.edu $
!===============================================================================
module advance_helper_module

! Description:
!   This module contains helper methods for the advance_* modules.
!------------------------------------------------------------------------

  implicit none

  public :: &
    set_boundary_conditions_lhs, &
    set_boundary_conditions_rhs, &
    calc_stability_correction

  private ! Set Default Scope

  contains

  !---------------------------------------------------------------------------
  subroutine set_boundary_conditions_lhs( diag_index, low_bound, high_bound, lhs, &
                                      diag_index2, low_bound2, high_bound2 )

  ! Description:
  !   Sets the boundary conditions for a left-hand side LAPACK matrix.
  !
  ! References:
  !   none
  !---------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Exernal 
    intrinsic :: present

    ! Input Variables
    integer, intent(in) :: &
      diag_index, low_bound, high_bound ! boundary indexes for the first variable

    ! Input / Output Variables
    real( kind = core_rknd ), dimension(:,:), intent(inout) :: &
      lhs ! left hand side of the LAPACK matrix equation

    ! Optional Input Variables
    integer, intent(in), optional :: &
      diag_index2, low_bound2, high_bound2 ! boundary indexes for the second variable

    ! --------------------- BEGIN CODE ----------------------

    if ( ( present( low_bound2 ) .or. present( high_bound2 ) ) .and. &
         ( .not. present( diag_index2 ) ) ) then

      stop "Boundary index provided without diag_index."

    end if

    ! Set the lower boundaries for the first variable
    lhs(:,low_bound) = 0.0_core_rknd
    lhs(diag_index,low_bound) = 1.0_core_rknd

    ! Set the upper boundaries for the first variable
    lhs(:,high_bound) = 0.0_core_rknd
    lhs(diag_index,high_bound) = 1.0_core_rknd

    ! Set the lower boundaries for the second variable, if it is provided
    if ( present( low_bound2 ) ) then

      lhs(:,low_bound2) = 0.0_core_rknd
      lhs(diag_index2,low_bound2) = 1.0_core_rknd

    end if

    ! Set the upper boundaries for the second variable, if it is provided
    if ( present( high_bound2 ) ) then

      lhs(:,high_bound2) = 0.0_core_rknd
      lhs(diag_index2,high_bound2) = 1.0_core_rknd

    end if

    return
  end subroutine set_boundary_conditions_lhs

  !--------------------------------------------------------------------------
  subroutine set_boundary_conditions_rhs( &
               low_value, low_bound, high_value, high_bound, &
               rhs, &
               low_value2, low_bound2, high_value2, high_bound2 )

  ! Description:
  !   Sets the boundary conditions for a right-hand side LAPACK vector.
  !
  ! References:
  !   none
  !---------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Exernal 
    intrinsic :: present

    ! Input Variables

    ! The values for the first variable
    real( kind = core_rknd ), intent(in) :: low_value, high_value

    ! The bounds for the first variable
    integer, intent(in) :: low_bound, high_bound

    ! Input / Output Variables

    ! The right-hand side vector
    real( kind = core_rknd ), dimension(:), intent(inout) :: rhs

    ! Optional Input Variables

    ! The values for the second variable
    real( kind = core_rknd ), intent(in), optional :: low_value2, high_value2

    ! The bounds for the second variable
    integer, intent(in), optional :: low_bound2, high_bound2


    ! -------------------- BEGIN CODE ------------------------

    ! Stop execution if a boundary was provided without a value
    if ( (present( low_bound2 ) .and. (.not. present( low_value2 ))) .or. &
         (present( high_bound2 ) .and. (.not. present( high_value2 ))) ) then

      stop "Boundary condition provided without value."

    end if

    ! Set the lower and upper bounds for the first variable
    rhs(low_bound) = low_value
    rhs(high_bound) = high_value

    ! If a lower bound was given for the second variable, set it
    if ( present( low_bound2 ) ) then
      rhs(low_bound2) = low_value2
    end if

    ! If an upper bound was given for the second variable, set it
    if ( present( high_bound2 ) ) then
      rhs(high_bound2) = high_value2
    end if

    return
  end subroutine set_boundary_conditions_rhs

  !===============================================================================
  function calc_stability_correction( thlm, Lscale, em ) &
    result ( stability_correction )
      !
      ! Description:
      !   Stability Factor
      !
      ! References:
      !
      !--------------------------------------------------------------------

      use parameters_model, only: &
        T0 ! Variables(s)

      use constants_clubb, only: &
        zero, & ! Constant(s)
        grav

      use grid_class, only:  &
        gr, & ! Variable(s)
        zt2zm, & ! Procedure(s)
        ddzt

      use clubb_precision, only:  &
        core_rknd ! Variable(s)

      implicit none

      ! Input Variables
      real( kind = core_rknd ), intent(in), dimension(gr%nz) :: &
        Lscale,          & ! Turbulent mixing length                   [m]
        em,              & ! Turbulent Kinetic Energy (TKE)            [m^2/s^2]
        thlm               ! th_l (thermo. levels)                     [K]

      ! Result
      real( kind = core_rknd ), dimension(gr%nz) :: &
        stability_correction

      ! Local Variables
      real( kind = core_rknd ) :: &
        lambda0_stability_coef ! []

      real( kind = core_rknd ), dimension(gr%nz) :: &
        brunt_vaisala_freq, & !  []
        lambda0_stability

      !------------ Begin Code --------------
      ! lambda0_stability_coef = 0.025_core_rknd
      ! changed to 0.030 to provide a simulation similar to track02 simulation
      lambda0_stability_coef = 0.030_core_rknd
      brunt_vaisala_freq = ( grav / T0 ) * ddzt( thlm )
      lambda0_stability = merge( lambda0_stability_coef, zero, brunt_vaisala_freq > zero )

      stability_correction = 1.0_core_rknd &
        + min( lambda0_stability * brunt_vaisala_freq * zt2zm( Lscale )**2 / em, 3.0_core_rknd )

      return
    end function calc_stability_correction

end module advance_helper_module
