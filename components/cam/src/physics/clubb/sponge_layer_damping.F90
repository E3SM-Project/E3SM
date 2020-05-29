!---------------------------------------------------------------------------
! $Id$
!===============================================================================
module sponge_layer_damping

  ! Description:
  ! This module is used for damping variables in upper altitudes of the grid.
  !
  ! References:
  !   None
  !-------------------------------------------------------------------------

  use clubb_precision, only: &
      core_rknd ! Variable(s)

  implicit none

  public :: sponge_damp_xm,             & ! Procedure(s)
            sponge_damp_xp2,            &
            sponge_damp_xp3,            &
            initialize_tau_sponge_damp, &
            finalize_tau_sponge_damp,   &
            sponge_damp_settings,       & ! Variable type(s)
            sponge_damp_profile


  type sponge_damp_settings

    real( kind = core_rknd ) :: &
      tau_sponge_damp_min, & ! Minimum damping time scale (model top)        [s]
      tau_sponge_damp_max, & ! Maximum damping time scale (damp layer base)  [s]
      sponge_damp_depth      ! Damping depth as a fraction of domain height  [-]

    logical :: &
      l_sponge_damping       ! True if damping is being used

  end type sponge_damp_settings

  type sponge_damp_profile

    real( kind = core_rknd ), allocatable, dimension(:) :: &
      tau_sponge_damp    ! Damping time scale    [1/s]

    real( kind = core_rknd ) :: &
      sponge_layer_depth    ! Depth of sponge damping layer  [m]

  end type sponge_damp_profile


  type(sponge_damp_settings), public :: &
    thlm_sponge_damp_settings,    & ! Variable(s)
    rtm_sponge_damp_settings,     &
    uv_sponge_damp_settings,      &
    wp2_sponge_damp_settings,     &
    wp3_sponge_damp_settings,     &
    up2_vp2_sponge_damp_settings
!$omp threadprivate( thlm_sponge_damp_settings, rtm_sponge_damp_settings, &
!$omp                uv_sponge_damp_settings, wp2_sponge_damp_settings, &
!$omp                wp3_sponge_damp_settings, up2_vp2_sponge_damp_settings )

  type(sponge_damp_profile), public :: &
    thlm_sponge_damp_profile,    & ! Variable(s)
    rtm_sponge_damp_profile,     &
    uv_sponge_damp_profile,      &
    wp2_sponge_damp_profile,     &
    wp3_sponge_damp_profile,     &
    up2_vp2_sponge_damp_profile
!$omp threadprivate( thlm_sponge_damp_profile, rtm_sponge_damp_profile, &
!$omp                uv_sponge_damp_profile, wp2_sponge_damp_profile, &
!$omp                wp3_sponge_damp_profile, up2_vp2_sponge_damp_profile )


  private

  contains

  !=============================================================================
  function sponge_damp_xm( dt, z, xm_ref, xm, damping_profile ) result( xm_p )

    ! Description:
    ! Damps specified mean field toward a reference profile.  The module must be
    ! initialized for this function to work.  Otherwise a stop is issued.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    !  "Sponge"-layer damping at the domain top region

    use grid_class, only: &
        gr    ! Variable(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! External
    intrinsic :: allocated

    ! Input Variable(s)
    real( kind = core_rknd ), intent(in) :: &
      dt    ! Model Timestep  [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      z,      & ! Height of model grid levels                [m]
      xm_ref, & ! Reference profile of x to damp xm towards  [units vary]
      xm        ! Mean field being damped                    [units vary]

    type(sponge_damp_profile), intent(in) :: &
      damping_profile

    ! Output Variable
    real( kind = core_rknd ), dimension(gr%nz) :: &
      xm_p   ! Damped value of xm  [units_vary]

    ! Local Variable(s)
    real( kind = core_rknd ) :: &
      dt_on_tau    ! Ratio of timestep to damping timescale  [-]

    integer :: k

    ! ---- Begin Code ----

    if ( allocated( damping_profile%tau_sponge_damp ) ) then

       xm_p = xm
     
       do k = gr%nz, 1, -1

          ! The height of the model top is gr%zm(gr%nz).
          if ( gr%zm(gr%nz) - z(k) < damping_profile%sponge_layer_depth ) then

             ! Vince Larson used implicit discretization in order to 
             ! reduce noise in rtm in cloud_feedback_s12 (CGILS) 
             !xm_p(k) = xm(k) - real( ( ( xm(k) - xm_ref(k) ) / & 
             !                  damping_profile%tau_sponge_damp(k) ) * dt )
             dt_on_tau = dt / damping_profile%tau_sponge_damp(k)

             ! Really, we should be using xm_ref at time n+1 rather than n.
             ! However, for steady profiles of xm_ref, it won't matter.        
             xm_p(k) = ( xm(k) + dt_on_tau * xm_ref(k) ) / &
                             ( 1.0_core_rknd + dt_on_tau )
             ! End Vince Larson's change

          else ! gr%zm(gr%nz) - z(k) >= damping_profile%sponge_layer_depth

             ! Below sponge damping layer; exit loop.
             exit

          endif ! gr%zm(gr%nz) - z(k) < damping_profile%sponge_layer_depth


       enddo ! k = gr%nz, 1, -1

    else

       stop "tau_sponge_damp in sponge_damp_xm used before initialization"

    endif


    return

  end function sponge_damp_xm

  !=============================================================================
  function sponge_damp_xp2( dt, z, xp2, x_tol_sqd, damping_profile ) &
  result( xp2_damped )

    ! Description:
    ! Calculates the effects of "sponge"-layer damping on the variance of x,
    ! xp2.
    !
    ! Sponge damping on a local value of x is given by the equation:
    !
    ! x_d = x - ( delta_t / tau ) * ( x - <x> ),
    !
    ! where x is the local value prior to damping, x_d is the damped local value
    ! of x, <x> is the grid-level mean value of x, delta_t is the model time
    ! step duration, and tau is the damping time scale.  Since delta_t / tau has
    ! the same value everywhere at a grid level, the grid-level mean of x does
    ! not change as a result of damping.
    !
    ! Subtracting <x> from both sides:
    !
    ! x_d - <x> = ( x - <x> ) - ( delta_t / tau ) * ( x - <x> ),
    !
    ! which results in:
    !
    ! x_d - <x> = ( 1 - delta_t / tau ) * ( x - <x> ).
    !
    ! Squaring both sides:
    !
    ! ( x_d - <x> )^2 = ( 1 - delta_t / tau )^2 * ( x - <x> )^2.
    !
    ! After averaging both sides, the damped value of xp2 is:
    !
    ! < x_d'^2 > = ( 1 - delta_t / tau )^2 * < x'^2 >.
    !
    ! Any sponge damping is applied to (predictive) xp2 after the field has been
    ! advanced in time.  This allows sponge damping to be applied in an implicit
    ! manner.  The damped value of xp2 must also be limited at a minimum value
    ! of x_tol^2.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable(s)

    use constants_clubb, only: &
        one    ! Constant(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variable(s)
    real( kind = core_rknd ), intent(in) :: &
      dt    ! Model Timestep  [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      z,   & ! Height of model grid levels               [m]
      xp2    ! Variance of x, <x'^2>, prior to damping   [units vary]

    real( kind = core_rknd ), intent(in) :: &
      x_tol_sqd    ! Square of the tolerance value of x    [units vary]

    type(sponge_damp_profile), intent(in) :: &
      damping_profile

    ! Output Variable
    real( kind = core_rknd ), dimension(gr%nz) :: &
      xp2_damped    ! Variance of x, <x'^2>, after damping   [units vary]

    ! Local Variable(s)
    real( kind = core_rknd ) :: &
      dt_on_tau    ! Ratio of model time step duration to damping timescale  [-]

    integer :: &
      k    ! Loop index


    if ( allocated( damping_profile%tau_sponge_damp ) ) then

       ! Set the entire profile of <x'^2> after damping to the profile of <x'^2>
       ! before damping.  The values of <x'^2> after damping will be overwritten
       ! at any levels where "sponge"-layer damping occurs.
       xp2_damped = xp2
     
       do k = gr%nz, 1, -1

          ! The height of the model top is gr%zm(gr%nz).
          if ( gr%zm(gr%nz) - z(k) < damping_profile%sponge_layer_depth ) then

             ! Calculate the value of delta_t / tau at the grid level.
             dt_on_tau = dt / damping_profile%tau_sponge_damp(k)

             ! Calculate the damped value of <x'^2>.
             xp2_damped(k) = ( one - dt_on_tau )**2 * xp2(k)

             ! The damped value of <x'^2> needs to be greater than or equal to
             ! x_tol^2.
             xp2_damped(k) = max( xp2_damped(k), x_tol_sqd )

          else ! gr%zm(gr%nz) - z(k) >= damping_profile%sponge_layer_depth

             ! Below sponge damping layer; exit loop.
             exit

          endif ! gr%zm(gr%nz) - z(k) < damping_profile%sponge_layer_depth


       enddo ! k = gr%nz, 1, -1

    else

       stop "tau_sponge_damp in sponge_damp_xp2 used before initialization"

    endif


    return

  end function sponge_damp_xp2

  !=============================================================================
  function sponge_damp_xp3( dt, z, xp3, damping_profile ) &
  result( xp3_damped )

    ! Description:
    ! Calculates the effects of "sponge"-layer damping on xp3.
    !
    ! Sponge damping on a local value of x is given by the equation:
    !
    ! x_d = x - ( delta_t / tau ) * ( x - <x> ),
    !
    ! where x is the local value prior to damping, x_d is the damped local value
    ! of x, <x> is the grid-level mean value of x, delta_t is the model time
    ! step duration, and tau is the damping time scale.  Since delta_t / tau has
    ! the same value everywhere at a grid level, the grid-level mean of x does
    ! not change as a result of damping.
    !
    ! Subtracting <x> from both sides:
    !
    ! x_d - <x> = ( x - <x> ) - ( delta_t / tau ) * ( x - <x> ),
    !
    ! which results in:
    !
    ! x_d - <x> = ( 1 - delta_t / tau ) * ( x - <x> ).
    !
    ! Taking both sides to the third power:
    !
    ! ( x_d - <x> )^3 = ( 1 - delta_t / tau )^3 * ( x - <x> )^3.
    !
    ! After averaging both sides, the damped value of xp3 is:
    !
    ! < x_d'^3 > = ( 1 - delta_t / tau )^3 * < x'^3 >.
    !
    ! Any sponge damping is applied to (predictive) xp3 after the field has been
    ! advanced in time.  This allows sponge damping to be applied in an implicit
    ! manner.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable(s)

    use constants_clubb, only: &
        one    ! Constant(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variable(s)
    real( kind = core_rknd ), intent(in) :: &
      dt    ! Model Timestep  [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      z,   & ! Height of model grid levels     [m]
      xp3    ! <x'^3> prior to damping         [units vary]

    type(sponge_damp_profile), intent(in) :: &
      damping_profile

    ! Output Variable
    real( kind = core_rknd ), dimension(gr%nz) :: &
      xp3_damped    ! <x'^3> after damping   [units vary]

    ! Local Variable(s)
    real( kind = core_rknd ) :: &
      dt_on_tau    ! Ratio of model time step duration to damping timescale  [-]

    integer :: &
      k    ! Loop index


    if ( allocated( damping_profile%tau_sponge_damp ) ) then

       ! Set the entire profile of <x'^3> after damping to the profile of <x'^3>
       ! before damping.  The values of <x'^3> after damping will be overwritten
       ! at any levels where "sponge"-layer damping occurs.
       xp3_damped = xp3
     
       do k = gr%nz, 1, -1

          ! The height of the model top is gr%zm(gr%nz).
          if ( gr%zm(gr%nz) - z(k) < damping_profile%sponge_layer_depth ) then

             ! Calculate the value of delta_t / tau at the grid level.
             dt_on_tau = dt / damping_profile%tau_sponge_damp(k)

             ! Calculate the damped value of <x'^3>.
             xp3_damped(k) = ( one - dt_on_tau )**3 * xp3(k)

          else ! gr%zm(gr%nz) - z(k) >= damping_profile%sponge_layer_depth

             ! Below sponge damping layer; exit loop.
             exit

          endif ! gr%zm(gr%nz) - z(k) < damping_profile%sponge_layer_depth


       enddo ! k = gr%nz, 1, -1

    else

       stop "tau_sponge_damp in sponge_damp_xp3 used before initialization"

    endif


    return

  end function sponge_damp_xp3

  !=============================================================================
  subroutine initialize_tau_sponge_damp( dt, z, settings, damping_profile )

    ! Description:
    ! Initialize time scale, tau_sponge_damp, used for damping.  The time scale
    ! attains its maximum value, tau_sponge_damp_max, at the bottom of the
    ! "sponge" damping layer, which results in minimal damping.  Likewise, the
    ! time scale attains its minimum value, tau_sponge_damp_min, at the top of
    ! the model, which results in maximum damping.  At levels in-between the top
    ! of the model and the base of the sponge damping layer, the value of
    ! tau_sponge_damp is in-between tau_sponge_damp_min and tau_sponge_damp_max,
    ! as calculated by an interpolation formula.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd    ! Variable(s)
    
    use constants_clubb, only: &
        two,     & ! Constant(s)
        fstderr

    use grid_class, only: &
        gr    ! Variable(s)

!    use interpolation, only: &
!        lin_interpolate_two_points    ! Procedure(s)

    implicit none

    ! Input Variable(s)
    real( kind = core_rknd ), intent(in) :: &
      dt    ! Model Timestep    [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      z    ! Height of model grid levels    [m]

    type(sponge_damp_settings), intent(in) :: &
      settings

    ! Output Variable(s)
    type(sponge_damp_profile), intent(out) :: &
      damping_profile

    ! Local Variable(s)
    real( kind = core_rknd ) :: &
      tau_sponge_damp_exponent  ! Exponent in calculation of tau_sponge_damp [-]

    integer :: &
      k    ! Loop index

    ! ---- Begin Code ----

    ! Allocate the damping time scale.
    allocate( damping_profile%tau_sponge_damp(1:gr%nz) )

    ! Calculate the depth of the sponge layer.
    ! The height of the model top is gr%zm(gr%nz).
    damping_profile%sponge_layer_depth &
    = settings%sponge_damp_depth * gr%zm(gr%nz)

    ! Check the value of tau_sponge_damp_min.
    if ( settings%tau_sponge_damp_min < two * dt ) then
       write(fstderr,*) "Error:  tau_sponge_damp_min is too small!"
       write(fstderr,*) "It must be at least 2.0 * dt"
       stop
    endif

    ! Calculate the value of the damping time scale, tau_sponge_damp, at levels
    ! that are within the sponge damping layer.
    do k = gr%nz, 1, -1

       ! The height of the model top is gr%zm(gr%nz).
       if ( gr%zm(gr%nz) - z(k) < damping_profile%sponge_layer_depth ) then

          ! Vince Larson added code to use standard linear interpolation.
          ! Brian Griffin reverted the linear interpolation in order to use code
          ! that is similar to what is found in SAM LES.

          tau_sponge_damp_exponent &
          = ( gr%zm(gr%nz) - z(k) ) / damping_profile%sponge_layer_depth

          damping_profile%tau_sponge_damp(k) &
          = settings%tau_sponge_damp_min &
            * ( settings%tau_sponge_damp_max &
                / settings%tau_sponge_damp_min )**tau_sponge_damp_exponent

          !damping_profile%tau_sponge_damp(k) &
          != lin_interpolate_two_points( z(k), gr%zm(gr%nz), &
          !                              gr%zm(gr%nz) &
          !                              - damping_profile%sponge_layer_depth, &
          !                              settings%tau_sponge_damp_min, &
          !                              settings%tau_sponge_damp_max )

          ! End Vince Larson's change
          ! End Brian Griffin's rebellious reversion.

       else ! gr%zm(gr%nz) - z(k) >= damping_profile%sponge_layer_depth

          ! Below sponge damping layer; exit loop.
          exit

       endif ! gr%zm(gr%nz) - z(k) < damping_profile%sponge_layer_depth

    enddo ! k = gr%nz, 1, -1


    return

  end subroutine initialize_tau_sponge_damp

  !=============================================================================
  subroutine finalize_tau_sponge_damp( damping_profile )

    ! Description:
    ! Frees memory allocated in initialize_tau_sponge_damp
    ! 
    ! References:
    ! None
    !-----------------------------------------------------------------------

    implicit none

    ! Input/Output Variable(s)
    type(sponge_damp_profile), intent(inout) :: &
      damping_profile ! Information for damping the profile

    ! ---- Begin Code ----

    ! Deallocate the damping time scale.
    deallocate( damping_profile%tau_sponge_damp )


    return

  end subroutine finalize_tau_sponge_damp

!===============================================================================

end module sponge_layer_damping
