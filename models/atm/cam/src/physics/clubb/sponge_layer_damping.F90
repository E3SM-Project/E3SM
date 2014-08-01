!$Id: sponge_layer_damping.F90 5623 2012-01-17 17:55:26Z connork@uwm.edu $
module sponge_layer_damping
! Description:
!   This module is used for damping variables in upper altitudes of the grid.
!
! References:
!   None
!---------------------------------------------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  public :: sponge_damp_xm, initialize_tau_sponge_damp, finalize_tau_sponge_damp, &
            sponge_damp_settings, sponge_damp_profile


  type sponge_damp_settings

    real( kind = core_rknd ) :: &
      tau_sponge_damp_min, & ! Minimum damping time-scale (at the top) [s]
      tau_sponge_damp_max, & ! Maximum damping time-scale (base of damping layer) [s]
      sponge_damp_depth      ! damping depth as a fraction of domain height [-]

    logical :: &
      l_sponge_damping       ! True if damping is being used

  end type sponge_damp_settings

  type sponge_damp_profile
    real( kind = core_rknd ), pointer, dimension(:) :: &
      tau_sponge_damp ! Damping factor

    integer :: &
      n_sponge_damp ! Number of levels damped

  end type sponge_damp_profile


  type(sponge_damp_settings), public :: &
    thlm_sponge_damp_settings, &
    rtm_sponge_damp_settings, &
    uv_sponge_damp_settings
!$omp threadprivate(thlm_sponge_damp_settings, rtm_sponge_damp_settings, uv_sponge_damp_settings)

  type(sponge_damp_profile), public :: &
    thlm_sponge_damp_profile, &
    rtm_sponge_damp_profile, &
    uv_sponge_damp_profile
!$omp threadprivate(thlm_sponge_damp_profile, rtm_sponge_damp_profile, uv_sponge_damp_profile)


  private

  contains

  !---------------------------------------------------------------------------------------------
  function sponge_damp_xm( dt, xm_ref, xm, damping_profile ) result( xm_p )
    !
    ! Description:
    !   Damps specified variable. The module must be initialized for
    !   this function to work. Otherwise a stop is issued.
    !
    ! References:
    !   None
    !-------------------------------------------------------------------------------------------

    !  "Sponge"-layer damping at the domain top region

    use grid_class, only: gr ! Variable(s)

    use clubb_precision, only: time_precision, core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: associated

    ! Input Variable(s)
    real(kind=time_precision), intent(in) :: dt ! Model Timestep

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      xm_ref ! Reference to damp to [-]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      xm ! Variable being damped [-]

    type(sponge_damp_profile), intent(in) :: &
        damping_profile

    ! Output Variable(s)
    real( kind = core_rknd ), dimension(gr%nz) :: xm_p ! Variable damped [-]

    real( kind = core_rknd ) :: dt_on_tau ! Ratio of timestep to damping timescale [-]

    integer :: k

    ! ---- Begin Code ----

    if ( associated( damping_profile%tau_sponge_damp ) ) then

      xm_p = xm

      do k = gr%nz, gr%nz-damping_profile%n_sponge_damp, -1

! Vince Larson used implicit discretization in order to 
! reduce noise in rtm in cloud_feedback_s12 (CGILS) 
!        xm_p(k) = xm(k) - real( ( ( xm(k) - xm_ref(k) ) / & 
!                        damping_profile%tau_sponge_damp(k) ) * dt )
        dt_on_tau = real( dt, kind = core_rknd ) / damping_profile%tau_sponge_damp(k)

! Really, we should be using xm_ref at time n+1 rather than n.
! However, for steady profiles of xm_ref, it won't matter.        
        xm_p(k) = ( xm(k) + dt_on_tau * xm_ref(k) ) / &
                        ( 1.0_core_rknd + dt_on_tau )
! End Vince Larson's change
      end do ! k

    else

      stop "tau_sponge_damp in damping used before initialization"

    end if

    return
  end function sponge_damp_xm

  !---------------------------------------------------------------------------------------------
  subroutine initialize_tau_sponge_damp( dt, settings, damping_profile )
    !
    ! Description:
    !   Initialize tau_sponge_damp used for damping
    !
    ! References:
    !   None
    !-------------------------------------------------------------------------------------------
    use clubb_precision, only: time_precision, core_rknd ! Variable(s)
    
    use constants_clubb, only: fstderr ! Constant(s)

    use grid_class, only: gr ! Variable(s)

    use interpolation, only: lin_int ! function

    implicit none

    ! Input Variable(s)
    real(kind=time_precision), intent(in) :: dt ! Model Timestep [s]

    type(sponge_damp_settings), intent(in) :: &
        settings

    type(sponge_damp_profile), intent(out) :: &
      damping_profile

    integer :: k ! Loop iterator

    ! ---- Begin Code ----

    allocate( damping_profile%tau_sponge_damp(1:gr%nz))

    if( settings%tau_sponge_damp_min < 2._core_rknd * real( dt, kind = core_rknd ) ) then
      write(fstderr,*) 'Error: in damping() tau_sponge_damp_min is too small!'
      stop
    end if

    do k=gr%nz,1,-1
      if(gr%zt(gr%nz)-gr%zt(k) < settings%sponge_damp_depth*gr%zt(gr%nz)) then
        damping_profile%n_sponge_damp=gr%nz-k+1
      endif
    end do

    do k=gr%nz,gr%nz-damping_profile%n_sponge_damp,-1
! Vince Larson added code to use standard linear interpolation.
!      damping_profile%tau_sponge_damp(k) = settings%tau_sponge_damp_min *&
!        (settings%tau_sponge_damp_max/settings%tau_sponge_damp_min)** &
!        ( ( gr%zt(gr%nz)-gr%zt(k) ) / &
!          (gr%zt(gr%nz) - gr%zt( gr%nz-damping_profile%n_sponge_damp ) ) )
      damping_profile%tau_sponge_damp(k) =                                     &
        lin_int( gr%zt(k), gr%zt(gr%nz),                                     &
          gr%zt(gr%nz) - gr%zt( gr%nz-damping_profile%n_sponge_damp ) ,    &
          settings%tau_sponge_damp_min, settings%tau_sponge_damp_max )         
! End Vince Larson's change
    end do

    return
  end subroutine initialize_tau_sponge_damp

  !---------------------------------------------------------------------------------------------
  subroutine finalize_tau_sponge_damp( damping_profile )
    !
    ! Description:
    !   Frees memory allocated in initialize_tau_sponge_damp
    ! 
    ! References:
    !   None
    !-------------------------------------------------------------------------------------------
    implicit none

    ! Input/Output Variable(s)
    type(sponge_damp_profile), intent(inout) :: &
      damping_profile ! Information for damping the profile

    ! ---- Begin Code ----

    deallocate( damping_profile%tau_sponge_damp )

    return
  end subroutine finalize_tau_sponge_damp


end module sponge_layer_damping
