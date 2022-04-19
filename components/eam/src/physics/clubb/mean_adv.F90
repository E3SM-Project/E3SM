!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module mean_adv

  ! Description:
  ! Module mean_adv computes the mean advection terms for all of the
  ! time-tendency (prognostic) equations in the CLUBB parameterization.  All of
  ! the mean advection terms are solved for completely implicitly, and therefore
  ! become part of the left-hand side of their respective equations.
  !
  ! Function term_ma_zt_lhs handles the mean advection terms for the variables
  ! located at thermodynamic grid levels.  These variables are:  rtm, thlm, wp3,
  ! all hydrometeor species, and sclrm.
  !
  ! Function term_ma_zm_lhs handles the mean advection terms for the variables
  ! located at momentum grid levels.  The variables are:  wprtp, wpthlp, wp2,
  ! rtp2, thlp2, rtpthlp, up2, vp2, wpsclrp, sclrprtp, sclrpthlp, and sclrp2.

  implicit none

  private ! Default scope

  public :: term_ma_zt_lhs, & 
            term_ma_zm_lhs

  contains

  !=============================================================================
  pure subroutine term_ma_zt_lhs( gr, wm_zt, invrs_dzt, invrs_dzm, & ! Intent(in)
                                  l_upwind_xm_ma,              & ! Intent(in)
                                  lhs_ma                       ) ! Intent(out)

    ! Description:
    ! Mean advection of var_zt:  implicit portion of the code.
    !
    ! The variable "var_zt" stands for a variable that is located at
    ! thermodynamic grid levels.
    !
    ! The d(var_zt)/dt equation contains a mean advection term:
    !
    ! - w * d(var_zt)/dz.
    !
    ! This term is solved for completely implicitly, such that:
    !
    ! - w * d( var_zt(t+1) )/dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign
    !        is reversed and the leading "-" in front of the term is changed to
    !        a "+".
    !
    ! The timestep index (t+1) means that the value of var_zt being used is from
    ! the next timestep, which is being advanced to in solving the d(var_zt)/dt
    ! equation.
    !
    ! This term is discretized as follows:
    !
    ! The values of var_zt are found on the thermodynamic levels, as are the
    ! values of wm_zt (mean vertical velocity on thermodynamic levels).  The
    ! variable var_zt is interpolated to the intermediate momentum levels.  The
    ! derivative of the interpolated values is taken over the central
    ! thermodynamic level.  The derivative is multiplied by wm_zt at the central
    ! thermodynamic level to get the desired result.
    !
    ! -----var_zt---------------------------------------------- t(k+1)
    !
    ! =================var_zt(interp)========================== m(k)
    !
    ! -----var_zt---------------------d(var_zt)/dz-----wm_zt--- t(k)
    !
    ! =================var_zt(interp)========================== m(k-1)
    !
    ! -----var_zt---------------------------------------------- t(k-1)
    !
    ! The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1) correspond
    ! with altitudes zt(k+1), zm(k), zt(k), zm(k-1), and zt(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is used
    ! for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) )
    !
    !
    ! Special discretization for upper boundary level:
    !
    ! Method 1:  Constant derivative method (or "one-sided" method).
    !
    ! The values of var_zt are found on the thermodynamic levels, as are the
    ! values of wm_zt (mean vertical velocity on the thermodynamic levels).  The
    ! variable var_zt is interpolated to momentum level gr%nz-1, based on
    ! the values of var_zt at thermodynamic levels gr%nz and gr%nz-1.
    ! However, the variable var_zt cannot be interpolated to momentum level
    ! gr%nz.  Rather, a linear extension is used to find the value of var_zt
    ! at momentum level gr%nz, based on the values of var_zt at thermodynamic
    ! levels gr%nz and gr%nz-1.  The derivative of the extended and
    ! interpolated values, d(var_zt)/dz, is taken over the central thermodynamic
    ! level.  Of course, this derivative will be the same as the derivative of
    ! var_zt between thermodynamic levels gr%nz and gr%nz-1.  The derivative
    ! is multiplied by wm_zt at the central thermodynamic level to get the
    ! desired result.
    !
    ! For the following diagram, k = gr%nz, which is the uppermost level of
    ! the model:
    !
    ! =================var_zt(extend)========================== m(k)   Boundary
    !
    ! -----var_zt---------------------d(var_zt)/dz-----wm_zt--- t(k)
    !
    ! =================var_zt(interp)========================== m(k-1)
    !
    ! -----var_zt---------------------------------------------- t(k-1)
    !
    !
    ! Method 2:  Zero derivative method:
    !            the derivative d(var_zt)/dz over the model top is set to 0.
    !
    ! This method corresponds with the "zero-flux" boundary condition option
    ! for eddy diffusion, where d(var_zt)/dz is set to 0 across the upper
    ! boundary.
    !
    ! In order to discretize the upper boundary condition, consider a new level
    ! outside the model (thermodynamic level gr%nz+1) just above the upper
    ! boundary level (thermodynamic level gr%nz).  The value of var_zt at the
    ! level just outside the model is defined to be the same as the value of
    ! var_zt at thermodynamic level gr%nz.  Therefore, the value of
    ! d(var_zt)/dz between the level just outside the model and the uppermost
    ! thermodynamic level is 0, staying consistent with the zero-flux boundary
    ! condition option for the eddy diffusion portion of the code.  Therefore,
    ! the value of var_zt at momentum level gr%nz, which is the upper boundary
    ! of the model, would be the same as the value of var_zt at the uppermost
    ! thermodynamic level.
    !
    ! The values of var_zt are found on the thermodynamic levels, as are the
    ! values of wm_zt (mean vertical velocity on the thermodynamic levels).  The
    ! variable var_zt is interpolated to momentum level gr%nz-1, based on
    ! the values of var_zt at thermodynamic levels gr%nz and gr%nz-1.  The
    ! value of var_zt at momentum level gr%nz is set equal to the value of
    ! var_zt at thermodynamic level gr%nz, as described above.  The derivative
    ! of the set and interpolated values, d(var_zt)/dz, is taken over the
    ! central thermodynamic level.  The derivative is multiplied by wm_zt at the
    ! central thermodynamic level to get the desired result.
    !
    ! For the following diagram, k = gr%nz, which is the uppermost level of
    ! the model:
    !
    ! --[var_zt(k+1) = var_zt(k)]----(level outside model)----- t(k+1)
    !
    ! ==[var_zt(top) = var_zt(k)]===[d(var_zt)/dz|_(top) = 0]== m(k)   Boundary
    !
    ! -----var_zt(k)------------------d(var_zt)/dz-----wm_zt--- t(k)
    !
    ! =================var_zt(interp)========================== m(k-1)
    !
    ! -----var_zt(k-1)----------------------------------------- t(k-1)
    !
    ! where (top) stands for the grid index of momentum level k = gr%nz, which
    ! is the upper boundary of the model.
    !
    ! This method of boundary discretization is also similar to the method
    ! currently employed at the lower boundary for most thermodynamic-level
    ! variables.  Since thermodynamic level k = 1 is below the model bottom,
    ! mean advection is not applied.  Thus, thermodynamic level k = 2 becomes
    ! the lower boundary level.  Now, the mean advection term at thermodynamic
    ! level 2 takes into account var_zt from levels 1, 2, and 3.  However, in
    ! most cases, the value of var_zt(1) is set equal to var_zt(2) after the
    ! matrix of equations has been solved.  Therefore, the derivative,
    ! d(var_zt)/dz, over the model bottom (momentum level k = 1) becomes 0.
    ! Thus, the method of setting d(var_zt)/dz to 0 over the model top keeps
    ! the way the upper and lower boundaries are handled consistent with each
    ! other.

    ! References:
    !   None
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        grid ! Type

    use constants_clubb, only: &
        one,  & ! Constant(s)
        zero

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1, & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2, & ! Thermodynamic main diagonal index.
      km1_tdiag = 3    ! Thermodynamic subdiagonal index.

    integer, parameter :: & 
      t_above = 1, & ! Index for upper thermodynamic level grid weight.
      t_below = 2    ! Index for lower thermodynamic level grid weight.

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      wm_zt,     & ! wm_zt                        [m/s]
      invrs_dzt, & ! Inverse of grid spacing      [1/m]
      invrs_dzm    ! Inverse of grid spacing      [1/m]

    logical, intent(in) :: &
      l_upwind_xm_ma    ! This flag determines whether we want to use an upwind
                        ! differencing approximation rather than a centered
                        ! differencing for turbulent or mean advection terms.
                        ! It affects rtm, thlm, sclrm, um and vm.

    ! Return Variable
    real( kind = core_rknd ), dimension(3,gr%nz), intent(out) :: &
      lhs_ma    ! Mean advection contributions to lhs    [1/s]

    ! Local Variables

    integer :: k    ! Vertical level index

    ! Set lower boundary array to 0
    lhs_ma(:,1) = 0.0_core_rknd


    if( .not. l_upwind_xm_ma ) then  ! Use centered differencing

       ! Most of the interior model; normal conditions.
       do k = 2, gr%nz, 1

          ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
          lhs_ma(kp1_tdiag,k) & 
          = + wm_zt(k) * invrs_dzt(k) * gr%weights_zt2zm(t_above,k)

          ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs_ma(k_tdiag,k) & 
          = + wm_zt(k) * invrs_dzt(k) * ( gr%weights_zt2zm(t_below,k) & 
                                          - gr%weights_zt2zm(t_above,k-1) )

          ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
          lhs_ma(km1_tdiag,k) & 
          = - wm_zt(k) * invrs_dzt(k) * gr%weights_zt2zm(t_below,k-1)

       enddo ! k = 2, gr%nz, 1

       ! Upper Boundary

        ! Special discretization for zero derivative method, where the
        ! derivative d(var_zt)/dz over the model top is set to 0, in order
        ! to stay consistent with the zero-flux boundary condition option
        ! in the eddy diffusion code.
        ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
        lhs_ma(kp1_tdiag,gr%nz) & 
        = zero
        ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
        lhs_ma(k_tdiag,gr%nz) & 
        = + wm_zt(gr%nz) &
            * invrs_dzt(gr%nz) * ( one - gr%weights_zt2zm(t_above,gr%nz-1) )

        ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
        lhs_ma(km1_tdiag,gr%nz) & 
        = - wm_zt(gr%nz) &
            * invrs_dzt(gr%nz) * gr%weights_zt2zm(t_below,gr%nz-1)


    else ! l_upwind_xm_ma == .true.; use "upwind" differencing

       ! Most of the interior model; normal conditions.
       do k = 2, gr%nz, 1

          if ( wm_zt(k) >= zero ) then  ! Mean wind is in upward direction

             ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
             lhs_ma(kp1_tdiag,k) &
             = zero

             ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
             lhs_ma(k_tdiag,k) &
             = + wm_zt(k) * invrs_dzm(k-1)

             ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
             lhs_ma(km1_tdiag,k) &
             = - wm_zt(k) * invrs_dzm(k-1)

          else  ! wm_zt < 0; Mean wind is in downward direction

             ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
             lhs_ma(kp1_tdiag,k) &
             = + wm_zt(k) * invrs_dzm(k)

             ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
             lhs_ma(k_tdiag,k) &
             = - wm_zt(k) * invrs_dzm(k)

             ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
             lhs_ma(km1_tdiag,k) &
             = zero

          endif ! wm_zt > 0

       enddo ! k = 2, gr%nz, 1

       ! Upper Boundary
       if ( wm_zt(gr%nz) >= zero ) then  ! Mean wind is in upward direction

          ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
          lhs_ma(kp1_tdiag,gr%nz) &
          = zero

          ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs_ma(k_tdiag,gr%nz) &
          = + wm_zt(gr%nz) * invrs_dzm(gr%nz-1)

          ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
          lhs_ma(km1_tdiag,gr%nz) &
          = - wm_zt(gr%nz) * invrs_dzm(gr%nz-1)

       else  ! wm_zt < 0; Mean wind is in downward direction

          ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
          lhs_ma(kp1_tdiag,gr%nz) &
          = zero

          ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs_ma(k_tdiag,gr%nz) &
          = zero

          ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
          lhs_ma(km1_tdiag,gr%nz) &
          = zero

       endif ! wm_zt > 0

    endif ! l_upwind_xm_ma


    return

  end subroutine term_ma_zt_lhs

  !=============================================================================
  pure subroutine term_ma_zm_lhs( gr, wm_zm, invrs_dzm, & 
                                  lhs_ma )

    ! Description:
    ! Mean advection of var_zm:  implicit portion of the code.
    !
    ! The variable "var_zm" stands for a variable that is located at momentum
    ! grid levels.
    !
    ! The d(var_zm)/dt equation contains a mean advection term:
    !
    ! - w * d(var_zm)/dz.
    !
    ! This term is solved for completely implicitly, such that:
    !
    ! - w * d( var_zm(t+1) )/dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign
    !        is reversed and the leading "-" in front of the term is changed to
    !        a "+".
    !
    ! The timestep index (t+1) means that the value of var_zm being used is from
    ! the next timestep, which is being advanced to in solving the d(var_zm)/dt
    ! equation.
    !
    ! This term is discretized as follows:
    !
    ! The values of var_zm are found on the momentum levels, as are the values
    ! of wm_zm (mean vertical velocity on momentum levels).  The variable var_zm
    ! is interpolated to the intermediate thermodynamic levels.  The derivative
    ! of the interpolated values is taken over the central momentum level.  The
    ! derivative is multiplied by wm_zm at the central momentum level to get the
    ! desired result.
    !
    ! =====var_zm============================================== m(k+1)
    !
    ! -----------------var_zm(interp)-------------------------- t(k+1)
    !
    ! =====var_zm=====================d(var_zm)/dz=====wm_zm=== m(k)
    !
    ! -----------------var_zm(interp)-------------------------- t(k)
    !
    ! =====var_zm============================================== m(k-1)
    !
    ! The vertical indices m(k+1), t(k+1), m(k), t(k), and m(k-1) correspond
    ! with altitudes zm(k+1), zt(k+1), zm(k), zt(k), and zm(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is used
    ! for momentum levels.
    !
    ! invrs_dzm(k) = 1 / ( zt(k+1) - zt(k) )

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        grid ! Type

    use constants_clubb, only: &
        zero ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Constant parameters
    integer, parameter :: & 
      kp1_mdiag = 1, & ! Momentum superdiagonal index.
      k_mdiag   = 2, & ! Momentum main diagonal index.
      km1_mdiag = 3    ! Momentum subdiagonal index.

    integer, parameter :: & 
      m_above = 1, & ! Index for upper momentum level grid weight.
      m_below = 2    ! Index for lower momentum level grid weight.

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      wm_zm,     & ! wm_zm                        [m/s]
      invrs_dzm    ! Inverse of grid spacing      [1/m]

    ! Return Variable
    real( kind = core_rknd ), dimension(3,gr%nz), intent(out) :: &
      lhs_ma    ! Mean advection contributions to lhs  [1/s]

    integer :: k    ! Vertical level index


    ! Set lower boundary array to 0
    lhs_ma(:,1) = zero

    ! Most of the interior model; normal conditions.
    do k = 2, gr%nz-1, 1

       ! Momentum superdiagonal: [ x var_zm(k+1,<t+1>) ]
       lhs_ma(kp1_mdiag,k) & 
       = + wm_zm(k) * invrs_dzm(k) * gr%weights_zm2zt(m_above,k+1)

       ! Momentum main diagonal: [ x var_zm(k,<t+1>) ]
       lhs_ma(k_mdiag,k) & 
       = + wm_zm(k) * invrs_dzm(k) * ( gr%weights_zm2zt(m_below,k+1) & 
                                       - gr%weights_zm2zt(m_above,k) )

       ! Momentum subdiagonal: [ x var_zm(k-1,<t+1>) ]
       lhs_ma(km1_mdiag,k) & 
       = - wm_zm(k) * invrs_dzm(k) * gr%weights_zm2zt(m_below,k)

    enddo ! k = 2, gr%nz-1, 1

    ! Set upper boundary array to 0
    lhs_ma(:,gr%nz) = zero


    return

  end subroutine term_ma_zm_lhs

  !=============================================================================

end module mean_adv
