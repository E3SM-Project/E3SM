! $Id: diffusion.F90 5623 2012-01-17 17:55:26Z connork@uwm.edu $
!===============================================================================
module diffusion

  ! Description:
  ! Module diffusion computes the eddy diffusion terms for all of the
  ! time-tendency (prognostic) equations in the CLUBB parameterization.  Most of
  ! the eddy diffusion terms are solved for completely implicitly, and therefore
  ! become part of the left-hand side of their respective equations.  However,
  ! wp2 and wp3 have an option to use a Crank-Nicholson eddy diffusion scheme,
  ! which has both implicit and explicit components.
  !
  ! Function diffusion_zt_lhs handles the eddy diffusion terms for the variables
  ! located at thermodynamic grid levels.  These variables are:  wp3 and all
  ! hydrometeor species.  The variables um and vm also use the Crank-Nicholson
  ! eddy-diffusion scheme for their turbulent advection term.
  !
  ! Function diffusion_zm_lhs handles the eddy diffusion terms for the variables
  ! located at momentum grid levels.  The variables are:  wprtp, wpthlp, wp2,
  ! rtp2, thlp2, rtpthlp, up2, vp2, wpsclrp, sclrprtp, sclrpthlp, and sclrp2.

  implicit none

  private ! Default Scope

  public :: diffusion_zt_lhs, & 
            diffusion_cloud_frac_zt_lhs, & 
            diffusion_zm_lhs

  contains

  !=============================================================================
  pure function diffusion_zt_lhs( K_zm, K_zmm1, nu,  &
                                  invrs_dzmm1, invrs_dzm,  &
                                  invrs_dzt, level )  &
  result( lhs )

    ! Description:
    ! Vertical eddy diffusion of var_zt:  implicit portion of the code.
    !
    ! The variable "var_zt" stands for a variable that is located at
    ! thermodynamic grid levels.
    !
    ! The d(var_zt)/dt equation contains an eddy diffusion term:
    !
    ! + d [ ( K_zm + nu ) * d(var_zt)/dz ] / dz.
    !
    ! This term is usually solved for completely implicitly, such that:
    !
    ! + d [ ( K_zm + nu ) * d( var_zt(t+1) )/dz ] / dz.
    !
    ! However, when a Crank-Nicholson scheme is used, the eddy diffusion term
    ! has both implicit and explicit components, such that:
    !
    ! + (1/2) * d [ ( K_zm + nu ) * d( var_zt(t+1) )/dz ] / dz
    !    + (1/2) * d [ ( K_zm + nu ) * d( var_zt(t) )/dz ] / dz;
    !
    ! for which the implicit component is:
    !
    ! + (1/2) * d [ ( K_zm + nu ) * d( var_zt(t+1) )/dz ] / dz.
    !
    ! Note:  When the implicit term is brought over to the left-hand side,
    !        the  sign is reversed and the leading "+" in front of the term
    !        is changed to a "-".
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is
    ! being advanced to in solving the d(var_zt)/dt equation.
    !
    ! The implicit portion of this term is discretized as follows:
    !
    ! The values of var_zt are found on the thermodynamic levels, while the
    ! values of K_zm are found on the momentum levels.  The derivatives (d/dz)
    ! of var_zt are taken over the intermediate momentum levels.  At the
    ! intermediate momentum levels, d(var_zt)/dz is multiplied by ( K_zm + nu ).
    ! Then, the derivative of the whole mathematical expression is taken over
    ! the central thermodynamic level, which yields the desired result.
    !
    ! --var_ztp1----------------------------------------------- t(k+1)
    !
    ! ==========d(var_zt)/dz==(K_zm+nu)======================== m(k)
    !
    ! --var_zt-------------------d[(K_zm+nu)*d(var_zt)/dz]/dz-- t(k)
    !
    ! ==========d(var_zt)/dz==(K_zmm1+nu)====================== m(k-1)
    !
    ! --var_ztm1----------------------------------------------- t(k-1)
    !
    ! The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1) correspond
    ! with altitudes zt(k+1), zm(k), zt(k), zm(k-1), and zt(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is used
    ! for momentum levels.
    !
    ! invrs_dzt(k)   = 1 / ( zm(k) - zm(k-1) )
    ! invrs_dzm(k)   = 1 / ( zt(k+1) - zt(k) )
    ! invrs_dzm(k-1) = 1 / ( zt(k) - zt(k-1) )
    !
    ! Note:  This function only computes the general implicit form:
    !        + d [ ( K_zm + nu ) * d( var_zt(t+1) )/dz ] / dz.
    !        For a Crank-Nicholson scheme, the left-hand side result of this
    !        function will have to be multiplied by (1/2).  For a
    !        Crank-Nicholson scheme, the right-hand side (explicit) component
    !        needs to be computed by multiplying the left-hand side results by
    !        (1/2), reversing the sign on each left-hand side element, and then
    !        multiplying each element by the appropriate var_zt(t) value from
    !        the appropriate vertical level.
    !
    !
    ! Boundary Conditions:
    !
    ! 1) Zero-flux boundary conditions.
    !    This function is set up to use zero-flux boundary conditions at both
    !    the lower boundary level and the upper boundary level.  The flux, F,
    !    is the amount of var_zt flowing normal through the boundary per unit
    !    time per unit surface area.  The derivative of the flux effects the
    !    time-tendency of var_zt, such that:
    !
    !    d(var_zt)/dt = -dF/dz.
    !
    !    For the 2nd-order eddy-diffusion term, +d[(K_zm+nu)*d(var_zt)/dz]/dz,
    !    the flux is:
    !
    !    F = -(K_zm+nu)*d(var_zt)/dz.
    !
    !    In order to have zero-flux boundary conditions, the derivative of
    !    var_zt, d(var_zt)/dz, needs to equal 0 at both the lower boundary and
    !    the upper boundary.
    !
    !    In order to discretize the lower boundary condition, consider a new
    !    level outside the model (thermodynamic level 0) just below the lower
    !    boundary level (thermodynamic level 1).  The value of var_zt at the
    !    level just outside the model is defined to be the same as the value of
    !    var_zt at the lower boundary level.  Therefore, the value of
    !    d(var_zt)/dz between the level just outside the model and the lower
    !    boundary level is 0, satisfying the zero-flux boundary condition.  The
    !    other value for d(var_zt)/dz (between thermodynamic level 2 and
    !    thermodynamic level 1) is taken over the intermediate momentum level
    !    (momentum level 1), where it is multiplied by the factor
    !    ( K_zm(1) + nu ).  Then, the derivative of the whole expression is
    !    taken over the central thermodynamic level.
    !
    !    -var_zt(2)-------------------------------------------- t(2)
    !
    !    ==========d(var_zt)/dz==(K_zm(1)+nu)================== m(1)
    !
    !    -var_zt(1)---------------d[(K_zm+nu)*d(var_zt)/dz]/dz- t(1) Boundary
    !
    !              [d(var_zt)/dz = 0]
    !
    !    -[var_zt(0) = var_zt(1)]-----(level outside model)---- t(0)
    !
    !    The result is dependent only on values of K_zm found at momentum
    !    level 1 and values of var_zt found at thermodynamic levels 1 and 2.
    !    Thus, it only affects 2 diagonals on the left-hand side matrix.
    !
    !    The same method can be used to discretize the upper boundary by
    !    considering a new level outside the model just above the upper boundary
    !    level.
    !
    ! 2) Fixed-point boundary conditions.
    !    Many equations in the model use fixed-point boundary conditions rather
    !    than zero-flux boundary conditions.  This means that the value of
    !    var_zt stays the same over the course of the timestep at the lower
    !    boundary, as well as at the upper boundary.
    !
    !    In order to discretize the boundary conditions for equations requiring
    !    fixed-point boundary conditions, either:
    !    a) in the parent subroutine or function (that calls this function),
    !       loop over all vertical levels from the second-lowest to the
    !       second-highest, ignoring the boundary levels.  Then set the values
    !       at the boundary levels in the parent subroutine; or
    !    b) in the parent subroutine or function, loop over all vertical levels
    !       and then overwrite the results at the boundary levels.
    !
    !    Either way, at the boundary levels, an array with a value of 1 at the
    !    main diagonal on the left-hand side and with values of 0 at all other
    !    diagonals on the left-hand side will preserve the right-hand side value
    !    at that level, thus satisfying the fixed-point boundary conditions.
    !
    !
    ! Conservation Properties:
    !
    ! When zero-flux boundary conditions are used, this technique of
    ! discretizing the eddy diffusion term leads to conservative differencing.
    ! When conservative differencing is in place, the column totals for each
    ! column in the left-hand side matrix (for the eddy diffusion term) should
    ! be equal to 0.  This ensures that the total amount of the quantity var_zt
    ! over the entire vertical domain is being conserved, meaning that nothing
    ! is lost due to diffusional effects.
    !
    ! To see that this conservation law is satisfied, compute the eddy diffusion
    ! of var_zt and integrate vertically.  In discretized matrix notation (where
    ! "i" stands for the matrix column and "j" stands for the matrix row):
    !
    !  0 = Sum_j Sum_i ( 1/invrs_dzt )_i 
    !                     ( invrs_dzt * ((K_zm+nu)*invrs_dzm) )_ij (var_zt)_j.
    !
    ! The left-hand side matrix, ( invrs_dzt * ((K_zm+nu)*invrs_dzm) )_ij, is
    ! partially written below.  The sum over i in the above equation removes
    ! invrs_dzt everywhere from the matrix below.  The sum over j leaves the
    ! column totals that are desired.
    !
    ! Left-hand side matrix contributions from eddy diffusion term; first four
    ! vertical levels:
    !
    !     -------------------------------------------------------------------------->
    !k=1 | +invrs_dzt(k)          -invrs_dzt(k)                          0
    !    |   *(K_zm(k)+nu)          *(K_zm(k)+nu)
    !    |     *invrs_dzm(k)          *invrs_dzm(k)
    !    |
    !k=2 | -invrs_dzt(k)          +invrs_dzt(k)              -invrs_dzt(k)
    !    |   *(K_zm(k-1)+nu)        *[ (K_zm(k)+nu)            *(K_zm(k)+nu)
    !    |     *invrs_dzm(k-1)          *invrs_dzm(k)            *invrs_dzm(k)
    !    |                            +(K_zm(k-1)+nu)
    !    |                              *invrs_dzm(k-1) ]
    !    |
    !k=3 |         0              -invrs_dzt(k)              +invrs_dzt(k)
    !    |                          *(K_zm(k-1)+nu)            *[ (K_zm(k)+nu)
    !    |                            *invrs_dzm(k-1)              *invrs_dzm(k)
    !    |                                                       +(K_zm(k-1)+nu)
    !    |                                                         *invrs_dzm(k-1) ]
    !    |
    !k=4 |         0                        0                -invrs_dzt(k)
    !    |                                                     *(K_zm(k-1)+nu)
    !    |                                                       *invrs_dzm(k-1)
    !   \ /
    !
    ! Note:  The superdiagonal term from level 3 and both the main diagonal and
    !        superdiagonal terms from level 4 are not shown on this diagram.
    !
    ! Note:  The matrix shown is a tridiagonal matrix.  For a band diagonal
    !        matrix (with 5 diagonals), there would be an extra row between each
    !        of the rows shown and an extra column between each of the columns
    !        shown.  However, for the purposes of the var_zt eddy diffusion
    !        term, those extra row and column values are all 0, and the
    !        conservation properties of the matrix aren't effected.
    !
    ! If fixed-point boundary conditions are used, the matrix entries at
    ! level 1 (k=1) read:  1   0   0; which means that conservative differencing
    ! is not in play.  The total amount of var_zt over the entire vertical
    ! domain is not being conserved, as amounts of var_zt may be fluxed out
    ! through the upper boundary or lower boundary through the effects of
    ! diffusion.
    !
    ! Brian Griffin.  April 26, 2008.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2,    & ! Thermodynamic main diagonal index.
      km1_tdiag = 3       ! Thermodynamic subdiagonal index.

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      K_zm,        & ! Coef. of eddy diffusivity at momentum level (k)   [m^2/s]
      K_zmm1,      & ! Coef. of eddy diffusivity at momentum level (k-1) [m^2/s 
      invrs_dzt,   & ! Inverse of grid spacing over thermo. level (k)    [1/m]
      invrs_dzm,   & ! Inverse of grid spacing over momentum level (k)   [1/m]
      invrs_dzmm1    ! Inverse of grid spacing over momentum level (k-1) [1/m]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      nu             ! Background constant coef. of eddy diffusivity     [m^2/s]

    integer, intent(in) ::  & 
      level     ! Thermodynamic level where calculation occurs.          [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(3) :: lhs

    if ( level == 1 ) then

      ! k = 1 (bottom level); lower boundary level.
      ! Only relevant if zero-flux boundary conditions are used.

      ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
      lhs(kp1_tdiag) = - invrs_dzt * (K_zm+nu(1)) * invrs_dzm

      ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
      lhs(k_tdiag)   = + invrs_dzt * (K_zm+nu(1)) * invrs_dzm

      ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
      lhs(km1_tdiag) = 0.0_core_rknd


    elseif ( level > 1 .and. level < gr%nz ) then

      ! Most of the interior model; normal conditions.

      ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
      lhs(kp1_tdiag) = - invrs_dzt * (K_zm+nu(level)) * invrs_dzm

      ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
      lhs(k_tdiag)   = + invrs_dzt * (  (K_zm+nu(level))*invrs_dzm  &
                                      + (K_zmm1+nu(level))*invrs_dzmm1 )

      ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
      lhs(km1_tdiag) = - invrs_dzt * (K_zmm1+nu(level)) * invrs_dzmm1

    elseif ( level == gr%nz ) then

      ! k = gr%nz (top level); upper boundary level.
      ! Only relevant if zero-flux boundary conditions are used.

      ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
      lhs(kp1_tdiag) = 0.0_core_rknd

      ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
      lhs(k_tdiag)   = + invrs_dzt * (K_zmm1+nu(gr%nz)) * invrs_dzmm1

      ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
      lhs(km1_tdiag) = - invrs_dzt * (K_zmm1+nu(gr%nz)) * invrs_dzmm1


    endif

  end function diffusion_zt_lhs

  !=============================================================================
  pure function diffusion_cloud_frac_zt_lhs &
                ( K_zm, K_zmm1, cloud_frac_zt, cloud_frac_ztm1, &
                  cloud_frac_ztp1, cloud_frac_zm, &
                  cloud_frac_zmm1, &
                  nu, invrs_dzmm1, invrs_dzm, invrs_dzt, level ) &
  result( lhs )

  ! Description:
  !   This function adds a weight of cloud fraction to the existing diffusion
  !   function for number concentration variables (e.g. cloud droplet number
  !   concentration).  This code should be considered experimental and may
  !   contain bugs.
  ! References:
  !   This algorithm uses equations derived from Guo, et al. 2009.
  !-----------------------------------------------------------------------------

   use grid_class, only: & 
     gr ! Variable(s)

   use clubb_precision, only: &
     core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: min

    ! Constant parameters
    real( kind = core_rknd ), parameter :: &
      cf_ratio = 10._core_rknd ! Maximum cloud-fraction coefficient applied to Kh_zm

    integer, parameter :: & 
      kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2,    & ! Thermodynamic main diagonal index.
      km1_tdiag = 3       ! Thermodynamic subdiagonal index.

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      K_zm,            & ! Coef. of eddy diffusivity at mom. level (k)   [m^2/s]
      K_zmm1,          & ! Coef. of eddy diffusivity at mom. level (k-1) [m^2/s]
      cloud_frac_zt,   & ! Cloud fraction at the thermo. level (k)       [-]
      cloud_frac_ztm1, & ! Cloud fraction at the thermo. level (k-1)     [-]
      cloud_frac_ztp1, & ! Cloud fraction at the thermo. level (k+1)     [-]
      cloud_frac_zm,   & ! Cloud fraction at the momentum level (k)      [-]
      cloud_frac_zmm1, & ! Cloud fraction at the momentum level (k-1)    [-]
      invrs_dzt,       & ! Inverse of grid spacing over thermo. lev. (k) [1/m]
      invrs_dzm,       & ! Inverse of grid spacing over mom. level (k)   [1/m]
      invrs_dzmm1        ! Inverse of grid spacing over mom. level (k-1) [1/m]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      nu                 ! Background constant coef. of eddy diffusivity [m^2/s]

    integer, intent(in) ::  & 
      level     ! Thermodynamic level where calculation occurs.           [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(3) :: lhs

    ! ---- Begin Code ----

    if ( level == 1 ) then

      ! k = 1 (bottom level); lower boundary level.
      ! Only relevant if zero-flux boundary conditions are used.

      ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
!     lhs(kp1_tdiag) = - invrs_dzt &
!                        * (K_zm+nu) &
!                           * ( cloud_frac_zm / cloud_frac_ztp1 ) * invrs_dzm
      lhs(kp1_tdiag) = - invrs_dzt &
                         * (K_zm &
                            * min( cloud_frac_zm / cloud_frac_ztp1, cf_ratio ) &
                            + nu(1)) * invrs_dzm

      ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
!     lhs(k_tdiag)   = + invrs_dzt &
!                        * (K_zm+nu) &
!                           * ( cloud_frac_zm / cloud_frac_ztp1 ) * invrs_dzm
      lhs(k_tdiag)   = + invrs_dzt &
                         * (K_zm &
                            * min( cloud_frac_zm / cloud_frac_ztp1, cf_ratio ) &
                            + nu(1)) * invrs_dzm

      ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
      lhs(km1_tdiag) = 0.0_core_rknd


    else if ( level > 1 .and. level < gr%nz ) then

      ! Most of the interior model; normal conditions.

      ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
!     lhs(kp1_tdiag) = - invrs_dzt &
!                        * (K_zm+nu) &
!                           * ( cloud_frac_zm / cloud_frac_ztp1 ) * invrs_dzm
!     lhs(kp1_tdiag) = - invrs_dzt &
!                        * (K_zm &
!                           * ( cloud_frac_zm / cloud_frac_ztp1 ) &
!                           + nu ) * invrs_dzm
      lhs(kp1_tdiag) = - invrs_dzt &
                         * (K_zm &
                            * min( cloud_frac_zm / cloud_frac_ztp1, cf_ratio ) &
                            + nu(level) ) * invrs_dzm

      ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
!     lhs(k_tdiag)   = + invrs_dzt &
!                        * (  ((K_zm+nu)*cloud_frac_zm)*invrs_dzm &
!                           + ((K_zmm1+nu)*cloud_frac_zmm1)*invrs_dzmm1 ) &
!                        / cloud_frac_zt
!     lhs(k_tdiag)   = + invrs_dzt &
!                        * ( nu*(invrs_dzm+invrs_dzmm1) + &
!                               ( ((K_zm*cloud_frac_zm)*invrs_dzm + 
!                                  (K_zmm1*cloud_frac_zmm1)*invrs_dzmm1)&
!                                / cloud_frac_zt &
!                               ) &
!                          )
      lhs(k_tdiag)   = + invrs_dzt &
                          * ( nu(level)*(invrs_dzm+invrs_dzmm1) + &
                              (   K_zm*invrs_dzm* &
                                     min( cloud_frac_zm / cloud_frac_zt, &
                                          cf_ratio ) &
                                + K_zmm1*invrs_dzmm1* &
                                     min( cloud_frac_zmm1 / cloud_frac_zt, &
                                          cf_ratio ) &
                              ) &
                            )

      ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
!     lhs(km1_tdiag) = - invrs_dzt * (K_zmm1+nu) * &
!                        ( cloud_frac_zmm1 / cloud_frac_ztm1 ) * invrs_dzmm1
      lhs(km1_tdiag) = - invrs_dzt &
                          * (K_zmm1 &
                             * min( cloud_frac_zmm1 / cloud_frac_ztm1, &
                                    cf_ratio ) &
                             + nu(level) ) * invrs_dzmm1

    else if ( level == gr%nz ) then

      ! k = gr%nz (top level); upper boundary level.
      ! Only relevant if zero-flux boundary conditions are used.

      ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
      lhs(kp1_tdiag) = 0.0_core_rknd

      ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
!     lhs(k_tdiag)   = + invrs_dzt &
!                         *(K_zmm1+nu) &
!                           *( cloud_frac_zmm1 / cloud_frac_ztm1 ) * invrs_dzmm1
      lhs(k_tdiag)   = + invrs_dzt &
                          * (K_zmm1 &
                             * min( cloud_frac_zmm1 / cloud_frac_ztm1, &
                                    cf_ratio ) &
                             + nu(gr%nz)) * invrs_dzmm1

      ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
!     lhs(km1_tdiag) = - invrs_dzt * (K_zmm1+nu) * &
!                        ( cloud_frac_zmm1 / cloud_frac_ztm1 ) * invrs_dzmm1
      lhs(km1_tdiag) = - invrs_dzt &
                          * (K_zmm1 &
                             * min( cloud_frac_zmm1 / cloud_frac_ztm1, &
                                    cf_ratio ) &
                             + nu(gr%nz)) * invrs_dzmm1

    end if

    return
  end function diffusion_cloud_frac_zt_lhs

  !=============================================================================
  pure function diffusion_zm_lhs( K_zt, K_ztp1, nu,  & 
                                  invrs_dztp1, invrs_dzt,  &
                                  invrs_dzm, level )  &
  result( lhs )

    ! Description:
    ! Vertical eddy diffusion of var_zm:  implicit portion of the code.
    !
    ! The variable "var_zm" stands for a variable that is located at momentum
    ! grid levels.
    !
    ! The d(var_zm)/dt equation contains an eddy diffusion term:
    !
    ! + d [ ( K_zt + nu ) * d(var_zm)/dz ] / dz.
    !
    ! This term is usually solved for completely implicitly, such that:
    !
    ! + d [ ( K_zt + nu ) * d( var_zm(t+1) )/dz ] / dz.
    !
    ! However, when a Crank-Nicholson scheme is used, the eddy diffusion term
    ! has both implicit and explicit components, such that:
    !
    ! + (1/2) * d [ ( K_zt + nu ) * d( var_zm(t+1) )/dz ] / dz
    !    + (1/2) * d [ ( K_zt + nu ) * d( var_zm(t) )/dz ] / dz;
    !
    ! for which the implicit component is:
    !
    ! + (1/2) * d [ ( K_zt + nu ) * d( var_zm(t+1) )/dz ] / dz.
    !
    ! Note:  When the implicit term is brought over to the left-hand side,
    !        the sign is reversed and the leading "+" in front of the term
    !        is changed to a "-".
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is
    ! being advanced to in solving the d(var_zm)/dt equation.
    !
    ! The implicit portion of this term is discretized as follows:
    !
    ! The values of var_zm are found on the momentum levels, while the values of
    ! K_zt are found on the thermodynamic levels.  The derivatives (d/dz) of
    ! var_zm are taken over the intermediate thermodynamic levels.  At the
    ! intermediate thermodynamic levels, d(var_zm)/dz is multiplied by
    ! ( K_zt + nu ).  Then, the derivative of the whole mathematical expression
    ! is taken over the central momentum level, which yields the desired result.
    !
    ! ==var_zmp1=============================================== m(k+1)
    !
    ! ----------d(var_zm)/dz--(K_ztp1+nu)---------------------- t(k+1)
    !
    ! ==var_zm===================d[(K_zt+nu)*d(var_zm)/dz]/dz== m(k)
    !
    ! ----------d(var_zm)/dz--(K_zt+nu)------------------------ t(k)
    !
    ! ==var_zmm1=============================================== m(k-1)
    !
    ! The vertical indices m(k+1), t(k+1), m(k), t(k), and m(k-1) correspond
    ! with altitudes zm(k+1), zt(k+1), zm(k), zt(k), and zm(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is used
    ! for momentum levels.
    !
    ! invrs_dzm(k)   = 1 / ( zt(k+1) - zt(k) )
    ! invrs_dzt(k+1) = 1 / ( zm(k+1) - zm(k) )
    ! invrs_dzt(k)   = 1 / ( zm(k) - zm(k-1) )
    !
    ! Note:  This function only computes the general implicit form:
    !        + d [ ( K_zt + nu ) * d( var_zm(t+1) )/dz ] / dz.
    !        For a Crank-Nicholson scheme, the left-hand side result of this
    !        function will have to be multiplied by (1/2).  For a
    !        Crank-Nicholson scheme, the right-hand side (explicit) component
    !        needs to be computed by multiplying the left-hand side results by
    !        (1/2), reversing the sign on each left-hand side element, and then
    !        multiplying each element by the appropriate var_zm(t) value from
    !        the appropriate vertical level.
    !
    !
    ! Boundary Conditions:
    !
    ! 1) Zero-flux boundary conditions.
    !    This function is set up to use zero-flux boundary conditions at both
    !    the lower boundary level and the upper boundary level.  The flux, F,
    !    is the amount of var_zm flowing normal through the boundary per unit
    !    time per unit surface area.  The derivative of the flux effects the
    !    time-tendency of var_zm, such that:
    !
    !    d(var_zm)/dt = -dF/dz.
    !
    !    For the 2nd-order eddy-diffusion term, +d[(K_zt+nu)*d(var_zm)/dz]/dz,
    !    the flux is:
    !
    !    F = -(K_zt+nu)*d(var_zm)/dz.
    !
    !    In order to have zero-flux boundary conditions, the derivative of
    !    var_zm, d(var_zm)/dz, needs to equal 0 at both the lower boundary and
    !    the upper boundary.
    !
    !    In order to discretize the lower boundary condition, consider a new
    !    level outside the model (momentum level 0) just below the lower
    !    boundary level (momentum level 1).  The value of var_zm at the level
    !    just outside the model is defined to be the same as the value of var_zm
    !    at the lower boundary level.  Therefore, the value of d(var_zm)/dz
    !    between the level just outside the model and the lower boundary level
    !    is 0, satisfying the zero-flux boundary condition.  The other value for
    !    d(var_zm)/dz (between momentum level 2 and momentum level 1) is taken
    !    over the intermediate thermodynamic level (thermodynamic level 2),
    !    where it is multiplied by the factor ( K_zt(2) + nu ).  Then, the
    !    derivative of the whole expression is taken over the central momentum
    !    level.
    !
    !    =var_zm(2)============================================ m(2)
    !
    !    ----------d(var_zm)/dz==(K_zt(2)+nu)------------------ t(2)
    !
    !    =var_zm(1)===============d[(K_zt+nu)*d(var_zm)/dz]/dz= m(1) Boundary
    !
    !    ----------[d(var_zm)/dz = 0]-------------------------- t(1)
    !
    !    =[var_zm(0) = var_zm(1)]=====(level outside model)==== m(0)
    !
    !    The result is dependent only on values of K_zt found at thermodynamic
    !    level 2 and values of var_zm found at momentum levels 1 and 2.  Thus,
    !    it only affects 2 diagonals on the left-hand side matrix.
    !
    !    The same method can be used to discretize the upper boundary by
    !    considering a new level outside the model just above the upper boundary
    !    level.
    !
    ! 2) Fixed-point boundary conditions.
    !    Many equations in the model use fixed-point boundary conditions rather
    !    than zero-flux boundary conditions.  This means that the value of
    !    var_zm stays the same over the course of the timestep at the lower
    !    boundary, as well as at the upper boundary.
    !
    !    In order to discretize the boundary conditions for equations requiring
    !    fixed-point boundary conditions, either:
    !    a) in the parent subroutine or function (that calls this function),
    !       loop over all vertical levels from the second-lowest to the
    !       second-highest, ignoring the boundary levels.  Then set the values
    !       at the boundary levels in the parent subroutine; or
    !    b) in the parent subroutine or function, loop over all vertical levels
    !       and then overwrite the results at the boundary levels.
    !
    !    Either way, at the boundary levels, an array with a value of 1 at the
    !    main diagonal on the left-hand side and with values of 0 at all other
    !    diagonals on the left-hand side will preserve the right-hand side value
    !    at that level, thus satisfying the fixed-point boundary conditions.
    !
    !
    ! Conservation Properties:
    !
    ! When zero-flux boundary conditions are used, this technique of
    ! discretizing the eddy diffusion term leads to conservative differencing.
    ! When conservative differencing is in place, the column totals for each
    ! column in the left-hand side matrix (for the eddy diffusion term) should
    ! be equal to 0.  This ensures that the total amount of the quantity var_zm
    ! over the entire vertical domain is being conserved, meaning that nothing
    ! is lost due to diffusional effects.
    !
    ! To see that this conservation law is satisfied, compute the eddy diffusion
    ! of var_zm and integrate vertically.  In discretized matrix notation (where
    ! "i" stands for the matrix column and "j" stands for the matrix row):
    !
    !  0 = Sum_j Sum_i ( 1/invrs_dzm )_i
    !                     ( invrs_dzm * ((K_zt+nu)*invrs_dzt) )_ij (var_zm)_j.
    !
    ! The left-hand side matrix, ( invrs_dzm * ((K_zt+nu)*invrs_dzt) )_ij, is
    ! partially written below.  The sum over i in the above equation removes
    ! invrs_dzm everywhere from the matrix below.  The sum over j leaves the
    ! column totals that are desired.
    !
    ! Left-hand side matrix contributions from eddy diffusion term; first four
    ! vertical levels:
    !
    !     ---------------------------------------------------------------------->
    !k=1 | +invrs_dzm(k)          -invrs_dzm(k)                      0
    !    |   *(K_zt(k+1)+nu)        *(K_zt(k+1)+nu)
    !    |     *invrs_dzt(k+1)        *invrs_dzt(k+1)
    !    |
    !k=2 | -invrs_dzm(k)          +invrs_dzm(k)            -invrs_dzm(k)
    !    |   *(K_zt(k)+nu)          *[ (K_zt(k+1)+nu)        *(K_zt(k+1)+nu)
    !    |     *invrs_dzt(k)            *invrs_dzt(k+1)        *invrs_dzt(k+1)
    !    |                            +(K_zt(k)+nu)
    !    |                              *invrs_dzt(k) ]
    !    |
    !k=3 |          0             -invrs_dzm(k)            +invrs_dzm(k)
    !    |                          *(K_zt(k)+nu)            *[ (K_zt(k+1)+nu)
    !    |                            *invrs_dzt(k)              *invrs_dzt(k+1)
    !    |                                                     +(K_zt(k)+nu)
    !    |                                                       *invrs_dzt(k) ]
    !    |
    !k=4 |          0                       0              -invrs_dzm(k)
    !    |                                                   *(K_zt(k)+nu)
    !    |                                                     *invrs_dzt(k)
    !   \ /
    !
    ! Note:  The superdiagonal term from level 3 and both the main diagonal and
    !        superdiagonal terms from level 4 are not shown on this diagram.
    !
    ! Note:  The matrix shown is a tridiagonal matrix.  For a band diagonal
    !        matrix (with 5 diagonals), there would be an extra row between each
    !        of the rows shown and an extra column between each of the columns
    !        shown.  However, for the purposes of the var_zm eddy diffusion
    !        term, those extra row and column values are all 0, and the
    !        conservation properties of the matrix aren't effected.
    !
    ! If fixed-point boundary conditions are used, the matrix entries at
    ! level 1 (k=1) read:  1   0   0; which means that conservative differencing
    ! is not in play.  The total amount of var_zm over the entire vertical
    ! domain is not being conserved, as amounts of var_zm may be fluxed out
    ! through the upper boundary or lower boundary through the effects of
    ! diffusion.
    !
    ! Brian Griffin.  April 26, 2008.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        gr       ! Variable(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_mdiag = 1,    & ! Momentum superdiagonal index.
      k_mdiag   = 2,    & ! Momentum main diagonal index.
      km1_mdiag = 3       ! Momentum subdiagonal index.

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      K_zt,        & ! Coef. of eddy diffusivity at thermo. level (k)   [m^2/s]
      K_ztp1,      & ! Coef. of eddy diffusivity at thermo. level (k+1) [m^2/s]
      invrs_dzm,   & ! Inverse of grid spacing over momentum level (k)  [1/m]
      invrs_dzt,   & ! Inverse of grid spacing over thermo. level (k)   [1/m]
      invrs_dztp1    ! Inverse of grid spacing over thermo. level (k+1) [1/m]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      nu             ! Background constant coef. of eddy diffusivity    [m^2/s]

    integer, intent(in) ::  & 
      level     ! Momentum level where calculation occurs.              [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(3) :: lhs

    if ( level == 1 ) then

      ! k = 1; lower boundary level at surface.
      ! Only relevant if zero-flux boundary conditions are used.

      ! Momentum superdiagonal: [ x var_zm(k+1,<t+1>) ]
      lhs(kp1_mdiag) = - invrs_dzm * (K_ztp1+nu(2)) * invrs_dztp1

      ! Momentum main diagonal: [ x var_zm(k,<t+1>) ]
      lhs(k_mdiag)   = + invrs_dzm * (K_ztp1+nu(2)) * invrs_dztp1

      ! Momentum subdiagonal: [ x var_zm(k-1,<t+1>) ]
      lhs(km1_mdiag) = 0.0_core_rknd


    elseif ( level > 1 .and. level < gr%nz ) then

      ! Most of the interior model; normal conditions.

      ! Momentum superdiagonal: [ x var_zm(k+1,<t+1>) ]
      lhs(kp1_mdiag) = - invrs_dzm * (K_ztp1+nu(level+1)) * invrs_dztp1

      ! Momentum main diagonal: [ x var_zm(k,<t+1>) ]
      lhs(k_mdiag)   = + invrs_dzm * (   (K_ztp1+nu(level+1))*invrs_dztp1  &
                                       + (K_zt+nu(level))*invrs_dzt )

      ! Momentum subdiagonal: [ x var_zm(k-1,<t+1>) ]
      lhs(km1_mdiag) = - invrs_dzm * (K_zt+nu(level)) * invrs_dzt


    elseif ( level == gr%nz ) then

      ! k = gr%nz (top level); upper boundary level.
      ! Only relevant if zero-flux boundary conditions are used.

      ! Momentum superdiagonal: [ x var_zm(k+1,<t+1>) ]
      lhs(kp1_mdiag) = 0.0_core_rknd

      ! Momentum main diagonal: [ x var_zm(k,<t+1>) ]
      lhs(k_mdiag)   = + invrs_dzm * (K_zt+nu(gr%nz)) * invrs_dzt

      ! Momentum subdiagonal: [ x var_zm(k-1,<t+1>) ]
      lhs(km1_mdiag) = - invrs_dzm * (K_zt+nu(gr%nz)) * invrs_dzt


    endif

  end function diffusion_zm_lhs

!===============================================================================

end module diffusion
