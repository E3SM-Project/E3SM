!-----------------------------------------------------------------------
! $Id: clip_semi_implicit.F90 7140 2014-07-31 19:14:05Z betlej@uwm.edu $
!===============================================================================
module clip_semi_implicit

  ! Description of the semi-implicit clipping code:
  ! The semi-implicit clipping code is based on an upper threshold and/or a
  ! lower threshold value for variable f.
  !
  ! The semi-implicit clipping code is used when the value of variable f should
  ! not exceed the designated threshold(s) when it is advanced to timestep
  ! index (t+1).
  !
  !
  ! Clipping at an Upper Threshold:
  !
  ! When there is an upper threshold to be applied, the equation for the clipped
  ! value of the variable f, f_clipped, is:
  !
  ! f_clipped(t+1) = MIN( f_unclipped(t+1), upper_threshold )
  !                = ( f_unclipped(t+1) - upper_threshold )
  !                    * H(upper_threshold-f_unclipped(t+1))
  !                  + upper_threshold;
  !
  ! where f_unclipped is the value of the variable f without clipping, and
  ! H(upper_threshold-f_unclipped(t+1)) is the Heaviside Step function.  The
  ! clipping term is turned into a time tendency term, such that:
  !
  ! (df/dt)_clipping = (1/dt_clip)
  !                    * ( f_clipped(t+1) - f_unclipped(t+1) );
  !
  ! where dt_clip is the time scale for the clipping term.  The difference
  ! between the threshold value and f_unclipped is defined as f_diff:
  !
  ! f_diff = upper_threshold - f_unclipped.
  !
  ! The clipping time tendency is now simplified as:
  !
  ! (df/dt)_clipping = + (1/dt_clip)
  !                      * { f_diff(t+1) * [ 1 - H(f_diff(t+1)) ] }.
  !
  ! Function R(f_diff) is defined as:
  !
  ! R(f_diff) = { f_diff * [ 1 - H(f_diff) ] }.
  !
  ! The clipping time tendency is now written as:
  !                
  ! (df/dt)_clipping = + (1/dt_clip) * R(f_diff(t+1)).
  !
  ! In order to solve for f_unclipped (and f_diff) at timestep index (t+1), the
  ! clipping term must be linearized.  A Taylor Series expansion (truncated
  ! after the first derivative term) of R(f_diff) around f_diff = f_diff(t) is
  ! used to linearize the term.  However, the Heaviside Step function,
  ! H(f_diff), is not differentiable when f_diff(t) = 0, as the function jumps
  ! at that point.  Likewise, the function R(f_diff) is not differentiable when
  ! f_diff(t) = 0, as the function has a corner at that point.  Therefore, a new
  ! function, F_R(f_diff) is used as an approximation of R(f_diff).  Function
  ! F_R(f_diff) is a three-piece function that has the exact same value as
  ! R(f_diff) when f_diff <= -sigma or f_diff >= sigma (sigma is an arbitrarily
  ! declared value).  However, when -sigma < f_diff < sigma, a parabolic
  ! function is used to approximate the corner found in R(f_diff).  The
  ! parabolic function needs to have the same values at f_diff = -sigma and
  ! f_diff = sigma as does R(f_diff).  Furthermore, the derivative of the
  ! parabolic function (with respect to f_diff) needs to have the same values at
  ! f_diff = -sigma and f_diff = sigma as does d(R)/d(f_diff).  The parabolic
  ! function that satisfies these properities is:
  ! f_diff - (sigma/4) * [ 1 + (f_diff/sigma) ]^2.
  ! Therefore:
  !
  !               | f_diff;   where f_diff <= -sigma
  !               |
  ! F_R(f_diff) = | f_diff - (sigma/4) * [ 1 + (f_diff/sigma) ]^2;
  !               |    where -sigma < f_diff < sigma
  !               |
  !               | 0;   where f_diff >= sigma;   and
  !
  !                        | 1;   where f_diff <= -sigma
  !                        |
  ! ( d F_R / d f_diff ) = | 1 - (1/2) * [ 1 + (f_diff/sigma) ];
  !                        |    where -sigma < f_diff < sigma
  !                        |
  !                        | 0;   where f_diff >= sigma.
  !
  ! Since, R(f_diff(t+1)) approx.= F_R(f_diff(t+1)), the Taylor Series expansion
  ! is done for F_R(f_diff) around f_diff = f_diff(t) in order to linearize the
  ! term:
  !
  ! F_R(f_diff(t+1)) approx.=
  !    A_fnc + B_fnc * ( f_diff(t+1) - f_diff(t) );
  !
  ! where A_fnc is defined as F_R(f_diff(t)) and B_fnc is defined as
  ! ( d F_R / d f_diff )|_(f_diff=f_diff(t)).
  !
  ! The approximation is substituted into the (df/dt)_clipping equation.  The
  ! rate of change of variable f due to clipping with the upper threshold is:
  !
  ! (df/dt)_clipping
  ! = + (1/dt_clip)
  !     * {   A_fnc - B_fnc * f_diff(t)
  !         + B_fnc * upper_threshold - B_fnc * f_unclipped(t+1) }.
  !
  ! The implicit (LHS) portion of the equation for clipping with the upper
  ! threshold is:
  !
  ! - (1/dt_clip) * B_fnc * f_unclipped(t+1).
  !
  ! Note:  When the term is brought over to the left-hand side, the sign
  !        is reversed and the leading "-" in front of the term is changed
  !        to a "+".
  !
  ! The explicit (RHS) portion of the equation for clipping with the upper
  ! threshold is:
  !
  ! + (1/dt_clip)
  !   * { A_fnc - B_fnc * f_diff(t) + B_fnc * upper_threshold }.
  !
  ! Timestep index (t) stands for the index of the current timestep, while
  ! timestep index (t+1) stands for the index of the next timestep, which is
  ! being advanced to in solving the d(f)/dt equation.
  !
  !
  ! Clipping at a Lower Threshold:
  !
  ! When there is a lower threshold to be applied, the equation for the clipped
  ! value of the variable f, f_clipped, is:
  !
  ! f_clipped(t+1) = MAX( f_unclipped(t+1), lower_threshold )
  !                = ( f_unclipped(t+1) - lower_threshold )
  !                    * H(f_unclipped(t+1)-lower_threshold)
  !                  + lower_threshold;
  !
  ! where f_unclipped is the value of the variable f without clipping, and
  ! H(f_unclipped(t+1)-lower_threshold) is the Heaviside Step function.  The
  ! clipping term is turned into a time tendency term, such that:
  !
  ! (df/dt)_clipping = (1/dt_clip)
  !                    * ( f_clipped(t+1) - f_unclipped(t+1) );
  !
  ! where dt_clip is the time scale for the clipping term.  The difference
  ! between f_unclipped and the threshold value is defined as f_diff:
  !
  ! f_diff = f_unclipped - lower_threshold.
  !
  ! The clipping time tendency is now simplified as:
  !
  ! (df/dt)_clipping = - (1/dt_clip)
  !                      * { f_diff(t+1) * [ 1 - H(f_diff(t+1)) ] }.
  !
  ! Function R(f_diff) is defined as:
  !
  ! R(f_diff) = { f_diff * [ 1 - H(f_diff) ] }.
  !
  ! The clipping time tendency is now written as:
  !
  ! (df/dt)_clipping = - (1/dt_clip) * R(f_diff(t+1)).
  !
  ! The linearization process is the same for the lower threshold as it is for
  ! the upper threshold.  The formulas for A_fnc and B_fnc are the same, but the
  ! values (based on a different f_diff) are different.  The rate of change of
  ! variable f due to clipping with the lower threshold is:
  !
  ! (df/dt)_clipping
  ! = - (1/dt_clip)
  !     * {   A_fnc - B_fnc * f_diff(t)
  !         - B_fnc * lower_threshold + B_fnc * f_unclipped(t+1) }.
  !
  ! The implicit (LHS) portion of the equation for clipping with the lower
  ! threshold is:
  !
  ! - (1/dt_clip) * B_fnc * f_unclipped(t+1).
  !
  ! Note:  When the term is brought over to the left-hand side, the sign
  !        is reversed and the leading "-" in front of the term is changed
  !        to a "+".
  !
  ! The explicit (RHS) portion of the equation for clipping with the lower
  ! threshold is:
  !
  ! - (1/dt_clip)
  !   * { A_fnc - B_fnc * f_diff(t) - B_fnc * lower_threshold }.
  !
  ! All variables in these equations are on the same vertical levels as the
  ! variable f.
  !
  !
  ! Adjustable parameters:
  !
  ! sigma:  sigma is the amount on either side of the threshold value to which
  !         the parabolic function portion of F_R(f_diff) is applied.  The value
  !         of sigma must be greater than 0.  A proportionally larger value of
  !         sigma can be used to effect values of f that are near the threshold,
  !         but not to it or over it.  The close-to-threshold values will be
  !         nudged away from the threshold.
  !
  ! dt_clip:  dt_clip is the clipping time scale.  It can be set equal to the
  !           model timestep, dt, but it doesn't have to be.  Smaller values of
  !           dt_clip produce a greater effect on the clipping term.

  ! References:
  !-----------------------------------------------------------------------

  use clubb_precision, only:  & 
    core_rknd ! Variable(s)

  implicit none

  private

  public :: clip_semi_imp_lhs, & 
            clip_semi_imp_rhs

  private :: compute_clip_lhs, & 
             compute_fncts_A_B

  ! Constant parameters.

  ! sigma coefficient:  A coefficient with dimensionless units that must have a
  !                     value greater than 0.  The value should be kept below 1.
  !                     The larger the value of sigma_coef, the larger the value
  !                     of sigma, and the larger the range of close-to-threshold
  !                     values that will be effected (nudged away from the
  !                     threshold) by the semi-implicit clipping.
  real( kind = core_rknd ), parameter :: sigma_coef = 0.15_core_rknd

  ! dt_clip coefficient:  A coefficient with dimensionless units that must have
  !                       a value greater than 0.  A value of 1 will set the
  !                       clipping time scale, dt_clip, equal to the model
  !                       timestep, dt.  The smaller the value of dt_clip_coef,
  !                       the smaller the value of dt_clip, and the larger the
  !                       magnitude of (df/dt)_clipping.
  real(kind=core_rknd), parameter :: dt_clip_coef = 1.0_core_rknd

  contains

  !=============================================================================
  function clip_semi_imp_lhs( dt, f_unclipped,  & 
                              l_upper_thresh, upper_threshold,  & 
                              l_lower_thresh, lower_threshold ) & 
  result( lhs )

    ! Description:
    ! The implicit portion of the semi-implicit clipping code.
    !
    ! The implicit (LHS) portion of the equation for clipping with the upper
    ! threshold is:
    !
    ! - (1/dt_clip) * B_fnc * f_unclipped(t+1).
    !
    ! The implicit (LHS) portion of the equation for clipping with the lower
    ! threshold is:
    !
    ! - (1/dt_clip) * B_fnc * f_unclipped(t+1).
    !
    ! Note:  When either term is brought over to the left-hand side, the sign
    !        is reversed and the leading "-" in front of either term is changed
    !        to a "+".
    !
    ! The timestep index (t+1) means that the value of f being used is from the
    ! next timestep, which is being advanced to in solving the d(f)/dt equation.
    !
    ! While the formulas are the same for both the upper threshold and the lower
    ! threshold, the values of A_fnc, B_fnc, and f_diff will differ between the
    ! two thresholds.
    !
    ! The overall implicit (LHS) portion for the clipping term is the sum of the
    ! implicit portion from the upper threshold and the implicit portion from
    ! the lower threshold.

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)implicit none

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      dt  ! Model timestep.  [s]

    real( kind = core_rknd ), intent(in) :: & 
      f_unclipped,        & ! The unclipped value of variable f at timestep (t). [f units]
      upper_threshold,    & ! Greatest allowable value of variable f.            [f units]
      lower_threshold       ! Smallest allowable value of variable f.            [f units]

    logical, intent(in) :: & 
      l_upper_thresh,    & ! Flag for having an upper threshold value.
      l_lower_thresh       ! Flag for having a lower threshold value.

    ! Return Variable
    real( kind = core_rknd ) :: lhs

    ! Local Variables
    real( kind = core_rknd ) ::  & 
      dt_clip       ! Time scale for semi-implicit clipping term.                [s]

    real( kind = core_rknd ) :: & 
      f_diff,     & ! Difference between the threshold value and f_unclipped.    [f units]
      A_fnc,      & ! Function that approximates { f_diff * [ 1 - H(f_diff) ] }. [f units]
      B_fnc,      & ! Derivative w/ respect to f_diff of function A_fnc.         []
      lhs_upper,  & ! Contribution of upper threshold to implicit portion (LHS). [s^-1]
      lhs_lower     ! Contribution of lower threshold to implicit portion (LHS). [s^-1]


    ! Compute the clipping time scale, dt_clip.
    dt_clip = dt_clip_coef * dt


    ! Upper Threshold
    if ( l_upper_thresh ) then

      ! f_diff is the difference between the threshold value and f_unclipped.
      ! In regards to the upper threshold, it is defined as
      ! upper_threshold - f_unclipped.
      f_diff = upper_threshold - f_unclipped

      ! Compute the values of functions A_fnc and B_fnc evaluated at f_diff(t)
      ! for the upper threshold.
      call compute_fncts_A_B( l_upper_thresh, upper_threshold,  & 
                              l_lower_thresh, lower_threshold,  & 
                              f_diff, A_fnc, B_fnc )

      ! Compute the implicit (LHS) contribution from clipping for the upper
      ! threshold.
      lhs_upper = compute_clip_lhs( dt_clip, B_fnc )

    else

      lhs_upper = 0.0_core_rknd

    endif


    ! Lower Threshold
    if ( l_lower_thresh ) then

      ! f_diff is the difference between the threshold value and f_unclipped.
      ! In regards to the lower threshold, it is defined as
      ! f_unclipped - lower_threshold.
      f_diff = f_unclipped - lower_threshold

      ! Compute the values of functions A_fnc and B_fnc evaluated at f_diff(t)
      ! for the lower threshold.
      call compute_fncts_A_B( l_upper_thresh, upper_threshold,  & 
                              l_lower_thresh, lower_threshold,  & 
                              f_diff, A_fnc, B_fnc )

      ! Compute the implicit (LHS) contribution from clipping for the lower
      ! threshold.
      lhs_lower = compute_clip_lhs( dt_clip, B_fnc )

    else

      lhs_lower = 0.0_core_rknd

    endif


    ! Total implicit (LHS) contribution to clipping.
    ! Main diagonal: [ x f_unclipped(k,<t+1>) ]
    lhs = lhs_upper + lhs_lower


  end function clip_semi_imp_lhs

  !=============================================================================
  pure function compute_clip_lhs( dt_clip, B_fnc ) & 
  result( lhs_contribution )

    ! Description:
    ! Calculation of the implicit portion of the semi-implicit clipping term.
    !
    ! The implicit portion of the semi-implicit clipping term is:
    !
    ! - (1/dt_clip) * B_fnc * f_unclipped(t+1).
    !
    ! Note:  When the term is brought over to the left-hand side, the sign
    !        is reversed and the leading "-" in front of the term is changed
    !        to a "+".
    !
    ! The timestep index (t+1) means that the value of f being used is from the
    ! next timestep, which is being advanced to in solving the d(f)/dt equation.

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd), intent(in) ::  & 
      dt_clip ! Time scale for semi-implicit clipping term.        [s]

    real( kind = core_rknd ), intent(in) :: & 
      B_fnc   ! Derivative w/ respect to f_diff of function A_fnc. []

    ! Return Variable
    real( kind = core_rknd ) :: lhs_contribution


    ! Main diagonal: [ x f_unclipped(k,<t+1>) ]
    lhs_contribution & 
    = + (1.0_core_rknd/dt_clip * B_fnc )


  end function compute_clip_lhs

  !=============================================================================
  function clip_semi_imp_rhs( dt, f_unclipped,  & 
                             l_upper_thresh, upper_threshold,  & 
                             l_lower_thresh, lower_threshold ) & 
  result( rhs )

    ! Description:
    ! The explicit portion of the semi-implicit clipping code.
    !
    ! The explicit (RHS) portion of the equation for clipping with the upper
    ! threshold is:
    !
    ! + (1/dt_clip)
    !   * { A_fnc - B_fnc * f_diff(t) + B_fnc * upper_threshold }.
    !
    ! The explicit (RHS) portion of the equation for clipping with the lower
    ! threshold is:
    !
    ! - (1/dt_clip)
    !   * { A_fnc - B_fnc * f_diff(t) - B_fnc * lower_threshold }.
    !
    ! Timestep index (t) stands for the index of the current timestep.
    !
    ! The values of A_fnc, B_fnc, and f_diff will differ between the two
    ! thresholds.
    !
    ! The overall explicit (RHS) portion for the clipping term is the sum of the
    ! explicit portion from the upper threshold and the explicit portion from
    ! the lower threshold.

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      dt                 ! Model timestep.                                    [s]

    real( kind = core_rknd ), intent(in) :: & 
      f_unclipped,        & ! The unclipped value of variable f at timestep (t). [f units]
      upper_threshold,    & ! Greatest allowable value of variable f.            [f units]
      lower_threshold       ! Smallest allowable value of variable f.            [f units]

    logical, intent(in) :: & 
      l_upper_thresh,    & ! Flag for having an upper threshold value.
      l_lower_thresh       ! Flag for having a lower threshold value.

    ! Return Variable
    real( kind = core_rknd ) :: rhs

    ! Local Variables
    real( kind = core_rknd) ::  & 
      dt_clip       ! Time scale for semi-implicit clipping term.                [s]

    real( kind = core_rknd ) :: & 
      f_diff,     & ! Difference between the threshold value and f_unclipped.    [f units]
      A_fnc,      & ! Function that approximates { f_diff * [ 1 - H(f_diff) ] }. [f units]
      B_fnc,      & ! Derivative w/ respect to f_diff of function A_fnc.         []
      rhs_upper,  & ! Contribution of upper threshold to explicit portion (RHS). [s^-1]
      rhs_lower     ! Contribution of lower threshold to explicit portion (RHS). [s^-1]


    ! Compute the clipping time scale, dt_clip.
    dt_clip = dt_clip_coef * dt


    ! Upper Threshold
    if ( l_upper_thresh ) then

      ! f_diff is the difference between the threshold value and f_unclipped.
      ! In regards to the upper threshold, it is defined as
      ! upper_threshold - f_unclipped.
      f_diff = upper_threshold - f_unclipped

      ! Compute the values of functions A_fnc and B_fnc evaluated at f_diff(t)
      ! for the upper threshold.
      call compute_fncts_A_B( l_upper_thresh, upper_threshold,  & 
                              l_lower_thresh, lower_threshold,  & 
                              f_diff, A_fnc, B_fnc )

      ! Compute the explicit (RHS) contribution from clipping for the upper
      ! threshold.
      rhs_upper  & 
      = + (1.0_core_rknd/dt_clip  & 
          * ( A_fnc - B_fnc * f_diff + B_fnc * upper_threshold ) )

    else

      rhs_upper = 0.0_core_rknd

    endif


    ! Lower Threshold
    if ( l_lower_thresh ) then

      ! f_diff is the difference between the threshold value and f_unclipped.
      ! In regards to the lower threshold, it is defined as
      ! f_unclipped - lower_threshold.
      f_diff = f_unclipped - lower_threshold

      ! Compute the values of functions A_fnc and B_fnc evaluated at f_diff(t)
      ! for the lower threshold.
      call compute_fncts_A_B( l_upper_thresh, upper_threshold,  & 
                              l_lower_thresh, lower_threshold,  & 
                              f_diff, A_fnc, B_fnc )

      ! Compute the explicit (RHS) contribution from clipping for the lower
      ! threshold.
      rhs_lower  & 
      = - (1.0_core_rknd/ dt_clip)  & 
          * ( A_fnc - B_fnc * f_diff - B_fnc * lower_threshold )

    else

      rhs_lower = 0.0_core_rknd

    endif


    ! Total explicit (RHS) contribution to clipping.
    rhs = rhs_upper + rhs_lower


  end function clip_semi_imp_rhs

  !=============================================================================
  subroutine compute_fncts_A_B( l_upper_thresh, upper_threshold,  & 
                                l_lower_thresh, lower_threshold,  & 
                                f_diff, A_fnc, B_fnc )

    ! Description:
    ! This subroutine computes the values of two functions used in semi-implicit
    ! clipping.  Both of the functions are based on the values of f_diff(t) and
    ! the parameter sigma.  One function is A_fnc, which is F_R(f_diff)
    ! evaluated at f_diff = f_diff(t).  F_R(f_diff) is a three-piece function
    ! that is used to approximate function R(f_diff).  The other function is
    ! B_fnc, the derivative with respect to f_diff of function A_fnc.  In other
    ! words, B_fnc is ( d F_R / d f_diff ) evaluated at f_diff = f_diff(t).
    !
    ! The equation for A_fnc is:
    !
    !         | f_diff(t);   where f_diff(t) <= -sigma
    !         |
    ! A_fnc = | f_diff(t) - (sigma/4) * [ 1 + (f_diff(t)/sigma) ]^2;
    !         |    where -sigma < f_diff(t) < sigma
    !         |
    !         | 0;   where f_diff(t) >= sigma;
    !
    ! while the equation for B_fnc is:
    !
    !         | 1;   where f_diff(t) <= -sigma
    !         |
    ! B_fnc = | 1 - (1/2) * [ 1 + (f_diff(t)/sigma) ];
    !         |    where -sigma < f_diff(t) < sigma
    !         |
    !         | 0;   where f_diff(t) >= sigma;
    !
    ! where timestep index (t) stands for the index of the current timestep.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: eps ! Variable(s)

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variable
    real( kind = core_rknd ), intent(in) :: & 
      f_diff,    & ! Difference between the threshold value and f_unclipped.  [f units]
      upper_threshold,    & ! Greatest allowable value of variable f.         [f units]
      lower_threshold       ! Smallest allowable value of variable f.         [f units]

    logical, intent(in) :: & 
      l_upper_thresh,    & ! Flag for having an upper threshold value.
      l_lower_thresh       ! Flag for having a lower threshold value.

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: & 
      A_fnc,   & ! Function that approximates { f_diff * [ 1 - H(f_diff) ] }. [f units]
      B_fnc      ! Derivative w/ respect to f_diff of function A_fnc.         []

    ! Local Variables
    real( kind = core_rknd ) :: sigma_val,       & ! Value of parameter sigma. [f units]
            thresh_avg_mag     ! Average magnitude of threshold(s).         [f units]

    thresh_avg_mag = 0.0_core_rknd ! Default Initialization

    ! Find the average magnitude of the threshold.
    ! In cases where only one threshold applies, the average magnitude of the
    ! threshold must be greater than 0.
    ! Note:  The constant eps is there in case only one threshold applies, and
    !        it has a value of 0 (or very close to 0).  However, eps is a very
    !        small number, and therefore it will not start curbing values until
    !        they get extremely close to the threshold.  A larger constant value
    !        may work better.
    if ( l_upper_thresh .and. l_lower_thresh ) then
      ! Both thresholds apply.
      thresh_avg_mag = 0.5_core_rknd * (   abs(upper_threshold)  & 
                               + abs(lower_threshold)   )
    elseif ( l_upper_thresh ) then
      ! Only the upper threshold applies.
      thresh_avg_mag = max( abs(upper_threshold), eps )
    elseif ( l_lower_thresh ) then
      ! Only the lower threshold applies.
      thresh_avg_mag = max( abs(lower_threshold), eps )
    endif

    ! Compute the value of sigma based on the magnitude of the threshold(s) for
    ! variable f and the sigma coefficient.  The value of sigma must always be
    ! positive.
    sigma_val = sigma_coef * thresh_avg_mag

    ! A_fnc is a three-piece function that approximates function
    ! R(f_diff(t)) = { f_diff(t) * [ 1 - H(f_diff(t)) ] }.  This is needed
    ! because the R(f_diff(t)) is not differentiable at point f_diff(t) = 0, as
    ! the function has a corner at that point.  Function A_fnc is differentiable
    ! at all points.  It is evaluated for f_diff at timestep index (t).
    if ( f_diff <= -sigma_val ) then
      A_fnc = f_diff
    elseif ( f_diff >= sigma_val ) then
      A_fnc = 0.0_core_rknd
    else ! -sigma_val < f_diff < sigma_val
      A_fnc = f_diff - ( (sigma_val/4.0_core_rknd)  & 
                         * ( 1.0_core_rknd + f_diff/sigma_val )**2 )
    endif

    ! B_fnc is the derivative with respect to f_diff of function A_fnc.  It is
    ! evaluated for f_diff at timestep index (t).
    if ( f_diff <= -sigma_val ) then
      B_fnc = 1.0_core_rknd
    elseif ( f_diff >= sigma_val ) then
      B_fnc = 0.0_core_rknd
    else ! -sigma_val < f_diff < sigma_val
      B_fnc = 1.0_core_rknd - (1.0_core_rknd/2.0_core_rknd)*( 1.0_core_rknd + f_diff/sigma_val )
    endif


  end subroutine compute_fncts_A_B

!===============================================================================

end module clip_semi_implicit
