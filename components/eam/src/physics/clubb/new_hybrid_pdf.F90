! $Id$
!===============================================================================
module new_hybrid_pdf

  ! Description:
  ! The portion of CLUBB's multivariate, two-component PDF that is the
  ! multivariate, two-component normal PDF of vertical velocity (w), total water
  ! mixing ratio (rt), liquid water potential temperature (thl), and optionally,
  ! the west-east horizontal wind component (u), the south-north horizontal wind
  ! component (v), and passive scalars (sclr).

  ! References:
  ! Griffin and Larson (2020)
  !-------------------------------------------------------------------------

  implicit none

  public :: calculate_mixture_fraction,  & ! Procedure(s)
            calculate_w_params,          &
            calculate_responder_params,  &
            calculate_coef_wp4_implicit, &
            calc_coef_wp2xp_implicit,    &
            calc_coefs_wpxp2_semiimpl,   &
            calc_coefs_wpxpyp_semiimpl

  private

  contains

  !=============================================================================
  !
  ! DESCRIPTION OF THE METHOD FOR THE VARIABLE THAT SETS THE MIXTURE FRACTION
  ! =========================================================================
  !
  ! The variable that sets the mixture fraction for the PDF is w.  There are
  ! five PDF parameters that need to be calculated, which are mu_w_1 (the mean
  ! of w is the 1st PDF component), mu_w_2 (the mean of w in the 2nd PDF
  ! component), sigma_w_1 (the standard deviation of w in the 1st PDF
  ! component), sigma_w_2 (the standard deviation of w in the 2nd PDF
  ! component), and mixt_frac (the mixture fraction, which is the weight of the
  ! 1st PDF component).  In order to solve for these five parameters, five
  ! equations are needed.  These five equations are the equations for <w>,
  ! <w'^2>, and <w'^3> as found by integrating over the PDF.  Additionally, two
  ! more equations, which involve tunable parameters F_w and zeta_w, and which
  ! are used to help control the spread of the PDF component means and the size
  ! of the PDF component standard deviations compared to each other,
  ! respectively, are used in this equation set.  The five equations are:
  !
  ! <w> = mixt_frac * mu_w_1 + ( 1 - mixt_frac ) * mu_w_2;
  !
  ! <w'^2> = mixt_frac * ( ( mu_w_1 - <w> )^2 + sigma_w_1^2 )
  !          + ( 1 - mixt_frac ) * ( ( mu_w_2 - <w> )^2 + sigma_w_2^2 );
  !
  ! <w'^3> = mixt_frac * ( mu_w_1 - <w> )
  !                    * ( ( mu_w_1 - <w> )^2 + 3 * sigma_w_1^2 )
  !          + ( 1 - mixt_frac ) * ( mu_w_2 - <w> )
  !                              * ( ( mu_w_2 - <w> )^2 + 3 * sigma_w_2^2 );
  !
  ! mu_w_1 - <w> = sqrt(F_w) * ( sqrt( 1 - mixt_frac ) / sqrt( mixt_frac ) )
  !                * sqrt( <w'^2> );
  !
  ! where 0 <= F_w <= 1; and
  !
  ! 1 + zeta_w = ( mixt_frac * sigma_w_1^2 )
  !              / ( ( 1 - mixt_frac ) * sigma_w_2^2 );
  !
  ! where zeta_w > -1.
  !
  ! Following convention for w, mu_w_1 is defined to be greater than or equal to
  ! mu_w_2 (and is also greater than or equal to <w>, while mu_w_2 is less than
  ! or equal to <w>).  This is relationship is found in the mu_w_1 - <w>
  ! equation above.
  !
  ! The resulting equations for the five PDF parameters are:
  !
  ! mixt_frac
  ! = ( 4 * F_w^3
  !     + 18 * F_w * ( zeta_w + 1 ) * ( 1 - F_w ) / ( zeta_w + 2 )
  !     + 6 * F_w^2 * ( 1 - F_w ) / ( zeta_w + 2 )
  !     + Skw^2
  !     - Skw * sqrt( 4 * F_w^3
  !                   + 12 * F_w^2 * ( 1 - F_w )
  !                   + 36 * F_w * ( zeta_w + 1 ) * ( 1 - F_w )^2
  !                     / ( zeta_w + 2 )^2
  !                   + Skw^2 ) )
  !   / ( 2 * F_w * ( F_w - 3 )^2 + 2 * Skw^2 );
  !
  ! mu_w_1 = <w> + sqrt( F_w * ( ( 1 - mixt_frac ) / mixt_frac ) * <w'^2> );
  !
  ! mu_w_2 = <w> - ( mixt_frac / ( 1 - mixt_frac ) ) * ( mu_w_1 - <w> );
  !
  ! sigma_w_1 = sqrt( ( ( zeta_w + 1 ) * ( 1 - F_w ) )
  !                   / ( ( zeta_w + 2 ) * mixt_frac ) * <w'^2> ); and
  !
  ! sigma_w_2 = sqrt( ( mixt_frac * sigma_w_1^2 )
  !                   / ( ( 1 - mixt_frac ) * ( 1 + zeta_w ) ) );
  !
  ! where Skw is the skewness of w, and Skw = <w'^3> / <w'^2>^(3/2).
  !
  ! This method works for all values of F_w (where 0 <= F_w <= 1) and zeta_w
  ! (where zeta_w > -1).
  !
  !
  ! Special case:
  !
  ! When Skw = 0 and F_w = 0, the equation for mixt_frac is undefined.  The
  ! equation for mixture fraction in this scenario can be derived by using the
  ! above equation for mixture fraction and then setting Skw = 0.  The resulting
  ! equation becomes:
  !
  ! mixt_frac
  ! = ( 4 * F_w^3
  !     + 18 * F_w * ( zeta_w + 1 ) * ( 1 - F_w ) / ( zeta_w + 2 )
  !     + 6 * F_w^2 * ( 1 - F_w ) / ( zeta_w + 2 ) )
  !   / ( 2 * F_w * ( F_w - 3 )^2 ).
  !
  ! All of the terms in the numerator and denominator contain a F_w, so this
  ! equation can be rewritten as:
  !
  ! mixt_frac
  ! = ( 4 * F_w^2
  !     + 18 * ( zeta_w + 1 ) * ( 1 - F_w ) / ( zeta_w + 2 )
  !     + 6 * F_w * ( 1 - F_w ) / ( zeta_w + 2 ) )
  !   / ( 2 * ( F_w - 3 )^2 ).
  !
  ! Now setting F_w = 0, the equation becomes:
  !
  ! mixt_frac = ( 18 * ( zeta_w + 1 ) / ( zeta_w + 2 ) ) / 18;
  !
  ! which can be rewritten as:
  !
  ! mixt_frac = ( zeta_w + 1 ) / ( zeta_w + 2 ).
  !
  ! When F_w = 0, Skw must have a value of 0 in order for the PDF to function
  ! correctly.  When F_w = 0, mu_w_1 = mu_w_2.  When the two PDF component means
  ! are equal to each other (and to the overall mean, <w>), the only value of
  ! Skw that can be represented is a value of 0.  In the equation for mixture
  ! fraction, when F_w is set to 0, but | Skw | > 0, mixt_frac will either have
  ! a value of 0 or 1, depending on whether Skw is positive or negative,
  ! respectively.
  !
  ! The value of F_w should be set as a function of Skw.  The value F_w should
  ! go toward 0 as | Skw | (or Skw^2) goes toward 0.  The value of F_w should
  ! go toward 1 as | Skw | (or Skw^2) goes to infinity.
  !
  !
  ! Tunable parameters:
  !
  ! 1) F_w:  This parameter controls the spread of the PDF component means.  The
  !          range of this parameter is 0 <= F_w <= 1.  When F_w = 0, the two
  !          PDF component means (mu_w_1 and mu_w_2) are equal to each other
  !          (and Skw must equal 0).  All of the variance of w is accounted for
  !          by the PDF component standard deviations (sigma_w_1 and sigma_w_2).
  !          When F_w = 1, mu_w_1 and mu_w_2 are spread as far apart as they can
  !          be.  Both PDF component standard deviations (sigma_w_1 and
  !          sigma_w_2) are equal to 0, and all of the variance of w is
  !          accounted for by the spread of the PDF component means.
  !
  !          When sigma_w_1 = sigma_w_2 = 0, the equation for <w'^2> becomes:
  !
  !          <w'^2> = mixt_frac * ( mu_w_1 - <w> )^2
  !                   + ( 1 - mixt_frac ) * ( mu_w_2 - <w> )^2.
  !
  !          Substituting the equation for <w> into the above equation for
  !          mu_w_2 - <w>, the above equation becomes:
  !
  !          <w'^2> = ( mixt_frac / ( 1 - mixt_frac ) ) * ( mu_w_1 - <w> )^2;
  !
  !          which can be rewritten as:
  !
  !          ( mu_w_1 - <w> )^2 = ( ( 1 - mixt_frac ) / mixt_frac ) * <w'^2>.
  !
  !          Taking the square root of the above equation:
  !
  !          mu_w_1 - <w> = +/- ( sqrt( 1 - mixt_frac ) / sqrt(mixt_frac) )
  !                             * sqrt( <w'^2> ).
  !
  !          This equation can be compared to the equation for mu_w_1 - <w> in
  !          the set of 5 equations, which is:
  !
  !          mu_w_1 - <w>
  !          = sqrt(F_w) * ( sqrt( 1 - mixt_frac ) / sqrt( mixt_frac ) )
  !            * sqrt( <w'^2> ).
  !
  !          The above equations give another example of the meaning of F_w.
  !          The value of sqrt(F_w) is ratio of mu_w_1 - <w> to its maximum
  !          value, which is:
  !
  !          sqrt( ( ( 1 - mixt_frac ) / mixt_frac ) * <w'^2> ).
  !
  !
  ! 2) zeta_w:  This parameter controls the size of the PDF component standard
  !             deviations compared to each other.  The equation for zeta_w is:
  !
  !             1 + zeta_w = ( mixt_frac * sigma_w_1^2 )
  !                          / ( ( 1 - mixt_frac ) * sigma_w_2^2 ).
  !
  !             When zeta_w > 0, mixt_frac * sigma_w_1^2 increases at the
  !             expense of ( 1 - mixt_frac ) * sigma_w_2^2, which decreases in
  !             this variance-preserving equation set.  When zeta_w = 0, then
  !             mixt_frac * sigma_w_1^2 = ( 1 - mixt_frac ) * sigma_w_2^2.
  !             When -1 < zeta_w < 0, ( 1 - mixt_frac ) * sigma_w_2^2 increases
  !             at the expense of mixt_frac * sigma_w_1^2, which decreases.  As
  !             a result, greater values of zeta_w cause the 1st PDF component
  !             to become broader while the 2nd PDF component becomes narrower,
  !             and smaller values of zeta_w cause the 1st PDF component to
  !             become narrower while the 2nd PDF component becomes broader.
  !
  !             Symmetry
  !
  !             When zeta_w = 0, the PDF is always symmetric.  In other words,
  !             the PDF at any positive value of Skw (for example, Skw = 2.5)
  !             will look like a mirror-image (reflection across the y-axis)
  !             of the PDF at a negative value of Skw of the same magnitude (in
  !             this example, Skw = -2.5).  However, when zeta_w /= 0, the PDF
  !             loses this quality and is not symmetric.
  !
  !             When symmetry is desired at values of zeta_w besides zeta_w = 0,
  !             the solution is to turn zeta_w into a function of Skw.  A basic
  !             example of a zeta_w skewness equation that produces a symmetric
  !             PDF for values of zeta_w other than 0 is:
  !
  !             zeta_w = | zeta_w_in,                      when Skw >= 0;
  !                      | ( 1 / ( 1 + zeta_w_in ) ) - 1,  when Skw < 0.
  !
  !
  ! Notes:
  !
  ! When F_w = 0 (which can only happen when Skw = 0), mu_w_1 = mu_w_2, and
  ! mixt_frac = ( zeta_w + 1 ) / ( zeta_w + 2 ).  When these equations are
  ! substituted into the equations for sigma_w_1 and sigma_w_2, the result is
  ! sigma_w_1 = sigma_w_2 = sqrt( <w'^2> ).  This means that the distribution
  ! becomes a single Gaussian when F_w = 0 (and Skw = 0).  This happens
  ! regardless of the value of zeta_w.
  !
  ! The equations for the PDF component means and standard deviations can also
  ! be written as:
  !
  ! mu_w_1 = <w> + sqrt( F_w * ( ( 1 - mixt_frac ) / mixt_frac ) * <w'^2> );
  !
  ! mu_w_2 = <w> - sqrt( F_w * ( mixt_frac / ( 1 - mixt_frac ) ) * <w'^2> );
  !
  ! sigma_w_1 = sqrt( coef_sigma_w_1_sqd * <w'^2> ); and
  !
  ! sigma_w_2 = sqrt( coef_sigma_w_2_sqd * <w'^2> ); where
  !
  ! coef_sigma_w_1_sqd = ( ( zeta_w + 1 ) * ( 1 - F_w ) )
  !                      / ( ( zeta_w + 2 ) * mixt_frac ); and
  !
  ! coef_sigma_w_2_sqd = ( 1 - F_w ) / ( ( zeta_w + 2 ) * ( 1 - mixt_frac ) ).
  !
  ! The above equations can be substituted into an equation for a variable that
  ! has been derived by integrating over the PDF.  Many variables like this are
  ! used in parts of the predictive equation set.  These substitutions allow
  ! some terms to solved implicitly or semi-implicitly in the predictive
  ! equations.
  !
  !
  ! Brian Griffin; September 2017.
  !
  !=============================================================================
  elemental function calculate_mixture_fraction( Skw, F_w, zeta_w ) &
  result( mixt_frac )

    ! Description:
    ! Calculates mixture fraction.

    ! References:
    ! Griffin and Larson (2020)
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        thirty_six, & ! Constant(s)
        eighteen,   &
        twelve,     &
        six,        &
        four,       &
        three,      &
        two,        &
        one,        &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Skw,    & ! Skewness of w                                            [-]
      F_w,    & ! Parameter for the spread of the PDF component means of w [-]
      zeta_w    ! Parameter for the PDF component variances of w           [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      mixt_frac    ! Mixture fraction    [-]


    ! Calculate mixture fraction, which is the weight of the 1st PDF component.
    ! The 2nd PDF component has a weight of 1 - mixt_frac.
    if ( F_w > zero ) then

       mixt_frac &
       = ( four * F_w**3 &
           + eighteen * F_w &
             * ( zeta_w + one ) * ( one - F_w ) / ( zeta_w + two ) &
           + six * F_w**2 * ( one - F_w ) / ( zeta_w + two ) &
           + Skw**2 &
           - Skw * sqrt( four * F_w**3 &
                         + twelve * F_w**2 * ( one - F_w ) &
                         + thirty_six * F_w &
                           * ( zeta_w + one ) * ( one - F_w )**2 &
                           / ( zeta_w + two )**2 &
                         + Skw**2 ) ) &
         / ( two * F_w * ( F_w - three )**2 + two * Skw**2 )

    else ! F_w = 0

       if ( abs( Skw ) > zero ) then

          ! When F_w = 0, | Skw | must have a value of 0.  In a scenario where
          ! F_w = 0 and | Skw | > 0, the mixture fraction (and the rest of the
          ! PDF parameters) can't be calculated.  Since mixture fraction can
          ! only have values 0 < mixt_frac < 1, set mixt_frac to -1 in this
          ! scenario.
          mixt_frac = -one

       else ! Skw = 0

          mixt_frac = ( zeta_w + one ) / ( zeta_w + two )

       endif ! | Skw | > 0

    endif ! F_w > 0


    return

  end function calculate_mixture_fraction

  !=============================================================================
  elemental subroutine calculate_w_params( wm, wp2, Skw, F_w, zeta_w, & ! In
                                           mu_w_1, mu_w_2, sigma_w_1, & ! Out
                                           sigma_w_2, mixt_frac,      & ! Out
                                           coef_sigma_w_1_sqd,        & ! Out
                                           coef_sigma_w_2_sqd         ) ! Out

    ! Description:
    ! Calculates the PDF component means, the PDF component standard deviations,
    ! and the mixture fraction for the variable that sets the PDF.

    ! References:
    ! Griffin and Larson (2020)
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two,  & ! Variable(s)
        one,  &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      wm,     & ! Mean of w (overall)                                     [m/s]
      wp2,    & ! Variance of w (overall)                             [m^2/s^2]
      Skw,    & ! Skewness of w                                             [-]
      F_w,    & ! Parameter for the spread of the PDF component means of w  [-]
      zeta_w    ! Parameter for the PDF component variances of w            [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      mu_w_1,    & ! Mean of w (1st PDF component)                  [m/s]
      mu_w_2,    & ! Mean of w (2nd PDF component)                  [m/s]
      sigma_w_1, & ! Standard deviation of w (1st PDF component)    [m/s]
      sigma_w_2, & ! Standard deviation of w (2nd PDF component)    [m/s]
      mixt_frac    ! Mixture fraction                               [-]

    real( kind = core_rknd ), intent(out) :: &
      coef_sigma_w_1_sqd, & ! sigma_w_1^2 = coef_sigma_w_1_sqd * <w'^2>     [-]
      coef_sigma_w_2_sqd    ! sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>     [-]


    ! Calculate the mixture fraction.
    mixt_frac = calculate_mixture_fraction( Skw, F_w, zeta_w )

    if ( mixt_frac > zero .and. mixt_frac < one ) then

       ! Calculate the mean of w in the 1st PDF component.
       mu_w_1 = wm + sqrt( F_w * ( ( one - mixt_frac ) / mixt_frac ) * wp2 )

       ! Calculate the mean of w in the 2nd PDF component.
       mu_w_2 = wm - ( mixt_frac / ( one - mixt_frac ) ) * ( mu_w_1 - wm )

       ! Calculate the standard deviation of w in the 1st PDF component.
       ! sigma_w_1 = sqrt( ( ( zeta_w + 1 ) * ( 1 - F_w ) )
       !                   / ( ( zeta_w + 2 ) * mixt_frac ) * <w'^2> )
       coef_sigma_w_1_sqd = ( ( zeta_w + one ) * ( one - F_w ) ) &
                            / ( ( zeta_w + two ) * mixt_frac )

       sigma_w_1 = sqrt( coef_sigma_w_1_sqd * wp2 )

       ! Calculate the standard deviation of w in the 2nd PDF component.
       ! sigma_w_2 = sqrt( ( mixt_frac * sigma_w_1^2 )
       !                   / ( ( 1 - mixt_frac ) * ( 1 + zeta_w ) ) )
       !           = sqrt( ( 1 - F_w )
       !                   / ( ( zeta_w + 2 ) * ( 1 - mixt_frac ) ) * <w'^2> )
       coef_sigma_w_2_sqd = ( one - F_w ) &
                            / ( ( zeta_w + two ) * ( one - mixt_frac ) )

       sigma_w_2 = sqrt( coef_sigma_w_2_sqd * wp2 )

    else ! mixt_frac <= 0 or mixt_frac >= 1

       ! The mixture fraction produced is invalid.  This should only happen in
       ! the scenario where F_w = 0 and | Skw | > 0, where the value of
       ! mixt_frac has been set to -1.  Set all output variables to 0 in this
       ! scenario.  Since F_w is a function of skewness, the mixture fraction
       ! and the PDF should always be valid, and this section of code shouldn't
       ! be entered.
       mu_w_1 = zero
       mu_w_2 = zero
       sigma_w_1 = zero
       sigma_w_2 = zero
       coef_sigma_w_1_sqd = zero
       coef_sigma_w_2_sqd = zero

    endif ! 0 < mixt_frac < 1


    return

  end subroutine calculate_w_params

  !=============================================================================
  !
  ! DESCRIPTION OF THE METHOD FOR EACH RESPONDING VARIABLE
  ! ======================================================
  !
  ! In order to find equations for the four PDF parameters for each responding
  ! variable, which are mu_x_1, mu_x_2, sigma_x_1, and sigma_x_2 (where x stands
  ! for a responding variable here), four equations are needed.  These four
  ! equations are the equations for <x>, <x'^2>, <x'^3>, and <w'x'> as found by
  ! integrating over the PDF.  The four equations are:
  !
  ! <x> = mixt_frac * mu_x_1 + ( 1 - mixt_frac ) * mu_x_2;
  !
  ! <x'^2> = mixt_frac * ( ( mu_x_1 - <x> )^2 + sigma_x_1^2 )
  !          + ( 1 - mixt_frac ) * ( ( mu_x_2 - <x> )^2 + sigma_x_2^2 );
  !
  ! <x'^3> = mixt_frac * ( mu_x_1 - <x> )
  !                    * ( ( mu_x_1 - <x> )^2 + 3 * sigma_x_1^2 )
  !          + ( 1 - mixt_frac ) * ( mu_x_2 - <x> )
  !                              * ( ( mu_x_2 - <x> )^2 + 3 * sigma_x_2^2 ); and
  !
  ! <w'x'> = mixt_frac * ( mu_w_1 - <w> ) * ( mu_x_1 - <x> )
  !          + ( 1 - mixt_frac ) * ( mu_w_2 - <w> ) * ( mu_x_2 - <x> );
  !
  ! where the correlations that are normally found in the <w'x'> equation,
  ! corr_w_x_1 and corr_w_x_2, have both been set to 0.
  !
  ! The equations for mu_w_1 - <w> and mu_w_2 - <w> are:
  !
  ! mu_w_1 - <w> = sqrt( F_w * ( ( 1 - mixt_frac ) / mixt_frac ) * <w'^2> );
  !
  ! mu_w_2 - <w> = - sqrt( F_w * ( mixt_frac / ( 1 - mixt_frac ) ) * <w'^2> );
  !
  ! The resulting equations for the four PDF parameters are:
  !
  ! mu_x_1 = <x> + sqrt( ( 1 - mixt_frac ) / mixt_frac )
  !                * <w'x'> / sqrt( F_w * <w'^2> );
  !
  ! mu_x_2 = <x> - sqrt( mixt_frac / ( 1 - mixt_frac ) )
  !                * <w'x'> / sqrt( F_w * <w'^2> );
  !
  ! sigma_x_1^2 = ( 1 + sqrt( ( 1 - mixt_frac ) / mixt_frac )
  !                     * Skx * sqrt( F_w * <w'^2> * <x'^2> ) / ( 3 * <w'x'> )
  !                 - ( ( 1 + mixt_frac ) / mixt_frac )
  !                   * <w'x'>^2 / ( 3 * F_w * <w'^2> * <x'^2> ) )
  !               * <x'^2>; and
  !
  ! sigma_x_2^2 = ( 1 - sqrt( mixt_frac / ( 1 - mixt_frac ) )
  !                     * Skx * sqrt( F_w * <w'^2> * <x'^2> ) / ( 3 * <w'x'> )
  !                 + ( ( mixt_frac - 2 ) / ( 1 - mixt_frac ) )
  !                   * <w'x'>^2 / ( 3 * F_w * <w'^2> * <x'^2> ) )
  !               * <x'^2>;
  !
  ! where Skx is the skewness of x, and Skx = <x'^3> / <x'^2>^(3/2).
  !
  !
  ! Limits on F_w:
  !
  ! The only limits placed on the value of F_w from the w equation set itself
  ! are 0 <= F_w <= 1.  However, use of the above equation set for responder
  ! variable x forces an additional limit to be placed on the value of F_w.
  ! That additional limit restricts the range of F_w to:
  !
  ! <w'x'>^2 / ( <w'^2> * <x'^2> ) <= F_w <= 1.
  !
  ! Furthermore, when there is more than one responder variable, F_w is limited
  ! by the most restrictive cases, such that:
  !
  ! max( <w'x'>^2 / ( <w'^2> * <x'^2> ), for all variables x ) <= F_w <= 1.
  !
  !
  ! Limits on Skx:
  !
  ! Since the PDF parameters for this variable need to work with the mixture
  ! fraction that has been provided by the setting variable, which is w, the
  ! method does not work for all values of Skx.  However, the limits of Skx can
  ! always be calculated.  The limits on Skw given by:
  !
  ! when <w'x'> > 0:
  !
  ! ( 1 + mixt_frac ) / sqrt( mixt_frac * ( 1 - mixt_frac ) )
  ! * <w'x'>^3 / ( F_w * <w'^2> * <x'^2> )^(3/2)
  ! - sqrt( mixt_frac / ( 1 - mixt_frac ) )
  !   * 3 * <w'x'> / sqrt( F_w * <w'^2> * <x'^2> )
  ! <= Skx <=
  ! ( mixt_frac - 2 ) / sqrt( mixt_frac * ( 1 - mixt_frac ) )
  ! * <w'x'>^3 / ( F_w * <w'^2> * <x'^2> )^(3/2)
  ! + sqrt( ( 1 - mixt_frac ) / mixt_frac )
  !   * 3 * <w'x'> / sqrt( F_w * <w'^2> * <x'^2> );
  !
  ! when <w'x'> < 0:
  !
  ! ( mixt_frac - 2 ) / sqrt( mixt_frac * ( 1 - mixt_frac ) )
  ! * <w'x'>^3 / ( F_w * <w'^2> * <x'^2> )^(3/2)
  ! + sqrt( ( 1 - mixt_frac ) / mixt_frac )
  !   * 3 * <w'x'> / sqrt( F_w * <w'^2> * <x'^2> )
  ! <= Skx <=
  ! ( 1 + mixt_frac ) / sqrt( mixt_frac * ( 1 - mixt_frac ) )
  ! * <w'x'>^3 / ( F_w * <w'^2> * <x'^2> )^(3/2)
  ! - sqrt( mixt_frac / ( 1 - mixt_frac ) )
  !   * 3 * <w'x'> / sqrt( F_w * <w'^2> * <x'^2> );
  !
  ! and when <w'x'> = 0, Skx = 0.
  !
  !
  ! Special cases:
  !
  ! When <w'x'> = 0, mu_x_1 = mu_x_2 = <x>, and the value of Skx must be 0.
  ! Since both <w'x'> = 0 and Skx = 0, the equations for sigma_x_1^2 and
  ! sigma_x_2^2 are both undefined.  In this situation, the equations for the
  ! PDF parameters of x are:
  !
  ! mu_x_1 = mu_x_2 = <x>; and
  ! sigma_x_1^2 = sigma_x_2^2 = <x'^2>.
  !
  ! The value of F_w is allowed to be 0 only when <w'x'> = 0 (for all variables
  ! x).  When <w'^2> = 0 and/or <x'^2> = 0, this means that <w'x'> = 0, as well.
  ! In all these situations, the equation set for the situation when <w'x'> = 0
  ! is used.  This means that the distribution becomes a single Gaussian when
  ! <w'x'> = 0 (and Skx = 0).
  !
  !
  ! Notes:
  !
  ! The equations for the PDF component means and standard deviations can also
  ! be written as:
  !
  ! mu_x_1 = <x> + sqrt( ( 1 - mixt_frac ) / mixt_frac )
  !                * <w'x'> / sqrt( F_w * <w'^2> );
  !
  ! mu_x_2 = <x> - sqrt( mixt_frac / ( 1 - mixt_frac ) )
  !                * <w'x'> / sqrt( F_w * <w'^2> );
  !
  ! sigma_x_1 = sqrt( coef_sigma_x_1_sqd * <x'^2> ); and
  !
  ! sigma_x_2 = sqrt( coef_sigma_x_2_sqd * <x'^2> ); where
  !
  ! coef_sigma_x_1_sqd
  ! = 1 + sqrt( ( 1 - mixt_frac ) / mixt_frac )
  !       * Skx * sqrt( F_w * <w'^2> * <x'^2> ) / ( 3 * <w'x'> )
  !   - ( ( 1 + mixt_frac ) / mixt_frac )
  !     * <w'x'>^2 / ( 3 * F_w * <w'^2> * <x'^2> ); and
  !
  ! coef_sigma_x_2_sqd
  ! = 1 - sqrt( mixt_frac / ( 1 - mixt_frac ) )
  !       * Skx * sqrt( F_w * <w'^2> * <x'^2> ) / ( 3 * <w'x'> )
  !   + ( ( mixt_frac - 2 ) / ( 1 - mixt_frac ) )
  !     * <w'x'>^2 / ( 3 * F_w * <w'^2> * <x'^2> ).
  !
  ! The above equations can be substituted into an equation for a variable that
  ! has been derived by integrating over the PDF.  Many variables like this are
  ! used in parts of the predictive equation set.  These substitutions allow
  ! some terms to solved implicitly or semi-implicitly in the predictive
  ! equations.
  !
  !
  ! Brian Griffin; September 2019.
  !
  !=============================================================================
  elemental subroutine calculate_responder_params( xm, xp2, Skx, wpxp,  & ! In
                                                   wp2, F_w, mixt_frac, & ! In
                                                   mu_x_1, mu_x_2,      & ! Out
                                                   sigma_x_1_sqd,       & ! Out
                                                   sigma_x_2_sqd,       & ! Out
                                                   coef_sigma_x_1_sqd,  & ! Out
                                                   coef_sigma_x_2_sqd   ) ! Out

    ! Description:
    ! Calculates the PDF component means and the PDF component standard
    ! deviations for a responding variable (a variable that is not used to set
    ! the mixture fraction).

    ! References:
    ! Griffin and Larson (2020)
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        three, & ! Variable(s)
        two,   &
        one,   &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      xm,       & ! Mean of x (overall)                             [units vary]
      xp2,      & ! Variance of x (overall)                     [(units vary)^2]
      Skx,      & ! Skewness of x                                            [-]
      wpxp,     & ! Covariance of w and x (overall)            [m/s(units vary)]
      wp2,      & ! Variance of w (overall)                            [m^2/s^2]
      F_w,      & ! Parameter for the spread of the PDF component means of w [-]
      mixt_frac   ! Mixture fraction                                         [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      mu_x_1,        & ! Mean of x (1st PDF component)              [units vary]
      mu_x_2,        & ! Mean of x (2nd PDF component)              [units vary]
      sigma_x_1_sqd, & ! Variance of x (1st PDF component)      [(units vary)^2]
      sigma_x_2_sqd    ! Variance of x (2nd PDF component)      [(units vary)^2]

    real( kind = core_rknd ), intent(out) :: &
      coef_sigma_x_1_sqd, & ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>      [-]
      coef_sigma_x_2_sqd    ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>      [-]


    if ( abs( wpxp ) > zero ) then

       ! Note:  when |<w'x'>| > 0, F_w, <w'^2>, and <x'^2> must all have values
       !        greater than 0.

       ! Calculate the mean of x in the 1st PDF component.
       mu_x_1 = xm + sqrt( ( one - mixt_frac ) / mixt_frac ) &
                     * wpxp / sqrt( F_w * wp2 )

       ! Calculate the mean of x in the 2nd PDF component.
       mu_x_2 = xm - sqrt( mixt_frac / ( one - mixt_frac ) ) &
                     * wpxp / sqrt( F_w * wp2 )

       ! Calculate the variance of x in the 1st PDF component.
       ! sigma_x_1^2
       ! = ( 1 + sqrt( ( 1 - mixt_frac ) / mixt_frac )
       !         * Skx * sqrt( F_w * <w'^2> * <x'^2> ) / ( 3 * <w'x'> )
       !     - ( ( 1 + mixt_frac ) / mixt_frac )
       !       * <w'x'>^2 / ( 3 * F_w * <w'^2> * <x'^2> ) )
       !   * <x'^2>
       coef_sigma_x_1_sqd &
       = one + sqrt( ( one - mixt_frac ) / mixt_frac ) &
               * Skx * sqrt( F_w * wp2 * xp2 ) / ( three * wpxp ) &
         - ( ( one + mixt_frac ) / mixt_frac ) &
           * wpxp**2 / ( three * F_w * wp2 * xp2 )

       ! Mathematically, the value of coef_sigma_x_1_sqd cannot be less than 0.
       ! Numerically, this can happen when numerical round off error causes an
       ! epsilon-sized negative value.  When this happens, reset the value of
       ! coef_sigma_x_1_sqd to 0.
       coef_sigma_x_1_sqd = max( coef_sigma_x_1_sqd, zero )

       sigma_x_1_sqd = coef_sigma_x_1_sqd * xp2

       ! Calculate the variance of x in the 2nd PDF component.
       ! sigma_x_2^2
       ! = ( 1 - sqrt( mixt_frac / ( 1 - mixt_frac ) )
       !         * Skx * sqrt( F_w * <w'^2> * <x'^2> ) / ( 3 * <w'x'> )
       !     + ( ( mixt_frac - 2 ) / ( 1 - mixt_frac ) )
       !       * <w'x'>^2 / ( 3 * F_w * <w'^2> * <x'^2> ) )
       !   * <x'^2>
       coef_sigma_x_2_sqd &
       = one - sqrt( mixt_frac / ( one - mixt_frac ) ) &
               * Skx * sqrt( F_w * wp2 * xp2 ) / ( three * wpxp ) &
         + ( ( mixt_frac - two ) / ( one - mixt_frac ) ) &
           * wpxp**2 / ( three * F_w * wp2 * xp2 )

       ! Mathematically, the value of coef_sigma_x_2_sqd cannot be less than 0.
       ! Numerically, this can happen when numerical round off error causes an
       ! epsilon-sized negative value.  When this happens, reset the value of
       ! coef_sigma_x_2_sqd to 0.
       coef_sigma_x_2_sqd = max( coef_sigma_x_2_sqd, zero )

       sigma_x_2_sqd = coef_sigma_x_2_sqd * xp2

    else ! | <w'x'> | = 0

       ! When <w'x'> has a value of 0, the PDF becomes a single Gaussian.  This
       ! only works when Skx = 0.  However, when Skx /= 0, the value of min_F_x
       ! is greater than 0, preventing a problem where F_x = 0 but | Skx | > 0.
       mu_x_1 = xm
       mu_x_2 = xm
       sigma_x_1_sqd = xp2
       sigma_x_2_sqd = xp2
       coef_sigma_x_1_sqd = one
       coef_sigma_x_2_sqd = one

    endif ! | <w'x'> | > 0


    return

  end subroutine calculate_responder_params

  !=============================================================================
  elemental function calculate_coef_wp4_implicit( mixt_frac, F_w, &
                                                  coef_sigma_w_1_sqd, &
                                                  coef_sigma_w_2_sqd ) &
  result( coef_wp4_implicit )

    ! Description:
    ! The predictive equation for <w'^3> contains a turbulent advection term of
    ! the form:
    !
    ! - ( 1 / rho_ds ) * d ( rho_ds * <w'^4> ) / dz;
    !
    ! where z is height, rho_ds is the dry, base-state density, and <w'^4> is
    ! calculated by integrating over the PDF.  The equation for <w'^4> is:
    !
    ! <w'^4> = mixt_frac * ( 3 * sigma_w_1^4
    !                        + 6 * ( mu_w_1 - <w> )^2 * sigma_w_1^2
    !                        + ( mu_w_1 - <w> )^4 )
    !          + ( 1 - mixt_frac ) * ( 3 * sigma_w_2^4
    !                                  + 6 * ( mu_w_2 - <w> )^2 * sigma_w_2^2
    !                                  + ( mu_w_2 - <w> )^4 ).
    !
    ! The following substitutions are made into the above equation:
    !
    ! mu_w_1 - <w> = sqrt(F_w) * sqrt( ( 1 - mixt_frac ) / mixt_frac )
    !                * sqrt( <w'^2> );
    !
    ! mu_w_2 - <w> = - sqrt(F_w) * sqrt( mixt_frac / ( 1 - mixt_frac ) )
    !                  * sqrt( <w'^2> );
    !
    ! sigma_w_1 = sqrt( coef_sigma_w_1_sqd * <w'^2> ); and
    !
    ! sigma_w_2 = sqrt( coef_sigma_w_2_sqd * <w'^2> ).
    !
    ! The equations for coef_sigma_w_1_sqd and coef_sigma_w_2_sqd are:
    !
    ! coef_sigma_w_1_sqd = ( ( zeta_w + 1 ) * ( 1 - F_w ) )
    !                      / ( ( zeta_w + 2 ) * mixt_frac ); and
    !
    ! coef_sigma_w_2_sqd = ( 1 - F_w ) / ( ( zeta_w + 2 ) * ( 1 - mixt_frac ) ).
    !
    ! The equation for <w'4> becomes:
    !
    ! <w'^4> = ( 3 * mixt_frac * coef_sigma_w_1_sqd^2
    !            + 6 * F_w * ( 1 - mixt_frac ) * coef_sigma_w_1_sqd
    !            + F_w^2 * ( 1 - mixt_frac )^2 / mixt_frac
    !            + 3 * ( 1 - mixt_frac ) * coef_sigma_w_2_sqd^2
    !            + 6 * F_w * mixt_frac * coef_sigma_w_2_sqd
    !            + F_w^2 * mixt_frac^2 / ( 1 - mixt_frac ) ) * <w'^2>^2.
    !
    ! This equation is of the form:
    !
    ! <w'^4> = coef_wp4_implicit * <w'^2>^2;
    !
    ! where:
    !
    ! coef_wp4_implicit = 3 * mixt_frac * coef_sigma_w_1_sqd^2
    !                     + 6 * F_w * ( 1 - mixt_frac ) * coef_sigma_w_1_sqd
    !                     + F_w^2 * ( 1 - mixt_frac )^2 / mixt_frac
    !                     + 3 * ( 1 - mixt_frac ) * coef_sigma_w_2_sqd^2
    !                     + 6 * F_w * mixt_frac * coef_sigma_w_2_sqd
    !                     + F_w^2 * mixt_frac^2 / ( 1 - mixt_frac ).
    !
    ! While the <w'^4> term is found in the <w'^3> predictive equation and not
    ! the <w'^2> predictive equation, the <w'^3> and <w'^2> predictive equations
    ! are solved together.  This allows the term containing <w'^4> to be solved
    ! implicitly.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        six,   & ! Variable(s)
        three, &
        one

    use clubb_precision, only: &
        core_rknd    ! Procedure(s)

    implicit none

    ! Input Variables
    real ( kind = core_rknd ), intent(in) :: &
      mixt_frac,          & ! Mixture fraction                               [-]
      F_w,                & ! Parameter: spread of the PDF comp. means of w  [-]
      coef_sigma_w_1_sqd, & ! sigma_w_1^2 = coef_sigma_w_1_sqd * <w'^2>      [-]
      coef_sigma_w_2_sqd    ! sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>      [-]

    ! Return Variable
    real ( kind = core_rknd ) :: &
      coef_wp4_implicit    ! Coef.: <w'^4> = coef_wp4_implicit * <w'^2>^2    [-]


    ! Calculate coef_wp4_implicit.
    coef_wp4_implicit = three * mixt_frac * coef_sigma_w_1_sqd**2 &
                        + six * F_w * ( one - mixt_frac ) * coef_sigma_w_1_sqd &
                        + F_w**2 * ( one - mixt_frac )**2 / mixt_frac &
                        + three * ( one - mixt_frac ) * coef_sigma_w_2_sqd**2 &
                        + six * F_w * mixt_frac * coef_sigma_w_2_sqd &
                        + F_w**2 * mixt_frac**2 / ( one - mixt_frac )


    return

  end function calculate_coef_wp4_implicit

  !=============================================================================
  elemental function calc_coef_wp2xp_implicit( wp2, mixt_frac, F_w, &
                                               coef_sigma_w_1_sqd,  &
                                               coef_sigma_w_2_sqd   ) &
  result( coef_wp2xp_implicit )

    ! Description:
    ! The predictive equation for <w'x'> contains a turbulent advection term of
    ! the form:
    !
    ! - ( 1 / rho_ds ) * d ( rho_ds * <w'^2 x'> ) / dz;
    !
    ! where z is height, rho_ds is the dry, base-state density, and <w'^2 x'> is
    ! calculated by integrating over the PDF.  The equation for <w'^2 x'> is:
    !
    ! <w'^2 x'>
    ! = mixt_frac * ( ( mu_x_1 - <x> ) * ( ( mu_w_1 - <w> )^2 + sigma_w_1^2 )
    !                 + 2 * corr_w_x_1 * sigma_w_1 * sigma_x_1
    !                   * ( mu_w_1 - <w> ) )
    !   + ( 1 - mixt_frac ) * ( ( mu_x_2 - <x> )
    !                           * ( ( mu_w_2 - <w> )^2 + sigma_w_2^2 )
    !                           + 2 * corr_w_x_2 * sigma_w_2 * sigma_x_2
    !                             * ( mu_w_2 - <w> ) ).
    !
    ! The following substitutions are made into the above equation:
    !
    ! mu_w_1 - <w> = sqrt(F_w) * sqrt( ( 1 - mixt_frac ) / mixt_frac )
    !                * sqrt( <w'^2> );
    !
    ! mu_w_2 - <w> = - sqrt(F_w) * sqrt( mixt_frac / ( 1 - mixt_frac ) )
    !                  * sqrt( <w'^2> );
    !
    ! mu_x_1 - <x> = sqrt( ( 1 - mixt_frac ) / mixt_frac )
    !                * <w'x'> / sqrt( F_w * <w'^2> );
    !
    ! mu_x_2 - <x> = - sqrt( mixt_frac / ( 1 - mixt_frac ) )
    !                  * <w'x'> / sqrt( F_w * <w'^2> );
    !
    ! sigma_w_1 = sqrt( coef_sigma_w_1_sqd * <w'^2> );
    !
    ! sigma_w_2 = sqrt( coef_sigma_w_2_sqd * <w'^2> );
    !
    ! sigma_x_1 = sqrt( coef_sigma_x_1_sqd * <x'^2> ); and
    !
    ! sigma_x_2 = sqrt( coef_sigma_x_2_sqd * <x'^2> ).
    !
    ! The equations for coef_sigma_w_1_sqd and coef_sigma_w_2_sqd are:
    !
    ! coef_sigma_w_1_sqd = ( ( zeta_w + 1 ) * ( 1 - F_w ) )
    !                      / ( ( zeta_w + 2 ) * mixt_frac ); and
    !
    ! coef_sigma_w_2_sqd = ( 1 - F_w ) / ( ( zeta_w + 2 ) * ( 1 - mixt_frac ) ).
    !
    ! The equations for coef_sigma_x_1_sqd and coef_sigma_x_2_sqd are:
    !
    ! coef_sigma_x_1_sqd
    ! = 1 + sqrt( ( 1 - mixt_frac ) / mixt_frac )
    !       * Skx * sqrt( F_w * <w'^2> * <x'^2> ) / ( 3 * <w'x'> )
    !   - ( ( 1 + mixt_frac ) / mixt_frac )
    !     * <w'x'>^2 / ( 3 * F_w * <w'^2> * <x'^2> ); and
    !
    ! coef_sigma_x_2_sqd
    ! = 1 - sqrt( mixt_frac / ( 1 - mixt_frac ) )
    !       * Skx * sqrt( F_w * <w'^2> * <x'^2> ) / ( 3 * <w'x'> )
    !   + ( ( mixt_frac - 2 ) / ( 1 - mixt_frac ) )
    !     * <w'x'>^2 / ( 3 * F_w * <w'^2> * <x'^2> ).
    !
    ! Additionally, corr_w_x_1 = corr_w_x_2 = 0.
    !
    ! The equation for <w'^2 x'> becomes:
    !
    ! <w'^2 x'> = sqrt( mixt_frac * ( 1 - mixt_frac ) )
    !             * ( F_w * ( ( 1 - mixt_frac ) / mixt_frac
    !                         - mixt_frac / ( 1 - mixt_frac ) )
    !                 + coef_sigma_w_1_sqd - coef_sigma_w_2_sqd )
    !             * sqrt( <w'^2> / F_w ) * <w'x'>.
    !
    ! This equation is of the form:
    !
    ! <w'^2 x'> = coef_wp2xp_implicit * <w'x'>;
    !
    ! where:
    !
    ! coef_wp2xp_implicit = sqrt( mixt_frac * ( 1 - mixt_frac ) )
    !                       * ( F_w * ( ( 1 - mixt_frac ) / mixt_frac
    !                                   - mixt_frac / ( 1 - mixt_frac ) )
    !                           + coef_sigma_w_1_sqd - coef_sigma_w_2_sqd )
    !                       * sqrt( <w'^2> / F_w ).
    !
    ! In the special case that F_w = 0, <w'x'> must have a value of 0, in which
    ! case mu_x_1 - <x> = mu_x_2 - <x> = 0, and <w'^2 x'> = 0.  The value of
    ! coef_wp2xp_implicit when F_w = 0 and <w'x'> = 0 can be calculated since
    ! Skw must also have a value of 0 when F_w = 0.  When F_w = 0 and Skw = 0,
    ! mixt_frac = ( zeta_w + 1 ) / ( zeta_w + 2 ).  When this happens,
    ! coef_sigma_w_1_sqd - coef_sigma_w_2_sqd = 0.  The equation for
    ! coef_wp2xp_implicit becomes:
    !
    ! coef_wp2xp_implicit = sqrt( mixt_frac * ( 1 - mixt_frac ) )
    !                       * ( ( 1 - mixt_frac ) / mixt_frac
    !                           - mixt_frac / ( 1 - mixt_frac ) )
    !                       * sqrt( F_w * <w'^2> );
    !
    ! and since F_w = 0, coef_wp2xp_implicit = 0.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,  & ! Variable(s)
        zero

    use clubb_precision, only: &
        core_rknd    ! Procedure(s)

    implicit none

    ! Input Variables
    real ( kind = core_rknd ), intent(in) :: &
      wp2,                & ! Variance of w (overall)                  [m^2/s^2]
      mixt_frac,          & ! Mixture fraction                               [-]
      F_w,                & ! Parameter: spread of the PDF comp. means of w  [-]
      coef_sigma_w_1_sqd, & ! sigma_w_1^2 = coef_sigma_w_1_sqd * <w'^2>      [-]
      coef_sigma_w_2_sqd    ! sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>      [-]

    ! Return Variable
    ! Coefficient: <w'^2 x'> = coef_wp2xp_implicit * <w'x'>
    real ( kind = core_rknd ) :: &
      coef_wp2xp_implicit    ! Coefficient that is multiplied by <w'x'>    [m/s]


    ! Calculate coef_wp2xp_implicit.
    if ( F_w > 0 ) then

       coef_wp2xp_implicit &
       = sqrt( mixt_frac * ( one - mixt_frac ) ) &
         * ( F_w * ( ( one - mixt_frac ) / mixt_frac &
                     - mixt_frac / ( one - mixt_frac ) ) &
             + coef_sigma_w_1_sqd - coef_sigma_w_2_sqd ) &
         * sqrt( wp2 / F_w )

    else ! F_w = 0

       coef_wp2xp_implicit = zero

    endif


    return

  end function calc_coef_wp2xp_implicit

  !=============================================================================
  elemental subroutine calc_coefs_wpxp2_semiimpl( wp2, wpxp,           & ! In
                                                  mixt_frac, F_w,      & ! In
                                                  coef_sigma_x_1_sqd,  & ! In
                                                  coef_sigma_x_2_sqd,  & ! In
                                                  coef_wpxp2_implicit, & ! Out
                                                  term_wpxp2_explicit )  ! Out

    ! Description:
    ! The predictive equation for <x'^2> contains a turbulent advection term of
    ! the form:
    !
    ! - ( 1 / rho_ds ) * d ( rho_ds * <w'x'^2> ) / dz;
    !
    ! where z is height, rho_ds is the dry, base-state density, and <w'x'^2> is
    ! calculated by integrating over the PDF.  The equation for <w'x'^2> is:
    !
    ! <w'x'^2>
    ! = mixt_frac * ( ( mu_w_1 - <w> ) * ( ( mu_x_1 - <x> )^2 + sigma_x_1^2 )
    !                 + 2 * corr_w_x_1 * sigma_w_1 * sigma_x_1
    !                   * ( mu_x_1 - <x> ) )
    !   + ( 1 - mixt_frac ) * ( ( mu_w_2 - <w> )
    !                           * ( ( mu_x_2 - <x> )^2 + sigma_x_2^2 )
    !                           + 2 * corr_w_x_2 * sigma_w_2 * sigma_x_2
    !                             * ( mu_x_2 - <x> ) ).
    !
    ! The following substitutions are made into the above equation:
    !
    ! mu_w_1 - <w> = sqrt(F_w) * sqrt( ( 1 - mixt_frac ) / mixt_frac )
    !                * sqrt( <w'^2> );
    !
    ! mu_w_2 - <w> = - sqrt(F_w) * sqrt( mixt_frac / ( 1 - mixt_frac ) )
    !                  * sqrt( <w'^2> );
    !
    ! mu_x_1 - <x> = sqrt( ( 1 - mixt_frac ) / mixt_frac )
    !                * <w'x'> / sqrt( F_w * <w'^2> );
    !
    ! mu_x_2 - <x> = - sqrt( mixt_frac / ( 1 - mixt_frac ) )
    !                  * <w'x'> / sqrt( F_w * <w'^2> );
    !
    ! sigma_w_1 = sqrt( coef_sigma_w_1_sqd * <w'^2> );
    !
    ! sigma_w_2 = sqrt( coef_sigma_w_2_sqd * <w'^2> );
    !
    ! sigma_x_1 = sqrt( coef_sigma_x_1_sqd * <x'^2> ); and
    !
    ! sigma_x_2 = sqrt( coef_sigma_x_2_sqd * <x'^2> ).
    !
    ! The equations for coef_sigma_w_1_sqd and coef_sigma_w_2_sqd are:
    !
    ! coef_sigma_w_1_sqd = ( ( zeta_w + 1 ) * ( 1 - F_w ) )
    !                      / ( ( zeta_w + 2 ) * mixt_frac ); and
    !
    ! coef_sigma_w_2_sqd = ( 1 - F_w ) / ( ( zeta_w + 2 ) * ( 1 - mixt_frac ) ).
    !
    ! The equations for coef_sigma_x_1_sqd and coef_sigma_x_2_sqd are:
    !
    ! coef_sigma_x_1_sqd
    ! = 1 + sqrt( ( 1 - mixt_frac ) / mixt_frac )
    !       * Skx * sqrt( F_w * <w'^2> * <x'^2> ) / ( 3 * <w'x'> )
    !   - ( ( 1 + mixt_frac ) / mixt_frac )
    !     * <w'x'>^2 / ( 3 * F_w * <w'^2> * <x'^2> ); and
    !
    ! coef_sigma_x_2_sqd
    ! = 1 - sqrt( mixt_frac / ( 1 - mixt_frac ) )
    !       * Skx * sqrt( F_w * <w'^2> * <x'^2> ) / ( 3 * <w'x'> )
    !   + ( ( mixt_frac - 2 ) / ( 1 - mixt_frac ) )
    !     * <w'x'>^2 / ( 3 * F_w * <w'^2> * <x'^2> ).
    !
    ! Additionally, corr_w_x_1 = corr_w_x_2 = 0.
    !
    ! The equation for <w'x'^2> becomes:
    !
    ! <w'x'^2>
    ! = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( F_w * <w'^2> )
    !   * ( coef_sigma_x_1_sqd - coef_sigma_x_2_sqd ) * <x'^2>
    !   + sqrt( mixt_frac * ( 1 - mixt_frac ) )
    !     * <w'x'>^2 / sqrt( F_w * <w'^2> )
    !     * ( ( 1 - mixt_frac ) / mixt_frac  - mixt_frac / ( 1 - mixt_frac ) ).
    !
    ! This equation is of the form:
    !
    ! <w'x'^2> = coef_wpxp2_implicit * <x'^2> + term_wpxp2_explicit;
    !
    ! where:
    !
    ! coef_wpxp2_implicit
    ! = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( F_w * <w'^2> )
    !   * ( coef_sigma_x_1_sqd - coef_sigma_x_2_sqd ); and
    !
    ! term_wpxp2_explicit
    ! = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * <w'x'>^2 / sqrt( F_w * <w'^2> )
    !   * ( ( 1 - mixt_frac ) / mixt_frac  - mixt_frac / ( 1 - mixt_frac ) ).
    !
    ! In the special case that F_w = 0, mu_w_1 - <w> = mu_w_2 - <w> = 0, and
    ! <w'x'^2> = 0.  The value of coef_wp2xp_implicit = 0 when F_w = 0.
    ! Likewise, the value of <w'x'> = 0 when F_w = 0, which makes
    ! term_wp2xp_explicit undefined in this scenario.  However, since both
    ! <w'x'^2> = 0 and coef_wp2xp_implicit = 0, term_wp2xp_explicit = 0 in this
    ! situation.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,  & ! Variable(s)
        zero

    use clubb_precision, only: &
        core_rknd    ! Procedure(s)

    implicit none

    ! Input Variables
    real ( kind = core_rknd ), intent(in) :: &
      wp2,                & ! Variance of w (overall)                  [m^2/s^2]
      wpxp,               & ! Covariance of w and x           [m/s (units vary)]
      mixt_frac,          & ! Mixture fraction                               [-]
      F_w,                & ! Parameter: spread of the PDF comp. means of w  [-]
      coef_sigma_x_1_sqd, & ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>      [-]
      coef_sigma_x_2_sqd    ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>      [-]

    ! Output Variables
    ! Coefs.: <w'x'^2> = coef_wpxp2_implicit * <x'^2> + term_wpxp2_explicit
    real ( kind = core_rknd ), intent(out) :: &
      coef_wpxp2_implicit, & ! Coefficient that is multiplied by <x'^2>  [m/s]
      term_wpxp2_explicit    ! Term that is on the RHS    [m/s (units vary)^2]


    ! Calculate coef_wpxp2_implicit and term_wpxp2_explicit.
    if ( F_w > 0 .and. wp2 > 0 ) then

       coef_wpxp2_implicit &
       = sqrt( mixt_frac * ( one - mixt_frac ) ) * sqrt( F_w * wp2 ) &
         * ( coef_sigma_x_1_sqd - coef_sigma_x_2_sqd )

       term_wpxp2_explicit &
       = sqrt( mixt_frac * ( one - mixt_frac ) ) * wpxp**2 / sqrt( F_w * wp2 ) &
         * ( ( one - mixt_frac ) / mixt_frac - mixt_frac / ( one - mixt_frac ) )

    else ! F_w = 0 or wp2 = 0

       coef_wpxp2_implicit = zero
       term_wpxp2_explicit = zero

    endif ! F_w > 0 and wp2 > 0


    return

  end subroutine calc_coefs_wpxp2_semiimpl

  !=============================================================================
  elemental subroutine calc_coefs_wpxpyp_semiimpl( wp2, wpxp, wpyp,      & ! In
                                                   mixt_frac, F_w,       & ! In
                                                   coef_sigma_x_1_sqd,   & ! In
                                                   coef_sigma_x_2_sqd,   & ! In
                                                   coef_sigma_y_1_sqd,   & ! In
                                                   coef_sigma_y_2_sqd,   & ! In
                                                   coef_wpxpyp_implicit, & ! Out
                                                   term_wpxpyp_explicit  ) ! Out

    ! Description:
    ! The predictive equation for <w'x'y'> contains a turbulent advection term
    ! of the form:
    !
    ! - ( 1 / rho_ds ) * d ( rho_ds * <w'x'y'> ) / dz;
    !
    ! where z is height, rho_ds is the dry, base-state density, and <w'x'y'> is
    ! calculated by integrating over the PDF.  The equation for <w'x'y'> is:
    !
    ! <w'x'y'>
    ! = mixt_frac
    !   * ( ( mu_w_1 - <w> ) * ( mu_x_1 - <x> ) * ( mu_y_1 - <y> )
    !       + corr_x_y_1 * sigma_x_1 * sigma_y_1 * ( mu_w_1 - <w> )
    !       + corr_w_y_1 * sigma_w_1 * sigma_y_1 * ( mu_x_1 - <x> )
    !       + corr_w_x_1 * sigma_w_1 * sigma_x_1 * ( mu_y_1 - <y> ) )
    !   + ( 1 - mixt_frac )
    !     * ( ( mu_w_2 - <w> ) * ( mu_x_2 - <x> ) * ( mu_y_2 - <y> )
    !         + corr_x_y_2 * sigma_x_2 * sigma_y_2 * ( mu_w_2 - <w> )
    !         + corr_w_y_2 * sigma_w_2 * sigma_y_2 * ( mu_x_2 - <x> )
    !         + corr_w_x_2 * sigma_w_2 * sigma_x_2 * ( mu_y_2 - <y> ) ).
    !
    ! The following substitutions are made into the above equation:
    !
    ! mu_w_1 - <w> = sqrt(F_w) * sqrt( ( 1 - mixt_frac ) / mixt_frac )
    !                * sqrt( <w'^2> );
    !
    ! mu_w_2 - <w> = - sqrt(F_w) * sqrt( mixt_frac / ( 1 - mixt_frac ) )
    !                  * sqrt( <w'^2> );
    !
    ! mu_x_1 - <x> = sqrt( ( 1 - mixt_frac ) / mixt_frac )
    !                * <w'x'> / sqrt( F_w * <w'^2> );
    !
    ! mu_x_2 - <x> = - sqrt( mixt_frac / ( 1 - mixt_frac ) )
    !                  * <w'x'> / sqrt( F_w * <w'^2> );
    !
    ! mu_y_1 - <y> = sqrt( ( 1 - mixt_frac ) / mixt_frac )
    !                * <w'y'> / sqrt( F_w * <w'^2> );
    !
    ! mu_y_2 - <y> = - sqrt( mixt_frac / ( 1 - mixt_frac ) )
    !                  * <w'y'> / sqrt( F_w * <w'^2> );
    !
    ! sigma_w_1 = sqrt( coef_sigma_w_1_sqd * <w'^2> );
    !
    ! sigma_w_2 = sqrt( coef_sigma_w_2_sqd * <w'^2> );
    !
    ! sigma_x_1 = sqrt( coef_sigma_x_1_sqd * <x'^2> );
    !
    ! sigma_x_2 = sqrt( coef_sigma_x_2_sqd * <x'^2> );
    !
    ! sigma_y_1 = sqrt( coef_sigma_y_1_sqd * <y'^2> ); and
    !
    ! sigma_y_2 = sqrt( coef_sigma_y_2_sqd * <y'^2> ).
    !
    ! The equations for coef_sigma_w_1_sqd and coef_sigma_w_2_sqd are:
    !
    ! coef_sigma_w_1_sqd = ( ( zeta_w + 1 ) * ( 1 - F_w ) )
    !                      / ( ( zeta_w + 2 ) * mixt_frac ); and
    !
    ! coef_sigma_w_2_sqd = ( 1 - F_w ) / ( ( zeta_w + 2 ) * ( 1 - mixt_frac ) ).
    !
    ! The equations for coef_sigma_x_1_sqd and coef_sigma_x_2_sqd are:
    !
    ! coef_sigma_x_1_sqd
    ! = 1 + sqrt( ( 1 - mixt_frac ) / mixt_frac )
    !       * Skx * sqrt( F_w * <w'^2> * <x'^2> ) / ( 3 * <w'x'> )
    !   - ( ( 1 + mixt_frac ) / mixt_frac )
    !     * <w'x'>^2 / ( 3 * F_w * <w'^2> * <x'^2> ); and
    !
    ! coef_sigma_x_2_sqd
    ! = 1 - sqrt( mixt_frac / ( 1 - mixt_frac ) )
    !       * Skx * sqrt( F_w * <w'^2> * <x'^2> ) / ( 3 * <w'x'> )
    !   + ( ( mixt_frac - 2 ) / ( 1 - mixt_frac ) )
    !     * <w'x'>^2 / ( 3 * F_w * <w'^2> * <x'^2> ).
    !
    ! The equations for coef_sigma_y_1_sqd and coef_sigma_y_2_sqd are:
    !
    ! coef_sigma_y_1_sqd
    ! = 1 + sqrt( ( 1 - mixt_frac ) / mixt_frac )
    !       * Sky * sqrt( F_w * <w'^2> * <y'^2> ) / ( 3 * <w'y'> )
    !   - ( ( 1 + mixt_frac ) / mixt_frac )
    !     * <w'y'>^2 / ( 3 * F_w * <w'^2> * <y'^2> ); and
    !
    ! coef_sigma_y_2_sqd
    ! = 1 - sqrt( mixt_frac / ( 1 - mixt_frac ) )
    !       * Sky * sqrt( F_w * <w'^2> * <y'^2> ) / ( 3 * <w'y'> )
    !   + ( ( mixt_frac - 2 ) / ( 1 - mixt_frac ) )
    !     * <w'y'>^2 / ( 3 * F_w * <w'^2> * <y'^2> ).
    !
    ! Additionally, corr_w_x_1 = corr_w_x_2 = corr_w_y_1 = corr_w_y_2 = 0; and:
    !
    ! corr_x_y_1 = corr_x_y_2
    ! = ( <x'y'> - mixt_frac * ( mu_x_1 - <x> ) * ( mu_y_1 - <y> )
    !            - ( 1 - mixt_frac ) * ( mu_x_2 - <x> ) * ( mu_y_2 - <y> ) )
    !   / ( mixt_frac * sigma_x_1 * sigma_y_1
    !       + ( 1 - mixt_frac ) * sigma_x_2 * sigma_y_2 );
    !
    ! where -1 <= corr_x_y_1 = corr_x_y_2 <= 1.  This equation can be rewritten
    ! as:
    !
    ! corr_x_y_1 = corr_x_y_2
    ! = ( <x'y'>
    !     - sqrt( F_x ) * sqrt( F_y ) * sgn( <w'x'> ) * sgn( <w'y'> )
    !       * sqrt( <x'^2> ) * sqrt( <y'^2 > ) )
    !   / ( ( mixt_frac * sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
    !         + ( 1 - mixt_frac )
    !           * sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
    !       * sqrt( <x'^2> ) * sqrt( <y'^2> ) ).
    !
    ! The equation for <w'x'y'> becomes:
    !
    ! <w'x'y'>
    ! = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( F_w * <w'^2> )
    !   * ( sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
    !             - sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
    !     / ( mixt_frac * sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
    !         + ( 1 - mixt_frac )
    !           * sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
    !   * <x'y'>
    !   + sqrt( mixt_frac * ( 1 - mixt_frac ) )
    !     * <w'x'> * <w'y'> / sqrt( F_w * <w'^2> )
    !     * ( ( 1 - mixt_frac ) / mixt_frac - mixt_frac / ( 1 - mixt_frac )
    !         - ( sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
    !             - sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
    !           / ( mixt_frac * sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
    !               + ( 1 - mixt_frac )
    !                 * sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) ) )
    !
    ! This equation is of the form:
    !
    ! <w'x'y'> = coef_wpxpyp_implicit * <x'y'> + term_wpxpyp_explicit;
    !
    ! where:
    !
    ! coef_wpxpyp_implicit
    ! = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( F_w * <w'^2> )
    !   * ( sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
    !             - sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
    !     / ( mixt_frac * sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
    !         + ( 1 - mixt_frac )
    !           * sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) ); and
    !
    ! term_wpxpyp_explicit
    ! = sqrt( mixt_frac * ( 1 - mixt_frac ) )
    !   * <w'x'> * <w'y'> / sqrt( F_w * <w'^2> )
    !   * ( ( 1 - mixt_frac ) / mixt_frac - mixt_frac / ( 1 - mixt_frac )
    !       - ( sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
    !           - sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
    !         / ( mixt_frac * sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
    !             + ( 1 - mixt_frac )
    !               * sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) ) ).
    !
    ! There are also special cases for the above equations.  In the scenario
    ! that sigma_x_1 * sigma_y_1 = sigma_x_2 * sigma_y_2 = 0, and equations for
    ! coef_wpxpyp_implicit and term_wpxpyp_explicit become:
    !
    ! coef_wpxpyp_implicit
    ! = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( F_w * <w'^2> )
    !   * ( ( 1 - mixt_frac ) / mixt_frac - mixt_frac / ( 1 - mixt_frac ) ); and
    !
    ! term_wpxpyp_explicit = 0.
    !
    ! In the scenario where F_w = 0, mu_w_1 - <w> = mu_w_2 - <w> = 0, and
    ! <w'x'> = <w'y'> = 0, which means mu_x_1 - <x> = mu_x_2 - <x> = 0, as well
    ! as mu_y_1 - <y> = mu_y_2 - <y> = 0, and as a result, <w'x'y'> = 0.  When
    ! F_w = 0, coef_wpxpyp_implicit = 0.  Since <w'x'y'> = 0 and
    ! coef_wpxpyp_implicit = 0 when F_w = 0, term_wpxpyp_explicit = 0.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,  & ! Variable(s)
        zero

    use clubb_precision, only: &
        core_rknd    ! Procedure(s)

    implicit none

    ! Input Variables
    real ( kind = core_rknd ), intent(in) :: &
      wp2,                & ! Variance of w (overall)                  [m^2/s^2]
      wpxp,               & ! Covariance of w and x (overall)     [m/s(x units)]
      wpyp,               & ! Covariance of w and y (overall)     [m/s(y units)]
      mixt_frac,          & ! Mixture fraction                               [-]
      F_w,                & ! Parameter: spread of the PDF comp. means of w  [-]
      coef_sigma_x_1_sqd, & ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>      [-]
      coef_sigma_x_2_sqd, & ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>      [-]
      coef_sigma_y_1_sqd, & ! sigma_y_1^2 = coef_sigma_y_1_sqd * <y'^2>      [-]
      coef_sigma_y_2_sqd    ! sigma_y_2^2 = coef_sigma_y_2_sqd * <y'^2>      [-]

    ! Output Variables
    ! Coefs.: <w'x'y'> = coef_wpxpyp_implicit * <x'y'> + term_wpxpyp_explicit
    real ( kind = core_rknd ), intent(out) :: &
      coef_wpxpyp_implicit, & ! Coefficient that is multiplied by <x'y'>   [m/s]
      term_wpxpyp_explicit    ! Term that is on the RHS  [m/s(x units)(y units)]

    ! Local Variables
    real ( kind = core_rknd ) :: &
      coefs_factor_xy    ! Factor involving coef_sigma_... x and y coefs     [-]


    ! Calculate coef_wpxpyp_implicit and term_wpxpyp_explicit.
    if ( ( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd > zero &
           .or. coef_sigma_x_2_sqd * coef_sigma_y_2_sqd > zero ) &
          .and. F_w > zero .and. wp2 > zero ) then

       ! coefs_factor_xy
       ! = ( sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
       !     - sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
       !   / ( mixt_frac * sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
       !       + ( 1 - mixt_frac )
       !         * sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
       coefs_factor_xy &
       = ( sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd ) &
           - sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) ) &
         / ( mixt_frac * sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd ) &
             + ( one - mixt_frac ) &
               * sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )

       coef_wpxpyp_implicit &
       = sqrt( mixt_frac * ( 1 - mixt_frac ) ) &
         * sqrt( F_w * wp2 ) * coefs_factor_xy

       term_wpxpyp_explicit &
       = sqrt( mixt_frac * ( one - mixt_frac ) ) &
         * wpxp * wpyp / sqrt( F_w * wp2 ) &
         * ( ( one - mixt_frac ) / mixt_frac &
             - mixt_frac / ( one - mixt_frac ) &
             - coefs_factor_xy )

    else ! ( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd = 0
         !   and coef_sigma_x_2_sqd * coef_sigma_y_2_sqd = 0 )
         ! or F_w = 0 or wp2 = 0

       if ( F_w > 0 ) then

          coef_wpxpyp_implicit &
          = sqrt( mixt_frac * ( one - mixt_frac ) ) * sqrt( F_w * wp2 ) &
            * ( ( one - mixt_frac ) / mixt_frac &
                - mixt_frac / ( one - mixt_frac ) )

       else ! F_w = 0

          coef_wpxpyp_implicit = zero

       endif

       term_wpxpyp_explicit = zero

    endif


    return

  end subroutine calc_coefs_wpxpyp_semiimpl

  !=============================================================================

end module new_hybrid_pdf
