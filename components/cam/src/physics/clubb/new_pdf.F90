! $Id$
!===============================================================================
module new_pdf

  ! Description:
  ! The portion of CLUBB's multivariate, two-component PDF that is the
  ! trivariate, two-component normal PDF of vertical velocity (w), total water
  ! mixing ratio (rt), and liquid water potential temperature (thl).

  ! References:
  ! Griffin and Larson (2018)
  !-------------------------------------------------------------------------

  implicit none

  public :: calc_mixture_fraction,      & ! Procedure(s)
            calc_setter_var_params,     &
            calc_responder_params,      &
            calc_limits_F_x_responder,  &
            calc_coef_wp4_implicit,     &
            calc_coef_wpxp2_implicit,   &
            calc_coefs_wp2xp_semiimpl,  &
            calc_coefs_wpxpyp_semiimpl

  private :: sort_roots    ! Procedure(s)

  private

  contains

  !=============================================================================
  !
  ! DESCRIPTION OF THE METHOD FOR THE VARIABLE THAT SETS THE MIXTURE FRACTION
  ! =========================================================================
  !
  ! Many times, w has been used as the variable that sets the mixture fraction
  ! for the PDF.  There are five PDF parameters that need to be calculated,
  ! which are mu_w_1 (the mean of w is the 1st PDF component), mu_w_2 (the mean
  ! of w in the 2nd PDF component), sigma_w_1 (the standard deviation of w in
  ! the 1st PDF component), sigma_w_2 (the standard deviation of w in the 2nd
  ! PDF component), and mixt_frac (the mixture fraction, which is the weight of
  ! the 1st PDF component).  In order to solve for these five parameters, five
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
  ! sigma_w_1
  ! = sqrt( ( ( zeta_w + 1 ) * ( 1 - F_w ) )
  !         / ( ( zeta_w + 2 ) * mixt_frac ) * <w'^2> ); and
  !
  ! sigma_w_2
  ! = sqrt( ( mixt_frac * sigma_w_1^2 )
  !         / ( ( 1 - mixt_frac ) * ( 1 + zeta_w ) ) );
  !
  ! where Skw is the skewness of w, and Skw = <w'^3> / <w'^2>^(3/2).
  !
  ! This method works for all values of F_w (where 0 <= F_w <= 1) and zeta_w
  ! (where zeta_w > -1).
  !
  !
  ! Generalized equations for any variable, x, that sets the mixture fraction:
  !
  ! A slight alteration is made to the above equations in order to have any
  ! variable, x, set the mixture fraction.  The same five PDF parameters need to
  ! be calculated for the setting variable, which are mu_x_1 (the mean of x in
  ! the 1st PDF component), mu_x_2 (the mean of x in the 2nd PDF component),
  ! sigma_x_1 (the standard deviation of x in the 1st PDF component), sigma_x_2
  ! (the standard deviation of x in the 2nd PDF component), and mixt_frac (the
  ! mixture fraction).  Again, five equations are needed, and they are the
  ! equations for <x>, <x'^2>, and <x'^3> as found by integrating over the PDF,
  ! as well as the equations that involve tunable parameters F_x and zeta_x.
  ! However, the equation for F_x is multiplied by a new variable,
  ! sgn( <w'x'> ), where <w'x'> is the covariance of w and x, and sgn( <w'x'> )
  ! is given by:
  !
  ! sgn( <w'x'> ) = |  1, when <w'x'> >= 0;
  !                 | -1, when <w'x'> < 0.
  !
  ! The five equations are:
  !
  ! <x> = mixt_frac * mu_x_1 + ( 1 - mixt_frac ) * mu_x_2;
  !
  ! <x'^2> = mixt_frac * ( ( mu_x_1 - <x> )^2 + sigma_x_1^2 )
  !          + ( 1 - mixt_frac ) * ( ( mu_x_2 - <x> )^2 + sigma_x_2^2 );
  !
  ! <x'^3> = mixt_frac * ( mu_x_1 - <x> )
  !                    * ( ( mu_x_1 - <x> )^2 + 3 * sigma_x_1^2 )
  !          + ( 1 - mixt_frac ) * ( mu_x_2 - <x> )
  !                              * ( ( mu_x_2 - <x> )^2 + 3 * sigma_x_2^2 );
  !
  ! mu_x_1 - <x> = sqrt(F_x) * ( sqrt( 1 - mixt_frac ) / sqrt( mixt_frac ) )
  !                * sqrt( <x'^2> ) * sgn( <w'x'> );
  !
  ! where 0 <= F_x <= 1; and
  !
  ! 1 + zeta_x = ( mixt_frac * sigma_x_1^2 )
  !              / ( ( 1 - mixt_frac ) * sigma_x_2^2 );
  !
  ! where zeta_x > -1.
  !
  ! The only equations that are altered are the equation for mu_x_1 and the
  ! equation for mixt_frac, which now both contain a sgn( <w'x'> ).  The mu_x_2
  ! equation is not altered, but the sign of mu_x_2 - <x> will be the opposite
  ! of the sign of mu_x_1 - <x>.  The resulting equations for the five PDF
  ! parameters are:
  !
  ! mixt_frac
  ! = ( 4 * F_x^3
  !     + 18 * F_x * ( zeta_x + 1 ) * ( 1 - F_x ) / ( zeta_x + 2 )
  !     + 6 * F_x^2 * ( 1 - F_x ) / ( zeta_x + 2 )
  !     + Skx^2
  !     - Skx * sgn( <w'x'> )
  !           * sqrt( 4 * F_x^3
  !                   + 12 * F_x^2 * ( 1 - F_x )
  !                   + 36 * F_x * ( zeta_x + 1 ) * ( 1 - F_x )^2
  !                     / ( zeta_x + 2 )^2
  !                   + Skx^2 ) )
  !   / ( 2 * F_x * ( F_x - 3 )^2 + 2 * Skx^2 );
  !
  ! mu_x_1 = <x> + sqrt( F_x * ( ( 1 - mixt_frac ) / mixt_frac ) * <x'^2> )
  !                * sgn( <w'x'> );
  !
  ! mu_x_2 = <x> - ( mixt_frac / ( 1 - mixt_frac ) ) * ( mu_x_1 - <x> );
  !
  ! sigma_x_1
  ! = sqrt( ( ( zeta_x + 1 ) * ( 1 - F_x ) )
  !         / ( ( zeta_x + 2 ) * mixt_frac ) * <x'^2> ); and
  !
  ! sigma_x_2
  ! = sqrt( ( mixt_frac * sigma_x_1^2 )
  !         / ( ( 1 - mixt_frac ) * ( 1 + zeta_x ) ) );
  !
  ! where Skx is the skewness of x, and Skx = <x'^3> / <x'^2>^(3/2).
  !
  ! This method works for all values of F_x (where 0 <= F_x <= 1) and zeta_x
  ! (where zeta_x > -1).
  !
  ! When the generalized form is solved for w (x = w), sgn( <w'^2> ) = 1, and
  ! the equations are unaltered from the equations listed above for w.
  !
  !
  ! Special case:
  !
  ! When Skx = 0 and F_x = 0, the equation for mixt_frac is undefined.  The
  ! equation for mixture fraction in this scenario can be derived by using the
  ! above equation for mixture fraction and then setting Skx = 0.  The resulting
  ! equation becomes:
  !
  ! mixt_frac
  ! = ( 4 * F_x^3
  !     + 18 * F_x * ( zeta_x + 1 ) * ( 1 - F_x ) / ( zeta_x + 2 )
  !     + 6 * F_x^2 * ( 1 - F_x ) / ( zeta_x + 2 ) )
  !   / ( 2 * F_x * ( F_x - 3 )^2 ).
  !
  ! All of the terms in the numerator and denominator contain a F_x, so this
  ! equation can be rewritten as:
  !
  ! mixt_frac
  ! = ( 4 * F_x^2
  !     + 18 * ( zeta_x + 1 ) * ( 1 - F_x ) / ( zeta_x + 2 )
  !     + 6 * F_x * ( 1 - F_x ) / ( zeta_x + 2 ) )
  !   / ( 2 * ( F_x - 3 )^2 ).
  !
  ! Now setting F_x = 0, the equation becomes:
  !
  ! mixt_frac = ( 18 * ( zeta_x + 1 ) / ( zeta_x + 2 ) ) / 18;
  !
  ! which can be rewritten as:
  !
  ! mixt_frac = ( zeta_x + 1 ) / ( zeta_x + 2 ).
  !
  ! When F_x = 0, Skx must have a value of 0 in order for the PDF to function
  ! correctly.  When F_x = 0, mu_x_1 = mu_x_2.  When the two PDF component means
  ! are equal to each other (and to the overall mean, <x>), the only value of
  ! Skx that can be represented is a value of 0.  In the equation for mixture
  ! fraction, when F_x is set to 0, but | Skx | > 0, mixt_frac will either have
  ! a value of 0 or 1, depending on whether Skx is positive or negative,
  ! respectively.
  !
  ! The value of F_x should be set as a function of Skx.  The value F_x should
  ! go toward 0 as | Skx | (or Skx^2) goes toward 0.  The value of F_x should
  ! go toward 1 as | Skx | (or Skx^2) goes to infinity.
  !
  !
  ! Tunable parameters:
  !
  ! 1) F_x:  This parameter controls the spread of the PDF component means.  The
  !          range of this parameter is 0 <= F_x <= 1.  When F_x = 0, the two
  !          PDF component means (mu_x_1 and mu_x_2) are equal to each other
  !          (and Skx must equal 0).  All of the variance of x is accounted for
  !          by the PDF component standard deviations (sigma_x_1 and sigma_x_2).
  !          When F_x = 1, mu_x_1 and mu_x_2 are spread as far apart as they can
  !          be.  Both PDF component standard deviations (sigma_x_1 and
  !          sigma_x_2) are equal to 0, and all of the variance of x is
  !          accounted for by the spread of the PDF component means.
  !
  !          When sigma_x_1 = sigma_x_2 = 0, the equation for <x'^2> becomes:
  !
  !          <x'^2> = mixt_frac * ( mu_x_1 - <x> )^2
  !                   + ( 1 - mixt_frac ) * ( mu_x_2 - <x> )^2.
  !
  !          Substituting the equation for <x> into the above equation for
  !          mu_x_2 - <x>, the above equation becomes:
  !
  !          <x'^2> = ( mixt_frac / ( 1 - mixt_frac ) ) * ( mu_x_1 - <x> )^2;
  !
  !          which can be rewritten as:
  !
  !          ( mu_x_1 - <x> )^2 = ( ( 1 - mixt_frac ) / mixt_frac ) * <x'^2>.
  !
  !          Taking the square root of the above equation:
  !
  !          mu_x_1 - <x> = +/- ( sqrt( 1 - mixt_frac ) / sqrt(mixt_frac) )
  !                             * sqrt( <x'^2> ).
  !
  !          This equation can be compared to the equation for mu_x_1 - <x> in
  !          the set of 5 equations, which is:
  !
  !          mu_x_1 - <x>
  !          = sqrt(F_x) * ( sqrt( 1 - mixt_frac ) / sqrt( mixt_frac ) )
  !            * sqrt( <x'^2> ) * sgn( <w'x'> ).
  !
  !          The above equations give another example of the meaning of F_x.
  !          The value of sqrt(F_x) is ratio of mu_x_1 - <x> to its maximum
  !          value (or minimum value, depending on sgn( <w'x'> )), which is:
  !
  !          sqrt( ( ( 1 - mixt_frac ) / mixt_frac ) * <x'^2> ) * sgn( <w'x'> ).
  !
  !
  ! 2) zeta_x:  This parameter controls the size of the PDF component standard
  !             deviations compared to each other.  The equation for zeta_x is:
  !
  !             1 + zeta_x = ( mixt_frac * sigma_x_1^2 )
  !                          / ( ( 1 - mixt_frac ) * sigma_x_2^2 ).
  !
  !             When zeta_x > 0, mixt_frac * sigma_x_1^2 increases at the
  !             expense of ( 1 - mixt_frac ) * sigma_x_2^2, which decreases in
  !             this variance-preserving equation set.  When zeta_x = 0, then
  !             mixt_frac * sigma_x_1^2 = ( 1 - mixt_frac ) * sigma_x_2^2.
  !             When -1 < zeta_x < 0, ( 1 - mixt_frac ) * sigma_x_2^2 increases
  !             at the expense of mixt_frac * sigma_x_1^2, which decreases.  As
  !             a result, greater values of zeta_x cause the 1st PDF component
  !             to become broader while the 2nd PDF component becomes narrower,
  !             and smaller values of zeta_x cause the 1st PDF component to
  !             become narrower while the 2nd PDF component becomes broader.
  !
  !             Symmetry
  !
  !             When zeta_x = 0, the PDF is always symmetric.  In other words,
  !             the PDF at any positive value of Skx (for example, Skx = 2.5)
  !             will look like a mirror-image (reflection across the y-axis)
  !             of the PDF at a negative value of Skx of the same magnitude (in
  !             this example, Skx = -2.5).  However, when zeta_x /= 0, the PDF
  !             loses this quality and is not symmetric.
  !
  !             When symmetry is desired at values of zeta_x besides zeta_x = 0,
  !             the solution is to turn zeta_x into a function of Skx.  A basic
  !             example of a zeta_x skewness equation that produces a symmetric
  !             PDF for values of zeta_x other than 0 is:
  !
  !             zeta_x = | zeta_x_in,                      when Skx >= 0;
  !                      | ( 1 / ( 1 + zeta_x_in ) ) - 1,  when Skx < 0.
  !
  !
  ! Notes:
  !
  ! When F_x = 0 (which can only happen when Skx = 0), mu_x_1 = mu_x_2, and
  ! mixt_frac = ( zeta_x + 1 ) / ( zeta_x + 2 ).  When these equations are
  ! substituted into the equations for sigma_x_1 and sigma_x_2, the result is
  ! sigma_x_1 = sigma_x_2 = sqrt( <x'^2> ).  This means that the distribution
  ! becomes a single Gaussian when F_x = 0 (and Skx = 0).  This happens
  ! regardless of the value of zeta_x.
  !
  ! The equations for the PDF component means and standard deviations can also
  ! be written as:
  !
  ! mu_x_1 = <x> + sqrt( F_x * ( ( 1 - mixt_frac ) / mixt_frac ) * <x'^2> )
  !                * sgn( <w'x'> );
  !
  ! mu_x_2 = <x> - sqrt( F_x * ( mixt_frac / ( 1 - mixt_frac ) ) * <x'^2> )
  !                * sgn( <w'x'> );
  !
  ! sigma_x_1 = sqrt( coef_sigma_x_1_sqd * <x'^2> ); and
  !
  ! sigma_x_2 = sqrt( coef_sigma_x_2_sqd * <x'^2> ); where
  !
  ! coef_sigma_x_1_sqd = ( ( zeta_x + 1 ) * ( 1 - F_x ) )
  !                      / ( ( zeta_x + 2 ) * mixt_frac ); and
  !
  ! coef_sigma_x_2_sqd = ( 1 - F_x ) / ( ( zeta_x + 2 ) * ( 1 - mixt_frac ) ).
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
  function calc_mixture_fraction( Skx, F_x, zeta_x, sgn_wpxp ) &
  result( mixt_frac )

    ! Description:
    ! Calculates mixture fraction.

    ! References:
    ! Griffin and Larson (2018)
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        thirty_six, & ! Constant(s)
        eighteen,   &
        twelve,     &
        six,        &
        four,       &
        three,      &
        two,        &
        one,        &
        zero,       &
        fstderr

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      Skx,      & ! Skewness of x                                            [-]
      F_x,      & ! Parameter for the spread of the PDF component means of x [-]
      zeta_x,   & ! Parameter for the PDF component variances of x           [-]
      sgn_wpxp    ! Sign of the covariance of w and x                        [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz) :: &
      mixt_frac    ! Mixture fraction    [-]

    ! Local Variable
    ! Flag that turns off when conditions aren't right for calculating mixt_frac
    logical, dimension(gr%nz) :: &
      l_calc_mixt_frac


    ! Initialize l_calc_mixt_frac
    l_calc_mixt_frac = .true.

    ! Calculate mixture fraction, which is the weight of the 1st PDF component.
    ! The 2nd PDF component has a weight of 1 - mixt_frac.
    where ( F_x > zero )

       mixt_frac &
       = ( four * F_x**3 &
           + eighteen * F_x &
             * ( zeta_x + one ) * ( one - F_x ) / ( zeta_x + two ) &
           + six * F_x**2 * ( one - F_x ) / ( zeta_x + two ) &
           + Skx**2 &
           - Skx * sgn_wpxp * sqrt( four * F_x**3 &
                                    + twelve * F_x**2 * ( one - F_x ) &
                                    + thirty_six * F_x &
                                      * ( zeta_x + one ) * ( one - F_x )**2 &
                                      / ( zeta_x + two )**2 &
                                    + Skx**2 ) ) &
         / ( two * F_x * ( F_x - three )**2 + two * Skx**2 )

    elsewhere ! F_x = 0

       where ( abs( Skx ) > zero )

          l_calc_mixt_frac = .false.

       elsewhere ! Skx = 0

          mixt_frac = ( zeta_x + one ) / ( zeta_x + two )

       endwhere ! | Skx | > 0

    endwhere ! F_x > 0


    if ( any( .not. l_calc_mixt_frac ) ) then
       write(fstderr,*) "Mixture fraction cannot be calculated."
       write(fstderr,*) "The value of F_x must be greater than 0 when " &
                        // "| Skx | > 0."
       stop
    endif ! any( .not. l_valid_mixt_frac )


    return

  end function calc_mixture_fraction

  !=============================================================================
  subroutine calc_setter_var_params( xm, xp2, Skx, sgn_wpxp,    & ! In
                                     F_x, zeta_x,               & ! In
                                     mu_x_1, mu_x_2, sigma_x_1, & ! Out
                                     sigma_x_2, mixt_frac,      & ! Out
                                     coef_sigma_x_1_sqd,        & ! Out
                                     coef_sigma_x_2_sqd         ) ! Out

    ! Description:
    ! Calculates the PDF component means, the PDF component standard deviations,
    ! and the mixture fraction for the variable that sets the PDF.

    ! References:
    ! Griffin and Larson (2018)
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        two, & ! Variable(s)
        one

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      xm,       & ! Mean of x (overall)                             [units vary]
      xp2,      & ! Variance of x (overall)                     [(units vary)^2]
      Skx,      & ! Skewness of x                                            [-]
      sgn_wpxp, & ! Sign of the covariance of w and x (overall)              [-]
      F_x,      & ! Parameter for the spread of the PDF component means of x [-]
      zeta_x      ! Parameter for the PDF component variances of x           [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      mu_x_1,    & ! Mean of x (1st PDF component)                  [units vary]
      mu_x_2,    & ! Mean of x (2nd PDF component)                  [units vary]
      sigma_x_1, & ! Standard deviation of x (1st PDF component)    [units vary]
      sigma_x_2, & ! Standard deviation of x (2nd PDF component)    [units vary]
      mixt_frac    ! Mixture fraction                                        [-]

    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      coef_sigma_x_1_sqd, & ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>      [-]
      coef_sigma_x_2_sqd    ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>      [-]


    ! Calculate the mixture fraction.
    mixt_frac = calc_mixture_fraction( Skx, F_x, zeta_x, sgn_wpxp )

    ! Calculate the mean of x in the 1st PDF component.
    mu_x_1 = xm + sqrt( F_x * ( ( one - mixt_frac ) / mixt_frac ) * xp2 ) &
                  * sgn_wpxp

    ! Calculate the mean of x in the 2nd PDF component.
    mu_x_2 = xm - ( mixt_frac / ( one - mixt_frac ) ) * ( mu_x_1 - xm )

    ! Calculate the standard deviation of x in the 1st PDF component.
    ! sigma_x_1 = sqrt( ( ( zeta_x + 1 ) * ( 1 - F_x ) )
    !                   / ( ( zeta_x + 2 ) * mixt_frac ) * <x'^2> )
    coef_sigma_x_1_sqd = ( ( zeta_x + one ) * ( one - F_x ) ) &
                         / ( ( zeta_x + two ) * mixt_frac )

    sigma_x_1 = sqrt( coef_sigma_x_1_sqd * xp2 )

    ! Calculate the standard deviation of x in the 2nd PDF component.
    ! sigma_x_2 = sqrt( ( mixt_frac * sigma_x_1^2 )
    !                   / ( ( 1 - mixt_frac ) * ( 1 + zeta_x ) ) )
    !           = sqrt( ( 1 - F_x )
    !                   / ( ( zeta_x + 2 ) * ( 1 - mixt_frac ) ) * <x'^2> )
    coef_sigma_x_2_sqd = ( one - F_x ) &
                         / ( ( zeta_x + two ) * ( one - mixt_frac ) )

    sigma_x_2 = sqrt( coef_sigma_x_2_sqd * xp2 )


    return

  end subroutine calc_setter_var_params

  !=============================================================================
  !
  ! DESCRIPTION OF THE METHOD FOR EACH RESPONDING VARIABLE
  ! ======================================================
  !
  ! In order to find equations for the four PDF parameters for each responding
  ! variable, which are mu_x_1, mu_x_2, sigma_x_1, and sigma_x_2 (where x stands
  ! for a responding variable here), four equations are needed.  These four
  ! equations are the equations for <x>, <x'^2>, and <x'^3> as found by
  ! integrating over the PDF.  Additionally, one more equation, which involves
  ! a tunable parameter F_x, and which is used to help control the spread of the
  ! PDF component means, is used in this equation set.  The four equations are:
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
  ! mu_x_1 - <x> = sqrt(F_x) * ( sqrt( 1 - mixt_frac ) / sqrt( mixt_frac ) )
  !                * sqrt( <x'^2> ) * sgn( <w'x'> );
  !
  ! where 0 <= F_x <= 1, and where sgn( <w'x'> ) is given by:
  !
  ! sgn( <w'x'> ) = |  1, when <w'x'> >= 0;
  !                 | -1, when <w'x'> < 0.
  !
  ! The resulting equations for the four PDF parameters are:
  !
  ! mu_x_1 = <x> + sqrt( F_x * ( ( 1 - mixt_frac ) / mixt_frac ) * <x'^2> )
  !                * sgn( <w'x'> );
  !
  ! mu_x_2 = <x> - ( mixt_frac / ( 1 - mixt_frac ) ) * ( mu_x_1 - <x> );
  !
  ! sigma_x_1^2
  ! = ( ( sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> )
  !       - ( 1 + mixt_frac ) * F_x^1.5 + 3 * mixt_frac * sqrt( F_x ) )
  !     / ( 3 * mixt_frac * sqrt( F_x ) ) )
  !   * <x'^2>; and
  !
  ! sigma_x_2^2 = ( ( 1 - F_x ) / ( 1 - mixt_frac ) ) * <x'^2>
  !               - ( mixt_frac / ( 1 - mixt_frac ) ) * sigma_x_1^2;
  !
  ! where Skx is the skewness of x, and Skx = <x'^3> / <x'^2>^(3/2).
  !
  !
  ! Special case:
  !
  ! When Skx = 0 and F_x = 0, the equations for sigma_x_1^2 and sigma_x_2^2 are
  ! both undefined.  The equations for sigma_x_1^2 and sigma_x_2^2 in this
  ! scenario can be derived by using the above equations for sigma_x_1^2 and
  ! sigma_x_2^2 and then setting Skx = 0.  The resulting equation for
  ! sigma_x_1^2 becomes:
  !
  ! sigma_x_1^2
  ! = ( ( - ( 1 + mixt_frac ) * F_x^1.5 + 3 * mixt_frac * sqrt( F_x ) )
  !     / ( 3 * mixt_frac * sqrt( F_x ) ) )
  !   * <x'^2>.
  !
  ! All of the terms in the numerator and denominator contain a sqrt( F_x ),
  ! so this equation can be rewritten as:
  !
  ! sigma_x_1^2
  ! = ( ( - ( 1 + mixt_frac ) * F_x + 3 * mixt_frac ) / ( 3 * mixt_frac ) )
  !   * <x'^2>.
  !
  ! Now setting F_x = 0, the equation becomes:
  !
  ! sigma_x_1^2 = ( ( 3 * mixt_frac ) / ( 3 * mixt_frac ) ) * <x'^2>;
  !
  ! which can be rewritten as:
  !
  ! sigma_x_1^2 = <x'^2>.
  !
  ! Substituting the equation for sigma_x_1^2 into the equation for sigma_x_2^2,
  ! and also setting F_x = 0, the equation for sigma_x_2^2 becomes:
  !
  ! sigma_x_2^2 = ( 1 / ( 1 - mixt_frac ) ) * <x'^2>
  !               - ( mixt_frac / ( 1 - mixt_frac ) ) * <x'^2>;
  !
  ! which can be rewritten as:
  !
  ! sigma_x_2^2
  ! = ( ( 1 / ( 1 - mixt_frac ) ) - ( mixt_frac / ( 1 - mixt_frac ) ) )
  !   * <x'^2>.
  !
  ! This equation becomes:
  !
  ! sigma_x_2^2 = ( ( 1 - mixt_frac ) / ( 1 - mixt_frac ) ) * <x'^2>;
  !
  ! which can be rewritten as:
  !
  ! sigma_x_2^2 = <x'^2>.
  !
  ! When F_x = 0, Skx must have a value of 0 in order for the PDF to function
  ! correctly.  When F_x = 0, mu_x_1 = mu_x_2.  When the two PDF component means
  ! are equal to each other (and to the overall mean, <x>), the only value of
  ! Skx that can be represented is a value of 0.  The equations that place
  ! limits on F_x for a responding variable (below) calculate the minimum
  ! allowable value of F_x to be greater than 0 when | Skx | > 0.
  !
  ! The value of F_x should be set as a function of Skx.  The value F_x should
  ! go toward 0 as | Skx | (or Skx^2) goes toward 0.  The value of F_x should
  ! go toward 1 as | Skx | (or Skx^2) goes to infinity.  However, the value of
  ! F_x must also be between the minimum and maximum allowable values of F_x for
  ! a responding variable (below).
  !
  !
  ! Tunable parameter:
  !
  ! F_x:  This parameter controls the spread of the PDF component means.  The
  !       range of this parameter is 0 <= F_x <= 1.  When F_x = 0, the two PDF
  !       component means (mu_x_1 and mu_x_2) are equal to each other (and Skx
  !       must equal 0).  All of the variance of x is accounted for by the PDF
  !       component standard deviations (sigma_x_1 and sigma_x_2).  When
  !       F_x = 1, mu_x_1 and mu_x_2 are spread as far apart as they can be.
  !       Both PDF component standard deviations (sigma_x_1 and sigma_x_2) are
  !       equal to 0, and all of the variance of x is accounted for by the
  !       spread of the PDF component means.
  !
  !
  ! Limits on F_x:
  !
  ! Since the PDF parameters for this variable need to work with the mixture
  ! fraction that has been provided by the setting variable, the method does
  ! not work for all values of F_x and Skx.  However, the limits of Skx and F_x
  ! can always be calculated.  The limits are based on keeping the values of
  ! sigma_x_1 and sigma_x_2 greater than or equal to 0.  The equation for
  ! keeping the value of sigma_x_1 greater than or equal to 0 is:
  !
  ! - ( 1 + mixt_frac ) * sqrt( F_x )^3 + 3 * mixt_frac * sqrt( F_x )
  ! + sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> ) >= 0.
  !
  ! The roots of sqrt( F_x ) can be solved by an equation of the form:
  !
  ! A * sqrt( F_x )^3 + B * sqrt( F_x )^2 + C * sqrt( F_x ) + D = 0;
  !
  ! where:
  !
  ! A = - ( 1 + mixt_frac );
  ! B = 0;
  ! C = 3 * mixt_frac; and
  ! D = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> ).
  !
  ! The equation for keeping the value of sigma_x_2 greater than or equal to 0
  ! is:
  !
  ! - ( 2 - mixt_frac ) * sqrt( F_x )^3 + 3 * ( 1 - mixt_frac ) * sqrt( F_x )
  ! - sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> ) >= 0.
  !
  ! The roots of sqrt( F_x ) can be solved by an equation of the form:
  !
  ! A * sqrt( F_x )^3 + B * sqrt( F_x )^2 + C * sqrt( F_x ) + D = 0;
  !
  ! where:
  !
  ! A = - ( 2 - mixt_frac );
  ! B = 0;
  ! C = 3 * ( 1 - mixt_frac ); and
  ! D = - sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> ).
  !
  ! After careful analysis of the above equations, the following properties
  ! emerge:
  !
  ! When Skx * sgn( <w'x'> ) >= 0,
  !    Skx^2 < 4 * ( 1 - mixt_frac )^2 / ( mixt_frac * ( 2 - mixt_frac ) )
  !    is required; and
  ! when Skx * sgn( <w'x'> ) < 0,
  !    Skx^2 < 4 * mixt_frac^2 / ( 1 - mixt_frac^2 ) is required.
  !
  ! Whenever Skx^2 exceeds these limits, Skx must be limited (preserving its
  ! sign) in order to have any value of F_x that will work in the equation set.
  !
  ! When Skx is found to be within the above limits (or after it has been
  ! limited to fall within its limits), the range of valid values of F_x can be
  ! found according to the following:
  !
  ! When Skx * sgn( <w'x'> ) >= 0:
  !
  !     When  4 * mixt_frac^2 / ( 1 - mixt_frac^2 )  <  Skx^2
  !           <  4 * ( 1 - mixt_frac )^2 / ( mixt_frac * ( 2 - mixt_frac ) ):
  !
  !          Minimum sqrt( F_x ):  2nd root (middle-valued root; also smallest
  !                                positive) of the second equation (sigma_x_2
  !                                based).
  !
  !          Maximum sqrt( F_x ):  Minimum of the largest root of the second
  !                                equation (sigma_x_2 based) and the only* root
  !                                of the first equation (sigma_x_1 based).
  !
  !     When  Skx^2  <=  4 * mixt_frac^2 / ( 1 - mixt_frac^2 ):
  !
  !          Minimum sqrt( F_x ):  2nd root (middle-valued root; also smallest
  !                                positive) of the second equation (sigma_x_2
  !                                based).
  !
  !          Maximum sqrt( F_x ):  Minimum of the largest root of the second
  !                                equation (sigma_x_2 based) and the largest
  !                                root of the first equation (sigma_x_1 based).
  !
  ! When Skx * sgn( <w'x'> ) < 0:
  !
  !     When  4 * ( 1 - mixt_frac )^2 / ( mixt_frac * ( 2 - mixt_frac ) )
  !           <  Skx^2  <  4 * mixt_frac^2 / ( 1 - mixt_frac^2 ):
  !
  !          Minimum sqrt( F_x ):  2nd root (middle-valued root; also smallest
  !                                positive) of the first equation (sigma_x_1
  !                                based).
  !
  !          Maximum sqrt( F_x ):  Minimum of the largest root of the first
  !                                equation (sigma_x_1 based) and the only* root
  !                                of the second equation (sigma_x_2 based).
  !
  !     When  Skx^2
  !           <=  4 * ( 1 - mixt_frac )^2 / ( mixt_frac * ( 2 - mixt_frac ) ):
  !
  !          Minimum sqrt( F_x ):  2nd root (middle-valued root; also smallest
  !                                positive) of the first equation (sigma_x_1
  !                                based).
  !
  !          Maximum sqrt( F_x ):  Minimum of the largest root of the first
  !                                equation (sigma_x_1 based) and the largest
  !                                root of the second equation (sigma_x_2
  !                                based).
  !
  ! Here, "only* root" means the the only root that isn't a complex root.
  !
  ! The value of sqrt( F_x ) is also limited with a minimum of 0 and a maximum
  ! of 1.  The minimum and maximum allowable values of F_x are found by taking
  ! the square of the minimum and maximum allowable values of sqrt( F_x ),
  ! respectively.
  !
  !
  ! Notes:
  !
  ! When F_x = 0 (which can only happen when Skx = 0), mu_x_1 = mu_x_2, and
  ! sigma_x_1 = sigma_x_2 = sqrt( <x'^2> ).  This means that the distribution
  ! becomes a single Gaussian when F_x = 0 (and Skx = 0).
  !
  ! The equations for the PDF component means and standard deviations can also
  ! be written as:
  !
  ! mu_x_1 = <x> + sqrt( F_x * ( ( 1 - mixt_frac ) / mixt_frac ) * <x'^2> )
  !                * sgn( <w'x'> );
  !
  ! mu_x_2 = <x> - sqrt( F_x * ( mixt_frac / ( 1 - mixt_frac ) ) * <x'^2> )
  !                * sgn( <w'x'> );
  !
  ! sigma_x_1 = sqrt( coef_sigma_x_1_sqd * <x'^2> ); and
  !
  ! sigma_x_2 = sqrt( coef_sigma_x_2_sqd * <x'^2> ); where
  !
  ! coef_sigma_x_1_sqd
  ! = ( sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> )
  !     - ( 1 + mixt_frac ) * F_x^1.5 + 3 * mixt_frac * sqrt( F_x ) )
  !   / ( 3 * mixt_frac * sqrt( F_x ) )
  ! = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> )
  !   / ( 3 * mixt_frac * sqrt( F_x ) )
  !   - ( 1 + mixt_frac ) * F_x / ( 3 * mixt_frac )
  !   + 1; and
  !
  ! coef_sigma_x_2_sqd
  ! = ( 1 - F_x ) / ( 1 - mixt_frac )
  !   - mixt_frac / ( 1 - mixt_frac )
  !     * ( sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> )
  !         / ( 3 * mixt_frac * sqrt( F_x ) )
  !         - ( 1 + mixt_frac ) * F_x / ( 3 * mixt_frac )
  !         + 1 )
  ! = ( ( 1 - F_x ) - mixt_frac * coef_sigma_x_1_sqd ) / ( 1 - mixt_frac ).
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
  subroutine calc_responder_params( xm, xp2, Skx, sgn_wpxp,       & ! In
                                    F_x, mixt_frac,               & ! In
                                    mu_x_1, mu_x_2,               & ! Out
                                    sigma_x_1_sqd, sigma_x_2_sqd, & ! Out
                                    coef_sigma_x_1_sqd,           & ! Out
                                    coef_sigma_x_2_sqd            ) ! Out

    ! Description:
    ! Calculates the PDF component means and the PDF component standard
    ! deviations for a responding variable (a variable that is not used to set
    ! the mixture fraction).

    ! References:
    ! Griffin and Larson (2018)
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        three, & ! Variable(s)
        one,   &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      xm,       & ! Mean of x (overall)                             [units vary]
      xp2,      & ! Variance of x (overall)                     [(units vary)^2]
      Skx,      & ! Skewness of x                                            [-]
      sgn_wpxp, & ! Sign of the covariance of w and x (overall)              [-]
      F_x,      & ! Parameter for the spread of the PDF component means of x [-]
      mixt_frac   ! Mixture fraction                                         [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      mu_x_1,        & ! Mean of x (1st PDF component)              [units vary]
      mu_x_2,        & ! Mean of x (2nd PDF component)              [units vary]
      sigma_x_1_sqd, & ! Variance of x (1st PDF component)      [(units vary)^2]
      sigma_x_2_sqd    ! Variance of x (2nd PDF component)      [(units vary)^2]

    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      coef_sigma_x_1_sqd, & ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>      [-]
      coef_sigma_x_2_sqd    ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>      [-]


    where ( F_x > zero )

       ! Calculate the mean of x in the 1st PDF component.
       mu_x_1 = xm + sqrt( F_x * ( ( one - mixt_frac ) / mixt_frac ) * xp2 ) &
                     * sgn_wpxp

       ! Calculate the mean of x in the 2nd PDF component.
       mu_x_2 = xm - ( mixt_frac / ( one - mixt_frac ) ) * ( mu_x_1 - xm )

       ! Calculate the variance of x in the 1st PDF component.
       ! sigma_x_1^2
       ! = ( ( sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> )
       !       - ( 1 + mixt_frac ) * F_x^1.5 + 3 * mixt_frac * sqrt( F_x ) )
       !     / ( 3 * mixt_frac * sqrt( F_x ) ) ) * <x'^2>
       ! = ( sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> )
       !     / ( 3 * mixt_frac * sqrt( F_x ) )
       !     - ( 1 + mixt_frac ) * F_x / ( 3 * mixt_frac )
       !     + 1 ) * <x'^2>
       coef_sigma_x_1_sqd &
       = sqrt( mixt_frac * ( one - mixt_frac ) ) * Skx * sgn_wpxp &
         / ( three * mixt_frac * sqrt( F_x ) ) &
         - ( one + mixt_frac ) * F_x / ( three * mixt_frac ) &
         + one

       sigma_x_1_sqd = coef_sigma_x_1_sqd * xp2

       ! Calculate the variance of x in the 2nd PDF component.
       ! sigma_x_2^2
       ! = ( ( 1 - F_x ) / ( 1 - mixt_frac )
       !     - mixt_frac / ( 1 - mixt_frac )
       !       * ( sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> )
       !           / ( 3 * mixt_frac * sqrt( F_x ) )
       !           - ( 1 + mixt_frac ) * F_x / ( 3 * mixt_frac )
       !           + 1 ) ) * <x'^2>
       ! = ( ( ( 1 - F_x ) - mixt_frac * coef_sigma_x_1_sqd )
       !     / ( 1 - mixt_frac ) ) * <x'^2>
       coef_sigma_x_2_sqd &
       = ( ( one - F_x ) - mixt_frac * coef_sigma_x_1_sqd ) &
         / ( one - mixt_frac )

       sigma_x_2_sqd = coef_sigma_x_2_sqd * xp2

    elsewhere ! F_x = 0

       ! When F_x has a value of 0, the PDF becomes a single Gaussian.  This
       ! only works when Skx = 0.  However, when Skx /= 0, the value of min_F_x
       ! is greater than 0, preventing a problem where F_x = 0 but | Skx | > 0.
       mu_x_1 = xm
       mu_x_2 = xm
       sigma_x_1_sqd = xp2
       sigma_x_2_sqd = xp2
       coef_sigma_x_1_sqd = one
       coef_sigma_x_2_sqd = one

    endwhere ! F_x > 0


    return

  end subroutine calc_responder_params

  !=============================================================================
  subroutine calc_limits_F_x_responder( mixt_frac, Skx, sgn_wpxp,  & ! In
                                        max_Skx2_pos_Skx_sgn_wpxp, & ! In
                                        max_Skx2_neg_Skx_sgn_wpxp, & ! In
                                        min_F_x, max_F_x )           ! Out

    ! Description:
    ! Calculates the minimum and maximum allowable values for F_x for a
    ! responding variable.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        three, & ! Variable(s)
        two,   &
        one,   &
        zero

    use calc_roots, only: &
        cubic_solve    ! Procedure(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      mixt_frac, & ! Mixture fraction                 [-]
      Skx,       & ! Skewness of x                    [-]
      sgn_wpxp     ! Sign of covariance of w and x    [-]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      max_Skx2_pos_Skx_sgn_wpxp, & ! Maximum Skx^2 when Skx*sgn(<w'x'>) >= 0 [-]
      max_Skx2_neg_Skx_sgn_wpxp    ! Maximum Skx^2 when Skx*sgn(<w'x'>) < 0  [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      min_F_x, & ! Minimum allowable value of F_x    [-]
      max_F_x    ! Maximum allowable value of F_x    [-]

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      coef_A_1, & ! Coef. A in Ax^3 + Bx^2 + Cx + D = 0 (1st PDF comp. lim.) [-]
      coef_B_1, & ! Coef. B in Ax^3 + Bx^2 + Cx + D = 0 (1st PDF comp. lim.) [-]
      coef_C_1, & ! Coef. C in Ax^3 + Bx^2 + Cx + D = 0 (1st PDF comp. lim.) [-]
      coef_D_1, & ! Coef. D in Ax^3 + Bx^2 + Cx + D = 0 (1st PDF comp. lim.) [-]
      coef_A_2, & ! Coef. A in Ax^3 + Bx^2 + Cx + D = 0 (2nd PDF comp. lim.) [-]
      coef_B_2, & ! Coef. B in Ax^3 + Bx^2 + Cx + D = 0 (2nd PDF comp. lim.) [-]
      coef_C_2, & ! Coef. C in Ax^3 + Bx^2 + Cx + D = 0 (2nd PDF comp. lim.) [-]
      coef_D_2    ! Coef. D in Ax^3 + Bx^2 + Cx + D = 0 (2nd PDF comp. lim.) [-]

    complex( kind = core_rknd ), dimension(gr%nz,3) :: &
      sqrt_F_x_roots_1, & ! Roots of sqrt(F_x) for the sigma_x_1 term    [-]
      sqrt_F_x_roots_2    ! Roots of sqrt(F_x) for the sigma_x_2 term    [-]

    real( kind = core_rknd ), dimension(gr%nz,3) :: &
      sqrt_F_x_roots_1_sorted, & ! Sorted roots of sqrt(F_x): sigma_x_1 term [-]
      sqrt_F_x_roots_2_sorted    ! Sorted roots of sqrt(F_x): sigma_x_2 term [-]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      min_sqrt_F_x, & ! Minimum allowable value of sqrt(F_x)    [-]
      max_sqrt_F_x    ! Maximum allowable value of sqrt(F_x)    [-]

 
    ! Set up the coefficients in the equation for the limit of sqrt(F_x) based
    ! on the 1st PDF component standard deviation (sigma_x_1) being greater than
    ! or equal to 0.  This equation has the form:
    ! A * sqrt(F_x)^3 + B * sqrt(F_x)^2 + C * sqrt(F_x) + D = 0.
    coef_A_1 = -( one + mixt_frac )
    coef_B_1 = zero
    coef_C_1 = three * mixt_frac
    coef_D_1 = sqrt( mixt_frac * ( one - mixt_frac ) ) * Skx * sgn_wpxp

    ! Solve for the roots (values of sqrt(F_x)) that satisfy the above equation.
    sqrt_F_x_roots_1 &
    = cubic_solve( gr%nz, coef_A_1, coef_B_1, coef_C_1, coef_D_1 )

    ! Sort the values of the roots (values of sqrt(F_x)) from smallest to
    ! largest.  Ignore any complex component of the roots.  The code below that
    ! uses sqrt_F_x_roots_1_sorted already factors the appropriate roots to use
    ! into account.
    sqrt_F_x_roots_1_sorted &
    = sort_roots( real( sqrt_F_x_roots_1, kind = core_rknd ) )

    ! Set up the coefficients in the equation for the limit of sqrt(F_x) based
    ! on the 2nd PDF component standard deviation (sigma_x_2) being greater than
    ! or equal to 0.  This equation has the form:
    ! A * sqrt(F_x)^3 + B * sqrt(F_x)^2 + C * sqrt(F_x) + D = 0.
    coef_A_2 = -( two - mixt_frac )
    coef_B_2 = zero
    coef_C_2 = three * ( one - mixt_frac )
    coef_D_2 = -sqrt( mixt_frac * ( one - mixt_frac ) ) * Skx * sgn_wpxp

    ! Solve for the roots (values of sqrt(F_x)) that satisfy the above equation.
    sqrt_F_x_roots_2 &
    = cubic_solve( gr%nz, coef_A_2, coef_B_2, coef_C_2, coef_D_2 )

    ! Sort the values of the roots (values of sqrt(F_x)) from smallest to
    ! largest.  Ignore any complex component of the roots.  The code below that
    ! uses sqrt_F_x_roots_2_sorted already factors the appropriate roots to use
    ! into account.
    sqrt_F_x_roots_2_sorted &
    = sort_roots( real( sqrt_F_x_roots_2, kind = core_rknd ) )


    ! Find the minimum and maximum allowable values of sqrt(F_x) based on Skx
    ! and sgn( <w'x'> ).
    where ( Skx * sgn_wpxp >= zero )

       where ( Skx**2 > max_Skx2_neg_Skx_sgn_wpxp )

          min_sqrt_F_x = sqrt_F_x_roots_2_sorted(:,2)
          max_sqrt_F_x = min( real( sqrt_F_x_roots_1(:,1), kind = core_rknd ), &
                              sqrt_F_x_roots_2_sorted(:,3) )

       elsewhere ! Skx^2 <= max_Skx2_neg_Skx_sgn_wpxp

          min_sqrt_F_x = sqrt_F_x_roots_2_sorted(:,2)
          max_sqrt_F_x = min( sqrt_F_x_roots_1_sorted(:,3), &
                              sqrt_F_x_roots_2_sorted(:,3) )

       endwhere ! Skx**2 > max_Skx2_neg_Skx_sgn_wpxp

    elsewhere ! Skx * sgn( <w'x'> ) < 0 

       where ( Skx**2 > max_Skx2_pos_Skx_sgn_wpxp )

          min_sqrt_F_x = sqrt_F_x_roots_1_sorted(:,2)
          max_sqrt_F_x = min( real( sqrt_F_x_roots_2(:,1), kind = core_rknd ), &
                              sqrt_F_x_roots_1_sorted(:,3) )

       elsewhere ! Skx^2 <= max_Skx2_pos_Skx_sgn_wpxp

          min_sqrt_F_x = sqrt_F_x_roots_1_sorted(:,2)
          max_sqrt_F_x = min( sqrt_F_x_roots_1_sorted(:,3), &
                              sqrt_F_x_roots_2_sorted(:,3) )

       endwhere ! Skx**2 > max_Skx2_pos_Skx_sgn_wpxp

    endwhere ! Skx * sgn( <w'x'> ) >= 0


    ! The minimum and maximum are also limited by 0 and 1, respectively.
    min_sqrt_F_x = max( min_sqrt_F_x, zero )
    max_sqrt_F_x = min( max_sqrt_F_x, one )

    ! The minimum and maximum allowable values for F_x are the squares of the
    ! minimum and maximum allowable values for sqrt(F_x).
    min_F_x = min_sqrt_F_x**2
    max_F_x = max_sqrt_F_x**2


    return

  end subroutine calc_limits_F_x_responder

  !=============================================================================
  function sort_roots( roots ) &
  result ( roots_sorted )

    ! Description:
    ! Sorts roots from smallest to largest.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variable
    real( kind = core_rknd ), dimension(gr%nz,3), intent(in) :: &
      roots    ! Roots    [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz,3) :: &
      roots_sorted    ! Roots sorted from smallest to largest    [-]


    where ( roots(:,1) <= roots(:,2) .and. roots(:,1) <= roots(:,3) )

       ! The value of roots(1) is the smallest root.
       roots_sorted(:,1) = roots(:,1)

       where ( roots(:,2) <= roots(:,3) )

          ! The value of roots(2) is the middle-valued root and the value of
          ! roots(3) is the largest root.
          roots_sorted(:,2) = roots(:,2)
          roots_sorted(:,3) = roots(:,3)

       elsewhere ! roots(3) < roots(2)

          ! The value of roots(3) is the middle-valued root and the value of
          ! roots(2) is the largest root.
          roots_sorted(:,2) = roots(:,3)
          roots_sorted(:,3) = roots(:,2)

       endwhere ! roots(2) <= roots(3)

    elsewhere ( roots(:,2) < roots(:,1) .and. roots(:,2) <= roots(:,3) )

       ! The value of roots(2) is the smallest root.
       roots_sorted(:,1) = roots(:,2)

       where ( roots(:,1) <= roots(:,3) )

          ! The value of roots(1) is the middle-valued root and the value of
          ! roots(3) is the largest root.
          roots_sorted(:,2) = roots(:,1)
          roots_sorted(:,3) = roots(:,3)

       elsewhere ! roots(3) < roots(1)

          ! The value of roots(3) is the middle-valued root and the value of
          ! roots(1) is the largest root.
          roots_sorted(:,2) = roots(:,3)
          roots_sorted(:,3) = roots(:,1)

       endwhere ! roots(1) <= roots(3)

    elsewhere ! roots(3) < roots(1) .and. roots(3) < roots(2)

       ! The value of roots(3) is the smallest root.
       roots_sorted(:,1) = roots(:,3)

       where ( roots(:,1) <= roots(:,2) )

          ! The value of roots(1) is the middle-valued root and the value of
          ! roots(2) is the largest root.
          roots_sorted(:,2) = roots(:,1)
          roots_sorted(:,3) = roots(:,2)

       elsewhere ! roots(2) < roots(1)

          ! The value of roots(2) is the middle-valued root and the value of
          ! roots(1) is the largest root.
          roots_sorted(:,2) = roots(:,2)
          roots_sorted(:,3) = roots(:,1)

       endwhere ! roots(1) <= roots(2)

    endwhere ! roots(1) <= roots(2) .and. roots(1) <= roots(3)


    return

  end function sort_roots

  !=============================================================================
  function calc_coef_wp4_implicit( mixt_frac, F_w, &
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
    ! When w is the setting variable, coef_sigma_w_1_sqd and coef_sigma_w_2_sqd
    ! are given by:
    !
    ! coef_sigma_w_1_sqd = ( ( zeta_w + 1 ) * ( 1 - F_w ) )
    !                      / ( ( zeta_w + 2 ) * mixt_frac ); and
    !
    ! coef_sigma_w_2_sqd = ( 1 - F_w ) / ( ( zeta_w + 2 ) * ( 1 - mixt_frac ) ).
    !
    ! When w is a responding variable, coef_sigma_w_1_sqd and coef_sigma_w_2_sqd
    ! are given by:
    !
    ! coef_sigma_w_1_sqd = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skw
    !                      / ( 3 * mixt_frac * sqrt( F_w ) )
    !                      - ( 1 + mixt_frac ) * F_w / ( 3 * mixt_frac )
    !                      + 1; and
    !
    ! coef_sigma_w_2_sqd = ( ( 1 - F_w ) - mixt_frac * coef_sigma_w_1_sqd )
    !                      / ( 1 - mixt_frac ).
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

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        six,   & ! Variable(s)
        three, &
        one

    use clubb_precision, only: &
        core_rknd    ! Procedure(s)

    implicit none

    ! Input Variables
    real ( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      mixt_frac,          & ! Mixture fraction                               [-]
      F_w,                & ! Parameter: spread of the PDF comp. means of w  [-]
      coef_sigma_w_1_sqd, & ! sigma_w_1^2 = coef_sigma_w_1_sqd * <w'^2>      [-]
      coef_sigma_w_2_sqd    ! sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>      [-]

    ! Return Variable
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wp4_implicit    ! Coef.: <w'^4> = coef_wp4_implicit * <w'^2>^2    [-]


    ! Calculate coef_wp4_implicit.
    coef_wp4_implicit = three * mixt_frac * coef_sigma_w_1_sqd**2 &
                        + six * F_w * ( one - mixt_frac ) * coef_sigma_w_1_sqd &
                        + F_w**2 * ( one - mixt_frac )**2 / mixt_frac &
                        + three * ( one - mixt_frac ) * coef_sigma_w_2_sqd**2 &
                        + six * F_w * mixt_frac * coef_sigma_w_2_sqd &
                        + F_w**2 * mixt_frac**2 / ( one - mixt_frac )


    return

  end function calc_coef_wp4_implicit

  !=============================================================================
  function calc_coef_wpxp2_implicit( wp2, xp2, wpxp, sgn_wpxp, &
                                     mixt_frac, F_w, F_x, &
                                     coef_sigma_w_1_sqd, &
                                     coef_sigma_w_2_sqd, &
                                     coef_sigma_x_1_sqd, &
                                     coef_sigma_x_2_sqd  ) &
  result( coef_wpxp2_implicit )

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
    ! mu_x_1 - <x> = sqrt(F_x) * sqrt( ( 1 - mixt_frac ) / mixt_frac )
    !                * sqrt( <x'^2> ) * sgn( <w'x'> );
    !
    ! mu_x_2 - <x> = - sqrt(F_x) * sqrt( mixt_frac / ( 1 - mixt_frac ) )
    !                  * sqrt( <x'^2> ) * sgn( <w'x'> );
    !
    ! sigma_w_1 = sqrt( coef_sigma_w_1_sqd * <w'^2> );
    !
    ! sigma_w_2 = sqrt( coef_sigma_w_2_sqd * <w'^2> );
    !
    ! sigma_x_1 = sqrt( coef_sigma_x_1_sqd * <x'^2> ); and
    !
    ! sigma_x_2 = sqrt( coef_sigma_x_2_sqd * <x'^2> ).
    !
    ! Either w can be the setting variable and x can be a responding variable,
    ! x can be the setting variable and w can be a responding variable, or both
    ! w and x can be responding variables.
    !
    ! When w is the setting variable, coef_sigma_w_1_sqd and coef_sigma_w_2_sqd
    ! are given by:
    !
    ! coef_sigma_w_1_sqd = ( ( zeta_w + 1 ) * ( 1 - F_w ) )
    !                      / ( ( zeta_w + 2 ) * mixt_frac ); and
    !
    ! coef_sigma_w_2_sqd = ( 1 - F_w ) / ( ( zeta_w + 2 ) * ( 1 - mixt_frac ) ).
    !
    ! When w is a responding variable, coef_sigma_w_1_sqd and coef_sigma_w_2_sqd
    ! are given by:
    !
    ! coef_sigma_w_1_sqd = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skw
    !                      / ( 3 * mixt_frac * sqrt( F_w ) )
    !                      - ( 1 + mixt_frac ) * F_w / ( 3 * mixt_frac )
    !                      + 1; and
    !
    ! coef_sigma_w_2_sqd = ( ( 1 - F_w ) - mixt_frac * coef_sigma_w_1_sqd )
    !                      / ( 1 - mixt_frac ).
    !
    ! When x is the setting variable, coef_sigma_x_1_sqd and coef_sigma_x_2_sqd
    ! are given by:
    !
    ! coef_sigma_x_1_sqd = ( ( zeta_x + 1 ) * ( 1 - F_x ) )
    !                      / ( ( zeta_x + 2 ) * mixt_frac ); and
    !
    ! coef_sigma_x_2_sqd = ( 1 - F_x ) / ( ( zeta_x + 2 ) * ( 1 - mixt_frac ) ).
    !
    ! When x is a responding variable, coef_sigma_x_1_sqd and coef_sigma_x_2_sqd
    ! are given by:
    !
    ! coef_sigma_x_1_sqd = sqrt( mixt_frac * ( 1 - mixt_frac ) )
    !                      * Skx * sgn( <w'x'> )
    !                      / ( 3 * mixt_frac * sqrt( F_x ) )
    !                      - ( 1 + mixt_frac ) * F_x / ( 3 * mixt_frac )
    !                      + 1; and
    !
    ! coef_sigma_x_2_sqd = ( ( 1 - F_x ) - mixt_frac * coef_sigma_x_1_sqd )
    !                      / ( 1 - mixt_frac ).
    !
    ! Additionally:
    !
    ! corr_w_x_1 = corr_w_x_2
    ! = ( <w'x'> - mixt_frac * ( mu_w_1 - <w> ) * ( mu_x_1 - <x> )
    !            - ( 1 - mixt_frac ) * ( mu_w_2 - <w> ) * ( mu_x_2 - <x> ) )
    !   / ( mixt_frac * sigma_w_1 * sigma_x_1
    !       + ( 1 - mixt_frac ) * sigma_w_2 * sigma_x_2 );
    !
    ! where -1 <= corr_w_x_1 = corr_w_x_2 <= 1.  This equation can be rewritten
    ! as:
    !
    ! corr_w_x_1 = corr_w_x_2
    ! = ( <w'x'>
    !     - sqrt( F_w ) * sqrt( F_x ) * sgn( <w'x'> )
    !       * sqrt( <w'^2> ) * sqrt( <x'^2 > ) )
    !   / ( ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !         + ( 1 - mixt_frac )
    !           * sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) )
    !       * sqrt( <w'^2> ) * sqrt( <x'^2> ) ).
    !
    ! The equation for <w'x'^2> becomes:
    !
    ! <w'x'^2>
    ! = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( <w'^2> )
    !   * ( sqrt( F_w ) * F_x
    !       * ( ( 1 - mixt_frac ) / mixt_frac - mixt_frac / ( 1 - mixt_frac ) )
    !       + sqrt( F_w ) * ( coef_sigma_x_1_sqd - coef_sigma_x_2_sqd )
    !       + ( 2 * sqrt( F_x ) * sgn( <w'x'> ) * <w'x'>
    !           * ( sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !               - sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) ) )
    !         / ( ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !               + ( 1 - mixt_frac )
    !                 * sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) )
    !             * sqrt( <w'^2> ) * sqrt( <x'^2> ) )
    !       - ( 2 * sqrt( F_w ) * F_x 
    !           * ( sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !               - sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) ) )
    !         / ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !             + ( 1 - mixt_frac )
    !               * sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) ) )
    !   * <x'^2>
    !
    ! This equation is of the form:
    !
    ! <w'x'^2> = coef_wpxp2_implicit * <x'^2>;
    !
    ! where:
    !
    ! coef_wpxp2_implicit
    ! = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( <w'^2> )
    !   * ( sqrt( F_w ) * F_x
    !       * ( ( 1 - mixt_frac ) / mixt_frac - mixt_frac / ( 1 - mixt_frac ) )
    !       + sqrt( F_w ) * ( coef_sigma_x_1_sqd - coef_sigma_x_2_sqd )
    !       + ( 2 * sqrt( F_x ) * sgn( <w'x'> ) * <w'x'>
    !           * ( sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !               - sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) ) )
    !         / ( ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !               + ( 1 - mixt_frac )
    !                 * sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) )
    !             * sqrt( <w'^2> ) * sqrt( <x'^2> ) )
    !       - ( 2 * sqrt( F_w ) * F_x 
    !           * ( sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !               - sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) ) )
    !         / ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !             + ( 1 - mixt_frac )
    !               * sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) ) ).
    !
    ! In the special case that coef_sigma_w_1_sqd * coef_sigma_x_1_sqd = 0 and
    ! coef_sigma_w_2_sqd * coef_sigma_x_2_sqd = 0, the above equation is
    ! undefined.  However, the equation for this special case can be derived by
    ! taking the original equation for <w'x'^2> and setting both
    ! sigma_w_1 * sigma_x_1 = 0 and sigma_w_2 * sigma_x_2 = 0.  The equation
    ! becomes:
    !
    ! <w'x'^2>
    ! = mixt_frac * ( ( mu_w_1 - <w> ) * ( ( mu_x_1 - <x> )^2 + sigma_x_1^2 )
    !   + ( 1 - mixt_frac ) * ( ( mu_w_2 - <w> )
    !                           * ( ( mu_x_2 - <x> )^2 + sigma_x_2^2 );
    !
    ! and making the same substitutions as before, it can be rewritten as:
    !
    ! <w'x'^2>
    ! = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( <w'^2> ) * sqrt( F_w )
    !   * ( F_x * ( ( 1 - mixt_frac ) / mixt_frac
    !               - mixt_frac / ( 1 - mixt_frac ) )
    !       + ( coef_sigma_x_1_sqd - coef_sigma_x_2_sqd ) ) * <x'^2>.
    !
    ! The coefficient in this special case is:
    !
    ! coef_wpxp2_implicit
    ! = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( <w'^2> ) * sqrt( F_w )
    !   * ( F_x * ( ( 1 - mixt_frac ) / mixt_frac
    !               - mixt_frac / ( 1 - mixt_frac ) )
    !       + ( coef_sigma_x_1_sqd - coef_sigma_x_2_sqd ) ).

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        two,  & ! Variable(s)
        one,  &
        zero

    use clubb_precision, only: &
        core_rknd    ! Procedure(s)

    implicit none

    ! Input Variables
    real ( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      wp2,                & ! Variance of w (overall)                  [m^2/s^2]
      xp2,                & ! Variance of x (overall)           [(units vary)^2]
      wpxp,               & ! Covariance of w and x           [m/s (units vary)]
      sgn_wpxp,           & ! Sign of the covariance of w and x              [-]
      mixt_frac,          & ! Mixture fraction                               [-]
      F_w,                & ! Parameter: spread of the PDF comp. means of w  [-]
      F_x,                & ! Parameter: spread of the PDF comp. means of x  [-]
      coef_sigma_w_1_sqd, & ! sigma_w_1^2 = coef_sigma_w_1_sqd * <w'^2>      [-]
      coef_sigma_w_2_sqd, & ! sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>      [-]
      coef_sigma_x_1_sqd, & ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>      [-]
      coef_sigma_x_2_sqd    ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>      [-]

    ! Return Variable
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wpxp2_implicit ! Coef.: <w'x'^2> = coef_wpxp2_implicit * <x'^2> [m/s]

    ! Local Variable
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coefs_factor    ! Factor involving coef_sigma_... coefficients         [-]


    ! Calculate coef_wpxp2_implicit.
    where ( ( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd > zero &
              .or. coef_sigma_w_2_sqd * coef_sigma_x_2_sqd > zero ) &
            .and. ( wp2 * xp2 > zero ) )

       coefs_factor &
       = ( sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd ) &
           - sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) ) &
         / ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd ) &
             + ( one - mixt_frac ) &
               * sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) )

       coef_wpxp2_implicit &
       = sqrt( mixt_frac * ( one - mixt_frac ) ) * sqrt( wp2 ) &
         * ( sqrt( F_w ) * F_x &
             * ( ( one - mixt_frac ) / mixt_frac &
                 - mixt_frac / ( one - mixt_frac ) ) &
             + sqrt( F_w ) * ( coef_sigma_x_1_sqd - coef_sigma_x_2_sqd ) &
             + two * sqrt( F_x ) * coefs_factor * sgn_wpxp * wpxp &
               / ( sqrt( wp2 ) * sqrt( xp2 ) ) &
             - two * sqrt( F_w ) * F_x * coefs_factor )

    elsewhere ! ( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd = 0
              !   and coef_sigma_w_2_sqd * coef_sigma_x_2_sqd = 0 )
              ! or wp2 * xp2 = 0

       coef_wpxp2_implicit &
       = sqrt( mixt_frac * ( one - mixt_frac ) ) * sqrt( wp2 ) * sqrt( F_w ) &
         * ( F_x * ( ( one - mixt_frac ) / mixt_frac &
                     - mixt_frac / ( one - mixt_frac ) ) &
             + ( coef_sigma_x_1_sqd - coef_sigma_x_2_sqd ) )

    endwhere


    return

  end function calc_coef_wpxp2_implicit

  !=============================================================================
  subroutine calc_coefs_wp2xp_semiimpl( wp2, xp2, sgn_wpxp,  & ! In
                                        mixt_frac, F_w, F_x, & ! In
                                        coef_sigma_w_1_sqd,  & ! In
                                        coef_sigma_w_2_sqd,  & ! In
                                        coef_sigma_x_1_sqd,  & ! In
                                        coef_sigma_x_2_sqd,  & ! In
                                        coef_wp2xp_implicit, & ! Out
                                        term_wp2xp_explicit  ) ! Out

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
    ! mu_x_1 - <x> = sqrt(F_x) * sqrt( ( 1 - mixt_frac ) / mixt_frac )
    !                * sqrt( <x'^2> ) * sgn( <w'x'> );
    !
    ! mu_x_2 - <x> = - sqrt(F_x) * sqrt( mixt_frac / ( 1 - mixt_frac ) )
    !                  * sqrt( <x'^2> ) * sgn( <w'x'> );
    !
    ! sigma_w_1 = sqrt( coef_sigma_w_1_sqd * <w'^2> );
    !
    ! sigma_w_2 = sqrt( coef_sigma_w_2_sqd * <w'^2> );
    !
    ! sigma_x_1 = sqrt( coef_sigma_x_1_sqd * <x'^2> ); and
    !
    ! sigma_x_2 = sqrt( coef_sigma_x_2_sqd * <x'^2> ).
    !
    ! Either w can be the setting variable and x can be a responding variable,
    ! x can be the setting variable and w can be a responding variable, or both
    ! w and x can be responding variables.
    !
    ! When w is the setting variable, coef_sigma_w_1_sqd and coef_sigma_w_2_sqd
    ! are given by:
    !
    ! coef_sigma_w_1_sqd = ( ( zeta_w + 1 ) * ( 1 - F_w ) )
    !                      / ( ( zeta_w + 2 ) * mixt_frac ); and
    !
    ! coef_sigma_w_2_sqd = ( 1 - F_w ) / ( ( zeta_w + 2 ) * ( 1 - mixt_frac ) ).
    !
    ! When w is a responding variable, coef_sigma_w_1_sqd and coef_sigma_w_2_sqd
    ! are given by:
    !
    ! coef_sigma_w_1_sqd = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skw
    !                      / ( 3 * mixt_frac * sqrt( F_w ) )
    !                      - ( 1 + mixt_frac ) * F_w / ( 3 * mixt_frac )
    !                      + 1; and
    !
    ! coef_sigma_w_2_sqd = ( ( 1 - F_w ) - mixt_frac * coef_sigma_w_1_sqd )
    !                      / ( 1 - mixt_frac ).
    !
    ! When x is the setting variable, coef_sigma_x_1_sqd and coef_sigma_x_2_sqd
    ! are given by:
    !
    ! coef_sigma_x_1_sqd = ( ( zeta_x + 1 ) * ( 1 - F_x ) )
    !                      / ( ( zeta_x + 2 ) * mixt_frac ); and
    !
    ! coef_sigma_x_2_sqd = ( 1 - F_x ) / ( ( zeta_x + 2 ) * ( 1 - mixt_frac ) ).
    !
    ! When x is a responding variable, coef_sigma_x_1_sqd and coef_sigma_x_2_sqd
    ! are given by:
    !
    ! coef_sigma_x_1_sqd = sqrt( mixt_frac * ( 1 - mixt_frac ) )
    !                      * Skx * sgn( <w'x'> )
    !                      / ( 3 * mixt_frac * sqrt( F_x ) )
    !                      - ( 1 + mixt_frac ) * F_x / ( 3 * mixt_frac )
    !                      + 1; and
    !
    ! coef_sigma_x_2_sqd = ( ( 1 - F_x ) - mixt_frac * coef_sigma_x_1_sqd )
    !                      / ( 1 - mixt_frac ).
    !
    ! Additionally:
    !
    ! corr_w_x_1 = corr_w_x_2
    ! = ( <w'x'> - mixt_frac * ( mu_w_1 - <w> ) * ( mu_x_1 - <x> )
    !            - ( 1 - mixt_frac ) * ( mu_w_2 - <w> ) * ( mu_x_2 - <x> ) )
    !   / ( mixt_frac * sigma_w_1 * sigma_x_1
    !       + ( 1 - mixt_frac ) * sigma_w_2 * sigma_x_2 );
    !
    ! where -1 <= corr_w_x_1 = corr_w_x_2 <= 1.  This equation can be rewritten
    ! as:
    !
    ! corr_w_x_1 = corr_w_x_2
    ! = ( <w'x'>
    !     - sqrt( F_w ) * sqrt( F_x ) * sgn( <w'x'> )
    !       * sqrt( <w'^2> ) * sqrt( <x'^2 > ) )
    !   / ( ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !         + ( 1 - mixt_frac )
    !           * sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) )
    !       * sqrt( <w'^2> ) * sqrt( <x'^2> ) ).
    !
    ! The equation for <w'^2 x'> becomes:
    !
    ! <w'^2 x'>
    ! = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( F_x )
    !   * sqrt( <x'^2> ) * <w'^2> * sgn( <w'x'> )
    !   * ( F_w * ( ( 1 - mixt_frac ) / mixt_frac
    !               - mixt_frac / ( 1 - mixt_frac ) )
    !       + ( coef_sigma_w_1_sqd - coef_sigma_w_2_sqd )
    !       - ( 2 * F_w * ( sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !                       - sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) )
    !         / ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !             + ( 1 - mixt_frac )
    !               * sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) ) ) )
    !   + ( sqrt( mixt_frac * ( 1 - mixt_frac ) )
    !       * 2 * sqrt( F_w ) * sqrt( <w'^2> )
    !       * ( sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !           - sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) )
    !         / ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !             + ( 1 - mixt_frac )
    !               * sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) ) )
    !     * <w'x'>
    !
    ! This equation is of the form:
    !
    ! <w'^2 x'> = coef_wp2xp_implicit * <w'x'> + term_wp2xp_explicit;
    !
    ! where:
    !
    ! coef_wp2xp_implicit
    ! = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * 2 * sqrt( F_w ) * sqrt( <w'^2> )
    !   * ( sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !       - sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) )
    !     / ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !         + ( 1 - mixt_frac )
    !           * sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) ); and
    !
    ! term_wp2xp_explicit
    ! = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( F_x )
    !   * sqrt( <x'^2> ) * <w'^2> * sgn( <w'x'> )
    !   * ( F_w * ( ( 1 - mixt_frac ) / mixt_frac
    !               - mixt_frac / ( 1 - mixt_frac ) )
    !       + ( coef_sigma_w_1_sqd - coef_sigma_w_2_sqd )
    !       - ( 2 * F_w * ( sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !                       - sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) )
    !         / ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !             + ( 1 - mixt_frac )
    !               * sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) ) ) ).
    !
    ! In the special case that coef_sigma_w_1_sqd * coef_sigma_x_1_sqd = 0 and
    ! coef_sigma_w_2_sqd * coef_sigma_x_2_sqd = 0, the above equation is
    ! undefined.  However, the equation for this special case can be derived by
    ! taking the original equation for <w'^2 x'> and setting both
    ! sigma_w_1 * sigma_x_1 = 0 and sigma_w_2 * sigma_x_2 = 0.  The equation
    ! becomes:
    !
    ! <w'^2 x'>
    ! = mixt_frac * ( ( mu_x_1 - <x> ) * ( ( mu_w_1 - <w> )^2 + sigma_w_1^2 )
    !   + ( 1 - mixt_frac ) * ( ( mu_x_2 - <x> )
    !                           * ( ( mu_w_2 - <w> )^2 + sigma_w_2^2 );
    !
    ! and making the same substitutions as before, it can be rewritten as:
    !
    ! <w'^2 x'>
    ! = sqrt( mixt_frac * ( 1 - mixt_frac ) )
    !   * sqrt( F_x ) * sqrt( <x'^2> ) * <w'^2> * sgn( <w'x'> )
    !   * ( F_w * ( ( 1 - mixt_frac ) / mixt_frac
    !               - mixt_frac / ( 1 - mixt_frac ) )
    !       + ( coef_sigma_w_1_sqd - coef_sigma_w_2_sqd ) ).
    !
    ! Likewise, the equation for <w'x'> in this special case becomes:
    !
    ! <w'x'> = mixt_frac * ( mu_w_1 - <w> ) * ( mu_x_1 - <x> )
    !          + ( 1 - mixt_frac ) * ( mu_w_2 - <w> ) * ( mu_x_2 - <x> );
    !
    ! and making the same substitutions as before, it can be rewritten as:
    !
    ! <w'x'> = sqrt( F_w ) * sqrt( F_x )
    !          * sqrt( <w'^2> ) * sqrt( <x'^2> ) * sgn( <w'x'> ).
    !
    ! The equation for <w'^2 x'> can be rewritten as:
    !
    ! <w'^2 x'>
    ! = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( F_w ) * sqrt( <w'^2> )
    !   * ( ( 1 - mixt_frac ) / mixt_frac - mixt_frac / ( 1 - mixt_frac ) )
    !   * <w'x'>
    !   + sqrt( mixt_frac * ( 1 - mixt_frac ) )
    !     * sqrt( F_x ) * sqrt( <x'^2> ) * <w'^2> * sgn( <w'x'> )
    !     * ( coef_sigma_w_1_sqd - coef_sigma_w_2_sqd ).
    !
    ! The coefficients in this special case are:
    !
    ! coef_wp2xp_implicit
    ! = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( F_w ) * sqrt( <w'^2> )
    !   * ( ( 1 - mixt_frac ) / mixt_frac - mixt_frac / ( 1 - mixt_frac ) ); and
    !
    ! term_wp2xp_explicit
    ! = sqrt( mixt_frac * ( 1 - mixt_frac ) )
    !   * sqrt( F_x ) * sqrt( <x'^2> ) * <w'^2> * sgn( <w'x'> )
    !   * ( coef_sigma_w_1_sqd - coef_sigma_w_2_sqd ).

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        two,  & ! Variable(s)
        one,  &
        zero

    use clubb_precision, only: &
        core_rknd    ! Procedure(s)

    implicit none

    ! Input Variables
    real ( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      wp2,                & ! Variance of w (overall)                  [m^2/s^2]
      xp2,                & ! Variance of x (overall)           [(units vary)^2]
      sgn_wpxp,           & ! Sign of the covariance of w and x              [-]
      mixt_frac,          & ! Mixture fraction                               [-]
      F_w,                & ! Parameter: spread of the PDF comp. means of w  [-]
      F_x,                & ! Parameter: spread of the PDF comp. means of x  [-]
      coef_sigma_w_1_sqd, & ! sigma_w_1^2 = coef_sigma_w_1_sqd * <w'^2>      [-]
      coef_sigma_w_2_sqd, & ! sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>      [-]
      coef_sigma_x_1_sqd, & ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>      [-]
      coef_sigma_x_2_sqd    ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>      [-]

    ! Output Variables
    ! Coefs.: <w'^2 x'> = coef_wp2xp_implicit * <w'x'> + term_wp2xp_explicit
    real ( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      coef_wp2xp_implicit, & ! Coefficient that is multiplied by <w'x'>    [m/s]
      term_wp2xp_explicit    ! Term that is on the RHS    [m^2/s^2 (units vary)]

    ! Local Variable
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coefs_factor    ! Factor involving coef_sigma_... coefficients         [-]


    ! Calculate coef_wp2xp_implicit and term_wp2xp_explicit.
    where ( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd > zero &
            .or. coef_sigma_w_2_sqd * coef_sigma_x_2_sqd > zero )

       coefs_factor &
       = ( sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd ) &
           - sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) ) &
         / ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd ) &
             + ( one - mixt_frac ) &
               * sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) )

       coef_wp2xp_implicit &
       = sqrt( mixt_frac * ( one - mixt_frac ) ) &
         * two * sqrt( F_w ) * sqrt( wp2 ) * coefs_factor

       term_wp2xp_explicit &
       = sqrt( mixt_frac * ( one - mixt_frac ) ) &
         * sqrt( F_x ) * sqrt( xp2 ) * wp2 * sgn_wpxp &
         * ( F_w * ( ( one - mixt_frac ) / mixt_frac &
                     - mixt_frac / ( one - mixt_frac ) ) &
             + ( coef_sigma_w_1_sqd - coef_sigma_w_2_sqd ) &
             - two * F_w * coefs_factor )

    elsewhere ! coef_sigma_w_1_sqd * coef_sigma_x_1_sqd = 0
              ! and coef_sigma_w_2_sqd * coef_sigma_x_2_sqd = 0

       coef_wp2xp_implicit &
       = sqrt( mixt_frac * ( one - mixt_frac ) ) &
         * sqrt( F_w ) * sqrt( wp2 ) &
         * ( ( one - mixt_frac ) / mixt_frac &
             - mixt_frac / ( one - mixt_frac ) )

       term_wp2xp_explicit &
       = sqrt( mixt_frac * ( one - mixt_frac ) ) &
         * sqrt( F_x ) * sqrt( xp2 ) * wp2 * sgn_wpxp &
         * ( coef_sigma_w_1_sqd - coef_sigma_w_2_sqd )

    endwhere


    return

  end subroutine calc_coefs_wp2xp_semiimpl

  !=============================================================================
  subroutine calc_coefs_wpxpyp_semiimpl( wp2, xp2, yp2, wpxp,      & ! In
                                         wpyp, sgn_wpxp, sgn_wpyp, & ! In
                                         mixt_frac, F_w, F_x, F_y, & ! In
                                         coef_sigma_w_1_sqd,       & ! In
                                         coef_sigma_w_2_sqd,       & ! In
                                         coef_sigma_x_1_sqd,       & ! In
                                         coef_sigma_x_2_sqd,       & ! In
                                         coef_sigma_y_1_sqd,       & ! In
                                         coef_sigma_y_2_sqd,       & ! In
                                         coef_wpxpyp_implicit,     & ! Out
                                         term_wpxpyp_explicit      ) ! Out

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
    ! mu_x_1 - <x> = sqrt(F_x) * sqrt( ( 1 - mixt_frac ) / mixt_frac )
    !                * sqrt( <x'^2> ) * sgn( <w'x'> );
    !
    ! mu_x_2 - <x> = - sqrt(F_x) * sqrt( mixt_frac / ( 1 - mixt_frac ) )
    !                  * sqrt( <x'^2> ) * sgn( <w'x'> );
    !
    ! mu_y_1 - <y> = sqrt(F_y) * sqrt( ( 1 - mixt_frac ) / mixt_frac )
    !                * sqrt( <y'^2> ) * sgn( <w'y'> );
    !
    ! mu_y_2 - <y> = - sqrt(F_y) * sqrt( mixt_frac / ( 1 - mixt_frac ) )
    !                  * sqrt( <y'^2> ) * sgn( <w'y'> );
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
    ! Either w can be the setting variable and both x and y can be responding
    ! variables, x can be the setting variable and both w and y can be
    ! responding variables, y can be the setting variable and both w and x can
    ! be responding variables, or all of w, x, and y can be responding
    ! variables.
    !
    ! When w is the setting variable, coef_sigma_w_1_sqd and coef_sigma_w_2_sqd
    ! are given by:
    !
    ! coef_sigma_w_1_sqd = ( ( zeta_w + 1 ) * ( 1 - F_w ) )
    !                      / ( ( zeta_w + 2 ) * mixt_frac ); and
    !
    ! coef_sigma_w_2_sqd = ( 1 - F_w ) / ( ( zeta_w + 2 ) * ( 1 - mixt_frac ) ).
    !
    ! When w is a responding variable, coef_sigma_w_1_sqd and coef_sigma_w_2_sqd
    ! are given by:
    !
    ! coef_sigma_w_1_sqd = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skw
    !                      / ( 3 * mixt_frac * sqrt( F_w ) )
    !                      - ( 1 + mixt_frac ) * F_w / ( 3 * mixt_frac )
    !                      + 1; and
    !
    ! coef_sigma_w_2_sqd = ( ( 1 - F_w ) - mixt_frac * coef_sigma_w_1_sqd )
    !                      / ( 1 - mixt_frac ).
    !
    ! When x is the setting variable, coef_sigma_x_1_sqd and coef_sigma_x_2_sqd
    ! are given by:
    !
    ! coef_sigma_x_1_sqd = ( ( zeta_x + 1 ) * ( 1 - F_x ) )
    !                      / ( ( zeta_x + 2 ) * mixt_frac ); and
    !
    ! coef_sigma_x_2_sqd = ( 1 - F_x ) / ( ( zeta_x + 2 ) * ( 1 - mixt_frac ) ).
    !
    ! When x is a responding variable, coef_sigma_x_1_sqd and coef_sigma_x_2_sqd
    ! are given by:
    !
    ! coef_sigma_x_1_sqd = sqrt( mixt_frac * ( 1 - mixt_frac ) )
    !                      * Skx * sgn( <w'x'> )
    !                      / ( 3 * mixt_frac * sqrt( F_x ) )
    !                      - ( 1 + mixt_frac ) * F_x / ( 3 * mixt_frac )
    !                      + 1; and
    !
    ! coef_sigma_x_2_sqd = ( ( 1 - F_x ) - mixt_frac * coef_sigma_x_1_sqd )
    !                      / ( 1 - mixt_frac ).
    !
    ! When y is the setting variable, coef_sigma_y_1_sqd and coef_sigma_y_2_sqd
    ! are given by:
    !
    ! coef_sigma_y_1_sqd = ( ( zeta_y + 1 ) * ( 1 - F_y ) )
    !                      / ( ( zeta_y + 2 ) * mixt_frac ); and
    !
    ! coef_sigma_y_2_sqd = ( 1 - F_y ) / ( ( zeta_y + 2 ) * ( 1 - mixt_frac ) ).
    !
    ! When y is a responding variable, coef_sigma_y_1_sqd and coef_sigma_y_2_sqd
    ! are given by:
    !
    ! coef_sigma_y_1_sqd = sqrt( mixt_frac * ( 1 - mixt_frac ) )
    !                      * Sky * sgn( <w'y'> )
    !                      / ( 3 * mixt_frac * sqrt( F_y ) )
    !                      - ( 1 + mixt_frac ) * F_y / ( 3 * mixt_frac )
    !                      + 1; and
    !
    ! coef_sigma_y_2_sqd = ( ( 1 - F_y ) - mixt_frac * coef_sigma_y_1_sqd )
    !                      / ( 1 - mixt_frac ).
    !
    ! Additionally:
    !
    ! corr_w_x_1 = corr_w_x_2
    ! = ( <w'x'> - mixt_frac * ( mu_w_1 - <w> ) * ( mu_x_1 - <x> )
    !            - ( 1 - mixt_frac ) * ( mu_w_2 - <w> ) * ( mu_x_2 - <x> ) )
    !   / ( mixt_frac * sigma_w_1 * sigma_x_1
    !       + ( 1 - mixt_frac ) * sigma_w_2 * sigma_x_2 );
    !
    ! where -1 <= corr_w_x_1 = corr_w_x_2 <= 1.  This equation can be rewritten
    ! as:
    !
    ! corr_w_x_1 = corr_w_x_2
    ! = ( <w'x'>
    !     - sqrt( F_w ) * sqrt( F_x ) * sgn( <w'x'> )
    !       * sqrt( <w'^2> ) * sqrt( <x'^2 > ) )
    !   / ( ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !         + ( 1 - mixt_frac )
    !           * sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) )
    !       * sqrt( <w'^2> ) * sqrt( <x'^2> ) ).
    !
    ! Likewise:
    !
    ! corr_w_y_1 = corr_w_y_2
    ! = ( <w'y'>
    !     - sqrt( F_w ) * sqrt( F_y ) * sgn( <w'y'> )
    !       * sqrt( <w'^2> ) * sqrt( <y'^2 > ) )
    !   / ( ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_y_1_sqd )
    !         + ( 1 - mixt_frac )
    !           * sqrt( coef_sigma_w_2_sqd * coef_sigma_y_2_sqd ) )
    !       * sqrt( <w'^2> ) * sqrt( <y'^2> ) ); and
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
    ! = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( F_w ) * sqrt( <w'^2> )
    !   * ( sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
    !             - sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
    !     / ( mixt_frac * sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
    !         + ( 1 - mixt_frac )
    !           * sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
    !   * <x'y'>
    !   + sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( F_w ) * sqrt( <w'^2> )
    !     * sqrt( F_x ) * sqrt( <x'^2> ) * sgn( <w'x'> )
    !     * sqrt( F_y ) * sqrt( <y'^2> ) * sgn( <w'y'> )
    !     * ( ( 1 - mixt_frac ) / mixt_frac - mixt_frac / ( 1 - mixt_frac )
    !         - ( sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
    !             - sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
    !           / ( mixt_frac * sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
    !               + ( 1 - mixt_frac )
    !                 * sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
    !         - ( sqrt( coef_sigma_w_1_sqd * coef_sigma_y_1_sqd )
    !             - sqrt( coef_sigma_w_2_sqd * coef_sigma_y_2_sqd ) )
    !           / ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_y_1_sqd )
    !               + ( 1 - mixt_frac )
    !                 * sqrt( coef_sigma_w_2_sqd * coef_sigma_y_2_sqd ) )
    !         - ( sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !             - sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) )
    !           / ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !               + ( 1 - mixt_frac )
    !                 * sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) ) )
    !   + sqrt( mixt_frac * ( 1 - mixt_frac ) )
    !     * sqrt( F_x ) * sqrt( <x'^2> ) * sgn( <w'x'> )
    !     * ( sqrt( coef_sigma_w_1_sqd * coef_sigma_y_1_sqd )
    !         - sqrt( coef_sigma_w_2_sqd * coef_sigma_y_2_sqd ) )
    !       / ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_y_1_sqd )
    !           + ( 1 - mixt_frac )
    !             * sqrt( coef_sigma_w_2_sqd * coef_sigma_y_2_sqd ) )
    !     * <w'y'>
    !   + sqrt( mixt_frac * ( 1 - mixt_frac ) )
    !     * sqrt( F_y ) * sqrt( <y'^2> ) * sgn( <w'y'> )
    !     * ( sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !         - sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) )
    !       / ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !           + ( 1 - mixt_frac )
    !             * sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) )
    !     * <w'x'>.
    !
    ! This equation is of the form:
    !
    ! <w'^2 x'> = coef_wp2xp_implicit * <w'x'> + term_wp2xp_explicit;
    !
    ! where:
    !
    ! coef_wpxpyp_implicit
    ! = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( F_w ) * sqrt( <w'^2> )
    !   * ( sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
    !             - sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
    !     / ( mixt_frac * sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
    !         + ( 1 - mixt_frac )
    !           * sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) ); and
    !
    ! term_wpxpyp_explicit
    ! = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( F_w ) * sqrt( <w'^2> )
    !   * sqrt( F_x ) * sqrt( <x'^2> ) * sgn( <w'x'> )
    !   * sqrt( F_y ) * sqrt( <y'^2> ) * sgn( <w'y'> )
    !   * ( ( 1 - mixt_frac ) / mixt_frac - mixt_frac / ( 1 - mixt_frac )
    !       - ( sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
    !           - sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
    !         / ( mixt_frac * sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
    !             + ( 1 - mixt_frac )
    !               * sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
    !       - ( sqrt( coef_sigma_w_1_sqd * coef_sigma_y_1_sqd )
    !           - sqrt( coef_sigma_w_2_sqd * coef_sigma_y_2_sqd ) )
    !         / ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_y_1_sqd )
    !             + ( 1 - mixt_frac )
    !               * sqrt( coef_sigma_w_2_sqd * coef_sigma_y_2_sqd ) )
    !       - ( sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !           - sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) )
    !         / ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !             + ( 1 - mixt_frac )
    !               * sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) ) )
    !   + sqrt( mixt_frac * ( 1 - mixt_frac ) )
    !     * sqrt( F_x ) * sqrt( <x'^2> ) * sgn( <w'x'> )
    !     * ( sqrt( coef_sigma_w_1_sqd * coef_sigma_y_1_sqd )
    !         - sqrt( coef_sigma_w_2_sqd * coef_sigma_y_2_sqd ) )
    !       / ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_y_1_sqd )
    !           + ( 1 - mixt_frac )
    !             * sqrt( coef_sigma_w_2_sqd * coef_sigma_y_2_sqd ) )
    !     * <w'y'>
    !   + sqrt( mixt_frac * ( 1 - mixt_frac ) )
    !     * sqrt( F_y ) * sqrt( <y'^2> ) * sgn( <w'y'> )
    !     * ( sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !         - sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) )
    !       / ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
    !           + ( 1 - mixt_frac )
    !             * sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) )
    !     * <w'x'>.
    !
    ! There are also special cases for the above equations.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        one,  & ! Variable(s)
        zero

    use clubb_precision, only: &
        core_rknd    ! Procedure(s)

    implicit none

    ! Input Variables
    real ( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      wp2,                & ! Variance of w (overall)                  [m^2/s^2]
      xp2,                & ! Variance of x (overall)              [(x units)^2]
      yp2,                & ! Variance of y (overall)              [(y units)^2]
      wpxp,               & ! Covariance of w and x (overall)     [m/s(x units)]
      wpyp,               & ! Covariance of w and y (overall)     [m/s(y units)]
      sgn_wpxp,           & ! Sign of the covariance of w and x              [-]
      sgn_wpyp,           & ! Sign of the covariance of w and y              [-]
      mixt_frac,          & ! Mixture fraction                               [-]
      F_w,                & ! Parameter: spread of the PDF comp. means of w  [-]
      F_x,                & ! Parameter: spread of the PDF comp. means of x  [-]
      F_y,                & ! Parameter: spread of the PDF comp. means of y  [-]
      coef_sigma_w_1_sqd, & ! sigma_w_1^2 = coef_sigma_w_1_sqd * <w'^2>      [-]
      coef_sigma_w_2_sqd, & ! sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>      [-]
      coef_sigma_x_1_sqd, & ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>      [-]
      coef_sigma_x_2_sqd, & ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>      [-]
      coef_sigma_y_1_sqd, & ! sigma_y_1^2 = coef_sigma_y_1_sqd * <y'^2>      [-]
      coef_sigma_y_2_sqd    ! sigma_y_2^2 = coef_sigma_y_2_sqd * <y'^2>      [-]

    ! Output Variables
    ! Coefs.: <w'x'y'> = coef_wpxpyp_implicit * <x'y'> + term_wpxpyp_explicit
    real ( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      coef_wpxpyp_implicit, & ! Coefficient that is multiplied by <x'y'>   [m/s]
      term_wpxpyp_explicit    ! Term that is on the RHS  [m/s(x units)(y units)]

    ! Local Variables
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coefs_factor_wx, & ! Factor involving coef_sigma_... w and x coefs     [-]
      coefs_factor_wy, & ! Factor involving coef_sigma_... w and y coefs     [-]
      coefs_factor_xy    ! Factor involving coef_sigma_... x and y coefs     [-]


    ! Calculate coefs_factor_wx.
    where ( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd > zero &
            .or. coef_sigma_w_2_sqd * coef_sigma_x_2_sqd > zero )

       ! coefs_factor_wx
       ! = ( sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
       !     - sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) )
       !   / ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd )
       !       + ( 1 - mixt_frac )
       !         * sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) )
       coefs_factor_wx &
       = ( sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd ) &
           - sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) ) &
         / ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd ) &
             + ( one - mixt_frac ) &
               * sqrt( coef_sigma_w_2_sqd * coef_sigma_x_2_sqd ) )

    elsewhere ! coef_sigma_w_1_sqd * coef_sigma_x_1_sqd = 0
              ! and coef_sigma_w_2_sqd * coef_sigma_x_2_sqd = 0

       ! When coef_sigma_w_1_sqd * coef_sigma_x_1_sqd = 0 and
       ! coef_sigma_w_2_sqd * coef_sigma_x_2_sqd = 0, the value of
       ! coefs_factor_wx is undefined.  However, setting coefs_factor_wx to a
       ! value of 0 in this scenario allows for the use of general form
       ! equations below for coef_wpxpyp_implicit and term_wpxpyp_explicit.
       coefs_factor_wx = zero

    endwhere

    ! Calculate coefs_factor_wy.
    where ( coef_sigma_w_1_sqd * coef_sigma_y_1_sqd > zero &
            .or. coef_sigma_w_2_sqd * coef_sigma_y_2_sqd > zero )

       ! coefs_factor_wy
       ! = ( sqrt( coef_sigma_w_1_sqd * coef_sigma_y_1_sqd )
       !     - sqrt( coef_sigma_w_2_sqd * coef_sigma_y_2_sqd ) )
       !   / ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_y_1_sqd )
       !       + ( 1 - mixt_frac )
       !         * sqrt( coef_sigma_w_2_sqd * coef_sigma_y_2_sqd ) )
       coefs_factor_wy &
       = ( sqrt( coef_sigma_w_1_sqd * coef_sigma_y_1_sqd ) &
           - sqrt( coef_sigma_w_2_sqd * coef_sigma_y_2_sqd ) ) &
         / ( mixt_frac * sqrt( coef_sigma_w_1_sqd * coef_sigma_y_1_sqd ) &
             + ( one - mixt_frac ) &
               * sqrt( coef_sigma_w_2_sqd * coef_sigma_y_2_sqd ) )

    elsewhere ! coef_sigma_w_1_sqd * coef_sigma_y_1_sqd = 0
              ! and coef_sigma_w_2_sqd * coef_sigma_y_2_sqd = 0

       ! When coef_sigma_w_1_sqd * coef_sigma_y_1_sqd = 0 and
       ! coef_sigma_w_2_sqd * coef_sigma_y_2_sqd = 0, the value of
       ! coefs_factor_wy is undefined.  However, setting coefs_factor_wy to a
       ! value of 0 in this scenario allows for the use of general form
       ! equations below for coef_wpxpyp_implicit and term_wpxpyp_explicit.
       coefs_factor_wy = zero

    endwhere

    ! Calculate coefs_factor_xy.
    where ( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd > zero &
            .or. coef_sigma_x_2_sqd * coef_sigma_y_2_sqd > zero )

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

    elsewhere ! coef_sigma_x_1_sqd * coef_sigma_y_1_sqd = 0
              ! and coef_sigma_x_2_sqd * coef_sigma_y_2_sqd = 0

       ! When coef_sigma_x_1_sqd * coef_sigma_y_1_sqd = 0 and
       ! coef_sigma_x_2_sqd * coef_sigma_y_2_sqd = 0, the value of
       ! coefs_factor_xy is undefined.  However, setting coefs_factor_xy to a
       ! value of 0 in this scenario allows for the use of general form
       ! equations below for coef_wpxpyp_implicit and term_wpxpyp_explicit.
       coefs_factor_xy = zero

    endwhere


    ! Calculate coef_wpxpyp_implicit and term_wpxpyp_explicit.
    where ( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd > zero &
            .or. coef_sigma_x_2_sqd * coef_sigma_y_2_sqd > zero )

       coef_wpxpyp_implicit &
       = sqrt( mixt_frac * ( one - mixt_frac ) ) &
         * sqrt( F_w ) * sqrt( wp2 ) * coefs_factor_xy

       term_wpxpyp_explicit &
       = sqrt( mixt_frac * ( one - mixt_frac ) ) * sqrt( F_w ) * sqrt( wp2 ) &
         * sqrt( F_x ) * sqrt( xp2 ) * sgn_wpxp &
         * sqrt( F_y ) * sqrt( yp2 ) * sgn_wpyp &
         * ( ( one - mixt_frac ) / mixt_frac - mixt_frac / ( one - mixt_frac ) &
             - coefs_factor_xy - coefs_factor_wy - coefs_factor_wx ) &
         + sqrt( mixt_frac * ( one - mixt_frac ) ) &
           * sqrt( F_x ) * sqrt( xp2 ) * sgn_wpxp * coefs_factor_wy * wpyp &
         + sqrt( mixt_frac * ( one - mixt_frac ) ) &
           * sqrt( F_y ) * sqrt( yp2 ) * sgn_wpyp * coefs_factor_wx * wpxp

    elsewhere ! coef_sigma_x_1_sqd * coef_sigma_y_1_sqd = 0
              ! and coef_sigma_x_2_sqd * coef_sigma_y_2_sqd = 0

       coef_wpxpyp_implicit &
       = sqrt( mixt_frac * ( one - mixt_frac ) ) * sqrt( F_w ) * sqrt( wp2 ) &
         * ( ( one - mixt_frac ) / mixt_frac - mixt_frac / ( one - mixt_frac ) &
             - coefs_factor_wy - coefs_factor_wx )

       term_wpxpyp_explicit &
       = sqrt( mixt_frac * ( one - mixt_frac ) ) &
         * sqrt( F_x ) * sqrt( xp2 ) * sgn_wpxp * coefs_factor_wy * wpyp &
         + sqrt( mixt_frac * ( one - mixt_frac ) ) &
           * sqrt( F_y ) * sqrt( yp2 ) * sgn_wpyp * coefs_factor_wx * wpxp

    endwhere


    return

  end subroutine calc_coefs_wpxpyp_semiimpl

  !=============================================================================

end module new_pdf
