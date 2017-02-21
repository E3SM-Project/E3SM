!---------------------------------------------------------------------------
! $Id: calc_roots.F90 7827 2015-07-14 18:38:11Z bmg2@uwm.edu $
!===============================================================================
module calc_roots

  implicit none

  public :: cubic_solve, &
            quadratic_solve, &
            cube_root

  private ! Set Default Scope

  contains

  !=============================================================================
  pure function cubic_solve( a_coef, b_coef, c_coef, d_coef ) &
  result( roots )

    ! Description:
    ! Solve for the roots of x in a cubic equation.
    !
    ! The cubic equation has the form:
    !
    ! f(x) = a*x^3 + b*x^2 + c*x + d;
    !
    ! where a /= 0.  When f(x) = 0, the cubic formula is used to solve:
    !
    ! a*x^3 + b*x^2 + c*x + d = 0.
    !
    ! The cubic formula is also called Cardano's Formula.
    !
    ! The three solutions for x are:
    !
    ! x(1) = -(1/3)*(b/a) + ( S + T );
    ! x(2) = -(1/3)*(b/a) - (1/2) * ( S + T ) + (1/2)i * sqrt(3) * ( S - T );
    ! x(3) = -(1/3)*(b/a) - (1/2) * ( S + T ) - (1/2)i * sqrt(3) * ( S - T );
    !
    ! where:
    !
    ! S = ( R + sqrt( D ) )^(1/3); and
    ! T = ( R - sqrt( D ) )^(1/3).
    !
    ! The determinant, D, is given by:
    !
    ! D = R^2 + Q^3.
    !
    ! The values of R and Q relate back to the a, b, c, and d coefficients:
    !
    ! Q = ( 3*(c/a) - (b/a)^2 ) / 9; and
    ! R = ( 9*(b/a)*(c/a) - 27*(d/a) - 2*(b/a)^3 ) / 54.
    !
    ! When D < 0, there are three unique, real-valued roots.  When D = 0, there
    ! are three real-valued roots, but one root is a double root or a triple
    ! root.  When D > 0, there is one real-valued root and there are two roots
    ! that are complex conjugates.

    ! References:
    ! http://mathworld.wolfram.com/CubicFormula.html
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        three,     & ! Constant(s)
        two,       &
        one_half,  &
        one_third, &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      a_coef, & ! Coefficient a (of x^3) in a*x^3 + b*x^2 + c^x + d = 0    [-]
      b_coef, & ! Coefficient b (of x^2) in a*x^3 + b*x^2 + c^x + d = 0    [-]
      c_coef, & ! Coefficient c (of x) in a*x^3 + b*x^2 + c^x + d = 0      [-]
      d_coef    ! Coefficient d in a*x^3 + b*x^2 + c^x + d = 0             [-]

    ! Return Variables
    complex( kind = core_rknd ), dimension(3) :: &
      roots    ! Roots of x that satisfy a*x^3 + b*x^2 + c*x + d = 0       [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      cap_Q_coef,  & ! Coefficient Q in cubic formula     [-]
      cap_R_coef,  & ! Coefficient R in cubic formula     [-]
      determinant    ! Determinant D in cubic formula     [-]

    complex( kind = core_rknd ) :: &
      sqrt_det,   & ! Square root of determinant D in cubic formula     [-]
      cap_S_coef, & ! Coefficient S in cubic formula                    [-]
      cap_T_coef    ! Coefficient T in cubic formula                    [-]

    complex( kind = core_rknd ), parameter :: &
      i_cmplx = ( 0.0_core_rknd, 1.0_core_rknd )  ! i = sqrt(-1)

    complex( kind = core_rknd ) :: &
      sqrt_3,          & ! Sqrt 3 (complex data type)
      one_half_cmplx,  & ! 1/2 (complex data type)
      one_third_cmplx    ! 1/3 (complex data type)


    ! Declare some constants as complex data types in order to prevent
    ! data-type conversion warning messages.
    sqrt_3 = cmplx( sqrt( three ), kind = core_rknd )
    one_half_cmplx = cmplx( one_half, kind = core_rknd )
    one_third_cmplx = cmplx( one_third, kind = core_rknd )

    ! Find the value of the coefficient Q; where
    ! Q = ( 3*(c/a) - (b/a)^2 ) / 9.
    cap_Q_coef = ( three * (c_coef/a_coef) - (b_coef/a_coef)**2 ) &
                 / 9.0_core_rknd

    ! Find the value of the coefficient R; where
    ! R = ( 9*(b/a)*(c/a) - 27*(d/a) - 2*(b/a)^3 ) / 54.
    cap_R_coef = ( 9.0_core_rknd * (b_coef/a_coef) * (c_coef/a_coef) &
                   - 27.0_core_rknd * (d_coef/a_coef) &
                   - two * (b_coef/a_coef)**3 ) / 54.0_core_rknd

    ! Find the value of the determinant D; where
    ! D = R^2 + Q^3.
    determinant = cap_Q_coef**3 + cap_R_coef**2

    if ( determinant < zero ) then

       ! Calculate the square root of the determinant.  This will be a complex
       ! number.
       sqrt_det = sqrt( cmplx( determinant, kind = core_rknd ) )

       ! Find the value of the coefficient S; where
       ! S = ( R + sqrt( D ) )^(1/3).
       cap_S_coef &
       = ( cmplx( cap_R_coef, kind = core_rknd ) + sqrt_det )**one_third_cmplx

       ! Find the value of the coefficient T; where
       ! T = ( R - sqrt( D ) )^(1/3).
       cap_T_coef &
       = ( cmplx( cap_R_coef, kind = core_rknd ) - sqrt_det )**one_third_cmplx

    else ! determinant >= 0

       ! Find the value of the coefficient S; where
       ! S = ( R + sqrt( D ) )^(1/3).
       cap_S_coef &
       = cmplx( cube_root( cap_R_coef + sqrt( determinant ) ), &
                kind = core_rknd )

       ! Find the value of the coefficient T; where
       ! T = ( R - sqrt( D ) )^(1/3).
       cap_T_coef &
       = cmplx( cube_root( cap_R_coef - sqrt( determinant ) ), &
                kind = core_rknd )

    endif ! determinant < 0

    ! Find the values of the roots.
    ! This root is always real-valued.
    ! x(1) = -(1/3)*(b/a) + ( S + T ).
    roots(1) = - one_third_cmplx * cmplx( b_coef/a_coef, kind = core_rknd ) &
               + ( cap_S_coef + cap_T_coef )

    ! This root is real-valued when D < 0 (even though the square root of the
    ! determinant is a complex number), as well as when D = 0 (when it is part
    ! of a double or triple root).  When D > 0, this root is a complex number.
    ! It is the complex conjugate of roots(3).
    ! x(2) = -(1/3)*(b/a) - (1/2) * ( S + T ) + (1/2)i * sqrt(3) * ( S - T ).
    roots(2) = - one_third_cmplx * cmplx( b_coef/a_coef, kind = core_rknd ) &
               - one_half_cmplx * ( cap_S_coef + cap_T_coef ) &
               + one_half_cmplx * i_cmplx * sqrt_3 * ( cap_S_coef - cap_T_coef )

    ! This root is real-valued when D < 0 (even though the square root of the
    ! determinant is a complex number), as well as when D = 0 (when it is part
    ! of a double or triple root).  When D > 0, this root is a complex number.
    ! It is the complex conjugate of roots(2).
    ! x(3) = -(1/3)*(b/a) - (1/2) * ( S + T ) - (1/2)i * sqrt(3) * ( S - T ).
    roots(3) = - one_third_cmplx * cmplx( b_coef/a_coef, kind = core_rknd ) &
               - one_half_cmplx * ( cap_S_coef + cap_T_coef ) &
               - one_half_cmplx * i_cmplx * sqrt_3 * ( cap_S_coef - cap_T_coef )


    return

  end function cubic_solve

  !=============================================================================
  pure function quadratic_solve( a_coef, b_coef, c_coef ) &
  result( roots )

    ! Description:
    ! Solve for the roots of x in a quadratic equation.
    !
    ! The equation has the form:
    !
    ! f(x) = a*x^2 + b*x + c;
    !
    ! where a /= 0.  When f(x) = 0, the quadratic formula is used to solve:
    !
    ! a*x^2 + b*x + c = 0.
    !
    ! The two solutions for x are:
    !
    ! x(1) = ( -b + sqrt( b^2 - 4*a*c ) ) / (2*a); and
    ! x(2) = ( -b - sqrt( b^2 - 4*a*c ) ) / (2*a).
    !
    ! The determinant, D, is given by:
    !
    ! D = b^2 - 4*a*c.
    !
    ! When D > 0, there are two unique, real-valued roots.  When D = 0, there
    ! are two real-valued roots, but they are a double root.  When D < 0, there
    ! there are two roots that are complex conjugates.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        four, & ! Constant(s)
        two,  &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      a_coef, & ! Coefficient a (of x^2) in a*x^2 + b*x + c = 0    [-]
      b_coef, & ! Coefficient b (of x) in a*x^2 + b*x + c = 0      [-]
      c_coef    ! Coefficient c in a*x^2 + b*x + c = 0             [-]

    ! Return Variables
    complex( kind = core_rknd ), dimension(2) :: &
      roots    ! Roots of x that satisfy a*x^2 + b*x + c = 0       [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      determinant    ! Determinant D in quadratic formula     [-]

    complex( kind = core_rknd ) :: &
      sqrt_det    ! Square root of determinant D in quadratic formula     [-]


    ! Find the value of the determinant D; where
    ! D = b^2 - 4*a*c.
    determinant = b_coef**2 - four * a_coef * c_coef

    if ( determinant >= zero ) then

       ! Calculate the square root of the determinant.
       sqrt_det = cmplx( sqrt( determinant ), kind = core_rknd )

    else ! determinant < 0

       ! Calculate the square root of the determinant.  This will be a complex
       ! number.
       sqrt_det = sqrt( cmplx( determinant, kind = core_rknd ) )

    endif ! determinant >= 0

    ! Find the values of the roots.
    ! This root is real-valued when D > 0, as well as when D = 0 (when it is
    ! part of a double root).  When D < 0, this root is a complex number.  It is
    ! the complex conjugate of roots(2).
    ! x(1) = ( -b + sqrt( b^2 - 4*a*c ) ) / (2*a); and
    roots(1) = ( -cmplx( b_coef, kind = core_rknd ) + sqrt_det ) &
               / cmplx( two * a_coef, kind = core_rknd )

    ! This root is real-valued when D > 0, as well as when D = 0 (when it is
    ! part of a double root).  When D < 0, this root is a complex number.  It is
    ! the complex conjugate of roots(1).
    ! x(2) = ( -b - sqrt( b^2 - 4*a*c ) ) / (2*a).
    roots(2) = ( -cmplx( b_coef, kind = core_rknd ) - sqrt_det ) &
               / cmplx( two * a_coef, kind = core_rknd )


    return

  end function quadratic_solve

  !=============================================================================
  pure function cube_root( x )

    ! Description:
    ! Calculates the cube root of x.
    !
    ! When x >= 0, this code simply calculates x^(1/3).  When x < 0, this code
    ! uses x^(1/3) = -|x|^(1/3).  This eliminates numerical errors when the
    ! exponent of 1/3 is not treated as exactly 1/3, which would sometimes
    ! result in values of NaN.
    !
    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one_third, & ! Constant(s)
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      x  ! Variable x

    ! Return Variables
    real( kind = core_rknd ) :: &
      cube_root  ! Cube root of x


    if ( x >= zero ) then
       cube_root = x**one_third
    else ! x < 0
       cube_root = -abs(x)**one_third
    endif ! x >= 0


    return

  end function cube_root

!===============================================================================

end module calc_roots
