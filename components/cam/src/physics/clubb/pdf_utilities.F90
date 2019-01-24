!-------------------------------------------------------------------------
! $Id$
!===============================================================================
module pdf_utilities

  implicit none

  private ! Set default scope to private

  public :: mean_L2N,                  &
            mean_L2N_dp,               &
            stdev_L2N,                 &
            stdev_L2N_1lev,            &
            stdev_L2N_dp,              &
            corr_NL2NN,                &
            corr_NL2NN_dp,             &
            corr_NN2NL,                &
            corr_NN2NL_1lev,           &
            corr_LL2NN,                &
            corr_LL2NN_dp,             &
            corr_NN2LL,                &
            compute_mean_binormal,     &
            compute_variance_binormal, &
            calc_comp_corrs_binormal,  &
            calc_corr_chi_x,           &
            calc_corr_eta_x,           &
            calc_corr_rt_x,            &
            calc_corr_thl_x,           &
            calc_xp2

  contains

  !=============================================================================
  pure function mean_L2N( nz, mu_x, sigma2_on_mu2 )  &
  result( mu_x_n )
  
    ! Description:
    ! For a lognormally-distributed variable x, this function finds the mean of
    ! ln x (mu_x_n) for the ith component of the PDF, given the mean of x (mu_x)
    ! and the variance of x (sigma_sqd_x) for the ith component of the PDF. The
    ! value ln x is distributed normally when x is distributed lognormally.

    ! References:
    !  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
    !    Marcel Dekker, 401 pp.
    !  -- App. B.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one  ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz    ! Number of vertical levels                                    [-]

    real( kind = core_rknd ), dimension(nz), intent(in) ::  &
      mu_x,          & ! Mean of x (ith PDF component)                     [-]
      sigma2_on_mu2    ! Ratio:  sigma_x^2 / mu_x^2 (ith PDF component)    [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(nz) ::  &
      mu_x_n  ! Mean of ln x (ith PDF component)           [-]


    ! Find the mean of ln x for the ith component of the PDF.
    ! The max( mu_x / sqrt( 1 + sigma_x^2 / mu_x^2 ), tiny( mu_x ) ) statement
    ! is used to prevent taking ln 0, which will produce a result of -infinity.
    ! This would happen when mu_x is 0.  However, this code should not be
    ! entered when mu_x has a value of 0.
    mu_x_n = log( max( mu_x / sqrt( one + sigma2_on_mu2 ), tiny( mu_x ) ) )


    return

  end function mean_L2N

  !=============================================================================
  pure function mean_L2N_dp( mu_x, sigma2_on_mu2 )  &
  result( mu_x_n )
  
    ! Description:
    ! For a lognormally-distributed variable x, this function finds the mean of
    ! ln x (mu_x_n) for the ith component of the PDF, given the mean of x (mu_x)
    ! and the variance of x (sigma_sqd_x) for the ith component of the PDF. The
    ! value ln x is distributed normally when x is distributed lognormally.
    ! This function uses double precision variables.

    ! References:
    !  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
    !    Marcel Dekker, 401 pp.
    !  -- App. B.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one_dp  ! Constant(s)

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) ::  &
      mu_x,          & ! Mean of x (ith PDF component)                     [-]
      sigma2_on_mu2    ! Ratio:  sigma_x^2 / mu_x^2 (ith PDF component)    [-]

    ! Return Variable
    real( kind = dp ) ::  &
      mu_x_n  ! Mean of ln x (ith PDF component)           [-]


    ! Find the mean of ln x for the ith component of the PDF.
    mu_x_n = log( mu_x / sqrt( one_dp + sigma2_on_mu2 ) )


    return

  end function mean_L2N_dp

  !=============================================================================
  pure function stdev_L2N( nz, sigma2_on_mu2 )  &
  result( sigma_x_n )

    ! Description:
    ! For a lognormally-distributed variable x, this function finds the standard
    ! deviation of ln x (sigma_x_n) for the ith component of the PDF, given the
    ! mean of x (mu_x) and the variance of x (sigma_sqd_x) for the ith component
    ! of the PDF.  The value ln x is distributed normally when x is distributed
    ! lognormally.

    ! References:
    !  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
    !    Marcel Dekker, 401 pp.
    !  -- App. B.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one  ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz    ! Number of vertical levels                                    [-]

    real( kind = core_rknd ), dimension(nz), intent(in) ::  &
      sigma2_on_mu2    ! Ratio:  sigma_x^2 / mu_x^2 (ith PDF component)    [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(nz) ::  &
      sigma_x_n  ! Standard deviation of ln x (ith PDF component)   [-]


    ! Find the standard deviation of ln x for the ith component of the PDF.
    sigma_x_n = sqrt( log( one + sigma2_on_mu2 ) )


    return

  end function stdev_L2N

  !=============================================================================
  pure function stdev_L2N_1lev( sigma2_on_mu2 )  &
  result( sigma_x_n )

    ! Description:
    ! For a lognormally-distributed variable x, this function finds the standard
    ! deviation of ln x (sigma_x_n) for the ith component of the PDF, given the
    ! mean of x (mu_x) and the variance of x (sigma_sqd_x) for the ith component
    ! of the PDF.  The value ln x is distributed normally when x is distributed
    ! lognormally.

    ! References:
    !  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
    !    Marcel Dekker, 401 pp.
    !  -- App. B.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one  ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  &
      sigma2_on_mu2    ! Ratio:  sigma_x^2 / mu_x^2 (ith PDF component)    [-]

    ! Return Variable
    real( kind = core_rknd ) ::  &
      sigma_x_n  ! Standard deviation of ln x (ith PDF component)   [-]


    ! Find the standard deviation of ln x for the ith component of the PDF.
    sigma_x_n = sqrt( log( one + sigma2_on_mu2 ) )


    return

  end function stdev_L2N_1lev

  !=============================================================================
  pure function stdev_L2N_dp( sigma2_on_mu2 )  &
  result( sigma_x_n )

    ! Description:
    ! For a lognormally-distributed variable x, this function finds the standard
    ! deviation of ln x (sigma_x_n) for the ith component of the PDF, given the
    ! mean of x (mu_x) and the variance of x (sigma_sqd_x) for the ith component
    ! of the PDF.  The value ln x is distributed normally when x is distributed
    ! lognormally.
    ! This function uses double precision variables.

    ! References:
    !  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
    !    Marcel Dekker, 401 pp.
    !  -- App. B.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one_dp  ! Constant(s)

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) ::  &
      sigma2_on_mu2    ! Ratio:  sigma_x^2 / mu_x^2 (ith PDF component)    [-]

    ! Return Variable
    real( kind = dp ) ::  &
      sigma_x_n  ! Standard deviation of ln x (ith PDF component)   [-]


    ! Find the standard deviation of ln x for the ith component of the PDF.
    sigma_x_n = sqrt( log( one_dp + sigma2_on_mu2 ) )


    return

  end function stdev_L2N_dp

  !=============================================================================
  pure function corr_NL2NN( corr_x_y, sigma_y_n, y_sigma2_on_mu2 )  &
  result( corr_x_y_n )

    ! Description:
    ! For a normally-distributed variable x and a lognormally-distributed
    ! variable y, this function finds the correlation of x and ln y (corr_x_y_n)
    ! for the ith component of the PDF, given the correlation of x and y
    ! (corr_x_y) and the standard deviation of ln y (sigma_y_n) for the ith
    ! component of the PDF.  The value ln y is distributed normally when y is
    ! distributed lognormally.

    ! References:
    !  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
    !    Marcel Dekker, 401 pp.
    !  -- Eq. B-1.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        max_mag_correlation, & ! Constant(s)
        zero

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      corr_x_y,        & ! Correlation of x and y (ith PDF component)       [-]
      sigma_y_n,       & ! Standard deviation of ln y (ith PDF component)   [-]
      y_sigma2_on_mu2    ! Ratio:  sigma_y^2 / mu_y^2 (ith PDF component)   [-]

    ! Return Variable
    real( kind = core_rknd ) ::  &
      corr_x_y_n  ! Correlation of x and ln y (ith PDF component) [-]


    ! Find the correlation of x and ln y for the ith component of the PDF.
    ! When sigma_y = 0 and mu_y > 0, y_sigma2_on_mu2 = 0.  This results in
    ! sigma_y_n = 0.  The resulting corr_x_y_n is undefined.  However, the
    ! divide-by-zero problem needs to be addressed in the code.
    if ( sigma_y_n > zero ) then
       corr_x_y_n = corr_x_y * sqrt( y_sigma2_on_mu2 ) / sigma_y_n
    else ! sigma_y_n = 0
       ! The value of sqrt( y_sigma2_on_mu2 ) / sigma_y_n can be rewritten as:
       ! sqrt( y_sigma2_on_mu2 ) / sqrt( ln( 1 + y_sigma2_on_mu2 ) ).
       ! This can be further rewritten as:
       ! sqrt( y_sigma2_on_mu2 / ln( 1 + y_sigma2_on_mu2 ) ),
       ! which has a limit of 1 as y_sigma2_on_mu2 approaches 0 from the right.
       ! When sigma_y_n = 0, the value of corr_x_y_n is undefined, so set it
       ! to corr_x_y.
       corr_x_y_n = corr_x_y
    endif ! sigma_y_n > 0

    ! Clip the magnitude of the correlation of x and ln y in the ith PDF
    ! component, just in case the correlation (ith PDF component) of x and y and
    ! the standard deviation (ith PDF component) of ln y are inconsistent,
    ! resulting in an unrealizable value for corr_x_y_n.
    if ( corr_x_y_n > max_mag_correlation ) then
       corr_x_y_n = max_mag_correlation
    elseif ( corr_x_y_n < -max_mag_correlation ) then
       corr_x_y_n = -max_mag_correlation
    endif


    return

  end function corr_NL2NN

  !=============================================================================
  pure function corr_NL2NN_dp( corr_x_y, sigma_y_n, y_sigma2_on_mu2 )  &
  result( corr_x_y_n )

    ! Description:
    ! For a normally-distributed variable x and a lognormally-distributed
    ! variable y, this function finds the correlation of x and ln y (corr_x_y_n)
    ! for the ith component of the PDF, given the correlation of x and y
    ! (corr_x_y) and the standard deviation of ln y (sigma_y_n) for the ith
    ! component of the PDF.  The value ln y is distributed normally when y is
    ! distributed lognormally.
    ! This function uses double precision variables.

    ! References:
    !  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
    !    Marcel Dekker, 401 pp.
    !  -- Eq. B-1.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        max_mag_correlation, & ! Constant(s)
        zero_dp

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      corr_x_y,        & ! Correlation of x and y (ith PDF component)       [-]
      sigma_y_n,       & ! Standard deviation of ln y (ith PDF component)   [-]
      y_sigma2_on_mu2    ! Ratio:  sigma_y^2 / mu_y^2 (ith PDF component)   [-]

    ! Return Variable
    real( kind = dp ) ::  &
      corr_x_y_n  ! Correlation of x and ln y (ith PDF component) [-]


    ! Find the correlation of x and ln y for the ith component of the PDF.
    ! When sigma_y = 0 and mu_y > 0, y_sigma2_on_mu2 = 0.  This results in
    ! sigma_y_n = 0.  The resulting corr_x_y_n is undefined.  However, the
    ! divide-by-zero problem needs to be addressed in the code.
    if ( sigma_y_n > zero_dp ) then
       corr_x_y_n = corr_x_y * sqrt( y_sigma2_on_mu2 ) / sigma_y_n
    else ! sigma_y_n = 0
       ! The value of sqrt( y_sigma2_on_mu2 ) / sigma_y_n can be rewritten as:
       ! sqrt( y_sigma2_on_mu2 ) / sqrt( ln( 1 + y_sigma2_on_mu2 ) ).
       ! This can be further rewritten as:
       ! sqrt( y_sigma2_on_mu2 / ln( 1 + y_sigma2_on_mu2 ) ),
       ! which has a limit of 1 as y_sigma2_on_mu2 approaches 0 from the right.
       ! When sigma_y_n = 0, the value of corr_x_y_n is undefined, so set it
       ! to corr_x_y.
       corr_x_y_n = corr_x_y
    endif ! sigma_y_n > 0

    ! Clip the magnitude of the correlation of x and ln y in the ith PDF
    ! component, just in case the correlation (ith PDF component) of x and y and
    ! the standard deviation (ith PDF component) of ln y are inconsistent,
    ! resulting in an unrealizable value for corr_x_y_n.
    if ( corr_x_y_n > real( max_mag_correlation, kind = dp ) ) then
       corr_x_y_n = real( max_mag_correlation, kind = dp )
    elseif ( corr_x_y_n < -real( max_mag_correlation, kind = dp ) ) then
       corr_x_y_n = -real( max_mag_correlation, kind = dp )
    endif


    return

  end function corr_NL2NN_dp

  !=============================================================================
  pure function corr_NN2NL( nz, corr_x_y_n, sigma_y_n, y_sigma2_on_mu2 )  &
  result( corr_x_y )

    ! Description:
    ! For a normally-distributed variable x and a lognormally-distributed
    ! variable y, this function finds the correlation of x and y (corr_x_y) for
    ! the ith component of the PDF, given the correlation of x and ln y
    ! (corr_x_y_n) and the standard deviation of ln y (sigma_y_n) for the ith
    ! component of the PDF.  The value ln y is distributed normally when y is
    ! distributed lognormally.

    ! References:
    !  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
    !    Marcel Dekker, 401 pp.
    !  -- Eq. B-1.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        max_mag_correlation, & ! Constant(s)
        zero

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz    ! Number of vertical levels                                    [-]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      corr_x_y_n,      & ! Correlation of x and ln y (ith PDF component)    [-]
      sigma_y_n,       & ! Standard deviation of ln y (ith PDF component)   [-]
      y_sigma2_on_mu2    ! Ratio:  sigma_y^2 / mu_y^2 (ith PDF component)   [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(nz) ::  &
      corr_x_y  ! Correlation of x and y (ith PDF component) [-]


    ! Find the correlation of x and y for the ith component of the PDF.
    ! When sigma_y = 0 and mu_y > 0, y_sigma2_on_mu2 = 0.  This results in
    ! sigma_y_n = 0.  The resulting corr_x_y and corr_x_y_n are undefined.
    ! However, the divide-by-zero problem needs to be addressed in the code.
    where ( sigma_y_n > zero )
       ! Use the maximum of y_sigma2_on_mu2 and tiny( y_sigma2_on_mu2 ) instead
       ! of just y_sigma2_on_mu2.  The value of y_sigma2_on_mu2 must already be
       ! greater than 0 in order for this block of code to be entered (when
       ! y_sigma2_on_mu2 = 0, sigma_y_n = 0, and this block of code is not
       ! entered).  However, since a "where" statement is used here, this block
       ! of code may be erroneously entered when sigma_y_n and y_sigma2_on_mu2
       ! are equal to 0.  While these erroneous results are thrown away, they
       ! may result in a floating point error that can cause the run to stop.
       corr_x_y = corr_x_y_n * sigma_y_n &
                  / sqrt( max( y_sigma2_on_mu2, tiny( y_sigma2_on_mu2 ) ) )
    elsewhere ! sigma_y_n = 0
       ! The value of sigma_y_n / sqrt( y_sigma2_on_mu2 ) can be rewritten as:
       ! sqrt( ln( 1 + y_sigma2_on_mu2 ) ) / sqrt( y_sigma2_on_mu2 ).
       ! This can be further rewritten as:
       ! sqrt( ln( 1 + y_sigma2_on_mu2 ) / y_sigma2_on_mu2 ),
       ! which has a limit of 1 as y_sigma2_on_mu2 approaches 0 from the right.
       ! When sigma_y_n = 0, the value of corr_x_y is undefined, so set it
       ! to corr_x_y_n.
       corr_x_y = corr_x_y_n
    endwhere ! sigma_y_n > 0

    ! Clip the magnitude of the correlation of x and y in the ith PDF component,
    ! just in case the correlation (ith PDF component) of x and ln y and the
    ! standard deviation (ith PDF component) of ln y are inconsistent, resulting
    ! in an unrealizable value for corr_x_y.
    where ( corr_x_y > max_mag_correlation )
       corr_x_y = max_mag_correlation
    elsewhere ( corr_x_y < -max_mag_correlation )
       corr_x_y = -max_mag_correlation
    endwhere


    return

  end function corr_NN2NL

  !=============================================================================
  pure function corr_NN2NL_1lev( corr_x_y_n, sigma_y_n, y_sigma2_on_mu2 )  &
  result( corr_x_y )

    ! Description:
    ! For a normally-distributed variable x and a lognormally-distributed
    ! variable y, this function finds the correlation of x and y (corr_x_y) for
    ! the ith component of the PDF, given the correlation of x and ln y
    ! (corr_x_y_n) and the standard deviation of ln y (sigma_y_n) for the ith
    ! component of the PDF.  The value ln y is distributed normally when y is
    ! distributed lognormally.

    ! References:
    !  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
    !    Marcel Dekker, 401 pp.
    !  -- Eq. B-1.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        max_mag_correlation, & ! Constant(s)
        zero

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      corr_x_y_n,      & ! Correlation of x and ln y (ith PDF component)    [-]
      sigma_y_n,       & ! Standard deviation of ln y (ith PDF component)   [-]
      y_sigma2_on_mu2    ! Ratio:  sigma_y^2 / mu_y^2 (ith PDF component)   [-]

    ! Return Variable
    real( kind = core_rknd ) ::  &
      corr_x_y  ! Correlation of x and y (ith PDF component) [-]


    ! Find the correlation of x and y for the ith component of the PDF.
    ! When sigma_y = 0 and mu_y > 0, y_sigma2_on_mu2 = 0.  This results in
    ! sigma_y_n = 0.  The resulting corr_x_y and corr_x_y_n are undefined.
    ! However, the divide-by-zero problem needs to be addressed in the code.
    if ( sigma_y_n > zero ) then
       corr_x_y = corr_x_y_n * sigma_y_n / sqrt( y_sigma2_on_mu2 )
    else ! sigma_y_n = 0
       ! The value of sigma_y_n / sqrt( y_sigma2_on_mu2 ) can be rewritten as:
       ! sqrt( ln( 1 + y_sigma2_on_mu2 ) ) / sqrt( y_sigma2_on_mu2 ).
       ! This can be further rewritten as:
       ! sqrt( ln( 1 + y_sigma2_on_mu2 ) / y_sigma2_on_mu2 ),
       ! which has a limit of 1 as y_sigma2_on_mu2 approaches 0 from the right.
       ! When sigma_y_n = 0, the value of corr_x_y is undefined, so set it
       ! to corr_x_y_n.
       corr_x_y = corr_x_y_n
    endif ! sigma_y_n > 0

    ! Clip the magnitude of the correlation of x and y in the ith PDF component,
    ! just in case the correlation (ith PDF component) of x and ln y and the
    ! standard deviation (ith PDF component) of ln y are inconsistent, resulting
    ! in an unrealizable value for corr_x_y.
    if ( corr_x_y > max_mag_correlation ) then
       corr_x_y = max_mag_correlation
    elseif ( corr_x_y < -max_mag_correlation ) then
       corr_x_y = -max_mag_correlation
    endif


    return

  end function corr_NN2NL_1lev

  !=============================================================================
  pure function corr_LL2NN( corr_x_y, sigma_x_n, sigma_y_n, &
                            x_sigma2_on_mu2, y_sigma2_on_mu2 )  &
  result( corr_x_y_n )

    ! Description:
    ! For lognormally-distributed variables x and y, this function finds the
    ! correlation of ln x and ln y (corr_x_y_n) for the ith component of the
    ! PDF, given the correlation of x and y (corr_x_y), the standard deviation
    ! of ln x (sigma_x_n), and the standard deviation of ln y (sigma_y_n) for
    ! the ith component of the PDF.  The value of ln x (or ln y) is distributed
    ! normally when x (or y) is distributed lognormally.

    ! References:
    !  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
    !    Marcel Dekker, 401 pp.
    !  -- Eq. C-3.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,                 & ! Constant(s)
        zero,                &
        max_mag_correlation

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  &
      corr_x_y,        & ! Correlation of x and y (ith PDF component)       [-]
      sigma_x_n,       & ! Standard deviation of ln x (ith PDF component)   [-]
      sigma_y_n,       & ! Standard deviation of ln y (ith PDF component)   [-]
      x_sigma2_on_mu2, & ! Ratio:  sigma_x^2 / mu_x^2 (ith PDF component)   [-]
      y_sigma2_on_mu2    ! Ratio:  sigma_y^2 / mu_y^2 (ith PDF component)   [-]

    ! Return Variable
    real( kind = core_rknd ) ::  &
      corr_x_y_n  ! Correlation of ln x and ln y (ith PDF component)  [-]

    ! Local Variable
    real( kind = core_rknd ) ::  &
      log_arg    ! Input into the ln function    [-]


    ! Find the correlation of ln x and ln y for the ith component of the PDF.
    ! When sigma_x = 0 and mu_x > 0, x_sigma2_on_mu2 = 0.  This results in
    ! sigma_x_n = 0.  The resulting corr_x_y_n is undefined.  The same holds
    ! true when sigma_y = 0 and mu_y > 0.  However, the divide-by-zero problem
    ! needs to be addressed in the code.
    if ( sigma_x_n > zero .and. sigma_y_n > zero ) then
!      corr_x_y_n = log( one + corr_x_y * sqrt( exp( sigma_x_n**2 ) - one ) &
!                                       * sqrt( exp( sigma_y_n**2 ) - one )  ) &
!                   / ( sigma_x_n * sigma_y_n )
       log_arg = one + corr_x_y * sqrt( x_sigma2_on_mu2 * y_sigma2_on_mu2 )
       corr_x_y_n = log( log_arg ) / ( sigma_x_n * sigma_y_n )
    else ! sigma_x_n = 0 or sigma_y_n = 0
       ! The value of corr_x_y_n is undefined, so set it to corr_x_y.
       corr_x_y_n = corr_x_y
    endif ! sigma_x_n > 0 and sigma_y_n > 0

    ! Clip the magnitude of the correlation of ln x and ln y in the ith PDF
    ! component, just in case the correlation (ith PDF component) of x and y,
    ! the standard deviation (ith PDF component) of ln x, and the standard
    ! deviation (ith PDF component) of ln y are inconsistent, resulting in an
    ! unrealizable value for corr_x_y_n.
    if ( corr_x_y_n > max_mag_correlation ) then
       corr_x_y_n = max_mag_correlation
    elseif ( corr_x_y_n < -max_mag_correlation ) then
       corr_x_y_n = -max_mag_correlation
    endif


    return

  end function corr_LL2NN

  !=============================================================================
  pure function corr_LL2NN_dp( corr_x_y, sigma_x_n, sigma_y_n, &
                               x_sigma2_on_mu2, y_sigma2_on_mu2 )  &
  result( corr_x_y_n )

    ! Description:
    ! For lognormally-distributed variables x and y, this function finds the
    ! correlation of ln x and ln y (corr_x_y_n) for the ith component of the
    ! PDF, given the correlation of x and y (corr_x_y), the standard deviation
    ! of ln x (sigma_x_n), and the standard deviation of ln y (sigma_y_n) for
    ! the ith component of the PDF.  The value of ln x (or ln y) is distributed
    ! normally when x (or y) is distributed lognormally.
    ! This function uses double precision variables.

    ! References:
    !  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
    !    Marcel Dekker, 401 pp.
    !  -- Eq. C-3.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one_dp,              & ! Constant(s)
        zero_dp,             &
        max_mag_correlation

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! Input Variables
    real( kind = dp ), intent(in) ::  &
      corr_x_y,        & ! Correlation of x and y (ith PDF component)       [-]
      sigma_x_n,       & ! Standard deviation of ln x (ith PDF component)   [-]
      sigma_y_n,       & ! Standard deviation of ln y (ith PDF component)   [-]
      x_sigma2_on_mu2, & ! Ratio:  sigma_x^2 / mu_x^2 (ith PDF component)   [-]
      y_sigma2_on_mu2    ! Ratio:  sigma_y^2 / mu_y^2 (ith PDF component)   [-]


    ! Return Variable
    real( kind = dp ) ::  &
      corr_x_y_n  ! Correlation of ln x and ln y (ith PDF component)  [-]


    ! Find the correlation of ln x and ln y for the ith component of the PDF.
    ! When sigma_x = 0 and mu_x > 0, x_sigma2_on_mu2 = 0.  This results in
    ! sigma_x_n = 0.  The resulting corr_x_y_n is undefined.  The same holds
    ! true when sigma_y = 0 and mu_y > 0.  However, the divide-by-zero problem
    ! needs to be addressed in the code.
    if ( sigma_x_n > zero_dp .and. sigma_y_n > zero_dp ) then
       corr_x_y_n &
       = log( one_dp + corr_x_y * sqrt( x_sigma2_on_mu2 * y_sigma2_on_mu2 ) ) &
         / ( sigma_x_n * sigma_y_n )
    else ! sigma_x_n = 0 or sigma_y_n = 0
       ! The value of corr_x_y_n is undefined, so set it to corr_x_y.
       corr_x_y_n = corr_x_y
    endif ! sigma_x_n > 0 and sigma_y_n > 0

    ! Clip the magnitude of the correlation of ln x and ln y in the ith PDF
    ! component, just in case the correlation (ith PDF component) of x and y,
    ! the standard deviation (ith PDF component) of ln x, and the standard
    ! deviation (ith PDF component) of ln y are inconsistent, resulting in an
    ! unrealizable value for corr_x_y_n.
    if ( corr_x_y_n > real( max_mag_correlation, kind = dp ) ) then
       corr_x_y_n = real( max_mag_correlation, kind = dp )
    elseif ( corr_x_y_n < -real( max_mag_correlation, kind = dp ) ) then
       corr_x_y_n = -real( max_mag_correlation, kind = dp )
    endif


    return

  end function corr_LL2NN_dp

  !=============================================================================
  pure function corr_NN2LL( nz, corr_x_y_n, sigma_x_n, sigma_y_n, &
                            x_sigma2_on_mu2, y_sigma2_on_mu2 )  &
  result( corr_x_y )

    ! Description:
    ! For lognormally-distributed variables x and y, this function finds the
    ! correlation of x and y (corr_x_y) for the ith component of the PDF, given
    ! the correlation of ln x and ln y (corr_x_y_n), the standard deviation of
    ! ln x (sigma_x_n), and the standard deviation of ln y (sigma_y_n) for
    ! the ith component of the PDF.  The value of ln x (or ln y) is distributed
    ! normally when x (or y) is distributed lognormally.

    ! References:
    !  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
    !    Marcel Dekker, 401 pp.
    !  -- Eq. C-3.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,                 & ! Constant(s)
        zero,                &
        max_mag_correlation

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz    ! Number of vertical levels                                    [-]

    real( kind = core_rknd ), dimension(nz), intent(in) ::  &
      corr_x_y_n,      & ! Correlation of ln x and ln y (ith PDF component) [-]
      sigma_x_n,       & ! Standard deviation of ln x (ith PDF component)   [-]
      sigma_y_n,       & ! Standard deviation of ln y (ith PDF component)   [-]
      x_sigma2_on_mu2, & ! Ratio:  sigma_x^2 / mu_x^2 (ith PDF component)   [-]
      y_sigma2_on_mu2    ! Ratio:  sigma_y^2 / mu_y^2 (ith PDF component)   [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(nz) ::  &
      corr_x_y  ! Correlation of x and y (ith PDF component)  [-]


    ! Find the correlation of x and y for the ith component of the PDF.
    ! When sigma_x = 0 and mu_x > 0, x_sigma2_on_mu2 = 0.  This results in
    ! sigma_x_n = 0.  The resulting corr_x_y and corr_x_y_n are undefined.  The
    ! same holds true when sigma_y = 0 and mu_y > 0.  However, the
    ! divide-by-zero problem needs to be addressed in the code.
    where ( sigma_x_n > zero .and. sigma_y_n > zero )
!       corr_x_y = ( exp( sigma_x_n * sigma_y_n * corr_x_y_n ) - one ) &
!                  / ( sqrt( exp( sigma_x_n**2 ) - one ) &
!                      * sqrt( exp( sigma_y_n**2 ) - one ) )
       corr_x_y = ( exp( sigma_x_n * sigma_y_n * corr_x_y_n ) - one ) &
                  / sqrt( x_sigma2_on_mu2 * y_sigma2_on_mu2 )
    elsewhere ! sigma_x_n = 0 or sigma_y_n = 0
       ! The value of corr_x_y is undefined, so set it to corr_x_y_n.
       corr_x_y = corr_x_y_n
    endwhere ! sigma_x_n > 0 and sigma_y_n > 0

    ! Clip the magnitude of the correlation of x and y in the ith PDF component,
    ! just in case the correlation (ith PDF component) of ln x and ln y, the
    ! standard deviation (ith PDF component) of ln x, and the standard deviation
    ! (ith PDF component) of ln y are inconsistent, resulting in an unrealizable
    ! value for corr_x_y.
    where ( corr_x_y > max_mag_correlation )
       corr_x_y = max_mag_correlation
    elsewhere ( corr_x_y < -max_mag_correlation )
       corr_x_y = -max_mag_correlation
    endwhere


    return

  end function corr_NN2LL

  !=============================================================================
  elemental function compute_mean_binormal( mu_x_1, mu_x_2, mixt_frac ) &
  result( xm )

    ! Description:
    ! Computes the overall grid-box mean of a binormal distribution from the
    ! mean of each component

    ! References:
    !   None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Constant

    use constants_clubb, only: &
        one ! Constant

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_x_1,    & ! First PDF component mean of 'x'                       [?]
      mu_x_2,    & ! Second PDF component mean of 'x'                      [?]
      mixt_frac    ! Weight of the first PDF component                     [-]

    ! Output Variables
    real( kind = core_rknd ) :: &
      xm           ! Mean of 'x' (overall)                                 [?]

    !-----------------------------------------------------------------------

    !----- Begin Code -----
    xm = mixt_frac * mu_x_1 + ( one - mixt_frac ) * mu_x_2


    return

  end function compute_mean_binormal

  !=============================================================================
  elemental function compute_variance_binormal( xm, mu_x_1, mu_x_2, &
                                                stdev_x_1, stdev_x_2, &
                                                mixt_frac ) &
  result( xp2 )

    ! Description:
    ! Computes the overall grid-box variance of a binormal distribution from the
    ! variance of each component.

    ! References:
    !   None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Constant

    use constants_clubb, only: &
        one ! Constant

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      xm,        & ! Overall mean of 'x'                                   [?]
      mu_x_1,    & ! First PDF component mean of 'x'                       [?]
      mu_x_2,    & ! Second PDF component mean of 'x'                      [?]
      stdev_x_1, & ! Standard deviation of 'x' in the first PDF component  [?]
      stdev_x_2, & ! Standard deviation of 'x' in the second PDF component [?]
      mixt_frac    ! Weight of the first PDF component                     [-]

    ! Output Variables
    real( kind = core_rknd ) :: &
      xp2          ! Variance of 'x' (overall)                             [?^2]

    !-----------------------------------------------------------------------

    !----- Begin Code -----
    xp2 = mixt_frac * ( ( mu_x_1 - xm )**2 + stdev_x_1**2 ) &
          + ( one - mixt_frac ) * ( ( mu_x_2 - xm )**2 + stdev_x_2**2 )


    return

  end function compute_variance_binormal

  !=============================================================================
  subroutine calc_comp_corrs_binormal( xpyp, xm, ym, mu_x_1, mu_x_2,  & ! In
                                       mu_y_1, mu_y_2, sigma_x_1_sqd, & ! In
                                       sigma_x_2_sqd, sigma_y_1_sqd,  & ! In
                                       sigma_y_2_sqd, mixt_frac,      & ! In
                                       corr_x_y_1, corr_x_y_2         ) ! Out

    ! Description:
    ! Calculates the PDF component correlations of variables x and y, where
    ! x and y are both distributed as two-component normals (or binormals).
    ! The PDF component correlations are set equal to each other.
    !
    ! The overall covariance of x and y, <x'y'>, can be expressed in terms of
    ! PDF parameters by integrating over the PDF:
    !
    ! <x'y'> = INT(-inf:inf) INT(-inf:inf) ( x - <x> ) ( y - <y> ) P(x,y) dy dx;
    !
    ! where <x> is the overall mean of x, <y> is the overall mean of y, and
    ! P(x,y) is the equation for the two-component normal PDF of x and y.
    !
    ! The integral is evaluated, and the equation for <x'y'> is:
    !
    ! <x'y'> = mixt_frac * ( ( mu_x_1 - <x> ) * ( mu_y_1 - <y> )
    !                        + corr_x_y_1 * sigma_x_1 * sigma_y_1 )
    !          + ( 1 - mixt_frac ) * ( ( mu_x_2 - <x> ) * ( mu_y_2 - <y> )
    !                                  + corr_x_y_2 * sigma_x_2 * sigma_y_2 );
    !
    ! where mu_x_1 is the mean of x in the 1st PDF component, mu_x_2 is the mean
    ! of x in the 2nd PDF component, mu_y_1 is the mean of y in the 1st PDF
    ! component, mu_y_2 is the mean of y in the 2nd PDF component, sigma_x_1 is
    ! the standard deviation of x in the 1st PDF component, sigma_x_2 is the
    ! standard deviation of x in the 2nd PDF component, sigma_y_1 is the
    ! standard deviation of y in the 1st PDF component, sigma_y_2 is the
    ! standard deviation of y in the 2nd PDF component, corr_x_y_1 is the
    ! correlation of x and y in the 1st PDF component, corr_x_y_2 is the
    ! correlation of x and y in the 2nd PDF component, and mixt_frac is the
    ! mixture fraction (weight of the 1st PDF component).
    !
    ! This equation can be rewritten as:
    !
    ! <x'y'> = mixt_frac * ( mu_x_1 - <x> ) * ( mu_y_1 - <y> )
    !          + mixt_frac * corr_x_y_1 * sigma_x_1 * sigma_y_1
    !          + ( 1 - mixt_frac ) * ( mu_x_2 - <x> ) * ( mu_y_2 - <y> )
    !          + ( 1 - mixt_frac ) * corr_x_y_2 * sigma_x_2 * sigma_y_2.
    !
    ! Setting the two PDF component correlations equal to each other
    ! (corr_x_y_1 = corr_x_y_2), the equation can be solved for the PDF
    ! component correlations:
    !
    ! corr_x_y_1 = corr_x_y_2
    ! = ( <x'y'> - mixt_frac * ( mu_x_1 - <x> ) * ( mu_y_1 - <y> )
    !            - ( 1 - mixt_frac ) * ( mu_x_2 - <x> ) * ( mu_y_2 - <y> ) )
    !   / ( mixt_frac * sigma_x_1 * sigma_y_1
    !       + ( 1 - mixt_frac ) * sigma_x_2 * sigma_y_2 );
    !
    ! where -1 <= corr_x_y_1 = corr_x_y_2 <= 1.
    !
    ! When sigma_x_1 * sigma_y_1 = 0 and sigma_x_2 * sigma_y_2 = 0, at least one
    ! of x or y are constant within each PDF component, and both PDF component
    ! correlations are undefined.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        max_mag_correlation, & ! Variable(s)
        one, &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real ( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      xpyp,          & ! Covariance of x and y (overall)    [(x units)(y units)]
      xm,            & ! Mean of x (overall)                [x units]
      ym,            & ! Mean of y (overall)                [y units]
      mu_x_1,        & ! Mean of x (1st PDF component)      [x units]
      mu_x_2,        & ! Mean of x (2nd PDF component)      [x units]
      mu_y_1,        & ! Mean of y (1st PDF component)      [y units]
      mu_y_2,        & ! Mean of y (2nd PDF component)      [y units]
      sigma_x_1_sqd, & ! Variance of x (1st PDF component)  [(x units)^2]
      sigma_x_2_sqd, & ! Variance of x (2nd PDF component)  [(x units)^2]
      sigma_y_1_sqd, & ! Variance of y (1st PDF component)  [(y units)^2]
      sigma_y_2_sqd, & ! Variance of y (2nd PDF component)  [(y units)^2]
      mixt_frac        ! Mixture fraction                   [-]

    ! Output Variables
    real ( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      corr_x_y_1, & ! Correlation of x and y (1st PDF component)    [-]
      corr_x_y_2    ! Correlation of x and y (2nd PDF component)    [-]


    where ( sigma_x_1_sqd * sigma_y_1_sqd > zero &
            .or. sigma_x_2_sqd * sigma_y_2_sqd > zero )

       ! Calculate corr_x_y_1 (which also equals corr_x_y_2).
       corr_x_y_1 &
       = ( xpyp &
           - mixt_frac * ( mu_x_1 - xm ) * ( mu_y_1 - ym ) &
           - ( one - mixt_frac ) * ( mu_x_2 - xm ) * ( mu_y_2 - ym ) ) &
         / ( mixt_frac * sqrt( sigma_x_1_sqd * sigma_y_1_sqd ) &
             + ( one - mixt_frac ) * sqrt( sigma_x_2_sqd * sigma_y_2_sqd ) )

       ! The correlation must fall within the bounds of
       ! -max_mag_correlation < corr_x_y_1 (= corr_x_y_2) < max_mag_correlation
       corr_x_y_1 = max( -max_mag_correlation, min( max_mag_correlation, corr_x_y_1 ) )

    elsewhere ! sigma_x_1^2 * sigma_y_1^2 = 0 and sigma_x_2^2 * sigma_y_2^2 = 0.

       ! The correlation is undefined (output as 0).
       corr_x_y_1 = zero

    endwhere

    ! Set corr_x_y_2 equal to corr_x_y_1.
    corr_x_y_2 = corr_x_y_1


    return

  end subroutine calc_comp_corrs_binormal

  !=============================================================================
  pure function calc_corr_chi_x( crt_i, cthl_i, sigma_rt_i, sigma_thl_i,  &
                                 sigma_chi_i, corr_rt_x_i, corr_thl_x_i )  &
  result( corr_chi_x_i )

    ! Description:
    ! This function calculates the correlation of extended liquid water mixing
    ! ratio, chi (old s), and a generic variable x, within the ith component of
    ! the PDF.  The variable chi can be split into mean and turbulent
    ! components, such that:
    !
    ! chi = <chi> + chi';
    !
    ! where < > denotes a mean field an ' denotes a turbulent component.
    !
    ! The linearized equation for chi' is given in Larson et al. (2001), where
    ! within the ith component of the PDF:
    !
    ! chi_(i)' = Coef_rt(i) * r_t(i)' - Coef_thl(i) * th_l(i)'.
    !
    ! The equation for chi' can be multiplied by x'.  The equation becomes:
    !
    ! chi'x'_(i) = Coef_rt(i) * r_t'x'_(i) - Coef_thl(i) * th_l'x'_(i).
    !
    ! Averaging both sides, the covariance <chi'x'> is given by the equation:
    !
    ! <chi'x'_(i)> = Coef_rt(i) * <r_t'x'_(i)> - Coef_thl(i) * <th_l'x'_(i)>.
    !
    ! This equation can be rewritten as:
    !
    ! sigma_chi(i) * sigma_x(i) * corr_chi_x(i)
    !   = Coef_rt(i) * sigma_rt(i) * sigma_x(i) * corr_rt_x(i)
    !     - Coef_thl(i) * sigma_thl(i) * sigma_x(i) * corr_thl_x(i).
    !
    ! This equation can be solved for corr_chi_x(i):
    !
    ! corr_chi_x(i)
    ! = Coef_rt(i) * ( sigma_rt(i) / sigma_chi(i) ) * corr_rt_x(i)
    !   - Coef_thl(i) * ( sigma_thl(i) / sigma_chi(i) ) * corr_thl_x(i).
    !
    ! The correlation of chi and x within the ith component of the PDF is
    ! calculated.

    ! References:
    ! Eq. (13) and Eq. (14) of Larson, V. E., R. Wood, P. R. Field, J.-C. Golaz,
    ! T. H. Vonder Haar, W. R. Cotton, 2001:  Systematic Biases in the
    ! Microphysics and Thermodynamics of Numerical Models That Ignore
    ! Subgrid-Scale Variability. J. Atmos. Sci., 58, 1117--1128,
    ! doi:https://doi.org/10.1175/1520-0469(2001)058%3C1117:SBITMA%3E2.0.CO;2.
    !
    ! Eq. (A29) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        zero,                & ! Constant(s)
        chi_tol,             &
        max_mag_correlation

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      crt_i,        & ! Coefficient of r_t for chi (old s) (ith PDF comp.)   [-]
      cthl_i,       & ! Coefficient of th_l for chi (ith PDF comp.)  [(kg/kg)/K]
      sigma_rt_i,   & ! Standard deviation of r_t (ith PDF component)    [kg/kg]
      sigma_thl_i,  & ! Standard deviation of th_l (ith PDF component)       [K]
      sigma_chi_i,  & ! Standard deviation of chi (ith PDF component)    [kg/kg]
      corr_rt_x_i,  & ! Correlation of r_t and x (ith PDF component)         [-]
      corr_thl_x_i    ! Correlation of th_l and x (ith PDF component)        [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz) :: &
      corr_chi_x_i  ! Correlation of chi and x (ith PDF component)   [-]


    ! Calculate the correlation of chi and x in the ith PDF component.
    where ( sigma_chi_i > chi_tol )

       corr_chi_x_i = crt_i * ( sigma_rt_i / sigma_chi_i ) * corr_rt_x_i  &
                      - cthl_i * ( sigma_thl_i / sigma_chi_i ) * corr_thl_x_i

    elsewhere  ! sigma_chi_i = 0

       ! The standard deviation of chi in the ith PDF component is 0.  This
       ! means that chi is constant within the ith PDF component, and the ith
       ! PDF component covariance of chi and x is also 0.  The correlation of
       ! chi and x is undefined in the ith PDF component, so a value of 0 will
       ! be used.
       corr_chi_x_i = zero

    endwhere

    ! Clip the magnitude of the correlation of chi and x in the ith PDF
    ! component, just in case the correlations and standard deviations used in
    ! calculating it are inconsistent, resulting in an unrealizable value for
    ! corr_chi_x_i.
    where ( corr_chi_x_i > max_mag_correlation )
       corr_chi_x_i = max_mag_correlation
    elsewhere ( corr_chi_x_i < -max_mag_correlation )
       corr_chi_x_i = -max_mag_correlation
    endwhere


    return

  end function calc_corr_chi_x

  !=============================================================================
  pure function calc_corr_eta_x( crt_i, cthl_i, sigma_rt_i, sigma_thl_i,  &
                                 sigma_eta_i, corr_rt_x_i, corr_thl_x_i )  &
  result( corr_eta_x_i )

    ! Description:
    ! This function calculates the correlation of the variable that is
    ! orthogonal to extended liquid water mixing ratio in a PDF transformation,
    ! eta (old t), and a generic variable x, within the ith component of
    ! the PDF.

    ! References:
    ! Eq. (A30) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    ! Representation of Hydrometeors and Microphysical Feedback Effects Using a
    ! Multivariate PDF.  Doctoral dissertation, University of
    ! Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    ! http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        zero,                & ! Constant(s)
        eta_tol,             &
        max_mag_correlation

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      crt_i,        & ! Coefficient of r_t for chi (old s) (ith PDF comp.)   [-]
      cthl_i,       & ! Coefficient of th_l for chi (ith PDF comp.)  [(kg/kg)/K]
      sigma_rt_i,   & ! Standard deviation of r_t (ith PDF component)    [kg/kg]
      sigma_thl_i,  & ! Standard deviation of th_l (ith PDF component)       [K]
      sigma_eta_i,  & ! Standard deviation of eta (ith PDF component)    [kg/kg]
      corr_rt_x_i,  & ! Correlation of r_t and x (ith PDF component)         [-]
      corr_thl_x_i    ! Correlation of th_l and x (ith PDF component)        [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz) :: &
      corr_eta_x_i  ! Correlation of eta and x (ith PDF component)   [-]


    ! Calculate the correlation of eta and x in the ith PDF component.
    where ( sigma_eta_i > eta_tol )

       corr_eta_x_i = crt_i * ( sigma_rt_i / sigma_eta_i ) * corr_rt_x_i  &
                      + cthl_i * ( sigma_thl_i / sigma_eta_i ) * corr_thl_x_i

    elsewhere  ! sigma_eta_i = 0

       ! The standard deviation of eta in the ith PDF component is 0.  This
       ! means that eta is constant within the ith PDF component, and the ith
       ! PDF component covariance of eta and x is also 0.  The correlation of
       ! eta and x is undefined in the ith PDF component, so a value of 0 will
       ! be used.
       corr_eta_x_i = zero

    endwhere

    ! Clip the magnitude of the correlation of eta and x in the ith PDF
    ! component, just in case the correlations and standard deviations used in
    ! calculating it are inconsistent, resulting in an unrealizable value for
    ! corr_eta_x_i.
    where ( corr_eta_x_i > max_mag_correlation )
       corr_eta_x_i = max_mag_correlation
    elsewhere ( corr_eta_x_i < -max_mag_correlation )
       corr_eta_x_i = -max_mag_correlation
    endwhere


    return

  end function calc_corr_eta_x

  !=============================================================================
  pure function calc_corr_rt_x( crt_i, sigma_rt_i, sigma_chi_i, &
                                sigma_eta_i, corr_chi_x_i, corr_eta_x_i )  &
  result( corr_rt_x_i )

    ! Description:
    ! This function calculates the correlation of rt and x based on the
    ! correlation of chi and x and the correlation of eta and x.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two,                 & ! Constant(s)
        zero,                &
        rt_tol,              &
        max_mag_correlation

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      crt_i,        & ! Coef. of r_t in chi/eta eqns. (ith PDF component)    [-]
      sigma_rt_i,   & ! Standard deviation of r_t (ith PDF component)    [kg/kg]
      sigma_chi_i,  & ! Standard deviation of chi (ith PDF component)    [kg/kg]
      sigma_eta_i,  & ! Standard deviation of eta (ith PDF component)    [kg/kg]
      corr_chi_x_i, & ! Correlation of chi and x (ith PDF component)         [-]
      corr_eta_x_i    ! Correlation of eta and x (ith PDF component)         [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      corr_rt_x_i   ! Correlation of rt and x (ith PDF component)            [-]


    ! Calculate the correlation of rt and x in the ith PDF component.
    if ( sigma_rt_i > rt_tol ) then

       corr_rt_x_i = ( sigma_eta_i * corr_eta_x_i &
                       + sigma_chi_i * corr_chi_x_i ) &
                     / ( two * crt_i * sigma_rt_i )

    else  ! sigma_rt_i = 0

       ! The standard deviation of rt in the ith PDF component is 0.  This means
       ! that rt is constant within the ith PDF component, and the ith PDF
       ! component covariance of rt and x is also 0.  The correlation of rt and
       ! x is undefined in the ith PDF component, so a value of 0 will be used.
       corr_rt_x_i = zero

    endif

    ! Clip the magnitude of the correlation of rt and x in the ith PDF
    ! component, just in case the correlations and standard deviations used in
    ! calculating it are inconsistent, resulting in an unrealizable value for
    ! corr_rt_x_i.
    if ( corr_rt_x_i > max_mag_correlation ) then
       corr_rt_x_i = max_mag_correlation
    elseif ( corr_rt_x_i < -max_mag_correlation ) then
       corr_rt_x_i = -max_mag_correlation
    endif


    return

  end function calc_corr_rt_x

  !=============================================================================
  pure function calc_corr_thl_x( cthl_i, sigma_thl_i, sigma_chi_i, &
                                 sigma_eta_i, corr_chi_x_i, corr_eta_x_i )  &
  result( corr_thl_x_i )

    ! Description:
    ! This function calculates the correlation of thl and x based on the
    ! correlation of chi and x and the correlation of eta and x.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two,                 & ! Constant(s)
        zero,                &
        thl_tol,             &
        max_mag_correlation

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      cthl_i,       & ! Coef. of thl:  chi/eta eqns. (ith PDF comp.) [(kg/kg)/K]
      sigma_thl_i,  & ! Standard deviation of thl (ith PDF component)        [K]
      sigma_chi_i,  & ! Standard deviation of chi (ith PDF component)    [kg/kg]
      sigma_eta_i,  & ! Standard deviation of eta (ith PDF component)    [kg/kg]
      corr_chi_x_i, & ! Correlation of chi and x (ith PDF component)         [-]
      corr_eta_x_i    ! Correlation of eta and x (ith PDF component)         [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      corr_thl_x_i    ! Correlation of thl and x (ith PDF component)         [-]


    ! Calculate the correlation of thl and x in the ith PDF component.
    if ( sigma_thl_i > thl_tol ) then

       corr_thl_x_i = ( sigma_eta_i * corr_eta_x_i &
                        - sigma_chi_i * corr_chi_x_i ) &
                      / ( two * cthl_i * sigma_thl_i )

    else  ! sigma_thl_i = 0

       ! The standard deviation of thl in the ith PDF component is 0.  This
       ! means that thl is constant within the ith PDF component, and the ith
       ! PDF component covariance of thl and x is also 0.  The correlation of
       ! thl and x is undefined in the ith PDF component, so a value of 0 will
       ! be used.
       corr_thl_x_i = zero

    endif

    ! Clip the magnitude of the correlation of thl and x in the ith PDF
    ! component, just in case the correlations and standard deviations used in
    ! calculating it are inconsistent, resulting in an unrealizable value for
    ! corr_thl_x_i.
    if ( corr_thl_x_i > max_mag_correlation ) then
       corr_thl_x_i = max_mag_correlation
    elseif ( corr_thl_x_i < -max_mag_correlation ) then
       corr_thl_x_i = -max_mag_correlation
    endif


    return

  end function calc_corr_thl_x

  !=============================================================================
  pure function calc_xp2( mu_x_1, mu_x_2, &
                          mu_x_1_n, mu_x_2_n, &
                          sigma_x_1, sigma_x_2, &
                          sigma_x_1_n, sigma_x_2_n, &
                          mixt_frac, x_frac_1, x_frac_2, &
                          x_mean )  &
  result( xp2 )

    ! Description:
    ! Calculates the overall variance of x, <x'^2>, where the distribution of x
    ! is a combination of a lognormal distribution and/or 0 in each PDF
    ! component.  The fraction of each component where x is lognormally
    ! distributed (amd greater than 0) is x_frac_i (x_frac_1 and x_frac_2 for
    ! PDF components 1 and 2, respectively).  The fraction of each component
    ! where x has a value of 0 is ( 1 - x_frac_i ).  This function should be
    ! called to calculate the total variance for x when <x'^2> is not provided
    ! by a predictive (or other) equation.
    !    
    ! This function is used to calculate the overall variance for rain water
    ! mixing ratio, <r_r'^2>, and the overall variance for rain drop
    ! concentration, <N_r'^2>.  The ratio of variance to mean-value-squared is
    ! specified for the in-precip values of r_r and N_r within each PDF
    ! component, allowing for the calculation of sigma_rr_i and sigma_Nr_i,
    ! as well as sigma_rr_i_n and sigma_Nr_i_n.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two,  & ! Constant(s)
        one,  &
        zero, &
        eps

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_x_1,      & ! Mean of x (1st PDF comp.) in x_frac                  [-]
      mu_x_2,      & ! Mean of x (2nd PDF comp.) in x_frac                  [-]
      mu_x_1_n,    & ! Mean of ln x (1st PDF comp.) in x_frac               [-]
      mu_x_2_n,    & ! Mean of ln x (2nd PDF comp.) in x_frac               [-]
      sigma_x_1,   & ! Standard deviation of x (1st PDF comp.) in x_frac    [-]
      sigma_x_2,   & ! Standard deviation of x (2nd PDF comp.) in x_frac    [-]
      sigma_x_1_n, & ! Standard deviation of ln x (1st PDF comp.) in x_frac [-]
      sigma_x_2_n, & ! Standard deviation of ln x (2nd PDF comp.) in x_frac [-]
      mixt_frac,   & ! Mixture fraction                                     [-]
      x_frac_1,    & ! Fraction: x distributed lognormally (1st PDF comp.)  [-]
      x_frac_2,    & ! Fraction: x distributed lognormally (2nd PDF comp.)  [-]
      x_mean         ! Overall mean value of x                              [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      xp2            ! Overall variance of x                                [-]


    ! Calculate overall variance of x, <x'^2>.
    if ( abs(sigma_x_1) < eps .and. abs(sigma_x_2) < eps ) then

       ! The value of x is constant within both PDF components.
       xp2 = ( mixt_frac * x_frac_1 * mu_x_1**2 &
             + ( one - mixt_frac ) * x_frac_2 * mu_x_2**2 &
             ) &
             - x_mean**2


    elseif ( abs(sigma_x_1) < eps ) then

       ! The value of x is constant within the 1st PDF component.
       xp2 = ( mixt_frac * x_frac_1 * mu_x_1**2 &
             + ( one - mixt_frac ) * x_frac_2 &
               * exp( two * mu_x_2_n + two * sigma_x_2_n**2 ) &
             ) &
             - x_mean**2


    elseif ( abs(sigma_x_2) < eps ) then

       ! The value of x is constant within the 2nd PDF component.
       xp2 = ( mixt_frac * x_frac_1 &
               * exp( two * mu_x_1_n + two * sigma_x_1_n**2 ) &
             + ( one - mixt_frac ) * x_frac_2 * mu_x_2**2 &
             ) &
             - x_mean**2


    else  ! sigma_x_1 and sigma_x_2 > 0

       ! The value of x varies within both PDF component.
       xp2 = ( mixt_frac * x_frac_1 &
               * exp( two * mu_x_1_n + two * sigma_x_1_n**2 ) &
             + ( one - mixt_frac ) * x_frac_2 &
               * exp( two * mu_x_2_n + two * sigma_x_2_n**2 ) &
             ) &
             - x_mean**2


    endif


    ! As a check, prevent negative values for hydrometeor variances due to
    ! numerical loss of precision error.
    if ( xp2 < zero ) then
       xp2 = zero
    endif


    return

  end function calc_xp2

!===============================================================================

end module pdf_utilities
