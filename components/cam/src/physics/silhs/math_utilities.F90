!-----------------------------------------------------------------------
!$Id: math_utilities.F90 7467 2015-01-28 00:32:39Z raut@uwm.edu $
!===============================================================================
module math_utilities
!-----------------------------------------------------------------------
! Various mathematical utilities
!-----------------------------------------------------------------------
  implicit none

  public :: compute_sample_mean, compute_sample_variance, compute_sample_covariance, &
            rand_integer_in_range

  private

  contains

!-----------------------------------------------------------------------
  pure function compute_sample_mean( n_levels, n_samples, weight, x_sample ) &
    result( mean )
! Description:
!   Find the mean of a set of sample points

! References:
!   None
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: real, sum

    ! Input Varibles
    integer, intent(in) :: &
      n_levels, &
      n_samples

    real( kind = core_rknd ), dimension(n_samples), intent(in) :: &
      weight   ! Weights for individual points of the vector

    real( kind = core_rknd ),dimension(n_levels,n_samples), intent(in) :: &
      x_sample ! Collection of sample points    [units vary]

    ! Return type
    real( kind = core_rknd ), dimension(n_levels) :: mean

    integer :: k

    ! ---- Begin Code ----

    ! Get rid of an annoying compiler warning.
    k = 1
    k = k

    forall( k = 1:n_levels )
      mean(k) = sum( weight(1:n_samples) * x_sample(k,1:n_samples) ) &
              / real( n_samples, kind=core_rknd )
    end forall


    return

  end function compute_sample_mean

!-----------------------------------------------------------------------
  pure function compute_sample_variance( n_levels, n_samples, x_sample, weight, x_mean ) &
    result( variance )

! Description:
!   Compute the variance of a set of sample points

! References:
!   None
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      n_levels, & ! Number of sample levels in the mean / variance
      n_samples   ! Number of sample points compute the variance of

    real( kind = core_rknd ),dimension(n_levels,n_samples), intent(in) :: &
      x_sample ! Collection of sample points    [units vary]

    real( kind = core_rknd ),dimension(n_samples), intent(in) :: &
      weight ! Coefficient to weight the nth sample point by [-]

    real( kind = core_rknd ),dimension(n_levels), intent(in) :: &
      x_mean ! Mean sample points [units vary]

    ! Output Variable
    real( kind = core_rknd ),dimension(n_levels) :: &
      variance ! Variance of x [(units vary)^2]

    ! Local Variable(s)
    integer :: sample ! Loop iterator

    ! ---- Begin Code ----

    variance(1:n_levels) = 0.0_core_rknd

    do sample=1, n_samples
      variance(1:n_levels) = variance(1:n_levels) &
        + weight(sample) * ( x_sample(1:n_levels,sample) - x_mean(1:n_levels) )**2
    end do

    variance(1:n_levels) = variance(1:n_levels) / real( n_samples, kind=core_rknd )

    return
  end function compute_sample_variance

!-----------------------------------------------------------------------
  pure function compute_sample_covariance( n_levels, n_samples, weight, &
                   x_sample, x_mean, y_sample, y_mean ) &
    result( covariance )

! Description:
!   Compute the covariance of a set of sample points of 2 variables
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      n_levels, & ! Number of sample levels in the mean / variance
      n_samples   ! Number of sample points compute the variance of

    real( kind = core_rknd ),dimension(n_levels,n_samples), intent(in) :: &
      x_sample, & ! Collection of sample points    [units vary]
      y_sample

    real( kind = core_rknd ),dimension(n_samples), intent(in) :: &
      weight ! Coefficient to weight the nth sample point by [-]

    real( kind = core_rknd ),dimension(n_levels), intent(in) :: &
      x_mean, & ! Mean sample points [units vary]
      y_mean

    ! Output Variable
    real( kind = core_rknd ),dimension(n_levels) :: &
      covariance ! Coariance of x and y [(units vary)^2]

    ! Local Variable(s)
    integer :: sample ! Loop iterator

    ! ---- Begin Code ----

    covariance(1:n_levels) = 0.0_core_rknd

    do sample=1, n_samples
      covariance(1:n_levels) = covariance(1:n_levels) &
        + weight(sample) * ( x_sample(1:n_levels,sample) - x_mean(1:n_levels) ) &
           * ( y_sample(1:n_levels,sample) - y_mean(1:n_levels) )
    end do

    covariance(1:n_levels) = covariance(1:n_levels) / real( n_samples, kind=core_rknd )

    return
  end function compute_sample_covariance

  !-----------------------------------------------------------------------
  function rand_integer_in_range(low, high)

  ! Description:
  !   Returns a uniformly distributed integer in the range [low,high]
  !   using the Mersenne Twister PRNG library.
  !
  !   The integers returned from this function are actually not quite
  !   evenly distributed because of the use of MOD. Smaller numbers are
  !   slightly more likely than larger ones. This could be fixed someday.

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use mt95, only: &
      genrand_intg, & ! Constant
      genrand_int32   ! Procedure

    implicit none

    ! Local Constants

    ! Input Variables
    integer, intent(in) :: &
      low,   &      ! Lowest possible returned value
      high          ! Highest possible returned value

    ! Output Variable
    integer :: &
      rand_integer_in_range  ! Random integer in range [low,high]

    ! Local Variables
    integer( kind = genrand_intg ) :: &
      rand_32                ! Random integer in range[-2^31, +2^31-1]

    integer :: &
      range_width

  !-----------------------------------------------------------------------
    !----- Begin Code -----

    range_width = high - low + 1
    call genrand_int32( rand_32 )
    rand_integer_in_range = abs( mod( rand_32, range_width ) ) + low

    return
  end function rand_integer_in_range
  !-----------------------------------------------------------------------


end module math_utilities
