!-------------------------------------------------------------------------------
!$Id: est_kessler_microphys_module.F90 7655 2015-04-29 19:22:26Z raut@uwm.edu $
!===============================================================================

module est_kessler_microphys_module

  implicit none

  public :: est_kessler_microphys

  private :: autoconv_estimate, rc_estimate

  private ! Default Scope

  contains

!------------------------------------------------------------------------

  subroutine est_kessler_microphys &
             ( nz, num_samples, d_variables, &
               X_nl_all_levs, pdf_params, rcm, cloud_frac, &
               X_mixt_comp_all_levs, lh_sample_point_weights, &
               lh_AKm, AKm, AKstd, AKstd_cld, &
               AKm_rcm, AKm_rcc, lh_rcm_avg )
! Description:
!   This subroutine computes microphysical grid box averages of the
!   Kessler autoconversion scheme, using both Latin hypercube sampling
!   and analytic integration, given a Latin Hypercube sample.
! References:
!   None
!------------------------------------------------------------------------

    use constants_clubb, only:  &
      pi,  & ! Variables(s)
      chi_tol, &
      zero_threshold, &
      zero

    use anl_erf, only:  &
      erf ! Procedure(s)

    use pdf_parameter_module, only:  &
      pdf_parameter  ! Type

    use clubb_precision, only: &
      core_rknd

    implicit none

    ! Input Variables

    integer, intent(in) :: &
      nz, &          ! Number of vertical levels
      num_samples, & ! Number of sample points
      d_variables    ! Number of variates

    real( kind = core_rknd ), dimension(nz,num_samples,d_variables), intent(in) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      cloud_frac    ! Cloud fraction           [-]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rcm          ! Liquid water mixing ratio                [kg/kg]

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params ! PDF parameters       [units vary]

    integer, dimension(nz,num_samples), intent(in) :: &
      X_mixt_comp_all_levs ! Whether we're in mixture component 1 or 2

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      lh_sample_point_weights ! Weight for cloud weighted sampling

    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      lh_AKm,    & ! Monte Carlo estimate of Kessler autoconversion [kg/kg/s]
      AKm,       & ! Exact Kessler autoconversion, AKm,             [kg/kg/s]
      AKstd,     & ! Exact standard deviation of gba Kessler        [kg/kg/s]
      AKstd_cld, & ! Exact w/in cloud std of gba Kessler            [kg/kg/s]
      AKm_rcm,   & ! Exact local gba Kessler auto based on rcm      [kg/kg/s]
      AKm_rcc      ! Exact local gba Kessler based on w/in cloud rc [kg/kg/s]

    ! For comparison, estimate kth liquid water using Monte Carlo
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      lh_rcm_avg ! LH estimate of grid box avg liquid water [kg/kg]

    ! Local Variables

    ! Level on which calculations are occuring
    integer :: level

    ! PDF parameters
    real( kind = core_rknd ) :: mixt_frac
!   real( kind = core_rknd ) :: w1, w2
!   real( kind = core_rknd ) :: sw1, sw2
!   real( kind = core_rknd ) :: thl1, thl2, sthl1, sthl2
!   real( kind = core_rknd ) :: rt1,rt2
!   real( kind = core_rknd ) :: srt1, srt2
    real( kind = core_rknd ) :: stdev_chi_1, stdev_chi_2, chi_1, chi_2
    real( kind = core_rknd ) :: cloud_frac_1, cloud_frac_2
!   real( kind = core_rknd ) :: rc1, rc2

    ! Cloud fraction 0<cloud_frac<1, mean liquid water mix ratio [kg/kg]
!   real( kind = core_rknd ) :: cloud_frac, rcm

    real( kind = core_rknd ), dimension(num_samples) :: &
      rcm_sample ! Sample points of rcm         [kg/kg]

    ! Variables needed for exact Kessler autoconversion, AKm
    real( kind = core_rknd ) :: r_crit, K_one
    real( kind = core_rknd ) :: chi_n_1_crit, cloud_frac_1_crit, chi_n_2_crit, cloud_frac_2_crit
    real( kind = core_rknd ) :: AK1, AK2

    ! Variables needed for exact std of Kessler autoconversion, AKstd
    !      and within cloud standard deviation, AKstd_cld
    real( kind = core_rknd ) :: AK1var, AK2var

    ! For comparison, compute within-cloud vertical velocity analytically.
    !real C_w_cld1, C_w_cld2, w_cld_avg

    ! ---- Begin Code ----

    ! Boundary condition
    lh_AKm(1)     = 0.0_core_rknd
    AKm(1)        = 0.0_core_rknd
    AKm_rcm(1)    = 0.0_core_rknd
    AKm_rcc(1)    = 0.0_core_rknd
    lh_rcm_avg(1) = 0.0_core_rknd
    AKstd(1)      = 0.0_core_rknd
    AKstd_cld(1)  = 0.0_core_rknd

    do level = 2, nz, 1
      ! Extract PDF parameters

      !w1         = pdf_params(level)%w1
      !w2         = pdf_params(level)%w2
      !sw1        = pdf_params(level)%sw1
      !sw2        = pdf_params(level)%sw2
      !rt1        = pdf_params(level)%rt1
      !rt2        = pdf_params(level)%rt2
      !srt1       = pdf_params(level)%srt1
      !srt2       = pdf_params(level)%srt2
      !thl1       = pdf_params(level)%thl1
      !thl2       = pdf_params(level)%thl2
      !sthl1      = pdf_params(level)%sthl1
      !sthl2      = pdf_params(level)%sthl2
      mixt_frac   = pdf_params(level)%mixt_frac
!     rc1         = pdf_params(level)%rc1
!     rc2         = pdf_params(level)%rc2
!     cloud_frac_1 = pdf_params(level)%cloud_frac_1
!     cloud_frac_2 = pdf_params(level)%cloud_frac_2
      cloud_frac_1 = 1.0_core_rknd ! For in and out of cloud sampling -dschanen 30 Jul 09
      cloud_frac_2 = 1.0_core_rknd !     "    "
      chi_1          = pdf_params(level)%chi_1
      chi_2          = pdf_params(level)%chi_2
      stdev_chi_1    = pdf_params(level)%stdev_chi_1
      stdev_chi_2    = pdf_params(level)%stdev_chi_2

      ! Compute mean cloud fraction and cloud water

!     cloud_frac = mixt_frac * cloud_frac_1 + (1-mixt_frac) * cloud_frac_2
!     rcm        = mixt_frac * rc1 + (1-mixt_frac) * rc2

      !------------------------------------------------------------------------
      ! Call Kessler autoconversion microphysics using Latin hypercube sample
      ! This acts as an interface between the boundary layer scheme
      !    and the microphysics.  
      ! Then we compute Kessler ac analytically.
      !------------------------------------------------------------------------

      ! We prognose rt-thl-w,
      !    but we set means, covariance of hydromet mixing ratio's and number
      !    concentrations to constants.

      rcm_sample(1:num_samples) = max( X_nl_all_levs(level,1:num_samples,1), zero)

      ! Call microphysics, i.e. Kessler autoconversion.
      ! A_K = (1e-3/s)*(rc-0.5kg/kg)*H(rc-0.5kg/kg)
      call autoconv_estimate &
           ( num_samples, mixt_frac, &
             cloud_frac_1, cloud_frac_2, &
             rcm_sample, & 
             !X_nl(1:n,3), X_nl(1:n,4), X_nl(1:n,5),
             X_mixt_comp_all_levs(level,:), lh_sample_point_weights, lh_AKm(level) )

      ! Compute Monte Carlo estimate of liquid for test purposes.
      call rc_estimate &
           ( num_samples, mixt_frac, cloud_frac_1, &
             cloud_frac_2, rcm_sample, &
             ! X_nl(1:n,3), X_nl(1:n,4), X_nl(1:n,5),
             X_mixt_comp_all_levs(level,: ), lh_rcm_avg(level) )

      ! A test of the scheme:
      ! Compare exact (rcm) and Monte Carlo estimates (lh_rcm_avg) of
      !     specific liq water content, rcm.
      !      print*, 'rcm=', rcm
      !       print*, 'lh_rcm_avg=', lh_rcm_avg

      ! Exact Kessler autoconversion in units of (kg/kg)/s
      !        r_crit = 0.3e-3_core_rknd
      !        r_crit = 0.7e-3_core_rknd
      r_crit            = 0.2e-3_core_rknd
      K_one             = 1.e-3_core_rknd
      chi_n_1_crit          = (chi_1-r_crit)/max( stdev_chi_1, chi_tol )
      cloud_frac_1_crit  = 0.5_core_rknd*(1._core_rknd+erf(chi_n_1_crit/sqrt(2.0_core_rknd)))
      AK1               = K_one * ( (chi_1-r_crit)*cloud_frac_1_crit  & 
                         + stdev_chi_1*exp(-0.5_core_rknd*chi_n_1_crit**2)/(sqrt(2._core_rknd*pi)) )
      chi_n_2_crit          = (chi_2-r_crit)/max( stdev_chi_2, chi_tol )
      cloud_frac_2_crit  = 0.5_core_rknd*(1._core_rknd+erf(chi_n_2_crit/sqrt(2.0_core_rknd)))
      AK2               = K_one * ( (chi_2-r_crit)*cloud_frac_2_crit  & 
                         + stdev_chi_2*exp(-0.5_core_rknd*chi_n_2_crit**2)/(sqrt(2._core_rknd*pi)) )
      AKm(level)        = mixt_frac * AK1 + (1._core_rknd-mixt_frac) * AK2

      ! Exact Kessler standard deviation in units of (kg/kg)/s
      ! For some reason, sometimes AK1var, AK2var are negative
      AK1var   = max( zero_threshold, K_one * (chi_1-r_crit) * AK1  & 
               + K_one * K_one * (stdev_chi_1**2) * cloud_frac_1_crit  & 
               - AK1**2  )
      AK2var   = max( zero_threshold, K_one * (chi_2-r_crit) * AK2  & 
               + K_one * K_one * (stdev_chi_2**2) * cloud_frac_2_crit  & 
               - AK2**2  )
      ! This formula is for a grid box average:
      AKstd(level)  = sqrt( mixt_frac * ( (AK1-AKm(level))**2 + AK1var ) & 
                  + (1._core_rknd-mixt_frac) * ( (AK2-AKm(level))**2 + AK2var ) &
                  )
      ! This formula is for a within-cloud average:
      if ( cloud_frac(level) > 0._core_rknd ) then
        AKstd_cld(level) = sqrt( max( real( zero_threshold, kind=core_rknd ),   &
                  (1._core_rknd/cloud_frac(level)) * ( mixt_frac * ( AK1**2 + AK1var ) &
                            + (1._core_rknd-mixt_frac) * ( AK2**2 + AK2var )  &
                            ) & 
                 - (AKm(level)/cloud_frac(level))**2  ) & 
                        )
      else
        AKstd_cld(level) = zero_threshold
      end if

      ! Kessler autoconversion, using grid box avg liquid, rcm, as input
      AKm_rcm(level) = K_one * max( zero_threshold, rcm(level)-r_crit )

      ! Kessler ac, using within cloud liquid, rcm/cloud_frac, as input
      ! We found that for small values of cloud_frac this formula
      ! can still produce NaN values and therefore added this
      ! threshold of 0.001 here. -dschanen 3 June 2009
      if ( cloud_frac(level) > 0.001_core_rknd ) then
        AKm_rcc(level) = cloud_frac(level) * K_one * &
                         max( zero_threshold, rcm(level)/cloud_frac(level)-r_crit )
      else
        AKm_rcc(level) = zero_threshold
      end if

    end do ! level = 2, nz

    return
  end subroutine est_kessler_microphys
!-----------------------------------------------------------------------
  subroutine autoconv_estimate( num_samples, mixt_frac, &
                              cloud_frac_1, cloud_frac_2, rc, &
                             !w, Nc, rr, &
                              X_mixt_comp_one_lev, lh_sample_point_weights, ac_m )
! Description:
!   Compute Kessler grid box avg autoconversion (kg/kg)/s.
! References:
!   None
!-----------------------------------------------------------------------

    use constants_clubb, only:  &
      fstderr, & ! Constant(s)
      g_per_kg, &
      zero, &
      one

!   use error_code, only:  &
!     clubb_at_least_debug_level  ! Procedure(s)

    use parameters_silhs, only: &
      l_lh_importance_sampling ! Variable(s)

    use clubb_precision, only: &
      core_rknd

    implicit none

    ! External
    intrinsic :: epsilon

    ! Constant parameters
    logical, parameter :: &
      l_cloud_weighted_averaging = .false.

    real(kind=core_rknd), parameter :: &
      r_crit_g_kg = 0.2_core_rknd

    ! Input Variables

    integer, intent(in) :: &
      num_samples  ! Number of calls to microphysics (normally=2)

    real( kind = core_rknd ), intent(in) :: &
      mixt_frac,               & ! Mixture fraction of Gaussians
      cloud_frac_1, cloud_frac_2   ! Cloud fraction associated w/ 1st, 2nd mixture component

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      rc !, & ! n in-cloud values of spec liq water content (when positive) [kg/kg].
!     w,  & ! n in-cloud values of vertical velocity (m/s)
!     Nc, & ! n in-cloud values of droplet number (#/mg air)
!     rr    ! n in-cloud values of specific rain content (g/kg)

    integer, dimension(num_samples), intent(in) :: &
      X_mixt_comp_one_lev ! Whether we're in the first or second mixture component

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
       lh_sample_point_weights ! Weight for cloud weighted sampling

    ! Output Variables

    ! a scalar representing grid box average autoconversion;
    ! has same units as rc; divide by total cloud fraction to obtain
    ! within-cloud autoconversion
    real( kind = core_rknd ), intent(out) :: &
      ac_m

    ! Local Variables

    integer :: sample
    integer :: n1, n2
    real( kind = core_rknd ) :: ac_m1, ac_m2
    real( kind = core_rknd ) :: coeff, r_crit

    ! ---- Begin Code ----

    ! Handle some possible errors re: proper ranges of mixt_frac, cloud_frac_1, cloud_frac_2.
    if ( mixt_frac > one .or. mixt_frac < zero ) then
      write(fstderr,*) 'Error in autoconv_estimate:  ',  &
                       'mixture fraction, mixt_frac, does not lie in [0,1].'
      write(fstderr,*) 'mixt_frac = ', mixt_frac
      stop
    end if
    if ( cloud_frac_1 > one .or. cloud_frac_1 < zero ) then
      write(fstderr,*) 'Error in autoconv_estimate:  ',  &
                       'cloud fraction 1, cloud_frac_1, does not lie in [0,1].'
      write(fstderr,*) 'cloud_frac_1 = ', cloud_frac_1
      stop
    end if
    if ( cloud_frac_2 > one .or. cloud_frac_2 < zero ) then
      write(fstderr,*) 'Error in autoconv_estimate:  ',  &
                       'cloud fraction 2, cloud_frac_2, does not lie in [0,1].'
      write(fstderr,*) 'cloud_frac_2 = ', cloud_frac_2
      stop
    end if

    ! Autoconversion formula prefactor and exponent.
    ! These are for Kessler autoconversion in (kg/kg)/s.
    coeff  = 1.e-3_core_rknd
    r_crit = r_crit_g_kg / g_per_kg

    ! Initialize autoconversion in each mixture component
    ac_m1 = zero
    ac_m2 = zero

    ! Initialize numbers of sample points corresponding
    !    to each mixture component
    n1 = 0
    n2 = 0

    do sample = 1, num_samples

      ! Choose which mixture fraction we are in.
      ! Account for cloud fraction.
      ! Follow M. E. Johnson (1987), p. 56.
!     fraction_1 = mixt_frac*cloud_frac_1 / &
!                    max( mixt_frac*cloud_frac_1+(1-mixt_frac)*cloud_frac_2, epsilon( mixt_frac ) )
!          print*, 'fraction_1= ', fraction_1

! V. Larson change to try to fix sampling
!          if ( in_mixt_frac_1( X_u_one_lev(sample,d_variables+1), fraction_1 ) ) then
!          print*, '-1+2*int((sample+1)/2)= ', -1+2*int((sample+1)/2)
!          print*, '-1+2*int((sample+1)/2)= ', int(sample)
      if ( X_mixt_comp_one_lev(sample) == 1 ) then
! End of V. Larson fix

! Use an idealized formula to compute autoconversion
!      in mixture comp. 1
! A_K = (1e-3/s)*(rc-0.5g/kg)*H(rc-0.5g/kg)
! This is the first of two lines where
!      a user must add a new microphysics scheme.
        if ( l_lh_importance_sampling ) then
          ac_m1 = ac_m1 + coeff*max(zero,rc(sample)-r_crit)&
                  * lh_sample_point_weights(sample)
        else
          ac_m1 = ac_m1 + coeff*max(zero,rc(sample)-r_crit)
        end if
        n1 = n1 + 1
      else
! Use an idealized formula to compute autoconversion
!      in mixture comp. 2
! A_K = (1e-3/s)*(rc-0.5g/kg)*H(rc-0.5g/kg)
! This is the second and last line where
!      a user must add a new microphysics scheme.

        if ( l_lh_importance_sampling ) then
          ac_m2 = ac_m2 + coeff*max(zero,rc(sample)-r_crit) &
                  * lh_sample_point_weights(sample)
        else
          ac_m2 = ac_m2 + coeff*max(zero,rc(sample)-r_crit)
        end if

        n2 = n2 + 1
      end if

      ! Loop to get new sample
    end do ! sample = 1, num_samples

!! Convert sums to averages.
!! Old code that underestimates if a plume has no sample points
!       if (n1 .eq. 0) then
!          ac_m1 = 0._core_rknd
!       else
!          ac_m1 = ac_m1/n1
!       end if

!       if (n2 .eq. 0) then
!         ac_m2 = 0._core_rknd
!       else
!          ac_m2 = ac_m2/n2
!       end if

    if ( n1 == 0 .and. n2 == 0 ) then
      stop 'Error:  no sample points in autoconv_estimate'
    end if

    if ( l_cloud_weighted_averaging ) then
      ! Convert sums to averages.
      ! If we have no sample points for a certain plume,
      !    then we estimate the plume liquid water by the
      !    other plume's value.
      if ( .not. (n1 == 0) ) then
        ac_m1 = ac_m1/ real( n1, kind=core_rknd )
      end if

      if ( .not. (n2 == 0) ) then
        ac_m2 = ac_m2/ real( n2, kind=core_rknd )
      end if

      if ( n1 == 0 ) then
        ac_m1 = ac_m2
      end if

      if ( n2 == 0 ) then
        ac_m2 = ac_m1
      end if

      ! Grid box average.
      ac_m = mixt_frac*cloud_frac_1*ac_m1 + (one-mixt_frac)*cloud_frac_2*ac_m2

    else
      ac_m = ( ac_m1 + ac_m2 ) / real( num_samples, kind=core_rknd )

    end if

!   print*, 'autoconv_estimate: acm=', ac_m

    return
  end subroutine autoconv_estimate

!----------------------------------------------------------------------
  subroutine rc_estimate( num_samples, mixt_frac, C1, C2, rc, & ! w,   & 
                         !N_pts, rr, 
                           X_mixt_comp_one_lev, rc_m )
! Description:
!   Compute an Monte Carlo estimate of grid box avg liquid water.
! References:
!   None
!---------------------------------------------------------------

    use constants_clubb, only:  &
      fstderr      ! Constant(s)

    use clubb_precision, only: &
      core_rknd    ! Core precision

    use constants_clubb, only: &
      one, &       ! Constant(s)
      zero

!   use error_code, only:  &
!       clubb_at_least_debug_level  ! Procedure(s)

    implicit none

    ! Constant parameters
    logical, parameter :: &
      l_cloud_weighted_averaging   = .false.

    ! Input Variables
    integer, intent(in) :: &
      num_samples ! Number of calls to microphysics (normally=2)

    real( kind = core_rknd ), intent(in) :: &
      mixt_frac, & ! Mixture fraction of Gaussians
      C1, C2       ! Cloud fraction associated w/ 1st, 2nd mixture component

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      rc !, & ! n in-cloud values of spec liq water content [kg/kg].
!     w,  & ! n in-cloud values of vertical velocity (m/s)
!     Npts, & ! n in-cloud values of droplet number (#/kg air)
!     rr    ! n in-cloud values of specific rain content (kg/kg)

    integer, dimension(num_samples), intent(in) :: &
      X_mixt_comp_one_lev ! Whether we're in the first or second mixture component

    ! Output Variables

    ! A scalar representing grid box avg specific liquid water;
    ! divide by total cloud fraction to obtain within-cloud liquid water
    real( kind = core_rknd ), intent(out) :: rc_m

    ! Local Variables

    integer :: sample
    integer :: n1, n2
    real( kind = core_rknd ) :: rc_m1, rc_m2
    real( kind = core_rknd ) :: coeff, expn

    ! ---- Begin Code ----

    ! Handle some possible errors re: proper ranges of mixt_frac, C1, C2.
    if ( mixt_frac > one .or. mixt_frac < zero ) then
      write(fstderr,*) 'Error in rc_estimate:  ',  &
                       'mixture fraction, mixt_frac, does not lie in [0,1].'
      stop
    end if
    if ( C1 > one .or. C1 < zero ) then
      write(fstderr,*) 'Error in rc_estimate:  ',  &
                       'cloud fraction 1, C1, does not lie in [0,1].'
      stop
    end if
    if ( C2 > one .or. C2 < zero ) then
      write(fstderr,*) 'Error in rc_estimate:  ',  &
                       'cloud fraction 2, C2, does not lie in [0,1].'
      stop
    end if

    ! To compute liquid water, need to set coeff=expn=1.
    coeff = one
    expn  = one

    ! Initialize liquid in each mixture component
    rc_m1 = zero
    rc_m2 = zero

    ! Initialize numbers of sample points corresponding
    !    to each mixture component
    n1    = 0
    n2    = 0

    do sample = 1, num_samples

      ! Choose which mixture fraction we are in.
      ! Account for cloud fraction.
      ! Follow M. E. Johnson (1987), p. 56.
!     fraction_1 = mixt_frac*C1/max( (mixt_frac*C1+(1-mixt_frac)*C2), epsilon( mixt_frac ) )
      if ( X_mixt_comp_one_lev(sample) == 1 ) then
        ! Use an idealized formula to compute liquid
        !      in mixture comp. 1
        rc_m1 = rc_m1 + coeff*(rc(sample))**expn
        n1    = n1 + 1
      else
        ! Use an idealized formula to compute liquid
        !      in mixture comp. 2
        rc_m2 = rc_m2 + coeff*(rc(sample))**expn
        n2    = n2 + 1
      end if

      ! Loop to get new sample
    end do

!! Convert sums to averages.
!! Old code that underestimates if n1 or n2 = 0.
!   if ( n1 == 0 ) then
!     rc_m1 = zero
!   else
!     rc_m1 = rc_m1/n1
!   end if

!   if ( n2 == 0 ) then
!     rc_m2 = zero
!   else
!     rc_m2 = rc_m2/n2
!   end if


    ! Convert sums to averages.
    ! If we have no sample points for a certain plume,
    !    then we estimate the plume liquid water by the
    !    other plume's value.
    if (n1 == 0 .and. n2 == 0) then
      stop 'Error:  no sample points in rc_estimate'
    end if

    if ( l_cloud_weighted_averaging ) then
      if ( .not. (n1 == 0) ) then
        rc_m1 = rc_m1/real( n1, kind=core_rknd )
      end if

      if ( .not. (n2 == 0) ) then
        rc_m2 = rc_m2/real( n2, kind=core_rknd )
      end if

      if (n1 == 0) then
        rc_m1 = rc_m2
      end if

      if (n2 == 0) then
        rc_m2 = rc_m1
      end if
      ! Grid box average.
      rc_m = mixt_frac*C1*rc_m1 + (one-mixt_frac)*C2*rc_m2

    end if ! l_cloud_weighted_averaging

    ! Grid box average.
    rc_m = ( rc_m1 + rc_m2 ) / real(num_samples, kind = core_rknd)

    return
  end subroutine rc_estimate
!---------------------------------------------------------------

end module est_kessler_microphys_module
