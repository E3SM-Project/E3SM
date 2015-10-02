! $Id: surface_varnce_module.F90 5623 2012-01-17 17:55:26Z connork@uwm.edu $
!===============================================================================
module surface_varnce_module

  implicit none

  private ! Default to private

  public :: surface_varnce

  contains

!=============================================================================
  subroutine surface_varnce( upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc, & 
                             um_sfc, vm_sfc, wpsclrp_sfc,  & 
                             wp2_sfc, up2_sfc, vp2_sfc, & 
                             thlp2_sfc, rtp2_sfc, rtpthlp_sfc, err_code, & 
                             sclrp2_sfc, & 
                             sclrprtp_sfc,  & 
                             sclrpthlp_sfc )

! Description:
!   This subroutine computes estimate of the surface thermodynamic
!   second order moments.

! References:
!   None
!-------------------------------------------------------------------------------

    use parameters_model, only:  & 
      T0 ! Variable(s)

    use constants_clubb, only: & 
      grav,  & ! Variable(s)
      eps,   &
      fstderr

    use parameters_model, only: & 
      sclr_dim  ! Variable(s)

    use numerical_check, only: & 
      surface_varnce_check ! Procedure

    use error_code, only:  & 
      clubb_var_equals_NaN,  & ! Variable(s)
      clubb_at_least_debug_level, &
      clubb_no_error  ! Constant

    use array_index, only: &
      iisclr_rt, & ! Index for a scalar emulating rt
      iisclr_thl   ! Index for a scalar emulating thetal

    use stats_type, only: & 
      stat_end_update_pt, & ! Procedure(s)
      stat_update_var_pt

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: sqrt, max

    ! Constant Parameters

    ! Logical for Andre et al., 1978 parameterization.
    logical, parameter :: l_andre_1978 = .false.

    real( kind = core_rknd ), parameter ::  & 
      a_const = 1.8_core_rknd, & 
      z_const = 1.0_core_rknd, & 
      ! Vince Larson increased ufmin to stabilize arm_97.  24 Jul 2007
!      ufmin = 0.0001_core_rknd, &
      ufmin = 0.01_core_rknd, & 
      ! End Vince Larson's change.
      ! Vince Larson changed in order to make correlations between [-1,1].  31 Jan 2008.
!      sclr_var_coef = 0.25_core_rknd, & ! This value is made up! - Vince Larson 12 Jul 2005
      sclr_var_coef = 0.4_core_rknd,  & ! This value is made up! - Vince Larson 12 Jul 2005
      ! End Vince Larson's change
      ! Vince Larson reduced surface spike in scalar variances associated
      ! w/ Andre et al. 1978 scheme
      reduce_coef   = 0.2_core_rknd

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      upwp_sfc,     & ! Surface u momentum flux   [m^2/s^2]
      vpwp_sfc,     & ! Surface v momentum flux   [m^2/s^2]
      wpthlp_sfc,   & ! Surface thetal flux       [K m/s]
      wprtp_sfc,    & ! Surface moisture flux     [kg/kg m/s]
      um_sfc,       & ! Surface u wind component  [m/s]
      vm_sfc          ! Surface v wind component  [m/s]

    real( kind = core_rknd ), intent(in), dimension(sclr_dim) ::  & 
      wpsclrp_sfc ! Passive scalar flux       [units m/s]

    ! Output Variables
    real( kind = core_rknd ), intent(out) ::  & 
      wp2_sfc,     & ! Vertical velocity variance        [m^2/s^2]
      up2_sfc,     & ! u'^2                              [m^2/s^2]
      vp2_sfc,     & ! v'^2                              [m^2/s^2]
      thlp2_sfc,   & ! thetal variance                   [K^2]
      rtp2_sfc,    & ! rt variance                       [(kg/kg)^2]
      rtpthlp_sfc    ! thetal rt covariance              [kg K/kg]

    integer, intent(out) :: & 
      err_code

    real( kind = core_rknd ), intent(out), dimension(sclr_dim) ::  & 
      sclrp2_sfc,    & ! Passive scalar variance                 [units^2]
      sclrprtp_sfc,  & ! Passive scalar r_t covariance           [units kg/kg]
      sclrpthlp_sfc    ! Passive scalar theta_l covariance       [units K]

    ! Local Variables
    real( kind = core_rknd ) :: ustar2, wstar
    real( kind = core_rknd ) :: uf

    ! Variables for Andre et al., 1978 parameterization.
    real( kind = core_rknd ) :: &
      um_sfc_sqd, & ! Surface value of <u>^2                           [m^2/s^2]
      vm_sfc_sqd, & ! Surface value of <v>^2                           [m^2/s^2]
      usp2_sfc,   & ! u_s (vector oriented w/ mean sfc. wind) variance [m^2/s^2]
      vsp2_sfc      ! v_s (vector perpen. to mean sfc. wind) variance  [m^2/s^2]

    real( kind = core_rknd ) :: ustar
    real( kind = core_rknd ) :: Lngth
    real( kind = core_rknd ) :: zeta

    integer :: i ! Loop index

    ! ---- Begin Code ----

    err_code = clubb_no_error

    if ( l_andre_1978 ) then

      ! Calculate <u>^2 and <v>^2.
      um_sfc_sqd = um_sfc**2
      vm_sfc_sqd = vm_sfc**2

      ! Calculate surface friction velocity, u*.
      ustar = MAX( ( upwp_sfc**2 + vpwp_sfc**2 )**(1.0_core_rknd/4.0_core_rknd), ufmin )

      ! Find Monin-Obukhov Length (Andre et al., 1978, p. 1866).
      Lngth = - ( ustar**3 ) /  & 
                ( 0.35_core_rknd * (1.0_core_rknd/T0) * grav * wpthlp_sfc ) ! Known magic number

      ! Find the value of dimensionless height zeta
      ! (Andre et al., 1978, p. 1866).
      zeta = z_const / Lngth

      ! Andre et al, 1978, Eq. 29.
      ! Notes:  1) "reduce_coef" is a reduction coefficient intended to make
      !            the values of rtp2, thlp2, and rtpthlp smaller at the
      !            surface.
      !         2) With the reduction coefficient having a value of 0.2, the
      !            surface correlations of both w & rt and w & thl have a value
      !            of about 0.845.  These correlations are greater if zeta < 0.
      !            The correlations have a value greater than 1 if
      !            zeta <= -0.212.
      !         3) The surface correlation of rt & thl is 1.
      ! Brian Griffin; February 2, 2008.
      if ( zeta < 0.0_core_rknd ) then
        thlp2_sfc   = reduce_coef  & 
                     * ( wpthlp_sfc**2 / ustar**2 ) & 
                     * 4.0_core_rknd * ( 1.0_core_rknd - 8.3_core_rknd*zeta )**&
                     (-2.0_core_rknd/3.0_core_rknd) ! Known magic number
        rtp2_sfc    = reduce_coef  & 
                     * ( wprtp_sfc**2 / ustar**2 ) & 
                     * 4.0_core_rknd * ( 1.0_core_rknd - 8.3_core_rknd*zeta )**&
                     (-2.0_core_rknd/3.0_core_rknd) ! Known magic number
        rtpthlp_sfc = reduce_coef  & 
                     * ( wprtp_sfc*wpthlp_sfc / ustar**2 ) & 
                     * 4.0_core_rknd * ( 1.0_core_rknd - 8.3_core_rknd*zeta )**&
                     (-2.0_core_rknd/3.0_core_rknd) ! Known magic number
        wp2_sfc     = ( ustar**2 ) & 
                     * ( 1.75_core_rknd + 2.0_core_rknd*(-zeta)**&
                     (2.0_core_rknd/3.0_core_rknd) ) ! Known magic number
      else
        thlp2_sfc   = reduce_coef  & 
                     * 4.0_core_rknd * ( wpthlp_sfc**2 / ustar**2 ) ! Known magic number
        rtp2_sfc    = reduce_coef  & 
                     * 4.0_core_rknd * ( wprtp_sfc**2 / ustar**2 ) ! Known magic number
        rtpthlp_sfc = reduce_coef  & 
                     * 4.0_core_rknd * ( wprtp_sfc*wpthlp_sfc / ustar**2 ) ! Known magic number
        wp2_sfc     = 1.75_core_rknd * ustar**2 ! Known magic number
      end if

      ! Calculate wstar following Andre et al., 1978, p. 1866.
      wstar = ( (1.0_core_rknd/T0) * grav * wpthlp_sfc * z_const )**(1.0_core_rknd/3.0_core_rknd)

      ! Andre et al., 1978, Eq. 29.
      ! Andre et al. (1978) defines horizontal wind surface variances in terms
      ! of orientation with the mean surface wind.  The vector u_s is the wind
      ! vector oriented with the mean surface wind.  The vector v_s is the wind
      ! vector oriented perpendicular to the mean surface wind.  Thus, <u_s> is
      ! equal to the mean surface wind (both in speed and direction), and <v_s>
      ! is 0.  Equation 29 gives the formula for the variance of u_s, which is
      ! <u_s'^2> (usp2_sfc in the code), and the formula for the variance of
      ! v_s, which is <v_s'^2> (vsp2_sfc in the code).
      if ( wpthlp_sfc > 0.0_core_rknd ) then
        usp2_sfc = 4.0_core_rknd * ustar**2 + 0.3_core_rknd * wstar**2 ! Known magic number
        vsp2_sfc = 1.75_core_rknd * ustar**2 + 0.3_core_rknd * wstar**2 ! Known magic number
      else
        usp2_sfc = 4.0_core_rknd * ustar**2 ! Known magic number
        vsp2_sfc = 1.75_core_rknd * ustar**2 ! Known magic number
      end if

      ! Variance of u, <u'^2>, at the surface can be found from <u_s'^2>,
      ! <v_s'^2>, and mean winds (at the surface) <u> and <v>, such that:
      !    <u'^2>|_sfc = <u_s'^2> * [ <u>^2 / ( <u>^2 + <v>^2 ) ]
      !                  + <v_s'^2> * [ <v>^2 / ( <u>^2 + <v>^2 ) ];
      ! where <u>^2 + <v>^2 /= 0.
      up2_sfc  &
         = usp2_sfc * ( um_sfc_sqd / max( um_sfc_sqd + vm_sfc_sqd , eps ) )  &
           + vsp2_sfc * ( vm_sfc_sqd / max( um_sfc_sqd + vm_sfc_sqd , eps ) )

      ! Variance of v, <v'^2>, at the surface can be found from <u_s'^2>,
      ! <v_s'^2>, and mean winds (at the surface) <u> and <v>, such that:
      !    <v'^2>|_sfc = <v_s'^2> * [ <u>^2 / ( <u>^2 + <v>^2 ) ]
      !                  + <u_s'^2> * [ <v>^2 / ( <u>^2 + <v>^2 ) ];
      ! where <u>^2 + <v>^2 /= 0.
      vp2_sfc  &
         = vsp2_sfc * ( um_sfc_sqd / max( um_sfc_sqd + vm_sfc_sqd , eps ) )  &
           + usp2_sfc * ( vm_sfc_sqd / max( um_sfc_sqd + vm_sfc_sqd , eps ) )

      ! Passive scalars
      if ( sclr_dim > 0 ) then
        do i = 1, sclr_dim
          ! Notes:  1) "reduce_coef" is a reduction coefficient intended to
          !            make the values of sclrprtp, sclrpthlp, and sclrp2
          !            smaller at the surface.
          !         2) With the reduction coefficient having a value of 0.2,
          !            the surface correlation of w & sclr has a value of
          !            about 0.845.  The correlation is greater if zeta < 0.
          !            The correlation has a value greater than 1 if
          !            zeta <= -0.212.
          !         3) The surface correlations of both rt & sclr and
          !            thl & sclr are 1.
          ! Brian Griffin; February 2, 2008.
          if ( zeta < 0.0_core_rknd ) then
            sclrprtp_sfc(i)  & 
            = reduce_coef  & 
             * ( wpsclrp_sfc(i)*wprtp_sfc / ustar**2 ) & 
             * 4.0_core_rknd * ( 1.0_core_rknd - 8.3_core_rknd*zeta )**&
             (-2.0_core_rknd/3.0_core_rknd) ! Known magic number
            sclrpthlp_sfc(i)  & 
            = reduce_coef  & 
             * ( wpsclrp_sfc(i)*wpthlp_sfc / ustar**2 ) & 
             * 4.0_core_rknd * ( 1.0_core_rknd - 8.3_core_rknd*zeta )**&
             (-2.0_core_rknd/3.0_core_rknd) ! Known magic number
            sclrp2_sfc(i)  & 
            = reduce_coef   & 
             * ( wpsclrp_sfc(i)**2 / ustar**2 ) & 
             * 4.0_core_rknd * ( 1.0_core_rknd - 8.3_core_rknd*zeta )**&
             (-2.0_core_rknd/3.0_core_rknd) ! Known magic number
          else
            sclrprtp_sfc(i)  & 
            = reduce_coef  & 
             * 4.0_core_rknd * ( wpsclrp_sfc(i)*wprtp_sfc / ustar**2 ) ! Known magic number
            sclrpthlp_sfc(i)  & 
            = reduce_coef  & 
             * 4.0_core_rknd * ( wpsclrp_sfc(i)*wpthlp_sfc / ustar**2 ) ! Known magic number
            sclrp2_sfc(i)  & 
            = reduce_coef & 
             * 4.0_core_rknd * ( wpsclrp_sfc(i)**2 / ustar**2 ) ! Known magic number
          end if
        end do ! 1,...sclr_dim
      end if

    else ! Previous code.

      ! Compute ustar^2

      ustar2 = sqrt( upwp_sfc * upwp_sfc + vpwp_sfc * vpwp_sfc )

      ! Compute wstar following Andre et al., 1976

      if ( wpthlp_sfc > 0._core_rknd ) then
        wstar = ( 1.0_core_rknd/T0 * grav * wpthlp_sfc * z_const ) ** (1._core_rknd/3._core_rknd)
      else
        wstar = 0._core_rknd
      end if

      ! Surface friction velocity following Andre et al. 1978

      uf = sqrt( ustar2 + 0.3_core_rknd * wstar * wstar ) ! Known magic number
      uf = max( ufmin, uf )

      ! Compute estimate for surface second order moments

      wp2_sfc     =  a_const * uf**2
      up2_sfc     =  2.0_core_rknd * a_const * uf**2  ! From Andre, et al. 1978
      vp2_sfc     =  2.0_core_rknd * a_const * uf**2  ! "  "
      ! Vince Larson changed to make correlations between [-1,1]  31 Jan 2008
!       thlp2_sfc   = 0.1 * a * ( wpthlp_sfc / uf )**2
!       rtp2_sfc    = 0.4 * a * ( wprtp_sfc / uf )**2
!       rtpthlp_sfc = a * ( wpthlp_sfc / uf ) * ( wprtp_sfc / uf )
      ! Notes:  1) With "a" having a value of 1.8, the surface correlations of
      !            both w & rt and w & thl have a value of about 0.878.
      !         2) The surface correlation of rt & thl is 0.5.
      ! Brian Griffin; February 2, 2008.

      thlp2_sfc = 0.4_core_rknd * a_const * ( wpthlp_sfc / uf )**2 ! Known magic number

      rtp2_sfc = 0.4_core_rknd * a_const * ( wprtp_sfc / uf )**2 ! Known magic number

      rtpthlp_sfc = 0.2_core_rknd * a_const * ( wpthlp_sfc / uf ) &
        * ( wprtp_sfc / uf )! Known magic number

      ! End Vince Larson's change.

      ! Passive scalars
      if ( sclr_dim > 0 ) then
        do i=1, sclr_dim
          ! Vince Larson changed coeffs to make correlations between [-1,1]. 31 Jan 2008
!             sclrprtp_sfc(i) &
!             = a * (wprtp_sfc / uf) * (wpsclrp_sfc(i) / uf)
!             sclrpthlp_sfc(i) &
!             = a * (wpthlp_sfc / uf) * (wpsclrp_sfc(i) / uf)
!             sclrp2_sfc(i) &
!             = sclr_var_coef * a * ( wpsclrp_sfc(i) / uf )**2
          ! Notes:  1) With "a" having a value of 1.8 and "sclr_var_coef"
          !            having a value of 0.4, the surface correlation of
          !            w & sclr has a value of about 0.878.
          !         2) With "sclr_var_coef" having a value of 0.4, the
          !            surface correlations of both rt & sclr and
          !            thl & sclr are 0.5.
          ! Brian Griffin; February 2, 2008.

          ! We use the following if..then's to make sclr_rt and sclr_thl close to
          ! the actual thlp2/rtp2 at the surface. -dschanen 25 Sep 08
          if ( i == iisclr_rt ) then
            ! If we are trying to emulate rt with the scalar, then we
            ! use the variance coefficient from above
            sclrprtp_sfc(i) = 0.4_core_rknd * a_const * (wprtp_sfc / uf) * &
               (wpsclrp_sfc(i) / uf)!Known magic number
          else
            sclrprtp_sfc(i) = 0.2_core_rknd * a_const * (wprtp_sfc / uf) * &
               (wpsclrp_sfc(i) / uf)!Known magic number
          end if

          if ( i == iisclr_thl ) then
            ! As above, but for thetal
            sclrpthlp_sfc(i) = 0.4_core_rknd * a_const * (wpthlp_sfc / uf) &
                                * (wpsclrp_sfc(i) / uf) ! Known magic number
          else
            sclrpthlp_sfc(i) = 0.2_core_rknd * a_const * (wpthlp_sfc / uf) &
                                * (wpsclrp_sfc(i) / uf) ! Known magic number
          end if

          sclrp2_sfc(i) = sclr_var_coef * a_const * ( wpsclrp_sfc(i) / uf )**2

          ! End Vince Larson's change.

        end do ! 1,...sclr_dim
      end if ! sclr_dim > 0

    end if

    if ( clubb_at_least_debug_level( 2 ) ) then

      call surface_varnce_check( wp2_sfc, up2_sfc, vp2_sfc,  & 
                                 thlp2_sfc, rtp2_sfc, rtpthlp_sfc, & 
                                 err_code, & 
                                 sclrp2_sfc, sclrprtp_sfc, sclrpthlp_sfc )

!       Error reporting
!       Joshua Fasching February 2008
      if ( err_code == clubb_var_equals_NaN ) then

        write(fstderr,*) "Error in surface_varnce"
        write(fstderr,*) "Intent(in)"

        write(fstderr,*) "upwp_sfc = ", upwp_sfc
        write(fstderr,*) "vpwp_sfc = ", vpwp_sfc
        write(fstderr,*) "wpthlp_sfc = ", wpthlp_sfc
        write(fstderr,*) "wprtp_sfc = ", wprtp_sfc

        if ( sclr_dim > 0 ) then
          write(fstderr,*) "wpsclrp_sfc = ", wpsclrp_sfc
        end if

        write(fstderr,*) "Intent(out)"

        write(fstderr,*) "wp2_sfc = ", wp2_sfc
        write(fstderr,*) "up2_sfc = ", up2_sfc
        write(fstderr,*) "vp2_sfc = ", vp2_sfc
        write(fstderr,*) "thlp2_sfc = ", thlp2_sfc
        write(fstderr,*) "rtp2_sfc = ", rtp2_sfc
        write(fstderr,*) "rtpthlp_sfc = ", rtpthlp_sfc

        if ( sclr_dim > 0 ) then
          write(fstderr,*) "sclrp2_sfc = ", sclrp2_sfc
          write(fstderr,*) "sclrprtp_sfc = ", sclrprtp_sfc
          write(fstderr,*) "sclrpthlp_sfc = ", sclrpthlp_sfc
        end if

      end if ! err_code == clubb_var_equals_NaN

    end if ! clubb_at_least_debug_level ( 2 )

    return

  end subroutine surface_varnce

!===============================================================================

end module surface_varnce_module
