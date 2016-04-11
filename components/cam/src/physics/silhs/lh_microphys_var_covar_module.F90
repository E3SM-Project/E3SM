!-------------------------------------------------------------------------------
! $Id: lh_microphys_var_covar_module.F90 7461 2015-01-17 02:44:34Z raut@uwm.edu $
!===============================================================================
module lh_microphys_var_covar_module

  implicit none

  public :: lh_microphys_var_covar_driver

  private ! Default scope

  contains

  !-----------------------------------------------------------------------
  subroutine lh_microphys_var_covar_driver &
             ( nz, num_samples, dt, lh_sample_point_weights, &
               lh_rt_all, lh_thl_all, lh_w_all, &
               lh_rcm_mc_all, lh_rvm_mc_all, lh_thlm_mc_all, &
               lh_rtp2_mc, lh_thlp2_mc, lh_wprtp_mc, &
               lh_wpthlp_mc, lh_rtpthlp_mc )

  ! Description:
  !   Computes the effect of microphysics on gridbox variances and covariances

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd ! Constant

    implicit none

    ! Input Variables!
    integer, intent(in) :: &
      nz,           &                  ! Number of vertical levels
      num_samples                      ! Number of SILHS sample points

    real( kind = core_rknd ), intent(in) :: &
      dt                               ! Model time step                             [s]

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      lh_sample_point_weights          ! Weight of SILHS sample points

    real( kind = core_rknd ), dimension(nz,num_samples), intent(in) :: &
      lh_rt_all, &                     ! SILHS samples of total water                [kg/kg]
      lh_thl_all, &                    ! SILHS samples of potential temperature      [K]
      lh_w_all, &                      ! SILHS samples of vertical velocity          [m/s]
      lh_rcm_mc_all, &                 ! SILHS microphys. tendency of rcm            [kg/kg/s]
      lh_rvm_mc_all, &                 ! SILHS microphys. tendency of rvm            [kg/kg/s]
      lh_thlm_mc_all                   ! SILHS microphys. tendency of thlm           [K/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      lh_rtp2_mc,     &                ! SILHS microphys. est. tendency of <rt'^2>   [(kg/kg)^2/s]
      lh_thlp2_mc,    &                ! SILHS microphys. est. tendency of <thl'^2>  [K^2/s]
      lh_wprtp_mc,    &                ! SILHS microphys. est. tendency of <w'rt'>   [m*(kg/kg)/s^2]
      lh_wpthlp_mc,   &                ! SILHS microphys. est. tendency of <w'thl'>  [m*K/s^2]
      lh_rtpthlp_mc                    ! SILHS microphys. est. tendency of <rt'thl'> [K*(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ), dimension(nz) :: &
      lh_rtp2_before_microphys,    & ! <rt'^2> before microphys_sub                  [(kg/kg)^2]
      lh_rtp2_after_microphys,     & ! <rt'^2> after microphys_sub                   [(kg/kg)^2]
      lh_thlp2_before_microphys,   & ! <thl'^2> before microphys_sub                 [K^2]
      lh_thlp2_after_microphys,    & ! <thl'^2> after microphys_sub                  [K^2]
      lh_wprtp_before_microphys,   & ! <w'rt'> before microphys_sub                  [m*(kg/kg)/s]
      lh_wprtp_after_microphys,    & ! <w'rt'> after microphys_sub                   [m*(kg/kg)/s]
      lh_wpthlp_before_microphys,  & ! <w'thl'> before microphys_sub                 [m*K/s]
      lh_wpthlp_after_microphys,   & ! <w'thl'> after microphys_sub                  [m*K/s]
      lh_rtpthlp_before_microphys, & ! <rt'thl'> before microphys_sub                [K*(kg/kg)/s]
      lh_rtpthlp_after_microphys     ! <rt'thl'> after microphys_sub                 [K*(kg/kg)/s]


  !-----------------------------------------------------------------------

    !----- Begin Code -----
    call lh_moments ( num_samples, lh_sample_point_weights, nz, &              ! Intent (in)
                      lh_rt_all, lh_thl_all, lh_w_all, &                       ! Intent (in)
                      lh_rtp2_before_microphys, lh_thlp2_before_microphys, &   ! Intent (out)
                      lh_wprtp_before_microphys, lh_wpthlp_before_microphys, & ! Intent (out)
                      lh_rtpthlp_before_microphys )                            ! Intent (out)

    call lh_moments ( num_samples, lh_sample_point_weights, nz, &              ! Intent (in)
                      lh_rt_all + dt * ( lh_rcm_mc_all + lh_rvm_mc_all ), &    ! Intent (in)
                      lh_thl_all + dt * lh_thlm_mc_all, &                      ! Intent (in)
                      lh_w_all, &                                              ! Intent (in)
                      lh_rtp2_after_microphys, lh_thlp2_after_microphys, &     ! Intent (out)
                      lh_wprtp_after_microphys, lh_wpthlp_after_microphys, &   ! Intent (out)
                      lh_rtpthlp_after_microphys )                             ! Intent (out)

    lh_wpthlp_mc = ( lh_wpthlp_after_microphys - lh_wpthlp_before_microphys ) / dt
    lh_wprtp_mc = ( lh_wprtp_after_microphys - lh_wprtp_before_microphys ) / dt
    lh_rtp2_mc = ( lh_rtp2_after_microphys - lh_rtp2_before_microphys ) / dt
    lh_thlp2_mc = ( lh_thlp2_after_microphys - lh_thlp2_before_microphys) / dt
    lh_rtpthlp_mc = ( lh_rtpthlp_after_microphys - lh_rtpthlp_before_microphys) / dt

    return
  end subroutine lh_microphys_var_covar_driver

  !-----------------------------------------------------------------------
  subroutine lh_moments ( n_samples, lh_weights, nz, & 
                 rt_all_samples, thl_all_samples, w_all_samples, &
                 lh_rtp2, lh_thlp2, &
                 lh_wprtp, lh_wpthlp, &
                 lh_rtpthlp )

  ! Description:
  !   Calculates variances and covariances using LH sample columns
  !
  ! References:
  !   None
  !
  ! TODO: This code assumes nz == gr%nnzp since it references zt2zm;  this is
  ! not necessarily the case.
  !-----------------------------------------------------------------------------

    use grid_class, only: &
      zt2zm    ! Procedures

    use math_utilities, only: &
      compute_sample_mean, & ! functions
      compute_sample_variance, &
      compute_sample_covariance

    use clubb_precision, only: &
      core_rknd ! Constants

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,            & ! Number of vertical levels
      n_samples        ! Number of sample columns from latin hypercube

    real( kind = core_rknd ), dimension(n_samples), intent(in) :: &
      lh_weights   ! Sample  point weights                   [-]

    real( kind = core_rknd ), dimension(nz,n_samples), intent(in) :: &
      rt_all_samples, &  ! rt columns from latin hypercube   [kg/kg]
      thl_all_samples, & ! thl columns from latin hypercube  [K]
      w_all_samples      ! w columns from latin hypercube    [m/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      lh_rtp2, &    ! Latin hypercube estimate of <rt'^2>     [(kg/kg)^2]
      lh_thlp2, &   ! Latin hypercube estimate of <thl'^2>    [K^2]
      lh_wprtp, &   ! Latin hypercube estimate of <w'rt'>     [m*(kg/kg)/s]
      lh_wpthlp, &  ! Latin hypercube estimate of <w'thl'>    [m*K/s]
      lh_rtpthlp    ! Latin hypercube estimate of <rt'thl'>   [K*(kg/kg)]

    ! Local variables
    real( kind = core_rknd ), dimension(nz) :: &  
      rt_mean, &    ! Latin hypercube estimate of rtm         [kg/kg]
      thl_mean, &   ! Latin hypercube estimate of thlm        [K]
      w_mean, &     ! Latin hypercube estimate of wm          [m/s]
      rtp2_zt, &    ! Estimate of <rt'^2> on the zt grid      [(kg/kg)^2]
      thlp2_zt, &   ! Estimate of <thl'^2> on the zt grid     [K^2]
      wprtp_zt, &   ! Estimate of <w'rt'> on the zt grid      [m*(kg/kg)/s]
      wpthlp_zt, &  ! Estimate of <w'thl'> on the zt grid     [m*K/s]
      rtpthlp_zt    ! Estimate of <rt'thl'> on the zt grid    [K*(kg/kg)]


    ! ---- Begin code ----

    rt_mean = compute_sample_mean( nz, n_samples, lh_weights, rt_all_samples )
    thl_mean = compute_sample_mean( nz, n_samples, lh_weights, thl_all_samples )
    w_mean = compute_sample_mean( nz, n_samples, lh_weights, w_all_samples )

    rtp2_zt = compute_sample_variance( nz, n_samples, rt_all_samples, lh_weights, rt_mean )
    thlp2_zt = compute_sample_variance( nz, n_samples, thl_all_samples, lh_weights, thl_mean )
  
    wprtp_zt = compute_sample_covariance( nz, n_samples, lh_weights, &
                   w_all_samples, w_mean, rt_all_samples, rt_mean ) 
    wpthlp_zt = compute_sample_covariance( nz, n_samples, lh_weights, &
                   w_all_samples, w_mean, thl_all_samples, thl_mean )
    rtpthlp_zt = compute_sample_covariance( nz, n_samples, lh_weights, &
                   rt_all_samples, rt_mean, thl_all_samples, thl_mean ) 

    lh_rtp2 = zt2zm( rtp2_zt )
    lh_thlp2 = zt2zm( thlp2_zt )
    lh_wprtp = zt2zm( wprtp_zt )
    lh_wpthlp = zt2zm( wpthlp_zt )
    lh_rtpthlp = zt2zm( rtpthlp_zt )

    return
  end subroutine lh_moments
  !-----------------------------------------------------------------------------
end module lh_microphys_var_covar_module
