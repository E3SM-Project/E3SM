! $Id$
!===============================================================================
module adg1_adg2_3d_luhar_pdf

  ! Description:
  ! Contains code to calculate the PDF parameters using the Analytic Double
  ! Gaussian 1 (ADG1) PDF closure, the Analytic Double Gaussian 2 (ADG2) PDF
  ! closure, or the 3D Luhar PDF closure.  The ADG1 PDF is the original PDF used
  ! by CLUBB in interactive runs.

  ! References:
  !-------------------------------------------------------------------------

  implicit none

  public :: ADG1_pdf_driver,     & ! Procedure(s)
            ADG2_pdf_driver,     &
            Luhar_3D_pdf_driver, &
            ADG1_w_closure,      &
            calc_Luhar_params,   &
            close_Luhar_pdf

  private :: ADG1_ADG2_responder_params, & ! Procedure(s)
             backsolve_Luhar_params,     &
             max_cubic_root

  private

  contains

  !=============================================================================
  subroutine ADG1_pdf_driver( wm, rtm, thlm, um, vm,                   & ! In
                              wp2, rtp2, thlp2, up2, vp2,              & ! In
                              Skw, wprtp, wpthlp, upwp, vpwp, sqrt_wp2,& ! In
                              sigma_sqd_w, mixt_frac_max_mag,          & ! In
                              sclrm, sclrp2, wpsclrp, l_scalar_calc,   & ! In
                              w_1, w_2, rt_1, rt_2, thl_1, thl_2,      & ! Out
                              u_1, u_2, v_1, v_2,                      & ! Out
                              varnce_w_1, varnce_w_2, varnce_rt_1,     & ! Out
                              varnce_rt_2, varnce_thl_1, varnce_thl_2, & ! Out
                              varnce_u_1, varnce_u_2,                  & ! Out
                              varnce_v_1, varnce_v_2,                  & ! Out
                              mixt_frac, alpha_rt, alpha_thl,          & ! Out
                              alpha_u, alpha_v,                        & ! Out
                              sclr_1, sclr_2, varnce_sclr_1,           & ! Out
                              varnce_sclr_2, alpha_sclr )                ! Out

    ! Calculates the PDF parameters using the Analytic Double Gaussian 1 (ADG1)
    ! PDF closure.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        rt_tol,  & ! Constant(s)
        thl_tol

    use parameters_model, only: &
        sclr_dim, & ! Variable(s)
        sclr_tol

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  &
      wm,          & ! Mean of w-wind comp. (vert. vel.)    [m/s] 
      rtm,         & ! Mean of total water mixing ratio     [kg/kg]
      thlm,        & ! Mean of liquid water potential temp. [K]
      um,          & ! Mean of eastward wind                [m/s]
      vm,          & ! Mean of northward wind               [m/s]
      wp2,         & ! Variance of w (overall)              [m^2/s^2] 
      rtp2,        & ! Variance of r_t (overall)            [(kg/kg)^2]
      thlp2,       & ! Variance of th_l (overall)           [K^2]
      up2,         & ! Variance of eastward wind (overall)  [m^2/s^2]
      vp2,         & ! Variance of northward wind (overall) [m^2/s^2]
      Skw,         & ! Skewness of w                        [-]
      wprtp,       & ! Covariance of w and r_t              [(kg/kg)(m/s)]
      wpthlp,      & ! Covariance of w and th_l             [K(m/s)]
      upwp,        & ! Covariance of u and w                [m^2/s^2]
      vpwp,        & ! Covariance of v and w                [m^2/s^2]
      sqrt_wp2,    & ! Square root of variance of w         [m/s]
      sigma_sqd_w    ! Width of individual w plumes         [-]

    real( kind = core_rknd ), intent(in) ::  &
      mixt_frac_max_mag    ! Maximum allowable mag. of mixt_frac  [-]

    real( kind = core_rknd ), dimension(gr%nz, sclr_dim), intent(in) ::  &
      sclrm,   & ! Mean of passive scalar (overall)           [units vary]
      sclrp2,  & ! Variance of passive scalar (overall)       [(units vary)^2]
      wpsclrp    ! Covariance of w and passive scalar         [m/s (units vary)]

    logical, intent(in) :: &
      l_scalar_calc    ! Flag to perform calculations for passive scalars

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) ::  &
      w_1,          & ! Mean of w (1st PDF component)                      [m/s]
      w_2,          & ! Mean of w (2nd PDF component)                      [m/s]
      rt_1,         & ! Mean of r_t (1st PDF component)                  [kg/kg]
      rt_2,         & ! Mean of r_t (2nd PDF component)                  [kg/kg]
      thl_1,        & ! Mean of th_l (1st PDF component)                     [K]
      thl_2,        & ! Mean of th_l (2nd PDF component)                     [K]
      u_1,          & ! Mean of u (eastward wind) (1st PDF component)      [m/s]
      u_2,          & ! Mean of u (eastward wind) (2nd PDF component)      [m/s]
      v_1,          & ! Mean of v (northward wind) (1st PDF component)     [m/s]
      v_2,          & ! Mean of v (northward wind) (2nd PDF component)     [m/s]
      varnce_w_1,   & ! Variance of w (1st PDF component)              [m^2/s^2]
      varnce_w_2,   & ! Variance of w (2nd PDF component)              [m^2/s^2]
      varnce_rt_1,  & ! Variance of r_t (1st PDF component)          [kg^2/kg^2]
      varnce_rt_2,  & ! Variance of r_t (2nd PDF component)          [kg^2/kg^2]
      varnce_thl_1, & ! Variance of th_l (1st PDF component)               [K^2]
      varnce_thl_2, & ! Variance of th_l (2nd PDF component)               [K^2]
      varnce_u_1,   & ! Variance of u wind (1st PDF component)         [m^2/s^2]
      varnce_u_2,   & ! Variance of u wind (2nd PDF component)         [m^2/s^2]
      varnce_v_1,   & ! Variance of v wind (1st PDF component)         [m^2/s^2]
      varnce_v_2,   & ! Variance of v wind (2nd PDF component)         [m^2/s^2]
      mixt_frac,    & ! Mixture fraction (weight of 1st PDF component)       [-]
      alpha_thl,    & ! Factor relating to normalized variance for th_l      [-]
      alpha_rt,     & ! Factor relating to normalized variance for r_t       [-]
      alpha_u,      & ! Factor relating to normalized variance for u wind    [-]
      alpha_v         ! Factor relating to normalized variance for v wind    [-]

    real( kind = core_rknd ), dimension(gr%nz, sclr_dim), intent(out) ::  &
      sclr_1,        & ! Mean of passive scalar (1st PDF component) [units vary]
      sclr_2,        & ! Mean of passive scalar (2nd PDF component) [units vary]
      varnce_sclr_1, & ! Variance of pass. sclr (1st PDF comp.) [(units vary)^2]
      varnce_sclr_2, & ! Variance of pass. sclr (2nd PDF comp.) [(units vary)^2]
      alpha_sclr       ! Factor relating to normalized variance for sclr     [-]

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) ::  &
      w_1_n, & ! Normalized mean of w (1st PDF component)     [-]
      w_2_n    ! Normalized mean of w (2nd PDF component)     [-]

    integer :: i  ! Loop index


    ! Calculate the mixture fraction and the PDF component means and variances
    ! of w.
    call ADG1_w_closure( wm, wp2, Skw, sigma_sqd_w, &        ! In
                         sqrt_wp2, mixt_frac_max_mag, &      ! In
                         w_1, w_2, w_1_n, w_2_n, &           ! Out
                         varnce_w_1, varnce_w_2, mixt_frac ) ! Out

    ! Calculate the PDF component means and variances of rt.
    call ADG1_ADG2_responder_params( rtm, rtp2, wp2, sqrt_wp2, &       ! In
                                     wprtp, w_1_n, w_2_n, mixt_frac, & ! In
                                     sigma_sqd_w, rt_tol, &            ! In
                                     rt_1, rt_2, varnce_rt_1, &        ! Out
                                     varnce_rt_2, alpha_rt )           ! Out

    ! Calculate the PDF component means and variances of thl.
    call ADG1_ADG2_responder_params( thlm, thlp2, wp2, sqrt_wp2, &      ! In
                                     wpthlp, w_1_n, w_2_n, mixt_frac, & ! In
                                     sigma_sqd_w, thl_tol, &            ! In
                                     thl_1, thl_2, varnce_thl_1, &      ! Out
                                     varnce_thl_2, alpha_thl )          ! Out

    ! Calculate the PDF component means and variances of u wind.
    call ADG1_ADG2_responder_params( um, up2, wp2, sqrt_wp2, &        ! In
                                     upwp, w_1_n, w_2_n, mixt_frac, & ! In
                                     sigma_sqd_w, thl_tol, &          ! In
                                     u_1, u_2, varnce_u_1, &          ! Out
                                     varnce_u_2, alpha_u )            ! Out

    ! Calculate the PDF component means and variances of v wind.
    call ADG1_ADG2_responder_params( vm, vp2, wp2, sqrt_wp2, &        ! In
                                     vpwp, w_1_n, w_2_n, mixt_frac, & ! In
                                     sigma_sqd_w, thl_tol, &          ! In
                                     v_1, v_2, varnce_v_1, &          ! Out
                                     varnce_v_2, alpha_v )            ! Out

    ! Calculate the PDF component means and variances of passive scalars.
    if ( l_scalar_calc ) then
       do i = 1, sclr_dim
          call ADG1_ADG2_responder_params( sclrm(:,i), sclrp2(:,i),     & ! In
                                           wp2, sqrt_wp2, wpsclrp(:,i), & ! In
                                           w_1_n, w_2_n, mixt_frac,     & ! In
                                           sigma_sqd_w, sclr_tol(i),    & ! In
                                           sclr_1(:,i), sclr_2(:,i),    & ! Out
                                           varnce_sclr_1(:,i),          & ! Out
                                           varnce_sclr_2(:,i),          & ! Out
                                           alpha_sclr(:,i) )              ! Out
       enddo ! i=1, sclr_dim
    endif ! l_scalar_calc


    return

  end subroutine ADG1_pdf_driver

  !=============================================================================
  subroutine ADG2_pdf_driver( wm, rtm, thlm, wp2, rtp2, thlp2,         & ! In
                              Skw, wprtp, wpthlp, sqrt_wp2,            & ! In
                              sclrm, sclrp2, wpsclrp, l_scalar_calc,   & ! In
                              w_1, w_2, rt_1, rt_2, thl_1, thl_2,      & ! Out
                              varnce_w_1, varnce_w_2, varnce_rt_1,     & ! Out
                              varnce_rt_2, varnce_thl_1, varnce_thl_2, & ! Out
                              mixt_frac, alpha_rt, alpha_thl,          & ! Out
                              sigma_sqd_w, sclr_1, sclr_2,             & ! Out
                              varnce_sclr_1, varnce_sclr_2, alpha_sclr ) ! Out

    ! Calculates the PDF parameters using the Analytic Double Gaussian 2 (ADG2)
    ! PDF closure.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        one,       & ! Constant(s)
        w_tol_sqd, &
        rt_tol,    &
        thl_tol

    use parameters_model, only: &
        sclr_dim, & ! Variable(s)
        sclr_tol

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      wm,       & ! Mean of w-wind component (vertical velocity)  [m/s] 
      rtm,      & ! Mean of total water mixing ratio              [kg/kg]
      thlm,     & ! Mean of liquid water potential temperature    [K]
      wp2,      & ! Variance of w (overall)                       [m^2/s^2] 
      rtp2,     & ! Variance of r_t (overall)                     [(kg/kg)^2]
      thlp2,    & ! Variance of th_l (overall)                    [K^2]
      Skw,      & ! Skewness of w                                 [-]
      wprtp,    & ! Covariance of w and r_t                       [(kg/kg)(m/s)]
      wpthlp,   & ! Covariance of w and th_l                      [K(m/s)]
      sqrt_wp2    ! Square root of variance of w                  [m/s]

    real( kind = core_rknd ), dimension(gr%nz, sclr_dim), intent(in) ::  &
      sclrm,   & ! Mean of passive scalar (overall)           [units vary]
      sclrp2,  & ! Variance of passive scalar (overall)       [(units vary)^2]
      wpsclrp    ! Covariance of w and passive scalar         [m/s (units vary)]

    logical, intent(in) :: &
      l_scalar_calc    ! Flag to perform calculations for passive scalars

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) ::  &
      w_1,          & ! Mean of w (1st PDF component)                      [m/s]
      w_2,          & ! Mean of w (2nd PDF component)                      [m/s]
      rt_1,         & ! Mean of r_t (1st PDF component)                  [kg/kg]
      rt_2,         & ! Mean of r_t (2nd PDF component)                  [kg/kg]
      thl_1,        & ! Mean of th_l (1st PDF component)                     [K]
      thl_2,        & ! Mean of th_l (2nd PDF component)                     [K]
      varnce_w_1,   & ! Variance of w (1st PDF component)              [m^2/s^2]
      varnce_w_2,   & ! Variance of w (2nd PDF component)              [m^2/s^2]
      varnce_rt_1,  & ! Variance of r_t (1st PDF component)          [kg^2/kg^2]
      varnce_rt_2,  & ! Variance of r_t (2nd PDF component)          [kg^2/kg^2]
      varnce_thl_1, & ! Variance of th_l (1st PDF component)               [K^2]
      varnce_thl_2, & ! Variance of th_l (2nd PDF component)               [K^2]
      mixt_frac,    & ! Mixture fraction (weight of 1st PDF component)       [-]
      alpha_thl,    & ! Factor relating to normalized variance for th_l      [-]
      alpha_rt,     & ! Factor relating to normalized variance for r_t       [-]
      sigma_sqd_w     ! Width of individual w plumes (ADG1 parameter)        [-]

    real( kind = core_rknd ), dimension(gr%nz, sclr_dim), intent(out) ::  &
      sclr_1,        & ! Mean of passive scalar (1st PDF component) [units vary]
      sclr_2,        & ! Mean of passive scalar (2nd PDF component) [units vary]
      varnce_sclr_1, & ! Variance of pass. sclr (1st PDF comp.) [(units vary)^2]
      varnce_sclr_2, & ! Variance of pass. sclr (2nd PDF comp.) [(units vary)^2]
      alpha_sclr       ! Factor relating to normalized variance for sclr     [-]

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) ::  &
      w_1_n,         & ! Normalized mean of w (1st PDF component)            [-]
      w_2_n,         & ! Normalized mean of w (2nd PDF component)            [-]
      small_m_w,     & ! Luhar's m (tunable parameter)                       [-]
      big_m_w,       & ! Luhar's M (a function of Luhar's m)                 [-]
      sigma_sqd_w_1, & ! Normalized width parameter of w (1st PDF component) [-]
      sigma_sqd_w_2    ! Normalized width parameter of w (2nd PDF component) [-]

    integer :: i  ! Loop index


    ! Calculate the mixture fraction and the PDF component means and variances
    ! of w.

    ! Reproduce ADG2_w_closure using separate functions
    call calc_Luhar_params( gr%nz, Skw, wp2, &              ! intent(in)
                            wp2, w_tol_sqd, &               ! intent(in)
                            mixt_frac, big_m_w, small_m_w ) ! intent(out)

    call close_Luhar_pdf( gr%nz, wm, wp2, mixt_frac,    & ! intent(in)
                          small_m_w, wp2, w_tol_sqd,    & ! intent(in)
                          sigma_sqd_w_1, sigma_sqd_w_2, & ! intent(out)
                          varnce_w_1, varnce_w_2,       & ! intent(out)
                          w_1_n, w_2_n, w_1, w_2 )        ! intent(out)

    ! Overwrite sigma_sqd_w for consistency with ADG1
    sigma_sqd_w = min( one / ( one + small_m_w**2 ), 0.99_core_rknd )

    ! Calculate the PDF component means and variances of rt.
    call ADG1_ADG2_responder_params( rtm, rtp2, wp2, sqrt_wp2, &       ! In
                                     wprtp, w_1_n, w_2_n, mixt_frac, & ! In
                                     sigma_sqd_w, rt_tol, &            ! In
                                     rt_1, rt_2, varnce_rt_1, &        ! Out
                                     varnce_rt_2, alpha_rt )           ! Out

    ! Calculate the PDF component means and variances of thl.
    call ADG1_ADG2_responder_params( thlm, thlp2, wp2, sqrt_wp2, &      ! In
                                     wpthlp, w_1_n, w_2_n, mixt_frac, & ! In
                                     sigma_sqd_w, thl_tol, &            ! In
                                     thl_1, thl_2, varnce_thl_1, &      ! Out
                                     varnce_thl_2, alpha_thl )          ! Out

    ! Calculate the PDF component means and variances of passive scalars.
    if ( l_scalar_calc ) then
       do i = 1, sclr_dim
          call ADG1_ADG2_responder_params( sclrm(:,i), sclrp2(:,i),     & ! In
                                           wp2, sqrt_wp2, wpsclrp(:,i), & ! In
                                           w_1_n, w_2_n, mixt_frac,     & ! In
                                           sigma_sqd_w, sclr_tol(i),    & ! In
                                           sclr_1(:,i), sclr_2(:,i),    & ! Out
                                           varnce_sclr_1(:,i),          & ! Out
                                           varnce_sclr_2(:,i),          & ! Out
                                           alpha_sclr(:,i) )              ! Out
       enddo ! i=1, sclr_dim
    endif ! l_scalar_calc


    return

  end subroutine ADG2_pdf_driver

  !=============================================================================
  subroutine Luhar_3D_pdf_driver( wm, rtm, thlm, wp2, rtp2, thlp2,      & ! In
                                  Skw, Skrt, Skthl, wprtp, wpthlp,      & ! In
                                  w_1, w_2, rt_1, rt_2, thl_1, thl_2,   & ! Out
                                  varnce_w_1, varnce_w_2, varnce_rt_1,  & ! Out
                                  varnce_rt_2, varnce_thl_1,            & ! Out
                                  varnce_thl_2, mixt_frac )               ! Out

    ! Calculates the PDF parameters using the 3D Luhar PDF closure.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        w_tol_sqd, & ! Constant(s)
        rt_tol,    &
        thl_tol

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      wm,     & ! Mean of w-wind component (vertical velocity)  [m/s] 
      rtm,    & ! Mean of total water mixing ratio              [kg/kg]
      thlm,   & ! Mean of liquid water potential temperature    [K]
      wp2,    & ! Variance of w (overall)                       [m^2/s^2] 
      rtp2,   & ! Variance of r_t (overall)                     [(kg/kg)^2]
      thlp2,  & ! Variance of th_l (overall)                    [K^2]
      Skw,    & ! Skewness of w                                 [-]
      Skrt,   & ! Skewness of rt                                [-]
      Skthl,  & ! Skewness of th_l                              [-]
      wprtp,  & ! Covariance of w and r_t                       [(kg/kg)(m/s)]
      wpthlp    ! Covariance of w and th_l                      [K(m/s)]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) ::  &
      w_1,          & ! Mean of w (1st PDF component)                      [m/s]
      w_2,          & ! Mean of w (2nd PDF component)                      [m/s]
      rt_1,         & ! Mean of r_t (1st PDF component)                  [kg/kg]
      rt_2,         & ! Mean of r_t (2nd PDF component)                  [kg/kg]
      thl_1,        & ! Mean of th_l (1st PDF component)                     [K]
      thl_2,        & ! Mean of th_l (2nd PDF component)                     [K]
      varnce_w_1,   & ! Variance of w (1st PDF component)              [m^2/s^2]
      varnce_w_2,   & ! Variance of w (2nd PDF component)              [m^2/s^2]
      varnce_rt_1,  & ! Variance of r_t (1st PDF component)          [kg^2/kg^2]
      varnce_rt_2,  & ! Variance of r_t (2nd PDF component)          [kg^2/kg^2]
      varnce_thl_1, & ! Variance of th_l (1st PDF component)               [K^2]
      varnce_thl_2, & ! Variance of th_l (2nd PDF component)               [K^2]
      mixt_frac       ! Mixture fraction (weight of 1st PDF component)       [-]

    real( kind = core_rknd ), dimension(gr%nz) ::  &
      w_1_n,           & ! Normalized mean of w (1st PDF component)          [-]
      w_2_n,           & ! Normalized mean of w (2nd PDF component)          [-]
      rt_1_n,          & ! Normalized mean of rt (1st PDF component)         [-]
      rt_2_n,          & ! Normalized mean of rt (2nd PDF component)         [-]
      thl_1_n,         & ! Normalized mean of thl (1st PDF component)        [-]
      thl_2_n,         & ! Normalized mean of thl (2nd PDF component)        [-]
      small_m_w,       & ! Luhar's m (tunable parameter) for w               [-]
      big_m_w,         & ! Luhar's M (a function of Luhar's m) for w         [-]
      small_m_rt,      & ! Luhar's m (tunable parameter) for rt              [-]
      big_m_rt,        & ! Luhar's M (a function of Luhar's m) for rt        [-]
      small_m_thl,     & ! Luhar's m (tunable parameter) for thl             [-]
      big_m_thl,       & ! Luhar's M (a function of Luhar's m) for thl       [-]
      sigma_sqd_w_1,   & ! Normalized width parameter of w (1st PDF comp.)   [-]
      sigma_sqd_w_2,   & ! Normalized width parameter of w (2nd PDF comp.)   [-]
      sigma_sqd_rt_1,  & ! Normalized width parameter of rt (1st PDF comp.)  [-]
      sigma_sqd_rt_2,  & ! Normalized width parameter of rt (2nd PDF comp.)  [-]
      sigma_sqd_thl_1, & ! Normalized width parameter of thl (1st PDF comp.) [-]
      sigma_sqd_thl_2    ! Normalized width parameter of thl (2nd PDF comp.) [-]

    integer :: k    ! Vertical loop index


    do k = 1, gr%nz, 1

       if ( ( abs( Skw(k) ) >= abs( Skthl(k) ) ) &
            .and. ( abs( Skw(k) ) >= abs( Skrt(k) ) ) ) then

          ! w has the greatest magnitude of skewness.

          ! Solve for the w PDF
          call calc_Luhar_params( 1, Skw(k), wp2(k), &                     ! In
                                  wp2(k), w_tol_sqd, &                     ! In
                                  mixt_frac(k), big_m_w(k), small_m_w(k) ) ! Out

          call close_Luhar_pdf( 1, wm(k), wp2(k), mixt_frac(k),     & ! In
                                small_m_w(k), wp2(k), w_tol_sqd,    & ! In
                                sigma_sqd_w_1(k), sigma_sqd_w_2(k), & ! Out
                                varnce_w_1(k), varnce_w_2(k),       & ! Out
                                w_1_n(k), w_2_n(k), w_1(k), w_2(k)  ) ! Out

          ! Solve for the thl PDF
          call backsolve_Luhar_params( Skw(k), Skthl(k),            & ! In
                                       big_m_w(k), mixt_frac(k),    & ! In
                                       big_m_thl(k), small_m_thl(k) ) ! Out

          call close_Luhar_pdf( 1, thlm(k), thlp2(k), mixt_frac(k),     & ! In
                                small_m_thl(k), wpthlp(k), thl_tol**2,  & ! In
                                sigma_sqd_thl_1(k), sigma_sqd_thl_2(k), & ! Out
                                varnce_thl_1(k), varnce_thl_2(k),       & ! Out
                                thl_1_n(k), thl_2_n(k),                 & ! Out
                                thl_1(k), thl_2(k)                      ) ! Out

          ! Solve for the rt PDF
          call backsolve_Luhar_params( Skw(k), Skrt(k),           & ! In
                                       big_m_w(k), mixt_frac(k),  & ! In 
                                       big_m_rt(k), small_m_rt(k) ) ! Out

          call close_Luhar_pdf( 1, rtm(k), rtp2(k), mixt_frac(k),      & ! In
                                small_m_rt(k), wprtp(k), rt_tol**2,    & ! In
                                sigma_sqd_rt_1(k), sigma_sqd_rt_2(k),  & ! Out
                                varnce_rt_1(k), varnce_rt_2(k),        & ! Out
                                rt_1_n(k), rt_2_n(k), rt_1(k), rt_2(k) ) ! Out

       elseif ( ( abs( Skthl(k) ) > abs( Skw(k) ) ) &
                  .and. ( abs( Skthl(k) ) >= abs( Skrt(k) ) ) ) then

          ! theta-l has the greatest magnitude of skewness.

          ! Solve for the thl PDF
          call calc_Luhar_params( 1, Skthl(k), wpthlp(k),     & ! In
                                  thlp2(k), thl_tol**2,       & ! In
                                  mixt_frac(k), big_m_thl(k), & ! Out
                                  small_m_thl(k)              ) ! Out

          ! Solve for the thl PDF
          call close_Luhar_pdf( 1, thlm(k), thlp2(k), mixt_frac(k),     & ! In
                                small_m_thl(k), wpthlp(k), thl_tol**2,  & ! In
                                sigma_sqd_thl_1(k), sigma_sqd_thl_2(k), & ! Out
                                varnce_thl_1(k), varnce_thl_2(k),       & ! Out
                                thl_1_n(k), thl_2_n(k),                 & ! Out
                                thl_1(k), thl_2(k)                      ) ! Out

          ! Solve for the w PDF
          call backsolve_Luhar_params( Skthl(k), Skw(k),           & ! In
                                       big_m_thl(k), mixt_frac(k), & ! In
                                       big_m_w(k), small_m_w(k) )    ! Out

          call close_Luhar_pdf( 1, wm(k), wp2(k), mixt_frac(k),     & ! In
                                small_m_w(k), wp2(k), w_tol_sqd,    & ! In
                                sigma_sqd_w_1(k), sigma_sqd_w_2(k), & ! Out
                                varnce_w_1(k), varnce_w_2(k),       & ! Out
                                w_1_n(k), w_2_n(k), w_1(k), w_2(k)  ) ! Out

          ! Solve for the rt PDF
          call backsolve_Luhar_params( Skthl(k), Skrt(k),           & ! In
                                       big_m_thl(k), mixt_frac(k),  & ! In
                                       big_m_rt(k), small_m_rt(k) )   ! Out

          call close_Luhar_pdf( 1, rtm(k), rtp2(k), mixt_frac(k),      & ! In
                                small_m_rt(k), wprtp(k), rt_tol**2,    & ! In
                                sigma_sqd_rt_1(k), sigma_sqd_rt_2(k),  & ! Out
                                varnce_rt_1(k), varnce_rt_2(k),        & ! Out
                                rt_1_n(k), rt_2_n(k), rt_1(k), rt_2(k) ) ! Out

       else

          ! rt has the greatest magnitude of skewness.

          ! Solve for the rt PDF
          call calc_Luhar_params( 1, Skrt(k), wprtp(k),      & ! In
                                  rtp2(k), rt_tol**2,        & ! In
                                  mixt_frac(k), big_m_rt(k), & ! Out
                                  small_m_rt(k)              ) ! Out

          ! Solve for the rt PDF
          call close_Luhar_pdf( 1, rtm(k), rtp2(k), mixt_frac(k),      & ! In
                                small_m_rt(k), wprtp(k), rt_tol**2,    & ! In
                                sigma_sqd_rt_1(k), sigma_sqd_rt_2(k),  & ! Out
                                varnce_rt_1(k), varnce_rt_2(k),        & ! Out
                                rt_1_n(k), rt_2_n(k), rt_1(k), rt_2(k) ) ! Out

          ! Solve for the w PDF
          call backsolve_Luhar_params( Skrt(k), Skw(k),           & ! In
                                       big_m_rt(k), mixt_frac(k), & ! In
                                       big_m_w(k), small_m_w(k) )   ! Out

          call close_Luhar_pdf( 1, wm(k), wp2(k), mixt_frac(k),     & ! In
                                small_m_w(k), wp2(k), w_tol_sqd,    & ! In
                                sigma_sqd_w_1(k), sigma_sqd_w_2(k), & ! Out
                                varnce_w_1(k), varnce_w_2(k),       & ! Out
                                w_1_n(k), w_2_n(k), w_1(k), w_2(k)  ) ! OUt

          ! Solve for the thl PDF
          call backsolve_Luhar_params( Skrt(k), Skthl(k),           & ! In
                                       big_m_rt(k), mixt_frac(k),   & ! In
                                       big_m_thl(k), small_m_thl(k) ) ! Out

          call close_Luhar_pdf( 1, thlm(k), thlp2(k), mixt_frac(k),     & ! In
                                small_m_thl(k), wpthlp(k), thl_tol**2,  & ! In
                                sigma_sqd_thl_1(k), sigma_sqd_thl_2(k), & ! Out
                                varnce_thl_1(k), varnce_thl_2(k),       & ! Out
                                thl_1_n(k), thl_2_n(k),                 & ! Out
                                thl_1(k), thl_2(k)                      ) ! Out

       endif

    enddo ! k = 1, gr%nz, 1


    return

  end subroutine Luhar_3D_pdf_driver

  !=============================================================================
  subroutine ADG1_w_closure( wm, wp2, Skw, sigma_sqd_w, &        ! In
                             sqrt_wp2, mixt_frac_max_mag, &      ! In
                             w_1, w_2, w_1_n, w_2_n, &           ! Out
                             varnce_w_1, varnce_w_2, mixt_frac ) ! Out

    ! Description:
    ! Calculates the mixture fraction, the PDF component means of w, and the PDF
    ! component variances of w for the Analytic Double Gaussian 1 (ADG1)
    ! closure.  It sets the widths of both w Gaussians to be the same, and
    ! furthermore, relates them to the overall variance of w (<w'^2>) by a
    ! parameter, sigma_sqd_w.  The equation is:
    !
    ! sigma_w_1^2 = sigma_w_2^2 = sigma_sqd_w * <w'^2>;
    !
    ! where sigma_w_1^2 is the variance of w in the 1st PDF component,
    ! sigma_w_2^2 is the variance of w in the 2nd PDF component, and parameter
    ! sigma_sqd_w must have a value within the range 0 <= sigma_sqd_w <= 1.
    !
    ! References:
    ! Golaz, J-C., V. E. Larson, and W. R. Cotton, 2002a: A PDF-based model for
    ! boundary layer clouds. Part I: Method and model description. J. Atmos.
    ! Sci., 59, 3540–3551.
    !
    ! Vincent E. Larson and Jean-Christophe Golaz, 2005: Using Probability
    ! Density Functions to Derive Consistent Closure Relationships among
    ! Higher-Order Moments. Mon. Wea. Rev., 133, 1023–1042.
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        four,      & ! Constant(s)
        one,       &
        one_half,  &
        zero,      &
        w_tol_sqd

    use clubb_precision, only: &
        core_rknd     ! Precision

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      wm,          & ! Mean of w (overall)                   [m/s]
      wp2,         & ! Variance of w (overall)               [m^2/s^2]
      Skw,         & ! Skewness of w                         [-]
      sigma_sqd_w, & ! Widths of each w Gaussian             [-]
      sqrt_wp2       ! Square root of the variance of w      [m/s]

    real( kind = core_rknd ), intent(in) :: &
      mixt_frac_max_mag    ! Maximum allowable mag. of mixt_frac   [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      w_1,        & ! Mean of w (1st PDF component)                [m/s]
      w_2,        & ! Mean of w (2nd PDF component)                [m/s]
      w_1_n,      & ! Normalized mean of w (1st PDF component)     [-]
      w_2_n,      & ! Normalized mean of w (2nd PDF component)     [-]
      varnce_w_1, & ! Variance of w (1st PDF component)            [m^2/s^2]
      varnce_w_2, & ! Variance of w (2nd PDF component)            [m^2/s^2]
      mixt_frac     ! Mixture fraction                             [-]


    !----- Begin Code -----

    where ( wp2 > w_tol_sqd )

       ! Width (standard deviation) parameters are non-zero

       ! The variable "mixt_frac" is the weight of the 1st PDF component.  The
       ! weight of the 2nd PDF component is "1-mixt_frac".  If there isn't any
       ! skewness of w (Sk_w = 0 because w'^3 = 0), mixt_frac = 0.5, and both
       ! PDF components are equally weighted.  If there is positive skewness of
       ! w (Sk_w > 0 because w'^3 > 0), 0 < mixt_frac < 0.5, and the 2nd PDF
       ! component has greater weight than does the 1st PDF component.  If there
       ! is negative skewness of w (Sk_w < 0 because w'^3 < 0),
       ! 0.5 < mixt_frac < 1, and the 1st PDF component has greater weight than
       ! does the 2nd PDF component.
       where ( abs( Skw ) <= 1.0e-5_core_rknd )
          mixt_frac = one_half
       elsewhere
          mixt_frac &
          = one_half &
            * ( one - Skw / sqrt( four * ( one - sigma_sqd_w )**3 + Skw**2 ) )
       endwhere

       ! Clip mixt_frac, and 1 - mixt_frac, to avoid dividing by a small number.
       ! Formula for mixt_frac_max_mag =
       ! 1 - ( 1/2 * ( 1 - Skw_max
       !                   / sqrt( 4*( 1 - sigma_sqd_w )^3 + Skw_max^2 ) ) ),
       ! where sigma_sqd_w is fixed at 0.4 for this calculation.
       mixt_frac = min( max( mixt_frac, one - mixt_frac_max_mag ), &
                        mixt_frac_max_mag )

       ! The normalized mean of w for Gaussian "plume" 1 is w_1_n.  It's value
       ! will always be greater than 0.  As an example, a value of 1.0 would
       ! indicate that the actual mean of w for Gaussian "plume" 1 is found 1.0
       ! standard deviation above the overall mean for w.
       w_1_n &
       = sqrt( ( ( one - mixt_frac ) / mixt_frac ) * ( one - sigma_sqd_w ) )
       ! The normalized mean of w for Gaussian "plume" 2 is w_2_n.  It's value
       ! will always be less than 0.  As an example, a value of -0.5 would
       ! indicate that the actual mean of w for Gaussian "plume" 2 is found 0.5
       ! standard deviations below the overall mean for w.
       w_2_n &
       = -sqrt( ( mixt_frac / ( one - mixt_frac ) ) * ( one - sigma_sqd_w ) )
       ! The mean of w for Gaussian "plume" 1 is w_1.
       w_1 = wm + sqrt_wp2 * w_1_n
       ! The mean of w for Gaussian "plume" 2 is w_2.
       w_2 = wm + sqrt_wp2 * w_2_n

       ! The variance of w for Gaussian "plume" 1 for varnce_w_1.
       varnce_w_1 = sigma_sqd_w * wp2
       ! The variance of w for Gaussian "plume" 2 for varnce_w_2.
       ! The variance in both Gaussian "plumes" is defined to be the same.
       varnce_w_2 = sigma_sqd_w * wp2

    elsewhere

       ! Vertical velocity doesn't vary.
       mixt_frac  = one_half
       w_1_n      = sqrt( one - sigma_sqd_w )
       w_2_n      = -sqrt( one - sigma_sqd_w )
       w_1        = wm
       w_2        = wm
       varnce_w_1 = zero
       varnce_w_2 = zero

    endwhere  ! Widths non-zero


    return

  end subroutine ADG1_w_closure

  !=============================================================================
  subroutine calc_Luhar_params( nz, Skx, wpxp,            & ! In
                                xp2, x_tol_sqd,           & ! In
                                mixt_frac, big_m, small_m ) ! Out

    ! Description:
    ! For the Luhar closure, this subroutine takes Skx (and Skw) as input and
    ! outputs the mixture fraction, big_m, and small_m. This code was written
    ! using the equations and nomenclature of Larson et al. (2002) Appendix
    ! section e.
    !
    ! The relationship between skewness of x (Skx), mixture fraction (a), and
    ! Luhar's small m (m) is given by:
    !
    ! Skx^2 = ( m^2 * ( m^2 + 3 )^2 / ( m^2 + 1 )^3 )
    !         * ( 1 - 2*a )^2 / ( a * ( 1 - a ) ).
    !
    ! Luhar's large M (M) is used to more easily express the factor involving
    ! the m's:
    !
    ! M = ( m^2 + 1 )^3 / ( m^2 * ( m^2 + 3 )^2 ).
    !
    ! The equation involving skewness of x becomes:
    !
    ! Skx^2 = ( 1 / M ) * ( 1 - 2*a )^2 / ( a * ( 1 - a ) );
    !
    ! or:
    !
    ! M * Skx^2 = ( 1 - 2*a )^2 / ( a * ( 1 - a ) ).
    !
    ! This equation can be rewritten as:
    !
    ! ( a * ( 1 - a ) ) * M * Skx^2 = ( 1 - 2*a )^2;
    !
    ! as well as:
    !
    ! ( a - a^2 ) * M * Skx^2 = 1 - 4*a + 4*a^2;
    !
    ! and eventually as:
    !
    ! ( 4 + M * Skx^2 ) * a^2 - ( 4 + M * Skx^2 ) * a + 1 = 0.
    !
    ! Solving the quadratic equation for a:
    !
    ! a = (1/2) * ( 1 +- Skx * sqrt( 1 / ( 4/M + Skx^2 ) ) ).
    !
    ! Since by definition, mu_w_1 >= mu_w_2, a < 0.5 when Skw > 0, the equation
    ! for mixture fraction is:
    !
    ! a = (1/2) * ( 1 - Skx * sqrt( 1 / ( 4/M + Skx^2 ) ) ).
    !
    ! For 3-D Luhar, the variable (w, rt, or theta-l) with the greatest
    ! magnitude of skewness is used to calculate mixture fraction.  Since it is
    ! desirable to still have a < 0.5 when Skw > 0 and a > 0.5 when Skw < 0, the
    ! sign function is used.  The value of Skx is replaced by:
    !
    ! Skx|_adj = sgn( <w'x'> ) * Skx;
    !
    ! where
    !
    ! sgn( <w'x'> ) = | 1 when <w'x'> >= 0
    !                 | -1 when <w'x'> < 0.
    !
    ! Since Skx|_adj^2 = ( sgn( <w'x'> ) * Skx )^2 = ( sgn( <w'x'> ) )^2 * Skx^2
    ! = Skx^2, the equation for mixture fraction is:
    !
    ! a = (1/2) * ( 1 - sgn( <w'x'> ) * Skx * sqrt( 1 / ( 4/M + Skx^2 ) ) ).
    !
    ! When using the ADG2 closure or when using the 3-D Luhar closure when the
    ! variable with the greatest magnitude of skewness is w, Skw = Skx and
    ! sgn( <w'^2> ) is always equal to 1, reducing the equation to its previous
    ! form.

    ! References:
    ! Vincent E. Larson, Jean-Christophe Golaz, and William R. Cotton, 2002:
    ! Small-Scale and Mesoscale Variability in Cloudy Boundary Layers: Joint
    ! Probability Density Functions. J. Atmos. Sci., 59, 3519–3539.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        four,       & ! Constant(s)
        three,      &
        one,        &
        two_thirds, &
        one_half,   &
        one_third,  &
        zero

    use clubb_precision, only: &
        core_rknd     ! Precision

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz    ! Number of vertical levels

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      Skx,  & ! Skewness of x ( <x'^3> / <x'2>^(3/2) )      [-]
      wpxp, & ! Covariance of w and x                       [m/s (x units)]
      xp2     ! Variance of x                               [(x units)^2]

    real( kind = core_rknd ), intent(in) :: &
      x_tol_sqd    ! Tolerance value of x squared           [(x units)^2]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      mixt_frac, & ! Mixture fraction                       [-]
      big_m,     & ! Luhar's M                              [-]
      small_m      ! Luhar's m                              [-]

    ! Local Variables
    real( kind = core_rknd ), dimension(nz) :: &
      small_m_sqd, & ! Luhar's m^2                                           [-]
      sgn_wpxp       ! Sgn(<w'x'>); 1 when <w'x'> >= 0 or -1 when <w'x'> < 0 [-]


    where ( xp2 > x_tol_sqd )

       ! Width (standard deviation) parameters are non-zero

       ! Calculate Luhar's m (small m).
       ! If Skx is very small, then small_m will tend to zero which risks
       ! divide-by-zero.  To ameliorate this problem, we enforce abs( x_1_n )
       ! and abs( x_2_n ) > 0.05.
       ! Note:  Luhar's small_m (m) is the only tunable parameter in the Luhar
       !        closure, so this equation can be changed.  However, the value
       !        of m should go toward 0 as Skx goes toward 0 so that the double
       !        Gaussian reduces to a single Gaussian when the distribution is
       !        unskewed.
       small_m = max( two_thirds * abs( Skx )**one_third, 0.05_core_rknd )

       ! Calculate m^2.
       small_m_sqd = small_m**2

       ! Calculate Luhar's M (big M).
       big_m = ( one + small_m_sqd )**3 &
               / ( ( three + small_m_sqd )**2 * small_m_sqd )

       ! Calculate sgn( <w'x'> ).
       where ( wpxp >= zero )
          sgn_wpxp = one
       elsewhere ! <w'x'> < 0
          sgn_wpxp = -one
       endwhere ! <w'x'> >= 0

       ! Calculate mixture fraction.
       mixt_frac = one_half &
                   * ( one - sgn_wpxp * Skx &
                             * sqrt( one / ( ( four / big_m ) + Skx**2 ) ) )

    elsewhere

       ! Variable x doesn't vary.
       mixt_frac = one_half

       ! For output purposes, set small_m and big_m to 0.
       small_m = zero
       big_m   = zero

    endwhere  ! Widths non-zero


    return

  end subroutine calc_Luhar_params

  !=============================================================================
  subroutine close_Luhar_pdf( nz, xm, xp2, mixt_frac,       & ! In
                              small_m, wpxp, x_tol_sqd,     & ! In
                              sigma_sqd_x_1, sigma_sqd_x_2, & ! Out
                              varnce_x_1, varnce_x_2,       & ! Out
                              x_1_n, x_2_n, x_1, x_2 )        ! Out

    ! Description:
    ! For the Luhar closure, this subroutine takes Skx, xm, xp2, and mixt_frac,
    ! big_m, and small_m (calculated in calc_Luhar_params) as input and outputs
    ! the PDF component means and variances of a variable x in the joint-PDF
    ! according to Luhar et al. (1996).  This code was written using the
    ! equations and nomenclature of Larson et al. (2002) Appendix section e.

    ! References:
    !    Vincent E. Larson, Jean-Christophe Golaz, and William R. Cotton, 2002:
    !    Small-Scale and Mesoscale Variability in Cloudy Boundary Layers: Joint
    !    Probability Density Functions. J. Atmos. Sci., 59, 3519–3539.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,  & ! Constant(s)
        zero

    use clubb_precision, only: &
        core_rknd     ! Precision

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz    ! Number of vertical levels

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      xm,        & ! Mean (overall) of x, <x>                    [(x units)]
      xp2,       & ! Variance (overall) of x, <x'^2>             [(x units)^2]
      mixt_frac, & ! Mixture fraction                            [-]
      small_m,   & ! Luhar's small m                             [-]
      wpxp         ! Covariance of w and x, <w'x'>               [m/s (x units)]

    real( kind = core_rknd ), intent(in) :: &
      x_tol_sqd    ! Tolerance value of x squared                [(x units)^2]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      sigma_sqd_x_1, & ! Normalized width parameter of x (1st PDF component) [-]
      sigma_sqd_x_2, & ! Normalized width parameter of x (1st PDF component) [-]
      varnce_x_1,    & ! Variance of x (1st PDF component)         [(x units)^2]
      varnce_x_2,    & ! Variance of x (2nd PDF component)         [(x units)^2]
      x_1_n,         & ! Normalized mean of x (1st PDF component)            [-]
      x_2_n,         & ! Normalized mean of x (2nd PDF component)            [-]
      x_1,           & ! Mean of x (1st PDF component)               [(x units)]
      x_2              ! Mean of x (2nd PDF component)               [(x units)]

    ! Local Variables
    real( kind = core_rknd), dimension(nz) :: &
      sqrt_xp2, & ! Square root of the variance of x                 [(x units)]
      sgn_wpxp    ! Sgn( <w'x'> ); 1 when <w'x'> >= 0 or -1 when <w'x'> < 0  [-]


    where ( xp2 > x_tol_sqd )

       ! Width (standard deviation) parameters are non-zero

       ! Calculate sgn( <w'x'> ).
       where ( wpxp >= zero )
          sgn_wpxp = one
       elsewhere ! <w'x'> < 0
          sgn_wpxp = -one
       endwhere ! <w'x'> >= 0

       ! Calculate the square root of the overall variance of x.
       sqrt_xp2 = sqrt( xp2 )

       ! Normalized width parameter of x in the 1st PDF component.
       sigma_sqd_x_1 &
       = ( one - mixt_frac ) / ( mixt_frac * ( one + small_m**2 ) )

       ! The variance of x in the 1st PDF component.
       varnce_x_1 = sigma_sqd_x_1 * xp2

       ! Normalized width parameter of x in the 2nd PDF component.
       sigma_sqd_x_2 &
       = mixt_frac / ( ( one - mixt_frac ) * ( one + small_m**2 ) )

       ! The variance of x in the 2nd PDF component.
       varnce_x_2 = sigma_sqd_x_2 * xp2

       ! Normalized mean of x in the 1st PDF component.
       x_1_n = sgn_wpxp * small_m * sqrt( sigma_sqd_x_1 )

       ! Normalized mean of x in the 2nd PDF component.
       x_2_n = -sgn_wpxp * small_m * sqrt( sigma_sqd_x_2 )

       ! The mean of x in the 1st PDF component.
       x_1 = xm + sqrt_xp2 * x_1_n

       ! The mean of x in the 2nd PDF component.
       x_2 = xm + sqrt_xp2 * x_2_n

    elsewhere

       ! Variable x doesn't vary.
       sigma_sqd_x_1 = ( one / ( one + small_m**2 ) )
       sigma_sqd_x_2 = ( one / ( one + small_m**2 ) )
       varnce_x_1    = zero
       varnce_x_2    = zero
       x_1_n         = sgn_wpxp * small_m * sqrt( sigma_sqd_x_1 )
       x_2_n         = -sgn_wpxp * small_m * sqrt( sigma_sqd_x_2 )
       x_1           = xm
       x_2           = xm

    endwhere  ! Widths non-zero


    return

  end subroutine close_Luhar_pdf

  !=============================================================================
  subroutine ADG1_ADG2_responder_params( xm, xp2, wp2, sqrt_wp2, &        ! In
                                         wpxp, w_1_n, w_2_n, mixt_frac, & ! In
                                         sigma_sqd_w, x_tol, &            ! In
                                         x_1, x_2, varnce_x_1, &          ! Out
                                         varnce_x_2, alpha_x )            ! Out

    ! Description:
    ! Calculates the PDF component means and PDF component variances for thl,
    ! rt, or sclr when either the ADG1 PDF or ADG2 PDF are used.
    !
    ! The normalized variance for thl, rt, and sclr for "plume" 1 is:
    !
    ! { 1 - [1/(1-sigma_sqd_w)]*[ <w'x'>^2 / (<w'^2> * <x'^2>) ] / mixt_frac }
    ! * { (1/3)*beta + mixt_frac*( 1 - (2/3)*beta ) };
    !
    ! where "x" stands for thl, rt, or sclr; "mixt_frac" is the weight of
    ! Gaussian "plume" 1, and 0 <= beta <= 3.
    !
    ! The factor { (1/3)*beta + mixt_frac*( 1 - (2/3)*beta ) } does not depend
    ! on which varable "x" stands for.  The factor is multiplied by 2 and
    ! defined as width_factor_1.
    !
    ! The factor
    ! { 1 - [1/(1-sigma_sqd_w)]*[ (w'x')^2 / (w'^2 * x'^2) ] / mixt_frac }
    ! depends on which variable "x" stands for.  It is multiplied by 1/2 and
    ! defined as alpha_x, where "x" stands for thl, rt, or sclr.

    ! Vince Larson added a dimensionless factor so that the
    ! width of plumes in theta_l, rt can vary.
    ! beta is a constant defined in module parameters_tunable
    !   Set 0<beta<3.
    ! beta=1.5_core_rknd recovers Chris Golaz' simplified formula.
    ! 3 Nov 2003

    ! References:

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        two,            & ! Constant(s)
        one,            &
        two_thirds,     &
        one_half,       &
        zero,           &
        zero_threshold, &
        w_tol_sqd

    use parameters_tunable, only: &
        beta    ! Variable(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real ( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      xm,          & ! Mean of x (overall)                    [units vary]
      xp2,         & ! Variance of x (overall)                [(units vary)^2]
      wp2,         & ! Variance of w (overall)                [m^2/s^2]
      sqrt_wp2,    & ! Square root of the variance of w       [m/s]
      wpxp,        & ! Covariance of w and x                  [m/s (units vary)]
      w_1_n,       & ! Normalized mean of w (1st PDF comp.)   [-]
      w_2_n,       & ! Normalized mean of w (2nd PDF comp.)   [-]
      mixt_frac,   & ! Mixture fraction                       [-]
      sigma_sqd_w    ! Width of individual w plumes           [-]

    real ( kind = core_rknd ), intent(in) :: &
      x_tol          ! Tolerance value for x                  [units vary]

    ! Output Variables
    real ( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      x_1,        & ! Mean of x (1st PDF component)             [units vary]
      x_2,        & ! Mean of x (2nd PDF component)             [units vary]
      varnce_x_1, & ! Variance of x (1st PDF component)         [(units vary)^2]
      varnce_x_2, & ! Variance of x (2nd PDF component)         [(units vary)^2]
      alpha_x       ! Factor relating to normalized variance for x           [-]

    ! Local Variables
    ! variables for a generalization of Chris Golaz' closure
    ! varies width of plumes in theta_l, rt
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      width_factor_1, & ! Width factor relating to PDF component 1    [-]
      width_factor_2    ! Width factor relating to PDF component 2    [-]


    where ( wp2 > w_tol_sqd )

       ! Width (standard deviation) parameters of w are non-zero

       width_factor_1 = two_thirds * beta &
                        + two * mixt_frac * ( one - two_thirds * beta )

       width_factor_2 = two - width_factor_1

       where ( xp2 > x_tol**2 )

!          x_1_n = - ( wpxp / ( sqrt( wp2 ) * sqrt( xp2 ) ) ) / w_2_n
!          x_2_n = - ( wpxp / ( sqrt( wp2 ) * sqrt( xp2 ) ) ) / w_1_n

          x_1 = xm - ( wpxp / sqrt_wp2 ) / w_2_n
          x_2 = xm - ( wpxp / sqrt_wp2 ) / w_1_n

          alpha_x = one_half &
                    * ( one - wpxp * wpxp &
                              / ( ( one - sigma_sqd_w ) * wp2 * xp2 ) )

          alpha_x = max( min( alpha_x, one ), zero_threshold )

          ! Vince Larson multiplied original expressions by width_factor_1,2
          !   to generalize scalar skewnesses.  05 Nov 03
          varnce_x_1 = ( alpha_x / mixt_frac * xp2 ) * width_factor_1
          varnce_x_2 = ( alpha_x / ( one - mixt_frac ) * xp2 ) &
                       * width_factor_2

       elsewhere

          ! Variable x doesn't vary.
          x_1        = xm
          x_2        = xm
          varnce_x_1 = zero
          varnce_x_2 = zero
          alpha_x    = one_half

       endwhere ! xp2 > x_tol**2

    elsewhere

       ! Vertical velocity doesn't vary.  Variable x is treated as a single
       ! Gaussian in this situation.
       x_1        = xm
       x_2        = xm
       varnce_x_1 = xp2
       varnce_x_2 = xp2
       alpha_x    = one_half

    endwhere  ! Widths of w non-zero


    return

  end subroutine ADG1_ADG2_responder_params

  !=============================================================================
  subroutine backsolve_Luhar_params( Sk_max, Skx,          & ! In
                                     big_m_max, mixt_frac, & ! In
                                     big_m_x, small_m_x    ) ! Out

    ! Description:
    ! This subroutine calculates Luhar's big_m and small_m for the variate 'x'
    ! consistent with the mixture fraction of the variate with the largest
    ! skewness.
    !
    ! The relationship between skewness of x (Skx), mixture fraction (a), and
    ! Luhar's small m (m) is given by:
    !
    ! Skx^2 = ( m^2 * ( m^2 + 3 )^2 / ( m^2 + 1 )^3 )
    !         * ( 1 - 2*a )^2 / ( a * ( 1 - a ) ).
    !
    ! Moving the factor involving mixture fraction to the right-hand side:
    !
    ! ( ( a * ( 1 - a ) ) / ( 1 - 2*a )^2 ) * Skx^2
    ! = m^2 * ( m^2 + 3 )^2 / ( m^2 + 1 )^3.
    !
    ! This can be rewritten as:
    !
    ! ( ( a * ( 1 - a ) ) / ( 1 - 2*a )^2 ) * Skx^2
    ! = ( m^6 + 6*m^4 + 9*m^2 ) / ( m^6 + 3*m^4 + 3*m^2 + 1 ).
    !
    ! Setting alpha = ( ( a * ( 1 - a ) ) / ( 1 - 2*a )^2 ) * Skx^2, the
    ! equation can be rewritten as:
    !
    ! ( m^6 + 3*m^4 + 3*m^2 + 1 ) * alpha = m^6 + 6*m^4 + 9*m^2.
    !
    ! This can be rearranged and rewritten as:
    !
    ! ( alpha - 1 ) * m^6 + ( 3 * alpha - 6 ) * m^4
    ! + ( 3 * alpha - 9 ) * m^2 + alpha = 0.
    !
    ! This can be rewritten again as:
    !
    ! ( alpha - 1 ) * (m^2)^3 + ( 3 * alpha - 6 ) * (m^2)^2
    ! + ( 3 * alpha - 9 ) * (m^2) + alpha = 0.
    !
    ! The goal is to solve for m^2, and then take the square root of m^2 to
    ! solve for m.  This can be accomplished by using the cubic formula (with
    ! the l_use_cubic_backsolve option), or else by a quadratic approximation.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        nine,     & ! Constant(s)
        six,      &
        five,     &
        four,     &
        three,    &
        two,      &
        one,      &
        one_half, &
        zero,     &
        eps

    use clubb_precision, only: &
        core_rknd     ! Precision

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Sk_max,    & ! Maximum skewness                                      [-]
      Skx,       & ! Skewness of the variate solving small_m and big_m for [-]
      big_m_max, & ! Luhar's big_m of the variate with maximum skewness    [-]
      mixt_frac    ! Mixture fraction                                      [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      big_m_x,   & ! Luhar's big_m for the variate being solved for    [-]
      small_m_x    ! Luhar's small_m for the variate being solved for  [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha,     & ! 1 / big_m_x
      a,         & ! For readability, quadratic equation
      b,         &
      c,         &
      alpha_upr, &
      alpha_low, &
      discrim

    real( kind = core_rknd ), dimension(1) :: &
      a_coef_in, & ! Coefficient a (of x^3) in a*x^3 + b*x^2 + c^x + d = 0   [-]
      b_coef_in, & ! Coefficient b (of x^2) in a*x^3 + b*x^2 + c^x + d = 0   [-]
      c_coef_in, & ! Coefficient c (of x) in a*x^3 + b*x^2 + c^x + d = 0     [-]
      d_coef_in    ! Coefficient d in a*x^3 + b*x^2 + c^x + d = 0            [-]

    ! Flag to backsolve for m^2 using cubic formula
    logical, parameter :: &
      l_use_cubic_backsolve = .true.


    if ( l_use_cubic_backsolve ) then

       if ( abs( mixt_frac - one_half ) < 0.001_core_rknd ) then

          ! When mixture fraction = 0.5 (based on the variable with the largest
          ! magnitude of skewness), all variables must have a skewness of 0.
          ! Set m to the minimum threshold of 0.05.
          small_m_x = 0.05_core_rknd

          ! Calculate the corresponding value of big_m_x.
          big_m_x = ( one + small_m_x**2 )**3 &
                    / ( ( three + small_m_x**2 )**2 * small_m_x**2 )

       elseif ( abs(Skx) < eps ) then

          ! Mixture fraction /= 0.5 because the variable with the largest
          ! magnitude of skewness has a skewness /= 0.  However, variable x has
          ! a skewness of 0.  In order to reproduce the correct skewness for
          ! variable x, set m to 0 (regardless of minimum thresholds used in
          ! other parts of the code).
          small_m_x = zero

          ! The value of big_m_x should be inf.  Set it to huge.  This is not
          ! used in any calculation, anyway.
          big_m_x = huge( big_m_x )

       else  ! mixt_frac /= 0.5 and Skx /= 0

          ! Backsolve for m, given mixt_frac and Skx.

          ! alpha = 1/M is given by:
          ! [ mixt_frac * ( 1 - mixt_frac ) / ( 1 - 2 * mixt_frac )^2 ] * Skx^2.
          alpha = ( mixt_frac * ( one - mixt_frac ) &
                    / ( one - two * mixt_frac )**2 ) * Skx**2

          ! Calculate big_m_x.
          big_m_x = one / alpha

          ! Solve the cubic equation for m^2:
          ! ( alpha - 1 ) * (m^2)^3 + ( 3 * alpha - 6 ) * (m^2)^2
          ! + ( 3 * alpha - 9 ) * (m^2) + alpha = 0.
          ! The largest root is preferred.
          a_coef_in(1) = alpha - one
          b_coef_in(1) = three * alpha - six
          c_coef_in(1) = three * alpha - nine
          d_coef_in(1) = alpha
          small_m_x &
          = sqrt( max( max_cubic_root( a_coef_in, b_coef_in, &
                                       c_coef_in, d_coef_in ), &
                       0.05_core_rknd**2 ) )

       endif

    else ! original formualation

       alpha = ( Skx**2 / (max(Sk_max**2, eps) * big_m_max) )  ! 1 / big_m_x

       ! This limit keeps the discriminant >= 0
       alpha_upr = two * sqrt( 13.0_core_rknd ) - five

       alpha_low = eps

       ! For this approximation, alpha must be less than 2*sqrt(13) - 5 to get a
       ! real ans.
       alpha = min( alpha, alpha_upr )

       ! For testing, eliminate possibility of divide by zero
       alpha = max( alpha, alpha_low )

       ! Use a piece-wise approximation
       if ( alpha < one ) then

          a = max( three * alpha - six, eps) ! Prevent divide by zero
          b = three * alpha - nine
          c = alpha

          discrim = b**2 - four * a * c
          small_m_x = sqrt( ( - b - sqrt( discrim ) ) / ( two * a ) )

       else ! alpha >= 1

          ! For this approximation, alpha must be less than 2*sqrt(13) - 5
          ! to get a real ans.
          alpha = min( alpha, two )

          a = max( six * alpha - nine, eps) ! Prevent divide by zero
          b = -six
          c = two * alpha - one

          discrim = b**2 - four * a * c
          small_m_x = sqrt( ( - b - sqrt( discrim ) ) / ( two * a ) )

       endif ! alpha < 1

       ! Clip consistently with subroutine calc_Luhar_params
       small_m_x = max( 5e-2_core_rknd, small_m_x)

       big_m_x = one / alpha

    endif ! l_use_cubic_backsolve


    return

  end subroutine backsolve_Luhar_params

  !=============================================================================
  pure function max_cubic_root( a_coef, b_coef, c_coef, d_coef ) &
  result( max_root )

    ! Description:
    ! Calculates the largest root that results from solving a cubic equation of
    ! the form a*x^3 + b*x^2 + c*x + d = 0.
    !
    ! This is done to backsolve for m^2 for the 3-D Luhar closure, given the
    ! values of mixt_frac and Skx.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        eps                ! Constant(s)
        
    use calc_roots, only: &
        cubic_solve,     & ! Procedure(s)
        quadratic_solve

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(1), intent(in) :: &
      a_coef, & ! Coefficient a (of x^3) in a*x^3 + b*x^2 + c^x + d = 0    [-]
      b_coef, & ! Coefficient b (of x^2) in a*x^3 + b*x^2 + c^x + d = 0    [-]
      c_coef, & ! Coefficient c (of x) in a*x^3 + b*x^2 + c^x + d = 0      [-]
      d_coef    ! Coefficient d in a*x^3 + b*x^2 + c^x + d = 0             [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      max_root    ! Maximum root that solves the cubic equation            [-]

    ! Local Variables
    complex( kind = core_rknd ), dimension(1,3) :: &
      cubic_roots    ! Roots of x that satisfy a*x^3 + b*x^2 + c*x + d = 0 [-]

    complex( kind = core_rknd ), dimension(1,2) :: &
      quadratic_roots    ! Roots of x that satisfy b*x^2 + c*x + d = 0     [-]

    real( kind = core_rknd ) :: &
      a_coef_thresh, & ! Minimum threshold of |a| to use cubic solver      [-]
      b_coef_thresh    ! Minimum threshold of |b| to use quadratic solver  [-]


    ! Calculate a minimum threshold for |a| to call this a cubic equation.
    a_coef_thresh = 0.001_core_rknd &
                    * max( abs(b_coef(1)), abs(c_coef(1)), abs(d_coef(1)) )

    ! Calculate a minimum threshold for |b| to call this a quadratic equation.
    ! This only matters when |a| <= a_coef_thresh.
    b_coef_thresh = 0.001_core_rknd * max( abs(c_coef(1)), abs(d_coef(1)) )

    if ( abs( a_coef(1) ) > a_coef_thresh ) then

       ! The equation is a cubic equation.
       cubic_roots = cubic_solve( 1, a_coef, b_coef, c_coef, d_coef )

       if ( abs(aimag( cubic_roots(1,2) )) < eps .and.  &
            abs(aimag( cubic_roots(1,3) )) < eps ) then

          ! Find the maximum root of the three roots.
          max_root = max( real( cubic_roots(1,1), kind = core_rknd ), &
                          real( cubic_roots(1,2), kind = core_rknd ), &
                          real( cubic_roots(1,3), kind = core_rknd ) )

       else  ! cubic_roots(2) and cubic_roots(3) are complex.

          max_root = real( cubic_roots(1,1), kind = core_rknd )

       endif

    elseif ( abs( b_coef(1) ) > b_coef_thresh ) then

       ! The equation is a quadratic equation, since a = 0, but b /= 0.
       ! This should very rarely occur for 3-D Luhar.  When it does, the result
       ! will always be two real-valued roots.
       quadratic_roots = quadratic_solve( 1, b_coef, c_coef, d_coef )

       ! Find the maximum root of the two roots.
       max_root = max( real( quadratic_roots(1,1), kind = core_rknd ), &
                       real( quadratic_roots(1,2), kind = core_rknd ) )

    else ! |a| = 0 and |b| = 0

       ! The equation is a linear equation.
       ! This won't happen for 3-D Luhar.
       max_root = - d_coef(1) / c_coef(1)

    endif ! |a| > 0


    return

  end function max_cubic_root

  !=============================================================================

end module adg1_adg2_3d_luhar_pdf
