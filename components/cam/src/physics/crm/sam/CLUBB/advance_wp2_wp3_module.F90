!------------------------------------------------------------------------
! $Id: advance_wp2_wp3_module.F90 6146 2013-04-05 18:02:22Z raut@uwm.edu $
!===============================================================================
module advance_wp2_wp3_module

  implicit none

  private ! Default Scope

  public :: advance_wp2_wp3

  private :: wp23_solve, & 
             wp23_lhs, & 
             wp23_rhs, & 
             wp2_term_ta_lhs, & 
             wp2_terms_ac_pr2_lhs, & 
             wp2_term_dp1_lhs, & 
             wp2_term_pr1_lhs, & 
             wp2_terms_bp_pr2_rhs, & 
             wp2_term_dp1_rhs, &
             wp2_term_pr3_rhs, & 
             wp2_term_pr1_rhs, & 
             wp3_terms_ta_tp_lhs, & 
             wp3_terms_ac_pr2_lhs, & 
             wp3_term_pr1_lhs, & 
             wp3_terms_bp1_pr2_rhs, & 
             wp3_term_pr1_rhs, &
             wp3_term_bp2_rhs

! private :: wp3_terms_ta_tp_rhs

  ! Private named constants to avoid string comparisons
  integer, parameter, private :: &
    clip_wp2 = 12 ! Named constant for wp2 clipping.
                  ! NOTE: This must be the same as the clip_wp2 declared in
                  ! clip_explicit!

  contains

  !=============================================================================
  subroutine advance_wp2_wp3( dt, sfc_elevation, sigma_sqd_w, wm_zm, wm_zt, &
                              a3, a3_zt, wp3_on_wp2, &
                              wpthvp, wp2thvp, um, vm, upwp, vpwp, &
                              up2, vp2, Kh_zm, Kh_zt, tau_zm, tau_zt, &
                              Skw_zm, Skw_zt, rho_ds_zm, rho_ds_zt, &
                              invrs_rho_ds_zm, invrs_rho_ds_zt, radf, &
                              thv_ds_zm, thv_ds_zt, mixt_frac, &
                              wp2, wp3, wp3_zm, wp2_zt, err_code )

    ! Description:
    ! Advance w'^2 and w'^3 one timestep.

    ! References:
    ! Eqn. 12 & 18 on p. 3545--3546 of
    ! ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    !   Method and Model Description'' Golaz, et al. (2002)
    ! JAS, Vol. 59, pp. 3540--3551.

    ! See also
    ! ``Equations for CLUBB'', Section 6:
    ! /Implict solution for the vertical velocity moments/
    !------------------------------------------------------------------------

    use grid_class, only:  & 
        gr,     & ! Variable(s)
        zt2zm,  & ! Procedure(s)
        zm2zt

    use parameters_tunable, only:  & 
        C11c,  & ! Variable(s)
        C11b,  & 
        C11,  & 
        C1c,  & 
        C1b,  & 
        C1,  & 
        c_K1,  & 
        c_K8

    use stats_type, only: & 
        stat_update_var

    use stats_variables, only: &
        iC1_Skw_fnc, &
        iC11_Skw_fnc, &
        zm, &
        zt, &
        l_stats_samp

    use constants_clubb, only:  & 
        fstderr    ! Variable(s)

    use model_flags, only: &
        l_hyper_dfsn ! Variable(s)

    use clubb_precision, only:  & 
        time_precision, & ! Variable(s)
        core_rknd

    use error_code, only:  & 
        fatal_error,  & ! Procedure(s)
        clubb_at_least_debug_level

    use error_code, only: &
      clubb_var_out_of_range ! Constant(s)

    implicit none

    intrinsic :: exp

    ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
      dt                 ! Model timestep                            [s]

    real( kind = core_rknd ), intent(in) ::  &
      sfc_elevation      ! Elevation of ground level                 [m AMSL]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  & 
      sigma_sqd_w,     & ! sigma_sqd_w (momentum levels)             [-]
      wm_zm,           & ! w wind component on momentum levels       [m/s]
      wm_zt,           & ! w wind component on thermodynamic levels  [m/s]
      a3,              & ! a_3 (momentum levels); See eqn. 25 in `Equations for CLUBB' [-]
      a3_zt,           & ! a_3 interpolated to thermodynamic levels  [-]
      wp3_on_wp2,      & ! Smoothed version of wp3 / wp2             [m/s]
      wpthvp,          & ! w'th_v' (momentum levels)                 [K m/s]
      wp2thvp,         & ! w'^2th_v' (thermodynamic levels)          [K m^2/s^2]
      um,              & ! u wind component (thermodynamic levels)   [m/s]
      vm,              & ! v wind component (thermodynamic levels)   [m/s]
      upwp,            & ! u'w' (momentum levels)                    [m^2/s^2]
      vpwp,            & ! v'w' (momentum levels)                    [m^2/s^2]
      up2,             & ! u'^2 (momentum levels)                    [m^2/s^2]
      vp2,             & ! v'^2 (momentum levels)                    [m^2/s^2]
      Kh_zm,           & ! Eddy diffusivity on momentum levels       [m^2/s]
      Kh_zt,           & ! Eddy diffusivity on thermodynamic levels  [m^2/s]
      tau_zm,          & ! Time-scale tau on momentum levels         [s]
      tau_zt,          & ! Time-scale tau on thermodynamic levels    [s]
      Skw_zm,          & ! Skewness of w on momentum levels          [-]
      Skw_zt,          & ! Skewness of w on thermodynamic levels     [-]
      rho_ds_zm,       & ! Dry, static density on momentum levels    [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ momentum levs. [m^3/kg]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. levs.  [m^3/kg]
      radf,            & ! Buoyancy production at the CL top         [m^2/s^3]
      thv_ds_zm,       & ! Dry, base-state theta_v on momentum levs. [K]
      thv_ds_zt,       & ! Dry, base-state theta_v on thermo. levs.  [K]
      mixt_frac          ! Weight of 1st normal distribution         [-]

    ! Input/Output
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) ::  & 
      wp2,  & ! w'^2 (momentum levels)                    [m^2/s^2]
      wp3,  & ! w'^3 (thermodynamic levels)               [m^3/s^3]
      wp3_zm  ! w'^3 interpolated to momentum levels      [m^3/s^3]

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) ::  &
      wp2_zt  ! w'^2 interpolated to thermodyamic levels  [m^2/s^2]

    integer, intent(inout) :: err_code ! Diagnostic

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) ::  & 
      tauw3t  ! Currently just tau_zt                           [s]

    ! Eddy Diffusion for w'^2 and w'^3.
    real( kind = core_rknd ), dimension(gr%nz) :: Kw1    ! w'^2 coef. eddy diff.  [m^2/s]
    real( kind = core_rknd ), dimension(gr%nz) :: Kw8    ! w'^3 coef. eddy diff.  [m^2/s]

    ! Internal variables for C11 function, Vince Larson 13 Mar 2005
    ! Brian added C1 function.
    real( kind = core_rknd ), dimension(gr%nz) ::  & 
      C1_Skw_fnc,  & ! C_1 parameter with Sk_w applied              [-]
      C11_Skw_fnc    ! C_11 parameter with Sk_w applied             [-]
    ! End Vince Larson's addition.

    integer :: &
      nsub,   & ! Number of subdiagonals in the LHS matrix.
      nsup      ! Number of superdiagonals in the LHS matrix.

    integer :: k ! Array indices

    integer :: wp2_wp3_err_code ! Error code from solving for wp2/wp3


    !-----------------------------------------------------------------------



!       Define tauw

!        tauw3t = tau_zt
!     .           / ( 1.
!     .                   + 3.0_core_rknd * max(
!     .                     min(1.-(mixt_frac-0.01_core_rknd)/(0.05_core_rknd-0.01_core_rknd)
!     .                         ,1.)
!     .                         ,0.)
!     .                   + 3.0_core_rknd * max(
!     .                     min(1.-(mixt_frac-0.99_core_rknd)/(0.95_core_rknd-0.99_core_rknd)
!     .                         ,1.)
!     .                         ,0.)
!     .              )

!        do k=1,gr%nz
!
!          Skw = abs( wp3(k)/max(wp2(k),1.e-8)**1.5_core_rknd )
!          Skw = min( 5.0_core_rknd, Skw )
!          tauw3t(k) = tau_zt(k) / ( 0.005_core_rknd*Skw**4 + 1.0_core_rknd )
!
!        end do

    tauw3t = tau_zt

    ! Vince Larson added code to make C11 function of Skw. 13 Mar 2005
    ! If this code is used, C11 is no longer relevant, i.e. constants
    !    are hardwired.

    ! Calculate C_{1} and C_{11} as functions of skewness of w.
    ! The if..then here is only for computational efficiency -dschanen 2 Sept 08
    if ( C11 /= C11b ) then
      C11_Skw_fnc(1:gr%nz) =  & 
        C11b + (C11-C11b)*EXP( -(1.0_core_rknd/2.0_core_rknd) * (Skw_zt(1:gr%nz)/C11c)**2 )
    else
      C11_Skw_fnc(1:gr%nz) = C11b
    end if

    ! The if..then here is only for computational efficiency -dschanen 2 Sept 08
    if ( C1 /= C1b ) then
      C1_Skw_fnc(1:gr%nz) =  & 
        C1b + (C1-C1b)*EXP( -(1.0_core_rknd/2.0_core_rknd) * (Skw_zm(1:gr%nz)/C1c)**2 )
    else
      C1_Skw_fnc(1:gr%nz) = C1b 
    end if

    !C11_Skw_fnc = C11
    !C1_Skw_fnc = C1

    if ( clubb_at_least_debug_level( 2 ) ) then
      ! Assertion check for C11_Skw_fnc
      if ( any( C11_Skw_fnc(:) > 1._core_rknd ) .or. any( C11_Skw_fnc(:) < 0._core_rknd ) ) then
        write(fstderr,*) "The C11_Skw_fnc is outside the valid range for this variable"
        err_code = clubb_var_out_of_range
        return
      end if
    end if

    if ( l_stats_samp ) then
      call stat_update_var( iC11_Skw_fnc, C11_Skw_fnc, zt )
      call stat_update_var( iC1_Skw_fnc, C1_Skw_fnc, zm )
    endif

    ! Define the Coefficent of Eddy Diffusivity for the wp2 and wp3.
    do k = 1, gr%nz, 1

      ! Kw1 is used for wp2, which is located on momentum levels.
      ! Kw1 is located on thermodynamic levels.
      ! Kw1 = c_K1 * Kh_zt
      Kw1(k) = c_K1 * Kh_zt(k)

      ! Kw8 is used for wp3, which is located on thermodynamic levels.
      ! Kw8 is located on momentum levels.
      ! Note: Kw8 is usually defined to be 1/2 of Kh_zm.
      ! Kw8 = c_K8 * Kh_zm
      Kw8(k) = c_K8 * Kh_zm(k)

    enddo

    ! Declare the number of subdiagonals and superdiagonals in the LHS matrix.
    if ( l_hyper_dfsn ) then
       ! There are nine overall diagonals (including four subdiagonals
       ! and four superdiagonals).
       nsub = 4
       nsup = 4
    else
       ! There are five overall diagonals (including two subdiagonals
       ! and two superdiagonals).
       nsub = 2
       nsup = 2
    endif

    ! Solve semi-implicitly
    call wp23_solve( dt, sfc_elevation, sigma_sqd_w, wm_zm, wm_zt, & ! Intent(in)
                     a3, a3_zt, wp3_on_wp2, &  ! Intent(in)
                     wpthvp, wp2thvp, um, vm, upwp, vpwp,    & ! Intent(in)
                     up2, vp2, Kw1, Kw8, Kh_zt, Skw_zt, tau_zm, tauw3t,    & ! Intent(in)
                     C1_Skw_fnc, C11_Skw_fnc, rho_ds_zm, rho_ds_zt, & ! Intent(in)
                     invrs_rho_ds_zm, invrs_rho_ds_zt, radf, thv_ds_zm,   & ! Intent(in)
                     thv_ds_zt, nsub, nsup,                         & ! Intent(in)
                     wp2, wp3, wp3_zm, wp2_zt, wp2_wp3_err_code     ) ! Intent(inout)

!       Error output
!       Joshua Fasching Feb 2008
    if ( fatal_error( wp2_wp3_err_code ) ) then  
     
      if ( clubb_at_least_debug_level( 1 ) ) then
        write(fstderr,*) "Errors in advance_wp2_wp3"

        write(fstderr,*) "Intent(in)"

        write(fstderr,*) "dt = ", dt
        write(fstderr,*) "sfc_elevation = ", sfc_elevation
        write(fstderr,*) "sigma_sqd_w = ", sigma_sqd_w
        write(fstderr,*) "wm_zm = ", wm_zm
        write(fstderr,*) "wm_zt = ", wm_zt
        write(fstderr,*) "wpthvp = ", wpthvp
        write(fstderr,*) "wp2thvp = ", wp2thvp
        write(fstderr,*) "um = ", um
        write(fstderr,*) "vm = ", vm
        write(fstderr,*) "upwp = ", upwp
        write(fstderr,*) "vpwp = ", vpwp
        write(fstderr,*) "up2 = ", up2
        write(fstderr,*) "vp2 = ", vp2
        write(fstderr,*) "Kh_zm = ", Kh_zm
        write(fstderr,*) "Kh_zt = ", Kh_zt
        write(fstderr,*) "tau_zm = ", tau_zm
        write(fstderr,*) "tau_zt = ", tau_zt
        write(fstderr,*) "Skw_zm = ", Skw_zm
        write(fstderr,*) "Skw_zt = ", Skw_zt
        write(fstderr,*) "mixt_frac = ", mixt_frac
        write(fstderr,*) "wp2zt = ", wp2_zt

        write(fstderr,*) "Intent(in/out)"

        write(fstderr,*) "wp2 = ", wp2
        write(fstderr,*) "wp3 = ", wp3

      end if

      err_code = wp2_wp3_err_code
    end if ! fatal error

    return

  end subroutine advance_wp2_wp3

  !=============================================================================
  subroutine wp23_solve( dt, sfc_elevation, sigma_sqd_w, wm_zm, wm_zt, &
                         a3, a3_zt, wp3_on_wp2, &
                         wpthvp, wp2thvp, um, vm, upwp, vpwp, &
                         up2, vp2, Kw1, Kw8, Kh_zt, Skw_zt, tau1m, tauw3t, &
                         C1_Skw_fnc, C11_Skw_fnc, rho_ds_zm, rho_ds_zt, &
                         invrs_rho_ds_zm, invrs_rho_ds_zt, radf, thv_ds_zm, &
                         thv_ds_zt, nsub, nsup, &
                         wp2, wp3, wp3_zm, wp2_zt, err_code )

    ! Description:
    ! Decompose, and back substitute the matrix for wp2/wp3

    ! References:
    ! _Equations for CLUBB_ section 6.3
    !------------------------------------------------------------------------

    use grid_class, only:  & 
        gr  ! Variable(s) 

    use grid_class, only:  & 
        zm2zt, & ! Function(s)
        zt2zm, & 
        ddzt

    use constants_clubb, only: & 
        w_tol_sqd,      & ! Variables(s)
        eps,           &
        zero_threshold, &
        fstderr

    use model_flags, only:  & 
        l_tke_aniso,  & ! Variable(s)
        l_hyper_dfsn, &
        l_hole_fill,  &
        l_gmres

    use clubb_precision, only:  & 
        time_precision, &  ! Variable(s)
        core_rknd

    use lapack_wrap, only:  & 
        band_solve,  & ! Procedure(s) 
        band_solvex

    use fill_holes, only: & 
        fill_holes_driver

    use clip_explicit, only: &
        clip_variance, & ! Procedure(s)
        clip_skewness

    use stats_type, only: & 
        stat_begin_update,  & ! Procedure(s)
        stat_update_var_pt, &
        stat_end_update,  &
        stat_end_update_pt

    use stats_variables, only:  & 
        zm,         & ! Variable(s)
        zt, & 
        sfc, & 
        l_stats_samp, & 
        iwp2_ta, & 
        iwp2_ma, & 
        iwp2_pd, & 
        iwp2_ac, & 
        iwp2_dp1, & 
        iwp2_dp2, & 
        iwp2_pr1, & 
        iwp2_pr2, &
        iwp2_4hd, &
        iwp3_ta, & 
        iwp3_ma, & 
        iwp3_tp, & 
        iwp3_ac, & 
        iwp3_dp1, & 
        iwp3_pr1, &
        iwp3_pr2, &
        iwp3_4hd, &
        iwp23_matrix_condt_num

    use stats_variables, only:  & 
        zmscr01, &
        zmscr02, &
        zmscr03, &
        zmscr04, &
        zmscr05, &
        zmscr06, &
        zmscr07, &
        zmscr08, &
        zmscr09, &
        zmscr10, &
        zmscr11, &
        zmscr12, &
        zmscr13, &
        zmscr14, &
        zmscr15, &
        zmscr16, &
        zmscr17, &
        ztscr01, &
        ztscr02

    use stats_variables, only: &
        ztscr03, &
        ztscr04, &
        ztscr05, &
        ztscr06, &
        ztscr07, &
        ztscr08, &
        ztscr09, &
        ztscr10, &
        ztscr11, &
        ztscr12, &
        ztscr13, &
        ztscr14, &
        ztscr15, &
        ztscr16, &
        ztscr17, &
        ztscr18, &
        ztscr19, &
        ztscr20, &
        ztscr21

    implicit none

    ! External
    intrinsic :: max, min, sqrt

    ! Parameter Constants
    integer, parameter :: & 
      nrhs = 1      ! Number of RHS vectors

    ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
      dt                 ! Timestep                                  [s]

    real( kind = core_rknd ), intent(in) ::  &
      sfc_elevation      ! Elevation of ground level                 [m AMSL]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  & 
      sigma_sqd_w,     & ! sigma_sqd_w (momentum levels)             [-]
      wm_zm,           & ! w wind component on momentum levels       [m/s]
      wm_zt,           & ! w wind component on thermodynamic levels  [m/s]
      a3,              & ! a_3 (momentum levels); See eqn. 25 in `Equations for CLUBB' [-]
      a3_zt,           & ! a_3 interpolated to thermodynamic levels  [-]
      wp3_on_wp2,      & ! Smoothed version of wp3 / wp2             [m/s]
      wpthvp,          & ! w'th_v' (momentum levels)                 [K m/s]
      wp2thvp,         & ! w'^2th_v' (thermodynamic levels)          [K m^2/s^2]
      um,              & ! u wind component (thermodynamic levels)   [m/s]
      vm,              & ! v wind component (thermodynamic levels)   [m/s]
      upwp,            & ! u'w' (momentum levels)                    [m^2/s^2]
      vpwp,            & ! v'w' (momentum levels)                    [m^2/s^2]
      up2,             & ! u'^2 (momentum levels)                    [m^2/s^2]
      vp2,             & ! v'^2 (momentum levels)                    [m^2/s^2]
      Kw1,             & ! Coefficient of eddy diffusivity for w'^2  [m^2/s]
      Kw8,             & ! Coefficient of eddy diffusivity for w'^3  [m^2/s]
      Kh_zt,           & ! Eddy diffusivity on thermodynamic levels  [m^2/s]
      Skw_zt,          & ! Skewness of w on thermodynamic levels     [-]
      tau1m,           & ! Time-scale tau on momentum levels         [s]
      tauw3t,          & ! Time-scale tau on thermodynamic levels    [s]
      C1_Skw_fnc,      & ! C_1 parameter with Sk_w applied           [-]
      C11_Skw_fnc,     & ! C_11 parameter with Sk_w applied          [-]
      rho_ds_zm,       & ! Dry, static density on momentum levels    [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ momentum levs. [m^3/kg]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. levs.  [m^3/kg]
      radf,            & ! Buoyancy production at CL top             [m^2/s^3]
      thv_ds_zm,       & ! Dry, base-state theta_v on momentum levs. [K]
      thv_ds_zt          ! Dry, base-state theta_v on thermo. levs.  [K]

    integer, intent(in) :: &
      nsub,   & ! Number of subdiagonals in the LHS matrix.
      nsup      ! Number of superdiagonals in the LHS matrix.

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) ::  & 
      wp2,  & ! w'^2 (momentum levels)                            [m^2/s^2]
      wp3,  & ! w'^3 (thermodynamic levels)                       [m^3/s^3]
      wp3_zm  ! w'^3 interpolated to momentum levels      [m^3/s^3]

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) ::  &
      wp2_zt  ! w'^2 interpolated to thermodyamic levels          [m^2/s^2]

    integer, intent(inout) :: err_code ! Have any errors occured?

    ! Local Variables
    real( kind = core_rknd ), dimension(nsup+nsub+1,2*gr%nz) ::  & 
      lhs ! Implicit contributions to wp2/wp3 (band diag. matrix)

    real( kind = core_rknd ), dimension(2*gr%nz) ::  & 
      rhs   ! RHS of band matrix

!        real, target, dimension(2*gr%nz) ::
    real( kind = core_rknd ), dimension(2*gr%nz) ::  & 
      solut ! Solution to band diagonal system.

    real( kind = core_rknd ), dimension(gr%nz) ::  & 
      a1,   & ! a_1 (momentum levels); See eqn. 23 in `Equations for CLUBB' [-]
      a1_zt   ! a_1 interpolated to thermodynamic levels                    [-]

!      real, dimension(gr%nz) ::  &
!        wp2_n ! w'^2 at the previous timestep           [m^2/s^2]

    real( kind = core_rknd ) ::  & 
      rcond  ! Est. of the reciprocal of the condition #

    ! Array indices
    integer :: k, km1, km2, kp1, kp2, k_wp2, k_wp3

    ! Set logical to true for Crank-Nicholson diffusion scheme
    ! or to false for completely implicit diffusion scheme.
    ! Note:  Although Crank-Nicholson diffusion has usually been used for wp2
    !        and wp3 in the past, we found that using completely implicit
    !        diffusion stabilized the deep convective cases more while having
    !        almost no effect on the boundary layer cases.  Brian; 1/4/2008.
!    logical, parameter :: l_crank_nich_diff = .true.
    logical, parameter :: l_crank_nich_diff = .false.

    ! Define a_1 and a_3 (both are located on momentum levels).
    ! They are variables that are both functions of sigma_sqd_w (where
    ! sigma_sqd_w is located on momentum levels).

    a1 = 1.0_core_rknd / ( 1.0_core_rknd - sigma_sqd_w )

    ! Interpolate a_1 from momentum levels to thermodynamic
    ! levels.  This will be used for the w'^3 turbulent advection
    ! (ta) and turbulent production (tp) combined term.
    a1_zt  = max( zm2zt( a1 ), zero_threshold )   ! Positive definite quantity

    ! Compute the explicit portion of the w'^2 and w'^3 equations.
    ! Build the right-hand side vector.
    call wp23_rhs( dt, wp2, wp3, a1, a1_zt, &
                   a3, a3_zt, wp3_on_wp2, wpthvp, wp2thvp, um, vm,  & 
                   upwp, vpwp, up2, vp2, Kw1, Kw8, Kh_zt,  & 
                   Skw_zt, tau1m, tauw3t, C1_Skw_fnc, &
                   C11_Skw_fnc, rho_ds_zm, invrs_rho_ds_zt, radf, &
                   thv_ds_zm, thv_ds_zt, l_crank_nich_diff, &
                   rhs )

    if (l_gmres) then
      call wp23_gmres( dt, wp2, wm_zm, wm_zt, a1, a1_zt, a3, a3_zt, &
                       wp3_on_wp2, &
                       Kw1, Kw8, Skw_zt, tau1m, tauw3t, C1_Skw_fnc, &
                       C11_Skw_fnc, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
                       invrs_rho_ds_zt, l_crank_nich_diff, nsup, nsub, nrhs, &
                       rhs, &
                       solut, err_code )
    else
      ! Compute the implicit portion of the w'^2 and w'^3 equations.
      ! Build the left-hand side matrix.
      call wp23_lhs( dt, wp2, wm_zm, wm_zt, a1, a1_zt, a3, a3_zt,  &
                     wp3_on_wp2, &
                     Kw1, Kw8, Skw_zt, tau1m, tauw3t, C1_Skw_fnc, &
                     C11_Skw_fnc, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
                     invrs_rho_ds_zt, l_crank_nich_diff, nsub, nsup,  & 
                     lhs )

      ! Solve the system with LAPACK
      if ( l_stats_samp .and. iwp23_matrix_condt_num > 0 ) then

        ! Perform LU decomp and solve system (LAPACK with diagnostics)
        ! Note that this can change the answer slightly
        call band_solvex( "wp2_wp3", nsup, nsub, 2*gr%nz, nrhs, & 
                          lhs, rhs, solut, rcond, err_code )

        ! Est. of the condition number of the w'^2/w^3 LHS matrix
        call stat_update_var_pt( iwp23_matrix_condt_num, 1, 1.0_core_rknd / rcond, sfc )

      else
        ! Perform LU decomp and solve system (LAPACK)
        call band_solve( "wp2_wp3", nsup, nsub, 2*gr%nz, nrhs, & 
                         lhs, rhs, solut, err_code )
      end if

    end if ! l_gmres

    ! Copy result into output arrays and clip

    do k = 1, gr%nz

      km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nz )

      k_wp3 = 2*k - 1
      k_wp2 = 2*k

      ! wp2_n(k) = wp2(k) ! For the positive definite scheme

      wp2(k) = solut(k_wp2)
      wp3(k) = solut(k_wp3)

    end do

    if (l_stats_samp) then

      ! Finalize implicit contributions for wp2

      do k = 2, gr%nz-1

        km1 = max( k-1, 1 )
        km2 = max( k-2, 1 )
        kp1 = min( k+1, gr%nz )
        kp2 = min( k+2, gr%nz )

        ! w'^2 term dp1 has both implicit and explicit components;
        ! call stat_end_update_pt.
        call stat_end_update_pt( iwp2_dp1, k, & 
           zmscr01(k) * wp2(k), zm )

        ! w'^2 term dp2 has both implicit and explicit components (if the
        ! Crank-Nicholson scheme is selected); call stat_end_update_pt.  
        ! If Crank-Nicholson diffusion is not selected, then w'^3 term dp1 is 
        ! completely implicit; call stat_update_var_pt.
        if ( l_crank_nich_diff ) then
           call stat_end_update_pt( iwp2_dp2, k, &
              zmscr02(k) * wp2(km1) & 
            + zmscr03(k) * wp2(k) & 
            + zmscr04(k) * wp2(kp1), zm )
        else
           call stat_update_var_pt( iwp2_dp2, k, &
              zmscr02(k) * wp2(km1) & 
            + zmscr03(k) * wp2(k) & 
            + zmscr04(k) * wp2(kp1), zm )
        endif

        ! w'^2 term ta is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwp2_ta, k, & 
           zmscr05(k) * wp3(k) & 
         + zmscr06(k) * wp3(kp1), zm )

        ! w'^2 term ma is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwp2_ma, k, & 
           zmscr07(k) * wp2(km1) & 
         + zmscr08(k) * wp2(k) & 
         + zmscr09(k) * wp2(kp1), zm )

        ! w'^2 term ac is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwp2_ac, k,  & 
           zmscr10(k) * wp2(k), zm )

        ! w'^2 term pr1 has both implicit and explicit components;
        ! call stat_end_update_pt.
        if ( l_tke_aniso ) then
          call stat_end_update_pt( iwp2_pr1, k, & 
             zmscr12(k) * wp2(k), zm )
        endif

        ! w'^2 term pr2 has both implicit and explicit components;
        ! call stat_end_update_pt.
        call stat_end_update_pt( iwp2_pr2, k, & 
           zmscr11(k) * wp2(k), zm )

        ! w'^2 term 4hd is completely implicit; call stat_update_var_pt.
        if ( l_hyper_dfsn ) then
           call stat_update_var_pt( iwp2_4hd, k, &
              zmscr13(k) * wp2(km2) &
            + zmscr14(k) * wp2(km1) &
            + zmscr15(k) * wp2(k) &
            + zmscr16(k) * wp2(kp1) &
            + zmscr17(k) * wp2(kp2), zm )
        endif
      enddo

      ! Finalize implicit contributions for wp3

      do k = 2, gr%nz-1, 1

        km1 = max( k-1, 1 )
        km2 = max( k-2, 1 )
        kp1 = min( k+1, gr%nz )
        kp2 = min( k+2, gr%nz )

        ! w'^3 term pr1 has both implicit and explicit components; 
        ! call stat_end_update_pt.
        call stat_end_update_pt( iwp3_pr1, k, & 
           ztscr01(k) * wp3(k), zt )

        ! w'^3 term dp1 has both implicit and explicit components (if the
        ! Crank-Nicholson scheme is selected); call stat_end_update_pt.  
        ! If Crank-Nicholson diffusion is not selected, then w'^3 term dp1 is 
        ! completely implicit; call stat_update_var_pt.
        if ( l_crank_nich_diff ) then
           call stat_end_update_pt( iwp3_dp1, k, & 
              ztscr02(k) * wp3(km1) & 
            + ztscr03(k) * wp3(k) & 
            + ztscr04(k) * wp3(kp1), zt )
        else
           call stat_update_var_pt( iwp3_dp1, k, & 
              ztscr02(k) * wp3(km1) & 
            + ztscr03(k) * wp3(k) & 
            + ztscr04(k) * wp3(kp1), zt )
        endif

        ! w'^3 term ta has both implicit and explicit components; 
        ! call stat_end_update_pt.
        call stat_end_update_pt( iwp3_ta, k, & 
           ztscr05(k) * wp3(km1) & 
         + ztscr06(k) * wp2(km1) & 
         + ztscr07(k) * wp3(k) & 
         + ztscr08(k) * wp2(k) & 
         + ztscr09(k) * wp3(kp1), zt )

        ! w'^3 term tp has both implicit and explicit components; 
        ! call stat_end_update_pt.
        call stat_end_update_pt( iwp3_tp, k,  & 
           ztscr10(k) * wp2(km1) & 
         + ztscr11(k) * wp2(k), zt )

        ! w'^3 term ma is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwp3_ma, k, & 
           ztscr12(k) * wp3(km1) & 
         + ztscr13(k) * wp3(k) & 
         + ztscr14(k) * wp3(kp1), zt )

        ! w'^3 term ac is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwp3_ac, k, & 
           ztscr15(k) * wp3(k), zt )

        ! w'^3 term pr2 has both implicit and explicit components; 
        ! call stat_end_update_pt.
        call stat_end_update_pt( iwp3_pr2, k, & 
           ztscr16(k) * wp3(k), zt )

        ! w'^3 term 4hd is completely implicit; call stat_update_var_pt.
        if ( l_hyper_dfsn ) then
           call stat_update_var_pt( iwp3_4hd, k, &
              ztscr17(k) * wp3(km2) &
            + ztscr18(k) * wp3(km1) &
            + ztscr19(k) * wp3(k) &
            + ztscr20(k) * wp3(kp1) &
            + ztscr21(k) * wp3(kp2), zt )
        endif
      enddo

    endif ! l_stats_samp


    if ( l_stats_samp ) then
      ! Store previous value for effect of the positive definite scheme
      call stat_begin_update( iwp2_pd, wp2 / real( dt, kind = core_rknd ), zm )
    endif

    if ( l_hole_fill .and. any( wp2 < w_tol_sqd ) ) then

      ! Use a simple hole filling algorithm
      call fill_holes_driver( 2, w_tol_sqd, "zm", &
                              rho_ds_zt, rho_ds_zm, &
                              wp2 )

    endif ! wp2

    ! Here we attempt to clip extreme values of wp2 to prevent a crash of the
    ! type found on the Climate Process Team ticket #49.  Chris Golaz found that
    ! instability caused by large wp2 in CLUBB led unrealistic results in AM3.
    ! -dschanen 11 Apr 2011
    where ( wp2 > 1000._core_rknd ) wp2 = 1000._core_rknd

    if ( l_stats_samp ) then
      ! Store updated value for effect of the positive definite scheme
      call stat_end_update( iwp2_pd, wp2 / real( dt, kind = core_rknd ), zm )
    endif


    ! Clip w'^2 at a minimum threshold.
    call clip_variance( clip_wp2, dt, w_tol_sqd, wp2 )

    ! Interpolate w'^2 from momentum levels to thermodynamic levels.
    ! This is used for the clipping of w'^3 according to the value
    ! of Sk_w now that w'^2 and w'^3 have been advanced one timestep.
    wp2_zt = max( zm2zt( wp2 ), w_tol_sqd )   ! Positive definite quantity

    ! Clip w'^3 by limiting skewness.
    call clip_skewness( dt, sfc_elevation, wp2_zt, wp3 )

    ! Compute wp3_zm for output purposes
    wp3_zm = zt2zm( wp3 )

    return
  end subroutine wp23_solve

  subroutine wp23_gmres( dt, wp2, wm_zm, wm_zt, a1, a1_zt, a3, a3_zt, &
                         wp3_on_wp2, &
                         Kw1, Kw8, Skw_zt, tau1m, tauw3t, C1_Skw_fnc, &
                         C11_Skw_fnc, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
                         invrs_rho_ds_zt, l_crank_nich_diff, nsup, nsub, nrhs, &
                         rhs, &
                         solut, err_code )
    ! Description:
    ! Perform all GMRES-specific matrix generation and solving for the
    ! wp2/wp3 matrices.
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        gr  ! Variable(s) 

    use clubb_precision, only:  & 
        time_precision, &  ! Variable(s)
        core_rknd

#ifdef MKL
    use error_code, only: &
      fatal_error ! Procedure(s)

    use stats_variables, only:  & 
        iwp23_matrix_condt_num, & ! Variable(s)
        l_stats_samp, & 
        sfc

    use constants_clubb, only: & 
        fstderr         ! Variable(s)

    use lapack_wrap, only:  & 
        band_solve,  & ! Procedure(s) 
        band_solvex

    use stats_type, only: & 
        stat_update_var_pt ! Procedure(s)

    use csr_matrix_class, only: &
        csr_intlc_5b_5b_ia, & ! Variables
        csr_intlc_5b_5b_ja, &
        intlc_5d_5d_ja_size

    use gmres_wrap, only: &
        gmres_solve   ! Subroutine

    use gmres_cache, only: &
        gmres_cache_soln, & ! Subroutine
        gmres_prev_soln, &        ! Variables
        gmres_prev_precond_a, &
        l_gmres_soln_ok, &
        gmres_idx_wp2wp3, &
        gmres_temp_intlc, &
        gmres_tempsize_intlc
#endif /* MKL */

    implicit none

    ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
      dt                 ! Timestep                                  [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      wp2                ! w'^2 (momentum levels)                    [m^2/s^2]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  & 
      wm_zm,           & ! w wind component on momentum levels       [m/s]
      wm_zt,           & ! w wind component on thermodynamic levels  [m/s]
      a1,              & ! a_1 (momentum levels); See eqn. 23 in `Equations for CLUBB' [-]
      a1_zt,           & ! a_1 interpolated to thermodynamic levels                    [-]
      a3,              & ! a_3 (momentum levels); See eqn. 25 in `Equations for CLUBB' [-]
      a3_zt,           & ! a_3 interpolated to thermodynamic levels  [-]
      wp3_on_wp2,      & ! Smoothed version of wp3 / wp2             [m/s]
      Kw1,             & ! Coefficient of eddy diffusivity for w'^2  [m^2/s]
      Kw8,             & ! Coefficient of eddy diffusivity for w'^3  [m^2/s]
      Skw_zt,          & ! Skewness of w on thermodynamic levels     [-]
      tau1m,           & ! Time-scale tau on momentum levels         [s]
      tauw3t,          & ! Time-scale tau on thermodynamic levels    [s]
      C1_Skw_fnc,      & ! C_1 parameter with Sk_w applied           [-]
      C11_Skw_fnc,     & ! C_11 parameter with Sk_w applied          [-]
      rho_ds_zm,       & ! Dry, static density on momentum levels    [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ momentum levs. [m^3/kg]
      invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs.  [m^3/kg]

    logical, intent(in) :: & 
      l_crank_nich_diff  ! Turns on/off Crank-Nicholson diffusion.

    integer, intent(in) :: &
      nsub,   & ! Number of subdiagonals in the LHS matrix.
      nsup,   & ! Number of superdiagonals in the LHS matrix.
      nrhs      ! Number of right-hand side vectors
                ! (GMRES currently only supports 1)

    ! Input/Output variables
    real( kind = core_rknd ), dimension(2*gr%nz), intent(inout) :: &
      rhs       ! Right hand side vector

    ! Output variables
    real( kind = core_rknd ), dimension(2*gr%nz), intent(out) :: &
      solut     ! Solution to band diagonal system

    integer, intent(out) :: err_code ! Have any errors occured?

#ifdef MKL
    ! Local variables
    real( kind = core_rknd ), dimension(nsup+nsub+1,2*gr%nz) :: &
      lhs, &    ! Implicit contributions to wp2/wp3 (band diag. matrix)
      lhs_cache ! Backup cache of LHS matrix

    real( kind = core_rknd ), dimension(intlc_5d_5d_ja_size) :: &
      lhs_a_csr ! Implicit contributions to wp2/wp3 (CSR format)

    real( kind = core_rknd ), dimension(2*gr%nz) :: &
      rhs_cache ! Backup cache of RHS vector

    real( kind = core_rknd )::  & 
      rcond  ! Est. of the reciprocal of the condition #

    ! Begin code

    if (nsup > 2) then
      write (fstderr, *) "WARNING: CSR-format solvers currently do not", &
                         "support solving with hyper diffusion", &
                         "at this time. l_hyper_dfsn ignored."
    end if
    call wp23_lhs_csr( dt, wp2, wm_zm, wm_zt, a1, a1_zt, a3, a3_zt,  &
                       wp3_on_wp2, &
                       Kw1, Kw8, Skw_zt, tau1m, tauw3t, C1_Skw_fnc, &
                       C11_Skw_fnc, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
                       invrs_rho_ds_zt, l_crank_nich_diff, & 
                       lhs_a_csr )

    if ( .not. l_gmres_soln_ok(gmres_idx_wp2wp3) ) then
      call wp23_lhs( dt, wp2, wm_zm, wm_zt, a1, a1_zt, a3, a3_zt,  &
                     wp3_on_wp2, &
                     Kw1, Kw8, Skw_zt, tau1m, tauw3t, C1_Skw_fnc, &
                     C11_Skw_fnc, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
                     invrs_rho_ds_zt, l_crank_nich_diff, nsub, nsup,  & 
                     lhs )

      ! Solve system with LAPACK to give us our first solution vector
        lhs_cache = lhs
        rhs_cache = rhs
        call band_solve( "wp2_wp3", nsup, nsub, 2*gr%nz, nrhs, &
                         lhs, rhs, solut, err_code )

        ! Use gmres_cache_wp2wp3_soln to set cache this solution for GMRES
        call gmres_cache_soln( gr%nz * 2, gmres_idx_wp2wp3, solut )
        lhs = lhs_cache
        rhs = rhs_cache
    end if ! .not. l_gmres_soln_ok(gmres_idx_wp2wp3)

    call gmres_solve( intlc_5d_5d_ja_size, (gr%nz * 2), &
                      lhs_a_csr, csr_intlc_5b_5b_ia, csr_intlc_5b_5b_ja, &
                      gmres_tempsize_intlc, &
                      gmres_prev_soln(:,gmres_idx_wp2wp3), &
                      gmres_prev_precond_a(:,gmres_idx_wp2wp3), rhs, &
                      gmres_temp_intlc, &
                      solut, err_code )
    ! Fall back to LAPACK if GMRES returned any errors
    if ( fatal_error( err_code ) ) then
      write(fstderr,*) "Errors encountered in GMRES solve."
      write(fstderr,*) "Falling back to LAPACK solver."

      ! Generate the LHS in LAPACK format
      call wp23_lhs( dt, wp2, wm_zm, wm_zt, a1, a1_zt, a3, a3_zt,  &
                     wp3_on_wp2, &
                     Kw1, Kw8, Skw_zt, tau1m, tauw3t, C1_Skw_fnc, &
                     C11_Skw_fnc, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
                     invrs_rho_ds_zt, l_crank_nich_diff, nsub, nsup,  & 
                     lhs )

      ! Note: The RHS does not need to be re-generated.

      ! Solve the system with LAPACK as a fall-back.
      if ( l_stats_samp .and. iwp23_matrix_condt_num > 0 ) then

        ! Perform LU decomp and solve system (LAPACK with diagnostics)
        ! Note that this can change the answer slightly
        call band_solvex( "wp2_wp3", nsup, nsub, 2*gr%nz, nrhs, & 
                          lhs, rhs, solut, rcond, err_code )

        ! Est. of the condition number of the w'^2/w^3 LHS matrix
        call stat_update_var_pt( iwp23_matrix_condt_num, 1, 1.0_core_rknd / rcond, sfc )

      else
        ! Perform LU decomp and solve system (LAPACK)
        call band_solve( "wp2_wp3", nsup, nsub, 2*gr%nz, nrhs, & 
                         lhs, rhs, solut, err_code )
      end if

    end if ! fatal_error

#else
    stop "This build was not compiled with PARDISO/GMRES support."

    ! These prevent compiler warnings when -DMKL not set.
    if ( l_crank_nich_diff .or. .true. ) print *, "This should be unreachable"
    solut = rhs
    solut(1:gr%nz) = a1
    solut(1:gr%nz) = a1_zt
    solut(1:gr%nz) = a3
    solut(1:gr%nz) = a3_zt
    solut(1:gr%nz) = C11_Skw_fnc
    solut(1:gr%nz) = C1_Skw_fnc
    solut(1:gr%nz) = invrs_rho_ds_zm
    solut(1:gr%nz) = invrs_rho_ds_zt
    solut(1:gr%nz) = rho_ds_zm
    solut(1:gr%nz) = rho_ds_zt
    solut(1:gr%nz) = Kw1
    solut(1:gr%nz) = Kw8
    solut(1:gr%nz) = Skw_zt
    solut(1:gr%nz) = tau1m
    solut(1:gr%nz) = tauw3t
    solut(1:gr%nz) = wm_zt
    solut(1:gr%nz) = wm_zm
    solut(1:gr%nz) = wp2
    solut(1:gr%nz) = wp3_on_wp2
    err_code = int( dt )
    err_code = nsup
    err_code = nsub
    err_code = nrhs

#endif /* MKL */

  end subroutine wp23_gmres

  !=============================================================================
  subroutine wp23_lhs( dt, wp2, wm_zm, wm_zt, a1, a1_zt, a3, a3_zt,  &
                       wp3_on_wp2, &
                       Kw1, Kw8, Skw_zt, tau1m, tauw3t, C1_Skw_fnc, &
                       C11_Skw_fnc, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
                       invrs_rho_ds_zt, l_crank_nich_diff, nsub, nsup,  & 
                       lhs )

    ! Description:
    ! Compute LHS band diagonal matrix for w'^2 and w'^3.
    ! This subroutine computes the implicit portion 
    ! of the w'^2 and w'^3 equations.
    !
    ! NOTE: If changes are made to this subroutine, ensure that the CSR
    !   version of the subroutine is updated as well! If the two are different,
    !   the results will be inconsistent between LAPACK and PARDISO/GMRES!

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        gr ! Variable

    use parameters_tunable, only:  & 
        C4,  & ! Variables
        C5,  & 
        C8,  & 
        C8b, & 
        C12, & 
        nu1_vert_res_dep, & 
        nu8_vert_res_dep, &
        nu_hd_vert_res_dep

    use constants_clubb, only:  & 
        eps,          & ! Variable(s)
        three_halves, &
        gamma_over_implicit_ts

    use model_flags, only: & 
        l_tke_aniso, & ! Variable(s)
        l_hyper_dfsn

    use diffusion, only: & 
        diffusion_zm_lhs,  & ! Procedures
        diffusion_zt_lhs

    use mean_adv, only: & 
        term_ma_zm_lhs,  & ! Procedures
        term_ma_zt_lhs

    use hyper_diffusion_4th_ord, only:  &
        hyper_dfsn_4th_ord_zm_lhs,  &
        hyper_dfsn_4th_ord_zt_lhs

    use clubb_precision, only: &
        time_precision, &
        core_rknd

    use stats_variables, only: & 
        zmscr01,    &
        zmscr02,    &
        zmscr03,    &
        zmscr04,    &
        zmscr05,    &
        zmscr06,    &
        zmscr07,    &
        zmscr08,    &
        zmscr09,    &
        zmscr11,    & 
        zmscr10,    & 
        zmscr12,    &
        zmscr13,    &
        zmscr14,    &
        zmscr15,    &
        zmscr16,    &
        zmscr17,    &
        ztscr01,    &
        ztscr02

    use stats_variables, only: &
        ztscr03,    &
        ztscr04,    &
        ztscr05,    &
        ztscr06,    &
        ztscr07,    &
        ztscr08,    &
        ztscr09,    &
        ztscr10,    &
        ztscr11,    &
        ztscr12,    &
        ztscr13,    &
        ztscr14,    &
        ztscr15,    &
        ztscr16,    &
        ztscr17,    &
        ztscr18,    &
        ztscr19,    &
        ztscr20,    &
        ztscr21

    use stats_variables, only: & 
        l_stats_samp, & 
        iwp2_dp1, & 
        iwp2_dp2, & 
        iwp2_ta, & 
        iwp2_ma, & 
        iwp2_ac, & 
        iwp2_pr2, & 
        iwp2_pr1, &
        iwp2_4hd, &
        iwp3_ta, & 
        iwp3_tp, & 
        iwp3_ma, & 
        iwp3_ac, & 
        iwp3_pr2, & 
        iwp3_pr1, & 
        iwp3_dp1, &
        iwp3_4hd

    use advance_helper_module, only: set_boundary_conditions_lhs ! Procedure(s)

    implicit none

    ! Parameter Constants
    ! Left-hand side matrix diagonal identifiers for
    ! momentum-level variable, w'^2.
    integer, parameter ::  &
      m_kp2_mdiag = 1, & ! Momentum super-super diagonal index for w'^2.
     !m_kp2_tdiag = 2, & ! Thermodynamic super-super diagonal index for w'^2.
      m_kp1_mdiag = 3, & ! Momentum super diagonal index for w'^2.
      m_kp1_tdiag = 4, & ! Thermodynamic super diagonal index for w'^2.
      m_k_mdiag   = 5, & ! Momentum main diagonal index for w'^2.
      m_k_tdiag   = 6, & ! Thermodynamic sub diagonal index for w'^2.
      m_km1_mdiag = 7, & ! Momentum sub diagonal index for w'^2.
     !m_km1_tdiag = 8, & ! Thermodynamic sub-sub diagonal index for w'^2.
      m_km2_mdiag = 9    ! Momentum sub-sub diagonal index for w'^2.

    ! Left-hand side matrix diagonal identifiers for
    ! thermodynamic-level variable, w'^3.
    integer, parameter ::  &
      t_kp2_tdiag = 1, & ! Thermodynamic super-super diagonal index for w'^3.
     !t_kp1_mdiag = 2, & ! Momentum super-super diagonal index for w'^3.
      t_kp1_tdiag = 3, & ! Thermodynamic super diagonal index for w'^3.
     !t_k_mdiag   = 4, & ! Momentum super diagonal index for w'^3.
      t_k_tdiag   = 5, & ! Thermodynamic main diagonal index for w'^3.
     !t_km1_mdiag = 6, & ! Momentum sub diagonal index for w'^3.
      t_km1_tdiag = 7, & ! Thermodynamic sub diagonal index for w'^3.
     !t_km2_mdiag = 8, & ! Momentum sub-sub diagonal index for w'^3.
      t_km2_tdiag = 9    ! Thermodynamic sub-sub diagonal index for w'^3.

    ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
      dt                 ! Timestep length                            [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      wp2,             & ! w'^2 (momentum levels)                     [m^2/s^2]
      wm_zm,           & ! w wind component on momentum levels        [m/s]
      wm_zt,           & ! w wind component on thermodynamic levels   [m/s]
      a1,              & ! sigma_sqd_w term a_1 (momentum levels)     [-]
      a1_zt,           & ! a_1 interpolated to thermodynamic levels   [-]
      a3,              & ! sigma_sqd_w term a_3 (momentum levels)     [-]
      a3_zt,           & ! a_3 interpolated to thermodynamic levels   [-]
      wp3_on_wp2,      & ! Smoothed version of wp3 / wp2              [m/s]
      Kw1,             & ! Coefficient of eddy diffusivity for w'^2   [m^2/s]
      Kw8,             & ! Coefficient of eddy diffusivity for w'^3   [m^2/s]
      Skw_zt,          & ! Skewness of w on thermodynamic levels      [-]
      tau1m,           & ! Time-scale tau on momentum levels          [s]
      tauw3t,          & ! Time-scale tau on thermodynamic levels     [s]
      C1_Skw_fnc,      & ! C_1 parameter with Sk_w applied            [-]
      C11_Skw_fnc,     & ! C_11 parameter with Sk_w applied           [-]
      rho_ds_zm,       & ! Dry, static density on momentum levels     [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels      [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ momentum levs.  [m^3/kg]
      invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs.   [m^3/kg]

    logical, intent(in) :: & 
      l_crank_nich_diff  ! Turns on/off Crank-Nicholson diffusion.

    integer, intent(in) :: &
      nsub,   & ! Number of subdiagonals in the LHS matrix.
      nsup      ! Number of superdiagonals in the LHS matrix.

    ! Output Variable
    real( kind = core_rknd ), dimension(5-nsup:5+nsub,2*gr%nz), intent(out) ::  & 
      lhs ! Implicit contributions to wp2/wp3 (band diag. matrix)

    ! Local Variables

    ! Array indices
    integer :: k, km1, km2, kp1, kp2, k_wp2, k_wp3, k_wp2_low, k_wp2_high, &
               k_wp3_low, k_wp3_high

    real( kind = core_rknd ), dimension(5) :: tmp


    ! Initialize the left-hand side matrix to 0.
    lhs = 0.0_core_rknd

    do k = 2, gr%nz-1, 1

      ! Define indices

      km1 = max( k-1, 1 )
      km2 = max( k-2, 1 )
      kp1 = min( k+1, gr%nz )
      kp2 = min( k+2, gr%nz )

      k_wp3 = 2*k - 1
      k_wp2 = 2*k


      !!!!!***** w'^2 *****!!!!!

      ! w'^2: Left-hand side (implicit w'^2 portion of the code).
      !
      ! Momentum sub-sub diagonal (lhs index: m_km2_mdiag)
      !         [ x wp2(k-2,<t+1>) ]
      ! Thermodynamic sub-sub diagonal (lhs index: m_km1_tdiag)
      !         [ x wp3(k-1,<t+1>) ]
      ! Momentum sub diagonal (lhs index: m_km1_mdiag)
      !         [ x wp2(k-1,<t+1>) ]
      ! Thermodynamic sub diagonal (lhs index: m_k_tdiag)
      !         [ x wp3(k,<t+1>) ]
      ! Momentum main diagonal (lhs index: m_k_mdiag)
      !         [ x wp2(k,<t+1>) ]
      ! Thermodynamic super diagonal (lhs index: m_kp1_tdiag)
      !         [ x wp3(k+1,<t+1>) ]
      ! Momentum super diagonal (lhs index: m_kp1_mdiag)
      !         [ x wp2(k+1,<t+1>) ]
      ! Thermodynamic super-super diagonal (lhs index: m_kp2_tdiag)
      !         [ x wp3(k+2,<t+1>) ]
      ! Momentum super-super diagonal (lhs index: m_kp2_mdiag)
      !         [ x wp2(k+2,<t+1>) ]

      ! LHS time tendency.
      lhs(m_k_mdiag,k_wp2) & 
      = + 1.0_core_rknd / real( dt, kind = core_rknd )

      ! LHS mean advection (ma) term.
      lhs((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/),k_wp2) & 
      = lhs((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/),k_wp2) & 
      + term_ma_zm_lhs( wm_zm(k), gr%invrs_dzm(k), k )

      ! LHS turbulent advection (ta) term.
      lhs((/m_kp1_tdiag,m_k_tdiag/),k_wp2) & 
      = lhs((/m_kp1_tdiag,m_k_tdiag/),k_wp2) & 
      + wp2_term_ta_lhs( rho_ds_zt(kp1), rho_ds_zt(k), &
                         invrs_rho_ds_zm(k), gr%invrs_dzm(k) )

      ! LHS accumulation (ac) term and pressure term 2 (pr2).
      lhs(m_k_mdiag,k_wp2) & 
      = lhs(m_k_mdiag,k_wp2) & 
      + wp2_terms_ac_pr2_lhs( C5, wm_zt(kp1), wm_zt(k), gr%invrs_dzm(k)  )

      ! LHS dissipation term 1 (dp1).
      ! Note:  An "over-implicit" weighted time step is applied to this term.
      !        A weighting factor of greater than 1 may be used to make the term
      !        more numerically stable (see note below for w'^3 LHS turbulent
      !        advection (ta) and turbulent production (tp) terms).
      lhs(m_k_mdiag,k_wp2)  & 
      = lhs(m_k_mdiag,k_wp2)  &
      + gamma_over_implicit_ts  & 
      * wp2_term_dp1_lhs( C1_Skw_fnc(k), tau1m(k) )

      ! LHS eddy diffusion term: dissipation term 2 (dp2).
      if ( l_crank_nich_diff ) then
        ! Eddy diffusion for wp2 using a Crank-Nicholson time step.
        lhs((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/),k_wp2) & 
        = lhs((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/),k_wp2) & 
        + (1.0_core_rknd/2.0_core_rknd) & 
        * diffusion_zm_lhs( Kw1(k), Kw1(kp1), nu1_vert_res_dep, & 
                            gr%invrs_dzt(kp1), gr%invrs_dzt(k), &
                            gr%invrs_dzm(k), k )
      else
        ! Eddy diffusion for wp2 using a completely implicit time step.
        lhs((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/),k_wp2) & 
        = lhs((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/),k_wp2) & 
        + diffusion_zm_lhs( Kw1(k), Kw1(kp1), nu1_vert_res_dep, & 
                            gr%invrs_dzt(kp1), gr%invrs_dzt(k), &
                            gr%invrs_dzm(k), k )
      endif

      ! LHS pressure term 1 (pr1).
      ! Note:  An "over-implicit" weighted time step is applied to this term.
      !        A weighting factor of greater than 1 may be used to make the term
      !        more numerically stable (see note below for w'^3 LHS turbulent
      !        advection (ta) and turbulent production (tp) terms).
      if ( l_tke_aniso ) then
        ! Add in this term if we're not assuming tke = 1.5 * wp2
        lhs(m_k_mdiag,k_wp2)  & 
        = lhs(m_k_mdiag,k_wp2)  &
        + gamma_over_implicit_ts  & 
        * wp2_term_pr1_lhs( C4, tau1m(k) )
      endif

      ! LHS 4th-order hyper-diffusion (4hd).
      if ( l_hyper_dfsn ) then
         ! Note:  w'^2 uses fixed-point boundary conditions.
         lhs( (/m_kp2_mdiag,m_kp1_mdiag,m_k_mdiag,m_km1_mdiag,m_km2_mdiag/), &
              k_wp2 )  &
         = lhs( (/m_kp2_mdiag,m_kp1_mdiag,m_k_mdiag,m_km1_mdiag,m_km2_mdiag/), &
                k_wp2 )  &
         + hyper_dfsn_4th_ord_zm_lhs( 'fixed-point', nu_hd_vert_res_dep, gr%invrs_dzm(k),  &
                                      gr%invrs_dzt(kp1), gr%invrs_dzt(k),     &
                                      gr%invrs_dzm(kp1), gr%invrs_dzm(km1),   &
                                      gr%invrs_dzt(kp2), gr%invrs_dzt(km1), k )
      endif

      if ( l_stats_samp ) then

        ! Statistics: implicit contributions for wp2.

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note below for w'^3 LHS
        !        turbulent advection (ta) and turbulent production (tp) terms).
        if ( iwp2_dp1 > 0 ) then
          zmscr01(k)  &
          = - gamma_over_implicit_ts  &
            * wp2_term_dp1_lhs( C1_Skw_fnc(k), tau1m(k) )
        endif

        if ( iwp2_dp2 > 0 ) then
          if ( l_crank_nich_diff ) then
            ! Eddy diffusion for wp2 using a Crank-Nicholson time step.
            tmp(1:3) & 
            = (1.0_core_rknd/2.0_core_rknd) & 
            * diffusion_zm_lhs( Kw1(k), Kw1(kp1), nu1_vert_res_dep, & 
                                gr%invrs_dzt(kp1), gr%invrs_dzt(k), &
                                gr%invrs_dzm(k), k )
          else
            ! Eddy diffusion for wp2 using a completely implicit time step.
            tmp(1:3) & 
            = diffusion_zm_lhs( Kw1(k), Kw1(kp1), nu1_vert_res_dep, & 
                                gr%invrs_dzt(kp1), gr%invrs_dzt(k), &
                                gr%invrs_dzm(k), k )
          endif

          zmscr02(k) = -tmp(3)
          zmscr03(k) = -tmp(2)
          zmscr04(k) = -tmp(1)

        endif

        if ( iwp2_ta > 0 ) then
          tmp(1:2) =  & 
          + wp2_term_ta_lhs( rho_ds_zt(kp1), rho_ds_zt(k), &
                             invrs_rho_ds_zm(k), gr%invrs_dzm(k) )
          zmscr05(k) = -tmp(2)
          zmscr06(k) = -tmp(1)
        endif

        if ( iwp2_ma > 0 ) then
          tmp(1:3) = & 
          + term_ma_zm_lhs( wm_zm(k), gr%invrs_dzm(k), k )
          zmscr07(k) = -tmp(3)
          zmscr08(k) = -tmp(2)
          zmscr09(k) = -tmp(1)
        endif

        ! Note:  To find the contribution of w'^2 term ac, substitute 0 for the
        !        C_5 input to function wp2_terms_ac_pr2_lhs.
        if ( iwp2_ac > 0 ) then
          zmscr10(k) =  & 
          - wp2_terms_ac_pr2_lhs( 0.0_core_rknd, wm_zt(kp1), wm_zt(k), gr%invrs_dzm(k)  )
        endif

        ! Note:  To find the contribution of w'^2 term pr2, add 1 to the
        !        C_5 input to function wp2_terms_ac_pr2_lhs.
        if ( iwp2_pr2 > 0 ) then
          zmscr11(k) =  & 
          - wp2_terms_ac_pr2_lhs( (1.0_core_rknd+C5), wm_zt(kp1), wm_zt(k),  & 
                                  gr%invrs_dzm(k)  )
        endif

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note below for w'^3 LHS
        !        turbulent advection (ta) and turbulent production (tp) terms).
        if ( iwp2_pr1 > 0 .and. l_tke_aniso ) then
          zmscr12(k)  &
          = - gamma_over_implicit_ts  &
            * wp2_term_pr1_lhs( C4, tau1m(k) )
        endif

        if ( iwp2_4hd > 0 .and. l_hyper_dfsn ) then
          tmp(1:5) = &
          hyper_dfsn_4th_ord_zm_lhs( 'fixed-point', nu_hd_vert_res_dep, gr%invrs_dzm(k),  &
                                     gr%invrs_dzt(kp1), gr%invrs_dzt(k),     &
                                     gr%invrs_dzm(kp1), gr%invrs_dzm(km1),   &
                                     gr%invrs_dzt(kp2), gr%invrs_dzt(km1), k )
          zmscr13(k) = -tmp(5)
          zmscr14(k) = -tmp(4)
          zmscr15(k) = -tmp(3)
          zmscr16(k) = -tmp(2)
          zmscr17(k) = -tmp(1)
        endif

      endif



      !!!!!***** w'^3 *****!!!!!

      ! w'^3: Left-hand side (implicit w'^3 portion of the code).
      !
      ! Thermodynamic sub-sub diagonal (lhs index: t_km2_tdiag)
      !         [ x wp3(k-2,<t+1>) ]
      ! Momentum sub-sub diagonal (lhs index: t_km2_mdiag)
      !         [ x wp2(k-2,<t+1>) ]
      ! Thermodynamic sub diagonal (lhs index: t_km1_tdiag)
      !         [ x wp3(k-1,<t+1>) ]
      ! Momentum sub diagonal (lhs index: t_km1_mdiag)
      !         [ x wp2(k-1,<t+1>) ]
      ! Thermodynamic main diagonal (lhs index: t_k_tdiag)
      !         [ x wp3(k,<t+1>) ]
      ! Momentum super diagonal (lhs index: t_k_mdiag)
      !         [ x wp2(k,<t+1>) ]
      ! Thermodynamic super diagonal (lhs index: t_kp1_tdiag)
      !         [ x wp3(k+1,<t+1>) ]
      ! Momentum super-super diagonal (lhs index: t_kp1_mdiag)
      !         [ x wp2(k+1,<t+1>) ]
      ! Thermodynamic super-super diagonal (lhs index: t_kp2_tdiag)
      !         [ x wp3(k+2,<t+1>) ]

      ! LHS time tendency.
      lhs(t_k_tdiag,k_wp3) & 
      =  + 1.0_core_rknd / real( dt, kind = core_rknd )

      ! LHS mean advection (ma) term.
      lhs((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/),k_wp3) & 
      = lhs((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/),k_wp3) & 
      + term_ma_zt_lhs( wm_zt(k), gr%invrs_dzt(k), k, gr%invrs_dzm(k), gr%invrs_dzm(k-1) )

      ! LHS turbulent advection (ta) and turbulent production (tp) terms.
      ! Note:  An "over-implicit" weighted time step is applied to these terms.
      !        The weight of the implicit portion of these terms is controlled
      !        by the factor gamma_over_implicit_ts (abbreviated "gamma" in the
      !        expression below).  A factor is added to the right-hand side of
      !        the equation in order to balance a weight that is not equal to 1,
      !        such that:
      !             -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
      !        where X is the variable that is being solved for in a predictive
      !        equation (w'^3 in this case), y(t) is the linearized portion of
      !        the terms that gets treated implicitly, and RHS is the portion of
      !        the terms that is always treated explicitly.  A weight of greater
      !        than 1 can be applied to make the terms more numerically stable.
      lhs(t_kp1_tdiag:t_km1_tdiag,k_wp3)  & 
      = lhs(t_kp1_tdiag:t_km1_tdiag,k_wp3)  &
      + gamma_over_implicit_ts  &
      * wp3_terms_ta_tp_lhs( wp2(k), wp2(km1),  &
                             a1(k), a1_zt(k), a1(km1),  &
                             a3(k), a3_zt(k), a3(km1),  &
                             wp3_on_wp2(k), wp3_on_wp2(km1), &
                             rho_ds_zm(k), rho_ds_zm(km1),  &
                             invrs_rho_ds_zt(k),  &
                             three_halves,  &
                             gr%invrs_dzt(k), k )

      ! LHS accumulation (ac) term and pressure term 2 (pr2).
      lhs(t_k_tdiag,k_wp3) & 
      = lhs(t_k_tdiag,k_wp3) & 
      + wp3_terms_ac_pr2_lhs( C11_Skw_fnc(k), & 
                              wm_zm(k), wm_zm(km1), gr%invrs_dzt(k) )

      ! LHS pressure term 1 (pr1).
      ! Note:  An "over-implicit" weighted time step is applied to this term.
      lhs(t_k_tdiag,k_wp3)  &
      = lhs(t_k_tdiag,k_wp3)  &
      + gamma_over_implicit_ts  &
      * wp3_term_pr1_lhs( C8, C8b, tauw3t(k), Skw_zt(k) )

      ! LHS eddy diffusion term: dissipation term 1 (dp1).
      !  Added a new constant, C12.
      !  Initially, this new constant will be set to 1.0 -dschanen 9/19/05
      if ( l_crank_nich_diff ) then
        ! Eddy diffusion for wp3 using a Crank-Nicholson time step.
        lhs((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/),k_wp3) & 
        = lhs((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/),k_wp3) & 
        + C12 * (1.0_core_rknd/2.0_core_rknd) & 
        * diffusion_zt_lhs( Kw8(k), Kw8(km1), nu8_vert_res_dep, & 
                            gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                            gr%invrs_dzt(k), k )
      else
        ! Eddy diffusion for wp3 using a completely implicit time step.
        lhs((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/),k_wp3) & 
        = lhs((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/),k_wp3) & 
        + C12  & 
        * diffusion_zt_lhs( Kw8(k), Kw8(km1), nu8_vert_res_dep, & 
                            gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                            gr%invrs_dzt(k), k )
      endif

      ! LHS 4th-order hyper-diffusion (4hd).
      if ( l_hyper_dfsn ) then
         ! Note:  w'^3 uses fixed-point boundary conditions.
         lhs( (/t_kp2_tdiag,t_kp1_tdiag,t_k_tdiag,t_km1_tdiag,t_km2_tdiag/), &
              k_wp3 )  &
         = lhs( (/t_kp2_tdiag,t_kp1_tdiag,t_k_tdiag,t_km1_tdiag,t_km2_tdiag/), &
                k_wp3 )  &
         + hyper_dfsn_4th_ord_zt_lhs( 'fixed-point', nu_hd_vert_res_dep, gr%invrs_dzt(k),  &
                                      gr%invrs_dzm(k), gr%invrs_dzm(km1),     &
                                      gr%invrs_dzt(kp1), gr%invrs_dzt(km1),   &
                                      gr%invrs_dzm(kp1), gr%invrs_dzm(km2), k )
      endif

      if ( l_stats_samp ) then

        ! Statistics: implicit contributions for wp3.

        ! Note:  To find the contribution of w'^3 term ta, add 3 to all of 
        !        the a_3 inputs and substitute 0 for the three_halves input to
        !        function wp3_terms_ta_tp_lhs.
        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note above for LHS turbulent
        !        advection (ta) and turbulent production (tp) terms).
        if ( iwp3_ta > 0 ) then
          tmp(1:5)  &
          = gamma_over_implicit_ts  &
          * wp3_terms_ta_tp_lhs( wp2(k), wp2(km1),  &
                                 a1(k), a1_zt(k), a1(km1),  &
                                 a3(k)+3.0_core_rknd, a3_zt(k)+3.0_core_rknd, &
                                 a3(km1)+3.0_core_rknd,  &
                                 wp3_on_wp2(k), wp3_on_wp2(km1), &
                                 rho_ds_zm(k), rho_ds_zm(km1),  &
                                 invrs_rho_ds_zt(k),  &
                                 0.0_core_rknd,  &
                                 gr%invrs_dzt(k), k )
          ztscr05(k) = -tmp(5)
          ztscr06(k) = -tmp(4)
          ztscr07(k) = -tmp(3)
          ztscr08(k) = -tmp(2)
          ztscr09(k) = -tmp(1)
        endif

        ! Note:  To find the contribution of w'^3 term tp, substitute 0 for all
        !        of the a_1 and a_3 inputs and subtract 3 from all of the a_3
        !        inputs to function wp3_terms_ta_tp_lhs.
        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note above for LHS turbulent
        !        advection (ta) and turbulent production (tp) terms).
        if ( iwp3_tp > 0 ) then
          tmp(1:5)  &
          = gamma_over_implicit_ts  &
          * wp3_terms_ta_tp_lhs( wp2(k), wp2(km1),  &
                                 0.0_core_rknd, 0.0_core_rknd, 0.0_core_rknd,  &
                                 0.0_core_rknd-3.0_core_rknd, 0.0_core_rknd-3.0_core_rknd, &
                                 0.0_core_rknd-3.0_core_rknd,  &
                                 0.0_core_rknd, 0.0_core_rknd, &
                                 rho_ds_zm(k), rho_ds_zm(km1),  &
                                 invrs_rho_ds_zt(k),  &
                                 three_halves,  &
                                 gr%invrs_dzt(k), k )
          ztscr10(k) = -tmp(4)
          ztscr11(k) = -tmp(2)
        endif

        if ( iwp3_ma > 0 ) then
          tmp(1:3) = & 
          term_ma_zt_lhs( wm_zt(k), gr%invrs_dzt(k), k, gr%invrs_dzm(k), gr%invrs_dzm(km1) )
          ztscr12(k) = -tmp(3)
          ztscr13(k) = -tmp(2)
          ztscr14(k) = -tmp(1)
        endif

        ! Note:  To find the contribution of w'^3 term ac, substitute 0 for the
        !        C_ll skewness function input to function wp3_terms_ac_pr2_lhs.
        if ( iwp3_ac > 0 ) then
          ztscr15(k) =  & 
          - wp3_terms_ac_pr2_lhs( 0.0_core_rknd, & 
                                  wm_zm(k), wm_zm(km1), gr%invrs_dzt(k) )
        endif

        ! Note:  To find the contribution of w'^3 term pr2, add 1 to the
        !        C_ll skewness function input to function wp3_terms_ac_pr2_lhs.
        if ( iwp3_pr2 > 0 ) then
          ztscr16(k) = & 
          - wp3_terms_ac_pr2_lhs( (1.0_core_rknd+C11_Skw_fnc(k)), & 
                                  wm_zm(k), wm_zm(km1), gr%invrs_dzt(k) )
        endif

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note above for LHS turbulent
        !        advection (ta) and turbulent production (tp) terms).
        if ( iwp3_pr1 > 0 ) then
          ztscr01(k)  &
          = - gamma_over_implicit_ts  &
            * wp3_term_pr1_lhs( C8, C8b, tauw3t(k), Skw_zt(k) )
        endif

        if ( iwp3_dp1 > 0 ) then
          if ( l_crank_nich_diff ) then
            ! Eddy diffusion for wp3 using a Crank-Nicholson time step.
            tmp(1:3) & 
            = C12 * (1.0_core_rknd/2.0_core_rknd) & 
            * diffusion_zt_lhs( Kw8(k), Kw8(km1), nu8_vert_res_dep, & 
                                gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                                gr%invrs_dzt(k), k )
          else
            ! Eddy diffusion for wp3 using a completely implicit time step.
            tmp(1:3) & 
            = C12  & 
            * diffusion_zt_lhs( Kw8(k), Kw8(km1), nu8_vert_res_dep, & 
                                gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                                gr%invrs_dzt(k), k )
          endif

          ztscr02(k) = -tmp(3)
          ztscr03(k) = -tmp(2)
          ztscr04(k) = -tmp(1)

        endif

        if ( iwp3_4hd > 0 .and. l_hyper_dfsn ) then
          tmp(1:5) = &
          hyper_dfsn_4th_ord_zt_lhs( 'fixed-point', nu_hd_vert_res_dep, gr%invrs_dzt(k),  &
                                     gr%invrs_dzm(k), gr%invrs_dzm(km1),     &
                                     gr%invrs_dzt(kp1), gr%invrs_dzt(km1),   &
                                     gr%invrs_dzm(kp1), gr%invrs_dzm(km2), k )
          ztscr17(k) = -tmp(5)
          ztscr18(k) = -tmp(4)
          ztscr19(k) = -tmp(3)
          ztscr20(k) = -tmp(2)
          ztscr21(k) = -tmp(1)
        endif

      endif

    enddo ! k = 2, gr%nz-1, 1


    ! Boundary conditions

    ! Both wp2 and wp3 used fixed-point boundary conditions.
    ! Therefore, anything set in the above loop at both the upper
    ! and lower boundaries would be overwritten here.  However, the
    ! above loop does not extend to the boundary levels.  An array
    ! with a value of 1 at the main diagonal on the left-hand side
    ! and with values of 0 at all other diagonals on the left-hand
    ! side will preserve the right-hand side value at that level.
    !
    !   wp3(1)  wp2(1) ... wp3(nzmax) wp2(nzmax)
    ! [  0.0     0.0         0.0     0.0  ]
    ! [  0.0     0.0         0.0     0.0  ]
    ! [  1.0     1.0   ...   1.0     1.0  ]
    ! [  0.0     0.0         0.0     0.0  ]
    ! [  0.0     0.0         0.0     0.0  ]

    ! Lower boundary
    k = 1
    k_wp3_low = 2*k - 1
    k_wp2_low = 2*k

    ! Upper boundary
    k = gr%nz
    k_wp3_high = 2*k - 1
    k_wp2_high = 2*k

    ! t_k_tdiag and m_k_mdiag need to be adjusted because the dimensions of lhs
    ! are offset
    call set_boundary_conditions_lhs( t_k_tdiag - nsup, k_wp3_low, k_wp3_high, lhs, &
                                  m_k_mdiag - nsup, k_wp2_low, k_wp2_high)

    return

  end subroutine wp23_lhs

#ifdef MKL
  !=============================================================================
  subroutine wp23_lhs_csr( dt, wp2, wm_zm, wm_zt, a1, a1_zt, a3, a3_zt,  &
                           wp3_on_wp2, &
                           Kw1, Kw8, Skw_zt, tau1m, tauw3t, C1_Skw_fnc, &
                           C11_Skw_fnc, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
                           invrs_rho_ds_zt, l_crank_nich_diff, & 
                           lhs_a_csr )

    ! Description:
    ! Compute LHS band diagonal matrix for w'^2 and w'^3.
    ! This subroutine computes the implicit portion 
    ! of the w'^2 and w'^3 equations.
    !
    ! This version of the subroutine computes the LHS in CSR (compressed
    !   sparse row) format.
    ! NOTE: This subroutine must be kept up to date with the non CSR version
    !   of the subroutine! If the two are different, the results will be
    !   inconsistent between LAPACK and PARDISO/GMRES results!

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        gr ! Variable

    use parameters_tunable, only:  & 
        C4,  & ! Variables
        C5,  & 
        C8,  & 
        C8b, & 
        C12, & 
        nu1_vert_res_dep, & 
        nu8_vert_res_dep, &
        nu_hd_vert_res_dep

    use constants_clubb, only:  & 
        eps,          & ! Variable(s)
        three_halves, &
        gamma_over_implicit_ts

    use model_flags, only: & 
        l_tke_aniso, & ! Variable(s)
        l_hyper_dfsn

    use diffusion, only: & 
        diffusion_zm_lhs,  & ! Procedures
        diffusion_zt_lhs

    use mean_adv, only: & 
        term_ma_zm_lhs,  & ! Procedures
        term_ma_zt_lhs

    use hyper_diffusion_4th_ord, only:  &
        hyper_dfsn_4th_ord_zm_lhs,  &
        hyper_dfsn_4th_ord_zt_lhs

    use clubb_precision, only: &
        time_precision, &
        core_rknd

    use stats_variables, only: & 
        zmscr01,    &
        zmscr02,    &
        zmscr03,    &
        zmscr04,    &
        zmscr05,    &
        zmscr06,    &
        zmscr07,    &
        zmscr08,    &
        zmscr09,    &
        zmscr11,    & 
        zmscr10,    & 
        zmscr12,    &
        zmscr13,    &
        zmscr14,    &
        zmscr15,    &
        zmscr16,    &
        zmscr17,    &
        ztscr01,    &
        ztscr02

    use stats_variables, only: &
        ztscr03,    &
        ztscr04,    &
        ztscr05,    &
        ztscr06,    &
        ztscr07,    &
        ztscr08,    &
        ztscr09,    &
        ztscr10,    &
        ztscr11,    &
        ztscr12,    &
        ztscr13,    &
        ztscr14,    &
        ztscr15,    &
        ztscr16,    &
        ztscr17,    &
        ztscr18,    &
        ztscr19,    &
        ztscr20,    &
        ztscr21

    use stats_variables, only: & 
        l_stats_samp, & 
        iwp2_dp1, & 
        iwp2_dp2, & 
        iwp2_ta, & 
        iwp2_ma, & 
        iwp2_ac, & 
        iwp2_pr2, & 
        iwp2_pr1, &
        iwp2_4hd, &
        iwp3_ta, & 
        iwp3_tp, & 
        iwp3_ma, & 
        iwp3_ac, & 
        iwp3_pr2, & 
        iwp3_pr1, & 
        iwp3_dp1, &
        iwp3_4hd

    use csr_matrix_class, only: &
        intlc_5d_5d_ja_size ! Variable

    implicit none

    ! Left-hand side matrix diagonal identifiers for
    ! momentum-level variable, w'^2.
    ! These are updated for each diagonal of the matrix as the
    ! LHS of the matrix is created.
    integer ::  &
     !m_kp2_mdiag, & ! Momentum super-super diagonal index for w'^2.
     !m_kp2_tdiag, & ! Thermodynamic super-super diagonal index for w'^2.
      m_kp1_mdiag, & ! Momentum super diagonal index for w'^2.
      m_kp1_tdiag, & ! Thermodynamic super diagonal index for w'^2.
      m_k_mdiag  , & ! Momentum main diagonal index for w'^2.
      m_k_tdiag  , & ! Thermodynamic sub diagonal index for w'^2.
      m_km1_mdiag    ! Momentum sub diagonal index for w'^2.
     !m_km1_tdiag, & ! Thermodynamic sub-sub diagonal index for w'^2.
     !m_km2_mdiag    ! Momentum sub-sub diagonal index for w'^2.

    ! Left-hand side matrix diagonal identifiers for
    ! thermodynamic-level variable, w'^3.
    ! These are updated for each diagonal of the matrix as the
    ! LHS of the matrix is created
    integer ::  &
     !t_kp2_tdiag, & ! Thermodynamic super-super diagonal index for w'^3.
     !t_kp1_mdiag, & ! Momentum super-super diagonal index for w'^3.
      t_kp1_tdiag, & ! Thermodynamic super diagonal index for w'^3.
     !t_k_mdiag  , & ! Momentum super diagonal index for w'^3.
      t_k_tdiag  , & ! Thermodynamic main diagonal index for w'^3.
     !t_km1_mdiag, & ! Momentum sub diagonal index for w'^3.
      t_km1_tdiag    ! Thermodynamic sub diagonal index for w'^3.
     !t_km2_mdiag, & ! Momentum sub-sub diagonal index for w'^3.
     !t_km2_tdiag    ! Thermodynamic sub-sub diagonal index for w'^3.

    ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
      dt                 ! Timestep length                            [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      wp2,             & ! w'^2 (momentum levels)                     [m^2/s^2]
      wm_zm,           & ! w wind component on momentum levels        [m/s]
      wm_zt,           & ! w wind component on thermodynamic levels   [m/s]
      a1,              & ! sigma_sqd_w term a_1 (momentum levels)     [-]
      a1_zt,           & ! a_1 interpolated to thermodynamic levels   [-]
      a3,              & ! sigma_sqd_w term a_3 (momentum levels)     [-]
      a3_zt,           & ! a_3 interpolated to thermodynamic levels   [-]
      wp3_on_wp2,      & ! Smoothed version of wp3 / wp2              [m/s]
      Kw1,             & ! Coefficient of eddy diffusivity for w'^2   [m^2/s]
      Kw8,             & ! Coefficient of eddy diffusivity for w'^3   [m^2/s]
      Skw_zt,          & ! Skewness of w on thermodynamic levels      [-]
      tau1m,           & ! Time-scale tau on momentum levels          [s]
      tauw3t,          & ! Time-scale tau on thermodynamic levels     [s]
      C1_Skw_fnc,      & ! C_1 parameter with Sk_w applied            [-]
      C11_Skw_fnc,     & ! C_11 parameter with Sk_w applied           [-]
      rho_ds_zm,       & ! Dry, static density on momentum levels     [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels      [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ momentum levs.  [m^3/kg]
      invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs.   [m^3/kg]

    logical, intent(in) :: & 
      l_crank_nich_diff  ! Turns on/off Crank-Nicholson diffusion.

!    integer, intent(in) :: &
!      nsub,   & ! Number of subdiagonals in the LHS matrix.
!      nsup      ! Number of superdiagonals in the LHS matrix.

    ! Output Variable
    real( kind = core_rknd ), dimension(intlc_5d_5d_ja_size), intent(out) ::  & 
      lhs_a_csr ! Implicit contributions to wp2/wp3 (band diag. matrix)

    ! Local Variables

    ! Array indices
    integer :: k, km1, km2, kp1, kp2, k_wp2, k_wp3, wp2_cur_row, wp3_cur_row

    real( kind = core_rknd ), dimension(5) :: tmp


    ! Initialize the left-hand side matrix to 0.
    lhs_a_csr = 0.0_core_rknd

    do k = 2, gr%nz-1, 1

      ! Define indices

      km1 = max( k-1, 1 )
      km2 = max( k-2, 1 )
      kp1 = min( k+1, gr%nz )
      kp2 = min( k+2, gr%nz )

      k_wp3 = 2*k - 1
      k_wp2 = 2*k

      wp2_cur_row = ((k_wp2 - 3) * 5) + 8
      wp3_cur_row = ((k_wp3 - 3) * 5) + 8

      !!!!!***** w'^2 *****!!!!!

      ! w'^2: Left-hand side (implicit w'^2 portion of the code).
      !
      ! Momentum sub-sub diagonal (lhs index: m_km2_mdiag)
      !         [ x wp2(k-2,<t+1>) ]
      ! Thermodynamic sub-sub diagonal (lhs index: m_km1_tdiag)
      !         [ x wp3(k-1,<t+1>) ]
      ! Momentum sub diagonal (lhs index: m_km1_mdiag)
      !         [ x wp2(k-1,<t+1>) ]
      ! Thermodynamic sub diagonal (lhs index: m_k_tdiag)
      !         [ x wp3(k,<t+1>) ]
      ! Momentum main diagonal (lhs index: m_k_mdiag)
      !         [ x wp2(k,<t+1>) ]
      ! Thermodynamic super diagonal (lhs index: m_kp1_tdiag)
      !         [ x wp3(k+1,<t+1>) ]
      ! Momentum super diagonal (lhs index: m_kp1_mdiag)
      !         [ x wp2(k+1,<t+1>) ]
      ! Thermodynamic super-super diagonal (lhs index: m_kp2_tdiag)
      !         [ x wp3(k+2,<t+1>) ]
      ! Momentum super-super diagonal (lhs index: m_kp2_mdiag)
      !         [ x wp2(k+2,<t+1>) ]

      ! NOTES FOR CSR-FORMAT MATRICES
      ! The various diagonals are referenced through the following
      ! array indices:
      ! (m_kp1_mdiag, k_wp2) ==> (wp2_cur_row + 4)
      ! (m_kp1_tdiag, k_wp2) ==> (wp2_cur_row + 3)
      ! (m_k_mdiag, k_wp2) ==> (wp2_cur_row + 2)
      ! (m_k_tdiag, k_wp2) ==> (wp2_cur_row + 1)
      ! (m_km1_mdiag, k_wp2) ==> (wp2_cur_row)
      ! For readability, these values are updated here.
      ! This means that to update the CSR version of the LHS subroutine,
      ! all that must be done is remove the ,k_wp2 from the array indices,
      ! as the CSR-format matrix is one-dimensional.

      ! NOTE: All references to lhs will need to be changed to lhs_a_csr

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! WARNING: If you have array indices that go from m_kp1_mdiag to
      ! m_km1_mdiag, you will need to set it to span by -1. This is because
      ! in the CSR-format arrays, the indices descend as you go from m_kp1_mdiag
      ! to m_km1_mdiag!
      !
      ! EXAMPLE: lhs((m_kp1_mdiag:m_km1_mdiag),wp2) would become
      !          lhs_a_csr((m_kp1_mdiag:m_km1_mdiag:-1))
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      m_kp1_mdiag = wp2_cur_row + 4
      m_kp1_tdiag = wp2_cur_row + 3
      m_k_mdiag = wp2_cur_row + 2
      m_k_tdiag = wp2_cur_row + 1
      m_km1_mdiag = wp2_cur_row

      ! LHS time tendency.
      lhs_a_csr(m_k_mdiag) & 
      = real( + 1.0_core_rknd / dt )

      ! LHS mean advection (ma) term.
      lhs_a_csr((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/)) & 
      = lhs_a_csr((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/)) & 
      + term_ma_zm_lhs( wm_zm(k), gr%invrs_dzm(k), k )

      ! LHS turbulent advection (ta) term.
      lhs_a_csr((/m_kp1_tdiag,m_k_tdiag/)) & 
      = lhs_a_csr((/m_kp1_tdiag,m_k_tdiag/)) & 
      + wp2_term_ta_lhs( rho_ds_zt(kp1), rho_ds_zt(k), &
                         invrs_rho_ds_zm(k), gr%invrs_dzm(k) )

      ! LHS accumulation (ac) term and pressure term 2 (pr2).
      lhs_a_csr(m_k_mdiag) & 
      = lhs_a_csr(m_k_mdiag) & 
      + wp2_terms_ac_pr2_lhs( C5, wm_zt(kp1), wm_zt(k), gr%invrs_dzm(k)  )

      ! LHS dissipation term 1 (dp1).
      ! Note:  An "over-implicit" weighted time step is applied to this term.
      !        A weighting factor of greater than 1 may be used to make the term
      !        more numerically stable (see note below for w'^3 LHS turbulent
      !        advection (ta) and turbulent production (tp) terms).
      lhs_a_csr(m_k_mdiag)  & 
      = lhs_a_csr(m_k_mdiag)  &
      + gamma_over_implicit_ts  & 
      * wp2_term_dp1_lhs( C1_Skw_fnc(k), tau1m(k) )

      ! LHS eddy diffusion term: dissipation term 2 (dp2).
      if ( l_crank_nich_diff ) then
        ! Eddy diffusion for wp2 using a Crank-Nicholson time step.
        lhs_a_csr((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/)) & 
        = lhs_a_csr((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/)) & 
        + (1.0_core_rknd/2.0_core_rknd) & 
        * diffusion_zm_lhs( Kw1(k), Kw1(kp1), nu1_vert_res_dep, & 
                            gr%invrs_dzt(kp1), gr%invrs_dzt(k), &
                            gr%invrs_dzm(k), k )
      else
        ! Eddy diffusion for wp2 using a completely implicit time step.
        lhs_a_csr((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/)) & 
        = lhs_a_csr((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/)) & 
        + diffusion_zm_lhs( Kw1(k), Kw1(kp1), nu1_vert_res_dep, & 
                            gr%invrs_dzt(kp1), gr%invrs_dzt(k), &
                            gr%invrs_dzm(k), k )
      endif

      ! LHS pressure term 1 (pr1).
      ! Note:  An "over-implicit" weighted time step is applied to this term.
      !        A weighting factor of greater than 1 may be used to make the term
      !        more numerically stable (see note below for w'^3 LHS turbulent
      !        advection (ta) and turbulent production (tp) terms).
      if ( l_tke_aniso ) then
        ! Add in this term if we're not assuming tke = 1.5 * wp2
        lhs_a_csr(m_k_mdiag)  & 
        = lhs_a_csr(m_k_mdiag)  &
        + gamma_over_implicit_ts  & 
        * wp2_term_pr1_lhs( C4, tau1m(k) )
      endif

      ! LHS 4th-order hyper-diffusion (4hd).
      ! NOTE: 4th-order hyper-diffusion is not yet supported in CSR-format.
      ! As such, this needs to remain commented out.
      !if ( l_hyper_dfsn ) then
      !   ! Note:  w'^2 uses fixed-point boundary conditions.
      !   lhs( (/m_kp2_mdiag,m_kp1_mdiag,m_k_mdiag,m_km1_mdiag,m_km2_mdiag/), &
      !        k_wp2) &
      !   = lhs( (/m_kp2_mdiag,m_kp1_mdiag,m_k_mdiag,m_km1_mdiag,m_km2_mdiag/), &
      !          k_wp2) &
      !   + hyper_dfsn_4th_ord_zm_lhs( 'fixed-point', nu_hd_vert_res_dep, gr%invrs_dzm(k),  &
      !                                gr%invrs_dzt(kp1), gr%invrs_dzt(k),     &
      !                                gr%invrs_dzm(kp1), gr%invrs_dzm(km1),   &
      !                                gr%invrs_dzt(kp2), gr%invrs_dzt(km1), k )
      !endif

      if ( l_stats_samp ) then

        ! Statistics: implicit contributions for wp2.

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note below for w'^3 LHS
        !        turbulent advection (ta) and turbulent production (tp) terms).
        if ( iwp2_dp1 > 0 ) then
          zmscr01(k)  &
          = - gamma_over_implicit_ts  &
            * wp2_term_dp1_lhs( C1_Skw_fnc(k), tau1m(k) )
        endif

        if ( iwp2_dp2 > 0 ) then
          if ( l_crank_nich_diff ) then
            ! Eddy diffusion for wp2 using a Crank-Nicholson time step.
            tmp(1:3) & 
            = (1.0_core_rknd/2.0_core_rknd) & 
            * diffusion_zm_lhs( Kw1(k), Kw1(kp1), nu1_vert_res_dep, & 
                              gr%invrs_dzt(kp1), gr%invrs_dzt(k), &
                              gr%invrs_dzm(k), k )
          else
            ! Eddy diffusion for wp2 using a completely implicit time step.
            tmp(1:3) & 
            = diffusion_zm_lhs( Kw1(k), Kw1(kp1), nu1_vert_res_dep, & 
                                gr%invrs_dzt(kp1), gr%invrs_dzt(k), &
                                gr%invrs_dzm(k), k )
          endif

          zmscr02(k) = -tmp(3)
          zmscr03(k) = -tmp(2)
          zmscr04(k) = -tmp(1)

        endif

        if ( iwp2_ta > 0 ) then
          tmp(1:2) =  & 
          + wp2_term_ta_lhs( rho_ds_zt(kp1), rho_ds_zt(k), &
                             invrs_rho_ds_zm(k), gr%invrs_dzm(k) )
          zmscr05(k) = -tmp(2)
          zmscr06(k) = -tmp(1)
        endif

        if ( iwp2_ma > 0 ) then
          tmp(1:3) = & 
          + term_ma_zm_lhs( wm_zm(k), gr%invrs_dzm(k), k )
          zmscr07(k) = -tmp(3)
          zmscr08(k) = -tmp(2)
          zmscr09(k) = -tmp(1)
        endif

        ! Note:  To find the contribution of w'^2 term ac, substitute 0 for the
        !        C_5 input to function wp2_terms_ac_pr2_lhs.
        if ( iwp2_ac > 0 ) then
          zmscr10(k) =  & 
          - wp2_terms_ac_pr2_lhs( 0.0_core_rknd, wm_zt(kp1), wm_zt(k), gr%invrs_dzm(k)  )
        endif

        ! Note:  To find the contribution of w'^2 term pr2, add 1 to the
        !        C_5 input to function wp2_terms_ac_pr2_lhs.
        if ( iwp2_pr2 > 0 ) then
          zmscr11(k) =  & 
          - wp2_terms_ac_pr2_lhs( (1.0_core_rknd+C5), wm_zt(kp1), wm_zt(k),  & 
                                  gr%invrs_dzm(k)  )
        endif

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note below for w'^3 LHS
        !        turbulent advection (ta) and turbulent production (tp) terms).
        if ( iwp2_pr1 > 0 .and. l_tke_aniso ) then
          zmscr12(k)  &
          = - gamma_over_implicit_ts  &
            * wp2_term_pr1_lhs( C4, tau1m(k) )
        endif

        if ( iwp2_4hd > 0 .and. l_hyper_dfsn ) then
          tmp(1:5) = &
          hyper_dfsn_4th_ord_zm_lhs( 'fixed-point', nu_hd_vert_res_dep, gr%invrs_dzm(k),  &
                                     gr%invrs_dzt(kp1), gr%invrs_dzt(k), &
                                     gr%invrs_dzm(kp1), gr%invrs_dzm(km1), &
                                     gr%invrs_dzt(kp2), gr%invrs_dzt(km1), k )
          zmscr13(k) = -tmp(5)
          zmscr14(k) = -tmp(4)
          zmscr15(k) = -tmp(3)
          zmscr16(k) = -tmp(2)
          zmscr17(k) = -tmp(1)
        endif

      endif



      !!!!!***** w'^3 *****!!!!!

      ! w'^3: Left-hand side (implicit w'^3 portion of the code).
      !
      ! Thermodynamic sub-sub diagonal (lhs index: t_km2_tdiag)
      !         [ x wp3(k-2,<t+1>) ]
      ! Momentum sub-sub diagonal (lhs index: t_km2_mdiag)
      !         [ x wp2(k-2,<t+1>) ]
      ! Thermodynamic sub diagonal (lhs index: t_km1_tdiag)
      !         [ x wp3(k-1,<t+1>) ]
      ! Momentum sub diagonal (lhs index: t_km1_mdiag)
      !         [ x wp2(k-1,<t+1>) ]
      ! Thermodynamic main diagonal (lhs index: t_k_tdiag)
      !         [ x wp3(k,<t+1>) ]
      ! Momentum super diagonal (lhs index: t_k_mdiag)
      !         [ x wp2(k,<t+1>) ]
      ! Thermodynamic super diagonal (lhs index: t_kp1_tdiag)
      !         [ x wp3(k+1,<t+1>) ]
      ! Momentum super-super diagonal (lhs index: t_kp1_mdiag)
      !         [ x wp2(k+1,<t+1>) ]
      ! Thermodynamic super-super diagonal (lhs index: t_kp2_tdiag)
      !         [ x wp3(k+2,<t+1>) ]

      ! NOTES FOR CSR-FORMAT MATRICES
      ! The various diagonals are referenced through the following
      ! array indices:
      ! (t_kp1_tdiag, k_wp3) ==> (wp3_cur_row + 4)
      ! (t_kp1_mdiag, k_wp3) ==> (wp3_cur_row + 3)
      ! (t_k_tdiag, k_wp3) ==> (wp3_cur_row + 2)
      ! (t_k_mdiag, k_wp3) ==> (wp3_cur_row + 1)
      ! (t_km1_tdiag, k_wp3) ==> (wp3_cur_row)
      ! For readability, these values are updated here.
      ! This means that to update the CSR version of the LHS subroutine,
      ! all that must be done is remove the ,k_wp2 from the array indices,
      ! as the CSR-format matrix is one-dimensional.

      ! NOTE: All references to lhs will need to be changed to lhs_a_csr

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! WARNING: If you have array indices that go from t_kp1_tdiag to
      ! t_km1_tdiag, you will need to set it to span by -1. This is because
      ! in the CSR-format arrays, the indices descend as you go from t_kp1_tdiag
      ! to t_km1_tdiag!
      !
      ! EXAMPLE: lhs((t_kp1_tdiag:t_km1_tdiag),wp3) would become
      !          lhs_a_csr((t_kp1_tdiag:t_km1_tdiag:-1))
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      t_kp1_tdiag = wp3_cur_row + 4
      !t_kp1_mdiag = wp3_cur_row + 3
      t_k_tdiag = wp3_cur_row + 2
      !t_k_mdiag = wp3_cur_row + 1
      t_km1_tdiag = wp3_cur_row

      ! LHS time tendency.
      lhs_a_csr(t_k_tdiag) & 
      = real( + 1.0_core_rknd / dt )

      ! LHS mean advection (ma) term.
      lhs_a_csr((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/)) & 
      = lhs_a_csr((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/)) & 
      + term_ma_zt_lhs( wm_zt(k), gr%invrs_dzt(k), k, gr%invrs_dzm(k), gr%invrs_dzm(km1) )

      ! LHS turbulent advection (ta) and turbulent production (tp) terms.
      ! Note:  An "over-implicit" weighted time step is applied to these terms.
      !        The weight of the implicit portion of these terms is controlled
      !        by the factor gamma_over_implicit_ts (abbreviated "gamma" in the
      !        expression below).  A factor is added to the right-hand side of
      !        the equation in order to balance a weight that is not equal to 1,
      !        such that:
      !             -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
      !        where X is the variable that is being solved for in a predictive
      !        equation (w'^3 in this case), y(t) is the linearized portion of
      !        the terms that gets treated implicitly, and RHS is the portion of
      !        the terms that is always treated explicitly.  A weight of greater
      !        than 1 can be applied to make the terms more numerically stable.
      lhs_a_csr(t_kp1_tdiag:t_km1_tdiag:-1)  & 
      = lhs_a_csr(t_kp1_tdiag:t_km1_tdiag:-1)  &
      + gamma_over_implicit_ts  &
      * wp3_terms_ta_tp_lhs( wp2(k), wp2(km1),  &
                             a1(k), a1_zt(k), a1(km1),  &
                             a3(k), a3_zt(k), a3(km1),  &
                             wp3_on_wp2(k), wp3_on_wp2(km1), &
                             rho_ds_zm(k), rho_ds_zm(km1),  &
                             invrs_rho_ds_zt(k),  &
                             three_halves,  &
                             gr%invrs_dzt(k), k )

      ! LHS accumulation (ac) term and pressure term 2 (pr2).
      lhs_a_csr(t_k_tdiag) & 
      = lhs_a_csr(t_k_tdiag) & 
      + wp3_terms_ac_pr2_lhs( C11_Skw_fnc(k), & 
                              wm_zm(k), wm_zm(km1), gr%invrs_dzt(k) )

      ! LHS pressure term 1 (pr1).
      ! Note:  An "over-implicit" weighted time step is applied to this term.
      lhs_a_csr(t_k_tdiag)  &
      = lhs_a_csr(t_k_tdiag)  &
      + gamma_over_implicit_ts  &
      * wp3_term_pr1_lhs( C8, C8b, tauw3t(k), Skw_zt(k) )

      ! LHS eddy diffusion term: dissipation term 1 (dp1).
      !  Added a new constant, C12.
      !  Initially, this new constant will be set to 1.0 -dschanen 9/19/05
      if ( l_crank_nich_diff ) then
        ! Eddy diffusion for wp3 using a Crank-Nicholson time step.
        lhs_a_csr((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/)) & 
        = lhs_a_csr((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/)) & 
        + C12 * (1.0_core_rknd/2.0_core_rknd) & 
        * diffusion_zt_lhs( Kw8(k), Kw8(km1), nu8_vert_res_dep, & 
                            gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                            gr%invrs_dzt(k), k )
      else
        ! Eddy diffusion for wp3 using a completely implicit time step.
        lhs_a_csr((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/)) & 
        = lhs_a_csr((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/)) & 
        + C12  & 
        * diffusion_zt_lhs( Kw8(k), Kw8(km1), nu8_vert_res_dep, & 
                            gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                            gr%invrs_dzt(k), k )
      endif

      ! LHS 4th-order hyper-diffusion (4hd).
      ! NOTE: 4th-order hyper-diffusion is not yet supported in CSR-format.
      ! As such, this needs to remain commented out.
      !if ( l_hyper_dfsn ) then
      !   ! Note:  w'^3 uses fixed-point boundary conditions.
      !   lhs( (/t_kp2_tdiag,t_kp1_tdiag,t_k_tdiag,t_km1_tdiag,t_km2_tdiag/), &
      !        k_wp3) &
      !   = lhs( (/t_kp2_tdiag,t_kp1_tdiag,t_k_tdiag,t_km1_tdiag,t_km2_tdiag/), &
      !          k_wp3) &
      !   + hyper_dfsn_4th_ord_zt_lhs( 'fixed-point', nu_hd_vert_res_dep, gr%invrs_dzt(k),  &
      !                                gr%invrs_dzm(k), gr%invrs_dzm(km1),     &
      !                                gr%invrs_dzt(kp1), gr%invrs_dzt(km1),   &
      !                                gr%invrs_dzm(kp1), gr%invrs_dzm(km2), k )
      !endif

      if (l_stats_samp) then

        ! Statistics: implicit contributions for wp3.

        ! Note:  To find the contribution of w'^3 term ta, add 3 to all of 
        !        the a_3 inputs and substitute 0 for the three_halves input to
        !        function wp3_terms_ta_tp_lhs.
        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note above for LHS turbulent
        !        advection (ta) and turbulent production (tp) terms).
        if ( iwp3_ta > 0 ) then
          tmp(1:5)  &
          = gamma_over_implicit_ts  &
          * wp3_terms_ta_tp_lhs( wp2(k), wp2(km1),  &
                                 a1(k), a1_zt(k), a1(km1),  &
                                 a3(k)+3.0_core_rknd, a3_zt(k)+3.0_core_rknd, &
                                 a3(km1)+3.0_core_rknd,  &
                                 wp3_on_wp2(k), wp3_on_wp2(km1), &
                                 rho_ds_zm(k), rho_ds_zm(km1),  &
                                 invrs_rho_ds_zt(k),  &
                                 0.0_core_rknd,  &
                                 gr%invrs_dzt(k), k )
          ztscr05(k) = -tmp(5)
          ztscr06(k) = -tmp(4)
          ztscr07(k) = -tmp(3)
          ztscr08(k) = -tmp(2)
          ztscr09(k) = -tmp(1)
        endif

        ! Note:  To find the contribution of w'^3 term tp, substitute 0 for all
        !        of the a_1 and a_3 inputs and subtract 3 from all of the a_3
        !        inputs to function wp3_terms_ta_tp_lhs.
        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note above for LHS turbulent
        !        advection (ta) and turbulent production (tp) terms).
        if ( iwp3_tp > 0 ) then
          tmp(1:5)  &
          = gamma_over_implicit_ts  &
          * wp3_terms_ta_tp_lhs( wp2(k), wp2(km1),  &
                                 0.0_core_rknd, 0.0_core_rknd, 0.0_core_rknd,  &
                                 0.0_core_rknd-3.0_core_rknd, 0.0_core_rknd-3.0_core_rknd, &
                                 0.0_core_rknd-3.0_core_rknd,  &
                                 0.0_core_rknd, 0.0_core_rknd, &
                                 rho_ds_zm(k), rho_ds_zm(km1),  &
                                 invrs_rho_ds_zt(k),  &
                                 three_halves,  &
                                 gr%invrs_dzt(k), k )
          ztscr10(k) = -tmp(4)
          ztscr11(k) = -tmp(2)
        endif

        if ( iwp3_ma > 0 ) then
          tmp(1:3) = & 
          term_ma_zt_lhs( wm_zt(k), gr%invrs_dzt(k), k, gr%invrs_dzm(k), gr%invrs_dzm(km1) )
          ztscr12(k) = -tmp(3)
          ztscr13(k) = -tmp(2)
          ztscr14(k) = -tmp(1)
        endif

        ! Note:  To find the contribution of w'^3 term ac, substitute 0 for the
        !        C_ll skewness function input to function wp3_terms_ac_pr2_lhs.
        if ( iwp3_ac > 0 ) then
          ztscr15(k) =  & 
          - wp3_terms_ac_pr2_lhs( 0.0_core_rknd, & 
                                  wm_zm(k), wm_zm(km1), gr%invrs_dzt(k) )
        endif

        ! Note:  To find the contribution of w'^3 term pr2, add 1 to the
        !        C_ll skewness function input to function wp3_terms_ac_pr2_lhs.
        if ( iwp3_pr2 > 0 ) then
          ztscr16(k) = & 
          - wp3_terms_ac_pr2_lhs( (1.0_core_rknd+C11_Skw_fnc(k)), & 
                                  wm_zm(k), wm_zm(km1), gr%invrs_dzt(k) )
        endif

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note above for LHS turbulent
        !        advection (ta) and turbulent production (tp) terms).
        if ( iwp3_pr1 > 0 ) then
          ztscr01(k)  &
          = - gamma_over_implicit_ts  &
            * wp3_term_pr1_lhs( C8, C8b, tauw3t(k), Skw_zt(k) )
        endif

        if ( iwp3_dp1 > 0 ) then
          if ( l_crank_nich_diff ) then
            ! Eddy diffusion for wp3 using a Crank-Nicholson time step.
            tmp(1:3) & 
            = C12 * (1.0_core_rknd/2.0_core_rknd) & 
            * diffusion_zt_lhs( Kw8(k), Kw8(km1), nu8_vert_res_dep, & 
                                gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                                gr%invrs_dzt(k), k )
          else
            ! Eddy diffusion for wp3 using a completely implicit time step.
            tmp(1:3) & 
            = C12  & 
            * diffusion_zt_lhs( Kw8(k), Kw8(km1), nu8_vert_res_dep, & 
                                gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                                gr%invrs_dzt(k), k )
          endif

          ztscr02(k) = -tmp(3)
          ztscr03(k) = -tmp(2)
          ztscr04(k) = -tmp(1)

        endif

        if ( iwp3_4hd > 0 .and. l_hyper_dfsn ) then
          tmp(1:5) = &
          hyper_dfsn_4th_ord_zt_lhs( 'fixed-point', nu_hd_vert_res_dep, gr%invrs_dzt(k),  &
                                     gr%invrs_dzm(k), gr%invrs_dzm(km1),     &
                                     gr%invrs_dzt(kp1), gr%invrs_dzt(km1),   &
                                     gr%invrs_dzm(kp1), gr%invrs_dzm(km2), k )
          ztscr17(k) = -tmp(5)
          ztscr18(k) = -tmp(4)
          ztscr19(k) = -tmp(3)
          ztscr20(k) = -tmp(2)
          ztscr21(k) = -tmp(1)
        endif

      endif

    enddo ! k = 2, gr%nz-1, 1


    ! Boundary conditions

    ! Both wp2 and wp3 used fixed-point boundary conditions.
    ! Therefore, anything set in the above loop at both the upper
    ! and lower boundaries would be overwritten here.  However, the
    ! above loop does not extend to the boundary levels.  An array
    ! with a value of 1 at the main diagonal on the left-hand side
    ! and with values of 0 at all other diagonals on the left-hand
    ! side will preserve the right-hand side value at that level.
    !
    !   wp3(1)  wp2(1) ... wp3(nzmax) wp2(nzmax)
    ! [  0.0     0.0         0.0     0.0  ]
    ! [  0.0     0.0         0.0     0.0  ]
    ! [  1.0     1.0   ...   1.0     1.0  ]
    ! [  0.0     0.0         0.0     0.0  ]
    ! [  0.0     0.0         0.0     0.0  ]

    ! Lower boundary
    k = 1
    k_wp3 = 2*k - 1
    k_wp2 = 2*k

    wp3_cur_row = 1
    wp2_cur_row = 4

    ! w'^2
    lhs_a_csr(wp2_cur_row:wp2_cur_row + 3) = 0.0_core_rknd
    lhs_a_csr(wp2_cur_row + 1) = 1.0_core_rknd

    ! w'^3
    lhs_a_csr(wp3_cur_row:wp3_cur_row + 2) = 0.0_core_rknd
    lhs_a_csr(wp3_cur_row) = 1.0_core_rknd

    ! w'^2
    !lhs(:,k_wp2)         = 0.0_core_rknd
    !lhs(m_k_mdiag,k_wp2) = 1.0_core_rknd
    ! w'^3
    !lhs(:,k_wp3)         = 0.0_core_rknd
    !lhs(t_k_tdiag,k_wp3) = 1.0_core_rknd

    ! Upper boundary
    k = gr%nz
    k_wp3 = 2*k - 1
    k_wp2 = 2*k

    ! w'^2
    lhs_a_csr(intlc_5d_5d_ja_size - 2:intlc_5d_5d_ja_size) = 0.0_core_rknd
    lhs_a_csr(intlc_5d_5d_ja_size) = 1.0_core_rknd

    ! w'^3
    lhs_a_csr(intlc_5d_5d_ja_size - 6:intlc_5d_5d_ja_size - 3) = 0.0_core_rknd
    lhs_a_csr(intlc_5d_5d_ja_size - 4) = 1.0_core_rknd

    ! w'^2
    !lhs(:,k_wp2)         = 0.0_core_rknd
    !lhs(m_k_mdiag,k_wp2) = 1.0_core_rknd
    ! w'^3
    !lhs(:,k_wp3)         = 0.0_core_rknd
    !lhs(t_k_tdiag,k_wp3) = 1.0_core_rknd


    return
  end subroutine wp23_lhs_csr
#endif /* MKL */

  !=============================================================================
  subroutine wp23_rhs( dt, wp2, wp3, a1, a1_zt, &
                       a3, a3_zt, wp3_on_wp2, wpthvp, wp2thvp, um, vm,  & 
                       upwp, vpwp, up2, vp2, Kw1, Kw8, Kh_zt, & 
                       Skw_zt, tau1m, tauw3t, C1_Skw_fnc, &
                       C11_Skw_fnc, rho_ds_zm, invrs_rho_ds_zt, radf, &
                       thv_ds_zm, thv_ds_zt, l_crank_nich_diff, &
                       rhs )

    ! Description:
    ! Compute RHS vector for w'^2 and w'^3.
    ! This subroutine computes the explicit portion of 
    ! the w'^2 and w'^3 equations.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        gr ! Variable

    use grid_class, only:  & 
        ddzt ! Procedure

    use parameters_tunable, only:  & 
        C4,  & ! Variables
        C5,  & 
        C8,  & 
        C8b, & 
        C12, & 
        C15, & 
        nu1_vert_res_dep, & 
        nu8_vert_res_dep

    use constants_clubb, only: & 
        w_tol_sqd,     & ! Variable(s)
        eps,          &
        three_halves, &
        gamma_over_implicit_ts

    use model_flags, only:  & 
        l_tke_aniso ! Variable

    use diffusion, only: & 
        diffusion_zm_lhs,  & ! Procedures
        diffusion_zt_lhs

    use clubb_precision, only:  & 
        time_precision, & ! Variable
        core_rknd

    use stats_variables, only:  & 
        l_stats_samp, iwp2_dp1, iwp2_dp2, zm, iwp2_bp,   & ! Variable(s)
        iwp2_pr1, iwp2_pr2, iwp2_pr3, iwp3_ta, zt, & 
        iwp3_tp, iwp3_bp1, iwp3_pr2, iwp3_pr1, iwp3_dp1, iwp3_bp2

    use stats_type, only:  &
        stat_update_var_pt,  & ! Procedure(s)
        stat_begin_update_pt,  &
        stat_modify_pt

    use advance_helper_module, only: set_boundary_conditions_rhs


    implicit none

    ! Constant parameters
    logical, parameter :: &
      l_wp3_2nd_buoyancy_term = .true.

    ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
      dt                 ! Timestep length                           [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      wp2,             & ! w'^2 (momentum levels)                    [m^2/s^2]
      wp3,             & ! w'^3 (thermodynamic levels)               [m^3/s^3]
      a1,              & ! sigma_sqd_w term a_1 (momentum levels)    [-]
      a1_zt,           & ! a_1 interpolated to thermodynamic levels  [-]
      a3,              & ! sigma_sqd_w term a_3 (momentum levels)    [-]
      a3_zt,           & ! a_3 interpolated to thermodynamic levels  [-]
      wp3_on_wp2,      & ! Smoothed version of wp3 / wp2             [m/s]
      wpthvp,          & ! w'th_v' (momentum levels)                 [K m/s]
      wp2thvp,         & ! w'^2th_v' (thermodynamic levels)          [K m^2/s^2]
      um,              & ! u wind component (thermodynamic levels)   [m/s]
      vm,              & ! v wind component (thermodynamic levels)   [m/s]
      upwp,            & ! u'w' (momentum levels)                    [m^2/s^2]
      vpwp,            & ! v'w' (momentum levels)                    [m^2/s^2]
      up2,             & ! u'^2 (momentum levels)                    [m^2/s^2]
      vp2,             & ! v'^2 (momentum levels)                    [m^2/s^2]
      Kw1,             & ! Coefficient of eddy diffusivity for w'^2  [m^2/s]
      Kw8,             & ! Coefficient of eddy diffusivity for w'^3  [m^2/s]
      Kh_zt,           & ! Eddy diffusivity on thermodynamic levels  [m^2/s]
      Skw_zt,          & ! Skewness of w on thermodynamic levels     [-]
      tau1m,           & ! Time-scale tau on momentum levels         [s]
      tauw3t,          & ! Time-scale tau on thermodynamic levels    [s]
      C1_Skw_fnc,      & ! C_1 parameter with Sk_w applied           [-]
      C11_Skw_fnc,     & ! C_11 parameter with Sk_w applied          [-]
      rho_ds_zm,       & ! Dry, static density on momentum levels    [kg/m^3]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. levs.  [m^3/kg]
      radf,            & ! Buoyancy production at the CL top         [m^2/s^3]
      thv_ds_zm,       & ! Dry, base-state theta_v on momentum levs. [K]
      thv_ds_zt          ! Dry, base-state theta_v on thermo. levs.  [K]

    logical, intent(in) :: & 
      l_crank_nich_diff   ! Turns on/off Crank-Nicholson diffusion.

    ! Output Variable
    real( kind = core_rknd ), dimension(2*gr%nz), intent(out) :: & 
      rhs   ! RHS of band matrix

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: & 
      dum_dz, dvm_dz ! Vertical derivatives of um and vm

    ! Array indices
    integer :: k, km1, kp1, k_wp2, k_wp3, k_wp2_low, k_wp2_high, &
               k_wp3_low, k_wp3_high

    ! For "over-implicit" weighted time step.
    ! This vector holds output from the LHS (implicit) portion of a term at a
    ! given vertical level.  This output is weighted and applied to the RHS.
    ! This is used if the implicit portion of the term is "over-implicit", which
    ! means that the LHS contribution is given extra weight (>1) in order to
    ! increase numerical stability.  A weighted factor must then be applied to
    ! the RHS in order to balance the weight.
    real( kind = core_rknd ), dimension(5) :: lhs_fnc_output

    real( kind = core_rknd ), dimension(3) :: &
      rhs_diff ! For use in Crank-Nicholson eddy diffusion.

    real( kind = core_rknd ) :: temp


    ! Initialize the right-hand side vector to 0.
    rhs = 0.0_core_rknd

    if ( l_wp3_2nd_buoyancy_term ) then
      ! Compute the vertical derivative of the u and v winds
      dum_dz = ddzt( um )
      dvm_dz = ddzt( vm )
    else
      dum_dz = -999._core_rknd
      dvm_dz = -999._core_rknd
    end if

    do k = 2, gr%nz-1, 1


      ! Define indices

      km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nz )

      k_wp3 = 2*k - 1
      k_wp2 = 2*k


      !!!!!***** w'^2 *****!!!!!

      ! w'^2: Right-hand side (explicit w'^2 portion of the code).

      ! RHS time tendency.
      rhs(k_wp2) & 
      = + ( 1.0_core_rknd / real( dt, kind = core_rknd ) ) * wp2(k)

      ! RHS buoyancy production (bp) term and pressure term 2 (pr2).
      rhs(k_wp2) & 
      = rhs(k_wp2) & 
      + wp2_terms_bp_pr2_rhs( C5, thv_ds_zm(k), wpthvp(k) )

      ! RHS buoyancy production at CL top due to LW radiative cooling
      rhs(k_wp2) = rhs(k_wp2) + radf(k) 

      ! RHS pressure term 3 (pr3).
      rhs(k_wp2) & 
      = rhs(k_wp2) & 
      + wp2_term_pr3_rhs( C5, thv_ds_zm(k), wpthvp(k), upwp(k), um(kp1), &
                          um(k), vpwp(k), vm(kp1), vm(k), gr%invrs_dzm(k) )

      ! RHS dissipation term 1 (dp1).
      rhs(k_wp2) &
      = rhs(k_wp2) &
      + wp2_term_dp1_rhs( C1_Skw_fnc(k), tau1m(k), w_tol_sqd )

      ! RHS contribution from "over-implicit" weighted time step
      ! for LHS dissipation term 1 (dp1).
      !
      ! Note:  An "over-implicit" weighted time step is applied to this term.
      !        A weighting factor of greater than 1 may be used to make the term
      !        more numerically stable (see note below for w'^3 RHS turbulent
      !        advection (ta) and turbulent production (tp) terms).
      lhs_fnc_output(1)  &
      = wp2_term_dp1_lhs( C1_Skw_fnc(k), tau1m(k) )
      rhs(k_wp2)  &
      = rhs(k_wp2)  &
      + ( 1.0_core_rknd - gamma_over_implicit_ts )  &
      * ( - lhs_fnc_output(1) * wp2(k) )

      ! RHS eddy diffusion term: dissipation term 2 (dp2).
      if ( l_crank_nich_diff ) then
        ! These lines are for the diffusional term with a Crank-Nicholson
        ! time step.  They are not used for completely implicit diffusion.
        rhs_diff(1:3) & 
        = (1.0_core_rknd/2.0_core_rknd) & 
        * diffusion_zm_lhs( Kw1(k), Kw1(kp1), nu1_vert_res_dep, & 
                            gr%invrs_dzt(kp1), gr%invrs_dzt(k), &
                            gr%invrs_dzm(k), k )
        rhs(k_wp2)   =   rhs(k_wp2) & 
                       - rhs_diff(3) * wp2(km1) & 
                       - rhs_diff(2) * wp2(k) & 
                       - rhs_diff(1) * wp2(kp1)
      endif

      ! RHS pressure term 1 (pr1).
      if ( l_tke_aniso ) then

        rhs(k_wp2) & 
        = rhs(k_wp2) & 
        + wp2_term_pr1_rhs( C4, up2(k), vp2(k), tau1m(k) )

        ! RHS contribution from "over-implicit" weighted time step
        ! for LHS dissipation term 1 (dp1).
        !
        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note below for w'^3 RHS
        !        turbulent advection (ta) and turbulent production (tp) terms).
        lhs_fnc_output(1)  &
        = wp2_term_pr1_lhs( C4, tau1m(k) )
        rhs(k_wp2)  &
        = rhs(k_wp2)  &
        + ( 1.0_core_rknd - gamma_over_implicit_ts )  &
        * ( - lhs_fnc_output(1) * wp2(k) )

      endif

      if ( l_stats_samp ) then

        ! Statistics: explicit contributions for wp2.

        ! w'^2 term dp2 has both implicit and explicit components (if the
        ! Crank-Nicholson scheme is selected); call stat_begin_update_pt.  
        ! Since stat_begin_update_pt automatically subtracts the value sent in, 
        ! reverse the sign on right-hand side diffusion component.  If 
        ! Crank-Nicholson diffusion is not selected, the stat_begin_update_pt 
        ! will not be called.
        if ( l_crank_nich_diff ) then
          call stat_begin_update_pt( iwp2_dp2, k, & 
            rhs_diff(3) * wp2(km1) & 
          + rhs_diff(2) * wp2(k) & 
          + rhs_diff(1) * wp2(kp1), zm )
        endif

        ! w'^2 term bp is completely explicit; call stat_update_var_pt.
        ! Note:  To find the contribution of w'^2 term bp, substitute 0 for the
        !        C_5 input to function wp2_terms_bp_pr2_rhs.
        call stat_update_var_pt( iwp2_bp, k, & 
          wp2_terms_bp_pr2_rhs( 0.0_core_rknd, thv_ds_zm(k), wpthvp(k) ), zm )

        ! w'^2 term pr1 has both implicit and explicit components; call
        ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
        ! subtracts the value sent in, reverse the sign on wp2_term_pr1_rhs.
        if ( l_tke_aniso ) then
          call stat_begin_update_pt( iwp2_pr1, k, & 
            -wp2_term_pr1_rhs( C4, up2(k), vp2(k), tau1m(k) ), zm )

          ! Note:  An "over-implicit" weighted time step is applied to this
          !        term.  A weighting factor of greater than 1 may be used to
          !        make the term more numerically stable (see note below for
          !        w'^3 RHS turbulent advection (ta) and turbulent
          !        production (tp) terms).
          lhs_fnc_output(1)  &
          = wp2_term_pr1_lhs( C4, tau1m(k) )
          call stat_modify_pt( iwp2_pr1, k, &
                               + ( 1.0_core_rknd - gamma_over_implicit_ts )  &
                               * ( - lhs_fnc_output(1) * wp2(k) ), zm )
        endif

        ! w'^2 term pr2 has both implicit and explicit components; call
        ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
        ! subtracts the value sent in, reverse the sign on wp2_terms_bp_pr2_rhs.
        ! Note:  To find the contribution of w'^2 term pr2, add 1 to the
        !        C_5 input to function wp2_terms_bp_pr2_rhs.
        call stat_begin_update_pt( iwp2_pr2, k, & 
          -wp2_terms_bp_pr2_rhs( (1.0_core_rknd+C5), thv_ds_zm(k), wpthvp(k) ), zm )

        ! w'^2 term dp1 has both implicit and explicit components; call
        ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
        ! subtracts the value sent in, reverse the sign on wp2_term_dp1_rhs.
        call stat_begin_update_pt( iwp2_dp1, k, &
          -wp2_term_dp1_rhs( C1_Skw_fnc(k), tau1m(k), w_tol_sqd ), zm )

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note below for w'^3 RHS
        !        turbulent advection (ta) and turbulent production (tp) terms).
        lhs_fnc_output(1)  &
        = wp2_term_dp1_lhs( C1_Skw_fnc(k), tau1m(k) )
        call stat_modify_pt( iwp2_dp1, k, &
                             + ( 1.0_core_rknd - gamma_over_implicit_ts )  &
                             * ( - lhs_fnc_output(1) * wp2(k) ), zm )

        ! w'^2 term pr3 is completely explicit; call stat_update_var_pt.
        call stat_update_var_pt( iwp2_pr3, k, & 
          wp2_term_pr3_rhs( C5, thv_ds_zm(k), wpthvp(k), upwp(k), um(kp1), &
                            um(k), vpwp(k), vm(kp1), vm(k), gr%invrs_dzm(k) ), &
                                 zm )

      endif



      !!!!!***** w'^3 *****!!!!!

      ! w'^3: Right-hand side (explicit w'^3 portion of the code).

      ! RHS time tendency.
      rhs(k_wp3) = & 
      + ( 1.0_core_rknd / real( dt, kind = core_rknd ) * wp3(k) )

      ! RHS turbulent advection (ta) and turbulent production (tp) terms.
!     rhs(k_wp3)  & 
!     = rhs(k_wp3)  & 
!     + wp3_terms_ta_tp_rhs( wp3_zm(k), wp3_zm(km1),  &
!                            wp2(k), wp2(km1),  &
!                            a1(k), a1_zt(k), a1(km1),  &
!                            a3(k), a3_zt(k), a3(km1),  &
!                            wp3_on_wp2(k), wp3_on_wp2(km1), &
!                            rho_ds_zm(k), rho_ds_zm(km1),  &
!                            invrs_rho_ds_zt(k),  &
!                            three_halves,  &
!                            gr%invrs_dzt(k) )

      ! RHS contribution from "over-implicit" weighted time step
      ! for LHS turbulent advection (ta) and turbulent production (tp) terms.
      !
      ! Note:  An "over-implicit" weighted time step is applied to these terms.
      !        The weight of the implicit portion of these terms is controlled
      !        by the factor gamma_over_implicit_ts (abbreviated "gamma" in the
      !        expression below).  A factor is added to the right-hand side of
      !        the equation in order to balance a weight that is not equal to 1,
      !        such that:
      !             -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
      !        where X is the variable that is being solved for in a predictive
      !        equation (w'^3 in this case), y(t) is the linearized portion of
      !        the terms that gets treated implicitly, and RHS is the portion of
      !        the terms that is always treated explicitly.  A weight of greater
      !        than 1 can be applied to make the terms more numerically stable.
      lhs_fnc_output(1:5)  &
      = wp3_terms_ta_tp_lhs( wp2(k), wp2(km1),  &
                             a1(k), a1_zt(k), a1(km1),  &
                             a3(k), a3_zt(k), a3(km1),  &
                             wp3_on_wp2(k), wp3_on_wp2(km1), &
                             rho_ds_zm(k), rho_ds_zm(km1),  &
                             invrs_rho_ds_zt(k),  &
                             three_halves,  &
                             gr%invrs_dzt(k), k )
      rhs(k_wp3)  & 
      = rhs(k_wp3)  &
      + ( 1.0_core_rknd - gamma_over_implicit_ts )  &
      * ( - lhs_fnc_output(1) * wp3(kp1)  &
          - lhs_fnc_output(2) * wp2(k)  &
          - lhs_fnc_output(3) * wp3(k)  &
          - lhs_fnc_output(4) * wp2(km1)  &
          - lhs_fnc_output(5) * wp3(km1) )

      ! RHS buoyancy production (bp) term and pressure term 2 (pr2).
      rhs(k_wp3) & 
      = rhs(k_wp3) & 
      + wp3_terms_bp1_pr2_rhs( C11_Skw_fnc(k), thv_ds_zt(k), wp2thvp(k) )

      ! RHS pressure term 1 (pr1).
      rhs(k_wp3) & 
      = rhs(k_wp3) & 
      + wp3_term_pr1_rhs( C8, C8b, tauw3t(k), Skw_zt(k), wp3(k) )

      ! RHS contribution from "over-implicit" weighted time step
      ! for LHS pressure term 1 (pr1).
      !
      ! Note:  An "over-implicit" weighted time step is applied to this term.
      lhs_fnc_output(1)  &
      = wp3_term_pr1_lhs( C8, C8b, tauw3t(k), Skw_zt(k) )
      rhs(k_wp3)  & 
      = rhs(k_wp3)  &
      + ( 1.0_core_rknd - gamma_over_implicit_ts )  &
      * ( - lhs_fnc_output(1) * wp3(k) )

      ! RHS eddy diffusion term: dissipation term 1 (dp1).
      if ( l_crank_nich_diff ) then
        ! These lines are for the diffusional term with a Crank-Nicholson
        ! time step.  They are not used for completely implicit diffusion.
        rhs_diff(1:3) & 
        = C12 * (1.0_core_rknd/2.0_core_rknd) & 
        * diffusion_zt_lhs( Kw8(k), Kw8(km1), nu8_vert_res_dep, & 
                            gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                            gr%invrs_dzt(k), k )
        rhs(k_wp3)   =   rhs(k_wp3) & 
                       - rhs_diff(3) * wp3(km1) & 
                       - rhs_diff(2) * wp3(k) & 
                       - rhs_diff(1) * wp3(kp1)
      endif

      if ( l_wp3_2nd_buoyancy_term ) then
        ! RHS 2nd bouyancy term
        rhs(k_wp3) = rhs(k_wp3) &
                   + wp3_term_bp2_rhs( C15, Kh_zt(k), wpthvp(k), wpthvp(km1), &
                                       dum_dz(k), dum_dz(km1), dvm_dz(k), dvm_dz(km1), &
                                       upwp(k), upwp(km1), vpwp(k), vpwp(km1), &
                                       thv_ds_zt(k), gr%invrs_dzt(k) )
      end if

      if ( l_stats_samp ) then

        ! Statistics: explicit contributions for wp3.

        ! w'^3 term ta has both implicit and explicit components; call 
        ! stat_begin_update_pt.  Since stat_begin_update_pt automatically 
        ! subtracts the value sent in, reverse the sign on wp3_terms_ta_tp_rhs.
        ! Note:  To find the contribution of w'^3 term ta, add 3 to all of the
        !        a_3 inputs and substitute 0 for the three_halves input to
        !        function wp3_terms_ta_tp_rhs.
!       call stat_begin_update_pt( iwp3_ta, k, &
!         -wp3_terms_ta_tp_rhs( wp3_zm(k), wp3_zm(km1),  &
!                               wp2(k), wp2(km1),  &
!                               a1(k), a1_zt(k), a1(km1),  &
!                               a3(k)+3.0_core_rknd, a3_zt(k)+3.0_core_rknd, 
!                               a3(km1)+3.0_core_rknd,  &
!                               wp3_on_wp2(k), wp3_on_wp2(km1), &
!                               rho_ds_zm(k), rho_ds_zm(km1),  &
!                               invrs_rho_ds_zt(k),  &
!                               0.0_core_rknd,  &
!                               gr%invrs_dzt(k) ),  &
!                                  zt )
        call stat_begin_update_pt( iwp3_ta, k, 0.0_core_rknd, zt )

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note above for RHS turbulent
        !        advection (ta) and turbulent production (tp) terms).
        lhs_fnc_output(1:5)  &
        = wp3_terms_ta_tp_lhs( wp2(k), wp2(km1),  &
                               a1(k), a1_zt(k), a1(km1),  &
                               a3(k)+3.0_core_rknd, a3_zt(k)+3.0_core_rknd, &
                               a3(km1)+3.0_core_rknd,  &
                               wp3_on_wp2(k), wp3_on_wp2(km1), &
                               rho_ds_zm(k), rho_ds_zm(km1),  &
                               invrs_rho_ds_zt(k),  &
                               0.0_core_rknd,  &
                               gr%invrs_dzt(k), k )
        call stat_modify_pt( iwp3_ta, k,  &
                             + ( 1.0_core_rknd - gamma_over_implicit_ts )  &
                             * ( - lhs_fnc_output(1) * wp3(kp1)  &
                                 - lhs_fnc_output(2) * wp2(k)  &
                                 - lhs_fnc_output(3) * wp3(k)  &
                                 - lhs_fnc_output(4) * wp2(km1)  &
                                 - lhs_fnc_output(5) * wp3(km1) ), zt )

        ! w'^3 term tp has both implicit and explicit components; call 
        ! stat_begin_update_pt.  Since stat_begin_update_pt automatically 
        ! subtracts the value sent in, reverse the sign on wp3_terms_ta_tp_rhs.
        ! Note:  To find the contribution of w'^3 term tp, substitute 0 for all
        !        of the a_1 and a_3 inputs and subtract 3 from all of the a_3
        !        inputs to function wp3_terms_ta_tp_rhs.
!       call stat_begin_update_pt( iwp3_tp, k,  &
!         -wp3_terms_ta_tp_rhs( wp3_zm(k), wp3_zm(km1),  &
!                               wp2(k), wp2(km1),  &
!                               0.0_core_rknd, 0.0_core_rknd, 0.0_core_rknd,  &
!                               0.0_core_rknd-3.0_core_rknd, 0.0_core_rknd-3.0_core_rknd, 
!                               0.0_core_rknd-3.0_core_rknd,  &
!                               0.0_core_rknd, 0.0_core_rknd, &
!                               rho_ds_zm(k), rho_ds_zm(km1),  &
!                               invrs_rho_ds_zt(k),  &
!                               three_halves,  &
!                               gr%invrs_dzt(k) ),  &
!                                  zt )
        call stat_begin_update_pt( iwp3_tp, k,  0.0_core_rknd, zt )

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note above for RHS turbulent
        !        advection (ta) and turbulent production (tp) terms).
        lhs_fnc_output(1:5)  &
        = wp3_terms_ta_tp_lhs( wp2(k), wp2(km1),  &
                               0.0_core_rknd, 0.0_core_rknd, 0.0_core_rknd,  &
                               0.0_core_rknd-3.0_core_rknd, 0.0_core_rknd-3.0_core_rknd, &
                               0.0_core_rknd-3.0_core_rknd,  &
                               0.0_core_rknd, 0.0_core_rknd, &
                               rho_ds_zm(k), rho_ds_zm(km1), &
                               invrs_rho_ds_zt(k), &
                               three_halves, &
                               gr%invrs_dzt(k), k )
        call stat_modify_pt( iwp3_tp, k,  &
                             + ( 1.0_core_rknd - gamma_over_implicit_ts )  &
                             * ( - lhs_fnc_output(2) * wp2(k)  &
                                 - lhs_fnc_output(4) * wp2(km1) ), zt )

        ! w'^3 term bp is completely explicit; call stat_update_var_pt.
        ! Note:  To find the contribution of w'^3 term bp, substitute 0 for the
        !        C_11 skewness function input to function wp3_terms_bp1_pr2_rhs.
        call stat_update_var_pt( iwp3_bp1, k, & 
          wp3_terms_bp1_pr2_rhs( 0.0_core_rknd, thv_ds_zt(k), wp2thvp(k) ), zt )

        ! w'^3 term pr2 has both implicit and explicit components; call
        ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
        ! subtracts the value sent in, reverse the sign on wp3_terms_bp1_pr2_rhs.
        ! Note:  To find the contribution of w'^3 term pr2, add 1 to the
        !        C_11 skewness function input to function wp3_terms_bp1_pr2_rhs.
        call stat_begin_update_pt( iwp3_pr2, k, & 
          -wp3_terms_bp1_pr2_rhs( (1.0_core_rknd+C11_Skw_fnc(k)), thv_ds_zt(k), &
                                 wp2thvp(k) ), & 
                                   zt )

        ! w'^3 term pr1 has both implicit and explicit components; call 
        ! stat_begin_update_pt.  Since stat_begin_update_pt automatically 
        ! subtracts the value sent in, reverse the sign on wp3_term_pr1_rhs.
        call stat_begin_update_pt( iwp3_pr1, k, & 
          -wp3_term_pr1_rhs( C8, C8b, tauw3t(k), Skw_zt(k), wp3(k) ), & 
                                   zt )

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note above for RHS turbulent
        !        advection (ta) and turbulent production (tp) terms).
        lhs_fnc_output(1)  &
        = wp3_term_pr1_lhs( C8, C8b, tauw3t(k), Skw_zt(k) )
        call stat_modify_pt( iwp3_pr1, k,  &
                             + ( 1.0_core_rknd - gamma_over_implicit_ts )  &
                             * ( - lhs_fnc_output(1) * wp3(k) ), zt )

        ! w'^3 term dp1 has both implicit and explicit components (if the
        ! Crank-Nicholson scheme is selected); call stat_begin_update_pt.  
        ! Since stat_begin_update_pt automatically subtracts the value sent in, 
        ! reverse the sign on right-hand side diffusion component.  If 
        ! Crank-Nicholson diffusion is not selected, the stat_begin_update_pt 
        ! will not be called.
        if ( l_crank_nich_diff ) then
          call stat_begin_update_pt( iwp3_dp1, k, & 
              rhs_diff(3) * wp3(km1) & 
            + rhs_diff(2) * wp3(k) & 
            + rhs_diff(1) * wp3(kp1), zt )
        endif
                  
        if ( l_wp3_2nd_buoyancy_term ) then
          temp = wp3_term_bp2_rhs( C15, Kh_zt(k), wpthvp(k), wpthvp(km1), &
                                   dum_dz(k), dum_dz(km1), dvm_dz(k), dvm_dz(km1), &
                                   upwp(k), upwp(km1), vpwp(k), vpwp(km1), &
                                   thv_ds_zt(k), gr%invrs_dzt(k) )
          call stat_update_var_pt( iwp3_bp2, k, temp, zt )
        end if

      endif ! l_stats_samp

    enddo ! k = 2..gr%nz-1


    ! Boundary conditions

    ! Both wp2 and wp3 used fixed-point boundary conditions.
    ! Therefore, anything set in the above loop at both the upper
    ! and lower boundaries would be overwritten here.  However, the
    ! above loop does not extend to the boundary levels.  An array
    ! with a value of 1 at the main diagonal on the left-hand side
    ! and with values of 0 at all other diagonals on the left-hand
    ! side will preserve the right-hand side value at that level.

    ! Lower boundary
    k = 1
    k_wp3_low = 2*k - 1
    k_wp2_low = 2*k

    ! Upper boundary
    k = gr%nz
    k_wp3_high = 2*k - 1
    k_wp2_high = 2*k


    ! The value of w'^2 at the lower boundary will remain the same.
    ! When the lower boundary is at the surface, the surface value of
    ! w'^2 is set in subroutine surface_varnce (surface_varnce_module.F).

    ! The value of w'^3 at the lower boundary will be 0.
 
    ! The value of w'^2 at the upper boundary will be set to the threshold
    ! minimum value of w_tol_sqd.

    ! The value of w'^3 at the upper boundary will be set to 0.
    call set_boundary_conditions_rhs( &
            wp2(1), k_wp2_low, w_tol_sqd, k_wp2_high, & ! Intent(in)
            rhs, & ! Intent(inout)
            0.0_core_rknd, k_wp3_low, 0.0_core_rknd, k_wp3_high )

    return

  end subroutine wp23_rhs

  !=============================================================================
  pure function wp2_term_ta_lhs( rho_ds_ztp1, rho_ds_zt, &
                                 invrs_rho_ds_zm, invrs_dzm ) &
  result( lhs )

    ! Description:
    ! Turbulent advection term for w'^2:  implicit portion of the code.
    !
    ! The d(w'^2)/dt equation contains a turbulent advection term:
    !
    ! - (1/rho_ds) * d( rho_ds * w'^3 )/dz.
    !
    ! The term is solved for completely implicitly, such that:
    !
    ! - (1/rho_ds) * d( rho_ds * w'^3(t+1) )/dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign 
    !        is reversed and the leading "-" in front of the term is changed 
    !        to a "+".
    !
    ! The timestep index (t+1) means that the value of w'^3 being used is from 
    ! the next timestep, which is being advanced to in solving the d(w'^2)/dt 
    ! and d(w'^3)/dt equations.
    !
    ! This term is discretized as follows:
    !
    ! While the values of w'^2 are found on the momentum levels, the values of 
    ! w'^3 are found on the thermodynamic levels.  Additionally, the values of
    ! rho_ds_zt are found on the thermodynamic levels, and the values of
    ! invrs_rho_ds_zm are found on the momentum levels.  On the thermodynamic
    ! levels, the values of rho_ds_zt are multiplied by the values of w'^3.  The
    ! derivative of (rho_ds_zt * w'^3) is taken over the intermediate (central)
    ! momentum level, where it is multiplied by invrs_rho_ds_zm, yielding the
    ! desired results.
    !
    ! -----rho_ds_ztp1--------wp3p1---------------------------- t(k+1)
    !
    ! ========invrs_rho_ds_zm==========d(rho_ds*wp3)/dz======== m(k)
    !
    ! -----rho_ds_zt----------wp3------------------------------ t(k)
    !
    ! The vertical indices t(k+1), m(k), and t(k) correspond with altitudes 
    ! zt(k+1), zm(k), and zt(k), respectively.  The letter "t" is used for 
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzm(k) = 1 / ( zt(k+1) - zt(k) )

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2       ! Thermodynamic subdiagonal index.

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      rho_ds_ztp1,     & ! Dry, static density at thermo. level (k+1)  [kg/m^3]
      rho_ds_zt,       & ! Dry, static density at thermo. level (k)    [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ moment. lev. (k) [m^3/kg]
      invrs_dzm          ! Inverse of grid spacing (k)                 [1/m]

    ! Return Variable
    real( kind = core_rknd ), dimension(2) :: lhs

    ! Thermodynamic superdiagonal: [ x wp3(k+1,<t+1>) ]
    lhs(kp1_tdiag) & 
    = + invrs_rho_ds_zm * invrs_dzm * rho_ds_ztp1

    ! Thermodynamic subdiagonal: [ x wp3(k,<t+1>) ]
    lhs(k_tdiag) & 
    = - invrs_rho_ds_zm * invrs_dzm * rho_ds_zt

    return

  end function wp2_term_ta_lhs

  !=============================================================================
  pure function wp2_terms_ac_pr2_lhs( C5, wm_ztp1, wm_zt, invrs_dzm ) & 
  result( lhs )

    ! Description:
    ! Accumulation of w'^2 and w'^2 pressure term 2:  implicit portion of the 
    ! code.
    !
    ! The d(w'^2)/dt equation contains an accumulation term:
    !
    ! - 2 w'^2 dw/dz;
    !
    ! and pressure term 2:
    !
    ! - C_5 ( -2 w'^2 dw/dz + 2 (g/th_0) w'th_v' ).
    !
    ! The w'^2 accumulation term is completely implicit, while w'^2 pressure 
    ! term 2 has both implicit and explicit components.  The accumulation term 
    ! and the implicit portion of pressure term 2 are combined and solved 
    ! together as:
    !
    ! + ( 1 - C_5 ) ( -2 w'^2(t+1) dw/dz ).
    !
    ! Note:  When the term is brought over to the left-hand side, the sign 
    !        is reversed and the leading "-" in front of the "2" is changed 
    !        to a "+".
    !
    ! The timestep index (t+1) means that the value of w'^2 being used is from 
    ! the next timestep, which is being advanced to in solving the d(w'^2)/dt 
    ! equation.
    !
    ! The terms are discretized as follows:
    !
    ! The values of w'^2 are found on the momentum levels, while the values of 
    ! wm_zt (mean vertical velocity on thermodynamic levels) are found on the 
    ! thermodynamic levels.  The vertical derivative of wm_zt is taken over the 
    ! intermediate (central) momentum level.  It is then multiplied by w'^2 
    ! (implicitly calculated at timestep (t+1)) and the coefficients to yield 
    ! the desired results.
    !
    ! -------wm_ztp1------------------------------------------- t(k+1)
    !
    ! ===============d(wm_zt)/dz============wp2================ m(k)
    !
    ! -------wm_zt--------------------------------------------- t(k)
    !
    ! The vertical indices t(k+1), m(k), and t(k) correspond with altitudes 
    ! zt(k+1), zm(k), and zt(k), respectively.  The letter "t" is used for 
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzm(k) = 1 / ( zt(k+1) - zt(k) )

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C5,        & ! Model parameter C_5                            [-]
      wm_ztp1,   & ! w wind component at t:hermodynamic levels (k+1) [m/s]
      wm_zt,     & ! w wind component at thermodynamic levels (k)   [m/s]
      invrs_dzm    ! Inverse of grid spacing (k)                    [1/m]

    ! Return Variable
    real( kind = core_rknd ) :: lhs

    ! Momentum main diagonal: [ x wp2(k,<t+1>) ]
    lhs & 
    = + ( 1.0_core_rknd - C5 ) * 2.0_core_rknd * invrs_dzm * ( wm_ztp1 - wm_zt )

    return

  end function wp2_terms_ac_pr2_lhs

  !=============================================================================
  pure function wp2_term_dp1_lhs( C1_Skw_fnc, tau1m ) & 
  result( lhs )

    ! Description:
    ! Dissipation term 1 for w'^2:  implicit portion of the code.
    !
    ! The d(w'^2)/dt equation contains dissipation term 1:
    !
    ! - ( C_1 / tau_1m ) w'^2.
    !
    ! Since w'^2 has a minimum threshold, the term should be damped only to that
    ! threshold.  The term becomes:
    !
    ! - ( C_1 / tau_1m ) * ( w'^2 - threshold ).
    !
    ! This term is broken into implicit and explicit portions.  The implicit 
    ! portion of this term is:
    !
    ! - ( C_1 / tau_1m ) w'^2(t+1).
    !
    ! Note:  When the implicit term is brought over to the left-hand side, the
    !        sign is reversed and the leading "-" in front of the term is 
    !        changed to a "+".
    !
    ! The timestep index (t+1) means that the value of w'^2 being used is from 
    ! the next timestep, which is being advanced to in solving the d(w'^2)/dt 
    ! equation.
    !
    ! The values of w'^2 are found on the momentum levels.  The values of the 
    ! C_1 skewness function and time-scale tau1m are also found on the momentum 
    ! levels.

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C1_Skw_fnc,  & ! C_1 parameter with Sk_w applied (k)   [-]
      tau1m          ! Time-scale tau at momentum levels (k) [s]

    ! Return Variable
    real( kind = core_rknd ) :: lhs

    ! Momentum main diagonal: [ x wp2(k,<t+1>) ]
    lhs & 
    = + C1_Skw_fnc / tau1m

    return
  end function wp2_term_dp1_lhs

  !=============================================================================
  pure function wp2_term_pr1_lhs( C4, tau1m ) & 
  result( lhs )

    ! Description
    ! Pressure term 1 for w'^2:  implicit portion of the code.
    !
    ! The d(w'^2)/dt equation contains pressure term 1:
    !
    ! - ( C_4 / tau_1m ) * ( w'^2 - (2/3)*em ),
    !
    ! where em = (1/2) * ( w'^2 + u'^2 + v'^2 ).
    !
    ! This simplifies to:
    !
    ! - ( C_4 / tau_1m ) * (2/3) * w'^2
    !    + ( C_4 / tau_1m ) * (1/3) * ( u'^2 + v'^2 ).
    !
    ! Pressure term 1 has both implicit and explicit components.  The implicit
    ! portion is:
    !
    ! - ( C_4 / tau_1m ) * (2/3) * w'^2(t+1);
    !
    ! and is computed in this function.
    !
    ! Note:  When the implicit term is brought over to the left-hand side, the
    !        sign is reversed and the leading "-" in front of the term is 
    !        changed to a "+".
    !
    ! The timestep index (t+1) means that the value of w'^2 being used is from 
    ! the next timestep, which is being advanced to in solving the d(w'^2)/dt 
    ! equation.
    !
    ! The values of w'^2 are found on momentum levels, as are the values of tau1m.

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C4,   & ! Model parameter C_4                   [-]
      tau1m   ! Time-scale tau at momentum levels (k) [s]

    ! Return Variable
    real( kind = core_rknd ) :: lhs

    ! Momentum main diagonal: [ x wp2(k,<t+1>) ]
    lhs & 
    = + ( 2.0_core_rknd * C4 ) / ( 3.0_core_rknd * tau1m )

    return
  end function wp2_term_pr1_lhs

  !=============================================================================
  pure function wp2_terms_bp_pr2_rhs( C5, thv_ds_zm, wpthvp ) & 
  result( rhs )

    ! Description:
    ! Buoyancy production of w'^2 and w'^2 pressure term 2:  explicit portion of
    ! the code.
    !
    ! The d(w'^2)/dt equation contains a buoyancy production term:
    !
    ! + 2 (g/thv_ds) w'th_v';
    !
    ! and pressure term 2:
    !
    ! - C_5 ( -2 w'^2 dw/dz + 2 (g/thv_ds) w'th_v' ).
    !
    ! The w'^2 buoyancy production term is completely explicit, while w'^2 
    ! pressure term 2 has both implicit and explicit components.  The buoyancy 
    ! production term and the explicit portion of pressure term 2 are combined 
    ! and solved together as:
    !
    ! + ( 1 - C_5 ) ( 2 (g/thv_ds) w'th_v' ).

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use constants_clubb, only:  & 
    ! Variable(s)        
        grav ! Gravitational acceleration [m/s^2]

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C5,        & ! Model parameter C_5                             [-]
      thv_ds_zm, & ! Dry, base-state theta_v at momentum level (k)   [K]
      wpthvp       ! w'th_v'(k)                                      [K m/s]

    ! Return Variable
    real( kind = core_rknd ) :: rhs

    rhs & 
    = + ( 1.0_core_rknd - C5 ) * 2.0_core_rknd * ( grav / thv_ds_zm ) * wpthvp

    return
  end function wp2_terms_bp_pr2_rhs

  !=============================================================================
  pure function wp2_term_dp1_rhs( C1_Skw_fnc, tau1m, threshold ) & 
  result( rhs )

    ! Description:
    ! Dissipation term 1 for w'^2:  explicit portion of the code.
    !
    ! The d(w'^2)/dt equation contains dissipation term 1:
    !
    ! - ( C_1 / tau_1m ) w'^2.
    !
    ! Since w'^2 has a minimum threshold, the term should be damped only to that
    ! threshold.  The term becomes:
    !
    ! - ( C_1 / tau_1m ) * ( w'^2 - threshold ).
    !
    ! This term is broken into implicit and explicit portions.  The explicit 
    ! portion of this term is:
    !
    ! + ( C_1 / tau_1m ) * threshold.
    !
    ! The values of the C_1 skewness function, time-scale tau1m, and the 
    ! threshold are found on the momentum levels.

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C1_Skw_fnc,  & ! C_1 parameter with Sk_w applied (k)   [-]
      tau1m,       & ! Time-scale tau at momentum levels (k) [s]
      threshold      ! Minimum allowable value of w'^2       [m^2/s^2]

    ! Return Variable
    real( kind = core_rknd ) :: rhs

    rhs & 
    = + ( C1_Skw_fnc / tau1m ) * threshold

    return
  end function wp2_term_dp1_rhs

  !=============================================================================
  pure function wp2_term_pr3_rhs( C5, thv_ds_zm, wpthvp, upwp, ump1, &
                                  um, vpwp, vmp1, vm, invrs_dzm ) &
  result( rhs )

    ! Description:
    ! Pressure term 3 for w'^2:  explicit portion of the code.
    !
    ! The d(w'^2)/dt equation contains pressure term 3:
    !
    ! + (2/3) C_5 [ (g/thv_ds) w'th_v' - u'w' du/dz - v'w' dv/dz ].
    !
    ! This term is solved for completely explicitly and is discretized as 
    ! follows:
    !
    ! The values of w'th_v', u'w', and v'w' are found on the momentum levels,
    ! whereas the values of um and vm are found on the thermodynamic levels.
    ! Additionally, the values of thv_ds_zm are found on the momentum levels.
    ! The derivatives of both um and vm are taken over the intermediate
    ! (central) momentum level.  All the remaining mathematical operations take
    ! place at the central momentum level, yielding the desired result.
    !
    ! -----ump1------------vmp1-------------------------------------- t(k+1)
    !
    ! =upwp====d(um)/dz========d(vm)/dz==vpwp===thv_ds_zm==wpthvp==== m(k)
    !
    ! -----um--------------vm---------------------------------------- t(k)
    !
    ! The vertical indices t(k+1), m(k), and t(k) correspond with altitudes 
    ! zt(k+1), zm(k), and zt(k), respectively.  The letter "t" is used for 
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzm(k) = 1 / ( zt(k+1) - zt(k) )

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use constants_clubb, only: & ! Variables 
        grav, & ! Gravitational acceleration [m/s^2]
        zero_threshold

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C5,        & ! Model parameter C_5                            [-]
      thv_ds_zm, & ! Dry, base-state theta_v at momentum level (k)  [K]
      wpthvp,    & ! w'th_v'(k)                                     [K m/s]
      upwp,      & ! u'w'(k)                                        [m^2/s^2]
      ump1,      & ! um(k+1)                                        [m/s]
      um,        & ! um(k)                                          [m/s]
      vpwp,      & ! v'w'(k)                                        [m^2/s^2]
      vmp1,      & ! vm(k+1)                                        [m/s]
      vm,        & ! vm(k)                                          [m/s]
      invrs_dzm    ! Inverse of grid spacing (k)                    [1/m]

    ! Return Variable
    real( kind = core_rknd ) :: rhs

    rhs & 
    ! Michael Falk, 2 August 2007
    ! Use the following code for standard mixing, with c_k=0.548:
    = + (2.0_core_rknd/3.0_core_rknd) * C5 & 
                  * ( ( grav / thv_ds_zm ) * wpthvp & 
                      - upwp * invrs_dzm * ( ump1 - um ) & 
                      - vpwp * invrs_dzm * ( vmp1 - vm ) & 
                    )
     ! Use the following code for alternate mixing, with c_k=0.1 or 0.2
!    = + (2.0_core_rknd/3.0_core_rknd) * C5 &
!                  * ( ( grav / thv_ds_zm ) * wpthvp &
!                      - 0. * upwp * invrs_dzm * ( ump1 - um ) &
!                      - 0. * vpwp * invrs_dzm * ( vmp1 - vm ) &
!                    )
!    eMFc

    ! Added by dschanen for ticket #36
    ! We have found that when shear generation is zero this term will only be
    ! offset by hole-filling (wp2_pd) and reduces turbulence 
    ! unrealistically at lower altitudes to make up the difference.
    rhs = max( rhs, zero_threshold )

    return
  end function wp2_term_pr3_rhs

  !=============================================================================
  pure function wp2_term_pr1_rhs( C4, up2, vp2, tau1m ) & 
  result( rhs )

    ! Description:
    ! Pressure term 1 for w'^2:  explicit portion of the code.
    !
    ! The d(w'^2)/dt equation contains pressure term 1:
    !
    ! - ( C_4 / tau_1m ) * ( w'^2 - (2/3)*em );
    !
    ! where em = (1/2) * ( w'^2 + u'^2 + v'^2 ).
    !
    ! This simplifies to:
    !
    ! - ( C_4 / tau_1m ) * (2/3) * w'^2
    !    + ( C_4 / tau_1m ) * (1/3) * ( u'^2 + v'^2 ).
    !
    ! Pressure term 1 has both implicit and explicit components.
    ! The explicit portion is:
    !
    ! + ( C_4 / tau_1m ) * (1/3) * ( u'^2 + v'^2 );
    !
    ! and is computed in this function.
    !
    ! The values of u'^2 and v'^2 are found on momentum levels, as are the 
    ! values of tau1m.

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C4,   & ! Model parameter C_4                   [-]
      up2,  & ! u'^2(k)                               [m^2/s^2]
      vp2,  & ! v'^2(k)                               [m^2/s^2]
      tau1m   ! Time-scale tau at momentum levels (k) [s]

    ! Return Variable
    real( kind = core_rknd ) :: rhs

    rhs & 
    = + ( C4 * ( up2 + vp2 ) ) / ( 3.0_core_rknd * tau1m )

    return
  end function wp2_term_pr1_rhs

  !=============================================================================
  pure function wp3_terms_ta_tp_lhs( wp2, wp2m1,  &
                                     a1, a1_zt, a1m1,  &
                                     a3, a3_zt, a3m1,  &
                                     wp3_on_wp2, wp3_on_wp2_m1, &
                                     rho_ds_zm, rho_ds_zmm1,  &
                                     invrs_rho_ds_zt,  &
                                     const_three_halves,  &
                                     invrs_dzt, level )  &
  result( lhs )

    ! Description:
    ! Turbulent advection and turbulent production of w'^3:  implicit portion of
    ! the code.
    !
    ! The d(w'^3)/dt equation contains a turbulent advection term:
    !
    ! - (1/rho_ds) * d( rho_ds * w'^4 )/dz;
    !
    ! and a turbulent production term:
    !
    ! + 3 * ( w'^2 / rho_ds ) * d( rho_ds * w'^2 )/dz.
    !
    ! A substitution is made in order to close the turbulent advection term, 
    ! such that:
    !
    ! w'^4 = coef_sig_sqd_w * (w'^2)^2  +  a_1 * ( (w'^3)^2 / w'^2 );
    !
    ! where both a_1 and coef_sig_sqd_w are variables that are functions of
    ! sigma_sqd_w, such that: 
    !
    ! coef_sig_sqd_w = 3*(sigma_sqd_w)^2 + 6*(1 - sigma_sqd_w)*sigma_sqd_w
    !                  + (1 - sigma_sqd_w)^2; and
    !
    ! a_1 = 1 / (1 - sigma_sqd_w).
    !
    ! Since the turbulent advection and turbulent production terms are being
    ! combined, a further substitution is made, such that:
    !
    ! a_3 = coef_sig_sqd_w - 3;
    !
    ! and thus:
    !
    ! w'^4 = (a_3 + 3) * (w'^2)^2  +  a_1 * ( (w'^3)^2 / w'^2 ).
    !
    ! The turbulent production term is rewritten as:
    !
    ! + 3 * ( w'^2 / rho_ds ) * d[ rho_ds * w'^2 ]/dz
    ! = + (3/rho_ds) * d[ rho_ds * (w'^2)^2 ]/dz - (3/2) * d[ (w'^2)^2 ]/dz.
    !
    ! The turbulent advection and turbulent production terms are combined as:
    !
    ! - (1/rho_ds) * d [ rho_ds * a_3 * (w'^2)^2 ] / dz
    ! - (1/rho_ds) * d [ rho_ds * a_1 * ( (w'^3)^2 / w'^2 ) ] / dz
    ! - (3/2) * d [ (w'^2)^2 ] / dz.
    !
    ! The (w'^2)^2 and (w'^3)^2 terms are both linearized, such that:
    !
    ! ( w'^2(t+1) )^2 = - ( w'^2(t) )^2  +  2 * w'^2(t) * w'^2(t+1);
    ! ( w'^3(t+1) )^2 = - ( w'^3(t) )^2  +  2 * w'^3(t) * w'^3(t+1);
    !
    ! which produces implicit and explicit portions of these terms.  The 
    ! implicit portion of these terms is:
    !
    ! - (1/rho_ds) * d [ rho_ds * a_3 * 2 * w'^2(t) * w'^2(t+1) ] / dz
    ! - (1/rho_ds) * d [ rho_ds * a_1 
    !                    * ( 2 * w'^3(t) * w'^3(t+1) ) / w'^2(t) ] / dz
    ! - (3/2) * d [ 2 * w'^2(t) * w'^2(t+1) ] /dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign is
    !        reversed and the leading "-" in front of all d[ ] / dz terms is
    !        changed to a "+".
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is 
    ! being advanced to in solving the d(w'^3)/dt and d(w'^2)/dt equations.
    !
    ! The implicit portion of these terms is discretized as follows:
    !
    ! The values of w'^3 are found on the thermodynamic levels, while the values
    ! of w'^2, a_1, and a_3 are found on the momentum levels.  Additionally, the
    ! values of rho_ds_zm are found on the momentum levels, and the values of
    ! invrs_rho_ds_zt are found on the thermodynamic levels.  The variable w'^3
    ! is interpolated to the intermediate momentum levels.  The values of the
    ! mathematical expressions (called F, G, and H here) within the dF/dz,
    ! dG/dz, and dH/dz terms are computed on the momentum levels.  Then, the
    ! derivatives (d/dz) of the expressions (F, G, and H) are taken over the
    ! central thermodynamic level, where dF/dz and dG/dz are multiplied by
    ! invrs_rho_ds_zt, and where dH/dz is multiplied by 3/2.  This yields the
    ! desired results.  In this function, the values of F, G, and H are as
    ! follows:
    !
    ! F = rho_ds_zm * a_3(t) * 2 * w'^2(t) * w'^2(t+1);
    !
    ! G = rho_ds_zm * a_1(t) * ( 2 * w'^3(t) * w'^3(t+1) ) / w'^2(t); and
    !
    ! H = 2 * w'^2(t) * w'^2(t+1).
    !
    !
    ! ------------------------------------------------wp3p1-------------- t(k+1)
    !
    ! ===a3====wp2====rho_ds_zm====a1======================wp3(interp)=== m(k)
    !
    ! ---dH/dz---dF/dz----invrs_rho_ds_zt----dG/dz----wp3---------------- t(k)
    !
    ! ===a3m1==wp2m1==rho_ds_zmm1==a1m1====================wp3(interp)=== m(k-1)
    !
    ! ------------------------------------------------wp3m1-------------- t(k-1)
    !
    ! The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1) correspond 
    ! with altitudes zt(k+1), zm(k), zt(k), zm(k-1), and zt(k-1), respectively. 
    ! The letter "t" is used for thermodynamic levels and the letter "m" is 
    ! used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) )

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use grid_class, only:  &
        gr ! Variable gr%weights_zt2zm

    use constants_clubb, only:  &
        w_tol_sqd

    use model_flags, only:  &
        l_standard_term_ta

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
      k_mdiag   = 2,    & ! Momentum superdiagonal index.
      k_tdiag   = 3,    & ! Thermodynamic main diagonal index.
      km1_mdiag = 4,    & ! Momentum subdiagonal index. 
      km1_tdiag = 5       ! Thermodynamic subdiagonal index.

    integer, parameter :: & 
      t_above = 1,    & ! Index for upper thermodynamic level grid weight.
      t_below = 2       ! Index for lower thermodynamic level grid weight.

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      wp2,                & ! w'^2(k)                                  [m^2/s^2]
      wp2m1,              & ! w'^2(k-1)                                [m^2/s^2]
      a1,                 & ! a_1(k)                                   [-]
      a1_zt,              & ! a_1 interpolated to thermo. level (k)    [-]
      a1m1,               & ! a_1(k-1)                                 [-]
      a3,                 & ! a_3(k)                                   [-]
      a3_zt,              & ! a_3 interpolated to thermo. level (k)    [-]
      a3m1,               & ! a_3(k-1)                                 [-]
      wp3_on_wp2,         & ! wp3 / wp2 (k)                            [m/s]
      wp3_on_wp2_m1,      & ! wp3 / wp2 (k-1)                          [m/s]
      rho_ds_zm,          & ! Dry, static density at moment. lev (k)   [kg/m^3]
      rho_ds_zmm1,        & ! Dry, static density at moment. lev (k-1) [kg/m^3]
      invrs_rho_ds_zt,    & ! Inv dry, static density @ thermo lev (k) [m^3/kg]
      const_three_halves, & ! "3/2" ("0" is sent in for wp3_ta budget) [-]
      invrs_dzt             ! Inverse of grid spacing (k)              [1/m]

    integer, intent(in) :: & 
      level ! Central thermodynamic level (on which calculation occurs).

    ! Return Variable
    real( kind = core_rknd ), dimension(5) :: lhs

    ! Local Variables
    integer :: & 
      mk,    & ! Momentum level directly above central thermodynamic level.
      mkm1     ! Momentum level directly below central thermodynamic level.


    ! Momentum level (k) is between thermodynamic level (k+1)
    ! and thermodynamic level (k).
    mk = level

    ! Momentum level (k-1) is between thermodynamic level (k)
    ! and thermodynamic level (k-1).
    mkm1 = level - 1

    if ( l_standard_term_ta ) then

       ! The turbulent advection term is discretized normally, in accordance
       ! with the model equations found in the documentation and the description
       ! listed above.

       ! Thermodynamic superdiagonal: [ x wp3(k+1,<t+1>) ]
       lhs(kp1_tdiag) &
       = + invrs_rho_ds_zt &
           * invrs_dzt &
             * rho_ds_zm * a1 &
             * wp3_on_wp2 &
             * gr%weights_zt2zm(t_above,mk)

       ! Momentum superdiagonal: [ x wp2(k,<t+1>) ]
       lhs(k_mdiag) &
       = + invrs_rho_ds_zt &
           * invrs_dzt * rho_ds_zm * a3 * wp2 &
         + const_three_halves &
           * invrs_dzt * wp2

       ! Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
       lhs(k_tdiag) &
       = + invrs_rho_ds_zt &
           * invrs_dzt &
             * (   rho_ds_zm * a1 &
                   * wp3_on_wp2 &
                   * gr%weights_zt2zm(t_below,mk) &
                 - rho_ds_zmm1 * a1m1 &
                   * wp3_on_wp2_m1 &
                   * gr%weights_zt2zm(t_above,mkm1) &
               )

       ! Momentum subdiagonal: [ x wp2(k-1,<t+1>) ]
       lhs(km1_mdiag) &
       = - invrs_rho_ds_zt &
           * invrs_dzt * rho_ds_zmm1 * a3m1 * wp2m1 &
         - const_three_halves &
           * invrs_dzt * wp2m1

       ! Thermodynamic subdiagonal: [ x wp3(k-1,<t+1>) ]
       lhs(km1_tdiag) &
       = - invrs_rho_ds_zt &
           * invrs_dzt &
             * rho_ds_zmm1 * a1m1 &
             * wp3_on_wp2_m1 &
             * gr%weights_zt2zm(t_below,mkm1)

    else

       ! Brian tried a new discretization for the turbulent advection term, 
       ! which contains the term:
       !  - (1/rho_ds) * d [ rho_ds * a_1 * (w'^3)^2 / w'^2 ] / dz.  In order
       ! to help stabilize w'^3, a_1 has been pulled outside of the derivative. 
       ! On the left-hand side of the equation, this effects the thermodynamic 
       ! superdiagonal (kp1_tdiag), the thermodynamic main diagonal (k_tdiag), 
       ! and the thermodynamic subdiagonal (km1_tdiag).

       ! Additionally, the discretization of the turbulent advection term, which
       ! contains the term:
       !  - (1/rho_ds) * d [ rho_ds * (a_3 + 3) * (w'^2)^2 ] / dz, has been 
       ! altered to pull (a_3 + 3) outside of the derivative.  This was done in
       ! order to help stabilize w'^3.  On the left-hand side of the equation,
       ! this effects the momentum superdiagonal (k_mdiag) and the momentum 
       ! subdiagonal (km1_mdiag).

       ! Thermodynamic superdiagonal: [ x wp3(k+1,<t+1>) ]
       lhs(kp1_tdiag) & 
       = + invrs_rho_ds_zt &
           * a1_zt * invrs_dzt &
             * rho_ds_zm &
             * wp3_on_wp2 &
             * gr%weights_zt2zm(t_above,mk)

       ! Momentum superdiagonal: [ x wp2(k,<t+1>) ]
       lhs(k_mdiag) & 
       = + invrs_rho_ds_zt &
           * a3_zt * invrs_dzt * rho_ds_zm * wp2 &
         + const_three_halves &
           * invrs_dzt * wp2

       ! Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
       lhs(k_tdiag) & 
       = + invrs_rho_ds_zt &
           * a1_zt * invrs_dzt & 
             * (   rho_ds_zm &
                   * wp3_on_wp2 & 
                   * gr%weights_zt2zm(t_below,mk) & 
                 - rho_ds_zmm1 &
                   * wp3_on_wp2_m1 & 
                   * gr%weights_zt2zm(t_above,mkm1) & 
               )

       ! Momentum subdiagonal: [ x wp2(k-1,<t+1>) ]
       lhs(km1_mdiag) & 
       = - invrs_rho_ds_zt &
           * a3_zt * invrs_dzt * rho_ds_zmm1 * wp2m1 &
         - const_three_halves &
           * invrs_dzt * wp2m1

       ! Thermodynamic subdiagonal: [ x wp3(k-1,<t+1>) ]
       lhs(km1_tdiag) & 
       = - invrs_rho_ds_zt &
           * a1_zt * invrs_dzt &
             * rho_ds_zmm1 &
             * wp3_on_wp2_m1 & 
             * gr%weights_zt2zm(t_below,mkm1)

       ! End of code that pulls out a3.
       ! End of Brian's a1 change.  Feb. 14, 2008.

    end if ! l_standard_term_ta


    return
  end function wp3_terms_ta_tp_lhs

  !=============================================================================
  pure function wp3_terms_ac_pr2_lhs( C11_Skw_fnc,  & 
                                      wm_zm, wm_zmm1, invrs_dzt ) & 
  result( lhs )

    ! Description:
    ! Accumulation of w'^3 and w'^3 pressure term 2:  implicit portion of the 
    ! code.
    !
    ! The d(w'^3)/dt equation contains an accumulation term:
    !
    ! - 3 w'^3 dw/dz;
    !
    ! and pressure term 2:
    !
    ! - C_11 ( -3 w'^3 dw/dz + 3 (g/th_0) w'^2th_v' ).
    !
    ! The w'^3 accumulation term is completely implicit, while w'^3 pressure 
    ! term 2 has both implicit and explicit components.  The accumulation term 
    ! and the implicit portion of pressure term 2 are combined and solved 
    ! together as:
    !
    ! + ( 1 - C_11 ) ( -3 w'^3(t+1) dw/dz ).
    !
    ! Note:  When the term is brought over to the left-hand side, the sign 
    !        is reversed and the leading "-" in front of the "3" is changed 
    !        to a "+".
    !
    ! The timestep index (t+1) means that the value of w'^3 being used is from 
    ! the next timestep, which is being advanced to in solving the d(w'^3)/dt 
    ! equation.
    !
    ! The terms are discretized as follows:
    !
    ! The values of w'^3 are found on thermodynamic levels, while the values of
    ! wm_zm (mean vertical velocity on momentum levels) are found on momentum
    ! levels.  The vertical derivative of wm_zm is taken over the intermediate
    ! (central) thermodynamic level.  It is then multiplied by w'^3 (implicitly
    ! calculated at timestep (t+1)) and the coefficients to yield the desired
    ! results.
    !
    ! =======wm_zm============================================= m(k)
    !
    ! ---------------d(wm_zm)/dz------------wp3---------------- t(k)
    !
    ! =======wm_zmm1=========================================== m(k-1)
    !
    ! The vertical indices m(k), t(k), and m(k-1) correspond with altitudes 
    ! zm(k), zt(k), and zm(k-1), respectively.  The letter "t" is used for 
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) )

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C11_Skw_fnc,  & ! C_11 parameter with Sk_w applied (k)      [-]
      wm_zm,        & ! w wind component at momentum levels (k)   [m/s]
      wm_zmm1,      & ! w wind component at momentum levels (k-1) [m/s]
      invrs_dzt       ! Inverse of grid spacing (k)               [1/m]

    ! Return Variable
    real( kind = core_rknd ) :: lhs

    ! Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
    lhs & 
    = + ( 1.0_core_rknd - C11_Skw_fnc ) & 
        * 3.0_core_rknd * invrs_dzt * ( wm_zm - wm_zmm1 )

    return
  end function wp3_terms_ac_pr2_lhs

  !=============================================================================
  pure function wp3_term_pr1_lhs( C8, C8b, tauw3t, Skw_zt ) & 
  result( lhs )

    ! Description:
    ! Pressure term 1 for w'^3:  implicit portion of the code.
    !
    ! Pressure term 1 is the term:
    !
    ! - (C_8/tau_w3t) * ( C_8b * Sk_wt^4 + 1 ) * w'^3;
    !
    ! where Sk_wt = w'^3 / (w'^2)^(3/2).
    !
    ! This term needs to be linearized, so function L(w'^3) is defined to be 
    ! equal to this term (pressure term 1), such that:
    !
    ! L(w'^3) = - (C_8/tau_w3t) * ( C_8b * (w'^3)^5 / (w'^2)^6 + w'^3 ).
    !
    ! A Taylor Series expansion (truncated after the first derivative term) of
    ! L(w'^3) around w'^3 = w'^3(t) is used to linearize pressure term 1.
    ! Evaluating L(w'^3) at w'^3(t+1):
    !
    ! L( w'^3(t+1) ) = L( w'^3(t) )
    !                  + ( d L(w'^3) / d w'^3 )|_(w'^3=w'^3(t))
    !                    * ( w'^3(t+1) - w'^3(t) ).
    !
    ! After evaluating the expression above, the term has become linearized.  It
    ! is broken down into implicit (LHS) and explicit (RHS) components.
    ! The implicit portion is:
    !
    ! - (C_8/tau_w3t) * ( 5 * C_8b * Sk_wt^4 + 1 ) * w'^3(t+1).
    !
    ! Note:  When the term is brought over to the left-hand side, the sign 
    !        is reversed and the leading "-" in front of the term is changed 
    !        to a "+".
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is 
    ! being advanced to in solving the d(w'^3)/dt equation.
    !
    ! The values of w'^3 are found on the thermodynamic levels, as are the 
    ! values of tau_w3t and Sk_wt (in Sk_wt, w'^3 is found on thermodynamic 
    ! levels and w'^2 is interpolated to thermodynamic levels).

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C8,      & ! Model parameter C_8                        [-]
      C8b,     & ! Model parameter C_8b                       [-]
      tauw3t,  & ! Time-scale tau at thermodynamic levels (k) [s]
      Skw_zt     ! Skewness of w at thermodynamic levels (k)  [-]

    ! Return Variable
    real( kind = core_rknd ) :: lhs

    ! Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
    lhs & 
    = + ( C8 / tauw3t ) * ( 5.0_core_rknd * C8b * Skw_zt**4 + 1.0_core_rknd )

    return
  end function wp3_term_pr1_lhs

  !=============================================================================
! pure function wp3_terms_ta_tp_rhs( wp3_zm, wp3_zmm1,  &
!                                    wp2, wp2m1,  &
!                                    a1, a1_zt, a1m1,  &
!                                    a3, a3_zt, a3m1,  &
!                                    wp3_on_wp2, wp3_on_wp2_m1, &
!                                    rho_ds_zm, rho_ds_zmm1,  &
!                                    invrs_rho_ds_zt,  &
!                                    const_three_halves,  &
!                                    invrs_dzt )  &
! result( rhs )

    ! Description:
    ! Turbulent advection and turbulent production of wp3:  explicit portion of 
    ! the code.
    !
    ! The d(w'^3)/dt equation contains a turbulent advection term:
    !
    ! - (1/rho_ds) * d( rho_ds * w'^4 )/dz;
    !
    ! and a turbulent production term:
    !
    ! + 3 * ( w'^2 / rho_ds ) * d( rho_ds * w'^2 )/dz.
    !
    ! A substitution is made in order to close the turbulent advection term, 
    ! such that:
    !
    ! w'^4 = coef_sig_sqd_w * (w'^2)^2  +  a_1 * ( (w'^3)^2 / w'^2 );
    !
    ! where both a_1 and coef_sig_sqd_w are variables that are functions of
    ! sigma_sqd_w, such that: 
    !
    ! coef_sig_sqd_w = 3*(sigma_sqd_w)^2 + 6*(1 - sigma_sqd_w)*sigma_sqd_w
    !                  + (1 - sigma_sqd_w)^2; and
    !
    ! a_1 = 1 / (1 - sigma_sqd_w).
    !
    ! Since the turbulent advection and turbulent production terms are being
    ! combined, a further substitution is made, such that:
    !
    ! a_3 = coef_sig_sqd_w - 3;
    !
    ! and thus:
    !
    ! w'^4 = (a_3 + 3) * (w'^2)^2  +  a_1 * ( (w'^3)^2 / w'^2 ).
    !
    ! The turbulent production term is rewritten as:
    !
    ! + 3 * ( w'^2 / rho_ds ) * d[ rho_ds * w'^2 ]/dz
    ! = + (3/rho_ds) * d[ rho_ds * (w'^2)^2 ]/dz - (3/2) * d[ (w'^2)^2 ]/dz.
    !
    ! The turbulent advection and turbulent production terms are combined as:
    !
    ! - (1/rho_ds) * d [ rho_ds * a_3 * (w'^2)^2 ] / dz
    ! - (1/rho_ds) * d [ rho_ds * a_1 * ( (w'^3)^2 / w'^2 ) ] / dz
    ! - (3/2) * d [ (w'^2)^2 ] / dz.
    !
    ! The (w'^2)^2 and (w'^3)^2 terms are both linearized, such that:
    !
    ! ( w'^2(t+1) )^2 = - ( w'^2(t) )^2  +  2 * w'^2(t) * w'^2(t+1);
    ! ( w'^3(t+1) )^2 = - ( w'^3(t) )^2  +  2 * w'^3(t) * w'^3(t+1);
    !
    ! which produces implicit and explicit portions of these terms.  The 
    ! explicit portion of these terms is:
    !
    ! + (1/rho_ds) * d [ rho_ds * a_3 * ( w'^2(t) )^2 ] / dz
    ! + (1/rho_ds) * d [ rho_ds * a_1 * ( w'^3(t) )^2 / w'^2(t) ] / dz
    ! + (3/2) * d [ ( w'^2(t) )^2 ] / dz.
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is 
    ! being advanced to in solving the d(w'^3)/dt and d(w'^2)/dt equations.
    !
    ! The explicit portion of these terms is discretized as follows:
    !
    ! The values of w'^3 are found on the thermodynamic levels, while the values
    ! of w'^2, a_1, and a_3 are found on the momentum levels.  Additionally, the
    ! values of rho_ds_zm are found on the momentum levels, and the values of
    ! invrs_rho_ds_zt are found on the thermodynamic levels.  The variable w'^3
    ! is interpolated to the intermediate momentum levels.  The values of the
    ! mathematical expressions (called F, G, and H here) within the dF/dz,
    ! dG/dz, and dH/dz terms are computed on the momentum levels.  Then, the
    ! derivatives (d/dz) of the expressions (F, G, and H) are taken over the
    ! central thermodynamic level, where dF/dz and dG/dz are multiplied by
    ! invrs_rho_ds_zt, and where dH/dz is multiplied by 3/2.  This yields the
    ! desired results.  In this function, the values of F, G, and H are as
    ! follows:
    !
    ! F = rho_ds_zm * a_3(t) * ( w'^2(t) )^2;
    !
    ! G = rho_ds_zm * a_1(t) * ( w'^3(t) )^2 / w'^2(t); and
    !
    ! H = ( w'^2(t) )^2.
    !
    !
    ! ------------------------------------------------wp3p1-------------- t(k+1)
    !
    ! ===a3====wp2====rho_ds_zm====a1======================wp3(interp)=== m(k)
    !
    ! ---dH/dz---dF/dz----invrs_rho_ds_zt----dG/dz----wp3---------------- t(k)
    !
    ! ===a3m1==wp2m1==rho_ds_zmm1==a1m1====================wp3(interp)=== m(k-1)
    !
    ! ------------------------------------------------wp3m1-------------- t(k-1)
    !
    ! The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1) correspond 
    ! with altitudes zt(k+1), zm(k), zt(k), zm(k-1), and zt(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is used
    ! for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) )

    ! References:
    !-----------------------------------------------------------------------

!   use constants_clubb, only:  &
!       w_tol_sqd

!   use model_flags, only:  &
!       l_standard_term_ta

!   implicit none

    ! Input Variables
!   real, intent(in) ::  & 
!     wp3_zm,             & ! w'^3 interpolated to momentum lev. (k)   [m^3/s^3]
!     wp3_zmm1,           & ! w'^3 interpolated to momentum lev. (k-1) [m^3/s^3]
!     wp2,                & ! w'^2(k)                                  [m^2/s^2]
!     wp2m1,              & ! w'^2(k-1)                                [m^2/s^2]
!     a1,                 & ! a_1(k)                                   [-]
!     a1_zt,              & ! a_1 interpolated to thermo. level (k)    [-]
!     a1m1,               & ! a_1(k-1)                                 [-]
!     a3,                 & ! a_3(k)                                   [-]
!     a3_zt,              & ! a_3 interpolated to thermo. level (k)    [-]
!     a3m1,               & ! a_3(k-1)                                 [-]
!     wp3_on_wp2,         & ! (k) [m/s]
!     wp3_on_wp2_m1,      & ! (k-1)                  [m/s]
!     rho_ds_zm,          & ! Dry, static density at moment. lev (k)   [kg/m^3]
!     rho_ds_zmm1,        & ! Dry, static density at moment. lev (k-1) [kg/m^3]
!     invrs_rho_ds_zt,    & ! Inv dry, static density @ thermo lev (k) [m^3/kg]
!     const_three_halves, & ! "3/2" ("0" is sent in for wp3_ta budget) [-]
!     invrs_dzt             ! Inverse of grid spacing (k)              [1/m]

    ! Return Variable
!   real :: rhs


!   if ( l_standard_term_ta ) then

       ! The turbulent advection term is discretized normally, in accordance
       ! with the model equations found in the documentation and the description
       ! listed above.

!      rhs & 
!      = + invrs_rho_ds_zt &
!          * invrs_dzt &
!            * (   rho_ds_zm * a3 * wp2**2 &
!                - rho_ds_zmm1 * a3m1 * wp2m1**2 &
!              ) &
!        + invrs_rho_ds_zt &
!          * invrs_dzt &
!            * (   rho_ds_zm * a1 &
!                  * wp3_zm * wp3_on_wp2 &
!                - rho_ds_zmm1 * a1m1 &
!                  * wp3_zmm1 * wp3_on_wp2_m1 &
!              ) &
!        + const_three_halves &
!          * invrs_dzt * ( wp2**2 - wp2m1**2 )

!   else

       ! Brian tried a new discretization for the turbulent advection term, 
       ! which contains the term:
       !  - (1/rho_ds) * d [ rho_ds * a_1 * (w'^3)^2 / w'^2 ] / dz.  In order
       ! to help stabilize w'^3, a_1 has been pulled outside of the derivative.
       ! This effects the right-hand side of the equation, as well as the 
       ! left-hand side.

       ! Additionally, the discretization of the turbulent advection term, which
       ! contains the term:
       !  - (1/rho_ds) * d [ rho_ds * (a_3 + 3) * (w'^2)^2 ] / dz, has been 
       ! altered to pull (a_3 + 3) outside of the derivative.  This was done in
       ! order to help stabilize w'^3.  This effects the right-hand side of the 
       ! equation, as well as the left-hand side.

!      rhs & 
!      = + invrs_rho_ds_zt &
!          * a3_zt * invrs_dzt &
!            * (   rho_ds_zm * wp2**2 &
!                - rho_ds_zmm1 * wp2m1**2 ) &
!        + invrs_rho_ds_zt &
!          * a1_zt * invrs_dzt & 
!            * (   rho_ds_zm &
!                  * ( wp3_zm * wp3_on_wp2 ) & 
!                - rho_ds_zmm1 &
!                  * ( wp3_zmm1 * wp3_on_wp2_m1 ) & 
!              ) &
!        + const_three_halves &
!          * invrs_dzt * ( wp2**2 - wp2m1**2 )

       ! End of code that pulls out a3.
       ! End of Brian's a1 change.  Feb. 14, 2008.

!   endif ! l_standard_term_ta


!   return
! end function wp3_terms_ta_tp_rhs

  !=============================================================================
  pure function wp3_terms_bp1_pr2_rhs( C11_Skw_fnc, thv_ds_zt, wp2thvp ) & 
  result( rhs )

    ! Description:
    ! Buoyancy production of w'^3 and w'^3 pressure term 2:  explicit portion of
    ! the code.
    !
    ! The d(w'^3)/dt equation contains a buoyancy production term:
    !
    ! + 3 (g/thv_ds) w'^2th_v';
    !
    ! and pressure term 2:
    !
    ! - C_11 ( -3 w'^3 dw/dz + 3 (g/thv_ds) w'^2th_v' ).
    !
    ! The w'^3 buoyancy production term is completely explicit, while w'^3 
    ! pressure term 2 has both implicit and explicit components.  The buoyancy 
    ! production term and the explicit portion of pressure term 2 are combined 
    ! and solved together as:
    !
    ! + ( 1 - C_ll ) ( 3 (g/thv_ds) w'^2th_v' ).

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use constants_clubb, only: & ! Constant(s) 
        grav ! Gravitational acceleration [m/s^2]

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C11_Skw_fnc, & ! C_11 parameter with Sk_w applied (k)        [-]
      thv_ds_zt,   & ! Dry, base-state theta_v at thermo. lev. (k) [K]
      wp2thvp        ! w'^2th_v'(k)                                [K m^2/s^2]

    ! Return Variable
    real( kind = core_rknd ) :: rhs

    rhs & 
    = + ( 1.0_core_rknd - C11_Skw_fnc ) * 3.0_core_rknd * ( grav / thv_ds_zt ) * wp2thvp

    return
  end function wp3_terms_bp1_pr2_rhs

  !=============================================================================
  pure function wp3_term_bp2_rhs( C15, Kh_zt, wpthvp, wpthvp_m1, &
                                  dum_dz, dum_dz_m1, dvm_dz, dvm_dz_m1, &
                                  upwp, upwp_m1, vpwp, vpwp_m1, &
                                  thv_ds_zt, invrs_dzt ) & 
    result( rhs )

    ! Description:
    !   Experimental term from CLUBB TRAC ticket #411. The derivative here is of
    !   the form:
    !   - C_15 * Kh * ∂{ grav / thv_ds * [w'th_v'(k) - w'th_v'(k-1)] 
    !                    -[ u'w'(k) * ∂u(k)/∂z - u'w'(k-1) * ∂u(k-1)/∂z ]
    !                    -[ v'w'(k) * ∂v(k)/∂z - v'w'(k-1) * ∂v(k-1)/∂z ]  }/∂z.
    !
    !   This does not appear in Andre et al. 1976 or Bougeault et al. 1981, but
    !   is based on experiments in matching LES data.
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use constants_clubb, only: & ! Constant(s) 
        grav ! Gravitational acceleration [m/s^2]

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      C15,       & ! Model parameter C15                [-]
      Kh_zt,     & ! Eddy-diffusivity on moment. levels [m^2/s]
      wpthvp,    & ! w'th_v'(k)                         [K m/s]
      wpthvp_m1, & ! w'th_v'(k-1)                       [K m/s]
      dum_dz,    & ! d u wind dz (k)                    [m/s]
      dvm_dz,    & ! d v wind dz (k)                    [m/s]
      dum_dz_m1, & ! d u wind dz (k-1)                  [m/s]
      dvm_dz_m1, & ! d v wind dz (k-1)                  [m/s]
      upwp,      & ! u'v'(k)                            [m^2/s^2]
      upwp_m1,   & ! u'v'(k-1)				[m^2/s^2]
      vpwp,      & ! v'w'(k)                            [m^2/s^2]
      vpwp_m1,   & ! v'w'(k-1)                          [m^2/s^2]
      thv_ds_zt, & ! Dry, base-state theta_v at thermo. lev. (k) [K]
      invrs_dzt    ! Inverse of grid spacing (k)        [1/m]

    ! Return Variable
    real( kind = core_rknd ) :: rhs

    ! ---- Begin Code ----

!   rhs =  - C15 * Kh_zt * invrs_dzt * grav / thv_ds_zt * ( wpthvp - wpthvp_m1 )

    rhs =  - C15 * Kh_zt * invrs_dzt * &
      ( grav / thv_ds_zt * ( wpthvp - wpthvp_m1 ) &
      - ( upwp * dum_dz - upwp_m1 * dum_dz_m1 ) &
      - ( vpwp * dvm_dz - vpwp_m1 * dvm_dz_m1 ) )

    return
  end function wp3_term_bp2_rhs


  !=============================================================================
  pure function wp3_term_pr1_rhs( C8, C8b, tauw3t, Skw_zt, wp3 ) & 
  result( rhs )

    ! Description:
    ! Pressure term 1 for w'^3:  explicit portion of the code.
    !
    ! Pressure term 1 is the term:
    !
    ! - (C_8/tau_w3t) * ( C_8b * Sk_wt^4 + 1 ) * w'^3;
    !
    ! where Sk_wt = w'^3 / (w'^2)^(3/2).
    !
    ! This term needs to be linearized, so function L(w'^3) is defined to be 
    ! equal to this term (pressure term 1), such that:
    !
    ! L(w'^3) = - (C_8/tau_w3t) * ( C_8b * (w'^3)^5 / (w'^2)^6 + w'^3 ).
    !
    ! A Taylor Series expansion (truncated after the first derivative term) of
    ! L(w'^3) around w'^3 = w'^3(t) is used to linearize pressure term 1.
    ! Evaluating L(w'^3) at w'^3(t+1):
    !
    ! L( w'^3(t+1) ) = L( w'^3(t) )
    !                  + ( d L(w'^3) / d w'^3 )|_(w'^3=w'^3(t))
    !                    * ( w'^3(t+1) - w'^3(t) ).
    !
    ! After evaluating the expression above, the term has become linearized.  It
    ! is broken down into implicit (LHS) and explicit (RHS) components.
    ! The explicit portion is:
    !
    ! + (C_8/tau_w3t) * ( 4 * C_8b * Sk_wt^4 + 1 ) * w'^3(t).
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is 
    ! being advanced to in solving the d(w'^3)/dt equation.
    !
    ! The values of w'^3 are found on the thermodynamic levels, as are the 
    ! values of tau_w3t and Sk_wt (in Sk_wt, w'^3 is found on thermodynamic 
    ! levels and w'^2 is interpolated to thermodynamic levels).

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C8,      & ! Model parameter C_8                        [-]
      C8b,     & ! Model parameter C_8b                       [-]
      tauw3t,  & ! Time-scale tau at thermodynamic levels (k) [s]
      Skw_zt,  & ! Skewness of w at thermodynamic levels (k)  [-]
      wp3        ! w'^3(k)                                    [m^3/s^3]

    ! Return Variable
    real( kind = core_rknd ) :: rhs

    rhs & 
    = + ( C8 / tauw3t ) * ( 4.0_core_rknd * C8b * Skw_zt**4 ) * wp3

    return
  end function wp3_term_pr1_rhs

!===============================================================================

end module advance_wp2_wp3_module
