!------------------------------------------------------------------------
! $Id$
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
             wp3_term_ta_new_pdf_lhs, &
             wp3_term_ta_ADG1_lhs, & 
             wp3_term_tp_lhs, & 
             wp3_terms_ac_pr2_lhs, & 
             wp3_term_pr1_lhs, & 
             wp3_term_ta_explicit_rhs, &
             wp3_terms_bp1_pr2_rhs, & 
             wp3_term_pr1_rhs, &
             wp3_term_bp2_rhs, &
             wp2_term_ta_lhs_all, & 
             wp2_terms_ac_pr2_lhs_all, & 
             wp2_term_dp1_lhs_all, & 
             wp2_term_pr1_lhs_all, & 
             wp2_terms_bp_pr2_rhs_all, & 
             wp2_term_dp1_rhs_all, &
             wp2_term_pr3_rhs_all, & 
             wp2_term_pr1_rhs_all, & 
             wp3_term_ta_new_pdf_lhs_all, &
             wp3_term_ta_ADG1_lhs_all, & 
             wp3_term_tp_lhs_all, & 
             wp3_terms_ac_pr2_lhs_all, & 
             wp3_term_pr1_lhs_all, & 
             wp3_term_ta_explicit_rhs_all, &
             wp3_terms_bp1_pr2_rhs_all, & 
             wp3_term_pr1_rhs_all, &
             wp3_term_bp2_rhs_all

  ! Private named constants to avoid string comparisons
  integer, parameter, private :: &
    clip_wp2 = 12 ! Named constant for wp2 clipping.
                  ! NOTE: This must be the same as the clip_wp2 declared in
                  ! clip_explicit!

  contains

  !=============================================================================
  subroutine advance_wp2_wp3( dt, sfc_elevation, sigma_sqd_w, wm_zm,   & ! In
                              wm_zt, a3, a3_zt, wp3_on_wp2, wp4,       & ! In
                              wpthvp, wp2thvp, um, vm, upwp, vpwp,     & ! In
                              up2, vp2, Kh_zm, Kh_zt, tau_zm, tau_zt,  & ! In
                              tau_C1_zm, Skw_zm, Skw_zt, rho_ds_zm,    & ! In
                              rho_ds_zt, invrs_rho_ds_zm,              & ! In
                              invrs_rho_ds_zt, radf, thv_ds_zm,        & ! In
                              thv_ds_zt, mixt_frac, Cx_fnc_Richardson, & ! In
                              wp2_splat, wp3_splat,                    & ! intent(in)
                              pdf_implicit_coefs_terms,                & ! In
                              wprtp, wpthlp, rtp2, thlp2,              & ! In
                              wp2, wp3, wp3_zm, wp2_zt )                 ! Inout

    ! Description:
    ! Advance w'^2 and w'^3 one timestep.

    ! References:
    ! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:wp2_wp3_eqns
    !
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

    use sponge_layer_damping, only: &
        wp2_sponge_damp_settings, & ! Variable(s)
        wp3_sponge_damp_settings, &
        wp2_sponge_damp_profile,  &
        wp3_sponge_damp_profile,  &
        sponge_damp_xp2, & ! Procedure(s)
        sponge_damp_xp3

    use stats_type_utilities, only: &
        stat_begin_update, & ! Procedure(s)
        stat_end_update, &
        stat_update_var

    use stats_variables, only: &
        iC1_Skw_fnc, &  ! Variable(s)
        iC11_Skw_fnc, &
        iwp2_sdmp, &
        iwp3_sdmp, &
        stats_zm, &
        stats_zt, &
        l_stats_samp

    use constants_clubb, only:  & 
        fstderr,   & ! Variables
        one,       &
        one_half,  &
        one_third, &
        w_tol_sqd, &
        eps

    use pdf_parameter_module, only: &
        implicit_coefs_terms    ! Variable Type

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_fatal_error              ! Constant

    use model_flags, only: &
        l_damp_wp2_using_em, &  ! Logical(s)
        l_use_C11_Richardson

    implicit none

    intrinsic :: exp

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      dt                 ! Model timestep                            [s]

    real( kind = core_rknd ), intent(in) ::  &
      sfc_elevation      ! Elevation of ground level                 [m AMSL]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  & 
      sigma_sqd_w,       & ! sigma_sqd_w (momentum levels)             [-]
      wm_zm,             & ! w wind component on momentum levels       [m/s]
      wm_zt,             & ! w wind component on thermodynamic levels  [m/s]
      a3,                & ! a_3 (momentum levels); See eqn. 25 in `Equations for CLUBB' [-]
      a3_zt,             & ! a_3 interpolated to thermodynamic levels  [-]
      wp3_on_wp2,        & ! Smoothed version of wp3 / wp2             [m/s]
      wp4,               & ! w'^4 (momentum levels)                    [m^4/s^4]
      wpthvp,            & ! w'th_v' (momentum levels)                 [K m/s]
      wp2thvp,           & ! w'^2th_v' (thermodynamic levels)          [K m^2/s^2]
      um,                & ! u wind component (thermodynamic levels)   [m/s]
      vm,                & ! v wind component (thermodynamic levels)   [m/s]
      upwp,              & ! u'w' (momentum levels)                    [m^2/s^2]
      vpwp,              & ! v'w' (momentum levels)                    [m^2/s^2]
      up2,               & ! u'^2 (momentum levels)                    [m^2/s^2]
      vp2,               & ! v'^2 (momentum levels)                    [m^2/s^2]
      Kh_zm,             & ! Eddy diffusivity on momentum levels       [m^2/s]
      Kh_zt,             & ! Eddy diffusivity on thermodynamic levels  [m^2/s]
      tau_zm,            & ! Time-scale tau on momentum levels         [s]
      tau_zt,            & ! Time-scale tau on thermodynamic levels    [s]
      tau_C1_zm,         & ! Tau values used for the C1 (dp1) term in wp2 [s]
      Skw_zm,            & ! Skewness of w on momentum levels          [-]
      Skw_zt,            & ! Skewness of w on thermodynamic levels     [-]
      rho_ds_zm,         & ! Dry, static density on momentum levels    [kg/m^3]
      rho_ds_zt,         & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zm,   & ! Inv. dry, static density @ momentum levs. [m^3/kg]
      invrs_rho_ds_zt,   & ! Inv. dry, static density @ thermo. levs.  [m^3/kg]
      radf,              & ! Buoyancy production at the CL top         [m^2/s^3]
      thv_ds_zm,         & ! Dry, base-state theta_v on momentum levs. [K]
      thv_ds_zt,         & ! Dry, base-state theta_v on thermo. levs.  [K]
      mixt_frac,         & ! Weight of 1st normal distribution         [-]
      wprtp,             & ! Flux of total water mixing ratio          [m/s kg/kg]
      wpthlp,            & ! Flux of liquid water potential temp.      [m/s K]
      rtp2,              & ! Variance of rt (overall)                  [kg^2/kg^2]
      thlp2,             & ! Variance of thl (overall)                 [K^2]
      Cx_fnc_Richardson, & ! Cx_fnc from Richardson_num                [-]
      wp2_splat,         & ! Tendency of <w'2> due to vertical compression of eddies [m^2/s^3]
      wp3_splat            ! Tendency of <w'3> due to vertical compression of eddies [m^3/s^4]

    type(implicit_coefs_terms), dimension(gr%nz), intent(in) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    ! Input/Output
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) ::  & 
      wp2,  & ! w'^2 (momentum levels)                    [m^2/s^2]
      wp3,  & ! w'^3 (thermodynamic levels)               [m^3/s^3]
      wp3_zm  ! w'^3 interpolated to momentum levels      [m^3/s^3]

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) ::  &
      wp2_zt  ! w'^2 interpolated to thermodyamic levels  [m^2/s^2]

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
      C11_Skw_fnc, & ! C_11 parameter with Sk_w applied             [-]
    ! End Vince Larson's addition.
      C16_fnc        ! C_16 parameter                               [-]

    integer :: k ! Array indices

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
!          tauw3t(k) = tau_zt(k) / ( 0.005_core_rknd*Skw**4 + one )
!
!        end do

    tauw3t = tau_zt

    ! Vince Larson added code to make C11 function of Skw. 13 Mar 2005
    ! If this code is used, C11 is no longer relevant, i.e. constants
    !    are hardwired.

    if ( l_use_C11_Richardson ) then
      C11_Skw_fnc = Cx_fnc_Richardson
    else
      ! Calculate C_{1} and C_{11} as functions of skewness of w.
      ! The if..then here is only for computational efficiency -dschanen 2 Sept 08
      if ( abs(C11-C11b) > abs(C11+C11b)*eps/2 ) then
        C11_Skw_fnc(1:gr%nz) =  & 
          C11b + (C11-C11b)*EXP( -one_half * (Skw_zt(1:gr%nz)/C11c)**2 )
      else
        C11_Skw_fnc(1:gr%nz) = C11b
      end if
    end if ! l_use_C11_Richardson

    ! The if..then here is only for computational efficiency -dschanen 2 Sept 08
    if ( abs(C1-C1b) > abs(C1+C1b)*eps/2 ) then
      C1_Skw_fnc(1:gr%nz) =  & 
        C1b + (C1-C1b)*EXP( -one_half * (Skw_zm(1:gr%nz)/C1c)**2 )
    else
      C1_Skw_fnc(1:gr%nz) = C1b 
    end if

    if ( l_damp_wp2_using_em ) then
      ! Insert 1/3 here to account for the fact that in the dissipation term, 
      !   (2/3)*em = (2/3)*(1/2)*(wp2+up2+vp2).  Then we can insert wp2, up2,
      !   and vp2 directly into the dissipation subroutines without prefixing them by (1/3).
      C1_Skw_fnc(1:gr%nz) = one_third * C1_Skw_fnc(1:gr%nz)
    end if
 
    !C11_Skw_fnc = C11
    !C1_Skw_fnc = C1

    ! Set C16_fnc based on Richardson_num
    C16_fnc = Cx_fnc_Richardson

    if ( clubb_at_least_debug_level( 0 ) ) then
      ! Assertion check for C11_Skw_fnc
      if ( any( C11_Skw_fnc(:) > one ) .or. any( C11_Skw_fnc(:) < 0._core_rknd ) ) then
        write(fstderr,*) "The C11_Skw_fnc is outside the valid range for this variable"
        err_code = clubb_fatal_error
        return
      end if

      if ( any( C16_fnc(:) > one ) .or. any( C16_fnc(:) < 0._core_rknd ) ) then
        write(fstderr,*) "The C16_fnc is outside the valid range for this variable"
        err_code = clubb_fatal_error
        return
      end if
    end if

    if ( l_stats_samp ) then
      call stat_update_var( iC11_Skw_fnc, C11_Skw_fnc, stats_zt )
      call stat_update_var( iC1_Skw_fnc, C1_Skw_fnc, stats_zm )
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

    ! Solve semi-implicitly
    call wp23_solve( dt, sfc_elevation, sigma_sqd_w, wm_zm, & ! Intent(in)
                     wm_zt, a3, a3_zt, wp3_on_wp2, wp4,     & ! Intent(in)
                     wpthvp, wp2thvp, um, vm, upwp, vpwp,   & ! Intent(in)
                     up2, vp2, Kw1, Kw8, Kh_zt, Skw_zt,     & ! Intent(in)
                     tau_zm, tauw3t, tau_C1_zm, C1_Skw_fnc, & ! Intent(in)
                     C11_Skw_fnc, C16_fnc, rho_ds_zm,       & ! Intent(in)
                     rho_ds_zt, invrs_rho_ds_zm,            & ! Intent(in)
                     invrs_rho_ds_zt, radf,                 & ! Intent(in)
                     thv_ds_zm, thv_ds_zt,                  & ! Intent(in)
                     wp2_splat, wp3_splat,                  & ! Intent(in)
                     pdf_implicit_coefs_terms,              & ! Intent(in)
                     wprtp, wpthlp, rtp2, thlp2,            & ! Intent(in)
                     wp2, wp3, wp3_zm, wp2_zt )               ! Intent(inout)

    ! When selected, apply sponge damping after wp2 and wp3 have been advanced.
    if ( wp2_sponge_damp_settings%l_sponge_damping ) then

       if ( l_stats_samp ) then
          call stat_begin_update( iwp2_sdmp, wp2 / dt, stats_zm )
       endif

       wp2 = sponge_damp_xp2( dt, gr%zm, wp2, w_tol_sqd, &
                              wp2_sponge_damp_profile )

       if ( l_stats_samp ) then
          call stat_end_update( iwp2_sdmp, wp2 / dt, stats_zm )
       endif

    endif ! wp2_sponge_damp_settings%l_sponge_damping

    if ( wp3_sponge_damp_settings%l_sponge_damping ) then

       if ( l_stats_samp ) then
          call stat_begin_update( iwp3_sdmp, wp3 / dt, stats_zt )
       endif

       wp3 = sponge_damp_xp3( dt, gr%zt, wp3, wp3_sponge_damp_profile )

       if ( l_stats_samp ) then
          call stat_end_update( iwp3_sdmp, wp3 / dt, stats_zt )
       endif

    endif ! wp3_sponge_damp_settings%l_sponge_damping

    if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code == clubb_fatal_error ) then  

            write(fstderr,*) "Error in advance_wp2_wp3"

            write(fstderr,*) "Intent(in)"

            write(fstderr,*) "gr%zt = ", gr%zt, new_line('c')
            write(fstderr,*) "dt = ", dt, new_line('c')
            write(fstderr,*) "sfc_elevation = ", sfc_elevation, new_line('c')
            write(fstderr,*) "sigma_sqd_w = ", sigma_sqd_w, new_line('c')
            write(fstderr,*) "wm_zm = ", wm_zm, new_line('c')
            write(fstderr,*) "wm_zt = ", wm_zt, new_line('c')
            write(fstderr,*) "wp4 = ", wp4, new_line('c')
            write(fstderr,*) "wpthvp = ", wpthvp, new_line('c')
            write(fstderr,*) "wp2thvp = ", wp2thvp, new_line('c')
            write(fstderr,*) "um = ", um, new_line('c')
            write(fstderr,*) "vm = ", vm, new_line('c')
            write(fstderr,*) "upwp = ", upwp, new_line('c')
            write(fstderr,*) "vpwp = ", vpwp, new_line('c')
            write(fstderr,*) "up2 = ", up2, new_line('c')
            write(fstderr,*) "vp2 = ", vp2, new_line('c')
            write(fstderr,*) "Kh_zm = ", Kh_zm, new_line('c')
            write(fstderr,*) "Kh_zt = ", Kh_zt, new_line('c')
            write(fstderr,*) "tau_zm = ", tau_zm, new_line('c')
            write(fstderr,*) "tau_zt = ", tau_zt, new_line('c')
            write(fstderr,*) "Skw_zm = ", Skw_zm, new_line('c')
            write(fstderr,*) "Skw_zt = ", Skw_zt, new_line('c')
            write(fstderr,*) "mixt_frac = ", mixt_frac, new_line('c')
            write(fstderr,*) "a3 = ", a3, new_line('c')
            write(fstderr,*) "a3_zt = ", a3_zt, new_line('c')
            write(fstderr,*) "wp3_on_wp2 = ", wp3_on_wp2, new_line('c')
            write(fstderr,*) "tau_C1_zm = ", tau_C1_zm, new_line('c')
            write(fstderr,*) "rho_ds_zm = ", rho_ds_zm, new_line('c')
            write(fstderr,*) "rho_ds_zt = ", rho_ds_zt, new_line('c')
            write(fstderr,*) "invrs_rho_ds_zm = ", invrs_rho_ds_zm, new_line('c')
            write(fstderr,*) "invrs_rho_ds_zt = ", invrs_rho_ds_zt, new_line('c')
            write(fstderr,*) "radf = ", radf, new_line('c')
            write(fstderr,*) "thv_ds_zm = ", thv_ds_zm, new_line('c')
            write(fstderr,*) "thv_ds_zt = ", thv_ds_zt, new_line('c')
            write(fstderr,*) "Cx_fnc_Richardson = ", Cx_fnc_Richardson, new_line('c')
            write(fstderr,*) "pdf_implicit_coefs_terms = ", pdf_implicit_coefs_terms
            write(fstderr,*) new_line('c')

            write(fstderr,*) "Intent(in/out)"

            write(fstderr,*) "wp2_zt = ", wp2_zt, new_line('c')
            write(fstderr,*) "wp3_zm = ", wp3_zm, new_line('c')
            write(fstderr,*) "wp2 = ", wp2, new_line('c')
            write(fstderr,*) "wp3 = ", wp3, new_line('c')

        end if ! fatal error
    end if

    return

  end subroutine advance_wp2_wp3

  !=============================================================================
  subroutine wp23_solve( dt, sfc_elevation, sigma_sqd_w, wm_zm, & ! Intent(in)
                         wm_zt, a3, a3_zt, wp3_on_wp2, wp4,     & ! Intent(in)
                         wpthvp, wp2thvp, um, vm, upwp, vpwp,   & ! Intent(in)
                         up2, vp2, Kw1, Kw8, Kh_zt, Skw_zt,     & ! Intent(in)
                         tau1m, tauw3t, tau_C1_zm, C1_Skw_fnc,  & ! Intent(in)
                         C11_Skw_fnc, C16_fnc, rho_ds_zm,       & ! Intent(in)
                         rho_ds_zt, invrs_rho_ds_zm,            & ! Intent(in)
                         invrs_rho_ds_zt, radf,                 & ! Intent(in)
                         thv_ds_zm, thv_ds_zt,                  & ! Intent(in)
                         wp2_splat, wp3_splat,                  & ! Intent(in)
                         pdf_implicit_coefs_terms,              & ! Intent(in)
                         wprtp, wpthlp, rtp2, thlp2,            & ! Intent(in)
                         wp2, wp3, wp3_zm, wp2_zt )               ! Intent(inout)

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
        w_tol_sqd,                & ! Variables(s)
        max_mag_correlation_flux, &
        one,                      &
        zero,                     &
        zero_threshold,           &
        fstderr

    use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_fatal_error              ! Constants

    use model_flags, only:  & 
        l_tke_aniso,                  & ! Variable(s)
        l_hole_fill,                  &
        l_explicit_turbulent_adv_wp3, &
        l_gmres,                      &
        l_min_wp2_from_corr_wx

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use lapack_wrap, only:  & 
        band_solve,  & ! Procedure(s) 
        band_solvex

    use fill_holes, only: & 
        fill_holes_vertical

    use clip_explicit, only: &
        clip_variance, & ! Procedure(s)
        clip_variance_level, &
        clip_skewness

    use pdf_closure_module, only: &
        iiPDF_ADG1, & ! Variable(s)
        iiPDF_new,  &
        iiPDF_type

    use pdf_parameter_module, only: &
        implicit_coefs_terms    ! Variable Type

    use stats_type_utilities, only: & 
        stat_begin_update, & ! Procedure(s)
        stat_update_var, &
        stat_update_var_pt, &
        stat_end_update, &
        stat_end_update_pt

    use stats_variables, only:  & 
        stats_zm, & ! Variable(s)
        stats_zt, & 
        stats_sfc, & 
        l_stats_samp, &
        icoef_wp4_implicit, &
        iwp2_ta, & 
        iwp2_ma, & 
        iwp2_pd, & 
        iwp2_ac, & 
        iwp2_dp1, & 
        iwp2_dp2, & 
        iwp2_pr1, & 
        iwp2_pr2, &
        iwp3_ta, & 
        iwp3_ma, & 
        iwp3_tp, & 
        iwp3_ac, & 
        iwp3_dp1, & 
        iwp3_pr1, &
        iwp3_pr2, &
        iwp3_pr3, &
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
        ztscr16

    implicit none

    ! External
    intrinsic :: max, min, sqrt

    ! Parameter Constants
    integer, parameter :: & 
      nrhs = 1      ! Number of RHS vectors

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
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
      wp4,             & ! w'^4 (momentum levels)                    [m^4/s^4]
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
      tau_C1_zm,       & ! Tau values used for the C1 (dp1) term in wp2 [s]
      C1_Skw_fnc,      & ! C_1 parameter with Sk_w applied           [-]
      C11_Skw_fnc,     & ! C_11 parameter with Sk_w applied          [-]
      C16_fnc,         & ! C_16 parameter                            [-]
      rho_ds_zm,       & ! Dry, static density on momentum levels    [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ momentum levs. [m^3/kg]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. levs.  [m^3/kg]
      radf,            & ! Buoyancy production at CL top             [m^2/s^3]
      thv_ds_zm,       & ! Dry, base-state theta_v on momentum levs. [K]
      thv_ds_zt,       & ! Dry, base-state theta_v on thermo. levs.  [K]
      wprtp,           & ! Flux of total water mixing ratio          [m/s kg/kg]
      wpthlp,          & ! Flux of liquid water potential temp.      [m/s K]
      rtp2,            & ! Variance of rt (overall)                  [kg^2/kg^2]
      thlp2,           & ! Variance of thl (overall)                 [K^2]
      wp2_splat,       & ! Tendency of <w'2> due to vertical compression of eddies  [m^2/s^3]
      wp3_splat          ! Tendency of <w'3> due to vertical compression of eddies  [m^3/s^4]

    type(implicit_coefs_terms), dimension(gr%nz), intent(in) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) ::  & 
      wp2,  & ! w'^2 (momentum levels)                            [m^2/s^2]
      wp3,  & ! w'^3 (thermodynamic levels)                       [m^3/s^3]
      wp3_zm  ! w'^3 interpolated to momentum levels      [m^3/s^3]

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) ::  &
      wp2_zt  ! w'^2 interpolated to thermodyamic levels          [m^2/s^2]

    ! Local Variables
    real( kind = core_rknd ), dimension(5,2*gr%nz) ::  & 
      lhs ! Implicit contributions to wp2/wp3 (band diag. matrix)

    real( kind = core_rknd ), dimension(2*gr%nz) ::  & 
      rhs,      & ! RHS of band matrix
      rhs_save    ! Saved RHS of band matrix

!        real, target, dimension(2*gr%nz) ::
    real( kind = core_rknd ), dimension(2*gr%nz) ::  & 
      solut ! Solution to band diagonal system.

    real( kind = core_rknd ), dimension(gr%nz) :: & 
      coef_wp4_implicit_zt, & ! <w'^4>|_zt=coef_wp4_implicit_zt*<w'^2>|_zt^2 [-]
      coef_wp4_implicit       ! <w'^4> = coef_wp4_implicit * <w'^2>^2        [-]

    real( kind = core_rknd ), dimension(gr%nz) ::  & 
      a1,   & ! a_1 (momentum levels); See eqn. 23 in `Equations for CLUBB' [-]
      a1_zt   ! a_1 interpolated to thermodynamic levels                    [-]

!      real, dimension(gr%nz) ::  &
!        wp2_n ! w'^2 at the previous timestep           [m^2/s^2]

    real( kind = core_rknd ) ::  & 
      rcond  ! Est. of the reciprocal of the condition #

    real( kind = core_rknd ), dimension(gr%nz,5) :: &
      wp3_pr3_lhs ! wp3_pr3 (implicit) contribution to lhs

    real( kind = core_rknd ) :: &
      threshold    ! Minimum value for wp2    [m^2/s^2]

    ! Array indices
    integer :: k, km1, kp1, k_wp2, k_wp3

    ! Set logical to true for Crank-Nicholson diffusion scheme
    ! or to false for completely implicit diffusion scheme.
    ! Note:  Although Crank-Nicholson diffusion has usually been used for wp2
    !        and wp3 in the past, we found that using completely implicit
    !        diffusion stabilized the deep convective cases more while having
    !        almost no effect on the boundary layer cases.  Brian; 1/4/2008.
!    logical, parameter :: l_crank_nich_diff = .true.
    logical, parameter :: l_crank_nich_diff = .false.

  !-----------------------------------------------------------------------
    !----- Begin Code -----

    if ( .not. l_explicit_turbulent_adv_wp3 ) then

       if ( iiPDF_type == iiPDF_new ) then

          ! Unpack coef_wp4_implicit from pdf_implicit_coefs_terms.
          ! Since PDF parameters and the resulting implicit coefficients and
          ! explicit terms are calculated on thermodynamic levels, the <w'^4>
          ! implicit coefficient needs to be unpacked as coef_wp4_implicit_zt.
          coef_wp4_implicit_zt = pdf_implicit_coefs_terms%coef_wp4_implicit

          ! The values of <w'^4> are located on momentum levels.  Interpolate
          ! coef_wp4_implicit_zt to momentum levels as coef_wp4_implicit.  The
          ! discretization diagram is found in the description section of
          ! function wp3_term_ta_new_pdf_lhs below.  These values are always
          ! positive.
          coef_wp4_implicit = max( zt2zm( coef_wp4_implicit_zt ), &
                                   zero_threshold )

          ! Set the value of coef_wp4_implicit to 0 at the lower boundary and at
          ! the upper boundary.  This sets the value of <w'^4> to 0 at the lower
          ! and upper boundaries.
          coef_wp4_implicit(1) = zero
          coef_wp4_implicit(gr%nz) = zero

          if ( l_stats_samp ) then
             call stat_update_var( icoef_wp4_implicit, coef_wp4_implicit, &
                                   stats_zm )
          endif ! l_stats_samp

       elseif ( iiPDF_type == iiPDF_ADG1 ) then

          ! Define a_1 and a_3 (both are located on momentum levels).
          ! They are variables that are both functions of sigma_sqd_w (where
          ! sigma_sqd_w is located on momentum levels).
          a1 = one / ( one - sigma_sqd_w )

          ! Interpolate a_1 from momentum levels to thermodynamic levels.  This
          ! will be used for the w'^3 turbulent advection (ta) term.
          a1_zt = max( zm2zt( a1 ), zero_threshold )  ! Positive def. quantity

       endif ! iiPDF_type

    endif ! .not. l_explicit_turbulent_adv_wp3

    ! Compute the explicit portion of the w'^2 and w'^3 equations.
    ! Build the right-hand side vector.
    call wp23_rhs( dt, wp2, wp3, a1, a1_zt, a3, a3_zt, wp3_on_wp2, &           ! intent(in)
                   coef_wp4_implicit, wp4, wpthvp, wp2thvp, um, vm, &          ! intent(in)
                   upwp, vpwp, up2, vp2, Kw1, Kw8, Kh_zt,  &                   ! intent(in)
                   Skw_zt, tau1m, tauw3t, tau_C1_zm, C1_Skw_fnc, &             ! intent(in)
                   C11_Skw_fnc, C16_fnc, rho_ds_zm, invrs_rho_ds_zt, radf, &   ! intent(in)
                   thv_ds_zm, thv_ds_zt, wp2_splat, wp3_splat, &               ! intent(in)
                   l_crank_nich_diff, &                                        ! intent(in)
                   rhs )                                                       ! intent(out)

    ! Save the value of rhs, which will be overwritten with the solution as
    ! part of the solving routine.
    rhs_save = rhs

    if (l_gmres) then

       call wp23_gmres( dt, wp2, wm_zm, wm_zt, a1, a1_zt, a3, a3_zt, &
                        wp3_on_wp2, coef_wp4_implicit, &
                        Kw1, Kw8, Skw_zt, tau1m, tauw3t, tau_C1_zm, &
                        C1_Skw_fnc, C11_Skw_fnc, C16_fnc, rho_ds_zm, &
                        rho_ds_zt, invrs_rho_ds_zm, &
                        invrs_rho_ds_zt, l_crank_nich_diff, nrhs, &
                        rhs, &
                        solut, wp3_pr3_lhs )

    else

       ! Compute the implicit portion of the w'^2 and w'^3 equations.
       ! Build the left-hand side matrix.
       call wp23_lhs( dt, wp2, wm_zm, wm_zt, a1, a1_zt, a3, a3_zt,  &
                      wp3_on_wp2, coef_wp4_implicit, &
                      Kw1, Kw8, Skw_zt, tau1m, tauw3t, tau_C1_zm, C1_Skw_fnc, &
                      C11_Skw_fnc, C16_fnc, rho_ds_zm, rho_ds_zt, &
                      invrs_rho_ds_zm, invrs_rho_ds_zt, l_crank_nich_diff, & 
                      lhs, wp3_pr3_lhs )

       ! Solve the system with LAPACK
       if ( l_stats_samp .and. iwp23_matrix_condt_num > 0 ) then

            ! Perform LU decomp and solve system (LAPACK with diagnostics)
            ! Note that this can change the answer slightly
            call band_solvex( "wp2_wp3", 2, 2, 2*gr%nz, nrhs, & 
                              lhs, rhs, solut, rcond )
            
            if ( clubb_at_least_debug_level( 0 ) ) then
               if ( err_code == clubb_fatal_error ) then
                  write(fstderr,*) "Error in wp23_solve calling band_solvex for wp2_wp3"
                  write(fstderr,*) "wp2 & wp3 LU decomp. failed"
                  write(fstderr,*) "wp2 and wp3 LHS"
                  do k = 1, gr%nz
                     write(fstderr,*) "zt level = ", k, "height [m] = ", &
                                      gr%zt(k), "LHS = ", lhs(1:5,2*k-1)
                     write(fstderr,*) "zm level = ", k, "height [m] = ", &
                                      gr%zm(k), "LHS = ", lhs(1:5,2*k)
                  enddo ! k = 1, gr%nz
                  write(fstderr,*) "wp2 and wp3 RHS"
                  do k = 1, gr%nz
                     write(fstderr,*) "zt level = ", k, "height [m] = ", &
                                      gr%zt(k), "RHS = ", rhs_save(2*k-1)
                     write(fstderr,*) "zm level = ", k, "height [m] = ", &
                                      gr%zm(k), "RHS = ", rhs_save(2*k)
                  enddo ! k = 1, gr%nz
                  return
               endif
            endif

          ! Est. of the condition number of the w'^2/w^3 LHS matrix
          call stat_update_var_pt( iwp23_matrix_condt_num, 1, one / rcond, stats_sfc )

       else

            ! Perform LU decomp and solve system (LAPACK)
            call band_solve( "wp2_wp3", 2, 2, 2*gr%nz, nrhs, & 
                             lhs, rhs, solut )

            if ( clubb_at_least_debug_level( 0 ) ) then
               if ( err_code == clubb_fatal_error ) then
                  write(fstderr,*) "Error in wp23_solve calling band_solve for wp2_wp3"
                  write(fstderr,*) "wp2 & wp3 LU decomp. failed"
                  write(fstderr,*) "wp2 and wp3 LHS"
                  do k = 1, gr%nz
                     write(fstderr,*) "zt level = ", k, "height [m] = ", &
                                      gr%zt(k), "LHS = ", lhs(1:5,2*k-1)
                     write(fstderr,*) "zm level = ", k, "height [m] = ", &
                                      gr%zm(k), "LHS = ", lhs(1:5,2*k)
                  enddo ! k = 1, gr%nz
                  write(fstderr,*) "wp2 and wp3 RHS"
                  do k = 1, gr%nz
                     write(fstderr,*) "zt level = ", k, "height [m] = ", &
                                      gr%zt(k), "RHS = ", rhs_save(2*k-1)
                     write(fstderr,*) "zm level = ", k, "height [m] = ", &
                                      gr%zm(k), "RHS = ", rhs_save(2*k)
                  enddo ! k = 1, gr%nz
                  return
               endif
            endif

       endif

    endif ! l_gmres

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

    if ( l_stats_samp ) then

      ! Finalize implicit contributions for wp2

      do k = 2, gr%nz-1

        km1 = max( k-1, 1 )
        kp1 = min( k+1, gr%nz )

        ! w'^2 term dp1 has both implicit and explicit components;
        ! call stat_end_update_pt.
        call stat_end_update_pt( iwp2_dp1, k, & 
           zmscr01(k) * wp2(k), stats_zm )

        ! w'^2 term dp2 has both implicit and explicit components (if the
        ! Crank-Nicholson scheme is selected); call stat_end_update_pt.  
        ! If Crank-Nicholson diffusion is not selected, then w'^3 term dp1 is 
        ! completely implicit; call stat_update_var_pt.
        if ( l_crank_nich_diff ) then
           call stat_end_update_pt( iwp2_dp2, k, &
              zmscr02(k) * wp2(km1) & 
            + zmscr03(k) * wp2(k) & 
            + zmscr04(k) * wp2(kp1), stats_zm )
        else
           call stat_update_var_pt( iwp2_dp2, k, &
              zmscr02(k) * wp2(km1) & 
            + zmscr03(k) * wp2(k) & 
            + zmscr04(k) * wp2(kp1), stats_zm )
        endif

        ! w'^2 term ta is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwp2_ta, k, & 
           zmscr05(k) * wp3(k) & 
         + zmscr06(k) * wp3(kp1), stats_zm )

        ! w'^2 term ma is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwp2_ma, k, & 
           zmscr07(k) * wp2(km1) & 
         + zmscr08(k) * wp2(k) & 
         + zmscr09(k) * wp2(kp1), stats_zm )

        ! w'^2 term ac is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwp2_ac, k,  & 
           zmscr10(k) * wp2(k), stats_zm )

        ! w'^2 term pr1 has both implicit and explicit components;
        ! call stat_end_update_pt.
        if ( l_tke_aniso ) then
          call stat_end_update_pt( iwp2_pr1, k, & 
             zmscr12(k) * wp2(k), stats_zm )
        endif

        ! w'^2 term pr2 has both implicit and explicit components;
        ! call stat_end_update_pt.
        call stat_end_update_pt( iwp2_pr2, k, & 
           zmscr11(k) * wp2(k), stats_zm )

      enddo

      ! Finalize implicit contributions for wp3

      do k = 2, gr%nz-1, 1

        km1 = max( k-1, 1 )
        kp1 = min( k+1, gr%nz )

        ! w'^3 term pr1 has both implicit and explicit components; 
        ! call stat_end_update_pt.
        call stat_end_update_pt( iwp3_pr1, k, & 
           ztscr01(k) * wp3(k), stats_zt )

        ! w'^3 term dp1 has both implicit and explicit components (if the
        ! Crank-Nicholson scheme is selected); call stat_end_update_pt.  
        ! If Crank-Nicholson diffusion is not selected, then w'^3 term dp1 is 
        ! completely implicit; call stat_update_var_pt.
        if ( l_crank_nich_diff ) then
           call stat_end_update_pt( iwp3_dp1, k, & 
              ztscr02(k) * wp3(km1) & 
            + ztscr03(k) * wp3(k) & 
            + ztscr04(k) * wp3(kp1), stats_zt )
        else
           call stat_update_var_pt( iwp3_dp1, k, & 
              ztscr02(k) * wp3(km1) & 
            + ztscr03(k) * wp3(k) & 
            + ztscr04(k) * wp3(kp1), stats_zt )
        endif

        ! w'^3 term ta has both implicit and explicit components; 
        ! call stat_end_update_pt.
        call stat_end_update_pt( iwp3_ta, k, & 
           ztscr05(k) * wp3(km1) & 
         + ztscr06(k) * wp2(km1) & 
         + ztscr07(k) * wp3(k) & 
         + ztscr08(k) * wp2(k) & 
         + ztscr09(k) * wp3(kp1), stats_zt )

        ! w'^3 term tp has both implicit and explicit components; 
        ! call stat_end_update_pt.
        call stat_end_update_pt( iwp3_tp, k,  & 
           ztscr10(k) * wp2(km1) & 
         + ztscr11(k) * wp2(k), stats_zt )

        ! w'^3 pressure term 3 (pr3) has both implicit and explicit components;
        ! call stat_end_update_pt
        call stat_end_update_pt( iwp3_pr3, k, &
         - wp3_pr3_lhs(k,5) * wp3(km1) &
         - wp3_pr3_lhs(k,4) * wp2(km1) &
         - wp3_pr3_lhs(k,3) * wp3(k) &
         - wp3_pr3_lhs(k,2) * wp2(k) &
         - wp3_pr3_lhs(k,1) * wp3(kp1), stats_zt )

        ! w'^3 term ma is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwp3_ma, k, & 
           ztscr12(k) * wp3(km1) & 
         + ztscr13(k) * wp3(k) & 
         + ztscr14(k) * wp3(kp1), stats_zt )

        ! w'^3 term ac is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwp3_ac, k, & 
           ztscr15(k) * wp3(k), stats_zt )

        ! w'^3 term pr2 has both implicit and explicit components; 
        ! call stat_end_update_pt.
        call stat_end_update_pt( iwp3_pr2, k, & 
           ztscr16(k) * wp3(k), stats_zt )

      enddo

    endif ! l_stats_samp


    if ( l_stats_samp ) then
      ! Store previous value for effect of the positive definite scheme
      call stat_begin_update( iwp2_pd, wp2 / dt, stats_zm )
    endif

    if ( l_hole_fill .and. any( wp2 < w_tol_sqd ) ) then

      ! Use a simple hole filling algorithm
      call fill_holes_vertical( 2, w_tol_sqd, "zm", &
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
      call stat_end_update( iwp2_pd, wp2 / dt, stats_zm )
    endif


    ! Clip <w'^2> at a minimum threshold.

    ! The value of <w'^2> is not allowed to become smaller than the threshold
    ! value of w_tol^2.  Additionally, that threshold value may be boosted at
    ! any grid level in order to keep the overall correlation of w and rt or
    ! the overall correlation of w and theta-l between the values of
    ! -max_mag_correlation_flux and max_mag_correlation_flux by boosting <w'^2>
    ! rather than by limiting the magnitude of <w'rt'> or <w'thl'>.
    if ( l_min_wp2_from_corr_wx ) then

       ! The overall correlation of w and rt is:
       !
       ! corr_w_rt = wprtp / ( sqrt( wp2 ) * sqrt( rtp2 ) );
       !
       ! and the overall correlation of w and thl is:
       !
       ! corr_w_thl = wpthlp / ( sqrt( wp2 ) * sqrt( thlp2 ) ).
       !
       ! Squaring both sides, the equations becomes:
       !
       ! corr_w_rt^2 = wprtp^2 / ( wp2 * rtp2 ); and
       !
       ! corr_w_thl^2 = wpthlp^2 / ( wp2 * thlp2 ).
       !
       ! Using max_mag_correlation_flux for the correlation and then solving for
       ! the minimum of wp2, the equation becomes:
       !
       ! wp2|_min = max( wprtp^2 / ( rtp2 * max_mag_correlation_flux^2 ),
       !                 wpthlp^2 / ( thlp2 * max_mag_correlation_flux^2 ) ).
       do k = 1, gr%nz, 1

          threshold &
          = max( w_tol_sqd, &
                 wprtp(k)**2 / ( rtp2(k) * max_mag_correlation_flux**2 ), &
                 wpthlp(k)**2 / ( thlp2(k) * max_mag_correlation_flux**2 ) )

          call clip_variance_level( clip_wp2, dt, threshold, k, & ! In
                                    wp2(k) )                      ! In/out

       enddo ! k = 1, gr%nz, 1

    else

       ! Consider only the minimum tolerance threshold value for wp2.
       threshold = w_tol_sqd

       call clip_variance( clip_wp2, dt, threshold, & ! Intent(in)
                           wp2 )                      ! Intent(inout)

    endif ! l_min_wp2_from_corr_wx

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
                         wp3_on_wp2, coef_wp4_implicit, &
                         Kw1, Kw8, Skw_zt, tau1m, tauw3t, tau_C1_zm, &
                         C1_Skw_fnc, C11_Skw_fnc, C16_fnc, rho_ds_zm, &
                         rho_ds_zt, invrs_rho_ds_zm, &
                         invrs_rho_ds_zt, l_crank_nich_diff, nrhs, &
                         rhs, &
                         solut, wp3_pr3_lhs )
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
        core_rknd ! Variable(s)

#ifdef MKL

    use stats_variables, only:  & 
        iwp23_matrix_condt_num, & ! Variable(s)
        l_stats_samp, & 
        stats_sfc

    use constants_clubb, only: & 
        fstderr, & ! Variable(s)
        one

    use lapack_wrap, only:  & 
        band_solve,  & ! Procedure(s) 
        band_solvex

    use stats_type_utilities, only: & 
        stat_update_var_pt ! Procedure(s)

    use csr_matrix_module, only: &
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
    
    use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_no_error                 ! Constant

#endif /* MKL */

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      dt                 ! Timestep                                  [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      wp2                ! w'^2 (momentum levels)                    [m^2/s^2]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  & 
      wm_zm,             & ! w wind component on momentum levels       [m/s]
      wm_zt,             & ! w wind component on thermodynamic levels  [m/s]
      a1,                & ! a_1 (momentum levels); See eqn. 23 in `Equations for CLUBB' [-]
      a1_zt,             & ! a_1 interpolated to thermodynamic levels                    [-]
      a3,                & ! a_3 (momentum levels); See eqn. 25 in `Equations for CLUBB' [-]
      a3_zt,             & ! a_3 interpolated to thermodynamic levels  [-]
      wp3_on_wp2,        & ! Smoothed version of wp3 / wp2             [m/s]
      coef_wp4_implicit, & ! <w'^4> = coef_wp4_implicit * <w'^2>^2     [-]
      Kw1,               & ! Coefficient of eddy diffusivity for w'^2  [m^2/s]
      Kw8,               & ! Coefficient of eddy diffusivity for w'^3  [m^2/s]
      Skw_zt,            & ! Skewness of w on thermodynamic levels     [-]
      tau1m,             & ! Time-scale tau on momentum levels         [s]
      tauw3t,            & ! Time-scale tau on thermodynamic levels    [s]
      tau_C1_zm,         & ! Tau values used for the C1 (dp1) term in wp2 [s]
      C1_Skw_fnc,        & ! C_1 parameter with Sk_w applied           [-]
      C11_Skw_fnc,       & ! C_11 parameter with Sk_w applied          [-]
      C16_fnc,           & ! C_16 parameter                            [-]
      rho_ds_zm,         & ! Dry, static density on momentum levels    [kg/m^3]
      rho_ds_zt,         & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zm,   & ! Inv. dry, static density @ momentum levs. [m^3/kg]
      invrs_rho_ds_zt      ! Inv. dry, static density @ thermo. levs.  [m^3/kg]

    logical, intent(in) :: & 
      l_crank_nich_diff  ! Turns on/off Crank-Nicholson diffusion.

    integer, intent(in) :: &
      nrhs      ! Number of right-hand side vectors
                ! (GMRES currently only supports 1)

    ! Input/Output variables
    real( kind = core_rknd ), dimension(2*gr%nz), intent(inout) :: &
      rhs       ! Right hand side vector

    ! Output variables
    real( kind = core_rknd ), dimension(2*gr%nz), intent(out) :: &
      solut     ! Solution to band diagonal system

    real( kind = core_rknd ), dimension(gr%nz,5), intent(out) :: &
      wp3_pr3_lhs ! w'^3 pressure term 3 (pr3) lhs contribution

#ifdef MKL
    ! Local variables
    real( kind = core_rknd ), dimension(5,2*gr%nz) :: &
      lhs, &    ! Implicit contributions to wp2/wp3 (band diag. matrix)
      lhs_cache ! Backup cache of LHS matrix

    real( kind = core_rknd ), dimension(intlc_5d_5d_ja_size) :: &
      lhs_a_csr ! Implicit contributions to wp2/wp3 (CSR format)

    real( kind = core_rknd ), dimension(2*gr%nz) :: &
      rhs_cache ! Backup cache of RHS vector

    real( kind = core_rknd )::  & 
      rcond  ! Est. of the reciprocal of the condition #

    ! Begin code

    call wp23_lhs_csr( dt, wp2, wm_zm, wm_zt, a1, a1_zt, a3, a3_zt,  &
                       wp3_on_wp2, coef_wp4_implicit, &
                       Kw1, Kw8, Skw_zt, tau1m, tauw3t, tau_C1_zm, &
                       C1_Skw_fnc, C11_Skw_fnc, C16_fnc, rho_ds_zm, &
                       rho_ds_zt, invrs_rho_ds_zm, &
                       invrs_rho_ds_zt, l_crank_nich_diff, & 
                       lhs_a_csr, wp3_pr3_lhs )

    if ( .not. l_gmres_soln_ok(gmres_idx_wp2wp3) ) then
      call wp23_lhs( dt, wp2, wm_zm, wm_zt, a1, a1_zt, a3, a3_zt,  &
                     wp3_on_wp2, coef_wp4_implicit, &
                     Kw1, Kw8, Skw_zt, tau1m, tauw3t, tau_C1_zm, C1_Skw_fnc, &
                     C11_Skw_fnc, C16_fnc, rho_ds_zm, rho_ds_zt, &
                     invrs_rho_ds_zm, invrs_rho_ds_zt, l_crank_nich_diff, & 
                     lhs, wp3_pr3_lhs )

      ! Solve system with LAPACK to give us our first solution vector
        lhs_cache = lhs
        rhs_cache = rhs
        call band_solve( "wp2_wp3", 2, 2, 2*gr%nz, nrhs, &
                         lhs, rhs, solut )
        if ( clubb_at_least_debug_level( 0 ) ) then
            if ( err_code == clubb_fatal_error ) then
                write(fstderr,*) "in wp23_solve calling band_solve for wp2_wp3"
                return
            end if
        end if

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
                      solut )

    ! Fall back to LAPACK if GMRES returned any errors
    if ( err_code == clubb_fatal_error ) then
      write(fstderr,*) "Errors encountered in GMRES solve."
      write(fstderr,*) "Falling back to LAPACK solver."
      err_code = clubb_no_error

      ! Generate the LHS in LAPACK format
      call wp23_lhs( dt, wp2, wm_zm, wm_zt, a1, a1_zt, a3, a3_zt,  &
                     wp3_on_wp2, coef_wp4_implicit, &
                     Kw1, Kw8, Skw_zt, tau1m, tauw3t, tau_C1_zm, C1_Skw_fnc, &
                     C11_Skw_fnc, C16_fnc, rho_ds_zm, rho_ds_zt, &
                     invrs_rho_ds_zm, invrs_rho_ds_zt, l_crank_nich_diff, & 
                     lhs, wp3_pr3_lhs )

      ! Note: The RHS does not need to be re-generated.

      ! Solve the system with LAPACK as a fall-back.
      if ( l_stats_samp .and. iwp23_matrix_condt_num > 0 ) then

        ! Perform LU decomp and solve system (LAPACK with diagnostics)
        ! Note that this can change the answer slightly
        call band_solvex( "wp2_wp3", 2, 2, 2*gr%nz, nrhs, & 
                          lhs, rhs, solut, rcond )
        
        if ( clubb_at_least_debug_level( 0 ) ) then
            if ( err_code == clubb_fatal_error ) then
                write(fstderr,*) "in wp23_solve calling band_solvex for wp2_wp3"
                return
            end if
        end if

        ! Est. of the condition number of the w'^2/w^3 LHS matrix
        call stat_update_var_pt( iwp23_matrix_condt_num, 1, one / rcond, stats_sfc )

      else
        ! Perform LU decomp and solve system (LAPACK)
        call band_solve( "wp2_wp3", 2, 2, 2*gr%nz, nrhs, & 
                         lhs, rhs, solut )

        if ( clubb_at_least_debug_level( 0 ) ) then
            if ( err_code == clubb_fatal_error ) then
                write(fstderr,*) "in wp23_solve calling band_solve for wp2_wp3"
                return
            end if
        end if

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
    solut(1:gr%nz) = C11_Skw_fnc + C16_fnc
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
    solut(1:gr%nz) = tau_C1_zm
    solut(1:gr%nz) = wm_zt
    solut(1:gr%nz) = wm_zm
    solut(1:gr%nz) = wp2
    solut(1:gr%nz) = wp3_on_wp2
    wp3_pr3_lhs = -9999._core_rknd

#endif /* MKL */

  end subroutine wp23_gmres

  !=================================================================================
  subroutine wp23_lhs( dt, wp2, wm_zm, wm_zt, a1, a1_zt, a3, a3_zt,  &
                       wp3_on_wp2, coef_wp4_implicit, &
                       Kw1, Kw8, Skw_zt, tau1m, tauw3t, tau_C1_zm, C1_Skw_fnc, &
                       C11_Skw_fnc, C16_fnc, rho_ds_zm, rho_ds_zt, &
                       invrs_rho_ds_zm, invrs_rho_ds_zt, l_crank_nich_diff, & 
                       lhs, wp3_pr3_lhs )
    ! Description:
    ! Compute LHS band diagonal matrix for w'^2 and w'^3.
    ! This subroutine computes the implicit portion 
    ! of the w'^2 and w'^3 equations.
    !
    ! NOTE: If changes are made to this subroutine, ensure that the CSR
    !   version of the subroutine is updated as well! If the two are different,
    !   the results will be inconsistent between LAPACK and PARDISO/GMRES!
    ! 
    ! 
    ! Boundary conditions
    ! 
    !   Both wp2 and wp3 used fixed-point boundary conditions.
    !   Therefore, anything set in the above loop at both the upper
    !   and lower boundaries would be overwritten here.  However, the
    !   above loop does not extend to the boundary levels.  An array
    !   with a value of 1 at the main diagonal on the left-hand side
    !   and with values of 0 at all other diagonals on the left-hand
    !   side will preserve the right-hand side value at that level.
    !
    !      wp3(1)  wp2(1)  ... wp3(nzmax) wp2(nzmax)
    !     [  0.0     0.0          0.0       0.0  ]
    !     [  0.0     0.0          0.0       0.0  ]
    !     [  1.0     1.0   ...    1.0       1.0  ]
    !     [  0.0     0.0          0.0       0.0  ]
    !     [  0.0     0.0          0.0       0.0  ]
    ! 
    ! 
    !  WARNING: This subroutine has been optimized. Significant changes could
    !           noticeably  impact computational efficiency. See clubb:ticket:834
    !-------------------------------------------------------------------------------

    use grid_class, only:  & 
        gr ! Variable

    use parameters_tunable, only:  & 
        C4,  & ! Variables
        C5,  & 
        C8,  & 
        C8b, & 
        C12, & 
        nu1_vert_res_dep, & 
        nu8_vert_res_dep

    use constants_clubb, only:  & 
        one, &
        one_half, &
        gamma_over_implicit_ts, &
        zero

    use model_flags, only: & 
        l_tke_aniso,                  & ! Variable(s)
        l_explicit_turbulent_adv_wp3, &
        l_use_wp3_pr3

    use diffusion, only: & 
        diffusion_zm_lhs,  & ! Procedures
        diffusion_zm_lhs_all, &
        diffusion_zt_lhs, &
        diffusion_zt_lhs_all

    use mean_adv, only: & 
        term_ma_zm_lhs,  & ! Procedures
        term_ma_zm_lhs_all, &
        term_ma_zt_lhs, &
        term_ma_zt_lhs_all

    use pdf_closure_module, only: &
        iiPDF_ADG1, & ! Variable(s)
        iiPDF_new,  &
        iiPDF_type

    use clubb_precision, only: &
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
        ztscr16

    use stats_variables, only: & 
        l_stats_samp, & 
        iwp2_dp1, & 
        iwp2_dp2, & 
        iwp2_ta, & 
        iwp2_ma, & 
        iwp2_ac, & 
        iwp2_pr2, & 
        iwp2_pr1, &
        iwp3_ta, & 
        iwp3_tp, & 
        iwp3_ma, & 
        iwp3_ac, & 
        iwp3_pr2, & 
        iwp3_pr1, & 
        iwp3_dp1

    use advance_helper_module, only: set_boundary_conditions_lhs ! Procedure(s)

    implicit none


    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      dt                 ! Timestep length                            [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      wp2,               & ! w'^2 (momentum levels)                    [m^2/s^2]
      wm_zm,             & ! w wind component on momentum levels       [m/s]
      wm_zt,             & ! w wind component on thermodynamic levels  [m/s]
      a1,                & ! sigma_sqd_w term a_1 (momentum levels)    [-]
      a1_zt,             & ! a_1 interpolated to thermodynamic levels  [-]
      a3,                & ! sigma_sqd_w term a_3 (momentum levels)    [-]
      a3_zt,             & ! a_3 interpolated to thermodynamic levels  [-]
      wp3_on_wp2,        & ! Smoothed version of wp3 / wp2             [m/s]
      coef_wp4_implicit, & ! <w'^4> = coef_wp4_implicit * <w'^2>^2     [-]
      Kw1,               & ! Coefficient of eddy diffusivity for w'^2  [m^2/s]
      Kw8,               & ! Coefficient of eddy diffusivity for w'^3  [m^2/s]
      Skw_zt,            & ! Skewness of w on thermodynamic levels     [-]
      tau1m,             & ! Time-scale tau on momentum levels         [s]
      tauw3t,            & ! Time-scale tau on thermodynamic levels    [s]
      tau_C1_zm,         & ! Tau values used for the C1 (dp1) term in wp2 [s]
      C1_Skw_fnc,        & ! C_1 parameter with Sk_w applied           [-]
      C11_Skw_fnc,       & ! C_11 parameter with Sk_w applied          [-]
      C16_fnc,           & ! C_16 parameter                            [-]
      rho_ds_zm,         & ! Dry, static density on momentum levels    [kg/m^3]
      rho_ds_zt,         & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zm,   & ! Inv. dry, static density @ momentum levs. [m^3/kg]
      invrs_rho_ds_zt      ! Inv. dry, static density @ thermo. levs.  [m^3/kg]

    logical, intent(in) :: & 
      l_crank_nich_diff  ! Turns on/off Crank-Nicholson diffusion.

    ! Output Variable
    real( kind = core_rknd ), dimension(5,2*gr%nz), intent(out) ::  & 
      lhs ! Implicit contributions to wp2/wp3 (band diag. matrix)

    real( kind = core_rknd ), dimension(gr%nz,5), intent(out) :: &
      wp3_pr3_lhs

    ! Local Variables

    ! Loop Variable
    integer :: k, k_wp2, k_wp3

    real( kind = core_rknd ), dimension(5,gr%nz) :: &
      wp3_term_ta_lhs_result

    real( kind = core_rknd ), dimension(3,gr%nz) :: &
        lhs_diff_zm, &  ! Completely implicit diffusion term for w'2
        lhs_diff_zt, &  ! Completely implicit diffusion term for w'3
        lhs_ma_zm, &    ! Mean advection term for w'2
        lhs_ma_zt       ! Mean advection term for w'3

    real( kind = core_rknd ), dimension(2,gr%nz) :: &
        lhs_ta_wp2, &   ! Turbulent advection terms for wp2
        lhs_ta_wp3, &   ! Turbulent advection terms for wp3
        lhs_tp_wp3      ! Turbulent production terms of w'^3

    real( kind = core_rknd ), dimension(gr%nz) :: &
        lhs_ac_pr2_wp2, &   ! Accumulation terms of w'^2 and w'^2 pressure term 2
        lhs_ac_pr2_wp3, &   ! Accumulation terms of w'^3 and w'^3 pressure term 2
        lhs_dp1_wp2, &      ! Dissipation terms 1 for w'^2
        lhs_pr1_wp3, &      ! Dissipation terms 1 for w'^3
        lhs_pr1_wp2         ! Pressure term 1 for w'2

    real( kind = core_rknd) :: &
        invrs_dt        ! Inverse of dt, 1/dt, used for computational efficiency

    !---------------------- Being Code ----------------------


    ! Initialize arrays to 0 and calculate invrs_dt
    lhs = 0.0_core_rknd
    wp3_pr3_lhs = 0.0_core_rknd
    wp3_term_ta_lhs_result = 0.0_core_rknd
    invrs_dt = 1.0_core_rknd / dt


    ! Calculated mean advection term for w'2
    call term_ma_zm_lhs_all( wm_zm(:), gr%invrs_dzm(:), &
                             lhs_ma_zm(:,:) )


    ! Calculated mean advection term for w'3
    call term_ma_zt_lhs_all( wm_zt(:), gr%invrs_dzt(:), gr%invrs_dzm(:), &
                             lhs_ma_zt(:,:) )


    ! Calculate diffusion term for w'2 using a completely implicit time step
    call diffusion_zm_lhs_all( Kw1(:), nu1_vert_res_dep(:), & 
                               gr%invrs_dzt(:), gr%invrs_dzm(:), &
                               lhs_diff_zm(:,:) )


    ! Calculate diffusion term for w'3 using a completely implicit time step
    call diffusion_zt_lhs_all( Kw8(:), nu8_vert_res_dep(:), & 
                               gr%invrs_dzm(:), gr%invrs_dzt(:), &
                               lhs_diff_zt(:,:) )

    lhs_diff_zt(:,:) = lhs_diff_zt(:,:) * C12

    if ( l_crank_nich_diff ) then

        ! Using a Crank-Nicholson time step for diffusion terms
        ! Modify diffusion terms
        do k = 2, gr%nz - 1

            lhs_diff_zm(1,k) = lhs_diff_zm(1,k) * 0.5_core_rknd
            lhs_diff_zm(2,k) = lhs_diff_zm(2,k) * 0.5_core_rknd
            lhs_diff_zm(3,k) = lhs_diff_zm(3,k) * 0.5_core_rknd

            lhs_diff_zt(1,k) = lhs_diff_zt(1,k) * 0.5_core_rknd
            lhs_diff_zt(2,k) = lhs_diff_zt(2,k) * 0.5_core_rknd
            lhs_diff_zt(3,k) = lhs_diff_zt(3,k) * 0.5_core_rknd

        end do

    end if


    ! Calculate turbulent advection terms for wp2
    call wp2_term_ta_lhs_all( rho_ds_zt(:), &
                              invrs_rho_ds_zm(:), &
                              gr%invrs_dzm(:), &
                              lhs_ta_wp2(:,:) )


    ! Calculate accumulation terms of w'^2 and w'^2 pressure term 2
    call wp2_terms_ac_pr2_lhs_all( C5, wm_zt(:), gr%invrs_dzm(:), &
                                   lhs_ac_pr2_wp2(:) )


    ! Calculate dissipation terms 1 for w'^2
    call wp2_term_dp1_lhs_all( C1_Skw_fnc(:), tau_C1_zm(:), &
                               lhs_dp1_wp2(:) )


    ! Calculate turbulent production terms of w'^3
    call wp3_term_tp_lhs_all( wp2(:), &
                              rho_ds_zm(:), &
                              invrs_rho_ds_zt(:), &
                              gr%invrs_dzt(:), &
                              lhs_tp_wp3(:,:) )


    ! Calculate accumulation terms of w'^3 and w'^3 pressure terms 2
    call wp3_terms_ac_pr2_lhs_all( C11_Skw_fnc(:), wm_zm(:), gr%invrs_dzt(:), &
                                   lhs_ac_pr2_wp3(:) )


    ! Calculate pressure terms 1 for w'^3
    call wp3_term_pr1_lhs_all( C8, C8b, tauw3t(:), Skw_zt(:), &
                               lhs_pr1_wp3(:) )


    ! Lower boundary for w'3
    lhs(1,1) = 0.0_core_rknd
    lhs(2,1) = 0.0_core_rknd
    lhs(3,1) = 1.0_core_rknd
    lhs(4,1) = 0.0_core_rknd
    lhs(5,1) = 0.0_core_rknd

    ! Lower boundary for w'2
    lhs(1,2) = 0.0_core_rknd
    lhs(2,2) = 0.0_core_rknd
    lhs(3,2) = 1.0_core_rknd
    lhs(4,2) = 0.0_core_rknd
    lhs(5,2) = 0.0_core_rknd

    ! Combine terms to calculate non-boundary lhs values
    do k = 2, gr%nz-1, 1

        k_wp3 = 2*k - 1
        k_wp2 = 2*k

        ! ------ w'3 ------

        ! LHS mean advection (ma) and diffusion (diff) terms
        lhs(1,k_wp3) = lhs(1,k_wp3) + lhs_ma_zt(1,k) + lhs_diff_zt(1,k)

        ! LHS turbulent production (tp) term.
        ! Note:  An "over-implicit" weighted time step is applied to this term.
        lhs(2,k_wp3) = lhs(2,k_wp3) + gamma_over_implicit_ts * lhs_tp_wp3(1,k)

        ! LHS mean advection (ma) and diffusion (diff) terms
        lhs(3,k_wp3) = lhs(3,k_wp3) + lhs_ma_zt(2,k) + lhs_diff_zt(2,k)
                                    
        ! LHS accumulation (ac) term and pressure term 2 (pr2).
        lhs(3,k_wp3) = lhs(3,k_wp3) + lhs_ac_pr2_wp3(k)

        ! LHS pressure term 1 (pr1).
        ! Note:  An "over-implicit" weighted time step is applied to this term.
        lhs(3,k_wp3) = lhs(3,k_wp3) + gamma_over_implicit_ts * lhs_pr1_wp3(k)

        ! LHS time tendency.
        lhs(3,k_wp3) = lhs(3,k_wp3) + invrs_dt

        ! LHS turbulent production (tp) term.
        ! Note:  An "over-implicit" weighted time step is applied to this term.
        lhs(4,k_wp3) = lhs(4,k_wp3) + gamma_over_implicit_ts * lhs_tp_wp3(2,k)

        ! LHS mean advection (ma) and diffusion (diff) terms
        lhs(5,k_wp3) = lhs(5,k_wp3) + lhs_ma_zt(3,k) + lhs_diff_zt(3,k)


        ! ------ w'2 ------

        ! LHS mean advection (ma) and diffusion (diff) terms
        lhs(1,k_wp2) = lhs(1,k_wp2) + lhs_ma_zm(1,k) + lhs_diff_zm(1,k)

        ! LHS turbulent advection (ta) term.
        lhs(2,k_wp2) = lhs(2,k_wp2) + lhs_ta_wp2(1,k)

        ! LHS mean advection (ma) and diffusion (diff) terms
        lhs(3,k_wp2) = lhs(3,k_wp2) + lhs_ma_zm(2,k) + lhs_diff_zm(2,k) 
                                    
        ! LHS accumulation (ac) term and pressure term 2 (pr2).
        lhs(3,k_wp2) = lhs(3,k_wp2) + lhs_ac_pr2_wp2(k)

        ! LHS dissipation term 1 (dp1).
        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the term
        !        more numerically stable (see note below for w'^3 LHS turbulent
        !        advection (ta) term).
        lhs(3,k_wp2) = lhs(3,k_wp2) + gamma_over_implicit_ts  * lhs_dp1_wp2(k)

        ! LHS time tendency.
        lhs(3,k_wp2) = lhs(3,k_wp2) + invrs_dt

        ! LHS turbulent advection (ta) term.
        lhs(4,k_wp2) = lhs(4,k_wp2) + lhs_ta_wp2(2,k)

        ! LHS mean advection (ma) and diffusion (diff) terms
        lhs(5,k_wp2) = lhs(5,k_wp2) + lhs_ma_zm(3,k) + lhs_diff_zm(3,k)

    enddo

    ! Upper boundary for w'3
    lhs(1,2*gr%nz-1) = 0.0_core_rknd
    lhs(2,2*gr%nz-1) = 0.0_core_rknd
    lhs(3,2*gr%nz-1) = 1.0_core_rknd
    lhs(4,2*gr%nz-1) = 0.0_core_rknd
    lhs(5,2*gr%nz-1) = 0.0_core_rknd

    ! Upper boundary for w'2
    lhs(1,2*gr%nz) = 0.0_core_rknd
    lhs(2,2*gr%nz) = 0.0_core_rknd
    lhs(3,2*gr%nz) = 1.0_core_rknd
    lhs(4,2*gr%nz) = 0.0_core_rknd
    lhs(5,2*gr%nz) = 0.0_core_rknd


    ! LHS pressure term 1 (pr1) for wp2
    if ( l_tke_aniso ) then

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the term
        !        more numerically stable (see note below for w'^3 LHS turbulent
        !        advection (ta) term).
        ! Reference:
        ! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:wp2_pr 

        ! Calculate terms
        call wp2_term_pr1_lhs_all( C4, tau1m(:), &
                                   lhs_pr1_wp2(:) )

        ! Add terms to lhs
        do k = 2, gr%nz-1

            k_wp2 = 2*k

            lhs(3,k_wp2) = lhs(3,k_wp2) + gamma_over_implicit_ts * lhs_pr1_wp2(k)

        end do

    endif


    ! LHS turbulent advection (ta) term for wp3
    if ( .not. l_explicit_turbulent_adv_wp3 ) then
        
        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        The weight of the implicit portion of this term is controlled
        !        by the factor gamma_over_implicit_ts (abbreviated "gamma" in
        !        the expression below).  A factor is added to the right-hand
        !        side of the equation in order to balance a weight that is not
        !        equal to 1, such that:
        !             -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
        !        where X is the variable that is being solved for in a
        !        predictive equation (w'^3 in this case), y(t) is the
        !        linearized portion of the term that gets treated implicitly,
        !        and RHS is the portion of the term that is always treated
        !        explicitly (in the case of the w'^3 turbulent advection term,
        !        RHS = 0).  A weight of greater than 1 can be applied to make
        !        the term more numerically stable.

        if ( iiPDF_type == iiPDF_ADG1 ) then

            ! The ADG1 PDF is used.

            ! Calculate terms
            call wp3_term_ta_ADG1_lhs_all( wp2(:), &
                                           a1(:), a1_zt(:), &
                                           a3(:), a3_zt(:), &
                                           wp3_on_wp2(:), &
                                           rho_ds_zm(:), &
                                           invrs_rho_ds_zt(:), &
                                           gr%invrs_dzt(:), &
                                           wp3_term_ta_lhs_result(:,:) )

        elseif ( iiPDF_type == iiPDF_new ) then

            ! The new PDF is used.

            ! Calculate terms
            call wp3_term_ta_new_pdf_lhs_all( coef_wp4_implicit(:), &
                                                 wp2(:), rho_ds_zm(:), &
                                                 invrs_rho_ds_zt(:), &
                                                 gr%invrs_dzt(:), &
                                                 lhs_ta_wp3(:,:) )

            ! Save terms in wp3_term_ta_lhs_result
            wp3_term_ta_lhs_result((/2,4/),:) = lhs_ta_wp3(:,:)

        endif

        ! Add terms to lhs
        do k = 2, gr%nz-1

            k_wp3 = 2*k - 1

            lhs(:,k_wp3) = lhs(:,k_wp3) + gamma_over_implicit_ts * wp3_term_ta_lhs_result(:,k)

        end do

    endif


    ! LHS pressure term 3 (pr3) for wp3
    if ( l_use_wp3_pr3 ) then

        do k = 2, gr%nz-1

            k_wp3 = 2*k - 1

            wp3_pr3_lhs(k,1) = - gamma_over_implicit_ts * C16_fnc(k) &
                               * wp3_term_ta_lhs_result(1,k)

            wp3_pr3_lhs(k,2) = - gamma_over_implicit_ts * C16_fnc(k) &
                               * ( wp3_term_ta_lhs_result(2,k) + lhs_tp_wp3(1,k) )

            wp3_pr3_lhs(k,3) = - gamma_over_implicit_ts * C16_fnc(k) &
                               * wp3_term_ta_lhs_result(3,k)

            wp3_pr3_lhs(k,4) = - gamma_over_implicit_ts * C16_fnc(k) &
                               * ( wp3_term_ta_lhs_result(4,k) + lhs_tp_wp3(2,k) )

            wp3_pr3_lhs(k,5) = - gamma_over_implicit_ts * C16_fnc(k) &
                               * wp3_term_ta_lhs_result(5,k)

            lhs(:,k_wp3) = lhs(:,k_wp3) + wp3_pr3_lhs(k,:)

        end do

    endif


    ! --------- Statistics output ---------
    if ( l_stats_samp ) then

        do k = 2, gr%nz-1

            !!!!!***** w'^2 *****!!!!!

            ! Note:  An "over-implicit" weighted time step is applied to this term.
            !        A weighting factor of greater than 1 may be used to make the
            !        term more numerically stable (see note below for w'^3 LHS
            !        turbulent advection (ta) term).
            if ( iwp2_dp1 > 0 ) then
                zmscr01(k) = - gamma_over_implicit_ts  * lhs_dp1_wp2(k)
            endif

            ! Eddy diffusion for wp2
            if ( iwp2_dp2 > 0 ) then
                zmscr02(k) = - lhs_diff_zm(3,k)
                zmscr03(k) = - lhs_diff_zm(2,k)
                zmscr04(k) = - lhs_diff_zm(1,k)
            endif

            ! Turbulent advection for wp2
            if ( iwp2_ta > 0 ) then
                zmscr05(k) = - lhs_ta_wp2(2,k)
                zmscr06(k) = - lhs_ta_wp2(1,k)
            endif

            ! Mean advection for wp2
            if ( iwp2_ma > 0 ) then
                zmscr07(k) = - lhs_ma_zm(3,k)
                zmscr08(k) = - lhs_ma_zm(2,k)
                zmscr09(k) = - lhs_ma_zm(1,k)
            endif

            ! Note:  To find the contribution of w'^2 term ac, substitute 0 for the
            !        C_5 input to function wp2_terms_ac_pr2_lhs.
            if ( iwp2_ac > 0 ) then
                zmscr10(k) = - wp2_terms_ac_pr2_lhs( 0.0_core_rknd, wm_zt(k+1), &
                                                     wm_zt(k), gr%invrs_dzm(k)  )
            endif

            ! Note:  To find the contribution of w'^2 term pr2, add 1 to the
            !        C_5 input to function wp2_terms_ac_pr2_lhs.
            if ( iwp2_pr2 > 0 ) then
                zmscr11(k) = - wp2_terms_ac_pr2_lhs( (one+C5), wm_zt(k+1), wm_zt(k),  & 
                                                      gr%invrs_dzm(k)  )
            endif

            ! Note:  An "over-implicit" weighted time step is applied to this term.
            !        A weighting factor of greater than 1 may be used to make the
            !        term more numerically stable (see note below for w'^3 LHS
            !        turbulent advection (ta) term).
            if ( iwp2_pr1 > 0 .and. l_tke_aniso ) then
                zmscr12(k) = - gamma_over_implicit_ts * lhs_pr1_wp2(k)
            endif



            !!!!!***** w'^3 *****!!!!!

            ! Turbulent advection for wp3
            if ( iwp3_ta > 0 ) then
                ztscr05(k) = - gamma_over_implicit_ts * wp3_term_ta_lhs_result(5,k)
                ztscr06(k) = - gamma_over_implicit_ts * wp3_term_ta_lhs_result(4,k)
                ztscr07(k) = - gamma_over_implicit_ts * wp3_term_ta_lhs_result(3,k)
                ztscr08(k) = - gamma_over_implicit_ts * wp3_term_ta_lhs_result(2,k)
                ztscr09(k) = - gamma_over_implicit_ts * wp3_term_ta_lhs_result(1,k)
            endif

            ! Note:  An "over-implicit" weighted time step is applied to this term.
            !        A weighting factor of greater than 1 may be used to make the
            !        term more numerically stable (see note above for LHS turbulent
            !        advection (ta) term).
            if ( iwp3_tp > 0 ) then
                ztscr10(k) = - gamma_over_implicit_ts * lhs_tp_wp3(2,k)
                ztscr11(k) = - gamma_over_implicit_ts * lhs_tp_wp3(1,k)
            endif

            ! Mean advection for wp2
            if ( iwp3_ma > 0 ) then
                ztscr12(k) = - lhs_ma_zt(3,k)
                ztscr13(k) = - lhs_ma_zt(2,k)
                ztscr14(k) = - lhs_ma_zt(1,k)
            endif

            ! Note:  To find the contribution of w'^3 term ac, substitute 0 for the
            !        C_ll skewness function input to function wp3_terms_ac_pr2_lhs.
            if ( iwp3_ac > 0 ) then
                ztscr15(k) = - wp3_terms_ac_pr2_lhs( 0.0_core_rknd, wm_zm(k), &
                                                     wm_zm(k-1), gr%invrs_dzt(k) )
            endif

            ! Note:  To find the contribution of w'^3 term pr2, add 1 to the
            !        C_ll skewness function input to function wp3_terms_ac_pr2_lhs.
            if ( iwp3_pr2 > 0 ) then
                ztscr16(k) = - wp3_terms_ac_pr2_lhs( (one+C11_Skw_fnc(k)), wm_zm(k), &
                                                     wm_zm(k-1), gr%invrs_dzt(k) )
            endif

            ! Note:  An "over-implicit" weighted time step is applied to this term.
            !        A weighting factor of greater than 1 may be used to make the
            !        term more numerically stable (see note above for LHS turbulent
            !        advection (ta) term).
            if ( iwp3_pr1 > 0 ) then
                ztscr01(k) = - gamma_over_implicit_ts  * lhs_pr1_wp3(k)
            endif

            ! Eddy diffusion for wp3 
            if ( iwp3_dp1 > 0 ) then
                ztscr02(k) = - lhs_diff_zt(3,k)
                ztscr03(k) = - lhs_diff_zt(2,k)
                ztscr04(k) = - lhs_diff_zt(1,k)
            endif

        end do

    end if

    return

  end subroutine wp23_lhs

#ifdef MKL
  !=============================================================================
  subroutine wp23_lhs_csr( dt, wp2, wm_zm, wm_zt, a1, a1_zt, a3, a3_zt,  &
                           wp3_on_wp2, coef_wp4_implicit, &
                           Kw1, Kw8, Skw_zt, tau1m, tauw3t, tau_C1_zm, &
                           C1_Skw_fnc, C11_Skw_fnc, C16_fnc, rho_ds_zm, &
                           rho_ds_zt, invrs_rho_ds_zm, &
                           invrs_rho_ds_zt, l_crank_nich_diff, & 
                           lhs_a_csr, wp3_pr3_lhs )

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
        nu8_vert_res_dep

    use constants_clubb, only:  & 
        eps,          & ! Variable(s)
        one,          &
        one_half,     &
        gamma_over_implicit_ts

    use model_flags, only: & 
        l_tke_aniso,                  & ! Variable(s)
        l_explicit_turbulent_adv_wp3, &
        l_use_wp3_pr3

    use diffusion, only: & 
        diffusion_zm_lhs,  & ! Procedures
        diffusion_zt_lhs

    use mean_adv, only: & 
        term_ma_zm_lhs,  & ! Procedures
        term_ma_zt_lhs

    use pdf_closure_module, only: &
        iiPDF_ADG1, & ! Variable(s)
        iiPDF_new,  &
        iiPDF_type

    use clubb_precision, only: &
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
        ztscr16

    use stats_variables, only: & 
        l_stats_samp, & 
        iwp2_dp1, & 
        iwp2_dp2, & 
        iwp2_ta, & 
        iwp2_ma, & 
        iwp2_ac, & 
        iwp2_pr2, & 
        iwp2_pr1, &
        iwp3_ta, & 
        iwp3_tp, & 
        iwp3_ma, & 
        iwp3_ac, & 
        iwp3_pr2, & 
        iwp3_pr1, & 
        iwp3_dp1

    use csr_matrix_module, only: &
        intlc_5d_5d_ja_size ! Variable

    implicit none

    ! Left-hand side matrix diagonal identifiers for
    ! momentum-level variable, w'^2.
    ! These are updated for each diagonal of the matrix as the
    ! LHS of the matrix is created.
    integer ::  &
      m_kp1_mdiag, & ! Momentum super diagonal index for w'^2.
      m_kp1_tdiag, & ! Thermodynamic super diagonal index for w'^2.
      m_k_mdiag  , & ! Momentum main diagonal index for w'^2.
      m_k_tdiag  , & ! Thermodynamic sub diagonal index for w'^2.
      m_km1_mdiag    ! Momentum sub diagonal index for w'^2.

    ! Left-hand side matrix diagonal identifiers for
    ! thermodynamic-level variable, w'^3.
    ! These are updated for each diagonal of the matrix as the
    ! LHS of the matrix is created
    integer ::  &
      t_kp1_tdiag, & ! Thermodynamic super diagonal index for w'^3.
      t_k_mdiag  , & ! Momentum super diagonal index for w'^3.
      t_k_tdiag  , & ! Thermodynamic main diagonal index for w'^3.
      t_km1_mdiag, & ! Momentum sub diagonal index for w'^3.
      t_km1_tdiag    ! Thermodynamic sub diagonal index for w'^3.

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      dt                 ! Timestep length                            [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      wp2,               & ! w'^2 (momentum levels)                    [m^2/s^2]
      wm_zm,             & ! w wind component on momentum levels       [m/s]
      wm_zt,             & ! w wind component on thermodynamic levels  [m/s]
      a1,                & ! sigma_sqd_w term a_1 (momentum levels)    [-]
      a1_zt,             & ! a_1 interpolated to thermodynamic levels  [-]
      a3,                & ! sigma_sqd_w term a_3 (momentum levels)    [-]
      a3_zt,             & ! a_3 interpolated to thermodynamic levels  [-]
      wp3_on_wp2,        & ! Smoothed version of wp3 / wp2             [m/s]
      coef_wp4_implicit, & ! <w'^4> = coef_wp4_implicit * <w'^2>^2     [-]
      Kw1,               & ! Coefficient of eddy diffusivity for w'^2  [m^2/s]
      Kw8,               & ! Coefficient of eddy diffusivity for w'^3  [m^2/s]
      Skw_zt,            & ! Skewness of w on thermodynamic levels     [-]
      tau1m,             & ! Time-scale tau on momentum levels         [s]
      tauw3t,            & ! Time-scale tau on thermodynamic levels    [s]
      tau_C1_zm,         & ! Tau values used for the C1 (dp1) term in wp2 [s]
      C1_Skw_fnc,        & ! C_1 parameter with Sk_w applied           [-]
      C11_Skw_fnc,       & ! C_11 parameter with Sk_w applied          [-]
      rho_ds_zm,         & ! Dry, static density on momentum levels    [kg/m^3]
      rho_ds_zt,         & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zm,   & ! Inv. dry, static density @ momentum levs. [m^3/kg]
      invrs_rho_ds_zt      ! Inv. dry, static density @ thermo. levs.  [m^3/kg]

    logical, intent(in) :: & 
      l_crank_nich_diff  ! Turns on/off Crank-Nicholson diffusion.

    ! Output Variable
    real( kind = core_rknd ), dimension(intlc_5d_5d_ja_size), intent(out) ::  & 
      lhs_a_csr ! Implicit contributions to wp2/wp3 (band diag. matrix)

    real( kind = core_rknd ), dimension(gr%nz,5), intent(out) :: &
      wp3_pr3_lhs

    ! Local Variables

    ! Array indices
    integer :: k, km1, kp1, k_wp2, k_wp3, wp2_cur_row, wp3_cur_row

    real( kind = core_rknd ), dimension(5) :: &
      tmp, &
      wp3_terms_ta_tp_lhs_result, &
      wp3_term_ta_lhs_result, &
      wp3_term_tp_lhs_result


    ! Initialize the left-hand side matrix to 0.
    lhs_a_csr = 0.0_core_rknd

    ! Initialize values to 0.
    wp3_term_ta_lhs_result = zero
    wp3_term_tp_lhs_result = zero
    wp3_terms_ta_tp_lhs_result = zero

    do k = 2, gr%nz-1, 1

      ! Define indices

      km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nz )

      k_wp3 = 2*k - 1
      k_wp2 = 2*k

      wp2_cur_row = ((k_wp2 - 3) * 5) + 8
      wp3_cur_row = ((k_wp3 - 3) * 5) + 8

      !!!!!***** w'^2 *****!!!!!

      ! w'^2: Left-hand side (implicit w'^2 portion of the code).
      !
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
      = real( + one / dt )

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
      !        advection (ta) term).
      lhs_a_csr(m_k_mdiag)  & 
      = lhs_a_csr(m_k_mdiag)  &
      + gamma_over_implicit_ts  & 
      * wp2_term_dp1_lhs( C1_Skw_fnc(k), tau_C1_zm(k) )

      ! LHS eddy diffusion term: dissipation term 2 (dp2).
      if ( l_crank_nich_diff ) then
        ! Eddy diffusion for wp2 using a Crank-Nicholson time step.
        lhs_a_csr((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/)) & 
        = lhs_a_csr((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/)) & 
        + one_half & 
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
      !        advection (ta) term).
      if ( l_tke_aniso ) then
        ! Add in this term if we're not assuming tke = 1.5 * wp2
        lhs_a_csr(m_k_mdiag)  & 
        = lhs_a_csr(m_k_mdiag)  &
        + gamma_over_implicit_ts  & 
        * wp2_term_pr1_lhs( C4, tau1m(k) )
      endif

      if ( l_stats_samp ) then

        ! Statistics: implicit contributions for wp2.

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note below for w'^3 LHS
        !        turbulent advection (ta) term).
        if ( iwp2_dp1 > 0 ) then
          zmscr01(k)  &
          = - gamma_over_implicit_ts  &
            * wp2_term_dp1_lhs( C1_Skw_fnc(k), tau_C1_zm(k) )
        endif

        if ( iwp2_dp2 > 0 ) then
          if ( l_crank_nich_diff ) then
            ! Eddy diffusion for wp2 using a Crank-Nicholson time step.
            tmp(1:3) & 
            = one_half & 
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
          - wp2_terms_ac_pr2_lhs( (one+C5), wm_zt(kp1), wm_zt(k),  & 
                                  gr%invrs_dzm(k)  )
        endif

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note below for w'^3 LHS
        !        turbulent advection (ta) term).
        if ( iwp2_pr1 > 0 .and. l_tke_aniso ) then
          zmscr12(k)  &
          = - gamma_over_implicit_ts  &
            * wp2_term_pr1_lhs( C4, tau1m(k) )
        endif

      endif



      !!!!!***** w'^3 *****!!!!!

      ! w'^3: Left-hand side (implicit w'^3 portion of the code).
      !
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
      = real( + one / dt )

      ! LHS mean advection (ma) term.
      lhs_a_csr((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/)) & 
      = lhs_a_csr((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/)) & 
      + term_ma_zt_lhs( wm_zt(k), gr%invrs_dzt(k), k, gr%invrs_dzm(k), gr%invrs_dzm(km1) )

      ! LHS turbulent advection (ta) term.
      if ( .not. l_explicit_turbulent_adv_wp3 ) then

         ! Note:  An "over-implicit" weighted time step is applied to this term.
         !        The weight of the implicit portion of this term is controlled
         !        by the factor gamma_over_implicit_ts (abbreviated "gamma" in
         !        the expression below).  A factor is added to the right-hand
         !        side of the equation in order to balance a weight that is not
         !        equal to 1, such that:
         !             -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
         !        where X is the variable that is being solved for in a
         !        predictive equation (w'^3 in this case), y(t) is the
         !        linearized portion of the term that gets treated implicitly,
         !        and RHS is the portion of the term that is always treated
         !        explicitly (in the case of the w'^3 turbulent advection term,
         !        RHS = 0).  A weight of greater than 1 can be applied to make
         !        the term more numerically stable.
         if ( iiPDF_type == iiPDF_ADG1 ) then

            ! The ADG1 PDF is used.
            wp3_term_ta_lhs_result(t_kp1_tdiag:t_km1_tdiag:-1) &
            = wp3_term_ta_ADG1_lhs( wp2(k), wp2(km1),  &
                                    a1(k), a1_zt(k), a1(km1),  &
                                    a3(k), a3_zt(k), a3(km1),  &
                                    wp3_on_wp2(k), wp3_on_wp2(km1), &
                                    rho_ds_zm(k), rho_ds_zm(km1),  &
                                    invrs_rho_ds_zt(k),  &
                                    gr%invrs_dzt(k), k )

            lhs_a_csr(t_kp1_tdiag:t_km1_tdiag:-1) & 
            = lhs_a_csr(t_kp1_tdiag:t_km1_tdiag:-1) &
              + gamma_over_implicit_ts * wp3_term_ta_lhs_result

         elseif ( iiPDF_type == iiPDF_new ) then

            ! The new PDF is used.
            wp3_term_ta_lhs_result((/t_k_mdiag,t_km1_mdiag/)) &
            = wp3_term_ta_new_pdf_lhs( coef_wp4_implicit(k), &
                                       coef_wp4_implicit(km1), &
                                       wp2(k), wp2(km1), rho_ds_zm(k), &
                                       rho_ds_zm(km1), invrs_rho_ds_zt(k), &
                                       gr%invrs_dzt(k) )

            lhs_a_csr((/t_k_mdiag,t_km1_mdiag/)) & 
            = lhs_a_csr((/t_k_mdiag,t_km1_mdiag/)) &
              + gamma_over_implicit_ts &
                * wp3_term_ta_lhs_result((/t_k_mdiag,t_km1_mdiag/))

         endif ! iiPDF_type

      else

         ! The turbulent advection term is being solved explicitly.
         wp3_term_ta_lhs_result(t_kp1_tdiag:t_km1_tdiag:-1) = zero

      endif ! .not. l_explicit_turbulent_adv_wp3

      ! LHS turbulent production (tp) term.
      ! Note:  An "over-implicit" weighted time step is applied to this term.
      wp3_term_tp_lhs_result((/t_k_mdiag,t_km1_mdiag/)) &
      = wp3_term_tp_lhs( wp2(k), wp2(km1), &
                         rho_ds_zm(k), rho_ds_zm(km1), &
                         invrs_rho_ds_zt(k), &
                         gr%invrs_dzt(k) )

      lhs_a_csr((/t_k_mdiag,t_km1_mdiag/)) & 
      = lhs_a_csr((/t_k_mdiag,t_km1_mdiag/)) & 
        + gamma_over_implicit_ts &
          * wp3_term_tp_lhs_result((/t_k_mdiag,t_km1_mdiag/))

      ! LHS pressure term 3 (pr3)
      if ( l_use_wp3_pr3 ) then

         wp3_terms_ta_tp_lhs_result &
         = wp3_term_ta_lhs_result + wp3_term_tp_lhs_result

         wp3_pr3_lhs(k,:) &
         = - gamma_over_implicit_ts * C16_fnc(k) * wp3_terms_ta_tp_lhs_result

         lhs(t_kp1_tdiag:t_km1_tdiag:-1)  &
         = lhs(t_kp1_tdiag:t_km1_tdiag:-1) + wp3_pr3_lhs(k,:)

      else

         wp3_pr3_lhs(k,:) = zero

      endif

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
        + C12 * one_half & 
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


      if (l_stats_samp) then

        ! Statistics: implicit contributions for wp3.

        if ( iwp3_ta > 0 ) then
          if ( .not. l_explicit_turbulent_adv_wp3 ) then
            ! Note:  An "over-implicit" weighted time step is applied to this
            !        term.  A weighting factor of greater than 1 may be used to
            !        make the term more numerically stable (see note above for
            !        LHS turbulent advection (ta) term).
            if ( iiPDF_type == iiPDF_ADG1 ) then
              tmp(1:5) &
              = gamma_over_implicit_ts &
                * wp3_term_ta_ADG1_lhs( wp2(k), wp2(km1), &
                                        a1(k), a1_zt(k), a1(km1), &
                                        a3(k), a3_zt(k), a3(km1), &
                                        wp3_on_wp2(k), wp3_on_wp2(km1), &
                                        rho_ds_zm(k), rho_ds_zm(km1), &
                                        invrs_rho_ds_zt(k), &
                                        gr%invrs_dzt(k), k )
              ztscr05(k) = -tmp(5)
              ztscr06(k) = -tmp(4)
              ztscr07(k) = -tmp(3)
              ztscr08(k) = -tmp(2)
              ztscr09(k) = -tmp(1)
            elseif ( iiPDF_type == iiPDF_new ) then
              tmp(1:2) &
              = gamma_over_implicit_ts &
                * wp3_term_ta_new_pdf_lhs( coef_wp4_implicit(k), &
                                           coef_wp4_implicit(km1), &
                                           wp2(k), wp2(km1), rho_ds_zm(k), &
                                           rho_ds_zm(km1), invrs_rho_ds_zt(k), &
                                           gr%invrs_dzt(k) )
              ztscr05(k) = zero
              ztscr06(k) = -tmp(2)
              ztscr07(k) = zero
              ztscr08(k) = -tmp(1)
              ztscr09(k) = zero
            endif ! iiPDF_type
          else
            ! The turbulent advection term is being solved explicitly.
            ztscr05(k) = zero
            ztscr06(k) = zero
            ztscr07(k) = zero
            ztscr08(k) = zero
            ztscr09(k) = zero
          endif ! .not. l_explicit_turbulent_adv_wp3
        endif ! iwp3_ta > 0

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note above for LHS turbulent
        !        advection (ta) term).
        if ( iwp3_tp > 0 ) then
          tmp(1:2)  &
          = gamma_over_implicit_ts  &
            * wp3_term_tp_lhs( wp2(k), wp2(km1), &
                               rho_ds_zm(k), rho_ds_zm(km1), &
                               invrs_rho_ds_zt(k), &
                               gr%invrs_dzt(k) )
          ztscr10(k) = -tmp(2)
          ztscr11(k) = -tmp(1)
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
          - wp3_terms_ac_pr2_lhs( (one+C11_Skw_fnc(k)), & 
                                  wm_zm(k), wm_zm(km1), gr%invrs_dzt(k) )
        endif

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note above for LHS turbulent
        !        advection (ta) term).
        if ( iwp3_pr1 > 0 ) then
          ztscr01(k)  &
          = - gamma_over_implicit_ts  &
            * wp3_term_pr1_lhs( C8, C8b, tauw3t(k), Skw_zt(k) )
        endif

        if ( iwp3_dp1 > 0 ) then
          if ( l_crank_nich_diff ) then
            ! Eddy diffusion for wp3 using a Crank-Nicholson time step.
            tmp(1:3) & 
            = C12 * one_half & 
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
    lhs_a_csr(wp2_cur_row + 1) = one

    ! w'^3
    lhs_a_csr(wp3_cur_row:wp3_cur_row + 2) = 0.0_core_rknd
    lhs_a_csr(wp3_cur_row) = one

    ! w'^2
    !lhs(:,k_wp2)         = 0.0_core_rknd
    !lhs(m_k_mdiag,k_wp2) = one
    ! w'^3
    !lhs(:,k_wp3)         = 0.0_core_rknd
    !lhs(t_k_tdiag,k_wp3) = one

    ! Upper boundary
    k = gr%nz
    k_wp3 = 2*k - 1
    k_wp2 = 2*k

    ! w'^2
    lhs_a_csr(intlc_5d_5d_ja_size - 2:intlc_5d_5d_ja_size) = 0.0_core_rknd
    lhs_a_csr(intlc_5d_5d_ja_size) = one

    ! w'^3
    lhs_a_csr(intlc_5d_5d_ja_size - 6:intlc_5d_5d_ja_size - 3) = 0.0_core_rknd
    lhs_a_csr(intlc_5d_5d_ja_size - 4) = one

    ! w'^2
    !lhs(:,k_wp2)         = 0.0_core_rknd
    !lhs(m_k_mdiag,k_wp2) = one
    ! w'^3
    !lhs(:,k_wp3)         = 0.0_core_rknd
    !lhs(t_k_tdiag,k_wp3) = one


    return
  end subroutine wp23_lhs_csr
#endif /* MKL */

  !=================================================================================
  subroutine wp23_rhs( dt, wp2, wp3, a1, a1_zt, a3, a3_zt, wp3_on_wp2, &
                       coef_wp4_implicit, wp4, wpthvp, wp2thvp, um, vm, & 
                       upwp, vpwp, up2, vp2, Kw1, Kw8, Kh_zt, & 
                       Skw_zt, tau1m, tauw3t, tau_C1_zm, C1_Skw_fnc, &
                       C11_Skw_fnc, C16_fnc, rho_ds_zm, invrs_rho_ds_zt, radf, &
                       thv_ds_zm, thv_ds_zt, wp2_splat, wp3_splat, & 
                       l_crank_nich_diff, &
                       rhs )

    ! Description:
    !   Compute RHS vector for w'^2 and w'^3.
    !   This subroutine computes the explicit portion of 
    !   the w'^2 and w'^3 equations.
    ! 
    !   Notes: 
    !        For LHS turbulent advection (ta) term.
    !           An "over-implicit" weighted time step is applied to this term.
    !           The weight of the implicit portion of this term is controlled
    !           by the factor gamma_over_implicit_ts (abbreviated "gamma" in
    !           the expression below).  A factor is added to the right-hand
    !           side of the equation in order to balance a weight that is not
    !           equal to 1, such that:
    !                -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
    !           where X is the variable that is being solved for in a
    !           predictive equation (w'^3 in this case), y(t) is the
    !           linearized portion of the term that gets treated implicitly,
    !           and RHS is the portion of the term that is always treated
    !           explicitly (in the case of the w'^3 turbulent advection term,
    !           RHS = 0).  A weight of greater than 1 can be applied to make
    !           the term more numerically stable.
    ! 
    ! 
    !  WARNING: This subroutine has been optimized. Significant changes could
    !           noticeably  impact computational efficiency. See clubb:ticket:834
    !-------------------------------------------------------------------------------

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
        one,           &
        one_half,      &
        zero,          &
        gamma_over_implicit_ts

    use model_flags, only:  & 
        l_tke_aniso,                  & ! Variable(s)
        l_explicit_turbulent_adv_wp3, &
        l_use_wp3_pr3

    use diffusion, only: & 
        diffusion_zm_lhs,  & ! Procedures
        diffusion_zm_lhs_all,  &
        diffusion_zt_lhs, &
        diffusion_zt_lhs_all

    use pdf_closure_module, only: &
        iiPDF_ADG1, & ! Variable(s)
        iiPDF_new,  &
        iiPDF_type

    use clubb_precision, only:  & 
        core_rknd ! Variable

    use stats_variables, only:  & 
        l_stats_samp, iwp2_dp1, iwp2_dp2, stats_zm, iwp2_bp,   & ! Variable(s)
        iwp2_pr1, iwp2_pr2, iwp2_pr3, iwp2_splat, iwp3_splat, &
        iwp3_ta, stats_zt, & 
        iwp3_tp, iwp3_bp1, iwp3_pr2, iwp3_pr1, iwp3_dp1, iwp3_bp2, iwp3_pr3
        
    use stats_type_utilities, only:  &
        stat_update_var_pt,  & ! Procedure(s)
        stat_begin_update_pt,  &
        stat_modify_pt

    use advance_helper_module, only: set_boundary_conditions_rhs


    implicit none

    ! Constant parameters
    logical, parameter :: &
      l_wp3_2nd_buoyancy_term = .true.

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      dt                 ! Timestep length                           [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      wp2,               & ! w'^2 (momentum levels)                    [m^2/s^2]
      wp3,               & ! w'^3 (thermodynamic levels)               [m^3/s^3]
      a1,                & ! sigma_sqd_w term a_1 (momentum levels)    [-]
      a1_zt,             & ! a_1 interpolated to thermodynamic levels  [-]
      a3,                & ! sigma_sqd_w term a_3 (momentum levels)    [-]
      a3_zt,             & ! a_3 interpolated to thermodynamic levels  [-]
      wp3_on_wp2,        & ! Smoothed version of wp3 / wp2             [m/s]
      coef_wp4_implicit, & ! <w'^4> = coef_wp4_implicit * <w'^2>^2     [-]
      wp4,               & ! w'^4 (momentum levels)                    [m^4/s^4]
      wpthvp,            & ! w'th_v' (momentum levels)                 [K m/s]
      wp2thvp,           & ! w'^2th_v' (thermodynamic levels)        [K m^2/s^2]
      um,                & ! u wind component (thermodynamic levels)   [m/s]
      vm,                & ! v wind component (thermodynamic levels)   [m/s]
      upwp,              & ! u'w' (momentum levels)                    [m^2/s^2]
      vpwp,              & ! v'w' (momentum levels)                    [m^2/s^2]
      up2,               & ! u'^2 (momentum levels)                    [m^2/s^2]
      vp2,               & ! v'^2 (momentum levels)                    [m^2/s^2]
      Kw1,               & ! Coefficient of eddy diffusivity for w'^2  [m^2/s]
      Kw8,               & ! Coefficient of eddy diffusivity for w'^3  [m^2/s]
      Kh_zt,             & ! Eddy diffusivity on thermodynamic levels  [m^2/s]
      Skw_zt,            & ! Skewness of w on thermodynamic levels     [-]
      tau1m,             & ! Time-scale tau on momentum levels         [s]
      tauw3t,            & ! Time-scale tau on thermodynamic levels    [s]
      tau_C1_zm,         & ! Tau values used for the C1 (dp1) term in wp2 [s]
      C1_Skw_fnc,        & ! C_1 parameter with Sk_w applied           [-]
      C11_Skw_fnc,       & ! C_11 parameter with Sk_w applied          [-]
      C16_fnc,           & ! C_16 parameter                            [-]
      rho_ds_zm,         & ! Dry, static density on momentum levels    [kg/m^3]
      invrs_rho_ds_zt,   & ! Inv. dry, static density @ thermo. levs.  [m^3/kg]
      radf,              & ! Buoyancy production at the CL top         [m^2/s^3]
      thv_ds_zm,         & ! Dry, base-state theta_v on momentum levs. [K]
      thv_ds_zt,         & ! Dry, base-state theta_v on thermo. levs.  [K]
      wp2_splat,         & ! Tendency of <w'^2> due to vertical compression of eddies [m^2/s^3]
      wp3_splat            ! Tendency of <w'^3> due to vertical compression of eddies [m^3/s^4]

    logical, intent(in) :: & 
      l_crank_nich_diff   ! Turns on/off Crank-Nicholson diffusion.

    ! Output Variable
    real( kind = core_rknd ), dimension(2*gr%nz), intent(out) :: & 
      rhs   ! RHS of band matrix

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: & 
      dum_dz, dvm_dz ! Vertical derivatives of um and vm

    ! Array indices
    integer :: k, k_wp2, k_wp3

    real( kind = core_rknd ), dimension(5,gr%nz) :: &
        wp3_term_ta_lhs_result

    real( kind = core_rknd ), dimension(3,gr%nz) :: &
        rhs_diff_zm, &
        rhs_diff_zt

    real( kind = core_rknd ), dimension(2,gr%nz) :: &
        lhs_tp_wp3, &
        lhs_ta_wp3

    real( kind = core_rknd ), dimension(gr%nz) :: &
        lhs_dp1_wp2, &          ! wp2 "over-implicit" dissipation term
        rhs_dp1_wp2, &          ! wp2 rhs dissipation term
        lhs_pr1_wp2, &          ! wp2 "over-implicit" pressure term 1
        rhs_pr1_wp2, &          ! wp2 rhs pressure term 1
        lhs_pr1_wp3, &          ! wp3 "over-implicit" pressure term 1
        rhs_pr1_wp3, &          ! wp3 rhs pressure term 1
        rhs_bp_pr2_wp2, &       ! wp2 bouyancy production and pressure term 2
        rhs_bp1_pr2_wp3, &      ! wp3 bouyancy production 1 and pressure term 2
        rhs_pr3_wp2, &          ! wp2 pressure term 3
        rhs_pr3_wp3, &          ! wp3 pressure term 3
        rhs_ta_wp3, &           ! wp3 turbulent advection term
        rhs_bp2_wp3             ! wp3 bouyancy production term 2 !--EXPERIMENTAL--!

    
    real( kind = core_rknd ) :: &
        invrs_dt        ! Inverse of dt, 1/dt, used for computational efficiency

    ! --------------- Begin Code ---------------
        

    ! Initialize arrays to 0 and calculate invers_dt
    invrs_dt = 1.0_core_rknd / dt
    rhs = 0.0_core_rknd
    wp3_term_ta_lhs_result = zero


    ! Experimental term from CLUBB TRAC ticket #411
    if ( l_wp3_2nd_buoyancy_term ) then

        ! Compute the vertical derivative of the u and v winds
          dum_dz = ddzt( um )
          dvm_dz = ddzt( vm )

        ! Calculate term
        call wp3_term_bp2_rhs_all( C15, Kh_zt(:), wpthvp(:), &
                                   dum_dz(:), dvm_dz(:), &
                                   upwp(:), vpwp(:), &
                                   thv_ds_zt(:), gr%invrs_dzt(:), &
                                   rhs_bp2_wp3(:) )
        ! Add term
        do k = 2, gr%nz-1

            k_wp3 = 2*k - 1

            rhs(k_wp3) = rhs(k_wp3) + rhs_bp2_wp3(k)

        end do

    end if



    ! These lines are for the diffusional term with a Crank-Nicholson
    ! time step.  They are not used for completely implicit diffusion.
    if ( l_crank_nich_diff ) then

        ! Calculate RHS eddy diffusion terms for w'2 and w'3
        
        call diffusion_zm_lhs_all( Kw1(:), nu1_vert_res_dep(:), & 
                                  gr%invrs_dzt(:), gr%invrs_dzm(:), &
                                  rhs_diff_zm(:,:) )

        call diffusion_zt_lhs_all( Kw8(:), nu8_vert_res_dep(:), & 
                                   gr%invrs_dzm(:), gr%invrs_dzt(:), &
                                   rhs_diff_zt(:,:) )
        ! Add diffusion terms
        do k = 2, gr%nz-1

            k_wp3 = 2*k - 1
            k_wp2 = 2*k

            rhs_diff_zm(1,k) = rhs_diff_zm(1,k) * one_half
            rhs_diff_zm(2,k) = rhs_diff_zm(2,k) * one_half
            rhs_diff_zm(3,k) = rhs_diff_zm(3,k) * one_half

            rhs_diff_zt(1,k) = rhs_diff_zt(1,k) * C12 * one_half
            rhs_diff_zt(2,k) = rhs_diff_zt(2,k) * C12 * one_half
            rhs_diff_zt(3,k) = rhs_diff_zt(3,k) * C12 * one_half
        
            rhs(k_wp2) = rhs(k_wp2) & 
                         - rhs_diff_zm(3,k) * wp2(k-1) & 
                         - rhs_diff_zm(2,k) * wp2(k) & 
                         - rhs_diff_zm(1,k) * wp2(k+1)

            rhs(k_wp3) = rhs(k_wp3) & 
                         - rhs_diff_zt(3,k) * wp3(k-1) & 
                         - rhs_diff_zt(2,k) * wp3(k) & 
                         - rhs_diff_zt(1,k) * wp3(k+1)
        end do

    endif
  

    if ( l_tke_aniso ) then

        ! Calculate "over-implicit" pressure terms for w'2 and w'3

        call wp2_term_pr1_rhs_all( C4, up2(:), vp2(:), tau1m(:), &
                                   rhs_pr1_wp2(:) )

        ! Note:  An "over-implicit" weighted time step is applied to the  term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note below for w'^3 RHS
        !        turbulent advection (ta) term).
        call wp2_term_pr1_lhs_all( C4, tau1m(:), &
                                   lhs_pr1_wp2(:) )

        ! Add pressure terms and splat terms
        do k = 2, gr%nz-1

            k_wp2 = 2*k

            rhs(k_wp2) = rhs(k_wp2) + rhs_pr1_wp2(k)

            rhs(k_wp2) = rhs(k_wp2) + ( one - gamma_over_implicit_ts ) &
                                    * ( - lhs_pr1_wp2(k) * wp2(k) )

            ! Effect of vertical compression of eddies
            rhs(k_wp2) = rhs(k_wp2) + wp2_splat(k)
        
        end do

    endif

    ! Calculate turbulent production terms of w'^3 
    call wp3_term_tp_lhs_all( wp2(:), &
                              rho_ds_zm(:), &
                              invrs_rho_ds_zt(:), &
                              gr%invrs_dzt(:), &
                              lhs_tp_wp3(:,:) )

    ! Calculate pressure terms 1 for w'^3
    call wp3_term_pr1_lhs_all( C8, C8b, tauw3t(:), Skw_zt(:), &
                               lhs_pr1_wp3(:) )

    ! Calculate dissipation terms 1 for w'^2
    call wp2_term_dp1_lhs_all( C1_Skw_fnc(:), tau_C1_zm(:), &
                               lhs_dp1_wp2(:) )

    ! Calculate buoyancy production of w'^2 and w'^2 pressure term 2
    call wp2_terms_bp_pr2_rhs_all( C5, thv_ds_zm(:), wpthvp(:), &
                                   rhs_bp_pr2_wp2(:) )

    ! Calculate pressure terms 3 for w'^2
    call wp2_term_pr3_rhs_all( C5, thv_ds_zm(:), wpthvp(:), upwp(:), &
                               um(:), vpwp(:), vm(:), gr%invrs_dzm(:), &
                               rhs_pr3_wp2(:) )

    ! Calculate dissipation terms 1 for w'^2
    call wp2_term_dp1_rhs_all( C1_Skw_fnc(:), tau_C1_zm(:), w_tol_sqd, up2(:), vp2(:), &
                               rhs_dp1_wp2(:) )

    ! Calculate buoyancy production of w'^3 and w'^3 pressure term 2
    call wp3_terms_bp1_pr2_rhs_all( C11_Skw_fnc(:), thv_ds_zt(:), wp2thvp(:), &
                                    rhs_bp1_pr2_wp3(:) )

    ! Calculate pressure terms 1 for w'^3
    call wp3_term_pr1_rhs_all( C8, C8b, tauw3t(:), Skw_zt(:), wp3(:), &
                               rhs_pr1_wp3(:) )


    ! Combine terms
    do k = 2, gr%nz-1

        k_wp3 = 2*k - 1
        k_wp2 = 2*k


        ! ------ Combine terms for 3rd moment of vertical velocity, <w'^3> ------ !

        ! RHS time tendency.
        rhs(k_wp3) = rhs(k_wp3) + invrs_dt * wp3(k)

        ! RHS contribution from "over-implicit" turbulent production (tp) term.
        rhs(k_wp3) = rhs(k_wp3) + ( one - gamma_over_implicit_ts )  &
                              * ( - lhs_tp_wp3(1,k) * wp2(k)  &
                                  - lhs_tp_wp3(2,k) * wp2(k-1) )

        ! RHS buoyancy production (bp) term and pressure term 2 (pr2).
        rhs(k_wp3) = rhs(k_wp3) + rhs_bp1_pr2_wp3(k)

        ! RHS term for vertical compression of eddies (w'^3 splat)
        rhs(k_wp3) = rhs(k_wp3) + wp3_splat(k) 

        ! RHS pressure term 1
        rhs(k_wp3) = rhs(k_wp3) + rhs_pr1_wp3(k)

        ! RHS "over implicit" pressure term 1 (pr1).
        rhs(k_wp3)  = rhs(k_wp3) + ( one - gamma_over_implicit_ts ) * ( - lhs_pr1_wp3(k) * wp3(k) )


        ! ------ Combine terms for 2nd moment of vertical velocity, <w'^2> ------ !

        ! RHS time tendency.
        rhs(k_wp2) = rhs(k_wp2) + invrs_dt * wp2(k)

        ! RHS buoyancy production (bp) term and pressure term 2 (pr2).
        rhs(k_wp2) = rhs(k_wp2) + rhs_bp_pr2_wp2(k)

        ! RHS buoyancy production at CL top due to LW radiative cooling
        rhs(k_wp2) = rhs(k_wp2) + radf(k) 

        ! RHS pressure term 3 (pr3).
        rhs(k_wp2) = rhs(k_wp2) + rhs_pr3_wp2(k)

        ! RHS dissipation term 1 (dp1).
        rhs(k_wp2) = rhs(k_wp2) + rhs_dp1_wp2(k)

        ! RHS "over implicit" pressure term 1 (pr1).
        rhs(k_wp2) = rhs(k_wp2) + ( one - gamma_over_implicit_ts ) * ( - lhs_dp1_wp2(k) * wp2(k) )

    enddo


    if ( l_explicit_turbulent_adv_wp3 ) then

        ! The turbulent advection term is being solved explicitly.

        call wp3_term_ta_explicit_rhs_all( wp4(:), &
                                           rho_ds_zm(:), &
                                           invrs_rho_ds_zt(:), &
                                           gr%invrs_dzt(:), &
                                           rhs_ta_wp3(:) )

        ! Add RHS turbulent advection (ta) terms
        do k = 2, gr%nz-1

            k_wp3 = 2*k - 1

            rhs(k_wp3) = rhs(k_wp3) + rhs_ta_wp3(k)

        end do

    else

        ! The turbulent advection term is being solved implicitly. See note above

        if ( iiPDF_type == iiPDF_ADG1 ) then

            ! The ADG1 PDF is used.

            ! Calculate terms
            call wp3_term_ta_ADG1_lhs_all( wp2(:), &
                                           a1(:), a1_zt(:), &
                                           a3(:), a3_zt(:), &
                                           wp3_on_wp2(:), &
                                           rho_ds_zm(:), &
                                           invrs_rho_ds_zt(:), &
                                           gr%invrs_dzt(:), &
                                           wp3_term_ta_lhs_result(:,:) )
            ! Add terms
            do k = 2, gr%nz-1

                k_wp3 = 2*k - 1

                rhs(k_wp3) = rhs(k_wp3) + ( one - gamma_over_implicit_ts ) &
                                        * ( - wp3_term_ta_lhs_result(1,k) * wp3(k+1) &
                                            - wp3_term_ta_lhs_result(2,k) * wp2(k) &
                                            - wp3_term_ta_lhs_result(3,k) * wp3(k) &
                                            - wp3_term_ta_lhs_result(4,k) * wp2(k-1) &
                                            - wp3_term_ta_lhs_result(5,k) * wp3(k-1) )
            end do

        elseif ( iiPDF_type == iiPDF_new ) then

            ! The new PDF is used.

            ! Calculate terms
            call wp3_term_ta_new_pdf_lhs_all( coef_wp4_implicit(:), &
                                                 wp2(:), rho_ds_zm(:), &
                                                 invrs_rho_ds_zt(:), &
                                                 gr%invrs_dzt(:), &
                                                 lhs_ta_wp3(:,:) )
            ! Add terms
            do k = 2, gr%nz-1

                k_wp3 = 2*k - 1

                wp3_term_ta_lhs_result(2,k) = lhs_ta_wp3(1,k)
                wp3_term_ta_lhs_result(4,k) = lhs_ta_wp3(2,k)

                rhs(k_wp3) = rhs(k_wp3) + ( one - gamma_over_implicit_ts ) &
                                        * ( - lhs_ta_wp3(1,k) * wp2(k) &
                                            - lhs_ta_wp3(2,k) * wp2(k-1) )
            end do

        endif ! iiPDF_type

    endif ! l_explicit_turbulent_adv_wp3



    if ( l_use_wp3_pr3 ) then

        ! Using pressure term 3 for w'3

        ! Calculate pressure term and add to rhs
        do k = 2, gr%nz-1

            k_wp3 = 2*k - 1

            rhs_pr3_wp3(k) = - ( one - gamma_over_implicit_ts ) * C16_fnc(k) &
                             * ( -   wp3_term_ta_lhs_result(1,k)                     * wp3(k+1) &
                                 - ( wp3_term_ta_lhs_result(2,k) + lhs_tp_wp3(1,k) ) * wp2(k) &
                                 -   wp3_term_ta_lhs_result(3,k)                     * wp3(k) &
                                 - ( wp3_term_ta_lhs_result(4,k) + lhs_tp_wp3(2,k) ) * wp2(k-1) &
                                 -   wp3_term_ta_lhs_result(5,k)                     * wp3(k-1) )

            rhs(k_wp3) = rhs(k_wp3) + rhs_pr3_wp3(k)

        end do

    else

        ! Not using pressure term, set to 0
        rhs_pr3_wp3 = zero

    endif

    ! --------- Boundary Conditions ---------

    ! Both wp2 and wp3 used fixed-point boundary conditions.
    ! Therefore, anything set in the above loop at both the upper
    ! and lower boundaries would be overwritten here.  However, the
    ! above loop does not extend to the boundary levels.  An array
    ! with a value of 1 at the main diagonal on the left-hand side
    ! and with values of 0 at all other diagonals on the left-hand
    ! side will preserve the right-hand side value at that level.

    ! The value of w'^2 at the lower boundary will remain the same.
    ! When the lower boundary is at the surface, the surface value of
    ! w'^2 is set in subroutine calc_surface_varnce (surface_varnce_module.F).

    ! The value of w'^3 at the lower boundary will be 0.
 
    ! The value of w'^2 at the upper boundary will be set to the threshold
    ! minimum value of w_tol_sqd.

    ! The value of w'^3 at the upper boundary will be set to 0.
    rhs(1) = 0.0_core_rknd
    rhs(2) = wp2(1)

    rhs(2*gr%nz-1) = 0.0_core_rknd
    rhs(2*gr%nz) = w_tol_sqd


    ! --------- Statistics output ---------
    if ( l_stats_samp ) then

        do k = 2, gr%nz-1

            ! ----------- w'2 -----------

            ! w'^2 term dp2 has both implicit and explicit components (if the
            ! Crank-Nicholson scheme is selected); call stat_begin_update_pt.  
            ! Since stat_begin_update_pt automatically subtracts the value sent in, 
            ! reverse the sign on right-hand side diffusion component.  If 
            ! Crank-Nicholson diffusion is not selected, the stat_begin_update_pt 
            ! will not be called.
            if ( l_crank_nich_diff ) then
              call stat_begin_update_pt( iwp2_dp2, k, & 
                rhs_diff_zm(3,k) * wp2(k-1) & 
              + rhs_diff_zm(2,k) * wp2(k) & 
              + rhs_diff_zm(1,k) * wp2(k+1), stats_zm )
            endif


            ! w'^2 term bp is completely explicit; call stat_update_var_pt.
            ! Note:  To find the contribution of w'^2 term bp, substitute 0 for the
            !        C_5 input to function wp2_terms_bp_pr2_rhs.
            call stat_update_var_pt( iwp2_bp, k, & 
              wp2_terms_bp_pr2_rhs( 0.0_core_rknd, thv_ds_zm(k), wpthvp(k) ), stats_zm )


            ! Include effect of vertical compression of eddies in wp2 budget
            call stat_update_var_pt( iwp2_splat, k, wp2_splat(k), stats_zm )


            if ( l_tke_aniso ) then

                ! w'^2 term pr1 has both implicit and explicit components; call
                ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
                ! subtracts the value sent in, reverse the sign on wp2_term_pr1_rhs.
                call stat_begin_update_pt( iwp2_pr1, k, -rhs_pr1_wp2(k), stats_zm )

                ! Note:  An "over-implicit" weighted time step is applied to this
                !        term.  A weighting factor of greater than 1 may be used to
                !        make the term more numerically stable (see note below for
                !        w'^3 RHS turbulent advection (ta) term).
                call stat_modify_pt( iwp2_pr1, k, &
                                   + ( one - gamma_over_implicit_ts )  &
                                   * ( - lhs_pr1_wp2(k) * wp2(k) ), stats_zm )
            endif

            ! w'^2 term pr2 has both implicit and explicit components; call
            ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
            ! subtracts the value sent in, reverse the sign on wp2_terms_bp_pr2_rhs.
            ! Note:  To find the contribution of w'^2 term pr2, add 1 to the
            !        C_5 input to function wp2_terms_bp_pr2_rhs.
            call stat_begin_update_pt( iwp2_pr2, k, & 
              -wp2_terms_bp_pr2_rhs( (one+C5), thv_ds_zm(k), wpthvp(k) ), stats_zm )

            ! w'^2 term dp1 has both implicit and explicit components; call
            ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
            ! subtracts the value sent in, reverse the sign on wp2_term_dp1_rhs.
            call stat_begin_update_pt( iwp2_dp1, k, -rhs_dp1_wp2(k), stats_zm )


            ! Note:  An "over-implicit" weighted time step is applied to this term.
            !        A weighting factor of greater than 1 may be used to make the
            !        term more numerically stable (see note below for w'^3 RHS
            !        turbulent advection (ta) term).
            call stat_modify_pt( iwp2_dp1, k, &
                                 + ( one - gamma_over_implicit_ts )  &
                                 * ( - lhs_dp1_wp2(k) * wp2(k) ), stats_zm )

            ! w'^2 term pr3 is completely explicit; call stat_update_var_pt.
            call stat_update_var_pt( iwp2_pr3, k, rhs_pr3_wp2(k), stats_zm )


            ! ----------- w'3 -----------

            if ( l_explicit_turbulent_adv_wp3 ) then !l_explicit_turbulent_adv_wp3

                ! The turbulent advection term is being solved explicitly.
                ! 
                ! The turbulent advection stats code is still set up in two parts,
                ! so call stat_begin_update_pt.  The implicit portion of the stat,
                ! which has a value of 0, will still be called later.  Since
                ! stat_begin_update_pt automatically subtracts the value sent in,
                ! reverse the sign on the input value.
                call stat_begin_update_pt( iwp3_ta, k, &
                                           -wp3_term_ta_explicit_rhs( wp4(k), wp4(k-1), &
                                                               rho_ds_zm(k), rho_ds_zm(k-1), &
                                                               invrs_rho_ds_zt(k), &
                                                               gr%invrs_dzt(k) ), &
                                           stats_zt )
            else

                ! The turbulent advection term is being solved implicitly.
                ! 
                ! Note:  An "over-implicit" weighted time step is applied to this
                !        term.  A weighting factor of greater than 1 may be used to
                !        make the term more numerically stable (see note above for
                !        RHS turbulent advection (ta) term).
                !        Call stat_begin_update_pt.  Since stat_begin_update_pt
                !        automatically subtracts the value sent in, reverse the sign
                !        on the input value.

                if ( iiPDF_type == iiPDF_ADG1 ) then

                    ! The ADG1 PDF is used.

                    call stat_begin_update_pt( iwp3_ta, k, &
                                                - ( one - gamma_over_implicit_ts )  &
                                                * ( - wp3_term_ta_lhs_result(1,k) * wp3(k+1)  &
                                                    - wp3_term_ta_lhs_result(2,k) * wp2(k)  &
                                                    - wp3_term_ta_lhs_result(3,k) * wp3(k)  &
                                                    - wp3_term_ta_lhs_result(4,k) * wp2(k-1)  &
                                                    - wp3_term_ta_lhs_result(5,k) * wp3(k-1) ), &
                                               stats_zt )

                elseif ( iiPDF_type == iiPDF_new ) then

                    ! The new PDF is used.

                    call stat_begin_update_pt( iwp3_ta, k, &
                                               - ( one - gamma_over_implicit_ts )  &
                                                 * ( - lhs_ta_wp3(1,k) * wp2(k)  &
                                                     - lhs_ta_wp3(2,k) * wp2(k-1) ), &
                                               stats_zt )
                endif

            endif

            ! Note:  An "over-implicit" weighted time step is applied to this term.
            !        A weighting factor of greater than 1 may be used to make the
            !        term more numerically stable (see note above for RHS turbulent
            !        production (tp) term).  Call stat_begin_update_pt.  Since
            !        stat_begin_update_pt automatically subtracts the value sent in,
            !        reverse the sign on the input value.
            call stat_begin_update_pt( iwp3_tp, k, &
                                       - ( one - gamma_over_implicit_ts )  &
                                         * ( - lhs_tp_wp3(1,k) * wp2(k)  &
                                             - lhs_tp_wp3(2,k) * wp2(k-1) ), &
                                       stats_zt )


            ! w'^3 pressure term 3 (pr3) explicit (rhs) contribution
            call stat_begin_update_pt( iwp3_pr3, k, rhs_pr3_wp3(k), stats_zt )


            ! w'^3 term bp is completely explicit; call stat_update_var_pt.
            ! Note:  To find the contribution of w'^3 term bp, substitute 0 for the
            !        C_11 skewness function input to function wp3_terms_bp1_pr2_rhs.
            call stat_update_var_pt( iwp3_bp1, k, & 
              wp3_terms_bp1_pr2_rhs( 0.0_core_rknd, thv_ds_zt(k), wp2thvp(k) ), stats_zt )


            ! w'^3 term pr2 has both implicit and explicit components; call
            ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
            ! subtracts the value sent in, reverse the sign on wp3_terms_bp1_pr2_rhs.
            ! Note:  To find the contribution of w'^3 term pr2, add 1 to the
            !        C_11 skewness function input to function wp3_terms_bp1_pr2_rhs.
            call stat_begin_update_pt( iwp3_pr2, k, & 
                                       -wp3_terms_bp1_pr2_rhs( (one+C11_Skw_fnc(k)), &
                                       thv_ds_zt(k), wp2thvp(k) ), & 
                                       stats_zt )

            ! w'^3 term pr1 has both implicit and explicit components; call 
            ! stat_begin_update_pt.  Since stat_begin_update_pt automatically 
            ! subtracts the value sent in, reverse the sign on wp3_term_pr1_rhs.
            call stat_begin_update_pt( iwp3_pr1, k, -rhs_pr1_wp3(k), stats_zt )


            ! Note:  An "over-implicit" weighted time step is applied to this term.
            !        A weighting factor of greater than 1 may be used to make the
            !        term more numerically stable (see note above for RHS turbulent
            !        advection (ta) term).
            call stat_modify_pt( iwp3_pr1, k,  &
                                 + ( one - gamma_over_implicit_ts )  &
                                 * ( - lhs_pr1_wp3(k) * wp3(k) ), stats_zt )


            ! Include effect of vertical compression of eddies in wp2 budget
            call stat_update_var_pt( iwp3_splat, k, wp3_splat(k), stats_zt )


            if ( l_crank_nich_diff ) then

                ! w'^3 term dp1 has both implicit and explicit components (if the
                ! Crank-Nicholson scheme is selected); call stat_begin_update_pt.  
                ! Since stat_begin_update_pt automatically subtracts the value sent in, 
                ! reverse the sign on right-hand side diffusion component.  If 
                ! Crank-Nicholson diffusion is not selected, the stat_begin_update_pt 
                ! will not be called.
                call stat_begin_update_pt( iwp3_dp1, k, & 
                                           rhs_diff_zt(3,k) * wp3(k-1) & 
                                         + rhs_diff_zt(2,k) * wp3(k) & 
                                         + rhs_diff_zt(1,k) * wp3(k+1), stats_zt )
            endif
                      
            ! Experimental bouyancy term
            if ( l_wp3_2nd_buoyancy_term ) then
                call stat_update_var_pt( iwp3_bp2, k, rhs_bp2_wp3(k), stats_zt )
            end if

        end do

    endif

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
    pure subroutine wp2_term_ta_lhs_all( rho_ds_zt, &
                                         invrs_rho_ds_zm, &
                                         invrs_dzm, &
                                         lhs_ta_wp2 )
    ! Description:
    !     This subroutine serves the same function as wp2_term_ta_lhs (above), but
    !     calculates terms for all grid levels at once rather than one at a time.
    !     This was done so that this code could be vectorized and thereby sped up
    !     by the compiler. See clubb:ticket:834 for more information.
    ! 
    !-----------------------------------------------------------------------------

        use clubb_precision, only: &
          core_rknd ! Variable(s)

        use grid_class, only:  & 
            gr ! Variable

        implicit none


        ! Input Variables
        real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
          rho_ds_zt,       & ! Dry, static density at thermo. level (k)    [kg/m^3]
          invrs_rho_ds_zm, & ! Inv. dry, static density @ moment. lev. (k) [m^3/kg]
          invrs_dzm          ! Inverse of grid spacing (k)                 [1/m]

        ! Return Variable
        real( kind = core_rknd ), dimension(2,gr%nz), intent(out) :: lhs_ta_wp2

        ! Loop variable
        integer :: k

        ! Set lower boundary to 0
        lhs_ta_wp2(1,1) = 0.0_core_rknd
        lhs_ta_wp2(2,1) = 0.0_core_rknd

        ! Calculate non-boundary terms
        do k = 2, gr%nz-1 

            ! Thermodynamic superdiagonal: [ x wp3(k+1,<t+1>) ]
            lhs_ta_wp2(1,k) = + invrs_rho_ds_zm(k) * invrs_dzm(k) * rho_ds_zt(k+1)

            ! Thermodynamic subdiagonal: [ x wp3(k,<t+1>) ]
            lhs_ta_wp2(2,k) = - invrs_rho_ds_zm(k) * invrs_dzm(k) * rho_ds_zt(k)

        end do

        ! Set upper boundary to 0
        lhs_ta_wp2(1,gr%nz) = 0.0_core_rknd
        lhs_ta_wp2(2,gr%nz) = 0.0_core_rknd

        return

    end subroutine wp2_term_ta_lhs_all

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

    use constants_clubb, only: &
        two, & ! Variable(s)
        one

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
    = + ( one - C5 ) * two * invrs_dzm * ( wm_ztp1 - wm_zt )

    return

  end function wp2_terms_ac_pr2_lhs

    !==================================================================================
    pure subroutine wp2_terms_ac_pr2_lhs_all( C5, wm_zt, invrs_dzm, &
                                              lhs_ac_pr2_wp2 )
    ! Description:
    !     This subroutine serves the same function as wp2_terms_ac_pr2_lhs (above), but
    !     calculates terms for all grid levels at once rather than one at a time.
    !     This was done so that this code could be vectorized and thereby sped up
    !     by the compiler. See clubb:ticket:834 for more information.
    ! 
    !----------------------------------------------------------------------------------

        use clubb_precision, only: &
            core_rknd ! Variable(s)

        use grid_class, only:  & 
            gr ! Variable

        implicit none

        ! Input Variables
        real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
          wm_zt,     & ! w wind component at thermodynamic levels (k)   [m/s]
          invrs_dzm    ! Inverse of grid spacing (k)                    [1/m]

        real( kind = core_rknd ), intent(in) :: & 
          C5

        ! Return Variable
        real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
            lhs_ac_pr2_wp2

        integer :: k

        ! Set lower boundary to 0
        lhs_ac_pr2_wp2(1) = 0.0_core_rknd

        ! Calculate non-boundary values
        do k = 2, gr%nz-1

            ! Momentum main diagonal: [ x wp2(k,<t+1>) ]
            lhs_ac_pr2_wp2(k) = + ( 1.0_core_rknd - C5 ) * 2.0_core_rknd &
                                * invrs_dzm(k) * ( wm_zt(k+1) - wm_zt(k) )

        end do

        ! Set upper boundary to 0
        lhs_ac_pr2_wp2(gr%nz) = 0.0_core_rknd

        return

    end subroutine wp2_terms_ac_pr2_lhs_all

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

    !==================================================================================
    pure subroutine wp2_term_dp1_lhs_all( C1_Skw_fnc, tau1m, &
                                          lhs_dp1_wp2 )
    ! Description:
    !     This subroutine serves the same function as wp2_term_dp1_lhs (above), but
    !     calculates terms for all grid levels at once rather than one at a time.
    !     This was done so that this code could be vectorized and thereby sped up
    !     by the compiler. See clubb:ticket:834 for more information.
    ! 
    !----------------------------------------------------------------------------------

        use clubb_precision, only: &
            core_rknd ! Variable(s)

        use grid_class, only:  & 
            gr ! Variable

        implicit none

        ! Input Variables
        real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
          C1_Skw_fnc,  & ! C_1 parameter with Sk_w applied (k)   [-]
          tau1m          ! Time-scale tau at momentum levels (k) [s]

        ! Return Variable
        real( kind = core_rknd ), dimension(gr%nz), intent(out) :: lhs_dp1_wp2

        ! Loop variable
        integer :: k

        ! Set lower boundary to 0
        lhs_dp1_wp2(1) = 0.0_core_rknd

        ! Calculate non-boundary values
        do k = 2, gr%nz-1

            ! Momentum main diagonal: [ x wp2(k,<t+1>) ]
            lhs_dp1_wp2(k) = + C1_Skw_fnc(k) / tau1m(k)

        end do

        ! Set upper boundary to 0
        lhs_dp1_wp2(gr%nz) = 0.0_core_rknd

        return
    end subroutine wp2_term_dp1_lhs_all

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

    use constants_clubb, only: &
        three, & ! Variable(s)
        two

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
    = + ( two * C4 ) / ( three * tau1m )

    return
  end function wp2_term_pr1_lhs

    !==================================================================================
    pure subroutine wp2_term_pr1_lhs_all( C4, tau1m, &
                                          lhs_pr1_wp2 )
    ! Description:
    !     This subroutine serves the same function as wp2_term_pr1_lhs (above), but
    !     calculates terms for all grid levels at once rather than one at a time.
    !     This was done so that this code could be vectorized and thereby sped up
    !     by the compiler. See clubb:ticket:834 for more information.
    ! 
    !----------------------------------------------------------------------------------

        use constants_clubb, only: &
            three, & ! Variable(s)
            two

        use clubb_precision, only: &
            core_rknd ! Variable(s)

        use grid_class, only:  &
            gr      ! Variable 

        implicit none

        ! Input Variables
        real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
          tau1m   ! Time-scale tau at momentum levels (k) [s]

        real( kind = core_rknd ), intent(in) :: & 
          C4      ! Model parameter C_4                   [-]

        ! Return Variable
        real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
            lhs_pr1_wp2
    
        ! Loop variable
        integer :: k

        ! Set lower boundary to 0
        lhs_pr1_wp2(1) = 0.0_core_rknd

        ! Calculate non-boundary values
        do k = 2, gr%nz-1

            ! Momentum main diagonal: [ x wp2(k,<t+1>) ]
            lhs_pr1_wp2(k) = + ( two * C4 ) / ( three * tau1m(k) )
    
        end do

        ! Set upper boundary to 0
        lhs_pr1_wp2(gr%nz) = 0.0_core_rknd

        return

    end subroutine wp2_term_pr1_lhs_all

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

    use constants_clubb, only:  & ! Variable(s)        
        grav, & ! Gravitational acceleration [m/s^2]
        two,  &
        one

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C5,        & ! Model parameter C_5                             [-]
      thv_ds_zm, & ! Dry, base-state theta_v at momentum level (k)   [K]
      wpthvp       ! w'th_v'(k)                                      [K m/s]

    ! Return Variable
    real( kind = core_rknd ) :: rhs

    rhs & 
    = + ( one - C5 ) * two * ( grav / thv_ds_zm ) * wpthvp

    return
  end function wp2_terms_bp_pr2_rhs

    !==================================================================================
    pure subroutine wp2_terms_bp_pr2_rhs_all( C5, thv_ds_zm, wpthvp, &
                                              rhs_bp_pr2_wp2 )
    ! Description:
    !     This subroutine serves the same function as wp2_terms_bp_pr2_rhs (above), but
    !     calculates terms for all grid levels at once rather than one at a time.
    !     This was done so that this code could be vectorized and thereby sped up
    !     by the compiler. See clubb:ticket:834 for more information.
    ! 
    !----------------------------------------------------------------------------------

        use clubb_precision, only: &
            core_rknd ! Variable(s)

        use grid_class, only: &
            gr

        use constants_clubb, only:  & ! Variable(s)        
            grav, & ! Gravitational acceleration [m/s^2]
            two,  &
            one

        implicit none

        ! Input Variables
        real( kind = core_rknd ), intent(in) :: & 
          C5           ! Model parameter C_5                             [-]

        real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
          thv_ds_zm, & ! Dry, base-state theta_v at momentum level (k)   [K]
          wpthvp       ! w'th_v'(k)                                      [K m/s]

        ! Return Variable
        real( kind = core_rknd ), dimension(gr%nz), intent(out) :: rhs_bp_pr2_wp2

        ! Loop variable
        integer :: k

        ! Set lower boundary to 0
        rhs_bp_pr2_wp2(1) = 0.0_core_rknd

        ! Calculate non-boundary values
        do k = 2, gr%nz-1
            rhs_bp_pr2_wp2(k) = + ( one - C5 ) * two * ( grav / thv_ds_zm(k) ) * wpthvp(k)
        end do

        ! Set upper boundary to 0
        rhs_bp_pr2_wp2(gr%nz) = 0.0_core_rknd

        return
    end subroutine wp2_terms_bp_pr2_rhs_all

  !=============================================================================
  pure function wp2_term_dp1_rhs( C1_Skw_fnc, tau1m, threshold, up2, vp2 ) & 
  result( rhs )

    ! Description:
    ! When l_damp_wp2_using_em == .false., then
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

    ! if l_damp_wp2_using_em == .true., then
    ! we damp wp2 using a more standard turbulence closure, -(2/3)*em/tau
    ! This only works if C1=C14 and l_stability_correct_tau_zm =.false.
    ! A factor of (1/3) is absorbed into C1.
    ! The threshold is implicitly set to 0.


    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use model_flags, only: &
        l_damp_wp2_using_em ! Logical

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C1_Skw_fnc,  & ! C_1 parameter with Sk_w applied (k)   [-]
      tau1m,       & ! Time-scale tau at momentum levels (k) [s]
      threshold,   & ! Minimum allowable value of w'^2       [m^2/s^2]
      up2,         & ! Horizontal (east-west) velocity variance, u'^2 [m^2/s^2]
      vp2            ! Horizontal (north-south) velocity variance, v'^2 [m^2/s^2]

    ! Return Variable
    real( kind = core_rknd ) :: rhs


    if ( l_damp_wp2_using_em ) then

      rhs & 
      = - ( C1_Skw_fnc / tau1m ) * ( up2 + vp2 )

    else

      rhs & 
      = + ( C1_Skw_fnc / tau1m ) * threshold

    end if

    return
  end function wp2_term_dp1_rhs

    !==================================================================================
    pure subroutine wp2_term_dp1_rhs_all( C1_Skw_fnc, tau1m, threshold, up2, vp2, &
                                          rhs_dp1_wp2 )
    ! Description:
    !     This subroutine serves the same function as wp2_term_dp1_rhs (above), but
    !     calculates terms for all grid levels at once rather than one at a time.
    !     This was done so that this code could be vectorized and thereby sped up
    !     by the compiler. See clubb:ticket:834 for more information.
    ! 
    !----------------------------------------------------------------------------------

        use clubb_precision, only: &
            core_rknd ! Variable(s)

        use model_flags, only: &
            l_damp_wp2_using_em ! Logical

        use grid_class, only: &
            gr

        implicit none

        ! Input Variables
        real( kind = core_rknd ), intent(in) :: & 
          threshold      ! Minimum allowable value of w'^2       [m^2/s^2]

        real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
          C1_Skw_fnc,  & ! C_1 parameter with Sk_w applied (k)   [-]
          tau1m,       & ! Time-scale tau at momentum levels (k) [s]
          up2,         & ! Horizontal (east-west) velocity variance, u'^2 [m^2/s^2]
          vp2            ! Horizontal (north-south) velocity variance, v'^2 [m^2/s^2]

        ! Return Variable
        real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
            rhs_dp1_wp2

        ! Loop variable
        integer :: k

        ! Set lower boundary to 0
        rhs_dp1_wp2(1) = 0.0_core_rknd

        ! Calculate non-boundary values
        if ( l_damp_wp2_using_em ) then

            do k = 2, gr%nz-1
                rhs_dp1_wp2(k) = - ( C1_Skw_fnc(k) / tau1m(k) ) * ( up2(k) + vp2(k) )
            end do

        else

            do k = 2, gr%nz-1
                rhs_dp1_wp2(k) = + ( C1_Skw_fnc(k) / tau1m(k) ) * threshold
            end do

        end if

        ! Set upper boundary to 0
        rhs_dp1_wp2(gr%nz) = 0.0_core_rknd

        return
    end subroutine wp2_term_dp1_rhs_all

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
        grav,       & ! Gravitational acceleration [m/s^2]
        two_thirds, &
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
    = + two_thirds * C5 & 
                   * ( ( grav / thv_ds_zm ) * wpthvp & 
                       - upwp * invrs_dzm * ( ump1 - um ) & 
                       - vpwp * invrs_dzm * ( vmp1 - vm ) & 
                     )
     ! Use the following code for alternate mixing, with c_k=0.1 or 0.2
!    = + two_thirds * C5 &
!                   * ( ( grav / thv_ds_zm ) * wpthvp &
!                       - 0. * upwp * invrs_dzm * ( ump1 - um ) &
!                       - 0. * vpwp * invrs_dzm * ( vmp1 - vm ) &
!                     )
!    eMFc

    ! Added by dschanen for ticket #36
    ! We have found that when shear generation is zero this term will only be
    ! offset by hole-filling (wp2_pd) and reduces turbulence 
    ! unrealistically at lower altitudes to make up the difference.
    rhs = max( rhs, zero_threshold )

    return
  end function wp2_term_pr3_rhs

    !==================================================================================
    pure subroutine wp2_term_pr3_rhs_all( C5, thv_ds_zm, wpthvp, upwp, &
                                          um, vpwp, vm, invrs_dzm, &
                                          rhs_pr3_wp2 )
    ! Description:
    !     This subroutine serves the same function as wp2_term_pr3_rhs (above), but
    !     calculates terms for all grid levels at once rather than one at a time.
    !     This was done so that this code could be vectorized and thereby sped up
    !     by the compiler. See clubb:ticket:834 for more information.
    ! 
    !--------------------------------------------------------------------------------------

        use clubb_precision, only: &
            core_rknd ! Variable(s)

        use constants_clubb, only: & ! Variables 
            grav,       & ! Gravitational acceleration [m/s^2]
            two_thirds, &
            zero_threshold

        use grid_class, only: &
            gr

        implicit none

        
        ! Input Variables
        real( kind = core_rknd ), intent(in) :: & 
          C5           ! Model parameter C_5                            [-]

        real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
          thv_ds_zm, & ! Dry, base-state theta_v at momentum level (k)  [K]
          wpthvp,    & ! w'th_v'(k)                                     [K m/s]
          upwp,      & ! u'w'(k)                                        [m^2/s^2]
          um,        & ! um(k)                                          [m/s]
          vpwp,      & ! v'w'(k)                                        [m^2/s^2]
          vm,        & ! vm(k)                                          [m/s]
          invrs_dzm    ! Inverse of grid spacing (k)                    [1/m]

        ! Return Variable
        real( kind = core_rknd ), dimension(gr%nz), intent(out) :: rhs_pr3_wp2

        ! Loop variable
        integer :: k

        ! Set lower boundary to 0
        rhs_pr3_wp2(1) = 0.0_core_rknd

        ! Calculate non-boundary value
        do k = 2, gr%nz-1

            rhs_pr3_wp2(k) = + two_thirds * C5 & 
                             * ( ( grav / thv_ds_zm(k) ) * wpthvp(k) & 
                                 - upwp(k) * invrs_dzm(k) * ( um(k+1) - um(k) ) & 
                                 - vpwp(k) * invrs_dzm(k) * ( vm(k+1) - vm(k) ) )
        end do

        ! Set upper boundary to 0
        rhs_pr3_wp2(gr%nz) = 0.0_core_rknd

        ! Added by dschanen for ticket #36
        ! We have found that when shear generation is zero this term will only be
        ! offset by hole-filling (wp2_pd) and reduces turbulence 
        ! unrealistically at lower altitudes to make up the difference.
        rhs_pr3_wp2 = max( rhs_pr3_wp2, zero_threshold )

        return
    end subroutine wp2_term_pr3_rhs_all

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

    use constants_clubb, only: &
        three    ! Variable(s)

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
    = + ( C4 * ( up2 + vp2 ) ) / ( three * tau1m )

    return
  end function wp2_term_pr1_rhs

    !==================================================================================
    pure subroutine wp2_term_pr1_rhs_all( C4, up2, vp2, tau1m, &
                                          rhs_pr1_wp2 )
    ! Description:
    !     This subroutine serves the same function as wp2_term_pr1_rhs (above), but
    !     calculates terms for all grid levels at once rather than one at a time.
    !     This was done so that this code could be vectorized and thereby sped up
    !     by the compiler. See clubb:ticket:834 for more information.
    ! 
    !--------------------------------------------------------------------------------------

        use constants_clubb, only: &
            three    ! Variable(s)

        use clubb_precision, only: &
            core_rknd ! Variable(s)

        use grid_class, only: &
            gr

        implicit none

        ! Input Variables
        real( kind = core_rknd ), intent(in) :: & 
          C4     ! Model parameter C_4                   [-]

        real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
          up2,  & ! u'^2(k)                               [m^2/s^2]
          vp2,  & ! v'^2(k)                               [m^2/s^2]
          tau1m   ! Time-scale tau at momentum levels (k) [s]

        ! Return Variable
        real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
            rhs_pr1_wp2

        ! Loop Variable
        integer :: k

        ! Set lower bounadry to 0
        rhs_pr1_wp2(1) = 0.0_core_rknd

        ! Calculate non-boundary values
        do k = 2, gr%nz-1

            rhs_pr1_wp2(k) = + ( C4 * ( up2(k) + vp2(k) ) ) / ( three * tau1m(k) )

        end do

        ! Set upper boundary to 0
        rhs_pr1_wp2(gr%nz) = 0.0_core_rknd

        return
    end subroutine wp2_term_pr1_rhs_all

  !=============================================================================
  pure function wp3_term_ta_new_pdf_lhs( coef_wp4_implicit, &
                                         coef_wp4_implicitm1, &
                                         wp2, wp2m1, rho_ds_zm, &
                                         rho_ds_zmm1, invrs_rho_ds_zt, &
                                         invrs_dzt ) &
  result( lhs )

    ! Description:
    ! Turbulent advection of <w'^3>:  implicit portion of the code.
    !
    ! This implicit discretization is specifically for the new PDF.
    !
    ! The d<w'^3>/dt equation contains a turbulent advection term:
    !
    ! - (1/rho_ds) * d( rho_ds * <w'^4> )/dz.
    !
    ! A substitution, which is specific to the new PDF, is made in order to
    ! close the turbulent advection term, such that:
    !
    ! <w'^4> = coef_wp4_implicit * <w'^2>^2.
    !
    ! The calculation of coef_wp4_implicit is detailed in function
    ! calc_coef_wp4_implicit, which is found in module new_pdf in new_pdf.F90.
    !
    ! The turbulent advection term is rewritten as:
    !
    ! - (1/rho_ds) * d( rho_ds * coef_wp4_implicit * <w'^2>^2 )/dz.
    !
    ! The <w'^2>^2 term is timestep split so that it can be expressed linearly
    ! in terms of <w'^2> at the (t+1) timestep, such that:
    !
    ! <w'^2>^2 = <w'^2>(t) * <w'^2>(t+1);
    !
    ! which allows the turbulent advection term to be expressed implicitly as:
    !
    ! - (1/rho_ds)
    !   * d( rho_ds * coef_wp4_implicit * <w'^2>(t) * <w'^2>(t+1) )/dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign is
    !        reversed and the leading "-" in front of all d[ ] / dz terms is
    !        changed to a "+".
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is
    ! being advanced to in solving the d<w'^3>/dt and d<w'^2>/dt equations.
    !
    ! The implicit discretization of this term is as follows:
    !
    ! The values of <w'^3> are found on the thermodynamic levels, while the
    ! values of <w'^2> are found on the momentum levels.  The values of
    ! coef_wp4_implicit_zt are originally calculated by the PDF on the
    ! thermodynamic levels.  They are interpolated to the intermediate momentum
    ! levels as coef_wp4_implicit.  Additionally, the values of rho_ds_zm are
    ! found on the momentum levels, and the values of invrs_rho_ds_zt are found
    ! on the thermodynamic levels.  At the intermediate momentum levels, the
    ! values of coef_wp4_implicit are multiplied by <w'^2>(t) * <w'^2>(t+1), and
    ! the resulting product is also multiplied by rho_ds_zm.  This product is
    ! referred to as G below.  Then, the derivative (d/dz) of that expression is
    ! taken over the central thermodynamic level, where it is multiplied by
    ! -invrs_rho_ds_zt.  This yields the desired result.  In this function,
    ! the values of G are as follows:
    !
    ! G = rho_ds_zm * coef_wp4_implicit * <w'^2>(t) * <w'^2>(t+1).
    !
    ! -------coef_wp4_implicit_zt---------------------------------------- t(k+1)
    !
    ! =======coef_wp4_implicit(interp)=======wp2=========rho_ds_zm======= m(k)
    !
    ! -------coef_wp4_implicit_zt-----dG/dz-----invrs_rho_ds_zt----wp3--- t(k)
    !
    ! =======coef_wp4_implicitm1(interp)=====wp2m1=======rho_ds_zmm1===== m(k-1)
    !
    ! -------coef_wp4_implicit_zt---------------------------------------- t(k-1)
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

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      k_mdiag   = 1, & ! Momentum superdiagonal index.
      km1_mdiag = 2    ! Momentum subdiagonal index. 

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      coef_wp4_implicit,   & ! <w'^4>=coef_wp4_implicit*<w'^2>^2; m-lev(k)   [-]
      coef_wp4_implicitm1, & ! <w'^4>=coef_wp4_implicit*<w'^2>^2; m-lev(k-1) [-]
      wp2,                 & ! w'^2(k)                                 [m^2/s^2]
      wp2m1,               & ! w'^2(k-1)                               [m^2/s^2]
      rho_ds_zm,           & ! Dry, static density at mom lev (k)       [kg/m^3]
      rho_ds_zmm1,         & ! Dry, static density at mom lev (k-1)     [kg/m^3]
      invrs_rho_ds_zt,     & ! Inv dry, static density @ thermo lev (k) [m^3/kg]
      invrs_dzt              ! Inverse of grid spacing (k)              [1/m]

    ! Return Variable
    real( kind = core_rknd ), dimension(2) :: lhs


    ! Momentum superdiagonal: [ x wp2(k,<t+1>) ]
    lhs(k_mdiag) &
    = + invrs_rho_ds_zt * invrs_dzt * rho_ds_zm * coef_wp4_implicit * wp2

    ! Momentum subdiagonal: [ x wp2(k-1,<t+1>) ]
    lhs(km1_mdiag) &
    = - invrs_rho_ds_zt * invrs_dzt * rho_ds_zmm1 * coef_wp4_implicitm1 * wp2m1


    return

  end function wp3_term_ta_new_pdf_lhs

    !======================================================================================
    pure subroutine wp3_term_ta_new_pdf_lhs_all( coef_wp4_implicit, &
                                                 wp2, rho_ds_zm, &
                                                 invrs_rho_ds_zt, &
                                                 invrs_dzt, &
                                                 lhs_ta_wp3 )
    ! Description:
    !     This subroutine serves the same function as wp3_term_ta_new_pdf_lhs (above), but
    !     calculates terms for all grid levels at once rather than one at a time.
    !     This was done so that this code could be vectorized and thereby sped up
    !     by the compiler. See clubb:ticket:834 for more information.
    ! 
    !--------------------------------------------------------------------------------------

        use clubb_precision, only: &
            core_rknd ! Variable(s)

        use grid_class, only:  &
            gr      ! Variable

        implicit none

        ! Constant parameters
        integer, parameter :: & 
          k_mdiag   = 1, & ! Momentum superdiagonal index.
          km1_mdiag = 2    ! Momentum subdiagonal index. 

        ! Input Variables
        real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
          coef_wp4_implicit,   & ! <w'^4>=coef_wp4_implicit*<w'^2>^2; m-lev(k)   [-]
          wp2,                 & ! w'^2(k)                                 [m^2/s^2]
          rho_ds_zm,           & ! Dry, static density at mom lev (k)       [kg/m^3]
          invrs_rho_ds_zt,     & ! Inv dry, static density @ thermo lev (k) [m^3/kg]
          invrs_dzt              ! Inverse of grid spacing (k)              [1/m]

        ! Return Variable
        real( kind = core_rknd ), dimension(2,gr%nz), intent(out) :: lhs_ta_wp3

        ! Loop variable
        integer :: k

        ! Set lower boundary to 0
        lhs_ta_wp3(1,1) = 0.0_core_rknd
        lhs_ta_wp3(2,1) = 0.0_core_rknd

        ! Calculate non-boundary terms
        do k = 2, gr%nz-1

            ! Momentum superdiagonal: [ x wp2(k,<t+1>) ]
            lhs_ta_wp3(1,k) = + invrs_rho_ds_zt(k) * invrs_dzt(k) * rho_ds_zm(k) &
                              * coef_wp4_implicit(k) * wp2(k)

            ! Momentum subdiagonal: [ x wp2(k-1,<t+1>) ]
            lhs_ta_wp3(2,k) = - invrs_rho_ds_zt(k) * invrs_dzt(k) * rho_ds_zm(k-1) &
                              * coef_wp4_implicit(k-1) * wp2(k-1)

        end do

        ! Set upper boundary to 0
        lhs_ta_wp3(1,gr%nz) = 0.0_core_rknd
        lhs_ta_wp3(2,gr%nz) = 0.0_core_rknd

        return

    end subroutine wp3_term_ta_new_pdf_lhs_all

  !=============================================================================
  pure function wp3_term_ta_ADG1_lhs( wp2, wp2m1, &
                                      a1, a1_zt, a1m1, &
                                      a3, a3_zt, a3m1, &
                                      wp3_on_wp2, wp3_on_wp2_m1, &
                                      rho_ds_zm, rho_ds_zmm1, &
                                      invrs_rho_ds_zt, &
                                      invrs_dzt, level ) &
  result( lhs )

    ! Description:
    ! Turbulent advection of w'^3:  implicit portion of the code.
    !
    ! This implicit discretization is specifically for the ADG1 PDF.
    !
    ! The d(w'^3)/dt equation contains a turbulent advection term:
    !
    ! - (1/rho_ds) * d( rho_ds * w'^4 )/dz.
    !
    ! A substitution, which is specific to ADG1, is made in order to close the
    ! turbulent advection term, such that:
    !
    ! w'^4 = a_3 * (w'^2)^2  +  a_1 * ( (w'^3)^2 / w'^2 );
    !
    ! where both a_1 and a_3 are variables that are functions of sigma_sqd_w,
    ! such that: 
    !
    ! a_1 = 1 / (1 - sigma_sqd_w); and
    !
    ! a_3 = 3*(sigma_sqd_w)^2 + 6*(1 - sigma_sqd_w)*sigma_sqd_w
    !       + (1 - sigma_sqd_w)^2.
    !
    ! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:wp4_diagnosis
    !
    ! The turbulent advection term is rewritten as:
    !
    ! - (1/rho_ds) * d [ rho_ds * a_3 * (w'^2)^2 ] / dz
    ! - (1/rho_ds) * d [ rho_ds * a_1 * ( (w'^3)^2 / w'^2 ) ] / dz.
    !
    ! The (w'^2)^2 and (w'^3)^2 terms are both timestep split so that they can
    ! be expressed linearly in terms of w'^2 and w'^3, respectively, at the
    ! (t+1) timestep, such that:
    !
    ! (w'^2)^2 = w'^2(t) * w'^2(t+1);
    ! (w'^3)^2 = w'^3(t) * w'^3(t+1);
    !
    ! which allows these terms to be expressed implicitly as:
    !
    ! - (1/rho_ds) * d [ rho_ds * a_3 * w'^2(t) * w'^2(t+1) ] / dz
    ! - (1/rho_ds) * d [ rho_ds * a_1 * w'^3(t) * w'^3(t+1) / w'^2(t) ] / dz.
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
    ! mathematical expressions (called F and G here) within the dF/dz and dG/dz
    ! terms are computed on the momentum levels.  Then, the derivatives (d/dz)
    ! of the expressions (F and G) are taken over the central thermodynamic
    ! level, where dF/dz and dG/dz are multiplied by -invrs_rho_ds_zt.  This
    ! yields the desired results.  In this function, the values of F and G are
    ! as follows:
    !
    ! F = rho_ds_zm * a_3(t) * w'^2(t) * w'^2(t+1); and
    !
    ! G = rho_ds_zm * a_1(t) * w'^3(t) * w'^3(t+1) / w'^2(t).
    !
    ! ------------------------------------------------wp3p1-------------- t(k+1)
    !
    ! ===a3====wp2====rho_ds_zm====a1======================wp3(interp)=== m(k)
    !
    ! -----------dF/dz----invrs_rho_ds_zt----dG/dz----wp3---------------- t(k)
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

    use grid_class, only:  &
        gr ! Variable gr%weights_zt2zm

    use model_flags, only:  &
        l_standard_term_ta

    use clubb_precision, only: &
        core_rknd ! Variable(s)

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
      wp2,             & ! w'^2(k)                                     [m^2/s^2]
      wp2m1,           & ! w'^2(k-1)                                   [m^2/s^2]
      a1,              & ! a_1(k)                                      [-]
      a1_zt,           & ! a_1 interpolated to thermodynamic level (k) [-]
      a1m1,            & ! a_1(k-1)                                    [-]
      a3,              & ! a_3(k)                                      [-]
      a3_zt,           & ! a_3 interpolated to thermodynamic level (k) [-]
      a3m1,            & ! a_3(k-1)                                    [-]
      wp3_on_wp2,      & ! w'^3 / w'^2 at momentum level (k)           [m/s]
      wp3_on_wp2_m1,   & ! w'^3 / w'^2 at momentum level (k-1)         [m/s]
      rho_ds_zm,       & ! Dry, static density at momentum level (k)   [kg/m^3]
      rho_ds_zmm1,     & ! Dry, static density at momentum level (k-1) [kg/m^3]
      invrs_rho_ds_zt, & ! Inv dry, static density at thermo level (k) [m^3/kg]
      invrs_dzt          ! Inverse of grid spacing (k)                 [1/m]

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
             * rho_ds_zm * a1 * wp3_on_wp2 &
             * gr%weights_zt2zm(t_above,mk)

       ! Momentum superdiagonal: [ x wp2(k,<t+1>) ]
       lhs(k_mdiag) &
       = + invrs_rho_ds_zt * invrs_dzt * rho_ds_zm * a3 * wp2

       ! Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
       lhs(k_tdiag) &
       = + invrs_rho_ds_zt &
           * invrs_dzt &
             * (   rho_ds_zm * a1 * wp3_on_wp2 &
                   * gr%weights_zt2zm(t_below,mk) &
                 - rho_ds_zmm1 * a1m1 * wp3_on_wp2_m1 &
                   * gr%weights_zt2zm(t_above,mkm1) &
               )

       ! Momentum subdiagonal: [ x wp2(k-1,<t+1>) ]
       lhs(km1_mdiag) &
       = - invrs_rho_ds_zt * invrs_dzt * rho_ds_zmm1 * a3m1 * wp2m1

       ! Thermodynamic subdiagonal: [ x wp3(k-1,<t+1>) ]
       lhs(km1_tdiag) &
       = - invrs_rho_ds_zt &
           * invrs_dzt &
             * rho_ds_zmm1 * a1m1 * wp3_on_wp2_m1 &
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
       !  - (1/rho_ds) * d [ rho_ds * a_3 * (w'^2)^2 ] / dz, has been altered to
       ! pull a_3 outside of the derivative.  This was done in order to help
       ! stabilize w'^3.  On the left-hand side of the equation, this effects
       ! the momentum superdiagonal (k_mdiag) and the momentum subdiagonal
       ! (km1_mdiag).

       ! Thermodynamic superdiagonal: [ x wp3(k+1,<t+1>) ]
       lhs(kp1_tdiag) &
       = + invrs_rho_ds_zt &
           * a1_zt * invrs_dzt &
             * rho_ds_zm * wp3_on_wp2 &
             * gr%weights_zt2zm(t_above,mk)

       ! Momentum superdiagonal: [ x wp2(k,<t+1>) ]
       lhs(k_mdiag) &
       = + invrs_rho_ds_zt * a3_zt * invrs_dzt * rho_ds_zm * wp2

       ! Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
       lhs(k_tdiag) &
       = + invrs_rho_ds_zt &
           * a1_zt * invrs_dzt & 
             * (   rho_ds_zm * wp3_on_wp2 & 
                   * gr%weights_zt2zm(t_below,mk) & 
                 - rho_ds_zmm1 * wp3_on_wp2_m1 & 
                   * gr%weights_zt2zm(t_above,mkm1) & 
               )

       ! Momentum subdiagonal: [ x wp2(k-1,<t+1>) ]
       lhs(km1_mdiag) &
       = - invrs_rho_ds_zt * a3_zt * invrs_dzt * rho_ds_zmm1 * wp2m1

       ! Thermodynamic subdiagonal: [ x wp3(k-1,<t+1>) ]
       lhs(km1_tdiag) &
       = - invrs_rho_ds_zt &
           * a1_zt * invrs_dzt &
             * rho_ds_zmm1 * wp3_on_wp2_m1 & 
             * gr%weights_zt2zm(t_below,mkm1)

       ! End of code that pulls out a3.
       ! End of Brian's a1 change.  Feb. 14, 2008.

    end if ! l_standard_term_ta


    return

  end function wp3_term_ta_ADG1_lhs

    !=============================================================================
    pure subroutine wp3_term_ta_ADG1_lhs_all( wp2, &
                                              a1, a1_zt, &
                                              a3, a3_zt, &
                                              wp3_on_wp2, &
                                              rho_ds_zm, &
                                              invrs_rho_ds_zt, &
                                              invrs_dzt, &
                                              lhs_ta_wp3 )
    ! Description:
    !     This subroutine serves the same function as wp3_term_ta_ADG1_lhs (above), but
    !     calculates terms for all grid levels at once rather than one at a time.
    !     This was done so that this code could be vectorized and thereby sped up
    !     by the compiler. See clubb:ticket:834 for more information.
    ! 
    !--------------------------------------------------------------------------------------

        use grid_class, only:  &
            gr ! Variable gr%weights_zt2zm

        use model_flags, only:  &
            l_standard_term_ta

        use clubb_precision, only: &
            core_rknd ! Variable(s)

        implicit none

        ! Input Variables
        real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
          wp2,             & ! w'^2(k)                                     [m^2/s^2]
          a1,              & ! a_1(k)                                      [-]
          a1_zt,           & ! a_1 interpolated to thermodynamic level (k) [-]
          a3,              & ! a_3(k)                                      [-]
          a3_zt,           & ! a_3 interpolated to thermodynamic level (k) [-]
          wp3_on_wp2,      & ! w'^3 / w'^2 at momentum level (k)           [m/s]
          rho_ds_zm,       & ! Dry, static density at momentum level (k)   [kg/m^3]
          invrs_rho_ds_zt, & ! Inv dry, static density at thermo level (k) [m^3/kg]
          invrs_dzt          ! Inverse of grid spacing (k)                 [1/m]

        ! Return Variable
        real( kind = core_rknd ), dimension(5,gr%nz), intent(out) :: lhs_ta_wp3

        ! Loop variable
        integer :: k

        ! Set lower boundary to 0
        lhs_ta_wp3(:,1) = 0.0_core_rknd


        if ( l_standard_term_ta ) then

            do k = 2, gr%nz-1

                ! Thermodynamic superdiagonal: [ x wp3(k+1,<t+1>) ]
                lhs_ta_wp3(1,k) = + invrs_rho_ds_zt(k) * invrs_dzt(k) * rho_ds_zm(k) &
                                  * a1(k) * wp3_on_wp2(k) * gr%weights_zt2zm(1,k)

                ! Momentum superdiagonal: [ x wp2(k,<t+1>) ]
                lhs_ta_wp3(2,k) = + invrs_rho_ds_zt(k) * invrs_dzt(k) &
                                  * rho_ds_zm(k) * a3(k) * wp2(k)

                ! Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
                lhs_ta_wp3(3,k) = + invrs_rho_ds_zt(k) * invrs_dzt(k) * ( rho_ds_zm(k) &
                                      * a1(k) * wp3_on_wp2(k) * gr%weights_zt2zm(2,k) &
                                      - rho_ds_zm(k-1) * a1(k-1) * wp3_on_wp2(k-1) &
                                      * gr%weights_zt2zm(1,k-1) )

                ! Momentum subdiagonal: [ x wp2(k-1,<t+1>) ]
                lhs_ta_wp3(4,k) = - invrs_rho_ds_zt(k) * invrs_dzt(k) &
                                  * rho_ds_zm(k-1) * a3(k-1) * wp2(k-1)

                ! Thermodynamic subdiagonal: [ x wp3(k-1,<t+1>) ]
                lhs_ta_wp3(5,k) = - invrs_rho_ds_zt(k) * invrs_dzt(k) * rho_ds_zm(k-1) &
                                  * a1(k-1) * wp3_on_wp2(k-1) * gr%weights_zt2zm(2,k-1)

            end do

        else

            do k = 2, gr%nz-1

                ! Thermodynamic superdiagonal: [ x wp3(k+1,<t+1>) ]
                lhs_ta_wp3(1,k) = + invrs_rho_ds_zt(k) * a1_zt(k) * invrs_dzt(k) &
                                  * rho_ds_zm(k) * wp3_on_wp2(k) * gr%weights_zt2zm(1,k)

                ! Momentum superdiagonal: [ x wp2(k,<t+1>) ]
                lhs_ta_wp3(2,k) = + invrs_rho_ds_zt(k) * a3_zt(k) * invrs_dzt(k) &
                                  * rho_ds_zm(k) * wp2(k)

                ! Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
                lhs_ta_wp3(3,k) = + invrs_rho_ds_zt(k) * a1_zt(k) * invrs_dzt(k) & 
                                  * ( rho_ds_zm(k) * wp3_on_wp2(k) * gr%weights_zt2zm(2,k) &
                                    - rho_ds_zm(k-1) * wp3_on_wp2(k-1) * gr%weights_zt2zm(1,k-1) )

                ! Momentum subdiagonal: [ x wp2(k-1,<t+1>) ]
                lhs_ta_wp3(4,k) = - invrs_rho_ds_zt(k) * a3_zt(k) * invrs_dzt(k) &
                                  * rho_ds_zm(k-1) * wp2(k-1)

                ! Thermodynamic subdiagonal: [ x wp3(k-1,<t+1>) ]
                lhs_ta_wp3(5,k) = - invrs_rho_ds_zt(k) * a1_zt(k) * invrs_dzt(k) &
                                  * rho_ds_zm(k-1) * wp3_on_wp2(k-1) * gr%weights_zt2zm(2,k-1)

            end do


        end if ! l_standard_term_ta

        ! Set lower boundary to 0
        lhs_ta_wp3(:,gr%nz) = 0.0_core_rknd


        return

    end subroutine wp3_term_ta_ADG1_lhs_all

  !=============================================================================
  pure function wp3_term_tp_lhs( wp2, wp2m1, &
                                 rho_ds_zm, rho_ds_zmm1, &
                                 invrs_rho_ds_zt, &
                                 invrs_dzt ) &
  result( lhs )

    ! Description:
    ! Turbulent production of w'^3:  implicit portion of the code.
    !
    ! The d(w'^3)/dt equation contains a turbulent production term:
    !
    ! + 3 * ( w'^2 / rho_ds ) * d( rho_ds * w'^2 )/dz.
    !
    ! The turbulent production term is rewritten as:
    !
    ! + 3 * ( w'^2 / rho_ds ) * d[ rho_ds * w'^2 ]/dz
    ! = + (3/rho_ds) * d[ rho_ds * (w'^2)^2 ]/dz - (3/2) * d[ (w'^2)^2 ]/dz.
    !
    ! The (w'^2)^2 terms are timestep split so that they can be expressed
    ! linearly in terms of w'^2 at the (t+1) timestep, such that:
    !
    ! (w'^2)^2 = w'^2(t) * w'^2(t+1).
    !
    ! The term can now be expressed implicitly as:
    !
    ! + (3/rho_ds) * d [ rho_ds * w'^2(t) * w'^2(t+1) ] / dz
    ! - (3/2) * d [ w'^2(t) * w'^2(t+1) ] /dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign is
    !        reversed and the leading "-" in front of a d[ ] / dz term is
    !        changed to a "+".  Likewise, the leading "+" in front of a
    !        d[ ] / dz term is changed to a "-".
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is 
    ! being advanced to in solving the d(w'^3)/dt and d(w'^2)/dt equations.
    !
    ! The implicit portion of these terms is discretized as follows:
    !
    ! While the values of w'^3 are found on the thermodynamic levels, the values
    ! of w'^2 are found on the momentum levels.  Additionally, the values of
    ! rho_ds_zm are found on the momentum levels, and the values of
    ! invrs_rho_ds_zt are found on the thermodynamic levels.  The values of the
    ! mathematical expressions (called F and G below) within the dF/dz and dG/dz
    ! terms are computed on the momentum levels.  Then, the derivatives (d/dz)
    ! of the expressions (F and G) are taken over the central thermodynamic
    ! level, where dF/dz and dG/dz are multiplied by -3 * invrs_rho_ds_zt and
    ! 3/2, respectively, yielding the desired results.  In this function, the
    ! values of F and G are as follows:
    !
    ! F = rho_ds_zm * w'^2(t) * w'^2(t+1);
    !
    ! G = w'^2(t) * w'^2(t+1).
    !
    ! ====wp2=========rho_ds_zm========================================== m(k)
    !
    ! -----------dF/dz----invrs_rho_ds_zt----dG/dz----wp3---------------- t(k)
    !
    ! ====wp2m1=======rho_ds_zmm1======================================== m(k-1)
    !
    ! The vertical indices m(k), t(k), and m(k-1) correspond with altitudes
    ! zm(k), zt(k), and zm(k-1), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) )

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        three,        & ! Constant(s)
        three_halves

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      k_mdiag   = 1, & ! Momentum superdiagonal index.
      km1_mdiag = 2    ! Momentum subdiagonal index. 

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      wp2,             & ! w'^2(k)                                     [m^2/s^2]
      wp2m1,           & ! w'^2(k-1)                                   [m^2/s^2]
      rho_ds_zm,       & ! Dry, static density at momentum level (k)   [kg/m^3]
      rho_ds_zmm1,     & ! Dry, static density at momentum level (k-1) [kg/m^3]
      invrs_rho_ds_zt, & ! Inv dry, static density at thermo level (k) [m^3/kg]
      invrs_dzt          ! Inverse of grid spacing (k)                 [1/m]

    ! Return Variable
    real( kind = core_rknd ), dimension(2) :: lhs


    ! Momentum superdiagonal: [ x wp2(k,<t+1>) ]
    lhs(k_mdiag) &
    = - three * invrs_rho_ds_zt * invrs_dzt * rho_ds_zm * wp2 &
      + three_halves * invrs_dzt * wp2

    ! Momentum subdiagonal: [ x wp2(k-1,<t+1>) ]
    lhs(km1_mdiag) &
    = + three * invrs_rho_ds_zt * invrs_dzt * rho_ds_zmm1 * wp2m1 &
      - three_halves * invrs_dzt * wp2m1


    return

  end function wp3_term_tp_lhs

    !==================================================================================
    pure subroutine wp3_term_tp_lhs_all( wp2, &
                                       rho_ds_zm, &
                                       invrs_rho_ds_zt, &
                                       invrs_dzt, &
                                       lhs_tp_wp3 )
    ! Description:
    !     This subroutine serves the same function as wp3_term_tp_lhs (above), but
    !     calculates terms for all grid levels at once rather than one at a time.
    !     This was done so that this code could be vectorized and thereby sped up
    !     by the compiler. See clubb:ticket:834 for more information.
    ! 
    !----------------------------------------------------------------------------------

        use constants_clubb, only:  &
            three,        & ! Constant(s)
            three_halves

        use clubb_precision, only: &
            core_rknd ! Variable(s)

        use grid_class, only:  & 
            gr       ! Variable(s)

        implicit none
     

        ! Input Variables
        real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
          wp2,             & ! w'^2(k)                                     [m^2/s^2]
          rho_ds_zm,       & ! Dry, static density at momentum level (k)   [kg/m^3]
          invrs_rho_ds_zt, & ! Inv dry, static density at thermo level (k) [m^3/kg]
          invrs_dzt          ! Inverse of grid spacing (k)                 [1/m]

        ! Return Variable
        real( kind = core_rknd ), dimension(2,gr%nz), intent(out) :: lhs_tp_wp3

        ! Loop variable
        integer :: k

        ! Set lower boundary to 0
        lhs_tp_wp3(1,1) = 0.0_core_rknd
        lhs_tp_wp3(2,1) = 0.0_core_rknd
        
        ! Calculate non-boundary values
        do k = 2, gr%nz-1

            ! Momentum superdiagonal: [ x wp2(k,<t+1>) ]
            lhs_tp_wp3(1,k) = - three * invrs_rho_ds_zt(k) * invrs_dzt(k) &
                              * rho_ds_zm(k) * wp2(k) + three_halves * invrs_dzt(k) * wp2(k)

            ! Momentum subdiagonal: [ x wp2(k-1,<t+1>) ]
            lhs_tp_wp3(2,k) = + three * invrs_rho_ds_zt(k) * invrs_dzt(k) &
                              * rho_ds_zm(k-1) * wp2(k-1) - three_halves * invrs_dzt(k) * wp2(k-1)

        end do

        ! Set upper boundary to 0
        lhs_tp_wp3(1,gr%nz) = 0.0_core_rknd
        lhs_tp_wp3(2,gr%nz) = 0.0_core_rknd


        return

    end subroutine wp3_term_tp_lhs_all

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

    use constants_clubb, only: &
        three, & ! Variable(s)
        one

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

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
    = + ( one - C11_Skw_fnc ) * three * invrs_dzt * ( wm_zm - wm_zmm1 )

    return
  end function wp3_terms_ac_pr2_lhs

    !==================================================================================
    pure subroutine wp3_terms_ac_pr2_lhs_all( C11_Skw_fnc, wm_zm, invrs_dzt, &
                                              lhs_ac_pr2_wp3 )
    ! Description:
    !     This subroutine serves the same function as wp3_term_tp_lhs (above), but
    !     calculates terms for all grid levels at once rather than one at a time.
    !     This was done so that this code could be vectorized and thereby sped up
    !     by the compiler. See clubb:ticket:834 for more information.
    ! 
    !----------------------------------------------------------------------------------

        use constants_clubb, only: &
            three, & ! Variable(s)
            one

        use clubb_precision, only: &
            core_rknd    ! Variable(s)

        use grid_class, only:  & 
            gr       ! Variable(s)

        implicit none

        ! Input Variables
        real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
          C11_Skw_fnc,  & ! C_11 parameter with Sk_w applied (k)      [-]
          wm_zm,        & ! w wind component at momentum levels (k)   [m/s]
          invrs_dzt       ! Inverse of grid spacing (k)               [1/m]

        ! Return Variable
        real( kind = core_rknd ), dimension(gr%nz), intent(out) :: & 
            lhs_ac_pr2_wp3

        ! Loop variable
        integer :: k

        ! Set lower boundary to 0
        lhs_ac_pr2_wp3(1) = 0.0_core_rknd

        ! Calculate non-boundary terms
        do k = 2, gr%nz-1

            ! Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
            lhs_ac_pr2_wp3(k) = + ( one - C11_Skw_fnc(k) ) * three &
                                * invrs_dzt(k) * ( wm_zm(k) - wm_zm(k-1) )

        end do

        ! Set upper boundary to 0
        lhs_ac_pr2_wp3(gr%nz) = 0.0_core_rknd

        return

    end subroutine wp3_terms_ac_pr2_lhs_all

  !=============================================================================
  pure function wp3_term_pr1_lhs( C8, C8b, tauw3t, Skw_zt ) & 
  result( lhs )

    ! Description:
    ! Pressure term 1 for w'^3:  implicit portion of the code.
    !
    ! Pressure term 1 is the term:
    !
    ! - (C_8/tau_w3t) * ( C_8b * Sk_wt^2 + 1 ) * w'^3;
    !
    ! where Sk_wt = w'^3 / (w'^2)^(3/2).
    !
    ! This term needs to be linearized, so function L(w'^3) is defined to be 
    ! equal to this term (pressure term 1), such that:
    !
    ! L(w'^3) = - (C_8/tau_w3t) * ( C_8b * (w'^3)^3 / (w'^2)^3 + w'^3 ).
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
    ! - (C_8/tau_w3t) * ( 3 * C_8b * Sk_wt^2 + 1 ) * w'^3(t+1).
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

    use constants_clubb, only: &
        one, & ! Variable(s)
        three, &
        five

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use model_flags, only: &
        l_damp_wp3_Skw_squared

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
    if ( l_damp_wp3_Skw_squared ) then 
        lhs & 
        = + ( C8 / tauw3t ) * ( three * C8b * Skw_zt**2 + one )
    else
        lhs & 
        = + ( C8 / tauw3t ) * ( five * C8b * Skw_zt**4 + one )
    end if

    return
  end function wp3_term_pr1_lhs

    !=============================================================================
    pure subroutine wp3_term_pr1_lhs_all( C8, C8b, tauw3t, Skw_zt, &
                                          lhs_pr1_wp3 )
    ! Description:
    !     This subroutine serves the same function as wp3_term_pr1_lhs (above), but
    !     calculates terms for all grid levels at once rather than one at a time.
    !     This was done so that this code could be vectorized and thereby sped up
    !     by the compiler. See clubb:ticket:834 for more information.
    ! 
    !----------------------------------------------------------------------------------

        use constants_clubb, only: &
            one, & ! Variable(s)
            three, &
            five

        use clubb_precision, only: &
            core_rknd ! Variable(s)

        use model_flags, only: &
            l_damp_wp3_Skw_squared

        use grid_class, only:  & 
            gr      ! Variable(s)

        implicit none

        ! Input Variables
        real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
          tauw3t,  & ! Time-scale tau at thermodynamic levels (k) [s]
          Skw_zt     ! Skewness of w at thermodynamic levels (k)  [-]

        real( kind = core_rknd ), intent(in) :: & 
          C8,      & ! Model parameter C_8                        [-]
          C8b        ! Model parameter C_8b                       [-]

        ! Return Variable
        real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
            lhs_pr1_wp3

        ! Loop variable
        integer :: k

        ! Set lower boundary to 0
        lhs_pr1_wp3(1) = 0.0_core_rknd

        
        if ( l_damp_wp3_Skw_squared ) then 

            ! Calculate non-boundary values using Skw_zt^2
            do k = 2, gr%nz-1

                ! Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
                lhs_pr1_wp3(k) = + ( C8 / tauw3t(k) ) * ( three * C8b * Skw_zt(k)**2 + one )

            end do

        else
            
            ! Calculate non-boundary values using Skw_zt^4
            do k = 2, gr%nz-1

                ! Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
                lhs_pr1_wp3(k) = + ( C8 / tauw3t(k) ) * ( five * C8b * Skw_zt(k)**4 + one )

            end do

        end if
        
        ! Set upper boundary to 0
        lhs_pr1_wp3(gr%nz) = 0.0_core_rknd

        return

    end subroutine wp3_term_pr1_lhs_all

  !=============================================================================
  pure function wp3_term_ta_explicit_rhs( wp4, wp4m1, &
                                          rho_ds_zm, rho_ds_zmm1, &
                                          invrs_rho_ds_zt, &
                                          invrs_dzt ) &
  result( rhs )

    ! Description:
    ! Turbulent advection of <w'^3>:  explicit portion of the code.
    !
    ! This explicit discretization works generally for any PDF.
    !
    ! The d<w'^3>/dt equation contains a turbulent advection term:
    !
    ! - (1/rho_ds) * d( rho_ds * <w'^4> )/dz.
    !
    ! The value of <w'^4> is found by integrating over the PDF of w, as detailed
    ! in function calc_wp4_pdf, which is found in module pdf_closure_module in
    ! pdf_closure_module.F90.
    !
    ! The explicit discretization of this term is as follows:
    !
    ! The values of <w'^3> are found on the thermodynamic levels, while the
    ! values of <w'^4> are found on the momentum levels.  The values of
    ! <w'^4>|_zt are originally calculated by the PDF on the thermodynamic
    ! levels.  They are interpolated to the intermediate momentum levels as
    ! <w'^4>.  Additionally, the values of rho_ds_zm are found on the momentum
    ! levels, and the values of invrs_rho_ds_zt are found on the thermodynamic
    ! levels.  At the intermediate momentum levels, the values of <w'^4> are
    ! multiplied by rho_ds_zm.  Then, the derivative (d/dz) of that expression
    ! is taken over the central thermodynamic level, where it is multiplied by
    ! -invrs_rho_ds_zt.  This yields the desired result.
    !
    ! ---------wp4_zt---------------------------------------------------- t(k+1)
    !
    ! =========wp4(interp)===========rho_ds_zm=========================== m(k)
    !
    ! ---------wp4_zt-----d( rho_ds_zm * wp4 )/dz-----invrs_rho_ds_zt---- t(k)
    !
    ! =========wp4m1(interp)=========rho_ds_zmm1========================= m(k-1)
    !
    ! ---------wp4_zt---------------------------------------------------- t(k-1)
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

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      wp4,             & ! <w'^4>(k)                                   [m^4/s^4]
      wp4m1,           & ! <w'^4>(k-1)                                 [m^4/s^4]
      rho_ds_zm,       & ! Dry, static density at momentum level (k)   [kg/m^3]
      rho_ds_zmm1,     & ! Dry, static density at momentum level (k-1) [kg/m^3]
      invrs_rho_ds_zt, & ! Inv dry, static density at thermo level (k) [m^3/kg]
      invrs_dzt          ! Inverse of grid spacing (k)                 [1/m]

    ! Return Variable
    real( kind = core_rknd ) :: rhs


    rhs &
    = - invrs_rho_ds_zt * invrs_dzt * ( rho_ds_zm * wp4 - rho_ds_zmm1 * wp4m1 )


    return

  end function wp3_term_ta_explicit_rhs

    !==================================================================================
    pure subroutine wp3_term_ta_explicit_rhs_all( wp4, &
                                                  rho_ds_zm, &
                                                  invrs_rho_ds_zt, &
                                                  invrs_dzt, &
                                                  rhs_ta_wp3 )
    ! Description:
    !     This subroutine serves the same function as wp3_term_ta_explicit_rhs (above), but
    !     calculates terms for all grid levels at once rather than one at a time.
    !     This was done so that this code could be vectorized and thereby sped up
    !     by the compiler. See clubb:ticket:834 for more information.
    ! 
    !----------------------------------------------------------------------------------

        use clubb_precision, only: &
            core_rknd ! Variable(s)

        use grid_class, only: &
            gr

        implicit none

        ! Input Variables
        real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
          wp4,             & ! <w'^4>(k)                                   [m^4/s^4]
          rho_ds_zm,       & ! Dry, static density at momentum level (k)   [kg/m^3]
          invrs_rho_ds_zt, & ! Inv dry, static density at thermo level (k) [m^3/kg]
          invrs_dzt          ! Inverse of grid spacing (k)                 [1/m]

        ! Return Variable
        real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
            rhs_ta_wp3

        ! Loop variable
        integer :: k

        ! Set lower boundary to 0
        rhs_ta_wp3(1) = 0.0_core_rknd

        ! Calculate non-boundary values
        do k = 2, gr%nz
            rhs_ta_wp3(k) = - invrs_rho_ds_zt(k) * invrs_dzt(k) &
                            * ( rho_ds_zm(k) * wp4(k) - rho_ds_zm(k-1) * wp4(k-1) )
        end do

        ! Set upper boundary to 0
        rhs_ta_wp3(gr%nz) = 0.0_core_rknd


        return

    end subroutine wp3_term_ta_explicit_rhs_all

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
        grav,  & ! Gravitational acceleration [m/s^2]
        three, &
        one

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C11_Skw_fnc, & ! C_11 parameter with Sk_w applied (k)        [-]
      thv_ds_zt,   & ! Dry, base-state theta_v at thermo. lev. (k) [K]
      wp2thvp        ! w'^2th_v'(k)                                [K m^2/s^2]

    ! Return Variable
    real( kind = core_rknd ) :: rhs

    rhs & 
    = + ( one - C11_Skw_fnc ) * three * ( grav / thv_ds_zt ) * wp2thvp

    return
  end function wp3_terms_bp1_pr2_rhs

    !==================================================================================
    pure subroutine wp3_terms_bp1_pr2_rhs_all( C11_Skw_fnc, thv_ds_zt, wp2thvp, &
                                               rhs_bp1_pr2_wp3 )
    ! Description:
    !     This subroutine serves the same function as wp3_terms_bp1_pr2_rhs (above), but
    !     calculates terms for all grid levels at once rather than one at a time.
    !     This was done so that this code could be vectorized and thereby sped up
    !     by the compiler. See clubb:ticket:834 for more information.
    ! 
    !----------------------------------------------------------------------------------

        use clubb_precision, only: &
            core_rknd ! Variable(s)

        use constants_clubb, only: & ! Constant(s) 
            grav,  & ! Gravitational acceleration [m/s^2]
            three, &
            one

        use grid_class, only: &
            gr

        implicit none

        ! Input Variables
        real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
          C11_Skw_fnc, & ! C_11 parameter with Sk_w applied (k)        [-]
          thv_ds_zt,   & ! Dry, base-state theta_v at thermo. lev. (k) [K]
          wp2thvp        ! w'^2th_v'(k)                                [K m^2/s^2]

        ! Return Variable
        real( kind = core_rknd ),  dimension(gr%nz), intent(out) :: &
            rhs_bp1_pr2_wp3

        ! Loop variable
        integer :: k

        ! Set lower boundary to 0
        rhs_bp1_pr2_wp3(1) = 0.0_core_rknd

        ! Calculate non-boundary values
        do k = 2, gr%nz-1
            rhs_bp1_pr2_wp3(k) = + ( one - C11_Skw_fnc(k) ) * three &
                                 * ( grav / thv_ds_zt(k) ) * wp2thvp(k)
        end do

        ! Set upper boundary to 0
        rhs_bp1_pr2_wp3(gr%nz) = 0.0_core_rknd

        return
    end subroutine wp3_terms_bp1_pr2_rhs_all

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

    !==================================================================================
    pure subroutine wp3_term_bp2_rhs_all( C15, Kh_zt, wpthvp, &
                                          dum_dz, dvm_dz, &
                                          upwp, vpwp, &
                                          thv_ds_zt, invrs_dzt, &
                                          rhs_bp2_wp3 )
    ! Description:
    !     This subroutine serves the same function as wp3_term_bp2_rhs (above), but
    !     calculates terms for all grid levels at once rather than one at a time.
    !     This was done so that this code could be vectorized and thereby sped up
    !     by the compiler. See clubb:ticket:834 for more information.
    ! 
    !----------------------------------------------------------------------------------


        use clubb_precision, only: &
            core_rknd ! Variable(s)

        use constants_clubb, only: & ! Constant(s) 
            grav ! Gravitational acceleration [m/s^2]

        use grid_class, only: &
            gr

        implicit none

        
        ! Input Variables
        real( kind = core_rknd ), intent(in) :: &
          C15          ! Model parameter C15                [-]

        real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
          Kh_zt,     & ! Eddy-diffusivity on moment. levels [m^2/s]
          wpthvp,    & ! w'th_v'(k)                         [K m/s]
          dum_dz,    & ! d u wind dz (k)                    [m/s]
          dvm_dz,    & ! d v wind dz (k)                    [m/s]
          upwp,      & ! u'v'(k)                            [m^2/s^2]
          vpwp,      & ! v'w'(k)                            [m^2/s^2]
          thv_ds_zt, & ! Dry, base-state theta_v at thermo. lev. (k) [K]
          invrs_dzt    ! Inverse of grid spacing (k)        [1/m]

        ! Return Variable
        real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
            rhs_bp2_wp3

        ! Loop variable
        integer :: k

        ! ---- Begin Code ----

        ! Set lower boundary to 0
        rhs_bp2_wp3(1) = 0.0_core_rknd

        ! Calculate non-boundary values
        do k = 2, gr%nz-1

            rhs_bp2_wp3(k) = - C15 * Kh_zt(k) * invrs_dzt(k) &
                               * ( grav / thv_ds_zt(k) * ( wpthvp(k) - wpthvp(k-1) ) &
                                 - ( upwp(k) * dum_dz(k) - upwp(k-1) * dum_dz(k-1) ) &
                                 - ( vpwp(k) * dvm_dz(k) - vpwp(k-1) * dvm_dz(k-1) ) )

        end do

        ! Set upper boundary to 0
        rhs_bp2_wp3(gr%nz) = 0.0_core_rknd

        return
    end subroutine wp3_term_bp2_rhs_all


  !=============================================================================
  pure function wp3_term_pr1_rhs( C8, C8b, tauw3t, Skw_zt, wp3 ) & 
  result( rhs )

    ! Description:
    ! Pressure term 1 for w'^3:  explicit portion of the code.
    !
    ! Pressure term 1 is the term:
    !
    ! - (C_8/tau_w3t) * ( C_8b * Sk_wt^2 + 1 ) * w'^3;
    !
    ! where Sk_wt = w'^3 / (w'^2)^(3/2).
    !
    ! This term needs to be linearized, so function L(w'^3) is defined to be 
    ! equal to this term (pressure term 1), such that:
    !
    ! L(w'^3) = - (C_8/tau_w3t) * ( C_8b * (w'^3)^3 / (w'^2)^3 + w'^3 ).
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
    ! + (C_8/tau_w3t) * ( 2 * C_8b * Sk_wt^2 + 1 ) * w'^3(t).
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

    use constants_clubb, only: &
        two, &
        four

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use model_flags, only: &
        l_damp_wp3_Skw_squared

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

    if ( l_damp_wp3_Skw_squared ) then 
        rhs & 
        = + ( C8 / tauw3t ) * ( two * C8b * Skw_zt**2 ) * wp3
    else 
        rhs & 
        = + ( C8 / tauw3t ) * ( four * C8b * Skw_zt**4 ) * wp3
    end if

    return
  end function wp3_term_pr1_rhs

    !==================================================================================
    pure subroutine wp3_term_pr1_rhs_all( C8, C8b, tauw3t, Skw_zt, wp3, &
                                          rhs_pr1_wp3 )
    ! Description:
    !     This subroutine serves the same function as wp3_term_pr1_rhs (above), but
    !     calculates terms for all grid levels at once rather than one at a time.
    !     This was done so that this code could be vectorized and thereby sped up
    !     by the compiler. See clubb:ticket:834 for more information.
    ! 
    !----------------------------------------------------------------------------------

        use constants_clubb, only: &
            two, &
            four

        use clubb_precision, only: &
            core_rknd ! Variable(s)

        use model_flags, only: &
            l_damp_wp3_Skw_squared

        use grid_class, only: &
            gr

        implicit none

        ! Input Variables
        real( kind = core_rknd ), intent(in) :: & 
          C8,      & ! Model parameter C_8                        [-]
          C8b        ! Model parameter C_8b                       [-]

        ! Input Variables
        real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
          tauw3t,  & ! Time-scale tau at thermodynamic levels (k) [s]
          Skw_zt,  & ! Skewness of w at thermodynamic levels (k)  [-]
          wp3        ! w'^3(k)                                    [m^3/s^3]

        ! Return Variable
        real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
            rhs_pr1_wp3

        ! Loop variable
        integer :: k

        ! Set lower boundary to 0
        rhs_pr1_wp3(1) = 0.0_core_rknd
        

        ! Calculate non-boundary values
        if ( l_damp_wp3_Skw_squared ) then 

            do k = 2, gr%nz-1
                rhs_pr1_wp3(k) = + ( C8 / tauw3t(k) ) * ( two * C8b * Skw_zt(k)**2 ) * wp3(k)
            end do

        else 

            do k = 2, gr%nz-1
                rhs_pr1_wp3(k) = + ( C8 / tauw3t(k) ) * ( four * C8b * Skw_zt(k)**4 ) * wp3(k)
            end do

        end if

        ! Set upper boundary to 0
        rhs_pr1_wp3(gr%nz) = 0.0_core_rknd

        return
    end subroutine wp3_term_pr1_rhs_all

!===============================================================================

end module advance_wp2_wp3_module
