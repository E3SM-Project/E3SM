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
             wp3_term_pr_turb_rhs, &
             wp3_term_pr_dfsn_rhs


  ! Private named constants to avoid string comparisons
  integer, parameter, private :: &
    clip_wp2 = 12 ! Named constant for wp2 clipping.
                  ! NOTE: This must be the same as the clip_wp2 declared in
                  ! clip_explicit!

  contains

  !=============================================================================
  subroutine advance_wp2_wp3( gr, dt, sfc_elevation, sigma_sqd_w, wm_zm,     & ! In
                              wm_zt, a3, a3_zt, wp3_on_wp2,                  & ! In
                              wpup2, wpvp2, wp2up2, wp2vp2, wp4,             & ! In
                              wpthvp, wp2thvp, um, vm, upwp, vpwp,           & ! In
                              up2, vp2, em, Kh_zm, Kh_zt, invrs_tau_C4_zm,   & ! In
                              invrs_tau_wp3_zt, invrs_tau_C1_zm, Skw_zm,     & ! In
                              Skw_zt, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, & ! In
                              invrs_rho_ds_zt, radf, thv_ds_zm,              & ! In
                              thv_ds_zt, mixt_frac, Cx_fnc_Richardson,       & ! In
                              wp2_splat, wp3_splat,                          & ! In
                              pdf_implicit_coefs_terms,                      & ! In
                              wprtp, wpthlp, rtp2, thlp2,                    & ! In
                              clubb_params, nu_vert_res_dep,                 & ! In
                              iiPDF_type,                                    & ! In
                              l_min_wp2_from_corr_wx,                        & ! In
                              l_upwind_xm_ma,                                & ! In
                              l_tke_aniso,                                   & ! In
                              l_standard_term_ta,                            & ! In
                              l_partial_upwind_wp3,                          & ! In
                              l_damp_wp2_using_em,                           & ! In
                              l_use_C11_Richardson,                          & ! In
                              l_damp_wp3_Skw_squared,                        & ! In
                              l_lmm_stepping,                                & ! In
                              l_use_tke_in_wp3_pr_turb_term,                 & ! In
                              l_use_tke_in_wp2_wp3_K_dfsn,                   & ! In
                              stats_zt, stats_zm, stats_sfc,                 & ! intent(inout)
                              wp2, wp3, wp3_zm, wp2_zt )                       ! Inout

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
        grid, & ! Type
        zt2zm,  & ! Procedure(s)
        zm2zt

    use parameter_indices, only: &
        nparams, & ! Variable(s)
        iC11c,   &
        iC11b,   & 
        iC11,    & 
        iC1c,    & 
        iC1b,    & 
        iC1,     & 
        ic_K1,   & 
        ic_K8

    use parameters_tunable, only: &
        nu_vertical_res_dep    ! Type(s)

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

    use stats_type, only: stats ! Type

    implicit none

    type (stats), target, intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc

    type (grid), target, intent(in) :: gr

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
      wpup2,             & ! w'u'^2 (thermodynamic levels)             [m^3/s^3]
      wpvp2,             & ! w'v'^2 (thermodynamic levels)             [m^3/s^3]
      wp2up2,            & ! w'^2u'^2 (momentum levels)                [m^4/s^4]
      wp2vp2,            & ! w'^2v'^2 (momentum levels)                [m^4/s^4]
      wp4,               & ! w'^4 (momentum levels)                    [m^4/s^4]
      wpthvp,            & ! w'th_v' (momentum levels)                 [K m/s]
      wp2thvp,           & ! w'^2th_v' (thermodynamic levels)          [K m^2/s^2]
      um,                & ! u wind component (thermodynamic levels)   [m/s]
      vm,                & ! v wind component (thermodynamic levels)   [m/s]
      upwp,              & ! u'w' (momentum levels)                    [m^2/s^2]
      vpwp,              & ! v'w' (momentum levels)                    [m^2/s^2]
      up2,               & ! u'^2 (momentum levels)                    [m^2/s^2]
      vp2,               & ! v'^2 (momentum levels)                    [m^2/s^2]
      em,                & ! Turbulence kinetic energy                 [m^2/s^2]
      Kh_zm,             & ! Eddy diffusivity on momentum levels       [m^2/s]
      Kh_zt,             & ! Eddy diffusivity on thermodynamic levels  [m^2/s]
      invrs_tau_C4_zm,   & ! Inverse time-scale tau on momentum levels         [1/s]
      invrs_tau_wp3_zt,  & ! Inverse time-scale tau on thermodynamic levels    [1/s]
      invrs_tau_C1_zm,   & ! Inverse tau values used for the C1 (dp1) term in wp2 [1/s]
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

    type(implicit_coefs_terms), intent(in) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]
 
    type(nu_vertical_res_dep), intent(in) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    integer, intent(in) :: &
      iiPDF_type    ! Selected option for the two-component normal (double
                    ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                    ! w, chi, and eta) portion of CLUBB's multivariate,
                    ! two-component PDF.

    logical, intent(in) :: &
      l_min_wp2_from_corr_wx,     & ! Flag to base the threshold minimum value of wp2 on keeping the
                                    ! overall correlation of w and x (w and rt, as well as w and
                                    ! theta-l) within the limits of -max_mag_correlation_flux to
                                    ! max_mag_correlation_flux.
      l_upwind_xm_ma,             & ! This flag determines whether we want to use an upwind
                                    ! differencing approximation rather than a centered differencing
                                    ! for turbulent or mean advection terms. It affects rtm, thlm,
                                    ! sclrm, um and vm.
      l_tke_aniso,                & ! For anisotropic turbulent kinetic energy, i.e. TKE = 1/2
                                    ! (u'^2 + v'^2 + w'^2)
      l_standard_term_ta,         & ! Use the standard discretization for the turbulent advection
                                    ! terms. Setting to .false. means that a_1 and a_3 are pulled
                                    ! outside of the derivative in advance_wp2_wp3_module.F90 and in
                                    ! advance_xp2_xpyp_module.F90.
      l_partial_upwind_wp3,       & ! Flag to use an "upwind" discretization rather
                                    ! than a centered discretization for the portion
                                    ! of the wp3 turbulent advection term for ADG1
                                    ! that is linearized in terms of wp3<t+1>.
                                    ! (Requires ADG1 PDF and l_standard_term_ta).
      l_damp_wp2_using_em,        & ! In wp2 equation, use a dissipation formula of 
                                    ! -(2/3)*em/tau_zm,
                                    ! as in Bougeault (1981)
      l_use_C11_Richardson,       & ! Parameterize C16 based on Richardson number
      l_damp_wp3_Skw_squared,     & ! Set damping on wp3 to use Skw^2 rather than Skw^4
      l_lmm_stepping,             & ! Apply Linear Multistep Method (LMM) Stepping
      l_use_tke_in_wp3_pr_turb_term, & ! Use TKE formulation for wp3 pr_turb term
      l_use_tke_in_wp2_wp3_K_dfsn   ! Use TKE in eddy diffusion for wp2 and wp3

    ! Input/Output
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) ::  & 
      wp2,  & ! w'^2 (momentum levels)                    [m^2/s^2]
      wp3,  & ! w'^3 (thermodynamic levels)               [m^3/s^3]
      wp3_zm  ! w'^3 interpolated to momentum levels      [m^3/s^3]

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) ::  &
      wp2_zt  ! w'^2 interpolated to thermodyamic levels  [m^2/s^2]

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) ::  &
      wp2_old, & ! w'^2 (momentum levels)                 [m^2/s^2]
      wp3_old    ! w'^3 (thermodynamic levels)            [m^3/s^3] 

    ! Eddy Diffusion for w'^2 and w'^3.
    real( kind = core_rknd ), dimension(gr%nz) :: Kw1    ! w'^2 coef. eddy diff.  [m^2/s]
    real( kind = core_rknd ), dimension(gr%nz) :: Kw8    ! w'^3 coef. eddy diff.  [m^2/s]

    ! Internal variables for C11 function, Vince Larson 13 Mar 2005
    real( kind = core_rknd ), dimension(gr%nz) ::  & 
      C1_Skw_fnc,  & ! C_1 parameter with Sk_w applied              [-]
      C11_Skw_fnc, & ! C_11 parameter with Sk_w applied             [-]
    ! End Vince Larson's addition.
      C16_fnc        ! C_16 parameter                               [-]

    real( kind = core_rknd ) ::  &
      C1,   & ! CLUBB tunable parameter C1
      C1b,  & ! CLUBB tunable parameter C1b
      C1c,  & ! CLUBB tunable parameter C1c
      C11,  & ! CLUBB tunable parameter C11
      C11b, & ! CLUBB tunable parameter C11b
      C11c, & ! CLUBB tunable parameter C11c
      c_K1, & ! CLUBB tunable parameter c_K1 
      c_K8    ! CLUBB tunable parameter c_K8

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

    ! Vince Larson added code to make C11 function of Skw. 13 Mar 2005
    ! If this code is used, C11 is no longer relevant, i.e. constants
    !    are hardwired.

    if ( l_use_C11_Richardson ) then

      C11_Skw_fnc = Cx_fnc_Richardson

    else

      ! Unpack CLUBB tunable parameters
      C11 = clubb_params(iC11)
      C11b = clubb_params(iC11b)
      C11c = clubb_params(iC11c)

      ! Calculate C_{1} and C_{11} as functions of skewness of w.
      ! The if..then here is only for computational efficiency -dschanen 2 Sept 08
      if ( abs(C11-C11b) > abs(C11+C11b)*eps/2 ) then
        C11_Skw_fnc(1:gr%nz) =  & 
          C11b + (C11-C11b)*EXP( -one_half * (Skw_zt(1:gr%nz)/C11c)**2 )
      else
        C11_Skw_fnc(1:gr%nz) = C11b
      end if

    end if ! l_use_C11_Richardson

    ! Unpack CLUBB tunable parameters
    C1 = clubb_params(iC1)
    C1b = clubb_params(iC1b)
    C1c = clubb_params(iC1c)

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
      call stat_update_var( iC11_Skw_fnc, C11_Skw_fnc, & ! intent(in)
                            stats_zt )                   ! intent(inout)
      call stat_update_var( iC1_Skw_fnc, C1_Skw_fnc, &   ! intent(in)
                            stats_zm )                   ! intent(inout)
    endif

    ! Unpack CLUBB tunable parameters
    c_K1 = clubb_params(ic_K1)
    c_K8 = clubb_params(ic_K8)

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

    if ( l_lmm_stepping ) then
       wp2_old=wp2
       wp3_old=wp3
    endif ! l_lmm_stepping

    ! Solve semi-implicitly
    call wp23_solve( gr, dt, sfc_elevation, sigma_sqd_w, wm_zm,               & ! Intent(in)
                     wm_zt, a3, a3_zt, wp3_on_wp2,                            & ! Intent(in)
                     wpup2, wpvp2, wp2up2, wp2vp2, wp4,                       & ! Intent(in)
                     wpthvp, wp2thvp, um, vm, upwp, vpwp,                     & ! Intent(in)
                     up2, vp2, em, Kw1, Kw8, Kh_zt, Skw_zt,                   & ! Intent(in)
                     invrs_tau_C4_zm, invrs_tau_wp3_zt,                       & ! Intent(in)
                     invrs_tau_C1_zm, C1_Skw_fnc,                             & ! Intent(in)
                     C11_Skw_fnc, rho_ds_zm,                                  & ! Intent(in)
                     rho_ds_zt, invrs_rho_ds_zm,                              & ! Intent(in)
                     invrs_rho_ds_zt, radf,                                   & ! Intent(in)
                     thv_ds_zm, thv_ds_zt,                                    & ! Intent(in)
                     wp2_splat, wp3_splat,                                    & ! Intent(in)
                     pdf_implicit_coefs_terms,                                & ! Intent(in)
                     wprtp, wpthlp, rtp2, thlp2,                              & ! Intent(in)
                     clubb_params, nu_vert_res_dep,                           & ! Intent(in)
                     iiPDF_type,                                              & ! Intent(in)
                     l_min_wp2_from_corr_wx,                                  & ! Intent(in)
                     l_upwind_xm_ma,                                          & ! Intent(in)
                     l_tke_aniso,                                             & ! Intent(in)
                     l_standard_term_ta,                                      & ! Intent(in)
                     l_partial_upwind_wp3,                                    & ! Intent(in)
                     l_damp_wp2_using_em,                                     & ! Intent(in)
                     l_damp_wp3_Skw_squared,                                  & ! Intent(in)
                     l_use_tke_in_wp3_pr_turb_term,                           & ! Intent(in)
                     l_use_tke_in_wp2_wp3_K_dfsn,                             & ! Intent(in)
                     stats_zt, stats_zm, stats_sfc,                           & ! intent(inout)
                     wp2, wp3, wp3_zm, wp2_zt )                                 ! Intent(inout)

    if ( l_lmm_stepping ) then
       wp2 = one_half * ( wp2_old + wp2 )
       wp3 = one_half * ( wp3_old + wp3 )
    endif ! l_lmm_stepping

    ! When selected, apply sponge damping after wp2 and wp3 have been advanced.
    if ( wp2_sponge_damp_settings%l_sponge_damping ) then

       if ( l_stats_samp ) then
          call stat_begin_update( gr, iwp2_sdmp, wp2 / dt, & ! intent(in)
                                  stats_zm )             ! intent(inout)
       endif

       wp2 = sponge_damp_xp2( gr, dt, gr%zm, wp2, w_tol_sqd, &
                              wp2_sponge_damp_profile )

       if ( l_stats_samp ) then
          call stat_end_update( gr, iwp2_sdmp, wp2 / dt, & ! intent(in)
                                stats_zm )             ! intent(inout)
       endif

    endif ! wp2_sponge_damp_settings%l_sponge_damping

    if ( wp3_sponge_damp_settings%l_sponge_damping ) then

       if ( l_stats_samp ) then
          call stat_begin_update( gr, iwp3_sdmp, wp3 / dt, & ! intent(in)
                                  stats_zt )             ! intent(inout)
       endif

       wp3 = sponge_damp_xp3( gr, dt, gr%zt, wp3, wp3_sponge_damp_profile )

       if ( l_stats_samp ) then
          call stat_end_update( gr, iwp3_sdmp, wp3 / dt, & ! intent(in) 
                                stats_zt )             ! intent(inout)
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
            write(fstderr,*) "wpup2 = ", wpup2, new_line('c')
            write(fstderr,*) "wpvp2 = ", wpvp2, new_line('c')
            write(fstderr,*) "wp2up2 = ", wp2up2, new_line('c')
            write(fstderr,*) "wp2vp2 = ", wp2vp2, new_line('c')
            write(fstderr,*) "wp4 = ", wp4, new_line('c')
            write(fstderr,*) "wpthvp = ", wpthvp, new_line('c')
            write(fstderr,*) "wp2thvp = ", wp2thvp, new_line('c')
            write(fstderr,*) "um = ", um, new_line('c')
            write(fstderr,*) "vm = ", vm, new_line('c')
            write(fstderr,*) "upwp = ", upwp, new_line('c')
            write(fstderr,*) "vpwp = ", vpwp, new_line('c')
            write(fstderr,*) "up2 = ", up2, new_line('c')
            write(fstderr,*) "vp2 = ", vp2, new_line('c')
            write(fstderr,*) "em = ", em, new_line('c')
            write(fstderr,*) "Kh_zm = ", Kh_zm, new_line('c')
            write(fstderr,*) "Kh_zt = ", Kh_zt, new_line('c')
            write(fstderr,*) "invrs_tau_C4 zm = ", invrs_tau_C4_zm, new_line('c')
            write(fstderr,*) "invrs_tau_wp3_zt = ", invrs_tau_wp3_zt, new_line('c')
            write(fstderr,*) "Skw_zm = ", Skw_zm, new_line('c')
            write(fstderr,*) "Skw_zt = ", Skw_zt, new_line('c')
            write(fstderr,*) "mixt_frac = ", mixt_frac, new_line('c')
            write(fstderr,*) "a3 = ", a3, new_line('c')
            write(fstderr,*) "a3_zt = ", a3_zt, new_line('c')
            write(fstderr,*) "wp3_on_wp2 = ", wp3_on_wp2, new_line('c')
            write(fstderr,*) "invrs_tau_C1_zm = ", invrs_tau_C1_zm, new_line('c')
            write(fstderr,*) "rho_ds_zm = ", rho_ds_zm, new_line('c')
            write(fstderr,*) "rho_ds_zt = ", rho_ds_zt, new_line('c')
            write(fstderr,*) "invrs_rho_ds_zm = ", invrs_rho_ds_zm, new_line('c')
            write(fstderr,*) "invrs_rho_ds_zt = ", invrs_rho_ds_zt, new_line('c')
            write(fstderr,*) "radf = ", radf, new_line('c')
            write(fstderr,*) "thv_ds_zm = ", thv_ds_zm, new_line('c')
            write(fstderr,*) "thv_ds_zt = ", thv_ds_zt, new_line('c')
            write(fstderr,*) "Cx_fnc_Richardson = ", Cx_fnc_Richardson, new_line('c')
            write(fstderr,*) "wp2_splat = ", wp2_splat, new_line('c')
            write(fstderr,*) "wp3_splat = ", wp3_splat, new_line('c')
            write(fstderr,*) "wprtp = ", wprtp, new_line('c')
            write(fstderr,*) "wpthlp = ", wpthlp, new_line('c')
            write(fstderr,*) "rtp2 = ", rtp2, new_line('c')
            write(fstderr,*) "thlp2 = ", thlp2, new_line('c')
            write(fstderr,*) "pdf_implicit_coefs_terms%coef_wp4_implicit = ", &
                             pdf_implicit_coefs_terms%coef_wp4_implicit
            write(fstderr,*) new_line('c')

            write(fstderr,*) "Intent(in/out)"

            write(fstderr,*) "wp2_zt = ", wp2_zt, new_line('c')
            write(fstderr,*) "wp3_zm = ", wp3_zm, new_line('c')
            if ( l_lmm_stepping ) &
               write(fstderr,*) "wp2 (pre-solve) = ", wp2_old, new_line('c')
            write(fstderr,*) "wp2 = ", wp2, new_line('c')
            if ( l_lmm_stepping ) &
               write(fstderr,*) "wp3 (pre-solve) = ", wp3_old, new_line('c')
            write(fstderr,*) "wp3 = ", wp3, new_line('c')

        end if ! fatal error
    end if

    return

  end subroutine advance_wp2_wp3

  !=============================================================================
  subroutine wp23_solve( gr, dt, sfc_elevation, sigma_sqd_w, wm_zm,               & ! Intent(in)
                         wm_zt, a3, a3_zt, wp3_on_wp2,                            & ! Intent(in)
                         wpup2, wpvp2, wp2up2, wp2vp2, wp4,                       & ! Intent(in)
                         wpthvp, wp2thvp, um, vm, upwp, vpwp,                     & ! Intent(in)
                         up2, vp2, em, Kw1, Kw8, Kh_zt, Skw_zt,                   & ! Intent(in)
                         invrs_tau_C4_zm, invrs_tau_wp3_zt,                       & ! Intent(in)
                         invrs_tau_C1_zm, C1_Skw_fnc,                             & ! Intent(in)
                         C11_Skw_fnc, rho_ds_zm,                                  & ! Intent(in)
                         rho_ds_zt, invrs_rho_ds_zm,                              & ! Intent(in)
                         invrs_rho_ds_zt, radf,                                   & ! Intent(in)
                         thv_ds_zm, thv_ds_zt,                                    & ! Intent(in)
                         wp2_splat, wp3_splat,                                    & ! Intent(in)
                         pdf_implicit_coefs_terms,                                & ! Intent(in)
                         wprtp, wpthlp, rtp2, thlp2,                              & ! Intent(in)
                         clubb_params, nu_vert_res_dep,                           & ! Intent(in)
                         iiPDF_type,                                              & ! Intent(in)
                         l_min_wp2_from_corr_wx,                                  & ! Intent(in)
                         l_upwind_xm_ma,                                          & ! Intent(in)
                         l_tke_aniso,                                             & ! Intent(in)
                         l_standard_term_ta,                                      & ! Intent(in)
                         l_partial_upwind_wp3,                                    & ! Intent(in)
                         l_damp_wp2_using_em,                                     & ! Intent(in)
                         l_damp_wp3_Skw_squared,                                  & ! Intent(in)
                         l_use_tke_in_wp3_pr_turb_term,                           & ! Intent(in)
                         l_use_tke_in_wp2_wp3_K_dfsn,                             & ! Intent(in)
                         stats_zt, stats_zm, stats_sfc,                           & ! intent(inout)
                         wp2, wp3, wp3_zm, wp2_zt )                                 ! Intent(inout)

    ! Description:
    ! Decompose, and back substitute the matrix for wp2/wp3

    ! References:
    ! _Equations for CLUBB_ section 6.3
    !------------------------------------------------------------------------

    use grid_class, only:  & 
        grid ! Type

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
        iiPDF_ADG1,                   & ! Variable(s)
        iiPDF_new,                    &
        iiPDF_new_hybrid,             &
        l_hole_fill,                  &
        l_explicit_turbulent_adv_wp3

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use lapack_wrap, only:  & 
        band_solve,  & ! Procedure(s) 
        band_solvex

    use parameter_indices, only: &
        nparams, & ! Variable(s)
        iSkw_max_mag

    use parameters_tunable, only: &
        nu_vertical_res_dep    ! Type(s)

    use fill_holes, only: & 
        fill_holes_vertical

    use clip_explicit, only: &
        clip_variance, & ! Procedure(s)
        clip_variance_level, &
        clip_skewness

    use pdf_parameter_module, only: &
        implicit_coefs_terms    ! Variable Type

    use stats_type_utilities, only: & 
        stat_begin_update, & ! Procedure(s)
        stat_update_var, &
        stat_update_var_pt, &
        stat_end_update, &
        stat_end_update_pt

    use stats_variables, only:  & 
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
        iwp3_pr_tp, &
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
        ztscr16, &
        ztscr17, &
        ztscr18

    use stats_type, only: stats ! Type

    implicit none

    type (stats), target, intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc

    type (grid), target, intent(in) :: gr

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
      wpup2,           & ! w'u'^2 (thermodynamic levels)             [m^3/s^3]
      wpvp2,           & ! w'v'^2 (thermodynamic levels)             [m^3/s^3]
      wp2up2,          & ! w'^2u'^2 (momentum levels)                [m^4/s^4]
      wp2vp2,          & ! w'^2v'^2 (momentum levels)                [m^4/s^4]
      wp4,             & ! w'^4 (momentum levels)                    [m^4/s^4]
      wpthvp,          & ! w'th_v' (momentum levels)                 [K m/s]
      wp2thvp,         & ! w'^2th_v' (thermodynamic levels)          [K m^2/s^2]
      um,              & ! u wind component (thermodynamic levels)   [m/s]
      vm,              & ! v wind component (thermodynamic levels)   [m/s]
      upwp,            & ! u'w' (momentum levels)                    [m^2/s^2]
      vpwp,            & ! v'w' (momentum levels)                    [m^2/s^2]
      up2,             & ! u'^2 (momentum levels)                    [m^2/s^2]
      vp2,             & ! v'^2 (momentum levels)                    [m^2/s^2]
      em,              & ! Turbulence kinetic energy                 [m^2/s^2]
      Kw1,             & ! Coefficient of eddy diffusivity for w'^2  [m^2/s]
      Kw8,             & ! Coefficient of eddy diffusivity for w'^3  [m^2/s]
      Kh_zt,           & ! Eddy diffusivity on thermodynamic levels  [m^2/s]
      Skw_zt,          & ! Skewness of w on thermodynamic levels     [-]
      invrs_tau_C4_zm, & ! Inverse time-scale tau on momentum levels         [1/s]
      invrs_tau_wp3_zt, & ! Inverse time-scale tau on thermodynamic levels    [1/s]
      invrs_tau_C1_zm, & ! Inverse tau values used for the C1 (dp1) term in wp2 [1/s]
      C1_Skw_fnc,      & ! C_1 parameter with Sk_w applied           [-]
      C11_Skw_fnc,     & ! C_11 parameter with Sk_w applied          [-]
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

    type(implicit_coefs_terms), intent(in) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    type(nu_vertical_res_dep), intent(in) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    integer, intent(in) :: &
      iiPDF_type    ! Selected option for the two-component normal (double
                    ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                    ! w, chi, and eta) portion of CLUBB's multivariate,
                    ! two-component PDF.

    logical, intent(in) :: &
      l_min_wp2_from_corr_wx,     & ! Flag to base the threshold minimum value of wp2 on keeping the
                                    ! overall correlation of w and x (w and rt, as well as w and
                                    ! theta-l) within the limits of -max_mag_correlation_flux to
                                    ! max_mag_correlation_flux.
      l_upwind_xm_ma,             & ! This flag determines whether we want to use an upwind
                                    ! differencing approximation rather than a centered differencing
                                    ! for turbulent or mean advection terms. It affects rtm, thlm,
                                    ! sclrm, um and vm.
      l_tke_aniso,                & ! For anisotropic turbulent kinetic energy, i.e. TKE = 1/2
                                    ! (u'^2 + v'^2 + w'^2)
      l_standard_term_ta,         & ! Use the standard discretization for the turbulent advection
                                    ! terms. Setting to .false. means that a_1 and a_3 are pulled
                                    ! outside of the derivative in advance_wp2_wp3_module.F90 and in
                                    ! advance_xp2_xpyp_module.F90.
      l_partial_upwind_wp3,       & ! Flag to use an "upwind" discretization rather
                                    ! than a centered discretization for the portion
                                    ! of the wp3 turbulent advection term for ADG1
                                    ! that is linearized in terms of wp3<t+1>.
                                    ! (Requires ADG1 PDF and l_standard_term_ta).
      l_damp_wp2_using_em,        & ! In wp2 equation, use a dissipation formula of 
                                    ! -(2/3)*em/tau_zm,
                                    ! as in Bougeault (1981)
      l_damp_wp3_Skw_squared,     & ! Set damping on wp3 to use Skw^2 rather than Skw^4
      l_use_tke_in_wp3_pr_turb_term, & ! Use TKE formulation for wp3 pr_turb term
      l_use_tke_in_wp2_wp3_K_dfsn   ! Use TKE in eddy diffusion for wp2 and wp3

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


    real( kind = core_rknd ), dimension(gr%nz) :: &
      threshold_array ! Minimum values for wp2 [m^2/s^2]

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

    if ( l_crank_nich_diff .and. l_use_tke_in_wp2_wp3_K_dfsn ) then
      write(fstderr,*) "The l_crank_nich_diff flag and l_use_tke_in_wp2_wp3_K_dfsn ", &
                       "flags cannot currently be used together."
      err_code = clubb_fatal_error
      return
    end if

    if ( .not. l_explicit_turbulent_adv_wp3 ) then

       if ( iiPDF_type == iiPDF_new .or. iiPDF_type == iiPDF_new_hybrid ) then

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
          coef_wp4_implicit = max( zt2zm( gr, coef_wp4_implicit_zt ), &
                                   zero_threshold )

          ! Set the value of coef_wp4_implicit to 0 at the lower boundary and at
          ! the upper boundary.  This sets the value of <w'^4> to 0 at the lower
          ! and upper boundaries.
          coef_wp4_implicit(1) = zero
          coef_wp4_implicit(gr%nz) = zero

          if ( l_stats_samp ) then
             call stat_update_var( icoef_wp4_implicit, coef_wp4_implicit, & ! intent(in)
                                   stats_zm )                               ! intent(inout)
          endif ! l_stats_samp

       elseif ( iiPDF_type == iiPDF_ADG1 ) then

          ! Define a_1 and a_3 (both are located on momentum levels).
          ! They are variables that are both functions of sigma_sqd_w (where
          ! sigma_sqd_w is located on momentum levels).
          a1 = one / ( one - sigma_sqd_w )

          ! Interpolate a_1 from momentum levels to thermodynamic levels.  This
          ! will be used for the w'^3 turbulent advection (ta) term.
          a1_zt = max( zm2zt( gr, a1 ), zero_threshold )  ! Positive def. quantity

       endif ! iiPDF_type

    endif ! .not. l_explicit_turbulent_adv_wp3

    ! Compute the explicit portion of the w'^2 and w'^3 equations.
    ! Build the right-hand side vector.
    call wp23_rhs( gr, dt, wp2, wp3, a1, a1_zt, a3, a3_zt, wp3_on_wp2, &             ! intent(in)
                   coef_wp4_implicit, wpup2, wpvp2, wp2up2, wp2vp2, wp4, &           ! intent(in)
                   wpthvp, wp2thvp, um, vm, &                                        ! intent(in)
                   upwp, vpwp, up2, vp2, em, Kw1, Kw8, Kh_zt,  &                     ! intent(in)
                   Skw_zt, invrs_tau_C4_zm, invrs_tau_wp3_zt, &                      ! intent(in)
                   invrs_tau_C1_zm, C1_Skw_fnc, &                                    ! intent(in)
                   C11_Skw_fnc, rho_ds_zm, rho_ds_zt, &                              ! intent(in)
                   invrs_rho_ds_zm, invrs_rho_ds_zt, radf, thv_ds_zm, thv_ds_zt, &   ! intent(in)
                   wp2_splat, wp3_splat, &                                           ! intent(in)
                   l_crank_nich_diff, &                                              ! intent(in)
                   clubb_params, nu_vert_res_dep, &                                  ! intent(in)
                   iiPDF_type, &                                                     ! intent(in)
                   l_tke_aniso, &                                                    ! intent(in)
                   l_standard_term_ta, &                                             ! intent(in)
                   l_partial_upwind_wp3, &                                           ! intent(in)
                   l_damp_wp2_using_em, &                                            ! intent(in)
                   l_damp_wp3_Skw_squared, &                                         ! intent(in)
                   l_use_tke_in_wp3_pr_turb_term, &                                  ! intent(in)
                   l_use_tke_in_wp2_wp3_K_dfsn, &                                    ! intent(in)
                   stats_zt, stats_zm, &                                             ! intent(inout)
                   rhs )                                                             ! intent(out)

    ! Save the value of rhs, which will be overwritten with the solution as
    ! part of the solving routine.
    rhs_save = rhs

    ! Compute the implicit portion of the w'^2 and w'^3 equations.
    ! Build the left-hand side matrix.
    call wp23_lhs( gr, dt, wp2, wm_zm, wm_zt, a1, a1_zt, a3, a3_zt,  & ! intent(in)
                   wp3_on_wp2, coef_wp4_implicit, & ! intent(in)
                   Kw1, Kw8, Skw_zt, & ! intent(in) 
                   invrs_tau_C4_zm, invrs_tau_wp3_zt, invrs_tau_C1_zm, C1_Skw_fnc, & !intent(in)
                   C11_Skw_fnc, rho_ds_zm, rho_ds_zt, & ! intent(in)
                   invrs_rho_ds_zm, invrs_rho_ds_zt, l_crank_nich_diff, & ! intent(in)
                   clubb_params, nu_vert_res_dep, & ! intent(in)
                   iiPDF_type, & ! intent(in)
                   l_upwind_xm_ma, & ! intent(in)
                   l_tke_aniso, & ! intent(in)
                   l_standard_term_ta, & ! intent(in)
                   l_partial_upwind_wp3, & ! intent(in)
                   l_damp_wp3_Skw_squared, & ! intent(in)
                   lhs, wp3_pr3_lhs ) ! intent(out)

    ! Solve the system with LAPACK
    if ( l_stats_samp .and. iwp23_matrix_condt_num > 0 ) then

        ! Perform LU decomp and solve system (LAPACK with diagnostics)
        ! Note that this can change the answer slightly
        call band_solvex( "wp2_wp3", 2, 2, 2*gr%nz, nrhs, & ! intent(in)
                          lhs, rhs, &                       ! intent(inout)
                          solut, rcond )                    ! intent(out)

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
        call stat_update_var_pt( iwp23_matrix_condt_num, 1, one / rcond, & ! intent(in) 
                                 stats_sfc )                               ! intent(inout)

    else

        ! Perform LU decomp and solve system (LAPACK)
        call band_solve( "wp2_wp3", 2, 2, 2*gr%nz, nrhs, & ! intent(in)
                         lhs, &                            ! intent(in)
                         rhs, &                            ! intent(inout)
                         solut )                           ! intent(out)

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
        call stat_end_update_pt( iwp2_dp1, k, & ! intent(in)
           zmscr01(k) * wp2(k), &               ! intent(in)
           stats_zm )                           ! intent(inout)

        ! w'^2 term dp2 has both implicit and explicit components (if the
        ! Crank-Nicholson scheme is selected or if l_use_tke_in_wp2_wp3_K_dfsn is true);
        ! call stat_end_update_pt.  
        ! If neither of these flags is true, then w'^2 term dp2 is 
        ! completely implicit; call stat_update_var_pt.
        if ( l_crank_nich_diff .or. l_use_tke_in_wp2_wp3_K_dfsn ) then
           call stat_end_update_pt( iwp2_dp2, k, & ! intent(in)
              zmscr02(k) * wp2(km1) & 
            + zmscr03(k) * wp2(k) & 
            + zmscr04(k) * wp2(kp1), &             ! intent(in)
              stats_zm )                           ! intent(inout)
        else
           call stat_update_var_pt( iwp2_dp2, k, & ! intent(in)
              zmscr02(k) * wp2(km1) & 
            + zmscr03(k) * wp2(k) & 
            + zmscr04(k) * wp2(kp1), &             ! intent(in)
              stats_zm )                           ! intent(inout)
        endif

        ! w'^2 term ta is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwp2_ta, k, & ! intent(in)
           zmscr05(k) * wp3(k) & 
         + zmscr06(k) * wp3(kp1), &            ! intent(in)
           stats_zm )                          ! intent(inout)

        ! w'^2 term ma is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwp2_ma, k, & ! intent(in)
           zmscr07(k) * wp2(km1) & 
         + zmscr08(k) * wp2(k) & 
         + zmscr09(k) * wp2(kp1), &            ! intent(in)
           stats_zm )                          ! intent(inout)

        ! w'^2 term ac is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwp2_ac, k,  & ! intent(in) 
           zmscr10(k) * wp2(k), &               ! intent(in)
           stats_zm )                           ! intent(inout)

        ! w'^2 term pr1 has both implicit and explicit components;
        ! call stat_end_update_pt.
        if ( l_tke_aniso ) then
          call stat_end_update_pt( iwp2_pr1, k, & ! intent(in)
             zmscr12(k) * wp2(k), &               ! intent(in)
             stats_zm )                           ! intent(inout)
        endif

        ! w'^2 term pr2 has both implicit and explicit components;
        ! call stat_end_update_pt.
        call stat_end_update_pt( iwp2_pr2, k, & ! intent(in) 
           zmscr11(k) * wp2(k), &               ! intent(in)
           stats_zm )                           ! intent(inout)

      enddo

      ! Finalize implicit contributions for wp3

      do k = 2, gr%nz-1, 1

        km1 = max( k-1, 1 )
        kp1 = min( k+1, gr%nz )

        ! w'^3 term pr1 has both implicit and explicit components; 
        ! call stat_end_update_pt.
        call stat_end_update_pt( iwp3_pr1, k, & ! intent(in) 
           ztscr01(k) * wp3(k), &               ! intent(in) 
           stats_zt )                           ! intent(inout)

        ! w'^3 term dp1 has both implicit and explicit components (if the
        ! Crank-Nicholson scheme is selected or l_use_tke_in_wp2_wp3_K_dfsn is true);
        ! call stat_end_update_pt.  
        ! If neither of these flags is true, then w'^3 term dp1 is 
        ! completely implicit; call stat_update_var_pt.
        if ( l_crank_nich_diff .or. l_use_tke_in_wp2_wp3_K_dfsn ) then
           call stat_end_update_pt( iwp3_dp1, k, & ! intent(in) 
              ztscr02(k) * wp3(km1) & 
            + ztscr03(k) * wp3(k) & 
            + ztscr04(k) * wp3(kp1), &             ! intent(in)
              stats_zt )                           ! intent(inout)
        else
           call stat_update_var_pt( iwp3_dp1, k, & ! intent(in)
              ztscr02(k) * wp3(km1) & 
            + ztscr03(k) * wp3(k) & 
            + ztscr04(k) * wp3(kp1), &             ! intent(in)
              stats_zt )                           ! intent(inout)
        endif

        ! w'^3 term ta has both implicit and explicit components; 
        ! call stat_end_update_pt.
        call stat_end_update_pt( iwp3_ta, k, & ! intent(in)
           ztscr05(k) * wp3(km1) & 
         + ztscr06(k) * wp2(km1) & 
         + ztscr07(k) * wp3(k) & 
         + ztscr08(k) * wp2(k) & 
         + ztscr09(k) * wp3(kp1), &            ! intent(in)
           stats_zt )                          ! intent(inout)

        ! w'^3 term tp has both implicit and explicit components; 
        ! call stat_end_update_pt.
        call stat_end_update_pt( iwp3_tp, k,  & ! intent(in)
           ztscr10(k) * wp2(km1) & 
         + ztscr11(k) * wp2(k), &               ! intent(in)
           stats_zt )                           ! intent(inout)

        ! w'^3 term pr_tp same as above tp term but opposite sign.
        call stat_end_update_pt( iwp3_pr_tp, k,  & ! intent(in)
           ztscr17(k) * wp2(km1) &
         + ztscr18(k) * wp2(k), & ! intent(in)
           stats_zt )             ! intent(inout)

        ! w'^3 pressure term 3 (pr3) has both implicit and explicit components;
        ! call stat_end_update_pt
        call stat_end_update_pt( iwp3_pr3, k, & ! intent(in)
         - wp3_pr3_lhs(k,5) * wp3(km1) &
         - wp3_pr3_lhs(k,4) * wp2(km1) &
         - wp3_pr3_lhs(k,3) * wp3(k) &
         - wp3_pr3_lhs(k,2) * wp2(k) &
         - wp3_pr3_lhs(k,1) * wp3(kp1), &       ! intent(in)
           stats_zt )                           ! intent(inout)

        ! w'^3 term ma is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwp3_ma, k, & ! intent(in)
           ztscr12(k) * wp3(km1) & 
         + ztscr13(k) * wp3(k) & 
         + ztscr14(k) * wp3(kp1), &            ! intent(in)
           stats_zt )                          ! intent(inout)

        ! w'^3 term ac is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwp3_ac, k, & ! intent(in) 
           ztscr15(k) * wp3(k), &              ! intent(in)
           stats_zt )                          ! intent(inout)

        ! w'^3 term pr2 has both implicit and explicit components; 
        ! call stat_end_update_pt.
        call stat_end_update_pt( iwp3_pr2, k, & ! intent(in) 
           ztscr16(k) * wp3(k), &               ! intent(in)
           stats_zt )                           ! intent(inout)

      enddo

    endif ! l_stats_samp


    if ( l_stats_samp ) then
      ! Store previous value for effect of the positive definite scheme
      call stat_begin_update( gr, iwp2_pd, wp2 / dt, & ! intent(in)
                              stats_zm )           ! intent(inout)
    endif

    if ( l_hole_fill .and. any( wp2 < w_tol_sqd ) ) then

      ! Use a simple hole filling algorithm
      call fill_holes_vertical( gr, 2, w_tol_sqd, "zm", &   ! intent(in)
                                rho_ds_zt, rho_ds_zm, & ! intent(in)
                                wp2 )                   ! intent(inout)

    endif ! wp2

    ! Here we attempt to clip extreme values of wp2 to prevent a crash of the
    ! type found on the Climate Process Team ticket #49.  Chris Golaz found that
    ! instability caused by large wp2 in CLUBB led unrealistic results in AM3.
    ! -dschanen 11 Apr 2011
    where ( wp2 > 1000._core_rknd ) wp2 = 1000._core_rknd

    if ( l_stats_samp ) then
      ! Store updated value for effect of the positive definite scheme
      call stat_end_update( gr, iwp2_pd, wp2 / dt, & ! intent(in)
                            stats_zm )           ! intent(inout)
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

          threshold_array(k) &
          = max( w_tol_sqd, &
                 wprtp(k)**2 / ( rtp2(k) * max_mag_correlation_flux**2 ), &
                 wpthlp(k)**2 / ( thlp2(k) * max_mag_correlation_flux**2 ) )

       enddo ! k = 1, gr%nz, 1

       call clip_variance_level( gr, clip_wp2, dt, threshold_array, & ! intent(in)
                                 stats_zm, &                   ! intent(inout)
                                 wp2 )                         ! intent(inout)

    else

       ! Consider only the minimum tolerance threshold value for wp2.
       threshold = w_tol_sqd

       call clip_variance( gr, clip_wp2, dt, threshold, & ! Intent(in)
                           stats_zm, &                ! intent(inout)
                           wp2 )                      ! Intent(inout)

    endif ! l_min_wp2_from_corr_wx

    ! Interpolate w'^2 from momentum levels to thermodynamic levels.
    ! This is used for the clipping of w'^3 according to the value
    ! of Sk_w now that w'^2 and w'^3 have been advanced one timestep.
    wp2_zt = max( zm2zt( gr, wp2 ), w_tol_sqd )   ! Positive definite quantity

    ! Clip w'^3 by limiting skewness.
    call clip_skewness( gr, dt, sfc_elevation, &              ! intent(in)
                        clubb_params(iSkw_max_mag), wp2_zt, & ! intent(in)
                        stats_zt, &                           ! intent(inout)
                        wp3 )                                 ! intent(inout)

    ! Compute wp3_zm for output purposes
    wp3_zm = zt2zm( gr, wp3 )

    return
  end subroutine wp23_solve

  !=================================================================================
  subroutine wp23_lhs( gr, dt, wp2, wm_zm, wm_zt, a1, a1_zt, a3, a3_zt,  &
                       wp3_on_wp2, coef_wp4_implicit, &
                       Kw1, Kw8, Skw_zt, &
                       invrs_tau_C4_zm, invrs_tau_wp3_zt, invrs_tau_C1_zm, C1_Skw_fnc, &
                       C11_Skw_fnc, rho_ds_zm, rho_ds_zt, &
                       invrs_rho_ds_zm, invrs_rho_ds_zt, l_crank_nich_diff, &
                       clubb_params, nu_vert_res_dep, &
                       iiPDF_type, &
                       l_upwind_xm_ma, &
                       l_tke_aniso, &
                       l_standard_term_ta, &
                       l_partial_upwind_wp3, &
                       l_damp_wp3_Skw_squared, &
                       lhs, wp3_pr3_lhs )
    ! Description:
    ! Compute LHS band diagonal matrix for w'^2 and w'^3.
    ! This subroutine computes the implicit portion 
    ! of the w'^2 and w'^3 equations.
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
        grid, & ! Type
        zm2zt, & ! Procedure(s)
        zt2zm

    use parameter_indices, only: &
        nparams, & ! Variable(s)
        iC4, &
        iC_uu_shr, & 
        iC8, & 
        iC8b, & 
        iC12, &
        iC_wp3_pr_tp

    use parameters_tunable, only: &
        nu_vertical_res_dep    ! Type(s)

    use constants_clubb, only:  & 
        one, & ! Constant(s)
        zero, &
        gamma_over_implicit_ts

    use model_flags, only: &
        iiPDF_ADG1,                   & ! Variable(s)
        iiPDF_new,                    &
        iiPDF_new_hybrid,             &
        l_explicit_turbulent_adv_wp3

    use diffusion, only: & 
        diffusion_zm_lhs, & ! Procedures
        diffusion_zt_lhs

    use mean_adv, only: & 
        term_ma_zm_lhs, & ! Procedures
        term_ma_zt_lhs

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
        ztscr16,    &
        ztscr17,    &
        ztscr18

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
        iwp3_pr_tp, & 
        iwp3_ma, & 
        iwp3_ac, & 
        iwp3_pr2, & 
        iwp3_pr1, & 
        iwp3_dp1

    use advance_helper_module, only: set_boundary_conditions_lhs ! Procedure(s)

    implicit none

    type (grid), target, intent(in) :: gr


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
      invrs_tau_C4_zm,   & ! Inverse time-scale tau on momentum levels         [1/s]
      invrs_tau_wp3_zt,  & ! Inverse time-scale tau on thermodynamic levels    [1/s]
      invrs_tau_C1_zm,   & ! Inverse tau values used for the C1 (dp1) term in wp2 [1/s]
      C1_Skw_fnc,        & ! C_1 parameter with Sk_w applied           [-]
      C11_Skw_fnc,       & ! C_11 parameter with Sk_w applied          [-]
      rho_ds_zm,         & ! Dry, static density on momentum levels    [kg/m^3]
      rho_ds_zt,         & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zm,   & ! Inv. dry, static density @ momentum levs. [m^3/kg]
      invrs_rho_ds_zt      ! Inv. dry, static density @ thermo. levs.  [m^3/kg]

    logical, intent(in) :: & 
      l_crank_nich_diff  ! Turns on/off Crank-Nicholson diffusion.

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    type(nu_vertical_res_dep), intent(in) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    integer, intent(in) :: &
      iiPDF_type    ! Selected option for the two-component normal (double
                    ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                    ! w, chi, and eta) portion of CLUBB's multivariate,
                    ! two-component PDF.

    logical, intent(in) :: &
      l_upwind_xm_ma,         & ! This flag determines whether we want to use an upwind
                                ! differencing approximation rather than a centered differencing
                                ! for turbulent or mean advection terms. It affects rtm, thlm,
                                ! sclrm, um and vm.
      l_tke_aniso,            & ! For anisotropic turbulent kinetic energy, i.e. TKE = 1/2
                                ! (u'^2 + v'^2 + w'^2)
      l_standard_term_ta,     & ! Use the standard discretization for the turbulent advection
                                ! terms. Setting to .false. means that a_1 and a_3 are pulled
                                ! outside of the derivative in advance_wp2_wp3_module.F90 and in
                                ! advance_xp2_xpyp_module.F90.
      l_partial_upwind_wp3,   & ! Flag to use an "upwind" discretization rather
                                ! than a centered discretization for the portion
                                ! of the wp3 turbulent advection term for ADG1
                                ! that is linearized in terms of wp3<t+1>.
                                ! (Requires ADG1 PDF and l_standard_term_ta).
      l_damp_wp3_Skw_squared    ! Set damping on wp3 to use Skw^2 rather than Skw^4

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
      lhs_ta_wp2,     & ! Turbulent advection terms for wp2
      lhs_ta_wp3,     & ! Turbulent advection terms for wp3
      lhs_tp_wp3,     & ! Turbulent production terms of w'^3
      lhs_adv_tp_wp3, & ! Turbulent production terms of w'^3 (for stats)
      lhs_pr_tp_wp3     ! Pressure scrambling terms for turbulent production of w'^3 (for stats)

    real( kind = core_rknd ), dimension(gr%nz) :: &
      lhs_ac_pr2_wp2, &   ! Accumulation terms of w'^2 and w'^2 pressure term 2
      lhs_ac_pr2_wp3, &   ! Accumulation terms of w'^3 and w'^3 pressure term 2
      lhs_dp1_wp2, &      ! Dissipation terms 1 for w'^2
      lhs_pr1_wp3, &      ! Dissipation terms 1 for w'^3
      lhs_pr1_wp2         ! Pressure term 1 for w'2

    real( kind = core_rknd) :: &
      invrs_dt        ! Inverse of dt, 1/dt, used for computational efficiency

    real( kind = core_rknd ), dimension(gr%nz) :: &
      Kw1_zm, &     ! Eddy diffusivity coefficient, momentum levels [m2/s]
      Kw8_zt        ! Eddy diffusivity coefficient, thermo. levels [m2/s]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      zero_vector    ! Vector of 0s

    real( kind = core_rknd ) ::  &
      C_uu_shr,    & ! CLUBB tunable parameter C_uu_shr
      C12,         & ! CLUBB tunable parameter C12
      C_wp3_pr_tp    ! CLUBB tunable parameter C_wp3_pr_tp

    !---------------------- Begin Code ----------------------


    ! Initialize arrays to 0 and calculate invrs_dt
    lhs = 0.0_core_rknd
    wp3_pr3_lhs = 0.0_core_rknd
    wp3_term_ta_lhs_result = 0.0_core_rknd
    invrs_dt = 1.0_core_rknd / dt

    lhs_tp_wp3 = zero
    lhs_adv_tp_wp3 = zero
    lhs_pr_tp_wp3 = zero

    C_uu_shr = clubb_params(iC_uu_shr)
    C12 = clubb_params(iC12)
    C_wp3_pr_tp = clubb_params(iC_wp3_pr_tp)

    ! Interpolate variables used for diffusion
    Kw1_zm = max( zt2zm( gr, Kw1 ), zero )
    Kw8_zt = max( zm2zt( gr, Kw8 ), zero )

    ! Calculated mean advection term for w'2
    call term_ma_zm_lhs( gr, wm_zm(:), gr%invrs_dzm(:), &                  ! intent(in)
                         lhs_ma_zm(:,:) )                              ! intent(out)


    ! Calculated mean advection term for w'3
    call term_ma_zt_lhs( gr, wm_zt(:), gr%invrs_dzt(:), gr%invrs_dzm(:), & ! intent(in)
                         l_upwind_xm_ma, &                             ! intent(in)
                         lhs_ma_zt(:,:) )                              ! intent(out)


    ! Calculate diffusion term for w'2 using a completely implicit time step
    call diffusion_zm_lhs( gr, Kw1(:), Kw1_zm(:), nu_vert_res_dep%nu1, & ! intent(in)
                           gr%invrs_dzt(:), gr%invrs_dzm(:), & ! intent(in)
                           invrs_rho_ds_zm(:), rho_ds_zt(:), & ! intent(in)
                           lhs_diff_zm(:,:) ) ! intent(out)


    ! Calculate diffusion term for w'3 using a completely implicit time step
    call diffusion_zt_lhs( gr, Kw8(:), Kw8_zt(:), nu_vert_res_dep%nu8, & ! intent(in)
                           gr%invrs_dzm(:), gr%invrs_dzt(:), & ! intent(in)
                           invrs_rho_ds_zt(:), rho_ds_zm(:), & ! intent(in)
                           lhs_diff_zt(:,:) ) ! intent(out)

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
    call wp2_term_ta_lhs( gr, rho_ds_zt(:), invrs_rho_ds_zm(:), gr%invrs_dzm(:), & ! intent(in)
                          lhs_ta_wp2(:,:) )                                    ! intent(out)

    ! Calculate accumulation terms of w'^2 and w'^2 pressure term 2
    call wp2_terms_ac_pr2_lhs( gr, C_uu_shr, wm_zt(:), gr%invrs_dzm(:), & ! intent(in)
                               lhs_ac_pr2_wp2(:) )                    ! intent(out)

    ! Calculate dissipation terms 1 for w'^2
    call wp2_term_dp1_lhs( gr, C1_Skw_fnc(:), invrs_tau_C1_zm(:), & ! intent(in)
                           lhs_dp1_wp2(:) )                     ! intent(out)

    ! Calculate turbulent production terms of w'^3
    call wp3_term_tp_lhs( gr, one, wp2(:), &         ! intent(in)
                          rho_ds_zm(:), &            ! intent(in)
                          invrs_rho_ds_zt(:), &      ! intent(in)
                          gr%invrs_dzt(:), &         ! intent(in)
                          lhs_adv_tp_wp3(:,:) )      ! intent(out)

    ! Calculate pressure damping of turbulent production of w'^3
    call wp3_term_tp_lhs( gr, -1*C_wp3_pr_tp, wp2(:),  & ! intent(in)
                          rho_ds_zm(:), &                ! intent(in)
                          invrs_rho_ds_zt(:), &          ! intent(in)
                          gr%invrs_dzt(:), &             ! intent(in)
                          lhs_pr_tp_wp3(:,:) )           ! intent(out)

    ! Sum contributions to turbulent production from standard term & damping
    lhs_tp_wp3 = lhs_adv_tp_wp3 + lhs_pr_tp_wp3

    ! Calculate accumulation terms of w'^3 and w'^3 pressure terms 2
    call wp3_terms_ac_pr2_lhs( gr, C11_Skw_fnc(:), wm_zm(:), gr%invrs_dzt(:), & ! intent(in)
                               lhs_ac_pr2_wp3(:) )                          ! intent(out)

    ! Calculate pressure terms 1 for w'^3
    call wp3_term_pr1_lhs( gr, clubb_params(iC8), clubb_params(iC8b), & ! intent(in)
                           invrs_tau_wp3_zt(:), Skw_zt(:), &            ! intent(in)
                           l_damp_wp3_Skw_squared, &                    ! intent(in)
                           lhs_pr1_wp3(:) )                             ! intent(out)

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
        call wp2_term_pr1_lhs( gr, clubb_params(iC4), invrs_tau_C4_zm(:), & ! intent(in)
                               lhs_pr1_wp2(:) )      ! intent(out)

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
            call wp3_term_ta_ADG1_lhs( gr, wp2, a1, a1_zt, a3, a3_zt, &        ! intent(in)
                                       wp3_on_wp2, rho_ds_zm, &            ! intent(in)
                                       rho_ds_zt, invrs_rho_ds_zt, &       ! intent(in)
                                       gr%invrs_dzt, l_standard_term_ta, & ! intent(in)
                                       l_partial_upwind_wp3, &             ! intent(in)
                                       wp3_term_ta_lhs_result )            ! intent(out)

        elseif ( iiPDF_type == iiPDF_new &
                 .or. iiPDF_type == iiPDF_new_hybrid ) then

            ! The new PDF or the new hybrid PDF is used.

            ! Calculate terms
            call wp3_term_ta_new_pdf_lhs( gr, coef_wp4_implicit(:), wp2(:), &     ! intent(in)
                                          rho_ds_zm(:), invrs_rho_ds_zt(:), & ! intent(in)
                                          gr%invrs_dzt(:), &                  ! intent(in)
                                          lhs_ta_wp3(:,:) )                   ! intent(out)

            ! Save terms in wp3_term_ta_lhs_result
            wp3_term_ta_lhs_result((/2,4/),:) = lhs_ta_wp3(:,:)

        endif

        ! Add terms to lhs
        do k = 2, gr%nz-1

            k_wp3 = 2*k - 1

            lhs(:,k_wp3) = lhs(:,k_wp3) + gamma_over_implicit_ts * wp3_term_ta_lhs_result(:,k)

        end do

    endif


    ! --------- Statistics output ---------
    if ( l_stats_samp ) then

        zero_vector = zero

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
                ztscr10(k) = - gamma_over_implicit_ts * lhs_adv_tp_wp3(2,k)
                ztscr11(k) = - gamma_over_implicit_ts * lhs_adv_tp_wp3(1,k)
            endif

            if ( iwp3_pr_tp > 0 ) then
                ztscr17(k) = - gamma_over_implicit_ts * lhs_pr_tp_wp3(2,k)
                ztscr18(k) = - gamma_over_implicit_ts * lhs_pr_tp_wp3(1,k)
            endif

            ! Mean advection for wp2
            if ( iwp3_ma > 0 ) then
                ztscr12(k) = - lhs_ma_zt(3,k)
                ztscr13(k) = - lhs_ma_zt(2,k)
                ztscr14(k) = - lhs_ma_zt(1,k)
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

        ! Note:  To find the contribution of w'^2 term ac, substitute 0 for the
        !        C_uu_shr input to function wp2_terms_ac_pr2_lhs.
        if ( iwp2_ac > 0 ) then
           call wp2_terms_ac_pr2_lhs( gr, zero, wm_zt, gr%invrs_dzm, & ! intent(in)
                                      zmscr10  )                   ! intent(out)
           zmscr10 = - zmscr10
        endif

        ! Note:  To find the contribution of w'^2 term pr2, add 1 to the
        !        C_uu_shr input to function wp2_terms_ac_pr2_lhs.
        if ( iwp2_pr2 > 0 ) then
           call wp2_terms_ac_pr2_lhs( gr, (one+C_uu_shr), wm_zt, gr%invrs_dzm, & ! intent(in)
                                       zmscr11 )                             ! intent(out)
           zmscr11 = - zmscr11
        endif

        ! Note:  To find the contribution of w'^3 term ac, substitute 0 for the
        !        C_ll skewness function input to function wp3_terms_ac_pr2_lhs.
        if ( iwp3_ac > 0 ) then
           call wp3_terms_ac_pr2_lhs( gr, zero_vector, wm_zm, gr%invrs_dzt, & ! intent(in)
                                      ztscr15 )                           ! intent(out)
           ztscr15 = - ztscr15
        endif

        ! Note:  To find the contribution of w'^3 term pr2, add 1 to the
        !        C_ll skewness function input to function wp3_terms_ac_pr2_lhs.
        if ( iwp3_pr2 > 0 ) then
           call wp3_terms_ac_pr2_lhs( gr, (one+C11_Skw_fnc), wm_zm, gr%invrs_dzt, & ! intent(in)
                                      ztscr16 )                                 ! intent(out)
           ztscr16 = - ztscr16
        endif

    end if

    return

  end subroutine wp23_lhs

  !=================================================================================
  subroutine wp23_rhs( gr, dt, wp2, wp3, a1, a1_zt, a3, a3_zt, wp3_on_wp2, &
                       coef_wp4_implicit, wpup2, wpvp2, wp2up2, wp2vp2, wp4, &
                       wpthvp, wp2thvp, um, vm, &
                       upwp, vpwp, up2, vp2, em, Kw1, Kw8, Kh_zt, & 
                       Skw_zt, invrs_tau_C4_zm, invrs_tau_wp3_zt, &
                       invrs_tau_C1_zm, C1_Skw_fnc, &
                       C11_Skw_fnc, rho_ds_zm, rho_ds_zt, &
                       invrs_rho_ds_zm, invrs_rho_ds_zt, radf, thv_ds_zm, thv_ds_zt, &
                       wp2_splat, wp3_splat, & 
                       l_crank_nich_diff, &
                       clubb_params, nu_vert_res_dep, &
                       iiPDF_type, &
                       l_tke_aniso, &
                       l_standard_term_ta, &
                       l_partial_upwind_wp3, &
                       l_damp_wp2_using_em, &
                       l_damp_wp3_Skw_squared, &
                       l_use_tke_in_wp3_pr_turb_term, &
                       l_use_tke_in_wp2_wp3_K_dfsn, &
                       stats_zt, stats_zm, &
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
        grid    ! Variable

    use grid_class, only:  & 
        ddzt, & ! Procedure
        zm2zt, & 
        zt2zm

    use parameter_indices, only: &
        nparams, & ! Variable(s)
        iC4, &
        iC_uu_shr, &
        iC_uu_buoy, & 
        iC8, & 
        iC8b, & 
        iC12, &
        iC_wp2_pr_dfsn, & 
        iC_wp3_pr_tp, &
        iC_wp3_pr_turb, &
        iC_wp3_pr_dfsn

    use parameters_tunable, only: &
        nu_vertical_res_dep    ! Type(s)

    use constants_clubb, only: & 
        w_tol_sqd,     & ! Variable(s)
        one,           &
        one_half,      &
        zero,          &
        gamma_over_implicit_ts

    use model_flags, only:  &
        iiPDF_ADG1,                   & ! Variable(s)
        iiPDF_new,                    &
        iiPDF_new_hybrid,             &
        l_explicit_turbulent_adv_wp3

    use diffusion, only: & 
        diffusion_zm_lhs,  & ! Procedures
        diffusion_zt_lhs

    use clubb_precision, only:  & 
        core_rknd ! Variable

    use stats_variables, only:  & 
        l_stats_samp, iwp2_dp1, iwp2_dp2, iwp2_bp,   & ! Variable(s)
        iwp2_pr1, iwp2_pr2, iwp2_pr3, iwp2_pr_dfsn, iwp2_splat, iwp3_splat, &
        iwp3_ta, & 
        iwp3_tp, iwp3_pr_tp, iwp3_bp1, iwp3_pr2, iwp3_pr1, iwp3_dp1, iwp3_pr_turb, &
        iwp3_pr_dfsn, iwp3_pr3
        
    use stats_type_utilities, only:  &
        stat_update_var_pt,  & ! Procedure(s)
        stat_begin_update_pt,  &
        stat_modify_pt

    use advance_helper_module, only: set_boundary_conditions_rhs


    use stats_type, only: stats ! Type

    implicit none

    type (stats), target, intent(inout) :: &
      stats_zt, &
      stats_zm

    type (grid), target, intent(in) :: gr

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
      wpup2,             & ! w'u'^2 (thermodynamic levels)             [m^3/s^3]
      wpvp2,             & ! w'v'^2 (thermodynamic levels)             [m^3/s^3]
      wp2up2,            & ! w'^2u'^2 (momentum levels)                [m^4/s^4]
      wp2vp2,            & ! w'^2v'^2 (momentum levels)                [m^4/s^4]
      wp4,               & ! w'^4 (momentum levels)                    [m^4/s^4]
      wpthvp,            & ! w'th_v' (momentum levels)                 [K m/s]
      wp2thvp,           & ! w'^2th_v' (thermodynamic levels)        [K m^2/s^2]
      um,                & ! u wind component (thermodynamic levels)   [m/s]
      vm,                & ! v wind component (thermodynamic levels)   [m/s]
      upwp,              & ! u'w' (momentum levels)                    [m^2/s^2]
      vpwp,              & ! v'w' (momentum levels)                    [m^2/s^2]
      up2,               & ! u'^2 (momentum levels)                    [m^2/s^2]
      vp2,               & ! v'^2 (momentum levels)                    [m^2/s^2]
      em,                & ! Turbulence kinetic energy (momentum levels)   [m^2/s^2]
      Kw1,               & ! Coefficient of eddy diffusivity for w'^2  [m^2/s]
      Kw8,               & ! Coefficient of eddy diffusivity for w'^3  [m^2/s]
      Kh_zt,             & ! Eddy diffusivity on thermodynamic levels  [m^2/s]
      Skw_zt,            & ! Skewness of w on thermodynamic levels     [-]
      invrs_tau_C4_zm,   & ! Inverse time-scale tau on momentum levels         [1/s]
      invrs_tau_wp3_zt,  & ! Inverse time-scale tau on thermodynamic levels    [1/s]
      invrs_tau_C1_zm,   & ! Inverse tau values used for the C1 (dp1) term in wp2 [1/s]
      C1_Skw_fnc,        & ! C_1 parameter with Sk_w applied           [-]
      C11_Skw_fnc,       & ! C_11 parameter with Sk_w applied          [-]
      rho_ds_zm,         & ! Dry, static density on momentum levels    [kg/m^3]
      rho_ds_zt,         & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zm,   & ! Inv. dry, static density @ mom. levs.     [m^3/kg]
      invrs_rho_ds_zt,   & ! Inv. dry, static density @ thermo. levs.  [m^3/kg]
      radf,              & ! Buoyancy production at the CL top         [m^2/s^3]
      thv_ds_zm,         & ! Dry, base-state theta_v on momentum levs. [K]
      thv_ds_zt,         & ! Dry, base-state theta_v on thermo. levs.  [K]
      wp2_splat,         & ! Tendency of <w'^2> due to vertical compression of eddies [m^2/s^3]
      wp3_splat            ! Tendency of <w'^3> due to vertical compression of eddies [m^3/s^4]

    logical, intent(in) :: & 
      l_crank_nich_diff   ! Turns on/off Crank-Nicholson diffusion.

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    type(nu_vertical_res_dep), intent(in) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    integer, intent(in) :: &
      iiPDF_type    ! Selected option for the two-component normal (double
                    ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                    ! w, chi, and eta) portion of CLUBB's multivariate,
                    ! two-component PDF.

    logical, intent(in) :: &
      l_tke_aniso,          &        ! For anisotropic turbulent kinetic energy, i.e. TKE = 1/2
                                     ! (u'^2 + v'^2 + w'^2)
      l_standard_term_ta,   &        ! Use the standard discretization for the 
                                     ! turbulent advection terms.
                                     ! Setting to .false. means that a_1 and a_3 are pulled outside 
                                     ! of the derivative in advance_wp2_wp3_module.F90 and in
                                     ! advance_xp2_xpyp_module.F90.
      l_partial_upwind_wp3, &        ! Flag to use an "upwind" discretization rather
                                     ! than a centered discretization for the portion
                                     ! of the wp3 turbulent advection term for ADG1
                                     ! that is linearized in terms of wp3<t+1>.
                                     ! (Requires ADG1 PDF and l_standard_term_ta).
      l_damp_wp2_using_em,  &        ! In wp2 equation, use a dissipation formula of 
                                     ! -(2/3)*em/tau_zm,
                                     ! as in Bougeault (1981)
      l_damp_wp3_Skw_squared,     &  ! Set damping on wp3 to use Skw^2 rather than Skw^4
      l_use_tke_in_wp3_pr_turb_term, &  ! Use TKE formulation for wp3 pr_turb term
      l_use_tke_in_wp2_wp3_K_dfsn    ! Use TKE in eddy diffusion for wp2 and wp3

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
      lhs_tp_wp3,      & ! Turbulent production terms of w'^3
      lhs_adv_tp_wp3,  & ! Turbulent production terms of w'^3 (for stats)
      lhs_pr_tp_wp3,   & ! Pressure scrambling terms for turbulent production of w'^3 (for stats) 
      lhs_ta_wp3         ! Turbulent advection terms for wp3

    real( kind = core_rknd ), dimension(gr%nz) :: &
      lhs_dp1_wp2, &          ! wp2 "over-implicit" dissipation term
      rhs_dp1_wp2, &          ! wp2 rhs dissipation term
      lhs_pr1_wp2, &          ! wp2 "over-implicit" pressure term 1
      rhs_pr1_wp2, &          ! wp2 rhs pressure term 1
      lhs_pr1_wp3, &          ! wp3 "over-implicit" pressure term 1
      rhs_pr1_wp3, &          ! wp3 rhs pressure term 1
      rhs_bp_pr2_wp2, &       ! wp2 bouyancy production and pressure term 2
      rhs_pr_dfsn_wp2, &      ! wp2 pressure diffusion term
      rhs_bp1_pr2_wp3, &      ! wp3 bouyancy production 1 and pressure term 2
      rhs_pr3_wp2, &          ! wp2 pressure term 3
      rhs_pr3_wp3, &          ! wp3 pressure term 3
      rhs_ta_wp3, &           ! wp3 turbulent advection term
      rhs_pr_turb_wp3, &      ! wp3 pressure-turbulence correlation term !--EXPERIMENTAL--!
      rhs_pr_dfsn_wp3         ! wp3 pressure diffusion term

    real( kind = core_rknd ), dimension(gr%nz) :: &
      rhs_bp_wp2, &  ! wp2 bouyancy production (stats only)
      rhs_pr2_wp2, & ! wp2 pressure term 2 (stats only)
      rhs_bp1_wp3, & ! wp3 bouyancy production 1 (stats only)
      rhs_pr2_wp3    ! wp3 pressure term 2 (stats only)
    
    real( kind = core_rknd ) :: &
      invrs_dt        ! Inverse of dt, 1/dt, used for computational efficiency

    real( kind = core_rknd ), dimension(gr%nz) :: &
      Kw1_zm, &    ! Eddy diffusivity coefficient, momentum levels [m2/s]
      Kw8_zt       ! Eddy diffusivity coefficient, thermo. levels [m2/s]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      zero_vector    ! Vector of 0s

    real( kind = core_rknd ) ::  &
      C4,          & ! CLUBB tunable parameter C4
      C_uu_shr,    & ! CLUBB tunable parameter C_uu_shr
      C_uu_buoy,   & ! CLUBB tunable parameter C_uu_buoy
      C8,          & ! CLUBB tunable parameter C8
      C8b,         & ! CLUBB tunable parameter C8b
      C12,         & ! CLUBB tunable parameter C12
      C_wp3_pr_tp    ! CLUBB tunable parameter C_wp3_pr_tp

    ! --------------- Begin Code ---------------

    ! Initialize arrays to 0 and calculate invers_dt
    invrs_dt = 1.0_core_rknd / dt
    rhs = 0.0_core_rknd
    wp3_term_ta_lhs_result = zero

    lhs_tp_wp3 = zero
    lhs_adv_tp_wp3 = zero
    lhs_pr_tp_wp3 = zero

    C_uu_shr = clubb_params(iC_uu_shr)
    C_uu_buoy = clubb_params(iC_uu_buoy)
    C8 = clubb_params(iC8)
    C8b = clubb_params(iC8b)
    C12 = clubb_params(iC12)
    C_wp3_pr_tp = clubb_params(iC_wp3_pr_tp)

    Kw1_zm = max( zt2zm( gr, Kw1 ), zero )
    Kw8_zt = max( zm2zt( gr, Kw8 ), zero )

    ! Experimental term from CLUBB TRAC ticket #411

      ! Compute the vertical derivative of the u and v winds
        dum_dz = ddzt( gr, um )
        dvm_dz = ddzt( gr, vm )

      ! Calculate term
      call wp3_term_pr_turb_rhs( gr, clubb_params(iC_wp3_pr_turb),   & ! intent(in)
                                 Kh_zt(:), wpthvp(:),                & ! intent(in)
                                 dum_dz(:), dvm_dz(:),               & ! intent(in)
                                 upwp(:), vpwp(:),                   & ! intent(in)
                                 thv_ds_zt(:), gr%invrs_dzt(:),      & ! intent(in)
                                 rho_ds_zm(:), invrs_rho_ds_zt(:),   & !intent(in)
                                 em(:), wp2(:),                      & ! intent(in)
                                 rhs_pr_turb_wp3(:),                 & ! intent(out)
                                 l_use_tke_in_wp3_pr_turb_term )       ! intent(in)

      call wp3_term_pr_dfsn_rhs( gr, clubb_params(iC_wp3_pr_dfsn),   & ! intent(in)
                                 rho_ds_zm(:), invrs_rho_ds_zt(:),   & ! intent(in)
                                 wp2up2(:), wp2vp2(:), wp4(:),       & ! intent(in)
                                 up2(:), vp2(:), wp2(:),             & ! intent(in)
                                 rhs_pr_dfsn_wp3(:) )                  ! intent(out)

      ! Add term
      do k = 2, gr%nz-1

          k_wp3 = 2*k - 1

          rhs(k_wp3) = rhs(k_wp3) + rhs_pr_turb_wp3(k) + rhs_pr_dfsn_wp3(k)

      end do


    call wp2_term_pr_dfsn_rhs( gr, clubb_params(iC_wp2_pr_dfsn), &
                               rho_ds_zt, invrs_rho_ds_zm, &
                               wpup2, wpvp2, wp3, &
                               rhs_pr_dfsn_wp2 )

    do k = 2, gr%nz-1

      k_wp2 = 2*k

      rhs(k_wp2) = rhs(k_wp2) + rhs_pr_dfsn_wp2(k)

    enddo

    ! These lines are for the diffusional term with a Crank-Nicholson
    ! time step.  They are not used for completely implicit diffusion.
    if ( l_crank_nich_diff ) then

        ! Calculate RHS eddy diffusion terms for w'2 and w'3
        
        call diffusion_zm_lhs( gr, Kw1(:), Kw1_zm(:), nu_vert_res_dep%nu1, & ! intent(in) 
                               gr%invrs_dzt(:), gr%invrs_dzm(:), & ! intent(in)
                               invrs_rho_ds_zm(:), rho_ds_zt(:), & ! intent(in)
                               rhs_diff_zm(:,:) )                  ! inetnt(out)

        call diffusion_zt_lhs( gr, Kw8(:), Kw8_zt(:), nu_vert_res_dep%nu8, & ! inetnt(in) 
                               gr%invrs_dzm(:), gr%invrs_dzt(:), & ! intent(in)
                               invrs_rho_ds_zt(:), rho_ds_zm(:), & ! intent(in)
                               rhs_diff_zt(:,:) )                  ! intent(out)
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
 
    ! This code block adds terms to the right-hand side so that TKE is being
    ! used in eddy diffusion instead of just wp2 or wp3.  For example, in the
    ! wp2 equation, if this flag is false, the eddy diffusion term would
    ! normally be completely implicit (hence no right-hand side contribution),
    ! and equal to +d/dz((K+nu)d/dz(wp2)).  With this flag set to true, the eddy
    ! diffusion term will be +d/dz((K+nu)d/dz(up2+vp2+wp2)), but the up2 and vp2
    ! parts are added on here as if they were right-hand side terms. For the wp3
    ! equation, with this flag false, the eddy diffusion term is
    ! +d/dz((K+nu)d/dz(wp3)), but with this flag true, it will be
    ! +d/dz((K+nu)d/dz(wpup2+wpvp2+wp3)).
    if ( l_use_tke_in_wp2_wp3_K_dfsn ) then

      ! This part handles the wp2 equation terms.
      call diffusion_zm_lhs( gr, Kw1(:), Kw1_zm(:), nu_vert_res_dep%nu1, & ! intent(in) 
                             gr%invrs_dzt(:), gr%invrs_dzm(:), &           ! intent(in)
                             invrs_rho_ds_zm(:), rho_ds_zt(:), &           ! intent(in)
                             rhs_diff_zm(:,:) )                            ! intent(out)

      do k = 2, gr%nz-1
        k_wp2 = 2*k
        rhs(k_wp2) = rhs(k_wp2) &
                     - rhs_diff_zm(3,k) * ( up2(k-1) + vp2(k-1) ) &
                     - rhs_diff_zm(2,k) * ( up2(k) + vp2(k) ) &
                     - rhs_diff_zm(1,k) * ( up2(k+1) + vp2(k+1) )
      end do

      ! This part handles the wp3 equation terms.
      call diffusion_zt_lhs( gr, Kw8(:), Kw8_zt(:), nu_vert_res_dep%nu8, & ! intent(in) 
                             gr%invrs_dzm(:), gr%invrs_dzt(:), &           ! intent(in)
                             invrs_rho_ds_zt(:), rho_ds_zm(:), &           ! intent(in)
                             rhs_diff_zt(:,:) )                            ! intent(out)

      do k = 2, gr%nz-1
        k_wp3 = 2*k - 1
        rhs(k_wp3) = rhs(k_wp3) &
                     - rhs_diff_zt(3,k) * ( wpup2(k-1) + wpvp2(k-1) ) &
                     - rhs_diff_zt(2,k) * ( wpup2(k) + wpvp2(k) ) &
                     - rhs_diff_zt(1,k) * ( wpup2(k+1) + wpvp2(k+1) )
      end do

    end if

    if ( l_tke_aniso ) then

        C4 = clubb_params(iC4)

        ! Calculate "over-implicit" pressure terms for w'2 and w'3

        call wp2_term_pr1_rhs( gr, C4, up2(:), vp2(:), invrs_tau_C4_zm(:), & ! intent(in)
                               rhs_pr1_wp2(:) )                      ! intent(out)

        ! Note:  An "over-implicit" weighted time step is applied to the  term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note below for w'^3 RHS
        !        turbulent advection (ta) term).
        call wp2_term_pr1_lhs( gr, C4, invrs_tau_C4_zm(:), & ! intent(in)
                               lhs_pr1_wp2(:) )      ! intent(out)

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
    call wp3_term_tp_lhs( gr, one, wp2(:),      & ! intent(in)
                          rho_ds_zm(:),         & ! intent(in)
                          invrs_rho_ds_zt(:),   & ! intent(in)
                          gr%invrs_dzt(:),      & ! intent(in)
                          lhs_adv_tp_wp3(:,:) )   ! intent(out)

    ! Calculate pressure damping of turbulent production of w'^3
    call wp3_term_tp_lhs( gr, -1*C_wp3_pr_tp, wp2(:),  & ! intent(in)
                          rho_ds_zm(:),                & ! intent(in)
                          invrs_rho_ds_zt(:),          & ! intent(in)
                          gr%invrs_dzt(:),             & ! intent(in)
                          lhs_pr_tp_wp3(:,:) )           ! intent(out)

    ! Sum contributions to turbulent production from standard term & damping
    lhs_tp_wp3 = lhs_adv_tp_wp3 + lhs_pr_tp_wp3

    ! Calculate pressure terms 1 for w'^3
    call wp3_term_pr1_lhs( gr, C8, C8b, invrs_tau_wp3_zt(:), Skw_zt(:), & ! intent(in)
                           l_damp_wp3_Skw_squared,              & ! intent(in)
                           lhs_pr1_wp3(:) )                       ! intent(out)

    ! Calculate dissipation terms 1 for w'^2
    call wp2_term_dp1_lhs( gr, C1_Skw_fnc(:), invrs_tau_C1_zm(:), & ! intent(in)
                           lhs_dp1_wp2(:) )                     ! intent(out)

    ! Calculate buoyancy production of w'^2 and w'^2 pressure term 2
    call wp2_terms_bp_pr2_rhs( gr, C_uu_buoy , thv_ds_zm(:), wpthvp(:), & ! intent(in)
                               rhs_bp_pr2_wp2(:) )                    ! intent(out)

    ! Calculate pressure terms 3 for w'^2
    call wp2_term_pr3_rhs( gr, C_uu_shr, C_uu_buoy, thv_ds_zm(:), wpthvp(:), upwp(:), & ! intent(in)
                           um(:), vpwp(:), vm(:), gr%invrs_dzm(:),                & ! intent(in)
                           rhs_pr3_wp2(:) )                                         ! intent(out)

    ! Calculate dissipation terms 1 for w'^2
    call wp2_term_dp1_rhs( gr, C1_Skw_fnc(:), invrs_tau_C1_zm(:), & ! intent(in)
                           w_tol_sqd, up2(:), vp2(:), & ! intent(in)
                           l_damp_wp2_using_em, & ! intent(in)
                           rhs_dp1_wp2(:) ) ! intent(out)

    ! Calculate buoyancy production of w'^3 and w'^3 pressure term 2
    call wp3_terms_bp1_pr2_rhs( gr, C11_Skw_fnc(:), thv_ds_zt(:), wp2thvp(:), & ! intent(in)
                                rhs_bp1_pr2_wp3(:) )                        ! intent(out)

    ! Calculate pressure terms 1 for w'^3
    call wp3_term_pr1_rhs( gr, C8, C8b, invrs_tau_wp3_zt(:), Skw_zt(:), wp3(:), & ! intent(in)
                           l_damp_wp3_Skw_squared, &                      ! intent(in)
                           rhs_pr1_wp3(:) )                               ! intent(out)


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

        call wp3_term_ta_explicit_rhs( gr, wp4(:),             & ! intent(in)
                                       rho_ds_zm(:),       & ! intent(in)
                                       invrs_rho_ds_zt(:), & ! intent(in)
                                       gr%invrs_dzt(:),    & ! intent(in)
                                       rhs_ta_wp3(:) )       ! intent(out)

        ! Add RHS turbulent advection (ta) terms
        do k = 2, gr%nz-1

            k_wp3 = 2*k - 1

            rhs(k_wp3) = rhs(k_wp3) + rhs_ta_wp3(k)

        end do

    else

        ! The turbulent advection term is being solved implicitly. See note above

        if ( iiPDF_type == iiPDF_ADG1 ) then

            ! The ADG1 PDF is used.
            call wp3_term_ta_ADG1_lhs( gr, wp2, a1, a1_zt, a3, a3_zt,        & ! intent(in)
                                       wp3_on_wp2, rho_ds_zm,            & ! intent(in)
                                       rho_ds_zt, invrs_rho_ds_zt,       & ! intent(in)
                                       gr%invrs_dzt, l_standard_term_ta, & ! intent(in)
                                       l_partial_upwind_wp3,             & ! intent(in)
                                       wp3_term_ta_lhs_result )            ! intent(out)

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

        elseif ( iiPDF_type == iiPDF_new &
                 .or. iiPDF_type == iiPDF_new_hybrid ) then

            ! The new PDF or the new hybrid PDF is used.

            ! Calculate terms
            call wp3_term_ta_new_pdf_lhs( gr, coef_wp4_implicit(:), wp2(:),     & ! intent(in)
                                          rho_ds_zm(:), invrs_rho_ds_zt(:), & ! intent(in)
                                          gr%invrs_dzt(:),                  & ! intent(in)
                                          lhs_ta_wp3(:,:) )                   ! intent(out)
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

      ! Not using pressure term, set to 0
      rhs_pr3_wp3 = zero

    

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

        zero_vector = zero

        ! w'^2 term bp is completely explicit; call stat_update_var_pt.
        ! Note:  To find the contribution of w'^2 term bp, substitute 0 for the
        !        C_uu_buoy input to function wp2_terms_bp_pr2_rhs.
        call wp2_terms_bp_pr2_rhs( gr, zero, thv_ds_zm(:), wpthvp(:), & ! intent(in)
                                   rhs_bp_wp2(:) )                  ! intent(out)

        ! w'^2 term pr2 has both implicit and explicit components; call
        ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
        ! subtracts the value sent in, reverse the sign on wp2_terms_bp_pr2_rhs.
        ! Note:  To find the contribution of w'^2 term pr2, add 1 to the
        !        C_uu_buoy input to function wp2_terms_bp_pr2_rhs.
        call wp2_terms_bp_pr2_rhs( gr, (one+C_uu_buoy), thv_ds_zm(:), wpthvp(:), & ! intent(in)
                                   rhs_pr2_wp2(:) )                            ! intent(out)

        if ( l_explicit_turbulent_adv_wp3 ) then

           ! The w'^3 turbulent advection term is being solved explicitly.
           !
           ! The turbulent advection stats code is still set up in two parts,
           ! so call stat_begin_update_pt.  The implicit portion of the stat,
           ! which has a value of 0, will still be called later.  Since
           ! stat_begin_update_pt automatically subtracts the value sent in,
           ! reverse the sign on the input value.
           call wp3_term_ta_explicit_rhs( gr, wp4(:),             & ! intent(in)
                                          rho_ds_zm(:),       & ! intent(in)
                                          invrs_rho_ds_zt(:), & ! intent(in)
                                          gr%invrs_dzt(:),    & ! intent(in)
                                          rhs_ta_wp3(:) )       ! intent(out)

        endif ! l_explicit_turbulent_adv_wp3

        ! w'^3 term bp is completely explicit; call stat_update_var_pt.
        ! Note:  To find the contribution of w'^3 term bp, substitute 0 for the
        !        C_11 skewness function input to function wp3_terms_bp1_pr2_rhs.
        call wp3_terms_bp1_pr2_rhs( gr, zero_vector, thv_ds_zt(:), wp2thvp(:), & ! intent(in) 
                                    rhs_bp1_wp3(:) )                         ! intent(out)

        ! w'^3 term pr2 has both implicit and explicit components; call
        ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
        ! subtracts the value sent in, reverse the sign on wp3_terms_bp1_pr2_rhs.
        ! Note:  To find the contribution of w'^3 term pr2, add 1 to the
        !        C_11 skewness function input to function wp3_terms_bp1_pr2_rhs.
        call wp3_terms_bp1_pr2_rhs( gr, ( one + C11_Skw_fnc(:) ), thv_ds_zt(:), & ! intent(in)
                                    wp2thvp(:), &                             ! intent(in)
                                    rhs_pr2_wp3(:) )                          ! intent(out)

        do k = 2, gr%nz-1
 
            ! ----------- w'2 -----------

            ! w'^2 term dp2 has both implicit and explicit components (if the
            ! Crank-Nicholson scheme is selected); call stat_begin_update_pt.  
            ! Since stat_begin_update_pt automatically subtracts the value sent in, 
            ! reverse the sign on right-hand side diffusion component.  If 
            ! Crank-Nicholson diffusion is not selected, the stat_begin_update_pt 
            ! will not be called.
            if ( l_crank_nich_diff ) then
              call stat_begin_update_pt( iwp2_dp2, k, & ! intent(in) 
                rhs_diff_zm(3,k) * wp2(k-1)           &  
              + rhs_diff_zm(2,k) * wp2(k)             & 
              + rhs_diff_zm(1,k) * wp2(k+1),          & ! intent(in)
                stats_zm )                              ! intent(inout)
            endif

            ! w'^2 term dp2 and w'^3 term dp1 have both implicit and explicit 
            ! components (if the l_use_tke_in_wp2_wp3_K_dfsn flag is true;
            ! call stat_begin_update_pt.
            if ( l_use_tke_in_wp2_wp3_K_dfsn ) then
              call stat_begin_update_pt( iwp2_dp2, k, &
                         + rhs_diff_zm(3,k) * ( up2(k-1) + vp2(k-1) )  &
                         + rhs_diff_zm(2,k) * ( up2(k)   + vp2(k)   )  &
                         + rhs_diff_zm(1,k) * ( up2(k+1) + vp2(k+1) ), &
                           stats_zm )
            endif

            ! w'^2 term bp is completely explicit; call stat_update_var_pt.
            ! Note:  To find the contribution of w'^2 term bp, substitute 0 for the
            !        C_uu_buoy input to function wp2_terms_bp_pr2_rhs.
            call stat_update_var_pt( iwp2_bp, k, rhs_bp_wp2(k), & ! intent(in)
                                     stats_zm )                   ! intent(inout)


            call stat_update_var_pt( iwp2_pr_dfsn, k, rhs_pr_dfsn_wp2(k), & ! intent(in)
                                     stats_zm )                             ! intent(out)


            ! Include effect of vertical compression of eddies in wp2 budget
            call stat_update_var_pt( iwp2_splat, k, wp2_splat(k), & ! intent(in)
                                     stats_zm )                     ! intent(inout)


            if ( l_tke_aniso ) then

                ! w'^2 term pr1 has both implicit and explicit components; call
                ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
                ! subtracts the value sent in, reverse the sign on wp2_term_pr1_rhs.
                call stat_begin_update_pt( iwp2_pr1, k, -rhs_pr1_wp2(k), & ! intent(in)
                                           stats_zm )                      ! intent(inout)

                ! Note:  An "over-implicit" weighted time step is applied to this
                !        term.  A weighting factor of greater than 1 may be used to
                !        make the term more numerically stable (see note below for
                !        w'^3 RHS turbulent advection (ta) term).
                call stat_modify_pt( iwp2_pr1, k,                      & ! intent(in)       
                                   + ( one - gamma_over_implicit_ts )  &
                                   * ( - lhs_pr1_wp2(k) * wp2(k) ),    & ! intent(in)
                                     stats_zm )                          ! intent(inout)
            endif

            ! w'^2 term pr2 has both implicit and explicit components; call
            ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
            ! subtracts the value sent in, reverse the sign on wp2_terms_bp_pr2_rhs.
            ! Note:  To find the contribution of w'^2 term pr2, add 1 to the
            !        C_uu_buoy input to function wp2_terms_bp_pr2_rhs.
            call stat_begin_update_pt( iwp2_pr2, k, -rhs_pr2_wp2(k), & ! intent(in)
                                       stats_zm )                      ! intent(inout)

            ! w'^2 term dp1 has both implicit and explicit components; call
            ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
            ! subtracts the value sent in, reverse the sign on wp2_term_dp1_rhs.
            call stat_begin_update_pt( iwp2_dp1, k, -rhs_dp1_wp2(k), & ! intent(in)
                                       stats_zm )                      ! intent(inout)


            ! Note:  An "over-implicit" weighted time step is applied to this term.
            !        A weighting factor of greater than 1 may be used to make the
            !        term more numerically stable (see note below for w'^3 RHS
            !        turbulent advection (ta) term).
            call stat_modify_pt( iwp2_dp1, k,                        & ! intent(in)
                                 + ( one - gamma_over_implicit_ts )  &
                                 * ( - lhs_dp1_wp2(k) * wp2(k) ),    & ! intent(in)
                                 stats_zm )                            ! intent(inout)

            ! w'^2 term pr3 is completely explicit; call stat_update_var_pt.
            call stat_update_var_pt( iwp2_pr3, k, rhs_pr3_wp2(k), & ! intent(in)
                                     stats_zm )                     ! intent(inout)


            ! ----------- w'3 -----------

            if ( l_explicit_turbulent_adv_wp3 ) then !l_explicit_turbulent_adv_wp3

                ! The turbulent advection term is being solved explicitly.
                ! 
                ! The turbulent advection stats code is still set up in two parts,
                ! so call stat_begin_update_pt.  The implicit portion of the stat,
                ! which has a value of 0, will still be called later.  Since
                ! stat_begin_update_pt automatically subtracts the value sent in,
                ! reverse the sign on the input value.
                call stat_begin_update_pt( iwp3_ta, k, -rhs_ta_wp3(k), & ! intent(in)
                                           stats_zt )                    ! intent(inout)
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

                    call stat_begin_update_pt( iwp3_ta, k, & ! intent(in)
                                                - ( one - gamma_over_implicit_ts ) & ! intent(in)
                                                * ( - wp3_term_ta_lhs_result(1,k) * wp3(k+1) &
                                                    - wp3_term_ta_lhs_result(2,k) * wp2(k) &
                                                    - wp3_term_ta_lhs_result(3,k) * wp3(k) &
                                                    - wp3_term_ta_lhs_result(4,k) * wp2(k-1) &
                                                    - wp3_term_ta_lhs_result(5,k) * wp3(k-1) ), &
                                               stats_zt ) ! intent(inout)

                elseif ( iiPDF_type == iiPDF_new &
                         .or. iiPDF_type == iiPDF_new_hybrid ) then

                    ! The new PDF or the new hybrid PDF is used.

                    call stat_begin_update_pt( iwp3_ta, k, & ! intent(in)
                                               - ( one - gamma_over_implicit_ts )  &
                                                 * ( - lhs_ta_wp3(1,k) * wp2(k)  &
                                                     - lhs_ta_wp3(2,k) * wp2(k-1) ), & ! intent(in)
                                               stats_zt ) ! intent(inout)
                endif

            endif

            ! Note:  An "over-implicit" weighted time step is applied to this term.
            !        A weighting factor of greater than 1 may be used to make the
            !        term more numerically stable (see note above for RHS turbulent
            !        production (tp) term).  Call stat_begin_update_pt.  Since
            !        stat_begin_update_pt automatically subtracts the value sent in,
            !        reverse the sign on the input value.
            call stat_begin_update_pt( iwp3_tp, k, & ! intent(in)
                                       - ( one - gamma_over_implicit_ts )  &
                                         * ( - lhs_adv_tp_wp3(1,k) * wp2(k)  &
                                             - lhs_adv_tp_wp3(2,k) * wp2(k-1) ), & ! intent(in)
                                       stats_zt ) ! intent(inout)

            call stat_begin_update_pt( iwp3_pr_tp, k, & ! intent(in)
                                       - ( one - gamma_over_implicit_ts )  &
                                         * ( - lhs_pr_tp_wp3(1,k) * wp2(k)  &
                                             - lhs_pr_tp_wp3(2,k) * wp2(k-1) ), & ! intent(in)
                                       stats_zt ) ! intent(inout)


            ! w'^3 pressure term 3 (pr3) explicit (rhs) contribution
            call stat_begin_update_pt( iwp3_pr3, k, rhs_pr3_wp3(k), & ! intent(in)
                                       stats_zt )                     ! intent(inout)


            ! w'^3 term bp is completely explicit; call stat_update_var_pt.
            ! Note:  To find the contribution of w'^3 term bp, substitute 0 for the
            !        C_11 skewness function input to function wp3_terms_bp1_pr2_rhs.
            call stat_update_var_pt( iwp3_bp1, k, rhs_bp1_wp3(k), & ! intent(in)
                                     stats_zt )                     ! intent(inout)


            ! w'^3 term pr2 has both implicit and explicit components; call
            ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
            ! subtracts the value sent in, reverse the sign on wp3_terms_bp1_pr2_rhs.
            ! Note:  To find the contribution of w'^3 term pr2, add 1 to the
            !        C_11 skewness function input to function wp3_terms_bp1_pr2_rhs.
            call stat_begin_update_pt( iwp3_pr2, k, -rhs_pr2_wp3(k), & ! intent(in)
                                       stats_zt )                      ! intent(inout)

            ! w'^3 term pr1 has both implicit and explicit components; call 
            ! stat_begin_update_pt.  Since stat_begin_update_pt automatically 
            ! subtracts the value sent in, reverse the sign on wp3_term_pr1_rhs.
            call stat_begin_update_pt( iwp3_pr1, k, -rhs_pr1_wp3(k), & ! intent(in)
                                       stats_zt )                      ! intent(inout)


            ! Note:  An "over-implicit" weighted time step is applied to this term.
            !        A weighting factor of greater than 1 may be used to make the
            !        term more numerically stable (see note above for RHS turbulent
            !        advection (ta) term).
            call stat_modify_pt( iwp3_pr1, k,                        & ! intent(in)
                                 + ( one - gamma_over_implicit_ts )  &
                                 * ( - lhs_pr1_wp3(k) * wp3(k) ),    & ! intent(in)
                                 stats_zt )                            ! intent(inout)


            ! Include effect of vertical compression of eddies in wp2 budget
            call stat_update_var_pt( iwp3_splat, k, wp3_splat(k), & ! intent(in)
                                     stats_zt )                     ! intent(inout)


            if ( l_crank_nich_diff ) then

                ! w'^3 term dp1 has both implicit and explicit components (if the
                ! Crank-Nicholson scheme is selected); call stat_begin_update_pt.  
                ! Since stat_begin_update_pt automatically subtracts the value sent in, 
                ! reverse the sign on right-hand side diffusion component.  If 
                ! Crank-Nicholson diffusion is not selected, the stat_begin_update_pt 
                ! will not be called.
                call stat_begin_update_pt( iwp3_dp1, k, & ! intent(in) 
                                           rhs_diff_zt(3,k) * wp3(k-1) & 
                                         + rhs_diff_zt(2,k) * wp3(k) & 
                                         + rhs_diff_zt(1,k) * wp3(k+1), & ! intent(in)
                                           stats_zt ) ! intent(inout)
            endif

            ! w'^2 term dp2 and w'^3 term dp1 have both implicit and explicit 
            ! components (if the l_use_tke_in_wp2_wp3_K_dfsn flag is true;
            ! call stat_begin_update_pt.
            if ( l_use_tke_in_wp2_wp3_K_dfsn ) then
              call stat_begin_update_pt( iwp3_dp1, k, &
                         + rhs_diff_zt(3,k) * ( wpup2(k-1) + wpvp2(k-1) ) &
                         + rhs_diff_zt(2,k) * ( wpup2(k)   + wpvp2(k)   ) &
                         + rhs_diff_zt(1,k) * ( wpup2(k+1) + wpvp2(k+1) ), &
                           stats_zt )
            endif
                      
            ! Experimental bouyancy term
              call stat_update_var_pt( iwp3_pr_turb, k, rhs_pr_turb_wp3(k), & ! intent(in)
                                       stats_zt )                             ! intent(inout)
              call stat_update_var_pt( iwp3_pr_dfsn, k, rhs_pr_dfsn_wp3(k), & ! intent(in)
                                       stats_zt )                             ! intent(inout)
        end do

    endif

    return

  end subroutine wp23_rhs

  !=============================================================================
  pure subroutine wp2_term_ta_lhs( gr, rho_ds_zt, invrs_rho_ds_zm, invrs_dzm, &
                                   lhs_ta_wp2 )

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
    ! -----rho_ds_zt----------wp3------------------------------ t(k+1)
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

    use constants_clubb, only: &
        zero    ! Constant(s)

    use grid_class, only: & 
        grid ! Type

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2       ! Thermodynamic subdiagonal index.

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      rho_ds_zt,       & ! Dry, static density at thermodynamic levels  [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density at momentum levels  [m^3/kg]
      invrs_dzm          ! Inverse of grid spacing                      [1/m]

    ! Return Variable
    real( kind = core_rknd ), dimension(2,gr%nz), intent(out) :: &
      lhs_ta_wp2    ! LHS coefficient of wp2 turbulent advection  [1/m]

    ! Local variables
    integer :: k    ! Vertical level index


    ! Set lower boundary to 0
    lhs_ta_wp2(kp1_tdiag,1) = zero
    lhs_ta_wp2(k_tdiag,1) = zero

    ! Calculate term at all interior grid levels.
    do k = 2, gr%nz-1 

       ! Thermodynamic superdiagonal: [ x wp3(k+1,<t+1>) ]
       lhs_ta_wp2(kp1_tdiag,k) &
       = + invrs_rho_ds_zm(k) * invrs_dzm(k) * rho_ds_zt(k+1)

       ! Thermodynamic subdiagonal: [ x wp3(k,<t+1>) ]
       lhs_ta_wp2(k_tdiag,k) &
       = - invrs_rho_ds_zm(k) * invrs_dzm(k) * rho_ds_zt(k)

    enddo

    ! Set upper boundary to 0
    lhs_ta_wp2(kp1_tdiag,gr%nz) = zero
    lhs_ta_wp2(k_tdiag,gr%nz) = zero


    return

  end subroutine wp2_term_ta_lhs

  !=============================================================================
  pure subroutine wp2_terms_ac_pr2_lhs( gr, C_uu_shr, wm_zt, invrs_dzm, &
                                        lhs_ac_pr2_wp2 )

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
    ! + ( 1 - C_uu_shr ) ( -2 w'^2(t+1) dw/dz ).
    !
    ! Note 1:  When the term is brought over to the left-hand side, the sign 
    !          is reversed and the leading "-" in front of the "2" is changed 
    !          to a "+".
    ! Note 2:  We have broken C5 up into C_uu_shr for this term
    !          and C_uu_buoy for the buoyancy term.
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
    ! -------wm_zt--------------------------------------------- t(k+1)
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

    use grid_class, only: &
        grid ! Type

    use constants_clubb, only: &
        two,  & ! Variable(s)
        one,  &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      wm_zt,     & ! w wind component at thermodynamic levels    [m/s]
      invrs_dzm    ! Inverse of grid spacing                     [1/m]

    real( kind = core_rknd ), intent(in) :: & 
      C_uu_shr    ! Model parameter C_uu_shr                       [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      lhs_ac_pr2_wp2    ! LHS coefficient of wp2 ac and pr2 terms [1/s]

    ! Local Variables
    integer :: k ! Vertical level index


    ! Set lower boundary to 0
    lhs_ac_pr2_wp2(1) = zero

    ! Calculate term at all interior grid levels.
    do k = 2, gr%nz-1

       ! Momentum main diagonal: [ x wp2(k,<t+1>) ]
       lhs_ac_pr2_wp2(k) &
       = + ( one - C_uu_shr ) * two * invrs_dzm(k) * ( wm_zt(k+1) - wm_zt(k) )

    enddo ! k = 2, gr%nz-1

    ! Set upper boundary to 0
    lhs_ac_pr2_wp2(gr%nz) = zero


    return

  end subroutine wp2_terms_ac_pr2_lhs

  !=============================================================================
  pure subroutine wp2_term_dp1_lhs( gr, C1_Skw_fnc, invrs_tau1m, &
                                    lhs_dp1_wp2 )

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

    use grid_class, only:  & 
        grid ! Type

    use constants_clubb, only: &
        zero    ! Constant(s) 

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      C1_Skw_fnc, & ! C_1 parameter with Sk_w applied    [-]
      invrs_tau1m   ! Inverse time-scale tau at momentum levels  [1/s]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      lhs_dp1_wp2    ! LHS coefficient of wp2 dissipation term 1  [1/s]

    ! Local Variable
    integer :: k    ! Vertical level index


    ! Set lower boundary to 0
    lhs_dp1_wp2(1) = zero

    ! Calculate term at all interior grid levels.
    do k = 2, gr%nz-1

       ! Momentum main diagonal: [ x wp2(k,<t+1>) ]
       lhs_dp1_wp2(k) = + C1_Skw_fnc(k) * invrs_tau1m(k)

    enddo ! k = 2, gr%nz-1

    ! Set upper boundary to 0
    lhs_dp1_wp2(gr%nz) = zero


    return

  end subroutine wp2_term_dp1_lhs

  !=============================================================================
  pure subroutine wp2_term_pr1_lhs( gr, C4, invrs_tau_C4_zm, &
                                    lhs_pr1_wp2 )

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
    ! The values of w'^2 are found on momentum levels, as are the values of
    ! tau1m.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only:  &
        grid ! Type

    use constants_clubb, only: &
        three, & ! Variable(s)
        two,   &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      invrs_tau_C4_zm    ! Inverse time-scale tau at momentum levels  [1/s]

    real( kind = core_rknd ), intent(in) :: & 
      C4    ! Model parameter C_4                [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      lhs_pr1_wp2    ! LHS coefficient of wp2 pressure term 1  [1/s]
    
    ! Local Variables
    integer :: k


    ! Set lower boundary to 0
    lhs_pr1_wp2(1) = zero

    ! Calculate term at all interior grid levels.
    do k = 2, gr%nz-1

        ! Momentum main diagonal: [ x wp2(k,<t+1>) ]
        lhs_pr1_wp2(k) = + ( two * C4 * invrs_tau_C4_zm(k) ) / three
    
    enddo ! k = 2, gr%nz-1

    ! Set upper boundary to 0
    lhs_pr1_wp2(gr%nz) = zero


    return

  end subroutine wp2_term_pr1_lhs

  !=============================================================================
  pure subroutine wp2_terms_bp_pr2_rhs( gr, C_uu_buoy, thv_ds_zm, wpthvp, &
                                        rhs_bp_pr2_wp2 )

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
    ! + ( 1 - C_uu_buoy ) ( 2 (g/thv_ds) w'th_v' ).
    !
    ! Note:  We have broken C5 up into C_uu_shr for the accumulation term
    !        and C_uu_buoy for the buoyancy term.
    !
    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid ! Type

    use constants_clubb, only:  & ! Variable(s)        
        grav, & ! Gravitational acceleration [m/s^2]
        two,  &
        one,  &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C_uu_buoy    ! Model parameter C_uu_buoy                         [-]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      thv_ds_zm, & ! Dry, base-state theta_v at momentum levels   [K]
      wpthvp       ! w'th_v'                                      [K m/s]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      rhs_bp_pr2_wp2    ! RHS portion of wp2 from terms bp and pr2  [m^2/s^3]

    ! Local variables
    integer :: k    ! Vertical level index


    ! Set lower boundary to 0
    rhs_bp_pr2_wp2(1) = zero

    ! Calculate term at all interior grid levels.
    do k = 2, gr%nz-1

       rhs_bp_pr2_wp2(k) &
       = + ( one - C_uu_buoy ) * two * ( grav / thv_ds_zm(k) ) * wpthvp(k)

    enddo ! k = 2, gr%nz-1

    ! Set upper boundary to 0
    rhs_bp_pr2_wp2(gr%nz) = zero


    return

  end subroutine wp2_terms_bp_pr2_rhs

  !=============================================================================
  pure subroutine wp2_term_dp1_rhs( gr, C1_Skw_fnc, invrs_tau1m, &
                                    threshold, up2, vp2, &
                                    l_damp_wp2_using_em, &
                                    rhs_dp1_wp2 )

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

    use grid_class, only: &
        grid ! Type

    use constants_clubb, only: &
        zero    ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      C1_Skw_fnc,  & ! C_1 parameter with Sk_w applied                  [-]
      invrs_tau1m, & ! Inverse time-scale tau at momentum levels        [1/s]
      up2,         & ! Horizontal (east-west) velocity variance, u'^2   [m^2/s^2]
      vp2            ! Horizontal (north-south) velocity variance, v'^2 [m^2/s^2]

    real( kind = core_rknd ), intent(in) :: & 
      threshold    ! Minimum allowable value of w'^2       [m^2/s^2]

    logical, intent(in) :: &
      l_damp_wp2_using_em ! In wp2 equation, use a dissipation formula of -(2/3)*em/tau_zm,
                          ! as in Bougeault (1981)

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      rhs_dp1_wp2    ! RHS portion of wp2 from dissipation term 1  [m^2/s^3]

    ! Local variables
    integer :: k    ! Vertical level index


    ! Set lower boundary to 0
    rhs_dp1_wp2(1) = zero

    ! Calculate term at all interior grid levels.
    if ( l_damp_wp2_using_em ) then

       do k = 2, gr%nz-1
          rhs_dp1_wp2(k) = - ( C1_Skw_fnc(k) * invrs_tau1m(k) ) * ( up2(k) + vp2(k) )
       enddo ! k = 2, gr%nz-1

    else

       do k = 2, gr%nz-1
          rhs_dp1_wp2(k) = + ( C1_Skw_fnc(k) * invrs_tau1m(k) ) * threshold
       enddo ! k = 2, gr%nz-1

    endif ! l_damp_wp2_using_em

    ! Set upper boundary to 0
    rhs_dp1_wp2(gr%nz) = zero


    return

  end subroutine wp2_term_dp1_rhs

  !=============================================================================
  pure subroutine wp2_term_pr3_rhs( gr, C_uu_shr, C_uu_buoy, thv_ds_zm, wpthvp, upwp, &
                                    um, vpwp, vm, invrs_dzm, &
                                    rhs_pr3_wp2 )

    ! Description:
    ! Pressure term 3 for w'^2:  explicit portion of the code.
    !
    ! The d(w'^2)/dt equation contains pressure term 3:
    !
    ! + (2/3) C_5 [ (g/thv_ds) w'th_v' - u'w' du/dz - v'w' dv/dz ].
    !
    ! Note that below we have broken up C5 into C_uu_shr for shear terms and 
    ! C_uu_buoy for buoyancy terms.
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
    ! -----um--------------vm---------------------------------------- t(k+1)
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

    use grid_class, only: &
        grid ! Type

    use constants_clubb, only: & ! Variables 
        grav,           & ! Gravitational acceleration [m/s^2]
        two_thirds,     &
        zero,           &
        zero_threshold

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C_uu_shr,  & ! Model parameter                            [-]
      C_uu_buoy    ! Model parameter                            [-]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      thv_ds_zm, & ! Dry, base-state theta_v at momentum level (k)  [K]
      wpthvp,    & ! w'th_v'(k)                                     [K m/s]
      upwp,      & ! u'w'(k)                                        [m^2/s^2]
      um,        & ! um(k)                                          [m/s]
      vpwp,      & ! v'w'(k)                                        [m^2/s^2]
      vm,        & ! vm(k)                                          [m/s]
      invrs_dzm    ! Inverse of grid spacing (k)                    [1/m]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      rhs_pr3_wp2    ! RHS portion of wp2 from pressure term 3  [m^2/s^3]

    ! Local variables
    integer :: k    ! Vertical level index


    ! Set lower boundary to 0
    rhs_pr3_wp2(1) = zero

    ! Calculate term at all interior grid levels.
    do k = 2, gr%nz-1

       rhs_pr3_wp2(k) &
       ! Michael Falk, 2 August 2007
       ! Use the following code for standard mixing, with c_k=0.548:
       = + two_thirds * &
                      ( C_uu_buoy &
                        * ( grav / thv_ds_zm(k) ) * wpthvp(k) &
                      + C_uu_shr &
                        * ( - upwp(k) * invrs_dzm(k) * ( um(k+1) - um(k) ) &
                            - vpwp(k) * invrs_dzm(k) * ( vm(k+1) - vm(k) ) &
                          ) &
                      )
        ! Use the following code for alternate mixing, with c_k=0.1 or 0.2
!       = + two_thirds * C_uu_shr &
!                      * ( ( grav / thv_ds_zm(k) ) * wpthvp(k) &
!                          - 0. * upwp(k) * invrs_dzm(k) * ( um(k+1) - um(k) ) &
!                          - 0. * vpwp(k) * invrs_dzm(k) * ( vm(k+1) - vm(k) ) &
!                        )
!       eMFc

    enddo ! k = 2, gr%nz-1

    ! Set upper boundary to 0
    rhs_pr3_wp2(gr%nz) = zero

    ! Added by dschanen for ticket #36
    ! We have found that when shear generation is zero this term will only be
    ! offset by hole-filling (wp2_pd) and reduces turbulence 
    ! unrealistically at lower altitudes to make up the difference.
    rhs_pr3_wp2 = max( rhs_pr3_wp2, zero_threshold )


    return

  end subroutine wp2_term_pr3_rhs

  !=============================================================================
  pure subroutine wp2_term_pr1_rhs( gr, C4, up2, vp2, invrs_tau_C4_zm, &
                                    rhs_pr1_wp2 )

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

    use grid_class, only: &
        grid ! Type

    use constants_clubb, only: &
        three, & ! Constant9(s)
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C4    ! Model parameter C_4                      [-]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      up2,        & ! u'^2(k)                               [m^2/s^2]
      vp2,        & ! v'^2(k)                               [m^2/s^2]
      invrs_tau_C4_zm   ! Inverse time-scale tau at momentum levels [1/s]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      rhs_pr1_wp2    ! RHS portion of wp2 from pressure term 1  [m^2/s^3]

    ! Local Variables
    integer :: k    ! Vertical level index


    ! Set lower bounadry to 0
    rhs_pr1_wp2(1) = zero

    ! Calculate term at all interior grid levels.
    do k = 2, gr%nz-1

      rhs_pr1_wp2(k) = + ( C4 * ( up2(k) + vp2(k) ) * invrs_tau_C4_zm(k) ) / three

    enddo ! k = 2, gr%nz-1

    ! Set upper boundary to 0
    rhs_pr1_wp2(gr%nz) = zero


    return

  end subroutine wp2_term_pr1_rhs

  !=============================================================================
  pure subroutine wp2_term_pr_dfsn_rhs( gr, C_wp2_pr_dfsn, &
                                        rho_ds_zt, invrs_rho_ds_zm, &
                                        wpup2, wpvp2, wp3, &
                                        rhs_pr_dfsn_wp2 )

    ! Description:
    !
    ! This term is intended to represent the "diffusion" part of the wp2 
    ! pressure correlation.  The total pressure diffusion term, 
    ! 
    !   -1 / rho * ( d( <u_k'p'> )/dx_i + d( <u_i'p'> )/dx_k  )
    !
    ! becomes 
    !
    !   -2 / rho * d( <w'p'> )/dz
    !
    ! for the w'^2 equation.  The factor of 2 is replaced with a tunable
    ! parameter, C_wp2_pr_dfsn, and p' is replaced with 
    !
    !   p' ~ - rho * ( u_i*u_i - <u_i*u_i> ),
    !
    ! following Lumley 1978.  The wp2 pressure diffusion term becomes
    !
    !   + C_wp2_pr_dfsn / rho * ( d( rho*<w'u_iu_i> )/dz )
    !
    ! References:
    !   Lumley 1978, p. 170.  See eq. 6.47 and accompanying discussion.
    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid  ! Type

    use constants_clubb, only: &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    type (grid), target, intent(in) :: gr

    real( kind = core_rknd ), intent(in) :: &
      C_wp2_pr_dfsn      ! Model parameter C_wp2_pr_dfsn                [-]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      invrs_rho_ds_zm, & ! Inverse dry static density (thermo levels) [kg/m^3] 
      rho_ds_zt,       & ! Dry static density on mom. levels       [kg/m^3]
      wpup2,           & ! w'u'^2 on thermodynamic levels          [m^4/s^4]
      wpvp2,           & ! w'v'^2 on thermodynamic levels          [m^4/s^4]
      wp3                ! w'^3 on thermo levels                   [m^4/s^4]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      rhs_pr_dfsn_wp2    ! RHS portion of wp2 from pressure-diffusion correlation [m^3/s^4]

    ! Local Variables
    integer :: k   ! Vertical level index 

    real( kind = core_rknd ), dimension(gr%nz) :: &
      wpuip2            ! 4th-order moment sum <w'u_i'u_i'>     [m^4/s^4]

    ! ---- Begin Code ----

    wpuip2 = wpup2 + wpvp2 + wp3

    do k = 2, gr%nz-1
      rhs_pr_dfsn_wp2(k) &
       = + C_wp2_pr_dfsn * invrs_rho_ds_zm(k) * gr%invrs_dzm(k) &
         * ( rho_ds_zt(k+1) * wpuip2(k+1) - rho_ds_zt(k) * wpuip2(k) )
    enddo ! k = 2, gr%nz-1

    ! Set lower boundary condition
    rhs_pr_dfsn_wp2(1) = rhs_pr_dfsn_wp2(2)

    ! Set upper boundary to 0
    rhs_pr_dfsn_wp2(gr%nz) = zero


    return

  end subroutine wp2_term_pr_dfsn_rhs

  !=============================================================================
  pure subroutine wp3_term_ta_new_pdf_lhs( gr, coef_wp4_implicit, wp2, rho_ds_zm, &
                                           invrs_rho_ds_zt, invrs_dzt, &
                                           lhs_ta_wp3 )

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
    ! =======coef_wp4_implicit(interp)=======wp2=========rho_ds_zm======= m(k-1)
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

    use grid_class, only:  &
        grid ! Type

    use constants_clubb, only: &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Constant parameters
    integer, parameter :: & 
      k_mdiag   = 1, & ! Momentum superdiagonal index.
      km1_mdiag = 2    ! Momentum subdiagonal index. 

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      coef_wp4_implicit, & ! <w'^4>=coef_wp4_implicit*<w'^2>^2; m-levs [-]
      wp2,               & ! <w'^2>                                    [m^2/s^2]
      rho_ds_zm,         & ! Dry, static density at momentum levels    [kg/m^3]
      invrs_rho_ds_zt,   & ! Inv dry, static density at thermo levels  [m^3/kg]
      invrs_dzt            ! Inverse of grid spacing                   [1/m]

    ! Return Variable
    real( kind = core_rknd ), dimension(2,gr%nz), intent(out) :: &
      lhs_ta_wp3   ! LHS coefficient of wp3 turbulent advection  [m/s^2]

    ! Local Variable
    integer :: k    ! Vertical index


    ! Set term at lower boundary to 0
    lhs_ta_wp3(k_mdiag,1) = zero
    lhs_ta_wp3(km1_mdiag,1) = zero

    ! Calculate term at all interior grid levels.
    do k = 2, gr%nz-1

       ! Momentum superdiagonal: [ x wp2(k,<t+1>) ]
       lhs_ta_wp3(k_mdiag,k) &
       = + invrs_rho_ds_zt(k) * invrs_dzt(k) * rho_ds_zm(k) &
           * coef_wp4_implicit(k) * wp2(k)

       ! Momentum subdiagonal: [ x wp2(k-1,<t+1>) ]
       lhs_ta_wp3(km1_mdiag,k) &
       = - invrs_rho_ds_zt(k) * invrs_dzt(k) * rho_ds_zm(k-1) &
           * coef_wp4_implicit(k-1) * wp2(k-1)

    enddo ! k = 2, gr%nz-1

    ! Set term at upper boundary to 0
    lhs_ta_wp3(k_mdiag,gr%nz) = zero
    lhs_ta_wp3(km1_mdiag,gr%nz) = zero


    return

  end subroutine wp3_term_ta_new_pdf_lhs

  !=============================================================================
  pure subroutine wp3_term_ta_ADG1_lhs( gr, wp2, a1, a1_zt, a3, a3_zt, &
                                        wp3_on_wp2, rho_ds_zm, &
                                        rho_ds_zt, invrs_rho_ds_zt, &
                                        invrs_dzt, l_standard_term_ta, &
                                        l_partial_upwind_wp3, &
                                        lhs_ta_wp3 )

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
    ! ------------------------------------------------wp3---------------- t(k+1)
    !
    ! ===a3====wp2====rho_ds_zm====a1======================wp3(interp)=== m(k)
    !
    ! -----------dF/dz----invrs_rho_ds_zt----dG/dz----wp3---------------- t(k)
    !
    ! ===a3====wp2====rho_ds_zm====a1======================wp3(interp)=== m(k-1)
    !
    ! ------------------------------------------------wp3---------------- t(k-1)
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
        grid ! Type

    use constants_clubb, only: &
        zero

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

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
    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      wp2,             & ! w'^2                                     [m^2/s^2]
      a1,              & ! a_1                                      [-]
      a1_zt,           & ! a_1 interpolated to thermodynamic levels [-]
      a3,              & ! a_3                                      [-]
      a3_zt,           & ! a_3 interpolated to thermodynamic levels [-]
      wp3_on_wp2,      & ! w'^3 / w'^2 at momentum levels           [m/s]
      rho_ds_zm,       & ! Dry, static density at momentum levels   [kg/m^3]
      rho_ds_zt,       & ! Dry, static density at thermo. levels    [kg/m^3]
      invrs_rho_ds_zt, & ! Inv dry, static density at thermo levels [m^3/kg]
      invrs_dzt          ! Inverse of grid spacing                  [1/m]

    logical, intent(in) :: &
      l_standard_term_ta,   & ! Use the standard discretization for the
                              ! turbulent advection terms.  Setting to .false.
                              ! means that a_1 and a_3 are pulled outside of the
                              ! derivative in advance_wp2_wp3_module.F90 and in
                              ! advance_xp2_xpyp_module.F90.
      l_partial_upwind_wp3    ! Flag to use an "upwind" discretization rather
                              ! than a centered discretization for the portion
                              ! of the wp3 turbulent advection term for ADG1
                              ! that is linearized in terms of wp3<t+1>.
                              ! (Requires ADG1 PDF and l_standard_term_ta).

    ! Output Variable
    real( kind = core_rknd ), dimension(5,gr%nz), intent(out) :: &
      lhs_ta_wp3    ! LHS coefficient of wp3 turbulent advection

    ! Local variables
    integer :: k    ! Vertical level index


    ! Set lower boundary to 0
    lhs_ta_wp3(:,1) = zero

    ! Calculate term at all interior grid levels.
    if ( l_standard_term_ta ) then

       ! The turbulent advection term is discretized normally, in accordance
       ! with the model equations found in the documentation and the description
       ! listed above.

       if ( .not. l_partial_upwind_wp3 ) then

          ! All portions of the wp3 turbulent advection term for ADG1 use
          ! centered discretization in accordance with description and diagram
          ! shown above.

          do k = 2, gr%nz-1, 1

             ! Thermodynamic superdiagonal: [ x wp3(k+1,<t+1>) ]
             lhs_ta_wp3(kp1_tdiag,k) &
             = + invrs_rho_ds_zt(k) &
                 * invrs_dzt(k) &
                   * rho_ds_zm(k) * a1(k) * wp3_on_wp2(k) &
                   * gr%weights_zt2zm(t_above,k)

             ! Momentum superdiagonal: [ x wp2(k,<t+1>) ]
             lhs_ta_wp3(k_mdiag,k) &
             = + invrs_rho_ds_zt(k) * invrs_dzt(k) &
                 * rho_ds_zm(k) * a3(k) * wp2(k)

             ! Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
             lhs_ta_wp3(k_tdiag,k) &
             = + invrs_rho_ds_zt(k) &
                 * invrs_dzt(k) &
                   * ( rho_ds_zm(k) * a1(k) * wp3_on_wp2(k) &
                       * gr%weights_zt2zm(t_below,k) &
                       - rho_ds_zm(k-1) * a1(k-1) * wp3_on_wp2(k-1) &
                         * gr%weights_zt2zm(t_above,k-1) )

             ! Momentum subdiagonal: [ x wp2(k-1,<t+1>) ]
             lhs_ta_wp3(km1_mdiag,k) &
             = - invrs_rho_ds_zt(k) * invrs_dzt(k) &
                 * rho_ds_zm(k-1) * a3(k-1) * wp2(k-1)

             ! Thermodynamic subdiagonal: [ x wp3(k-1,<t+1>) ]
             lhs_ta_wp3(km1_tdiag,k) &
             = - invrs_rho_ds_zt(k) &
                 * invrs_dzt(k) &
                   * rho_ds_zm(k-1) * a1(k-1) * wp3_on_wp2(k-1) &
                   * gr%weights_zt2zm(t_below,k-1)

          enddo ! k = 2, gr%nz-1, 1

       else ! l_partial_upwind_wp3

          ! Partial upwinding of the wp3 turbulent advection term, where the
          ! portion of the wp3 turbulent advection term that is linearized in
          ! terms of wp2<t+1> is still handled using centered discretization,
          ! but the portion of the term that is linearized in terms of wp3<t+1>
          ! is handled using an "upwind" discretization that also takes into
          ! "winds" that converge or diverge around the central thermodynamic
          ! grid level.  Provided by Chris Vogl and Shixuan Zhang.

          do k = 2, gr%nz-1, 1

             ! Thermodynamic superdiagonal: [ x wp3(k+1,<t+1>) ]
             lhs_ta_wp3(kp1_tdiag,k) &
             = + invrs_rho_ds_zt(k) &
                 * invrs_dzt(k) * rho_ds_zt(k+1) &
                 * min( a1(k) * wp3_on_wp2(k), zero )

             ! Momentum superdiagonal: [ x wp2(k,<t+1>) ]
             lhs_ta_wp3(k_mdiag,k) &
             = + invrs_rho_ds_zt(k) * invrs_dzt(k) &
                 * rho_ds_zm(k) * a3(k) * wp2(k)

             ! Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
             lhs_ta_wp3(k_tdiag,k) &
             = + invrs_rho_ds_zt(k) &
                 * invrs_dzt(k) * rho_ds_zt(k) &
                 * ( max( a1(k) * wp3_on_wp2(k), zero ) &
                     - min( a1(k-1) * wp3_on_wp2(k-1), zero ) )

             ! Momentum subdiagonal: [ x wp2(k-1,<t+1>) ]
             lhs_ta_wp3(km1_mdiag,k) &
             = - invrs_rho_ds_zt(k) * invrs_dzt(k) &
                 * rho_ds_zm(k-1) * a3(k-1) * wp2(k-1)

             ! Thermodynamic subdiagonal: [ x wp3(k-1,<t+1>) ]
             lhs_ta_wp3(km1_tdiag,k) &
             = - invrs_rho_ds_zt(k) &
                 * invrs_dzt(k) * rho_ds_zt(k-1) &
                 * max( a1(k-1) * wp3_on_wp2(k-1), zero )

          enddo ! k = 2, gr%nz-1, 1

       endif ! .not. l_partial_upwind_wp3

    else

       ! Alternate discretization for the turbulent advection term, which
       ! contains the term:
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

       do k = 2, gr%nz-1

          ! Thermodynamic superdiagonal: [ x wp3(k+1,<t+1>) ]
          lhs_ta_wp3(kp1_tdiag,k) &
          = + invrs_rho_ds_zt(k) &
              * a1_zt(k) * invrs_dzt(k) &
              * rho_ds_zm(k) * wp3_on_wp2(k) &
              * gr%weights_zt2zm(t_above,k)

          ! Momentum superdiagonal: [ x wp2(k,<t+1>) ]
          lhs_ta_wp3(k_mdiag,k) &
          = + invrs_rho_ds_zt(k) * a3_zt(k) * invrs_dzt(k) &
              * rho_ds_zm(k) * wp2(k)

          ! Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
          lhs_ta_wp3(k_tdiag,k) &
          = + invrs_rho_ds_zt(k) &
              * a1_zt(k) * invrs_dzt(k) &
                * ( rho_ds_zm(k) * wp3_on_wp2(k) &
                    * gr%weights_zt2zm(t_below,k) &
                    - rho_ds_zm(k-1) * wp3_on_wp2(k-1) &
                      * gr%weights_zt2zm(t_above,k-1) )

          ! Momentum subdiagonal: [ x wp2(k-1,<t+1>) ]
          lhs_ta_wp3(km1_mdiag,k) &
          = - invrs_rho_ds_zt(k) * a3_zt(k) * invrs_dzt(k) &
              * rho_ds_zm(k-1) * wp2(k-1)

          ! Thermodynamic subdiagonal: [ x wp3(k-1,<t+1>) ]
          lhs_ta_wp3(km1_tdiag,k) &
          = - invrs_rho_ds_zt(k) &
              * a1_zt(k) * invrs_dzt(k) &
              * rho_ds_zm(k-1) * wp3_on_wp2(k-1) &
              * gr%weights_zt2zm(t_below,k-1)

       enddo ! k = 2, gr%nz-1

    endif ! l_standard_term_ta

    ! Set upper boundary to 0
    lhs_ta_wp3(:,gr%nz) = zero


    return

  end subroutine wp3_term_ta_ADG1_lhs

  !=============================================================================
  pure subroutine wp3_term_tp_lhs( gr, coef_wp3_tp, wp2, &
                                   rho_ds_zm, &
                                   invrs_rho_ds_zt, &
                                   invrs_dzt, &
                                   lhs_tp_wp3 )

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
    ! ====wp2=========rho_ds_zm========================================== m(k-1)
    !
    ! The vertical indices m(k), t(k), and m(k-1) correspond with altitudes
    ! zm(k), zt(k), and zm(k-1), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) )

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        grid ! Type

    use constants_clubb, only:  &
        three,        & ! Constant(s)
        three_halves, &
        zero

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Constant parameters
    integer, parameter :: & 
      k_mdiag   = 1, & ! Momentum superdiagonal index.
      km1_mdiag = 2    ! Momentum subdiagonal index. 

    ! Input Variables
   real( kind = core_rknd ), intent(in) :: &
      coef_wp3_tp      ! Coefficient for tp pressure scrambling term   [-]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      wp2,             & ! w'^2                                        [m^2/s^2]
      rho_ds_zm,       & ! Dry, static density at momentum levels      [kg/m^3]
      invrs_rho_ds_zt, & ! Inv dry, static density at thermo levels    [m^3/kg]
      invrs_dzt          ! Inverse of grid spacing                     [1/m]

    ! Return Variable
    real( kind = core_rknd ), dimension(2,gr%nz), intent(out) :: &
      lhs_tp_wp3    ! LHS coefficient of wp3 turbulent production  [1/s]

    ! Local variables
    integer :: k    ! Vertical level index


    ! Set lower boundary to 0
    lhs_tp_wp3(k_mdiag,1) = zero
    lhs_tp_wp3(km1_mdiag,1) = zero

    ! Calculate term at all interior grid levels.
    do k = 2, gr%nz-1

       ! Momentum superdiagonal: [ x wp2(k,<t+1>) ]
       lhs_tp_wp3(k_mdiag,k) &
       = - coef_wp3_tp * three * invrs_rho_ds_zt(k) * invrs_dzt(k) &
           * rho_ds_zm(k) * wp2(k) &
         + coef_wp3_tp * three_halves * invrs_dzt(k) * wp2(k)

       ! Momentum subdiagonal: [ x wp2(k-1,<t+1>) ]
       lhs_tp_wp3(km1_mdiag,k) &
       = + coef_wp3_tp * three * invrs_rho_ds_zt(k) * invrs_dzt(k) &
           * rho_ds_zm(k-1) * wp2(k-1) &
         - coef_wp3_tp * three_halves * invrs_dzt(k) * wp2(k-1)

    enddo ! k = 2, gr%nz-1

    ! Set upper boundary to 0
    lhs_tp_wp3(k_mdiag,gr%nz) = zero
    lhs_tp_wp3(km1_mdiag,gr%nz) = zero


    return

  end subroutine wp3_term_tp_lhs

  !=============================================================================
  pure subroutine wp3_terms_ac_pr2_lhs( gr, C11_Skw_fnc, wm_zm, invrs_dzt, &
                                        lhs_ac_pr2_wp3 )

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
    ! =======wm_zm============================================= m(k-1)
    !
    ! The vertical indices m(k), t(k), and m(k-1) correspond with altitudes 
    ! zm(k), zt(k), and zm(k-1), respectively.  The letter "t" is used for 
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) )

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid ! Type

    use constants_clubb, only: &
        three, & ! Variable(s)
        one,   &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      C11_Skw_fnc,  & ! C_11 parameter with Sk_w applied       [-]
      wm_zm,        & ! w wind component at momentum levels    [m/s]
      invrs_dzt       ! Inverse of grid spacing                [1/m]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: & 
      lhs_ac_pr2_wp3     ! LHS coefficient of wp3 from terms ac and pr2 [1/s]

    ! Local variable
    integer :: k    ! Vertical level index


    ! Set lower boundary to 0
    lhs_ac_pr2_wp3(1) = zero

    ! Calculate term at all interior grid levels.
    do k = 2, gr%nz-1

       ! Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
       lhs_ac_pr2_wp3(k) &
       = + ( one - C11_Skw_fnc(k) ) &
           * three * invrs_dzt(k) * ( wm_zm(k) - wm_zm(k-1) )

    enddo ! k = 2, gr%nz-1

    ! Set upper boundary to 0
    lhs_ac_pr2_wp3(gr%nz) = zero


    return

  end subroutine wp3_terms_ac_pr2_lhs

  !=============================================================================
  pure subroutine wp3_term_pr1_lhs( gr, C8, C8b, invrs_tau_wp3_zt, Skw_zt, &
                                    l_damp_wp3_Skw_squared, &
                                    lhs_pr1_wp3 )

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

    use grid_class, only: &
        grid ! Type

    use constants_clubb, only: &
        one, & ! Variable(s)
        three, &
        five, &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      invrs_tau_wp3_zt,  & ! Inverse time-scale tau at thermodynamic levels  [1/s]
      Skw_zt               ! Skewness of w at thermodynamic levels   [-]

    real( kind = core_rknd ), intent(in) :: & 
      C8,      & ! Model parameter C_8                     [-]
      C8b        ! Model parameter C_8b                    [-]

    logical, intent(in) :: &
      l_damp_wp3_Skw_squared ! Set damping on wp3 to use Skw^2 rather than Skw^4

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      lhs_pr1_wp3    ! LHS coefficient of wp3 from pressure term 1  [1/s]

    ! Local variables
    integer :: k    ! Vertical level index


    ! Set lower boundary to 0
    lhs_pr1_wp3(1) = zero


    ! Calculate term at all interior grid levels.
    if ( l_damp_wp3_Skw_squared ) then
 
       do k = 2, gr%nz-1

          ! Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
          lhs_pr1_wp3(k) &
          = + ( C8 * invrs_tau_wp3_zt(k) ) * ( three * C8b * Skw_zt(k)**2 + one )

       enddo ! k = 2, gr%nz-1

    else

       do k = 2, gr%nz-1

          ! Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
          lhs_pr1_wp3(k) &
          = + ( C8 * invrs_tau_wp3_zt(k) ) * ( five * C8b * Skw_zt(k)**4 + one )

       enddo ! k = 2, gr%nz-1

    endif ! l_damp_wp3_Skw_squared

    ! Set upper boundary to 0
    lhs_pr1_wp3(gr%nz) = zero


    return

  end subroutine wp3_term_pr1_lhs

  !=============================================================================
  pure subroutine wp3_term_ta_explicit_rhs( gr, wp4, &
                                            rho_ds_zm, &
                                            invrs_rho_ds_zt, &
                                            invrs_dzt, &
                                            rhs_ta_wp3 )

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
    ! =========wp4(interp)===========rho_ds_zm=========================== m(k-1)
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

    use grid_class, only: &
        grid ! Type

    use constants_clubb, only: &
        zero    ! Constant(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      wp4,             & ! <w'^4>                                   [m^4/s^4]
      rho_ds_zm,       & ! Dry, static density at momentum level    [kg/m^3]
      invrs_rho_ds_zt, & ! Inv dry, static density at thermo level  [m^3/kg]
      invrs_dzt          ! Inverse of grid spacing                  [1/m]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      rhs_ta_wp3    ! Rate of change of wp3 from turbulent advection  [m^3/s^4]

    ! Local variables
    integer :: k    ! Vertical level index


    ! Set lower boundary to 0
    rhs_ta_wp3(1) = zero

    ! Calculate term at all interior grid levels.
    do k = 2, gr%nz

       rhs_ta_wp3(k) &
       = - invrs_rho_ds_zt(k) * invrs_dzt(k) &
           * ( rho_ds_zm(k) * wp4(k) - rho_ds_zm(k-1) * wp4(k-1) )

    enddo ! k = 2, gr%nz

    ! Set upper boundary to 0
    rhs_ta_wp3(gr%nz) = zero


    return

  end subroutine wp3_term_ta_explicit_rhs

  !=============================================================================
  pure subroutine wp3_terms_bp1_pr2_rhs( gr, C11_Skw_fnc, thv_ds_zt, wp2thvp, &
                                         rhs_bp1_pr2_wp3 )

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

    use grid_class, only: &
        grid ! Type

    use constants_clubb, only: & ! Constant(s) 
        grav,  & ! Gravitational acceleration [m/s^2]
        three, &
        one,   &
        zero

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      C11_Skw_fnc, & ! C_11 parameter with Sk_w applied         [-]
      thv_ds_zt,   & ! Dry, base-state theta_v at thermo. levs  [K]
      wp2thvp        ! w'^2 th_v'                               [K m^2/s^2]

    ! Return Variable
    real( kind = core_rknd ),  dimension(gr%nz), intent(out) :: &
      rhs_bp1_pr2_wp3   ! RHS portion of wp3 from terms bp1 and pr2 [m^3/s^4]

    ! Local Variables
    integer :: k    ! Vertical loop index


    ! Set lower boundary to 0
    rhs_bp1_pr2_wp3(1) = zero

    ! Calculate term at all interior grid levels.
    do k = 2, gr%nz-1

       rhs_bp1_pr2_wp3(k) &
       = + ( one - C11_Skw_fnc(k) ) &
           * three * ( grav / thv_ds_zt(k) ) * wp2thvp(k)

    enddo ! k = 2, gr%nz-1

    ! Set upper boundary to 0
    rhs_bp1_pr2_wp3(gr%nz) = zero


    return

  end subroutine wp3_terms_bp1_pr2_rhs

  !=============================================================================
  pure subroutine wp3_term_pr_turb_rhs( gr, C_wp3_pr_turb, Kh_zt, wpthvp, &
                                        dum_dz, dvm_dz, &
                                        upwp, vpwp, &
                                        thv_ds_zt, invrs_dzt, &
                                        rho_ds_zm, invrs_rho_ds_zt,  &
                                        em, wp2, &
                                        rhs_pr_turb_wp3, &
                                        l_use_tke_in_wp3_pr_turb_term )
    use grid_class, only: grid

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

    use grid_class, only: &
        grid, &
        zm2zt    ! Variable type(s)

    use constants_clubb, only: & ! Constant(s) 
        grav, & ! Gravitational acceleration [m/s^2]
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      C_wp3_pr_turb         ! Model parameter C_wp3_pr_turb                [-]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      Kh_zt,           & ! Eddy-diffusivity on moment. levels      [m^2/s]
      wpthvp,          & ! w'th_v'                                 [K m/s]
      dum_dz,          & ! derivative of u wind with respect to z  [m/s]
      dvm_dz,          & ! derivative of v wind with respect to z  [m/s]
      upwp,            & ! u'v'                                    [m^2/s^2]
      vpwp,            & ! v'w'                                    [m^2/s^2]
      thv_ds_zt,       & ! Dry, base-state theta_v at thermo. levs [K]
      invrs_dzt,       & ! Inverse of grid spacing                 [1/m]
      rho_ds_zm,       & ! Dry static density on mom. levels       [kg/m^3]
      invrs_rho_ds_zt, & ! Inverse dry static density on thermo. levs [kg/m^3]
      wp2,             & ! w'^2                                    [m^2/s^2]
      em                 ! Turbulence kinetic energy               [m^2/s^2]

    logical, intent(in) :: &
      l_use_tke_in_wp3_pr_turb_term  ! Use TKE formulation for wp3 pr_turb term

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      rhs_pr_turb_wp3    ! RHS portion of wp3 from pressure-turbulence correlation [m^3/s^4]

    ! Local Variables
    integer :: k   ! Vertical level index 

    ! ---- Begin Code ----

    ! Set lower boundary to 0
    rhs_pr_turb_wp3(1) = zero

    do k = 2, gr%nz-1

      if ( .not. l_use_tke_in_wp3_pr_turb_term ) then

        rhs_pr_turb_wp3(k) &
        = - C_wp3_pr_turb * Kh_zt(k) * invrs_dzt(k) &
            * ( grav / thv_ds_zt(k) * ( wpthvp(k) - wpthvp(k-1) ) &
                - ( upwp(k) * dum_dz(k) - upwp(k-1) * dum_dz(k-1) ) &
                - ( vpwp(k) * dvm_dz(k) - vpwp(k-1) * dvm_dz(k-1) ) )

      else

        rhs_pr_turb_wp3(k) &
        = - C_wp3_pr_turb * invrs_rho_ds_zt(k) * invrs_dzt(k) &
            * ( rho_ds_zm(k) * wp2(k) * em(k) - rho_ds_zm(k-1) * wp2(k-1) * em(k-1) )

      endif

    enddo ! k = 2, gr%nz-1

    ! Set upper boundary to 0
    rhs_pr_turb_wp3(gr%nz) = zero


    return

  end subroutine wp3_term_pr_turb_rhs

  !=============================================================================
  pure subroutine wp3_term_pr_dfsn_rhs( gr, C_wp3_pr_dfsn, &
                                        rho_ds_zm, invrs_rho_ds_zt, &
                                        wp2up2, wp2vp2, wp4, &
                                        up2, vp2, wp2, &
                                        rhs_pr_dfsn_wp3 )

    ! Description:
    !
    ! This term is intended to represent the "diffusion" part of the total wp3 
    ! pressure correlation.  The total wp3 pressure term, -3w'^2/rho*dp'/dz, can be
    ! split into
    ! 
    !   -3w'^2/rho*dp'/dz = + 3p'/rho*d(w'^2)/dz - 3/rho*d(w'^2p')/dz 
    !
    ! using the product rule.  The second term on the RHS we consider to be the
    ! diffusion part, calculated by this subroutine.  We replace the factor of 3
    ! with a tunable parameter, C_wp3_pr_dfsn, and we replace p' with
    !
    !   p' ~ - rho * ( u_i*u_i - <u_i*u_i> ),
    !
    ! following Lumley 1978.  The wp3 pressure diffusion term then becomes
    !
    !   + C_wp3_pr_dfsn / rho * ( d( rho*( <w'^2u_i'u_i'> - <w'^2>*<u_i'u_i'> ) )/dz )
    !
    ! References:
    !   Lumley 1978, p. 170.  See eq. 6.47 and accompanying discussion.
    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid ! Type

    use constants_clubb, only: &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      C_wp3_pr_dfsn      ! Model parameter C_wp3_pr_dfsn              [-]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      invrs_rho_ds_zt, & ! Inverse dry static density (thermo levels) [kg/m^3] 
      rho_ds_zm,       & ! Dry static density on mom. levels          [kg/m^3]
      wp2up2,          & ! w'^2u'^2 on momentum levels                [m^4/s^4]
      wp2vp2,          & ! w'^2v'^2 on momentum levels                [m^4/s^4]
      wp4,             & ! w'^4 on momentum levels                    [m^4/s^4]
      up2,             & ! u'^2 on momentum levels                    [m^2/s^2]
      vp2,             & ! v'^2 on momentum levels                    [m^2/s^2]
      wp2                ! w'^2 on momentum levels                    [m^2/s^2]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      rhs_pr_dfsn_wp3    ! RHS portion of wp3 from pressure-diffusion correlation [m^3/s^4]

    ! Local Variables
    integer :: k   ! Vertical level index 

    real( kind = core_rknd ), dimension(gr%nz) :: &
      wp2uip2,   & ! 4th-order moment sum <w'^2u_i'u_i'>     [m^4/s^4]
      wp2_uip2     ! 2nd-order moment sum <w'^2>*<u_i'u_i'>  [m^4/s^4]

    ! ---- Begin Code ----

    wp2uip2 = wp2up2 + wp2vp2 + wp4
    wp2_uip2 = wp2*up2 + wp2*vp2 + wp2*wp2

    do k = 2, gr%nz-1
      rhs_pr_dfsn_wp3(k) &
       = + C_wp3_pr_dfsn * invrs_rho_ds_zt(k) * gr%invrs_dzt(k) &
         * ( rho_ds_zm(k) * ( wp2uip2(k) - wp2_uip2(k) ) &
           - rho_ds_zm(k-1) * ( wp2uip2(k-1) - wp2_uip2(k-1) ) )
    enddo ! k = 2, gr%nz-1

    ! Set lower boundary condition
    rhs_pr_dfsn_wp3(1) = zero

    ! Set upper boundary to 0
    rhs_pr_dfsn_wp3(gr%nz) = zero

    return

  end subroutine wp3_term_pr_dfsn_rhs

  !=============================================================================
  pure subroutine wp3_term_pr1_rhs( gr, C8, C8b, invrs_tau_wp3_zt, Skw_zt, wp3, &
                                    l_damp_wp3_Skw_squared, &
                                    rhs_pr1_wp3 )

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

    use grid_class, only: &
        grid ! Type

    use constants_clubb, only: &
        two,  & ! Constant(s)
        four, &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C8,  & ! Model parameter C_8                        [-]
      C8b    ! Model parameter C_8b                       [-]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      invrs_tau_wp3_zt, & ! Inverse time-scale tau at thermodynamic levels  [1/s]
      Skw_zt,           & ! Skewness of w at thermodynamic levels      [-]
      wp3                 ! w'^3                                       [m^3/s^3]

    logical, intent(in) :: &
      l_damp_wp3_Skw_squared ! Set damping on wp3 to use Skw^2 rather than Skw^4

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      rhs_pr1_wp3    ! RHS portion of wp3 from pressure term 1  [m^3/s^4]

    ! Local Variables
    integer :: k    ! Vertical level index


    ! Set lower boundary to 0
    rhs_pr1_wp3(1) = zero

    ! Calculate term at all interior grid levels.
    if ( l_damp_wp3_Skw_squared ) then

       do k = 2, gr%nz-1

          rhs_pr1_wp3(k) &
          = + ( C8 * invrs_tau_wp3_zt(k) ) * ( two * C8b * Skw_zt(k)**2 ) * wp3(k)

       enddo ! k = 2, gr%nz-1

    else

       do k = 2, gr%nz-1

          rhs_pr1_wp3(k) &
          = + ( C8 * invrs_tau_wp3_zt(k) ) * ( four * C8b * Skw_zt(k)**4 ) * wp3(k)

       enddo ! k = 2, gr%nz-1

    endif ! l_damp_wp3_Skw_squared

    ! Set upper boundary to 0
    rhs_pr1_wp3(gr%nz) = zero


    return

  end subroutine wp3_term_pr1_rhs

!===============================================================================

end module advance_wp2_wp3_module
