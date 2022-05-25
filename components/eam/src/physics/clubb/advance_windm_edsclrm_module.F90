!------------------------------------------------------------------------
! $Id$
!===============================================================================
module advance_windm_edsclrm_module

  implicit none

  private ! Set Default Scope

  public :: advance_windm_edsclrm

  private :: windm_edsclrm_solve, &
    compute_uv_tndcy,  &
    windm_edsclrm_lhs, &
    windm_edsclrm_rhs
    

  ! Private named constants to avoid string comparisons
  integer, parameter, private :: &
    windm_edsclrm_um = 1, &     ! Named constant to handle um solves
    windm_edsclrm_vm = 2, &     ! Named constant to handle vm solves
    windm_edsclrm_scalar = 3, & ! Named constant to handle scalar solves
    clip_upwp = 10, &           ! Named constant for upwp clipping
                                ! NOTE: This must be the same as the clip_upwp
                                ! declared in clip_explicit!
    clip_vpwp = 11              ! Named constant for vpwp clipping
                                ! NOTE: This must be the same as the clip_vpwp
                                ! declared in clip_explicit!

  contains

  !=============================================================================
  subroutine advance_windm_edsclrm( nz, ngrdcol, gr, dt, &
                                    wm_zt, Km_zm, Kmh_zm, &
                                    ug, vg, um_ref, vm_ref, &
                                    wp2, up2, vp2, um_forcing, vm_forcing, &
                                    edsclrm_forcing, &
                                    rho_ds_zm, invrs_rho_ds_zt, &
                                    fcor, l_implemented, &
                                    nu_vert_res_dep, &
                                    l_predict_upwp_vpwp, &
                                    l_upwind_xm_ma, &
                                    l_uv_nudge, &
                                    l_tke_aniso, &
                                    l_lmm_stepping, &
                                    l_linearize_pbl_winds, &
                                    order_xp2_xpyp, order_wp2_wp3, order_windm, &
                                    stats_zt, stats_zm, stats_sfc, & 
                                    um, vm, edsclrm, &
                                    upwp, vpwp, wpedsclrp, &
                                    um_pert, vm_pert, upwp_pert, vpwp_pert )

    ! Description:
    ! Solves for both mean horizontal wind components, um and vm, and for the
    ! eddy-scalars (passive scalars that don't use the high-order closure).

    ! Uses the LAPACK tridiagonal solver subroutine with 2 + # of scalar(s)
    ! back substitutions (since the left hand side matrix is the same for all
    ! input variables).

    ! References:
    ! Eqn. 8 & 9 on p. 3545 of
    ! ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    ! Method and Model Description'' Golaz, et al. (2002)
    ! JAS, Vol. 59, pp. 3540--3551.
    !-----------------------------------------------------------------------

    use grid_class, only:  &
        grid, &  ! Type
        zm2zt

    use parameters_model, only:  &
        ts_nudge,  & ! Variable(s)
        edsclr_dim

    use parameters_tunable, only: &
        nu_vertical_res_dep    ! Type(s)

    use clubb_precision, only:  &
        core_rknd ! Variable(s)

    use stats_type_utilities, only: &
        stat_begin_update, & ! Subroutines
        stat_end_update, &
        stat_update_var

    use stats_variables, only: &
        ium_ref, & ! Variables
        ivm_ref, &
        ium_sdmp, &
        ivm_sdmp, &
        ium_ndg, &
        ivm_ndg, &
        iwindm_matrix_condt_num, &
        l_stats_samp

    use clip_explicit, only:  &
        clip_covar  ! Procedure(s)

    use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_fatal_error              ! Constant

    use constants_clubb, only:  &
        one_half, & ! Constant(s)
        zero,     &
        fstderr,  &
        eps

    use sponge_layer_damping, only: &
        uv_sponge_damp_settings, &
        uv_sponge_damp_profile, &
        sponge_damp_xm     ! Procedure(s)

    use stats_type, only: stats ! Type

    use diffusion, only:  & 
        diffusion_zt_lhs   ! Procedure(s)

    use mean_adv, only: & 
        term_ma_zt_lhs    ! Procedures
        
    use model_flags, only: &
        l_upwind_Kh_dp_term
        
    use advance_helper_module, only: &
        calc_xpwp

    implicit none

    ! ------------------------ Input Variables ------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, dimension(ngrdcol), intent(in) :: &
      gr
  
    real( kind = core_rknd ), intent(in) ::  &
      dt                 ! Model timestep                             [s]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
      wm_zt,           & ! w wind component on thermodynamic levels    [m/s]
      Km_zm,           & ! Eddy diffusivity of winds on momentum levs. [m^2/s]
      Kmh_zm,          & ! Eddy diffusivity of themo on momentum levs. [m^s/s]
      ug,              & ! u (west-to-east) geostrophic wind comp.     [m/s]
      vg,              & ! v (south-to-north) geostrophic wind comp.   [m/s]
      um_ref,          & ! Reference u wind component for nudging      [m/s]
      vm_ref,          & ! Reference v wind component for nudging      [m/s]
      wp2,             & ! w'^2 (momentum levels)                      [m^2/s^2]
      up2,             & ! u'^2 (momentum levels)                      [m^2/s^2]
      vp2,             & ! v'^2 (momentum levels)                      [m^2/s^2]
      um_forcing,      & ! u forcing                                   [m/s/s]
      vm_forcing,      & ! v forcing                                   [m/s/s]
      rho_ds_zm,       & ! Dry, static density on momentum levels      [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density at thermo. levels  [m^3/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nz,edsclr_dim), intent(in) ::  &
      edsclrm_forcing  ! Eddy scalar large-scale forcing        [{units vary}/s]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
      fcor           ! Coriolis parameter                            [s^-1]

    logical, intent(in) ::  &
      l_implemented  ! Flag for CLUBB being implemented in a larger model.

    type(nu_vertical_res_dep), dimension(ngrdcol), intent(in) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    logical, intent(in) :: &
      l_predict_upwp_vpwp,   & ! Flag to predict <u'w'> and <v'w'> along with <u> and <v> alongside
                               ! the advancement of <rt>, <w'rt'>, <thl>, <wpthlp>, <sclr>, and
                               ! <w'sclr'> in subroutine advance_xm_wpxp.  Otherwise, <u'w'> and
                               ! <v'w'> are still approximated by eddy diffusivity when <u> and <v>
                               ! are advanced in subroutine advance_windm_edsclrm.
      l_upwind_xm_ma,        & ! This flag determines whether we want to use an upwind differencing
                               ! approximation rather than a centered differencing for turbulent or
                               ! mean advection terms. It affects rtm, thlm, sclrm, um and vm.
      l_uv_nudge,            & ! For wind speed nudging
      l_tke_aniso,           & ! For anisotropic turbulent kinetic energy, i.e. TKE = 1/2
                               ! (u'^2 + v'^2 + w'^2)
      l_lmm_stepping,        & ! Apply Linear Multistep Method (LMM) Stepping
      l_linearize_pbl_winds    ! Flag (used by E3SM) to linearize PBL winds

    integer, intent(in) :: &
      order_xp2_xpyp, &
      order_wp2_wp3, &
      order_windm

    ! ------------------------ Input/Output Variables ------------------------
    type (stats), dimension(ngrdcol), target, intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc
      
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(inout) ::  &
      um,   & ! Mean u (west-to-east) wind component         [m/s]
      vm,   & ! Mean v (south-to-north) wind component       [m/s]
      upwp, & ! <u'w'> (momentum levels)                     [m^2/s^2]
      vpwp    ! <v'w'> (momentum levels)                     [m^2/s^2]

    ! Input/Output Variable for eddy-scalars
    real( kind = core_rknd ), dimension(ngrdcol,nz,edsclr_dim), intent(inout) ::  &
      edsclrm        ! Mean eddy scalar quantity             [units vary]

    ! Output Variable for eddy-scalars
    real( kind = core_rknd ), dimension(ngrdcol,nz,edsclr_dim), intent(inout) ::  &
      wpedsclrp      ! w'edsclr' (momentum levels)           [m/s {units vary}]

    ! Variables used to track perturbed version of winds.
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(inout) :: &
      um_pert,   & ! perturbed <u>    [m/s]
      vm_pert,   & ! perturbed <v>    [m/s]
      upwp_pert, & ! perturbed <u'w'> [m^2/s^2]
      vpwp_pert    ! perturbed <v'w'> [m^2/s^2]

    ! ------------------------ Local Variables ------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz) ::  &
      um_old, & ! Saved value of mean u (west-to-east) wind component    [m/s]
      vm_old    ! Saved value of Mean v (south-to-north) wind component  [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,nz,edsclr_dim) ::  &
      edsclrm_old    ! Saved value of mean eddy scalar quantity   [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nz) ::  &
      um_tndcy,    & ! u wind component tendency                    [m/s^2]
      vm_tndcy       ! v wind component tendency                    [m/s^2]

    real( kind = core_rknd ), dimension(ngrdcol,nz) ::  &
      upwp_chnge,  & ! Net change of u'w' due to clipping           [m^2/s^2]
      vpwp_chnge     ! Net change of v'w' due to clipping           [m^2/s^2]

    real( kind = core_rknd ), dimension(3,ngrdcol,nz) :: &
      lhs ! The implicit part of the tridiagonal matrix             [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nz,max(2,edsclr_dim)) :: &
      rhs,     &! The explicit part of the tridiagonal matrix       [units vary]
      solution  ! The solution to the tridiagonal matrix            [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      wind_speed,      & ! wind speed; sqrt(u^2 + v^2)              [m/s]
      wind_speed_pert    ! perturbed wind speed; sprt(u^2 + v^2)    [m/s]

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      u_star_sqd,      & ! Surface friction velocity, u_star, squared      [m/s]
      u_star_sqd_pert    ! perturbed u_star, squared                       [m/s]

    logical :: &
      l_imp_sfc_momentum_flux  ! Flag for implicit momentum surface fluxes.

    integer :: nrhs  ! Number of right hand side terms

    integer :: i, k, edsclr  ! Array index

    logical :: l_first_clip_ts, l_last_clip_ts ! flags for clip_covar
    
    real( kind = core_rknd ), dimension(ngrdcol) :: &
      nu10, &
      nu_zero
      
    real( kind = core_rknd ), dimension(3,ngrdcol,nz) :: &
      lhs_diff, & ! LHS diffustion terms
      lhs_ma_zt   ! LHS mean advection terms
        
    real( kind = core_rknd ), dimension(ngrdcol,nz) ::  &
      Km_zt,          & ! Eddy diffusivity of winds on momentum levs. [m^2/s]
      Kmh_zt,         & ! Eddy diffusivity of themo on momentum levs. [m^s/s]
      Km_zm_p_nu10,   & ! Km_zm plus nu_vert_res_dep%nu10
      xpwp              ! x'w' for arbitrary x
      
    integer, parameter :: &
      kp1_tdiag = 1, & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2, & ! Thermodynamic main diagonal index.
      km1_tdiag = 3    ! Thermodynamic subdiagonal index.

    ! Whether perturbed winds are being solved.
    logical :: l_perturbed_wind

    ! ------------------------ Begin Code ------------------------
    do i = 1, ngrdcol
      nu10(i) = nu_vert_res_dep(i)%nu10
      nu_zero(i) = zero
    end do
    
    do k = 1, nz
      do i = 1, ngrdcol
        Km_zm_p_nu10(i,:) = Km_zm(i,:) + nu10(i)
      end do
    end do

    l_perturbed_wind = ( .not. l_predict_upwp_vpwp ) .and. l_linearize_pbl_winds
    
    if ( .not. l_implemented ) then
      do i = 1, ngrdcol
        call term_ma_zt_lhs( gr(i), wm_zt(i,:),                       & ! intent(in)
                             gr(i)%invrs_dzt(:), gr(i)%invrs_dzm(:),  & ! intent(in)
                             l_upwind_xm_ma,                          & ! intent(in)
                             lhs_ma_zt(:,i,:) )                         ! intent(out)
      end do
    else
      lhs_ma_zt(:,:,:) = zero
    end if

    if ( .not. l_predict_upwp_vpwp ) then
      
      Km_zt(:,:) = max( zm2zt( nz, ngrdcol, gr, Km_zm(:,:) ), zero )

      ! Calculate diffusion terms
      do i = 1, ngrdcol
        call diffusion_zt_lhs( gr(i), Km_zm(i,:), Km_zt(i,:), nu10(i),  & ! intent(in)
                               gr(i)%invrs_dzm(:), gr(i)%invrs_dzt(:),  & ! intent(in)
                               invrs_rho_ds_zt(i,:), rho_ds_zm(i,:),    & ! intent(in)
                               lhs_diff(:,i,:) )                          ! intent(out)
      end do
      
      ! The lower boundary condition needs to be applied here at level 2.
      if ( .not. l_upwind_Kh_dp_term ) then 
        
        do i = 1, ngrdcol
          ! Thermodynamic superdiagonal: [ x xm(k+1,<t+1>) ]
          lhs_diff(kp1_tdiag,i,2) = - gr(i)%invrs_dzt(2) * invrs_rho_ds_zt(i,2) &
                                    * ( Km_zm(i,2) + nu10(i) ) * rho_ds_zm(i,2) * gr(i)%invrs_dzm(2)

          ! Thermodynamic main diagonal: [ x xm(k,<t+1>) ]
          lhs_diff(k_tdiag,i,2) = + gr(i)%invrs_dzt(2) * invrs_rho_ds_zt(i,2) &
                                  * ( Km_zm(i,2) + nu10(i) ) * rho_ds_zm(i,2) * gr(i)%invrs_dzm(2)

          ! Thermodynamic subdiagonal: [ x xm(k-1,<t+1>) ]
          lhs_diff(km1_tdiag,i,2) = zero
        end do

      else

        do i = 1, ngrdcol
          ! Thermodynamic superdiagonal: [ x xm(k+1,<t+1>) ]
          lhs_diff(kp1_tdiag,i,2) = - gr(i)%invrs_dzt(2) &
                                      * ( Km_zt(i,2) + nu10(i) ) * gr(i)%invrs_dzm(2)

          ! Thermodynamic main diagonal: [ x xm(k,<t+1>) ]
          lhs_diff(k_tdiag,i,2) = + gr(i)%invrs_dzt(2) &
                                    * ( Km_zt(i,2) + nu10(i) ) * gr(i)%invrs_dzm(2)

          ! Thermodynamic subdiagonal: [ x xm(k-1,<t+1>) ]
          lhs_diff(km1_tdiag,i,2) = zero
        end do

      end if

      if ( l_lmm_stepping ) then
        um_old(:,:) = um(:,:)
        vm_old(:,:) = vm(:,:)
      end if ! l_lmm_stepping

      !----------------------------------------------------------------
      ! Prepare tridiagonal system for horizontal winds, um and vm
      !----------------------------------------------------------------

      ! Compute Coriolis, geostrophic, and other prescribed wind forcings for um.
      call compute_uv_tndcy( nz, ngrdcol, gr, windm_edsclrm_um, & ! intent(in)
                             fcor, vm, vg,                      & ! intent(in)
                             um_forcing, l_implemented,         & ! intent(in)
                             stats_zt,                          & ! intent(inout)
                             um_tndcy )                           ! intent(out)

      ! Compute Coriolis, geostrophic, and other prescribed wind forcings for vm.
      call compute_uv_tndcy( nz, ngrdcol, gr, windm_edsclrm_vm, & ! intent(in)
                             fcor, um, ug,                      & ! intent(in)
                             vm_forcing, l_implemented,         & ! intent(in)
                             stats_zt,                          & ! intent(inout) 
                             vm_tndcy )                           ! intent(out)

      ! Momentum surface fluxes, u'w'|_sfc and v'w'|_sfc, are applied through
      ! an implicit method, such that:
      !    x'w'|_sfc = - ( u_star(t)^2 / wind_speed(t) ) * xm(t+1).
      l_imp_sfc_momentum_flux = .true.
      ! Compute wind speed (use threshold "eps" to prevent divide-by-zero
      ! error).
      wind_speed(:,:) = max( sqrt( um(:,:)**2 + vm(:,:)**2 ), eps )
      ! Compute u_star_sqd according to the definition of u_star.
      u_star_sqd(:) = sqrt( upwp(:,1)**2 + vpwp(:,1)**2 )

      ! Compute the explicit portion of the um equation.
      ! Build the right-hand side vector.
      call windm_edsclrm_rhs( nz, ngrdcol, gr, windm_edsclrm_um, dt,  & ! intent(in)
                              lhs_diff, um, um_tndcy,                 & ! intent(in)
                              rho_ds_zm, invrs_rho_ds_zt,             & ! intent(in)
                              l_imp_sfc_momentum_flux, upwp(:,1),     & ! intent(in)
                              stats_zt,                               & ! intent(inout)
                              rhs(:,:,windm_edsclrm_um) )               ! intent(out)

      ! Compute the explicit portion of the vm equation.
      ! Build the right-hand side vector.
      call windm_edsclrm_rhs( nz, ngrdcol, gr, windm_edsclrm_vm, dt,  & ! intent(in)
                              lhs_diff, vm, vm_tndcy,                 & ! intent(in)
                              rho_ds_zm, invrs_rho_ds_zt,             & ! intent(in)
                              l_imp_sfc_momentum_flux, vpwp(:,1),     & ! intent(in)
                              stats_zt,                               & ! intent(inout)
                              rhs(:,:,windm_edsclrm_vm) )               ! intent(out)

      ! Store momentum flux (explicit component)

      ! The surface flux, x'w'(1) = x'w'|_sfc, is set elsewhere in the model.
      ! upwp(1) = upwp_sfc
      ! vpwp(1) = vpwp_sfc
      
      call calc_xpwp( nz, ngrdcol, gr, &
                      Km_zm_p_nu10, um, &
                      xpwp )

      ! Solve for x'w' at all intermediate model levels.
      ! A Crank-Nicholson timestep is used.
      upwp(:,2:nz-1) = -one_half * xpwp(:,2:nz-1)
      
      call calc_xpwp( nz, ngrdcol, gr, &
                      Km_zm_p_nu10, vm, &
                      xpwp )
                                             
      vpwp(:,2:nz-1) = -one_half * xpwp(:,2:nz-1)

      ! A zero-flux boundary condition at the top of the model, d(xm)/dz = 0,
      ! means that x'w' at the top model level is 0,
      ! since x'w' = - K_zm * d(xm)/dz.
      upwp(:,nz) = zero
      vpwp(:,nz) = zero
      

      ! Compute the implicit portion of the um and vm equations.
      ! Build the left-hand side matrix.
      call windm_edsclrm_lhs( nz, ngrdcol, gr, dt,                    & ! intent(in)
                              lhs_ma_zt, lhs_diff,                    & ! intent(in)
                              wind_speed, u_star_sqd,                 & ! intent(in)
                              rho_ds_zm, invrs_rho_ds_zt,             & ! intent(in)
                              l_implemented, l_imp_sfc_momentum_flux, & ! intent(in)
                              lhs )                                     ! intent(out)

      ! Decompose and back substitute for um and vm
      nrhs = 2
      call windm_edsclrm_solve( nz, ngrdcol, gr, nrhs, iwindm_matrix_condt_num, & ! intent(in)
                                stats_sfc, &                                      ! intent(inout)
                                lhs, rhs, &                                       ! intent(inout)
                                solution )                                        ! intent(out)

      ! Check for singular matrices and bad LAPACK arguments
      if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code == clubb_fatal_error ) then
          write(fstderr,*) "Fatal error solving for um/vm"
          return
        end if
      end if

      !----------------------------------------------------------------
      ! Update zonal (west-to-east) component of mean wind, um
      !----------------------------------------------------------------
      um(:,1:nz) = solution(:,1:nz,windm_edsclrm_um)

      !----------------------------------------------------------------
      ! Update meridional (south-to-north) component of mean wind, vm
      !----------------------------------------------------------------
      vm(:,1:nz) = solution(:,1:nz,windm_edsclrm_vm)

      if ( l_stats_samp ) then
        do i = 1, ngrdcol
          ! Implicit contributions to um and vm
          call windm_edsclrm_implicit_stats( gr(i), windm_edsclrm_um, um(i,:),    & ! intent(in)
                                             lhs_diff(:,i,:), lhs_ma_zt(:,i,:),   & ! intent(in)
                                             invrs_rho_ds_zt(i,:), u_star_sqd(i), & ! intent(in)
                                             rho_ds_zm(i,:), wind_speed(i,:),     & ! intent(in)
                                             l_imp_sfc_momentum_flux,             & ! intent(in)
                                             stats_zt(i) )                          ! intent(inout)

          call windm_edsclrm_implicit_stats( gr(i), windm_edsclrm_vm, vm(i,:),    & ! intent(in)
                                             lhs_diff(:,i,:), lhs_ma_zt(:,i,:),   & ! intent(in)
                                             invrs_rho_ds_zt(i,:), u_star_sqd(i), & ! intent(in)
                                             rho_ds_zm(i,:), wind_speed(i,:),     & ! intent(in)
                                             l_imp_sfc_momentum_flux,             & ! intent(in)
                                             stats_zt(i) )                          ! intent(inout)
        end do
      end if ! l_stats_samp
  
      ! The values of um(1) and vm(1) are located below the model surface and
      ! do not affect the rest of the model.  The values of um(1) or vm(1) are
      ! simply set to the values of um(2) and vm(2), respectively, after the
      ! equation matrices has been solved.  Even though um and vm would sharply
      ! decrease to a value of 0 at the surface, this is done to avoid
      ! confusion on plots of the vertical profiles of um and vm.
      um(:,1) = um(:,2)
      vm(:,1) = vm(:,2)

      if ( l_lmm_stepping ) then
        um(:,:) = one_half * ( um_old(:,:) + um(:,:) )
        vm(:,:) = one_half * ( vm_old(:,:) + vm(:,:) )
      endif ! l_lmm_stepping ) then

      if ( uv_sponge_damp_settings%l_sponge_damping ) then
        
        ! _sponge_damp_settings and _sponge_damp_profile could potentially vary
        ! from column to column, but there is no column index in those variables.
        ! Thus this code is potentially unsafe when implemented in a host model, 
        ! which is indicated by l_implemented = T
        if ( l_implemented ) then
          write(fstderr,*) "l_sponge_damping = T and l_implemented = T &
                            -- this is likely unsafe and considered fatal"
          err_code = clubb_fatal_error
          return
        end if
          
        if ( l_stats_samp ) then
          do i = 1, ngrdcol
            call stat_begin_update( gr(i), ium_sdmp, um(i,:) / dt, & ! intent(in)
                                    stats_zt(i) )           ! intent(inout)
            call stat_begin_update( gr(i), ivm_sdmp, vm(i,:) / dt, & ! intent(in)
                                    stats_zt(i) )           ! intent(inout)
          end do
        end if

        do i = 1, ngrdcol
          um(i,1:nz) = sponge_damp_xm( gr(i), dt, gr(i)%zt, um_ref(i,1:nz), &
                                        um(i,1:nz), uv_sponge_damp_profile )
        end do

        do i = 1, ngrdcol
          vm(i,1:nz) = sponge_damp_xm( gr(i), dt, gr(i)%zt, vm_ref(i,1:nz), &
                                        vm(i,1:nz), uv_sponge_damp_profile )
        end do

        if ( l_stats_samp ) then
          do i = 1, ngrdcol
            call stat_end_update( gr(i), ium_sdmp, um(i,:) / dt, & ! intent(in) 
                                  stats_zt(i) )           ! intent(inout)
            call stat_end_update( gr(i), ivm_sdmp, vm(i,:) / dt, & ! intent(in)
                                  stats_zt(i) )           ! intent(inout)
          end do
        end if

      end if ! uv_sponge_damp_settings%l_sponge_damping

      ! Second part of momentum (implicit component)

      ! Solve for x'w' at all intermediate model levels.
      ! A Crank-Nicholson timestep is used.
      call calc_xpwp( nz, ngrdcol, gr, &
                      Km_zm_p_nu10, um, &
                      xpwp )
                      
      upwp(:,2:nz-1) = upwp(:,2:nz-1) - one_half * xpwp(:,2:nz-1)
                                                             
      call calc_xpwp( nz, ngrdcol, gr, &
                      Km_zm_p_nu10, vm, &
                      xpwp )
                      
      vpwp(:,2:nz-1) = vpwp(:,2:nz-1) - one_half * xpwp(:,2:nz-1)

      ! Adjust um and vm if nudging is turned on.
      if ( l_uv_nudge ) then

        ! Reflect nudging in budget
        if ( l_stats_samp ) then
          do i = 1, ngrdcol
            call stat_begin_update( gr(i), ium_ndg, um(i,:) / dt, & ! intent(in)
                                    stats_zt(i) )          ! intent(inout)
            call stat_begin_update( gr(i), ivm_ndg, vm(i,:) / dt, & ! intent(in)
                                    stats_zt(i) )          ! intent(inout)
          end do
        end if
    
        um(:,:) = um(:,:) - ( ( um(:,:) - um_ref(:,:) ) * (dt/ts_nudge) )
        vm(:,:) = vm(:,:) - ( ( vm(:,:) - vm_ref(:,:) ) * (dt/ts_nudge) )

        if ( l_stats_samp ) then
          do i = 1, ngrdcol
            call stat_end_update( gr(i), ium_ndg, um(i,:) / dt, & ! intent(in)
                                  stats_zt(i) )          ! intent(inout)
            call stat_end_update( gr(i), ivm_ndg, vm(i,:) / dt, & ! intent(in)
                                  stats_zt(i) )          ! intent(inout)
          end do
        end if
    
      end if ! l_uv_nudge

      if ( l_stats_samp ) then
        do i = 1, ngrdcol
          call stat_update_var( ium_ref, um_ref(i,:), & ! intent(in)
                                stats_zt(i) )         ! intent(inout)
          call stat_update_var( ivm_ref, vm_ref(i,:), & ! intent(in)
                                stats_zt(i) )         ! intent(inout)
        end do
      end if

      if ( order_windm < order_wp2_wp3 &
          .and. order_windm < order_xp2_xpyp ) then
        l_first_clip_ts = .true.
        l_last_clip_ts = .false.
      elseif ( order_windm > order_wp2_wp3 &
              .and. order_windm > order_xp2_xpyp ) then
        l_first_clip_ts = .false.
        l_last_clip_ts = .true.
      else
        l_first_clip_ts = .false.
        l_last_clip_ts = .false.
      endif

      if ( l_tke_aniso ) then

        ! Clipping for u'w'
        !
        ! Clipping u'w' at each vertical level, based on the
        ! correlation of u and w at each vertical level, such that:
        ! corr_(u,w) = u'w' / [ sqrt(u'^2) * sqrt(w'^2) ];
        ! -1 <= corr_(u,w) <= 1.
        !
        ! Since u'^2, w'^2, and u'w' are each advanced in different
        ! subroutines from each other in advance_clubb_core, clipping for u'w'
        ! has to be done three times during each timestep (once after each
        ! variable has been updated).
        ! This is the third instance of u'w' clipping.
        !l_first_clip_ts = .false.
        !l_last_clip_ts = .true.
        do i = 1, ngrdcol
          call clip_covar( gr(i), clip_upwp, l_first_clip_ts,       & ! intent(in)
                           l_last_clip_ts, dt, wp2(i,:), up2(i,:),  & ! intent(in)
                           l_predict_upwp_vpwp,                     & ! intent(in)
                           stats_zm(i),                             & ! intent(inout)
                           upwp(i,:), upwp_chnge(i,:) )               ! intent(inout)
        end do

        ! Clipping for v'w'
        !
        ! Clipping v'w' at each vertical level, based on the
        ! correlation of v and w at each vertical level, such that:
        ! corr_(v,w) = v'w' / [ sqrt(v'^2) * sqrt(w'^2) ];
        ! -1 <= corr_(v,w) <= 1.
        !
        ! Since v'^2, w'^2, and v'w' are each advanced in different
        ! subroutines from each other in advance_clubb_core, clipping for v'w'
        ! has to be done three times during each timestep (once after each
        ! variable has been updated).
        ! This is the third instance of v'w' clipping.
        !l_first_clip_ts = .false.
        !l_last_clip_ts = .true.
        do i = 1, ngrdcol
          call clip_covar( gr(i), clip_vpwp, l_first_clip_ts,       & ! intent(in)
                           l_last_clip_ts, dt, wp2(i,:), vp2(i,:),  & ! intent(in)
                           l_predict_upwp_vpwp,                     & ! intent(in)
                           stats_zm(i),                             & ! intent(inout)
                           vpwp(i,:), vpwp_chnge(i,:) )               ! intent(inout)
        end do
      else

        ! intent(in) this case, it is assumed that
        !   u'^2 == v'^2 == w'^2, and the variables `up2' and `vp2' do not
        ! interact with any other variables.
        !l_first_clip_ts = .false.
        !l_last_clip_ts = .true.
        do i = 1, ngrdcol
          call clip_covar( gr(i), clip_upwp, l_first_clip_ts,       & ! intent(in)
                           l_last_clip_ts, dt, wp2(i,:), wp2(i,:),  & ! intent(in)
                           l_predict_upwp_vpwp,                     & ! intent(in)
                           stats_zm(i),                             & ! intent(inout)
                           upwp(i,:), upwp_chnge(i,:) )               ! intent(inout)
        end do

        do i = 1, ngrdcol
          call clip_covar( gr(i), clip_vpwp, l_first_clip_ts,       & ! intent(in)
                           l_last_clip_ts, dt, wp2(i,:), wp2(i,:),  & ! intent(in)
                           l_predict_upwp_vpwp,                     & ! intent(in)
                           stats_zm(i),                             & ! intent(inout)
                           vpwp(i,:), vpwp_chnge(i,:) )               ! intent(inout)
        end do
      endif ! l_tke_aniso

    endif ! .not. l_predict_upwp_vpwp

    
    if ( l_perturbed_wind ) then

      !----------------------------------------------------------------
      ! Prepare tridiagonal system for horizontal winds, um and vm
      !----------------------------------------------------------------

      ! Momentum surface fluxes, u'w'|_sfc and v'w'|_sfc, are applied through
      ! an implicit method, such that:
      !    x'w'|_sfc = - ( u_star(t)^2 / wind_speed(t) ) * xm(t+1).
      l_imp_sfc_momentum_flux = .true.
      ! Compute wind speed (use threshold "eps" to prevent divide-by-zero
      ! error).
      wind_speed_pert(:,:) = max( sqrt( (um_pert(:,:))**2 + (vm_pert(:,:))**2 ), eps )
      ! Compute u_star_sqd according to the definition of u_star.
      u_star_sqd_pert(:) = sqrt( upwp_pert(:,1)**2 + vpwp_pert(:,1)**2 )

      ! Compute the explicit portion of the um equation.
      ! Build the right-hand side vector.
      call windm_edsclrm_rhs( nz, ngrdcol, gr, windm_edsclrm_um, dt,    & ! intent(in)
                              lhs_diff, um_pert, um_tndcy,              & ! intent(in)
                              rho_ds_zm, invrs_rho_ds_zt,               & ! intent(in)
                              l_imp_sfc_momentum_flux, upwp_pert(:,1),  & ! intent(in)
                              stats_zt,                                 & ! intent(inout)
                              rhs(:,:,windm_edsclrm_um) )                 ! intent(out)
      
      ! Compute the explicit portion of the vm equation.
      ! Build the right-hand side vector.
      call windm_edsclrm_rhs( nz, ngrdcol, gr, windm_edsclrm_vm, dt,    & ! intent(in)
                              lhs_diff, vm_pert, vm_tndcy,              & ! intent(in)
                              rho_ds_zm, invrs_rho_ds_zt,               & ! intent(in)
                              l_imp_sfc_momentum_flux, vpwp_pert(:,1),  & ! intent(in)
                              stats_zt,                                 & ! intent(inout)
                              rhs(:,:,windm_edsclrm_vm) )                 ! intent(out)

      ! Store momentum flux (explicit component)

      ! The surface flux, x'w'(1) = x'w'|_sfc, is set elsewhere in the model.
      !      upwp(1) = upwp_sfc
      !      vpwp(1) = vpwp_sfc

      ! Solve for x'w' at all intermediate model levels.
      ! A Crank-Nicholson timestep is used.
      call calc_xpwp( nz, ngrdcol, gr, &
                      Km_zm_p_nu10, um_pert, &
                      xpwp )
                      
      upwp_pert(:,2:nz-1) = -one_half * xpwp(:,2:nz-1)
                                                  
      call calc_xpwp( nz, ngrdcol, gr, &
                      Km_zm_p_nu10, vm_pert, &
                      xpwp )
                                                  
      vpwp_pert(:,2:nz-1) = -one_half * xpwp(:,2:nz-1)
      
      ! A zero-flux boundary condition at the top of the model, d(xm)/dz = 0,
      ! means that x'w' at the top model level is 0,
      ! since x'w' = - K_zm * d(xm)/dz.
      upwp_pert(:,nz) = zero
      vpwp_pert(:,nz) = zero

      ! Compute the implicit portion of the um and vm equations.
      ! Build the left-hand side matrix.
      call windm_edsclrm_lhs( nz, ngrdcol, gr, dt,                        & ! intent(in)
                              lhs_ma_zt, lhs_diff,                        & ! intent(in)
                              wind_speed_pert, u_star_sqd_pert,           & ! intent(in)
                              rho_ds_zm, invrs_rho_ds_zt,                 & ! intent(in)
                              l_implemented, l_imp_sfc_momentum_flux,     & ! intent(in)
                              lhs )                                         ! intent(out)
      
      ! Decompose and back substitute for um and vm
      nrhs = 2
      call windm_edsclrm_solve( nz, ngrdcol, gr, nrhs, iwindm_matrix_condt_num, & ! intent(in)
                                stats_sfc, &                                      ! intent(in)
                                lhs, rhs, &                                       ! intent(inout)
                                solution )                                        ! intent(out)
      
      ! Check for singular matrices and bad LAPACK arguments
      if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code == clubb_fatal_error ) then
          write(fstderr,*) "Fatal error solving for um_pert/vm_pert"
          return
        endif
      endif

      !----------------------------------------------------------------
      ! Update zonal (west-to-east) component of mean wind, um
      !----------------------------------------------------------------
      um_pert(:,1:nz) = solution(:,1:nz,windm_edsclrm_um)

      !----------------------------------------------------------------
      ! Update meridional (south-to-north) component of mean wind, vm
      !----------------------------------------------------------------
      vm_pert(:,1:nz) = solution(:,1:nz,windm_edsclrm_vm)
  
      ! The values of um(1) and vm(1) are located below the model surface and
      ! do not affect the rest of the model.  The values of um(1) or vm(1) are
      ! simply set to the values of um(2) and vm(2), respectively, after the
      ! equation matrices has been solved.  Even though um and vm would sharply
      ! decrease to a value of 0 at the surface, this is done to avoid
      ! confusion on plots of the vertical profiles of um and vm.
      um_pert(:,1) = um_pert(:,2)
      vm_pert(:,1) = vm_pert(:,2)

      ! Second part of momentum (implicit component)

      ! Solve for x'w' at all intermediate model levels.
      ! A Crank-Nicholson timestep is used.
      call calc_xpwp( nz, ngrdcol, gr, &
                      Km_zm_p_nu10, um_pert, &
                      xpwp )
                      
      upwp_pert(:,2:nz-1) = upwp_pert(:,2:nz-1) - one_half * xpwp(:,2:nz-1)
                                                          
      call calc_xpwp( nz, ngrdcol, gr, &
                      Km_zm_p_nu10, vm_pert, &
                      xpwp )
                         
      vpwp_pert(:,2:nz-1) = vpwp_pert(:,2:nz-1) - one_half * xpwp(:,2:nz-1)
                            
      if ( l_tke_aniso ) then

        ! Clipping for u'w'
        !
        ! Clipping u'w' at each vertical level, based on the
        ! correlation of u and w at each vertical level, such that:
        ! corr_(u,w) = u'w' / [ sqrt(u'^2) * sqrt(w'^2) ];
        ! -1 <= corr_(u,w) <= 1.
        !
        ! Since u'^2, w'^2, and u'w' are each advanced in different
        ! subroutines from each other in advance_clubb_core, clipping for u'w'
        ! has to be done three times during each timestep (once after each
        ! variable has been updated).
        ! This is the third instance of u'w' clipping.
        l_first_clip_ts = .false.
        l_last_clip_ts = .true.
        do i = 1, ngrdcol
          call clip_covar( gr(i), clip_upwp, l_first_clip_ts,       & ! intent(in)
                           l_last_clip_ts, dt, wp2(i,:), up2(i,:),  & ! intent(in)
                           l_predict_upwp_vpwp,                     & ! intent(in)
                           stats_zm(i),                             & ! intent(inout)
                           upwp_pert(i,:), upwp_chnge(i,:) )          ! intent(inout)
        end do
        
        ! Clipping for v'w'
        !
        ! Clipping v'w' at each vertical level, based on the
        ! correlation of v and w at each vertical level, such that:
        ! corr_(v,w) = v'w' / [ sqrt(v'^2) * sqrt(w'^2) ];
        ! -1 <= corr_(v,w) <= 1.
        !
        ! Since v'^2, w'^2, and v'w' are each advanced in different
        ! subroutines from each other in advance_clubb_core, clipping for v'w'
        ! has to be done three times during each timestep (once after each
        ! variable has been updated).
        ! This is the third instance of v'w' clipping.
        l_first_clip_ts = .false.
        l_last_clip_ts = .true.
        do i = 1, ngrdcol
          call clip_covar( gr(i), clip_vpwp, l_first_clip_ts,       & ! intent(in)
                           l_last_clip_ts, dt, wp2(i,:), vp2(i,:),  & ! intent(in)
                           l_predict_upwp_vpwp,                     & ! intent(in)
                           stats_zm(i),                             & ! intent(inout)
                           vpwp_pert(i,:), vpwp_chnge(i,:) )          ! intent(inout)
        end do
      else

        ! intent(in) this case, it is assumed that
        !   u'^2 == v'^2 == w'^2, and the variables `up2' and `vp2' do not
        ! interact with any other variables.
        l_first_clip_ts = .false.
        l_last_clip_ts = .true.
        do i = 1, ngrdcol
          call clip_covar( gr(i), clip_upwp, l_first_clip_ts,         & ! intent(in)
                           l_last_clip_ts, dt, wp2(i,:), wp2(i,:),    & ! intent(in)
                           l_predict_upwp_vpwp,                       & ! intent(in)
                           stats_zm(i),                               & ! intent(inout)
                           upwp_pert(i,:), upwp_chnge(i,:) )            ! intent(inout)
        end do
        
        do i = 1, ngrdcol
          call clip_covar( gr(i), clip_vpwp, l_first_clip_ts,       & ! intent(in)
                           l_last_clip_ts, dt, wp2(i,:), wp2(i,:),  & ! intent(in)
                           l_predict_upwp_vpwp,                     & ! intent(in)
                           stats_zm(i),                             & ! intent(inout)
                           vpwp_pert(i,:), vpwp_chnge(i,:) )          ! intent(inout)
        end do
        
      end if ! l_tke_aniso
    end if ! .not. l_predict_upwp_vpwp

    !----------------------------------------------------------------
    ! Prepare tridiagonal system for eddy-scalars
    !----------------------------------------------------------------

    if ( edsclr_dim > 0 ) then
      
      Kmh_zt(:,:) = max( zm2zt( nz, ngrdcol, gr, Kmh_zm(:,:) ), zero )

      ! Calculate diffusion terms
      do i = 1, ngrdcol
        call diffusion_zt_lhs( gr(i), Kmh_zm(i,:), Kmh_zt(i,:), nu_zero(i), & ! intent(in)
                               gr(i)%invrs_dzm(:), gr(i)%invrs_dzt(:),      & ! intent(in)
                               invrs_rho_ds_zt(i,:), rho_ds_zm(i,:),        & ! intent(in)
                               lhs_diff(:,i,:) )                              ! intent(out)
      end do
      
      ! The lower boundary condition needs to be applied here at level 2.
      if ( .not. l_upwind_Kh_dp_term ) then 
        
        do i = 1, ngrdcol
          ! Thermodynamic superdiagonal: [ x xm(k+1,<t+1>) ]
          lhs_diff(kp1_tdiag,i,2) = - gr(i)%invrs_dzt(2) * invrs_rho_ds_zt(i,2) &
                                    * ( Kmh_zm(i,2) + nu_zero(i) ) &
                                    * rho_ds_zm(i,2) * gr(i)%invrs_dzm(2)

          ! Thermodynamic main diagonal: [ x xm(k,<t+1>) ]
          lhs_diff(k_tdiag,i,2) = + gr(i)%invrs_dzt(2) * invrs_rho_ds_zt(i,2) &
                                  * ( Kmh_zm(i,2) + nu_zero(i) ) &
                                  * rho_ds_zm(i,2) * gr(i)%invrs_dzm(2)

          ! Thermodynamic subdiagonal: [ x xm(k-1,<t+1>) ]
          lhs_diff(km1_tdiag,i,2) = zero
        end do

      else
        do i = 1, ngrdcol
          ! Thermodynamic superdiagonal: [ x xm(k+1,<t+1>) ]
          lhs_diff(kp1_tdiag,i,2) = - gr(i)%invrs_dzt(2) &
                                      * ( Kmh_zt(i,2) + nu_zero(i) ) * gr(i)%invrs_dzm(2)

          ! Thermodynamic main diagonal: [ x xm(k,<t+1>) ]
          lhs_diff(k_tdiag,i,2) = + gr(i)%invrs_dzt(2) &
                                    * ( Kmh_zt(i,2) + nu_zero(i) ) * gr(i)%invrs_dzm(2)

          ! Thermodynamic subdiagonal: [ x xm(k-1,<t+1>) ]
          lhs_diff(km1_tdiag,i,2) = zero
        end do

      end if

      if ( l_lmm_stepping ) then
        edsclrm_old(:,:,:) = edsclrm(:,:,:)
      endif ! l_lmm_stepping

      ! Eddy-scalar surface fluxes, x'w'|_sfc, are applied through an explicit
      ! method.
      l_imp_sfc_momentum_flux = .false.

      ! Compute the explicit portion of eddy scalar equation.
      ! Build the right-hand side vector.
      ! Because of statistics, we have to use a DO rather than a FORALL here
      ! -dschanen 7 Oct 2008
      do edsclr = 1, edsclr_dim
        call windm_edsclrm_rhs( nz, ngrdcol, gr, windm_edsclrm_scalar, dt,      & ! intent(in)
                                lhs_diff, edsclrm(:,:,edsclr),                  & ! intent(in)
                                edsclrm_forcing(:,:,edsclr),                    & ! intent(in)
                                rho_ds_zm, invrs_rho_ds_zt,                     & ! intent(in)
                                l_imp_sfc_momentum_flux, wpedsclrp(:,1,edsclr), & ! intent(in)
                                stats_zt,                                       & ! intent(inout)
                                rhs(:,:,edsclr) )                                 ! intent(out)
      enddo


      ! Store momentum flux (explicit component)

      ! The surface flux, x'w'(1) = x'w'|_sfc, is set elsewhere in the model.
      ! wpedsclrp(1,1:edsclr_dim) =  wpedsclrp_sfc(1:edsclr_dim)

      ! Solve for x'w' at all intermediate model levels.
      ! A Crank-Nicholson timestep is used
      do edsclr = 1, edsclr_dim
        
        call calc_xpwp( nz, ngrdcol, gr, &
                        Km_zm_p_nu10, edsclrm(:,:,edsclr), &
                        xpwp )
                        
        wpedsclrp(:,2:nz-1,edsclr) = -one_half * xpwp(:,2:nz-1)
      end do

      ! A zero-flux boundary condition at the top of the model, d(xm)/dz = 0,
      ! means that x'w' at the top model level is 0,
      ! since x'w' = - K_zm * d(xm)/dz.
      wpedsclrp(:,nz,1:edsclr_dim) = zero

      ! Compute the implicit portion of the xm (eddy-scalar) equations.
      ! Build the left-hand side matrix.
      call windm_edsclrm_lhs( nz, ngrdcol, gr, dt,                    & ! intent(in)
                              lhs_ma_zt, lhs_diff,                    & ! intent(in)
                              wind_speed, u_star_sqd,                 & ! intent(in)
                              rho_ds_zm, invrs_rho_ds_zt,             & ! intent(in)
                              l_implemented, l_imp_sfc_momentum_flux, & ! intent(in)
                              lhs )                                     ! intent(out)
                                    
      ! Decompose and back substitute for all eddy-scalar variables
      call windm_edsclrm_solve( nz, ngrdcol, gr, edsclr_dim, 0, & ! intent(in)
                                stats_sfc,                      & ! intent(inout)
                                lhs, rhs,                       & ! intent(inout)
                                solution )                        ! intent(out)
      
      if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code == clubb_fatal_error ) then
          write(fstderr,*) "Fatal error solving for eddsclrm"
        end if
      end if

      !----------------------------------------------------------------
      ! Update Eddy-diff. Passive Scalars
      !----------------------------------------------------------------
      edsclrm(:,:,:) = solution(:,:,:)

      ! The value of edsclrm(1) is located below the model surface and does not
      ! effect the rest of the model.  The value of edsclrm(1) is simply set to
      ! the value of edsclrm(2) after the equation matrix has been solved.
      do edsclr = 1, edsclr_dim
        do i = 1, ngrdcol
          edsclrm(i,1,edsclr) = edsclrm(i,2,edsclr)
        end do
      end do

      if ( l_lmm_stepping ) then
        edsclrm(:,:,:) = one_half * ( edsclrm_old(:,:,:) + edsclrm(:,:,:) )
      endif ! l_lmm_stepping

      ! Second part of momentum (implicit component)

      ! Solve for x'w' at all intermediate model levels.
      ! A Crank-Nicholson timestep is used.
      do edsclr = 1, edsclr_dim
        
        call calc_xpwp( nz, ngrdcol, gr, &
                        Kmh_zm, edsclrm(:,:,edsclr), &
                        xpwp )
                        
        wpedsclrp(:,2:nz-1,edsclr) = -one_half * xpwp(:,2:nz-1)
      end do

      ! Note that the w'edsclr' terms are not clipped, since we don't compute
      ! the variance of edsclr anywhere. -dschanen 7 Oct 2008

    endif
    
    if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code == clubb_fatal_error ) then

          write(fstderr,*) "Error in advance_windm_edsclrm"

          write(fstderr,*) "intent(in)"

          write(fstderr,*) "dt = ", dt
          write(fstderr,*) "wm_zt = ", wm_zt
          write(fstderr,*) "Km_zm = ", Km_zm
          write(fstderr,*) "ug = ", ug
          write(fstderr,*) "vg = ", vg
          write(fstderr,*) "um_ref = ", um_ref
          write(fstderr,*) "vm_ref = ", vm_ref
          write(fstderr,*) "wp2 = ", wp2
          write(fstderr,*) "up2 = ", up2
          write(fstderr,*) "vp2 = ", vp2
          write(fstderr,*) "um_forcing = ", um_forcing
          write(fstderr,*) "vm_forcing = ", vm_forcing
          do edsclr = 1, edsclr_dim
            write(fstderr,*) "edsclrm_forcing # = ",edsclr, edsclrm_forcing
          end do
          write(fstderr,*) "fcor = ", fcor
          write(fstderr,*) "l_implemented = ", l_implemented

          write(fstderr,*) "intent(inout)"

          if ( l_lmm_stepping ) &
             write(fstderr,*) "um (pre-solve) = ", um_old
          write(fstderr,*) "um = ", um
          if ( l_lmm_stepping ) &
             write(fstderr,*) "vm (pre-solve) = ", vm_old
          write(fstderr,*) "vm = ", vm
          do edsclr = 1, edsclr_dim
            if ( l_lmm_stepping ) &
               write(fstderr,*) "edsclrm (pre-solve) # ", edsclr, "=", edsclrm_old(:,:,edsclr)
            write(fstderr,*) "edsclrm # ", edsclr, "=", edsclrm(:,:,edsclr)
          end do
          write(fstderr,*) "upwp = ", upwp
          write(fstderr,*) "vpwp = ", vpwp
          write(fstderr,*) "wpedsclrp = ", wpedsclrp

          return
        end if
    end if

    return

  end subroutine advance_windm_edsclrm

  !=============================================================================
  subroutine windm_edsclrm_solve( nz, ngrdcol, gr, nrhs, ixm_matrix_condt_num, &
                                  stats_sfc, & 
                                  lhs, rhs, solution )

    ! Note:  In the "Description" section of this subroutine, the variable
    !        "invrs_dzm" will be written as simply "dzm", and the variable
    !        "invrs_dzt" will be written as simply "dzt".  This is being done as
    !        as device to save space and to make some parts of the description
    !        more readable.  This change does not pertain to the actual code.

    ! Description:
    ! Solves the horizontal wind or eddy-scalar time-tendency equation, and
    ! diagnoses the turbulent flux.  A Crank-Nicholson time-stepping algorithm
    ! is used in solving the turbulent advection term and in diagnosing the
    ! turbulent flux.
    !
    ! The rate of change of an eddy-scalar quantity, xm, is:
    !
    ! d(xm)/dt = - w * d(xm)/dz - (1/rho_ds) * d( rho_ds * x'w' )/dz 
    !            + xm_forcings.
    !
    !
    ! The Turbulent Advection Term
    ! ----------------------------
    !
    ! The above equation contains a turbulent advection term:
    !
    ! - (1/rho_ds) * d( rho_ds * x'w' )/dz;
    !
    ! where the momentum flux, x'w', is closed using a down gradient approach:
    !
    ! x'w' = - K_zm * d(xm)/dz.
    !
    ! The turbulent advection term becomes:
    !
    ! + (1/rho_ds) * d [ rho_ds * K_zm * d(xm)/dz ] / dz;
    !
    ! which is the same as a standard eddy-diffusion term (if "rho_ds * K_zm" in
    ! the term above is substituted for "K_zm" in a standard eddy-diffusion
    ! term, and if the standard eddy-diffusion term is multiplied by
    ! "1/rho_ds").  Thus, the turbulent advection term is treated and solved in
    ! the same way that a standard eddy-diffusion term would be solved.  The
    ! term is discretized as follows:
    !
    ! The values of xm are found on the thermodynamic levels, while the values
    ! of K_zm are found on the momentum levels.  Additionally, the values of
    ! rho_ds_zm are found on the momentum levels, and the values of
    ! invrs_rho_ds_zt are found on the thermodynamic levels.  The
    ! derivatives (d/dz) of xm are taken over the intermediate momentum levels.
    ! At the intermediate momentum levels, d(xm)/dz is multiplied by K_zm and by
    ! rho_ds_zm.  Then, the derivative of the whole mathematical expression is
    ! taken over the central thermodynamic level, where it is multiplied by
    ! invrs_rho_ds_zt, which yields the desired result.
    !
    ! ---xm(kp1)----------------------------------------------------- t(k+1)
    !
    ! ===========d(xm)/dz===K_zm(k)=====rho_ds_zm(k)================= m(k)
    !
    ! ---xm(k)---invrs_rho_ds_zt---d[rho_ds_zm*K_zm*d(xm)/dz]/dz----- t(k)
    !
    ! ===========d(xm)/dz===K_zm(km1)===rho_ds_zm(km1)=============== m(k-1)
    !
    ! ---xm(km1)----------------------------------------------------- t(k-1)
    !
    ! The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1) correspond
    ! with altitudes zt(k+1), zm(k), zt(k), zm(k-1), and zt(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is used
    ! for momentum levels.
    !
    ! dzt(k)   = 1 / ( zm(k) - zm(k-1) )
    ! dzm(k)   = 1 / ( zt(k+1) - zt(k) )
    ! dzm(k-1) = 1 / ( zt(k) - zt(k-1) )
    !
    ! The vertically discretized form of the turbulent advection term (treated
    ! as an eddy diffusion term) is written out as:
    !
    ! + invrs_rho_ds_zt(k)
    !   * dzt(k)
    !     * [   rho_ds_zm(k) * K_zm(k) * dzm(k) * ( xm(k+1) - xm(k) )
    !         - rho_ds_zm(k-1) * K_zm(k-1) * dzm(k-1) * ( xm(k) - xm(k-1) ) ].
    !
    ! For this equation, a Crank-Nicholson (semi-implicit) diffusion scheme is
    ! used to solve the (1/rho_ds) * d [ rho_ds * K_zm * d(xm)/dz ] / dz
    ! eddy-diffusion term.  The discretized implicit portion of the term is
    ! written out as:
    !
    ! + (1/2) * invrs_rho_ds_zt(k)
    !   * dzt(k)
    !     * [   rho_ds_zm(k) * K_zm(k)
    !           * dzm(k) * ( xm(k+1,<t+1>) - xm(k,<t+1>) )
    !         - rho_ds_zm(k-1) * K_zm(k-1)
    !           * dzm(k-1) * ( xm(k,<t+1>) - xm(k-1,<t+1>) ) ].
    !
    ! Note:  When the implicit term is brought over to the left-hand side,
    !        the sign is reversed and the leading "+" in front of the term
    !        is changed to a "-".
    !
    ! The discretized explicit portion of the term is written out as:
    !
    ! + (1/2) * invrs_rho_ds_zt(k)
    !   * dzt(k)
    !     * [   rho_ds_zm(k) * K_zm(k)
    !           * dzm(k) * ( xm(k+1,<t>) - xm(k,<t>) )
    !         - rho_ds_zm(k-1) * K_zm(k-1)
    !           * dzm(k-1) * ( xm(k,<t>) - xm(k-1,<t>) ) ].
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is
    ! being advanced to in solving the d(xm)/dt equation.
    !
    !
    ! Boundary Conditions:
    !
    ! An eddy-scalar quantity is not allowed to flux out the upper boundary.
    ! Thus, a zero-flux boundary condition is used for the upper boundary in the
    ! eddy-diffusion equation.
    !
    ! The lower boundary condition is much more complicated.  It is neither a
    ! zero-flux nor a fixed-point boundary condition.  Rather, it is a
    ! fixed-flux boundary condition.  This term is a turbulent advection term,
    ! but with the eddy-scalars, the only value of x'w' relevant in solving the
    ! d(xm)/dt equation is the value of x'w' at the surface (the first momentum
    ! level), which is written as x'w'|_sfc.
    !
    ! 1) x'w' surface flux; generalized explicit form
    !
    !    The x'w' surface flux is applied to the d(xm)/dt equation through the
    !    turbulent advection term, which is:
    !
    !    - (1/rho_ds) * d( rho_ds * x'w' )/dz.
    !
    !    At most vertical levels, a substitution can be made for x'w', such
    !    that:
    !
    !    x'w' = - K_zm * d(xm)/dz.
    !
    !    However, the same substitution cannot be made at the surface (momentum
    !    level 1), as x'w'|_sfc is a surface flux that is explicitly computed
    !    elsewhere in the model code.
    !
    !    The lower boundary condition, which in this case needs to be applied to
    !    the d(xm)/dt equation at level 2, is discretized as follows:
    !
    !    --xm(3)------------------------------------------------------- t(3)
    !
    !    ========[x'w'(2) = -K_zm(2)*d(xm)/dz]===rho_ds_zm(2)========== m(2)
    !
    !    --xm(2)---invrs_rho_ds_zt(2)---d[rho_ds_zm*K_zm*d(xm)/dz]/dz-- t(2)
    !
    !    ========[x'w'|_sfc]=====================rho_ds_zm(1)========== m(1) sfc
    !
    !    --xm(1)-------(below surface; not applicable)----------------- t(1)
    !
    !    where "sfc" is the level of the model surface or lower boundary.
    !
    !    The vertically discretized form of the turbulent advection term
    !    (treated as an eddy diffusion term), with the explicit surface flux,
    !    x'w'|_sfc, in place, is written out as:
    !
    !    - invrs_rho_ds_zt(2)
    !      * dzt(2) * [ rho_ds_zm(2) * x'w'(2) - rho_ds_zm(1) * x'w'|_sfc ];
    !
    !    which can be re-written as:
    !
    !    + invrs_rho_ds_zt(2)
    !      * dzt(2)
    !        * [   rho_ds_zm(2) * K_zm(2) * dzm(2) * ( xm(3) - xm(2) )
    !            + rho_ds_zm(1) * x'w'|_sfc ];
    !
    !    which can be re-written again as:
    !
    !    + invrs_rho_ds_zt(2)
    !      * dzt(2)
    !        * rho_ds_zm(2) * K_zm(2) * dzm(2) * ( xm(3) - xm(2) )
    !    + invrs_rho_ds_zt(2)
    !      * dzt(2)
    !        * rho_ds_zm(1) * x'w'|_sfc.
    !
    !    For this equation, a Crank-Nicholson (semi-implicit) diffusion scheme
    !    is used to solve the (1/rho_ds) * d [ rho_ds * K_zm * d(xm)/dz ] / dz
    !    eddy-diffusion term.  The discretized implicit portion of the term is
    !    written out as:
    !
    !    + (1/2) * invrs_rho_ds_zt(2)
    !      * dzt(2)
    !        * [ rho_ds_zm(2) * K_zm(2)
    !            * dzm(2) * ( xm(3,<t+1>) - xm(2,<t+1>) ) ].
    !
    !    Note:  When the implicit term is brought over to the left-hand side,
    !           the sign is reversed and the leading "+" in front of the term
    !           is changed to a "-".
    !
    !    The discretized explicit portion of the term is written out as:
    !
    !    + (1/2) * invrs_rho_ds_zt(2)
    !      * dzt(2)
    !        * [ rho_ds_zm(2) * K_zm(2)
    !            * dzm(2) * ( xm(3,<t>) - xm(2,<t>) ) ]
    !    + invrs_rho_ds_zt(2)
    !      * dzt(2)
    !        * rho_ds_zm(1) * x'w'|_sfc.
    !
    !    Note:  The x'w'|_sfc portion of the term written above has been pulled
    !           away from the rest of the explicit form written above because
    !           the (1/2) factor due to Crank-Nicholson time_stepping does not
    !           apply to it, as there isn't an implicit portion for x'w'|_sfc.
    !
    !    Timestep index (t) stands for the index of the current timestep, while
    !    timestep index (t+1) stands for the index of the next timestep, which
    !    is being advanced to in solving the d(xm)/dt equation.
    !
    ! 2) x'w' surface flux; implicit form for momentum fluxes u'w' and v'w'
    !
    !    The x'w' surface flux is applied to the d(xm)/dt equation through the
    !    turbulent advection term, which is:
    !
    !    - (1/rho_ds) * d( rho_ds * x'w' )/dz.
    !
    !    At most vertical levels, a substitution can be made for x'w', such
    !    that:
    !
    !    x'w' = - K_zm * d(xm)/dz.
    !
    !    However, the same substitution cannot be made at the surface (momentum
    !    level 1), as x'w'|_sfc is a surface momentum flux that is found by the
    !    following equation:
    !
    !    x'w'|_sfc = - [ u_star^2 / sqrt( um^2 + vm^2 ) ] * xm;
    !
    !    where x'w'|_sfc and xm are either u'w'|_sfc and um, respectively, or
    !    v'w'|_sfc and vm, respectively (um and vm are located at the first
    !    thermodynamic level above the surface, which is thermodynamic level 2),
    !    sqrt( um^2 + vm^2 ) is the wind speed (also at thermodynamic level 2),
    !    and u_star is defined as:
    !
    !    u_star = ( u'w'|_sfc^2 + v'w'|_sfc^2 )^(1/4);
    !
    !    and thus u_star^2 is defined as:
    !
    !    u_star^2 = sqrt( u'w'|_sfc^2 + v'w'|_sfc^2 ).
    !
    !    The value of u_star is either set to a constant value or computed
    !    (through function diag_ustar) based on the surface wind speed, the
    !    height above surface of the surface wind speed (as compared to the
    !    roughness height), and the buoyancy flux at the surface.  Either way,
    !    u_star is computed elsewhere in the model, and the values of u'w'|_sfc
    !    and v'w'|_sfc are based on it and computed along with it.  The values
    !    of u'w'|_sfc and v'w'|_sfc are then passed into advance_clubb_core,
    !    and are eventually passed into advance_windm_edsclrm.  In subroutine
    !    advance_windm_edsclrm, the value of u_star_sqd is then recomputed
    !    based on u'w'|_sfc and v'w'|_sfc.  The value of sqrt( u_star_sqd ) is
    !    consistent with the value of the original computation of u_star.
    !
    !    The equation listed above is substituted for x'w'|_sfc.  The lower
    !    boundary condition, which in this case needs to be applied to the
    !    d(xm)/dt equation at level 2, is discretized as follows:
    !
    !    --xm(3)------------------------------------------------------- t(3)
    !
    !    ===[x'w'(2) = -K_zm(2)*d(xm)/dz]=================rho_ds_zm(2)= m(2)
    !
    !    --xm(2)---invrs_rho_ds_zt(2)---d[rho_ds_zm*K_zm*d(xm)/dz]/dz-- t(2)
    !
    !    ===[x'w'|_sfc = -[u_star^2/sqrt(um^2+vm^2)]*xm]==rho_ds_zm(1)= m(1) sfc
    !
    !    --xm(1)-------(below surface; not applicable)----------------- t(1)
    !
    !    where "sfc" is the level of the model surface or lower boundary.
    !
    !    The vertically discretized form of the turbulent advection term
    !    (treated as an eddy diffusion term), with the implicit surface momentum
    !    flux in place, is written out as:
    !
    !    - invrs_rho_ds_zt(2)
    !      * dzt(2) * [ rho_ds_zm(2) * x'w'(2) - rho_ds_zm(1) * x'w'|_sfc ];
    !
    !    which can be re-written as:
    !
    !    - invrs_rho_ds_zt(2)
    !      * dzt(2)
    !        * [   rho_ds_zm(2)
    !              * { - K_zm(2) * dzm(2) * ( xm(3) - xm(2) ) }
    !            - rho_ds_zm(1)
    !              * { - [ u_star^2 / sqrt( um(2)^2 + vm(2)^2 ) ] * xm(2) } ];
    !
    !    which can be re-written as:
    !
    !    + invrs_rho_ds_zt(2)
    !      * dzt(2)
    !        * rho_ds_zm(2) * K_zm(2) * dzm(2) * ( xm(3) - xm(2) )
    !    - invrs_rho_ds_zt(2)
    !      * dzt(2)
    !        * rho_ds_zm(1) * [ u_star^2 / sqrt( um(2)^2 + vm(2)^2 ) ] * xm(2).
    !
    !    For this equation, a Crank-Nicholson (semi-implicit) diffusion scheme
    !    is used to solve the (1/rho_ds) * d [ rho_ds * K_zm * d(xm)/dz ] / dz
    !    eddy-diffusion term.  The discretized implicit portion of the term is
    !    written out as:
    !
    !    + (1/2) * invrs_rho_ds_zt(2)
    !      * dzt(2)
    !        * [ rho_ds_zm(2) * K_zm(2)
    !            * dzm(2) * ( xm(3,<t+1>) - xm(2,<t+1>) ) ]
    !    - invrs_rho_ds_zt(2)
    !      * dzt(2)
    !        * rho_ds_zm(1) 
    !        * [u_star^2/sqrt( um(2,<t>)^2 + vm(2,<t>)^2 )] * xm(2,<t+1>).
    !
    !    Note:  When the implicit term is brought over to the left-hand side,
    !           the signs are reversed and the leading "+" in front of the first
    !           part of the term is changed to a "-", while the leading "-" in
    !           front of the second part of the term is changed to a "+".
    !
    !    Note:  The x'w'|_sfc portion of the term written above has been pulled
    !           away from the rest of the implicit form written above because
    !           the (1/2) factor due to Crank-Nicholson time_stepping does not
    !           apply to it.  The x'w'|_sfc portion of the term is treated
    !           completely implicitly in order to enhance numerical stability.
    !
    !    The discretized explicit portion of the term is written out as:
    !
    !    + (1/2) * invrs_rho_ds_zt(2)
    !      * dzt(2)
    !        * [ rho_ds_zm(2) * K_zm(2)
    !            * dzm(2) * ( xm(3,<t>) - xm(2,<t>) ) ].
    !
    !    Timestep index (t) stands for the index of the current timestep, while
    !    timestep index (t+1) stands for the index of the next timestep, which
    !    is being advanced to in solving the d(xm)/dt equation.
    !
    !
    ! The lower boundary condition for the implicit and explicit portions of the
    ! turbulent advection term, without the x'w'|_sfc portion of the term, can
    ! easily be invoked by using the zero-flux boundary conditions found in the
    ! generalized diffusion function (function diffusion_zt_lhs), which is used
    ! for many other equations in this model.  Either the generalized explicit
    ! surface flux needs to be added onto the explicit term after the diffusion
    ! function has been called from subroutine windm_edsclrm_rhs, or the
    ! implicit momentum surface flux needs to be added onto the implicit term
    ! after the diffusion function has been called from subroutine
    ! windm_edsclrm_lhs.  However, all other equations in this model that use
    ! zero-flux diffusion have level 1 as the level to which the lower boundary
    ! condition needs to be applied.  Thus, an adjuster will have to be used at
    ! level 2 to call diffusion_zt_lhs with level 1 as the input level (the last
    ! variable being passed in during the function call).  However, the other
    ! variables passed in (rho_ds_zm*K_zm, gr%dzt, and gr%dzm variables) will
    ! have to be passed in as solving for level 2.
    !
    ! The value of xm(1) is located below the model surface and does not effect
    ! the rest of the model.  Since xm can be either a horizontal wind component
    ! or a generic eddy scalar quantity, the value of xm(1) is simply set to the
    ! value of xm(2) after the equation matrix has been solved.
    !
    !
    ! Conservation Properties:
    !
    ! When a fixed-flux lower boundary condition is used (combined with a
    ! zero-flux upper boundary condition), this technique of discretizing the
    ! turbulent advection term (treated as an eddy-diffusion term) leads to
    ! conservative differencing.  When the implicit momentum surface flux is
    ! either zero or not used, the column totals for each column in the
    ! left-hand side matrix (for the turbulent advection term) should be equal
    ! to 0.  Otherwise, the column total for the second column will be equal to
    ! rho_ds_zm(1) * x'w'|_sfc<t+1>.   When the generalized explicit surface
    ! flux is either zero or not used, the column total for the right-hand side
    ! vector (for the turbulent advection term) should be equal to 0.
    ! Otherwise, the column total for the right-hand side vector (for the
    ! turbulent advection term) will be equal to rho_ds_zm(1) * x'w'|_sfc<t>.
    ! This ensures that the total amount of quantity xm over the entire vertical
    ! domain is only changed by the surface flux (neglecting any forcing terms).
    ! The total amount of change is equal to rho_ds_zm(1) * x'w'|_sfc.
    !
    ! To see that this conservation law is satisfied by the left-hand side
    ! matrix, compute the turbulent advection (treated as eddy diffusion) of xm,
    ! neglecting any implicit momentum surface flux, multiply by rho_ds_zt, and
    ! integrate vertically.  In discretized matrix notation (where "i" stands
    ! for the matrix column and "j" stands for the matrix row):
    !
    !  0 = Sum_j Sum_i
    !       (rho_ds_zt)_i ( 1/dzt )_i
    !       ( 0.5_core_rknd * (1/rho_ds_zt) * dzt * (rho_ds_zm*K_zm*dzm) )_ij (xm<t+1>)_j.
    !
    ! The left-hand side matrix,
    ! ( 0.5_core_rknd * (1/rho_ds_zt) * dzt * (rho_ds_zm*K_zm*dzm) )_ij, is partially
    ! written below.  The sum over i in the above equation removes (1/rho_ds_zt)
    ! and dzt everywhere from the matrix below.  The sum over j leaves the
    ! column totals that are desired, which are 0.
    !
    ! Left-hand side matrix contributions from the turbulent advection term
    ! (treated as an eddy-diffusion term using a Crank-Nicholson timestep);
    ! first five vertical levels:
    !
    !     ------------------------------------------------------------------------------->
    !k=1 |  0             0                        0                          0
    !    |
    !k=2 |  0   +0.5*                  -0.5*                                  0
    !    |        (1/rho_ds_zt(k))*      (1/rho_ds_zt(k))*
    !    |        dzt(k)*                dzt(k)*
    !    |        rho_ds_zm(k)*          rho_ds_zm(k)*
    !    |        K_zm(k)*dzm(k)         K_zm(k)*dzm(k)
    !    |
    !k=3 |  0   -0.5*                  +0.5*                      -0.5*
    !    |        (1/rho_ds_zt(k))*      (1/rho_ds_zt(k))*          (1/rho_ds_zt(k))*
    !    |        dzt(k)*                dzt(k)*                    dzt(k)*
    !    |        rho_ds_zm(k-1)*        [ rho_ds_zm(k)*            rho_ds_zm(k)*
    !    |        K_zm(k-1)*dzm(k-1)       K_zm(k)*dzm(k)           K_zm(k)*dzm(k)
    !    |                                +rho_ds_zm(k-1)*
    !    |                                 K_zm(k-1)*dzm(k-1) ]
    !    |
    !k=4 |  0             0            -0.5*                      +0.5*
    !    |                               (1/rho_ds_zt(k))*          (1/rho_ds_zt(k))*
    !    |                               dzt(k)*                    dzt(k)*
    !    |                               rho_ds_zm(k-1)*            [ rho_ds_zm(k)*
    !    |                               K_zm(k-1)*dzm(k-1)           K_zm(k)*dzm(k)
    !    |                                                           +rho_ds_zm(k-1)*
    !    |                                                            K_zm(k-1)*dzm(k-1) ]
    !    |
    !k=5 |  0             0                        0              -0.5*
    !    |                                                          (1/rho_ds_zt(k))*
    !    |                                                          dzt(k)*
    !    |                                                          rho_ds_zm(k-1)*
    !    |                                                          K_zm(k-1)*dzm(k-1)
    !   \ /
    !
    ! Note:  The superdiagonal term from level 4 and both the main diagonal and
    !        superdiagonal terms from level 5 are not shown on this diagram.
    !
    ! Note:  If an implicit momentum surface flux is used, an additional term,
    !        + (1/rho_ds_zt(2)) * dzt(2) * rho_ds_zm(1)
    !          * [ u_star^2 / sqrt( um(2,<t>)^2 + vm(2,<t>)^2 ) ], is added to
    !        row 2 (k=2), column 2.
    !
    ! To see that the above conservation law is satisfied by the right-hand side
    ! vector, compute the turbulent advection (treated as eddy diffusion) of xm,
    ! neglecting any generalized explicit surface flux, multiply by rho_ds_zt,
    ! and integrate vertically.  In discretized matrix notation (where "i"
    ! stands for the matrix column and "j" stands for the matrix row):
    !
    !  0 = Sum_j Sum_i (rho_ds_zt)_i ( 1/dzt )_i ( rhs_vector )_j.
    !
    ! The right-hand side vector, ( rhs_vector )_j, is partially written below.
    ! The sum over i in the above equation removes (1/rho_ds_zt) and dzt
    ! everywhere from the vector below.  The sum over j leaves the column total
    ! that is desired, which is 0.
    !
    ! Right-hand side vector contributions from the turbulent advection term
    ! (treated as an eddy-diffusion term using a Crank-Nicholson timestep);
    ! first five vertical levels:
    !
    !     --------------------------------------------
    !k=1 |                      0                     |
    !    |                                            |
    !    |                                            |
    !k=2 | +0.5*(1/rho_ds_zt(k))*                     |
    !    |      dzt(k)*                               |
    !    |       [ rho_ds_zm(k)*K_zm(k)*              |
    !    |         dzm(k)*(xm(k+1,<t>)-xm(k,<t>)) ]   |
    !    |                                            |
    !k=3 | +0.5*(1/rho_ds_zt(k))*                     |
    !    |      dzt(k)*                               |
    !    |       [ rho_ds_zm(k)*K_zm(k)*              |
    !    |         dzm(k)*(xm(k+1,<t>)-xm(k,<t>))     |
    !    |        -rho_ds_zm(k-1)*K_zm(k-1)*          |
    !    |         dzm(k-1)*(xm(k,<t>)-xm(k-1,<t>)) ] |
    !    |                                            |
    !k=4 | +0.5*(1/rho_ds_zt(k))*                     |
    !    |      dzt(k)*                               |
    !    |       [ rho_ds_zm(k)*K_zm(k)*              |
    !    |         dzm(k)*(xm(k+1,<t>)-xm(k,<t>))     |
    !    |        -rho_ds_zm(k-1)*K_zm(k-1)*          |
    !    |         dzm(k-1)*(xm(k,<t>)-xm(k-1,<t>)) ] |
    !    |                                            |
    !k=5 | +0.5*(1/rho_ds_zt(k))*                     |
    !    |      dzt(k)*                               |
    !    |       [ rho_ds_zm(k)*K_zm(k)*              |
    !    |         dzm(k)*(xm(k+1,<t>)-xm(k,<t>))     |
    !    |        -rho_ds_zm(k-1)*K_zm(k-1)*          |
    !    |         dzm(k-1)*(xm(k,<t>)-xm(k-1,<t>)) ] |
    !   \ /                                          \ /
    !
    ! Note:  If a generalized explicit surface flux is used, an additional term,
    !        + (1/rho_ds_zt(2)) * dzt(2) * rho_ds_zm(1) * x'w'|_sfc, is added to
    !        row 2 (k=2).
    !
    ! Note:  Only the contributions by the turbulent advection term are shown
    !        for both the left-hand side matrix and the right-hand side vector.
    !        There are more terms in the equation, and thus more factors to be
    !        added to both the left-hand side matrix (such as time tendency and
    !        mean advection) and the right-hand side vector (such as xm
    !        forcings).  The left-hand side matrix is set-up so that a singular
    !        matrix is not encountered.

    ! References:
    ! Eqn. 8 & 9 on p. 3545 of
    ! ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    !   Method and Model Description'' Golaz, et al. (2002)
    ! JAS, Vol. 59, pp. 3540--3551.
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        grid ! Type

    use lapack_wrap, only:  & 
        tridag_solve, & ! Procedure(s)
        tridag_solvex

    use stats_variables, only: & 
        l_stats_samp

    use stats_type_utilities, only:  &
        stat_update_var_pt  ! Subroutine

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use stats_type, only: stats ! Type

    implicit none

    ! Constant parameters

    integer, parameter :: &
      kp1_tdiag = 1, & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2, & ! Thermodynamic main diagonal index.
      km1_tdiag = 3    ! Thermodynamic subdiagonal index.

    ! ------------------------ Input Variables ------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol
    
    type (grid), target, dimension(ngrdcol), intent(in) :: gr

    integer, intent(in) :: &
      nrhs ! Number of right-hand side (explicit) vectors & Number of solution vectors.

    integer, intent(in) :: &
      ixm_matrix_condt_num  ! Stats index of the condition numbers

    ! ------------------------ Inout variables ------------------------
    type (stats), target, dimension(ngrdcol), intent(inout) :: &
      stats_sfc
      
    real( kind = core_rknd ), dimension(3,ngrdcol,nz), intent(inout) :: &
      lhs    ! Implicit contributions to um, vm, and eddy scalars  [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nz,nrhs), intent(inout) :: &
      rhs    ! Right-hand side (explicit) contributions.

    ! ------------------------ Output variables ------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz,nrhs), intent(out) :: &
      solution ! Solution to the system of equations    [units vary]

    ! ------------------------ Local variables ------------------------
    real( kind = core_rknd ), dimension(ngrdcol) :: &
      rcond ! Estimate of the reciprocal of the condition number on the LHS matrix

    integer :: i
    
    ! ------------------------ Begin Code ------------------------

    ! Solve tridiagonal system for xm.
    if ( l_stats_samp .and. ixm_matrix_condt_num > 0 ) then
      
      do i = 1, ngrdcol
        call tridag_solvex( & 
               "windm_edsclrm", nz, nrhs, &                            ! intent(in) 
               lhs(kp1_tdiag,i,:), lhs(k_tdiag,i,:), lhs(km1_tdiag,i,:), rhs(i,:,:), & ! intent(inout)
               solution(i,:,:), rcond(i) )                                          ! intent(out)
      end do

      ! Est. of the condition number of the variance LHS matrix
      do i = 1, ngrdcol
        call stat_update_var_pt( ixm_matrix_condt_num, 1, 1.0_core_rknd/rcond(i), &  ! intent(in)
                                 stats_sfc(i) )                                      ! intent(inout)
      end do
    else

      do i = 1, ngrdcol
        call tridag_solve( "windm_edsclrm", nz, nrhs, &                             ! intent(in)
                           lhs(kp1_tdiag,i,:),  lhs(k_tdiag,i,:), lhs(km1_tdiag,i,:), rhs(i,:,:), & ! intent(inout)
                           solution(i,:,:) )                                                  ! intent(out)
      end do
    end if

    return
  end subroutine windm_edsclrm_solve

  !=============================================================================
  subroutine windm_edsclrm_implicit_stats( gr, solve_type, xm, & !intent(in)
                                           lhs_diff, lhs_ma_zt, & 
                                           invrs_rho_ds_zt, u_star_sqd,&
                                           rho_ds_zm, wind_speed, &
                                           l_imp_sfc_momentum_flux, &
                                           stats_zt ) ! intent(inout)

    ! Description:
    ! Compute implicit contributions to um and vm

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use stats_variables, only: & 
        ium_ma,  & ! Variables
        ium_ta,  & 
        ivm_ma,  &
        ivm_ta
        
    use constants_clubb, only: &
      zero

    use stats_type_utilities, only:  &
        stat_end_update_pt,  & ! Subroutines
        stat_update_var_pt

    use clubb_precision, only:  & 
        core_rknd

    use grid_class, only: &
        grid ! Type

    use stats_type, only: stats ! Type

    implicit none

    type (stats), target, intent(inout) :: &
      stats_zt

    type (grid), target, intent(in) :: gr

    ! Input variables
    integer, intent(in) :: & 
      solve_type     ! Desc. of what is being solved for

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      xm !  Computed value um or vm at <t+1>    [m/s]
      
    real( kind = core_rknd ), dimension(3,gr%nz), intent(in) :: &
      lhs_diff, & ! LHS diffustion terms
      lhs_ma_zt   ! LHS mean advection terms
      
    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  &
      wind_speed,      & ! wind speed; sqrt(u^2 + v^2)              [m/s]
      rho_ds_zm,       & ! Dry, static density on momentum levels      [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density at thermo. levels  [m^3/kg]
      
    real( kind = core_rknd ), intent(in) :: &
      u_star_sqd   ! Surface friction velocity, u_star, squared      [m/s]
      
    logical, intent(in) :: &
      l_imp_sfc_momentum_flux  ! Flag for implicit momentum surface fluxes.

    ! Local variables
    integer :: k, kp1, km1 ! Array indices
    
    real( kind = core_rknd ), dimension(gr%nz) :: &
      imp_sfc_flux

    ! Budget indices
    integer :: ixm_ma, ixm_ta

    select case ( solve_type )
    case ( windm_edsclrm_um )
      ixm_ma = ium_ma
      ixm_ta = ium_ta

    case ( windm_edsclrm_vm )
      ixm_ma = ivm_ma
      ixm_ta = ivm_ta

    case default
      ixm_ma = 0
      ixm_ta = 0

    end select
    
    imp_sfc_flux(:) = zero
    
    if ( l_imp_sfc_momentum_flux ) then

      ! Statistics:  implicit contributions for um or vm.

      ! xm term ta is modified at level 2 to include the effects of the
      ! surface flux.  In this case, this effects the implicit portion of
      ! the term, which handles the main diagonal for the turbulent advection term 
      if ( ium_ta + ivm_ta > 0 ) then
        imp_sfc_flux(2) =  - invrs_rho_ds_zt(2) * gr%invrs_dzt(2) &
                             * rho_ds_zm(1) * ( u_star_sqd / wind_speed(2) )
      endif

    endif ! l_imp_sfc_momentum_flux

    ! Finalize implicit contributions for xm

    do k = 2, gr%nz-1, 1

      km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nz )

      ! xm mean advection
      ! xm term ma is completely implicit; call stat_update_var_pt.
      call stat_update_var_pt( ixm_ma, k, & ! intent(in)
             - lhs_ma_zt(3,k) * xm(km1)   &
             - lhs_ma_zt(2,k) * xm(k)     &
             - lhs_ma_zt(1,k) * xm(kp1),  & ! intent(in)
              stats_zt )                    ! intent(inout)

      ! xm turbulent transport (implicit component)
      ! xm term ta has both implicit and explicit components;
      ! call stat_end_update_pt.
      call stat_end_update_pt( ixm_ta, k,               & ! intent(in)
             - 0.5_core_rknd * lhs_diff(3,k) * xm(km1)  &
             + ( - 0.5_core_rknd * lhs_diff(2,k)        &
                 + imp_sfc_flux(k) ) * xm(k)               &         
             - 0.5_core_rknd * lhs_diff(1,k) * xm(kp1), & ! intent(in) 
           stats_zt )                                     ! intent(inout)

    enddo


    ! Upper boundary conditions
    k   = gr%nz
    km1 = max( k-1, 1 )

    ! xm mean advection
    ! xm term ma is completely implicit; call stat_update_var_pt.
    call stat_update_var_pt( ixm_ma, k, & ! intent(in)
           - lhs_ma_zt(3,k) * xm(km1)   &
           - lhs_ma_zt(2,k) * xm(k),    & ! intent(in)
            stats_zt )                    ! intent(inout)

    ! xm turbulent transport (implicit component)
    ! xm term ta has both implicit and explicit components;
    ! call stat_end_update_pt.
    call stat_end_update_pt( ixm_ta, k,               & ! intent(in)
           - 0.5_core_rknd * lhs_diff(3,k) * xm(km1)  &
           + ( - 0.5_core_rknd * lhs_diff(2,k)        &
               + imp_sfc_flux(k) ) * xm(k),              & ! intent(in)
           stats_zt )                                   ! intent(inout)


    return
  end subroutine windm_edsclrm_implicit_stats

  !=============================================================================
  subroutine compute_uv_tndcy( nz, ngrdcol, gr, solve_type, &
                               fcor, perp_wind_m, perp_wind_g, &
                               xm_forcing, l_implemented, &
                               stats_zt, & 
                               xm_tndcy )

    ! Description:
    ! Computes the explicit tendency for the um and vm wind components.
    !
    ! The only explicit tendency that is involved in the d(um)/dt or d(vm)/dt
    ! equations is the Coriolis tendency.
    !
    ! The d(um)/dt equation contains the term:
    !
    ! - f * ( v_g - vm );
    !
    ! where f is the Coriolis parameter and v_g is the v component of the
    ! geostrophic wind.
    !
    ! Likewise, the d(vm)/dt equation contains the term:
    !
    ! + f * ( u_g - um );
    !
    ! where u_g is the u component of the geostrophic wind.
    !
    ! This term is treated completely explicitly.  The values of um, vm, u_g,
    ! and v_g are all found on the thermodynamic levels.
    !
    ! Wind forcing from the GCSS cases is also added here.
    !
    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        grid ! Type

    use stats_type_utilities, only: & 
        stat_update_var

    use stats_variables, only:      &
        ium_gf, & 
        ium_cf, & 
        ivm_gf, & 
        ivm_cf, & 
        ium_f,  &
        ivm_f,  &
        l_stats_samp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use stats_type, only: stats ! Type

    implicit none

    ! -------------------------- Input Variables --------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, dimension(ngrdcol), intent(in) :: gr
    
    integer, intent(in) ::  &
      solve_type      ! Description of what is being solved for

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  & 
      fcor            ! Coriolis parameter     [s^-1]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: & 
      perp_wind_m,  & ! Perpendicular component of the mean wind (e.g. v, for the u-eqn) [m/s]
      perp_wind_g,  & ! Perpendicular component of the geostropic wind (e.g. vg)         [m/s]
      xm_forcing      ! Prescribed wind forcing                                          [m/s/s]

    logical, intent(in) :: & 
      l_implemented   ! Flag for CLUBB being implemented in a larger model.
      
    ! -------------------------- Output Variables --------------------------
    type (stats), target, dimension(ngrdcol), intent(inout) :: &
      stats_zt

    ! -------------------------- Output Variables --------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) ::  &
      xm_tndcy        ! xm tendency            [m/s^2]

    ! -------------------------- Local Variables --------------------------
    integer :: & 
      ixm_gf, & 
      ixm_cf, &
      ixm_f

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: & 
      xm_gf, & 
      xm_cf
      
    integer :: i, k

    ! -------------------------- Begin Code --------------------------

    if ( .not. l_implemented ) then
      ! Only compute the Coriolis term if the model is running on it's own,
      ! and is not part of a larger, host model.

      select case ( solve_type )

      case ( windm_edsclrm_um )

        ixm_gf = ium_gf
        ixm_cf = ium_cf
        ixm_f  = ium_f
        
        do k = 1, nz
          do i = 1, ngrdcol
            xm_gf(i,k) = - fcor(i) * perp_wind_g(i,k)
          end do
        end do
          
        do k = 1, nz
          do i = 1, ngrdcol
            xm_cf(i,k) = fcor(i) * perp_wind_m(i,k)
          end do
        end do

      case ( windm_edsclrm_vm )

        ixm_gf = ivm_gf
        ixm_cf = ivm_cf
        ixm_f  = ivm_f

        do k = 1, nz
          do i = 1, ngrdcol
            xm_gf(i,k) = fcor(i) * perp_wind_g(i,k)
          end do
        end do

        do k = 1, nz
          do i = 1, ngrdcol
            xm_cf(i,k) = -fcor(i) * perp_wind_m(i,k)
          end do
        end do

      case default

        ixm_gf = 0
        ixm_cf = 0
        ixm_f = 0

        xm_gf = 0._core_rknd
        xm_cf = 0._core_rknd

      end select

      do k = 1, nz
        do i = 1, ngrdcol
          xm_tndcy(i,k) = xm_gf(i,k) + xm_cf(i,k) + xm_forcing(i,k)
        end do
      end do

      if ( l_stats_samp ) then
        do i = 1, ngrdcol
          ! xm term gf is completely explicit; call stat_update_var.
          call stat_update_var( ixm_gf, xm_gf(i,:), & ! intent(in)
                                stats_zt(i) )         ! intent(inout)

          ! xm term cf is completely explicit; call stat_update_var.
          call stat_update_var( ixm_cf, xm_cf(i,:), & ! intent(in)
                                stats_zt(i) )         ! intent(inout)

          ! xm term F
          call stat_update_var( ixm_f, xm_forcing(i,:), & ! intent(in)
                                stats_zt(i) )             ! intent(inout)
        end do
      endif

    else   ! implemented in a host model.

      xm_tndcy(:,:) = 0.0_core_rknd

    endif

    return
  end subroutine compute_uv_tndcy

!======================================================================================
  subroutine windm_edsclrm_lhs( nz, ngrdcol, gr, dt, &
                                lhs_ma_zt, lhs_diff, &
                                wind_speed, u_star_sqd,  &
                                rho_ds_zm, invrs_rho_ds_zt,  &
                                l_implemented, l_imp_sfc_momentum_flux,  &
                                lhs )
    ! Description:
    ! Calculate the implicit portion of the horizontal wind or eddy-scalar
    ! time-tendency equation.  See the description in subroutine
    ! windm_edsclrm_solve for more details.
    ! 
    ! Notes: 
    !   Lower Boundary:
    !       The lower boundary condition is a fixed-flux boundary condition, which
    !       gets added into the time-tendency equation at level 2.
    !       The value of xm(1) is located below the model surface and does not effect
    !       the rest of the model.  Since xm can be either a horizontal wind component
    !       or a generic eddy scalar quantity, the value of xm(1) is simply set to the
    !       value of xm(2) after the equation matrix has been solved.
    !
    !   --- THIS SUBROUTINE HAS BEEN OPTIMIZED ---
    !   Simple changes to this procedure may adversely affect computational speed
    !       - Gunther Huebler, Aug. 2018, clubb:ticket:834
    !----------------------------------------------------------------------------------

    use grid_class, only:  & 
        grid   ! Type

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    implicit none

    ! ----------------------- Input Variables -----------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, dimension(ngrdcol), intent(in) :: gr
    
    real( kind = core_rknd ), intent(in) :: & 
      dt                 ! Model timestep                             [s]
      
    real( kind = core_rknd ), intent(in), dimension(3,ngrdcol,nz) :: &
      lhs_diff, & ! LHS diffustion terms
      lhs_ma_zt   ! LHS mean advection terms

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      wind_speed,      & ! wind speed; sqrt( u^2 + v^2 )              [m/s]
      rho_ds_zm,       & ! Dry, static density on momentum levels     [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density at thermo. levels [m^3/kg]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      u_star_sqd    ! Surface friction velocity, u_*, squared  [m/s]

    logical, intent(in) ::  & 
      l_implemented, & ! Flag for CLUBB being implemented in a larger model.
      l_imp_sfc_momentum_flux  ! Flag for implicit momentum surface fluxes.

    ! ----------------------- Output Variable -----------------------
    real( kind = core_rknd ), dimension(3,ngrdcol,nz), intent(out) :: &
      lhs           ! Implicit contributions to xm (tridiagonal matrix)

    ! ----------------------- Local Variables -----------------------
    real( kind = core_rknd ) :: &
        invrs_dt    ! Inverse of dt, 1/dt, used for computational efficiency
        
    integer :: k, i ! Loop variable
    
    ! ----------------------- Begin Code -----------------------

    ! Calculate coefs of eddy diffusivity and inverse of dt
    invrs_dt = 1.0_core_rknd / dt

    ! Set lower boundary, see notes 
    do i = 1, ngrdcol
      lhs(1,i,1) = 0.0_core_rknd
      lhs(2,i,1) = 1.0_core_rknd
      lhs(3,i,1) = 0.0_core_rknd
    end do

    ! Add terms to lhs
    do k = 2, nz
      do i = 1, ngrdcol
        lhs(1,i,k) = 0.5_core_rknd * lhs_diff(1,i,k)

        lhs(2,i,k) = 0.5_core_rknd * lhs_diff(2,i,k)
        
        ! LHS time tendency.
        lhs(2,i,k) = lhs(2,i,k) + invrs_dt

        lhs(3,i,k) = 0.5_core_rknd * lhs_diff(3,i,k)
      end do
    end do 

    ! LHS mean advection term.
    if ( .not. l_implemented ) then

      do k = 2, nz-1
        do i = 1, ngrdcol
          lhs(1:3,i,k) = lhs(1:3,i,k) + lhs_ma_zt(:,i,k)
        end do
      end do

    endif

    if ( l_imp_sfc_momentum_flux ) then

      ! LHS momentum surface flux.
      do i = 1, ngrdcol
        lhs(2,i,2) = lhs(2,i,2) + invrs_rho_ds_zt(i,2) * gr(i)%invrs_dzt(2) &
                                  * rho_ds_zm(i,1) * ( u_star_sqd(i) / wind_speed(i,2) )
      end do
      
    end if ! l_imp_sfc_momentum_flux

    return
  end subroutine windm_edsclrm_lhs

  !=============================================================================
  subroutine windm_edsclrm_rhs( nz, ngrdcol, gr, solve_type, dt, &
                                lhs_diff, xm, xm_tndcy,  &
                                rho_ds_zm, invrs_rho_ds_zt,  &
                                l_imp_sfc_momentum_flux, xpwp_sfc, &
                                stats_zt, &
                                rhs )
    ! Description:
    !   Calculate the explicit portion of the horizontal wind or eddy-scalar
    !   time-tendency equation.  See the description in subroutine
    !   windm_edsclrm_solve for more details.
    ! 
    ! References:
    !   None
    ! 
    ! Notes:
    !   The lower boundary condition needs to be applied here at level 2.
    !   The lower boundary condition is a "fixed flux" boundary condition.
    !   The coding is the same as for a zero-flux boundary condition, but with
    !   an extra term added on the right-hand side at the boundary level.  For
    !   the rest of the model code, a zero-flux boundary condition is applied
    !   at level 1, and thus subroutine diffusion_zt_lhs is set-up to do that.
    !   In order to apply the same boundary condition code here at level 2, an
    !   adjuster needs to be used to tell diffusion_zt_lhs to use the code at
    !   level 2 that it normally uses at level 1.
    ! 
    !   --- THIS SUBROUTINE HAS BEEN OPTIMIZED ---
    !   Simple changes to this procedure may adversely affect computational speed
    !       - Gunther Huebler, Aug. 2018, clubb:ticket:834
    !----------------------------------------------------------------------------------

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use diffusion, only:  & 
        diffusion_zt_lhs    ! Procedure(s)

    use stats_variables, only: &
        ium_ta,  & ! Variable(s)
        ivm_ta,  &
        l_stats_samp

    use stats_type_utilities, only: &
        stat_begin_update_pt,  & ! Procedure(s)
        stat_modify_pt

    use grid_class, only:  & 
        grid ! Type

    use stats_type, only: stats ! Type

    implicit none

    !------------------- Input Variables -------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol
    
    type (grid), target, dimension(ngrdcol), intent(in) :: gr
    
    integer, intent(in) :: &
      solve_type ! Description of what is being solved for

    real( kind = core_rknd ), intent(in) :: & 
      dt                 ! Model timestep                             [s]
      
    real( kind = core_rknd ), dimension(3,ngrdcol,nz), intent(in) :: &
      lhs_diff   ! LHS diffustion terms

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      xm,              & ! Eddy-scalar variable, xm (thermo. levels)  [units vary]
      xm_tndcy,        & ! The explicit time-tendency acting on xm    [units vary]
      rho_ds_zm,       & ! Dry, static density on momentum levels     [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density at thermo. levels [m^3/kg]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      xpwp_sfc     ! x'w' at the surface                              [units vary]

    logical, intent(in) :: &
      l_imp_sfc_momentum_flux  ! Flag for implicit momentum surface fluxes.

    !------------------- Inout Variable -------------------
    type (stats), target, dimension(ngrdcol), intent(inout) :: &
      stats_zt

    !------------------- Output Variable -------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      rhs          ! Right-hand side (explicit) contributions.

    !------------------- Local Variables -------------------
    integer :: i, k    ! Loop variable

    real( kind = core_rknd ) :: invrs_dt

    integer :: ixm_ta

    !------------------- Begin Code -------------------

    invrs_dt = 1.0_core_rknd / dt   ! Precalculate 1.0/dt to avoid redoing the divide

    select case ( solve_type )
      case ( windm_edsclrm_um )
        ixm_ta = ium_ta
      case ( windm_edsclrm_vm )
        ixm_ta = ivm_ta
      case default  ! Eddy scalars
        ixm_ta = 0
    end select

    ! For purposes of the matrix equation, rhs(1) is simply set to 0.
    rhs(:,1) = 0.0_core_rknd

    ! Non-boundary rhs calculation, this is a highly vectorized loop
    do k = 2, nz-1
      do i = 1, ngrdcol
        rhs(i,k) = 0.5_core_rknd  & 
                   * ( - lhs_diff(3,i,k) * xm(i,k-1)      &
                       - lhs_diff(2,i,k) * xm(i,k)        &
                       - lhs_diff(1,i,k) * xm(i,k+1) )    &
                   + xm_tndcy(i,k)                        & ! RHS forcings
                   + invrs_dt * xm(i,k)                     ! RHS time tendency
      end do
    end do

    ! Upper boundary calculation
    do i = 1, ngrdcol
      rhs(i,nz) = 0.5_core_rknd  & 
                  * ( - lhs_diff(3,i,nz) * xm(i,nz-1)  &
                      - lhs_diff(2,i,nz) * xm(i,nz) )  &
                  + xm_tndcy(i,nz)                     & ! RHS forcings
                  + invrs_dt * xm(i,nz)                  ! RHS time tendency
    end do

    if ( l_stats_samp .and. ixm_ta > 0 ) then
      
      do i = 1, ngrdcol

        ! Statistics:  explicit contributions for um or vm.

        ! xm term ta has both implicit and explicit components; call
        ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
        ! subtracts the value sent in, reverse the sign on right-hand side
        ! turbulent advection component.
        do k = 2, nz-1
          
          call stat_begin_update_pt( ixm_ta, k, &                         ! intent(in)
                                     0.5_core_rknd  &
                                   * ( lhs_diff(3,i,k) * xm(i,k-1) &
                                   +   lhs_diff(2,i,k) * xm(i,k)   &      ! intent(in)
                                   +   lhs_diff(1,i,k) * xm(i,k+1) ), &
                                       stats_zt(i) )                      ! intent(inout)
        end do

        ! Upper boundary
        call stat_begin_update_pt( ixm_ta, nz, &
                                   0.5_core_rknd  &                       ! intent(in)
                                 * ( lhs_diff(3,i,nz) * xm(i,nz-1) &
                                 +   lhs_diff(2,i,nz) * xm(i,nz) ), &     ! intent(in)
                                     stats_zt(i) )                        ! intent(inout)
      end do
    endif

    if ( .not. l_imp_sfc_momentum_flux ) then

      ! RHS generalized surface flux.
      do i = 1, ngrdcol
        rhs(i,2) = rhs(i,2) + invrs_rho_ds_zt(i,2)  &
                              * gr(i)%invrs_dzt(2)  &
                              * rho_ds_zm(i,1) * xpwp_sfc(i)
      end do

      if ( l_stats_samp .and. ixm_ta > 0 ) then
        do i = 1, ngrdcol
          call stat_modify_pt( ixm_ta, 2,  &                    ! intent(in)  
                             + invrs_rho_ds_zt(i,2)  &  
                             * gr(i)%invrs_dzt(2)  &
                             * rho_ds_zm(i,1) * xpwp_sfc(i),  & ! intent(in)
                               stats_zt(i) )                    ! intent(inout)
        end do
      end if

    endif ! l_imp_sfc_momentum_flux

    return

  end subroutine windm_edsclrm_rhs

end module advance_windm_edsclrm_module
