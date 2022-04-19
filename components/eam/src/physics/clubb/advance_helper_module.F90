!-------------------------------------------------------------------------
! $Id: advance_helper_module.F90 8769 2018-08-11 22:54:36Z vlarson@uwm.edu $
!===============================================================================
module advance_helper_module

! Description:
!   This module contains helper methods for the advance_* modules.
!------------------------------------------------------------------------

  implicit none

  public :: &
    set_boundary_conditions_lhs, &
    set_boundary_conditions_rhs, &
    calc_stability_correction,   &
    calc_brunt_vaisala_freq_sqd, &
    compute_Cx_fnc_Richardson, &
    term_wp2_splat, term_wp3_splat, &
    smooth_min, smooth_max, &
    smooth_heaviside_peskin

  private ! Set Default Scope

!===============================================================================
  interface smooth_min

    ! These functions wrap the intrinsic Fortran 'min' function in 
    ! zt2zm(zm2zt()), thus smoothing the result for a variable defined on the
    ! momentum grid ("zm").  These functions can accept two 1d arrays,
    ! or a scalar and a 1d array in either order. In the case of an 
    ! array being compared with a scalar (eg zero), an additional
    ! 'min' is applied to guarantee that the smoothing did not violate
    ! the original min requirement.

    module procedure smooth_min_sclr_array
    module procedure smooth_min_array_sclr
    module procedure smooth_min_arrays

  end interface

!===============================================================================
  interface smooth_max

    ! These functions wrap the intrinsic Fortran 'max' function in 
    ! zt2zm(zm2zt()), thus smoothing the result for a variable defined on the
    ! momentum grid ("zm").  These functions can accept two 1d arrays,
    ! or a scalar and a 1d array in either order. In the case of an  
    ! array being compared with a scalar (eg zero), an additional
    ! 'max' is applied to guarantee that the smoothing did not violate
    ! the original max requirement.

    module procedure smooth_max_sclr_array
    module procedure smooth_max_array_sclr
    module procedure smooth_max_arrays

  end interface

!===============================================================================
  contains

  !---------------------------------------------------------------------------
  subroutine set_boundary_conditions_lhs( diag_index, low_bound, high_bound, &
                                          lhs, &
                                          diag_index2, low_bound2, high_bound2 )

  ! Description:
  !   Sets the boundary conditions for a left-hand side LAPACK matrix.
  !
  ! References:
  !   none
  !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)
        
    use constants_clubb, only: &
        one, zero

    implicit none

    ! Exernal 
    intrinsic :: present

    ! Input Variables
    integer, intent(in) :: &
      diag_index, low_bound, high_bound ! boundary indexes for the first variable

    ! Input / Output Variables
    real( kind = core_rknd ), dimension(:,:), intent(inout) :: &
      lhs ! left hand side of the LAPACK matrix equation

    ! Optional Input Variables
    integer, intent(in), optional :: &
      diag_index2, low_bound2, high_bound2 ! boundary indexes for the second variable

    ! --------------------- BEGIN CODE ----------------------

    if ( ( present( low_bound2 ) .or. present( high_bound2 ) ) .and. &
         ( .not. present( diag_index2 ) ) ) then

      error stop "Boundary index provided without diag_index."

    end if

    ! Set the lower boundaries for the first variable
    lhs(:,low_bound) = zero
    lhs(diag_index,low_bound) = one

    ! Set the upper boundaries for the first variable
    lhs(:,high_bound) = zero
    lhs(diag_index,high_bound) = one

    ! Set the lower boundaries for the second variable, if it is provided
    if ( present( low_bound2 ) ) then

      lhs(:,low_bound2) = zero
      lhs(diag_index2,low_bound2) = one

    end if

    ! Set the upper boundaries for the second variable, if it is provided
    if ( present( high_bound2 ) ) then

      lhs(:,high_bound2) = zero
      lhs(diag_index2,high_bound2) = one

    end if

    return
  end subroutine set_boundary_conditions_lhs

  !--------------------------------------------------------------------------
  subroutine set_boundary_conditions_rhs( &
               low_value, low_bound, high_value, high_bound, &
               rhs, &
               low_value2, low_bound2, high_value2, high_bound2 )

  ! Description:
  !   Sets the boundary conditions for a right-hand side LAPACK vector.
  !
  ! References:
  !   none
  !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Exernal 
    intrinsic :: present

    ! Input Variables

    ! The values for the first variable
    real( kind = core_rknd ), intent(in) :: low_value, high_value

    ! The bounds for the first variable
    integer, intent(in) :: low_bound, high_bound

    ! Input / Output Variables

    ! The right-hand side vector
    real( kind = core_rknd ), dimension(:), intent(inout) :: rhs

    ! Optional Input Variables

    ! The values for the second variable
    real( kind = core_rknd ), intent(in), optional :: low_value2, high_value2

    ! The bounds for the second variable
    integer, intent(in), optional :: low_bound2, high_bound2


    ! -------------------- BEGIN CODE ------------------------

    ! Stop execution if a boundary was provided without a value
    if ( (present( low_bound2 ) .and. (.not. present( low_value2 ))) .or. &
         (present( high_bound2 ) .and. (.not. present( high_value2 ))) ) then

      error stop "Boundary condition provided without value."

    end if

    ! Set the lower and upper bounds for the first variable
    rhs(low_bound) = low_value
    rhs(high_bound) = high_value

    ! If a lower bound was given for the second variable, set it
    if ( present( low_bound2 ) ) then
      rhs(low_bound2) = low_value2
    end if

    ! If an upper bound was given for the second variable, set it
    if ( present( high_bound2 ) ) then
      rhs(high_bound2) = high_value2
    end if

    return
  end subroutine set_boundary_conditions_rhs

  !===============================================================================
  function calc_stability_correction( gr, thlm, Lscale, em, &
                                      exner, rtm, rcm, &
                                      p_in_Pa, thvm, ice_supersat_frac, &
                                      lambda0_stability_coef, &
                                      l_brunt_vaisala_freq_moist, &
                                      l_use_thvm_in_bv_freq ) &
    result ( stability_correction )
  !
  ! Description:
  !   Stability Factor
  !
  ! References:
  !
  !--------------------------------------------------------------------

    use constants_clubb, only: &
        zero, one, three    ! Constant(s)

    use grid_class, only:  &
        grid, & ! Type
        zt2zm    ! Procedure(s)

    use clubb_precision, only:  &
        core_rknd ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: &
      Lscale,          & ! Turbulent mixing length                   [m]
      em,              & ! Turbulent Kinetic Energy (TKE)            [m^2/s^2]
      thlm,            & ! th_l (thermo. levels)                     [K]
      exner,           & ! Exner function                            [-]
      rtm,             & ! total water mixing ratio, r_t             [kg/kg]
      rcm,             & ! cloud water mixing ratio, r_c             [kg/kg]
      p_in_Pa,         & ! Air pressure                              [Pa]
      thvm,            & ! Virtual potential temperature             [K]
      ice_supersat_frac

    real( kind = core_rknd ), intent(in) :: &
      lambda0_stability_coef    ! CLUBB tunable parameter lambda0_stability_coef

    logical, intent(in) :: &
      l_brunt_vaisala_freq_moist, & ! Use a different formula for the Brunt-Vaisala frequency in
                                    ! saturated atmospheres (from Durran and Klemp, 1982)
      l_use_thvm_in_bv_freq         ! Use thvm in the calculation of Brunt-Vaisala frequency

    ! Result
    real( kind = core_rknd ), dimension(gr%nz) :: &
      stability_correction

   real( kind = core_rknd ), dimension(gr%nz) :: &
      brunt_vaisala_freq_sqd, & !  []
      brunt_vaisala_freq_sqd_mixed, &
      brunt_vaisala_freq_sqd_dry, & !  []
      brunt_vaisala_freq_sqd_moist, &
      brunt_vaisala_freq_sqd_plus, &
      lambda0_stability
      
    ! Locals
    type (grid), target, dimension(1) :: gr_col
    
    ! Input Variables
    real( kind = core_rknd ), dimension(1,gr%nz) :: &
      Lscale_col,          & ! Turbulent mixing length                   [m]
      em_col,              & ! Turbulent Kinetic Energy (TKE)            [m^2/s^2]
      thlm_col,            & ! th_l (thermo. levels)                     [K]
      exner_col,           & ! Exner function                            [-]
      rtm_col,             & ! total water mixing ratio, r_t             [kg/kg]
      rcm_col,             & ! cloud water mixing ratio, r_c             [kg/kg]
      p_in_Pa_col,         & ! Air pressure                              [Pa]
      thvm_col,            & ! Virtual potential temperature             [K]
      ice_supersat_frac_col
      
    ! Result
    real( kind = core_rknd ), dimension(1,gr%nz) :: &
      stability_correction_col

    real( kind = core_rknd ), dimension(1,gr%nz) :: &
      brunt_vaisala_freq_sqd_col, & !  []
      brunt_vaisala_freq_sqd_mixed_col, &
      brunt_vaisala_freq_sqd_dry_col, & !  []
      brunt_vaisala_freq_sqd_moist_col, &
      brunt_vaisala_freq_sqd_plus_col, &
      lambda0_stability_col

    !------------ Begin Code --------------
    gr_col(1) = gr
    thlm_col(1,:) = thlm
    exner_col(1,:) = exner
    rtm_col(1,:) = rtm
    rcm_col(1,:) = rcm
    p_in_Pa_col (1,:) = p_in_Pa
    thvm_col(1,:) = thvm
    ice_supersat_frac_col(1,:) = ice_supersat_frac
    
    call calc_brunt_vaisala_freq_sqd( gr%nz, 1, gr_col, thlm_col, &          ! intent(in)
                                      exner_col, rtm_col, rcm_col, p_in_Pa_col, thvm_col, & ! intent(in)
                                      ice_supersat_frac_col, &              ! intent(in)
                                      l_brunt_vaisala_freq_moist, &     ! intent(in)
                                      l_use_thvm_in_bv_freq, &          ! intent(in)
                                      brunt_vaisala_freq_sqd_col, &         ! intent(out)
                                      brunt_vaisala_freq_sqd_mixed_col,&    ! intent(out)
                                      brunt_vaisala_freq_sqd_dry_col, &     ! intent(out)
                                      brunt_vaisala_freq_sqd_moist_col, &   ! intent(out)
                                      brunt_vaisala_freq_sqd_plus_col )     ! intent(out)
 
   brunt_vaisala_freq_sqd = brunt_vaisala_freq_sqd_col(1,:)
   brunt_vaisala_freq_sqd_mixed = brunt_vaisala_freq_sqd_mixed_col(1,:)
   brunt_vaisala_freq_sqd_dry = brunt_vaisala_freq_sqd_dry_col(1,:)
   brunt_vaisala_freq_sqd_moist = brunt_vaisala_freq_sqd_moist_col(1,:)
   brunt_vaisala_freq_sqd_plus = brunt_vaisala_freq_sqd_plus_col(1,:)

    lambda0_stability = merge( lambda0_stability_coef, zero, brunt_vaisala_freq_sqd > zero )

    stability_correction = one &
    + min( lambda0_stability * brunt_vaisala_freq_sqd * zt2zm(gr, Lscale)**2 / em, three )

    return
  end function calc_stability_correction

  !===============================================================================
  subroutine calc_brunt_vaisala_freq_sqd(  nz, ngrdcol, gr, thlm, &
                                           exner, rtm, rcm, p_in_Pa, thvm, &
                                           ice_supersat_frac, &
                                           l_brunt_vaisala_freq_moist, &
                                           l_use_thvm_in_bv_freq, &
                                           brunt_vaisala_freq_sqd, &
                                           brunt_vaisala_freq_sqd_mixed,&
                                           brunt_vaisala_freq_sqd_dry, &
                                           brunt_vaisala_freq_sqd_moist, &
                                           brunt_vaisala_freq_sqd_plus )

  ! Description:
  !   Calculate the Brunt-Vaisala frequency squared, N^2.

  ! References:
  !   ?
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Konstant

    use constants_clubb, only: &
        grav, & ! Constant(s)
        Lv, Cp, Rd, ep, &
        one

    use parameters_model, only: & 
        T0 ! Variable! 

    use grid_class, only: &
        grid, & ! Type
        ddzt,   &  ! Procedure(s)
        zt2zm

    use T_in_K_module, only: &
        thlm2T_in_K ! Procedure

    use saturation, only: &
        sat_mixrat_liq ! Procedure

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, dimension(ngrdcol), intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      thlm,    &  ! th_l (thermo. levels)              [K]
      exner,   &  ! Exner function                     [-]
      rtm,     &  ! total water mixing ratio, r_t      [kg/kg]
      rcm,     &  ! cloud water mixing ratio, r_c      [kg/kg]
      p_in_Pa, &  ! Air pressure                       [Pa]
      thvm,    &  ! Virtual potential temperature      [K]
      ice_supersat_frac

    logical, intent(in) :: &
      l_brunt_vaisala_freq_moist, & ! Use a different formula for the Brunt-Vaisala frequency in
                                    ! saturated atmospheres (from Durran and Klemp, 1982)
      l_use_thvm_in_bv_freq         ! Use thvm in the calculation of Brunt-Vaisala frequency

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      brunt_vaisala_freq_sqd, & ! Brunt-Vaisala frequency squared, N^2 [1/s^2]
      brunt_vaisala_freq_sqd_mixed, &
      brunt_vaisala_freq_sqd_dry,&
      brunt_vaisala_freq_sqd_moist, &
      brunt_vaisala_freq_sqd_plus

    ! Local Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      T_in_K, T_in_K_zm, rsat, rsat_zm, thm, thm_zm, ddzt_thlm, &
      ddzt_thm, ddzt_rsat, ddzt_rtm, thvm_zm, ddzt_thvm

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      stat_dry, stat_liq, ddzt_stat_liq, ddzt_stat_liq_zm, &
      stat_dry_virtual, stat_dry_virtual_zm,  ddzt_rtm_zm


    integer :: i, k

  !---------------------------------------------------------------------
    !----- Begin Code -----
    ddzt_thlm = ddzt( nz, ngrdcol, gr, thlm )
    thvm_zm = zt2zm( nz, ngrdcol, gr, thvm )
    ddzt_thvm = ddzt( nz, ngrdcol, gr, thvm )

    if ( .not. l_brunt_vaisala_freq_moist ) then

        ! Dry Brunt-Vaisala frequency
        if ( l_use_thvm_in_bv_freq ) then

          brunt_vaisala_freq_sqd(:,:) = ( grav / thvm_zm(:,:) ) * ddzt_thvm(:,:)

        else

          brunt_vaisala_freq_sqd(:,:) = ( grav / T0 ) * ddzt_thlm(:,:)

        end if

        T_in_K = thlm2T_in_K( thlm, exner, rcm )
        T_in_K_zm = zt2zm( nz, ngrdcol, gr, T_in_K )
        rsat = sat_mixrat_liq( p_in_Pa, T_in_K )
        rsat_zm = zt2zm( nz, ngrdcol, gr, rsat )
        ddzt_rsat = ddzt( nz, ngrdcol, gr, rsat )
        thm = thlm + Lv/(Cp*exner) * rcm
        thm_zm = zt2zm( nz, ngrdcol, gr, thm )
        ddzt_thm = ddzt( nz, ngrdcol, gr, thm )
        ddzt_rtm = ddzt( nz, ngrdcol, gr, rtm )

        do k = 1, nz
          do i = 1, ngrdcol
            stat_dry(i,k)  =  Cp * T_in_K(i,k) + grav * gr(i)%zt(k)
          end do
        end do
        
        stat_liq  =  stat_dry -Lv * rcm
        ddzt_stat_liq       = ddzt( nz, ngrdcol, gr, stat_liq )
        ddzt_stat_liq_zm    = zt2zm( nz, ngrdcol, gr, ddzt_stat_liq)
        stat_dry_virtual    = stat_dry + Cp * T_in_K *(0.608*(rtm-rcm)- rcm)
        stat_dry_virtual_zm = zt2zm( nz, ngrdcol, gr, stat_dry_virtual)
        ddzt_rtm_zm         = zt2zm( nz, ngrdcol, gr, ddzt_rtm )

        brunt_vaisala_freq_sqd_dry = ( grav / thm_zm)* ddzt_thm

        do k=1, nz
          do i = 1, ngrdcol
            brunt_vaisala_freq_sqd_plus(i,k) = grav/stat_dry_virtual(i,k) &
                      * ( ( ice_supersat_frac(i,k) * 0.5 + (1-ice_supersat_frac(i,k))) &
                          * ddzt_stat_liq_zm(i,k) &
                          + ( ice_supersat_frac(i,k) * Lv - (1-ice_supersat_frac(i,k)) *0.608*Cp) &
                            * ddzt_rtm_zm(i,k) )
          end do
        end do ! k=1, gr%nz

        do k=1, nz
          do i = 1, ngrdcol
            ! In-cloud Brunt-Vaisala frequency. This is Eq. (36) of Durran and
            ! Klemp (1982)
            brunt_vaisala_freq_sqd_moist(i,k) = &
              grav * ( ((one + Lv*rsat_zm(i,k) / (Rd*T_in_K_zm(i,k))) / &
              (one + ep*(Lv**2)*rsat_zm(i,k)/(Cp*Rd*T_in_K_zm(i,k)**2))) * &
              ( (one/thm_zm(i,k) * ddzt_thm(i,k)) + (Lv/(Cp*T_in_K_zm(i,k)))*ddzt_rsat(i,k)) - &
              ddzt_rtm(i,k) )
          end do
        end do ! k=1, gr%nz

        brunt_vaisala_freq_sqd_mixed =  &
               merge (brunt_vaisala_freq_sqd_moist,brunt_vaisala_freq_sqd_dry,&
               ice_supersat_frac > 0 )

    else ! l_brunt_vaisala_freq_moist

        T_in_K = thlm2T_in_K( thlm, exner, rcm )
        T_in_K_zm = zt2zm( nz, ngrdcol, gr, T_in_K )
        rsat = sat_mixrat_liq( p_in_Pa, T_in_K )
        rsat_zm = zt2zm( nz, ngrdcol, gr, rsat )
        ddzt_rsat = ddzt( nz, ngrdcol, gr, rsat )
        thm = thlm + Lv/(Cp*exner) * rcm
        thm_zm = zt2zm( nz, ngrdcol, gr, thm )
        ddzt_thm = ddzt( nz, ngrdcol, gr, thm )
        ddzt_rtm = ddzt( nz, ngrdcol, gr, rtm )

        do k=1, nz
          do i = 1, ngrdcol

            ! In-cloud Brunt-Vaisala frequency. This is Eq. (36) of Durran and
            ! Klemp (1982)
            brunt_vaisala_freq_sqd(i,k) = &
              grav * ( ((one + Lv*rsat_zm(i,k) / (Rd*T_in_K_zm(i,k))) / &
                (one + ep*(Lv**2)*rsat_zm(i,k)/(Cp*Rd*T_in_K_zm(i,k)**2))) * &
                ( (one/thm_zm(i,k) * ddzt_thm(i,k)) +(Lv/(Cp*T_in_K_zm(i,k)))*ddzt_rsat(i,k)) - &
                ddzt_rtm(i,k) )
                
          end do
        end do ! k=1, gr%nz


    end if ! .not. l_brunt_vaisala_freq_moist

    return

  end subroutine calc_brunt_vaisala_freq_sqd

!===============================================================================
  subroutine compute_Cx_fnc_Richardson( nz, ngrdcol, gr, &
                                        thlm, um, vm, em, Lscale, exner, rtm, &
                                        rcm, p_in_Pa, thvm, rho_ds_zm, &
                                        ice_supersat_frac, &
                                        clubb_params, &
                                        l_brunt_vaisala_freq_moist, &
                                        l_use_thvm_in_bv_freq, &
                                        l_use_shear_Richardson, &
                                        stats_zm, & 
                                        Cx_fnc_Richardson )

  ! Description:
  !   Compute Cx as a function of the Richardson number

  ! References:
  !   cam:ticket:59
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd  ! Konstant

    use grid_class, only: &
        grid, & ! Type
        ddzt, & ! Procedure(s)
        zt2zm

    use constants_clubb, only: &
        one, zero

    use interpolation, only: &
        linear_interp_factor ! Procedure

    use parameter_indices, only: &
        nparams,             & ! Variable(s)
        iCx_min,             &
        iCx_max,             &
        iRichardson_num_min, &
        iRichardson_num_max

    use stats_variables, only: &
        iRichardson_num, &    ! Variable(s)
        ishear_sqd, &
        l_stats_samp

    use stats_type_utilities, only: &
        stat_update_var      ! Procedure

    use stats_type, only: stats ! Type

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (stats), target, dimension(ngrdcol), intent(inout) :: &
      stats_zm

    type (grid), target, dimension(ngrdcol), intent(in) :: gr

    ! Constant Parameters
    real( kind = core_rknd ), parameter :: &
      Richardson_num_divisor_threshold = 1.0e-6_core_rknd, &
      Cx_fnc_Richardson_below_ground_value = one

    logical, parameter :: &
      l_Cx_fnc_Richardson_vert_avg = .false.,& ! Vertically average Cx_fnc_Richardson over a
                                        !  distance of Lscale
      l_Richardson_vert_avg = .false. , & ! Vertically average Richardson_num over a
                                         !  distance of Lscale
      l_use_shear_turb_freq_sqd = .false.! Use turb_freq_sqd and shear_sqd in denominator of
                                         !  Richardson_num

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      thlm,      & ! th_l (liquid water potential temperature)      [K]
      um,        & ! u mean wind component (thermodynamic levels)   [m/s]
      vm,        & ! v mean wind component (thermodynamic levels)   [m/s]
      em,        & ! Turbulent Kinetic Energy (TKE)                 [m^2/s^2]
      Lscale,    & ! Turbulent mixing length                        [m]
      exner,     & ! Exner function                                 [-]
      rtm,       & ! total water mixing ratio, r_t                  [kg/kg]
      rcm,       & ! cloud water mixing ratio, r_c                  [kg/kg]
      p_in_Pa,   & ! Air pressure                                   [Pa]
      thvm,      & ! Virtual potential temperature                  [K]
      rho_ds_zm, &  ! Dry static density on momentum levels          [kg/m^3]
      ice_supersat_frac  ! ice cloud fraction

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    logical, intent(in) :: &
      l_brunt_vaisala_freq_moist, & ! Use a different formula for the Brunt-Vaisala frequency in
                                    ! saturated atmospheres (from Durran and Klemp, 1982)
      l_use_thvm_in_bv_freq,      & ! Use thvm in the calculation of Brunt-Vaisala frequency
      l_use_shear_Richardson        ! Use shear in the calculation of Richardson number

    ! Output Variable
    real( kind = core_rknd), dimension(ngrdcol,nz), intent(out) :: &
      Cx_fnc_Richardson

    ! Local Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      brunt_vaisala_freq_sqd, &
      brunt_vaisala_freq_sqd_mixed,&
      brunt_vaisala_freq_sqd_dry, &
      brunt_vaisala_freq_sqd_moist, &
      brunt_vaisala_freq_sqd_plus, &
      Richardson_num, &
      Ri_zm, &
      dum_dz, dvm_dz, &
      shear_sqd, &
      turb_freq_sqd, &
      Lscale_zm

    real ( kind = core_rknd ) :: &
      invrs_min_max_diff, &
      invrs_num_div_thresh

    real( kind = core_rknd ) :: &
      Richardson_num_max, & ! CLUBB tunable parameter Richardson_num_max
      Richardson_num_min    ! CLUBB tunable parameter Richardson_num_min
      
    integer :: i

  !-----------------------------------------------------------------------

    !----- Begin Code -----
    call calc_brunt_vaisala_freq_sqd( nz, ngrdcol, gr, thlm, &          ! intent(in)
                                      exner, rtm, rcm, p_in_Pa, thvm, & ! intent(in)
                                      ice_supersat_frac, &              ! intent(in)
                                      l_brunt_vaisala_freq_moist, &     ! intent(in)
                                      l_use_thvm_in_bv_freq, &          ! intent(in)
                                      brunt_vaisala_freq_sqd, &         ! intent(out)
                                      brunt_vaisala_freq_sqd_mixed,&    ! intent(out)
                                      brunt_vaisala_freq_sqd_dry, &     ! intent(out)
                                      brunt_vaisala_freq_sqd_moist, &   ! intent(out)
                                      brunt_vaisala_freq_sqd_plus )     ! intent(out)

    Richardson_num_max = clubb_params(iRichardson_num_max)
    Richardson_num_min = clubb_params(iRichardson_num_min)

    invrs_min_max_diff = one / ( Richardson_num_max - Richardson_num_min )
    invrs_num_div_thresh = one / Richardson_num_divisor_threshold

    Lscale_zm = zt2zm( nz, ngrdcol, gr, Lscale )

    if ( l_use_shear_turb_freq_sqd ) then
      ! Calculate shear_sqd
      dum_dz = ddzt( nz, ngrdcol, gr, um )
      dvm_dz = ddzt( nz, ngrdcol, gr, vm )
      shear_sqd = dum_dz**2 + dvm_dz**2

      turb_freq_sqd = em / Lscale_zm**2
      Richardson_num = brunt_vaisala_freq_sqd / max( shear_sqd, turb_freq_sqd, &
                                                     Richardson_num_divisor_threshold )

      if ( l_stats_samp ) then
        do i = 1, ngrdcol
          call stat_update_var( ishear_sqd, shear_sqd(i,:), & ! intent(in)
                                stats_zm(i) )               ! intent(inout)
        end do
      end if
      
    else

      if ( l_use_shear_Richardson ) then
         Richardson_num = brunt_vaisala_freq_sqd_mixed * invrs_num_div_thresh
         Ri_zm &
         = max( 1.0e-7_core_rknd, brunt_vaisala_freq_sqd_mixed ) &
           / max( ( ddzt( nz, ngrdcol, gr, um)**2 &
                    + ddzt( nz, ngrdcol, gr, vm )**2 ), 1.0e-7_core_rknd )
      else
         Richardson_num = brunt_vaisala_freq_sqd * invrs_num_div_thresh
         Ri_zm = Richardson_num
      endif

    end if

    if ( l_Richardson_vert_avg ) then
      ! Clip below-min values of Richardson_num
      Richardson_num = max( Richardson_num, Richardson_num_min )

      do i = 1, ngrdcol
        Richardson_num(i,:) = Lscale_width_vert_avg( gr(i), Richardson_num(i,:), Lscale_zm(i,:), rho_ds_zm(i,:), &
                                                Richardson_num_max )
      end do
    end if

    ! Cx_fnc_Richardson is interpolated based on the value of Richardson_num
    ! The min function ensures that Cx does not exceed Cx_max, regardless of the
    !     value of Richardson_num_max.
    Cx_fnc_Richardson = linear_interp_factor( &
                    ( max(min(Richardson_num_max,Ri_zm),Richardson_num_min) &
                    - Richardson_num_min )  * invrs_min_max_diff, &
                                              clubb_params(iCx_max), clubb_params(iCx_min) )

    if ( l_Cx_fnc_Richardson_vert_avg ) then
      do i = 1, ngrdcol
        Cx_fnc_Richardson(i,:) = Lscale_width_vert_avg( gr(i), Cx_fnc_Richardson(i,:), Lscale_zm(i,:), rho_ds_zm(i,:), &
                                                   Cx_fnc_Richardson_below_ground_value )
      end do
    end if

    ! On some compilers, roundoff error can result in Cx_fnc_Richardson being
    ! slightly outside the range [0,1]. Thus, it is clipped here.
    Cx_fnc_Richardson = max( zero, min( one, Cx_fnc_Richardson ) )

    ! Stats sampling
    if ( l_stats_samp ) then
      do i = 1, ngrdcol
        call stat_update_var( iRichardson_num, Richardson_num(i,:), & ! intent(in)
                              stats_zm(i) )                      ! intent(inout)
      end do
    end if

  end subroutine compute_Cx_fnc_Richardson
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  function Lscale_width_vert_avg( gr, var_profile, Lscale_zm, rho_ds_zm, var_below_ground_value )&
  result (Lscale_width_vert_avg_output)

  ! Description:
  !   Averages a profile with a running mean of width Lscale_zm

  ! References:
  !   cam:ticket:59

    use clubb_precision, only: &
        core_rknd ! Precision

    use grid_class, only: &
        grid ! Type
        
    use constants_clubb, only: &
        zero

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      var_profile, &      ! Profile on momentum levels
      Lscale_zm, &        ! Lscale on momentum levels
      rho_ds_zm           ! Dry static energy on momentum levels!

    real( kind = core_rknd ), intent(in) :: &
      var_below_ground_value ! Value to use below ground

    ! Result Variable
    real( kind = core_rknd ), dimension(gr%nz) :: &
      Lscale_width_vert_avg_output ! Vertically averaged profile (on momentum levels)

    ! Local Variables
    integer :: &
        k, i,        & ! Loop variable
        k_avg_lower, &
        k_avg_upper

    real( kind = core_rknd ), dimension(gr%nz) :: &
      one_half_avg_width, &
      numer_terms, &
      denom_terms

    integer :: n_below_ground_levels

    real( kind = core_rknd ) :: & 
      numer_integral, & ! Integral in the numerator (see description)
      denom_integral    ! Integral in the denominator (see description)

  !----------------------------------------------------------------------

    !----- Begin Code -----

    one_half_avg_width = max( Lscale_zm, 500.0_core_rknd )

    ! Pre calculate numerator and denominator terms
    do k=1, gr%nz
        numer_terms(k) = rho_ds_zm(k) * gr%dzm(k) * var_profile(k)
        denom_terms(k) = rho_ds_zm(k) * gr%dzm(k)
    end do

    k_avg_upper = 2
    k_avg_lower = 1

    ! For every grid level
    do k=1, gr%nz

        !-----------------------------------------------------------------------
        ! Hunt down all vertical levels with one_half_avg_width(k) of gr%zm(k).
        ! 
        !     k_avg_upper and k_avg_lower are saved each loop iteration, this 
        !     improves computational efficiency since their values are likely
        !     within one or two grid levels of where they were last found to 
        !     be. This is because one_half_avg_width does not change drastically
        !     from one grid level to the next. Thus, less searching is required
        !     by allowing the search to start at a value that is close to the
        !     desired value and allowing each value to increment or decrement 
        !     as needed. 
        !-----------------------------------------------------------------------


        ! Determine if k_avg_upper needs to increment or decrement
        if ( gr%zm(k_avg_upper) - gr%zm(k) > one_half_avg_width(k) ) then

            ! k_avg_upper is too large, decrement it
            do while ( gr%zm(k_avg_upper) - gr%zm(k) > one_half_avg_width(k) )
                k_avg_upper = k_avg_upper - 1
            end do

        elseif ( k_avg_upper < gr%nz ) then

            ! k_avg_upper is too small, increment it
            do while ( gr%zm(k_avg_upper+1) - gr%zm(k) <= one_half_avg_width(k) )

                k_avg_upper = k_avg_upper + 1

                if ( k_avg_upper == gr%nz ) exit

            end do

        end if


        ! Determine if k_avg_lower needs to increment or decrement
        if ( gr%zm(k) - gr%zm(k_avg_lower) > one_half_avg_width(k) ) then

            ! k_avg_lower is too small, increment it
            do while ( gr%zm(k) - gr%zm(k_avg_lower) > one_half_avg_width(k) )

                k_avg_lower = k_avg_lower + 1

            end do

        elseif ( k_avg_lower > 1 ) then

            ! k_avg_lower is too large, decrement it
            do while ( gr%zm(k) - gr%zm(k_avg_lower-1) <= one_half_avg_width(k) )

                k_avg_lower = k_avg_lower - 1

                if ( k_avg_lower == 1 ) exit

            end do 

        end if


        ! Compute the number of levels below ground to include.
        if ( k_avg_lower > 1 ) then

            ! k=1, the lowest "real" level, is not included in the average, so no
            ! below-ground levels should be included.
            n_below_ground_levels = 0

            numer_integral = zero
            denom_integral = zero

        else

            ! The number of below-ground levels included is equal to the distance
            ! below the lowest level spanned by one_half_avg_width(k)
            ! divided by the distance between vertical levels below ground; the
            ! latter is assumed to be the same as the distance between the first and
            ! second vertical levels.
            n_below_ground_levels = int( ( one_half_avg_width(k)-(gr%zm(k)-gr%zm(1)) ) / &
                                        ( gr%zm(2)-gr%zm(1) ) )

            numer_integral = n_below_ground_levels * denom_terms(1) * var_below_ground_value
            denom_integral = n_below_ground_levels * denom_terms(1)

        end if

            
        ! Add numerator and denominator terms for all above-ground levels
        do i = k_avg_lower, k_avg_upper

            numer_integral = numer_integral + numer_terms(i)
            denom_integral = denom_integral + denom_terms(i)

        end do

        Lscale_width_vert_avg_output(k) = numer_integral / denom_integral

    end do

    return

  end function Lscale_width_vert_avg

 !============================================================================
  subroutine term_wp2_splat( nz, ngrdcol, gr, C_wp2_splat, dt, &
                             wp2, wp2_zt, tau_zm, &
                             wp2_splat )


  ! Description:
  !   This subroutine computes the (negative) tendency of wp2 due
  !   to "splatting" of eddies, e.g., near the ground or a Sc inversion.
  !   term_splat is intended to be added to the right-hand side of 
  !   the wp2 equation, and -0.5*term_splat is intended to be added to each 
  !   of the up2 and vp2 equations.  The functional form of term splat is
  !
  !   term_splat \propto - w'2 * (turbulent time scale) * ( d/dz( sqrt(w'2) ) )^2

    ! Included Modules
    use grid_class, only: &
        ddzt,  &   ! Procedure(s)
        grid

    use clubb_precision, only: & 
        core_rknd 

    use constants_clubb, only: &
        five  ! Constant(s)

    implicit none 
    
    ! Input Variables
    integer, intent(in) :: &
      nz, &
      ngrdcol
      
    type(grid), target, dimension(ngrdcol), intent(in) :: gr

    real( kind = core_rknd ), intent(in) :: & 
      C_wp2_splat, &          ! Tuning parameter                    [-]
      dt                      ! CLUBB computational time step       [s]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: & 
      wp2,     &  ! Variance of vertical velocity on the momentum grid [m^2/s^2]
      wp2_zt,  &  ! Variance of vertical velocity on the thermodynamic grid [m^2/s^2]
      tau_zm     ! Turbulent time scale on the momentum grid               [s]

    ! Output Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: & 
      wp2_splat     ! Tendency of <w'^2> due to splatting of eddies (on zm grid) [m^2/s^3]

    ! Local Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: & 
      d_sqrt_wp2_dz    ! d/dz( sqrt( w'2 ) )                 [1/s]

    integer :: i, k

    ! ---- Begin Code ----

    d_sqrt_wp2_dz(:,:) = ddzt( nz, ngrdcol, gr, sqrt( wp2_zt ) )
    
    ! The splatting term is clipped so that the incremental change doesn't exceed 5 times the
    !   value of wp2 itself.  This prevents spikes in wp2 from being propagated to up2 and vp2.
    !   However, it does introduce undesired dependence on the time step.
    !   Someday we may wish to treat this term using a semi-implicit discretization.
    do k = 1, nz
      do i = 1, ngrdcol
        wp2_splat(i,k) = - wp2(i,k) * min( five/dt, C_wp2_splat * tau_zm(i,k) &
                                                    * d_sqrt_wp2_dz(i,k)**2 )
      end do
    end do

  end subroutine term_wp2_splat

  !============================================================================
  subroutine term_wp3_splat( nz, ngrdcol, gr, C_wp2_splat, dt, &
                             wp2, wp3, tau_zt, &
                             wp3_splat )

  ! Description:
  !   This subroutine computes the damping of wp3 due
  !   to "splatting" of eddies, e.g., as they approach the ground or a Sc inversion.
  !   term_wp3_splat is intended to be added to the right-hand side of 
  !   the wp3 equation.  The functional form of wp3_splat is
  !
  !   wp3_splat \propto - w'3 * (turbulent time scale) * ( d/dz( sqrt(w'2) ) )^2
  !
  !   If the coefficient on wp3_splat is at least 1.5 times greater than the 
  !   coefficient on wp2_splat, then skewness will be damped, promoting
  !   more stratiform layers.

    ! Included Modules
    use grid_class, only: &
        ddzm, &   ! Procedure(s)
        grid    

    use clubb_precision, only: & 
        core_rknd 

    use constants_clubb, only: &
        three, &  ! Constant(s)
        five

    implicit none
   
    ! Input Variables
    integer, intent(in) :: &
      nz, &
      ngrdcol
      
    type(grid), target, dimension(ngrdcol), intent(in) :: gr

    real( kind = core_rknd ), intent(in) :: & 
      C_wp2_splat, &          ! Tuning parameter                    [-]
      dt                      ! CLUBB computational time step       [s]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: & 
      wp2,     &  ! Variance of vertical velocity on the momentum grid [m^2/s^2]
      wp3,     &  ! Third moment of vertical velocity on the momentum grid [m^3/s^3]
      tau_zt      ! Turbulent time scale on the thermal grid               [s]

    ! Output Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: & 
      wp3_splat     ! Tendency of <w'^3> due to splatting of eddies (on zt grid) [m^3/s^4]

    ! Local Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: & 
      d_sqrt_wp2_dz    ! d/dz( sqrt( w'2 ) )                 [1/s]
      
    integer :: i, k

    ! ---- Begin Code ----

    d_sqrt_wp2_dz(:,:) = ddzm( nz, ngrdcol, gr, sqrt( wp2 ) )
    
    ! The splatting term is clipped so that the incremental change doesn't exceed 5 times the
    ! value of wp3 itself. Someday we may wish to treat this term using a semi-implicit 
    ! discretization.
    do k = 1, nz
      do i = 1, ngrdcol
        wp3_splat(i,k) = - wp3(i,k) &
                           * min( five/dt, three * C_wp2_splat * tau_zt(i,k) &
                                           * d_sqrt_wp2_dz(i,k)**2 )
      end do
    end do

  end subroutine term_wp3_splat

!===============================================================================
  function smooth_min_sclr_array( nz, ngrdcol, input_var1, input_var2, smth_coef ) &
  result( output_var )

  ! Description:
  !   Computes a smoothed version of the min function using zt2zm/zm2zt, using
  !   one scalar and one 1d array as inputs.

  ! References:
  !   See clubb:ticket:894, updated version: 965
  !----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd                     ! Constant(s)
        
    use constants_clubb, only: &
        one_half

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

  ! Input Variables
    real ( kind = core_rknd ), intent(in) :: &
      input_var1, &       ! Units vary
      smth_coef          ! smoothing happens on interval [-smth_range, +smth_range]

    real ( kind = core_rknd ), dimension(ngrdcol, nz), intent(in) :: &
      input_var2          ! Units vary

  ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol, nz) :: &
      output_var          ! Units vary

  !----------------------------------------------------------------------

    output_var = one_half * ( (input_var1+input_var2) - &
                              sqrt((input_var1-input_var2)**2 + smth_coef**2) )

    return
  end function smooth_min_sclr_array

!===============================================================================
  function smooth_min_array_sclr( nz, ngrdcol, input_var1, input_var2, smth_coef ) &
  result( output_var )

  ! Description:
  !   Computes a smoothed version of the min function using zt2zm/zm2zt, using
  !   one scalar and one 1d array as inputs.

  ! References:
  !   See clubb:ticket:894, updated version: 965
  !----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd                     ! Constant(s)
        
    use constants_clubb, only: &
        one_half

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

  ! Input Variables
    real ( kind = core_rknd ), dimension(ngrdcol, nz), intent(in) :: &
      input_var1          ! Units vary

    real ( kind = core_rknd ), intent(in) :: &
      input_var2, &       ! Units vary
      smth_coef          ! smoothing happens on interval [-smth_range, +smth_range]

  ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol, nz) :: &
      output_var          ! Units vary

  !----------------------------------------------------------------------

    output_var = one_half * ( (input_var1+input_var2) - &
                              sqrt((input_var1-input_var2)**2 + smth_coef**2) )

    return
  end function smooth_min_array_sclr

!===============================================================================
  function smooth_min_arrays( nz, ngrdcol, input_var1, input_var2, smth_coef ) &
  result( output_var )

  ! Description:
  !   Computes a smoothed version of the min function using zt2zm/zm2zt, using
  !   two 1d arrays as inputs.

  ! References:
  !   See clubb:ticket:894, updated version: 965
  !----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd                     ! Constant(s)
        
    use constants_clubb, only: &
        one_half

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

  ! Input Variables
    real ( kind = core_rknd ), dimension(ngrdcol, nz), intent(in) :: &
      input_var1, &       ! Units vary
      input_var2          ! Units vary
      
    real ( kind = core_rknd ), intent(in) :: &
      smth_coef          ! smoothing happens on interval [-smth_range, +smth_range]

  ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol, nz) :: &
      output_var          ! Units vary

  !----------------------------------------------------------------------

    output_var = one_half * ( (input_var1+input_var2) - &
                              sqrt((input_var1-input_var2)**2 + smth_coef**2) )

    return
  end function smooth_min_arrays

!===============================================================================
  function smooth_max_sclr_array( nz, ngrdcol, input_var1, input_var2, smth_coef ) &
  result( output_var )

  ! Description:
  !   Computes a smoothed version of the max function using zt2zm/zm2zt, 
  !   using one scalar and one 1d array as inputs.

  ! References:
  !   See clubb:ticket:894, updated version: 965
  !----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd                     ! Constant(s)
        
    use constants_clubb, only: &
        one_half

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

  ! Input Variables
    real ( kind = core_rknd ), intent(in) :: &
      input_var1, &       ! Units vary
      smth_coef          ! smoothing happens on interval [-smth_range, +smth_range]

    real ( kind = core_rknd ), dimension(ngrdcol, nz), intent(in) :: &
      input_var2          ! Units vary

  ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol, nz) :: &
      output_var          ! Units vary

  !----------------------------------------------------------------------
  
    output_var = one_half * ( (input_var1+input_var2) + &
                              sqrt((input_var1-input_var2)**2 + smth_coef**2) )

    return
  end function smooth_max_sclr_array

!===============================================================================
  function smooth_max_array_sclr( nz, ngrdcol, input_var1, input_var2, smth_coef ) &
  result( output_var )

  ! Description:
  !   Computes a smoothed version of the max function using zt2zm/zm2zt, 
  !   using one scalar and one 1d array as inputs.

  ! References:
  !   See clubb:ticket:894, updated version: 965
  !----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd                     ! Constant(s)
        
    use constants_clubb, only: &
        one_half

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

  ! Input Variables
    real ( kind = core_rknd ), dimension(ngrdcol, nz), intent(in) :: &
      input_var1          ! Units vary

    real ( kind = core_rknd ), intent(in) :: &
      input_var2, &       ! Units vary
      smth_coef          ! smoothing happens on interval [-smth_range, +smth_range]

  ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol, nz) :: &
      output_var          ! Units vary

  !----------------------------------------------------------------------

    output_var = one_half * ( (input_var1+input_var2) + &
                              sqrt((input_var1-input_var2)**2 + smth_coef**2) )

    return
  end function smooth_max_array_sclr

!===============================================================================
  function smooth_max_arrays( nz, ngrdcol, input_var1, input_var2, smth_coef ) &
  result( output_var )

  ! Description:
  !   Computes a smoothed version of the max function using zt2zm/zm2zt, using
  !   two 1d arrays as inputs.

  ! References:
  !   See clubb:ticket:894, updated version: 965
  !----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd                     ! Constant(s)
        
    use constants_clubb, only: &
        one_half

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

  ! Input Variables
    real ( kind = core_rknd ), dimension(ngrdcol, nz), intent(in) :: &
      input_var1, &       ! Units vary
      input_var2          ! Units vary
      
    real( kind = core_rknd ), intent(in) :: &
      smth_coef          ! smoothing happens on interval [-smth_range, +smth_range]

  ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol, nz) :: &
      output_var          ! Units vary

  !----------------------------------------------------------------------

    output_var = one_half * ( (input_var1+input_var2) + &
                              sqrt((input_var1-input_var2)**2 + smth_coef**2) )

    return
  end function smooth_max_arrays
  
  elemental function smooth_heaviside_peskin( input, smth_range ) &
    result( smth_output )
    
  ! Description:
  !   Computes a smoothed heaviside function as in 
  !   [Lin, Lee et al., 2005, A level set characteristic Galerkin finite element 
  !   method for free surface flows], equation (2)
  
  ! References:
  !   See clubb:ticket:965
  !----------------------------------------------------------------------
  
    use clubb_precision, only: &
        core_rknd                     ! Constant(s)
        
    use constants_clubb, only: &
        pi, invrs_pi, one, one_half, zero

    implicit none
    
    ! Input Variables      
    real ( kind = core_rknd ), intent(in) :: &
      input, &    ! Units vary
      smth_range  ! Outside of [-smth_range, smth_range], smooth Heaviside = Heaviside
    
    ! Local Variables
    real ( kind = core_rknd ) :: &
      input_over_smth_range  ! input divided by smth_range

    ! Output Variables
    real( kind = core_rknd ) :: &
      smth_output    ! Units vary
      
  !----------------------------------------------------------------------
    if (input < -smth_range ) then 
      smth_output = zero
    elseif ( input > smth_range ) then
       smth_output = one
    else 
      ! Note that this case will only ever be reached if smth_range != 0,
      ! so this division is fine and should not cause any issues
      input_over_smth_range = input / smth_range
      smth_output = one_half &
                    * (one + input_over_smth_range &
                       + invrs_pi * sin(pi * input_over_smth_range))
    end if
    
    return
  end function smooth_heaviside_peskin

end module advance_helper_module
