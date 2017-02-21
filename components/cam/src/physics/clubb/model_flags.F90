!-----------------------------------------------------------------------
! $Id: model_flags.F90 8220 2016-07-21 18:48:32Z raut@uwm.edu $
!===============================================================================
module model_flags

! Description:
!   Various model options that can be toggled off and on as desired.

! References:
!   None
!-------------------------------------------------------------------------------

  implicit none

  public :: setup_model_flags, read_model_flags_from_file, setup_configurable_model_flags, &
    get_configurable_model_flags, write_model_flags_to_file

  private ! Default Scope

  logical, parameter, public ::  & 
    l_pos_def                     = .false., & ! Flux limiting positive definite scheme on rtm
    l_hole_fill                   = .true.,  & ! Hole filling pos def scheme on wp2,up2,rtp2,etc
    l_clip_semi_implicit          = .false., & ! Semi-implicit clipping scheme on wpthlp and wprtp
    l_clip_turb_adv               = .false., & ! Corrects thlm/rtm when w'th_l'/w'r_t' is clipped
    l_gmres                       = .false., & ! Use GMRES iterative solver rather than LAPACK
    l_sat_mixrat_lookup           = .false.    ! Use a lookup table for mixing length
                                      ! saturation vapor pressure calculations

  logical, parameter, public :: &
#ifdef BYTESWAP_IO
    l_byteswap_io   = .true.,  & ! Don't use the native byte ordering in GrADS output
#else
    l_byteswap_io   = .false., & ! Use the native byte ordering in GrADS output
#endif
    l_gamma_Skw     = .true.      ! Use a Skw dependent gamma parameter

  logical, parameter, public :: &
    l_use_boussinesq = .false.  ! Flag to use the Boussinesq form of the
                                ! predictive equations.  The predictive
                                ! equations are anelastic by default.

  logical, public :: &
    l_use_precip_frac = .true.   ! Flag to use precipitation fraction in KK
                                 ! microphysics.  The precipitation fraction
                                 ! is automatically set to 1 when this flag
                                 ! is turned off.

!$omp threadprivate( l_use_precip_frac )

  ! These flags determine whether or not we want CLUBB to do diffusion
  !   on thlm and rtm and if a stability correction is applied
  logical, public :: &
    l_diffuse_rtm_and_thlm        = .false., & ! Diffuses rtm and thlm
    l_stability_correct_Kh_N2_zm  = .false.    ! Divides Kh_N2_zm by a stability factor

!$omp threadprivate( l_diffuse_rtm_and_thlm, l_stability_correct_Kh_N2_zm )

  logical, parameter, public :: &
    l_morr_xp2_mc = .false. !Flag to include the effects of rain evaporation
                                  !on rtp2 and thlp2.  The moister (rt_1 or rt_2)
                                  !and colder (thl_1 or thl_2) will be fed into
                                  !the morrison microphys, and rain evaporation will
                                  !be allowed to increase variances

  logical, parameter, public :: &
    l_evaporate_cold_rcm = .false. ! Flag to evaporate cloud water at temperatures
                                   ! colder than -37C.  This is to be used for 
                                   ! Morrison microphysics, to prevent excess ice

  logical, parameter, public :: &
    l_cubic_interp = .false.      ! Flag to convert grid points with cubic monotonic
                                  ! spline interpolation as opposed to linear interpolation.

  ! See clubb:ticket:632 for details
  logical, public :: &
    l_calc_thlp2_rad = .true.         ! Include the contribution of radiation to thlp2
!$omp threadprivate( l_calc_thlp2_rad )

  ! These are the integer constants that represent the various saturation
  ! formulas. To add a new formula, add an additional constant here,
  ! add the logic to check the strings for the new formula in clubb_core and
  ! this module, and add logic in saturation to call the proper function--
  ! the control logic will be based on these named constants.

  integer, parameter, public :: &
    saturation_bolton = 1, & ! Constant for Bolton approximations of saturation
    saturation_gfdl   = 2, & ! Constant for the GFDL approximation of saturation
    saturation_flatau = 3    ! Constant for Flatau approximations of saturation

  !-----------------------------------------------------------------------------
  ! Options that can be changed at runtime 
  ! The default values are chosen below and overwritten if desired by the user
  !-----------------------------------------------------------------------------

  ! These flags determine whether we want to use an upwind differencing approximation 
  ! rather than a centered differencing for turbulent or mean advection terms.
  ! wpxp_ta affects wprtp, wpthlp, & wpsclrp
  ! xpyp_ta affects rtp2, thlp2, up2, vp2, sclrp2, rtpthlp, sclrprtp, & sclrpthlp
  ! xm_ma affects rtm, thlm, sclrm, um and vm.
  logical, public :: & 
    l_upwind_wpxp_ta = .false., & 
    l_upwind_xpyp_ta = .true.,  &
    l_upwind_xm_ma   = .true.

!$omp threadprivate(l_upwind_wpxp_ta, l_upwind_xpyp_ta, l_upwind_xm_ma)

  logical, public :: & 
    l_quintic_poly_interp = .false. ! Use a quintic polynomial in mono_cubic_interp

!$omp threadprivate(l_quintic_poly_interp)


  logical, public :: & 
    l_uv_nudge = .false.,  & ! For wind speed nudging. - Michael Falk
    l_rtm_nudge = .false., & ! For rtm nudging
    l_tke_aniso = .true.     ! For anisotropic turbulent kinetic energy, 
                             ! i.e. TKE = 1/2 (u'^2 + v'^2 + w'^2)
!$omp threadprivate(l_uv_nudge, l_tke_aniso, l_rtm_nudge)

  ! Use 2 calls to pdf_closure and the trapezoidal rule to  compute the 
  ! varibles that are output from high order closure
  logical, private :: &
    l_vert_avg_closure  = .false.
!$omp threadprivate(l_vert_avg_closure)

  ! These are currently set based on l_vert_avg_closure
  logical, public :: &
    l_trapezoidal_rule_zt = .false.,& ! If true, the trapezoidal rule is called for
                                         ! the thermodynamic-level variables output 
                                         ! from pdf_closure.  
    l_trapezoidal_rule_zm = .false.,& ! If true, the trapezoidal rule is called for
                                         ! three momentum-level variables - wpthvp,
                                       ! thlpthvp, and rtpthvp - output from pdf_closure.
    l_call_pdf_closure_twice = .false.,& ! This logical flag determines whether or not to
    ! call subroutine pdf_closure twice.  If true,
    ! pdf_closure is called first on thermodynamic levels
    ! and then on momentum levels so that each variable is
    ! computed on its native level.  If false, pdf_closure
    ! is only called on thermodynamic levels, and variables
    ! which belong on momentum levels are interpolated.
    l_single_C2_Skw = .false. ! Use a single Skewness dependent C2 for rtp2, thlp2, and rtpthlp

!$omp threadprivate(l_trapezoidal_rule_zt, l_trapezoidal_rule_zm, &
!$omp   l_call_pdf_closure_twice, l_single_C2_Skw)

  logical, public :: &
    l_standard_term_ta = .false.    ! Use the standard discretization for the
  ! turbulent advection terms.  Setting to
  ! .false. means that a_1 and a_3 are pulled
  ! outside of the derivative in advance_wp2_wp3_module.F90
  ! and in advance_xp2_xpyp_module.F90.
!$omp threadprivate(l_standard_term_ta)

  ! Use to determine whether a host model has already applied the surface flux,
  ! to avoid double counting.
  logical, public :: &
    l_host_applies_sfc_fluxes = .false.

!$omp threadprivate(l_host_applies_sfc_fluxes)

  ! Use cloud_cover and rcm_in_layer to help boost cloud_frac and rcm to help increase cloudiness
  ! at coarser grid resolutions.
  logical, public :: &
    l_use_cloud_cover = .false.
!$omp threadprivate(l_use_cloud_cover)

  integer, public :: &
    saturation_formula = saturation_flatau ! Integer that stores the saturation formula to be used

!$omp threadprivate(saturation_formula)

  ! See clubb:ticket:514 for details
  logical, public :: &
    l_diagnose_correlations = .false., & ! Diagnose correlations instead of using fixed ones
    l_calc_w_corr = .false.    ! Calculate the correlations between w and the hydrometeors
!$omp threadprivate(l_diagnose_correlations, l_calc_w_corr)

  logical, parameter, public :: &
    l_silhs_rad = .false.    ! Resolve radiation over subcolumns using SILHS

  logical, public :: &
    l_const_Nc_in_cloud = .false.,      & ! Use a constant cloud droplet conc. within cloud (K&K)
    l_fix_chi_eta_correlations = .true.   ! Use a fixed correlation for s and t Mellor(chi/eta) 
!$omp threadprivate( l_const_Nc_in_cloud, l_fix_chi_eta_correlations )

  logical, public :: &
    l_use_ADG2 = .false.    ! Use Luhar et al. (2002) to close the w Gaussians.
                            ! Allows for each w Gaussian to have a different
                            ! width
!$omp threadprivate(l_use_ADG2)

  logical, public :: &
    l_use_3D_closure = .false., & ! Use Luhar et al. (2002) to close the w, thl, and rt Gaussians.
    l_stability_correct_tau_zm = .true., & ! Use tau_N2_zm instead of tau_zm in wpxp_pr1
                                           !  stability correction
    l_damp_wp2_using_em = .false., & ! In wp2 equation, use a dissipation
                                     ! formula of -(2/3)*em/tau_zm, as in Bougeault (1981)
    l_do_expldiff_rtm_thlm = .false., & ! Diffuse rtm and thlm explicitly
    l_Lscale_plume_centered = .false., &   ! Alternate that uses the PDF to
                                           !  compute the perturbed values
    l_use_ice_latent = .false., &          ! Includes the effects of ice latent heating in
                                           !  turbulence terms
    l_use_C7_Richardson = .false., &       ! Parameterize C7 based on Richardson number
    l_use_C11_Richardson = .false., &      ! Parameterize C16 based on Richardson number
    l_brunt_vaisala_freq_moist = .false., &! Use a different formula for the Brunt-Vaisala 
                                           !  frequency in saturated atmospheres
                                           !  (from Durran and Klemp, 1982)
    l_use_thvm_in_bv_freq = .false., &     ! Use thvm in the calculation of Brunt-Vaisala frequency
    l_use_wp3_pr3 = .false., &             ! Include pressure term 3 (pr3) in wp3
    l_rcm_supersat_adj = .true.            ! Add excess supersaturated vapor to cloud water

!$omp threadprivate( l_use_3D_closure, l_stability_correct_tau_zm, l_damp_wp2_using_em, &
!$omp                l_do_expldiff_rtm_thlm, &
!$omp                l_Lscale_plume_centered, l_use_ice_latent, l_use_C7_Richardson, &
!$omp                l_use_C11_Richardson, l_brunt_vaisala_freq_moist, l_use_thvm_in_bv_freq, &
!$omp                l_use_wp3_pr3, l_rcm_supersat_adj )

#ifdef GFDL
  logical, public :: &
     I_sat_sphum       ! h1g, 2010-06-15
!$omp threadprivate( I_sat_sphum )
#endif

  namelist /configurable_clubb_flags_nl/ &
    l_upwind_wpxp_ta, l_upwind_xpyp_ta, l_upwind_xm_ma, l_quintic_poly_interp, &
    l_tke_aniso, l_vert_avg_closure, l_single_C2_Skw, l_standard_term_ta, &
    l_use_cloud_cover, l_calc_thlp2_rad, l_use_ADG2, l_use_3D_closure, l_rcm_supersat_adj

  contains

!===============================================================================
  subroutine setup_model_flags & 
             ( l_host_applies_sfc_fluxes_in, & 
               l_uv_nudge_in, saturation_formula_in &
#ifdef GFDL
               ,  I_sat_sphum_in   &  ! h1g, 2010-06-15
#endif
                )

! Description:
!   Setup flags that influence the numerics, etc. of CLUBB core

! References:
!   None
!-------------------------------------------------------------------------------

    implicit none

    ! External
    intrinsic :: trim

    ! Input Variables
    logical, intent(in) ::  & 
      l_host_applies_sfc_fluxes_in, &
      l_uv_nudge_in

    character(len=*), intent(in) :: &
      saturation_formula_in

#ifdef GFDL
    logical, intent(in) ::  & 
      I_sat_sphum_in           ! h1g, 2010-06-15
#endif

    !---- Begin Code ----

    ! Logicals

    l_uv_nudge  = l_uv_nudge_in

    l_host_applies_sfc_fluxes = l_host_applies_sfc_fluxes_in

    ! Integers

    ! Set up the saturation formula value
    select case ( trim( saturation_formula_in ) )
    case ( "bolton", "Bolton" )
      saturation_formula = saturation_bolton

    case ( "flatau", "Flatau" )
      saturation_formula = saturation_flatau

    case ( "gfdl", "GFDL" )
      saturation_formula = saturation_gfdl

      ! Add new saturation formulas after this.
    end select

#ifdef GFDL
    I_sat_sphum = I_sat_sphum_in  ! h1g, 2010-06-15
#endif
    return
  end subroutine setup_model_flags

!===============================================================================
  subroutine read_model_flags_from_file( iunit, filename )

! Description:
!   Read in some of the model flags of interest from a namelist file. If the
!   variable isn't in the file it will just be the default value.
!
! References:
!   None
!-------------------------------------------------------------------------------

    implicit none

    integer, intent(in) :: &
      iunit ! File I/O unit to use

    character(len=*), intent(in) :: &
      filename ! Name of the file with the namelist

   ! Read the namelist
    open(unit=iunit, file=filename, status='old', action='read')

    read(unit=iunit, nml=configurable_clubb_flags_nl)

    close(unit=iunit)

    if ( l_vert_avg_closure ) then
      l_trapezoidal_rule_zt    = .true.
      l_trapezoidal_rule_zm    = .true.
      l_call_pdf_closure_twice = .true.
    else
      l_trapezoidal_rule_zt    = .false.
      l_trapezoidal_rule_zm    = .false.
      l_call_pdf_closure_twice = .false.
    end if

    return
  end subroutine read_model_flags_from_file

!===============================================================================
  subroutine write_model_flags_to_file( iunit, filename )

! Description:
!   Write a new namelist for the configurable model flags
!
! References:
!   None
!-------------------------------------------------------------------------------

    implicit none

    integer, intent(in) :: &
      iunit ! File I/O unit to use

    character(len=*), intent(in) :: &
      filename ! Name of the file with the namelist

   ! Read the namelist
    open(unit=iunit, file=filename, status='unknown', action='write')

    write(unit=iunit, nml=configurable_clubb_flags_nl)

    close(unit=iunit)

    return
  end subroutine write_model_flags_to_file
!===============================================================================
  subroutine setup_configurable_model_flags &
             ( l_upwind_wpxp_ta_in, l_upwind_xpyp_ta_in, & 
               l_upwind_xm_ma_in, l_quintic_poly_interp_in, &
               l_vert_avg_closure_in, &
               l_single_C2_Skw_in, l_standard_term_ta_in, &
               l_tke_aniso_in, l_use_cloud_cover_in, l_use_ADG2_in, &
               l_use_3D_closure_in, l_rcm_supersat_adj_in )

! Description:
!   Set a model flag based on the input arguments for the purposes of trying
!   all possible combinations in the clubb_tuner.
!
! References:
!   None
!-------------------------------------------------------------------------------

    implicit none

    ! Input Variables
    logical, intent(in) :: &
      l_upwind_wpxp_ta_in, & ! Model flags
      l_upwind_xpyp_ta_in, & 
      l_upwind_xm_ma_in, &
      l_quintic_poly_interp_in, &
      l_vert_avg_closure_in, &
      l_single_C2_Skw_in, &
      l_standard_term_ta_in, &
      l_tke_aniso_in, &
      l_use_cloud_cover_in, &
      l_use_ADG2_in, &
      l_use_3D_closure_in, &
      l_rcm_supersat_adj_in
    ! ---- Begin Code ----

    l_upwind_wpxp_ta = l_upwind_wpxp_ta_in
    l_upwind_xpyp_ta = l_upwind_xpyp_ta_in
    l_upwind_xm_ma = l_upwind_xm_ma_in
    l_quintic_poly_interp = l_quintic_poly_interp_in
    l_vert_avg_closure = l_vert_avg_closure_in
    l_single_C2_Skw = l_single_C2_Skw_in
    l_standard_term_ta = l_standard_term_ta_in
    l_tke_aniso = l_tke_aniso_in
    l_use_cloud_cover = l_use_cloud_cover_in
    l_use_ADG2 = l_use_ADG2_in
    l_use_3D_closure = l_use_3D_closure_in
    l_rcm_supersat_adj = l_rcm_supersat_adj_in

    if ( l_vert_avg_closure ) then
      l_trapezoidal_rule_zt    = .true.
      l_trapezoidal_rule_zm    = .true.
      l_call_pdf_closure_twice = .true.
    else
      l_trapezoidal_rule_zt    = .false.
      l_trapezoidal_rule_zm    = .false.
      l_call_pdf_closure_twice = .false.
    end if

    return
  end subroutine setup_configurable_model_flags

!===============================================================================
  subroutine get_configurable_model_flags &
             ( l_upwind_wpxp_ta_out, l_upwind_xpyp_ta_out, & 
               l_upwind_xm_ma_out, l_quintic_poly_interp_out, &
               l_vert_avg_closure_out, &
               l_single_C2_Skw_out, l_standard_term_ta_out, &
               l_tke_aniso_out, l_use_cloud_cover_out, l_use_ADG2_out, &
               l_use_3D_closure_out, l_rcm_supersat_adj_out )

! Description:
!   Get the current model flags.
!
! References:
!   None
!-------------------------------------------------------------------------------

    implicit none

    ! Input Variables
    logical, intent(out) :: &
      l_upwind_wpxp_ta_out, & ! Model flags
      l_upwind_xpyp_ta_out, & 
      l_upwind_xm_ma_out, &
      l_quintic_poly_interp_out, &
      l_vert_avg_closure_out, &
      l_single_C2_Skw_out, &
      l_standard_term_ta_out, &
      l_tke_aniso_out, &
      l_use_cloud_cover_out, &
      l_use_ADG2_out, &
      l_use_3D_closure_out, &
      l_rcm_supersat_adj_out
    ! ---- Begin Code ----

    l_upwind_wpxp_ta_out = l_upwind_wpxp_ta
    l_upwind_xpyp_ta_out = l_upwind_xpyp_ta
    l_upwind_xm_ma_out = l_upwind_xm_ma
    l_quintic_poly_interp_out = l_quintic_poly_interp
    l_vert_avg_closure_out = l_vert_avg_closure
    l_single_C2_Skw_out = l_single_C2_Skw
    l_standard_term_ta_out = l_standard_term_ta
    l_tke_aniso_out = l_tke_aniso
    l_use_cloud_cover_out = l_use_cloud_cover
    l_use_ADG2_out = l_use_ADG2
    l_use_3D_closure_out =  l_use_3D_closure
    l_rcm_supersat_adj_out = l_rcm_supersat_adj
    return
  end subroutine get_configurable_model_flags

end module model_flags
