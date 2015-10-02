! $Id: variables_diagnostic_module.F90 5623 2012-01-17 17:55:26Z connork@uwm.edu $
module variables_diagnostic_module

! Description:
!   This module contains definitions of all diagnostic
!   arrays used in the single column model, as well as subroutines
!   to allocate, deallocate and initialize them.

!   Note that while these are all same dimension, there is a
!   thermodynamic and momentum grid and they have different levels
!-----------------------------------------------------------------------

  use pdf_parameter_module, only: &
    pdf_parameter ! derived type

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  private ! Set default scope

  public :: setup_diagnostic_variables, & 
            cleanup_diagnostic_variables


  ! Diagnostic variables

  real( kind = core_rknd ), target, allocatable, dimension(:), public :: &
    sigma_sqd_w_zt, & ! PDF width parameter interpolated to t-levs.  [-]
    Skw_zm,         & ! Skewness of w on momentum levels             [-]
    Skw_zt,         & ! Skewness of w on thermodynamic levels        [-]
    ug,             & ! u geostrophic wind                           [m/s]
    vg,             & ! v geostrophic wind                           [m/s]
    um_ref,         & ! Initial u wind; Michael Falk                 [m/s]
    vm_ref,         & ! Initial v wind; Michael Falk                 [m/s]
    thlm_ref,       & ! Initial liquid water potential temperature   [K]
    rtm_ref,        & ! Initial total water mixing ratio             [kg/kg]
    thvm              ! Virtual potential temperature                [K]

!!! Important Note !!!
! Do not indent the omp comments, they need to be in the first 4 columns
!!! End Important Note !!!
!$omp threadprivate(sigma_sqd_w_zt, Skw_zm, Skw_zt, ug, vg, &
!$omp   um_ref, vm_ref, thlm_ref, rtm_ref, thvm )

  real( kind = core_rknd ), target, allocatable, dimension(:), public :: & 
    rsat ! Saturation mixing ratio  ! Brian

!$omp threadprivate(rsat)

  type(pdf_parameter), allocatable, dimension(:), target, public :: &
    pdf_params_zm ! pdf_params on momentum levels  [units vary]

!$omp threadprivate(pdf_params_zm)

  real( kind = core_rknd ), target, allocatable, dimension(:), public :: & 
    Frad,         & ! Radiative flux (momentum point)   [W/m^2]
    radht,        & ! SW + LW heating rate              [K/s]
    Frad_SW_up,   & ! SW radiative upwelling flux       [W/m^2]
    Frad_LW_up,   & ! LW radiative upwelling flux       [W/m^2]
    Frad_SW_down, & ! SW radiative downwelling flux     [W/m^2]
    Frad_LW_down ! LW radiative downwelling flux        [W/m^2]

!$omp threadprivate(Frad, radht, Frad_SW_up, Frad_SW_down, Frad_LW_up, Frad_LW_down)

! Second order moments
  real( kind = core_rknd ), target, allocatable, dimension(:), public :: & 
    thlprcp,  & ! thl'rc'              [K kg/kg]
    rtprcp,   & ! rt'rc'               [kg^2/kg^2]
    rcp2        ! rc'^2                [kg^2/kg^2]

!$omp threadprivate(thlprcp, rtprcp, rcp2)

! Third order moments
  real( kind = core_rknd ), target, allocatable, dimension(:), public :: & 
    wpthlp2,   & ! w'thl'^2    [m K^2/s]
    wp2thlp,   & ! w'^2 thl'   [m^2 K/s^2]
    wprtp2,    & ! w'rt'^2     [m kg^2/kg^2]
    wp2rtp,    & ! w'^2rt'     [m^2 kg/kg]
    wprtpthlp, & ! w'rt'thl'   [m kg K/kg s]
    wp2rcp,    & ! w'^2 rc'    [m^2 kg/kg s^2]
    wp3_zm       ! w'^3        [m^3/s^3]

!$omp threadprivate(wpthlp2, wp2thlp, wprtp2, wp2rtp, &
!$omp   wprtpthlp, wp2rcp, wp3_zm )

! Fourth order moments
  real( kind = core_rknd ), target, allocatable, dimension(:), public :: & 
    wp4 ! w'^4      [m^4/s^4]

!$omp threadprivate(wp4)

! Buoyancy related moments
  real( kind = core_rknd ), target, allocatable, dimension(:), public :: & 
    rtpthvp,  & ! rt'thv'     [K kg/kg]
    thlpthvp, & ! thl'thv'    [K^2]
    wpthvp,   & ! w'thv'      [K m/s]
    wp2thvp     ! w'^2thv'    [K m^2/s^2]

!$omp threadprivate(rtpthvp, thlpthvp, wpthvp, wp2thvp)

  real( kind = core_rknd ), target, allocatable, dimension(:), public :: &
    Kh_zt, & ! Eddy diffusivity coefficient on thermodynamic levels   [m^2/s]
    Kh_zm    ! Eddy diffusivity coefficient on momentum levels        [m^2/s]

!$omp threadprivate(Kh_zt, Kh_zm)

! Mixing lengths
  real( kind = core_rknd ), target, allocatable, dimension(:), public :: & 
    Lscale, Lscale_up, Lscale_down ! [m]

!$omp threadprivate(Lscale, Lscale_up, Lscale_down)

  real( kind = core_rknd ), target, allocatable, dimension(:), public :: & 
    em,     & ! Turbulent Kinetic Energy (TKE)                        [m^2/s^2]
    tau_zm, & ! Eddy dissipation time scale on momentum levels        [s]
    tau_zt    ! Eddy dissipation time scale on thermodynamic levels   [s]

!$omp threadprivate(em, tau_zm, tau_zt)

! hydrometeors variable array
  real( kind = core_rknd ), allocatable, dimension(:,:), public :: hydromet
!$omp threadprivate(hydromet)

  real( kind = core_rknd ), target, allocatable, dimension(:), public :: & 
    Ncnm     ! Cloud nuclei number concentration       [num/m^3]
!$omp threadprivate(Ncnm)


! Surface data
  real( kind = core_rknd ), public  :: ustar ! Average value of friction velocity [m/s]

  real( kind = core_rknd ), public :: soil_heat_flux    ! Soil Heat Flux [W/m^2]
!$omp threadprivate(ustar, soil_heat_flux)

! Passive scalar variables

  real( kind = core_rknd ), target, allocatable, dimension(:,:), public :: & 
    wpedsclrp   ! w'edsclr'
!$omp threadprivate(wpedsclrp)

  real( kind = core_rknd ), target, allocatable, dimension(:,:), public :: & 
    sclrpthvp,   & ! sclr'th_v'
    sclrprcp,    & ! sclr'rc'
    wp2sclrp,    & ! w'^2 sclr'
    wpsclrp2,    & ! w'sclr'^2
    wpsclrprtp,  & ! w'sclr'rt'
    wpsclrpthlp    ! w'sclr'thl'

!$omp threadprivate(sclrpthvp, sclrprcp, &
!$omp   wp2sclrp, wpsclrp2, wpsclrprtp, wpsclrpthlp )

! Interpolated variables for tuning
!
  real( kind = core_rknd ), target, allocatable, dimension(:), public :: & 
    wp2_zt,     & ! w'^2 on thermo. grid     [m^2/s^2]
    thlp2_zt,   & ! thl'^2 on thermo. grid   [K^2]
    wpthlp_zt,  & ! w'thl' on thermo. grid   [m K/s]
    wprtp_zt,   & ! w'rt' on thermo. grid    [m kg/(kg s)]
    rtp2_zt,    & ! rt'^2 on therm. grid     [(kg/kg)^2]
    rtpthlp_zt, & ! rt'thl' on thermo. grid  [kg K/kg]
    up2_zt,     & ! u'^2 on thermo. grid     [m^2/s^2]
    vp2_zt,     & ! v'^2 on thermo. grid     [m^2/s^2]
    upwp_zt,    & ! u'w' on thermo. grid     [m^2/s^2]
    vpwp_zt       ! v'w' on thermo. grid     [m^2/s^2]

!$omp threadprivate(wp2_zt, thlp2_zt, wpthlp_zt, wprtp_zt, &
!$omp   rtp2_zt, rtpthlp_zt, &
!$omp   up2_zt, vp2_zt, upwp_zt, vpwp_zt)


! Latin Hypercube arrays.  Vince Larson 22 May 2005
  real( kind = core_rknd ), target, allocatable, dimension(:), public :: & 
    lh_AKm,   & ! Kessler ac estimate                 [kg/kg/s]
    AKm,       & ! Exact Kessler ac                    [kg/kg/s]
    AKstd,     & ! St dev of exact Kessler ac          [kg/kg/s]
    AKstd_cld, & ! Stdev of exact w/in cloud ac        [kg/kg/s]
    lh_rcm_avg,   & ! Monte Carlo rcm estimate            [kg/kg]
    AKm_rcm,   & ! Kessler ac based on rcm             [kg/kg/s]
    AKm_rcc      ! Kessler ac based on rcm/cloud_frac  [kg/kg/s]

!$omp threadprivate(lh_AKm, AKm, AKstd, AKstd_cld, lh_rcm_avg, AKm_rcm, &
!$omp   AKm_rcc)

  ! Diagnostics from the pdf_closure subroutine
  real( kind = core_rknd ), target, allocatable, dimension(:), public :: & 
    sptp_mellor_1, sptp_mellor_2, &      ! Covariance of s and t[(kg/kg)^2] 
    tp2_mellor_1, tp2_mellor_2,   &      ! Variance of t [(kg/kg)^2]
    corr_st_mellor1, corr_st_mellor2 ! Correlation between s and t [-]
!$omp threadprivate(sptp_mellor_1, sptp_mellor_2, tp2_mellor_1, tp2_mellor_2, &
!$omp   corr_st_mellor1, corr_st_mellor2 )

  real( kind = core_rknd ), target, allocatable, dimension(:), public :: & 
    Skw_velocity, & ! Skewness velocity    [m/s]
    a3_coef,      & ! The a3 coefficient from CLUBB eqns                [-]
    a3_coef_zt      ! The a3 coefficient interpolated to the zt grid    [-]

!$omp threadprivate(Skw_velocity, a3_coef, a3_coef_zt)

  real( kind = core_rknd ), target, allocatable, dimension(:), public :: & 
    wp3_on_wp2,   &  ! w'^3 / w'^2 on the zm grid [m/s]
    wp3_on_wp2_zt    ! w'^3 / w'^2 on the zt grid [m/s]

!$omp threadprivate(wp3_on_wp2, wp3_on_wp2_zt)

  contains

!-----------------------------------------------------------------------
  subroutine setup_diagnostic_variables( nzmax )
! Description:
!   Allocates and initializes prognostic scalar and array variables
!   for the CLUBB model code
!-----------------------------------------------------------------------

    use constants_clubb, only:  & 
      em_min ! Variables

    use parameters_model, only: & 
      hydromet_dim, & ! Variables
      sclr_dim, &
      edsclr_dim

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: nzmax ! Nunber of grid levels [-]

    ! Local Variables
    integer :: i

!   --- Allocation ---

    ! Diagnostic variables

    allocate( sigma_sqd_w_zt(1:nzmax) ) ! PDF width parameter interp. to t-levs.
    allocate( Skw_zm(1:nzmax) )         ! Skewness of w on momentum levels
    allocate( Skw_zt(1:nzmax) )         ! Skewness of w on thermodynamic levels
    allocate( ug(1:nzmax) )             ! u geostrophic wind
    allocate( vg(1:nzmax) )             ! v geostrophic wind
    allocate( um_ref(1:nzmax) )         ! Reference u wind for nudging; Michael Falk, 17 Oct 2007
    allocate( vm_ref(1:nzmax) )         ! Reference v wind for nudging; Michael Falk, 17 Oct 2007
    allocate( thlm_ref(1:nzmax) )       ! Reference liquid water potential for nudging
    allocate( rtm_ref(1:nzmax) )        ! Reference total water mixing ratio for nudging
    allocate( thvm(1:nzmax) )           ! Virtual potential temperature

    allocate( rsat(1:nzmax) )       ! Saturation mixing ratio  ! Brian

    allocate( Frad(1:nzmax) )      ! radiative flux (momentum point)
    allocate( Frad_SW_up(1:nzmax) )
    allocate( Frad_LW_up(1:nzmax) )
    allocate( Frad_SW_down(1:nzmax) )
    allocate( Frad_LW_down(1:nzmax) )

    allocate( radht(1:nzmax) )     ! SW + LW heating rate

    ! pdf_params on momentum levels
    allocate( pdf_params_zm(1:nzmax) )

    ! Second order moments

    allocate( thlprcp(1:nzmax) )   ! thl'rc'
    allocate( rtprcp(1:nzmax) )    ! rt'rc'
    allocate( rcp2(1:nzmax) )      ! rc'^2

    ! Third order moments

    allocate( wpthlp2(1:nzmax) )   ! w'thl'^2
    allocate( wp2thlp(1:nzmax) )   ! w'^2thl'
    allocate( wprtp2(1:nzmax) )    ! w'rt'^2
    allocate( wp2rtp(1:nzmax) )    ! w'^2rt'
    allocate( wprtpthlp(1:nzmax) ) ! w'rt'thl'
    allocate( wp2rcp(1:nzmax) )    ! w'^2rc'

    allocate( wp3_zm(1:nzmax) )    ! w'^3

    ! Fourth order moments

    allocate( wp4(1:nzmax) )

    ! Buoyancy related moments

    allocate( rtpthvp(1:nzmax) )  ! rt'thv'
    allocate( thlpthvp(1:nzmax) ) ! thl'thv'
    allocate( wpthvp(1:nzmax) )   ! w'thv'
    allocate( wp2thvp(1:nzmax) )  ! w'^2thv'

    allocate( Kh_zt(1:nzmax) )  ! Eddy diffusivity coefficient: thermo. levels
    allocate( Kh_zm(1:nzmax) )  ! Eddy diffusivity coefficient: momentum levels

    allocate( em(1:nzmax) )
    allocate( Lscale(1:nzmax) )
    allocate( Lscale_up(1:nzmax) )
    allocate( Lscale_down(1:nzmax) )

    allocate( tau_zm(1:nzmax) ) ! Eddy dissipation time scale: momentum levels
    allocate( tau_zt(1:nzmax) ) ! Eddy dissipation time scale: thermo. levels


    ! Interpolated Variables
    allocate( wp2_zt(1:nzmax) )     ! w'^2 on thermo. grid
    allocate( thlp2_zt(1:nzmax) )   ! thl'^2 on thermo. grid
    allocate( wpthlp_zt(1:nzmax) )  ! w'thl' on thermo. grid
    allocate( wprtp_zt(1:nzmax) )   ! w'rt' on thermo. grid
    allocate( rtp2_zt(1:nzmax) )    ! rt'^2 on thermo. grid
    allocate( rtpthlp_zt(1:nzmax) ) ! rt'thl' on thermo. grid
    allocate( up2_zt(1:nzmax) )     ! u'^2 on thermo. grid
    allocate( vp2_zt(1:nzmax) )     ! v'^2 on thermo. grid
    allocate( upwp_zt(1:nzmax) )    ! u'w' on thermo. grid
    allocate( vpwp_zt(1:nzmax) )    ! v'w' on thermo. grid


    ! Microphysics Variables
    allocate( Ncnm(1:nzmax) )
    allocate( hydromet(1:nzmax,1:hydromet_dim) ) ! All hydrometeor fields

    ! Variables for Latin hypercube microphysics.  Vince Larson 22 May 2005
    allocate( lh_AKm(1:nzmax) )    ! Kessler ac estimate
    allocate( AKm(1:nzmax) )        ! Exact Kessler ac
    allocate( AKstd(1:nzmax) )      ! St dev of exact Kessler ac
    allocate( AKstd_cld(1:nzmax) )  ! St dev of exact w/in cloud Kessler ac
    allocate( lh_rcm_avg(1:nzmax) )      ! Monte Carlo rcm estimate
    allocate( AKm_rcm(1:nzmax) )      ! Kessler ac based on rcm
    allocate( AKm_rcc(1:nzmax) )      ! Kessler ac based on rcm/cloud_frac
    ! End of variables for Latin hypercube.

    ! High-order passive scalars
    allocate( sclrpthvp(1:nzmax, 1:sclr_dim) )
    allocate( sclrprcp(1:nzmax, 1:sclr_dim) )

    allocate( wp2sclrp(1:nzmax, 1:sclr_dim) )
    allocate( wpsclrp2(1:nzmax, 1:sclr_dim) )
    allocate( wpsclrprtp(1:nzmax, 1:sclr_dim) )
    allocate( wpsclrpthlp(1:nzmax, 1:sclr_dim) )

    ! Eddy Diff. Scalars
    allocate( wpedsclrp(1:nzmax, 1:edsclr_dim) )

    ! Diagnostics for s and t Mellor
    allocate( sptp_mellor_1(1:nzmax) )
    allocate( sptp_mellor_2(1:nzmax) )
    allocate( tp2_mellor_1(1:nzmax) )
    allocate( tp2_mellor_2(1:nzmax) )
    allocate( corr_st_mellor1(1:nzmax) )
    allocate( corr_st_mellor2(1:nzmax) )

    allocate( Skw_velocity(1:nzmax) )

    allocate( a3_coef(1:nzmax) )
    allocate( a3_coef_zt(1:nzmax) )

    allocate( wp3_on_wp2(1:nzmax) )
    allocate( wp3_on_wp2_zt(1:nzmax) )

    !   --- Initializaton ---

    ! Diagnostic variables

    sigma_sqd_w_zt = 0.0_core_rknd ! PDF width parameter interp. to t-levs.
    Skw_zm         = 0.0_core_rknd ! Skewness of w on momentum levels
    Skw_zt         = 0.0_core_rknd ! Skewness of w on thermodynamic levels
    ug             = 0.0_core_rknd ! u geostrophic wind
    vg             = 0.0_core_rknd ! v geostrophic wind
    um_ref         = 0.0_core_rknd
    vm_ref         = 0.0_core_rknd
    thlm_ref       = 0.0_core_rknd
    rtm_ref        = 0.0_core_rknd

    thvm = 0.0_core_rknd  ! Virtual potential temperature
    rsat  = 0.0_core_rknd  ! Saturation mixing ratio  ! Brian

    radht = 0.0_core_rknd ! Heating rate
    Frad  = 0.0_core_rknd ! Radiative flux
    Frad_SW_up = 0.0_core_rknd
    Frad_LW_up = 0.0_core_rknd
    Frad_SW_down = 0.0_core_rknd
    Frad_LW_down = 0.0_core_rknd


    ! pdf_params on momentum levels
    pdf_params_zm%w1          = 0.0_core_rknd
    pdf_params_zm%w2          = 0.0_core_rknd
    pdf_params_zm%varnce_w1   = 0.0_core_rknd
    pdf_params_zm%varnce_w2   = 0.0_core_rknd
    pdf_params_zm%rt1         = 0.0_core_rknd
    pdf_params_zm%rt2         = 0.0_core_rknd
    pdf_params_zm%varnce_rt1  = 0.0_core_rknd
    pdf_params_zm%varnce_rt2  = 0.0_core_rknd
    pdf_params_zm%thl1        = 0.0_core_rknd
    pdf_params_zm%thl2        = 0.0_core_rknd
    pdf_params_zm%varnce_thl1 = 0.0_core_rknd
    pdf_params_zm%varnce_thl2 = 0.0_core_rknd
    pdf_params_zm%mixt_frac   = 0.0_core_rknd
    pdf_params_zm%rc1         = 0.0_core_rknd
    pdf_params_zm%rc2         = 0.0_core_rknd
    pdf_params_zm%rsl1        = 0.0_core_rknd
    pdf_params_zm%rsl2        = 0.0_core_rknd
    pdf_params_zm%cloud_frac1 = 0.0_core_rknd
    pdf_params_zm%cloud_frac2 = 0.0_core_rknd
    pdf_params_zm%s1          = 0.0_core_rknd
    pdf_params_zm%s2          = 0.0_core_rknd
    pdf_params_zm%stdev_s1    = 0.0_core_rknd
    pdf_params_zm%stdev_s2    = 0.0_core_rknd
    pdf_params_zm%rrtthl      = 0.0_core_rknd
    pdf_params_zm%alpha_thl   = 0.0_core_rknd
    pdf_params_zm%alpha_rt    = 0.0_core_rknd
    pdf_params_zm%crt1        = 0.0_core_rknd
    pdf_params_zm%crt2        = 0.0_core_rknd
    pdf_params_zm%cthl1       = 0.0_core_rknd
    pdf_params_zm%cthl2       = 0.0_core_rknd


    ! Second order moments
    thlprcp = 0.0_core_rknd
    rtprcp  = 0.0_core_rknd
    rcp2    = 0.0_core_rknd

    ! Third order moments
    wpthlp2   = 0.0_core_rknd
    wp2thlp   = 0.0_core_rknd
    wprtp2    = 0.0_core_rknd
    wp2rtp    = 0.0_core_rknd
    wp2rcp    = 0.0_core_rknd
    wprtpthlp = 0.0_core_rknd

    wp3_zm    = 0.0_core_rknd

    ! Fourth order moments
    wp4 = 0.0_core_rknd

    ! Buoyancy related moments
    rtpthvp  = 0.0_core_rknd ! rt'thv'
    thlpthvp = 0.0_core_rknd ! thl'thv'
    wpthvp   = 0.0_core_rknd ! w'thv'
    wp2thvp  = 0.0_core_rknd ! w'^2thv'

    ! Eddy diffusivity
    Kh_zt = 0.0_core_rknd  ! Eddy diffusivity coefficient: thermo. levels
    Kh_zm = 0.0_core_rknd  ! Eddy diffusivity coefficient: momentum levels

    ! TKE
    em       = em_min

    ! Length scale
    Lscale   = 0.0_core_rknd
    Lscale_up      = 0.0_core_rknd
    Lscale_down    = 0.0_core_rknd

    ! Dissipation time
    tau_zm = 0.0_core_rknd ! Eddy dissipation time scale: momentum levels
    tau_zt = 0.0_core_rknd ! Eddy dissipation time scale: thermo. levels

    ! Hydrometer types
    Ncnm(1:nzmax) = 0.0_core_rknd ! Cloud nuclei number concentration (COAMPS)

    do i = 1, hydromet_dim, 1
      hydromet(1:nzmax,i) = 0.0_core_rknd
    end do


    ! Variables for Latin hypercube microphysics.  Vince Larson 22 May 2005
    lh_AKm   = 0.0_core_rknd  ! Kessler ac estimate
    AKm       = 0.0_core_rknd  ! Exact Kessler ac
    AKstd     = 0.0_core_rknd  ! St dev of exact Kessler ac
    AKstd_cld = 0.0_core_rknd  ! St dev of exact w/in cloud Kessler ac
    lh_rcm_avg   = 0.0_core_rknd  ! Monte Carlo rcm estimate
    AKm_rcm   = 0.0_core_rknd  ! Kessler ac based on rcm
    AKm_rcc   = 0.0_core_rknd  ! Kessler ac based on rcm/cloud_frac

    ! Passive scalars
    if ( sclr_dim > 0 ) then
      sclrpthvp(:,:)     = 0.0_core_rknd
      sclrprcp(:,:)      = 0.0_core_rknd

      wp2sclrp(:,:)      = 0.0_core_rknd
      wpsclrp2(:,:)      = 0.0_core_rknd
      wpsclrprtp(:,:)    = 0.0_core_rknd
      wpsclrpthlp(:,:)   = 0.0_core_rknd

    end if

    if ( edsclr_dim > 0 ) then
      wpedsclrp(:,:)     = 0.0_core_rknd
    end if

    sptp_mellor_1 = 0.0_core_rknd
    sptp_mellor_2 = 0.0_core_rknd
    tp2_mellor_1  = 0.0_core_rknd
    tp2_mellor_2  = 0.0_core_rknd

    corr_st_mellor1 = 0.0_core_rknd
    corr_st_mellor2 = 0.0_core_rknd

    Skw_velocity = 0.0_core_rknd

    a3_coef    = 0.0_core_rknd
    a3_coef_zt = 0.0_core_rknd

    wp3_on_wp2    = 0.0_core_rknd
    wp3_on_wp2_zt = 0.0_core_rknd

    return
  end subroutine setup_diagnostic_variables

!------------------------------------------------------------------------
  subroutine cleanup_diagnostic_variables( )

! Description:
!   Subroutine to deallocate variables defined in module global
!------------------------------------------------------------------------

    implicit none


    ! --- Deallocate ---

    deallocate( sigma_sqd_w_zt ) ! PDF width parameter interp. to t-levs.
    deallocate( Skw_zm )         ! Skewness of w on momentum levels
    deallocate( Skw_zt )         ! Skewness of w on thermodynamic levels
    deallocate( ug )             ! u geostrophic wind
    deallocate( vg )             ! v geostrophic wind
    deallocate( um_ref )         ! u initial
    deallocate( vm_ref )         ! v initial
    deallocate( thlm_ref )
    deallocate( rtm_ref )

    deallocate( thvm )      ! virtual potential temperature
    deallocate( rsat )      ! saturation mixing ratio  ! Brian

    deallocate( Frad )      ! radiative flux (momentum point)

    deallocate( Frad_SW_up ) ! upwelling shortwave radiative flux
    deallocate( Frad_LW_up ) ! upwelling longwave radiative flux
    deallocate( Frad_SW_down ) ! downwelling shortwave radiative flux
    deallocate( Frad_LW_down ) ! downwelling longwave radiative flux

    deallocate( radht )     ! SW + LW heating rate

    deallocate( pdf_params_zm )

    ! Second order moments

    deallocate( thlprcp )   ! thl'rc'
    deallocate( rtprcp )    ! rt'rc'
    deallocate( rcp2 )      ! rc'^2

    ! Third order moments

    deallocate( wpthlp2 )   ! w'thl'^2
    deallocate( wp2thlp )   ! w'^2thl'
    deallocate( wprtp2 )    ! w'rt'^2
    deallocate( wp2rtp )    ! w'^2rt'
    deallocate( wprtpthlp ) ! w'rt'thl'
    deallocate( wp2rcp )    ! w'^2rc'

    deallocate( wp3_zm )

    ! Fourth order moments

    deallocate( wp4 )

    ! Buoyancy related moments

    deallocate( rtpthvp )  ! rt'thv'
    deallocate( thlpthvp ) ! thl'thv'
    deallocate( wpthvp )   ! w'thv'
    deallocate( wp2thvp )  ! w'^2thv'

    deallocate( Kh_zt )  ! Eddy diffusivity coefficient: thermo. levels
    deallocate( Kh_zm )  ! Eddy diffusivity coefficient: momentum levels

    deallocate( em )
    deallocate( Lscale )
    deallocate( Lscale_up )
    deallocate( Lscale_down )
    deallocate( tau_zm ) ! Eddy dissipation time scale: momentum levels
    deallocate( tau_zt ) ! Eddy dissipation time scale: thermo. levels

    ! Cloud water variables

    deallocate( Ncnm )

    deallocate( hydromet )  ! Hydrometeor fields


    ! Interpolated variables for tuning
    deallocate( wp2_zt )     ! w'^2 on thermo. grid
    deallocate( thlp2_zt )   ! th_l'^2 on thermo. grid
    deallocate( wpthlp_zt )  ! w'th_l' on thermo. grid
    deallocate( wprtp_zt )   ! w'rt' on thermo. grid
    deallocate( rtp2_zt )    ! rt'^2 on thermo. grid
    deallocate( rtpthlp_zt ) ! rt'th_l' on thermo. grid
    deallocate( up2_zt )     ! u'^2 on thermo. grid
    deallocate( vp2_zt )     ! v'^2 on thermo. grid
    deallocate( upwp_zt )    ! u'w' on thermo. grid
    deallocate( vpwp_zt )    ! v'w' on thermo. grid

    ! Variables for Latin hypercube microphysics.  Vince Larson 22 May 2005
    deallocate( lh_AKm )   ! Kessler ac estimate
    deallocate( AKm )       ! Exact Kessler ac
    deallocate( AKstd )     ! St dev of exact Kessler ac
    deallocate( AKstd_cld ) ! St dev of exact w/in cloud Kessler ac
    deallocate( lh_rcm_avg )   ! Monte Carlo rcm estimate
    deallocate( AKm_rcm )   ! Kessler ac based on rcm
    deallocate( AKm_rcc )   ! Kessler ac based on rcm/cloud_frac

    ! Passive scalars
    deallocate( sclrpthvp )
    deallocate( sclrprcp )

    deallocate( wp2sclrp )
    deallocate( wpsclrp2 )
    deallocate( wpsclrprtp )
    deallocate( wpsclrpthlp )

    deallocate( wpedsclrp )

    ! Diagnostics for s and t Mellor
    deallocate( sptp_mellor_1 )
    deallocate( sptp_mellor_2 )
    deallocate( tp2_mellor_1 )
    deallocate( tp2_mellor_2 )
    deallocate( corr_st_mellor1 )
    deallocate( corr_st_mellor2 )

    deallocate( Skw_velocity )

    deallocate( a3_coef )
    deallocate( a3_coef_zt )

    deallocate( wp3_on_wp2 )
    deallocate( wp3_on_wp2_zt )

    return
  end subroutine cleanup_diagnostic_variables

end module variables_diagnostic_module
