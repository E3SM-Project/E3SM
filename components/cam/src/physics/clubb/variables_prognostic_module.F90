!-----------------------------------------------------------------------
! $Id: variables_prognostic_module.F90 8020 2016-03-21 22:28:37Z raut@uwm.edu $
!===============================================================================
module variables_prognostic_module

!       This module contains definitions of all prognostic
!       arrays used in the single column model, as well as subroutines
!       to allocate, deallocate and initialize them.

!       Note that while these are all same dimension, there is a
!       thermodynamic grid and a momentum grid, and the grids have
!       different points.
!-----------------------------------------------------------------------
  use pdf_parameter_module, only: &
    pdf_parameter ! Derived type

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  private ! Set Default Scoping

  public ::  & 
    setup_prognostic_variables,  & 
    cleanup_prognostic_variables

  ! Prognostic variables
! ---> h1g, 2010-06-16
#ifdef GFDL
  real( kind = core_rknd ), target, allocatable, dimension(:), public :: & 
    um,      & ! u wind                        [m/s]
    vm,      & ! v wind                        [m/s]
    upwp,    & ! vertical u momentum flux      [m^2/s^2]
    vpwp,    & ! vertical v momentum flux      [m^2/s^2]
    up2,     & ! u'^2                          [m^2/s^2]
    vp2,     & ! v'^2                          [m^2/s^2]
    thlm,    & ! liquid potential temperature  [K]
!---> h1g
    temp_clubb, & ! air temperature [K]
!<--- h1g
    rtm,     & ! total water mixing ratio      [kg/kg]
    wprtp,   & ! w'rt'                         [(kg/kg) m/s]
    wpthlp,  & ! w'thl'                        [m K/s]
    wprcp,   & ! w'rc'                         [(kg/kg) m/s]
    wp2,     & ! w'^2                          [m^2/s^2]
    wp3,     & ! w'^3                          [m^3/s^3]
    rtp2,    & ! rt'^2                         [(kg/kg)^2]
    thlp2,   & ! thl'^2                        [K^2]
    rtpthlp    ! rt'thl'                       [kg/kg K]
!$omp threadprivate( temp_clubb )
#else
  real( kind = core_rknd ), target, allocatable, dimension(:), public :: & 
    um,      & ! u wind                        [m/s]
    vm,      & ! v wind                        [m/s]
    upwp,    & ! vertical u momentum flux      [m^2/s^2]
    vpwp,    & ! vertical v momentum flux      [m^2/s^2]
    up2,     & ! u'^2                          [m^2/s^2]
    vp2,     & ! v'^2                          [m^2/s^2]
    thlm,    & ! liquid potential temperature  [K]
    rtm,     & ! total water mixing ratio      [kg/kg]
    wprtp,   & ! w'rt'                         [(kg/kg) m/s]
    wpthlp,  & ! w'thl'                        [m K/s]
    wprcp,   & ! w'rc'                         [(kg/kg) m/s]
    wp2,     & ! w'^2                          [m^2/s^2]
    wp3,     & ! w'^3                          [m^3/s^3]
    rtp2,    & ! rt'^2                         [(kg/kg)^2]
    thlp2,   & ! thl'^2                        [K^2]
    rtpthlp    ! rt'thl'                       [kg/kg K]
#endif
! <--- h1g, 2010-06-16

!$omp   threadprivate(um, vm, upwp, vpwp, up2, vp2)
!$omp   threadprivate(thlm, rtm, wprtp, wpthlp, wprcp)
!$omp   threadprivate(wp2, wp3, rtp2, thlp2, rtpthlp)

  real( kind = core_rknd ), target, allocatable, dimension(:), public :: & 
    p_in_Pa,         & ! Pressure (Pa) (thermodynamic levels)          [Pa]
    exner,           & ! Exner function = ( p / p0 ) ** kappa          [-]
    rho,             & ! Density (thermodynamic levels)                [kg/m^3]
    rho_zm,          & ! Density on momentum levels                    [kg/m^3]
    rho_ds_zm,       & ! Dry, static density (momentum levels)         [kg/m^3]
    rho_ds_zt,       & ! Dry, static density (thermodynamic levels)    [kg/m^3]
    invrs_rho_ds_zm, & ! Inverse dry, static density (momentum levs.)  [m^3/kg]
    invrs_rho_ds_zt, & ! Inverse dry, static density (thermo. levs.)   [m^3/kg]
    thv_ds_zm,       & ! Dry, base-state theta_v (momentum levels)     [K]
    thv_ds_zt,       & ! Dry, base-state theta_v (thermodynamic levs.) [K]
    thlm_forcing,    & ! thlm large-scale forcing                      [K/s]
    rtm_forcing,     & ! rtm large-scale forcing                       [kg/kg/s]
    um_forcing,      & ! u wind forcing                                [m/s/s] 
    vm_forcing,      & ! v wind forcing                                [m/s/s]
    wprtp_forcing,   & ! <w'r_t'> forcing (momentum levels)      [m*K/s^2]
    wpthlp_forcing,  & ! <w'th_l'> forcing (momentum levels)     [m*(kg/kg)/s^2]
    rtp2_forcing,    & ! <r_t'^2> forcing (momentum levels)      [(kg/kg)^2/s]
    thlp2_forcing,   & ! <th_l'^2> forcing (momentum levels)     [K^2/s]
    rtpthlp_forcing    ! <r_t'th_l'> forcing (momentum levels)   [K*(kg/kg)/s]

!$omp   threadprivate( p_in_Pa, exner, rho, rho_zm, rho_ds_zm, rho_ds_zt, &
!$omp     invrs_rho_ds_zm, invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &
!$omp     thlm_forcing, rtm_forcing, um_forcing, vm_forcing, wprtp_forcing, &
!$omp     wpthlp_forcing, rtp2_forcing, thlp2_forcing, rtpthlp_forcing )

  ! Imposed large scale w
  real( kind = core_rknd ), target, allocatable, dimension(:), public :: & 
    wm_zm, & ! w on momentum levels              [m/s]
    wm_zt    ! w on thermodynamic levels         [m/s]

!$omp   threadprivate(wm_zm, wm_zt)

  ! Cloud water variables
  real( kind = core_rknd ), target, allocatable, dimension(:), public :: & 
    rcm,                & ! Cloud water mixing ratio                 [kg/kg]
    cloud_frac,         & ! Cloud fraction                           [-]
    ice_supersat_frac,  & ! Ice cloud fraction                       [-]
    rcm_in_layer,       & ! Cloud water mixing ratio in cloud layer  [kg/kg]
    cloud_cover           ! Cloud cover                              [-]

!$omp   threadprivate(rcm, cloud_frac, ice_supersat_frac, rcm_in_layer, cloud_cover)

  ! Surface fluxes
  real( kind = core_rknd ), public ::  & 
    wpthlp_sfc,        & ! w'thl'      [m K/s]
    wprtp_sfc,         & ! w'rt'       [m kg/(kg s)]
    upwp_sfc, vpwp_sfc   ! u'w' & v'w' [m^2/s^2]

!$omp   threadprivate(wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc)

  ! Surface fluxes for passive scalars
  real( kind = core_rknd ), dimension(:), allocatable, public :: & 
    wpsclrp_sfc,     & ! w'sclr' at surface    [units m/s]
    wpedsclrp_sfc      ! w'edsclr' at surface  [units m/s]

!$omp   threadprivate(wpsclrp_sfc, wpedsclrp_sfc)

  ! More surface data
  real( kind = core_rknd ), public ::  & 
    T_sfc,  & ! surface temperature     [K]
    p_sfc,  & ! surface pressure        [Pa]
    sens_ht,    & ! sensible heat flux      [K m/s]
    latent_ht       ! latent heat flux        [m/s]

!$omp   threadprivate(T_sfc, p_sfc, sens_ht, latent_ht)

  ! Passive scalars
  real( kind = core_rknd ), target, allocatable, dimension(:,:), public :: & 
    sclrm,           & ! Mean passive scalars           [units vary]
    sclrp2,          & ! sclr'^2                        [units^2]
    sclrprtp,        & ! sclr'rt'                       [units kg/kg]
    sclrpthlp,       & ! sclr'th_l'                     [units K]
    sclrm_forcing,   & ! Scalars' forcing               [units/s]
    edsclrm,         & ! Mean eddy-diffusivity scalars  [units vary]
    edsclrm_forcing, & ! Eddy-diff. scalars forcing     [units/s]
    wpsclrp            ! w'sclr'                        [units vary m/s]

!---> h1g, 2010-06-16
#ifdef GFDL
  real( kind = core_rknd ), target, allocatable, dimension( : , : , : ), public :: & 
    RH_crit  ! critical relative humidity for droplet and ice nucleation
!$omp threadprivate( RH_crit )
#endif
!<--- h1g, 2010-06-16

!$omp   threadprivate(sclrm, sclrp2, sclrprtp, sclrpthlp, sclrm_forcing, &
!$omp     edsclrm, edsclrm_forcing, wpsclrp)

  ! PDF parameters
  real( kind = core_rknd ), target, allocatable, dimension(:), public :: &
    sigma_sqd_w    ! PDF width parameter (momentum levels)   [-]

!$omp threadprivate(sigma_sqd_w)

  type(pdf_parameter), target, allocatable, dimension(:), public :: &
    pdf_params, &
    pdf_params_frz !for use when l_use_ice_latent = .true.

!$omp threadprivate(pdf_params, pdf_params_frz)

  contains
!-----------------------------------------------------------------------
  subroutine setup_prognostic_variables( nz )

! Description:
!   Allocates and Initializes prognostic scalar and array variables
!   for the  CLUBB parameterization.  Variables contained within this module
!   will be arguments to the advance_clubb_core subroutine rather than brought
!   in through a use statement.

! References:
!   None
!-----------------------------------------------------------------------
    use constants_clubb, only:  & 
        rt_tol,    & ! Constant(s)
        thl_tol,   &
        w_tol_sqd, &
        zero

    use parameters_model, only: & 
        sclr_dim,  & ! Variable(s)
        edsclr_dim

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    integer, intent(in) :: nz ! Number of grid levels [-]

    integer :: i

!   --- Allocation ---

! Prognostic variables

    allocate( um(1:nz) )        ! u wind
    allocate( vm(1:nz) )        ! v wind

    allocate( upwp(1:nz) )      ! vertical u momentum flux
    allocate( vpwp(1:nz) )      ! vertical v momentum flux

    allocate( up2(1:nz) )
    allocate( vp2(1:nz) )

    allocate( thlm(1:nz) )      ! liquid potential temperature
!---> h1g, 2010-06-16
#ifdef GFDL
    allocate( temp_clubb(1:nz) )      ! air temperature
#endif
!<--- h1g, 2010-06-16

    allocate( rtm(1:nz) )       ! total water mixing ratio
    allocate( wprtp(1:nz) )     ! w'rt'
    allocate( wpthlp(1:nz) )    ! w'thl'
    allocate( wprcp(1:nz) )     ! w'rc'
    allocate( wp2(1:nz) )       ! w'^2
    allocate( wp3(1:nz) )       ! w'^3
    allocate( rtp2(1:nz) )      ! rt'^2
    allocate( thlp2(1:nz) )     ! thl'^2
    allocate( rtpthlp(1:nz) )   ! rt'thlp'

    allocate( p_in_Pa(1:nz) )         ! pressure (pascals)
    allocate( exner(1:nz) )           ! exner function
    allocate( rho(1:nz) )             ! density: t points
    allocate( rho_zm(1:nz) )          ! density: m points
    allocate( rho_ds_zm(1:nz) )       ! dry, static density: m-levs
    allocate( rho_ds_zt(1:nz) )       ! dry, static density: t-levs
    allocate( invrs_rho_ds_zm(1:nz) ) ! inv. dry, static density: m-levs
    allocate( invrs_rho_ds_zt(1:nz) ) ! inv. dry, static density: t-levs
    allocate( thv_ds_zm(1:nz) )       ! dry, base-state theta_v: m-levs
    allocate( thv_ds_zt(1:nz) )       ! dry, base-state theta_v: t-levs

    allocate( thlm_forcing(1:nz) )    ! thlm ls forcing
    allocate( rtm_forcing(1:nz) )     ! rtm ls forcing
    allocate( um_forcing(1:nz) )      ! u forcing
    allocate( vm_forcing(1:nz) )      ! v forcing
    allocate( wprtp_forcing(1:nz) )   ! <w'r_t'> forcing (microphysics)
    allocate( wpthlp_forcing(1:nz) )  ! <w'th_l'> forcing (microphysics)
    allocate( rtp2_forcing(1:nz) )    ! <r_t'^2> forcing (microphysics)
    allocate( thlp2_forcing(1:nz) )   ! <th_l'^2> forcing (microphysics)
    allocate( rtpthlp_forcing(1:nz) ) ! <r_t'th_l'> forcing (microphysics)

    ! Imposed large scale w

    allocate( wm_zm(1:nz) )       ! momentum levels
    allocate( wm_zt(1:nz) )       ! thermodynamic levels

    ! Cloud water variables

    allocate( rcm(1:nz) )
    allocate( cloud_frac(1:nz) )
    allocate( ice_supersat_frac(1:nz) )
    allocate( rcm_in_layer(1:nz) )
    allocate( cloud_cover(1:nz) )

    ! Passive scalar variables
    ! Note that sclr_dim can be 0
    allocate( wpsclrp_sfc(1:sclr_dim) )
    allocate( sclrm(1:nz, 1:sclr_dim) )
    allocate( sclrp2(1:nz, 1:sclr_dim) )
    allocate( sclrm_forcing(1:nz, 1:sclr_dim) )
    allocate( sclrprtp(1:nz, 1:sclr_dim) )
    allocate( sclrpthlp(1:nz, 1:sclr_dim) )

    allocate( wpedsclrp_sfc(1:edsclr_dim) )
    allocate( edsclrm_forcing(1:nz, 1:edsclr_dim) )

    allocate( edsclrm(1:nz, 1:edsclr_dim) )
    allocate( wpsclrp(1:nz, 1:sclr_dim) )

!---> h1g, 2010-06-16
#ifdef GFDL
    allocate( RH_crit(1:nz, 1:min(1,sclr_dim), 2) )
#endif
!<--- h1g, 2010-06-16

    allocate( sigma_sqd_w(1:nz) )    ! PDF width parameter (momentum levels)

    ! Variables for pdf closure scheme
    allocate( pdf_params(1:nz) )
    allocate( pdf_params_frz(1:nz) )

!--------- Set initial values for array variables ---------

    ! Prognostic variables

    um(1:nz)      = 0.0_core_rknd     ! u wind
    vm (1:nz)     = 0.0_core_rknd     ! v wind

    upwp(1:nz)    = 0.0_core_rknd     ! vertical u momentum flux
    vpwp(1:nz)    = 0.0_core_rknd     ! vertical v momentum flux

    up2(1:nz)     = w_tol_sqd ! u'^2
    vp2(1:nz)     = w_tol_sqd ! v'^2
    wp2(1:nz)     = w_tol_sqd ! w'^2

    thlm(1:nz)    = 0.0_core_rknd         ! liquid potential temperature
    rtm(1:nz)     = 0.0_core_rknd         ! total water mixing ratio
    wprtp(1:nz)   = 0.0_core_rknd         ! w'rt'
    wpthlp(1:nz)  = 0.0_core_rknd         ! w'thl'
    wprcp(1:nz)   = 0.0_core_rknd         ! w'rc'
    wp3(1:nz)     = 0.0_core_rknd         ! w'^3
    rtp2(1:nz)    = rt_tol**2    ! rt'^2
    thlp2(1:nz)   = thl_tol**2   ! thl'^2
    rtpthlp(1:nz) = 0.0_core_rknd         ! rt'thl'

    p_in_Pa(1:nz)= 0.0_core_rknd           ! pressure (Pa)
    exner(1:nz) = 0.0_core_rknd            ! exner
    rho(1:nz)  = 0.0_core_rknd             ! density on thermo. levels
    rho_zm(1:nz)  = 0.0_core_rknd          ! density on moment. levels
    rho_ds_zm(1:nz) = 0.0_core_rknd        ! dry, static density: m-levs
    rho_ds_zt(1:nz) = 0.0_core_rknd        ! dry, static density: t-levs
    invrs_rho_ds_zm(1:nz) = 0.0_core_rknd  ! inv. dry, static density: m-levs
    invrs_rho_ds_zt(1:nz) = 0.0_core_rknd  ! inv. dry, static density: t-levs
    thv_ds_zm(1:nz) = 0.0_core_rknd        ! dry, base-state theta_v: m-levs
    thv_ds_zt(1:nz) = 0.0_core_rknd        ! dry, base-state theta_v: t-levs

    thlm_forcing(1:nz)    = zero  ! thlm large-scale forcing
    rtm_forcing(1:nz)     = zero  ! rtm large-scale forcing
    um_forcing(1:nz)      = zero  ! u forcing
    vm_forcing(1:nz)      = zero  ! v forcing
    wprtp_forcing(1:nz)   = zero  ! <w'r_t'> forcing (microphysics)
    wpthlp_forcing(1:nz)  = zero  ! <w'th_l'> forcing (microphysics)
    rtp2_forcing(1:nz)    = zero  ! <r_t'^2> forcing (microphysics)
    thlp2_forcing(1:nz)   = zero  ! <th_l'^2> forcing (microphysics)
    rtpthlp_forcing(1:nz) = zero  ! <r_t'th_l'> forcing (microphysics)

    ! Imposed large scale w

    wm_zm(1:nz) = 0.0_core_rknd      ! Momentum levels
    wm_zt(1:nz) = 0.0_core_rknd      ! Thermodynamic levels

    ! Cloud water variables

    rcm(1:nz)               = 0.0_core_rknd
    cloud_frac(1:nz)        = 0.0_core_rknd
    ice_supersat_frac(1:nz) = 0.0_core_rknd
    rcm_in_layer(1:nz)      = 0.0_core_rknd
    cloud_cover(1:nz)       = 0.0_core_rknd

    sigma_sqd_w           = 0.0_core_rknd ! PDF width parameter (momentum levels)

    ! Variables for PDF closure scheme
    pdf_params(:)%w_1                 = zero
    pdf_params(:)%w_2                 = zero
    pdf_params(:)%varnce_w_1          = zero
    pdf_params(:)%varnce_w_2          = zero
    pdf_params(:)%rt_1                = zero
    pdf_params(:)%rt_2                = zero
    pdf_params(:)%varnce_rt_1         = zero
    pdf_params(:)%varnce_rt_2         = zero
    pdf_params(:)%thl_1               = zero
    pdf_params(:)%thl_2               = zero
    pdf_params(:)%varnce_thl_1        = zero
    pdf_params(:)%varnce_thl_2        = zero
    pdf_params(:)%rrtthl              = zero
    pdf_params(:)%alpha_thl           = zero
    pdf_params(:)%alpha_rt            = zero
    pdf_params(:)%crt_1               = zero
    pdf_params(:)%crt_2               = zero
    pdf_params(:)%cthl_1              = zero
    pdf_params(:)%cthl_2              = zero
    pdf_params(:)%chi_1               = zero
    pdf_params(:)%chi_2               = zero
    pdf_params(:)%stdev_chi_1         = zero
    pdf_params(:)%stdev_chi_2         = zero
    pdf_params(:)%stdev_eta_1         = zero
    pdf_params(:)%stdev_eta_2         = zero
    pdf_params(:)%covar_chi_eta_1     = zero
    pdf_params(:)%covar_chi_eta_2     = zero
    pdf_params(:)%corr_chi_eta_1      = zero
    pdf_params(:)%corr_chi_eta_2      = zero
    pdf_params(:)%rsatl_1             = zero
    pdf_params(:)%rsatl_2             = zero
    pdf_params(:)%rc_1                = zero
    pdf_params(:)%rc_2                = zero
    pdf_params(:)%cloud_frac_1        = zero
    pdf_params(:)%cloud_frac_2        = zero
    pdf_params(:)%mixt_frac           = zero
    pdf_params(:)%ice_supersat_frac_1 = zero
    pdf_params(:)%ice_supersat_frac_2 = zero

    pdf_params_frz(:)%w_1                 = zero
    pdf_params_frz(:)%w_2                 = zero
    pdf_params_frz(:)%varnce_w_1          = zero
    pdf_params_frz(:)%varnce_w_2          = zero
    pdf_params_frz(:)%rt_1                = zero
    pdf_params_frz(:)%rt_2                = zero
    pdf_params_frz(:)%varnce_rt_1         = zero
    pdf_params_frz(:)%varnce_rt_2         = zero
    pdf_params_frz(:)%thl_1               = zero
    pdf_params_frz(:)%thl_2               = zero
    pdf_params_frz(:)%varnce_thl_1        = zero
    pdf_params_frz(:)%varnce_thl_2        = zero
    pdf_params_frz(:)%rrtthl              = zero
    pdf_params_frz(:)%alpha_thl           = zero
    pdf_params_frz(:)%alpha_rt            = zero
    pdf_params_frz(:)%crt_1               = zero
    pdf_params_frz(:)%crt_2               = zero
    pdf_params_frz(:)%cthl_1              = zero
    pdf_params_frz(:)%cthl_2              = zero
    pdf_params_frz(:)%chi_1               = zero
    pdf_params_frz(:)%chi_2               = zero
    pdf_params_frz(:)%stdev_chi_1         = zero
    pdf_params_frz(:)%stdev_chi_2         = zero
    pdf_params_frz(:)%stdev_eta_1         = zero
    pdf_params_frz(:)%stdev_eta_2         = zero
    pdf_params_frz(:)%covar_chi_eta_1     = zero
    pdf_params_frz(:)%covar_chi_eta_2     = zero
    pdf_params_frz(:)%corr_chi_eta_1      = zero
    pdf_params_frz(:)%corr_chi_eta_2      = zero
    pdf_params_frz(:)%rsatl_1             = zero
    pdf_params_frz(:)%rsatl_2             = zero
    pdf_params_frz(:)%rc_1                = zero
    pdf_params_frz(:)%rc_2                = zero
    pdf_params_frz(:)%cloud_frac_1        = zero
    pdf_params_frz(:)%cloud_frac_2        = zero
    pdf_params_frz(:)%mixt_frac           = zero
    pdf_params_frz(:)%ice_supersat_frac_1 = zero
    pdf_params_frz(:)%ice_supersat_frac_2 = zero

    ! Surface fluxes
    wpthlp_sfc = 0.0_core_rknd
    wprtp_sfc  = 0.0_core_rknd
    upwp_sfc   = 0.0_core_rknd
    vpwp_sfc   = 0.0_core_rknd

! ---> h1g, 2010-06-16
! initialize critical relative humidity for liquid and ice nucleation
#ifdef GFDL
    RH_crit = 1.0_core_rknd
#endif
!<--- h1g, 2010-06-16

    ! Passive scalars
    do i = 1, sclr_dim, 1
      wpsclrp_sfc(i)   = 0.0_core_rknd

      sclrm(1:nz,i)         = 0.0_core_rknd
      sclrp2(1:nz,i)        = 0.0_core_rknd
      sclrprtp(1:nz,i)      = 0.0_core_rknd
      sclrpthlp(1:nz,i)     = 0.0_core_rknd
      sclrm_forcing(1:nz,i) = 0.0_core_rknd
      wpsclrp(1:nz,i)         = 0.0_core_rknd
    end do

    do i = 1, edsclr_dim, 1
      wpedsclrp_sfc(i) = 0.0_core_rknd

      edsclrm(1:nz,i)         = 0.0_core_rknd
      edsclrm_forcing(1:nz,i) = 0.0_core_rknd
    end do

    return
  end subroutine setup_prognostic_variables
!-----------------------------------------------------------------------
  subroutine cleanup_prognostic_variables
    implicit none

    ! Prognostic variables
    ! TODO: use a more appropriate condition
    if (allocated(um)) then
      deallocate( um )        ! u wind
      deallocate( vm )        ! v wind

      deallocate( upwp )      ! vertical u momentum flux
      deallocate( vpwp )      ! vertical v momentum flux

      deallocate( up2, vp2 )

      deallocate( thlm )      ! liquid potential temperature

!---> h1g, 2010-06-16
#ifdef GFDL
      deallocate( temp_clubb )
#endif
!<--- h1g, 2010-06-16

      deallocate( rtm )       ! total water mixing ratio
      deallocate( wprtp )     ! w'rt'
      deallocate( wpthlp )    ! w'thl'
      deallocate( wprcp )     ! w'rc'
      deallocate( wp2 )       ! w'^2
      deallocate( wp3 )       ! w'^3
      deallocate( rtp2 )      ! rt'^2
      deallocate( thlp2 )     ! thl'^2
      deallocate( rtpthlp )   ! rt'thl'

      deallocate( p_in_Pa )         ! pressure
      deallocate( exner )           ! exner
      deallocate( rho )             ! density: t points
      deallocate( rho_zm )          ! density: m points
      deallocate( rho_ds_zm )       ! dry, static density: m-levs
      deallocate( rho_ds_zt )       ! dry, static density: t-levs
      deallocate( invrs_rho_ds_zm ) ! inv. dry, static density: m-levs
      deallocate( invrs_rho_ds_zt ) ! inv. dry, static density: t-levs
      deallocate( thv_ds_zm )       ! dry, base-state theta_v: m-levs
      deallocate( thv_ds_zt )       ! dry, base-state theta_v: t-levs

      deallocate( thlm_forcing )    ! thlm large-scale forcing
      deallocate( rtm_forcing )     ! rtm large-scale forcing
      deallocate( um_forcing )      ! u forcing
      deallocate( vm_forcing )      ! v forcing
      deallocate( wprtp_forcing )   ! <w'r_t'> forcing (microphysics)
      deallocate( wpthlp_forcing )  ! <w'th_l'> forcing (microphysics)
      deallocate( rtp2_forcing )    ! <r_t'^2> forcing (microphysics)
      deallocate( thlp2_forcing )   ! <th_l'^2> forcing (microphysics)
      deallocate( rtpthlp_forcing ) ! <r_t'th_l'> forcing (microphysics)

      ! Imposed large scale w

      deallocate( wm_zm )     ! momentum levels
      deallocate( wm_zt )     ! thermodynamic levels

      ! Cloud water variables

      deallocate( rcm )
      deallocate( cloud_frac )
      deallocate( ice_supersat_frac )
      deallocate( rcm_in_layer )
      deallocate( cloud_cover )

      deallocate( sigma_sqd_w )    ! PDF width parameter (momentum levels)

      ! Variable for pdf closure scheme
      deallocate( pdf_params )
      deallocate( pdf_params_frz )

      ! Passive scalars
      deallocate( wpsclrp_sfc, wpedsclrp_sfc )
      deallocate( sclrm )
      deallocate( sclrp2 )
      deallocate( sclrprtp )
      deallocate( sclrpthlp )
      deallocate( sclrm_forcing )
      deallocate( wpsclrp )

      deallocate( edsclrm )
      deallocate( edsclrm_forcing )

!---> h1g, 2010-06-16
#ifdef GFDL
      deallocate( RH_crit )
#endif
! <--- h1g, 2010-06-16
    end if

    return
  end subroutine cleanup_prognostic_variables

end module variables_prognostic_module
