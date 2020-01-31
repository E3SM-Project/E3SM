! $Id: clubbvars.F90 1103 2013-05-14 18:35:02Z minghuai.wang@pnl.gov $
module clubbvars
#ifdef CLUBB_CRM
! Description:
! This module contains variables that exist in CLUBB but not in SAM

  use grid, only: &
    ntracers, &
    nx, &
    ny,&
    nz,&
    nzm,&
    dimx1_s,&
    dimx2_s,&
    dimy1_s,&
    dimy2_s,&
    nxp1,&
    nyp1,&
    YES3D

  use clubb_precision, only: &
    core_rknd ! CLUBB core real kind

  implicit none

  private ! Default Scope

  intrinsic :: selected_real_kind, max

  ! Determines whether to use CLUBB's eddy scalar or high order scalar code on
  ! a tracer in SAM
  ! To enable the passive scalars, set enable_<type> to 1,
  ! and the dimensions for edsclr or sclr will be 1*ntracers.
  integer, private, parameter :: &
    enable_eddy_scalars       = 0, &
    enable_high_order_scalars = 0

  integer, public, parameter :: &
    edsclr_dim = enable_eddy_scalars*ntracers,       & ! Number of eddy scalars
    sclr_dim   = enable_high_order_scalars*ntracers    ! Number of high order scalars

  integer, parameter, public :: &
    tndcy_precision = selected_real_kind( p=12 )

  real(kind = core_rknd), public, dimension(nx, ny, nz) :: &
    upwp,        &! u'w'.                         [m^2/s^2]
    vpwp,        &! u'w'.                         [m^2/s^2]
    up2,         &! u'^2                          [m^2/s^2]
    vp2,         &! v'^2                          [m^2/s^2]
    wprtp,       &! w' r_t'.                      [(m kg)/(s kg)]
    wpthlp,      &! w' th_l'.                     [(m K)/s]
    wprcp,       &! w' r_c'.                      [(kg/kg) m/s]
    wp2,         &! w'^2.                         [m^2/s^2]
    rtp2,        &! r_t'^2.                       [(kg/kg)^2]
    thlp2,       &! th_l'^2.                      [K^2]
    rtpthlp,     &! r_t' th_l'.                   [(kg K)/kg]
    rcm,         &! Cloud water                   [kg/kg]
    cloud_frac,  &! Cloud Fraction.               [-]
    rcm_in_layer,&! rcm in cloud layer            [kg/kg]
    cloud_cover   ! Cloud cover                   [-]

  real, public, dimension(0:nxp1, 1-YES3D:nyp1, nzm) :: &
    khzm,        &! eddy diffusivity on momentum grids                  [m^2/s]
    khzt,        &! eddy diffusivity on thermo grids                    [m^2/s]
    qclvarg,      &! cloud water variance                                [kg^2/kg^2]
    relvarg,     &! relative cloud water variance                        
    accre_enhang  ! accretion enhancement 

    
  real(kind=core_rknd), public, dimension(nx, ny) :: &
     rtm_spurious_source, & ! Spurious source of total water        [kg/kg/s]
     thlm_spurious_source   ! Spurious source of liquid pot. temp.  [K/s]

  ! w'^3 is requires additional ghost points on the x and y dimension,
  ! for the purposes of horizontal advection.  The variables um and vm
  ! require them for the purposes of horizontal interpolation.
  real(kind=core_rknd), public, dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nz) :: &
    wp3,&    ! w'^3.                      [m^3/s^3]
    um, &    ! x-wind                     [m/s]
    vm       ! y-wind                     [m/s]

  real(tndcy_precision), public, dimension(nx, ny, nzm) :: &
    t_tndcy,  & ! CLUBB contribution to moist static energy  [K/s]
    qc_tndcy, & ! CLUBB contribution to liquid water         [kg/kg/s]
    qv_tndcy, & ! CLUBB contribution to vapor water          [kg/kg/s]
    u_tndcy,  & ! CLUBB contribution to x-wind               [m/s^2]
    v_tndcy     ! CLUBB contribution to y-wind               [m/s^2]

  real(tndcy_precision), public, dimension(nx, ny, nzm, ntracers) :: &
    tracer_tndcy ! CLUBB contribution to the tracers [{units vary}/s]

  real(kind=core_rknd), public, dimension(nx,ny,nz,sclr_dim) :: &
    sclrp2,      & ! Passive scalar variance.       [{units vary}^2]
    sclrpthlp,   & ! Passive scalar covariance.     [{units vary} K]
    sclrprtp,    & ! Passive scalar covariance.     [{units vary} kg/kg]
    wpsclrp        ! w'sclr'                        [units vary m/s]

  real(kind=core_rknd), public, dimension(sclr_dim) :: &
    sclr_tol ! Tolerance on passive scalar [units vary]

  real(kind=core_rknd), public, dimension(nz) :: &
    rho_ds_zm,       & ! Dry, static density on momentum levels      [kg/m^3]
    rho_ds_zt,       & ! Dry, static density on thermodynamic levels [kg/m^3]
    invrs_rho_ds_zm, & ! Inv. dry, static density on momentum levels [m^3/kg]
    invrs_rho_ds_zt, & ! Inv. dry, static density on thermo. levels  [m^3/kg]
    thv_ds_zm,       & ! Dry, base-state theta_v on momentum levels  [K]
    thv_ds_zt          ! Dry, base-state theta_v on thermo. levels   [K]

  logical, public :: l_stats_samgrid  ! Stats on sam grid enabled (T/F)

#ifdef CRM
  logical, public :: lrestart_clubb = .false.
#endif
#endif /*CLUBB_CRM*/
end module clubbvars
