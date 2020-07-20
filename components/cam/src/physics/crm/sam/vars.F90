
module vars
  use grid
  use params, only: crm_rknd
#ifdef CRM
#ifdef MODAL_AERO
  use modal_aero_data,   only: ntot_amode
#endif
#endif

  implicit none
  public
  !--------------------------------------------------------------------
  ! prognostic variables:
  real(crm_rknd), allocatable :: u   (:,:,:,:) ! x-wind
  real(crm_rknd), allocatable :: v   (:,:,:,:) ! y-wind
  real(crm_rknd), allocatable :: w   (:,:,:,:) ! z-wind
  real(crm_rknd), allocatable :: t   (:,:,:,:) ! liquid/ice water static energy

  !--------------------------------------------------------------------
  ! diagnostic variables:
  real(crm_rknd), allocatable :: p       (:,:,:,:)     ! perturbation pressure (from Poison eq)
  real(crm_rknd), allocatable :: tabs    (:,:,:,:)                 ! temperature
  real(crm_rknd), allocatable :: qv      (:,:,:,:)                ! water vapor
  real(crm_rknd), allocatable :: qcl     (:,:,:,:)                ! liquid water  (condensate)
  real(crm_rknd), allocatable :: qpl     (:,:,:,:)                ! liquid water  (precipitation)
  real(crm_rknd), allocatable :: qci     (:,:,:,:)                ! ice water  (condensate)
  real(crm_rknd), allocatable :: qpi     (:,:,:,:)                ! ice water  (precipitation)
  real(crm_rknd), allocatable :: tke2    (:,:,:,:)   ! SGS TKE
  real(crm_rknd), allocatable :: tk2     (:,:,:,:) ! SGS eddyviscosity

  !--------------------------------------------------------------------
  ! time-tendencies for prognostic variables

  real(crm_rknd), allocatable :: dudt   (:,:,:,:,:)
  real(crm_rknd), allocatable :: dvdt   (:,:,:,:,:)
  real(crm_rknd), allocatable :: dwdt   (:,:,:,:,:)

  !----------------------------------------------------------------
  ! Temporary storage array:
  real(crm_rknd), allocatable :: misc(:,:,:,:)

  !------------------------------------------------------------------
  ! fluxes at the top and bottom of the domain:
  real(crm_rknd), allocatable :: fluxbu  (:,:,:)
  real(crm_rknd), allocatable :: fluxbv  (:,:,:)
  real(crm_rknd), allocatable :: fluxbt  (:,:,:)
  real(crm_rknd), allocatable :: fluxbq  (:,:,:)
  real(crm_rknd), allocatable :: fluxtu  (:,:,:)
  real(crm_rknd), allocatable :: fluxtv  (:,:,:)
  real(crm_rknd), allocatable :: fluxtt  (:,:,:)
  real(crm_rknd), allocatable :: fluxtq  (:,:,:)
  real(crm_rknd), allocatable :: fzero   (:,:,:)
  real(crm_rknd), allocatable :: precsfc (:,:,:) ! surface precip. rate
  real(crm_rknd), allocatable :: precssfc(:,:,:) ! surface ice precip. rate

  !-----------------------------------------------------------------
  ! profiles
  real(crm_rknd), allocatable :: t0   (:,:)
  real(crm_rknd), allocatable :: q0   (:,:)
  real(crm_rknd), allocatable :: qv0  (:,:)
  real(crm_rknd), allocatable :: tabs0(:,:)
  real(crm_rknd), allocatable :: tv0  (:,:)
  real(crm_rknd), allocatable :: u0   (:,:)
  real(crm_rknd), allocatable :: v0   (:,:)
  real(crm_rknd), allocatable :: tg0  (:,:)
  real(crm_rknd), allocatable :: qg0  (:,:)
  real(crm_rknd), allocatable :: ug0  (:,:)
  real(crm_rknd), allocatable :: vg0  (:,:)
  real(crm_rknd), allocatable :: p0   (:,:)
  real(crm_rknd), allocatable :: tke0 (:,:)
  real(crm_rknd), allocatable :: t01  (:,:)
  real(crm_rknd), allocatable :: q01  (:,:)
  real(crm_rknd), allocatable :: qp0  (:,:)
  real(crm_rknd), allocatable :: qn0  (:,:)

  !-----------------------------------------------------------------
  ! reference vertical profiles:
  real(crm_rknd), allocatable :: prespot(:,:)  ! (1000./pres)**R/cp
  real(crm_rknd), allocatable :: rho    (:,:)   ! air density at pressure levels,kg/m3
  real(crm_rknd), allocatable :: rhow   (:,:)   ! air density at vertical velocity levels,kg/m3
  real(crm_rknd), allocatable :: bet    (:,:)   ! = ggr/tv0
  real(crm_rknd), allocatable :: gamaz  (:,:) ! ggr/cp*z
  real(crm_rknd), allocatable :: wsub   (:,:)   ! Large-scale subsidence velocity,m/s
  real(crm_rknd), allocatable :: qtend  (:,:) ! Large-scale tendency for total water
  real(crm_rknd), allocatable :: ttend  (:,:) ! Large-scale tendency for temp.
  real(crm_rknd), allocatable :: utend  (:,:) ! Large-scale tendency for u
  real(crm_rknd), allocatable :: vtend  (:,:) ! Large-scale tendency for v

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !---------------------------------------------------------------------
  !  Horizontally varying stuff (as a function of xy)
  !
  real(crm_rknd), allocatable :: sstxy    (:,:,:) !  surface temperature xy-distribution
  real(crm_rknd), allocatable :: fcory    (:,:)      !  Coriolis parameter xy-distribution
  real(crm_rknd), allocatable :: fcorzy   (:,:)      !  z-Coriolis parameter xy-distribution
  real(crm_rknd), allocatable :: latitude (:,:,:)      ! latitude (degrees,:)
  real(crm_rknd), allocatable :: longitude(:,:,:)      ! longitude(degrees,:)
  real(crm_rknd), allocatable :: prec_xy  (:,:,:) ! mean precip. rate for outout
  real(crm_rknd), allocatable :: pw_xy    (:,:,:)   ! precipitable water
  real(crm_rknd), allocatable :: cw_xy    (:,:,:)   ! cloud water path
  real(crm_rknd), allocatable :: iw_xy    (:,:,:)   ! ice water path
  real(crm_rknd), allocatable :: cld_xy   (:,:,:)   ! cloud frequency
  real(crm_rknd), allocatable :: u200_xy  (:,:,:) ! u-wind at 200 mb
  real(crm_rknd), allocatable :: usfc_xy  (:,:,:) ! u-wind at at the surface
  real(crm_rknd), allocatable :: v200_xy  (:,:,:) ! v-wind at 200 mb
  real(crm_rknd), allocatable :: vsfc_xy  (:,:,:) ! v-wind at the surface
  real(crm_rknd), allocatable :: w500_xy  (:,:,:) ! w at 500 mb

  !----------------------------------------------------------------------
  ! Vertical profiles of quantities sampled for statitistics purposes:
  real(crm_rknd), allocatable :: w_max(:)
  real(crm_rknd), allocatable :: u_max(:)

  real(crm_rknd), allocatable :: twsb(:,:)
  real(crm_rknd), allocatable :: precflux(:,:)
  real(crm_rknd), allocatable :: uwle(:,:)
  real(crm_rknd), allocatable :: uwsb(:,:)
  real(crm_rknd), allocatable :: vwle(:,:)
  real(crm_rknd), allocatable :: vwsb(:,:)
  real(crm_rknd), allocatable :: tkelediss(:,:)
  real(crm_rknd), allocatable :: tdiff(:,:)
  real(crm_rknd), allocatable :: tlat(:,:)
  real(crm_rknd), allocatable :: tlatqi(:,:)
  real(crm_rknd), allocatable :: qifall(:,:)
  real(crm_rknd), allocatable :: qpfall(:,:)

  ! energy conservation diagnostics:
  real(8), allocatable :: total_water_evap(:)
  real(8), allocatable :: total_water_prec(:)

  real(crm_rknd), allocatable :: CF3D(:,:,:,:)  ! Cloud fraction
  ! =1.0 when there is no fractional cloudiness scheme
  ! = cloud fraction produced by fractioal cloudiness scheme when avaiable

  ! 850 mbar horizontal winds
  real(crm_rknd), allocatable :: u850_xy(:,:,:) ! zonal velocity at 850 mb
  real(crm_rknd), allocatable :: v850_xy(:,:,:) ! meridional velocity at 850 mb

  ! Surface pressure
  real(crm_rknd), allocatable :: psfc_xy(:,:,:) ! pressure (in millibar) at lowest grid point

  ! Saturated water vapor path, useful for computing column relative humidity
  real(crm_rknd), allocatable :: swvp_xy(:,:,:)  ! saturated water vapor path (wrt water)

  ! Cloud and echo top heights, and cloud top temperature (instantaneous)
  real(crm_rknd), allocatable :: cloudtopheight(:,:,:)
  real(crm_rknd), allocatable :: echotopheight (:,:,:)
  real(crm_rknd), allocatable :: cloudtoptemp  (:,:,:)

  ! END UW ADDITIONS
  !===========================================================================
#if (defined CRM && defined MODAL_AERO)
  real(crm_rknd), allocatable :: naer (:,:,:)     ! Aerosol number concentration [/m3]
  real(crm_rknd), allocatable :: vaer (:,:,:)     ! aerosol volume concentration [m3/m3]
  real(crm_rknd), allocatable :: hgaer(:,:,:)    ! hygroscopicity of aerosol mode
#endif

  public :: allocate_vars
  public :: deallocate_vars
#if defined(_OPENMP)
  public :: update_device_vars
  public :: update_host_vars
#endif
contains
  subroutine allocate_vars(ncrms)
#if defined(_OPENACC)
    use openacc_utils
#endif
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) :: zero

    allocate( u(ncrms,dimx1_u:dimx2_u,dimy1_u:dimy2_u,nzm)  )
    allocate( v(ncrms,dimx1_v:dimx2_v,dimy1_v:dimy2_v,nzm)  )
    allocate( w(ncrms,dimx1_w:dimx2_w,dimy1_w:dimy2_w,nz )  )
    allocate( t(ncrms,dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm)  )
    allocate( p(ncrms,0:nx, (1-YES3D):ny, nzm)      )
    allocate( tabs(ncrms,nx, ny, nzm)                  )
    allocate( qv(ncrms,nx, ny, nzm)                 )
    allocate( qcl(ncrms,nx, ny, nzm)                 )
    allocate( qpl(ncrms,nx, ny, nzm)                 )
    allocate( qci(ncrms,nx, ny, nzm)                 )
    allocate( qpi(ncrms,nx, ny, nzm)                 )
    allocate( tke2(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)    )
    allocate( tk2  (ncrms,0:nxp1, (1-YES3D):nyp1, nzm)  )
    allocate( dudt(ncrms,nxp1, ny, nzm, 3) )
    allocate( dvdt(ncrms,nx, nyp1, nzm, 3) )
    allocate( dwdt(ncrms,nx, ny  , nz,  3) )
    allocate( misc(ncrms,nx, ny, nz) )
    allocate( fluxbu(ncrms,nx,ny) )
    allocate( fluxbv(ncrms,nx,ny) )
    allocate( fluxbt(ncrms,nx,ny) )
    allocate( fluxbq(ncrms,nx,ny) )
    allocate( fluxtu(ncrms,nx,ny) )
    allocate( fluxtv(ncrms,nx,ny) )
    allocate( fluxtt(ncrms,nx,ny) )
    allocate( fluxtq(ncrms,nx,ny) )
    allocate( fzero(ncrms,nx,ny) )
    allocate( precsfc(ncrms,nx,ny)  )
    allocate( precssfc(ncrms,nx,ny)  )
    allocate( t0(ncrms,nzm) )
    allocate( q0(ncrms,nzm) )
    allocate( qv0(ncrms,nzm) )
    allocate( tabs0(ncrms,nzm) )
    allocate( tv0(ncrms,nzm) )
    allocate( u0(ncrms,nzm) )
    allocate( v0(ncrms,nzm) )
    allocate( tg0(ncrms,nzm) )
    allocate( qg0(ncrms,nzm) )
    allocate( ug0(ncrms,nzm) )
    allocate( vg0(ncrms,nzm) )
    allocate( p0(ncrms,nzm) )
    allocate( tke0(ncrms,nzm) )
    allocate( t01(ncrms,nzm) )
    allocate( q01(ncrms,nzm) )
    allocate( qp0(ncrms,nzm) )
    allocate( qn0(ncrms,nzm) )
    allocate( prespot(ncrms,nzm)   )
    allocate( rho(ncrms,nzm)     )
    allocate( rhow(ncrms,nz )    )
    allocate( bet(ncrms,nzm)     )
    allocate( gamaz(ncrms,nzm)  )
    allocate( wsub(ncrms,nz )    )
    allocate( qtend(ncrms,nzm)  )
    allocate( ttend(ncrms,nzm)  )
    allocate( utend(ncrms,nzm)  )
    allocate( vtend(ncrms,nzm)  )
    allocate( sstxy    (ncrms,0:nx,(1-YES3D):ny)   )
    allocate( fcory(ncrms,0:ny)       )
    allocate( fcorzy(ncrms,ny)       )
    allocate( latitude(ncrms,nx,ny)        )
    allocate( longitude(ncrms,nx,ny)        )
    allocate( prec_xy(ncrms,nx,ny)  )
    allocate( pw_xy(ncrms,nx,ny)    )
    allocate( cw_xy(ncrms,nx,ny)    )
    allocate( iw_xy(ncrms,nx,ny)    )
    allocate( cld_xy(ncrms,nx,ny)    )
    allocate( u200_xy(ncrms,nx,ny)  )
    allocate( usfc_xy(ncrms,nx,ny)  )
    allocate( v200_xy(ncrms,nx,ny)  )
    allocate( vsfc_xy(ncrms,nx,ny)  )
    allocate( w500_xy(ncrms,nx,ny)  )
    allocate( twsb(ncrms,nz) )
    allocate( precflux(ncrms,nz) )
    allocate( uwle(ncrms,nz) )
    allocate( uwsb(ncrms,nz) )
    allocate( vwle(ncrms,nz) )
    allocate( vwsb(ncrms,nz) )
    allocate( tkelediss(ncrms,nz) )
    allocate( tdiff    (ncrms,nz) )
    allocate( tlat(ncrms,nz) )
    allocate( tlatqi(ncrms,nz) )
    allocate( qifall(ncrms,nz) )
    allocate( qpfall(ncrms,nz) )
    allocate( cf3d(ncrms,1:nx, 1:ny, 1:nzm)   )
    allocate( u850_xy(ncrms,nx,ny)  )
    allocate( v850_xy(ncrms,nx,ny)  )
    allocate( psfc_xy(ncrms,nx,ny)  )
    allocate( swvp_xy(ncrms,nx,ny)   )
    allocate( cloudtopheight(ncrms,nx,ny) )
    allocate( echotopheight(ncrms,nx,ny) )
    allocate( cloudtoptemp(ncrms,nx,ny) )
    allocate( u_max(ncrms) )
    allocate( w_max(ncrms) )
    allocate( total_water_evap(ncrms) )
    allocate( total_water_prec(ncrms) )
#if (defined CRM && defined MODAL_AERO)
    allocate( naer(ncrms,nzm, ntot_amode) )
    allocate( vaer(ncrms,nzm, ntot_amode) )
    allocate( hgaer(ncrms,nzm, ntot_amode) )
#endif
#if defined(_OPENACC)
    call prefetch( u )
    call prefetch( v )
    call prefetch( w )
    call prefetch( t )
    call prefetch( p )
    call prefetch( tabs )
    call prefetch( qv )
    call prefetch( qcl )
    call prefetch( qpl )
    call prefetch( qci )
    call prefetch( qpi )
    call prefetch( tke2 )
    call prefetch( tk2 )
    call prefetch( dudt )
    call prefetch( dvdt )
    call prefetch( dwdt )
    call prefetch( misc )
    call prefetch( fluxbu )
    call prefetch( fluxbv )
    call prefetch( fluxbt )
    call prefetch( fluxbq )
    call prefetch( fluxtu )
    call prefetch( fluxtv )
    call prefetch( fluxtt )
    call prefetch( fluxtq )
    call prefetch( fzero )
    call prefetch( precsfc )
    call prefetch( precssfc )
    call prefetch( t0 )
    call prefetch( q0 )
    call prefetch( qv0 )
    call prefetch( tabs0 )
    call prefetch( tv0 )
    call prefetch( u0 )
    call prefetch( v0 )
    call prefetch( tg0 )
    call prefetch( qg0 )
    call prefetch( ug0 )
    call prefetch( vg0 )
    call prefetch( p0 )
    call prefetch( tke0 )
    call prefetch( t01 )
    call prefetch( q01 )
    call prefetch( qp0 )
    call prefetch( qn0 )
    call prefetch( prespot )
    call prefetch( rho )
    call prefetch( rhow )
    call prefetch( bet )
    call prefetch( gamaz )
    call prefetch( wsub )
    call prefetch( qtend )
    call prefetch( ttend )
    call prefetch( utend )
    call prefetch( vtend )
    call prefetch( sstxy )
    call prefetch( fcory )
    call prefetch( fcorzy )
    call prefetch( latitude )
    call prefetch( longitude )
    call prefetch( prec_xy )
    call prefetch( pw_xy )
    call prefetch( cw_xy )
    call prefetch( iw_xy )
    call prefetch( cld_xy )
    call prefetch( u200_xy )
    call prefetch( usfc_xy )
    call prefetch( v200_xy )
    call prefetch( vsfc_xy )
    call prefetch( w500_xy )
    call prefetch( twsb )
    call prefetch( precflux )
    call prefetch( uwle )
    call prefetch( uwsb )
    call prefetch( vwle )
    call prefetch( vwsb )
    call prefetch( tkelediss )
    call prefetch( tdiff )
    call prefetch( tlat )
    call prefetch( tlatqi )
    call prefetch( qifall )
    call prefetch( qpfall )
    call prefetch( CF3D )
    call prefetch( u850_xy )
    call prefetch( v850_xy )
    call prefetch( psfc_xy )
    call prefetch( swvp_xy )
    call prefetch( cloudtopheight )
    call prefetch( echotopheight )
    call prefetch( cloudtoptemp )
    call prefetch( u_max )
    call prefetch( w_max )
    call prefetch( total_water_evap )
    call prefetch( total_water_prec )
#if (defined CRM && defined MODAL_AERO)
    call prefetch( naer )
    call prefetch( vaer )
    call prefetch( hgaer  )
#endif
#elif defined(_OPENMP)
    !$omp target enter data map(alloc: u )
    !$omp target enter data map(alloc: v )
    !$omp target enter data map(alloc: w )
    !$omp target enter data map(alloc: t )
    !$omp target enter data map(alloc: p )
    !$omp target enter data map(alloc: tabs )
    !$omp target enter data map(alloc: qv )
    !$omp target enter data map(alloc: qcl )
    !$omp target enter data map(alloc: qpl )
    !$omp target enter data map(alloc: qci )
    !$omp target enter data map(alloc: qpi )
    !$omp target enter data map(alloc: tke2 )
    !$omp target enter data map(alloc: tk2 )
    !$omp target enter data map(alloc: dudt )
    !$omp target enter data map(alloc: dvdt )
    !$omp target enter data map(alloc: dwdt )
    !$omp target enter data map(alloc: misc )
    !$omp target enter data map(alloc: fluxbu ) 
    !$omp target enter data map(alloc: fluxbv )
    !$omp target enter data map(alloc: fluxbt ) 
    !$omp target enter data map(alloc: fluxbq )
    !$omp target enter data map(alloc: fluxtu )
    !$omp target enter data map(alloc: fluxtv )
    !$omp target enter data map(alloc: fluxtt )
    !$omp target enter data map(alloc: fluxtq )
    !$omp target enter data map(alloc: fzero )
    !$omp target enter data map(alloc: precsfc )
    !$omp target enter data map(alloc: precssfc )
    !$omp target enter data map(alloc: t0 )
    !$omp target enter data map(alloc: q0 )
    !$omp target enter data map(alloc: qv0 )
    !$omp target enter data map(alloc: tabs0 )
    !$omp target enter data map(alloc: tv0 )
    !$omp target enter data map(alloc: u0 )
    !$omp target enter data map(alloc: v0 )
    !$omp target enter data map(alloc: tg0 )
    !$omp target enter data map(alloc: qg0 )
    !$omp target enter data map(alloc: ug0 )
    !$omp target enter data map(alloc: vg0 )
    !$omp target enter data map(alloc: p0 )
    !$omp target enter data map(alloc: tke0 )
    !$omp target enter data map(alloc: t01 )
    !$omp target enter data map(alloc: q01 )
    !$omp target enter data map(alloc: qp0 )
    !$omp target enter data map(alloc: qn0 )
    !$omp target enter data map(alloc: prespot )
    !$omp target enter data map(alloc: rho )
    !$omp target enter data map(alloc: rhow )
    !$omp target enter data map(alloc: bet )
    !$omp target enter data map(alloc: gamaz )
    !$omp target enter data map(alloc: wsub )
    !$omp target enter data map(alloc: qtend )
    !$omp target enter data map(alloc: ttend )
    !$omp target enter data map(alloc: utend )
    !$omp target enter data map(alloc: vtend )
    !$omp target enter data map(alloc: sstxy )
    !$omp target enter data map(alloc: fcory )
    !$omp target enter data map(alloc: fcorzy )
    !$omp target enter data map(alloc: latitude )
    !$omp target enter data map(alloc: prec_xy )
    !$omp target enter data map(alloc: pw_xy )
    !$omp target enter data map(alloc: cw_xy )
    !$omp target enter data map(alloc: iw_xy )
    !$omp target enter data map(alloc: cld_xy )
    !$omp target enter data map(alloc: u200_xy )
    !$omp target enter data map(alloc: usfc_xy )
    !$omp target enter data map(alloc: v200_xy )
    !$omp target enter data map(alloc: vsfc_xy )
    !$omp target enter data map(alloc: w500_xy )
    !$omp target enter data map(alloc: twsb )
    !$omp target enter data map(alloc: precflux )
    !$omp target enter data map(alloc: uwle )
    !$omp target enter data map(alloc: uwsb )
    !$omp target enter data map(alloc: vwle )
    !$omp target enter data map(alloc: vwsb )
    !$omp target enter data map(alloc: tkelediss )
    !$omp target enter data map(alloc: tdiff )
    !$omp target enter data map(alloc: tlat )
    !$omp target enter data map(alloc: tlatqi )
    !$omp target enter data map(alloc: qifall )
    !$omp target enter data map(alloc: qpfall )
    !$omp target enter data map(alloc: cf3d )
    !$omp target enter data map(alloc: u850_xy )
    !$omp target enter data map(alloc: v850_xy )
    !$omp target enter data map(alloc: psfc_xy )
    !$omp target enter data map(alloc: swvp_xy )
    !$omp target enter data map(alloc: cloudtopheight )
    !$omp target enter data map(alloc: echotopheight )
    !$omp target enter data map(alloc: cloudtoptemp )
    !$omp target enter data map(alloc: u_max )
    !$omp target enter data map(alloc: w_max )
    !$omp target enter data map(alloc: total_water_evap )
    !$omp target enter data map(alloc: total_water_prec )
#if (defined CRM && defined MODAL_AERO)
    !$omp target enter data map(alloc: naer )
    !$omp target enter data map(alloc: vaer )
    !$omp target enter data map(alloc: hgaer )
#endif
#endif
    zero = 0.0_crm_rknd
    u = zero
    v = zero
    w = zero
    t = zero
    p = zero
    tabs = zero
    qv = zero
    qcl = zero
    qpl = zero
    qci = zero
    qpi = zero
    tke2 = zero
    tk2 = zero
    dudt = zero
    dvdt = zero
    dwdt = zero
    misc = zero
    fluxbu = zero
    fluxbv = zero
    fluxbt = zero
    fluxbq = zero
    fluxtu = zero
    fluxtv = zero
    fluxtt = zero
    fluxtq = zero
    fzero = zero
    precsfc = zero
    precssfc = zero
    t0 = zero
    q0 = zero
    qv0 = zero
    tabs0 = zero
    tv0 = zero
    u0 = zero
    v0 = zero
    tg0 = zero
    qg0 = zero
    ug0 = zero
    vg0 = zero
    p0 = zero
    tke0 = zero
    t01 = zero
    q01 = zero
    qp0 = zero
    qn0 = zero
    prespot = zero
    rho = zero
    rhow = zero
    bet = zero
    gamaz = zero
    wsub = zero
    qtend = zero
    ttend = zero
    utend = zero
    vtend = zero
    sstxy = zero
    fcory = zero
    fcorzy = zero
    latitude = zero
    longitude = zero
    prec_xy = zero
    pw_xy = zero
    cw_xy = zero
    iw_xy = zero
    cld_xy = zero
    u200_xy = zero
    usfc_xy = zero
    v200_xy = zero
    vsfc_xy = zero
    w500_xy = zero
    twsb = zero
    precflux = zero
    uwle = zero
    uwsb = zero
    vwle = zero
    vwsb = zero
    tkelediss = zero
    tdiff = zero
    tlat = zero
    tlatqi = zero
    qifall = zero
    qpfall = zero
    CF3D = 1.
    u850_xy = zero
    v850_xy = zero
    psfc_xy = zero
    swvp_xy = zero
    cloudtopheight = zero
    echotopheight = zero
    cloudtoptemp = zero
    u_max = zero
    w_max = zero
    total_water_evap = zero
    total_water_prec = zero
#if (defined CRM && defined MODAL_AERO)
    naer = zero
    vaer = zero
    hgaer = zero
#endif
  end subroutine allocate_vars

#if defined(_OPENMP)
  subroutine update_device_vars()
    implicit none
    !$omp target update to( u )
    !$omp target update to( v )
    !$omp target update to( w )
    !$omp target update to( t )
    !$omp target update to( p )
    !$omp target update to( tabs )
    !$omp target update to( qv )
    !$omp target update to( qcl )
    !$omp target update to( qpl )
    !$omp target update to( qci )
    !$omp target update to( qpi )
    !$omp target update to( tke2 )
    !$omp target update to( tk2 )
    !$omp target update to( dudt )
    !$omp target update to( dvdt )
    !$omp target update to( dwdt )
    !$omp target update to( misc )
    !$omp target update to( fluxbu )
    !$omp target update to( fluxbv )
    !$omp target update to( fluxbt )
    !$omp target update to( fluxbq )
    !$omp target update to( fluxtu )
    !$omp target update to( fluxtv )
    !$omp target update to( fluxtt )
    !$omp target update to( fluxtq )
    !$omp target update to( fzero )
    !$omp target update to( precsfc )
    !$omp target update to( precssfc )
    !$omp target update to( t0 )
    !$omp target update to( q0 )
    !$omp target update to( qv0 )
    !$omp target update to( tabs0 )
    !$omp target update to( tv0 )
    !$omp target update to( u0 )
    !$omp target update to( v0 )
    !$omp target update to( tg0 )
    !$omp target update to( qg0 )
    !$omp target update to( ug0 )
    !$omp target update to( vg0 )
    !$omp target update to( p0 )
    !$omp target update to( tke0 )
    !$omp target update to( t01 )
    !$omp target update to( q01 )
    !$omp target update to( qp0 )
    !$omp target update to( qn0 )
    !$omp target update to( prespot )
    !$omp target update to( rho )
    !$omp target update to( rhow )
    !$omp target update to( bet )
    !$omp target update to( gamaz )
    !$omp target update to( wsub )
    !$omp target update to( qtend )
    !$omp target update to( ttend )
    !$omp target update to( utend )
    !$omp target update to( vtend )
    !$omp target update to( sstxy )
    !$omp target update to( fcory )
    !$omp target update to( fcorzy )
    !$omp target update to( latitude )
    !$omp target update to( longitude )
    !$omp target update to( prec_xy )
    !$omp target update to( pw_xy )
    !$omp target update to( cw_xy )
    !$omp target update to( iw_xy )
    !$omp target update to( cld_xy )
    !$omp target update to( u200_xy )
    !$omp target update to( vsfc_xy )
    !$omp target update to( w500_xy )
    !$omp target update to( twsb )
    !$omp target update to( precflux )
    !$omp target update to( uwle )
    !$omp target update to( uwsb )
    !$omp target update to( vwle )
    !$omp target update to( vwsb )
    !$omp target update to( tkelediss )
    !$omp target update to( tdiff )
    !$omp target update to( tlat )
    !$omp target update to( tlatqi )
    !$omp target update to( qifall )
    !$omp target update to( qpfall )
    !$omp target update to( CF3D )
    !$omp target update to( u850_xy )
    !$omp target update to( v850_xy )
    !$omp target update to( psfc_xy )
    !$omp target update to( swvp_xy )
    !$omp target update to( cloudtopheight )
    !$omp target update to( echotopheight )
    !$omp target update to( cloudtoptemp )
    !$omp target update to( u_max )
    !$omp target update to( w_max )
    !$omp target update to( total_water_evap )
    !$omp target update to( total_water_prec )
#if (defined CRM && defined MODAL_AERO)
    !$omp target update to( naer )
    !$omp target update to( vaer )
    !$omp target update to( hgaer )
#endif
  end subroutine update_device_vars

  subroutine update_host_vars()
    implicit none
    !$omp target update from( u )
    !$omp target update from( v )
    !$omp target update from( w )
    !$omp target update from( t )
    !$omp target update from( p )
    !$omp target update from( tabs )
    !$omp target update from( qv )
    !$omp target update from( qcl )
    !$omp target update from( qpl )
    !$omp target update from( qci )
    !$omp target update from( qpi )
    !$omp target update from( tke2 )
    !$omp target update from( tk2 )
    !$omp target update from( dudt )
    !$omp target update from( dvdt )
    !$omp target update from( dwdt )
    !$omp target update from( misc )
    !$omp target update from( fluxbu )
    !$omp target update from( fluxbv )
    !$omp target update from( fluxbt )
    !$omp target update from( fluxbq )
    !$omp target update from( fluxtu )
    !$omp target update from( fluxtv )
    !$omp target update from( fluxtt )
    !$omp target update from( fluxtq )
    !$omp target update from( fzero )
    !$omp target update from( precsfc )
    !$omp target update from( precssfc )
    !$omp target update from( t0 )
    !$omp target update from( q0 )
    !$omp target update from( qv0 )
    !$omp target update from( tabs0 )
    !$omp target update from( tv0 )
    !$omp target update from( u0 )
    !$omp target update from( v0 )
    !$omp target update from( tg0 )
    !$omp target update from( qg0 )
    !$omp target update from( ug0 )
    !$omp target update from( vg0 )
    !$omp target update from( p0 )
    !$omp target update from( tke0 )
    !$omp target update from( t01 )
    !$omp target update from( q01 )
    !$omp target update from( qp0 )
    !$omp target update from( qn0 )
    !$omp target update from( prespot )
    !$omp target update from( rho )
    !$omp target update from( rhow )
    !$omp target update from( bet )
    !$omp target update from( gamaz )
    !$omp target update from( wsub )
    !$omp target update from( qtend )
    !$omp target update from( ttend )
    !$omp target update from( utend )
    !$omp target update from( vtend )
    !$omp target update from( sstxy )
    !$omp target update from( fcory )
    !$omp target update from( fcorzy )
    !$omp target update from( latitude )
    !$omp target update from( longitude )
    !$omp target update from( prec_xy )
    !$omp target update from( pw_xy )
    !$omp target update from( cw_xy )
    !$omp target update from( iw_xy )
    !$omp target update from( cld_xy )
    !$omp target update from( u200_xy )
    !$omp target update from( vsfc_xy )
    !$omp target update from( w500_xy )
    !$omp target update from( twsb )
    !$omp target update from( precflux )
    !$omp target update from( uwle )
    !$omp target update from( uwsb )
    !$omp target update from( vwle )
    !$omp target update from( vwsb )
    !$omp target update from( tkelediss )
    !$omp target update from( tdiff )
    !$omp target update from( tlat )
    !$omp target update from( tlatqi )
    !$omp target update from( qifall )
    !$omp target update from( qpfall )
    !$omp target update from( CF3D )
    !$omp target update from( u850_xy )
    !$omp target update from( v850_xy )
    !$omp target update from( psfc_xy )
    !$omp target update from( swvp_xy )
    !$omp target update from( cloudtopheight )
    !$omp target update from( echotopheight )
    !$omp target update from( cloudtoptemp )
    !$omp target update from( u_max )
    !$omp target update from( w_max )
    !$omp target update from( total_water_evap )
    !$omp target update from( total_water_prec )
#if (defined CRM && defined MODAL_AERO)
    !$omp target update from( naer )
    !$omp target update from( vaer )
    !$omp target update from( hgaer )
#endif
  end subroutine update_host_vars
#endif

  subroutine deallocate_vars()
    implicit none
#if defined(_OPENMP)
    !$omp target exit data map(delete: u )
    !$omp target exit data map(delete: v )
    !$omp target exit data map(delete: w )
    !$omp target exit data map(delete: t )
    !$omp target exit data map(delete: p )
    !$omp target exit data map(delete: tabs )
    !$omp target exit data map(delete: qv )
    !$omp target exit data map(delete: qcl )
    !$omp target exit data map(delete: qpl )
    !$omp target exit data map(delete: qci )
    !$omp target exit data map(delete: qpi )
    !$omp target exit data map(delete: tke2 )
    !$omp target exit data map(delete: tk2 )
    !$omp target exit data map(delete: dudt )
    !$omp target exit data map(delete: dvdt )
    !$omp target exit data map(delete: dwdt )
    !$omp target exit data map(delete: misc )
    !$omp target exit data map(delete: fluxbu )
    !$omp target exit data map(delete: fluxbv )
    !$omp target exit data map(delete: fluxbt )
    !$omp target exit data map(delete: fluxbq )
    !$omp target exit data map(delete: fluxtu )
    !$omp target exit data map(delete: fluxtv )
    !$omp target exit data map(delete: fluxtt )
    !$omp target exit data map(delete: fluxtq )
    !$omp target exit data map(delete: fzero )
    !$omp target exit data map(delete: precsfc )
    !$omp target exit data map(delete: precssfc )
    !$omp target exit data map(delete: t0 )
    !$omp target exit data map(delete: q0 )
    !$omp target exit data map(delete: qv0 )
    !$omp target exit data map(delete: tabs0 )
    !$omp target exit data map(delete: tv0 )
    !$omp target exit data map(delete: u0 )
    !$omp target exit data map(delete: v0 )
    !$omp target exit data map(delete: tg0 )
    !$omp target exit data map(delete: qg0 )
    !$omp target exit data map(delete: ug0 )
    !$omp target exit data map(delete: vg0 )
    !$omp target exit data map(delete: p0 )
    !$omp target exit data map(delete: tke0 )
    !$omp target exit data map(delete: t01 )
    !$omp target exit data map(delete: q01 )
    !$omp target exit data map(delete: qp0 )
    !$omp target exit data map(delete: qn0 )
    !$omp target exit data map(delete: prespot )
    !$omp target exit data map(delete: rho )
    !$omp target exit data map(delete: rhow )
    !$omp target exit data map(delete: bet )
    !$omp target exit data map(delete: gamaz )
    !$omp target exit data map(delete: wsub )
    !$omp target exit data map(delete: qtend )
    !$omp target exit data map(delete: ttend )
    !$omp target exit data map(delete: utend )
    !$omp target exit data map(delete: vtend )
    !$omp target exit data map(delete: sstxy )
    !$omp target exit data map(delete: fcory )
    !$omp target exit data map(delete: fcorzy )
    !$omp target exit data map(delete: latitude )
    !$omp target exit data map(delete: longitude )
    !$omp target exit data map(delete: prec_xy )
    !$omp target exit data map(delete: pw_xy )
    !$omp target exit data map(delete: cw_xy )
    !$omp target exit data map(delete: iw_xy )
    !$omp target exit data map(delete: cld_xy )
    !$omp target exit data map(delete: u200_xy )
    !$omp target exit data map(delete: usfc_xy )
    !$omp target exit data map(delete: v200_xy )
    !$omp target exit data map(delete: vsfc_xy )
    !$omp target exit data map(delete: w500_xy )
    !$omp target exit data map(delete: twsb )
    !$omp target exit data map(delete: precflux )
    !$omp target exit data map(delete: uwle )
    !$omp target exit data map(delete: uwsb )
    !$omp target exit data map(delete: vwle )
    !$omp target exit data map(delete: vwsb )
    !$omp target exit data map(delete: tkelediss )
    !$omp target exit data map(delete: tdiff )
    !$omp target exit data map(delete: tlat )
    !$omp target exit data map(delete: tlatqi )
    !$omp target exit data map(delete: qifall )
    !$omp target exit data map(delete: qpfall )
    !$omp target exit data map(delete: CF3D )
    !$omp target exit data map(delete: u850_xy )
    !$omp target exit data map(delete: v850_xy )
    !$omp target exit data map(delete: psfc_xy )
    !$omp target exit data map(delete: swvp_xy )
    !$omp target exit data map(delete: cloudtopheight )
    !$omp target exit data map(delete: echotopheight )
    !$omp target exit data map(delete: cloudtoptemp )
    !$omp target exit data map(delete: u_max )
    !$omp target exit data map(delete: w_max )
    !$omp target exit data map(delete: total_water_evap )
    !$omp target exit data map(delete: total_water_prec )
#if (defined CRM && defined MODAL_AERO)
    !$omp target exit data map(delete: naer )
    !$omp target exit data map(delete: vaer )
    !$omp target exit data map(delete: hgaer )
#endif
#endif
    deallocate( u )
    deallocate( v )
    deallocate( w )
    deallocate( t )
    deallocate( p )
    deallocate( tabs )
    deallocate( qv )
    deallocate( qcl )
    deallocate( qpl )
    deallocate( qci )
    deallocate( qpi )
    deallocate( tke2 )
    deallocate( tk2 )
    deallocate( dudt )
    deallocate( dvdt )
    deallocate( dwdt )
    deallocate( misc )
    deallocate( fluxbu )
    deallocate( fluxbv )
    deallocate( fluxbt )
    deallocate( fluxbq )
    deallocate( fluxtu )
    deallocate( fluxtv )
    deallocate( fluxtt )
    deallocate( fluxtq )
    deallocate( fzero )
    deallocate( precsfc )
    deallocate( precssfc )
    deallocate( t0 )
    deallocate( q0 )
    deallocate( qv0 )
    deallocate( tabs0 )
    deallocate( tv0 )
    deallocate( u0 )
    deallocate( v0 )
    deallocate( tg0 )
    deallocate( qg0 )
    deallocate( ug0 )
    deallocate( vg0 )
    deallocate( p0 )
    deallocate( tke0 )
    deallocate( t01 )
    deallocate( q01 )
    deallocate( qp0 )
    deallocate( qn0 )
    deallocate( prespot )
    deallocate( rho )
    deallocate( rhow )
    deallocate( bet )
    deallocate( gamaz )
    deallocate( wsub )
    deallocate( qtend )
    deallocate( ttend )
    deallocate( utend )
    deallocate( vtend )
    deallocate( sstxy )
    deallocate( fcory )
    deallocate( fcorzy )
    deallocate( latitude )
    deallocate( longitude )
    deallocate( prec_xy )
    deallocate( pw_xy )
    deallocate( cw_xy )
    deallocate( iw_xy )
    deallocate( cld_xy )
    deallocate( u200_xy )
    deallocate( usfc_xy )
    deallocate( v200_xy )
    deallocate( vsfc_xy )
    deallocate( w500_xy )
    deallocate( twsb )
    deallocate( precflux )
    deallocate( uwle )
    deallocate( uwsb )
    deallocate( vwle )
    deallocate( vwsb )
    deallocate( tkelediss )
    deallocate( tdiff )
    deallocate( tlat )
    deallocate( tlatqi )
    deallocate( qifall )
    deallocate( qpfall )
    deallocate( CF3D )
    deallocate( u850_xy )
    deallocate( v850_xy )
    deallocate( psfc_xy )
    deallocate( swvp_xy )
    deallocate( cloudtopheight )
    deallocate( echotopheight )
    deallocate( cloudtoptemp )
    deallocate( u_max )
    deallocate( w_max )
    deallocate( total_water_evap )
    deallocate( total_water_prec )
#if (defined CRM && defined MODAL_AERO)
    deallocate( naer )
    deallocate( vaer )
    deallocate( hgaer  )
#endif
end subroutine deallocate_vars
end module vars
