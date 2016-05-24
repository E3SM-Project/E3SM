
subroutine scandyn(ztodt   ,detam   ,cwava   ,etamid  ,          &
                    grlps1  ,grlps2  ,grt1    ,grt2    ,grq1    , &
                    grq2    ,grfu1   ,grfu2   ,grfv1   ,grfv2   , &
                    grfu    ,grfv    ,t2      ,flx_net , &
                    vcour   ,vmax2d  ,vmax2dt ,adv_state )
!-----------------------------------------------------------------------
!
! Purpose:
! Driving routine for semi-lagrangian transport and SLD dynamics.
! Set up  FFT and combine terms in preparation for Fourier -> spectral
! quadrature.
! 
! The latitude loop in this routine is multitasked.
! The naming convention is as follows:
!  - prefix gr contains grid point values before FFT and Fourier
!    coefficients after
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plev, plat, plon, beglat, endlat
  use prognostics,  only: ps, u3, v3, q3, t3, etadot, div, parrsld, &
                          tarrsld, n3, n3m1, prhs, trhs, vrhs, urhs, &
                          omga, phism, phisl, phis, dpsm, dpsl, ql, &
                          tl, ed1, qm, tm
  use rgrid,        only: nlon
  use comspe,       only: maxm
  use scanslt,      only: advection_state, slt_run, slt_run_setup
  use perf_mod

  implicit none

!------------------------------Arguments--------------------------------
!
  real(r8), intent(in)   :: ztodt                     ! twice the timestep unless nstep=0
  real(r8), intent(in)   :: detam  (plev)             ! delta eta at levels
  real(r8), intent(in)   :: cwava  (plat)             ! weight for global water vapor int.
  real(r8), intent(in)   :: etamid (plev)             ! eta at levels

  real(r8), intent(out)   :: grlps1(2*maxm,plat/2)      ! ------------------------------
  real(r8), intent(out)   :: grlps2(2*maxm,plat/2)      ! |
  real(r8), intent(out)   :: grt1  (2*maxm,plev,plat/2) ! |
  real(r8), intent(out)   :: grt2  (2*maxm,plev,plat/2) ! |
  real(r8), intent(out)   :: grq1  (2*maxm,plev,plat/2) ! |
  real(r8), intent(out)   :: grq2  (2*maxm,plev,plat/2) ! |- see quad for definitions
  real(r8), intent(out)   :: grfu1 (2*maxm,plev,plat/2) ! |
  real(r8), intent(out)   :: grfu2 (2*maxm,plev,plat/2) ! |
  real(r8), intent(out)   :: grfv1 (2*maxm,plev,plat/2) ! |
  real(r8), intent(out)   :: grfv2 (2*maxm,plev,plat/2) ! ------------------------------
  real(r8), intent(inout) :: grfu  (plon,plev,beglat:endlat)   ! nonlinear term - u momentum eqn
  real(r8), intent(inout) :: grfv  (plon,plev,beglat:endlat)   ! nonlinear term - v momentum eqn
  real(r8), intent(inout) :: t2    (plon,plev,beglat:endlat)    ! tot dT/dt to to physics
  real(r8), intent(in)   :: flx_net(plon,beglat:endlat)         ! net flx from physics
  real(r8), intent(out)  :: vcour  (plev,plat)          ! maximum Courant number in vert.
  real(r8), intent(out)  :: vmax2d (plev,plat)          ! max. wind at each level, latitude
  real(r8), intent(out)  :: vmax2dt(plev,plat)          ! max. truncated wind at each lvl,lat
  type(advection_state), intent(inout) :: adv_state
!
!---------------------------Local workspace-----------------------------
!
  integer :: lat        ! latitude index
  real(r8) :: lnpssld(plon,plev,beglat:endlat)
  real(r8) :: prhssld(plon,plev,beglat:endlat)
!
!-----------------------------------------------------------------------
!
  call t_startf ('slt_run_setup')
  call slt_run_setup( adv_state )
  call t_stopf  ('slt_run_setup')

  call t_startf ('linemsdyn')
!$OMP PARALLEL DO PRIVATE (LAT)

  do lat=beglat,endlat
     call linemsdyn(lat                     ,ps     (1,lat,n3)   ,u3   (1,1,lat,n3)  , &
                       u3     (1,1,lat,n3m1) ,v3     (1,1,lat,n3)  ,                    &
                    v3        (1,1,lat,n3m1) ,t3     (1,1  ,lat,n3),t3   (1,1,lat,n3m1), &
                       q3     (1,1,1,lat,n3) ,etadot (1,1,lat)  ,                    &
                                              etamid              ,ztodt             , &
                       vcour  (1,lat)       ,vmax2d(1,lat)       ,vmax2dt(1,lat)    , &
                       detam               ,                                          &
                    ed1       (1,1,lat)     ,grfu   (1,1,lat)    ,grfv (1,1,lat)    , &
                       lnpssld(1,1,lat)      ,prhssld(1,1,lat)     ,                    &
                    tarrsld   (1,1,lat)     ,parrsld(1,1,lat)    ,t2   (1,1,lat)    , &
                       div    (1,1,lat,n3)  ,tl     (1,1,lat)    ,                    &
                    tm        (1,1,lat)     ,ql     (1,1,lat)    ,qm   (1,1,lat)    , &
                       dpsl   (1,lat)       ,dpsm   (1,lat)      ,                    &
                    phis      (1,lat)       ,phisl  (1,lat)      ,phism(1,lat)      , &
                       omga   (1,1,lat)     ,                                         &
                       urhs   (1,1,lat)     ,vrhs   (1,1,lat)    ,                    &
                    trhs      (1,1,lat)     ,prhs   (1,1,lat)    ,nlon (lat),         &
                       cwava(lat), flx_net(1,lat))

  end do
  call t_stopf ('linemsdyn')

  call t_startf ('sltrun')
  call slt_run( ztodt   ,detam   ,cwava   ,etamid  ,          &
                grlps1  ,grlps2  ,grt1    ,grt2    ,grq1    , &
                grq2    ,grfu1   ,grfu2   ,grfv1   ,grfv2   , &
                grfu    ,grfv    ,lnpssld ,prhssld ,adv_state )
  call t_stopf  ('sltrun')
!
!
  return
end subroutine scandyn
