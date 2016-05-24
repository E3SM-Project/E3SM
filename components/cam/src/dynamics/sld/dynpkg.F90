
subroutine dynpkg(t2      ,fu      ,fv      ,etamid  , &
                  cwava   ,detam   ,flx_net , ztodt  ,adv_state  )
!-----------------------------------------------------------------------
!
! Purpose:
! driving routines dynamics,and transport
!
! Note that this routine has many "#if ..." constructs.  
! The message-passing model needs to test SPMD since many
! arrays have their space allocated dynamically.
!
! Original version:  CCM3
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plon, plev, plevp, beglat, endlat, plat
  use comspe, only: maxm
  use scanslt,      only: advection_state, plond, platd
  use perf_mod
  implicit none

!------------------------------Arguments--------------------------------
!
  real(r8), intent(inout):: t2(plon,plev,beglat:endlat)  ! temperature tendency
  real(r8), intent(inout):: fu(plon,plev,beglat:endlat)  ! u wind tendency
  real(r8), intent(inout):: fv(plon,plev,beglat:endlat)  ! v wind tendency
  real(r8), intent(in)   :: etamid (plev)        ! vertical coords at midpoints 
  real(r8), intent(in)   :: cwava  (plat)        ! weight applied to global integrals
  real(r8), intent(in)   :: detam  (plev)        ! intervals between vert full levs.
  real(r8), intent(in)   :: flx_net(plon,beglat:endlat)  ! net flux from physics
  real(r8), intent(in)   :: ztodt                ! time step (s)
  type(advection_state), intent(inout):: adv_state
!
!---------------------------Local workspace-----------------------------
!
! Fourier coefficient arrays which have a latitude index on them for
! multitasking. These arrays are defined in LINEMSAC and used in QUAD
! to compute spectral coefficients. They contain a latitude index so
! that the sums over latitude can be performed in a specified order.
!
  real(r8) grlps1(2*maxm,plat/2)      ! ------------------------------
  real(r8) grlps2(2*maxm,plat/2)      ! |
  real(r8) grt1  (2*maxm,plev,plat/2) ! |
  real(r8) grt2  (2*maxm,plev,plat/2) ! |
  real(r8) grq1  (2*maxm,plev,plat/2) ! |
  real(r8) grq2  (2*maxm,plev,plat/2) ! |
  real(r8) grz1  (2*maxm,plev,plat/2) ! |
  real(r8) grz2  (2*maxm,plev,plat/2) ! |
  real(r8) grd1  (2*maxm,plev,plat/2) ! |
  real(r8) grd2  (2*maxm,plev,plat/2) ! |
  real(r8) grfu1 (2*maxm,plev,plat/2) ! |- see linemsac and quad for
  real(r8) grfu2 (2*maxm,plev,plat/2) ! |  definitions
  real(r8) grfv1 (2*maxm,plev,plat/2) ! |
  real(r8) grfv2 (2*maxm,plev,plat/2) ! |
  real(r8) vcour(plev,plat)            ! maximum Courant number in slice
  real(r8) vmax2d (plev,plat)          ! max. wind at each level, latitude
  real(r8) vmax2dt(plev,plat)          ! max. truncated wind at each lvl,lat
!
!----------------------------------------------------------
! SCANDYN Dynamics scan
!----------------------------------------------------------
!
  call t_startf ('scandyn')
  call scandyn(ztodt   ,detam   ,cwava   ,etamid  ,          &
                grlps1  ,grlps2  ,grt1    ,grt2    ,grq1    , &
                grq2    ,grfu1   ,grfu2   ,grfv1   ,grfv2   , &
                fu      ,fv      ,t2      ,flx_net ,          &
                vcour   ,vmax2d  ,vmax2dt ,adv_state )
  call t_stopf  ('scandyn')
!
!----------------------------------------------------------
! Accumulate spectral coefficients
!----------------------------------------------------------
!
  call t_startf ('dyndrv')
  call dyndrv(grlps1  ,grt1    ,grq1    ,grz1    ,grd1    , &
              grfu1   ,grfv1   ,grlps2  ,grt2    ,grq2    , &
              grz2    ,grd2    ,grfu2   ,grfv2   ,vmax2d  , &
              vmax2dt ,vcour   )
  call t_stopf  ('dyndrv')
!
!----------------------------------------------------------
! Second gaussian scan (spectral -> grid)
!----------------------------------------------------------
!
  call t_startf ('scan2')
  call scan2 (ztodt   ,cwava   ,etamid  )
  call t_stopf  ('scan2')

  return
end subroutine dynpkg
