
subroutine scandyn (ztodt,   etadot,  etamid,  grlps1,  grt1,   &
                    grz1,    grd1,    grfu1,   grfv1,   grut1,  &
                    grvt1,   grrh1,   grlps2,  grt2,    grz2,   &
                    grd2,    grfu2,   grfv2,   grut2,   grvt2,  &
                    grrh2,   vcour,   vmax2d,  vmax2dt, detam,  &
                    cwava,   flx_net, t2,      fu,      fv)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! 
! Method: 
! "After coupling" gaussian latitude scan for which some of the physics
! and nonlinear dynamics calculations are completed.  The main loop over
! latitude in this routine is multitasked.
!
! Note: the "ifdef" constructs in this routine are associated with the
! message-passing version of CAM.  Messages are sent  which
! have no relevance to the shared-memory case.  
! 
! Author: 
! Original version:  CCM3
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plat, plev, beglat, endlat, plevp
   use prognostics,  only: u3, v3, q3, t3, div, vort, phis, omga, dpsl, &
                           dpsm, ps, n3m1, n3, n3m2, qminus, pdeld
   use constituents, only: pcnst
   use scanslt,      only: hw1lat
   use rgrid,        only: nlon
   use comspe, only: maxm
   use linemsdyn,    only: linemsdyn_bft, linemsdyn_fft, linemsdyn_aft, &
                     plondfft
   use commap,       only: w
   use qmassa,       only: qmassarun
   use perf_mod
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! Input arguments
!
   real(r8), intent(in) :: ztodt                       ! two delta t unless nstep =0
   real(r8), intent(inout) :: etadot(plon,plevp,beglat:endlat)     ! vertical motion (slt)
   real(r8), intent(in) :: etamid(plev)                ! hybrd coord value at levels
   real(r8), intent(in) :: detam(plev)      
!
! Fourier coefficient arrays which have a latitude index on them for
! multitasking. These arrays are defined in LINEMSDYN and and used in QUAD
! to compute spectral coefficients. They contain a latitude index so
! that the sums over latitude can be performed in a specified order.
!
   real(r8), intent(in)    :: cwava(plat)           ! weight applied to global integrals
   real(r8), intent(in)    :: flx_net(plon,beglat:endlat)         ! net flx from physics
   real(r8), intent(inout) :: t2(plon,plev,beglat:endlat)         ! tot dT/dt to to physics
   real(r8), intent(inout) :: fu(plon,plev,beglat:endlat)         ! u wind tend
   real(r8), intent(inout) :: fv(plon,plev,beglat:endlat)         ! v wind tend
!
! Output arguments
!
   real(r8), intent(out) :: grlps1(2*maxm,(plat+1)/2)      ! sym. undiff. term in lnps eqn.
   real(r8), intent(out) :: grlps2(2*maxm,(plat+1)/2)      ! antisym undiff. term in lnps eqn.
   real(r8), intent(out) :: grt1(2*maxm,plev,(plat+1)/2)   ! sym. undiff. term in t eqn.
   real(r8), intent(out) :: grt2(2*maxm,plev,(plat+1)/2)   ! antisym. undiff. term in t eqn.
   real(r8), intent(out) :: grz1(2*maxm,plev,(plat+1)/2)   ! sym. undiff. term in z eqn.
   real(r8), intent(out) :: grz2(2*maxm,plev,(plat+1)/2)   ! antisym. undiff. term in z eqn.
   real(r8), intent(out) :: grd1(2*maxm,plev,(plat+1)/2)   ! sym. undiff. term in d eqn.
   real(r8), intent(out) :: grd2(2*maxm,plev,(plat+1)/2)   ! antisym. undiff. term in d eqn.
   real(r8), intent(out) :: grfu1(2*maxm,plev,(plat+1)/2)  ! sym. nonlinear terms in u eqn.
   real(r8), intent(out) :: grfu2(2*maxm,plev,(plat+1)/2)  ! antisym. nonlinear terms in u eqn.
   real(r8), intent(out) :: grfv1(2*maxm,plev,(plat+1)/2)  ! sym. nonlinear terms in v eqn.
   real(r8), intent(out) :: grfv2(2*maxm,plev,(plat+1)/2)  ! antisym. nonlinear terms in v eqn.
   real(r8), intent(out) :: grut1(2*maxm,plev,(plat+1)/2)  ! sym. lambda deriv. term in t eqn.
   real(r8), intent(out) :: grut2(2*maxm,plev,(plat+1)/2)  ! antisym. lambda deriv. term in t eqn.
   real(r8), intent(out) :: grvt1(2*maxm,plev,(plat+1)/2)  ! sym. mu derivative term in t eqn.
   real(r8), intent(out) :: grvt2(2*maxm,plev,(plat+1)/2)  ! antisym. mu deriv. term in t eqn.
   real(r8), intent(out) :: grrh1(2*maxm,plev,(plat+1)/2)  ! sym. del**2 term in d eqn.
   real(r8), intent(out) :: grrh2(2*maxm,plev,(plat+1)/2)  ! antisym. del**2 term in d eqn.
   real(r8), intent(out) :: vcour(plev,plat)           ! maximum Courant number in vert.
   real(r8), intent(out) :: vmax2d(plev,plat)          ! max. wind at each level, latitude
   real(r8), intent(out) :: vmax2dt(plev,plat)         ! max. truncated wind at each lvl,lat

! Local variables

   integer irow              ! latitude pair index
   integer lat,latn,lats     ! latitude indices
   integer nlon_fft_in       ! FFT work array inner dimension
   integer nlon_fft_out      ! FFT work array inner dimension
   real(r8) pmid(plon,plev)  ! pressure at model levels
   real(r8) pint(plon,plevp) ! pressure at interfaces
   real(r8) pdel(plon,plev)  ! pressure difference between
   integer :: m              ! constituent index
!
! FFT buffers
!
   real(r8), allocatable:: fftbuf_in(:,:,:,:)          ! fftbuf_in(nlon_fft_in,9,plev,beglat:endlat) 
   real(r8), allocatable:: fftbuf_out(:,:,:,:)         ! fftbuf_out(nlon_fft_out,9,plev,plat)
!
   call t_startf ('scandyn_alloc')
   nlon_fft_in = plondfft
   allocate(fftbuf_in(nlon_fft_in,9,plev,beglat:endlat))

#if ( defined SPMD )
#ifdef NEC_SX
   nlon_fft_out = 2*maxm + 1
#else
   nlon_fft_out = 2*maxm
#endif
   allocate(fftbuf_out(nlon_fft_out,9,plev,plat))
#else
   nlon_fft_out = 1
   allocate(fftbuf_out(1,1,1,1))
#endif
   call t_stopf ('scandyn_alloc')
!
   call t_startf ('linemsdyn_bft')
#ifdef OUTER_OMP
!$OMP PARALLEL DO PRIVATE (LAT)
#endif
   do lat=beglat,endlat

      call linemsdyn_bft (lat, nlon(lat), nlon_fft_in, &
                      ps(1,lat,n3m1), ps(1,lat,n3m2), u3(1,1,lat,n3m1), &
                      u3(1,1,lat,n3m2), v3(1,1,lat,n3m1), v3(1,1,lat,n3m2), t3(1,1,lat,n3m1), t3(1,1,lat,n3m2), &
                      q3(1,1,1,lat,n3m1), etadot(1,1,lat), etamid, &
                      ztodt, vcour(1,lat), vmax2d(1,lat), vmax2dt(1,lat),       &
                      detam, t2(1,1,lat), fu(1,1,lat), fv(1,1,lat),                     &
                      div(1,1,lat,n3m1), vort(1,1,lat,n3m2), div(1,1,lat,n3m2), vort(1,1,lat,n3m1), &
                      phis(1,lat), dpsl(1,lat), dpsm(1,lat), omga(1,1,lat), &
                      cwava(lat), flx_net(1,lat), fftbuf_in(1,1,1,lat) )
   end do
   call t_stopf ('linemsdyn_bft')

   call t_startf ('linemsdyn_fft')
   call linemsdyn_fft (nlon_fft_in,nlon_fft_out,fftbuf_in,fftbuf_out)
   call t_stopf ('linemsdyn_fft')

   call t_startf ('linemsdyn_aft')
!$OMP PARALLEL DO PRIVATE (IROW, LATN, LATS)
   do irow=1,plat/2

      lats = irow
      latn = plat - irow + 1
#if ( defined SPMD )
      call linemsdyn_aft (irow, nlon_fft_out, fftbuf_out(1,1,1,lats), fftbuf_out(1,1,1,latn), &
                      grlps1(1,irow), grt1(1,1,irow), grz1(1,1,irow), grd1(1,1,irow), &
                      grfu1(1,1,irow),  grfv1(1,1,irow),   &
                      grut1(1,1,irow), grvt1(1,1,irow), grrh1(1,1,irow), grlps2(1,irow),grt2(1,1,irow),    &
                      grz2(1,1,irow), grd2(1,1,irow), grfu2(1,1,irow), grfv2(1,1,irow),  grut2(1,1,irow),  &
                      grvt2(1,1,irow), grrh2(1,1,irow) )
#else
      call linemsdyn_aft (irow, nlon_fft_in, fftbuf_in(1,1,1,lats), fftbuf_in(1,1,1,latn), &
                      grlps1(1,irow), grt1(1,1,irow), grz1(1,1,irow), grd1(1,1,irow), &
                      grfu1(1,1,irow),  grfv1(1,1,irow),   &
                      grut1(1,1,irow), grvt1(1,1,irow), grrh1(1,1,irow), grlps2(1,irow),grt2(1,1,irow),    &
                      grz2(1,1,irow), grd2(1,1,irow), grfu2(1,1,irow), grfv2(1,1,irow),  grut2(1,1,irow),  &
                      grvt2(1,1,irow), grrh2(1,1,irow) )
#endif
   end do
   call t_stopf ('linemsdyn_aft')
!
   call t_startf ('scandyn_dealloc')
   deallocate(fftbuf_in)
   deallocate(fftbuf_out)
   call t_stopf ('scandyn_dealloc')

!
   call t_startf ('moisture_mass')
!
! Initialize moisture mass integrals.
!
   hw1lat = 0.0_r8
!
! Calculate total mass of moisture in fields advected
!
#ifdef OUTER_OMP
!$OMP PARALLEL DO PRIVATE (LAT, IROW)
#endif
   do lat=beglat,endlat
      if(lat.le.plat/2) then
         irow = lat
      else
         irow = plat + 1 - lat
      end if
!
! Only pdel is needed pint and pmid are not.
!
      call plevs0 (nlon(lat),plon,plev,ps(1,lat,n3m2), pint, pmid, pdel)
!
! Calculate mass of moisture in field being advected
!

!  q3     is plon,plev,pcnst,beglat:endlat,ptimelevs
!  qminus is plon,plev,pcnst,beglat:endlat
      call qmassarun (cwava(lat),w(irow) ,qminus(1,1,1,lat),pdel    , &
                   hw1lat(1,lat),nlon(lat), q3(1,1,1,lat,n3m2), lat, &
		   pdeld(:,:,lat,n3m2 ))
   end do
   call t_stopf ('moisture_mass')

   return
end subroutine scandyn

