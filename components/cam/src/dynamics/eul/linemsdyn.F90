
module linemsdyn

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Control non-linear dynamical terms, FFT and combine terms
! in preparation for Fourier -> spectral quadrature.
! 
! Method: 
! The naming convention is as follows:
!  - prefix gr contains grid point values before FFT and Fourier
!     coefficients after
!  - t, q, d, z and ps refer to temperature, specific humidity,
!     divergence, vorticity and surface pressure
!  - "1" suffix to an array => symmetric component current latitude pair
!  - "2" suffix to an array => antisymmetric component.
!
! Author: 
! Original version:  CCM3
! Modified:          P. Worley, October 2002
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plev, plevp, plat, beglat, endlat
   use spmd_utils,   only: iam
   use perf_mod
   implicit none

   private
!
! Public interfaces
!
   public linemsdyn_bft    ! Before FFT
   public linemsdyn_fft    ! FFT
   public linemsdyn_aft    ! After FFT
!
! Public data
!
   integer, public, parameter :: plondfft  = plon + 2      ! Length needed for FFT
   integer, public, parameter :: plndlvfft = plondfft*plev ! Length of multilevel 3-d field slice

!
!-----------------------------------------------------------------------
!

contains

!-----------------------------------------------------------------------

subroutine linemsdyn_bft(                                          &
                     lat     ,nlon    ,nlon_fft,                   &
                     psm1    ,psm2    ,u3m1    , &
                     u3m2    ,v3m1    ,v3m2    ,t3m1    ,t3m2    , &
                     q3m1    ,etadot  ,etamid  ,                   &
                     ztodt   , vcour   ,vmax   ,vmaxt   ,          &
                     detam   ,t2      ,fu      ,fv      ,          &
                     divm1   ,vortm2  ,divm2   ,vortm1  ,phis    , &
                     dpsl    ,dpsm    ,omga    ,cwava   ,flx_net , &
                     fftbuf             )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Control non-linear dynamical terms and fill FFT buffer 
! in preparation for Fourier -> spectral quadrature.
! 
! Author: 
! Original version:  CCM3
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$

   use constituents, only: pcnst
   use pspect,       only: ptrm, ptrn
   use scanslt,      only: engy1lat
   use commap,       only: clat, tau, w
   use cam_history,  only: outfld
   use time_manager, only: get_step_size
   use hycoef,       only : hypd, hypi
   use cam_control_mod, only : adiabatic
   use eul_control_mod, only : eul_nsplit
!
! Input arguments
!     
   integer lat                              ! latitude index for S->N storage
   integer nlon
   integer, intent(in) :: nlon_fft          ! first dimension of FFT work array

   real(r8), intent(in) :: psm1(plon)       ! surface pressure (time n)
   real(r8), intent(in) :: psm2(plon)       ! surface pressure (time n-1)
   real(r8), intent(in) :: u3m1(plon,plev)  ! u-wind (time n)
   real(r8), intent(in) :: u3m2(plon,plev)  ! u-wind (time n-1)
   real(r8), intent(in) :: v3m1(plon,plev)  ! v-wind (time n)
   real(r8), intent(in) :: v3m2(plon,plev)  ! v-wind (time n-1)
   real(r8), intent(in) :: t3m1(plon,plev)  ! temperature (time n)
   real(r8), intent(in) :: q3m1(plon,plev,pcnst)   ! constituent conc(time n: h2o first)
   real(r8), intent(inout) :: etadot(plon,plevp) ! vertical motion (3-d used by slt)
   real(r8), intent(in) :: etamid(plev)     ! midpoint values of eta (a+b)
   real(r8), intent(in) :: ztodt            ! 2*timestep unless nstep = 0
   real(r8), intent(in) :: detam(plev)      ! maximum Courant number in vert.
!     
! Input/Output arguments
!     
   real(r8), intent(inout) :: t2(plon,plev)    ! t tend
   real(r8), intent(inout) :: fu(plon,plev)    ! nonlinear term - u momentum eqn.
   real(r8), intent(inout) :: fv(plon,plev)    ! nonlinear term - v momentum eqn.
   real(r8), intent(inout) :: divm1(plon,plev)
   real(r8), intent(inout) :: vortm2(plon,plev)
   real(r8), intent(inout) :: divm2(plon,plev)
   real(r8), intent(inout) :: vortm1(plon,plev)
   real(r8), intent(inout) :: phis(plon)
   real(r8), intent(inout) :: dpsl(plon)
   real(r8), intent(inout) :: dpsm(plon)
   real(r8), intent(inout) :: omga(plon,plev)
   real(r8), intent(inout) :: t3m2(plon,plev)  ! temperature (time n-1)
   real(r8), intent(in)    :: cwava                  ! weight for global water vapor int.
   real(r8), intent(in)    :: flx_net(plon)          ! net flux from physics
!     
! Output arguments
!     
   real(r8), intent(out) :: fftbuf(nlon_fft,9,plev) ! buffer used for in-place FFTs
   real(r8), intent(out) :: vcour(plev)      ! maximum Courant number in vert.
   real(r8), intent(out) :: vmax(plev)       ! maximum wind speed squared (m^2/s^2)
   real(r8), intent(out) :: vmaxt(plev)      ! maximum truncated wind speed (m^2/s^2)
!     
!---------------------------Local workspace-----------------------------
!     
   real(r8) :: dtime          ! timestep size
   real(r8) :: bpstr(plon)    ! 
   real(r8) pmid(plon,plev)   ! pressure at model levels (time n)
   real(r8) rpmid(plon,plev)  ! 1./pmid
   real(r8) pint(plon,plevp)  ! pressure at model interfaces (n  )
   real(r8) pdel(plon,plev)   ! pdel(k)   = pint  (k+1)-pint  (k)
   real(r8) rpdel(plon,plev)  ! 1./pdel
   real(r8) tdyn(plon,plev)   ! temperature for dynamics
   real(r8) logpsm1(plon)     ! log(psm1)
   real(r8) logpsm2(plon)     ! log(psm2)
   real(r8) engy(plon,plev)   ! kinetic energy
   real(r8) vat  (plon,plev)  ! Vertical advection of temperature
   real(r8) ktoop(plon,plev)  ! (Kappa*T)*(omega/P)
   real(r8) ut(plon,plev)     ! (u*T) - heat flux - zonal
   real(r8) vt(plon,plev)     ! (v*T) - heat flux - meridional
   real(r8) drhs(plon,plev)   ! RHS of divergence eqn. (del^2 term)
   real(r8) lvcour            ! local vertical courant number
   real(r8) dtdz              ! dt/detam(k)
   real(r8) ddivdt(plon,plev) ! temporary workspace
   real(r8) ddpn(plon)        ! complete sum of d*delta p
   real(r8) vpdsn(plon)       ! complete sum V dot grad(ln(ps)) delta b
   real(r8) dpslat(plon,plev) ! Pressure gradient term 
   real(r8) dpslon(plon,plev) ! Pressure gradient term 
   real(r8) coslat            ! cosine(latitude)
   real(r8) rcoslat           ! 1./cosine(latitude)
   real(r8) rhypi             ! 1./hypi(plevp)

   real(r8) wind              ! u**2 + v**2 (m/s)
   real(r8) utfac             ! asymmetric truncation factor for courant calculation
   real(r8) vtfac             ! asymmetric truncation factor for courant calculation

   real(r8) tmp               ! accumulator
   integer i,k,kk             ! longitude,level,constituent indices
   integer, parameter :: tdyndex = 1     ! indices into fftbuf 
   integer, parameter :: fudex = 2
   integer, parameter :: fvdex = 3
   integer, parameter :: utdex = 4
   integer, parameter :: vtdex = 5
   integer, parameter :: drhsdex = 6
   integer, parameter :: vortdyndex = 7
   integer, parameter :: divdyndex = 8
   integer, parameter :: bpstrdex = 9
!
! This group of arrays are glued together via equivalence to exbuf for
! communication from LINEMSBC.
!
!
!-----------------------------------------------------------------------
!
!
! Compute maximum wind speed this latitude (used in Courant number estimate)
!
   if (ptrm .lt. ptrn) then
      utfac = real(ptrm,r8)/real(ptrn,r8)
      vtfac = 1._r8
   else if (ptrn .lt. ptrm) then
      utfac = 1._r8
      vtfac = real(ptrn,r8)/real(ptrm,r8) 
   else if (ptrn .eq. ptrm) then
      utfac = 1._r8
      vtfac = 1._r8
   end if

!$OMP PARALLEL DO PRIVATE (K, I, WIND)
   do k=1,plev
      vmax(k) = 0._r8
      vmaxt(k) = 0._r8
      do i=1,nlon
         wind = u3m2(i,k)**2 + v3m2(i,k)**2
         vmax(k) = max(wind,vmax(k))
!
! Change to Courant limiter for non-triangular truncations.
!
         wind = utfac*u3m2(i,k)**2 + vtfac*v3m2(i,k)**2
         vmaxt(k) = max(wind,vmaxt(k))
      end do
   end do
!
! Variables needed in tphysac
!
   coslat = cos(clat(lat))
   rcoslat = 1._r8/coslat
!
! Set current time pressure arrays for model levels etc.
!
   call plevs0(nlon,plon,plev,psm1,pint,pmid,pdel)
!$OMP PARALLEL DO PRIVATE (K, I)
   do k=1,plev
      do i=1,nlon
         rpmid(i,k) = 1._r8/pmid(i,k)
         rpdel(i,k) = 1._r8/pdel(i,k)
      end do
   end do
!
! Accumulate statistics for diagnostic print
!
   call stats(lat,     pint,    pdel,      psm1,   &
              vortm1,  divm1,   t3m1,      q3m1(:,:,1), nlon  )
!
! Compute log(surface pressure) for use by grmult and when adding tendency.
!
!$OMP PARALLEL DO PRIVATE (I)
   do i=1,nlon
      logpsm1(i) = log(psm1(i))
      logpsm2(i) = log(psm2(i))
   end do
!     
! Compute integrals
!     
   call plevs0(nlon,plon,plev,psm2,pint,pmid,pdel)
   call engy_te (cwava,w(lat),t3m2,u3m2,v3m2,phis    ,pdel, tmp  ,nlon)
   engy1lat(lat) = tmp 
   call plevs0(nlon,plon,plev,psm1,pint,pmid,pdel)
!
! Include top/bottom flux integral to energy integral
!
   call flxint  (w(lat) ,flx_net ,tmp  ,nlon )
   engy1lat(lat) = engy1lat(lat) + tmp *ztodt
!
! Calculate non-linear terms in tendencies
!
   if (adiabatic) t2(:,:) = 0._r8
   call outfld('FU      ',fu    ,plon,lat)
   call outfld('FV      ',fv    ,plon,lat)
   call grmult(rcoslat ,divm1     ,q3m1(1,1,1),t3m1   ,u3m1    , &
               v3m1    ,vortm1    ,t3m2    ,phis    ,dpsl    , &
               dpsm    ,omga    ,pdel    ,pint(1,plevp),logpsm2, &
               logpsm1 ,rpmid   ,rpdel   ,fu      ,fv      , &
               t2      ,ut      ,vt      ,drhs    ,pmid    , &
               etadot  ,etamid  ,engy    ,ddpn    ,vpdsn   , &
               dpslon  ,dpslat  ,vat     ,ktoop   ,nlon    )
!
! Add tendencies to previous timestep values of surface pressure,
! temperature, and (if spectral transport) moisture.  Store *log* surface
! pressure in bpstr array for transform to spectral space.
!
   rhypi = 1._r8/hypi(plevp)
!$OMP PARALLEL DO PRIVATE (K, I)
   do k=1,plev
      do i=1,nlon
         ddivdt(i,k) = ztodt*(0.5_r8*divm2(i,k) - divm1(i,k))
         tdyn(i,k) = t3m2(i,k) + ztodt*t2(i,k)
      end do
   end do

!$OMP PARALLEL DO PRIVATE (I, K)
   do i=1,nlon
      bpstr(i) = logpsm2(i) - ztodt*(vpdsn(i)+ddpn(i))/psm1(i)
      do k=1,plev
         bpstr(i) = bpstr(i) - ddivdt(i,k)*hypd(k)*rhypi
      end do
   end do

!$OMP PARALLEL DO PRIVATE (K, KK, I)
   do k=1,plev
      do kk=1,plev
         do i=1,nlon
            tdyn(i,k) = tdyn(i,k) - ddivdt(i,kk)*tau(kk,k)
         end do
      end do
   end do

!
! Compute maximum vertical Courant number this latitude.
!
   dtime = get_step_size()/eul_nsplit
   vcour(:) = 0._r8
!$OMP PARALLEL DO PRIVATE (K, DTDZ, I, LVCOUR)
   do k=2,plev
      dtdz = dtime/detam(k-1)
      do i=1,nlon
         lvcour = abs(etadot(i,k))*dtdz
         vcour(k) = max(lvcour,vcour(k))
      end do
   end do

   call outfld('ETADOT  ',etadot,plon,lat)
   call outfld('VAT     ',vat   ,plon,lat)
   call outfld('KTOOP   ',ktoop ,plon,lat)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Apply cos(lat) to momentum terms before fft
!
!$OMP PARALLEL DO PRIVATE (K, I)
   do k=1,plev
      do i=1,nlon
         fu(i,k) = coslat*fu(i,k)
         fv(i,k) = coslat*fv(i,k)
         ut(i,k) = coslat*ut(i,k)
         vt(i,k) = coslat*vt(i,k)
      end do
   end do

!
! Copy fields into FFT buffer
!
!$OMP PARALLEL DO PRIVATE (K, I)
   do k=1,plev
      do i=1,nlon
!
! undifferentiated terms
         fftbuf(i,tdyndex,k) = tdyn(i,k)
! longitudinally and latitudinally differentiated terms
         fftbuf(i,fudex,k)   = fu(i,k)
         fftbuf(i,fvdex,k)   = fv(i,k)
         fftbuf(i,utdex,k)   = ut(i,k)
         fftbuf(i,vtdex,k)   = vt(i,k)
         fftbuf(i,drhsdex,k) = drhs(i,k)
! vort,div
         fftbuf(i,vortdyndex,k) = vortm2(i,k)
         fftbuf(i,divdyndex,k)  = divm2(i,k)
!
      enddo
   enddo
! ps
   do i=1,nlon
      fftbuf(i,bpstrdex,1) = bpstr(i)
   enddo

   return
end subroutine linemsdyn_bft

!-----------------------------------------------------------------------

subroutine linemsdyn_fft(nlon_fft,nlon_fft2,fftbuf,fftbuf2)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute FFT of non-linear dynamical terms
! in preparation for Fourier -> spectral quadrature.
! 
! Author: 
! Original version:  CCM3
! Modified:          P. Worley, September 2002
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$

   use pmgrid,  only: plon, plat
   use rgrid,   only: nlon
   use eul_control_mod, only : trig, ifax
#if (defined SPMD)
   use mpishorthand, only: mpicom
#endif

!     
! Input arguments
!     
   integer, intent(in) :: nlon_fft         ! first dimension of first FFT work array
   integer, intent(in) :: nlon_fft2        ! first dimension of second FFT work array
!     
! Input/Output arguments
!     
   real(r8), intent(inout) :: fftbuf(nlon_fft,9,plev,beglat:endlat) 
                            ! buffer used for in-place FFTs
!     
! Output arguments
!     
#if (defined SPMD)
   real(r8), intent(out) :: fftbuf2(nlon_fft2,9,plev,plat) 
                            ! buffer for returning reorderd Fourier coefficients
#else
   real(r8), intent(in) :: fftbuf2(1) 
                            ! buffer unused
#endif
!     
!---------------------------Local workspace-----------------------------
!     
! The "work" array has a different size requirement depending upon whether
! the proprietary Cray assembly language version of the FFT library
! routines, or the all-Fortran version, is being used.
!     
#if ( ! defined USEFFTLIB )
   real(r8) work((plon+1)*plev*9)
#else 
   real(r8) work((plon+1)*pcray) ! workspace array for fft991
#endif
   integer lat               ! latitude index
   integer inc               ! increment for fft991
   integer isign             ! flag indicates transform direction
   integer ntr               ! number of transforms to perform
   integer k                 ! vertical level index
!
   inc = 1
   isign = -1
#ifdef OUTER_OMP
!$OMP PARALLEL DO PRIVATE (LAT, NTR, K, WORK)
#endif
   do lat=beglat,endlat
      ntr = 8
!$OMP PARALLEL DO PRIVATE (K, WORK)
      do k=1,plev
         fftbuf(nlon(lat)+1:nlon_fft,:,k,lat) = 0.0_r8
         call fft991(fftbuf(1,1,k,lat)     ,work    ,trig(1,lat),ifax(1,lat),inc     ,&
                     nlon_fft ,nlon(lat)   ,ntr     ,isign   )
      enddo
      ntr = 1
      fftbuf(nlon(lat)+1:nlon_fft,9,1,lat) = 0.0_r8
      call fft991(fftbuf(1,9,1,lat)     ,work    ,trig(1,lat),ifax(1,lat),inc     ,&
                  nlon_fft ,nlon(lat)   ,ntr     ,isign   )
   enddo
!
#if ( defined SPMD )
!
!  reorder Fourier coefficients
!
   call t_barrierf ('sync_realloc4a', mpicom)
   call t_startf('realloc4a')
   call realloc4a(nlon_fft, nlon_fft2, fftbuf, fftbuf2)
   call t_stopf('realloc4a')
#endif

   return
end subroutine linemsdyn_fft

!-----------------------------------------------------------------------

subroutine linemsdyn_aft(                                          &
                     irow    ,nlon_fft,fftbufs ,fftbufn ,          &
                     grlps1  ,grt1    ,grz1    ,grd1    ,          &
                     grfu1   ,grfv1   ,grut1   ,grvt1   ,grrh1   , &
                     grlps2  ,grt2    ,grz2    ,grd2    ,grfu2   , &
                     grfv2   ,grut2   ,grvt2   ,grrh2              )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Combine terms in preparation for Fourier -> spectral quadrature.
! 
! Author: 
! Original version:  CCM3
! Modified:          P. Worley, September 2002
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$

#if (defined SPMD)
   use comspe, only: numm, maxm
#else
   use comspe, only: maxm
   use rgrid, only: nmmax
#endif
! Input arguments
!     
   integer, intent(in) :: irow                ! latitude pair index
   integer, intent(in) :: nlon_fft            ! first dimension of FFT work arrays

   real(r8), intent(in) :: fftbufs(nlon_fft,9,plev) ! southern latitude Fourier coefficients
   real(r8), intent(in) :: fftbufn(nlon_fft,9,plev) ! northern latitude Fourier coefficients
!     
! Output arguments
!     
   real(r8), intent(out) :: grlps1(2*maxm)  ! sym. undiff. term in lnps eqn.
   real(r8), intent(out) :: grlps2(2*maxm)  ! antisym undiff. term in lnps eqn.
   real(r8), intent(out) :: grt1(2*maxm,plev) ! sym. undiff. term in t eqn.
   real(r8), intent(out) :: grt2(2*maxm,plev) ! antisym. undiff. term in t eqn.
   real(r8), intent(out) :: grz1(2*maxm,plev) ! sym. undiff. term in z eqn.
   real(r8), intent(out) :: grz2(2*maxm,plev) ! antisym. undiff. term in z eqn.
   real(r8), intent(out) :: grd1(2*maxm,plev) ! sym. undiff. term in d eqn.
   real(r8), intent(out) :: grd2(2*maxm,plev) ! antisym. undiff. term in d eqn.
   real(r8), intent(out) :: grfu1(2*maxm,plev) ! sym. nonlinear terms in u eqn.
   real(r8), intent(out) :: grfu2(2*maxm,plev) ! antisym. nonlinear terms in u eqn.
   real(r8), intent(out) :: grfv1(2*maxm,plev) ! sym. nonlinear terms in v eqn.
   real(r8), intent(out) :: grfv2(2*maxm,plev) ! antisym. nonlinear terms in v eqn.
   real(r8), intent(out) :: grut1(2*maxm,plev) ! sym. lambda deriv. term in t eqn.
   real(r8), intent(out) :: grut2(2*maxm,plev) ! antisym. lambda deriv. term in t eqn.
   real(r8), intent(out) :: grvt1(2*maxm,plev) ! sym. mu derivative term in t eqn.
   real(r8), intent(out) :: grvt2(2*maxm,plev) ! antisym. mu deriv. term in t eqn.
   real(r8), intent(out) :: grrh1(2*maxm,plev) ! sym. del**2 term in d eqn.
   real(r8), intent(out) :: grrh2(2*maxm,plev) ! antisym. del**2 term in d eqn.
!     
!---------------------------Local workspace-----------------------------
!     
   integer i,k            ! longitude,level indices
   integer mlength        ! number of wavenumbers
   integer, parameter :: tdyndex = 1     ! indices into fftbuf 
   integer, parameter :: fudex = 2
   integer, parameter :: fvdex = 3
   integer, parameter :: utdex = 4
   integer, parameter :: vtdex = 5
   integer, parameter :: drhsdex = 6
   integer, parameter :: vortdyndex = 7
   integer, parameter :: divdyndex = 8
   integer, parameter :: bpstrdex = 9
!
#if (defined SPMD)
   mlength = numm(iam)
#else
   mlength = nmmax(irow)
#endif
   do k=1,plev
!cdir loopchg
      do i=1,2*mlength

         grt1(i,k) = 0.5_r8*(fftbufn(i,tdyndex,k)+fftbufs(i,tdyndex,k))
         grt2(i,k) = 0.5_r8*(fftbufn(i,tdyndex,k)-fftbufs(i,tdyndex,k))

         grz1(i,k) = 0.5_r8*(fftbufn(i,vortdyndex,k)+fftbufs(i,vortdyndex,k))
         grz2(i,k) = 0.5_r8*(fftbufn(i,vortdyndex,k)-fftbufs(i,vortdyndex,k))

         grd1(i,k) = 0.5_r8*(fftbufn(i,divdyndex,k)+fftbufs(i,divdyndex,k))
         grd2(i,k) = 0.5_r8*(fftbufn(i,divdyndex,k)-fftbufs(i,divdyndex,k))

         grfu1(i,k) = 0.5_r8*(fftbufn(i,fudex,k)+fftbufs(i,fudex,k))
         grfu2(i,k) = 0.5_r8*(fftbufn(i,fudex,k)-fftbufs(i,fudex,k))

         grfv1(i,k) = 0.5_r8*(fftbufn(i,fvdex,k)+fftbufs(i,fvdex,k))
         grfv2(i,k) = 0.5_r8*(fftbufn(i,fvdex,k)-fftbufs(i,fvdex,k))

         grut1(i,k) = 0.5_r8*(fftbufn(i,utdex,k)+fftbufs(i,utdex,k))
         grut2(i,k) = 0.5_r8*(fftbufn(i,utdex,k)-fftbufs(i,utdex,k))

         grvt1(i,k) = 0.5_r8*(fftbufn(i,vtdex,k)+fftbufs(i,vtdex,k))
         grvt2(i,k) = 0.5_r8*(fftbufn(i,vtdex,k)-fftbufs(i,vtdex,k))

         grrh1(i,k) = 0.5_r8*(fftbufn(i,drhsdex,k)+fftbufs(i,drhsdex,k))
         grrh2(i,k) = 0.5_r8*(fftbufn(i,drhsdex,k)-fftbufs(i,drhsdex,k))

      end do
   end do

!cdir altcode=(loopcnt)
   do i=1,2*mlength
      grlps1(i) = 0.5_r8*(fftbufn(i,bpstrdex,1)+fftbufs(i,bpstrdex,1))
      grlps2(i) = 0.5_r8*(fftbufn(i,bpstrdex,1)-fftbufs(i,bpstrdex,1))
   end do

   return
end subroutine linemsdyn_aft

!-----------------------------------------------------------------------

end module linemsdyn
