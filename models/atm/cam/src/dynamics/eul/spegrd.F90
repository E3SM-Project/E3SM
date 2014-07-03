
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Transfrom variables from spherical harmonic coefficients 
! to grid point values during second gaussian latitude scan (scan2)
! 
! Method: 
! Assemble northern and southern hemisphere grid values from the
! symmetric and antisymmetric fourier coefficients. 
! 1. Determine the fourier coefficients for the northern or southern
!    hemisphere latitude. 
! 2. Transform to gridpoint values
! 3. Clean up
!
! Author: 
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, J. Hack, August 1992
! Reviewed:          B. Boville, April 1996
! Modified:          P. Worley, September 2002
! 
!-----------------------------------------------------------------------
!

subroutine spegrd_bft (lat     ,nlon_fft, &
                       grdps   ,grzs    ,grds    ,gruhs   ,grvhs   , &
                       grths   ,grpss   ,grus    ,grvs    ,grts    , &
                       grpls   ,grpms   ,grdpa   ,grza    ,grda    , &
                       gruha   ,grvha   ,grtha   ,grpsa   ,grua    , &
                       grva    ,grta    ,grpla   ,grpma   ,fftbuf  )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Preparation for transform of variables from spherical harmonic 
! coefficients to grid point values during second gaussian latitude scan 
! (scan2)
! 
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, J. Hack, August 1992
! Reviewed:          B. Boville, April 1996
! Modified:          P. Worley, September 2002
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plat, plev, plevp
   use spmd_utils,   only:  iam
   use comspe, only: maxm, numm
!-----------------------------------------------------------------------
   implicit none
!---------------------------------------------------------------------
!
! Arguments
!
   integer, intent(in) :: lat                  ! latitude index
   integer, intent(in) :: nlon_fft             ! first dimension of FFT work array
!
! Symmetric fourier coefficient arrays for all variables transformed
! from spherical harmonics (see grcalc)
!
   real(r8), intent(in) :: grdps(2*maxm)       ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt*lnps(n,m)*P(n,m)
   real(r8), intent(in) :: grzs(2*maxm,plev)   ! sum(n) of z(n,m)*P(n,m) 
   real(r8), intent(in) :: grds(2*maxm,plev)   ! sum(n) of d(n,m)*P(n,m)
   real(r8), intent(in) :: gruhs(2*maxm,plev)  ! sum(n) of K(2i)*z(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(in) :: grvhs(2*maxm,plev)  ! sum(n) of K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(in) :: grths(2*maxm,plev)  ! sum(n) of K(2i)*t(n,m)*P(n,m)
   real(r8), intent(in) :: grpss(2*maxm)       ! sum(n) of lnps(n,m)*P(n,m)
   real(r8), intent(in) :: grus(2*maxm,plev)   ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(in) :: grvs(2*maxm,plev)   ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(in) :: grts(2*maxm,plev)   ! sum(n) of t(n,m)*P(n,m)
   real(r8), intent(in) :: grpls(2*maxm)       ! sum(n) of lnps(n,m)*P(n,m)*m/a
   real(r8), intent(in) :: grpms(2*maxm)       ! sum(n) of lnps(n,m)*H(n,m)
!
! Antisymmetric fourier coefficient arrays for all variables transformed
! from spherical harmonics (see grcalc)
!
   real(r8), intent(in) :: grdpa(2*maxm)       ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt*lnps(n,m)*P(n,m)
   real(r8), intent(in) :: grza(2*maxm,plev)   ! sum(n) of z(n,m)*P(n,m)
   real(r8), intent(in) :: grda(2*maxm,plev)   ! sum(n) of d(n,m)*P(n,m)
   real(r8), intent(in) :: gruha(2*maxm,plev)  ! sum(n)K(2i)*z(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(in) :: grvha(2*maxm,plev)  ! sum(n)K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(in) :: grtha(2*maxm,plev)  ! sum(n) of K(2i)*t(n,m)*P(n,m)
   real(r8), intent(in) :: grpsa(2*maxm)       ! sum(n) of lnps(n,m)*P(n,m)
   real(r8), intent(in) :: grua(2*maxm,plev)   ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(in) :: grva(2*maxm,plev)   ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(in) :: grta(2*maxm,plev)   ! sum(n) of t(n,m)*P(n,m)
   real(r8), intent(in) :: grpla(2*maxm)       ! sum(n) of lnps(n,m)*P(n,m)*m/a
   real(r8), intent(in) :: grpma(2*maxm)       ! sum(n) of lnps(n,m)*H(n,m)

   real(r8), intent(out) :: fftbuf(nlon_fft,8,plevp) ! buffer used for in-place FFTs

!
!---------------------------Local workspace-----------------------------
!
   integer i,k                           ! longitude, level indices
   integer rmlength                      ! twice number of local wavenumbers
   integer, parameter :: vortdex = 1     ! indices into fftbuf 
   integer, parameter :: divdex = 2
   integer, parameter :: duhdex = 3
   integer, parameter :: dvhdex = 4
   integer, parameter :: dthdex = 5
   integer, parameter :: u3dex = 6
   integer, parameter :: v3dex = 7
   integer, parameter :: t3dex = 8
   integer, parameter :: dpsdex = 1
   integer, parameter :: psdex = 2
   integer, parameter :: dpsldex = 3
   integer, parameter :: dpsmdex = 4
!
!-----------------------------------------------------------------------
!
! Assemble northern and southern hemisphere grid values from the
! symmetric and antisymmetric fourier coefficients: pre-FFT
!
!DIR$ NOSTREAM
   rmlength = 2*numm(iam)
   if (lat > plat/2) then                       ! Northern hemisphere
!DIR$ PREFERVECTOR
      do k=1,plev
!cdir loopchg
         do i=1,rmlength
            fftbuf(i,vortdex,k) = grzs(i,k) + grza(i,k)
            fftbuf(i,divdex,k) = grds(i,k) + grda(i,k)
            fftbuf(i,duhdex,k) = gruhs(i,k) + gruha(i,k)
            fftbuf(i,dvhdex,k) = grvhs(i,k) + grvha(i,k)
            fftbuf(i,dthdex,k) = grths(i,k) + grtha(i,k)
            fftbuf(i,u3dex,k) = grus(i,k) + grua(i,k)
            fftbuf(i,v3dex,k) = grvs(i,k) + grva(i,k)
            fftbuf(i,t3dex,k) = grts(i,k) + grta(i,k)
         end do
      end do
!
!cdir altcode=(loopcnt)
      do i=1,rmlength
         fftbuf(i,dpsdex,plevp) = grdps(i) + grdpa(i)
         fftbuf(i,psdex,plevp) = grpss(i) + grpsa(i)
         fftbuf(i,dpsldex,plevp) = grpls(i) + grpla(i)
         fftbuf(i,dpsmdex,plevp) = grpms(i) + grpma(i)
      end do

   else                                          ! Southern hemisphere

!DIR$ PREFERVECTOR
      do k=1,plev
!cdir loopchg
         do i=1,rmlength
            fftbuf(i,vortdex,k) = grzs(i,k) - grza(i,k)
            fftbuf(i,divdex,k) = grds(i,k) - grda(i,k)
            fftbuf(i,duhdex,k) = gruhs(i,k) - gruha(i,k)
            fftbuf(i,dvhdex,k) = grvhs(i,k) - grvha(i,k)
            fftbuf(i,dthdex,k) = grths(i,k) - grtha(i,k)
            fftbuf(i,u3dex,k) = grus(i,k) - grua(i,k)
            fftbuf(i,v3dex,k) = grvs(i,k) - grva(i,k)
            fftbuf(i,t3dex,k) = grts(i,k) - grta(i,k)
         end do
      end do

!cdir altcode=(loopcnt)
      do i=1,rmlength
         fftbuf(i,dpsdex,plevp) = grdps(i) - grdpa(i)
         fftbuf(i,psdex,plevp) = grpss(i) - grpsa(i)
         fftbuf(i,dpsldex,plevp) = grpls(i) - grpla(i)
         fftbuf(i,dpsmdex,plevp) = grpms(i) - grpma(i)
      end do

   end if

   return
end subroutine spegrd_bft

subroutine spegrd_ift (nlon_fft_in, nlon_fft_out, fftbuf_in, fftbuf_out)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Inverse Fourier transform of variables from spherical harmonic 
! coefficients to grid point values during second gaussian latitude scan 
! (scan2)
! 
! Author:             P. Worley, September 2002
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plat, plevp, beglat, endlat, plev
   use rgrid
   use comspe, only: maxm
#if ( defined SPMD )
   use mpishorthand
#endif
   use eul_control_mod, only : trig, ifax, pcray
   use perf_mod
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------

!---------------------------------------------------------------------
!
! Arguments
!
!     
! Input arguments
!     
   integer, intent(in) :: nlon_fft_in      ! first dimension of first FFT work array
   integer, intent(in) :: nlon_fft_out     ! first dimension of second FFT work array
#if (defined SPMD)
   real(r8), intent(in) :: fftbuf_in(nlon_fft_in,8,plevp,plat) 
                            ! buffer containing fields dcomposed over wavenumbers
#else
   real(r8), intent(in) :: fftbuf_in(1,1,1,1) 
                            ! buffer unused
#endif
!     
! Input/Output arguments
!     
   real(r8), intent(inout) :: fftbuf_out(nlon_fft_out,8,plevp,beglat:endlat) 
                            ! buffer used for in-place FFTs
!
!---------------------------Local workspace-----------------------------
!
#if ( ! defined USEFFTLIB )
   real(r8) work((plon+1)*8*plevp)
#else
   real(r8) work((plon+1)*pcray) ! workspace needed by fft991
#endif
   integer lat             ! latitude index
   integer isign           ! +1 => transform spectral to grid
   integer ntr             ! number of transforms to perform
   integer inc             ! distance between transform elements
   integer begtrm          ! (real) location of first truncated wavenumber
   integer k, ifld         ! level and field indices
!
!-----------------------------------------------------------------------
!
!
#if ( defined SPMD )
!
!  reorder Fourier coefficients
!
   call t_barrierf ('sync_realloc4b', mpicom)
   call t_startf('realloc4b')
   call realloc4b(nlon_fft_in, nlon_fft_out, fftbuf_in, fftbuf_out)
   call t_stopf('realloc4b')
#endif
!
! Zero elements corresponding to truncated wavenumbers, then
! transform from fourier coefficients to gridpoint values.
! ps,vort,div,duh,dvh,dth,dpsl,dpsm,dps,
! u,v,t (SLT) [If you want to do spectral transport, do q as well]
!
   begtrm = 2*pmmax+1
   inc = 1
   isign = +1
#ifdef OUTER_OMP
!$OMP PARALLEL DO PRIVATE (LAT, NTR, K, IFLD, WORK)
#endif
   do lat=beglat,endlat
      ntr = 8
!$OMP PARALLEL DO PRIVATE (K, WORK)
      do k=1,plev
         fftbuf_out(begtrm:nlon_fft_out,:,k,lat) = 0.0_r8
         call fft991 (fftbuf_out(1,1,k,lat), work, trig(1,lat), ifax(1,lat), inc, &
                      nlon_fft_out, nlon(lat), ntr, isign)
      enddo
      ntr = 1
!$OMP PARALLEL DO PRIVATE (IFLD, WORK)
      do ifld=1,4
         fftbuf_out(begtrm:nlon_fft_out,ifld,plevp,lat) = 0.0_r8
         call fft991 (fftbuf_out(1,ifld,plevp,lat), work, trig(1,lat), ifax(1,lat), inc, &
                      nlon_fft_out, nlon(lat), ntr, isign)
      enddo
   enddo
!
   return
end subroutine spegrd_ift

subroutine spegrd_aft (ztodt   ,lat     ,nlon    ,nlon_fft, &
                       cwava   ,qfcst   , &
                       etamid  ,ps      ,u3      ,v3      ,t3      , &
                       qminus  ,vort    ,div     ,hw2al   ,hw2bl   , &
                       hw3al   ,hw3bl   ,hwxal   ,hwxbl   ,q3m1    , &
                       dps     ,dpsl    ,dpsm    ,t3m2    ,engy2alat, &
                       engy2blat,difftalat, difftblat,phis,fftbuf  )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Completion of transformation of variables from spherical harmonic 
! coefficients to grid point values during second gaussian latitude scan 
! (scan2)
! 
! Method: 
! 
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, J. Hack, August 1992
! Reviewed:          B. Boville, April 1996
! Modified:          P. Worley, September 2002
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plat, plev, plevp
   use pspect
   use commap
   use cam_history,  only: outfld
   use physconst,    only: rga
   use constituents, only: pcnst
   use eul_control_mod
   use hycoef,       only: nprlev
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! Arguments
!
   integer, intent(in) :: lat                 ! latitude index
   integer, intent(in) :: nlon                ! number of longitudes
   integer, intent(in) :: nlon_fft            ! first dimension of FFT work arrays

   real(r8), intent(in) :: ztodt              ! twice the timestep unles nstep=0
   real(r8), intent(in) :: cwava              ! normalization factor (1/g*plon)
   real(r8), intent(in) :: qfcst(plon,plev,pcnst)
   real(r8), intent(in) :: qminus(plon,plev,pcnst)
   real(r8), intent(in) :: etamid(plev)       ! vertical coords at midpoints

   real(r8), intent(inout) :: ps(plon)  
   real(r8), intent(inout) :: u3(plon,plev)  
   real(r8), intent(inout) :: v3(plon,plev)  
   real(r8), intent(inout) :: t3(plon,plev)  
   real(r8), intent(inout) :: vort(plon,plev)
   real(r8), intent(inout) :: div(plon,plev)
   real(r8), intent(inout) :: q3m1(plon,plev,pcnst)

   real(r8), intent(out) :: hw2al(pcnst)  ! -
   real(r8), intent(out) :: hw2bl(pcnst)  !  | lat contributions to components
   real(r8), intent(out) :: hw3al(pcnst)  !  | of slt global mass integrals 
   real(r8), intent(out) :: hw3bl(pcnst)  ! -
   real(r8), intent(out) :: hwxal(pcnst,4)
   real(r8), intent(out) :: hwxbl(pcnst,4)

   real(r8), intent(out) :: dps(plon)
   real(r8), intent(out) :: dpsl(plon)
   real(r8), intent(out) :: dpsm(plon)
   real(r8), intent(in)  :: t3m2(plon,plev) ! temperature
   real(r8), intent(out) :: engy2alat
   real(r8), intent(out) :: engy2blat
   real(r8), intent(out) :: difftalat
   real(r8), intent(out) :: difftblat
   real(r8), intent(in)  :: phis(plon)
   real(r8), intent(in) :: fftbuf(nlon_fft,8,plevp) ! buffer used for in-place FFTs
!
!---------------------------Local workspace-----------------------------
!
   real(r8) :: duh(plon,plev)  ! 
   real(r8) :: dvh(plon,plev)  ! 
   real(r8) :: dth(plon,plev)  ! 

   real(r8) pmid(plon,plev)    ! pressure at model levels
   real(r8) pint(plon,plevp)   ! pressure at model interfaces
   real(r8) pdel(plon,plev)    ! pdel(k) = pint(k+1) - pint(k)
   real(r8) pdelb(plon,plev)   ! pressure diff bet intfcs (press defined using the "B" part 
                               ! of the hybrid grid only)
   real(r8) hcwavaw            ! 0.5*cwava*w(lat)
   real(r8) sum
!
   real(r8) rcoslat            ! 1./cosine(latitude)
   real(r8) dotproda           ! dot product
   real(r8) dotprodb           ! dot product
   integer i,k,m               ! longitude, level, constituent indices
   integer klev                ! top level where hybrid coordinates apply
   integer, parameter :: vortdex = 1     ! indices into fftbuf 
   integer, parameter :: divdex = 2
   integer, parameter :: duhdex = 3
   integer, parameter :: dvhdex = 4
   integer, parameter :: dthdex = 5
   integer, parameter :: u3dex = 6
   integer, parameter :: v3dex = 7
   integer, parameter :: t3dex = 8
   integer, parameter :: dpsdex = 1
   integer, parameter :: psdex = 2
   integer, parameter :: dpsldex = 3
   integer, parameter :: dpsmdex = 4
!
!-----------------------------------------------------------------------
!
! Copy 3D fields out of FFT buffer, removing cosine(latitude) from momentum variables
!
   rcoslat = 1._r8/cos(clat(lat))
!$OMP PARALLEL DO PRIVATE (K, I)
   do k=1,plev
      do i=1,nlon
         vort(i,k) = fftbuf(i,vortdex,k)
         div(i,k) = fftbuf(i,divdex,k)
         duh(i,k) = fftbuf(i,duhdex,k)*rcoslat
         dvh(i,k) = fftbuf(i,dvhdex,k)*rcoslat
         dth(i,k) = fftbuf(i,dthdex,k)
         u3(i,k) = fftbuf(i,u3dex,k)*rcoslat
         v3(i,k) = fftbuf(i,v3dex,k)*rcoslat
         t3(i,k) = fftbuf(i,t3dex,k)
      end do
   end do
!
! Copy 2D fields out of FFT buffer, converting
! log(ps) to ps.
!
!$OMP PARALLEL DO PRIVATE (I)
   do i=1,nlon
      dps(i) = fftbuf(i,dpsdex,plevp)
      dpsl(i) = fftbuf(i,dpsldex,plevp)
      dpsm(i) = fftbuf(i,dpsmdex,plevp)
      ps(i) = exp(fftbuf(i,psdex,plevp))
   end do

!
! Diagnose pressure arrays needed by DIFCOR
!
   call plevs0 (nlon, plon, plev, ps, pint, pmid, pdel)
   call pdelb0 (ps, pdelb, nlon)
!
! Accumulate mass integrals
!
   sum = 0._r8
   do i=1,nlon
      sum = sum + ps(i)
   end do
   tmass(lat) = w(lat)*rga*sum/nlon
!
! Finish horizontal diffusion: add pressure surface correction term to t and
! q diffusions; add kinetic energy dissipation to internal energy (temperature)
!
   klev = max(kmnhd4,nprlev)
   call difcor (klev,        ztodt,  dps,    u3,     v3, &
                q3m1(1,1,1), pdel,   pint,   t3,     dth, &
                duh,         dvh,    nlon)
!
! Calculate SLT moisture, constituent, energy, and temperature integrals
!
   hcwavaw   = 0.5_r8*cwava*w(lat)
   engy2alat = 0._r8
   engy2blat = 0._r8
   difftalat = 0._r8
   difftblat = 0._r8
!$OMP PARALLEL DO PRIVATE (M, K, DOTPRODA, DOTPRODB, I)
   do m=1,pcnst
      hw2al(m) = 0._r8
      hw2bl(m) = 0._r8
      hw3al(m) = 0._r8
      hw3bl(m) = 0._r8
      hwxal(m,1) = 0._r8
      hwxal(m,2) = 0._r8
      hwxal(m,3) = 0._r8
      hwxal(m,4) = 0._r8
      hwxbl(m,1) = 0._r8
      hwxbl(m,2) = 0._r8
      hwxbl(m,3) = 0._r8
      hwxbl(m,4) = 0._r8
      do k=1,plev
         dotproda = 0._r8
         dotprodb = 0._r8
         do i=1,nlon
            dotproda = dotproda + qfcst(i,k,m)*pdela(i,k)
            dotprodb = dotprodb + qfcst(i,k,m)*pdelb(i,k)
         end do
         hw2al(m) = hw2al(m) + hcwavaw*dotproda
         hw2bl(m) = hw2bl(m) + hcwavaw*dotprodb
      end do
   end do

! using do loop and select to enable functional parallelism with OpenMP
!$OMP PARALLEL DO PRIVATE (I)
   do i=1,6
      select case (i)
      case (1)
         call engy_te  (cwava ,w(lat) ,t3  ,u3  ,v3 ,phis    ,pdela, engy2alat ,nlon)
      case (2)
         call engy_te  (cwava ,w(lat) ,t3  ,u3  ,v3 ,phis    ,pdelb, engy2blat ,nlon)
      case (3)
         call engy_tdif(cwava ,w(lat) ,t3  ,t3m2             ,pdela, difftalat ,nlon)
      case (4)
         call engy_tdif(cwava ,w(lat) ,t3  ,t3m2             ,pdelb, difftblat ,nlon)
      case (5)
         call qmassd (cwava, etamid, w(lat), qminus, qfcst, &
                      pdela, hw3al, nlon)
      case (6)
         call qmassd (cwava, etamid, w(lat), qminus, qfcst, &
                      pdelb, hw3bl, nlon)
      end select
   end do

   if (pcnst.gt.1) then
      call xqmass (cwava, etamid, w(lat), qminus, qfcst, &
                   qminus, qfcst, pdela, pdelb, hwxal, &
                   hwxbl, nlon)
   end if

   call outfld ('DTH     ',dth     ,plon   ,lat     )

   return
end subroutine spegrd_aft


