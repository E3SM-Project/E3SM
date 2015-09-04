
module spegrd
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
! Original version:  J. Rosinski
! Modified:          P. Worley, October 2002
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plev, plat, plevp
   use perf_mod

   implicit none

   private
!
! Public interfaces
!
   public spegrd_bft   ! Before FFT
   public spegrd_ift   ! FFT
   public spegrd_aft   ! After FFT

contains

!
!-----------------------------------------------------------------------
!

subroutine spegrd_bft (lat     ,nlon_fft, &
                       grts    ,grqs    ,grths   , &
                       grds    ,grus    ,gruhs   ,grvs    ,grvhs   , &
                       grpss   ,grdps   ,grpms   ,grpls   ,grtms   , &
                       grtls   ,grqms   ,grqls   ,grta    ,grqa    , &
                       grtha   ,grda    ,grua    ,gruha   ,grva    , &
                       grvha   ,grpsa   ,grdpa   ,grpma   ,grpla   , &
                       grtma   ,grtla   ,grqma   ,grqla   ,fftbuf   )
!-----------------------------------------------------------------------
!
! Purpose:
! Preparation for transform of variables from spherical harmonic 
! coefficients to grid point values during second gaussian latitude scan 
! (scan2)
!
! Method: 
! 
! Original version:  J. Rosinski
! Modified:          P. Worley, October 2002
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

   use pmgrid
   use comspe, only: maxm, numm
   use spmd_utils, only : iam
!
! Arguments
!
   integer, intent(in) :: lat                   ! latitude index
   integer, intent(in) :: nlon_fft             ! first dimension of FFT work array
!
! Symmetric fourier coefficient arrays for all variables transformed 
! from spherical harmonics (see subroutine grcalc)
!                                
   real(r8), intent(in) :: grdps(2*maxm)         ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt*lnps(n,m)*P(n,m)
   real(r8), intent(in) :: grds (2*maxm,plev)    ! sum(n) of d(n,m)*P(n,m)
   real(r8), intent(in) :: gruhs(2*maxm,plev)    ! sum(n) of K(2i)*z(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(in) :: grvhs(2*maxm,plev)    ! sum(n) of K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(in) :: grths(2*maxm,plev)    ! sum(n) of K(2i)*t(n,m)*P(n,m)
   real(r8), intent(in) :: grpss(2*maxm)         ! sum(n) of lnps(n,m)*P(n,m)
   real(r8), intent(in) :: grus (2*maxm,plev)    ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(in) :: grvs (2*maxm,plev)    ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(in) :: grts (2*maxm,plev)    ! sum(n) of t(n,m)*P(n,m)
   real(r8), intent(in) :: grqs (2*maxm,plev)    ! sum(n) of q(n,m)*P(n,m)
   real(r8), intent(in) :: grpls(2*maxm)         ! sum(n) of lnps(n,m)*P(n,m)*m/a
   real(r8), intent(in) :: grtms(2*maxm,plev)
   real(r8), intent(in) :: grtls(2*maxm,plev)
   real(r8), intent(in) :: grqms(2*maxm,plev)
   real(r8), intent(in) :: grqls(2*maxm,plev)
   real(r8), intent(in) :: grpms(2*maxm)         ! sum(n) of lnps(n,m)*H(n,m)
!
! Antisymmetric fourier coefficient arrays for all variables
! transformed from spherical harmonics (see grcalc)
!
   real(r8), intent(in) :: grdpa(2*maxm)         ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt*lnps(n,m)*P(n,m)
   real(r8), intent(in) :: grda (2*maxm,plev)    ! sum(n) of d(n,m)*P(n,m)
   real(r8), intent(in) :: gruha(2*maxm,plev)    ! sum(n)K(2i)*z(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(in) :: grvha(2*maxm,plev)    ! sum(n)K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(in) :: grtha(2*maxm,plev)    ! sum(n) of K(2i)*t(n,m)*P(n,m)
   real(r8), intent(in) :: grpsa(2*maxm)         ! sum(n) of lnps(n,m)*P(n,m)
   real(r8), intent(in) :: grua (2*maxm,plev)    ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(in) :: grva (2*maxm,plev)    ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(in) :: grta (2*maxm,plev)    ! sum(n) of t(n,m)*P(n,m)
   real(r8), intent(in) :: grqa (2*maxm,plev)    ! sum(n) of q(n,m)*P(n,m)
   real(r8), intent(in) :: grpla(2*maxm)         ! sum(n) of lnps(n,m)*P(n,m)*m/a
   real(r8), intent(in) :: grtma(2*maxm,plev)
   real(r8), intent(in) :: grtla(2*maxm,plev)
   real(r8), intent(in) :: grqma(2*maxm,plev)
   real(r8), intent(in) :: grqla(2*maxm,plev)
   real(r8), intent(in) :: grpma(2*maxm)         ! sum(n) of lnps(n,m)*H(n,m)
!
   real(r8), intent(out) :: fftbuf(nlon_fft,11,plevp) ! buffer used for in-place FFTs

!
!---------------------------Local workspace-----------------------------
!
!
! Local workspace
!
   integer i,k                           ! longitude, level
   integer rmlength                      ! twice number of local wavenumbers
   integer, parameter :: divdex = 1      ! indices into fftbuf 
   integer, parameter :: duhdex = 2
   integer, parameter :: dvhdex = 3
   integer, parameter :: dthdex = 4
   integer, parameter :: tldex = 5
   integer, parameter :: tmdex = 6
   integer, parameter :: qldex = 7
   integer, parameter :: qmdex = 8
   integer, parameter :: u3dex = 9
   integer, parameter :: v3dex = 10
   integer, parameter :: t3dex = 11
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
   rmlength = 2*numm(iam)
   if (lat > plat/2) then                       ! Northern hemisphere
      do k=1,plev
         do i=1,rmlength
            fftbuf(i,divdex,k) = grds(i,k) + grda(i,k)
            fftbuf(i,duhdex,k) = gruhs(i,k) + gruha(i,k)
            fftbuf(i,dvhdex,k) = grvhs(i,k) + grvha(i,k)
            fftbuf(i,dthdex,k) = grths(i,k) + grtha(i,k)
            fftbuf(i,tldex,k)  = grtls(i,k) + grtla(i,k)
            fftbuf(i,tmdex,k)  = grtms(i,k) + grtma(i,k)
            fftbuf(i,qldex,k)  = grqls(i,k) + grqla(i,k)
            fftbuf(i,qmdex,k)  = grqms(i,k) + grqma(i,k)
            fftbuf(i,u3dex,k)  = grus(i,k) + grua(i,k)
            fftbuf(i,v3dex,k)  = grvs(i,k) + grva(i,k)
            fftbuf(i,t3dex,k)  = grts(i,k) + grta(i,k)
         end do
      end do

      do i=1,rmlength
         fftbuf(i,dpsdex,plevp)  = grdps(i) + grdpa(i)
         fftbuf(i,psdex,plevp)   = grpss(i) + grpsa(i)
         fftbuf(i,dpsldex,plevp) = grpls(i) + grpla(i)
         fftbuf(i,dpsmdex,plevp) = grpms(i) + grpma(i)
      end do

   else                                          ! Southern hemisphere

      do k=1,plev
         do i=1,rmlength
            fftbuf(i,divdex,k) = grds(i,k) - grda(i,k)
            fftbuf(i,duhdex,k) = gruhs(i,k) - gruha(i,k)
            fftbuf(i,dvhdex,k) = grvhs(i,k) - grvha(i,k)
            fftbuf(i,dthdex,k) = grths(i,k) - grtha(i,k)
            fftbuf(i,tldex,k)  = grtls(i,k) - grtla(i,k)
            fftbuf(i,tmdex,k)  = grtms(i,k) - grtma(i,k)
            fftbuf(i,qldex,k)  = grqls(i,k) - grqla(i,k)
            fftbuf(i,qmdex,k)  = grqms(i,k) - grqma(i,k)
            fftbuf(i,u3dex,k)  = grus(i,k) - grua(i,k)
            fftbuf(i,v3dex,k)  = grvs(i,k) - grva(i,k)
            fftbuf(i,t3dex,k)  = grts(i,k) - grta(i,k)
         end do
      end do

      do i=1,rmlength
         fftbuf(i,dpsdex,plevp)  = grdps(i) - grdpa(i)
         fftbuf(i,psdex,plevp)   = grpss(i) - grpsa(i)
         fftbuf(i,dpsldex,plevp) = grpls(i) - grpla(i)
         fftbuf(i,dpsmdex,plevp) = grpms(i) - grpma(i)
      end do
   end if
!
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
! Method: 
! 
! Original version:  J. Rosinski
! Modified:          P. Worley, October 2002
!
!-----------------------------------------------------------------------
   use rgrid,  only: nlon, pmmax
   use pmgrid, only: plat, plon, beglat, endlat
   use comspe, only: maxm
#if ( defined SPMD )
   use mpishorthand
#endif
   use sld_control_mod, only : trig, ifax, pcray
!-----------------------------------------------------------------------
!---------------------------------------------------------------------
!
! Arguments
!
   integer, intent(in) :: nlon_fft_in      ! first dimension of first FFT work array
   integer, intent(in) :: nlon_fft_out     ! first dimension of second FFT work array
#if (defined SPMD)
   real(r8), intent(in) :: fftbuf_in(nlon_fft_in,11,plevp,plat) 
                            ! buffer containing fields dcomposed over wavenumbers
#else
   real(r8), intent(in) :: fftbuf_in(1,1,1,1) 
                            ! buffer unused
#endif
!     
! Input/Output arguments
!     
   real(r8), intent(inout) :: fftbuf_out(nlon_fft_out,11,plevp,beglat:endlat) 
                            ! buffer used for in-place FFTs
!
!---------------------------Local workspace-----------------------------
!
#if ( ! defined USEFFTLIB )
   real(r8) work((plon+1)*11*plevp)
#else
   real(r8) work((plon+1)*pcray) ! workspace needed by fft991
#endif
   integer lat               ! latitude index
   integer isign           ! +1 => transform spectral to grid
   integer ntr             ! number of transforms to perform
   integer inc             ! distance between transform elements
   integer begtrm          ! (real) location of first truncated wavenumber
!
!-----------------------------------------------------------------------
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
! transform from fourier coefficients to gridpoint values:
! ps,div,tl,tm,dpsl,dpsm,ql,qm,dth,duh,dvh,dps,u,v,t
!
   begtrm = 2*pmmax+1
   inc = 1
   isign = +1
   ntr = 11*plev + 4
   fftbuf_out(begtrm:nlon_fft_out,:,:,:) = 0.0_r8
!$OMP PARALLEL DO PRIVATE (LAT, WORK)
   do lat=beglat,endlat
      call fft991 (fftbuf_out(1,1,1,lat), work, trig(1,lat), ifax(1,lat), inc, &
                   nlon_fft_out, nlon(lat), ntr, isign)
   enddo
!
   return
end subroutine spegrd_ift

subroutine spegrd_aft (ztodt   ,lat     ,nlon    ,nlon_fft, &
                       cwava   ,qfcst   ,q3      , &
                       etamid  ,ps      ,u3      ,v3      ,t3      , &
                       div     ,hw2al   ,hw2bl   ,hw3al   ,hw3bl   , &
                       hwxal   ,hwxbl   ,dps     , &
                       dpsl    ,dpsm    ,tl      ,tm      ,ql      , &
                       qm      ,t3m1    ,engy2alat,engy2blat,difftalat, &
                       difftblat,phis   ,fftbuf  )
!-----------------------------------------------------------------------
!
! Purpose:
! Completion of transformation of variables from spherical harmonic 
! coefficients to grid point values during second gaussian latitude scan 
! (scan2)
!
! Method: 
! 
! Author: 
! 
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------
   use constituents,    only: pcnst
   use pspect
   use comspe
   use commap
   use cam_history,     only: outfld
   use physconst,       only: rga
   use hycoef,          only: nprlev
   use sld_control_mod, only: tmass, pdela, kmnhd4
!
! Arguments
!
   integer, intent(in) :: lat                     ! latitude index
   integer, intent(in) :: nlon                    ! number of longitudes
   integer, intent(in) :: nlon_fft                ! first dimension of FFT work arrays
!
   real(r8), intent(in)   :: ztodt                ! timestep
   real(r8), intent(in)   :: cwava                ! normalization factor (1/g*plon)
   real(r8), intent(in)   :: qfcst(plon,plev,pcnst)       ! fcst q + consts
   real(r8), intent(in)   :: q3(plon,plev,pcnst)  ! q + consts
   real(r8), intent(in)   :: etamid(plev)                     ! vertical coords at midpts 
   real(r8), intent(inout)  :: ps(plon)       ! surface pressure
   real(r8), intent(inout)  :: u3(plon,plev)  ! u-wind
   real(r8), intent(inout)  :: v3(plon,plev)  ! v-wind
   real(r8), intent(inout)  :: t3(plon,plev)  ! temperature
   real(r8), intent(inout) :: div(plon,plev)  ! divergence

   real(r8), intent(out)  :: hw2al(pcnst)               ! -
   real(r8), intent(out)  :: hw2bl(pcnst)               !  | lat contributions to
   real(r8), intent(out)  :: hw3al(pcnst)               !  | components of slt global
   real(r8), intent(out)  :: hw3bl(pcnst)               !  | mass integrals
   real(r8), intent(out)  :: hwxal(pcnst,4)             !  |
   real(r8), intent(out)  :: hwxbl(pcnst,4)             ! -

   real(r8), intent(out) :: dps(plon)
   real(r8), intent(out) :: dpsl(plon)
   real(r8), intent(out) :: dpsm(plon)

   real(r8), intent(out) :: tl(plon,plev)
   real(r8), intent(out) :: tm(plon,plev)
   real(r8), intent(out) :: ql(plon,plev)
   real(r8), intent(out) :: qm(plon,plev)
   real(r8), intent(in)  :: t3m1(plon,plev) ! temperature
   real(r8), intent(out) :: engy2alat
   real(r8), intent(out) :: engy2blat
   real(r8), intent(out) :: difftalat
   real(r8), intent(out) :: difftblat
   real(r8), intent(in)  :: phis(plon)
!
   real(r8), intent(in) :: fftbuf(nlon_fft,11,plevp) ! buffer used for in-place FFTs
!
!---------------------------Local workspace-----------------------------
!
   real(r8) :: duh(plon,plev) ! 
   real(r8) :: dvh(plon,plev) ! 
   real(r8) :: dth(plon,plev) ! 

   real(r8) pmid (plon,plev)     ! pressure at model levels
   real(r8) pint (plon,plevp)    ! pressure at model interfaces
   real(r8) pdel (plon,plev)     ! pdel(k) = pint(k+1) - pint(k)
   real(r8) pdelb(plon,plev)     ! pressure diff between interfaces
!                                ! (press defined using the "B" part 
!                                ! of the hybrid grid only)
   real(r8) hcwavaw              ! 0.5*cwava*w(irow)
   real(r8) sum

   real(r8) rcoslat              ! 1./cosine(latitude)
   real(r8) dotproda             ! dot product
   real(r8) dotprodb             ! dot product

   integer i,k,m                 ! longitude, level, constituent indices
   integer ihem                  ! hemisphere index
   integer klev                  ! top level where hybrid coordinates apply
   
   integer, parameter :: divdex = 1      ! indices into fftbuf 
   integer, parameter :: duhdex = 2
   integer, parameter :: dvhdex = 3
   integer, parameter :: dthdex = 4
   integer, parameter :: tldex = 5
   integer, parameter :: tmdex = 6
   integer, parameter :: qldex = 7
   integer, parameter :: qmdex = 8
   integer, parameter :: u3dex = 9
   integer, parameter :: v3dex = 10
   integer, parameter :: t3dex = 11
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
   do k=1,plev
      do i=1,nlon
         div(i,k) = fftbuf(i,divdex,k)
         duh(i,k) = fftbuf(i,duhdex,k)*rcoslat
         dvh(i,k) = fftbuf(i,dvhdex,k)*rcoslat
         dth(i,k) = fftbuf(i,dthdex,k)
         tl (i,k) = fftbuf(i,tldex,k)
         tm (i,k) = fftbuf(i,tmdex,k)
         ql (i,k) = fftbuf(i,qldex,k)
         qm (i,k) = fftbuf(i,qmdex,k)
         u3(i,k)  = fftbuf(i,u3dex,k)*rcoslat
         v3(i,k)  = fftbuf(i,v3dex,k)*rcoslat
         t3(i,k)  = fftbuf(i,t3dex,k)
      end do
   end do
!
! Copy 2D fields out of FFT buffer, converting
! log(ps) to ps.
!
   do i=1,nlon
      dps(i)  = fftbuf(i,dpsdex,plevp)
      ps(i)   = exp(fftbuf(i,psdex,plevp))
      dpsl(i) = fftbuf(i,dpsldex,plevp)
      dpsm(i) = fftbuf(i,dpsmdex,plevp)
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
! Finish horizontal diffusion: add pressure surface correction term to
! t and q diffusions; add kinetic energy dissipation to internal energy
! (temperature)
!
   klev = max(kmnhd4,nprlev)

   call difcor (klev,      ztodt,  dps,    u3,     v3, &
                q3(1,1,1), pdel,   pint,   t3,     dth, &
                duh,       dvh,    nlon)
!
! Calculate SLT moisture, constituent, energy, and temperature integrals
!
   hcwavaw   = 0.5_r8*cwava*w(lat)
   engy2alat = 0._r8
   engy2blat = 0._r8
   difftalat = 0._r8
   difftblat = 0._r8
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

   call engy_te  (cwava ,w(lat) ,t3  ,u3  ,v3 ,phis    ,pdela, engy2alat ,nlon)
   call engy_te  (cwava ,w(lat) ,t3  ,u3  ,v3 ,phis    ,pdelb, engy2blat ,nlon)
   call engy_tdif(cwava ,w(lat) ,t3  ,t3m1             ,pdela, difftalat ,nlon)
   call engy_tdif(cwava ,w(lat) ,t3  ,t3m1             ,pdelb, difftblat ,nlon)

   call qmassd (cwava, etamid, w(lat), q3(1,1,1), qfcst(1,1,1), &
                pdela, hw3al, nlon)

   call qmassd (cwava, etamid, w(lat), q3(1,1,1), qfcst(1,1,1), &
                pdelb, hw3bl, nlon)

   if (pcnst.gt.1) then
      call xqmass (cwava, etamid, w(lat), q3(1,1,1), qfcst(1,1,1), &
                   q3(1,1,1), qfcst(1,1,1), pdela, pdelb, hwxal, &
                   hwxbl, nlon)
   end if

   call outfld ('DTH     ',dth     ,plon   ,lat     )

   return
end subroutine spegrd_aft

end module spegrd
