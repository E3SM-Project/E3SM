!-----------------------------------------------------------------------
module scan2
!----------------------------------------------------------------------- 
! 
! Purpose: Module for second gaussian latitude scan, to convert from
! spectral coefficients to grid point values.
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plat, plev, plon, beglat, endlat, plevp
   use constituents, only: pcnst
   use perf_mod
!-----------------------------------------------------------------------
   implicit none
!
! By default everything is private to this module
!
   private
!
! Public interfaces
!
   public scan2run     ! Public run method

!
! Private module data
!
   integer, parameter :: plondfft = plon + 2

!----------------------------------------------------------------------- 
contains
!----------------------------------------------------------------------- 

!
!----------------------------------------------------------------------- 
!

subroutine scan2run (ztodt, cwava, etamid,t2      ,fu      ,fv    )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Second gaussian latitude scan, converts from spectral coefficients to 
! grid point values, from poles to equator, with read/calculate/write cycle.
! 
! Method: 
! The latitude pair loop in this routine is multitasked.
!
! The grid point values of ps, t, u, v, z (vorticity), and d (divergence)
! are calculated and stored for each latitude from the spectral coefficients.
! In addition, the pressure-surface corrections to the horizontal diffusion
! are applied and the global integrals of the constituent fields are 
! computed for the mass fixer.
!
! Author: 
! Original version:  CCM1
!
!-----------------------------------------------------------------------
   use prognostics,  only: ps, u3, v3, q3, t3, dps, dpsl, dpsm, vort, &
                           qminus, div, n3, n3m1, n3m2, phis, omga,   &
                           shift_time_indices, hadv, pdeld
   use comspe,       only: maxm
   use rgrid,        only: nlon
   use scanslt,      only: hw1lat, engy1lat, qfcst
#ifdef SPMD
   use mpishorthand, only: mpicom, mpir8
#endif
   use physconst,    only: cpair
   use scamMod,      only: fixmascam,alphacam,betacam,use_iop, single_column
   use pspect,       only: pnmax
   use tfilt_massfix, only: tfilt_massfixrun
   use massfix,      only: hw1,hw2,hw3,alpha
   use cam_control_mod, only: ideal_phys, adiabatic
   use eul_control_mod, only: qmassf, tmass, tmass0, fixmas, tmassf

!-----------------------------------------------------------------------
!
! Input arguments
!
   real(r8), intent(in) :: ztodt                ! twice the timestep unless nstep = 0
   real(r8), intent(in) :: cwava(plat)          ! weight applied to global integrals
   real(r8), intent(in) :: etamid(plev)         ! vertical coords at midpoints 
   real(r8), optional, intent(inout) :: t2(plon,plev,beglat:endlat)         ! tot dT/dt to to physics
   real(r8), optional, intent(inout) :: fu(plon,plev,beglat:endlat)         ! u wind tend
   real(r8), optional, intent(inout) :: fv(plon,plev,beglat:endlat)         ! v wind tend
!
!---------------------------Local workspace-----------------------------
!
   real(r8) engy1         ! component of global energy integral (for time step n)
   real(r8) engy2         ! component of global energy integral (for time step n+1)
   real(r8) engy2a        ! component of global energy integral (for time step n+1)
   real(r8) engy2b        ! component of global energy integral (for time step n+1)
   real(r8) difft         ! component of global delta-temp integral ( (n+1) - n )
   real(r8) diffta        ! component of global delta-temp integral ( (n+1) - n )
   real(r8) difftb        ! component of global delta-temp integral ( (n+1) - n )
   real(r8) hw2a(pcnst)   ! component of constituent global mass integral (mass weighting is 
                          ! based upon the "A" portion of the hybrid grid)
   real(r8) hw2b(pcnst)   ! component of constituent global mass integral (mass weighting is 
                          ! based upon the "B" portion of the hybrid grid)
   real(r8) hw3a(pcnst)   ! component of constituent global mass integral (mass weighting is 
                          ! based upon the "A" portion of the hybrid grid)
   real(r8) hw3b(pcnst)   ! component of constituent global mass integral (mass weighting is 
                          ! based upon the "B" portion of the hybrid grid)
   real(r8) hwxa(pcnst,4)
   real(r8) hwxb(pcnst,4)
   real(r8) engy2alat(plat)     ! lat contribution to total energy integral
   real(r8) engy2blat(plat)     ! lat contribution to total energy integral
   real(r8) difftalat(plat)     ! lat contribution to delta-temperature integral
   real(r8) difftblat(plat)     ! lat contribution to delta-temperature integral
   real(r8) hw2al(pcnst,plat)   ! |------------------------------------
   real(r8) hw2bl(pcnst,plat)   ! | latitudinal contributions to the
   real(r8) hw3al(pcnst,plat)   ! | components of global mass integrals
   real(r8) hw3bl(pcnst,plat)   ! |
   real(r8) hwxal(pcnst,4,plat) ! |
   real(r8) hwxbl(pcnst,4,plat) ! |-----------------------------------
!
! Symmetric fourier coefficient arrays for all variables transformed 
! from spherical harmonics (see subroutine grcalc)
!                                
   real(r8) grdpss(2*maxm,(plat+1)/2)       ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt*lnps(n,m)*P(n,m)
   real(r8) grzs(2*maxm,plev,(plat+1)/2)   ! sum(n) of z(n,m)*P(n,m) 
   real(r8) grds(2*maxm,plev,(plat+1)/2)   ! sum(n) of d(n,m)*P(n,m)
   real(r8) gruhs(2*maxm,plev,(plat+1)/2)  ! sum(n) of K(2i)*z(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grvhs(2*maxm,plev,(plat+1)/2)  ! sum(n) of K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grths(2*maxm,plev,(plat+1)/2)  ! sum(n) of K(2i)*t(n,m)*P(n,m)
   real(r8) grpss(2*maxm,(plat+1)/2)       ! sum(n) of lnps(n,m)*P(n,m)
   real(r8) grus(2*maxm,plev,(plat+1)/2)   ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grvs(2*maxm,plev,(plat+1)/2)   ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grts(2*maxm,plev,(plat+1)/2)   ! sum(n) of t(n,m)*P(n,m)
   real(r8) grpls(2*maxm,(plat+1)/2)       ! sum(n) of lnps(n,m)*P(n,m)*m/a
   real(r8) grpms(2*maxm,(plat+1)/2)       ! sum(n) of lnps(n,m)*H(n,m)
!
! Antisymmetric fourier coefficient arrays for all variables transformed
! from spherical harmonics (see grcalc)
!
   real(r8) grdpsa(2*maxm,(plat+1)/2)       ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt*lnps(n,m)*P(n,m)
   real(r8) grza(2*maxm,plev,(plat+1)/2)   ! sum(n) of z(n,m)*P(n,m)
   real(r8) grda(2*maxm,plev,(plat+1)/2)   ! sum(n) of d(n,m)*P(n,m)
   real(r8) gruha(2*maxm,plev,(plat+1)/2)  ! sum(n)K(2i)*z(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grvha(2*maxm,plev,(plat+1)/2)  ! sum(n)K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grtha(2*maxm,plev,(plat+1)/2)  ! sum(n) of K(2i)*t(n,m)*P(n,m)
   real(r8) grpsa(2*maxm,(plat+1)/2)       ! sum(n) of lnps(n,m)*P(n,m)
   real(r8) grua(2*maxm,plev,(plat+1)/2)   ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grva(2*maxm,plev,(plat+1)/2)   ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
   real(r8) grta(2*maxm,plev,(plat+1)/2)   ! sum(n) of t(n,m)*P(n,m)
   real(r8) grpla(2*maxm,(plat+1)/2)       ! sum(n) of lnps(n,m)*P(n,m)*m/a
   real(r8) grpma(2*maxm,(plat+1)/2)       ! sum(n) of lnps(n,m)*H(n,m)
   real(r8) residual                   ! residual energy integral
   real(r8) beta                       ! energy fixer coefficient
!
   integer m,n                         ! indices
   integer lat,j,irow                  ! latitude indices
   integer nlon_fft_in                 ! FFT work array inner dimension
   integer nlon_fft_out                ! FFT work array inner dimension
!
! FFT buffers
!
   real(r8), allocatable:: fftbuf_in(:,:,:,:)  ! fftbuf_in(nlon_fft_in,8,plevp,plat) 
   real(r8), allocatable:: fftbuf_out(:,:,:,:) ! fftbuf_out(nlon_fft_out,8,plevp,beglat:endlat) 
!
! Temporal space for rearranged spectral coeffs. The rearrangement will
! be made in prepGRcalc and the rearranged coeffs will be transformed
! to Fourier coeffs in grcalca and grcalcs.
! 
   real(r8) tmpSPEcoef(plev*24,pnmax,maxm) 

!
!-----------------------------------------------------------------------
   if (.not. single_column) then

      call t_startf ('grcalc')

      call prepGRcalc(tmpSPEcoef)

#if ( defined SPMD )

!$OMP PARALLEL DO PRIVATE (J)
      do j=1,plat/2
         call grcalcs (j, ztodt, grts(1,1,j), grths(1,1,j), grds(1,1,j), &
              grzs(1,1,j), grus(1,1,j), gruhs(1,1,j), grvs(1,1,j), grvhs(1,1,j), &
              grpss(1,j), grdpss(1,j), grpms(1,j), grpls(1,j), tmpSPEcoef)
         
         call grcalca (j, ztodt, grta(1,1,j), grtha(1,1,j), grda(1,1,j), &
              grza(1,1,j), grua(1,1,j), gruha(1,1,j), grva(1,1,j), grvha(1,1,j), &
              grpsa(1,j), grdpsa(1,j), grpma(1,j), grpla(1,j), tmpSPEcoef)
      end do

#else

!$OMP PARALLEL DO PRIVATE (LAT, J)
   do lat=beglat,endlat
      if (lat > plat/2) then
         j = plat - lat + 1
         call grcalcs (j, ztodt, grts(1,1,j), grths(1,1,j), grds(1,1,j), &
                       grzs(1,1,j), grus(1,1,j), gruhs(1,1,j), grvs(1,1,j), grvhs(1,1,j), &
                       grpss(1,j), grdpss(1,j), grpms(1,j), grpls(1,j), tmpSPEcoef)
      else
         j = lat
         call grcalca (j, ztodt, grta(1,1,j), grtha(1,1,j), grda(1,1,j), &
                       grza(1,1,j), grua(1,1,j), gruha(1,1,j), grva(1,1,j), grvha(1,1,j), &
                       grpsa(1,j), grdpsa(1,j), grpma(1,j), grpla(1,j), tmpSPEcoef)
      end if
   end do

#endif

   call t_stopf ('grcalc')

   call t_startf('spegrd_alloc')
#if ( defined SPMD )
   nlon_fft_in = 2*maxm
   allocate(fftbuf_in(nlon_fft_in,8,plevp,plat))
#else
   nlon_fft_in = 1
   allocate(fftbuf_in(1,1,1,1))
#endif

   nlon_fft_out = plondfft
   allocate(fftbuf_out(nlon_fft_out,8,plevp,beglat:endlat))
   call t_stopf('spegrd_alloc')
!
   call t_startf('spegrd_bft')
!$OMP PARALLEL DO PRIVATE (LAT, IROW)
   do lat=1,plat
      irow = lat
      if (lat > plat/2) irow = plat - lat + 1
#if ( defined SPMD )
      call spegrd_bft (lat, nlon_fft_in, &
                       grdpss(1,irow), grzs(1,1,irow), grds(1,1,irow), gruhs(1,1,irow), grvhs(1,1,irow), &
                       grths(1,1,irow), grpss(1,irow), grus(1,1,irow), grvs(1,1,irow), grts(1,1,irow), &
                       grpls(1,irow), grpms(1,irow), grdpsa(1,irow), grza(1,1,irow), grda(1,1,irow), &
                       gruha(1,1,irow), grvha(1,1,irow), grtha(1,1,irow), grpsa(1,irow), grua(1,1,irow), &
                       grva(1,1,irow), grta(1,1,irow), grpla(1,irow), grpma(1,irow), fftbuf_in(1,1,1,lat) )
#else
      call spegrd_bft (lat, nlon_fft_out, &
                       grdpss(1,irow), grzs(1,1,irow), grds(1,1,irow), gruhs(1,1,irow), grvhs(1,1,irow), &
                       grths(1,1,irow), grpss(1,irow), grus(1,1,irow), grvs(1,1,irow), grts(1,1,irow), &
                       grpls(1,irow), grpms(1,irow), grdpsa(1,irow), grza(1,1,irow), grda(1,1,irow), &
                       gruha(1,1,irow), grvha(1,1,irow), grtha(1,1,irow), grpsa(1,irow), grua(1,1,irow), &
                       grva(1,1,irow), grta(1,1,irow), grpla(1,irow), grpma(1,irow), fftbuf_out(1,1,1,lat) )
#endif
   end do
   call t_stopf('spegrd_bft')

   call t_startf('spegrd_ift')
   call spegrd_ift ( nlon_fft_in, nlon_fft_out, fftbuf_in, fftbuf_out )
   call t_stopf('spegrd_ift')
                   
   call t_startf('spegrd_aft')
#ifdef OUTER_OMP
!$OMP PARALLEL DO PRIVATE (LAT)
#endif
   do lat=beglat,endlat
      call spegrd_aft (ztodt, lat, nlon(lat), nlon_fft_out, &
                   cwava(lat), qfcst(1,1,1,lat), etamid, ps(1,lat,n3), &
                   u3(1,1,lat,n3), v3(1,1,lat,n3), t3(1,1,lat,n3), &
                   qminus(1,1,1,lat), vort(1,1,lat,n3), div(1,1,lat,n3), hw2al(1,lat), hw2bl(1,lat), &
                   hw3al(1,lat), hw3bl(1,lat), hwxal(1,1,lat), hwxbl(1,1,lat), q3(1,1,1,lat,n3m1), &
                   dps(1,lat), dpsl(1,lat), dpsm(1,lat), t3(1,1,lat,n3m2) ,engy2alat(lat), engy2blat(lat), &
                   difftalat(lat), difftblat(lat), phis(1,lat), fftbuf_out(1,1,1,lat) )
                   
   end do
   call t_stopf('spegrd_aft')
!
   call t_startf('spegrd_dealloc')
   deallocate(fftbuf_in)
   deallocate(fftbuf_out)
   call t_stopf('spegrd_dealloc')
!
#ifdef SPMD
   call t_barrierf ('sync_realloc5', mpicom)
   call t_startf('realloc5')
   call realloc5 (hw2al   ,hw2bl   ,hw3al   ,hw3bl   ,tmass    , &
                  hw1lat  ,hwxal   ,hwxbl   ,engy1lat,engy2alat, &
                  engy2blat, difftalat, difftblat)
   call t_stopf('realloc5')
#endif

else

  do lat=beglat,endlat
      j = lat
      irow = lat
      if (lat > plat/2) irow = plat - lat + 1
      call forecast(lat, ps(1,lat,n3m1), ps(1,lat,n3m2), ps(1,lat,n3), &
                      u3(1,1,j,n3),u3(1,1,j,n3m1), u3(1,1,j,n3m2), &
		      v3(1,1,j,n3), v3(1,1,j,n3m1), v3(1,1,j,n3m2), &
		      t3(1,1,j,n3),t3(1,1,j,n3m1), t3(1,1,j,n3m2), &
                      q3(1,1,1,j,n3), q3(1,1,1,j,n3m1), q3(1,1,1,j,n3m2), &
	   	      ztodt, t2(1,1,lat), &
		      fu(1,1,lat), fv(1,1,lat), qfcst(1,1,1,lat), etamid,cwava(lat), &
                      qminus(1,1,1,j), hw2al(1,lat), hw2bl(1,lat), &
                      hw3al(1,lat), hw3bl(1,lat), hwxal(1,1,lat), &
                      hwxbl(1,1,lat), &
                      nlon(lat)) 

   end do

!
! Initialize fixer variables for routines not called in scam version of
! model
!
	engy2alat=0._r8
	engy2blat=0._r8
	difftalat=0._r8
	difftblat=0._r8
        engy2b=0._r8

endif ! if not SCAM
!
! Accumulate and normalize global integrals for mass fixer (dry mass of
! atmosphere is held constant).
!
   call t_startf ('scan2_single')
   tmassf = 0._r8
   do lat=1,plat
      tmassf = tmassf + tmass(lat)
   end do
   tmassf = tmassf*.5_r8
!
! Initialize moisture, mass, energy, and temperature integrals
!
   hw1(1) = 0._r8
   engy1  = 0._r8
   engy2a = 0._r8
   engy2b = 0._r8
   diffta = 0._r8
   difftb = 0._r8
   do m=1,pcnst
      hw2a(m) = 0._r8
      hw2b(m) = 0._r8
      hw3a(m) = 0._r8
      hw3b(m) = 0._r8
      do n=1,4
         hwxa(m,n) = 0._r8
         hwxb(m,n) = 0._r8
      end do
   end do
!
! Sum water and energy integrals over latitudes
!
   do lat=1,plat
      engy1   = engy1   + engy1lat (lat)
      engy2a  = engy2a  + engy2alat(lat)
      engy2b  = engy2b  + engy2blat(lat)
      diffta  = diffta  + difftalat(lat)
      difftb  = difftb  + difftblat(lat)
      hw1(1)  = hw1(1)  + hw1lat(1,lat)
      hw2a(1) = hw2a(1) + hw2al(1,lat)
      hw2b(1) = hw2b(1) + hw2bl(1,lat)
      hw3a(1) = hw3a(1) + hw3al(1,lat)
      hw3b(1) = hw3b(1) + hw3bl(1,lat)
   end do
!
! Compute atmospheric mass fixer coefficient
!
   qmassf     = hw1(1)
   if (adiabatic .or. ideal_phys) then
      fixmas = tmass0/tmassf
   else
      fixmas = (tmass0 + qmassf)/tmassf
   end if
!
! Compute alpha for water ONLY
!
   hw2(1)    = hw2a(1) + fixmas*hw2b(1)
   hw3(1)    = hw3a(1) + fixmas*hw3b(1)
   if(hw3(1) .ne. 0._r8) then
      alpha(1)  = ( hw1(1) - hw2(1) )/hw3(1)
   else
      alpha(1)  = 1._r8
   endif
!
! Compute beta for energy
!
   engy2    = engy2a + fixmas*engy2b
   difft    = diffta + fixmas*difftb
   residual = (engy2 - engy1)/ztodt
   if(difft .ne. 0._r8) then
     beta = -residual*ztodt/(cpair*difft)
   else
     beta = 0._r8
   endif
!!   write(iulog,125) residual,beta
!!125 format('      resid, beta      = ',25x,2f25.15)
!
! Compute alpha for non-water constituents
!
   do m = 2,pcnst
      hw1(m) = 0._r8
      do lat=1,plat
         hw1(m) = hw1(m) + hw1lat(m,lat)
      end do
      do n = 1,4
         do lat=1,plat
            hwxa(m,n) = hwxa(m,n) + hwxal(m,n,lat)
            hwxb(m,n) = hwxb(m,n) + hwxbl(m,n,lat)
         end do
      end do
      hw2a(m) = hwxa(m,1) - alpha(1)*hwxa(m,2)
      hw2b(m) = hwxb(m,1) - alpha(1)*hwxb(m,2)
      hw3a(m) = hwxa(m,3) - alpha(1)*hwxa(m,4)
      hw3b(m) = hwxb(m,3) - alpha(1)*hwxb(m,4)
      hw2 (m) = hw2a(m) + fixmas*hw2b(m)
      hw3 (m) = hw3a(m) + fixmas*hw3b(m)
      if(hw3(m) .ne. 0._r8) then
         alpha(m)  = ( hw1(m) - hw2(m) )/hw3(m)
      else
         alpha(m)  = 1._r8
      end if
   end do

   call t_stopf ('scan2_single')

   call t_startf ('tfilt_massfix')

   if (single_column) then
!
! read in fixer for scam
!
      if ( use_iop ) then
         fixmas=fixmascam
         beta=betacam
         do m = 1, pcnst
            alpha(m)=alphacam(m)
         end do
      else
         fixmas=1._r8
         beta=0._r8
         alpha(:)=0._r8
      end if
   endif

#ifdef OUTER_OMP
!$OMP PARALLEL DO PRIVATE (LAT)
#endif
   do lat=beglat,endlat

      call tfilt_massfixrun (ztodt,       lat, u3(1,1,lat,n3m1),u3(1,1,lat,n3), &
                          v3(1,1,lat,n3m1), v3(1,1,lat,n3), t3(1,1,lat,n3m1), t3(1,1,lat,n3), &
                             q3(1,1,1,lat,n3m1), &
                          q3(1,1,1,lat,n3), ps(1,lat,n3m1), ps(1,lat,n3), alpha, &
                          etamid, qfcst(1,1,1,lat), vort(1,1,lat,n3), div(1,1,lat,n3), &
                             vort(1,1,lat,n3m2), &
                          div(1,1,lat,n3m2), qminus(1,1,1,lat), ps(1,lat,n3m2), &
                             u3(1,1,lat,n3m2), &
                          v3(1,1,lat,n3m2), t3(1,1,lat,n3m2), q3(1,1,1,lat,n3m2), vort(1,1,lat,n3m1), &
                             div(1,1,lat,n3m1), &
                          omga(1,1,lat), dpsl(1,lat), dpsm(1,lat), beta, hadv(1,1,1,lat) ,nlon(lat), &
                          pdeld(:,:,lat,n3), pdeld(:,:,lat,n3m1), pdeld(:,:,lat,n3m2))
                          
   end do
   call t_stopf ('tfilt_massfix')
!
! Shift time pointers
!
   call shift_time_indices ()

   return
end subroutine scan2run

!
!-----------------------------------------------------------------------
!

#ifdef SPMD
subroutine realloc5 (hw2al   ,hw2bl   ,hw3al   ,hw3bl   ,tmass    , &
                     hw1lat  ,hwxal   ,hwxbl   ,engy1lat,engy2alat, &
                     engy2blat,difftalat,difftblat      )
!-----------------------------------------------------------------------
!
! Purpose: Reallocation routine for slt variables.
!
! Method: MPI_Allgatherv (or point-to-point implementation)
! 
! Author:  J. Rosinski
! Standardized:      J. Rosinski, Oct 1995
!                    J. Truesdale, Feb. 1996
! Modified: P. Worley, December 2003, October 2004
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------
   use pmgrid, only: numlats, plat
   use mpishorthand, only: mpicom, mpir8
   use spmd_dyn
   use spmd_utils, only : iam, npes, altalltoallv
!---------------------------------Parameters----------------------------------
   integer, parameter :: msgtag  = 5000
!---------------------------------Commons-------------------------------------
#include <comsta.h>
!-----------------------------------------------------------------------
!
! Input arguments
!
   real(r8), intent(inout) :: hw2al(pcnst,plat)
   real(r8), intent(inout) :: hw2bl(pcnst,plat)
   real(r8), intent(inout) :: hw3al(pcnst,plat)
   real(r8), intent(inout) :: hw3bl(pcnst,plat)
   real(r8), intent(inout) :: tmass (plat)
   real(r8), intent(inout) :: hw1lat(pcnst,plat)
   real(r8), intent(inout) :: hwxal(pcnst,4,plat)
   real(r8), intent(inout) :: hwxbl(pcnst,4,plat)
!                                                ! -
   real(r8), intent(inout)   :: engy1lat (plat)  ! lat contribution to total energy (n)
   real(r8), intent(inout)   :: engy2alat(plat)  ! lat contribution to total energy (n+1)
   real(r8), intent(inout)   :: engy2blat(plat)  ! lat contribution to total energy (n+1)
   real(r8), intent(inout)   :: difftalat(plat)  ! lat contribution to delta-T integral
   real(r8), intent(inout)   :: difftblat(plat)  ! lat contribution to delta-T integral
!
!---------------------------Local workspace-----------------------------
!
   integer procid
   integer bufpos
   integer procj
   integer step, i, j, m, jstrt
   integer beglat_p, endlat_p, numlats_p, jstrt_p
!
   logical, save :: first = .true.
   integer, save :: sndcnt
   integer, allocatable, save :: sndcnts(:), sdispls(:)
   integer, allocatable, save :: rcvcnts(:), rdispls(:)
   integer, allocatable, save :: pdispls(:)
!-----------------------------------------------------------------------
   if (first) then
! Compute send/recv/put counts and displacements
      allocate(sndcnts(0:npes-1))
      allocate(sdispls(0:npes-1))
      allocate(rcvcnts(0:npes-1))
      allocate(rdispls(0:npes-1))
      allocate(pdispls(0:npes-1))
!
! Compute send count
      sndcnt = (pcnst*(5 + 2*4) + 6)*numlats
      sndcnts(:) = 0
      do step=1,allgather_steps
         procid = allgather_proc(step)
         sndcnts(procid) = sndcnt
      enddo
!   
      sdispls(0) = 0
      do procid=1,npes-1
        sdispls(procid) = 0
      enddo
!
! Compute recv counts and displacements
      rcvcnts(:) = 0
      do step=1,allgather_steps
         procid = allgather_proc(step)
         rcvcnts(procid) = (pcnst*(5 + 2*4) + 6)*nlat_p(procid)
      enddo
      rcvcnts(iam) = (pcnst*(5 + 2*4) + 6)*numlats
!   
      rdispls(0) = 0
      do procid=1,npes-1
        rdispls(procid) = rdispls(procid-1) + rcvcnts(procid-1)
      enddo
!
      pdispls(:) = 0
      call mpialltoallint(rdispls, 1, pdispls, 1, mpicom)
!
      first = .false.
   endif
!
! Fill send buffer
   jstrt = beglat - 1
   bufpos = 0
! tmass
   do j=1,numlats
      buf1(bufpos+j) = tmass(jstrt+j)
   enddo
   bufpos = bufpos + numlats
! engy1lat
   do j=1,numlats
      buf1(bufpos+j) = engy1lat(jstrt+j)
   enddo
   bufpos = bufpos + numlats
! engy2alat
   do j=1,numlats
      buf1(bufpos+j) = engy2alat(jstrt+j)
   enddo
   bufpos = bufpos + numlats
! engy2blat
   do j=1,numlats
      buf1(bufpos+j) = engy2blat(jstrt+j)
   enddo
   bufpos = bufpos + numlats
! difftalat
   do j=1,numlats
      buf1(bufpos+j) = difftalat(jstrt+j)
   enddo
   bufpos = bufpos + numlats
! difftblat
   do j=1,numlats
      buf1(bufpos+j) = difftblat(jstrt+j)
   enddo
   bufpos = bufpos + numlats
!hw1lat
   do j=beglat,endlat
      do m=1,pcnst
         buf1(bufpos+m) = hw1lat(m,j)
      enddo
      bufpos = bufpos + pcnst
   enddo
!hw2al
   do j=beglat,endlat
      do m=1,pcnst
         buf1(bufpos+m) = hw2al(m,j)
      enddo
      bufpos = bufpos + pcnst
   enddo
!hw2bl
   do j=beglat,endlat
      do m=1,pcnst
         buf1(bufpos+m) = hw2bl(m,j)
      enddo
      bufpos = bufpos + pcnst
   enddo
!hw3al
   do j=beglat,endlat
      do m=1,pcnst
         buf1(bufpos+m) = hw3al(m,j)
      enddo
      bufpos = bufpos + pcnst
   enddo
!hw3bl
   do j=beglat,endlat
      do m=1,pcnst
         buf1(bufpos+m) = hw3bl(m,j)
      enddo
      bufpos = bufpos + pcnst
   enddo
!hwxal
   do j=beglat,endlat
      do i=1,4
         do m=1,pcnst
            buf1(bufpos+m) = hwxal(m,i,j)
         enddo
         bufpos = bufpos + pcnst
      enddo
   enddo
!hwxbl
   do j=beglat,endlat
      do i=1,4
         do m=1,pcnst
            buf1(bufpos+m) = hwxbl(m,i,j)
         enddo
         bufpos = bufpos + pcnst
      enddo
   enddo
!
! Gather the data
!
   if (dyn_allgather .eq. 0) then
      call mpiallgatherv(buf1, sndcnt, mpir8, &
                         buf2, rcvcnts, rdispls, mpir8, &
                         mpicom)
   else
      call altalltoallv(dyn_allgather, iam, npes, &
                        allgather_steps, allgather_proc, &
                        buf1, spmdbuf_siz, sndcnts, sdispls, mpir8, &
                        buf2, spmdbuf_siz, rcvcnts, rdispls, mpir8, &
                        msgtag, pdispls, mpir8, buf2win, mpicom)
   endif
!
! Copy out of message buffers
!
!$OMP PARALLEL DO PRIVATE (STEP, PROCID, BEGLAT_P, ENDLAT_P, NUMLATS_P, BUFPOS, JSTRT_P, I, J, M)
   do step=1,allgather_steps
      procid = allgather_proc(step)
      beglat_p = cut(1,procid)
      endlat_p = cut(2,procid)
      numlats_p = nlat_p(procid)
      bufpos = rdispls(procid)
! tmass
      jstrt_p  = beglat_p - 1
      do j=1,numlats_p
         tmass(jstrt_p+j) = buf2(bufpos+j)
      enddo
      bufpos = bufpos + numlats_p
! engy1lat
      jstrt_p  = beglat_p - 1
      do j=1,numlats_p
         engy1lat(jstrt_p+j) = buf2(bufpos+j)
      enddo
      bufpos = bufpos + numlats_p
! engy2alat
      jstrt_p  = beglat_p - 1
      do j=1,numlats_p
         engy2alat(jstrt_p+j) = buf2(bufpos+j)
      enddo
      bufpos = bufpos + numlats_p
! engy2blat
      jstrt_p  = beglat_p - 1
      do j=1,numlats_p
         engy2blat(jstrt_p+j) = buf2(bufpos+j)
      enddo
      bufpos = bufpos + numlats_p
! difftalat
      jstrt_p  = beglat_p - 1
      do j=1,numlats_p
         difftalat(jstrt_p+j) = buf2(bufpos+j)
      enddo
      bufpos = bufpos + numlats_p
! difftblat
      jstrt_p  = beglat_p - 1
      do j=1,numlats_p
         difftblat(jstrt_p+j) = buf2(bufpos+j)
      enddo
      bufpos = bufpos + numlats_p
! hw1lat
      do j=beglat_p,endlat_p
         do m=1,pcnst
            hw1lat(m,j) = buf2(bufpos+m)
         enddo
         bufpos = bufpos + pcnst
      enddo
! hw2al
      do j=beglat_p,endlat_p
         do m=1,pcnst
            hw2al(m,j) = buf2(bufpos+m)
         enddo
         bufpos = bufpos + pcnst
      enddo
! hw2bl
      do j=beglat_p,endlat_p
         do m=1,pcnst
            hw2bl(m,j) = buf2(bufpos+m)
         enddo
         bufpos = bufpos + pcnst
      enddo
! hw3al
      do j=beglat_p,endlat_p
         do m=1,pcnst
            hw3al(m,j) = buf2(bufpos+m)
         enddo
         bufpos = bufpos + pcnst
      enddo
! hw3bl
      do j=beglat_p,endlat_p
         do m=1,pcnst
            hw3bl(m,j) = buf2(bufpos+m)
         enddo
         bufpos = bufpos + pcnst
      enddo
! hwxal
      do j=beglat_p,endlat_p
         do i=1,4
            do m=1,pcnst
               hwxal(m,i,j) = buf2(bufpos+m)
            enddo
            bufpos = bufpos + pcnst
         enddo
      enddo
! hwxbl
      do j=beglat_p,endlat_p
         do i=1,4
            do m=1,pcnst
               hwxbl(m,i,j) = buf2(bufpos+m)
            enddo
            bufpos = bufpos + pcnst
         enddo
      enddo
!
   end do
!
   return
end subroutine realloc5
#endif

!
!-----------------------------------------------------------------------
!


end module scan2
