
module spetru

!----------------------------------------------------------------------- 
! 
! Purpose: Spectrally truncate initial data fields.
!
! Method: Truncate one or a few fields at a time, to minimize the 
!         memory  requirements
! 
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, J. Hack, August 1992
! Modified to implement processing of subsets of fields: P. Worley, May 2003
! 
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plon, plev, plat
  use pspect,       only: psp, pspt, ptrn, pmmax
  use comspe,       only: alp, nlen, nstart, dalp
  use rgrid,        only: nlon, nmmax
  use commap,       only: w, xm
  use physconst,    only: rearth, ra
  use eul_control_mod, only: trig, ifax, pcray
  implicit none
!
! By default make data and interfaces to this module private
!
  private

!
! Public interfaces
!
  public spetru_phis      ! Spectrally truncate PHIS
  public spetru_ps        ! Spectrally truncate PS
  public spetru_3d_scalar ! Spectrally truncate 3D scalar fields
  public spetru_uv        ! Spectrally truncate winds (U and V)
!
! Private module data
!
  integer, parameter :: plondfft = plon + 2   ! Size of longitude needed for FFT's

!
!======================================================================= 
contains

!************************************************************************
subroutine spetru_phis  (phis, phis_hires, phisl, phism, phi_out)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! 
! Method: 
! Spectrally truncate PHIS input field.
! 
! Author: 
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, J. Hack, August 1992
! Modified:          P. Worley, May 2003
! Modified:          J. Olson, Apr 2004
!
!-----------------------------------------------------------------------

   use pmgrid,   only: plon, plat

!
! Input/Output arguments
!
   real(r8), intent(inout) :: phis(plon,plat)            ! Fourier -> spec. coeffs. for sfc geo.
   logical,  intent(in)    :: phis_hires                 ! true => PHIS came from hi res topo file
   real(r8), intent(out), optional :: phisl(plon,plat)   ! Spectrally trunc d(phis)/d(longitude)
   real(r8), intent(out), optional :: phism(plon,plat)   ! Spectrally trunc d(phis)/d(latitude)
   real(r8), intent(out), optional :: phi_out(2,psp/2)   ! used in spectral truncation of phis
!
!---------------------------Local workspace-----------------------------
!
   real(r8), pointer :: phis_tmp(:,:)      ! Temporary to compute Phis of size needed for FFT
   real(r8), pointer :: phisl_tmp(:,:)     ! Temporary to compute phisl of size needed for FFT
   real(r8), pointer :: phism_tmp(:,:)     ! Temporary to compute phism of size needed for FFT
   real(r8) tmp1               ! vector temporary
   real(r8) tmp2               ! vector temporary
   real(r8) phialpr,phialpi    ! phi*alp (real and imaginary)
   real(r8) phdalpr,phdalpi    ! phi*dalp (real and imaginary)
   real(r8) zwalp              ! zw*alp
   real(r8) zw                 ! w**2
   real(r8) filtlim            ! filter function
   real(r8) ft                 ! filter multiplier for spectral coefficients
   real(r8) phi(2,psp/2)       ! used in spectral truncation of phis
#if ( ! defined USEFFTLIB )
   real(r8) work((plon+1)*plev)  ! Workspace for fft
#else
   real(r8) work((plon+1)*pcray)   ! Workspace for fft
#endif

   integer i                   ! longitude index
   integer irow                ! latitude pair index
   integer latm,latp           ! symmetric latitude indices
   integer lat
   integer m                   ! longitudinal wavenumber index
   integer n                   ! latitudinal wavenumber index
   integer nspec
   integer mr                  ! spectral indices
!
!-----------------------------------------------------------------------
!
! Zero spectral array
!
   phi(:,:) = 0._r8
!
! Transform grid -> fourier
!
   allocate(phis_tmp(plondfft,plat))
   phis_tmp(:plon,:) = phis(:plon,:)
   do lat=1,plat
      irow = lat
      if (lat.gt.plat/2) irow = plat - lat + 1
      call fft991(phis_tmp(1,lat),work,trig(1,irow),ifax(1,irow),1,plondfft, &
                  nlon(lat),1,-1)
   end do                    ! lat=1,plat
!
! Loop over latitude pairs
!
   do irow=1,plat/2
      latp = irow
      latm = plat - irow + 1
      zw = w(irow)*2._r8
      do i=1,2*nmmax(irow)
!
! Compute symmetric and antisymmetric components
!
         tmp1 = 0.5_r8*(phis_tmp(i,latm) - phis_tmp(i,latp))
         tmp2 = 0.5_r8*(phis_tmp(i,latm) + phis_tmp(i,latp))
         phis_tmp(i,latm) = tmp1
         phis_tmp(i,latp) = tmp2
      end do
!     
! Compute phi*mn
!
      do m=1,nmmax(irow)
         mr = nstart(m)
         do n=1,nlen(m),2
            zwalp = zw*alp(mr+n,irow)
            phi(1,mr+n) = phi(1,mr+n) + zwalp*phis_tmp(2*m-1,latp)
            phi(2,mr+n) = phi(2,mr+n) + zwalp*phis_tmp(2*m  ,latp)
         end do

         do n=2,nlen(m),2
            zwalp = zw*alp(mr+n,irow)
            phi(1,mr+n) = phi(1,mr+n) + zwalp*phis_tmp(2*m-1,latm)
            phi(2,mr+n) = phi(2,mr+n) + zwalp*phis_tmp(2*m  ,latm)
         end do
      end do
   enddo                  ! irow=1,plat/2
!
   if (phis_hires) then
!
! Apply spectral filter to phis
!     filter is a function of n 
!        if n < filter limit then
!           spectral_coeff = spectral_coeff * (1. - (real(n,r8)/filtlim)**2)
!        else         
!           spectral_coeff = 0.
!        endif
!     where filter limit = 1.4*PTRN
!     
      filtlim = real(int(1.4_r8*real(ptrn,r8)),r8)
      do m=1,pmmax
         mr = nstart(m)
         do n=1,nlen(m)
            nspec=m-1+n
            ft = 1._r8 - (real(nspec,r8)/filtlim)**2
            if (real(nspec,r8) .ge. filtlim) ft = 0._r8
            phi(1,mr+n) = phi(1,mr+n)*ft 
            phi(2,mr+n) = phi(2,mr+n)*ft 
         end do
      end do
      call hordif1(rearth,phi)
   end if
!
! Compute grid point values of phi*.
!
   do irow=1,plat/2
      latp = irow
      latm = plat - irow + 1
!
! Zero fourier fields
!
      phis_tmp(:,latm) = 0._r8
      phis_tmp(:,latp) = 0._r8
!
! Compute(phi*)m
!
      do m=1,nmmax(irow)
         mr = nstart(m)
         do n=1,nlen(m),2
            phialpr = phi(1,mr+n)*alp(mr+n,irow)
            phialpi = phi(2,mr+n)*alp(mr+n,irow)
            phis_tmp(2*m-1,latm) = phis_tmp(2*m-1,latm) + phialpr
            phis_tmp(2*m  ,latm) = phis_tmp(2*m  ,latm) + phialpi
         end do
      end do

      do m=1,nmmax(irow)
         mr = nstart(m)
         do n=2,nlen(m),2
            phialpr = phi(1,mr+n)*alp(mr+n,irow)
            phialpi = phi(2,mr+n)*alp(mr+n,irow)
            phis_tmp(2*m-1,latp) = phis_tmp(2*m-1,latp) + phialpr
            phis_tmp(2*m  ,latp) = phis_tmp(2*m  ,latp) + phialpi
         end do
      end do
!
! Recompute real fields from symmetric and antisymmetric parts
!
      do i=1,nlon(latm)+2
         tmp1 = phis_tmp(i,latm) + phis_tmp(i,latp)
         tmp2 = phis_tmp(i,latm) - phis_tmp(i,latp)
         phis_tmp(i,latm) = tmp1
         phis_tmp(i,latp) = tmp2
      end do

   enddo                 ! irow=1,plat/2

   if(present(phisl)) then
      allocate(phisl_tmp(plondfft,plat))
      do irow=1,plat/2
         latp = irow
         latm = plat - irow + 1
!
! Zero fourier fields
!
         phisl_tmp(:,latm) = 0._r8
         phisl_tmp(:,latp) = 0._r8
!
! Compute(phi*)m
!
         do m=1,nmmax(irow)
            mr = nstart(m)
            do n=1,nlen(m),2
               phialpr = phi(1,mr+n)*alp(mr+n,irow)
               phialpi = phi(2,mr+n)*alp(mr+n,irow)
               phisl_tmp(2*m-1,latm) = phisl_tmp(2*m-1,latm) - phialpi*ra
               phisl_tmp(2*m  ,latm) = phisl_tmp(2*m  ,latm) + phialpr*ra
            end do
         end do

         do m=1,nmmax(irow)
            mr = nstart(m)
            do n=2,nlen(m),2
               phialpr = phi(1,mr+n)*alp(mr+n,irow)
               phialpi = phi(2,mr+n)*alp(mr+n,irow)
               phisl_tmp(2*m-1,latp) = phisl_tmp(2*m-1,latp) - phialpi*ra
               phisl_tmp(2*m  ,latp) = phisl_tmp(2*m  ,latp) + phialpr*ra
            end do
         end do
!
! d(Phi)/d(lamda)
!
         do m=1,nmmax(irow)
            phisl_tmp(2*m-1,latm) = xm(m)*phisl_tmp(2*m-1,latm)
            phisl_tmp(2*m  ,latm) = xm(m)*phisl_tmp(2*m  ,latm)
            phisl_tmp(2*m-1,latp) = xm(m)*phisl_tmp(2*m-1,latp)
            phisl_tmp(2*m  ,latp) = xm(m)*phisl_tmp(2*m  ,latp)
         end do
!
! Recompute real fields from symmetric and antisymmetric parts
!
         do i=1,nlon(latm)+2
            tmp1 = phisl_tmp(i,latm) + phisl_tmp(i,latp)
            tmp2 = phisl_tmp(i,latm) - phisl_tmp(i,latp)
            phisl_tmp(i,latm) = tmp1
            phisl_tmp(i,latp) = tmp2
         end do
      enddo                 ! irow=1,plat/2
   end if

   if(present(phism)) then
      allocate(phism_tmp(plondfft,plat))
      do irow=1,plat/2
         latp = irow
         latm = plat - irow + 1
!
! Zero fourier fields
!
         phism_tmp(:,latm) = 0._r8
         phism_tmp(:,latp) = 0._r8
!
! Compute(phi*)m
!
         do m=1,nmmax(irow)
            mr = nstart(m)
            do n=1,nlen(m),2
               phdalpr = phi(1,mr+n)*dalp(mr+n,irow)
               phdalpi = phi(2,mr+n)*dalp(mr+n,irow)
               phism_tmp(2*m-1,latp) = phism_tmp(2*m-1,latp) + phdalpr*ra
               phism_tmp(2*m  ,latp) = phism_tmp(2*m  ,latp) + phdalpi*ra
            end do
         end do

         do m=1,nmmax(irow)
            mr = nstart(m)
            do n=2,nlen(m),2
               phdalpr = phi(1,mr+n)*dalp(mr+n,irow)
               phdalpi = phi(2,mr+n)*dalp(mr+n,irow)
               phism_tmp(2*m-1,latm) = phism_tmp(2*m-1,latm) + phdalpr*ra
               phism_tmp(2*m  ,latm) = phism_tmp(2*m  ,latm) + phdalpi*ra
            end do
         end do
!
! Recompute real fields from symmetric and antisymmetric parts
!
         do i=1,nlon(latm)+2
            tmp1 = phism_tmp(i,latm) + phism_tmp(i,latp)
            tmp2 = phism_tmp(i,latm) - phism_tmp(i,latp)
            phism_tmp(i,latm) = tmp1
            phism_tmp(i,latp) = tmp2
         end do
      enddo                 ! irow=1,plat/2
   end if
!
   do lat=1,plat
!     
! Transform Fourier -> grid, obtaining spectrally truncated
! grid point values.
!
      irow = lat
      if (lat.gt.plat/2) irow = plat - lat + 1

      call fft991(phis_tmp(1,lat),work,trig(1,irow),ifax(1,irow),1,plondfft, &
                  nlon(lat),1,+1)
      phis(:nlon(lat),lat) = phis_tmp(:nlon(lat),lat)
      if(present(phisl)) then
         call fft991 (phisl_tmp(1,lat),work     ,trig(1,irow),ifax(1,irow),1       , &
                      plondfft       ,nlon(lat),1           ,+1      )
         phisl(:nlon(lat),lat) = phisl_tmp(:nlon(lat),lat)
      end if
      if(present(phism)) then
         call fft991 (phism_tmp(1,lat),work     ,trig(1,irow),ifax(1,irow),1       , &
                      plondfft       ,nlon(lat),1           ,+1      )
         phism(:nlon(lat),lat) = phism_tmp(:nlon(lat),lat)
      end if
   enddo
   deallocate( phis_tmp )
   if ( present(phisl) ) deallocate( phisl_tmp )
   if ( present(phism) ) deallocate( phism_tmp )

   if(present(phi_out)) then
      phi_out(:,:) = phi(:,:)
   end if

   return
end subroutine spetru_phis

!************************************************************************
subroutine spetru_ps(ps      ,dpsl    ,dpsm)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! 
! Method: 
! Spectrally truncate PS input field.
! 
! Author: 
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, J. Hack, August 1992
! Modified:          P. Worley, May 2003
! Modified:          J. Olson, Apr 2004
!
!-----------------------------------------------------------------------

   use pmgrid,   only: plon, plat

!
! Input/Output arguments
!
   real(r8), intent(inout) :: ps(plon,plat)          ! Fourier -> spec. coeffs. for ln(ps)
!
! Output arguments
!
   real(r8), intent(out) :: dpsl(plon,plat)          ! Spectrally trunc d(ln(ps))/d(longitude)
   real(r8), intent(out) :: dpsm(plon,plat)          ! Spectrally trunc d(ln(ps))/d(latitude)

!
!---------------------------Local workspace-----------------------------
!
   real(r8), pointer :: log_ps(:,:)    ! log(ps)
   real(r8), pointer :: dpsm_tmp(:,:)  ! Temporary to compute dpsm of size needed for FFT
   real(r8), pointer :: dpsl_tmp(:,:)  ! Temporary to compute dpsl of size needed for FFT
   real(r8) alps_tmp(psp)      ! used in spectral truncation of phis
   real(r8) tmp1               ! vector temporary
   real(r8) tmp2               ! vector temporary
   real(r8) zwalp              ! zw*alp
   real(r8) psdalpr,psdalpi    ! alps (real and imaginary)*dalp
   real(r8) psalpr,psalpi      ! alps (real and imaginary)*alp
   real(r8) zw                 ! w**2
#if ( ! defined USEFFTLIB )
   real(r8) work((plon+1)*plev)  ! Workspace for fft
#else
   real(r8) work((plon+1)*pcray)   ! Workspace for fft
#endif

   integer ir,ii               ! indices complex coeffs. of spec. arrs.
   integer i,k                 ! longitude, level indices
   integer irow                ! latitude pair index
   integer latm,latp           ! symmetric latitude indices
   integer lat
   integer m                   ! longitudinal wavenumber index
   integer n                   ! latitudinal wavenumber index
   integer nspec
   integer mr,mc               ! spectral indices
!
!-----------------------------------------------------------------------
!
! Zero spectral array
!
   alps_tmp(:) = 0._r8
!
! Compute the 2D quantities which are transformed to spectral space:
!  ps= ln(ps). 
!
   allocate( log_ps(plondfft,plat) )
   do lat=1,plat
      irow = lat
      if (lat.gt.plat/2) irow = plat - lat + 1
      do i=1,nlon(lat)
         log_ps(i,lat) = log(ps(i,lat))
      end do
!
! Transform grid -> fourier
!
      call fft991(log_ps(1,lat),work,trig(1,irow),ifax(1,irow),1,plondfft, &
                  nlon(lat),1,-1)

   end do                    ! lat=1,plat
   allocate( dpsl_tmp(plondfft,plat) )
   allocate( dpsm_tmp(plondfft,plat) )
!
! Loop over latitude pairs
!
   do irow=1,plat/2
      latp = irow
      latm = plat - irow + 1
      zw = w(irow)*2._r8
      do i=1,2*nmmax(irow)
!
! Compute symmetric and antisymmetric components
!
         tmp1 = 0.5_r8*(log_ps(i,latm) - log_ps(i,latp))
         tmp2 = 0.5_r8*(log_ps(i,latm) + log_ps(i,latp))
         log_ps(i,latm) = tmp1
         log_ps(i,latp) = tmp2

      end do
!     
! Compute ln(p*)mn
!
      do m=1,nmmax(irow)
         mr = nstart(m)
         mc = 2*mr
         do n=1,nlen(m),2
            zwalp = zw*alp(mr+n,irow)
            ir = mc + 2*n - 1
            ii = ir + 1
            alps_tmp(ir) = alps_tmp(ir) + zwalp*log_ps(2*m-1,latp)
            alps_tmp(ii) = alps_tmp(ii) + zwalp*log_ps(2*m  ,latp)
         end do

         do n=2,nlen(m),2
            zwalp = zw*alp(mr+n,irow)
            ir = mc + 2*n - 1
            ii = ir + 1
            alps_tmp(ir) = alps_tmp(ir) + zwalp*log_ps(2*m-1,latm)
            alps_tmp(ii) = alps_tmp(ii) + zwalp*log_ps(2*m  ,latm)
         end do
      end do
   enddo                  ! irow=1,plat/2
!
! Compute grid point values of:ln(p*) and grad(ln(p*)).
!
   do irow=1,plat/2
      latp = irow
      latm = plat - irow + 1
!
! Zero fourier fields
!
      log_ps(:,latm) = 0._r8
      log_ps(:,latp) = 0._r8

      dpsl_tmp(:,latm) = 0._r8
      dpsl_tmp(:,latp) = 0._r8

      dpsm_tmp(:,latm) = 0._r8
      dpsm_tmp(:,latp) = 0._r8

!
! Compute(ln(p*),grad(ln(p*)))m
!
      do m=1,nmmax(irow)
         mr = nstart(m)
         mc = 2*mr
         do n=1,nlen(m),2
            ir = mc + 2*n - 1
            ii = ir + 1
            psalpr = alps_tmp(ir)*alp(mr+n,irow)
            psalpi = alps_tmp(ii)*alp(mr+n,irow)
!     
            log_ps(2*m-1,latm) = log_ps(2*m-1,latm) + psalpr
            log_ps(2*m  ,latm) = log_ps(2*m  ,latm) + psalpi
            dpsl_tmp(2*m-1,latm) = dpsl_tmp(2*m-1,latm) - psalpi*ra
            dpsl_tmp(2*m  ,latm) = dpsl_tmp(2*m  ,latm) + psalpr*ra
!
            psdalpr = alps_tmp(ir)*dalp(mr+n,irow)
            psdalpi = alps_tmp(ii)*dalp(mr+n,irow)
!
            dpsm_tmp(2*m-1,latp) = dpsm_tmp(2*m-1,latp) + psdalpr*ra
            dpsm_tmp(2*m  ,latp) = dpsm_tmp(2*m  ,latp) + psdalpi*ra
         end do
      end do

      do m=1,nmmax(irow)
         mr = nstart(m)
         mc = 2*mr
         do n=2,nlen(m),2
            ir = mc + 2*n - 1
            ii = ir + 1
            psalpr = alps_tmp(ir)*alp(mr+n,irow)
            psalpi = alps_tmp(ii)*alp(mr+n,irow)
!     
            log_ps(2*m-1,latp) = log_ps(2*m-1,latp) + psalpr
            log_ps(2*m  ,latp) = log_ps(2*m  ,latp) + psalpi
            dpsl_tmp(2*m-1,latp) = dpsl_tmp(2*m-1,latp) - psalpi*ra
            dpsl_tmp(2*m  ,latp) = dpsl_tmp(2*m  ,latp) + psalpr*ra
!
            psdalpr = alps_tmp(ir)*dalp(mr+n,irow)
            psdalpi = alps_tmp(ii)*dalp(mr+n,irow)
!
            dpsm_tmp(2*m-1,latm) = dpsm_tmp(2*m-1,latm) + psdalpr*ra
            dpsm_tmp(2*m  ,latm) = dpsm_tmp(2*m  ,latm) + psdalpi*ra
         end do
      end do

      do m=1,nmmax(irow)
         dpsl_tmp(2*m-1,latm) = xm(m)*dpsl_tmp(2*m-1,latm)
         dpsl_tmp(2*m  ,latm) = xm(m)*dpsl_tmp(2*m  ,latm)
         dpsl_tmp(2*m-1,latp) = xm(m)*dpsl_tmp(2*m-1,latp)
         dpsl_tmp(2*m  ,latp) = xm(m)*dpsl_tmp(2*m  ,latp)
      end do
!
! Recompute real fields from symmetric and antisymmetric parts
!
      do i=1,nlon(latm)+2
!
         tmp1 = log_ps(i,latm) + log_ps(i,latp)
         tmp2 = log_ps(i,latm) - log_ps(i,latp)
         log_ps(i,latm) = tmp1
         log_ps(i,latp) = tmp2
!
         tmp1 = dpsl_tmp(i,latm) + dpsl_tmp(i,latp)
         tmp2 = dpsl_tmp(i,latm) - dpsl_tmp(i,latp)
         dpsl_tmp(i,latm) = tmp1
         dpsl_tmp(i,latp) = tmp2
!
         tmp1 = dpsm_tmp(i,latm) + dpsm_tmp(i,latp)
         tmp2 = dpsm_tmp(i,latm) - dpsm_tmp(i,latp)
         dpsm_tmp(i,latm) = tmp1
         dpsm_tmp(i,latp) = tmp2
      end do
!
   enddo                 ! irow=1,plat/2
!
   do lat=1,plat
!     
! Transform Fourier -> grid, obtaining spectrally truncated
! grid point values.
!
      irow = lat
      if (lat.gt.plat/2) irow = plat - lat + 1

      call fft991(log_ps(1,lat),work,trig(1,irow),ifax(1,irow),1,plondfft, &
                  nlon(lat),1,+1)
      call fft991(dpsl_tmp(1,lat),work,trig(1,irow),ifax(1,irow),1,plondfft, &
                  nlon(lat),1,+1)
      call fft991(dpsm_tmp(1,lat),work,trig(1,irow),ifax(1,irow),1,plondfft, &
                  nlon(lat),1,+1)
!
! Convert from ln(ps) to ps, copy temporaries to input arrays
!
      do i=1,nlon(lat)
         ps(i,lat) = exp(log_ps(i,lat))
         dpsl(i,lat) = dpsl_tmp(i,lat)
         dpsm(i,lat) = dpsm_tmp(i,lat)
      end do
!
   enddo
   deallocate( log_ps )
   deallocate( dpsm_tmp )
   deallocate( dpsl_tmp )

   return
end subroutine spetru_ps

!************************************************************************

subroutine spetru_3d_scalar(x3, dl, dm)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! 
! Method: 
! Spectrally truncate 3-D scalar field.
! 
! Author: 
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, J. Hack, August 1992
! Modified:          P. Worley, May 2003
! Modified:          J. Olson, Apr 2004
!
!-----------------------------------------------------------------------

   use pmgrid,   only: plon, plat

!
! Input/Output arguments
!
   real(r8), intent(inout) :: x3(plon,plev,plat)          ! Fourier -> spec. coeffs. for X
   real(r8), intent(out), optional :: dl(plon,plev,plat)  ! Spectrally trunc d(X)/d(longitude)
   real(r8), intent(out), optional :: dm(plon,plev,plat)  ! Spectrally trunc d(X)/d(latitude)
!
!---------------------------Local workspace-----------------------------
!
   real(r8), pointer :: x3_tmp(:,:,:)   ! Temporary to compute x3 of size needed for FFT
   real(r8), pointer :: dl_tmp(:,:,:)   ! Temporary to compute dl of size needed for FFT
   real(r8), pointer :: dm_tmp(:,:,:)   ! Temporary to compute dm of size needed for FFT
   real(r8) t_tmp(psp)         ! used in spectral truncation of t
   real(r8) tmp1               ! vector temporary
   real(r8) tmp2               ! vector temporary
   real(r8) tmpr               ! vector temporary (real)
   real(r8) tmpi               ! vector temporary (imaginary)
   real(r8) zwalp              ! zw*alp
   real(r8) zw                 ! w**2
#if ( ! defined USEFFTLIB )
   real(r8) work((plon+1)*plev)  ! Workspace for fft
#else
   real(r8) work((plon+1)*pcray)   ! Workspace for fft
#endif

   integer ir,ii               ! indices complex coeffs. of spec. arrs.
   integer i,k                 ! longitude, level indices
   integer irow                ! latitude pair index
   integer latm,latp           ! symmetric latitude indices
   integer lat
   integer m                   ! longitudinal wavenumber index
   integer n                   ! latitudinal wavenumber index
   integer nspec
   integer mr,mc               ! spectral indices
!
!-----------------------------------------------------------------------
!
! Transform grid -> fourier
!
   allocate( x3_tmp(plondfft,plev,plat) )
   if(present(dm)) allocate( dm_tmp(plondfft,plev,plat) )
   if(present(dl)) allocate( dl_tmp(plondfft,plev,plat) )
   do lat=1,plat
      irow = lat
      if (lat.gt.plat/2) irow = plat - lat + 1
      x3_tmp(:nlon(lat),:,lat) = x3(:nlon(lat),:,lat)
      call fft991(x3_tmp(1,1,lat),work,trig(1,irow),ifax(1,irow),1,plondfft, &
                  nlon(lat),plev,-1)
   end do                    ! lat=1,plat
!
! Loop over vertical levels
!
   do k=1,plev
!
! Zero spectral array
!
      t_tmp(:) = 0._r8
!
! Loop over latitude pairs
!
      do irow=1,plat/2
         latp = irow
         latm = plat - irow + 1
         zw = w(irow)*2._r8
!
! Multi-level field: T
!
         do i=1,2*nmmax(irow)
            tmp1 = 0.5_r8*(x3_tmp(i,k,latm) - x3_tmp(i,k,latp))
            tmp2 = 0.5_r8*(x3_tmp(i,k,latm) + x3_tmp(i,k,latp))
            x3_tmp(i,k,latm) = tmp1
            x3_tmp(i,k,latp) = tmp2
         end do
!     
! Compute tmn
!
         do m=1,nmmax(irow)
            mr = nstart(m)
            mc = 2*mr
            do n=1,nlen(m),2
               zwalp  = zw*alp (mr+n,irow)
               ir = mc + 2*n - 1
               ii = ir + 1
               t_tmp(ir) = t_tmp(ir) + zwalp*x3_tmp(2*m-1,k,latp)
               t_tmp(ii) = t_tmp(ii) + zwalp*x3_tmp(2*m  ,k,latp)
            end do
         end do

         do m=1,nmmax(irow)
            mr = nstart(m)
            mc = 2*mr
            do n=2,nlen(m),2
               zwalp  = zw*alp (mr+n,irow)
               ir = mc + 2*n - 1
               ii = ir + 1
               t_tmp(ir) = t_tmp(ir) + zwalp*x3_tmp(2*m-1,k,latm)
               t_tmp(ii) = t_tmp(ii) + zwalp*x3_tmp(2*m  ,k,latm)
            end do
         end do
      enddo                ! irow=1,plat/2
!
! Compute grid point values of:t.
!
      do irow=1,plat/2
         latp = irow
         latm = plat - irow + 1
!
! Zero fourier fields
!
         x3_tmp(:,k,latm) = 0._r8
         x3_tmp(:,k,latp) = 0._r8

         do m=1,nmmax(irow)
            mr = nstart(m)
            mc = 2*mr
            do n=1,nlen(m),2
               ir = mc + 2*n - 1
               ii = ir + 1
               tmpr = t_tmp(ir)*alp(mr+n,irow)
               tmpi = t_tmp(ii)*alp(mr+n,irow)
               x3_tmp(2*m-1,k,latm) = x3_tmp(2*m-1,k,latm) + tmpr
               x3_tmp(2*m  ,k,latm) = x3_tmp(2*m  ,k,latm) + tmpi
            end do
         end do

         do m=1,nmmax(irow)
            mr = nstart(m)
            mc = 2*mr
            do n=2,nlen(m),2
               ir = mc + 2*n - 1
               ii = ir + 1
               tmpr = t_tmp(ir)*alp(mr+n,irow)
               tmpi = t_tmp(ii)*alp(mr+n,irow)
               x3_tmp(2*m-1,k,latp) = x3_tmp(2*m-1,k,latp) + tmpr
               x3_tmp(2*m  ,k,latp) = x3_tmp(2*m  ,k,latp) + tmpi
            end do
         end do
!
! Recompute real fields from symmetric and antisymmetric parts
!
         do i=1,nlon(latm)+2
            tmp1 = x3_tmp(i,k,latm) + x3_tmp(i,k,latp)
            tmp2 = x3_tmp(i,k,latm) - x3_tmp(i,k,latp)
            x3_tmp(i,k,latm) = tmp1
            x3_tmp(i,k,latp) = tmp2
         end do
      enddo                ! irow=1,plat/2

      if(present(dl)) then
         do irow=1,plat/2
            latp = irow
            latm = plat - irow + 1
!
! Zero fourier fields
!
            dl_tmp(:,k,latm) = 0._r8
            dl_tmp(:,k,latp) = 0._r8

            do m=1,nmmax(irow)
               mr = nstart(m)
               mc = 2*mr
               do n=1,nlen(m),2
                  ir = mc + 2*n - 1
                  ii = ir + 1
                  tmpr = t_tmp(ir)*alp(mr+n,irow)
                  tmpi = t_tmp(ii)*alp(mr+n,irow)
                  dl_tmp(2*m-1,k,latm) = dl_tmp(2*m-1,k,latm) - tmpi*ra
                  dl_tmp(2*m  ,k,latm) = dl_tmp(2*m  ,k,latm) + tmpr*ra
               end do
            end do

            do m=1,nmmax(irow)
               mr = nstart(m)
               mc = 2*mr
               do n=2,nlen(m),2
                  ir = mc + 2*n - 1
                  ii = ir + 1
                  tmpr = t_tmp(ir)*alp(mr+n,irow)
                  tmpi = t_tmp(ii)*alp(mr+n,irow)
                  dl_tmp(2*m-1,k,latp) = dl_tmp(2*m-1,k,latp) - tmpi*ra
                  dl_tmp(2*m  ,k,latp) = dl_tmp(2*m  ,k,latp) + tmpr*ra
               end do
            end do
!
! d(T)/d(lamda)
!
            do m=1,nmmax(irow)
               dl_tmp(2*m-1,k,latm) = xm(m)*dl_tmp(2*m-1,k,latm)
               dl_tmp(2*m  ,k,latm) = xm(m)*dl_tmp(2*m  ,k,latm)
               dl_tmp(2*m-1,k,latp) = xm(m)*dl_tmp(2*m-1,k,latp)
               dl_tmp(2*m  ,k,latp) = xm(m)*dl_tmp(2*m  ,k,latp)
            end do
!
! Recompute real fields from symmetric and antisymmetric parts
!
            do i=1,nlon(latm)+2
               tmp1 = dl_tmp(i,k,latm) + dl_tmp(i,k,latp)
               tmp2 = dl_tmp(i,k,latm) - dl_tmp(i,k,latp)
               dl_tmp(i,k,latm) = tmp1
               dl_tmp(i,k,latp) = tmp2
            end do
         enddo                ! irow=1,plat/2
      end if

      if(present(dm)) then
         do irow=1,plat/2
            latp = irow
            latm = plat - irow + 1
!
! Zero fourier fields
!
            dm_tmp(:,k,latm) = 0._r8
            dm_tmp(:,k,latp) = 0._r8

            do m=1,nmmax(irow)
               mr = nstart(m)
               mc = 2*mr
               do n=1,nlen(m),2
                  ir = mc + 2*n - 1
                  ii = ir + 1
                  tmpr = t_tmp(ir)*dalp(mr+n,irow)
                  tmpi = t_tmp(ii)*dalp(mr+n,irow)
                  dm_tmp(2*m-1,k,latp) = dm_tmp(2*m-1,k,latp) + tmpr*ra
                  dm_tmp(2*m  ,k,latp) = dm_tmp(2*m  ,k,latp) + tmpi*ra
               end do
            end do

            do m=1,nmmax(irow)
               mr = nstart(m)
               mc = 2*mr
               do n=2,nlen(m),2
                  ir = mc + 2*n - 1
                  ii = ir + 1
                  tmpr = t_tmp(ir)*dalp(mr+n,irow)
                  tmpi = t_tmp(ii)*dalp(mr+n,irow)
                  dm_tmp(2*m-1,k,latm) = dm_tmp(2*m-1,k,latm) + tmpr*ra
                  dm_tmp(2*m  ,k,latm) = dm_tmp(2*m  ,k,latm) + tmpi*ra
               end do
            end do
!
! Recompute real fields from symmetric and antisymmetric parts
!
            do i=1,nlon(latm)+2
               tmp1 = dm_tmp(i,k,latm) + dm_tmp(i,k,latp)
               tmp2 = dm_tmp(i,k,latm) - dm_tmp(i,k,latp)
               dm_tmp(i,k,latm) = tmp1
               dm_tmp(i,k,latp) = tmp2
            end do
         enddo                ! irow=1,plat/2
      end if

   enddo                   ! k=1,plev
!
   do lat=1,plat
!     
! Transform Fourier -> grid, obtaining spectrally truncated
! grid point values.

      irow = lat
      if (lat.gt.plat/2) irow = plat - lat + 1

      call fft991(x3_tmp(1,1,lat) ,work     ,trig(1,irow),ifax(1,irow),1       , &
                  plondfft       ,nlon(lat),plev        ,+1)
      x3(:nlon(lat),:,lat) = x3_tmp(:nlon(lat),:,lat)
      if(present(dl)) then
         call fft991(dl_tmp(1,1,lat) ,work     ,trig(1,irow),ifax(1,irow),1       , &
                     plondfft       ,nlon(lat),plev        ,+1      )
         dl(:nlon(lat),:,lat) = dl_tmp(:nlon(lat),:,lat)
      end if
      if(present(dm)) then
         call fft991(dm_tmp(1,1,lat) ,work     ,trig(1,irow),ifax(1,irow),1       , &
                     plondfft       ,nlon(lat),plev        ,+1      )
         dm(:nlon(lat),:,lat) = dm_tmp(:nlon(lat),:,lat)
      end if
   end do
   deallocate( x3_tmp )
   if ( present(dm) ) deallocate( dm_tmp )
   if ( present(dl) ) deallocate( dl_tmp )

   return
end subroutine spetru_3d_scalar

!***********************************************************************

subroutine spetru_uv(u3      ,v3      ,div     ,vort    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! 
! Method: 
! Spectrally truncate U, V input fields.
! 
! Author: 
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, J. Hack, August 1992
! Modified:          P. Worley, May 2003
! Modified:          J. Olson, Apr 2004
!
!-----------------------------------------------------------------------

   use pmgrid,   only: plon, plat
   use commap,   only: rsq, cs
   use physconst,only: ez

!
! Input/Output arguments
!
   real(r8), intent(inout) :: u3(plon,plev,plat)     ! Fourier -> spec. coeffs. for u-wind
   real(r8), intent(inout) :: v3(plon,plev,plat)     ! Fourier -> spec. coeffs. for v-wind
!
! Output arguments
!
   real(r8), intent(out), optional :: div (plon,plev,plat)  ! Spectrally truncated divergence
   real(r8), intent(out), optional :: vort(plon,plev,plat)  ! Spectrally truncated vorticity

!
!---------------------------Local workspace-----------------------------
!
   real(r8), pointer :: u_cosphi(:,:,:)  ! u3*cos(phi)
   real(r8), pointer :: v_cosphi(:,:,:)  ! v3*cos(phi)
   real(r8), pointer :: div_tmp(:,:,:)   ! Temporary to compute div of size needed for FFT
   real(r8), pointer :: vort_tmp(:,:,:)  ! Temporary to compute vort of size needed for FFT
   real(r8) d_tmp(psp)         ! used in spectral truncation of div
   real(r8) vz_tmp(psp)        ! used in spectral truncation of vort
   real(r8) alpn(pspt)         ! alp*rsq*xm*ra
   real(r8) dalpn(pspt)        ! dalp*rsq*ra
   real(r8) tmp1               ! vector temporary
   real(r8) tmp2               ! vector temporary
   real(r8) tmpr               ! vector temporary (real)
   real(r8) tmpi               ! vector temporary (imaginary)
   real(r8) zcor               ! correction for absolute vorticity
   real(r8) zwalp              ! zw*alp
   real(r8) zwdalp             ! zw*dalp
   real(r8) zrcsj              ! ra/(cos**2 latitude)
   real(r8) zw                 ! w**2
#if ( ! defined USEFFTLIB )
   real(r8) work((plon+1)*plev)  ! Workspace for fft
#else
   real(r8) work((plon+1)*pcray)   ! Workspace for fft
#endif
   real(r8) zsqcs

   integer ir,ii               ! indices complex coeffs. of spec. arrs.
   integer i,k                 ! longitude, level indices
   integer irow                ! latitude pair index
   integer latm,latp           ! symmetric latitude indices
   integer lat
   integer m                   ! longitudinal wavenumber index
   integer n                   ! latitudinal wavenumber index
   integer nspec
   integer mr,mc               ! spectral indices

!
!-----------------------------------------------------------------------
!
! Compute the quantities which are transformed to spectral space:
!   1. u = u*sqrt(1-mu*mu),   u * cos(phi)
!   2. v = v*sqrt(1-mu*mu),   v * cos(phi)
!
   allocate( u_cosphi(plondfft,plev,plat) )
   allocate( v_cosphi(plondfft,plev,plat) )
   do lat=1,plat
      irow = lat
      if (lat.gt.plat/2) irow = plat - lat + 1
      zsqcs = sqrt(cs(irow))
      do k=1,plev
         do i=1,nlon(lat)
            u_cosphi(i,k,lat) = u3(i,k,lat)*zsqcs
            v_cosphi(i,k,lat) = v3(i,k,lat)*zsqcs
         end do
      end do
!
! Transform grid -> fourier
! 1st transform: U,V,T: note contiguity assumptions
! 2nd transform: LN(PS).  3rd transform: surface geopotential
!
      call fft991(u_cosphi(1,1,lat),work,trig(1,irow),ifax(1,irow),1,plondfft, &
                  nlon(lat),plev,-1)
      call fft991(v_cosphi(1,1,lat),work,trig(1,irow),ifax(1,irow),1,plondfft, &
                  nlon(lat),plev,-1)

   end do                    ! lat=1,plat
!
! Multi-level fields: U, V
!
   if ( present(div) ) allocate( div_tmp(plondfft,plev,plat) )
   if ( present(vort) ) allocate( vort_tmp(plondfft,plev,plat) )
   do k=1,plev
!
! Zero spectral arrays
!
      vz_tmp(:) = 0._r8
      d_tmp(:) = 0._r8
!
! Loop over latitude pairs
!
      do irow=1,plat/2
         latp = irow
         latm = plat - irow + 1
         zrcsj = ra/cs(irow)
         zw = w(irow)*2._r8
         do i=1,2*nmmax(irow)

            tmp1 = 0.5_r8*(u_cosphi(i,k,latm) - u_cosphi(i,k,latp))
            tmp2 = 0.5_r8*(u_cosphi(i,k,latm) + u_cosphi(i,k,latp))
            u_cosphi(i,k,latm) = tmp1
            u_cosphi(i,k,latp) = tmp2

            tmp1 = 0.5_r8*(v_cosphi(i,k,latm) - v_cosphi(i,k,latp))
            tmp2 = 0.5_r8*(v_cosphi(i,k,latm) + v_cosphi(i,k,latp))
            v_cosphi(i,k,latm) = tmp1
            v_cosphi(i,k,latp) = tmp2

         end do
!     
! Compute vzmn and dmn
!
         do m=1,nmmax(irow)
            mr = nstart(m)
            mc = 2*mr
            do n=1,nlen(m),2
               zwdalp = zw*dalp(mr+n,irow)
               zwalp  = zw*alp (mr+n,irow)
               ir = mc + 2*n - 1
               ii = ir + 1
               d_tmp(ir) = d_tmp(ir) - (zwdalp*v_cosphi(2*m-1,k,latm) + &
                  xm(m)*zwalp*u_cosphi(2*m  ,k,latp))*zrcsj
               d_tmp(ii) = d_tmp(ii) - (zwdalp*v_cosphi(2*m  ,k,latm) - &
                  xm(m)*zwalp*u_cosphi(2*m-1,k,latp))*zrcsj
               vz_tmp(ir) = vz_tmp(ir) + (zwdalp*u_cosphi(2*m-1,k,latm) - &
                  xm(m)*zwalp*v_cosphi(2*m  ,k,latp))*zrcsj
               vz_tmp(ii) = vz_tmp(ii) + (zwdalp*u_cosphi(2*m  ,k,latm) + &
                  xm(m)*zwalp*v_cosphi(2*m-1,k,latp))*zrcsj
            end do
         end do

         do m=1,nmmax(irow)
            mr = nstart(m)
            mc = 2*mr
            do n=2,nlen(m),2
               zwdalp = zw*dalp(mr+n,irow)
               zwalp  = zw*alp (mr+n,irow)
               ir = mc + 2*n - 1
               ii = ir + 1
               d_tmp(ir) = d_tmp(ir) - (zwdalp*v_cosphi(2*m-1,k,latp) + &
                  xm(m)*zwalp*u_cosphi(2*m  ,k,latm))*zrcsj
               d_tmp(ii) = d_tmp(ii) - (zwdalp*v_cosphi(2*m  ,k,latp) - &
                  xm(m)*zwalp*u_cosphi(2*m-1,k,latm))*zrcsj
               vz_tmp(ir) = vz_tmp(ir) + (zwdalp*u_cosphi(2*m-1,k,latp) - &
                  xm(m)*zwalp*v_cosphi(2*m  ,k,latm))*zrcsj
               vz_tmp(ii) = vz_tmp(ii) + (zwdalp*u_cosphi(2*m  ,k,latp) + &
                  xm(m)*zwalp*v_cosphi(2*m-1,k,latm))*zrcsj
            end do
         end do
      enddo               ! irow=1,plat/2
!
! Compute grid point values of:u,v,vz, and d.
!
      do irow=1,plat/2
         latp = irow
         latm = plat - irow + 1
         zcor = ez*alp(2,irow)
!
! Compute(u,v,vz,d)m
!
         do m=1,nmmax(irow)
            mr = nstart(m)
            do n=1,nlen(m)
!
! These statements will likely not be bfb since xm*ra is now a scalar
!
               alpn (mr+n) =  alp(mr+n,irow)*rsq(n+m-1)*xm(m)*ra
               dalpn(mr+n) = dalp(mr+n,irow)*rsq(n+m-1)      *ra
            end do
         end do
!
! Zero fourier fields
!
         u_cosphi(:,k,latm) = 0._r8
         u_cosphi(:,k,latp) = 0._r8

         v_cosphi(:,k,latm) = 0._r8
         v_cosphi(:,k,latp) = 0._r8

         if(present(vort)) then
            vort_tmp(:,k,latm) = 0._r8
            vort_tmp(:,k,latp) = 0._r8
         end if

         if(present(div)) then
            div_tmp(:,k,latm) = 0._r8
            div_tmp(:,k,latp) = 0._r8
         end if

         do m=1,nmmax(irow)
            mr = nstart(m)
            mc = 2*mr
            do n=1,nlen(m),2
               ir = mc + 2*n - 1
               ii = ir + 1
!
               tmpr = d_tmp(ir)*alpn(mr+n)
               tmpi = d_tmp(ii)*alpn(mr+n)
               u_cosphi(2*m-1,k,latm) = u_cosphi(2*m-1,k,latm) + tmpi
               u_cosphi(2*m  ,k,latm) = u_cosphi(2*m  ,k,latm) - tmpr
!
               tmpr = d_tmp(ir)*dalpn(mr+n)
               tmpi = d_tmp(ii)*dalpn(mr+n)
               v_cosphi(2*m-1,k,latp) = v_cosphi(2*m-1,k,latp) - tmpr
               v_cosphi(2*m  ,k,latp) = v_cosphi(2*m  ,k,latp) - tmpi
!
               tmpr = vz_tmp(ir)*dalpn(mr+n)
               tmpi = vz_tmp(ii)*dalpn(mr+n)
               u_cosphi(2*m-1,k,latp) = u_cosphi(2*m-1,k,latp) + tmpr
               u_cosphi(2*m  ,k,latp) = u_cosphi(2*m  ,k,latp) + tmpi
!
               tmpr = vz_tmp(ir)*alpn(mr+n)
               tmpi = vz_tmp(ii)*alpn(mr+n)
               v_cosphi(2*m-1,k,latm) = v_cosphi(2*m-1,k,latm) + tmpi
               v_cosphi(2*m  ,k,latm) = v_cosphi(2*m  ,k,latm) - tmpr
!
               if(present(div)) then
                  tmpr = d_tmp(ir)*alp(mr+n,irow)
                  tmpi = d_tmp(ii)*alp(mr+n,irow)
                  div_tmp(2*m-1,k,latm) = div_tmp(2*m-1,k,latm) + tmpr
                  div_tmp(2*m  ,k,latm) = div_tmp(2*m  ,k,latm) + tmpi
               end if
!
               if(present(vort)) then
                  tmpr = vz_tmp(ir)*alp(mr+n,irow)
                  tmpi = vz_tmp(ii)*alp(mr+n,irow)
                  vort_tmp(2*m-1,k,latm) = vort_tmp(2*m-1,k,latm) + tmpr
                  vort_tmp(2*m  ,k,latm) = vort_tmp(2*m  ,k,latm) + tmpi
               end if
            end do
         end do

         do m=1,nmmax(irow)
            mr = nstart(m)
            mc = 2*mr
            do n=2,nlen(m),2
               ir = mc + 2*n - 1
               ii = ir + 1
!     
               tmpr = d_tmp(ir)*alpn(mr+n)
               tmpi = d_tmp(ii)*alpn(mr+n)
               u_cosphi(2*m-1,k,latp) = u_cosphi(2*m-1,k,latp) + tmpi
               u_cosphi(2*m  ,k,latp) = u_cosphi(2*m  ,k,latp) - tmpr
!
               tmpr = d_tmp(ir)*dalpn(mr+n)
               tmpi = d_tmp(ii)*dalpn(mr+n)
               v_cosphi(2*m-1,k,latm) = v_cosphi(2*m-1,k,latm) - tmpr
               v_cosphi(2*m  ,k,latm) = v_cosphi(2*m  ,k,latm) - tmpi
!
               tmpr = vz_tmp(ir)*dalpn(mr+n)
               tmpi = vz_tmp(ii)*dalpn(mr+n)
               u_cosphi(2*m-1,k,latm) = u_cosphi(2*m-1,k,latm) + tmpr
               u_cosphi(2*m  ,k,latm) = u_cosphi(2*m  ,k,latm) + tmpi
!
               tmpr = vz_tmp(ir)*alpn(mr+n)
               tmpi = vz_tmp(ii)*alpn(mr+n)
               v_cosphi(2*m-1,k,latp) = v_cosphi(2*m-1,k,latp) + tmpi
               v_cosphi(2*m  ,k,latp) = v_cosphi(2*m  ,k,latp) - tmpr
!
               if(present(div)) then
                  tmpr = d_tmp(ir)*alp(mr+n,irow)
                  tmpi = d_tmp(ii)*alp(mr+n,irow)
                  div_tmp(2*m-1,k,latp) = div_tmp(2*m-1,k,latp) + tmpr
                  div_tmp(2*m  ,k,latp) = div_tmp(2*m  ,k,latp) + tmpi
               end if
!
               if(present(vort)) then
                  tmpr = vz_tmp(ir)*alp(mr+n,irow)
                  tmpi = vz_tmp(ii)*alp(mr+n,irow)
                  vort_tmp(2*m-1,k,latp) = vort_tmp(2*m-1,k,latp) + tmpr
                  vort_tmp(2*m  ,k,latp) = vort_tmp(2*m  ,k,latp) + tmpi
               end if
            end do
         end do
!
! Correction to get the absolute vorticity.
!     
         if(present(vort)) then
            vort_tmp(1,k,latp) = vort_tmp(1,k,latp) + zcor
         end if
!
! Recompute real fields from symmetric and antisymmetric parts
!
         do i=1,nlon(latm)+2
            tmp1 = u_cosphi(i,k,latm) + u_cosphi(i,k,latp)
            tmp2 = u_cosphi(i,k,latm) - u_cosphi(i,k,latp)
            u_cosphi(i,k,latm) = tmp1
            u_cosphi(i,k,latp) = tmp2
!
            tmp1 = v_cosphi(i,k,latm) + v_cosphi(i,k,latp)
            tmp2 = v_cosphi(i,k,latm) - v_cosphi(i,k,latp)
            v_cosphi(i,k,latm) = tmp1
            v_cosphi(i,k,latp) = tmp2
!
            if(present(vort)) then
               tmp1 = vort_tmp(i,k,latm) + vort_tmp(i,k,latp)
               tmp2 = vort_tmp(i,k,latm) - vort_tmp(i,k,latp)
               vort_tmp(i,k,latm) = tmp1
               vort_tmp(i,k,latp) = tmp2
            end if
!
            if(present(div)) then
               tmp1 = div_tmp(i,k,latm) + div_tmp(i,k,latp)
               tmp2 = div_tmp(i,k,latm) - div_tmp(i,k,latp)
               div_tmp(i,k,latm) = tmp1
               div_tmp(i,k,latp) = tmp2
            end if
         end do
      enddo               ! irow=1,plat/2
   enddo                  ! k=1,plev
!
   do lat=1,plat
!     
! Transform Fourier -> grid, obtaining spectrally truncated
! grid point values.
!
      irow = lat
      if (lat.gt.plat/2) irow = plat - lat + 1

      call fft991(u_cosphi(1,1,lat),work,trig(1,irow),ifax(1,irow),1,plondfft, &
                  nlon(lat),plev,+1)
      call fft991(v_cosphi(1,1,lat),work,trig(1,irow),ifax(1,irow),1,plondfft, &
                  nlon(lat),plev,+1)
      if(present(vort)) then
         call fft991(vort_tmp(1,1,lat),work,trig(1,irow),ifax(1,irow),1, &
                     plondfft,nlon(lat),plev,+1)
         vort(:nlon(lat),:,lat) = vort_tmp(:nlon(lat),:,lat)
      end if
      if(present(div)) then
         call fft991(div_tmp(1,1,lat),work,trig(1,irow),ifax(1,irow),1,plondfft, &
                     nlon(lat),plev,+1)
         div(:nlon(lat),:,lat) = div_tmp(:nlon(lat),:,lat)
      end if
!
! Convert U,V to u,v
!
      zsqcs = sqrt(cs(irow))
      do k=1,plev
         do i=1,nlon(lat)
            u3(i,k,lat) = u_cosphi(i,k,lat)/zsqcs
            v3(i,k,lat) = v_cosphi(i,k,lat)/zsqcs
         end do
      end do
   enddo
   deallocate( u_cosphi )
   deallocate( v_cosphi )
   if ( present(div) ) deallocate( div_tmp )
   if ( present(vort) ) deallocate( vort_tmp )

   return
end subroutine spetru_uv

end module spetru
