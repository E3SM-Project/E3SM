
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Reallocation routines for the Fourier coefficients
! 
! Method: 
!   1) After FFT preceding Legendre analysis, reallocate fftbuf
!      to decompose over wavenumber, recombining latitudes.
!   2) Before FFT following Legendre synthesis, reallocate fftbuf
!      to recombine wavenumbers, decomposing over latitude.
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------
subroutine realloc4a(nlon_fft_in, nlon_fft_out, fftbuf_in, fftbuf_out )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Reallocation routines for the Fourier coefficients
! 
! Method: 
!   After FFT preceding Legendre analysis, reallocate fftbuf
!   to decompose over wavenumber, combining latitudes.
! 
! Author: 
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, Oct 1995
!                    J. Truesdale, Feb. 1996
! Modified:          P. Worley, September 2002, December 2003, 
!                               October 2004, April 2007
! 
!-----------------------------------------------------------------------

#ifdef SPMD

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use comspe
   use spmd_dyn
   use mpishorthand
   use spmd_utils, only : iam, npes, altalltoallv
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <comsta.h>
!------------------------------Parameters-------------------------------
!
  integer, parameter :: msgtag  = 1000
!---------------------------Input arguments-----------------------------
!
   integer, intent(in) :: nlon_fft_in      ! first dimension of input array
   integer, intent(in) :: nlon_fft_out     ! first dimension of output array
   real(r8), intent(in)  :: fftbuf_in(nlon_fft_in,9,plev,beglat:endlat) 
                            ! buffer used for in-place FFTs
   real(r8), intent(out) :: fftbuf_out(nlon_fft_out,9,plev,plat) 
                            ! buffer used for reordered Fourier coefficients
!
!---------------------------Local workspace-----------------------------
!
! xxx_l: local decomposition
! xxx_r: remote decomposition
   integer :: procid
   integer :: length_r, length_l
   integer :: bpos
   integer :: step, ifld, k, i
   integer :: lat_l, lat_r, beglat_r, endlat_r
!
   logical, save :: first = .true.
   integer, allocatable, save :: sndcnts(:), sdispls(:)
   integer, allocatable, save :: rcvcnts(:), rdispls(:)
   integer, allocatable, save :: sndcnts_act(:), sdispls_act(:)
   integer, allocatable, save :: rcvcnts_act(:), rdispls_act(:)
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
      sndcnts(:) = 0
      do step=1,realloc4_steps
         procid = realloc4_proc(step)
         length_r = 2*numm(procid)
         sndcnts(procid) = length_r*(plev*8 + 1)*numlats
      enddo
!   
      sdispls(0) = 0
      do procid=1,npes-1
        sdispls(procid) = sdispls(procid-1) + sndcnts(procid-1)
      enddo
!
      length_l = 2*numm(iam)
      rcvcnts(:) = 0
      do step=1,realloc4_steps
         procid = realloc4_proc(step)
         rcvcnts(procid) = length_l*(plev*8 + 1)*nlat_p(procid)
      enddo
!   
      rdispls(0) = 0
      do procid=1,npes-1
        rdispls(procid) = rdispls(procid-1) + rcvcnts(procid-1)
      enddo
!
      pdispls(:) = 0
      call mpialltoallint(rdispls, 1, pdispls, 1, mpicom)
!
      allocate(sndcnts_act(0:dyn_npes-1))
      allocate(sdispls_act(0:dyn_npes-1))
      allocate(rcvcnts_act(0:dyn_npes-1))
      allocate(rdispls_act(0:dyn_npes-1))
!   
      do procid=0,dyn_npes-1
         sndcnts_act(procid) = sndcnts(procid*dyn_npes_stride)
         sdispls_act(procid) = sdispls(procid*dyn_npes_stride)
      enddo
!   
      do procid=0,dyn_npes-1
         rcvcnts_act(procid) = rcvcnts(procid*dyn_npes_stride)
         rdispls_act(procid) = rdispls(procid*dyn_npes_stride)
      enddo
!
      first = .false.
   endif
!     
! Copy local data to new location
   length_l = 2*numm(iam)
   do lat_l=beglat,endlat
!$OMP PARALLEL DO PRIVATE(K, IFLD, I)
      do k=1,plev
         do ifld=1,8
            do i=1,length_l
               fftbuf_out(i,ifld,k,lat_l) = fftbuf_in(locrm(i,iam),ifld,k,lat_l)
            enddo
         enddo
      enddo
      do i=1,length_l
         fftbuf_out(i,9,1,lat_l) = fftbuf_in(locrm(i,iam),9,1,lat_l)
      enddo
   enddo
!
! Fill message buffer
!$OMP PARALLEL DO PRIVATE (STEP, PROCID, LENGTH_R, BPOS, LAT_L, IFLD, K, I)
   do step=1,realloc4_steps
      procid = realloc4_proc(step)
      length_r = 2*numm(procid)
!
      bpos = sdispls(procid)
      do lat_l=beglat,endlat
         do k=1,plev
            do ifld=1,8
               do i=1,length_r
                  buf1(bpos+i) = fftbuf_in(locrm(i,procid),ifld,k,lat_l)
               enddo
               bpos = bpos+length_r
            enddo
         enddo
         do i=1,length_r
            buf1(bpos+i) = fftbuf_in(locrm(i,procid),9,1,lat_l)
         enddo
         bpos = bpos+length_r
      enddo
   enddo
!
! Get remote data
!
   if (dyn_alltoall .eq. 0) then
      if (beglat <= endlat) then
         call mpialltoallv(buf1, sndcnts_act, sdispls_act, mpir8, &
                           buf2, rcvcnts_act, rdispls_act, mpir8, &
                           mpicom_dyn_active)
      endif
   else
      call altalltoallv(dyn_alltoall, iam, npes,     &
                        realloc4_steps, realloc4_proc, &
                        buf1, spmdbuf_siz, sndcnts, sdispls, mpir8, &
                        buf2, spmdbuf_siz, rcvcnts, rdispls, mpir8, &
                        msgtag, pdispls, mpir8, buf2win, mpicom)
   endif
!
! Copy out of message buffers
!
!$OMP PARALLEL DO PRIVATE (STEP, PROCID, BEGLAT_R, ENDLAT_R, BPOS, LAT_R, IFLD, K, I)
   do step=1,realloc4_steps
      procid = realloc4_proc(step)
      beglat_r = cut(1,procid)
      endlat_r = cut(2,procid)
      bpos = rdispls(procid)
      do lat_r=beglat_r,endlat_r
         do k=1,plev
            do ifld=1,8
               do i=1,length_l
                  fftbuf_out(i,ifld,k,lat_r) = buf2(bpos+i)
               enddo
               bpos = bpos+length_l
            enddo
         enddo
         do i=1,length_l
            fftbuf_out(i,9,1,lat_r) = buf2(bpos+i)
         enddo
         bpos = bpos+length_l
      enddo
!
   end do
#endif
   return
   end subroutine realloc4a

subroutine realloc4b(nlon_fft_in, nlon_fft_out, fftbuf_in, fftbuf_out )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Reallocation routines for the Fourier coefficients
! 
! Method: 
!   Before FFT following Legendre synthesis, reallocate fftbuf
!   to combine wavenumbers, decomposing over latitude.
! 
! Author:  P. Worley, September 2002
! Modified: P. Worley, December 2003, October 2004
! 
!-----------------------------------------------------------------------

#ifdef SPMD

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use comspe
   use spmd_dyn
   use mpishorthand
   use spmd_utils, only : iam, npes, altalltoallv

!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <comsta.h>
!------------------------------Parameters-------------------------------
!
  integer, parameter :: msgtag  = 2000
!---------------------------Input arguments--------------------------
!
   integer, intent(in) :: nlon_fft_in      ! first dimension of input array
   integer, intent(in) :: nlon_fft_out     ! first dimension of output array
   real(r8), intent(in)  :: fftbuf_in(nlon_fft_in,8,plevp,plat) 
                            ! buffer of Fourier coefficients to be reordered
   real(r8), intent(out) :: fftbuf_out(nlon_fft_out,8,plevp,beglat:endlat) 
                            ! buffer used for in-place FFTs
!
!---------------------------Local workspace-----------------------------
!
! xxx_l: local decomposition
! xxx_r: remote decomposition
   integer :: procid
   integer :: length_r, length_l
   integer :: bpos
   integer :: step, ifld, k, i
   integer :: lat_l, lat_r
   integer :: beglat_r, endlat_r
!
   logical, save :: first = .true.
   integer, allocatable, save :: sndcnts(:), sdispls(:)
   integer, allocatable, save :: rcvcnts(:), rdispls(:)
   integer, allocatable, save :: sndcnts_act(:), sdispls_act(:)
   integer, allocatable, save :: rcvcnts_act(:), rdispls_act(:)
   integer, allocatable, save :: pdispls(:)
!-----------------------------------------------------------------------
   if (first) then
! Compute send/recv counts and displacements
      allocate(sndcnts(0:npes-1))
      allocate(sdispls(0:npes-1))
      allocate(rcvcnts(0:npes-1))
      allocate(rdispls(0:npes-1))
      allocate(pdispls(0:npes-1))
!
      length_l = 2*numm(iam)
      sndcnts(:) = 0
      do step=1,realloc4_steps
         procid = realloc4_proc(step)
         sndcnts(procid) = length_l*(8*plev + 4)*nlat_p(procid)
      enddo
!   
      sdispls(0) = 0
      do procid=1,npes-1
        sdispls(procid) = sdispls(procid-1) + sndcnts(procid-1)
      enddo
!
      rcvcnts(:) = 0
      do step=1,realloc4_steps
         procid = realloc4_proc(step)
         length_r = 2*numm(procid)
         rcvcnts(procid) = length_r*(8*plev + 4)*numlats
      enddo
!   
      rdispls(0) = 0
      do procid=1,npes-1
        rdispls(procid) = rdispls(procid-1) + rcvcnts(procid-1)
      enddo
!
      pdispls(:) = 0
      call mpialltoallint(rdispls, 1, pdispls, 1, mpicom)
!
      allocate(sndcnts_act(0:dyn_npes-1))
      allocate(sdispls_act(0:dyn_npes-1))
      allocate(rcvcnts_act(0:dyn_npes-1))
      allocate(rdispls_act(0:dyn_npes-1))
!   
      do procid=0,dyn_npes-1
         sndcnts_act(procid) = sndcnts(procid*dyn_npes_stride)
         sdispls_act(procid) = sdispls(procid*dyn_npes_stride)
      enddo
!   
      do procid=0,dyn_npes-1
         rcvcnts_act(procid) = rcvcnts(procid*dyn_npes_stride)
         rdispls_act(procid) = rdispls(procid*dyn_npes_stride)
      enddo
!
      first = .false.
   endif
!
! Copy local data to new location
   length_l = 2*numm(iam)
   do lat_l=beglat,endlat
!$OMP PARALLEL DO PRIVATE(K, IFLD, I)
      do k=1,plev
         do ifld=1,8
            do i=1,length_l
               fftbuf_out(locrm(i,iam),ifld,k,lat_l) = fftbuf_in(i,ifld,k,lat_l)
            enddo
         enddo
      enddo
!
!$OMP PARALLEL DO PRIVATE(IFLD, I)
      do ifld=1,4
         do i=1,length_l
            fftbuf_out(locrm(i,iam),ifld,plevp,lat_l) = fftbuf_in(i,ifld,plevp,lat_l)
         enddo
      enddo
   enddo
!
! Fill message buffer
!$OMP PARALLEL DO PRIVATE (STEP, PROCID, BEGLAT_R, ENDLAT_R, BPOS, LAT_R, K, IFLD, I)
   do step=1,realloc4_steps
      procid = realloc4_proc(step)
      beglat_r = cut(1,procid)
      endlat_r = cut(2,procid)
      bpos = sdispls(procid)
!
      do lat_r=beglat_r,endlat_r
         do k=1,plev
            do ifld=1,8
               do i=1,length_l
                  buf1(bpos+i) = fftbuf_in(i,ifld,k,lat_r)
               enddo
               bpos = bpos+length_l
            enddo
         enddo
         do ifld=1,4
            do i=1,length_l
               buf1(bpos+i) = fftbuf_in(i,ifld,plevp,lat_r)
            enddo
            bpos = bpos+length_l
         enddo
      enddo
   enddo
!
! Get remote data
!
   if (dyn_alltoall .eq. 0) then
      if (beglat <= endlat) then
         call mpialltoallv(buf1, sndcnts_act, sdispls_act, mpir8, &
                           buf2, rcvcnts_act, rdispls_act, mpir8, &
                           mpicom_dyn_active)
      endif
   else
      call altalltoallv(dyn_alltoall, iam, npes,       &
                        realloc4_steps, realloc4_proc, &
                        buf1, spmdbuf_siz, sndcnts, sdispls, mpir8, &
                        buf2, spmdbuf_siz, rcvcnts, rdispls, mpir8, &
                        msgtag, pdispls, mpir8, buf2win, mpicom)
   endif
!
! Copy out of message buffers
!
!$OMP PARALLEL DO PRIVATE (STEP, PROCID, LENGTH_R, BPOS, LAT_L, K, IFLD, I)
   do step=1,realloc4_steps
      procid = realloc4_proc(step)
      length_r = 2*numm(procid)
      bpos = rdispls(procid)

      do lat_l=beglat,endlat
         do k=1,plev
            do ifld=1,8
               do i=1,length_r
                  fftbuf_out(locrm(i,procid),ifld,k,lat_l) = buf2(bpos+i)
               enddo
               bpos = bpos+length_r
            enddo
         enddo

         do ifld=1,4
            do i=1,length_r
               fftbuf_out(locrm(i,procid),ifld,plevp,lat_l) = buf2(bpos+i)
            enddo
            bpos = bpos+length_r
         enddo

      enddo
!
   end do
#endif
   return
   end subroutine realloc4b

