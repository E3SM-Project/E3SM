
subroutine realloc7 (vmax2d, vmax2dt, vcour)

!----------------------------------------------------------------------- 
! 
! Purpose: Reallocation routine for energy and log stats
! 
! Method: MPI_Allgatherv (or point-to-point implementation)
! 
! Author: J. Rosinski
! Modified: P. Worley, September 2002, December 2003, October 2004
! 
!-----------------------------------------------------------------------

#ifdef SPMD
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plat, plev, numlats, beglat, endlat
   use mpishorthand
   use spmd_dyn
   use spmd_utils, only : iam, npes, altalltoallv
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <comsta.h>
!------------------------------Parameters-------------------------------
!
   integer, parameter :: msgtag  = 3000
!---------------------------Input arguments-----------------------------
!
   real(r8), intent(inout) :: vmax2d(plev,plat)   ! Max. wind at each lvl, lat
   real(r8), intent(inout) :: vmax2dt(plev,plat)  ! Max. truncated wind at each lvl, lat
   real(r8), intent(inout) :: vcour(plev,plat)    ! Max. Courant number at each lvl, lat
!
!---------------------------Local workspace-----------------------------
!
   integer procid
   integer bufpos
   integer procj
   integer step, j, k, jstrt
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
      sndcnt = (plev*3 + 5)*numlats
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
         rcvcnts(procid) = (plev*3 + 5)*nlat_p(procid)
      enddo
      rcvcnts(iam) = (plev*3 + 5)*numlats
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
! psurf
   do j=1,numlats
      buf1(bufpos+j) = psurf(jstrt+j)
   enddo
   bufpos = bufpos + numlats
! stq
   do j=1,numlats
      buf1(bufpos+j) = stq(jstrt+j)
   enddo
   bufpos = bufpos + numlats
! rmst
   do j=1,numlats
      buf1(bufpos+j) = rmst(jstrt+j)
   enddo
   bufpos = bufpos + numlats
! rmsd
   do j=1,numlats
      buf1(bufpos+j) = rmsd(jstrt+j)
   enddo
   bufpos = bufpos + numlats
! rmsz
   do j=1,numlats
      buf1(bufpos+j) = rmsz(jstrt+j)
   enddo
   bufpos = bufpos + numlats
!vmax2d
   do j=beglat,endlat
      do k=1,plev
         buf1(bufpos+k) = vmax2d(k,j)
      enddo
      bufpos = bufpos + plev
   enddo
! vmax2dt
   do j=beglat,endlat
      do k=1,plev
         buf1(bufpos+k) = vmax2dt(k,j)
      enddo
      bufpos = bufpos + plev
   enddo
! vcour
   do j=beglat,endlat
      do k=1,plev
         buf1(bufpos+k) = vcour(k,j)
      enddo
      bufpos = bufpos + plev
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
!$OMP PARALLEL DO PRIVATE (STEP, PROCID, BEGLAT_P, ENDLAT_P, NUMLATS_P, JSTRT_P, BUFPOS, J, K)
   do step=1,allgather_steps
      procid = allgather_proc(step)
      beglat_p = cut(1,procid)
      endlat_p = cut(2,procid)
      numlats_p = nlat_p(procid)
      bufpos = rdispls(procid)
! psurf
      jstrt_p  = beglat_p - 1
      do j=1,numlats_p
         psurf(jstrt_p+j) = buf2(bufpos+j)
      enddo
      bufpos = bufpos + numlats_p
! stq
      do j=1,numlats_p
         stq(jstrt_p+j) = buf2(bufpos+j)
      enddo
      bufpos = bufpos + numlats_p
! rmst
      do j=1,numlats_p
         rmst(jstrt_p+j) = buf2(bufpos+j)
      enddo
      bufpos = bufpos + numlats_p
! rmsd
      do j=1,numlats_p
         rmsd(jstrt_p+j) = buf2(bufpos+j) 
      enddo
      bufpos = bufpos + numlats_p
! rmsz
      do j=1,numlats_p
         rmsz(jstrt_p+j) = buf2(bufpos+j) 
      enddo
      bufpos = bufpos + numlats_p
! vmax2d
      do j=beglat_p,endlat_p
         do k=1,plev
            vmax2d(k,j) = buf2(bufpos+k)
         enddo
         bufpos = bufpos + plev
      enddo
! vmax2dt
      do j=beglat_p,endlat_p
         do k=1,plev
            vmax2dt(k,j) = buf2(bufpos+k)
         enddo
         bufpos = bufpos + plev
      enddo
! vcour
      do j=beglat_p,endlat_p
         do k=1,plev
            vcour(k,j) = buf2(bufpos+k)
         enddo
         bufpos = bufpos + plev
      enddo
!
   enddo
#endif
   return
end subroutine realloc7

