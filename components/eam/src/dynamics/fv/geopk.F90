!-----------------------------------------------------------------------
!BOP
! !ROUTINE: geopk --- Calculate geopotential to the kappa
!
!-----------------------------------------------------------------------
! There are three versions of geopk below. The first is the standard
! version and is typically used with transposes between yz and xy
! space. The second (called geopk16) operates in yz space and performs
! semi-global communication in the z direction (to avoid transposes).
! It also can use 16-byte reals to preserve accuracy through round-off;
! this is accomplished by toggling DSIZE to 16 immediately below.
! The third version (called geopk_d) also operates in yz space
! and implements a ring-pipeline algorithm in the z direction.
! Numerics are identical with the first version without requiring
! 16-byte arithmetic. While less parallel, communication costs are
! smaller, and this is often the fastest option.
!
! Note that the interfaces to the first, second, and third versions are 
! slightly different. Also, geopk (the standard version with transposes) 
! is called for the D-grid during the last two small timesteps in cd_core.
! Geopk16 uses mod_comm communication calls; one can activate the old
! Pilgrim calls (for debugging) by activating PaREXCH immediately below.

!#define PAREXCH
!#define DSIZE 16
#define DSIZE 8

#if (DSIZE == 16)
# define DTWO 2
#else
# define DTWO 1
#endif
!-----------------------------------------------------------------------
!
! !INTERFACE:
      subroutine geopk(grid, pe, delp, pk, wz, hs, pt, cp, akap, nx)

      use shr_kind_mod, only: r8 => shr_kind_r8
      use dynamics_vars, only: T_FVDYCORE_GRID

      implicit none

! !INPUT PARAMETERS:

      type (T_FVDYCORE_GRID), intent(in) :: grid
      integer nx                        ! # of pieces in longitude direction
      real(r8)    akap, cp
      real(r8)    hs(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy)
      real(r8)    pt(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)
      real(r8)  delp(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)

! !OUTPUT PARAMETERS:
      real(r8) wz(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km+1)  ! space N*1 S*1
      real(r8) pk(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km+1)  ! space N*1 S*1
      real(r8) pe(grid%ifirstxy:grid%ilastxy,grid%km+1,grid%jfirstxy:grid%jlastxy)

! !DESCRIPTION:
!     Calculates geopotential and pressure to the kappa.  This is an expensive
!     operation and several out arrays are kept around for future use.
!
! !REVISION HISTORY:
!
!  WS  99.10.22: MPIed SJ's original SMP version
!  SJL 00.01.01: Merged C-core and D-core computation
!                SMP "decmposition" in E-W by combining i and j loops
!  WS  00.12.01: Replaced MPI_ON with SPMD; hs now distributed
!  AAM 01.06.27: Generalize for 2D decomposition
!  AAM 01.07.24: Removed dpcheck
!  WS  04.10.07: Simplified interface using Grid as input argument
!  WS  05.05.25: Merged CAM and GEOS5 versions (mostly CAM)
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local:
      real(r8), parameter ::  D0_0                    =  0.0_r8
      integer :: im, jm, km, jfirst, jlast, ifirst, ilast
      real(r8) :: ptop

      integer i, j, k
      integer ixj, jp, it, i1, i2, nxu, itot
      real(r8) delpp(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)

      ptop = grid%ptop
      im = grid%im
      jm = grid%jm
      km = grid%km
      ifirst = grid%ifirstxy
      ilast  = grid%ilastxy
      jfirst = grid%jfirstxy
      jlast  = grid%jlastxy

      itot = ilast - ifirst + 1
!     nxu = nx
      nxu = 1
      it = itot / nxu
      jp = nxu * ( jlast - jfirst + 1 )

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i1, i2, ixj, i, j, k )

!     do 2000 j=jfirst,jlast
      do 2000 ixj=1, jp

         j  = jfirst + (ixj-1)/nxu
         i1 = ifirst + it * mod(ixj-1, nxu)
         i2 = i1 + it - 1

         do i=i1,i2
            pe(i,1,j) = D0_0
            wz(i,j,km+1) = D0_0
         enddo

! Top down
         do k=2,km+1
            do i=i1,i2
               pe(i,k,j)  = pe(i,k-1,j) + delp(i,j,k-1)
            enddo
         enddo
         do k=1,km+1
            do i=i1,i2
               pe(i,k,j)  = pe(i,k,j) + ptop
               pk(i,j,k) = pe(i,k,j)**akap
            enddo
         enddo

! Bottom up
         do k=1,km
            do i=i1,i2
               delpp(i,j,k) = cp*pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))
            enddo
         enddo
         do k=km,1,-1
            do i=i1,i2
               wz(i,j,k) = wz(i,j,k+1)+delpp(i,j,k)
            enddo
         enddo
         do k=1,km+1
            do i=i1,i2
               wz(i,j,k) = wz(i,j,k)+hs(i,j)
            enddo
         enddo
2000  continue

      return
!EOC
      end subroutine geopk
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !ROUTINE: geopk16 --- Calculate geopotential to the kappa
!
! !INTERFACE:
      subroutine geopk16(grid, pe, delp, pk, wz, hs, pt, ng, cp, akap )

      use shr_kind_mod,  only : r8 => shr_kind_r8, i8 => shr_kind_i8
      use decompmodule,  only : decomptype
      use dynamics_vars, only : T_FVDYCORE_GRID

#if defined( SPMD )
      use parutilitiesmodule, only : parexchangevector
      use mod_comm, only : blockdescriptor, get_partneroffset,      &
                           mp_sendirr, mp_recvirr, max_nparcels
      use spmd_dyn, only: npes_yz
#endif

      implicit none

#if defined ( SPMD )
#include "mpif.h"
#endif

! !INPUT PARAMETERS:

      type (T_FVDYCORE_GRID), intent(in) :: grid
      integer, intent(in)  :: ng      ! Halo size (not always = ng_d)

      real(r8)    akap, cp
      real(r8)    hs(1:grid%im,grid%jfirst:grid%jlast)

! !INPUT PARAMETERS:
      real(r8)    pt(1:grid%im,grid%jfirst-ng:grid%jlast+ng,grid%kfirst:grid%klast) 
      real(r8)  delp(1:grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)

! !OUTPUT PARAMETERS:
      real(r8) wz(1:grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast+1)  ! space N*1 S*1
      real(r8) pk(1:grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast+1)  ! space N*1 S*1
      real(r8) pe(1:grid%im,grid%kfirst:grid%klast+1,grid%jfirst:grid%jlast)  ! temporary variable

! !DESCRIPTION:
!     Calculates geopotential and pressure to the kappa.  This is an expensive
!     operation and several out arrays are kept around for future use.
!     To preserve accuracy through round-off, 16-byte reals are used
!     for some intermediate data.
!
! !REVISION HISTORY:
!
!  AAM 00.12.18: Original version
!  AAM 03.01.21: Use mod_comm
!  WS  03.11.19: Merged latest CAM version (by AAM)
!  WS  04.10.07: Simplified interface using Grid as input argument
!  WS  05.05.17: Merged CAM and GEOS5 versions
!
!EOP
!---------------------------------------------------------------------
!BOC

#ifndef NO_CRAY_POINTERS

! Local:
      integer :: i, j, k, nk, ijtot, ierror, ione

      integer :: im,jm,km, ifirst, ilast, jfirst, jlast, kfirst, klast
      real(r8):: ptop

      integer :: npr_y, npr_z, myid_y, myid_z
      integer :: twod_decomp, mod_geopk

#if (DSIZE == 16)
#ifdef NO_R16
      integer,parameter :: r16= selected_real_kind(12) ! 8 byte real
#else
      integer,parameter :: r16= selected_real_kind(24) ! 16 byte real
#endif
      real(r16), parameter ::  DP0_0                    =  0.0_r16
      real(r16)  delp16(1:grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)
      real(r16)  pe16(1:grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast+1)
      real(r16)  inbuf(1:grid%im,grid%jfirst:grid%jlast,0:grid%npr_z-1)
      real(r16)  outbuf(1:grid%im,grid%jfirst:grid%jlast,0:grid%npr_z-1)
#else
      real (r8), parameter ::  DP0_0                    =  0.0_r8
      real (r8) delp16(1:grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)
      real (r8) pe16(1:grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast+1)
      real (r8) inbuf(1:grid%im,grid%jfirst:grid%jlast,0:grid%npr_z-1)
      real (r8) outbuf(1:grid%im,grid%jfirst:grid%jlast,0:grid%npr_z-1)
#endif
      integer sendcount(0:grid%npr_z-1), recvcount(0:grid%npr_z-1)

#if defined(SPMD)
!
! data structures for mp_sendirr, mp_recvirr
!
      type (blockdescriptor), allocatable, save :: sendbl1(:), recvbl1(:)
      type (blockdescriptor), allocatable, save :: sendbl2(:), recvbl2(:)

#endif

      integer first_time_through
      data first_time_through / 0 /

! Arrays inbuf8 and outbuf8 are created to fool the compiler
!  into accepting them as calling arguments for parexchangevector.
!  The trickery below equivalences them to inbuf and outbuf
      real (r8) inbuf8(1), outbuf8(1)
      pointer (ptr_inbuf8, inbuf8)
      pointer (ptr_outbuf8, outbuf8)
      integer (i8) locinbuf, locoutbuf

!
! Initialize variables from Grid
!
      ptop = grid%ptop

      im       = grid%im
      jm       = grid%jm
      km       = grid%km

      ifirst   = 1               ! 2004.10.04 (WS): Now hardwired for 1..im
      ilast    = grid%im         ! Code was always used in this mode 
      jfirst   = grid%jfirst
      jlast    = grid%jlast
      kfirst   = grid%kfirst
      klast    = grid%klast

      myid_y = grid%myid_y
      myid_z = grid%myid_z

      npr_y = grid%npr_y
      npr_z = grid%npr_z

      twod_decomp = grid%twod_decomp
      mod_geopk   = grid%mod_geopk

      ijtot = (jlast-jfirst+1) * (ilast-ifirst+1)

#if defined (SPMD)
      if (first_time_through .eq. 0) then
       first_time_through = 1
       ione = 1
       if (npr_z .gt. 1) then
        allocate( sendbl1(0:npes_yz-1) )
        allocate( recvbl1(0:npes_yz-1) )
        allocate( sendbl2(0:npes_yz-1) )
        allocate( recvbl2(0:npes_yz-1) )

        do nk = 0,npes_yz-1

          sendbl1(nk)%method = mod_geopk
          sendbl2(nk)%method = mod_geopk
          recvbl1(nk)%method = mod_geopk
          recvbl2(nk)%method = mod_geopk

! allocate for either method (safety)
          allocate( sendbl1(nk)%blocksizes(1) )
          allocate( sendbl1(nk)%displacements(1) )
          allocate( recvbl1(nk)%blocksizes(1) )
          allocate( recvbl1(nk)%displacements(1) )
          allocate( sendbl2(nk)%blocksizes(1) )
          allocate( sendbl2(nk)%displacements(1) )
          allocate( recvbl2(nk)%blocksizes(1) )
          allocate( recvbl2(nk)%displacements(1) )

          sendbl1(nk)%type = MPI_DATATYPE_NULL

          if ( (nk/npr_y) > myid_z .and. mod(nk,npr_y) == myid_y ) then 

             if (mod_geopk .ne. 0) then
                call MPI_TYPE_INDEXED(ione, DTWO*ijtot,   &
                     DTWO*ijtot*(klast-kfirst+1), MPI_DOUBLE_PRECISION, &
                     sendbl1(nk)%type, ierror)
                call MPI_TYPE_COMMIT(sendbl1(nk)%type, ierror)
             endif

             sendbl1(nk)%blocksizes(1) = DTWO*ijtot
             sendbl1(nk)%displacements(1) = DTWO*ijtot*(klast-kfirst+1)
             sendbl1(nk)%partneroffset = myid_z * ijtot * DTWO

          else

             sendbl1(nk)%blocksizes(1) = 0
             sendbl1(nk)%displacements(1) = 0
             sendbl1(nk)%partneroffset = 0

          endif
          sendbl1(nk)%nparcels = size(sendbl1(nk)%displacements)
          sendbl1(nk)%tot_size = sum(sendbl1(nk)%blocksizes)
          max_nparcels = max(max_nparcels, sendbl1(nk)%nparcels)

          recvbl1(nk)%type = MPI_DATATYPE_NULL

          if ( (nk/npr_y) < myid_z .and. mod(nk,npr_y) == myid_y ) then

             if (mod_geopk .ne. 0) then
                call MPI_TYPE_INDEXED(ione, DTWO*ijtot,   &
                     nk/npr_y * ijtot * DTWO, MPI_DOUBLE_PRECISION, &
                     recvbl1(nk)%type, ierror)
                call MPI_TYPE_COMMIT(recvbl1(nk)%type, ierror)
             endif

             recvbl1(nk)%blocksizes(1) = DTWO*ijtot
             recvbl1(nk)%displacements(1) = nk/npr_y * ijtot * DTWO
             recvbl1(nk)%partneroffset = 0

          else

             recvbl1(nk)%blocksizes(1) = 0
             recvbl1(nk)%displacements(1) = 0
             recvbl1(nk)%partneroffset = 0

          endif
          recvbl1(nk)%nparcels = size(recvbl1(nk)%displacements)
          recvbl1(nk)%tot_size = sum(recvbl1(nk)%blocksizes)
          max_nparcels = max(max_nparcels, recvbl1(nk)%nparcels)

          if ( (nk/npr_y) < myid_z .and. mod(nk,npr_y) == myid_y ) then 

             call MPI_TYPE_INDEXED(ione, DTWO*ijtot,   &
                  0, MPI_DOUBLE_PRECISION, &
                  sendbl2(nk)%type, ierror)
             call MPI_TYPE_COMMIT(sendbl2(nk)%type, ierror)

             sendbl2(nk)%blocksizes(1) = DTWO*ijtot
             sendbl2(nk)%displacements(1) = 0
             sendbl2(nk)%partneroffset = (myid_z-nk/npr_y-1) * ijtot * DTWO

          else

             sendbl2(nk)%type = MPI_DATATYPE_NULL

             sendbl2(nk)%blocksizes(1) = 0
             sendbl2(nk)%displacements(1) = 0
             sendbl2(nk)%partneroffset = 0

          endif
          sendbl2(nk)%nparcels = size(sendbl2(nk)%displacements)
          sendbl2(nk)%tot_size = sum(sendbl2(nk)%blocksizes)
          max_nparcels = max(max_nparcels, sendbl2(nk)%nparcels)

          if ( (nk/npr_y) > myid_z .and. mod(nk,npr_y) == myid_y ) then

             call MPI_TYPE_INDEXED(ione, DTWO*ijtot,   &
                  nk/npr_y * ijtot * DTWO, MPI_DOUBLE_PRECISION, &
                  recvbl2(nk)%type, ierror)
             call MPI_TYPE_COMMIT(recvbl2(nk)%type, ierror)

             recvbl2(nk)%blocksizes(1) = DTWO*ijtot
             recvbl2(nk)%displacements(1) = nk/npr_y * ijtot * DTWO
             recvbl2(nk)%partneroffset = 0

          else

             recvbl2(nk)%type = MPI_DATATYPE_NULL

             recvbl2(nk)%blocksizes(1) = 0
             recvbl2(nk)%displacements(1) = 0
             recvbl2(nk)%partneroffset = 0

          endif
          recvbl2(nk)%nparcels = size(recvbl2(nk)%displacements)
          recvbl2(nk)%tot_size = sum(recvbl2(nk)%blocksizes)
          max_nparcels = max(max_nparcels, recvbl2(nk)%nparcels)
        enddo

        call get_partneroffset(grid%commyz, sendbl1, recvbl1)
        call get_partneroffset(grid%commyz, sendbl2, recvbl2)

       endif
      endif

#if (!defined PAREXCH)
      locinbuf = loc(pe16)
#else
      locinbuf = loc(inbuf)
#endif
      locoutbuf = loc(outbuf)
      ptr_inbuf8 = locinbuf
      ptr_outbuf8 = locoutbuf
#endif

! Top down

#if (DSIZE == 16)
!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j, k)
      do k = kfirst,klast
      do j = jfirst,jlast
      do i = ifirst,ilast
         delp16(i,j,k) = delp(i,j,k)
      enddo
      enddo
      enddo
#endif

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j)
      do j = jfirst,jlast
      do i = ifirst,ilast
        pe16(i,j,kfirst) = DP0_0
      enddo
      enddo

! compute partial sums

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j, k)
      do j = jfirst,jlast
        do k = kfirst+1,klast+1
        do i = ifirst,ilast
#if (DSIZE == 16)
          pe16(i,j,k) = pe16(i,j,k-1) + delp16(i,j,k-1)
#else
          pe16(i,j,k) = pe16(i,j,k-1) + delp(i,j,k-1)
#endif
        enddo
        enddo
      enddo

#if defined( SPMD )
      if (npr_z .gt. 1) then

! communicate upward

# if !defined (PAREXCH)
        call mp_sendirr(grid%commyz, sendbl1, recvbl1, inbuf8, outbuf8,            &
                        modc=grid%modc_cdcore )
        call mp_recvirr(grid%commyz, sendbl1, recvbl1, inbuf8, outbuf8,            &
                        modc=grid%modc_cdcore )
# else

        do nk = 0, npr_z-1
          sendcount(nk) = 0
          recvcount(nk) = 0
        enddo

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j, nk)
        do nk = myid_z+1, npr_z-1
          do j = jfirst,jlast
          do i = ifirst,ilast
            inbuf(i,j,nk-myid_z-1) = pe16(i,j,klast+1)
          enddo
          enddo
! Double sendcount since quantities are 16-bytes long
          sendcount(nk) = DTWO*ijtot
        enddo

        call parexchangevector(grid%comm_z, sendcount, inbuf8, recvcount, outbuf8)

# endif

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j, k, nk)
        do k = kfirst,klast+1
          do nk = 0, myid_z-1
          do j = jfirst,jlast
          do i = ifirst,ilast
             pe16(i,j,k) = pe16(i,j,k) + outbuf(i,j,nk)
          enddo
          enddo
          enddo
        enddo

      endif
#endif

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j, k)
      do k = kfirst,klast+1
      do j = jfirst,jlast
      do i = ifirst,ilast
        pe(i,k,j) = pe16(i,j,k) + ptop
        pk(i,j,k) = pe(i,k,j) ** akap
      enddo
      enddo
      enddo

! Bottom up

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j, k)
      do k = kfirst,klast
      do j = jfirst,jlast
      do i = ifirst,ilast
        delp16(i,j,k) = cp*pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))
      enddo
      enddo
      enddo

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j)
      do j = jfirst,jlast
      do i = ifirst,ilast
        pe16(i,j,klast+1) = DP0_0
      enddo
      enddo

! compute partial sums

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j, k)
      do j = jfirst,jlast
        do k = klast,kfirst,-1
        do i = ifirst,ilast
          pe16(i,j,k) = pe16(i,j,k+1) + delp16(i,j,k)
        enddo
        enddo
      enddo

#if defined( SPMD )
      if (npr_z .gt. 1) then

! communicate downward

# if !defined (PAREXCH)
        call mp_sendirr(grid%commyz, sendbl2, recvbl2, inbuf8, outbuf8,            &
                        modc=grid%modc_cdcore )
        call mp_recvirr(grid%commyz, sendbl2, recvbl2, inbuf8, outbuf8,            &
                        modc=grid%modc_cdcore )
# else

        do nk = 0, npr_z-1
          sendcount(nk) = 0
          recvcount(nk) = 0
        enddo

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j, nk)
        do nk = 0, myid_z-1
          do j = jfirst,jlast
          do i = ifirst,ilast
            inbuf(i,j,nk) = pe16(i,j,kfirst)
          enddo
          enddo
! Double sendcount since quantities are 16-bytes long
          sendcount(nk) = DTWO*ijtot
        enddo

        call parexchangevector(grid%comm_z, sendcount, inbuf8, recvcount, outbuf8)

# endif

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j, k, nk)
        do k = kfirst,klast+1
          do nk = myid_z+1, npr_z-1
          do j = jfirst,jlast
          do i = ifirst,ilast
# if !defined (PAREXCH)
            pe16(i,j,k) = pe16(i,j,k) + outbuf(i,j,nk)
# else
            pe16(i,j,k) = pe16(i,j,k) + outbuf(i,j,nk-myid_z-1)
# endif
          enddo
          enddo
          enddo
        enddo

      endif
#endif

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(i, j, k)
      do k = kfirst,klast+1
      do j = jfirst,jlast
      do i = ifirst,ilast
        wz(i,j,k) = pe16(i,j,k) + hs(i,j)
      enddo
      enddo
      enddo

      return
! endif for NO_CRAY_POINTERS
#endif
!EOC
      end subroutine geopk16
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: geopk_d --- Calculate geopotential to the kappa
!
! !INTERFACE:
      subroutine geopk_d( grid, pe, delp, pk, wz, hs, pt, ng, cp, akap )

      use shr_kind_mod,  only : r8 => shr_kind_r8, i8 => shr_kind_i8
      use dynamics_vars, only : T_FVDYCORE_GRID

      implicit none

#if defined ( SPMD )
#include "mpif.h"
#endif

! !INPUT PARAMETERS:

      type (T_FVDYCORE_GRID), intent(in) :: grid
      integer, intent(in)  :: ng      ! Halo size (not always = ng_d)

      real(r8)    akap, cp
      real(r8)    hs(1:grid%im,grid%jfirst:grid%jlast)

! !INPUT PARAMETERS:
      real(r8)    pt(1:grid%im,grid%jfirst-ng:grid%jlast+ng,grid%kfirst:grid%klast) 
      real(r8)  delp(1:grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)

! !OUTPUT PARAMETERS:
      real(r8) wz(1:grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast+1)  ! space N*1 S*1
      real(r8) pk(1:grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast+1)  ! space N*1 S*1
      real(r8) pe(1:grid%im,grid%kfirst:grid%klast+1,grid%jfirst:grid%jlast)  ! temporary variable

! !DESCRIPTION:
!     Calculates geopotential and pressure to the kappa.  This is an expensive
!     operation and several out arrays are kept around for future use.
!     To preserve reproducibility, ordering of transposed-based geopk algorithm
!     is preserved at the cost of a serialization of computation in the Z-direction.
!
! !REVISION HISTORY:
!
!  PW  08.06.27: Original: simple ring r8 version of geopk16 - 
!                  serialized Z-direction but minimized communication overhead
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local:
      real(r8):: ptop
      integer :: km, ifirst, ilast, jfirst, jlast, kfirst, klast
      integer :: npr_z

      integer :: itot, jtot

      integer :: n_blocks
      logical :: sendd   

      integer :: i, il, ilmax, ib, iblksiz
      integer :: j, jl, jlmax, jb, jblksiz
      integer :: k, block, ierror

      integer :: klimits(2), klimits_all(2,0:grid%npr_z-1)
      integer, save :: k_succ_pid, k_pred_pid

      integer, allocatable :: rcvreq(:), sndreq(:)

      real (r8), parameter ::  DP0_0                    =  0.0_r8
      real (r8) l_pe(1:grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast+1)
      real (r8) l_delp(1:grid%im,grid%jfirst:grid%jlast,grid%kfirst:grid%klast)

      real (r8) inbuf(1:grid%im,grid%jfirst:grid%jlast,0:grid%npr_z-1)
      real (r8) outbuf(1:grid%im,grid%jfirst:grid%jlast,0:grid%npr_z-1)
      integer sendcount(0:grid%npr_z-1), recvcount(0:grid%npr_z-1)

      integer first_time_through
      data first_time_through / 0 /

#if defined ( SPMD )
      integer status (MPI_STATUS_SIZE) ! Status of message
#endif

!
! Initialize variables from Grid
!
      ptop     = grid%ptop

      km       = grid%km

      ifirst   = 1               ! 2004.10.04 (WS): Now hardwired for 1..im
      ilast    = grid%im         ! Code was always used in this mode 
      jfirst   = grid%jfirst
      jlast    = grid%jlast
      kfirst   = grid%kfirst
      klast    = grid%klast

      npr_z = grid%npr_z

      itot  = (ilast-ifirst+1)
      jtot  = (jlast-jfirst+1)

      if (grid%modc_cdcore(3) .eq. 1) then
         sendd = .true.
      else
         sendd = .false.
      endif

      n_blocks = max(1,grid%geopkblocks)

      if (n_blocks < jtot) then
         jblksiz = ceiling(float(jtot)/float(n_blocks))
         iblksiz = itot
      else 
         jblksiz = 1
         iblksiz = ceiling(float(itot*jtot)/float(n_blocks))
      endif

      block = 0
      do j=jfirst,jlast,jblksiz
         do i=ifirst,ilast,iblksiz
            block = block + 1
         enddo
      enddo

      allocate( sndreq(block) )
      allocate( rcvreq(block) )

      if (first_time_through .eq. 0) then
         first_time_through = 1
         k_pred_pid = -1
         k_succ_pid = -1
#if defined (SPMD)
         klimits(1) = kfirst
         klimits(2) = klast
         call mpi_allgather (klimits, 2, mpi_integer, &
                             klimits_all, 2, mpi_integer, &
                             grid%comm_z, ierror)
         do i=0,npr_z-1
            if (klimits_all(2,i) == kfirst-1) k_pred_pid = i
            if (klimits_all(1,i) == klast+1)  k_succ_pid = i
         enddo
#endif
      endif

! Top down

! prepost first set of receive requests
#if defined (SPMD)
      if (k_pred_pid /= -1) then
         block = 0
         do j=jfirst,jlast,jblksiz
            if (j+jblksiz > jlast) then
               jb = jlast-j+1
            else
               jb = jblksiz
            endif

            do i=ifirst,ilast,iblksiz
               if (i+iblksiz > ilast) then
                  ib = ilast-i+1
               else
                  ib = iblksiz
               endif

               block = block + 1
               call mpi_irecv (l_pe(i,j,kfirst), jb*ib, &
                               mpi_real8, k_pred_pid, block, &
                               grid%comm_z, rcvreq(block), ierror)
            enddo
         enddo
      endif
#endif

      block = 0
      do j=jfirst,jlast,jblksiz

         if (j+jblksiz > jlast) then
            jb = jlast-j+1
         else
            jb = jblksiz
         endif
         jlmax = j+jb-1

         do i=ifirst,ilast,iblksiz
            if (i+iblksiz > ilast) then
               ib = ilast-i+1
            else
               ib = iblksiz
            endif
            ilmax = i+ib-1

            block = block + 1

! get data from k predecessor
            if (k_pred_pid /= -1) then
#if defined (SPMD)
               call mpi_wait (rcvreq(block), status, ierror)
#endif
            else
               do jl=j,jlmax
                  do il = i,ilmax
                     l_pe(il,jl,kfirst) = DP0_0
                  enddo
               enddo
            endif

! compute partial sums (note that can not thread over k-loop)
            do k = kfirst+1,klast+1
               do jl=j,jlmax
                  do il = i,ilmax
                     l_pe(il,jl,k) = l_pe(il,jl,k-1) + delp(il,jl,k-1)
                  enddo
               enddo
            enddo

! send results to k successor
#if defined (SPMD)
            if (k_succ_pid /= -1) then
               if (sendd) then
                  call mpi_send  (l_pe(i,j,klast+1), jb*ib, mpi_real8, &
                                  k_succ_pid, block, grid%comm_z, &
                                  ierror)
               else
                  call mpi_isend (l_pe(i,j,klast+1), jb*ib, mpi_real8, &
                                  k_succ_pid, block, grid%comm_z, &
                                  sndreq(block), ierror)
               endif
            endif
#endif
!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(il, jl, k)
            do k = kfirst,klast+1
               do jl = j,jlmax
                  do il = i,ilmax
                     pe(il,k,jl) = l_pe(il,jl,k) + ptop
                     pk(il,jl,k) = pe(il,k,jl) ** akap
                  enddo
               enddo
            enddo

#if defined (SPMD)
            if (k_succ_pid /= -1) then
               if (.not. sendd) then
                  call mpi_wait (sndreq(block), status, ierror)
               endif
            endif
#endif
         enddo
      enddo

! Bottom up

! prepost second set of receive requests
#if defined (SPMD)
      if (k_succ_pid /= -1) then
         block = 0
         do j=jfirst,jlast,jblksiz
            if (j+jblksiz > jlast) then
               jb = jlast-j+1
            else
               jb = jblksiz
            endif

            do i=ifirst,ilast,iblksiz
               if (i+iblksiz > ilast) then
                  ib = ilast-i+1
               else
                  ib = iblksiz
               endif

               block = block + 1
               call mpi_irecv (l_pe(i,j,klast+1), jb*ib, &
                               mpi_real8, k_succ_pid, block, &
                               grid%comm_z, rcvreq(block), ierror)
            enddo
         enddo
      endif
#endif

      block = 0
      do j=jfirst,jlast,jblksiz

         if (j+jblksiz > jlast) then
            jb = jlast-j+1
         else
            jb = jblksiz
         endif
         jlmax = j+jb-1

         do i=ifirst,ilast,iblksiz
            if (i+iblksiz > ilast) then
               ib = ilast-i+1
            else
               ib = iblksiz
            endif
            ilmax = i+ib-1

            block = block + 1

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(il, jl, k)
            do k = kfirst,klast
               do jl=j,jlmax
                  do il = i,ilmax
                     l_delp(il,jl,k) = &
                        cp*pt(il,jl,k)*(pk(il,jl,k+1)-pk(il,jl,k))
                  enddo
               enddo
            enddo

! get data from k predecessor
            if (k_succ_pid /= -1) then
#if defined (SPMD)
               call mpi_wait (rcvreq(block), status, ierror)
#endif
            else
               do jl=j,jlmax
                  do il = i,ilmax
                     l_pe(il,jl,klast+1) = DP0_0
                  enddo
               enddo
            endif

! compute partial sums (note that can not thread over k-loop)
            do k = klast,kfirst,-1
               do jl=j,jlmax
                  do il = i,ilmax
                     l_pe(il,jl,k) = l_pe(il,jl,k+1) + l_delp(il,jl,k)
                  enddo
               enddo
            enddo

! send results to k predecessor
#if defined (SPMD)
            if (k_pred_pid /= -1) then
               if (sendd) then
                  call mpi_send  (l_pe(i,j,kfirst), jb*ib, mpi_real8, &
                                  k_pred_pid, block, &
                                  grid%comm_z, ierror)
               else
                  call mpi_isend (l_pe(i,j,kfirst), jb*ib, mpi_real8, &
                                  k_pred_pid, block, &
                                  grid%comm_z, sndreq(block), ierror)
               endif
            endif
#endif

!$omp  parallel do      &
!$omp  default(shared)  &
!$omp  private(il, jl, k)
            do k = kfirst,klast+1
               do jl=j,jlmax
                  do il = i,ilmax
                     wz(il,jl,k) = l_pe(il,jl,k) + hs(il,jl)
                  enddo
               enddo
            enddo

#if defined (SPMD)
            if (k_pred_pid /= -1) then
               if (.not. sendd) then
                  call mpi_wait (sndreq(block), status, ierror)
               endif
            endif
#endif

         enddo

      enddo

      deallocate( sndreq )
      deallocate( rcvreq )

      return
!EOC
      end subroutine geopk_d
!-----------------------------------------------------------------------
