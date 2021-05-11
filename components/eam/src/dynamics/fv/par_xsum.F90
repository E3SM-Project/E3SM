!-----------------------------------------------------------------------
!BOP
! !ROUTINE: par_xsum --- Calculate x-sum bit-wise consistently
!
! !INTERFACE:
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine par_xsum(grid, a, ltot, sum)
!****6***0*********0*********0*********0*********0*********0**********72
!
! !USES:
#if defined ( SPMD )
      use parutilitiesmodule, only : parexchangevector
#endif
      use dynamics_vars, only : T_FVDYCORE_GRID
      use shr_kind_mod, only: r8 => shr_kind_r8
      use shr_reprosum_mod, only : shr_reprosum_calc, shr_reprosum_tolExceeded, &
                                shr_reprosum_reldiffmax, &
                                shr_reprosum_recompute
      use cam_logfile,   only : iulog
      use FVperf_module, only : FVstartclock, FVstopclock

      implicit none

! !INPUT PARAMETERS:
      type (T_FVDYCORE_GRID), intent(in) :: grid
      integer, intent(in) :: ltot       ! number of quantities to be summed
      ! input vector to be summed
      real (r8), intent(in) :: a(grid%ifirstxy:grid%ilastxy,ltot)

! !OUTPUT PARAMETERS:
      real (r8) sum(ltot)               ! sum of all vector entries

! !DESCRIPTION:
!     This subroutine calculates the sum of "a" in a reproducible
!     (sequentialized) fashion which should give bit-wise identical
!     results irrespective of the number of MPI processes.
!
! !CALLED FROM:
!     te_map
!
! !REVISION HISTORY:
!
!     AAM 00.11.01 : Created
!     WS  03.10.22 : pmgrid removed (now spmd_dyn)
!     WS  04.10.04 : added grid as an argument; removed spmd_dyn
!     WS  05.05.25 : removed ifirst, ilast, im as arguments (in grid)
!     PW  08.06.25 : added fixed point reproducible sum
!
!EOP
!---------------------------------------------------------------------
!BOC
 
! !Local
      real(r8), parameter ::  D0_0                    =  0.0_r8

      real(r8) :: rel_diff(2,ltot)
      real(r8),allocatable :: l_a(:)
      real(r8),allocatable :: a_tmp(:)

      integer :: i,ipe,l,im,lim,nprxy_x,offset
      integer :: sendcount(grid%nprxy_x)
      integer :: recvcount(grid%nprxy_x)

      logical :: write_warning

      im  = grid%im
      lim = (grid%ilastxy-grid%ifirstxy) + 1
      nprxy_x = grid%nprxy_x
      offset  = grid%ifirstxy - 1

      call FVstartclock(grid,'xsum_reprosum')
      call shr_reprosum_calc(a, sum, lim, lim, ltot, gbl_count=im, &
                     commid=grid%commxy_x, rel_diff=rel_diff)
      call FVstopclock(grid,'xsum_reprosum')

      ! check that "fast" reproducible sum is accurate enough. If not, calculate
      ! using old method
      write_warning = .false.
      if (grid%myidxy_x == 0) write_warning = .true.
      if ( shr_reprosum_tolExceeded('par_xsum', ltot, write_warning, &
           iulog, rel_diff) ) then
         if ( shr_reprosum_recompute ) then
            call FVstartclock(grid,'xsum_sumfix')
            allocate( l_a(lim*nprxy_x) )
            allocate( a_tmp(im) )
            sendcount(:) = lim

            do l=1,ltot
               if (rel_diff(1,l) > shr_reprosum_reldiffmax) then
                  sum(l) = D0_0
#if defined ( SPMD )
                  do ipe=1,nprxy_x
                     do i=1,lim
                        l_a(i+(ipe-1)*lim) = a(i+offset,l)
                     enddo
                  enddo
                  call parexchangevector( grid%commxy_x, sendcount, l_a, &
                                          recvcount, a_tmp )
                  do i=1,im
                     sum(l) = sum(l) + a_tmp(i)
                  enddo
#else
                  do i=1,im
                     sum(l) = sum(l) + a(i,l)
                  enddo
#endif
               endif

            enddo

            deallocate( a_tmp )
            deallocate( l_a )
            call FVstopclock(grid,'xsum_sumfix')
         endif
      endif

      return
!EOC
      end subroutine par_xsum
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: par_xsum_r4 --- Calculate x-sum bit-wise consistently (real4)
!
! !INTERFACE:
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine par_xsum_r4(grid, a, ltot, sum)
!****6***0*********0*********0*********0*********0*********0**********72
!
! !USES:
#if defined ( SPMD )
      use parutilitiesmodule, only : parexchangevector
#endif
      use dynamics_vars, only : T_FVDYCORE_GRID
      use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
      use shr_reprosum_mod, only : shr_reprosum_calc, shr_reprosum_tolExceeded, &
                                shr_reprosum_reldiffmax, &
                                shr_reprosum_recompute
      use cam_logfile,   only : iulog
      use FVperf_module, only : FVstartclock, FVstopclock

      implicit none

! !INPUT PARAMETERS:
      type (T_FVDYCORE_GRID), intent(in) :: grid
      integer, intent(in) :: ltot       ! number of quantities to be summed
      real (r4) a(grid%ifirstxy:grid%ilastxy,ltot)    ! input vector to be summed

! !OUTPUT PARAMETERS:
      real (r8) sum(ltot)               ! sum of all vector entries

! !DESCRIPTION:
!     This subroutine calculates the sum of "a" in a reproducible
!     (sequentialized) fashion which should give bit-wise identical
!     results irrespective of the number of MPI processes.
!
! !REVISION HISTORY:
!
!     WS  05.04.08 : Created from par_xsum
!     WS  05.05.25 : removed ifirst, ilast, im as arguments (in grid)
!     WS  06.06.28 : Fixed bug in sequential version
!     PW  08.06.25 : added fixed point reproducible sum
!
!EOP
!---------------------------------------------------------------------
!BOC
 
! !Local
      real(r8), parameter ::  D0_0                    =  0.0_r8

      real(r8) :: a8(grid%ifirstxy:grid%ilastxy,ltot)
      real(r8) :: rel_diff(2,ltot)
      real(r4),allocatable :: l_a(:)
      real(r4),allocatable :: a_tmp(:)

      integer i,ipe,l,im,lim,nprxy_x,offset
      integer sendcount(grid%nprxy_x)
      integer recvcount(grid%nprxy_x)

      logical :: write_warning

      im  = grid%im
      lim = (grid%ilastxy-grid%ifirstxy) + 1
      nprxy_x = grid%nprxy_x
      offset  = grid%ifirstxy - 1

      call FVstartclock(grid,'xsum_r4_reprosum')
      a8(:,:) = a(:,:)
      call shr_reprosum_calc(a8, sum, lim, lim, ltot, gbl_count=im, &
                     commid=grid%commxy_x, rel_diff=rel_diff)
      call FVstopclock(grid,'xsum_r4_reprosum')

      ! check that "fast" reproducible sum is accurate enough. If not, calculate
      ! using old method
      write_warning = .false.
      if (grid%myidxy_x == 0) write_warning = .true.
      if ( shr_reprosum_tolExceeded('par_xsum_r4', ltot, write_warning, &
           iulog, rel_diff) ) then
         if ( shr_reprosum_recompute ) then
            call FVstartclock(grid,'xsum_r4_sumfix')
            allocate( l_a(lim*nprxy_x) )
            allocate( a_tmp(im) )
            sendcount(:) = lim

            do l=1,ltot
               if (rel_diff(1,l) > shr_reprosum_reldiffmax) then
                  sum(l) = D0_0
#if defined ( SPMD )
                  do ipe=1,nprxy_x
                     do i=1,lim
                        l_a(i+(ipe-1)*lim) = a(i+offset,l)
                     enddo
                  enddo
                  call parexchangevector( grid%commxy_x, sendcount, l_a, &
                                          recvcount, a_tmp )
                  do i=1,im
                     sum(l) = sum(l) + a_tmp(i)
                  enddo
#else
                  do i=1,im
                     sum(l) = sum(l) + a(i,l)
                  enddo
#endif
               endif

            enddo

            deallocate( a_tmp )
            deallocate( l_a )
            call FVstopclock(grid,'xsum_r4_sumfix')
         endif
      endif

      return
!EOC
      end subroutine par_xsum_r4
!-----------------------------------------------------------------------
