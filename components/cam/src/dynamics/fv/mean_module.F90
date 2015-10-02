module mean_module

  use shr_kind_mod,  only : r8 => shr_kind_r8
  use shr_reprosum_mod, only : shr_reprosum_calc, shr_reprosum_tolExceeded, &
                            shr_reprosum_recompute
  use perf_mod
  use cam_logfile,   only : iulog

  public gmean, gmeanxy 

  private
  real(r8), parameter ::  D0_0                    =  0.0_r8

contains
!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  gmean --- Calculate the mean of a 2D field
!
! !INTERFACE:

subroutine gmean(grid, q, qmean)

! !USES:
  use commap, only : w
  use dynamics_vars, only : T_FVDYCORE_GRID

#if defined( SPMD )
  use parutilitiesmodule, only : parcollective, sumop
#endif

  implicit none

! !INPUT PARAMETERS:

  type (T_FVDYCORE_GRID), intent(in) :: grid                ! Grid information
  real(r8), intent(in) :: q(grid%im,grid%jfirst:grid%jlast) ! 2D field 

  real(r8) qmean

! !DESCRIPTION:
!     Calculate the mean of a 2D field
!
! !REVISION HISTORY:
!   00.08.01   Lin     Creation
!   01.01.10   Lin     Revised
!   01.06.27   Mirin   Use y communicator
!   05.07.12   Sawyer  Simplified interface with grid argument
!
!EOP
!-----------------------------------------------------------------------
!BOC

  real(r8) :: xsum(grid%jm)
  integer  :: i, j, im, jm, jfirst, jlast

  im      = grid%im
  jm      = grid%jm
  jfirst  = grid%jfirst
  jlast   = grid%jlast

  do j=1,jm
     xsum(j) = D0_0
  enddo
  do j=jfirst,jlast
     do i=1,im
        xsum(j) = xsum(j) + q(i,j)
     enddo
     xsum(j) = xsum(j)*w(j)
  enddo

#if defined( SPMD )
  if (grid%npr_y .ne. 1) then
     call parcollective( grid%comm_y, sumop, jm, xsum )
  endif
#endif

  qmean = D0_0
  do j=1,jm
     qmean = qmean + xsum(j)
  enddo
  qmean = qmean / (2*im)

  return
!EOC
end subroutine gmean
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  gmeanxy --- Calculate the mean of a 2D field (XY decomp)
!
! !INTERFACE:

subroutine gmeanxy(grid, q, qmean)

! !USES:
  use commap, only : w
  use dynamics_vars, only : T_FVDYCORE_GRID

#if defined( SPMD )
  use parutilitiesmodule, only : parcollective, sumop
#endif

  implicit none

! !INPUT PARAMETERS:

  type (T_FVDYCORE_GRID), intent(in) :: grid            ! Grid information
  real(r8), intent(in), target :: q(grid%ifirstxy:grid%ilastxy,            &
                                    grid%jfirstxy:grid%jlastxy) ! 2D field 

  real(r8) qmean

! !DESCRIPTION:
!     Calculate the mean of a 2D field on an XY decomposition
!     This is inefficiently programmed (global collective operation),
!     and is therefore only intended for initialization phase.
!
!     PW: gmeanxy is called in fv_prints, so replaced inefficient algorithm
!     with shr_reprosum_calc.
!
! !REVISION HISTORY:
!   00.08.01   Lin     Creation
!   01.01.10   Lin     Revised
!   01.06.27   Mirin   Use y communicator
!   05.07.12   Sawyer  Simplified interface with grid argument
!   05.08.26   Sawyer  Modified for XY decomposition
!   08.07.03   Worley  Introduced repro_sum logic
!   12.10.29   Santos  repro_sum is now shr_reprosum
!
!EOP
!-----------------------------------------------------------------------
!BOC

  real(r8) :: q_tmp(grid%ifirstxy:grid%ilastxy, &
                    grid%jfirstxy:grid%jlastxy)
  real(r8) :: rel_diff(2), qmean_tmp(1), xsum
  real(r8), allocatable :: q_global(:,:)

  integer  :: i, j, im, jm, ifirstxy, ilastxy, jfirstxy, jlastxy
  integer  :: lim, ljm, lijm

  logical  :: write_warning

  im        = grid%im
  jm        = grid%jm
  ifirstxy  = grid%ifirstxy
  ilastxy   = grid%ilastxy
  jfirstxy  = grid%jfirstxy
  jlastxy   = grid%jlastxy

  lim = ilastxy - ifirstxy + 1
  ljm = jlastxy - jfirstxy + 1
  lijm = lim*ljm

  do j=jfirstxy,jlastxy
     do i=ifirstxy,ilastxy
        q_tmp(i,j) = q(i,j)*w(j)
     enddo
  enddo

  call t_startf("gmeanxy_reprosum")
  call shr_reprosum_calc(q_tmp, qmean_tmp, lijm, lijm, 1, gbl_count=im*jm, &
                 commid=grid%commxy, rel_diff=rel_diff)
  qmean = qmean_tmp(1)
  call t_stopf("gmeanxy_reprosum")

  ! check that "fast" reproducible sum is accurate enough. If not, calculate
  ! using old method
  write_warning = .false.
  if (grid%iam == 0) write_warning = .true.
  if ( shr_reprosum_tolExceeded('gmeanxy', 1, write_warning, &
                              iulog, rel_diff) ) then
     if ( shr_reprosum_recompute ) then
        call t_startf("gmeanxy_sumfix")
        allocate( q_global(im,jm) )
        q_global = D0_0
        do j=jfirstxy,jlastxy
           do i=ifirstxy,ilastxy
              q_global(i,j) = q_tmp(i,j)
           enddo
        enddo

#if defined( SPMD )
        call parcollective( grid%commxy, sumop, im, jm, q_global )
#endif
        qmean = D0_0
        do j=1,jm
           xsum = D0_0
           do i=1,im
              xsum = xsum + q_global(i,j)
           enddo
           qmean = qmean + xsum
        enddo

        deallocate( q_global )
        call t_stopf("gmeanxy_sumfix")
     endif
  endif

  qmean = qmean / (2*im)

  return
!EOC
end subroutine gmeanxy
!-----------------------------------------------------------------------

end module mean_module
