module BandDiagonalMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Band Diagonal matrix solution
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use elm_varctl     , only : iulog
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: BandDiagonal
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine BandDiagonal(numf, lbj, ubj, jtop, jbot, b, r, u)
    !
    ! !DESCRIPTION:
    ! Band diagonal matrix solution for a batch of columns. The arrays are
    ! compressed onto a filter: entry fc holds the system of the fc-th filtered
    ! column, whose unknowns are levels jtop(fc):jbot(fc).
    ! All the systems are solved at once by dgbsv_batch (one thread per system).
    !
    ! !USES:
    use lapack_acc_seq, only : dgbsv_batch
    !
    ! !ARGUMENTS:
    implicit none
    integer , parameter     :: nband = 5                       ! band width
    integer , intent(in)    :: numf                            ! filter dimension
    integer , intent(in)    :: lbj, ubj                        ! lbinning and ubing level indices
    integer , intent(in)    :: jtop(numf)                      ! top level for each column [col]
    integer , intent(in)    :: jbot(numf)                      ! bottom level for each column [col]
    real(r8), intent(in)    :: b(numf, nband, lbj:ubj)         ! compact band matrix [col, nband, j]
    real(r8), intent(in)    :: r(numf, lbj:ubj)                ! "r" rhs of linear system [col, j]
    real(r8), intent(inout) :: u(numf, lbj:ubj)                ! solution [col, j]
    !
    ! ! LOCAL VARIABLES:
    integer , parameter :: kl = (nband-1)/2, ku = kl           ! number of sub/super diagonals
    integer , parameter :: m = 2*kl+ku+1                       ! number of rows of the band storage (dgbsv LDAB)
    integer  :: fc, j, k, n, nmax, ierr                        ! indices
    integer  :: nlev(numf)                                     ! order of each system (jbot-jtop+1)
    integer  :: info(numf)                                     ! dgbsv status of each system
    integer  :: ipiv(numf, ubj-lbj+1)                          ! pivot indices
    real(r8) :: ab(numf, m, ubj-lbj+1)                         ! compact storage array
    real(r8) :: rhs(numf, ubj-lbj+1, 1)                        ! rhs on input, solution on output
    !-----------------------------------------------------------------------

    !Set up input matrix AB
    !An m-by-n band matrix with kl subdiagonals and ku superdiagonals
    !may be stored compactly in a two-dimensional array with
    !kl+ku+1 rows and n columns
    !AB(KL+KU+1+i-j,j) = A(i,j)

    ! the largest system spans every level
    nmax = ubj - lbj + 1

    !$acc data create(nlev, info, ipiv, ab, rhs)

    !$acc parallel loop gang vector collapse(3) default(present)
    do j = 1, nmax
       do k = 1, m
          do fc = 1, numf
             ab(fc,k,j) = 0._r8
          end do
       end do
    end do

    !$acc parallel loop gang vector default(present) private(n)
    do fc = 1, numf
       ! n is the number of levels (snow/soil)
       n = jbot(fc)-jtop(fc)+1
       nlev(fc) = n

       ! on input rhs is the rhs, on output the solution vector
       !$acc loop seq
       do j = 1, n
          rhs(fc,j,1) = r(fc,jtop(fc)+j-1)
       end do

       !$acc loop seq
       do j = 3, n
          ab(fc,kl+ku-1,j) = b(fc,1,jtop(fc)+j-3)   ! 2nd superdiagonal
       end do
       !$acc loop seq
       do j = 2, n
          ab(fc,kl+ku+0,j) = b(fc,2,jtop(fc)+j-2)   ! 1st superdiagonal
       end do
       !$acc loop seq
       do j = 1, n
          ab(fc,kl+ku+1,j) = b(fc,3,jtop(fc)+j-1)   ! diagonal
       end do
       !$acc loop seq
       do j = 1, n-1
          ab(fc,kl+ku+2,j) = b(fc,4,jtop(fc)+j)     ! 1st subdiagonal
       end do
       !$acc loop seq
       do j = 1, n-2
          ab(fc,kl+ku+3,j) = b(fc,5,jtop(fc)+j+1)   ! 2nd subdiagonal
       end do
    end do

    ! DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO ) for every column
    call dgbsv_batch(numf, nmax, kl, ku, 1, nlev, ab, rhs, ipiv, info)

    ierr = 0
    !$acc parallel loop gang vector default(present) reduction(max:ierr)
    do fc = 1, numf
       ierr = max(ierr, abs(info(fc)))
    end do

    if (ierr /= 0) then
       !$acc update self(info, nlev)
       do fc = 1, numf
          if (info(fc) /= 0) then
             write(iulog,*)'filter index: ', fc
             write(iulog,*)'n,kl,ku,m ',nlev(fc),kl,ku,m
             write(iulog,*)'dgbsv info: ',fc,info(fc)
          end if
       end do
       stop
    end if

    !$acc parallel loop gang vector default(present)
    do fc = 1, numf
       !$acc loop seq
       do j = 1, nlev(fc)
          u(fc,jtop(fc)+j-1) = rhs(fc,j,1)
       end do
    end do

    !$acc end data

  end subroutine BandDiagonal

end module BandDiagonalMod
