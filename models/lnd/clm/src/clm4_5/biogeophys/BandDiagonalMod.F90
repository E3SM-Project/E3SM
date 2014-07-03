module BandDiagonalMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Band Diagonal matrix solution
  !
  !
  ! !USES:
  use decompMod   , only : bounds_type
  use abortutils  , only : endrun
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl  , only : iulog
  use shr_log_mod , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: BandDiagonal
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine BandDiagonal(bounds, lbj, ubj, jtop, jbot, numf, filter, nband, b, r, u)
    !
    ! !DESCRIPTION:
    ! Tridiagonal matrix solution
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                    ! bounds
    integer , intent(in)    :: lbj, ubj                        ! lbinning and ubing level indices
    integer , intent(in)    :: jtop( bounds%begc: )            ! top level for each column [col]
    integer , intent(in)    :: jbot( bounds%begc: )            ! bottom level for each column [col]
    integer , intent(in)    :: numf                            ! filter dimension
    integer , intent(in)    :: nband                           ! band width
    integer , intent(in)    :: filter(:)                       ! filter
    real(r8), intent(in)    :: b( bounds%begc: , 1:   , lbj: ) ! compact band matrix [col, nband, j]
    real(r8), intent(in)    :: r( bounds%begc: , lbj: )        ! "r" rhs of linear system [col, j]
    real(r8), intent(inout) :: u( bounds%begc: , lbj: )        ! solution [col, j]
    !
    ! ! LOCAL VARIABLES:
    integer  :: j,ci,fc,info,m,n              !indices
    integer  :: kl,ku                         !number of sub/super diagonals
    integer, allocatable :: ipiv(:)           !temporary
    real(r8),allocatable :: ab(:,:),temp(:,:) !compact storage array
    real(r8),allocatable :: result(:)

    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(jtop) == (/bounds%endc/)),             errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(jbot) == (/bounds%endc/)),             errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(b)    == (/bounds%endc, nband, ubj/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(r)    == (/bounds%endc, ubj/)),        errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(u)    == (/bounds%endc, ubj/)),        errMsg(__FILE__, __LINE__))


!!$     SUBROUTINE SGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
!!$*
!!$*  -- LAPACK driver routine (version 3.1) --
!!$*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!!$*     November 2006
!!$*
!!$*     .. Scalar Arguments ..
!!$      INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
!!$*     ..
!!$*     .. Array Arguments ..
!!$      INTEGER            IPIV( * )
!!$      REAL               AB( LDAB, * ), B( LDB, * )
!!$*     ..
!!$*
!!$*  Purpose
!!$*  =======
!!$*
!!$*  SGBSV computes the solution to a real system of linear equations
!!$*  A * X = B, where A is a band matrix of order N with KL subdiagonals
!!$*  and KU superdiagonals, and X and B are N-by-NRHS matrices.
!!$*
!!$*  The LU decomposition with partial pivoting and row interchanges is
!!$*  used to factor A as A = L * U, where L is a product of permutation
!!$*  and unit lower triangular matrices with KL subdiagonals, and U is
!!$*  upper triangular with KL+KU superdiagonals.  The factored form of A
!!$*  is then used to solve the system of equations A * X = B.
!!$*
!!$*  Arguments
!!$*  =========
!!$*
!!$*  N       (input) INTEGER
!!$*          The number of linear equations, i.e., the order of the
!!$*          matrix A.  N >= 0.
!!$*
!!$*  KL      (input) INTEGER
!!$*          The number of subdiagonals within the band of A.  KL >= 0.
!!$*
!!$*  KU      (input) INTEGER
!!$*          The number of superdiagonals within the band of A.  KU >= 0.
!!$*
!!$*  NRHS    (input) INTEGER
!!$*          The number of right hand sides, i.e., the number of columns
!!$*          of the matrix B.  NRHS >= 0.
!!$*
!!$*  AB      (input/output) REAL array, dimension (LDAB,N)
!!$*          On entry, the matrix A in band storage, in rows KL+1 to
!!$*          2*KL+KU+1; rows 1 to KL of the array need not be set.
!!$*          The j-th column of A is stored in the j-th column of the
!!$*          array AB as follows:
!!$*          AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)
!!$*          On exit, details of the factorization: U is stored as an
!!$*          upper triangular band matrix with KL+KU superdiagonals in
!!$*          rows 1 to KL+KU+1, and the multipliers used during the
!!$*          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!!$*          See below for further details.
!!$*
!!$*  LDAB    (input) INTEGER
!!$*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!!$*
!!$*  IPIV    (output) INTEGER array, dimension (N)
!!$*          The pivot indices that define the permutation matrix P;
!!$*          row i of the matrix was interchanged with row IPIV(i).
!!$*
!!$*  B       (input/output) REAL array, dimension (LDB,NRHS)
!!$*          On entry, the N-by-NRHS right hand side matrix B.
!!$*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!!$*
!!$*  LDB     (input) INTEGER
!!$*          The leading dimension of the array B.  LDB >= max(1,N).
!!$*
!!$*  INFO    (output) INTEGER
!!$*          = 0:  successful exit
!!$*          < 0:  if INFO = -i, the i-th argument had an illegal value
!!$*          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
!!$*                has been completed, but the factor U is exactly
!!$*                singular, and the solution has not been computed.
!!$*
!!$*  Further Details
!!$*  ===============
!!$*
!!$*  The band storage scheme is illustrated by the following example, when
!!$*  M = N = 6, KL = 2, KU = 1:
!!$*
!!$*  On entry:                       On exit:
!!$*
!!$*      *    *    *    +    +    +       *    *    *   u14  u25  u36
!!$*      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!!$*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!!$*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!!$*     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!!$*     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!!$*
!!$*  Array elements marked * are not used by the routine; elements marked
!!$*  + need not be set on entry, but are required by the routine to store
!!$*  elements of U because of fill-in resulting from the row interchanges.


!Set up input matrix AB
!An m-by-n band matrix with kl subdiagonals and ku superdiagonals 
!may be stored compactly in a two-dimensional array with 
!kl+ku+1 rows and n columns
!AB(KL+KU+1+i-j,j) = A(i,j)

    do fc = 1,numf
       ci = filter(fc)

       kl=(nband-1)/2
       ku=kl
! m is the number of rows required for storage space by dgbsv
       m=2*kl+ku+1
! n is the number of levels (snow/soil)
!scs: replace ubj with jbot
       n=jbot(ci)-jtop(ci)+1

       allocate(ab(m,n))
       ab=0.0

       ab(kl+ku-1,3:n)=b(ci,1,jtop(ci):jbot(ci)-2)   ! 2nd superdiagonal
       ab(kl+ku+0,2:n)=b(ci,2,jtop(ci):jbot(ci)-1)   ! 1st superdiagonal
       ab(kl+ku+1,1:n)=b(ci,3,jtop(ci):jbot(ci))     ! diagonal
       ab(kl+ku+2,1:n-1)=b(ci,4,jtop(ci)+1:jbot(ci)) ! 1st subdiagonal
       ab(kl+ku+3,1:n-2)=b(ci,5,jtop(ci)+2:jbot(ci)) ! 2nd subdiagonal

       allocate(temp(m,n))
       temp=ab

       allocate(ipiv(n))
       allocate(result(n))

! on input result is rhs, on output result is solution vector
       result(:)=r(ci,jtop(ci):jbot(ci))

!       DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
       call dgbsv( n, kl, ku, 1, ab, m, ipiv, result, n, info )
       u(ci,jtop(ci):jbot(ci))=result(:)

       if(info /= 0) then 
          write(iulog,*)'index: ', ci
          write(iulog,*)'n,kl,ku,m ',n,kl,ku,m
          write(iulog,*)'dgbsv info: ',ci,info
          
          write(iulog,*) ''
          write(iulog,*) 'ab matrix'
          do j=1,n
             !             write(iulog,'(i2,7f18.7)') j,temp(:,j)
             write(iulog,'(i2,5f18.7)') j,temp(3:7,j)
          enddo
          write(iulog,*) ''
          stop
       endif
       deallocate(temp)

       deallocate(ab)
       deallocate(ipiv)
       deallocate(result)
    end do

  end subroutine BandDiagonal

end module BandDiagonalMod
