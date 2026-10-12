module lapack_acc_seq
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  !! This module contains the needed LAPACK/BLAS routines that need
  !! to be called on the GPU.  As of CUDA 10.1, cuBLAS being called in
  !! device code is no longer supported

contains

INTEGER FUNCTION idamax(N,DX,INCX)
  !$acc routine seq
  !*
  !*  -- Reference BLAS level1 routine (version 3.8.0) --
  !*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
  !*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  !*     November 2017
  !*
  !*     .. Scalar Arguments ..
  integer :: incx,n
  !*     ..
  !*     .. Array Arguments ..
  real(r8) :: dx(*)
  !*     ..
  !*
  !*  =====================================================================
  !*
  !*     .. Local Scalars ..
  real(r8) :: dmax
  integer :: i,ix
  !*     ..
  !*     .. Intrinsic Functions ..
  !*     ..
  idamax = 0
  IF (n.LT.1 .OR. incx.LE.0) RETURN
  idamax = 1
  IF (n.EQ.1) RETURN
  IF (incx.EQ.1) THEN
    !
    !        code for increment equal to 1
    !
     dmax = dabs(dx(1))
     DO i = 2,n
        IF (dabs(dx(i)).GT.dmax) THEN
           idamax = i
           dmax = dabs(dx(i))
        END IF
     END DO
  ELSE
    !
    !        code for increment not equal to 1
    !
     ix = 1
     dmax = dabs(dx(1))
     ix = ix + incx
     DO i = 2,n
        IF (dabs(dx(ix)).GT.dmax) THEN
           idamax = i
           dmax = dabs(dx(ix))
        END IF
        ix = ix + incx
     END DO
  END IF
  RETURN
END FUNCTION idamax

 subroutine dswap_oacc(N,DX,INCX,DY,INCY)
  !$acc routine seq
  !*
  !*  -- Reference BLAS level1 routine (version 3.8.0) --
  !*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
  !*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  !*     November 2017
  !*
  !*     .. Scalar Arguments ..
  integer :: INCX,INCY,N
  !*     ..
  !*     .. Array Arguments ..
  real(r8) :: DX(*),DY(*)
  !*     ..
  !*
  !*  =====================================================================
  !*
  !*     .. Local Scalars ..
  real(r8) :: DTEMP
  integer :: I,IX,IY,M,MP1
  !*     ..
  !*     .. Intrinsic Functions ..
  !*     ..
  IF (n.LE.0) RETURN
  IF (incx.EQ.1 .AND. incy.EQ.1) THEN
    !*
    !*       code for both increments equal to 1
    !*
    !*
    !*       clean-up loop
    !*
     m = mod(n,3)
     IF (m.NE.0) THEN
        DO i = 1,m
           dtemp = dx(i)
           dx(i) = dy(i)
           dy(i) = dtemp
        END DO
        IF (n.LT.3) RETURN
     END IF
     mp1 = m + 1
     DO i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i+1)
        dx(i+1) = dy(i+1)
        dy(i+1) = dtemp
        dtemp = dx(i+2)
        dx(i+2) = dy(i+2)
        dy(i+2) = dtemp
     END DO
  ELSE
    !*
    !*       code for unequal increments or equal increments not equal
    !*         to 1
    !*
     ix = 1
     iy = 1
     IF (incx.LT.0) ix = (-n+1)*incx + 1
     IF (incy.LT.0) iy = (-n+1)*incy + 1
     DO i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
     END DO
  END IF
  RETURN
 end subroutine dswap_oacc

 subroutine dscal_oacc(N,DA,DX,INCX)
  !$acc routine seq
  !*
  !*  -- Reference BLAS level1 routine (version 3.8.0) --
  !*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
  !*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  !*     November 2017
  !*
  !*     .. Scalar Arguments ..
 real(r8) :: DA
 integer :: INCX,N
 !*     ..
 !*     .. Array Arguments ..
 real(r8) :: DX(*)
 !*     ..
 !*
 !*  =====================================================================
 !*
 !*     .. Local Scalars ..
 integer :: I,M,MP1,NINCX
 !*     ..
 !*     .. Intrinsic Functions ..
 !*     ..
 IF (n.LE.0 .OR. incx.LE.0) RETURN
 IF (incx.EQ.1) THEN
   !*
   !*        code for increment equal to 1
   !*
   !*
   !*        clean-up loop
   !*
    m = mod(n,5)
    IF (m.NE.0) THEN
       DO i = 1,m
          dx(i) = da*dx(i)
       END DO
       IF (n.LT.5) RETURN
    END IF
    mp1 = m + 1
    DO i = mp1,n,5
       dx(i) = da*dx(i)
       dx(i+1) = da*dx(i+1)
       dx(i+2) = da*dx(i+2)
       dx(i+3) = da*dx(i+3)
       dx(i+4) = da*dx(i+4)
    END DO
 ELSE
   !*
   !*        code for increment not equal to 1
   !*
    nincx = n*incx
    DO i = 1,nincx,incx
       dx(i) = da*dx(i)
    END DO
  END IF
  RETURN
 end subroutine dscal_oacc

 subroutine dger_oacc(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
  !$acc routine seq
  !*
  !*  -- Reference BLAS level2 routine (version 3.7.0) --
  !*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
  !*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  !*     December 2016
  !*
  !*     .. Scalar Arguments ..
 real(r8) :: ALPHA
 integer :: INCX,INCY,LDA,M,N
 !*     ..
 !*     .. Array Arguments ..
 real(r8) :: A(LDA,*),X(*),Y(*)
 !*     ..
 !*
 !*  =====================================================================
 !*
 !*     .. Parameters ..
 real(r8), parameter :: zero = 0.0_r8
 !*     ..
 !*     .. Local Scalars ..
 real(r8) :: TEMP
 integer :: I,INFO,IX,J,JY,KX
 !*     ..
 !*     .. Intrinsic Functions ..
 !*     ..
 !*
 !*     Test the input parameters.
 !*
 info = 0
 IF (m.LT.0) THEN
     info = 1
 ELSE IF (n.LT.0) THEN
     info = 2
 ELSE IF (incx.EQ.0) THEN
     info = 5
 ELSE IF (incy.EQ.0) THEN
     info = 7
 ELSE IF (lda.LT.max(1,m)) THEN
     info = 9
 END IF
 IF (info.NE.0) THEN
     print *,"ERROR",info
     RETURN
 END IF
 !*
 !*     Quick return if possible.
 !*
 IF ((m.EQ.0) .OR. (n.EQ.0) .OR. (alpha.EQ.zero)) RETURN
 !*
 !*     Start the operations. In this version the elements of A are
 !*     accessed sequentially with one pass through A.
 !*
 IF (incy.GT.0) THEN
     jy = 1
 ELSE
     jy = 1 - (n-1)*incy
 END IF
 IF (incx.EQ.1) THEN
     DO j = 1,n
         IF (y(jy).NE.zero) THEN
             temp = alpha*y(jy)
             DO i = 1,m
                 a(i,j) = a(i,j) + x(i)*temp
             END DO
         END IF
         jy = jy + incy
     END DO
 ELSE
     IF (incx.GT.0) THEN
         kx = 1
     ELSE
         kx = 1 - (m-1)*incx
     END IF
     DO j = 1,n
         IF (y(jy).NE.zero) THEN
             temp = alpha*y(jy)
             ix = kx
             DO i = 1,m
                 a(i,j) = a(i,j) + x(ix)*temp
                 ix = ix + incx
             END DO
         END IF
         jy = jy + incy
     END DO
  END IF
  !*
  RETURN
  !*
  !*     End of DGER  .
  !*
  end subroutine dger_oacc

SUBROUTINE dgbtf2_oacc( M, N, KL, KU, AB, LDAB, IPIV, INFO )
  !$acc routine seq
  !  -- LAPACK computational routine (version 3.7.0) --
  !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  !     December 2016
  !
  !     .. Scalar Arguments ..
 integer :: INFO, KL, KU, LDAB, M, N
 !     ..
 !     .. Array Arguments ..
 integer :: IPIV( * )
 real(r8) :: AB( LDAB, * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
 real(r8), parameter :: one = 1.0_r8, zero = 0.0_r8
 !     ..
 !     .. Local Scalars ..
 integer :: I, J, JP, JU, KM, KV
 !     ..
 !     .. External Functions ..
 !INTEGER            IDAMAX
 !EXTERNAL           idamax
 !     ..
 !     .. External Subroutines ..
 !EXTERNAL           dger, dscal, dswap
 !     ..
 !     .. Intrinsic Functions ..
 !     ..
 !     .. Executable Statements ..
 !
 !     KV is the number of superdiagonals in the factor U, allowing for
 !     fill-in.
 !
 kv = ku + kl
 !
 !     Test the input parameters.
 !
 info = 0
 IF( m.LT.0 ) THEN
    info = -1
 ELSE IF( n.LT.0 ) THEN
    info = -2
 ELSE IF( kl.LT.0 ) THEN
    info = -3
 ELSE IF( ku.LT.0 ) THEN
    info = -4
 ELSE IF( ldab.LT.kl+kv+1 ) THEN
    info = -6
 END IF
 IF( info.NE.0 ) THEN
    print *, info, "error in "
    RETURN
 END IF
 !
 !     Quick return if possible
 !
 IF( m.EQ.0 .OR. n.EQ.0 ) RETURN
 !
 !     Gaussian elimination with partial pivoting
 !
 !     Set fill-in elements in columns KU+2 to KV to zero.
 !
 DO j = ku + 2, min( kv, n )
    DO i = kv - j + 2, kl
       ab( i, j ) = zero
    END DO
 END DO
  !
  ! JU is the index of the last column affected by the current stage
  ! of the factorization.
  !
 ju = 1
 !
 DO j = 1, min( m, n )
   !
   !  Set fill-in elements in column J+KV to zero.
   !
    IF( j+kv.LE.n ) THEN
       DO i = 1, kl
          ab( i, j+kv ) = zero
       END DO
    END IF
    !
    ! Find pivot and test for singularity. KM is the number of
    ! subdiagonal elements in the current column.
    !
    km = min( kl, m-j )
    jp = idamax( km+1, ab( kv+1, j ), 1 )
    ipiv( j ) = jp + j - 1
    IF( ab( kv+jp, j ).NE.zero ) THEN
       ju = max( ju, min( j+ku+jp-1, n ) )
       !
       ! Apply interchange to columns J to JU.
       !
       IF( jp.NE.1 ) &
         CALL dswap_oacc( ju-j+1, ab( kv+jp, j ), ldab-1, &
                     ab( kv+1, j ), ldab-1 )

       IF( km.GT.0 ) THEN
         !
         !   Compute multipliers.
         !
          CALL dscal_oacc( km, one / ab( kv+1, j ), ab( kv+2, j ), 1 )
          !
          !  Update trailing submatrix within the band.
          !
          IF( ju.GT.j ) &
            CALL dger_oacc( km, ju-j, -one, ab( kv+2, j ), 1, &
                       ab( kv, j+1 ), ldab-1, ab( kv+1, j+1 ), &
                       ldab-1 )
       END IF
    ELSE
      !
      !  If pivot is zero, set INFO to the index of the pivot
      !  unless a zero pivot has already been found.
      !
       IF( info.EQ.0 )  info = j

    END IF
 END DO
 RETURN
!
!     End of DGBTF2
!
END SUBROUTINE DGBTF2_oacc

SUBROUTINE dgbtrf_oacc( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!$acc routine seq
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
integer :: INFO, KL, KU, LDAB, M, N
!     ..
!     .. Array Arguments ..
integer :: IPIV( * )
real(r8) :: AB( LDAB, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
real(r8), parameter :: one = 1.0_r8, zero = 0.0_r8
integer, parameter :: nbmax = 64, ldwork = nbmax+1
!     ..
!     .. Local Scalars ..
integer :: I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP, &
                   JU, K2, KM, KV, NB, NW
real(r8) :: TEMP
!     ..
!     .. Local Arrays ..
real(r8) :: WORK13( LDWORK, NBMAX ), &
                   WORK31( LDWORK, NBMAX )
!     ..
!     .. External Functions ..
!INTEGER            IDAMAX !, ILAENV
!EXTERNAL           idamax !, ilaenv
!     ..
!     .. External Subroutines ..
!EXTERNAL           dgbtf2, dger, dscal, dswap
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in
!
kv = ku + kl
!
!     Test the input parameters.
!
info = 0
IF( m.LT.0 ) THEN
   info = -1
ELSE IF( n.LT.0 ) THEN
   info = -2
ELSE IF( kl.LT.0 ) THEN
   info = -3
ELSE IF( ku.LT.0 ) THEN
   info = -4
ELSE IF( ldab.LT.kl+kv+1 ) THEN
   info = -6
END IF
IF( info.NE.0 ) THEN
   print *, "error in dgbtrf", info
   RETURN
END IF
!
!     Quick return if possible
!
IF( m.EQ.0 .OR. n.EQ.0 )   RETURN
!
!     Determine the block size for this environment
!
!nb = ilaenv( 1, 'DGBTRF', ' ', m, n, kl, ku )
!
!     The block size must not exceed the limit set by the size of the
!     local arrays WORK13 and WORK31.
!
!
!        Use unblocked code
!
   CALL dgbtf2_oacc( m, n, kl, ku, ab, ldab, ipiv, info )

!
RETURN
!
!     End of DGBTRF
!
END subroutine dgbtrf_oacc


SUBROUTINE dtbsv_oacc(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
  !$acc routine seq
  !NOTE : FOR DGBTRS WE HARDCODE UPLO='UPPER',TRANS='NO TRANSPOSE',DIAG='NON-UNIT'
  !*
  !*  -- Reference BLAS level2 routine (version 3.7.0) --
  !*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
  !*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  !*     December 2016
  !*
  !*     .. Scalar Arguments ..
 integer :: INCX,K,LDA,N
 integer :: DIAG,TRANS,UPLO
 !*     ..
 !*     .. Array Arguments ..
 real(r8) :: A(LDA,*),X(*)
 !*     ..
 !*
 !*  =====================================================================
 !*
 !*     .. Parameters ..
 real(r8), parameter :: zero = 0.0_r8
 !*     ..
 !*     .. Local Scalars ..
 real(r8) :: TEMP
 integer :: I,INFO,IX,J,JX,KPLUS1,KX,L
 logical :: NOUNIT
 !*     ..
 !*     .. External Subroutines ..
 !*     ..
 !*     .. Intrinsic Functions ..
 !*     ..
 !*
 !*     Test the input parameters.
 !*
 info = 0
 IF (uplo .ne. 1 ) THEN
     info = 1
 ELSE IF (trans .ne. 0 ) THEN
     info = 2
 ELSE IF (diag .ne. 0) THEN
     info = 3
 ELSE IF (n.LT.0) THEN
     info = 4
 ELSE IF (k.LT.0) THEN
     info = 5
 ELSE IF (lda.LT. (k+1)) THEN
     info = 7
 ELSE IF (incx.EQ.0) THEN
     info = 9
 END IF
 IF (info.NE.0) THEN
     print *, "ERRROR in dtbsv", info
     RETURN
 END IF
 !*
 !*     Quick return if possible.
 !*
 IF (n.EQ.0) RETURN
 !*
 nounit = .true.
 !*
 !*     Set up the start point in X if the increment is not unity. This
 !*     will be  ( N - 1 )*INCX  too small for descending loops.
 !*
 IF (incx.LE.0) THEN
     kx = 1 - (n-1)*incx
 ELSE IF (incx.NE.1) THEN
     kx = 1
 END IF
 !*
 !*     Start the operations. In this version the elements of A are
 !*     accessed by sequentially with one pass through A.
 !*
   !*
   !*        Form  x := inv( A )*x.
   !*
         kplus1 = k + 1
         IF (incx.EQ.1) THEN
             DO j = n,1,-1
                 IF (x(j).NE.zero) THEN
                     l = kplus1 - j
                     IF (nounit) x(j) = x(j)/a(kplus1,j)
                     temp = x(j)
                     DO i = j - 1,max(1,j-k),-1
                         x(i) = x(i) - temp*a(l+i,j)
                     END DO
                 END IF
             END DO
         ELSE
             kx = kx + (n-1)*incx
             jx = kx
             DO j = n,1,-1
                 kx = kx - incx
                 IF (x(jx).NE.zero) THEN
                     ix = kx
                     l = kplus1 - j
                     IF (nounit) x(jx) = x(jx)/a(kplus1,j)
                     temp = x(jx)
                     DO i = j - 1,max(1,j-k),-1
                         x(ix) = x(ix) - temp*a(l+i,j)
                         ix = ix - incx
                     END DO
                 END IF
                 jx = jx - incx
             END DO
         END IF
         !*
 RETURN
 !*
 !*     End of DTBSV .
 !*
END subroutine DTBSV_oacc

SUBROUTINE dgbtrs_oacc( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
  !$acc routine seq
  !  -- LAPACK computational routine (version 3.7.0) --
  !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  !     December 2016
  !
  !     .. Scalar Arguments ..
   integer :: TRANS
   integer :: INFO, KL, KU, LDAB, LDB, N, NRHS
   !     ..
   !     .. Array Arguments ..
   integer :: IPIV( * )
   real(r8) :: AB( LDAB, * ), B( LDB, * )
   !     ..
   !
   !  =====================================================================
   !
   !     .. Parameters ..
   real(r8), parameter :: one = 1.0_r8
   !     ..
   !     .. Local Scalars ..
   logical :: LNOTI, NOTRAN
   integer :: I, J, KD, L, LM
   !*     ..
   !*     .. External Functions ..
   !LOGICAL            LSAME
   !EXTERNAL           lsame
   !*     ..
   !*     .. External Subroutines ..
   !EXTERNAL           dgemv, dger, dswap, dtbsv
   !*     ..
   !*     .. Intrinsic Functions ..
   !*     ..
   !*     .. Executable Statements ..
   !*
   !*     Test the input parameters.
   !*
   info = 0
   notran = .true.
   IF( .NOT.notran ) THEN
      info = -1
   ELSE IF( n.LT.0 ) THEN
      info = -2
   ELSE IF( kl.LT.0 ) THEN
      info = -3
   ELSE IF( ku.LT.0 ) THEN
      info = -4
   ELSE IF( nrhs.LT.0 ) THEN
      info = -5
   ELSE IF( ldab.LT.( 2*kl+ku+1 ) ) THEN
      info = -7
   ELSE IF( ldb.LT.max( 1, n ) ) THEN
      info = -10
   END IF
   IF( info.NE.0 ) THEN
      print *,"error in solving Ax=b",info
      RETURN
   END IF
   !
   !     Quick return if possible
   !
   IF( n.EQ.0 .OR. nrhs.EQ.0 ) RETURN
   !
   kd = ku + kl + 1
   lnoti = kl.GT.0
   !
   IF( notran ) THEN
     !*
     !*        Solve  A*X = B.
     !*
     !*        Solve L*X = B, overwriting B with X.
     !*
     !*        L is represented as a product of permutations and unit lower
     !*        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
     !*        where each transformation L(i) is a rank-one modification of
     !*        the identity matrix.
     !*
      IF( lnoti ) THEN
         DO j = 1, n - 1
            lm = min( kl, n-j )
            l = ipiv( j )
            IF( l.NE.j )  CALL dswap_oacc( nrhs, b( l, 1 ), ldb, b( j, 1 ), ldb )

            CALL dger_oacc( lm, nrhs, -one, ab( kd+1, j ), 1, b( j, 1 ), &
                      ldb, b( j+1, 1 ), ldb )
       END DO
      END IF
      !
      DO i = 1, nrhs
        !*
        !*           Solve U*X = B, overwriting B with X.
        !*
         CALL dtbsv_oacc( 1, 0, 0, n, kl+ku, &
                    ab, ldab, b( 1, i ), 1 )
    END DO
  !!
   END IF
   RETURN
   !*
   !*     End of DGBTRS
   !*
 END subroutine DGBTRS_oacc



subroutine dgbsv_oacc( N, KL, KU, NRHS, AB, LDAB, B, LDB, INFO )
!$acc routine seq
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
  integer :: INFO, KL, KU, LDAB, LDB, N, NRHS
!*     ..
!*     .. Array Arguments ..
  real(r8) :: AB( LDAB, * ), B( LDB, * )
!     ..
!
!  =====================================================================
!
!     .. External Subroutines ..
!  EXTERNAL           dgbtrf, dgbtrs
!     ..
!     .. Intrinsic Functions ..
  INTEGER   :: IPIV(1:n)
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  info = 0
  IF( n.LT.0 ) THEN
     info = -1
  ELSE IF( kl.LT.0 ) THEN
     info = -2
  ELSE IF( ku.LT.0 ) THEN
     info = -3
  ELSE IF( nrhs.LT.0 ) THEN
     info = -4
  ELSE IF( ldab.LT.2*kl+ku+1 ) THEN
     info = -6
  ELSE IF( ldb.LT.max( n, 1 ) ) THEN
     info = -9
  END IF
  IF( info.NE.0 ) THEN
      print *,"error in dgbsv",info
     RETURN
  END IF
!
!     Compute the LU factorization of the band matrix A.
!
  CALL dgbtrf_oacc( n, n, kl, ku, ab, ldab, ipiv, info )
  IF( info.EQ.0 ) THEN
!
!        Solve the system A*X = B, overwriting B with X.
!   I changed the first argument from 'No transpose' to 0
     CALL dgbtrs_oacc( 0, n, kl, ku, nrhs, ab, ldab, ipiv, &
                 b, ldb, info )
  else
    print *,"error: dgbsv",info

  END IF
  RETURN
!
!     End of DGBSV
!
END subroutine dgbsv_oacc

subroutine dgbsv_batch(nbatch, nmax, kl, ku, nrhs, n, ab, b, ipiv, info)
  ! Batched DGBSV: solves A_s * X_s = B_s for s = 1..nbatch, where each A_s is
  ! an n(s) x n(s) band matrix with kl sub- and ku super-diagonals.
  !
  ! Same algorithm and operation order as reference LAPACK DGBSV (unblocked
  ! DGBTF2 + DGBTRS), but the batch index is the leading (fastest) dimension so
  ! that one thread per system gives coalesced memory access on the GPU:
  !   ab(s, kl+ku+1+i-j, j) = A_s(i,j)  for max(1,j-ku) <= i <= min(n(s),j+kl)
  ! Rows 1:kl of ab are workspace for fill-in.  On exit ab holds the LU factors,
  ! ipiv the pivots, and b the solution (only entries 1:n(s) are referenced).
  ! info(s) = 0 success, -i illegal i-th argument, > 0 U(info,info) is exactly zero.
  ! ab, b, ipiv, info must be present on the device when compiled with OpenACC.
  integer , intent(in)    :: nbatch, nmax, kl, ku, nrhs
  integer , intent(in)    :: n(nbatch)
  real(r8), intent(inout) :: ab(nbatch, 2*kl+ku+1, nmax)
  real(r8), intent(inout) :: b(nbatch, nmax, nrhs)
  integer , intent(out)   :: ipiv(nbatch, nmax)
  integer , intent(out)   :: info(nbatch)

  real(r8), parameter :: zero = 0.0_r8, one = 1.0_r8
  integer  :: s, ns, i, j, jj, k, kv, kd, km, jp, ju, l, lm
  real(r8) :: dmax, temp, recip

  kv = ku + kl
  kd = kv + 1

  !$acc parallel loop gang vector default(present) &
  !$acc private(ns, i, j, jj, k, km, jp, ju, l, lm, dmax, temp, recip)
  do s = 1, nbatch
     ns = n(s)
     info(s) = 0
     if (ns < 0 .or. ns > nmax) then
        info(s) = -1
     else if (kl < 0) then
        info(s) = -2
     else if (ku < 0) then
        info(s) = -3
     else if (nrhs < 0) then
        info(s) = -4
     end if
     if (info(s) /= 0 .or. ns == 0) cycle

     ! ---- LU factorization (DGBTF2) ----
     ! Zero fill-in elements in columns ku+2 to kv
     do j = ku + 2, min(kv, ns)
        do i = kv - j + 2, kl
           ab(s,i,j) = zero
        end do
     end do

     ju = 1
     do j = 1, ns
        if (j + kv <= ns) then
           do i = 1, kl
              ab(s,i,j+kv) = zero
           end do
        end if

        ! Pivot search (IDAMAX: first index of max |.|)
        km = min(kl, ns - j)
        jp = 1
        dmax = abs(ab(s,kv+1,j))
        do i = 2, km + 1
           if (abs(ab(s,kv+i,j)) > dmax) then
              jp = i
              dmax = abs(ab(s,kv+i,j))
           end if
        end do
        ipiv(s,j) = jp + j - 1

        if (ab(s,kv+jp,j) /= zero) then
           ju = max(ju, min(j + ku + jp - 1, ns))

           ! Interchange rows j and j+jp-1 in columns j:ju
           if (jp /= 1) then
              do jj = 0, ju - j
                 temp                 = ab(s,kv+jp-jj,j+jj)
                 ab(s,kv+jp-jj,j+jj)  = ab(s,kv+1-jj,j+jj)
                 ab(s,kv+1-jj,j+jj)   = temp
              end do
           end if

           if (km > 0) then
              ! Multipliers
              recip = one / ab(s,kv+1,j)
              do i = 1, km
                 ab(s,kv+1+i,j) = recip * ab(s,kv+1+i,j)
              end do
              ! Rank-1 update of the trailing band
              do jj = 1, ju - j
                 if (ab(s,kv+1-jj,j+jj) /= zero) then
                    temp = -one * ab(s,kv+1-jj,j+jj)
                    do i = 1, km
                       ab(s,kv+1+i-jj,j+jj) = ab(s,kv+1+i-jj,j+jj) + ab(s,kv+1+i,j)*temp
                    end do
                 end if
              end do
           end if
        else if (info(s) == 0) then
           info(s) = j
        end if
     end do

     if (info(s) /= 0) cycle

     ! ---- Solve (DGBTRS, no transpose) ----
     ! L*X = B
     if (kl > 0) then
        do j = 1, ns - 1
           lm = min(kl, ns - j)
           l  = ipiv(s,j)
           if (l /= j) then
              do k = 1, nrhs
                 temp     = b(s,l,k)
                 b(s,l,k) = b(s,j,k)
                 b(s,j,k) = temp
              end do
           end if
           do k = 1, nrhs
              if (b(s,j,k) /= zero) then
                 temp = -one * b(s,j,k)
                 do i = 1, lm
                    b(s,j+i,k) = b(s,j+i,k) + ab(s,kd+i,j)*temp
                 end do
              end if
           end do
        end do
     end if

     ! U*X = B (DTBSV upper, non-unit, bandwidth kl+ku)
     do k = 1, nrhs
        do j = ns, 1, -1
           if (b(s,j,k) /= zero) then
              b(s,j,k) = b(s,j,k) / ab(s,kd,j)
              temp = b(s,j,k)
              do i = j - 1, max(1, j - kv), -1
                 b(s,i,k) = b(s,i,k) - temp*ab(s,kd-j+i,j)
              end do
           end if
        end do
     end do
  end do

end subroutine dgbsv_batch


end module lapack_acc_seq
