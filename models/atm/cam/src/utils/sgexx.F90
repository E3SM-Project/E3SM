module sgexx

! A few LINPACK/LAPACK and supporting BLAS routines

! The routines have all been converted to use the shr_kind_r8 type which
! is set to double precision.  The LINPACK naming convention isn't 
! followed by the 's' routines (which should contain 'd's).

use shr_kind_mod,    only: r8 => shr_kind_r8
use cam_logfile,     only: iulog
use abortutils,      only: endrun

implicit none

public :: &
   dgeco,  &
   dgedi,  &
   sgefa,  &
   sdot,   &
   isamax, &
   sasum,  &
   saxpy,  &
   sscal,  &
   sswap,  &
   sgeev,  &
   sgebak, &
   sgebal, &
   sgehd2, &
   sgehrd, &
   shseqr, &
   slabad, &
   slacpy, &
   sladiv, &
   slahqr, &
   slahrd, &
   slaln2, &
   slange, &
   slanhs, &
   slanv2, &
   slapy2, &
   slarf,  &
   slarfb, &
   slarfg, &
   slarft, &
   slarfx, &
   slartg, &
   slascl, &
   slaset, &
   slassq, &
   sorg2r, &
   sorghr, &
   sorgqr, &
   strevc, &
   ilaenv, &
   lsame,  &
   slamch, &
   slamc1, &
   slamc2, &
   slamc3, &
   slamc4, &
   slamc5, &
   xerbla, &
   scopy,  &
   sger,   &
   snrm2,  &
   srot,   &
   sgemm,  &
   sgemv,  &
   strmm,  &
   strmv

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

      SUBROUTINE DGECO (A,LDA,N,IPVT,RCOND,Z)
!

      INTEGER LDA,N,IPVT(N)
      REAL(r8) A(LDA,*),Z(*)
      REAL(r8) RCOND
!
!     SGECO FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION
!     AND ESTIMATES THE CONDITION OF THE MATRIX.
!
!     IF  RCOND  IS NOT NEEDED, SGEFA IS SLIGHTLY FASTER.
!     TO SOLVE  A*X = B , FOLLOW SGECO BY SGESL.
!     TO COMPUTE  INVERSE(A)*C , FOLLOW SGECO BY SGESL.
!     TO COMPUTE  DETERMINANT(A) , FOLLOW SGECO BY SGEDI.
!     TO COMPUTE  INVERSE(A) , FOLLOW SGECO BY SGEDI.
!
!     ON ENTRY
!
!        A       REAL(LDA, N)
!                THE MATRIX TO BE FACTORED.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!     ON RETURN
!
!        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
!                WHICH WERE USED TO OBTAIN IT.
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
!                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
!                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
!
!        IPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.
!
!        RCOND   REAL
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
!                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
!                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
!                           1.0 + RCOND .EQ. 1.0
!                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
!                UNDERFLOWS.
!
!        Z       REAL(N)
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
!                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     LINPACK SGEFA
!     BLAS SAXPY,SDOT,SSCAL,SASUM
!     FORTRAN ABS,MAX,SIGN
!
!     INTERNAL VARIABLES
!
      REAL(r8) EK,T,WK,WKM
      REAL(r8) ANORM,S,SM,YNORM
      INTEGER INFO,J,K,KB,KP1,L
!
!
!
!     COMPUTE 1-NORM OF A
!
      ANORM = 0.0E0_r8
      DO 10 J = 1, N
         ANORM = MAX(ANORM,SASUM(N,A(1,J),1))
   10 CONTINUE
!
!     FACTOR
!
      CALL SGEFA(A,LDA,N,IPVT,INFO)
!
!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
!     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
!     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
!     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
!     OVERFLOW.
!
!     SOLVE TRANS(U)*W = E
!
      EK = 1.0E0_r8
      DO 20 J = 1, N
         Z(J) = 0.0E0_r8
   20 CONTINUE
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0E0_r8) EK = SIGN(EK,-Z(K))
         IF (ABS(EK-Z(K)) .LE. ABS(A(K,K))) GO TO 30
            S = ABS(A(K,K))/ABS(EK-Z(K))
            CALL SSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = ABS(WK)
         SM = ABS(WKM)
         IF (A(K,K) .EQ. 0.0E0_r8) GO TO 40
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0E0_r8
            WKM = 1.0E0_r8
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               SM = SM + ABS(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + ABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0E0_r8/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
!
!     SOLVE TRANS(L)*Y = W
!
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + SDOT(N-K,A(K+1,K),1,Z(K+1),1)
         IF (ABS(Z(K)) .LE. 1.0E0_r8) GO TO 110
            S = 1.0E0_r8/ABS(Z(K))
            CALL SSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0E0_r8/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
!
      YNORM = 1.0E0_r8
!
!     SOLVE L*V = Y
!
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL SAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
         IF (ABS(Z(K)) .LE. 1.0E0_r8) GO TO 130
            S = 1.0E0_r8/ABS(Z(K))
            CALL SSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0E0_r8/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
      YNORM = S*YNORM
!
!     SOLVE  U*Z = V
!
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (ABS(Z(K)) .LE. ABS(A(K,K))) GO TO 150
            S = ABS(A(K,K))/ABS(Z(K))
            CALL SSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (A(K,K) .NE. 0.0E0_r8) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0E0_r8) Z(K) = 1.0E0_r8
         T = -Z(K)
         CALL SAXPY(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
!     MAKE ZNORM = 1.0
      S = 1.0E0_r8/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
      YNORM = S*YNORM
!
      IF (ANORM .NE. 0.0E0_r8) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0_r8) RCOND = 0.0E0_r8
      RETURN
      END SUBROUTINE DGECO
 
!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE DGEDI (A,LDA,N,IPVT,DET,WORK,JOB)

      INTEGER LDA,N,IPVT(N),JOB
      REAL(r8) A(LDA,*),DET(2),WORK(*)
!
!     SGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX
!     USING THE FACTORS COMPUTED BY SGECO OR SGEFA.
!
!     ON ENTRY
!
!        A       REAL(LDA, N)
!                THE OUTPUT FROM SGECO OR SGEFA.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!        IPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM SGECO OR SGEFA.
!
!        WORK    REAL(N)
!                WORK VECTOR.  CONTENTS DESTROYED.
!
!        JOB     INTEGER
!                = 11   BOTH DETERMINANT AND INVERSE.
!                = 01   INVERSE ONLY.
!                = 10   DETERMINANT ONLY.
!
!     ON RETURN
!
!        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
!                OTHERWISE UNCHANGED.
!
!        DET     REAL(2)
!                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
!                OTHERWISE NOT REFERENCED.
!                DETERMINANT = DET(1) * 10.0**DET(2)
!                WITH  1.0 .LE. ABS(DET(1)) .LT. 10.0
!                OR  DET(1) .EQ. 0.0 .
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
!        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
!        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
!        AND IF SGECO HAS SET RCOND .GT. 0.0 OR SGEFA HAS SET
!        INFO .EQ. 0 .
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS SAXPY,SSCAL,SSWAP
!     FORTRAN ABS,MOD
!
!     INTERNAL VARIABLES
!
      REAL(r8) T
      REAL(r8) TEN
      INTEGER I,J,K,KB,KP1,L,NM1
!
!
!
!     COMPUTE DETERMINANT
!
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0E0_r8
         DET(2) = 0.0E0_r8
         TEN = 10.0E0_r8
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
!        ...EXIT
            IF (DET(1) .EQ. 0.0E0_r8) GO TO 60
   10       IF (ABS(DET(1)) .GE. 1.0E0_r8) GO TO 20
               DET(1) = TEN*DET(1)
               DET(2) = DET(2) - 1.0E0_r8
            GO TO 10
   20       CONTINUE
   30       IF (ABS(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/TEN
               DET(2) = DET(2) + 1.0E0_r8
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
!
!     COMPUTE INVERSE(U)
!
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = 1.0E0_r8/A(K,K)
            T = -A(K,K)
            CALL SSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0E0_r8
               CALL SAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
!
!        FORM INVERSE(U)*INVERSE(L)
!
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = 0.0E0_r8
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
               CALL SAXPY(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) CALL SSWAP(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END SUBROUTINE DGEDI

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SGEFA (A,LDA,N,IPVT,INFO)

      INTEGER LDA,N,IPVT(N),INFO
      REAL(r8) A(LDA,*)
!
!     SGEFA FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.
!
!     SGEFA IS USUALLY CALLED BY SGECO, BUT IT CAN BE CALLED
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
!     (TIME FOR SGECO) = (1 + 9/N)*(TIME FOR SGEFA) .
!
!     ON ENTRY
!
!        A       REAL(LDA, N)
!                THE MATRIX TO BE FACTORED.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!     ON RETURN
!
!        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
!                WHICH WERE USED TO OBTAIN IT.
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
!                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
!                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
!
!        IPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.
!
!        INFO    INTEGER
!                = 0  NORMAL VALUE.
!                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
!                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
!                     INDICATE THAT SGESL OR SGEDI WILL DIVIDE BY ZERO
!                     IF CALLED.  USE  RCOND  IN SGECO FOR A RELIABLE
!                     INDICATION OF SINGULARITY.
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS SAXPY,SSCAL,ISAMAX
!
!     INTERNAL VARIABLES
!
      REAL(r8) T
      INTEGER J,K,KP1,L,NM1
!
!
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
!
!        FIND L = PIVOT INDEX
!
         L = ISAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
!
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!
         IF (A(L,K) .EQ. 0.0E0_r8) GO TO 40
!
!           INTERCHANGE IF NECESSARY
!
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
!
!           COMPUTE MULTIPLIERS
!
            T = -1.0E0_r8/A(K,K)
            CALL SSCAL(N-K,T,A(K+1,K),1)
!
!           ROW ELIMINATION WITH COLUMN INDEXING
!
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL SAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0E0_r8) INFO = N
      RETURN
      END SUBROUTINE SGEFA

!--------1---------2---------3---------4---------5---------6---------7---------8
 
      FUNCTION SDOT(N,SX,INCX,SY,INCY)
!
!     FORMS THE DOT PRODUCT OF TWO VECTORS.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
!
      real(r8) sdot

      REAL(r8) SX(*),SY(*),STEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
!
      STEMP = 0.0E0_r8
      SDOT = 0.0E0_r8
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
!
!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        STEMP = STEMP + SX(IX)*SY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      SDOT = STEMP
      RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        STEMP = STEMP + SX(I)*SY(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        STEMP = STEMP + SX(I)*SY(I) + SX(I + 1)*SY(I + 1) + &
         SX(I + 2)*SY(I + 2) + SX(I + 3)*SY(I + 3) + SX(I + 4)*SY(I + 4)
   50 CONTINUE
   60 SDOT = STEMP
      RETURN
      END FUNCTION SDOT

!--------1---------2---------3---------4---------5---------6---------7---------8

      INTEGER FUNCTION ISAMAX(N,SX,INCX)
!
!     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
      REAL(r8) SX(*),SMAX
      INTEGER I,INCX,IX,N
!
      ISAMAX = 0
      IF( N .LT. 1 ) RETURN
      ISAMAX = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
!
!        CODE FOR INCREMENT NOT EQUAL TO 1
!
      IX = 1
      SMAX = ABS(SX(1))
      IX = IX + INCX
      DO 10 I = 2,N
         IF(ABS(SX(IX)).LE.SMAX) GO TO 5
         ISAMAX = I
         SMAX = ABS(SX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
!
!        CODE FOR INCREMENT EQUAL TO 1
!
   20 SMAX = ABS(SX(1))
      DO 30 I = 2,N
         IF(ABS(SX(I)).LE.SMAX) GO TO 30
         ISAMAX = I
         SMAX = ABS(SX(I))
   30 CONTINUE
      RETURN
      END FUNCTION ISAMAX
 
!--------1---------2---------3---------4---------5---------6---------7---------8

      FUNCTION SASUM(N,SX,INCX)
!
!     TAKES THE SUM OF THE ABSOLUTE VALUES.
!     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
      real(r8) sasum

      REAL(r8) SX(*),STEMP
      INTEGER I,INCX,M,MP1,N,NINCX
!
      SASUM = 0.0E0_r8
      STEMP = 0.0E0_r8
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
!
!        CODE FOR INCREMENT NOT EQUAL TO 1
!
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        STEMP = STEMP + ABS(SX(I))
   10 CONTINUE
      SASUM = STEMP
      RETURN
!
!        CODE FOR INCREMENT EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        STEMP = STEMP + ABS(SX(I))
   30 CONTINUE
      IF( N .LT. 6 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        STEMP = STEMP + ABS(SX(I)) + ABS(SX(I + 1)) + ABS(SX(I + 2)) &
        + ABS(SX(I + 3)) + ABS(SX(I + 4)) + ABS(SX(I + 5))
   50 CONTINUE
   60 SASUM = STEMP
      RETURN
      END FUNCTION SASUM
 
!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
!
!     CONSTANT TIMES A VECTOR PLUS A VECTOR.
!     USES UNROLLED LOOP FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
      REAL(r8) SX(*),SY(*),SA
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
!
      IF(N.LE.0)RETURN
      IF (SA .EQ. 0.0_r8) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
!
!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SY(IY) + SA*SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        SY(I) = SY(I) + SA*SX(I)
        SY(I + 1) = SY(I + 1) + SA*SX(I + 1)
        SY(I + 2) = SY(I + 2) + SA*SX(I + 2)
        SY(I + 3) = SY(I + 3) + SA*SX(I + 3)
   50 CONTINUE
      RETURN
      END SUBROUTINE SAXPY
 
!--------1---------2---------3---------4---------5---------6---------7---------8

      subroutine sscal(n,sa,sx,incx)
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to 1.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      integer i,incx,m,mp1,n,nincx
      real(r8) sa,sx(*)
!
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      do 10 i = 1,nincx,incx
        sx(i) = sa*sx(i)
   10 continue
      return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sx(i) = sa*sx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        sx(i) = sa*sx(i)
        sx(i + 1) = sa*sx(i + 1)
        sx(i + 2) = sa*sx(i + 2)
        sx(i + 3) = sa*sx(i + 3)
        sx(i + 4) = sa*sx(i + 4)
   50 continue
      return
      end subroutine sscal

!--------1---------2---------3---------4---------5---------6---------7---------8

      subroutine sswap (n,sx,incx,sy,incy)
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal to 1.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      real(r8) sx(*),sy(*),stemp
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = sx(ix)
        sx(ix) = sy(iy)
        sy(iy) = stemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
        stemp = sx(i + 1)
        sx(i + 1) = sy(i + 1)
        sy(i + 1) = stemp
        stemp = sx(i + 2)
        sx(i + 2) = sy(i + 2)
        sy(i + 2) = stemp
   50 continue
      return
      end subroutine sswap
 
!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, &
                        LDVR, WORK, LWORK, INFO )

!
!  -- LAPACK driver routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
!     ..
!     .. Array Arguments ..
      REAL(r8)               A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), &
                         WI( * ), WORK( * ), WR( * )
!     ..
!
!  Purpose
!  =======
!
!  SGEEV computes for an N-by-N real nonsymmetric matrix A, the
!  eigenvalues and, optionally, the left and/or right eigenvectors.
!
!  The right eigenvector v(j) of A satisfies
!                   A * v(j) = lambda(j) * v(j)
!  where lambda(j) is its eigenvalue.
!  The left eigenvector u(j) of A satisfies
!                u(j)**H * A = lambda(j) * u(j)**H
!  where u(j)**H denotes the conjugate transpose of u(j).
!
!  The computed eigenvectors are normalized to have Euclidean norm
!  equal to 1 and largest component real.
!
!  Arguments
!  =========
!
!  JOBVL   (input) CHARACTER*1
!          = 'N': left eigenvectors of A are not computed;
!          = 'V': left eigenvectors of A are computed.
!
!  JOBVR   (input) CHARACTER*1
!          = 'N': right eigenvectors of A are not computed;
!          = 'V': right eigenvectors of A are computed.
!
!  N       (input) INTEGER
!          The order of the matrix A. N >= 0.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the N-by-N matrix A.
!          On exit, A has been overwritten.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  WR      (output) REAL array, dimension (N)
!  WI      (output) REAL array, dimension (N)
!          WR and WI contain the real and imaginary parts,
!          respectively, of the computed eigenvalues.  Complex
!          conjugate pairs of eigenvalues appear consecutively
!          with the eigenvalue having the positive imaginary part
!          first.
!
!  VL      (output) REAL array, dimension (LDVL,N)
!          If JOBVL = 'V', the left eigenvectors u(j) are stored one
!          after another in the columns of VL, in the same order
!          as their eigenvalues.
!          If JOBVL = 'N', VL is not referenced.
!          If the j-th eigenvalue is real, then u(j) = VL(:,j),
!          the j-th column of VL.
!          If the j-th and (j+1)-st eigenvalues form a complex
!          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
!          u(j+1) = VL(:,j) - i*VL(:,j+1).
!
!  LDVL    (input) INTEGER
!          The leading dimension of the array VL.  LDVL >= 1; if
!          JOBVL = 'V', LDVL >= N.
!
!  VR      (output) REAL array, dimension (LDVR,N)
!          If JOBVR = 'V', the right eigenvectors v(j) are stored one
!          after another in the columns of VR, in the same order
!          as their eigenvalues.
!          If JOBVR = 'N', VR is not referenced.
!          If the j-th eigenvalue is real, then v(j) = VR(:,j),
!          the j-th column of VR.
!          If the j-th and (j+1)-st eigenvalues form a complex
!          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
!          v(j+1) = VR(:,j) - i*VR(:,j+1).
!
!  LDVR    (input) INTEGER
!          The leading dimension of the array VR.  LDVR >= 1; if
!          JOBVR = 'V', LDVR >= N.
!
!  WORK    (workspace/output) REAL array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,3*N), and
!          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good
!          performance, LWORK must generally be larger.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = i, the QR algorithm failed to compute all the
!                eigenvalues, and no eigenvectors have been computed;
!                elements i+1:N of WR and WI contain eigenvalues which
!                have converged.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0_r8, ONE = 1.0E0_r8 )
!     ..
!     .. Local Scalars ..
      LOGICAL            SCALEA, WANTVL, WANTVR
      CHARACTER          SIDE
      INTEGER            HSWORK, I, IBAL, IERR, IHI, ILO, ITAU, IWRK, K, &
                         MAXB, MAXWRK, MINWRK, NOUT
      REAL(r8)               ANRM, BIGNUM, CS, CSCALE, EPS, R, SCL, SMLNUM, &
                         SN
!     ..
!     .. Local Arrays ..
      LOGICAL            SELECT( 1 )
      REAL(r8)               DUM( 1 )
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      WANTVL = LSAME( JOBVL, 'V' )
      WANTVR = LSAME( JOBVR, 'V' )
      IF( ( .NOT.WANTVL ) .AND. ( .NOT.LSAME( JOBVL, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( ( .NOT.WANTVR ) .AND. ( .NOT.LSAME( JOBVR, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDVL.LT.1 .OR. ( WANTVL .AND. LDVL.LT.N ) ) THEN
         INFO = -9
      ELSE IF( LDVR.LT.1 .OR. ( WANTVR .AND. LDVR.LT.N ) ) THEN
         INFO = -11
      END IF
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of workspace needed at that point in the code,
!       as well as the preferred amount for good performance.
!       NB refers to the optimal block size for the immediately
!       following subroutine, as returned by ILAENV.
!       HSWORK refers to the workspace preferred by SHSEQR, as
!       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
!       the worst case.)
!
      MINWRK = 1
      IF( INFO.EQ.0 .AND. LWORK.GE.1 ) THEN
         MAXWRK = 2*N + N*ILAENV( 1, 'SGEHRD', ' ', N, 1, N, 0 )
         IF( ( .NOT.WANTVL ) .AND. ( .NOT.WANTVR ) ) THEN
            MINWRK = MAX( 1, 3*N )
            MAXB = MAX( ILAENV( 8, 'SHSEQR', 'EN', N, 1, N, -1 ), 2 )
            K = MIN( MAXB, N, MAX( 2, ILAENV( 4, 'SHSEQR', 'EN', N, 1, &
                N, -1 ) ) )
            HSWORK = MAX( K*( K+2 ), 2*N )
            MAXWRK = MAX( MAXWRK, N+1, N+HSWORK )
         ELSE
            MINWRK = MAX( 1, 4*N )
            MAXWRK = MAX( MAXWRK, 2*N+( N-1 )* &
                     ILAENV( 1, 'SORGHR', ' ', N, 1, N, -1 ) )
            MAXB = MAX( ILAENV( 8, 'SHSEQR', 'SV', N, 1, N, -1 ), 2 )
            K = MIN( MAXB, N, MAX( 2, ILAENV( 4, 'SHSEQR', 'SV', N, 1, &
                N, -1 ) ) )
            HSWORK = MAX( K*( K+2 ), 2*N )
            MAXWRK = MAX( MAXWRK, N+1, N+HSWORK )
            MAXWRK = MAX( MAXWRK, 4*N )
         END IF
         WORK( 1 ) = MAXWRK
      END IF
      IF( LWORK.LT.MINWRK ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEEV ', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
!     Get machine constants
!
      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL SLABAD( SMLNUM, BIGNUM )
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      ANRM = SLANGE( 'M', N, N, A, LDA, DUM )
      SCALEA = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      END IF
      IF( SCALEA ) &
         CALL SLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )
!
!     Balance the matrix
!     (Workspace: need N)
!
      IBAL = 1
      CALL SGEBAL( 'B', N, A, LDA, ILO, IHI, WORK( IBAL ), IERR )
!
!     Reduce to upper Hessenberg form
!     (Workspace: need 3*N, prefer 2*N+N*NB)
!
      ITAU = IBAL + N
      IWRK = ITAU + N
      CALL SGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), &
                   LWORK-IWRK+1, IERR )
!
      IF( WANTVL ) THEN
!
!        Want left eigenvectors
!        Copy Householder vectors to VL
!
         SIDE = 'L'
         CALL SLACPY( 'L', N, N, A, LDA, VL, LDVL )
!
!        Generate orthogonal matrix in VL
!        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!
         CALL SORGHR( N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ), &
                      LWORK-IWRK+1, IERR )
!
!        Perform QR iteration, accumulating Schur vectors in VL
!        (Workspace: need N+1, prefer N+HSWORK (see comments) )
!
         IWRK = ITAU
         CALL SHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VL, LDVL, &
                      WORK( IWRK ), LWORK-IWRK+1, INFO )
!
         IF( WANTVR ) THEN
!
!           Want left and right eigenvectors
!           Copy Schur vectors to VR
!
            SIDE = 'B'
            CALL SLACPY( 'F', N, N, VL, LDVL, VR, LDVR )
         END IF
!
      ELSE IF( WANTVR ) THEN
!
!        Want right eigenvectors
!        Copy Householder vectors to VR
!
         SIDE = 'R'
         CALL SLACPY( 'L', N, N, A, LDA, VR, LDVR )
!
!        Generate orthogonal matrix in VR
!        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!
         CALL SORGHR( N, ILO, IHI, VR, LDVR, WORK( ITAU ), WORK( IWRK ), &
                      LWORK-IWRK+1, IERR )
!
!        Perform QR iteration, accumulating Schur vectors in VR
!        (Workspace: need N+1, prefer N+HSWORK (see comments) )
!
         IWRK = ITAU
         CALL SHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR, &
                      WORK( IWRK ), LWORK-IWRK+1, INFO )
!
      ELSE
!
!        Compute eigenvalues only
!        (Workspace: need N+1, prefer N+HSWORK (see comments) )
!
         IWRK = ITAU
         CALL SHSEQR( 'E', 'N', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR, &
                      WORK( IWRK ), LWORK-IWRK+1, INFO )
      END IF
!
!     If INFO > 0 from SHSEQR, then quit
!
      IF( INFO.GT.0 ) &
         GO TO 50
!
      IF( WANTVL .OR. WANTVR ) THEN
!
!        Compute left and/or right eigenvectors
!        (Workspace: need 4*N)
!
         CALL STREVC( SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, &
                      N, NOUT, WORK( IWRK ), IERR )
      END IF
!
      IF( WANTVL ) THEN
!
!        Undo balancing of left eigenvectors
!        (Workspace: need N)
!
         CALL SGEBAK( 'B', 'L', N, ILO, IHI, WORK( IBAL ), N, VL, LDVL, &
                      IERR )
!
!        Normalize left eigenvectors and make largest component real
!
         DO 20 I = 1, N
            IF( WI( I ).EQ.ZERO ) THEN
               SCL = ONE / SNRM2( N, VL( 1, I ), 1 )
               CALL SSCAL( N, SCL, VL( 1, I ), 1 )
            ELSE IF( WI( I ).GT.ZERO ) THEN
               SCL = ONE / SLAPY2( SNRM2( N, VL( 1, I ), 1 ), &
                     SNRM2( N, VL( 1, I+1 ), 1 ) )
               CALL SSCAL( N, SCL, VL( 1, I ), 1 )
               CALL SSCAL( N, SCL, VL( 1, I+1 ), 1 )
               DO 10 K = 1, N
                  WORK( IWRK+K-1 ) = VL( K, I )**2 + VL( K, I+1 )**2
   10          CONTINUE
               K = ISAMAX( N, WORK( IWRK ), 1 )
               CALL SLARTG( VL( K, I ), VL( K, I+1 ), CS, SN, R )
               CALL SROT( N, VL( 1, I ), 1, VL( 1, I+1 ), 1, CS, SN )
               VL( K, I+1 ) = ZERO
            END IF
   20    CONTINUE
      END IF
!
      IF( WANTVR ) THEN
!
!        Undo balancing of right eigenvectors
!        (Workspace: need N)
!
         CALL SGEBAK( 'B', 'R', N, ILO, IHI, WORK( IBAL ), N, VR, LDVR, &
                      IERR )
!
!        Normalize right eigenvectors and make largest component real
!
         DO 40 I = 1, N
            IF( WI( I ).EQ.ZERO ) THEN
               SCL = ONE / SNRM2( N, VR( 1, I ), 1 )
               CALL SSCAL( N, SCL, VR( 1, I ), 1 )
            ELSE IF( WI( I ).GT.ZERO ) THEN
               SCL = ONE / SLAPY2( SNRM2( N, VR( 1, I ), 1 ), &
                     SNRM2( N, VR( 1, I+1 ), 1 ) )
               CALL SSCAL( N, SCL, VR( 1, I ), 1 )
               CALL SSCAL( N, SCL, VR( 1, I+1 ), 1 )
               DO 30 K = 1, N
                  WORK( IWRK+K-1 ) = VR( K, I )**2 + VR( K, I+1 )**2
   30          CONTINUE
               K = ISAMAX( N, WORK( IWRK ), 1 )
               CALL SLARTG( VR( K, I ), VR( K, I+1 ), CS, SN, R )
               CALL SROT( N, VR( 1, I ), 1, VR( 1, I+1 ), 1, CS, SN )
               VR( K, I+1 ) = ZERO
            END IF
   40    CONTINUE
      END IF
!
!     Undo scaling if necessary
!
   50 CONTINUE
      IF( SCALEA ) THEN
         CALL SLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WR( INFO+1 ), &
                      MAX( N-INFO, 1 ), IERR )
         CALL SLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WI( INFO+1 ), &
                      MAX( N-INFO, 1 ), IERR )
         IF( INFO.GT.0 ) THEN
            CALL SLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WR, N, &
                         IERR )
            CALL SLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WI, N, &
                         IERR )
         END IF
      END IF
!
      WORK( 1 ) = MAXWRK
      RETURN
!
!     End of SGEEV
!
      END SUBROUTINE SGEEV
 
!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV, &
                         INFO )

!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          JOB, SIDE
      INTEGER            IHI, ILO, INFO, LDV, M, N
!     ..
!     .. Array Arguments ..
      REAL(r8)               V( LDV, * ), SCALE( * )
!     ..
!
!  Purpose
!  =======
!
!  SGEBAK forms the right or left eigenvectors of a real general matrix
!  by backward transformation on the computed eigenvectors of the
!  balanced matrix output by SGEBAL.
!
!  Arguments
!  =========
!
!  JOB     (input) CHARACTER*1
!          Specifies the type of backward transformation required:
!          = 'N', do nothing, return immediately;
!          = 'P', do backward transformation for permutation only;
!          = 'S', do backward transformation for scaling only;
!          = 'B', do backward transformations for both permutation and
!                 scaling.
!          JOB must be the same as the argument JOB supplied to SGEBAL.
!
!  SIDE    (input) CHARACTER*1
!          = 'R':  V contains right eigenvectors;
!          = 'L':  V contains left eigenvectors.
!
!  N       (input) INTEGER
!          The number of rows of the matrix V.  N >= 0.
!
!  ILO     (input) INTEGER
!  IHI     (input) INTEGER
!          The integers ILO and IHI determined by SGEBAL.
!          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!
!  SCALE   (input) REAL array, dimension (N)
!          Details of the permutation and scaling factors, as returned
!          by SGEBAL.
!
!  M       (input) INTEGER
!          The number of columns of the matrix V.  M >= 0.
!
!  V       (input/output) REAL array, dimension (LDV,M)
!          On entry, the matrix of right or left eigenvectors to be
!          transformed, as returned by SHSEIN or STREVC.
!          On exit, V is overwritten by the transformed eigenvectors.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V. LDV >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ONE
      PARAMETER          ( ONE = 1.0E+0_r8 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFTV, RIGHTV
      INTEGER            I, II, K
      REAL(r8)               S
!     ..
!     .. Executable Statements ..
!
!     Decode and Test the input parameters
!
      RIGHTV = LSAME( SIDE, 'R' )
      LEFTV = LSAME( SIDE, 'L' )
!
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND. &
          .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -7
      ELSE IF( LDV.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEBAK', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
      IF( M.EQ.0 ) &
         RETURN
      IF( LSAME( JOB, 'N' ) ) &
         RETURN
!
      IF( ILO.EQ.IHI ) &
         GO TO 30
!
!     Backward balance
!
      IF( LSAME( JOB, 'S' ) .OR. LSAME( JOB, 'B' ) ) THEN
!
         IF( RIGHTV ) THEN
            DO 10 I = ILO, IHI
               S = SCALE( I )
               CALL SSCAL( M, S, V( I, 1 ), LDV )
   10       CONTINUE
         END IF
!
         IF( LEFTV ) THEN
            DO 20 I = ILO, IHI
               S = ONE / SCALE( I )
               CALL SSCAL( M, S, V( I, 1 ), LDV )
   20       CONTINUE
         END IF
!
      END IF
!
!     Backward permutation
!
!     For  I = ILO-1 step -1 until 1,
!              IHI+1 step 1 until N do --
!
   30 CONTINUE
      IF( LSAME( JOB, 'P' ) .OR. LSAME( JOB, 'B' ) ) THEN
         IF( RIGHTV ) THEN
            DO 40 II = 1, N
               I = II
               IF( I.GE.ILO .AND. I.LE.IHI ) &
                  GO TO 40
               IF( I.LT.ILO ) &
                  I = ILO - II
               K = SCALE( I )
               IF( K.EQ.I ) &
                  GO TO 40
               CALL SSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   40       CONTINUE
         END IF
!
         IF( LEFTV ) THEN
            DO 50 II = 1, N
               I = II
               IF( I.GE.ILO .AND. I.LE.IHI ) &
                  GO TO 50
               IF( I.LT.ILO ) &
                  I = ILO - II
               K = SCALE( I )
               IF( K.EQ.I ) &
                  GO TO 50
               CALL SSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   50       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of SGEBAK 
!
      END SUBROUTINE SGEBAK
 
!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )

!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            IHI, ILO, INFO, LDA, N
!     ..
!     .. Array Arguments ..
      REAL(r8)               A( LDA, * ), SCALE( * )
!     ..
!
!  Purpose
!  =======
!
!  SGEBAL balances a general real matrix A.  This involves, first,
!  permuting A by a similarity transformation to isolate eigenvalues
!  in the first 1 to ILO-1 and last IHI+1 to N elements on the
!  diagonal; and second, applying a diagonal similarity transformation
!  to rows and columns ILO to IHI to make the rows and columns as
!  close in norm as possible.  Both steps are optional.
!
!  Balancing may reduce the 1-norm of the matrix, and improve the
!  accuracy of the computed eigenvalues and/or eigenvectors.
!
!  Arguments
!  =========
!
!  JOB     (input) CHARACTER*1
!          Specifies the operations to be performed on A:
!          = 'N':  none:  simply set ILO = 1, IHI = N, SCALE(I) = 1.0
!                  for i = 1,...,N;
!          = 'P':  permute only;
!          = 'S':  scale only;
!          = 'B':  both permute and scale.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the input matrix A.
!          On exit,  A is overwritten by the balanced matrix.
!          If JOB = 'N', A is not referenced.
!          See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  ILO     (output) INTEGER
!  IHI     (output) INTEGER
!          ILO and IHI are set to integers such that on exit
!          A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N.
!          If JOB = 'N' or 'S', ILO = 1 and IHI = N.
!
!  SCALE   (output) REAL array, dimension (N)
!          Details of the permutations and scaling factors applied to
!          A.  If P(j) is the index of the row and column interchanged
!          with row and column j and D(j) is the scaling factor
!          applied to row and column j, then
!          SCALE(j) = P(j)    for j = 1,...,ILO-1
!                   = D(j)    for j = ILO,...,IHI
!                   = P(j)    for j = IHI+1,...,N.
!          The order in which the interchanges are made is N to IHI+1,
!          then 1 to ILO-1.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The permutations consist of row and column interchanges which put
!  the matrix in the form
!
!             ( T1   X   Y  )
!     P A P = (  0   B   Z  )
!             (  0   0   T2 )
!
!  where T1 and T2 are upper triangular matrices whose eigenvalues lie
!  along the diagonal.  The column indices ILO and IHI mark the starting
!  and ending columns of the submatrix B. Balancing consists of applying
!  a diagonal similarity transformation inv(D) * B * D to make the
!  1-norms of each row of B and its corresponding column nearly equal.
!  The output matrix is
!
!     ( T1     X*D          Y    )
!     (  0  inv(D)*B*D  inv(D)*Z ).
!     (  0      0           T2   )
!
!  Information about the permutations P and the diagonal matrix D is
!  returned in the vector SCALE.
!
!  This subroutine is based on the EISPACK routine BALANC.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0_r8, ONE = 1.0E+0_r8 )
      REAL(r8)               SCLFAC
      PARAMETER          ( SCLFAC = 1.0E+1_r8 )
      REAL(r8)               FACTOR
      PARAMETER          ( FACTOR = 0.95E+0_r8 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOCONV
      INTEGER            I, ICA, IEXC, IRA, J, K, L, M
      REAL(r8)               C, CA, F, G, R, RA, S, SFMAX1, SFMAX2, SFMIN1, &
                         SFMIN2
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND. &
          .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEBAL', -INFO )
         RETURN
      END IF
!
      K = 1
      L = N
!
      IF( N.EQ.0 ) &
         GO TO 210
!
      IF( LSAME( JOB, 'N' ) ) THEN
         DO 10 I = 1, N
            SCALE( I ) = ONE
   10    CONTINUE
         GO TO 210
      END IF
!
      IF( LSAME( JOB, 'S' ) ) &
         GO TO 120
!
!     Permutation to isolate eigenvalues if possible
!
      GO TO 50
!
!     Row and column exchange.
!
   20 CONTINUE
      SCALE( M ) = J
      IF( J.EQ.M ) &
         GO TO 30
!
      CALL SSWAP( L, A( 1, J ), 1, A( 1, M ), 1 )
      CALL SSWAP( N-K+1, A( J, K ), LDA, A( M, K ), LDA )
!
   30 CONTINUE
      GO TO ( 40, 80 )IEXC
!
!     Search for rows isolating an eigenvalue and push them down.
!
   40 CONTINUE
      IF( L.EQ.1 ) &
         GO TO 210
      L = L - 1
!
   50 CONTINUE
      DO 70 J = L, 1, -1
!
         DO 60 I = 1, L
            IF( I.EQ.J ) &
               GO TO 60
            IF( A( J, I ).NE.ZERO ) &
               GO TO 70
   60    CONTINUE
!
         M = L
         IEXC = 1
         GO TO 20
   70 CONTINUE
!
      GO TO 90
!
!     Search for columns isolating an eigenvalue and push them left.
!
   80 CONTINUE
      K = K + 1
!
   90 CONTINUE
      DO 110 J = K, L
!
         DO 100 I = K, L
            IF( I.EQ.J ) &
               GO TO 100
            IF( A( I, J ).NE.ZERO ) &
               GO TO 110
  100    CONTINUE
!
         M = K
         IEXC = 2
         GO TO 20
  110 CONTINUE
!
  120 CONTINUE
      DO 130 I = K, L
         SCALE( I ) = ONE
  130 CONTINUE
!
      IF( LSAME( JOB, 'P' ) ) &
         GO TO 210
!
!     Balance the submatrix in rows K to L.
!
!     Iterative loop for norm reduction
!
      SFMIN1 = SLAMCH( 'S' ) / SLAMCH( 'P' )
      SFMAX1 = ONE / SFMIN1
      SFMIN2 = SFMIN1*SCLFAC
      SFMAX2 = ONE / SFMIN2
  140 CONTINUE
      NOCONV = .FALSE.
!
      DO 200 I = K, L
         C = ZERO
         R = ZERO
!
         DO 150 J = K, L
            IF( J.EQ.I ) &
               GO TO 150
            C = C + ABS( A( J, I ) )
            R = R + ABS( A( I, J ) )
  150    CONTINUE
         ICA = ISAMAX( L, A( 1, I ), 1 )
         CA = ABS( A( ICA, I ) )
         IRA = ISAMAX( N-K+1, A( I, K ), LDA )
         RA = ABS( A( I, IRA+K-1 ) )
!
!        Guard against zero C or R due to underflow.
!
         IF( C.EQ.ZERO .OR. R.EQ.ZERO ) &
            GO TO 200
         G = R / SCLFAC
         F = ONE
         S = C + R
  160    CONTINUE
         IF( C.GE.G .OR. MAX( F, C, CA ).GE.SFMAX2 .OR. &
             MIN( R, G, RA ).LE.SFMIN2 )GO TO 170
         F = F*SCLFAC
         C = C*SCLFAC
         CA = CA*SCLFAC
         R = R / SCLFAC
         G = G / SCLFAC
         RA = RA / SCLFAC
         GO TO 160
!
  170    CONTINUE
         G = C / SCLFAC
  180    CONTINUE
         IF( G.LT.R .OR. MAX( R, RA ).GE.SFMAX2 .OR. &
             MIN( F, C, G, CA ).LE.SFMIN2 )GO TO 190
         F = F / SCLFAC
         C = C / SCLFAC
         G = G / SCLFAC
         CA = CA / SCLFAC
         R = R*SCLFAC
         RA = RA*SCLFAC
         GO TO 180
!
!        Now balance.
!
  190    CONTINUE
         IF( ( C+R ).GE.FACTOR*S ) &
            GO TO 200
         IF( F.LT.ONE .AND. SCALE( I ).LT.ONE ) THEN
            IF( F*SCALE( I ).LE.SFMIN1 ) &
               GO TO 200
         END IF
         IF( F.GT.ONE .AND. SCALE( I ).GT.ONE ) THEN
            IF( SCALE( I ).GE.SFMAX1 / F ) &
               GO TO 200
         END IF
         G = ONE / F
         SCALE( I ) = SCALE( I )*F
         NOCONV = .TRUE.
!
         CALL SSCAL( N-K+1, G, A( I, K ), LDA )
         CALL SSCAL( L, F, A( 1, I ), 1 )
!
  200 CONTINUE
!
      IF( NOCONV ) &
         GO TO 140
!
  210 CONTINUE
      ILO = K
      IHI = L
!
      RETURN
!
!     End of SGEBAL
!
      END SUBROUTINE SGEBAL
 
!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )

!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, N
!     ..
!     .. Array Arguments ..
      REAL(r8)               A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  SGEHD2 reduces a real general matrix A to upper Hessenberg form H by
!  an orthogonal similarity transformation:  Q' * A * Q = H .
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  ILO     (input) INTEGER
!  IHI     (input) INTEGER
!          It is assumed that A is already upper triangular in rows
!          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
!          set by a previous call to SGEBAL; otherwise they should be
!          set to 1 and N respectively. See Further Details.
!          1 <= ILO <= IHI <= max(1,N).
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the n by n general matrix to be reduced.
!          On exit, the upper triangle and the first subdiagonal of A
!          are overwritten with the upper Hessenberg matrix H, and the
!          elements below the first subdiagonal, with the array TAU,
!          represent the orthogonal matrix Q as a product of elementary
!          reflectors. See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  TAU     (output) REAL array, dimension (N-1)
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace) REAL array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of (ihi-ilo) elementary
!  reflectors
!
!     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
!  exit in A(i+2:ihi,i), and tau in TAU(i).
!
!  The contents of A are illustrated by the following example, with
!  n = 7, ilo = 2 and ihi = 6:
!
!  on entry,                        on exit,
!
!  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
!  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
!  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
!  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
!  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
!  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
!  (                         a )    (                          a )
!
!  where a denotes an element of the original matrix A, h denotes a
!  modified element of the upper Hessenberg matrix H, and vi denotes an
!  element of the vector defining H(i).
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ONE
      PARAMETER          ( ONE = 1.0E+0_r8 )
!     ..
!     .. Local Scalars ..
      INTEGER            I
      REAL(r8)               AII
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEHD2', -INFO )
         RETURN
      END IF
!
      DO 10 I = ILO, IHI - 1
!
!        Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
!
         CALL SLARFG( IHI-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1, &
                      TAU( I ) )
         AII = A( I+1, I )
         A( I+1, I ) = ONE
!
!        Apply H(i) to A(1:ihi,i+1:ihi) from the right
!
         CALL SLARF( 'Right', IHI, IHI-I, A( I+1, I ), 1, TAU( I ), &
                     A( 1, I+1 ), LDA, WORK )
!
!        Apply H(i) to A(i+1:ihi,i+1:n) from the left
!
         CALL SLARF( 'Left', IHI-I, N-I, A( I+1, I ), 1, TAU( I ), &
                     A( I+1, I+1 ), LDA, WORK )
!
         A( I+1, I ) = AII
   10 CONTINUE
!
      RETURN
!
!     End of SGEHD2
!
      END SUBROUTINE SGEHD2

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )

!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      REAL(r8)               A( LDA, * ), TAU( * ), WORK( LWORK )
!     ..
!
!  Purpose
!  =======
!
!  SGEHRD reduces a real general matrix A to upper Hessenberg form H by
!  an orthogonal similarity transformation:  Q' * A * Q = H .
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  ILO     (input) INTEGER
!  IHI     (input) INTEGER
!          It is assumed that A is already upper triangular in rows
!          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
!          set by a previous call to SGEBAL; otherwise they should be
!          set to 1 and N respectively. See Further Details.
!          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the N-by-N general matrix to be reduced.
!          On exit, the upper triangle and the first subdiagonal of A
!          are overwritten with the upper Hessenberg matrix H, and the
!          elements below the first subdiagonal, with the array TAU,
!          represent the orthogonal matrix Q as a product of elementary
!          reflectors. See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  TAU     (output) REAL array, dimension (N-1)
!          The scalar factors of the elementary reflectors (see Further
!          Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to
!          zero.
!
!  WORK    (workspace/output) REAL array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The length of the array WORK.  LWORK >= max(1,N).
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of (ihi-ilo) elementary
!  reflectors
!
!     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
!  exit in A(i+2:ihi,i), and tau in TAU(i).
!
!  The contents of A are illustrated by the following example, with
!  n = 7, ilo = 2 and ihi = 6:
!
!  on entry,                        on exit,
!
!  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
!  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
!  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
!  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
!  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
!  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
!  (                         a )    (                          a )
!
!  where a denotes an element of the original matrix A, h denotes a
!  modified element of the upper Hessenberg matrix H, and vi denotes an
!  element of the vector defining H(i).
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
      REAL(r8)               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0_r8, ONE = 1.0E+0_r8 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IB, IINFO, IWS, LDWORK, NB, NBMIN, NH, NX
      REAL(r8)               EI
!     ..
!     .. Local Arrays ..
      REAL(r8)               T( LDT, NBMAX )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEHRD', -INFO )
         RETURN
      END IF
!
!     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero
!
      DO 10 I = 1, ILO - 1
         TAU( I ) = ZERO
   10 CONTINUE
      DO 20 I = MAX( 1, IHI ), N - 1
         TAU( I ) = ZERO
   20 CONTINUE
!
!     Quick return if possible
!
      NH = IHI - ILO + 1
      IF( NH.LE.1 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
!     Determine the block size.
!
      NB = MIN( NBMAX, ILAENV( 1, 'SGEHRD', ' ', N, ILO, IHI, -1 ) )
      NBMIN = 2
      IWS = 1
      IF( NB.GT.1 .AND. NB.LT.NH ) THEN
!
!        Determine when to cross over from blocked to unblocked code
!        (last block is always handled by unblocked code).
!
         NX = MAX( NB, ILAENV( 3, 'SGEHRD', ' ', N, ILO, IHI, -1 ) )
         IF( NX.LT.NH ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            IWS = N*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  determine the
!              minimum value of NB, and reduce NB or force use of
!              unblocked code.
!
               NBMIN = MAX( 2, ILAENV( 2, 'SGEHRD', ' ', N, ILO, IHI, &
                       -1 ) )
               IF( LWORK.GE.N*NBMIN ) THEN
                  NB = LWORK / N
               ELSE
                  NB = 1
               END IF
            END IF
         END IF
      END IF
      LDWORK = N
!
      IF( NB.LT.NBMIN .OR. NB.GE.NH ) THEN
!
!        Use unblocked code below
!
         I = ILO
!
      ELSE
!
!        Use blocked code
!
         DO 30 I = ILO, IHI - 1 - NX, NB
            IB = MIN( NB, IHI-I )
!
!           Reduce columns i:i+ib-1 to Hessenberg form, returning the
!           matrices V and T of the block reflector H = I - V*T*V'
!           which performs the reduction, and also the matrix Y = A*V*T
!
            CALL SLAHRD( IHI, I, IB, A( 1, I ), LDA, TAU( I ), T, LDT, &
                         WORK, LDWORK )
!
!           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
!           right, computing  A := A - Y * V'. V(i+ib,ib-1) must be set
!           to 1.
!
            EI = A( I+IB, I+IB-1 )
            A( I+IB, I+IB-1 ) = ONE
            CALL SGEMM( 'No transpose', 'Transpose', IHI, IHI-I-IB+1, &
                        IB, -ONE, WORK, LDWORK, A( I+IB, I ), LDA, ONE, &
                        A( 1, I+IB ), LDA )
            A( I+IB, I+IB-1 ) = EI
!
!           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
!           left
!
            CALL SLARFB( 'Left', 'Transpose', 'Forward', 'Columnwise', &
                         IHI-I, N-I-IB+1, IB, A( I+1, I ), LDA, T, LDT, &
                         A( I+1, I+IB ), LDA, WORK, LDWORK )
   30    CONTINUE
      END IF
!
!     Use unblocked code to reduce the rest of the matrix
!
      CALL SGEHD2( N, I, IHI, A, LDA, TAU, WORK, IINFO )
      WORK( 1 ) = IWS
!
      RETURN
!
!     End of SGEHRD
!
      END SUBROUTINE SGEHRD

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, Z, &
                         LDZ, WORK, LWORK, INFO )

!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          COMPZ, JOB
      INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N
!     ..
!     .. Array Arguments ..
      REAL(r8)               H( LDH, * ), WI( * ), WORK( * ), WR( * ), &
                         Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  SHSEQR computes the eigenvalues of a real upper Hessenberg matrix H
!  and, optionally, the matrices T and Z from the Schur decomposition
!  H = Z T Z**T, where T is an upper quasi-triangular matrix (the Schur
!  form), and Z is the orthogonal matrix of Schur vectors.
!
!  Optionally Z may be postmultiplied into an input orthogonal matrix Q,
!  so that this routine can give the Schur factorization of a matrix A
!  which has been reduced to the Hessenberg form H by the orthogonal
!  matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.
!
!  Arguments
!  =========
!
!  JOB     (input) CHARACTER*1
!          = 'E':  compute eigenvalues only;
!          = 'S':  compute eigenvalues and the Schur form T.
!
!  COMPZ   (input) CHARACTER*1
!          = 'N':  no Schur vectors are computed;
!          = 'I':  Z is initialized to the unit matrix and the matrix Z
!                  of Schur vectors of H is returned;
!          = 'V':  Z must contain an orthogonal matrix Q on entry, and
!                  the product Q*Z is returned.
!
!  N       (input) INTEGER
!          The order of the matrix H.  N >= 0.
!
!  ILO     (input) INTEGER
!  IHI     (input) INTEGER
!          It is assumed that H is already upper triangular in rows
!          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
!          set by a previous call to SGEBAL, and then passed to SGEHRD
!          when the matrix output by SGEBAL is reduced to Hessenberg
!          form. Otherwise ILO and IHI should be set to 1 and N
!          respectively.
!          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!
!  H       (input/output) REAL array, dimension (LDH,N)
!          On entry, the upper Hessenberg matrix H.
!          On exit, if JOB = 'S', H contains the upper quasi-triangular
!          matrix T from the Schur decomposition (the Schur form);
!          2-by-2 diagonal blocks (corresponding to complex conjugate
!          pairs of eigenvalues) are returned in standard form, with
!          H(i,i) = H(i+1,i+1) and H(i+1,i)*H(i,i+1) < 0. If JOB = 'E',
!          the contents of H are unspecified on exit.
!
!  LDH     (input) INTEGER
!          The leading dimension of the array H. LDH >= max(1,N).
!
!  WR      (output) REAL array, dimension (N)
!  WI      (output) REAL array, dimension (N)
!          The real and imaginary parts, respectively, of the computed
!          eigenvalues. If two eigenvalues are computed as a complex
!          conjugate pair, they are stored in consecutive elements of
!          WR and WI, say the i-th and (i+1)th, with WI(i) > 0 and
!          WI(i+1) < 0. If JOB = 'S', the eigenvalues are stored in the
!          same order as on the diagonal of the Schur form returned in
!          H, with WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2
!          diagonal block, WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and
!          WI(i+1) = -WI(i).
!
!  Z       (input/output) REAL array, dimension (LDZ,N)
!          If COMPZ = 'N': Z is not referenced.
!          If COMPZ = 'I': on entry, Z need not be set, and on exit, Z
!          contains the orthogonal matrix Z of the Schur vectors of H.
!          If COMPZ = 'V': on entry Z must contain an N-by-N matrix Q,
!          which is assumed to be equal to the unit matrix except for
!          the submatrix Z(ILO:IHI,ILO:IHI); on exit Z contains Q*Z.
!          Normally Q is the orthogonal matrix generated by SORGHR after
!          the call to SGEHRD which formed the Hessenberg matrix H.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.
!          LDZ >= max(1,N) if COMPZ = 'I' or 'V'; LDZ >= 1 otherwise.
!
!  WORK    (workspace) REAL array, dimension (N)
!
!  LWORK   (input) INTEGER
!          This argument is currently redundant.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, SHSEQR failed to compute all of the
!                eigenvalues in a total of 30*(IHI-ILO+1) iterations;
!                elements 1:ilo-1 and i+1:n of WR and WI contain those
!                eigenvalues which have been successfully computed.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0E+0_r8, ONE = 1.0E+0_r8, TWO = 2.0E+0_r8 )
      REAL(r8)               CONST
      PARAMETER          ( CONST = 1.5E+0_r8 )
      INTEGER            NSMAX, LDS
      PARAMETER          ( NSMAX = 15, LDS = NSMAX )
!     ..
!     .. Local Scalars ..
      LOGICAL            INITZ, WANTT, WANTZ
      INTEGER            I, I1, I2, IERR, II, ITEMP, ITN, ITS, J, K, L, &
                         MAXB, NH, NR, NS, NV
      REAL(r8)               ABSW, OVFL, SMLNUM, TAU, TEMP, TST1, ULP, UNFL
!     ..
!     .. Local Arrays ..
      REAL(r8)               S( LDS, NSMAX ), V( NSMAX+1 ), VV( NSMAX+1 )
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters
!
      WANTT = LSAME( JOB, 'S' )
      INITZ = LSAME( COMPZ, 'I' )
      WANTZ = INITZ .OR. LSAME( COMPZ, 'V' )
!
      INFO = 0
      IF( .NOT.LSAME( JOB, 'E' ) .AND. .NOT.WANTT ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( COMPZ, 'N' ) .AND. .NOT.WANTZ ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -5
      ELSE IF( LDH.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDZ.LT.1 .OR. WANTZ .AND. LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SHSEQR', -INFO )
         RETURN
      END IF
!
!     Initialize Z, if necessary
!
      IF( INITZ ) &
         CALL SLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
!
!     Store the eigenvalues isolated by SGEBAL.
!
      DO 10 I = 1, ILO - 1
         WR( I ) = H( I, I )
         WI( I ) = ZERO
   10 CONTINUE
      DO 20 I = IHI + 1, N
         WR( I ) = H( I, I )
         WI( I ) = ZERO
   20 CONTINUE
!
!     Quick return if possible.
!
      IF( N.EQ.0 ) &
         RETURN
      IF( ILO.EQ.IHI ) THEN
         WR( ILO ) = H( ILO, ILO )
         WI( ILO ) = ZERO
         RETURN
      END IF
!
!     Set rows and columns ILO to IHI to zero below the first
!     subdiagonal.
!
      DO 40 J = ILO, IHI - 2
         DO 30 I = J + 2, N
            H( I, J ) = ZERO
   30    CONTINUE
   40 CONTINUE
      NH = IHI - ILO + 1
!
!     Determine the order of the multi-shift QR algorithm to be used.
!
      NS = ILAENV( 4, 'SHSEQR', JOB // COMPZ, N, ILO, IHI, -1 )
      MAXB = ILAENV( 8, 'SHSEQR', JOB // COMPZ, N, ILO, IHI, -1 )
      IF( NS.LE.2 .OR. NS.GT.NH .OR. MAXB.GE.NH ) THEN
!
!        Use the standard double-shift algorithm
!
         CALL SLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, &
                      IHI, Z, LDZ, INFO )
         RETURN
      END IF
      MAXB = MAX( 3, MAXB )
      NS = MIN( NS, MAXB, NSMAX )
!
!     Now 2 < NS <= MAXB < NH.
!
!     Set machine-dependent constants for the stopping criterion.
!     If norm(H) <= sqrt(OVFL), overflow should not occur.
!
      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      CALL SLABAD( UNFL, OVFL )
      ULP = SLAMCH( 'Precision' )
      SMLNUM = UNFL*( NH / ULP )
!
!     I1 and I2 are the indices of the first row and last column of H
!     to which transformations must be applied. If eigenvalues only are
!     being computed, I1 and I2 are set inside the main loop.
!
      IF( WANTT ) THEN
         I1 = 1
         I2 = N
      END IF
!
!     ITN is the total number of multiple-shift QR iterations allowed.
!
      ITN = 30*NH
!
!     The main loop begins here. I is the loop index and decreases from
!     IHI to ILO in steps of at most MAXB. Each iteration of the loop
!     works with the active submatrix in rows and columns L to I.
!     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
!     H(L,L-1) is negligible so that the matrix splits.
!
      I = IHI
   50 CONTINUE
      L = ILO
      IF( I.LT.ILO ) &
         GO TO 170
!
!     Perform multiple-shift QR iterations on rows and columns ILO to I
!     until a submatrix of order at most MAXB splits off at the bottom
!     because a subdiagonal element has become negligible.
!
      DO 150 ITS = 0, ITN
!
!        Look for a single small subdiagonal element.
!
         DO 60 K = I, L + 1, -1
            TST1 = ABS( H( K-1, K-1 ) ) + ABS( H( K, K ) )
            IF( TST1.EQ.ZERO ) &
               TST1 = SLANHS( '1', I-L+1, H( L, L ), LDH, WORK )
            IF( ABS( H( K, K-1 ) ).LE.MAX( ULP*TST1, SMLNUM ) ) &
               GO TO 70
   60    CONTINUE
   70    CONTINUE
         L = K
         IF( L.GT.ILO ) THEN
!
!           H(L,L-1) is negligible.
!
            H( L, L-1 ) = ZERO
         END IF
!
!        Exit from loop if a submatrix of order <= MAXB has split off.
!
         IF( L.GE.I-MAXB+1 ) &
            GO TO 160
!
!        Now the active submatrix is in rows and columns L to I. If
!        eigenvalues only are being computed, only the active submatrix
!        need be transformed.
!
         IF( .NOT.WANTT ) THEN
            I1 = L
            I2 = I
         END IF
!
         IF( ITS.EQ.20 .OR. ITS.EQ.30 ) THEN
!
!           Exceptional shifts.
!
            DO 80 II = I - NS + 1, I
               WR( II ) = CONST*( ABS( H( II, II-1 ) )+ &
                          ABS( H( II, II ) ) )
               WI( II ) = ZERO
   80       CONTINUE
         ELSE
!
!           Use eigenvalues of trailing submatrix of order NS as shifts.
!
            CALL SLACPY( 'Full', NS, NS, H( I-NS+1, I-NS+1 ), LDH, S, &
                         LDS )
            CALL SLAHQR( .FALSE., .FALSE., NS, 1, NS, S, LDS, &
                         WR( I-NS+1 ), WI( I-NS+1 ), 1, NS, Z, LDZ, &
                         IERR )
            IF( IERR.GT.0 ) THEN
!
!              If SLAHQR failed to compute all NS eigenvalues, use the
!              unconverged diagonal elements as the remaining shifts.
!
               DO 90 II = 1, IERR
                  WR( I-NS+II ) = S( II, II )
                  WI( I-NS+II ) = ZERO
   90          CONTINUE
            END IF
         END IF
!
!        Form the first column of (G-w(1)) (G-w(2)) . . . (G-w(ns))
!        where G is the Hessenberg submatrix H(L:I,L:I) and w is
!        the vector of shifts (stored in WR and WI). The result is
!        stored in the local array V.
!
         V( 1 ) = ONE
         DO 100 II = 2, NS + 1
            V( II ) = ZERO
  100    CONTINUE
         NV = 1
         DO 120 J = I - NS + 1, I
            IF( WI( J ).GE.ZERO ) THEN
               IF( WI( J ).EQ.ZERO ) THEN
!
!                 real shift
!
                  CALL SCOPY( NV+1, V, 1, VV, 1 )
                  CALL SGEMV( 'No transpose', NV+1, NV, ONE, H( L, L ), &
                              LDH, VV, 1, -WR( J ), V, 1 )
                  NV = NV + 1
               ELSE IF( WI( J ).GT.ZERO ) THEN
!
!                 complex conjugate pair of shifts
!
                  CALL SCOPY( NV+1, V, 1, VV, 1 )
                  CALL SGEMV( 'No transpose', NV+1, NV, ONE, H( L, L ), &
                              LDH, V, 1, -TWO*WR( J ), VV, 1 )
                  ITEMP = ISAMAX( NV+1, VV, 1 )
                  TEMP = ONE / MAX( ABS( VV( ITEMP ) ), SMLNUM )
                  CALL SSCAL( NV+1, TEMP, VV, 1 )
                  ABSW = SLAPY2( WR( J ), WI( J ) )
                  TEMP = ( TEMP*ABSW )*ABSW
                  CALL SGEMV( 'No transpose', NV+2, NV+1, ONE, &
                              H( L, L ), LDH, VV, 1, TEMP, V, 1 )
                  NV = NV + 2
               END IF
!
!              Scale V(1:NV) so that max(abs(V(i))) = 1. If V is zero,
!              reset it to the unit vector.
!
               ITEMP = ISAMAX( NV, V, 1 )
               TEMP = ABS( V( ITEMP ) )
               IF( TEMP.EQ.ZERO ) THEN
                  V( 1 ) = ONE
                  DO 110 II = 2, NV
                     V( II ) = ZERO
  110             CONTINUE
               ELSE
                  TEMP = MAX( TEMP, SMLNUM )
                  CALL SSCAL( NV, ONE / TEMP, V, 1 )
               END IF
            END IF
  120    CONTINUE
!
!        Multiple-shift QR step
!
         DO 140 K = L, I - 1
!
!           The first iteration of this loop determines a reflection G
!           from the vector V and applies it from left and right to H,
!           thus creating a nonzero bulge below the subdiagonal.
!
!           Each subsequent iteration determines a reflection G to
!           restore the Hessenberg form in the (K-1)th column, and thus
!           chases the bulge one step toward the bottom of the active
!           submatrix. NR is the order of G.
!
            NR = MIN( NS+1, I-K+1 )
            IF( K.GT.L ) &
               CALL SCOPY( NR, H( K, K-1 ), 1, V, 1 )
            CALL SLARFG( NR, V( 1 ), V( 2 ), 1, TAU )
            IF( K.GT.L ) THEN
               H( K, K-1 ) = V( 1 )
               DO 130 II = K + 1, I
                  H( II, K-1 ) = ZERO
  130          CONTINUE
            END IF
            V( 1 ) = ONE
!
!           Apply G from the left to transform the rows of the matrix in
!           columns K to I2.
!
            CALL SLARFX( 'Left', NR, I2-K+1, V, TAU, H( K, K ), LDH, &
                         WORK )
!
!           Apply G from the right to transform the columns of the
!           matrix in rows I1 to min(K+NR,I).
!
            CALL SLARFX( 'Right', MIN( K+NR, I )-I1+1, NR, V, TAU, &
                         H( I1, K ), LDH, WORK )
!
            IF( WANTZ ) THEN
!
!              Accumulate transformations in the matrix Z
!
               CALL SLARFX( 'Right', NH, NR, V, TAU, Z( ILO, K ), LDZ, &
                            WORK )
            END IF
  140    CONTINUE
!
  150 CONTINUE
!
!     Failure to converge in remaining number of iterations
!
      INFO = I
      RETURN
!
  160 CONTINUE
!
!     A submatrix of order <= MAXB in rows and columns L to I has split
!     off. Use the double-shift QR algorithm to handle it.
!
      CALL SLAHQR( WANTT, WANTZ, N, L, I, H, LDH, WR, WI, ILO, IHI, Z, &
                   LDZ, INFO )
      IF( INFO.GT.0 ) &
         RETURN
!
!     Decrement number of remaining iterations, and return to start of
!     the main loop with a new value of I.
!
      ITN = ITN - ITS
      I = L - 1
      GO TO 50
!
  170 CONTINUE
      RETURN
!
!     End of SHSEQR
!
      END SUBROUTINE SHSEQR

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SLABAD( SMALL, LARGE )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      REAL(r8)               LARGE, SMALL
!     ..
!
!  Purpose
!  =======
!
!  SLABAD takes as input the values computed by SLAMCH for underflow and
!  overflow, and returns the square root of each of these values if the
!  log of LARGE is sufficiently large.  This subroutine is intended to
!  identify machines with a large exponent range, such as the Crays, and
!  redefine the underflow and overflow limits to be the square roots of
!  the values computed by SLAMCH.  This subroutine is needed because
!  SLAMCH does not compensate for poor arithmetic in the upper half of
!  the exponent range, as is found on a Cray.
!
!  Arguments
!  =========
!
!  SMALL   (input/output) REAL
!          On entry, the underflow threshold as computed by SLAMCH.
!          On exit, if LOG10(LARGE) is sufficiently large, the square
!          root of SMALL, otherwise unchanged.
!
!  LARGE   (input/output) REAL
!          On entry, the overflow threshold as computed by SLAMCH.
!          On exit, if LOG10(LARGE) is sufficiently large, the square
!          root of LARGE, otherwise unchanged.
!
!  =====================================================================
!
!     ..
!     .. Executable Statements ..
!
!     If it looks like we're on a Cray, take the square root of
!     SMALL and LARGE to avoid overflow and underflow problems.
!
      IF( LOG10( LARGE ).GT.2000._r8 ) THEN
         SMALL = SQRT( SMALL )
         LARGE = SQRT( LARGE )
      END IF
!
      RETURN
!
!     End of SLABAD
!
      END SUBROUTINE SLABAD

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SLACPY( UPLO, M, N, A, LDA, B, LDB )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
!     ..
!     .. Array Arguments ..
      REAL(r8)               A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  SLACPY copies all or part of a two-dimensional matrix A to another
!  matrix B.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies the part of the matrix A to be copied to B.
!          = 'U':      Upper triangular part
!          = 'L':      Lower triangular part
!          Otherwise:  All of the matrix A
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input) REAL array, dimension (LDA,N)
!          The m by n matrix A.  If UPLO = 'U', only the upper triangle
!          or trapezoid is accessed; if UPLO = 'L', only the lower
!          triangle or trapezoid is accessed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  B       (output) REAL array, dimension (LDB,N)
!          On exit, B = A in the locations specified by UPLO.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,M).
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
      RETURN
!
!     End of SLACPY
!
      END SUBROUTINE SLACPY

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SLADIV( A, B, C, D, P, Q )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      REAL(r8)               A, B, C, D, P, Q
!     ..
!
!  Purpose
!  =======
!
!  SLADIV performs complex division in  real arithmetic
!
!                        a + i*b
!             p + i*q = ---------
!                        c + i*d
!
!  The algorithm is due to Robert L. Smith and can be found
!  in D. Knuth, The art of Computer Programming, Vol.2, p.195
!
!  Arguments
!  =========
!
!  A       (input) REAL
!  B       (input) REAL
!  C       (input) REAL
!  D       (input) REAL
!          The scalars a, b, c, and d in the above expression.
!
!  P       (output) REAL
!  Q       (output) REAL
!          The scalars p and q in the above expression.
!
!  =====================================================================
!
!     .. Local Scalars ..
      REAL(r8)               E, F
!     ..
!     .. Executable Statements ..
!
      IF( ABS( D ).LT.ABS( C ) ) THEN
         E = D / C
         F = C + D*E
         P = ( A+B*E ) / F
         Q = ( B-A*E ) / F
      ELSE
         E = C / D
         F = D + C*E
         P = ( B+A*E ) / F
         Q = ( -A+B*E ) / F
      END IF
!
      RETURN
!
!     End of SLADIV
!
      END SUBROUTINE SLADIV

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, &
                         ILOZ, IHIZ, Z, LDZ, INFO )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            WANTT, WANTZ
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
!     ..
!     .. Array Arguments ..
      REAL(r8)               H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  SLAHQR is an auxiliary routine called by SHSEQR to update the
!  eigenvalues and Schur decomposition already computed by SHSEQR, by
!  dealing with the Hessenberg submatrix in rows and columns ILO to IHI.
!
!  Arguments
!  =========
!
!  WANTT   (input) LOGICAL
!          = .TRUE. : the full Schur form T is required;
!          = .FALSE.: only eigenvalues are required.
!
!  WANTZ   (input) LOGICAL
!          = .TRUE. : the matrix of Schur vectors Z is required;
!          = .FALSE.: Schur vectors are not required.
!
!  N       (input) INTEGER
!          The order of the matrix H.  N >= 0.
!
!  ILO     (input) INTEGER
!  IHI     (input) INTEGER
!          It is assumed that H is already upper quasi-triangular in
!          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless
!          ILO = 1). SLAHQR works primarily with the Hessenberg
!          submatrix in rows and columns ILO to IHI, but applies
!          transformations to all of H if WANTT is .TRUE..
!          1 <= ILO <= max(1,IHI); IHI <= N.
!
!  H       (input/output) REAL array, dimension (LDH,N)
!          On entry, the upper Hessenberg matrix H.
!          On exit, if WANTT is .TRUE., H is upper quasi-triangular in
!          rows and columns ILO:IHI, with any 2-by-2 diagonal blocks in
!          standard form. If WANTT is .FALSE., the contents of H are
!          unspecified on exit.
!
!  LDH     (input) INTEGER
!          The leading dimension of the array H. LDH >= max(1,N).
!
!  WR      (output) REAL array, dimension (N)
!  WI      (output) REAL array, dimension (N)
!          The real and imaginary parts, respectively, of the computed
!          eigenvalues ILO to IHI are stored in the corresponding
!          elements of WR and WI. If two eigenvalues are computed as a
!          complex conjugate pair, they are stored in consecutive
!          elements of WR and WI, say the i-th and (i+1)th, with
!          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
!          eigenvalues are stored in the same order as on the diagonal
!          of the Schur form returned in H, with WR(i) = H(i,i), and, if
!          H(i:i+1,i:i+1) is a 2-by-2 diagonal block,
!          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).
!
!  ILOZ    (input) INTEGER
!  IHIZ    (input) INTEGER
!          Specify the rows of Z to which transformations must be
!          applied if WANTZ is .TRUE..
!          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
!
!  Z       (input/output) REAL array, dimension (LDZ,N)
!          If WANTZ is .TRUE., on entry Z must contain the current
!          matrix Z of transformations accumulated by SHSEQR, and on
!          exit Z has been updated; transformations are applied only to
!          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
!          If WANTZ is .FALSE., Z is not referenced.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z. LDZ >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          > 0: SLAHQR failed to compute all the eigenvalues ILO to IHI
!               in a total of 30*(IHI-ILO+1) iterations; if INFO = i,
!               elements i+1:ihi of WR and WI contain those eigenvalues
!               which have been successfully computed.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0_r8, ONE = 1.0E+0_r8 )
      REAL(r8)               DAT1, DAT2
      PARAMETER          ( DAT1 = 0.75E+0_r8, DAT2 = -0.4375E+0_r8 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I1, I2, ITN, ITS, J, K, L, M, NH, NR, NZ
      REAL(r8)               CS, H00, H10, H11, H12, H21, H22, H33, H33S, &
                         H43H34, H44, H44S, OVFL, S, SMLNUM, SN, SUM, &
                         T1, T2, T3, TST1, ULP, UNFL, V1, V2, V3
!     ..
!     .. Local Arrays ..
      REAL(r8)               V( 3 ), WORK( 1 )
!     ..
!     .. Executable Statements ..
!
      INFO = 0
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
      IF( ILO.EQ.IHI ) THEN
         WR( ILO ) = H( ILO, ILO )
         WI( ILO ) = ZERO
         RETURN
      END IF
!
      NH = IHI - ILO + 1
      NZ = IHIZ - ILOZ + 1
!
!     Set machine-dependent constants for the stopping criterion.
!     If norm(H) <= sqrt(OVFL), overflow should not occur.
!
      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      CALL SLABAD( UNFL, OVFL )
      ULP = SLAMCH( 'Precision' )
      SMLNUM = UNFL*( NH / ULP )
!
!     I1 and I2 are the indices of the first row and last column of H
!     to which transformations must be applied. If eigenvalues only are
!     being computed, I1 and I2 are set inside the main loop.
!
      IF( WANTT ) THEN
         I1 = 1
         I2 = N
      END IF
!
!     ITN is the total number of QR iterations allowed.
!
      ITN = 30*NH
!
!     The main loop begins here. I is the loop index and decreases from
!     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
!     with the active submatrix in rows and columns L to I.
!     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
!     H(L,L-1) is negligible so that the matrix splits.
!
      I = IHI
   10 CONTINUE
      L = ILO
      IF( I.LT.ILO ) &
         GO TO 150
!
!     Perform QR iterations on rows and columns ILO to I until a
!     submatrix of order 1 or 2 splits off at the bottom because a
!     subdiagonal element has become negligible.
!
      DO 130 ITS = 0, ITN
!
!        Look for a single small subdiagonal element.
!
         DO 20 K = I, L + 1, -1
            TST1 = ABS( H( K-1, K-1 ) ) + ABS( H( K, K ) )
            IF( TST1.EQ.ZERO ) &
               TST1 = SLANHS( '1', I-L+1, H( L, L ), LDH, WORK )
            IF( ABS( H( K, K-1 ) ).LE.MAX( ULP*TST1, SMLNUM ) ) &
               GO TO 30
   20    CONTINUE
   30    CONTINUE
         L = K
         IF( L.GT.ILO ) THEN
!
!           H(L,L-1) is negligible
!
            H( L, L-1 ) = ZERO
         END IF
!
!        Exit from loop if a submatrix of order 1 or 2 has split off.
!
         IF( L.GE.I-1 ) &
            GO TO 140
!
!        Now the active submatrix is in rows and columns L to I. If
!        eigenvalues only are being computed, only the active submatrix
!        need be transformed.
!
         IF( .NOT.WANTT ) THEN
            I1 = L
            I2 = I
         END IF
!
         IF( ITS.EQ.10 .OR. ITS.EQ.20 ) THEN
!
!           Exceptional shift.
!
            S = ABS( H( I, I-1 ) ) + ABS( H( I-1, I-2 ) )
            H44 = DAT1*S
            H33 = H44
            H43H34 = DAT2*S*S
         ELSE
!
!           Prepare to use Wilkinson's double shift
!
            H44 = H( I, I )
            H33 = H( I-1, I-1 )
            H43H34 = H( I, I-1 )*H( I-1, I )
         END IF
!
!        Look for two consecutive small subdiagonal elements.
!
         DO 40 M = I - 2, L, -1
!
!           Determine the effect of starting the double-shift QR
!           iteration at row M, and see if this would make H(M,M-1)
!           negligible.
!
            H11 = H( M, M )
            H22 = H( M+1, M+1 )
            H21 = H( M+1, M )
            H12 = H( M, M+1 )
            H44S = H44 - H11
            H33S = H33 - H11
            V1 = ( H33S*H44S-H43H34 ) / H21 + H12
            V2 = H22 - H11 - H33S - H44S
            V3 = H( M+2, M+1 )
            S = ABS( V1 ) + ABS( V2 ) + ABS( V3 )
            V1 = V1 / S
            V2 = V2 / S
            V3 = V3 / S
            V( 1 ) = V1
            V( 2 ) = V2
            V( 3 ) = V3
            IF( M.EQ.L ) &
               GO TO 50
            H00 = H( M-1, M-1 )
            H10 = H( M, M-1 )
            TST1 = ABS( V1 )*( ABS( H00 )+ABS( H11 )+ABS( H22 ) )
            IF( ABS( H10 )*( ABS( V2 )+ABS( V3 ) ).LE.ULP*TST1 ) &
               GO TO 50
   40    CONTINUE
   50    CONTINUE
!
!        Double-shift QR step
!
         DO 120 K = M, I - 1
!
!           The first iteration of this loop determines a reflection G
!           from the vector V and applies it from left and right to H,
!           thus creating a nonzero bulge below the subdiagonal.
!
!           Each subsequent iteration determines a reflection G to
!           restore the Hessenberg form in the (K-1)th column, and thus
!           chases the bulge one step toward the bottom of the active
!           submatrix. NR is the order of G.
!
            NR = MIN( 3, I-K+1 )
            IF( K.GT.M ) &
               CALL SCOPY( NR, H( K, K-1 ), 1, V, 1 )
            CALL SLARFG( NR, V( 1 ), V( 2 ), 1, T1 )
            IF( K.GT.M ) THEN
               H( K, K-1 ) = V( 1 )
               H( K+1, K-1 ) = ZERO
               IF( K.LT.I-1 ) &
                  H( K+2, K-1 ) = ZERO
            ELSE IF( M.GT.L ) THEN
               H( K, K-1 ) = -H( K, K-1 )
            END IF
            V2 = V( 2 )
            T2 = T1*V2
            IF( NR.EQ.3 ) THEN
               V3 = V( 3 )
               T3 = T1*V3
!
!              Apply G from the left to transform the rows of the matrix
!              in columns K to I2.
!
               DO 60 J = K, I2
                  SUM = H( K, J ) + V2*H( K+1, J ) + V3*H( K+2, J )
                  H( K, J ) = H( K, J ) - SUM*T1
                  H( K+1, J ) = H( K+1, J ) - SUM*T2
                  H( K+2, J ) = H( K+2, J ) - SUM*T3
   60          CONTINUE
!
!              Apply G from the right to transform the columns of the
!              matrix in rows I1 to min(K+3,I).
!
               DO 70 J = I1, MIN( K+3, I )
                  SUM = H( J, K ) + V2*H( J, K+1 ) + V3*H( J, K+2 )
                  H( J, K ) = H( J, K ) - SUM*T1
                  H( J, K+1 ) = H( J, K+1 ) - SUM*T2
                  H( J, K+2 ) = H( J, K+2 ) - SUM*T3
   70          CONTINUE
!
               IF( WANTZ ) THEN
!
!                 Accumulate transformations in the matrix Z
!
                  DO 80 J = ILOZ, IHIZ
                     SUM = Z( J, K ) + V2*Z( J, K+1 ) + V3*Z( J, K+2 )
                     Z( J, K ) = Z( J, K ) - SUM*T1
                     Z( J, K+1 ) = Z( J, K+1 ) - SUM*T2
                     Z( J, K+2 ) = Z( J, K+2 ) - SUM*T3
   80             CONTINUE
               END IF
            ELSE IF( NR.EQ.2 ) THEN
!
!              Apply G from the left to transform the rows of the matrix
!              in columns K to I2.
!
               DO 90 J = K, I2
                  SUM = H( K, J ) + V2*H( K+1, J )
                  H( K, J ) = H( K, J ) - SUM*T1
                  H( K+1, J ) = H( K+1, J ) - SUM*T2
   90          CONTINUE
!
!              Apply G from the right to transform the columns of the
!              matrix in rows I1 to min(K+3,I).
!
               DO 100 J = I1, I
                  SUM = H( J, K ) + V2*H( J, K+1 )
                  H( J, K ) = H( J, K ) - SUM*T1
                  H( J, K+1 ) = H( J, K+1 ) - SUM*T2
  100          CONTINUE
!
               IF( WANTZ ) THEN
!
!                 Accumulate transformations in the matrix Z
!
                  DO 110 J = ILOZ, IHIZ
                     SUM = Z( J, K ) + V2*Z( J, K+1 )
                     Z( J, K ) = Z( J, K ) - SUM*T1
                     Z( J, K+1 ) = Z( J, K+1 ) - SUM*T2
  110             CONTINUE
               END IF
            END IF
  120    CONTINUE
!
  130 CONTINUE
!
!     Failure to converge in remaining number of iterations
!
      INFO = I
      RETURN
!
  140 CONTINUE
!
      IF( L.EQ.I ) THEN
!
!        H(I,I-1) is negligible: one eigenvalue has converged.
!
         WR( I ) = H( I, I )
         WI( I ) = ZERO
      ELSE IF( L.EQ.I-1 ) THEN
!
!        H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
!
!        Transform the 2-by-2 submatrix to standard Schur form,
!        and compute and store the eigenvalues.
!
         CALL SLANV2( H( I-1, I-1 ), H( I-1, I ), H( I, I-1 ), &
                      H( I, I ), WR( I-1 ), WI( I-1 ), WR( I ), WI( I ), &
                      CS, SN )
!
         IF( WANTT ) THEN
!
!           Apply the transformation to the rest of H.
!
            IF( I2.GT.I ) &
               CALL SROT( I2-I, H( I-1, I+1 ), LDH, H( I, I+1 ), LDH, &
                          CS, SN )
            CALL SROT( I-I1-1, H( I1, I-1 ), 1, H( I1, I ), 1, CS, SN )
         END IF
         IF( WANTZ ) THEN
!
!           Apply the transformation to Z.
!
            CALL SROT( NZ, Z( ILOZ, I-1 ), 1, Z( ILOZ, I ), 1, CS, SN )
         END IF
      END IF
!
!     Decrement number of remaining iterations, and return to start of
!     the main loop with new value of I.
!
      ITN = ITN - ITS
      I = L - 1
      GO TO 10
!
  150 CONTINUE
      RETURN
!
!     End of SLAHQR
!
      END SUBROUTINE SLAHQR

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SLAHRD( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      INTEGER            K, LDA, LDT, LDY, N, NB
!     ..
!     .. Array Arguments ..
      REAL(r8)               A( LDA, * ), T( LDT, NB ), TAU( NB ), &
                         Y( LDY, NB )
!     ..
!
!  Purpose
!  =======
!
!  SLAHRD reduces the first NB columns of a real general n-by-(n-k+1)
!  matrix A so that elements below the k-th subdiagonal are zero. The
!  reduction is performed by an orthogonal similarity transformation
!  Q' * A * Q. The routine returns the matrices V and T which determine
!  Q as a block reflector I - V*T*V', and also the matrix Y = A * V * T.
!
!  This is an auxiliary routine called by SGEHRD.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.
!
!  K       (input) INTEGER
!          The offset for the reduction. Elements below the k-th
!          subdiagonal in the first NB columns are reduced to zero.
!
!  NB      (input) INTEGER
!          The number of columns to be reduced.
!
!  A       (input/output) REAL array, dimension (LDA,N-K+1)
!          On entry, the n-by-(n-k+1) general matrix A.
!          On exit, the elements on and above the k-th subdiagonal in
!          the first NB columns are overwritten with the corresponding
!          elements of the reduced matrix; the elements below the k-th
!          subdiagonal, with the array TAU, represent the matrix Q as a
!          product of elementary reflectors. The other columns of A are
!          unchanged. See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  TAU     (output) REAL array, dimension (NB)
!          The scalar factors of the elementary reflectors. See Further
!          Details.
!
!  T       (output) REAL array, dimension (NB,NB)
!          The upper triangular matrix T.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T.  LDT >= NB.
!
!  Y       (output) REAL array, dimension (LDY,NB)
!          The n-by-nb matrix Y.
!
!  LDY     (input) INTEGER
!          The leading dimension of the array Y. LDY >= N.
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of nb elementary reflectors
!
!     Q = H(1) H(2) . . . H(nb).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in
!  A(i+k+1:n,i), and tau in TAU(i).
!
!  The elements of the vectors v together form the (n-k+1)-by-nb matrix
!  V which is needed, with T and Y, to apply the transformation to the
!  unreduced part of the matrix, using an update of the form:
!  A := (I - V*T*V') * (A - Y*V').
!
!  The contents of A on exit are illustrated by the following example
!  with n = 7, k = 3 and nb = 2:
!
!     ( a   h   a   a   a )
!     ( a   h   a   a   a )
!     ( a   h   a   a   a )
!     ( h   h   a   a   a )
!     ( v1  h   a   a   a )
!     ( v1  v2  a   a   a )
!     ( v1  v2  a   a   a )
!
!  where a denotes an element of the original matrix A, h denotes a
!  modified element of the upper Hessenberg matrix H, and vi denotes an
!  element of the vector defining H(i).
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0_r8, ONE = 1.0E+0_r8 )
!     ..
!     .. Local Scalars ..
      INTEGER            I
      REAL(r8)               EI
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.LE.1 ) &
         RETURN
!
      DO 10 I = 1, NB
         IF( I.GT.1 ) THEN
!
!           Update A(1:n,i)
!
!           Compute i-th column of A - Y * V'
!
            CALL SGEMV( 'No transpose', N, I-1, -ONE, Y, LDY, &
                        A( K+I-1, 1 ), LDA, ONE, A( 1, I ), 1 )
!
!           Apply I - V * T' * V' to this column (call it b) from the
!           left, using the last column of T as workspace
!
!           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
!                    ( V2 )             ( b2 )
!
!           where V1 is unit lower triangular
!
!           w := V1' * b1
!
            CALL SCOPY( I-1, A( K+1, I ), 1, T( 1, NB ), 1 )
            CALL STRMV( 'Lower', 'Transpose', 'Unit', I-1, A( K+1, 1 ), &
                        LDA, T( 1, NB ), 1 )
!
!           w := w + V2'*b2
!
            CALL SGEMV( 'Transpose', N-K-I+1, I-1, ONE, A( K+I, 1 ), &
                        LDA, A( K+I, I ), 1, ONE, T( 1, NB ), 1 )
!
!           w := T'*w
!
            CALL STRMV( 'Upper', 'Transpose', 'Non-unit', I-1, T, LDT, &
                        T( 1, NB ), 1 )
!
!           b2 := b2 - V2*w
!
            CALL SGEMV( 'No transpose', N-K-I+1, I-1, -ONE, A( K+I, 1 ), &
                        LDA, T( 1, NB ), 1, ONE, A( K+I, I ), 1 )
!
!           b1 := b1 - V1*w
!
            CALL STRMV( 'Lower', 'No transpose', 'Unit', I-1, &
                        A( K+1, 1 ), LDA, T( 1, NB ), 1 )
            CALL SAXPY( I-1, -ONE, T( 1, NB ), 1, A( K+1, I ), 1 )
!
            A( K+I-1, I-1 ) = EI
         END IF
!
!        Generate the elementary reflector H(i) to annihilate
!        A(k+i+1:n,i)
!
         CALL SLARFG( N-K-I+1, A( K+I, I ), A( MIN( K+I+1, N ), I ), 1, &
                      TAU( I ) )
         EI = A( K+I, I )
         A( K+I, I ) = ONE
!
!        Compute  Y(1:n,i)
!
         CALL SGEMV( 'No transpose', N, N-K-I+1, ONE, A( 1, I+1 ), LDA, &
                     A( K+I, I ), 1, ZERO, Y( 1, I ), 1 )
         CALL SGEMV( 'Transpose', N-K-I+1, I-1, ONE, A( K+I, 1 ), LDA, &
                     A( K+I, I ), 1, ZERO, T( 1, I ), 1 )
         CALL SGEMV( 'No transpose', N, I-1, -ONE, Y, LDY, T( 1, I ), 1, &
                     ONE, Y( 1, I ), 1 )
         CALL SSCAL( N, TAU( I ), Y( 1, I ), 1 )
!
!        Compute T(1:i,i)
!
         CALL SSCAL( I-1, -TAU( I ), T( 1, I ), 1 )
         CALL STRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, LDT, &
                     T( 1, I ), 1 )
         T( I, I ) = TAU( I )
!
   10 CONTINUE
      A( K+NB, NB ) = EI
!
      RETURN
!
!     End of SLAHRD
!
      END SUBROUTINE SLAHRD

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SLALN2( LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, B, &
                         LDB, WR, WI, X, LDX, SCALE, XNORM, INFO )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            LTRANS
      INTEGER            INFO, LDA, LDB, LDX, NA, NW
      REAL(r8)               CA, D1, D2, SCALE, SMIN, WI, WR, XNORM
!     ..
!     .. Array Arguments ..
      REAL(r8)               A( LDA, * ), B( LDB, * ), X( LDX, * )
!     ..
!
!  Purpose
!  =======
!
!  SLALN2 solves a system of the form  (ca A - w D ) X = s B
!  or (ca A' - w D) X = s B   with possible scaling ("s") and
!  perturbation of A.  (A' means A-transpose.)
!
!  A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA
!  real diagonal matrix, w is a real or complex value, and X and B are
!  NA x 1 matrices -- real if w is real, complex if w is complex.  NA
!  may be 1 or 2.
!
!  If w is complex, X and B are represented as NA x 2 matrices,
!  the first column of each being the real part and the second
!  being the imaginary part.
!
!  "s" is a scaling factor (.LE. 1), computed by SLALN2, which is
!  so chosen that X can be computed without overflow.  X is further
!  scaled if necessary to assure that norm(ca A - w D)*norm(X) is less
!  than overflow.
!
!  If both singular values of (ca A - w D) are less than SMIN,
!  SMIN*identity will be used instead of (ca A - w D).  If only one
!  singular value is less than SMIN, one element of (ca A - w D) will be
!  perturbed enough to make the smallest singular value roughly SMIN.
!  If both singular values are at least SMIN, (ca A - w D) will not be
!  perturbed.  In any case, the perturbation will be at most some small
!  multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values
!  are computed by infinity-norm approximations, and thus will only be
!  correct to a factor of 2 or so.
!
!  Note: all input quantities are assumed to be smaller than overflow
!  by a reasonable factor.  (See BIGNUM.)
!
!  Arguments
!  ==========
!
!  LTRANS  (input) LOGICAL
!          =.TRUE.:  A-transpose will be used.
!          =.FALSE.: A will be used (not transposed.)
!
!  NA      (input) INTEGER
!          The size of the matrix A.  It may (only) be 1 or 2.
!
!  NW      (input) INTEGER
!          1 if "w" is real, 2 if "w" is complex.  It may only be 1
!          or 2.
!
!  SMIN    (input) REAL
!          The desired lower bound on the singular values of A.  This
!          should be a safe distance away from underflow or overflow,
!          say, between (underflow/machine precision) and  (machine
!          precision * overflow ).  (See BIGNUM and ULP.)
!
!  CA      (input) REAL
!          The coefficient c, which A is multiplied by.
!
!  A       (input) REAL array, dimension (LDA,NA)
!          The NA x NA matrix A.
!
!  LDA     (input) INTEGER
!          The leading dimension of A.  It must be at least NA.
!
!  D1      (input) REAL
!          The 1,1 element in the diagonal matrix D.
!
!  D2      (input) REAL
!          The 2,2 element in the diagonal matrix D.  Not used if NW=1.
!
!  B       (input) REAL array, dimension (LDB,NW)
!          The NA x NW matrix B (right-hand side).  If NW=2 ("w" is
!          complex), column 1 contains the real part of B and column 2
!          contains the imaginary part.
!
!  LDB     (input) INTEGER
!          The leading dimension of B.  It must be at least NA.
!
!  WR      (input) REAL
!          The real part of the scalar "w".
!
!  WI      (input) REAL
!          The imaginary part of the scalar "w".  Not used if NW=1.
!
!  X       (output) REAL array, dimension (LDX,NW)
!          The NA x NW matrix X (unknowns), as computed by SLALN2.
!          If NW=2 ("w" is complex), on exit, column 1 will contain
!          the real part of X and column 2 will contain the imaginary
!          part.
!
!  LDX     (input) INTEGER
!          The leading dimension of X.  It must be at least NA.
!
!  SCALE   (output) REAL
!          The scale factor that B must be multiplied by to insure
!          that overflow does not occur when computing X.  Thus,
!          (ca A - w D) X  will be SCALE*B, not B (ignoring
!          perturbations of A.)  It will be at most 1.
!
!  XNORM   (output) REAL
!          The infinity-norm of X, when X is regarded as an NA x NW
!          real matrix.
!
!  INFO    (output) INTEGER
!          An error flag.  It will be set to zero if no error occurs,
!          a negative number if an argument is in error, or a positive
!          number if  ca A - w D  had to be perturbed.
!          The possible values are:
!          = 0: No error occurred, and (ca A - w D) did not have to be
!                 perturbed.
!          = 1: (ca A - w D) had to be perturbed to make its smallest
!               (or only) singular value greater than SMIN.
!          NOTE: In the interests of speed, this routine does not
!                check the inputs for errors.
!
! =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0_r8, ONE = 1.0E0_r8 )
      REAL(r8)               TWO
      PARAMETER          ( TWO = 2.0E0_r8 )
!     ..
!     .. Local Scalars ..
      INTEGER            ICMAX, J
      REAL(r8)               BBND, BI1, BI2, BIGNUM, BNORM, BR1, BR2, CI21, &
                         CI22, CMAX, CNORM, CR21, CR22, CSI, CSR, LI21, &
                         LR21, SMINI, SMLNUM, TEMP, U22ABS, UI11, UI11R, &
                         UI12, UI12S, UI22, UR11, UR11R, UR12, UR12S, &
                         UR22, XI1, XI2, XR1, XR2
!     ..
!     .. Local Arrays ..
      LOGICAL            CSWAP( 4 ), RSWAP( 4 )
      INTEGER            IPIVOT( 4, 4 )
      REAL(r8)               CI( 2, 2 ), CIV( 4 ), CR( 2, 2 ), CRV( 4 )
!     ..
!     .. Equivalences ..
      EQUIVALENCE        ( CI( 1, 1 ), CIV( 1 ) ), &
                         ( CR( 1, 1 ), CRV( 1 ) )
!     ..
!     .. Data statements ..
      DATA               CSWAP / .FALSE., .FALSE., .TRUE., .TRUE. /
      DATA               RSWAP / .FALSE., .TRUE., .FALSE., .TRUE. /
      DATA               IPIVOT / 1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4, &
                         3, 2, 1 /
!     ..
!     .. Executable Statements ..
!
!     Compute BIGNUM
!
      SMLNUM = TWO*SLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      SMINI = MAX( SMIN, SMLNUM )
!
!     Don't check for input errors
!
      INFO = 0
!
!     Standard Initializations
!
      SCALE = ONE
!
      IF( NA.EQ.1 ) THEN
!
!        1 x 1  (i.e., scalar) system   C X = B
!
         IF( NW.EQ.1 ) THEN
!
!           Real 1x1 system.
!
!           C = ca A - w D
!
            CSR = CA*A( 1, 1 ) - WR*D1
            CNORM = ABS( CSR )
!
!           If | C | < SMINI, use C = SMINI
!
            IF( CNORM.LT.SMINI ) THEN
               CSR = SMINI
               CNORM = SMINI
               INFO = 1
            END IF
!
!           Check scaling for  X = B / C
!
            BNORM = ABS( B( 1, 1 ) )
            IF( CNORM.LT.ONE .AND. BNORM.GT.ONE ) THEN
               IF( BNORM.GT.BIGNUM*CNORM ) &
                  SCALE = ONE / BNORM
            END IF
!
!           Compute X
!
            X( 1, 1 ) = ( B( 1, 1 )*SCALE ) / CSR
            XNORM = ABS( X( 1, 1 ) )
         ELSE
!
!           Complex 1x1 system (w is complex)
!
!           C = ca A - w D
!
            CSR = CA*A( 1, 1 ) - WR*D1
            CSI = -WI*D1
            CNORM = ABS( CSR ) + ABS( CSI )
!
!           If | C | < SMINI, use C = SMINI
!
            IF( CNORM.LT.SMINI ) THEN
               CSR = SMINI
               CSI = ZERO
               CNORM = SMINI
               INFO = 1
            END IF
!
!           Check scaling for  X = B / C
!
            BNORM = ABS( B( 1, 1 ) ) + ABS( B( 1, 2 ) )
            IF( CNORM.LT.ONE .AND. BNORM.GT.ONE ) THEN
               IF( BNORM.GT.BIGNUM*CNORM ) &
                  SCALE = ONE / BNORM
            END IF
!
!           Compute X
!
            CALL SLADIV( SCALE*B( 1, 1 ), SCALE*B( 1, 2 ), CSR, CSI, &
                         X( 1, 1 ), X( 1, 2 ) )
            XNORM = ABS( X( 1, 1 ) ) + ABS( X( 1, 2 ) )
         END IF
!
      ELSE
!
!        2x2 System
!
!        Compute the real part of  C = ca A - w D  (or  ca A' - w D )
!
         CR( 1, 1 ) = CA*A( 1, 1 ) - WR*D1
         CR( 2, 2 ) = CA*A( 2, 2 ) - WR*D2
         IF( LTRANS ) THEN
            CR( 1, 2 ) = CA*A( 2, 1 )
            CR( 2, 1 ) = CA*A( 1, 2 )
         ELSE
            CR( 2, 1 ) = CA*A( 2, 1 )
            CR( 1, 2 ) = CA*A( 1, 2 )
         END IF
!
         IF( NW.EQ.1 ) THEN
!
!           Real 2x2 system  (w is real)
!
!           Find the largest element in C
!
            CMAX = ZERO
            ICMAX = 0
!
            DO 10 J = 1, 4
               IF( ABS( CRV( J ) ).GT.CMAX ) THEN
                  CMAX = ABS( CRV( J ) )
                  ICMAX = J
               END IF
   10       CONTINUE
!
!           If norm(C) < SMINI, use SMINI*identity.
!
            IF( CMAX.LT.SMINI ) THEN
               BNORM = MAX( ABS( B( 1, 1 ) ), ABS( B( 2, 1 ) ) )
               IF( SMINI.LT.ONE .AND. BNORM.GT.ONE ) THEN
                  IF( BNORM.GT.BIGNUM*SMINI ) &
                     SCALE = ONE / BNORM
               END IF
               TEMP = SCALE / SMINI
               X( 1, 1 ) = TEMP*B( 1, 1 )
               X( 2, 1 ) = TEMP*B( 2, 1 )
               XNORM = TEMP*BNORM
               INFO = 1
               RETURN
            END IF
!
!           Gaussian elimination with complete pivoting.
!
            UR11 = CRV( ICMAX )
            CR21 = CRV( IPIVOT( 2, ICMAX ) )
            UR12 = CRV( IPIVOT( 3, ICMAX ) )
            CR22 = CRV( IPIVOT( 4, ICMAX ) )
            UR11R = ONE / UR11
            LR21 = UR11R*CR21
            UR22 = CR22 - UR12*LR21
!
!           If smaller pivot < SMINI, use SMINI
!
            IF( ABS( UR22 ).LT.SMINI ) THEN
               UR22 = SMINI
               INFO = 1
            END IF
            IF( RSWAP( ICMAX ) ) THEN
               BR1 = B( 2, 1 )
               BR2 = B( 1, 1 )
            ELSE
               BR1 = B( 1, 1 )
               BR2 = B( 2, 1 )
            END IF
            BR2 = BR2 - LR21*BR1
            BBND = MAX( ABS( BR1*( UR22*UR11R ) ), ABS( BR2 ) )
            IF( BBND.GT.ONE .AND. ABS( UR22 ).LT.ONE ) THEN
               IF( BBND.GE.BIGNUM*ABS( UR22 ) ) &
                  SCALE = ONE / BBND
            END IF
!
            XR2 = ( BR2*SCALE ) / UR22
            XR1 = ( SCALE*BR1 )*UR11R - XR2*( UR11R*UR12 )
            IF( CSWAP( ICMAX ) ) THEN
               X( 1, 1 ) = XR2
               X( 2, 1 ) = XR1
            ELSE
               X( 1, 1 ) = XR1
               X( 2, 1 ) = XR2
            END IF
            XNORM = MAX( ABS( XR1 ), ABS( XR2 ) )
!
!           Further scaling if  norm(A) norm(X) > overflow
!
            IF( XNORM.GT.ONE .AND. CMAX.GT.ONE ) THEN
               IF( XNORM.GT.BIGNUM / CMAX ) THEN
                  TEMP = CMAX / BIGNUM
                  X( 1, 1 ) = TEMP*X( 1, 1 )
                  X( 2, 1 ) = TEMP*X( 2, 1 )
                  XNORM = TEMP*XNORM
                  SCALE = TEMP*SCALE
               END IF
            END IF
         ELSE
!
!           Complex 2x2 system  (w is complex)
!
!           Find the largest element in C
!
            CI( 1, 1 ) = -WI*D1
            CI( 2, 1 ) = ZERO
            CI( 1, 2 ) = ZERO
            CI( 2, 2 ) = -WI*D2
            CMAX = ZERO
            ICMAX = 0
!
            DO 20 J = 1, 4
               IF( ABS( CRV( J ) )+ABS( CIV( J ) ).GT.CMAX ) THEN
                  CMAX = ABS( CRV( J ) ) + ABS( CIV( J ) )
                  ICMAX = J
               END IF
   20       CONTINUE
!
!           If norm(C) < SMINI, use SMINI*identity.
!
            IF( CMAX.LT.SMINI ) THEN
               BNORM = MAX( ABS( B( 1, 1 ) )+ABS( B( 1, 2 ) ), &
                       ABS( B( 2, 1 ) )+ABS( B( 2, 2 ) ) )
               IF( SMINI.LT.ONE .AND. BNORM.GT.ONE ) THEN
                  IF( BNORM.GT.BIGNUM*SMINI ) &
                     SCALE = ONE / BNORM
               END IF
               TEMP = SCALE / SMINI
               X( 1, 1 ) = TEMP*B( 1, 1 )
               X( 2, 1 ) = TEMP*B( 2, 1 )
               X( 1, 2 ) = TEMP*B( 1, 2 )
               X( 2, 2 ) = TEMP*B( 2, 2 )
               XNORM = TEMP*BNORM
               INFO = 1
               RETURN
            END IF
!
!           Gaussian elimination with complete pivoting.
!
            UR11 = CRV( ICMAX )
            UI11 = CIV( ICMAX )
            CR21 = CRV( IPIVOT( 2, ICMAX ) )
            CI21 = CIV( IPIVOT( 2, ICMAX ) )
            UR12 = CRV( IPIVOT( 3, ICMAX ) )
            UI12 = CIV( IPIVOT( 3, ICMAX ) )
            CR22 = CRV( IPIVOT( 4, ICMAX ) )
            CI22 = CIV( IPIVOT( 4, ICMAX ) )
            IF( ICMAX.EQ.1 .OR. ICMAX.EQ.4 ) THEN
!
!              Code when off-diagonals of pivoted C are real
!
               IF( ABS( UR11 ).GT.ABS( UI11 ) ) THEN
                  TEMP = UI11 / UR11
                  UR11R = ONE / ( UR11*( ONE+TEMP**2 ) )
                  UI11R = -TEMP*UR11R
               ELSE
                  TEMP = UR11 / UI11
                  UI11R = -ONE / ( UI11*( ONE+TEMP**2 ) )
                  UR11R = -TEMP*UI11R
               END IF
               LR21 = CR21*UR11R
               LI21 = CR21*UI11R
               UR12S = UR12*UR11R
               UI12S = UR12*UI11R
               UR22 = CR22 - UR12*LR21
               UI22 = CI22 - UR12*LI21
            ELSE
!
!              Code when diagonals of pivoted C are real
!
               UR11R = ONE / UR11
               UI11R = ZERO
               LR21 = CR21*UR11R
               LI21 = CI21*UR11R
               UR12S = UR12*UR11R
               UI12S = UI12*UR11R
               UR22 = CR22 - UR12*LR21 + UI12*LI21
               UI22 = -UR12*LI21 - UI12*LR21
            END IF
            U22ABS = ABS( UR22 ) + ABS( UI22 )
!
!           If smaller pivot < SMINI, use SMINI
!
            IF( U22ABS.LT.SMINI ) THEN
               UR22 = SMINI
               UI22 = ZERO
               INFO = 1
            END IF
            IF( RSWAP( ICMAX ) ) THEN
               BR2 = B( 1, 1 )
               BR1 = B( 2, 1 )
               BI2 = B( 1, 2 )
               BI1 = B( 2, 2 )
            ELSE
               BR1 = B( 1, 1 )
               BR2 = B( 2, 1 )
               BI1 = B( 1, 2 )
               BI2 = B( 2, 2 )
            END IF
            BR2 = BR2 - LR21*BR1 + LI21*BI1
            BI2 = BI2 - LI21*BR1 - LR21*BI1
            BBND = MAX( ( ABS( BR1 )+ABS( BI1 ) )* &
                   ( U22ABS*( ABS( UR11R )+ABS( UI11R ) ) ), &
                   ABS( BR2 )+ABS( BI2 ) )
            IF( BBND.GT.ONE .AND. U22ABS.LT.ONE ) THEN
               IF( BBND.GE.BIGNUM*U22ABS ) THEN
                  SCALE = ONE / BBND
                  BR1 = SCALE*BR1
                  BI1 = SCALE*BI1
                  BR2 = SCALE*BR2
                  BI2 = SCALE*BI2
               END IF
            END IF
!
            CALL SLADIV( BR2, BI2, UR22, UI22, XR2, XI2 )
            XR1 = UR11R*BR1 - UI11R*BI1 - UR12S*XR2 + UI12S*XI2
            XI1 = UI11R*BR1 + UR11R*BI1 - UI12S*XR2 - UR12S*XI2
            IF( CSWAP( ICMAX ) ) THEN
               X( 1, 1 ) = XR2
               X( 2, 1 ) = XR1
               X( 1, 2 ) = XI2
               X( 2, 2 ) = XI1
            ELSE
               X( 1, 1 ) = XR1
               X( 2, 1 ) = XR2
               X( 1, 2 ) = XI1
               X( 2, 2 ) = XI2
            END IF
            XNORM = MAX( ABS( XR1 )+ABS( XI1 ), ABS( XR2 )+ABS( XI2 ) )
!
!           Further scaling if  norm(A) norm(X) > overflow
!
            IF( XNORM.GT.ONE .AND. CMAX.GT.ONE ) THEN
               IF( XNORM.GT.BIGNUM / CMAX ) THEN
                  TEMP = CMAX / BIGNUM
                  X( 1, 1 ) = TEMP*X( 1, 1 )
                  X( 2, 1 ) = TEMP*X( 2, 1 )
                  X( 1, 2 ) = TEMP*X( 1, 2 )
                  X( 2, 2 ) = TEMP*X( 2, 2 )
                  XNORM = TEMP*XNORM
                  SCALE = TEMP*SCALE
               END IF
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of SLALN2
!
      END SUBROUTINE SLALN2

!--------1---------2---------3---------4---------5---------6---------7---------8

      FUNCTION SLANGE( NORM, M, N, A, LDA, WORK )

      real(r8) slange
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, M, N
!     ..
!     .. Array Arguments ..
      REAL(r8)               A( LDA, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  SLANGE  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  real matrix A.
!
!  Description
!  ===========
!
!  SLANGE returns the value
!
!     SLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in SLANGE as described
!          above.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.  When M = 0,
!          SLANGE is set to zero.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.  When N = 0,
!          SLANGE is set to zero.
!
!  A       (input) REAL array, dimension (LDA,N)
!          The m by n matrix A.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(M,1).
!
!  WORK    (workspace) REAL array, dimension (LWORK),
!          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
!          referenced.
!
! =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0_r8, ZERO = 0.0E+0_r8 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      REAL(r8)               SCALE, SUM, VALUE
!     ..
!     .. Executable Statements ..
!
      IF( MIN( M, N ).EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
         VALUE = ZERO
         DO 20 J = 1, N
            DO 10 I = 1, M
               VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
!
!        Find norm1(A).
!
         VALUE = ZERO
         DO 40 J = 1, N
            SUM = ZERO
            DO 30 I = 1, M
               SUM = SUM + ABS( A( I, J ) )
   30       CONTINUE
            VALUE = MAX( VALUE, SUM )
   40    CONTINUE
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
!
!        Find normI(A).
!
         DO 50 I = 1, M
            WORK( I ) = ZERO
   50    CONTINUE
         DO 70 J = 1, N
            DO 60 I = 1, M
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = ZERO
         DO 80 I = 1, M
            VALUE = MAX( VALUE, WORK( I ) )
   80    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         DO 90 J = 1, N
            CALL SLASSQ( M, A( 1, J ), 1, SCALE, SUM )
   90    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
!
      SLANGE = VALUE
      RETURN
!
!     End of SLANGE
!
      END FUNCTION SLANGE

!--------1---------2---------3---------4---------5---------6---------7---------8

      FUNCTION SLANHS( NORM, N, A, LDA, WORK )

      real(r8) slanhs

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, N
!     ..
!     .. Array Arguments ..
      REAL(r8)               A( LDA, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  SLANHS  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  Hessenberg matrix A.
!
!  Description
!  ===========
!
!  SLANHS returns the value
!
!     SLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in SLANHS as described
!          above.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, SLANHS is
!          set to zero.
!
!  A       (input) REAL array, dimension (LDA,N)
!          The n by n upper Hessenberg matrix A; the part of A below the
!          first sub-diagonal is not referenced.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(N,1).
!
!  WORK    (workspace) REAL array, dimension (LWORK),
!          where LWORK >= N when NORM = 'I'; otherwise, WORK is not
!          referenced.
!
! =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0_r8, ZERO = 0.0E+0_r8 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      REAL(r8)               SCALE, SUM, VALUE
!     ..
!     .. Executable Statements ..
!
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
         VALUE = ZERO
         DO 20 J = 1, N
            DO 10 I = 1, MIN( N, J+1 )
               VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
!
!        Find norm1(A).
!
         VALUE = ZERO
         DO 40 J = 1, N
            SUM = ZERO
            DO 30 I = 1, MIN( N, J+1 )
               SUM = SUM + ABS( A( I, J ) )
   30       CONTINUE
            VALUE = MAX( VALUE, SUM )
   40    CONTINUE
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
!
!        Find normI(A).
!
         DO 50 I = 1, N
            WORK( I ) = ZERO
   50    CONTINUE
         DO 70 J = 1, N
            DO 60 I = 1, MIN( N, J+1 )
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = ZERO
         DO 80 I = 1, N
            VALUE = MAX( VALUE, WORK( I ) )
   80    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         DO 90 J = 1, N
            CALL SLASSQ( MIN( N, J+1 ), A( 1, J ), 1, SCALE, SUM )
   90    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF 
!
      SLANHS = VALUE
      RETURN
!
!     End of SLANHS
!
      END FUNCTION SLANHS

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      REAL(r8)               A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN
!     ..
!
!  Purpose
!  =======
!
!  SLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric
!  matrix in standard form:
!
!       [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
!       [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]
!
!  where either
!  1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or
!  2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex
!  conjugate eigenvalues.
!
!  Arguments
!  =========
!
!  A       (input/output) REAL
!  B       (input/output) REAL
!  C       (input/output) REAL
!  D       (input/output) REAL
!          On entry, the elements of the input matrix.
!          On exit, they are overwritten by the elements of the
!          standardised Schur form.
!
!  RT1R    (output) REAL
!  RT1I    (output) REAL
!  RT2R    (output) REAL
!  RT2I    (output) REAL
!          The real and imaginary parts of the eigenvalues. If the
!          eigenvalues are both real, abs(RT1R) >= abs(RT2R); if the
!          eigenvalues are a complex conjugate pair, RT1I > 0.
!
!  CS      (output) REAL
!  SN      (output) REAL
!          Parameters of the rotation matrix.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0E+0_r8, HALF = 0.5E+0_r8, ONE = 1.0E+0_r8 )
!     ..
!     .. Local Scalars ..
      REAL(r8)               AA, BB, CC, CS1, DD, P, SAB, SAC, SIGMA, SN1, &
                         TAU, TEMP
!     ..
!     .. Executable Statements ..
!
!     Initialize CS and SN
!
      CS = ONE
      SN = ZERO
!
      IF( C.EQ.ZERO ) THEN
         GO TO 10
!
      ELSE IF( B.EQ.ZERO ) THEN
!
!        Swap rows and columns
!
         CS = ZERO
         SN = ONE
         TEMP = D
         D = A
         A = TEMP
         B = -C
         C = ZERO
         GO TO 10
      ELSE IF( (A-D).EQ.ZERO .AND. SIGN( ONE, B ).NE. &
         SIGN( ONE, C ) ) THEN
         GO TO 10
      ELSE
!
!        Make diagonal elements equal
!
         TEMP = A - D
         P = HALF*TEMP
         SIGMA = B + C
         TAU = SLAPY2( SIGMA, TEMP )
         CS1 = SQRT( HALF*( ONE+ABS( SIGMA ) / TAU ) )
         SN1 = -( P / ( TAU*CS1 ) )*SIGN( ONE, SIGMA )
!
!        Compute [ AA  BB ] = [ A  B ] [ CS1 -SN1 ]
!                [ CC  DD ]   [ C  D ] [ SN1  CS1 ]
!
         AA = A*CS1 + B*SN1
         BB = -A*SN1 + B*CS1
         CC = C*CS1 + D*SN1
         DD = -C*SN1 + D*CS1
!
!        Compute [ A  B ] = [ CS1  SN1 ] [ AA  BB ]
!                [ C  D ]   [-SN1  CS1 ] [ CC  DD ]
!
         A = AA*CS1 + CC*SN1
         B = BB*CS1 + DD*SN1
         C = -AA*SN1 + CC*CS1
         D = -BB*SN1 + DD*CS1
!
!        Accumulate transformation
!
         TEMP = CS*CS1 - SN*SN1
         SN = CS*SN1 + SN*CS1
         CS = TEMP
!
         TEMP = HALF*( A+D )
         A = TEMP
         D = TEMP
!
         IF( C.NE.ZERO ) THEN
            IF ( B.NE.ZERO ) THEN
               IF( SIGN( ONE, B ).EQ.SIGN( ONE, C ) ) THEN
!
!                 Real eigenvalues: reduce to upper triangular form
!
                  SAB = SQRT( ABS( B ) )
                  SAC = SQRT( ABS( C ) )
                  P = SIGN( SAB*SAC, C )
                  TAU = ONE / SQRT( ABS( B+C ) )
                  A = TEMP + P
                  D = TEMP - P
                  B = B - C
                  C = ZERO
                  CS1 = SAB*TAU
                  SN1 = SAC*TAU
                  TEMP = CS*CS1 - SN*SN1
                  SN = CS*SN1 + SN*CS1
                  CS = TEMP
               END IF
            ELSE
               B = -C
               C = ZERO
               TEMP = CS
               CS = -SN
               SN = TEMP
            ENDIF
         ENDIF
      END IF
!
   10 CONTINUE
!
!     Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).
!
      RT1R = A
      RT2R = D
      IF( C.EQ.ZERO ) THEN
         RT1I = ZERO
         RT2I = ZERO
      ELSE
         RT1I = SQRT( ABS( B ) )*SQRT( ABS( C ) )
         RT2I = -RT1I
      END IF
      RETURN
!
!     End of SLANV2
!
      END SUBROUTINE SLANV2

!--------1---------2---------3---------4---------5---------6---------7---------8

      FUNCTION SLAPY2( X, Y )

      real(r8) slapy2

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      REAL(r8)               X, Y
!     ..
!
!  Purpose
!  =======
!
!  SLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
!  overflow.
!
!  Arguments
!  =========
!
!  X       (input) REAL
!  Y       (input) REAL
!          X and Y specify the values x and y.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ZERO
      PARAMETER          ( ZERO = 0.0E0_r8 )
      REAL(r8)               ONE
      PARAMETER          ( ONE = 1.0E0_r8 )
!     ..
!     .. Local Scalars ..
      REAL(r8)               W, XABS, YABS, Z
!     ..
!     .. Executable Statements ..
!
      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z.EQ.ZERO ) THEN
         SLAPY2 = W
      ELSE
         SLAPY2 = W*SQRT( ONE+( Z / W )**2 )
      END IF
      RETURN
!
!     End of SLAPY2
!
      END FUNCTION SLAPY2

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      REAL(r8)               TAU
!     ..
!     .. Array Arguments ..
      REAL(r8)               C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  SLARF applies a real elementary reflector H to a real m by n matrix
!  C, from either the left or the right. H is represented in the form
!
!        H = I - tau * v * v'
!
!  where tau is a real scalar and v is a real vector.
!
!  If tau = 0, then H is taken to be the unit matrix.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': form  H * C
!          = 'R': form  C * H
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  V       (input) REAL array, dimension
!                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!          The vector v in the representation of H. V is not used if
!          TAU = 0.
!
!  INCV    (input) INTEGER
!          The increment between elements of v. INCV <> 0.
!
!  TAU     (input) REAL
!          The value tau in the representation of H.
!
!  C       (input/output) REAL array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!          or C * H if SIDE = 'R'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) REAL array, dimension
!                         (N) if SIDE = 'L'
!                      or (M) if SIDE = 'R'
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0_r8, ZERO = 0.0E+0_r8 )
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( SIDE, 'L' ) ) THEN
!
!        Form  H * C
!
         IF( TAU.NE.ZERO ) THEN
!
!           w := C' * v
!
            CALL SGEMV( 'Transpose', M, N, ONE, C, LDC, V, INCV, ZERO, &
                        WORK, 1 )
!
!           C := C - v * w'
!
            CALL SGER( M, N, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
!
!        Form  C * H
!
         IF( TAU.NE.ZERO ) THEN
!
!           w := C * v
!
            CALL SGEMV( 'No transpose', M, N, ONE, C, LDC, V, INCV, &
                        ZERO, WORK, 1 )
!
!           C := C - w * v'
!
            CALL SGER( M, N, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
!
!     End of SLARF
!
      END SUBROUTINE SLARF

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, &
                         T, LDT, C, LDC, WORK, LDWORK )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
!     ..
!     .. Array Arguments ..
      REAL(r8)               C( LDC, * ), T( LDT, * ), V( LDV, * ), &
                         WORK( LDWORK, * )
!     ..
!
!  Purpose
!  =======
!
!  SLARFB applies a real block reflector H or its transpose H' to a
!  real m by n matrix C, from either the left or the right.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply H or H' from the Left
!          = 'R': apply H or H' from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply H (No transpose)
!          = 'T': apply H' (Transpose)
!
!  DIRECT  (input) CHARACTER*1
!          Indicates how H is formed from a product of elementary
!          reflectors
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
!  STOREV  (input) CHARACTER*1
!          Indicates how the vectors which define the elementary
!          reflectors are stored:
!          = 'C': Columnwise
!          = 'R': Rowwise
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  K       (input) INTEGER
!          The order of the matrix T (= the number of elementary
!          reflectors whose product defines the block reflector).
!
!  V       (input) REAL array, dimension
!                                (LDV,K) if STOREV = 'C'
!                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!          The matrix V. See further details.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
!          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
!          if STOREV = 'R', LDV >= K.
!
!  T       (input) REAL array, dimension (LDT,K)
!          The triangular k by k matrix T in the representation of the
!          block reflector.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= K.
!
!  C       (input/output) REAL array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDA >= max(1,M).
!
!  WORK    (workspace) REAL array, dimension (LDWORK,K)
!
!  LDWORK  (input) INTEGER
!          The leading dimension of the array WORK.
!          If SIDE = 'L', LDWORK >= max(1,N);
!          if SIDE = 'R', LDWORK >= max(1,M).
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ONE
      PARAMETER          ( ONE = 1.0E+0_r8 )
!     ..
!     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( M.LE.0 .OR. N.LE.0 ) &
         RETURN
!
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      END IF
!
      IF( LSAME( STOREV, 'C' ) ) THEN
!
         IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1 )    (first K rows)
!                     ( V2 )
!           where  V1  is unit lower triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
!
!              W := C1'
!
               DO 10 J = 1, K
                  CALL SCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
   10          CONTINUE
!
!              W := W * V1
!
               CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                           K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C2'*V2
!
                  CALL SGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                              ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV, &
                              ONE, WORK, LDWORK )
               END IF
!
!              W := W * T'  or  W * T
!
               CALL STRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W'
!
               IF( M.GT.K ) THEN
!
!                 C2 := C2 - V2 * W'
!
                  CALL SGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                              -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE, &
                              C( K+1, 1 ), LDC )
               END IF
!
!              W := W * V1'
!
               CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
                           ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W'
!
               DO 30 J = 1, K
                  DO 20 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
   20             CONTINUE
   30          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C1
!
               DO 40 J = 1, K
                  CALL SCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
!
!              W := W * V1
!
               CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                           K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C2 * V2
!
                  CALL SGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                              ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV, &
                              ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T'
!
               CALL STRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V'
!
               IF( N.GT.K ) THEN
!
!                 C2 := C2 - W * V2'
!
                  CALL SGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                              -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE, &
                              C( 1, K+1 ), LDC )
               END IF
!
!              W := W * V1'
!
               CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                           ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 60 J = 1, K
                  DO 50 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
!
         ELSE
!
!           Let  V =  ( V1 )
!                     ( V2 )    (last K rows)
!           where  V2  is unit upper triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
!
!              W := C2'
!
               DO 70 J = 1, K
                  CALL SCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
   70          CONTINUE
!
!              W := W * V2
!
               CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                           K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C1'*V1
!
                  CALL SGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                              ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T'  or  W * T
!
               CALL STRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W'
!
               IF( M.GT.K ) THEN
!
!                 C1 := C1 - V1 * W'
!
                  CALL SGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                              -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
!
!              W := W * V2'
!
               CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                           ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W'
!
               DO 90 J = 1, K
                  DO 80 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
   80             CONTINUE
   90          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C2
!
               DO 100 J = 1, K
                  CALL SCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
!
!              W := W * V2
!
               CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                           K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C1 * V1
!
                  CALL SGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                              ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T'
!
               CALL STRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V'
!
               IF( N.GT.K ) THEN
!
!                 C1 := C1 - W * V1'
!
                  CALL SGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                              -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
!
!              W := W * V2'
!
               CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                           ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W
!
               DO 120 J = 1, K
                  DO 110 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
!
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
!
         IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1  V2 )    (V1: first K columns)
!           where  V1  is unit upper triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
!
!              W := C1'
!
               DO 130 J = 1, K
                  CALL SCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
  130          CONTINUE
!
!              W := W * V1'
!
               CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                           ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C2'*V2'
!
                  CALL SGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                              C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE, &
                              WORK, LDWORK )
               END IF
!
!              W := W * T'  or  W * T
!
               CALL STRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V' * W'
!
               IF( M.GT.K ) THEN
!
!                 C2 := C2 - V2' * W'
!
                  CALL SGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                              V( 1, K+1 ), LDV, WORK, LDWORK, ONE, &
                              C( K+1, 1 ), LDC )
               END IF
!
!              W := W * V1
!
               CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                           K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W'
!
               DO 150 J = 1, K
                  DO 140 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
  140             CONTINUE
  150          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
!
!              W := C1
!
               DO 160 J = 1, K
                  CALL SCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
!
!              W := W * V1'
!
               CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                           ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C2 * V2'
!
                  CALL SGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                              ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV, &
                              ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T'
!
               CALL STRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
               IF( N.GT.K ) THEN
!
!                 C2 := C2 - W * V2
!
                  CALL SGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                              -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE, &
                              C( 1, K+1 ), LDC )
               END IF
!
!              W := W * V1
!
               CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                           K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 180 J = 1, K
                  DO 170 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
!
            END IF
!
         ELSE
!
!           Let  V =  ( V1  V2 )    (V2: last K columns)
!           where  V2  is unit lower triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
!
!              W := C2'
!
               DO 190 J = 1, K
                  CALL SCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
  190          CONTINUE
!
!              W := W * V2'
!
               CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
                           ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C1'*V1'
!
                  CALL SGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                              C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T'  or  W * T
!
               CALL STRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V' * W'
!
               IF( M.GT.K ) THEN
!
!                 C1 := C1 - V1' * W'
!
                  CALL SGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                              V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
!
!              W := W * V2
!
               CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                           K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W'
!
               DO 210 J = 1, K
                  DO 200 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
  200             CONTINUE
  210          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
!
!              W := C2
!
               DO 220 J = 1, K
                  CALL SCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
!
!              W := W * V2'
!
               CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                           ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C1 * V1'
!
                  CALL SGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                              ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T'
!
               CALL STRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                           ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
               IF( N.GT.K ) THEN
!
!                 C1 := C1 - W * V1
!
                  CALL SGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                              -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
!
!              W := W * V2
!
               CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                           K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 240 J = 1, K
                  DO 230 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
!
            END IF
!
         END IF
      END IF
!
      RETURN
!
!     End of SLARFB
!
      END SUBROUTINE SLARFB

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SLARFG( N, ALPHA, X, INCX, TAU )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      REAL(r8)               ALPHA, TAU
!     ..
!     .. Array Arguments ..
      REAL(r8)               X( * )
!     ..
!
!  Purpose
!  =======
!
!  SLARFG generates a real elementary reflector H of order n, such
!  that
!
!        H * ( alpha ) = ( beta ),   H' * H = I.
!            (   x   )   (   0  )
!
!  where alpha and beta are scalars, and x is an (n-1)-element real
!  vector. H is represented in the form
!
!        H = I - tau * ( 1 ) * ( 1 v' ) ,
!                      ( v )
!
!  where tau is a real scalar and v is a real (n-1)-element
!  vector.
!
!  If the elements of x are all zero, then tau = 0 and H is taken to be
!  the unit matrix.
!
!  Otherwise  1 <= tau <= 2.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the elementary reflector.
!
!  ALPHA   (input/output) REAL
!          On entry, the value alpha.
!          On exit, it is overwritten with the value beta.
!
!  X       (input/output) REAL array, dimension
!                         (1+(N-2)*abs(INCX))
!          On entry, the vector x.
!          On exit, it is overwritten with the vector v.
!
!  INCX    (input) INTEGER
!          The increment between elements of X. INCX > 0.
!
!  TAU     (output) REAL
!          The value tau.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0_r8, ZERO = 0.0E+0_r8 )
!     ..
!     .. Local Scalars ..
      INTEGER            J, KNT
      REAL(r8)               BETA, RSAFMN, SAFMIN, XNORM
!     ..
!     .. Executable Statements ..
!
      IF( N.LE.1 ) THEN
         TAU = ZERO
         RETURN
      END IF
!
      XNORM = SNRM2( N-1, X, INCX )
!
      IF( XNORM.EQ.ZERO ) THEN
!
!        H  =  I
!
         TAU = ZERO
      ELSE
!
!        general case
!
         BETA = -SIGN( SLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = SLAMCH( 'S' ) / SLAMCH( 'E' )
         IF( ABS( BETA ).LT.SAFMIN ) THEN
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
            RSAFMN = ONE / SAFMIN
            KNT = 0
   10       CONTINUE
            KNT = KNT + 1
            CALL SSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN ) &
               GO TO 10
!
!           New BETA is at most 1, at least SAFMIN
!
            XNORM = SNRM2( N-1, X, INCX )
            BETA = -SIGN( SLAPY2( ALPHA, XNORM ), ALPHA )
            TAU = ( BETA-ALPHA ) / BETA
            CALL SSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
!
!           If ALPHA is subnormal, it may lose relative accuracy
!
            ALPHA = BETA
            DO 20 J = 1, KNT
               ALPHA = ALPHA*SAFMIN
   20       CONTINUE
         ELSE
            TAU = ( BETA-ALPHA ) / BETA
            CALL SSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
            ALPHA = BETA
         END IF
      END IF
!
      RETURN
!
!     End of SLARFG
!
      END SUBROUTINE SLARFG

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
!     ..
!     .. Array Arguments ..
      REAL(r8)               T( LDT, * ), TAU( * ), V( LDV, * )
!     ..
!
!  Purpose
!  =======
!
!  SLARFT forms the triangular factor T of a real block reflector H
!  of order n, which is defined as a product of k elementary reflectors.
!
!  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
!
!  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
!
!  If STOREV = 'C', the vector which defines the elementary reflector
!  H(i) is stored in the i-th column of the array V, and
!
!     H  =  I - V * T * V'
!
!  If STOREV = 'R', the vector which defines the elementary reflector
!  H(i) is stored in the i-th row of the array V, and
!
!     H  =  I - V' * T * V
!
!  Arguments
!  =========
!
!  DIRECT  (input) CHARACTER*1
!          Specifies the order in which the elementary reflectors are
!          multiplied to form the block reflector:
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
!  STOREV  (input) CHARACTER*1
!          Specifies how the vectors which define the elementary
!          reflectors are stored (see also Further Details):
!          = 'C': columnwise
!          = 'R': rowwise
!
!  N       (input) INTEGER
!          The order of the block reflector H. N >= 0.
!
!  K       (input) INTEGER
!          The order of the triangular factor T (= the number of
!          elementary reflectors). K >= 1.
!
!  V       (input/output) REAL array, dimension
!                               (LDV,K) if STOREV = 'C'
!                               (LDV,N) if STOREV = 'R'
!          The matrix V. See further details.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
!
!  TAU     (input) REAL array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i).
!
!  T       (output) REAL array, dimension (LDT,K)
!          The k by k triangular factor T of the block reflector.
!          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
!          lower triangular. The rest of the array is not used.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= K.
!
!  Further Details
!  ===============
!
!  The shape of the matrix V and the storage of the vectors which define
!  the H(i) is best illustrated by the following example with n = 5 and
!  k = 3. The elements equal to 1 are not stored; the corresponding
!  array elements are modified but restored on exit. The rest of the
!  array is not used.
!
!  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!
!               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!                   ( v1  1    )                     (     1 v2 v2 v2 )
!                   ( v1 v2  1 )                     (        1 v3 v3 )
!                   ( v1 v2 v3 )
!                   ( v1 v2 v3 )
!
!  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!
!               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!                   (     1 v3 )
!                   (        1 )
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0_r8, ZERO = 0.0E+0_r8 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      REAL(r8)               VII
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
      IF( LSAME( DIRECT, 'F' ) ) THEN
         DO 20 I = 1, K
            IF( TAU( I ).EQ.ZERO ) THEN
!
!              H(i)  =  I
!
               DO 10 J = 1, I
                  T( J, I ) = ZERO
   10          CONTINUE
            ELSE
!
!              general case
!
               VII = V( I, I )
               V( I, I ) = ONE
               IF( LSAME( STOREV, 'C' ) ) THEN
!
!                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i)
!
                  CALL SGEMV( 'Transpose', N-I+1, I-1, -TAU( I ), &
                              V( I, 1 ), LDV, V( I, I ), 1, ZERO, &
                              T( 1, I ), 1 )
               ELSE
!
!                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)'
!
                  CALL SGEMV( 'No transpose', I-1, N-I+1, -TAU( I ), &
                              V( 1, I ), LDV, V( I, I ), LDV, ZERO, &
                              T( 1, I ), 1 )
               END IF
               V( I, I ) = VII
!
!              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
!
               CALL STRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, &
                           LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
            END IF
   20    CONTINUE
      ELSE
         DO 40 I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
!
!              H(i)  =  I
!
               DO 30 J = I, K
                  T( J, I ) = ZERO
   30          CONTINUE
            ELSE
!
!              general case
!
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
                     VII = V( N-K+I, I )
                     V( N-K+I, I ) = ONE
!
!                    T(i+1:k,i) :=
!                            - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i)
!
                     CALL SGEMV( 'Transpose', N-K+I, K-I, -TAU( I ), &
                                 V( 1, I+1 ), LDV, V( 1, I ), 1, ZERO, &
                                 T( I+1, I ), 1 )
                     V( N-K+I, I ) = VII
                  ELSE
                     VII = V( I, N-K+I )
                     V( I, N-K+I ) = ONE
!
!                    T(i+1:k,i) :=
!                            - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)'
!
                     CALL SGEMV( 'No transpose', K-I, N-K+I, -TAU( I ), &
                                 V( I+1, 1 ), LDV, V( I, 1 ), LDV, ZERO, &
                                 T( I+1, I ), 1 )
                     V( I, N-K+I ) = VII
                  END IF
!
!                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
!
                  CALL STRMV( 'Lower', 'No transpose', 'Non-unit', K-I, &
                              T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
               END IF
               T( I, I ) = TAU( I )
            END IF
   40    CONTINUE
      END IF
      RETURN
!
!     End of SLARFT
!
      END SUBROUTINE SLARFT

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SLARFX( SIDE, M, N, V, TAU, C, LDC, WORK )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            LDC, M, N
      REAL(r8)               TAU
!     ..
!     .. Array Arguments ..
      REAL(r8)               C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  SLARFX applies a real elementary reflector H to a real m by n
!  matrix C, from either the left or the right. H is represented in the
!  form
!
!        H = I - tau * v * v'
!
!  where tau is a real scalar and v is a real vector.
!
!  If tau = 0, then H is taken to be the unit matrix
!
!  This version uses inline code if H has order < 11.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': form  H * C
!          = 'R': form  C * H
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  V       (input) REAL array, dimension (M) if SIDE = 'L'
!                                     or (N) if SIDE = 'R'
!          The vector v in the representation of H.
!
!  TAU     (input) REAL
!          The value tau in the representation of H.
!
!  C       (input/output) REAL array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!          or C * H if SIDE = 'R'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDA >= (1,M).
!
!  WORK    (workspace) REAL array, dimension
!                      (N) if SIDE = 'L'
!                      or (M) if SIDE = 'R'
!          WORK is not referenced if H has order < 11.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0_r8, ONE = 1.0E+0_r8 )
!     ..
!     .. Local Scalars ..
      INTEGER            J
      REAL(r8)               SUM, T1, T10, T2, T3, T4, T5, T6, T7, T8, T9, &
                         V1, V10, V2, V3, V4, V5, V6, V7, V8, V9
!     ..
!     .. Executable Statements ..
!
      IF( TAU.EQ.ZERO ) &
         RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
!
!        Form  H * C, where H has order m.
!
         GO TO ( 10, 30, 50, 70, 90, 110, 130, 150, &
                 170, 190 )M
!
!        Code for general M
!
!        w := C'*v
!
         CALL SGEMV( 'Transpose', M, N, ONE, C, LDC, V, 1, ZERO, WORK, &
                     1 )
!
!        C := C - tau * v * w'
!
         CALL SGER( M, N, -TAU, V, 1, WORK, 1, C, LDC )
         GO TO 410
   10    CONTINUE
!
!        Special code for 1 x 1 Householder
!
         T1 = ONE - TAU*V( 1 )*V( 1 )
         DO 20 J = 1, N
            C( 1, J ) = T1*C( 1, J )
   20    CONTINUE
         GO TO 410
   30    CONTINUE
!
!        Special code for 2 x 2 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         DO 40 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
   40    CONTINUE
         GO TO 410
   50    CONTINUE
!
!        Special code for 3 x 3 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         DO 60 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
   60    CONTINUE
         GO TO 410
   70    CONTINUE
!
!        Special code for 4 x 4 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         DO 80 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
   80    CONTINUE
         GO TO 410
   90    CONTINUE
!
!        Special code for 5 x 5 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         DO 100 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
  100    CONTINUE
         GO TO 410
  110    CONTINUE
!
!        Special code for 6 x 6 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         DO 120 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
  120    CONTINUE
         GO TO 410
  130    CONTINUE
!
!        Special code for 7 x 7 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         DO 140 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + &
                  V7*C( 7, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
            C( 7, J ) = C( 7, J ) - SUM*T7
  140    CONTINUE
         GO TO 410
  150    CONTINUE
!
!        Special code for 8 x 8 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         DO 160 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + &
                  V7*C( 7, J ) + V8*C( 8, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
            C( 7, J ) = C( 7, J ) - SUM*T7
            C( 8, J ) = C( 8, J ) - SUM*T8
  160    CONTINUE
         GO TO 410
  170    CONTINUE
!
!        Special code for 9 x 9 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         V9 = V( 9 )
         T9 = TAU*V9
         DO 180 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + &
                  V7*C( 7, J ) + V8*C( 8, J ) + V9*C( 9, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
            C( 7, J ) = C( 7, J ) - SUM*T7
            C( 8, J ) = C( 8, J ) - SUM*T8
            C( 9, J ) = C( 9, J ) - SUM*T9
  180    CONTINUE
         GO TO 410
  190    CONTINUE
!
!        Special code for 10 x 10 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         V9 = V( 9 )
         T9 = TAU*V9
         V10 = V( 10 )
         T10 = TAU*V10
         DO 200 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + &
                  V7*C( 7, J ) + V8*C( 8, J ) + V9*C( 9, J ) + &
                  V10*C( 10, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
            C( 7, J ) = C( 7, J ) - SUM*T7
            C( 8, J ) = C( 8, J ) - SUM*T8
            C( 9, J ) = C( 9, J ) - SUM*T9
            C( 10, J ) = C( 10, J ) - SUM*T10
  200    CONTINUE
         GO TO 410
      ELSE
!
!        Form  C * H, where H has order n.
!
         GO TO ( 210, 230, 250, 270, 290, 310, 330, 350, &
                 370, 390 )N
!
!        Code for general N
!
!        w := C * v
!
         CALL SGEMV( 'No transpose', M, N, ONE, C, LDC, V, 1, ZERO, &
                     WORK, 1 )
!
!        C := C - tau * w * v'
!
         CALL SGER( M, N, -TAU, WORK, 1, V, 1, C, LDC )
         GO TO 410
  210    CONTINUE
!
!        Special code for 1 x 1 Householder
!
         T1 = ONE - TAU*V( 1 )*V( 1 )
         DO 220 J = 1, M
            C( J, 1 ) = T1*C( J, 1 )
  220    CONTINUE
         GO TO 410
  230    CONTINUE
!
!        Special code for 2 x 2 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         DO 240 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
  240    CONTINUE
         GO TO 410
  250    CONTINUE
!
!        Special code for 3 x 3 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         DO 260 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
  260    CONTINUE
         GO TO 410
  270    CONTINUE
!
!        Special code for 4 x 4 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         DO 280 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
  280    CONTINUE
         GO TO 410
  290    CONTINUE
!
!        Special code for 5 x 5 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         DO 300 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
  300    CONTINUE
         GO TO 410
  310    CONTINUE
!
!        Special code for 6 x 6 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         DO 320 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
  320    CONTINUE
         GO TO 410
  330    CONTINUE
!
!        Special code for 7 x 7 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         DO 340 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + &
                  V7*C( J, 7 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
            C( J, 7 ) = C( J, 7 ) - SUM*T7
  340    CONTINUE
         GO TO 410
  350    CONTINUE
!
!        Special code for 8 x 8 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         DO 360 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + &
                  V7*C( J, 7 ) + V8*C( J, 8 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
            C( J, 7 ) = C( J, 7 ) - SUM*T7
            C( J, 8 ) = C( J, 8 ) - SUM*T8
  360    CONTINUE
         GO TO 410
  370    CONTINUE
!
!        Special code for 9 x 9 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         V9 = V( 9 )
         T9 = TAU*V9
         DO 380 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + &
                  V7*C( J, 7 ) + V8*C( J, 8 ) + V9*C( J, 9 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
            C( J, 7 ) = C( J, 7 ) - SUM*T7
            C( J, 8 ) = C( J, 8 ) - SUM*T8
            C( J, 9 ) = C( J, 9 ) - SUM*T9
  380    CONTINUE
         GO TO 410
  390    CONTINUE
!
!        Special code for 10 x 10 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         V9 = V( 9 )
         T9 = TAU*V9
         V10 = V( 10 )
         T10 = TAU*V10
         DO 400 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + &
                  V7*C( J, 7 ) + V8*C( J, 8 ) + V9*C( J, 9 ) + &
                  V10*C( J, 10 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
            C( J, 7 ) = C( J, 7 ) - SUM*T7
            C( J, 8 ) = C( J, 8 ) - SUM*T8
            C( J, 9 ) = C( J, 9 ) - SUM*T9
            C( J, 10 ) = C( J, 10 ) - SUM*T10
  400    CONTINUE
         GO TO 410
      END IF
  410 RETURN
!
!     End of SLARFX
!
      END SUBROUTINE SLARFX

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SLARTG( F, G, CS, SN, R )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      REAL(r8)               CS, F, G, R, SN
!     ..
!
!  Purpose
!  =======
!
!  SLARTG generate a plane rotation so that
!
!     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
!     [ -SN  CS  ]     [ G ]     [ 0 ]
!
!  This is a slower, more accurate version of the BLAS1 routine SROTG,
!  with the following other differences:
!     F and G are unchanged on return.
!     If G=0, then CS=1 and SN=0.
!     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
!        floating point operations (saves work in SBDSQR when
!        there are zeros on the diagonal).
!
!  If F exceeds G in magnitude, CS will be positive.
!
!  Arguments
!  =========
!
!  F       (input) REAL
!          The first component of vector to be rotated.
!
!  G       (input) REAL
!          The second component of vector to be rotated.
!
!  CS      (output) REAL
!          The cosine of the rotation.
!
!  SN      (output) REAL
!          The sine of the rotation.
!
!  R       (output) REAL
!          The nonzero component of the rotated vector.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ZERO
      PARAMETER          ( ZERO = 0.0E0_r8 )
      REAL(r8)               ONE
      PARAMETER          ( ONE = 1.0E0_r8 )
      REAL(r8)               TWO
      PARAMETER          ( TWO = 2.0E0_r8 )
!     ..
!     .. Local Scalars ..
      LOGICAL            FIRST
      INTEGER            COUNT, I
      REAL(r8)               EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
!     ..
!     .. Save statement ..
      SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         SAFMIN = SLAMCH( 'S' )
         EPS = SLAMCH( 'E' )
         SAFMN2 = SLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) / &
                  LOG( SLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
      END IF
      IF( G.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF( F.EQ.ZERO ) THEN
         CS = ZERO
         SN = ONE
         R = G
      ELSE
         F1 = F
         G1 = G
         SCALE = MAX( ABS( F1 ), ABS( G1 ) )
         IF( SCALE.GE.SAFMX2 ) THEN
            COUNT = 0
   10       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.GE.SAFMX2 ) &
               GO TO 10
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 20 I = 1, COUNT
               R = R*SAFMX2
   20       CONTINUE
         ELSE IF( SCALE.LE.SAFMN2 ) THEN
            COUNT = 0
   30       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.LE.SAFMN2 ) &
               GO TO 30
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 40 I = 1, COUNT
               R = R*SAFMN2
   40       CONTINUE
         ELSE
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
         END IF
         IF( ABS( F ).GT.ABS( G ) .AND. CS.LT.ZERO ) THEN
            CS = -CS
            SN = -SN
            R = -R
         END IF
      END IF
      RETURN
!
!     End of SLARTG
!
      END SUBROUTINE SLARTG

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      REAL(r8)               CFROM, CTO
!     ..
!     .. Array Arguments ..
      REAL(r8)               A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  SLASCL multiplies the M by N real matrix A by the real scalar
!  CTO/CFROM.  This is done without over/underflow as long as the final
!  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
!  A may be full, upper triangular, lower triangular, upper Hessenberg,
!  or banded.
!
!  Arguments
!  =========
!
!  TYPE    (input) CHARACTER*1
!          TYPE indices the storage type of the input matrix.
!          = 'G':  A is a full matrix.
!          = 'L':  A is a lower triangular matrix.
!          = 'U':  A is an upper triangular matrix.
!          = 'H':  A is an upper Hessenberg matrix.
!          = 'B':  A is a symmetric band matrix with lower bandwidth KL
!                  and upper bandwidth KU and with the only the lower
!                  half stored.
!          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
!                  and upper bandwidth KU and with the only the upper
!                  half stored.
!          = 'Z':  A is a band matrix with lower bandwidth KL and upper
!                  bandwidth KU.
!
!  KL      (input) INTEGER
!          The lower bandwidth of A.  Referenced only if TYPE = 'B',
!          'Q' or 'Z'.
!
!  KU      (input) INTEGER
!          The upper bandwidth of A.  Referenced only if TYPE = 'B',
!          'Q' or 'Z'.
!
!  CFROM   (input) REAL
!  CTO     (input) REAL
!          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
!          without over/underflow if the final result CTO*A(I,J)/CFROM
!          can be represented without over/underflow.  CFROM must be
!          nonzero.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) REAL array, dimension (LDA,M)
!          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
!          storage type.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  INFO    (output) INTEGER
!          0  - successful exit
!          <0 - if INFO = -i, the i-th argument had an illegal value.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0_r8, ONE = 1.0E0_r8 )
!     ..
!     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      REAL(r8)               BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
!
      IF( LSAME( TYPE, 'G' ) ) THEN
         ITYPE = 0
      ELSE IF( LSAME( TYPE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
         ITYPE = 3
      ELSE IF( LSAME( TYPE, 'B' ) ) THEN
         ITYPE = 4
      ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
         ITYPE = 5
      ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
         ITYPE = 6
      ELSE
         ITYPE = -1
      END IF
!
      IF( ITYPE.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( CFROM.EQ.ZERO ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR. &
               ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR. &
                  ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) ) &
                   THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR. &
                  ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR. &
                  ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLASCL', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. M.EQ.0 ) &
         RETURN
!
!     Get machine parameters
!
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
!
      CFROMC = CFROM
      CTOC = CTO
!
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      CTO1 = CTOC / BIGNUM
      IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
         MUL = SMLNUM
         DONE = .FALSE.
         CFROMC = CFROM1
      ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
         MUL = BIGNUM
         DONE = .FALSE.
         CTOC = CTO1
      ELSE
         MUL = CTOC / CFROMC
         DONE = .TRUE.
      END IF
!
      IF( ITYPE.EQ.0 ) THEN
!
!        Full matrix
!
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
!
      ELSE IF( ITYPE.EQ.1 ) THEN
!
!        Lower triangular matrix
!
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
!
      ELSE IF( ITYPE.EQ.2 ) THEN
!
!        Upper triangular matrix
!
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
!
      ELSE IF( ITYPE.EQ.3 ) THEN
!
!        Upper Hessenberg matrix
!
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
!
      ELSE IF( ITYPE.EQ.4 ) THEN
!
!        Lower half of a symmetric band matrix
!
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
!
      ELSE IF( ITYPE.EQ.5 ) THEN
!
!        Upper half of a symmetric band matrix
!
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
!
      ELSE IF( ITYPE.EQ.6 ) THEN
!
!        Band matrix
!
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
!
      END IF
!
      IF( .NOT.DONE ) &
         GO TO 10
!
      RETURN
!
!     End of SLASCL
!
      END SUBROUTINE SLASCL

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SLASET( UPLO, M, N, ALPHA, BETA, A, LDA )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      REAL(r8)               ALPHA, BETA
!     ..
!     .. Array Arguments ..
      REAL(r8)               A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  SLASET initializes an m-by-n matrix A to BETA on the diagonal and
!  ALPHA on the offdiagonals.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies the part of the matrix A to be set.
!          = 'U':      Upper triangular part is set; the strictly lower
!                      triangular part of A is not changed.
!          = 'L':      Lower triangular part is set; the strictly upper
!                      triangular part of A is not changed.
!          Otherwise:  All of the matrix A is set.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  ALPHA   (input) REAL
!          The constant to which the offdiagonal elements are to be set.
!
!  BETA    (input) REAL
!          The constant to which the diagonal elements are to be set.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On exit, the leading m-by-n submatrix of A is set as follows:
!
!          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
!          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
!          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
!
!          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Set the strictly upper triangular or trapezoidal part of the
!        array to ALPHA.
!
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
!
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
!
!        Set the strictly lower triangular or trapezoidal part of the
!        array to ALPHA.
!
         DO 40 J = 1, MIN( M, N )
            DO 30 I = J + 1, M
               A( I, J ) = ALPHA
   30       CONTINUE
   40    CONTINUE
!
      ELSE
!
!        Set the leading m-by-n submatrix to ALPHA.
!
         DO 60 J = 1, N
            DO 50 I = 1, M
               A( I, J ) = ALPHA
   50       CONTINUE
   60    CONTINUE
      END IF
!
!     Set the first min(M,N) diagonal elements to BETA.
!
      DO 70 I = 1, MIN( M, N )
         A( I, I ) = BETA
   70 CONTINUE
!
      RETURN
!
!     End of SLASET
!
      END SUBROUTINE SLASET

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SLASSQ( N, X, INCX, SCALE, SUMSQ )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      REAL(r8)               SCALE, SUMSQ
!     ..
!     .. Array Arguments ..
      REAL(r8)               X( * )
!     ..
!
!  Purpose
!  =======
!
!  SLASSQ  returns the values  scl  and  smsq  such that
!
!     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
!
!  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
!  assumed to be non-negative and  scl  returns the value
!
!     scl = max( scale, abs( x( i ) ) ).
!
!  scale and sumsq must be supplied in SCALE and SUMSQ and
!  scl and smsq are overwritten on SCALE and SUMSQ respectively.
!
!  The routine makes only one pass through the vector x.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements to be used from the vector X.
!
!  X       (input) REAL
!          The vector for which a scaled sum of squares is computed.
!             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
!
!  INCX    (input) INTEGER
!          The increment between successive values of the vector X.
!          INCX > 0.
!
!  SCALE   (input/output) REAL
!          On entry, the value  scale  in the equation above.
!          On exit, SCALE is overwritten with  scl , the scaling factor
!          for the sum of squares.
!
!  SUMSQ   (input/output) REAL
!          On entry, the value  sumsq  in the equation above.
!          On exit, SUMSQ is overwritten with  smsq , the basic sum of
!          squares from which  scl  has been factored out.
!
! =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ZERO
      PARAMETER          ( ZERO = 0.0E+0_r8 )
!     ..
!     .. Local Scalars ..
      INTEGER            IX
      REAL(r8)               ABSXI
!     ..
!     .. Executable Statements ..
!
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( X( IX ).NE.ZERO ) THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
!
!     End of SLASSQ
!
      END SUBROUTINE SLASSQ

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SORG2R( M, N, K, A, LDA, TAU, WORK, INFO )

!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
!     ..
!     .. Array Arguments ..
      REAL(r8)               A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  SORG2R generates an m by n real matrix Q with orthonormal columns,
!  which is defined as the first n columns of a product of k elementary
!  reflectors of order m
!
!        Q  =  H(1) H(2) . . . H(k)
!
!  as returned by SGEQRF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the i-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by SGEQRF in the first k columns of its array
!          argument A.
!          On exit, the m-by-n matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) REAL array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by SGEQRF.
!
!  WORK    (workspace) REAL array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0_r8, ZERO = 0.0E+0_r8 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, L
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SORG2R', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) &
         RETURN
!
!     Initialise columns k+1:n to columns of the unit matrix
!
      DO 20 J = K + 1, N
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( J, J ) = ONE
   20 CONTINUE
!
      DO 40 I = K, 1, -1
!
!        Apply H(i) to A(i:m,i:n) from the left
!
         IF( I.LT.N ) THEN
            A( I, I ) = ONE
            CALL SLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), &
                        A( I, I+1 ), LDA, WORK )
         END IF
         IF( I.LT.M ) &
            CALL SSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
         A( I, I ) = ONE - TAU( I )
!
!        Set A(1:i-1,i) to zero
!
         DO 30 L = 1, I - 1
            A( L, I ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
!
!     End of SORG2R
!
      END SUBROUTINE SORG2R

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SORGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )

!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      REAL(r8)               A( LDA, * ), TAU( * ), WORK( LWORK )
!     ..
!
!  Purpose
!  =======
!
!  SORGHR generates a real orthogonal matrix Q which is defined as the
!  product of IHI-ILO elementary reflectors of order N, as returned by
!  SGEHRD:
!
!  Q = H(ilo) H(ilo+1) . . . H(ihi-1).
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix Q. N >= 0.
!
!  ILO     (input) INTEGER
!  IHI     (input) INTEGER
!          ILO and IHI must have the same values as in the previous call
!          of SGEHRD. Q is equal to the unit matrix except in the
!          submatrix Q(ilo+1:ihi,ilo+1:ihi).
!          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the vectors which define the elementary reflectors,
!          as returned by SGEHRD.
!          On exit, the N-by-N orthogonal matrix Q.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,N).
!
!  TAU     (input) REAL array, dimension (N-1)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by SGEHRD.
!
!  WORK    (workspace/output) REAL array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= IHI-ILO.
!          For optimum performance LWORK >= (IHI-ILO)*NB, where NB is
!          the optimal blocksize.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0_r8, ONE = 1.0E+0_r8 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IINFO, J, NH
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, IHI-ILO ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SORGHR', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
!     Shift the vectors which define the elementary reflectors one
!     column to the right, and set the first ilo and the last n-ihi
!     rows and columns to those of the unit matrix
!
      DO 40 J = IHI, ILO + 1, -1
         DO 10 I = 1, J - 1
            A( I, J ) = ZERO
   10    CONTINUE
         DO 20 I = J + 1, IHI
            A( I, J ) = A( I, J-1 )
   20    CONTINUE
         DO 30 I = IHI + 1, N
            A( I, J ) = ZERO
   30    CONTINUE
   40 CONTINUE
      DO 60 J = 1, ILO
         DO 50 I = 1, N
            A( I, J ) = ZERO
   50    CONTINUE
         A( J, J ) = ONE
   60 CONTINUE
      DO 80 J = IHI + 1, N
         DO 70 I = 1, N
            A( I, J ) = ZERO
   70    CONTINUE
         A( J, J ) = ONE
   80 CONTINUE
!
      NH = IHI - ILO
      IF( NH.GT.0 ) THEN
!
!        Generate Q(ilo+1:ihi,ilo+1:ihi)
!
         CALL SORGQR( NH, NH, NH, A( ILO+1, ILO+1 ), LDA, TAU( ILO ), &
                      WORK, LWORK, IINFO )
      END IF
      RETURN
!
!     End of SORGHR
!
      END SUBROUTINE SORGHR

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )

!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      REAL(r8)               A( LDA, * ), TAU( * ), WORK( LWORK )
!     ..
!
!  Purpose
!  =======
!
!  SORGQR generates an M-by-N real matrix Q with orthonormal columns,
!  which is defined as the first N columns of a product of K elementary
!  reflectors of order M
!
!        Q  =  H(1) H(2) . . . H(k)
!
!  as returned by SGEQRF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the i-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by SGEQRF in the first k columns of its array
!          argument A.
!          On exit, the M-by-N matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) REAL array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by SGEQRF.
!
!  WORK    (workspace/output) REAL array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= max(1,N).
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ZERO
      PARAMETER          ( ZERO = 0.0E+0_r8 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, NB, &
                         NBMIN, NX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SORGQR', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
!     Determine the block size.
!
      NB = ILAENV( 1, 'SORGQR', ' ', M, N, K, -1 )
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
         NX = MAX( 0, ILAENV( 3, 'SORGQR', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'SORGQR', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
!
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!
!        Use blocked code after the last block.
!        The first kk columns are handled by the block method.
!
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
!
!        Set A(1:kk,kk+1:n) to zero.
!
         DO 20 J = KK + 1, N
            DO 10 I = 1, KK
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
!
!     Use unblocked code for the last or only block.
!
      IF( KK.LT.N ) &
         CALL SORG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA, &
                      TAU( KK+1 ), WORK, IINFO )
!
      IF( KK.GT.0 ) THEN
!
!        Use blocked code
!
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.N ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
               CALL SLARFT( 'Forward', 'Columnwise', M-I+1, IB, &
                            A( I, I ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H to A(i:m,i+ib:n) from the left
!
               CALL SLARFB( 'Left', 'No transpose', 'Forward', &
                            'Columnwise', M-I+1, N-I-IB+1, IB, &
                            A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ), &
                            LDA, WORK( IB+1 ), LDWORK )
            END IF
!
!           Apply H to rows i:m of current block
!
            CALL SORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK, &
                         IINFO )
!
!           Set rows 1:i-1 of current block to zero
!
            DO 40 J = I, I + IB - 1
               DO 30 L = 1, I - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of SORGQR
!
      END SUBROUTINE SORGQR

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE STREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, &
                         LDVR, MM, M, WORK, INFO )

!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          HOWMNY, SIDE
      INTEGER            INFO, LDT, LDVL, LDVR, M, MM, N
!     ..
!     .. Array Arguments ..
      LOGICAL            SELECT( * )
      REAL(r8)               T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), &
                         WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  STREVC computes some or all of the right and/or left eigenvectors of
!  a real upper quasi-triangular matrix T.
!
!  The right eigenvector x and the left eigenvector y of T corresponding
!  to an eigenvalue w are defined by:
!
!               T*x = w*x,     y'*T = w*y'
!
!  where y' denotes the conjugate transpose of the vector y.
!
!  If all eigenvectors are requested, the routine may either return the
!  matrices X and/or Y of right or left eigenvectors of T, or the
!  products Q*X and/or Q*Y, where Q is an input orthogonal
!  matrix. If T was obtained from the real-Schur factorization of an
!  original matrix A = Q*T*Q', then Q*X and Q*Y are the matrices of
!  right or left eigenvectors of A.
!
!  T must be in Schur canonical form (as returned by SHSEQR), that is,
!  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
!  2-by-2 diagonal block has its diagonal elements equal and its
!  off-diagonal elements of opposite sign.  Corresponding to each 2-by-2
!  diagonal block is a complex conjugate pair of eigenvalues and
!  eigenvectors; only one eigenvector of the pair is computed, namely
!  the one corresponding to the eigenvalue with positive imaginary part.

!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'R':  compute right eigenvectors only;
!          = 'L':  compute left eigenvectors only;
!          = 'B':  compute both right and left eigenvectors.
!
!  HOWMNY  (input) CHARACTER*1
!          = 'A':  compute all right and/or left eigenvectors;
!          = 'B':  compute all right and/or left eigenvectors,
!                  and backtransform them using the input matrices
!                  supplied in VR and/or VL;
!          = 'S':  compute selected right and/or left eigenvectors,
!                  specified by the logical array SELECT.
!
!  SELECT  (input/output) LOGICAL array, dimension (N)
!          If HOWMNY = 'S', SELECT specifies the eigenvectors to be
!          computed.
!          If HOWMNY = 'A' or 'B', SELECT is not referenced.
!          To select the real eigenvector corresponding to a real
!          eigenvalue w(j), SELECT(j) must be set to .TRUE..  To select
!          the complex eigenvector corresponding to a complex conjugate
!          pair w(j) and w(j+1), either SELECT(j) or SELECT(j+1) must be
!          set to .TRUE.; then on exit SELECT(j) is .TRUE. and
!          SELECT(j+1) is .FALSE..
!
!  N       (input) INTEGER
!          The order of the matrix T. N >= 0.
!
!  T       (input) REAL array, dimension (LDT,N)
!          The upper quasi-triangular matrix T in Schur canonical form.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= max(1,N).
!
!  VL      (input/output) REAL array, dimension (LDVL,MM)
!          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
!          contain an N-by-N matrix Q (usually the orthogonal matrix Q
!          of Schur vectors returned by SHSEQR).
!          On exit, if SIDE = 'L' or 'B', VL contains:
!          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
!          if HOWMNY = 'B', the matrix Q*Y;
!          if HOWMNY = 'S', the left eigenvectors of T specified by
!                           SELECT, stored consecutively in the columns
!                           of VL, in the same order as their
!                           eigenvalues.
!          A complex eigenvector corresponding to a complex eigenvalue
!          is stored in two consecutive columns, the first holding the
!          real part, and the second the imaginary part.
!          If SIDE = 'R', VL is not referenced.
!
!  LDVL    (input) INTEGER
!          The leading dimension of the array VL.  LDVL >= max(1,N) if
!          SIDE = 'L' or 'B'; LDVL >= 1 otherwise.
!
!  VR      (input/output) REAL array, dimension (LDVR,MM)
!          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
!          contain an N-by-N matrix Q (usually the orthogonal matrix Q
!          of Schur vectors returned by SHSEQR).
!          On exit, if SIDE = 'R' or 'B', VR contains:
!          if HOWMNY = 'A', the matrix X of right eigenvectors of T;
!          if HOWMNY = 'B', the matrix Q*X;
!          if HOWMNY = 'S', the right eigenvectors of T specified by
!                           SELECT, stored consecutively in the columns
!                           of VR, in the same order as their
!                           eigenvalues.
!          A complex eigenvector corresponding to a complex eigenvalue
!          is stored in two consecutive columns, the first holding the
!          real part and the second the imaginary part.
!          If SIDE = 'L', VR is not referenced.
!
!  LDVR    (input) INTEGER
!          The leading dimension of the array VR.  LDVR >= max(1,N) if
!          SIDE = 'R' or 'B'; LDVR >= 1 otherwise.
!
!  MM      (input) INTEGER
!          The number of columns in the arrays VL and/or VR. MM >= M.
!
!  M       (output) INTEGER
!          The number of columns in the arrays VL and/or VR actually
!          used to store the eigenvectors.
!          If HOWMNY = 'A' or 'B', M is set to N.
!          Each selected real eigenvector occupies one column and each
!          selected complex eigenvector occupies two columns.
!
!  WORK    (workspace) REAL array, dimension (3*N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  The algorithm used in this program is basically backward (forward)
!  substitution, with scaling to make the the code robust against
!  possible overflow.
!
!  Each eigenvector is normalized so that the element of largest
!  magnitude has magnitude 1; here the magnitude of a complex number
!  (x,y) is taken to be |x| + |y|.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0_r8, ONE = 1.0E+0_r8 )
!     ..
!     .. Local Scalars ..
      LOGICAL            ALLV, BOTHV, LEFTV, OVER, PAIR, RIGHTV, SOMEV
      INTEGER            I, IERR, II, IP, IS, J, J1, J2, JNXT, K, KI, N2
      REAL(r8)               BETA, BIGNUM, EMAX, OVFL, REC, REMAX, SCALE, &
                         SMIN, SMLNUM, ULP, UNFL, VCRIT, VMAX, WI, WR, &
                         XNORM
!     ..
!     .. Local Arrays ..
      REAL(r8)               X( 2, 2 )
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters
!
      BOTHV = LSAME( SIDE, 'B' )
      RIGHTV = LSAME( SIDE, 'R' ) .OR. BOTHV
      LEFTV = LSAME( SIDE, 'L' ) .OR. BOTHV
!
      ALLV = LSAME( HOWMNY, 'A' )
      OVER = LSAME( HOWMNY, 'B' ) .OR. LSAME( HOWMNY, 'O' )
      SOMEV = LSAME( HOWMNY, 'S' )
!
      INFO = 0
      IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
         INFO = -1
      ELSE IF( .NOT.ALLV .AND. .NOT.OVER .AND. .NOT.SOMEV ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDVL.LT.1 .OR. ( LEFTV .AND. LDVL.LT.N ) ) THEN
         INFO = -8
      ELSE IF( LDVR.LT.1 .OR. ( RIGHTV .AND. LDVR.LT.N ) ) THEN
         INFO = -10
      ELSE
!
!        Set M to the number of columns required to store the selected
!        eigenvectors, standardize the array SELECT if necessary, and
!        test MM.
!
         IF( SOMEV ) THEN
            M = 0
            PAIR = .FALSE.
            DO 10 J = 1, N
               IF( PAIR ) THEN
                  PAIR = .FALSE.
                  SELECT( J ) = .FALSE.
               ELSE
                  IF( J.LT.N ) THEN
                     IF( T( J+1, J ).EQ.ZERO ) THEN
                        IF( SELECT( J ) ) &
                           M = M + 1
                     ELSE
                        PAIR = .TRUE.
                        IF( SELECT( J ) .OR. SELECT( J+1 ) ) THEN
                           SELECT( J ) = .TRUE.
                           M = M + 2
                        END IF
                     END IF
                  ELSE
                     IF( SELECT( N ) ) &
                        M = M + 1
                  END IF
               END IF
   10       CONTINUE
         ELSE
            M = N
         END IF
!
         IF( MM.LT.M ) THEN
            INFO = -11
         END IF
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'STREVC', -INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 ) &
         RETURN
!
!     Set the constants to control overflow.
!
      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      CALL SLABAD( UNFL, OVFL )
      ULP = SLAMCH( 'Precision' )
      SMLNUM = UNFL*( N / ULP )
      BIGNUM = ( ONE-ULP ) / SMLNUM
!
!     Compute 1-norm of each column of strictly upper triangular
!     part of T to control overflow in triangular solver.
!
      WORK( 1 ) = ZERO
      DO 30 J = 2, N
         WORK( J ) = ZERO
         DO 20 I = 1, J - 1
            WORK( J ) = WORK( J ) + ABS( T( I, J ) )
   20    CONTINUE
   30 CONTINUE
!
!     Index IP is used to specify the real or complex eigenvalue:
!       IP = 0, real eigenvalue,
!            1, first of conjugate complex pair: (wr,wi)
!           -1, second of conjugate complex pair: (wr,wi)
!
      N2 = 2*N
!
      IF( RIGHTV ) THEN
!
!        Compute right eigenvectors.
!
         IP = 0
         IS = M
         DO 140 KI = N, 1, -1
!
            IF( IP.EQ.1 ) &
               GO TO 130
            IF( KI.EQ.1 ) &
               GO TO 40
            IF( T( KI, KI-1 ).EQ.ZERO ) &
               GO TO 40
            IP = -1
!
   40       CONTINUE
            IF( SOMEV ) THEN
               IF( IP.EQ.0 ) THEN
                  IF( .NOT.SELECT( KI ) ) &
                     GO TO 130
               ELSE
                  IF( .NOT.SELECT( KI-1 ) ) &
                     GO TO 130
               END IF
            END IF
!
!           Compute the KI-th eigenvalue (WR,WI).
!
            WR = T( KI, KI )
            WI = ZERO
            IF( IP.NE.0 ) &
               WI = SQRT( ABS( T( KI, KI-1 ) ) )* &
                    SQRT( ABS( T( KI-1, KI ) ) )
            SMIN = MAX( ULP*( ABS( WR )+ABS( WI ) ), SMLNUM )
!
            IF( IP.EQ.0 ) THEN
!
!              Real right eigenvector
!
               WORK( KI+N ) = ONE
!
!              Form right-hand side
!
               DO 50 K = 1, KI - 1
                  WORK( K+N ) = -T( K, KI )
   50          CONTINUE
!
!              Solve the upper quasi-triangular system:
!                 (T(1:KI-1,1:KI-1) - WR)*X = SCALE*WORK.
!
               JNXT = KI - 1
               DO 60 J = KI - 1, 1, -1
                  IF( J.GT.JNXT ) &
                     GO TO 60
                  J1 = J
                  J2 = J
                  JNXT = J - 1
                  IF( J.GT.1 ) THEN
                     IF( T( J, J-1 ).NE.ZERO ) THEN
                        J1 = J - 1
                        JNXT = J - 2
                     END IF
                  END IF
!
                  IF( J1.EQ.J2 ) THEN
!
!                    1-by-1 diagonal block
!
                     CALL SLALN2( .FALSE., 1, 1, SMIN, ONE, T( J, J ), &
                                  LDT, ONE, ONE, WORK( J+N ), N, WR, &
                                  ZERO, X, 2, SCALE, XNORM, IERR )
!
!                    Scale X(1,1) to avoid overflow when updating
!                    the right-hand side.
!
                     IF( XNORM.GT.ONE ) THEN
                        IF( WORK( J ).GT.BIGNUM / XNORM ) THEN
                           X( 1, 1 ) = X( 1, 1 ) / XNORM
                           SCALE = SCALE / XNORM
                        END IF
                     END IF
!
!                    Scale if necessary
!
                     IF( SCALE.NE.ONE ) &
                        CALL SSCAL( KI, SCALE, WORK( 1+N ), 1 )
                     WORK( J+N ) = X( 1, 1 )
!
!                    Update right-hand side
!
                     CALL SAXPY( J-1, -X( 1, 1 ), T( 1, J ), 1, &
                                 WORK( 1+N ), 1 )
!
                  ELSE
!
!                    2-by-2 diagonal block
!
                     CALL SLALN2( .FALSE., 2, 1, SMIN, ONE, &
                                  T( J-1, J-1 ), LDT, ONE, ONE, &
                                  WORK( J-1+N ), N, WR, ZERO, X, 2, &
                                  SCALE, XNORM, IERR )
!
!                    Scale X(1,1) and X(2,1) to avoid overflow when
!                    updating the right-hand side.
!
                     IF( XNORM.GT.ONE ) THEN
                        BETA = MAX( WORK( J-1 ), WORK( J ) )
                        IF( BETA.GT.BIGNUM / XNORM ) THEN
                           X( 1, 1 ) = X( 1, 1 ) / XNORM
                           X( 2, 1 ) = X( 2, 1 ) / XNORM
                           SCALE = SCALE / XNORM
                        END IF
                     END IF
!
!                    Scale if necessary
!
                     IF( SCALE.NE.ONE ) &
                        CALL SSCAL( KI, SCALE, WORK( 1+N ), 1 )
                     WORK( J-1+N ) = X( 1, 1 )
                     WORK( J+N ) = X( 2, 1 )
!
!                    Update right-hand side
!
                     CALL SAXPY( J-2, -X( 1, 1 ), T( 1, J-1 ), 1, &
                                 WORK( 1+N ), 1 )
                     CALL SAXPY( J-2, -X( 2, 1 ), T( 1, J ), 1, &
                                 WORK( 1+N ), 1 )
                  END IF
   60          CONTINUE
!
!              Copy the vector x or Q*x to VR and normalize.
!
               IF( .NOT.OVER ) THEN
                  CALL SCOPY( KI, WORK( 1+N ), 1, VR( 1, IS ), 1 )
!
                  II = ISAMAX( KI, VR( 1, IS ), 1 )
                  REMAX = ONE / ABS( VR( II, IS ) )
                  CALL SSCAL( KI, REMAX, VR( 1, IS ), 1 )
!
                  DO 70 K = KI + 1, N
                     VR( K, IS ) = ZERO
   70             CONTINUE
               ELSE
                  IF( KI.GT.1 ) &
                     CALL SGEMV( 'N', N, KI-1, ONE, VR, LDVR, &
                                 WORK( 1+N ), 1, WORK( KI+N ), &
                                 VR( 1, KI ), 1 )
!
                  II = ISAMAX( N, VR( 1, KI ), 1 )
                  REMAX = ONE / ABS( VR( II, KI ) )
                  CALL SSCAL( N, REMAX, VR( 1, KI ), 1 )
               END IF
!
            ELSE
!
!              Complex right eigenvector.
!
!              Initial solve
!                [ (T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I* WI)]*X = 0.
!                [ (T(KI,KI-1)   T(KI,KI)   )               ]
!
               IF( ABS( T( KI-1, KI ) ).GE.ABS( T( KI, KI-1 ) ) ) THEN
                  WORK( KI-1+N ) = ONE
                  WORK( KI+N2 ) = WI / T( KI-1, KI )
               ELSE
                  WORK( KI-1+N ) = -WI / T( KI, KI-1 )
                  WORK( KI+N2 ) = ONE
               END IF
               WORK( KI+N ) = ZERO
               WORK( KI-1+N2 ) = ZERO
!
!              Form right-hand side
!
               DO 80 K = 1, KI - 2
                  WORK( K+N ) = -WORK( KI-1+N )*T( K, KI-1 )
                  WORK( K+N2 ) = -WORK( KI+N2 )*T( K, KI )
   80          CONTINUE
!
!              Solve upper quasi-triangular system:
!              (T(1:KI-2,1:KI-2) - (WR+i*WI))*X = SCALE*(WORK+i*WORK2)
!
               JNXT = KI - 2
               DO 90 J = KI - 2, 1, -1
                  IF( J.GT.JNXT ) &
                     GO TO 90
                  J1 = J
                  J2 = J
                  JNXT = J - 1
                  IF( J.GT.1 ) THEN
                     IF( T( J, J-1 ).NE.ZERO ) THEN
                        J1 = J - 1
                        JNXT = J - 2
                     END IF
                  END IF
!
                  IF( J1.EQ.J2 ) THEN
!
!                    1-by-1 diagonal block
!
                     CALL SLALN2( .FALSE., 1, 2, SMIN, ONE, T( J, J ), &
                                  LDT, ONE, ONE, WORK( J+N ), N, WR, WI, &
                                  X, 2, SCALE, XNORM, IERR )
!
!                    Scale X(1,1) and X(1,2) to avoid overflow when
!                    updating the right-hand side.
!
                     IF( XNORM.GT.ONE ) THEN
                        IF( WORK( J ).GT.BIGNUM / XNORM ) THEN
                           X( 1, 1 ) = X( 1, 1 ) / XNORM
                           X( 1, 2 ) = X( 1, 2 ) / XNORM
                           SCALE = SCALE / XNORM
                        END IF
                     END IF
!
!                    Scale if necessary
!
                     IF( SCALE.NE.ONE ) THEN
                        CALL SSCAL( KI, SCALE, WORK( 1+N ), 1 )
                        CALL SSCAL( KI, SCALE, WORK( 1+N2 ), 1 )
                     END IF
                     WORK( J+N ) = X( 1, 1 )
                     WORK( J+N2 ) = X( 1, 2 )
!
!                    Update the right-hand side
!
                     CALL SAXPY( J-1, -X( 1, 1 ), T( 1, J ), 1, &
                                 WORK( 1+N ), 1 )
                     CALL SAXPY( J-1, -X( 1, 2 ), T( 1, J ), 1, &
                                 WORK( 1+N2 ), 1 )
!
                  ELSE
!
!                    2-by-2 diagonal block
!
                     CALL SLALN2( .FALSE., 2, 2, SMIN, ONE, &
                                  T( J-1, J-1 ), LDT, ONE, ONE, &
                                  WORK( J-1+N ), N, WR, WI, X, 2, SCALE, &
                                  XNORM, IERR )
!
!                    Scale X to avoid overflow when updating
!                    the right-hand side.
!
                     IF( XNORM.GT.ONE ) THEN
                        BETA = MAX( WORK( J-1 ), WORK( J ) )
                        IF( BETA.GT.BIGNUM / XNORM ) THEN
                           REC = ONE / XNORM
                           X( 1, 1 ) = X( 1, 1 )*REC
                           X( 1, 2 ) = X( 1, 2 )*REC
                           X( 2, 1 ) = X( 2, 1 )*REC
                           X( 2, 2 ) = X( 2, 2 )*REC
                           SCALE = SCALE*REC
                        END IF
                     END IF
!
!                    Scale if necessary
!
                     IF( SCALE.NE.ONE ) THEN
                        CALL SSCAL( KI, SCALE, WORK( 1+N ), 1 )
                        CALL SSCAL( KI, SCALE, WORK( 1+N2 ), 1 )
                     END IF
                     WORK( J-1+N ) = X( 1, 1 )
                     WORK( J+N ) = X( 2, 1 )
                     WORK( J-1+N2 ) = X( 1, 2 )
                     WORK( J+N2 ) = X( 2, 2 )
!
!                    Update the right-hand side
!
                     CALL SAXPY( J-2, -X( 1, 1 ), T( 1, J-1 ), 1, &
                                 WORK( 1+N ), 1 )
                     CALL SAXPY( J-2, -X( 2, 1 ), T( 1, J ), 1, &
                                 WORK( 1+N ), 1 )
                     CALL SAXPY( J-2, -X( 1, 2 ), T( 1, J-1 ), 1, &
                                 WORK( 1+N2 ), 1 )
                     CALL SAXPY( J-2, -X( 2, 2 ), T( 1, J ), 1, &
                                 WORK( 1+N2 ), 1 )
                  END IF
   90          CONTINUE
!
!              Copy the vector x or Q*x to VR and normalize.
!
               IF( .NOT.OVER ) THEN
                  CALL SCOPY( KI, WORK( 1+N ), 1, VR( 1, IS-1 ), 1 )
                  CALL SCOPY( KI, WORK( 1+N2 ), 1, VR( 1, IS ), 1 )
!
                  EMAX = ZERO
                  DO 100 K = 1, KI
                     EMAX = MAX( EMAX, ABS( VR( K, IS-1 ) )+ &
                            ABS( VR( K, IS ) ) )
  100             CONTINUE
!
                  REMAX = ONE / EMAX
                  CALL SSCAL( KI, REMAX, VR( 1, IS-1 ), 1 )
                  CALL SSCAL( KI, REMAX, VR( 1, IS ), 1 )
!
                  DO 110 K = KI + 1, N
                     VR( K, IS-1 ) = ZERO
                     VR( K, IS ) = ZERO
  110             CONTINUE
!
               ELSE
!
                  IF( KI.GT.2 ) THEN
                     CALL SGEMV( 'N', N, KI-2, ONE, VR, LDVR, &
                                 WORK( 1+N ), 1, WORK( KI-1+N ), &
                                 VR( 1, KI-1 ), 1 )
                     CALL SGEMV( 'N', N, KI-2, ONE, VR, LDVR, &
                                 WORK( 1+N2 ), 1, WORK( KI+N2 ), &
                                 VR( 1, KI ), 1 )
                  ELSE
                     CALL SSCAL( N, WORK( KI-1+N ), VR( 1, KI-1 ), 1 )
                     CALL SSCAL( N, WORK( KI+N2 ), VR( 1, KI ), 1 )
                  END IF
!
                  EMAX = ZERO
                  DO 120 K = 1, N
                     EMAX = MAX( EMAX, ABS( VR( K, KI-1 ) )+ &
                            ABS( VR( K, KI ) ) )
  120             CONTINUE
                  REMAX = ONE / EMAX
                  CALL SSCAL( N, REMAX, VR( 1, KI-1 ), 1 )
                  CALL SSCAL( N, REMAX, VR( 1, KI ), 1 )
               END IF
            END IF
!
            IS = IS - 1
            IF( IP.NE.0 ) &
               IS = IS - 1
  130       CONTINUE
            IF( IP.EQ.1 ) &
               IP = 0
            IF( IP.EQ.-1 ) &
               IP = 1
  140    CONTINUE
      END IF
!
      IF( LEFTV ) THEN
!
!        Compute left eigenvectors.
!
         IP = 0
         IS = 1
         DO 260 KI = 1, N
!
            IF( IP.EQ.-1 ) &
               GO TO 250
            IF( KI.EQ.N ) &
               GO TO 150
            IF( T( KI+1, KI ).EQ.ZERO ) &
               GO TO 150
            IP = 1
!
  150       CONTINUE
            IF( SOMEV ) THEN
               IF( .NOT.SELECT( KI ) ) &
                  GO TO 250
            END IF
!
!           Compute the KI-th eigenvalue (WR,WI).
!
            WR = T( KI, KI )
            WI = ZERO
            IF( IP.NE.0 ) &
               WI = SQRT( ABS( T( KI, KI+1 ) ) )* &
                    SQRT( ABS( T( KI+1, KI ) ) )
            SMIN = MAX( ULP*( ABS( WR )+ABS( WI ) ), SMLNUM )
!
            IF( IP.EQ.0 ) THEN
!
!              Real left eigenvector.
!
               WORK( KI+N ) = ONE
!
!              Form right-hand side
!
               DO 160 K = KI + 1, N
                  WORK( K+N ) = -T( KI, K )
  160          CONTINUE
!
!              Solve the quasi-triangular system:
!                 (T(KI+1:N,KI+1:N) - WR)'*X = SCALE*WORK
!
               VMAX = ONE
               VCRIT = BIGNUM
!
               JNXT = KI + 1
               DO 170 J = KI + 1, N
                  IF( J.LT.JNXT ) &
                     GO TO 170
                  J1 = J
                  J2 = J
                  JNXT = J + 1
                  IF( J.LT.N ) THEN
                     IF( T( J+1, J ).NE.ZERO ) THEN
                        J2 = J + 1
                        JNXT = J + 2
                     END IF
                  END IF
!
                  IF( J1.EQ.J2 ) THEN
!
!                    1-by-1 diagonal block
!
!                    Scale if necessary to avoid overflow when forming
!                    the right-hand side.
!
                     IF( WORK( J ).GT.VCRIT ) THEN
                        REC = ONE / VMAX
                        CALL SSCAL( N-KI+1, REC, WORK( KI+N ), 1 )
                        VMAX = ONE
                        VCRIT = BIGNUM
                     END IF
!
                     WORK( J+N ) = WORK( J+N ) - &
                                   SDOT( J-KI-1, T( KI+1, J ), 1, &
                                   WORK( KI+1+N ), 1 )
!
!                    Solve (T(J,J)-WR)'*X = WORK
!
                     CALL SLALN2( .FALSE., 1, 1, SMIN, ONE, T( J, J ), &
                                  LDT, ONE, ONE, WORK( J+N ), N, WR, &
                                  ZERO, X, 2, SCALE, XNORM, IERR )
!
!                    Scale if necessary
!
                     IF( SCALE.NE.ONE ) &
                        CALL SSCAL( N-KI+1, SCALE, WORK( KI+N ), 1 )
                     WORK( J+N ) = X( 1, 1 )
                     VMAX = MAX( ABS( WORK( J+N ) ), VMAX )
                     VCRIT = BIGNUM / VMAX
!
                  ELSE
!
!                    2-by-2 diagonal block
!
!                    Scale if necessary to avoid overflow when forming
!                    the right-hand side.
!
                     BETA = MAX( WORK( J ), WORK( J+1 ) )
                     IF( BETA.GT.VCRIT ) THEN
                        REC = ONE / VMAX
                        CALL SSCAL( N-KI+1, REC, WORK( KI+N ), 1 )
                        VMAX = ONE
                        VCRIT = BIGNUM
                     END IF
!
                     WORK( J+N ) = WORK( J+N ) - &
                                   SDOT( J-KI-1, T( KI+1, J ), 1, &
                                   WORK( KI+1+N ), 1 )
!
                     WORK( J+1+N ) = WORK( J+1+N ) - &
                                     SDOT( J-KI-1, T( KI+1, J+1 ), 1, &
                                     WORK( KI+1+N ), 1 )
!
!                    Solve
!                      [T(J,J)-WR   T(J,J+1)     ]'* X = SCALE*( WORK1 )
!                      [T(J+1,J)    T(J+1,J+1)-WR]             ( WORK2 )
!
                     CALL SLALN2( .TRUE., 2, 1, SMIN, ONE, T( J, J ), &
                                  LDT, ONE, ONE, WORK( J+N ), N, WR, &
                                  ZERO, X, 2, SCALE, XNORM, IERR )
!
!                    Scale if necessary
!
                     IF( SCALE.NE.ONE ) &
                        CALL SSCAL( N-KI+1, SCALE, WORK( KI+N ), 1 )
                     WORK( J+N ) = X( 1, 1 )
                     WORK( J+1+N ) = X( 2, 1 )
!
                     VMAX = MAX( ABS( WORK( J+N ) ), &
                            ABS( WORK( J+1+N ) ), VMAX )
                     VCRIT = BIGNUM / VMAX
!
                  END IF
  170          CONTINUE
!
!              Copy the vector x or Q*x to VL and normalize.
!
               IF( .NOT.OVER ) THEN
                  CALL SCOPY( N-KI+1, WORK( KI+N ), 1, VL( KI, IS ), 1 )
!
                  II = ISAMAX( N-KI+1, VL( KI, IS ), 1 ) + KI - 1
                  REMAX = ONE / ABS( VL( II, IS ) )
                  CALL SSCAL( N-KI+1, REMAX, VL( KI, IS ), 1 )
!
                  DO 180 K = 1, KI - 1
                     VL( K, IS ) = ZERO
  180             CONTINUE
!
               ELSE
!
                  IF( KI.LT.N ) &
                     CALL SGEMV( 'N', N, N-KI, ONE, VL( 1, KI+1 ), LDVL, &
                                 WORK( KI+1+N ), 1, WORK( KI+N ), &
                                 VL( 1, KI ), 1 )
!
                  II = ISAMAX( N, VL( 1, KI ), 1 )
                  REMAX = ONE / ABS( VL( II, KI ) )
                  CALL SSCAL( N, REMAX, VL( 1, KI ), 1 )
!
               END IF
!
            ELSE
!
!              Complex left eigenvector.
!
!               Initial solve:
!                 ((T(KI,KI)    T(KI,KI+1) )' - (WR - I* WI))*X = 0.
!                 ((T(KI+1,KI) T(KI+1,KI+1))                )
!
               IF( ABS( T( KI, KI+1 ) ).GE.ABS( T( KI+1, KI ) ) ) THEN
                  WORK( KI+N ) = WI / T( KI, KI+1 )
                  WORK( KI+1+N2 ) = ONE
               ELSE
                  WORK( KI+N ) = ONE
                  WORK( KI+1+N2 ) = -WI / T( KI+1, KI )
               END IF
               WORK( KI+1+N ) = ZERO
               WORK( KI+N2 ) = ZERO
!
!              Form right-hand side
!
               DO 190 K = KI + 2, N
                  WORK( K+N ) = -WORK( KI+N )*T( KI, K )
                  WORK( K+N2 ) = -WORK( KI+1+N2 )*T( KI+1, K )
  190          CONTINUE
!
!              Solve complex quasi-triangular system:
!              ( T(KI+2,N:KI+2,N) - (WR-i*WI) )*X = WORK1+i*WORK2
!
               VMAX = ONE
               VCRIT = BIGNUM
!
               JNXT = KI + 2
               DO 200 J = KI + 2, N
                  IF( J.LT.JNXT ) &
                     GO TO 200
                  J1 = J
                  J2 = J
                  JNXT = J + 1
                  IF( J.LT.N ) THEN
                     IF( T( J+1, J ).NE.ZERO ) THEN
                        J2 = J + 1
                        JNXT = J + 2
                     END IF
                  END IF
!
                  IF( J1.EQ.J2 ) THEN
!
!                    1-by-1 diagonal block
!
!                    Scale if necessary to avoid overflow when
!                    forming the right-hand side elements.
!
                     IF( WORK( J ).GT.VCRIT ) THEN
                        REC = ONE / VMAX
                        CALL SSCAL( N-KI+1, REC, WORK( KI+N ), 1 )
                        CALL SSCAL( N-KI+1, REC, WORK( KI+N2 ), 1 )
                        VMAX = ONE
                        VCRIT = BIGNUM
                     END IF
!
                     WORK( J+N ) = WORK( J+N ) - &
                                   SDOT( J-KI-2, T( KI+2, J ), 1, &
                                   WORK( KI+2+N ), 1 )
                     WORK( J+N2 ) = WORK( J+N2 ) - &
                                    SDOT( J-KI-2, T( KI+2, J ), 1, &
                                    WORK( KI+2+N2 ), 1 )
!
!                    Solve (T(J,J)-(WR-i*WI))*(X11+i*X12)= WK+I*WK2
!
                     CALL SLALN2( .FALSE., 1, 2, SMIN, ONE, T( J, J ), &
                                  LDT, ONE, ONE, WORK( J+N ), N, WR, &
                                  -WI, X, 2, SCALE, XNORM, IERR )
!
!                    Scale if necessary
!
                     IF( SCALE.NE.ONE ) THEN
                        CALL SSCAL( N-KI+1, SCALE, WORK( KI+N ), 1 )
                        CALL SSCAL( N-KI+1, SCALE, WORK( KI+N2 ), 1 )
                     END IF
                     WORK( J+N ) = X( 1, 1 )
                     WORK( J+N2 ) = X( 1, 2 )
                     VMAX = MAX( ABS( WORK( J+N ) ), &
                            ABS( WORK( J+N2 ) ), VMAX )
                     VCRIT = BIGNUM / VMAX
!
                  ELSE
!
!                    2-by-2 diagonal block
!
!                    Scale if necessary to avoid overflow when forming
!                    the right-hand side elements.
!
                     BETA = MAX( WORK( J ), WORK( J+1 ) )
                     IF( BETA.GT.VCRIT ) THEN
                        REC = ONE / VMAX
                        CALL SSCAL( N-KI+1, REC, WORK( KI+N ), 1 )
                        CALL SSCAL( N-KI+1, REC, WORK( KI+N2 ), 1 )
                        VMAX = ONE
                        VCRIT = BIGNUM
                     END IF
!
                     WORK( J+N ) = WORK( J+N ) - &
                                   SDOT( J-KI-2, T( KI+2, J ), 1, &
                                   WORK( KI+2+N ), 1 )
!
                     WORK( J+N2 ) = WORK( J+N2 ) - &
                                    SDOT( J-KI-2, T( KI+2, J ), 1, &
                                    WORK( KI+2+N2 ), 1 )
!
                     WORK( J+1+N ) = WORK( J+1+N ) - &
                                     SDOT( J-KI-2, T( KI+2, J+1 ), 1, &
                                     WORK( KI+2+N ), 1 )
!
                     WORK( J+1+N2 ) = WORK( J+1+N2 ) - &
                                      SDOT( J-KI-2, T( KI+2, J+1 ), 1, &
                                      WORK( KI+2+N2 ), 1 )
!
!                    Solve 2-by-2 complex linear equation
!                      ([T(j,j)   T(j,j+1)  ]'-(wr-i*wi)*I)*X = SCALE*B
!                      ([T(j+1,j) T(j+1,j+1)]             )
!
                     CALL SLALN2( .TRUE., 2, 2, SMIN, ONE, T( J, J ), &
                                  LDT, ONE, ONE, WORK( J+N ), N, WR, &
                                  -WI, X, 2, SCALE, XNORM, IERR )
!
!                    Scale if necessary
!
                     IF( SCALE.NE.ONE ) THEN
                        CALL SSCAL( N-KI+1, SCALE, WORK( KI+N ), 1 )
                        CALL SSCAL( N-KI+1, SCALE, WORK( KI+N2 ), 1 )
                     END IF
                     WORK( J+N ) = X( 1, 1 )
                     WORK( J+N2 ) = X( 1, 2 )
                     WORK( J+1+N ) = X( 2, 1 )
                     WORK( J+1+N2 ) = X( 2, 2 )
                     VMAX = MAX( ABS( X( 1, 1 ) ), ABS( X( 1, 2 ) ), &
                            ABS( X( 2, 1 ) ), ABS( X( 2, 2 ) ), VMAX )
                     VCRIT = BIGNUM / VMAX
!
                  END IF
  200          CONTINUE
!
!              Copy the vector x or Q*x to VL and normalize.
!
  210          CONTINUE
               IF( .NOT.OVER ) THEN
                  CALL SCOPY( N-KI+1, WORK( KI+N ), 1, VL( KI, IS ), 1 )
                  CALL SCOPY( N-KI+1, WORK( KI+N2 ), 1, VL( KI, IS+1 ), &
                              1 )
!
                  EMAX = ZERO
                  DO 220 K = KI, N
                     EMAX = MAX( EMAX, ABS( VL( K, IS ) )+ &
                            ABS( VL( K, IS+1 ) ) )
  220             CONTINUE
                  REMAX = ONE / EMAX
                  CALL SSCAL( N-KI+1, REMAX, VL( KI, IS ), 1 )
                  CALL SSCAL( N-KI+1, REMAX, VL( KI, IS+1 ), 1 )
!
                  DO 230 K = 1, KI - 1
                     VL( K, IS ) = ZERO
                     VL( K, IS+1 ) = ZERO
  230             CONTINUE
               ELSE
                  IF( KI.LT.N-1 ) THEN
                     CALL SGEMV( 'N', N, N-KI-1, ONE, VL( 1, KI+2 ), &
                                 LDVL, WORK( KI+2+N ), 1, WORK( KI+N ), &
                                 VL( 1, KI ), 1 )
                     CALL SGEMV( 'N', N, N-KI-1, ONE, VL( 1, KI+2 ), &
                                 LDVL, WORK( KI+2+N2 ), 1, &
                                 WORK( KI+1+N2 ), VL( 1, KI+1 ), 1 )
                  ELSE
                     CALL SSCAL( N, WORK( KI+N ), VL( 1, KI ), 1 )
                     CALL SSCAL( N, WORK( KI+1+N2 ), VL( 1, KI+1 ), 1 )
                  END IF
!
                  EMAX = ZERO
                  DO 240 K = 1, N
                     EMAX = MAX( EMAX, ABS( VL( K, KI ) )+ &
                            ABS( VL( K, KI+1 ) ) )
  240             CONTINUE
                  REMAX = ONE / EMAX
                  CALL SSCAL( N, REMAX, VL( 1, KI ), 1 )
                  CALL SSCAL( N, REMAX, VL( 1, KI+1 ), 1 )
!
               END IF
!
            END IF
!
            IS = IS + 1
            IF( IP.NE.0 ) &
               IS = IS + 1
  250       CONTINUE
            IF( IP.EQ.-1 ) &
               IP = 0
            IF( IP.EQ.1 ) &
               IP = -1
!
  260    CONTINUE
!
      END IF
!
      RETURN
!
!     End of STREVC
!
      END SUBROUTINE STREVC

!--------1---------2---------3---------4---------5---------6---------7---------8

      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, &
                               N4 )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
!     ..
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPEC for a description of
!  the parameters.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR and QZ methods
!               for nonsymmetric eigenvalue problems.
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
! (ILAENV) (output) INTEGER
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 ) NB = MAX( 1, N )
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
!     ..
!     .. Executable Statements ..
!
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
!
!     Invalid value for ISPEC
!
      ILAENV = -1
      RETURN
!
  100 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
!
!        ASCII character set
!
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 ) &
                  SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
!
!        EBCDIC character set
!
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
             ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
             ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
                   ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
                   ( IC.GE.162 .AND. IC.LE.169 ) ) &
                  SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
!
!        Prime machines:  ASCII+128
!
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 ) &
                  SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
!
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) ) &
         RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
!
      GO TO ( 110, 200, 300 ) ISPEC
!
  110 CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
      NB = 1
!
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
                  C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
!
  200 CONTINUE
!
!     ISPEC = 2:  minimum block size
!
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
             C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
!
  300 CONTINUE
!
!     ISPEC = 3:  crossover point
!
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
             C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
                C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
                C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
!
  400 CONTINUE
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
      ILAENV = 6
      RETURN
!
  500 CONTINUE
!
!     ISPEC = 5:  minimum column dimension (not used)
!
      ILAENV = 2
      RETURN
!
  600 CONTINUE 
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      ILAENV = INT( REAL( MIN( N1, N2 ) ,r8)*1.6E0_r8 )
      RETURN
!
  700 CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
      ILAENV = 1
      RETURN
!
  800 CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      ILAENV = 50
      RETURN
!
!     End of ILAENV
!
      END FUNCTION ILAENV

!--------1---------2---------3---------4---------5---------6---------7---------8

      LOGICAL FUNCTION LSAME( CA, CB )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          CA, CB
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
! =====================================================================
!
!     ..
!     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      LSAME = CA.EQ.CB
      IF( LSAME ) &
         RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
!
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR. &
             INTA.GE.145 .AND. INTA.LE.153 .OR. &
             INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR. &
             INTB.GE.145 .AND. INTB.LE.153 .OR. &
             INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
!
!     RETURN
!
!     End of LSAME
!
      END FUNCTION LSAME

!--------1---------2---------3---------4---------5---------6---------7---------8

      FUNCTION SLAMCH( CMACH )

      real(r8) slamch

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992 
!
!     .. Scalar Arguments ..
      CHARACTER          CMACH
!     ..
!
!  Purpose
!  =======
!
!  SLAMCH determines single precision machine parameters.
!
!  Arguments
!  =========
!
!  CMACH   (input) CHARACTER*1
!          Specifies the value to be returned by SLAMCH:
!          = 'E' or 'e',   SLAMCH := eps
!          = 'S' or 's ,   SLAMCH := sfmin
!          = 'B' or 'b',   SLAMCH := base
!          = 'P' or 'p',   SLAMCH := eps*base
!          = 'N' or 'n',   SLAMCH := t
!          = 'R' or 'r',   SLAMCH := rnd
!          = 'M' or 'm',   SLAMCH := emin
!          = 'U' or 'u',   SLAMCH := rmin
!          = 'L' or 'l',   SLAMCH := emax
!          = 'O' or 'o',   SLAMCH := rmax
!
!          where
!
!          eps   = relative machine precision
!          sfmin = safe minimum, such that 1/sfmin does not overflow
!          base  = base of the machine
!          prec  = eps*base
!          t     = number of (base) digits in the mantissa
!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!          emin  = minimum exponent before (gradual) underflow
!          rmin  = underflow threshold - base**(emin-1)
!          emax  = largest exponent before overflow
!          rmax  = overflow threshold  - (base**emax)*(1-eps)
!
! =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0_r8, ZERO = 0.0E+0_r8 )
!     ..
!     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      REAL(r8)               BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN, &
                         RND, SFMIN, SMALL, T
!     ..
!     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN, &
                         EMAX, RMAX, PREC
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         CALL SLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
!
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
!
      SLAMCH = RMACH
      RETURN
!
!     End of SLAMCH
!
      END FUNCTION SLAMCH
!
!***********************************************************************
!
      SUBROUTINE SLAMC1( BETA, T, RND, IEEE1 )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
!     ..
!
!  Purpose
!  =======
!
!  SLAMC1 determines the machine parameters given by BETA, T, RND, and
!  IEEE1.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  IEEE1   (output) LOGICAL
!          Specifies whether rounding appears to be done in the IEEE
!          'round to nearest' style.
!
!  Further Details
!  ===============
!
!  The routine is based on the routine  ENVRON  by Malcolm and
!  incorporates suggestions by Gentleman and Marovich. See
!
!     Malcolm M. A. (1972) Algorithms to reveal properties of
!        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
!
!     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
!        that reveal properties of floating point arithmetic units.
!        Comms. of the ACM, 17, 276-277.
!
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      REAL(r8)               A, B, C, F, ONE, QTR, SAVEC, T1, T2
!     ..
!     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1
!
!        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
!        IEEE1, T and RND.
!
!        Throughout this routine  we use the function  SLAMC3  to ensure
!        that relevant values are  stored and not held in registers,  or
!        are not affected by optimizers.
!
!        Compute  a = 2.0**m  with the  smallest positive integer m such
!        that
!
!           fl( a + 1.0 ) = a.
!
         A = 1
         C = 1
!
!+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = SLAMC3( A, ONE )
            C = SLAMC3( C, -A )
            GO TO 10
         END IF
!+       END WHILE
!
!        Now compute  b = 2.0**m  with the smallest positive integer m
!        such that
!
!           fl( a + b ) .gt. a.
!
         B = 1
         C = SLAMC3( A, B )
!
!+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = SLAMC3( A, B )
            GO TO 20
         END IF
!+       END WHILE
!
!        Now compute the base.  a and c  are neighbouring floating point
!        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
!        their difference is beta. Adding 0.25 to c is to ensure that it
!        is truncated to beta and not ( beta - 1 ).
!
         QTR = ONE / 4
         SAVEC = C
         C = SLAMC3( C, -A )
         LBETA = C + QTR
!
!        Now determine whether rounding or chopping occurs,  by adding a
!        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
!
         B = LBETA
         F = SLAMC3( B / 2, -B / 100 )
         C = SLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = SLAMC3( B / 2, B / 100 )
         C = SLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) ) &
            LRND = .FALSE.
!
!        Try and decide whether rounding is done in the  IEEE  'round to
!        nearest' style. B/2 is half a unit in the last place of the two
!        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
!        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
!        A, but adding B/2 to SAVEC should change SAVEC.
!
         T1 = SLAMC3( B / 2, A )
         T2 = SLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
!
!        Now find  the  mantissa, t.  It should  be the  integer part of
!        log to the base beta of a,  however it is safer to determine  t
!        by powering.  So we find t as the smallest positive integer for
!        which
!
!           fl( beta**t + 1.0 ) = 1.0.
!
         LT = 0
         A = 1
         C = 1
!
!+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = SLAMC3( A, ONE )
            C = SLAMC3( C, -A )
            GO TO 30
         END IF
!+       END WHILE
!
      END IF
!
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN
!
!     End of SLAMC1
!
      END SUBROUTINE SLAMC1
!
!***********************************************************************
!
      SUBROUTINE SLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      REAL(r8)               EPS, RMAX, RMIN
!     ..
!
!  Purpose
!  =======
!
!  SLAMC2 determines the machine parameters specified in its argument
!  list.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  EPS     (output) REAL
!          The smallest positive number such that
!
!             fl( 1.0 - EPS ) .LT. 1.0,
!
!          where fl denotes the computed value.
!
!  EMIN    (output) INTEGER
!          The minimum exponent before (gradual) underflow occurs.
!
!  RMIN    (output) REAL
!          The smallest normalized number for the machine, given by
!          BASE**( EMIN - 1 ), where  BASE  is the floating point value
!          of BETA.
!
!  EMAX    (output) INTEGER
!          The maximum exponent before overflow occurs.
!
!  RMAX    (output) REAL
!          The largest positive number for the machine, given by
!          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
!          value of BETA.
!
!  Further Details
!  ===============
!
!  The computation of  EPS  is based on a routine PARANOIA by
!  W. Kahan of the University of California at Berkeley.
!
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT, &
                         NGNMIN, NGPMIN
      REAL(r8)               A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE, &
                         SIXTH, SMALL, THIRD, TWO, ZERO
!     ..
!     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX, &
                         LRMIN, LT
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2
!
!        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
!        BETA, T, RND, EPS, EMIN and RMIN.
!
!        Throughout this routine  we use the function  SLAMC3  to ensure
!        that relevant values are stored  and not held in registers,  or
!        are not affected by optimizers.
!
!        SLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
!
         CALL SLAMC1( LBETA, LT, LRND, LIEEE1 )
!
!        Start to find EPS.
!
         B = LBETA
         A = B**( -LT )
         LEPS = A
!
!        Try some tricks to see whether or not this is the correct  EPS.
!
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = SLAMC3( B, -HALF )
         THIRD = SLAMC3( SIXTH, SIXTH )
         B = SLAMC3( THIRD, -HALF )
         B = SLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS ) &
            B = LEPS
!
         LEPS = 1
!
!+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = SLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = SLAMC3( HALF, -C )
            B = SLAMC3( HALF, C )
            C = SLAMC3( HALF, -B )
            B = SLAMC3( HALF, C )
            GO TO 10
         END IF
!+       END WHILE
!
         IF( A.LT.LEPS ) &
            LEPS = A
!
!        Computation of EPS complete.
!
!        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
!        Keep dividing  A by BETA until (gradual) underflow occurs. This
!        is detected when we cannot recover the previous A.
!
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = SLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = SLAMC3( ONE, SMALL )
         CALL SLAMC4( NGPMIN, ONE, LBETA )
         CALL SLAMC4( NGNMIN, -ONE, LBETA )
         CALL SLAMC4( GPMIN, A, LBETA )
         CALL SLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
!
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
!            ( Non twos-complement machines, no gradual underflow;
!              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
!            ( Non twos-complement machines, with gradual underflow;
!              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
!            ( Twos-complement machines, no gradual underflow;
!              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND. &
                  ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
!            ( Twos-complement machines with gradual underflow;
!              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
!         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
!**
! Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            write(iulog, FMT = 9999 )LEMIN
         END IF
!**
!
!        Assume IEEE arithmetic if we found denormalised  numbers above,
!        or if arithmetic seems to round in the  IEEE style,  determined
!        in routine SLAMC1. A true IEEE machine should have both  things
!        true; however, faulty machines may have one or the other.
!
         IEEE = IEEE .OR. LIEEE1
!
!        Compute  RMIN by successive division by  BETA. We could compute
!        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
!        this computation.
!
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = SLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
!
!        Finally, call SLAMC5 to compute EMAX and RMAX.
!
         CALL SLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
!
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
!
      RETURN
!
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-', &
            '  EMIN = ', I8, / &
            ' If, after inspection, the value EMIN looks', &
            ' acceptable please comment out ', &
            / ' the IF block as marked within the code of routine', &
            ' SLAMC2,', / ' otherwise supply EMIN explicitly.', / )
!
!     End of SLAMC2
!
      END SUBROUTINE SLAMC2
!
!***********************************************************************
!
      FUNCTION SLAMC3( A, B )

      real(r8) slamc3

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      REAL(r8)               A, B
!     ..
!
!  Purpose
!  =======
!
!  SLAMC3  is intended to force  A  and  B  to be stored prior to doing
!  the addition of  A  and  B ,  for use in situations where optimizers
!  might hold one of these in a register.
!
!  Arguments
!  =========
!
!  A, B    (input) REAL
!          The values A and B.
!
! =====================================================================
!
!     .. Executable Statements ..
!
      SLAMC3 = A + B
!
      RETURN
!
!     End of SLAMC3
!
      END FUNCTION SLAMC3
!
!***********************************************************************
!
      SUBROUTINE SLAMC4( EMIN, START, BASE )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      REAL(r8)               START
!     ..
!
!  Purpose
!  =======
!
!  SLAMC4 is a service routine for SLAMC2.
!
!  Arguments
!  =========
!
!  EMIN    (output) EMIN
!          The minimum exponent before (gradual) underflow, computed by
!          setting A = START and dividing by BASE until the previous A
!          can not be recovered.
!
!  START   (input) REAL
!          The starting point for determining EMIN.
!
!  BASE    (input) INTEGER
!          The base of the machine.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I
      REAL(r8)               A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
!     ..
!     .. Executable Statements ..
!
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = SLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
!+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
!    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND. &
          ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = SLAMC3( A / BASE, ZERO )
         C1 = SLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = SLAMC3( A*RBASE, ZERO )
         C2 = SLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
!+    END WHILE
!
      RETURN
!
!     End of SLAMC4
!
      END SUBROUTINE SLAMC4
!
!***********************************************************************
!
      SUBROUTINE SLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      REAL(r8)               RMAX
!     ..
!
!  Purpose
!  =======
!
!  SLAMC5 attempts to compute RMAX, the largest machine floating-point
!  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
!  approximately to a power of 2.  It will fail on machines where this
!  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
!  EMAX = 28718).  It will also fail if the value supplied for EMIN is
!  too large (i.e. too close to zero), probably with overflow.
!
!  Arguments
!  =========
!
!  BETA    (input) INTEGER
!          The base of floating-point arithmetic.
!
!  P       (input) INTEGER
!          The number of base BETA digits in the mantissa of a
!          floating-point value.
!
!  EMIN    (input) INTEGER
!          The minimum exponent before (gradual) underflow.
!
!  IEEE    (input) LOGICAL
!          A logical flag specifying whether or not the arithmetic
!          system is thought to comply with the IEEE standard.
!
!  EMAX    (output) INTEGER
!          The largest exponent before overflow
!
!  RMAX    (output) REAL
!          The largest machine floating-point number.
!
! =====================================================================
!
!     .. Parameters ..
      REAL(r8)               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0_r8, ONE = 1.0E0_r8 )
!     ..
!     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      REAL(r8)               OLDY, RECBAS, Y, Z
!     ..
!     .. Executable Statements ..
!
!     First compute LEXP and UEXP, two powers of 2 that bound
!     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
!     approximately to the bound that is closest to abs(EMIN).
!     (EMAX is the exponent of the required number RMAX).
!
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
!
!     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
!     than or equal to EMIN. EXBITS is the number of bits needed to
!     store the exponent.
!
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
!
!     EXPSUM is the exponent range, approximately equal to
!     EMAX - EMIN + 1 .
!
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
!
!     NBITS is the total number of bits needed to store a
!     floating-point number.
!
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
!
!        Either there are an odd number of bits used to store a
!        floating-point number, which is unlikely, or some bits are
!        not used in the representation of numbers, which is possible,
!        (e.g. Cray machines) or the mantissa has an implicit bit,
!        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
!        most likely. We have to assume the last alternative.
!        If this is true, then we need to reduce EMAX by one because
!        there must be some way of representing zero in an implicit-bit
!        system. On machines like Cray, we are reducing EMAX by one
!        unnecessarily.
!
         EMAX = EMAX - 1
      END IF
!
      IF( IEEE ) THEN
!
!        Assume we are on an IEEE machine which reserves one exponent
!        for infinity and NaN.
!
         EMAX = EMAX - 1
      END IF
!
!     Now create RMAX, the largest machine number, which should
!     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
!
!     First compute 1.0 - BETA**(-P), being careful that the
!     result is less than 1.0 .
!
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE ) &
            OLDY = Y
         Y = SLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE ) &
         Y = OLDY
!
!     Now multiply by BETA**EMAX to get RMAX.
!
      DO 30 I = 1, EMAX
         Y = SLAMC3( Y*BETA, ZERO )
   30 CONTINUE
!
      RMAX = Y
      RETURN
!
!     End of SLAMC5
!
      END SUBROUTINE SLAMC5

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE XERBLA( SRNAME, INFO )

!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
!     ..
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========
!
!  SRNAME  (input) CHARACTER*6
!          The name of the routine which called XERBLA.
!
!  INFO    (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
! =====================================================================
!
!     .. Executable Statements ..
!
      write(iulog, FMT = 9999 )SRNAME, INFO
!
      call endrun('XERBLA')
!
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ', &
            'an illegal value' )
!
!     End of XERBLA
!
      END SUBROUTINE XERBLA

!--------1---------2---------3---------4---------5---------6---------7---------8

      subroutine scopy(n,sx,incx,sy,incy)

!
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to 1.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      real(r8) sx(*),sy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        sy(i) = sx(i)
        sy(i + 1) = sx(i + 1)
        sy(i + 2) = sx(i + 2)
        sy(i + 3) = sx(i + 3)
        sy(i + 4) = sx(i + 4)
        sy(i + 5) = sx(i + 5)
        sy(i + 6) = sx(i + 6)
   50 continue
      return
      end subroutine scopy

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )

!     .. Scalar Arguments ..
      REAL(r8)               ALPHA
      INTEGER            INCX, INCY, LDA, M, N
!     .. Array Arguments ..
      REAL(r8)               A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  SGER   performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Parameters
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - REAL             array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - REAL             array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - REAL             array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      REAL(r8)               ZERO
      PARAMETER        ( ZERO = 0.0E+0_r8 )
!     .. Local Scalars ..
      REAL(r8)               TEMP
      INTEGER            I, INFO, IX, J, JY, KX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SGER  ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) ) &
         RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
!
      RETURN
!
!     End of SGER  .
!
      END SUBROUTINE SGER

!--------1---------2---------3---------4---------5---------6---------7---------8

      FUNCTION SNRM2 ( N, X, INCX )

      real(r8) snrm2

!     .. Scalar Arguments ..
      INTEGER                           INCX, N
!     .. Array Arguments ..
      REAL(r8)                              X( * )
!     ..
!
!  SNRM2 returns the euclidean norm of a vector via the function
!  name, so that
!
!     SNRM2 := sqrt( x'*x )
!
!
!
!  -- This version written on 25-October-1982.
!     Modified on 14-October-1993 to inline the call to SLASSQ.
!     Sven Hammarling, Nag Ltd.
!
!
!     .. Parameters ..
      REAL(r8)                  ONE         , ZERO
      PARAMETER           ( ONE = 1.0E+0_r8, ZERO = 0.0E+0_r8 )
!     .. Local Scalars ..
      INTEGER               IX
      REAL(r8)                  ABSXI, NORM, SCALE, SSQ
!     ..
!     .. Executable Statements ..
      IF( N.LT.1 .OR. INCX.LT.1 )THEN
         NORM  = ZERO
      ELSE IF( N.EQ.1 )THEN
         NORM  = ABS( X( 1 ) )
      ELSE
         SCALE = ZERO
         SSQ   = ONE
!        The following loop is equivalent to this call to the LAPACK
!        auxiliary routine:
!        CALL SLASSQ( N, X, INCX, SCALE, SSQ )
!
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( X( IX ).NE.ZERO )THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI )THEN
                  SSQ   = ONE   + SSQ*( SCALE/ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SSQ   = SSQ   +     ( ABSXI/SCALE )**2
               END IF
            END IF
   10    CONTINUE
         NORM  = SCALE * SQRT( SSQ )
      END IF
!
      SNRM2 = NORM
      RETURN
!
!     End of SNRM2.
!
      END FUNCTION SNRM2

!--------1---------2---------3---------4---------5---------6---------7---------8

      subroutine srot (n,sx,incx,sy,incy,c,s)
!
!     applies a plane rotation.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!

      real(r8) sx(*),sy(*),stemp,c,s
      integer i,incx,incy,ix,iy,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = c*sx(ix) + s*sy(iy)
        sy(iy) = c*sy(iy) - s*sx(ix)
        sx(ix) = stemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!       code for both increments equal to 1
!
   20 do 30 i = 1,n
        stemp = c*sx(i) + s*sy(i)
        sy(i) = c*sy(i) - s*sx(i)
        sx(i) = stemp
   30 continue
      return
      end subroutine srot

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
                         BETA, C, LDC )

!     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      REAL(r8)               ALPHA, BETA
!     .. Array Arguments ..
      REAL(r8)               A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  SGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - REAL            .
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - REAL             array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      REAL(r8)               TEMP
!     .. Parameters ..
      REAL(r8)               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0_r8, ZERO = 0.0E+0_r8 )
!     ..
!     .. Executable Statements ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
!
!     Test the input parameters.
!
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND. &
               ( .NOT.LSAME( TRANSA, 'C' ) ).AND. &
               ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND. &
               ( .NOT.LSAME( TRANSB, 'C' ) ).AND. &
               ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SGEMM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR. &
          ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) ) &
         RETURN
!
!     And if  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
!
!     Start the operations.
!
      IF( NOTB )THEN
         IF( NOTA )THEN
!
!           Form  C := alpha*A*B + beta*C.
!
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B + beta*C
!
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
!
!           Form  C := alpha*A*B' + beta*C
!
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B' + beta*C
!
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of SGEMM .
!
      END SUBROUTINE SGEMM

!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE SGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, &
                         BETA, Y, INCY )

!     .. Scalar Arguments ..
      REAL(r8)               ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
!     .. Array Arguments ..
      REAL(r8)               A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  SGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - REAL             array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - REAL             array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - REAL            .
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - REAL             array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      REAL(r8)               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0_r8, ZERO = 0.0E+0_r8 )
!     .. Local Scalars ..
      REAL(r8)               TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND. &
               .NOT.LSAME( TRANS, 'T' ).AND. &
               .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SGEMV ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR. &
          ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) ) &
         RETURN
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO ) &
         RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
!
!        Form  y := alpha*A*x + y.
!
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
!
!        Form  y := alpha*A'*x + y.
!
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of SGEMV .
!
      END SUBROUTINE SGEMV
 
!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE STRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
                         B, LDB )

!     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      REAL(r8)               ALPHA
!     .. Array Arguments ..
      REAL(r8)               A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  STRMM  performs one of the matrix-matrix operations
!
!     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
!
!  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry,  SIDE specifies whether  op( A ) multiplies B from
!           the left or right as follows:
!
!              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
!
!              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL            .
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - REAL             array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - REAL             array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain the matrix  B,  and  on exit  is overwritten  by the
!           transformed matrix.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      REAL(r8)               TEMP
!     .. Parameters ..
      REAL(r8)               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0_r8, ZERO = 0.0E+0_r8 )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
!
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND. &
               ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND. &
               ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND. &
               ( .NOT.LSAME( TRANSA, 'T' ) ).AND. &
               ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND. &
               ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'STRMM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 ) &
         RETURN
!
!     And when  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
!
!     Start the operations.
!
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*A*B.
!
            IF( UPPER )THEN
               DO 50, J = 1, N
                  DO 40, K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*B( K, J )
                        DO 30, I = 1, K - 1
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   30                   CONTINUE
                        IF( NOUNIT ) &
                           TEMP = TEMP*A( K, K )
                        B( K, J ) = TEMP
                     END IF
   40             CONTINUE
   50          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70 K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP      = ALPHA*B( K, J )
                        B( K, J ) = TEMP
                        IF( NOUNIT ) &
                           B( K, J ) = B( K, J )*A( K, K )
                        DO 60, I = K + 1, M
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   60                   CONTINUE
                     END IF
   70             CONTINUE
   80          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*A'*B.
!
            IF( UPPER )THEN
               DO 110, J = 1, N
                  DO 100, I = M, 1, -1
                     TEMP = B( I, J )
                     IF( NOUNIT ) &
                        TEMP = TEMP*A( I, I )
                     DO 90, K = 1, I - 1
                        TEMP = TEMP + A( K, I )*B( K, J )
   90                CONTINUE
                     B( I, J ) = ALPHA*TEMP
  100             CONTINUE
  110          CONTINUE
            ELSE
               DO 140, J = 1, N
                  DO 130, I = 1, M
                     TEMP = B( I, J )
                     IF( NOUNIT ) &
                        TEMP = TEMP*A( I, I )
                     DO 120, K = I + 1, M
                        TEMP = TEMP + A( K, I )*B( K, J )
  120                CONTINUE
                     B( I, J ) = ALPHA*TEMP
  130             CONTINUE
  140          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*B*A.
!
            IF( UPPER )THEN
               DO 180, J = N, 1, -1
                  TEMP = ALPHA
                  IF( NOUNIT ) &
                     TEMP = TEMP*A( J, J )
                  DO 150, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  150             CONTINUE
                  DO 170, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 160, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  160                   CONTINUE
                     END IF
  170             CONTINUE
  180          CONTINUE
            ELSE
               DO 220, J = 1, N
                  TEMP = ALPHA
                  IF( NOUNIT ) &
                     TEMP = TEMP*A( J, J )
                  DO 190, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  190             CONTINUE
                  DO 210, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 200, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  200                   CONTINUE
                     END IF
  210             CONTINUE
  220          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*B*A'.
!
            IF( UPPER )THEN
               DO 260, K = 1, N
                  DO 240, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT ) &
                     TEMP = TEMP*A( K, K )
                  IF( TEMP.NE.ONE )THEN
                     DO 250, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  250                CONTINUE
                  END IF
  260          CONTINUE
            ELSE
               DO 300, K = N, 1, -1
                  DO 280, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 270, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  270                   CONTINUE
                     END IF
  280             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT ) &
                     TEMP = TEMP*A( K, K )
                  IF( TEMP.NE.ONE )THEN
                     DO 290, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
                  END IF
  300          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of STRMM .
!
      END SUBROUTINE STRMM
 
!--------1---------2---------3---------4---------5---------6---------7---------8

      SUBROUTINE STRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )

!     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
!     .. Array Arguments ..
      REAL(r8)               A( LDA, * ), X( * )
!     ..
!
!  Purpose
!  =======
!
!  STRMV  performs one of the matrix-vector operations
!
!     x := A*x,   or   x := A'*x,
!
!  where x is an n element vector and  A is an n by n unit, or non-unit,
!  upper or lower triangular matrix.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   x := A*x.
!
!              TRANS = 'T' or 't'   x := A'*x.
!
!              TRANS = 'C' or 'c'   x := A'*x.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  A      - REAL             array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular matrix and the strictly lower triangular part of
!           A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular matrix and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!           A are not referenced either, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!  X      - REAL             array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x. On exit, X is overwritten with the
!           tranformed vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      REAL(r8)               ZERO
      PARAMETER        ( ZERO = 0.0E+0_r8 )
!     .. Local Scalars ..
      REAL(r8)               TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND. &
               .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND. &
               .NOT.LSAME( TRANS, 'T' ).AND. &
               .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND. &
               .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'STRMV ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 ) &
         RETURN
!
      NOUNIT = LSAME( DIAG, 'N' )
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF( LSAME( TRANS, 'N' ) )THEN
!
!        Form  x := A*x.
!
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 10, I = 1, J - 1
                        X( I ) = X( I ) + TEMP*A( I, J )
   10                CONTINUE
                     IF( NOUNIT ) &
                        X( J ) = X( J )*A( J, J )
                  END IF
   20          CONTINUE
            ELSE
               JX = KX
               DO 40, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 30, I = 1, J - 1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      + INCX
   30                CONTINUE
                     IF( NOUNIT ) &
                        X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX + INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 50, I = N, J + 1, -1
                        X( I ) = X( I ) + TEMP*A( I, J )
   50                CONTINUE
                     IF( NOUNIT ) &
                        X( J ) = X( J )*A( J, J )
                  END IF
   60          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 80, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 70, I = N, J + 1, -1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      - INCX
   70                CONTINUE
                     IF( NOUNIT ) &
                        X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX - INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
!
!        Form  x := A'*x.
!
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 100, J = N, 1, -1
                  TEMP = X( J )
                  IF( NOUNIT ) &
                     TEMP = TEMP*A( J, J )
                  DO 90, I = J - 1, 1, -1
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 120, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT ) &
                     TEMP = TEMP*A( J, J )
                  DO 110, I = J - 1, 1, -1
                     IX   = IX   - INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  110             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = 1, N
                  TEMP = X( J )
                  IF( NOUNIT ) &
                     TEMP = TEMP*A( J, J )
                  DO 130, I = J + 1, N
                     TEMP = TEMP + A( I, J )*X( I )
  130             CONTINUE
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               JX = KX
               DO 160, J = 1, N
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT ) &
                     TEMP = TEMP*A( J, J )
                  DO 150, I = J + 1, N
                     IX   = IX   + INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  150             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  160          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of STRMV .
!
      END SUBROUTINE STRMV
 
!--------1---------2---------3---------4---------5---------6---------7---------8

end module sgexx
