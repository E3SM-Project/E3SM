! FFT99F
!
! PURPOSE      PERFORMS MULTIPLE FAST FOURIER TRANSFORMS.  THIS PACKAGE
!              WILL PERFORM A NUMBER OF SIMULTANEOUS REAL/HALF-COMPLEX
!              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE
!              TRANSFORMS, I.E.  GIVEN A SET OF REAL DATA VECTORS, THE
!              PACKAGE RETURNS A SET OF 'HALF-COMPLEX' FOURIER
!              COEFFICIENT VECTORS, OR VICE VERSA.  THE LENGTH OF THE
!              TRANSFORMS MUST BE AN EVEN NUMBER GREATER THAN 4 THAT HAS
!              NO OTHER FACTORS EXCEPT POSSIBLY POWERS OF 2, 3, AND 5.
!              THIS IS AN ALL FORTRAN VERSION OF THE CRAYLIB PACKAGE
!              THAT IS MOSTLY WRITTEN IN CAL.
!
!              THE PACKAGE FFT99F CONTAINS SEVERAL USER-LEVEL ROUTINES:
!
!            SUBROUTINE SET99
!                AN INITIALIZATION ROUTINE THAT MUST BE CALLED ONCE
!                BEFORE A SEQUENCE OF CALLS TO THE FFT ROUTINES
!                (PROVIDED THAT N IS NOT CHANGED).
!
!            SUBROUTINES FFT99 AND FFT991
!                TWO FFT ROUTINES THAT RETURN SLIGHTLY DIFFERENT
!                ARRANGEMENTS OF THE DATA IN GRIDPOINT SPACE.
!
!
! ACCESS       THIS FORTRAN VERSION MAY BE ACCESSED WITH
!
!                   *FORTRAN,P=XLIB,SN=FFT99F
!
!              TO ACCESS THE CRAY OBJECT CODE, CALLING THE USER ENTRY
!              POINTS FROM A CRAY PROGRAM IS SUFFICIENT.  THE SOURCE
!              FORTRAN AND CAL CODE FOR THE CRAYLIB VERSION MAY BE
!              ACCESSED USING
!
!                   FETCH P=CRAYLIB,SN=FFT99
!                   FETCH P=CRAYLIB,SN=CAL99
!
! USAGE        LET N BE OF THE FORM 2**P * 3**Q * 5**R, WHERE P .GE. 1,
!              Q .GE. 0, AND R .GE. 0.  THEN A TYPICAL SEQUENCE OF
!              CALLS TO TRANSFORM A GIVEN SET OF REAL VECTORS OF LENGTH
!              N TO A SET OF 'HALF-COMPLEX' FOURIER COEFFICIENT VECTORS
!              OF LENGTH N IS
!
!                   DIMENSION IFAX(13),TRIGS(3*N/2+1),A(M*(N+2)),
!                  +          WORK(M*(N+1))
!
!                   CALL SET99 (TRIGS, IFAX, N)
!                   CALL FFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)
!
!              SEE THE INDIVIDUAL WRITE-UPS FOR SET99, FFT99, AND
!              FFT991 BELOW, FOR A DETAILED DESCRIPTION OF THE
!              ARGUMENTS.
!
! HISTORY      THE PACKAGE WAS WRITTEN BY CLIVE TEMPERTON AT ECMWF IN
!              NOVEMBER, 1978.  IT WAS MODIFIED, DOCUMENTED, AND TESTED
!              FOR NCAR BY RUSS REW IN SEPTEMBER, 1980.
!
!-----------------------------------------------------------------------
!
! SUBROUTINE SET99 (TRIGS, IFAX, N)
!
! PURPOSE      A SET-UP ROUTINE FOR FFT99 AND FFT991.  IT NEED ONLY BE
!              CALLED ONCE BEFORE A SEQUENCE OF CALLS TO THE FFT
!              ROUTINES (PROVIDED THAT N IS NOT CHANGED).
!
! ARGUMENT     IFAX(13),TRIGS(3*N/2+1)
! DIMENSIONS
!
! ARGUMENTS
!
! ON INPUT     TRIGS
!               A FLOATING POINT ARRAY OF DIMENSION 3*N/2 IF N/2 IS
!               EVEN, OR 3*N/2+1 IF N/2 IS ODD.
!
!              IFAX
!               AN INTEGER ARRAY.  THE NUMBER OF ELEMENTS ACTUALLY USED
!               WILL DEPEND ON THE FACTORIZATION OF N.  DIMENSIONING
!               IFAX FOR 13 SUFFICES FOR ALL N LESS THAN A MILLION.
!
!              N
!               AN EVEN NUMBER GREATER THAN 4 THAT HAS NO PRIME FACTOR
!               GREATER THAN 5.  N IS THE LENGTH OF THE TRANSFORMS (SEE
!               THE DOCUMENTATION FOR FFT99 AND FFT991 FOR THE
!               DEFINITIONS OF THE TRANSFORMS).
!
! ON OUTPUT    IFAX
!               CONTAINS THE FACTORIZATION OF N/2.  IFAX(1) IS THE
!               NUMBER OF FACTORS, AND THE FACTORS THEMSELVES ARE STORED
!               IN IFAX(2),IFAX(3),...  IF SET99 IS CALLED WITH N ODD,
!               OR IF N HAS ANY PRIME FACTORS GREATER THAN 5, IFAX(1)
!               IS SET TO -99.
!
!              TRIGS
!               AN ARRAY OF TRIGONOMETRIC FUNCTION VALUES SUBSEQUENTLY
!               USED BY THE FFT ROUTINES.
!
!-----------------------------------------------------------------------
!
! SUBROUTINE FFT991 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)
!                       AND
! SUBROUTINE FFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)
!
! PURPOSE      PERFORM A NUMBER OF SIMULTANEOUS REAL/HALF-COMPLEX
!              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE
!              TRANSFORMS, USING ORDINARY SPATIAL ORDER OF GRIDPOINT
!              VALUES (FFT991) OR EXPLICIT CYCLIC CONTINUITY IN THE
!              GRIDPOINT VALUES (FFT99).  GIVEN A SET
!              OF REAL DATA VECTORS, THE PACKAGE RETURNS A SET OF
!              'HALF-COMPLEX' FOURIER COEFFICIENT VECTORS, OR VICE
!              VERSA.  THE LENGTH OF THE TRANSFORMS MUST BE AN EVEN
!              NUMBER THAT HAS NO OTHER FACTORS EXCEPT POSSIBLY POWERS
!              OF 2, 3, AND 5.  THESE VERSION OF FFT991 AND FFT99 ARE
!              OPTIMIZED FOR USE ON THE CRAY-1.
!
! ARGUMENT     A(M*(N+2)), WORK(M*(N+1)), TRIGS(3*N/2+1), IFAX(13)
! DIMENSIONS
!
! ARGUMENTS
!
! ON INPUT     A
!               AN ARRAY OF LENGTH M*(N+2) CONTAINING THE INPUT DATA
!               OR COEFFICIENT VECTORS.  THIS ARRAY IS OVERWRITTEN BY
!               THE RESULTS.
!
!              WORK
!               A WORK ARRAY OF DIMENSION M*(N+1)
!
!              TRIGS
!               AN ARRAY SET UP BY SET99, WHICH MUST BE CALLED FIRST.
!
!              IFAX
!               AN ARRAY SET UP BY SET99, WHICH MUST BE CALLED FIRST.
!
!              INC
!               THE INCREMENT (IN WORDS) BETWEEN SUCCESSIVE ELEMENTS OF
!               EACH DATA OR COEFFICIENT VECTOR (E.G.  INC=1 FOR
!               CONSECUTIVELY STORED DATA).
!
!              JUMP
!               THE INCREMENT (IN WORDS) BETWEEN THE FIRST ELEMENTS OF
!               SUCCESSIVE DATA OR COEFFICIENT VECTORS.  ON THE CRAY-1,
!               TRY TO ARRANGE DATA SO THAT JUMP IS NOT A MULTIPLE OF 8
!               (TO AVOID MEMORY BANK CONFLICTS).  FOR CLARIFICATION OF
!               INC AND JUMP, SEE THE EXAMPLES BELOW.
!
!              N
!               THE LENGTH OF EACH TRANSFORM (SEE DEFINITION OF
!               TRANSFORMS, BELOW).
!
!              M
!               THE NUMBER OF TRANSFORMS TO BE DONE SIMULTANEOUSLY.
!
!              ISIGN
!               = +1 FOR A TRANSFORM FROM FOURIER COEFFICIENTS TO
!                    GRIDPOINT VALUES.
!               = -1 FOR A TRANSFORM FROM GRIDPOINT VALUES TO FOURIER
!                    COEFFICIENTS.
!
! ON OUTPUT    A
!               IF ISIGN = +1, AND M COEFFICIENT VECTORS ARE SUPPLIED
!               EACH CONTAINING THE SEQUENCE:
!
!               A(0),B(0),A(1),B(1),...,A(N/2),B(N/2)  (N+2 VALUES)
!
!               THEN THE RESULT CONSISTS OF M DATA VECTORS EACH
!               CONTAINING THE CORRESPONDING N+2 GRIDPOINT VALUES:
!
!               FOR FFT991, X(0), X(1), X(2),...,X(N-1),0,0.
!               FOR FFT99, X(N-1),X(0),X(1),X(2),...,X(N-1),X(0).
!                   (EXPLICIT CYCLIC CONTINUITY)
!
!               WHEN ISIGN = +1, THE TRANSFORM IS DEFINED BY:
!                 X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!                 WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!                 AND I=SQRT (-1)
!
!               IF ISIGN = -1, AND M DATA VECTORS ARE SUPPLIED EACH
!               CONTAINING A SEQUENCE OF GRIDPOINT VALUES X(J) AS
!               DEFINED ABOVE, THEN THE RESULT CONSISTS OF M VECTORS
!               EACH CONTAINING THE CORRESPONDING FOURIER COFFICIENTS
!               A(K), B(K), 0 .LE. K .LE N/2.
!
!               WHEN ISIGN = -1, THE INVERSE TRANSFORM IS DEFINED BY:
!                 C(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*EXP(-2*I*J*K*PI/N))
!                 WHERE C(K)=A(K)+I*B(K) AND I=SQRT(-1)
!
!               A CALL WITH ISIGN=+1 FOLLOWED BY A CALL WITH ISIGN=-1
!               (OR VICE VERSA) RETURNS THE ORIGINAL DATA.
!
!               NOTE: THE FACT THAT THE GRIDPOINT VALUES X(J) ARE REAL
!               IMPLIES THAT B(0)=B(N/2)=0.  FOR A CALL WITH ISIGN=+1,
!               IT IS NOT ACTUALLY NECESSARY TO SUPPLY THESE ZEROS.
!
! EXAMPLES      GIVEN 19 DATA VECTORS EACH OF LENGTH 64 (+2 FOR EXPLICIT
!               CYCLIC CONTINUITY), COMPUTE THE CORRESPONDING VECTORS OF
!               FOURIER COEFFICIENTS.  THE DATA MAY, FOR EXAMPLE, BE
!               ARRANGED LIKE THIS:
!
! FIRST DATA   A(1)=    . . .                A(66)=             A(70)
! VECTOR       X(63) X(0) X(1) X(2) ... X(63) X(0)  (4 EMPTY LOCATIONS)
!
! SECOND DATA  A(71)=   . . .                                  A(140)
! VECTOR       X(63) X(0) X(1) X(2) ... X(63) X(0)  (4 EMPTY LOCATIONS)
!
!               AND SO ON.  HERE INC=1, JUMP=70, N=64, M=19, ISIGN=-1,
!               AND FFT99 SHOULD BE USED (BECAUSE OF THE EXPLICIT CYCLIC
!               CONTINUITY).
!
!               ALTERNATIVELY THE DATA MAY BE ARRANGED LIKE THIS:
!
!                FIRST         SECOND                          LAST
!                DATA          DATA                            DATA
!                VECTOR        VECTOR                          VECTOR
!
!                 A(1)=         A(2)=                           A(19)=
!
!                 X(63)         X(63)       . . .               X(63)
!        A(20)=   X(0)          X(0)        . . .               X(0)
!        A(39)=   X(1)          X(1)        . . .               X(1)
!                  .             .                               .
!                  .             .                               .
!                  .             .                               .
!
!               IN WHICH CASE WE HAVE INC=19, JUMP=1, AND THE REMAINING
!               PARAMETERS ARE THE SAME AS BEFORE.  IN EITHER CASE, EACH
!               COEFFICIENT VECTOR OVERWRITES THE CORRESPONDING INPUT
!               DATA VECTOR.
!-----------------------------------------------------------------------
!
! $Id$
! $Author$

!================================================================================================

      SUBROUTINE FFT99(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
!
!-----------------------------------------------------------------------
!     SUBROUTINE "FFT99" - MULTIPLE FAST REAL PERIODIC TRANSFORM
!     CORRESPONDING TO OLD SCALAR ROUTINE FFT9
!     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM
!     IS GIVEN BY COOLEY, LEWIS AND WELCH (J. SOUND VIB., VOL. 12
!     (1970), 315-337)
!
!     A IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA
!     WORK IS AN AREA OF SIZE (N+1)*LOT
!     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2
!     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
!         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!     N IS THE LENGTH OF THE DATA VECTORS
!     LOT IS THE NUMBER OF DATA VECTORS
!     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
!
!     ORDERING OF COEFFICIENTS:
!         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!
!     ORDERING OF DATA:
!         X(N-1),X(0),X(1),X(2),...,X(N),X(0)
!         I.E. EXPLICIT CYCLIC CONTINUITY; (N+2) LOCATIONS REQUIRED
!
!     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
!     PARALLEL
!
!     *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER
!
!     DEFINITION OF TRANSFORMS:
!     -------------------------
!
!     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!
!     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
!               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
!
!-----------------------------------------------------------------------
!
      use shr_kind_mod,    only: r8 => shr_kind_r8
      implicit none
!
!------------------------------Arguments--------------------------------
!
      integer IFAX(13), inc, jump, n, lot, isign
      real(r8) A(LOT* (N+2) ), WORK(LOT*(N+1)), TRIGS(3*N/2+1)
!
!---------------------------Local variables-----------------------------
!
      integer i, j, k, l, m, ia, ib, la, ink, nh, nx, nfax
      integer ibase, jbase, igo
!
!-----------------------------------------------------------------------
!
      NFAX=IFAX(1)
      NX=N+1
      NH=N/2
      INK=INC+INC
      IF (ISIGN.EQ.+1) GO TO 30
!
!     IF NECESSARY, TRANSFER DATA TO WORK AREA
      IGO=50
      IF (MOD(NFAX,2).EQ.1) GOTO 40
      IBASE=INC+1
      JBASE=1
      DO 20 L=1,LOT
      I=IBASE
      J=JBASE
!cdir nodep
!DIR$ CONCURRENT
      DO 10 M=1,N
      WORK(J)=A(I)
      I=I+INC
      J=J+1
   10 CONTINUE
      IBASE=IBASE+JUMP
      JBASE=JBASE+NX
   20 CONTINUE
!
      IGO=60
      GO TO 40
!
!     PREPROCESSING (ISIGN=+1)
!     ------------------------
!
   30 CONTINUE
      CALL FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      IGO=60
!
!     COMPLEX TRANSFORM
!     -----------------
!
   40 CONTINUE
      IA=INC+1
      LA=1
      DO 80 K=1,NFAX
      IF (IGO.EQ.60) GO TO 60
   50 CONTINUE
      CALL VPASSM(A(IA),A(IA+INC),WORK(1),WORK(2),TRIGS, &
         INK,2,JUMP,NX,LOT,NH,IFAX(K+1),LA)
      IGO=60
      GO TO 70
   60 CONTINUE
      CALL VPASSM(WORK(1),WORK(2),A(IA),A(IA+INC),TRIGS, &
          2,INK,NX,JUMP,LOT,NH,IFAX(K+1),LA)
      IGO=50
   70 CONTINUE
      LA=LA*IFAX(K+1)
   80 CONTINUE
!
      IF (ISIGN.EQ.-1) GO TO 130
!
!     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      IF (MOD(NFAX,2).EQ.1) GO TO 110
      IBASE=1
      JBASE=IA
      DO 100 L=1,LOT
      I=IBASE
      J=JBASE
!cdir nodep
!DIR$ CONCURRENT
      DO 90 M=1,N
      A(J)=WORK(I)
      I=I+1
      J=J+INC
   90 CONTINUE
      IBASE=IBASE+NX
      JBASE=JBASE+JUMP
  100 CONTINUE
!
!     FILL IN CYCLIC BOUNDARY POINTS
  110 CONTINUE
      IA=1
      IB=N*INC+1
!cdir nodep
!DIR$ CONCURRENT
      DO 120 L=1,LOT
      A(IA)=A(IB)
      A(IB+INC)=A(IA+INC)
      IA=IA+JUMP
      IB=IB+JUMP
  120 CONTINUE
      GO TO 140
!
!     POSTPROCESSING (ISIGN=-1):
!     --------------------------
!
  130 CONTINUE
      CALL FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
!
  140 CONTINUE
      RETURN
      END SUBROUTINE FFT99

!================================================================================================

      SUBROUTINE FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
!
!-----------------------------------------------------------------------
!     SUBROUTINE FFT99A - PREPROCESSING STEP FOR FFT99, ISIGN=+1
!     (SPECTRAL TO GRIDPOINT TRANSFORM)
!-----------------------------------------------------------------------
!
      use shr_kind_mod,    only: r8 => shr_kind_r8
      implicit none
!
!------------------------------Arguments--------------------------------
!
      integer inc, jump, n, lot
      real(r8) A(*), WORK(*), TRIGS(*)
!
!---------------------------Local variables-----------------------------
!
      integer iabase, ibbase, jabase, jbbase, ia, ib, ink
      integer ja, jb, k, l, nh, nx
      real(r8) c, s
!
!-----------------------------------------------------------------------
!
      NH=N/2
      NX=N+1
      INK=INC+INC
!
!     A(0) AND A(N/2)
      IA=1
      IB=N*INC+1
      JA=1
      JB=2
!cdir nodep
!DIR$ CONCURRENT
      DO 10 L=1,LOT
      WORK(JA)=A(IA)+A(IB)
      WORK(JB)=A(IA)-A(IB)
      IA=IA+JUMP
      IB=IB+JUMP
      JA=JA+NX
      JB=JB+NX
   10 CONTINUE
!
!     REMAINING WAVENUMBERS
      IABASE=2*INC+1
      IBBASE=(N-2)*INC+1
      JABASE=3
      JBBASE=N-1
!
      DO 30 K=3,NH,2
      IA=IABASE
      IB=IBBASE
      JA=JABASE
      JB=JBBASE
      C=TRIGS(N+K)
      S=TRIGS(N+K+1)
!cdir nodep
!DIR$ CONCURRENT
      DO 20 L=1,LOT
      WORK(JA)=(A(IA)+A(IB))- &
          (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
      WORK(JB)=(A(IA)+A(IB))+ &
          (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
      WORK(JA+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))+ &
          (A(IA+INC)-A(IB+INC))
      WORK(JB+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))- &
          (A(IA+INC)-A(IB+INC))
      IA=IA+JUMP
      IB=IB+JUMP
      JA=JA+NX
      JB=JB+NX
   20 CONTINUE
      IABASE=IABASE+INK
      IBBASE=IBBASE-INK
      JABASE=JABASE+2
      JBBASE=JBBASE-2
   30 CONTINUE
!
      IF (IABASE.NE.IBBASE) GO TO 50
!     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE
!cdir nodep
!DIR$ CONCURRENT
      DO 40 L=1,LOT
      WORK(JA)=2.0_r8*A(IA)
      WORK(JA+1)=-2.0_r8*A(IA+INC)
      IA=IA+JUMP
      JA=JA+NX
   40 CONTINUE
!
   50 CONTINUE
      RETURN
      END SUBROUTINE FFT99A

!================================================================================================

      SUBROUTINE FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
!
!-----------------------------------------------------------------------
!     SUBROUTINE FFT99B - POSTPROCESSING STEP FOR FFT99, ISIGN=-1
!     (GRIDPOINT TO SPECTRAL TRANSFORM)
!-----------------------------------------------------------------------
!
      use shr_kind_mod,    only: r8 => shr_kind_r8
      implicit none
!
!------------------------------Arguments--------------------------------
!
      integer inc, jump, n, lot
      real(r8) WORK(*), A(*), TRIGS(*)
!
!---------------------------Local variables-----------------------------
!
      integer iabase, ibbase, jabase, jbbase, ia, ib, ink
      integer ja, jb, k, l, nh, nx
      real(r8) scale, s, c
!
!-----------------------------------------------------------------------
!
      NH=N/2
      NX=N+1
      INK=INC+INC
!
!     A(0) AND A(N/2)
      SCALE=1.0_r8/real(N,r8)
      IA=1
      IB=2
      JA=1
      JB=N*INC+1
!cdir nodep
!DIR$ CONCURRENT
      DO 10 L=1,LOT
      A(JA)=SCALE*(WORK(IA)+WORK(IB))
      A(JB)=SCALE*(WORK(IA)-WORK(IB))
      A(JA+INC)=0.0_r8
      A(JB+INC)=0.0_r8
      IA=IA+NX
      IB=IB+NX
      JA=JA+JUMP
      JB=JB+JUMP
   10 CONTINUE
!
!     REMAINING WAVENUMBERS
      SCALE=0.5_r8*SCALE
      IABASE=3
      IBBASE=N-1
      JABASE=2*INC+1
      JBBASE=(N-2)*INC+1
!
      DO 30 K=3,NH,2
      IA=IABASE
      IB=IBBASE
      JA=JABASE
      JB=JBBASE
      C=TRIGS(N+K)
      S=TRIGS(N+K+1)
!cdir nodep
!DIR$ CONCURRENT
      DO 20 L=1,LOT
      A(JA)=SCALE*((WORK(IA)+WORK(IB)) &
         +(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
      A(JB)=SCALE*((WORK(IA)+WORK(IB)) &
         -(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
      A(JA+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1))) &
          +(WORK(IB+1)-WORK(IA+1)))
      A(JB+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1))) &
          -(WORK(IB+1)-WORK(IA+1)))
      IA=IA+NX
      IB=IB+NX
      JA=JA+JUMP
      JB=JB+JUMP
   20 CONTINUE
      IABASE=IABASE+2
      IBBASE=IBBASE-2
      JABASE=JABASE+INK
      JBBASE=JBBASE-INK
   30 CONTINUE
!
      IF (IABASE.NE.IBBASE) GO TO 50
!     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE
      SCALE=2.0_r8*SCALE
!cdir nodep
!DIR$ CONCURRENT
      DO 40 L=1,LOT
      A(JA)=SCALE*WORK(IA)
      A(JA+INC)=-SCALE*WORK(IA+1)
      IA=IA+NX
      JA=JA+JUMP
   40 CONTINUE
!
   50 CONTINUE
      RETURN
      END SUBROUTINE FFT99B

!================================================================================================

      SUBROUTINE FFT991(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
!
!-----------------------------------------------------------------------
!     SUBROUTINE "FFT991" - MULTIPLE REAL/HALF-COMPLEX PERIODIC
!     FAST FOURIER TRANSFORM
!
!     SAME AS FFT99 EXCEPT THAT ORDERING OF DATA CORRESPONDS TO
!     THAT IN MRFFT2
!
!     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM
!     IS GIVEN BY COOLEY, LEWIS AND WELCH (J. SOUND VIB., VOL. 12
!     (1970), 315-337)
!
!     A IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA
!     WORK IS AN AREA OF SIZE (N+1)*LOT
!     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2
!     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
!         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!     N IS THE LENGTH OF THE DATA VECTORS
!     LOT IS THE NUMBER OF DATA VECTORS
!     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
!
!     ORDERING OF COEFFICIENTS:
!         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!
!     ORDERING OF DATA:
!         X(0),X(1),X(2),...,X(N-1)
!
!     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
!     PARALLEL
!
!     *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER
!
!     DEFINITION OF TRANSFORMS:
!     -------------------------
!
!     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!
!     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
!               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
!-----------------------------------------------------------------------
!
      use shr_kind_mod,    only: r8 => shr_kind_r8
      implicit none
!
!------------------------------Arguments--------------------------------
!
      integer IFAX(13), inc, jump, n, lot, isign
      real(r8) A(*), WORK(*), TRIGS(*)
!
!---------------------------Local variables-----------------------------
!
      integer ibase, jbase, i, j, ia, ib, igo, ink, k
      integer l, la, m, nh, nfax, nx
!
!-----------------------------------------------------------------------
!
      NFAX=IFAX(1)
      NX=N+1
      NH=N/2
      INK=INC+INC
      IF (ISIGN.EQ.+1) GO TO 30
!
!     IF NECESSARY, TRANSFER DATA TO WORK AREA
      IGO=50
      IF (MOD(NFAX,2).EQ.1) GOTO 40
      IBASE=1
      JBASE=1
      DO 20 L=1,LOT
      I=IBASE
      J=JBASE
!cdir nodep
!DIR$ CONCURRENT
      DO 10 M=1,N
      WORK(J)=A(I)
      I=I+INC
      J=J+1
   10 CONTINUE
      IBASE=IBASE+JUMP
      JBASE=JBASE+NX
   20 CONTINUE
!
      IGO=60
      GO TO 40
!
!     PREPROCESSING (ISIGN=+1)
!     ------------------------
!
   30 CONTINUE
      CALL FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      IGO=60
!
!     COMPLEX TRANSFORM
!     -----------------
!
   40 CONTINUE
      IA=1
      LA=1
      DO 80 K=1,NFAX
      IF (IGO.EQ.60) GO TO 60
   50 CONTINUE
      CALL VPASSM(A(IA),A(IA+INC),WORK(1),WORK(2),TRIGS, &
         INK,2,JUMP,NX,LOT,NH,IFAX(K+1),LA)
      IGO=60
      GO TO 70
   60 CONTINUE
      CALL VPASSM(WORK(1),WORK(2),A(IA),A(IA+INC),TRIGS, &
          2,INK,NX,JUMP,LOT,NH,IFAX(K+1),LA)
      IGO=50
   70 CONTINUE
      LA=LA*IFAX(K+1)
   80 CONTINUE
!
      IF (ISIGN.EQ.-1) GO TO 130
!
!     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      IF (MOD(NFAX,2).EQ.1) GO TO 110
      IBASE=1
      JBASE=1
      DO 100 L=1,LOT
      I=IBASE
      J=JBASE
!cdir nodep
!DIR$ CONCURRENT
      DO 90 M=1,N
      A(J)=WORK(I)
      I=I+1
      J=J+INC
   90 CONTINUE
      IBASE=IBASE+NX
      JBASE=JBASE+JUMP
  100 CONTINUE
!
!     FILL IN ZEROS AT END
  110 CONTINUE
      IB=N*INC+1
!cdir nodep
!DIR$ CONCURRENT
      DO 120 L=1,LOT
      A(IB)=0.0_r8
      A(IB+INC)=0.0_r8
      IB=IB+JUMP
  120 CONTINUE
      GO TO 140
!
!     POSTPROCESSING (ISIGN=-1):
!     --------------------------
!
  130 CONTINUE
      CALL FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
!
  140 CONTINUE
      RETURN
      END SUBROUTINE FFT991

!================================================================================================

      SUBROUTINE SET99 (TRIGS, IFAX, N)
!
!-----------------------------------------------------------------------
!
      use shr_kind_mod,    only: r8 => shr_kind_r8
      use cam_logfile,     only: iulog
      implicit none
!
!------------------------------Arguments--------------------------------
!
      integer n, IFAX(13)
      real(r8) TRIGS(*)
!
!---------------------------Local variables-----------------------------
!
      integer i
!
! MODE 3 IS USED FOR REAL/HALF-COMPLEX TRANSFORMS.  IT IS POSSIBLE
! TO DO COMPLEX/COMPLEX TRANSFORMS WITH OTHER VALUES OF MODE, BUT
! DOCUMENTATION OF THE DETAILS WERE NOT AVAILABLE WHEN THIS ROUTINE
! WAS WRITTEN.
!
      integer mode
      DATA MODE /3/
!
!-----------------------------------------------------------------------
!
      CALL FAX (IFAX, N, MODE)
      I = IFAX(1)
      IF (IFAX(I+1) .GT. 5 .OR. N .LE. 4) IFAX(1) = -99
      IF (IFAX(1) .LE. 0 ) THEN 
        write(iulog,*) ' SET99 -- INVALID N'
        STOP 'SET99'
      ENDIF
      CALL FFTRIG (TRIGS, N, MODE)
      RETURN
      END SUBROUTINE SET99

!================================================================================================

      SUBROUTINE FAX(IFAX,N,MODE)
!
!-----------------------------------------------------------------------
!
      use shr_kind_mod,    only: r8 => shr_kind_r8
      implicit none
!
!------------------------------Arguments--------------------------------
!
      integer IFAX(10), n, mode   
!
!---------------------------Local variables-----------------------------
!
      integer ii, nfax, inc, item, i, istop, l, k, nn
!
!-----------------------------------------------------------------------
!
      NN=N
      IF (IABS(MODE).EQ.1) GO TO 10
      IF (IABS(MODE).EQ.8) GO TO 10
      NN=N/2
      IF ((NN+NN).EQ.N) GO TO 10
      IFAX(1)=-99
      RETURN
   10 K=1
!     TEST FOR FACTORS OF 4
   20 IF (MOD(NN,4).NE.0) GO TO 30
      K=K+1
      IFAX(K)=4
      NN=NN/4
      IF (NN.EQ.1) GO TO 80
      GO TO 20
!     TEST FOR EXTRA FACTOR OF 2
   30 IF (MOD(NN,2).NE.0) GO TO 40
      K=K+1
      IFAX(K)=2
      NN=NN/2
      IF (NN.EQ.1) GO TO 80
!     TEST FOR FACTORS OF 3
   40 IF (MOD(NN,3).NE.0) GO TO 50
      K=K+1
      IFAX(K)=3
      NN=NN/3
      IF (NN.EQ.1) GO TO 80
      GO TO 40
!     NOW FIND REMAINING FACTORS
   50 L=5
      INC=2
!     INC ALTERNATELY TAKES ON VALUES 2 AND 4
   60 IF (MOD(NN,L).NE.0) GO TO 70
      K=K+1
      IFAX(K)=L
      NN=NN/L
      IF (NN.EQ.1) GO TO 80
      GO TO 60
   70 L=L+INC
      INC=6-INC
      GO TO 60
   80 IFAX(1)=K-1
!     IFAX(1) CONTAINS NUMBER OF FACTORS
      NFAX=IFAX(1)
!     SORT FACTORS INTO ASCENDING ORDER
      IF (NFAX.EQ.1) GO TO 110
      DO 100 II=2,NFAX
      ISTOP=NFAX+2-II
      DO 90 I=2,ISTOP
      IF (IFAX(I+1).GE.IFAX(I)) GO TO 90
      ITEM=IFAX(I)
      IFAX(I)=IFAX(I+1)
      IFAX(I+1)=ITEM
   90 CONTINUE
  100 CONTINUE
  110 CONTINUE
      RETURN
      END SUBROUTINE FAX

!================================================================================================

      SUBROUTINE FFTRIG(TRIGS,N,MODE)
!
!-----------------------------------------------------------------------
!
      use shr_kind_mod,    only: r8 => shr_kind_r8
      implicit none
!
!------------------------------Arguments--------------------------------
!
      integer n, mode
      real(r8) TRIGS(*)
!
!---------------------------Local variables-----------------------------
!
      integer i, l, la, nh, imode, nn
      real(r8) del, pi, angle
!
!-----------------------------------------------------------------------
!
      PI=2.0_r8*ASIN(1.0_r8)
      IMODE=IABS(MODE)
      NN=N
      IF (IMODE.GT.1.AND.IMODE.LT.6) NN=N/2
      DEL=(PI+PI)/real(NN,r8)
      L=NN+NN
      DO 10 I=1,L,2
      ANGLE=0.5_r8*real(I-1,r8)*DEL
      TRIGS(I)=COS(ANGLE)
      TRIGS(I+1)=SIN(ANGLE)
   10 CONTINUE
      IF (IMODE.EQ.1) RETURN
      IF (IMODE.EQ.8) RETURN
      DEL=0.5_r8*DEL
      NH=(NN+1)/2
      L=NH+NH
      LA=NN+NN
      DO 20 I=1,L,2
      ANGLE=0.5_r8*real(I-1,r8)*DEL
      TRIGS(LA+I)=COS(ANGLE)
      TRIGS(LA+I+1)=SIN(ANGLE)
   20 CONTINUE
      IF (IMODE.LE.3) RETURN
      DEL=0.5_r8*DEL
      LA=LA+NN
      IF (MODE.EQ.5) GO TO 40
      DO 30 I=2,NN
      ANGLE=real(I-1,r8)*DEL
      TRIGS(LA+I)=2.0_r8*SIN(ANGLE)
   30 CONTINUE
      RETURN
   40 CONTINUE
      DEL=0.5_r8*DEL
      DO 50 I=2,N
      ANGLE=real(I-1,r8)*DEL
      TRIGS(LA+I)=SIN(ANGLE)
   50 CONTINUE
      RETURN
      END SUBROUTINE FFTRIG

!================================================================================================

      SUBROUTINE VPASSM(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,LA)
!
!-----------------------------------------------------------------------
!     SUBROUTINE "VPASSM" - MULTIPLE VERSION OF "VPASSA"
!     PERFORMS ONE PASS THROUGH DATA
!     AS PART OF MULTIPLE COMPLEX FFT ROUTINE
!     A IS FIRST REAL INPUT VECTOR
!     B IS FIRST IMAGINARY INPUT VECTOR
!     C IS FIRST REAL OUTPUT VECTOR
!     D IS FIRST IMAGINARY OUTPUT VECTOR
!     TRIGS IS PRECALCULATED TABLE OF SINES " COSINES
!     INC1 IS ADDRESSING INCREMENT FOR A AND B
!     INC2 IS ADDRESSING INCREMENT FOR C AND D
!     INC3 IS ADDRESSING INCREMENT BETWEEN A"S & B"S
!     INC4 IS ADDRESSING INCREMENT BETWEEN C"S & D"S
!     LOT IS THE NUMBER OF VECTORS
!     N IS LENGTH OF VECTORS
!     IFAC IS CURRENT FACTOR OF N
!     LA IS PRODUCT OF PREVIOUS FACTORS
!-----------------------------------------------------------------------
!
      use shr_kind_mod,    only: r8 => shr_kind_r8
      implicit none
!
!------------------------------Arguments--------------------------------
!
      integer inc1, inc2, inc3, inc4, lot, n, ifac, la
      real(r8) A(*),B(*),C(*),D(*),TRIGS(*)
!
!---------------------------Local variables-----------------------------
!
      integer ie, je, ke, kd, ib, ja, ia, i, l, jb, igo, jink 
      integer iink, m, jbase, ibase, jump, j, kc, jc, jd, id
      integer ic, k, la1, ijk, kb

      real(r8) s3, c3, s4, c4, c2, s2, s1, c1

      real(r8) sin36, cos36, sin72, cos72, sin60
      DATA SIN36/0.587785252292473_r8/,COS36/0.809016994374947_r8/, &
           SIN72/0.951056516295154_r8/,COS72/0.309016994374947_r8/, &
           SIN60/0.866025403784437_r8/
!
!-----------------------------------------------------------------------
!
      M=N/IFAC
      IINK=M*INC1
      JINK=LA*INC2
      JUMP=(IFAC-1)*JINK
      IBASE=0
      JBASE=0
      IGO=IFAC-1
      IF (IGO.GT.4) RETURN
      GO TO (10,50,90,130),IGO
!
!     CODING FOR FACTOR 2
!
   10 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      DO 20 L=1,LA
      I=IBASE
      J=JBASE
!cdir nodep
!DIR$ CONCURRENT
      DO 15 IJK=1,LOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)+B(IB+I)
      C(JB+J)=A(IA+I)-A(IB+I)
      D(JB+J)=B(IA+I)-B(IB+I)
      I=I+INC3
      J=J+INC4
   15 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   20 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 40 K=LA1,M,LA
      KB=K+K-2
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      DO 30 L=1,LA
      I=IBASE
      J=JBASE
!cdir nodep
!DIR$ CONCURRENT
      DO 25 IJK=1,LOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)+B(IB+I)
      C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)-B(IB+I))
      D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)-B(IB+I))
      I=I+INC3
      J=J+INC4
   25 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   30 CONTINUE
      JBASE=JBASE+JUMP
   40 CONTINUE
      RETURN
!
!     CODING FOR FACTOR 3
!
   50 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      DO 60 L=1,LA
      I=IBASE
      J=JBASE
!cdir nodep
!DIR$ CONCURRENT
      DO 55 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
      C(JB+J)=(A(IA+I)-0.5_r8*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I)))
      C(JC+J)=(A(IA+I)-0.5_r8*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I)))
      D(JB+J)=(B(IA+I)-0.5_r8*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I)))
      D(JC+J)=(B(IA+I)-0.5_r8*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I)))
      I=I+INC3
      J=J+INC4
   55 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   60 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 80 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      DO 70 L=1,LA
      I=IBASE
      J=JBASE
!cdir nodep
!DIR$ CONCURRENT
      DO 65 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
      C(JB+J)= &
          C1*((A(IA+I)-0.5_r8*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I)))) &
         -S1*((B(IA+I)-0.5_r8*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
      D(JB+J)= &
          S1*((A(IA+I)-0.5_r8*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I)))) &
         +C1*((B(IA+I)-0.5_r8*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
      C(JC+J)= &
          C2*((A(IA+I)-0.5_r8*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I)))) &
         -S2*((B(IA+I)-0.5_r8*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
      D(JC+J)= &
          S2*((A(IA+I)-0.5_r8*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I)))) &
         +C2*((B(IA+I)-0.5_r8*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
      I=I+INC3
      J=J+INC4
   65 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   70 CONTINUE
      JBASE=JBASE+JUMP
   80 CONTINUE
      RETURN
!
!     CODING FOR FACTOR 4
!
   90 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      DO 100 L=1,LA
      I=IBASE
      J=JBASE
!cdir nodep
!DIR$ CONCURRENT
      DO 95 IJK=1,LOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      C(JC+J)=(A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
      D(JC+J)=(B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I))
      C(JB+J)=(A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I))
      C(JD+J)=(A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I))
      D(JB+J)=(B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I))
      D(JD+J)=(B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I))
      I=I+INC3
      J=J+INC4
   95 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  100 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 120 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      KD=KC+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      DO 110 L=1,LA
      I=IBASE
      J=JBASE
!cdir nodep
!DIR$ CONCURRENT
      DO 105 IJK=1,LOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
      C(JC+J)= &
          C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))) &
         -S2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
      D(JC+J)= &
          S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))) &
         +C2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
      C(JB+J)= &
          C1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I))) &
         -S1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
      D(JB+J)= &
          S1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I))) &
         +C1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
      C(JD+J)= &
          C3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I))) &
         -S3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
      D(JD+J)= &
          S3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I))) &
         +C3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  105 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  110 CONTINUE
      JBASE=JBASE+JUMP
  120 CONTINUE
      RETURN
!
!     CODING FOR FACTOR 5
!
  130 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      IE=ID+IINK
      JE=JD+JINK
      DO 140 L=1,LA
      I=IBASE
      J=JBASE
!cdir nodep
!DIR$ CONCURRENT
      DO 135 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
      C(JB+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I))) &
        -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      C(JE+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I))) &
        +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      D(JB+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I))) &
        +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      D(JE+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I))) &
        -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      C(JC+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I))) &
        -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      C(JD+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I))) &
        +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      D(JC+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I))) &
        +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
      D(JD+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I))) &
        -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  135 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  140 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 160 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      KD=KC+KB
      KE=KD+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      C4=TRIGS(KE+1)
      S4=TRIGS(KE+2)
      DO 150 L=1,LA
      I=IBASE
      J=JBASE
!cdir nodep
!DIR$ CONCURRENT
      DO 145 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
      C(JB+J)= &
          C1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I))) &
            -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))) &
         -S1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I))) &
            +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      D(JB+J)= &
          S1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I))) &
            -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))) &
         +C1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I))) &
            +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      C(JE+J)= &
          C4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I))) &
            +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))) &
         -S4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I))) &
            -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      D(JE+J)= &
          S4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I))) &
            +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))) &
         +C4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I))) &
            -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      C(JC+J)= &
          C2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I))) &
            -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))) &
         -S2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I))) &
            +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      D(JC+J)= &
          S2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I))) &
            -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))) &
         +C2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I))) &
            +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      C(JD+J)= &
          C3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I))) &
            +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))) &
         -S3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I))) &
            -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      D(JD+J)= &
          S3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I))) &
            +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))) &
         +C3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I))) &
            -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      I=I+INC3
      J=J+INC4
  145 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  150 CONTINUE
      JBASE=JBASE+JUMP
  160 CONTINUE
      RETURN
      END SUBROUTINE VPASSM

!===============================================================================

 
