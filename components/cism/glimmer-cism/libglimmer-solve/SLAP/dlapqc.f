      PROGRAM DLAPQC
C***BEGIN PROLOGUE  DLAPQC
C***SUBSIDIARY
C***PURPOSE  Driver for testing SLATEC Sparse Linear Algebra Package
C            (SLAP) Version 2.0.
C***LIBRARY   SLATEC(SLAP)
C***CATEGORY  D2A4, D2B4
C***TYPE      SINGLE (DLAPQC-S)
C***KEYWORDS  QUICK CHECK DRIVER, SLAP
C***AUTHOR  Mark K. Seager (LLNL)
C             seager@lll-crg.llnl.gov
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550
C             (415)423-3141
C***DESCRIPTION
C
C *Usage:
C     One input data record is required
C         READ (LIN,990) KPRINT
C     999 FORMAT (I1)
C
C *Arguments:
C     KPRINT = 0  Quick checks - No printing.
C                 Driver       - Short pass or fail message printed.
C              1  Quick checks - No message printed for passed tests,
C                                short message printed for failed tests.
C                 Driver       - Short pass or fail message printed.
C              2  Quick checks - Print short message for passed tests,
C                                fuller information for failed tests.
C                 Driver       - Pass or fail message printed.
C              3  Quick checks - Print complete quick check results.
C                 Driver       - Pass or fail message printed.
C              4  Quick checks - Print complete quick check results.  
C                                Prints matricies, etc.  Very verbose. 
C                 Driver       - Pass or fail message printed.
C
C *Description:
C                                 The
C                    Sparse Linear Algebra Package
C
C                @@@@@@@  @            @@@    @@@@@@@@
C               @       @ @           @   @   @       @
C               @         @          @     @  @       @
C                @@@@@@@  @         @       @ @@@@@@@@
C                       @ @         @@@@@@@@@ @
C               @       @ @         @       @ @        
C                @@@@@@@  @@@@@@@@@ @       @ @        
C                                        
C      @       @                            @@@@@@@        @@@@@
C      @       @                           @       @      @    @@ 
C      @       @  @@@@@@@  @ @@                    @     @    @  @
C      @       @ @       @ @@  @             @@@@@@      @   @   @
C       @     @  @@@@@@@@@ @                @            @  @    @
C        @   @   @         @               @         @@@  @@    @ 
C         @@@     @@@@@@@  @               @@@@@@@@@ @@@   @@@@@  
C
C----------------------------------------------------------------------
C                              Written By
C                        Mark K. Seager (LLNL)
C                   Lawrence Livermore National Lab.
C                          PO Box 808, L-300
C                         Livermore, CA 94550
C                            (415) 423-3141
C                       seager@lll-crg.llnl.gov
C----------------------------------------------------------------------
C         This is a SLATEC Quick Checks program to test the *SLAP* 
C         Version 2.0 package.  It generates a "random" matrix (See 
C         DRMGEN) and then runs all the various methods with all the 
C         various preconditoners and all the various stop tests.
C
C         It is assumed that the test is being run interactively and 
C         that STDIN (STANDARD INPUT) is Fortran I/O unit I1MACH(1) 
C         and STDOUT (STANDARD OUTPUT) is unit I1MACH(2).
C
C         *************************************************************
C         **** WARNING !!! WARNING !!! WARNING !!! WARNING !!! WARNING
C         *************************************************************
C         **** THIS PROGRAM WILL NOT FUNCTION PROPERLY IF THE FORTRAN
C         **** I/O UNITS I1MACH(1) and I1MACH(2) are not connected
C         **** to the program for I/O.
C         *************************************************************
C
C         All errors in the driver are handled with the SLATEC error
C         handler (revision date 851111).
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DS2Y, DCPPLT, DSJAC, DSGS, DSILUR, DSDCG, DSICCG,
C                    DSDCGN, DSLUCN, DSDBCG, DSLUBC, DSDCGS, DSLUCS, 
C                    DSDOMN, DSLUOM, DSDCMR, DSLUCM
C***REVISION HISTORY  (YYMMDD)
C   880601  DATE WRITTEN
C   881213  Revised to meet the new SLATEC prologue standards.
C***END PROLOGUE  DLAPQC
      PARAMETER(MAXN=441, MXNELT=50000, MAXIW=50000, MAXRW=50000)
C$$$      PARAMETER(MAXN=25, MXNELT=50000, MAXIW=50000, MAXRW=50000)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*72     MESG
      DOUBLE PRECISION A(MXNELT), F(MAXN), XITER(MAXN), RWORK(MAXRW)
      INTEGER          IA(MXNELT), JA(MXNELT), IWORK(MAXIW)
      COMMON /SOLBLK/ SOLN(MAXN)
C
C         The following lines is for the braindammaged Sun FPE handler.
C
C$$$      integer oldmode, fpmode
C$$$      oldmode = fpmode( 62464 )
C
C   READ KPRINT PARAMETER
C
C***FIRST EXECUTABLE STATEMENT  DLAPQC
C
      ISTDI = I1MACH(1)
      ISTDO = I1MACH(2)
      NFAIL = 0
      READ(ISTDI,990) KPRINT
  990 FORMAT(I1)
      CALL XSETUN(LUN)
      IF( KPRINT.LE.1 ) THEN
         CALL XSETF(0)
      ELSE
         CALL XSETF(1)
      ENDIF
      CALL XERMAX(1000)
C
C         Maximum problem sizes.
C
      NELTMX = MXNELT
      NMAX   = MAXN
      LENIW  = MAXIW
      LENW   = MAXRW
C
C        Set some input data.
C
      N      = NMAX
      ITMAX  = N
      IOUT   = KPRINT
      FACTOR = 1.2
      IF( IOUT.LT.3 ) THEN
         IUNIT = 0
      ELSE
         IUNIT = ISTDO
      ENDIF
C
C         Set the Error tolerance to depend on the machine epsilon.
C
      TOL = MAX(1.0D3*D1MACH(3),1.0D-6)
C         
C         Test routines using various convergence criteria.
C
      DO 10 KASE = 3, 3
         IF(KASE .EQ. 1 .OR. KASE .EQ. 2) ITOL = KASE
         IF(KASE .EQ. 3) ITOL = 11
C         
C         Test routines using nonsymmetric (ISYM=0) and symmetric
C         storage (ISYM=1).  For ISYM=0 a really non-symmetric matrix
C         is generated.  The amount of non-symmetry is controlled by
C         user.
C
         DO 20 ISYM = 0, 1
C
C         Set up a random matrix.
C
            CALL DRMGEN( NELTMX, FACTOR, IERR, N, NELT, 
     $           ISYM, IA, JA, A, F, SOLN, RWORK, IWORK, IWORK(N+1) )
            IF( IERR.NE.0 ) THEN
               MESG = 'DLAPQC -- Fatal error (i1) generating '//
     $              '*RANDOM* Matrix.'
               CALL XERRWV( MESG,LEN(MESG),IERR,2,1,IERR,0,
     $              0,0.0,0.0 )
            ENDIF
            IF( ISYM.EQ.0 ) THEN
               DENS = FLOAT(NELT)/FLOAT(N*N)
            ELSE
               DENS = FLOAT(2*NELT)/FLOAT(N*N)
            ENDIF
            IF( IOUT.GE.2 ) THEN
              WRITE(ISTDO,1020) N, NELT, DENS
              WRITE(ISTDO,1030) TOL
            ENDIF
C         
C         Convert to the SLAP-Column format and
C         write out matrix in SLAP-Column format, if desired.
C
            CALL DS2Y( N, NELT, IA, JA, A, ISYM )
            IF( IOUT.GE.4 ) THEN
               WRITE(ISTDO,1040) (K,IA(K),JA(K),A(K),K=1,NELT)
               CALL DCPPLT( N, NELT, IA, JA, A, ISYM, ISTDO )
            ENDIF
C
C**********************************************************************
C                    BEGINING OF SLAP QUICK TESTS
C**********************************************************************
C
C         * * * * * *   DSJAC   * * * * * *
C
            IF( IOUT.GE.3 ) THEN
              WRITE(ISTDO,1000) 'DSJAC ', ITOL, ISYM
            ENDIF
            CALL DFILL( N, XITER, 0.0D0 )
C
            CALL DSJAC(N, F, XITER, NELT, IA, JA, A, ISYM,
     $           ITOL, TOL, 2*ITMAX, ITER, ERR, IERR, IUNIT,
     $           RWORK, LENW, IWORK, LENIW )
C
            CALL DUTERR( 'DSJAC ',IERR,IOUT,NFAIL,ISTDO,ITER,ERR )
C         
C         * * * * *  DSGS  * * * * *
C         
            IF( IOUT.GE.3 ) THEN
              WRITE(ISTDO,1000) 'DSGS  ',ITOL,ISYM
            ENDIF
            CALL DFILL( N, XITER, 0.0D0 )
C
            CALL DSGS(N, F, XITER, NELT, IA, JA, A, ISYM,
     $           ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
     $           RWORK, LENW, IWORK, LENIW )
C
            CALL DUTERR( 'DSGS  ',IERR,IOUT,NFAIL,ISTDO,ITER,ERR )
C         
C         * * * * * *   DSILUR   * * * * * *
C         
            IF( IOUT.GE.3 ) THEN
              WRITE(ISTDO,1000) 'DSILUR',ITOL,ISYM
            ENDIF
            CALL DFILL( N, XITER, 0.0D0 )
C
            CALL DSILUR(N, F, XITER, NELT, IA, JA, A, ISYM,
     $           ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
     $           RWORK, LENW, IWORK, LENIW )
C
            CALL DUTERR( 'DSILUR',IERR,IOUT,NFAIL,ISTDO,ITER,ERR )
C         
C         * * * * * *   DSDCG    * * * * * *
C         
            IF( ISYM.EQ.1 ) THEN
               IF( IOUT.GE.3 ) THEN
                  WRITE(ISTDO,1000) 'DSDCG',ITOL,ISYM
               ENDIF
               CALL DFILL( N, XITER, 0.0D0 )
C
               CALL DSDCG(N, F, XITER, NELT, IA, JA, A, ISYM, 
     $              ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
     $              RWORK, LENW, IWORK, LENIW )
C
               CALL DUTERR( 'DSDCG ',IERR,IOUT,NFAIL,ISTDO,ITER,ERR )
            ENDIF
C         
C         * * * * * *    DSICCG    * * * * * *
C         
            IF( ISYM.EQ.1 ) THEN
               IF( IOUT.GE.3 ) THEN
                  WRITE(ISTDO,1000) 'DSICCG',ITOL,ISYM
               ENDIF
               CALL DFILL( N, XITER, 0.0D0 )
C
               CALL DSICCG(N, F, XITER, NELT, IA, JA, A, ISYM, 
     $              ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, 
     $              LENW, IWORK, LENIW )
C
               CALL DUTERR( 'DSICCG',IERR,IOUT,NFAIL,ISTDO,ITER,ERR )
            ENDIF
C         
C         * * * * * *    DSDCGN   * * * * * *
C         
            IF( IOUT.GE.3 ) THEN
               WRITE(ISTDO,1000) 'DSDCGN',ITOL,ISYM
            ENDIF
            CALL DFILL( N, XITER, 0.0D0 )
C
            CALL DSDCGN(N, F, XITER, NELT, IA, JA, A, ISYM, ITOL,
     $           TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW,
     $           IWORK, LENIW )
C
            CALL DUTERR( 'DSDCGN',IERR,IOUT,NFAIL,ISTDO,ITER,ERR )
C         
C         * * * * * *   DSLUCN   * * * * * *
C
            IF( IOUT.GE.3 ) THEN 
               WRITE(ISTDO,1000) 'DSLUCN',ITOL,ISYM
            ENDIF
            CALL DFILL( N, XITER, 0.0D0 )
C
            CALL DSLUCN(N, F, XITER, NELT, IA, JA, A, ISYM, ITOL,
     $           TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW,
     $           IWORK, LENIW )
C
            CALL DUTERR( 'DSLUCN',IERR,IOUT,NFAIL,ISTDO,ITER,ERR )
C         
C         * * * * * *    DSDBCG   * * * * * *
C         
            IF( IOUT.GE.3 ) THEN
               WRITE(ISTDO,1000) 'DSDBCG',ITOL,ISYM
            ENDIF
            CALL DFILL( N, XITER, 0.0D0 )
C
            CALL DSDBCG(N, F, XITER, NELT, IA, JA, A, ISYM, ITOL,
     $           TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW,
     $           IWORK, LENIW )
C
            CALL DUTERR( 'DSDBCG',IERR,IOUT,NFAIL,ISTDO,ITER,ERR )
C         
C         * * * * * *   DSLUBC   * * * * * *
C         
            IF( IOUT.GE.3 ) THEN
               WRITE(ISTDO,1000) 'DSLUBC',ITOL,ISYM
            ENDIF
            CALL DFILL( N, XITER, 0.0D0 )
C
            CALL DSLUBC(N, F, XITER, NELT, IA, JA, A, ISYM, 
     $           ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, 
     $           RWORK, LENW, IWORK, LENIW )
C
            CALL DUTERR( 'DSLUBC',IERR,IOUT,NFAIL,ISTDO,ITER,ERR )
C         
C         * * * * * *    DSDCGS   * * * * * *
C         
            IF( IOUT.GE.3 ) THEN
               WRITE(ISTDO,1000) 'DSDCGS',ITOL,ISYM
            ENDIF
            CALL DFILL( N, XITER, 0.0D0 )
C
            CALL DSDCGS(N, F, XITER, NELT, IA, JA, A, ISYM, ITOL,
     $           TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW,
     $           IWORK, LENIW )
C
            CALL DUTERR( 'DSDCGS',IERR,IOUT,NFAIL,ISTDO,ITER,ERR )
C         
C         * * * * * *   DSLUCS   * * * * * *
C         
            IF( IOUT.GE.3 ) THEN
               WRITE(ISTDO,1000) 'DSLUCS',ITOL,ISYM
            ENDIF
            CALL DFILL( N, XITER, 0.0D0 )
C
            CALL DSLUCS(N, F, XITER, NELT, IA, JA, A, ISYM, 
     $           ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, 
     $           RWORK, LENW, IWORK, LENIW )
C
            CALL DUTERR( 'DSLUCS',IERR,IOUT,NFAIL,ISTDO,ITER,ERR )
C         
C         * * * * * *    DSDOMN   * * * * * *
C         
CVD$ NOVECTOR
            DO 30 NSAVE = 0, 3
               IF( IOUT.GE.3 ) THEN
                  WRITE(ISTDO,1010) 'DSDOMN',ITOL, ISYM, NSAVE
               ENDIF
               CALL DFILL( N, XITER, 0.0D0 )
C
               CALL DSDOMN(N, F, XITER, NELT, IA, JA, A,
     $              ISYM, NSAVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, 
     $              IUNIT, RWORK, LENW, IWORK, LENIW )
C
               CALL DUTERR( 'DSDOMN',IERR,IOUT,NFAIL,ISTDO,ITER,ERR )
 30         CONTINUE
C         
C         * * * * * *   DSLUOM   * * * * * *
C         
CVD$ NOVECTOR
            DO 40 NSAVE=0,3
               IF( IOUT.GE.3 ) THEN
                  WRITE(ISTDO,1010) 'DSLUOM',ITOL, ISYM, NSAVE
               ENDIF
               CALL DFILL( N, XITER, 0.0D0 )
C
               CALL DSLUOM(N, F, XITER, NELT, IA, JA, A,
     $              ISYM, NSAVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, 
     $              IUNIT, RWORK, LENW, IWORK, LENIW )
C
               CALL DUTERR( 'DSLUOM',IERR,IOUT,NFAIL,ISTDO,ITER,ERR )
 40         CONTINUE
C         
C         * * * * * *   DSDGMR   * * * * * *
C         
CVD$ NOVECTOR
            DO 50 NSAVE = 5, 12
               IF( IOUT.GE.3 ) THEN
                  WRITE(ISTDO,1010) 'DSDGMR',ITOL, ISYM, NSAVE
               ENDIF
               CALL DFILL( N, XITER, 0.0D0 )
               ITOLGM = 0
C
               CALL DSDGMR(N, F, XITER, NELT, IA, JA, A,
     $              ISYM, NSAVE, ITOLGM, TOL, ITMAX, ITER, ERR, IERR, 
     $              IUNIT, RWORK, LENW, IWORK, LENIW )
C
               CALL DUTERR( 'DSDGMR',IERR,IOUT,NFAIL,ISTDO,ITER,ERR )
 50         CONTINUE
C         
C         * * * * * *   DSLUGM   * * * * * *
C         
CVD$ NOVECTOR
            DO 60 NSAVE = 5, 12
               IF( IOUT.GE.3 ) THEN
                  WRITE(ISTDO,1010) 'DSLUGM',ITOL, ISYM, NSAVE
               ENDIF
               CALL DFILL( N, XITER, 0.0D0 )
C
               CALL DSLUGM(N, F, XITER, NELT, IA, JA, A,
     $              ISYM, NSAVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, 
     $              IUNIT, RWORK, LENW, IWORK, LENIW )
C
               CALL DUTERR( 'DSLUGM',IERR,IOUT,NFAIL,ISTDO,ITER,ERR )
 60         CONTINUE
 20      CONTINUE
 10   CONTINUE
C         
      IF( NFAIL.EQ.0 ) THEN
         WRITE(ISTDO,1050)
      ELSE
         WRITE(ISTDO,1060) NFAIL
      ENDIF
C
      STOP 'All Done'
C
 1000 FORMAT(/1X,A6,' : ITOL = ',I2,'   ISYM = ',I1)
 1010 FORMAT(/1X,A6,' : ITOL = ',I2,'   ISYM = ',I1,' NSAVE = ',I2) 
 1020 FORMAT(/'                * RANDOM Matrix of size',I5,'*'
     $     /'                ',
     $     'Number of non-zeros & Density = ', I5,E16.7)
 1030 FORMAT('                Error tolerance = ',E16.7) 
 1040 FORMAT(/'  ***** SLAP Column Matrix *****'/
     $        ' Indx   ia   ja     a'/(1X,I4,1X,I4,1X,I4,1X,E16.7))
 1050 FORMAT(//
     $     '*******************************************************'/
     $     '**** All SLAP Double Precision Quick Checks Passed ****'/
     $     '****                 No Errors                     ****'/
     $     '*******************************************************')
 1060 FORMAT(//
     $     '************************************************'/
     $     '**     ===>',I3,' Failures detected <===      **'/
     $     '**     SLAP Double Precision Quick Checks     **'/
     $     '** Set KPRINT = 3 for DEBUG information and   **'/
     $     '** rerun the tests to determine the problem.  **'/
     $     '************************************************')
      END
*DECK DUTERR
      SUBROUTINE DUTERR( METHOD, IERR, IOUT, NFAIL, ISTDO, ITER, ERR )
C***BEGIN PROLOGUE  DUTERR
C***SUBSIDIARY
C***DATE WRITTEN   890404   (YYMMDD)
C***REVISION DATE  890404   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=DOUBLE PRECISION(DUTERR-D),
C             Linear system, Sparse, Iterative Precondition
C***AUTHOR  Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Output error messages for the SLAP Quick Check
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DUTERR
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*6 METHOD
      INTEGER IERR, IOUT, NFAIL, ISTDO, ITER
      DOUBLE PRECISION ERR
C
C***FIRST EXECUTABLE STATEMENT  DUTERR
      IF( IERR.NE.0 ) NFAIL = NFAIL+1    
      IF( IOUT.EQ.1 .AND. IERR.NE.0 ) THEN
         WRITE(ISTDO,1000) METHOD
      ENDIF
      IF( IOUT.EQ.2 ) THEN
         IF( IERR.EQ.0 ) THEN
            WRITE(ISTDO,1010) METHOD
         ELSE
            WRITE(ISTDO,1020) METHOD,IERR,ITER,ERR
         ENDIF
      ENDIF
      IF( IOUT.GE.3 ) THEN
         IF( IERR.EQ.0 ) THEN
            WRITE(ISTDO,1030) METHOD,IERR,ITER,ERR
         ELSE
            WRITE(ISTDO,1020) METHOD,IERR,ITER,ERR
         ENDIF
      ENDIF
      RETURN
 1000 FORMAT( 1X,A6,' : **** FAILURE ****')
 1010 FORMAT( 1X,A6,' : **** PASSED  ****')
 1020 FORMAT(' **************** WARNING ***********************'/
     $       ' **** ',A6,' Quick Test FAILED: IERR = ',I5,' ****'/
     $       ' **************** WARNING ***********************'/
     $       ' Iteration Count = ',I3,' Stop Test = ',E12.6)
 1030 FORMAT(' ***************** PASSED ***********************'/
     $       ' **** ',A6,' Quick Test PASSED: IERR = ',I5,' ****'/
     $       ' ***************** PASSED ***********************'/
     $       ' Iteration Count = ',I3,' Stop Test = ',E12.6)
C------------- LAST LINE OF DUTERR FOLLOWS ----------------------------
      END
*DECK DRMGEN
      SUBROUTINE DRMGEN( NELTMX, FACTOR, IERR, N, NELT, ISYM, 
     $     IA, JA, A, F, SOLN, DSUM, ITMP, IDIAG )
C***BEGIN PROLOGUE  DRMGEN
C***SUBSIDIARY
C***PURPOSE  This routine generates a "Random" symmetric or 
C            non-symmetric matrix of size N for use in the SLAP
C            Quick Checks.
C***LIBRARY   SLATEC(SLAP)
C***AUTHOR  Seager, Mark K., (LLNL)
C             seager@lll-crg.llnl.gov
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550
C             (415)423-3141
C***DESCRIPTION
C
C *Usage:
C       INTEGER NELTMX, IERR, N, NELT, ISYM, 
C       INTEGER IA(NELTMX), JA(NELTMX), ITMP(N), IDIAG(N)
C       DOUBLE PRECISION FACTOR, A(NELTMX), F(N), SOLN(N), DSUM(N)
C
C       CALL DRMGEN( NELTMX, FACTOR, IERR, N, NELT, ISYM, 
C      $     IA, JA, A, F, SOLN, DSUM, ITMP, IDIAG )
C
C *Arguments:
C  
C NELTMX :IN       Integer.
C         Maximum number of non-zeros that can be created by this
C         routine for storage in the IA, JA, A arrays,  see below.
C FACTOR :IN       Double Precision.
C         Non-zeros in the upper triangle are set to FACTOR times
C         the coresponding entry in the lower triangle when a non-
C         symmetric matrix is requested (See ISYM, below).
C IERR   :OUT      Integer.
C         Return error flag.  
C             IERR = 0 => everything went OK. 
C                  = 1 => Ran out of space trying to create matrix.
C                         Set NELTMX to something larger and retry.
C N      :IN       Integer.
C         Size of the linear system to generate (number of unknowns).
C NELT   :OUT      Integer.
C         Number of non-zeros stored in the IA, JA, A arrays, see below.
C ISYM   :IN       Integer.
C         Flag to indicate the type of matrix to generate:
C             ISYM = 0 => Non-Symmetric Matrix (See FACTOR, above).
C                  = 1 => Symmetric Matrix.
C IA     :OUT      Integer IA(NELTMX).
C         Stores the row indicies for the non-zeros.
C JA     :OUT      Integer JA(NELTMX).
C         Stores the column indicies for the non-zeros.
C A      :OUT      Double Precision A(NELTMX).
C         Stores the values of the non-zeros.
C F      :OUT      Double Precision F(N).
C         The right hand side of the linear system.  Obtained by mult-
C         iplying the matrix time SOLN, see below.
C SOLN   :OUT      Double Precision SOLN(N).
C         The true solution to the linear system.  Each component is
C         chosen at random (0.0<SOLN(I)<1.0, I=1,N)
C DSUM   :WORK     Double Precision DSUM(N).
C ITMP   :WORK     Integer ITMP(N).
C IDIAG  :WORK     Integer IDIAG(N).
C
C *Description
C         The matrix is generated by choosing a random number of 
C         entries for each column and then chosing negative random 
C         numbers for each off diagionals.   The diagionals elements 
C         are chosen to be positive and large enough so the matrix 
C         is slightly diagionally domainate.  The lower triangle of 
C         the matrix is generated and if isym.eq.0 (all matrix elements 
C         stored) the upper triangle elements are chosen so that they
C         are FACTOR times the coresponding lower triangular element.
C
C***ROUTINES CALLED  RAND, DMPL
C***REVISION HISTORY  (YYMMDD)
C   881120  DATE WRITTEN
C***END PROLOGUE  DRMGEN
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER NELTMX, IERR, N, NELT, ISYM
      INTEGER IA(NELTMX), JA(NELTMX)
      INTEGER ITMP(N), IDIAG(N)
      DOUBLE PRECISION FACTOR, A(NELTMX)
      DOUBLE PRECISION F(N), SOLN(N), DSUM(N)
      INTEGER DUMMY
      REAL, EXTERNAL :: RAND
C
C         Start by setting the random number generator seed.
C         This is done for reproducablility in debuggin.  Remove
C         the seed seeting call for production testing.
C
C***FIRST EXECUTABLE STATEMENT  DRMGEN
      DUMMY = 16381
      ISEED = RAND( DUMMY )
      IERR = 0
      DO 10 I = 1, N
         IDIAG(I) = 0
         DSUM(I) = -1.0
 10   CONTINUE
C
C         Set the matrix elements.
C         Loop over the columns.
      DUMMY = 0.0
      NELT = 0
CVD$ NOCONCUR
      DO 30 ICOL = 1, N
         NL = N+1-ICOL
C
C         To keep things sparse divide by two, three or four or ...
C
         INUM = (IFIX( RAND(DUMMY)*NL ) + 1)/3
         CALL DMPL( NL, INUM, ITMP )
C
C         Set up this column (and row, if non-sym structure).
CVD$ NOVECTOR
CVD$ NOCONCUR
         DO 20 IROW = 1, INUM
            NELT = NELT + 1
            IF( NELT.GT.NELTMX ) THEN
               IERR = 1
               RETURN
            ENDIF
            IA(NELT) = N+1-ITMP(IROW)
            JA(NELT) = ICOL
            IF( IA(NELT).EQ.ICOL ) THEN
               IDIAG(ICOL) = NELT
            ELSE
               A(NELT) = -RAND(DUMMY)
               DSUM(ICOL) = DSUM(ICOL) + A(NELT)
               IF( ISYM.EQ.0 ) THEN
C
C         Copy this element into upper triangle.
C
                  NELT = NELT + 1
                  IF( NELT.GT.NELTMX ) THEN
                     IERR = 1
                     RETURN
                  ENDIF
                  IA(NELT) = ICOL
                  JA(NELT) = IA(NELT-1)
                  A(NELT)  = A(NELT-1)*FACTOR
                  DSUM(JA(NELT)) = DSUM(JA(NELT)) + A(NELT)
               ELSE
                  DSUM(IA(NELT)) = DSUM(IA(NELT)) + A(NELT)
               ENDIF
            ENDIF
 20      CONTINUE
         IF( IDIAG(ICOL).EQ.0 ) THEN
C
C         Add a diagional to the column.
C
            NELT = NELT + 1
            IF( NELT.GT.NELTMX ) THEN
               IERR = 1
               RETURN
            ENDIF
            IDIAG(ICOL) = NELT
            A(NELT) = 0.0D0
            IA(NELT) = ICOL
            JA(NELT) = ICOL
         ENDIF
 30   CONTINUE
C
C         Clean up the diagionals.
C
CVD$ NODEPCHK
CLLL. OPTION ASSERT (NOHAZARD)
CDIR$ IVDEP
      DO 40 I = 1, N
         A(IDIAG(I)) = -1.0001*DSUM(I)
 40   CONTINUE
C
C         Set a random soln and determine the right-hand side.
C
CVD$ NOVECTOR
CVD$ NOCONCUR
      DO 50 I = 1, N
         SOLN(I) = RAND(DUMMY)
         F(I) = 0.0D0
 50   CONTINUE
C
CVD$ NOVECTOR
CVD$ NOCONCUR
      DO 60 K = 1, NELT
         F(IA(K)) = F(IA(K)) + A(K)*SOLN(JA(K))
         IF( ISYM.NE.0 .AND. IA(K).NE.JA(K) ) THEN
            F(JA(K)) = F(JA(K)) + A(K)*SOLN(IA(K))
         ENDIF
 60   CONTINUE
      RETURN
C------------- LAST LINE OF DRMGEN FOLLOWS ----------------------------
      END
*DECK DMPL
      SUBROUTINE DMPL( N, M, INDX )
C***BEGIN PROLOGUE  DMPL
C***SUBSIDIARY
C***DATE WRITTEN   890404   (YYMMDD)
C***REVISION DATE  890404   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=DOUBLE PRECISION(DMPL-D),
C             Linear system, Sparse, Iterative Precondition
C***AUTHOR  Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Picks m special random integers in the range 1 to n.
C            This routine picks m "random" integers in the range 1 to
C            n with out any repetitions.
C***ROUTINES CALLED  RAND
C***END PROLOGUE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL, EXTERNAL :: RAND
      INTEGER DUMMY
      INTEGER N, M, INDX(M)
C
C..       Check the input
C***FIRST EXECUTABLE STATEMENT  DMPL
      DUMMY = 0
      IF( N*M.LT.0 .OR. M.GT.N ) RETURN
C
C..       Set the indeicies.      
      INDX(1) = IFIX( RAND(DUMMY)*N ) + 1
CVD$ NOCONCUR
      DO 30 I = 2, M
 10      ID = IFIX( RAND(DUMMY)*N ) + 1
C
C..       Check to see if id has already been chosen.
CVD$ NOVECTOR
CVD$ NOCONCUR
         DO 20 J = 1, I-1
            IF( ID.EQ.INDX(J) ) GOTO 10
 20      CONTINUE
         INDX(I) = ID
 30   CONTINUE
      RETURN
C------------- LAST LINE OF DMPL FOLLOWS ------------------------------
      END
*DECK DFILL
      SUBROUTINE DFILL (N,V,VAL)
C***BEGIN PROLOGUE  DFILL
C***DATE WRITTEN   890404   (YYMMDD)
C***REVISION DATE  890404   (YYMMDD)
C***CATEGORY NO.  D2A4,
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=DOUBLE PRECISION(DFILL-D),
C             Linear system, Sparse, Iterative Precondition
C***AUTHOR  Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Fill a vector with a value.
C***DESCRIPTION
C *Usage:
C     INTEGER  N
C     DOUBLE PRECISION V(N), VAL
C
C     CALL DFILL( N, V, VAL )
C         
C *Arguments:
C N      :IN       Integer.
C         Length of the vector 
C V      :OUT      Double Precision V(N).
C         Vectored to be set.
C VAL    :IN       Double Precision.
C         Value to seed the vector with.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DFILL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N
      DOUBLE PRECISION V(*), VAL
C
C***FIRST EXECUTABLE STATEMENT  DFILL
      IF (N .LE. 0) RETURN
      NR=MOD(N,4)
C
C         The following construct assumes a zero pass do loop.
C
      IS=1
      GOTO(1,2,3,4), NR+1
    4   IS=4
        V(1)=VAL
        V(2)=VAL
        V(3)=VAL
        GOTO 1
    3   IS=3
        V(1)=VAL
        V(2)=VAL
        GOTO 1
    2   IS=2
        V(1)=VAL
    1 DO 10 I=IS,N,4
        V(I)  =VAL
        V(I+1)=VAL
        V(I+2)=VAL
        V(I+3)=VAL
 10   CONTINUE
      RETURN
C------------- LAST LINE OF DFILL FOLLOWS -----------------------------
      END
