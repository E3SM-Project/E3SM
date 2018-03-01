!/ ------------------------------------------------------------------- /
      MODULE W3INITMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    28-Dec-2004 : Origination (out of W3WAVEMD).      ( version 3.06 )
!/                  Multiple grid version.
!/    03-Jan-2005 : Add US2x to MPI communication.      ( version 3.06 )
!/    04-Jan-2005 : Add grid output flags to W3INIT.    ( version 3.06 )
!/    07-Feb-2005 : Combined vs. separate test output.  ( version 3.07 )
!/    04-May-2005 : Change to MPI_COMM_WAVE.            ( version 3.07 )
!/    21-Jul-2005 : Add output fields.                  ( version 3.07 )
!/    09-Nov-2005 : Drying out of points added.         ( version 3.08 )
!/    13-Jun-2006 : Splitting STORE in G/SSTORE.        ( version 3.09 )
!/    26-Jun-2006 : adding wiring for output type 6.    ( version 3.09 )
!/    27-Jun-2006 : Adding file name preamble.          ( version 3.09 )
!/    04-Jul-2006 : Consolidate stress arrays.          ( version 3.09 )
!/    02-Aug-2006 : Adding W3MPIP.                      ( version 3.10 )
!/    02-Nov-2006 : Adding partitioning options.        ( version 3.10 )
!/    11-Jan-2007 : Updating IAPPRO computation.        ( version 3.10 )
!/    02-Apr-2007 : Add partitioned field data.         ( version 3.11 )
!/                  Add user-defined field data.
!/    01-May-2007 : Move O7a output to W3IOPP.          ( version 3.11 )
!/    08-May-2007 : Starting from calm as an option.    ( version 3.11 )
!/    17-May-2007 : Adding NTPROC/NAPROC separation.    ( version 3.11 )
!/    21-Jun-2007 : Dedicated output processes.         ( version 3.11 )
!/    29-Feb-2008 : Add NEC compiler directives.        ( version 3.13 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!/    Note: Changes in version numbers not logged above.
!/
!  1. Purpose :
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      CRITOS    R.P.  Public   Critical percentage of resources used
!                               for output to trigger warning.
!      WWVER     C*10  Public   Model version number.
!     ----------------------------------------------------------------
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3INIT    Subr. Public   Wave model initialization.
!      W3MPII    Subr. Public   Initialize MPI data transpose.
!      W3MPIO    Subr. Public   Initialize MPI output gathering.
!      W3MPIP    Subr. Public   Initialize MPI point output gathering.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!     See subroutine documentation.
!
!  5. Remarks :
!
!  6. Switches :
!
!       !/SHRD  Switch for shared / distributed memory architecture.
!       !/DIST  Id.
!       !/MPI   Id.
!       !/LLG   Spherical grid.
!       !/XYG   Cartesian grid.
!
!       !/S     Enable subroutine tracing.
!       !/Tn    Enable test output.
!       !/MPIT  Enable test output (MPI).
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
      PUBLIC
!/
      REAL, PARAMETER              :: CRITOS = 15.
      CHARACTER*10, PARAMETER      :: WWVER  = '3.14 '
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3INIT ( IMOD, FEXT, MDS, MTRACE, ODAT,              &
                          FLGRD, NPT, XPT, YPT, PNAMES, IPRT, PRTFRM, &
                          MPI_COMM )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         21-Jun-2007 |
!/                  +-----------------------------------+
!/
!/    17-Mar-1999 : Distributed FORTRAN 77 version.     ( version 1.18 )
!/    13-Jan-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/                  Major changes to logistics.
!/    14-Feb-2000 : Exact-NL added.                     ( version 2.01 )
!/    24-Jan-2001 : Flat grid version.                  ( version 2.06 )
!/    24-Jan-2002 : Zero time step for data ass.        ( version 2.17 )
!/    18-Feb-2002 : Point output diagnostics added.     ( version 2.18 )
!/    13-Nov-2002 : Add stress vector.                  ( version 3.00 )
!/    20-Aug-2003 : Output server options added.        ( version 3.04 )
!/    28-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/                  Taken out of W3WAVE.
!/    04-Jan-2005 : Add grid output flags to par list.  ( version 3.06 )
!/    07-Feb-2005 : Combined vs. separate test output.  ( version 3.07 )
!/    04-May-2005 : Change to MPI_COMM_WAVE.            ( version 3.07 )
!/    09-Nov-2005 : Drying out of points added.         ( version 3.08 )
!/    26-Jun-2006 : adding wiring for output type 6.    ( version 3.09 )
!/    27-Jun-2006 : Adding file name preamble.          ( version 3.09 )
!/    02-Aug-2006 : Adding W3MPIP.                      ( version 3.10 )
!/    02-Nov-2006 : Adding partitioning options.        ( version 3.10 )
!/    11-Jan-2007 : Updating IAPPRO computation.        ( version 3.10 )
!/    01-May-2007 : Move O7a output to W3IOPP.          ( version 3.11 )
!/    08-May-2007 : Starting from calm as an option.    ( version 3.11 )
!/    17-May-2007 : Adding NTPROC/NAPROC separation.    ( version 3.11 )
!/    21-Jun-2007 : Dedicated output processes.         ( version 3.11 )
!/
!  1. Purpose :
!
!     Initialize WAVEWATCH III.
!
!  2. Method :
!
!     Initialize data structure and wave fields from data files.
!     Initialize grid from local and instantaneous data.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.   I   Model number.
!       FEXT    Char   I   Extension of data files.
!       MDS     I.A.   I   Array with dataset numbers (see below),
!                          saved as NDS in W3ODATMD.
!                           1: General output unit number ("log file").
!                           2: Error output unit number.
!                           3: Test output unit number.
!                           4: "screen", i.e., direct output location,
!                              can be the screen or the output file of
!                              the shell.
!                           5: Model definition file unit number.
!                           6: Restart file unit number.
!                           7: Grid output file unit number.
!                           8: Point output file unit number.
!                           9: Input boundary data file unit number.
!                          10: Output boundary data file unit number
!                              (first).
!                          11: Track information file unit number.
!                          12: Track output file unit number.
!       MTRACE  I.A.   I   Array with subroutine tracing information.
!                           1: Output unit number for trace.
!                           2: Maximum number of trace prints.
!       ODAT    I.A.   I   Output data, five parameters per output type
!                           1-5  Data for OTYPE = 1; gridded fields.
!                                1 YYYMMDD for first output.
!                                2 HHMMSS for first output.
!                                3 Output interval in seconds.
!                                4 YYYMMDD for last output.
!                                5 HHMMSS for last output.
!                           6-10 Id. for OTYPE = 2; point output.
!                          11-15 Id. for OTYPE = 3; track point output.
!                          16-20 Id. for OTYPE = 4; restart files.
!                          21-25 Id. for OTYPE = 5; boundary data.
!       FLGRD   L.A.   I   Flags for gridded output.
!       NPT     Int.   I   Number of output points
!       X/YPT   R.A.   I   Coordinates of output points.
!       PNAMES  C.A.   I   Output point names.
!       IPRT    I.A.   I   Partitioning grid info.
!       PRTFRM  I.A.   I   Partitioning format flag.
!       MPI_COMM Int.  I   MPI communicator to be used for model.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SETG    Subr. W3GDATMD Point to data structure.
!      W3SETW    Subr. W3WDATMD Point to data structure.
!      W3DIMW    Subr.   Id.    Set array sizes in data structure.
!      W3SETA    Subr. W3ADATMD Point to data structure.
!      W3DIMA    Subr.   Id.    Set array sizes in data structure.
!      W3SETI    Subr. W3IDATMD Point to data structure.
!      W3DIMI    Subr.   Id.    Set array sizes in data structure.
!      W3SETO    Subr. W3ODATMD Point to data structure.
!      W3DMO5    Subr.   Id.    Set array sizes in data structure.
!      ITRACE    Subr. W3SERVMD Subroutine tracing initialization.
!      STRACE    Subr.   Id.    Subroutine tracing.
!      EXTCDE    Subr.   Id.    Program abort.
!      WWDATE    Subr.   Id.    System date.
!      WWTIME    Subr.   Id.    System time.
!      DSEC21    Func. W3TIMEMD Compute time difference.
!      TICK21    Func.   Id.    Advance the clock.
!      STME21    Func.   Id.    Print the time readable.
!      PRTBLK    Func. W3ARRYMD Print plot of array.
!      W3IOGR    Subr. W3IOGRMD Read/write model definition file.
!      W3IORS    Subr. W3IORSMD Read/write restart file.
!      W3IOPP    Subr. W3IOPOMD Preprocess point output.
!      CALL MPI_COMM_SIZE, CALL MPI_COMM_RANK
!                Subr. mpif.h   Standard MPI routines.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!     Any program shell or integrated model which uses WAVEWATCH III.
!
!  6. Error messages :
!
!     On opening of log file only. Other error messages are generated
!     by W3IOGR and W3IORS.
!
!  7. Remarks :
!
!     - The log file is called 'log.FEXT', where FEXT is passed to
!       the routine.
!     - The test output file is called 'test.FEXT' in shared memory
!       version or testNNN.FEXT in distributed memory version.
!     - A water level and ice coverage are transferred with the
!       restart file. To assure consistency within the model, the
!       water level and ice coverage are re-evaluated at the 0th
!       time step in the actual wave model routine.
!
!  8. Structure :
!
!     ----------------------------------------------------
!      1.  Set-up of idata structures and I/O.
!        a Point to proper data structures.
!        b Number of processors and processor number.
!        c Open files.
!        d Dataset unit numbers
!        e Subroutine tracing
!        f Initial and test outputs
!      2.  Model definition.
!        a Read model definition file         ( W3IOGR )
!        b Save MAPSTA.
!        c MPP preparation
!      3.  Model initialization.
!        a Read restart file.                 ( W3IORS )
!        b Compare grid and restart MAPSTA.
!        c Initialize with winds if requested (set flag).
!        d Initialize calm conditions if requested.
!        e Preparations for prop. scheme.
!      4.  Set-up output times.
!        a Unpack ODAT.
!        b Check if output available.
!        c Get first time per output and overall.
!        d Prepare point output               ( W3IOPP )
!      5.  Define wavenumber grid.
!        a Calculate depth.
!        b Fill wavenumber and group velocity arrays.
!      6.  Initialize arrays.
!      7.  Write info to log file.
!      8.  Final MPI set up  ( W3MPII , W3MPIO , W3MPIP )
!     ----------------------------------------------------
!
!  9. Switches :
!
!       !/SHRD  Switch for shared / distributed memory architecture.
!       !/DIST  Id.
!       !/MPI   Id.
!       !/LLG   Spherical grid.
!       !/XYG   Cartesian grid.
!
!       !/S     Enable subroutine tracing.
!       !/Tn    Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3CONSTANTS
!/
      USE W3GDATMD, ONLY: W3SETG
      USE W3WDATMD, ONLY: W3SETW, W3DIMW
      USE W3ADATMD, ONLY: W3SETA, W3DIMA
      USE W3IDATMD, ONLY: W3SETI, W3DIMI
      USE W3ODATMD, ONLY: W3SETO, W3DMO5
      USE W3IOGRMD, ONLY: W3IOGR
      USE W3IORSMD, ONLY: W3IORS
      USE W3IOPOMD, ONLY: W3IOPP
      USE W3SERVMD, ONLY: ITRACE, EXTCDE, WWDATE, WWTIME
      USE W3TIMEMD, ONLY: DSEC21, TICK21, STME21
      USE W3ARRYMD, ONLY: PRTBLK
!/
      USE W3GDATMD, ONLY: NX, NY, NSEA, NSEAL, MAPSTA, MAPST2, MAPFS, &
                          GLOBAL, ZB, TRNX, TRNY, DMIN, DTCFL, DTMAX, &
                          FLCK, NK, NTH, NSPEC, SIG, GNAME,           &
                          X0, Y0, SX, SY
      USE W3WDATMD, ONLY: TIME, TLEV, TICE, WLV, UST, USTDIR, VA
      USE W3ODATMD, ONLY: NDSO, NDSE, NDST, SCREEN, NDS, NTPROC,      &
                          NAPROC, IAPROC, NAPLOG, NAPOUT, NAPERR,     &
                          NAPFLD, NAPPNT, NAPTRK, NAPRST, NAPBPT,     &
                          NAPPRT, TOFRST, DTOUT, TONEXT, TOLAST,      &
                          FLOUT, NOGRD, FLOGRD, FLBPO, NOPTS, PTNME,  &
                          PTLOC, IPTINT, PTIFAC, UNDEF, IDOUT,        &
                          OUTPTS, FNMPRE, IX0, IXN, IXS, IY0, IYN,    &
                          IYS, FLFORM, IOSTYP, UNIPTS, UPPROC
      USE W3ADATMD, ONLY: NSEALM, IAPPRO, FLCOLD, FLIWND, DW, CG, WN, &
                          UA, UD, U10, U10D, AS
      USE W3ADATMD, ONLY: MPI_COMM_WAVE, MPI_COMM_WCMP
      USE W3IDATMD, ONLY: FLLEV, FLCUR, FLWIND, FLICE
      USE W3DISPMD, ONLY: NAR1D, DFAC, N1MAX, ECG1, EWN1, DSIE
!/
      IMPLICIT NONE
!
      INCLUDE "mpif.h"
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)           :: IMOD, MDS(13), MTRACE(2),      &
                                       ODAT(30), NPT, IPRT(6), MPI_COMM
      REAL, INTENT(INOUT)           :: XPT(NPT), YPT(NPT)
      LOGICAL, INTENT(IN)           :: FLGRD(NOGRD), PRTFRM
      CHARACTER*(*), INTENT(IN)     :: FEXT
      CHARACTER(LEN=10), INTENT(IN) :: PNAMES(NPT)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IE, IFL, IFT, IERR, NTTOT, NTLOC,    &
                                 NTTARG, IK, IP, ITH, INTYPE, IX, IY, &
                                 J, J0, TOUT(2), TLST(2), ISEA, IS,   &
                                 I1, I2, JSEA, NTTMAX
      INTEGER                 :: ISTEP, ISP
      INTEGER                 :: IERR_MPI, BGROUP, LGROUP
      INTEGER, ALLOCATABLE    :: TMPRNK(:)
      INTEGER, ALLOCATABLE    :: NT(:), MAPTST(:,:)
      REAL                    :: DTTST, DEPTH, SQRTH, SIX, R1, R2,    &
                                 FRACOS
       REAL                   :: FACTOR = 1.
      LOGICAL                 :: OPENED
      CHARACTER(LEN=8)        :: STTIME
      CHARACTER(LEN=10)       :: STDATE
      CHARACTER(LEN=23)       :: DTME21
      CHARACTER(LEN=30)       :: LFILE, TFILE
!/
!/ ------------------------------------------------------------------- /
!
! 1.  Set-up of data structures and I/O  ----------------------------- /
! 1.a Point to proper data structures.
!
      CALL W3SETO ( IMOD, MDS(2), MDS(3) )
      CALL W3SETG ( IMOD, MDS(2), MDS(3) )
      CALL W3SETW ( IMOD, MDS(2), MDS(3) )
      CALL W3SETA ( IMOD, MDS(2), MDS(3) )
      CALL W3SETI ( IMOD, MDS(2), MDS(3) )
!
! 1.b Number of processors and processor number.
!     Overwrite some initializations from W3ODATMD.
!
!     *******************************************************
!     *** NOTE : OUTPUT PROCESSOR ASSIGNMENT NEEDS TO BE  ***
!     ***        CONSISTENT WITH ASSIGNMENT IN WMINIT.    ***
!     *******************************************************
!
      MPI_COMM_WAVE = MPI_COMM
      CALL MPI_COMM_SIZE ( MPI_COMM_WAVE, NTPROC, IERR_MPI )
      NAPROC = NTPROC
      CALL MPI_COMM_RANK ( MPI_COMM_WAVE, IAPROC, IERR_MPI )
      IAPROC = IAPROC + 1
!
      IF ( IOSTYP .LE. 1 ) THEN
!
          NAPFLD = MAX(1,NAPROC-1)
          NAPPNT = MAX(1,NAPROC-2)
          NAPTRK = MAX(1,NAPROC-5)
          NAPRST = NAPROC
          NAPBPT = MAX(1,NAPROC-3)
          NAPPRT = MAX(1,NAPROC-4)
!
        ELSE
!
          NAPPNT = NAPROC
          IF ( UNIPTS .AND. UPPROC ) NAPROC = MAX(1,NTPROC - 1)
          NAPFLD = NAPROC
          NAPRST = NAPROC
          NAPBPT = NAPROC
          NAPTRK = NAPROC
          NAPPRT = NAPROC
!
          IF ( IOSTYP .EQ. 2 ) THEN
              NAPROC = MAX(1,NAPROC-1)
            ELSE IF ( IOSTYP .EQ. 3 ) THEN
              IF ( ODAT( 3).GT.0 ) THEN
                  NAPFLD =       NAPROC
                  NAPROC = MAX(1,NAPROC-1)
                END IF
              IF ( ODAT(13).GT.0 ) THEN
                  NAPTRK =       NAPROC
                  NAPROC = MAX(1,NAPROC-1)
                END IF
              IF ( ODAT(28).GT.0 ) THEN
                  NAPPRT =       NAPROC
                  NAPROC = MAX(1,NAPROC-1)
                END IF
              IF ( ODAT( 8).GT.0 ) NAPPNT = NAPROC
              IF ( ODAT(18).GT.0 ) NAPRST = NAPROC
              IF ( ODAT(23).GT.0 ) NAPBPT = NAPROC
              IF ( ( ODAT( 8).GT.0 .OR. ODAT(18).GT.0 .OR.            &
                     ODAT(23).GT.0 ) ) NAPROC = MAX(1,NAPROC-1)
            END IF
        END IF
!
      FRACOS = 100. * REAL(NTPROC-NAPROC) / REAL(NTPROC)
      IF ( FRACOS.GT.CRITOS .AND. IAPROC.EQ.NAPERR )                  &
                                           WRITE (NDSE,8002) FRACOS
!
      IF ( NAPROC .EQ. NTPROC ) THEN
          MPI_COMM_WCMP = MPI_COMM_WAVE
        ELSE
          CALL MPI_COMM_GROUP ( MPI_COMM_WAVE, BGROUP, IERR_MPI )
          ALLOCATE ( TMPRNK(NAPROC) )
          DO J=1, NAPROC
            TMPRNK(J) = J - 1
            END DO
          CALL MPI_GROUP_INCL ( BGROUP, NAPROC, TMPRNK, LGROUP,  &
                                IERR_MPI )
          CALL MPI_COMM_CREATE ( MPI_COMM_WAVE, LGROUP,          &
                                 MPI_COMM_WCMP, IERR_MPI )
          CALL MPI_GROUP_FREE ( LGROUP, IERR_MPI )
          CALL MPI_GROUP_FREE ( BGROUP, IERR_MPI )
          DEALLOCATE ( TMPRNK )
        END IF
!
! 1.c Open files without unpacking MDS ,,,
!
      IE     = LEN_TRIM(FEXT)
      LFILE  = 'log.' // FEXT(:IE)
      IFL    = LEN_TRIM(LFILE)
       WRITE (TFILE,'(A4,I3.3,2A)') 'test',                     &
                    OUTPTS(IMOD)%IAPROC, '.', FEXT(:IE)
      IFT    = LEN_TRIM(TFILE)
      J      = LEN_TRIM(FNMPRE)
!
!MV      IF ( OUTPTS(IMOD)%IAPROC .EQ. OUTPTS(IMOD)%NAPLOG )             &
!MV          OPEN (MDS(1),FILE=FNMPRE(:J)//LFILE(:IFL),ERR=888,IOSTAT=IERR)
!
      IF ( MDS(3).NE.MDS(1) .AND. MDS(3).NE.MDS(4) ) THEN
          INQUIRE (MDS(3),OPENED=OPENED)
          IF ( .NOT. OPENED ) OPEN                                    &
               (MDS(3),FILE=FNMPRE(:J)//TFILE(:IFT),ERR=889,IOSTAT=IERR)
        END IF
! 1.d Dataset unit numbers
!
      NDS    = MDS
      NDSO   = NDS(1)
      NDSE   = NDS(2)
      NDST   = NDS(3)
      SCREEN = NDS(4)
!
! 1.e Subroutine tracing
!
      CALL ITRACE ( MTRACE(1), MTRACE(2) )
!
! 1.f Initial and test outputs
!
      IF ( IAPROC .EQ. NAPLOG ) THEN
          CALL WWDATE ( STDATE )
          CALL WWTIME ( STTIME )
          WRITE (NDSO,900) WWVER, STDATE, STTIME
        END IF
!
! 2.  Model defintition ---------------------------------------------- /
! 2.a Read model defintition file
!
      CALL W3IOGR ( 'READ', NDS(5), IMOD, FEXT )
      IF ( IAPROC .EQ. NAPLOG ) WRITE (NDSO,920)
!
! 2.b Save MAPSTA
!
      ALLOCATE ( MAPTST(NY,NX) )
      MAPTST  = MAPSTA
!
! 2.c MPP preparation
! 2.c.1 Set simple counters and variables
!
      IF ( IAPROC .LE. NAPROC ) THEN
          NSEAL  = 1 + (NSEA-IAPROC)/NAPROC
        ELSE
          NSEAL  = 0
        END IF
      NSEALM = 1 + (NSEA-1)/NAPROC
      IF ( NSEA .LT. NAPROC ) GOTO 820
      IF ( NSPEC .LT. NAPROC ) GOTO 821
!
! 2.c.2 Allocate arrays
!
      IF ( IAPROC .LE. NAPROC ) THEN
          CALL W3DIMW ( IMOD, NDSE, NDST )
        ELSE
          CALL W3DIMW ( IMOD, NDSE, NDST, .FALSE. )
        END IF
      CALL W3DIMA ( IMOD, NDSE, NDST )
      CALL W3DIMI ( IMOD, NDSE, NDST )
!
! 2.c.3 Calculated expected number of prop. calls per processor
!
      NTTOT  = 0
      DO IK=1, NK
        NTLOC  = 1 + INT(DTMAX/(DTCFL*SIG(IK)/SIG(1))-0.001)
        NTTOT  = NTTOT + NTLOC*NTH
        END DO
      NTTARG = 1 + (NTTOT-1)/NAPROC
      NTTARG = NTTARG + INT(DTMAX/(DTCFL*SIG(NK)/SIG(1))-0.001)
      NTTMAX = NTTARG + 5
!
! 2.c.4 Initialize IAPPRO
!
      IAPPRO = 1
      ALLOCATE ( NT(NSPEC) )
      NT     = NTTOT
!
      DO
!
! 2.c.5 First sweep filling IAPPRO
!
        DO IP=1, NAPROC
          ISTEP  = IP
          ISP    = 0
          NT(IP) = 0
          DO J=1, 1+NSPEC/NAPROC
            ISP    = ISP + ISTEP
            IF ( MOD(J,2) .EQ. 1 ) THEN
                ISTEP  = 2*(NAPROC-IP) + 1
              ELSE
                ISTEP  = 2*IP - 1
              END IF
            IF ( ISP .LE. NSPEC ) THEN
                IK     = 1 + (ISP-1)/NTH
                NTLOC  = 1 + INT(DTMAX/(DTCFL*SIG(IK)/SIG(1))-0.001)
                IF ( NT(IP)+NTLOC .LE. NTTARG ) THEN
                    IAPPRO(ISP) = IP
                    NT(IP)      = NT(IP) + NTLOC
                  ELSE
                    IAPPRO(ISP) = -1
                 END IF
              END IF
            END DO
          END DO
!
! 2.c.6 Second sweep filling IAPPRO
!
        DO IP=1, NAPROC
          IF ( NT(IP) .LT. NTTARG ) THEN
              DO ISP=1, NSPEC
                IF ( IAPPRO(ISP) .EQ. -1 ) THEN
                    IK     = 1 + (ISP-1)/NTH
                    NTLOC  = 1 + INT(DTMAX/(DTCFL*SIG(IK)/SIG(1))-0.001)
                  IF ( NT(IP)+NTLOC .LE. NTTARG ) THEN
                        IAPPRO(ISP) = IP
                        NT(IP)      = NT(IP) + NTLOC
                     END IF
                  END IF
                END DO
            END IF
          END DO
!
! 2.c.7 Check if all served
!
        IF ( MINVAL(IAPPRO(1:NSPEC)) .GT. 0 ) THEN
            EXIT
          ELSE
            NTTARG = NTTARG + 1
            IF ( NTTARG .GE. NTTMAX ) EXIT
            IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,8028)
          END IF
!
        END DO
!
! 2.c.8 Test output
!
! 2.c.9 Test if any spectral points are left out
!
      DO ISP=1, NSPEC
        IF ( IAPPRO(ISP) .EQ. -1. ) GOTO 829
        END DO
!
      DEALLOCATE ( NT )
!
! 3.  Model initialization ------------------------------------------- /
! 3.a Read restart file
!
      CALL W3IORS ( 'READ', NDS(6), SIG(NK), INTYPE, IMOD )
      FLCOLD = INTYPE.LE.1  .OR. INTYPE.EQ.4
      IF ( IAPROC .EQ. NAPLOG ) THEN
          IF (INTYPE.EQ.0) THEN
              WRITE (NDSO,930) 'cold start (idealized).'
            ELSE IF ( INTYPE .EQ. 1 ) THEN
              WRITE (NDSO,930) 'cold start (wind).'
            ELSE IF ( INTYPE .EQ. 4 ) THEN
              WRITE (NDSO,930) 'cold start (calm).'
            ELSE
              WRITE (NDSO,930) 'full restart.'
            END IF
        END IF
!
! 3.b Compare MAPSTA from grid and restart
!
      DO IX=1, NX
        DO IY=1, NY
          IF ( ABS(MAPSTA(IY,IX)).EQ.2 .OR.                           &
               ABS(MAPTST(IY,IX)).EQ.2 ) THEN
              MAPSTA(IY,IX) = SIGN ( MAPTST(IY,IX) , MAPSTA(IY,IX) )
            END IF
          END DO
        END DO
!
! 3.c Initialization from wind fields
!
      FLIWND = INTYPE.EQ.1
!
! 3.d Initialization with calm conditions
!
      IF ( INTYPE .EQ. 4 ) THEN
          VA     = 0.
        END IF
!
! 3.e Prepare propagation scheme
!
      IF ( .NOT. FLCUR ) FLCK = .FALSE.
!
! 4.  Set-up output times -------------------------------------------- *
! 4.a Unpack ODAT
!
      DO J=1, 6
        J0 = (J-1)*5
        TONEXT(1,J) =        ODAT(J0+1)
        TONEXT(2,J) =        ODAT(J0+2)
        DTOUT (  J) = REAL ( ODAT(J0+3) )
        TOLAST(1,J) =        ODAT(J0+4)
        TOLAST(2,J) =        ODAT(J0+5)
        END DO
!
! 4.b Check if output available
!
      FLOUT(1) = .FALSE.
      FLOGRD   = FLGRD
      DO J=1, NOGRD
        FLOUT(1) = FLOUT(1) .OR. FLOGRD(J)
        END DO
!
      FLOUT(2) = NPT .GT. 0
!
      FLOUT(3) = .TRUE.
!
      FLOUT(4) = .TRUE.
!
      FLOUT(5) = FLBPO
      IF ( FLBPO ) THEN
          CALL W3DMO5 ( IMOD, NDSE, NDST, 4 )
        ELSE
          DTOUT(5) = 0.
        END IF
!
      IX0    = MAX (  1, IPRT(1) )
      IXN    = MIN ( NX, IPRT(2) )
      IXS    = MAX (  1, IPRT(3) )
      IY0    = MAX (  1, IPRT(4) )
      IYN    = MIN ( NY, IPRT(5) )
      IYS    = MAX (  1, IPRT(6) )
      FLFORM = PRTFRM
      FLOUT(6) = IX0.LE.IXN .AND. IY0.LE.IYN
!
! 4.c Get first time per output and overall.
!
      TOFRST(1) = -1
      TOFRST(2) =  0
!
      DO J=1, 6
!
! ... check time step
!
        DTOUT(J) = MAX ( 0. , DTOUT(J) )
        FLOUT(J) = FLOUT(J) .AND. ( DTOUT(J) .GT. 0.5 )
!
! ... get first time
!
        IF ( FLOUT(J) ) THEN
            TOUT = TONEXT(:,J)
            TLST = TOLAST(:,J)
!
            DO
              DTTST   = DSEC21 ( TIME , TOUT )
              IF ( ( J.NE.4 .AND. DTTST.LT.0. ) .OR.                  &
                   ( J.EQ.4 .AND. DTTST.LE.0. ) ) THEN
                  CALL TICK21 ( TOUT, DTOUT(J) )
                ELSE
                  EXIT
                END IF
              END DO
!
! ... reset first time
!
            TONEXT(:,J) = TOUT
!
! ... check last time
!
            DTTST  = DSEC21 ( TOUT , TLST )
            IF ( DTTST.LT.0.) FLOUT(J) = .FALSE.
!
! ... check overall first time
!
            IF ( FLOUT(J) ) THEN
                IF ( TOFRST(1).EQ.-1 ) THEN
                    TOFRST = TOUT
                  ELSE
                    DTTST  = DSEC21 ( TOUT , TOFRST )
                    IF ( DTTST.GT.0.) THEN
                        TOFRST = TOUT
                      END IF
                  END IF
              END IF
!
          END IF
!
        END DO
!
! 4.d Preprocessing for point output.
!
      IF ( FLOUT(2) ) CALL W3IOPP ( NPT, XPT, YPT, PNAMES, IMOD )
!
! 5.  Define wavenumber grid ----------------------------------------- *
! 5.a Calculate depth
!
      MAPTST = MOD(MAPST2/2,2)
      MAPST2 = MAPST2 - 2*MAPTST
!
      DO IY=1, NY
        DO IX=1, NX
          ISEA   = MAPFS(IY,IX)
          IF ( ISEA .NE. 0) THEN
              DW(ISEA) = MAX ( 0. , WLV(ISEA)-ZB(ISEA) )
              IF ( WLV(ISEA)-ZB(ISEA) .LE.0. ) THEN
                  MAPTST(IY,IX) = 1
                  MAPSTA(IY,IX) = -ABS(MAPSTA(IY,IX))
                  IF ( MOD(ISEA-IAPROC,NAPROC) .EQ. 0 ) THEN
                      JSEA   = 1 + (ISEA-1)/NAPROC
                      VA(:,JSEA) = 0.
                    END IF
                END IF
            END IF
           END DO
         END DO
!
      MAPST2 = MAPST2 + 2*MAPTST
!
      DEALLOCATE ( MAPTST )
!
! 5.b Fill wavenumber and group velocity arrays.
!
      DO IS=1, NSEA
        DEPTH  = MAX ( DMIN , DW(IS) )
        SQRTH  = SQRT(DEPTH)
!
        DO IK=0, NK+1
!
! Inlined version of WAVNU1.
!
          SIX    = SIG(IK) * SQRTH
          I1     = INT(SIX/DSIE)
          IF (I1.LE.N1MAX) THEN
              I2 = I1 + 1
              R1 = SIX/DSIE - REAL(I1)
              R2 = 1. - R1
              WN(IK,IS) = ( R2*EWN1(I1) + R1*EWN1(I2) ) / DEPTH
              CG(IK,IS) = ( R2*ECG1(I1) + R1*ECG1(I2) ) * SQRTH
            ELSE
              WN(IK,IS) = SIG(IK)*SIG(IK)/GRAV
              CG(IK,IS) = 0.5 * GRAV / SIG(IK)
            END IF
!
          END DO
        END DO
!
      DO IK=1, NK
        CG(IK,0) = CG(IK,1)
        WN(IK,0) = WN(IK,1)
        END DO
!
! 6.  Initialize arrays ---------------------------------------------- /
!     Some initialized in W3IORS
!
      UA     = 0.
      UD     = 0.
      U10    = 0.
      U10D   = 0.
!
      AS     = UNDEF
!
      UST   (0) = 0.
      USTDIR(0) = 0.
      AS    (0) = 0.
      DW    (0) = 0.
!
! 7.  Write info to log file ----------------------------------------- /
!
      IF ( IAPROC .EQ. NAPLOG ) THEN
!
          WRITE (NDSO,970)GNAME
          WRITE (NDSO,*) '   Xaxis X0,SX : ',X0,SX
          WRITE (NDSO,*) '   Yaxis Y0,SY : ',Y0,SY
          IF (   FLLEV    ) WRITE (NDSO,971) 'Prescribed'
          IF (.NOT. FLLEV ) WRITE (NDSO,971) 'No'
          IF (   FLCUR    ) WRITE (NDSO,972) 'Prescribed'
          IF (.NOT. FLCUR ) WRITE (NDSO,972) 'No'
          IF (   FLWIND   ) WRITE (NDSO,973) 'Prescribed'
          IF (.NOT. FLWIND) WRITE (NDSO,973) 'No'
          IF (   FLICE    ) WRITE (NDSO,974) 'Prescribed'
          IF (.NOT. FLICE ) WRITE (NDSO,974) 'No'
!
          IF ( FLOUT(1) ) THEN
              WRITE (NDSO,975)
              DO J=1,NOGRD
                IF ( FLOGRD(J) ) WRITE (NDSO,976) IDOUT(J)
                END DO
            END IF
!
          IF ( FLOUT(2) ) THEN
              WRITE (NDSO,977) NOPTS
              IF ( NOPTS .EQ. 0 ) THEN
                  WRITE (NDSO,978)
                ELSE
                  WRITE (NDSO,979)
                  DO IP=1, NOPTS
                    WRITE (NDSO,980) IP, FACTOR*PTLOC(1,IP),          &
                                         FACTOR*PTLOC(2,IP), PTNME(IP)
                    END DO
                END IF
            END IF
!
          CALL STME21 ( TIME , DTME21 )
          WRITE (NDSO,981) DTME21
          IF (FLLEV) THEN
              CALL STME21 ( TLEV , DTME21 )
              WRITE (NDSO,982) DTME21
            END IF
          IF (FLICE) THEN
              CALL STME21 ( TICE , DTME21 )
              WRITE (NDSO,983) DTME21
            END IF
!
          WRITE (NDSO,984)
!
        END IF
!
      IF ( NOPTS .EQ. 0 ) FLOUT(2) = .FALSE.
!
! 8.  Final MPI set up ----------------------------------------------- /
!
      CALL W3MPII ( IMOD )
      CALL W3MPIO ( IMOD )
      IF ( FLOUT(2) ) CALL W3MPIP ( IMOD )
!
      RETURN
!
! Escape locations read errors :
!
  820 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,8020) NSEA, NAPROC
      CALL EXTCDE ( 820 )
!
  821 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,8021) NSPEC, NAPROC
      CALL EXTCDE ( 821 )
!
  829 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,8029)
      CALL EXTCDE ( 829 )
 
!
  888 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,8000) IERR
      CALL EXTCDE ( 1 )
!
  889 CONTINUE
! === no process number filtering for test file !!! ===
      WRITE (NDSE,8001) IERR
      CALL EXTCDE ( 2 )
!
! Formats
!
  900 FORMAT ( ' WAVEWATCH III log file            ',                 &
               '                     version ',A/                     &
               ' ==================================',                 &
               '==================================='/                 &
               50X,'date : ',A10/50X,'time :  ',A8)
  920 FORMAT (/' Model definition file read.')
  930 FORMAT ( ' Restart file read; ',A)
!
  970 FORMAT (/' Grid name : ',A)
  971 FORMAT (/' ',A,' water levels.')
  972 FORMAT ( ' ',A,' curents.')
  973 FORMAT ( ' ',A,' winds.')
  974 FORMAT ( ' ',A,' ice fields.')
!
  975 FORMAT (/' Gridded output fields : '/                           &
               '--------------------------------------------------')
  976 FORMAT ( '     ',A)
!
  977 FORMAT (/' Point output requested for',I4,' points : '/         &
               '------------------------------------------')
  978 FORMAT (/'      Point output disabled')
  979 FORMAT                                                     &
        (/'      point  |  longitude  |   latitude  |  name  '/  &
     '     --------|-------------|-------------|----------------')
  980 FORMAT ( 5X,I5,'   |',2(F10.2,'   |'),2X,A)
!
  981 FORMAT (/' Initial time     : ',A)
  982 FORMAT ( ' Water level time : ',A)
  983 FORMAT ( ' Ice field time   : ',A)
!
  984 FORMAT (//                                                      &
        43X,'  |    input    |    output   |'/                        &
        43X,'  |-------------|-------------|'/                        &
         8X,'   step | pass |    date      time   | b w l c i d |',   &
              ' g p t r b f |'/                                       &
         8X,'--------|------|---------------------|-------------|',   &
              '-------------|'/                                       &
         8X,'--------+------+---------------------+-------------+',   &
              '-------------+')
!
 8000 FORMAT (/' *** WAVEWATCH III ERROR IN W3INIT : '/               &
               '     ERROR IN OPENING LOG FILE'/                      &
               '     IOSTAT =',I5/)
 8001 FORMAT (/' *** WAVEWATCH III ERROR IN W3INIT : '/               &
               '     ERROR IN OPENING TEST FILE'/                     &
               '     IOSTAT =',I5/)
 8002 FORMAT (/' *** WAVEWATCH III WARNING IN W3INIT : '/             &
               '     SIGNIFICANT PART OF RESOURCES RESERVED FOR',     &
                   ' OUTPUT :',F6.1,'%'/)
 8020 FORMAT (/' *** WAVEWATCH III ERROR IN W3INIT : '/         &
         '     NUMBER OF SEA POINTS LESS THAN NUMBER OF PROC.'/ &
         '     NSEA, NAPROC =',2I8/)
 8021 FORMAT (/' *** WAVEWATCH III ERROR IN W3INIT : '/         &
    '     NUMBER OF SPECTRAL POINTS LESS THAN NUMBER OF PROC.'/ &
         '     NSEA, NAPROC =',2I8/)
 8028 FORMAT (/' *** WAVEWATCH III WARNING IN W3INIT : '/       &
         '     INCREASING TARGET IN MPP PROPAGATION MAP.'/      &
         '     IMBALANCE BETWEEN OVERALL AND CFL TIME STEPS'/)
 8029 FORMAT (/' *** WAVEWATCH III ERROR IN W3INIT : '/         &
         '     SOMETHING WRONG WITH MPP PROPAGATION MAP.'/      &
         '     CALL HENDRIK !!!'/)
!
!/
!/ End of W3INIT ----------------------------------------------------- /
!/
      END SUBROUTINE W3INIT
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3MPII ( IMOD )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         11-May-2007 |
!/                  +-----------------------------------+
!/
!/    04-Jan-1999 : Distributed FORTRAN 77 version.     ( version 1.18 )
!/    13-Jan-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    28-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/                  Taken out of W3WAVE.
!/    04-May-2005 : Change to MPI_COMM_WAVE.            ( version 3.07 )
!/    13-Jun-2006 : Splitting STORE in G/SSTORE.        ( version 3.09 )
!/    11-May-2007 : Adding NTPROC/NAPROC separation.    ( version 3.11 )
!/
!  1. Purpose :
!
!     Perform initializations for MPI version of model.
!     Data transpose only.
!
!  2. Method :
!
!     Some derived data types are defined.  All communiction in
!     W3GATH, W3SCAT and W3WAVE are initialized so that all
!     communication can be performed with single MPI_STARTALL,
!     MPI_TESTALL and MPI_WAITALL calls.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.   I   Model number.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!
!      MPI_TYPE_VECTOR, MPI_TYPE_COMMIT
!                Subr. mpif.h   MPI derived data type routines.
!
!      MPI_SEND_INIT, MPI_RECV_INIT
!                Subr. mpif.h   MPI persistent communication calls.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3INIT    Subr. W3INITMD Wave model initialization routine.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     None.
!
!  7. Remarks :
!
!     - Basic MPP set up partially performed in W3INIT.
!     - Each processor has to be able to send out individual error
!       messages in this routine !
!     - No testing on IMOD, since only called by W3INIT.
!     - In version 3.09 STORE was split into a send and receive
!       buffer, to avoid/reduce possible conflicts btween the FORTRAN
!       and MPI standards when a gather is posted in a given buffer
!       right after a send in completed.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!       !/SHRD  Switch for shared / distributed memory architecture.
!       !/DIST  Id.
!       !/MPI   MPI communication calls.
!
!       !/S     Subroutine tracing,
!       !/T     Test output, general.
!       !/MPIT  Test output, MPI communications details.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE W3GDATMD, ONLY: NSEA
      USE W3GDATMD, ONLY: NSPEC
      USE W3WDATMD, ONLY: VA
      USE W3ADATMD, ONLY: MPI_COMM_WAVE, WW3_FIELD_VEC,         &
                          WW3_SPEC_VEC, IAPPRO, WADATS,         &
                          NRQSG1, IRQSG1, NRQSG2, IRQSG2,       &
                          GSTORE, SSTORE, MPIBUF, BSTAT,        &
                          BISPL, ISPLOC, IBFLOC, NSPLOC
      USE W3ODATMD, ONLY: NDST, NAPROC, IAPROC
!/
      IMPLICIT NONE
!
      INCLUDE "mpif.h"
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: IMOD
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: NSEALM
      INTEGER                 :: IERR_MPI, ISP, IH, ITARG,       &
                                 IERR1, IERR2, IP
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Set up derived data types -------------------------------------- /
!
      NSEALM = 1 + (NSEA-1)/NAPROC
!
      CALL MPI_TYPE_VECTOR ( NSEALM, 1, NAPROC, MPI_REAL,        &
                             WW3_FIELD_VEC, IERR_MPI )
      CALL MPI_TYPE_VECTOR ( NSEALM, 1, NSPEC, MPI_REAL,         &
                             WW3_SPEC_VEC, IERR_MPI )
      CALL MPI_TYPE_COMMIT ( WW3_FIELD_VEC, IERR_MPI )
      CALL MPI_TYPE_COMMIT ( WW3_SPEC_VEC, IERR_MPI )
!
      IF( IAPROC .GT. NAPROC ) THEN
          NSPLOC = 0
          NRQSG1 = 0
          NRQSG2 = 0
          RETURN
        END IF
!
! 2.  Set up scatters and gathers for W3WAVE ------------------------- /
!     ( persistent communication calls )
!
      NSPLOC = 0
      DO ISP=1, NSPEC
        IF ( IAPPRO(ISP) .EQ. IAPROC ) NSPLOC = NSPLOC + 1
        END DO
!
      NRQSG1 = NSPEC - NSPLOC
      ALLOCATE ( WADATS(IMOD)%IRQSG1(MAX(1,NRQSG1),2) )
      IRQSG1 => WADATS(IMOD)%IRQSG1
      IH     = 0
!
      DO ISP=1, NSPEC
        IF ( IAPPRO(ISP) .NE. IAPROC ) THEN
            ITARG  = IAPPRO(ISP) - 1
            IH     = IH + 1
            CALL MPI_SEND_INIT ( VA(ISP,1), 1, WW3_SPEC_VEC,     &
                 ITARG, ISP, MPI_COMM_WAVE, IRQSG1(IH,1), IERR1 )
            CALL MPI_RECV_INIT ( VA(ISP,1), 1, WW3_SPEC_VEC,     &
                 ITARG, ISP, MPI_COMM_WAVE, IRQSG1(IH,2), IERR2 )
          END IF
        END DO
!
! 3.  Set up scatters and gathers for W3SCAT and W3GATH -------------- /
!     Also set up buffering of data.
!
      NRQSG2 = MAX( 1 , NAPROC-1 )
      ALLOCATE ( WADATS(IMOD)%IRQSG2(NRQSG2*NSPLOC,2),           &
                 WADATS(IMOD)%GSTORE(NAPROC*NSEALM,MPIBUF),      &
                 WADATS(IMOD)%SSTORE(NAPROC*NSEALM,MPIBUF) )
      NRQSG2 = NAPROC - 1
!
      IRQSG2 => WADATS(IMOD)%IRQSG2
      GSTORE => WADATS(IMOD)%GSTORE
      SSTORE => WADATS(IMOD)%SSTORE
!
      IH     = 0
      ISPLOC = 0
      IBFLOC = 0
      WADATS(IMOD)%GSTORE = 0.
      WADATS(IMOD)%SSTORE = 0.
!
! 3.a Loop over local spectral components
!
      DO ISP=1, NSPEC
        IF ( IAPPRO(ISP) .EQ. IAPROC ) THEN
!
            ISPLOC = ISPLOC + 1
            IBFLOC = IBFLOC + 1
            IF ( IBFLOC .GT. MPIBUF ) IBFLOC = 1
!
! 3.b Loop over non-local processes
!
            DO IP=1, NAPROC
              IF ( IP .NE. IAPROC ) THEN
!
                  ITARG  = IP - 1
                  IH     = IH + 1
!
                  CALL MPI_RECV_INIT                             &
                     ( WADATS(IMOD)%GSTORE(IP,IBFLOC), 1,        &
                       WW3_FIELD_VEC, ITARG, ISP, MPI_COMM_WAVE, &
                       IRQSG2(IH,1), IERR2 )
                  CALL MPI_SEND_INIT                             &
                     ( WADATS(IMOD)%SSTORE(IP,IBFLOC), 1,        &
                       WW3_FIELD_VEC, ITARG, ISP, MPI_COMM_WAVE, &
                       IRQSG2(IH,2), IERR2 )
!
! ... End of loops
!
                END IF
              END DO
!
          END IF
        END DO
!
! 4.  Initialize buffer management ----------------------------------- /
!
      BSTAT  = 0
      BISPL  = 0
      ISPLOC = 0
      IBFLOC = 0
!
      RETURN
!
! Format statements
!
!/
!/ End of W3MPII ----------------------------------------------------- /
!/
      END SUBROUTINE W3MPII
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3MPIO ( IMOD )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         21-Jun-2007 |
!/                  +-----------------------------------+
!/
!/    17-Mar-1999 : Distributed FORTRAN 77 version.     ( version 1.18 )
!/    11-Jan-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    20-Aug-2003 : Output server options added.        ( version 3.04 )
!/    28-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/                  Taken out of W3WAVE.
!/    03-Jan-2005 : Add US2x to MPI communication.      ( version 3.06 )
!/    04-May-2005 : Change to MPI_COMM_WAVE.            ( version 3.07 )
!/    21-Jul-2005 : Add output fields.                  ( version 3.07 )
!/    04-Jul-2006 : Consolidate stress arrays.          ( version 3.09 )
!/    02-Aug-2006 : W3MPIP split off.                   ( version 3.10 )
!/    02-Apr-2007 : Add partitioned field data.         ( version 3.11 )
!/                  Add user-defined field data.
!/    17-May-2007 : Adding NTPROC/NAPROC separation.    ( version 3.11 )
!/    21-Jun-2007 : Dedicated output processes.         ( version 3.11 )
!/
!  1. Purpose :
!
!     Prepare MPI persistent communication needed for WAVEWATCH I/O
!     routines.
!
!  2. Method :
!
!     Create handles as needed.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.   I   Model number.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!
!      MPI_SEND_INIT, MPI_RECV_INIT
!                Subr. mpif.h   MPI persistent communication calls.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3INIT    Subr. W3INITMD Wave model initialization routine.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!  7. Remarks :
!
!     - The communication as set up in W3MPII uses tags with number
!       ranging from 1 through NSPEC. New and unique tags for IO
!       related communication are assigned here dynamically.
!     - No testing on IMOD, since only called by W3INIT.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!       !/MPI   MPI communication calls.
!
!       !/S     Enable subroutine tracing.
!       !/MPIT  Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE W3SERVMD, ONLY: EXTCDE
!/
      USE W3GDATMD, ONLY: NSEA
      USE W3GDATMD, ONLY: NX, NSPEC, MAPFS
      USE W3WDATMD, ONLY: VA, UST, USTDIR, ASF, FPIS
      USE W3ADATMD, ONLY: MPI_COMM_WAVE, WW3_FIELD_VEC, HS, WLM, &
                          TMN, THM, THS, FP0, THP0, FP1, THP1,   &
                          DTDYN, FCUT, SPPNT, ABA, ABD, UBA, UBD,&
                          SXX, SYY, SXY, USERO, PHS, PTP, PLP,   &
                          PTH, PSI, PWS, PWST, PNR
      USE W3ODATMD, ONLY: NDST, IAPROC, NAPROC, NTPROC, FLOUT,   &
                          NOGRD, NAPFLD, NAPPNT, NAPRST, NAPBPT, &
                          NAPTRK
      USE W3ODATMD, ONLY: OUTPTS, NRQGO, IRQGO, FLOGRD, NRQPO,   &
                          NRQPO2, IRQPO1, IRQPO2, NOPTS, IPTINT, &
                          NRQRS, IRQRS, NBLKRS, RSBLKS, IRQRSS,  &
                          VAAUX, NRQBP, NRQBP2, IRQBP1, IRQBP2,  &
                          NFBPO, NBO2, ISBPO, ABPOS, NRQTR,      &
                          IRQTR, IT0PNT, IT0TRK, IT0PRT, NOSWLL, &
                          NOEXTR, NDSE, IOSTYP
!/
      IMPLICIT NONE
!
      INCLUDE "mpif.h"
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: IMOD
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: NSEALM
      INTEGER                 :: IH, IT0, IROOT, IT, IERR, I0,   &
                                 IFROM, IX(4), IY(4), IS(4),     &
                                 IP(4), I, J, JSEA, ITARG, IB,   &
                                 JSEA0, JSEAN, NSEAB, IBOFF,     &
                                 ISEA, ISPROC, K, NRQMAX
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Set-up for W3IOGO ---------------------------------------------- /
!
      NRQGO  = 0
      IH     = 0
      IT0    = NSPEC
      IROOT  = NAPFLD - 1
!
      IF ( FLOUT(1) ) THEN
          IF ( NTPROC .EQ. 1 ) THEN
              NRQMAX = 0
            ELSE
              NRQMAX = 29 + NOEXTR + 6*NOSWLL
              IF ( IAPROC .EQ. NAPFLD ) THEN
                  IF ( IAPROC .LE. NAPROC ) THEN
                      NRQMAX = NRQMAX * (NAPROC-1)
                    ELSE
                      NRQMAX = NRQMAX * NAPROC
                    END IF
                END IF
            END IF
          IF ( NRQMAX .GT. 0 ) ALLOCATE                          &
                           ( OUTPTS(IMOD)%OUT1%IRQGO(NRQMAX) )
          IRQGO  => OUTPTS(IMOD)%OUT1%IRQGO
!
      IF ( IAPROC.NE.NAPFLD .AND. IAPROC.LE.NAPROC ) THEN
          IT     = IT0
!
          IF ( FLOGRD( 5) ) THEN
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (UST  (IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (USTDIR(IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (ASF  (IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
            END IF
          IF ( FLOGRD( 6) ) THEN
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (HS   (IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
            END IF
          IF ( FLOGRD( 7) ) THEN
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (WLM  (IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
            END IF
          IF ( FLOGRD( 8) ) THEN
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (TMN  (IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
            END IF
          IF ( FLOGRD( 9) ) THEN
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (THM  (IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
            END IF
          IF ( FLOGRD(10) ) THEN
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (THS  (IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
            END IF
          IF ( FLOGRD(11) ) THEN
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (FP0  (IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
            END IF
          IF ( FLOGRD(12) ) THEN
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (THP0 (IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
            END IF
          IF ( FLOGRD(13) ) THEN
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (FP1  (IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
            END IF
          IF ( FLOGRD(14) ) THEN
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (THP1 (IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
            END IF
          IF ( FLOGRD(15) ) THEN
            DO K=0, NOSWLL
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (PHS(IAPROC,K),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
              END DO
            END IF
          IF ( FLOGRD(16) ) THEN
            DO K=0, NOSWLL
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (PTP(IAPROC,K),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
              END DO
            END IF
          IF ( FLOGRD(17) ) THEN
            DO K=0, NOSWLL
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (PLP(IAPROC,K),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
              END DO
            END IF
          IF ( FLOGRD(18) ) THEN
            DO K=0, NOSWLL
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (PTH(IAPROC,K),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
              END DO
            END IF
          IF ( FLOGRD(19) ) THEN
            DO K=0, NOSWLL
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (PSI(IAPROC,K),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
              END DO
            END IF
          IF ( FLOGRD(20) ) THEN
            DO K=0, NOSWLL
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (PWS(IAPROC,K),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
              END DO
            END IF
          IF ( FLOGRD(21) ) THEN
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (PWST (IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
            END IF
          IF ( FLOGRD(22) ) THEN
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (PNR  (IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
            END IF
          IF ( FLOGRD(23) ) THEN
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (DTDYN(IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
            END IF
          IF ( FLOGRD(24) ) THEN
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (FCUT (IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
            END IF
          IF ( FLOGRD(27) ) THEN
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (ABA  (IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (ABD  (IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
            END IF
          IF ( FLOGRD(28) ) THEN
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (UBA  (IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (UBD  (IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
            END IF
          IF ( FLOGRD(29) ) THEN
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (SXX  (IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (SYY  (IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (SXY  (IAPROC),1,WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
            END IF
          IF ( FLOGRD(30) ) THEN
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (USERO(IAPROC,1),1,WW3_FIELD_VEC,&
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
            END IF
          IF ( FLOGRD(31) ) THEN
              IH     = IH + 1
              IT     = IT + 1
              CALL MPI_SEND_INIT (USERO(IAPROC,2),1,WW3_FIELD_VEC,&
                       IROOT, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
            END IF
!
        ELSE IF ( IAPROC .EQ. NAPFLD ) THEN
!
          DO I0=1, NAPROC
            IT     = IT0
            IF ( I0 .NE. IAPROC ) THEN
                IFROM  = I0 - 1
                IF ( FLOGRD( 5) ) THEN
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (UST  (I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (USTDIR(I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (ASF  (I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                  END IF
                IF ( FLOGRD( 6) ) THEN
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (HS   (I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                  END IF
                IF ( FLOGRD( 7) ) THEN
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (WLM  (I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                  END IF
                IF ( FLOGRD( 8) ) THEN
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (TMN  (I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                  END IF
                IF ( FLOGRD( 9) ) THEN
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (THM  (I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                  END IF
                IF ( FLOGRD(10) ) THEN
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (THS  (I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                  END IF
                IF ( FLOGRD(11) ) THEN
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (FP0  (I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                  END IF
                IF ( FLOGRD(12) ) THEN
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (THP0 (I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                  END IF
                IF ( FLOGRD(13) ) THEN
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (FP1  (I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                  END IF
                IF ( FLOGRD(14) ) THEN
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (THP1 (I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                  END IF
                IF ( FLOGRD(15) ) THEN
                  DO K=0, NOSWLL
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (PHS(I0,K),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                    END DO
                  END IF
                IF ( FLOGRD(16) ) THEN
                  DO K=0, NOSWLL
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (PTP(I0,K),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                    END DO
                  END IF
                IF ( FLOGRD(17) ) THEN
                  DO K=0, NOSWLL
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (PLP(I0,K),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                    END DO
                  END IF
                IF ( FLOGRD(18) ) THEN
                  DO K=0, NOSWLL
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (PTH(I0,K),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                    END DO
                  END IF
                IF ( FLOGRD(19) ) THEN
                  DO K=0, NOSWLL
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (PSI(I0,K),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                    END DO
                  END IF
                IF ( FLOGRD(20) ) THEN
                  DO K=0, NOSWLL
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (PWS(I0,K),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                    END DO
                  END IF
                IF ( FLOGRD(21) ) THEN
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (PWST (I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                  END IF
                IF ( FLOGRD(22) ) THEN
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (PNR  (I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                  END IF
                IF ( FLOGRD(23) ) THEN
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (DTDYN(I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                  END IF
                IF ( FLOGRD(24) ) THEN
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (FCUT (I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                  END IF
                IF ( FLOGRD(27) ) THEN
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (ABA  (I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (ABD  (I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                  END IF
                IF ( FLOGRD(28) ) THEN
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (UBA  (I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (UBD  (I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                  END IF
                IF ( FLOGRD(29) ) THEN
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (SXX  (I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (SYY  (I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (SXY  (I0),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                  END IF
                IF ( FLOGRD(30) ) THEN
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (USERO(I0,1),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                  END IF
                IF ( FLOGRD(31) ) THEN
                    IH     = IH + 1
                    IT     = IT + 1
                    CALL MPI_RECV_INIT (USERO(I0,2),1,WW3_FIELD_VEC,&
                       IFROM, IT, MPI_COMM_WAVE, IRQGO(IH), IERR )
                  END IF
              END IF
            END DO
!
        END IF
!
      NRQGO  = IH
      IT0    = IT0 + ( 29 + NOEXTR + 6*NOSWLL ) * NAPROC
!
      IF ( NRQGO .GT. NRQMAX ) THEN
          WRITE (NDSE,1010)
          CALL EXTCDE (10)
        END IF
!
      END IF
!
! 2.  Set-up for W3IORS ---------------------------------------------- /
! 2.a General preparations
!
      NRQRS  = 0
      IH     = 0
      IROOT  = NAPRST - 1
!
      IF ( FLOUT(4) ) THEN
          ALLOCATE ( OUTPTS(IMOD)%OUT4%IRQRS(2*NSPEC) )
          IRQRS  => OUTPTS(IMOD)%OUT4%IRQRS
!
! 2.b Fields at end of file (allways)
!
          IF ( IAPROC.NE.NAPRST .AND. IAPROC.LE.NAPROC ) THEN
!
              IH     = IH + 1
              IT     = IT0 + 1
              CALL MPI_SEND_INIT (UST (IAPROC), 1, WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQRS(IH), IERR )
!
              IH     = IH + 1
              IT     = IT0 + 2
              CALL MPI_SEND_INIT (USTDIR(IAPROC), 1, WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQRS(IH), IERR )
!
              IH     = IH + 1
              IT     = IT0 + 3
              CALL MPI_SEND_INIT (FPIS(IAPROC), 1, WW3_FIELD_VEC, &
                       IROOT, IT, MPI_COMM_WAVE, IRQRS(IH), IERR )
!
            ELSE IF ( IAPROC .EQ. NAPRST ) THEN
              DO I0=1, NAPROC
                IFROM  = I0 - 1
                IF ( I0 .NE. IAPROC ) THEN
!
                    IH     = IH + 1
                    IT     = IT0 + 1
                    CALL MPI_RECV_INIT (UST (I0),1,WW3_FIELD_VEC, &
                       IFROM, IT, MPI_COMM_WAVE, IRQRS(IH), IERR )
!
                    IH     = IH + 1
                    IT     = IT0 + 2
                    CALL MPI_RECV_INIT (USTDIR(I0),1,WW3_FIELD_VEC, &
                       IFROM, IT, MPI_COMM_WAVE, IRQRS(IH), IERR )
!
                    IH     = IH + 1
                    IT     = IT0 + 3
                    CALL MPI_RECV_INIT (FPIS(I0),1,WW3_FIELD_VEC, &
                       IFROM, IT, MPI_COMM_WAVE, IRQRS(IH), IERR )
!
                  END IF
!
                END DO
            END IF
!
          NRQRS  = IH
          IT0    = IT0 + 3
!
! 2.c Data server mode
!
          IF ( IOSTYP .GT. 0 ) THEN
!
              NSEALM = 1 + (NSEA-1)/NAPROC
              NBLKRS = 10
              RSBLKS = MAX ( 5 , NSEALM/NBLKRS )
              IF ( NBLKRS*RSBLKS .LT. NSEALM ) RSBLKS = RSBLKS + 1
              NBLKRS = 1 + (NSEALM-1)/RSBLKS
!
              IH     = 0
!
              IF ( IAPROC .NE. NAPRST ) THEN
!
                  ALLOCATE ( OUTPTS(IMOD)%OUT4%IRQRSS(NBLKRS) )
                  IRQRSS => OUTPTS(IMOD)%OUT4%IRQRSS
!
                  DO IB=1, NBLKRS
                    IH     = IH + 1
                    IT     = IT0 + 3 + IB
                    JSEA0  = 1 + (IB-1)*RSBLKS
                    JSEAN  = MIN ( NSEALM , IB*RSBLKS )
                    NSEAB  = 1 + JSEAN - JSEA0
                    CALL MPI_SEND_INIT (VA(1,JSEA0), NSPEC*NSEAB,&
                         MPI_REAL, IROOT, IT, MPI_COMM_WAVE,    &
                         IRQRSS(IH), IERR )
                    END DO
!
                ELSE
!
                  ALLOCATE                                       &
                 ( OUTPTS(IMOD)%OUT4%IRQRSS(NAPROC*NBLKRS) ,     &
                   OUTPTS(IMOD)%OUT4%VAAUX(NSPEC,2*RSBLKS,NAPROC) )
!
                  IRQRSS => OUTPTS(IMOD)%OUT4%IRQRSS
                  VAAUX  => OUTPTS(IMOD)%OUT4%VAAUX
                  DO IB=1, NBLKRS
                    IT     = IT0 + 3 + IB
                    JSEA0  = 1 + (IB-1)*RSBLKS
                    JSEAN  = MIN ( NSEALM , IB*RSBLKS )
                    NSEAB  = 1 + JSEAN - JSEA0
                    DO I0=1, NAPROC
                      IF ( I0 .NE. NAPRST ) THEN
                          IH     = IH + 1
                          IFROM  = I0 - 1
                          IBOFF  = MOD(IB-1,2)*RSBLKS
                          CALL MPI_RECV_INIT (VAAUX(1,1+IBOFF,I0),&
                               NSPEC*NSEAB, MPI_REAL, IFROM, IT,  &
                               MPI_COMM_WAVE, IRQRSS(IH), IERR )
                        END IF
                      END DO
                    END DO
!
                END IF
!
              IT0    = IT0 + NBLKRS
!
            END IF
!
        END IF
!
! 3.  Set-up for W3IOBC ( SENDs ) ------------------------------------ /
!
      NRQBP  = 0
      NRQBP2 = 0
      IH     = 0
      IT     = IT0
      IROOT  = NAPBPT - 1
!
      IF ( FLOUT(5) ) THEN
          ALLOCATE ( OUTPTS(IMOD)%OUT5%IRQBP1(NBO2(NFBPO)),      &
                     OUTPTS(IMOD)%OUT5%IRQBP2(NBO2(NFBPO)) )
          IRQBP1 => OUTPTS(IMOD)%OUT5%IRQBP1
          IRQBP2 => OUTPTS(IMOD)%OUT5%IRQBP2
!
! 3.a Loops over files and points
!
          DO J=1, NFBPO
            DO I=NBO2(J-1)+1, NBO2(J)
!
               IT     = IT + 1
!
! 3.b Residence processor of point
!
              ISEA   = ISBPO(I)
              JSEA   = 1 + (ISEA-1)/NAPROC
              ISPROC = ISEA - (JSEA-1)*NAPROC
!
! 3.c If stored locally, send data
!
              IF ( IAPROC .EQ. ISPROC ) THEN
                  IH     = IH + 1
                  CALL MPI_SEND_INIT (VA(1,JSEA),NSPEC,MPI_REAL, &
                       IROOT, IT, MPI_COMM_WAVE, IRQBP1(IH), IERR)
                END IF
!
              END DO
            END DO
!
! ... End of loops 4.a
!
          NRQBP  = IH
!
! 3.d Set-up for W3IOBC ( RECVs ) ------------------------------------ /
!
          IF ( IAPROC .EQ. NAPBPT ) THEN
!
              IH     = 0
              IT     = IT0
!
! 3.e Loops over files and points
!
              DO J=1, NFBPO
                DO I=NBO2(J-1)+1, NBO2(J)
!
! 3.f Residence processor of point
!
                  ISEA   = ISBPO(I)
                  JSEA   = 1 + (ISEA-1)/NAPROC
                  ISPROC = ISEA - (JSEA-1)*NAPROC
!
! 3.g Receive in correct array
!
                  IH     = IH + 1
                  IT     = IT + 1
                  ITARG  = ISPROC - 1
                  CALL MPI_RECV_INIT (ABPOS(1,IH),NSPEC,MPI_REAL,&
                       ITARG, IT, MPI_COMM_WAVE, IRQBP2(IH), IERR)
!
                  END DO
                END DO
!
              NRQBP2 = IH
!
! ... End of loops 4.e
!
            END IF
!
          IT0    = IT0 + NBO2(NFBPO)
!
        END IF
!
! 4.  Set-up for W3IOTR ---------------------------------------------- /
!
      IH     = 0
      IROOT  = NAPTRK - 1
!
      IF ( FLOUT(3) ) THEN
!
! 4.a U*
!
          IF ( IAPROC .NE. NAPTRK ) THEN
              ALLOCATE ( OUTPTS(IMOD)%OUT3%IRQTR(2) )
              IRQTR  => OUTPTS(IMOD)%OUT3%IRQTR
              IH     = IH + 1
              IT     = IT0 + 1
              CALL MPI_SEND_INIT (UST   (IAPROC),1,WW3_FIELD_VEC,&
                   IROOT, IT, MPI_COMM_WAVE, IRQTR(IH), IERR )
              IH     = IH + 1
              IT     = IT0 + 2
              CALL MPI_SEND_INIT (USTDIR(IAPROC),1,WW3_FIELD_VEC,&
                   IROOT, IT, MPI_COMM_WAVE, IRQTR(IH), IERR )
            ELSE
              ALLOCATE ( OUTPTS(IMOD)%OUT3%IRQTR(2*NAPROC) )
              IRQTR  => OUTPTS(IMOD)%OUT3%IRQTR
              DO I0=1, NAPROC
                IFROM  = I0 - 1
                IF ( I0 .NE. IAPROC ) THEN
                    IH     = IH + 1
                    IT     = IT0 + 1
                    CALL MPI_RECV_INIT(UST   (I0),1,WW3_FIELD_VEC,&
                         IFROM,IT,MPI_COMM_WAVE, IRQTR(IH), IERR)
                    IH     = IH + 1
                    IT     = IT0 + 2
                    CALL MPI_RECV_INIT(USTDIR(I0),1,WW3_FIELD_VEC,&
                         IFROM,IT,MPI_COMM_WAVE, IRQTR(IH), IERR)
                  END IF
                END DO
            END IF
!
          NRQTR  = IH
          IT0    = IT0 + 2
!
        END IF
!
! 5.  Set-up remaining counters -------------------------------------- /
!
      IT0PRT = IT0
      IT0PNT = IT0PRT + 2*NAPROC
      IT0TRK = IT0PNT + 5000
!
      RETURN
!
!     Formats :
!
  1010 FORMAT (/' *** ERROR W3MPIO : ARRAY IRQGO TOO SMALL *** '/)
!
!/
!/ End of W3MPIO ----------------------------------------------------- /
!/
      END SUBROUTINE W3MPIO
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3MPIP ( IMOD )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         17-May-2007 |
!/                  +-----------------------------------+
!/
!/    02-Aug-2006 : Origination.                        ( version 3.10 )
!/    17-May-2007 : Adding NTPROC/NAPROC separation.    ( version 3.11 )
!/
!  1. Purpose :
!
!     Prepare MPI persistent communication needed for WAVEWATCH I/O
!     routines.
!
!  2. Method :
!
!     Create handles as needed.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.   I   Model number.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!
!      MPI_SEND_INIT, MPI_RECV_INIT
!                Subr. mpif.h   MPI persistent communication calls.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3INIT    Subr. W3INITMD Wave model initialization routine.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!       !/MPI   MPI communication calls.
!
!       !/S     Enable subroutine tracing.
!       !/MPIT  Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3SERVMD, ONLY: EXTCDE
!/
      USE W3GDATMD, ONLY: NX, NSPEC, MAPFS
      USE W3WDATMD, ONLY: VA
      USE W3ADATMD, ONLY: MPI_COMM_WAVE, SPPNT
      USE W3ODATMD, ONLY: NDST, NDSE, IAPROC, NAPROC, NAPPNT, FLOUT
      USE W3ODATMD, ONLY: OUTPTS, NRQPO, NRQPO2, IRQPO1, IRQPO2, &
                          NOPTS, IPTINT, IT0PNT, IT0TRK, O2IRQI
!/
      IMPLICIT NONE
!
      INCLUDE "mpif.h"
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: IMOD
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IH, IROOT, I, J, IT, IT0, JSEA, &
                                 IERR, ITARG, IX(4), IY(4),      &
                                 IS(4), IP(4)
!/
!/ ------------------------------------------------------------------- /
!/
!
      IF ( O2IRQI ) THEN
          WRITE (NDSE,1001)
          CALL EXTCDE (1)
        END IF
!
! 1.  Set-up for W3IOPE/O ( SENDs ) ---------------------------------- /
!
      NRQPO  = 0
      NRQPO2 = 0
      IH     = 0
      IT0    = IT0PNT
      IROOT  = NAPPNT - 1
!
      ALLOCATE ( OUTPTS(IMOD)%OUT2%IRQPO1(4*NOPTS),              &
                 OUTPTS(IMOD)%OUT2%IRQPO2(4*NOPTS) )
      IRQPO1 => OUTPTS(IMOD)%OUT2%IRQPO1
      IRQPO2 => OUTPTS(IMOD)%OUT2%IRQPO2
      O2IRQI = .TRUE.
!
! 1.a Loop over output locations
!
      DO I=1, NOPTS
!
        IX(1)  = IPTINT(1,I)
        IX(2)  = 1 + MOD(IPTINT(1,I),NX)
        IX(3)  = IPTINT(1,I)
        IX(4)  = IX(2)
!
        IY(1)  = IPTINT(2,I)
        IY(2)  = IPTINT(2,I)
        IY(3)  = IPTINT(2,I) + 1
        IY(4)  = IPTINT(2,I) + 1
!
! 1.b Loop over corner points
!
        DO J=1, 4
!
          IT     = IT0 + (I-1)*4 + J
          IS(J)  = MAPFS (IY(J),IX(J))
          IF ( IS(J) .EQ. 0 ) THEN
              JSEA   = 0
              IP(J)  = NAPPNT
            ELSE
              JSEA   = 1 + (IS(J)-1)/NAPROC
              IP(J)  = IS(J) - (JSEA-1)*NAPROC
            END IF
!
! 1.c Send if point is stored here
!
          IF ( IP(J) .EQ. IAPROC ) THEN
              IH     = IH + 1
              CALL MPI_SEND_INIT ( VA(1,JSEA), NSPEC, MPI_REAL, &
                   IROOT, IT, MPI_COMM_WAVE, IRQPO1(IH), IERR )
            END IF
!
! ... End of loop 1.b
!
          END DO
!
! ... End of loop 1.a
!
        END DO
!
      NRQPO  = IH
!
! 1.d Set-up for W3IOPE/O ( RECVs ) ---------------------------------- /
!
      IF ( IAPROC .EQ. NAPPNT ) THEN
!
          IH     = 0
!
! 2.e Loop over output locations
!
          DO I=1, NOPTS
!
            IX(1)  = IPTINT(1,I)
            IX(2)  = 1 + MOD(IPTINT(1,I),NX)
            IX(3)  = IPTINT(1,I)
            IX(4)  = IX(2)
!
            IY(1)  = IPTINT(2,I)
            IY(2)  = IPTINT(2,I)
            IY(3)  = IPTINT(2,I) + 1
            IY(4)  = IPTINT(2,I) + 1
!
            DO J=1, 4
!
              IT     = IT0 + (I-1)*4 + J
              IS(J)  = MAPFS (IY(J),IX(J))
              IF ( IS(J) .EQ. 0 ) THEN
                  JSEA   = 0
                  IP(J)  = NAPPNT
                ELSE
                  JSEA   = 1 + (IS(J)-1)/NAPROC
                  IP(J)  = IS(J) - (JSEA-1)*NAPROC
                END IF
!
! 1.g Receive in correct array
!
              IH     = IH + 1
              ITARG  = IP(J) - 1
              CALL MPI_RECV_INIT ( SPPNT(1,1,J), NSPEC, MPI_REAL, &
                   ITARG, IT, MPI_COMM_WAVE, IRQPO2(IH), IERR )
!
! ... End of loop 1.f
!
              END DO
!
! ... End of loop 1.e
!
            END DO
!
          NRQPO2 = NOPTS*4
!
        END IF
!
      IT0    = IT0 + 8*NOPTS
!
! 1.h Base tag number for track output
!
      IT0TRK = IT0
!
      RETURN
!
!     Formats :
!
  1001 FORMAT (/' *** ERROR W3MPIP : ARRAYS ALREADY ALLOCATED *** '/)
!
!/
!/ End of W3MPIP ----------------------------------------------------- /
!/
      END SUBROUTINE W3MPIP
!/
!/ End of module W3INITMD -------------------------------------------- /
!/
      END MODULE W3INITMD
