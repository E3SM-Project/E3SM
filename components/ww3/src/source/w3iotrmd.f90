!/ ------------------------------------------------------------------- /
      MODULE W3IOTRMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!
!/    See subroutine for update history.
!/
!  1. Purpose :
!
!     Generate track output.
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      VERTRK    C*10  Private  Version number of routine.
!      IDSTRI    C*34  Private  ID string input file.
!      IDSTRO    C*34  Private  ID string output file.
!     ----------------------------------------------------------------
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3IOTR    Subr. Public   Track output subroutine.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SETO    Subr. W3ODATMD Point to data structure.
!      W3SETG    Subr. W3GDATMD Point to data structure.
!      W3SETW    Subr. W3WDATMD Point to data structure.
!      W3SETA    Subr. W3ADATMD Point to data structure.
!      W3DMO3    Subr. W3ODATMD Allocate work arrays.
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!      TICK21    Subr. W3TIMEMD Increment time.
!      DSEC21    Func. W3TIMEMD Time difference.
!      MPI_SEND, MPI_RECV, MPI_STARTALL, MPI_WAITALL
!                               MPI send and recieve routines
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!  6. Switches :
!
!     See documentation of W3IOTR.
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
!/
!/ Private parameter statements (ID strings)
!/
      CHARACTER(LEN=10), PARAMETER, PRIVATE :: VERTRK = 'III  1.02 '
      CHARACTER(LEN=34), PARAMETER, PRIVATE ::                        &
                       IDSTRI = 'WAVEWATCH III TRACK LOCATIONS DATA', &
                       IDSTRO = 'WAVEWATCH III TRACK OUTPUT SPECTRA'
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3IOTR ( NDSINP, NDSOUT, A, IMOD )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    22-Dec-1998 : Final FORTRAN 77                    ( version 1.18 )
!/    27-Dec-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    24-Jan-2001 : Flat grid version                   ( version 2.06 )
!/    20-Aug-2003 : Output through NAPTRK, seq. file.   ( version 3.04 )
!/    24-Nov-2004 : Multiple grid version.              ( version 3.06 )
!/    04-May-2005 : Change to MPI_COMM_WAVE.            ( version 3.07 )
!/    27-Jun-2005 : Adding MAPST2,                      ( version 3.07 )
!/    27-Jun-2006 : Adding file name preamble.          ( version 3.09 )
!/    17-May-2007 : Adding NTPROC/NAPROC separation.    ( version 3.11 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Perform output of spectral information along provided tracks.
!
!  2. Method :
!
!     Time and location data for the track is read from the file
!     track_i.FILEXT, and output spectra ad adtional information are
!     written to track_o.FILEXT.
!
!     The spectrum dumped is the frequency-direction spectrum in
!     m**2/Hz/rad.
!
!     The output spectra are energy density spectra in terms of the
!     true frequency and a direction in radians. The corresponding
!     band widths are part of the file header.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NDSINP  Int.   I   Unit number of input file track_i.FILEXT
!                          If negative, file is unformatted and v.v.
!       NDSOUT  Int.   I   Unit number of outptu file track_o.FILEXT
!       A       R.A.   I   Spectra (shape conversion through par list).
!       IMOD    Int.   I   Model grid number.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     See module documentation.
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3WAVE    Subr. W3WAVEMD Actual wave model routine.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     - If input file not found, a warning is printed and output
!       type is disabled.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/SHRD  Switch for shared / distributed memory architecture.
!     !/DIST  Id.
!     !/MPI   MPI interface routines.
!
!     !/LLG   Spherical grid.
!     !/XYG   Cartesian grid.
!
!     !/S     Enable subroutine tracing.
!     !/T     General test output.
!     !/T1    Test output on track point status.
!     !/T2    Test output of mask arrays.
!     !/T3    Test output for writing of file.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3CONSTANTS
!/
      USE W3GDATMD, ONLY: W3SETG
      USE W3WDATMD, ONLY: W3SETW
      USE W3ADATMD, ONLY: W3SETA
      USE W3ODATMD, ONLY: W3SETO, W3DMO3
!/
      USE W3GDATMD, ONLY: NK, NTH, NSPEC, NSEA, NSEAL, NX, NY,        &
                          X0, Y0, SX, SY, GLOBAL, MAPSTA, MAPST2,     &
                          MAPFS, TH, DTH, SIG, DSIP, FILEXT
      USE W3WDATMD, ONLY: TIME, UST
      USE W3ADATMD, ONLY: CG, DW, CX, CY, UA, UD, AS
      USE W3ADATMD, ONLY: MPI_COMM_WAVE
      USE W3ODATMD, ONLY: NDST, NDSE, IAPROC, NAPROC, NAPTRK, NAPERR, &
                          IPASS => IPASS3, ATOLAST => TOLAST,         &
                          ADTOUT => DTOUT, O3INIT, STOP, MASK1,       &
                          MASK2, FNMPRE
      USE W3ODATMD, ONLY: IT0TRK, NRQTR, IRQTR
!/
      USE W3TIMEMD
!
      IMPLICIT NONE
!
      INCLUDE "mpif.h"
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: NDSINP, NDSOUT, IMOD
      REAL, INTENT(IN)        :: A(NTH,NK,0:NSEAL)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER, PARAMETER      :: OTYPE = 3
      INTEGER                 :: NDSTI, NDSTO, ISPROC, IERR,          &
                                 IK, ITH, IX, IY, TIMEB(2), TIMEE(2), &
                                 TTIME(2), IX1, IX2, IY1, IY2,        &
                                 IXX(4), IYY(4), I, J, ISEA, JSEA,    &
                                 TOLAST(2)
      INTEGER                 :: IT, IROOT, IFROM, IERR_MPI
      INTEGER, ALLOCATABLE    :: STATUS(:,:)
      REAL                    :: XN, YN, XT, YT, RD, X, Y, WX, WY,    &
                                 SPEC(NK,NTH), FACTOR, ASPTRK(NTH,NK),&
                                 DTOUT
      REAL, SAVE              :: RDCHCK = 0.05, RTCHCK = 0.05
      LOGICAL                 :: FORMI, FLAG1, FLAG2
      CHARACTER(LEN=3)        :: TSTSTR
      CHARACTER(LEN=34)       :: IDTST
!
      EQUIVALENCE                (IXX(1),IX1) , (IXX(2),IX2) ,        &
                                 (IYY(1),IY1) , (IYY(3),IY2)
!/
!/ ------------------------------------------------------------------- /
!/
!
      CALL W3SETO ( IMOD, NDSE, NDST )
      CALL W3SETG ( IMOD, NDSE, NDST )
      CALL W3SETA ( IMOD, NDSE, NDST )
      CALL W3SETW ( IMOD, NDSE, NDST )
!
      TOLAST = ATOLAST(:,OTYPE)
      DTOUT  = ADTOUT(OTYPE)
!
      IF ( .NOT. O3INIT ) CALL W3DMO3 ( IMOD, NDSE, NDST )
!
      FORMI  = NDSINP .GT. 0
      NDSTI  = ABS(NDSINP)
      NDSTO  = ABS(NDSOUT)
!
      XN     = X0 + REAL(NX-1)*SX
      YN     = Y0 + REAL(NY-1)*SY
!
      ISPROC = IAPROC
      IPASS  = IPASS + 1
!
      FACTOR = 1.
!
      IF ( NRQTR .GT. 0 ) THEN
          CALL MPI_STARTALL ( NRQTR, IRQTR, IERR_MPI )
          ALLOCATE ( STATUS(MPI_STATUS_SIZE,NRQTR) )
          CALL MPI_WAITALL ( NRQTR, IRQTR , STATUS, IERR_MPI )
          DEALLOCATE ( STATUS )
        END IF
!
! 1.  First pass through routine ------------------------------------- *
!
      IF ( IPASS .EQ. 1 ) THEN
!
          I      = LEN_TRIM(FILEXT)
          J      = LEN_TRIM(FNMPRE)
!
! 1.a Open input file
!
          IF ( FORMI ) THEN
              OPEN (NDSTI,FILE=FNMPRE(:J)//'track_i.'//FILEXT(:I),     &
                    STATUS='OLD',ERR=800,FORM='FORMATTED',IOSTAT=IERR)
              READ (NDSTI,'(A)',ERR=801,END=802,IOSTAT=IERR) IDTST
            ELSE
              OPEN (NDSTI,FILE=FNMPRE(:J)//'track_i.'//FILEXT(:I),     &
                    STATUS='OLD',ERR=800,FORM='UNFORMATTED',IOSTAT=IERR)
              READ (NDSTI,ERR=801,END=802,IOSTAT=IERR) IDTST
            END IF
!
          IF ( IDTST .NE. IDSTRI ) GOTO 803
!
! 1.b Open output file
!
          IF ( IAPROC .EQ. NAPTRK ) THEN
              OPEN (NDSTO,FILE=FNMPRE(:J)//'track_o.'//FILEXT(:I),     &
                    FORM='UNFORMATTED',ERR=810,IOSTAT=IERR)
              WRITE (NDSTO,ERR=811,IOSTAT=IERR) IDSTRO, NK, NTH
              WRITE (NDSTO,ERR=811,IOSTAT=IERR) 0.5*PI-TH(1), -DTH,    &
                    (SIG(IK)*TPIINV,IK=1,NK),                          &
                    (DSIP(IK)*TPIINV,IK=1,NK)
            END IF
!
! 1.c Initialize maps
!
          MASK2  = .FALSE.
!
        END IF
!
! 2.  Preparations --------------------------------------------------- *
! 2.a Shift mask arrays
!
      MASK1  = MASK2
      MASK2  = .FALSE.
!
! 2.b Set time frame
!
      TIMEB  = TIME
      TIMEE  = TIME
      CALL TICK21 ( TIMEE ,  DTOUT )
!
      IF ( DSEC21(TIMEE,TOLAST) .LT. 0. ) THEN
          TIMEE  = TOLAST
        END IF
!
! 3.  Loop over input points ----------------------------------------- *
!
! 3.a Read new track point (infinite loop)
!
      IF ( STOP ) THEN
          TOLAST = TIME
          GOTO 399
        END IF
!
      DO
!
        IF ( FORMI ) THEN
            READ (NDSTI,*,ERR=801,END=390,IOSTAT=IERR) TTIME, XT, YT
          ELSE
            READ (NDSTI , ERR=801,END=390,IOSTAT=IERR) TTIME, XT, YT
          END IF
!
! 3.b Point before time interval
!
        IF ( DSEC21(TIMEB,TTIME) .LT. 0. ) THEN
            CYCLE
          END IF
!
! 3.c Point after time interval
!
        IF ( DSEC21(TIMEE,TTIME) .GT. 0. ) THEN
            BACKSPACE (NDSTI)
            GOTO 399
          END IF
!
! 3.d Check time in interval
!
        FLAG1  = DSEC21(TTIME,TIMEE) .GT. RTCHCK*DTOUT
        FLAG2  = DSEC21(TIMEB,TTIME) .GT. RTCHCK*DTOUT
!
! 3.e Check point coordinates
!
        IF ( XT .GT. X0 ) THEN
            XT     = XT - 360.*REAL(INT((XT-X0)/360.))
          ELSE
            XT     = XT + 360.*REAL(1+INT((X0-XT)/360.))
          END IF
!
        IF ( .NOT.GLOBAL .AND. XT.GT.XN ) THEN
            CYCLE
          END IF
!
        IF ( YT.LT.Y0 .OR. YT.GT.YN ) THEN
            CYCLE
          END IF
!
        IX1    = 1 + INT((XT-X0)/SX)
        RD     = (XT-X0-REAL(IX1-1)*SX)/SX
        IF ( RD .LT. RDCHCK ) THEN
            IX2    = IX1
          ELSE IF ( RD .GT. 1.-RDCHCK ) THEN
            IX1    = IX1 + 1
            IX2    = IX1
          ELSE
            IX2    = IX1 + 1
          END IF
!
        IF ( GLOBAL .AND. IX1.GT.NX ) IX1 = IX1 - NX
        IF ( GLOBAL .AND. IX2.GT.NX ) IX2 = IX2 - NX
!
        IY1    = 1 + INT((YT-Y0)/SY)
        RD     = (YT-Y0-REAL(IY1-1)*SY)/SY
        IF ( RD .LT. RDCHCK ) THEN
            IY2    = IY1
          ELSE IF ( RD .GT. 1.-RDCHCK ) THEN
            IY1    = IY1 + 1
            IY2    = IY1
          ELSE
            IY2    = IY1 + 1
          END IF
!
        IXX(3) = IX1
        IXX(4) = IX2
        IYY(2) = IY1
        IYY(4) = IY2
!
! 3.f Mark the four corner points
!
        DO J=1, 4
!
          IX     = IXX(J)
          IY     = IYY(J)
          X      = X0 + REAL(IX-1)*SX
          Y      = Y0 + REAL(IY-1)*SY
          MASK1(IY,IX) = MASK1(IY,IX) .OR. FLAG1
          MASK2(IY,IX) = MASK2(IY,IX) .OR. FLAG2
!
          END DO
!
        END DO
!
! 3.g End of input file escape location
!
  390 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1004)
      STOP   = .TRUE.
!
! 3.h Read end escape location
!
  399 CONTINUE
!
! 3.h Mask test output
!
! 4.  Write data for flagged locations ------------------------------- *
!
      IT     = IT0TRK
      IROOT  = NAPTRK - 1
      ALLOCATE ( STATUS(MPI_STATUS_SIZE,1) )
!
      DO IY=1, NY
        DO IX=1, NX
          IF ( MASK1(IY,IX) ) THEN
!
              X      = X0 + REAL(IX-1)*SX
              Y      = Y0 + REAL(IY-1)*SY
              IT     = IT + 1
!
! 4.a Status of point
!
              IF ( MAPSTA(IY,IX) .EQ. 0 ) THEN
                  IF ( MAPST2(IY,IX) .EQ. 0 ) THEN
                      TSTSTR = 'LND'
                    ELSE
                      TSTSTR = 'XCL'
                    END IF
                ELSE IF ( MAPSTA(IY,IX) .LT. 0 ) THEN
                  IF ( MAPST2(IY,IX) .EQ. 1 ) THEN
                      TSTSTR = 'ICE'
                    ELSE IF ( MAPST2(IY,IX) .EQ. 2 ) THEN
                      TSTSTR = 'DRY'
                    ELSE
                      TSTSTR = 'DIS'
                    END IF
                ELSE
                  TSTSTR = 'SEA'
                END IF
!
! 4.b Determine where point is stored
!     ( land point assumed stored on IAPROC = NAPTRK
!       set to -99 in test output )
!
              ISEA   = MAPFS(IY,IX)
              IF ( ISEA .EQ. 0 ) THEN
                  ISPROC = NAPTRK
                ELSE
                  JSEA   = 1 + (ISEA-1)/NAPROC
                  ISPROC = ISEA - (JSEA-1)*NAPROC
                END IF
              IFROM  = ISPROC - 1
! 4.c Spectrum is at local processor, but this is not the NAPTRK
!     Send the spectrum to NAPTRK
 
              IF ( ISPROC.EQ.IAPROC .AND. IAPROC.NE.NAPTRK ) THEN
                  CALL MPI_SEND ( A(1,1,JSEA), NSPEC, MPI_REAL,  &
                           IROOT, IT, MPI_COMM_WAVE, IERR_MPI )
                END IF
!
! 4.d This is NAPTRK, perform all output
!
              IF ( IAPROC .EQ. NAPTRK ) THEN
!
! 4.e Sea point, prepare data
!
                  IF ( TSTSTR .EQ. 'SEA' ) THEN
!
                      WX     = UA(ISEA) * COS(UD(ISEA))
                      WY     = UA(ISEA) * SIN(UD(ISEA))
!
! ..... Local spectra
!
                      IF ( IAPROC .EQ. ISPROC ) THEN
                          DO IK=1, NK
                            DO ITH=1, NTH
                              SPEC(IK,ITH) =                          &
                               TPI*A(ITH,IK,JSEA)*SIG(IK)/CG(IK,ISEA)
                              END DO
                            END DO
!
! ..... Non-local spectra
!
                        ELSE
                          CALL MPI_RECV (ASPTRK, NSPEC, MPI_REAL,&
                                     IFROM, IT, MPI_COMM_WAVE,   &
                                     STATUS, IERR_MPI )
!
                          DO IK=1, NK
                            DO ITH=1, NTH
                              SPEC(IK,ITH) =                          &
                               TPI*ASPTRK(ITH,IK)*SIG(IK)/CG(IK,ISEA)
                              END DO
                            END DO
                        END IF
!
! 4.e Sea point, write general data + spectrum
!
                      WRITE (NDSTO,ERR=811,IOSTAT=IERR)               &
                         TIME, X, Y, TSTSTR
                      WRITE (NDSTO,ERR=811,IOSTAT=IERR)               &
                         DW(ISEA), CX(ISEA), CY(ISEA), WX, WY,        &
                         UST(ISEA), AS(ISEA), SPEC
!
! 4.f Non-sea point, write
!
                    ELSE
                      WRITE (NDSTO,ERR=811,IOSTAT=IERR)               &
                         TIME, X, Y, TSTSTR
!
! ..... Sea and non-sea points processed
!
                    END IF
!
! ..... End of action at NAPTRK
!
                END IF
!
! ..... Close IF for mask flag (top section 4)
!
            END IF
!
! ..... End of loop over map
!
          END DO
        END DO
!
      DEALLOCATE ( STATUS )
!
      GOTO 888
!
!     Error Escape Locations
!
  800 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1000) IERR
      GOTO 880
!
  801 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1001) IERR
      GOTO 880
!
  802 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1002)
      GOTO 880
!
  803 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1003) IDSTRI, IDTST
      GOTO 880
!
  810 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1010) IERR
      GOTO 880
!
  811 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1011) IERR
!
!     Disabeling output
!
  880 CONTINUE
      ATOLAST(:,3) = TIME
!
  888 CONTINUE
!
      RETURN
!
! Formats
!
 1000 FORMAT (/' *** WAVEWATCH III WARNING IN W3IOTR : '/             &
               '     INPUT FILE WITH TRACK DATA NOT FOUND ',          &
               '(FILE track_i.ww3, IOSTAT =',I6,')'/                  &
               '     TRACK OUTPUT DISABLED '/)
 1001 FORMAT (/' *** WAVEWATCH III WARNING IN W3IOTR : '/             &
               '     ERROR IN READING FILE track_i.ww3, IOSTAT =',I6/ &
               '     (ADITIONAL) TRACK OUTPUT DISABLED '/)
 1002 FORMAT (/' *** WAVEWATCH III WARNING IN W3IOTR : '/             &
               '     PREMATURE END OF FILE track_i.ww3'/              &
               '     TRACK OUTPUT DISABLED '/)
 1003 FORMAT (/' *** WAVEWATCH III WARNING IN W3IOTR : '/             &
               '     UNEXPECTED CONTENTS OF OF FILE track_i.ww3'/     &
               '       EXPECTED : ',A/                                &
               '       FOUND    : ',A/                                &
               '     TRACK OUTPUT DISABLED '/)
 1004 FORMAT (/' *** WAVEWATCH III WARNING IN W3IOTR : '/             &
               '     END OF FILE track_i.ww3'/)
!
 1010 FORMAT (/' *** WAVEWATCH III WARNING IN W3IOTR : '/             &
               '     ERROR IN OPENING OUTPUT FILE ',                  &
               '(FILE track_o.ww3, IOSTAT =',I6,')'/                  &
               '     TRACK OUTPUT DISABLED '/)
 1011 FORMAT (/' *** WAVEWATCH III WARNING IN W3IOTR : '/             &
               '     ERROR IN WRITING TO FILE track_o.ww3, IOSTAT =',I6/ &
               '     (ADITIONAL) TRACK OUTPUT DISABLED '/)
!
!/
!/ End of W3IOTR ----------------------------------------------------- /
!/
      END SUBROUTINE W3IOTR
!/
!/ End of module W3IOTRMD -------------------------------------------- /
!/
      END MODULE W3IOTRMD
