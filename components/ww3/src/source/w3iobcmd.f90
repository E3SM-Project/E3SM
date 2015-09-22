!/ ------------------------------------------------------------------- /
      MODULE W3IOBCMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    See subroutine for update log.
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Processing of boundary data output.
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      VERBPT    C*10  Private  Restart file version number.
!      IDSTR     C*32  Private  Restart file ID string.
!     ----------------------------------------------------------------
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3IOBC    Subr. Public   Boundary data IO.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SETO, W3SETG, W3SETW, W3SETA, W3DMO5
!                Subr. W3xDATMD Manage data structures.
!      W3CSPC    Subr. W3CSPCMD Spectral grid conversion.
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!      EXTCDE    Subr. W3SERVMD Abort program with exit code.
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!     None.
!
!  6. Switches :
!
!     See subroutine W3IOBC.
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
      PUBLIC
!/
!/ Private parameter statements (ID strings)
!/
      CHARACTER(LEN=10), PARAMETER, PRIVATE :: VERBPT = 'III  1.03 '
      CHARACTER(LEN=32), PARAMETER, PRIVATE ::                        &
                          IDSTR = 'WAVEWATCH III BOUNDARY DATA FILE'
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3IOBC ( INXOUT, NDSB, TIME1, TIME2, IOTST, IMOD )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    12-Jan-1999 : Distributed FORTRAN 77 version.     ( version 1.18 )
!/    20-May-1999 : Remove read bug for IPBP and RDBP   ( see web page )
!/    30-Dec-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/                  Major changes to logistics.
!/    13-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    19-Sep-2005 : Allow for change of spec. res.      ( version 3.08 )
!/                  (on read only).
!/    30-Sep-2005 : Add 'DUMP' option.                  ( version 3.08 )
!/    27-Jun-2006 : Adding file name preamble.          ( version 3.09 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!  1. Purpose :
!
!     Write/read boundary conditions file(s).
!
!  2. Method :
!
!     The file(s) are opened within the routine, the names are
!     pre-defined as nest.FILEXT for the input file and nest1.FILEXT
!     through nest9.FILEXT for up to 9 output files.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       INXOUT  C*(*)  I   Test string for read/write, valid are:
!                          'READ', 'WRITE' or 'DUMP'.
!       NDSB    Int.   I   Data set unit number.
!       TIME1   I.A.  I/O  Present time.                          (w)
!                          Time of first field.                   (r)
!       TIME2   I.A.   O   Time of second field.                  (r)
!       IOTST   Int.   O   Test indictor for reading.
!                           1 : File not found.
!                           0 : Fields read.
!                          -1 : Past end of file.
!       IMOD    Int.   I   Optional grid number, defaults to 1.
!     ----------------------------------------------------------------
!                                            (w) used for write only
!                                            (r) used for write only
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
!       Tests on INXOUT, file status and data present in file.
!
!  7. Remarks :
!
!     - Array dimensions are tested in W3IOGR.
!     - Spectra are stored as frequency (sigma) spectra to guarantee
!       conservation under grid transformation.
!     - At the moment it is mplicitly assumed that the number of
!       spectral components is larger that the number of spectra
!       per time step per file.
!     - Dump option used in multi-grid model.
!
!  8. Structure :
!
!       See source code.
!
!  9. Switches :
!
!     !/SHRD  Switch for shared / distributed memory architecture.
!     !/DIST  Id.
!
!     !/S     Enable subroutine tracing.
!     !/T     General test output.
!     !/T0    Point info test output.
!     !/T1    Wave heights at input/output points.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3CONSTANTS
!
      USE W3GDATMD, ONLY: W3SETG
      USE W3WDATMD, ONLY: W3SETW
      USE W3ADATMD, ONLY: W3SETA
      USE W3ODATMD, ONLY: W3SETO, W3DMO5
      USE W3CSPCMD, ONLY: W3CSPC
!
      USE W3GDATMD, ONLY: NK, NTH, NSPEC, NSEA, NSEAL, NX, NY,        &
                          X0, Y0, SX, SY, MAPSTA, MAPFS, MAPSF,       &
                          XFR, FR1, SIG2, TH, DTH, FILEXT, FACHFE
      USE W3WDATMD, ONLY: VA
      USE W3ADATMD, ONLY: CG
      USE W3ODATMD, ONLY: NDSE, NDST, IAPROC, NAPROC, NAPERR, NAPBPT, &
                          NBI, NBI2, NFBPO, NBO, NBO2, NDSL,          &
                          NKI, NTHI, XFRI, FR1I, TH1I,                &
                          IPBPI, ISBPI, XBPI, YBPI, RDBPI,            &
                          IPBPO, ISBPO, XBPO, YBPO, RDBPO,            &
                          ABPI0, ABPIN, ABPOS, FLBPI, FILER, FILEW,   &
                          FILED, SPCONV, FNMPRE
!
      USE W3SERVMD, ONLY: EXTCDE
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)           :: NDSB
      INTEGER, INTENT(INOUT)        :: TIME1(2)
      INTEGER, INTENT(OUT)          :: TIME2(2), IOTST
      INTEGER, INTENT(IN), OPTIONAL :: IMOD
      CHARACTER, INTENT(IN)         :: INXOUT*(*)
!/
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IFILE, IERR, I, J, IX, IY,&
                                 ISEA, IP, ISP, NPTS, ISOUT, IS, IGRD
      REAL, ALLOCATABLE       :: TMPSPC(:,:)
      LOGICAL                 :: FLOK
      CHARACTER(LEN=16)       :: FILEN
      CHARACTER(LEN=10)       :: VERTST
      CHARACTER(LEN=32)       :: IDTST
!/
!/ ------------------------------------------------------------------- /
!/
!
      IOTST  = 0
!
! test parameter list input ------------------------------------------ *
!
      IF ( PRESENT(IMOD) ) THEN
          IGRD   = IMOD
        ELSE
          IGRD   = 1
        END IF
!
      CALL W3SETO ( IGRD, NDSE, NDST )
      CALL W3SETG ( IGRD, NDSE, NDST )
      CALL W3SETW ( IGRD, NDSE, NDST )
      CALL W3SETA ( IGRD, NDSE, NDST )
!
      IF (INXOUT.NE.'READ' .AND. INXOUT.NE.'WRITE' .AND.              &
          INXOUT.NE.'DUMP' ) THEN
          IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,900) INXOUT
          CALL EXTCDE ( 1 )
        END IF
!
! open file ---------------------------------------------------------- *
!
      I      = LEN_TRIM(FILEXT)
      J      = LEN_TRIM(FNMPRE)
!
      IF ( INXOUT.EQ.'READ'  .AND. FILER ) THEN
          WRITE (FILEN,'(A5,A)') 'nest.', FILEXT(:I)
          OPEN (NDSB,FILE=FNMPRE(:J)//FILEN(:5+I),FORM='UNFORMATTED', &
                ERR=801,IOSTAT=IERR,STATUS='OLD')
        END IF
!
      IF ( INXOUT.EQ.'WRITE' .AND. FILEW ) THEN
          DO IFILE=1, NFBPO
            NDSL(IFILE) = NDSB + IFILE - 1
            WRITE (FILEN,'(A4,I1,A1,A)') 'nest', IFILE, '.',          &
                                         FILEXT(:I)
            OPEN (NDSL(IFILE),FILE=FNMPRE(:J)//FILEN(:6+I),           &
                  FORM='UNFORMATTED',ERR=800,IOSTAT=IERR)
            END DO
        END IF
!
      IF ( INXOUT.EQ.'DUMP'  .AND. FILED ) THEN
          WRITE (FILEN,'(A5,A)') 'nest.', FILEXT(:I)
          OPEN (NDSB,FILE=FNMPRE(:J)//FILEN(:5+I),FORM='UNFORMATTED', &
                ERR=800,IOSTAT=IERR)
        END IF
!
! test info ---------------------------------------------------------- *
! ( new files only )
! ... writing
!
      IF ( INXOUT.EQ.'WRITE' .AND. FILEW ) THEN
          IF ( IAPROC .EQ. NAPBPT ) THEN
              DO IFILE=1, NFBPO
                WRITE (NDSL(IFILE))                                   &
                      IDSTR, VERBPT, NK, NTH, XFR, FR1, TH(1),        &
                      NBO(IFILE)-NBO(IFILE-1)
!
                WRITE (NDSL(IFILE))                                   &
                     (XBPO(I),I=NBO(IFILE-1)+1,NBO(IFILE)),           &
                     (YBPO(I),I=NBO(IFILE-1)+1,NBO(IFILE)),           &
                     ((IPBPO(I,J),I=NBO(IFILE-1)+1,NBO(IFILE)),J=1,4),&
                     ((RDBPO(I,J),I=NBO(IFILE-1)+1,NBO(IFILE)),J=1,4)
!
                END DO
            END IF
        END IF
!
! ... dumping
!
      IF ( INXOUT.EQ.'DUMP' .AND. FILED ) THEN
          IF ( IAPROC .EQ. NAPBPT ) THEN
              WRITE (NDSB) IDSTR, VERBPT, NK, NTH, XFR, FR1, TH(1), NBI
!
              WRITE (NDSB) (XBPI(I),I=1,NBI), (YBPI(I),I=1,NBI),      &
                           ((IPBPI(I,J),I=1,NBI),J=1,4),              &
                           ((RDBPI(I,J),I=1,NBI),J=1,4)
!
            END IF
        END IF
!
! ... reading
!
      IF ( INXOUT.EQ.'READ' .AND. FILER ) THEN
!
          READ (NDSB,ERR=803,IOSTAT=IERR)                             &
                IDTST, VERTST, NKI, NTHI, XFRI, FR1I, TH1I, NBI
!
          IF ( IDTST .NE. IDSTR ) THEN
              IF ( IAPROC .EQ. NAPERR )                               &
                   WRITE (NDSE,901) IDTST, IDSTR
              CALL EXTCDE ( 10 )
            END IF
          IF ( VERTST .NE. VERBPT ) THEN
              IF ( IAPROC .EQ. NAPERR )                               &
                   WRITE (NDSE,902) VERTST, VERBPT
              CALL EXTCDE ( 11 )
            END IF
          SPCONV = NKI.NE.NK .OR. NTHI.NE.NTH .OR.                    &
                   ABS(XFRI/XFR-1.).GT.0.01 .OR.                      &
                   ABS(FR1I/FR1-1.).GT.0.01 .OR.                      &
                   ABS(TH1I-TH(1)).GT.0.01*DTH
!
          CALL W3DMO5 ( IGRD, NDSE, NDST, 1 )
!
          READ (NDSB,ERR=803,IOSTAT=IERR)                             &
              (XBPI(I),I=1,NBI), (YBPI(I),I=1,NBI),                   &
              ((IPBPI(I,J),I=1,NBI),J=1,4),                           &
              ((RDBPI(I,J),I=1,NBI),J=1,4)
!
          FLOK   = .TRUE.
          DO I=1, NBI
            IX     = 1 + NINT((XBPI(I)-X0)/SX)
            IY     = 1 + NINT((YBPI(I)-Y0)/SY)
            IF ( ABS(MAPSTA(IY,IX)) .NE. 2 ) THEN
                IF ( IAPROC .EQ. NAPERR )                             &
                    WRITE (NDSE,910) IX, IY, ABS(MAPSTA(IY,IX))
                FLOK   = .FALSE.
              END IF
            ISBPI(I) = MAPFS(IY,IX)
            END DO
!
          IF ( .NOT.FLOK ) CALL EXTCDE ( 20 )
!
          DO ISEA=1, NSEA
            IX     = MAPSF(ISEA,1)
            IY     = MAPSF(ISEA,2)
            IF ( ABS(MAPSTA(IY,IX)) .EQ. 2 ) THEN
                FLOK   = .FALSE.
                DO I=1, NBI
                  IF ( ISEA .EQ. ISBPI(I) ) FLOK = .TRUE.
                  END DO
                IF ( .NOT.FLOK .AND. IAPROC.EQ.NAPERR )               &
                    WRITE (NDSE,911) IX, IY
              END IF
            END DO
!
!     Read first time and allocate ABPI0/N
!
          READ (NDSB,END=810,ERR=810) TIME2, NBI2
          BACKSPACE (NDSB)
          CALL W3DMO5 ( IGRD, NDSE, NDST, 3 )
!
        END IF
!
! Save previous spectra on read -------------------------------------- *
!
      IF ( INXOUT.EQ.'READ' .AND. .NOT.FILER ) THEN
          TIME1  = TIME2
          ABPI0(:,1:NBI2) = ABPIN(:,1:NBI2)
        END IF
!
! TIME --------------------------------------------------------------- *
!
      IF ( INXOUT .EQ. 'WRITE'  ) THEN
          DO IFILE=1, NFBPO
            NPTS   = NBO2(IFILE) - NBO2(IFILE-1)
            WRITE (NDSL(IFILE)) TIME1, NPTS
            END DO
        END IF
!
      IF ( INXOUT .EQ. 'DUMP'  ) THEN
          WRITE (NDSB) TIME1, NBI2
        END IF
!
      IF ( INXOUT .EQ. 'READ'  ) THEN
          READ (NDSB,ERR=810,END=810) TIME2, NBI2
        END IF
!
! Spectra ------------------------------------------------------------ *
!
      IF ( INXOUT .EQ. 'WRITE' ) THEN
!
          DO IFILE=1, NFBPO
            DO ISOUT=NBO2(IFILE-1)+1, NBO2(IFILE)
!
              ISEA   = ISBPO(ISOUT)
!
! ... Shared memory version data gather
!
! ... Distributed memory version data gather
!   ( Array pre-filled in W3WAVE )
!
              DO IS=1, NSPEC
                ABPOS(IS,ISOUT) = ABPOS(IS,ISOUT) * SIG2(IS) /  &
                                    CG(1+(IS-1)/NTH,ISEA)
                END DO
!
              WRITE (NDSL(IFILE)) (ABPOS(IS,ISOUT),IS=1,NSPEC)
!
              END DO
            END DO
!
        END IF
!
      IF ( INXOUT .EQ. 'DUMP' ) THEN
          DO I=1, NBI2
            WRITE (NDSB) ABPIN(:,I)
            END DO
        END IF
!
      IF ( INXOUT .EQ. 'READ' ) THEN
!
          IF ( .NOT. SPCONV ) THEN
              DO IP=1, NBI2
                READ (NDSB,ERR=803,IOSTAT=IERR) ABPIN(:,IP)
                END DO
            ELSE
              ALLOCATE ( TMPSPC(NKI*NTHI,NBI2) )
              DO IP=1, NBI2
                READ (NDSB,ERR=803,IOSTAT=IERR) TMPSPC(:,IP)
                END DO
              CALL W3CSPC ( TMPSPC     , NKI, NTHI, XFRI, FR1I, TH1I, &
                            ABPIN(:,1:), NK , NTH , XFR , FR1, TH(1), &
                            NBI2, NDST, NDSE, FACHFE )
              DEALLOCATE ( TMPSPC )
            END IF
!
        END IF
!
! Set first spectra on first read ------------------------------------ *
!
      IF ( INXOUT.EQ.'READ' .AND. FILER ) THEN
          TIME1 = TIME2
          DO IP=1, NBI2
            ABPI0(:,IP) = ABPIN(:,IP)
            END DO
          ABPI0(:,0) = 0.
          ABPIN(:,0) = 0.
        END IF
!
! Reset flags -------------------------------------------------------- *
!
      IF ( INXOUT .EQ. 'WRITE' ) FILEW  = .FALSE.
      IF ( INXOUT .EQ. 'DUMP'  ) FILED  = .FALSE.
      IF ( INXOUT .EQ. 'READ'  ) FILER  = .FALSE.
!
      RETURN
!
! Escape locations IO errors
!
  800 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1000) FILEN, IERR
      CALL EXTCDE ( 40 )
!
  801 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1001)
      IOTST  = 1
      FLBPI  = .FALSE.
      RETURN
!
  802 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1002)
      CALL EXTCDE ( 41 )
!
  803 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1003) IERR
      CALL EXTCDE ( 42 )
!
  810 CONTINUE
      IF ( FILER ) THEN
          IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1010)
          CALL EXTCDE ( 43 )
        END IF
!
      TIME1(1) = TIME2(1)
      TIME1(2) = TIME2(2)
      DO 812, IP=0, NBI2
        DO 811, ISP=1, NSPEC
          ABPI0(ISP,IP) = ABPIN(ISP,IP)
  811     CONTINUE
  812   CONTINUE
!
      IOTST  = -1
      FLBPI  = .FALSE.
      RETURN
!
! Formats
!
  900 FORMAT (/' *** WAVEWATCH III ERROR IN W3IOBC :'/                &
               '     ILEGAL INXOUT VALUE: ',A/)
  901 FORMAT (/' *** WAVEWATCH III ERROR IN W3IOBC :'/                &
               '     ILEGAL IDSTR, READ : ',A/                        &
               '                  CHECK : ',A/)
  902 FORMAT (/' *** WAVEWATCH III ERROR IN W3IOBC :'/                &
               '     ILEGAL VEROGR, READ : ',A/                       &
               '                   CHECK : ',A/)
!
  910 FORMAT (/' *** WAVEWATCH III ERROR IN W3IOBC :'/                &
               '     POINT',2I4,' NOT ACTIVE SEA POINT (',I1,')')
  911 FORMAT ( ' *** WAVEWATCH III WARNING : POINT',2I4,              &
                   ' WILL NOT BE UPDATED')
!
 1000 FORMAT (/' *** WAVEWATCH III ERROR IN W3IOBC : '/               &
               '     ERROR IN OPENING FILE ',A/                       &
               '     IOSTAT =',I5/)
 1001 FORMAT (/' *** WAVEWATCH III WARNING IN W3IOBC : '/             &
               '     INPUT FILE WITH BOUNDARY CONDITIONS NOT FOUND'/  &
               '     BOUNDARY CONDITIONS WILL NOT BE UPDATED '/)
 1002 FORMAT (/' *** WAVEWATCH III ERROR IN W3IOBC : '/               &
               '     PREMATURE END OF FILE'/)
 1003 FORMAT (/' *** WAVEWATCH III ERROR IN W3IOBC : '/               &
               '     ERROR IN READING FROM FILE'/                     &
               '     IOSTAT =',I5/)
!
 1010 FORMAT (/' *** WAVEWATCH III ERROR IN W3IOBC : '/               &
               '     NO DATA IN INPUT FILE'/)
!
!/
!/ End of W3IOBC ----------------------------------------------------- /
!/
      END SUBROUTINE W3IOBC
!/
!/ End of module W3IOBCMD -------------------------------------------- /
!/
      END MODULE W3IOBCMD
