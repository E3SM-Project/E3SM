!/ ------------------------------------------------------------------- /
      MODULE W3IORSMD
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
!  1. Purpose :
!
!     Read/write restart files.
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      VERINI    C*10  Private  Restart file version number.
!      IDSTR     C*26  Private  Restart file UD string.
!     ----------------------------------------------------------------
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3IORS    Subr. Public   Read/write restart files.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SETO, W3SETG, W3SETW, W3DIMW
!                Subr. W3xDATMD Manage data structures.
!      STRACE    Subr. W3SERVMD Subroutine tracing.            (!/S)
!      EXTCDE    Subr. W3SERVMD Abort program with exit code.
!      MPI_STARTALL, MPI_WAITALL                              (!/MPI)
!                Subr.          MPI persistent communication routines
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!  6. Switches :
!
!     See also routine.
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
      PUBLIC
!/
!/ Private parameter statements (ID strings)
!/
      CHARACTER(LEN=10), PARAMETER, PRIVATE :: VERINI = 'III  3.01 '
      CHARACTER(LEN=26), PARAMETER, PRIVATE ::                        &
                               IDSTR = 'WAVEWATCH III RESTART FILE'
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3IORS ( INXOUT, NDSR, DUMFPI, INTYPE, IMOD )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    12-Jan-1999 : Final FORTRAN 77                    ( version 1.18 )
!/    27-Dec-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    30-Apr-2002 : Add ice for transparencies.         ( version 2.20 )
!/    13-Nov-2002 : Add stress as vector.               ( version 3.00 )
!/    19-Aug-2003 : Output server options added.        ( version 3.04 )
!/    09-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    24-Jun-2005 : Adding MAPST2.                      ( version 3.07 )
!/    27-Jun-2006 : Adding file name preamble.          ( version 3.09 )
!/    05-Jul-2006 : Consolidate stress arrays.          ( version 3.09 )
!/    08-May-2007 : Starting from calm as an option.    ( version 3.11 )
!/    17-May-2007 : Adding NTPROC/NAPROC separation.    ( version 3.11 )
!/    22-Jun-2007 : Dedicated output processes.         ( version 3.11 )
!/    15-Apr-2008 : Clean up for distribution.          ( version 3.14 )
!/    21-Apr-2008 : Remove PGI bug internal files.      ( version 3.14 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Reads/writes restart files.
!
!  2. Method :
!
!     The file is opened within the routine, the name is pre-defined
!     and the unit number is given in the parameter list. The restart
!     file is written using UNFORMATTED write statements. The routine
!     generates new names when called more than once. File names are :
!
!                                 restart.FILEXT
!                                 restart1.FILEXT
!                                 restart2.FILEXT etc.
!
!     The file to be read thus always is unnumbered, whereas all
!     written files are automatically numbered.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       INXOUT  C*(*)  I   Test string for read/write, valid are:
!                          'READ' Reading of a restart file.
!                          'HOT'  Writing a full restart from the model.
!                          'COLD' Writing a cold start file.
!                          'WIND' Initialize fields using first wind
!                                 field.
!                          'CALM' Starting from calm conditions.
!       NDSR    Int.  I/O  File unit number.
!       DUMFPI  Real   I   Dummy values for FPIS for cold start.
!       INTYPE  Int.   O   Type of input field,
!                           0 : cold start,
!                           1 : cold start with fetch-limited spectra,
!                           2 : full restart,
!                           3 : for writing file.
!                           4 : starting from calm.
!       IMOD    Int.   I   Optional grid number, defaults to 1.
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
!      W3INIT    Subr. W3INITMD Wave model initialization routine.
!      W3WAVE    Subr. W3WAVEMD Actual wave model routine.
!      WW3_STRT  Prog.   N/A    Initial conditions program.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       Tests on INXOUT, file status and on array dimensions.
!
!  7. Remarks :
!
!     - MAPSTA is dumped as it contains information on inactive points.
!       Note that the original MAPSTA is dumped in the model def. file
!       for use in the initial conditions (and output) programs.
!     - Note that MAPSTA and MAPST2 data is combinded in the file.
!     - The depth is recalculated in a write to avoid floting point
!       errors in W3STRT.
!     - Fields and field info read by all, written by las processor
!       only.
!     - The MPP version of the model will perform a gather here to
!       maximize hiding of communication with IO.
!
!  8. Structure :
!
!     +---------------------------------------------------------------+
!     | initialisations                                               |
!     | test INXOUT                                                   |
!     | open file                                                     |
!     +---------------------------------------------------------------|
!     |                             WRITE ?                           |
!     | Y                                                           N |
!     |-------------------------------|-------------------------------|
!     | Write identifiers and         | Write identifiers and         |
!     |   dimensions.                 |   dimensions.                 |
!     |                               | Check ident. and dimensions.  |
!     +-------------------------------+-------------------------------|
!     |                       Full restart ?                          |
!     | Y                                                           N |
!     |-------------------------------|-------------------------------|
!     | read/write/test time          |                               |
!     +-------------------------------+-------------------------------|
!     |                             WRITE ?                           |
!     | Y                                                           N |
!     |-------------------------------|-------------------------------|
!     |          TYPE = 'WIND' ?      |          TYPE = 'WIND' ?      |
!     | Y                           N | Y                           N |
!     |---------------|---------------|---------------|---------------|
!     | close file    | write spectra | gen. fetch-l. | read spectra  |
!     | RETURN        |               |   spectra.    |               |
!     |---------------+---------------+---------------+---------------|
!     |                             WRITE ?                           |
!     | Y                                                           N |
!     |-------------------------------|-------------------------------|
!     |          TYPE = 'FULL' ?      |          TYPE = 'FULL' ?      |
!     | Y                           N | Y                           N |
!     |---------------|---------------|---------------|---------------|
!     | write level & | ( prep. level | read level &  | initalize l.& |
!     |   (ice) map & |   for test    |   (ice) map.& |   times       |
!     |   times       |   output )    |   times       | ( no ice )    |
!     +---------------+---------------+---------------+-------------- +
!
!  9. Switches :
!
!     !/SEED  Linear input / seeding option.
!     !/LNx
!
!     !/SHRD  Switch for shared / distributed memory architecture.
!     !/DIST  Id.
!     !/MPI   Id.
!
!     !/LRBn  Word length in bites.
!
!     !/S     Enable subroutine tracing.
!     !/T     Enable test output
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3GDATMD, ONLY: W3SETG
      USE W3ODATMD, ONLY: W3SETO
!/
      USE W3GDATMD, ONLY: NX, NY, NSEA, NSEAL, NSPEC, MAPSTA, MAPST2, &
                          GNAME, FILEXT
      USE W3WDATMD
      USE W3ODATMD, ONLY: NDSE, NDST, IAPROC, NAPROC, NAPERR, NAPRST, &
                          IFILE => IFILE4, FNMPRE, NTPROC, IOSTYP
      USE W3ODATMD, ONLY: NRQRS, NBLKRS, RSBLKS, IRQRS, IRQRSS, VAAUX
!/
      USE W3SERVMD, ONLY: EXTCDE
!
      IMPLICIT NONE
!
      INCLUDE "mpif.h"
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)           :: NDSR
      INTEGER, INTENT(IN), OPTIONAL :: IMOD
      INTEGER, INTENT(OUT)          :: INTYPE
      REAL, INTENT(INOUT)           :: DUMFPI
      CHARACTER, INTENT(IN)         :: INXOUT*(*)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER, PARAMETER      :: LRB = 8
!
      INTEGER                 :: IGRD, I, J, LRECL, NSIZE, IYS, IERR, &
                                 NSEAT, MSPEC, TTIME(2), ISEA, JSEA,  &
                                 NREC, NPART, IPART, IX, IY, IYL, IP
      INTEGER, ALLOCATABLE    :: MAPTMP(:,:)
      INTEGER                 :: IERR_MPI, IH, IB, ISEA0, ISEAN, &
                                 NRQ
      INTEGER, ALLOCATABLE    :: STAT1(:,:), STAT2(:,:)
      LOGICAL                 :: WRITE, IOSFLG
      CHARACTER(LEN=4)        :: TYPE
      CHARACTER(LEN=10)       :: VERTST
      CHARACTER(LEN=19)       :: FNAME
      CHARACTER(LEN=26)       :: IDTST
      CHARACTER(LEN=30)       :: TNAME
!/
!/ ------------------------------------------------------------------- /
!/
!
      IOSFLG = IOSTYP .GT. 0
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
!
      IF (INXOUT.NE.'READ' .AND. INXOUT.NE.'HOT'  .AND.               &
          INXOUT.NE.'COLD' .AND. INXOUT.NE.'WIND' .AND.               &
          INXOUT.NE.'CALM' ) THEN
          IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,900) INXOUT
          CALL EXTCDE ( 1 )
        END IF
!
      WRITE = INXOUT .NE. 'READ'
      IF ( INXOUT .EQ. 'HOT' ) THEN
          TYPE   = 'FULL'
        ELSE
          TYPE   = INXOUT
        END IF
!
! initializations ---------------------------------------------------- *
!
      IF ( .NOT.DINIT ) THEN
          IF ( IAPROC .LE. NAPROC ) THEN
              CALL W3DIMW ( IMOD, NDSE, NDST )
            ELSE
              CALL W3DIMW ( IMOD, NDSE, NDST, .FALSE. )
            END IF
        END IF
!
      IF ( IAPROC .LE. NAPROC ) VA(:,0) = 0.
!
      LRECL  = MAX ( LRB*NSPEC , LRB*NX ,                             &
                     LRB*(6+(25/LRB)+(9/LRB)+(29/LRB)+(3/LRB)) )
      NSIZE  = LRECL / LRB
      IYS    = NSIZE / NX
!
! open file ---------------------------------------------------------- *
!
      I      = LEN_TRIM(FILEXT)
      J      = LEN_TRIM(FNMPRE)
!
      IF ( IFILE.EQ.0 ) THEN
          FNAME  = 'restart.'//FILEXT(:I)
        ELSE
          FNAME  = 'restartN.'//FILEXT(:I)
          IF ( WRITE .AND. IAPROC.EQ.NAPRST )                         &
               WRITE (FNAME(8:8),'(I1)') IFILE
        END IF
      IFILE  = IFILE + 1
!
      IF ( WRITE ) THEN
          IF ( .NOT.IOSFLG .OR. IAPROC.EQ.NAPRST )                    &
          OPEN (NDSR,FILE=FNMPRE(:J)//FNAME,FORM='UNFORMATTED',       &
                ACCESS='DIRECT',RECL=LRECL,ERR=800,IOSTAT=IERR)
        ELSE
          OPEN (NDSR,FILE=FNMPRE(:J)//FNAME,FORM='UNFORMATTED',       &
                ACCESS='DIRECT',RECL=LRECL,ERR=800,IOSTAT=IERR,       &
                STATUS='OLD')
        END IF
!
! test info ---------------------------------------------------------- *
!
      IF ( WRITE ) THEN
!
          IF ( IAPROC .EQ. NAPRST ) WRITE (NDSR,REC=1)                &
            IDSTR, VERINI, GNAME, TYPE, NSEA, NSPEC
          INTYPE = 3
!
        ELSE
          READ (NDSR,REC=1,ERR=802,IOSTAT=IERR)                       &
            IDTST, VERTST, TNAME, TYPE, NSEAT, MSPEC
!
          IF ( IDTST .NE. IDSTR ) THEN
              IF ( IAPROC .EQ. NAPERR )                               &
                  WRITE (NDSE,901) IDTST, IDSTR
              CALL EXTCDE ( 10 )
            END IF
          IF ( VERTST .NE. VERINI ) THEN
              IF ( IAPROC .EQ. NAPERR )                               &
                  WRITE (NDSE,902) VERTST, VERINI
              CALL EXTCDE ( 11 )
            END IF
          IF ( TNAME .NE. GNAME ) THEN
              IF ( IAPROC .EQ. NAPERR )                               &
                  WRITE (NDSE,903) TNAME, GNAME
            END IF
          IF (TYPE.NE.'FULL' .AND. TYPE.NE.'COLD' .AND.               &
              TYPE.NE.'WIND' .AND. TYPE.NE.'CALM' ) THEN
              IF ( IAPROC .EQ. NAPERR )                               &
                  WRITE (NDSE,904) TYPE
              CALL EXTCDE ( 12 )
            END IF
          IF (NSEAT.NE.NSEA .OR. NSPEC.NE.MSPEC) THEN
              IF ( IAPROC .EQ. NAPERR )                               &
                  WRITE (NDSE,905) MSPEC, NSEAT, NSPEC, NSEA
              CALL EXTCDE ( 13 )
            END IF
          IF (TYPE.EQ.'FULL') THEN
              INTYPE = 2
            ELSE IF (TYPE.EQ.'WIND') THEN
              INTYPE = 1
            ELSE IF (TYPE.EQ.'CALM') THEN
              INTYPE = 4
            ELSE
              INTYPE = 0
            END IF
!
        END IF
!
  100 CONTINUE
!
! TIME if required --------------------------------------------------- *
!
      IF (TYPE.EQ.'FULL') THEN
          IF ( WRITE ) THEN
              IF ( IAPROC .EQ. NAPRST ) WRITE (NDSR,REC=2) TIME
            ELSE
              READ (NDSR,REC=2,ERR=802,IOSTAT=IERR) TTIME
              IF (TIME(1).NE.TTIME(1) .OR. TIME(2).NE.TTIME(2)) THEN
                  IF ( IAPROC .EQ. NAPERR )                           &
                      WRITE (NDSE,906) TTIME, TIME
                  CALL EXTCDE ( 20 )
                END IF
            END IF
!
        END IF
!
! Spectra ------------------------------------------------------------ *
!          ( Bail out if write for TYPE.EQ.'WIND' )
!
      IF ( WRITE ) THEN
          IF ( TYPE.EQ.'WIND' .OR. TYPE.EQ.'CALM' ) THEN
              IF ( .NOT.IOSFLG .OR. IAPROC.EQ.NAPRST ) CLOSE ( NDSR )
              RETURN
            ELSE IF ( IAPROC.LE.NAPROC .OR. IAPROC.EQ. NAPRST ) THEN
!
! Original non-server version writing of spectra
!
              IF ( .NOT.IOSFLG .OR. (NAPROC.EQ.1.AND.NAPRST.EQ.1) ) THEN
                  DO JSEA=1, NSEAL
                    ISEA   = IAPROC + (JSEA-1)*NAPROC
                    NREC   = ISEA + 2
                    WRITE (NDSR,REC=NREC) (VA(I,JSEA),I=1,NSPEC)
                    END DO
!
! I/O server version writing of spectra ( !/MPI )
!
                ELSE
!
                  IF ( IAPROC .NE. NAPRST ) THEN
                      NRQ    = 1
                    ELSE IF ( NAPRST .LE. NAPROC ) THEN
                      NRQ    = NAPROC - 1
                    ELSE
                      NRQ    = NAPROC
                    END IF
!
                  ALLOCATE ( STAT1(MPI_STATUS_SIZE,NRQ) )
                  IF ( IAPROC .EQ. NAPRST ) CALL MPI_STARTALL    &
                                      ( NRQ, IRQRSS, IERR_MPI )
!
                  DO IB=1, NBLKRS
                    ISEA0  = 1 + (IB-1)*RSBLKS*NAPROC
                    ISEAN  = MIN ( NSEA , IB*RSBLKS*NAPROC )
!
                    IF ( IAPROC .EQ. NAPRST ) THEN
!
                        IH     = 1 + NRQ * (IB-1)
                        CALL MPI_WAITALL                         &
                           ( NRQ, IRQRSS(IH), STAT1, IERR_MPI )
                        IF ( IB .LT. NBLKRS ) THEN
                            IH     = 1 + NRQ * IB
                            CALL MPI_STARTALL                    &
                               ( NRQ, IRQRSS(IH), IERR_MPI )
                          END IF
!
                        DO ISEA=ISEA0, ISEAN
                          NREC   = ISEA + 2
                          JSEA   = 1 + (ISEA-1)/NAPROC
                          IP     = 1 + MOD(ISEA-1,NAPROC)
                          IF ( IP .EQ. NAPRST ) THEN
                              WRITE (NDSR,REC=NREC) VA(:,JSEA)
                            ELSE
                              JSEA   = JSEA - 2*((IB-1)/2)*RSBLKS
                              WRITE (NDSR,REC=NREC)              &
                                     VAAUX(:,JSEA,IP)
                            END IF
                          END DO
!
                      ELSE
!
                        CALL MPI_STARTALL                        &
                           ( 1, IRQRSS(IB), IERR_MPI )
                        CALL MPI_WAITALL                         &
                           ( 1, IRQRSS(IB), STAT1, IERR_MPI )
!
                      END IF
                    END DO
!
                  DEALLOCATE ( STAT1 )
!
                END IF
!
            END IF
        ELSE
!
! Reading spectra
!
          IF ( TYPE.EQ.'WIND' .OR. TYPE.EQ.'CALM' ) THEN
            ELSE
              DO JSEA=1, NSEAL
                ISEA   = IAPROC + (JSEA-1)*NAPROC
                NREC   = ISEA + 2
                READ (NDSR,REC=NREC,ERR=802,IOSTAT=IERR)              &
                         (VA(I,JSEA),I=1,NSPEC)
                END DO
            END IF
        END IF
!
! Water level etc. if required --------------------------------------- *
!     ( For cold start write test output and cold start initialize
!       water levels. Note that MAPSTA overwrites the one read from the
!       model definition file, so that it need not be initialized. )
!
      NREC   = NSEA + 3
      NPART  = 1 + (NSEA-1)/NSIZE
!
      IF ( WRITE ) THEN
!
          IF (TYPE.EQ.'FULL') THEN
!
              IF ( IAPROC .EQ. NAPRST ) THEN
!
                  ALLOCATE ( STAT2(MPI_STATUS_SIZE,NRQRS) )
                  CALL MPI_WAITALL                               &
                     ( NRQRS, IRQRS , STAT2, IERR_MPI )
                  DEALLOCATE ( STAT2 )
!
                  WRITE (NDSR,REC=NREC) TLEV, TICE
                  DO IPART=1,NPART
                    NREC  = NREC + 1
                    WRITE (NDSR,REC=NREC)                             &
                          (WLV(ISEA),ISEA=1+(IPART-1)*NSIZE,          &
                                          MIN(NSEA,IPART*NSIZE))
                    END DO
                  DO IPART=1,NPART
                    NREC  = NREC + 1
                    WRITE (NDSR,REC=NREC)                             &
                          (ICE(ISEA),ISEA=1+(IPART-1)*NSIZE,          &
                                          MIN(NSEA,IPART*NSIZE))
                    END DO
                  ALLOCATE ( MAPTMP(NY,NX) )
                  MAPTMP = MAPSTA + 8*MAPST2
                  DO IY=1, NY, IYS
                    NREC  = NREC + 1
                    WRITE (NDSR,REC=NREC)                             &
                          ((MAPTMP(IYL,IX),IX=1,NX),IYL=IY,           &
                                                    MIN(NY,IY+IYS-1))
                    END DO
                  DEALLOCATE ( MAPTMP )
                  DO IPART=1,NPART
                    NREC  = NREC + 1
                    WRITE (NDSR,REC=NREC)                             &
                          (UST(ISEA),ISEA=1+(IPART-1)*NSIZE,          &
                                          MIN(NSEA,IPART*NSIZE))
                    END DO
                  DO IPART=1,NPART
                    NREC  = NREC + 1
                    WRITE (NDSR,REC=NREC)                             &
                          (USTDIR(ISEA),ISEA=1+(IPART-1)*NSIZE,       &
                                          MIN(NSEA,IPART*NSIZE))
                    END DO
                  DO IPART=1,NPART
                    NREC  = NREC + 1
                    WRITE (NDSR,REC=NREC)                             &
                          (ASF(ISEA),ISEA=1+(IPART-1)*NSIZE,          &
                                          MIN(NSEA,IPART*NSIZE))
                    END DO
                  DO IPART=1,NPART
                    NREC  = NREC + 1
                    WRITE (NDSR,REC=NREC)                             &
                          (FPIS(ISEA),ISEA=1+(IPART-1)*NSIZE,         &
                                          MIN(NSEA,IPART*NSIZE))
                    END DO
                END IF
            END IF
        ELSE
          IF (TYPE.EQ.'FULL') THEN
              READ (NDSR,REC=NREC,ERR=802,IOSTAT=IERR) TLEV, TICE
              DO IPART=1,NPART
                NREC  = NREC + 1
                READ (NDSR,REC=NREC,ERR=802,IOSTAT=IERR)              &
                      (WLV(ISEA),ISEA=1+(IPART-1)*NSIZE,              &
                                      MIN(NSEA,IPART*NSIZE))
                END DO
              DO IPART=1,NPART
                NREC  = NREC + 1
                READ (NDSR,REC=NREC,ERR=802,IOSTAT=IERR)              &
                      (ICE(ISEA),ISEA=1+(IPART-1)*NSIZE,              &
                                      MIN(NSEA,IPART*NSIZE))
                END DO
              ALLOCATE ( MAPTMP(NY,NX) )
              DO IY=1, NY, IYS
                NREC  = NREC + 1
                READ (NDSR,REC=NREC,ERR=802,IOSTAT=IERR)              &
                      ((MAPTMP(IYL,IX),IX=1,NX),IYL=IY,               &
                                                MIN(NY,IY+IYS-1))
                END DO
              MAPSTA = MOD(MAPTMP+2,8) - 2
              MAPST2 = (MAPTMP-MAPSTA) / 8
              DEALLOCATE ( MAPTMP )
              DO IPART=1,NPART
                NREC  = NREC + 1
                READ (NDSR,REC=NREC,ERR=802,IOSTAT=IERR)              &
                      (UST(ISEA),ISEA=1+(IPART-1)*NSIZE,              &
                                      MIN(NSEA,IPART*NSIZE))
                END DO
              DO IPART=1,NPART
                NREC  = NREC + 1
                READ (NDSR,REC=NREC,ERR=802,IOSTAT=IERR)              &
                      (USTDIR(ISEA),ISEA=1+(IPART-1)*NSIZE,           &
                                      MIN(NSEA,IPART*NSIZE))
                END DO
              DO IPART=1,NPART
                NREC  = NREC + 1
                READ (NDSR,REC=NREC,ERR=802,IOSTAT=IERR)              &
                      (ASF(ISEA),ISEA=1+(IPART-1)*NSIZE,              &
                                      MIN(NSEA,IPART*NSIZE))
                END DO
              DO IPART=1,NPART
                NREC  = NREC + 1
                READ (NDSR,REC=NREC,ERR=802,IOSTAT=IERR)              &
                      (FPIS(ISEA),ISEA=1+(IPART-1)*NSIZE,             &
                                      MIN(NSEA,IPART*NSIZE))
                END DO
            ELSE
              TLEV(1) = -1
              TLEV(2) =  0
              TICE(1) = -1
              TICE(2) =  0
              WLV     =  0.
              ICE     =  0.
              UST     =  1.
              USTDIR  =  0.
              ASF     =  1.
              FPIS    =  DUMFPI
            END IF
        END IF
!
! Close file --------------------------------------------------------- *
!
      IF ( .NOT.IOSFLG .OR. IAPROC.EQ.NAPRST ) CLOSE ( NDSR )
!
      RETURN
!
! Escape locations read errors :
!
  800 CONTINUE
      TYPE   = 'CALM'
      INTYPE = 4
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,990) TYPE, IERR
      GOTO 100
!
  801 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,991)
      CALL EXTCDE ( 30 )
!
  802 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,992) IERR
      CALL EXTCDE ( 31 )
!
! Formats
!
  900 FORMAT (/' *** WAVEWATCH III ERROR IN W3IORS :'/                &
               '     ILLEGAL INXOUT VALUE: ',A/)
  901 FORMAT (/' *** WAVEWATCH III ERROR IN W3IORS :'/                &
               '     ILLEGAL IDSTR, READ : ',A/                       &
               '                   CHECK : ',A/)
  902 FORMAT (/' *** WAVEWATCH III ERROR IN W3IORS :'/                &
               '     ILLEGAL VERINI, READ : ',A/                      &
               '                    CHECK : ',A/)
  903 FORMAT (/' *** WAVEWATCH III WARNING IN W3IORS :'/              &
               '     ILLEGAL GNAME, READ : ',A/                       &
               '                   CHECK : ',A/)
  904 FORMAT (/' *** WAVEWATCH III ERROR IN W3IORS :'/                &
               '     ILLEGAL TYPE : ',A/)
  905 FORMAT (/' *** WAVEWATCH III ERROR IN W3IORS :'/                &
               '     CONFLICTING NSPEC, NSEA GRID : ',2I8/            &
               '                         EXPECTED : ',2I8/)
  906 FORMAT (/' *** WAVEWATCH III ERROR IN W3IORS :'/                &
               '     CONFLICTING TIMES: FILE : ',I10.8,I8.6/          &
               '                       MODEL : ',I10.8,I8.6/)
!
  990 FORMAT (/' *** WAVEWATCH III WARNING IN W3IORS : '/             &
               '     ERROR IN OPENING FILE, ',                        &
                    'INITIALIZE WITH ''',A,''' INSTEAD'/              &
               '     IOSTAT =',I5/)
  991 FORMAT (/' *** WAVEWATCH III ERROR IN W3IORS : '/               &
               '     PREMATURE END OF FILE'/)
  992 FORMAT (/' *** WAVEWATCH III ERROR IN W3IORS : '/               &
               '     ERROR IN READING FROM FILE'/                     &
               '     IOSTAT =',I5/)
!
!/
!/ End of W3IORS ----------------------------------------------------- /
!/
      END SUBROUTINE W3IORS
!/
!/ End of module W3IORSMD -------------------------------------------- /
!/
      END MODULE W3IORSMD
