!/ ------------------------------------------------------------------- /
      MODULE W3FLDSMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |            A. Chawla              |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    30-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    25-Jan-2002 : Data assimilation set up.           ( version 2.17 )
!/    26-Dec-2002 : Continuously moving grid.           ( version 3.02 )
!/    04-Sep-2003 : Bug fix W3FLHD.                     ( version 3.04 )
!/    27-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    05-Jul-2005 : Correct first level/ice.            ( version 3.07 )
!/    27-Jun-2006 : Adding file name preamble.          ( version 3.09 )
!/    09-Oct-2007 : Make file header optional.          ( version 3.13 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Gathers a set of routines to manage input fields of depth,
!     current, wind and ice concentration.
!
!  2. Variables and types :
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3FLDO    Subr. Public   Open data file.
!      W3FLDG    Subr. Public.  Read/write data file (fields).
!      W3FLDD    Subr. Public.  Read/write data file (data).
!      W3FLDP    Subr. Public.  Generic field interpolation.
!      W3FLDH    Subr. Public.  Process homogeneous fields.
!      W3FLDM    Subr. Public.  Process moving grid data.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.           ( !/S )
!      TICK21    Subr. W3TIMEMD Increment the clock.
!      DSEC21    R.F.  W3TIMEMD Calculate time differnces.
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!     - By design, these routines do not use the WAVEWATCH III data
!       structure. With this approach, they can be used in a straight-
!       forward way in other programs to generate WAVEWATCH III input
!       data sets directly from such programs.
!
!  6. Switches :
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
      PUBLIC
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3FLDO ( INXOUT, IDFLD, NDS, NDST, NDSE, NX, NY,     &
                          X0, Y0, SX, SY, IERR, FEXT, FPRE, FHDR )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |            A. Chawla              |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         09-Oct-2007 |
!/                  +-----------------------------------+
!/
!/    15-Jan-1999 : Final FORTRAN 77                    ( version 1.18 )
!/    30-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    24-Jan-2001 : Flat grid version (formats only)    ( version 2.06 )
!/    24-Jan-2002 : Assimilation data added.            ( version 2.17 )
!/    27-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    27-Jun-2006 : Adding file name preamble.          ( version 3.09 )
!/    09-Oct-2007 : Make file header optional.          ( version 3.13 )
!/
!  1. Purpose :
!
!     Open and prepare WAVEWATCH III field files as used by the
!     generic shell and the field preprocessor.
!
!  2. Method :
!
!     The file header contains a general WAVEWATCH III ID string,
!     a field ID string and the dimensions of the grid. If a file
!     is opened to be read, these parameters are all checked.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       INXOUT  C*(*) I  Test string for read/write, valid are:
!                        'READ' and 'WRITE'.
!       IDFLD   C*3  I/O ID string for field type, valid are:
!                        'LEV', 'CUR', 'WND', 'WNS', 'ICE' and 'DTn'.
!       NDS     Int.  I  Dataset number for fields file.
!       NDST    Int.  I  Dataset number for test output.
!       NDSE    Int.  I  Dataset number for error output.
!                        (No output if NDSE < 0).
!       NX, NY  Int.  I  Discrete grid dimensions.                 \
!       X0, Y0  Real  I  Coordinates (deg.) of grid point (1,1).   |a
!       SX, SY  Real  I  Grid increments (deg.).                   /
!       NX      Int. I/O Record length.                            \
!       X0      Real I/O Undefined value.                          /b
!       IERR    Int.  O  Error indicator.
!                         0 : No errors.
!                         1 : Illegal INXOUT.
!                         2 : Illegal ID.
!                         3 : Error in opening file.
!                         4 : Write error in file.
!                         5 : Read error in file.
!                         6 : Premature EOF in read.
!                         7 : Unexpected file identifier read.
!                         8 : Unexpected field identifier read.
!                         9 : Unexpected grid dimensions read.
!                        10 : Unexpected data info.
!     ----------------------------------------------------------------
!      a) for output fields.
!      b) for input data.
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      WW3_PREP  Prog.   N/A    Input data preprocessor.
!      WW3_SHEL  Prog.   N/A    Basic wave model driver.
!      ......    Prog.   N/A    Any other program that reads or
!                               writes WAVEWATCH III data files.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     See end of subroutine.
!
!  7. Remarks :
!
!     - On read, the ID 'WND' may be changed to 'WNS' (including
!       stability data).
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!     !/T  Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)                 :: NDS, NDST, NDSE, NY
      INTEGER, INTENT(INOUT)              :: NX
      INTEGER, INTENT(OUT)                :: IERR
      REAL, INTENT(IN)                    :: Y0, SX, SY
      REAL, INTENT(INOUT)                 :: X0
      CHARACTER(LEN=3), INTENT(INOUT)     :: IDFLD
      CHARACTER*(*), INTENT(IN)           :: INXOUT
      CHARACTER*(*), INTENT(IN), OPTIONAL :: FEXT, FPRE
      LOGICAL, INTENT(IN), OPTIONAL       :: FHDR
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: NXT, NYT, I
      REAL                    :: X0T, Y0T, SXT, SYT, TOL
      LOGICAL                 :: WRITE
      CHARACTER(LEN=3)        :: TSFLD
      CHARACTER(LEN=11)       :: FORM = 'UNFORMATTED'
      CHARACTER(LEN=13)       :: TSSTR, IDSTR = 'WAVEWATCH III'
      CHARACTER(LEN=20)       :: TEMPXT
      CHARACTER(LEN=30)       :: FNAME
      LOGICAL                 :: FDHDR = .TRUE.
!
! 'FORM' is used for initial testing of new files only.
!/
!/ ------------------------------------------------------------------- /
!/
!
! test input parameters ---------------------------------------------- *
!
      IF (INXOUT.NE.'READ' .AND. INXOUT.NE.'WRITE') GOTO 801
      IF ( IDFLD.NE.'LEV' .AND. IDFLD.NE.'CUR' .AND.                  &
           IDFLD.NE.'WND' .AND. IDFLD.NE.'WNS' .AND.                  &
           IDFLD.NE.'ICE' .AND. IDFLD.NE.'DT0' .AND.                  &
           IDFLD.NE.'DT1' .AND. IDFLD.NE.'DT2' )    GOTO 802
!
      IF ( PRESENT(FEXT) ) THEN
          TEMPXT = FEXT
          I      = LEN_TRIM(FEXT)
        ELSE
          TEMPXT = 'ww3'
          I      = 3
        END IF
!
      IF ( PRESENT(FHDR) ) THEN
          FDHDR = FHDR
        END IF
!
! Set internal variables --------------------------------------------- *
!
      IF ( IDFLD.EQ.'LEV' ) THEN
          FNAME = 'level.' // TEMPXT(:I)
          I     = I + 6
        ELSE IF ( IDFLD.EQ.'CUR' ) THEN
          FNAME = 'current.' // TEMPXT(:I)
          I     = I + 8
        ELSE IF ( IDFLD.EQ.'WND' .OR. IDFLD.EQ.'WNS' ) THEN
          FNAME = 'wind.' // TEMPXT(:I)
          I     = I + 5
        ELSE IF ( IDFLD.EQ.'ICE' ) THEN
          FNAME = 'ice.' // TEMPXT(:I)
          I     = I + 4
        ELSE IF ( IDFLD.EQ.'DT0' ) THEN
          FNAME = 'data0.' // TEMPXT(:I)
          I     = I + 6
        ELSE IF ( IDFLD.EQ.'DT1' ) THEN
          FNAME = 'data1.' // TEMPXT(:I)
          I     = I + 6
        ELSE
          FNAME = 'data2.' // TEMPXT(:I)
          I     = I + 6
        END IF
!
      WRITE  = INXOUT .EQ. 'WRITE'
!
! Open file ---------------------------------------------------------- *
!
      IF ( WRITE ) THEN
          IF ( PRESENT(FPRE) ) THEN
              OPEN (NDS,FILE=FPRE//FNAME(:I),FORM=FORM,ERR=803,       &
                    IOSTAT=IERR)
            ELSE
              OPEN (NDS,FILE=FNAME(:I),FORM=FORM,ERR=803,IOSTAT=IERR)
            END IF
        ELSE
          IF ( PRESENT(FPRE) ) THEN
              OPEN (NDS,FILE=FPRE//FNAME(:I),FORM=FORM,               &
                    STATUS='OLD',ERR=803,IOSTAT=IERR)
            ELSE
              OPEN (NDS,FILE=FNAME(:I),FORM=FORM,                     &
                    STATUS='OLD',ERR=803,IOSTAT=IERR)
            END IF
        END IF
!
! Process test data -------------------------------------------------- *
!
      IF ( WRITE ) THEN
          IF ( FDHDR ) THEN
             IF ( FORM .EQ. 'UNFORMATTED' ) THEN
                 WRITE (NDS,ERR=804,IOSTAT=IERR)                      &
                       IDSTR, IDFLD, NX, NY, X0, Y0, SX, SY
               ELSE
                 WRITE (NDS,900,ERR=804,IOSTAT=IERR)                  &
                       IDSTR, IDFLD, NX, NY, X0, Y0, SX, SY
               END IF
            END IF
        ELSE
          IF ( FORM .EQ. 'UNFORMATTED' ) THEN
              READ (NDS,END=806,ERR=805,IOSTAT=IERR)                  &
                    TSSTR, TSFLD, NXT, NYT, X0T, Y0T, SXT, SYT
            ELSE
              READ (NDS,900,END=806,ERR=805,IOSTAT=IERR)              &
                    TSSTR, TSFLD, NXT, NYT, X0T, Y0T, SXT, SYT
            END IF
          IF ( IDSTR .NE. TSSTR ) GOTO 807
          IF ( IDFLD.EQ.'WND' .AND. TSFLD.EQ.'WNS' ) THEN
              IDFLD  = TSFLD
            END IF
          IF ( IDFLD .NE. TSFLD ) GOTO 808
          IF ( IDFLD(1:2) .NE. 'DT' ) THEN
              TOL    = 1.E-4 * MIN(SX,SY)
              IF ( NX.NE.NXT .OR. NY.NE.NYT .OR.                      &
                   ABS(X0-X0T).GT.TOL .OR. ABS(Y0-Y0T).GT.TOL .OR.    &
                   ABS(SX-SXT).GT.TOL .OR. ABS(SY-SYT).GT.TOL ) GOTO 809
            ELSE
              NX     = NXT
              X0     = X0T
            END IF
        END IF
!
! File OK ------------------------------------------------------------ *
!
      IERR   = 0
      RETURN
!
! Error escape locations
!
  801 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1001) INXOUT
      IERR   = 1
      RETURN
!
  802 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1002) IDFLD
      IERR   = 2
      RETURN
!
  803 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1003) IDFLD, IERR
      IERR   = 3
      RETURN
!
  804 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1004) IDFLD, IERR
      IERR   = 4
      RETURN
!
  805 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1005) IDFLD, IERR
      IERR   = 5
      RETURN
!
  806 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1006) IDFLD
      IERR   = 6
      RETURN
!
  807 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1007) TSSTR, IDSTR
      IERR   = 7
      RETURN
!
  808 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1008) TSFLD, IDFLD
      IERR   = 8
      RETURN
!
  809 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1009)                      &
                                NXT, NYT, X0T, Y0T, SXT, SYT,   &
                                NX , NY , X0 , Y0 , SX , SY
      IERR   = 9
      RETURN
!
! Formats
!
  900 FORMAT (1X,A13,1X,A3,2I12,4E12.5)
!
 1001 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDO : '/         &
               '     ILLEGAL INXOUT STRING : ',A/)
 1002 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDO : '/         &
               '     ILLEGAL FIELD ID STRING : ',A/)
 1003 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDO : '/         &
               '     ERROR IN OPENING ',A,' FILE, IOSTAT =',I6/)
 1004 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDO : '/         &
               '     ERROR IN WRITING TO ',A,' FILE, IOSTAT =',I6/)
 1005 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDO : '/         &
               '     ERROR IN READING ',A,' FILE, IOSTAT =',I6/)
 
 1006 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDO : '/         &
               '     PREMATURE END OF ',A,' FILE'/)
 1007 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDO : '/         &
               '     ILLEGAL FILE ID STRING >',A,'<'/           &
               '                  SHOULD BE >',A,'<'/)
 1008 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDO : '/         &
               '     ILLEGAL FIELD ID STRING >',A,'<'/          &
               '                   SHOULD BE >',A,'<'/)
 1009 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDO : '/         &
         '     INCOMPATIBLE GRID DATA : ',2(1X,I4),4(1X,F11.3)/ &
         '                  SHOULD BE : ',2(1X,I4),4(1X,F11.3)/)
!
!/
!/ End of W3FLDO  ---------------------------------------------------- /
!/
      END SUBROUTINE W3FLDO
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3FLDG (INXOUT, IDFLD, NDS, NDST, NDSE, MX, MY,      &
                         NX, NY, T0, TN, TF0, FX0, FY0, FA0,          &
                         TFN, FXN, FYN, FAN, IERR)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         05-Jul-2005 |
!/                  +-----------------------------------+
!/
!/    15-Jan-1999 : Final FORTRAN 77                    ( version 1.18 )
!/    30-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    05-Jul-2005 : Correct first level/ice.            ( version 3.07 )
!/
!  1. Purpose :
!
!     Update input fields in the WAVEWATCH III generic shell from a
!     WAVEWATCH III shell data file or write from preprocessor.
!
!  2. Method :
!
!     Read from file opened by W3FLDO.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       INXOUT  C*(*)  I   Test string for read/write, valid are:
!                          'READ' and 'WRITE'.
!       IDFLD   C*3    I   ID string for field type, valid are:
!                          'LEV', 'CUR', 'WND', 'WNS'  and 'ICE'.
!       NDS     Int.   I   Dataset number for fields file.
!       NDST    Int.   I   Dataset number for test output.
!       NDSE    Int.   I   Dataset number for error output.
!                          (No error output if NDSE < 0 ).
!       MX,MY   Int.   I   Array dimensions output fields.
!       NX,NY   Int.   I   Discrete grid dimensions.
!       T0-N    I.A.   I   Time interval considered (dummy for write).
!       TF0-N   I.A.  I/O  Field times (TFN dummy for write).
!       Fxx     R.A.  I/O  Input fields (FxN dummy for write).
!       IERR    Int.   O   Error indicator,
!                          -1 Past last data
!                           0 OK,
!                           1 : Illegal INXOUT.
!                           2 : Illegal IDFLD.
!                           3 : Error in writing time.
!                           4 : Error in writing field.
!                           5 : Error in reading time.
!                           6 : Premature EOF reading field.
!                           7 : Error reading field.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr.   Id.    Subroutine tracing.
!      TICK21    Subr. W3TIMEMD Advance time.
!      DSEC21    Func.   Id.    Difference between times.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      WW3_PREP  Prog.   N/A    Input data preprocessor.
!      WW3_SHEL  Prog.   N/A    Basic wave model driver.
!      ......    Prog.   N/A    Any other program that reads or
!                               writes WAVEWATCH III data files.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     See end of subroutine.
!
!  7. Remarks :
!
!     - Saving of previous fields needed only for reading of 2-D fields.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!     !/T  Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      USE W3TIMEMD
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)          :: NDS, NDST, NDSE, MX, MY,        &
                                      NX, NY, T0(2), TN(2)
      INTEGER, INTENT(INOUT)       :: TF0(2), TFN(2)
      INTEGER, INTENT(OUT)         :: IERR
      REAL, INTENT(INOUT)          :: FX0(MX,MY), FY0(MX,MY),         &
                                      FXN(MX,MY), FYN(MX,MY),         &
                                      FA0(MX,MY), FAN(MX,MY)
      CHARACTER*(*), INTENT(IN)    :: INXOUT
      CHARACTER(LEN=3), INTENT(IN) :: IDFLD
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IX, IY, J, ISTAT
      REAL                    :: DTTST
      LOGICAL                 :: WRITE, FL2D, FLFRST, FLST
!/
!/ ------------------------------------------------------------------- /
!/
!/
      IERR   = 0
!
! test input parameters ---------------------------------------------- *
!
      IF (INXOUT.NE.'READ' .AND. INXOUT.NE.'WRITE') GOTO 801
      IF ( IDFLD.NE.'LEV' .AND. IDFLD.NE.'CUR' .AND.                  &
           IDFLD.NE.'WND' .AND. IDFLD.NE.'WNS' .AND.                  &
           IDFLD.NE.'ICE' )    GOTO 802
!
! Set internal variables --------------------------------------------- *
!
      WRITE  = INXOUT .EQ. 'WRITE'
      FL2D   = IDFLD.EQ.'CUR' .OR. IDFLD.EQ.'WND' .OR. IDFLD.EQ.'WNS'
      FLST   = IDFLD.EQ.'WNS'
      FLFRST = TFN(1) .EQ. -1
!
! Loop over times / fields ========================================== *
!
      DO
!
! Shift fields
!
        IF ( (.NOT.WRITE) .AND. FL2D ) THEN
!
            TF0(1) = TFN(1)
            TF0(2) = TFN(2)
            IF ( TFN(1) .NE. -1 ) THEN
                DO IX=1, NX
                  DO IY=1, NY
                    FX0(IX,IY) = FXN(IX,IY)
                    FY0(IX,IY) = FYN(IX,IY)
                    END DO
                  IF( FLST ) THEN
                      DO IY=1, NY
                        FA0(IX,IY) = FAN(IX,IY)
                        END DO
                    END IF
                  END DO
              END IF
!
          END IF
!
! Process fields, write --------------------------------------------- *
!
        IF ( WRITE ) THEN
!
            WRITE (NDS,ERR=803,IOSTAT=ISTAT) TF0
            IF ( .NOT. FL2D ) THEN
                J      = 1
                WRITE (NDS,ERR=804,IOSTAT=ISTAT)                      &
                           ((FA0(IX,IY),IX=1,NX),IY=1,NY)
              ELSE
                J      = 1
                WRITE (NDS,ERR=804,IOSTAT=ISTAT)                      &
                           ((FX0(IX,IY),IX=1,NX),IY=1,NY)
                J      = 2
                WRITE (NDS,ERR=804,IOSTAT=ISTAT)                      &
                           ((FY0(IX,IY),IX=1,NX),IY=1,NY)
                J      = 3
                IF ( FLST ) WRITE (NDS,ERR=804,IOSTAT=ISTAT)          &
                           ((FA0(IX,IY),IX=1,NX),IY=1,NY)
              END IF
!
            EXIT
!
! Process fields, read ---------------------------------------------- *
!
          ELSE
!
            READ (NDS,END=800,ERR=805,IOSTAT=ISTAT) TFN
            IF ( .NOT. FL2D ) THEN
                J      = 1
                READ (NDS,END=806,ERR=807,IOSTAT=ISTAT)               &
                           ((FAN(IX,IY),IX=1,NX),IY=1,NY)
              ELSE
                J      = 1
                READ (NDS,END=806,ERR=807,IOSTAT=ISTAT)               &
                           ((FXN(IX,IY),IX=1,NX),IY=1,NY)
                J      = 2
                READ (NDS,END=806,ERR=807,IOSTAT=ISTAT)               &
                           ((FYN(IX,IY),IX=1,NX),IY=1,NY)
                J      = 3
                IF ( FLST ) READ (NDS,END=806,ERR=807,IOSTAT=ISTAT)   &
                           ((FAN(IX,IY),IX=1,NX),IY=1,NY)
              END IF
!
! Check time, branch back if necessary
!
            DTTST  = DSEC21 ( T0 , TFN )
            IF ( .NOT.FL2D .AND. FLFRST .AND. DTTST .EQ. 0. ) EXIT
            IF ( DTTST .GT. 0. ) EXIT
!
          END IF
!
        END DO
!
! Branch point for EOF in current and wind
!
  300 CONTINUE
!
! Check first field
!
      IF ( .NOT.WRITE .AND. FL2D .AND. TF0(1) .EQ. -1 ) THEN
!
          TF0(1) = T0(1)
          TF0(2) = T0(2)
!
          DO IX=1, NX
            DO IY=1, NY
              FX0(IX,IY) = FXN(IX,IY)
              FY0(IX,IY) = FYN(IX,IY)
              END DO
            IF( FLST ) THEN
                DO IY=1, NY
                  FA0(IX,IY) = FAN(IX,IY)
                  END DO
              END IF
            END DO
!
        END IF
!
! Branch point for EOF in level and ice
!
  500 CONTINUE
!
! Process fields, end ----------------------------------------------- *
!
      RETURN
!
! EOF escape location
!
  800 CONTINUE
      IERR   = -1
!
      IF ( FL2D ) THEN
          TFN(1) = TN(1)
          TFN(2) = TN(2)
          CALL TICK21 ( TFN , 1. )
        END IF
!
      IF ( FL2D ) THEN
          GOTO 300
        ELSE
          GOTO 500
        END IF
!
! Error escape locations
!
  801 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1001) INXOUT
      IERR   = 1
      RETURN
!
  802 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1002) IDFLD
      IERR   = 2
      RETURN
!
  803 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1003) ISTAT
      IERR   = 3
      RETURN
!
  804 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1004) J, ISTAT
      IERR   = 4
      RETURN
!
  805 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1005) ISTAT
      IERR   = 5
      RETURN
!
  806 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1006) J, ISTAT
      IERR   = 6
      RETURN
!
  807 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1007) J, ISTAT
      IERR   = 7
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDG : '/               &
               '     ILLEGAL INXOUT STRING : ',A/)
 1002 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDG : '/               &
               '     ILLEGAL FIELD ID STRING : ',A/)
 1003 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDG : '/               &
               '     ERROR IN WRITING TIME, IOSTAT =',I6/)
 1004 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDG : '/               &
               '     ERROR IN WRITING FIELD ',I1,', IOSTAT =',I6/)
 1005 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDG : '/               &
               '     ERROR IN READING TIME, IOSTAT =',I6/)
 1006 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDG : '/               &
               '     PRMATURE EOF READING FIELD ',I1,', IOSTAT =',I6/)
 1007 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDG : '/               &
               '     ERROR IN READING FIELD ',I1,', IOSTAT =',I6/)
!
!/
!/ End of W3FLDG ----------------------------------------------------- /
!/
      END SUBROUTINE W3FLDG
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3FLDD (INXOUT, IDFLD, NDS, NDST, NDSE, TIME, TD,    &
                         NR, ND, NDOUT, DATA, IERR )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         24-Jan-2002 |
!/                  +-----------------------------------+
!/
!/    24-Jan-2002 : Origination.                        ( version 2.17 )
!/
!  1. Purpose :
!
!     Update assimilation data in the WAVEWATCH III generic shell from
!     a WAVEWATCH III shell data file or write from preprocessor.
!
!  2. Method :
!
!     Read from file opened by W3FLDO.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       INXOUT  C*(*)  I   Test string for read/write, valid are:
!                           'WRITE'  Write a data field to file.
!                           'SIZE'   Get the number of records of
!                                    next data set.
!                           'READ'   Read the data set found by
!                                    'SIZE' after allocating proper
!                                    data array.
!       IDFLD   C*3    I   ID string for field type, valid are:
!                          'DT0', 'DT1', and 'DT2'.
!       NDS     Int.   I   Dataset number for fields file.
!       NDST    Int.   I   Dataset number for test output.
!       NDSE    Int.   I   Dataset number for error output.
!                          (No error output if NDSE < 0 ).
!       TIME    I.A.   I   Minimum time for data.
!       TD      I.A.  I/O  Data time.
!       NR,ND   Int.   I   Array dimensions.
!       NDOUT   Int.   O   Number of data to be read next.
!       DATA    R.A.  I/O  Data array.
!       IERR    Int.   O   Error indicator,
!                          -1 Past last data
!                           0 OK,
!                           1 : Illegal INXOUT.
!                           2 : Illegal IDFLD.
!                           3 : Error in writing time.
!                           4 : Error in writing data.
!                           5 : Error in reading time.
!                           6 : Premature EOF reading data.
!                           7 : Error reading data.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr.   Id.    Subroutine tracing.
!      TICK21    Subr. W3TIMEMD Advance time.
!      DSEC21    Func.   Id.    Difference between times.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      WW3_PREP  Prog.   N/A    Input data preprocessor.
!      WW3_SHEL  Prog.   N/A    Basic wave model driver.
!      ......    Prog.   N/A    Any other program that reads or
!                               writes WAVEWATCH III data files.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     See end of subroutine.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!     !/T  Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      USE W3TIMEMD
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)          :: NDS, NDST, NDSE, TIME(2), NR, ND
      INTEGER, INTENT(INOUT)       :: TD(2), NDOUT
      INTEGER, INTENT(OUT)         :: IERR
      REAL, INTENT(INOUT)          :: DATA(NR,ND)
      CHARACTER*(*), INTENT(IN)    :: INXOUT
      CHARACTER(LEN=3), INTENT(IN) :: IDFLD
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: ISTAT, NRT
      REAL                    :: DTTST
      LOGICAL                 :: WRITE, SIZE
!/
!/ ------------------------------------------------------------------- /
!/
!/
      IERR   = 0
!
! test input parameters ---------------------------------------------- *
!
      IF ( INXOUT.NE.'READ' .AND. INXOUT.NE.'WRITE' .AND.             &
           INXOUT.NE.'SIZE' ) GOTO 801
      IF ( IDFLD.NE.'DT0' .AND. IDFLD.NE.'DT1' .AND.                  &
           IDFLD.NE.'DT2' )    GOTO 802
!
! Set internal variables --------------------------------------------- *
!
      WRITE  = INXOUT .EQ. 'WRITE'
      SIZE   = INXOUT .EQ. 'SIZE'
!
! Process fields, write --------------------------------------------- *
!
      IF ( WRITE ) THEN
!
          WRITE (NDS,ERR=803,IOSTAT=ISTAT) TD, ND
          WRITE (NDS,ERR=804,IOSTAT=ISTAT) DATA
!
! Process fields, read size ----------------------------------------- *
!
        ELSE IF ( SIZE ) THEN
!
  100     CONTINUE
          READ (NDS,END=800,ERR=805,IOSTAT=ISTAT) TD, NDOUT
!
! Check time, read and branch back if necessary
!
          DTTST  = DSEC21 ( TIME , TD )
          IF ( DTTST.LT.0. .OR. NDOUT.EQ.0 ) THEN
              IF (NDOUT.GT.0) READ (NDS,END=806,ERR=807,IOSTAT=ISTAT)
              GOTO 100
           END IF
!
! Process fields, read data ----------------------------------------- *
!
        ELSE
!
          READ (NDS,END=806,ERR=807,IOSTAT=ISTAT) DATA
        END IF
!
! Process fields, end ----------------------------------------------- *
!
      RETURN
!
! EOF escape location
!
  800 CONTINUE
      IERR   = -1
      RETURN
!
! Error escape locations
!
  801 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1001) INXOUT
      IERR   = 1
      RETURN
!
  802 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1002) IDFLD
      IERR   = 2
      RETURN
!
  803 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1003) ISTAT
      IERR   = 3
      RETURN
!
  804 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1004) ISTAT
      IERR   = 4
      RETURN
!
  805 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1005) ISTAT
      IERR   = 5
      RETURN
!
  806 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1006) ISTAT
      IERR   = 6
      RETURN
!
  807 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1007) ISTAT
      IERR   = 7
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDD : '/               &
               '     ILLEGAL INXOUT STRING : ',A/)
 1002 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDD : '/               &
               '     ILLEGAL FIELD ID STRING : ',A/)
 1003 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDD : '/               &
               '     ERROR IN WRITING TIME, IOSTAT =',I6/)
 1004 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDD : '/               &
               '     ERROR IN WRITING DATA, IOSTAT =',I6/)
 1005 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDD : '/               &
               '     ERROR IN READING TIME, IOSTAT =',I6/)
 1006 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDD : '/               &
               '     PRMATURE EOF READING DATA, IOSTAT =',I6/)
 1007 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDD : '/               &
               '     ERROR IN READING DATA, IOSTAT =',I6/)
!
!/
!/ End of W3FLDD ----------------------------------------------------- /
!/
      END SUBROUTINE W3FLDD
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3FLDP ( NDSM, NDST, NDSE, IERR, MX, MY, NX, NY,     &
                          X0, SX, Y0, SY, MAPOVR, ILAND, MXI, MYI,    &
                          NXI, NYI, CLOSED, ALAT, ALON, MASK,         &
                          RD11, RD21, RD12, RD22, IX1, IX2, IY1, IY2 )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         30-Nov-1999 |
!/                  +-----------------------------------+
!/
!/    08-Feb-1999 : Final FORTRAN 77                    ( version 1.18 )
!/    30-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/
!  1. Purpose :
!
!     General purpose routine for interpolating data of and irregular
!     lat-long grid given by ALA AND ALO to a regular lat-long grid
!     given by X0, Y0, SX and SY.
!
!  2. Method :
!
!     Steepest descent search.
!     Bi-linear interpolation.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NDSM    Int.  I  Unit number message  output (disabled if 0).
!       NDST    Int.  I  Unit number test output.
!       NDSE    Int.  I  Unit number error output.
!       IERR    Int.  O  Error indicator (number of lost points due
!                        to ap conflicts).
!       MX,MY   Int.  I  Array dimensions for output type arrays.
!       NX,NY   Int.  I  Id. actual field syze.
!       X0,Y0   Real  I  (long,lat) for array point (1,1) of output
!                        arrays.
!       SX,SY   Real  I  Longitude and latitude grid increments.
!       MAPOVR  I.A. I/O Overlay map, the value of a grid point is
!                        incremeted by 1 of the corresponding grid
!                        point of the output grid is covered by the
!                        input grid. Land points are masked out by
!                        setting them to ILAND.
!       ILAND   Int.  I  Value for land points in MAPOVR (typically<0)
!       MXI,MYI Int.  I  Array dimensions for input fields.
!       NXI,NYI Int.  I  Id. actual field sizes.
!       CLOSED  Log.  I  Flag for closed longitude range in input.
!       ALAT    R.A.  I  Latitude of input grid (degr.).
!       ALON    R.A. I/O Longitudes of input grid (degr.).
!                        (will be modified if CLOSED)
!       MASK    I.A.  I  Land-sea mask for input field (0=land).
!       RDnn    R.A.  O  Interpolation factors (see below).
!       IXn,IYn I.A.  O  Interpolation addresses (see below).
!     ----------------------------------------------------------------
!
!                             RD12|          |RD22
!                     IY2       --+----------+--
!                                 |          |
!                                 |          |
!                                 |          |
!                                 |          |
!                     IY1       --+----------+--
!                             RD11|          |RD21
!
!                                IX1        IX2
!
!     Internal parameters
!     ----------------------------------------------------------------
!       IPMAX   Int.  Maximum number of passes in descent.    ( DATA )
!       TOL     Real  "Tolerance" for misfit after search.    ( DATA )
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr.   Id.    Subroutine tracing.
!      TICK21    Subr. W3TIMEMD Advance time.
!      DSEC21    Func.   Id.    Difference between times.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      WW3_PREP  Prog.   N/A    Input data preprocessor.
!      ......    Prog.   N/A    Any other program that reads or
!                               writes WAVEWATCH III data files.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!  7. Remarks :
!
!     - The present method will break down for output grids including
!       a pole, but should work for imput grids including a pole.
!     - Land points in the input grid are taken out of the interp.
!       algorithm. If this results in zero weight factors through the
!       interpolation box in the input grid, the closest 2 sea point
!       for an extended 4x4 grid are used for interpolation, weighted
!       threi inverse distance.
!
!  8. Structure :
!
!     -----------------------------------------------------------------
!      1.  Initializations.
!        a Initialize counters and factors.
!        b Get range of input grid.
!        c Determine maximum covered range of output grid.
!      2.  Do for all longitudes
!      --------------------------------------------------------------
!          a Calculate longitude
!          b Shift longitude array if necessary
!          c Find suitable start point along boundary of input grid
!          d Loop over latitudes
!          --------------------------------------------------------
!            1 Check if sea point
!            2 Calculate latitude
!            3 Steepest-descent search
!            ----------------------------------------------------
!              a Aux. data for interpolation grid box
!              b Distances within grid box
!              c Height above lines in input grid
!              d Check if value in grid box
!              e Find new grid box
!              f Check if new itteration required / allowed
!            ----------------------------------------------------
!            4 If point found, calculate/store interpolation data
!            5 Correct interpolation data for input mask
!     -----------------------------------------------------------------
!
!  9. Switches :
!
!     !/S   Enable subroutine tracing.
!
!     !/T   Enable test output.
!     !/T1  Overlay map before and after.
!     !/T2  Test output longitude loop.
!     !/T3  Test output latitude loop.
!     !/T4  Test output steepest descent.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: NDSM, NDST, NDSE, MX, MY, NX, NY,    &
                                 MXI, MYI, NXI, NYI, MASK(MXI,MYI)
      INTEGER, INTENT(INOUT)  :: MAPOVR(MX,MY), ILAND
      INTEGER, INTENT(OUT)    :: IERR, IX1(MX,MY), IX2(MX,MY),        &
                                       IY1(MX,MY), IY2(MX,MY)
      REAL, INTENT(IN)        :: X0, Y0, SX, SY, ALAT(MXI,MYI)
      REAL, INTENT(INOUT)     :: ALON(MXI,MYI)
      REAL, INTENT(OUT)       :: RD11(MX,MY), RD12(MX,MY),            &
                                  RD21(MX,MY), RD22(MX,MY)
      LOGICAL, INTENT(IN)     :: CLOSED
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IX, IY, IY0, IYN, JX, JY, IPASS, JXO,&
                                 JYO, INMASK, J, I, ILOC, II(16),     &
                                 JJ(16), JX1, JX2, JY1, JY2, IFOUND,  &
                                 IMASK, ICOR1
      INTEGER, SAVE           :: IPMAX = 10
      REAL                    :: ALTMIN, ALTMAX, ALNMIN, ALNMAX,      &
                                 ALNAVG, XMIN, XMAX, XAVG, XADD,      &
                                 RDSUM, LO1, LO2, LO3, LO4, LA1, LA2, &
                                 LA3, LA4,   D12, D13, D24, D34,      &
                                 D12X34, D13X24, HT12, HT34, HT13,    &
                                 HT24, RDB1, RDB2, RR(16), RRMAX1,    &
                                 RRMAX2, X, Y
      REAL, SAVE              :: TOL = 0.01
      LOGICAL                 :: F12X34, F13X24, FL1234
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Initializations ------------------------------------------------ *
! 1.a Initialize counters and factors
!
      IERR   = 0
      IFOUND = 0
      IMASK  = 0
      ICOR1  = 0
!
      DO 110, IX=1, NX
        DO 100, IY=1, NY
          RD11(IX,IY) = 0.
          RD12(IX,IY) = 0.
          RD21(IX,IY) = 0.
          RD22(IX,IY) = 0.
          IX1(IX,IY)  = 1
          IX2(IX,IY)  = 1
          IY1(IX,IY)  = 1
          IY2(IX,IY)  = 1
  100     CONTINUE
  110   CONTINUE
!
! 1.b Get range of input grid
!
      ALTMIN = ALAT(1,1)
      ALTMAX = ALTMIN
      ALNMIN = ALON(1,1)
      ALNMAX = ALNMIN
!
      DO 130, IX=1, NXI
        DO 120, IY=1, NYI
          ALTMIN = MIN ( ALTMIN , ALAT(IX,IY) )
          ALTMAX = MAX ( ALTMAX , ALAT(IX,IY) )
          ALNMIN = MIN ( ALNMIN , ALON(IX,IY) )
          ALNMAX = MAX ( ALNMAX , ALON(IX,IY) )
  120     CONTINUE
  130   CONTINUE
!
      ALNAVG = 0.5 * ( ALNMIN + ALNMAX )
      XAVG   = X0 + 0.5*SX*REAL(NX-1)
      XADD   = -360. * REAL(NINT((ALNAVG-XAVG)/360.))
!
      IF ( XADD .NE. 0. ) THEN
          DO 150, IX=1, NXI
            DO 140, IY=1, NYI
              ALON(IX,IY) = ALON(IX,IY) + XADD
  140         CONTINUE
  150       CONTINUE
          ALNMIN = ALNMIN + XADD
          ALNMAX = ALNMAX + XADD
        ENDIF
!
! 1.c Determine maximum covered latitude range
!
      IY0    = NINT( (ALTMIN-Y0)/SY )
      IYN    = NINT( (ALTMAX-Y0)/SY ) + 2
      IY0    = MAX ( 1 , MIN(NY,IY0) )
      IYN    = MAX ( 1 , MIN(NY,IYN) )
!
! 2.  Loop over longitudes ------------------------------------------- *
!
      DO 500, IX=1, NX
!
! 2.a Calculate longitude
!
        X      = X0 + REAL(IX-1)*SX
!
! 2.b Shift longitude (array) if necessary
!
        IF ( CLOSED ) THEN
!
            XMIN   = X - 180.
            XMAX   = X + 180.
            DO 210, JX=1, NXI
              DO 200, JY=1, NYI
                IF (ALON(JX,JY).LT.XMIN) ALON(JX,JY) = ALON(JX,JY)+360.
                IF (ALON(JX,JY).GT.XMAX) ALON(JX,JY) = ALON(JX,JY)-360.
  200           CONTINUE
  210         CONTINUE
!
          ELSE
!
            IF ( X+0.001*SX .LT. ALNMIN ) THEN
                X      = X + 360.*REAL(1+INT((ALNMIN-X)/360.+0.01))
                IF ( X .GT. ALNMAX ) THEN
                  GOTO 500
                  ENDIF
              ENDIF
            IF ( X-0.001*SX .GT. ALNMAX ) THEN
                X      = X - 360.*REAL(1+INT((X-ALNMAX)/360.+0.01))
                IF ( X .LT. ALNMIN ) THEN
                  GOTO 500
                  ENDIF
              ENDIF
!
          ENDIF
!
! 2.c Find suitable start point along boundary of input grid
!
      DO 270, JX=1, NXI-1
        JY     =  1
        IF ((ALON(JX,JY)-X)*(ALON(JX+1,JY)-X).LE.0.                   &
             .AND. ABS(ALON(JX,JY)-X).LT.45. ) GOTO 290
        JY     = NYI
        IF ((ALON(JX,JY)-X)*(ALON(JX+1,JY)-X).LE.0.                   &
             .AND. ABS(ALON(JX,JY)-X).LT.45. ) GOTO 290
  270   CONTINUE
!
      DO 280, JY=1, NYI-1
        JX     =  1
        IF ((ALON(JX,JY)-X)*(ALON(JX,JY+1)-X).LE.0.                   &
             .AND. ABS(ALON(JX,JY)-X).LT.45. ) GOTO 290
        JX     = NXI
        IF ((ALON(JX,JY)-X)*(ALON(JX,JY+1)-X).LE.0.                   &
             .AND. ABS(ALON(JX,JY)-X).LT.45. ) GOTO 290
  280   CONTINUE
!
  290 CONTINUE
      JX     = MIN ( JX , NXI-1 )
      JY     = MIN ( JY , NYI-1 )
!
! 2.d Loop over latitudes  - - - - - - - - - - - - - - - - - - - - - - -
!
        DO 400, IY=IY0, IYN
!
! 2.d.1 Check if sea point
!
          IF ( MAPOVR(IX,IY) .NE. ILAND ) THEN
!
! 2.d.2 Calculate latitude
!
              Y      = Y0 + REAL(IY-1)*SY
!
! 2.d.3 Steepest-descent search                  3      4
!                                                  +--+
!                         counters aux. var.       |  |
!                                                  +--+
!                                                1      2
              IPASS  = 0
!
! ...   Search re-entry point
!
  300         CONTINUE
              IPASS  = IPASS + 1
!
! 2.d.3.a Aux. data for interpolation grid box
!
              LO1    = ALON(JX  ,JY  )
              LO2    = ALON(JX+1,JY  )
              LO3    = ALON(JX  ,JY+1)
              LO4    = ALON(JX+1,JY+1)
              LA1    = ALAT(JX  ,JY  )
              LA2    = ALAT(JX+1,JY  )
              LA3    = ALAT(JX  ,JY+1)
              LA4    = ALAT(JX+1,JY+1)
!
! 2.d.3.b Distances within grid box
!
              D12    = SQRT ( (LO2-LO1)**2 + (LA2-LA1)**2 )
              D13    = SQRT ( (LO3-LO1)**2 + (LA3-LA1)**2 )
              D24    = SQRT ( (LO4-LO2)**2 + (LA4-LA2)**2 )
              D34    = SQRT ( (LO4-LO3)**2 + (LA4-LA3)**2 )
              D12X34 = 0.5 * ( D13 + D24 )
              D13X24 = 0.5 * ( D12 + D34 )
!
! 2.d.3.c Height above lines in input grid
!
              HT12   = ( - (LA2-LA1)*(X-LO1) + (LO2-LO1)*(Y-LA1) ) / &
                                     ( D12 * D12X34 )
              HT34   = ( - (LA4-LA3)*(X-LO3) + (LO4-LO3)*(Y-LA3) ) / &
                                     ( D34 * D12X34 )
              HT13   = (   (LA3-LA1)*(X-LO1) - (LO3-LO1)*(Y-LA1) ) / &
                                     ( D13 * D13X24 )
              HT24   = (   (LA4-LA2)*(X-LO2) - (LO4-LO2)*(Y-LA2) ) / &
                                     ( D24 * D13X24 )
!
! 2.d.3.d Check if value in grid box
!
              F12X34 = HT12*HT34 .LE. 0.
              F13X24 = HT13*HT24 .LE. 0.
              FL1234 = F12X34 .AND. F13X24
!
! 2.d.3.e Find new grid box
!
              JXO    = JX
              JYO    = JY
!
              IF (.NOT.F12X34) THEN
                  IF (HT12.GT.0.) THEN
                      JY = JY + INT (  HT34 ) + 1
                      JY = MIN ( JY , NYI-1 )
                    ELSE
                      JY = JY - INT ( -HT12 ) - 1
                      JY = MAX ( JY , 1 )
                    ENDIF
                ENDIF
!
              IF (.NOT.F13X24) THEN
                  IF (HT13.GT.0.) THEN
                      JX = JX + INT (  HT24 ) + 1
                      JX = MIN ( JX , NXI-1 )
                    ELSE
                      JX = JX - INT ( -HT13 ) - 1
                      JX = MAX ( JX , 1 )
                    ENDIF
                ENDIF
!
! 2.d.3.f Check if new itteration required / allowed
!
              IF ( FL1234 ) GOTO 320
              IF ( JXO.EQ.JX .AND. JYO.EQ.JY ) GOTO 310
              IF ( IPASS .EQ. IPMAX ) GOTO 310
              GOTO 300
!
  310         CONTINUE
              F12X34 = ABS(HT12).LT.TOL .OR. ABS(HT34).LT.TOL
              F13X24 = ABS(HT13).LT.TOL .OR. ABS(HT24).LT.TOL
              FL1234 = F12X34 .AND. F13X24
  320         CONTINUE
!
! 2.d.4 If point found, calculate/store interpolation data
!
              IF ( FL1234 ) THEN
!
                  IX1 (IX,IY) = JX
                  IX2 (IX,IY) = JX + 1
                  IY1 (IX,IY) = JY
                  IY2 (IX,IY) = JY + 1
                  RDB1      = MIN ( 1. , MAX(0.,ABS(HT24)) )
                  RDB2      = MIN ( 1. , MAX(0.,ABS(HT34)) )
                  RD11(IX,IY) =   RDB1    *   RDB2
                  RD12(IX,IY) =   RDB1    * (1.-RDB2)
                  RD21(IX,IY) = (1.-RDB1) *   RDB2
                  RD22(IX,IY) = (1.-RDB1) * (1.-RDB2)
!
! 2.d.5 Correct interpolation data for input mask
!
                  INMASK = 0
!
! ..... Number of points in mask
!
                  IF ( MASK( JX , JY ) .EQ. 0 ) THEN
                      INMASK = INMASK + 1
                      RD11(IX,IY) = 0.
                    ENDIF
!
                  IF ( MASK( JX ,JY+1) .EQ. 0 ) THEN
                      INMASK = INMASK + 1
                      RD12(IX,IY) = 0.
                    ENDIF
!
                  IF ( MASK(JX+1, JY ) .EQ. 0 ) THEN
                      INMASK = INMASK + 1
                      RD21(IX,IY) = 0.
                    ENDIF
!
                  IF ( MASK(JX+1,JY+1) .EQ. 0 ) THEN
                      INMASK = INMASK + 1
                      RD22(IX,IY) = 0.
                    ENDIF
!
! ..... Correction required
!
                  IF ( INMASK .NE. 0 ) THEN
                      IMASK  = IMASK + 1
!
! ..... Summate remaining interpolation factors
!
                      RDSUM  = RD11(IX,IY) + RD12(IX,IY)        &
                             + RD21(IX,IY) + RD22(IX,IY)
!
! ..... Point partially on land mask, correct interpolation factors
!
                      IF ( RDSUM .GT. 1.E-5 ) THEN
!
                          RD11(IX,IY) = RD11(IX,IY) / RDSUM
                          RD12(IX,IY) = RD12(IX,IY) / RDSUM
                          RD21(IX,IY) = RD21(IX,IY) / RDSUM
                          RD22(IX,IY) = RD22(IX,IY) / RDSUM
!
! ..... SUM RDnn = 0, reset interpolation factors
!
                        ELSE
                          RD11(IX,IY) = 0.
                          RD12(IX,IY) = 0.
                          RD21(IX,IY) = 0.
                          RD22(IX,IY) = 0.
!
! ..... SUM RDnn = 0, determine distance to sea points in 4x4 grid
!
                          ILOC        = 0
                          DO 340, I=MAX(1,JX-1), MIN(NXI,JX+2)
                            DO 330, J=MAX(1,JY-1), MIN(NYI,JY+2)
                              IF ( MASK(I,J).NE.0 ) THEN
                                  ILOC     = ILOC + 1
                                  RR(ILOC) = SQRT ( (ALON(I,J)-X)**2 &
                                                  + (ALAT(I,J)-Y)**2 )
                                  II(ILOC) = I
                                  JJ(ILOC) = J
                                ENDIF
  330                         CONTINUE
  340                       CONTINUE
!
                          IF ( ILOC .GT. 0 ) ICOR1  = ICOR1 + 1
!
! ..... One point found, set interpolation data
!
                          IF ( ILOC .EQ. 1 ) THEN
                              RD11(IX,IY) = 1.
                              IX1 (IX,IY) = II(1)
                              IX2 (IX,IY) = II(1)
                              IY1 (IX,IY) = JJ(1)
                              IY2 (IX,IY) = JJ(1)
!
! ..... Several points found, use closest 2
!
                            ELSE IF ( ILOC .GT. 1 ) THEN
                              RRMAX1 = 1.E10
                              RRMAX2 = 1.E10
                              JX1    = 0
                              JY1    = 0
                              DO 350, I=1, ILOC
                                IF ( RR(I) .LT. RRMAX1 ) THEN
                                     RRMAX2 = RRMAX1
                                     JX2    = JX1
                                     JY2    = JY1
                                     RRMAX1 = RR(I)
                                     JX1    = II(I)
                                     JY1    = JJ(I)
                                   ELSE IF ( RR(I) .LT. RRMAX2 ) THEN
                                     RRMAX2 = RR(I)
                                     JX2    = II(I)
                                     JY2    = JJ(I)
                                   ENDIF
  350                           CONTINUE
                              IX1 (IX,IY) = JX1
                              IY1 (IX,IY) = JY1
                              IX2 (IX,IY) = JX2
                              IY2 (IX,IY) = JY2
                              RD11(IX,IY) = RRMAX2 / ( RRMAX1 + RRMAX2 )
                              RD22(IX,IY) = RRMAX1 / ( RRMAX1 + RRMAX2 )
!
! ..... No luck, uncorrectable error
!
                            ELSE
                              IERR   = IERR + 1
                              WRITE (NDSE,910) IX, IY, X, Y,    &
                                     JX, JX+1, JY, JY+1
                            ENDIF
!
                        ENDIF
!
                    ENDIF
!
! 2.d.6 Update overlay map
!
                  MAPOVR(IX,IY) = MAPOVR(IX,IY) + 1
                  IFOUND  = IFOUND + 1
!
                ENDIF
            ENDIF
!
! ... End of loop over latitudes - - - - - - - - - - - - - - - - - - - -
!
  400     CONTINUE
!
! ... End loop over longitudes --------------------------------------- *
!
  500   CONTINUE
!
!     Final output :
!
      IF (NDSM.NE.0) WRITE (NDSM,900) IFOUND, IMASK, ICOR1, IERR
!
      RETURN
!
! Formats
!
  900 FORMAT (/' *** MESSAGE W3FLDP: FINAL SEA POINT COUNT     :',I8/ &
               '                     INTERPOLATION ACROSS SHORE:',I8/ &
               '                     CORRECTED COASTAL POINTS  :',I8/ &
               '                     UNCORRECTABLE C. POINTS   :',I8/)
!
  910 FORMAT ( ' *** WARNING W3FLDP : SEA POINT ON LAND MASK ', &
                    '(COULD NOT BE CORRECTED)'/                 &
               '     COORDINATES IN OUTPUT GRID :',2I4,2F8.2/   &
               '     X-COUNTERS IN INPUT GRID   :',2I4/         &
               '     Y-COUNTERS IN INPUT GRID   :',2I4)
!
 
!
!/
!/ End of W3FLDP ----------------------------------------------------- /
!/
      END SUBROUTINE W3FLDP
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3FLDH (J, NDST, NDSE, MX, MY, NX, NY, T0, TN,       &
                         NH, NHM, THO, HA, HD, HS, TF0, FX0, FY0, FS0,&
                         TFN, FXN, FYN, FSN, IERR)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         05-Jul-2005 |
!/                  +-----------------------------------+
!/
!/    15-Jan-1999 : Final FORTRAN 77                    ( version 1.18 )
!/    30-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    04-Sep-2003 : Bug fix par. list declaration.      ( version 3.04 )
!/    05-Jul-2005 : Correct first level/ice.            ( version 3.07 )
!/
!  1. Purpose :
!
!     Update homogeneous input fields for the WAVEWATCH III generic
!     shell.
!
!  2. Method :
!
!     Variables defining the homogeneous fields are transfered through
!     the parameter list (see section 3).
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       J       Int    I   Field number of input field as in shell.
!                           1 : water levels
!                           2 : currents
!                           3 : winds
!       NDST    Int.   I   Unit number test output.
!       NDSE    Int.   I   Unit number error messages.
!                          (No output if NDSE < 0).
!       MX,MY   Int.   I   Array dimensions output fields.
!       NX,NY   Int.   I   Field dimensions output fields.
!       T0-N    I.A.   I   Time interval considered.
!       NH      Int.  I/O  Number of homogeneous fields J.
!       NHM     Int.   I   Array dimension corresponding to NH.
!       THO     I.A.  I/O  Times for all homogeneous fields left.
!       HA      R.A.  I/O  Id. amplitude.
!       HD      R.A.  I/O  Id. direction (degr., Naut.).
!       HS      R.A.  I/O  Id. air-sea temperature difference (degr.).
!       TF0-N   I.A.  I/O  Times of input fields
!       Fxx     R.A.  I/O  Input fields (X, Y, Scalar)
!       IERR    Int.   O   Error indicator,
!                           0 OK,
!                           1 Illegal field number
!                          -1 Past last data
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr.   Id.    Subroutine tracing.
!      TICK21    Subr. W3TIMEMD Advance time.
!      DSEC21    Func.   Id.    Difference between times.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      WW3_SHEL  Prog.   N/A    Basic wave model driver.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     - See end of subroutine.
!     - Array dimensions not checked.
!
!  7. Remarks :
!
!     - No homogeneous ice fields available.
!     - Previous fields needed only for 2-D fields.
!
!  8. Structure :
!
!       See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!     !/T  Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      USE W3TIMEMD
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: J, NDST, NDSE, MX, MY, NX, NY,       &
                                 T0(2), TN(2), NHM
      INTEGER, INTENT(INOUT)  :: NH, THO(2,4,NHM), TF0(2), TFN(2)
      INTEGER, INTENT(OUT)    :: IERR
      REAL, INTENT(INOUT)     :: HA(NHM,4), HD(NHM,4), HS(NHM,4),     &
                                 FX0(MX,MY), FY0(MX,MY), FS0(MX,MY),  &
                                 FXN(MX,MY), FYN(MX,MY), FSN(MX,MY)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IX, IY, I
      REAL                    :: X, Y, DIR, DTTST, DERA
      LOGICAL                 :: FLFRST
!/
!/ ------------------------------------------------------------------- /
!/
!
      IERR   = 0
      DERA   = ATAN(1.)/45.
!
! Test field ID number for validity
!
      IF ( J.LE.0 .OR. J .GT.3 ) GOTO 801
      FLFRST = TFN(1) .EQ. -1
!
! Loop over times / fields ========================================== *
!
      DO
!
! Shift fields
!
        TF0(1) = TFN(1)
        TF0(2) = TFN(2)
        IF ( TFN(1) .NE. -1 ) THEN
            IF ( J .EQ. 2 ) THEN
                DO IX=1, NX
                  DO IY=1, NY
                    FX0(IX,IY) = FXN(IX,IY)
                    FY0(IX,IY) = FYN(IX,IY)
                    END DO
                  END DO
               ELSE IF ( J .EQ. 3 ) THEN
                DO IX=1, NX
                  DO IY=1, NY
                    FX0(IX,IY) = FXN(IX,IY)
                    FY0(IX,IY) = FYN(IX,IY)
                    FS0(IX,IY) = FSN(IX,IY)
                    END DO
                  END DO
               END IF
          END IF
!
! New field
!
        IF ( NH .NE. 0. ) THEN
            TFN(1) = THO(1,J,1)
            TFN(2) = THO(2,J,1)
!
            IF ( J .EQ. 1 ) THEN
                DO IX=1, NX
                  DO IY=1, NY
                    FSN(IX,IY) = HS(1,J)
                    END DO
                  END DO
              END IF
!
            IF ( J .EQ. 2 ) THEN
                DIR    = ( 270. - HD(1,J) ) * DERA
                X      = HA(1,J) * COS(DIR)
                Y      = HA(1,J) * SIN(DIR)
                DO IX=1, NX
                  DO IY=1, NY
                    FXN(IX,IY) = X
                    FYN(IX,IY) = Y
                    END DO
                  END DO
              END IF
!
            IF ( J .EQ. 3 ) THEN
                DIR    = ( 270. - HD(1,J) ) * DERA
                X      = HA(1,J) * COS(DIR)
                Y      = HA(1,J) * SIN(DIR)
                DO IX=1, NX
                  DO IY=1, NY
                    FXN(IX,IY) = X
                    FYN(IX,IY) = Y
                    FSN(IX,IY) = HS(1,J)
                    END DO
                  END DO
              END IF
!
! Shift data arrays
!
            DO I=1, NH-1
              THO(1,J,I) = THO(1,J,I+1)
              THO(2,J,I) = THO(2,J,I+1)
              HA(I,J)    = HA(I+1,J)
              HD(I,J)    = HD(I+1,J)
              HS(I,J)    = HS(I+1,J)
              END DO
            NH      = NH - 1
!
          ELSE
!
            TFN(1) = TN(1)
            TFN(2) = TN(2)
            CALL TICK21 ( TFN , 1. )
            IERR   = -1
!
          END IF
!
! Check time
!
        DTTST  = DSEC21 ( T0 , TFN )
        IF ( J.EQ.1 .AND. FLFRST .AND. DTTST.EQ.0. ) EXIT
        IF ( DTTST .GT. 0. ) EXIT
        END DO
!
! Check if first field
!
      IF ( J.NE.1 .AND. TFN(1) .EQ. -1 ) THEN
          TF0(1) = T0(1)
          TF0(2) = T0(2)
!
          DO IX=1, NX
            DO IY=1, NY
              FX0(IX,IY) = FXN(IX,IY)
              FY0(IX,IY) = FYN(IX,IY)
              FS0(IX,IY) = FSN(IX,IY)
              END DO
            END DO
        END IF
!
      RETURN
!
! Error escape locations
!
  801 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1001) J
      IERR   = 1
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDH : '/               &
               '     ILLEGAL FIELD ID NR : ',I4/)
!
!/
!/ End of W3FLDH ----------------------------------------------------- /
!/
      END SUBROUTINE W3FLDH
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3FLDM (J, NDST, NDSE, T0, TN, NH, NHM, THO, HA, HD, &
                         TF0, A0, D0, TFN, AN, DN, IERR)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         26-Dec-2002 |
!/                  +-----------------------------------+
!/
!/    26-Dec-2002 : Origination.                        ( version 3.02 )
!/
!  1. Purpose :
!
!     Update moving grid info for the WAVEWATCH III generic
!     shell.
!
!  2. Method :
!
!     Variables defining the homogeneous fields are transfered through
!     the parameter list (see section 3).
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       J       Int    I   Field number, should be 4.
!       NDST    Int.   I   Unit number test output.
!       NDSE    Int.   I   Unit number error messages.
!                          (No output if NDSE < 0).
!       T0-N    I.A.   I   Time interval considered.
!       NH      Int.  I/O  Number of homogeneous fields J.
!       NHM     Int.   I   Array dimension corresponding to NH.
!       THO     I.A.  I/O  Times for all homogeneous fields left.
!       HA      R.A.  I/O  Id. amplitude.
!       HD      R.A.  I/O  Id. direction (degr., Naut.).
!       TF0-N   I.A.  I/O  Times of input fields
!       A/D0/N  R.A.  I/O  Input data.
!       IERR    Int.   O   Error indicator,
!                           0 OK,
!                           1 Illegal field number
!                          -1 Past last data
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr.   Id.    Subroutine tracing.
!      TICK21    Subr. W3TIMEMD Advance time.
!      DSEC21    Func.   Id.    Difference between times.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      WW3_SHEL  Prog.   N/A    Basic wave model driver.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     - See end of subroutine.
!     - Array dimensions not checked.
!
!  7. Remarks :
!
!  8. Structure :
!
!       See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!     !/T  Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      USE W3TIMEMD
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: J, NDST, NDSE, T0(2), TN(2), NHM
      INTEGER, INTENT(INOUT)  :: NH, THO(2,4,NHM), TF0(2), TFN(2)
      INTEGER, INTENT(OUT)    :: IERR
      REAL, INTENT(INOUT)     :: HA(NHM,4), HD(NHM,4), A0, AN, D0, DN
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: I
      REAL                    :: DTTST, DERA
      LOGICAL                 :: FLFRST
!/
!/ ------------------------------------------------------------------- /
!/
!
      IERR   = 0
      DERA   = ATAN(1.)/45.
!
! Test field ID number for validity
!
      IF ( J .NE. 4 ) GOTO 801
      FLFRST = TFN(1) .EQ. -1
!
! Backward branch point ============================================= *
!
  100 CONTINUE
!
! Shift data
!
      TF0(1) = TFN(1)
      TF0(2) = TFN(2)
      IF ( TFN(1) .NE. -1 ) THEN
          A0     = AN
          D0     = DN
        END IF
!
! New field
!
      IF ( NH .NE. 0. ) THEN
          TFN(1) = THO(1,J,1)
          TFN(2) = THO(2,J,1)
          AN     = HA(1,J)
          DN     = ( 90. - HD(1,J) ) * DERA
!
! Shift data arrays
!
          DO I=1, NH-1
            THO(1,J,I) = THO(1,J,I+1)
            THO(2,J,I) = THO(2,J,I+1)
            HA(I,J)    = HA(I+1,J)
            HD(I,J)    = HD(I+1,J)
            END DO
          NH      = NH - 1
!
        ELSE
!
          TFN(1) = TN(1)
          TFN(2) = TN(2)
          CALL TICK21 ( TFN , 1. )
          IERR   = -1
!
        END IF
!
! Check time
!
      DTTST  = DSEC21 ( T0 , TFN )
      IF ( DTTST .LE. 0. ) GOTO 100
!
! Check if first field
!
      IF ( TF0(1).EQ.-1 ) THEN
          TF0(1) = T0(1)
          TF0(2) = T0(2)
          A0     = AN
          D0     = DN
        END IF
!
      RETURN
!
! Error escape locations
!
  801 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,1001) J
      IERR   = 1
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLDM : '/               &
               '     ILLEGAL FIELD ID NR : ',I4/)
!
!/
!/ End of W3FLDM ----------------------------------------------------- /
!/
      END SUBROUTINE W3FLDM
!/
!/ End of module W3FLDSMD -------------------------------------------- /
!/
      END MODULE W3FLDSMD
