!/ ------------------------------------------------------------------- /
      MODULE W3ARRYMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     In this module all service routines for in and output (binary
!     and test) of arrays are gathered.
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      ICOL      Int.  Private  Number of collums four array output
!                               (if not 80, 132 assumed).
!      NFRMAX    Int.  Private  Max number of frequencies in 1D
!                               print plots of spectra.
!     ----------------------------------------------------------------
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      INA2R     Subr. Public   Read 2D real array.
!      INA2I     Subr. Public   Read 2D integer array.
!      OUTA2R    Subr. Public   Write 2D real array.
!      OUTA2I    Subr. Public   Write 2D integer array.
!      OUTREA    Subr. Public   Print out 1D real array.
!      OUTINT    Subr. Public   Print out 1D integer array.
!      OUTMAT    Subr. Public   Print out 2D real array.
!      PRTBLK    Subr. Public   Print a block-type table of a 2D
!                               real array.
!      PRT1DS    Subr. Public   Print plot of 1D spectrum.
!      PRT1DM    Subr. Public   Print plot of 1D spectra.
!      PRT2DS    Subr. Public   Print plot of 2D spectrum.
!      ANGSTR    Subr. PRT2DS   Convert direction to string.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.            ( !/S )
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!  6. Switches :
!
!       !/S    Enable subroutine tracing troughout module.
!       !/T    Switch on test output for INA2R/I and OUTA2R/I.
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
      PUBLIC
!
      INTEGER, PARAMETER, PRIVATE :: ICOL   = 80
      INTEGER, PARAMETER, PRIVATE :: NFRMAX = 50
      INTEGER, PARAMETER, PRIVATE :: NFM2   = NFRMAX+1
!
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE INA2R  (ARRAY, MX, MY, LX, HX, LY, HY,         &
                         NDS, NDST, NDSE, IDFM, RFORM, IDLA, VSC)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-Nov-1999 |
!/                  +-----------------------------------+
!/                                  Based on INAR2D by N.Booij, DUT.
!/
!/    31-Mar-1993 : Final FORTRAN 77                    ( version 1.18 )
!/    29-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/
!  1. Purpose :
!
!     Reads 2-D array of pre-described layout and format.
!
!  3. Parameter list
!     ----------------------------------------------------------------
!       ARRAY   R.A.   O   Array to be read.
!       MX,MY   Int.   I   Declared size of array.
!       LX,HX   Int.   I   Range of x-counters to be read.
!       LY,HY   Int.   I   Range of y-counters to be read.
!       NDS     Int.   I   Unit number for dataset with array.
!       NDST    Int.   I   Unit number for test output.
!       NDSE    Int.   I   Unit number for error messages.
!       IDFM    Int.   I   Format indicator.
!                           IDFM = 1 : Free format.
!                           IDFM = 2 : Fixed format RFORM.
!                           IDFM = 3 : Unformatted.
!       RFORM   C*(*)  I   Format, if IDFM = 2
!       IDLA    Int.   I   Lay out indicator.
!                           IDLA = 1 : Read for IY=LY-HY, IX=LX-HX,
!                                      IX line by IX line.
!                           IDLA = 2 : Idem, one read statement.
!                           IDLA = 3 : Read for IY=HY-LY, IX=LX,HX,
!                                      IX line by IX line.
!                           IDLA = 4 : Idem, one read statement.
!       VSC     Real   I   Scaling factor (multiplication).
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     See mudule documentation.
!
!  5. Called by :
!
!     Any.
!
!  6. Error messages :
!
!     See error escape locations at end of routine.
!
!  8. Structure :
!
!     See comments in code.
!
!  9. Switches :
!
!     !/S   Enable subroutine tracing.
!     !/T   Dump of input parameters in parameter list.
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
      INTEGER, INTENT(IN)     :: MX, MY, LX, HX, LY, HY, NDS, NDST,  &
                                 NDSE, IDFM, IDLA
      REAL, INTENT(IN)        :: VSC
      CHARACTER, INTENT(IN)   :: RFORM*(*)
      REAL, INTENT(OUT)       :: ARRAY(MX,MY)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IIDFM, IIDLA, IX, IY, ISTAT
!/
!/ ------------------------------------------------------------------- /
!/
!
      IF (IDFM.LT.1 .OR. IDFM.GT.3) THEN
          IIDFM = 1
        ELSE
          IIDFM = IDFM
        END IF
      IF (IDLA.LT.1 .OR. IDLA.GT.4) THEN
          IIDLA = 1
        ELSE
          IIDLA = IDLA
        END IF
!
! Free format read :
!
      IF (IIDFM.EQ.1) THEN
          IF (IIDLA.EQ.1) THEN
                DO IY=LY, HY
                  READ (NDS,*,END=800,ERR=801,IOSTAT=ISTAT)     &
                       (ARRAY(IX,IY),IX=LX,HX)
                  END DO
            ELSE IF (IIDLA.EQ.2) THEN
                  READ (NDS,*,END=800,ERR=801,IOSTAT=ISTAT)     &
                       ((ARRAY(IX,IY),IX=LX,HX),IY=LY,HY)
            ELSE IF (IIDLA.EQ.3) THEN
                DO IY=HY, LY, -1
                  READ (NDS,*,END=800,ERR=801,IOSTAT=ISTAT)     &
                       (ARRAY(IX,IY),IX=LX,HX)
                  END DO
            ELSE
                  READ (NDS,*,END=800,ERR=801,IOSTAT=ISTAT)     &
                       ((ARRAY(IX,IY),IX=LX,HX),IY=HY,LY,-1)
            END IF
!
! Fixed format read :
!
        ELSE IF (IIDFM.EQ.2) THEN
          IF (IIDLA.EQ.1) THEN
                DO IY=LY, HY
                  READ (NDS,RFORM,END=800,ERR=801,IOSTAT=ISTAT) &
                       (ARRAY(IX,IY),IX=LX,HX)
                  END DO
            ELSE IF (IIDLA.EQ.2) THEN
                  READ (NDS,RFORM,END=800,ERR=801,IOSTAT=ISTAT) &
                       ((ARRAY(IX,IY),IX=LX,HX),IY=LY,HY)
            ELSE IF (IIDLA.EQ.3) THEN
                DO IY=HY, LY, -1
                  READ (NDS,RFORM,END=800,ERR=801,IOSTAT=ISTAT) &
                       (ARRAY(IX,IY),IX=LX,HX)
                  END DO
            ELSE
                  READ (NDS,RFORM,END=800,ERR=801,IOSTAT=ISTAT) &
                       ((ARRAY(IX,IY),IX=LX,HX),IY=HY,LY,-1)
            END IF
!
! Unformat read :
!
        ELSE
          IF (IIDLA.EQ.1) THEN
                DO IY=LY, HY
                  READ (NDS,END=800,ERR=801,IOSTAT=ISTAT)       &
                       (ARRAY(IX,IY),IX=LX,HX)
                  END DO
            ELSE IF (IIDLA.EQ.2) THEN
                  READ (NDS,END=800,ERR=801,IOSTAT=ISTAT)       &
                       ((ARRAY(IX,IY),IX=LX,HX),IY=LY,HY)
            ELSE IF (IIDLA.EQ.3) THEN
                DO IY=HY, LY, -1
                  READ (NDS,END=800,ERR=801,IOSTAT=ISTAT)       &
                       (ARRAY(IX,IY),IX=LX,HX)
                  END DO
            ELSE
                  READ (NDS,END=800,ERR=801,IOSTAT=ISTAT)       &
                       ((ARRAY(IX,IY),IX=LX,HX),IY=HY,LY,-1)
            END IF
        END IF
!
! Scaling :
!
      DO IX=LX, HX
        DO IY=LY, HY
          ARRAY(IX,IY) = VSC * ARRAY(IX,IY)
          END DO
        END DO
!
      RETURN
!
! Escape locations read errors :
!
  800 CONTINUE
      WRITE (NDSE,900)
      STOP
!
  801 CONTINUE
      WRITE (NDSE,901) ISTAT
      STOP
!
! Formats
!
  900 FORMAT (/' *** ERROR INA2R : '/                           &
               '     PREMATURE END OF FILE'/)
  901 FORMAT (/' *** ERROR INA2R : '/                           &
               '     ERROR IN READING FROM FILE'/               &
               '     IOSTAT =',I5/)
!
!/
!/ End of INA2R  ----------------------------------------------------- /
!/
      END SUBROUTINE INA2R
!/ ------------------------------------------------------------------- /
      SUBROUTINE INA2I  (ARRAY, MX, MY, LX, HX, LY, HY,         &
                         NDS, NDST, NDSE, IDFM, RFORM, IDLA, VSC)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-Nov-1999 |
!/                  +-----------------------------------+
!/                                  Based on INAR2D by N.Booij, DUT.
!/
!/    31-Mar-1993 : Final FORTRAN 77                    ( version 1.18 )
!/    29-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/
!  1. Purpose :
!
!     Like INA2R , integer ARRAY and VSC, see INA2R .
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
      INTEGER, INTENT(IN)     :: MX, MY, LX, HX, LY, HY, NDS, NDST,  &
                                 NDSE, IDFM, IDLA, VSC
      INTEGER, INTENT(OUT)    :: ARRAY(MX,MY)
      CHARACTER, INTENT(IN)   :: RFORM*(*)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IIDFM, IIDLA, IX, IY, ISTAT
!/
!/ ------------------------------------------------------------------- /
!/
!
      IF (IDFM.LT.1 .OR. IDFM.GT.3) THEN
          IIDFM = 1
        ELSE
          IIDFM = IDFM
        END IF
      IF (IDLA.LT.1 .OR. IDLA.GT.4)THEN
          IIDLA = 1
        ELSE
          IIDLA = IDLA
        END IF
!
! Free format read :
!
      IF (IIDFM.EQ.1) THEN
          IF (IIDLA.EQ.1) THEN
                DO IY=LY, HY
                  READ (NDS,*,END=800,ERR=801,IOSTAT=ISTAT)     &
                       (ARRAY(IX,IY),IX=LX,HX)
                  END DO
            ELSE IF (IIDLA.EQ.2) THEN
                  READ (NDS,*,END=800,ERR=801,IOSTAT=ISTAT)     &
                       ((ARRAY(IX,IY),IX=LX,HX),IY=LY,HY)
            ELSE IF (IIDLA.EQ.3) THEN
                DO IY=HY, LY, -1
                  READ (NDS,*,END=800,ERR=801,IOSTAT=ISTAT)     &
                       (ARRAY(IX,IY),IX=LX,HX)
                  END DO
            ELSE
                  READ (NDS,*,END=800,ERR=801,IOSTAT=ISTAT)     &
                       ((ARRAY(IX,IY),IX=LX,HX),IY=HY,LY,-1)
            END IF
!
! Fixed format read :
!
        ELSE IF (IIDFM.EQ.2) THEN
          IF (IIDLA.EQ.1) THEN
                DO IY=LY, HY
                  READ (NDS,RFORM,END=800,ERR=801,IOSTAT=ISTAT) &
                       (ARRAY(IX,IY),IX=LX,HX)
                  END DO
            ELSE IF (IIDLA.EQ.2) THEN
                  READ (NDS,RFORM,END=800,ERR=801,IOSTAT=ISTAT) &
                       ((ARRAY(IX,IY),IX=LX,HX),IY=LY,HY)
            ELSE IF (IIDLA.EQ.3) THEN
                DO IY=HY, LY, -1
                  READ (NDS,RFORM,END=800,ERR=801,IOSTAT=ISTAT) &
                       (ARRAY(IX,IY),IX=LX,HX)
                  END DO
            ELSE
                  READ (NDS,RFORM,END=800,ERR=801,IOSTAT=ISTAT) &
                       ((ARRAY(IX,IY),IX=LX,HX),IY=HY,LY,-1)
            END IF
!
! Unformat read :
!
        ELSE
          IF (IIDLA.EQ.1) THEN
                DO IY=LY, HY
                  READ (NDS,END=800,ERR=801,IOSTAT=ISTAT)       &
                       (ARRAY(IX,IY),IX=LX,HX)
                  END DO
            ELSE IF (IIDLA.EQ.2) THEN
                  READ (NDS,END=800,ERR=801,IOSTAT=ISTAT)       &
                       ((ARRAY(IX,IY),IX=LX,HX),IY=LY,HY)
            ELSE IF (IIDLA.EQ.3) THEN
                DO IY=HY, LY, -1
                  READ (NDS,END=800,ERR=801,IOSTAT=ISTAT)       &
                       (ARRAY(IX,IY),IX=LX,HX)
                  END DO
            ELSE
                  READ (NDS,END=800,ERR=801,IOSTAT=ISTAT)       &
                       ((ARRAY(IX,IY),IX=LX,HX),IY=HY,LY,-1)
            END IF
        END IF
!
! Scaling :
!
      DO IX=LX, HX
        DO IY=LY, HY
          ARRAY(IX,IY) = VSC * ARRAY(IX,IY)
          END DO
        END DO
!
      RETURN
!
! Escape locations read errors :
!
  800 CONTINUE
      WRITE (NDSE,900)
      STOP
!
  801 CONTINUE
      WRITE (NDSE,901) ISTAT
      STOP
!
! Formats
!
  900 FORMAT (/' *** ERROR INA2I : '/                           &
               '     PREMATURE END OF FILE'/)
  901 FORMAT (/' *** ERROR INA2I : '/                           &
               '     ERROR IN READING FROM FILE'/               &
               '     IOSTAT =',I5/)
!
!/
!/ End of INA2I  ----------------------------------------------------- /
!/
      END SUBROUTINE INA2I
!/ ------------------------------------------------------------------- /
      SUBROUTINE OUTA2R (ARRAY, MX, MY, LX, HX, LY, HY,         &
                         NDS, NDST, NDSE, IDFM, RFORM, IDLA, VSC)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         21-Feb-2008 |
!/                  +-----------------------------------+
!/
!/    31-Mar-1993 : Final FORTRAN 77                    ( version 1.18 )
!/    29-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    21-Feb-2008 ; Bug fix IDFM=1, IDLA=2 writing      ( version 3.13 )
!/
!  1. Purpose :
!
!     Writes 2-D array of pre-described layout and format. "Inverse"
!     version of INA2R . For documentation see INA2R .
!
!     N.B. - Array is devided by VSC !
!          - No error trapping on write.
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
      INTEGER, INTENT(IN)     :: MX, MY, LX, HX, LY, HY, NDS, NDST,  &
                                 NDSE, IDFM, IDLA
      REAL, INTENT(IN)        :: VSC, ARRAY(MX,MY)
      CHARACTER, INTENT(IN)   :: RFORM*(*)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IIDFM, IIDLA, IX, IY, ISTAT
!/
!/ ------------------------------------------------------------------- /
!/
!
      IF (IDFM.LT.1 .OR. IDFM.GT.3) THEN
          IIDFM = 1
        ELSE
          IIDFM = IDFM
        END IF
      IF (IDLA.LT.1 .OR. IDLA.GT.4) THEN
          IIDLA = 1
        ELSE
          IIDLA = IDLA
        END IF
!
! Free format write :
!
      IF (IIDFM.EQ.1) THEN
          IF (IIDLA.EQ.1) THEN
              DO IY=LY, HY
                WRITE (NDS,*,ERR=800,IOSTAT=ISTAT)              &
                      (ARRAY(IX,IY)/VSC,IX=LX,HX)
                END DO
            ELSE IF (IIDLA.EQ.2) THEN
              WRITE (NDS,*,ERR=800,IOSTAT=ISTAT)                &
                   ((ARRAY(IX,IY)/VSC,IX=LX,HX),IY=LY,HY)
            ELSE IF (IIDLA.EQ.3) THEN
              DO IY=HY, LY, -1
                WRITE (NDS,*,ERR=800,IOSTAT=ISTAT)              &
                      (ARRAY(IX,IY)/VSC,IX=LX,HX)
                END DO
            ELSE
              WRITE (NDS,*,ERR=800,IOSTAT=ISTAT)                &
                   ((ARRAY(IX,IY)/VSC,IX=LX,HX),IY=HY,LY,-1)
            END IF
!
! Fixed format write :
!
        ELSE IF (IIDFM.EQ.2) THEN
          IF (IIDLA.EQ.1) THEN
              DO IY=LY, HY
                WRITE (NDS,RFORM,ERR=800,IOSTAT=ISTAT)          &
                      (ARRAY(IX,IY)/VSC,IX=LX,HX)
                END DO
            ELSE IF (IIDLA.EQ.2) THEN
              WRITE (NDS,RFORM,ERR=800,IOSTAT=ISTAT)            &
                   ((ARRAY(IX,IY)/VSC,IX=LX,HX),IY=LY,HY)
            ELSE IF (IIDLA.EQ.3) THEN
              DO IY=HY, LY, -1
                WRITE (NDS,RFORM,ERR=800,IOSTAT=ISTAT)          &
                      (ARRAY(IX,IY)/VSC,IX=LX,HX)
                END DO
            ELSE
              WRITE (NDS,RFORM,ERR=800,IOSTAT=ISTAT)            &
                   ((ARRAY(IX,IY)/VSC,IX=LX,HX),IY=HY,LY,-1)
            END IF
!
! Unformat write :
!
        ELSE
          IF (IIDLA.EQ.1) THEN
              DO IY=LY, HY
                WRITE (NDS,ERR=800,IOSTAT=ISTAT)                &
                      (ARRAY(IX,IY)/VSC,IX=LX,HX)
                END DO
            ELSE IF (IIDLA.EQ.2) THEN
              WRITE (NDS,ERR=800,IOSTAT=ISTAT)                  &
                   ((ARRAY(IX,IY)/VSC,IX=LX,HX),IY=LY,HY)
            ELSE IF (IIDLA.EQ.3) THEN
              DO IY=HY, LY, -1
                WRITE (NDS,ERR=800,IOSTAT=ISTAT)                &
                      (ARRAY(IX,IY)/VSC,IX=LX,HX)
                END DO
            ELSE
              WRITE (NDS,ERR=800,IOSTAT=ISTAT)                  &
                    ((ARRAY(IX,IY)/VSC,IX=LX,HX),IY=HY,LY,-1)
            END IF
        END IF
!
      RETURN
!
! Escape locations write errors :
!
  800 CONTINUE
      WRITE (NDSE,900) ISTAT
      STOP
!
! Formats
!
  900 FORMAT (/' *** ERROR OUTA2R : '/                          &
               '     ERROR IN WRITING TO FILE'/                 &
               '     IOSTAT =',I5/)
!
!/
!/ End of OUTA2R ----------------------------------------------------- /
!/
      END SUBROUTINE OUTA2R
!/ ------------------------------------------------------------------- /
      SUBROUTINE OUTA2I (ARRAY, MX, MY, LX, HX, LY, HY,         &
                         NDS, NDST, NDSE, IDFM, RFORM, IDLA, VSC)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-Nov-1999 |
!/                  +-----------------------------------+
!/
!/    31-Mar-1993 : Final FORTRAN 77                    ( version 1.18 )
!/    29-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/
!  1. Purpose :
!
!     Like OUTA2R, integer ARRAY and VSC, see OUTA2R.
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
      INTEGER, INTENT(IN)     :: MX, MY, LX, HX, LY, HY, NDS, NDST,  &
                                 NDSE, IDFM, IDLA, VSC, ARRAY(MX,MY)
      CHARACTER, INTENT(IN)   :: RFORM*(*)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IIDFM, IIDLA, IX, IY, ISTAT
!/
!/ ------------------------------------------------------------------- /
!/
!
      IF (IDFM.LT.1 .OR. IDFM.GT.3) THEN
          IIDFM = 1
        ELSE
          IIDFM = IDFM
        END IF
      IF (IDLA.LT.1 .OR. IDLA.GT.4) THEN
          IIDLA = 1
        ELSE
          IIDLA = IDLA
        END IF
!
! Free format write :
!
      IF (IIDFM.EQ.1) THEN
          IF (IIDLA.EQ.1) THEN
              DO IY=LY, HY
                WRITE (NDS,*,ERR=800,IOSTAT=ISTAT)              &
                      (ARRAY(IX,IY)/VSC,IX=LX,HX)
                END DO
            ELSE IF (IIDLA.EQ.2) THEN
              WRITE (NDS,*,ERR=800,IOSTAT=ISTAT)                &
                   ((ARRAY(IX,IY)/VSC,IX=LX,HX),IY=LY,HY)
            ELSE IF (IIDLA.EQ.3) THEN
              DO IY=HY, LY, -1
                WRITE (NDS,*,ERR=800,IOSTAT=ISTAT)              &
                      (ARRAY(IX,IY)/VSC,IX=LX,HX)
                END DO
            ELSE
              WRITE (NDS,*,ERR=800,IOSTAT=ISTAT)                &
                   ((ARRAY(IX,IY)/VSC,IX=LX,HX),IY=HY,LY,-1)
            END IF
!
! Fixed format write :
!
        ELSE IF (IIDFM.EQ.2) THEN
          IF (IIDLA.EQ.1) THEN
              DO IY=LY, HY
                WRITE (NDS,RFORM,ERR=800,IOSTAT=ISTAT)          &
                      (ARRAY(IX,IY)/VSC,IX=LX,HX)
                END DO
            ELSE IF (IIDLA.EQ.2) THEN
              WRITE (NDS,RFORM,ERR=800,IOSTAT=ISTAT)            &
                   ((ARRAY(IX,IY)/VSC,IX=LX,HX),IY=LY,HY)
            ELSE IF (IIDLA.EQ.3) THEN
              DO IY=HY, LY, -1
                WRITE (NDS,RFORM,ERR=800,IOSTAT=ISTAT)          &
                      (ARRAY(IX,IY)/VSC,IX=LX,HX)
                END DO
            ELSE
              WRITE (NDS,RFORM,ERR=800,IOSTAT=ISTAT)            &
                   ((ARRAY(IX,IY)/VSC,IX=LX,HX),IY=HY,LY,-1)
            END IF
!
! Unformat write :
!
        ELSE
          IF (IIDLA.EQ.1) THEN
              DO IY=LY, HY
                WRITE (NDS,ERR=800,IOSTAT=ISTAT)                &
                      (ARRAY(IX,IY)/VSC,IX=LX,HX)
                END DO
            ELSE IF (IIDLA.EQ.2) THEN
              WRITE (NDS,ERR=800,IOSTAT=ISTAT)                  &
                   ((ARRAY(IX,IY)/VSC,IX=LX,HX),IY=LY,HY)
            ELSE IF (IIDLA.EQ.3) THEN
              DO IY=HY, LY, -1
                WRITE (NDS,ERR=800,IOSTAT=ISTAT)                &
                      (ARRAY(IX,IY)/VSC,IX=LX,HX)
                END DO
            ELSE
              WRITE (NDS,ERR=800,IOSTAT=ISTAT)                  &
                   ((ARRAY(IX,IY)/VSC,IX=LX,HX),IY=HY,LY,-1)
            END IF
        END IF
!
      RETURN
!
! Escape locations write errors :
!
  800 CONTINUE
      WRITE (NDSE,900) ISTAT
      STOP
!
! Formats
!
  900 FORMAT (/' *** ERROR OUTA2I : '/                          &
               '     ERROR IN WRITING TO FILE'/                 &
               '     IOSTAT =',I5/)
!
!/
!/ End of OUTA2I ----------------------------------------------------- /
!/
      END SUBROUTINE OUTA2I
!/ ------------------------------------------------------------------- /
      SUBROUTINE OUTREA (NDS,ARRAY,DIM,ANAME)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-Nov-1999 |
!/                  +-----------------------------------+
!/                        Original versions G. Ph. van Vledder
!/                                          P. H. Willems
!/
!/    29-Mar-1993 : Final FORTRAN 77                    ( version 1.18 )
!/    29-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/
!  1. Purpose :
!
!     Print contents of a 1-D real array, see OUTINT.
!
!/ ------------------------------------------------------------------- /
!/
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: NDS, DIM
      REAL, INTENT(IN)        :: ARRAY(DIM)
      CHARACTER, INTENT(IN)   :: ANAME*(*)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: I, K
!/
!/ ------------------------------------------------------------------- /
!/
!
      WRITE (NDS,8000) ANAME
!
      IF (ICOL.EQ.80) THEN
!
          WRITE (NDS,8005) (I, I=1, 5)
          WRITE (NDS,8010)
          DO K=0, DIM, 5
          IF (DIM-K.GE.5) THEN
              WRITE (NDS,'(1X,I4,A,5E12.4,A)')                  &
                K,'  |',(ARRAY(I),I= K+1, K+5),'  |'
            ELSE
              WRITE (NDS,'(1X,T71,''|'',T2,I4,A,5E12.4)')       &
                K,'  |',(ARRAY(I),I= K+1, DIM)
             END IF
          END DO
          WRITE (NDS,8010)
!
        ELSE
!
          WRITE (NDS,9005) (I, I=1, 10)
          WRITE (NDS,9010)
          DO K=0, DIM, 10
          IF (DIM-K.GE.10) THEN
              WRITE (NDS,'(1X,I4,A,10E12.4,A)')                 &
                K,'  |',(ARRAY(I),I= K+1, K+10),'  |'
            ELSE
              WRITE (NDS,'(1X,T131,''|'',T2,I4,A,10E12.4)')     &
                K,'  |',(ARRAY(I),I= K+1, DIM)
             END IF
          END DO
          WRITE (NDS,9010)
         END IF
!
      RETURN
!
 8000 FORMAT (/,1X,'A R R A Y   D U M P  (REAL) / NAME: ',A)
 8005 FORMAT (8X,5I12)
 8010 FORMAT (7X,'+',62('-'),'+')
 9005 FORMAT (8X,10I12)
 9010 FORMAT (7X,'+',122('-'),'+')
!/
!/ End of OUTREA ----------------------------------------------------- /
!/
      END SUBROUTINE OUTREA
!/ ------------------------------------------------------------------- /
      SUBROUTINE OUTINT ( NDS, IARRAY, DIM, ANAME )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-Mar-1993 |
!/                  +-----------------------------------+
!/                        Original versions G. Ph. van Vledder
!/                                          P. H. Willems
!/
!/    29-Mar-1993 : Final FORTRAN 77                    ( version 1.18 )
!/    29-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/
!  1. Purpose :
!
!     Print contents of a 1-D integer array.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NDS     Int.   I   Output unit number.
!       IARRAY  I.A.   I   Array to be printed.
!       DIM     Int.   I   Number of elements to be printed.
!       ANAME   C*(*)  I   Name of array.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     See mudule documentation.
!
!  5. Called by :
!
!       Anny routine or program.
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
      INTEGER, INTENT(IN)     :: NDS, DIM, IARRAY(DIM)
      CHARACTER, INTENT(IN)   :: ANAME*(*)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: I, K
!/
!/ ------------------------------------------------------------------- /
!/
!
      WRITE (NDS,8000) ANAME
!
!  ------- 80 COLUMNS -----
!
      IF (ICOL.EQ.80) THEN
          WRITE (NDS,8005) (I, I=1, 5)
          WRITE (NDS,8010)
          DO K=0, DIM, 5
          IF (DIM-K.GE.5) THEN
              WRITE (NDS,'(1X,I4,A,5I12,A)')                    &
                K,'  |',(IARRAY(I),I= K+1, K+5),'  |'
            ELSE
              WRITE (NDS,'(1X,T71,''|'',T2,I4,A,5I12)')         &
                K,'  |',(IARRAY(I),I= K+1, DIM)
            END IF
          END DO
          WRITE (NDS,8010)
      ELSE
!
!    ---- 132 COLUMNS ----
!
          WRITE (NDS,9005) (I, I=1, 10)
          WRITE (NDS,9010)
          DO K=0, DIM, 10
          IF (DIM-K.GE.10) THEN
              WRITE (NDS,'(1X,I4,A,10I12,A)')                   &
                K,'  |',(IARRAY(I),I= K+1, K+10),'  |'
            ELSE
              WRITE (NDS,'(1X,T131,''|'',T2,I4,A,10I12)')       &
                K,'  |',(IARRAY(I),I= K+1, DIM)
             END IF
          END DO
          WRITE (NDS,9010)
         END IF
!
      RETURN
!
 8000 FORMAT (/,1X,'A R R A Y   D U M P  (INTEGER) / NAME: ',A)
 8005 FORMAT (8X,5I12)
 8010 FORMAT (7X,'+',62('-'),'+')
 9005 FORMAT (8X,10I12)
 9010 FORMAT (7X,'+',122('-'),'+')
!/
!/ End of OUTINT ----------------------------------------------------- /
!/
      END SUBROUTINE OUTINT
!/ ------------------------------------------------------------------- /
      SUBROUTINE OUTMAT (NDS,A,MX,NX,NY,MNAME)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-Nov-1999 |
!/                  +-----------------------------------+
!/                        Original versions G. Ph. van Vledder
!/
!/    29-Mar-1993 : Final FORTRAN 77                    ( version 1.18 )
!/    29-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/
!  1. Purpose :
!
!     Print contents of a 2-D real array.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NDS     Int.   I   Output unit number.
!       A       R.A.   I   Matrix to be printed.
!       MX      Int.   I   Dimension of first index.
!       NX      Int.   I   Number of points for first index.
!       NY      Int.   I   Number of points for scond index.
!       MNAME   C*(*)  I   Name of matrix.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     See mudule documentation.
!
!  5. Called by :
!
!       Anny routine or program.
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
      INTEGER, INTENT(IN)     :: NDS, MX, NX, NY
      REAL, INTENT(IN)        :: A(MX,NY)
      CHARACTER, INTENT(IN)   :: MNAME*(*)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: LBLOK, NBLOK, IBLOK, IX, IX1, IX2, IY
!/
!/ ------------------------------------------------------------------- /
!/
!
      WRITE(NDS,8000) MNAME
!
!  ------ 80 COLUMNS -----
!
      IF(ICOL.EQ.80) THEN
          LBLOK = 6
          NBLOK = (NX-1)/LBLOK + 1
          DO IBLOK = 1,NBLOK
            IX1 = (IBLOK-1)*LBLOK + 1
            IX2 = IX1 + LBLOK - 1
            IF(IX2.GT.NX) IX2 = NX
            WRITE(NDS,8001) (IX,IX = IX1,IX2)
            WRITE(NDS,8002)
            DO IY = 1,NY
              WRITE(NDS,8003) IY,(A(IX,IY),IX = IX1,IX2)
              END DO
            WRITE(NDS,8002)
            END DO
        ELSE
!
!   ---- 132 COLUMNS ----
!
          LBLOK = 12
          NBLOK = (NX-1)/LBLOK + 1
          DO IBLOK = 1,NBLOK
            IX1 = (IBLOK-1)*LBLOK + 1
            IX2 = IX1 + LBLOK - 1
            IF(IX2.GT.NX) IX2 = NX
            WRITE(NDS,9001) (IX,IX = IX1,IX2)
            WRITE(NDS,9002)
            DO IY = 1,NY
              WRITE(NDS,9003) IY,(A(IX,IY),IX = IX1,IX2)
              END DO
            WRITE(NDS,9002)
            END DO
         END IF
!
      RETURN
!
! Formats
!
 8000 FORMAT(/,1X,' M A T R I X   D U M P  (REAL) / NAME: ',A)
 8001 FORMAT(9X,6I10)
 8002 FORMAT(1X,6X,'+',62('-'),'+')
 8003 FORMAT(1X,T71,'|',T2,I5,' | ',12E10.3)
 9001 FORMAT(9X,12I10)
 9002 FORMAT(1X,6X,'+',122('-'),'+')
 9003 FORMAT(1X,T131,'|',T2,I5,' | ',12E10.3)
!/
!/ End of OUTMAT ----------------------------------------------------- /
!/
      END SUBROUTINE OUTMAT
!/ ------------------------------------------------------------------- /
      SUBROUTINE PRTBLK (NDS, NX, NY, MX, F, MAP, MAP0, FSC,    &
                         IX1, IX2, IX3, IY1, IY2, IY3, PRVAR, PRUNIT)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-Nov-1999 |
!/                  +-----------------------------------+
!/
!/    04-Jun-1996 : Final FORTRAN 77                    ( version 1.18 )
!/    29-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/
!  1. Purpose :
!
!     Print a block-type table of a two-dimensional field using a
!     land-sea array.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NDS     Int.   I   File unit number.
!       NX, NY  Int.   I   X and Y range of arrays.
!       MY      Int.   I   Actual X size of arrays.
!       F       R.A.   I   Array to pr presented.
!       MAP     I.A.   I   Map array for land points.
!       MAP0    Int.   I   Map value for land points in MAP.
!       FSC     Real   I   Scaling factor.
!       IX1-3   Int.   I   Firts, last, increment grid points in X
!                          direction.
!       IY1-3   Int.   I   Id. Y direction.
!       PRVAR   C*(*)  I   Name of variable.
!       PRUNIT  C*(*)  I   Units of spectrum.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     See mudule documentation.
!
!  5. Called by :
!
!     Any program.
!
!  6. Error messages :
!
!     None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     ------------------------------------------------
!       Check if automatic scaling
!       If automatic scaling : get extermata
!       Print heading
!       Print table
!       Print ending
!     ------------------------------------------------
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing using STRACE.
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
      INTEGER, INTENT(IN)     :: NDS, NX, NY, MX, MAP(MX,NY), MAP0,  &
                                 IX1, IX2, IX3, IY1, IY2, IY3
      REAL, INTENT(IN)        :: F(MX,NY), FSC
      CHARACTER, INTENT(IN)   :: PRVAR*(*), PRUNIT*(*)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IX, IY, JJ, JM, K1, LX, I
      REAL                    :: FMAX, RR
      LOGICAL                 :: FLSCLE
      CHARACTER               :: PNUM*5, STRA*5, PNUM2*2, STRA3*3
      DIMENSION               :: PNUM(25), PNUM2(61)
!/
!/ ------------------------------------------------------------------- /
!/
!
! Check scaling
!
      FLSCLE = (FSC.LE.0.)
!
! Extremata
!
      IF (FLSCLE) THEN
          FMAX   = 1.E-15
          DO IX=1, NX
            DO IY=1, NY
              IF ( MAP(IX,IY) .NE. MAP0 )                       &
                   FMAX   = MAX ( FMAX , ABS(F(IX,IY)) )
              END DO
            END DO
        END IF
!
! Normalized print plot -----------------------------------------------
!
      IF (FLSCLE) THEN
!
! Heading
!
          WRITE (NDS,901) PRVAR, FMAX, PRUNIT
!
          STRA   = '     '
          JJ     = 0
          DO IX = IX1, IX2, IX3
            JJ = JJ + 1
            END DO
          LX = JJ
          WRITE (NDS,911)
          WRITE (NDS,912) (IX,IX=IX1,IX2,2*IX3)
          PNUM2(1) = '--'
          WRITE (NDS,910) STRA, ' +', (PNUM2(1), I=1, LX), '-+'
!
! Write table
!
          JM = 0
          DO IY = IY2, IY1, IY3*(-1)
!
            JJ = 0
            DO IX = IX1, IX2, IX3
              JJ = JJ + 1
              IF (MAP(IX,IY).EQ.MAP0) THEN
                  PNUM2(JJ) = '  '
                ELSE
                  RR = 10.*F(IX,IY)/FMAX
                  WRITE (STRA, FMT='(I2,3X)') INT(RR*1.000001)
                  PNUM2(JJ) = STRA(1:2)
                  IF (PNUM2(JJ).EQ.'10' .OR. PNUM2(JJ).EQ.'**' .OR. &
                      F(IX,IY).EQ.FMAX) THEN
                      IF ( RR .LT. 0. ) THEN
                          PNUM2(JJ) = '-*'
                        ELSE
                          PNUM2(JJ) = ' *'
                        END IF
                    END IF
                END IF
              END DO
!
            IF (JM.EQ.0) THEN
                WRITE (STRA, FMT='(I5)') IY
                JM   = 2
              ELSE
                STRA = '     '
                JM   = JM-1
              END IF
!
            LX = JJ
            WRITE (NDS,910) STRA, ' |', (PNUM2(I), I=1, LX), ' |'
            END DO
!
          STRA     = '     '
          PNUM2(1) = '--'
          WRITE (NDS,910) STRA, ' +', (PNUM2(1), I=1, LX), '-+'
          WRITE (NDS,912) (IX,IX=IX1,IX2,2*IX3)
          WRITE (NDS,911)
!
! Non-normalized print plot -------------------------------------------
!
        ELSE
!
! Heading
!
          WRITE (NDS,900) PRVAR, FSC, PRUNIT
!
          JJ = 0
          PNUM(1) = '     '
          DO IX = IX1, IX2, IX3
            JJ = JJ + 1
            END DO
          LX = JJ
          WRITE (NDS,921)
          WRITE (NDS,922) (IX,IX=IX1,IX2,IX3)
          STRA3   = '   '
          PNUM(1) = '-----'
          WRITE (NDS,920) STRA3, ' +', (PNUM(1), I=1, LX), '-+   '
!
! Write table
!
          JM = 0
          DO IY = IY2, IY1, IY3*(-1)
            IF (JM.EQ.0) THEN
                WRITE (STRA3, FMT='(I3)') IY
                JM = 2
              ELSE
                STRA3  = '   '
                JM = JM-1
              END IF
!
            JJ = 0
            DO IX = IX1, IX2, IX3
              JJ = JJ + 1
              IF (MAP(IX,IY).EQ.MAP0) THEN
                  PNUM(JJ) = '     '
                ELSE
                  RR     = F(IX,IY)
                  K1 = NINT (RR / FSC)
                  WRITE (STRA, FMT='(I5)') K1
                  PNUM(JJ) = STRA
                END IF
              END DO
!
            LX = JJ
            WRITE (NDS,920) STRA3, ' |', (PNUM(I), I=1, LX), ' |   '
            END DO
!
          STRA3   = '   '
          PNUM(1) = '-----'
          WRITE (NDS,920) STRA3, ' +', (PNUM(1), I=1, LX), '-+   '
          WRITE (NDS,922) (IX,IX=IX1,IX2,IX3)
          WRITE (NDS,921)
!
        END IF
!
      RETURN
!
! Formats
!
  900 FORMAT (/, ' Variable: ',A,' Units: ',E10.3,1X,A)
  901 FORMAT (/, ' Variable: ',A,' Max.: ',E10.3,1X,A)
!
  910 FORMAT (1X,A5,63A2)
  911 FORMAT (' ')
  912 FORMAT (6X,36I4)
!
  920 FORMAT (1X,A3,A2,25A5)
  921 FORMAT (' ')
  922 FORMAT (6X,25I5)
!/
!/ End of PRTBLK ----------------------------------------------------- /
!/
      END SUBROUTINE PRTBLK
!/ ------------------------------------------------------------------- /
      SUBROUTINE PRT1DS (NDS, NFR, E, FR, UFR, NLINES, FTOPI,    &
                         PRVAR, PRUNIT, PNTNME)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-Nov-1999 |
!/                  +-----------------------------------+
!/
!/    10-Mar-1992 : Final FORTRAN 77                    ( version 1.18 )
!/    29-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/
!  1. Purpose :
!
!     Produces a print plot of a 1-D spectrum.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NDS     Int.   I   File unit number.
!       NFR     Int.   I   Number of frequencies.
!       E       R.A.   I   Spectral densities.
!       FR      R.A.   I   Frequencies.
!       UFR     C*(*)  I   If 'HZ', frequencies in Hz, otherwise in
!                          rad/s (N.B., does not re-scale spectrum).
!       NLINES  Int.   I   Hight of plot in lines.
!       FTOPI   Real   I   Highest value of density in plot,
!                          if FTOPI.LE.0., automatic scaling.
!       PRVAR   C*(*)  I   Name of variable.
!       PRUNIT  C*(*)  I   Units of spectrum.
!       PNTNME  C*(*)  I   Name of location.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     See mudule documentation.
!
!  5. Called by :
!
!       Any routine.
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!     - Paperwidth is "set" by NFRMAX.
!
!  8. Structure :
!
!     ------------------------------------------------
!       Initializations and preparations.
!       Determine maximum of spectra.
!       Scaling / normalization.
!       Printing of spectrum
!       ----------------------------------------------
!         Print ID
!         Print heading
!         Print table
!         Print ending
!     ------------------------------------------------
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing using STRACE.
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
      INTEGER, INTENT(IN)     :: NDS, NFR, NLINES
      REAL, INTENT(IN)        :: FTOPI, E(NFR), FR(NFR)
      CHARACTER, INTENT(IN)   :: PRVAR*(*), PRUNIT*(*), PNTNME*(*),  &
                                 UFR*(*)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: NFRB, IFR, IL, IL0
      REAL, SAVE              :: TOPFAC = 1.1
      REAL                    :: FTOP, RLINES, FACFR, FSC, FLINE,    &
                                 EMAX, EMIN, EXTR, FLOC
      LOGICAL                 :: FLSCLE
      CHARACTER               :: STRA*10, STRA2*2, PNUM2*2
      DIMENSION               :: PNUM2(NFM2)
!/
!/ ------------------------------------------------------------------- /
!/
!
      FTOP   = FTOPI
!
      NFRB   = MIN (NFR,50)
      RLINES = REAL(NLINES)
      FLSCLE = FTOP.LE.0.
!
      IF (UFR.EQ.'HZ') THEN
         FACFR  = 1.
       ELSE
         FACFR  = 0.159155
       END IF
!
! Maximum of 1-D spectrum
!
      EMAX   = 0.
      EMIN   = 0.
!
      DO IFR=1, NFR
        EMAX = MAX ( EMAX , E(IFR) )
        EMIN = MIN ( EMIN , E(IFR) )
        END DO
!
      IF (EMAX.EQ.0. .AND. EMIN.EQ.0.) THEN
          EMAX   =  1.E-20
          EMIN   = -1.E-20
        END IF
!
      IF (EMAX.GT.ABS(EMIN)) THEN
          EXTR   = EMAX
        ELSE
          EXTR   = EMIN
        END IF
!
! Scaling / Normalization
!
      IF (FLSCLE) THEN
          IF (EMAX.GT.ABS(EMIN)) THEN
              FLOC   = EMAX * TOPFAC
              FSC    = FLOC / REAL(NINT(EMAX/(EMAX-EMIN)*RLINES))
            ELSE
              FLOC   = EMIN * TOPFAC
              FSC    = FLOC / REAL(NINT(EMIN/(EMAX-EMIN)*RLINES))
              FLOC   = FTOP + RLINES*FSC
              IF (EMAX.LT.0.01*FSC) FTOP = 0.
            END IF
        ELSE
          FLOC   = FTOP
          FSC    = FLOC  / RLINES
          IF (EMAX*EMIN.LT.0) FSC = 2.*FSC
          IF (EMAX.LT.0.01*FSC) FLOC = 0.
        END IF
!
      IL0   = MOD ( NINT(FLOC/FSC) , 2 ) + 1
!
! Print ID
!
      WRITE (NDS,900) PNTNME, PRVAR, EXTR, PRUNIT
!
! Print heading
!
      FLINE  = FLOC
      IF (MOD(NLINES+IL0,2).EQ.0) THEN
          WRITE (STRA, FMT='(E10.3)') FLINE
        ELSE
          STRA=  '          '
        END IF
!
      DO IFR=1, NFRB
        IF ( NINT( (E(IFR)-FLINE)/FSC ) .EQ.0) THEN
            PNUM2(IFR) = '-*'
          ELSE
            PNUM2(IFR) = '--'
          END IF
        END DO
!
      PNUM2(NFRB+1) = '-+'
      STRA2 = ' +'
      WRITE (NDS,910) STRA, STRA2, (PNUM2(IFR),IFR=1, NFRB+1)
!
! Print table
!
      DO IL = 1, NLINES-1
        FLINE  = FLOC - FSC * REAL(IL)
        IF (ABS(FLINE).LT.0.01*FSC) FLINE = 0.
        IF (MOD(NLINES+IL0-IL,2).EQ.0) THEN
            WRITE (STRA, FMT='(E10.3)') FLINE
            STRA2 =  ' +'
          ELSE
            STRA  =  '          '
            STRA2 =  ' |'
          END IF
        DO IFR=1, NFRB
          IF (ABS(FLINE).LT.0.1*FSC) THEN
              PNUM2(NFRB+1) = '-|'
              IF ( NINT( (E(IFR)-FLINE)/FSC ) .EQ.0) THEN
                  PNUM2(IFR) = '-*'
                ELSE
                  PNUM2(IFR) = '--'
                END IF
            ELSE
              PNUM2(NFRB+1) = ' |'
              IF ( NINT( (E(IFR)-FLINE)/FSC ) .EQ.0) THEN
                  PNUM2(IFR) = ' *'
                ELSE
                  PNUM2(IFR) = '  '
                END IF
            END IF
          END DO
        WRITE (NDS,910) STRA, STRA2, (PNUM2(IFR),IFR=1, NFRB+1)
        END DO
!
! write ending
!
      FLINE  = FLOC - FSC * REAL(IL)
      IF (ABS(FLINE).LT.0.01*FSC) FLINE = 0.
      WRITE (STRA, FMT='(E10.3)') FLINE
        IF (MOD(IL0,2).EQ.0) THEN
            WRITE (STRA, FMT='(E10.3)') FLINE
          ELSE
            STRA  =  '          '
          END IF
      STRA2         = ' +'
      PNUM2(NFRB+1) = '-+'
!
      DO IFR=1, NFRB
        IF ( NINT( (E(IFR)-FLINE)/FSC ) .EQ.0) THEN
            PNUM2(IFR) = '-*'
          ELSE IF ( MOD (IFR-2,4) .EQ. 0 ) THEN
            PNUM2(IFR) = '-|'
          ELSE
            PNUM2(IFR) = '--'
          END IF
        END DO
!
      WRITE (NDS,910) STRA, STRA2, (PNUM2(IFR),IFR=1, NFRB+1)
      WRITE (NDS,911) (FR(IFR)*FACFR,IFR=2,NFRB,4)
      WRITE (NDS,920)
!
      RETURN
!
! Formats
!
  900 FORMAT (/'  Location : ',A                                &
              /'  Spectrum : ',A,'  Extreme value : ',E10.3,1X,A/)
!
  910 FORMAT (A10,A2,60A2)
  911 FORMAT (10X,15F8.3)
!
  920 FORMAT (' ')
!/
!/ End of PRT1DS ----------------------------------------------------- /
!/
      END SUBROUTINE PRT1DS
!/ ------------------------------------------------------------------- /
      SUBROUTINE PRT1DM (NDS, NFR, NE, E, FR, UFR, NLINES, FTOPI, &
                         PRVAR, PRUNIT, PNTNME)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         17-Apr-1992 |
!/                  +-----------------------------------+
!/
!/    17-Apr-1992 : Final FORTRAN 77                    ( version 1.18 )
!/    29-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/
!  1. Purpose :
!
!     Produces a print plot of several 1-D spectra.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NDS     Int.   I   File unit number.
!       NFR     Int.   I   Number of frequencies.
!       NE      Int.   I   Number of spectra.
!       E       R.A.   I   Spectral densities.
!       FR      R.A.   I   Frequencies.
!       UFR     C*     I   If 'HZ', frequencies in Hz, otherwise in
!                          rad/s
!       NLINES  Int.   I   Hight of plot in lines.
!       FTOPI   Real   I   Highest value of density in plot,
!                          if FTOP.LE.0., automatic scaling.
!       PRVAR   C*(*)  I   Name of variable.
!       PRUNIT  C*(*)  I   Units of spectrum.
!       PNTNME  C*(*)  I   Name of location.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     See mudule documentation.
!
!  5. Called by :
!
!       Any routine.
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!     - Paperwidth is "set" by NFRMAX.
!
!  8. Structure :
!
!     ------------------------------------------------
!       Initializations and preparations.
!       Determine maximum of spectrum.
!       Scaling / normalization.
!       Printing of spectrum
!       ----------------------------------------------
!         Print ID
!         Print heading
!         Print table
!         Print ending
!     ------------------------------------------------
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing using STRACE.
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
      INTEGER, INTENT(IN)     :: NDS, NFR, NE, NLINES
      REAL, INTENT(IN)        :: FTOPI, E(NFR,NE), FR(NFR)
      CHARACTER, INTENT(IN)   :: PRVAR*(*), PRUNIT*(*), PNTNME*(*),  &
                                 UFR*(*)
      DIMENSION               :: PRVAR(NE)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER, PARAMETER      :: NFRMAX = 50
      INTEGER, PARAMETER      :: NFM2   = NFRMAX+1
      INTEGER                 :: NFRB, IFR, IE, IL
      REAL, SAVE              :: TOPFAC = 1.1
      REAL                    :: FTOP, RLINES, FACFR, FSC, FLINE,    &
                                 EMAX, EMIN, EXTR, FLOC
      LOGICAL                 :: FLSCLE
      CHARACTER               :: STRA*10, STRA2*2, STRAX*2, PNUM2*2
      DIMENSION               :: PNUM2(NFM2)
!/
!/ ------------------------------------------------------------------- /
!/
!
! Test output, echo input
!
      FTOP   = FTOPI
      NFRB   = MIN (NFR,50)
      RLINES = REAL(NLINES)
      FLSCLE = FTOP.LE.0.
!
      IF (UFR.EQ.'HZ') THEN
         FACFR  = 1.
       ELSE
         FACFR  = 0.159155
       END IF
!
! Maximum of 1-D spectrum
!
      EMAX   = 0.
      EMIN   = 0.
!
      DO IE=1, NE
        DO IFR=1, NFR
          EMAX = MAX ( EMAX , E(IFR,IE) )
          EMIN = MIN ( EMIN , E(IFR,IE) )
          END DO
        END DO
!
      IF (EMAX.EQ.0. .AND. EMIN.EQ.0.) THEN
          EMAX   =  1.E-20
          EMIN   = -1.E-20
        END IF
!
      IF (EMAX.GT.ABS(EMIN)) THEN
          EXTR   = EMAX
        ELSE
          EXTR   = EMIN
        END IF
!
! Scaling / Normalization
!
      IF (FLSCLE) THEN
          IF (EMAX.GT.ABS(EMIN)) THEN
              FTOP   = EMAX * TOPFAC
              FSC    = FTOP / REAL(NINT(EMAX/(EMAX-EMIN)*RLINES))
            ELSE
              FTOP   = EMIN * TOPFAC
              FSC    = FTOP / REAL(NINT(EMIN/(EMAX-EMIN)*RLINES))
              FTOP   = FTOP + RLINES*FSC
              IF (ABS(FTOP).LT.0.01*FSC) FTOP = 0.
            END IF
        ELSE
          FSC    = FTOP  / RLINES
          IF (EMAX*EMIN.LT.0) FSC = 2.*FSC
          IF (EMAX.EQ.0.) FTOP = 0.
        END IF
!
! Print ID
!
      WRITE (NDS,900) PNTNME, EXTR, PRUNIT
!
! Print heading
!
      FLINE  = FTOP
      IF (MOD(NLINES,2).EQ.0) THEN
          WRITE (STRA, FMT='(E10.3)') FLINE
        ELSE
          STRA=  '          '
        END IF
!
      DO IFR=1, NFRB
        PNUM2(IFR) = '--'
        DO IE=1, NE
          IF ( NINT( (E(IFR,IE)-FLINE)/FSC ) .EQ.0) THEN
              IF (IE.LT.10) THEN
                  WRITE (STRAX,'(A1,I1)') '-', IE
                ELSE
                  WRITE (STRAX,'(I2)') IE
                END IF
              PNUM2(IFR) = STRAX
            END IF
          END DO
        END DO
!
      PNUM2(NFRB+1) = '-+'
      STRA2 = ' +'
      WRITE (NDS,910) STRA, STRA2, (PNUM2(IFR),IFR=1, NFRB+1)
!
! Print table
!
      PNUM2(NFRB+1) = ' |'
!
      DO IL = 1, NLINES-1
        FLINE  = FTOP - FSC * REAL(IL)
        IF (ABS(FLINE).LT.0.01*FSC) FLINE = 0.
        IF (MOD(NLINES-IL,2).EQ.0) THEN
            WRITE (STRA, FMT='(E10.3)') FLINE
            STRA2 =  ' +'
          ELSE
            STRA  =  '          '
            STRA2 =  ' |'
          END IF
        DO IFR=1, NFRB
          PNUM2(NFRB+1) = ' |'
          IF (ABS(FLINE).LT.0.1*FSC) THEN
              PNUM2(IFR) = '--'
              PNUM2(NFRB+1) = '-+'
              DO IE=1, NE
                IF ( NINT( (E(IFR,IE)-FLINE)/FSC ) .EQ.0) THEN
                    IF (IE.LT.10) THEN
                        WRITE (STRAX,'(A1,I1)') '-', IE
                      ELSE
                        WRITE (STRAX,'(I2)') IE
                      END IF
                    PNUM2(IFR) = STRAX
                  END IF
                END DO
            ELSE
              PNUM2(IFR) = '  '
              DO IE=1, NE
                IF ( NINT( (E(IFR,IE)-FLINE)/FSC ) .EQ.0) THEN
                    WRITE (STRAX,'(I2)') IE
                    PNUM2(IFR) = STRAX
                  END IF
                END DO
            END IF
          END DO
        WRITE (NDS,910) STRA, STRA2, (PNUM2(IFR),IFR=1, NFRB+1)
        END DO
!
! write ending
!
      FLINE  = FTOP - FSC * REAL(IL)
      IF (ABS(FLINE).LT.0.01*FSC) FLINE = 0.
      WRITE (STRA, FMT='(E10.3)') FLINE
      STRA2         = ' +'
      PNUM2(NFRB+1) = '-+'
!
      DO IFR=1, NFRB
        IF ( MOD (IFR-2,4) .EQ. 0 ) THEN
            PNUM2(IFR) = '-|'
          ELSE
            PNUM2(IFR) = '--'
          END IF
        DO IE=1, NE
          IF ( NINT( (E(IFR,IE)-FLINE)/FSC ) .EQ.0) THEN
              IF (IE.LT.10) THEN
                  WRITE (STRAX,'(A1,I1)') '-', IE
                ELSE
                  WRITE (STRAX,'(I2)') IE
                END IF
              PNUM2(IFR) = STRAX
            END IF
          END DO
        END DO
!
      WRITE (NDS,910) STRA, STRA2, (PNUM2(IFR),IFR=1, NFRB+1)
      WRITE (NDS,911) (FR(IFR)*FACFR,IFR=2,NFRB,4)
      WRITE (NDS,920)
      WRITE (NDS,921) (PRVAR(IE),IE=1,NE)
      WRITE (NDS,920)
      IF (FLSCLE) FTOP = 0.
!
      RETURN
!
! Formats
!
  900 FORMAT (/'  Location : ',A                                &
              /'  Extreme value : ',E10.3,1X,A/)
!
  910 FORMAT (A10,A2,60A2)
  911 FORMAT (10X,15F8.3)
!
  920 FORMAT (' ')
  921 FORMAT (10X,'spectra : ',10(A,'  ')/)
!/
!/ End of PRT1DM ----------------------------------------------------- /
!/
      END SUBROUTINE PRT1DM
!/ ------------------------------------------------------------------- /
      SUBROUTINE PRT2DS (NDS, NFR0, NFR, NTH, E, FR, UFR, FACSP, FSC, &
                         RRCUT, PRVAR, PRUNIT, PNTNME)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-Nov-1999 |
!/                  +-----------------------------------+
!/
!/    07-Jun-1996 : Final FORTRAN 77                    ( version 1.18 )
!/    29-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/
!  1. Purpose :
!
!     Prints a block type table of a 2-D spectrum. Input considers
!     cartesian directions, output according to meteorological
!     conventions (compass direction where waves come from).
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NDS     Int.   I   File unit number.
!       NFR0    Int.   I   Array size for freq.
!       NFR     Int.   I   Number of frequencies.
!       NTH     Int.   I   Number of frequencies.
!       E       R.A.   I   Spectral densities.
!       FR      R.A.   I   Frequencies.
!       UFR     C*(*)  I   If 'HZ', frequencies in Hz, otherwise in
!                          rad/s
!       FACSP   Real   I   Conversion factor to obtain (Hz,degr)
!                          spectrum from E
!       FSC     Real   I   Scale factor, if FSC.eq.0. automatic
!                          scaling for "compressed" block.
!       RRCUT   Real   I   Relative cut-off for printing.
!       PRVAR   C*(*)  I   Name of variable.
!       PRUNIT  C*(*)  I   Units of spectrum.
!       PNTNME  C*(*)  I   Name of location.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       ANGSTR (Internal)
!
!  5. Called by :
!
!       Any program.
!
!  6. Error messages :
!
!       None.
!
!  8. Structure :
!
!     ------------------------------------------------
!       Initializations and preparations.
!       Determine maximum of spectrum.
!       Scaling / normalization.
!       Do for normalized or non-norm. spectrum
!       ----------------------------------------------
!         Print ID
!         Print heading
!         Print table
!         Print ending
!     ------------------------------------------------
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing using STRACE.
!     !/T  Diagnostic test output.
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
      INTEGER, INTENT(IN)     :: NDS, NFR0, NFR, NTH
      REAL, INTENT(IN)        :: E(NFR0,*), FR(*), FACSP, FSC, RRCUT
      CHARACTER, INTENT(IN)   :: PRVAR*(*), PRUNIT*(*), PNTNME*(*),  &
                                 UFR*(*)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IFR, ITH, NFRB, INTANG, ITHSEC
      LOGICAL                 :: FLSCLE
      REAL                    :: FACFR, EMAX, EMIN, DTHDEG, RR, RRC
      CHARACTER               :: PNUM*5, STRA*5, STRANG*5, PNUM2*2,  &
                                 STRA2*2
      DIMENSION               :: PNUM(25), PNUM2(51)
!/
!/ ------------------------------------------------------------------- /
!/
!
! initialisations
!
      FLSCLE = .FALSE.
      IF (FSC.EQ.0.) THEN
          FLSCLE = .TRUE.
          RRC    = RRCUT * 10.
        END IF
!
      IF (UFR.EQ.'HZ') THEN
           FACFR  = 1.
         ELSE
           FACFR  = 0.159155
         END IF
!
! Maximum of spectrum
!
      EMAX   = 1.E-20
      EMIN   = 0.
!
      DO IFR=1, NFR
        DO ITH=1, NTH
          EMAX = MAX ( EMAX , E(IFR,ITH) )
          EMIN = MIN ( EMIN , E(IFR,ITH) )
          END DO
        END DO
!
      EMAX = MAX (EMAX, ABS(EMIN) )
!
      DTHDEG = 360. / REAL(NTH)
!
! Normalized spectra :  = = = = = = = = = = = = = = = = = = = = = =
!
      IF (FLSCLE) THEN
!
! Write ID
!
          WRITE (NDS,900) PNTNME, PRVAR, EMAX*FACSP, PRUNIT
!
! Write Head
!
          NFRB  = MIN (NFR,50)
          WRITE (NDS,910) (FR(IFR)*FACFR,IFR=2,NFRB,4)
!
          DO IFR=1, NFR
            IF ( MOD((IFR-2),4) .EQ. 0) THEN
                PNUM2(IFR) = '-|'
              ELSE
                PNUM2(IFR) = '--'
              END IF
            END DO
!
          PNUM2(NFRB+1) = '-+'
          WRITE (NDS,920) (PNUM2(IFR),IFR=1, NFRB+1)
!
! Write table
!
          ITHSEC = NTH + 1
!
          DO ITH= NTH, 1, -1
            INTANG = 270 - NINT (DTHDEG*REAL(ITH-1))
            IF (INTANG.LT.0) THEN
                ITHSEC = ITH
                CYCLE
              END IF
            CALL ANGSTR (INTANG, STRANG, 4, 2)
            DO IFR=1, NFRB
              RR     = E(IFR,ITH)/EMAX
              IF (E(IFR,ITH).EQ.EMAX .OR. RR.GE.1.) THEN
                  PNUM2(IFR) = ' *'
                ELSE IF (-E(IFR,ITH).EQ.EMAX .OR. RR.LE.-1.) THEN
                  PNUM2(IFR) = ' #'
                ELSE IF (ABS(RR).LT.RRC) THEN
                  PNUM2(IFR) = '  '
                ELSE IF ((RR*10.).LT.0. .AND. (RR*10.).GT.-1.) THEN
                  PNUM2(IFR) = '-0'
                ELSE
                  WRITE (STRA2, FMT='(I2)') INT (RR*10.)
                  PNUM2(IFR) = STRA2
                END IF
              END DO
            PNUM2(NFRB+1) = ' |'
            WRITE (NDS,930) STRANG, (PNUM2(IFR),IFR=1, NFRB+1)
            END DO
!
          DO ITH= NTH, ITHSEC, -1
            INTANG = 630 - NINT (DTHDEG*REAL(ITH-1))
            CALL ANGSTR (INTANG, STRANG, 4, 2)
            DO IFR=1, NFRB
              RR     = E(IFR,ITH)/EMAX
              IF (E(IFR,ITH).EQ.EMAX .OR. RR.GE.1.) THEN
                  PNUM2(IFR) = ' *'
                ELSE IF (-E(IFR,ITH).EQ.EMAX .OR. RR.LE.-1.) THEN
                  PNUM2(IFR) = ' #'
                ELSE IF (ABS(RR).LT.RRC) THEN
                  PNUM2(IFR) = '  '
                ELSE IF ((RR*10.).LT.0. .AND. (RR*10.).GT.-1.) THEN
                  PNUM2(IFR) = '-0'
                ELSE
                  WRITE (STRA2, FMT='(I2)') INT (RR*10.)
                  PNUM2(IFR) = STRA2
                END IF
              END DO
            PNUM2(NFRB+1) = ' |'
            WRITE (NDS,930) STRANG, (PNUM2(IFR),IFR=1, NFRB+1)
            END DO
!
! Write ending:
!
          PNUM2(1) = '--'
          PNUM2(2) = '-+'
          WRITE (NDS,920) (PNUM2(1),IFR=1, NFRB), PNUM2(2)
          WRITE (NDS,950)
!
! Scaled spectra :  = = = = = = = = = = = = = = = = = = = = = = = =
!
        ELSE
!
! Write ID
!
          WRITE (NDS,901) PNTNME, PRVAR, FSC, PRUNIT,           &
                          EMAX*FACSP, PRUNIT
!
! Write heading
!
          NFRB  = MIN (NFR,25)
!
          WRITE (NDS,911) (FR(IFR)*FACFR,IFR=2,NFRB,2)
          PNUM(1) = '-----'
          PNUM(2) = '--   '
!
          IF (NFRB.LT.25) THEN
              WRITE (NDS,921) (PNUM(1),IFR=1, NFRB), PNUM(2)
            ELSE
              WRITE (NDS,921) (PNUM(1),IFR=1, NFRB)
            END IF
!
!     write table :
!
          ITHSEC = NTH + 1
!
          DO ITH= NTH, 1, -1
            INTANG = 270 - NINT (DTHDEG*REAL(ITH-1))
            IF (INTANG.LT.0) THEN
                ITHSEC = ITH
                CYCLE
              END IF
            CALL ANGSTR (INTANG, STRANG, 4, 2)
            DO IFR=1, NFRB
              RR = E(IFR,ITH)
              IF (ABS(RR/EMAX).LT.RRCUT) THEN
                  PNUM(IFR) = '     '
                ELSE
                  WRITE (STRA, FMT='(I5)') NINT (RR*FACSP/FSC)
                  PNUM(IFR) = STRA
                END IF
              END DO
            WRITE (NDS,931) STRANG, (PNUM(IFR),IFR=1, NFRB)
            END DO
!
          DO ITH= NTH, ITHSEC, -1
            INTANG = 630 - NINT (DTHDEG*REAL(ITH-1))
            CALL ANGSTR (INTANG, STRANG, 4, 2)
            DO IFR=1, NFRB
              RR = E(IFR,ITH)
              IF (ABS(RR/EMAX).LT.RRCUT) THEN
                  PNUM(IFR) = '     '
                ELSE
                  WRITE (STRA, FMT='(I5)') NINT (RR*FACSP/FSC)
                  PNUM(IFR) = STRA
                END IF
              END DO
            WRITE (NDS,931) STRANG, (PNUM(IFR),IFR=1, NFRB)
            END DO
!
!     write ending :
!
          PNUM(1) = '-----'
          PNUM(2) = '--   '
          IF (NFRB.LT.25) THEN
              WRITE (NDS,921) (PNUM(1),IFR=1, NFRB), PNUM(2)
            ELSE
              WRITE (NDS,921) (PNUM(1),IFR=1, NFRB)
            END IF
          WRITE (NDS,950)
!
        END IF
!
      RETURN
!
! Formats
!
  900 FORMAT (/'  Location : ',A/                               &
               '  Spectrum : ',A,' (Normalized) ',              &
               '  Maximum value : ',E8.3,1X,A/)
  901 FORMAT (/'  Location : ',A/                               &
               '  Spectrum : ',A,'  Units : ',E8.3,1X,A,        &
               '  Maximum value : ',E8.3,1X,A/)
!
  910 FORMAT (5X,'  ang.|  frequencies (Hz) '/                  &
              5X,'  deg.|',F6.3,15F8.3)
  920 FORMAT (5X,'  ----+',60A2)
  930 FORMAT (5X,' ',A4,' |',60A2)
!
  911 FORMAT ('  ang.|  frequencies (Hz) '/                     &
              '  deg.|',12F10.3)
  921 FORMAT ('  ----|',25A5)
  931 FORMAT (' ',A4,' |',25A5)
!
  950 FORMAT (' ')
!
!/
!/ Internal subroutine ANGSTR ---------------------------------------- /
!/
      CONTAINS
!/
!/ ------------------------------------------------------------------- /
      SUBROUTINE ANGSTR (IANG, SANG, ILEN, INUM)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-Nov-1999 |
!/                  +-----------------------------------+
!/
!/    10-Mar-1992 : Final FORTRAN 77                    ( version 1.18 )
!/    29-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!
!     INPUT  : IANG --> INTEGER ANGLE (DEGREES)
!              ILEN --> STRING LENGTH
!              INUM --> <1 : ONLY FOUR MAIN DIRECTIONS
!                        1 : N,E,S,W AND NUMERICAL OUTPUT
!                        2 : EIGHT MAIN DIRECTIONS
!                       >2 : EIGHT DIRECTIONS + NUMERICAL OUTPUT
!     OUTPUT : SANG --> STRING
!
!/ ------------------------------------------------------------------- /
!/
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: IANG, ILEN, INUM
      CHARACTER, INTENT(OUT)  :: SANG*(*)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: I, J
      CHARACTER               :: SAUX*4
!/
!/ ------------------------------------------------------------------- /
!/
!     numerical :
!
      IF (INUM.EQ.1 .OR. INUM.GE.3) THEN
          WRITE (SAUX, FMT='(I4)') IANG
        ELSE
          SAUX = '    '
        END IF
!
!     string :
!
      IF (IANG.EQ.0) THEN
          SAUX = '   N'
        ELSE IF (IANG.EQ.90) THEN
          SAUX = '   E'
        ELSE IF (IANG.EQ.180) THEN
          SAUX = '   S'
        ELSE IF (IANG.EQ.270) THEN
          SAUX = '   W'
        ELSE IF (INUM.GE.2) THEN
          IF (IANG.EQ.45) THEN
              SAUX = '  NE'
            ELSE IF (IANG.EQ.135) THEN
              SAUX = '  SE'
            ELSE IF (IANG.EQ.225) THEN
              SAUX = '  SW'
            ELSE IF (IANG.EQ.315) THEN
              SAUX = '  NW'
            END IF
        END IF
!
!     Auxilary string to output :
!
      DO I=1, ILEN-4
        SANG = ' '
        END DO
      J = 0
      DO I=ILEN-3, ILEN
        J = J + 1
        SANG(I:I) = SAUX(J:J)
        END DO
      RETURN
!/
!/ End of ANGSTR ----------------------------------------------------- /
!/
      END SUBROUTINE ANGSTR
!/
!/ End of PRT2DS ----------------------------------------------------- /
!/
      END SUBROUTINE PRT2DS
!/
!/ End of module W3ARRYMD -------------------------------------------- /
!/
      END MODULE W3ARRYMD
