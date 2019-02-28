#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      PROGRAM W3TRCK
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         05-Mar-2014 |
!/                  +-----------------------------------+
!/
!/    14-Jan-1999 : Final FORTRAN 77                    ( version 1.18 )
!/    21-Jan-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    25-Jan-2001 : Flat grid version                   ( version 2.06 )
!/    20-Aug-2003 : Sequential file version             ( version 3.04 )
!/    29-Jun-2006 : Adding file name preamble.          ( version 3.09 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/    30-Oct-2009 : Implement run-time grid selection.  ( version 3.14 )
!/                  (W. E. Rogers & T. J. Campbell, NRL)
!/    05-Mar-2014 : Now calls W3SETG for pointer def.   ( version 4.18 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS. 
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Convert direct access track output file to free-format
!     readable sequential file.
!
!  2. Method :
!
!     Info read from track_o.ww3, written to track.ww3.
!
!  3. Parameters :
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3NMOD    Subr. W3GDATMD Set number of model.
!      W3NOUT    Subr. W3ODATMD Set number of model for output.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!     None, stand-alone program.
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
!       !/S    Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3GDATMD, ONLY : W3NMOD, W3SETG, FLAGLL, XFR
      USE W3ODATMD, ONLY : W3NOUT, W3SETO, FNMPRE
      USE W3SERVMD, ONLY : ITRACE, NEXTLN, EXTCDE
!/S      USE W3SERVMD, ONLY : STRACE
      USE W3TIMEMD, ONLY : STME21
!
      USE W3ODATMD, ONLY: NDSO, NDSE, NDST
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      CHARACTER*34, PARAMETER ::                                      &
                       IDTST  = 'WAVEWATCH III TRACK OUTPUT SPECTRA'
!
      INTEGER                 :: NDSI, NDSINP,                        &
                                 NDSOUT, NDSTRC, NTRACE, NK, NTH,     &
                                 NSPEC, IERR, MK, MTH,                &
                                 NREC, ILOC, ISPEC, TIME(2), TTST(2), &
                                 ILAST, NZERO, IK, ITH, IWZERO, ICH,  &
                                 IWDTH, J
!/S      INTEGER, SAVE           :: IENT   = 0
      INTEGER                 :: LINELN = 81
      REAL                    :: TH1, DTH, X, Y, DW, CX, CY, WX, WY,  &
                                 UST, AS, VALUE
      REAL                    :: SCALE  = 0.001
      REAL                    :: FACTOR
      REAL, ALLOCATABLE       :: SIG(:), DSIP(:), SPEC(:,:)
      CHARACTER               :: COMSTR*1, IDSTR*34, TSTSTR*3,        &
                                 STIME*23, STRING*81, EMPTY*81,       &
                                 PART*9, ZEROS*9, TRCKID*32
!
      DATA EMPTY(01:40) / '                                        ' /
      DATA EMPTY(41:81) / '                                         ' /
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.a Initialize data structure
!
      CALL W3NMOD ( 1, 6, 6 )
      CALL W3SETG ( 1, 6, 6 )
      CALL W3NOUT (    6, 6 )
      CALL W3SETO ( 1, 6, 6 )
!
! 1.b IO set-up.
!
      NDSI   = 10
      NDSINP = 11
      NDSOUT = 51
!
      NDSTRC =  6
      NTRACE = 10
      CALL ITRACE ( NDSTRC, NTRACE )
!
!/S      CALL STRACE ( IENT, 'W3TRCK' )
!
      WRITE (NDSO,900)
!
      J      = LEN_TRIM(FNMPRE)
      OPEN (NDSI,FILE=FNMPRE(:J)//'ww3_trck.inp',STATUS='OLD',        &
            ERR=805,IOSTAT=IERR)
      READ (NDSI,'(A)',END=806,ERR=807) COMSTR
      IF (COMSTR.EQ.' ') COMSTR = '$'
      WRITE (NDSO,901) COMSTR
!
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
      READ (NDSI,*,END=806,ERR=807) NK, NTH
      NSPEC  = NK * NTH
      WRITE (NDSO,902) NK, NTH
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 2.  Open and test input data file
!
      WRITE (NDSO,920)
!
      OPEN (NDSINP,FILE=FNMPRE(:J)//'track_o.ww3',FORM='UNFORMATTED', &
            STATUS='OLD',ERR=800,IOSTAT=IERR)
      READ (NDSINP,ERR=801,IOSTAT=IERR) IDSTR, FLAGLL, MK, MTH, XFR
!
      IF ( FLAGLL ) THEN
          FACTOR  = 1.
        ELSE
          FACTOR  = 1.E-3
        END IF
!
      IF ( IDSTR .NE. IDTST ) GOTO 810
      IF ( NK.NE.MK .OR. NTH.NE.MTH ) GOTO 811

      ALLOCATE ( SIG(MK), DSIP(MK), SPEC(MK,MTH) )
!
      READ (NDSINP,ERR=801,IOSTAT=IERR) TH1, DTH, SIG, DSIP
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 3.  Open output file and prepare
!
      WRITE (NDSO,930)
!
      OPEN (NDSOUT,FILE=FNMPRE(:J)//'track.ww3',                      &
            FORM='FORMATTED',ERR=802,IOSTAT=IERR)
!
      WRITE (NDSOUT,980,ERR=803,IOSTAT=IERR) IDSTR
      WRITE (NDSOUT,981,ERR=803,IOSTAT=IERR) MK, MTH, TH1, DTH
      WRITE (NDSOUT,982,ERR=803,IOSTAT=IERR) SIG
      WRITE (NDSOUT,983,ERR=803,IOSTAT=IERR) DSIP
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 4.  Process data
!
      ILOC    = 0
      ISPEC   = 0
      READ (NDSINP,END=444, ERR=801,IOSTAT=IERR) TTST
      BACKSPACE (NDSINP)
      WRITE (NDSO,940)
!
  400 CONTINUE
!
! 4.a Read/write basic data
!
      READ (NDSINP,END=444, ERR=801,IOSTAT=IERR) TIME, X, Y, TSTSTR,  &
                                                 TRCKID
      IF ( FLAGLL ) THEN
          WRITE (NDSOUT,984,ERR=803,IOSTAT=IERR)                      &
                 TIME, FACTOR*X, FACTOR*Y, TSTSTR, TRCKID
        ELSE
          WRITE (NDSOUT,974,ERR=803,IOSTAT=IERR)                      &
                 TIME, FACTOR*X, FACTOR*Y, TSTSTR, TRCKID
        END IF
!
      IF ( TIME(1).EQ.TTST(1) .AND. TIME(2).EQ.TTST(2) ) THEN
          ILOC = ILOC + 1
          IF ( TSTSTR .EQ. 'SEA' ) ISPEC = ISPEC + 1
        ENDIF
      IF ( TIME(1).NE.TTST(1) .OR. TIME(2).NE.TTST(2) ) THEN
          CALL STME21 ( TTST , STIME )
          WRITE (NDSO,941) STIME, ILOC, ISPEC
          ILOC    = 1
          ISPEC   = 0
          IF ( TSTSTR .EQ. 'SEA' ) ISPEC = ISPEC + 1
          TTST(1) = TIME(1)
          TTST(2) = TIME(2)
        ENDIF
!
! 4.b Check if sea point
!
      IF ( TSTSTR .NE. 'SEA' ) GOTO 400
!
! 4.c Read all data
!
      READ (NDSINP,ERR=801,IOSTAT=IERR) DW, CX, CY, WX, WY, UST, AS,  &
                                        SPEC
      IF ( UST .LT. 0. ) UST = -1.0
!
! 4.d Write the basic stuff
!
      WRITE (NDSOUT,985,ERR=803,IOSTAT=IERR)                          &
            DW, CX, CY, WX, WY, UST, AS, SCALE
!
! 4.e Start of integer packing
!
      STRING = EMPTY
      ILAST  = 0
      NZERO  = 0
!
! 4.e.1 Loop over spectrum
!
      DO IK=1, NK
        DO ITH=1, NTH
          VALUE  = MAX ( 0.1 , 1.1*SPEC(IK,ITH)/SCALE )
          IWDTH  = 2 + MAX( 0 , INT( ALOG10(VALUE) ) )
!
! 4.e.2 Put value in string and test overflow
!
          IF ( IWDTH .GT. 9 ) THEN
              IWDTH   = 9
              PART    = ' 99999999'
            ELSE
              WRITE (PART,987) NINT(SPEC(IK,ITH)/SCALE)
              IF ( PART(11-IWDTH:11-IWDTH) .EQ. ' ' )                 &
                   IWDTH   = IWDTH - 1
            ENDIF
!
! 4.e.3 It's a zero, wait with writing
!
          IF ( PART(8:9) .EQ. ' 0' ) THEN
              NZERO  = NZERO + 1
            ELSE
!
! 4.e.4 It's not a zero, write unwritten zeros
!
              IF ( NZERO .NE. 0 ) THEN
                  IF ( NZERO .EQ. 1 ) THEN
                      ZEROS  = '        0'
                      IWZERO = 2
                    ELSE
                      WRITE (ZEROS,'(I7,A2)') NZERO, '*0'
                      IWZERO = 4
                      DO
                        ICH    = 10 - IWZERO
                        IF ( ZEROS(ICH:ICH) .NE. ' ' ) THEN
                            IWZERO = IWZERO + 1
                          ELSE
                            EXIT
                          ENDIF
                        END DO
                    ENDIF
                  IF ( ILAST+IWZERO .GT. LINELN ) THEN
                      WRITE (NDSOUT,986,ERR=803,IOSTAT=IERR)          &
                                   STRING(2:ILAST)
                      STRING = EMPTY
                      ILAST  = 0
                    ENDIF
                  STRING(ILAST+1:ILAST+IWZERO) =                      &
                                   ZEROS(10-IWZERO:9)
                  ILAST  = ILAST + IWZERO
                  NZERO  = 0
                ENDIF
!
! 4.e.5 It's not a zero, put in string
!
              IF ( ILAST+IWDTH .GT. LINELN ) THEN
                  WRITE (NDSOUT,986,ERR=803,IOSTAT=IERR)              &
                               STRING(2:ILAST)
                  STRING = EMPTY
                  ILAST  = 0
                ENDIF
!
              STRING(ILAST+1:ILAST+IWDTH) = PART(10-IWDTH:9)
              ILAST  = ILAST + IWDTH
!
            ENDIF
!
          END DO
        END DO
!
! ..... End of loop over spectrum (4.e.1)
!
! 4.e.6 Write trailing zeros
!
      IF ( NZERO .NE. 0 ) THEN
          IF ( NZERO .EQ. 1 ) THEN
              ZEROS  = '        0'
              IWZERO = 2
            ELSE
              WRITE (ZEROS,'(I7,A2)') NZERO, '*0'
              IWZERO = 4
              DO
                ICH    = 10 - IWZERO
                IF ( ZEROS(ICH:ICH) .NE. ' ' ) THEN
                    IWZERO = IWZERO + 1
                  ELSE
                    EXIT
                  ENDIF
                END DO
            ENDIF
          IF ( ILAST+IWZERO .GT. LINELN ) THEN
              WRITE (NDSOUT,986,ERR=803,IOSTAT=IERR)                  &
                           STRING(2:ILAST)
              STRING = EMPTY
              ILAST  = 0
            ENDIF
          STRING(ILAST+1:ILAST+IWZERO) = ZEROS(10-IWZERO:9)
          ILAST  = ILAST + IWZERO
          NZERO  = 0
        ENDIF
!
! 4.e.7 Write last line
!
      IF ( ILAST .NE. 0 ) THEN
          WRITE (NDSOUT,986,ERR=803,IOSTAT=IERR) STRING(2:ILAST)
        ENDIF
!
! ... Loop back to top
!
      GOTO 400
!
! 4.f All data done, write last batch info
!
  444 CONTINUE
!
      CALL STME21 ( TTST , STIME )
      WRITE (NDSO,941) STIME, ILOC, ISPEC
!
      GOTO 888
!
! Escape locations read errors :
!
  800 CONTINUE
      WRITE (NDSE,1000) IERR
      CALL EXTCDE ( 1 )
!
  801 CONTINUE
      WRITE (NDSE,1001) IERR
      CALL EXTCDE ( 2 )
!
  802 CONTINUE
      WRITE (NDSE,1002) IERR
      CALL EXTCDE ( 3 )
!
  803 CONTINUE
      WRITE (NDSE,1003) IERR
      CALL EXTCDE ( 4 )
!
  805 CONTINUE
      WRITE (NDSE,1004) IERR
      CALL EXTCDE ( 5 )
!
  806 CONTINUE
      WRITE (NDSE,1005) IERR
      CALL EXTCDE ( 6 )
!
  807 CONTINUE
      WRITE (NDSE,1006) IERR
      CALL EXTCDE ( 7 )
!
  810 CONTINUE
      WRITE (NDSE,1010) IDSTR, IDTST
      CALL EXTCDE ( 5 )
!
  811 CONTINUE
      WRITE (NDSE,1011) MK, MTH, NK, NTH
      CALL EXTCDE ( 6 )
!
  888 CONTINUE
!
      WRITE (NDSO,999)
!
! Formats
!
  900 FORMAT (/15X,'    *** WAVEWATCH III Track output post.***    '/ &
               15X,'==============================================='/)
  901 FORMAT ( '  Comment character is ''',A,''''/)
  902 FORMAT ( '  Spectral grid size is ',I3,' by ',I3//              &
                   '  Opening files : '/                              &
                   ' -----------------------------------------------')
  920 FORMAT ( '     Input file ...')
  930 FORMAT ( '     Output file ...')
  940 FORMAT (/'  Processing data : '/                                &
               ' -----------------------------------------------')
  941 FORMAT ( '     ',A,' :',I6,' points and',I6,'  spectra.')
!
  980 FORMAT (A)
  981 FORMAT (2I6,2E13.5)
  982 FORMAT (7E11.4)
  983 FORMAT (7E11.4)
  984 FORMAT (I8.8,I7.6,2F9.3,2X,A3,2X,A32)
  974 FORMAT (I8.8,I7.6,2(F9.2,'E3'),2X,A3,2X,A32)
  985 FORMAT (F8.1,2F6.2,2F8.2,f9.5,f7.2,E12.5)
  986 FORMAT (A)
  987 FORMAT (I9)
!
  999 FORMAT (/'  End of program '/                                   &
               ' ========================================='/          &
               '         WAVEWATCH III Track output '/)
!
 1000 FORMAT (/' *** WAVEWATCH III ERROR IN W3TRCK : '/               &
               '     ERROR IN OPENING INPUT DATA FILE'/               &
               '     IOSTAT =',I5/)
!
 1001 FORMAT (/' *** WAVEWATCH III ERROR IN W3TRCK : '/               &
               '     ERROR IN READING FROM INPUT DATA FILE'/          &
               '     IOSTAT =',I5/)
!
 1002 FORMAT (/' *** WAVEWATCH III ERROR IN W3TRCK : '/               &
               '     ERROR IN OPENING OUTPUT DATA FILE'/              &
               '     IOSTAT =',I5/)
!
 1003 FORMAT (/' *** WAVEWATCH III ERROR IN W3TRCK : '/               &
               '     ERROR IN WRITING TO OUTPUT FILE'/                &
               '     IOSTAT =',I5/)
!
 1004 FORMAT (/' *** WAVEWATCH III ERROR IN W3TRCK : '/               &
               '     ERROR IN OPENING INPUT FILE'/                    &
               '     IOSTAT =',I5/)
!
 1005 FORMAT (/' *** WAVEWATCH III ERROR IN W3TRCK : '/               &
               '     ERROR IN READING FROM INPUT FILE'/               &
               '     IOSTAT =',I5/)
!
 1006 FORMAT (/' *** WAVEWATCH III ERROR IN W3TRCK : '/               &
               '     ERROR IN OPENING OUTPUT FILE'/                   &
               '     IOSTAT =',I5/)
!
 1010 FORMAT (/' *** WAVEWATCH III ERROR IN W3TRCK : '/               &
               '     UNEXPECTED ID STRING IN INPUT : ',A/             &
               '                         SHOULD BE : ',A/)
!
 1011 FORMAT (/' *** WAVEWATCH III ERROR IN W3TRCK : '/               &
               '     UNEXPECTED SPECTRAL DIMENSIONS : ',2I4/          &
               '                          SHOULD BE : ',2I4/)
!/
!/ End of W3TRCK ----------------------------------------------------- /
!/
      END PROGRAM W3TRCK
