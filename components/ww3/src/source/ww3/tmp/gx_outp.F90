#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      PROGRAM GXOUTP
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |            J.H. Alves             |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         27-Aug-2015 |
!/                  +-----------------------------------+
!/
!/    30-Jun-1999 : Final FORTRAN 77                    ( version 1.18 )
!/    24-Jan-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    14-Feb-2000 : Exact nonlinear interactions        ( version 2.01 )
!/    25-Jan-2001 : Cartesian grid version              ( version 2.06 )
!/    02-Feb-2001 : Xnl version 3.0                     ( version 2.07 )
!/    13-Nov-2002 : Add stress vector.                  ( version 3.00 )
!/    27-Nov-2002 : First version of VDIA and MDIA.     ( version 3.01 )
!/    01-Aug-2003 : Fix format for SH output points.    ( version 3.03 )
!/    24-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    23-Jun-2006 : Linear input added.                 ( version 3.09 )
!/    29-Jun-2006 : Adding file name preamble.          ( version 3.09 )
!/    03-Jul-2006 : Separate flux modules.              ( version 3.09 )
!/    25-Jul-2006 : Grid ID for each point.             ( version 3.10 )
!/    25-Apr-2007 : EMEAN in W3SPR2 par list.           ( version 3.11 )
!/    09-Oct-2007 : WAM 4+ Sin and Sds added.           ( version 3.13 )
!/                  (F. Ardhuin)
!/    09-Oct-2007 : Experimental Sbs (BS1) added.       ( version 3.13 )
!/                  (F. Ardhuin)
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/    30-Oct-2009 : Implement run-time grid selection.  ( version 3.14 )
!/                  (W. E. Rogers & T. J. Campbell, NRL)
!/    30-Aug-2010 : Adding ST4                          ( version 3.14 )
!/    20-Apr-2010 : Fix initialization of USTAR.      ( version 3.14.1 )
!/    23-Aug-2012 : Adding movable bed friction BT4     ( version 4.07 )
!/    16-Jul-2012 : Move GMD (SNL3) and nonlinear filter (SNLS)
!/                  from 3.15 (HLT).                    ( version 4.08 )
!/    26-Dec-2012 : Modified obsolete declarations.     ( version 4.11 )
!/    27-Aug-2015 : Sice add as additional output       ( version 5.10 )
!/                  (in source terms)
!/
!/    Copyright 2009-2012 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Post-processing of point output for GrADS post-processing.
!
!  2. Method :
!
!     In order to be able to plot spectra and source terms as
!     fields, spectral data is written as if it is fields data.
!     The spectral direction becomes the longitude, 90.-FREQ
!     become the latitude. This way, polar plots can be made
!     using the GrADS 'NPS' map option. The level or z coordinate
!     is used to store spectra and source terms for separate
!     output points. The name of the output point is stored in
!     the control file as the 'description' of the field.
!     Also written is a separate file with mean input and wave
!     parameters. This file contains per level and per time a
!     single line containing :
!
!        Station ID, Longitude, Latitude, Depth, , Wind speed.
!           U and V components, Air-Sea Temperature difference,
!           Current velocity, U and V components, Significant
!           wave height.
!
!     The files generated are :
!
!       ww3.spec.ctl     GrADS control file.
!       ww3.spec.grads   GrADS data file.
!       ww3.mean.grads   File with additional input and wave
!                        parameters.
!
!     The first direction set to 90 degr. Grads NPS plot should
!     therefore have 'set lon -180 180' for oceanographic directional
!     convention.
!
!     Examples of using the three files can be found in spec.gs and
!     source.gs.
!
!  3. Parameters :
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3NMOD    Subr. W3GDATMD Set number of model.
!      W3SETG    Subr.   Id.    Point to selected model.
!      W3NDAT    Subr. W3WDATMD Set number of model for wave data.
!      W3SETW    Subr.   Id.    Point to selected model for wave data.
!      W3NAUX    Subr. W3ADATMD Set number of model for aux data.
!      W3SETA    Subr.   Id.    Point to selected model for aux data.
!      W3NOUT    Subr. W3ODATMD Set number of model for output.
!      W3SETO    Subr.   Id.    Point to selected model for output.
!      ITRACE    Subr. W3SERVMD Subroutine tracing initialization.
!      STRACE    Subr.   Id.    Subroutine tracing.
!      NEXTLN    Subr.   Id.    Get next line from input filw
!      EXTCDE    Subr.   Id.    Abort program as graceful as possible.
!      STME21    Subr. W3TIMEMD Convert time to string.
!      TICK21    Subr.   Id.    Advance time.
!      DSEC21    Func.   Id.    Difference between times.
!      W3IOGR    Subr. W3IOGRMD Reading/writing model definition file.
!      W3IOPO    Subr. W3IOPOMD Reading/writing raw point output file.
!      GXEXPO    Subr. Internal Execute point output.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!     None, stand-alone program.
!
!  6. Error messages :
!
!     Checks on input, checks in W3IOxx.
!     Check on grid type.
!
!  7. Remarks :
!
!     - Curvilinear grids currently not supported.
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
      USE CONSTANTS
!/
!     USE W3GDATMD, ONLY: W3NMOD, W3SETG
      USE W3WDATMD, ONLY: W3SETW, W3NDAT
      USE W3ADATMD, ONLY: W3SETA, W3NAUX
      USE W3ODATMD, ONLY: W3SETO, W3NOUT
      USE W3IOGRMD, ONLY: W3IOGR
      USE W3IOPOMD, ONLY: W3IOPO
      USE W3SERVMD, ONLY : ITRACE, NEXTLN, EXTCDE
      USE W3TIMEMD, ONLY: STME21, TICK21, DSEC21
!/
      USE W3GDATMD
      USE W3WDATMD, ONLY: TIME
      USE W3ODATMD, ONLY: NDSE, NDST, NDSO, NOPTS, PTLOC, PTNME,      &
                          DPO, WAO, WDO, ASO, CAO, CDO, SPCO, FNMPRE, &
                          GRDID, ICEO, ICEHO, ICEFO
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: NDSI, NDSM, NDSOP, NDSGRD, NDSPNT,   &
                                 NDSCGR, NDSTRC, NTRACE, IERR,        &
                                 IOTEST, I, TOUT(2), NOUT, TDUM(2),   &
                                 NREQ, IPOINT, NLEV, IOUT, TIME0(2),  &
                                 IH0, IM0, ID0, IID, IJ0, IINC, IK,   &
                                 IREQ, TIMEN(2), J
      REAL                    :: DTREQ, DTEST
      REAL                    :: UNDEFP = -99.E20
      REAL                    :: FACT
      LOGICAL                 :: FLSRCE(7)
      LOGICAL, ALLOCATABLE    :: FLREQ(:)
      CHARACTER               :: COMSTR*1, IDTIME*23, IDDDAY*11,      &
                                 CINC*2
      CHARACTER(LEN=3)        :: MNTH(12)
      CHARACTER(LEN=25)       :: IDSRCE(7)
!/
!/ ------------------------------------------------------------------- /
!/
      DATA IDSRCE / 'Spectrum                 ' ,                     &
                    'Wind-wave interactions   ' ,                     &
                    'Nonlinear interactions   ' ,                     &
                    'Dissipation              ' ,                     &
                    'Wave-bottom interactions ' ,                     &
                    'Wave-ice interactions ' ,                     &
                    'Sum of selected sources  ' /
      DATA FLSRCE / .FALSE. , .FALSE. , .FALSE. ,                     &
                    .FALSE. , .FALSE. , .FALSE., .FALSE. /
      DATA TIME0  / -1, 0 /
      DATA MNTH   / 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',         &
                    'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 1.  IO set-up.
!
      CALL W3NMOD ( 1, 6, 6 )
      CALL W3SETG ( 1, 6, 6 )
      CALL W3NDAT (    6, 6 )
      CALL W3SETW ( 1, 6, 6 )
      CALL W3NAUX (    6, 6 )
      CALL W3SETA ( 1, 6, 6 )
      CALL W3NOUT (    6, 6 )
      CALL W3SETO ( 1, 6, 6 )
!
      NDSI   = 10
      NDSM   = 20
      NDSOP  = 20
      NDSGRD = 30
      NDSPNT = 31
      NDSCGR = 32
!
      NDSTRC =  6
      NTRACE =  0
!
      WRITE (NDSO,900)
!
      CALL ITRACE ( NDSTRC, NTRACE )
!
      J      = LEN_TRIM(FNMPRE)
      OPEN (NDSI,FILE=FNMPRE(:J)//'gx_outp.inp',STATUS='OLD',         &
            ERR=800,IOSTAT=IERR)
      READ (NDSI,'(A)',END=801,ERR=802) COMSTR
      IF (COMSTR.EQ.' ') COMSTR = '$'
      WRITE (NDSO,901) COMSTR
!
      OPEN (NDSGRD,FILE=FNMPRE(:J)//'ww3.spec.grads',                 &
            FORM='UNFORMATTED', CONVERT='big_endian')
      OPEN (NDSPNT,FILE=FNMPRE(:J)//'ww3.mean.grads',FORM='FORMATTED')
      OPEN (NDSCGR,FILE=FNMPRE(:J)//'ww3.spec.ctl',FORM='FORMATTED')
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 2.  Read model definition file.
!
      CALL W3IOGR ( 'READ', NDSM )
      WRITE (NDSO,920) GNAME
      IF ( FLAGLL ) THEN
          FACT   = 1.
        ELSE
          FACT   = 1.E-3
        END IF
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 3.  Read general data and first fields from file
!
      CALL W3IOPO ( 'READ', NDSOP, IOTEST )
      ALLOCATE ( FLREQ(NOPTS) )
!
      WRITE (NDSO,930)
      DO I=1, NOPTS
        IF ( FLAGLL ) THEN
            WRITE (NDSO,931) PTNME(I), FACT*PTLOC(1,I), FACT*PTLOC(2,I)
          ELSE
            WRITE (NDSO,932) PTNME(I), FACT*PTLOC(1,I), FACT*PTLOC(2,I)
          END IF
        END DO
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 4.  Read requests from input file.
!     Output times
!
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
      READ (NDSI,*,END=801,ERR=802) TOUT, DTREQ, NOUT
      DTREQ  = MAX ( 0. , DTREQ )
      IF ( DTREQ.EQ.0 ) NOUT = 1
      NOUT   = MAX ( 1 , NOUT )
!
      CALL STME21 ( TOUT , IDTIME )
      WRITE (NDSO,940) IDTIME
!
      TDUM = 0
      CALL TICK21 ( TDUM , DTREQ )
      CALL STME21 ( TDUM , IDTIME )
      IF ( DTREQ .GE. 86400. ) THEN
          WRITE (IDDDAY,'(I10,1X)') INT(DTREQ/86400.)
        ELSE
          IDDDAY = '           '
        END IF
      IDTIME(1:11) = IDDDAY
      IDTIME(21:23) = '   '
      WRITE (NDSO,941) IDTIME, NOUT
!
! ... Output points
!
      FLREQ = .FALSE.
      NREQ   = 0
!
      DO
        CALL NEXTLN ( COMSTR , NDSI , NDSE )
        READ (NDSI,*,END=801,ERR=802) IPOINT
        IF ( IPOINT .GT. 0 ) THEN
            IF ( IPOINT .LE. NOPTS ) THEN
                IF ( .NOT. FLREQ(IPOINT) ) NREQ   = NREQ + 1
                FLREQ(IPOINT) = .TRUE.
              END IF
          ELSE
            EXIT
          END IF
        END DO
!
! ... Output of output points
!
      WRITE (NDSO,950) NREQ
      DO I=1, NOPTS
        IF (FLREQ(I)) THEN
            IF ( FLAGLL ) THEN
                WRITE (NDSO,951) PTNME(I), FACT*PTLOC(1,I), &
                                           FACT*PTLOC(2,I)
              ELSE
                WRITE (NDSO,956) PTNME(I), FACT*PTLOC(1,I), &
                                           FACT*PTLOC(2,I)
              END IF
          END IF
        END DO
!
! ... Output of output points
!
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
      READ (NDSI,*,END=801,ERR=802) FLSRCE
      WRITE (NDSO,952)
      NLEV   = 0
      DO I=1, 7
        IF ( FLSRCE(I) ) THEN
            WRITE (NDST,953) IDSRCE(I)
            NLEV   = NLEV + 1
          END IF
        END DO
!
      WRITE (NDSO,955)
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 5.  Time management.
!
      IOUT   = 0
!
      DO
        DTEST  = DSEC21 ( TIME , TOUT )
        IF ( DTEST .GT. 0. ) THEN
            CALL W3IOPO ( 'READ', NDSOP, IOTEST )
            IF ( IOTEST .EQ. -1 ) THEN
                WRITE (NDSO,998)
                EXIT
              END IF
            CYCLE
          END IF
        IF ( DTEST .LT. 0. ) THEN
            CALL TICK21 ( TOUT , DTREQ )
            CYCLE
          END IF
!
        IOUT   = IOUT + 1
        CALL STME21 ( TOUT , IDTIME )
!
        CALL GXEXPO
        TIMEN  = TOUT
!
        IF ( TIME0(1) .EQ. -1 ) TIME0 = TIME
!
        CALL TICK21 ( TOUT , DTREQ )
        IF ( IOUT .GE. NOUT ) EXIT
        END DO
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 6.  Close data file and write control file
! 6.a Close data sets
!
      WRITE (NDSO,960)
!
      WRITE (NDSO,961)
      CLOSE (NDSGRD)
      CLOSE (NDSPNT)
!
      WRITE (NDSO,962)
!
! 6.b Set up timing info
!
      IH0    = TIME0(2)/10000
      IM0    = MOD(TIME0(2)/100,100)
      ID0    = MOD(TIME0(1),100)
      IID    = MOD(TIME0(1)/100,100)
      IJ0    = TIME0(1)/10000
!
      IF ( IOUT .GT. 1 ) DTREQ  = DSEC21 ( TIME0, TIMEN ) / REAL(IOUT-1)
      IF ( IOUT .EQ. 1 ) DTREQ  = 3600.
      IF ( DTREQ .GT. 3599. ) THEN
          CINC   = 'HR'
          IINC   = NINT(DTREQ/3600.)
          IF ( MOD(NINT(DTREQ),3600) .NE. 0 ) GOTO 820
        ELSE
          CINC   = 'MN'
          IINC   = NINT(DTREQ/60.)
        END IF
!
      WRITE (NDSO,963) IOUT, IH0, IM0, ID0, MNTH(IID), IJ0, IINC, CINC
!
! 6.c Write control file for spectral data
!
      WRITE (NDSO,964)
!
      WRITE (NDSCGR,970) UNDEFP, NTH, 90.+TH(1)*RADE, DTH*RADE,       &
             NK, (90.-TPIINV*SIG(IK),IK=NK,MAX(1,NK-4),-1)
      WRITE (NDSCGR,971) (90.-TPIINV*SIG(IK),IK=NK-5,1,-1)
      WRITE (NDSCGR,972) NLEV, 1., 1.,                                &
                         IOUT, IH0, IM0, ID0, MNTH(IID), IJ0,         &
                         IINC, CINC, NREQ
!
      IREQ    = 0
      DO I=1, NOPTS
        IF ( FLREQ(I) ) THEN
            IREQ   = IREQ + 1
            WRITE (NDSCGR,973) IREQ, NLEV, 99, PTNME(I)
          END IF
        END DO
!
      WRITE (NDSCGR,974)
!
      GOTO 888
!
! Escape locations read errors :
!
  800 CONTINUE
      WRITE (NDSE,1000) IERR
      CALL EXTCDE ( 10 )
!
  801 CONTINUE
      WRITE (NDSE,1001)
      CALL EXTCDE ( 11 )
!
  802 CONTINUE
      WRITE (NDSE,1002) IERR
      CALL EXTCDE ( 12 )
!
  820 CONTINUE
      WRITE (NDSE,1020) DTREQ
      CALL EXTCDE ( 20 )
!
  821 CONTINUE
      WRITE (NDSE,1021)
      CALL EXTCDE ( 21 )
!
  888 CONTINUE
!
      WRITE (NDSO,999)
!
! Formats
!
  900 FORMAT (/12X,'    *** WAVEWATCH III GrADS point output post.***    '/ &
               12X,'====================================================='/)
  901 FORMAT ( '  Comment character is ''',A,''''/)
!
  920 FORMAT ( '  Grid name : ',A/)
!
  930 FORMAT ( '  Points in file : '/                                 &
               ' ------------------------------------')
 
  931 FORMAT ( '      ',A,2F10.2)
 
  932 FORMAT ( '      ',A,2(F8.1,'E3'))
!
  940 FORMAT (/'  Output time data : '/                               &
               ' --------------------------------------------------'/ &
               '      First time         : ',A)
  941 FORMAT ( '      Interval           : ',A/                       &
               '      Number of requests : ',I4)
!
  950 FORMAT (/'  Requested output for',I3,' points : '/              &
               ' --------------------------------------------------')
 
  951 FORMAT ( '      ',A,2F10.2)
 
  956 FORMAT ( '      ',A,2(F8.1,'E3'))
 
  952 FORMAT (/'  Requested output fields :'/                         &
               ' --------------------------------------------------')
  953 FORMAT ( '      ',A)
  955 FORMAT (/'  Output times :'/                                    &
               ' --------------------------------------------------')
!
  960 FORMAT (//'  Final file management '/                           &
               ' -----------------------------------------------------')
  961 FORMAT ( '      Closing file ww3.spec.grads'/                   &
               '      Closing file ww3.mean.grads')
  962 FORMAT ( '      Preparing control files :')
  963 FORMAT ( '         Number of times : ',I6/                      &
           '         Initial time ID : ',I2.2,':',I2.2,'Z',I2.2,A3,I4/ &
           '         Time step ID    : ',I2,A2)
  964 FORMAT ( '      Writing ww3.spec.ctl'/)
!
  970 FORMAT ('DSET      ww3.spec.grads'/                             &
              'TITLE     WAVEWATCH III spectra and source terms'/     &
              'OPTIONS   sequential'/                                 &
              'OPTIONS   big_endian'/                                 &
              'UNDEF    ',E10.2/                                      &
              'XDEF     ',I4,'  LINEAR ',2F8.2/                       &
              'YDEF     ',I4,'  LEVELS ',5F8.4)
  971 FORMAT (22X,5F8.4)
  972 FORMAT ('ZDEF     ',I4,'  LINEAR ',2F8.2/                       &
              'TDEF     ',I4,'  LINEAR ',I6.2,':',I2.2,'Z',I2.2,A3,I4, &
               2x,I2,A2/                                              &
              'VARS     ',I4)
  973 FORMAT ('LOC',I3.3,2I4,2X,A)
  974 FORMAT ('ENDVARS')
!
  998 FORMAT (/'      End of file reached '/)
!
  999 FORMAT (/'  End of program '/                                   &
               ' ========================================='/          &
               '         WAVEWATCH III GrADS point output '/)
!
 1000 FORMAT (/' *** WAVEWATCH III ERROR IN GXOUTP : '/               &
               '     ERROR IN OPENING INPUT FILE'/                    &
               '     IOSTAT =',I5/)
!
 1001 FORMAT (/' *** WAVEWATCH III ERROR IN GXOUTP : '/               &
               '     PREMATURE END OF INPUT FILE'/)
!
 1002 FORMAT (/' *** WAVEWATCH III ERROR IN GXOUTP : '/               &
               '     ERROR IN READING FROM INPUT FILE'/               &
               '     IOSTAT =',I5/)
!
 1020 FORMAT (/' *** WAVEWATCH III ERROR IN GXOUTF : '/               &
               '     FIELD INCREMENT > 1HR BUT NOT MULTIPLE',F10.0/)
!
 1021 FORMAT (/' *** WAVEWATCH III ERROR IN GXOUTF : '/               &
               '     UPDATE PARS IN LOOP 610 !!!'/)
!/
!/ Internal subroutine GXEXPO ---------------------------------------- /
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE GXEXPO
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         16-Jul-2012 |
!/                  +-----------------------------------+
!/
!/    30-Jun-1999 : Final FORTRAN 77                    ( version 1.18 )
!/    24-Jan-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/                  Massive changes to logistics
!/    25-Jan-2001 : Cartesian grid version              ( version 2.06 )
!/    02-Feb-2001 : Xnl version 5                       ( version 2.07 )
!/    01-Aug-2003 : Fix format for SH output points.    ( version 3.03 )
!/    24-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    23-Jun-2006 : Linear input added.                 ( version 3.09 )
!/    03-Jul-2006 : Separate flux modules.              ( version 3.09 )
!/    25-Jul-2006 : Grid ID for each point.             ( version 3.10 )
!/    25-Apr-2007 : EMEAN in W3SPR2 par list.           ( version 3.11 )
!/    09-Oct-2007 : WAM 4+ Sin and Sds added.           ( version 3.13 )
!/                  (F. Ardhuin)
!/    09-Oct-2007 : Experimental Sbs (BS1) added.       ( version 3.13 )
!/                  (F. Ardhuin)
!/    16-Jul-2012 : Move GMD (SNL3) and nonlinear filter (SNLS)
!/                  from 3.15 (HLT).                    ( version 4.08 )
!/
!  1. Purpose :
!
!     Perform actual point output.
!
!  3. Parameters :
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SPRn    Subr. W3SRCnMD Mean wave parameters for use in
!                               source terms.
!      W3FLXn    Subr. W3FLXnMD Flux/stress computation.
!      W3SLNn    Subr. W3SLNnMD Linear input.
!      W3SINn    Subr. W3SRCnMD Input source term.
!      W3SDSn    Subr. W3SRCnMD Whitecapping source term
!      W3SNLn    Subr. W3SNLnMD Nonlinear interactions.
!      W3SBTn    Subr. W3SBTnMD Bottom friction source term.
!      W3SDBn    Subr. W3SBTnMD Depth induced breaking source term.
!      W3STRn    Subr. W3STRnMD Triad interaction source term.
!      W3SBSn    Subr. W3SBSnMD Bottom scattering source term.
!      W3SXXn    Subr. W3SXXnMD Unclassified source term.
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!      STME21    Subr. W3TIMEMD Convert time to string.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!     Program in which it is contained.
!
!  6. Error messages :
!
!     None.
!
!  7. Remarks :
!
!     - Spectra are relative frequency energy spectra.
!     - Note that arrays CX and CY of the main program now contain
!       the absolute current speed and direction respectively.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!       !/S      Enable subroutine tracing.
!       !/T      Enable test output.
!
!       !/FLXx   Flux/stress computation.
!       !/LNx    Linear input package
!       !/STx    Source term package
!       !/NLx    Nonlinear interaction package
!       !/BTx    Bottom friction package
!       !/ICx    Ice source term package
!       !/DBx    Depth-induced breaking package
!       !/TRx    Triad interaction package
!       !/BSx    Bottom scattering package
!       !/XXx    Arbitrary adittional source term package
!
!       !/STAB2  Stability correction for !/ST2
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3SLN1MD
      USE W3SRC3MD
      USE W3SNL1MD
      USE W3SBT1MD
      USE W3SDB1MD
!/
      USE W3DISPMD, ONLY: NAR1D, DFAC, N1MAX, ECG1, EWN1, DSIE
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: J, I1, I2, IK, ITH, ISPEC, IKM, IKL, &
                                 IKH, ITT, IX, IY, ISEA
      REAL                    :: XL, XH, XL2, XH2, DEPTH, SQRTH, UDIR,&
                                 UDIRR, UABS, CDIR, SIX, R1, R2, ET,  &
                                 EWN, ETR, ETX, ETY, EBND, EBX, EBY,  &
                                 HSIG, WLEN, TMEAN, THMEAN, THSPRD,   &
                                 EMAX, EL, EH, DENOM, FP, THP, SPP,   &
                                 FACTOR, CD, USTAR, FHIGH, ZWND, ICE, &
                                 USTD, Z0, CHARN, EMEAN, FMEAN, WNMEAN,&
                                 ICETHICK, ICECON
      REAL                    :: FMEANS, FMEANWS, TAUWX, TAUWY, AMAX, &
                                 TAUWNX, TAUWNY
      REAL                    :: HSMIN = 0.05
      REAL                    :: WN(NK), CG(NK), E(NK,NTH), E1(NK),   &
                                 APM(NK), THBND(NK), SPBND(NK),       &
                                 A(NTH,NK), WN2(NTH,NK),WN_R(NK),     &
                                 ALPHA_LIU(NK), CG_ICE(NK)
      REAL                    :: DIA(NTH,NK), SWI(NK,NTH), SNL(NK,NTH),&
                                 SDS(NK,NTH), SBT(NK,NTH), SIS(NK,NTH),&
                                 STT(NK,NTH), DIA2(NK,NTH)
      REAL                    :: XLN(NTH,NK), XWI(NTH,NK), XNL(NTH,NK),&
                                 XTR(NTH,NK), XDS(NTH,NK), XDB(NTH,NK),&
                                 XBT(NTH,NK), XBS(NTH,NK), XXX(NTH,NK),&
                                 XWL(NTH,NK), XIS(NTH,NK)
      LOGICAL                :: LLWS(NTH,NK)
      CHARACTER               :: DTME21*23
!/
!/ ------------------------------------------------------------------- /
!/
!
      XL     = 1./XFR - 1.
      XH     =  XFR - 1.
      XL2    = XL**2
      XH2    = XH**2
      ICE = 0.
!
      XLN = 0.
      XWI = 0.
      XNL = 0.
      XTR = 0.
      XDS = 0.
      XDB = 0.
      XBT = 0.
      XBS = 0.
      XWL = 0.
      XIS = 0.
      XXX = 0.
!
!     Output of time
!
      CALL STME21 ( TIME , DTME21 )
      WRITE (NDSO,905) DTME21
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Loop over output points.
!
      DO J=1, NOPTS
        IF ( FLREQ(J) ) THEN
!
! 2. Calculate grid parameters using and inlined version of WAVNU1.
!
          DEPTH  = MAX ( DMIN, DPO(J) )
          SQRTH  = SQRT ( DEPTH )
          UDIR   = MOD ( 270. - WDO(J)*RADE , 360. )
          UDIRR  = WDO(J)
          UABS   = MAX ( 0.001 , WAO(J) )
          CDIR   = MOD ( 270. - CDO(J)*RADE , 360. )
                ICETHICK = MAX (0., ICEHO(J))
                ICECON = MAX (0., ICEO(J))
!
            DO IK=1, NK
              SIX    = SIG(IK) * SQRTH
              I1     = INT(SIX/DSIE)
              IF (I1.LE.N1MAX) THEN
                  I2 = I1 + 1
                  R1 = SIX/DSIE - REAL(I1)
                  R2 = 1. - R1
                  WN(IK) = ( R2*EWN1(I1) + R1*EWN1(I2) ) / DEPTH
                  CG(IK) = ( R2*ECG1(I1) + R1*ECG1(I2) ) * SQRTH
                ELSE
                  WN(IK) = SIG(IK)*SIG(IK)/GRAV
                  CG(IK) = 0.5 * GRAV / SIG(IK)
                END IF
!
              END DO
!
! 3.  Prepare spectra etc.
! 3.a Mean wave parameters.
!
            ET     = 0.
            EWN    = 0.
            ETR    = 0.
            ETX    = 0.
            ETY    = 0.
            DO IK=1, NK
              EBND   = 0.
              EBX    = 0.
              EBY    = 0.
              DO ITH=1, NTH
                ISPEC  = ITH + (IK-1)*NTH
                E(IK,ITH) = SPCO(ISPEC,J)
                EBND   = EBND + SPCO(ISPEC,J)
                EBX    = EBX  + SPCO(ISPEC,J)*ECOS(ITH)
                EBY    = EBY  + SPCO(ISPEC,J)*ESIN(ITH)
                END DO
              E1(IK) = EBND * DTH
              APM(IK)= E1(IK) / ( TPI * GRAV**2 / SIG(IK)**5  )
              IF ( E1(IK) .GT. 1.E-5) THEN
                  THBND(IK) = MOD(630.- RADE*ATAN2(EBY,EBX),360.)
                  SPBND(IK) = RADE * SQRT ( MAX ( 0. , 2.*( 1. -      &
                    SQRT( MAX(0.,(EBX**2+EBY**2)/EBND**2) ) ) ) )
                ELSE
                  THBND(IK) = -999.9
                  SPBND(IK) = -999.9
                END IF
              EBND   = E1(IK) * DSII(IK) * TPIINV
              ET     = ET  + EBND
              EWN    = EWN + EBND / WN(IK)
              ETR    = ETR + EBND / SIG(IK)
              ETX    = ETX + EBX * DSII(IK)
              ETY    = ETY + EBY * DSII(IK)
              END DO
!
! tail factors for radian action etc ...!
!
            EBND   = E1(NK) * TPIINV / ( SIG(NK) * DTH )
            ET     = ET  + FTE *EBND
            EWN    = EWN + FTWL*EBND
            ETR    = ETR + FTTR*EBND
            ETX    = DTH*ETX*TPIINV + FTE*EBX*TPIINV/SIG(NK)
            ETY    = DTH*ETY*TPIINV + FTE*EBY*TPIINV/SIG(NK)
!
            HSIG   = 4. * SQRT ( ET )
            IF ( HSIG .GT. HSMIN ) THEN
                WLEN   = EWN / ET * TPI
                TMEAN  = ETR / ET * TPI
                THMEAN = MOD ( 630. - RADE*ATAN2(ETY,ETX) , 360. )
                THSPRD = RADE * SQRT ( MAX ( 0. , 2.*( 1. - SQRT(     &
                           MAX(0.,(ETX**2+ETY**2)/ET**2) ) ) ) )
              ELSE
                WLEN   = 0.
                TMEAN  = 0.
                THMEAN = 0.
                THSPRD = 0.
                DO IK=1, NK
                  E1(IK) = 0.
                  DO ITH=1, NTH
                    E(IK,ITH) = 0.
                    END DO
                  END DO
              END IF
!
! peak frequency
!
            EMAX   = E1(NK)
            IKM    = NK
!
            DO IK=NK-1, 1, -1
              IF ( E1(IK) .GT. EMAX ) THEN
                  EMAX   = E1(IK)
                  IKM    = IK
                END IF
              END DO
!
            IKL    = MAX (  1 , IKM-1 )
            IKH    = MIN ( NK , IKM+1 )
            EL     = E1(IKL) - E1(IKM)
            EH     = E1(IKH) - E1(IKM)
            DENOM  = XL*EH - XH*EL
!
            IF ( HSIG .GE. HSMIN ) THEN
                FP     = SIG(IKM) * ( 1. + 0.5 * ( XL2*EH - XH2*EL )  &
                            / SIGN ( MAX(ABS(DENOM),1.E-15) , DENOM ) )
                THP    = THBND(IKM)
                SPP    = SPBND(IKM)
              ELSE
                FP     = 0.
                THP    = 0.
                SPP    = 0.
              END IF
!
! 3.4 source terms
!
            DO IK=1, NK
              FACTOR = TPIINV * CG(IK) / SIG(IK)
              DO ITH=1, NTH
                ISPEC  = ITH + (IK-1)*NTH
                A(ITH,IK)   = FACTOR * SPCO(ISPEC,J)
                WN2(ITH,IK) = WN(IK)
                END DO
              END DO
!
            ZWND   = ZZWND
            TAUWX  = 0.
            TAUWY  = 0.
            LLWS(:,:)  = .TRUE.
                 USTAR  = 1.
!
            CALL W3SPR3 (A, CG, WN, EMEAN, FMEAN, FMEANS,        &
                         WNMEAN, AMAX, UABS, UDIRR, USTAR, USTD, &
                         TAUWX, TAUWY, CD, Z0, CHARN, LLWS, FMEANWS)
!
            DO ITT=1, 3
              CALL W3SIN3 (A, CG, WN2, UABS, USTAR, DAIR/DWAT,   &
                           ASO(J), UDIRR, Z0, CD, TAUWX, TAUWY,  &
                           TAUWNX, TAUWNY,                       &
                           ICE, XWI, DIA, LLWS, IX, IY )
            CALL W3SPR3 (A, CG, WN, EMEAN, FMEAN, FMEANS,        &
                         WNMEAN, AMAX, UABS, UDIRR, USTAR, USTD, &
                         TAUWX, TAUWY, CD, Z0, CHARN, LLWS, FMEANWS)
              END DO
!
            IF ( FLSRCE(2) ) THEN
                CALL W3SLN1 ( WN, FHIGH, USTAR, UDIRR, XLN )
!
                CALL W3SIN3 (A, CG, WN2, UABS, USTAR, DAIR/DWAT, &
                             ASO(J), UDIRR, Z0, CD,              &
                             TAUWX, TAUWY, TAUWNX, TAUWNY,       &
                             ICE, XWI, DIA, LLWS, IX, IY )
              END IF
            IF ( FLSRCE(3) ) THEN
                CALL W3SNL1 ( A, CG, WNMEAN*DEPTH,      XNL, DIA )
!!/NLS                CALL W3SNLS ( A, CG, WN, DEPTH, UABS, 900.,      &
!!/NLS                                                 SNL=XNL, AA=DIA )
!
              END IF
            IF ( FLSRCE(4) ) THEN
                CALL W3SDS3 ( A, WN, CG, EMEAN, FMEANS, WNMEAN,         &
                                        USTAR, USTD, DEPTH, XDS, DIA, IX, IY )
!
                CALL W3SDB1 ( A, WN2, DEPTH, EMEAN, FMEAN, WNMEAN,  &
                                                            XDB, DIA )
!
              END IF
            IF ( FLSRCE(5) ) THEN
 
                CALL W3SBT1 ( A, CG, WN, DEPTH,   XBT, DIA )
 
 
 
!
 
 
!
              END IF
           IF ( FLSRCE(6) ) THEN
           END IF
!
            DO IK=1, NK
              FACTOR = TPI / CG(IK) * SIG(IK)
              DO ITH=1, NTH
                ISPEC       = ITH + (IK-1)*NTH
                E  (IK,ITH) = SPCO(ISPEC,J)
                SWI(IK,ITH) = ( XWI(ITH,IK) + XLN(ITH,IK) ) * FACTOR
                SNL(IK,ITH) = ( XNL(ITH,IK) + XTR(ITH,IK) ) * FACTOR
                SDS(IK,ITH) = ( XDS(ITH,IK) + XDB(ITH,IK) ) * FACTOR
                SBT(IK,ITH) = ( XBT(ITH,IK) + XBS(ITH,IK) ) * FACTOR
                SIS(IK,ITH) = XIS(ITH,IK) * FACTOR
                STT(IK,ITH) = XXX(ITH,IK) * FACTOR
                END DO
              END DO
            STT    = STT + SWI + SNL + SDS + SBT + SIS
 
!
! 4.a Perform output
!
            IF ( FLSRCE(1) ) WRITE (NDSGRD)                           &
               ((E  (IK,ITH),ITH=1,NTH),IK=NK,1,-1)
            IF ( FLSRCE(2) ) WRITE (NDSGRD)                           &
               ((SWI(IK,ITH),ITH=1,NTH),IK=NK,1,-1)
            IF ( FLSRCE(3) ) WRITE (NDSGRD)                           &
               ((SNL(IK,ITH),ITH=1,NTH),IK=NK,1,-1)
            IF ( FLSRCE(4) ) WRITE (NDSGRD)                           &
               ((SDS(IK,ITH),ITH=1,NTH),IK=NK,1,-1)
            IF ( FLSRCE(5) ) WRITE (NDSGRD)                           &
               ((SBT(IK,ITH),ITH=1,NTH),IK=NK,1,-1)
            IF ( FLSRCE(6) ) WRITE (NDSGRD)                           &
               ((SIS(IK,ITH),ITH=1,NTH),IK=NK,1,-1)
            IF ( FLSRCE(7) ) WRITE (NDSGRD)                           &
               ((STT(IK,ITH),ITH=1,NTH),IK=NK,1,-1)
!
            IF ( FLAGLL ) THEN
                WRITE (NDSPNT,940) PTNME(J),                          &
                   FACT*PTLOC(1,J), FACT*PTLOC(2,J), DPO(J), WAO(J),  &
                   WAO(J)*COS(WDO(J)), WAO(J)*SIN(WDO(J)), ASO(J),    &
                   CAO(J), CAO(J)*COS(CDO(J)), CAO(J)*SIN(CDO(J)),    &
                   HSIG, GRDID(J)
              ELSE
                WRITE (NDSPNT,941) PTNME(J),                          &
                   FACT*PTLOC(1,J), FACT*PTLOC(2,J), DPO(J), WAO(J),  &
                   WAO(J)*COS(WDO(J)), WAO(J)*SIN(WDO(J)), ASO(J),    &
                   CAO(J), CAO(J)*COS(CDO(J)), CAO(J)*SIN(CDO(J)),    &
                   HSIG, GRDID(J)
              END IF
!
! ... End of points loop
!
          END IF
        END DO
!
      RETURN
!
! Formats
!
  905 FORMAT (9X,A)
 
  940 FORMAT (A10,1X,2F6.1,f7.1,3F7.1,F8.2,3F7.2,F6.2,2X,A)
 
  941 FORMAT (A10,1X,2F8.1,f7.1,3F7.1,F8.2,3F7.2,F6.2,2X,A)
 
!
!/
!/ End of GXEXPO ----------------------------------------------------- /
!/
      END SUBROUTINE GXEXPO
!/
!/ End of GXOUTP ----------------------------------------------------- /
!/
      END PROGRAM GXOUTP
