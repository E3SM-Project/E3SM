#include "w3macros.h"
!/ ------------------------------------------------------------------- /
     PROGRAM W3SHEL
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         11-May-2015 |
!/                  +-----------------------------------+
!/
!/    19-Jan-1999 : Final FORTRAN 77                    ( version 1.18 )
!/    19-Jan-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    08-Mar-2000 : Fix time managament bug.            ( version 2.04 )
!/    09-Jan-2001 : Fix FOUT allocation bug.            ( version 2.05 )
!/    24-Jan-2001 : Flat grid version.                  ( version 2.06 )
!/    25-Jan-2002 : Data assimilation set up.           ( version 2.17 )
!/    08-May-2002 : Clean up for timers.                ( version 2.21 )
!/    26-Aug-2002 : Generalizing timer.                 ( version 2.22 )
!/    26-Dec-2002 : Continuously moving grid.           ( version 3.02 )
!/    01-Aug-2003 : Continuously moving grid, input.    ( version 3.03 )
!/    07-Oct-2003 : Fixed NHMAX test.                   ( version 3.05 )
!/    05-Jan-2005 : Multiple grid version.              ( version 3.06 )
!/    04-May-2005 : Change to MPI_COMM[_WAVE.           ( version 3.07 )
!/    26-Jun-2006 : Add wiring for output type 6.       ( version 3.07 )
!/    28-Jun-2006 : Adding file name preamble.          ( version 3.09 )
!/    28-Oct-2006 : Adding partitioning options.        ( version 3.10 )
!/    21-Jun-2007 : Dedicated output processes.         ( version 3.11 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/    30-Oct-2009 : Fix format statement 2945.          ( version 3.14 )
!/                  (T. J. Campbell, NRL)
!/    30-Oct-2009 : Implement run-time grid selection.  ( version 3.14 )
!/                  (W. E. Rogers & T. J. Campbell, NRL)
!/    13-Sep-2009 : Add coupling option                 ( version 3.14_SHOM )
!/    30-Oct-2009 : Implement curvilinear grid type.    ( version 3.14 )
!/                  (W. E. Rogers & T. J. Campbell, NRL)
!/    29-Oct-2010 : Implement unstructured grids        ( version 3.14.4 )
!/                  (A. Roland and F. Ardhuin)
!/    23-Nov-2011 : Comments clean up                   ( version 4.04 )
!/    06-Mar-2012 : Repairing test output.              ( version 4.07 )
!/    03-Sep-2012 : Output initialization time.         ( version 4.10 )
!/    27-Sep-2012 : Implement use of tidal constituents ( version 4.08 )
!/    04-Feb-2014 : Switched clock to DATE_AND_TIME     ( version 4.18 )
!/                  (A. Chawla and Mark Szyszka)
!/    23-Apr-2015 : Adding NCEP Coupler                 ( version 5.06 )
!/                  (A. Chawla and Dmitry Sheinin)
!/    24-Apr-2015 : Adding OASIS coupling calls         ( version 5.07 )
!/                  (M. Accensi & F. Ardhuin, IFREMER)
!/    11-May-2015 : Checks dates for output types       ( version 5.08 )
!/
!/    Copyright 2009-2012 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     A generic shell for WAVEWATCH III, using preformatted
!     input fields.
!
!  2. Method :
!
!     Driver for the actual wave model (W3WAVE).
!
!     Files : ww3_shel.inp  Input commands for shell.
!             level.ww3     Water level fields (optional).
!             current.ww3   Current fields (optional).
!             wind.ww3      Wind fields (optional).
!             muddens.ww3   Mud parameter (optional)
!             mudthk.ww3    Mud parameter (optional)
!             mudvisc.ww3   Mud parameter (optional)
!             ice(n).ww3    Ice parameters (n=1 to 5) (optional)
!             ice.ww3       ice concentration fields (optional).
!             data0.ww3     Files with assimilation data (optional).
!             data1.ww3
!             data2.ww3
!
!     The file names of the input files are set in W3FLDO
!
!  3. Parameters :
!
!     Local parameters.
!     ----------------------------------------------------------------
!       NHMAX   I.P.  Maximum number of homogeneous fields.
!
!       NDSI    Int.  General input unit number (shell only).
!       NDSS    Int.  Scratch file.
!       NDSO    Int.  General output unit number (shell only).
!       NDSE    Int.  Error output unit number (shell only).
!       NDST    Int.  Test output unit number (shell only).
!       NDSF    I.A.  Field files unit numbers (shell only).
!       FLH     L.A.  Flags for homogeneous fields.
!       NH      I.A.  Number of times for homogeneous fields.
!       THO     I.A.  Times of homogeneous fields.
!       TIME0   I.A.  Starting time.
!       TIMEN   I.A.  Ending time.
!     ----------------------------------------------------------------
!
!       NDS, NTRACE, ..., see W3WAVE
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3NMOD    Subr. W3GDATMD Set nummber of data structures
!      W3SETG    Subr.   Id.    Point to data structure.
!      W3NDAT    Subr. W3WDATMD Set nummber of data structures
!      W3SETW    Subr.   Id.    Point to data structure.
!      W3NMOD    Subr. W3ADATMD Set nummber of data structures
!      W3NAUX    Subr.   Id.    Point to data structure.
!      W3NOUT    Subr. W3ODATMD Set nummber of data structures
!      W3SETO    Subr.   Id.    Point to data structure.
!      W3NINP    Subr. W3IDATMD Set nummber of data structures
!      W3SETI    Subr.   Id.    Point to data structure.
!
!      NEXTLN    Subr. W3SERVMD Skip to next input line.
!      STME21    Subr. W3TIMEMD Print date and time readable.
!      DSEC21    Func.   Id.    Difference between times.
!      TICK21    Subr.   Id.    Increment time.
!
!      W3FLDO    Subr. W3FLDSMD Opens and checks input files.
!      W3FLDG    Subr.   Id.    Reads from input files.
!      W3FLDD    Subr.   Id.    Reads from data files.
!      W3FLDH    Subr.   Id.    Udates homogeneous fields.
!
!      W3INIT    Subr. W3INITMD Wave model initialization.
!      W3READFLGRD Subr. W3IOGOMD Reading output fields flags.
!      W3WAVE    Subr. W3WAVEMD Wave model.
!      W3WDAS    Subr. W3WDASMD Data assimilation interface.
!
!      MPI_INIT, MPI_COMM_SIZE, MPI_COMM_RANK, MPI_BARRIER,
!         MPI_FINALIZE
!                Subr.          Standard MPI routines.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!     None, stand-alone program.
!
!  6. Error messages :
!
!     - Checks on I-O.
!     - Check on time interval.
!
!  7. Remarks :
!
!     - A rigourous input check is made in W3INIT.
!     - See W3WDAS for documentation on the set-up of the data
!       assimilation.
!     - in "7.a.2 Check if update is needed"
!       Field is updated when compute time is past old input time, and
!       (in case of homogeneous input field),  grabs field value at next
!       input time, which may in fact be far in the future from current
!       compute time. Example: user says
!       field=1   on 19680101 000000 and
!       field=100 on 20160101 000000
!       then on if 7.a.2 is reached on 19680101 010000, WW3 will set
!       field to 100.
!
!  8. Structure :
!
!     ----------------------------------------------------------------
!        0.   Set up data structures.                ( W3NMOD, etc. )
!        1.   I-O setup.
!          a  For shell.
!          b  For WAVEWATCH III.
!          c  Local parameters.
!        2.   Define input fields
!        3.   Set time frame.
!        4.   Define output
!          a  Loop over types, do
!        +--------------------------------------------------------+
!        | b    Process standard line                             |
!        | c    If type 1: fields of mean wave parameters         |
!        | d    If type 2: point output                           |
!        | e    If type 3: track output                           |
!        | f    If type 4: restart files                          |
!        | g    If type 5: boundary output                        |
!        | h    If type 6: separated wave fields                  |
!        | i    If type 7: coupling fields                        |
!        +--------------------------------------------------------+
!        5.   Initialzations
!          a  Wave model.                              ( W3INIT )
!          b  Read homogeneous field data.
!          c  Prepare input files.                     ( W3FLDO )
!          d  Set field times.
!        6.   If no input fields required, run model in a single
!             sweep and exit.                          ( W3WAVE )
!        7.   Run model with input
!             Do until end time is reached
!        +--------------------------------------------------------+
!        | a  Determine next time interval and input fields.      |
!        |   1  Preparation                                       |
!        |      Loop over input fields                            |
!        | +------------------------------------------------------|
!        | | 2  Check if update is needed                         |
!        | | 3  Update time and fields                 ( W3FLDG ) |
!        | |                                           ( W3FLDH ) |
!        | | 4  Update next ending time                           |
!        | +------------------------------------------------------|
!        | b  Run wave model.                          ( W3WAVE ) |
!        | c  If requested, data assimilation.         ( W3WDAS ) |
!        | d  Final output if needed.                  ( W3WAVE ) |
!        | e  Check time                                          |
!        +--------------------------------------------------------+
!     ----------------------------------------------------------------
!
!  9. Switches :
!
!       !/SHRD  Switch for shared / distributed memory architecture.
!       !/DIST  Id.
!       !/MPI   Id.
!
!       !/MGW   Moving grid wind correction.
!       !/MGP   Moving grid propagation correction.
!
!       !/T     Enable test output.
!       !/O7    Echo input homogeneous fields.
!
!       !/NCO   NCEP NCO modifications for operational implementation.
!
!       !/F90   Timer function included for F90.
!
!       !/NCC   Ncep Coupler
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3GDATMD
      USE W3WDATMD, ONLY: TIME, W3NDAT, W3DIMW, W3SETW
      USE W3ADATMD, ONLY: W3NAUX, W3DIMA, W3SETA
      USE W3IDATMD
      USE W3ODATMD, ONLY: W3NOUT, W3SETO
      USE W3ODATMD, ONLY: NAPROC, IAPROC, NAPOUT, NAPERR, NOGRP,      &
                          NGRPP, IDOUT, FNMPRE, IOSTYP, NOTYPE, NOGE
!/
      USE W3FLDSMD
      USE W3INITMD
      USE W3WAVEMD
      USE W3WDASMD
!/
      USE W3IOGRMD, ONLY: W3IOGR
      USE W3IOGOMD, ONLY: W3READFLGRD, FLDOUT
      USE W3IOPOMD
      USE W3SERVMD, ONLY : NEXTLN, EXTCDE
      USE W3TIMEMD
!
      IMPLICIT NONE
!
      INCLUDE "mpif.h"
!/
!/ ------------------------------------------------------------------- /
!/ Local PARAMETER statements
!/
      INTEGER, PARAMETER  :: NHMAX =    200
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER             :: NDSI, NDSI2, NDSS, NDSO, NDSE, NDST,     &
                             NDSF(-7:7), NDSEN, NDS(13), NTRACE(2),   &
                             TIME0(2), TIMEN(2), TTIME(2), TTT(2),    &
                             IERR, J, I, ODAT(35), ILOOP, NPTS,       &
                             NH(-7:4), THO(2,-7:4,NHMAX), RCLD(5:7),  &
                             NDT(5:7), NDTNEW, MPI_COMM = -99, JJ,    &
                             IPRT(6) = 0, IFI, IFJ
      INTEGER             :: CLKDT1(8), CLKDT2(8), CLKDT3(8)
      INTEGER             :: IERR_MPI
      INTEGER             :: NODATA(5:7), FLAGTIDE
      INTEGER             :: COUPL_COMM
      REAL                :: FACTOR, DTTST, XX, YY,                   &
                             HA(NHMAX,-7:4), HD(NHMAX,-7:4),          &
                             HS(NHMAX,-7:4)
 
      REAL                :: CLKFIN, CLKFEL
      REAL, ALLOCATABLE   :: X(:), Y(:), XXX(:,:), DATA0(:,:),        &
                             DATA1(:,:), DATA2(:,:)
      LOGICAL             :: FLLSTL, FLLSTI, FLH(-7:8), FLFLG, FLHOM, &
                             TFLAGI, FLGDAS(3), FLGRD(NOGRP,NGRPP),   &
                             FLT, FLGD(NOGRP)
      LOGICAL             :: FLGR2(NOGRP,NGRPP), FLG2(NOGRP)
      LOGICAL             :: FLAGSTIDE(4)
      LOGICAL             :: PRTFRM
      LOGICAL             :: FLLST_ALL(-7:8)
      LOGICAL             :: DEBUG_NCC = .FALSE.
 
      CHARACTER(LEN=1)    :: COMSTR,FLAGTFC
      CHARACTER(LEN=3)    :: IDSTR(-7:8), IDTST
      CHARACTER(LEN=6)    :: YESXNO
      CHARACTER(LEN=10)   :: PN
      CHARACTER(LEN=10),                                              &
              ALLOCATABLE :: PNAMES(:)
      CHARACTER(LEN=13)   :: IDFLDS(-7:8)
      CHARACTER(LEN=20)   :: STRNG
      CHARACTER(LEN=23)   :: DTME21
      CHARACTER(LEN=30)   :: IDOTYP(7)
      CHARACTER(LEN=80)   :: LINE
      CHARACTER(LEN=1024) :: FLDIN
!/
!/ ------------------------------------------------------------------- /
!/
      DATA IDFLDS / 'ice param. 1 ' , 'ice param. 2 ' ,               &
                    'ice param. 3 ' , 'ice param. 4 ' ,               &
                    'ice param. 5 ' ,                                 &
                    'mud density  ' , 'mud thkness  ' ,               &
                    'mud viscos.  ' ,                                 &
                    'water levels ' , 'currents     ' ,               &
                    'winds        ' , 'ice fields   ' ,               &
                    'mean param.  ' , '1D spectra   ' ,               &
                    '2D spectra   ' , 'moving grid  ' /
      DATA IDOTYP / 'Fields of mean wave parameters' ,                &
                    'Point output                  ' ,                &
                    'Track point output            ' ,                &
                    'Restart files                 ' ,                &
                    'Nesting data                  ' ,                &
                    'Partitioned wave field data   ' ,                &
                    'Fields for coupling           ' /
      DATA IDSTR  / 'IC1', 'IC2', 'IC3', 'IC4', 'IC5', 'MDN', 'MTH',  &
                    'MVS', 'LEV', 'CUR', 'WND', 'ICE', 'DT0', 'DT1',  &
                    'DT2', 'MOV' /
!
!     IF (FLAGLL) THEN
!         FACTOR = 1.
!       ELSE
!         FACTOR = 1.E-3
!       END IF
!
      FLAGSTIDE(:) = .FALSE.
      FLH(:)       = .FALSE.
!
      CALL DATE_AND_TIME ( VALUES=CLKDT1 )
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 0.  Set up data structures
!
      CALL W3NMOD ( 1, 6, 6 )
      CALL W3NDAT (    6, 6 )
      CALL W3NAUX (    6, 6 )
      CALL W3NOUT (    6, 6 )
      CALL W3NINP (    6, 6 )
!
      CALL W3SETG ( 1, 6, 6 )
      CALL W3SETW ( 1, 6, 6 )
      CALL W3SETA ( 1, 6, 6 )
      CALL W3SETO ( 1, 6, 6 )
      CALL W3SETI ( 1, 6, 6 )
!
      CALL MPI_INIT      ( IERR_MPI )
      MPI_COMM = MPI_COMM_WORLD
!
      CALL MPI_COMM_SIZE ( MPI_COMM, NAPROC, IERR_MPI )
      CALL MPI_COMM_RANK ( MPI_COMM, IAPROC, IERR_MPI )
      IAPROC = IAPROC + 1
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 1.  IO set-up
! 1.a For shell
!
      NDSI   = 10
      NDSS   = 90
      NDSO   =  6
      NDSE   =  6
      NDST   =  6
 
 
      NDSF(-7)  = 1008
      NDSF(-6)  = 1009
      NDSF(-5)  = 1010
      NDSF(-4)  = 1011
      NDSF(-3)  = 1012
      NDSF(-2)  = 1013
      NDSF(-1)  = 1014
      NDSF(0)   = 1015
 
      NDSF(1)  = 11
      NDSF(2)  = 12
      NDSF(3)  = 13
      NDSF(4)  = 14
      NDSF(5)  = 15
      NDSF(6)  = 16
      NDSF(7)  = 17
!
      NAPOUT = 1
      NAPERR = 1
!
 
      IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,900)
!
      IF ( IAPROC .EQ. NAPERR ) THEN
          NDSEN  = NDSE
        ELSE
          NDSEN  = -1
        END IF
!
      JJ     = LEN_TRIM(FNMPRE)
      OPEN (NDSI,FILE=FNMPRE(:JJ)//'ww3_shel.inp',STATUS='OLD',       &
            ERR=2000,IOSTAT=IERR)
      REWIND (NDSI)
      READ (NDSI,'(A)',END=2001,ERR=2002) COMSTR
      IF (COMSTR.EQ.' ') COMSTR = '$'
!
      IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,901) COMSTR
!
! 1.b For WAVEWATCH III (See W3INIT)
!
      NDS( 1) = 20
      NDS( 2) =  6
!     NDS( 3) =  6
      NDS( 3) = 21
      NDS( 4) =  6
      NDS( 5) = 30
      NDS( 6) = 30
      NDS( 7) = 31
      NDS( 8) = 32
      NDS( 9) = 33
      NDS(10) = 35
      NDS(11) = 22
      NDS(12) = 23
      NDS(13) = 34
!
      NTRACE(1) =  NDS(3)
      NTRACE(2) =  10
!
! 1.c Local parameters
!
! inferred from context: these flags (FL) are to indicate that the last (LST)
!   field has been read from a file.
      FLLSTL = .FALSE. ! This is associated with J.EQ.1 (wlev)
      FLLSTI = .FALSE. ! This is associated with J.EQ.4 (ice)
      FLLST_ALL = .FALSE. ! For all
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 2.  Define input fields
!
      IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,920)
!
 
! If using experimental mud or ice physics, additional lines will
!  be read in from ww3_shel.inp and applied, so JFIRST is changed from
!  its initialization setting "JFIRST=1" to some lower value.
 
      DO J=JFIRST, 7
        CALL NEXTLN ( COMSTR , NDSI , NDSEN )
        IF ( J .LT. 4 ) THEN
            READ (NDSI,*,END=2001,ERR=2002) FLAGTFC, FLH(J)
            IF (FLAGTFC.EQ.'T') THEN
              INFLAGS1(J)=.TRUE.
              FLAGSC(J)=.FALSE.
            END IF
            IF (FLAGTFC.EQ.'F') THEN
              INFLAGS1(J)=.FALSE.
              FLAGSC(J)=.FALSE.
            END IF
            IF (FLAGTFC.EQ.'C') THEN
              INFLAGS1(J)=.TRUE.
              FLAGSC(J)=.TRUE.
            END IF
            FLH(J) = FLH(J) .AND. INFLAGS1(J)
          ELSE
            READ (NDSI,*,END=2001,ERR=2002) INFLAGS1(J)
            FLH(J) = .FALSE.
          END IF
        IF ( INFLAGS1(J) ) THEN
            YESXNO = 'YES/--'
          ELSE
            YESXNO = '---/NO'
          END IF
        IF ( FLH(J) ) THEN
            STRNG  = '(homogeneous field) '
          ELSE IF ( FLAGSC(J) ) THEN
            STRNG  = '(coupling field) '
          ELSE
            STRNG  = '                    '
          END IF
        IF ( IAPROC .EQ. NAPOUT )                                     &
              WRITE (NDSO,921) IDFLDS(J), YESXNO, STRNG
        END DO
!
      INFLAGS1(8) = .FALSE.
      FLH(8)   = .FALSE.
      IF ( INFLAGS1(8) .AND. IAPROC.EQ.NAPOUT )                          &
           WRITE (NDSO,921) IDFLDS(8), 'YES/--', ' '
!
      FLFLG  = INFLAGS1(-7) .OR. INFLAGS1(-6) .OR. INFLAGS1(-5) .OR. INFLAGS1(-4) &
               .OR. INFLAGS1(-3) .OR. INFLAGS1(-2) .OR. INFLAGS1(-1)           &
               .OR. INFLAGS1(0)  .OR. INFLAGS1(1)  .OR. INFLAGS1(2)            &
               .OR. INFLAGS1(3)  .OR. INFLAGS1(4)  .OR. INFLAGS1(5)            &
               .OR. INFLAGS1(6)  .OR. INFLAGS1(7)
      FLHOM  = FLH(-7) .OR. FLH(-6) .OR. FLH(-5) .OR. INFLAGS1(-4)       &
               .OR. FLH(-3) .OR. FLH(-2) .OR. FLH(-1) .OR. INFLAGS1(0)   &
               .OR. FLH(1) .OR. FLH(2) .OR. FLH(3) .OR. INFLAGS1(8)
!
      IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,922)
!
!     INFLAGS2 is just "initial value of INFLAGS1", i.e. does *not* get
!        changed when model reads last record of ice.ww3
      INFLAGS2=INFLAGS1
 
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 3.  Set time frame
!
      IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,930)
!
      CALL NEXTLN ( COMSTR , NDSI , NDSEN )
      READ (NDSI,*,END=2001,ERR=2002) TIME0
      CALL STME21 ( TIME0 , DTME21 )
      IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,931) DTME21
      TIME = TIME0
!
      CALL NEXTLN ( COMSTR , NDSI , NDSEN )
      READ (NDSI,*,END=2001,ERR=2002) TIMEN
      CALL STME21 ( TIMEN , DTME21 )
      IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,932) DTME21
!
      DTTST  = DSEC21 ( TIME0 , TIMEN )
      IF ( DTTST .LE. 0. ) GOTO 2003
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 4.  Define output
!
      CALL W3IOGR ( 'GRID', NDSF(5) )
      IF ( FLAGLL ) THEN
          FACTOR = 1.
        ELSE
          FACTOR = 1.E-3
        END IF
!
      CALL NEXTLN ( COMSTR , NDSI , NDSEN )
      READ (NDSI,*,END=2001,ERR=2002) IOSTYP
      IOSTYP = MAX ( 0 , MIN ( 3 , IOSTYP ) )
!
      IF ( IAPROC .EQ. NAPOUT ) THEN
          IF ( IOSTYP .EQ. 0 ) THEN
              WRITE (NDSO,940) 'No dedicated output process, ' //   &
                               'parallel file system required.'
          ELSE IF ( IOSTYP .EQ. 1 ) THEN
              WRITE (NDSO,940) 'No dedicated output process, ' //   &
                               'any file system.'
          ELSE IF ( IOSTYP .EQ. 2 ) THEN
              WRITE (NDSO,940) 'Single dedicated output process.'
          ELSE IF ( IOSTYP .EQ. 3 ) THEN
              WRITE (NDSO,940) 'Multiple dedicated output processes.'
          ELSE
              WRITE (NDSO,940) 'IOSTYP NOT RECOGNIZED'
          END IF
        END IF
!
! 4.a Loop over types
!
      NPTS   = 0
!
      NOTYPE = 6
      DO J = 1, NOTYPE
!
! 4.b Process standard line
!
        CALL NEXTLN ( COMSTR , NDSI , NDSEN )
        READ (NDSI,*,END=2001,ERR=2002) (ODAT(I),I=5*(J-1)+1,5*J)
        ODAT(5*(J-1)+3) = MAX ( 0 , ODAT(5*(J-1)+3) )
!
        IF ( ODAT(5*(J-1)+3) .NE. 0 ) THEN
            IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,941) J, IDOTYP(J)
            TTIME(1) = ODAT(5*(J-1)+1)
            TTIME(2) = ODAT(5*(J-1)+2)
            CALL STME21 ( TTIME , DTME21 )
            IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,942) DTME21
            TTIME(1) = ODAT(5*(J-1)+4)
            TTIME(2) = ODAT(5*(J-1)+5)
            CALL STME21 ( TTIME , DTME21 )
            IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,943) DTME21
            TTIME(1) = 0
            TTIME(2) = 0
            DTTST    = REAL ( ODAT(5*(J-1)+3) )
            CALL TICK21 ( TTIME , DTTST  )
            CALL STME21 ( TTIME , DTME21 )
            IF ( ( ODAT(5*(J-1)+1) .NE. ODAT(5*(J-1)+4) .OR.          &
                   ODAT(5*(J-1)+2) .NE. ODAT(5*(J-1)+5) ) .AND.       &
                   IAPROC .EQ. NAPOUT ) THEN
                IF ( DTME21(9:9) .NE. '0' ) THEN
                    WRITE (NDSO,1944) DTME21( 9:19)
                  ELSE IF ( DTME21(10:10) .NE. '0' ) THEN
                    WRITE (NDSO,2944) DTME21(10:19)
                  ELSE
                    WRITE (NDSO,3944) DTME21(12:19)
                  END IF
              END IF
!
            IF ( J .EQ. 1 ) THEN
!
! 4.c Type 1: fields of mean wave parameters
!
              CALL W3READFLGRD ( NDSI, NDSO, 9, NDSEN, COMSTR, FLGD,   &
                                 FLGRD, IAPROC, NAPOUT, IERR )
              IF ( IERR .NE. 0 ) GOTO 2222
!
              ELSE IF ( J .EQ. 2 ) THEN
!
! 4.d Type 2: point output
!
                DO ILOOP=1, 2
                  JJ     = LEN_TRIM(FNMPRE)
                  IF ( ILOOP .EQ. 1 ) THEN
                      NDSI2  = NDSI
                      IF ( IAPROC .EQ. 1 ) OPEN                       &
                          (NDSS,FILE=FNMPRE(:JJ)//'ww3_shel.scratch')
                    ELSE
                      NDSI2  = NDSS
                      CALL MPI_BARRIER (MPI_COMM,IERR_MPI)
                      OPEN (NDSS,FILE=FNMPRE(:JJ)//'ww3_shel.scratch')
                      REWIND (NDSS)
!
                      IF (NPTS.GT.0) THEN
                         ALLOCATE ( X(NPTS), Y(NPTS), PNAMES(NPTS) )
                      ELSE
                         GOTO 2004
                      END IF
                    END IF
!
                  NPTS   = 0
                  DO
                    CALL NEXTLN ( COMSTR , NDSI , NDSEN )
                    READ (NDSI2,*,END=2001,ERR=2002) XX, YY, PN
                    IF ( ILOOP.EQ.1 .AND. IAPROC.EQ.1 ) THEN
                        BACKSPACE (NDSI)
                        READ (NDSI,'(A)') LINE
                        WRITE (NDSS,'(A)') LINE
                      END IF
                    IF ( PN .EQ. 'STOPSTRING' ) EXIT
                    NPTS   = NPTS + 1
                    IF ( ILOOP .EQ. 1 ) CYCLE
                    X(NPTS)      = XX
                    Y(NPTS)      = YY
                    PNAMES(NPTS) = PN
                    IF ( IAPROC .EQ. NAPOUT ) THEN
                        IF ( FLAGLL ) THEN
                            IF ( NPTS .EQ. 1 ) THEN
                                WRITE (NDSO,2945)                     &
                                              FACTOR*XX, FACTOR*YY, PN
                              ELSE
                                WRITE (NDSO,2946) NPTS,               &
                                              FACTOR*XX, FACTOR*YY, PN
                              END IF
                          ELSE
                            IF ( NPTS .EQ. 1 ) THEN
                                WRITE (NDSO,2955)                     &
                                              FACTOR*XX, FACTOR*YY, PN
                              ELSE
                                WRITE (NDSO,2956) NPTS,               &
                                              FACTOR*XX, FACTOR*YY, PN
                              END IF
                          END IF
                      END IF
                    END DO
!
                  IF ( IAPROC.EQ.1 .AND. ILOOP.EQ.1 ) CLOSE (NDSS)
                  END DO
!
                IF ( NPTS.EQ.0 .AND. IAPROC.EQ.NAPOUT )               &
                     WRITE (NDSO,2947)
                IF ( IAPROC .EQ. 1 ) THEN
                    CALL MPI_BARRIER ( MPI_COMM, IERR_MPI )
                    CLOSE (NDSS,STATUS='DELETE')
                  ELSE
                    CLOSE (NDSS)
                    CALL MPI_BARRIER ( MPI_COMM, IERR_MPI )
                  END IF
!
              ELSE IF ( J .EQ. 3 ) THEN
!
! 4.e Type 3: track output
!
                CALL NEXTLN ( COMSTR , NDSI , NDSEN )
                READ (NDSI,*,END=2001,ERR=2002) TFLAGI
!
                IF ( .NOT. TFLAGI ) NDS(11) = -NDS(11)
                IF ( IAPROC .EQ. NAPOUT ) THEN
                    IF ( .NOT. TFLAGI ) THEN
                        WRITE (NDSO,3945) 'input', 'UNFORMATTED'
                      ELSE
                        WRITE (NDSO,3945) 'input', 'FORMATTED'
                      END IF
                  END IF
!
              ELSE IF ( J .EQ. 4 ) THEN
!
! 4.f Type 4: restart files (no additional data)
!
              ELSE IF ( J .EQ. 5 ) THEN
!
! 4.g Type 5: nesting data (no additional data)
!
              ELSE IF ( J .EQ. 6 ) THEN
!
! 4.h Type 6: partitioning
!
!             IPRT: IX0, IXN, IXS, IY0, IYN, IYS
!
                CALL NEXTLN ( COMSTR , NDSI , NDSEN )
                READ (NDSI,*,END=2001,ERR=2002) IPRT, PRTFRM
!
                IF ( IAPROC .EQ. NAPOUT ) THEN
                    IF ( PRTFRM ) THEN
                        YESXNO = 'YES/--'
                      ELSE
                        YESXNO = '---/NO'
                      END IF
                    WRITE (NDSO,6945) IPRT, YESXNO
                  END IF
!
              END IF
!
          END IF
!
        END DO
!
      IF ( NPTS.EQ.0 ) ALLOCATE ( X(1), Y(1), PNAMES(1) )
!
! ... End loop over output types
!
! For outputs with non-zero time step, check dates :
! If output ends before run start OR output starts after run end,
! deactivate output cleanly with output time step = 0
! This is usefull for IOSTYP=3 (Multiple dedicated output processes)
! to avoid the definition of dedicated proc. for unused output.
!
      DO J = 1, NOTYPE
        DTTST  = DSEC21 ( TIME0 , ODAT(5*(J-1)+4:5*(J-1)+5) )
        IF ( DTTST .LT. 0 ) THEN
          ODAT(5*(J-1)+3) = 0
          IF ( IAPROC .EQ. NAPOUT )  WRITE (NDSO,8945) trim(IDOTYP(J))
          CONTINUE
        ENDIF
        DTTST  = DSEC21 ( ODAT(5*(J-1)+1:5*(J-1)+2), TIMEN )
        IF ( DTTST .LT. 0 ) THEN
          ODAT(5*(J-1)+3) = 0
          IF ( IAPROC .EQ. NAPOUT )  WRITE (NDSO,8945) trim(IDOTYP(J))
          CONTINUE
        ENDIF
      ENDDO
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 5.  Initializations
!
      IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,950)
 
!
! 5.a Opening field and data files
!
      IF ( FLFLG ) THEN
          IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,951)                  &
                                          'Preparing input files ...'
!
 
          DO J=JFIRST, 4
            IF ( INFLAGS1(J) .AND. .NOT. FLAGSC(J)) THEN
                IF ( FLH(J) ) THEN
                    IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,954) IDFLDS(J)
                  ELSE
                    JJ     = LEN_TRIM(FNMPRE)
                    FLAGTIDE = 0
                    CALL W3FLDO ('READ', IDSTR(J), NDSF(J), NDST,     &
                                  NDSEN, NX, NY, GTYPE,               &
                                  IERR, FPRE=FNMPRE(:JJ), TIDEFLAGIN=FLAGTIDE )
                    IF ( IERR .NE. 0 ) GOTO 2222
                    IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,955) IDFLDS(J)
                  END IF
              ELSE
                IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,954) IDFLDS(J)
              END IF
          END DO
!
          DO J=5, 7
            IF ( INFLAGS1(J) .AND. .NOT. FLAGSC(J)) THEN
                CALL W3FLDO ('READ', IDSTR(J), NDSF(J), NDST, NDSEN, &
                             RCLD(J), NY, NODATA(J),                 &
                             IERR, FPRE=FNMPRE(:JJ) )
                IF ( IERR .NE. 0 ) GOTO 2222
                IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,956) IDFLDS(J),&
                             RCLD(J), NODATA(J)
              ELSE
                IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,954) IDFLDS(J)
              END IF
          END DO
!
        END IF
 
!
! 5.b Wave model
!
      IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,951) 'Wave model ...'
!
     CALL W3INIT ( 1, 'ww3', NDS, NTRACE, ODAT, FLGRD, FLGR2, FLGD,    &
                   FLG2, NPTS, X, Y, PNAMES, IPRT, PRTFRM, MPI_COMM,   &
                   FLAGSTIDEIN=FLAGSTIDE )
!
      ALLOCATE ( XXX(NX,NY) )
!
! 5.c Homogeneous field data
!
      IF ( FLHOM ) THEN
          IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,951)                  &
                        'Homogeneous field data (and moving grid) ...'
          NH     = 0
!
! ... Start of loop.
!
          DO
            CALL NEXTLN ( COMSTR , NDSI , NDSEN )
            READ (NDSI,*,END=2001,ERR=2002) IDTST
 
            IF ( IDTST.NE.IDSTR(-7) .AND. IDTST.NE.IDSTR(-6) .AND.   &
                 IDTST.NE.IDSTR(-5) .AND. IDTST.NE.IDSTR(-4) .AND.   &
                 IDTST.NE.IDSTR(-3) .AND. IDTST.NE.IDSTR(-2) .AND.   &
                 IDTST.NE.IDSTR(-1) .AND. IDTST.NE.IDSTR(0)  .AND.   &
                 IDTST.NE.IDSTR(1)  .AND. IDTST.NE.IDSTR(2)  .AND.   &
                 IDTST.NE.IDSTR(3)  .AND. IDTST.NE.IDSTR(8)  .AND.   &
                 IDTST.NE.'STP' ) GOTO 2005
 
!
! ... Stop conditions
!
            IF ( IDTST .EQ. 'STP' ) THEN
                EXIT
              ELSE
                BACKSPACE ( NDSI )
              END IF
!
! ... Store data
!
            DO J=LBOUND(IDSTR,1), 4
              I      = J
              IF ( J .EQ. 4 ) I = 8
              IF ( IDTST .EQ. IDSTR(I) ) THEN
                  NH(J)    = NH(J) + 1
                  IF ( NH(J) .GT. NHMAX ) GOTO 2006
                  IF ( J .LE. 1  ) THEN ! water levels, etc. : get HS
                     READ (NDSI,*,END=2001,ERR=2002) IDTST,           &
                           THO(1,J,NH(J)), THO(2,J,NH(J)),            &
                           HS(NH(J),J)
                    ELSE IF ( J .EQ. 2 ) THEN ! currents: get HA and HD
                     READ (NDSI,*,END=2001,ERR=2002) IDTST,           &
                           THO(1,J,NH(J)), THO(2,J,NH(J)),            &
                           HA(NH(J),J), HD(NH(J),J)
                    ELSE IF ( J .EQ. 3 ) THEN ! wind: get HA HD and HS
                     READ (NDSI,*,END=2001,ERR=2002) IDTST,           &
                           THO(1,J,NH(J)), THO(2,J,NH(J)),            &
                           HA(NH(J),J), HD(NH(J),J), HS(NH(J),J)
                    ELSE IF ( J .EQ. 4 ) THEN ! ice: HA and HD
                     READ (NDSI,*,END=2001,ERR=2002) IDTST,           &
                           THO(1,J,NH(J)), THO(2,J,NH(J)),            &
                           HA(NH(J),J), HD(NH(J),J)
                    END IF
                END IF
              END DO
!
            END DO
!
! ... End of loop, output
!
          DO J=JFIRST, 3
            IF ( FLH(J) .AND. IAPROC.EQ.NAPOUT ) THEN
                WRITE (NDSO,952) NH(J), IDFLDS(J)
                DO I=1, NH(J)
                  IF ( J .LE. 1 ) THEN
                      WRITE (NDSO,953) I, THO(1,J,I), THO(2,J,I), &
                                          HS(I,J)
                    ELSE IF ( J .EQ. 2 ) THEN
                      WRITE (NDSO,953) I, THO(1,J,I), THO(2,J,I), &
                                          HA(I,J), HD(I,J)
                    ELSE IF ( J .EQ. 3 ) THEN
                      WRITE (NDSO,953) I, THO(1,J,I), THO(2,J,I), &
                                          HA(I,J), HD(I,J), HS(I,J)
                    END IF
                  END DO
              END IF
            END DO
!
          IF ( INFLAGS1(8) .AND. IAPROC.EQ.NAPOUT ) THEN
              WRITE (NDSO,952) NH(4), IDFLDS(8)
              DO I=1, NH(J)
                WRITE (NDSO,953) I, THO(1,4,I), THO(2,4,I),       &
                                    HA(I,4), HD(I,4)
                END DO
              END IF
!
          IF ( ( FLH(-7) .AND. (NH(-7).EQ.0) ) .OR.                     &
               ( FLH(-6) .AND. (NH(-6).EQ.0) ) .OR.                     &
               ( FLH(-5) .AND. (NH(-5).EQ.0) ) .OR.                     &
               ( FLH(-4) .AND. (NH(-4).EQ.0) ) .OR.                     &
               ( FLH(-3) .AND. (NH(-3).EQ.0) ) .OR.                     &
               ( FLH(-2) .AND. (NH(-2).EQ.0) ) .OR.                     &
               ( FLH(-1) .AND. (NH(-1).EQ.0) ) .OR.                     &
               ( FLH(0)  .AND. (NH(0).EQ.0)  ) .OR.                     &
               ( FLH(1)  .AND. (NH(1).EQ.0)  ) .OR.                     &
               ( FLH(2)  .AND. (NH(2).EQ.0)  ) .OR.                     &
               ( FLH(3)  .AND. (NH(3).EQ.0)  ) .OR.                     &
               ( INFLAGS1(8) .AND. (NH(4).EQ.0) ) ) GOTO 2007
!
        END IF
!
      CALL MPI_BARRIER ( MPI_COMM, IERR_MPI )
!
      IF ( IAPROC .EQ. NAPOUT ) THEN
          CALL DATE_AND_TIME ( VALUES=CLKDT2 )
        END IF
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
 
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 6.  Model without input
!
      IF ( .NOT. FLFLG ) THEN
!
          IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,960)
          CALL W3WAVE ( 1, TIMEN                      &
                      )
!
          GOTO 2222
!
        END IF
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 7.  Model with input
!
      IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,970)
!
  700 CONTINUE
!
! 7.a Determine next time interval and input fields
! 7.a.1 Preparation
!
      TTIME  = TIMEN
!
      CALL STME21 ( TIME0 , DTME21 )
      IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,971) DTME21
!
 
      DO J=JFIRST,8
!
        IF ( INFLAGS1(J) ) THEN
!
! 7.a.2 Check if update is needed
            IF (.NOT.FLAGSC(J)) THEN
              TTT(1) = TFN(1,J)
              TTT(2) = TFN(2,J)
              IF ( TTT(1) .EQ. -1 ) THEN
                DTTST  = 0.
              ELSE
                DTTST  = DSEC21 ( TIME0 , TTT )
              END IF
            ELSE
             END IF
!
! 7.a.3 Update time and fields / data
!
            IF ( DTTST .LE. 0. ) THEN
                IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,972) IDFLDS(J)
!
! IC1 : (in context of IC3, this is ice thickness)
                IF ( J .EQ. -7 ) THEN
                   IF ( FLH(J) ) THEN
                      CALL W3FLDH (J, NDST, NDSEN, NX, NY, NX, NY,    &
                           TIME0, TIMEN, NH(J), NHMAX, THO, HA, HD, HS,&
                           TTT, XXX, XXX, XXX, TI1, XXX, XXX, ICEP1, IERR)
                   ELSE
                      CALL W3FLDG ('READ', IDSTR(J), NDSF(J),         &
                           NDST, NDSEN, NX, NY, NX, NY, TIME0, TIMEN, &
                           TTT, XXX, XXX, XXX, TI1, XXX, XXX, ICEP1, IERR)
                   END IF
                   IF ( IERR .LT. 0 )FLLST_ALL(J) = .TRUE.
! IC2 : (in context of IC3, this is ice viscosity)
                ELSE IF ( J .EQ. -6 ) THEN
                   IF ( FLH(J) ) THEN
                      CALL W3FLDH (J, NDST, NDSEN, NX, NY, NX, NY,    &
                           TIME0, TIMEN, NH(J), NHMAX, THO, HA, HD, HS,&
                           TTT, XXX, XXX, XXX, TI2, XXX, XXX, ICEP2, IERR)
                   ELSE
                      CALL W3FLDG ('READ', IDSTR(J), NDSF(J),         &
                           NDST, NDSEN, NX, NY, NX, NY, TIME0, TIMEN, &
                           TTT, XXX, XXX, XXX, TI2, XXX, XXX, ICEP2, IERR)
                   END IF
                   IF ( IERR .LT. 0 )FLLST_ALL(J) = .TRUE.
! IC3 : (in context of IC3, this is ice density)
                ELSE IF ( J .EQ. -5 ) THEN
                   IF ( FLH(J) ) THEN
                      CALL W3FLDH (J, NDST, NDSEN, NX, NY, NX, NY,    &
                           TIME0, TIMEN, NH(J), NHMAX, THO, HA, HD, HS,&
                           TTT, XXX, XXX, XXX, TI3, XXX, XXX, ICEP3, IERR)
                   ELSE
                      CALL W3FLDG ('READ', IDSTR(J), NDSF(J),         &
                           NDST, NDSEN, NX, NY, NX, NY, TIME0, TIMEN, &
                           TTT, XXX, XXX, XXX, TI3, XXX, XXX, ICEP3, IERR)
                   END IF
                   IF ( IERR .LT. 0 )FLLST_ALL(J) = .TRUE.
! IC4 : (in context of IC3, this is ice modulus)
                ELSE IF ( J .EQ. -4 ) THEN
                   IF ( FLH(J) ) THEN
                      CALL W3FLDH (J, NDST, NDSEN, NX, NY, NX, NY,    &
                           TIME0, TIMEN, NH(J), NHMAX, THO, HA, HD, HS,&
                           TTT, XXX, XXX, XXX, TI4, XXX, XXX, ICEP4, IERR)
                   ELSE
                      CALL W3FLDG ('READ', IDSTR(J), NDSF(J),         &
                           NDST, NDSEN, NX, NY, NX, NY, TIME0, TIMEN, &
                           TTT, XXX, XXX, XXX, TI4, XXX, XXX, ICEP4, IERR)
                   END IF
                   IF ( IERR .LT. 0 )FLLST_ALL(J) = .TRUE.
! IC5 : ice flow diam.
                ELSE IF ( J .EQ. -3 ) THEN
                   IF ( FLH(J) ) THEN
                      CALL W3FLDH (J, NDST, NDSEN, NX, NY, NX, NY,    &
                           TIME0, TIMEN, NH(J), NHMAX, THO, HA, HD, HS,&
                           TTT, XXX, XXX, XXX, TI5, XXX, XXX, ICEP5, IERR)
                   ELSE
                      CALL W3FLDG ('READ', IDSTR(J), NDSF(J),         &
                           NDST, NDSEN, NX, NY, NX, NY, TIME0, TIMEN, &
                           TTT, XXX, XXX, XXX, TI5, XXX, XXX, ICEP5, IERR)
                   END IF
                   IF ( IERR .LT. 0 )FLLST_ALL(J) = .TRUE.
! MUD1 : mud density
                ELSE IF ( J .EQ. -2 ) THEN
                   IF ( FLH(J) ) THEN
                      CALL W3FLDH (J, NDST, NDSEN, NX, NY, NX, NY,    &
                           TIME0, TIMEN, NH(J), NHMAX, THO, HA, HD, HS,&
                           TTT, XXX, XXX, XXX, TZN, XXX, XXX, MUDD, IERR)
                   ELSE
                      CALL W3FLDG ('READ', IDSTR(J), NDSF(J),         &
                           NDST, NDSEN, NX, NY, NX, NY, TIME0, TIMEN, &
                           TTT, XXX, XXX, XXX, TZN, XXX, XXX, MUDD, IERR)
                   END IF
                   IF ( IERR .LT. 0 )FLLST_ALL(J) = .TRUE.
! MUD2 : mud thickness
                ELSE IF ( J .EQ. -1 ) THEN
                   IF ( FLH(J) ) THEN
                      CALL W3FLDH (J, NDST, NDSEN, NX, NY, NX, NY,    &
                           TIME0, TIMEN, NH(J), NHMAX, THO, HA, HD, HS,&
                           TTT, XXX, XXX, XXX, TTN, XXX, XXX, MUDT, IERR)
                   ELSE
                      CALL W3FLDG ('READ', IDSTR(J), NDSF(J),         &
                           NDST, NDSEN, NX, NY, NX, NY, TIME0, TIMEN, &
                           TTT, XXX, XXX, XXX, TTN, XXX, XXX, MUDT, IERR)
                   END IF
                   IF ( IERR .LT. 0 )FLLST_ALL(J) = .TRUE.
! MUD3 : mud viscosity
                ELSE IF ( J .EQ. 0 ) THEN
                   IF ( FLH(J) ) THEN
                      CALL W3FLDH (J, NDST, NDSEN, NX, NY, NX, NY,    &
                           TIME0, TIMEN, NH(J), NHMAX, THO, HA, HD, HS,&
                           TTT, XXX, XXX, XXX, TVN, XXX, XXX, MUDV, IERR)
                   ELSE
                      CALL W3FLDG ('READ', IDSTR(J), NDSF(J),         &
                           NDST, NDSEN, NX, NY, NX, NY, TIME0, TIMEN, &
                           TTT, XXX, XXX, XXX, TVN, XXX, XXX, MUDV, IERR)
                   END IF
                   IF ( IERR .LT. 0 )FLLST_ALL(J) = .TRUE.
! LEV : water levels
                ELSE IF ( J .EQ. 1 ) THEN
                   IF ( FLH(J) ) THEN
                      CALL W3FLDH (J, NDST, NDSEN, NX, NY, NX, NY,    &
                           TIME0, TIMEN, NH(J), NHMAX, THO, HA, HD, HS,&
                           TTT, XXX, XXX, XXX, TLN, XXX, XXX, WLEV, IERR)
                   ELSE
                      CALL W3FLDG ('READ', IDSTR(J), NDSF(J),         &
                           NDST, NDSEN, NX, NY, NX, NY, TIME0, TIMEN, &
                           TTT, XXX, XXX, XXX, TLN, XXX, XXX, WLEV,   &
                           IERR                                       &
                           )
                   END IF
                   IF ( IERR .LT. 0 ) FLLSTL = .TRUE.
!could use this:   IF ( IERR .LT. 0 ) FLLST_ALL(J) = .TRUE.
! CUR : currents
                ELSE IF ( J .EQ. 2 ) THEN
                   IF ( FLH(J) ) THEN
                      CALL W3FLDH (J, NDST, NDSEN, NX, NY, NX, NY,    &
                           TIME0, TIMEN, NH(J), NHMAX, THO, HA, HD, HS,&
                           TC0, CX0, CY0, XXX, TCN, CXN, CYN, XXX, IERR)
                   ELSE
 
                      CALL W3FLDG ('READ', IDSTR(J), NDSF(J),         &
                           NDST, NDSEN, NX, NY, NX, NY, TIME0, TIMEN, &
                           TC0, CX0, CY0, XXX, TCN, CXN, CYN, XXX,    &
                           IERR                                       &
                           )
                   END IF
! WND : winds
                ELSE IF ( J .EQ. 3 ) THEN
                   IF ( FLH(J) ) THEN
                      CALL W3FLDH (J, NDST, NDSEN, NX, NY, NX, NY,    &
                           TIME0, TIMEN, NH(J), NHMAX, THO, HA, HD, HS,&
                           TW0, WX0, WY0, DT0, TWN, WXN, WYN, DTN, IERR)
                   ELSE
                      CALL W3FLDG ('READ', IDSTR(J), NDSF(J),         &
                           NDST, NDSEN, NX, NY, NX, NY, TIME0, TIMEN, &
                           TW0, WX0, WY0, DT0, TWN, WXN, WYN, DTN,    &
                           IERR                                       &
                           )
                      END IF
! ICE : ice conc.
                ELSE IF ( J .EQ. 4 ) THEN
                   CALL W3FLDG ('READ', IDSTR(J), NDSF(J),            &
                        NDST, NDSEN, NX, NY, NX, NY, TIME0, TIMEN,    &
                        TTT, XXX, XXX, XXX, TIN, XXX, BERGI, ICEI, IERR)
                   IF ( IERR .LT. 0 ) FLLSTI = .TRUE.
!could use this:   IF ( IERR .LT. 0 ) FLLST_ALL(J) = .TRUE.
! Assim data
                ELSE IF ( J .EQ. 5 ) THEN
                   CALL W3FLDD ('SIZE', IDSTR(J), NDSF(J), NDST,      &
                        NDSEN, TIME0, T0N, RCLD(J), NDT(J),           &
                        NDTNEW, DATA0, IERR )
                   IF ( IERR .LT. 0 ) THEN
                        INFLAGS1(J) = .FALSE.
                        IF ( ALLOCATED(DATA0) ) DEALLOCATE(DATA0)
                   ELSE
                        NDT(J) = NDTNEW
                        IF ( ALLOCATED(DATA0) ) DEALLOCATE(DATA0)
                        ALLOCATE ( DATA0(RCLD(J),NDT(J)) )
                        CALL W3FLDD ('READ', IDSTR(J), NDSF(J), NDST, &
                             NDSEN, TIME0, T0N, RCLD(J), NDT(J),      &
                             NDTNEW, DATA0, IERR )
                   END IF
! Assim data
                ELSE IF ( J .EQ. 6 ) THEN
                   CALL W3FLDD ('SIZE', IDSTR(J), NDSF(J), NDST,      &
                        NDSEN, TIME0, T1N, RCLD(J), NDT(J),           &
                        NDTNEW, DATA1, IERR )
                   IF ( IERR .LT. 0 ) THEN
                        INFLAGS1(J) = .FALSE.
                        IF ( ALLOCATED(DATA1) ) DEALLOCATE(DATA1)
                   ELSE
                        NDT(J) = NDTNEW
                        IF ( ALLOCATED(DATA1) ) DEALLOCATE(DATA1)
                        ALLOCATE ( DATA1(RCLD(J),NDT(J)) )
                        CALL W3FLDD ('READ', IDSTR(J), NDSF(J), NDST, &
                             NDSEN, TIME0, T1N, RCLD(J), NDT(J),      &
                             NDTNEW, DATA1, IERR )
                   END IF
! Assim data
                ELSE IF ( J .EQ. 7 ) THEN
                   CALL W3FLDD ('SIZE', IDSTR(J), NDSF(J), NDST,      &
                        NDSEN, TIME0, T2N, RCLD(J), NDT(J),           &
                        NDTNEW, DATA2, IERR )
                   IF ( IERR .LT. 0 ) THEN
                        INFLAGS1(J) = .FALSE.
                        IF ( ALLOCATED(DATA2) ) DEALLOCATE(DATA2)
                   ELSE
                        NDT(J) = NDTNEW
                        IF ( ALLOCATED(DATA2) ) DEALLOCATE(DATA2)
                        ALLOCATE ( DATA2(RCLD(J),NDT(J)) )
                        CALL W3FLDD ('READ', IDSTR(J), NDSF(J), NDST, &
                             NDSEN, TIME0, T2N, RCLD(J), NDT(J),      &
                             NDTNEW, DATA2, IERR )
                   END IF
! Track
                ELSE IF ( J .EQ. 8 ) THEN
                   CALL W3FLDM (4, NDST, NDSEN, TIME0, TIMEN, NH(4),  &
                           NHMAX, THO, HA, HD, TG0, GA0, GD0,         &
                           TGN, GAN, GDN, IERR)
                END IF
!
                IF ( IERR.GT.0 ) GOTO 2222
                IF ( IERR.LT.0 .AND. IAPROC.EQ.NAPOUT )               &
                                 WRITE (NDSO,973) IDFLDS(J)
!
              END IF
!
! 7.a.4 Update next ending time
!
            IF ( INFLAGS1(J) ) THEN
                TTT    = TFN(:,J)
                DTTST  = DSEC21 ( TTT , TTIME )
                IF ( DTTST.GT.0. .AND. .NOT.                          &
                       ( (FLLSTL .AND. J.EQ.1) .OR.                   &
                         (FLLST_ALL(J) .AND. J.EQ.-7) .OR.            &
                         (FLLST_ALL(J) .AND. J.EQ.-6) .OR.            &
                         (FLLST_ALL(J) .AND. J.EQ.-5) .OR.            &
                         (FLLST_ALL(J) .AND. J.EQ.-4) .OR.            &
                         (FLLST_ALL(J) .AND. J.EQ.-3) .OR.            &
                         (FLLST_ALL(J) .AND. J.EQ.-2) .OR.            &
                         (FLLST_ALL(J) .AND. J.EQ.-1) .OR.            &
                         (FLLST_ALL(J) .AND. J.EQ.0 ) .OR.            &
                         (FLLSTI .AND. J.EQ.4) ) ) THEN
                    TTIME  = TTT
! notes: if model has run out beyond field input, then this line should not
!    be reached.
                  END IF
              END IF
!
          END IF
!
        END DO ! J=JFIRST,8
!
! update the next assimilation data time
 
      TDN = TTIME
      CALL TICK21 ( TDN, 1. )
      DO J=5, 7
        IF ( INFLAGS1(J) ) THEN
            TTT    = TFN(:,J)
            DTTST  = DSEC21 ( TTT , TDN )
            IF ( DTTST.GT.0. ) TDN = TTT
          END IF
        END DO
!
      IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,*) ' '
!
! 7.b Run the wave model for the given interval
!
      TIME0  = TTIME
!
      CALL W3WAVE ( 1, TIME0                                          &
                  )
!
      ! The following lines prevents us from trying to read past the end
      ! of the files. This feature existed in v3.14.
      ! "1" is for water levels
      ! "4" is for ice concentration:
      IF ( FLLSTL ) INFLAGS1(1) = .FALSE.
      IF ( FLLSTI ) INFLAGS1(4) = .FALSE.
 
      ! We include something like this for mud and ice parameters also:
      DO J=-7,0
         IF (FLLST_ALL(J))THEN
            INFLAGS1(J)=.FALSE.
         END IF
      END DO
 
!
! 7.c Run data assimilation at ending time
!
      DTTST  = DSEC21 ( TIME , TDN )
      IF ( DTTST .EQ. 0 ) THEN
          CALL STME21 ( TIME0 , DTME21 )
          IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,975) DTME21
!
          FLGDAS(1) = DSEC21(TIME,T0N) .EQ. 0.
          FLGDAS(2) = DSEC21(TIME,T1N) .EQ. 0.
          FLGDAS(3) = DSEC21(TIME,T2N) .EQ. 0.
!
          CALL W3WDAS ( FLGDAS, RCLD, NDT, DATA0, DATA1, DATA2 )
!
! 7.d Call wave model again after data assimilation for output only
!
          DTTST  = DSEC21 ( TIME , TIMEN )
 
          IF ( DTTST .EQ. 0. ) THEN
              IF ( IAPROC .EQ. NAPOUT ) WRITE (NDSO,*) ' '
              CALL W3WAVE ( 1, TIME0                                  &
                          )
            END IF
        END IF
!
! 7.e Check times
!
      DTTST  = DSEC21 ( TIME0 , TIMEN )
      IF ( DTTST .GT. 0. ) GOTO 700
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     End of shell
!
      GOTO 2222
!
! Error escape locations
!
 2000 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1000) IERR
      GOTO 2222
!
 2001 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1001)
      GOTO 2222
!
 2002 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1002) IERR
      GOTO 2222
!
 2003 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1003)
      GOTO 2222
!
 2004 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1004)
      GOTO 2222
!
 2005 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1005) IDTST
      GOTO 2222
!
 2006 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1006) IDTST, NH(J)
      GOTO 2222
!
 2007 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1007)
      GOTO 2222
!
 2008 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1008) IERR
      GOTO 2222
!
 2222 CONTINUE
!
      CALL MPI_BARRIER ( MPI_COMM, IERR_MPI )
!
      IF ( IAPROC .EQ. NAPOUT ) THEN
          CALL DATE_AND_TIME ( VALUES=CLKDT3 )
          CLKFIN = TDIFF ( CLKDT1,CLKDT2 )
          CLKFEL = TDIFF ( CLKDT1,CLKDT3 )
          WRITE (NDSO,997) CLKFIN
          WRITE (NDSO,998) CLKFEL
          IF ( NDSO .NE. NDS(1) ) THEN
              WRITE (NDS(1),997) CLKFIN
              WRITE (NDS(1),998) CLKFEL
            END IF
          WRITE (NDSO,999)
        END IF
!
      CALL MPI_FINALIZE  ( IERR_MPI )
!
! Formats
!
  900 FORMAT (/15X,'      *** WAVEWATCH III Program shell ***      '/ &
               15X,'==============================================='/)
  901 FORMAT ( '  Comment character is ''',A,''''/)
!
  920 FORMAT (/'  Input fields : '/                                   &
               ' --------------------------------------------------')
  921 FORMAT ( '       ',A,2X,A,2X,A)
  922 FORMAT ( ' ' )
!
  930 FORMAT (/'  Time interval : '/                                  &
               ' --------------------------------------------------')
  931 FORMAT ( '       Starting time : ',A)
  932 FORMAT ( '       Ending time   : ',A/)
!
  940 FORMAT (/'  Output requests : '/                                &
               ' --------------------------------------------------'/ &
               '       ',A)
  941 FORMAT (/'       Type',I2,' : ',A/                              &
               '      -----------------------------------------')
  942 FORMAT ( '            From     : ',A)
  943 FORMAT ( '            To       : ',A)
 1944 FORMAT ( '            Interval : ', 8X,A11/)
 2944 FORMAT ( '            Interval : ', 9X,A10/)
 3944 FORMAT ( '            Interval : ',11X,A8/)
 1945 FORMAT ( '            Fields   : ',A)
 2945 FORMAT ( '            Point  1 : ',2F8.2,2X,A)
 2955 FORMAT ( '            Point  1 : ',2(F8.1,'E3'),2X,A)
 2946 FORMAT ( '              ',I6,' : ',2F8.2,2X,A)
 2956 FORMAT ( '              ',I6,' : ',2(F8.1,'E3'),2X,A)
 2947 FORMAT ( '            No points defined')
 3945 FORMAT ( '            The file with ',A,' data is ',A,'.')
 6945 FORMAT ( '            IX first,last,inc :',3I5/                 &
               '            IY first,last,inc :',3I5/                 &
               '            Formatted file    :    ',A)
 8945 FORMAT ( '            Dates out of run dates : output ', A,     &
               ' deactivated')
!
  950 FORMAT (/'  Initializations :'/                                 &
               ' --------------------------------------------------')
  951 FORMAT ( '       ',A)
  952 FORMAT ( '       ',I6,2X,A)
  953 FORMAT ( '          ',I6,I11.8,I7.6,3E12.4)
  954 FORMAT ( '            ',A,': file not needed')
  955 FORMAT ( '            ',A,': file OK')
  956 FORMAT ( '            ',A,': file OK, recl =',I3,               &
               '  undef = ',E10.3)
!
  960 FORMAT (/'  Running model without input fields'/                &
               ' --------------------------------------------------'/)
!
  970 FORMAT (/'  Running model with input fields'/                   &
               ' --------------------------------------------------')
  971 FORMAT (/'  Updating input at ',A)
  972 FORMAT ( '     Updating ',A)
  973 FORMAT ( '        Past last ',A)
  975 FORMAT (/'  Data assimmilation at ',A)
!
  997 FORMAT (/'  Initialization time :',F10.2,' s')
  998 FORMAT ( '  Elapsed time        :',F10.2,' s')
!
  999 FORMAT(//'  End of program '/                                   &
               ' ===================================='/               &
               '         WAVEWATCH III Program shell '/)
!
 1000 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : '/               &
               '     ERROR IN OPENING INPUT FILE'/                    &
               '     IOSTAT =',I5/)
!
 1001 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : '/               &
               '     PREMATURE END OF INPUT FILE'/)
!
 1002 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : '/               &
               '     ERROR IN READING FROM INPUT FILE'/               &
               '     IOSTAT =',I5/)
!
 1003 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : '/               &
               '     ILLEGAL TIME INTERVAL'/)
!
 1004 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : '/               &
               '     POINT OUTPUT ACTIVATED, BUT NO POINTS DEFINED'/)
!
 1005 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : '/               &
               '     ILLEGAL ID STRING HOMOGENEOUS FIELD : ',A/)
!
 1006 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : '/               &
               '     TOO MANY HOMOGENEOUS FIELDS : ',A,1X,I4/)
!
 1007 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : '/               &
               '     INSUFFICIENT DATA FOR HOMOGENEOUS FIELDS'/)
!
 1008 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : '/               &
               '     ERROR IN OPENING OUTPUT FILE'/                   &
               '     IOSTAT =',I5/)
!
!/
!/ End of W3SHEL ----------------------------------------------------- /
!/
     END PROGRAM W3SHEL
