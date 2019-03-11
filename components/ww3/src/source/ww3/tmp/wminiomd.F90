#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE WMINIOMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         28-Sep-2016 |
!/                  +-----------------------------------+
!/
!/    29-May-2006 : Origination.                        ( version 3.09 )
!/    21-Dec-2006 : VTIME change in WMIOHx and WMIOEx.  ( version 3.10 )
!/    22-Jan-2007 : Adding NAVMAX in WMIOEG.            ( version 3.10 )
!/    30-Jan-2007 : Fix memory leak WMIOBS.             ( version 3.10 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/    28-Sep-2016 : Add error traps for MPI tags.       ( version 5.15 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Internal IO routines for the multi-grid model.
!
!  2. Variables and types :
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      WMIOBS    Subr. Public   Stage internal boundary data.
!      WMIOBG    Subr. Public   Gather internal boundary data.
!      WMIOBF    Subr. Public   Finalize WMIOBS.            ( !/MPI )
!      WMIOHS    Subr. Public   Stage internal high to low data.
!      WMIOHG    Subr. Public   Gather internal high to low data.
!      WMIOHF    Subr. Public   Finalize WMIOHS.            ( !/MPI )
!      WMIOES    Subr. Public   Stage internal same rank data.
!      WMIOEG    Subr. Public   Gather internal same rank data.
!      WMIOEF    Subr. Public   Finalize WMIOES.            ( !/MPI )
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SETG, W3SETW, W3SETA, W3SETO, WMSETM
!                Subr. WxxDATMD Manage data structures.
!      W3UBPT    Subr. W3UBPTMD Update internal bounday spectra.
!      W3IOBC    Subr  W3IOBCMD I/O of boundary data.
!      W3CSPC    Subr. W3CSPCMD Spectral grid conversion.
!      STRACE    Sur.  W3SERVMD Subroutine tracing.
!
!      MPI_ISEND, MPI_IRECV, MPI_TESTALL, MPI_WAITALL
!                Subr.  mpif.h  MPI routines.
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!       !/SHRD   Shared/distributed memory models.
!       !/DIST
!       !/MPI
!
!       !/S      Enable subroutine tracing.
!       !/T      Enable test output
!       !/MPIT
!
!  6. Switches :
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
      PUBLIC
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE WMIOBS ( IMOD )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         28-Sep-2016 !
!/                  +-----------------------------------+
!/
!/    06-Oct-2005 : Origination.                        ( version 3.08 )
!/    29-May-2006 : Adding buffering for MPI.           ( version 3.09 )
!/    30-Jan-2007 : Fix memory leak.                    ( version 3.10 )
!/    28-Sep-2016 : Add error traps for MPI tags.       ( version 5.15 )
!/
!  1. Purpose :
!
!     Stage internal boundary data in the data structure BPSTGE.
!
!  2. Method :
!
!     For the shared memory version, arrays are initialized and the
!     data are copied. For the distributed memory version, the data
!     are moved using a non-blocking send. in this case, the arrays
!     are dimensioned on the recieving side.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.   I   Model number of grid from which data is to
!                          be staged.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SETG, W3SETW, W3SETA, W3SETO, WMSETM
!                Subr. WxxDATMD Manage data structures.
!      W3CSPC    Subr. W3CSPCMD Spectral grid conversion.
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!      EXTCDE    Sur.    Id.    Program abort.
!
!      MPI_ISEND
!                Subr. mpif.h   MPI routines.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      WMINIT    Subr  WMINITMD Multi-grid model initialization.
!      WMWAVE    Subr  WMWAVEMD Multi-grid wave model.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     See FORMAT label 1001.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!       !/SHRD   Shared/distributed memory models.
!       !/DIST
!       !/MPI
!
!       !/S      Enable subroutine tracing.
!       !/T      Enable test output
!       !/MPIT
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE W3GDATMD
      USE W3WDATMD
      USE W3ADATMD
      USE W3ODATMD
      USE WMMDATMD
!
      USE W3CSPCMD, ONLY: W3CSPC
      USE W3SERVMD, ONLY: EXTCDE
!
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
      INTEGER                 :: J, I, IOFF, ISEA, JSEA, IS
      INTEGER                 :: ISPROC
      INTEGER                 :: IP, IT0, ITAG, IERR_MPI
      INTEGER, POINTER        :: NRQ, IRQ(:)
      REAL, POINTER           :: SBPI(:,:), TSTORE(:,:)
!/
!
! -------------------------------------------------------------------- /
! 0.  Initializations
!
      IF ( SUM(NBI2G(:,IMOD)) .EQ. 0 ) RETURN
!
      CALL W3SETO ( IMOD, MDSE, MDST )
      CALL W3SETG ( IMOD, MDSE, MDST )
      CALL W3SETW ( IMOD, MDSE, MDST )
      CALL W3SETA ( IMOD, MDSE, MDST )
!
! -------------------------------------------------------------------- /
! 1.  Loop over grids
!
      DO J=1, NRGRD
!
        IF ( NBI2G(J,IMOD) .EQ. 0 ) CYCLE
!
        CALL WMSETM (   J , MDSE, MDST )
!
        IF ( IMOD .EQ. 1 ) THEN
            IOFF   = 0
          ELSE
            IOFF   = SUM(NBI2G(J,1:IMOD-1))
          END IF
!
! -------------------------------------------------------------------- /
! 2.  Allocate arrays
!
        NAPROC => OUTPTS(J)%NAPROC
        ALLOCATE ( IRQ(NBI2G(J,IMOD)*NAPROC+NAPROC) )
        ALLOCATE ( BPSTGE(J,IMOD)%TSTORE(NSPEC,NBI2G(J,IMOD)) )
        NAPROC => OUTPTS(IMOD)%NAPROC
!
        NRQ    => BPSTGE(J,IMOD)%NRQBPS
        SBPI   => BPSTGE(J,IMOD)%TSTORE
!
        NRQ    = 0
        IRQ    = 0
!
! -------------------------------------------------------------------- /
! 3.  Set the time
!     Note that with MPI the send needs to be posted to the local
!     processor too to make time management possible.
!
        IF ( IAPROC .EQ. 1 ) THEN
            BPSTGE(J,IMOD)%STIME = TIME
            ITAG   = MTAG0 + IMOD + (J-1)*NRGRD
            IF ( ITAG .GT. MTAG1 ) THEN
                WRITE (MDSE,1001)
                CALL EXTCDE (1001)
              END IF
            DO IP=1, NMPROC
              IF ( ALLPRC(IP,J) .NE. 0 .AND.                 &
                   ALLPRC(IP,J) .LE. OUTPTS(J)%NAPROC ) THEN
                  NRQ    = NRQ + 1
                  CALL MPI_ISEND ( BPSTGE(J,IMOD)%STIME, 2,  &
                                   MPI_INTEGER, IP-1, ITAG,  &
                                   MPI_COMM_MWAVE, IRQ(NRQ), &
                                   IERR_MPI )
                END IF
              END DO
          END IF
!
! -------------------------------------------------------------------- /
! 4.  Stage the spectral data
!
        DO I=1, NBI2G(J,IMOD)
!
          ISEA   = NBI2S(IOFF+I,2)
          JSEA   = 1 + (ISEA-1)/NAPROC
          ISPROC = ISEA - (JSEA-1)*NAPROC
          IF ( ISPROC .NE. IAPROC ) CYCLE
          IT0    = MTAG0 + NRGRD**2 + SUM(NBI2G(1:J-1,:)) +      &
                                      SUM(NBI2G(J,1:IMOD-1))
!
          DO IS=1, NSPEC
            SBPI(IS,I) = VA(IS,JSEA) * SIG2(IS) / CG(1+(IS-1)/NTH,ISEA)
            END DO
!
          DO IP=1, NMPROC
            IF ( ALLPRC(IP,J) .NE. 0 .AND.                   &
                 ALLPRC(IP,J) .LE. OUTPTS(J)%NAPROC ) THEN
                NRQ    = NRQ + 1
                ITAG   = IT0 + I
                IF ( ITAG .GT. MTAG1 ) THEN
                    WRITE (MDSE,1001)
                    CALL EXTCDE (1001)
                  END IF
                CALL MPI_ISEND ( SBPI(1,I), NSPEC, MPI_REAL, &
                                 IP-1, ITAG, MPI_COMM_MWAVE, &
                                 IRQ(NRQ), IERR_MPI )
              END IF
            END DO
!
          END DO
!
        IF ( NRQ .GT. 0 ) THEN
            ALLOCATE ( BPSTGE(J,IMOD)%IRQBPS(NRQ) )
            BPSTGE(J,IMOD)%IRQBPS = IRQ(:NRQ)
          ELSE
            DEALLOCATE ( BPSTGE(J,IMOD)%TSTORE )
          END IF
!
        DEALLOCATE ( IRQ )
!
! -------------------------------------------------------------------- /
! 5.  Convert spectra ( !/SHRD only )
!
! ... End of loop over grids
!
        END DO
!
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** ERROR WMIOBS : REQUESTED MPI TAG EXCEEDS', &
                                    ' UPPER BOUND (MTAG1) ***')
!
!/
!/ End of WMIOBS ----------------------------------------------------- /
!/
      END SUBROUTINE WMIOBS
!/ ------------------------------------------------------------------- /
      SUBROUTINE WMIOBG ( IMOD, DONE )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2006 !
!/                  +-----------------------------------+
!/
!/    18-Oct-2005 : Origination.                        ( version 3.08 )
!/    29-May-2006 : Adding buffering for MPI.           ( version 3.09 )
!/
!  1. Purpose :
!
!     Gather internal boundary data for a given model.
!
!  2. Method :
!
!     For the shared memory version, datat are gathered from the data
!     structure BPSTGE. For the distributed memeory version, the
!     gathering of thee data are finished first.
!
!     Gathering of data is triggered by the time stamp of the data
!     that is presently in the storage arrays.
!
!     This routine preempts the data flow normally executed by
!     W3IOBC and W3UBPT, and hence bypasses both routines in W3WAVE.
!
!  2. Method :
!
!     Using storage array BPSTAGE and time stamps.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.   I   Model number of grid from which data is to
!                          be gathered.
!       DONE    Log.   O   Flag for completion of operation (opt).
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SETG, W3SETW, W3SETA, W3SETO, WMSETM
!                Subr. WxxDATMD Manage data structures.
!      W3CSPC    Subr. W3CSPCMD Spectral grid conversion.
!      W3UBPT    Subr. W3UBPTMD Update internal bounday spectra.
!      W3IOBC    Subr  W3IOBCMD I/O of boundary data.
!      STRACE    Sur.  W3SERVMD Subroutine tracing.
!      EXTCDE    Sur.    Id.    Program abort.
!      DSEC21    Func. W3TIMEMD Difference between times.
!
!      MPI_IRECV, MPI_TESTALL, MPI_WAITALL
!                Subr.  mpif.h  MPI routines.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      WMINIT    Subr  WMINITMD Multi-grid model initialization.
!      WMWAVE    Subr  WMWAVEMD Multi-grid wave model.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     See FORMAT labels 1001-1002.
!
!  7. Remarks :
!
!  8. Structure :
!
!  9. Switches :
!
!       !/SHRD   Shared/distributed memory models.
!       !/DIST
!       !/MPI
!
!       !/S      Enable subroutine tracing.
!       !/T      Enable test output
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE W3GDATMD
      USE W3WDATMD
      USE W3ADATMD
      USE W3ODATMD
      USE WMMDATMD
!
      USE W3CSPCMD, ONLY: W3CSPC
      USE W3TIMEMD, ONLY: DSEC21
      USE W3UPDTMD, ONLY: W3UBPT
      USE W3IOBCMD, ONLY: W3IOBC
      USE W3SERVMD, ONLY: EXTCDE
!
      IMPLICIT NONE
!
      INCLUDE "mpif.h"
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)            :: IMOD
      LOGICAL, INTENT(OUT), OPTIONAL :: DONE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: J, I, IOFF, TTEST(2), ITEST
      INTEGER                 :: IERR_MPI, IT0, ITAG, IFROM,     &
                                 ISEA, JSEA, ISPROC
      INTEGER, POINTER        :: VTIME(:)
      INTEGER, POINTER        :: NRQ, IRQ(:)
      INTEGER, ALLOCATABLE    :: STATUS(:,:)
      REAL                    :: DTTST, DT1, DT2, W1, W2
      REAL, POINTER           :: SBPI(:,:)
      REAL, ALLOCATABLE       :: TSTORE(:,:)
      LOGICAL                 :: FLAGOK
!/
!
! -------------------------------------------------------------------- /
! 0.  Initializations
!
      IF ( PRESENT(DONE) ) DONE = .FALSE.
!
      CALL W3SETO ( IMOD, MDSE, MDST )
!
      IF ( IAPROC .GT. NAPROC ) THEN
          IF ( PRESENT(DONE) ) DONE = .TRUE.
          RETURN
        END IF
!
      IF ( SUM(NBI2G(IMOD,:)) .EQ. 0 ) THEN
          IF ( PRESENT(DONE) ) DONE = .TRUE.
          RETURN
        END IF
!
      CALL W3SETG ( IMOD, MDSE, MDST )
      CALL W3SETW ( IMOD, MDSE, MDST )
      CALL W3SETA ( IMOD, MDSE, MDST )
!
      IF ( TBPIN(1) .NE. -1 ) THEN
          IF ( DSEC21(TIME,TBPIN) .GT. 0. ) THEN
              IF ( PRESENT(DONE) ) DONE = .TRUE.
              RETURN
            END IF
        END IF
!
! -------------------------------------------------------------------- /
! 1.  Testing / gathering data in staging arrays
!
! 1.a Shared memory version, test valid times. - - - - - - - - - - - - /
!
! 1.b Distributed memory version - - - - - - - - - - - - - - - - - - - /
!
! 1.b.1 NBISTA = 0
!       Check if staging arrays are initialized.
!       Post the proper receives.
!
      IF ( NBISTA(IMOD) .EQ. 0 ) THEN
!
          NRQ    => MDATAS(IMOD)%NRQBPG
          NRQ    = NRGRD + SUM(NBI2G(IMOD,:))
          ALLOCATE ( MDATAS(IMOD)%IRQBPG(NRQ) )
          IRQ    => MDATAS(IMOD)%IRQBPG
          IRQ    = 0
          NRQ    = 0
!
          DO J=1, NRGRD
            IF ( NBI2G(IMOD,J) .EQ. 0 ) CYCLE
!
! ..... Staging arrays
!
            IF ( BPSTGE(IMOD,J)%INIT ) THEN
                IF ( RESPEC(IMOD,J) ) THEN
                    DEALLOCATE ( BPSTGE(IMOD,J)%SBPI )
                    BPSTGE(IMOD,J)%INIT  = .FALSE.
                  ELSE
                    IF ( SIZE(BPSTGE(IMOD,J)%SBPI(:,1)) .NE.     &
                                             SGRDS(J)%NSPEC .OR. &
                         SIZE(BPSTGE(IMOD,J)%SBPI(1,:)) .NE.     &
                                             NBI2G(IMOD,J) ) THEN
                        IF ( IMPROC .EQ. NMPERR ) WRITE (MDSE,1003)
                        CALL EXTCDE (1003)
                      END IF
                  END IF
              END IF
!
            IF ( .NOT. BPSTGE(IMOD,J)%INIT ) THEN
                NSPEC  => SGRDS(J)%NSPEC
                ALLOCATE (BPSTGE(IMOD,J)%SBPI(NSPEC,NBI2G(IMOD,J)))
                NSPEC  => SGRDS(IMOD)%NSPEC
                BPSTGE(IMOD,J)%INIT  = .TRUE.
              END IF
!
! ..... Check valid time to determine staging.
!
            VTIME  => BPSTGE(IMOD,J)%VTIME
            IF ( VTIME(1) .EQ. -1 ) THEN
                DTTST  = 0.
              ELSE
                DTTST  = DSEC21 ( TIME, VTIME )
              END IF
!
! ..... Post receives for data gather
!
            IF ( DTTST .LE. 0. ) THEN
!
! ..... Time
!
                ITAG   = MTAG0 + J + (IMOD-1)*NRGRD
                IFROM  = MDATAS(J)%CROOT - 1
                NRQ    = NRQ + 1
                CALL MPI_IRECV ( BPSTGE(IMOD,J)%VTIME, 2,        &
                                 MPI_INTEGER, IFROM, ITAG,       &
                                 MPI_COMM_MWAVE, IRQ(NRQ),       &
                                 IERR_MPI )
!
! ..... Spectra
!
                IF ( J .EQ. 1 ) THEN
                    IOFF   = 0
                  ELSE
                    IOFF   = SUM(NBI2G(IMOD,1:J-1))
                  END IF
!
                IT0 = MTAG0 + NRGRD**2 + SUM(NBI2G(1:IMOD-1,:))  &
                                       + SUM(NBI2G(IMOD,1:J-1))
!
                SBPI  => BPSTGE(IMOD,J)%SBPI
!
                NAPROC => OUTPTS(J)%NAPROC
                NSPEC  => SGRDS(J)%NSPEC
                DO I=1, NBI2G(IMOD,J)
                  ISEA   = NBI2S(IOFF+I,2)
                  JSEA   = 1 + (ISEA-1)/NAPROC
                  ISPROC = MDATAS(J)%CROOT - 1 + ISEA -          &
                                             (JSEA-1)*NAPROC
                  NRQ    = NRQ + 1
                  ITAG   = IT0 + I
                  CALL MPI_IRECV ( SBPI(1,I), NSPEC,             &
                                   MPI_REAL, ISPROC-1,           &
                                   ITAG, MPI_COMM_MWAVE,         &
                                   IRQ(NRQ), IERR_MPI )
                  END DO
                NSPEC  => SGRDS(IMOD)%NSPEC
                NAPROC => OUTPTS(IMOD)%NAPROC
!
! ..... End IF for posting receives 1.b.1
!
              END IF
!
! ..... End grid loop J in 1.b.1
!
            END DO
!
! ..... Reset status
!       NOTE: if NBI.EQ.0 all times are already OK, skip to section 2
!
          IF ( NBI .GT. 0 ) THEN
              NBISTA(IMOD) = 1
            END IF
!
! ..... End IF in 1.b.1
!
        END IF
!
! 1.b.2 NBISTA = 1
!       Wait for communication to finish.
!       If DONE defined, check if done, otherwise wait.
!
      IF ( NBISTA(IMOD) .EQ. 1 ) THEN
!
          NRQ    => MDATAS(IMOD)%NRQBPG
          IRQ    => MDATAS(IMOD)%IRQBPG
          ALLOCATE ( STATUS(MPI_STATUS_SIZE,NRQ) )
!
! ..... Test communication if DONE is present, wait otherwise
!
          IF ( PRESENT(DONE) ) THEN
!
              CALL MPI_TESTALL ( NRQ, IRQ, FLAGOK, STATUS,       &
                                 IERR_MPI )
!
            ELSE
!
              CALL MPI_WAITALL ( NRQ, IRQ, STATUS, IERR_MPI )
              FLAGOK = .TRUE.
!
            END IF
!
! ..... Go on based on FLAGOK
!
          IF ( FLAGOK ) THEN
              DEALLOCATE ( STATUS, MDATAS(IMOD)%IRQBPG )
              NRQ    = 0
            ELSE
              RETURN
            END IF
!
          NBISTA(IMOD) = 2
!
! 1.b.3 Convert spectra if needed
!
          DO J=1, NRGRD
!
            IF ( RESPEC(IMOD,J) .AND. NBI2G(IMOD,J).NE.0 ) THEN
!
                NSPEC  => SGRDS(J)%NSPEC
                ALLOCATE ( TSTORE(NSPEC,NBI2G(IMOD,J)))
                NSPEC  => SGRDS(IMOD)%NSPEC
                TSTORE = BPSTGE(IMOD,J)%SBPI
                DEALLOCATE ( BPSTGE(IMOD,J)%SBPI )
                ALLOCATE (BPSTGE(IMOD,J)%SBPI(NSPEC,NBI2G(IMOD,J)))
!
                SBPI   => BPSTGE(IMOD,J)%SBPI
                CALL W3CSPC ( TSTORE, SGRDS(J)%NK, SGRDS(J)%NTH, &
                     SGRDS(J)%XFR, SGRDS(J)%FR1, SGRDS(J)%TH(1), &
                     SBPI, NK, NTH, XFR, FR1, TH(1),             &
                     NBI2G(IMOD,J), MDST, MDSE, SGRDS(IMOD)%FACHFE)
!
                DEALLOCATE ( TSTORE )
!
              END IF
!
            END DO
!
          NBISTA(IMOD) = 0
!
        END IF
!
! -------------------------------------------------------------------- /
! 2.  Update arrays ABPI0/N and data times
!
! 2.a Determine next valid time
!
      TTEST  = -1
      DO J=1, NRGRD
        IF ( NBI2G(IMOD,J) .EQ. 0 ) CYCLE
        VTIME  => BPSTGE(IMOD,J)%VTIME
        IF ( TTEST(1) .EQ. -1 ) THEN
            TTEST  = VTIME
          ELSE
            DTTST  = DSEC21(VTIME,TTEST)
            IF ( DTTST .GT. 0. ) TTEST  = VTIME
          END IF
        END DO
!
! 2.b Shift data
!
      IF ( TBPIN(1) .EQ. -1 ) THEN
          DTTST  = DSEC21(TTEST,TIME)
          IF ( DTTST .NE. 0. ) THEN
              IF ( NMPROC .EQ. NMPERR ) WRITE (MDSE,1002)
              CALL EXTCDE(1002)
            END IF
          ABPI0  = 0.
        ELSE
          TBPI0  = TBPIN
          ABPI0  = ABPIN
        END IF
!
! 2.c Loop over grids for new spectra
!
      DO J=1, NRGRD
!
        IF ( NBI2G(IMOD,J) .EQ. 0 ) CYCLE
        VTIME  => BPSTGE(IMOD,J)%VTIME
        SBPI   => BPSTGE(IMOD,J)%SBPI
!
        IF ( J .EQ. 1 ) THEN
            IOFF   = 0
          ELSE
            IOFF   = SUM(NBI2G(IMOD,1:J-1))
          END IF
!
        IF ( TBPIN(1) .EQ. -1 ) THEN
            W1     = 0.
            W2     = 1.
          ELSE
            DT1    = DSEC21(TBPI0,VTIME)
            DT2    = DSEC21(TBPI0,TTEST)
            W2     = DT2 / DT1
            W1     = 1. - W2
          END IF
!
        ABPIN(:,IOFF+1:IOFF+NBI2G(IMOD,J)) =                          &
                    W1 * ABPI0(:,IOFF+1:IOFF+NBI2G(IMOD,J)) +         &
                    W2 * SBPI(:,1:NBI2G(IMOD,J))
!
        END DO
!
! 2.d New time
!
      TBPIN  = TTEST
!
! -------------------------------------------------------------------- /
! 3.  Dump data to file if requested
!
      IF ( IAPROC.EQ.NAPBPT .AND. BCDUMP(IMOD) ) THEN
          CALL W3IOBC ( 'DUMP', NDS(9), TBPIN, TBPIN, ITEST, IMOD )
        END IF
!
! -------------------------------------------------------------------- /
! 4.  Update arrays BBPI0/N
!
      CALL W3UBPT
!
! -------------------------------------------------------------------- /
! 5.  Successful update
!
      IF ( PRESENT(DONE) ) DONE = .TRUE.
!
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** ERROR WMIOBG : NO DATA IN STAGING ARRAY ***'/    &
               '                    CALL WMIOBS FIRST '/)
 1002 FORMAT (/' *** ERROR WMIOBG : INITIAL DATA NOT AT INITAL ',     &
                                   'MODEL TIME ***'/)
 1003 FORMAT (/' *** ERROR WMIOBG : UNEXPECTED SIZE OF STAGING', &
                                   ' ARRAY ***')
!
!/
!/ End of WMIOBG ----------------------------------------------------- /
!/
      END SUBROUTINE WMIOBG
!/ ------------------------------------------------------------------- /
      SUBROUTINE WMIOBF ( IMOD )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2006 !
!/                  +-----------------------------------+
!/
!/    18-Oct-2005 : Origination.                        ( version 3.08 )
!/    29-May-2006 : Adding buffering for MPI.           ( version 3.09 )
!/
!  1. Purpose :
!
!     Finalize staging of  internal boundary data in the data
!     structure BPSTGE (MPI only).
!
!  2. Method :
!
!     Post appropriate 'wait' functions to assure that the
!     communication has finished.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.   I   Model number of grid from which data has
!                          been staged.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!
!      MPI_WAITALL
!                Subr. mpif.h   MPI routines.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      WMINIT    Subr  WMINITMD Multi-grid model initialization.
!      WMWAVE    Subr  WMWAVEMD Multi-grid wave model.
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
!       !/SHRD   Shared/distributed memory models.
!       !/DIST
!       !/MPI
!
!       !/S      Enable subroutine tracing.
!       !/T      Test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE WMMDATMD
!
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
      INTEGER                 :: J
      INTEGER                 :: IERR_MPI
      INTEGER, POINTER        :: NRQ, IRQ(:)
      INTEGER, ALLOCATABLE    :: STATUS(:,:)
!/
!
! -------------------------------------------------------------------- /
! 0.  Initializations
!
! -------------------------------------------------------------------- /
! 1.  Loop over grids
!
      DO J=1, NRGRD
!
        NRQ    => BPSTGE(J,IMOD)%NRQBPS
!
! 1.a Nothing to finalize
!
        IF ( NRQ .EQ. 0 ) CYCLE
        IRQ    => BPSTGE(J,IMOD)%IRQBPS
!
! 1.b Wait for communication to end
!
        ALLOCATE ( STATUS(MPI_STATUS_SIZE,NRQ) )
        CALL MPI_WAITALL ( NRQ, IRQ, STATUS, IERR_MPI )
        DEALLOCATE ( STATUS )
!
! 1.c Reset arrays and counter
!
        NRQ    = 0
        DEALLOCATE ( BPSTGE(J,IMOD)%IRQBPS ,                     &
                     BPSTGE(J,IMOD)%TSTORE )
!
        END DO
!
      RETURN
!
! Formats
!
!/
!/ End of WMIOBF ----------------------------------------------------- /
!/
      END SUBROUTINE WMIOBF
!/ ------------------------------------------------------------------- /
      SUBROUTINE WMIOHS ( IMOD )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         28-Sep-2016 !
!/                  +-----------------------------------+
!/
!/    27-Jan-2006 : Origination.                        ( version 3.08 )
!/    20-Dec-2006 : Remove VTIME from MPI comm.         ( version 3.10 )
!/    28-Sep-2016 : Add error traps for MPI tags.       ( version 5.15 )
!/
!  1. Purpose :
!
!     Stage internal high-to-low data in the data structure HGSTGE.
!
!  2. Method :
!
!     Directly fill staging arrays in shared memory version, or post
!     the corresponding sends in distributed memory version.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.   I   Model number of grid from which data is to
!                          be staged.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SETG, W3SETW, W3SETA, W3SETO, WMSETM
!                Subr. WxxDATMD Manage data structures.
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!      EXTCDE    Sur.    Id.    Program abort.
!      DSEC21    Func. W3TIMEMD Difference between times.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      WMWAVE    Subr  WMWAVEMD Multi-grid wave model.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     See FORMAT label 1001.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!       !/SHRD   Shared/distributed memory models.
!       !/DIST
!       !/MPI
!
!       !/S      Enable subroutine tracing.
!       !/T      Enable test output
!       !/MPIT
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE W3GDATMD
      USE W3WDATMD
      USE W3ADATMD
      USE W3ODATMD
      USE WMMDATMD
!
      USE W3SERVMD, ONLY: EXTCDE
      USE W3TIMEMD, ONLY: DSEC21
!
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
      INTEGER                 :: J, NR, I, JSEA, ISEA, IS
      INTEGER                 :: ITAG, IP, IT0, IERR_MPI
      INTEGER                 :: I1, I2
      INTEGER, POINTER        :: NRQ, IRQ(:), NRQOUT, OUTDAT(:,:)
      REAL                    :: DTOUTP
      REAL, POINTER           :: SHGH(:,:)
!/
!
! -------------------------------------------------------------------- /
! 0.  Initializations
!
      IF ( .NOT. FLGHG1 ) THEN
          IF ( SUM(HGSTGE(:,IMOD)%NSND) .EQ. 0 ) RETURN
        ELSE
          IF ( SUM(HGSTGE(:,IMOD)%NSN1) .EQ. 0 ) RETURN
        END IF
!
      CALL W3SETO ( IMOD, MDSE, MDST )
      CALL W3SETG ( IMOD, MDSE, MDST )
      CALL W3SETW ( IMOD, MDSE, MDST )
      CALL W3SETA ( IMOD, MDSE, MDST )
!
! -------------------------------------------------------------------- /
! 1.  Loop over grids
!
      DO J=1, NRGRD
!
        IF ( J .EQ. IMOD ) CYCLE
!
        IF ( .NOT. FLGHG1 ) THEN
            NR     = HGSTGE(J,IMOD)%NSND
          ELSE IF ( FLGHG2 ) THEN
            NR     = HGSTGE(J,IMOD)%NSN1
          ELSE
            IF ( TOUTP(1,J) .EQ. -1 ) THEN
                DTOUTP = 1.
              ELSE
                DTOUTP = DSEC21(TIME,TOUTP(:,J))
              END IF
            IF ( DTOUTP .EQ. 0. ) THEN
                NR     = HGSTGE(J,IMOD)%NSND
              ELSE
                NR     = HGSTGE(J,IMOD)%NSN1
              END IF
          END IF
!
        IF ( NR .EQ. 0 ) CYCLE
        IF ( DSEC21(TIME,TSYNC(:,J)) .NE. 0. ) CYCLE
!
! -------------------------------------------------------------------- /
! 2.  Allocate arrays and/or point pointers
!
        ALLOCATE ( HGSTGE(J,IMOD)%TSTORE(NSPEC,NR) )
        SHGH   => HGSTGE(J,IMOD)%TSTORE
!
        ALLOCATE ( HGSTGE(J,IMOD)%IRQHGS(NR) )
        ALLOCATE ( HGSTGE(J,IMOD)%OUTDAT(NR,3) )
!
        NRQ    => HGSTGE(J,IMOD)%NRQHGS
        NRQOUT => HGSTGE(J,IMOD)%NRQOUT
        IRQ    => HGSTGE(J,IMOD)%IRQHGS
        OUTDAT => HGSTGE(J,IMOD)%OUTDAT
        NRQ    = 0
        NRQOUT = 0
        IRQ    = 0
!
! -------------------------------------------------------------------- /
! 3.  Set the time
!     !/SHRD only.
!
! -------------------------------------------------------------------- /
! 4.  Stage the spectral data
!
        IT0    = MTAG1 + 1
!
        DO I=1, NR
!
          JSEA   = HGSTGE(J,IMOD)%ISEND(I,1)
          ISEA   = IAPROC + NAPROC*(JSEA-1)
          IP     = HGSTGE(J,IMOD)%ISEND(I,2)
          I1     = HGSTGE(J,IMOD)%ISEND(I,3)
          I2     = HGSTGE(J,IMOD)%ISEND(I,4)
          ITAG   = HGSTGE(J,IMOD)%ISEND(I,5) + IT0
          IF ( ITAG .GT. MTAG2 ) THEN
              WRITE (MDSE,1001)
              CALL EXTCDE (1001)
            END IF
!
          DO IS=1, NSPEC
            SHGH(  IS,I  ) = VA(IS,JSEA) * SIG2(IS)             &
                                 / CG(1+(IS-1)/NTH,ISEA)
            END DO
!
          IF ( IP .NE. IMPROC ) THEN
              NRQ    = NRQ + 1
              CALL MPI_ISEND ( SHGH(1,I), NSPEC, MPI_REAL, IP-1, &
                       ITAG, MPI_COMM_MWAVE, IRQ(NRQ), IERR_MPI )
            ELSE
              NRQOUT = NRQOUT + 1
              OUTDAT(NRQOUT,1) = I
              OUTDAT(NRQOUT,2) = I2
              OUTDAT(NRQOUT,3) = I1
            END IF
!
          END DO
!
        END DO
!
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** ERROR WMIOHS : REQUESTED MPI TAG EXCEEDS', &
                                    ' UPPER BOUND (MTAG2) ***')
!
!/
!/ End of WMIOHS ----------------------------------------------------- /
!/
      END SUBROUTINE WMIOHS
!/ ------------------------------------------------------------------- /
      SUBROUTINE WMIOHG ( IMOD, DONE )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         20-Dec-2006 !
!/                  +-----------------------------------+
!/
!/    27-Jan-2006 : Origination.                        ( version 3.08 )
!/    20-Dec-2006 : Remove VTIME from MPI comm.         ( version 3.10 )
!/
!  1. Purpose :
!
!     Gather internal high-to-low data for a given model.
!
!  2. Method :
!
!     For distributed memory version first receive all staged data.
!     After staged data is present, average, convert as necessary,
!     and store in basic spatral arrays.
!
!  2. Method :
!
!     Using storage array HGSTAGE and time stamps.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.   I   Model number of grid from which data is to
!                          be gathered.
!       DONE    Log.   O   Flag for completion of operation (opt).
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SETG, W3SETW, W3SETA, W3SETO
!                Subr. WxxDATMD Manage data structures.
!      W3CSPC    Subr. W3CSPCMD Spectral grid conversion.
!      STRACE    Sur.  W3SERVMD Subroutine tracing.
!      DSEC21    Func. W3TIMEMD Difference between times.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      WMWAVE    Subr  WMWAVEMD Multi-grid wave model.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     See FORMAT labels 1001-1002.
!
!  7. Remarks :
!
!  8. Structure :
!
!  9. Switches :
!
!       !/SHRD   Shared/distributed memory models.
!       !/DIST
!       !/MPI
!
!       !/S      Enable subroutine tracing.
!       !/T      Enable test output
!       !/MPIT
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE W3GDATMD
      USE W3WDATMD
      USE W3ADATMD
      USE W3ODATMD
      USE WMMDATMD
!
      USE W3CSPCMD, ONLY: W3CSPC
      USE W3TIMEMD, ONLY: DSEC21
!     USE W3SERVMD, ONLY: EXTCDE
!
      IMPLICIT NONE
!
      INCLUDE "mpif.h"
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)            :: IMOD
      LOGICAL, INTENT(OUT), OPTIONAL :: DONE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: NTOT, J, IS, NA, IA, JSEA, ISEA, I
      INTEGER                 :: ITAG, IT0, IFROM, ILOC, NLOC,   &
                                 ISPROC, IERR_MPI, ICOUNT,       &
                                 I0, I1, I2
      INTEGER, POINTER        :: VTIME(:)
      INTEGER, POINTER        :: NRQ, IRQ(:), STATUS(:,:)
      REAL                    :: DTTST, WGTH
      REAL, POINTER           :: SPEC1(:,:), SPEC2(:,:), SPEC(:,:)
      REAL, POINTER           :: SHGH(:,:,:)
      LOGICAL                 :: FLGALL
      LOGICAL                 :: FLAGOK
!/
!
! -------------------------------------------------------------------- /
! 0.  Initializations
!
      IF ( TOUTP(1,IMOD) .EQ. -1 ) THEN
          DTTST  = 1.
        ELSE
          DTTST  = DSEC21 ( WDATAS(IMOD)%TIME , TOUTP(:,IMOD) )
        END IF
!
      IF ( .NOT. FLGHG1 ) THEN
          FLGALL = .TRUE.
        ELSE IF ( FLGHG2 ) THEN
          FLGALL = .FALSE.
        ELSE IF ( DTTST .EQ. 0. ) THEN
          FLGALL = .TRUE.
        ELSE
          FLGALL = .FALSE.
       END IF
!
      IF ( FLGALL ) THEN
          NTOT   = SUM(HGSTGE(IMOD,:)%NREC)
        ELSE
          NTOT   = SUM(HGSTGE(IMOD,:)%NRC1)
        END IF
!
      IF ( PRESENT(DONE) ) DONE = .FALSE.
!
      IF ( NTOT .EQ. 0 ) THEN
          IF ( PRESENT(DONE) ) DONE = .TRUE.
          RETURN
        END IF
!
      CALL W3SETO ( IMOD, MDSE, MDST )
      CALL W3SETG ( IMOD, MDSE, MDST )
      CALL W3SETW ( IMOD, MDSE, MDST )
      CALL W3SETA ( IMOD, MDSE, MDST )
!
! -------------------------------------------------------------------- /
! 1.  Testing / gathering data in staging arrays
!
! 1.a Shared memory version, test valid times. - - - - - - - - - - - - /
!
! 1.b Distributed memory version - - - - - - - - - - - - - - - - - - - /
!
! 1.b.1 HGHSTA = 0
!       Check if staging arrays are initialized.
!       Post the proper receives.
!
      IF ( HGHSTA(IMOD) .EQ. 0 ) THEN
!
          NRQ    => MDATAS(IMOD)%NRQHGG
          NRQ    = 0
          DO J=1, NRGRD
            IF ( FLGALL ) THEN
                NRQ    = NRQ + HGSTGE(IMOD,J)%NREC *             &
                               HGSTGE(IMOD,J)%NSMX
              ELSE
                NRQ    = NRQ + HGSTGE(IMOD,J)%NRC1 *             &
                               HGSTGE(IMOD,J)%NSMX
              END IF
            END DO
          NRQ    = MAX(1,NRQ)
          ALLOCATE ( IRQ(NRQ) )
          IRQ    = 0
          NRQ    = 0
!
          DO J=1, NRGRD
            IF ( HGSTGE(IMOD,J)%NTOT .EQ. 0 ) CYCLE
!
! ..... Check valid time to determine staging.
!
            VTIME  => HGSTGE(IMOD,J)%VTIME
            IF ( VTIME(1) .EQ. -1 ) THEN
                DTTST  = 1.
              ELSE
                DTTST  = DSEC21 ( TIME, VTIME )
              END IF
!
! ..... Post receives for data gather
!
            IF ( DTTST .NE. 0. ) THEN
!
! ..... Spectra
!
                IT0 = MTAG1 + 1
                SHGH  => HGSTGE(IMOD,J)%SHGH
!
                IF ( FLGALL ) THEN
                    NTOT   = HGSTGE(IMOD,J)%NREC
                  ELSE
                    NTOT   = HGSTGE(IMOD,J)%NRC1
                  END IF
!
                DO I=1, NTOT
                  NLOC   = HGSTGE(IMOD,J)%NRAVG(I)
                  DO ILOC=1, NLOC
                    ISPROC = HGSTGE(IMOD,J)%IMPSRC(I,ILOC)
                    ITAG   = HGSTGE(IMOD,J)%ITAG(I,ILOC) + IT0
                    IF ( ISPROC .NE. IMPROC ) THEN
                        NRQ    = NRQ + 1
                        CALL MPI_IRECV ( SHGH(1,ILOC,I),         &
                             SGRDS(J)%NSPEC, MPI_REAL,           &
                             ISPROC-1, ITAG, MPI_COMM_MWAVE,     &
                             IRQ(NRQ), IERR_MPI )
                      END IF
                    END DO
                  END DO
!
! ..... End IF for posting receives 1.b.1
!
              END IF
!
! ..... End grid loop J in 1.b.1
!
            END DO
!
          ALLOCATE ( MDATAS(IMOD)%IRQHGG(NRQ) )
          MDATAS(IMOD)%IRQHGG = IRQ(1:NRQ)
          DEALLOCATE ( IRQ )
!
! ..... Reset status
!
          IF ( NRQ .GT. 0 ) THEN
              HGHSTA(IMOD) = 1
            END IF
!
! ..... End IF in 1.b.1
!
        END IF
!
! 1.b.2 HGHSTA = 1
!       Wait for communication to finish.
!       If DONE defined, check if done, otherwise wait.
!
      IF ( HGHSTA(IMOD) .EQ. 1 ) THEN
!
          NRQ    => MDATAS(IMOD)%NRQHGG
          IRQ    => MDATAS(IMOD)%IRQHGG
          ALLOCATE ( STATUS(MPI_STATUS_SIZE,NRQ) )
!
! ..... Test communication if DONE is present, wait otherwise
!
          IF ( PRESENT(DONE) ) THEN
!
              CALL MPI_TESTALL ( NRQ, IRQ, FLAGOK, STATUS,       &
                                 IERR_MPI )
!
            ELSE
!
              CALL MPI_WAITALL ( NRQ, IRQ, STATUS, IERR_MPI )
              FLAGOK = .TRUE.
!
            END IF
!
          DEALLOCATE ( STATUS )
!
! ..... Go on based on FLAGOK
!
          IF ( FLAGOK ) THEN
              NRQ    = 0
              DEALLOCATE ( MDATAS(IMOD)%IRQHGG )
            ELSE
              RETURN
            END IF
!
          HGHSTA(IMOD) = 0
!
        END IF
!
! ..... process locally stored data
!
      DO J=1, NRGRD
        HGSTGE(IMOD,J)%VTIME = TIME
        IF ( J .EQ. IMOD ) CYCLE
        DO IS=1, HGSTGE(IMOD,J)%NRQOUT
          I0     = HGSTGE(IMOD,J)%OUTDAT(IS,1)
          I2     = HGSTGE(IMOD,J)%OUTDAT(IS,2)
          I1     = HGSTGE(IMOD,J)%OUTDAT(IS,3)
          HGSTGE(IMOD,J)%SHGH(:,I2,I1) = HGSTGE(IMOD,J)%TSTORE(:,I0)
          END DO
      END DO
!
! -------------------------------------------------------------------- /
! 2.  Data available, process grid by grid
!
! 2.a Loop over grids
!
      DO J=1, NRGRD
!
        IF ( FLGALL ) THEN
            NTOT   = HGSTGE(IMOD,J)%NREC
          ELSE
            NTOT   = HGSTGE(IMOD,J)%NRC1
          END IF
        IF ( NTOT .EQ. 0 ) CYCLE
!
! 2.b Set up temp data structures
!
        IF ( RESPEC(IMOD,J) ) THEN
            ALLOCATE ( SPEC1(SGRDS(J)%NSPEC,NTOT), SPEC2(NSPEC,NTOT) )
            SPEC   => SPEC1
          ELSE
            ALLOCATE ( SPEC2(NSPEC,NTOT) )
            SPEC   => SPEC2
          END IF
!
! 2.c Average spectra to temp storage
!
        DO IS=1, NTOT
          NA     = HGSTGE(IMOD,J)%NRAVG(IS)
          WGTH   = HGSTGE(IMOD,J)%WGTH(IS,1)
          SPEC(:,IS) = WGTH * HGSTGE(IMOD,J)%SHGH(:,1,IS)
          DO IA=2, NA
            WGTH   = HGSTGE(IMOD,J)%WGTH(IS,IA)
            SPEC(:,IS) = SPEC(:,IS) + WGTH*HGSTGE(IMOD,J)%SHGH(:,IA,IS)
            END DO
          END DO
!
! 2.d Convert spectral grid as needed
!
        IF ( RESPEC(IMOD,J) ) THEN
!
            CALL W3CSPC ( SPEC1, SGRDS(J)%NK, SGRDS(J)%NTH,           &
                          SGRDS(J)%XFR, SGRDS(J)%FR1, SGRDS(J)%TH(1), &
                          SPEC2 , NK, NTH, XFR, FR1, TH(1),           &
                          NTOT, MDST, MDSE, FACHFE)
            DEALLOCATE ( SPEC1 )
!
          END IF
!
! 2.e Move spectra to model
!
        DO IS=1, NTOT
          JSEA   = HGSTGE(IMOD,J)%LJSEA(IS)
          ISEA   = IAPROC + NAPROC*(JSEA-1)
          DO I=1, NSPEC
            VA(I,JSEA) = SPEC2(I,IS) / SIG2(I) * CG(1+(I-1)/NTH,ISEA)
            END DO
          END DO
!
        DEALLOCATE ( SPEC2 )
!
        END DO
!
! -------------------------------------------------------------------- /
! 3.  Set flag if reqeusted
!
      IF ( PRESENT(DONE) ) DONE = .TRUE.
!
      RETURN
!
! Formats
!
!/
!/ End of WMIOHG ----------------------------------------------------- /
!/
      END SUBROUTINE WMIOHG
!/ ------------------------------------------------------------------- /
      SUBROUTINE WMIOHF ( IMOD )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         16-Jan-2006 !
!/                  +-----------------------------------+
!/
!/    16-Jan-2006 : Origination.                        ( version 3.08 )
!/
!  1. Purpose :
!
!     Finalize staging of internal high-to-low data in the data
!     structure HGSTGE (MPI only).
!
!  2. Method :
!
!     Post appropriate 'wait' functions to assure that the
!     communication has finished.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.   I   Model number of grid from which data has
!                          been staged.
!     ----------------------------------------------------------------
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
!      WMWAVE    Subr  WMWAVEMD Multi-grid wave model.
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
!       !/SHRD   Shared/distributed memory models.
!       !/DIST
!       !/MPI
!
!       !/S      Enable subroutine tracing.
!       !/T      Test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE WMMDATMD
!
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
      INTEGER                 :: J
      INTEGER                 :: IERR_MPI
      INTEGER, POINTER        :: NRQ, IRQ(:)
      INTEGER, ALLOCATABLE    :: STATUS(:,:)
!/
!
! -------------------------------------------------------------------- /
! 0.  Initializations
!
! -------------------------------------------------------------------- /
! 1.  Loop over grids
!
      DO J=1, NRGRD
!
        NRQ    => HGSTGE(J,IMOD)%NRQHGS
!
! 1.a Nothing to finalize
!
        IF ( NRQ .EQ. 0 ) CYCLE
        IRQ    => HGSTGE(J,IMOD)%IRQHGS
!
! 1.b Wait for communication to end
!
        ALLOCATE ( STATUS(MPI_STATUS_SIZE,NRQ) )
        CALL MPI_WAITALL ( NRQ, IRQ, STATUS, IERR_MPI )
        DEALLOCATE ( STATUS )
!
! 1.c Reset arrays and counter
!
        NRQ    = 0
        DEALLOCATE ( HGSTGE(J,IMOD)%IRQHGS,                      &
                     HGSTGE(J,IMOD)%TSTORE,                      &
                     HGSTGE(J,IMOD)%OUTDAT )
!
        END DO
!
      RETURN
!
! Formats
!
!/
!/ End of WMIOHF ----------------------------------------------------- /
!/
      END SUBROUTINE WMIOHF
!/ ------------------------------------------------------------------- /
      SUBROUTINE WMIOES ( IMOD )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         28-Sep-2016 !
!/                  +-----------------------------------+
!/
!/    25-May-2006 : Origination.                        ( version 3.09 )
!/    21-Dec-2006 : Remove VTIME from MPI comm.         ( version 3.10 )
!/    28-Sep-2016 : Add error traps for MPI tags.       ( version 5.15 )
!/
!  1. Purpose :
!
!     Stage internal same-rank data in the data structure EQSTGE.
!
!  2. Method :
!
!     Directly fill staging arrays in shared memory version, or post
!     the corresponding sends in distributed memory version.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.   I   Model number of grid from which data is to
!                          be staged.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SETG, W3SETW, W3SETA, W3SETO, WMSETM
!                Subr. WxxDATMD Manage data structures.
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!      EXTCDE    Sur.    Id.    Program abort.
!      DSEC21    Func. W3TIMEMD Difference between times.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      WMWAVE    Subr  WMWAVEMD Multi-grid wave model.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     See FORMAT label 1001.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!       !/SHRD   Shared/distributed memory models.
!       !/DIST
!       !/MPI
!
!       !/S      Enable subroutine tracing.
!       !/T      Enable test output
!       !/MPIT
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE W3GDATMD
      USE W3WDATMD
      USE W3ADATMD
      USE W3ODATMD
      USE WMMDATMD
!
      USE W3SERVMD, ONLY: EXTCDE
      USE W3TIMEMD, ONLY: DSEC21
!
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
      INTEGER                 :: J, NR, I, ISEA, JSEA, IS, I1, I2
      INTEGER                 :: IT0, ITAG, IP, IERR_MPI
      INTEGER, POINTER        :: NRQ, IRQ(:), NRQOUT, OUTDAT(:,:)
      REAL, POINTER           :: SEQL(:,:)
!/
!
! -------------------------------------------------------------------- /
! 0.  Initializations
!
      CALL W3SETO ( IMOD, MDSE, MDST )
      CALL W3SETG ( IMOD, MDSE, MDST )
      CALL W3SETW ( IMOD, MDSE, MDST )
      CALL W3SETA ( IMOD, MDSE, MDST )
!
! -------------------------------------------------------------------- /
! 1.  Loop over grids
!
      DO J=1, NRGRD
!
        IF ( J .EQ. IMOD ) CYCLE
        NR     = EQSTGE(J,IMOD)%NSND
!
        IF ( NR .EQ. 0 ) CYCLE
        IF ( DSEC21(TIME,TSYNC(:,J)) .NE. 0. ) STOP
!
! -------------------------------------------------------------------- /
! 2.  Allocate arrays and/or point pointers
!
        ALLOCATE ( EQSTGE(J,IMOD)%TSTORE(NSPEC,NR) )
        SEQL   => EQSTGE(J,IMOD)%TSTORE
!
        ALLOCATE ( EQSTGE(J,IMOD)%IRQEQS(NR)   ,                 &
                   EQSTGE(J,IMOD)%OUTDAT(NR,3) )
!
        NRQ    => EQSTGE(J,IMOD)%NRQEQS
        NRQOUT => EQSTGE(J,IMOD)%NRQOUT
        IRQ    => EQSTGE(J,IMOD)%IRQEQS
        OUTDAT => EQSTGE(J,IMOD)%OUTDAT
        NRQ    = 0
        NRQOUT = 0
        IRQ    = 0
!
! -------------------------------------------------------------------- /
! 3.  Set the time
!     Note that with MPI the send needs to be posted to the local
!     processor too to make time management possible.
!
! -------------------------------------------------------------------- /
! 4.  Stage the spectral data
!
        IT0 = MTAG2 + 1
!
        DO I=1, NR
!
          ISEA   = EQSTGE(J,IMOD)%SIS(I)
          JSEA   = EQSTGE(J,IMOD)%SJS(I)
          I1     = EQSTGE(J,IMOD)%SI1(I)
          I2     = EQSTGE(J,IMOD)%SI2(I)
          IP     = EQSTGE(J,IMOD)%SIP(I)
          ITAG   = EQSTGE(J,IMOD)%STG(I) + IT0
          IF ( ITAG .GT. MTAG_UB ) THEN
              WRITE (MDSE,1001)
              CALL EXTCDE (1001)
            END IF
!
          DO IS=1, NSPEC
            SEQL(  IS,I  ) = VA(IS,JSEA) * SIG2(IS)             &
                                 / CG(1+(IS-1)/NTH,ISEA)
            END DO
!
          IF ( IP .NE. IMPROC ) THEN
              NRQ    = NRQ + 1
              CALL MPI_ISEND ( SEQL(1,I), NSPEC, MPI_REAL, IP-1, &
                       ITAG, MPI_COMM_MWAVE, IRQ(NRQ), IERR_MPI )
            ELSE
              NRQOUT = NRQOUT + 1
              OUTDAT(NRQOUT,1) = I
              OUTDAT(NRQOUT,2) = I1
              OUTDAT(NRQOUT,3) = I2
            END IF
!
          END DO
!
        END DO
!
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** ERROR WMIOES : REQUESTED MPI TAG EXCEEDS', &
                                  ' UPPER BOUND (MTAG_UB) ***')
!
!/
!/
!/ End of WMIOES ----------------------------------------------------- /
!/
      END SUBROUTINE WMIOES
!/ ------------------------------------------------------------------- /
      SUBROUTINE WMIOEG ( IMOD, DONE )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         22-Jan-2007 !
!/                  +-----------------------------------+
!/
!/    25-May-2006 : Origination.                        ( version 3.09 )
!/    21-Dec-2006 : Remove VTIME from MPI comm.         ( version 3.10 )
!/    22-Jan-2007 : Adding NAVMAX.                      ( version 3.10 )
!/
!  1. Purpose :
!
!     Gather internal same-rank data for a given model.
!
!  2. Method :
!
!     For distributed memory version first receive all staged data.
!     After staged data is present, average, convert as necessary,
!     and store in basic spatral arrays.
!
!  2. Method :
!
!     Using storage array EQSTGE and time stamps.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.   I   Model number of grid from which data is to
!                          be gathered.
!       DONE    Log.   O   Flag for completion of operation (opt).
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SETG, W3SETW, W3SETA, W3SETO
!                Subr. WxxDATMD Manage data structures.
!      W3CSPC    Subr. W3CSPCMD Spectral grid conversion.
!      STRACE    Sur.  W3SERVMD Subroutine tracing.
!      DSEC21    Func. W3TIMEMD Difference between times.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      WMWAVE    Subr  WMWAVEMD Multi-grid wave model.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     See FORMAT labels 1001-1002.
!
!  7. Remarks :
!
!  8. Structure :
!
!  9. Switches :
!
!       !/SHRD   Shared/distributed memory models.
!       !/DIST
!       !/MPI
!
!       !/S      Enable subroutine tracing.
!       !/T      Enable test output
!       !/MPIT
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE W3GDATMD
      USE W3WDATMD
      USE W3ADATMD
      USE W3ODATMD
      USE WMMDATMD
!
      USE W3CSPCMD, ONLY: W3CSPC
      USE W3TIMEMD, ONLY: DSEC21
      USE W3SERVMD, ONLY: EXTCDE
!
      IMPLICIT NONE
!
      INCLUDE "mpif.h"
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)            :: IMOD
      LOGICAL, INTENT(OUT), OPTIONAL :: DONE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: J, I, ISEA, JSEA, IA, IS
      INTEGER                 :: IT0, ITAG, IFROM, IERR_MPI,     &
                                 NA, IP, I1, I2
      INTEGER, POINTER        :: VTIME(:)
      INTEGER, POINTER        :: NRQ, IRQ(:), STATUS(:,:)
      REAL                    :: DTTST, WGHT
      REAL, POINTER           :: SPEC1(:,:), SPEC2(:,:), SPEC(:,:)
      REAL, POINTER           :: SEQL(:,:,:)
      LOGICAL                 :: FLAGOK
      LOGICAL                 :: FLAG
!/
!
! -------------------------------------------------------------------- /
! 0.  Initializations
!
      IF ( PRESENT(DONE) ) DONE = .FALSE.
!
      IF ( EQSTGE(IMOD,IMOD)%NREC .EQ. 0 ) THEN
          IF ( PRESENT(DONE) ) DONE = .TRUE.
          RETURN
        END IF
!
      CALL W3SETO ( IMOD, MDSE, MDST )
      CALL W3SETG ( IMOD, MDSE, MDST )
      CALL W3SETW ( IMOD, MDSE, MDST )
      CALL W3SETA ( IMOD, MDSE, MDST )
!
! -------------------------------------------------------------------- /
! 1.  Testing / gathering data in staging arrays
!
! 1.a Shared memory version, test valid times. - - - - - - - - - - - - /
!
! 1.b Distributed memory version - - - - - - - - - - - - - - - - - - - /
!
! 1.b.1 EQLSTA = 0
!       Check if staging arrays are initialized.
!       Post the proper receives.
!
      IF ( EQLSTA(IMOD) .EQ. 0 ) THEN
!
          NRQ    => MDATAS(IMOD)%NRQEQG
          NRQ    = 0
          DO J=1, NRGRD
            IF ( J .EQ. IMOD ) CYCLE
            NRQ    = NRQ + EQSTGE(IMOD,J)%NREC *                 &
                           EQSTGE(IMOD,J)%NAVMAX
            END DO
          ALLOCATE ( IRQ(NRQ) )
          IRQ    = 0
          NRQ    = 0
!
          DO J=1, NRGRD
            IF ( IMOD .EQ. J ) CYCLE
            IF ( EQSTGE(IMOD,J)%NREC .EQ. 0 ) CYCLE
!
! ..... Check valid time to determine staging.
!
            VTIME  => EQSTGE(IMOD,J)%VTIME
            IF ( VTIME(1) .EQ. -1 ) THEN
                DTTST  = 1.
              ELSE
                DTTST  = DSEC21 ( TIME, VTIME )
              END IF
!
! ..... Post receives for data gather
!
            IF ( DTTST .NE. 0. ) THEN
!
! ..... Spectra
!
                IT0 = MTAG2 + 1
                SEQL  => EQSTGE(IMOD,J)%SEQL
!
                DO I=1, EQSTGE(IMOD,J)%NREC
                  JSEA   = EQSTGE(IMOD,J)%JSEA(I)
                  NA     = EQSTGE(IMOD,J)%NAVG(I)
                  DO IA=1, NA
                    IP     = EQSTGE(IMOD,J)%RIP(I,IA)
                    ITAG   = EQSTGE(IMOD,J)%RTG(I,IA) + IT0
                    IF ( IP .NE. IMPROC ) THEN
                        NRQ    = NRQ + 1
                        CALL MPI_IRECV ( SEQL(1,I,IA),           &
                             SGRDS(J)%NSPEC, MPI_REAL,           &
                             IP-1, ITAG, MPI_COMM_MWAVE,         &
                             IRQ(NRQ), IERR_MPI )
                      END IF
                    END DO
                  END DO
!
! ..... End IF for posting receives 1.b.1
!
              END IF
!
! ..... End grid loop J in 1.b.1
!
            END DO
!
          IF ( NRQ .NE. 0 ) THEN
              ALLOCATE ( MDATAS(IMOD)%IRQEQG(NRQ) )
              MDATAS(IMOD)%IRQEQG = IRQ(1:NRQ)
            END IF
!
          DEALLOCATE ( IRQ )
!
! ..... Reset status
!
          IF ( NRQ .GT. 0 ) THEN
              EQLSTA(IMOD) = 1
            END IF
!
! ..... End IF in 1.b.1
!
        END IF
!
! 1.b.2 EQLSTA = 1
!       Wait for communication to finish.
!       If DONE defined, check if done, otherwise wait.
!
      IF ( EQLSTA(IMOD) .EQ. 1 ) THEN
!
          NRQ    => MDATAS(IMOD)%NRQEQG
          IRQ    => MDATAS(IMOD)%IRQEQG
          ALLOCATE ( STATUS(MPI_STATUS_SIZE,NRQ) )
!
! ..... Test communication if DONE is present, wait otherwise
!
          IF ( PRESENT(DONE) ) THEN
!
              CALL MPI_TESTALL ( NRQ, IRQ, FLAGOK, STATUS,       &
                                 IERR_MPI )
!
            ELSE
!
              CALL MPI_WAITALL ( NRQ, IRQ, STATUS, IERR_MPI )
              FLAGOK = .TRUE.
!
            END IF
!
          DEALLOCATE ( STATUS )
!
! ..... Go on based on FLAGOK
!
          IF ( FLAGOK ) THEN
              IF ( NRQ.NE.0 ) DEALLOCATE ( MDATAS(IMOD)%IRQEQG )
              NRQ    = 0
            ELSE
              RETURN
            END IF
!
          EQLSTA(IMOD) = 0
!
        END IF
!
! ..... process locally stored data
!
      DO J=1, NRGRD
        EQSTGE(IMOD,J)%VTIME = TIME
        IF ( J .EQ. IMOD ) CYCLE
        DO IS=1, EQSTGE(IMOD,J)%NRQOUT
          I      = EQSTGE(IMOD,J)%OUTDAT(IS,1)
          I1     = EQSTGE(IMOD,J)%OUTDAT(IS,2)
          I2     = EQSTGE(IMOD,J)%OUTDAT(IS,3)
          EQSTGE(IMOD,J)%SEQL(:,I1,I2) = EQSTGE(IMOD,J)%TSTORE(:,I)
          END DO
      END DO
!
! -------------------------------------------------------------------- /
! 2.  Data available, process grid by grid
!
! 2.a Do 'native' grid IMOD
!
      DO I=1, EQSTGE(IMOD,IMOD)%NREC
        JSEA   = EQSTGE(IMOD,IMOD)%JSEA(I)
        WGHT   = EQSTGE(IMOD,IMOD)%WGHT(I)
        VA(:,JSEA) = WGHT * VA(:,JSEA)
        END DO
!
! 2.b Loop over other grids
!
      DO J=1, NRGRD
        IF ( IMOD.EQ.J .OR. EQSTGE(IMOD,J)%NREC.EQ.0 ) CYCLE
!
! 2.c Average spectra
!
        ALLOCATE ( SPEC1(SGRDS(J)%NSPEC,EQSTGE(IMOD,J)%NREC) )
        SPEC1  = 0.
!
        DO I=1, EQSTGE(IMOD,J)%NREC
          DO IA=1, EQSTGE(IMOD,J)%NAVG(I)
            SPEC1(:,I) = SPEC1(:,I) + EQSTGE(IMOD,J)%SEQL(:,I,IA) *   &
                                       EQSTGE(IMOD,J)%WAVG(I,IA)
            END DO
          END DO
!
! 2.d Convert spectra
!
        IF ( RESPEC(IMOD,J) ) THEN
            ALLOCATE ( SPEC2(NSPEC,EQSTGE(IMOD,J)%NREC) )
!
            CALL W3CSPC ( SPEC1, SGRDS(J)%NK, SGRDS(J)%NTH,           &
                          SGRDS(J)%XFR, SGRDS(J)%FR1, SGRDS(J)%TH(1), &
                          SPEC2 , NK, NTH, XFR, FR1, TH(1),           &
                          EQSTGE(IMOD,J)%NREC, MDST, MDSE, FACHFE)
!
            SPEC   => SPEC2
          ELSE
            SPEC   => SPEC1
          END IF
!
! 2.e Apply to native grid
!
        DO I=1, EQSTGE(IMOD,J)%NREC
          ISEA   = EQSTGE(IMOD,J)%ISEA(I)
          JSEA   = EQSTGE(IMOD,J)%JSEA(I)
          WGHT   = EQSTGE(IMOD,J)%WGHT(I)
          DO IS=1, NSPEC
            VA(IS,JSEA) = VA(IS,JSEA) + WGHT *                        &
                 SPEC(IS,I) / SIG2(IS) * CG(1+(IS-1)/NTH,ISEA)
            END DO
          END DO
!
! 2.f Final clean up
!
        DEALLOCATE ( SPEC1 )
        IF ( RESPEC(IMOD,J) ) DEALLOCATE ( SPEC2 )
        END DO
!
! -------------------------------------------------------------------- /
! 3.  Set flag if requested
!
      IF ( PRESENT(DONE) ) DONE = .TRUE.
!
      RETURN
!
! Formats
!
!/
!/ End of WMIOEG ----------------------------------------------------- /
!/
      END SUBROUTINE WMIOEG
!/ ------------------------------------------------------------------- /
      SUBROUTINE WMIOEF ( IMOD )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         25-May-2006 !
!/                  +-----------------------------------+
!/
!/    25-May-2006 : Origination.                        ( version 3.09 )
!/
!  1. Purpose :
!
!     Finalize staging of internal same-rank data in the data
!     structure EQSTGE (MPI only).
!
!  2. Method :
!
!     Post appropriate 'wait' functions to assure that the
!     communication has finished.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.   I   Model number of grid from which data has
!                          been staged.
!     ----------------------------------------------------------------
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
!      WMWAVE    Subr  WMWAVEMD Multi-grid wave model.
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
!       !/SHRD   Shared/distributed memory models.
!       !/DIST
!       !/MPI
!
!       !/S      Enable subroutine tracing.
!       !/T      Test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE WMMDATMD
!
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
      INTEGER                 :: J
      INTEGER                 :: IERR_MPI
      INTEGER, POINTER        :: NRQ, IRQ(:)
      INTEGER, ALLOCATABLE    :: STATUS(:,:)
!/
!
! -------------------------------------------------------------------- /
! 0.  Initializations
!
! -------------------------------------------------------------------- /
! 1.  Loop over grids
!
      DO J=1, NRGRD
!
        NRQ    => EQSTGE(J,IMOD)%NRQEQS
!
! 1.a Nothing to finalize
!
        IF ( NRQ .EQ. 0 ) CYCLE
        IRQ    => EQSTGE(J,IMOD)%IRQEQS
!
! 1.b Wait for communication to end
!
        ALLOCATE ( STATUS(MPI_STATUS_SIZE,NRQ) )
        CALL MPI_WAITALL ( NRQ, IRQ, STATUS, IERR_MPI )
        DEALLOCATE ( STATUS )
!
! 1.c Reset arrays and counter
!
        DEALLOCATE ( EQSTGE(J,IMOD)%IRQEQS,                      &
                     EQSTGE(J,IMOD)%TSTORE,                      &
                     EQSTGE(J,IMOD)%OUTDAT )
        NRQ    = 0
!
        END DO
!
      RETURN
!
! Formats
!
!/
!/ End of WMIOEF ----------------------------------------------------- /
!/
      END SUBROUTINE WMIOEF
!/
!/ End of module WMINIOMD -------------------------------------------- /
!/
      END MODULE WMINIOMD
