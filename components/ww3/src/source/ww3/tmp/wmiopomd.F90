#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE WMIOPOMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         06-Jun-2012 |
!/                  +-----------------------------------+
!/
!/    09-Aug-2006 : Origination.                        ( version 3.10 )
!/    01-May-2007 : Addd diagnostic output O7a/b.       ( version 3.11 )
!/    21-Jun-2007 : Dedicated output processes.         ( version 3.11 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/    30-Oct-2009 : Implement run-time grid selection.  ( version 3.14 )
!/                  (W. E. Rogers & T. J. Campbell, NRL)
!/    06-Dec-2010 : Change from GLOBAL (logical) to ICLOSE (integer) to
!/                  specify index closure for a grid.   ( version 3.14 )
!/                  (T. J. Campbell, NRL)
!/    06-Mar-2012 : Using MPI_COMM_NULL in checks.      ( version 4.07 )
!/    06-Jun-2012 : Porting bugfixes from 3.14 to 4.07  ( version 4.07 )
!/
!/    Copyright 2009-2012 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Module for generating a single point output file for a multi-
!     grid model implementation.
!
!  2. Variables and types :
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      WMIOPP    Subr  Public   Initialization routine.
!      WMIOPO    Subr  Public   Gather and write routine.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SETG    Subr  W3GDATMD Point to model grid.
!      W3SETW    Subr  W3WDATMD Point to model grid.
!      W3SETA    Subr  W3ADATMD Point to model grid.
!      W3SETO    Subr  W3ODATMD Point to model grid.
!      W3DMO2    Subr     Id.   Dimention model grids output 2.
!      WMSETM    Subr  WMMDATMD Point to model grid.
!      W3MPIP    Subr  W3INITMD Model intiailization.
!      W3IOPP    Sunr  W3IOPOMD Prepare point output for single model.
!      W3IOPO    Sunr     Id.   Point output for single model.
!      W3CSPC    Subr. W3CSPCMD Spectral grid conversion.
!      STRACE    Subr  W3SERVMD Subroutine tracing.
!      EXTCDE    Subr     Id.   Program abort.
!      MPI_SEND, MPI_RECV
!                Subr.  mpif.h  Standard MPI library routines.
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!  6. Switches :
!
!       !/SHRD Distributed memory model.
!       !/MPI
!
!       !O7a   Disgnostic output to NMPSCR.
!       !O7b
!
!       !/S    Enable subroutine tracing.
!       !/T    Enable test output
!       !/MPIT
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
      PUBLIC
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE WMIOPP ( NPT, XPT, YPT, PNAMES )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Sep-2012 !
!/                  +-----------------------------------+
!/
!/    09-Aug-2006 : Origination.                        ( version 3.10 )
!/    01-May-2007 : Addd diagnostic output O7a,b        ( version 3.11 )
!/    30-Oct-2009 : Implement run-time grid selection.  ( version 3.14 )
!/                  (W. E. Rogers & T. J. Campbell, NRL)
!/    06-Dec-2010 : Change from GLOBAL (logical) to ICLOSE (integer) to
!/                  specify index closure for a grid.   ( version 3.14 )
!/                  (T. J. Campbell, NRL)
!/    16-Mar-2012 : Using MPI_COMM_NULL in checks.      ( version 4.07 )
!/    06-Jun-2012 : Porting bugfixes from 3.14 to 4.07  ( version 4.07 )
!/    01-Sep-2012 : Added tests for unstructured grid   ( version 4.07 )
!/                  (M. Dutour Sikiric, IRB & Aron Roland, Z&P)
!/
!  1. Purpose :
!
!     Initialization for unified point output.
!
!  2. Method :
!
!     Find highest resolution grid for each point.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NPT     Int.   I   Number of output points in input.
!       XPT     R.A.   I   X (longitude) coordinates of output points.
!       YPT     R.A.   I   Id. Y.
!       PNAMES  C*10   I   Names of output points.
!     ----------------------------------------------------------------
!       Note: all are optional, and should be given on the first call
!             only, will be taken from storage after that.
!             NPT needs to be ginve always, but can be dummy after
!             first call.
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SETG    Subr  W3GDATMD Point to model grid.
!      W3SETW    Subr  W3WDATMD Point to model grid.
!      W3SETA    Subr  W3ADATMD Point to model grid.
!      W3SETO    Subr  W3ODATMD Point to model grid.
!      W3DMO2    Subr     Id.   Dimension model grids output 2.
!      WMSETM    Subr  WMMDATMD Point to model grid.
!      W3MPIP    Subr  W3INITMD Model intiailization.
!      W3IOPP    Sunr  W3IOPOMD Point output for single model.
!      STRACE    Subr  W3SERVMD Subroutine tracing.
!      EXTCDE    Subr     Id.   Program abort.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      WMINIT    Subr. WMINITMD Wave model initialization routine.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!  7. Remarks :
!
!     - The algorithm used to decide if the pont is in the grid needs
!       to be strictly consistent with W3IOPP.
!     - MPI communication is set up separately from W3MPIO to assure
!       that data are gathered in a single processor even if this
!       procesor is not part of the communicator of the individual
!       model.
!     - In section 2.b the soring of the grids by rand is utilized.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!       !/SHRD Distributed memory model.
!       !/MPI
!
!       !O7a   Disgnostic output to NMPSCR.
!       !O7b
!
!       !/S    Enable subroutine tracing.
!       !/T    Enable test output
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE W3GSRUMD
      USE W3GDATMD, ONLY: W3SETG
      USE W3ADATMD, ONLY: W3SETA
      USE W3WDATMD, ONLY: W3SETW
      USE W3ODATMD, ONLY: W3SETO, W3DMO2
      USE WMMDATMD, ONLY: WMSETM
      USE W3INITMD, ONLY: W3MPIP
      USE W3IOPOMD, ONLY: W3IOPP
      USE W3SERVMD, ONLY: EXTCDE
!
      USE W3GDATMD, ONLY: NX, NY, X0, Y0, SX, MAPSTA, GRIDS,            &
                          FLAGLL, ICLOSE, ICLOSE_NONE, GTYPE, UNGTYPE,  &
                          CLGTYPE, GSU
      USE W3GDATMD, ONLY: XYB, TRIGP, MAXX, MAXY, DXYMAX  ! unstructured grids
      USE W3ODATMD, ONLY: O2INIT, NOPTS, PTLOC, PTNME, GRDID, OUTPTS
      USE W3ODATMD, ONLY: O2IRQI
      USE WMMDATMD, ONLY: MDSE, MDST, NRGRD, MDATAS, IMPROC, NMPSCR,  &
                          NMPERR, MDSS
      USE W3TRIAMD
      USE WMMDATMD, ONLY: MPI_COMM_GRD, MPI_COMM_MWAVE
!
      IMPLICIT NONE
!
      INCLUDE "mpif.h"
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)                    :: NPT
      REAL, INTENT(IN), OPTIONAL             :: XPT(NPT), YPT(NPT)
      CHARACTER(LEN=10),INTENT(IN), OPTIONAL :: PNAMES(NPT)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IPT, J, II
      INTEGER                 :: IX(4), IY(4)         ! created by w3grmp
      REAL                    :: RD(4)                ! created by w3grmp
      INTEGER                 :: itout, I1, I2, I3    ! unstructured grids
      INTEGER                 :: IERR_MPI
      REAL                    :: RX, RY, RDX, RDY
      REAL, PARAMETER         :: ACC = 0.05
      REAL, ALLOCATABLE       :: XP(:), YP(:)
      REAL                    :: FACTOR
      LOGICAL, ALLOCATABLE    :: INGRID(:,:)
      LOGICAL, SAVE           :: SETUP = .FALSE., FLGO7a = .FALSE.
      CHARACTER(LEN=10), ALLOCATABLE :: PN(:)
!/
!
! -------------------------------------------------------------------- /
! 0.  Initializations
!
      CALL W3SETO ( 0, MDSE, MDST )
!
! -------------------------------------------------------------------- /
! 1.  Initialize if necessary and possible
!
      IF ( .NOT. O2INIT ) THEN
!
          IF ( .NOT.PRESENT(XPT) .OR. .NOT.PRESENT(YPT) .OR.          &
               .NOT.PRESENT(PNAMES) ) THEN
              WRITE (MDSE,1000)
              CALL EXTCDE (1)
            END IF
!
          CALL W3DMO2 ( 0, MDSE, MDST, NPT )
!
          NOPTS      = NPT
          PTLOC(1,:) = XPT
          PTLOC(2,:) = YPT
          PTNME      = PNAMES
          GRDID      = 'none'
!
        END IF
!
! -------------------------------------------------------------------- /
! 2.  Locate points in grids
! 2.a Check all points for all grids
!
      IF ( FLAGLL ) THEN
          FACTOR = 1.
        ELSE
          FACTOR = 1.E-3
        END IF
!
      ALLOCATE ( INGRID(NRGRD,NOPTS), XP(NOPTS), YP(NOPTS) )
!
      INGRID = .FALSE.
      XP     = PTLOC(1,:)
      YP     = PTLOC(2,:)
!
      DO J=1, NRGRD
!
        CALL W3SETG ( J, MDSE, MDST )
!
        IF ( GTYPE .EQ. CLGTYPE ) THEN
          IF ( IMPROC .EQ. NMPERR ) WRITE (MDSE,'(/4A)')               &
                ' *** WARNING WMIOPP: ',                               &
                'CURVILINEAR GRID SUPPORT IN THIS ROUTINE ',           &
                'HAS BEEN IMPLEMENTED BUT NOT YET VERIFIED ',          &
                '(curvilinear + unified point output) ***'
          END IF
!
! Loop over output points
!
! notes.....Here, we have pulled coding for UNGTYPE and CLGTYPE from w3iopomd.ftn
! ..........in w3iopomd.ftn, it is "DO IPT=1, NPT" but otherwise very similar
        DO IPT=1, NOPTS
!
!     Check if point within grid
!
          IF (GTYPE .NE. UNGTYPE) THEN
            INGRID(J,IPT) = W3GRMP( GSU, XPT(IPT), YPT(IPT), IX, IY, RD )
            IF ( .NOT.INGRID(J,IPT) ) THEN
              CYCLE
              END IF
          ELSE
            CALL IS_IN_UNGRID(J, XPT(IPT), YPT(IPT), itout, IX, IY, RD )
            IF (itout.eq.0) THEN
              INGRID(J,IPT)=.FALSE.
              END IF
            END IF
!
!     Check if point not on land
!
            IF ( MAPSTA(IY(1),IX(1)) .EQ. 0 .AND. &
              MAPSTA(IY(2),IX(2)) .EQ. 0 .AND. &
              MAPSTA(IY(3),IX(3)) .EQ. 0 .AND. &
              MAPSTA(IY(4),IX(4)) .EQ. 0 ) THEN
              INGRID(J,IPT) = .FALSE.
              CYCLE
              END IF
 
!.........If we've gotten to this point, then we are satisfied that
!................the point is in this grid.
 
        END DO !        DO IPT=1, NOPTS
!
      END DO !      DO J=1, NRGRD
!
      DEALLOCATE ( XP, YP )
!
! 2.b Select a grid for each point
!     start from last, which is supposedly higher resolution
!
      MDATAS(:)%NRUPTS = 0
!
      DO IPT=1, NOPTS
        GRDID(IPT) = '...none...'
        DO J= NRGRD, 1, -1
          IF ( INGRID(J,IPT) ) THEN
            GRDID(IPT) = GRIDS(J)%FILEXT
            MDATAS(J)%NRUPTS = MDATAS(J)%NRUPTS + 1
            EXIT
            END IF
          END DO
        END DO
!
! 2.c Diagnostic output
!
! 2.d Test output
!
      DEALLOCATE ( INGRID )
!
! -------------------------------------------------------------------- /
! 3.  Initialize individual grids
! 3.a Loop over grids
!
      DO J=1, NRGRD
!
! 3.b (De)allocate map arrays
!
        IPT      = MAX ( 1 , MDATAS(J)%NRUPTS )
        IF ( SETUP ) DEALLOCATE ( MDATAS(J)%UPTMAP )
        ALLOCATE ( MDATAS(J)%UPTMAP(IPT) )
!
        IF ( MDATAS(J)%NRUPTS .EQ. 0 ) CYCLE
!
        ALLOCATE ( XP(IPT), YP(IPT), PN(IPT) )
!
! 3.c Set up mapping and point arrays
!
        IPT      = 0
        DO II=1, NOPTS
          IF ( GRDID(II) .NE. GRIDS(J)%FILEXT ) CYCLE
          IPT      = IPT + 1
          MDATAS(J)%UPTMAP(IPT) = II
          XP(IPT)  = PTLOC(1,II)
          YP(IPT)  = PTLOC(2,II)
          PN(IPT)  = PTNME(II)
          END DO
!
        IF ( FLGO7a ) CALL MPI_BARRIER ( MPI_COMM_MWAVE, IERR_MPI )
!
! 3.d Preprocessing for output
!
! 3.d.1 Shared memory version
!
! 3.d.2 Distributed memory version
!
        CALL WMSETM ( J, MDSE, MDST )
!
        IF ( MPI_COMM_GRD .NE. MPI_COMM_NULL ) THEN
!
            CALL W3SETO ( J, MDSE, MDST )
            CALL W3SETG ( J, MDSE, MDST )
            CALL W3SETA ( J, MDSE, MDST )
            CALL W3SETW ( J, MDSE, MDST )
!
            IF ( O2INIT ) THEN
                DEALLOCATE ( OUTPTS(J)%OUT2%IPTINT,               &
                    OUTPTS(J)%OUT2%IL   , OUTPTS(J)%OUT2%IW    ,  &
                    OUTPTS(J)%OUT2%II   , OUTPTS(J)%OUT2%PTIFAC,  &
                    OUTPTS(J)%OUT2%PTNME, OUTPTS(J)%OUT2%GRDID ,  &
                    OUTPTS(J)%OUT2%DPO  , OUTPTS(J)%OUT2%WAO   ,  &
                    OUTPTS(J)%OUT2%WDO  , OUTPTS(J)%OUT2%ASO   ,  &
                    OUTPTS(J)%OUT2%CAO  , OUTPTS(J)%OUT2%CDO   ,  &
                    OUTPTS(J)%OUT2%SPCO , OUTPTS(J)%OUT2%PTLOC )
                O2INIT = .FALSE.
              END IF
!
            CALL W3IOPP ( MDATAS(J)%NRUPTS, XP, YP, PN, J )
!
            IF ( O2IRQI ) THEN
                DEALLOCATE (OUTPTS(J)%OUT2%IRQPO1,                &
                            OUTPTS(J)%OUT2%IRQPO2 )
                O2IRQI = .FALSE.
              END IF
!
            CALL W3MPIP ( J )
!
          END IF
!
! This barrier is needed to straighten out output.
!
! 3.e Reset pointers and clean up
!
        CALL W3SETO ( 0, MDSE, MDST )
        DEALLOCATE ( XP, YP, PN )
!
        END DO
!
      IF ( FLGO7a ) CALL MPI_BARRIER ( MPI_COMM_MWAVE, IERR_MPI )
!
! -------------------------------------------------------------------- /
! 4.  Finalize
!
      SETUP  = .TRUE.
!
      RETURN
!
! Formats
!
 1000 FORMAT (/' *** ERROR WMIOPP : INITALIZATION DATA NOT',          &
               ' AVAILABLE *** '/)
!
!/
!/ End of WMIOPP ----------------------------------------------------- /
!/
      END SUBROUTINE WMIOPP
!/ ------------------------------------------------------------------- /
      SUBROUTINE WMIOPO ( TOUT )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         16-Mar-2012 !
!/                  +-----------------------------------+
!/
!/    09-Aug-2006 : Origination.                        ( version 3.10 )
!/    21-Jun-2007 : Dedicated output processes.         ( version 3.11 )
!/    16-Mar-2012 : Using MPI_COMM_NULL in checks.      ( version 3.14 )
!/
!  1. Purpose :
!
!     Gather and write unified point output.
!
!  2. Method :
!
!     Per-grid point output is already gathered. All data are gathered
!     in the porper storage, and writen using the standard W3IOPO
!     routint from grid number 0.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       TOUT    I.A.   I   Time for output file.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SETG    Subr. W3GDATMD Point to model grid.
!      W3SETW    Subr. W3WDATMD Point to model grid.
!      W3SETO    Subr. W3ODATMD Point to model grid.
!      WMSETM    Subr. WMMDATMD Point to model grid.
!      W3CSPC    Subr. W3CSPCMD Spectral grid conversion.
!      W3IOPO    Subr. W3IOPOMD Point output for single model.
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!      MPI_SEND, MPI_RECV
!                Subr.  mpif.h  Standard MPI library routines.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      WMWAVE    Prog. WMWAVEMD Multi-grid wave model routine.
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
!       !/MPI  Distributed memory model.
!
!       !/S    Enable subroutine tracing.
!       !/T    Enable test output
!       !/MPIT
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!     USE CONSTANTS
!
      USE W3GDATMD, ONLY: W3SETG
      USE W3WDATMD, ONLY: W3SETW
      USE W3ODATMD, ONLY: W3SETO
      USE WMMDATMD, ONLY: WMSETM
      USE W3CSPCMD, ONLY: W3CSPC
      USE W3IOPOMD, ONLY: W3IOPO
!
      USE W3GDATMD, ONLY: NK, NTH, NSPEC, XFR, FR1, TH, SGRDS
      USE W3WDATMD, ONLY: TIME
      USE W3ODATMD, ONLY: IAPROC, NAPROC, NAPPNT, NOPTS, SPCO, DPO,   &
                          WAO, WDO, ASO, CAO, CDO, OUTPTS
      USE WMMDATMD, ONLY: MDST, MDSE, IMPROC, NMPROC, NMPUPT, NRGRD,  &
                          RESPEC, UPTMAP, MDSUP
      USE WMMDATMD, ONLY: MPI_COMM_MWAVE, MPI_COMM_GRD, ALLPRC,  &
                          MTAG0
!
      IMPLICIT NONE
!
      INCLUDE "mpif.h"
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: TOUT(2)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: J, I, II, IT0, IT, ITARG, IFROM
      INTEGER                 :: IERR_MPI, NMPPNT
      INTEGER, ALLOCATABLE    :: STATUS(:,:)
      REAL, POINTER           :: SPEC(:,:)
      REAL, POINTER           :: SPCR(:,:), DPR(:), WAR(:),      &
                                 WDR(:), ASR(:), CAR(:), CDR(:)
!/
!
! -------------------------------------------------------------------- /
! 0.  Initializations
!
      IF ( IMPROC .EQ. NMPUPT ) THEN
          OUTPTS(0)%OUT2%SPCO  = 0.
          OUTPTS(0)%OUT2%DPO   = 1.
          OUTPTS(0)%OUT2%WAO   = 0.
          OUTPTS(0)%OUT2%WDO   = 0.
          OUTPTS(0)%OUT2%ASO   = 0.
          OUTPTS(0)%OUT2%CAO   = 0.
          OUTPTS(0)%OUT2%CDO   = 0.
        END IF
!
! -------------------------------------------------------------------- /
! 1.  Loop over grids for processing local data
!
      DO J=1, NRGRD
!
! 1.a Set up loop
!
        CALL W3SETO ( J, MDSE, MDST )
        CALL W3SETG ( J, MDSE, MDST )
        CALL WMSETM ( J, MDSE, MDST )
!
! 1.b Determine if action
!
        IF ( MPI_COMM_GRD .EQ. MPI_COMM_NULL ) THEN
            CYCLE
          END IF
!
        IF ( NOPTS .EQ. 0 ) THEN
            CYCLE
          END IF
!
        IF ( IAPROC .NE. NAPPNT ) THEN
            CYCLE
          END IF
!
! 1.c Data here, and to remain on present processor.
!
        IF ( IMPROC .EQ. NMPUPT ) THEN
!
! 1.c.1 Spectral conversion if needed
!
            IF ( RESPEC(0,J) ) THEN
                ALLOCATE ( SPEC(SGRDS(0)%NSPEC,NOPTS) )
                CALL W3CSPC ( SPCO, NK, NTH, XFR, FR1, TH(1), SPEC,   &
                     SGRDS(0)%NK, SGRDS(0)%NTH, SGRDS(0)%XFR,         &
                     SGRDS(0)%FR1, SGRDS(0)%TH(1), NOPTS, MDST, MDSE, &
                     SGRDS(0)%FACHFE )
!
! 1.c.2 Spectral conversion not needed
!
              ELSE
                SPEC   => SPCO
              END IF
!
! 1.d Store data at grid 0
!
            DO I=1, NOPTS
              II     = UPTMAP(I)
              OUTPTS(0)%OUT2%SPCO(:,II)  = SPEC(:,I)
              OUTPTS(0)%OUT2%DPO(II)     = DPO(I)
              OUTPTS(0)%OUT2%WAO(II)     = WAO(I)
              OUTPTS(0)%OUT2%WDO(II)     = WDO(I)
              OUTPTS(0)%OUT2%ASO(II)     = ASO(I)
              OUTPTS(0)%OUT2%CAO(II)     = CAO(I)
              OUTPTS(0)%OUT2%CDO(II)     = CDO(I)
              END DO
!
            IF ( RESPEC(0,J) ) DEALLOCATE ( SPEC )
!
! 1.e Data here, and to be sent to other processor.
!
          ELSE
!
            IT0    = MTAG0 - 7*NRGRD - 1
            IT     = IT0 + (J-1)*7
            ITARG  = NMPUPT - 1
!
            IT     = IT + 1
            CALL MPI_SEND ( SPCO(1,1), NSPEC*NOPTS, MPI_REAL,    &
                            ITARG, IT, MPI_COMM_MWAVE, IERR_MPI )
            IT     = IT + 1
            CALL MPI_SEND ( DPO(1), NOPTS, MPI_REAL, ITARG, IT,  &
                            MPI_COMM_MWAVE, IERR_MPI )
            IT     = IT + 1
            CALL MPI_SEND ( WAO(1), NOPTS, MPI_REAL, ITARG, IT,  &
                            MPI_COMM_MWAVE, IERR_MPI )
            IT     = IT + 1
            CALL MPI_SEND ( WDO(1), NOPTS, MPI_REAL, ITARG, IT,  &
                            MPI_COMM_MWAVE, IERR_MPI )
            IT     = IT + 1
            CALL MPI_SEND ( ASO(1), NOPTS, MPI_REAL, ITARG, IT,  &
                            MPI_COMM_MWAVE, IERR_MPI )
            IT     = IT + 1
            CALL MPI_SEND ( CAO(1), NOPTS, MPI_REAL, ITARG, IT,  &
                            MPI_COMM_MWAVE, IERR_MPI )
            IT     = IT + 1
            CALL MPI_SEND ( CDO(1), NOPTS, MPI_REAL, ITARG, IT,  &
                            MPI_COMM_MWAVE, IERR_MPI )
!
          END IF
!
        END DO
!
! -------------------------------------------------------------------- /
! 2.  Check if this is output processor, otherwise exit
!
      IF ( IMPROC .NE. NMPUPT ) THEN
          RETURN
       END IF
!
! -------------------------------------------------------------------- /
! 3.  Loop over grids for processing remote data
!
! 3.a Loop setup
!
      DO J=1, NRGRD
!
        CALL W3SETO ( J, MDSE, MDST )
        CALL W3SETG ( J, MDSE, MDST )
        CALL WMSETM ( J, MDSE, MDST )
!
        DO NMPPNT= NMPROC, 1, -1
          IF ( ALLPRC(NMPPNT,J) .EQ. NAPPNT ) EXIT
          END DO
!
        IF ( NMPPNT.EQ.NMPUPT .OR. NOPTS.EQ.0 ) THEN
            CYCLE
          END IF
!
! 3.b Receive data
!
        IT0    = MTAG0 - 7*NRGRD - 1
        IT     = IT0 + (J-1)*7
        IFROM  = NMPPNT - 1
        ALLOCATE ( SPCR(NSPEC,NOPTS), STATUS(MPI_STATUS_SIZE,1),  &
                   DPR(NOPTS), WAR(NOPTS), WDR(NOPTS), ASR(NOPTS),&
                   CAR(NOPTS), CDR(NOPTS) )
!
        IT     = IT + 1
        CALL MPI_RECV ( SPCR(1,1), NSPEC*NOPTS, MPI_REAL, IFROM,  &
                        IT, MPI_COMM_MWAVE, STATUS, IERR_MPI )
        IT     = IT + 1
        CALL MPI_RECV ( DPR(1), NSPEC*NOPTS, MPI_REAL, IFROM,     &
                        IT, MPI_COMM_MWAVE, STATUS, IERR_MPI )
        IT     = IT + 1
        CALL MPI_RECV ( WAR(1), NSPEC*NOPTS, MPI_REAL, IFROM,     &
                        IT, MPI_COMM_MWAVE, STATUS, IERR_MPI )
        IT     = IT + 1
        CALL MPI_RECV ( WDR(1), NSPEC*NOPTS, MPI_REAL, IFROM,     &
                        IT, MPI_COMM_MWAVE, STATUS, IERR_MPI )
        IT     = IT + 1
        CALL MPI_RECV ( ASR(1), NSPEC*NOPTS, MPI_REAL, IFROM,     &
                        IT, MPI_COMM_MWAVE, STATUS, IERR_MPI )
        IT     = IT + 1
        CALL MPI_RECV ( CAR(1), NSPEC*NOPTS, MPI_REAL, IFROM,     &
                        IT, MPI_COMM_MWAVE, STATUS, IERR_MPI )
        IT     = IT + 1
        CALL MPI_RECV ( CDR(1), NSPEC*NOPTS, MPI_REAL, IFROM,     &
                        IT, MPI_COMM_MWAVE, STATUS, IERR_MPI )
!
! 3.c Convert if necessary
!
        IF ( RESPEC(0,J) ) THEN
            ALLOCATE ( SPEC(SGRDS(0)%NSPEC,NOPTS) )
            CALL W3CSPC ( SPCR, NK, NTH, XFR, FR1, TH(1), SPEC,   &
                 SGRDS(0)%NK, SGRDS(0)%NTH, SGRDS(0)%XFR,         &
                 SGRDS(0)%FR1, SGRDS(0)%TH(1), NOPTS, MDST, MDSE, &
                 SGRDS(0)%FACHFE )
          ELSE
            SPEC   => SPCR
          END IF
!
! 3.d Store data at grid 0
!
        DO I=1, NOPTS
          II     = UPTMAP(I)
          OUTPTS(0)%OUT2%SPCO(:,II)  = SPEC(:,I)
          OUTPTS(0)%OUT2%DPO(II)     = DPR(I)
          OUTPTS(0)%OUT2%WAO(II)     = WAR(I)
          OUTPTS(0)%OUT2%WDO(II)     = WDR(I)
          OUTPTS(0)%OUT2%ASO(II)     = ASR(I)
          OUTPTS(0)%OUT2%CAO(II)     = CAR(I)
          OUTPTS(0)%OUT2%CDO(II)     = CDR(I)
          END DO
!
        IF ( RESPEC(0,J) ) DEALLOCATE ( SPEC )
        DEALLOCATE ( SPCR, DPR, WAR, WDR, ASR, CAR, CDR, STATUS )
!
        END DO
!
! -------------------------------------------------------------------- /
! 4.  Output data
!
      CALL W3SETO ( 0, MDSE, MDST )
      CALL W3SETG ( 0, MDSE, MDST )
      CALL W3SETW ( 0, MDSE, MDST )
!
      TIME   = TOUT
!
      CALL W3IOPO ( 'WRITE', MDSUP, II, 0 )
!
      RETURN
!
! Formats
!
 9040 FORMAT ( ' TEST WMIOPO : PERFORM OUTPUT')
!/
!/ End of WMIOPO ----------------------------------------------------- /
!/
      END SUBROUTINE WMIOPO
!/
!/ End of module WMIOPOMD -------------------------------------------- /
!/
      END MODULE WMIOPOMD
