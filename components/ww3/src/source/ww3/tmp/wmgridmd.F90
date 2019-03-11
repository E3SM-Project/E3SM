#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE WMGRIDMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |           W. E. Rogers            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         10-Dec-2014 |
!/                  +-----------------------------------+
!/
!/    28-Dec-2005 : Origination WMGLOW, WMGHGH, WMRSPC. ( version 3.08 )
!/    09-Mar-2006 : Carry land mask in WMGHGH.          ( version 3.09 )
!/    24-Apr-2006 : Origination WMGEQL.                 ( version 3.09 )
!/    25-Jul-2006 : Point output grid in WMRSPC.        ( version 3.10 )
!/    23-Dec-2006 : Adding group test in WMGEQL.        ( version 3.10 )
!/    28-Dec-2006 : Simplify NIT for partial comm.      ( version 3.10 )
!/    22-Jan-2007 : Add saving of NAVMAX in WMGEQL.     ( version 3.10 )
!/    02-Feb-2007 : Setting FLAGST in WMGEQL.           ( version 3.10 )
!/    07-Feb-2007 : Setting FLAGST in WMGHGH.           ( version 3.10 )
!/    15-Feb-2007 : Tweaking MAPODI algorithm in WMGEQL.( version 3.10 )
!/    11-Apr-2008 : Bug fix active edges WMGEQL.        ( version 3.13 )
!/    14-Apr-2008 : Bug fix for global grids WMGEQL.    ( version 3.13 )
!/    26-Mar-2009 : Adding test output !/T9 to WMGLOW.  ( version 3.14 )
!/    20-May-2009 : Linking FLAGST and FLGHG1.          ( version 3.14 )
!/    26-May-2009 : Fix erroneous cyclic upd in WMGHGH. ( version 3.14 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/    30-Oct-2009 : Implement run-time grid selection.  ( version 3.14 )
!/                  (W. E. Rogers & T. J. Campbell, NRL)
!/    06-Dec-2010 : Change from GLOBAL (logical) to ICLOSE (integer) to
!/                  specify index closure for a grid.   ( version 3.14 )
!/                  (T. J. Campbell, NRL)
!/    23-Dec-2010 : Fix HPFAC and HQFAC by including the COS(YGRD)
!/                  factor with DXDP and DXDQ terms.    ( version 3.14 )
!/                  (T. J. Campbell, NRL)
!/    12-Mar-2012 : Use MPI_COMM_NULL in checks.        ( version 3.14 )
!/    06-Jun-2012 : Porting bugfixes from 3.14 to 4.07  ( version 4.07 )
!/    05-Sep-2012 : Implementation of UNGTYPE with SCRIP
!/                       (Mathieu Dutour Sikiric, IRB; Aron Roland, Z&P)
!/    21-Sep-2012 : Modify WMGHGH to support SCRIP remap( version 4.11 )
!/                  write/read capabilities (K. Lind, NRL)
!/    05-Aug-2013 : Change PR2/3 to UQ/UNO in distances.( version 4.12 )
!/    10-Dec-2014 : Add checks for allocate status      ( version 5.04 )
!/
!/    Copyright 2009-2013 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Routines to determine and process grid dependencies in the
!     multi-grid wave model.
!
!  2. Variables and types :
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      WMGLOW    Subr. Public   Dependencies to lower ranked grids.
!      WMGHGH    Subr. Public   Dependencies to higher ranked grids.
!      WMGEQL    Subr. Public   Dependencies to same ranked grids.
!      WMRSPC    Subr. Public   Make map of flags for spectral
!                               conversion between grids.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SETO, W3SETG, W3DMO5, WMSETM
!                Subr. W3xDATMD Manage data structures.
!
!      STRACE    Sur.  W3SERVMD Subroutine tracing.
!      EXTCDE    Subr.   Id.    Program abort.
!
!      MPI_BCAST, MPI_BARRIER
!                Subr. mpif.h   Comunication routines.
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!      - WMGLOW and WMGHGH need to be run in this order to
!        assure proper resolving of cross-dependencies.
!      - WMGLOW and WMGEQL, idem.
!
!  6. Switches :
!
!     !/PRn  propagation scheme.
!     !/UQ   propagation scheme.
!     !/UNO  propagation scheme.
!
!     !/SHRD Distributed memory approach
!     !/DIST
!     !/MPI
!
!     !/O12  Removed boundary points output WMGEQL (central).
!     !/O13  Removed boundary points output WMGEQL (edge).
!
!     !/S    Enable subroutine tracing.
!     !/Tn   Enable test output.
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
!/
!/ Specify default accessibility
!/
      PUBLIC
!/
!/ Module private variable for checking error returns
!/
      INTEGER, PRIVATE        :: ISTAT
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE WMGLOW ( FLRBPI )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |           W. E. Rogers            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         10-Dec-2014 !
!/                  +-----------------------------------+
!/
!/    06-Oct-2005 : Origination.                        ( version 3.08 )
!/    10-Feb-2006 : Add test on grid resolution.        ( version 3.09 )
!/    26-Mar-2009 : Adding test output !/T9.            ( version 3.14 )
!/    30-Oct-2009 : Implement run-time grid selection.  ( version 3.14 )
!/                  (W. E. Rogers & T. J. Campbell, NRL)
!/    22-Dec-2010 : Adapt for use with irregular grids  ( version 3.14 )
!/                  (W. E. Rogers, NRL)
!/    12-Mar-2012 : Use MPI_COMM_NULL in checks.        ( version 4.07 )
!/    06-Jun-2012 : Porting bugfixes from 3.14 to 4.07  ( version 4.07 )
!/    10-Dec-2014 : Add checks for allocate status      ( version 5.04 )
!/
!  1. Purpose :
!
!     Determine relations to lower ranked grids for each grid.
!     On the fly, the oposite relations are also saved.
!
!  2. Method :
!
!     Map active boundary points to lower ranked grids.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       FLRBPI  L.A.   O   Array with flags for external file use.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SETO, W3SETG, W3DMO5
!                Subr. W3xDATMD Manage data structures.
!
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!      EXTCDE    Subr.   Id.    Program abort.
!
!      MPI_BCAST, MPI_BARRIER
!                Subr. mpif.h   Comunication routines.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      WMINIT    Subr  WMINITMD Multi-grid model initialization.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!  7. Remarks :
!
!     - For MPI version it is assumed that NX, NY, NSEA, and NSEAL are
!       properly initialized even if the grid is not run on the local
!       process.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/MPI  Distribbuted memory management.
!
!     !/S    Enable subroutine tracing.
!     !/T    Enable test output.
!     !/T1   Test output for individual boundary points
!     !/T2   Test output cross-reference table
!     !/T9   Test output of map of boundary origine.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE W3SERVMD, ONLY: EXTCDE
!
      USE W3GDATMD
      USE W3ODATMD
      USE W3TRIAMD
      USE WMMDATMD
!
      IMPLICIT NONE
!
      INCLUDE "mpif.h"
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      LOGICAL, INTENT(OUT), OPTIONAL :: FLRBPI(NRGRD)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: I, IBI, IX, IY, JS, J,       &
                                 JTOT, I1, J1, I2, J2
      INTEGER                 :: NXYG, IERR_MPI
      INTEGER, ALLOCATABLE    :: TSTORE(:,:)
      LOGICAL                 :: FLBARR
      REAL                    :: XA, YA
      REAL                    :: FACTOR
      LOGICAL                 :: GRIDD(NRGRD,NRGRD) ! indicates grid-to-grid
                                                    ! dependency
      LOGICAL                 :: RFILE(NRGRD),  FLAGOK
      LOGICAL                 :: INGRID ! indicates whether boundary point
                               ! is in lower rank grid
      INTEGER                 :: IVER(4),JVER(4) ! (I,J) indices of vertices
                               ! of cell (in lower rank grid J) enclosing
                               ! boundary point (in higher rank grid I)
      REAL                    :: RW(4) ! Array of interpolation weights.
      INTEGER                 :: KVER ! counter for 4 vertices
 
      REAL                    :: DX_MIN_GRIDI,DY_MIN_GRIDI,DX_MAX_GRIDI, &
                                 DY_MAX_GRIDI
      REAL                    :: DX_MIN_GRIDJ,DY_MIN_GRIDJ,DX_MAX_GRIDJ, &
                                 DY_MAX_GRIDJ
      INTEGER  :: ITRI, IM1, IM2, IT, JT, ISFIRST, ITOUT, NBRELEVANT
      REAL :: DIST_MIN, DIST_MAX, EDIST
      LOGICAL RESOL_CHECK
!
!/
!
! -------------------------------------------------------------------- /
! 1.  Test grid, Initialize and synchronize grids as needed ( !/MPI )
!
      FLBARR = .FALSE.
!
      DO I=1, NRGRD
!
        IF ( .NOT. GRIDS(I)%GINIT ) THEN
            IF ( IMPROC .EQ. NMPERR ) WRITE (MDSE,1000) I
            CALL EXTCDE ( 1000 )
          END IF
 
        CALL W3SETO ( I, MDSE, MDST )
        CALL W3SETG ( I, MDSE, MDST )
!
        FLBARR = FLBARR .OR. MDATAS(I)%FBCAST
        IF ( MDATAS(I)%FBCAST .AND.                              &
             MDATAS(I)%MPI_COMM_BCT.NE.MPI_COMM_NULL ) THEN
             NXYG   = GRIDS(I)%NX * GRIDS(I)%NY
             CALL MPI_BCAST ( GRIDS(I)%MAPSTA(1,1), NXYG,        &
                              MPI_INTEGER, 0,                    &
                              MDATAS(I)%MPI_COMM_BCT, IERR_MPI )
             CALL MPI_BCAST ( GRIDS(I)%MAPST2(1,1), NXYG,        &
                              MPI_INTEGER, 0,                    &
                              MDATAS(I)%MPI_COMM_BCT, IERR_MPI )
             CALL MPI_BCAST ( GRIDS(I)%MAPFS (1,1), NXYG,        &
                              MPI_INTEGER, 0,                    &
                              MDATAS(I)%MPI_COMM_BCT, IERR_MPI )
             NXYG   = 3*GRIDS(I)%NSEA
             CALL MPI_BCAST ( GRIDS(I)%MAPSF (1,1), NXYG,        &
                              MPI_INTEGER, 0,                    &
                              MDATAS(I)%MPI_COMM_BCT, IERR_MPI )
             CALL MPI_BCAST ( GRIDS(I)%CLATIS(1), NSEA, MPI_REAL, 0,&
                              MDATAS(I)%MPI_COMM_BCT, IERR_MPI )
             CALL MPI_BCAST ( SGRDS(I)%SIG(0), NK+2, MPI_REAL, 0,&
                              MDATAS(I)%MPI_COMM_BCT, IERR_MPI )
           END IF
!
        END DO
!
      IF (FLBARR) CALL MPI_BARRIER (MPI_COMM_MWAVE,IERR_MPI)
!
! -------------------------------------------------------------------- /
! 2.  Process grids
!
      IF ( FLAGLL ) THEN
          FACTOR = 1.
        ELSE
          FACTOR = 1.E-3
        END IF
!
      GRIDD = .FALSE.
      RFILE = .FALSE.
!
      IF ( .NOT. ALLOCATED(NBI2G) ) THEN
          ALLOCATE ( NBI2G(NRGRD,NRGRD), STAT=ISTAT )
          CHECK_ALLOC_STATUS ( ISTAT )
        END IF
      NBI2G  = 0
!
      DO I=1, NRGRD
!
! 2.a Test for input boundary points
!
        IF ( OUTPTS(I)%OUT5%NBI .EQ. 0 ) THEN
            CYCLE
          END IF
!
! 2.b Test for lowest rank
!
        IF ( GRANK(I) .EQ. 1 ) THEN
            RFILE(I) = .TRUE.
            CYCLE
          END IF
!
! 2.c Search for input boundary points
!
 
        IBI    = 0
!
! ... Set up data structure for grid
!
        CALL W3SETO ( I, MDSE, MDST )
        CALL W3SETG ( I, MDSE, MDST )
        CALL W3DMO5 ( I, MDSE, MDST, 1 )
        ALLOCATE ( TSTORE(NBI,0:4), STAT=ISTAT )
        CHECK_ALLOC_STATUS ( ISTAT )
!
! ... Set up loop structure for grid
!
        DO IY=1, NY
          DO IX=1, NX
 
!notes : MAPSTA refers to GRIDS(I)%MAPSTA ...this is set in W3SETG
            IF ( ABS(MAPSTA(IY,IX)) .EQ. 2 ) THEN
                XA     = XGRD(IY,IX) !old code: X0 + REAL(IX-1)*SX
                YA     = YGRD(IY,IX) !old code: Y0 + REAL(IY-1)*SY
!
! ... Loop over previous (lower ranked) grids, going in order from highest
!          of lower ranked grids (I-1) to lowest of lower ranked grids (1)
!
                JS     = 0
!
                DO J=I-1, 1, -1
!
                  IF ( GRANK(J) .GE. GRANK(I) ) CYCLE
!
! ... Check if in grid
 
! notes:
! old version (v4.00):
!        if in grid, return location in grid: a) JX, JY
!                                   (lower left indices of cell),
!                                             b) RX, RY
!                                   (normalized location in cell)
!        in not in grid, cycle (search next grid)
! new version (v4.01):
!        Check if point within grid and compute interpolation weights using GSU
!        in not in grid, cycle (search next grid)
!
        IF (GRIDS(J)%GTYPE .EQ. UNGTYPE) THEN
          CALL IS_IN_UNGRID(J, XA, YA, ITOUT, IVER, JVER, RW)
          IF (ITOUT.EQ.0) THEN
            INGRID=.FALSE.
          ELSE
            INGRID=.TRUE.
            FLAGOK =( ABS(GRIDS(J)%MAPSTA(JVER(1),IVER(1))).GE.1 .OR.  &
                      RW(1).LT.0.05 ) .AND.      &
                    ( ABS(GRIDS(J)%MAPSTA(JVER(2),IVER(2))).GE.1 .OR.  &
                      RW(2).LT.0.05 ) .AND.      &
                    ( ABS(GRIDS(J)%MAPSTA(JVER(3),IVER(3))).GE.1 .OR.  &
                      RW(3)   .LT.0.05 )
          END IF
          NbRelevant=3
        ELSE
          INGRID = W3GRMP( GRIDS(J)%GSU,  XA,  YA, IVER , JVER, RW )
!          Print *, 'J=', J, 'IX=', IX, 'IY=', IY
!          Print *, 'IN=', INGRID, 'XA=', XA, 'YA=', YA
!          Print *, '    1: IVER=', IVER(1), 'JVER=', JVER(1), 'RW=', RW(1)
!          Print *, '    2: IVER=', IVER(2), 'JVER=', JVER(2), 'RW=', RW(2)
!          Print *, '    3: IVER=', IVER(3), 'JVER=', JVER(3), 'RW=', RW(3)
!          Print *, '    4: IVER=', IVER(4), 'JVER=', JVER(4), 'RW=', RW(4)
          IF (INGRID) THEN
            FLAGOK =( ABS(GRIDS(J)%MAPSTA(JVER(1),IVER(1))).GE.1 .OR.  &
                      RW(1).LT.0.05 ) .AND.      &
                    ( ABS(GRIDS(J)%MAPSTA(JVER(2),IVER(2))).GE.1 .OR.  &
                      RW(2).LT.0.05 ) .AND.      &
                    ( ABS(GRIDS(J)%MAPSTA(JVER(4),IVER(4))).GE.1 .OR.  &
                      RW(4)   .LT.0.05 ) .AND.      &
                    ( ABS(GRIDS(J)%MAPSTA(JVER(3),IVER(3))).GE.1 .OR.  &
                      RW(3)   .LT.0.05 )
          END IF
          NbRelevant=4
        END IF
!  internal name=                 GSU XTIN YTIN    IS     JS  RW   (notes)
!  role=out                       in    in   in   out    out out
!  size=                          ---    1    1     4      4   4
!
! notes:
!      - organization of IVER(4),JVER(4),RW(4) as returned by W3GRMP are
!        as follows:
!              Point 1 : lower i , lower j (JY1,JX1)
!              Point 2 : upper i , lower j (JY1,JX2)
!              Point 3 : upper i , upper j (JY2,JX2)
!              Point 4 : lower i , upper j (JY2,JX1)
!              (counter-clockwise starting from lower i, lower j)
!
!  ... if not in grid, warning message and cycle (search next grid)
        IF ( .NOT.INGRID ) THEN
           CYCLE
        END IF
 
!
! ... Check against MAPSTA
!
 
! Notes:
! Old code        | becomes | New code
!-----------------| --------| -------
! (1.-RX)*(1.-RY) | becomes | RW(1)
! RX*(1.-RY)      | becomes | RW(2)
! (1.-RX)*RY      | becomes | RW(4)
! RX*RY           | becomes | RW(3)
! JX1             | becomes | IVER(1)
! JY1             | becomes | JVER(1)
! JX2             | becomes | IVER(3)
! JY2             | becomes | JVER(3)
 
! Notes:
! IVER(1)=IVER(4), IVER(2)=IVER(3)
! JVER(1)=JVER(2), JVER(3)=JVER(4)
 
 ! point 1:
                  FLAGOK = ( ABS(GRIDS(J)%MAPSTA(JVER(1),IVER(1))).GE.1 .OR.  &
                                 RW(1).LT.0.05 ) .AND.      &
 ! point 2:
                           ( ABS(GRIDS(J)%MAPSTA(JVER(2),IVER(2))).GE.1 .OR.  &
                                    RW(2).LT.0.05 ) .AND.      &
 ! point 4:
                           ( ABS(GRIDS(J)%MAPSTA(JVER(4),IVER(4))).GE.1 .OR.  &
                                 RW(4)   .LT.0.05 ) .AND.      &
 ! point 3:
                           ( ABS(GRIDS(J)%MAPSTA(JVER(3),IVER(3))).GE.1 .OR.  &
                                    RW(3)   .LT.0.05 )
!
                  IF ( .NOT.FLAGOK ) CYCLE
!
! ... We found interpolation data !
!
                  JS     = J
                  IBI    = IBI + 1
                  GRIDD(I,JS) = .TRUE.
!
                  XBPI(IBI)  = XA
                  YBPI(IBI)  = YA
                  ISBPI(IBI) = MAPFS(IY,IX)
!
                  TSTORE(IBI, 0) = JS
!
! notes:
!   To maintain perfect consistency with old code, we would make code such that:
!   -  point 1 in GSU goes to point 1 in RDBPI, TSTORE
!   -  point 2 in GSU goes to point 2 in RDBPI, TSTORE
!   -  point 4 in GSU goes to point 3 in RDBPI, TSTORE
!   -  point 3 in GSU goes to point 4 in RDBPI, TSTORE
!   Instead, here, we map point 4 in GSU goes to point 4 in RDBPI, TSTORE, etc.
!   Thus the ordering of RDBPI, TSTORE has changed.
!   I have no reason to believe that the ordering in RDBPI, TSTORE is important.
!   I have gone through test case mww3_test_02 for gridsets a,b,c,d and found
!      no change in result vs v4.00.
 
                  DO KVER=1,4
                    IF (KVER .LE. NbRelevant) THEN
                      IF ( ABS(GRIDS(J)%MAPSTA(JVER(KVER),IVER(KVER))).GE.1 &
                          .AND. RW(KVER)   .GT.0.05 ) THEN
                        RDBPI (IBI,KVER) = RW(KVER)
                        TSTORE(IBI,KVER) = GRIDS(J)%MAPFS(JVER(KVER),IVER(KVER))
                      ELSE
                        RDBPI (IBI,KVER) = 0.
                        TSTORE(IBI,KVER) = 0
                      END IF
                    ELSE
                      RDBPI (IBI,KVER) = 0.
                      TSTORE(IBI,KVER) = 0
                    END IF
 
                  END DO
 
!
! .....normalize weights to give sum(R)=1
                  RDBPI(IBI,:) = RDBPI(IBI,:) / SUM(RDBPI(IBI,:))
!
! Search was successful, so no need to search through other grids, so exit loop
                  EXIT
                END DO ! "DO J=..."
!
                IF ( JS.EQ.0 .AND. IMPROC.EQ.NMPERR )                 &
                                 WRITE (MDSE,1020) I, IX, IY, XA, YA
!
              END IF ! If a boundary point...
 
            END DO ! "DO IX=..."
          END DO ! "DO IY=..."
 
!
! 2.d Error checks
!
        IF ( IBI .EQ. 0 ) THEN
            RFILE(I) = .TRUE.
            IF ( IMPROC .EQ. NMPERR ) WRITE (MDSE,1021)
            DEALLOCATE ( OUTPTS(I)%OUT5%IPBPI, OUTPTS(I)%OUT5%ISBPI,  &
                         OUTPTS(I)%OUT5%XBPI,  OUTPTS(I)%OUT5%YBPI,   &
                         OUTPTS(I)%OUT5%RDBPI, STAT=ISTAT )
            CHECK_DEALLOC_STATUS ( ISTAT )
            CYCLE
          ELSE IF ( IBI .NE. OUTPTS(I)%OUT5%NBI ) THEN
             CALL EXTCDE ( 1020 )
          ENDIF
!
! 2.e Sort spectra by grid, fill IPBPI, and get NBI2 and ....
!
 
        IPBPI  = 0
        NBI2   = 0
!
        DO J=1, NRGRD
          DO I1=1, NBI
            IF ( TSTORE(I1,0) .NE. J ) CYCLE
            DO J1=1, 4
              IF ( TSTORE(I1,J1).NE.0 .AND. IPBPI(I1,J1).EQ.0 ) THEN
                  NBI2         = NBI2 + 1
                  IPBPI(I1,J1) = NBI2
                  DO I2=I1, NBI
                    IF ( TSTORE(I2,0) .NE. J ) CYCLE
                    DO J2=1, 4
                      IF ( TSTORE(I2,J2) .EQ. TSTORE(I1,J1) )         &
                           IPBPI(I2,J2) = NBI2
                      END DO
                    END DO
                END IF
              END DO
            END DO
          END DO
!
! 2.f Set up spectral storage and cross-grid mapping
!
        CALL W3DMO5 ( I, MDSE, MDST, 3 )
!
        ALLOCATE ( MDATAS(I)%NBI2S(NBI2,2), STAT=ISTAT )
        CHECK_ALLOC_STATUS ( ISTAT )
        NBI2S  =>  MDATAS(I)%NBI2S
!
        DO I1=1, NBI
          DO J1=1, 4
            IF ( IPBPI(I1,J1) .NE. 0 ) THEN
                NBI2S(IPBPI(I1,J1),1) = TSTORE(I1,0)
                NBI2S(IPBPI(I1,J1),2) = TSTORE(I1,J1)
              END IF
            END DO
          END DO
!
        DO I1=1, NBI2
          NBI2G(I,NBI2S(I1,1)) = NBI2G(I,NBI2S(I1,1)) + 1
          END DO
!
! 2.g Test output
!
        DEALLOCATE ( TSTORE, STAT=ISTAT )
        CHECK_DEALLOC_STATUS ( ISTAT )
!
        END DO
!
! -------------------------------------------------------------------- /
! 3.  Finalyze grid dependencies in GRDLOW
! 3.a Get size of array and dimension
!
 
! notes:
!    GRIDD(I,J) indicates whether grid I is dependent on lower ranked grid J
!    JS counts the number of grids J that grid I is dependent on
!    GRDLOW is sized to accomodate the grid with the largest JS
 
      JTOT   = 0
      DO I=1, NRGRD
        JS     = 0
        DO J=1, NRGRD
          IF ( GRIDD(I,J) ) JS = JS + 1
          END DO
        JTOT   = MAX ( JTOT , JS )
        END DO
!
      IF ( ALLOCATED(GRDLOW) ) THEN
          DEALLOCATE ( GRDLOW, STAT=ISTAT )
          CHECK_DEALLOC_STATUS ( ISTAT )
        END IF
      ALLOCATE ( GRDLOW(NRGRD,0:JTOT), STAT=ISTAT )
      CHECK_ALLOC_STATUS ( ISTAT )
      GRDLOW = 0
!
! 3.b Fill array
!
      FLAGOK = .TRUE.
!
      DO I=1, NRGRD
        JTOT   = 0
        DO J=1, NRGRD
          IF ( GRIDD(I,J) ) THEN
            JTOT   = JTOT + 1
            GRDLOW(I,JTOT) = J
! ... error checking: catch situation where ranks are inconsistent with
!                     resolution
 
! notes:
! old code:   SXJ=GRIDS(J)%SX
!             SXI=GRIDS(I)%SX
!             SYJ=GRIDS(J)%SY
!             SYI=GRIDS(I)%SY
!             also, old code did not need to check both min and max,
!             since they were the same
! new code:
!       SXI(:,:) ==> GRIDS(I)%HPFAC ! resolution in higher rank grid I
!                    (approximate in case of irregular grids)
!       SYI(:,:) ==> GRIDS(I)%HQFAC ! viz.
!       SXJ(:,:) ==> GRIDS(J)%HPFAC ! resolution in lower rank grid J
!                    (approximate in case of irregular grids)
!       SYJ(:,:) ==> GRIDS(J)%HQFAC ! viz.
 
! notes:
! for irregular grids, we require
!       1) smallest cell in low rank grid is larger than smallest cell
!          in high rank grid
!       2) largest cell in low rank grid is larger than largest cell
!          in high rank grid
! Each dimension (along i/p and j/q axes) is checked separately,
!     making 4 checks total.
! This is strict, and may generate "false positives" in error checking
!     here. In this case, the user may wish to disable this error checking.
! For case of regular grids, we cannot use HPFAC, since it goes to zero
!     at pole. We instead use good ol' SX and SY
 
            IF ( GRIDS(I)%GTYPE .EQ. CLGTYPE ) THEN
               DX_MIN_GRIDI=MINVAL(GRIDS(I)%HPFAC)
               DY_MIN_GRIDI=MINVAL(GRIDS(I)%HQFAC)
               DX_MAX_GRIDI=MAXVAL(GRIDS(I)%HPFAC)
               DY_MAX_GRIDI=MAXVAL(GRIDS(I)%HQFAC)
            ELSEIF ( GRIDS(I)%GTYPE .EQ. RLGTYPE ) THEN
               DX_MIN_GRIDI=GRIDS(I)%SX
               DY_MIN_GRIDI=GRIDS(I)%SY
               DX_MAX_GRIDI=GRIDS(I)%SX
               DY_MAX_GRIDI=GRIDS(I)%SY
            ELSEIF ( GRIDS(I)%GTYPE .EQ. UNGTYPE ) THEN
               ISFIRST=1
               DIST_MAX=0
               DIST_MIN=0
               DO ITRI=1,GRIDS(I)%NTRI
                 DO IT=1,3
                   IF (IT.EQ.3) THEN
                     JT=1
                   ELSE
                     JT=IT+1
                   END IF
                   IM1=GRIDS(I)%TRIGP(ITRI,IT)
                   IM2=GRIDS(I)%TRIGP(ITRI,JT)
                   EDIST=W3DIST(FLAGLL, GRIDS(I)%XYB(IM1,1), &
                         GRIDS(I)%XYB(IM1,2), GRIDS(I)%XYB(IM2,1), &
                         GRIDS(I)%XYB(IM2,2))
                   IF (ISFIRST.EQ.1) THEN
                     DIST_MAX=EDIST
                     DIST_MIN=EDIST
                     ISFIRST=0
                   ELSE
                     IF (EDIST.GT.DIST_MAX) THEN
                       DIST_MAX=EDIST
                     END IF
                     IF (EDIST.LT.DIST_MIN) THEN
                       DIST_MIN=EDIST
                     END IF
                   END IF
                 END DO
               END DO
               DX_MIN_GRIDI=DIST_MIN
               DY_MIN_GRIDI=DIST_MIN
               DX_MAX_GRIDI=DIST_MAX
               DY_MAX_GRIDI=DIST_MAX
            ELSE
               CALL EXTCDE ( 601 )
            END IF
 
            IF ( GRIDS(J)%GTYPE .EQ. CLGTYPE ) THEN
               DX_MIN_GRIDJ=MINVAL(GRIDS(J)%HPFAC)
               DY_MIN_GRIDJ=MINVAL(GRIDS(J)%HQFAC)
               DX_MAX_GRIDJ=MAXVAL(GRIDS(J)%HPFAC)
               DY_MAX_GRIDJ=MAXVAL(GRIDS(J)%HQFAC)
            ELSEIF ( GRIDS(J)%GTYPE .EQ. RLGTYPE ) THEN
               DX_MIN_GRIDJ=GRIDS(J)%SX
               DY_MIN_GRIDJ=GRIDS(J)%SY
               DX_MAX_GRIDJ=GRIDS(J)%SX
               DY_MAX_GRIDJ=GRIDS(J)%SY
            ELSEIF ( GRIDS(J)%GTYPE .EQ. UNGTYPE ) THEN
               ISFIRST=1
               DIST_MAX=0
               DIST_MIN=0
               DO ITRI=1,GRIDS(J)%NTRI
                 DO IT=1,3
                   IF (IT.EQ.3) THEN
                     JT=1
                   ELSE
                     JT=IT+1
                   END IF
                   IM1=GRIDS(J)%TRIGP(ITRI,IT)
                   IM2=GRIDS(J)%TRIGP(ITRI,JT)
                   EDIST=W3DIST(FLAGLL, GRIDS(J)%XYB(IM1,1), &
                         GRIDS(J)%XYB(IM1,2), GRIDS(J)%XYB(IM2,1), &
                         GRIDS(J)%XYB(IM2,2))
                   IF (ISFIRST.EQ.1) THEN
                     DIST_MAX=EDIST
                     DIST_MIN=EDIST
                     ISFIRST=0
                   ELSE
                     IF (EDIST.GT.DIST_MAX) THEN
                       DIST_MAX=EDIST
                     END IF
                     IF (EDIST.LT.DIST_MIN) THEN
                       DIST_MIN=EDIST
                     END IF
                   END IF
                 END DO
               END DO
               DX_MIN_GRIDJ=DIST_MIN
               DY_MIN_GRIDJ=DIST_MIN
               DX_MAX_GRIDJ=DIST_MAX
               DY_MAX_GRIDJ=DIST_MAX
           ELSE
               CALL EXTCDE ( 602 )
           END IF
 
           RESOL_CHECK=.FALSE.
           IF (RESOL_CHECK) THEN
              IF ( DX_MIN_GRIDJ .LT. 0.99*DX_MIN_GRIDI .OR.    &
                   DY_MIN_GRIDJ .LT. 0.99*DY_MIN_GRIDI .OR.    &
                   DX_MAX_GRIDJ .LT. 0.99*DX_MAX_GRIDI .OR.    &
                   DY_MAX_GRIDJ .LT. 0.99*DY_MAX_GRIDI  ) THEN
                 Print *, 'DX_MIN_GRID I=', DX_MIN_GRIDI, ' J=', DX_MIN_GRIDJ
                 Print *, 'DX_MAX_GRID I=', DX_MAX_GRIDI, ' J=', DX_MAX_GRIDJ
                 IF ( IMPROC.EQ.NMPERR ) WRITE (MDSE,1030)  &
                      J, GRANK(J), DX_MIN_GRIDJ, DY_MIN_GRIDJ, &
                      DX_MAX_GRIDJ, DY_MAX_GRIDJ, &
                      I, GRANK(I), DX_MIN_GRIDI, DY_MIN_GRIDI, &
                      DX_MAX_GRIDI, DY_MAX_GRIDI
                 FLAGOK = .FALSE.
              END IF
           END IF
 
          END IF ! IF ( GRIDD(I,J) ) THEN
 
        END DO ! DO J=...
        GRDLOW(I,0) = JTOT
      END DO ! DO I=...
!
      IF ( .NOT. FLAGOK ) CALL EXTCDE ( 1030 )
!
! -------------------------------------------------------------------- /
! 4.  Finalyze grid dependencies in GRDHGH
! 4.a Get size of array and dimension
!
      JTOT   = 0
      DO I=1, NRGRD
        JS     = 0
        DO J=1, NRGRD
          IF ( GRIDD(J,I) ) JS = JS + 1
          END DO
        JTOT   = MAX ( JTOT , JS )
        END DO
!
      IF ( ALLOCATED(GRDHGH) ) THEN
          DEALLOCATE ( GRDHGH, STAT=ISTAT )
          CHECK_DEALLOC_STATUS ( ISTAT )
        END IF
      ALLOCATE ( GRDHGH(NRGRD,0:JTOT), STAT=ISTAT )
      CHECK_ALLOC_STATUS ( ISTAT )
      GRDHGH = 0
!
! 4.b Fill array
!
      DO I=1, NRGRD ! low rank grid
        JTOT   = 0
        DO J=1, NRGRD
          IF ( GRIDD(J,I) ) THEN ! grid j is of higher rank than grid i
                                 ! *and* there is dependency
              JTOT   = JTOT + 1  ! count the number of grids of higher
                                 ! rank than grid i
              GRDHGH(I,JTOT) = J ! save the grid number of the higher rank grid
            END IF
          END DO
        GRDHGH(I,0) = JTOT ! save the count of higher ranked grids
        END DO
!
! -------------------------------------------------------------------- /
! 5.  Export file flags
!
      IF ( PRESENT(FLRBPI) ) FLRBPI = RFILE
!
      RETURN
!
! Formats
!
 1000 FORMAT (/' *** WAVEWATCH III ERROR IN WMGLOW : *** '/        &
               '     GRID NOT INITIALIZED, GRID NR',I4 /)
!
 1020 FORMAT (/' *** WAVEWATCH III ERROR IN WMGLOW : *** '/       &
               '     CANNOT FIND SOURCE FOR BOUNDARY DATA '/          &
               '     GRID, IX, IY, X, Y:',3I6,2E12.4/)
!
 1021 FORMAT (/' *** WAVEWATCH III ERROR IN WMGLOW : *** '/       &
               '     NONE OF BOUNDARY POINTS CAN BE MAPPED'/          &
               '     READING FROM FILE INSTEAD'/)
!
 1030 FORMAT (/' *** WAVEWATCH III ERROR IN WMGLOW : *** '/       &
               '     RANKS AND RESOLUTIONS INCONSISTENT'/             &
               '        GRID',I4,' RANK',I4,' RESOLUTION :',4E10.3/   &
               '        GRID',I4,' RANK',I4,' RESOLUTION :',4E10.3/)
!
 2000   FORMAT (/' *** WAVEWATCH-III WARNING : BOUNDARY POINT'/  &
                 '     NOT FOUND IN LOWER RANK GRID : ',2F10.3/  &
                 '     POINT SKIPPED '/)
!
 2001   FORMAT (/' *** WAVEWATCH-III WARNING : BOUNDARY POINT'/  &
                 '     NOT FOUND IN LOWER RANK GRID : ',2E10.3/  &
                 '     POINT SKIPPED '/)
!
!/
!/ End of WMGLOW ----------------------------------------------------- /
!/
      END SUBROUTINE WMGLOW
 
!/ ------------------------------------------------------------------- /
      SUBROUTINE WMGHGH
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |           W. E. Rogers            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         10-Dec-2014 !
!/                  +-----------------------------------+
!/
!/    28-Dec-2005 : Origination.                        ( version 3.08 )
!/    09-Mar-2006 : Carry over land mask.               ( version 3.09 )
!/    28-Dec-2006 : Simplify NIT for partial comm.      ( version 3.10 )
!/    07-Feb-2007 : Setting FLAGST.                     ( version 3.10 )
!/    20-May-2009 : Linking FLAGST and FLGHG1.          ( version 3.14 )
!/    26-May-2009 : Fix erroneous cyclic updating.      ( version 3.14 )
!/    30-Oct-2009 : Implement run-time grid selection.  ( version 3.14 )
!/                  (W. E. Rogers & T. J. Campbell, NRL)
!/    06-Dec-2010 : Change from GLOBAL (logical) to ICLOSE (integer) to
!/                  specify index closure for a grid.   ( version 3.14 )
!/                  (T. J. Campbell, NRL)
!/    23-Dec-2010 : Fix HPFAC and HQFAC by including the COS(YGRD)
!/                  factor with DXDP and DXDQ terms.    ( version 3.14 )
!/                  (T. J. Campbell, NRL)
!/    07-Jul-2011 : Bug fix for IX bounds with wrapping ( version 3.14+)
!/                  grids (see use of "IDSTLA" below)
!/                  (W. E. Rogers, NRL)
!/    02-Aug-2011 : Adapted for use with irregular      ( version 3.14+)
!/                  grids (W. E. Rogers, NRL)
!/    21-Sep-2012 : Modified to implement SCRIP remap   ( version 4.11 )
!/                  file read and write option
!/                  (K. R. Lind, NRL)
!/    05-Aug-2013 : Change PR2/3 to UQ/UNO in distances.( version 4.12 )
!/    10-Dec-2014 : Add checks for allocate status      ( version 5.04 )
!/
!  1. Purpose :
!
!     Determine relation to higher ranked grids for each grid.
!     Base map set in WMGLOW, supplemental data computed here.
!
!  2. Method :
!
!     Map averaging information for higher ranked grid to lower
!     ranked grid.
!
!  3. Parameters :
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SETO, W3SETG, W3DMO5, WMSETM
!                Subr. W3xDATMD Manage data structures.
!      STRACE    Sur.  W3SERVMD Subroutine tracing.
!      EXTCDE    Sur.    Id.    Program abort.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!  6. Error messages :
!
!  7. Remarks :
!
!     Regarding the map of distances to the boundary :
!      - v4.00 : the map of distances to the boundary was intentionally
!                not an accurate characteristic distance. It was felt that
!                it was more important that it be 'safe' and quick to compute.
!                An iterative method was used to compute distance by starting
!                at boundary and working inwards one grid row layer at a time,
!                incrementing distance by dx etc. until the distance map was
!                filled in. This was characterized as "local increment solution
!                only."
!      - v4.01 : conversion to work with irregular grids. Author could not
!                think of any way to retain "local increment solution" method
!                for situation of irregular grids. Therefore method has been
!                changed to compute accurate distances. New method is also
!                more transparent and simpler with much less code, thus
!                easier to modify or debug. It is expected that this method
!                could be more expensive to compute. Isolated timings were
!                not performed. Since the iteration step has been removed,
!                it is hoped that the expense is at least offset somewhat.
!
!     Regarding method of calculating weights :
!              o If SCRIP software is not compiled into WW3 by user
!                (i.e. if SCRIP switch is not set, then original method
!                (denoted "_OM") will be used.
!              o If SCRIP is activated by user, and all grids are
!                regular and specified in terms of meters (cartesian),
!                then WMGHGH will calculate weights using both methods,
!                and then compare the two, producing an error message
!                if they do not match (built-in regression testing)
!                For more info, see Section 0a below.
!
!     re: Inconsistent RANK vs NBI (warning message) in Section 1.d :
!         This was an error, but has been changed to a warning to allow
!         more flexibility, e.g. having two outer grids with different
!         rank (latter to avoid handling via WMGEQL).
!         Change made July 2016.
!         Old system:
!           * grid rank > 1 and NBI>0 : do computations
!           * grid rank > 1 and NBI=0 : error message
!           * grid rank = 1 and NBI>0 : do nothing
!           * grid rank = 1 and NBI=0 : do nothing
!         New system:
!           * grid rank > 1 and NBI>0 : do computations
!           * grid rank > 1 and NBI=0 : do nothing w/ warning message
!           * grid rank = 1 and NBI>0 : do nothing
!           * grid rank = 1 and NBI=0 : do nothing
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/SHRD Distributed memory approach
!     !/DIST
!
!     !/PRn  propagation scheme.
!
!     !/S    Enable subroutine tracing.
!     !/T    Enable test output.
!     !/T3   Test output for received spectra.
!     !/T4   Test output for sent spectra.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE CONSTANTS
      USE W3SERVMD, ONLY: EXTCDE
      USE W3GSRUMD, ONLY: W3DIST
!
      USE W3GDATMD
      USE W3ODATMD
      USE WMMDATMD
!/
      IMPLICIT NONE
!
      INCLUDE "mpif.h"
!
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
 
! notes re: variable names: During the extension for irregular grids,
!    some variable were renamed to make the code more readable:
!         JX==> ISRC
!         JY==> JSRC
!         IX==> IDST
!         IY==> JDST
!         grid I ==> grid GDST
!         grid J ==> grid GSRC
 
      INTEGER                 :: GDST, IJ, IDST, JDST, GSRC, JJ, IB, ISEA,    &
                                 JSEA, IDSTLA, IDSTHA, JDSTLA, JDSTHA,        &
                                 ISRC, JSRC, ISRCL, ISRCH, JSRCL, JSRCH, NIT, &
                                 NRTOT, NROK, JF, JR, NLMAX, ISPROC, ISPRO2,  &
                                 IREC, ISND, ITMP,ILOC
 
      INTEGER                 :: LTAG0
       INTEGER                 :: IERR_MPI
 
      INTEGER, ALLOCATABLE    :: IDSTL(:), IDSTH(:), JDSTL(:), JDSTH(:),      &
                                 MAPTST(:,:),                                 &
                                 I1(:,:), I2(:,:), I3(:), I4(:),              &
                                 INFLND(:,:)
      INTEGER, ALLOCATABLE    :: NX_BEG(:), NX_END(:)
 INTEGER                 :: IM, NX_REM, TAG, NRQ
 
      INTEGER, ALLOCATABLE    :: TMPINT_OM(:,:),TMPINT(:,:)
      REAL, ALLOCATABLE       :: TMPRL_OM(:,:) ,TMPRL(:,:)
      REAL, ALLOCATABLE       :: BDIST_OM(:)   ,BDIST(:)
      INTEGER                 :: NR0   , NR1   , NR2   , NRL   , NLOC
      INTEGER                 :: NR0_OM, NR1_OM, NR2_OM, NRL_OM, NLOC_OM
 
      INTEGER, ALLOCATABLE    :: LTAG(:)
 
      REAL                    :: FACTOR, STX, STY, STXY, NEWVAL,      &
                                 XL, XH, YL, YH, XA, YA, DXC, JD,     &
                                 WX, WY, WTOT
 
      LOGICAL                 :: FLGREC
 
      LOGICAL, ALLOCATABLE    :: GRIDOK(:),                           &
                                 STMASK(:,:), MASKI(:,:), TMPLOG(:)
 
      INTEGER                 :: JBND,IBND ! counter for boundary points
      REAL                    :: DD  ! distance to boundary point
                                     ! (temporary variable)
      REAL                    :: XDST,YDST
      REAL                    :: XSRC,YSRC
      REAL                    :: WXWY
      INTEGER                 :: NJDST,NIDST,KDST
      INTEGER                 :: NJSRC,NISRC,KSRC
      INTEGER                 :: IPNT,ICOUNT,IPNT2
      INTEGER                 :: DST_GRID_SIZE,ISTOP,JTMP
 
      REAL                    :: DX_MAX_GDST,DY_MAX_GDST
      REAL                    :: DX_MIN_GSRC,DY_MIN_GSRC
 
 
      LOGICAL                 :: LSCRIP=.FALSE.      ! true if SCRIP switch is set,
                                                     ! indicates that SCRIP code has
                                                     ! been compiled into WW3
      LOGICAL                 :: LSCRIPNC=.FALSE.    ! true if SCRIPNC switch is set,
                                                     ! indicates that SCRIP code has
                                                     ! been compiled with netCDF
                                                     ! into WW3
      LOGICAL                 :: T38=.FALSE.         ! true if T38 switch is set.
                                                     ! This logical is necessary
                                                     ! since it isn't possible to
                                                     ! have two switches disabling
                                                     ! the same line of code.
      LOGICAL                 :: ALL_REGULAR=.TRUE.  ! true if all grids are
                                                     ! regular grids
      LOGICAL                 :: DO_CHECKING=.FALSE. ! true if we will be
                                                     ! checking "old method" of
                                                     ! computing weights vs.
                                                     ! SCRIP method of computing
                                                     ! weights.
      LOGICAL                 :: OLD_METHOD=.FALSE.  ! true if we will compute
                                                     ! using "old method" (does
                                                     ! not necessarily mean
                                                     ! that this solution is
                                                     ! utilized)
      LOGICAL                 :: LMPIBDI=.FALSE.     ! true if MPIBDI switch is set
      LOGICAL                 :: CALLED_SCRIP=.FALSE.! true if SCRIP has been
                                                     ! called for this processor
 
 
      INTEGER  :: ITRI, IM1, IM2, IT, JT, IsFirst
      REAL :: DIST_MIN, DIST_MAX, eDist
 
!/
 
!
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR_MPI)
 
 
! -------------------------------------------------------------------- /
! 0.  Initializations / tests
!
      IF ( .NOT. ALLOCATED(GRDHGH) ) THEN
         IF ( IMPROC.EQ.NMPERR ) WRITE(MDSE,1000)
         CALL EXTCDE (1000)
      END IF
 
!
!KRL  Allocate helper arrays to enable bottleneck loop parallelization
      ALLOCATE ( NX_BEG(NMPROC), NX_END(NMPROC), STAT=ISTAT )
      CHECK_ALLOC_STATUS ( ISTAT )
!
!!HT:
!!HT: Set up and initialize storage data structures ....
!!HT:
      DO GDST=1, NRGRD
         DO GSRC=1, NRGRD
            IF ( HGSTGE(GDST,GSRC)%INIT ) THEN
               IF ( HGSTGE(GDST,GSRC)%NREC .NE. 0 ) THEN
                   DEALLOCATE (                                         &
                     HGSTGE(GDST,GSRC)%LJSEA , HGSTGE(GDST,GSRC)%NRAVG, &
                     HGSTGE(GDST,GSRC)%IMPSRC, HGSTGE(GDST,GSRC)%ITAG , &
                     HGSTGE(GDST,GSRC)%WGTH  , HGSTGE(GDST,GSRC)%SHGH , &
                     STAT=ISTAT )
                   CHECK_DEALLOC_STATUS ( ISTAT )
                 END IF
               IF ( HGSTGE(GDST,GSRC)%NSND .NE. 0 ) THEN
                   DEALLOCATE (                                         &
                     HGSTGE(GDST,GSRC)%ISEND ,                          &
                     STAT=ISTAT )
                   CHECK_DEALLOC_STATUS ( ISTAT )
                 END IF
               HGSTGE(GDST,GSRC)%NTOT = 0
               HGSTGE(GDST,GSRC)%NREC = 0
               HGSTGE(GDST,GSRC)%NRC1 = 0
               HGSTGE(GDST,GSRC)%NSND = 0
               HGSTGE(GDST,GSRC)%NSN1 = 0
               HGSTGE(GDST,GSRC)%NSMX = 0
               HGSTGE(GDST,GSRC)%INIT = .FALSE.
            END IF
         END DO
      END DO
      GDST=-999 ! unset grid
      GSRC=-999 ! unset grid
 
 
! -------------------------------------------------------------------- /
! 0.a Plan future behavior by setting logical variables.
 
 
      DO GDST=1, NRGRD
         IF ( GRIDS(GDST)%GTYPE .NE. RLGTYPE ) THEN
            ALL_REGULAR=.FALSE.
         END IF
      END DO
 
! Notes re: FLAGLL case: Old method calculates overlap area based on deg lat
!       and deg lon. New method (SCRIP) calculates overlap area based on real
!       distances. Therefore weights will not match for FLAGLL case, so we
!       do not perform checking for FLAGLL case.
 
      IF ( (.NOT.FLAGLL) .AND. ALL_REGULAR .AND. LSCRIP ) THEN
         IF ( IMPROC.EQ.NMPERR ) &
         WRITE (MDSE,'(/2A)')'We will check SCRIP calculations ', &
                    'against old method of calculating weights.'
         DO_CHECKING=.TRUE.
      END IF
 
      IF (DO_CHECKING .OR. (.NOT.LSCRIP)) OLD_METHOD=.TRUE.
 
! -------------------------------------------------------------------- /
! 0.b Check solution method
 
      IF ( (.NOT.LSCRIP) .AND. (.NOT.ALL_REGULAR) .AND. &
                                      (NRGRD.GT.1) ) THEN
         IF ( IMPROC.EQ.NMPERR ) &
              WRITE (MDSE,'(/3A)') ' *** ERROR WMGHGH: ', &
              'IRREGULAR or UNSTRUCTURED grid detected: this requires ', &
              'SCRIP switch.'
         CALL EXTCDE ( 999 )
      END IF
 
!
! -------------------------------------------------------------------- /
! 1.  Set boundary distance maps
! 1.a Check if needed
!
!!HT: FLGBDI is a flag set in WMMDATMD to .FALSE. and is used to identify
!!HT:        if the boundary distance maps have been initialized
!!HT:
!!HT: For each individual grid a map is generated identifying the distance
!!HT: to open boundaries (MAPBDI, saved in structure MDATA in WMMDATMD).
!!HT: This map is used later to choose if more that 1 high-res grids
!!HT: could provide data to a low-res grid. The high-res grid with data
!!HT: furthest away from its own open boundary will be used.
 
      IF ( .NOT. FLGBDI ) THEN
!
         IF ( FLAGLL ) THEN
            FACTOR = RADIUS * DERA
!notes: was FACTOR = RADIUS / 360. (bug fix)
         ELSE
            FACTOR = 1.
         END IF
!
! 1.b Loop over grids
!
         DO GDST=1, NRGRD
 
 
            CALL W3SETO ( GDST, MDSE, MDST )
            CALL W3SETG ( GDST, MDSE, MDST )
            CALL WMSETM ( GDST, MDSE, MDST )
 
!            IF ( GTYPE .EQ. UNGTYPE ) THEN
!               IF ( IMPROC.EQ.NMPERR ) &
!                    WRITE (MDSE,'(/2A)') ' *** ERROR WMGHGH: ', &
!                    'UNSTRUCTURED GRID SUPPORT NOT YET IMPLEMENTED ***'
!               CALL EXTCDE ( 999 )
!            END IF
 
!
! -------------------------------------------------------------------- /
! Inconsistent RANK vs NBI (warning message)
!   This was an error, now changed to a warning (see notes section)
            IF ( (GRANK(GDST).NE.1) .AND. (NBI.EQ.0) ) THEN
               IF ( IMPROC.EQ.NMPERR ) &
                  WRITE (MDSE,'(/2A)') ' WARNING in WMGHGH:   ', &
                  'NBI=0 AND RANK > 1 '
            END IF
 
! -------------------------------------------------------------------- /
! 1.c NBI=0, so computations not needed (test output only)
!
            IF ( (NBI.EQ.0) .OR. (GRANK(GDST).EQ.1) ) THEN
! (then do nothing except test output)
 
! -------------------------------------------------------------------- /
! 1.d NBI>0, Generate map with distances to boundary.
 
!!HT: Initialize MAPBDI
!!HT:   0. for active boundary points
!!HT:  -1. for points that are not considered at all (rescaled for test
!!HT:      output only, only negative value is essentially later).
!!HT:  -2. for points that still need to be processed.
 
            ELSE
 
               IF(IMPROC.EQ.NMPERR)WRITE(MDSE,'(A)') &
                   '  Generating map with distances to boundary.'
!  for purposes of screen output, would be useful to wait for other processors to catch up here...(if mpibdi switch used)
               ALLOCATE ( MDATAS(GDST)%MAPBDI(NY,NX), STAT=ISTAT )
               CHECK_ALLOC_STATUS ( ISTAT )
               MAPBDI => MDATAS(GDST)%MAPBDI
!
!KRL           Set up ranges for X. If not MPIBDI, just 1 to NX
               NX_BEG(IMPROC) = 1
               NX_END(IMPROC) = NX
!KRL           Setup complete
!
! -------------------------------------------------------------------- /
! Loop to determine MAPBDI
! -------------------------------------------------------------------- /
 
               IF(IMPROC.EQ.NMPERR)WRITE(MDSE,'(A)') &
                   'Starting MAPBDI 1st loop.'
 
               DO IDST=NX_BEG(IMPROC), NX_END(IMPROC)
                  IF(MOD(IDST,250).EQ.0)THEN
                     IF(LMPIBDI)THEN
                        WRITE(MDSE,'(4x,3(A,I5))')&
                             'processing column ',IDST,' out of ',NX, &
                             ' on processor ',IMPROC
                     ELSEIF(IMPROC.EQ.NMPERR)THEN
                        WRITE(MDSE,'(4x,2(A,I5))')&
                             'processing column ',IDST,' out of ',NX
                     ENDIF
                  ENDIF
                  DO JDST=1, NY
                     IF ( MAPSTA(JDST,IDST) .EQ. 0 ) THEN ! (excluded point)
                        MAPBDI(JDST,IDST) = -1. / SIG(1) * DTMAX ! new (bug fix)
                     ELSE IF ( ABS(MAPSTA(JDST,IDST)) .EQ. 2 ) THEN
                        ! (boundary point)
                        MAPBDI(JDST,IDST) =  0.
                     ELSE ! ABS(MAPSTA)=1 (sea point)
                        MAPBDI(JDST,IDST) = 1.0E+10
                     ENDIF ! IF MAPSTA
                  END DO ! DO JDST...
               END DO ! DO IDST...
 
! -------------------------------------------------------------------- /
 
               IF(IMPROC.EQ.NMPERR)WRITE(MDSE,'(A)') &
                   'Starting MAPBDI 2nd loop.'
 
               DO IBND=1,NX
                  IF ( (MOD(IBND,25).EQ.0) .AND. &
                       (IMPROC.EQ.NMPERR) ) THEN
                           WRITE(MDSE,'(4x,2(A,I5))') &
                          'bnd. point ',IBND,' out of ',NX
                  ENDIF
                  DO JBND=1,NY
                     IF ( ABS(MAPSTA(JBND,IBND)) .EQ. 2 ) THEN
                        ! (boundary point)
                        DO IDST=NX_BEG(IMPROC), NX_END(IMPROC)
                           DO JDST=1, NY
                              IF (ABS(MAPSTA(JDST,IDST)) .EQ. 1) THEN
                                 !....find distance to this boundary point.
                                 DD=FACTOR*W3DIST(FLAGLL,XGRD(JDST,IDST), &
                                    YGRD(JDST,IDST),XGRD(JBND,IBND), &
                                    YGRD(JBND,IBND))
 
! Notes: The origin of "0.58 * GRAV" is to translate from distance (in meters)
! to time (in seconds) required for a wave to travel from the boundary to point
! JDST,IDST based on a specific group velocity 0.58*grav would be the group
! velocity of a 7.3 s wave in deep water. Significance of T=7.3 s is explained
! in notes by HT below.
 
                                 DD=DD/ ( 0.58 * GRAV )
                                 MAPBDI(JDST,IDST)=MIN(MAPBDI(JDST,IDST),DD)
                              ENDIF
                           END DO ! DO JDST
                        END DO ! DO IDST
                     ENDIF ! (if BND point)
 
                  END DO ! DO JBND
               END DO ! DO IBND
 
               IF(IMPROC.EQ.NMPERR)WRITE(MDSE,'(A)') &
                   'Finished MAPBDI 2nd loop.'
 
! -------------------------------------------------------------------- /
 
 
               IF(IMPROC.EQ.NMPERR)WRITE(MDSE,'(A/)') &
               '  Finished generating map with distances to boundary.'
 
!...notes regarding old method of doing what we just did
!!HT:
!!HT: (1)
!!HT:
!!HT: CHANGE array is used to identify grid points that still need to
!!HT: be processed, and that are adjacent to points that have been
!!HT: processed. Only those points can be updated in this step of the
!!HT: loop started above here. The two loops below set the CHANGE array.
!!HT:
!!HT: (2)
!!HT:
!!HT: CHANGD identify if more points have been updated
!!HT:
!!HT: STX and STY are partial normalized distances, defined as the
!!HT: physical distance Delta Y ( FACTOR * SY ) and Delta X
!!HT: ( FACTOR * SX * XLAT(JDST) ), devided by the sistance traveled,
!!HT: which is CgMAX * DTMAX. CgMAX is approximately 1.15 * CgDEEP,
!!HT: or 1.15 * 0.5 * C_DEEP = 0.58 * GRAV / SIG(1). Since SIG(1) and
!!HT: DTMAX may vary, these two factors are not included in MAPBDI.
!!HT:
!!HT: This defines MAPBDI similar to an inverse CFL number.
!!HT:
!!HT: (3)
!!HT:
!!HT: ERROR : Should be CLAT(JDST), not CLATI(JDST) : "STX    = FACTOR * SX * CLATI(JDST) / ( 0.58 * GRAV )"
 
! 1.e Test output
!
!!HT: Note that SIG(1) and DTMAX are included here so that the map defines
!!HT: how many time steps DTMAX it takes to reach this place.
 
!
            END IF
         END DO
         FLGBDI = .TRUE.
      END IF
!
! -------------------------------------------------------------------- /
! 2.  Data sources for reconcilliation
! 2.a Loop over grids, processing check
!
 
!!HT: GRDHGH(GDST,0) was set in WMGLOW to identify how many grids may
!!HT: contribute from higher ranks to the present grid (GDST).
 
      ALLOCATE ( I1(NRGRD,NMPROC), I2(NRGRD,NMPROC),                  &
                 I3(NRGRD), I4(NRGRD), STAT=ISTAT )
      CHECK_ALLOC_STATUS ( ISTAT )
 
      LTAG0  = 0
 
              IF (LSCRIPNC) CALL MPI_BARRIER(MPI_COMM_WORLD, IERR_MPI)
 
      LOWRANK_GRID : DO GDST=1, NRGRD
 
 
! Test output
 
!
         IF ( GRDHGH(GDST,0) .EQ. 0 ) THEN ! no grids of higher rank than this
                                           ! one.
         ELSE ! processing required
 
!
! 2.b Process grid
! 2.b.1 Preparations
!
!!HT: Grid I has higher rank grids covering it, we now set up MAPTST
!!HT:   MAPTST shows from which gid the data is averages.
!!HT:   INFLND inferred land points based on land in high-res grids.
!!HT:
            CALL W3SETO ( GDST, MDSE, MDST )
            CALL W3SETG ( GDST, MDSE, MDST )
            CALL WMSETM ( GDST, MDSE, MDST )
 
!  W3SETG set ICLOSE for us, and we have determined that there is
!   interaction between high and low rank. So this is a good point
!   to check the closure type.
 
            IF ( ICLOSE .EQ. ICLOSE_TRPL ) THEN
              IF ( IMPROC.EQ.NMPERR ) WRITE(MDSE,*)'SUBROUTINE WMGHGH IS'// &
              ' NOT YET ADAPTED FOR TRIPOLE GRIDS. STOPPING NOW.'
              CALL EXTCDE ( 1 )
            END IF
 
            ALLOCATE ( MAPTST(NY,NX), INFLND(NY,NX), STAT=ISTAT )
            CHECK_ALLOC_STATUS ( ISTAT )
            MAPTST = 0
            INFLND = 0
 
!################################################################
! Start new block of code: Calculate weights by calling SCRIP interface
!################################################################
 
! Notes on grid variables:
! GRIDS(GSRC)%{grid variable} (src grid, high rank, high resolution grid)
! GRIDS(GDST)%{grid variable} (dst grid, low rank, low resolution grid)
 
!   At this point, we are working on a particular low rank (dst) grid.
!   We will save our weight information in the structure "ALLWGTS".
!   For this dst grid, it is possible to have many src grids. That is
!      why we store it this way.
!   First, we ALLOCATE ALLWGTS from 1 up to the largest value of all
!      the possible source grids. This will be referenced as "GSRC"
!      Not every value of GSRC will be filled (e.g. "1" usually isn't filled)
!      but since we are doing this as a derived data type, we are still
!      efficient in terms of memory usage.
!   Inside SCRIP interface, we have:
!      type weight_data
!         integer (kind=int_kind)              :: n    ! number of weights for
!                                               dst cell, formerly npnts(:)
!         real    (kind=dbl_kind), allocatable :: w(:) ! weights, sized by n,
!                                               formerly wxwy(:,:)
!         integer (kind=int_kind), allocatable :: k(:) ! source grid cells,
!                                               sized by n, formerly KSRC(:,:)
!      end type weight_data
!   ....
!     type(weight_data), allocatable :: WGTDATA(:)
!   ....
!     ALLOCATE ( WGTDATA(grid2_size), STAT=ISTAT ) ! grid2=destination grid
!     CHECK_ALLOC_STATUS ( ISTAT )
 
 
! Next, we loop through the src grids for the dst grid that we are working on.
 
 
! Next, we call SCRIP for this src grid.
! Conditions for calling SCRIP are:
! 1) Not using L_STOP: in this case, all processes need all the weight
!    information, so all processes need to call SCRIP for all grid pairs
! OR
! 2) Using L_STOP, writing .nc files and not reading .nc files. With
!    L_STOP, different processors are doing different things, and so
!    have different settings for L_WRITE. L_READ is the same for all
!    processors, since it is simply based on whether the file already
!    exists.
 
 
 
 
 
 
! SCRIP has now created the data strucure "WGTDATA" and stored the weights
! in it. However, this is only for the present src grid. We want to store the
! data for all the src grids. Thus, we use a new data structure of type
! "ALLWGT" to store this data. First though, we need to ALLOCATE it:
! (note: "k" is equivalent to isea, but includes *all* points)
 
 
! Now that we have it allocated, we can just copy WGTDATA into ALLWGTS
 
! Notes re: short and long way to do this:
! pgf90 on IBM Opteron, gfortran, g95, xlf, all tested ok with "short method"
! pgf90 on our linux workstations (Intel) requires the "long method"
! (possible compiler bug)
! ALLWGTS(GSRC)%WGTDATA = WGTDATA                               !short method
 
! BEGIN long method for filling derived data type "ALLWGTS"
 
 
! END long method for filling derived data type "ALLWGTS"
 
! We're done with WGTDATA, so we can DEALLOCATE it. This is important,
! since it will be allocated again the next time SCRIP is called.
 
 
! Here's a "test output" block of code to demonstrate how the weights can
! be called up from ALLWGTS...and to verify that the data is stored properly.
! (again note that "k" is equivalent to isea, but includes *all* points)
 
 
 
!################################################################
! End new block of code: Calculate weights by calling SCRIP interface
!################################################################
 
! 2.b.2 Find points used for boundary data in higher ranked grids
!
!!HT: These points are marked in MAPTST as negative values to assure
!!HT: that the grid poits used for boundary data are not getting
!!HT: averaged values from high-reswolution grids as that will result
!!HT: in cyclic, possibly non-conservative updating.
!!HT:
!!HT: NBI2S has all necessary data set in WMGLOW as called before WMGHGH.
!!HT:
!!HT: JJ loop goes over grids that reviously have been identified as
!!HT: getting data from the grid presently cousidered.
!
! notes: The purpose of this is to identify points
!        that should not be used in the averaging procedure. It is
!        related to statement in Tolman (OM, 2008):  "Second, Eq (7) is not
!        applied to grid points in the low resolution grid that contribute
!        to boundary data for the high resolution grid. This avoids cyclic
!        updating of data between grids.
 
! notes:   GRDHGH(GDST,0) is number of grids of higher rank than the present
!                         grid (GDST)
!          GRDHGH(GDST,1...etc.) is the grid number
 
! notes: Setting MAPTST=negative here is probably overkill, since it means
!        we will not have weights for this point. However, to change this,
!        we would need a new variable to use in its place, since we need
!        to mark the point for use in STMASK determination.
 
            DO JJ=1, GRDHGH(GDST,0)
               GSRC   = GRDHGH(GDST,JJ)
               DO IB=1, SIZE(MDATAS(GSRC)%NBI2S(:,1))
                  IF ( MDATAS(GSRC)%NBI2S(IB,1) .EQ. GDST ) THEN
                     IDST     = MAPSF(MDATAS(GSRC)%NBI2S(IB,2),1)
                     JDST     = MAPSF(MDATAS(GSRC)%NBI2S(IB,2),2)
                     MAPTST(JDST,IDST) = - GSRC
                  END IF
               END DO
            END DO
            GSRC   = -999 ! unset grid
!
! 2.b.3 Range of coverage per grid
 
!!HT:
!!HT: In this JJ loop, we go over all higher resolution grids to find
!!HT: ranges that can be averaged to replace data in the 'I' (GDST) grid.!
!!HT:
 
            ALLOCATE ( IDSTL(GRDHGH(GDST,0)), IDSTH(GRDHGH(GDST,0)),  &
                       JDSTL(GRDHGH(GDST,0)), JDSTH(GRDHGH(GDST,0)),  &
                       GRIDOK(GRDHGH(GDST,0)),BDIST(GRDHGH(GDST,0)),  &
                       STAT=ISTAT                                     )
            CHECK_ALLOC_STATUS ( ISTAT )
 
            IF (OLD_METHOD) THEN
                ALLOCATE (BDIST_OM(GRDHGH(GDST,0)), STAT=ISTAT )
                CHECK_ALLOC_STATUS ( ISTAT )
              END IF
 
!
!       Notes: For case of lower ranked grid GDST being irregular, grid indices
!              i and j do not correspond to x and y, so optimization
!              by limiting search in manner of pre-curvilinear versions of
!              WW3 is not appropriate.
!
            IF ( (GTYPE .EQ. CLGTYPE).or.(GTYPE .EQ. UNGTYPE) ) THEN
 
               IDSTLA   = 1
               IDSTHA   = NX
               JDSTLA   = 1
               JDSTHA   = NY
 
            ELSE
 
! loop through higher ranked grids GSRC
 
               DO JJ=1, GRDHGH(GDST,0)
                  GSRC   = GRDHGH(GDST,JJ)
!!HT:
!!HT: XL,XH, and YL,YH and the low and high (X,Y) values of the grid box
!!HT: in the grid 'I' fro which the high-res data needs to be averaged.
!!HT: To be efficient, we compute a range of high-res grid point that
!!HT: could be considered, rather than looking through the whole grid.
!!HT: This will work only for the old grids, not for the newer curvilinear
!!HT: and unstructured grids.
!!HT:
!!HT: This sets the range in the low-res grid to consider.
!
!       Notes (HLT): outer edges already taken off here ...
!              will not work in a simple way for spherical grids,
!              so we don't even try ....
!
!       Notes: SX and SY are only used in cases where GTYPE .NE. CLGTYPE,
!              i.e. regular grids. In case of regular grids, SX and SY
!              can be replaced with HPFAC HQFAC, if desired.
!
! find upper and lower bounds of higher ranks grids
 
                  IF ( (GRIDS(GSRC)%GTYPE .EQ. CLGTYPE) .OR. &
                     (GRIDS(GSRC)%GTYPE .EQ. UNGTYPE) ) THEN
 
!       Notes: in case of irregular grids, there is no obvious way to
!              offset by dx/2 dy/2, so we omit that sliver (thus we increase
!              search area slightly).
 
                     XL=MINVAL(GRIDS(GSRC)%XGRD)
                     YL=MINVAL(GRIDS(GSRC)%YGRD)
                     XH=MAXVAL(GRIDS(GSRC)%XGRD)
                     YH=MAXVAL(GRIDS(GSRC)%YGRD)
 
                  ELSE
 
                     XL     = GRIDS(GSRC)%X0 + 0.5 * GRIDS(GSRC)%SX
                     XH     = GRIDS(GSRC)%X0 + ( REAL(GRIDS(GSRC)%NX) - 1.5 ) &
                            * GRIDS(GSRC)%SX
                     YL     = GRIDS(GSRC)%Y0 + 0.5 * GRIDS(GSRC)%SY
                     YH     = GRIDS(GSRC)%Y0 + ( REAL(GRIDS(GSRC)%NY) - 1.5 ) &
                            * GRIDS(GSRC)%SY
 
                  ENDIF ! IF ( GRIDS(GSRC)%GTYPE .EQ. CLGTYPE )
 
!
! find where this falls in the current (low) ranked grid
 
                  IF ( FLAGLL ) THEN
                     IDSTL(JJ) = 1
                     IDSTH(JJ) = NX
                  ELSE
 
! Notes: from check "IF ( GTYPE .EQ. CLGTYPE ) THEN" above, we know that
!         GTYPE .NE CLGTYPE....so it is safe to use SX SY etc here
 
                     IDSTL(JJ) = 2 + INT( (XL-X0)/SX + 0.49 )
                     IDSTH(JJ) = 1 + INT( (XH-X0)/SX - 0.49 )
                  END IF
 
                  JDSTL(JJ) = 2 + INT( (YL-Y0)/SY + 0.49 )
                  JDSTH(JJ) = 1 + INT( (YH-Y0)/SY - 0.49 )
 
                  IDSTL(JJ) = MAX (  1 , IDSTL(JJ) )
                  IDSTH(JJ) = MIN ( NX , IDSTH(JJ) )
                  JDSTL(JJ) = MAX (  1 , JDSTL(JJ) )
                  JDSTH(JJ) = MIN ( NY , JDSTH(JJ) )
!
               END DO  ! end loop through higher ranked grids
               GSRC   = -999 ! unset grid
!
 
! save the extremities of that set of high-ranked grids
               IDSTLA   = MINVAL(IDSTL)
               IDSTHA   = MAXVAL(IDSTH)
               JDSTLA   = MINVAL(JDSTL)
               JDSTHA   = MAXVAL(JDSTH)
 
            ENDIF ! IF ( (GTYPE .EQ. CLGTYPE ) .or. (GTYPE .EQ. UNGTYPE))
 
! loop through higher ranked grids
 
!
! 2.b.4 Point by point check
!
! Notes: We loop through all grids of higher rank
!        GSRC=the grid number of the higher rank grid.
!        NLMAX is used for dimensioning purposes.
!        It is apparently using the ratio between the resolution
!        of the low rank grid (GDST) and high rank grid (GSRC)
!        Obviously, we cannot use this calculation for irregular grids.
 
            NLMAX  = 0
            DO JJ=1, GRDHGH(GDST,0)
               GSRC   = GRDHGH(GDST,JJ)
 
! Notes: NLMAX is used to dimension TMPINT,TMPRL, and to set ITAG and LTAG
!        (MPI case).
!        As we remove more of the older code, it may turn out that
!        NLMAX is no longer needed, in which case we can remove this
!        block of code. For example, the weights data structure is introduced
!        to WW3 already dimensioned.
 
               IF ( GRIDS(GDST)%GTYPE .EQ. CLGTYPE ) THEN
                  DX_MAX_GDST=MAXVAL(GRIDS(GDST)%HPFAC)
                  DY_MAX_GDST=MAXVAL(GRIDS(GDST)%HQFAC)
               ELSEIF ( GRIDS(GDST)%GTYPE .EQ. RLGTYPE ) THEN
                  DX_MAX_GDST=GRIDS(GDST)%SX
                  DY_MAX_GDST=GRIDS(GDST)%SY
               ELSE
                 IsFirst=1
                 DIST_MAX=0
                 DIST_MIN=0
                 DO ITRI=1,GRIDS(GDST)%NTRI
                   DO IT=1,3
                     IF (IT.eq.3) THEN
                       JT=1
                     ELSE
                       JT=IT+1
                     END IF
                     IM1=GRIDS(GDST)%TRIGP(ITRI,IT)
                     IM2=GRIDS(GDST)%TRIGP(ITRI,JT)
                     eDist=W3DIST(FLAGLL, GRIDS(GDST)%XYB(IM1,1), &
                           GRIDS(GDST)%XYB(IM1,2), &
                           GRIDS(GDST)%XYB(IM2,1), GRIDS(GDST)%XYB(IM2,2))
                     IF (IsFirst.eq.1) THEN
                       DIST_MAX=eDist
                       DIST_MIN=eDist
                       IsFirst=0
                     ELSE
                       IF (eDist.gt.DIST_MAX) THEN
                         DIST_MAX=eDist
                       END IF
                       IF (eDist.lt.DIST_MIN) THEN
                         DIST_MIN=eDist
                       END IF
                     END IF
                   END DO
                 END DO
                 DX_MAX_GDST=DIST_MAX
                 DY_MAX_GDST=DIST_MAX
               END IF
 
               IF ( GRIDS(GSRC)%GTYPE .EQ. CLGTYPE ) THEN
                  DX_MIN_GSRC=MINVAL(GRIDS(GSRC)%HPFAC)
                  DY_MIN_GSRC=MINVAL(GRIDS(GSRC)%HQFAC)
               ELSEIF ( GRIDS(GSRC)%GTYPE .EQ. RLGTYPE ) THEN
                  DX_MIN_GSRC=GRIDS(GSRC)%SX
                  DY_MIN_GSRC=GRIDS(GSRC)%SY
               ELSE
                 IsFirst=1
                 DIST_MAX=0
                 DIST_MIN=0
                 DO ITRI=1,GRIDS(GSRC)%NTRI
                   DO IT=1,3
                     IF (IT.eq.3) THEN
                       JT=1
                     ELSE
                       JT=IT+1
                     END IF
                     IM1=GRIDS(GSRC)%TRIGP(ITRI,IT)
                     IM2=GRIDS(GSRC)%TRIGP(ITRI,JT)
                     eDist=W3DIST(FLAGLL, GRIDS(GSRC)%XYB(IM1,1), &
                           GRIDS(GSRC)%XYB(IM1,2), &
                           GRIDS(GSRC)%XYB(IM2,1), GRIDS(GSRC)%XYB(IM2,2))
                     IF (IsFirst.eq.1) THEN
                       DIST_MAX=eDist
                       DIST_MIN=eDist
                       IsFirst=0
                     ELSE
                       IF (eDist.gt.DIST_MAX) THEN
                         DIST_MAX=eDist
                       END IF
                       IF (eDist.lt.DIST_MIN) THEN
                         DIST_MIN=eDist
                       END IF
                     END IF
                   END DO
                 END DO
                 DX_MIN_GSRC=DIST_MIN
                 DY_MIN_GSRC=DIST_MIN
               END IF
 
! notes: original code was much simpler:
!              NLMAX  = MAX ( NLMAX , (2+INT(SX/GRIDS(J)%SX+0.001)) *  &
!                                     (2+INT(SY/GRIDS(J)%SY+0.001)) )
 
               NLMAX = MAX ( NLMAX ,                             &
                       (2+INT(DX_MAX_GDST/DX_MIN_GSRC+0.001)) *  &
                       (2+INT(DY_MAX_GDST/DY_MIN_GSRC+0.001)) )
 
 
 
            END DO !  DO JJ=1, GRDHGH(GDST,0)
            GSRC=-999 ! unset grid
 
! Notes regarding 3 possible scenarios:
!        If only using SCRIP, then
!            * set NLMAX=NLMAX_SCRIP here.
!            * TMPRL_OM will not be created
!            * TMPRL will be calculated using SCRIP
!        If only using old method
!            * NLMAX is already set, and SCRIP switch does not exist, so
!               nothing is done here
!            * both TMPRL and TMPRL_OM will be dimensioned
!            * TMPRL_OM will be calculated
!            * TMPRL_OM will be copied to TMPRL for use
!        If using both methods ("DO_CHECKING")
!            * set NLMAX=MAX(NLMAX, NLMAX_SCRIP) here.
!            * both TMPRL_OM and TMPRL will be created
!            * both will be calculated using the 2 methods, and
!               checked against each other
!            * the SCRIP version of weights (TMPRL) will be the ones used.
 
 
            IF(NLMAX.GT.100)THEN
               WRITE(MDSE,'(/A,I8)') &
               'WARNING: unusually large value for NLMAX : ',NLMAX
            END IF
 
            NRTOT  = 0
            IF(OLD_METHOD)THEN
               ALLOCATE ( TMPINT_OM(NX*NY,-4:NLMAX), STAT=ISTAT )
               CHECK_ALLOC_STATUS ( ISTAT )
               ALLOCATE ( TMPRL_OM(NX*NY,0:NLMAX),   STAT=ISTAT )
               CHECK_ALLOC_STATUS ( ISTAT )
            ENDIF
            ALLOCATE ( TMPINT(NX*NY,-4:NLMAX),       STAT=ISTAT )
            CHECK_ALLOC_STATUS ( ISTAT )
            ALLOCATE ( TMPRL(NX*NY,0:NLMAX),         STAT=ISTAT )
            CHECK_ALLOC_STATUS ( ISTAT )
            ALLOCATE ( TMPLOG(NX*NY),                STAT=ISTAT )
            CHECK_ALLOC_STATUS ( ISTAT )
!
            ALLOCATE ( LTAG(NLMAX), STAT=ISTAT )
            CHECK_ALLOC_STATUS ( ISTAT )
            DO JJ=1, NLMAX
              LTAG(JJ) = JJ + LTAG0
              END DO
!
!!HT:
!!HT: After the search range is set, we are actually searching in the
!!HT: high-res grid. IDST, JDST are counters in the grid to which the
!!HT: averaged data is to go. XA and YA are center locatons of target
!!HT: grid. Necxt two loops over all relevant point in target grid.
!!HT:
 
! Notes: This is the start of the large loop through the individual
!        grid points of the low-rank grid.
!        The checks below for JDST.LT.JDSTLA , IDST.LT.IDSTLA etc are to save
!        time but will only be useful for the case of regular grids.
 
            LOWRANK_J : DO JDST=1, NY
               IF ( JDST.LT.JDSTLA .OR. JDST.GT.JDSTHA ) CYCLE
 
               LOWRANK_I : DO IDST=1, NX
                  IF ( IDST.LT.IDSTLA .OR. IDST.GT.IDSTHA ) CYCLE
 ! check that this is a sea point
                  IF ( ABS(MAPSTA(JDST,IDST)) .NE. 1 ) CYCLE
 ! MAPTST: see Section 2.b.2 above
                  IF ( MAPTST(JDST,IDST) .LT. 0 ) CYCLE
                  XA     = XGRD(JDST,IDST) ! old code: X0 + REAL(IDST-1)*SX
                  YA     = YGRD(JDST,IDST) ! old code: Y0 + REAL(JDST-1)*SY
 
!!HT: For each point in the target grid, loop over all relevant high-res
!!HT: grid (JJ loop).
 
                  NROK   = 0
 
! notes:   GRDHGH(GDST,0) is number of grids of higher rank than the present
!          grid (GDST)
!          GRDHGH(GDST,1...etc.) is the grid number
 
                  DO JJ=1, GRDHGH(GDST,0)
                     GSRC   = GRDHGH(GDST,JJ)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Start counting using old method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!       Note for LLG: Assumption is made that the higher ranked grid
!                     cannot be global.
!
!!HT: Set search range in [candidate] high-res grid.
 
! Notes: The quantities XL YL XH YH apear to be a bounding box in
!        index space for later searching within the high rank grid (or
!        otherwise making computations from the high rank grid). They
!        are the distance from the coarse grid point to the origin of
!        the high rank grid, measured in grid cells of the high rank
!        grid. So, they are the i and j values in the high rank
!        bounding the low rank grid cell.
 
                     IF (OLD_METHOD)THEN
 ! ...then  we do the counting using the old method
 
! Notes: Resulting "old method" variables are saved with "_OM" suffix.
 
                        IF ( FLAGLL ) THEN
                           DXC  = MOD ( 1080.+XA-GRIDS(GSRC)%X0 , 360. )
                           XL   = 1. + (DXC-0.5*SX)/GRIDS(GSRC)%SX
                           XH   = 1. + (DXC+0.5*SX)/GRIDS(GSRC)%SX
                        ELSE
                           XL   = 1. + (XA-GRIDS(GSRC)%X0-0.5*SX)/GRIDS(GSRC)%SX
                           XH   = 1. + (XA-GRIDS(GSRC)%X0+0.5*SX)/GRIDS(GSRC)%SX
                        END IF
                        YL     = 1. + (YA-GRIDS(GSRC)%Y0-0.5*SY)/GRIDS(GSRC)%SY
                        YH     = 1. + (YA-GRIDS(GSRC)%Y0+0.5*SY)/GRIDS(GSRC)%SY
 
                        ISRCL  = NINT(XL+0.01)
                        ISRCH  = NINT(XH-0.01)
                        JSRCL  = NINT(YL+0.01)
                        JSRCH  = NINT(YH-0.01)
 
                        IF ( ISRCL.LT.1 .OR. ISRCH.GT.GRIDS(GSRC)%NX .OR.          &
                             JSRCL.LT.1 .OR. JSRCH.GT.GRIDS(GSRC)%NY ) THEN
!                            dst point not in src grid, so go to next src grid
                           GRIDOK(JJ) = .FALSE.  ! does this get used anywhere?
                           CYCLE ! leave GSRC loop
                        END IF
 
!!HT: Loop over search range in high-res grid, ISRC and JSRC loops.
!!HT:  NR0_OM counts high-res grid points with MAPSTA=0, etc.
!!HT:  NRL_OM separately identifies explitcit land points.
!!HT:  BDIST_OM saves the boundary data from the source grid.
!!HT:
 
! Notes: We appear to be searching for the smallest boundary distance and
!        doing some counting
!        Purpose of counting is unknown (for dimensioning?)
 
! Initialize
                        NR0_OM    = 0 ! counter of MAPSTA=0 (indicates
                                      ! excluded point)
                        NRL_OM    = 0 ! counter of MAPSTA=0 (indicates
                                      ! excluded point) and MAPST2=0
                                      ! (indicates land)
                        NR1_OM    = 0 ! counter of |MAPSTA|=1
                                      ! (indicates sea point)
                        NR2_OM    = 0 ! counter of |MAPSTA|=2
                                      ! (indicates boundary point)
                        BDIST_OM(JJ) = 9.99E33
 
                        DO ISRC=ISRCL, ISRCH
                           DO JSRC=JSRCL, JSRCH
                              IF (GRIDS(GSRC)%MAPSTA(JSRC,ISRC).EQ.0) THEN
                                 ! excluded point
                                 NR0_OM    = NR0_OM + 1
 
! Notes: Q: What does MAPST2=0 indicate?
!        A: MAPST2 is the "second grid status map"
!           For disabled points (MAPSTA=0) , MAPST2 indicates land (0) or
!           excluded (1). For sea and active boundary points, MAPST2 indicates
!           a) ice coverage b) drying out of points c) land in moving grid or
!           inferred land in nesting and d) masked in two-way nesting
 
                                 IF (GRIDS(GSRC)%MAPST2(JSRC,ISRC).EQ.0) &
                                   NRL_OM = NRL_OM + 1
                              ELSE IF (ABS(GRIDS(GSRC)%MAPSTA(JSRC,ISRC)) &
                                      .EQ.1) THEN ! sea point
                                 NR1_OM    = NR1_OM + 1
 
! Notes: check against stored "distance to boundary point"
!        This BDIST_OM array will be used later, when we select
!        the high rank grid to average from.
 
                                 BDIST_OM(JJ) = MIN ( BDIST_OM(JJ) ,     &
                                      MDATAS(GSRC)%MAPBDI(JSRC,ISRC) )
                              ELSE IF (ABS(GRIDS(GSRC)%MAPSTA(JSRC,ISRC)) &
                                      .EQ.2) THEN ! bnd point
                                 NR2_OM    = NR2_OM + 1
                              END IF
                           END DO ! DO JSRC=JSRCL, JSRCH
                        END DO ! DO ISRC=ISRCL, ISRCH
 
                     END IF ! (if OLD_METHOD)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Done with counting using old method.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Start counting using new method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
! Initialize
 
 
! Notes on variables used here:
! IDST, JDST given by loop, NIDST set above, the rest we need to set here
 
 
                               ! sea point
 
! Pull NR0, etc. from storage...
!                             counter of MAPSTA=0 (indicates excluded point)
!                             counter of MAPSTA=0 (indicates excluded point)
!                             and MAPST2=0 (indicates land)
!                             counter of |MAPSTA|=1 (indicates sea point)
!                             counter of |MAPSTA|=2 (indicates boundary point)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finished counting using new method.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
! Compare results
                     IF(DO_CHECKING)THEN
! then it is OK to compare with the values that we got using the old method
                        IF(NR0_OM.NE.NR0)THEN
                           IF ( IMPROC.EQ.NMPERR )WRITE (MDSE,'(/1A,2(I8))') &
                                ' *** ERROR WMGHGH: NR0_OM,NR0 = ',NR0_OM,NR0
                           CALL EXTCDE ( 999 )
                        ENDIF
                        IF(NR1_OM.NE.NR1)THEN
                           IF ( IMPROC.EQ.NMPERR )WRITE (MDSE,'(/1A,2(I8))') &
                                ' *** ERROR WMGHGH: NR1_OM,NR1 = ',NR1_OM,NR1
                           CALL EXTCDE ( 999 )
                        ENDIF
                        IF(NR2_OM.NE.NR2)THEN
                           IF ( IMPROC.EQ.NMPERR )WRITE (MDSE,'(/1A,2(I8))') &
                                ' *** ERROR WMGHGH: NR2_OM,NR2 = ',NR2_OM,NR2
                           CALL EXTCDE ( 999 )
                        ENDIF
                        IF(NRL_OM.NE.NRL)THEN
                           IF ( IMPROC.EQ.NMPERR )WRITE (MDSE,'(/1A,2(I8))') &
                                ' *** ERROR WMGHGH: NRL_OM,NRL = ',NRL_OM,NRL
                           CALL EXTCDE ( 999 )
                        ENDIF
                        IF(BDIST_OM(JJ).NE.BDIST(JJ))THEN
                           IF ( IMPROC.EQ.NMPERR ) &
                                WRITE (MDSE,'(/2A,2(F12.5))') &
                                ' *** ERROR WMGHGH: ', &
                                ' BDIST_OM(JJ),BDIST(JJ) = ', &
                                BDIST_OM(JJ),BDIST(JJ)
                           CALL EXTCDE ( 999 )
                        ENDIF
                     END IF ! (if DO_CHECKING)
 
! Notes: We are done with the counting. If we didn't use SCRIP to get NR0,
!  etc., then we need to set them using the _OM variables.
 
                     IF(.NOT.LSCRIP)THEN
                        NR0=NR0_OM
                        NR1=NR1_OM
                        NR2=NR2_OM
                        NRL=NRL_OM
                        BDIST=BDIST_OM
                     END IF
 
! Notes: Potential future improvement: for irregular grids, it would make
! more sense to use the overlapped area, rather than simply counting cells
! to decide on "inferred land". However, since grid cell size it typically
! fairly uniform locally, the current approach will suffice for now.
 
! Notes: This is the only place that the "NRL" "NR0" "NR1" and "NR2" variables
!        are used directly. They affect MAPST2 indirectly below.
!        The calculation itself is essentially "50% or more of the grid
!        cells are land".
 
                     IF ( NRL .GT. (NR0+NR1+NR2)/2 ) THEN
 
! Notes: This is not considered an OK grid (NROK is not incremented) and
!        it is considered "inferred land"
 
                        INFLND(JDST,IDST) = 1
                     ELSE
                        GRIDOK(JJ) = NR1.GT.0 .AND. NR2.EQ.0
 
! Notes: for a grid cell to be considered "OK", we require that there is
! at least one sea point being used, and no boundary points being used
 
                        IF ( GRIDOK(JJ) ) NROK = NROK + 1
                     END IF
 
                  END DO !  GSRC loop
                  GSRC=-999 ! unset grid
 
                  IF ( NROK .EQ. 0 ) THEN
 
! Notes:  exit IDST loop since there are no "OK" source grid cells for this
! dst point. At this point, INFLND could be 1, but isn't necessarily 1
 
                     CYCLE
 
                  ELSE
 
! Notes: If any grids are OK for this dst point, then we override any prior
! setting of INFLD=1. Apparently this is for the situation of having some src
! grids giving INFLD=1 and another giving INFLD=0 for the same dst point.
! I wouldn't expect this to happen very often.
 
                     INFLND(JDST,IDST) = 0
 
                  END IF
!
! 2.b.5 Select source grid
!
! Notes: It appears that we are selecting the high rank grid from
!        which we will perform the averaging.
!        The code is written such that the first higher rank
!        grid that we find has the rank that we want, but isn't necessarily the
!        grid that we want.
!        Are grids necessarily in order of rank? If so, then we want the grid
!        that is higher rank but of nearest rank to the present grid.
!        Anyway, once we have decided on the grid rank that we want, we select
!        the specific grid according to criterion: larger distance to
!        boundary = better
!        Keep in mind that this grid is selected for *this* (IDST,JDST) and not
!        necesssarily for the next...
 
                  JF     = 0
 
!!HT: Another loop over all high-res grid to decide which grid will
!!HT: be used to average data. If more than 1 grids are available,
!!HT: the boundary distance in the high-res grid, stored in BDIST is
!!HT: used to make the choice.
 
!!ER: Old logic was to select a grid that is of the "next rank up,
!!ER: for which data is available". This was done by searching
!!ER: from 1 to GRDHGH (remember that the available source grids
!!ER: are in order of rank), and exiting when the rank increased.
!!ER: The problem with selecting the "lowest rank grid with rank
!!ER: greater than that of GDST" is that at this point, we have
!!ER: no knowledge of what will be masked in that SRC grid, since
!!ER: we haven't updated MAPSTA for that GSRC yet, based on what
!!ER: points are covered by higher rank grids (in case of masking
!!ER: computations). We can avoid this problem by reversing the
!!ER: order of GSRC search (from highest rank to lowest rank of
!!ER: higher rank). For example, grid 1 wants data from grid 2,
!!ER: but just gets zeros because grid 2 is masked there, because
!!ER: grid 2 is masked by grid 3 in Section 2.3.2 below. Going
!!ER: directly to GSRC=3 for GDST=1 (skipping GSRC=2) avoids
!!ER: this: the highest rank at that location will never be
!!ER: masked by a higher rank grid.
 
                  DO JJ=GRDHGH(GDST,0),1,-1
!old              DO JJ=1, GRDHGH(GDST,0) ! used before Aug 2014
 
                     GSRC   = GRDHGH(GDST,JJ)
 
                     IF ( GRIDOK(JJ) ) THEN
                        IF ( JF .EQ. 0 ) THEN ! we haven't already found a grid
                           JF     = GSRC ! now we have found a grid.
                           JR     = GRANK(GSRC)
! this is the rank that we want....the rank of the first grid that we find
                           JD     = BDIST(JJ) ! larger distance = better
                        ELSE
! we already found a grid, but maybe this one is better
                           IF ( GRANK(GSRC) .NE. JR ) EXIT
! this is not the rank that we want
                           IF ( BDIST(JJ) .GT. JD ) THEN
! we like this grid better
                              JF     = GSRC
                              JD     = BDIST(JJ)
                           END IF
                        END IF
                     END IF
                  END DO
                  GSRC=JF
 
!!HT: Data for grid point IDST,JDST in the low-res grid will be taken from
!!HT: high-res grid GSRC.
 
                  MAPTST(JDST,IDST) = GSRC
!
! 2.b.6 Store data (temp)
!
! Notes: This section is for calculations of weights for the
!        area-weighted averaging.
 
                  NRTOT  = NRTOT + 1
                  TMPINT(NRTOT,-4) = IDST
                  TMPINT(NRTOT,-3) = JDST
                  TMPINT(NRTOT,-2) = MAPFS(JDST,IDST)
                  TMPINT(NRTOT,-1) = GSRC
                  TMPRL (NRTOT, 0) = JD * SIG(1) / DTMAX
 
! Notes: Calculation for XL YL XH YH is same as it was in section 2.b.4, so
!        see notes in that section.
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!...Begin block of code for computing weights using old method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
                  IF (OLD_METHOD)THEN
! it is OK to do the counting using the old method
! (These variables are saved with "_OM" suffix)
 
                     DO ITMP=-4,-1
                        TMPINT_OM(NRTOT,ITMP)=TMPINT(NRTOT,ITMP)
                     END DO
                     TMPRL_OM(NRTOT,0)=TMPRL(NRTOT,0)
 
                     IF ( FLAGLL ) THEN
                        DXC    = MOD ( 1080.+XA-GRIDS(GSRC)%X0 , 360. )
                        XL     = 1. + (DXC-0.5*SX)/GRIDS(GSRC)%SX
                        XH     = 1. + (DXC+0.5*SX)/GRIDS(GSRC)%SX
                     ELSE
                        XL     = 1. + (XA-GRIDS(GSRC)%X0-0.5*SX)/GRIDS(GSRC)%SX
                        XH     = 1. + (XA-GRIDS(GSRC)%X0+0.5*SX)/GRIDS(GSRC)%SX
                     END IF
                     YL     = 1. + (YA-GRIDS(GSRC)%Y0-0.5*SY)/GRIDS(GSRC)%SY
                     YH     = 1. + (YA-GRIDS(GSRC)%Y0+0.5*SY)/GRIDS(GSRC)%SY
 
! Notes: Here, we save the search bounds. These bounds are given in terms of
!        index space of the high rank grid. "L" and "H" here do *not* refer
!        to grid rank! They refer to lower and upper bounds.
 
                     ISRCL    = NINT(XL+0.01)
                     ISRCH    = NINT(XH-0.01)
                     JSRCL    = NINT(YL+0.01)
                     JSRCH    = NINT(YH-0.01)
 
                     WTOT   = 0.
                     NLOC_OM   = 0
                     DO ISRC=ISRCL, ISRCH
                        WX     = MIN(XH,REAL(ISRC)+0.5) - MAX(XL,REAL(ISRC)-0.5)
                        DO JSRC=JSRCL, JSRCH
                           IF (ABS(GRIDS(GSRC)%MAPSTA(JSRC,ISRC)).EQ.1) THEN
! sea point
                              WY     = MIN(YH,REAL(JSRC)+0.5) -               &
                                   MAX(YL,REAL(JSRC)-0.5)
                              WTOT   = WTOT + WX*WY
                              NLOC_OM   = NLOC_OM + 1
! Notes: check here that we are sufficiently dimensioned.
                              IF ( NLOC_OM .GT. NLMAX ) THEN
                                 IF ( IMPROC.EQ.NMPERR ) WRITE (MDSE,1020)
                                 CALL EXTCDE(1020)
                              END IF
                              TMPINT_OM(NRTOT,NLOC_OM) =  &
                                 GRIDS(GSRC)%MAPFS(JSRC,ISRC)
                              TMPRL_OM (NRTOT,NLOC_OM) = WX*WY
                           END IF
                        END DO
                     END DO
                     TMPINT_OM(NRTOT,0) = NLOC_OM
                     TMPRL_OM (NRTOT,1:NLOC_OM) = TMPRL_OM(NRTOT,1:NLOC_OM) &
                        / WTOT
 
                  END IF ! (if OLD_METHOD)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!...End block of code for computing weights using old method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!...Begin block of code for "computing" weights using new method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
! Notes: Weights have already been computed by SCRIP.
!        We just need to transfer them to TMPINT and TMPRL
 
 
! Test output
 
! Notes: check here that we are sufficiently dimensioned.
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!...End block of code for "computing" weights using new method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!...Begin block of code that is just for testing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
                  IF (DO_CHECKING)THEN
! compare with the values that we got using the old method
                     if (NLOC.NE.NLOC_OM) THEN
                           IF ( IMPROC.EQ.NMPERR )WRITE (MDSE,'(/1A,2(I8))') &
                           ' *** ERROR WMGHGH: NLOC,NLOC_OM  = ',NLOC,NLOC_OM
                           CALL EXTCDE ( 999 )
                     END IF
                     ISTOP=0
                     ICOUNT=0
                     DO IPNT=1,NLOC
                        DO IPNT2=1,NLOC
                           IF (TMPINT_OM(NRTOT,IPNT).EQ.TMPINT(NRTOT,IPNT2))THEN
! we found our point
                              ICOUNT=ICOUNT+1
                              IF(ABS(TMPRL_OM(NRTOT,IPNT)-TMPRL(NRTOT,IPNT2)) &
                                 .GT.4.0e-5)then
                                 IF ( IMPROC.EQ.NMPERR )WRITE &
                                      (MDSE,'(/2A,2(F12.5))') &
                           ' *** ERROR WMGHGH: ', &
                           ' *** TMPRL_OM(NRTOT,IPNT),TMPRL(NRTOT,IPNT2) = ', &
                                      TMPRL_OM(NRTOT,IPNT),TMPRL(NRTOT,IPNT2)
                                 ISTOP=1
                              END IF
                           END IF
                        END DO
                     END DO
                     IF(ICOUNT.NE.NLOC)THEN
                        IF ( IMPROC.EQ.NMPERR )WRITE (MDSE,'(/1A,2(I8))') &
                         ' *** ERROR WMGHGH: ICOUNT,NLOC = ',ICOUNT,NLOC
                        ISTOP=1
                     END IF
                     IF(ISTOP.EQ.1)THEN
                       CALL EXTCDE ( 999 )
                     END IF
 
                  END IF ! (if both grids are regular grids)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!...End block of code that is just for testing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
                  NROK   = 0
 
               END DO LOWRANK_I !  DO IDST=1, NX
            END DO LOWRANK_J ! DO JDST=1, NY
 
 
! Notes: We are done with the counting. If we didn't use SCRIP to get
!     TMPINT, TMPRL, then we need to set them using the _OM variables.
 
            IF(.NOT.LSCRIP)THEN
               TMPINT=TMPINT_OM
               TMPRL=TMPRL_OM
            END IF
 
!
! 2.c Set up masks based on stencil width of scheme and inferred land
! 2.c.1 Inferred land
!
!!HT: Inferred land from INFLND is added to MAPSTA / MAPST2
!!HT:
            MAPST2 = MAPST2 - 4*MOD(MAPST2/4,2)
            MAPST2 = MAPST2 + 4*INFLND
            DO IDST=1, NX
               DO JDST=1, NY
                  IF ( MAPST2(JDST,IDST).GT.0 ) MAPSTA(JDST,IDST) =   &
                       - ABS(MAPSTA(JDST,IDST))
               END DO
            END DO
!
! 2.c.2 Masking
!
!!HT: This is masking in the low-res grid to identify where the grid
!!HT: is covered by high-res grids, and far enough away from the
!!HT: high-res grid edges so that no dynamic computations are needed
!!HT: in the low-res grid.
!!HT:
            ALLOCATE ( STMASK(NY,0:NX+1), MASKI(NY,NX), STAT=ISTAT )
            CHECK_ALLOC_STATUS ( ISTAT )
            IF ( MDATAS(GDST)%MSKINI ) THEN
                DEALLOCATE ( MDATAS(GDST)%MAPMSK, STAT=ISTAT )
                CHECK_DEALLOC_STATUS ( ISTAT )
              END IF
            ALLOCATE ( MDATAS(GDST)%MAPMSK(NY,NX), STAT=ISTAT )
            CHECK_ALLOC_STATUS ( ISTAT )
            MAPMSK => MDATAS(GDST)%MAPMSK
            MDATAS(GDST)%MSKINI = .TRUE.
 
            MAPMSK = MOD(MAPST2/8,2)
            MAPST2 = MAPST2 - 8*MAPMSK
 
!!HT: STMASK (logical) is used to start this up. We first use the point
!!HT:        MAPTST that have been marked as used for boundary points in
!!HT:        the corrsponding high-res grids.
!!HT: NIT sets the stencil width of the propagation scheme, used to see
!!HT:        how far we need to move in from the boundary points of
!!HT:        the high-res grid to reach the area in the low-res grid
!!HT:        where we do not need to compute.
 
            STMASK(:,1:NX) = MAPTST .LT. 0
            STMASK(:,0) = STMASK(:,NX)
            STMASK(:,NX+1) = STMASK(:,1)
 
            NIT    = ( 1 + INT(DTMAX/DTCFL-0.001) ) * 3
 
            IDSTLA=2
            IDSTHA=NX-1
 
! notes....bug fix: in official release 3.14, the if-then below
!  was missing. This would produce incorrect results for a global grid that
!  had a higher rank grid on the branch cut (180 to -180 or 360 to 0). See
!  treatment of STMASK after the if-then statement. There, it is clear that
!  it was intended that MASKI be available for i=1 and i=nx, ... but it wasn't
!  available. Symptoms of bug: when using "T T" for masking options, a strip
!  of land would be placed along the i-column just east of the branch cut.
!  This would be seen in the global (low rank) grid.
 
            IF ( ICLOSE.NE.ICLOSE_NONE ) THEN
               IDSTLA=1
               IDSTHA=NX
            END IF
 
            DO JTMP=1, NIT
               MASKI  = .FALSE.
               DO IDST=IDSTLA,IDSTHA
                  DO JDST=2, NY-1
                     IF ( .NOT. STMASK(JDST,IDST) .AND. (                     &
                        STMASK(JDST+1,IDST+1) .OR.  STMASK(JDST+1,IDST  ) .OR. &
                        STMASK(JDST+1,IDST-1) .OR.  STMASK(JDST  ,IDST-1) .OR. &
                        STMASK(JDST-1,IDST-1) .OR.  STMASK(JDST-1,IDST  ) .OR. &
                        STMASK(JDST-1,IDST+1) .OR.  STMASK(JDST  ,IDST+1) ) )  &
                          MASKI(JDST,IDST) = .TRUE.
                  END DO
               END DO
               STMASK(:,1:NX) = STMASK(:,1:NX) .OR. MASKI
               STMASK(:,0) = STMASK(:,NX)
               STMASK(:,NX+1) = STMASK(:,1)
            END DO
 
!!HT: Loop over all point from which low-res grid gets data from
!!HT: high-res grid(s). Comparing to STMASK shows which points can be
!!HT: masked out for computation.
!!HT:
!!HT: MAPMSK is stored in WMMDATMD for use in wave model.
 
            DO ILOC=1, NRTOT
               IDST     = TMPINT(ILOC,-4)
               JDST     = TMPINT(ILOC,-3)
               TMPLOG(ILOC) = STMASK(JDST,IDST)
               IF ( .NOT. STMASK(JDST,IDST) ) THEN
                  MAPMSK(JDST,IDST) = 1
                  IF ( FLGHG1 ) MAPSTA(JDST,IDST) = -ABS(MAPSTA(JDST,IDST))
                  MAPTST(JDST,IDST) = 99
               END IF
            END DO
 
            IF ( FLGHG1 ) MAPST2 = MAPST2 + 8*MAPMSK
 
            DEALLOCATE ( STMASK, MASKI, STAT=ISTAT )
            CHECK_DEALLOC_STATUS ( ISTAT )
 
!!HT: Now that all temporary data is stored, and all mosks are set, all
!!HT: can be put from temporaty storage in permanent storage.
!!HT:
!!HT: Should require no modifications for newer grids ...
!!HT:    ... unless more data is needed than for old grids .....
 
!
! 2.d Set up mapping for staging data
! 2.d.1 Set counters / required array sizes
!
            I1     = 0
            I2     = 0
            I3     = 0
            I4     = 0
 
            DO ILOC=1, NRTOT
 
               JJ     = TMPINT(ILOC,-1)
               HGSTGE(GDST,JJ)%NTOT = HGSTGE(GDST,JJ)%NTOT + 1
               ISEA   = TMPINT(ILOC,-2)
              JSEA   = 1 + (ISEA-1)/NAPROC
              ISPROC = ISEA - (JSEA-1)*NAPROC + CROOT - 1
!
               I1(JJ,ISPROC) = I1(JJ,ISPROC) + 1
               IF ( TMPLOG(ILOC) ) I2(JJ,ISPROC) = I2(JJ,ISPROC) + 1
               IF ( IMPROC .EQ. ISPROC ) THEN
                 HGSTGE(GDST,JJ)%NSMX = MAX(HGSTGE(GDST,JJ)%NSMX,TMPINT(ILOC,0))
               END IF
 
               DO JR=1, TMPINT(ILOC,0)
                  ISEA   = TMPINT(ILOC,JR)
                JSEA   = 1 + (ISEA-1)/OUTPTS(JJ)%NAPROC
                ISPRO2 = ISEA - (JSEA-1)*OUTPTS(JJ)%NAPROC +    &
                                MDATAS(JJ)%CROOT - 1
                  IF ( ISPRO2 .EQ. IMPROC ) THEN
                     HGSTGE(GDST,JJ)%NSND = HGSTGE(GDST,JJ)%NSND + 1
                     IF ( TMPLOG(ILOC) ) HGSTGE(GDST,JJ)%NSN1 =              &
                          HGSTGE(GDST,JJ)%NSN1 + 1
                  END IF
               END DO
!
            END DO
 
            HGSTGE(GDST,:)%NREC = I1(:,IMPROC)
            HGSTGE(GDST,:)%NRC1 = I2(:,IMPROC)
!
! 2.d.2 ALLOCATE (DEALLOCATE in section 0 as needed)
!
            DO GSRC=1, NRGRD
               IF ( HGSTGE(GDST,GSRC)%NREC .GT. 0 ) THEN
                   ALLOCATE (                                          &
                     HGSTGE(GDST,GSRC)%LJSEA (HGSTGE(GDST,GSRC)%NREC), &
                     HGSTGE(GDST,GSRC)%NRAVG (HGSTGE(GDST,GSRC)%NREC), &
                     HGSTGE(GDST,GSRC)%IMPSRC(HGSTGE(GDST,GSRC)%NREC,  &
                     HGSTGE(GDST,GSRC)%NSMX),                          &
                     HGSTGE(GDST,GSRC)%ITAG  (HGSTGE(GDST,GSRC)%NREC,  &
                     HGSTGE(GDST,GSRC)%NSMX),                          &
                     HGSTGE(GDST,GSRC)%WGTH  (HGSTGE(GDST,GSRC)%NREC,  &
                     HGSTGE(GDST,GSRC)%NSMX),                          &
                     HGSTGE(GDST,GSRC)%SHGH  (SGRDS(GSRC)%NSPEC,       &
                     HGSTGE(GDST,GSRC)%NSMX,                           &
                     HGSTGE(GDST,GSRC)%NREC), STAT=ISTAT               )
                   CHECK_ALLOC_STATUS ( ISTAT )
               END IF
               IF ( HGSTGE(GDST,GSRC)%NSND .GT. 0 ) THEN
                   ALLOCATE ( HGSTGE(GDST,GSRC)%ISEND (HGSTGE(GDST,GSRC)%NSND,5), &
                              STAT=ISTAT )
                   CHECK_ALLOC_STATUS ( ISTAT )
               END IF
               HGSTGE(GDST,GSRC)%INIT = .TRUE.
            END DO
!
! 2.d.3 Fill allocated arrays
!
            FLGREC = .TRUE.
            I2     = I1 + 1
            I1     = 0
            I4     = HGSTGE(GDST,:)%NSND + 1
            I3     = 0
 
            DO ILOC=1, NRTOT
 
               ISEA   = TMPINT(ILOC,-2)
               JJ     = TMPINT(ILOC,-1)
               NR0    = TMPINT(ILOC, 0)
              JSEA   = 1 + (ISEA-1)/NAPROC
              ISPROC = ISEA - (JSEA-1)*NAPROC + CROOT - 1
              FLGREC = ISPROC .EQ. IMPROC
!
               IF ( TMPLOG(ILOC) ) THEN
                  I1(JJ,ISPROC) = I1(JJ,ISPROC) + 1
                  IREC          = I1(JJ,ISPROC)
               ELSE
                  I2(JJ,ISPROC) = I2(JJ,ISPROC) - 1
                  IREC          = I2(JJ,ISPROC)
               END IF
 
               IF ( FLGREC ) THEN
                  HGSTGE(GDST,JJ)%LJSEA(IREC) = JSEA
                  HGSTGE(GDST,JJ)%NRAVG(IREC) = NR0
                  HGSTGE(GDST,JJ)%WGTH(IREC,:NR0) = TMPRL(ILOC,1:NR0)
                  HGSTGE(GDST,JJ)%ITAG(IREC,:NR0) = LTAG(:NR0)
               END IF
 
               DO IJ=1, NR0
 
                  ISEA   = TMPINT(ILOC,IJ)
                JSEA   = 1 + (ISEA-1)/OUTPTS(JJ)%NAPROC
                ISPRO2 = ISEA - (JSEA-1)*OUTPTS(JJ)%NAPROC      &
                              + MDATAS(JJ)%CROOT - 1
!
                  IF ( FLGREC ) HGSTGE(GDST,JJ)%IMPSRC(IREC,IJ) = ISPRO2
 
                  IF ( ISPRO2 .EQ. IMPROC ) THEN
                     IF ( TMPLOG(ILOC) ) THEN
                        I3(JJ) = I3(JJ) + 1
                        ISND   = I3(JJ)
                     ELSE
                        I4(JJ) = I4(JJ) - 1
                        ISND   = I4(JJ)
                     END IF
                     HGSTGE(GDST,JJ)%ISEND(ISND,1) = JSEA
                  HGSTGE(GDST,JJ)%ISEND(ISND,2) = ISPROC
                     HGSTGE(GDST,JJ)%ISEND(ISND,3) = IREC
                     HGSTGE(GDST,JJ)%ISEND(ISND,4) = IJ
                  HGSTGE(GDST,JJ)%ISEND(ISND,5) = LTAG(IJ)
                  END IF
 
               END DO
!
              LTAG   = LTAG + NR0
              LTAG0  = LTAG0 + NR0
!
            END DO
!
! 2.e Adjust FLAGST using MAPTST
!
            DO ISEA=1, NSEA
               IDST     = MAPSF(ISEA,1)
               JDST     = MAPSF(ISEA,2)
               IF ( MAPTST(JDST,IDST) .GT. 0 ) FLAGST(ISEA) = .NOT. FLGHG1
            END DO
!
! 2.f Test output map
!
            DEALLOCATE ( MAPTST, INFLND, STAT=ISTAT )
            CHECK_DEALLOC_STATUS ( ISTAT )
!
! 2.g Test output receiving
!
! 2.h Test output sending
!
! 2.i Final clean up
!
            DEALLOCATE ( IDSTL, IDSTH, JDSTL, JDSTH, GRIDOK,  BDIST, &
                         TMPINT, TMPRL, TMPLOG, STAT=ISTAT )
            CHECK_DEALLOC_STATUS ( ISTAT )
 
            IF (OLD_METHOD) THEN
                DEALLOCATE ( BDIST_OM, TMPINT_OM, TMPRL_OM, STAT=ISTAT )
                CHECK_DEALLOC_STATUS ( ISTAT )
              END IF
 
            DEALLOCATE ( LTAG, STAT=ISTAT )
            CHECK_DEALLOC_STATUS ( ISTAT )
!
 
! Notes: We are done with this dst (low rank) grid, so we deallocate ALLWGTS .
!        This is important because ALLWGTS will be allocated again for the next
!        dst grid.
 
 
         END IF !  IF ( GRDHGH(GDST,0) ...
      END DO LOWRANK_GRID
 
!
      DEALLOCATE ( I1, I2, I3, I4, STAT=ISTAT )
      CHECK_DEALLOC_STATUS ( ISTAT )
      DEALLOCATE ( NX_BEG, NX_END, STAT=ISTAT )
      CHECK_DEALLOC_STATUS ( ISTAT )
 
!
! 2.j Test output counters
!
 
      RETURN
!
! Formats
!
 1000 FORMAT (/' *** WAVEWATCH III ERROR IN WMGHGH : *** '/           &
               '     GRDHGH NOT YET ALLOCATED, CALL WMGLOW FIRST'/)
 1020 FORMAT (/' *** WAVEWATCH III ERROR IN WMGHGH : *** '/           &
               '     TMPINT AND TMPRL TOO SMALL (w/out SCRIP)'/)
 1021 FORMAT (/' *** WAVEWATCH III ERROR IN WMGHGH : *** '/           &
               '     TMPINT AND TMPRL TOO SMALL (w/SCRIP) '/)
!
!/
!
!/
!/ End of WMGHGH ----------------------------------------------------- /
!/
      END SUBROUTINE WMGHGH
!/ ------------------------------------------------------------------- /
      SUBROUTINE WMGEQL
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         10-Dec-2014 !
!/                  +-----------------------------------+
!/
!/    24-Apr-2006 : Origination.                        ( version 3.09 )
!/    23-Dec-2006 : Adding group test.                  ( version 3.10 )
!/    28-Dec-2006 : Simplify NIT for partial comm.      ( version 3.10 )
!/    22-Jan-2007 : Add saving og NAVMAX.               ( version 3.10 )
!/    02-Feb-2007 : Setting FLAGST for replaced points. ( version 3.10 )
!/    15-Feb-2007 : Tweaking MAPODI algorithm in WMGEQL.( version 3.10 )
!/    11-Apr-2008 : Big fix active edges (MAPSTA=2)     ( version 3.13 )
!/    14-Apr-2008 : Big fix for global grids.           ( version 3.13 )
!/    20-May-2009 : Linking FLAGST and FLGHG1.          ( version 3.14 )
!/    30-Oct-2009 : Implement run-time grid selection.  ( version 3.14 )
!/                  (W. E. Rogers & T. J. Campbell, NRL)
!/    06-Dec-2010 : Change from GLOBAL (logical) to ICLOSE (integer) to
!/                  specify index closure for a grid.   ( version 3.14 )
!/                  (T. J. Campbell, NRL)
!/    23-Dec-2010 : Fix HPFAC and HQFAC by including the COS(YGRD)
!/                  factor with DXDP and DXDQ terms.    ( version 3.14 )
!/                  (T. J. Campbell, NRL)
!/    05-Aug-2013 : Change PR2/3 to UQ/UNO in distances.( version 4.12 )
!/    10-Dec-2014 : Add checks for allocate status      ( version 5.04 )
!/
!  1. Purpose :
!
!     Determine relations to same ranked grids for each grid.
!
!  2. Method :
!
!     Cross mapping of grid points, determine boundary distance data
!     and interpolation weights.
!
!  3. Parameters :
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SETG, W3SETO, WMSETM
!                Subr. W3GDATMD Manage data structures.
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!      EXTCDE    Subr.   Id.    Program abort.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!  6. Error messages :
!
!  7. Remarks :
!
!     - In looking for compatable boundary points in overlapping grids
!       two assumptions hav been made.
!        a) No active boundary points exist in global grids.
!        b) For a lower resolution grid an expanded sewarch area is
!           required for corresponding active grid points. By limiting
!           the resolution ratio to 2, only one extra grid point needs
!           to be considered (JXL2 versus JXL etc.).
!
!  8. Structure :
!
!  9. Switches :
!
!     !/PRn  Propagation scheme.
!
!     !/O12  Removed boundary points output (central).
!     !/O13  Removed boundary points output (edge).
 
!     !/S    Enable subroutine tracing.
!     !/T    Enable test output.
!     !/T5   Detailed test output 'receiving'.
!     !/T6   Detailed test output 'sending'.
!     !/T7   Detailed test output all.
!
!     !/MPI  Distribbuted memory management.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE CONSTANTS
      USE W3GDATMD
      USE W3ODATMD
      USE W3ADATMD
      USE WMMDATMD
!
      USE W3SERVMD, ONLY: EXTCDE
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                  :: I, J, IX, IXL, IXH, IY, IYL, IYH,   &
                                  JX, JXL, JXH, JXL2, JXH2,           &
                                  JY, JYL, JYH, JYL2, JYH2,           &
                                  NR, NT, NA, NTL, JJ, NIT, NG, NOUT, &
                                  ISEA, JSEA, ISPROC, ITAG, TGRP,     &
                                  EXTRA, IP, NP
      INTEGER, ALLOCATABLE     :: MAP3D(:,:,:), NREC(:), NSND(:),     &
                                  NTPP(:), MAPOUT(:,:)
      REAL                     :: FACTOR, XSL, XSH, YSL, YSH, XA, YA, &
                                  XR, YR, RX(2), RY(2), STX, STY,     &
                                  STXY, NEWVAL, WGTH
      REAL, PARAMETER          :: TODO   = -9.99E25
      REAL, PARAMETER          :: ODIMAX = 25.
      REAL, PARAMETER          :: FACMAX = 2.001
      REAL, ALLOCATABLE        :: WGT3D(:,:,:)
      LOGICAL                  :: CHANGE, XEXPND, YEXPND
      LOGICAL, ALLOCATABLE     :: SHRANK(:,:), DOGRID(:),             &
                                  MASKA(:,:), MASKI(:,:)
!
      TYPE STORE
        INTEGER                :: NTOT, NFIN
        INTEGER, POINTER       :: IX(:), IY(:), NAV(:), ISS(:,:),     &
                                  JSS(:,:), IPS(:,:), ITG(:,:)
        REAL, POINTER          :: AWG(:,:)
        LOGICAL, POINTER       :: FLA(:)
        LOGICAL                :: INIT
      END TYPE STORE
!
      TYPE(STORE), ALLOCATABLE :: STORES(:,:)
!/
!
! -------------------------------------------------------------------- /
! 0.  Initializations
!
 
      ALLOCATE ( SHRANK(NRGRD,NRGRD), STORES(NRGRD,NRGRD),            &
                 DOGRID(NRGRD), STAT=ISTAT )
      CHECK_ALLOC_STATUS ( ISTAT )
!
      SHRANK = .FALSE.
!
      DO I=1, NRGRD
 
        DO J=1, NRGRD
          STORES(I,J)%INIT = .FALSE.
          STORES(I,J)%NTOT = 0
          STORES(I,J)%NFIN = 0
          END DO
        END DO
!
      IF ( FLAGLL ) THEN
        FACTOR = RADIUS * DERA
!notes: was FACTOR = RADIUS / 360. (I don't know where this came from....
!       ...maybe it was supposed to be CIRCUMFERENCE/360)
        ELSE
          FACTOR = 1.
        END IF
      ITAG   = 0
!
! -------------------------------------------------------------------- /
! 1.  Grid point relations and temp data storage
! 1.a Outer loop over all grids
!
      DO I=1, NRGRD
 
!
! 1.b Find grids with same rank
!
        NR     = 0
!
        DO J=1, NRGRD
 
          IF ( GRANK(I).NE.GRANK(J) .OR. I.EQ.J ) CYCLE
          SHRANK(I,J) = .TRUE.
          NR          = NR + 1
          END DO
!
        CALL W3SETG ( I, MDSE, MDST )
!
        DOGRID(I) = NR .GT. 0
 
!..notes: we will reach this point even if there are no equal rank grids
 
        IF ( NR .EQ. 0 ) CYCLE
 
!..notes: we will not reach this point if are no equal rank grids. that makes it a good place to check against grid type
 
        IF ( ICLOSE .EQ. ICLOSE_TRPL ) THEN
           IF ( IMPROC.EQ.NMPERR ) WRITE(MDSE,*)'SUBROUTINE WMGEQL IS'// &
           ' NOT YET ADAPTED FOR TRIPOLE GRIDS. STOPPING NOW.'
           CALL EXTCDE ( 1 )
        END IF
 
! Unresolved bug: this triggers even for 2 irregular grids that are not overlapping!
!   We should only be checking for cases of 2 irregular grids of equal rank that
!   are overlapping. Unfortunately, at this point in the routine, we don't know
!   whether they are overlapping...requires more code to do this, since all code
!   in this routine is for regular grids. Fortunately, there is really no
!   disadvantage to making the two irregular grids to be different rank using
!   ww3_multi.inp
 
        IF ( GTYPE .EQ. UNGTYPE ) THEN
           IF ( IMPROC.EQ.NMPERR )WRITE (MDSE,'(/3A)') ' *** ERROR ',  &
                'WMGEQL: UNSTRUCTURED GRID SUPPORT NOT YET ',          &
                'IMPLEMENTED ***'
           CALL EXTCDE ( 999 )
        END IF
        IF ( GTYPE .EQ. CLGTYPE ) THEN
           IF ( IMPROC.EQ.NMPERR )WRITE (MDSE,'(/3A)') ' *** ERROR ',  &
                'WMGEQL: CURVILINEAR GRID SUPPORT NOT IMPLEMENTED ',   &
                'FOR NRGRD > 1 ***'
           CALL EXTCDE ( 999 )
        END IF
 
!
! 1.c Fill TMPMAP with raw relational data
!
! 1.c.1 Loop over grids, select same rank
!
        DO J=1, NRGRD
 
          IF ( .NOT. SHRANK(I,J) ) CYCLE
!
! 1.c.2 Determine shared area
!       Don't even try for X in LLG
!
 
!       Note: Check is against FLAGLL. Would it be more appropriate
!                   to check against ICLOSE?
          IF ( FLAGLL ) THEN
              IXL    =  1
              IXH    = NX
            ELSE
              XSL    = ( GRIDS(J)%X0 - X0 ) / SX - 0.01
              XSH    = ( GRIDS(J)%X0 + GRIDS(J)%SX*(GRIDS(J)%NX-1)   &
                                     - X0 ) / SX + 0.01
              IXL    = MAX ( 1+NINT(XSL) ,  1 )
              IXH    = MIN ( 1+NINT(XSH) , NX )
            END IF
!
          YSL    = ( GRIDS(J)%Y0 - Y0 ) / SY - 0.01
          YSH    = ( GRIDS(J)%Y0 + GRIDS(J)%SY*(GRIDS(J)%NY-1)   &
                                 - Y0 ) / SY + 0.01
          IYL    = MAX ( 1+NINT(YSL) ,  1 )
          IYH    = MIN ( 1+NINT(YSH) , NY )
!
          NT     = (1+IXH-IXL) * (1+IYH-IYL)
          IF ( NT .EQ. 0 ) CYCLE
!
          STORES(I,J)%INIT = .TRUE.
          ALLOCATE ( STORES(I,J)%IX(NT)   , STORES(I,J)%IY(NT)   ,    &
                     STORES(I,J)%NAV(NT)  , STORES(I,J)%FLA(NT)  ,    &
                     STORES(I,J)%ISS(NT,4), STORES(I,J)%JSS(NT,4),    &
                     STORES(I,J)%IPS(NT,4), STORES(I,J)%ITG(NT,4),    &
                     STORES(I,J)%AWG(NT,4), STAT=ISTAT )
          CHECK_ALLOC_STATUS ( ISTAT )
          STORES(I,J)%NAV  = 0
          STORES(I,J)%FLA  = .FALSE.
          STORES(I,J)%ISS  = 0
          STORES(I,J)%JSS  = 0
          STORES(I,J)%IPS  = 0
          STORES(I,J)%ITG  = 0
          STORES(I,J)%AWG  = 0.
!
! 1.c.3 Loops over shared area
!
          NT     = 0
!
          XEXPND = SX .GT. GRIDS(J)%SX
          YEXPND = SY .GT. GRIDS(J)%SY
!
          DO IX=IXL, IXH
            XA     = X0 + REAL(IX-1)*SX
            IF ( FLAGLL ) THEN
                XR     = 1. + MOD (1080. + XA - GRIDS(J)%X0 , 360. ) &
                                  / GRIDS(J)%SX
              ELSE
                XR     = 1. + (XA-GRIDS(J)%X0) / GRIDS(J)%SX
              END IF
            JXL    = INT(XR)
            JXH    = JXL + 1
            RX(1)     = 1. - MOD(XR,1.)
            IF ( RX(1).GT.0.99 .OR. JXH.EQ.GRIDS(J)%NX+1 ) THEN
                JXH    = JXL
                RX(1)  = 1.
              END IF
            IF ( RX(1).LT.0.01 .OR. JXL.EQ.0 ) THEN
                JXL    = JXH
                RX(1)  = 1.
              END IF
            RX(2)  = 1. - RX(1)
!
            IF ( JXL.LT.1 .OR. JXH.GT.GRIDS(J)%NX ) CYCLE
!
            IF ( XEXPND ) THEN
                JXL2   = MAX ( 1 , JXL-1 )
                JXH2   = MIN ( GRIDS(J)%NX , JXH+1 )
              ELSE
                JXL2   = JXL
                JXH2   = JXH
              END IF
!
            DO IY=IYL, IYH
              YA     = Y0 + REAL(IY-1)*SY
              YR     = 1. + (YA-GRIDS(J)%Y0) / GRIDS(J)%SY
              JYL    = INT(YR)
              JYH    = JYL + 1
              RY(1)  = 1. - MOD(YR,1.)
              IF ( RY(1).GT.0.99 .OR. JYH.EQ.GRIDS(J)%NY+1 ) THEN
                  JYH    = JYL
                  RY(1)  = 1.
                END IF
              IF ( RY(1).LT.0.01 .OR. JYL.EQ.0 ) THEN
                  JYL    = JYH
                  RY(1)  = 1.
                END IF
              IF ( RY(1) .GT. 0.99 ) JYH = JYL
              RY(2)  = 1. - RY(1)
!
              IF ( JYL.LT.1 .OR. JYH.GT.GRIDS(J)%NY ) CYCLE
!
              IF ( YEXPND ) THEN
                  JYL2   = MAX ( 1 , JYL-1 )
                  JYH2   = MIN ( GRIDS(J)%NY , JYH+1 )
                ELSE
                  JYL2   = JYL
                  JYH2   = JYH
                END IF
!
! 1.c.4 Temp storage of raw data
!
              NT     = NT + 1
              NA     = 0
              STORES(I,J)%IX(NT)  = IX
              STORES(I,J)%IY(NT)  = IY
!
              DO JX = JXL, JXH
                DO JY = JYL, JYH
                  IF ( GRIDS(J)%MAPSTA(JY,JX) .NE. 0 ) THEN
                      NA     = NA + 1
                      ITAG   = ITAG + 1
                      WGTH   = RX(1+JX-JXL) * RY(1+JY-JYL)
                      ISEA   = GRIDS(J)%MAPFS(JY,JX)
                      IF ( ISEA .EQ. 0 ) THEN
                          JSEA   = 0
                          ISPROC = 1
                        ELSE
                          JSEA   = 1 + (ISEA-1)/OUTPTS(J)%NAPROC
                          ISPROC = ISEA -                        &
                                      (JSEA-1)*OUTPTS(J)%NAPROC  &
                                           + MDATAS(J)%CROOT - 1
                        END IF
                      STORES(I,J)%AWG(NT,NA) = WGTH
                      STORES(I,J)%ISS(NT,NA) = ISEA
                      STORES(I,J)%JSS(NT,NA) = JSEA
                      STORES(I,J)%IPS(NT,NA) = ISPROC
                      STORES(I,J)%ITG(NT,NA) = ITAG
                    END IF
                  END DO
                END DO
!
              DO JX = JXL2, JXH2
                DO JY = JYL2, JYH2
                  IF ( ABS(GRIDS(J)%MAPSTA(JY,JX)) .EQ. 2 )           &
                      STORES(I,J)%FLA(NT) = .TRUE.
                  END DO
                END DO
!
              WGTH   = SUM ( STORES(I,J)%AWG(NT,1:NA) )
              IF ( WGTH .LT. 0.499 ) THEN
                  NA    = 0
                ELSE
                  STORES(I,J)%AWG(NT,:) = STORES(I,J)%AWG(NT,:) / WGTH
                END IF
!
              STORES(I,J)%NAV(NT) = NA
!
! ... End of loops in 1.c
!
              END DO
            END DO
!
          STORES(I,J)%NTOT = NT
!
          END DO
!
! -------------------------------------------------------------------- /
! 2.  Generate open edge distance maps
! 2.a Base map based on MAPSTA only, time step not included.
!
        ALLOCATE ( MDATAS(I)%MAPODI(NY,NX), STAT=ISTAT )
        CHECK_ALLOC_STATUS ( ISTAT )
        MAPODI => MDATAS(I)%MAPODI
        MAPODI = 0.
!
        DO IX=1, NX
          DO IY=1, NY
            IF ( ABS(MAPSTA(IY,IX)) .EQ. 1 ) THEN
                MAPODI(IY,IX) = TODO
              ELSE IF ( ABS(MAPSTA(IY,IX)) .EQ. 2 ) THEN
                MAPODI(IY,IX) = -2. / SIG(1) * DTMAX
              ELSE
                MAPODI(IY,IX) = -1. / SIG(1) * DTMAX
              END IF
            END DO
          END DO
!
! 2.b Add in cross-grid information from STORES
!
        ALLOCATE ( MASKA(NY,NX), STAT=ISTAT )
        CHECK_ALLOC_STATUS ( ISTAT )
        MASKA  = .FALSE.
!
        DO J=1, NRGRD
          IF ( .NOT. SHRANK(I,J) ) CYCLE
          DO JJ=1, STORES(I,J)%NTOT
            IX     = STORES(I,J)%IX(JJ)
            IY     = STORES(I,J)%IY(JJ)
            IF ( IX.EQ.1 .OR. IX.EQ.NX .OR. IY.EQ.1 .OR. IY.EQ.NY ) THEN
                MASKA(IY,IX) = STORES(I,J)%FLA(JJ) .OR.               &
                               STORES(I,J)%NAV(JJ).EQ.0
                IF ( ABS(MAPSTA(IY,IX)).EQ.2 .AND.                    &
                     .NOT.STORES(I,J)%FLA(JJ) .AND.                   &
                     STORES(I,J)%NAV(JJ).GT.0 ) THEN
                     MAPODI(IY,IX) = 0.
                  END IF
              ELSE
                MASKA(IY,IX) = STORES(I,J)%FLA(JJ)
              END IF
            IF ( MAPSTA(IY,IX).EQ.0 .AND. MAPST2(IY,IX) .EQ.1 .AND.   &
                 STORES(I,J)%NAV(JJ).GT.0 ) MAPODI(IY,IX) = 0.
            END DO
          END DO
!
! 2.c Remove incompatable boundary points
!
        ALLOCATE ( MASKI(NY,NX), STAT=ISTAT )
        CHECK_ALLOC_STATUS ( ISTAT )
        MASKI  = .FALSE.
!
        DO IX=2, NX-1
          DO IY=2, NY-1
            IF ( ABS(MAPSTA(IY,IX)) .EQ. 2 .AND.                      &
                 .NOT. MASKA(IY,IX) .AND. (                           &
                          MAPODI(IY-1,IX  ) .GE. 0. .OR.              &
                          MAPODI(IY+1,IX  ) .GE. 0. .OR.              &
                          MAPODI(IY  ,IX-1) .GE. 0. .OR.              &
                          MAPODI(IY  ,IX+1) .GE. 0. ) ) THEN
                MASKI(IY,IX) = .TRUE.
              END IF
            END DO
          END DO
!
        DEALLOCATE ( MASKA, STAT=ISTAT )
        CHECK_DEALLOC_STATUS ( ISTAT )
!
        DO IX=1, NX
          DO IY=1, NY
            IF ( MASKI(IY,IX) ) MAPODI(IY,IX) = 0.
            END DO
          END DO
!
! 2.d Mask out influenced edge
!
        NIT    = ( 1 + INT(DTMAX/DTCFL-0.001) ) * 3
!
        IF ( ICLOSE.NE.ICLOSE_NONE ) THEN
            IXL    =  1
            IXH    = NX
          ELSE
            IXL    =  2
            IXH    = NX - 1
          END IF
!
        DO J=1, NIT
!
          MASKI  = .FALSE.
!
          DO IX=IXL, IXH
            IF ( IX .EQ. 1 ) THEN
                JXL    = NX
                JXH    =  2
              ELSE IF ( IX .EQ. NX ) THEN
                JXL    = NX - 1
                JXH    = 1
              ELSE
                JXL    = IX - 1
                JXH    = IX + 1
              END IF
!
            DO IY=2, NY-1
              IF ( MAPODI(IY,IX) .EQ. TODO .AND. (                    &
                   MAPODI(IY+1,IX ) .GE. 0. .OR.                      &
                   MAPODI(IY  ,JXL) .GE. 0. .OR.                      &
                   MAPODI(IY-1,IX ) .GE. 0. .OR.                      &
                   MAPODI(IY  ,JXH) .GE. 0. .OR.                      &
                 ( MAPODI(IY+1,JXH) .GE. 0. .AND. .NOT.               &
                 ( MAPSTA(IY+1,IX ) .NE. 1 .AND.                      &
                   MAPSTA(IY  ,JXH) .NE. 1 ) ) .OR.                   &
                 ( MAPODI(IY+1,JXL) .GE. 0. .AND. .NOT.               &
                 ( MAPSTA(IY+1,IX ) .NE. 1 .AND.                      &
                   MAPSTA(IY  ,JXL) .NE. 1 ) ) .OR.                   &
                 ( MAPODI(IY-1,JXL) .GE. 0. .AND. .NOT.               &
                 ( MAPSTA(IY-1,IX ) .NE. 1 .AND.                      &
                   MAPSTA(IY  ,JXL) .NE. 1 ) ) .OR.                   &
                 ( MAPODI(IY-1,JXH) .GE. 0. .AND. .NOT.               &
                 ( MAPSTA(IY-1,IX ) .NE. 1 .AND.                      &
                   MAPSTA(IY  ,JXH) .NE. 1 ) ) ) )                    &
                  MASKI(IY,IX) = .TRUE.
              END DO
!
            END DO
!
          DO IX=IXL, IXH
            DO IY=2, NY-1
              IF ( MASKI(IY,IX) ) MAPODI(IY,IX) = 0.
              END DO
            END DO
!
          END DO
!
! 2.e Compute distances
!
        DO
          MASKI  = .FALSE.
!
          DO IX=IXL, IXH
            IF ( IX .EQ. 1 ) THEN
                JXL    = NX
                JXH    =  2
              ELSE IF ( IX .EQ. NX ) THEN
                JXL    = NX - 1
                JXH    = 1
              ELSE
                JXL    = IX - 1
                JXH    = IX + 1
              END IF
            DO IY=2, NY-1
              IF ( MAPODI(IY,IX) .EQ. TODO .AND. (                    &
                   MAPODI(IY+1,IX ) .GE. 0. .OR.                      &
                   MAPODI(IY-1,IX ) .GE. 0. .OR.                      &
                   MAPODI(IY  ,JXH) .GE. 0. .OR.                      &
                   MAPODI(IY  ,JXL) .GE. 0. .OR.                      &
                 ( MAPODI(IY+1,JXH) .GE. 0. .AND. .NOT.               &
                 ( MAPSTA(IY+1,IX ) .NE. 1 .AND.                      &
                   MAPSTA(IY  ,JXH) .NE. 1 ) ) .OR.                   &
                 ( MAPODI(IY+1,JXL) .GE. 0. .AND. .NOT.               &
                 ( MAPSTA(IY+1,IX ) .NE. 1 .AND.                      &
                   MAPSTA(IY  ,JXL) .NE. 1 ) ) .OR.                   &
                 ( MAPODI(IY-1,JXL) .GE. 0. .AND. .NOT.               &
                 ( MAPSTA(IY-1,IX ) .NE. 1 .AND.                      &
                   MAPSTA(IY  ,JXL) .NE. 1 ) ) .OR.                   &
                 ( MAPODI(IY-1,JXH) .GE. 0. .AND. .NOT.               &
                 ( MAPSTA(IY-1,IX ) .NE. 1 .AND.                      &
                   MAPSTA(IY  ,JXH) .NE. 1 ) ) ) )                    &
                  MASKI(IY,IX) = .TRUE.
              END DO
            END DO
!
          CHANGE = .FALSE.
          DO IY=2, NY-1
            DO IX=IXL, IXH
              IF ( IX .EQ. 1 ) THEN
                  JXL    = NX
                  JXH    =  2
                ELSE IF ( IX .EQ. NX ) THEN
                  JXL    = NX - 1
                  JXH    = 1
                ELSE
                  JXL    = IX - 1
                  JXH    = IX + 1
                END IF
              ISEA   = MAPFS(IY,IX)
              STY    = FACTOR * HQFAC(IY,IX) / ( 0.58 * GRAV )
              STX    = FACTOR * HPFAC(IY,IX) &
                                             / ( 0.58 * GRAV )
              STXY   = SQRT ( STX**2 + STY**2 )
              IF ( MASKI(IY,IX) ) THEN
                  NEWVAL = ODIMAX / SIG(1) * DTMAX
                  IF ( MAPODI(IY+1,IX ).GE.0. .AND. .NOT.             &
                       MASKI (IY+1,IX ) ) NEWVAL = MIN (              &
                       NEWVAL , MAPODI(IY+1,IX )+STY )
                  IF ( MAPODI(IY-1,IX ).GE.0. .AND. .NOT.             &
                       MASKI (IY-1,IX ) ) NEWVAL = MIN (              &
                       NEWVAL , MAPODI(IY-1,IX )+STY )
                  IF ( MAPODI(IY  ,JXH).GE.0. .AND. .NOT.             &
                       MASKI (IY  ,JXH) ) NEWVAL = MIN (              &
                       NEWVAL , MAPODI(IY  ,JXH)+STX)
                  IF ( MAPODI(IY  ,JXL).GE.0. .AND. .NOT.             &
                       MASKI (IY  ,JXL) ) NEWVAL = MIN (              &
                       NEWVAL , MAPODI(IY  ,JXL)+STX)
                  IF ( MAPODI(IY+1,JXH).GE.0. .AND. .NOT.             &
                     ( MAPSTA(IY+1,IX ) .NE. 1 .AND.                  &
                       MAPSTA(IY  ,JXH) .NE. 1 ) ) NEWVAL =           &
                       MIN ( NEWVAL , MAPODI(IY+1,JXH)+STXY)
                  IF ( MAPODI(IY+1,JXL).GE.0. .AND. .NOT.             &
                     ( MAPSTA(IY+1,IX ) .NE. 1 .AND.                  &
                       MAPSTA(IY  ,JXL) .NE. 1 ) ) NEWVAL =           &
                       MIN ( NEWVAL , MAPODI(IY+1,JXL)+STXY)
                  IF ( MAPODI(IY-1,JXL).GE.0. .AND. .NOT.             &
                     ( MAPSTA(IY-1,IX ) .NE. 1 .AND.                  &
                       MAPSTA(IY  ,JXL) .NE. 1 ) ) NEWVAL =           &
                       MIN ( NEWVAL , MAPODI(IY-1,JXL)+STXY)
                  IF ( MAPODI(IY-1,JXH).GE.0. .AND. .NOT.             &
                     ( MAPSTA(IY-1,IX ) .NE. 1 .AND.                  &
                       MAPSTA(IY  ,JXH) .NE. 1 ) ) NEWVAL =           &
                       MIN ( NEWVAL , MAPODI(IY-1,JXH)+STXY)
                  MAPODI(IY,IX) = NEWVAL
                  CHANGE = .TRUE.
                END IF
              END DO
            END DO
!
          IF ( .NOT. CHANGE ) EXIT
          END DO
!
        DO IX=2, NX-1
          DO IY=2, NY-1
            IF ( MAPODI(IY,IX) .EQ. TODO )                            &
                 MAPODI(IY,IX) = 2. * ODIMAX / SIG(1) * DTMAX
            END DO
          END DO
!
        DEALLOCATE ( MASKI, STAT=ISTAT )
        CHECK_DEALLOC_STATUS ( ISTAT )
!
! 2.f Update FLAGST
!
        DO ISEA=1, NSEA
          IX     = MAPSF(ISEA,1)
          IY     = MAPSF(ISEA,2)
          IF ( MAPODI(IY,IX) .EQ. 0. ) FLAGST(ISEA) = .NOT. FLGHG1
          END DO
!
! 2.g Test output
!
! ... End of loop in 1.a
!
        END DO
!
! -------------------------------------------------------------------- /
! 3.  Final data base (full data base, scratched at end of routine)
! 3.a Loop over grids
!
      ALLOCATE ( NREC(NRGRD), NSND(NRGRD), NTPP(NMPROC), STAT=ISTAT )
      CHECK_ALLOC_STATUS ( ISTAT )
!
      DO I=1, NRGRD
        IF ( .NOT. DOGRID(I) ) CYCLE
!
        CALL W3SETG ( I, MDSE, MDST )
        CALL W3SETO ( I, MDSE, MDST )
        CALL WMSETM ( I, MDSE, MDST )
!
        ALLOCATE ( MAP3D(NY,NX,-4:NRGRD), WGT3D(NY,NX,0:NRGRD), STAT=ISTAT )
        CHECK_ALLOC_STATUS ( ISTAT )
        MAP3D  = 0
        WGT3D  = 0.
        NREC   = 0
        NSND   = 0
!
! 3.b Filling MAP3D and WGT3D, as well as NREC and NSND
!
        DO J=1, NRGRD
          IF ( .NOT. SHRANK(I,J) ) CYCLE
          MAPODI => MDATAS(J)%MAPODI
!
          DO JJ=1, STORES(I,J)%NTOT
            IX     = STORES(I,J)%IX(JJ)
            IY     = STORES(I,J)%IY(JJ)
            WGT3D(IY,IX,0) = MDATAS(I)%MAPODI(IY,IX)
            MAP3D(IY,IX,-2) = MAPFS(IY,IX)
            IF ( MAP3D(IY,IX,-2) .NE. 0 ) THEN
                MAP3D(IY,IX,-3) = 1 + (MAP3D(IY,IX,-2)-1)/NAPROC
                MAP3D(IY,IX,-4) = MAP3D(IY,IX,-2) -              &
                     (MAP3D(IY,IX,-3)-1)*NAPROC + CROOT - 1
              END IF
            IF ( WGT3D(IY,IX,0).GE.0. .AND. MAPSTA(IY,IX).NE.0. .AND. &
                 STORES(I,J)%NAV(JJ).GT.0 ) THEN
                WGT3D(IY,IX,J) = ODIMAX / SIG(1) * DTMAX
                DO NA=1, STORES(I,J)%NAV(JJ)
                  JX     = GRIDS(J)%MAPSF(STORES(I,J)%ISS(JJ,NA),1)
                  JY     = GRIDS(J)%MAPSF(STORES(I,J)%ISS(JJ,NA),2)
                  IF ( MAPODI(JY,JX) .GE. 0. ) WGT3D(IY,IX,J) =       &
                         MIN( WGT3D(IY,IX,J) , MAPODI(JY,JX) )
                  END DO
                IF ( WGT3D(IY,IX,J) .GT. 0. ) MAP3D(IY,IX,J) = 1
              END IF
            END DO
!
          STORES(I,J)%NFIN = SUM(MAP3D(:,:,J))
!
          END DO
!
        MAPODI => MDATAS(I)%MAPODI
        DO IX=1, NX
          DO IY=1, NY
            MAP3D(IY,IX, 0) = MAXVAL(MAP3D(IY,IX,1:))
            MAP3D(IY,IX,-1) = SUM(MAP3D(IY,IX,1:))
            IF ( MAP3D(IY,IX,-1) .GT. 0 ) THEN
                IF ( MAPODI(IY,IX)*SIG(1)/DTMAX .GT. 1.5*ODIMAX ) THEN
                    WGT3D(IY,IX, 0:) = 0.
                    MAP3D(IY,IX,-1:) = 0
                  ELSE
                    WGTH   = SUM(WGT3D(IY,IX,:))
                    IF ( WGTH .GT. 1.E-25 ) THEN
                        WGT3D(IY,IX,:) = WGT3D(IY,IX,:) / WGTH
                      ELSE
                        WGT3D(IY,IX,:) = 0.
                      END IF
                    IF ( MAP3D(IY,IX,-4) .EQ. IMPROC ) THEN
                        NREC(I) = NREC(I) + 1
                        DO JJ=1, NRGRD
                          IF ( MAP3D(IY,IX,JJ) .GT. 0 )               &
                               NREC(JJ) = NREC(JJ) + 1
                            END DO
                      END IF
                  END IF
              END IF
            END DO
          END DO
!
        DO J=1, NRGRD
          IF ( .NOT. SHRANK(I,J) ) CYCLE
          DO JJ=1, STORES(I,J)%NTOT
            IX     = STORES(I,J)%IX(JJ)
            IY     = STORES(I,J)%IY(JJ)
            IF ( MAP3D(IY,IX,J) .NE. 0 ) THEN
                DO NA=1, STORES(I,J)%NAV(JJ)
                  IF ( STORES(I,J)%IPS(JJ,NA) .EQ. IMPROC )           &
                      NSND(J) = NSND(J) + 1
                  END DO
              END IF
            END DO
          END DO
!
        NG     = MAXVAL(MAP3D(:,:,-1))
        NTL    = SUM(MAP3D(:,:,0))
!
! 3.c Check for points with all ODI = 0
!
        MAPODI => MDATAS(I)%MAPODI
        NOUT   = 0
!
        JXL    = NX
        JXH    =  1
        JYL    = NY
        JYH    =  1
!
        ALLOCATE ( MAPOUT(NY,NX), STAT=ISTAT )
        CHECK_ALLOC_STATUS ( ISTAT )
        MAPOUT = MAPSTA
!
        DO IX=1, NX
          DO IY=1, NY
            IF ( ABS(MAPSTA(IY,IX)).EQ. 1 .AND.                       &
                 MAPODI(IY,IX) .EQ. 0.    .AND.                       &
                 MAP3D(IY,IX,-1) .EQ. 0 ) THEN
                NOUT = NOUT + 1
                IF ( IMPROC.EQ.NMPERR .AND. NOUT.EQ. 1 )              &
                     WRITE(MDSE,*) ' '
                IF ( IMPROC.EQ.NMPERR .AND. NOUT.LE.25 )              &
                     WRITE(MDSE,1001) I, IX, IY
                IF ( IMPROC.EQ.NMPERR .AND. NOUT.EQ.25 )              &
                     WRITE(MDSE,1006)
                JXL    = MIN ( IX, JXL )
                JXH    = MAX ( IX, JXH )
                JYL    = MIN ( IY, JYL )
                JYH    = MAX ( IY, JYH )
                MAPOUT(IY,IX) = 999
              END IF
            END DO
          END DO
!
! 3.d Test and error output
!
        IF ( NOUT .GT. 0 ) THEN
            IF ( IMPROC.EQ.NMPERR ) THEN
                WRITE(MDSE,1000) I, NOUT
                EXTRA  = 2
                JXL = MAX (  1, JXL - EXTRA )
                JXH = MIN ( NX, JXH + EXTRA )
                JYL = MAX (  1, JYL - EXTRA )
                JYH = MIN ( NY, JYH + EXTRA )
                WRITE (MDSE,1002) JXL, JXH, JYL, JYH
                NP     = 1 + (JXH-JXL)/65
                DO IP=1, NP
                  IXL    = JXL + (IP-1)*65
                  IXH    = MIN( NX, IXL+64 )
                  WRITE (MDSE,1005) IXL, IXH
                  WRITE (MDSE,1003) 'STATUS MAP MAPSTA'
                  DO IY=JYH, JYL, -1
                    WRITE (MDSE,1004) MAPSTA(IY,IXL:IXH)
                    END DO
                  WRITE (MDSE,1003) 'MISSING POINTS IN MAPSTA (**)'
                  DO IY=JYH, JYL, -1
                    WRITE (MDSE,1004) MAPOUT(IY,IXL:IXH)
                    END DO
                  WRITE (MDSE,1003) 'OPEN BOUND. DISTANCE MAP MAPODI'
                  DO IY=JYH, JYL, -1
                   WRITE (MDSE,1004)                                 &
                          NINT(MAPODI(IY,IXL:IXH)*SIG(1)/DTMAX)
                    END DO
                  WRITE (MDSE,1003) 'GRID COVERAGE MAP MAP3D'
                  DO IY=JYH, JYL, -1
                    WRITE (MDSE,1004) MAP3D(IY,IXL:IXH,-1)
                    END DO
                  WRITE (MDSE,*)
                  END DO
              END IF
            CALL EXTCDE (1000)
          END IF
!
        DEALLOCATE ( MAPOUT, STAT=ISTAT )
        CHECK_DEALLOC_STATUS ( ISTAT )
!
! -------------------------------------------------------------------- /
! 4.  Save data base as needed in EQSTGE
!
! 4.a   ALLOCATE storage
! 4.a.1 Local counters, weights and sea counters (grid 'I')
!
        IF ( EQSTGE(I,I)%NREC .NE. 0 ) THEN
            DEALLOCATE (EQSTGE(I,I)%ISEA , EQSTGE(I,I)%JSEA ,         &
                        EQSTGE(I,I)%WGHT, STAT=ISTAT )
            CHECK_DEALLOC_STATUS ( ISTAT )
            EQSTGE(I,I)%NREC = 0
          END IF
!
        IF ( NREC(I) .GT. 0 ) THEN
            ALLOCATE ( EQSTGE(I,I)%ISEA(NREC(I))  ,                   &
                       EQSTGE(I,I)%JSEA(NREC(I))  ,                   &
                       EQSTGE(I,I)%WGHT(NREC(I)), STAT=ISTAT )
            CHECK_ALLOC_STATUS ( ISTAT )
            EQSTGE(I,I)%NREC = NREC(I)
          END IF
!
! 4.a.1 Local counters, arrays weights etc.  (grid 'J' receive)
!
        DO J=1, NRGRD
          IF ( I .EQ. J ) CYCLE
          EQSTGE(I,I)%NTOT = STORES(I,J)%NFIN
!
          IF ( EQSTGE(I,J)%NREC .NE. 0 ) THEN
              DEALLOCATE ( EQSTGE(I,J)%ISEA , EQSTGE(I,J)%JSEA ,      &
                           EQSTGE(I,J)%WGHT , EQSTGE(I,J)%SEQL ,      &
                           EQSTGE(I,J)%NAVG , EQSTGE(I,J)%WAVG ,      &
                           EQSTGE(I,J)%RIP  , EQSTGE(I,J)%RTG, STAT=ISTAT )
              CHECK_DEALLOC_STATUS ( ISTAT )
              EQSTGE(I,J)%NREC   = 0
              EQSTGE(I,J)%NAVMAX = 1
            END IF
!
          IF ( NREC(J) .GT. 0 ) THEN
              NA     =  MAXVAL ( STORES(I,J)%NAV(1:STORES(I,J)%NTOT) )
              EQSTGE(I,J)%NAVMAX = NA
              ALLOCATE ( EQSTGE(I,J)%ISEA(NREC(J))  ,                 &
                         EQSTGE(I,J)%JSEA(NREC(J))  ,                 &
                         EQSTGE(I,J)%WGHT(NREC(J))  ,                 &
                         EQSTGE(I,J)%SEQL(SGRDS(J)%NSPEC,NREC(J),NA), &
                         EQSTGE(I,J)%NAVG(NREC(J))  ,                 &
                         EQSTGE(I,J)%WAVG(NREC(J),NA),                &
                         EQSTGE(I,J)%RIP(NREC(J),NA),                 &
                         EQSTGE(I,J)%RTG(NREC(J),NA), STAT=ISTAT )
              CHECK_ALLOC_STATUS ( ISTAT )
              EQSTGE(I,J)%NREC = NREC(J)
            END IF
!
          IF ( EQSTGE(I,J)%NSND .NE. 0 ) THEN
              DEALLOCATE ( EQSTGE(I,J)%SIS , EQSTGE(I,J)%SJS ,        &
                           EQSTGE(I,J)%SI1 , EQSTGE(I,J)%SI2 ,        &
                           EQSTGE(I,J)%SIP , EQSTGE(I,J)%STG, STAT=ISTAT )
              CHECK_DEALLOC_STATUS ( ISTAT )
              EQSTGE(I,J)%NSND = 0
            END IF
!
          IF ( NSND(J) .GT. 0 ) THEN
              ALLOCATE ( EQSTGE(I,J)%SIS(NSND(J)) ,                   &
                         EQSTGE(I,J)%SJS(NSND(J)) ,                   &
                         EQSTGE(I,J)%SI1(NSND(J)) ,                   &
                         EQSTGE(I,J)%SI2(NSND(J)) ,                   &
                         EQSTGE(I,J)%SIP(NSND(J)) ,                   &
                         EQSTGE(I,J)%STG(NSND(J)), STAT=ISTAT )
              CHECK_ALLOC_STATUS ( ISTAT )
              EQSTGE(I,J)%NSND = NSND(J)
            END IF
!
          END DO
!
! 4.b   Store data in EQSTGE
! 4.b.1 Grid I (JSEA and weight only)
!
        IF ( EQSTGE(I,I)%NREC .GT. 0 ) THEN
            NTL    = 0
            DO IX=1, NX
              DO IY=1, NY
                IF ( MAP3D(IY,IX,0) .EQ. 0 ) CYCLE
                IF ( MAP3D(IY,IX,-4) .NE. IMPROC ) CYCLE
                NTL    = NTL + 1
                EQSTGE(I,I)%ISEA(NTL) = MAP3D(IY,IX,-2)
                EQSTGE(I,I)%JSEA(NTL) = MAP3D(IY,IX,-3)
                EQSTGE(I,I)%WGHT(NTL) = WGT3D(IY,IX,0)
                END DO
              END DO
          END IF
!
! 4.b.2 All other grids, info for receiving
!
        DO J=1, NRGRD
          IF ( .NOT. SHRANK(I,J) ) CYCLE
          IF ( EQSTGE(I,J)%NREC .EQ. 0 ) CYCLE
          NTL    = 0
!
          DO JJ=1, STORES(I,J)%NTOT
            IX     = STORES(I,J)%IX(JJ)
            IY     = STORES(I,J)%IY(JJ)
            IF ( MAP3D(IY,IX,J) .EQ. 0 ) CYCLE
            IF ( MAP3D(IY,IX,-4) .NE. IMPROC ) CYCLE
            NTL    = NTL + 1
            EQSTGE(I,J)%ISEA(NTL) = MAP3D(IY,IX,-2)
            EQSTGE(I,J)%JSEA(NTL) = MAP3D(IY,IX,-3)
            EQSTGE(I,J)%WGHT(NTL) = WGT3D(IY,IX,J)
            NA     = STORES(I,J)%NAV(JJ)
            EQSTGE(I,J)%NAVG(NTL) = NA
            EQSTGE(I,J)%WAVG(NTL,1:NA) = STORES(I,J)%AWG(JJ,1:NA)
            EQSTGE(I,J)%RIP (NTL,1:NA) = STORES(I,J)%IPS(JJ,1:NA)
            EQSTGE(I,J)%RTG (NTL,1:NA) = STORES(I,J)%ITG(JJ,1:NA)
            END DO
!
          END DO
!
! 4.b.3 All other grids, info for sending
!
        DO J=1, NRGRD
          IF ( .NOT. SHRANK(I,J) ) CYCLE
          IF ( EQSTGE(I,J)%NSND .EQ. 0 ) CYCLE
          NTPP   = 0
          NTL    = 0
!
          DO JJ=1, STORES(I,J)%NTOT
            IX     = STORES(I,J)%IX(JJ)
            IY     = STORES(I,J)%IY(JJ)
            IF ( MAP3D(IY,IX,J) .NE. 0 ) THEN
                NTPP(MAP3D(IY,IX,-4)) = NTPP(MAP3D(IY,IX,-4)) + 1
                DO NA=1, STORES(I,J)%NAV(JJ)
                  IF ( STORES(I,J)%IPS(JJ,NA) .EQ. IMPROC ) THEN
                      NTL    = NTL + 1
                      EQSTGE(I,J)%SIS(NTL) = STORES(I,J)%ISS(JJ,NA)
                      EQSTGE(I,J)%SJS(NTL) = STORES(I,J)%JSS(JJ,NA)
                      EQSTGE(I,J)%SI1(NTL) = NTPP(MAP3D(IY,IX,-4))
                      EQSTGE(I,J)%SI2(NTL) = NA
                      EQSTGE(I,J)%SIP(NTL) = MAP3D(IY,IX,-4)
                      EQSTGE(I,J)%STG(NTL) = STORES(I,J)%ITG(JJ,NA)
                    END IF
                  END DO
              END IF
            END DO
!
          END DO
!
! 4.c Detailed test output
!
! ... End of loop started in 3.a
!
        DEALLOCATE ( MAP3D, WGT3D, STAT=ISTAT )
        CHECK_DEALLOC_STATUS ( ISTAT )
        END DO
!
! -------------------------------------------------------------------- /
! 5.  Generate GRDEQL
! 5.a Size of array
!
      NREC   = 0
!
      DO I=1, NRGRD
        DO J=1, NRGRD
          IF ( I.EQ.J .OR. STORES(I,J)%NFIN.EQ.0 ) CYCLE
          NREC(I) = NREC(I) + 1
          END DO
        END DO
!
      NA     = MAXVAL(NREC)
      ALLOCATE ( GRDEQL(NRGRD,0:NA), STAT=ISTAT )
      CHECK_ALLOC_STATUS ( ISTAT )
      GRDEQL = 0
!
! 5.b Fill array
!
      DO I=1, NRGRD
        GRDEQL(I,0) = NREC(I)
        JJ          = 0
        DO J=1, NRGRD
          IF ( I.EQ.J .OR. STORES(I,J)%NFIN.EQ.0 ) CYCLE
          JJ           = JJ + 1
          GRDEQL(I,JJ) = J
          END DO
        END DO
!
! 5.c Resolution test
!
 
      IF ( FLAGLL ) THEN
          FACTOR = 1.
        ELSE
          FACTOR = 1.E-3
        END IF
!
! notes: This resolution test, with FACMAX=2, is pretty strict, so
!        it is not going to be appropriate for irregular grids.
!        We'll just have to trust the judgement of the user in the
!        case of irregular grids. But if we change our minds and do
!        some kind of check for irregular grids, we could make
!        a check against median(HPFAC) and median(HQFAC).
 
      DO I=1, NRGRD
         CALL W3SETG ( I, MDSE, MDST )
         IF ( GTYPE.EQ.RLGTYPE ) THEN
            DO JJ=1, GRDEQL(I,0)
               J      = GRDEQL(I,JJ)
               IF ( GRIDS(J)%GTYPE.EQ.RLGTYPE ) THEN
                  IF ( SX/GRIDS(J)%SX .GT. FACMAX     .OR.                   &
                       SX/GRIDS(J)%SX .LT. 1./FACMAX  .OR.                   &
                       SY/GRIDS(J)%SY .GT. FACMAX     .OR.                   &
                       SY/GRIDS(J)%SY .LT. 1./FACMAX ) THEN
                     IF ( IMPROC.EQ.NMPERR ) WRITE(MDSE,1050) I, FACTOR*SX,  &
                          FACTOR*SY, J, FACTOR*GRIDS(J)%SX, FACTOR*GRIDS(J)%SY
                     CALL EXTCDE ( 1050 )
                  END IF ! IF ( SX/GR ...
               END IF ! IF ( GRIDS(J)%GTYPE...
            END DO ! DO JJ=...
         END IF ! IF ( GTYPE....
      END DO ! DO I=...
!
! 5.d Group number test
!
      DO I=1, NRGRD
        IF ( GRDEQL(I,0) .GE. 2 ) THEN
            TGRP   = GRGRP(GRDEQL(I,1))
            DO J=2, GRDEQL(I,0)
              IF ( GRGRP(GRDEQL(I,J)) .NE. TGRP ) THEN
                  IF ( IMPROC .EQ. NMPERR ) WRITE(MDSE,1051)          &
                                   GRDEQL(I,1), GRGRP(GRDEQL(I,1)),   &
                                   GRDEQL(I,J), GRGRP(GRDEQL(I,J))
                  CALL EXTCDE ( 1051 )
                END IF
              END DO
          END IF
        END DO
!
! -------------------------------------------------------------------- /
! 6.  Final clean up
!
      DO I=1, NRGRD
        DO J=1, NRGRD
          IF ( STORES(I,J)%INIT ) THEN
              DEALLOCATE ( STORES(I,J)%IX  , STORES(I,J)%IY  ,        &
                           STORES(I,J)%NAV , STORES(I,J)%FLA ,        &
                           STORES(I,J)%ISS , STORES(I,J)%JSS ,        &
                           STORES(I,J)%IPS , STORES(I,J)%ITG ,        &
                           STORES(I,J)%AWG , STAT=ISTAT )
              CHECK_DEALLOC_STATUS ( ISTAT )
            END IF
          END DO
        END DO
!
      DEALLOCATE ( SHRANK, STORES, NREC, NSND, NTPP, STAT=ISTAT )
      CHECK_DEALLOC_STATUS ( ISTAT )
!
      RETURN
!
! Formats
!
 1000 FORMAT (/' *** WAVEWATCH III ERROR IN WMGEQL : *** '/           &
               '     UNCOVERED EDGE POINTS FOR GRID',I4,'  (',I6,')'/)
 1001 FORMAT ( '     GRID',I4,'  POINT',2I5,' NOT COVERED (WMGEQL)')
 1002 FORMAT ( '     DIAGNOSTICS IX AND IY RANGE:',4I6)
 1003 FORMAT (/'     SHOWING ',A/)
 1004 FORMAT (2X,65I2)
 1005 FORMAT (/'     SHOWING IX RANGE ',2I6)
 1006 FORMAT ( '        (WILL NOT PRINT ANY MORE UNCOVERED POINTS)')
!
 1020 FORMAT (/' *** WAVEWATCH III WARNING WMGEQL : *** '/            &
       '     REMOVED BOUNDARY POINT FROM OPEN EDGE DISTANCE MAP'/     &
       '     GRID, IX, IY :',3I6)
!
 1050 FORMAT (/' *** WAVEWATCH III ERROR IN WMGEQL : *** '/           &
               '     GRID INCREMENTS TOO DIFFERENT '/                 &
               '     GRID',I4,'   INCREMENTS ',2F8.2/                 &
               '     GRID',I4,'   INCREMENTS ',2F8.2/)
 1051 FORMAT (/' *** WAVEWATCH III ERROR IN WMGEQL : *** '/           &
               '     OVERLAPPING GRIDS NEED TO BE IN SAME GROUP '/    &
               '         GRID',I4,' IN GROUP',I4/                     &
               '         GRID',I4,' IN GROUP',I4/)
!
!/
!/ End of WMGEQL ----------------------------------------------------- /
!/
      END SUBROUTINE WMGEQL
!/ ------------------------------------------------------------------- /
      SUBROUTINE WMRSPC
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         10-Dec-2014 !
!/                  +-----------------------------------+
!/
!/    22-Sep-2005 : Origination.                        ( version 3.08 )
!/    25-Jul-2006 : Point output grid added.            ( version 3.10 )
!/    30-Oct-2009 : Implement run-time grid selection.  ( version 3.14 )
!/                  (W. E. Rogers & T. J. Campbell, NRL)
!/    10-Dec-2014 : Add checks for allocate status      ( version 5.04 )
!/
!  1. Purpose :
!
!     Generate map with flogs for need of spectral grid conversion
!     between models.
!
!  2. Method :
!
!     Test of parameters as introduced before in W3IOBC.
!
!  3. Parameters :
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Sur.  W3SERVMD Subroutine tracing.
!      EXTCDE    Subr.   Id.    Program abort.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      WMINIT    Subr  WMINITMD Multi-grid model initialization.
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
!     !/S    Enable subroutine tracing.
!     !/T    Enable test output
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE W3SERVMD, ONLY: EXTCDE
!
      USE W3GDATMD
      USE W3ODATMD, ONLY: UNIPTS
      USE WMMDATMD
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: I, J, LOW
!/
!
! -------------------------------------------------------------------- /
! 0.  Initializations
!
      IF ( UNIPTS ) THEN
          LOW    = 0
        ELSE
          LOW    = 1
        END IF
      IF ( .NOT. ALLOCATED(RESPEC) ) THEN
          ALLOCATE ( RESPEC(LOW:NRGRD,LOW:NRGRD), STAT=ISTAT )
          CHECK_ALLOC_STATUS ( ISTAT )
        END IF
      RESPEC = .FALSE.
!
! -------------------------------------------------------------------- /
! 1.  Fill map with flags
!
      DO I=LOW, NRGRD
        DO J=I+1, NRGRD
          RESPEC(I,J) = SGRDS(I)%NK    .NE. SGRDS(J)%NK     .OR.      &
                        SGRDS(I)%NTH   .NE. SGRDS(J)%NTH    .OR.      &
                        SGRDS(I)%XFR   .NE. SGRDS(J)%XFR    .OR.      &
                        SGRDS(I)%FR1   .NE. SGRDS(J)%FR1    .OR.      &
                        SGRDS(I)%TH(1) .NE. SGRDS(J)%TH(1)
          RESPEC(J,I) = RESPEC(I,J)
          END DO
        END DO
!
! -------------------------------------------------------------------- /
! 2.  Test output
!
      RETURN
!
! Formats
!
!/
!/ End of WMRSPC ----------------------------------------------------- /
!/
      END SUBROUTINE WMRSPC
!/
!/ End of module WMGRIDMD -------------------------------------------- /
!/
      END MODULE WMGRIDMD
