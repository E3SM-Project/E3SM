#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      PROGRAM W3GSPL
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         18-Nov-2013 |
!/                  +-----------------------------------+
!/
!/    24-Sep-2012 : Origination.                        ( version 4.10 )
!/    16-Jan-2013 : Add output of mask file (no halo).  ( version 4.10 )
!/    19-Jan-2013 : Tweaking the template file.         ( version 4.10 )
!/    24-Jan-2013 : Set up for minimum of 2 grids.      ( version 4.10 )
!/                  Add XOFF to grid origin in X.
!/                  Fix IDCLSE for partial grids.
!/                  Add FRFLAG option to disable side-by-side
!/                      running of grids in ww3_multi.
!/    29-Jan-2013 : Add error code on stop.             ( version 4.10 )
!/    31-Jan-2013 : Add routine GRLOST.                 ( version 4.10 )
!/    01-Feb-2013 : Speed up GRSEPA.                    ( version 4.10 )
!/                  Add dynamic trim range in GRTRIM.
!/                  Speed up GRFILL.
!/                  Add small grid merge (GRFSML) early in loop.
!/    04-Feb-2013 : Testing on zero grid size added.    ( version 4.10 )
!/                  Corner point in halo for GR1GRD.
!/    04-Mar-2013 : Adding GrADS output.                ( version 4.10 )
!/    05-Aug-2013 : Add UQ/UNO for distances.           ( version 4.12 )
!/    18-Nov-2013 : Add user-defined halo extension.    ( version 4.14 )
!/
!/    Copyright 2012-2013 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Take an existing grid and create from this the grid data for a set
!     of overlapping grids to be used in the ww3_multi code for hybid
!     paralellization.
!
!  2. Method :
!
!     See Section 8.
!
!  3. Parameters :
!
!     Local parameters.
!     ----------------------------------------------------------------
!       NDSI    Int.  Input unit number ("ww3_prep.inp").
!       NDSO    Int.  Output unit number.
!       NDSE    Int.  Error unit number.
!       NDST    Int.  Test output unit number.
!       NDSM    Int.  Unit number for mod_def file.
!       NG      Int.  Number of grids to be generated.
!       NITMAX  Int.  Maximum number of iterations on grid ref.
!       STARG   Real  std target in percent.
!       GLOBAL  Log.  Closure flag.
!       SEA     L.A.  Sea point map.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3NMOD    Subr. W3GDATMD Set number of model.
!      W3SETG    Subr.   Id.    Point to selected model.
!      W3NDAT    Subr. W3WDATMD Set number of model for wave data.
!      W3SETW    Subr.   Id.    Point to selected model for wave data.
!      W3NOUT    Subr. W3ODATMD Set number of model for output.
!      W3SETO    Subr.   Id.    Point to selected model for output.
!      ITRACE    Subr. W3SERVMD Subroutine tracing initialization.
!      STRACE    Subr.   Id.    Subroutine tracing.
!      NEXTLN    Subr.   Id.    Get next line from input filw
!      EXTCDE    Subr.   Id.    Abort program as graceful as possible.
!      W3IOGR    Subr. W3IOGRMD Reading/writing model definition file.
!
!      GRINFO    Subr. Internal Compile info on all grids.
!      GRTRIM    Subr. Internal Trim edges of grids.
!      GRFILL    Subr. Internal Fill unassigned space in grid.
!      GRLOST    Subr. Internal Assign "lost points".
!      GRSQRG    Subr. Internal Attempt to square-up grid.
!      GRSNGL    Subr. Internal Remove grid points that stick out.
!      GRSEPA    Subr. Internal Remove separated grid pieces.
!      GRFSML    Subr. Internal Deal with fixed minimum size.
!      GRFLRG    Subr. Internal Deal with fixed maximum size.
!      GR1GRD    Subr. Internal Extract single grid from map.
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
!     ----------------------------------------------------
!        1.a  Number of models.
!                   ( W3NMOD , W3NOUT , W3SETG , W3SETO )
!          b  I-O setup.
!          c  Print heading(s).
!        2.   Read model definition file.      ( W3IOGR )
!        3.   Read options from file.
!        4.  Generate first-guess map of sub-grids
!          a Set up array
!          b First cut with regular grid set up
!            1 Set up 'checkerboard'
!            2 Fill checkerboard
!            3 Remove smallest grids as necessary
!            4 Store first guess in MSPLIT
!        5.  Refine map of sub-grids (no halo).
!          a Set up loop.                      ( GRINFO )
!          b Remove small grids.               ( GRFSML )
!          c Trim edges of grids               ( GRTRIM )
!                                      ( GRFILL, GRLOST )
!          d Attempt to square-up grid         ( GRINFO )
!                                              ( GRSQRG )
!                                      ( GRFILL, GRLOST )
!          e Remove mid-sea points sticking out of grid
!                                              ( GRSNGL )
!          f Remove detached grid parts.       ( GRSEPA )
!          g Recompute stats                   ( GRINFO )
!          h Optional GrADS output.
!          i Test convergence
!            Check if stuck on min or max.     ( GRFSML )
!                                              ( GRFLRG )
!          j Test output
!        6.  Output info for all sub grids.
!          a Set up loop.                      ( GRINFO )
!          b Extract grid including halo.      ( GR1GRD )
!        7.  End of program.
!     ----------------------------------------------------
!
!  9. Switches :
!
!     !/PRn   Select propgation scheme.
!
!     !/O16   Generate GrADS output of grid partitioning.
!
!     !/S     Enable subroutine tracing.
!     !/T     Enable test output (main).
!     !/T1    Enable test output (GRINFO).
!     !/T2    Enable test output (GRFILL).
!     !/T3    Enable test output (GRSNGL).
!     !/T4    Enable test output (GRSEPA).
!     !/T5    Enable test output (GRFSML).
!     !/T6    Enable test output (GRFLRG).
!     !/T7    Enable test output (GR1GRD).
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS
!/
!     USE W3GDATMD, ONLY: W3NMOD, W3SETG
      USE W3ADATMD, ONLY: W3NAUX, W3SETA
      USE W3ODATMD, ONLY: W3NOUT, W3SETO
      USE W3SERVMD, ONLY : ITRACE, NEXTLN, EXTCDE
      USE W3ARRYMD, ONLY : OUTA2I, OUTA2R
      USE W3IOGRMD, ONLY: W3IOGR
!/
      USE W3GDATMD
      USE W3ODATMD, ONLY: NDSE, NDST, NDSO, FNMPRE
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: NDSI, NDSM, NDSTRC, NTRACE, J, IERR, &
                                 NG, IX, IY, NGB, NGX, NGY, IG, IGG,  &
                                 IGX, IGY, IGY0, IGYN, IGX0, IGXN,    &
                                 MINGRD, MINNR, MINNXT, MINNNR,       &
                                 NITMAX, IIT, INGMIN, INGMAX,         &
                                 INGMNC, INGMXC, INGLAG, JJ,          &
                                 NSTDLG, MSTDLG = 5, NSEAT, J1, J2,   &
                                 J3, J4, J5, IDFM1, IDFM2, IDFM3,     &
                                 IDLA1, IDLA2, IDLA3, VSC3, NHEXT
      INTEGER, ALLOCATABLE    :: MSPLIT(:,:), MTEMP(:,:), INGRD(:)
      REAL                    :: RATIO1, XMEAN, STARG, STDMIN,        &
                                 ZBDUM, ZBMIN, VSC1, VSC2, FRACL, FRACH
      LOGICAL                 :: GLOBAL, OK, DONE, FRFLAG
      LOGICAL, ALLOCATABLE    :: ISNEXT(:), SEA(:,:)
      CHARACTER(LEN=1)        :: COMSTR
      CHARACTER(LEN=3)        :: G0ID
      CHARACTER(LEN=4)        :: IDGRID, IDCLSE, PTCLSE
      CHARACTER(LEN=6)        :: NRFMT
      CHARACTER(LEN=11)       :: FEXT, AEXT
      CHARACTER(LEN=16)       :: RFORM1, RFORM2, RFORM3
      CHARACTER(LEN=20)       :: FNAME, INAME
!
      TYPE STATS_GRID
        LOGICAL               :: STRADLE, INSTAT
        INTEGER               :: NPTS, NYL, NYH, NXL, NXH
      END TYPE STATS_GRID
!
      TYPE STATS_MEAN
        INTEGER               :: NMIN, NMAX
        REAL                  :: RSTD
      END TYPE STATS_MEAN
!
      TYPE PART_GRID
        INTEGER               :: NX, NY, NSEA
        INTEGER, POINTER      :: MASK(:,:)
        REAL                  :: X0, Y0, SX, SY
        REAL, POINTER         :: ZBIN(:,:), OBSX(:,:), OBSY(:,:)
        LOGICAL               :: GLOBAL
      END TYPE PART_GRID
!
      TYPE(STATS_GRID), POINTER :: GSTATS(:), GSTOLD(:)
      TYPE(STATS_MEAN)          :: MSTATS   , MSTOLD
      TYPE(PART_GRID), POINTER  :: PGRID(:)
!/
!/ ------------------------------------------------------------------- /
!/
! 1.a  Set number of models
!
      CALL W3NMOD ( 1, 6, 6 )
      CALL W3SETG ( 1, 6, 6 )
      CALL W3NAUX (    6, 6 )
      CALL W3SETA ( 1, 6, 6 )
      CALL W3NOUT (    6, 6 )
      CALL W3SETO ( 1, 6, 6 )
!
! 1.b  IO set-up.
!
      NDSI   = 10
      NDSO   =  6
      NDSE   =  6
      NDST   =  6
      NDSM   = 11
!
      NDSTRC =  6
      NTRACE = 100
      CALL ITRACE ( NDSTRC, NTRACE )
!
! 1.c Print header
!
      WRITE (NDSO,900)
!
      J      = LEN_TRIM(FNMPRE)
      OPEN (NDSI,FILE=FNMPRE(:J)//'ww3_gspl.inp',STATUS='OLD',        &
            ERR=800,IOSTAT=IERR)
      REWIND (NDSI)
      READ (NDSI,'(A)',END=801,ERR=802,IOSTAT=IERR) COMSTR
      IF (COMSTR.EQ.' ') COMSTR = '$'
      WRITE (NDSO,901) COMSTR
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 2.  Read model definition file.
!
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
      READ (NDSI,*,END=801,ERR=802,IOSTAT=IERR) FEXT
!
      CALL W3IOGR ( 'READ', NDSM, 1, FEXT )
      CLOSE (NDSM)
!
      WRITE (NDSO,902) FEXT, GNAME
!
      SELECT CASE (GTYPE)
        CASE (RLGTYPE)
          WRITE ( NDSO,903) 'rectilinear'
          IDGRID = 'RECT'
        CASE (CLGTYPE)
          WRITE ( NDSO,903) 'curvictilinear'
          IDGRID = 'CURV'
        CASE (UNGTYPE)
          WRITE ( NDSO,903) 'unstructured'
          IDGRID = 'UNST'
          GOTO 820
        CASE DEFAULT
          WRITE ( NDSO,903) 'not recognized'
          GOTO 821
      END SELECT
!
      SELECT CASE (ICLOSE)
        CASE (ICLOSE_NONE)
          WRITE ( NDSO,904) 'none'
          IDCLSE = 'NONE'
          GLOBAL = .FALSE.
        CASE (ICLOSE_SMPL)
          WRITE ( NDSO,904) 'global (simple)'
          IDCLSE = 'SMPL'
          GLOBAL = .TRUE.
        CASE (ICLOSE_TRPL)
          WRITE ( NDSO,904) 'global (tripolar)'
          IDCLSE = 'TRPL'
          GLOBAL = .TRUE.
          GOTO 822
        CASE DEFAULT
          WRITE ( NDSO,904) 'not recognized'
          GOTO 823
      END SELECT
!
      WRITE (NDSO,905) NX, NY, NSEA
      IF ( NSEA .EQ. 0 ) GOTO 824
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 3.  Read options from input file.
!
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
      READ (NDSI,*,END=801,ERR=802,IOSTAT=IERR) NG, NITMAX, STARG, NHEXT
      NG     = MAX ( 2, NG )
      NITMAX = MAX ( 1, NITMAX )
      STARG  = MAX ( 0. , STARG )
      NHEXT  = MAX ( 0, NHEXT )
      WRITE (NDSO,930) NG, NITMAX, STARG, NHEXT
!
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
      READ (NDSI,*,END=801,ERR=802,IOSTAT=IERR) IDLA1, IDFM1,         &
                                                VSC1, RFORM1
      IF (IDLA1.LT.1 .OR. IDLA1.GT.4) IDLA1  = 1
      IF (IDFM1.LT.1 .OR. IDFM1.GT.3) IDFM1  = 1
      IF ( ABS(VSC1) .LT. 1.E-15 )    VSC1   = 1.
 
      WRITE (NDSO,931) IDLA1, IDFM1, VSC1, RFORM1
!
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
      READ (NDSI,*,END=801,ERR=802,IOSTAT=IERR) IDLA2, IDFM2,         &
                                                VSC2, RFORM2
      IF (IDLA2.LT.1 .OR. IDLA2.GT.4) IDLA2  = 1
      IF (IDFM2.LT.1 .OR. IDFM2.GT.3) IDFM2  = 1
      IF ( ABS(VSC2) .LT. 1.E-15 )    VSC2   = 1.
      IF ( TRFLAG .EQ. 0 ) THEN
          WRITE (NDSO,932)
        ELSE
          WRITE (NDSO,933) IDLA2, IDFM2, VSC2, RFORM2
        END IF
!
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
      READ (NDSI,*,END=801,ERR=802,IOSTAT=IERR) IDLA3, IDFM3,         &
                                                VSC3, RFORM3
      IF (IDLA3.LT.1 .OR. IDLA3.GT.4) IDLA3  = 1
      IF (IDFM3.LT.1 .OR. IDFM3.GT.3) IDFM3  = 1
      IF (     VSC3  .EQ. 0      )    VSC3   = 1
      WRITE (NDSO,934) IDLA3, IDFM3, VSC3, RFORM3
!
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
      READ (NDSI,*,END=801,ERR=802,IOSTAT=IERR) FRACL, FRACH, FRFLAG
      FRACL = MAX ( 0. , FRACL )
      FRACH = MIN ( 1. , FRACH )
      WRITE (NDSO,935) FRACL, FRACH
      IF ( FRACL .GT. FRACH ) GOTO 830
      IF ( .NOT. FRFLAG ) WRITE (NDSO,936)
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 4.  Generate map of sub-grids (no halo)
! 4.a Set up array
!
      ALLOCATE ( MSPLIT(NY,NX) , MTEMP(NY,NX), SEA(NY,NX) )
!
      DO IY=1, NY
        DO IX=1, NX
          IF ( MAPSTA(IY,IX) .EQ. 0 ) THEN
              MSPLIT(IY,IX) = 0
              SEA   (IY,IX) = .FALSE.
            ELSE
              MSPLIT(IY,IX) = -1
              SEA   (IY,IX) = .TRUE.
            END IF
          END DO
        END DO
 
!
! 4.b First cut with regular grid set up
! 4.b.1 Set up 'checkerboard'
!
      RATIO1 = REAL(NX) / REAL(NY)
!
      NGX    = 1
      NGY    = 1
!
      DO
        IF ( NGX*NGY .GE. NG ) EXIT
        IF ( REAL(NGX)/REAL(NGY) .GT. RATIO1 ) THEN
            NGY    = NGY + 1
          ELSE
            NGX    = NGX + 1
          END IF
        END DO
!
      IF ( NGX .GT. NGY ) THEN
          IF ( (NGY-1)*NGX .GE. NG ) NGY = NGY - 1
          IF ( (NGX-1)*NGY .GE. NG ) NGX = NGX - 1
        ELSE
          IF ( (NGY-1)*NGX .GE. NG ) NGY = NGY - 1
          IF ( (NGX-1)*NGY .GE. NG ) NGX = NGX - 1
        END IF
!
! 4.b.2 Fill checkerboard
!
      J      = 0
      DO
!
        MTEMP  = MSPLIT
        IG     = 1
        IGYN   = 0
        J      = J + 1
        ALLOCATE ( INGRD(NGX*NGY) )
        INGRD  = 0
!
        DO IGY=1, NGY
!
          IGY0   = IGYN + 1
          IF ( IGY .EQ. NGY ) THEN
              IGYN   = NY
            ELSE
              IGYN   = NINT ( REAL(NY*IGY) / REAL(NGY)  )
            END IF
          IGXN   = 0
!
          DO IGX=1, NGX
!
            IGX0   = IGXN + 1
            IF ( IGX .EQ. NGX ) THEN
                IGXN   = NX
              ELSE
                IGXN   = NINT ( REAL(NX*IGX) / REAL(NGX)  )
              END IF
!
            DO IX=IGX0, IGXN
              DO IY=IGY0, IGYN
                IF ( MTEMP(IY,IX) .EQ. -1 ) THEN
                    MTEMP(IY,IX) = IG
                    INGRD(IG) = INGRD(IG) + 1
                  END IF
                END DO
              END DO
!
            IF ( INGRD(IG) .GT. 0 ) THEN
                IG     = IG + 1
              END IF
!
            END DO
!
          END DO
!
        IG     = IG - 1
        IF ( IG .LT. NG ) THEN
            IF ( NGX .LT. NGY ) THEN
                NGY     = NGY + 1
              ELSE
                NGX     = NGX + 1
              END IF
            DEALLOCATE ( INGRD )
          ELSE
            EXIT
          END IF
!
        END DO
!
      MINGRD = 0
      DO J=1, IG
        MINGRD = MINGRD + INGRD(J)
        END DO
      IF ( MINGRD .NE. NSEA ) GOTO 825
!
! 4.b.3 Merge smallest grids as necessary
!
      IGG    = IG
!
      DO
!
        IF ( IGG .EQ. NG ) EXIT
!
        MINGRD = NSEA
        MINNR  = 0
        DO J=1, IG
          IF ( INGRD(J) .LT. MINGRD ) THEN
              MINGRD = INGRD(J)
              MINNR  = J
            END IF
          END DO
        INGRD(MINNR) = NSEA + 1
!
        ALLOCATE ( ISNEXT(0:IG) )
        ISNEXT = .FALSE.
!
        DO IY=1, NY-1
          DO IX=1, NX-1
            IF ( ( MTEMP(IY  ,IX  ) - MINNR ) *                       &
                 ( MTEMP(IY+1,IX  ) - MINNR ) *                       &
                 ( MTEMP(IY  ,IX+1) - MINNR ) *                       &
                 ( MTEMP(IY+1,IX+1) - MINNR ) .EQ. 0 ) THEN
                ISNEXT(MTEMP(IY  ,IX  )) = .TRUE.
                ISNEXT(MTEMP(IY+1,IX  )) = .TRUE.
                ISNEXT(MTEMP(IY  ,IX+1)) = .TRUE.
                ISNEXT(MTEMP(IY+1,IX+1)) = .TRUE.
              END IF
            END DO
          END DO
!
        IF ( GLOBAL ) THEN
            DO IY=1, NY-1
              IF ( ( MTEMP(IY  ,NX) - MINNR ) *                       &
                   ( MTEMP(IY+1,NX) - MINNR ) *                       &
                   ( MTEMP(IY  , 1) - MINNR ) *                       &
                   ( MTEMP(IY+1, 1) - MINNR ) .EQ. 0 ) THEN
                  ISNEXT(MTEMP(IY  ,NX)) = .TRUE.
                  ISNEXT(MTEMP(IY+1,NX)) = .TRUE.
                  ISNEXT(MTEMP(IY  , 1)) = .TRUE.
                  ISNEXT(MTEMP(IY+1, 1)) = .TRUE.
                END IF
              END DO
          END IF
!
        MINNXT = NSEA
        MINNNR = 0
        DO J=1, IG
          IF ( ISNEXT(J) .AND. ( INGRD(J) .LT. MINNXT ) ) THEN
              MINNXT = INGRD(J)
              MINNNR = J
            END IF
          END DO
!
        IF ( MINNNR .GT. 0 ) THEN
            DO IY=1, NY
              DO IX=1, NX
                IF ( MTEMP(IY,IX) .EQ. MINNR ) THEN
                    MTEMP(IY,IX) = MINNNR
                    INGRD(MINNNR) = INGRD(MINNNR) + 1
                  END IF
                END DO
              END DO
            IGG    = IGG - 1
          END IF
!
        DEALLOCATE ( ISNEXT)
!
        END DO
!
      DO J=1, IG
        IF ( INGRD(J) .GT. NSEA ) INGRD(J) = 0
        END DO
!
! 4.b.4 Store first guess in MSPLT
!
      IGG    = 0
      DO J=1, IG
        IF ( INGRD(J) .NE. 0 ) THEN
            IGG    = IGG + 1
            DO IY=1, NY
              DO IX=1, NX
                IF ( MTEMP(IY,IX) .EQ. J ) MSPLIT(IY,IX) = IGG
                END DO
              END DO
          END IF
        END DO
!
! 5.b.5 Optional GrADS output
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 5.  Refine grids
! 5.a Set up loop
!
      ALLOCATE ( GSTATS(NG), GSTOLD(NG), PGRID(NG) )
      GSTATS(:)%INSTAT = .TRUE.
      WRITE (NDSO,950)
      DONE   = .FALSE.
!
      CALL GRINFO
      WRITE (NDSO,951) 0, MSTATS%NMIN, MSTATS%NMAX,                   &
                          100.*MSTATS%RSTD/XMEAN
      G0ID   = '5.a'
      IF ( MSTATS%NMIN .EQ. 0 ) GOTO 850
      INGMIN = MSTATS%NMIN
      INGMAX = MSTATS%NMAX
      INGMNC = 0
      INGMXC = 0
      INGLAG = 3
      STDMIN = 100.*MSTATS%RSTD/XMEAN
      NSTDLG = 0
!
      DO IIT=1, NITMAX
!
        IF ( NG .EQ. 1 ) EXIT
!
        MSTOLD = MSTATS
        GSTOLD = GSTATS
!
! 5.b Small grid attempt to merge
!
        IF ( MSTATS%NMIN .LT. NINT(0.45*XMEAN) ) THEN
!
            CALL GRFSML
            CALL GRINFO
!
            G0ID   = '5.b'
            IF ( MSTOLD%NMIN .NE. MSTATS%NMIN ) THEN
                WRITE (NDSO,951) IIT, MSTATS%NMIN, MSTATS%NMAX,       &
                              100.*MSTATS%RSTD/XMEAN
                IF ( MSTATS%NMIN .EQ. 0 ) GOTO 850
 
                CYCLE
              ELSE
                WRITE (NDSO,952)      MSTATS%NMIN, MSTATS%NMAX,       &
                              100.*MSTATS%RSTD/XMEAN
                IF ( MSTATS%NMIN .EQ. 0 ) GOTO 850
              END IF
!
         END IF
!
! 5.c Trim edges of grids and reassign
!
        CALL GRTRIM
        CALL GRFILL ( 2 )
!
! 5.d Attempt to quare-up grid
!
        CALL GRINFO   ! call needed as GRSQRG uses grid ranges
        CALL GRSQRG
        CALL GRFILL ( 1 )
!
! 5.e Remove mid-sea points sticking out of grid
!     Call more than once to remove most .....
!
        OK     = .TRUE.
!
        DO JJ=1, 4
          CALL GRSNGL ( OK )
          END DO
!
! 5.f Remove parts of grid separated from main body, and attachable to
!     other grids.
!
        CALL GRSEPA ( OK , 0.10 )
        IF ( .NOT. OK ) THEN
            CALL GRFILL ( 1 )
            OK     = .TRUE.
          END IF
!
! 5.g Re-compute grid stats
!
        CALL GRINFO
        WRITE (NDSO,951) IIT, MSTATS%NMIN, MSTATS%NMAX,               &
                              100.*MSTATS%RSTD/XMEAN
!
        G0ID   = '5.g'
        IF ( MSTATS%NMIN .EQ. 0 ) GOTO 850
!
! 5.h Optional GrADS output
!
! 5.i Convergence tests
! ... The quick one
!
        IF ( 100.*MSTATS%RSTD/XMEAN .LE. STARG ) THEN
            WRITE (NDSO,959)
            EXIT
          END IF
!
! ... Monitoring convergence ....
!
        IF ( 100.*MSTATS%RSTD/XMEAN .LT. 1.0001*STDMIN ) THEN
            IF ( NSTDLG .LT. MSTDLG ) THEN
                NSTDLG = 0
              ELSE
                WRITE (NDSO,959)
                EXIT
              END IF
            STDMIN = 100.*MSTATS%RSTD/XMEAN
          ELSE
            NSTDLG = NSTDLG + 1
            IF ( NSTDLG .GT. MSTDLG ) STDMIN = 1.01*STDMIN
          END IF
!
! ... Check if stuck on min or max
!
        IF ( MSTATS%NMAX .LT. INGMAX ) THEN
            INGMAX = MSTATS%NMAX
            INGMXC = 0
          ELSE
            INGMXC = INGMXC + 1
          END IF
!
        IF ( MSTATS%NMIN .GT. INGMIN ) THEN
            INGMIN = MSTATS%NMIN
            INGMNC = 0
          ELSE
            INGMNC = INGMNC + 1
          END IF
!
! ... Stuck in min ...
!
        IF ( INGMNC .GE. INGLAG ) THEN
!
            IF ( REAL(INGMIN) .LT. 0.85*XMEAN ) THEN
!
                CALL GRFSML
                CALL GRINFO
                WRITE (NDSO,952) MSTATS%NMIN, MSTATS%NMAX,            &
                                 100.*MSTATS%RSTD/XMEAN
                INGMIN = MSTATS%NMIN
                INGMAX = MSTATS%NMAX
                INGMNC = 0
                INGMXC = 0
                IF ( DONE ) EXIT
!
              END IF
!
          END IF
!
! ... Stuck in max ...
!
        IF ( INGMXC .GE. INGLAG ) THEN
!
            IF ( REAL(INGMAX) .GT. 1.075*XMEAN ) THEN
!
                CALL GRINFO
                WRITE (NDSO,952) MSTATS%NMIN, MSTATS%NMAX,            &
                                 100.*MSTATS%RSTD/XMEAN
                INGMIN = MSTATS%NMIN
                INGMAX = MSTATS%NMAX
                INGMNC = 0
                INGMXC = 0
                IF ( DONE ) EXIT
!
              END IF
!
          END IF
!
        END DO
!
! 5.j Test output
!
      WRITE (NDSO,955)
      ALLOCATE ( ISNEXT(NG) )
      ISNEXT = .TRUE.
      CALL GRINFO
!
      DO JJ=1, NG
        MINNR  = NSEA + 1
        DO J=1, NG
          IF ( ISNEXT(J) .AND. GSTATS(J)%NPTS.LT.MINNR ) THEN
              MINNR  = GSTATS(J)%NPTS
              IG     = J
            END IF
          END DO
        ISNEXT(IG) = .FALSE.
        WRITE (NDST,956) IG, GSTATS(IG)%STRADLE, GSTATS(IG)%NPTS,     &
                             GSTATS(IG)%NXL, GSTATS(IG)%NXH,          &
                             GSTATS(IG)%NYL, GSTATS(IG)%NYH
        END DO
!
      DEALLOCATE ( ISNEXT )
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 6.  Generate output to make separate grids
! 6.a Set up loop
!
      WRITE (NDSO,960)
!
      ZBDUM  = 999.
      IF ( MAXVAL(ZB) .LT. -0.11 ) THEN
          ZBMIN  = -0.1
        ELSE
          ZBMIN  = MAXVAL(ZB) + 1.
          ZBDUM  = MAX ( ZBDUM , ZBMIN+1 )
        END IF
!
      J1     = LEN_TRIM(FEXT)
      J2     = 1 + INT(LOG10(REAL(NG)+0.5))
      WRITE (NRFMT,'(A2,I1,A1,I1,A1)') '(I', J2, '.', J2, ')'
!
      IF ( J1 + J2 + 2 .LE. 10 ) THEN
          FNAME = FEXT(:J1) // '_p'
          J3     = J1 + 3
        ELSE
          FNAME = 'part_'
          J3     = 6
        END IF
      J4     = J3 + J2 - 1
!
      NSEAT  = 0
!
      DO IG=1, NG
!
! 6.b Extract grid including halo
!
        WRITE (NDSO,961) IG
        CALL GR1GRD
        NSEAT  = NSEAT + PGRID(IG)%NSEA
!
        WRITE (AEXT,NRFMT) IG
        FNAME(J3:J4) = AEXT(:J2)
        J      = LEN_TRIM(FNMPRE)
!
! 6.c Writing bottom file
!
        J5    = J4 + 4
        FNAME(J4+1:J5) = '.bot'
        WRITE (NDSO,962) FNAME(:J5)
!
        IF ( IDFM1 .EQ. 3 ) THEN
            OPEN (NDSM,FILE=FNMPRE(:J)//FNAME(:J5),                   &
                  FORM='UNFORMATTED',ERR=860,IOSTAT=IERR)
          ELSE
            OPEN (NDSM,FILE=FNMPRE(:J)//FNAME(:J5), ERR=860,IOSTAT=IERR)
          END IF
        REWIND (NDSM)
        CALL OUTA2R ( PGRID(IG)%ZBIN, PGRID(IG)%NX, PGRID(IG)%NY,     &
                      1, PGRID(IG)%NX, 1, PGRID(IG)%NY, NDSM, NDST,   &
                      NDSE, IDFM1, RFORM1, IDLA1, VSC1, 0.0 )
        CLOSE (NDSM)
!
! 6.d Writing obstruction file
!
        J5    = J4 + 5
        FNAME(J4+1:J5) = '.obst'
!
        IF ( TRFLAG .EQ. 0 ) THEN
            WRITE (NDSO,963) FNAME(:J5)
          ELSE
            WRITE (NDSO,962) FNAME(:J5)
!
            IF ( IDFM2 .EQ. 3 ) THEN
                OPEN (NDSM,FILE=FNMPRE(:J)//FNAME(:J5),               &
                      FORM='UNFORMATTED',ERR=860,IOSTAT=IERR)
              ELSE
                OPEN (NDSM,FILE=FNMPRE(:J)//FNAME(:J5),               &
                      ERR=860,IOSTAT=IERR)
              END IF
            REWIND (NDSM)
            CALL OUTA2R ( PGRID(IG)%OBSX, PGRID(IG)%NX, PGRID(IG)%NY, &
                          1, PGRID(IG)%NX, 1, PGRID(IG)%NY, NDSM,     &
                          NDST, NDSE, IDFM2, RFORM2, IDLA2, VSC2, 0.0 )
            CALL OUTA2R ( PGRID(IG)%OBSY, PGRID(IG)%NX, PGRID(IG)%NY, &
                          1, PGRID(IG)%NX, 1, PGRID(IG)%NY, NDSM,     &
                          NDST, NDSE, IDFM2, RFORM2, IDLA2, VSC2, 0.0 )
            CLOSE (NDSM)
!
          END IF
!
! 6.e Writing mask file
!
        J5    = J4 + 5
        FNAME(J4+1:J5) = '.mask'
        WRITE (NDSO,962) FNAME(:J5)
!
        IF ( IDFM3 .EQ. 3 ) THEN
            OPEN (NDSM,FILE=FNMPRE(:J)//FNAME(:J5),                   &
                  FORM='UNFORMATTED',ERR=860,IOSTAT=IERR)
          ELSE
            OPEN (NDSM,FILE=FNMPRE(:J)//FNAME(:J5), ERR=860,IOSTAT=IERR)
          END IF
        REWIND (NDSM)
        CALL OUTA2I ( PGRID(IG)%MASK, PGRID(IG)%NX, PGRID(IG)%NY,     &
                      1, PGRID(IG)%NX, 1, PGRID(IG)%NY, NDSM, NDST,   &
                      NDSE, IDFM3, RFORM3, IDLA3, VSC3, 0 )
        CLOSE (NDSM)
!
! 6.f Writing input file
!
        J5    = J4 + 5
        FNAME(J4+1:J5) = '.tmpl'
        WRITE (NDSO,962) FNAME(:J5)
!
        OPEN (NDSM,FILE=FNMPRE(:J)//FNAME(:J5), ERR=860,IOSTAT=IERR)
!
        GNAME(31-J2:30) = AEXT
        GNAME(30-J2:30-J2) = 'p'
        WRITE (NDSM,965) GNAME, SIG(2)/SIG(1), TPIINV*SIG(1), NK,     &
                         NTH, TH(1)/DTH, FLDRY, FLCX, FLCY, FLCTH,    &
                         FLCK, FLSOU, DTMAX, DTCFL, DTCFLI, DTMIN
        J5     = LEN_TRIM(RFORM1)
        IF ( REAL(PGRID(IG)%NX) * PGRID(IG)%SX .LT. 359.9 ) THEN
            PTCLSE = 'NONE'
          ELSE
            PTCLSE = IDCLSE
          END IF
        WRITE (NDSM,966) IDGRID, FLAGLL, PTCLSE,                      &
                         PGRID(IG)%NX, PGRID(IG)%NY,                  &
                         PGRID(IG)%SX, PGRID(IG)%SY,                  &
                         PGRID(IG)%X0, PGRID(IG)%Y0,                  &
                         ZBMIN, DMIN, VSC1, IDLA1, IDFM1,             &
                         RFORM1(:J5), FNAME(:J4)//'.bot'
        IF ( TRFLAG .NE. 0 ) THEN
        J5     = LEN_TRIM(RFORM2)
            WRITE (NDSM,967) VSC2,IDLA2, IDFM2, RFORM2(:J5),          &
                             FNAME(:J4)//'.obst'
          END IF
        J5     = LEN_TRIM(RFORM3)
        WRITE (NDSM,968) IDLA3, IDFM3, RFORM3(:J5), FNAME(:J4)//'.mask'
        CLOSE (NDSM)
!
        END DO
!
      WRITE (NDSO,969) 100. * (REAL(NSEAT)/REAL(NSEA)-1.)
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 7.  Write part of ww3_multi.inp
!
      J5     = 11+J1+J2
      INAME(:J5) = 'ww3_multi.'//FEXT(:J1)//'.'//AEXT(:J2)
      OPEN (NDSM,FILE=FNMPRE(:J)//INAME(:J5), ERR=870,IOSTAT=IERR)
!
      DO IG=1, NG
        WRITE (AEXT,NRFMT) IG
        FNAME(J3:J4) = AEXT(:J2)
        IF ( FRFLAG ) THEN
            WRITE (NDSM,970) FNAME(:J4),                              &
                         FRACL + REAL(IG-1)*(FRACH-FRACL)/REAL(NG),   &
                         FRACL + REAL( IG )*(FRACH-FRACL)/REAL(NG)
          ELSE
            WRITE (NDSM,970) FNAME(:J4), FRACL, FRACH
          END IF
        END DO
!
      CLOSE (NDSM)
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 8.  Write mask file (no halo)
!
      J5     = 10+J1+J2
      INAME(:J5) = 'ww3_mask.'//FEXT(:J1)//'.'//AEXT(:J2)
      OPEN (NDSM,FILE=FNMPRE(:J)//INAME(:J5), ERR=870,IOSTAT=IERR)
!
      DO IY=1, NY
        WRITE (NDSM,980) MSPLIT(IY,:)
        END DO
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 9.  End of program
!
      GOTO 888
!
! Error escape locations
!
  800 CONTINUE
      WRITE (NDSE,1000) IERR
      CALL EXTCDE ( 40 )
!
  801 CONTINUE
      WRITE (NDSE,1001)
      CALL EXTCDE ( 41 )
!
  802 CONTINUE
      WRITE (NDSE,1002) IERR
      CALL EXTCDE ( 42 )
!
  820 CONTINUE
      WRITE (NDSE,1020) GTYPE
      CALL EXTCDE ( 20 )
!
  821 CONTINUE
      WRITE (NDSE,1021) GTYPE
      CALL EXTCDE ( 21 )
!
  822 CONTINUE
      WRITE (NDSE,1022) ICLOSE
      CALL EXTCDE ( 22 )
!
  823 CONTINUE
      WRITE (NDSE,1023) ICLOSE
      CALL EXTCDE ( 23 )
!
  824 CONTINUE
      WRITE (NDSE,1024)
      CALL EXTCDE ( 24 )
!
  825 CONTINUE
      WRITE (NDSE,1025) MINGRD, NSEA
      CALL EXTCDE ( 25 )
!
  830 CONTINUE
      WRITE (NDSE,1030)
      CALL EXTCDE ( 30 )
!
  850 CONTINUE
      WRITE (NDSE,1050) G0ID
      CALL EXTCDE ( 50 )
!
  860 CONTINUE
      WRITE (NDSE,1060) FNMPRE(:J)//FNAME(:J5), IERR
      CALL EXTCDE ( 60 )
!
  870 CONTINUE
      WRITE (NDSE,1070) FNMPRE(:J)//INAME(:J5), IERR
      CALL EXTCDE ( 70 )
!
  888 CONTINUE
      WRITE (NDSO,999)
!
! Formats
!
  900 FORMAT (/15X,'  *** WAVEWATCH III  Grid splitting ***  '/ &
               15X,'=========================================='/)
  901 FORMAT ( '  Comment character is ''',A,''''/)
  902 FORMAT ( '  Grid ID   : ',A/                                    &
               '  Grid name : ',A)
  903 FORMAT ( '  Grid type : ',A)
  904 FORMAT ( '  Closure   : ',A)
  905 FORMAT ( '  Grid size : ',I4,' x',I4,'   (',I8,')'/)
!
  930 FORMAT ( '  Generating ',I3,' grids'/                           &
               '  No more than',I4,' refinement iterations'/          &
               '  Grid point count std target (%) :',F6.2/            &
               '  Halo per sub grid extended by',I3,' grid point.')
  931 FORMAT ( '  Format info for bottom file      :',2I2,F12.4,2X,A)
  932 FORMAT ( '  Format info for obstruction file not used')
  933 FORMAT ( '  Format info for obstruction file :',2I2,F12.4,2X,A)
  934 FORMAT ( '  Format info for mask file        :',2I2,I7,7X,A)
  935 FORMAT ( '  Part of cummunicator to be used  :',2F7.4)
  936 FORMAT ( '  Not running grids side-by-side'/                    &
               '     *** NON CONVENTIONAL OPERATION ***'/)
!
  950 FORMAT (/'  Iterations:'/                                       &
               '     nr     min     max    std (%) '/                 &
               '  ---------------------------------')
  951 FORMAT (2X,I5,2I8,2F10.2)
  952 FORMAT (2X,5x,2I8,2F10.2)
  955 FORMAT (/'  Resulting grids:'/                                  &
               '     grid  stradle points  range X   range Y '/       &
               '  ---------------------------------------------')
  956 FORMAT ( '     ',I4,5X,L1,2X,I7,4I5)
  959 FORMAT ( '  Convergence reached')
!
  960 FORMAT (/'  Generating grid data:'/                             &
               '  ---------------------------------------------')
  961 FORMAT ( '     Extracting data for grid',I4)
  962 FORMAT ( '        Writing file ',A)
  963 FORMAT ( '        File ',A,' not requested')
!
  970 FORMAT ( ' ''',A,''' ''LEV'' ''CUR'' ''WND'' ''ICE''',          &
               ' ''D1'' ''D2'' ''D3'' RANK GROUP',2F10.7,' BFLAG')
!
  980 FORMAT (1X,360I2)
!
  965 FORMAT ( '$ -------------------------------------',             &
               '------------------------------- $'/                   &
               '$ WAVEWATCH III Grid preprocessor input',             &
               ' file                           $'/                   &
               '$ -------------------------------------',             &
               '------------------------------- $'/                   &
               '   ''',A,''''/'$'/                                    &
               '   ',F8.4,F10.6,2I6,F8.4/'   ',6L2/'   ',4F12.4/      &
               '$ NAMELISTS'/'$')
  966 FORMAT ( '   ''',A4,'''  ',L1,'  ''',A4,''''/1X,I8,I12/         &
               4X,2F12.6,'    1.0'/4X,2F12.6,'    1.0'/2F8.2,'  20',  &
               F12.6,2I2,' ''',A,''' ''NAME'' ''',A,'''')
  967 FORMAT ( 18X,'30',F12.6,2I2,' ''',A,''' ''NAME'' ''',A,'''' )
  968 FORMAT ( 18X,'40',12X,2I2,' ''',A,''' ''NAME'' ''',A,''''/'$'/  &
               '$ Note: cannot make output boundary points here'/'$'/ &
               ' 0. 0. 0. 0. 0'/                                      &
               '$ -------------------------------------',             &
               '------------------------------- $'/                   &
               '$ End of input file                    ',             &
               '                                $'/                   &
               '$ -------------------------------------',             &
               '------------------------------- $')
!
  969 FORMAT (/'  Grid point inflation',F7.2,'%')
!
  999 FORMAT(//'  End of program '/                                   &
               ' ========================================='/          &
               '             WAVEWATCH III Grid splitting '/)
!
 1000 FORMAT (/' *** WAVEWATCH III ERROR IN W3GSPL : '/               &
               '     ERROR IN OPENING INPUT FILE'/                    &
               '     IOSTAT =',I5/)
!
 1001 FORMAT (/' *** WAVEWATCH III ERROR IN W3GSPL : '/               &
               '     PREMATURE END OF INPUT FILE'/)
!
 1002 FORMAT (/' *** WAVEWATCH III ERROR IN W3GSPL : '/               &
               '     ERROR IN READING FROM INPUT FILE'/               &
               '     IOSTAT =',I5/)
!
 1020 FORMAT (/' *** WAVEWATCH III ERROR IN W3GSPL : '/               &
               '     SPLITTING NOT AVAILABLE FOR GRID TYPE'/          &
               '     GTYPE =',I5/)
!
 1021 FORMAT (/' *** WAVEWATCH III ERROR IN W3GSPL : '/               &
               '     GRID TYPE NOT RECOGNIZED'/                       &
               '     GTYPE =',I5/)
!
 1022 FORMAT (/' *** WAVEWATCH III ERROR IN W3GSPL : '/               &
               '     SPLITTING NOT AVAILABLE FOR CLOSURE TYPE'/       &
               '     ICLOSE =',I5/)
!
 1023 FORMAT (/' *** WAVEWATCH III ERROR IN W3GSPL : '/               &
               '     CLOSURE TYPE NOT RECOGNIZED'/                    &
               '     ICLOSE =',I5/)
!
 1024 FORMAT (/' *** WAVEWATCH III ERROR IN W3GSPL : '/               &
               '     NO ACTIVE SEA POINT IN GRID'/)
!
 1025 FORMAT (/' *** WAVEWATCH III ERROR IN W3GSPL : '/               &
               '     WRONG NUMBER OF SEA POINTS'/                     &
               '     MINGRD, NSEA =',2I7/)
!
 1030 FORMAT (/' *** WAVEWATCH III ERROR IN W3GSPL : '/               &
               '     ILLEGAL PART OF COMMUNICATOR REQUESTED'/)
!
 1050 FORMAT (/' *** WAVEWATCH III ERROR IN W3GSPL : '/               &
               '     SHOULD NOT HAVE ZERO GRID SIZE (',A,') ...'/)
!
 1060 FORMAT (/' *** WAVEWATCH III ERROR IN W3GSPL : '/               &
               '     ERROR IN OPENING FILE ',A/                       &
               '     IOSTAT =',I5/)
!
 1070 FORMAT (/' *** WAVEWATCH III ERROR IN W3GSPL : '/               &
               '     ERROR IN OPENING FILE ',A/                       &
               '     IOSTAT =',I5/)
!
 9052 FORMAT ( 'TEST W3GSPL: STUCK ON ',A,' GRID SIZE')
 9053 FORMAT ( '             OUT OF RANGE, PROCESSING (',F6.3,')')
 9054 FORMAT ( '             IN RANGE, NO ACTION')
!/
!/ Embedded subroutines ---------------------------------------------- /
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE GRINFO
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         13-Sep-2012 |
!/                  +-----------------------------------+
!/
!/    06-Sep-2012 : Origination.                        ( version 4.10 )
!/    13-Sep-2012 : Option to exclude grids from stats. ( version 4.10 )
!/
!  1. Purpose :
!
!     Compile statistical info on all sub grids (no halo).
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: NOCNT, NOCNTM, NOCNTL, NGC, NSEAC
      REAL                    :: SUMSQR
      LOGICAL                 :: LEFT, RIGHT, THERE
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Initialization ------------------------------------------------- *
!
      GSTATS(:)%STRADLE = .FALSE.
      GSTATS(:)%NPTS    = 0
      GSTATS(:)%NXL     = NX
      GSTATS(:)%NXH     = 1
      GSTATS(:)%NYL     = NY
      GSTATS(:)%NYH     = 1
!
! 2.  Get STRADLE, NGC ----------------------------------------------- *
!
      NGC    = 0
!
      DO IG=1, NG
        LEFT   = .FALSE.
        RIGHT  = .FALSE.
        IF ( GSTATS(IG)%INSTAT ) NGC = NGC + 1
        DO IY=1, NY
          IF ( MSPLIT(IY, 1) .EQ. IG ) LEFT   = .TRUE.
          IF ( MSPLIT(IY,NX) .EQ. IG ) RIGHT  = .TRUE.
          END DO
        GSTATS(IG)%STRADLE = LEFT .AND. RIGHT
        END DO
!
      IF ( NGC .EQ. 0 ) THEN
          NGC    = 1
          DONE   = .TRUE.
        END IF
!
! 3.  Run grid stats ------------------------------------------------- *
! 3.a General
!
      DO IY=1, NY
        DO IX=1, NX
          IG     = MSPLIT(IY,IX)
          IF ( MSPLIT(IY,IX) .GT. 0 ) THEN
              GSTATS(IG)%NPTS = GSTATS(IG)%NPTS + 1
              GSTATS(IG)%NXL  = MIN ( GSTATS(IG)%NXL , IX )
              GSTATS(IG)%NXH  = MAX ( GSTATS(IG)%NXH , IX )
              GSTATS(IG)%NYL  = MIN ( GSTATS(IG)%NYL , IY )
              GSTATS(IG)%NYH  = MAX ( GSTATS(IG)%NYH , IY )
            END IF
          END DO
        END DO
!
! 3.b Stradled grids
!
      IF ( NG .GT. 1) THEN
        DO IG=1, NG
          IF ( GSTATS(IG)%STRADLE ) THEN
              NOCNT  = 0
              NOCNTM = 0
              NOCNTL = 0
              DO IX=1, NX
                THERE  = .FALSE.
                DO IY=1, NY
                  IF ( MSPLIT(IY,IX) .EQ. IG ) THEN
                      THERE  = .TRUE.
                      EXIT
                    END IF
                  END DO
                IF ( THERE ) THEN
                    NOCNT = 0
                  ELSE
                    NOCNT = NOCNT + 1
                    IF ( NOCNT .GT. NOCNTM ) THEN
                        NOCNTM = NOCNT
                        NOCNTL = IX
                      END IF
                  END IF
                END DO
              GSTATS(IG)%NXL = NOCNTL + 1
              GSTATS(IG)%NXH = NOCNTL - NOCNTM
            END IF
          END DO
        ELSE
          GSTATS(1)%STRADLE = .FALSE.
        END IF
!
! 3.c Corrected NSEA
!
      NSEAC   = 0
!
      DO IG=1, NG
        IF ( GSTATS(IG)%INSTAT ) NSEAC = NSEAC + GSTATS(IG)%NPTS
        END DO
!
! 4.  Run overall stats ---------------------------------------------- *
!
      MSTATS%NMIN = NSEA + 1
      MSTATS%NMAX = 0
      XMEAN  = REAL(NSEAC) / REAL(NGC)
      SUMSQR = 0.
!
      DO IG=1, NG
        IF ( .NOT. GSTATS(IG)%INSTAT ) CYCLE
        MSTATS%NMIN = MIN ( MSTATS%NMIN , GSTATS(IG)%NPTS )
        MSTATS%NMAX = MAX ( MSTATS%NMAX , GSTATS(IG)%NPTS )
        SUMSQR      = SUMSQR + ( REAL(GSTATS(IG)%NPTS) - XMEAN )**2
        END DO
!
      MSTATS%RSTD = SQRT ( SUMSQR / REAL(NGC) )
!
! 5.  Test output ---------------------------------------------------- *
!
      RETURN
!
! Formats
!
!/ End of GRINFO ----------------------------------------------------- /
!/
      END SUBROUTINE GRINFO
!/ ------------------------------------------------------------------- /
      SUBROUTINE GRTRIM
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Feb-2013 |
!/                  +-----------------------------------+
!/
!/    07-Sep-2012 : Origination.                        ( version 4.10 )
!/    18-Sep-2012 : Include edge points of grid.        ( version 4.10 )
!/    01-Feb-2013 : Add dynamic trim range.             ( version 4.10 )
!/
!  1. Purpose :
!
!    Trim edges of all grids where they are next to another grid or next
!    to unassigned grid points. This is done in preparation for
!    reassigning edges of grids to smaller adjacent grids.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: ITARG, ITL, IPTS, MX, MY, ICIRC, NWDTH
      LOGICAL                 :: MASK(NY,NX)
!/
!/ ------------------------------------------------------------------- /
!/
!
      ITARG  = NSEA / NG
!
! 1.  Loop over grids ------------------------------------------------ *
!
      DO IG=1, NG
!
        IPTS   = GSTATS(IG)%NPTS
        MY     = 1 + GSTATS(IG)%NYH - GSTATS(IG)%NYL
        MX     = 1 + GSTATS(IG)%NXH - GSTATS(IG)%NXL
        IF ( GSTATS(IG)%STRADLE ) MX     = MX + NX
        ICIRC  = 2 * ( MX + MY )
!
        NWDTH  = 1
!
        ITL    = MIN ( ITARG , MAX ( ITARG-2*ICIRC , 3*ICIRC ) )
        IF ( IPTS .LT. ITL ) NWDTH  = 0
!
        IF ( IPTS.GT.ITARG ) THEN
            NWDTH = 1 +                                              &
                   MAX(0,+NINT((REAL((IPTS-ITARG))/REAL(ICIRC)-1.)/3.))
          ENDIF
!
        DO J=1, NWDTH
!
          MASK   = .FALSE.
!
! 2.  Mark points to be removed -------------------------------------- *
!
          DO IX=2, NX-1
            IF ( MSPLIT( 1,IX) .EQ. IG ) MASK( 1,IX) =                &
                     (SEA( 2,IX  ).AND.(MSPLIT( 2,IX  ).NE.IG))       &
                .OR. (SEA( 1,IX+1).AND.(MSPLIT( 1,IX+1).NE.IG))       &
                .OR. (SEA( 1,IX-1).AND.(MSPLIT( 1,IX-1).NE.IG))
            DO IY=2, NY-1
              IF ( MSPLIT(IY,IX) .EQ. IG ) MASK(IY,IX) =              &
                       (SEA(IY+1,IX  ).AND.(MSPLIT(IY+1,IX  ).NE.IG)) &
                  .OR. (SEA(IY-1,IX  ).AND.(MSPLIT(IY-1,IX  ).NE.IG)) &
                  .OR. (SEA(IY  ,IX+1).AND.(MSPLIT(IY  ,IX+1).NE.IG)) &
                  .OR. (SEA(IY  ,IX-1).AND.(MSPLIT(IY  ,IX-1).NE.IG))
              END DO
            IF ( MSPLIT(NY,IX) .EQ. IG ) MASK(NY,IX) =                &
                     (SEA(NY-1,IX  ).AND.(MSPLIT(NY-1,IX  ).NE.IG))   &
                .OR. (SEA(NY  ,IX+1).AND.(MSPLIT(NY  ,IX+1).NE.IG))   &
                .OR. (SEA(NY  ,IX-1).AND.(MSPLIT(NY  ,IX-1).NE.IG))
            END DO
!
          IF ( GLOBAL ) THEN
              IF ( MSPLIT( 1, 1) .EQ. IG ) MASK( 1, 1) =              &
                       (SEA( 2, 1).AND.(MSPLIT( 2, 1).NE.IG))         &
                  .OR. (SEA( 1, 2).AND.(MSPLIT( 1, 2).NE.IG))         &
                  .OR. (SEA( 1,NX).AND.(MSPLIT( 1,NX).NE.IG))
              IF ( MSPLIT( 1,NX) .EQ. IG ) MASK( 1,NX) =              &
                       (SEA( 2,NX  ).AND.(MSPLIT( 2,NX  ).NE.IG))     &
                  .OR. (SEA( 1, 1  ).AND.(MSPLIT( 1, 1  ).NE.IG))     &
                  .OR. (SEA( 1,NX-1).AND.(MSPLIT( 1,NX-1).NE.IG))
              DO IY=2, NY-1
                IF ( MSPLIT(IY, 1) .EQ. IG ) MASK(IY, 1) =            &
                         (SEA(IY+1, 1).AND.(MSPLIT(IY+1, 1).NE.IG))   &
                    .OR. (SEA(IY-1, 1).AND.(MSPLIT(IY-1, 1).NE.IG))   &
                    .OR. (SEA(IY  , 2).AND.(MSPLIT(IY  , 2).NE.IG))   &
                    .OR. (SEA(IY  ,NX).AND.(MSPLIT(IY  ,NX).NE.IG))
                IF ( MSPLIT(IY,NX) .EQ. IG ) MASK(IY,NX) =            &
                         (SEA(IY+1,NX).AND.(MSPLIT(IY+1,NX).NE.IG))   &
                    .OR. (SEA(IY-1,NX).AND.(MSPLIT(IY-1,NX).NE.IG))   &
                    .OR. (SEA(IY  , 1).AND.(MSPLIT(IY  , 1).NE.IG))   &
                    .OR. (SEA(IY,NX-1).AND.(MSPLIT(IY,NX-1).NE.IG))
                END DO
              IF ( MSPLIT(NY, 1) .EQ. IG ) MASK(NY, 1) =              &
                       (SEA(NY-1, 1).AND.(MSPLIT(NY-1, 1).NE.IG))     &
                  .OR. (SEA(NY  , 2).AND.(MSPLIT(NY  , 2).NE.IG))     &
                  .OR. (SEA(NY  ,NX).AND.(MSPLIT(NY  ,NX).NE.IG))
              IF ( MSPLIT(NY,NX) .EQ. IG ) MASK(NY,NX) =              &
                       (SEA(NY-1,NX).AND.(MSPLIT(NY-1,NX).NE.IG))     &
                  .OR. (SEA(NY  , 1).AND.(MSPLIT(NY  , 1).NE.IG))     &
                  .OR. (SEA(NY,NX-1).AND.(MSPLIT(NY,NX-1).NE.IG))
          ELSE
              IF ( MSPLIT( 1, 1) .EQ. IG ) MASK( 1, 1) =              &
                       (SEA( 2, 1).AND.(MSPLIT( 2, 1).NE.IG))         &
                  .OR. (SEA( 1, 2).AND.(MSPLIT( 1, 2).NE.IG))
              IF ( MSPLIT( 1,NX) .EQ. IG ) MASK( 1,NX) =              &
                       (SEA( 2,NX  ).AND.(MSPLIT( 2,NX  ).NE.IG))     &
                  .OR. (SEA( 1,NX-1).AND.(MSPLIT( 1,NX-1).NE.IG))
              DO IY=2, NY-1
                IF ( MSPLIT(IY, 1) .EQ. IG ) MASK(IY, 1) =            &
                         (SEA(IY+1, 1).AND.(MSPLIT(IY+1, 1).NE.IG))   &
                    .OR. (SEA(IY-1, 1).AND.(MSPLIT(IY-1, 1).NE.IG))   &
                    .OR. (SEA(IY  , 2).AND.(MSPLIT(IY  , 2).NE.IG))
                IF ( MSPLIT(IY,NX) .EQ. IG ) MASK(IY,NX) =            &
                         (SEA(IY+1,NX).AND.(MSPLIT(IY+1,NX).NE.IG))   &
                    .OR. (SEA(IY-1,NX).AND.(MSPLIT(IY-1,NX).NE.IG))   &
                    .OR. (SEA(IY,NX-1).AND.(MSPLIT(IY,NX-1).NE.IG))
                END DO
              IF ( MSPLIT(NY, 1) .EQ. IG ) MASK(NY, 1) =              &
                       (SEA(NY-1, 1).AND.(MSPLIT(NY-1, 1).NE.IG))     &
                  .OR. (SEA(NY  , 2).AND.(MSPLIT(NY  , 2).NE.IG))
              IF ( MSPLIT(NY,NX) .EQ. IG ) MASK(NY,NX) =              &
                       (SEA(NY-1,NX).AND.(MSPLIT(NY-1,NX).NE.IG))     &
                  .OR. (SEA(NY,NX-1).AND.(MSPLIT(NY,NX-1).NE.IG))
            END IF
!
! 3.  Remove marked points ------------------------------------------- *
!
          DO IX=1, NX
            DO IY=1, NY
              IF ( MASK(IY,IX) ) THEN
                  MSPLIT(IY,IX) = -1
                END IF
              END DO
            END DO
!
! ... End loops started in 1.
!
          END DO
        END DO
!
      RETURN
!
! Formats
!
!/ End of GRTRIM ----------------------------------------------------- /
!/
      END SUBROUTINE GRTRIM
!/ ------------------------------------------------------------------- /
      SUBROUTINE GRFILL ( ND )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Feb-2013 |
!/                  +-----------------------------------+
!/
!/    07-Sep-2012 : Origination.                        ( version 4.10 )
!/    18-Sep-2012 : Include edge points of grid.        ( version 4.10 )
!/                  Add convergence check.
!/    29-Jan-2013 : Add error code on stop.             ( version 4.10 )
!/    29-Jan-2013 : Add error test output.              ( version 4.10 )
!/    01-Feb-2013 : Loop over selected sea points only. ( version 4.10 )
!/
!  1. Purpose :
!
!     Reassign unassigned grid points to grids, starting with the
!     smallest grids.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       ND      Int.   I   Depth of halo for first sweep.
!     ----------------------------------------------------------------
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: ND
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: NMIN, I, NDEPTH, NITT, NADD, IXL, IXR,&
                                 NLEFT, NRIGHT, NXL, NXH, NYL, NYH
      INTEGER                 :: NXYOFF = 3
      INTEGER                 :: IIX(NSEA), IIY(NSEA), ISEA, NSEAL
      LOGICAL                 :: DONE(NG), MASK(NY,NX), FLOST(NG),     &
                                 XFL(NX), YFL(NY)
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Loop to assure all reassigned ---------------------------------- *
!
      NDEPTH = ND
      NITT   = 0
      NLEFT  = -1
      FLOST  = .FALSE.
!
      NSEAL  = 0
      DO IX=1, NX
        DO IY=1, NY
          IF ( MSPLIT(IY,IX) .EQ. -1 ) THEN
              NSEAL  = NSEAL + 1
              IIX(NSEAL) = IX
              IIY(NSEAL) = IY
            END IF
          END DO
        END DO
!
      DO
        NITT   = NITT + 1
!
! 2.  Loop over all grids -------------------------------------------- *
!
        DONE   = .FALSE.
!
        DO J=1, NG
!
! 3.  Find smallest unprocessed grid --------------------------------- *
!
          NMIN   = NSEA + 1
          IG     = 0
!
          DO I=1, NG
            IF ( .NOT.DONE(I) .AND. GSTATS(I)%NPTS.LT.NMIN ) THEN
                IG     = I
                NMIN   = GSTATS(I)%NPTS
              END IF
            END DO
!
          DONE(IG) = .TRUE.
!
! 4.  Loop for halos per grid ---------------------------------------- *
!
          DO, I=1, NDEPTH
!
            MASK   = .FALSE.
!
! 5.  Mark grid point for adding ------------------------------------- *
!
            DO ISEA=1, NSEAL
              IX     = IIX(ISEA)
              IY     = IIY(ISEA)
              IXL    = 1 + MOD(IX-2+NX,NX)
              IXR    = 1 + MOD(IX,NX)
              IF ( MSPLIT(IY,IX) .EQ. -1 ) MASK(IY,IX) =              &
                                      ( MSPLIT(IY+1,IX ) .EQ. IG )    &
                                .OR.  ( MSPLIT(IY-1,IX ) .EQ. IG )    &
                                .OR.  ( MSPLIT(IY  ,IXR) .EQ. IG )    &
                                .OR.  ( MSPLIT(IY  ,IXL) .EQ. IG )
              END DO
!
! 6.  Add marked grid point ------------------------------------------ *
!
            NADD   = 0
!
            DO ISEA=1, NSEAL
              IX     = IIX(ISEA)
              IY     = IIY(ISEA)
              IF ( MASK(IY,IX) ) THEN
                  MSPLIT(IY,IX) = IG
                  NADD          = NADD + 1
                END IF
              END DO
!
            IF ( NADD .EQ. 0 ) EXIT
!
! ... End loop started in 4.
!
            END DO
!
! ... End loop started in 2.
!
          END DO
!
        NDEPTH = 1
!
! 7.  Check convergence ---------------------------------------------- *
! 7.a Find number of points left
!
        NRIGHT = NLEFT
        NLEFT  = 0
!
        DO ISEA=1, NSEAL
          IX     = IIX(ISEA)
          IY     = IIY(ISEA)
          IF ( MSPLIT(IY,IX) .EQ. -1 ) NLEFT = NLEFT + 1
          END DO
!
! 7.b No point left, exit loop
!
        IF ( NLEFT .EQ. 0 ) EXIT
!
! 7.c Stuck with points left
!
        IF ( NRIGHT .GT. 0 ) THEN
            IF ( NLEFT .EQ. NRIGHT ) THEN
!
! 7.d Do lost point correction once
!
                IF ( .NOT. FLOST(IG) ) THEN
                    CALL GRLOST
                    FLOST(IG) = .TRUE.
                  ELSE
!
! 7.e Got stuck for good, error message and ouput
!
                    WRITE (NDSE,1000) IG, NITT, NLEFT
!
                    XFL    = .FALSE.
                    YFL    = .FALSE.
!
                    DO ISEA=1, NSEAL
                      IX     = IIX(ISEA)
                      IY     = IIY(ISEA)
                      IF ( MSPLIT(IY,IX) .EQ. -1 ) THEN
                          XFL(MAX(1,IX-NXYOFF):MIN(NX,IX+NXYOFF)) = .TRUE.
                          YFL(MAX(1,IY-NXYOFF):MIN(NY,IY+NXYOFF)) = .TRUE.
                          END IF
                      END DO
!
                    NXL   = 0
                    NXH   = 0
                    DO IX=1, NX
                      IF ( XFL(IX) .AND. NXL.EQ. 0 ) NXL = IX
                      IF ( XFL(IX) .AND. IX.EQ. NX ) NXH = IX
                      IF ( .NOT. XFL(IX) .AND. NXL.NE. 0 ) NXH = IX-1
                      IF ( NXH .NE. 0 ) THEN
                          NYL   = 0
                          NYH   = 0
                          DO IY=1, NY
                            IF ( YFL(IY) .AND. NYL.EQ. 0 ) NYL = IY
                            IF ( YFL(IY) .AND. IY.EQ. NY ) NYH = IY
                            IF ( .NOT. YFL(IY) .AND. NYL.NE. 0 )      &
                                                            NYH = IY-1
                            IF ( NYH .NE. 0 ) THEN
                                WRITE (NDST,1001) NXL, NXH, NYH, NYL
                                DO I=NYH, NYL, -1
                              WRITE (NDST,1002) MSPLIT(I,NXL:NXH)
                                  END DO
                                NYL    = 0
                                NYH    = 0
                              END IF
                            END DO
                          NXL   = 0
                          NXH   = 0
                        END IF
                      END DO
!
! ... Stop program with error output ...
!
                    STOP 01
                  ENDIF
!
              END IF
          END IF
!
! ... End loop started in 1.
!
        END DO
!
      RETURN
!
! Formats
!
 1000 FORMAT (/' *** ERROR GRFILL : NO MORE CONVERGENCE, ',           &
               'NITT, NLEFT:',2I8,' ***'/)
 1001 FORMAT ( ' MAP OUTPUT FOR GRID',I3,' AND X AND Y RANGE :',4I6/)
 1002 FORMAT ( ' ',60I2)
!
!/ End of GRFILL ----------------------------------------------------- /
!/
      END SUBROUTINE GRFILL
!/ ------------------------------------------------------------------- /
      SUBROUTINE GRLOST
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         .9-Jan-2013 |
!/                  +-----------------------------------+
!/
!/    31-Jan-2013 : Origination.                        ( version 4.10 )
!/
!  1. Purpose :
!
!     Reassign unassigned grid points to gridsR. Dealing with lost
!     point by finding clostst grids.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IX, IY, IOFF, JJX, JX, JY, IG, I
      INTEGER                 :: IFOUND(-1:NG)
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Loop over all grid points -------------------------------------- *
!
      DO IX=1, NX
        DO IY=1, NY
!
          IF ( MSPLIT(IY,IX) .EQ. -1 ) THEN
!
! 2.  Find nearest grid(s) ------------------------------------------- *
!
              IOFF   = 1
!
              DO
!
                IFOUND = 0
                DO JJX=IX-IOFF, IX+IOFF
                  IF ( GLOBAL ) THEN
                      JX     = 1 + MOD(JJX-1+2*NX,NX)
                    ELSE
                      JX     = JJX
                    END IF
                   IF ( JX.LT.1 .OR. JX.GT.NX ) CYCLE
                   DO JY=IY-IOFF, IY+IOFF
                     IF ( JY.LT.1 .OR. JY.GT.NY ) CYCLE
                     IFOUND(MSPLIT(JY,JX)) = IFOUND(MSPLIT(JY,JX)) + 1
                    END DO
                  END DO
!
                IG     = 0
                DO I=1, NG
                  IF ( IFOUND(I) .GT. 0 ) THEN
                      IG     = I
                      EXIT
                    END IF
                  END DO
!
                IF ( IG .NE. 0 ) THEN
                    MSPLIT(IY,IX) = IG
                    EXIT
                  END IF
!
                IOFF   = IOFF + 1
                IF ( IOFF .GT. NX .AND. IOFF.GT.NY ) EXIT
                END DO
!
! ... End of loops and logic started in 1.
!
            END IF
!
          END DO
        END DO
!
      RETURN
!
! Formats
!
!/ End of GRLOST ----------------------------------------------------- /
!/
      END SUBROUTINE GRLOST
!/ ------------------------------------------------------------------- /
      SUBROUTINE GRSQRG
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         07-Sep-2012 |
!/                  +-----------------------------------+
!/
!/    07-Sep-2012 : Origination.                        ( version 4.10 )
!/
!  1. Purpose :
!
!     Attemp to square-up grid by taking off grid point in outermost
!     grid point in X and Y only, after which GRFILL is to be run to
!     re-assign grid points,
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: MX, MY
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Loop over grids ------------------------------------------------ *
!
      DO IG=1, NG
!
        MY     = 1 + GSTATS(IG)%NYH - GSTATS(IG)%NYL
        MX     = 1 + GSTATS(IG)%NXH - GSTATS(IG)%NXL
        IF ( GSTATS(IG)%STRADLE ) MX     = MX + NX
!
! 2.  Top ------------------------------------------------------------ *
!
        IF ( MY .GE. 5 ) THEN
!
            DO IX=1, NX
              IF (MSPLIT(GSTATS(IG)%NYH,IX) .EQ. IG )                 &
                  MSPLIT(GSTATS(IG)%NYH,IX) = -1
              END DO
!
! 3.  Bottom --------------------------------------------------------- *
!
            DO IX=1, NX
              IF (MSPLIT(GSTATS(IG)%NYL,IX) .EQ. IG )                 &
                  MSPLIT(GSTATS(IG)%NYL,IX) = -1
              END DO
!
          END IF
!
! 4.  Left ----------------------------------------------------------- *
!
        IF ( MX .GE. 5 ) THEN
!
            DO IY=GSTATS(IG)%NYL, GSTATS(IG)%NYH
              IF (MSPLIT(IY,GSTATS(IG)%NXL) .EQ. IG )                 &
                  MSPLIT(IY,GSTATS(IG)%NXL) = -1
              END DO
!
! 5.  Right ---------------------------------------------------------- *
!
            DO IY=GSTATS(IG)%NYH, GSTATS(IG)%NYH
              IF (MSPLIT(IY,GSTATS(IG)%NXH) .EQ. IG )                 &
                  MSPLIT(IY,GSTATS(IG)%NXH) = -1
              END DO
!
          END IF
!
! ... End loop started in 1.
!
        END DO
!
      RETURN
!
! Formats
!
!/ End of GRSQRG ----------------------------------------------------- /
!/
      END SUBROUTINE GRSQRG
!/ ------------------------------------------------------------------- /
      SUBROUTINE GRSNGL ( OK )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         09-Sep-2012 |
!/                  +-----------------------------------+
!/
!/    09-Sep-2012 : Origination.                        ( version 4.10 )
!/
!  1. Purpose :
!
!     Remove points from a grid that are in the middle of the sea, but
!     that have omly one adjacent point in the same grid. Directly
!     select a new grid for this point rather than deactivate and use
!     GRFILL.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       OK      Log.  I/O  Flag for grid status, .F. if values of
!                          -1 are left in MSPLIT.
!     ----------------------------------------------------------------
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      LOGICAL, INTENT(INOUT)  :: OK
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: NX0, NXN, IXL, IXH, COUNT(-1:NG),    &
                                 INEW1, INEW2, INEW
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Set up looping ------------------------------------------------- *
!
      IF ( GLOBAL ) THEN
          NX0    =  1
          NXN    = NX
        ELSE
          NX0    =  2
          NXN    = NX-1
        END IF
!
! 2.  Loops over 2D grid --------------------------------------------- *
!
      DO IX=NX0, NXN
!
        IXL    = IX - 1
        IXH    = IX + 1
        IF ( IX .EQ.  1 ) IXL = NX
        IF ( IX .EQ. NX ) IXH =  1
!
        DO IY=2, NY-1
!
! 3.  Central sea points only ---------------------------------------- *
!
          IF ( SEA(IY,IX) .AND. SEA(IY-1,IX ) .AND. SEA(IY+1,IX )      &
                          .AND. SEA(IY  ,IXL) .AND. SEA(IY  ,IXH) ) THEN
!
! 4.  Check for 'lost points' ---------------------------------------- *
!
              COUNT  = 0
              IG     = MSPLIT(IY,IX)
!
              COUNT(IG) = 1
              COUNT(MSPLIT(IY-1,IX )) = COUNT(MSPLIT(IY-1,IX )) + 1
              COUNT(MSPLIT(IY+1,IX )) = COUNT(MSPLIT(IY+1,IX )) + 1
              COUNT(MSPLIT(IY  ,IXL)) = COUNT(MSPLIT(IY  ,IXL)) + 1
              COUNT(MSPLIT(IY  ,IXH)) = COUNT(MSPLIT(IY  ,IXH)) + 1
!
              IF ( COUNT(IG) .LE. 2 ) THEN
!
                  INEW1  = -1
                  INEW2  = -1
!
                  DO J=1, NG
                    IF ( COUNT(J) .GE. 2 ) THEN
                        IF ( INEW1 .EQ. -1 ) THEN
                            INEW1  = J
                          ELSE
                            INEW2  = J
                            EXIT
                          END IF
                      END IF
                    END DO
!
                 IF ( INEW1 .EQ. -1 ) THEN
                     INEW   = -1
                     OK     = .FALSE.
                   ELSE IF ( INEW2 .EQ. -1 ) THEN
                     INEW   = INEW1
                   ELSE
                     IF ( GSTATS(INEW1)%NPTS .GT.                     &
                          GSTATS(INEW2)%NPTS ) THEN
                         INEW   = INEW2
                       ELSE
                         INEW   = INEW1
                       END IF
                   END IF
!
                MSPLIT(IY,IX) = INEW
!
                END IF
!
            END IF
!
! ... End loops started in 2.
!
          END DO
!
        END DO
!
      RETURN
!
! Formats
!
!/ End of GRSNGL ----------------------------------------------------- /
!/
      END SUBROUTINE GRSNGL
!/ ------------------------------------------------------------------- /
      SUBROUTINE GRSEPA ( OK, FRAC )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Feb-2013 |
!/                  +-----------------------------------+
!/
!/    10-Sep-2012 : Origination.                        ( version 4.10 )
!/    18-Sep-2012 : Include edge points of grid.        ( version 4.10 )
!/    01-Feb-2013 : Much faster algorithms.             ( version 4.10 )
!/
!  1. Purpose :
!
!     Remove smller parts of a grid that are separated from the main
!     body, and that can be attached to other grids.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       OK      Log.  I/O  Flag for grid status, .F. if values of
!                          -1 are left in MSPLIT.
!       FRAC    Real   I   Fraction of average size used to remove grid
!                          part.
!     ----------------------------------------------------------------
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: FRAC
      LOGICAL, INTENT(INOUT)  :: OK
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IPAVG, IPCHCK, ID, IPTOT, IX, IY,    &
                                 IXL, IYL, IDL, JX, JY, KY, IPT,      &
                                 IXH, IYH, I, J, K, L, IMIN, LMIN
      INTEGER                 :: GMASK(NY,NX), IIX(NSEA), IIY(NSEA)
      INTEGER, ALLOCATABLE    :: PMAP(:), INGRD(:)
      LOGICAL                 :: PREV
      LOGICAL,ALLOCATABLE     :: FLNEXT(:), NEXTTO(:,:)
!/
!/ ------------------------------------------------------------------- /
!/
!
      IPAVG  = NINT ( REAL(NSEA) / REAL(NG) )
      IPCHCK = NINT ( FRAC * REAL(NSEA) / REAL(NG) )
!
! 1.  Loop over grids ------------------------------------------------ *
!
      DO IG=1, NG
!
        GMASK  = 0
        ID     = 0
!
! 2.  Find all parts ------------------------------------------------- *
! 2.a First loop, partial parts
!
        IPTOT  = 0
!
        DO IX=1, NX
!
          IXL    = 1 + MOD(IX-2+NX,NX)
          PREV   = .FALSE.
!
          DO IY=1, NY
            IF (MSPLIT(IY,IX) .EQ. IG ) THEN
                IPTOT  = IPTOT + 1
                IIX(IPTOT) = IX
                IIY(IPTOT) = IY
                IF ( .NOT. PREV) THEN
                    ID     = ID + 1
                    PREV   = .TRUE.
                  END IF
                GMASK(IY,IX) = ID
              ELSE IF ( PREV ) THEN
                PREV   = .FALSE.
                IDL    = 0
                DO JY=IY-1, 1, -1
                  IF ( GMASK(JY,IX) .EQ. 0 ) EXIT
                  IF ( GMASK(JY,IXL).NE.0 .AND. IDL.EQ.0 )            &
                                                 IDL = GMASK(JY,IXL)
                  END DO
                IF ( IDL .NE. 0 ) THEN
                    DO KY=JY+1, IY-1
                      IF ( GMASK(KY,IX).EQ.ID ) GMASK(KY,IX) = IDL
                      END DO
                    ID     = ID - 1
                  END IF
!
              END  IF
            END DO
          END DO
!
! 2.b Grid too small, do not cut
!
        IF ( IPTOT .LE. IPAVG ) THEN
            CYCLE
          END IF
!
! 2.c Neighbouring grid parts
!     Raw data
!
        ALLOCATE ( NEXTTO(0:ID,0:ID), PMAP(0:ID) )
        NEXTTO = .FALSE.
!
        DO IPT=1, IPTOT
          IX     = IIX(IPT)
          IY     = IIY(IPT)
          IXL    = 1 + MOD(IX-2+NX,NX)
          IYL    = IY - 1
          IXH    = 1 + MOD(IX,NX)
          IYH    = IY + 1
          NEXTTO( GMASK(IY,IX) , GMASK(IY ,IXL) ) = .TRUE.
          NEXTTO( GMASK(IY,IX) , GMASK(IY ,IXH) ) = .TRUE.
          NEXTTO( GMASK(IY,IX) , GMASK(IYL,IX ) ) = .TRUE.
          NEXTTO( GMASK(IY,IX) , GMASK(IYH,IX ) ) = .TRUE.
          END DO
!
!     Make symmetric
!
        DO I=1, ID
          DO J=1, ID
            NEXTTO(I,J) = NEXTTO(I,J) .OR. NEXTTO(J,I)
            END DO
          END DO
!
!     Connect accross neighbours
!
        DO I=1, ID
          DO J=1, ID
            IF ( NEXTTO(I,J) ) THEN
                DO K=1, ID
                  IF ( NEXTTO(K,J) ) THEN
                      NEXTTO(K,I) = .TRUE.
                      NEXTTO(I,K) = .TRUE.
                    END IF
                  END DO
              END IF
            END DO
          END DO
!
!     Map the parts
!
        IDL    = ID
        PMAP   = 0
        ID     = 0
!
        DO I=1, IDL
          IF ( PMAP(I) .EQ. 0 ) THEN
              ID     = ID + 1
              DO J=1, IDL
                IF ( NEXTTO(J,I) ) EXIT
                END DO
              IF ( J .GT. IDL ) THEN
                  PMAP(I) = ID
                ELSE
                  DO K=I, IDL
                    IF ( PMAP(K).EQ.0 .AND. NEXTTO(J,K) ) PMAP(K) = ID
                    END DO
                END IF
            END IF
          END DO
!
        DEALLOCATE ( NEXTTO )
!
! 3.  Grid is contiguous --------------------------------------------- *
!
        IF ( ID .EQ. 1 ) THEN
            DEALLOCATE ( PMAP )
            CYCLE
          END IF
!
! 4.  Grid is split, get stats --------------------------------------- *
!
! 4.a Construct final map for grid
!
        DO IPT=1, IPTOT
          IX     = IIX(IPT)
          IY     = IIY(IPT)
          GMASK(IY,IX) = PMAP(GMASK(IY,IX))
          END DO
!
        DEALLOCATE ( PMAP )
!
! 4.b Run stats
!
        ALLOCATE ( INGRD(ID), FLNEXT(ID) )
        INGRD  = 0
        FLNEXT = .FALSE.
        IPTOT  = 0
!
        DO JX=1, NX
          DO JY=1, NY
            IF ( GMASK(JY,JX) .GT. 0 ) THEN
                INGRD(GMASK(JY,JX))  = INGRD(GMASK(JY,JX)) + 1
                IPTOT                = IPTOT + 1
              END IF
            END DO
          END DO
!
        DO JX=1, NX
          DO JY=1, NY-1
            IF ( ( GMASK(JY  ,JX) .GT. 0 ) .AND.                      &
                 ( SEA(JY+1,JX) .AND. MSPLIT(JY+1,JX).NE.IG ) )       &
                FLNEXT(GMASK(JY  ,JX)) = .TRUE.
            IF ( ( GMASK(JY+1,JX) .GT. 0 ) .AND.                      &
                 ( SEA(JY  ,JX) .AND. MSPLIT(JY  ,JX).NE.IG ) )       &
                FLNEXT(GMASK(JY+1,JX)) = .TRUE.
            END DO
          END DO
!
        DO JY=1, NY
          DO JX=1, NX-1
            IF ( ( GMASK(JY,JX  ) .GT. 0 ) .AND.                      &
                 ( SEA(JY,JX+1) .AND. MSPLIT(JY,JX+1).NE.IG ) )       &
                FLNEXT(GMASK(JY,JX  )) = .TRUE.
            IF ( ( GMASK(JY,JX+1) .GT. 0 ) .AND.                      &
                 ( SEA(JY,JX  ) .AND. MSPLIT(JY,JX  ).NE.IG ) )       &
                FLNEXT(GMASK(JY,JX+1)) = .TRUE.
            END DO
          IF ( GLOBAL ) THEN
              IF ( ( GMASK(JY,NX) .GT. 0 ) .AND.                      &
                   ( SEA(JY, 1) .AND. MSPLIT(JY, 1).NE.IG ) )         &
                  FLNEXT(GMASK(JY,NX)) = .TRUE.
              IF ( ( GMASK(JY, 1) .GT. 0 ) .AND.                      &
                   ( SEA(JY,NX) .AND. MSPLIT(JY,NX).NE.IG ) )         &
                  FLNEXT(GMASK(JY, 1)) = .TRUE.
            END IF
          END DO
!
! 5.  Grid large enough, find smallest part -------------------------- *
!
        IMIN   = NSEA
        LMIN   = 0
!
        DO J=1, ID
          IF ( FLNEXT(J) .AND. INGRD(J).LT.IMIN ) THEN
              IMIN   = INGRD(J)
              LMIN   = J
            END IF
          END DO
!
        IF ( LMIN .EQ. 0 ) THEN
            DEALLOCATE ( INGRD, FLNEXT )
            CYCLE
          END IF
!
        IF ( IMIN .GT. IPCHCK ) THEN
            DEALLOCATE ( INGRD, FLNEXT )
            CYCLE
          END IF
!
! 6.  Part to cut has been identified -------------------------------- *
!
        DO JX=1, NX
          DO JY=1, NY
              IF ( GMASK(JY,JX) .EQ. LMIN ) MSPLIT(JY,JX) = -1
            END DO
          END DO
!
        DEALLOCATE ( INGRD, FLNEXT )
        OK     = .FALSE.
!
! ... End loops started in 1.
!
        END DO
!
      RETURN
!
! Formats
!
!/ End of GRSEPA ----------------------------------------------------- /
!/
      END SUBROUTINE GRSEPA
!/ ------------------------------------------------------------------- /
      SUBROUTINE GRFSML
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         04-Feb-2013 |
!/                  +-----------------------------------+
!/
!/    13-Sep-2012 : Origination.                        ( version 4.10 )
!/    04-Feb-2013 : Bug fix grid splitting.             ( version 4.10 )
!/
!  1. Purpose :
!
!     Subroutine called when lowest grid size is stuck. Attempting to
!     joint to neighbor grid, otherwise mark as accepted small grid.
!     note that small grid does not influence parallel scaling like a
!     big grid does .....
!
!     1-Feb-2013: Also used for early small-grid merging.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: NSMALL, IGMIN(NG), NNEXT, JG, IGADD, &
                                 IGTEST, FREE(NG), NFREE, NBIG, IGB,  &
                                 MX, MY, NX0, NXN, NY0, NYN, JX
      CHARACTER(LEN=1)        :: NEXTTO(0:NG,0:NG), TEMP(NG)
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Find small(s) -------------------------------------------------- *
!
      NSMALL = 0
      IGMIN  = 0
!
      DO IG=1,NG
        IF ( GSTATS(IG)%INSTAT .AND.                                  &
             GSTATS(IG)%NPTS .EQ. MSTATS%NMIN ) THEN
            NSMALL = NSMALL + 1
            IGMIN(NSMALL)  = IG
          END IF
        END DO
!
! 2.  Find neighbours ------------------------------------------------ *
!
      NEXTTO = '.'
!
      DO IX=1, NX-1
        DO IY=1, NY-1
          NEXTTO(MSPLIT(IY  ,IX  ),MSPLIT(IY+1,IX  )) = 'X'
          NEXTTO(MSPLIT(IY+1,IX  ),MSPLIT(IY  ,IX  )) = 'X'
          NEXTTO(MSPLIT(IY  ,IX+1),MSPLIT(IY  ,IX  )) = 'X'
          NEXTTO(MSPLIT(IY  ,IX  ),MSPLIT(IY  ,IX+1)) = 'X'
          END DO
        END DO
!
      IF ( GLOBAL ) THEN
          DO IY=1, NY-1
            NEXTTO(MSPLIT(IY  ,NX),MSPLIT(IY+1,NX)) = 'X'
            NEXTTO(MSPLIT(IY+1,NX),MSPLIT(IY  ,NX)) = 'X'
            NEXTTO(MSPLIT(IY  , 1),MSPLIT(IY  ,NX)) = 'X'
            NEXTTO(MSPLIT(IY  ,NX),MSPLIT(IY  , 1)) = 'X'
            END DO
        END IF
!
      DO IG=0,NG
        NEXTTO(IG,IG) = '-'
        END DO
!
! 3.  Loop over small grids ------------------------------------------ *
!
      FREE   = 0
      NFREE  = 0
!
      DO J=1, NSMALL
!
! 3.a Find neighbours
!
        IG     = IGMIN(J)
        IGADD   = 0
        IGTEST  = NSEA + 1
        NNEXT  = 0
        DO JG=1, NG
          IF ( NEXTTO(IG,JG) .EQ. 'X' ) THEN
              NNEXT  = NNEXT + 1
              IF ( GSTATS(JG)%NPTS .LT. IGTEST ) THEN
                  IGTEST = GSTATS(JG)%NPTS
                  IGADD  = JG
                END IF
            END IF
          END DO
!
! 3.b No neighbours found, mark as 'not to be processed further'
!
        IF ( NNEXT .EQ. 0 ) THEN
            GSTATS(IG)%INSTAT = .FALSE.
          ELSE
!
! 3.c Check smallest neighbor
!
            IF ( IGTEST + INGMIN .LT. NINT(XMEAN) ) THEN
!
! ... Merge grids
!
                DO IX=1, NX
                 DO IY=1, NY
                   IF ( MSPLIT(IY,IX) .EQ. IG ) MSPLIT(IY,IX) = IGADD
                   END DO
                 END DO
!
                NFREE  = NFREE + 1
                FREE(NFREE) = IG
!
              ELSE
!
! ... Remove grid(s) from stats
!
                GSTATS(IG)%INSTAT = .FALSE.
                NNEXT  = 0
                DO JG=1, NG
                  IF ( NEXTTO(IGADD,JG) .EQ. 'X' ) NNEXT  = NNEXT + 1
                  END DO
                IF ( NNEXT .EQ. 1 ) THEN
                    GSTATS(IGADD)%INSTAT = .FALSE.
                  END IF
!
              END IF
!
          END IF
!
        END DO
!
! 4.  Make new grids as needed --------------------------------------- *
!
      DO J=1, NFREE
!
! 4.a Find biggest grid
!
        NBIG   = 0
        IGB    = 0
!
        DO IG=1, NG
          IF ( GSTATS(IG)%NPTS .GT. NBIG ) THEN
              NBIG   = GSTATS(IG)%NPTS
              IGB    = IG
            END IF
          END DO
!
! 4.a Split biggest grid
!
        NX0    = GSTATS(IGB)%NXL
        NXN    = GSTATS(IGB)%NXH
        NY0    = GSTATS(IGB)%NYL
        NYN    = GSTATS(IGB)%NYH
!
        MY     = 1 + GSTATS(IGB)%NYH - GSTATS(IGB)%NYL
        MX     = 1 + GSTATS(IGB)%NXH - GSTATS(IGB)%NXL
        IF ( GSTATS(IGB)%STRADLE ) MX     = MX + NX
!
        IF ( MY .GE. MX ) THEN
            NYN    = NY0 + MY/2
          ELSE
            NXN    = NX0 + MX/2
          END IF
!
        DO IX=NX0, NXN
          JX     = 1 + MOD(IX-1,NX)
          DO IY=NY0, NYN
            IF ( MSPLIT(IY,JX) .EQ. IGB ) MSPLIT(IY,JX) = FREE(J)
            END DO
          END DO
!
        GSTATS(IGB)%NPTS = 0
        GSTATS(FREE(J))%NPTS = 0
!
        END DO
!
      RETURN
!
! Formats
!
!/ End of GRFSML ----------------------------------------------------- /
!/
      END SUBROUTINE GRFSML
!/ ------------------------------------------------------------------- /
      SUBROUTINE GRFLRG
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-Jan-2013 |
!/                  +-----------------------------------+
!/
!/    19-Sep-2012 : Origination.                        ( version 4.10 )
!/    29-Jan-2013 : Add error code on stop.             ( version 4.10 )
!/
!  1. Purpose :
!
!     Like GRFSML for largest grid ...
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: NBIG, IGMAX(NG), NNEXT, JG
!!!      INTEGER                 :: NSMALL, IGMIN(NG), NNEXT, JG, IGADD, &
!!!                                 IGTEST, FREE(NG), NFREE, NBIG, IGB,  &
!!!                                 MX, MY, NX0, NXN, NY0, NYN
!!!!/S      INTEGER, SAVE           :: IENT = 0
      CHARACTER(LEN=1)        :: NEXTTO(0:NG,0:NG), TEMP(NG)
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Find big(s) ---------------------------------------------------- *
!
      NBIG   = 0
      IGMAX  = 0
!
      DO IG=1,NG
        IF ( GSTATS(IG)%INSTAT .AND.                                  &
             GSTATS(IG)%NPTS .EQ. MSTATS%NMAX ) THEN
            NBIG   = NBIG + 1
            IGMAX(NBIG)  = IG
          END IF
        END DO
!
! 2.  Find neighbours ------------------------------------------------ *
!
      NEXTTO = '.'
!
      DO IX=1, NX-1
        DO IY=1, NY-1
          NEXTTO(MSPLIT(IY  ,IX  ),MSPLIT(IY+1,IX  )) = 'X'
          NEXTTO(MSPLIT(IY+1,IX  ),MSPLIT(IY  ,IX  )) = 'X'
          NEXTTO(MSPLIT(IY  ,IX+1),MSPLIT(IY  ,IX  )) = 'X'
          NEXTTO(MSPLIT(IY  ,IX  ),MSPLIT(IY  ,IX+1)) = 'X'
          END DO
        END DO
!
      IF ( GLOBAL ) THEN
          DO IY=1, NY-1
            NEXTTO(MSPLIT(IY  ,NX),MSPLIT(IY+1,NX)) = 'X'
            NEXTTO(MSPLIT(IY+1,NX),MSPLIT(IY  ,NX)) = 'X'
            NEXTTO(MSPLIT(IY  , 1),MSPLIT(IY  ,NX)) = 'X'
            NEXTTO(MSPLIT(IY  ,NX),MSPLIT(IY  , 1)) = 'X'
            END DO
        END IF
!
      DO IG=0,NG
        NEXTTO(IG,IG) = '-'
        END DO
!
! 3.  Loop over big grids -------------------------------------------- *
!
      DO J=1, NBIG
!
! 3.a Find neighbours
!
        IG     = IGMAX(J)
        NNEXT  = 0
        DO JG=1, NG
          IF ( NEXTTO(IG,JG) .EQ. 'X' ) NNEXT  = NNEXT + 1
          END DO
!
! 3.b Enough neighbours found, mark as 'not to be processed further'
!
        IF ( NNEXT .GE. 1 ) THEN
            GSTATS(IG)%INSTAT = .FALSE.
          ELSE
!
! 3.c Biggest grid is isolated, should split
!
            WRITE (NDSE,930)
            STOP 11
!
          END IF
!
        END DO
!
      RETURN
!
! Formats
!
  930 FORMAT ( ' *** ERROR GRFLRG: LARGEST GRID IS ISOLATED ***'      &
               '     SPLITTING NOT YET IMPLEMENTED '/)
!
!/ End of GRFLRG ----------------------------------------------------- /
!/
      END SUBROUTINE GRFLRG
!/ ------------------------------------------------------------------- /
      SUBROUTINE GR1GRD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         18-Nov-2012 |
!/                  +-----------------------------------+
!/
!/    23-Sep-2012 : Origination.                        ( version 4.10 )
!/    24-Jan-2013 : Correct X0 to be in range.          ( version 4.10 )
!/    04-Feb-2013 : Add corner point to halo.           ( version 4.10 )
!/    18-Nov-2012 : Add user-defined halo extension.    ( version 4.14 )
!/
!  1. Purpose :
!
!     Extract single grid from master map, including halo needed for
!     grid overlap in ww3_multi.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: NIT, IIT,  IXL, IXH, IYL, IYH, NOCNT,&
                                 NOCNTM, NOCNTL, JX, JY, ISEA, MX, MY
      INTEGER                 :: MTMP2(NY,NX)
      REAL                    :: XOFF
      LOGICAL                 :: MASK(NY,NX), LEFT, RIGHT, THERE
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Set up MTEMP with MAPSTA 0,1,3 for grid ------------------------ *
!
      DO IX=1, NX
        DO IY=1, NY
          IF ( MSPLIT(IY,IX) .EQ. IG ) THEN
              MTEMP(IY,IX) = 1
            ELSE IF ( MSPLIT(IY,IX) .GT. 0 ) THEN
              MTEMP(IY,IX) = 3
            ELSE
              MTEMP(IY,IX) = 0
            END IF
          END DO
        END DO
!
! 2.  Add ALL MAPSTA = 2 points to grid ------------------------------ *
!
      DO IX=1, NX
        DO IY=1, NY
          IF ( MAPSTA(IY,IX) .EQ. 2 ) THEN
              MTEMP(IY,IX) = 2
            END IF
          END DO
        END DO
!
! 3.  Add halo ------------------------------------------------------- *
! 3.a Set up halo width depending on scheme and time steps
!     NEEDED TO SET UP A LITTLE WIDER. NOT SURE WHY. NEED TO CHECK WITH
!     WMEQL SUBROUTINE.
!
      NIT    = 1 + NHEXT + ( 1 + INT(DTMAX/DTCFL-0.001) ) * 3
!
! 3.b Exand halo
!
      DO IIT=1, NIT
!
        MASK   = .FALSE.
!
        DO IX=1, NX
          IXL    = 1 + MOD(IX-2+NX,NX)
          IXH    = 1 + MOD(IX,NX)
          DO IY=2, NY-1
            IF ( MTEMP(IY,IX) .EQ. 3 ) MASK(IY,IX) =                  &
                              ( ( MTEMP(IY+1,IX ) .EQ. 1 ) .OR.       &
                                ( MTEMP(IY-1,IX ) .EQ. 1 ) .OR.       &
                                ( MTEMP(IY  ,IXH) .EQ. 1 ) .OR.       &
                                ( MTEMP(IY  ,IXL) .EQ. 1 ) )          &
                    .OR.    (   ( MTEMP(IY+1,IXL) .EQ. 1 ) .AND.      &
                              ( ( MTEMP(IY  ,IXL) .EQ. 1 ) .OR.       &
                                ( MTEMP(IY+1,IX ) .EQ. 1 ) ) )        &
                    .OR.    (   ( MTEMP(IY+1,IXH) .EQ. 1 ) .AND.      &
                              ( ( MTEMP(IY  ,IXH) .EQ. 1 ) .OR.       &
                                ( MTEMP(IY+1,IX ) .EQ. 1 ) ) )        &
                    .OR.    (   ( MTEMP(IY-1,IXH) .EQ. 1 ) .AND.      &
                              ( ( MTEMP(IY  ,IXH) .EQ. 1 ) .OR.       &
                                ( MTEMP(IY-1,IX ) .EQ. 1 ) ) )        &
                    .OR.    (   ( MTEMP(IY-1,IXL) .EQ. 1 ) .AND.      &
                              ( ( MTEMP(IY  ,IXL) .EQ. 1 ) .OR.       &
                                ( MTEMP(IY-1,IX ) .EQ. 1 ) ) )
            END DO
          END DO
!
        DO IX=1, NX
          DO IY=1, NY
            IF ( MASK(IY,IX) ) MTEMP(IY,IX) = 1
            END DO
          END DO
!
        END DO
!
! 3.c Contract halo
!
!     MTMP2 = MTEMP
!
!     DO IIT=1, NIT
!
!       MASK   = .FALSE.
!
!       DO IX=1, NX
!         IXL    = 1 + MOD(IX-2+NX,NX)
!         IXH    = 1 + MOD(IX,NX)
!         DO IY=2, NY-1
!           IF ( MTMP2(IY,IX) .EQ. 1 ) MASK(IY,IX) =                  &
!                             ( ( MTMP2(IY+1,IX ) .EQ. 3 ) .OR.       &
!                               ( MTMP2(IY-1,IX ) .EQ. 3 ) .OR.       &
!                               ( MTMP2(IY  ,IXH) .EQ. 3 ) .OR.       &
!                               ( MTMP2(IY  ,IXL) .EQ. 3 ) )          &
!                   .OR.    (   ( MTMP2(IY+1,IXL) .EQ. 3 ) .AND.      &
!                             ( ( MTMP2(IY  ,IXL) .EQ. 3 ) .OR.       &
!                               ( MTMP2(IY+1,IX ) .EQ. 3 ) ) )        &
!                   .OR.    (   ( MTMP2(IY+1,IXH) .EQ. 3 ) .AND.      &
!                             ( ( MTMP2(IY  ,IXH) .EQ. 3 ) .OR.       &
!                               ( MTMP2(IY+1,IX ) .EQ. 3 ) ) )        &
!                   .OR.    (   ( MTMP2(IY-1,IXH) .EQ. 3 ) .AND.      &
!                             ( ( MTMP2(IY  ,IXH) .EQ. 3 ) .OR.       &
!                               ( MTMP2(IY-1,IX ) .EQ. 3 ) ) )        &
!                   .OR.    (   ( MTMP2(IY-1,IXL) .EQ. 3 ) .AND.      &
!                             ( ( MTMP2(IY  ,IXL) .EQ. 3 ) .OR.       &
!                               ( MTMP2(IY-1,IX ) .EQ. 3 ) ) )
!           END DO
!         END DO
!
!       DO IX=1, NX
!         DO IY=1, NY
!           IF ( MASK(IY,IX) ) MTMP2(IY,IX) = 3
!           END DO
!         END DO
!
!       END DO
!
! 3.d Check if consistent .....
!
!     DO IX=1, NX
!       DO IY=1, NY
!         IF ( MSPLIT(IY,IX).EQ.IG .OR. MTMP2(IY,IX).EQ.1 ) THEN
!             IF ( MSPLIT(IY,IX).EQ.IG .AND. MTMP2(IY,IX).NE.1 ) THEN
!                 write (ndst,*) ix, iy, '  in grid, not in e-c grid'
!               END IF
!             IF ( MSPLIT(IY,IX).NE.IG .AND. MTMP2(IY,IX).EQ.1 ) THEN
!                 write (ndst,*) ix, iy, '  in e-c grid, not in grid'
!               END IF
!           END IF
!         END DO
!       END DO
!
! 4.  Remove extraeneous MAPSTA = 2 ---------------------------------- *
!
      DO IX=1, NX
!
        IF ( GLOBAL ) THEN
            IXL    = 1 + MOD(IX-2+NX,NX)
            IXH    = 1 + MOD(IX,NX)
          ELSE
            IXL    = MAX ( 1 , IX-1 )
            IXH    = MIN ( NX , IX+1 )
          END IF
!
        DO IY=1, NY
          IF ( MTEMP(IY,IX) .EQ. 2 ) THEN
              IYL    = MAX ( 1 , IY-1 )
              IYH    = MIN ( NY , IY+1 )
              IF ( .NOT. ( ( MTEMP(IYL,IX ) .EQ. 1 ) .OR.             &
                           ( MTEMP(IYH,IX ) .EQ. 1 ) .OR.             &
                           ( MTEMP(IY ,IXL) .EQ. 1 ) .OR.             &
                           ( MTEMP(IY ,IXH) .EQ. 1 ) ) )              &
                   MTEMP(IY,IX) = 3
            END IF
          END DO
!
        END DO
!
! 5.  Recompute grid range ------------------------------------------- *
!     Using GSTOLD to store info for modified grid
!
      GSTOLD(IG)%STRADLE = .FALSE.
      GSTOLD(IG)%NPTS    = 0
      GSTOLD(IG)%NXL     = NX
      GSTOLD(IG)%NXH     = 1
      GSTOLD(IG)%NYL     = NY
      GSTOLD(IG)%NYH     = 1
!
      IF ( GLOBAL ) THEN
!
          LEFT   = .FALSE.
          RIGHT  = .FALSE.
!
          DO IY=1, NY
            IF ( MTEMP(IY, 1).EQ.1 .OR.  MTEMP(IY, 1).EQ.2 ) LEFT   = .TRUE.
            IF ( MTEMP(IY,NX).EQ.1 .OR.  MTEMP(IY,NX).EQ.2 ) RIGHT  = .TRUE.
            END DO
          GSTOLD(IG)%STRADLE = LEFT .AND. RIGHT
!
        END IF
!
      DO IY=1, NY
        DO IX=1, NX
          IF ( MTEMP(IY,IX).EQ.1 .OR.  MTEMP(IY,IX).EQ.2 ) THEN
              GSTOLD(IG)%NPTS = GSTOLD(IG)%NPTS + 1
              GSTOLD(IG)%NXL  = MIN ( GSTOLD(IG)%NXL , IX )
              GSTOLD(IG)%NXH  = MAX ( GSTOLD(IG)%NXH , IX )
              GSTOLD(IG)%NYL  = MIN ( GSTOLD(IG)%NYL , IY )
              GSTOLD(IG)%NYH  = MAX ( GSTOLD(IG)%NYH , IY )
            END IF
          END DO
        END DO
!
      IF ( GSTOLD(IG)%STRADLE ) THEN
          NOCNT  = 0
          NOCNTM = 0
          NOCNTL = 0
          DO IX=1, NX
            THERE  = .FALSE.
            DO IY=1, NY
              IF ( MTEMP(IY,IX).EQ.1 .OR.  MTEMP(IY,IX).EQ.2 ) THEN
                  THERE  = .TRUE.
                  EXIT
                END IF
              END DO
            IF ( THERE ) THEN
                NOCNT = 0
              ELSE
                NOCNT = NOCNT + 1
                IF ( NOCNT .GT. NOCNTM ) THEN
                    NOCNTM = NOCNT
                    NOCNTL = IX
                  END IF
              END IF
            END DO
          GSTOLD(IG)%NXL = NOCNTL + 1
          GSTOLD(IG)%NXH = NOCNTL - NOCNTM
        END IF
!
! ... Make sure outside of grid is 2 or 3
!
      LEFT   = .FALSE.
      RIGHT  = .FALSE.
!
      DO IX=1, NX
        LEFT   = LEFT  .OR. ( MTEMP(GSTOLD(IG)%NYL,IX) .EQ. 1 )
        RIGHT  = RIGHT .OR. ( MTEMP(GSTOLD(IG)%NYH,IX) .EQ. 1 )
        END DO
!
      IF ( LEFT  ) GSTOLD(IG)%NYL = GSTOLD(IG)%NYL - 1
      IF ( RIGHT ) GSTOLD(IG)%NYH = GSTOLD(IG)%NYH + 1
!
      DO IY=1, NY
        LEFT   = LEFT  .OR. ( MTEMP(IY,GSTOLD(IG)%NXL) .EQ. 1 )
        RIGHT  = RIGHT .OR. ( MTEMP(IY,GSTOLD(IG)%NXH) .EQ. 1 )
        END DO
!
      IF ( LEFT  ) GSTOLD(IG)%NXL = GSTOLD(IG)%NXL - 1
      IF ( RIGHT ) GSTOLD(IG)%NXH = GSTOLD(IG)%NXH + 1
!
      IF ( GLOBAL .AND. GSTOLD(IG)%NXL.EQ.0 ) THEN
          GSTOLD(IG)%NXL     = NX
          GSTOLD(IG)%STRADLE = .TRUE.
        END IF
!
      IF ( GLOBAL .AND. GSTOLD(IG)%NXH.EQ.NX+1 ) THEN
          GSTOLD(IG)%NXH     = 1
          GSTOLD(IG)%STRADLE = .TRUE.
        END IF
!
! 6.  Extract reduced grid data -------------------------------------- *
!
      MY             = 1 + GSTOLD(IG)%NYH - GSTOLD(IG)%NYL
      MX             = 1 + GSTOLD(IG)%NXH - GSTOLD(IG)%NXL
      IF ( GSTOLD(IG)%STRADLE ) MX = MX + NX
      PGRID(IG)%NY   = MY
      PGRID(IG)%NX   = MX
      PGRID(IG)%NSEA = GSTOLD(IG)%NPTS
      PGRID(IG)%X0   = X0 + REAL(GSTOLD(IG)%NXL-1)*SX
      PGRID(IG)%Y0   = Y0 + REAL(GSTOLD(IG)%NYL-1)*SY
      PGRID(IG)%SX   = SX
      PGRID(IG)%SY   = SY
!
      XOFF   = 360. * REAL ( NINT((PGRID(IG)%X0+0.5*REAL(MX-1)*SX)/360.) )
      PGRID(IG)%X0   = PGRID(IG)%X0 - XOFF
!
      ALLOCATE ( PGRID(IG)%ZBIN(MX,MY) ,          &
                 PGRID(IG)%OBSX(MX,MY) ,          &
                 PGRID(IG)%OBSY(MX,MY) ,          &
                 PGRID(IG)%MASK(MX,MY) )
!
      PGRID(IG)%ZBIN   = ZBDUM
      PGRID(IG)%OBSX   = 0.
      PGRID(IG)%OBSY   = 0.
      PGRID(IG)%MASK   = 99
!
      DO IX=1, PGRID(IG)%NX
        JX     = 1 + MOD ( IX+GSTOLD(IG)%NXL-2 , NX )
        DO IY=1, PGRID(IG)%NY
          JY     = IY + GSTOLD(IG)%NYL - 1
          ISEA   = MAPFS(JY,JX)
          IF ( MTEMP(JY,JX) .NE. 0 ) THEN
              PGRID(IG)%ZBIN(IX,IY) = ZB(ISEA)
            END IF
          IF ( TRFLAG .NE. 0 ) THEN
              PGRID(IG)%OBSX(IX,IY) = 1. - TRNX(JY,JX)
              PGRID(IG)%OBSY(IX,IY) = 1. - TRNY(JY,JX)
            END IF
          PGRID(IG)%MASK(IX,IY) = MTEMP(JY,JX)
          END DO
        END DO
!
      RETURN
!
! Formats
!
!/ End of GR1GRD ----------------------------------------------------- /
!/
      END SUBROUTINE GR1GRD
!/
!/ End of W3GSPL ----------------------------------------------------- /
!/
      END PROGRAM W3GSPL
