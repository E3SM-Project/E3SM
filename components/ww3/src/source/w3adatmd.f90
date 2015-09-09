!/ ------------------------------------------------------------------- /
      MODULE W3ADATMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    28-Dec-2004 : Origination.                        ( version 3.06 )
!/    04-May-2005 : Adding MPI_COMM_WAVE.               ( version 3.07 )
!/    20-Jul-2005 : Adding output fields.               ( version 3.07 )
!/    09-Nov-2005 : Removing soft boundary option.      ( version 3.08 )
!/    13-Jun-2006 : Splitting STORE in G/SSTORE.        ( version 3.09 )
!/    04-Oct-2006 : Add filter to array pointers.       ( version 3.10 )
!/    28_Mar-2007 : Add partitioned data arrays.        ( version 3.11 )
!/                  Add aditional undefined arrays.
!/    22-Feb-2008 ; Modify MAPTH2 declaration.          ( version 3.13 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Define data structures to set up wave model auxiliary daya for
!     several models simultaneously.
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      NADATA    Int.  Public   Number of models in array dim.
!      IADATA    Int.  Public   Selected model for output, init. at -1.
!      MPIBUF    I.P.  Public   Number of buffer arrays for 'hidden'
!                               MPI communications (no hiding for
!                               MPIBUF = 1).
!      WADAT     TYPE  Public   Basic data structure.
!      WADATS    WADAT Public   Array of data structures.
!     ----------------------------------------------------------------
!
!     All elements of WADAT are aliased to pointers with the same
!     name. These pointers are defined as :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!     Internal model definition:
!
!      CG        R.A.  Public   Group velocities for all wave model
!                               sea points and frequencies.
!      WN        R.A.  Public   Idem, wavenumbers.
!
!     Aux. arrays for model input:
!
!      CA0-I     R.A.  Public   Absolute current velocity (initial
!                               and inc.) in W3UCUR.
!      CD0-I     R.A.  Public   Current direction (initial and
!                               increment) in W3UCUR.
!      UA0-I     R.A.  Public   Absolute wind speeds (initial and
!                               incr.) in W3UWND                (m/s)
!      UD0-I     R.A.  Public   Wind direction (initial and incr.)
!                               in W3UWND                       (rad)
!      AS0-I     R.A.  Public   Stability par. (initial and incr.)
!                               in W3UWND                      (degr)
!      ATRNX/Y   R.A.  Public   Actual transparency info.
!
!     Fields of mean wave parameters:
!
!      DW        R.A.  Public   Water depths.
!      UA        R.A.  Public   Absolute wind speeds.
!      UD        R.A.  Public   Absolute wind direction.
!      U10       R.A.  Public   Wind speed used.
!      U10D      R.A.  Public   Wind direction used.
!      AS        R.A.  Public   Stability parameter.
!      CX/Y      R.A.  Public   Current components.
!      EMN       R.A.  Public   Mean energy.
!      FMN       R.A.  Public   Mean frequency.
!      WNM       R.A.  Public   Mean wavenumber.
!      AMX       R.A.  Public   Spectral maximum.
!      CDS       R.A.  Public   Drag coefficient.
!      Z0S       R.A.  Public   Roughness parameter.
!      HS        R.A.  Public   Wave Height.
!      WLM       R.A.  Public   Mean wave length.
!      TMN       R.A.  Public   Mean wave period.
!      THM       R.A.  Public   Mean wave direction.
!      THS       R.A.  Public   Mean directional spread.
!      FP0       R.A.  Public   Peak frequency.
!      THP0      R.A.  Public   Peak direction.
!      FP1       R.A.  Public   Wind sea peak frequency.
!      THP1      R.A.  Public   Wind sea peak direction.
!      DTDYN     R.A.  Public   Mean dynamic time step (raw).
!      FCUT      R.A.  Public   Cut-off frequency for tail.
!      ABA       R.A.  Public   Near-bottom rms wave ex. apmplitude.
!      ABD       R.A.  Public   Corresponding direction.
!      UBA       R.A.  Public   Near-bottom rms wave velocity.
!      UBD       R.A.  Public   Corresponding direction.
!      Sxx       R.A.  Public   Radiation stresses.
!      DDDx      R.A.  Public   Spatial derivatives of the depth.
!      DCxDx     R.A.  Public   Spatial dirivatives of the current.
!
!    Mean parameters from partitiones spectra, 2D array with el.
!    0 holding wind sea data, and 1:NOSWLL holding swell fields.
!    Last two arrays are regular single-entry arrays.
!
!      PHS       R.A.  Public   Wave height of partition.
!      PTP       R.A.  Public   Peak period of partition.
!      PLP       R.A.  Public   Peak wave leingth of partition.
!      PTH       R.A.  Public   Mean direciotn of partition.
!      PSI       R.A.  Public   Mean spread of partition.
!      PWS       R.A.  Public   Wind sea fraction of partition.
!
!      PWST      R.A.  Public   Total wind sea fractio
!      PNR       R.A.  Public   Number of partitions found.
!
!     Empty dummy fields (NOEXTR)
!
!      USERO     R.A.  Public   Empty output arrays than can be
!                               used by users as a simple means to
!                               add output.
!
!     Map data for propagation schemes (1Up).
!
!      IS0/2     I.A.  Public   Spectral propagation maps.
!      FACVX/Y   R.A.  Public   Spatial propagation factor map.
!
!     Map data for propagation schemes (UQ).
!
!      NMXn      Int.  Public    Counters for MAPX2, see W3MAP3.
!      NMYn      Int.  Public
!      NMXY      Int.  Public    Dimension of MAPXY.
!      NACTn     Int.  Public    Dimension of MAPAXY.
!      NCENT     Int.  Public    Dimension of MAPAXY.
!      MAPX2     I.A.  Public    Map for prop. in 'x' (longitude) dir.
!      MAPY2     I.A.  Public    Idem in y' (latitude) direction.
!      MAPXY     I.A.  Public
!      MAPAXY    I.A.  Public    List of active points used in W3QCK1.
!      MAPCXY    I.A.  Public    List of central points used in avg.
!      MAPTH2    I.A.  Public    Like MAPX2 for refraction (rotated
!                                and shifted, see W3KTP3). Like MAPAXY.
!      MAPWN2    I.A.  Public    Like MAPX2 for wavenumber shift.
!      MAPTRN    L.A.  Public    Map to block out GSE mitigation in
!                                proper grid points.
!
!     Nonlinear interactions ( !/NL1 ) :
!
!      NFR       Int.  Public   Nuber of frequencies ( NFR = NK )
!      NFRHGH    Int.  Public   Auxiliary frequency counter.
!      NFRCHG    Int.  Public   Id.
!      NSPECX-Y  Int.  Public   Auxiliary spectral counter.
!      IPnn      I.A.  Public   Spectral address for Snl.
!      IMnn      I.A.  Public   Id.
!      ICnn      I.A.  Public   Id.
!      DALn      Real  Public   Lambda dependend weight factors.
!      AWGn      Real  Public   Interpolation weights for Snl.
!      SWGn      Real  Public   Interpolation weights for diag. term.
!      AF11      R.A.  Public   Scaling array (f**11)
!      NLINIT    Log.  Public   Flag for initialization.
!
!     MPP / MPI variables :
!
!      IAPPRO    I.A.  Public   Processor numbers for propagation calc.
!                               for each spectral component.
!      MPI_COMM_WAVE
!                Int.  Public   Communicator used in the wave model.
!      MPI_COMM_WCMP
!                Int.  Public   Idem, computational proc. only.
!      WW3_FIELD_VEC, WW3_SPEC_VEC
!                Int.  Public   MPI derived vecor types.
!      NRQSG1    Int.  Public   Number of handles in IRQSG1.
!      NRQSG2    Int.  Public   Number of handles in IRQSG2.
!      IBFLOC    Int.  Public   Present active buffer number.
!      ISPLOC    Int.  Public   Corresponding local spectral bin number
!                               (1,NSPLOC,1).
!      NSPLOC    Int.  Public   Total number of spectral bins for which
!                               prop. is performed on present CPU.
!      BSTAT     I.A.  Public   Status of buffer (size MPIBUF):
!                                 0: Inactive.
!                                 1: A --> STORE (active or finished).
!                                 2: STORE --> A (active or finished).
!      BISPL     I.A.  Public   Local spectral bin number for buffer
!                               (size MPIBUF).
!      IRQSG1    I.A.  Public   MPI request handles for scatters and
!                               gathers to A() (persistent).
!      IRQSG2    I.A.  Public   MPI request handles for gathers and
!                               scatters to STORE (persistent).
!    G/SSTORE    R.A.  Public   Communication buffer (NSEA,MPIBUF).
!      SPPNT     R.A.  Public   Point output buffer.
!
!     Other:
!
!      ITIME     Int.  Public   Discrete time step counter.
!      IPASS     Int.  Public   Pass counter for log file.
!      IDLAST    Int.  Public   Last day ID for log file.
!      NSEALM    Int.  Public   Maximum number of local sea points.
!      ALPHA     R.A.  Public   Phillips' alpha.
!      FLCOLD    Log.  Public   Flag for 'cold start' of model.
!      FLIWND    Log.  Public   Flag for initialization of model
!                               based on wind.
!      AINIT     Log.  Public   Flag for array initialization.
!      FL_ALL    Log.  Public   Flag for all/partial  initialization.
!     ----------------------------------------------------------------
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3NAUX    Subr. Public   Set number of grids/models.
!      W3DIMA    Subr. Public   Set dimensions of arrays.
!      W3DMNL    Subr. Public   Set dimensions of arrays.   ( !/NL1 )
!      W3SETA    Subr. Public   Point to selected grid / model.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SETG    Subr. W3GDATMD Point to proper model grid.
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!      EXTCDE    Subr. W3SERVMD Abort program with exit code.
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!     - The number of grids is taken from W3GDATMD, and needs to be
!       set first with W3DIMG.
!
!  6. Switches :
!
!     !/SHRD, !/DIST, !/MPI
!              Shared / distributed memory model
!
!     !/PRn    Propagation scheme selection.
!
!     !/S      Enable subroutine tracing.
!     !/T      Enable test output
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
      PUBLIC
!/
!/ Conventional declarations
!/
      INTEGER                 :: NADATA = -1, IADATA = -1
      INTEGER, PARAMETER      :: MPIBUF = 6
!/
!/ Data structure WADAT
!/
      TYPE WADAT
!
        REAL, POINTER         :: CG(:,:), WN(:,:)
!
        REAL, POINTER         :: CA0(:), CAI(:), CD0(:), CDI(:),      &
                                 UA0(:), UAI(:), UD0(:), UDI(:),      &
                                 AS0(:), ASI(:), ATRNX(:,:), ATRNY(:,:)
!
        REAL, POINTER         :: DW(:), UA(:), UD(:), U10(:), U10D(:),&
                                 AS(:), CX(:), CY(:), EMN(:), FMN(:), &
                                 WNM(:), AMX(:), CDS(:), Z0S(:),      &
                                 HS(:), WLM(:), TMN(:), THM(:),       &
                                 THS(:), FP0(:), THP0(:), FP1(:),     &
                                 THP1(:), DTDYN(:), FCUT(:),          &
                                 ABA(:), ABD(:), UBA(:), UBD(:),      &
                                 SXX(:), SYY(:), SXY(:), USERO(:,:)
        REAL, POINTER         :: PHS(:,:), PTP(:,:), PLP(:,:),        &
                                 PTH(:,:), PSI(:,:), PWS(:,:),        &
                                 PWST(:), PNR(:)
        REAL, POINTER         :: DDDX(:,:), DDDY(:,:), DCXDX(:,:),    &
                                 DCYDX(:,:), DCXDY(:,:), DCYDY(:,:)
!
        INTEGER               :: NMX0, NMX1, NMX2, NMY0, NMY1, NMY2,  &
                                 NACT, NCENT
        INTEGER, POINTER      :: MAPX2(:), MAPY2(:), MAPAXY(:),       &
                                 MAPCXY(:), MAPTH2(:), MAPWN2(:)
        LOGICAL, POINTER      :: MAPTRN(:)
!
        INTEGER               :: NFR, NFRHGH, NFRCHG, NSPECX, NSPECY
        INTEGER, POINTER      :: IP11(:), IP12(:), IP13(:), IP14(:),  &
                                 IM11(:), IM12(:), IM13(:), IM14(:),  &
                                 IP21(:), IP22(:), IP23(:), IP24(:),  &
                                 IM21(:), IM22(:), IM23(:), IM24(:),  &
                                 IC11(:), IC12(:), IC21(:), IC22(:),  &
                                 IC31(:), IC32(:), IC41(:), IC42(:),  &
                                 IC51(:), IC52(:), IC61(:), IC62(:),  &
                                 IC71(:), IC72(:), IC81(:), IC82(:)
        REAL                  :: DAL1, DAL2, DAL3,                    &
                                 AWG1, AWG2, AWG3, AWG4, AWG5, AWG6,  &
                                 AWG7, AWG8, SWG1, SWG2, SWG3, SWG4,  &
                                 SWG5, SWG6, SWG7, SWG8
        REAL, POINTER         :: AF11(:)
        LOGICAL               :: NLINIT
!
        INTEGER, POINTER      :: IAPPRO(:)
        INTEGER               :: MPI_COMM_WAVE, MPI_COMM_WCMP,        &
                                 WW3_FIELD_VEC, WW3_SPEC_VEC,         &
                                 NRQSG1, NRQSG2, IBFLOC, ISPLOC,      &
                                 NSPLOC
        INTEGER               :: BSTAT(MPIBUF), BISPL(MPIBUF)
        INTEGER, POINTER      :: IRQSG1(:,:), IRQSG2(:,:)
        REAL, POINTER         :: GSTORE(:,:), SSTORE(:,:)
        REAL, POINTER         :: SPPNT(:,:,:)
!
        INTEGER               :: ITIME, IPASS, IDLAST, NSEALM
        REAL, POINTER         :: ALPHA(:,:)
        LOGICAL               :: AINIT, FL_ALL, FLCOLD, FLIWND
!
      END TYPE WADAT
!/
!/ Data storage
!/
      TYPE(WADAT), TARGET, ALLOCATABLE :: WADATS(:)
!/
!/ Data aliasses for structure WADAT(S)
!/
      REAL, POINTER           :: CG(:,:), WN(:,:)
!
      REAL, POINTER           :: CA0(:), CAI(:), CD0(:), CDI(:),      &
                                 UA0(:), UAI(:), UD0(:), UDI(:),      &
                                 AS0(:), ASI(:), ATRNX(:,:), ATRNY(:,:)
!
      REAL, POINTER           :: DW(:), UA(:), UD(:), U10(:), U10D(:),&
                                 AS(:), CX(:), CY(:), EMN(:), FMN(:), &
                                 WNM(:), AMX(:), CDS(:), Z0S(:),      &
                                 HS(:), WLM(:), TMN(:), THM(:),       &
                                 THS(:), FP0(:), THP0(:), FP1(:),     &
                                 THP1(:), DTDYN(:), FCUT(:),          &
                                 ABA(:), ABD(:), UBA(:), UBD(:),      &
                                 SXX(:), SYY(:), SXY(:), USERO(:,:)
      REAL, POINTER           :: PHS(:,:), PTP(:,:), PLP(:,:),        &
                                 PTH(:,:), PSI(:,:), PWS(:,:),        &
                                 PWST(:), PNR(:)
      REAL, POINTER           :: DDDX(:,:), DDDY(:,:), DCXDX(:,:),    &
                                 DCYDX(:,:), DCXDY(:,:), DCYDY(:,:)
!
      INTEGER, POINTER        :: NMX0, NMX1, NMX2, NMY0, NMY1, NMY2,  &
                                 NACT, NCENT
      INTEGER, POINTER        :: MAPX2(:), MAPY2(:), MAPAXY(:),       &
                                 MAPCXY(:), MAPTH2(:), MAPWN2(:)
      LOGICAL, POINTER        :: MAPTRN(:)
!
      INTEGER, POINTER        :: NFR, NFRHGH, NFRCHG, NSPECX, NSPECY
      INTEGER, POINTER        :: IP11(:), IP12(:), IP13(:), IP14(:),  &
                                 IM11(:), IM12(:), IM13(:), IM14(:),  &
                                 IP21(:), IP22(:), IP23(:), IP24(:),  &
                                 IM21(:), IM22(:), IM23(:), IM24(:),  &
                                 IC11(:), IC12(:), IC21(:), IC22(:),  &
                                 IC31(:), IC32(:), IC41(:), IC42(:),  &
                                 IC51(:), IC52(:), IC61(:), IC62(:),  &
                                 IC71(:), IC72(:), IC81(:), IC82(:)
      REAL, POINTER           :: DAL1, DAL2, DAL3,                    &
                                 AWG1, AWG2, AWG3, AWG4, AWG5, AWG6,  &
                                 AWG7, AWG8, SWG1, SWG2, SWG3, SWG4,  &
                                 SWG5, SWG6, SWG7, SWG8
      REAL, POINTER           :: AF11(:)
      LOGICAL, POINTER        :: NLINIT
!
      INTEGER, POINTER        :: IAPPRO(:)
      INTEGER, POINTER        :: MPI_COMM_WAVE, MPI_COMM_WCMP,        &
                                 WW3_FIELD_VEC, WW3_SPEC_VEC,         &
                                 NRQSG1, NRQSG2, IBFLOC, ISPLOC,      &
                                 NSPLOC
      INTEGER, POINTER        :: BSTAT(:), BISPL(:)
      INTEGER, POINTER        :: IRQSG1(:,:), IRQSG2(:,:)
      REAL, POINTER           :: GSTORE(:,:), SSTORE(:,:)
      REAL, POINTER           :: SPPNT(:,:,:)
!
      INTEGER, POINTER        :: ITIME, IPASS, IDLAST, NSEALM
      REAL, POINTER           :: ALPHA(:,:)
      LOGICAL, POINTER        :: AINIT, FL_ALL, FLCOLD, FLIWND
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3NAUX ( NDSE, NDST )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         14-Dec-2004 !
!/                  +-----------------------------------+
!/
!/    14-Dec-2004 : Origination.                        ( version 3.06 )
!/    04-Oct-2006 : Add filter to array pointers.       ( version 3.10 )
!/
!  1. Purpose :
!
!     Set up the number of grids to be used.
!
!  2. Method :
!
!     Use data stored in NGRIDS in W3GDATMD.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NDSE    Int.   I   Error output unit number.
!       NDST    Int.   I   Test output unit number.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     See module documentation.
!
!  5. Called by :
!
!     Any program that uses this grid structure.
!
!  6. Error messages :
!
!     - Error checks on previous setting of variable NGRIDS.
!
!  7. Remarks :
!
!  8. Structure :
!
!  9. Switches :
!
!     !/S    Enable subroutine tracing.
!     !/T    Enable test output
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3GDATMD, ONLY: NGRIDS
      USE W3SERVMD, ONLY: EXTCDE
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: NDSE, NDST
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: I
!/
!
! -------------------------------------------------------------------- /
! 1.  Test input and module status
!
      IF ( NGRIDS .EQ. -1 ) THEN
          WRITE (NDSE,1001) NGRIDS
          CALL EXTCDE (1)
        END IF
!
! -------------------------------------------------------------------- /
! 2.  Set variable and allocate arrays
!
      ALLOCATE ( WADATS(NGRIDS) )
      NADATA = NGRIDS
!
! -------------------------------------------------------------------- /
! 3.  Initialize parameters
!
      DO I=1, NGRIDS
        WADATS(I)%ITIME  = 0
        WADATS(I)%IPASS  = 0
        WADATS(I)%IDLAST = 0
        WADATS(I)%AINIT  = .FALSE.
        WADATS(I)%FL_ALL = .FALSE.
        WADATS(I)%NLINIT  = .FALSE.
        END DO
!
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** ERROR W3NAUX : NGRIDS NOT YET SET *** '/         &
               '                    NGRIDS = ',I10/                   &
               '                    RUN W3NMOD FIRST'/)
!
!/
!/ End of W3NAUX ----------------------------------------------------- /
!/
      END SUBROUTINE W3NAUX
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3DIMA  ( IMOD, NDSE, NDST, D_ONLY )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         22-Feb-2008 !
!/                  +-----------------------------------+
!/
!/    28-Dec-2004 : Origination.                        ( version 3.06 )
!/    20-Jul-2005 : Adding output fields.               ( version 3.07 )
!/    04-Oct-2006 : Add filter to array pointers.       ( version 3.10 )
!/    28-Mar-2007 : Add partitioned data arrays.        ( version 3.11 )
!/                  Add aditional undefined arrays.
!/    22-Feb-2008 ; Modify MAPTH2 declaration.          ( version 3.14 )
!/
!  1. Purpose :
!
!     Initialize an individual data grid at the proper dimensions.
!
!  2. Method :
!
!     Allocate directly into the structure array. Note that
!     this cannot be done through the pointer alias!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.   I   Model number to point to.
!       NDSE    Int.   I   Error output unit number.
!       NDST    Int.   I   Test output unit number.
!       D_ONLY  L.O.   I   FLag for initializing data arrays only.
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
!      W3IOGO    Subr. W3IOGOMD Grid output IO routine.
!      WW3_SHEL  Prog.   N/A    Wave model driver.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     - Check on input parameters.
!     - Check on previous allocation.
!
!  7. Remarks :
!
!     - W3SETA needs to be called after allocation to point to
!       proper allocated arrays.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/SHRD, !/DIST
!              Shared / distributed memory model
!
!     !/PRn    Propagation scheme selection.
!
!     !/S      Enable subroutine tracing.
!     !/T    Enable test output
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3GDATMD, ONLY: NGRIDS, IGRID, W3SETG, NK, NX, NY, NSEA,    &
                          NSEAL, NSPEC, NTH
      USE W3ODATMD, ONLY: NAPROC, NOSWLL, NOEXTR, UNDEF
      USE W3IDATMD, ONLY: FLCUR, FLWIND
      USE W3SERVMD, ONLY: EXTCDE
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)           :: IMOD, NDSE, NDST
      LOGICAL, INTENT(IN), OPTIONAL :: D_ONLY
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: JGRID, NXXX
!/
!
! -------------------------------------------------------------------- /
! 1.  Test input and module status
!
      IF ( PRESENT(D_ONLY) ) THEN
          FL_ALL = .NOT. D_ONLY
        ELSE
          FL_ALL = .TRUE.
        END IF
!
      IF ( NGRIDS .EQ. -1 ) THEN
          WRITE (NDSE,1001)
          CALL EXTCDE (1)
        END IF
!
      IF ( IMOD.LT.1 .OR. IMOD.GT.NADATA ) THEN
          WRITE (NDSE,1002) IMOD, NADATA
          CALL EXTCDE (2)
        END IF
!
      IF ( WADATS(IMOD)%AINIT ) THEN
          WRITE (NDSE,1003)
          CALL EXTCDE (3)
        END IF
!
      JGRID  = IGRID
      IF ( JGRID .NE. IMOD ) CALL W3SETG ( IMOD, NDSE, NDST )
!
! -------------------------------------------------------------------- /
! 2.  Allocate arrays
!     Call W3SETA to assure of pointes FLCUR an FLWND
!
      CALL W3SETA ( IMOD, NDSE, NDST )
!
      NSEALM = 1 + (NSEA-1)/NAPROC
      NXXX   = NSEALM * NAPROC
!
      ALLOCATE ( WADATS(IMOD)%DW(0:NSEA) , WADATS(IMOD)%UA(0:NSEA)  , &
                 WADATS(IMOD)%UD(0:NSEA) , WADATS(IMOD)%U10(NSEA)   , &
                 WADATS(IMOD)%U10D(NSEA) , WADATS(IMOD)%AS(0:NSEA)  , &
                 WADATS(IMOD)%CX(0:NSEA) , WADATS(IMOD)%CY(0:NSEA)  , &
                 WADATS(IMOD)%EMN(NSEA)  , WADATS(IMOD)%FMN(NSEA)   , &
                 WADATS(IMOD)%WNM(NSEA)  , WADATS(IMOD)%AMX(NSEA)   , &
                 WADATS(IMOD)%CDS(NSEA)  , WADATS(IMOD)%Z0S(NSEA)   , &
                 WADATS(IMOD)%HS(NXXX)   , WADATS(IMOD)%WLM(NXXX)   , &
                 WADATS(IMOD)%TMN(NXXX)  , WADATS(IMOD)%THM(NXXX)   , &
                 WADATS(IMOD)%THS(NXXX)  , WADATS(IMOD)%FP0(NXXX)   , &
                 WADATS(IMOD)%THP0(NXXX) , WADATS(IMOD)%FP1(NXXX)   , &
                 WADATS(IMOD)%THP1(NXXX) , WADATS(IMOD)%DTDYN(NXXX) , &
                 WADATS(IMOD)%FCUT(NXXX) , WADATS(IMOD)%ABA(NXXX)   , &
                 WADATS(IMOD)%ABD(NXXX)  , WADATS(IMOD)%UBA(NXXX)   , &
                 WADATS(IMOD)%UBD(NXXX)  , WADATS(IMOD)%SXX(NXXX)   , &
                 WADATS(IMOD)%SYY(NXXX)  , WADATS(IMOD)%SXY(NXXX)   , &
                 WADATS(IMOD)%USERO(NXXX,NOEXTR) )
!
      WADATS(IMOD)%USERO = UNDEF
!
      ALLOCATE ( WADATS(IMOD)%PHS(NXXX,0:NOSWLL),                     &
                 WADATS(IMOD)%PTP(NXXX,0:NOSWLL),                     &
                 WADATS(IMOD)%PLP(NXXX,0:NOSWLL),                     &
                 WADATS(IMOD)%PTH(NXXX,0:NOSWLL),                     &
                 WADATS(IMOD)%PSI(NXXX,0:NOSWLL),                     &
                 WADATS(IMOD)%PWS(NXXX,0:NOSWLL),                     &
                 WADATS(IMOD)%PWST(NXXX), WADATS(IMOD)%PNR(NXXX) )
!
      IF ( FL_ALL ) THEN
!
          ALLOCATE ( WADATS(IMOD)%CG(0:NK+1,0:NSEA) ,                 &
                     WADATS(IMOD)%WN(0:NK+1,0:NSEA) )
!
          IF ( FLCUR  ) ALLOCATE ( WADATS(IMOD)%CA0(NSEA) ,           &
                                   WADATS(IMOD)%CAI(NSEA) ,           &
                                   WADATS(IMOD)%CD0(NSEA) ,           &
                                   WADATS(IMOD)%CDI(NSEA) )
!
          IF ( FLWIND ) ALLOCATE ( WADATS(IMOD)%UA0(NSEA) ,           &
                                   WADATS(IMOD)%UAI(NSEA) ,           &
                                   WADATS(IMOD)%UD0(NSEA) ,           &
                                   WADATS(IMOD)%UDI(NSEA) ,           &
                                   WADATS(IMOD)%AS0(NSEA) ,           &
                                   WADATS(IMOD)%ASI(NSEA) )
!
          ALLOCATE ( WADATS(IMOD)%ATRNX(NY*NX,-1:1) ,                 &
                     WADATS(IMOD)%ATRNY(NY*NX,-1:1) )
!
          ALLOCATE ( WADATS(IMOD)%DDDX(NY,NX)  ,                      &
                     WADATS(IMOD)%DDDY(NY,NX)  ,                      &
                     WADATS(IMOD)%DCXDX(NY,NX) ,                      &
                     WADATS(IMOD)%DCYDX(NY,NX) ,                      &
                     WADATS(IMOD)%DCXDY(NY,NX) ,                      &
                     WADATS(IMOD)%DCYDY(NY,NX) )
!
          ALLOCATE ( WADATS(IMOD)%ALPHA(NK,NSEAL) )
!
          ALLOCATE ( WADATS(IMOD)%MAPX2(NY*NX)       ,           &
                     WADATS(IMOD)%MAPY2(NY*NX)       ,           &
                     WADATS(IMOD)%MAPAXY(NY*NX)      ,           &
                     WADATS(IMOD)%MAPCXY(NSEA)       ,           &
                     WADATS(IMOD)%MAPTH2((NK+2)*NTH) ,           &
                     WADATS(IMOD)%MAPWN2(NSPEC+NTH)  ,           &
                     WADATS(IMOD)%MAPTRN(NY*NX) )
          WADATS(IMOD)%MAPTH2 = 0
!
          ALLOCATE ( WADATS(IMOD)%IAPPRO(NSPEC) ,                     &
                     WADATS(IMOD)%SPPNT(NTH,NK,4) )
!
        END IF
!
      WADATS(IMOD)%AINIT  = .TRUE.
!
! -------------------------------------------------------------------- /
! 3.  Point to allocated arrays
!
      CALL W3SETA ( IMOD, NDSE, NDST )
!
! -------------------------------------------------------------------- /
! 4.  Update counters in grid
!
! -------------------------------------------------------------------- /
! 5.  Restore previous grid setting if necessary
!
      IF ( JGRID .NE. IMOD ) CALL W3SETG ( JGRID, NDSE, NDST )
!
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** ERROR W3DIMA : GRIDS NOT INITIALIZED *** '/      &
               '                    RUN W3NMOD FIRST '/)
 1002 FORMAT (/' *** ERROR W3DIMA : ILLEGAL MODEL NUMBER *** '/       &
               '                    IMOD   = ',I10/                   &
               '                    NADATA = ',I10/)
 1003 FORMAT (/' *** ERROR W3DIMA : ARRAY(S) ALREADY ALLOCATED *** ')
!
!/
!/ End of W3DIMA ----------------------------------------------------- /
!/
      END SUBROUTINE W3DIMA
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3DMNL  ( IMOD, NDSE, NDST, NSP, NSPX )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         04-Oct-2006 !
!/                  +-----------------------------------+
!/
!/    24-Dec-2004 : Origination.                        ( version 3.06 )
!/    04-Oct-2006 : Add filter to array pointers.       ( version 3.10 )
!/
!  1. Purpose :
!
!     Initialize an individual data grid at the proper dimensions (DIA).
!
!  2. Method :
!
!     Allocate directly into the structure array. Note that
!     this cannot be done through the pointer alias!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.   I   Model number to point to.
!       NDSE    Int.   I   Error output unit number.
!       NDST    Int.   I   Test output unit number.
!       NSP(X)  Int.   I   Array dimensions.
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
!      INSNL1    Subr. W3SNL1MD Traditional DIA approach to Snl.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     - Check on input parameters.
!     - Check on previous allocation.
!
!  7. Remarks :
!
!     - W3SETA needs to be called after allocation to point to
!       proper allocated arrays.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S      Enable subroutine tracing.
!     !/T    Enable test output
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3GDATMD, ONLY: NGRIDS, IGRID, NK, NX, NY, NSEA, NSEAL,     &
                         NSPEC, NTH
      USE W3ODATMD, ONLY: NAPROC
      USE W3IDATMD, ONLY: FLCUR, FLWIND
      USE W3SERVMD, ONLY: EXTCDE
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)           :: IMOD, NDSE, NDST, NSP, NSPX
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
!/
!
! -------------------------------------------------------------------- /
! 1.  Test input and module status
!
      IF ( NGRIDS .EQ. -1 ) THEN
          WRITE (NDSE,1001)
          CALL EXTCDE (1)
        END IF
!
      IF ( IMOD.LT.1 .OR. IMOD.GT.NADATA ) THEN
          WRITE (NDSE,1002) IMOD, NADATA
          CALL EXTCDE (2)
        END IF
!
      IF ( WADATS(IMOD)%NLINIT ) THEN
          WRITE (NDSE,1003)
          CALL EXTCDE (3)
        END IF
!
! -------------------------------------------------------------------- /
! 2.  Allocate arrays
!
      ALLOCATE ( WADATS(IMOD)%IP11(NSPX),        &
                 WADATS(IMOD)%IP12(NSPX),        &
                 WADATS(IMOD)%IP13(NSPX),        &
                 WADATS(IMOD)%IP14(NSPX),        &
                 WADATS(IMOD)%IM11(NSPX),        &
                 WADATS(IMOD)%IM12(NSPX),        &
                 WADATS(IMOD)%IM13(NSPX),        &
                 WADATS(IMOD)%IM14(NSPX),        &
                 WADATS(IMOD)%IP21(NSPX),        &
                 WADATS(IMOD)%IP22(NSPX),        &
                 WADATS(IMOD)%IP23(NSPX),        &
                 WADATS(IMOD)%IP24(NSPX),        &
                 WADATS(IMOD)%IM21(NSPX),        &
                 WADATS(IMOD)%IM22(NSPX),        &
                 WADATS(IMOD)%IM23(NSPX),        &
                 WADATS(IMOD)%IM24(NSPX),        &
                 WADATS(IMOD)%IC11(NSP) ,        &
                 WADATS(IMOD)%IC12(NSP) ,        &
                 WADATS(IMOD)%IC21(NSP) ,        &
                 WADATS(IMOD)%IC22(NSP) ,        &
                 WADATS(IMOD)%IC31(NSP) ,        &
                 WADATS(IMOD)%IC32(NSP) ,        &
                 WADATS(IMOD)%IC41(NSP) ,        &
                 WADATS(IMOD)%IC42(NSP) ,        &
                 WADATS(IMOD)%IC51(NSP) ,        &
                 WADATS(IMOD)%IC52(NSP) ,        &
                 WADATS(IMOD)%IC61(NSP) ,        &
                 WADATS(IMOD)%IC62(NSP) ,        &
                 WADATS(IMOD)%IC71(NSP) ,        &
                 WADATS(IMOD)%IC72(NSP) ,        &
                 WADATS(IMOD)%IC81(NSP) ,        &
                 WADATS(IMOD)%IC82(NSP) ,        &
                 WADATS(IMOD)%AF11(NSPX) )
!
      WADATS(IMOD)%NLINIT = .TRUE.
!
! -------------------------------------------------------------------- /
! 3.  Point to allocated arrays
!
      CALL W3SETA ( IMOD, NDSE, NDST )
!
! -------------------------------------------------------------------- /
! 4.  Update counters in grid
!
      NSPECX = NSPX
!
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** ERROR W3DMNL : GRIDS NOT INITIALIZED *** '/      &
               '                    RUN W3NMOD FIRST '/)
 1002 FORMAT (/' *** ERROR W3DMNL : ILLEGAL MODEL NUMBER *** '/       &
               '                    IMOD   = ',I10/                   &
               '                    NADATA = ',I10/)
 1003 FORMAT (/' *** ERROR W3DMNL : ARRAY(S) ALREADY ALLOCATED *** ')
!
!/
!/ End of W3DMNL ----------------------------------------------------- /
!/
      END SUBROUTINE W3DMNL
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SETA ( IMOD, NDSE, NDST )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         28_Mar-2007 |
!/                  +-----------------------------------+
!/
!/    28-Dec-2004 : Origination.                        ( version 3.06 )
!/    04-May-2005 : Adding MPI_COMM_WAVE.               ( version 3.07 )
!/    20-Jul-2005 : Adding output fields.               ( version 3.07 )
!/    09-Nov-2005 : Removing soft boundary option.      ( version 3.08 )
!/    13-Jun-2006 : Splitting STORE in G/SSTORE.        ( version 3.09 )
!/    04-Oct-2006 : Add filter to array pointers.       ( version 3.10 )
!/    28_Mar-2007 : Add partitioned data arrays.        ( version 3.11 )
!/                  Add aditional undefined arrays.
!/
!  1. Purpose :
!
!     Select one of the WAVEWATCH III grids / models.
!
!  2. Method :
!
!     Point pointers to the proper variables in the proper element of
!     the GRIDS array.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.   I   Model number to point to.
!       NDSE    Int.   I   Error output unit number.
!       NDST    Int.   I   Test output unit number.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     See module documentation below.
!
!  5. Called by :
!
!     Many subroutines in the WAVEWATCH system.
!
!  6. Error messages :
!
!     Checks on parameter list IMOD.
!
!  7. Remarks :
!
!  8. Structure :
!
!  9. Switches :
!
!     !/MPI  Paralllel model environment.
!
!     !/PRn    Propagation scheme selection.
!
!     !/S    Enable subroutine tracing.
!     !/T    Enable test output
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE W3IDATMD, ONLY: INPUTS
!
      USE W3SERVMD, ONLY: EXTCDE
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: IMOD, NDSE, NDST
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
!/
!
! -------------------------------------------------------------------- /
! 1.  Test input and module status
!
      IF ( NADATA .EQ. -1 ) THEN
          WRITE (NDSE,1001)
          CALL EXTCDE (1)
        END IF
!
      IF ( IMOD.LT.1 .OR. IMOD.GT.NADATA ) THEN
          WRITE (NDSE,1002) IMOD, NADATA
          CALL EXTCDE (2)
        END IF
!
! -------------------------------------------------------------------- /
! 2.  Set model numbers
!
      IADATA = IMOD
!
! -------------------------------------------------------------------- /
! 3.  Set pointers
!
      ITIME  => WADATS(IMOD)%ITIME
      IPASS  => WADATS(IMOD)%IPASS
      IDLAST => WADATS(IMOD)%IDLAST
      NSEALM => WADATS(IMOD)%NSEALM
      FLCOLD => WADATS(IMOD)%FLCOLD
      FLIWND => WADATS(IMOD)%FLIWND
      AINIT  => WADATS(IMOD)%AINIT
      FL_ALL => WADATS(IMOD)%FL_ALL
!
      NMX0   => WADATS(IMOD)%NMX0
      NMX1   => WADATS(IMOD)%NMX1
      NMX2   => WADATS(IMOD)%NMX2
      NMY0   => WADATS(IMOD)%NMY0
      NMY1   => WADATS(IMOD)%NMY1
      NMY2   => WADATS(IMOD)%NMY2
      NACT   => WADATS(IMOD)%NACT
      NCENT  => WADATS(IMOD)%NCENT
!
      NFR    => WADATS(IMOD)%NFR
      NFRHGH => WADATS(IMOD)%NFRHGH
      NFRCHG => WADATS(IMOD)%NFRCHG
      NSPECX => WADATS(IMOD)%NSPECX
      NSPECY => WADATS(IMOD)%NSPECY
      DAL1   => WADATS(IMOD)%DAL1
      DAL2   => WADATS(IMOD)%DAL2
      DAL3   => WADATS(IMOD)%DAL3
      AWG1   => WADATS(IMOD)%AWG1
      AWG2   => WADATS(IMOD)%AWG2
      AWG3   => WADATS(IMOD)%AWG3
      AWG4   => WADATS(IMOD)%AWG4
      AWG5   => WADATS(IMOD)%AWG5
      AWG6   => WADATS(IMOD)%AWG6
      AWG7   => WADATS(IMOD)%AWG7
      AWG8   => WADATS(IMOD)%AWG8
      SWG1   => WADATS(IMOD)%SWG1
      SWG2   => WADATS(IMOD)%SWG2
      SWG3   => WADATS(IMOD)%SWG3
      SWG4   => WADATS(IMOD)%SWG4
      SWG5   => WADATS(IMOD)%SWG5
      SWG6   => WADATS(IMOD)%SWG6
      SWG7   => WADATS(IMOD)%SWG7
      SWG8   => WADATS(IMOD)%SWG8
      NLINIT => WADATS(IMOD)%NLINIT
!
      MPI_COMM_WAVE => WADATS(IMOD)%MPI_COMM_WAVE
      MPI_COMM_WCMP => WADATS(IMOD)%MPI_COMM_WCMP
      WW3_FIELD_VEC => WADATS(IMOD)%WW3_FIELD_VEC
      WW3_SPEC_VEC => WADATS(IMOD)%WW3_SPEC_VEC
      NRQSG1 => WADATS(IMOD)%NRQSG1
      NRQSG2 => WADATS(IMOD)%NRQSG2
      IBFLOC => WADATS(IMOD)%IBFLOC
      ISPLOC => WADATS(IMOD)%ISPLOC
      NSPLOC => WADATS(IMOD)%NSPLOC
      BSTAT  => WADATS(IMOD)%BSTAT
      BISPL  => WADATS(IMOD)%BISPL
!
      IF ( AINIT ) THEN
!
          DW     => WADATS(IMOD)%DW
          UA     => WADATS(IMOD)%UA
          UD     => WADATS(IMOD)%UD
          U10    => WADATS(IMOD)%U10
          U10D   => WADATS(IMOD)%U10D
          AS     => WADATS(IMOD)%AS
          CX     => WADATS(IMOD)%CX
          CY     => WADATS(IMOD)%CY
          EMN    => WADATS(IMOD)%EMN
          FMN    => WADATS(IMOD)%FMN
          WNM    => WADATS(IMOD)%WNM
          AMX    => WADATS(IMOD)%AMX
          CDS    => WADATS(IMOD)%CDS
          Z0S    => WADATS(IMOD)%Z0S
          HS     => WADATS(IMOD)%HS
          WLM    => WADATS(IMOD)%WLM
          TMN    => WADATS(IMOD)%TMN
          THM    => WADATS(IMOD)%THM
          THS    => WADATS(IMOD)%THS
          FP0    => WADATS(IMOD)%FP0
          THP0   => WADATS(IMOD)%THP0
          FP1    => WADATS(IMOD)%FP1
          THP1   => WADATS(IMOD)%THP1
          PHS    => WADATS(IMOD)%PHS
          PTP    => WADATS(IMOD)%PTP
          PLP    => WADATS(IMOD)%PLP
          PTH    => WADATS(IMOD)%PTH
          PSI    => WADATS(IMOD)%PSI
          PWS    => WADATS(IMOD)%PWS
          PWST   => WADATS(IMOD)%PWST
          PNR    => WADATS(IMOD)%PNR
          DTDYN  => WADATS(IMOD)%DTDYN
          FCUT   => WADATS(IMOD)%FCUT
          ABA    => WADATS(IMOD)%ABA
          ABD    => WADATS(IMOD)%ABD
          UBA    => WADATS(IMOD)%UBA
          UBD    => WADATS(IMOD)%UBD
          SXX    => WADATS(IMOD)%SXX
          SYY    => WADATS(IMOD)%SYY
          SXY    => WADATS(IMOD)%SXY
          USERO  => WADATS(IMOD)%USERO
!
          IF ( FL_ALL ) THEN
!
              CG     => WADATS(IMOD)%CG
              WN     => WADATS(IMOD)%WN
!
              ATRNX  => WADATS(IMOD)%ATRNX
              ATRNY  => WADATS(IMOD)%ATRNY
!
              DDDX   => WADATS(IMOD)%DDDX
              DDDY   => WADATS(IMOD)%DDDY
              DCXDX  => WADATS(IMOD)%DCXDX
              DCYDX  => WADATS(IMOD)%DCYDX
              DCXDY  => WADATS(IMOD)%DCXDY
              DCYDY  => WADATS(IMOD)%DCYDY
!
              ALPHA  => WADATS(IMOD)%ALPHA
!
              IF ( INPUTS(IMOD)%FLAGS(2) ) THEN
                  CA0    => WADATS(IMOD)%CA0
                  CAI    => WADATS(IMOD)%CAI
                  CD0    => WADATS(IMOD)%CD0
                  CDI    => WADATS(IMOD)%CDI
                END IF
!
              IF ( INPUTS(IMOD)%FLAGS(3) ) THEN
                  UA0    => WADATS(IMOD)%UA0
                  UAI    => WADATS(IMOD)%UAI
                  UD0    => WADATS(IMOD)%UD0
                  UDI    => WADATS(IMOD)%UDI
                  AS0    => WADATS(IMOD)%AS0
                  ASI    => WADATS(IMOD)%ASI
                END IF
!
              MAPX2  => WADATS(IMOD)%MAPX2
              MAPY2  => WADATS(IMOD)%MAPY2
              MAPAXY => WADATS(IMOD)%MAPAXY
              MAPCXY => WADATS(IMOD)%MAPCXY
              MAPTH2 => WADATS(IMOD)%MAPTH2
              MAPWN2 => WADATS(IMOD)%MAPWN2
              MAPTRN => WADATS(IMOD)%MAPTRN
!
              IAPPRO => WADATS(IMOD)%IAPPRO
              SPPNT  => WADATS(IMOD)%SPPNT
!
            END IF
!
        END IF
!
      IF ( NLINIT ) THEN
          IP11   => WADATS(IMOD)%IP11
          IP12   => WADATS(IMOD)%IP12
          IP13   => WADATS(IMOD)%IP13
          IP14   => WADATS(IMOD)%IP14
          IM11   => WADATS(IMOD)%IM11
          IM12   => WADATS(IMOD)%IM12
          IM13   => WADATS(IMOD)%IM13
          IM14   => WADATS(IMOD)%IM14
          IP21   => WADATS(IMOD)%IP21
          IP22   => WADATS(IMOD)%IP22
          IP23   => WADATS(IMOD)%IP23
          IP24   => WADATS(IMOD)%IP24
          IM21   => WADATS(IMOD)%IM21
          IM22   => WADATS(IMOD)%IM22
          IM23   => WADATS(IMOD)%IM23
          IM24   => WADATS(IMOD)%IM24
          IC11   => WADATS(IMOD)%IC11
          IC12   => WADATS(IMOD)%IC12
          IC21   => WADATS(IMOD)%IC21
          IC22   => WADATS(IMOD)%IC22
          IC31   => WADATS(IMOD)%IC31
          IC32   => WADATS(IMOD)%IC32
          IC41   => WADATS(IMOD)%IC41
          IC42   => WADATS(IMOD)%IC42
          IC51   => WADATS(IMOD)%IC51
          IC52   => WADATS(IMOD)%IC52
          IC61   => WADATS(IMOD)%IC61
          IC62   => WADATS(IMOD)%IC62
          IC71   => WADATS(IMOD)%IC71
          IC72   => WADATS(IMOD)%IC72
          IC81   => WADATS(IMOD)%IC81
          IC82   => WADATS(IMOD)%IC82
          AF11   => WADATS(IMOD)%AF11
        END IF
 
      IF ( NRQSG1 .NE. 0 ) THEN
          IRQSG1 => WADATS(IMOD)%IRQSG1
          IRQSG2 => WADATS(IMOD)%IRQSG2
        END IF
!
      GSTORE => WADATS(IMOD)%GSTORE
      SSTORE => WADATS(IMOD)%SSTORE
!
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** ERROR W3SETA : GRIDS NOT INITIALIZED *** '/      &
               '                    RUN W3NMOD FIRST '/)
 1002 FORMAT (/' *** ERROR W3SETA : ILLEGAL MODEL NUMBER *** '/       &
               '                    IMOD   = ',I10/                   &
               '                    NADATA = ',I10/)
!
!/
!/ End of W3SETA ----------------------------------------------------- /
!/
      END SUBROUTINE W3SETA
!/
!/ End of module W3ADATMD -------------------------------------------- /
!/
      END MODULE W3ADATMD
