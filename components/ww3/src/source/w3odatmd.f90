!/ ------------------------------------------------------------------- /
      MODULE W3ODATMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    13-Dec-2004 : Origination.                        ( version 3.06 )
!/    20-Jul-2005 : Adding output fields.               ( version 3.07 )
!/    29-Sep-2005 : Second storage for input bound. sp. ( version 3.08 )
!/                  Add FILED for the dump of data.
!/    26-Jun-2006 : Add output type 6, wave field sep.  ( version 3.09 )
!/                  Wiring of code only.
!/    27-Jun-2006 : Adding file name preamble.          ( version 3.09 )
!/    24-Jul-2006 : Adding unified point output storage.( version 3.10 )
!/    25-Jul-2006 : Originating grid ID for points.     ( version 3.10 )
!/    04-Oct-2006 : Add filter to array pointers.       ( version 3.10 )
!/    30-Oct-2006 : Add pars for partitioning.          ( version 3.10 )
!/    26-Mar-2007 : Add pars for partitioning.          ( version 3.11 )
!/    17-May-2007 : Adding NTPROC/NAPROC separation.    ( version 3.11 )
!/    21-Jun-2007 : Dedicated output processes.         ( version 3.11 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Define data structures to set up wave model grids and aliases
!     to use individual grids transparently. Also includes subroutines
!     to manage data structure and pointing to individual models.
!     This module considers the parameters required for model output.
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      NOUTP     Int.  Public   Number of models in array dim.
!      IOUTP     Int.  Public   Selected model for output, init. at -1.
!      IOSTYP    Int.  Public   Output data server type.
!      NOGRD     I.P.  Public   Number of output fields available.
!      NOSWLL    I.P.  Public   Number of swell fields from part.
!                               to be used in field output.
!      NOEXTR    I.P.  Public   Number of extra (user available)
!                               output fields.
!      IDOUT     C.A.  Public   ID strings for output fields.
!      FNMPRE    Char  Public   File name preamble.
!      UNDEF     Real  Public   Value for undefined parameters in
!                               gridded output fields.
!      UNIPTS    Log.  Public   Flag for unified point output (output
!                               to single file).
!      UPPROC    Log.  Public   FLag for dedicated point output proc.
!      OUTPUT    TYPE  Public   Data structure defining output.
!      OUTPTS    GRID  Public   Array of output for models.
!     ----------------------------------------------------------------
!
!     Elements of OUTPUT are aliased to pointers with the same
!     name. These pointers are defined as :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      NDSO      Int.  Public   General output unit number ("log
!                               file").
!      NDSE      Int.  Public   Error output unit number.
!      NDST      Int.  Public   Test output unit number.
!      SCREEN    Int.  Public   Unit for 'direct' output.
!      NTPROC    Int.  Public   Number of processors.
!      NAPROC    Int.  Public   Number of processors for computation.
!      IAPROC    Int.  Public   Actual processor number (base 1),
!      NAPLOG    Int.  Public   Proc. dealing with log output.
!      NAPOUT    Int.  Public   Proc. dealing with standard output.
!      NAPERR    Int.  Public   Proc. dealing with error output.
!      NAPFLD    Int.  Public   Proc. dealing with raw field output.
!      NAPPNT    Int.  Public   Proc. dealing with raw point output.
!      NAPTRK    Int.  Public   Proc. dealing with track output.
!      NAPRST    Int.  Public   Proc. dealing with restart output.
!      NAPBPT    Int.  Public   Proc. dealing with boundary output.
!      NAPPRT    Int.  Public   Proc. dealing with partition output.
!      TOFRST    I.A.  Public   Times for first output.
!      TONEXT    I.A.  Public   Times for next output.
!      TOLAST    I.A.  Public   Times for last output.
!      TBPI0     I.A   Public   Time of first set of input boundary
!                               spectra.
!      TBPIN     I.A   Public   Id. second set.
!      NDS       I.A.  Public   Data set numbers (see W3INIT).
!      DTOUT     R.A.  Public   Output intervals.
!      FLOUT     L.A.  Public   Output flags.
!      OUT1      TYPE  Public   Data structure of type OTYPE1 with
!                               suppl. data for output type 1.
!      OUT2      TYPE  Public   Data structure of type OTYPE2 with
!                               suppl. data for output type 2.
!      OUT3      TYPE  Public   Data structure of type OTYPE3 with
!                               suppl. data for output type 3.
!      OUT4      TYPE  Public   Data structure of type OTYPE4 with
!                               suppl. data for output type 4.
!      OUT5      TYPE  Public   Data structure of type OTYPE5 with
!                               suppl. data for output type 5.
!      OUT6      TYPE  Public   Data structure of type OTYPE6 with
!                               suppl. data for output type 6.
!     ----------------------------------------------------------------
!
!     Elements of OUT1 are aliased to pointers with the same
!     name. These pointers are defined as :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      IPASS1    Int.  Public   Pass counter for file management,
!                               renamed to IPASS in routine.
!      WRITE1    Int.  Public   Write flag for file management,
!                               renamed to WRITE in routine.
!      NRQGO     Int.  Public   Number of MPI handles W3IOGO.
!      IRQGO     I.A.  Public   Array with MPI handles W3IOGO.
!      FLOGRD    L.A.  Public   FLags for output fields.
!     ----------------------------------------------------------------
!
!     Elements of OUT2 are aliased to pointers with the same
!     name. These pointers are defined as :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      IPASS2    Int.  Public   Pass counter for file management,
!                               renamed to IPASS in routine.
!      NOPTS     Int.  Public   Number of output points.
!      NRQPO(2)  Int.  Public   Number of MPI handles IRQPOn. (!/MPI)
!      IPTINT    I.A.  Public   Id. interpolation counters.
!      IL        I.A.  Public   Number of land points in interpola-
!                               tion box for output point.
!      IW        I.A.  Public   Id. water.
!      II        I.A.  Public   Id. ice.
!      IRQPO1/2  I.A.  Public   Array with MPI handles.       (!/MPI)
!      PTLOC     R.A.  Public   Name of output locations.
!      PTIFAC    R.A.  Public   Id. weights.
!      DPO       R.A.  Public   Interpolated depths.
!      WAO       R.A.  Public   Interpolated wind speeds.
!      WDO       R.A.  Public   Interpolated wind directions.
!      ASO       R.A.  Public   Interpolated air-sea temp. diff.
!      CAO       R.A.  Public   Interpolated current speeds.
!      CDO       R.A.  Public   Interpolated current directions.
!      SPCO      R.A.  Public   Output spectra.
!      PTNME     C.A.  Public   Output locations.
!      GRDID     C.A.  Public   Originating grid ID.
!      O2INIT    Log.  Public   Flag for array initialization.
!      O2IRQI    Log.  Public   Flag for array initialization.
!     ----------------------------------------------------------------
!
!     Elements of OUT3 are aliased to pointers with the same
!     name. These pointers are defined as :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      IPASS3    Int.  Public   Pass counter for file management,
!                               renamed to IPASS in routine.
!      IT0PNT    Int.  Public   Base tag number of MPI communication.
!      IT0TRK    Int.  Public   Base tag number of MPI communication.
!      IT0PRT    Int.  Public   Base tag number of MPI communication.
!      NRQTR     Int.  Public   Number of handles in IRQTR.
!      IRQTR     I.A.  Public   Array with MPI handles.
!      O3INIT    Log.  Public   Flag for array initialization.
!      STOP      Log.  Public   Flag for end of output.
!      MASKn     L.A.  Public   Mask arrays for internal use.
!     ----------------------------------------------------------------
!
!     Elements of OUT4 are aliased to pointers with the same
!     name. These pointers are defined as :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      IFILE4    Int.  Public   File number for output files.
!      NBLKRS    Int.  Public   Number of blocks in communication of
!                               spectra.
!      RSBLKS    Int.  Public   Corresponding block size.
!      NRQSR     Int.  Public   Number of MPI handles.
!      IRQRS     I.A.  Public   Array with MPI handles.
!      IRQRSS    I.A.  Public   Array with MPI handles.
!      VAAUX     I.A.  Public   Aux. spectra storage.
!     ----------------------------------------------------------------
!
!     Elements of OUT5 are aliased to pointers with the same
!     name. These pointers are defined as :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      NBI(2)    Int.  Public   Number of input bound. points.
!      NFBPO     Int.  Public   Number of files for output bound. data.
!      NRQBP(2)  Int.  Public   Number of MPI handles.
!      NBO(2)    I.A.  Public   Number of output bound. pts. per file.
!      NDSL      I.A.  Public   Array with unit numbers.
!      IPBPI     I.A.  Public   Interpolation data input b.p.
!      ISBPI     I.A.  Public   Sea point counters for input b.p.
!      IRQBP1/2  I.A.  Public   Array with MPI handles.
!      X/YBPI    R.A.  Public   Location of input boundary points.
!      RDBPI     R.A.  Public   Interpolation factors input b.p.
!      ABPI0/N   R.A.  Public   Storage of spectra from which to
!                               interpolate b.d.
!      BBPI0/N   R.A.  Public   idem, secondary storage.
!      ABPOS     R.A.  Public   Temporarily storage for output b.d.
!      IPBPO, ISBPO, X/YBPO, RDBPO
!                Misc. Public   Id. for output b.p.
!      FLBPI     Log.  Public   Flag for input of boundary data.
!      FLBPO     Log.  Public   Flag for output of boundary data.
!      FILER/W/D Log.  Public   Read/write flags for file management.
!      SPCONV    Log.  Public   Flag for change of spectral res.
!      O5INIn    Log.  Public   Flag for array initializations.
!     ----------------------------------------------------------------
!
!     Elements of OUT6 are aliased to pointers with the same
!     name. These pointers are defined as :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      IPASS6    Int.  Public   Pass counter for file management,
!                               renamed to IPASS in routine.
!      IHMAX     Int.  Public   Number of discrete spectral levels.
!      IX0/N/S   Int.  Public   First-last-step IX counters.
!      IY0/N/S   Int.  Public   Idem IY counters.
!      HSPMIN    Real  Public   Minimum significant height per part.
!      WSMULT    Real  Public   Multiplier for wind sea boundary.
!      WSCUT     Real  Public   Cut-off wind factor for wind seas.
!      ICPRT     I.A.  Public   Counters for partitions.
!      DTPRT     R.A.  Public   Data from partitions.
!      FLCOMB    Log.  Public   Flag for combining wind seas.
!      FLFORM    Log.  Public   Flag for (un)formatted output
!      O6INIT    Log.  Public   Flag for array initializations.
!     ----------------------------------------------------------------
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3NOUT    Subr. Public   Set number of grids.
!      W3DMO2    Subr. Public   Allocate arrays output type 2.
!      W3DMO3    Subr. Public   Allocate arrays output type 3.
!      W3DMO5    Subr. Public   Allocate arrays output type 5.
!      W3SETO    Subr. Public   Point to selected grid / model.
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
!     !/MPI    MPI specific calls.
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
      INTEGER                 :: NOUTP = -1, IOUTP = -1, IOSTYP = 1
      INTEGER, PARAMETER      :: NOGRD = 31
      INTEGER, PARAMETER      :: NOSWLL=  2
      INTEGER, PARAMETER      :: NOEXTR=  2
      CHARACTER(LEN=20)       :: IDOUT(NOGRD)
      CHARACTER(LEN=80)       :: FNMPRE = './'
      REAL                    :: UNDEF  = -999.9
      LOGICAL                 :: UNIPTS = .FALSE., UPPROC = .FALSE.
!
      DATA IDOUT  / 'Water depth         ' , 'Current vel.        ' , &
                    'Wind speed          ' , 'Air-sea temp. dif.  ' , &
                    'Friction velocity   ' , 'Wave height         ' , &
                    'Mean wave length    ' , 'Mean wave period    ' , &
                    'Mean wave direction ' , 'Mean dir. spread    ' , &
                    'Peak frequency      ' , 'Peak direction      ' , &
                    'Wind sea pk. fr.    ' , 'Wind sea direction  ' , &
                    'Part. wave heigt    ' , 'Part. peak period   ' , &
                    'Part. peak wave l.  ' , 'Part. mean direction' , &
                    'Part. dir. spread   ' , 'Part. wind sea frac.' , &
                    'Total wind sea frac.' , 'Number of partitions' , &
                    'Avg. time step.     ' , 'Cut-off freq.       ' , &
                    'Ice concentration   ' , 'Water level         ' , &
                    'Bottom rms ampl.    ' , 'Bottom rms velocity ' , &
                    'Radiation stresses  ' , 'User defined #1     ' , &
                    'User defined #2     ' /
!/
!/ Data structures
!/
      TYPE OTYPE1
        INTEGER               :: IPASS1
        INTEGER               :: NRQGO
        INTEGER, POINTER      :: IRQGO(:)
        LOGICAL               :: FLOGRD(NOGRD), WRITE1
      END TYPE OTYPE1
!/
      TYPE OTYPE2
        INTEGER               :: IPASS2, NOPTS
        INTEGER               :: NRQPO, NRQPO2
        INTEGER, POINTER      :: IPTINT(:,:), IL(:), IW(:), II(:)
        INTEGER, POINTER      :: IRQPO1(:), IRQPO2(:)
        REAL, POINTER         :: PTLOC(:,:), PTIFAC(:,:),            &
                                 DPO(:), WAO(:), WDO(:), ASO(:),     &
                                 CAO(:), CDO(:), SPCO(:,:)
        CHARACTER(LEN=10), POINTER :: PTNME(:), GRDID(:)
        LOGICAL               :: O2INIT
        LOGICAL               :: O2IRQI
      END TYPE OTYPE2
!/
      TYPE OTYPE3
        INTEGER               :: IPASS3
        INTEGER               :: IT0PNT, IT0TRK, IT0PRT, NRQTR
        INTEGER, POINTER      :: IRQTR(:)
        LOGICAL               :: O3INIT, STOP
        LOGICAL, POINTER      :: MASK1(:,:), MASK2(:,:)
      END TYPE OTYPE3
!/
      TYPE OTYPE4
        INTEGER               :: IFILE4
        INTEGER               :: NRQRS, NBLKRS, RSBLKS
        INTEGER, POINTER      :: IRQRS(:), IRQRSS(:)
        REAL, POINTER         :: VAAUX(:,:,:)
      END TYPE OTYPE4
!/
      TYPE OTYPE5
        INTEGER               :: NBI, NBI2, NFBPO, NBO(0:9),          &
                                 NBO2(0:9), NDSL(9), NKI, NTHI
        INTEGER               :: NRQBP, NRQBP2
        INTEGER, POINTER      :: IPBPI(:,:), ISBPI(:),                &
                                 IPBPO(:,:), ISBPO(:)
        INTEGER, POINTER      :: IRQBP1(:), IRQBP2(:)
        REAL                  :: XFRI, FR1I, TH1I
        REAL, POINTER         :: XBPI(:), YBPI(:), RDBPI(:,:),        &
                                 XBPO(:), YBPO(:), RDBPO(:,:),        &
                                 ABPI0(:,:), ABPIN(:,:), ABPOS(:,:),  &
                                 BBPI0(:,:), BBPIN(:,:)
        LOGICAL               :: O5INI1, O5INI2, O5INI3, O5INI4
        LOGICAL               :: FLBPI, FLBPO, FILER, FILEW, FILED,   &
                                 SPCONV
      END TYPE OTYPE5
!/
      TYPE OTYPE6
        INTEGER               :: IPASS6, IHMAX, IX0, IXN, IXS,        &
                                 IY0, IYN, IYS
        INTEGER, POINTER      :: ICPRT(:,:)
        REAL                  :: HSPMIN, WSMULT, WSCUT
        REAL, POINTER         :: DTPRT(:,:)
        LOGICAL               :: FLFORM, FLCOMB, O6INIT
      END TYPE OTYPE6
!/
      TYPE OUTPUT
        INTEGER               :: NDSO, NDSE, NDST, SCREEN
        INTEGER               :: NTPROC, NAPROC, IAPROC, NAPLOG,      &
                                 NAPOUT, NAPERR, NAPFLD, NAPPNT,      &
                                 NAPTRK, NAPRST, NAPBPT, NAPPRT
        INTEGER               :: TOFRST(2), TONEXT(2,6), TOLAST(2,6), &
                                 TBPI0(2), TBPIN(2), NDS(13)
        REAL                  :: DTOUT(6)
        LOGICAL               :: FLOUT(6)
        TYPE(OTYPE1)          :: OUT1
        TYPE(OTYPE2)          :: OUT2
        TYPE(OTYPE3)          :: OUT3
        TYPE(OTYPE4)          :: OUT4
        TYPE(OTYPE5)          :: OUT5
        TYPE(OTYPE6)          :: OUT6
      END TYPE OUTPUT
!/
!/ Data storage
!/
      TYPE(OUTPUT), TARGET, ALLOCATABLE :: OUTPTS(:)
!/
!/ Data aliasses for structure OUTPUT
!/
      INTEGER, POINTER        :: NDSO, NDSE, NDST, SCREEN
      INTEGER, POINTER        :: NTPROC, NAPROC, IAPROC, NAPLOG,      &
                                 NAPOUT, NAPERR, NAPFLD, NAPPNT,      &
                                 NAPTRK, NAPRST, NAPBPT, NAPPRT
      INTEGER, POINTER        :: TOFRST(:), TONEXT(:,:), TOLAST(:,:), &
                                 TBPI0(:), TBPIN(:), NDS(:)
      REAL, POINTER           :: DTOUT(:)
      LOGICAL, POINTER        :: FLOUT(:)
!/
!/ Data aliasses for substructures for output types
!/ Type 1 ...
!/
      INTEGER, POINTER        :: IPASS1
      INTEGER, POINTER        :: NRQGO
      INTEGER, POINTER        :: IRQGO(:)
      LOGICAL, POINTER        :: FLOGRD(:), WRITE1
!/
!/ Type 2 ...
!/
      INTEGER, POINTER        :: IPASS2, NOPTS
      INTEGER, POINTER        :: NRQPO, NRQPO2
      INTEGER, POINTER        :: IPTINT(:,:), IL(:), IW(:), II(:)
      INTEGER, POINTER        :: IRQPO1(:), IRQPO2(:)
      REAL, POINTER           :: PTLOC(:,:), PTIFAC(:,:),            &
                                 DPO(:), WAO(:), WDO(:), ASO(:),     &
                                 CAO(:), CDO(:), SPCO(:,:)
      CHARACTER(LEN=10), POINTER :: PTNME(:), GRDID(:)
      LOGICAL, POINTER        :: O2INIT
        LOGICAL, POINTER      :: O2IRQI
!/
!/ Type 3 ...
!/
      INTEGER, POINTER        :: IPASS3
      INTEGER, POINTER        :: IT0PNT, IT0TRK, IT0PRT, NRQTR
      INTEGER, POINTER        :: IRQTR(:)
      LOGICAL, POINTER        :: O3INIT, STOP
      LOGICAL, POINTER        :: MASK1(:,:), MASK2(:,:)
!/
!/ Type 4 ...
!/
      INTEGER, POINTER        :: IFILE4
      INTEGER, POINTER        :: NRQRS, NBLKRS, RSBLKS
      INTEGER, POINTER        :: IRQRS(:), IRQRSS(:)
      REAL, POINTER           :: VAAUX(:,:,:)
!/
!/ Type 5 ...
!/
      INTEGER, POINTER        :: NBI, NBI2, NFBPO
      INTEGER, POINTER        :: NBO(:), NBO2(:), NDSL(:), NKI, NTHI
      INTEGER, POINTER        :: NRQBP, NRQBP2
      INTEGER, POINTER        :: IPBPI(:,:), ISBPI(:),                &
                                 IPBPO(:,:), ISBPO(:)
      INTEGER, POINTER        :: IRQBP1(:), IRQBP2(:)
      REAL, POINTER           :: XFRI, FR1I, TH1I
      REAL, POINTER           :: XBPI(:), YBPI(:), RDBPI(:,:),        &
                                 XBPO(:), YBPO(:), RDBPO(:,:),        &
                                 ABPI0(:,:), ABPIN(:,:), ABPOS(:,:),  &
                                 BBPI0(:,:), BBPIN(:,:)
      LOGICAL, POINTER        :: O5INI1, O5INI2, O5INI3, O5INI4
      LOGICAL, POINTER        :: FLBPI, FLBPO, FILER, FILEW, FILED,   &
                                 SPCONV
!/
!/ Type 6 ...
!/
      INTEGER, POINTER        :: IPASS6, IHMAX, IX0, IXN, IXS,        &
                                 IY0, IYN, IYS, ICPRT(:,:)
      REAL, POINTER           :: HSPMIN, WSMULT, WSCUT, DTPRT(:,:)
      LOGICAL, POINTER        :: FLFORM, FLCOMB, O6INIT
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3NOUT ( NDSERR, NDSTST )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         17-May-2007 !
!/                  +-----------------------------------+
!/
!/    13-Dec-2004 : Origination.                        ( version 3.06 )
!/    27-Jun-2006 : Adding file name preamble           ( version 3.09 )
!/    24-Jul-2006 : Adding unified point output storage.( version 3.10 )
!/    04-Oct-2006 : Add filter to array pointers.       ( version 3.10 )
!/    30-Oct-2006 : Add pars for partitioning.          ( version 3.10 )
!/    26-Mar-2007 : Add pars for partitioning.          ( version 3.11 )
!/    17-May-2007 : Adding NTPROC/NAPROC separation.    ( version 3.11 )
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
!       NDSERR  Int.   I   Error output unit number.
!       NDSTST  Int.   I   Test output unit number.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     See module documentation below.
!
!  5. Called by :
!
!     Any main program that uses this grid structure.
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
      USE W3GDATMD, ONLY: NGRIDS, NAUXGR
      USE W3SERVMD, ONLY: EXTCDE
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: NDSERR, NDSTST
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: I, NLOW
!/
!
! -------------------------------------------------------------------- /
! 1.  Test input and module status
!
      IF ( NGRIDS .EQ. -1 ) THEN
          WRITE (NDSERR,1001) NGRIDS
          CALL EXTCDE (1)
        END IF
!
! -------------------------------------------------------------------- /
! 2.  Set variable and allocate arrays
!
      NLOW   = MIN ( 0 , -NAUXGR )
      ALLOCATE ( OUTPTS(NLOW:NGRIDS) )
      NOUTP  = NGRIDS
!
! -------------------------------------------------------------------- /
! 3.  Initialize parameters
!
      DO I=NLOW, NGRIDS
!
        OUTPTS(I)%NDSO   = 6
        OUTPTS(I)%NDSE   = 6
        OUTPTS(I)%NDST   = 6
        OUTPTS(I)%SCREEN = 6
!
        OUTPTS(I)%NTPROC = 1
        OUTPTS(I)%NAPROC = 1
        OUTPTS(I)%IAPROC = 1
        OUTPTS(I)%NAPLOG = 1
        OUTPTS(I)%NAPOUT = 1
        OUTPTS(I)%NAPERR = 1
        OUTPTS(I)%NAPFLD = 1
        OUTPTS(I)%NAPPNT = 1
        OUTPTS(I)%NAPTRK = 1
        OUTPTS(I)%NAPRST = 1
        OUTPTS(I)%NAPBPT = 1
        OUTPTS(I)%NAPPRT = 1
!
        OUTPTS(I)%TBPI0 = (-1,0)
        OUTPTS(I)%TBPIN = (-1,0)
!
        OUTPTS(I)%OUT1%IPASS1 = 0
        OUTPTS(I)%OUT1%NRQGO  = 0
!
        OUTPTS(I)%OUT2%IPASS2 = 0
        OUTPTS(I)%OUT2%NOPTS  = 0
        OUTPTS(I)%OUT2%O2INIT = .FALSE.
        OUTPTS(I)%OUT2%O2IRQI = .FALSE.
!
        OUTPTS(I)%OUT3%IPASS3 = 0
        OUTPTS(I)%OUT3%O3INIT = .FALSE.
        OUTPTS(I)%OUT3%STOP   = .FALSE.
        OUTPTS(I)%OUT3%NRQTR  = 0
!
        OUTPTS(I)%OUT4%IFILE4 = 0
        OUTPTS(I)%OUT4%NRQRS  = 0
!
        OUTPTS(I)%OUT5%O5INI1 = .FALSE.
        OUTPTS(I)%OUT5%O5INI2 = .FALSE.
        OUTPTS(I)%OUT5%O5INI3 = .FALSE.
        OUTPTS(I)%OUT5%O5INI4 = .FALSE.
        OUTPTS(I)%OUT5%FILER  = .TRUE.
        OUTPTS(I)%OUT5%FILEW  = .TRUE.
        OUTPTS(I)%OUT5%FILED  = .TRUE.
!
        OUTPTS(I)%OUT6%IPASS6 = 0
        OUTPTS(I)%OUT6%O6INIT = .FALSE.
!
        END DO
!
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** ERROR W3NOUT : NGRIDS NOT YET SET *** '/         &
               '                    NGRIDS = ',I10/                   &
               '                    RUN W3NMOD FIRST'/)
!
!/
!/ End of W3NOUT ----------------------------------------------------- /
!/
      END SUBROUTINE W3NOUT
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3DMO2 ( IMOD, NDSE, NDST, NPT )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         04-Oct-2006 !
!/                  +-----------------------------------+
!/
!/    10-Nov-2004 : Origination.                        ( version 3.06 )
!/    24-Jul-2006 : Adding unified point output storage.( version 3.10 )
!/    25-Jul-2006 : Originating grid ID for points.     ( version 3.10 )
!/    04-Oct-2006 : Add filter to array pointers.       ( version 3.10 )
!/
!  1. Purpose :
!
!     Initialize an individual data storage for point output.
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
!       NPT     Int.   I   Array size.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     See module documentation below.
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3IOPO    Subr. W3IOPOMD Point output module.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     - Check on input parameters.
!     - Check on previous allocation.
!
!  7. Remarks :
!
!     - W3SETO needs to be called after allocation to point to
!       proper allocated arrays.
!     - Note that NOPTS is overwritten in W3IOPP.
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
      USE W3GDATMD, ONLY: W3SETG, NGRIDS, NAUXGR, IGRID, NSPEC
      USE W3SERVMD, ONLY: EXTCDE
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)           :: IMOD, NDSE, NDST, NPT
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: JGRID, NLOW
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
      NLOW   = MIN ( 0 , -NAUXGR )
      IF ( IMOD.LT.NLOW .OR. IMOD.GT.NOUTP ) THEN
          WRITE (NDSE,1002) IMOD, NLOW, NOUTP
          CALL EXTCDE (2)
        END IF
!
      IF ( OUTPTS(IMOD)%OUT2%O2INIT ) THEN
          WRITE (NDSE,1003)
          CALL EXTCDE (3)
        END IF
!
      JGRID  = IGRID
      IF ( JGRID .NE. IMOD ) CALL W3SETG ( IMOD, NDSE, NDST )
!
! -------------------------------------------------------------------- /
! 2.  Allocate arrays
!
      ALLOCATE ( OUTPTS(IMOD)%OUT2%IPTINT(2,NPT)   ,                  &
                 OUTPTS(IMOD)%OUT2%IL(NPT)         ,                  &
                 OUTPTS(IMOD)%OUT2%IW(NPT)         ,                  &
                 OUTPTS(IMOD)%OUT2%II(NPT)         ,                  &
                 OUTPTS(IMOD)%OUT2%PTIFAC(2,NPT)   ,                  &
                 OUTPTS(IMOD)%OUT2%PTNME(NPT)      ,                  &
                 OUTPTS(IMOD)%OUT2%GRDID(NPT)      ,                  &
                 OUTPTS(IMOD)%OUT2%DPO(NPT)        ,                  &
                 OUTPTS(IMOD)%OUT2%WAO(NPT)        ,                  &
                 OUTPTS(IMOD)%OUT2%WDO(NPT)        ,                  &
                 OUTPTS(IMOD)%OUT2%ASO(NPT)        ,                  &
                 OUTPTS(IMOD)%OUT2%CAO(NPT)        ,                  &
                 OUTPTS(IMOD)%OUT2%CDO(NPT)        ,                  &
                 OUTPTS(IMOD)%OUT2%SPCO(NSPEC,NPT) ,                  &
                 OUTPTS(IMOD)%OUT2%PTLOC(2,NPT)    )
!
      OUTPTS(IMOD)%OUT2%O2INIT = .TRUE.
!
! -------------------------------------------------------------------- /
! 3.  Point to allocated arrays
!
      CALL W3SETO ( IMOD, NDSE, NDST )
!
! -------------------------------------------------------------------- /
! 4.  Update counters in grid
!
      NOPTS  = NPT
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
 1001 FORMAT (/' *** ERROR W3DMO2 : GRIDS NOT INITIALIZED *** '/      &
               '                    RUN W3NMOD FIRST '/)
 1002 FORMAT (/' *** ERROR W3DMO2 : ILLEGAL MODEL NUMBER *** '/       &
               '                    IMOD   = ',I10/                   &
               '                    NLOW   = ',I10/                   &
               '                    NOUTP  = ',I10/)
 1003 FORMAT (/' *** ERROR W3DMO2 : ARRAY(S) ALREADY ALLOCATED *** ')
!
!/
!/ End of W3DMO2 ----------------------------------------------------- /
!/
      END SUBROUTINE W3DMO2
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3DMO3 ( IMOD, NDSE, NDST )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         04-Oct-2006 !
!/                  +-----------------------------------+
!/
!/    24-Nov-2004 : Origination.                        ( version 3.06 )
!/    04-Oct-2006 : Add filter to array pointers.       ( version 3.10 )
!/
!  1. Purpose :
!
!     Initialize an individual data storage for track output.
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
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     See module documentation below.
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3IOTR    Subr. W3IOTRMD Track output module.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     - Check on input parameters.
!     - Check on previous allocation.
!
!  7. Remarks :
!
!     - W3SETO needs to be called after allocation to point to
!       proper allocated arrays.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/SHRD, !/DIST, !/MPI
!              Shared / distributed memory model
!
!     !/S      Enable subroutine tracing.
!     !/T    Enable test output
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3GDATMD, ONLY: W3SETG, NGRIDS, IGRID, NX, NY
      USE W3SERVMD, ONLY: EXTCDE
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)           :: IMOD, NDSE, NDST
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: JGRID
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
      IF ( IMOD.LT.1 .OR. IMOD.GT.NOUTP ) THEN
          WRITE (NDSE,1002) IMOD, NOUTP
          CALL EXTCDE (2)
        END IF
!
      IF ( OUTPTS(IMOD)%OUT3%O3INIT ) THEN
          WRITE (NDSE,1003)
          CALL EXTCDE (3)
        END IF
!
      JGRID  = IGRID
      IF ( JGRID .NE. IMOD ) CALL W3SETG ( IMOD, NDSE, NDST )
!
! -------------------------------------------------------------------- /
! 2.  Allocate arrays
!
      ALLOCATE ( OUTPTS(IMOD)%OUT3%MASK1(NY,NX) ,                     &
                 OUTPTS(IMOD)%OUT3%MASK2(NY,NX) )
!
      OUTPTS(IMOD)%OUT3%O3INIT = .TRUE.
!
! -------------------------------------------------------------------- /
! 3.  Point to allocated arrays
!
      CALL W3SETO ( IMOD, NDSE, NDST )
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
 1001 FORMAT (/' *** ERROR W3DMO3 : GRIDS NOT INITIALIZED *** '/      &
               '                    RUN W3NMOD FIRST '/)
 1002 FORMAT (/' *** ERROR W3DMO3 : ILLEGAL MODEL NUMBER *** '/       &
               '                    IMOD   = ',I10/                   &
               '                    NOUTP  = ',I10/)
 1003 FORMAT (/' *** ERROR W3DMO3 : ARRAY(S) ALREADY ALLOCATED *** ')
!
!/
!/ End of W3DMO3 ----------------------------------------------------- /
!/
      END SUBROUTINE W3DMO3
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3DMO5 ( IMOD, NDSE, NDST, IBLOCK )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         04-Oct-2006 !
!/                  +-----------------------------------+
!/
!/    13-Dec-2004 : Origination.                        ( version 3.06 )
!/    06-Sep-2005 : Second storage for input bound. sp. ( version 3.08 )
!/    04-Oct-2006 : Add filter to array pointers.       ( version 3.10 )
!/
!  1. Purpose :
!
!     Initialize an individual data storage for boundary data.
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
!       IBLOCK  Int.   I   Select block to allocate.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     See module documentation below.
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3IOBC    Subr. W3IOBCMD Boundary data output module.
!      W3IOGR    Subr. W3IOGRMD Grid data output module.
!      W3WAVE    Subr. W3WAVEMD Actual wave model routine.
!      WW3_GRID  Prog.   N/A    Grid preprocessing program.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     - Check on input parameters.
!     - Check on previous allocation.
!
!  7. Remarks :
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
      USE W3GDATMD, ONLY: W3SETG, NGRIDS, IGRID, NX, NY, NSPEC
      USE W3SERVMD, ONLY: EXTCDE
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)           :: IMOD, NDSE, NDST, IBLOCK
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: JGRID
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
      IF ( IMOD.LT.1 .OR. IMOD.GT.NOUTP ) THEN
          WRITE (NDSE,1002) IMOD, NOUTP
          CALL EXTCDE (2)
        END IF
!
! -------------------------------------------------------------------- /
! 2.  Allocate arrays and reset pointers
!
      SELECT CASE (IBLOCK)
!
      CASE (1)
!
        ALLOCATE ( OUTPTS(IMOD)%OUT5%IPBPI(NBI,4),                    &
                   OUTPTS(IMOD)%OUT5%ISBPI(NBI)  ,                    &
                   OUTPTS(IMOD)%OUT5%XBPI(NBI)   ,                    &
                   OUTPTS(IMOD)%OUT5%YBPI(NBI)   ,                    &
                   OUTPTS(IMOD)%OUT5%RDBPI(NBI,4) )
!
        IPBPI  => OUTPTS(IMOD)%OUT5%IPBPI
        ISBPI  => OUTPTS(IMOD)%OUT5%ISBPI
        XBPI   => OUTPTS(IMOD)%OUT5%XBPI
        YBPI   => OUTPTS(IMOD)%OUT5%YBPI
        RDBPI  => OUTPTS(IMOD)%OUT5%RDBPI
!
        OUTPTS(IMOD)%OUT5%O5INI1 = .TRUE.
!
      CASE (2)
!
        ALLOCATE ( OUTPTS(IMOD)%OUT5%IPBPO(NBO(NFBPO),4),             &
                   OUTPTS(IMOD)%OUT5%ISBPO(4*NBO(NFBPO)),             &
                   OUTPTS(IMOD)%OUT5%XBPO(NBO(NFBPO))   ,             &
                   OUTPTS(IMOD)%OUT5%YBPO(NBO(NFBPO))   ,             &
                   OUTPTS(IMOD)%OUT5%RDBPO(NBO(NFBPO),4) )
!
        IPBPO  => OUTPTS(IMOD)%OUT5%IPBPO
        ISBPO  => OUTPTS(IMOD)%OUT5%ISBPO
        XBPO   => OUTPTS(IMOD)%OUT5%XBPO
        YBPO   => OUTPTS(IMOD)%OUT5%YBPO
        RDBPO  => OUTPTS(IMOD)%OUT5%RDBPO
!
        OUTPTS(IMOD)%OUT5%O5INI2 = .TRUE.
!
      CASE (3)
!
        ALLOCATE ( OUTPTS(IMOD)%OUT5%ABPI0(NSPEC,0:NBI2),             &
                   OUTPTS(IMOD)%OUT5%ABPIN(NSPEC,0:NBI2),             &
                   OUTPTS(IMOD)%OUT5%BBPI0(NSPEC,0:NBI),              &
                   OUTPTS(IMOD)%OUT5%BBPIN(NSPEC,0:NBI) )
!
        ABPI0  => OUTPTS(IMOD)%OUT5%ABPI0
        ABPIN  => OUTPTS(IMOD)%OUT5%ABPIN
        BBPI0  => OUTPTS(IMOD)%OUT5%BBPI0
        BBPIN  => OUTPTS(IMOD)%OUT5%BBPIN
        BBPI0 = -1.
!
        OUTPTS(IMOD)%OUT5%O5INI3 = .TRUE.
!
      CASE (4)
!
        ALLOCATE ( OUTPTS(IMOD)%OUT5%ABPOS(NSPEC,0:NBO2(NFBPO)) )
!
        ABPOS  => OUTPTS(IMOD)%OUT5%ABPOS
!
        OUTPTS(IMOD)%OUT5%O5INI4 = .TRUE.
!
      CASE DEFAULT
        WRITE (NDSE,1010)
        CALL EXTCDE (10)
!
      END SELECT
!
! -------------------------------------------------------------------- /
!
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** ERROR W3DMO5 : GRIDS NOT INITIALIZED *** '/      &
               '                    RUN W3NMOD FIRST '/)
 1002 FORMAT (/' *** ERROR W3DMO5 : ILLEGAL MODEL NUMBER *** '/       &
               '                    IMOD   = ',I10/                   &
               '                    NOUTP  = ',I10/)
 1010 FORMAT (/' *** ERROR W3DMO5 : ILLEGAL BLOCK NUMBER  *** '/      &
               '                    IBLOCK = ',I10/)
!
!/
!/ End of W3DMO5 ----------------------------------------------------- /
!/
      END SUBROUTINE W3DMO5
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SETO ( IMOD, NDSERR, NDSTST )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         17-May-2007 !
!/                  +-----------------------------------+
!/
!/    13-Dec-2004 : Origination.                        ( version 3.06 )
!/    06-Sep-2005 : Second storage for input bound. sp. ( version 3.08 )
!/    24-Jul-2006 : Adding unified point output storage.( version 3.10 )
!/    25-Jul-2006 : Originating grid ID for points.     ( version 3.10 )
!/    04-Oct-2006 : Add filter to array pointers.       ( version 3.10 )
!/    30-Oct-2006 : Add pars for partitioning.          ( version 3.10 )
!/    26-Mar-2007 : Add pars for partitioning.          ( version 3.11 )
!/    17-May-2007 : Adding NTPROC/NAPROC separation.    ( version 3.11 )
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
!       NDSERR  Int.   I   Error output unit number.
!       NDSTST  Int.   I   Test output unit number.
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
!     !/MPI  MPI specific calls.
!     !/S    Enable subroutine tracing.
!     !/T    Enable test output
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3GDATMD, ONLY: NAUXGR
      USE W3SERVMD, ONLY: EXTCDE
!
      IMPLICIT NONE
!
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: IMOD, NDSERR, NDSTST
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: NLOW
!/
!
! -------------------------------------------------------------------- /
! 1.  Test input and module status
!
      IF ( NOUTP .EQ. -1 ) THEN
          WRITE (NDSERR,1001)
          CALL EXTCDE (1)
        END IF
!
      NLOW   = MIN ( 0 , -NAUXGR )
      IF ( IMOD.LT.NLOW .OR. IMOD.GT.NOUTP ) THEN
          WRITE (NDSERR,1002) IMOD, NLOW, NOUTP
          CALL EXTCDE (2)
        END IF
!
! -------------------------------------------------------------------- /
! 2.  Set model number
!
      IOUTP  = IMOD
!
! -------------------------------------------------------------------- /
! 3.  Set pointers in structure OUTPUT
!
      NDSO   => OUTPTS(IMOD)%NDSO
      NDSE   => OUTPTS(IMOD)%NDSE
      NDST   => OUTPTS(IMOD)%NDST
      SCREEN => OUTPTS(IMOD)%SCREEN
!
      NTPROC => OUTPTS(IMOD)%NTPROC
      NAPROC => OUTPTS(IMOD)%NAPROC
      IAPROC => OUTPTS(IMOD)%IAPROC
      NAPLOG => OUTPTS(IMOD)%NAPLOG
      NAPOUT => OUTPTS(IMOD)%NAPOUT
      NAPERR => OUTPTS(IMOD)%NAPERR
      NAPFLD => OUTPTS(IMOD)%NAPFLD
      NAPPNT => OUTPTS(IMOD)%NAPPNT
      NAPTRK => OUTPTS(IMOD)%NAPTRK
      NAPRST => OUTPTS(IMOD)%NAPRST
      NAPBPT => OUTPTS(IMOD)%NAPBPT
      NAPPRT => OUTPTS(IMOD)%NAPPRT
!
      TOFRST => OUTPTS(IMOD)%TOFRST
      TONEXT => OUTPTS(IMOD)%TONEXT
      TOLAST => OUTPTS(IMOD)%TOLAST
      TBPI0  => OUTPTS(IMOD)%TBPI0
      TBPIN  => OUTPTS(IMOD)%TBPIN
      NDS    => OUTPTS(IMOD)%NDS
!
      DTOUT  => OUTPTS(IMOD)%DTOUT
      FLOUT  => OUTPTS(IMOD)%FLOUT
!
      IPASS1 => OUTPTS(IMOD)%OUT1%IPASS1
      WRITE1 => OUTPTS(IMOD)%OUT1%WRITE1
      NRQGO  => OUTPTS(IMOD)%OUT1%NRQGO
      IF ( NRQGO .NE. 0 ) IRQGO  => OUTPTS(IMOD)%OUT1%IRQGO
      FLOGRD => OUTPTS(IMOD)%OUT1%FLOGRD
!
      IPASS2 => OUTPTS(IMOD)%OUT2%IPASS2
      NOPTS  => OUTPTS(IMOD)%OUT2%NOPTS
      NRQPO  => OUTPTS(IMOD)%OUT2%NRQPO
      NRQPO2 => OUTPTS(IMOD)%OUT2%NRQPO2
      O2INIT => OUTPTS(IMOD)%OUT2%O2INIT
      O2IRQI => OUTPTS(IMOD)%OUT2%O2IRQI
!
      IF ( O2INIT ) THEN
          IPTINT => OUTPTS(IMOD)%OUT2%IPTINT
          IL     => OUTPTS(IMOD)%OUT2%IL
          IW     => OUTPTS(IMOD)%OUT2%IW
          II     => OUTPTS(IMOD)%OUT2%II
          PTLOC  => OUTPTS(IMOD)%OUT2%PTLOC
          PTIFAC => OUTPTS(IMOD)%OUT2%PTIFAC
          DPO    => OUTPTS(IMOD)%OUT2%DPO
          WAO    => OUTPTS(IMOD)%OUT2%WAO
          WDO    => OUTPTS(IMOD)%OUT2%WDO
          ASO    => OUTPTS(IMOD)%OUT2%ASO
          CAO    => OUTPTS(IMOD)%OUT2%CAO
          CDO    => OUTPTS(IMOD)%OUT2%CDO
          SPCO   => OUTPTS(IMOD)%OUT2%SPCO
          PTNME  => OUTPTS(IMOD)%OUT2%PTNME
          GRDID  => OUTPTS(IMOD)%OUT2%GRDID
        END IF
!
      IF ( O2IRQI ) THEN
          IRQPO1 => OUTPTS(IMOD)%OUT2%IRQPO1
          IRQPO2 => OUTPTS(IMOD)%OUT2%IRQPO2
        END IF
!
      IPASS3 => OUTPTS(IMOD)%OUT3%IPASS3
      IT0PNT => OUTPTS(IMOD)%OUT3%IT0PNT
      IT0TRK => OUTPTS(IMOD)%OUT3%IT0TRK
      IT0PRT => OUTPTS(IMOD)%OUT3%IT0PRT
      NRQTR  => OUTPTS(IMOD)%OUT3%NRQTR
      IF ( NRQTR .EQ. 0 ) IRQTR  => OUTPTS(IMOD)%OUT3%IRQTR
      O3INIT => OUTPTS(IMOD)%OUT3%O3INIT
      STOP   => OUTPTS(IMOD)%OUT3%STOP
!
      IF ( O3INIT ) THEN
          MASK1  => OUTPTS(IMOD)%OUT3%MASK1
          MASK2  => OUTPTS(IMOD)%OUT3%MASK2
        END IF
!
      IFILE4 => OUTPTS(IMOD)%OUT4%IFILE4
      NRQRS  => OUTPTS(IMOD)%OUT4%NRQRS
      NBLKRS => OUTPTS(IMOD)%OUT4%NBLKRS
      RSBLKS => OUTPTS(IMOD)%OUT4%RSBLKS
      IF ( NRQRS .NE. 0 ) THEN
          IRQRS  => OUTPTS(IMOD)%OUT4%IRQRS
          IRQRSS => OUTPTS(IMOD)%OUT4%IRQRSS
          VAAUX  => OUTPTS(IMOD)%OUT4%VAAUX
        END IF
!
      NBI    => OUTPTS(IMOD)%OUT5%NBI
      NBI2   => OUTPTS(IMOD)%OUT5%NBI2
      NFBPO  => OUTPTS(IMOD)%OUT5%NFBPO
      NRQBP  => OUTPTS(IMOD)%OUT5%NRQBP
      NRQBP2 => OUTPTS(IMOD)%OUT5%NRQBP2
      NBO    => OUTPTS(IMOD)%OUT5%NBO
      NBO2   => OUTPTS(IMOD)%OUT5%NBO2
      NDSL   => OUTPTS(IMOD)%OUT5%NDSL
      NKI    => OUTPTS(IMOD)%OUT5%NKI
      NTHI   => OUTPTS(IMOD)%OUT5%NTHI
      XFRI   => OUTPTS(IMOD)%OUT5%XFRI
      FR1I   => OUTPTS(IMOD)%OUT5%FR1I
      TH1I   => OUTPTS(IMOD)%OUT5%TH1I
      FLBPI  => OUTPTS(IMOD)%OUT5%FLBPI
      FLBPO  => OUTPTS(IMOD)%OUT5%FLBPO
      FILER  => OUTPTS(IMOD)%OUT5%FILER
      FILEW  => OUTPTS(IMOD)%OUT5%FILEW
      FILED  => OUTPTS(IMOD)%OUT5%FILED
      SPCONV => OUTPTS(IMOD)%OUT5%SPCONV
      O5INI1 => OUTPTS(IMOD)%OUT5%O5INI1
      O5INI2 => OUTPTS(IMOD)%OUT5%O5INI2
      O5INI3 => OUTPTS(IMOD)%OUT5%O5INI3
      O5INI4 => OUTPTS(IMOD)%OUT5%O5INI4
!
      IF ( O5INI1 ) THEN
          IPBPI  => OUTPTS(IMOD)%OUT5%IPBPI
          ISBPI  => OUTPTS(IMOD)%OUT5%ISBPI
          XBPI   => OUTPTS(IMOD)%OUT5%XBPI
          YBPI   => OUTPTS(IMOD)%OUT5%YBPI
          RDBPI  => OUTPTS(IMOD)%OUT5%RDBPI
        END IF
!
      IF ( O5INI2 ) THEN
          IPBPO  => OUTPTS(IMOD)%OUT5%IPBPO
          ISBPO  => OUTPTS(IMOD)%OUT5%ISBPO
          XBPO   => OUTPTS(IMOD)%OUT5%XBPO
          YBPO   => OUTPTS(IMOD)%OUT5%YBPO
          RDBPO  => OUTPTS(IMOD)%OUT5%RDBPO
        END IF
!
      IF ( O5INI3 ) THEN
          ABPI0  => OUTPTS(IMOD)%OUT5%ABPI0
          ABPIN  => OUTPTS(IMOD)%OUT5%ABPIN
          BBPI0  => OUTPTS(IMOD)%OUT5%BBPI0
          BBPIN  => OUTPTS(IMOD)%OUT5%BBPIN
        END IF
!
      IF ( O5INI4 ) THEN
          ABPOS  => OUTPTS(IMOD)%OUT5%ABPOS
        END IF
!
      IF ( NRQBP  .NE. 0 ) IRQBP1 => OUTPTS(IMOD)%OUT5%IRQBP1
      IF ( NRQBP2 .NE. 0 ) IRQBP2 => OUTPTS(IMOD)%OUT5%IRQBP2
!
      IPASS6 => OUTPTS(IMOD)%OUT6%IPASS6
      IHMAX  => OUTPTS(IMOD)%OUT6%IHMAX
      HSPMIN => OUTPTS(IMOD)%OUT6%HSPMIN
      WSMULT => OUTPTS(IMOD)%OUT6%WSMULT
      WSCUT  => OUTPTS(IMOD)%OUT6%WSCUT
      IX0    => OUTPTS(IMOD)%OUT6%IX0
      IXN    => OUTPTS(IMOD)%OUT6%IXN
      IXS    => OUTPTS(IMOD)%OUT6%IXS
      IY0    => OUTPTS(IMOD)%OUT6%IY0
      IYN    => OUTPTS(IMOD)%OUT6%IYN
      IYS    => OUTPTS(IMOD)%OUT6%IYS
      ICPRT  => OUTPTS(IMOD)%OUT6%ICPRT
      DTPRT  => OUTPTS(IMOD)%OUT6%DTPRT
      FLCOMB => OUTPTS(IMOD)%OUT6%FLCOMB
      FLFORM => OUTPTS(IMOD)%OUT6%FLFORM
      O6INIT => OUTPTS(IMOD)%OUT6%O6INIT
!
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** ERROR W3SETO : GRIDS NOT INITIALIZED *** '/      &
               '                    RUN W3NMOD FIRST '/)
 1002 FORMAT (/' *** ERROR W3SETO : ILLEGAL MODEL NUMBER *** '/       &
               '                    IMOD   = ',I10/                   &
               '                    NLOW   = ',I10/                   &
               '                    NOUTP  = ',I10/)
!
!/
!/ End of W3SETO ----------------------------------------------------- /
!/
      END SUBROUTINE W3SETO
!/
!/ End of module W3ODATMD -------------------------------------------- /
!/
      END MODULE W3ODATMD
