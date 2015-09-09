!/ ------------------------------------------------------------------- /
      MODULE W3WDATMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    22-Oct-2004 : Origination.                        ( version 3.06 )
!/    13-Jun-2006 : Allocate VA consistent with MPI     ( version 3.09 )
!/                  data types and initialize as needed.
!/    05-Jul-2006 : Consolidate stress vector.          ( version 3.09 )
!/    04-Oct-2006 : Add filter to array pointers.       ( version 3.10 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Define data structures to set up wave model dynamic data for
!     several models simultaneously.
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      NWDATA    Int.  Public   Number of models in array dim.
!      IWDATA    Int.  Public   Selected model for output, init. at -1.
!      WDATA     TYPE  Public   Basic data structure.
!      WDATAS    WDATA Public   Array of data structures.
!     ----------------------------------------------------------------
!
!     All elements of WDATA are aliased to pointers with the same
!     name. These pointers are defined as :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      TIME      I.A.  Public   Valid time for spectra.
!      TLEV      I.A.  Public   Valid time for water levels.
!      TICE      I.A.  Public   Valid time for ice.
!      VA        R.A.  Public   Storage array for spectra.
!      WLV       R.A.  Public   Water levels.
!      ICE       R.A.  Public   Ice coverage.
!      UST       R.A.  Public   Friction velocity (absolute).
!      USTDIR    R.A.  Public   Friction velocity direction.
!      ASF       R.A.  Public   Stability correction factor.
!      FPIS      R.A.  Public   Input peak frequencies.
!      DINIT     Log.  Public   Flag for array initialization.
!      FL_ALL    Log.  Public   Flag for initializing all arrays,
!                               otherwise VA is skipped.
!     ----------------------------------------------------------------
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3NDAT    Subr. Public   Set number of grids/models.
!      W3DIMW    Subr. Public   Set dimensions of arrays.
!      W3SETW    Subr. Public   Point to selected grid / model.
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
!     !/S    Enable subroutine tracing.
!     !/T    Enable test output
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
      PUBLIC
!/
!/ Conventional declarations
!/
      INTEGER                 :: NWDATA = -1, IWDATA = -1
!/
!/ Data structures
!/
      TYPE WDATA
        INTEGER               :: TIME(2), TLEV(2), TICE(2)
        REAL, POINTER         :: VA(:,:), WLV(:), ICE(:), UST(:),     &
                                 USTDIR(:), ASF(:), FPIS(:)
        LOGICAL               :: DINIT, FL_ALL
      END TYPE WDATA
!
!/
!/ Data storage
!/
      TYPE(WDATA), TARGET, ALLOCATABLE :: WDATAS(:)
!/
!/ Data aliasses for structure WDATA(S)
!/
      INTEGER, POINTER        :: TIME(:), TLEV(:), TICE(:)
      REAL, POINTER           :: VA(:,:), WLV(:), ICE(:), UST(:),     &
                                 USTDIR(:), ASF(:), FPIS(:)
      LOGICAL, POINTER        :: DINIT, FL_ALL
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3NDAT ( NDSE, NDST )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         31-Mar-2004 !
!/                  +-----------------------------------+
!/
!/    31-Mar-2004 : Origination.                        ( version 3.06 )
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
      ALLOCATE ( WDATAS(0:NGRIDS) )
      NWDATA = NGRIDS
!
! -------------------------------------------------------------------- /
! 3.  Initialize parameters
!
      DO I=0, NGRIDS
        WDATAS(I)%DINIT  = .FALSE.
        WDATAS(I)%FL_ALL = .FALSE.
        END DO
!
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** ERROR W3NDAT : NGRIDS NOT YET SET *** '/         &
               '                    NGRIDS = ',I10/                   &
               '                    RUN W3NMOD FIRST'/)
!
!/
!/ End of W3NDAT ----------------------------------------------------- /
!/
      END SUBROUTINE W3NDAT
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3DIMW  ( IMOD, NDSE, NDST, F_ONLY )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         04-Oct-2006 !
!/                  +-----------------------------------+
!/
!/    22-Oct-2004 : Origination.                        ( version 3.06 )
!/    13-Jun-2006 : Allocate VA consistent with MPI     ( version 3.09 )
!/                  data types and initialize as needed.
!/    05-Jul-2006 : Consolidate stress vector.          ( version 3.09 )
!/    04-Oct-2006 : Add filter to array pointers.       ( version 3.10 )
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
!       F_ONLY  L.O.   I   FLag for initializing field arrays only.
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
!      W3IORS    Subr. W3IORSMD Restart file IO routine.
!      WW3_SHEL  Prog.   N/A    Main wave model driver.
!      WW3_STRT  Prog.   N/A    Initial conditions program.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     - Check on input parameters.
!     - Check on previous allocation.
!
!  7. Remarks :
!
!     - W3SETW needs to be called after allocation to point to
!       proper allocated arrays.
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
      USE W3GDATMD, ONLY: NGRIDS, IGRID, W3SETG, NSPEC, NSEA, NSEAL
      USE W3ODATMD, ONLY: NAPROC
      USE W3SERVMD, ONLY: EXTCDE
!
      IMPLICIT NONE
!
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)           :: IMOD, NDSE, NDST
      LOGICAL, INTENT(IN), OPTIONAL :: F_ONLY
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: JGRID, NSEALM, NSEATM
!/
!
! -------------------------------------------------------------------- /
! 1.  Test input and module status
!
      IF ( PRESENT(F_ONLY) ) THEN
          FL_ALL = .NOT. F_ONLY
        ELSE
          FL_ALL = .TRUE.
        END IF
!
      IF ( NGRIDS .EQ. -1 ) THEN
          WRITE (NDSE,1001)
          CALL EXTCDE (1)
        END IF
!
      IF ( IMOD.LT.1 .OR. IMOD.GT.NWDATA ) THEN
          WRITE (NDSE,1002) IMOD, NWDATA
          CALL EXTCDE (2)
        END IF
!
      IF ( WDATAS(IMOD)%DINIT ) THEN
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
      NSEALM = 1 + (NSEA-1)/NAPROC
      NSEATM = NSEALM * NAPROC
!
      IF ( FL_ALL ) THEN
          ALLOCATE ( WDATAS(IMOD)%VA(NSPEC,0:NSEALM) )
          IF ( NSEAL .NE. NSEALM ) WDATAS(IMOD)%VA(:,NSEALM) = 0.
        END IF
!
      ALLOCATE ( WDATAS(IMOD)%WLV(NSEA),                              &
                 WDATAS(IMOD)%ICE(NSEA),                              &
                 WDATAS(IMOD)%UST(0:NSEATM),                          &
                 WDATAS(IMOD)%USTDIR(0:NSEATM),                       &
                 WDATAS(IMOD)%ASF(NSEATM),                            &
                 WDATAS(IMOD)%FPIS(NSEATM) )
!
      WDATAS(IMOD)%DINIT  = .TRUE.
      CALL W3SETW ( IMOD, NDSE, NDST )
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
 1001 FORMAT (/' *** ERROR W3DIMW : GRIDS NOT INITIALIZED *** '/      &
               '                    RUN W3NMOD FIRST '/)
 1002 FORMAT (/' *** ERROR W3DIMW : ILLEGAL MODEL NUMBER *** '/       &
               '                    IMOD   = ',I10/                   &
               '                    NWDATA = ',I10/)
 1003 FORMAT (/' *** ERROR W3DIMW : ARRAY(S) ALREADY ALLOCATED *** ')
!
! -------------------------------------------------------------------- /
! 3.  Point to allocated arrays
!
      CALL W3SETW ( IMOD, NDSE, NDST )
!
! -------------------------------------------------------------------- /
! 4.  Update counters in grid
!/
!/ End of W3DIMW ----------------------------------------------------- /
!/
      END SUBROUTINE W3DIMW
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SETW ( IMOD, NDSE, NDST )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         04-Oct-2006 !
!/                  +-----------------------------------+
!/
!/    31-Mar-2004 : Origination.                        ( version 3.06 )
!/    05-Jul-2006 : Consolidate stress vector.          ( version 3.09 )
!/    04-Oct-2006 : Add filter to array pointers.       ( version 3.10 )
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
!     See module documentation.
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
!     !/S    Enable subroutine tracing.
!     !/T    Enable test output
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
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
      IF ( NWDATA .EQ. -1 ) THEN
          WRITE (NDSE,1001)
          CALL EXTCDE (1)
        END IF
!
      IF ( IMOD.LT.0 .OR. IMOD.GT.NWDATA ) THEN
          WRITE (NDSE,1002) IMOD, NWDATA
          CALL EXTCDE (2)
        END IF
!
! -------------------------------------------------------------------- /
! 2.  Set model numbers
!
      IWDATA = IMOD
!
! -------------------------------------------------------------------- /
! 3.  Set pointers
!
      TIME   => WDATAS(IMOD)%TIME
      TLEV   => WDATAS(IMOD)%TLEV
      TICE   => WDATAS(IMOD)%TICE
      DINIT  => WDATAS(IMOD)%DINIT
      FL_ALL => WDATAS(IMOD)%FL_ALL
!
      IF ( DINIT ) THEN
          IF ( FL_ALL ) VA     => WDATAS(IMOD)%VA
          WLV    => WDATAS(IMOD)%WLV
          ICE    => WDATAS(IMOD)%ICE
          UST    => WDATAS(IMOD)%UST
          USTDIR => WDATAS(IMOD)%USTDIR
          ASF    => WDATAS(IMOD)%ASF
          FPIS   => WDATAS(IMOD)%FPIS
        END IF
!
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** ERROR W3SETW : GRIDS NOT INITIALIZED *** '/      &
               '                    RUN W3NMOD FIRST '/)
 1002 FORMAT (/' *** ERROR W3SETW : ILLEGAL MODEL NUMBER *** '/       &
               '                    IMOD   = ',I10/                   &
               '                    NWDATA = ',I10/)
!
!/
!/ End of W3SETW ----------------------------------------------------- /
!/
      END SUBROUTINE W3SETW
!/
!/ End of module W3WDATMD -------------------------------------------- /
!/
      END MODULE W3WDATMD
