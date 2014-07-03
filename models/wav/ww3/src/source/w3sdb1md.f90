!/ ------------------------------------------------------------------- /
      MODULE W3SDB1MD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           J. H. Alves             |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    25-Apr-2007 : Origination of module.              ( version 3.11 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Dummy slot for bottom friction source term.
!
!  2. Variables and types :
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3SDB1    Subr. Public   Battjes and Janssen depth-induced
!                               breaking.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!     See subroutine documentation.
!
!  5. Remarks :
!
!  6. Switches :
!
!     See subroutine documentation.
!
!  7. Source code :
!/
!/ ------------------------------------------------------------------- /
!/
      PUBLIC
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SDB1 (A, K, DEPTH, EMEAN, FMEAN, WNMEAN, S, D )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |                        FORTRAN 90 |
!/                  |           J. H. Alves             |
!/                  |           H. L. Tolman            |
!/                  | Last update :         25-Apr-2007 |
!/                  +-----------------------------------+
!/
!/    25-Apr-2007 : Origination of module.              ( version 3.11 )
!/
!  1. Purpose :
!
!     Compute depth-induced breaking using Battjes and Janssen bore
!     model approach
!
!  2. Method :
!
!       Sdb = - CDB * FMEAN * QB * B * B * SPEC
!
!     Where CDB   =  2.00
!           B     =  HM / HRMS
!           HM    =  GAMMA * DEP
!           GAMMA =  0.73 (average Battjes and Janssen value)
!
!     And QB is estimated by iterations using the nonlinear expression
!
!           1 - QB = HRMS**2
!           ------   -------
!            ln QB    HM**2
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       A       R.A.  I   Action density spectrum (1-D)
!       K       R.A.  I   Wavenumber for entire spectrum.          *)
!       EMEAN   Real  I   Mean wave energy.
!       FMEAN   Real  I   Mean wave frequency.
!       WNMEAN  Real  I   Mean wave number.
!       DEPTH   Real  I   Mean water depth.
!       S       R.A.  O   Source term (1-D version).
!       D       R.A.  O   Diagonal term of derivative (1-D version).
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       STRACE   Subroutine tracing (!/S switch).
!
!  5. Called by :
!
!       W3SRCE   Source term integration.
!       W3EXPO   Point output post-processor.
!       GXEXPO   GrADS point output post-processor.
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!     - Note that the Miche criterion con influence wave growth.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S   Enable subroutine tracing.
!     !/Tn  Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      USE W3CONSTANTS
      USE W3GDATMD, ONLY: NK, NTH, NSPEC, SDBC1, SDBC2, FDONLY
      USE W3ODATMD, ONLY: NDST
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: A(NSPEC), K(NSPEC),                  &
                                 EMEAN, FMEAN, WNMEAN, DEPTH
      REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IS
      REAL                    :: HM, BB, ARG, Q0, QB, B, CBJ, HRMS
!/
!/ ------------------------------------------------------------------- /
!/
!
! 0.  Initialzations ------------------------------------------------- /
!
! 1. Source term parameters ------------------------------------------ /
!
! 1.a. Maximum wave height
! 1.a.1. Simple limit
!
      IF ( FDONLY ) THEN
          HM     = SDBC2 * DEPTH
        ELSE
!
! 1.a.2. Miche style criterion
!
          HM     = SDBC2 / WNMEAN * TANH ( WNMEAN * MAX(DEPTH,0.) )
        END IF
!
! 1.b. Hrms and ratio Hrms / Hmax
!
      HRMS   = SQRT ( 8. * MAX(0.,EMEAN) )
      IF ( HM .GT. 1.E-20 ) THEN
          BB     = HRMS * HRMS / ( HM * HM )
          B      = SQRT(BB)
        ELSE
          B      = 0.
        END IF
!
! 2. Fraction of breaking waves -------------------------------------- /
! 2.a. First guess breaking fraction
!
      IF ( B .LE. 0.5 ) THEN
          Q0     = 0.
        ELSE IF ( B .LE. 1.0 ) THEN
          Q0     = ( 2. * B - 1. ) ** 2
        END IF
!
! 2.b. Iterate to obtain actual breaking fraction
!
      IF ( B .LE. 0.2 ) THEN
          QB     = 0.
        ELSE IF ( B .LT. 1.0 ) THEN
          ARG    = EXP  (( Q0 - 1. ) / BB )
          QB     = Q0 - BB * ( Q0 - ARG ) / ( BB - ARG )
          DO IS=1, 3
            QB     = EXP((QB-1.)/BB)
            END DO
        ELSE
          QB = 1.0 - 1.E-20
        END IF
!
! 3. Estimate the breaking coefficient ------------------------------- /
!
      IF ( ( BB .GT. 0. ) .AND. ( ABS ( BB - QB ) .GT. 0. ) ) THEN
          IF ( BB .LT. 1.0 ) THEN
              CBJ    = SDBC1 * QB * FMEAN * HM * HM / EMEAN
            ELSE
              CBJ    = SDBC1 * FMEAN / EMEAN
            END IF
        ELSE
          CBJ    = 0.
        END IF
!
! 4. Source and diagonal terms --------------------------------------- /
!
      D      = -1. * ( CBJ )
      S      = D * A
!
! ... Test output of arrays
!
      RETURN
!
! Formats
!
!/
!/ End of W3SDB1 ----------------------------------------------------- /
!/
      END SUBROUTINE W3SDB1
!/
!/ End of module W3SDB1MD -------------------------------------------- /
!/
      END MODULE W3SDB1MD
