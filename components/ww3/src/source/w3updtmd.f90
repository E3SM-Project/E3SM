!/ ------------------------------------------------------------------- /
      MODULE W3UPDTMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    21-Jan-2000 : Origination.                        ( version 2.00 )
!/    24-Jan-2001 : Flat grid version.                  ( version 2.06 )
!/    02-Apr-2001 : Adding sub-grid obstacles.          ( version 2.10 )
!/    18-May-2001 : Clean up and bug fixes.             ( version 2.11 )
!/    11-Jan-2002 : Sub-grid ice.                       ( version 2.15 )
!/    30-Apr-2002 : Water level fixes.                  ( version 2.20 )
!/    13-Nov-2002 : Add stress vector.                  ( version 3.00 )
!/    26-Dec-2002 : Moving grid wind correction.        ( version 3.02 )
!/    15-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    15-Jul-2005 : Adding MAPST2.                      ( version 3.07 )
!/    07-Sep-2005 : Upgrading W3UBPT.                   ( version 3.08 )
!/    04-Jul-2006 : Consolidate stress arrays.          ( version 3.09 )
!/    11-Jan-2007 : Clean-up W3UTRN boundary points.    ( version 3.10 )
!/    11-May-2007 : Adding NTPROC/NAPROC separation.    ( version 3.11 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Bundles all input updating routines for WAVEWATCH III.
!
!  2. Variables and types :
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3UCUR    Subr. Public   Update current fields.
!      W3UWND    Subr. Public   Update wind fields.
!      W3UINI    Subr. Public   Update initial conditions.
!      W3UBPT    Subr. Public   Update boundary conditions.
!      W3UICE    Subr. Public   Update ice concentrations.
!      W3ULEV    Subr. Public   Update water levels.
!      W3UTRN    Subr. Public   Update cell boundary transparancies.
!      W3DDXY    Subr. Public   Calculate depth dirivatives.
!      W3DCXY    Subr. Public   Calculate current dirivatives.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      DSEC21    Func. W3TIMEMD Difference in time.
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!      EXTCDE    Subr. W3SERVMD Exit program with error code.
!      PRTBLK    Subr. W3ARRYMD Print plot output.
!      PRT2DS    Subr. W3ARRYMD Print plot output.
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!  6. Switches :
!
!       !/SHRD   Switch for shared / distributed memory architecture.
!       !/DIST   Id.
!
!       !/OMP1   OpenMP compiler directives.
!
!       !/LLG    Spherical grid.
!       !/XYG    Cartesian grid.
!
!       !/CRT1   Linear current interpolation.
!       !/CRT2   Quasi-quadratic current interpolation.
!
!       !/WNT0   No wind interpolation.
!       !/WNT1   Linear wind interpolation.
!       !/WNT2   Energy conservation in wind interpolation.
!
!       !/RWND   Use wind speeds relative to currents.
!
!       !/STAB2  Calculate effective wind speed factor for stability
!                to be used with !/ST2.
!
!       !/S      Enable subroutine tracing.
!       !/Tn     Test output
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3CONSTANTS
      USE W3ODATMD, ONLY: NDSE, NDST, NAPROC, IAPROC, NAPERR
      USE W3TIMEMD, ONLY: DSEC21
!/
!/ ------------------------------------------------------------------- /
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3UCUR ( FLFRST )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         15-Dec-2004 |
!/                  +-----------------------------------+
!/
!/    09-Dec-1996 : Final FORTRAN 77                    ( version 1.18 )
!/    20-Dec-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    15-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/
!  1. Purpose :
!
!     Interpolate the current field to the present time.
!
!  2. Method :
!
!     Linear interpolation of speed and direction, with optionally
!     a correction to get approximate quadratic interpolation of speed
!     only.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       FLFRST  Log.  I   Flag for first pass through routine.
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
!      W3WAVE    Subr. W3WAVEMD Actual wave model routine.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     None.
!
!  7. Remarks :
!
!     - Only currents at sea points are considered.
!     - Time ranges checked in W3WAVE.
!     - Currents are stored by components to save on the use of
!       SIN and COS functions. The actual interpolations, however
!       are by absolute value and direction.
!
!  8. Structure :
!
!     --------------------------------------
!      1.  Prepare auxiliary arrays.
!      2.  Calculate interpolation factors.
!      3.  Get actual winds.
!     --------------------------------------
!
!  9. Switches :
!
!     !/CRT1  Linear current interpolation.
!     !/CRT2  Quasi-quadratic current interpolation.
!
!     !/S     Enable subroutine tracing.
!     !/T     Test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3GDATMD, ONLY: NX, NY, NSEA, MAPSF
      USE W3WDATMD, ONLY: TIME
      USE W3ADATMD, ONLY: CX, CY, CA0, CAI, CD0, CDI
      USE W3IDATMD, ONLY: TC0, CX0, CY0, TCN, CXN, CYN
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      LOGICAL, INTENT(IN)     :: FLFRST
!/
!/ ------------------------------------------------------------------- /
!/
      INTEGER                 :: ISEA, IX, IY
      REAL                    :: D0, DN, DD, DTT, DT0, RD, CABS, CDIR
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Prepare auxiliary arrays
!
      IF ( FLFRST ) THEN
          DO ISEA=1, NSEA
            IX        = MAPSF(ISEA,1)
            IY        = MAPSF(ISEA,2)
            CA0(ISEA) = SQRT ( CX0(IX,IY)**2 + CY0(IX,IY)**2 )
            CAI(ISEA) = SQRT ( CXN(IX,IY)**2 + CYN(IX,IY)**2 )
            IF ( CA0(ISEA) .GT. 1.E-7) THEN
                D0     = MOD ( TPI+ATAN2(CY0(IX,IY),CX0(IX,IY)) , TPI )
              ELSE
                D0     = 0
              END IF
            IF ( CAI(ISEA) .GT. 1.E-7) THEN
                DN     = MOD ( TPI+ATAN2(CYN(IX,IY),CXN(IX,IY)) , TPI )
              ELSE
                DN     = D0
              END IF
            IF ( CA0(ISEA) .GT. 1.E-7) THEN
                CD0(ISEA) = D0
              ELSE
                CD0(ISEA) = DN
              END IF
            DD     = DN - CD0(ISEA)
            IF (ABS(DD).GT.PI) DD = DD - TPI*SIGN(1.,DD)
            CDI(ISEA) = DD
            CAI(ISEA) = CAI(ISEA) - CA0(ISEA)
            END DO
        END IF
!
! 2.  Calculate interpolation factor
!
      DTT    = DSEC21 ( TC0, TCN )
      DT0    = DSEC21 ( TC0, TIME )
!
      RD     = DT0 / MAX ( 1.E-7 , DTT )
!
! 3.  Actual currents for all grid points
!
      DO ISEA=1, NSEA
!
        CABS    = CA0(ISEA) + RD * CAI(ISEA)
        CDIR    = CD0(ISEA) + RD * CDI(ISEA)
        CX(ISEA) = CABS * COS(CDIR)
        CY(ISEA) = CABS * SIN(CDIR)
!
        END DO
!
      RETURN
!
! Formats
!
!/
!/ End of W3UCUR ----------------------------------------------------- /
!/
      END SUBROUTINE W3UCUR
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3UWND ( FLFRST, VGX, VGY )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         04-Jul-2006 |
!/                  +-----------------------------------+
!/
!/    03-Dec-1998 : Final FORTRAN 77                    ( version 1.18 )
!/    20-Dec-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    13-Nov-2002 : Add stress vector.                  ( version 3.00 )
!/    26-Dec-2002 : Moving grid wind correction.        ( version 3.02 )
!/    15-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    04-Jul-2006 : Consolidate stress arrays.          ( version 3.09 )
!/
!  1. Purpose :
!
!     Interpolate wind fields to the given time.
!
!  2. Method :
!
!     Linear interpolation of wind speed and direction, with a simple
!     correction to obtain quasi-conservation of energy.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       FLFRST  Log.  I   Flag for first pass through routine.
!       VGX/Y   Real  I   Grid velocity                      (!/MGW)
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
!      W3WAVE    Subr. W3WAVEMD Actual wave model routine.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!     - Only winds over sea points are considered.
!     - Time ranges checked in W3WAVE.
!
!  8. Structure :
!
!     --------------------------------------
!      1.  Prepare auxiliary arrays.
!      2.  Calculate interpolation factors
!      3.  Get actual winds
!      4.  Correct for currents
!      5.  Convert to stresses
!      6.  Stability correction
!     --------------------------------------
!
!  9. Switches :
!
!     !/OMP1   OpenMP compiler directives.
!
!     !/WNT0   No wind interpolation.
!     !/WNT1   Linear wind interpolation.
!     !/WNT2   Energy conservation in wind interpolation.
!
!     !/RWND   Use wind speeds relative to currents.
!     !/MGW    Moving grid wind correction.
!
!     !/STAB2  Calculate effective wind speed factor for stability
!              to be used with !/ST2.
!
!     !/S      Enable subroutine tracing.
!     !/T      Test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3GDATMD, ONLY: NX, NY, NSEA, MAPSF
      USE W3WDATMD, ONLY: TIME, ASF
      USE W3ADATMD, ONLY: DW, CX, CY, UA, UD, U10, U10D, AS,          &
                          UA0, UAI, UD0, UDI, AS0, ASI
      USE W3IDATMD, ONLY: TW0, WX0, WY0, DT0, TWN, WXN, WYN, DTN, FLCUR
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: VGX, VGY
      LOGICAL, INTENT(IN)     :: FLFRST
!/
!/ ------------------------------------------------------------------- /
!/
      INTEGER                 :: ISEA, IX, IY
      REAL                    :: D0, DN, DD, DTT, DT00, RD, UI2,      &
                                 UXR, UYR
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Prepare auxiliary arrays
!
      IF ( FLFRST ) THEN
          DO ISEA=1, NSEA
            IX        = MAPSF(ISEA,1)
            IY        = MAPSF(ISEA,2)
            UA0(ISEA) = SQRT ( WX0(IX,IY)**2 + WY0(IX,IY)**2 )
            UAI(ISEA) = SQRT ( WXN(IX,IY)**2 + WYN(IX,IY)**2 )
            IF ( UA0(ISEA) .GT. 1.E-7) THEN
                D0     = MOD ( TPI+ATAN2(WY0(IX,IY),WX0(IX,IY)) , TPI )
              ELSE
                D0     = 0
              END IF
            IF ( UAI(ISEA) .GT. 1.E-7) THEN
                DN     = MOD ( TPI+ATAN2(WYN(IX,IY),WXN(IX,IY)) , TPI )
              ELSE
                DN     = D0
              END IF
            IF ( UA0(ISEA) .GT. 1.E-7) THEN
                UD0(ISEA) = D0
              ELSE
                UD0(ISEA) = DN
              END IF
            DD     = DN - UD0(ISEA)
            IF (ABS(DD).GT.PI) DD = DD - TPI*SIGN(1.,DD)
            UDI(ISEA) = DD
            UAI(ISEA) = UAI(ISEA) - UA0(ISEA)
            AS0(ISEA) = DT0(IX,IY)
            ASI(ISEA) = DTN(IX,IY) - DT0(IX,IY)
            END DO
        END IF
!
! 2.  Calculate interpolation factor
!
      DTT    = DSEC21 ( TW0, TWN )
      DT00   = DSEC21 ( TW0, TIME )
!
      RD     = DT00 / MAX ( 1.E-7 , DTT )
!
! 3.  Actual wind for all grid points
!
      DO ISEA=1, NSEA
!
        UA(ISEA) = UA0(ISEA) + RD * UAI(ISEA)
        UD(ISEA) = UD0(ISEA) + RD * UDI(ISEA)
        AS(ISEA) = AS0(ISEA) + RD * ASI(ISEA)
!
        END DO
!
! 4.  Correct for currents and grid motion
!
      IF ( FLCUR ) THEN
!
          DO ISEA=1, NSEA
            UXR        = UA(ISEA)*COS(UD(ISEA)) - CX(ISEA)
            UYR        = UA(ISEA)*SIN(UD(ISEA)) - CY(ISEA)
            U10 (ISEA) = MAX ( 0.001 , SQRT(UXR**2+UYR**2) )
            U10D(ISEA) = MOD ( TPI+ATAN2(UYR,UXR) , TPI )
            END DO
!
        ELSE
!
          DO ISEA=1, NSEA
            U10 (ISEA) = MAX ( UA(ISEA) , 0.001 )
            U10D(ISEA) = UD(ISEA)
            END DO
!
        END IF
!
! 5.  Stability correction ( !/STAB2 )
!     Original settings :
!
!     SHSTAB =    1.4
!     OFSTAB =   -0.01
!     CCNG   =   -0.1
!     CCPS   =    0.1
!     FFNG   = -150.
!     FFPS   =  150.
!
      RETURN
!
! Formats
!
!/
!/ End of W3UWND ----------------------------------------------------- /
!/
      END SUBROUTINE W3UWND
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3UINI ( A )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         15-Dec-2004 |
!/                  +-----------------------------------+
!/
!/    19-Oct-1998 : Final FORTRAN 77                    ( version 1.18 )
!/    20-Dec-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    24-Jan-2001 : Flat grid version.                  ( version 2.06 )
!/    18-May-2001 : Fix CG declaration.                 ( version 2.11 )
!/    15-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/
!  1. Purpose :
!
!     Initialize the wave field with fetch-limited spectra before the
!     actual calculation start. (Named as an update routine due to
!     placement in code.)
!
!  2. Method :
!
!     Fetch-limited JONSWAP spectra with a cosine^2 directional
!     distribution and a mean direction taken from the wind.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       A       R.A.   O   Action density spectra.
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
!      W3WAVE    Subr. W3WAVEMD Actual wave model routine.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!     - Wind speeds filtered by U10MIN and U10MAX (DATA statements)
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!       !/SHRD  Switch for shared / distributed memory architecture.
!       !/DIST  Id.
!
!       !/LLG    Spherical grid.
!       !/XYG    Cartesian grid.
!
!       !/S     Enable subroutine tracing.
!       !/T     General test output.
!       !/T1    Parameters at grid points.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3GDATMD, ONLY: NX, NY, NSEA, NSEAL, MAPSF, SX, SY,         &
                          NK, NTH, TH, SIG, DTH, DSIP
      USE W3ADATMD, ONLY: U10, U10D, CG
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(OUT)       :: A(NTH,NK,0:NSEAL)
!/
!/ ------------------------------------------------------------------- /
!/ Local variables
!/
      INTEGER                 :: ISEA, JSEA, IK, ITH
      REAL                    :: ALFA(NSEAL), FP(NSEAL), YLN(NSEAL),  &
                                 AA, BB, CC
      REAL                    :: XGR, U10C, U10DIR, XSTAR, FSTAR,     &
                                 GAMMA, FR, D1(NTH), D1INT, F1, F2
      REAL                    :: ETOT, E1I
      REAL                    :: U10MIN =  1.
      REAL                    :: U10MAX = 20.
!/
!/ ------------------------------------------------------------------- /
!/
!
      XGR    = 0.5 * RADIUS * SQRT(SX**2+SY**2) * DERA
!
! Pre-process JONSWAP data for all grid points ----------------------- *
!
      DO JSEA=1, NSEAL
        ISEA   = IAPROC + (JSEA-1)*NAPROC
!
        U10C   = MAX ( MIN(U10(ISEA),U10MAX) , U10MIN )
!
        XSTAR  = GRAV * XGR / U10C**2
        FSTAR  = 3.5 / XSTAR**(0.33)
        GAMMA  = MAX ( 1. , 7.0 / XSTAR**(0.143) )
!
        ALFA(JSEA) = 0.076 / XSTAR**(0.22)
        FP  (JSEA) = FSTAR * GRAV / U10C
        YLN (JSEA) = LOG ( GAMMA )
!
        END DO
!
! 1-D spectrum at location ITH = NTH --------------------------------- *
!
      DO IK=1, NK
        FR     = SIG(IK) * TPIINV
        DO JSEA=1, NSEAL
!
!/ ----- INLINED EJ5P (REDUCED) -------------------------------------- /
!
          AA     = ALFA(JSEA) * 0.06175/FR**5
          BB     = MAX( -50. , -1.25*(FP(JSEA)/FR)**4 )
          CC     = MAX( -50. , -0.5*((FR-FP(JSEA))/(0.07*FP(JSEA)))**2 )
          A(NTH,IK,JSEA)                                              &
                 = AA * EXP(BB + EXP(CC) * YLN(JSEA))
!
!/ ----- INLINED EJ5P (END) ------------------------------------------ /
!
          END DO
        END DO
!
! Apply directional distribution ------------------------------------- *
!
      DO JSEA=1, NSEAL
        ISEA   = IAPROC + (JSEA-1)*NAPROC
!
        U10DIR = U10D(ISEA)
        D1INT  = 0.
!
        DO ITH=1, NTH
          D1(ITH) = ( MAX ( 0. , COS(TH(ITH)-U10DIR) ) )**2
          D1INT   = D1INT + D1(ITH)
          END DO
!
        D1INT  = D1INT * DTH
        F1     = TPIINV / D1INT
!
        DO IK=1, NK
          F2     = F1 * A(NTH,IK,JSEA) * CG(IK,ISEA) / SIG(IK)
          DO ITH=1, NTH
            A(ITH,IK,JSEA) = F2 * D1(ITH)
            END DO
          END DO
!
        END DO
!
! Test output -------------------------------------------------------- *
!
      RETURN
!
! Formats
!
!/
!/ End of W3UINI ----------------------------------------------------- /
!/
      END SUBROUTINE W3UINI
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3UBPT
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         07-Sep-2005 |
!/                  +-----------------------------------+
!/
!/    19-Oct-1998 : Final FORTRAN 77                    ( version 1.18 )
!/    20-Dec-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    15-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    07-Sep-2005 : Moving update to end of time step.  ( version 3.08 )
!/
!  1. Purpose :
!
!     Update spectra at the active boundary points.
!
!  2. Method :
!
!     Spectra are read and interpolated in space and time from the
!     data read by W3IOBC.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
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
!      W3WAVE    Subr. W3WAVEMD Actual wave model routine.
!     ----------------------------------------------------------------
!       STRACE, DSEC21
!                Service routines.
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!     - The data arrays contain sigma spectra to assure conservation
!       when changing grids.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S     Enable subroutine tracing.
!     !/T0    Test output of wave heights.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3GDATMD, ONLY: NSPEC, MAPWN, SIG2, DDEN
      USE W3ADATMD, ONLY: CG
      USE W3ODATMD, ONLY: NBI, ABPI0, ABPIN, ISBPI, IPBPI, RDBPI,     &
                          BBPI0, BBPIN
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/ ------------------------------------------------------------------- /
!/
      INTEGER                 :: IBI, ISP, ISEA
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Process BBPI0 -------------------------------------------------- *
! 1.a First intialization
!
      IF ( BBPI0(1,0) .EQ. -1. ) THEN
!
          BBPI0(:,0) = 0.
          BBPIN(:,0) = 0.
          ABPI0(:,0) = 0.
          ABPIN(:,0) = 0.
!
          DO IBI=1, NBI
            ISEA   = ISBPI(IBI)
            DO ISP=1, NSPEC
              BBPI0(ISP,IBI) = CG(MAPWN(ISP),ISEA) / SIG2(ISP) *      &
                           ( RDBPI(IBI,1) * ABPI0(ISP,IPBPI(IBI,1))   &
                           + RDBPI(IBI,2) * ABPI0(ISP,IPBPI(IBI,2))   &
                           + RDBPI(IBI,3) * ABPI0(ISP,IPBPI(IBI,3))   &
                           + RDBPI(IBI,4) * ABPI0(ISP,IPBPI(IBI,4)) )
              END DO
            END DO
!
! 1.b Shift BBPIN
!
        ELSE
          BBPI0 = BBPIN
        END IF
!
! 2.  Process BBPIN -------------------------------------------------- *
!
      DO IBI=1, NBI
        ISEA   = ISBPI(IBI)
        DO ISP=1, NSPEC
          BBPIN(ISP,IBI) = CG(MAPWN(ISP),ISEA) / SIG2(ISP) *          &
                       ( RDBPI(IBI,1) * ABPIN(ISP,IPBPI(IBI,1))       &
                       + RDBPI(IBI,2) * ABPIN(ISP,IPBPI(IBI,2))       &
                       + RDBPI(IBI,3) * ABPIN(ISP,IPBPI(IBI,3))       &
                       + RDBPI(IBI,4) * ABPIN(ISP,IPBPI(IBI,4)) )
          END DO
        END DO
!
! 3.  Wave height test output ---------------------------------------- *
!
      RETURN
!
! Formats
!
 
!/
!/ End of W3UBPT ----------------------------------------------------- /
!/
      END SUBROUTINE W3UBPT
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3UICE ( A, VA )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         11-May-2007 |
!/                  +-----------------------------------+
!/
!/    19-Oct-1998 : Final FORTRAN 77                    ( version 1.18 )
!/    20-Dec-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    11-Jan-2002 : Sub-grid ice.                       ( version 2.15 )
!/    15-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    28-Jun-2005 : Adding MAPST2.                      ( version 3.07 )
!/                  Taking out initilization.
!/    11-May-2007 : Adding NTPROC/NAPROC separation.    ( version 3.11 )
!/
!  1. Purpose :
!
!     Update ice map in the wave model.
!
!  2. Method :
!
!     Points with an ice concentration larger than FICEN are removed
!     from the sea map in the wave model. Such points are identified
!     by negative numbers is the grid status map MAPSTA. For ice
!     points spectra are set to zero. Points from wich ice disappears
!     are initialized with a "small" JONSWAP spectrum, based on the
!     frequency SIG(NK-1) and the local wind direction.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!      (V)A     R.A.  I/O   Spectra in 1-D or 2-D representation
!                           (points to same address).
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
!      W3WAVE    Subr. W3WAVEMD Actual wave model routine.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/SHRD  Switch for shared / distributed memory architecture.
!     !/DIST  Id.
!
!     !/S  Enable subroutine tracing.
!     !/T  Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3GDATMD, ONLY: NX, NY, NSEA, NSEAL, MAPSF, MAPSTA, MAPST2, &
                          NTH, NK, NSPEC, SIG, TH, DTH, FICEN
      USE W3WDATMD, ONLY: TIME, TICE, ICE
!!    USE W3ADATMD, ONLY: U10, U10D, CG
      USE W3ADATMD, ONLY: CG
      USE W3IDATMD, ONLY: TIN, ICEI
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(INOUT)     :: A(NTH,NK,0:NSEAL), VA(NSPEC,0:NSEAL)
!/
!/ ------------------------------------------------------------------- /
!/
      INTEGER                 :: IK, ITH, ISEA, JSEA, IX, IY, ISP
      INTEGER                 :: MAPICE(NY,NX)
      LOGICAL                 :: LOCAL
!/
!/ ------------------------------------------------------------------- /
!/
!
      LOCAL   = IAPROC .LE. NAPROC
!
! 1.  Preparations --------------------------------------------------- *
! 1.a Update times
!
      TICE(1) = TIN(1)
      TICE(2) = TIN(2)
!
! 1.b Process maps
!
      MAPICE = MOD(MAPST2,2)
      MAPST2 = MAPST2 - MAPICE
!
! 2.  Main loop over sea points -------------------------------------- *
!
      DO ISEA=1, NSEA
!
! 2.a Get grid counters
!
        IX        = MAPSF(ISEA,1)
        IY        = MAPSF(ISEA,2)
        ICE(ISEA) = ICEI(IX,IY)
!
! 2.b Sea point to be de-activated..
!
        IF ( ICEI(IX,IY).GE.FICEN .AND. MAPICE(IY,IX).EQ.0 ) THEN
            MAPSTA(IY,IX) = - ABS(MAPSTA(IY,IX))
            MAPICE(IY,IX) = 1
            IF ( LOCAL .AND. MOD(ISEA-IAPROC,NAPROC) .EQ. 0 ) THEN
                JSEA   = 1 + (ISEA-1)/NAPROC
                VA(:,JSEA) = 0.
              END IF
!
          END IF
!
! 2.b Ice point to be re-activated.
!
        IF ( ICEI(IX,IY).LT.FICEN .AND. MAPICE(IY,IX).EQ.1 ) THEN
!
            MAPICE(IY,IX) = 0
!
            IF ( MAPST2(IY,IX) .EQ. 0 ) THEN
                MAPSTA(IY,IX) = ABS(MAPSTA(IY,IX))
!
                IF ( LOCAL .AND. MOD(ISEA-IAPROC,NAPROC) .EQ. 0 ) THEN
                    JSEA   = 1 + (ISEA-1)/NAPROC
                    VA(:,JSEA) = 0.
!
                  END IF
!
              END IF
!
          END IF
!
        END DO
!
! 3.  Update MAPST2 -------------------------------------------------- *
!
      MAPST2 = MAPST2 + MAPICE
!
      RETURN
!
! Formats
!
!/
!/ End of W3UICE ----------------------------------------------------- /
!/
      END SUBROUTINE W3UICE
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3ULEV ( A, VA )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         11-May-2007 |
!/                  +-----------------------------------+
!/
!/    15-Jan-1998 : Final FORTRAN 77                    ( version 1.18 )
!/    21-Jan-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    30-Apr-2002 : Water level fixes.                  ( version 2.20 )
!/    15-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    15-Jul-2005 : Adding drying out of points.        ( version 3.07 )
!/    11-May-2007 : Adding NTPROC/NAPROC separation.    ( version 3.11 )
!/
!  1. Purpose :
!
!     Update the water level.
!
!  2. Method :
!
!     The wavenumber grid is modified without modyfying the spectrum
!     (conservative linear interpolation to new grid).
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!      (V)A     R.A.  I/O  2-D and 1-D represetation of the spectra.
!     ----------------------------------------------------------------
!
!     Local variables
!     ----------------------------------------------------------------
!       KDMAX   Real   Deep water cut-off for kd.
!       WNO     R.A.   Old wavenumbers.
!       CGO     R.A.   Old group velocities.
!       OWN     R.A.   Old wavenumber band width.
!       DWN     R.A.   New wavenumber band width.
!       TA      R.A.   Auxiliary spectrum.
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
!      W3WAVE    Subr. W3WAVEMD Actual wave model routine.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     None.
!
!  7. Remarks :
!
!     - The grid is updated only if KDmin > KDMAX.
!     - The grid is updated for inactive points too.
!     - The local wavenumber bandwidth is DSIGMA/CG.
!     - The local spectrum is updated only if the grid is updated,
!       the grid point is not disabled (MAPST2) and if the change of
!       the lowest wavenumber exceeds RDKMIN times the band width.
!     - No spectral initialization for newly wet points.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S     Enable subroutine tracing.
!     !/T     Basic test output.
!     !/T2    Output of minimum relative depth per grid point.
!     !/T3    Spectra before and after
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3GDATMD, ONLY: NX, NY, NSEA, NSEAL, MAPSF, MAPSTA, MAPST2, &
                          ZB, DMIN, NK, NTH, NSPEC, SIG, DSIP,        &
                          MAPWN, MAPTH, FACHFA
      USE W3WDATMD, ONLY: TIME, TLEV, WLV
      USE W3ADATMD, ONLY: CG, WN, DW
      USE W3IDATMD, ONLY: TLN, WLEV
      USE W3SERVMD, ONLY: EXTCDE
      USE W3DISPMD, ONLY: NAR1D, DFAC, N1MAX, ECG1, EWN1, DSIE
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(INOUT)     :: A(NTH,NK,0:NSEAL), VA(NSPEC,0:NSEAL)
!/
!/ ------------------------------------------------------------------- /
!/
      INTEGER                 :: ISEA, JSEA, IX, IY, IK, I1, I2,      &
                                 ISPEC, IK0, ITH
      INTEGER                 :: MAPDRY(NY,NX)
      REAL                    :: DWO(NSEA), KDCHCK, WNO(0:NK+1),      &
                                 CGO(0:NK+1), DEPTH, SQRTH, SIX,      &
                                 R1, R2, RDK, RD1, RD2, TA(NTH,NK),   &
                                 OWN(NK), DWN(NK)
      REAL                    :: KDMAX = 4., RDKMIN = 0.05
      LOGICAL                 :: LOCAL
!/
!/ ------------------------------------------------------------------- /
!/
!
      LOCAL   = IAPROC .LE. NAPROC
!
! 1.  Preparations --------------------------------------------------- *
! 1.a Check NK
!
      IF ( NK .LT. 2 ) THEN
          IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1000)
          CALL EXTCDE ( 1 )
        END IF
!
! 1.b Update times
!
      TLEV = TLN
!
! 1.c Extract dry point map, and residual MAPST2
!
      MAPDRY = MOD(MAPST2/2,2)
      MAPST2 = MAPST2 - 2*MAPDRY
!
! 1.d Update water levels and save old
!
      DO ISEA=1, NSEA
        IX     = MAPSF(ISEA,1)
        IY     = MAPSF(ISEA,2)
        DWO(ISEA) = DW(ISEA)
        WLV(ISEA) = WLEV(IX,IY)
        DW (ISEA) = MAX ( 0. , WLV(ISEA)-ZB(ISEA) )
        END DO
!
! 2.  Loop over all sea points --------------------------------------- *
!
      DO ISEA=1, NSEA
!
        IX     = MAPSF(ISEA,1)
        IY     = MAPSF(ISEA,2)
!
! 2.a Check if deep water
!
        KDCHCK = WN(1,ISEA) * MIN( DWO(ISEA) , DW(ISEA) )
        IF ( KDCHCK .LT. KDMAX ) THEN
!
! 2.b Update grid and save old grid
!
            DEPTH  = MAX ( DMIN, DW(ISEA) )
            SQRTH  = SQRT(DEPTH)
!
            DO IK=0, NK+1
              WNO(IK) = WN(IK,ISEA)
              CGO(IK) = CG(IK,ISEA)
!
!/ --- INLINED WAVNU1 (START) ---------------------------------------- /
!
              SIX    = SIG(IK) * SQRTH
              I1     = INT(SIX/DSIE)
              IF (I1.LE.N1MAX) THEN
                  I2 = I1 + 1
                  R1 = SIX/DSIE - REAL(I1)
                  R2 = 1. - R1
                  WN(IK,ISEA) = ( R2*EWN1(I1) + R1*EWN1(I2) ) / DEPTH
                  CG(IK,ISEA) = ( R2*ECG1(I1) + R1*ECG1(I2) ) * SQRTH
                ELSE
                  WN(IK,ISEA) = SIG(IK)*SIG(IK)/GRAV
                  CG(IK,ISEA) = 0.5 * GRAV / SIG(IK)
                END IF
!
!/ --- INLINED WAVNU1 (END) ------------------------------------------ /
!
              END DO
!
            DO IK=1, NK
              OWN(IK) = DSIP(IK) / CGO(IK)
              DWN(IK) = DSIP(IK) / CG(IK,ISEA)
              END DO
!
! 2.c Process dry points
!
            IF ( WLV(ISEA)-ZB(ISEA) .LE.0. ) THEN
                IF ( MAPDRY(IY,IX) .EQ. 0 ) THEN
                    IF ( LOCAL .AND. MOD(ISEA-IAPROC,NAPROC).EQ.0 ) THEN
                        JSEA   = 1 + (ISEA-1)/NAPROC
                        VA(:,JSEA) = 0.
                      END IF
                    MAPDRY(IY,IX) = 1
                    MAPSTA(IY,IX) = -ABS(MAPSTA(IY,IX))
                  ENDIF
                CYCLE
             END IF
!
! 2.d Process new wet point
!
            IF (WLV(ISEA)-ZB(ISEA).GT.0. .AND. MAPDRY(IY,IX).EQ.1) THEN
                MAPDRY(IY,IX) = 0
                IF ( MAPST2(IY,IX) .EQ. 0 ) THEN
                    MAPSTA(IY,IX) = ABS(MAPSTA(IY,IX))
                  END IF
                CYCLE
             END IF
!
! 2.e Check if ice on grid point, or if grid changes negligible
!
            RDK    = ABS(WNO(1)-WN(1,ISEA)) / DWN(1)
!
            IF ( RDK.LT.RDKMIN .OR. MAPSTA(IY,IX).LT.0 ) CYCLE
            IF ( MOD(ISEA-IAPROC,NAPROC) .NE. 0 ) CYCLE
            JSEA   = 1 + (ISEA-1)/NAPROC
!
            IF ( .NOT. LOCAL ) CYCLE
!
! 2.d Save discrete actions and clean spectrum
!
            DO IK=1, NK
              DO ITH=1, NTH
                TA(ITH,IK) = A(ITH,IK,JSEA) * OWN(IK)
                END DO
              END DO
!
            VA(:,JSEA) = 0.
!
! 2.e Redistribute discrete action density
!
            IF ( WNO(1) .LT. WN(1,ISEA) ) THEN
                IK0    = 1
                I1     = 0
                I2     = 1
  220           CONTINUE
                IK0    = IK0 + 1
                IF ( IK0 .GT. NK+1 ) GOTO 251
                IF ( WNO(IK0) .GE. WN(1,ISEA) ) THEN
                     IK0    = IK0 - 1
                   ELSE
                     GOTO 220
                   END IF
              ELSE
                IK0    = 1
                I1     = 1
                I2     = 2
              END IF
!
            DO 250, IK=IK0, NK
!
  230         CONTINUE
              IF ( WNO(IK) .GT. WN(I2,ISEA) ) THEN
                  I1     = I1 + 1
                  IF ( I1 .GT. NK ) GOTO 250
                  I2     = I1 + 1
                  GOTO 230
                END IF
!
              IF ( I1 .EQ. 0 ) THEN
                  RD1    = ( WN(1,ISEA) - WNO(IK) ) / DWN(1)
                  RD2    = 1. - RD1
                ELSE
                  RD1    = ( WN(I2,ISEA) - WNO(IK) ) /                &
                           ( WN(I2,ISEA) - WN(I1,ISEA) )
                  RD2    = 1. - RD1
                END IF
!
                IF ( I1 .GE. 1 ) THEN
                    DO ITH=1, NTH
                      A(ITH,I1,JSEA) = A(ITH,I1,JSEA) + RD1*TA(ITH,IK)
                      END DO
                  END IF
!
                IF ( I2 .LE. NK ) THEN
                    DO ITH=1, NTH
                      A(ITH,I2,JSEA) = A(ITH,I2,JSEA) + RD2*TA(ITH,IK)
                      END DO
                  END IF
!
  250         CONTINUE
  251       CONTINUE
!
! 2.f Convert discrete action densities to spectrum
!
            DO ISPEC=1, NSPEC
              VA(ISPEC,JSEA) = VA(ISPEC,JSEA) / DWN(MAPWN(ISPEC))
              END DO
!
! 2.f Add tail if necessary
!
            IF ( I2.LE.NK .AND. RD2.LE.0.95 ) THEN
                DO IK=MAX(I2,2), NK
                  DO ITH=1, NTH
                    A(ITH,IK,JSEA) = FACHFA * A(ITH,IK-1,JSEA)
                    END DO
                  END DO
              END IF
!
          END IF
!
        END DO
!
! 3.  Reconstruct new MAPST2 ----------------------------------------- *
!
      MAPST2 = MAPST2 + 2*MAPDRY
!
      RETURN
!
! Formats
!
 1000 FORMAT (/' *** ERROR W3ULEV *** '/                              &
               '     THIS ROUTINE REQUIRES NK > 1 '/)
!
!/
!/ End of W3ULEV ----------------------------------------------------- /
!/
      END SUBROUTINE W3ULEV
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3UTRN ( TRNX, TRNY )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         11-Jan-2007 |
!/                  +-----------------------------------+
!/
!/    02-Apr-2001 : Origination.                        ( version 2.10 )
!/    11-Jan-2002 : Sub-grid ice.                       ( version 2.15 )
!/    30-Apr-2002 : Change to ICE on storage grid.      ( version 2.20 )
!/    15-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    11-Jan-2007 : Clean-up for boundary points.       ( version 3.10 )
!/
!  1. Purpose :
!
!     Update cell boundary transparencies for general use in propagation
!     routines.
!
!  2. Method :
!
!     Two arrays are generated with the size (NY*NX,-1:1). The value
!     at (IXY,-1) indicates the transparency to be used if the lower
!     or left boundary is an inflow boundary. (IXY,1) is used if the
!     upper or right boundary is an inflow boundary. (IXY,0) is used
!     for all other cases (by definition full transparency).
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       TRNX/Y  R.A.   I   Transparencies from model defintion file.
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
!      W3WAVE    Subr. W3WAVEMD Actual wave model routine.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!       !/S      Enable subroutine tracing.
!       !/T      Basic test output.
!
!       !/LLG    Spherical grid.
!       !/XYG    Cartesian grid.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3GDATMD, ONLY: NX, NY, NSEA, MAPSTA, MAPSF, SX, SY,        &
                          TRFLAG, CLAT, FICE0, FICEN
      USE W3WDATMD, ONLY: ICE
      USE W3ADATMD, ONLY: ATRNX, ATRNY
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: TRNX(NY*NX), TRNY(NY*NX)
!/
!/ ------------------------------------------------------------------- /
!/
      INTEGER                 :: ISEA, IX, IY, IXY, IXN, IXP, IYN, IYP
 
      REAL                    :: TRIX(NY*NX), TRIY(NY*NX), DX, DY,    &
                                 LICE0, LICEN
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Preparations --------------------------------------------------- *
!
      ATRNX = 1.
      ATRNY = 1.
!
! 2.  Filling arrays from TRNX/Y for obstructions -------------------- *
! 2.a TRFLAG = 0, no action needed
!
      IF ( TRFLAG .EQ. 0 ) THEN
          RETURN
!
! 2.b TRFLAG = 1,3: TRNX/Y defined at boundaries
!
        ELSE IF ( TRFLAG.EQ.1 .OR. TRFLAG.EQ.3 ) THEN
!
          DO ISEA=1, NSEA
!
            IX            = MAPSF(ISEA,1)
            IY            = MAPSF(ISEA,2)
            IXY           = MAPSF(ISEA,3)
            IF ( IX .EQ. 1 ) THEN
                ATRNX(IXY,-1) = TRNX(IY+(NX-1)*NY)
                ATRNX(IXY, 1) = TRNX(IXY)
              ELSE IF ( IX .EQ. NX ) THEN
                ATRNX(IXY,-1) = TRNX(IXY-NY)
                ATRNX(IXY, 1) = TRNX(IY)
              ELSE
                ATRNX(IXY,-1) = TRNX(IXY-NY)
                ATRNX(IXY, 1) = TRNX(IXY)
              END IF
            ATRNY(IXY,-1) = TRNY(IXY-1)
            ATRNY(IXY, 1) = TRNY(IXY)
!
            END DO
!
! 2.c TRFLAG = 2,4: TRNX/Y defined at cell centers
!
        ELSE
!
          DO ISEA=1, NSEA
!
            IX            = MAPSF(ISEA,1)
            IY            = MAPSF(ISEA,2)
            IXY           = MAPSF(ISEA,3)
!
            IF ( IX .EQ. 1 ) THEN
                IXN    = IY + (NX-1)*NY
                IXP    = IY +  IX   *NY
              ELSE IF ( IX .EQ. NX ) THEN
                IXN    = IY + (IX-2)*NY
                IXP    = IY
              ELSE
                IXN    = IY + (IX-2)*NY
                IXP    = IY +  IX   *NY
              END IF
!
            IF ( IY .EQ. 1 ) THEN
                IYN    = IXY
                IYP    = IXY + 1
              ELSE IF ( IY .EQ. NY ) THEN
                IYN    = IXY - 1
                IYP    = IXY
              ELSE
                IYN    = IXY - 1
                IYP    = IXY + 1
              END IF
!
! factors 0.5 in first term and 2. in second term cancel
!
            ATRNX(IXY,-1) = (1.+TRNX(IXY)) * TRNX(IXN)/(1.+TRNX(IXN))
            ATRNX(IXY, 1) = (1.+TRNX(IXY)) * TRNX(IXP)/(1.+TRNX(IXP))
            ATRNY(IXY,-1) = (1.+TRNY(IXY)) * TRNY(IYN)/(1.+TRNY(IYN))
            ATRNY(IXY, 1) = (1.+TRNY(IXY)) * TRNY(IYP)/(1.+TRNY(IYP))
!
            IF ( MAPSTA(IY,IX) .EQ. 2 ) THEN
                IF ( IX .EQ. 1  ) THEN
                    ATRNX(IXY,-1) = 1.
                  ELSE IF ( MAPSTA( IY ,IX-1) .LE. 0 ) THEN
                    ATRNX(IXY,-1) = 1.
                  END IF
                IF ( IX .EQ. NX ) THEN
                     ATRNX(IXY, 1) = 1.
                  ELSE IF ( MAPSTA( IY ,IX+1) .LE. 0 ) THEN
                     ATRNX(IXY, 1) = 1.
                  END IF
                IF ( IY .EQ. 1  ) THEN
                     ATRNY(IXY,-1) = 1.
                  ELSE IF ( MAPSTA(IY-1, IX ) .LE. 0 ) THEN
                     ATRNY(IXY,-1) = 1.
                  END IF
                IF ( IY .EQ. NY ) THEN
                     ATRNY(IXY, 1) = 1.
                  ELSE IF ( MAPSTA(IY+1, IX ) .LE. 0 ) THEN
                     ATRNY(IXY, 1) = 1.
                  END IF
              END IF
!
            END DO
        END IF
!
! 3.  Adding ice to obstructions ------------------------------------- *
! 3.a TRFLAG < 3, no action needed
!
      IF ( TRFLAG.LT.3 .OR. FICEN-FICE0.LT.1.E-6 ) THEN
          RETURN
!
! 3.b TRFLAG = 3,4: Calculate ice transparencies
!
        ELSE
          TRIX   = 1.
          TRIY   = 1.
!
          DO ISEA=1, NSEA
!
            IX     = MAPSF(ISEA,1)
            IY     = MAPSF(ISEA,2)
            IXY    = MAPSF(ISEA,3)
!
            DX     = SX * RADIUS * DERA * CLAT(IY)
            DY     = SY * RADIUS * DERA
!
            LICE0  = FICE0*MIN(DX,DY)
            LICEN  = FICEN*MIN(DX,DY)
!           TRIX(IXY) = MAX ( 0. , MIN ( 1. , 1.5-ICE(ISEA)*DX/LICE ) )
!           TRIY(IXY) = MAX ( 0. , MIN ( 1. , 1.5-ICE(ISEA)*DY/LICE ) )
            TRIX(IXY) = ( LICEN - ICE(ISEA)*DX ) / ( LICEN - LICE0 )
            TRIX(IXY) = MAX ( 0. , MIN ( 1. , TRIX(IXY) ) )
            TRIY(IXY) = ( LICEN - ICE(ISEA)*DY ) / ( LICEN - LICE0 )
            TRIY(IXY) = MAX ( 0. , MIN ( 1. , TRIY(IXY) ) )
!
          END DO
!
! 3.c Combine transparencies, ice always defined at cell center !
!
          DO ISEA=1, NSEA
!
            IX            = MAPSF(ISEA,1)
            IY            = MAPSF(ISEA,2)
            IXY           = MAPSF(ISEA,3)
!
            IF ( IX .EQ. 1 ) THEN
                IXN    = IY + (NX-1)*NY
                IXP    = IY +  IX   *NY
              ELSE IF ( IX .EQ. NX ) THEN
                IXN    = IY + (IX-2)*NY
                IXP    = IY
              ELSE
                IXN    = IY + (IX-2)*NY
                IXP    = IY +  IX   *NY
              END IF
!
            IF ( IY .EQ. 1 ) THEN
                IYN    = IXY
                IYP    = IXY + 1
              ELSE IF ( IY .EQ. NY ) THEN
                IYN    = IXY - 1
                IYP    = IXY
              ELSE
                IYN    = IXY - 1
                IYP    = IXY + 1
              END IF
!
            ATRNX(IXY,-1) = ATRNX(IXY,-1)                             &
                          * (1.+TRIX(IXY)) * TRIX(IXN)/(1.+TRIX(IXN))
            ATRNX(IXY, 1) = ATRNX(IXY, 1)                             &
                          * (1.+TRIX(IXY)) * TRIX(IXP)/(1.+TRIX(IXP))
            ATRNY(IXY,-1) = ATRNY(IXY,-1)                             &
                          * (1.+TRIY(IXY)) * TRIY(IYN)/(1.+TRIY(IYN))
            ATRNY(IXY, 1) = ATRNY(IXY, 1)                             &
                          * (1.+TRIY(IXY)) * TRIY(IYP)/(1.+TRIY(IYP))
!
            END DO
!
        END IF
!
      RETURN
!
! Formats
!
!/
!/ End of W3UTRN ----------------------------------------------------- /
!/
      END SUBROUTINE W3UTRN
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3DDXY
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         15-Dec-2004 |
!/                  +-----------------------------------+
!/
!/    17-Dec-1996 : Final FORTRAN 77                    ( version 1.18 )
!/    22-Dec-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    24-Jan-2001 : Flat grid version                   ( version 2.06 )
!/    15-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/
!  1. Purpose :
!
!     Calculate dirivatives of the depth field.
!
!  2. Method :
!
!     Derivaties are calculated in m/m from the longitude/latitude
!     grid, central in space for iternal points, one-sided for
!     coastal points.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
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
!      W3WAVE    Subr. W3WAVEMD Actual wave model routine.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  8. Structure :
!
!     ----------------------------------------
!      1.  Preparations
!        a Initialize arrays
!        b Set constants
!      2.  Derivatives in X-direction (W-E).
!      3.  Derivatives in Y-direction (S-N).
!     ----------------------------------------
!
!  9. Switches :
!
!       !/LLG    Spherical grid.
!       !/XYG    Cartesian grid.
!
!       !/S   Enable subroutine tracing.
!       !/T   Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3GDATMD, ONLY: NX, NY, NSEA, MAPSTA, MAPFS, SX, SY,        &
                          GLOBAL, CLATI
      USE W3ADATMD, ONLY: DW, DDDX, DDDY
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IX, IY, IXP, IXM, IYP, IYM
      REAL                    :: DX0I, DY0I, STX, STY
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Preparations --------------------------------------------------- *
! 1.a Initialize arrays
!
      DDDX = 0.
      DDDY = 0.
!
! 1.b Set constants
!
      DX0I   = 1. / ( SX * DERA * RADIUS )
      DY0I   = 1. / ( SY * DERA * RADIUS )
!
! 2.  Derivatives in X-direction (W-E) ------------------------------- *
!
      DO IY=2, NY-1
        DO IX=2, NX-1
          IF ( MAPSTA(IY,IX) .NE. 0 ) THEN
              STX    = 0.5
              IF (MAPSTA(IY,IX+1).EQ.0) THEN
                  IXP    = IX
                  STX    = 1.0
                ELSE
                  IXP    = IX + 1
                END IF
              IF (MAPSTA(IY,IX-1).EQ.0) THEN
                  IXM    = IX
                  STX    = 1.0
                ELSE
                  IXM    = IX - 1
                END IF
              DDDX(IY,IX) = STX * DX0I * CLATI(IY)                    &
                 * ( DW(MAPFS(IY,IXP)) - DW(MAPFS(IY,IXM)) )
            END IF
          END DO
        END DO
!
      IF ( GLOBAL ) THEN
!
          DO IY=2, NY-1
            IF ( MAPSTA(IY,1) .NE. 0 ) THEN
                STX    = 0.5
                IF (MAPSTA(IY,2).EQ.0) THEN
                    IXP    = 1
                    STX    = 1.0
                  ELSE
                    IXP    = 2
                  END IF
                IF (MAPSTA(IY,NX).EQ.0) THEN
                    IXM    = 1
                    STX    = 1.0
                  ELSE
                    IXM    = NX
                  END IF
                DDDX(IY,1) = STX * DX0I * CLATI(IY)                   &
                   * ( DW(MAPFS(IY,IXP)) - DW(MAPFS(IY,IXM)) )
              END IF
            END DO
!
          DO IY=2, NY-1
            IF ( MAPSTA(IY,NX) .NE. 0 ) THEN
                STX    = 0.5
                IF (MAPSTA(IY,1).EQ.0) THEN
                    IXP    = NX
                    STX    = 1.0
                  ELSE
                    IXP    = 1
                  END IF
                IF (MAPSTA(IY,NX-1).EQ.0) THEN
                    IXM    = NX
                    STX    = 1.0
                  ELSE
                    IXM    = NX-1
                  END IF
                DDDX(IY,NX) = STX * DX0I * CLATI(IY)                  &
                   * ( DW(MAPFS(IY,IXP)) - DW(MAPFS(IY,IXM)) )
              END IF
            END DO
!
        END IF
!
! 3.  Derivatives in Y-direction (S-N) ------------------------------- *
!
      DO IY=2, NY-1
        DO IX=1, NX
          IF ( MAPSTA(IY,IX) .NE. 0 ) THEN
              STY    = 0.5
              IF (MAPSTA(IY+1,IX).EQ.0) THEN
                  IYP    = IY
                  STY    = 1.0
                ELSE
                  IYP    = IY + 1
                END IF
              IF (MAPSTA(IY-1,IX).EQ.0) THEN
                  IYM    = IY
                  STY    = 1.0
                ELSE
                  IYM    = IY - 1
                END IF
!
              DDDY(IY,IX) = STY * DY0I                                &
                 * ( DW(MAPFS(IYP,IX)) - DW(MAPFS(IYM,IX)) )
!
            END IF
          END DO
        END DO
!
! 4.  Test output of fields ------------------------------------------ *
!
      RETURN
!
! Formats
!
!/
!/ End of W3DDXY ----------------------------------------------------- /
!/
      END SUBROUTINE W3DDXY
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3DCXY
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         15-Dec-2004 |
!/                  +-----------------------------------+
!/
!/    29-Aug-1997 : Final FORTRAN 77                    ( version 1.18 )
!/    22-Dec-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    24-Jan-2001 : Flat grid version                   ( version 2.06 )
!/    15-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/
!  1. Purpose :
!
!     Calculate dirivatives of the current field.
!
!  2. Method :
!
!     Derivaties are calculated in m/m from the longitude/latitude
!     grid, central in space for iternal points, one-sided for
!     coastal points.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
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
!      W3WAVE    Subr. W3WAVEMD Actual wave model routine.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  8. Structure :
!
!     ----------------------------------------
!      1.  Preparations
!        a Initialize arrays
!        b Set constants
!      2.  Derivatives in X-direction (W-E).
!      3.  Derivatives in Y-direction (S-N).
!     ----------------------------------------
!
!  9. Switches :
!
!       !/LLG    Spherical grid.
!       !/XYG    Cartesian grid.
!
!       !/S   Enable subroutine tracing.
!       !/T   Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3GDATMD, ONLY: NX, NY, NSEA, MAPSTA, MAPFS, SX, SY, CLATI, &
                          GLOBAL
      USE W3ADATMD, ONLY: CX, CY, DCXDX, DCXDY, DCYDX, DCYDY
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IX, IY, IXP, IXM, IYP, IYM
      INTEGER, SAVE           :: IPASS = 0
      REAL                    :: DX0I, DY0I, STX, STY
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Preparations --------------------------------------------------- *
! 1.a Initialize arrays
!
      IF ( IPASS .EQ. 0 ) THEN
          DCXDX = 0.
          DCXDY = 0.
          DCYDX = 0.
          DCYDY = 0.
          IPASS  = 1
        END IF
!
! 1.b Set constants
!
      DX0I   = 1. / ( SX * DERA * RADIUS )
      DY0I   = 1. / ( SY * DERA * RADIUS )
!
! 2.  Derivatives in X-direction (W-E) ------------------------------- *
!
      DO IY=2, NY-1
        DO IX=2, NX-1
          IF ( MAPSTA(IY,IX) .NE. 0 ) THEN
              STX    = 0.5
              IF (MAPSTA(IY,IX+1).EQ.0) THEN
                  IXP    = IX
                  STX    = 1.0
                ELSE
                  IXP    = IX + 1
                END IF
              IF (MAPSTA(IY,IX-1).EQ.0) THEN
                  IXM    = IX
                  STX    = 1.0
                ELSE
                  IXM    = IX - 1
                END IF
              DCXDX(IY,IX) = (CX(MAPFS(IY,IXP))-CX(MAPFS(IY,IXM))) *  &
                             STX * DX0I * CLATI(IY)
              DCYDX(IY,IX) = (CY(MAPFS(IY,IXP))-CY(MAPFS(IY,IXM))) *  &
                             STX * DX0I * CLATI(IY)
            END IF
          END DO
        END DO
!
      IF ( GLOBAL ) THEN
!
        DO IY=2, NY-1
          IF ( MAPSTA(IY,1) .NE. 0 ) THEN
              STX    = 0.5
              IF (MAPSTA(IY,2).EQ.0) THEN
                  IXP    = 1
                  STX    = 1.0
                ELSE
                  IXP    = 2
                END IF
              IF (MAPSTA(IY,NX).EQ.0) THEN
                  IXM    = 1
                  STX    = 1.0
                ELSE
                  IXM    = NX
                END IF
              DCXDX(IY,1) = (CX(MAPFS(IY,IXP))-CX(MAPFS(IY,IXM))) *   &
                            STX * DX0I * CLATI(IY)
              DCYDX(IY,1) = (CY(MAPFS(IY,IXP))-CY(MAPFS(IY,IXM))) *   &
                            STX * DX0I * CLATI(IY)
            END IF
          END DO
!
        DO IY=2, NY-1
          IF ( MAPSTA(IY,NX) .NE. 0 ) THEN
              STX    = 0.5
              IF (MAPSTA(IY,1).EQ.0) THEN
                  IXP    = NX
                  STX    = 1.0
                ELSE
                  IXP    = 1
                END IF
              IF (MAPSTA(IY,NX-1).EQ.0) THEN
                  IXM    = NX
                  STX    = 1.0
                ELSE
                  IXM    = NX-1
                END IF
              DCXDX(IY,NX) = (CX(MAPFS(IY,IXP))-CX(MAPFS(IY,IXM))) *  &
                             STX * DX0I * CLATI(IY)
              DCYDX(IY,NX) = (CY(MAPFS(IY,IXP))-CY(MAPFS(IY,IXM))) *  &
                             STX * DX0I * CLATI(IY)
            END IF
          END DO
!
        END IF
!
! 3.  Derivatives in Y-direction (S-N) ------------------------------- *
!
      DO IY=2, NY-1
        DO IX=1, NX
          IF ( MAPSTA(IY,IX) .NE. 0 ) THEN
              STY    = 0.5
              IF (MAPSTA(IY+1,IX).EQ.0) THEN
                  IYP    = IY
                  STY    = 1.0
                ELSE
                  IYP    = IY + 1
                END IF
              IF (MAPSTA(IY-1,IX).EQ.0) THEN
                  IYM    = IY
                  STY    = 1.0
                ELSE
                  IYM    = IY - 1
                END IF
!
              DCXDY(IY,IX) = (CX(MAPFS(IYP,IX))-CX(MAPFS(IYM,IX)))    &
                           * STY * DY0I
              DCYDY(IY,IX) = (CY(MAPFS(IYP,IX))-CY(MAPFS(IYM,IX)))    &
                           * STY * DY0I
!
            END IF
          END DO
        END DO
!
! 4.  Test output of fields ------------------------------------------ *
!
      RETURN
!
! Formats
!
!/
!/ End of W3DCXY ----------------------------------------------------- /
!/
      END SUBROUTINE W3DCXY
!/
!/ End of module W3UPDTMD -------------------------------------------- /
!/
      END MODULE W3UPDTMD
