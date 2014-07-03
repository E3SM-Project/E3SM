!/ ------------------------------------------------------------------- /
      MODULE W3PRO3MD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    27-Feb-2000 : Origination.                        ( version 2.08 )
!/    17-Sep-2000 : Clean-up.                           ( version 2.13 )
!/    10-Dec-2001 : Sub-grid obstructions.              ( version 2.14 )
!/    16-Oct-2002 : Change INTENT for ATRN in W3XYP3.   ( version 3.00 )
!/    26-Dec-2002 : Moving grid version.                ( version 3.02 )
!/    01-Aug-2003 : Moving grid GSE correction.         ( version 3.03 )
!/    17-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    07-Sep-2005 : Upgrade XY boundary conditions.     ( version 3.08 )
!/    09-Nov-2005 : Removing soft boundary option.      ( version 3.08 )
!/    05-Mar-2008 : Added NEC sxf90 compiler directives.
!/                  (Chris Bunney, UK Met Office)       ( version 3.13 )
!/    01-Apr-2008 : Bug fix W3MAP3 MAPSTA range check.  ( version 3.13 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Bundles routines for third order porpagation scheme in single
!     module.
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      TRNMIN    R.P.  Private   Minimum transparancy for local
!                                switching off of averaging.
!     ----------------------------------------------------------------
!
!     Also work arrays for W3KTP3 (private).
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3MAP3    Subr. Public   Set up auxiliary maps.
!      W3MAPT    Subr. Public   Set up transparency map for GSE.
!      W3XYP3    Subr. Public   Third order spatial propagation.
!      W3KTP3    Subr. Public   Third order spectral propagation.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!      W3QCK1    Subr. W3UQCKMD Regular grid UQ scheme.
!      W3QCK2    Subr. W3UQCKMD Irregular grid UQ scheme.
!      W3QCK3    Subr. W3UQCKMD Regular grid UQ scheme + obstructions.
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!     - The averaging is not performed around semi-transparent grid
!       points to avoid that leaking through barriers occurs.
!
!  6. Switches :
!
!       !/C90   Cray FORTRAN 90 compiler directives.
!       !/NEC   NEC SXF90 compiler directives.
!
!       !/MGP   Moving grid corrections.
!       !/MGG   Moving grid corrections.
!
!       !/S     Enable subroutine tracing.
!       !/Tn    Enable test output.
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
!/
!/ Public variables
!/
      PUBLIC
!/
!/ Private data
!/
      REAL, PRIVATE, PARAMETER:: TRNMIN = 0.95
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3MAP3
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Apr-2008 |
!/                  +-----------------------------------+
!/
!/    27-Feb-2000 : Origination.                        ( version 2.08 )
!/    10-Dec-2001 : Sub-grid obstructions.              ( version 2.14 )
!/                  (array allocation only.)
!/    17-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    09-Nov-2005 : Removing soft boundary option.      ( version 3.08 )
!/    01-Apr-2008 : Bug fix sec. 4 MAPSTA range check.  ( version 3.13 )
!/
!  1. Purpose :
!
!     Generate 'map' arrays for the ULTIMATE QUICKEST scheme.
!
!  2. Method :
!
!     MAPX2, MAPY2, MAPTH2 and MAPWN2 contain consecutive 1-D spatial
!     grid counters (e.g., IXY = (IX-1)*MY + IY). The arrays are
!     devided in three parts. For MAPX2, these ranges are :
!
!         1    - NMX0  Counters IXY for which grid point (IX,IY) and
!                      (IX+1,IY) both are active grid points.
!       NMX0+1 - NMX1  Id. only (IX,IY) active.
!       NMX1+1 - NMX2  Id. only (IX+1,IY) active.
!
!     The array MAPY2 has a similar layout varying IY instead of IX.
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
!      W3WAVE    Subr. W3WAVEMD Wave model routine.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!  7. Remarks :
!
!  8. Structure :
!
!     ------------------------------------------------------
!      1.   Map MAPX2
!        a  Range 1 to NMX0
!        b  Range NMX0+1 to NMX1
!        c  Range NMX1+1 to NMX2
!      2.   Map MAPY2
!        a  Range 1 to NMY0
!        b  Range NMY0+1 to NMY1
!        c  Range NMY1+1 to NMY2
!      3.   Map MAPAXY
!      4.   Map MAPCXY
!      5.   Maps for intra-spectral propagation
!        a  MAPTH2, MAPATK
!        b  MAPWN2
!     ------------------------------------------------------
!
!  9. Switches :
!
!     !/S   Enable subroutine tracing.
!     !/T   Enable test output.
!
! 10. Source code :
!/ ------------------------------------------------------------------- /
      USE W3GDATMD, ONLY: NK, NTH, NSPEC, NX, NY, NSEA, MAPSTA, MAPSF
      USE W3ADATMD, ONLY: NMX0, NMX1, NMX2, NMY0, NMY1, NMY2, NACT,   &
                          NCENT, MAPX2, MAPY2, MAPAXY, MAPCXY,        &
                          MAPTH2, MAPWN2
      USE W3ODATMD, ONLY: NDST
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IX, IY, IXY0, IX2, IY2, IX0, IY0,    &
                                 ISEA, IK, ITH, ISP, ISP0, ISP2, NCENTC
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Map MAPX2 ------------------------------------------------------ *
! 1.a Range 1 to NMX0
!
      NMX0   = 0
      DO IX=1, NX
        IXY0   = (IX-1)*NY
        IX2    = 1 + MOD(IX,NX)
        DO IY=2, NY-1
          IF ( MAPSTA(IY,IX).EQ.1 .AND. MAPSTA(IY,IX2).EQ.1 ) THEN
              NMX0   = NMX0 + 1
              MAPX2(NMX0) = IXY0 + IY
            END IF
          END DO
        END DO
!
! 1.b Range NMX0+1 to NMX1
!
      NMX1   = NMX0
      DO IX=1, NX
        IXY0   = (IX-1)*NY
        IX2    = 1 + MOD(IX,NX)
        DO IY=2, NY-1
          IF ( MAPSTA(IY,IX).EQ.1 .AND. MAPSTA(IY,IX2).NE.1 ) THEN
              NMX1   = NMX1 + 1
              MAPX2(NMX1) = IXY0 + IY
            END IF
          END DO
        END DO
!
! 1.c Range NMX1+1 to NMX2
!
      NMX2   = NMX1
      DO IX=1, NX
        IXY0   = (IX-1)*NY
        IX2    = 1 + MOD(IX,NX)
        DO IY=2, NY-1
          IF ( MAPSTA(IY,IX).NE.1 .AND. MAPSTA(IY,IX2).EQ.1 ) THEN
              NMX2   = NMX2 + 1
              MAPX2(NMX2) = IXY0 + IY
            END IF
          END DO
        END DO
!
! 2.  Map MAPY2 ------------------------------------------------------ *
! 2.a Range 1 to NMY0
!
      NMY0   = 0
      DO IX=1, NX
        IXY0   = (IX-1)*NY
        DO IY=1, NY-1
          IY2    = IY + 1
          IF ( MAPSTA(IY,IX).EQ.1 .AND. MAPSTA(IY2,IX).EQ.1 ) THEN
              NMY0   = NMY0 + 1
              MAPY2(NMY0) = IXY0 + IY
            END IF
          END DO
        END DO
!
! 2.b Range NMY0+1 to NMY1
!
      NMY1   = NMY0
      DO IX=1, NX
        IXY0   = (IX-1)*NY
        DO IY=1, NY-1
          IY2    = IY + 1
          IF ( MAPSTA(IY,IX).EQ.1 .AND. MAPSTA(IY2,IX).NE.1 ) THEN
              NMY1   = NMY1 + 1
              MAPY2(NMY1) = IXY0 + IY
            END IF
          END DO
        END DO
!
! 2.c Range NMY1+1 to NMY2
!
      NMY2   = NMY1
      DO IX=1, NX
        IXY0   = (IX-1)*NY
        DO IY=1, NY-1
          IY2    = IY + 1
          IF ( MAPSTA(IY,IX).NE.1 .AND. MAPSTA(IY2,IX).EQ.1 ) THEN
              NMY2   = NMY2 + 1
              MAPY2(NMY2) = IXY0 + IY
            END IF
           END DO
         END DO
!
! 3.  Map MAPAXY ----------------------------------------------------- *
!
      NACT   = 0
      DO IX=1, NX
        IY0    = (IX-1)*NY
        DO IY=2, NY-1
          IF ( MAPSTA(IY,IX).EQ.1 ) THEN
              NACT         = NACT + 1
              MAPAXY(NACT) = IY0 + IY
            END IF
          END DO
        END DO
!
! 4.  Map MAPCXY ----------------------------------------------------- *
!
      NCENT  = 0
      NCENTC = NSEA
      MAPCXY = 0
!
      DO ISEA=1,  NSEA
        IX      = MAPSF(ISEA,1)
        IX0    = IX-1
        IX2    = IX+1
        IY      = MAPSF(ISEA,2)
        IY0    = IY-1
        IY2    = IY+1
        IF ( IX .EQ. NX ) IX2 = 1
        IF ( IX .EQ. 1 ) IX0 = NX
        IF ( MAPSTA(IY,IX).EQ.2 .OR. MAPSTA(IY,IX).LT.0 ) THEN
            MAPCXY(NCENTC) = ISEA
            NCENTC = NCENTC - 1
          ELSE IF ( MAPSTA(IY0,IX0).GE.1 .AND.                     &
                    MAPSTA(IY0,IX ).GE.1 .AND.                     &
                    MAPSTA(IY0,IX2).GE.1 .AND.                     &
                    MAPSTA(IY ,IX0).GE.1 .AND.                     &
                    MAPSTA(IY ,IX2).GE.1 .AND.                     &
                    MAPSTA(IY2,IX0).GE.1 .AND.                     &
                    MAPSTA(IY2,IX ).GE.1 .AND.                     &
                    MAPSTA(IY2,IX2).GE.1 ) THEN
            NCENT  = NCENT + 1
            MAPCXY(NCENT) = ISEA
          ELSE
            MAPCXY(NCENTC) = ISEA
            NCENTC = NCENTC - 1
          END IF
        END DO
!
! 5.  Maps for intra-spectral propagation ---------------------------- *
!
      IF ( MAPTH2(1) .NE. 0 ) RETURN
!
! 5.a MAPTH2 and MAPBTK
!
      DO IK=1, NK
        DO ITH=1, NTH
          ISP    = ITH + (IK-1)*NTH
          ISP2   = (IK+1) + (ITH-1)*(NK+2)
          MAPTH2(ISP) = ISP2
          END DO
        END DO
!
! 5.b MAPWN2
!
      ISP0   = 0
      DO IK=1, NK-1
        DO ITH=1, NTH
          ISP0   = ISP0 + 1
          ISP2   = (IK+1) + (ITH-1)*(NK+2)
          MAPWN2(ISP0) = ISP2
          END DO
        END DO
!
      DO ITH=1, NTH
        ISP0   = ISP0 + 1
        ISP2   = NK+1 + (ITH-1)*(NK+2)
        MAPWN2(ISP0) = ISP2
        END DO
!
      DO ITH=1, NTH
        ISP0   = ISP0 + 1
        ISP2   = 1 + (ITH-1)*(NK+2)
        MAPWN2(ISP0) = ISP2
        END DO
!
      RETURN
!
! Formats
!
!/
!/ End of W3MAP3 ----------------------------------------------------- /
!/
      END SUBROUTINE W3MAP3
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3MAPT
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         17-Dec-2004 |
!/                  +-----------------------------------+
!/
!/    10-Dec-2001 : Origination.                        ( version 2.14 )
!/    17-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/
!  1. Purpose :
!
!     Generate 'map' arrays for the ULTIMATE QUICKEST scheme to combine
!     GSE alleviation with obstructions.
!
!  2. Method :
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
!      W3WAVE    Subr. W3WAVEMD Wave model routine.
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
!     !/S   Enable subroutine tracing.
!
! 10. Source code :
!/ ------------------------------------------------------------------- /
      USE W3GDATMD, ONLY: NX, NY, NSEA, MAPSF
      USE W3ADATMD, ONLY: ATRNX, ATRNY, MAPTRN
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: ISEA, IXY
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Map MAPTRN ----------------------------------------------------- *
!
      DO ISEA=1, NSEA
        IXY         = MAPSF(ISEA,3)
        MAPTRN(IXY) = MIN( ATRNX(IXY,1) ,ATRNY(IXY,-1) ,              &
                           ATRNY(IXY,1), ATRNY(IXY,-1) ) .LT. TRNMIN
        END DO
!
      RETURN
!
! Formats
!/
!/ End of W3MAPT ----------------------------------------------------- /
!/
      END SUBROUTINE W3MAPT
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3XYP3 ( ISP, FACX, FACY, DTG, MAPSTA, MAPFS, VQ,    &
                          VGX, VGY )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         05-Mar-2008 |
!/                  +-----------------------------------+
!/
!/    27-Feb-2000 : Origination.                        ( version 2.08 )
!/    17-Sep-2000 : Clean-up.                           ( version 2.13 )
!/    10-Dec-2001 : Sub-grid obstructions.              ( version 2.14 )
!/    16-Oct-2002 : Change INTENT for ATRNRX/Y.         ( version 3.00 )
!/    26-Dec-2002 : Moving grid version.                ( version 3.02 )
!/    01-Aug-2003 : Moving grid GSE correction.         ( version 3.03 )
!/    17-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    07-Sep-2005 : Upgrade XY boundary conditions.     ( version 3.08 )
!/    09-Nov-2005 : Removing soft boundary option.      ( version 3.08 )
!/    05-Mar-2008 : Added NEC sxf90 compiler directives.
!/                  (Chris Bunney, UK Met Office)       ( version 3.13 )
!/
!  1. Purpose :
!
!     Propagation in phyiscal space for a given spectral component.
!
!  2. Method :
!
!     Third-order ULTIMATE QUICKEST scheme with averaging.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       ISP     Int.   I   Number of spectral bin (IK-1)*NTH+ITH
!       FACX/Y  Real   I   Factor in propagation velocity.
!                          ( 1 or 0 * DT / DX )
!       DTG     Real   I   Total time step.
!       MAPSTA  I.A.   I   Grid point status map.
!       MAPFS   I.A.   I   Storage map.
!       VQ      R.A.  I/O  Field to propagate.
!       VGX/Y   Real   I   Speed of grid.
!     ----------------------------------------------------------------
!
!     Local variables.
!     ----------------------------------------------------------------
!       NTLOC   Int   Number of local time steps.
!       DTLOC   Real  Local propagation time step.
!       VCFL0X  R.A.  Local courant numbers for absolute group vel.
!                     using local X-grid step.
!       VCFL0Y  R.A.  Id. in Y.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       W3QCK3   Actual propagation algorithm
!
!       STRACE   Service routine.
!
!  5. Called by :
!
!       W3WAVE   Wave model routine.
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!     - Note that the ULTIMATE limiter does not guarantee non-zero
!       energies.
!     - The present scheme shows a strong distorsion when propaga-
!       ting a field under an angle with the grid in a truly 2-D
!       fashion. Propagation is therefore split along the two
!       axes.
!     - Two boundary treatments are available. The first uses real
!       boundaries in each space. In this case waves will not
!       penetrate in narrow straights under an angle with the grid.
!       This behavior is improved by using a 'soft' option, in
!       which the 'X' or 'Y' sweep allows for energy to go onto
!       the land. This improves the above behavior, but implies
!       that X-Y connenctions are required in barriers for them
!       to become inpenetrable.
!     - Note that unlike in W3XYP2, isotropic diffusion is never
!       used for growth.
!
!  8. Structure :
!
!     ---------------------------------------------
!       1.  Preparations
!         a Set constants
!         b Initialize arrays
!       2.  Prepare arrays
!         a Velocities and 'Q'
!       3.  Loop over sub-steps
!       ----------------------------------------
!         a Average
!         b Propagate
!         c Update boundaries.
!       ----------------------------------------
!       4.  Store Q field in spectra
!     ---------------------------------------------
!
!  9. Switches :
!
!       !/NEC   Enable NEC SXF90 compiler directives.
!
!       !/S     Enable subroutine tracing.
!
!       !/MGP   Moving grid corrections.
!       !/MGG   Moving grid corrections.
!
!       !/T     Enable general test output.
!       !/T0    Dump of precalcaulted interpolation data.
!       !/T1    Dump of input field and fluxes.
!       !/T2    Dump of output field.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3CONSTANTS
!
      USE W3TIMEMD, ONLY: DSEC21
!
      USE W3GDATMD, ONLY: NX, NY, NSEA, MAPSF, DTCFL, CLAT, CLATS,    &
                          GLOBAL, FLCX, FLCY, NK, NTH, DTH, XFR,      &
                          ECOS, ESIN, SIG, WDCG, WDTH, PFMOVE
      USE W3WDATMD, ONLY: TIME
      USE W3ADATMD, ONLY: NMX0, NMX1, NMX2, NMY0, NMY1, NMY2, NACT,   &
                          NCENT, MAPX2, MAPY2, MAPAXY, MAPCXY,        &
                          MAPTRN, CG, CX, CY, ATRNX, ATRNY, ITIME
      USE W3IDATMD, ONLY: FLCUR
      USE W3ODATMD, ONLY: NDSE, NDST, FLBPI, NBI, TBPI0, TBPIN,       &
                          ISBPI, BBPI0, BBPIN
      USE W3UQCKMD
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: ISP, MAPSTA(NY*NX), MAPFS(NY*NX)
      REAL, INTENT(IN)        :: FACX, FACY, DTG, VGX, VGY
      REAL, INTENT(INOUT)     :: VQ(1-NY:NY*(NX+2))
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: ITH, IK, NTLOC, ITLOC, ISEA, IXY, IP
      INTEGER                 :: IX, IY, IXC, IYC, IIXY1(NY),         &
                                 IIXY2(NY), IIXY3(NY), IIXY4(NY), IBI
      REAL                    :: CG0, CGA, CGN, CGX, CGY, CXC, CYC,   &
                                 CXMIN, CXMAX, CYMIN, CYMAX
      REAL                    :: CGC, FGSE = 1., FTH, FCG
      REAL                    :: DTLOC, CCOS, CSIN, CCURX, CCURY
      REAL                    :: DXCGN, DYCGN, DXCGS, DYCGS, DXCGC,   &
                                 DYCGC, RDI1(NY), RDI2(NY), RDI3(NY), &
                                 RDI4(NY), RD1, RD2, RD3, RD4
      LOGICAL                 :: YFIRST
!/
!/ Automatic work arrays
!/
      INTEGER                 :: MAPSTX(1-2*NY:NY*(NX+2))
      REAL                    :: VLCFLX((NX+1)*NY), VLCFLY((NX+1)*NY),&
                                 AQ(1-NY:NY*(NX+2))
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Preparations --------------------------------------------------- *
! 1.a Set constants
!
      ITH    = 1 + MOD(ISP-1,NTH)
      IK     = 1 + (ISP-1)/NTH
!
      CG0    = 0.575 * GRAV / SIG(1)
      CGA    = 0.575 * GRAV / SIG(IK)
      CGX    = CGA * ECOS(ITH)
      CGY    = CGA * ESIN(ITH)
      CGC    = SQRT ( CGX**2 + CGY**2 )
!
      IF ( FLCUR ) THEN
          CXMIN  = MINVAL ( CX(1:NSEA) )
          CXMAX  = MAXVAL ( CX(1:NSEA) )
          CYMIN  = MINVAL ( CY(1:NSEA) )
          CYMAX  = MAXVAL ( CY(1:NSEA) )
          IF ( ABS(CGX+CXMIN) .GT. ABS(CGX+CXMAX) ) THEN
              CGX    = CGX + CXMIN
            ELSE
              CGX    = CGX + CXMAX
            END IF
          IF ( ABS(CGY+CYMIN) .GT. ABS(CGY+CYMAX) ) THEN
              CGY    = CGY + CYMIN
            ELSE
              CGY    = CGY + CYMAX
            END IF
          CXC    = MAX ( ABS(CXMIN) , ABS(CXMAX) )
          CYC    = MAX ( ABS(CYMIN) , ABS(CYMAX) )
        ELSE
          CXC    = 0.
          CYC    = 0.
        END IF
!
      CGN    = MAX ( ABS(CGX) , ABS(CGY) , CXC, CYC, 0.001*CG0 )
!
      NTLOC  = 1 + INT(DTG/(DTCFL*CG0/CGN))
      DTLOC  = DTG / REAL(NTLOC)
!
      CCOS   = FACX * ECOS(ITH) / REAL(NTLOC)
      CSIN   = FACY * ESIN(ITH) / REAL(NTLOC)
!
      CCURX  = FACX / REAL(NTLOC)
      CCURY  = FACY / REAL(NTLOC)
!
      YFIRST = MOD(ITIME,2) .EQ. 0
!
! 1.b Initialize arrays
!
      VLCFLX = 0.
      VLCFLY = 0.
!
      MAPSTX(1:NX*NY) = MAPSTA(1:NX*NY)
!
      IF ( GLOBAL ) THEN
          MAPSTX(1-2*NY:0)            = MAPSTA((NX-2)*NY+1:NX*NY)
          MAPSTX(NX*NY+1:NX*NY+2*NY)  = MAPSTA(1:2*NY)
        ELSE
          MAPSTX(1-2*NY:0)            = 0
          MAPSTX(NX*NY+1:NX*NY+2*NY)  = 0
        END IF
!
! 1.c Pre-calculate interpolation info
!
      FTH    =   FGSE * WDTH * DTH / REAL(NTLOC)
      FCG    =   FGSE * WDCG * 0.5 * (XFR-1./XFR) / REAL(NTLOC)
!
      DO IY=1, NY
!
! 1.c.1 Normal and parallel width ...
!
        DXCGN  = - FTH * FACX * ESIN(ITH) / CLAT(IY)
        DYCGN  =   FTH * FACY * ECOS(ITH)
        DXCGS  =   FCG * CCOS / CLAT(IY)
        DYCGS  =   FCG * CSIN
!
! 1.c.2 "Sum" corner (and mirror image) ...
!
        DXCGC  = DXCGN + DXCGS
        DYCGC  = DYCGN + DYCGS
!
        IXC    = NY
        IF ( DXCGC .LT. 0. ) IXC = - IXC
        IYC    = 1
        IF ( DYCGC .LT. 0. ) IYC = - IYC
!
        IIXY1(IY) = IXC + IYC
        IF ( ABS(DXCGC) .GT. ABS(DYCGC) ) THEN
            IIXY2(IY) = IXC
            RDI1 (IY) = ABS(DYCGC/DXCGC)
            RDI2 (IY) = ABS(DXCGC)
          ELSE
            IIXY2(IY) = IYC
            IF ( ABS(DYCGC) .GT. 1.E-5 ) THEN
                RDI1(IY) = ABS(DXCGC/DYCGC)
              ELSE
                RDI1(IY) = 1.
              END IF
            RDI2(IY) = ABS(DYCGC)
          END IF
!
! 1.c.2 "Difference" corner (and mirror image) ...
!
        DXCGC  = DXCGN - DXCGS
        DYCGC  = DYCGN - DYCGS
!
        IXC    = NY
        IF ( DXCGC .LT. 0. ) IXC = - IXC
        IYC    = 1
        IF ( DYCGC .LT. 0. ) IYC = - IYC
!
        IIXY3(IY) = IXC + IYC
        IF ( ABS(DXCGC) .GT. ABS(DYCGC) ) THEN
            IIXY4(IY) = IXC
            RDI3 (IY) = ABS(DYCGC/DXCGC)
            RDI4 (IY) = ABS(DXCGC)
          ELSE
            IIXY4(IY) = IYC
            IF ( ABS(DYCGC) .GT. 1.E-5 ) THEN
                RDI3(IY) = ABS(DXCGC/DYCGC)
              ELSE
                RDI3(IY) = 1.
              END IF
            RDI4(IY) = ABS(DYCGC)
          END IF
!
        END DO
!
! 2.  Calculate velocities and diffusion coefficients ---------------- *
! 2.a Velocities
!
!     Q     = ( A / CG * CLATS )
!     LCFLX = ( COS*CG / CLATS ) * DT / DX
!     LCFLY = (     SIN*CG )     * DT / DY
!
      DO ISEA=1, NSEA
        IXY         = MAPSF(ISEA,3)
        VQ    (IXY) = VQ(IXY) / CG(IK,ISEA) * CLATS(ISEA)
        VLCFLX(IXY) = CCOS * CG(IK,ISEA) / CLATS(ISEA)
        VLCFLY(IXY) = CSIN * CG(IK,ISEA)
        END DO
!
      IF ( FLCUR ) THEN
          DO ISEA=1, NSEA
            IXY         = MAPSF(ISEA,3)
            VLCFLX(IXY) = VLCFLX(IXY) + CCURX*CX(ISEA)/CLATS(ISEA)
            VLCFLY(IXY) = VLCFLY(IXY) + CCURY*CY(ISEA)
           END DO
        END IF
!
! 3.  Loop over sub-steps -------------------------------------------- *
!
      DO ITLOC=1, NTLOC
!
! 3.a Average
!
        AQ     = VQ
        VQ     = 0.
!
! 3.a.1 Central points
!
        DO IP=1, NCENT
          ISEA    = MAPCXY(IP)
          IXY     = MAPSF(ISEA,3)
          IY      = MAPSF(ISEA,2)
          IF ( MAPTRN(IXY) ) THEN
              VQ(IXY) = AQ(IXY)
            ELSE
              RD1     = RDI1(IY)
              RD2     = MIN ( 1. , RDI2(IY) * CG(IK,ISEA) )
              RD3     = RDI3(IY)
              RD4     = MIN ( 1. , RDI4(IY) * CG(IK,ISEA) )
              VQ(IXY          ) = VQ(IXY          )                   &
                                   + AQ(IXY) * (3.-RD2-RD4)/3.
              VQ(IXY+IIXY1(IY)) = VQ(IXY+IIXY1(IY))                   &
                                   + AQ(IXY) * RD2*RD1/6.
              VQ(IXY+IIXY2(IY)) = VQ(IXY+IIXY2(IY))                   &
                                   + AQ(IXY) * (1.-RD1)*RD2/6.
              VQ(IXY+IIXY3(IY)) = VQ(IXY+IIXY3(IY))                   &
                                   + AQ(IXY) * RD4*RD3/6.
              VQ(IXY+IIXY4(IY)) = VQ(IXY+IIXY4(IY))                   &
                                   + AQ(IXY) * (1.-RD3)*RD4/6.
              VQ(IXY-IIXY1(IY)) = VQ(IXY-IIXY1(IY))                   &
                                   + AQ(IXY) * RD2*RD1/6.
              VQ(IXY-IIXY2(IY)) = VQ(IXY-IIXY2(IY))                   &
                                   + AQ(IXY) * (1.-RD1)*RD2/6.
              VQ(IXY-IIXY3(IY)) = VQ(IXY-IIXY3(IY))                   &
                                   + AQ(IXY) * RD4*RD3/6.
              VQ(IXY-IIXY4(IY)) = VQ(IXY-IIXY4(IY))                   &
                                   + AQ(IXY) * (1.-RD3)*RD4/6.
            END IF
          END DO
!
! 3.a.2 Near-coast points
!
        DO IP=NCENT+1, NSEA
          ISEA    = MAPCXY(IP)
          IX      = MAPSF(ISEA,1)
          IY      = MAPSF(ISEA,2)
          IXY     = MAPSF(ISEA,3)
          IF ( MAPSTA(IXY) .LE. 0 ) CYCLE
          IF ( MAPTRN(IXY) ) THEN
              VQ(IXY) = AQ(IXY)
            ELSE
              RD1     = RDI1(IY)
              RD3     = RDI3(IY)
              RD2     = MIN ( 1. , RDI2(IY) * CG(IK,ISEA) )
              RD4     = MIN ( 1. , RDI4(IY) * CG(IK,ISEA) )
              VQ(IXY          ) = VQ(IXY          )                   &
                                   + AQ(IXY) * (3.-RD2-RD4)/3.
!
              IXC    = SIGN(NY,IIXY1(IY))
              IYC    = IIXY1(IY) - IXC
              IF ( MAPSTX(IXY+IIXY1(IY)) .GE. 1 .AND.                 &
                   .NOT. ( MAPSTX(IXY+IXC).LE.0 .AND.                 &
                           MAPSTX(IXY+IYC).LE.0 ) ) THEN
                       VQ(IXY+IIXY1(IY)) = VQ(IXY+IIXY1(IY))          &
                                            + AQ(IXY) * RD2*RD1/6.
                ELSE
                       VQ(IXY          ) = VQ(IXY          )          &
                                            + AQ(IXY) * RD2*RD1/6.
                END IF
              IF ( MAPSTX(IXY-IIXY1(IY)) .GE. 1 .AND.                 &
                   .NOT. ( MAPSTX(IXY-IXC).LE.0 .AND.                 &
                           MAPSTX(IXY-IYC).LE.0 ) ) THEN
                       VQ(IXY-IIXY1(IY)) = VQ(IXY-IIXY1(IY))          &
                                            + AQ(IXY) * RD2*RD1/6.
                ELSE
                       VQ(IXY          ) = VQ(IXY          )          &
                                            + AQ(IXY) * RD2*RD1/6.
                END IF
 
              IF ( MAPSTX(IXY+IIXY2(IY)) .GE. 1 ) THEN
                       VQ(IXY+IIXY2(IY)) = VQ(IXY+IIXY2(IY))          &
                                            + AQ(IXY) * (1.-RD1)*RD2/6.
                ELSE
                       VQ(IXY          ) = VQ(IXY          )          &
                                            + AQ(IXY) * (1.-RD1)*RD2/6.
                END IF
              IF ( MAPSTX(IXY-IIXY2(IY)) .GE. 1 ) THEN
                       VQ(IXY-IIXY2(IY)) = VQ(IXY-IIXY2(IY))          &
                                            + AQ(IXY) * (1.-RD1)*RD2/6.
                ELSE
                       VQ(IXY          ) = VQ(IXY          )          &
                                            + AQ(IXY) * (1.-RD1)*RD2/6.
                END IF
!
              IXC    = SIGN(NY,IIXY3(IY))
              IYC    = IIXY3(IY) - IXC
              IF ( MAPSTX(IXY+IIXY3(IY)) .GE. 1 .AND.                 &
                   .NOT. ( MAPSTX(IXY+IXC).LE.0 .AND.                 &
                           MAPSTX(IXY+IYC).LE.0 ) ) THEN
                       VQ(IXY+IIXY3(IY)) = VQ(IXY+IIXY3(IY))          &
                                            + AQ(IXY) * RD4*RD3/6.
                ELSE
                       VQ(IXY          ) = VQ(IXY          )          &
                                            + AQ(IXY) * RD4*RD3/6.
                END IF
              IF ( MAPSTX(IXY-IIXY3(IY)) .GE. 1 .AND.                 &
                   .NOT. ( MAPSTX(IXY-IXC).LE.0 .AND.                 &
                           MAPSTX(IXY-IYC).LE.0 ) ) THEN
                       VQ(IXY-IIXY3(IY)) = VQ(IXY-IIXY3(IY))          &
                                            + AQ(IXY) * RD4*RD3/6.
                ELSE
                       VQ(IXY          ) = VQ(IXY          )          &
                                            + AQ(IXY) * RD4*RD3/6.
                END IF
!
              IF ( MAPSTX(IXY+IIXY4(IY)) .GE. 1 ) THEN
                       VQ(IXY+IIXY4(IY)) = VQ(IXY+IIXY4(IY))          &
                                            + AQ(IXY) * (1.-RD3)*RD4/6.
                ELSE
                       VQ(IXY          ) = VQ(IXY          )          &
                                            + AQ(IXY) * (1.-RD3)*RD4/6.
                END IF
              IF ( MAPSTX(IXY-IIXY4(IY)) .GE. 1 ) THEN
                       VQ(IXY-IIXY4(IY)) = VQ(IXY-IIXY4(IY))          &
                                            + AQ(IXY) * (1.-RD3)*RD4/6.
                ELSE
                       VQ(IXY          ) = VQ(IXY          )          &
                                            + AQ(IXY) * (1.-RD3)*RD4/6.
                END IF
!
            END IF
          END DO
!
! 3.a.3 Restore boundary data
!
        DO IXY=1, NX*NY
          IF ( MAPSTA(IXY).EQ.2 ) VQ(IXY) = AQ(IXY)
          END DO
!
! 3.a.4 Global closure
!
        IF ( GLOBAL ) THEN
            DO IY=1, NY
              VQ(IY          ) = VQ(IY          ) + VQ(NX*NY+IY)
              VQ((NX-1)*NY+IY) = VQ((NX-1)*NY+IY) + VQ(IY-NY)
              END DO
          END IF
!
! 3.b Propagate fields
!
        IF ( YFIRST ) THEN
            IF ( FLCY ) CALL W3QCK3                                   &
                        (NX, NY, NX, NY, VLCFLY, ATRNY, VQ,           &
                         .FALSE., 1, MAPAXY, NACT, MAPY2, NMY0,       &
                         NMY1, NMY2, NDSE, NDST )
            IF ( FLCX ) CALL W3QCK3                                   &
                        (NX, NY, NX, NY, VLCFLX, ATRNX, VQ,           &
                         GLOBAL, NY, MAPAXY, NACT, MAPX2, NMX0,       &
                         NMX1, NMX2, NDSE, NDST )
          ELSE
            IF ( FLCX ) CALL W3QCK3                                   &
                        (NX, NY, NX, NY, VLCFLX, ATRNX, VQ,           &
                         GLOBAL, NY, MAPAXY, NACT, MAPX2, NMX0,       &
                         NMX1, NMX2, NDSE, NDST )
            IF ( FLCY ) CALL W3QCK3                                   &
                        (NX, NY, NX, NY, VLCFLY, ATRNY, VQ,           &
                         .FALSE., 1, MAPAXY, NACT, MAPY2, NMY0,       &
                         NMY1, NMY2, NDSE, NDST )
          END IF
!
! 3.c Update boundaries
!
        IF ( FLBPI ) THEN
            RD1    = DSEC21 ( TBPI0, TIME ) - DTG *                   &
                                      REAL(NTLOC-ITLOC)/REAL(NTLOC)
            RD2    = DSEC21 ( TBPI0, TBPIN )
            IF ( RD2 .GT. 0.001 ) THEN
                 RD2    = MIN(1.,MAX(0.,RD1/RD2))
                 RD1    = 1. - RD2
              ELSE
                 RD1    = 0.
                 RD2    = 1.
              END IF
            DO IBI=1, NBI
              IXY     = MAPSF(ISBPI(IBI),3)
              VQ(IXY) = ( RD1*BBPI0(ISP,IBI) + RD2*BBPIN(ISP,IBI) )   &
                          / CG(IK,ISBPI(IBI)) * CLATS(ISBPI(IBI))
              END DO
          END IF
!
        YFIRST = .NOT. YFIRST
        END DO
!
! 4.  Store results in VQ in proper format --------------------------- *
!
      DO ISEA=1, NSEA
        IXY    = MAPSF(ISEA,3)
        IF ( MAPSTA(IXY) .GT. 0 ) THEN
            VQ(IXY) =  MAX ( 0. , CG(IK,ISEA)/CLATS(ISEA)*VQ(IXY) )
          END IF
        END DO
!
      RETURN
!
! Formats
!
!/
!/ End of W3XYP3 ----------------------------------------------------- /
!/
      END SUBROUTINE W3XYP3
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3KTP3 ( ISEA, FACTH, FACK, CTHG0, CG, WN, DEPTH,    &
                          DDDX, DDDY, CX, CY, DCXDX, DCXDY,           &
                          DCYDX, DCYDY, VA )
!/
!/    *** THIS ROUTINE SHOULD BE IDENTICAL TO W3KTP2 ***
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         17-Dec-2004 |
!/                  +-----------------------------------+
!/
!/    14-Feb-2000 : Origination.                        ( version 2.08 )
!/    17-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/
!  1. Purpose :
!
!     Propagation in spectral space.
!
!  2. Method :
!
!     Third order QUICKEST scheme with ULTIMATE limiter.
!
!     As with the spatial propagation, the two spaces are considered
!     independently, but the propagation is performed in a 2-D space.
!     Compared to the propagation in physical space, the directions
!     rerpesent a closed space and are therefore comparable to the
!     longitudinal or 'X' propagation. The wavenumber space has to be
!     extended to allow for boundary treatment. Using a simple first
!     order boundary treatment at both sided, two points need to
!     be added. This implies that the spectrum needs to be extended,
!     shifted and rotated, as is performed using MAPTH2 as set
!     in W3MAP3.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       ISEA    Int.   I   Number of sea point.
!       FACTH/K Real   I   Factor in propagation velocity.
!       CTHG0   Real   I   Factor in great circle refracftion term.
!       MAPxx2  I.A.   I   Propagation and storage maps.
!       CG      R.A.   I   Local group velocities.
!       WN      R.A.   I   Local wavenumbers.
!       DEPTH   R.A.   I   Depth.
!       DDDx    Real   I   Depth gradients.
!       CX/Y    Real   I   Current components.
!       DCxDx   Real   I   Current gradients.
!       VA      R.A.  I/O  Spectrum.
!     ----------------------------------------------------------------
!
!     Local variables.
!     ----------------------------------------------------------------
!       DSDD    R.A.  Partial derivative of sigma for depth.
!       FDD, FDU, FDG, FCD, FCU
!               R.A.  Directionally varying part of depth, current and
!                     great-circle refraction terms and of consit.
!                     of Ck term.
!       CFLT-K  R.A.  Propagation velocities of local fluxes.
!       DB      R.A.  Wavenumber band widths at cell centers.
!       DM      R.A.  Wavenumber band widths between cell centers and
!                     next cell center.
!       Q       R.A.  Extracted spectrum
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       W3QCK1   Actual propagation routine.
!       W3QCK2   Actual propagation routine.
!       STRACE   Service routine.
!
!  5. Called by :
!
!       W3WAVE   Wave model routine.
!
!  6. Error messages :
!
!       None.
!
!  8. Structure :
!
!     -----------------------------------------------------------------
!       1.  Preparations
!         a Initialize arrays
!         b Set constants and counters
!       2.  Point  preparations
!         a Calculate DSDD
!         b Extract spectrum
!       3.  Refraction velocities
!         a Filter level depth reffraction.
!         b Depth refratcion velocity.
!         c Current refraction velocity.
!       4.  Wavenumber shift velocities
!         a Prepare directional arrays
!         b Calcuate velocity.
!       5.  Propagate.
!       6.  Store results.
!     -----------------------------------------------------------------
!
!  9. Switches :
!
!       !/S     Enable subroutine tracing.
!       !/T     Enable general test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3CONSTANTS
      USE W3GDATMD, ONLY: NK, NK2, NTH, NSPEC, SIG, DSIP, ECOS, ESIN, &
                          EC2, ESC, ES2, FACHFA, MAPWN, FLCTH, FLCK,  &
                          CTMAX
      USE W3ADATMD, ONLY: MAPTH2, MAPWN2, ITIME
      USE W3IDATMD, ONLY: FLCUR
      USE W3ODATMD, ONLY: NDSE, NDST
      USE W3UQCKMD
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: ISEA
      REAL, INTENT(IN)        :: FACTH, FACK, CTHG0, CG(0:NK+1),      &
                                 WN(0:NK+1), DEPTH, DDDX, DDDY,       &
                                 CX, CY, DCXDX, DCXDY, DCYDX, DCYDY
      REAL, INTENT(INOUT)     :: VA(NSPEC)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: ITH, IK, ISP
      REAL                    :: FDDMAX, FDG, FKD, FKD0, DCYX,        &
                                 DCXXYY, DCXY, DCXX, DCXYYX, DCYY
      REAL                    :: DSDD(0:NK+1), FRK(NK), FRG(NK),      &
                                 FKC(NTH), VQ(-NK-1:NK2*(NTH+2)),     &
                                 DB(NK2,NTH+1), DM(NK2,0:NTH+1),      &
                                 VCFLT(NK2*(NTH+1)), CFLK(NK2,NTH)
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Preparations --------------------------------------------------- *
! 1.a Initialize arrays
!
      IF ( FLCK ) VQ = 0.
!
! 2.  Preparation for point ------------------------------------------ *
! 2.a Array with partial derivative of sigma versus depth
!
      DO IK=0, NK+1
        IF ( DEPTH*WN(IK) .LT. 5. ) THEN
            DSDD(IK) = MAX ( 0. ,                                     &
                CG(IK)*WN(IK)-0.5*SIG(IK) ) / DEPTH
          ELSE
            DSDD(IK) = 0.
          END IF
        END DO
!
! 2.b Extract spectrum
!
      DO ISP=1, NSPEC
        VQ(MAPTH2(ISP)) = VA(ISP)
        END DO
!
! 3.  Refraction velocities ------------------------------------------ *
!
      IF ( FLCTH ) THEN
!
! 3.a Set slope filter for depth refraction
!
          FDDMAX = 0.
          FDG    = FACTH * CTHG0
!
          DO ITH=1, NTH/2
            FDDMAX = MAX(FDDMAX,ABS(ESIN(ITH)*DDDX-ECOS(ITH)*DDDY))
            END DO
!
          DO IK=1, NK
            FRK(IK) = FACTH * DSDD(IK) / WN(IK)
            FRK(IK) = FRK(IK) / MAX ( 1. , FRK(IK)*FDDMAX/CTMAX )
            FRG(IK) = FDG * CG(IK)
            END DO
!
! 3.b Depth refraction and great-circle propagation
!
          DO ISP=1, NSPEC
            VCFLT(MAPTH2(ISP)) = FRG(MAPWN(ISP)) * ECOS(ISP)          &
             + FRK(MAPWN(ISP)) * ( ESIN(ISP)*DDDX - ECOS(ISP)*DDDY )
            END DO
!
! 3.c Current refraction
!
          IF ( FLCUR ) THEN
!
              DCYX   = FACTH *   DCYDX
              DCXXYY = FACTH * ( DCXDX - DCYDY )
              DCXY   = FACTH *   DCXDY
!
              DO ISP=1, NSPEC
                VCFLT(MAPTH2(ISP)) = VCFLT(MAPTH2(ISP)) +             &
                  ES2(ISP)*DCYX  + ESC(ISP)*DCXXYY - EC2(ISP)*DCXY
                  END DO
!
            END IF
!
        END IF
!
! 4.  Wavenumber shift velocities ------------------------------------ *
!     FACK is just the time step, which is accounted for in W3QCK2
!
      IF ( FLCK ) THEN
!
! 4.a Directionally dependent part
!
          DCXX   =  -   DCXDX
          DCXYYX =  - ( DCXDY + DCYDX )
          DCYY   =  -   DCYDY
          FKD    =    ( CX*DDDX + CY*DDDY )
!
          DO ITH=1, NTH
            FKC(ITH) = EC2(ITH)*DCXX +                                &
                       ESC(ITH)*DCXYYX + ES2(ITH)*DCYY
            END DO
!
! 4.b Velocities
!
          DO IK=0, NK+1
            FKD0   = FKD / CG(IK) * DSDD(IK)
            DO ITH=1, NTH
              CFLK(IK+1,ITH) = FKD0 + WN(IK)*FKC(ITH)
              END DO
            END DO
!
! 4.c Band widths
!
          DO IK=0, NK
            DB(IK+1,1) = DSIP(IK) / CG(IK)
            DM(IK+1,1) = WN(IK+1) - WN(IK)
            END DO
          DB(NK+2,1) = DSIP(NK+1) / CG(NK+1)
          DM(NK+2,1) = 0.
!
          DO ITH=2, NTH
            DO IK=1, NK+2
              DB(IK,ITH) = DB(IK,1)
              DM(IK,ITH) = DM(IK,1)
              END DO
            END DO
!
        END IF
!
! 5.  Propagate ------------------------------------------------------ *
!
      IF ( MOD(ITIME,2) .EQ. 0 ) THEN
          IF ( FLCK ) THEN
              DO ITH=1, NTH
                VQ(NK+2+(ITH-1)*NK2) = FACHFA * VQ(NK+1+(ITH-1)*NK2)
                END DO
              CALL W3QCK2 ( NTH, NK2, NTH, NK2, CFLK, FACK, DB, DM,   &
                            VQ, .FALSE., 1, MAPTH2, NSPEC,            &
                            MAPWN2, NSPEC-NTH, NSPEC, NSPEC+NTH,      &
                            NDSE, NDST )
            END IF
          IF ( FLCTH ) THEN
              CALL W3QCK1 ( NTH, NK2, NTH, NK2, VCFLT, VQ, .TRUE.,    &
                            NK2, MAPTH2, NSPEC, MAPTH2, NSPEC, NSPEC, &
                            NSPEC, NDSE, NDST )
            END IF
        ELSE
          IF ( FLCTH ) THEN
              CALL W3QCK1 ( NTH, NK2, NTH, NK2, VCFLT, VQ, .TRUE.,    &
                            NK2, MAPTH2, NSPEC, MAPTH2, NSPEC, NSPEC, &
                            NSPEC, NDSE, NDST )
            END IF
          IF ( FLCK )  THEN
              DO ITH=1, NTH
                VQ(NK+2+(ITH-1)*NK2) = FACHFA * VQ(NK+1+(ITH-1)*NK2)
                END DO
              CALL W3QCK2 ( NTH, NK2, NTH, NK2, CFLK, FACK, DB, DM,   &
                            VQ, .FALSE., 1, MAPTH2, NSPEC,            &
                            MAPWN2, NSPEC-NTH, NSPEC, NSPEC+NTH,      &
                            NDSE, NDST )
          END IF
        END IF
!
! 6.  Store reults --------------------------------------------------- *
!
      DO ISP=1, NSPEC
        VA(ISP) = VQ(MAPTH2(ISP))
        END DO
!
      RETURN
!
! Formats
!
!/
!/ End of W3KTP3 ----------------------------------------------------- /
!/
      END SUBROUTINE W3KTP3
!/
!/ End of module W3PRO3MD -------------------------------------------- /
!/
      END MODULE W3PRO3MD
