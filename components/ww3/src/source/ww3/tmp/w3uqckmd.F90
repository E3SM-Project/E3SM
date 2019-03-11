#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE W3UQCKMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         27-May-2014 |
!/                  +-----------------------------------+
!/
!/    08-Feb-2001 : Origination of module. Routines     ( version 2.08 )
!/                  taken out of w3pro2md.ftn
!/    13-Nov-2001 : Version with obstacles added.       ( version 2.14 )
!/    16-Oct-2002 : Fix par list W3QCK3.                ( version 3.00 )
!/    05-Mar-2008 : Added NEC sxf90 compiler directives.
!/                  (Chris Bunney, UK Met Office)       ( version 3.13 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/    30-Oct-2009 : Fixed a couple of doc lines.        ( version 3.14 )
!/                  (T. J. Campbell, NRL)
!/    27-May-2014 : Added OMPH switches in W3QCK3.      ( version 5.02 )
!/
!/    Copyright 2009-2014 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Portable ULTIMATE QUICKEST schemes.
!
!  2. Variables and types :
!
!     None.
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3QCK1    Subr. Public   Original ULTIMATE QUICKEST scheme.
!      W3QCK2    Subr. Public   UQ scheme for irregular grid.
!      W3QCK3    Subr. Public   Original ULTIMATE QUICKEST with obst.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!     - STRACE and !/S irrelevant for running code. The module is
!       therefore fully portable to any other model.
!
!  6. Switches :
!
!       !/C90   Cray FORTRAN 90 compiler directives.
!       !/NEC   NEC SXF90 compiler directives.
!
!       !/OMPH  Ading OMP directves for hybrid paralellization.
!
!       !/S     Enable subroutine tracing.
!       !/Tn    Enable test output.
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3QCK1 (MX, MY, NX, NY, CFLL, Q, CLOSE, INC,         &
                         MAPACT, NACT, MAPBOU, NB0, NB1, NB2,         &
                         NDSE, NDST )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         30-Oct-2009 |
!/                  +-----------------------------------+
!/
!/    11-Mar-1997 : Final FORTRAN 77                    ( version 1.18 )
!/    15-Dec-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    15-Feb-2001 : Unit numbers added to par list.     ( version 2.08 )
!/    05-Mar-2008 : Added NEC sxf90 compiler directives.
!/                  (Chris Bunney, UK Met Office)       ( version 3.13 )
!/    30-Oct-2009 : Fixed "Called by" doc line.         ( version 3.14 )
!/                  (T. J. Campbell, NRL)
!/
!  1. Purpose :
!
!     Preform one-dimensional propagation in a two-dimensional space
!     with irregular boundaries and regular grid.
!
!  2. Method :
!
!     ULTIMATE QUICKEST scheme (see manual).
!
!     Note that the check on monotonous behavior of QCN is performed
!     using weights CFAC, to avoid the need for IF statements.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       MX,MY   Int.   I   Field dimensions, if grid is 'closed' or
!                          circular, MX is the closed dimension.
!       NX,NY   Int.   I   Part of field actually used.
!       CFLL    R.A.   I   Local Courant numbers.         (MY,  MX+1)
!       Q       R.A.  I/O  Propagated quantity.           (MY,0:MX+2)
!       CLOSE   Log.   I   Flag for closed 'X' dimension'
!       INC     Int.   I   Increment in 1-D array corresponding to
!                          increment in 2-D space.
!       MAPACT  I.A.   I   List of active grid points.
!       NACT    Int.   I   Size of MAPACT.
!       MAPBOU  I.A.   I   Map with boundary information (see W3MAP2).
!       NBn     Int.   I   Counters in MAPBOU.
!       NDSE    Int.   I   Error output unit number.
!       NDST    Int.   I   Test output unit number.
!     ----------------------------------------------------------------
!           - CFLL amd Q need only bee filled in the (MY,MX) range,
!             extension is used internally for closure.
!           - CFLL and Q are defined as 1-D arrays internally.
!
!  4. Subroutines used :
!
!       STRACE   Service routine.
!
!  5. Called by :
!
!       W3KTP2   Propagation in spectral space
!
!  6. Error messages :
!
!     None.
!
!  7. Remarks :
!
!     - This routine can be used independently from WAVEWATCH III.
!
!  8. Structure :
!
!     ------------------------------------------------------
!       1. Initialize aux. array FLA.
!       2. Fluxes for central points (3rd order + limiter).
!       3. Fluxes boundary point above (1st order).
!       4. Fluxes boundary point below (1st order).
!       5. Closure of 'X' if required
!       6. Propagate.
!     ------------------------------------------------------
!
!  9. Switches :
!
!     !/S   Enable subroutine tracing.
!     !/T   Enable test output.
!     !/T0  Test output input/output fields.
!     !/T1  Test output fluxes.
!     !/T2  Test output integration.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: MX, MY, NX, NY, INC, MAPACT(MY*MX),  &
                                 NACT, MAPBOU(MY*MX), NB0, NB1, NB2,  &
                                 NDSE, NDST
      REAL, INTENT(INOUT)     :: CFLL(MY*(MX+1)), Q(1-MY:MY*(MX+2))
      LOGICAL, INTENT(IN)     :: CLOSE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IXY, IP, IXYC, IXYU, IXYD, IY, IX,   &
                                 IAD00, IAD02, IADN0, IADN1, IADN2
      REAL                    :: CFL, QB, DQ, DQNZ, QCN, QBN, QBR, CFAC
      REAL                    :: FLA(1-MY:MY*MX)
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Initialize aux. array FLA and closure ------------------------- *
!
      FLA = 0.
!
      IF ( CLOSE ) THEN
          IAD00  = -MY
          IAD02  =  MY
          IADN0  = IAD00 + MY*NX
          IADN1  =         MY*NX
          IADN2  = IAD02 + MY*NX
          DO IY=1, NY
            Q   (IY+IAD00) = Q   (IY+IADN0)
            Q   (IY+IADN1) = Q   (   IY   )
            Q   (IY+IADN2) = Q   (IY+IAD02)
            CFLL(IY+IADN1) = CFLL(   IY   )
            END DO
        END IF
!
! 2.  Fluxes for central points ------------------------------------- *
!     ( 3rd order + limiter )
!
      DO IP=1, NB0
!
        IXY    = MAPBOU(IP)
        CFL    = 0.5 * ( CFLL(IXY) + CFLL(IXY+INC) )
        IXYC   = IXY - INC * INT( MIN ( 0. , SIGN(1.1,CFL) ) )
        QB     = 0.5 * ( (1.-CFL)*Q(IXY+INC) + (1.+CFL)*Q(IXY) )      &
            - (1.-CFL**2)/6. * (Q(IXYC-INC)-2.*Q(IXYC)+Q(IXYC+INC))
!
        IXYU   = IXYC - INC * INT ( SIGN (1.1,CFL) )
        IXYD   = 2*IXYC - IXYU
        DQ     = Q(IXYD) - Q(IXYU)
        DQNZ   = SIGN ( MAX(1.E-15,ABS(DQ)) , DQ )
        QCN    = ( Q(IXYC) - Q(IXYU) ) / DQNZ
        QCN    = MIN ( 1.1, MAX ( -0.1 , QCN ) )
!
        QBN    = MAX ( (QB-Q(IXYU))/DQNZ , QCN )
        QBN    = MIN ( QBN , 1. , QCN/MAX(1.E-10,ABS(CFL)) )
        QBR    = Q(IXYU) + QBN*DQ
        CFAC   = REAL ( INT( 2. * ABS(QCN-0.5) ) )
        QB     = (1.-CFAC)*QBR + CFAC*Q(IXYC)
!
        FLA(IXY) = CFL * QB
!
        END DO
!
! 3.  Fluxes for points with boundary above ------------------------- *
!     ( 1st order without limiter )
!
      DO IP=NB0+1, NB1
        IXY    = MAPBOU(IP)
        CFL    = CFLL(IXY)
        IXYC   = IXY - INC * INT( MIN ( 0. , SIGN(1.1,CFL) ) )
        FLA(IXY) = CFL * Q(IXYC)
        END DO
!
! 4.  Fluxes for points with boundary below ------------------------- *
!     ( 1st order without limiter )
!
      DO IP=NB1+1, NB2
        IXY    = MAPBOU(IP)
        CFL    = CFLL(IXY+INC)
        IXYC   = IXY - INC * INT( MIN ( 0. , SIGN(1.1,CFL) ) )
        FLA(IXY) = CFL * Q(IXYC)
        END DO
!
! 5.  Global closure ----------------------------------------------- *
!
      IF ( CLOSE ) THEN
          DO IY=1, NY
            FLA (IY+IAD00) = FLA (IY+IADN0)
            END DO
        END IF
!
! 6.  Propagation -------------------------------------------------- *
!
      DO IP=1, NACT
        IXY    = MAPACT(IP)
        Q(IXY) = MAX ( 0. , Q(IXY) + FLA(IXY-INC) - FLA(IXY) )
        END DO
!
      RETURN
!
! Formats
!
!/
!/ End of W3QCK1 ----------------------------------------------------- /
!/
      END SUBROUTINE W3QCK1
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3QCK2 (MX, MY, NX, NY, VELO, DT, DX1, DX2, Q, CLOSE,&
                        INC,  MAPACT, NACT, MAPBOU, NB0, NB1, NB2,    &
                        NDSE, NDST )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         30-Oct-2009 |
!/                  +-----------------------------------+
!/
!/    07-Sep-1997 : Final FORTRAN 77                    ( version 1.18 )
!/    16-Dec-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    14-Feb-2001 : Unit numbers added to par list.     ( version 2.08 )
!/    05-Mar-2008 : Added NEC sxf90 compiler directives.
!/                  (Chris Bunney, UK Met Office)       ( version 3.13 )
!/    30-Oct-2009 : Fixed "Called by" doc line.         ( version 3.14 )
!/                  (T. J. Campbell, NRL)
!/
!  1. Purpose :
!
!     Like W3QCK1 with variable grid spacing.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       MX,MY   Int.   I   Field dimensions, if grid is 'closed' or
!                          circular, MX is the closed dimension.
!       NX,NY   Int.   I   Part of field actually used.
!       VELO    R.A.   I   Local velocities.              (MY,  MX+1)
!       DT      Real   I   Time step.
!       DX1     R.A.  I/O  Band width at points.          (MY,  MX+1)
!       DX2     R.A.  I/O  Band width between points.     (MY,0:MX+1)
!                          (local counter and counter+INC)
!       Q       R.A.  I/O  Propagated quantity.           (MY,0:MX+2)
!       CLOSE   Log.   I   Flag for closed 'X' dimension'
!       INC     Int.   I   Increment in 1-D array corresponding to
!                          increment in 2-D space.
!       MAPACT  I.A.   I   List of active grid points.
!       NACT    Int.   I   Size of MAPACT.
!       MAPBOU  I.A.   I   Map with boundary information (see W3MAP2).
!       NBn     Int.   I   Counters in MAPBOU.
!       NDSE    Int.   I   Error output unit number.
!       NDST    Int.   I   Test output unit number.
!     ----------------------------------------------------------------
!           - VELO amd Q need only bee filled in the (MY,MX) range,
!             extension is used internally for closure.
!           - VELO and Q are defined as 1-D arrays internally.
!
!  4. Subroutines used :
!
!       STRACE   Service routine.
!
!  5. Called by :
!
!       W3KTP2   Propagation in spectral space
!
!  6. Error messages :
!
!     None.
!
!  7. Remarks :
!
!     - This routine can be used independently from WAVEWATCH III.
!
!  8. Structure :
!
!     ------------------------------------------------------
!       1. Initialize aux. array FLA.
!       2. Fluxes for central points (3rd order + limiter).
!       3. Fluxes boundary point above (1st order).
!       4. Fluxes boundary point below (1st order).
!       5. Closure of 'X' if required
!       6. Propagate.
!     ------------------------------------------------------
!
!  9. Switches :
!
!     !/S   Enable subroutine tracing.
!     !/T   Enable test output.
!     !/T0  Test output input/output fields.
!     !/T1  Test output fluxes.
!     !/T2  Test output integration.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: MX, MY, NX, NY, INC, MAPACT(MY*MX),  &
                                 NACT, MAPBOU(MY*MX), NB0, NB1, NB2,  &
                                 NDSE, NDST
      REAL, INTENT(IN)        :: DT
      REAL, INTENT(INOUT)     :: VELO(MY*(MX+1)), DX1(MY*(MX+1)),     &
                                 DX2(1-MY:MY*(MX+1)), Q(1-MY:MY*(MX+2))
      LOGICAL, INTENT(IN)     :: CLOSE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IXY, IP, IXYC, IXYU, IXYD, IY, IX,   &
                                 IAD00, IAD02, IADN0, IADN1, IADN2
      REAL                    :: CFL, VEL, QB, DQ, DQNZ, QCN, QBN,    &
                                 QBR, CFAC, FLA(1-MY:MY*MX)
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Initialize aux. array FLA and closure ------------------------- *
!
      FLA = 0.
!
      IF ( CLOSE ) THEN
          IAD00  = -MY
          IAD02  =  MY
          IADN0  = IAD00 + MY*NX
          IADN1  =         MY*NX
          IADN2  = IAD02 + MY*NX
          DO IY=1, NY
            Q   (IY+IAD00) = Q   (IY+IADN0)
            Q   (IY+IADN1) = Q   (   IY   )
            Q   (IY+IADN2) = Q   (IY+IAD02)
            VELO(IY+IADN1) = VELO(   IY   )
            DX1 (IY+IADN1) = DX1 (   IY   )
            DX2 (IY+IAD00) = DX1 (IY+IADN0)
            DX2 (IY+IADN1) = DX1 (   IY   )
            END DO
        END IF
!
! 2.  Fluxes for central points ------------------------------------- *
!     ( 3rd order + limiter )
!
      DO IP=1, NB0
!
        IXY    = MAPBOU(IP)
        VEL    = 0.5 * ( VELO(IXY) + VELO(IXY+INC) )
        CFL    = DT *  VEL / DX2(IXY)
        IXYC   = IXY - INC * INT( MIN ( 0. , SIGN(1.1,CFL) ) )
        QB     = 0.5 * ( (1.-CFL)*Q(IXY+INC) + (1.+CFL)*Q(IXY) )      &
                    - DX2(IXY)**2 / DX1(IXYC) * (1.-CFL**2) / 6.      &
                        * ( (Q(IXYC+INC)-Q(IXYC))/DX2(IXYC)           &
                          - (Q(IXYC)-Q(IXYC-INC))/DX2(IXYC-INC) )
!
        IXYU   = IXYC - INC * INT ( SIGN (1.1,CFL) )
        IXYD   = 2*IXYC - IXYU
        DQ     = Q(IXYD) - Q(IXYU)
        DQNZ   = SIGN ( MAX(1.E-15,ABS(DQ)) , DQ )
        QCN    = ( Q(IXYC) - Q(IXYU) ) / DQNZ
        QCN    = MIN ( 1.1, MAX ( -0.1 , QCN ) )
!
        QBN    = MAX ( (QB-Q(IXYU))/DQNZ , QCN )
        QBN    = MIN ( QBN , 1. , QCN/MAX(1.E-10,ABS(CFL)) )
        QBR    = Q(IXYU) + QBN*DQ
        CFAC   = REAL ( INT( 2. * ABS(QCN-0.5) ) )
        QB     = (1.-CFAC)*QBR + CFAC*Q(IXYC)
!
        FLA(IXY) = VEL * QB
!
        END DO
!
! 3.  Fluxes for points with boundary above ------------------------- *
!     ( 1st order without limiter )
!
      DO IP=NB0+1, NB1
        IXY    = MAPBOU(IP)
        VEL    = VELO(IXY)
        IXYC   = IXY - INC * INT( MIN ( 0. , SIGN(1.1,VEL) ) )
        FLA(IXY) = VEL * Q(IXYC)
        END DO
!
! 4.  Fluxes for points with boundary below ------------------------- *
!     ( 1st order without limiter )
!
      DO IP=NB1+1, NB2
        IXY    = MAPBOU(IP)
        VEL    = VELO(IXY+INC)
        IXYC   = IXY - INC * INT( MIN ( 0. , SIGN(1.1,VEL) ) )
        FLA(IXY) = VEL * Q(IXYC)
        END DO
!
! 5.  Global closure ----------------------------------------------- *
!
      IF ( CLOSE ) THEN
          DO IY=1, NY
            FLA (IY+IAD00) = FLA (IY+IADN0)
            END DO
        END IF
!
! 6.  Propagation -------------------------------------------------- *
!
      DO IP=1, NACT
        IXY    = MAPACT(IP)
        Q(IXY) = MAX ( 0. , Q(IXY) + DT/DX1(IXY) *                    &
                              (FLA(IXY-INC)-FLA(IXY)) )
        END DO
!
      RETURN
!
! Formats
!
!/
!/ End of W3QCK2 ----------------------------------------------------- /
!/
      END SUBROUTINE W3QCK2
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3QCK3 (MX, MY, NX, NY, CFLL, TRANS, Q, CLOSE,       &
                         INC, MAPACT, NACT, MAPBOU, NB0, NB1, NB2,    &
                         NDSE, NDST )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         27-May-2014 |
!/                  +-----------------------------------+
!/
!/    13_nov-2001 : Origination.                        ( version 2.14 )
!/    16-Oct-2002 : Fix INTENT for TRANS.               ( version 3.00 )
!/    05-Mar-2008 : Added NEC sxf90 compiler directives.
!/                  (Chris Bunney, UK Met Office)       ( version 3.13 )
!/    27-May-2014 : Added OMPH switches in W3QCK3.      ( version 5.02 )
!/
!  1. Purpose :
!
!     Like W3QCK1 with cell transparencies added.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       MX,MY   Int.   I   Field dimensions, if grid is 'closed' or
!                          circular, MX is the closed dimension.
!       NX,NY   Int.   I   Part of field actually used.
!       CFLL    R.A.   I   Local Courant numbers.         (MY,  MX+1)
!       Q       R.A.  I/O  Propagated quantity.           (MY,0:MX+2)
!       CLOSE   Log.   I   Flag for closed 'X' dimension'
!       INC     Int.   I   Increment in 1-D array corresponding to
!                          increment in 2-D space.
!       MAPACT  I.A.   I   List of active grid points.
!       NACT    Int.   I   Size of MAPACT.
!       MAPBOU  I.A.   I   Map with boundary information (see W3MAP2).
!       NBn     Int.   I   Counters in MAPBOU.
!       NDSE    Int.   I   Error output unit number.
!       NDST    Int.   I   Test output unit number.
!     ----------------------------------------------------------------
!           - CFLL amd Q need only bee filled in the (MY,MX) range,
!             extension is used internally for closure.
!           - CFLL and Q are defined as 1-D arrays internally.
!
!  4. Subroutines used :
!
!       STRACE   Service routine.
!
!  5. Called by :
!
!       W3XYP2   Propagation in physical space
!
!  6. Error messages :
!
!     None.
!
!  7. Remarks :
!
!     - This routine can be used independently from WAVEWATCH III.
!
!  8. Structure :
!
!     ------------------------------------------------------
!       1. Initialize aux. array FLA.
!       2. Fluxes for central points (3rd order + limiter).
!       3. Fluxes boundary point above (1st order).
!       4. Fluxes boundary point below (1st order).
!       5. Closure of 'X' if required
!       6. Propagate.
!     ------------------------------------------------------
!
!  9. Switches :
!
!     !/OMPH  Ading OMP directves for hybrid paralellization.
!
!     !/S   Enable subroutine tracing.
!     !/T   Enable test output.
!     !/T0  Test output input/output fields.
!     !/T1  Test output fluxes.
!     !/T2  Test output integration.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: MX, MY, NX, NY, INC, MAPACT(MY*MX),  &
                                 NACT, MAPBOU(MY*MX), NB0, NB1, NB2,  &
                                 NDSE, NDST
      REAL, INTENT(IN)        :: TRANS(MY*MX,-1:1)
      REAL, INTENT(INOUT)     :: CFLL(MY*(MX+1)), Q(1-MY:MY*(MX+2))
      LOGICAL, INTENT(IN)     :: CLOSE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IXY, IP, IXYC, IXYU, IXYD, IY, IX,   &
                                 IAD00, IAD02, IADN0, IADN1, IADN2,   &
                                 JN, JP
      REAL                    :: CFL, QB, DQ, DQNZ, QCN, QBN, QBR, CFAC
      REAL                    :: FLA(1-MY:MY*MX)
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Initialize aux. array FLA and closure ------------------------- *
!
      FLA = 0.
!
      IF ( CLOSE ) THEN
          IAD00  = -MY
          IAD02  =  MY
          IADN0  = IAD00 + MY*NX
          IADN1  =         MY*NX
          IADN2  = IAD02 + MY*NX
!
          DO IY=1, NY
            Q   (IY+IAD00) = Q   (IY+IADN0) ! 1 ghost column to left
            Q   (IY+IADN1) = Q   (   IY   ) ! 1st ghost column to right
            Q   (IY+IADN2) = Q   (IY+IAD02) ! 2nd ghost column to right
            CFLL(IY+IADN1) = CFLL(   IY   ) ! as for Q above, 1st to rt
            END DO
!
        END IF
!
! 2.  Fluxes for central points ------------------------------------- *
!     ( 3rd order + limiter )
!
      DO IP=1, NB0
!
        IXY    = MAPBOU(IP)
        CFL    = 0.5 * ( CFLL(IXY) + CFLL(IXY+INC) )
        IXYC   = IXY - INC * INT( MIN ( 0. , SIGN(1.1,CFL) ) )
        QB     = 0.5 * ( (1.-CFL)*Q(IXY+INC) + (1.+CFL)*Q(IXY) )      &
            - (1.-CFL**2)/6. * (Q(IXYC-INC)-2.*Q(IXYC)+Q(IXYC+INC))
!
        IXYU   = IXYC - INC * INT ( SIGN (1.1,CFL) )
        IXYD   = 2*IXYC - IXYU
        DQ     = Q(IXYD) - Q(IXYU)
        DQNZ   = SIGN ( MAX(1.E-15,ABS(DQ)) , DQ )
        QCN    = ( Q(IXYC) - Q(IXYU) ) / DQNZ
        QCN    = MIN ( 1.1, MAX ( -0.1 , QCN ) )
!
        QBN    = MAX ( (QB-Q(IXYU))/DQNZ , QCN )
        QBN    = MIN ( QBN , 1. , QCN/MAX(1.E-10,ABS(CFL)) )
        QBR    = Q(IXYU) + QBN*DQ
        CFAC   = REAL ( INT( 2. * ABS(QCN-0.5) ) )
        QB     = (1.-CFAC)*QBR + CFAC*Q(IXYC)
!
        FLA(IXY) = CFL * QB
!
        END DO
!
! 3.  Fluxes for points with boundary above ------------------------- *
!     ( 1st order without limiter )
!
!!!/OMPH/!$OMP PARALLEL DO PRIVATE (IP, IXY, CFL, IXYC)
!!!
      DO IP=NB0+1, NB1
        IXY    = MAPBOU(IP)
        CFL    = CFLL(IXY)
        IXYC   = IXY - INC * INT( MIN ( 0. , SIGN(1.1,CFL) ) )
        FLA(IXY) = CFL * Q(IXYC)
        END DO
!!!
!!!/OMPH/!$OMP END PARALLEL DO
!
! 4.  Fluxes for points with boundary below ------------------------- *
!     ( 1st order without limiter )
!
!!!/OMPH/!$OMP PARALLEL DO PRIVATE (IP, IXY, CFL, IXYC)
!!!
      DO IP=NB1+1, NB2
        IXY    = MAPBOU(IP)
        CFL    = CFLL(IXY+INC)
        IXYC   = IXY - INC * INT( MIN ( 0. , SIGN(1.1,CFL) ) )
        FLA(IXY) = CFL * Q(IXYC)
        END DO
!
!!!/OMPH/!$OMP END PARALLEL DO
!
! 5.  Global closure ----------------------------------------------- *
!
      IF ( CLOSE ) THEN
          DO IY=1, NY
            FLA (IY+IAD00) = FLA (IY+IADN0)
            END DO
        END IF
!
! 6.  Propagation -------------------------------------------------- *
!
      DO IP=1, NACT
!
        IXY    = MAPACT(IP)
        IF ( FLA(IXY-INC) .GT. 0. ) THEN
            JN    = -1
          ELSE
            JN    =  0
          END IF
        IF ( FLA(IXY    ) .LT. 0. ) THEN
            JP    =  1
          ELSE
            JP    =  0
          END IF
!
        Q(IXY) = MAX ( 0. , Q(IXY) + TRANS(IXY,JN) * FLA(IXY-INC)     &
                                   - TRANS(IXY,JP) * FLA(IXY) )
 
        END DO
!
      RETURN
!
! Formats
!
!/
!/ End of W3QCK3 ----------------------------------------------------- /
!/
      END SUBROUTINE W3QCK3
!/
!/ End of module W3UQCKMD -------------------------------------------- /
!/
      END MODULE W3UQCKMD
