!/ ------------------------------------------------------------------- /
      MODULE W3PARTMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III          USACE/NOAA |
!/                  |          Barbara  Tracy           |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         15-Apr-2008 |
!/                  +-----------------------------------+
!/
!/    01-Nov-2006 : Origination.                        ( version 3.10 )
!/    02-Nov-2006 : Adding tail to integration.         ( version 3.10 )
!/    24-Mar-2007 : Bug fix IMI, adding overall field   ( version 3.11 )
!/                  and sorting.
!/    15-Apr-2008 : Clean up for distribution.          ( version 3.14 )
!/
!  1. Purpose :
!
!     Spectral partitioning according to the watershed method.
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      MK, MTH   Int.  Private  Dimensions of stored neighour array.
!      NEIGH     I.A.  Private  Nearest Neighbor array.
!     ----------------------------------------------------------------
!      Note: IHMAX, HSPMIN, WSMULT, WSCUT and FLCOMB used from W3ODATMD.
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3PART    Subr. Public   Interface to watershed routines.
!      PTSORT    Subr. Public   Sort discretized image.
!      PTNGHB    Subr. Public   Defeine nearest neighbours.
!      PT_FLD    Subr. Public   Incremental flooding algorithm.
!      FIFO_ADD, FIFO_EMPTY, FIFO_FIRST
!                Subr. PT_FLD   Queue management.
!      PTMEAN    Subr. Public   Compute mean parameters.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine traceing.
!      WAVNU1    Subr. W3DISPMD Wavenumber computation.
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!  6. Switches :
!
!     !/S    Enable subroutine tracing.
!     !/T    Enable test output
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE W3ODATMD, ONLY: IHMAX, HSPMIN, WSMULT
!
      PUBLIC
!
      INTEGER, PRIVATE              :: MK = -1, MTH = -1
      INTEGER, ALLOCATABLE, PRIVATE :: NEIGH(:,:)
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3PART ( SPEC, UABS, UDIR, DEPTH, WN, NP, XP, DIMXP )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III          USACE/NOAA |
!/                  |          Barbara  Tracy           |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         28-Oct-2006 !
!/                  +-----------------------------------+
!/
!/    28-Oct-2006 : Origination.                       ( version 3.10 )
!/
!  1. Purpose :
!
!     Interface to watershed partitioning routines.
!
!  2. Method :
!
!     Watershed Algorithm of Vincent and Soille, 1991, implemented by
!     Barbara Tracy (USACE/ERDC) for NOAA/NCEP.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       SPEC    R.A.   I   2-D spectrum E(f,theta).
!       UABS    Real   I   Wind speed.
!       UDIR    Real   I   Wind direction.
!       DEPTH   Real   I   Water depth.
!       WN      R.A.   I   Wavenumebers for each frequency.
!       NP      Int.   O   Number of partitions.
!                           -1 : Spectrum without minumum energy.
!                            0 : Spectrum with minumum energy.
!                                but no partitions.
!       XP      R.A.   O   Parameters describing partitions.
!                          Entry '0' contains entire spectrum.
!       DIMXP   Int.   I   Second dimension of XP.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Sur.  W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!  6. Error messages :
!
!  7. Remarks :
!
!     - To achieve minimum storage but guaranteed storage of all
!       partitions DIMXP = ((NK+1)/2) * ((NTH-1)/2)
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
!/
      USE W3CONSTANTS
!
      USE W3GDATMD, ONLY: NK, NTH, NSPEC
      USE W3ODATMD, ONLY: WSCUT, FLCOMB
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(OUT)          :: NP
      INTEGER, INTENT(IN)           :: DIMXP
      REAL, INTENT(IN)              :: SPEC(NK,NTH), WN(NK), UABS,    &
                                       UDIR, DEPTH
      REAL, INTENT(OUT)             :: XP(6,0:DIMXP)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: ITH, IMI(NSPEC), IMD(NSPEC),         &
                                 IMO(NSPEC), IND(NSPEC), NP_MAX,      &
                                 IP, IT(1), INDEX(DIMXP), NWS,        &
                                 IPW, IPT, ISP, PMAP(DIMXP)
      REAL                    :: ZP(NSPEC), ZMIN, ZMAX, Z(NSPEC),     &
                                 FACT, WSMAX, HSMAX
      REAL                    :: TP(6,DIMXP)
!/
!/ ------------------------------------------------------------------- /
! 0.  Initializations
!
      NP     = 0
      XP     = 0.
!
! -------------------------------------------------------------------- /
! 1.  Process input spectrum
! 1.a 2-D to 1-D spectrum
!
      DO ITH=1, NTH
        ZP(1+(ITH-1)*NK:ITH*NK) = SPEC(:,ITH)
        END DO
!
! 1.b Invert spectrum and 'digitize'
!
      ZMIN   = MINVAL ( ZP )
      ZMAX   = MAXVAL ( ZP )
      IF ( ZMAX-ZMIN .LT. 1.E-9 ) RETURN
!
      Z      = ZMAX - ZP
!
      FACT   = REAL(IHMAX-1) / ( ZMAX - ZMIN )
      IMI    = MAX ( 1 , MIN ( IHMAX , NINT ( 1. + Z*FACT ) ) )
!
! 1.c Sort digitized image
!
      CALL PTSORT ( IMI, IND, IHMAX )
!
! -------------------------------------------------------------------- /
! 2.  Perform partitioning
! 2.a Update nearest neighbor info as needed.
!
      CALL PTNGHB
!
! 2.b Incremental flooding
!
      CALL PT_FLD ( IMI, IND, IMO, ZP, NP_MAX )
!
! 2.c Compute parameters per partition
!     NP and NX initialized inside routine.
!
      CALL PTMEAN(NP_MAX,IMO,ZP,DEPTH,UABS,UDIR,WN,NP,XP,DIMXP,PMAP)
!
! -------------------------------------------------------------------- /
! 3.  Sort and recombine wind seas as needed
! 3.a Sort by wind sea fraction
!
      IF ( NP .LE. 1 ) RETURN
!
      TP(:,1:NP)  = XP(:,1:NP)
      XP(:,1:NP)  = 0.
      INDEX(1:NP) = 0
      NWS         = 0
!
      DO IP=1, NP
        IT          = MAXLOC(TP(6,1:NP))
        INDEX(IP)   = IT(1)
        XP(:,IP)    = TP(:,INDEX(IP))
        IF ( TP(6,IT(1)) .GE. WSCUT ) NWS = NWS + 1
        TP(6,IT(1)) = -1.
        END DO
!
! 3.b Combine wind seas as needed and resort
!
      IF ( NWS.GT.1 .AND. FLCOMB ) THEN
          IPW    = PMAP(INDEX(1))
          DO IP=2, NWS
            IPT    = PMAP(INDEX(IP))
            DO ISP=1, NSPEC
              IF ( IMO(ISP) .EQ. IPT ) IMO(ISP) = IPW
              END DO
            END DO
!
          CALL PTMEAN(NP_MAX,IMO,ZP,DEPTH,UABS,UDIR,WN,NP,XP,DIMXP,PMAP)
          IF ( NP .LE. 1 ) RETURN
!
          TP(:,1:NP)  = XP(:,1:NP)
          XP(:,1:NP)  = 0.
          INDEX(1:NP) = 0
          NWS         = 0
!
          DO IP=1, NP
            IT          = MAXLOC(TP(6,1:NP))
            INDEX(IP)   = IT(1)
            XP(:,IP)    = TP(:,INDEX(IP))
            IF ( TP(6,IT(1)) .GE. WSCUT ) NWS = NWS + 1
            TP(6,IT(1)) = -1.
            END DO
!
        END IF
!
! 3.c Sort remaining fields by wave height
!
      NWS    = MIN ( 1 , NWS )
!
      TP(:,1:NP)  = XP(:,1:NP)
      XP(:,1:NP)  = 0.
!
      IF ( NWS .GT. 0 ) THEN
          XP(:,1) = TP(:,1)
          TP(1,1) = -1.
          NWS     = 1
        END IF
!
      DO IP=NWS+1, NP
        IT          = MAXLOC(TP(1,1:NP))
        XP(:,IP)    = TP(:,IT(1))
        TP(1,IT(1)) = -1.
        END DO
!
! -------------------------------------------------------------------- /
! 4.  End of routine
!
      RETURN
!/
!/ End of W3PART ----------------------------------------------------- /
!/
      END SUBROUTINE W3PART
!/ ------------------------------------------------------------------- /
      SUBROUTINE PTSORT ( IMI, IND, IHMAX )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III          USACE/NOAA |
!/                  |          Barbara  Tracy           |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         19-Oct-2006 !
!/                  +-----------------------------------+
!/
!/    19-Oct-2006 : Origination.                       ( version 3.10 )
!/
!  1. Purpose :
!
!     This subroutine sorts the image data in ascending order.
!     This sort original to F.T.Tracy (2006)
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMI     I.A.   I   Input discretized spectrum.
!       IND     I.A.   O   Sorted data.
!       IHMAX   Int.   I   Number of integer levels.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Sur.  W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE W3GDATMD, ONLY: NSPEC
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)      :: IHMAX, IMI(NSPEC)
      INTEGER, INTENT(OUT)     :: IND(NSPEC)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: I, IN, IV
      INTEGER                 :: NUMV(IHMAX), IADDR(IHMAX),           &
                                 IORDER(NSPEC)
!/
!
! -------------------------------------------------------------------- /
! 1.  Occurences per height
!
      NUMV   = 0
      DO I=1, NSPEC
        NUMV(IMI(I)) = NUMV(IMI(I)) + 1
        END DO
!
! -------------------------------------------------------------------- /
! 2.  Starting address per height
!
      IADDR(1) = 1
      DO I=1, IHMAX-1
        IADDR(I+1) = IADDR(I) + NUMV(I)
      END DO
!
! -------------------------------------------------------------------- /
! 3.  Order points
!
      DO I=1, NSPEC
        IV        = IMI(I)
        IN        = IADDR(IV)
        IORDER(I) = IN
        IADDR(IV) = IN + 1
        END DO
!
! -------------------------------------------------------------------- /
! 4.  Sort points
!
      DO I=1, NSPEC
        IND(IORDER(I)) = I
        END DO
!
      RETURN
!/
!/ End of PTSORT ----------------------------------------------------- /
!/
      END SUBROUTINE PTSORT
!/ ------------------------------------------------------------------- /
      SUBROUTINE PTNGHB
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III          USACE/NOAA |
!/                  |          Barbara  Tracy           |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         20-Oct-2006 !
!/                  +-----------------------------------+
!/
!/    20-Oct-2006 : Origination.                       ( version 3.10 )
!/
!  1. Purpose :
!
!     This subroutine computes the nearest neighbors for each grid
!     point. Wrapping of directional distribution (0 to 360)is taken
!     care of using the nearest neighbor system
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMI     I.A.   I   Input discretized spectrum.
!       IMD     I.A.   O   Sorted data.
!       IHMAX   Int.   I   Number of integer levels.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Sur.  W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE W3GDATMD, ONLY: NK, NTH, NSPEC
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!     INTEGER, INTENT(IN)      :: IHMAX, IMI(NSPEC)
!     INTEGER, INTENT(IN)      :: IMD(NSPEC)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: N, J, I, K
!/
!
! -------------------------------------------------------------------- /
! 1.  Check on need of processing
!
      IF ( MK.EQ.NK .AND. MTH.EQ.NTH ) RETURN
!
      IF ( MK.GT.0 ) DEALLOCATE ( NEIGH )
      ALLOCATE ( NEIGH(9,NSPEC) )
      MK     = NK
      MTH    = NTH
!
! -------------------------------------------------------------------- /
! 2.  Build map
!
      NEIGH  = 0
!
! ... Base loop
!
      DO N = 1, NSPEC
!
        J      = (N-1) / NK + 1
        I      = N - (J-1) * NK
        K      = 0
!
! ... Point at the left(1)
!
        IF ( I .NE. 1 ) THEN
            K           = K + 1
            NEIGH(K, N) = N - 1
          END IF
!
! ... Point at the right (2)
!
        IF ( I .NE. NK ) THEN
            K           = K + 1
            NEIGH(K, N) = N + 1
          END IF
!
! ... Point at the bottom(3)
!
        IF ( J .NE. 1 ) THEN
            K           = K + 1
            NEIGH(K, N) = N - NK
          END IF
!
! ... ADD Point at bottom_wrap to top
!
        IF ( J .EQ. 1 ) THEN
            K          = K + 1
            NEIGH(K,N) = NSPEC - (NK-I)
          END IF
!
! ... Point at the top(4)
!
        IF ( J .NE. NTH ) THEN
            K           = K + 1
            NEIGH(K, N) = N + NK
          END IF
!
! ... ADD Point to top_wrap to bottom
!
         IF ( J .EQ. NTH ) THEN
             K          = K + 1
             NEIGH(K,N) = N - (NTH-1) * NK
            END IF
!
! ... Point at the bottom, left(5)
!
        IF ( (I.NE.1) .AND. (J.NE.1) ) THEN
            K           = K + 1
            NEIGH(K, N) = N - NK - 1
          END IF
!
! ... Point at the bottom, left with wrap.
!
         IF ( (I.NE.1) .AND. (J.EQ.1) ) THEN
             K          = K + 1
             NEIGH(K,N) = N - 1 + NK * (NTH-1)
           END IF
!
! ... Point at the bottom, right(6)
!
        IF ( (I.NE.NK) .AND. (J.NE.1) ) THEN
            K           = K + 1
            NEIGH(K, N) = N - NK + 1
          END IF
!
! ... Point at the bottom, right with wrap
!
        IF ( (I.NE.NK) .AND. (J.EQ.1) ) THEN
            K           = K + 1
            NEIGH(K,N) = N + 1 + NK * (NTH - 1)
          END  IF
!
! ... Point at the top, left(7)
!
        IF ( (I.NE.1) .AND. (J.NE.NTH) ) THEN
            K           = K + 1
            NEIGH(K, N) = N + NK - 1
          END IF
!
! ... Point at the top, left with wrap
!
         IF ( (I.NE.1) .AND. (J.EQ.NTH) ) THEN
             K           = K + 1
             NEIGH(K,N) = N - 1 - (NK) * (NTH-1)
           END IF
!
! ... Point at the top, right(8)
!
        IF ( (I.NE.NK) .AND. (J.NE.NTH) ) THEN
            K           = K + 1
            NEIGH(K, N) = N + NK + 1
          END IF
!
! ... Point at top, right with wrap
!
        IF ( (I.NE.NK) .AND. (J.EQ.NTH) ) THEN
            K           = K + 1
            NEIGH(K,N) = N + 1 - (NK) * (NTH-1)
          END IF
!
        NEIGH(9,N) = K
!
        END DO
!
      RETURN
!/
!/ End of PTNGHB ----------------------------------------------------- /
!/
      END SUBROUTINE PTNGHB
!/ ------------------------------------------------------------------- /
      SUBROUTINE PT_FLD ( IMI, IND, IMO, ZP, NPART )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Nov-2006 !
!/                  +-----------------------------------+
!/
!/    01-Nov-2006 : Origination.                       ( version 3.10 )
!/
!  1. Purpose :
!
!     This subroutine does incremental flooding of the image to
!     determine the watershed image.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMI     I.A.   I   Input discretized spectrum.
!       IND     I.A.   I   Sorted addresses.
!       IMO     I.A.   O   Output partitioned spectrum.
!       ZP      R.A.   I   Spectral array.
!       NPART   Int.   O   Number of partitions found.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Sur.  W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE W3GDATMD, ONLY: NSPEC
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: IMI(NSPEC), IND(NSPEC)
      INTEGER, INTENT(OUT)    :: IMO(NSPEC), NPART
      REAL, INTENT(IN)        :: ZP(NSPEC)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: MASK, INIT, IWSHED, IMD(NSPEC),      &
                                 IC_LABEL, IFICT_PIXEL, M, IH, MSAVE, &
                                 IP, I, IPP, IC_DIST, IEMPTY, IPPP,   &
                                 JL, JN, IPT, J
      INTEGER                 :: IQ(NSPEC), IQ_START, IQ_END
      REAL                    :: ZPMAX, EP1, DIFF
!/
!
! -------------------------------------------------------------------- /
! 0.  Initializations
!
      MASK        = -2
      INIT        = -1
      IWSHED      =  0
      IMO         = INIT
      IC_LABEL    =  0
      IMD         =  0
      IFICT_PIXEL = -100
!
      IQ_START    =  1
      IQ_END      =  1
!
      ZPMAX       = MAXVAL ( ZP )
!
! -------------------------------------------------------------------- /
! 1.  Loop over levels
!
      M      =  1
!
      DO IH=1, IHMAX
        MSAVE  = M
!
! 1.a Pixels at level IH
!
        DO
          IP     = IND(M)
          IF ( IMI(IP) .NE. IH ) EXIT
!
!     Flag the point, if it stays flagge, it is a separate minimum.
!
          IMO(IP) = MASK
!
!     Consider neighbors. If there is neighbor, set distance and add
!     to queue.
!
          DO I=1, NEIGH(9,IP)
            IPP    = NEIGH(I,IP)
            IF ( (IMO(IPP).GT.0) .OR. (IMO(IPP).EQ.IWSHED) ) THEN
                IMD(IP) = 1
                CALL FIFO_ADD (IP)
                EXIT
              END IF
            END DO
!
          IF ( M+1 .GT. NSPEC ) THEN
              EXIT
            ELSE
              M = M + 1
            END IF
!
          END DO
!
! 1.b Process the queue
!
        IC_DIST = 1
        CALL FIFO_ADD (IFICT_PIXEL)
!
        DO
          CALL FIFO_FIRST (IP)
!
!     Check for end of processing
!
          IF ( IP .EQ. IFICT_PIXEL ) THEN
              CALL FIFO_EMPTY (IEMPTY)
              IF ( IEMPTY .EQ. 1 ) THEN
                  EXIT
                ELSE
                  CALL FIFO_ADD (IFICT_PIXEL)
                  IC_DIST = IC_DIST + 1
                  CALL FIFO_FIRST (IP)
                END IF
            END IF
!
!     Process queue
!
          DO I=1, NEIGH(9,IP)
            IPP = NEIGH(I,IP)
!
!     Check for labeled watersheds or basins
!
            IF ( (IMD(IPP).LT.IC_DIST) .AND. ( (IMO(IPP).GT.0) .OR.  &
                 (IMO(IPP).EQ.IWSHED))) THEN
!
                IF ( IMO(IPP) .GT. 0 ) THEN
!
                    IF ((IMO(IP) .EQ. MASK) .OR. (IMO(IP) .EQ. &
                        IWSHED)) THEN
                        IMO(IP) = IMO(IPP)
                      ELSE IF (IMO(IP) .NE. IMO(IPP)) THEN
                        IMO(IP) = IWSHED
                      END IF
!
                  ELSE IF (IMO(IP) .EQ. MASK) THEN
!
                    IMO(IP) = IWSHED
!
                  END IF
!
              ELSE IF ( (IMO(IPP).EQ.MASK) .AND. (IMD(IPP).EQ.0) ) THEN
!
                 IMD(IPP) = IC_DIST + 1
                 CALL FIFO_ADD (IPP)
!
              END IF
!
            END DO
!
          END DO
!
! 1.c Check for mask values in IMO to identify new basins
!
        M = MSAVE
!
        DO
          IP     = IND(M)
          IF ( IMI(IP) .NE. IH ) EXIT
          IMD(IP) = 0
!
          IF (IMO(IP) .EQ. MASK) THEN
!
! ... New label for pixel
!
              IC_LABEL = IC_LABEL + 1
              CALL FIFO_ADD (IP)
              IMO(IP) = IC_LABEL
!
! ... and all connected to it ...
!
              DO
                CALL FIFO_EMPTY (IEMPTY)
                IF ( IEMPTY .EQ. 1 ) EXIT
                CALL FIFO_FIRST (IPP)
!
                DO I=1, NEIGH(9,IPP)
                  IPPP   = NEIGH(I,IPP)
                  IF ( IMO(IPPP) .EQ. MASK ) THEN
                      CALL FIFO_ADD (IPPP)
                      IMO(IPPP) = IC_LABEL
                    END IF
                  END DO
!
                END DO
!
            END IF
!
          IF ( M + 1 .GT. NSPEC ) THEN
              EXIT
            ELSE
              M = M + 1
            END IF
!
          END DO
!
        END DO
!
! -------------------------------------------------------------------- /
! 2.  Find nearest neighbor of 0 watershed points and replace
!     use original input to check which group to affiliate with 0
!     Soring changes first in IMD to assure symetry in adjustment.
!
      DO J=1, 5
        IMD    = IMO
        DO JL=1 , NSPEC
          IPT    = -1
          IF ( IMO(JL) .EQ. 0 ) THEN
              EP1    = ZPMAX
              DO JN=1, NEIGH (9,JL)
                DIFF   = ABS ( ZP(JL) - ZP(NEIGH(JN,JL)))
                IF ( (DIFF.LE.EP1) .AND. (IMO(NEIGH(JN,JL)).NE.0) ) THEN
                    EP1    = DIFF
                    IPT    = JN
                  END IF
                END DO
              IF ( IPT .GT. 0 ) IMD(JL) = IMO(NEIGH(IPT,JL))
            END IF
          END DO
        IMO    = IMD
        IF ( MINVAL(IMO) .GT. 0 ) EXIT
        END DO
!
      NPART = IC_LABEL
!
      RETURN
!
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE FIFO_ADD ( IV )
!
!     Add point to FIFO queue.
!
      INTEGER, INTENT(IN)      :: IV
!
      IQ(IQ_END) = IV
!
      IQ_END = IQ_END + 1
      IF ( IQ_END .GT. NSPEC ) IQ_END = 1
!
      RETURN
      END SUBROUTINE
!/ ------------------------------------------------------------------- /
      SUBROUTINE FIFO_EMPTY ( IEMPTY )
!
!     Check if queue is empty.
!
      INTEGER, INTENT(OUT)     :: IEMPTY
!
      IF ( IQ_START .NE. IQ_END ) THEN
        IEMPTY = 0
      ELSE
        IEMPTY = 1
      END IF
!
      RETURN
      END SUBROUTINE
!/ ------------------------------------------------------------------- /
      SUBROUTINE FIFO_FIRST ( IV )
!
!     Get point out of queue.
!
      INTEGER, INTENT(OUT)     :: IV
!
      IV = IQ(IQ_START)
!
      IQ_START = IQ_START + 1
      IF ( IQ_START .GT. NSPEC ) IQ_START = 1
!
      RETURN
      END SUBROUTINE
!/
!/ End of PT_FLD ----------------------------------------------------- /
!/
      END SUBROUTINE PT_FLD
!/ ------------------------------------------------------------------- /
      SUBROUTINE PTMEAN ( NPI, IMO, ZP, DEPTH, UABS, UDIR, WN,        &
                          NPO, XP, DIMXP, PMAP )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III          USACE/NOAA |
!/                  |          Barbara  Tracy           |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         24-Mar-2007 !
!/                  +-----------------------------------+
!/
!/    28-Oct-2006 : Origination.                       ( version 3.10 )
!/    02-Nov-2006 : Adding tail to integration.        ( version 3.10 )
!/    24-Mar-2007 : Adding overall field.              ( version 3.11 )
!/
!  1. Purpose :
!
!     Compute pean parameters per partition.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NPI     Int.   I   Number of partitions found.
!       IMO     I.A.   I   Partition map.
!       ZP      R.A.   I   Input spectrum.
!       DEPTH   Real   I   Water depth.
!       UABS    Real   I   Wind speed.
!       UDIR    Real   I   Wind direction.
!       WN      R.A.   I   Wavenumebers for each frequency.
!       NPO     Int.   O   Number of partitions with mean parameters.
!       XP      R.A.   O   Array with output parameters.
!       DIMXP   int.   I   Second dimesion of XP.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Sur.  W3SERVMD Subroutine tracing.
!      WAVNU1    Subr. W3DISPMD Wavenumber computation.
!     ----------------------------------------------------------------
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE W3CONSTANTS
      USE W3DISPMD, ONLY: WAVNU1
!
      USE W3GDATMD, ONLY: NK, NTH, NSPEC, DTH, SIG, DSII, DSIP,       &
                          ECOS, ESIN, XFR, FACHFE, TH, FTE
      USE W3ODATMD, ONLY: IAPROC, NAPERR, NDSE, NDST
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: NPI, IMO(NSPEC), DIMXP
      INTEGER, INTENT(OUT)    :: NPO
      REAL, INTENT(IN)        :: ZP(NSPEC), DEPTH, UABS, UDIR, WN(NK)
      REAL, INTENT(OUT)       :: XP(6,0:DIMXP)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IK, ITH, ISP, IP, IFPMAX(0:NPI)
      INTEGER, INTENT(OUT)    :: PMAP(DIMXP)
      REAL                    :: SUMF(0:NK+1,0:NPI), SUMFW(NK,0:NPI), &
                                 SUMFX(NK,0:NPI), SUMFY(NK,0:NPI),    &
                                 SUME(0:NPI), SUMEW(0:NPI),           &
                                 SUMEX(0:NPI), SUMEY(0:NPI),          &
                                 EFPMAX(0:NPI), FCDIR(NTH)
      REAL                    :: HS, XL, XH, XL2, XH2, EL, EH, DENOM, &
                                 SIGP, WNP, CGP, UPAR, C(NK), RD, FACT
!/
!
! -------------------------------------------------------------------- /
! 1.  Check on need of processing
!
      NPO    = 0
      XP     = 0.
!
      IF ( NPI .EQ. 0 ) RETURN
!
! -------------------------------------------------------------------- /
! 2.  Initialize arrays
!
      SUMF   = 0.
      SUMFW  = 0.
      SUMFX  = 0.
      SUMFY  = 0.
      SUME   = 0.
      SUMEW  = 0.
      SUMEX  = 0.
      SUMEY  = 0.
      IFPMAX = 0
      EFPMAX = 0.
!
      DO IK=1, NK
        C(IK)  = SIG(IK) / WN(IK)
        END DO
!
      DO ITH=1, NTH
        UPAR   = WSMULT * UABS * MAX(0.,COS(TH(ITH)-DERA*UDIR))
        IF ( UPAR .LT. C(NK) ) THEN
            FCDIR(ITH) = SIG(NK+1)
          ELSE
            DO IK=NK-1, 2, -1
              IF ( UPAR .LT. C(IK) ) EXIT
              END DO
            RD     = (C(IK)-UPAR) / (C(IK)-C(IK+1))
            IF ( RD .LT. 0 ) THEN
                IK     = 0
                RD     = MAX ( 0., RD+1. )
              END IF
            FCDIR(ITH) = RD*SIG(IK+1) + (1.-RD)*SIG(IK)
          END IF
        END DO
!
! -------------------------------------------------------------------- /
! 3.  Spectral integrals and preps
! 3.a Integrals
!     NOTE: Factor DTH only used in Hs computation.
!
      DO IK=1, NK
        DO ITH=1, NTH
          ISP    = IK + (ITH-1)*NK
          IP     = IMO(ISP)
          FACT   = MAX ( 0. , MIN ( 1. ,                              &
            1. - ( FCDIR(ITH) - 0.5*(SIG(IK-1)+SIG(IK)) ) / DSIP(IK) ) )
          SUMF (IK, 0) = SUMF (IK, 0) + ZP(ISP)
          SUMFW(IK, 0) = SUMFW(IK, 0) + ZP(ISP) * FACT
          SUMFX(IK, 0) = SUMFX(IK, 0) + ZP(ISP) * ECOS(ITH)
          SUMFY(IK, 0) = SUMFY(IK, 0) + ZP(ISP) * ESIN(ITH)
          IF ( IP .EQ. 0 ) CYCLE
          SUMF (IK,IP) = SUMF (IK,IP) + ZP(ISP)
          SUMFW(IK,IP) = SUMFW(IK,IP) + ZP(ISP) * FACT
          SUMFX(IK,IP) = SUMFX(IK,IP) + ZP(ISP) * ECOS(ITH)
          SUMFY(IK,IP) = SUMFY(IK,IP) + ZP(ISP) * ESIN(ITH)
          END DO
        END DO
      SUMF(NK+1,:) = SUMF(NK,:) * FACHFE
!
      DO IP=0, NPI
        DO IK=1, NK
          SUME (IP) = SUME (IP) + SUMF (IK,IP) * DSII(IK)
          SUMEW(IP) = SUMEW(IP) + SUMFW(IK,IP) * DSII(IK)
          SUMEX(IP) = SUMEX(IP) + SUMFX(IK,IP) * DSII(IK)
          SUMEY(IP) = SUMEY(IP) + SUMFY(IK,IP) * DSII(IK)
          IF ( SUMF(IK,IP) .GT. EFPMAX(IP) ) THEN
              IFPMAX(IP) = IK
              EFPMAX(IP) = SUMF(IK,IP)
            END IF
          END DO
        SUME (IP) = SUME (IP) + SUMF (NK,IP) * FTE
        SUMEW(IP) = SUMEW(IP) + SUMFW(NK,IP) * FTE
        SUMEX(IP) = SUMEX(IP) + SUMFX(NK,IP) * FTE
        SUMEY(IP) = SUMEY(IP) + SUMFY(NK,IP) * FTE
        END DO
!
! -------------------------------------------------------------------- /
! 4.  Compute pars
!
      NPO    = -1
!
      DO IP=0, NPI
!
        HS     = 4. * SQRT ( SUME(IP) * DTH * TPIINV )
        IF ( HS .LT. HSPMIN ) CYCLE
!
        XL     = 1./XFR - 1.
        XH     =  XFR - 1.
        XL2    = XL**2
        XH2    = XH**2
        EL     = SUMF(IFPMAX(IP)-1,IP) - SUMF(IFPMAX(IP),IP)
        EH     = SUMF(IFPMAX(IP)+1,IP) - SUMF(IFPMAX(IP),IP)
        DENOM  = XL*EH - XH*EL
        SIGP   = SIG(IFPMAX(IP)) * ( 1. + 0.5 * ( XL2*EH - XH2*EL )   &
                       / SIGN ( MAX(ABS(DENOM),1.E-15) , DENOM ) )
        CALL WAVNU1 ( SIGP, DEPTH, WNP, CGP )
!
        IF ( NPO .GE. DIMXP ) GOTO 2000
        NPO       = NPO + 1
        IF (IP.GT.0) PMAP(NPO) = IP
        XP(1,NPO) = HS
        XP(2,NPO) = TPI / SIGP
        XP(3,NPO) = TPI / WNP
        XP(4,NPO) = MOD( 630.-ATAN2(SUMEY(IP),SUMEX(IP))*RADE , 360. )
        XP(5,NPO) = RADE * SQRT ( MAX ( 0. , 2. * ( 1. - SQRT ( &
                MAX(0.,(SUMEX(IP)**2+SUMEY(IP)**2)/SUME(IP)**2) ) ) ) )
 
        XP(6,NPO) = SUMEW(IP) / SUME(IP)
!
        END DO
!
      RETURN
!
! Escape locations read errors --------------------------------------- *
!
 2000 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1000) NPO+1
      RETURN
!
! Formats
!
 1000 FORMAT (/' *** WAVEWATCH III ERROR IN PTMEAN :'/                &
               '     XP ARRAY TOO SMALL AT PARTITION',I6/)
!/
!/ End of PTMEAN ----------------------------------------------------- /
!/
      END SUBROUTINE PTMEAN
!/
!/ End of module W3PARTMD -------------------------------------------- /
!/
      END MODULE W3PARTMD
