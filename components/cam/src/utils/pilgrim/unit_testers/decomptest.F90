!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS
!-------------------------------------------------------------------------
!BOP
! !ROUTINE: DecompTest --- Unit tester for the decomposition utilities
!
! !INTERFACE:
      PROGRAM decomptest

! !USES:
      USE  decompmodule

#include "pilgrim.h"
#include "debug.h"

      IMPLICIT NONE

! !DESCRIPTION:
!
!    This main program tests the functionality of the DecompModule
!    It performs the following tests:
!
!    \begin{enumerate}
!      \item DecompRegular1D
!      \item DecompRegular2D
!      \item DecompGlobalToLocal
!      \item DecompLocalToGlobal
!    \end{enumerate}
!
!    Validation check: ./DecompTest
!
!    Should yield a single message (if -DDEBUG_ON is *not* defined):
!
!      Passed all tests
!
!    Be patient, it may take 2 minutes.
!
! !LOCAL VARIABLES:
      TYPE (DecompType)  :: Decomp1d, Decomp2d, Decomp1dPerm

! For the Observation decomposition
      INTEGER  ::  NPEsComp, BlockLen, I, J, Local, Global, Pe, Local2, Pe2
      INTEGER  ::  Nactual, NPEsMax, Nx, Ny, Iglobal, Jglobal, Kglobal, K
      PARAMETER (Nactual = 131, NPEsMax = 4, Nx = 72, Ny = 46 )

      LOGICAL :: Passed
      REAL (CPP_REAL8), ALLOCATABLE :: Rtmp(:)
      INTEGER , ALLOCATABLE :: itmp(:), ilocal(:), Dist(:), Tags(:)
      INTEGER , ALLOCATABLE :: Xdist(:), Ydist(:), Perm(:)

! !REVISION HISTORY:
!   98.03.20   Sawyer     Creation
!   98.05.11   Sawyer     Added test of DecompCopy, DecompPermute
!   99.03.05   Sawyer     Renovated for complete unit test concept
!   01.02.07   Sawyer     Removed DG2L 2D test, added DecompCreate tests
!   01.05.01   Sawyer     free-format
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
      Passed = .TRUE.
      NPEsComp = 1
      DO WHILE( Passed .AND. NPEsComp .LE. NPEsMax )
!
! Test 1 : Test DecompRegular1D
!          using a block-wise distribution.
!
        ALLOCATE( Dist( NPEsComp ) )
!
! Decomposition for Observations:  Block distribution with remainder
! on last PE.  Should be OK if #obs >> #PEs
!
        BlockLen = Nactual
        DO I = 1, NPEsComp-1
          Dist( I ) = BlockLen / 2
          BlockLen  = BlockLen - Dist(I)
        ENDDO
        Dist( NPEsComp ) = BlockLen
        IF ( SUM( Dist ) .ne. Nactual ) THEN
          print *, "Error: Dist contains ", SUM(Dist), " != ",Nactual
        ENDIF
        CALL DecompCreate( NPEsComp, Dist, Decomp1D )

        DEALLOCATE( Dist )

        DO J = 1, Nactual
          CALL DecompGlobalToLocal( Decomp1D, J, Local, Pe )
          CALL DecompLocalToGlobal( Decomp1D, Local, Pe, Global )
          IF ( J .NE. Global ) THEN
            PRINT *, "DecompTest failed: 1D Global<->Local mapping: "
            PRINT *, "GlobalIn ", J, " = ( ", Local, ",", Pe, ")"
            PRINT *, "But: (", Local, ",", Pe, ") = ", Global
            Passed = .FALSE.
          ENDIF
        ENDDO
        CALL DecompFree( Decomp1D )

!
! Test 2 : Test DecompRegular2D
!
        ALLOCATE( Xdist( NPEsComp ) )
        ALLOCATE( Ydist( NPEsComp ) )
!
        BlockLen = Nactual
        DO I = 1, NPEsComp-1
          Xdist( I ) = BlockLen / 2
          Ydist( I ) = Nactual / NPEsComp
          BlockLen  = BlockLen - Xdist(I)
        ENDDO
        Xdist( NPEsComp ) = BlockLen
        Ydist( NPEsComp ) = Nactual - (NPEsComp-1)*(Nactual/NPEsComp)
        CALL DecompCreate( NPEsComp, NPEsComp, Xdist, Ydist, Decomp2D )
        DO J = 1, Nactual
          DO I = 1, Nactual
            K = (J-1)*Nactual + I
            CALL DecompGlobalToLocal( Decomp2D, K, Local, Pe )
            CALL DecompLocalToGlobal( Decomp2D, Local, Pe, Kglobal )
            Iglobal = MOD( Kglobal - 1, Nactual ) + 1
            Jglobal = ( Kglobal - 1 ) / Nactual + 1
            IF ( I .NE. Iglobal .OR. J .NE. Jglobal ) THEN
              PRINT *, "DecompTest failed: 2D Global<->Local mapping: "
              PRINT *, "( ",I,J," ) != ( ", Iglobal, Jglobal, ")"
              Passed = .FALSE.
            ENDIF
          ENDDO
        ENDDO

        DEALLOCATE( Ydist )
        DEALLOCATE( Xdist )
        CALL DecompFree( Decomp2D )

!
! Test 3 : Test DecompPermute
!
        ALLOCATE( Dist( NPEsComp ) )
!
! Decomposition for Observations:  Block distribution with remainder
! on last PE.  Should be OK if #obs >> #PEs  Same as Test 1
!
        BlockLen = Nactual
        DO I = 1, NPEsComp-1
          Dist( I ) = BlockLen / 2
          BlockLen  = BlockLen - Dist(I)
        ENDDO
        Dist( NPEsComp ) = BlockLen
        IF ( SUM( Dist ) .ne. Nactual ) THEN
          print *, " Error: Dist contains ", SUM(Dist), " != ",Nactual
        ENDIF
        CALL DecompCreate( NPEsComp, Dist, Decomp1D )
        
        DEALLOCATE( Dist )
!
! Copy and permute decomposition
!
        CALL DecompCopy( Decomp1d, Decomp1dPerm )
        ALLOCATE( Perm( NPEsComp ) )
        DO I = 1, NPEsComp
          Perm( NPEsComp - I + 1 ) = I
        ENDDO
        CALL DecompPermute( Perm, Decomp1dPerm )

!
! Run a simple test of the permutation
!
        DO J = 1, Nactual
          CALL DecompGlobalToLocal( Decomp1D, J, Local, Pe )
          CALL DecompGlobalToLocal( Decomp1DPerm, J, Local2, Pe2 )
          IF ( (Pe+1) .NE. Perm( Pe2+1 ) .OR. Local .NE. Local2 ) THEN
            PRINT *, "DecompTest failed, 1D permuted decomposition"
            PRINT *, "GlobalIn ", J, " = ( ", Local, ",", Pe, ")"
            PRINT *, "But permuted: (", Local2, ",", Perm(Pe2+1)-1, ")"
            Passed = .FALSE.
          ENDIF
        ENDDO
        CALL DecompFree( Decomp1D )
        DEALLOCATE( Perm )

!
!
! Test 4 : Test DecompCreate
!
        ALLOCATE( Tags( Nactual ) )
        ALLOCATE( Dist( Nactual ) )
        ALLOCATE( Rtmp( Nactual ) )
        ALLOCATE( Perm( NPEsComp ) )

!
! A random PE assignment is by far the hardest test for the library
!
        CALL RANDOM_NUMBER( HARVEST = Rtmp )
        Dist = INT( NPesComp*Rtmp - 0.5_r8 )
!
! This is the simple version of an irregular decomposition
!
        CALL DecompCreate( NPEsComp, Dist, Nactual, Decomp1D )
!
! Now some tests: basically go through all the local index to see
! if every global tag is accounted for
!
        Perm = 0
        Tags = 0
        DO I = 1, Nactual
          Perm( Dist(I) + 1 ) = Perm( Dist(I) + 1 ) + 1
        ENDDO
        DO pe=1,NPEsComp
          DO Local=1,Perm(pe)
            CALL DecompLocalToGlobal( Decomp1D, Local, Pe-1, Global )
            IF ( Tags( Global ) .NE. 0 ) THEN
              print *, "Error: DecompCreate"
              print *, "Local index",Local, Pe-1, "maps to", Global
              print *, "but", Global, "is taken by another index"
              Passed = .FALSE.
            ENDIF
          ENDDO
        ENDDO
          
!
! Now get trickier: define a unique, but not contiguous set of tags,
! for example a subset of 1..Nactual.  
! 
        CALL RANDOM_NUMBER( HARVEST = Rtmp )
        global = 0
        DO I=1, Nactual
          IF ( Rtmp(I) .GE. 0.3333_r8 .AND. Rtmp(I) .LT. 0.6667_r8 ) THEN
            global = global + 1
            Tags( global ) = I
          ENDIF
        ENDDO
!
        CALL RANDOM_NUMBER( HARVEST = Rtmp )
        Dist = INT( NPesComp*Rtmp - 0.5_r8 )

!
! This is the esoteric version of an irregular decomposition
!
        CALL DecompCreate( NPEsComp, Dist, Global, Tags, Decomp1Dperm )
!
! Now check that each of the active tags is properly defined
!
        K = 0
        DO i=1, Nactual
          CALL DecompGlobalToLocal( Decomp1Dperm, i, Local, Pe )
          IF ( Pe .NE. -1 ) THEN
            K = K + 1
            IF ( Dist( K ) .NE. Pe ) THEN
              print *, "Error DecompCreate test"
              print *, "Element", I,"on", Pe, "instead of", Dist(K)
              Passed = .FALSE.
            ENDIF
          ENDIF
        ENDDO
        IF ( K .NE. Global ) THEN
          print *, "Error: DecompCreate test"
          print *, "Found", K, "unique tags", "not correct", Global
          Passed = .FALSE.
        ENDIF

        DEALLOCATE( Perm )
        DEALLOCATE( Rtmp )
        DEALLOCATE( Dist )
        DEALLOCATE( Tags )

!
! Next PE configuration
!
        NPEsComp = NPEsComp * 2
      ENDDO

!
! That's all folks
!
      IF ( Passed ) THEN
        PRINT *, "Passed DecompTest"
      ELSE
        PRINT *, "Failed DecompTest"
      ENDIF

!EOC
!-------------------------------------------------------------------------
      END PROGRAM decomptest
