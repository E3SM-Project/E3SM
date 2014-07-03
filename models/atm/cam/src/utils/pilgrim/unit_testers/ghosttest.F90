!------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS
!------------------------------------------------------------------------
!BOP
! !ROUTINE: GhostTest --- Unit tester for the decomposition utilities
!
! !INTERFACE:
      PROGRAM GhostTest

! !USES:
      USE  decompmodule, ONLY: DecompType, DecompFree, DecompRegular1D,  &
     &              DecompRegular2D, DecompRegular3D, DecompCreate,      &
     &              DecompLocalToGlobal
      USE  ghostmodule, ONLY: GhostType, GhostFree, GhostCreate,         &
     &              GhostCopy, GhostInfo
      USE  parutilitiesmodule, ONLY: CommGlobal, GID, GSize,             &
     &              ParPatternType, ParPatternCreate, ParPatternFree,    &
     &              ParInit, ParExit, ParBeginTransfer, ParEndTransfer
#if defined(TIMING)
      USE  perf_mod
#endif

#include "debug.h"
#include "pilgrim.h"

      IMPLICIT NONE
#if defined(TIMING)
#include "gptl.inc"
#endif

! !DESCRIPTION:
!
!    This main program tests the functionality of the GhostModule
!    It performs the following tests:
!
!    \begin{enumerate}
!      \item 1D ghost region of a 1D decomposition
!      \item 2D ghost region of a 2D decomposition
!      \item 3D ghost region of a 3D decomposition
!      \item irregular ghost region of an irregular decomposition
!    \end{enumerate}
!
!    Validation check: 
!
!      mpirun -np 7 GhostTest
!
!    Should yield a single message (if -DDEBUG_ON is *not* defined):
!
!      Passed all tests
!
!    Be patient, it tests many complex cases, so it could take a while
!
! !REVISION HISTORY:
!   01.02.07   Sawyer     Creation
!   01.05.01   Sawyer     Minor changes for CCM framework
!   02.08.14   Sawyer     No uses explicit precisions
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      TYPE (DecompType)     :: Decomp
      TYPE (GhostType)      :: Ghost
      TYPE (ParPatternType) :: Pattern, Pattern2d, Pattern3d

      INTEGER  ::  Nactual, GhostWidth, Nx, Ny
      PARAMETER (Nactual = 100, GhostWidth = 4, Nx = 72, Ny = 46 )

! For the Observation decomposition
      INTEGER  ::  BlockLen, I, J, K, Local, Global, Pe
      INTEGER  ::  PEsInX, PEsInY, PEsInZ, IamInX, IamInY, IamInZ
      INTEGER  ::  Xstart, Xend, Ystart, Yend, Zstart, Zend, Ytrue, Ztrue

      LOGICAL :: Passed
      REAL (CPP_REAL8), ALLOCATABLE :: Rtmp(:), Rtmp2d(:,:), Rtmp3d(:,:,:)
      INTEGER , ALLOCATABLE :: itmp(:), ilocal(:), Dist(:), &
                      Tags(:), Xdist(:), Ydist(:), Zdist(:), Perm(:)

!
! GhostModule is communication-free, but in this test a communication
! pattern is constructed for different ghost regions.  This makes it
! an SPMD code.
!
      CALL ParInit()
      Passed = .TRUE.

!
! Initialize timing library.  2nd arg 0 means disable, 1 means enable
!
#if defined(TIMING)
      call t_setoptionf (gptlcpu, 1)
      call t_initializef ()
#endif

!
! Test 1 : Test GhostRegular1D
!          using a block-wise distribution.
!
#if defined(TIMING)
      call t_startf('1D Ghosting Total')
#endif

      ALLOCATE( Xdist( GSize ) )
!
! Decomposition for Observations:  Block distribution with remainder
! on last PE.  Should be OK if #obs >> #PEs
!
      Global = Nactual*Gsize
      BlockLen = Global
      DO I = 1, GSize-1
        Xdist( I ) = BlockLen / 2
        BlockLen  = BlockLen - Xdist(I)
      ENDDO
      Xdist( GSize ) = BlockLen
      CALL DecompRegular1D ( GSize, Xdist, Decomp )

!
! Now define a ghost region (i.e., a subset of the entire domain)
!
      Xstart = 1
      IF (GID .GT. 0) Xstart = SUM( Xdist(1:GID) ) + 1
      Xend = Global
      IF (GID .LT. Gsize-1) Xend   = Xstart + Xdist(GID+1) - 1
      DEALLOCATE( Xdist )

!
! Define ghost region with GhostWidth overlap (and wrap-around)
!
      CALL GhostCreate( Decomp, Gid, Global,                             &
     &               Xstart-GhostWidth, Xend+GhostWidth, .TRUE., Ghost )

! Allocate the ghosted region itself
      ALLOCATE( Rtmp( Xstart-GhostWidth:Xend+GhostWidth ) ) 

!
! Put the correct global tag into entry of the array, but zero out ghost region
!
      Rtmp = 0.0_r8
      DO I=Xstart, Xend
        Rtmp(I) = I
      ENDDO

!
! Now create a communication pattern which interrelates all the 
! ghosted vectors
!
#if defined(TIMING)
      call t_startf('1D PatternCreate')
#endif
      CALL ParPatternCreate( CommGlobal, Ghost, Pattern )
#if defined(TIMING)
      call t_stopf('1D PatternCreate')
#endif


!
! Do a test with the communication pattern
!
#if defined(TIMING)
      call t_startf('1D Ghost Transfer')
#endif
      CALL ParBeginTransfer( Pattern, Rtmp, Rtmp )
      CALL ParEndTransfer( Pattern, Rtmp, Rtmp )
#if defined(TIMING)
      call t_stopf('1D Ghost Transfer')
#endif
      DO I=Xstart-GhostWidth, Xend+GhostWidth
        IF ( Rtmp(I) .NE. MODULO(I-1,Global)+1 ) THEN
          print *, "Error on PE", GID, "Rtmp(",I,")=",Rtmp(I)
          Passed = .FALSE.
        ENDIF
      ENDDO

!
! Free the communication pattern
!
      DEALLOCATE( Rtmp )
      CALL ParPatternFree( CommGlobal, Pattern )
      CALL GhostFree( Ghost )
      CALL DecompFree( Decomp )

#if defined(TIMING)
      call t_stopf('1D Ghosting Total')
#endif

!
! Test 2 : Test DecompRegular2D
!

#if defined(TIMING)
      call t_startf('2D Ghosting Total')
#endif
      IF ( Gsize .GT. 1 ) THEN
        PEsInX = 2
        DO WHILE ( MOD(Gsize,PEsInX) .NE. 0 )
          PEsInX = PEsInX + 1
        ENDDO
      ELSE
          PEsInX = 1
      ENDIF
!
! In the worst case PEsInX = Gsize, PEsInY=1
!
      PEsInY = Gsize / PEsInX

      IamInY = GID / PEsInX
      IamInX = MOD( GID, PEsInX )

      ALLOCATE( Xdist( PEsInX ) )
      ALLOCATE( Ydist( PEsInY ) )
!
      BlockLen = Nactual
      DO I = 1, PEsInX-1
        Xdist( I ) = BlockLen / 2
        BlockLen  = BlockLen - Xdist(I)
      ENDDO
      Xdist( PEsInX ) = BlockLen

      DO J = 1, PEsInY-1
        Ydist( J ) = Nactual / PEsInY
      ENDDO
      Ydist( PEsInY ) = Nactual - (PEsInY-1)*(Nactual/PEsInY)

      CALL DecompRegular2D( PEsInX, PEsInY, Xdist, Ydist, Decomp )

      Xstart = 1
      IF (IamInX .GT. 0) Xstart = SUM( Xdist(1:IamInX) ) + 1
      Xend = Nactual
      IF (IamInX .LT. PEsInX-1) Xend   = Xstart + Xdist(IamInX+1) - 1
      Ystart = 1
      IF (IamInY .GT. 0) Ystart = SUM( Ydist(1:IamInY) ) + 1
      Yend = Nactual
      IF (IamInY .LT. PEsInY-1) Yend   = Ystart + Ydist(IamInY+1) - 1

      DEALLOCATE( Ydist )
      DEALLOCATE( Xdist )

!
! Now define a ghost region (i.e., a subset of the entire domain)
!
      CALL GhostCreate( Decomp, Gid,                                     &
     &              Nactual, Xstart-GhostWidth, Xend+GhostWidth,.TRUE.,  &
     &              Nactual, Ystart-GhostWidth, Yend+GhostWidth,.FALSE., &
     &              Ghost)

!
! Allocated the corresponding ghosted array: Note that some ghost regions
! will not be used (there is no wrap around)
!
      ALLOCATE( Rtmp2d(Xstart-GhostWidth:Xend+GhostWidth,                &
     &                 Ystart-GhostWidth:Yend+GhostWidth) )


!
! Put the correct global tag into entry of the array, but zero out ghost region
!
      Rtmp2d = 0.0_r8
      DO J=Ystart, Yend
        DO I=Xstart, Xend
          Rtmp2d(I,J) = (J-1)*Nactual + I
        ENDDO
      ENDDO

#if defined(TIMING)
      call t_startf('2D PatternCreate')
#endif
      CALL ParPatternCreate( CommGlobal, Ghost, Pattern2d )
#if defined(TIMING)
      call t_stopf('2D PatternCreate')
#endif

!
! Do a test with the communication pattern
!
#if defined(TIMING)
      call t_startf('2D Ghost Transfer')
#endif
      CALL ParBeginTransfer( Pattern2d, Rtmp2d, Rtmp2d )
      CALL ParEndTransfer( Pattern2d, Rtmp2d, Rtmp2d )
#if defined(TIMING)
      call t_stopf('2D Ghost Transfer')
#endif

      DO J=Ystart, Yend
        Ytrue = MODULO(J-1,Nactual)
        DO I=Xstart-GhostWidth, Xend+GhostWidth
          Global = Ytrue*Nactual + MODULO(I-1,Nactual) + 1
          IF ( Rtmp2D(I,J) .NE. Global ) THEN
            print *, "Error on PE", GID, "Rtmp2d(",I,J,")=",Rtmp2d(I,J)
            Passed = .FALSE.
          ENDIF
        ENDDO
      ENDDO

!
! Free the communication pattern
!
      CALL ParPatternFree( CommGlobal, Pattern2d )
      DEALLOCATE( Rtmp2D )
      CALL GhostFree( Ghost )
      CALL DecompFree( Decomp )

#if defined(TIMING)
      call t_stopf('2D Ghosting Total')
#endif

#if 0
!
! Test 3 : Test DecompRegular3D
!
#if defined(TIMING)
      call t_startf('3D Ghosting Total')
#endif
!
! In the case of a prime: PEsInZ = Gsize, PEsInY=1, PEsInX=1
!
      IF ( Gsize .GT. 1 ) THEN
        PEsInZ = 2
        DO WHILE ( MOD(Gsize,PEsInZ) .NE. 0 )
          PEsInZ = PEsInZ + 1
        ENDDO
      ELSE
        PEsInZ = 1
      ENDIF
      Pe = Gsize / PEsInZ

      IF ( Pe .GT. 1 ) THEN
        PEsInY = 2
        DO WHILE ( MOD(Pe,PEsInY) .NE. 0 )
          PEsInY = PEsInY + 1
        ENDDO
      ELSE
        PEsInY = 1
      ENDIF
!
      PEsInX = Pe / PEsInY
!
      IamInX = MOD( GID, PEsInX )
      IamInY = MOD( GID/PEsInX, PEsInY )
      IamInZ = GID / Pe

      ALLOCATE( Xdist( PEsInX ) )
      ALLOCATE( Ydist( PEsInY ) )
      ALLOCATE( Zdist( PEsInZ ) )
!
      BlockLen = Nactual
      DO I = 1, PEsInX-1
        Xdist( I ) = BlockLen / 2
        BlockLen  = BlockLen - Xdist(I)
      ENDDO
      Xdist( PEsInX ) = BlockLen

      DO J = 1, PEsInY-1
        Ydist( J ) = Nactual / PEsInY
      ENDDO
      Ydist( PEsInY ) = Nactual - (PEsInY-1)*(Nactual/PEsInY)

      BlockLen = Nactual
      DO K = PEsInZ,2,-1
        Zdist( K ) = BlockLen / 2
        BlockLen  = BlockLen - Zdist(K)
      ENDDO
      Zdist( 1 ) = BlockLen

      CALL DecompRegular3D( PEsInX, PEsInY, PEsInZ,                      &
     &                      Xdist, Ydist, Zdist, Decomp )

      Xstart = 1
      IF (IamInX .GT. 0) Xstart = SUM( Xdist(1:IamInX) ) + 1
      Xend = Nactual
      IF (IamInX .LT. PEsInX-1) Xend   = Xstart + Xdist(IamInX+1) - 1
      Ystart = 1
      IF (IamInY .GT. 0) Ystart = SUM( Ydist(1:IamInY) ) + 1
      Yend = Nactual
      IF (IamInY .LT. PEsInY-1) Yend   = Ystart + Ydist(IamInY+1) - 1
      Zstart = 1
      IF (IamInZ .GT. 0) Zstart = SUM( Zdist(1:IamInZ) ) + 1
      Zend = Nactual
      IF (IamInZ .LT. PEsInZ-1) Zend   = Zstart + Zdist(IamInZ+1) - 1

      DEALLOCATE( Zdist )
      DEALLOCATE( Ydist )
      DEALLOCATE( Xdist )

#if defined(DEBUG_GHOSTTEST)
      print *, GID, "Xstart", Xstart, "Xend", Xend, "Ystart", Ystart,    &
     &         "Yend", Yend, "Zstart", Zstart, "Zend", Zend,             &
     &         "IamInX", IamInX, "IamInY", IamInY, "IamInZ", IamInZ
#endif

!
! Now define a ghost region (i.e., a subset of the entire domain)
!
      CALL GhostCreate( Decomp, Gid,                                     &
     &              Nactual, Xstart-GhostWidth, Xend+GhostWidth,.FALSE., &
     &              Nactual, Ystart-GhostWidth, Yend+GhostWidth,.FALSE., &
     &              Nactual, Zstart-GhostWidth, Zend+GhostWidth,.TRUE.,  &
     &              Ghost)

!
! Allocated the corresponding ghosted array: Note that some ghost regions
! will not be used (there is no wrap around)
!
      ALLOCATE( Rtmp3d(Xstart-GhostWidth:Xend+GhostWidth,                &
     &                 Ystart-GhostWidth:Yend+GhostWidth,                &
     &                 Zstart-GhostWidth:Zend+GhostWidth) )

!
! Put the correct global tag into entry of the array, but zero out ghost region
!
      Rtmp3d = 0.0_r8
      DO K=Zstart, Zend
        DO J=Ystart, Yend
          DO I=Xstart, Xend
            Rtmp3d(I,J,K) = (K-1)*Nactual*Nactual + (J-1)*Nactual + I
          ENDDO
        ENDDO
      ENDDO
#if defined(TIMING)
      call t_startf('3D PatternCreate')
#endif
      CALL ParPatternCreate( CommGlobal, Ghost, Pattern3d )
#if defined(TIMING)
      call t_stopf('3D PatternCreate')
#endif

!
! Do a test with the communication pattern
!
#if defined(TIMING)
      call t_startf('3D Ghost Transfer')
#endif
      CALL BeginTransfer( Pattern3d, Rtmp3d, Rtmp3d )
      CALL EndTransfer( Pattern3d, Rtmp3d, Rtmp3d )
#if defined(TIMING)
      call t_stopf('3D Ghost Transfer')
#endif

      DO K=Zstart-GhostWidth, Zend+GhostWidth
        Ztrue = MODULO(K-1,Nactual)
        DO J=Ystart, Yend
          Ytrue = MODULO(J-1,Nactual)
          DO I=Xstart, Xend
            Global = (Ztrue*Nactual+Ytrue)*Nactual+MODULO(I-1,Nactual)+1
            IF ( Rtmp3D(I,J,K) .NE. Global ) THEN
              print *, "Error on",GID,"Rtmp3d(",I,J,K,")=",Rtmp3d(I,J,K)
              Passed = .FALSE.
            ENDIF
          ENDDO
        ENDDO
      ENDDO

!
! Free the communication pattern
!
      CALL ParPatternFree( CommGlobal, Pattern3d )
      CALL GhostFree( Ghost )
      CALL DecompFree( Decomp )
#if defined(TIMING)
      call t_stopf('3D Ghosting Total')
#endif

!
! Test 4 : Test Irregular Decomposition
!
      ALLOCATE( Tags( Nactual ) )
      ALLOCATE( Dist( Nactual ) )
      ALLOCATE( Rtmp( Nactual ) )
      ALLOCATE( Perm( GSize ) )

!
! A random PE assignment is by far the hardest test for the library
!
      CALL RANDOM_NUMBER( HARVEST = Rtmp )
      Dist = INT( GSize*Rtmp - 0.5_r8 )
!
! This is the simple version of an irregular decomposition
!
      CALL DecompCreate( GSize, Dist, Nactual, Decomp )

!
! Define the Ghost region through an arbitrary set of unique tags
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
      Dist = INT( GSize*Rtmp - 0.5_r8 )

!
! This is the esoteric version of an irregular decomposition
!
      CALL GhostCreate( Decomp, Gid, Global, Tags, Ghost )
      DEALLOCATE( Perm )
      DEALLOCATE( Rtmp )
      DEALLOCATE( Dist )
      DEALLOCATE( Tags )

      CALL DecompFree( Decomp )

      CALL ParPatternCreate( CommGlobal, Ghost, Pattern )
      CALL GhostFree( Ghost )

!
! Do a test with the communication pattern
!

!
! Free the communication pattern
!
      CALL ParPatternFree( CommGlobal, Pattern )
#endif

!
! That's all folks
!
#if defined(TIMING)
      call t_prf(GID)
#endif

      IF ( Passed ) THEN
        PRINT *, "Passed GhostTest"
      ELSE
        PRINT *, "Failed GhostTest"
      ENDIF

      CALL ParExit()

!EOC
!-------------------------------------------------------------------------
      END PROGRAM GhostTest

