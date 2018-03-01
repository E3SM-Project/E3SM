!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS
!-------------------------------------------------------------------------
!BOP
! !ROUTINE: ParPatternTest --- Unit tester for the parutilities patterns
!
! !INTERFACE:
      PROGRAM parpatterntest

! !USES:
      USE  decompmodule, ONLY: DecompType, DecompFree, DecompPermute,      &
     &     DecompCreate
      USE  parutilitiesmodule, ONLY: Gsize, GID, CommGlobal,               &
     &     ParPatternType, ParInit, ParExit, ParScatter, ParGather,        &
     &     ParPatternCreate, ParPatternFree,                               &
     &     ParBeginTransfer, ParEndTransfer
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
!    This main program tests the functionality of the decompmodule 
!    and parutilitiesmodule which relates to communication patterns.
!
!      Test 1: DecompRegular1D, ParPatternCreate, ParBegin/EndTransfer
!
!    Validation check:
!
!         mpirun -np 7 parpatterntest
!
!    Should yield a single message (if -DDEBUG_ON is *not* defined):
!
!         Passed all tests
!
! !REVISION HISTORY:
!   01.06.03   Sawyer     Creation from RedistributeTest
!   02.08.14   Sawyer     Now using explicit precision from pilgrim.h
!
!EOP
!-------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER  :: GlobalPEs, GlobalRank, I, J, Ierror
      TYPE (DecompType)  :: DecompA, DecompB, ObsDecomp
      TYPE (ParPatternType) :: InterAB, InterBA
      INTEGER  :: BlockLen, Remainder, Ipe
      REAL (CPP_REAL8) :: time1, time2, time3
      LOGICAL Passed


! For the Observation decomposition
      INTEGER  ::   Nactual, Xdim, Ydim, Zdim
      PARAMETER ( Nactual = 7 )
      PARAMETER ( Xdim = 144 )
      PARAMETER ( Ydim =  91 )
      PARAMETER ( Zdim =  30 )

      REAL (CPP_REAL8), ALLOCATABLE :: Rglobal(:), Rtmp(:), RlocalA(:), RlocalB(:)
      REAL (CPP_REAL8), ALLOCATABLE :: R2dlocalA(:), R2dlocalB(:)
      REAL (CPP_REAL8), ALLOCATABLE :: R2dGlobal(:), R2dtmp(:)
!
! 3D Arrays are smashed down to 1D
!
      REAL (CPP_REAL8), ALLOCATABLE :: R3dlocalA(:), R3dlocalB(:)
      REAL (CPP_REAL8), ALLOCATABLE :: R3dGlobal(:), R3dTmp(:)

      INTEGER , ALLOCATABLE :: itmp(:), DistA(:), DistB(:)
      INTEGER , ALLOCATABLE :: Xdist(:), Ydist(:), Zdist(:), Perm(:)

      CALL ParInit( )
      Passed = .TRUE.

!
! Initialize timing library.  2nd arg 0 means disable, 1 means enable
!
#if defined(TIMING)
      call t_setoptionf (gptlcpu, 1)
      call t_initializef ()
#endif

!
! Test 1 : Test DecompRegular1D, DecompPermute, ParScatter1D, ParGather1D 
!          and Redistribution using a block-wise distribution.
!

#if defined(TIMING)
      call t_startf('1D Redist Total')
#endif

!
!     Set the global vector to random values
!
      ALLOCATE( Rtmp(Nactual) )
      rtmp = 1.0_r8
!
! Decomposition for Observations:  Block distribution with remainder
! on last PE.  Should be OK if #obs >> #PEs
!
!
      ALLOCATE( DistA( Gsize ) )
      ALLOCATE( DistB( Gsize ) )
      ALLOCATE( Perm( Gsize ) )
      BlockLen = Nactual
      DO I = 1, Gsize-1
        DistA( I ) = BlockLen / 2
        BlockLen  = BlockLen - DistA(I)
      ENDDO
      DistA( Gsize ) = BlockLen
      IF ( SUM( DistA ) .ne. Nactual ) THEN
        print *, " Error: DistA contains ", SUM(DistA), " != ",Nactual
      ENDIF
      DO I = 1, Gsize
        DistB( I ) = DistA( Gsize-I+1 )
      ENDDO

      CALL DecompCreate ( Gsize, DistA, DecompA )
      CALL DecompCreate ( Gsize, DistB, DecompB )

      DO I=1, Gsize-1
        Perm( I+1 ) = I
      ENDDO
      Perm( 1 ) = Gsize
      CALL DecompPermute( Perm, DecompB )

      ALLOCATE( RlocalA( DecompA%NumEntries(GID+1) ) )
      ALLOCATE( RlocalB( DecompB%NumEntries(GID+1) ) )
      ALLOCATE( Rglobal( Nactual ) )

!!!      IF ( GID .EQ. 0 ) THEN
!!!        CALL RANDOM_NUMBER( HARVEST = Rglobal )
!!!      ENDIF
      DO i=1, Nactual
         Rglobal(i) = REAL(i,r8)
      ENDDO

!
! Now scatter the arrays over all PEs
!
      CALL ParScatter( CommGlobal, 0, Rglobal, DecompA, RlocalA )

!
! Now redistribute the local arrays from one decomposition to another
!
#if defined(TIMING)
      call t_startf('1D Redist Create')
#endif
      CALL ParPatternCreate( CommGlobal, DecompA, DecompB, InterAB )
      CALL ParPatternCreate( CommGlobal, DecompB, DecompA, InterBA )
#if defined(TIMING)
      call t_stopf('1D Redist Create')
#endif

#if defined(TIMING)
      call t_startf('1D Redist Forward')
#endif
      CALL ParBeginTransfer( InterAB, RlocalA, RlocalB )
      CALL ParEndTransfer( InterAB, RlocalA, RlocalB )
#if defined(TIMING)
      call t_stopf('1D Redist Forward')
#endif
      RlocalA = 0.0_r8

#if defined(TIMING)
      call t_startf('1D Redist Back')
#endif
      CALL ParBeginTransfer( InterBA, RlocalB, RlocalA )
      CALL ParEndTransfer( InterBA, RlocalB, RlocalA )
#if defined(TIMING)
      call t_stopf('1D Redist Back')
#endif
      CALL ParPatternFree( CommGlobal, InterBA )
      CALL ParPatternFree( CommGlobal, InterAB )

      CALL ParGather( CommGlobal, 0, RlocalA, DecompA, Rtmp )

      IF ( GID .eq. 0 ) THEN
        Rtmp = Rtmp - Rglobal
        IF ( SUM(Rtmp) .ne. 0.0_r8 ) THEN
          PRINT *, "Redistribution failed: 1D Gathered ver. != Orig."
          Passed = .FALSE.
        ENDIF
      ENDIF

      CALL DecompFree( DecompB )
      CALL DecompFree( DecompA )

      DEALLOCATE( DistB )
      DEALLOCATE( DistA )

      DEALLOCATE( RlocalB )
      DEALLOCATE( RlocalA )
      DEALLOCATE( Rtmp )
      DEALLOCATE( Rglobal )

#if defined(TIMING)
      call t_stopf('1D Redist Total')
#endif

!
! Test 2 : Test DecompRegular2D, ParScatter2D and ParGather2D 
!          and Redistribute using a 2-D block-wise distribution.
!
#if defined(TIMING)
      call t_startf('2D Redist Total')
#endif

!
! Set the target vector to non-random values
!
!
! Make sure that the array is not square
!  
      ALLOCATE( R2dtmp( XDim*YDim ) )

! Set the global vector to random values
! Make sure that the array is not square
!
      ALLOCATE( R2dGlobal( XDim*YDim ) )
      IF ( GID .EQ. 0 ) THEN
        CALL RANDOM_NUMBER( HARVEST = R2dglobal )
      ENDIF

!     Decomposition for Observations:  Block distribution with remainder
!       on last PE.  Should be OK if #obs >> #PEs
!

      ALLOCATE( DistA( Gsize ) )
      ALLOCATE( DistB( Gsize ) )
      ALLOCATE( XDist( 1 ) )
      ALLOCATE( YDist( 1 ) )

      BlockLen   = Xdim
      DO I = 1, Gsize-1
        DistA( I ) = BlockLen / 2
        BlockLen   = BlockLen - DistA( I ) 
      ENDDO
      DistA( Gsize ) = BlockLen
      YDist( 1 ) = Ydim
      
      XDist( 1 ) = Xdim
      BlockLen   = Ydim
      DO J = 1, Gsize-1
        DistB( J ) = BlockLen / 2
        BlockLen   = BlockLen - DistB( J ) 
      ENDDO
      DistB( Gsize ) = BlockLen

!
! Row-major ordering
!
      ALLOCATE( R2dlocalA( DistA(GID+1)*YDist(1) ) )
      ALLOCATE( R2dlocalB( XDist(1)*DistB(GID+1) ) )

      CALL DecompCreate( Gsize, 1, DistA, YDist, DecompA )
      CALL DecompCreate( 1, Gsize, Xdist, DistB, DecompB )

!
!     Now scatter the arrays over all PEs
!
      CALL ParScatter( CommGlobal, 0, R2dglobal, DecompA, R2dlocalA )

#if defined(TIMING)
      call t_startf('2D Redist Create')
#endif
      CALL ParPatternCreate( CommGlobal, DecompA, DecompB, InterAB )
      CALL ParPatternCreate( CommGlobal, DecompB, DecompA, InterBA )

#if defined(TIMING)
      call t_stopf('2D Redist Create')
#endif
#if defined(TIMING)
      call t_startf('2D Redist Forward')
#endif
      CALL ParBeginTransfer( InterAB, R2dlocalA, R2dlocalB )
      CALL ParEndTransfer( InterAB, R2dlocalA, R2dlocalB )
#if defined(TIMING)
      call t_stopf('2D Redist Forward')
#endif
      R2dlocalA = 0.0_r8
#if defined(TIMING)
      call t_startf('2D Redist Back')
#endif
      CALL ParBeginTransfer( InterBA, R2dlocalB, R2dlocalA )
      CALL ParEndTransfer( InterBA, R2dlocalB, R2dlocalA )
#if defined(TIMING)
      call t_stopf('2D Redist Back')
#endif

      CALL ParPatternFree( CommGlobal, InterAB )
      CALL ParPatternFree( CommGlobal, InterBA )

      CALL ParGather( CommGlobal, 0, R2dlocalA, DecompA, R2dtmp )
      IF ( GID .eq. 0 ) THEN
        R2dtmp = R2dtmp - R2dglobal
        IF ( SUM(R2dtmp) .ne. 0.0_r8 ) THEN
          PRINT *,"RedistributeTest Failed: 2D Gathered ver. != Orig."
          Passed = .FALSE.
        ENDIF
      ENDIF

      CALL DecompFree( DecompB )
      CALL DecompFree( DecompA )

      DEALLOCATE( R2dlocalB )
      DEALLOCATE( R2dlocalA )

      DEALLOCATE( YDist )
      DEALLOCATE( XDist )

      DEALLOCATE( DistB )
      DEALLOCATE( DistA )

      DEALLOCATE( R2dtmp )
      DEALLOCATE( R2dglobal )

#if defined(TIMING)
      call t_stopf('2D Redist Total')
#endif


!
! Test 3 : Test 3-D redistribution
!
#if defined(TIMING)
      call t_startf('3D Redist Total')
#endif

!
! Set the target vector to non-random values
!
!
! Make sure that the array is not square
!  
      ALLOCATE( R3dTmp( Xdim*Ydim*ZDim ) )

! Set the global vector to random values
! Make sure that the array is not square
!
      ALLOCATE( R3dGlobal( XDim*YDim*ZDim  ) )
      IF ( GID .eq. 0 ) THEN
        CALL RANDOM_NUMBER( HARVEST = R3dglobal )
      ENDIF
      r3dtmp = 1.0_r8


!
! Now define the distribution
!
      ALLOCATE( DistA( Gsize ) )
      ALLOCATE( DistB( Gsize ) )

      ALLOCATE( XDist( 1 ) )
      ALLOCATE( YDist( 1 ) )
      ALLOCATE( ZDist( 1 ) )

      XDist( 1 ) = Xdim
      YDist( 1 ) = Ydim
      ZDist( 1 ) = Zdim

!
! Optimal distribution in Z
!
      BlockLen   = Zdim / Gsize
      Remainder  = MOD( Zdim, Gsize ) 
 
      IF ( Remainder .gt. 0 ) DistA( 1:Remainder ) = BlockLen+1
      DistA( Remainder+1 : Gsize ) = BlockLen
      CALL DecompCreate( 1,1,Gsize,XDist,YDist,DistA,DecompA )

!
! Optimal distribution in Y
!
      BlockLen   = Ydim / Gsize
      Remainder  = MOD( Ydim, Gsize ) 
      IF ( Remainder .gt. 0 ) DistB( 1:Remainder ) = BlockLen+1
      DistB( Remainder+1 : Gsize ) = BlockLen
      CALL DecompCreate( 1,Gsize,1,XDist,DistB,ZDist,DecompB )

      ALLOCATE( R3dlocalA( XDist(1)*YDist(1)*DistA(GID+1) ) )
      ALLOCATE( R3dlocalB( XDist(1)*DistB(GID+1)*ZDist(1) ) )

!
! Do all the stuff here
!
      CALL ParScatter( CommGlobal, 0, R3dglobal, DecompA, R3dlocalA )

#if defined( TIMING )
      call t_startf('3D Redist Create')
#endif

      CALL ParPatternCreate( CommGlobal, DecompA, DecompB, InterAB )
      CALL ParPatternCreate( CommGlobal, DecompB, DecompA, InterBA )

#if defined(TIMING)
      call t_stopf('3D Redist Create')
#endif
#if defined(TIMING)
      call t_startf('3D Redist Forward')
#endif
      CALL ParBeginTransfer( InterAB, R3dlocalA, R3dlocalB )
      CALL ParEndTransfer( InterAB, R3dlocalA, R3dlocalB )
#if defined(TIMING)
      call t_stopf('3D Redist Forward')
#endif
      R3dlocalA = 0.0_r8
#if defined(TIMING)
      call t_startf('3D Redist Back')
#endif
      CALL ParBeginTransfer( InterBA, R3dlocalB, R3dlocalA )
      CALL ParEndTransfer( InterBA, R3dlocalB, R3dlocalA )
#if defined(TIMING)
      call t_stopf('3D Redist Back')
#endif
      CALL ParPatternFree( CommGlobal, InterAB )
      CALL ParPatternFree( CommGlobal, InterBA )

      CALL ParGather( CommGlobal, 0, R3dlocalA, DecompA, R3dtmp )

      IF ( GID .eq. 0 ) THEN
        R3dtmp = R3dtmp - R3dglobal
        IF ( SUM(R3dtmp) .ne. 0.0_r8 ) THEN
          PRINT *, "RedistributeTest failed: 3d Gathered ver. != Orig."
          Passed = .FALSE.
        ENDIF
      ENDIF

      CALL DecompFree( DecompB )
      CALL DecompFree( DecompA )

      DEALLOCATE( R3dlocalB )
      DEALLOCATE( R3dlocalA )

      DEALLOCATE( ZDist )
      DEALLOCATE( YDist )
      DEALLOCATE( XDist )

      DEALLOCATE( DistB )
      DEALLOCATE( DistA )

      DEALLOCATE( R3dtmp )
      DEALLOCATE( R3dglobal )

#if defined(TIMING)
      call t_stopf('3D Redist Total')
#endif

!
! That's all folks
!
#if defined(TIMING)
      call t_prf(GID)
#endif
      IF ( gid == 0 ) THEN
        IF ( Passed ) THEN
          PRINT *, "Passed ParPatternTest"
        ELSE
          PRINT *, "Failed ParPatternTest" 
        ENDIF
      ENDIF

      CALL ParExit( )

!EOC
!-------------------------------------------------------------------------
      END PROGRAM parpatterntest
