!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS
!-------------------------------------------------------------------------
!BOP
! !ROUTINE: ParUtilitiesTest --- Unit tester for the parallel utilities
!
! !INTERFACE:
      PROGRAM parutilitiestest

! !USES:
#include "pilgrim.h"
      USE  decompmodule, ONLY: DecompType, DecompFree, DecompCreate
      USE  parutilitiesmodule

      IMPLICIT NONE

! !DESCRIPTION:
!
!    This main program tests the functionality of the ParUtilitites
!    module.  It performs the following tests:
!
!      Test 1: ParSplit, DecompRegular3D
!
!      Test 2: DecompRegular2D, ParScatter and ParGather
!         
!      Test 3: ParExchangeVector   
!
!    Validation check:
!
!         mpirun -np 7 ParUtilitiesTest
!
!    Should yield a single message (if -DDEBUG_ON is *not* defined):
!
!         Passed all tests
!
! !REVISION HISTORY:
!   00.07.20   Sawyer     Creation, from GEOS3_DAS_CORE version
!   00.08.21   Sawyer     Tests for ParCollective, ParGather/Scatter
!   01.05.01   Sawyer     free format, new decompmodule interfaces
!   02.08.14   Sawyer     Added explicit precisions from pilgrim.h
!
!EOP
!-------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER  :: GlobalPEs, GlobalRank, I, J, K, Ierror
      TYPE (DecompType)  :: Y3D, Z3D, YZ3D, XY2D, XY3D, ObsDecomp
      INTEGER  :: BlockLen, Ipe, Comm_1, Comm_2,npr_1,npr_2
      INTEGER  :: myid_1, myid_2, size_1, size_2, rank_1, rank_2
      LOGICAL Passed

! For the 2D and 3D decompositions
      INTEGER  ::   Im, Jm, Km
      PARAMETER ( Im = 72, Jm = 46, Km = 18 )
      REAL (CPP_REAL8)  ::     Tolerance

      REAL (CPP_REAL8), ALLOCATABLE :: Aglobal2d(:), Bglobal2d(:)
      REAL (CPP_REAL8), ALLOCATABLE :: Aglobal3d(:), Bglobal3d(:)
      REAL (CPP_REAL8), ALLOCATABLE :: Rsemiglobal3d(:)
      REAL (CPP_REAL8), ALLOCATABLE :: Rlocal2d(:), Rlocal3d(:)
      REAL (CPP_REAL8), ALLOCATABLE :: Rtmp(:), Rlocal(:), Rglobal(:)
      REAL (CPP_REAL8) :: Scalar, Q, Array1D(Im), Array2D(Im,Jm), Array3D(Im,Jm,Km)
      INTEGER  :: Inc, pe, IScalar, IArray1D(Im)
      INTEGER , ALLOCATABLE :: itmp(:),Dist(:),Xdist(:),Ydist(:),Zdist(:)


!
! Test 0 : Try to initialize the PEs
!
      CALL ParInit( )
      Passed = .TRUE.
      Tolerance = Gsize*1.0E-15_r8

!
! Create a virtual 2-D PE mesh.  Remember that this is not 
! inherently supported by PILGRIM which considers the PEs 
! as an agglommeration of Gsize processes
!
      npr_1 = 1
      npr_2 = 1
! 
! The following loop is guaranteed to terminate when npr_1 = Gsize
! and sooner if it finds a factorizaton of Gsize.  For best results,
! use Gsize = product of two primes
!
      DO WHILE ( npr_1 * npr_2 .LT. Gsize )
        npr_1 = npr_1 + 1
        npr_2 = Gsize / npr_1  
      ENDDO

      myid_2 = gid / npr_1
      myid_1 = gid - myid_2 * npr_1

#if !defined(USE_MLP)
!
! Test 1 : Test ParSplit, DecompRegular2D, DecompRegular3D, ParScatter,
!          and ParGather, as they will be used in the LLNL
!          2D decomposition FVCCM.  This test is supported in 
!          MPI PILGRIM but not in MLP (Jim Taft) since it requires
!          communicators and a group-wise barrier to be supported 
!          via shared memory.   This will come in due time.
!


!
! Split the communicators as LLNL needs
!
      call parsplit(commglobal, myid_2, gid, Comm_1, rank_1, size_1 )
      call parsplit(commglobal, myid_1, gid, Comm_2, rank_2, size_2 )

      IF ( myid_1 /= rank_1 .OR. myid_2 /= rank_2 .OR.                   &
     &     npr_1 /= size_1 .OR. npr_2 /= size_2 ) THEN
        print *, "ERROR in ParSplit: ranks or sizes are incorrect"
      ENDIF

      ALLOCATE( Dist( 1 ) )

!
! Create a latitude/level distribution which is intentionally unbalanced
! -- this is a tougher test of the software than a trivial case
!
      ALLOCATE( XDist( 1 ) )
      ALLOCATE( YDist( npr_1 ) )
      ALLOCATE( ZDist( npr_2 ) )

      Xdist(1) = Im
      BlockLen = Jm
      DO I = 1, npr_1-1
        YDist( I ) = BlockLen / 2
        BlockLen   = BlockLen - YDist( I ) 
      ENDDO
      YDist( npr_1 ) = BlockLen

      BlockLen = Km
      DO J = 1, npr_2-1
        ZDist( J ) = BlockLen / 2
        BlockLen   = BlockLen - ZDist( J ) 
      ENDDO
      ZDist( npr_2 ) = BlockLen

!
! Latitude/vertical decompositions 3D array
!
      ALLOCATE( Rlocal3D( im * ydist(myid_1+1) * zdist(myid_2+1) ) )
      call decompcreate( 1, npr_1, npr_2, xdist, ydist, zdist, YZ3D )

!
! Latitude strip decompositions 3D array
!
      dist(1) = km
      call decompcreate( 1, npr_1, 1, xdist, ydist, dist, Y3D )

!
! Now define a 3D decomposition local to the vertical column as 
! defined by LLNL.  This is slightly tricky, 
! since the decomposition is consistent only over
! all processors with the same myid_1
!

      ALLOCATE( Rsemiglobal3D( Im*ydist(myid_1+1)*Km ) )
      Dist(1) = ydist(myid_1+1)   ! DecompRegular3D requires arrays
      call decompcreate( 1,1,npr_2,xdist,dist,zdist, Z3D )

      DEALLOCATE( zdist )
      DEALLOCATE( ydist )
      DEALLOCATE( xdist )
      DEALLOCATE(  dist )

!
! Initialize the global arrays
!
      ALLOCATE( Bglobal3D( im*jm*km ) )
      ALLOCATE( Aglobal3D( im*jm*km ) )
      CALL RANDOM_NUMBER( HARVEST = Aglobal3D )  ! Only PE 0 is of interest

      CALL ParScatter( CommGlobal, 0, Aglobal3D, YZ3D, Rlocal3D )
      CALL ParGather(  Comm_2, 0, Rlocal3D, Z3D, Rsemiglobal3D )
!
! LLNL will want to do the following on only the myid_z == 0 PES
! This will cause trouble in the current version of MLP (Jim Taft)
! PILGRIM since all PEs meet at a barrier in all communication primitives 
! (since a group-wise barrier will take some time to implement)
!
      IF ( myid_2 == 0 ) THEN  
        CALL ParGather(  Comm_1, 0, Rsemiglobal3D, Y3D, Bglobal3D )
      ENDIF

      IF ( GID == 0 ) THEN
        IF ( SUM(Aglobal3D-Bglobal3D) /= 0.0_r8 ) THEN
          PRINT *, "ParUtilitiesTest failed: Scatter/Gather ver. != Orig."
          Passed = .FALSE.
        END IF
      ENDIF

      DEALLOCATE( Bglobal3D )
      DEALLOCATE( Aglobal3D )
      DEALLOCATE( Rsemiglobal3D )
      DEALLOCATE( Rlocal3D )

      CALL DecompFree( Z3D )
      CALL DecompFree( Y3D )
      CALL DecompFree( YZ3D )

      CALL ParFree( Comm_1 )
      CALL ParFree( Comm_2 )
#endif

!
! Test 2: A simple XY column distribution for 2D and 3D arrays
!
      ALLOCATE( XDist( npr_1 ) )
      ALLOCATE( YDist( npr_2 ) )
      ALLOCATE( ZDist( 1 ) )

!
! Consider a fairly irregular decomposition
!
      BlockLen = Im
      DO I = 1, npr_1-1
        XDist( I ) = BlockLen / 2
        BlockLen   = BlockLen - XDist( I ) 
      ENDDO
      XDist( npr_1 ) = BlockLen

      BlockLen = Jm
      DO J = 1, npr_2-1
        YDist( J ) = BlockLen / 2
        BlockLen   = BlockLen - YDist( J ) 
      ENDDO
      YDist( npr_2 ) = BlockLen

      Zdist(1) = Km

!
! Classical column distribution which could be used in the physics
!
      ALLOCATE( Rlocal2D( xdist(myid_1+1)*ydist(myid_2+1) ) )
      call decompcreate( npr_1, npr_2, xdist, ydist, XY2D )
      ALLOCATE( Rlocal3D( xdist(myid_1+1)*ydist(myid_2+1)*km ) )
      call decompcreate( npr_1, npr_2, 1, xdist,ydist,zdist, XY3D )


      DEALLOCATE( ZDist )
      DEALLOCATE( YDist )
      DEALLOCATE( XDist )

      ALLOCATE( Aglobal2D( im*jm ) )
      ALLOCATE( Aglobal3D( im*jm*km ) )
      ALLOCATE( Bglobal2D( im*jm ) )
      ALLOCATE( Bglobal3D( im*jm*km ) )

      CALL RANDOM_NUMBER( HARVEST = Aglobal2D )
      CALL RANDOM_NUMBER( HARVEST = Aglobal3D )

!
!     Now scatter the arrays over all PEs
!
      CALL ParScatter( CommGlobal, 0, Aglobal2D, XY2D, Rlocal2D )
      CALL ParGather( CommGlobal, 0, Rlocal2D, XY2D, Bglobal2D )
      CALL ParScatter( CommGlobal, 0, Aglobal3D, XY3D, Rlocal3D )
      CALL ParGather( CommGlobal, 0, Rlocal3D, XY3D, Bglobal3D )
      IF ( GID == 0 ) THEN
        IF ( SUM( Aglobal2d - Bglobal2d ) /= 0.0_r8 .OR.                     &
     &       SUM( Aglobal3d - Bglobal3d ) /= 0.0_r8 )  THEN
          PRINT *, "ParUtilitiesTest failed: Gather/scatter != Orig."
          Passed = .FALSE.
        END IF
      END IF

      DEALLOCATE( Rlocal3D )
      DEALLOCATE( Rlocal2D )
      DEALLOCATE( Bglobal3D )
      DEALLOCATE( Aglobal3D )
      DEALLOCATE( Bglobal2D )
      DEALLOCATE( Aglobal2D )

      CALL DecompFree( XY3D )
      CALL DecompFree( XY2D )

!
! Test 3 : Test ParExchangeVector by exchanging a vector twice
!
      ALLOCATE( Xdist( Gsize ) )
      ALLOCATE( Ydist( Gsize ) )

!
! Initialize seed to a different value on every PE, thus ensuring
! different values in the subsequent vectors
!
      CALL RANDOM_SEED( SIZE = BlockLen )
      ALLOCATE( Itmp( BlockLen ) )
      DO I=1, BlockLen
        Itmp(BlockLen) = Gid*(BlockLen-I+1)*Gid
      ENDDO
      CALL RANDOM_SEED( PUT = Itmp )
      DEALLOCATE( Itmp )


!
! Loop several times for better testing
!
      DO J = 1, 10
        ALLOCATE( Rtmp(Gsize) )
        CALL RANDOM_NUMBER( HARVEST = Rtmp )
!
! Determine a random destination pattern
!
        Xdist = INT( 10.0_r8 * Rtmp )
        DEALLOCATE( Rtmp )

        ALLOCATE( Rglobal( SUM(Xdist) ) )
        ALLOCATE( Rlocal( SUM(Xdist) ) )
        ALLOCATE( Rtmp( Gsize*SUM(Xdist) ) )   ! A maximum buffer size

        CALL RANDOM_NUMBER( HARVEST = Rglobal )
        CALL ParExchangeVector( CommGlobal,Xdist,Rglobal,Ydist,Rtmp )
        CALL ParExchangeVector( CommGlobal,Ydist,Rtmp,Xdist,Rlocal )
        IF ( SUM( Rglobal - Rlocal ) /= 0.0_r8 ) THEN
          PRINT *, "ParUtilitiesTest failed: ParExchangeVector"
          PRINT *, "Loop index ", J
          Passed = .FALSE.
          stop
        ENDIF

        DEALLOCATE( Rtmp )
        DEALLOCATE( Rlocal )
        DEALLOCATE( Rglobal )

      ENDDO

      DEALLOCATE( Ydist )
      DEALLOCATE( Xdist )

!
! Test 4: Parallel sums, ParBeginTransfer/ParEndTransfer
!
      ALLOCATE( Xdist( Gsize ) )
      ALLOCATE( Ydist( Gsize ) )

! Scalar sum test
      ALLOCATE( Rlocal(1) )    ! ParExchangeVector expects arrays
      ALLOCATE( Rglobal( Gsize ) )
      CALL RANDOM_NUMBER( HARVEST = Q )
      Rlocal(1) = Q
      Xdist = 0
      Xdist(1) = 1
      CALL ParExchangeVector( CommGlobal,Xdist,Rlocal,Ydist,Rglobal )
      CALL ParCollective( CommGlobal, SUMOP, Q )
      IF ( GID == 0 ) THEN
        Scalar = 0.0_r8
        DO pe=0,Gsize-1
          Scalar = Scalar + Rglobal( pe+1 )
        ENDDO
        IF ( ABS( Scalar - Q ) > Tolerance ) THEN
            print *, "Error in Scalar sum: ", Scalar-Q
        ENDIF
      ENDIF
      DEALLOCATE( Rlocal )
      DEALLOCATE( Rglobal )

! 1D Array sum test
      CALL RANDOM_NUMBER( HARVEST = Array1D )
      ALLOCATE( Rglobal( Im*Gsize ) )
      Xdist = 0
      Xdist(1) = Im
      CALL ParExchangeVector( CommGlobal,Xdist,Array1D,Ydist,Rglobal )
      CALL ParCollective( CommGlobal, SUMOP, Im, Array1D )
      IF ( GID == 0 ) THEN
        DO i=1, Im
          Scalar = 0.0_r8
          DO pe=0,Gsize-1
            Scalar = Scalar + Rglobal( Im*pe+I )
          ENDDO
          IF ( ABS( Scalar - Array1D(I) ) > Tolerance ) THEN
            print *, "Error in 1D Sum: ", Scalar-Array1D(I), " pos ",I
          ENDIF
        ENDDO
      ENDIF
      DEALLOCATE( Rglobal )

! 2D Array sum test
      CALL RANDOM_NUMBER( HARVEST = Array2D )

      ALLOCATE( Rglobal( Im*Jm*Gsize ) )
      ALLOCATE( Rlocal( Im*Jm ) )
      Xdist = 0
      Xdist(1) = Im*Jm

      inc = 0
      DO J=1,Jm
        DO I=1,Im
          inc = inc+1
          Rlocal(inc) = Array2D(i,j)
        ENDDO
      ENDDO
      CALL ParExchangeVector( CommGlobal,Xdist,Rlocal,Ydist,Rglobal )
      CALL ParCollective( CommGlobal, SUMOP, Im, Jm, Array2D )
      IF ( GID == 0 ) THEN
        inc = 0
        DO j=1, Jm
          DO i=1, Im
            inc = inc + 1
            Scalar = 0.0_r8
            DO pe=0,Gsize-1
              Scalar = Scalar + Rglobal( Jm*Im*pe+Inc )
            ENDDO
            IF ( ABS( Scalar - Array2D(I,J) ) > Tolerance ) THEN
              print *, "Error 2D Sum: ",Scalar-Array2D(I,J), "pos ",I,J
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      DEALLOCATE( Rglobal )
      DEALLOCATE( Rlocal )

! 3D Array sum test
      CALL RANDOM_NUMBER( HARVEST = Array3D )

      ALLOCATE( Rglobal( Im*Jm*Km*Gsize ) )
      ALLOCATE( Rlocal( Im*Jm*Km ) )
      Xdist = 0
      Xdist(1) = Im*Jm*Km

      inc = 0
      DO K=1,Km
        DO J=1,Jm
          DO I=1,Im
            inc = inc+1
            Rlocal(inc) = Array3D(i,j,k)
          ENDDO
        ENDDO
      ENDDO
      CALL ParExchangeVector( CommGlobal,Xdist,Rlocal,Ydist,Rglobal )
      CALL ParCollective( CommGlobal, SUMOP, Im, Jm, Km, Array3D )
      IF ( GID == 0 ) THEN
        inc = 0
        DO k=1, Km
          DO j=1, Jm
            DO i=1, Im
              inc = inc + 1
              Scalar = 0.0_r8
              DO pe=0,Gsize-1
                Scalar = Scalar + Rglobal( Km*Jm*Im*pe+inc )
              ENDDO
              IF ( ABS( Scalar - Array3D(I,J,K) ) > Tolerance ) THEN
              print *, "Error 3D Sum: ",Scalar-Array3D(I,J,K), I, J, K
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF


      DEALLOCATE( Rglobal )
      DEALLOCATE( Rlocal )

      DEALLOCATE( Ydist )
      DEALLOCATE( Xdist )

!
! That's all folks
!

      IF ( Gid == 0 ) THEN
        IF ( Passed ) THEN
          PRINT *, "Passed ParUtilitiesTest"
        ELSE
          PRINT *, "Failed ParUtilitiesTest"
        ENDIF
      END IF


      CALL ParExit( )

!EOC
!-------------------------------------------------------------------------
      END PROGRAM parutilitiestest
