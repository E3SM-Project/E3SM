#include "pilgrim.h"
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS
!-------------------------------------------------------------------------
      MODULE redistributemodule
#if defined( SPMD )
!BOP
!
! !MODULE: redistributemodule
!
! !USES:
#include "debug.h"
#if !defined(STAND_ALONE)
      use shr_kind_mod, only: r8 => shr_kind_r8
#endif
      IMPLICIT NONE

!
! !DESCRIPTION:
!
!
! !REVISION HISTORY:
!   99.01.18   Sawyer     Creation
!   99.11.17   Sawyer     Added RedistributeStart, RedistributeFinish
!   00.07.20   Sawyer     Minor cosmetic changes
!   00.08.28   Sawyer     Accommodated change to ParEndTranfer interface
!   01.02.12   Sawyer     Converted to free format
!
! !PUBLIC TYPES:
      PUBLIC RedistributeType
      PUBLIC RedistributeCreate, RedistributeFree, RedistributePerform
      PUBLIC RedistributeStart, RedistributeFinish
 
! Redistribution info

      TYPE RedistributeType
        INTEGER, POINTER     :: CountA(:) ! Per PE counts in Decomp A
        INTEGER, POINTER     :: CountB(:) ! Per PE counts in Decomp B
        INTEGER, POINTER     :: PermA(:)  ! Permutation in Decomp A
        INTEGER, POINTER     :: PermB(:)  ! Permutation in Decomp B
      END TYPE RedistributeType
!EOP
      REAL(CPP_REAL8), ALLOCATABLE, SAVE :: InStatic(:), OutStatic(:)

      CONTAINS

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: RedistributeCreate --- Create an inter-decomp. structure
!
! !INTERFACE:
      SUBROUTINE RedistributeCreate( DecompA, DecompB, Inter )
! !USES:
      USE decompmodule, ONLY:  DecompType, DecompGlobalToLocal
      USE parutilitiesmodule, ONLY:  GID, Gsize
      IMPLICIT NONE
#include "mpif.h"

!
! !INPUT PARAMETERS:
      TYPE(DecompType), INTENT( IN ) :: DecompA       ! Decomposition A
      TYPE(DecompType), INTENT( IN ) :: DecompB       ! Decomposition B

! !OUTPUT PARAMETERS:
      TYPE(RedistributeType), INTENT( OUT ) :: Inter  ! Inter info.

!
! !DESCRIPTION:
!
!     This routine constructs a RedistributeType structure which
!     can be efficiently used in the RedistributePerform routine.
!
! !SYSTEM ROUTINES:
!     ALLOCATE
!
! !REVISION HISTORY:
!   99.01.15   Sawyer     Creation
!
! !BUGS:
!    Currently untested.
!   
!EOP
!---------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER IndexA, IndexB, I, J, K, Tag, Local, Pe, Offsets( Gsize )
      LOGICAL Found, First, Search

      CPP_ENTER_PROCEDURE( "REDISTRIBUTECREATE" )
!
! Allocate the number of entries and list head arrays
!
      CPP_ASSERT_F90( SIZE( DecompA%NumEntries ).EQ. Gsize )

      CPP_ASSERT_F90( SIZE( DecompB%NumEntries ).EQ. Gsize )

      ALLOCATE( Inter%CountA( Gsize ) )
      ALLOCATE( Inter%CountB( Gsize ) )

      ALLOCATE( Inter%PermA( DecompA%NumEntries( GID+1 ) ) )
      ALLOCATE( Inter%PermB( DecompB%NumEntries( GID+1 ) ) )

      Inter%CountA = 0
      Inter%CountB = 0
      IndexA       = 0
      IndexB       = 0

      DO I = 1, Gsize
        DO J = 1, SIZE( DecompB%Head(I)%StartTags )
          First = .TRUE.
          DO Tag=DecompB%Head(I)%StartTags(J),DecompB%Head(I)%EndTags(J)

!
! CODE INLINED FOR PERFORMANCE
!
!!!            CALL DecompGlobalToLocal( DecompA, Tag, Local, Pe )

            Search = .TRUE.
            IF ( .NOT. First )                                           &
              Search = First .OR. Tag .GT. DecompA%Head(Pe+1)%EndTags(K)

            IF ( Search ) THEN
              First = .FALSE.

!
! Search over all the PEs
!
              Pe = -1
              Found = .FALSE.
              DO WHILE ( .NOT. Found )
!
! Copy the number of entries on each PE
!
                Pe = Pe + 1
!
! Search through the local data segment
!
                Local = 1
                K = 1
                DO WHILE ( .NOT. Found .AND.                             &
                    K .LE. SIZE( DecompA%Head(Pe+1)%StartTags ) )
                  IF ( Tag .GE. DecompA%Head(Pe+1)%StartTags(K) .AND.    &
                       Tag .LE. DecompA%Head(Pe+1)%EndTags(K) ) THEN
                    Local = Local+Tag - DecompA%Head(Pe+1)%StartTags(K)
                    Found = .TRUE.
                  ELSE
                    Local = Local + DecompA%Head(Pe+1)%EndTags(K) -      &
                            DecompA%Head(Pe+1)%StartTags(K) + 1
                    K = K+1
                  ENDIF
                ENDDO
!
! Emergency brake
!
                IF ( Pe.EQ.(SIZE(DecompA%Head)-1).AND. .NOT.Found ) THEN
                  Found = .TRUE.
                  Local = 0
                  Pe    = -1
                ENDIF
              ENDDO
          
!
! END OF INLINING
!
            ELSE
              Local = Local + 1
            ENDIF
!
! Calculate the sorting permutation for A
!
            IF ( Pe .EQ. GID ) THEN
              Inter%CountA( I ) = Inter%CountA( I ) + 1
              IndexA = IndexA + 1
              Inter%PermA( IndexA ) = local
            ENDIF
!
! Calculate the sorting permutation for B
!
            IF ( I-1 .EQ. GID ) THEN
              Inter%CountB( Pe+1 ) = Inter%CountB( Pe+1 ) + 1
              IndexB = IndexB + 1
              Inter%PermB( IndexB ) = Inter%CountB( Pe+1 )*Gsize + Pe
            ENDIF
          ENDDO
        ENDDO
      ENDDO

!
! Finally decode PermB and add in the proper offsets
!
      Offsets = 0
      DO I=1, Gsize-1
        Offsets( I+1 ) = Offsets( I ) + Inter%CountB( I )
      ENDDO
      DO I=1, IndexB
        Pe = MOD( Inter%PermB( I ), Gsize )
        Inter%PermB( I ) = Inter%PermB(I)/Gsize + Offsets( Pe+1 )
      ENDDO
        
      CPP_LEAVE_PROCEDURE( "REDISTRIBUTECREATE" )
      RETURN
!EOC
      END SUBROUTINE RedistributeCreate
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: RedistributePerform --- Perform the Redistribution
!
! !INTERFACE:
      SUBROUTINE RedistributePerform( Inter, Forward, Input, Output )
! !USES:
      USE parutilitiesmodule, ONLY : CommGlobal, Gsize,                  &
                                     ParExchangeVector,GID
      IMPLICIT NONE

!
! !INPUT PARAMETERS:
      TYPE(RedistributeType), INTENT( INOUT ) :: Inter ! Inter info.
      LOGICAL             :: Forward     ! True: A -> B  False: B -> A
      REAL(CPP_REAL8), INTENT( IN )  :: Input( * )  ! Input Array
! !OUTPUT PARAMETERS:
      REAL(CPP_REAL8), INTENT( OUT ) :: Output( * ) ! Output Array
!
! !DESCRIPTION:
!
!     This routine performs the redistribution of Input to Output
!     according to the RedistributionType data structure Inter.
!     The redistribution can be from A -> B ("forward") or B -> A
!     ("backward").  This feature has been added to avoid the
!     need of a separate Inter (which requires considerable
!     memory) to perform the backward redistribution.
!
! !SYSTEM ROUTINES:
!     ALLOCATE, DEALLOCATE
!
! !BUGS:
!     Currently limited to the global communicator.
!
! !REVISION HISTORY:
!   99.01.15   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC

! !LOCAL VARIABLES:
      INTEGER I, Ierr, LenOutBuf( Gsize )
      REAL(CPP_REAL8), ALLOCATABLE :: InBuf(:), OutBuf(:)

      CPP_ENTER_PROCEDURE( "REDISTRIBUTEPERFORM" )

      IF ( Forward ) THEN
!
! Forward redistribution
!
        ALLOCATE( InBuf( SUM( Inter%CountA ) ) )
        ALLOCATE( OutBuf( SUM( Inter%CountB ) ) )
        DO I = 1, SUM( Inter%CountA )
          InBuf( I ) = Input( Inter%PermA( I ) )
        ENDDO

        CALL ParExchangeVector( CommGlobal, Inter%CountA, InBuf,         &
                                LenOutBuf, OutBuf )
        DO I = 1, SUM( Inter%CountB )
          Output( I ) = OutBuf( Inter%PermB( I ) )
        ENDDO
        DEALLOCATE( OutBuf )
        DEALLOCATE( InBuf )

      ELSE
!
! Backward redistribution
!
        ALLOCATE( InBuf( SUM( Inter%CountB ) ) )
        ALLOCATE( OutBuf( SUM( Inter%CountA ) ) )
        DO I = 1, SUM( Inter%CountB )
          InBuf( Inter%PermB( I ) ) = Input( I )
        ENDDO

        CALL ParExchangeVector( CommGlobal, Inter%CountB, InBuf,         &
                                LenOutBuf, OutBuf )
        DO I = 1, SUM( Inter%CountA )
          Output( Inter%PermA( I ) ) = OutBuf( I )
        ENDDO
        DEALLOCATE( OutBuf )
        DEALLOCATE( InBuf )
     
      ENDIF
      CPP_LEAVE_PROCEDURE( "REDISTRIBUTEPERFORM" )
      RETURN
!EOC
      END SUBROUTINE RedistributePerform
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: RedistributeFree --- Free an inter-decomp. structure
!
! !INTERFACE:
      SUBROUTINE RedistributeFree( Inter )
! !USES:
      IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
      TYPE(RedistributeType), INTENT( INOUT ) :: Inter ! Inter info.
!
! !DESCRIPTION:
!
!     This routine frees a RedistributeType structure.
!
! !SYSTEM ROUTINES:
!     DEALLOCATE
!
! !REVISION HISTORY:
!   99.01.15   Sawyer     Creation
!
! !BUGS:
!    Currently untested.
!   
!EOP
!-----------------------------------------------------------------------
!BOC

! !LOCAL VARIABLES:
      INTEGER Ierr

      CPP_ENTER_PROCEDURE( "REDISTRIBUTEFREE" )

      DEALLOCATE( Inter%PermB  )
      DEALLOCATE( Inter%PermA  )
      DEALLOCATE( Inter%CountB )
      DEALLOCATE( Inter%CountA )

      CPP_LEAVE_PROCEDURE( "REDISTRIBUTEFREE" )
      RETURN
!EOC
      END SUBROUTINE RedistributeFree
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: RedistributeStart --- Perform Asynchronous Redistribution
!
! !INTERFACE:
      SUBROUTINE RedistributeStart( Inter, Forward, Input )
! !USES:
      USE parutilitiesmodule, ONLY : CommGlobal, Gsize,                  &
                                     ParBeginTransfer,GID
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      TYPE(RedistributeType), INTENT( INOUT ) :: Inter ! Inter info.
      LOGICAL             :: Forward     ! True: A -> B  False: B -> A
      REAL(CPP_REAL8), INTENT( IN )  :: Input( * )  ! Input Array
!
! !DESCRIPTION:
!
!     This routine starts an asynchronous redistribution of Input 
!     to Output according to the RedistributionType data structure Inter.
!     The redistribution can be from A -> B ("forward") or B -> A
!     ("backward").  This feature has been added to avoid the
!     need of a separate Inter (which requires considerable
!     memory) to perform the backward redistribution.
!
!     Beware: both RedistributeStart and RedistributeFinish *must*
!     be called with the same values of Inter and Forward.  Nesting
!     of asynchronous distributions is forbidden.  In addition, any
!     other communication in the between RedistributeStart and
!     RedistributeFinish cannot used the communicator "CommGlobal" 
!     provided by parutilitiesmodule.
!
! !SYSTEM ROUTINES:
!     ALLOCATE
!
! !REVISION HISTORY:
!   99.11.17   Sawyer     Creation from RedistributePerform
!
! !BUGS:
!     Currently limited to the global communicator.
!   
!EOP
!-----------------------------------------------------------------------
!BOC

! !LOCAL VARIABLES:
      INTEGER I, Ierr, Dest( Gsize ), Src( Gsize )

      CPP_ENTER_PROCEDURE( "REDISTRIBUTESTART" )

      DO I = 1, Gsize
        Dest( I ) = I-1
        Src( I )  = I-1
      ENDDO

      IF ( Forward ) THEN
!
! Forward redistribution
!
        ALLOCATE( InStatic( SUM( Inter%CountA ) ) )
        ALLOCATE( OutStatic( SUM( Inter%CountB ) ) )
        DO I = 1, SUM( Inter%CountA )
          InStatic( I ) = Input( Inter%PermA( I ) )
        ENDDO
        CALL ParBeginTransfer( CommGlobal, Gsize, Gsize, Dest, Src,      &
                               InStatic, Inter%CountA,                   &
                               OutStatic, Inter%CountB )

      ELSE
!
! Backward redistribution
!
        ALLOCATE( InStatic( SUM( Inter%CountB ) ) )
        ALLOCATE( OutStatic( SUM( Inter%CountA ) ) )
        DO I = 1, SUM( Inter%CountB )
          InStatic( Inter%PermB( I ) ) = Input( I )
        ENDDO
        CALL ParBeginTransfer( CommGlobal, Gsize, Gsize, Dest, Src,      &
                               InStatic, Inter%CountB,                   &
                               OutStatic, Inter%CountA )

      ENDIF
      CPP_LEAVE_PROCEDURE( "REDISTRIBUTESTART" )
      RETURN
!EOC
      END SUBROUTINE RedistributeStart
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: RedistributeFinish --- Complete Asynchronous Redistribution
!
! !INTERFACE:
      SUBROUTINE RedistributeFinish( Inter, Forward, Output )
! !USES:
      USE parutilitiesmodule, ONLY: CommGlobal,Gsize,ParEndTransfer,GID
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      TYPE(RedistributeType), INTENT( INOUT ) :: Inter ! Inter info.
      LOGICAL             :: Forward     ! True: A -> B  False: B -> A
! !OUTPUT PARAMETERS:
      REAL(CPP_REAL8), INTENT( OUT ) :: Output( * ) ! Output Array
!
! !DESCRIPTION:
!
!     This routine completes an asynchronous redistribution of Input 
!     to Output according to the RedistributionType data structure Inter.
!     The redistribution can be from A -> B ("forward") or B -> A
!     ("backward").  This feature has been added to avoid the
!     need of a separate Inter (which requires considerable
!     memory) to perform the backward redistribution.
!
!     See additional documentation in RedistributeStart.
!
! !SYSTEM ROUTINES:
!     DEALLOCATE
!
! !REVISION HISTORY:
!   99.11.17   Sawyer     Creation from RedistributePerform
!
! !BUGS:
!     Currently limited to the global communicator.
!   
!EOP
!-----------------------------------------------------------------------
!BOC

! !LOCAL VARIABLES:
      INTEGER I, Dest( Gsize ), Src( Gsize )

      CPP_ENTER_PROCEDURE( "REDISTRIBUTEFINISH" )

      DO I = 1, Gsize
        Dest( I ) = I-1
        Src( I )  = I-1
      ENDDO

      IF ( Forward ) THEN
        CALL ParEndTransfer( CommGlobal, Gsize, Gsize, Dest, Src,        &
                             InStatic, Inter%CountA,                     &
                             OutStatic, Inter%CountB )
        DO I = 1, SUM( Inter%CountB )
          Output( I ) = OutStatic( Inter%PermB( I ) )
        ENDDO
      ELSE
        CALL ParEndTransfer( CommGlobal, Gsize, Gsize, Dest, Src,        &
                             InStatic, Inter%CountB,                     &
                             OutStatic, Inter%CountA )
        DO I = 1, SUM( Inter%CountA )
          Output( Inter%PermA( I ) ) = OutStatic( I )
        ENDDO
      ENDIF
      DEALLOCATE( OutStatic )
      DEALLOCATE( InStatic )

      CPP_LEAVE_PROCEDURE( "REDISTRIBUTEFINISH" )
      RETURN
!EOC
      END SUBROUTINE RedistributeFinish
!-----------------------------------------------------------------------
#endif
      END MODULE redistributemodule
