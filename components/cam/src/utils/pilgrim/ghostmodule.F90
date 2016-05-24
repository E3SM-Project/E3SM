!------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS
!------------------------------------------------------------------------
      MODULE ghostmodule
!BOP
!
! !MODULE: ghostmodule
!
! !USES:
      USE decompmodule, ONLY : DecompType
#include "debug.h"
#include "pilgrim.h"
      IMPLICIT NONE

!
! !DESCRIPTION:
!
!      This module provides the basic support for "ghost regions".  In
!      reality the ghost region just subset of the global domain 
!      described by a decomposition (pro memoria: a decomposition 
!      describes a partition of a global index space over a number 
!      of PEs; this is inherently non-overlapping).  
!
!      It contains the following public types and routines.
!      \begin{center}
!      \begin{tabular}{|l|l|} \hline \hline
!         GhostType   & Type to describe ghosted local vector \\ \hline
!         GhostFree   & Destroy a ghost definition \\ \hline
!         GhostCreate & Copy ghost definition to newly created one\\ \hline 
!         GhostInfo   & Returns some information about the region \\ \hline 
!         \hline  \hline
!      \end{tabular}
!      \end{center}
!
!      GhostCreate is overloaded to support different types of domains:
!
!      \begin{center}
!      \begin{tabular}{|l|l|} \hline \hline
!         GhostCopy      & Copy a ghost region \\ \hline
!         GhostRegular1D & Define a subset of a 1D domain \\ \hline
!         GhostRegular2D & Define a subset of a 2D domain \\ \hline
!         GhostRegular3D & Define a subset of a 3D domain \\ \hline
!         GhostRegular4D & Define a subset of a 4D domain \\ \hline
!         GhostIrregular & Define a subset of an irregular domain \\ \hline
!         \hline  \hline
!      \end{tabular}
!      \end{center}
!
!      Generally one will use the GhostCreate routine which corresponds
!      to the underlying decomposition; e.g., if the decomposition was
!      defined with DecompRegular3D you would probably use GhostRegular3D
!      to define the ghost region.  But since decompositions and ghost
!      regions are generic, i.e., one-size-fits-all, this is not a requirement.
!      Be very careful if you use non-corresponding routines!
!
!      The ghost type contains a decomposition which describes the
!      {\it non-overlapping} distribution of the global domain 
!      (this is a replicated data structure with complete information
!      about all data structures on all PEs).  Its other components are
!      a list of the global indices of the on the boundary 
!      (not replicated), and a description of the mapping of the ghosted
!      local region to global indices.
!
!      This module is communication-free and is a foundation
!      for ParUtilitiesModule.  Since GhostType is local to the
!      PE, the modules routines can and should be called with 
!      non-replicated data structures.  Before boundary communication
!      takes place, the communication pattern derived from the ghost regions
!      must be created (see ParUtilitiesModule).  
!       
! !REVISION HISTORY:
!   00.11.10   Sawyer     Creation
!   01.02.07   Sawyer     Improvements; added Border to GhostType
!   01.02.12   Sawyer     Converted to free format
!   02.08.27   Zaslavsky  Changed intent from OUT to INOUT for objects of 
!                         GhostType
!   02.12.23   Sawyer     Added GhostRegular4D
!
! !PUBLIC TYPES:
      PUBLIC GhostType
      PUBLIC GhostFree
      PUBLIC GhostCreate
      PUBLIC GhostInfo

      INTERFACE     GhostCreate
        MODULE PROCEDURE GhostCopy
        MODULE PROCEDURE GhostIrregular
        MODULE PROCEDURE GhostRegular1D
        MODULE PROCEDURE GhostRegular2D
        MODULE PROCEDURE GhostRegular3D
        MODULE PROCEDURE GhostRegular4D
      END INTERFACE
 
! Decomposition info

       TYPE GhostType
         LOGICAL          :: Defined! Is it defined?
         TYPE(DecompType) :: Decomp ! Decomposition of global partition
         TYPE(DecompType) :: Local  ! Decomposition of local region
         TYPE(DecompType) :: Border ! Decomposition of local segment
       END TYPE GhostType

!EOP
      CONTAINS

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: GhostFree --- Free a ghosted region
!
! !INTERFACE:
      SUBROUTINE GhostFree ( Ghost )
! !USES:
      USE decompmodule, ONLY : DecompFree
      IMPLICIT NONE

! !INPUT/OUTPUT PARAMETERS:
      TYPE(GhostType), INTENT( INOUT ):: Ghost  ! Ghost information

!
! !DESCRIPTION:
!     Free the ghost decomposition -- deallocate the data structures.
!
! !SYSTEM ROUTINES:
!     ASSOCIATED, DEALLOCATE
!
! !REVISION HISTORY:
!   00.11.12   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
!
      CPP_ENTER_PROCEDURE( "GHOSTFREE" )

      IF ( Ghost%Defined ) THEN
         CALL DecompFree( Ghost%Border )
         CALL DecompFree( Ghost%Local  )
         CALL DecompFree( Ghost%Decomp ) 
         Ghost%Defined = .FALSE. 
      ENDIF

      CPP_LEAVE_PROCEDURE( "GHOSTFREE" )
      RETURN
!EOC
      END SUBROUTINE GhostFree
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: GhostDefined --- Is the ghost type de
!
! !INTERFACE:
      LOGICAL FUNCTION GhostDefined ( Ghost )
! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      TYPE(GhostType), INTENT( IN ):: Ghost  ! Ghost information

!
! !DESCRIPTION:
!     Returns true if Ghost has been created but not yet destroyed
!
! !REVISION HISTORY:
!   02.07.18   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
!
      CPP_ENTER_PROCEDURE( "GHOSTDEFINED" )
      GhostDefined = Ghost%Defined
      CPP_LEAVE_PROCEDURE( "GHOSTDEFINED" )

      RETURN
!EOC
      END FUNCTION GhostDefined
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: GhostCopy --- Copy one decomposition to another
!
! !INTERFACE:
      SUBROUTINE GhostCopy ( GhostIn, GhostOut )
! !USES:
      USE decompmodule, ONLY : DecompCopy
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      TYPE(GhostType), INTENT( IN )   :: GhostIn  ! Ghost information
!
! !OUTPUT PARAMETERS:
      TYPE(GhostType), INTENT( INOUT )  :: GhostOut ! Ghost information
!
! !DESCRIPTION:
!
!   Creates an output ghost definition and copies GhostIn to it 
!
! !SYSTEM ROUTINES:
!     ALLOCATE
!
! !REVISION HISTORY:
!   00.11.12   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER  :: I, Nsize

      CPP_ENTER_PROCEDURE( "GHOSTCOPY" )

      CALL DecompCopy( GhostIn%Decomp, GhostOut%Decomp )
      CALL DecompCopy( GhostIn%Local,  GhostOut%Local  )
      CALL DecompCopy( GhostIn%Border, GhostOut%Border )
      GhostOut%Defined = .TRUE.

      CPP_LEAVE_PROCEDURE( "GHOSTCOPY" )
      RETURN
!EOC
      END SUBROUTINE GhostCopy
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: GhostIrregular --- Create a ghost definition for 1-D grid
!
! !INTERFACE:
      SUBROUTINE GhostIrregular( Decomp, Id, LocalSize, Tags, Ghost )
! !USES:
      USE decompmodule, ONLY : DecompCreate, DecompCopy,                 &
                               DecompGlobalToLocal, DecompInfo
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      TYPE(DecompType), INTENT( IN ) :: Decomp     ! Decomp information
      INTEGER, INTENT( IN )          :: Id         ! Local PE identifer
      INTEGER, INTENT( IN )          :: LocalSize  ! Size of local segment
      INTEGER, INTENT( IN )          :: Tags(:)    ! Global tags
!
! !OUTPUT PARAMETERS:
      TYPE(GhostType), INTENT( INOUT ) :: Ghost  ! Ghost definition
!
!
! !DESCRIPTION:
!     Creates a ghost definition for a ghosted array given by
!     the PEs and Tags of the local points.  Note that none of the 
!     array bounds can be outside the global domain!
!
! !SYSTEM ROUTINES:
!     ALLOCATE, DEALLOCATE
!
! !REVISION HISTORY:
!   00.11.12   Sawyer     Creation
!
! !BUGS:
!   None of the array bounds can be outside of the global domain!
!   This is significant if the local region is on the edge of the
!   domain, and, in other words, the ghost region cannot cover
!   empty space.  This limitation may be relaxed in the future.
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER :: I, NPEs, GlobalSize, Local, Cnt, Ipe
      INTEGER, ALLOCATABLE :: Pe(:), Other(:)
!
!
      CPP_ENTER_PROCEDURE( "GHOSTIRREGULAR" )
!
! Allocate the basic data structures
!
      CALL DecompInfo( Decomp, Npes, GlobalSize )

      ALLOCATE( Pe( LocalSize ) )
      ALLOCATE( Other( LocalSize ) )

!
! Use decompmodule to create global and local portions of Ghost
! The local version is only on the local processor "0"
 
      Other = Id
      CALL DecompCreate( Npes, Other, LocalSize, Tags, Ghost%Local )

!
! Perform over all points local segment
!
      Cnt = 0
      DO I= 1, LocalSize
        CALL DecompGlobalToLocal( Decomp, Tags(I), Local, Ipe )
        CPP_ASSERT_F90( (Local .GT. 0) .AND. (ipe .GE. 0) )
        IF ( Ipe .ne. id ) THEN
          Cnt = Cnt + 1
          Other( Cnt ) = Tags(I)
          Pe( Cnt )    = Ipe
        ENDIF
      ENDDO

!
! Define the border regions.  Presumably Cnt << LocalSize
!

      CALL DecompCreate( Npes, Pe, Cnt, Other, Ghost%Border )

!
! Copy the decomposition too
!
      CALL DecompCopy( Decomp, Ghost%Decomp )

! Clean up

      DEALLOCATE( Pe )
      DEALLOCATE( Other )
      Ghost%Defined = .TRUE.

      CPP_LEAVE_PROCEDURE( "GHOSTIRREGULAR" )
      RETURN
!EOC
      END SUBROUTINE GhostIrregular
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: GhostRegular1D --- Create a ghost definition for 1-D grid
!
! !INTERFACE:
      SUBROUTINE GhostRegular1D( Decomp, Id, Xglobal, Xfrom, Xto, Xwrap, &
                                 Ghost )
! !USES:
      USE decompmodule, ONLY : DecompCreate, DecompCopy,                 &
                               DecompGlobalToLocal, DecompInfo
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      TYPE(DecompType), INTENT( IN )   :: Decomp ! Decomp information
      INTEGER, INTENT( IN )            :: Id     ! Local PE identifer
      INTEGER, INTENT( IN )            :: Xglobal! Total in X
      INTEGER, INTENT( IN )            :: Xfrom  ! Low index in X
      INTEGER, INTENT( IN )            :: Xto    ! High index in X
      LOGICAL, INTENT( IN )            :: Xwrap  ! Wrap in X?
!
! !OUTPUT PARAMETERS:p
      TYPE(GhostType), INTENT( INOUT )   :: Ghost  ! Ghost definition
!
!
! !DESCRIPTION:
!     Creates a ghost definition for a regular 1-D array with the
!     array bounds Xfrom:Xto.
!
!     If the array bounds are outside of the global domain they may
!     be wrapped around back into the global domain (variable Xwrap).  
!     If the region is not wrapped, it is advisable that the ghost 
!     region end at the boundary (which usually requires
!     special case treatment depending on the PE number). If 
!     it does not end at the boundary, undefined points are 
!     introduced.
!
! !SYSTEM ROUTINES:
!     ALLOCATE, DEALLOCATE
!
! !REVISION HISTORY:
!   00.11.12   Sawyer     Creation
!
! !BUGS:
!
!   There are certain limitations to ghost regions which can be
!   avoided by clean programming practices.  If the ghosted region
!   wraps back onto core regions of the same PE, problems can arise.  
!   The simple case -- a ghosted region on 1 PE -- is supported in
!   most cases.  However, if it wraps back onto the local PE 
!   in such a way that more than one ghost points is mapped to
!   one core domain global index, then the code may fail.  Note
!   that this is rarely the case if the ghost regions are small
!   and enough processors are used to avoid wrapping back on the
!   local one.
!   
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER :: I, L, NPEs, GlobalSize, LocalSize, Cnt, Local, Ipe
      INTEGER :: Global
      INTEGER, ALLOCATABLE :: Pe(:), Tags(:), Other(:)
!
!
      CPP_ENTER_PROCEDURE( "GHOSTREGULAR1D" )

!
! Allocate the basic data structures
!
      CALL DecompInfo( Decomp, NPEs, GlobalSize )
      CPP_ASSERT_F90( GlobalSize .EQ. Xglobal )

      LocalSize = Xto - Xfrom + 1
      CPP_ASSERT_F90( LocalSize .GE. 0 )

      ALLOCATE( Pe( LocalSize ) )
      ALLOCATE( Tags( LocalSize ) )
      ALLOCATE( Other( LocalSize ) )

!
! Perform over all points local segment
!
      Cnt = 0
      L = 0
      DO I = Xfrom, Xto
        L = L + 1
        Global = MODULO(I-1,Xglobal)+1  ! Wrap around condition
        IF (Xwrap .OR. Global==I) THEN
          Tags(L) = Global                 ! Global Tags
          CALL DecompGlobalToLocal( Decomp, Global, Local, Ipe )
          IF ( Ipe .ne. Id .AND. Ipe .GE. 0 ) THEN
            Cnt = Cnt + 1
            Other( Cnt ) = Global        ! Local Tags
            Pe( Cnt )    = Ipe
          ENDIF
!
! Special case: the domain wraps-around onto the same PE.  This is
! very tricky:  the ghost points are distinguished from their true
! local core domain counterparts by a minus sign.  This makes the
! address space in both Ghost%Border and Ghost%Local unique
!
          IF ( Ipe .eq. Id .AND. I .ne. Global ) THEN
            Cnt = Cnt + 1
            Other( Cnt ) = -Global       ! Local Tags
            Pe( Cnt )    = Ipe
            Tags(L) = -Global              ! Global Tags (mark ghost region!)
          ENDIF
        ELSE
          Tags(L) = 0
        ENDIF
      ENDDO

!
! Perform over all points local segment
!
      CALL DecompCreate( Npes, Pe, Cnt, Other, Ghost%Border )

!
! Use decompmodule to create global and local portions of Ghost
! The local version is only on the local PE
!
      Other = Id 
      CALL DecompCreate( Npes, Other, LocalSize, Tags, Ghost%Local )

!
! Copy the decomposition too
!
      CALL DecompCopy( Decomp, Ghost%Decomp )

! Clean up

      DEALLOCATE( Other    )
      DEALLOCATE( Tags     )
      DEALLOCATE( Pe       )

      Ghost%Defined = .TRUE.

      CPP_LEAVE_PROCEDURE( "GHOSTREGULAR1D" )
      RETURN
!EOC
      END SUBROUTINE GhostRegular1D
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: GhostRegular2D --- Create a ghost definition for 2-D grid
!
! !INTERFACE:
      SUBROUTINE GhostRegular2D( Decomp, Id, Xglobal, Xfrom, Xto, Xwrap, &
                                 Yglobal, Yfrom, Yto, Ywrap, Ghost )
! !USES:
      USE decompmodule, ONLY : DecompCreate, DecompCopy,                 &
                               DecompGlobalToLocal, DecompInfo
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      TYPE(DecompType), INTENT( IN )   :: Decomp ! Decomp information
      INTEGER, INTENT( IN )            :: Id     ! Local PE identifer
      INTEGER, INTENT( IN )            :: Xglobal! Total in X
      INTEGER, INTENT( IN )            :: Xfrom  ! Low index in X
      INTEGER, INTENT( IN )            :: Xto    ! High index in X
      LOGICAL, INTENT( IN )            :: Xwrap  ! Wrap in X?
      INTEGER, INTENT( IN )            :: Yglobal! Total in X
      INTEGER, INTENT( IN )            :: Yfrom  ! Distribution in X
      INTEGER, INTENT( IN )            :: Yto    ! Distribution in Y
      LOGICAL, INTENT( IN )            :: Ywrap  ! Wrap in Y?

!
! !OUTPUT PARAMETERS:
      TYPE(GhostType), INTENT( INOUT )   :: Ghost  ! Ghost definition
!
!
! !DESCRIPTION:
!     Creates a ghost definition for a regular 2-D array with the
!     array bounds Xfrom:Xto,Yfrom:Yto.
!
!     If the array bounds are outside of the global domain they may
!     be wrapped around back into the global domain (Xwrap, Ywrap).  
!     If the region is not wrapped, it is advisable that the ghost 
!     region end at the boundary (which usually requires
!     special case treatment depending on the PE number). If 
!     it does not end at the boundary, undefined points are 
!     introduced.
!
! !SYSTEM ROUTINES:
!     ALLOCATE, DEALLOCATE
!
! !REVISION HISTORY:
!   00.11.12   Sawyer     Creation
!
! !BUGS:
!
!   There are certain limitations to ghost regions which can be
!   avoided by clean programming practices.  If the ghosted region
!   wraps back onto core regions of the same PE, problems can arise.  
!   The simple case -- a ghosted region on 1 PE -- is supported in
!   most cases.  However, if it wraps back onto the local PE 
!   in such a way that more than one ghost points is mapped to
!   one core domain global index, then the code may fail.  Note
!   that this is rarely the case if the ghost regions are small
!   and enough processors are used to avoid wrapping back on the
!   local one.
!   
!   WARNING:  If the domain wraps around in both X and Y there is a 
!   the code should be run with at least 2 PEs so that in one of the
!   two dimensions there is no wrap-around onto the same PE.
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER :: I, J, L, Ipe, Npes, GlobalSize, LocalSize
      INTEGER :: Global, Cnt, Local, Xtrue, Ytrue
      INTEGER, ALLOCATABLE :: Pe(:), Tags(:), Other(:)
!
!
      CPP_ENTER_PROCEDURE( "GHOSTREGULAR2D" )

!
! Allocate the basic data structures
!
      CALL DecompInfo( Decomp, Npes, GlobalSize )
      CPP_ASSERT_F90( GlobalSize .EQ. Xglobal*Yglobal )

      LocalSize = (Xto - Xfrom + 1)*(Yto - Yfrom + 1)
      CPP_ASSERT_F90( LocalSize .GE. 0 )
      ALLOCATE( Pe( LocalSize ) )
      ALLOCATE( Tags( LocalSize ) )
      ALLOCATE( Other( LocalSize ) )
!
! Perform over all points local segment
!
      Cnt = 0
      L = 0
      DO J= Yfrom, Yto
        Ytrue = MODULO(J-1,Yglobal) + 1
        DO I= Xfrom, Xto
          Xtrue = MODULO(I-1,Xglobal) + 1
          L = L + 1
          Global = (Ytrue-1)*Xglobal + Xtrue
          IF ( (Xwrap.OR.(Xtrue==I)) .AND. (Ywrap.OR.(Ytrue==J)) ) THEN
            Tags( L ) = Global
            CALL DecompGlobalToLocal( Decomp, Global, Local, Ipe )
            IF ( Ipe .ne. Id .AND. Ipe .GE. 0 ) THEN
              Cnt = Cnt + 1
              Other( Cnt ) = Global     ! Local Tags
              Pe( Cnt )    = Ipe
            ENDIF
!
! Special case: the domain wraps-around onto the same PE.  This is
! very tricky:  the ghost points are distinguished from their true
! local core domain counterparts by a minus sign.  This makes the
! address space in both Ghost%Border and Ghost%Local unique
!
            IF ( Ipe.EQ.Id .AND. ( I.NE.Xtrue .OR. J.NE.Ytrue ) ) THEN
              Cnt = Cnt + 1
              Other( Cnt ) = -Global       ! Local Tags
              Pe( Cnt )    = Ipe
              Tags(L) = -Global              ! Global Tags (mark ghost region!)
            ENDIF
          ELSE
            Tags(L) = 0
          ENDIF
        ENDDO
      ENDDO

!
! Perform over all points local segment
!
      CALL DecompCreate( Npes, Pe, Cnt, Other, Ghost%Border )

!
! Use decompmodule to create global and local portions of Ghost
! The local version is only on the local PE
!
      Other = Id 
      CALL DecompCreate( Npes, Other, LocalSize, Tags, Ghost%Local )

!
! Copy the decomposition too
!
      CALL DecompCopy( Decomp, Ghost%Decomp )

! Clean up

      DEALLOCATE( Other )
      DEALLOCATE( Tags  )
      DEALLOCATE( Pe    )

      Ghost%Defined = .TRUE.

      CPP_LEAVE_PROCEDURE( "GHOSTREGULAR2D" )
      RETURN
!EOC
      END SUBROUTINE GhostRegular2D
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: GhostRegular3D --- Create a ghost definition for 3-D grid
!
! !INTERFACE:
      SUBROUTINE GhostRegular3D( Decomp, Id, Xglobal, Xfrom, Xto, Xwrap, &
                                 Yglobal, Yfrom, Yto, Ywrap,             &
                                 Zglobal, Zfrom, Zto, Zwrap, Ghost )
! !USES:
      USE decompmodule, ONLY : DecompCreate, DecompCopy,                 &
                               DecompGlobalToLocal, DecompInfo
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      TYPE(DecompType), INTENT( IN )   :: Decomp ! Decomp information
      INTEGER, INTENT( IN )            :: Id     ! Local PE identifer
      INTEGER, INTENT( IN )            :: Xglobal! Total in X
      INTEGER, INTENT( IN )            :: Xfrom  ! Low index in X
      INTEGER, INTENT( IN )            :: Xto    ! High index in X
      LOGICAL, INTENT( IN )            :: Xwrap  ! Wrap in X?
      INTEGER, INTENT( IN )            :: Yglobal! Total in Y
      INTEGER, INTENT( IN )            :: Yfrom  ! Distribution in Y
      INTEGER, INTENT( IN )            :: Yto    ! Distribution in Y
      LOGICAL, INTENT( IN )            :: Ywrap  ! Wrap in Y?
      INTEGER, INTENT( IN )            :: Zglobal! Total in Z
      INTEGER, INTENT( IN )            :: Zfrom  ! Distribution in Z
      INTEGER, INTENT( IN )            :: Zto    ! Distribution in Z
      LOGICAL, INTENT( IN )            :: Zwrap  ! Wrap in Z?
!
! !OUTPUT PARAMETERS:
      TYPE(GhostType), INTENT( INOUT )   :: Ghost  ! Ghost definition
!
!
! !DESCRIPTION:
!     Creates a ghost definition for a regular 3-D array with the
!     array bounds Xfrom:Xto,Yfrom:Yto,Zfrom:Zto. 
!
!     If the array bounds are outside of the global domain they may
!     be wrapped around back into the global domain (Xwrap, Ywrap).  
!     If the region is not wrapped, it is advisable that the ghost 
!     region end at the boundary (which usually requires
!     special case treatment depending on the PE number). If 
!     it does not end at the boundary, undefined points are 
!     introduced.
!
!
! !SYSTEM ROUTINES:
!     ALLOCATE, DEALLOCATE
!
! !REVISION HISTORY:
!   00.11.12   Sawyer     Creation
!
! !BUGS:
!   There are certain limitations to ghost regions which can be
!   avoided by clean programming practices.  If the ghosted region
!   wraps back onto core regions of the same PE, problems can arise.  
!   The simple case -- a ghosted region on 1 PE -- is supported in
!   most cases.  However, if it wraps back onto the local PE 
!   in such a way that more than one ghost points is mapped to
!   one core domain global index, then the code may fail.  Note
!   that this is rarely the case if the ghost regions are small
!   and enough processors are used to avoid wrapping back on the
!   local one.
!   
!   WARNING:  If the domain wraps around in two of the three dims 
!   the code should be run with at least 2 PEs so that in one of the
!   two dimensions there is no wrap-around onto the same PE.  If it
!   wraps around in all three dimensions it should be run on at least
!   4 PEs.  Note these are extremely rare toriodal cases.
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER :: I, J, K, L, Ipe, Npes, GlobalSize, LocalSize
      INTEGER :: Global, Cnt, Local, Xtrue, Ytrue, Ztrue
      LOGICAL :: IsX, IsY, IsZ
      INTEGER, ALLOCATABLE :: Pe(:), Tags(:), Other(:)
!
!
      CPP_ENTER_PROCEDURE( "GHOSTREGULAR3D" )

!
! Allocate the basic data structures
!
      CALL DecompInfo( Decomp, Npes, GlobalSize )
      CPP_ASSERT_F90( GlobalSize .EQ. Xglobal*Yglobal*Zglobal )

      LocalSize = (Xto-Xfrom+1) * (Yto-Yfrom+1) * (Zto-Zfrom+1)

      CPP_ASSERT_F90( LocalSize .GE. 0 )
      ALLOCATE( Pe( LocalSize ) )
      ALLOCATE( Tags( LocalSize ) )
      ALLOCATE( Other( LocalSize ) )
!
! Perform over all points local segment
!
      Cnt = 0
      L = 0
      DO K = Zfrom, Zto
        Ztrue = MODULO(K-1,Zglobal) + 1
        DO J = Yfrom, Yto
          Ytrue = MODULO(J-1,Yglobal) + 1
          DO I = Xfrom, Xto
            Xtrue = MODULO(I-1,Xglobal) + 1
            L = L + 1
            Global = ((Ztrue-1)*Yglobal+(Ytrue-1))*Xglobal+Xtrue
!
! Check to see if this is an defined global index
!
            CALL DecompGlobalToLocal( Decomp, Global, Local, Ipe )
            CPP_ASSERT_F90( (Local .GT. 0) .AND. (Ipe .GE. 0) )
!
! The wrapping case: mark as undefined

            IsX = Xtrue/=I
            IsY = Ytrue/=J
            IsZ = Ztrue/=K
            IF ( (.NOT.Xwrap.AND.IsX) .OR. (.NOT.Ywrap.AND.IsY)          &
                .OR. (.NOT.Zwrap.AND.IsZ) ) THEN
              Cnt = Cnt + 1
              Other( Cnt ) = 0           ! Local Tags
              Pe( Cnt )    = Ipe
              Tags( L )      = 0
            ELSE IF ( Ipe .ne. Id ) THEN
!
! Boundary case:  Global is in a ghost region not belonging
! to this PE.  Mark it in the border data structure (Arrays Other and Pe)
!
              Cnt = Cnt + 1
              Other( Cnt ) = Global      ! Local Tags
              Pe( Cnt )    = Ipe
              Tags( L )      = Global
            ELSE IF ( Ipe==Id .AND. (IsX.OR.IsY.OR.IsZ) ) THEN
!
! Special case: the domain wraps-around onto the same PE.  This is
! very tricky:  the ghost points are distinguished from their true
! local core domain counterparts by a minus sign.  This makes the
! address space in both Ghost%Border and Ghost%Local unique
!
              Cnt = Cnt + 1
              Other( Cnt ) = -Global     ! Local Tags
              Pe( Cnt )    = Ipe
              Tags(L)        = -Global     ! Global Tags (mark ghost region!)
            ELSE
              Tags( L ) = Global
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      CPP_ASSERT_F90( LocalSize==L )
!
! Perform over all points local segment
!
      CALL DecompCreate( Npes, Pe, Cnt, Other, Ghost%Border )

!
! Use decompmodule to create global and local portions of Ghost
! The local version is only on the local PE
!
      Other = Id
      CALL DecompCreate( Npes, Other, LocalSize, Tags, Ghost%Local )

!
! Copy the decomposition too
!
      CALL DecompCopy( Decomp, Ghost%Decomp )

! Clean up

      DEALLOCATE( Other )
      DEALLOCATE( Tags  )
      DEALLOCATE( Pe    )

      Ghost%Defined = .TRUE.

      CPP_LEAVE_PROCEDURE( "GHOSTREGULAR3D" )
      RETURN
!EOC
      END SUBROUTINE GhostRegular3D
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: GhostRegular4D --- Create a ghost definition for 4-D grid
!
! !INTERFACE:
      SUBROUTINE GhostRegular4D( Decomp, Id, Xglobal, Xfrom, Xto, Xwrap, &
                                 Yglobal, Yfrom, Yto, Ywrap,             &
                                 Zglobal, Zfrom, Zto, Zwrap,             &
                                 Tglobal, Tfrom, Tto, Twrap, Ghost )
! !USES:
      USE decompmodule, ONLY : DecompCreate, DecompCopy,                 &
                               DecompGlobalToLocal, DecompInfo
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      TYPE(DecompType), INTENT( IN )   :: Decomp ! Decomp information
      INTEGER, INTENT( IN )            :: Id     ! Local PE identifer
      INTEGER, INTENT( IN )            :: Xglobal! Total in X
      INTEGER, INTENT( IN )            :: Xfrom  ! Low index in X
      INTEGER, INTENT( IN )            :: Xto    ! High index in X
      LOGICAL, INTENT( IN )            :: Xwrap  ! Wrap in X?
      INTEGER, INTENT( IN )            :: Yglobal! Total in Y
      INTEGER, INTENT( IN )            :: Yfrom  ! Distribution in Y
      INTEGER, INTENT( IN )            :: Yto    ! Distribution in Y
      LOGICAL, INTENT( IN )            :: Ywrap  ! Wrap in Y?
      INTEGER, INTENT( IN )            :: Zglobal! Total in Z
      INTEGER, INTENT( IN )            :: Zfrom  ! Distribution in Z
      INTEGER, INTENT( IN )            :: Zto    ! Distribution in Z
      LOGICAL, INTENT( IN )            :: Zwrap  ! Wrap in Z?
      INTEGER, INTENT( IN )            :: Tglobal! Total in T
      INTEGER, INTENT( IN )            :: Tfrom  ! Distribution in T
      INTEGER, INTENT( IN )            :: Tto    ! Distribution in T
      LOGICAL, INTENT( IN )            :: Twrap  ! Wrap in T?
!
! !OUTPUT PARAMETERS:
      TYPE(GhostType), INTENT( INOUT )   :: Ghost  ! Ghost definition
!
!
! !DESCRIPTION:
!     Creates a ghost definition for a regular 3-D array with the
!     array bounds Xfrom:Xto,Yfrom:Yto,Zfrom:Zto,Tfrom:Tto.
!
!     If the array bounds are outside of the global domain they may
!     be wrapped around back into the global domain (Xwrap, Ywrap).  
!     If the region is not wrapped, it is advisable that the ghost 
!     region end at the boundary (which usually requires
!     special case treatment depending on the PE number). If 
!     it does not end at the boundary, undefined points are 
!     introduced.
!
! !SYSTEM ROUTINES:
!     ALLOCATE, DEALLOCATE
!
! !REVISION HISTORY:
!   02.12.23   Sawyer     Creation from GhostRegular3D
!
! !BUGS:
!   There are certain limitations to ghost regions which can be
!   avoided by clean programming practices.  If the ghosted region
!   wraps back onto core regions of the same PE, problems can arise.  
!   The simple case -- a ghosted region on 1 PE -- is supported in
!   most cases.  However, if it wraps back onto the local PE 
!   in such a way that more than one ghost points is mapped to
!   one core domain global index, then the code may fail.  Note
!   that this is rarely the case if the ghost regions are small
!   and enough processors are used to avoid wrapping back on the
!   local one.
!   
!   WARNING:  If the domain wraps around in two of the three dims 
!   the code should be run with at least 2 PEs so that in one of the
!   two dimensions there is no wrap-around onto the same PE.  If it
!   wraps around in all three dimensions it should be run on at least
!   4 PEs.  Note these are extremely rare toriodal cases.
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER :: I, J, K, L, M, Ipe, Npes, GlobalSize, LocalSize
      INTEGER :: Global, Cnt, Local, Xtrue, Ytrue, Ztrue, Ttrue
      LOGICAL :: IsX, IsY, IsZ, IsT
      INTEGER, ALLOCATABLE :: Pe(:), Tags(:), Other(:)
!
!
      CPP_ENTER_PROCEDURE( "GHOSTREGULAR4D" )

!
! Allocate the basic data structures
!
      CALL DecompInfo( Decomp, Npes, GlobalSize )
      CPP_ASSERT_F90( GlobalSize .EQ. Xglobal*Yglobal*Zglobal*Tglobal )

      LocalSize = (Xto-Xfrom+1)*(Yto-Yfrom+1)*(Zto-Zfrom+1)*(Tto-Tfrom+1)

      CPP_ASSERT_F90( LocalSize .GE. 0 )
      ALLOCATE( Pe( LocalSize ) )
      ALLOCATE( Tags( LocalSize ) )
      ALLOCATE( Other( LocalSize ) )
!
! Perform over all points local segment
!
      Cnt = 0
      M = 0
      DO L = Tfrom, Tto
        Ttrue = MODULO(L-1,Tglobal) + 1
        DO K = Zfrom, Zto
          Ztrue = MODULO(K-1,Zglobal) + 1
          DO J = Yfrom, Yto
            Ytrue = MODULO(J-1,Yglobal) + 1
            DO I = Xfrom, Xto
              Xtrue = MODULO(I-1,Xglobal) + 1
              M = M + 1
              Global = (((Ttrue-1)*Zglobal+(Ztrue-1))*Yglobal+(Ytrue-1)) &
                       *Xglobal+Xtrue
!
! Check to see if this is an defined global index
!
              CALL DecompGlobalToLocal( Decomp, Global, Local, Ipe )
              CPP_ASSERT_F90( (Local .GT. 0) .AND. (Ipe .GE. 0) )
!
! The wrapping case: mark as undefined

              IsX = Xtrue/=I
              IsY = Ytrue/=J
              IsZ = Ztrue/=K
              IsT = Ttrue/=L
              IF ( (.NOT.Xwrap.AND.IsX) .OR. (.NOT.Ywrap.AND.IsY)          &
                 .OR. (.NOT.Zwrap.AND.IsZ) .OR. (.NOT.Twrap.AND.IsT) ) THEN
                Cnt = Cnt + 1
                Other( Cnt ) = 0           ! Local Tags
                Pe( Cnt )    = Ipe
                Tags(M)      = 0
              ELSE IF ( Ipe .ne. Id ) THEN
!
! Boundary case:  Global is in a ghost region not belonging
! to this PE.  Mark it in the border data structure (Arrays Other and Pe)
!
                Cnt = Cnt + 1
                Other( Cnt ) = Global      ! Local Tags
                Pe( Cnt )    = Ipe
                Tags(M)      = Global
              ELSE IF ( Ipe==Id .AND. (IsX.OR.IsY.OR.IsZ.OR.IsT) ) THEN
!
! Special case: the domain wraps-around onto the same PE.  This is
! very tricky:  the ghost points are distinguished from their true
! local core domain counterparts by a minus sign.  This makes the
! address space in both Ghost%Border and Ghost%Local unique
!
                Cnt = Cnt + 1
                Other( Cnt ) = -Global     ! Local Tags
                Pe( Cnt )    = Ipe
                Tags(M)      = -Global     ! Global Tags (mark ghost region!)
              ELSE
                Tags(M) = Global
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CPP_ASSERT_F90( LocalSize==M )
!
! Perform over all points local segment
!
      CALL DecompCreate( Npes, Pe, Cnt, Other, Ghost%Border )

!
! Use decompmodule to create global and local portions of Ghost
! The local version is only on the local PE
!
      Other = Id
      CALL DecompCreate( Npes, Other, LocalSize, Tags, Ghost%Local )

!
! Copy the decomposition too
!
      CALL DecompCopy( Decomp, Ghost%Decomp )

! Clean up

      DEALLOCATE( Other )
      DEALLOCATE( Tags  )
      DEALLOCATE( Pe    )

      Ghost%Defined = .TRUE.

      CPP_LEAVE_PROCEDURE( "GHOSTREGULAR4D" )
      RETURN
!EOC
      END SUBROUTINE GhostRegular4D
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: GhostInfo --- Information about ghosted decompostion
!
! !INTERFACE:
      SUBROUTINE GhostInfo( Ghost, Npes,                                 &
                            GlobalSize, LocalSize, BorderSize )
! !USES:
      USE decompmodule, ONLY : DecompInfo
      IMPLICIT NONE

! !INPUT PARAMETERS:
      TYPE(GhostType), INTENT( IN ):: Ghost  ! Ghost information

! !INPUT PARAMETERS:
      INTEGER, INTENT( OUT )   :: Npes       ! Number of Pes
      INTEGER, INTENT( OUT )   :: GlobalSize ! Size of global domain
      INTEGER, INTENT( OUT )   :: LocalSize  ! Size of ghosted local region
      INTEGER, INTENT( OUT )   :: BorderSize ! Size of border
!
! !DESCRIPTION:
!     Return information about the ghosted region
!
! !SYSTEM ROUTINES:
!
! !REVISION HISTORY:
!   00.11.12   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
!
      CPP_ENTER_PROCEDURE( "GHOSTINFO" )

      CALL DecompInfo( Ghost%Decomp, Npes, GlobalSize )
      CALL DecompInfo( Ghost%Local,  Npes, LocalSize ) 
      CALL DecompInfo( Ghost%Border, Npes, BorderSize ) 

      CPP_LEAVE_PROCEDURE( "GHOSTINFO" )
      RETURN
!EOC
      END SUBROUTINE GhostInfo
!-----------------------------------------------------------------------

      END MODULE ghostmodule
