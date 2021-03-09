#if !defined(STAND_ALONE)
#endif
#define _SMEMORY 1
!-----------------------------------------------------------------------
!         Nasa/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS
!-----------------------------------------------------------------------
      MODULE parutilitiesmodule
#if defined( SPMD )
!BOP
!
! !MODULE: parutilitiesmodule
!
! !USES:
#if defined( STAND_ALONE )
# define iulog 6
#else
      use cam_logfile, only: iulog
#endif
#if !defined(STAND_ALONE)
      USE shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8, &
                              r4 => shr_kind_r4
#endif
      USE mod_comm, ONLY : commglobal, gid, numpro, blockdescriptor, max_nparcels
#include "debug.h"
      IMPLICIT NONE
#include "mpif.h"
#include "pilgrim.h"

!
! !PUBLIC DATA MEMBERS:
      PUBLIC     Gsize
      PUBLIC     INT4, REAL4, REAL8
      PUBLIC     SUMOP, MAXOP, MINOP, BCSTOP
       

      INTEGER,SAVE :: GSize        ! Size of communicator CommGlobal
                                   ! Equivalent to mod_comm::numpro
#define CPP_SUM_OP 101
#define CPP_MAX_OP 102
#define CPP_MIN_OP 103
#define CPP_BCST_OP 104

      INTEGER,SAVE :: INT4  = MPI_INTEGER
      INTEGER,SAVE :: REAL4 = MPI_REAL
      INTEGER,SAVE :: REAL8 = MPI_DOUBLE_PRECISION
      INTEGER,SAVE :: SUMOP = MPI_SUM
      INTEGER,SAVE :: MAXOP = MPI_MAX
      INTEGER,SAVE :: MINOP = MPI_MIN
      INTEGER,SAVE :: BCSTOP = CPP_BCST_OP

! !PUBLIC MEMBER FUNCTIONS:
      PUBLIC ParPatternType

      TYPE ParPatternType
        INTEGER ::     Comm                  ! Communicator
        INTEGER ::     Iam                   ! My rank in communicator
        INTEGER ::     Size                  ! Size of communicator
        TYPE(BlockDescriptor), POINTER :: SendDesc(:) ! Array of descriptors
        TYPE(BlockDescriptor), POINTER :: RecvDesc(:) ! Array of descriptors
      END TYPE ParPatternType 


#ifdef _SMEMORY
      TYPE ParInfoType
        INTEGER :: numRecvSeg               ! number of received segments
        INTEGER :: numSendSeg               ! number of send segments
        INTEGER :: maxNumSeg                ! maximum number of segments over all processors
        INTEGER :: numRecvNeigh             ! number of receive neighbors
        INTEGER :: numSendNeigh             ! number of send neighbors
      END TYPE ParInfoType
#endif

      PUBLIC     ParInit, ParSplit, ParFree, ParExit
      PUBLIC     ParScatter, ParGather
      PUBLIC     ParBeginTransfer, ParEndTransfer
      PUBLIC     ParExchangeVector, ParCollective
      PUBLIC     ParPatternCreate, ParPatternFree

      INTERFACE     ParPatternCreate
        MODULE PROCEDURE ParPatternCopy
        MODULE PROCEDURE ParPatternGhost
        MODULE PROCEDURE ParPatternDecompToDecomp
        MODULE PROCEDURE ParPatternDecompToGhost
        MODULE PROCEDURE ParPatternGhostToDecomp
        MODULE PROCEDURE ParPatternGhostToGhost
      END INTERFACE
 
      INTERFACE     ParScatter
        MODULE PROCEDURE ParScatterReal
        MODULE PROCEDURE ParScatterReal4
        MODULE PROCEDURE ParScatterInt
      END INTERFACE
 
      INTERFACE     ParGather
        MODULE PROCEDURE ParGatherReal
        MODULE PROCEDURE ParGatherReal4
        MODULE PROCEDURE ParGatherInt
      END INTERFACE

      INTERFACE     ParBeginTransfer
        MODULE PROCEDURE ParBeginTransferReal
        MODULE PROCEDURE ParBeginTransferPattern1D
        MODULE PROCEDURE ParBeginTransferPattern1Dint
        MODULE PROCEDURE ParBeginTransferPattern2D
        MODULE PROCEDURE ParBeginTransferPattern3D
        MODULE PROCEDURE ParBeginTransferPattern4D
!        MODULE PROCEDURE ParBeginTransferInt
      END INTERFACE

      INTERFACE     ParEndTransfer
        MODULE PROCEDURE ParEndTransferReal
        MODULE PROCEDURE ParEndTransferPattern1D
        MODULE PROCEDURE ParEndTransferPattern1Dint
        MODULE PROCEDURE ParEndTransferPattern2D
        MODULE PROCEDURE ParEndTransferPattern3D
        MODULE PROCEDURE ParEndTransferPattern4D
!        MODULE PROCEDURE ParEndTransferInt
      END INTERFACE

      INTERFACE     ParExchangeVector
        MODULE PROCEDURE ParExchangeVectorReal
        MODULE PROCEDURE ParExchangeVectorReal4
        MODULE PROCEDURE ParExchangeVectorInt
      END INTERFACE

      INTERFACE     ParCollective
        MODULE PROCEDURE ParCollectiveBarrier
        MODULE PROCEDURE ParCollective0D
        MODULE PROCEDURE ParCollective1D
        MODULE PROCEDURE ParCollective1DReal4
        MODULE PROCEDURE ParCollective2D
        MODULE PROCEDURE ParCollective2DReal4
        MODULE PROCEDURE ParCollective3D
        MODULE PROCEDURE ParCollective0DInt
        MODULE PROCEDURE ParCollective0DStr
        MODULE PROCEDURE ParCollective1DInt
        MODULE PROCEDURE ParCollective1DStr
        MODULE PROCEDURE ParCollective2DInt
      END INTERFACE

#ifdef _SMEMORY
      INTERFACE   ParCalcInfo
        MODULE PROCEDURE ParCalcInfoDecompToGhost
        MODULE PROCEDURE ParCalcInfoDecompToDecomp
        MODULE PROCEDURE ParCalcInfoGhostToGhost
        MODULE PROCEDURE ParCalcInfoGhostToDecomp
      END INTERFACE
#endif

!
! !DESCRIPTION:
!
!      This module provides the basic utilities to support parallelism
!      on a distributed or shared memory multiprocessor.
!
!      \begin{center}
!      \begin{tabular}{|l|l|} \hline \hline
!        ParInit           & Initialize the parallel system \\ \hline
!        ParExit           & Exit from the parallel system \\ \hline
!        ParSplit          & Create a Compute grid of PEs   \\ \hline
!        ParFree           & Free a split communicator \\ \hline
!        ParScatter        & Scatter global slice to local slices \\ \hline
!        ParGather         & Gather local slices to one global \\ \hline
!        ParBeginTransfer  & Initiate an all-to-all packet transfer \\ \hline
!        ParEndTransfer    & Complete an all-to-all packet transfer \\ \hline
!        ParExchangeVector & Complete an all-to-all packet transfer \\ \hline
!        ParCollective     & Collective operation across communicator \\ \hline
!      \end{tabular}
!      \end{center}
!      \vspace{2mm}
!
!      Other utilities can be added to this module as needs evolve.
!
!      Conceptually the intention is to aggregate as many of the
!      MPI communication calls as possible into a well-maintained
!      module.  This will help avoid the occurrence of MPI spaghetti 
!      code.  
!
!      This module is tailored to GEOS DAS and implements the 
!      design of Lucchesi/Mirin/Sawyer/Larson.
!
! !REVISION HISTORY:
!   97.02.01   Sawyer     Creation
!   97.07.22   Sawyer     Removal of DecompType related subroutines
!   97.08.13   Sawyer     Added ParScatter/Gather for Integers
!   97.09.26   Sawyer     Additions of Sparse communication primitives
!   97.12.01   Sawyer     Changed all MPI_SSEND to MPI_ISEND
!   97.12.23   Lucchesi   Added member variables IsIONode and InterComm
!   98.01.06   Sawyer     Additions from RL for I/O Nodes
!   98.02.02   Sawyer     Added the Cartesian data members
!   98.02.05   Sawyer     Removed the use of intercommunicators
!   98.02.23   Sawyer     Added ghosting utilities
!   98.02.25   Sawyer     Modified interface of BeginTransfer
!   98.03.03   Sawyer     Added Global ID number to public data members
!   98.03.25   Sawyer     Added documentation for walkthrough
!   98.04.16   Sawyer     Removed all use of MPI_CART (CommRow redefined)
!   98.07.23   Sawyer     Added ParGhost, ParPoleDot; ParBegin/EndGhost out
!   98.09.15   Sawyer     Added ParMerge, ParPoleGhost
!   98.09.17   Sawyer     Added ParSum, removed ParPoleDot
!   99.01.18   Sawyer     Minor cleaning
!   99.03.04   Sawyer     Revised SHMEM concept for Transfer
!   99.04.22   Sawyer     Removed COMMON for handles -- they are
!                         always used in same program unit.
!   99.05.21   Sawyer     Reintroduced barriers in Scatter/Gather
!   99.06.03   Sawyer     USE_SHMEM revisions
!   99.12.10   Sawyer     ParInit now sets GID, Gsize
!   99.12.13   Sawyer     Version slimmed down for FVCCM release
!   00.06.14   Sawyer     Precision module now used
!   00.07.07   Sawyer     Removed 2D scatter/gather; simplified API
!   00.07.30   Sawyer     Full implementation with shared memory
!   00.08.09   Sawyer     Replaced ParSum with ParCollective
!   00.08.28   Sawyer     Moved LLNL 2D data to LLNL2DModule; new MLP impl
!   01.02.04   Sawyer     Added PatternType and related routines
!   01.02.12   Sawyer     Converted to free format
!   02.10.30   Sawyer     Welded with mod_comm
!   03.03.06   Sawyer     Fix parpatterncreate for MPI2; use MPI_DATATYPE_NULL
!   05.10.12   Worley     Support for vectorization modifications in mod_comm
!   06.03.01   Sawyer     Merged CAM and GEOS5 versions
!   07.01.05   Mirin      Eliminated direct use of Gsize
!   07.09.04   Dennis     Reduced temporary memory usage
!
! !BUGS:
!   There are several MPI_Barriers at locations in the code.
!   These avoid potential race conditions which probably only occur
!   if the number of real processors is less than the number of
!   message passing processes.  Remove these barriers at your own risk
!
!EOP

      INTEGER, SAVE :: InHandle(MAX_PAX, MAX_SMP, MAX_TRF)
      INTEGER, SAVE :: OutHandle(MAX_PAX,MAX_SMP, MAX_TRF)
      INTEGER, SAVE :: BegTrf = 0  ! Ongoing overlapped begintransfer # 
      INTEGER, SAVE :: EndTrf = 0  ! Ongoing overlapped endtransfer #
      LOGICAL, SAVE :: Initialized = .FALSE. ! Flag for initialization of parutilitiesmodule.

      CONTAINS
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParInit --- Initialize the parallel execution
!
! !INTERFACE: 
      SUBROUTINE ParInit ( Comm, npryzxy, mod_method, mod_geopk, mod_gatscat, mod_maxirr )
!
! !USES:
      USE mod_comm, ONLY : mp_init
      IMPLICIT NONE
! !INPUT PARAMETERS:
      INTEGER, OPTIONAL  :: Comm
      INTEGER, OPTIONAL, INTENT(IN) :: npryzxy(4)      ! 2D decompositions
      INTEGER, OPTIONAL, INTENT(IN) :: mod_method      ! CAM optimization
      INTEGER, OPTIONAL, INTENT(IN) :: mod_geopk       ! CAM optimization
      INTEGER, OPTIONAL, INTENT(IN) :: mod_gatscat     ! CAM optimization
      INTEGER, OPTIONAL, INTENT(IN) :: mod_maxirr      ! CAM max simul. trsps.

!
! !DESCRIPTION:
!     Initializes the system.  In MPI mode, call MPI\_INIT if not done 
!     already. If the optional arguments are not provided, default 
!     values will be chosen.  But it is advisable to provide COMM
!     (main communicator) and NPRYZXY (internal 2D decomposition).
!
! !SYSTEM ROUTINES:
!     MPI_INITIALIZED, MPI_INIT
!
! !REVISION HISTORY:
!   97.03.20   Sawyer     Creation
!   97.04.16   Sawyer     Cleaned up for walk-through
!   97.07.03   Sawyer     Reformulated documentation
!   00.07.23   Sawyer     Added shared memory arena implementation
!   02.10.30   Sawyer     Now uses mp_init from mod_comm
!   06.06.15   Sawyer     Added CAM optimizations (passed to mod_comm)
!
!EOP
!-----------------------------------------------------------------------
!BOC

! Initialize mod_comm

      IF (.NOT. Initialized) THEN
         CALL mp_init( Comm, npryzxy, mod_method, mod_geopk, mod_gatscat, mod_maxirr )
         Gsize = numpro   !   Support PILGRIM's Gsize for now
         Initialized = .TRUE.
      ENDIF

      RETURN
!EOC
      END SUBROUTINE ParInit
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParExit --- Finalize the parallel execution
!
! !INTERFACE:
      SUBROUTINE ParExit ( Comm )

! !USES:
      USE mod_comm, ONLY: mp_exit
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, OPTIONAL  :: Comm

! !DESCRIPTION:
!     All PEs, compute nodes and IO nodes alike meet here to terminate
!     themselves.  If someone does not check in, everything will hang
!     here.
!
!     This routine is the very {\em last} thing which is executed!
!
! !LOCAL VARIABLES:
      INTEGER Ierror
!
! !SYSTEM ROUTINES:
!     MPI_BARRIER, MPI_FINALIZE
!
! !REVISION HISTORY:
!   97.03.20   Sawyer     Creation
!   97.04.16   Sawyer     Cleaned up for walk-through
!   97.07.03   Sawyer     Reformulated documentation
!   00.07.23   Sawyer     Added shared memory arena implementation
!   02.08.13   Sawyer     Incorporated mod_comm for low level comm.
!
!EOP
!-----------------------------------------------------------------------
!BOC
      CALL mp_exit(Comm)
      RETURN
!EOC
      END SUBROUTINE ParExit
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE:   ParSplit --- Split into group for I/O and computation
!
! !INTERFACE:
      SUBROUTINE ParSplit( InComm, Color, InID, Comm, MyID, Nprocs )
!
! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )     :: InComm    ! Communicator to split
      INTEGER, INTENT( IN )     :: Color     ! Group label
      INTEGER, INTENT( IN )     :: InID      ! Input ID

! !OUTPUT PARAMETERS:
      INTEGER, INTENT( OUT )    :: Comm      ! Split communicator
      INTEGER, INTENT( OUT )    :: MyID      ! Group label
      INTEGER, INTENT( OUT )    :: Nprocs    ! Number of PEs in my group
!
! !DESCRIPTION:
!     This routine splits the PEs into groups.  This is currently only
!     supported in MPI mode. Read the chapter on MPI\_COMM\_SPLIT 
!     thoroughly.  
!
! !SYSTEM ROUTINES:
!     MPI_COMM_SPLIT, MPI_COMM_SIZE, MPI_COMM_RANK
!
! !REVISION HISTORY:
!   97.03.20   Sawyer     Creation
!   97.04.16   Sawyer     Cleaned up for walk-through
!   97.07.03   Sawyer     Reformulated documentation
!   97.12.01   Sawyer     Xnodes and Ynodes are explicit arguments
!   97.12.23   Lucchesi   Added call to MPI_INTERCOMM_CREATE
!   98.01.06   Sawyer     Additions from RL for I/O Nodes
!   98.02.02   Sawyer     Added the Cartesian information
!   98.02.05   Sawyer     Removed the use of intercommunicators
!   98.04.16   Sawyer     Removed all use of MPI_CART (CommRow redefined)
!   99.01.10   Sawyer     CommRow now defined for all rows
!   00.07.09   Sawyer     Removed 2D computational mesh
!   00.08.08   Sawyer     Redefined as wrapper to mpi_comm_split
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER  Ierror

      CPP_ENTER_PROCEDURE( "PARSPLIT" )
!
!     Split the communicators
!
      CALL MPI_COMM_SPLIT( InComm, Color, InID, Comm, Ierror )
      IF ( Comm .ne. MPI_COMM_NULL ) THEN
        CALL MPI_COMM_RANK( Comm, MyID, Ierror )
        CALL MPI_COMM_SIZE( Comm, Nprocs, Ierror )
      ELSE
!
!     This PE does not participate: mark with impossible values
!
        MyID = -1
        Nprocs = -1
      ENDIF

      CPP_LEAVE_PROCEDURE( "PARSPLIT" )
      RETURN
!EOC
      END SUBROUTINE ParSplit
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE:   ParFree --- Free a communicator
!
! !INTERFACE:
      SUBROUTINE ParFree( InComm ) 
!
! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER InComm

!
! !DESCRIPTION:
!     This routine frees a communicator created with ParSplit
!
! !REVISION HISTORY:
!   97.09.11   Sawyer     Creation, to complement ParSplit
!   00.07.24   Sawyer     Revamped ParMerge into a free communicator 
!
! !LOCAL VARIABLES:
      INTEGER  Ierror
!
!EOP
!-----------------------------------------------------------------------
!BOC
      CPP_ENTER_PROCEDURE( "PARFREE" )
!
      CALL MPI_COMM_FREE( InComm, Ierror ) 
      CPP_LEAVE_PROCEDURE( "PARFREE" )
      RETURN
!EOC
      END SUBROUTINE ParFree
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE:   ParPatternCopy --- Duplicate/replicate a comm pattern
!
! !INTERFACE:
      SUBROUTINE ParPatternCopy( InComm, PatternIn, PatternOut, Multiplicity )
!
! !USES:
      USE mod_comm, ONLY : get_partneroffset
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER,  INTENT( IN )               :: InComm  ! # of PEs
      TYPE(ParPatternType), INTENT( IN )   :: PatternIn   ! Comm Pattern
      INTEGER, INTENT( IN ),  OPTIONAL     :: Multiplicity

! !OUTPUT PARAMETERS:
      TYPE(ParPatternType), INTENT( OUT )  :: PatternOut  ! Comm Pattern
!
! !DESCRIPTION:
!     This routine duplicates a given communication pattern. 
!
!     Optionally a multiplicity can be added.  This replicates the
!     communication pattern Mult times, that is for the case that
!     the data structures are replicated in the final dimension
!     Mult times.  A typical example is a pattern describing a 2D
!     array, e.g. a a lat-lon decomposition, which will be used
!     to copy a 3D lat-lon-lev array.  The strides (e.g. the number
!     of elements in one plane) of the source (send) and target 
!     (recv) arrays are now calculated internally.
!
! !SYSTEM ROUTINES:
!     MPI_TYPE_UB, MPI_TYPE_HVECTOR, MPI_TYPE_COMMIT
!
! !REVISION HISTORY:
!   03.03.20   Sawyer     Creation
!   03.06.26   Sawyer     Removed StrideSend/Recv from API
!
!EOP 
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER  Stride_S, Stride_R, Mult, Iam, GroupSize, Ipe, Ierror
      INTEGER  Disp, Length, I, J, ub, method

      CPP_ENTER_PROCEDURE( "PARPATTERNCOPY" )

      method = PatternIn%RecvDesc(1)%method

!
! Decide if this is a simple copy, or a multiple replication
!
      IF ( present(Multiplicity) ) THEN
          Mult = Multiplicity
      ELSE
          Mult = 1
      ENDIF

      CALL MPI_COMM_DUP( PatternIn%Comm, PatternOut%Comm, Ierror )
      CALL MPI_COMM_SIZE( PatternIn%Comm, GroupSize, Ierror )
      CALL MPI_COMM_RANK( PatternIn%Comm, Iam, Ierror )

      PatternOut%Iam  = Iam
      PatternOut%Size = GroupSize

      ALLOCATE( PatternOut%SendDesc( GroupSize ) )
      ALLOCATE( PatternOut%RecvDesc( GroupSize ) )

      PatternOut%SendDesc(:)%method = PatternIn%SendDesc(:)%method
      PatternOut%RecvDesc(:)%method = PatternIn%RecvDesc(:)%method
!
! Determine the strides which are by construction the maximum upper
! bound of all the derived types.  This is due to the fact that
! there are no 'holes' in the data types: even if one PE does not
! send to any other PEs, it will still have a data type for 'sending'
! data to itself.
!
        Stride_S = 0
        Stride_R = 0       
        DO Ipe=1, GroupSize
          IF ( PatternIn%SendDesc(Ipe)%type /= MPI_DATATYPE_NULL ) THEN
            CALL MPI_TYPE_UB( PatternIn%SendDesc(Ipe)%type, ub, ierror )
            Stride_S = max(Stride_S,ub)
          ENDIF
          IF ( PatternIn%RecvDesc(Ipe)%type /= MPI_DATATYPE_NULL ) THEN
            CALL MPI_TYPE_UB( PatternIn%RecvDesc(Ipe)%type, ub, ierror )
            Stride_R = max(Stride_R,ub)
          ENDIF
        ENDDO

!
! Determine the output data types
!
        DO Ipe=1, GroupSize
          IF ( PatternIn%SendDesc(ipe)%type /= MPI_DATATYPE_NULL ) THEN
            CALL MPI_TYPE_HVECTOR( Mult, 1, Stride_S, PatternIn%SendDesc(Ipe)%type,&
                                   PatternOut%SendDesc(Ipe)%type, Ierror )
            CALL MPI_TYPE_COMMIT( PatternOut%SendDesc(Ipe)%type, Ierror )
          ELSE
            PatternOut%SendDesc(ipe)%type = MPI_DATATYPE_NULL
          ENDIF
          IF ( PatternIn%RecvDesc(Ipe)%type /= MPI_DATATYPE_NULL ) THEN
            CALL MPI_TYPE_HVECTOR( Mult, 1, Stride_R, PatternIn%RecvDesc(Ipe)%type,&
                                   PatternOut%RecvDesc(Ipe)%type, Ierror )
            CALL MPI_TYPE_COMMIT( PatternOut%RecvDesc(Ipe)%type, Ierror )
          ELSE
            PatternOut%RecvDesc(ipe)%type = MPI_DATATYPE_NULL
          ENDIF
        ENDDO

!
! Determine the stride, which is the sum of all the blocksizes for all
! the derived types (there are no 'holes').
!      
        Stride_S = 0
        Stride_R = 0       
        DO Ipe=1, GroupSize
          Stride_S = Stride_S + sum( PatternIn%SendDesc(ipe)%BlockSizes(:) )
          Stride_R = Stride_R + sum( PatternIn%RecvDesc(ipe)%BlockSizes(:) )
        ENDDO

        DO ipe=1, GroupSize
          Length = SIZE(PatternIn%SendDesc(ipe)%BlockSizes) 
          ALLOCATE( PatternOut%SendDesc(ipe)%Displacements(Length*Mult) )
          ALLOCATE( PatternOut%SendDesc(ipe)%BlockSizes(Length*Mult) )
#if defined( DEBUG_PARPATTERNCOPY )
          write(iulog,*) "Multiplicity", Mult
          write(iulog,*) "Old send blocksizes", PatternIn%SendDesc(ipe)%BlockSizes
#endif
          DO i=1, Length
            Disp = PatternIn%SendDesc(ipe)%Displacements(i)
            DO j=1, Mult
              PatternOut%SendDesc(ipe)%BlockSizes(i+(j-1)*Length) =     &
                    PatternIn%SendDesc(ipe)%BlockSizes(i)
              PatternOut%SendDesc(ipe)%Displacements(i+(j-1)*Length) = Disp
              Disp = Disp + Stride_S
            ENDDO
          ENDDO
          PatternOut%SendDesc(ipe)%Nparcels  = &
            size (PatternOut%SendDesc(ipe)%Displacements)
          PatternOut%SendDesc(ipe)%Tot_Size = &
            sum  (PatternOut%SendDesc(ipe)%Blocksizes)
          Max_Nparcels = max (Max_Nparcels, PatternOut%SendDesc(ipe)%Nparcels)
#if defined( DEBUG_PARPATTERNCOPY )
          write(iulog,*) "Send blocksizes", PatternOut%SendDesc(ipe)%BlockSizes
          write(iulog,*) "Old recv blocksizes", PatternIn%RecvDesc(ipe)%BlockSizes
#endif
          Length = SIZE(PatternIn%RecvDesc(ipe)%BlockSizes) 
          ALLOCATE( PatternOut%RecvDesc(ipe)%Displacements(Length*Mult) )
          ALLOCATE( PatternOut%RecvDesc(ipe)%BlockSizes(Length*Mult) )
          DO i=1, Length
            Disp = PatternIn%RecvDesc(ipe)%Displacements(i)
            DO j=1, Mult
              PatternOut%RecvDesc(ipe)%BlockSizes(i+(j-1)*Length) =     &
                    PatternIn%RecvDesc(ipe)%BlockSizes(i)
              PatternOut%RecvDesc(ipe)%Displacements(i+(j-1)*Length) = Disp
              Disp = Disp + Stride_R
            ENDDO
          ENDDO
          PatternOut%RecvDesc(ipe)%Nparcels  = &
            size (PatternOut%RecvDesc(ipe)%Displacements)
          PatternOut%RecvDesc(ipe)%Tot_Size = &
            sum  (PatternOut%RecvDesc(ipe)%Blocksizes)
          Max_Nparcels = max (Max_Nparcels, PatternOut%RecvDesc(ipe)%Nparcels)
#if defined( DEBUG_PARPATTERNCOPY )
          write(iulog,*) "Recv blocksizes", PatternOut%RecvDesc(ipe)%BlockSizes
#endif
        ENDDO

        CALL get_partneroffset( InComm, PatternOut%SendDesc, PatternOut%RecvDesc )
      
      CPP_LEAVE_PROCEDURE( "PARPATTERNCOPY" )
      RETURN
!EOC
      END SUBROUTINE ParPatternCopy
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE:   ParPatternGhost --- Create pattern for given ghosting
!
! !INTERFACE:
      SUBROUTINE ParPatternGhost( InComm, Ghost, Pattern, mod_method, T )
!
! !USES:
      USE decompmodule, ONLY : DecompGlobalToLocal, DecompLocalToGlobal
      USE ghostmodule, ONLY : GhostType, GhostInfo
      USE mod_comm, ONLY : get_partneroffset
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER,  INTENT( IN )               :: InComm  ! # of PEs
      TYPE(GhostType),  INTENT( IN )       :: Ghost   ! # of PEs
      INTEGER,  INTENT( IN ), OPTIONAL     :: mod_method ! contiguous or derived type
      INTEGER, INTENT( IN ),  OPTIONAL     :: T       ! 

! !OUTPUT PARAMETERS:
      TYPE(ParPatternType), INTENT( OUT )  :: Pattern ! Comm Pattern
!
! !DESCRIPTION:
!     This routine contructs a communication pattern from the ghost
!     region definition.  That is, the resulting communication pattern
!     can be used in ParBegin/EndTransfer with the ghosted arrays as
!     inputs.  
!
! !SYSTEM ROUTINES:
!    MPI_COMM_SIZE,  MPI_COMM_RANK, MPI_COMM_DUP
!    MPI_TYPE_INDEXED, MPI_TYPE_COMMIT (depending on method)
!
! !REVISION HISTORY:
!   01.02.10   Sawyer     Creation
!   01.06.02   Sawyer     Renamed ParPatternGhost
!   02.06.27   Sawyer     Added data type "T" as optional argument
!   03.03.04   Sawyer     Set partneroffsets field
!   03.11.11   Mirin      Added optional argument mod_method
!
!EOP 
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER  i, j, ipe, pe, Iam, GroupSize, Num, Length, Ptr, Ierror
      INTEGER  Global, End, Local, GlobalSize, LocalSize, BorderSize
      INTEGER  DataType
      INTEGER, ALLOCATABLE :: InVector(:), OutVector(:)
      INTEGER, ALLOCATABLE :: LenInVector(:), LenOutVector(:)
      INTEGER              :: method

      CPP_ENTER_PROCEDURE( "PARPATTERNGHOST" )

      IF (present(T)) THEN
        DataType = T
      ELSE
        DataType = CPP_MPI_REAL8
      ENDIF

      IF (present(mod_method)) THEN
        method = mod_method
      ELSE
        method = 0     ! Default method - see mod_comm for description
      ENDIF
!
! First request the needed ghost values from other processors.
!
      CALL MPI_COMM_DUP( InComm, Pattern%Comm, Ierror )
      CALL MPI_COMM_SIZE( InComm, GroupSize, Ierror )
      CALL MPI_COMM_RANK( InComm, Iam, Ierror )

      Pattern%Iam  = Iam
      Pattern%Size = GroupSize

      ALLOCATE( Pattern%SendDesc( GroupSize ) )
      ALLOCATE( Pattern%RecvDesc( GroupSize ) )

      Pattern%SendDesc(:)%method = method
      Pattern%RecvDesc(:)%method = method

!
! Temporary variables
!
      ALLOCATE( LenInVector( GroupSize ) )
      ALLOCATE( LenOutVector( GroupSize ) )

      CALL GhostInfo( Ghost,GroupSize,GlobalSize,LocalSize,BorderSize )
      ALLOCATE( InVector( 2*BorderSize ) )
      ALLOCATE( OutVector( 2*LocalSize ) )

!
! A rather complicated loop to define the local ghost region.
! The concept is the following:  go through all the points in the
! border data structure.   It contains global indices of the points
! which have to be copied over from neighboring PEs.  These indices
! are collected into InVector for transmission to those PEs, in
! effect informing them of the local PEs requirements.
!
! A special case is supported:  if the ghost domain wraps around
! onto the domain of the local PE!  This is very tricky, because
! the index space in both Ghost%Border and Ghost%Local MUST be
! unique for DecompGlobalToLocal to work.   Solution:  ghost 
! points are marked with the negative value of the needed domain 
! value in both Ghost%Border and Ghost%Local.  These are "snapped 
! over" to the true global index with the ABS function, so that 
! they can be subsequently found in the true local domain.
!
      j = 1
      DO ipe=1, GroupSize
        Num = SIZE(Ghost%Border%Head(ipe)%StartTags)
        Length = 0
        DO i = 1, Num
          Global = Ghost%Border%Head(ipe)%StartTags(i)
          IF ( Global /= 0 ) THEN
            Length = Length + 1
            End    = Ghost%Border%Head(ipe)%EndTags(i)
            InVector(j) = ABS(Global)
            InVector(j+1) = ABS(End)
            CALL DecompGlobalToLocal( Ghost%Local, Global, Local, Pe )
            OutVector(Length) = Local-1                ! Zero-based address
            OutVector(Length+Num) = End - Global+1     ! Parcel size
            j = j + 2
          ENDIF
        ENDDO
        LenInVector(ipe) = 2*Length

!
! Set the receive buffer descriptor
!
#if defined(DEBUG_PARPATTERNGHOST)
        write(iulog,*) "Iam",Iam,"Pe",Ipe-1,"Lens",OutVector(Num+1:Num+Length), &
             "Displacements", OutVector(1:Length)
#endif

          IF ( Length > 0 .and. method > 0 ) THEN
            CALL MPI_TYPE_INDEXED( Length, OutVector(Num+1), OutVector,    &
                                   DataType, Ptr, Ierror )
            CALL MPI_TYPE_COMMIT( Ptr, Ierror )
            Pattern%RecvDesc(ipe)%type = Ptr
          ELSE
            Pattern%RecvDesc(ipe)%type = MPI_DATATYPE_NULL
          ENDIF

          ALLOCATE( Pattern%RecvDesc(ipe)%Displacements(Length) )
          ALLOCATE( Pattern%RecvDesc(ipe)%BlockSizes(Length) )
          DO i=1, Length
            Pattern%RecvDesc(ipe)%Displacements(i) = OutVector(i)
            Pattern%RecvDesc(ipe)%BlockSizes(i)    = OutVector(Num+i)
          ENDDO            
          Pattern%RecvDesc(ipe)%Nparcels  = &
            size (Pattern%RecvDesc(ipe)%Displacements)
          Pattern%RecvDesc(ipe)%Tot_Size = &
            sum  (Pattern%RecvDesc(ipe)%Blocksizes)
          Max_Nparcels = max (Max_Nparcels, Pattern%RecvDesc(ipe)%Nparcels)

      ENDDO

!
! Everybody exchanges the needed information
!
#if defined(DEBUG_PARPATTERNGHOST)
      write(iulog,*) "iam", iam, "In", LenInVector,                            &
                InVector( 1:SUM(LenInVector) )
#endif
      CALL ParExchangeVectorInt( InComm, LenInVector, InVector,          &
                                     LenOutVector, OutVector )
#if defined(DEBUG_PARPATTERNGHOST)
      write(iulog,*) "iam", iam, "Out", LenOutVector,                          &
                OutVector( 1:SUM(LenOutVector) )
#endif

!
! Now everyone has the segments which need to be sent to the 
! immediate neighbors.  Save these in PatternType.
!
      j = 1
      DO ipe = 1, GroupSize
        Num = LenOutVector(ipe) / 2
        DO i = 1, Num
          CALL DecompGlobalToLocal( Ghost%Local,OutVector(j),Local,pe )
          InVector(i) = Local-1
          InVector(i+Num) = OutVector(j+1) - OutVector(j) + 1
          j = j + 2
        ENDDO
#if defined(DEBUG_PARPATTERNGHOST)
        write(iulog,*) "Iam", Iam, "To", ipe-1, "InVector",                    &
              InVector(1:Num), "block size", InVector(Num+1:2*Num)
#endif

          IF ( Num > 0 .and. method > 0 ) THEN
            CALL MPI_TYPE_INDEXED( Num, InVector(Num+1), InVector,         &
                                   DataType, Ptr, Ierror )
            CALL MPI_TYPE_COMMIT( Ptr, Ierror )
            Pattern%SendDesc(ipe)%type = Ptr
          ELSE
            Pattern%SendDesc(ipe)%type = MPI_DATATYPE_NULL
          ENDIF

          ALLOCATE( Pattern%SendDesc(ipe)%Displacements(Num) )
          ALLOCATE( Pattern%SendDesc(ipe)%BlockSizes(Num) )
          DO i=1, Num
            Pattern%SendDesc(ipe)%Displacements(i) = InVector(i)
            Pattern%SendDesc(ipe)%BlockSizes(i)    = InVector(Num+i)
          ENDDO            
          Pattern%SendDesc(ipe)%Nparcels  = &
            size (Pattern%SendDesc(ipe)%Displacements)
          Pattern%SendDesc(ipe)%Tot_Size = &
            sum  (Pattern%SendDesc(ipe)%Blocksizes)
          Max_Nparcels = max (Max_Nparcels, Pattern%SendDesc(ipe)%Nparcels)

      ENDDO

      CALL get_partneroffset( InComm, Pattern%SendDesc, Pattern%RecvDesc )

!
! Clean up the locally allocate variables
!
      DEALLOCATE( OutVector )
      DEALLOCATE( InVector )
      DEALLOCATE( LenOutVector )
      DEALLOCATE( LenInVector )

      CPP_LEAVE_PROCEDURE( "PARPATTERNGHOST" )
      RETURN
!EOC
      END SUBROUTINE ParPatternGhost
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE:   ParPatternDecompToDecomp --- Create pattern between decomps
!
! !INTERFACE:
      SUBROUTINE ParPatternDecompToDecomp( InComm, DA, DB, Pattern, mod_method, T )
!
! !USES:
      USE decompmodule, ONLY : DecompType, DecompGlobalToLocal, DecompInfo
      USE mod_comm, ONLY : get_partneroffset
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER,  INTENT( IN )               :: InComm  ! # of PEs
      TYPE(DecompType),  INTENT( IN )      :: DA      ! Source Decomp Desc
      TYPE(DecompType),  INTENT( IN )      :: DB      ! Target Decomp Desc
      INTEGER,  INTENT( IN ), OPTIONAL     :: mod_method ! contiguous or derived type
      INTEGER, INTENT( IN ),  OPTIONAL     :: T       ! 

! !OUTPUT PARAMETERS:
      TYPE(ParPatternType), INTENT( OUT )  :: Pattern ! Comm Pattern
!
! !DESCRIPTION:
!     This routine contructs a communication pattern for a 
!     transformation from one decomposition to another, i.e., a 
!     so-called "transpose". The resulting communication pattern 
!     can be used in ParBegin/EndTransfer with the decomposed 
!     arrays as inputs.  
!
! !SYSTEM ROUTINES:
!    MPI_COMM_SIZE,  MPI_COMM_RANK, MPI_COMM_DUP
!    MPI_TYPE_INDEXED, MPI_TYPE_COMMIT (depending on method)
!
! !REVISION HISTORY:
!   01.05.29   Sawyer     Creation from RedistributeCreate
!   01.07.13   Sawyer     Rewritten to minimize DecompGlobalToLocal
!   02.07.16   Sawyer     Added data type T
!   03.11.11   Mirin      Added optional argument mod_method
!   07.03.11   Mirin      Generalized to different sized decompositions
!   07.09.04   Dennis     Reduced amount of temporary memory usage
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER I, J, Tag, Local, Pe, LenB, JB, Ipe, Num, Inc, Off
      INTEGER Ptr                                ! Pointer type
      INTEGER GroupSize, Iam, Ierror, DataType
      INTEGER OldPe, TotalPtsA, NpesA, TotalPtsB, NpesB
      INTEGER              :: method
      INTEGER              :: nCount,maxCount,ierr,sz
      INTEGER              :: lenBjmd,nNeigh,maxLenB,maxNeigh
#ifdef _SMEMORY
      TYPE (ParInfoType) :: Info
#endif

      INTEGER, ALLOCATABLE :: Count(:)           ! # segments for each recv PE
      INTEGER, ALLOCATABLE :: CountOut(:)        ! # segments for each send PE

      INTEGER, ALLOCATABLE :: DisplacementsA(:)  ! Generic displacements
      INTEGER, ALLOCATABLE :: BlockSizesA(:)     ! Generic block sizes
      INTEGER, ALLOCATABLE :: LocalA(:)          ! Generic Local indices

      INTEGER, ALLOCATABLE :: DisplacementsB(:)  ! Displacements for B
      INTEGER, ALLOCATABLE :: BlockSizesB(:)     ! Block sizes for B
      INTEGER, ALLOCATABLE :: LocalB(:)          ! Local indices for B
      INTEGER, ALLOCATABLE :: PeB(:)             ! Processor element numbers

      CPP_ENTER_PROCEDURE( "PARPATTERNDECOMPTODECOMP" )

      IF (present(T)) THEN
        DataType = T
      ELSE
        DataType = CPP_MPI_REAL8
      ENDIF

      IF (present(mod_method)) THEN
        method = mod_method
      ELSE
        method = 0     ! Default method - see mod_comm for description
      ENDIF

! Assume this routine is called by processes [ 0,max(NpesA,NpesB) )

      CALL DecompInfo( DA, NpesA, TotalPtsA )
      CALL DecompInfo( DB, NpesB, TotalPtsB )

      CALL MPI_COMM_SIZE( InComm, GroupSize, Ierror )
      CALL MPI_COMM_RANK( InComm, Iam, Ierror )
      CALL MPI_COMM_DUP( InComm, Pattern%Comm, Ierror )

#ifdef _SMEMORY
! Calculate info about the pattern 
      call ParCalcInfo(InComm,DA,DB, Info)
      TotalPtsA=Info%maxNumSeg
      TotalPtsB=Info%maxNumSeg
#endif

      Pattern%Size = GroupSize
      Pattern%Iam  = Iam
!
! Allocate the number of entries and list head arrays
!

!
! Allocate the patterns
!
      ALLOCATE( Pattern%SendDesc( NpesB ) )
      Pattern%SendDesc(:)%method = method
      if (iam .ge. NpesA) then
         do ipe = 1, NpesB
            ALLOCATE( Pattern%SendDesc(ipe)%Displacements(1) )
            ALLOCATE( Pattern%SendDesc(ipe)%BlockSizes(1) )
            Pattern%SendDesc(ipe)%Tot_Size = -1
            Pattern%SendDesc(ipe)%Nparcels = -1
            Pattern%SendDesc(ipe)%type = MPI_DATATYPE_NULL
            Pattern%SendDesc(ipe)%Displacements(1) = -1
            Pattern%SendDesc(ipe)%Blocksizes(1) = -1
         enddo
      endif

      ALLOCATE( Pattern%RecvDesc( NpesA ) )
      Pattern%RecvDesc(:)%method = method
      if (iam .ge. NpesB) then
         do ipe = 1, NpesA
            ALLOCATE( Pattern%RecvDesc(ipe)%Displacements(1) )
            ALLOCATE( Pattern%RecvDesc(ipe)%BlockSizes(1) )
            Pattern%RecvDesc(ipe)%Tot_Size = -1
            Pattern%RecvDesc(ipe)%Nparcels = -1
            Pattern%RecvDesc(ipe)%type = MPI_DATATYPE_NULL
            Pattern%RecvDesc(ipe)%Displacements(1) = -1
            Pattern%RecvDesc(ipe)%Blocksizes(1) = -1
         enddo
      endif

!
! Local allocations
!
      ALLOCATE( DisplacementsA( TotalPtsA ) )   ! Allocate for worst case
      ALLOCATE( BlockSizesA( TotalPtsA ) )      ! Allocate for worst case
      ALLOCATE( LocalA( TotalPtsA ) )           ! Allocate for worst case

      ALLOCATE( DisplacementsB( TotalPtsB ) )   ! Allocate for worst case
      ALLOCATE( BlockSizesB( TotalPtsB ) )      ! Allocate for worst case
      ALLOCATE( LocalB( TotalPtsB ) )           ! Allocate for worst case
      ALLOCATE( PeB( TotalPtsB ) )              ! Allocate for worst case

      ALLOCATE( Count( GroupSize ) )
      ALLOCATE( CountOut( GroupSize ) )

      JB        = 0
      Count     = 0
      LenB      = 0
      LocalA      = 0   !  (needed for parexchangevector later)
      BlocksizesA = 0   !  (needed for parexchangevector later)

      Num    = 0
      Inc    = 0

    if (iam .lt. NpesB) then

!
! Parse through all the tags in the local segment
      DO J = 1, SIZE( DB%Head(iam+1)%StartTags )
        OldPe     = -1         ! Set PE undefined
        DO Tag=DB%Head(iam+1)%StartTags(J), DB%Head(iam+1)%EndTags(J)
!
! Determine the index and PE of this entry on A. This might be inlined later
!
          CALL DecompGlobalToLocal( DA, Tag, Local, Pe )

!
! If ipe-1 is my id, then this is an entry ipe will receive from Pe
!
          IF ( Pe /= OldPe ) THEN
            OldPe   = Pe
            IF ( jb > 0 ) THEN
              BlockSizesB(jb) = LenB
              LenB = 0
            ENDIF
            jb = jb+1                     ! increment the segment index
            DisplacementsB(jb) = Inc      ! Zero-based offset of local segment
            LocalB(jb) = Local-1          ! The local index (zero-based)
            PeB(jb) = Pe                  ! Note the ID of the sender
            Count(Pe+1) = Count(Pe+1)+1 ! Increment counter of segments
          ENDIF
          LenB = LenB+1                   ! Good -- segment is getting longer
          Inc = Inc+1                     ! Increment local index
        ENDDO
      ENDDO
!
! Clean up
!
      IF ( jb>0 ) BlockSizesB(jb) = LenB
#if defined(DEBUG_PARPATTERNDECOMPTODECOMP)
      write(iulog,*) iam, "BlockSizes", BlockSizesB(1:jb), DisplacementsB(1:jb), PeB(1:jb), Count
#endif

      CPP_ASSERT_F90( JB .LE. TotalPtsB )
!
! Now create the pattern from the displacements and block sizes
!
      Inc = 0
      DO ipe = 1, NpesA
!
! Find the segments which are relevant for the sender ipe
! Make compact arrays BlockSizes and Displacements 
!
        DO j = 1, jb
          IF ( PeB(j) == ipe-1 ) THEN
            Inc = Inc + 1
            BlockSizesA(Inc) = BlockSizesB(j)
            DisplacementsA(Inc) = DisplacementsB(j)
            LocalA(Inc)      = LocalB(j)
          ENDIF
        ENDDO
      ENDDO
      CPP_ASSERT_F90( Inc .LE. TotalPtsA )

!
! Create the receiver communication pattern
!
      Off = 0
      DO ipe = 1, NpesA
        Num = Count(ipe)
        if(Num >0) then 
#if defined(DEBUG_PARPATTERNDECOMPTODECOMP)
        write(iulog,*) "Receiver Iam", Iam, "Ipe", Ipe-1, "Num", Num,         &
                 "Displacements", DisplacementsA(Off+1:Off+Num),        &
                 "BlockSizes", BlockSizesA(Off+1:Off+Num)
#endif
        endif
          IF ( Num > 0 .and. method > 0 ) THEN

            CALL MPI_TYPE_INDEXED( Num, BlockSizesA(Off+1),             &
                                   DisplacementsA(Off+1),               &
                                   DataType, Ptr, Ierror )
            CALL MPI_TYPE_COMMIT( Ptr, Ierror )
            Pattern%RecvDesc(ipe)%type = Ptr
          ELSE
            Pattern%RecvDesc(ipe)%type = MPI_DATATYPE_NULL
          ENDIF

          ALLOCATE( Pattern%RecvDesc(ipe)%Displacements(Num) )
          ALLOCATE( Pattern%RecvDesc(ipe)%BlockSizes(Num) )
          DO i=1, Num
            Pattern%RecvDesc(ipe)%Displacements(i) = DisplacementsA(i+Off)
            Pattern%RecvDesc(ipe)%BlockSizes(i)    = BlockSizesA(i+Off)
          ENDDO
          Pattern%RecvDesc(ipe)%Nparcels  = &
            size (Pattern%RecvDesc(ipe)%Displacements)
          Pattern%RecvDesc(ipe)%Tot_Size = &
            sum  (Pattern%RecvDesc(ipe)%Blocksizes)
          Max_Nparcels = max (Max_Nparcels, Pattern%RecvDesc(ipe)%Nparcels)

        Off = Off + Num
      ENDDO

    endif !  (iam .lt. NpesB)

!
! Now communicate what the receiver is expecting from the sender
!
      CALL ParExchangeVectorInt( InComm, Count, LocalA,                 &
                                 CountOut, DisplacementsB  )
      CALL ParExchangeVectorInt( InComm, Count, BlockSizesA,            &
                                 CountOut, BlockSizesB )

!
! Sender A: BlockSizes and Displacements can now be stored
!

    if (iam .lt. NpesA) then

      Off = 0
      DO ipe=1, NpesB
        Num = CountOut(ipe)
        if(Num>0) then 
#if defined(DEBUG_PARPATTERNDECOMPTODECOMP)
        write(iulog,*) "Sender Iam", Iam, "Ipe", Ipe-1, "Num", Num,           &
                 "Displacements", DisplacementsB(Off+1:Off+Num),        &
                 "BlockSizes", BlockSizesB(Off+1:Off+Num)
#endif
        endif
          IF ( Num > 0 .and. method > 0 ) THEN
            CALL MPI_TYPE_INDEXED( Num, BlockSizesB(Off+1),             &
                                   DisplacementsB(Off+1),               &
                                   DataType, Ptr, Ierror )
            CALL MPI_TYPE_COMMIT( Ptr, Ierror )
            Pattern%SendDesc(ipe)%type = Ptr
          ELSE
            Pattern%SendDesc(ipe)%type = MPI_DATATYPE_NULL
          ENDIF

          ALLOCATE( Pattern%SendDesc(ipe)%Displacements(Num) )
          ALLOCATE( Pattern%SendDesc(ipe)%BlockSizes(Num) )
          DO i=1, Num
            Pattern%SendDesc(ipe)%Displacements(i) = DisplacementsB(i+Off)
            Pattern%SendDesc(ipe)%BlockSizes(i)    = BlockSizesB(i+Off)
          ENDDO
          Pattern%SendDesc(ipe)%Nparcels  =  &
            size (Pattern%SendDesc(ipe)%Displacements)
          Pattern%SendDesc(ipe)%Tot_Size = &
            sum  (Pattern%SendDesc(ipe)%Blocksizes)
          Max_Nparcels = max (Max_Nparcels, Pattern%SendDesc(ipe)%Nparcels)

        Off = Off + Num
      ENDDO

    endif !  (iam .lt. NpesA)

      CALL get_partneroffset( InComm, Pattern%SendDesc, Pattern%RecvDesc )
      
      DEALLOCATE( CountOut )
      DEALLOCATE( Count )

      DEALLOCATE( PeB )
      DEALLOCATE( LocalB )
      DEALLOCATE( BlockSizesB )
      DEALLOCATE( DisplacementsB )

      DEALLOCATE( LocalA )
      DEALLOCATE( BlockSizesA )
      DEALLOCATE( DisplacementsA )

      CPP_LEAVE_PROCEDURE( "PARPATTERNDECOMPTODECOMP" )
      RETURN
!EOC
      END SUBROUTINE ParPatternDecompToDecomp
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE:   ParPatternDecompToGhost --- Create pattern decomp to ghost
!
! !INTERFACE:
      SUBROUTINE ParPatternDecompToGhost( InComm, DA, GB, Pattern, mod_method, T )
!
! !USES:
      USE decompmodule, ONLY : DecompType, DecompGlobalToLocal,         &
                               DecompInfo
      USE ghostmodule, ONLY : GhostType, GhostInfo
      USE mod_comm, ONLY : get_partneroffset
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER,  INTENT( IN )               :: InComm  ! # of PEs
      TYPE(DecompType),  INTENT( IN )      :: DA      ! Source Ghost Desc
      TYPE(GhostType),  INTENT( IN )       :: GB      ! Target Ghost Desc
      INTEGER,  INTENT( IN ), OPTIONAL     :: mod_method ! contiguous or derived type
      INTEGER, INTENT( IN ),  OPTIONAL     :: T       !

! !OUTPUT PARAMETERS:
      TYPE(ParPatternType), INTENT( OUT )  :: Pattern ! Comm Pattern
!
! !DESCRIPTION:
!     This routine contructs a communication pattern for a transformation
!     from decomposition to a ghosted decomposition, i.e., a so-called 
!     "transpose".  The resulting communication pattern can be used in 
!     ParBegin/EndTransfer with the decomposed arrays as inputs.  
!
! !SYSTEM ROUTINES:
!    MPI_COMM_SIZE,  MPI_COMM_RANK, MPI_COMM_DUP
!    MPI_TYPE_INDEXED, MPI_TYPE_COMMIT (depending on method)
!
! !REVISION HISTORY:
!   01.07.12   Sawyer     Creation from ParPatternDecompToDecomp
!   02.03.20   Sawyer     Bug fix: added OldLocal, increment Off
!   02.07.16   Sawyer     Added data type T
!   03.11.11   Mirin      Added optional argument mod_method
!   07.03.11   Mirin      Generalized to different sized decompositions
!   07.09.04   Dennis     Reduced amount of temporary memory usage
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER I, J, Tag, Local, Pe, LenB, JB, Ipe, Num, Inc, Off
      INTEGER Ptr                                ! Pointer type
      INTEGER GroupSize, Iam, Ierror
      INTEGER OldPe, OldLocal, TotalPtsA, NpesA
      INTEGER GlobalSizeB, LocalSizeB, BorderSizeB, NpesB
      INTEGER DataType
      INTEGER              :: method
      INTEGER              :: nCount, maxCount, ierr
#ifdef _SMEMORY
      TYPE (ParInfoType) :: Info
#endif

      INTEGER, ALLOCATABLE :: Count(:)           ! # segments for each recv PE
      INTEGER, ALLOCATABLE :: CountOut(:)        ! # segments for each send PE

      INTEGER, ALLOCATABLE :: DisplacementsA(:)  ! Generic displacements
      INTEGER, ALLOCATABLE :: BlockSizesA(:)     ! Generic block sizes
      INTEGER, ALLOCATABLE :: LocalA(:)          ! Generic Local indices

      INTEGER, ALLOCATABLE :: DisplacementsB(:)  ! Displacements for B
      INTEGER, ALLOCATABLE :: BlockSizesB(:)     ! Block sizes for B
      INTEGER, ALLOCATABLE :: LocalB(:)          ! Local indices for B
      INTEGER, ALLOCATABLE :: PeB(:)             ! Processor element numbers

      CPP_ENTER_PROCEDURE( "PARPATTERNDECOMPTOGHOST" )

      IF (present(T)) THEN
        DataType = T
      ELSE
        DataType = CPP_MPI_REAL8
      ENDIF

      IF (present(mod_method)) THEN
        method = mod_method
      ELSE
        method = 0     ! Default method - see mod_comm for description
      ENDIF

! Assume this routine is called by processes [ 0,max(NpesA,NpesB) )

      CALL DecompInfo( DA, NpesA, TotalPtsA )
      CALL GhostInfo( GB, NpesB, GlobalSizeB, LocalSizeB, BorderSizeB )

      CALL MPI_COMM_SIZE( InComm, GroupSize, Ierror )
      CALL MPI_COMM_RANK( InComm, Iam, Ierror )
      CALL MPI_COMM_DUP( InComm, Pattern%Comm, Ierror )

#ifdef _SMEMORY
! Calculate info about the pattern 
      call ParCalcInfo(InComm,DA,GB, Info)
      TotalPtsA=Info%maxNumSeg
      GlobalSizeB=Info%maxNumSeg
#endif

      Pattern%Size = GroupSize
      Pattern%Iam  = Iam
!
! Allocate the number of entries and list head arrays
!

!
! Allocate the patterns
!
      ALLOCATE( Pattern%SendDesc( NpesB ) )
      Pattern%SendDesc(:)%method = method
      if (iam .ge. NpesA) then
         do ipe = 1, NpesB
            ALLOCATE( Pattern%SendDesc(ipe)%Displacements(1) )
            ALLOCATE( Pattern%SendDesc(ipe)%BlockSizes(1) )
            Pattern%SendDesc(ipe)%Tot_Size = -1
            Pattern%SendDesc(ipe)%Nparcels = -1
            Pattern%SendDesc(ipe)%type = MPI_DATATYPE_NULL
            Pattern%SendDesc(ipe)%Displacements(1) = -1
            Pattern%SendDesc(ipe)%Blocksizes(1) = -1
         enddo
      endif

      ALLOCATE( Pattern%RecvDesc( NpesA ) )
      Pattern%RecvDesc(:)%method = method
      if (iam .ge. NpesB) then
         do ipe = 1, NpesA
            ALLOCATE( Pattern%RecvDesc(ipe)%Displacements(1) )
            ALLOCATE( Pattern%RecvDesc(ipe)%BlockSizes(1) )
            Pattern%RecvDesc(ipe)%Tot_Size = -1
            Pattern%RecvDesc(ipe)%Nparcels = -1
            Pattern%RecvDesc(ipe)%type = MPI_DATATYPE_NULL
            Pattern%RecvDesc(ipe)%Displacements(1) = -1
            Pattern%RecvDesc(ipe)%Blocksizes(1) = -1
         enddo
      endif

!
! Local allocations
!
      ALLOCATE( DisplacementsA( TotalPtsA ) )   ! Allocate for worst case
      ALLOCATE( BlockSizesA( TotalPtsA ) )      ! Allocate for worst case
      ALLOCATE( LocalA( TotalPtsA ) )           ! Allocate for worst case

      ALLOCATE( DisplacementsB( GlobalSizeB ) ) ! Allocate for worst case
      ALLOCATE( BlockSizesB( GlobalSizeB ) )    ! Allocate for worst case
      ALLOCATE( LocalB( GlobalSizeB ) )         ! Allocate for worst case
      ALLOCATE( PeB( GlobalSizeB ) )            ! Allocate for worst case

      ALLOCATE( Count( GroupSize ) )
      ALLOCATE( CountOut( GroupSize ) )

      JB        = 0
      Count     = 0
      LenB      = 0
      LocalA      = 0   !  (needed for parexchangevector later)
      BlocksizesA = 0   !  (needed for parexchangevector later)

      Num    = 0
      Inc    = 0

    if (iam .lt. NpesB) then

!
! Parse through all the tags in the local segment
      DO J = 1, SIZE( GB%Local%Head(iam+1)%StartTags )
        OldPe     = -1         ! Set PE undefined
        OldLocal  =  0         ! Set local index undefined
        DO Tag=GB%Local%Head(iam+1)%StartTags(J),                         &
                GB%Local%Head(iam+1)%EndTags(J)
          IF ( Tag > 0 ) THEN        ! Active point
!
! Determine the index and PE of this entry on A. This might be inlined later
!
            CALL DecompGlobalToLocal( DA, Tag, Local, Pe )

!
! If ipe-1 is my id, then this is an entry ipe will receive from Pe
!
            IF ( Pe /= OldPe .OR. Local /= OldLocal+1 ) THEN
              IF ( jb > 0 ) THEN
                BlockSizesB(jb) = LenB
                LenB = 0
              ENDIF
              jb = jb+1                     ! increment the segment index
              DisplacementsB(jb) = Inc      ! Zero-based offset of local segment
              LocalB(jb) = Local-1          ! Local indices (zero-based)
              PeB(jb) = Pe                  ! Note the ID of the sender
              Count(Pe+1) = Count(Pe+1)+1 ! Increment counter of segments
            ENDIF
            OldPe   = Pe                    ! Update PE
            OldLocal= Local                 ! Update local index
            LenB = LenB+1                   ! Good -- segment is getting longer
          ENDIF
          Inc = Inc+1                     ! Increment local index
        ENDDO
      ENDDO
!
! Clean up
!
      IF ( jb>0 ) BlockSizesB(jb) = LenB

      CPP_ASSERT_F90( JB .LE. GlobalSize )
!
! Now create the pattern from the displacements and block sizes
!
      Inc = 0
      DO ipe = 1, NpesA
!
! Find the segments which are relevant for the sender ipe
! Make compact arrays BlockSizes and Displacements 
!
        DO j = 1, jb
          IF ( PeB(j) == ipe-1 ) THEN
            Inc = Inc + 1
            BlockSizesA(Inc) = BlockSizesB(j)
            DisplacementsA(Inc) = DisplacementsB(j)
            LocalA(Inc)      = LocalB(j)
          ENDIF
        ENDDO
      ENDDO

      CPP_ASSERT_F90( Inc .LE. TotalPtsA )

      Off = 0
      DO ipe = 1, NpesA
        Num = Count(ipe)
#if defined( DEBUG_PARPATTERNDECOMPTOGHOST )
        write(iulog,*) "Receiver Iam", Iam, "Ipe", Ipe-1, "Num", Num, &
                 "Displacements", DisplacementsA(Off+1:Off+Num), &
                 "BlockSizes", BlockSizesA(Off+1:Off+Num)
#endif

!
! Create the receiver communication pattern
!
          IF ( Num > 0 .and. method > 0 ) THEN
            CALL MPI_TYPE_INDEXED( Num, BlockSizesA(Off+1),        &
                   DisplacementsA(Off+1), DataType, Ptr, Ierror )
            CALL MPI_TYPE_COMMIT( Ptr, Ierror )
            Pattern%RecvDesc(ipe)%type = Ptr
          ELSE
            Pattern%RecvDesc(ipe)%type = MPI_DATATYPE_NULL
          ENDIF

          ALLOCATE( Pattern%RecvDesc(ipe)%Displacements(Num) )
          ALLOCATE( Pattern%RecvDesc(ipe)%BlockSizes(Num) )
          DO i=1, Num
            Pattern%RecvDesc(ipe)%Displacements(i) = DisplacementsA(i+Off)
            Pattern%RecvDesc(ipe)%BlockSizes(i)    = BlockSizesA(i+Off)
          ENDDO
          Pattern%RecvDesc(ipe)%Nparcels  = &
            size (Pattern%RecvDesc(ipe)%Displacements)
          Pattern%RecvDesc(ipe)%Tot_Size = &
            sum  (Pattern%RecvDesc(ipe)%Blocksizes)
          Max_Nparcels = max (Max_Nparcels, Pattern%RecvDesc(ipe)%Nparcels)

        Off = Off + Num
      ENDDO

    endif !  (iam .lt. NpesB)

!
! Now communicate what the receiver is expecting to the sender
!
      CALL ParExchangeVectorInt( InComm, Count, LocalA,                 &
                                 CountOut, DisplacementsB  )
      CALL ParExchangeVectorInt( InComm, Count, BlockSizesA,            &
                                 CountOut, BlockSizesB )

!
! Sender A: BlockSizes and Displacements can now be stored
!

    if (iam .lt. NpesA) then

      Off = 0
      DO ipe=1, NpesB
        Num = CountOut(ipe)
#if defined( DEBUG_PARPATTERNDECOMPTOGHOST )
        write(iulog,*) "Sender Iam", Iam, "Ipe", Ipe-1, "Num", Num,           &
                 "Displacements", DisplacementsB(Off+1:Off+Num),        &
                 "BlockSizes", BlockSizesB(Off+1:Off+Num)
#endif

          IF ( Num > 0 .and. method > 0 ) THEN
            CALL MPI_TYPE_INDEXED( Num, BlockSizesB(Off+1),          &
                    DisplacementsB(Off+1), DataType, Ptr, Ierror )
            CALL MPI_TYPE_COMMIT( Ptr, Ierror )
            Pattern%SendDesc(ipe)%type = Ptr
          ELSE
            Pattern%SendDesc(ipe)%type = MPI_DATATYPE_NULL
          ENDIF

          ALLOCATE( Pattern%SendDesc(ipe)%Displacements(Num) )
          ALLOCATE( Pattern%SendDesc(ipe)%BlockSizes(Num) )
          DO i=1, Num
            Pattern%SendDesc(ipe)%Displacements(i) = DisplacementsB(i+Off)
            Pattern%SendDesc(ipe)%BlockSizes(i)    = BlockSizesB(i+Off)
          ENDDO
          Pattern%SendDesc(ipe)%Nparcels  = &
            size (Pattern%SendDesc(ipe)%Displacements)
          Pattern%SendDesc(ipe)%Tot_Size = &
            sum  (Pattern%SendDesc(ipe)%Blocksizes)
          Max_Nparcels = max (Max_Nparcels, Pattern%SendDesc(ipe)%Nparcels)

        Off = Off + Num
      ENDDO

    endif !  (iam .lt. NpesA)

      CALL get_partneroffset( InComm, Pattern%SendDesc, Pattern%RecvDesc )
      
      DEALLOCATE( CountOut )
      DEALLOCATE( Count )

      DEALLOCATE( PeB )
      DEALLOCATE( LocalB )
      DEALLOCATE( BlockSizesB )
      DEALLOCATE( DisplacementsB )

      DEALLOCATE( LocalA )
      DEALLOCATE( BlockSizesA )
      DEALLOCATE( DisplacementsA )

      CPP_LEAVE_PROCEDURE( "PARPATTERNDECOMPTOGHOST" )
      RETURN
!EOC
      END SUBROUTINE ParPatternDecompToGhost
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE:   ParPatternGhostToDecomp --- Create pattern between decomps
!
! !INTERFACE:
      SUBROUTINE ParPatternGhostToDecomp( InComm, GA, DB, Pattern, mod_method, T )
!
! !USES:
      USE decompmodule, ONLY : DecompType, DecompGlobalToLocal, DecompInfo
      USE ghostmodule, ONLY : GhostType, GhostInfo
      USE mod_comm, ONLY : get_partneroffset
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER,  INTENT( IN )               :: InComm  ! # of PEs
      TYPE(GhostType),   INTENT( IN )      :: GA      ! Source Decomp Desc
      TYPE(DecompType),  INTENT( IN )      :: DB      ! Target Decomp Desc
      INTEGER,  INTENT( IN ), OPTIONAL     :: mod_method ! contiguous or derived type
      INTEGER, INTENT( IN ),  OPTIONAL     :: T       !
! !OUTPUT PARAMETERS:
      TYPE(ParPatternType), INTENT( OUT )  :: Pattern ! Comm Pattern
!
! !DESCRIPTION:
!     This routine contructs a communication pattern for a 
!     transformation from one ghosted decomposition to partitioned
!     one, i.e., a so-called "transpose". The resulting communication 
!     pattern can be used in ParBegin/EndTransfer with the decomposed 
!     arrays as inputs.  
!
! !SYSTEM ROUTINES:
!    MPI_COMM_SIZE,  MPI_COMM_RANK, MPI_COMM_DUP
!    MPI_TYPE_INDEXED, MPI_TYPE_COMMIT (depending on method)
!
! !REVISION HISTORY:
!   02.01.10   Sawyer     Creation from DecompToDecomp
!   02.07.16   Sawyer     Added data type T
!   03.11.11   Mirin      Added optional argument mod_method
!   07.03.11   Mirin      Generalized to different sized decompositions
!   07.09.04   Dennis     Reduced amount of temporary memory usage
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER I, J, Tag, Local, Pe, Len, JA, Ipe, Num, Inc, Off
      INTEGER NpesA, GlobalSizeA, LocalSizeA, BorderSizeA
      INTEGER OldPe, OldLocal, TotalPtsB, NpesB
      INTEGER GroupSize, Iam, Ierror
      INTEGER Ptr                                ! Pointer type
      INTEGER  DataType
      INTEGER              :: method
      INTEGER              :: nCount, maxCount, ierr
#ifdef _SMEMORY
      TYPE (ParInfoType) :: Info
#endif

      INTEGER, ALLOCATABLE :: Count(:)           ! # segments for each recv PE
      INTEGER, ALLOCATABLE :: CountOut(:)        ! # segments for each send PE

      INTEGER, ALLOCATABLE :: DisplacementsA(:)  ! Generic displacements
      INTEGER, ALLOCATABLE :: BlockSizesA(:)     ! Generic block sizes
      INTEGER, ALLOCATABLE :: GlobalA(:)          ! Generic Local indices
      INTEGER, ALLOCATABLE :: PeA(:)             ! Processor element numbers

      INTEGER, ALLOCATABLE :: DisplacementsB(:)  ! Displacements for B
      INTEGER, ALLOCATABLE :: BlockSizesB(:)     ! Block sizes for B
      INTEGER, ALLOCATABLE :: GlobalB(:)         ! Global indices for B

      CPP_ENTER_PROCEDURE( "PARPATTERNGHOSTTODECOMP" )

      IF (present(T)) THEN
        DataType = T
      ELSE
        DataType = CPP_MPI_REAL8
      ENDIF

      IF (present(mod_method)) THEN
        method = mod_method
      ELSE
        method = 0     ! Default method - see mod_comm for description
      ENDIF

! Assume this routine is called by processes [ 0,max(NpesA,NpesB) )

      CALL GhostInfo( GA, NpesA, GlobalSizeA, LocalSizeA, BorderSizeA )
      CALL DecompInfo( DB, NpesB, TotalPtsB )

      CALL MPI_COMM_SIZE( InComm, GroupSize, Ierror )
      CALL MPI_COMM_RANK( InComm, Iam, Ierror )
      CALL MPI_COMM_DUP( InComm, Pattern%Comm, Ierror )

#ifdef _SMEMORY
! Calculate info about the pattern 
      call ParCalcInfo(InComm,GA,DB, Info)
      GlobalSizeA=Info%maxNumSeg
      TotalPtsB=Info%maxNumSeg
#endif

      Pattern%Size = GroupSize
      Pattern%Iam  = Iam
!
! Allocate the number of entries and list head arrays
!

!
! Allocate the patterns
!
      ALLOCATE( Pattern%SendDesc( NpesB ) )
      Pattern%SendDesc(:)%method = method
      if (iam .ge. NpesA) then
         do ipe = 1, NpesB
            ALLOCATE( Pattern%SendDesc(ipe)%Displacements(1) )
            ALLOCATE( Pattern%SendDesc(ipe)%BlockSizes(1) )
            Pattern%SendDesc(ipe)%Tot_Size = -1
            Pattern%SendDesc(ipe)%Nparcels = -1
            Pattern%SendDesc(ipe)%type = MPI_DATATYPE_NULL
            Pattern%SendDesc(ipe)%Displacements(1) = -1
            Pattern%SendDesc(ipe)%Blocksizes(1) = -1
         enddo
      endif

      ALLOCATE( Pattern%RecvDesc( NpesA ) )
      Pattern%RecvDesc(:)%method = method
      if (iam .ge. NpesB) then
         do ipe = 1, NpesA
            ALLOCATE( Pattern%RecvDesc(ipe)%Displacements(1) )
            ALLOCATE( Pattern%RecvDesc(ipe)%BlockSizes(1) )
            Pattern%RecvDesc(ipe)%Tot_Size = -1
            Pattern%RecvDesc(ipe)%Nparcels = -1
            Pattern%RecvDesc(ipe)%type = MPI_DATATYPE_NULL
            Pattern%RecvDesc(ipe)%Displacements(1) = -1
            Pattern%RecvDesc(ipe)%Blocksizes(1) = -1
         enddo
      endif

!
! Local allocations
!
      ALLOCATE( DisplacementsA( GlobalSizeA ) )   ! Allocate for worst case
      ALLOCATE( BlockSizesA( GlobalSizeA ) )      ! Allocate for worst case
      ALLOCATE( GlobalA( GlobalSizeA ) )          ! Allocate for worst case
      ALLOCATE( PeA( GlobalSizeA ) )              ! Allocate for worst case

      ALLOCATE( DisplacementsB( TotalPtsB ) )     ! Allocate for worst case
      ALLOCATE( BlockSizesB( TotalPtsB ) )        ! Allocate for worst case
      ALLOCATE( GlobalB( TotalPtsB ) )            ! Allocate for worst case

      ALLOCATE( Count( GroupSize ) )
      ALLOCATE( CountOut( GroupSize ) )

      JA     = 0
      Count  = 0
      Len    = 0
      GlobalB     = 0   !  (needed for parexchangevector later)
      BlockSizesB = 0   !  (needed for parexchangevector later)

      Num    = 0
      Inc    = 0

    if (iam .lt. NpesB) then

!
! Parse through all the tags in the local segment
      DO J = 1, SIZE( DB%Head(iam+1)%StartTags )
        OldPe     = -1         ! Set PE undefined
        OldLocal  = 0          ! Set index value undefined
        DO Tag=DB%Head(iam+1)%StartTags(J), DB%Head(iam+1)%EndTags(J)
!
! Determine the index and PE of this entry on A. This might be inlined later
!
          CALL DecompGlobalToLocal( GA%Decomp, Tag, Local, Pe )

!
! If ipe-1 is my id, then this is an entry ipe will receive from Pe
!
          IF ( Pe /= OldPe  .OR. Local /= OldLocal+1 ) THEN
            IF ( ja > 0 ) THEN
              BlockSizesA(ja) = Len
              Len = 0
            ENDIF
            ja = ja+1                     ! increment the segment index
            DisplacementsA(ja) = Inc      ! Zero-based offset of local segment
            GlobalA(ja) = Tag             ! The global tag of the desired datum
            PeA(ja) = Pe                  ! Note the ID of the sender
            Count(Pe+1) = Count(Pe+1)+1   ! Increment counter of segments
          ENDIF
          OldPe    = Pe                   ! Update old PE
          OldLocal = Local                ! Update old local index
          Len = Len+1                     ! Good -- segment is getting longer
          Inc = Inc+1                     ! Increment local index
        ENDDO
      ENDDO
!
! Clean up
!
      BlockSizesA(ja) = Len
      CPP_ASSERT_F90( JA .LE. GlobalSizeA )
!
! Now create the pattern from the displacements and block sizes
!
      Inc = 0
      DO ipe = 1, NpesA
!
! Find the segments which are relevant for the sender ipe
! Make compact arrays BlockSizes and Displacements 
!
        DO j = 1, ja
          IF ( PeA(j) == ipe-1 ) THEN
            Inc = Inc + 1
            BlockSizesB(Inc) = BlockSizesA(j)
            DisplacementsB(Inc) = DisplacementsA(j)
            GlobalB(Inc)      = GlobalA(j)
          ENDIF
        ENDDO
      ENDDO

     CPP_ASSERT_F90(Inc .LE. TotalPtsB)

!
! Create the receiver communication pattern
!
      Off = 0
      DO ipe = 1, NpesA
        Num = Count(ipe)
#if defined( DEBUG_PARPATTERNGHOSTTODECOMP )
        write(iulog,*) "Receiver Iam", Iam, "Ipe", Ipe-1, "Num", Num, &
                 "Displacements", DisplacementsB(Off+1:Off+Num), &
                 "BlockSizes", BlockSizesB(Off+1:Off+Num)
#endif

          IF ( Num > 0 .and. method > 0 ) THEN
            CALL MPI_TYPE_INDEXED( Num, BlockSizesB(Off+1),       &
               DisplacementsB(Off+1), DataType, Ptr, Ierror )
            CALL MPI_TYPE_COMMIT( Ptr, Ierror )
            Pattern%RecvDesc(ipe)%type = Ptr
          ELSE
            Pattern%RecvDesc(ipe)%type = MPI_DATATYPE_NULL
          ENDIF

          ALLOCATE( Pattern%RecvDesc(ipe)%Displacements(Num) )
          ALLOCATE( Pattern%RecvDesc(ipe)%BlockSizes(Num) )
          DO i=1, Num
            Pattern%RecvDesc(ipe)%Displacements(i) = DisplacementsB(i+Off)
            Pattern%RecvDesc(ipe)%BlockSizes(i)    = BlockSizesB(i+Off)
          ENDDO
          Pattern%RecvDesc(ipe)%Nparcels  = &
            size (Pattern%RecvDesc(ipe)%Displacements)
          Pattern%RecvDesc(ipe)%Tot_Size = &
            sum  (Pattern%RecvDesc(ipe)%Blocksizes)
          Max_Nparcels = max (Max_Nparcels, Pattern%RecvDesc(ipe)%Nparcels)

        Off = Off + Num
      ENDDO

    endif !  (iam .lt. NpesB)

!
! Now communicate what the receiver is expecting to the sender
!
      CALL ParExchangeVectorInt( InComm, Count, GlobalB,                 &
                                 CountOut, GlobalA  )
      CALL ParExchangeVectorInt( InComm, Count, BlockSizesB,            &
                                 CountOut, BlockSizesA )

    if (iam .lt. NpesA) then

!
! Sender A: BlockSizes and Displacements can now be stored
!
      Off = 0
      DO ipe=1, NpesB
        Num = CountOut(ipe)
        DO i=1, Num
          CALL DecompGlobalToLocal( GA%Local, GlobalA(i+Off), Local, Pe )
          DisplacementsA(i+Off) = Local-1    ! zero-based displacement
        ENDDO
#if defined( DEBUG_PARPATTERNGHOSTTODECOMP )
        write(iulog,*) "Sender Iam", Iam, "Ipe", Ipe-1, "Num", Num,  &
                 "Displacements", DisplacementsA(Off+1:Off+Num), &
                 "BlockSizes", BlockSizesA(Off+1:Off+Num)
#endif

          IF ( Num > 0 .and. method > 0 ) THEN
            CALL MPI_TYPE_INDEXED( Num, BlockSizesA(Off+1),        &
                    DisplacementsA(Off+1), DataType, Ptr, Ierror )
            CALL MPI_TYPE_COMMIT( Ptr, Ierror )
            Pattern%SendDesc(ipe)%type = Ptr
          ELSE
            Pattern%SendDesc(ipe)%type = MPI_DATATYPE_NULL
          ENDIF

          ALLOCATE( Pattern%SendDesc(ipe)%Displacements(Num) )
          ALLOCATE( Pattern%SendDesc(ipe)%BlockSizes(Num) )
          DO i=1, Num
            Pattern%SendDesc(ipe)%Displacements(i) = DisplacementsA(i+Off)
            Pattern%SendDesc(ipe)%BlockSizes(i)    = BlockSizesA(i+Off)
          ENDDO
          Pattern%SendDesc(ipe)%Nparcels  = &
            size (Pattern%SendDesc(ipe)%Displacements)
          Pattern%SendDesc(ipe)%Tot_Size = &
            sum  (Pattern%SendDesc(ipe)%Blocksizes)
          Max_Nparcels = max (Max_Nparcels, Pattern%SendDesc(ipe)%Nparcels)

        Off = Off + Num
      ENDDO

    endif !  (iam .lt. NpesA)

      CALL get_partneroffset( InComm, Pattern%SendDesc, Pattern%RecvDesc )
      
      DEALLOCATE( CountOut )
      DEALLOCATE( Count )

      DEALLOCATE( PeA )
      DEALLOCATE( GlobalA )
      DEALLOCATE( BlockSizesA )
      DEALLOCATE( DisplacementsA )

      DEALLOCATE( GlobalB )
      DEALLOCATE( BlockSizesB )
      DEALLOCATE( DisplacementsB )

      CPP_LEAVE_PROCEDURE( "PARPATTERNGHOSTTODECOMP" )
      RETURN
!EOC
      END SUBROUTINE ParPatternGhostToDecomp
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE:   ParPatternGhostToGhost --- Create pattern between decomps
!
! !INTERFACE:
      SUBROUTINE ParPatternGhostToGhost( InComm, GA, GB, Pattern, mod_method, T )
!
! !USES:
      USE decompmodule, ONLY : DecompGlobalToLocal
      USE ghostmodule, ONLY : GhostType, GhostInfo
      USE mod_comm, ONLY : get_partneroffset
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER,  INTENT( IN )               :: InComm  ! # of PEs
      TYPE(GhostType),   INTENT( IN )      :: GA      ! Source Ghost Decomp
      TYPE(GhostType),   INTENT( IN )      :: GB      ! Target Ghost Decomp
      INTEGER,  INTENT( IN ), OPTIONAL     :: mod_method ! contiguous or derived type
      INTEGER, INTENT( IN ),  OPTIONAL     :: T       !
! !OUTPUT PARAMETERS:
      TYPE(ParPatternType), INTENT( OUT )  :: Pattern ! Comm Pattern
!
! !DESCRIPTION:
!     This routine contructs a communication pattern for a 
!     transformation from one ghosted decomposition to partitioned
!     one, i.e., a so-called "transpose". The resulting communication 
!     pattern can be used in ParBegin/EndTransfer with the decomposed 
!     arrays as inputs.  
!
! !SYSTEM ROUTINES:
!    MPI_COMM_SIZE,  MPI_COMM_RANK, MPI_COMM_DUP
!    MPI_TYPE_INDEXED, MPI_TYPE_COMMIT (depending on method)
!
! !REVISION HISTORY:
!   02.01.10   Sawyer     Creation from DecompToDecomp
!   02.07.16   Sawyer     Added data type T
!   03.11.11   Mirin      Added optional argument mod_method
!   07.03.11   Mirin      Generalized to different sized decompositions
!   07.09.04   Dennis     Reduced amount of temporary memory usage
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER I, J, Tag, Local, Pe, Len, JA, Ipe, Num, Inc, Off
      INTEGER NpesA, GlobalSizeA, LocalSizeA, BorderSizeA
      INTEGER NpesB, GlobalSizeB, LocalSizeB, BorderSizeB
      INTEGER GroupSize, Iam, Ierror, OldPe, OldLocal 
      INTEGER Ptr                                ! Pointer type
      INTEGER  DataType
      INTEGER              :: method
      INTEGER              :: nCount, maxCount, ierr
#ifdef _SMEMORY
      TYPE (ParInfoType) :: Info
#endif

      INTEGER, ALLOCATABLE :: Count(:)           ! # segments for each recv PE
      INTEGER, ALLOCATABLE :: CountOut(:)        ! # segments for each send PE

      INTEGER, ALLOCATABLE :: DisplacementsA(:)  ! Generic displacements
      INTEGER, ALLOCATABLE :: BlockSizesA(:)     ! Generic block sizes
      INTEGER, ALLOCATABLE :: GlobalA(:)         ! Generic Local indices
      INTEGER, ALLOCATABLE :: PeA(:)             ! Processor element numbers

      INTEGER, ALLOCATABLE :: DisplacementsB(:)  ! Displacements for B
      INTEGER, ALLOCATABLE :: BlockSizesB(:)     ! Block sizes for B
      INTEGER, ALLOCATABLE :: GlobalB(:)         ! Global indices for B

      CPP_ENTER_PROCEDURE( "PARPATTERNGHOSTTOGHOST" )

      IF (present(T)) THEN
        DataType = T
      ELSE
        DataType = CPP_MPI_REAL8
      ENDIF

      IF (present(mod_method)) THEN
        method = mod_method
      ELSE
        method = 0     ! Default method - see mod_comm for description
      ENDIF

! Assume this routine is called by processes [ 0,max(NpesA,NpesB) )

      CALL GhostInfo( GA, NpesA, GlobalSizeA, LocalSizeA, BorderSizeA )
      CALL GhostInfo( GB, NpesB, GlobalSizeB, LocalSizeB, BorderSizeB )

      CALL MPI_COMM_SIZE( InComm, GroupSize, Ierror )
      CALL MPI_COMM_RANK( InComm, Iam, Ierror )
      CALL MPI_COMM_DUP( InComm, Pattern%Comm, Ierror )

#ifdef _SMEMORY
! Calculate info about the pattern 
      call ParCalcInfo(InComm,GA,GB, Info)
      GlobalSizeA=Info%maxNumSeg
      GlobalSizeB=Info%maxNumSeg
#endif

      Pattern%Size = GroupSize
      Pattern%Iam  = Iam
!
! Allocate the number of entries and list head arrays
!

!
! Allocate the patterns
!
      ALLOCATE( Pattern%SendDesc( NpesB ) )
      Pattern%SendDesc(:)%method = method
      if (iam .ge. NpesA) then
         do ipe = 1, NpesB
            ALLOCATE( Pattern%SendDesc(ipe)%Displacements(1) )
            ALLOCATE( Pattern%SendDesc(ipe)%BlockSizes(1) )
            Pattern%SendDesc(ipe)%Tot_Size = -1
            Pattern%SendDesc(ipe)%Nparcels = -1
            Pattern%SendDesc(ipe)%type = MPI_DATATYPE_NULL
            Pattern%SendDesc(ipe)%Displacements(1) = -1
            Pattern%SendDesc(ipe)%Blocksizes(1) = -1
         enddo
      endif

      ALLOCATE( Pattern%RecvDesc( NpesA ) )
      Pattern%RecvDesc(:)%method = method
      if (iam .ge. NpesB) then
         do ipe = 1, NpesA
            ALLOCATE( Pattern%RecvDesc(ipe)%Displacements(1) )
            ALLOCATE( Pattern%RecvDesc(ipe)%BlockSizes(1) )
            Pattern%RecvDesc(ipe)%Tot_Size = -1
            Pattern%RecvDesc(ipe)%Nparcels = -1
            Pattern%RecvDesc(ipe)%type = MPI_DATATYPE_NULL
            Pattern%RecvDesc(ipe)%Displacements(1) = -1
            Pattern%RecvDesc(ipe)%Blocksizes(1) = -1
         enddo
      endif

!
! Local allocations
!
      ALLOCATE( DisplacementsA( GlobalSizeA ) )   ! Allocate for worst case
      ALLOCATE( BlockSizesA( GlobalSizeA ) )      ! Allocate for worst case
      ALLOCATE( GlobalA( GlobalSizeA ) )          ! Allocate for worst case
      ALLOCATE( PeA( GlobalSizeA ) )              ! Allocate for worst case

      ALLOCATE( DisplacementsB( GlobalSizeB ) )   ! Allocate for worst case
      ALLOCATE( BlockSizesB( GlobalSizeB ) )      ! Allocate for worst case
      ALLOCATE( GlobalB( GlobalSizeB ) )          ! Allocate for worst case

      ALLOCATE( Count( GroupSize ) )
      ALLOCATE( CountOut( GroupSize ) )

      JA        = 0
      Count     = 0
      Len      = 0
      GlobalB     = 0   !  (needed for parexchangevector later)
      BlocksizesB = 0   !  (needed for parexchangevector later)

      Num    = 0
      Inc    = 0

    if (iam .lt. NpesB) then

!
! Parse through all the tags in the local segment
      DO J = 1, SIZE( GB%Local%Head(iam+1)%StartTags )
        OldPe     = -1         ! Set PE undefined
        OldLocal  = 0          ! Set index value undefined
        DO Tag=GB%Local%Head(iam+1)%StartTags(J), GB%Local%Head(iam+1)%EndTags(J)
          IF ( Tag > 0 ) THEN       ! Active point
!
! Determine the index and PE of this entry on A. This might be inlined later
!
            CALL DecompGlobalToLocal( GA%Decomp, Tag, Local, Pe )
!
! If ipe-1 is my id, then this is an entry ipe will receive from Pe
!
            IF ( Pe /= OldPe  .OR. Local /= OldLocal+1 ) THEN
              IF ( ja > 0 ) THEN
                BlockSizesA(ja) = Len
                Len = 0
              ENDIF
              ja = ja+1                     ! increment the segment index
              DisplacementsA(ja) = Inc      ! Zero-based offset of local segment
              GlobalA(ja) = Tag             ! The global tag of the desired datum
              PeA(ja) = Pe                  ! Note the ID of the sender
              Count(Pe+1) = Count(Pe+1)+1   ! Increment counter of segments
            ENDIF
            OldPe   = Pe                    ! Update old PE
            OldLocal = Local                ! Update old local index
            Len = Len+1                     ! Good -- segment is getting longer
          ENDIF
          Inc = Inc+1                       ! Increment local index
        ENDDO
      ENDDO
!
! Clean up
!
      BlockSizesA(ja) = Len

      CPP_ASSERT_F90( JA .LE. GlobalSizeA )

!
! Now create the pattern from the displacements and block sizes
!
      Inc = 0
      DO ipe = 1, NpesA
!
! Find the segments which are relevant for the sender ipe
! Make compact arrays BlockSizes and Displacements 
!
        DO j = 1, ja
          IF ( PeA(j) == ipe-1 ) THEN
            Inc = Inc + 1
            BlockSizesB(Inc) = BlockSizesA(j)
            DisplacementsB(Inc) = DisplacementsA(j)
            GlobalB(Inc)      = GlobalA(j)
          ENDIF
        ENDDO
      ENDDO
      CPP_ASSERT_F90( Inc .LE. GlobalSizeB )

!
! Create the receiver communication pattern
!
      Off = 0
      DO ipe = 1, NpesA
        Num = Count(ipe)
#if defined(DEBUG_PARPATTERNGHOSTTOGHOST)
        write(iulog,*) "Receiver Iam", Iam, "Ipe", Ipe-1, "Num", Num, &
                 "Displacements", DisplacementsB(Off+1:Off+Num), &
                 "BlockSizes", BlockSizesB(Off+1:Off+Num)
#endif

          IF ( Num > 0 .and. method > 0 ) THEN
            CALL MPI_TYPE_INDEXED( Num, BlockSizesB(Off+1),         &
                 DisplacementsB(Off+1), DataType, Ptr, Ierror )
            CALL MPI_TYPE_COMMIT( Ptr, Ierror )
            Pattern%RecvDesc(ipe)%type = Ptr
          ELSE
            Pattern%RecvDesc(ipe)%type = MPI_DATATYPE_NULL
          ENDIF

          ALLOCATE( Pattern%RecvDesc(ipe)%Displacements(Num) )
          ALLOCATE( Pattern%RecvDesc(ipe)%BlockSizes(Num) )
          DO i=1, Num
            Pattern%RecvDesc(ipe)%Displacements(i) = DisplacementsB(i+Off)
            Pattern%RecvDesc(ipe)%BlockSizes(i)    = BlockSizesB(i+Off)
          ENDDO
          Pattern%RecvDesc(ipe)%Nparcels  = &
            size (Pattern%RecvDesc(ipe)%Displacements)
          Pattern%RecvDesc(ipe)%Tot_Size = &
            sum  (Pattern%RecvDesc(ipe)%Blocksizes)
          Max_Nparcels = max (Max_Nparcels, Pattern%RecvDesc(ipe)%Nparcels)

        Off = Off + Num
      ENDDO

    endif !  (iam .lt. NpesB)

!
! Now communicate what the receiver is expecting to the sender
!
      CALL ParExchangeVectorInt( InComm, Count, GlobalB,                &
                                 CountOut, GlobalA  )
      CALL ParExchangeVectorInt( InComm, Count, BlockSizesB,            &
                                 CountOut, BlockSizesA )

    if (iam .lt. NpesA) then

!
! Sender A: BlockSizes and Displacements can now be stored
!
      Off = 0
      DO ipe=1, NpesB
        Num = CountOut(ipe)
        DO i=1, Num
          CALL DecompGlobalToLocal( GA%Local, GlobalA(i+Off), Local, Pe )
          DisplacementsA(i+Off) = Local-1    ! zero-based displacement
        ENDDO
#if defined(DEBUG_PARPATTERNGHOSTTOGHOST)
        write(iulog,*) "Sender Iam", Iam, "Ipe", Ipe-1, "Num", Num,  &
                 "Displacements", DisplacementsA(Off+1:Off+Num), &
                 "BlockSizes", BlockSizesA(Off+1:Off+Num)
#endif

          IF ( Num > 0 .and. method > 0 ) THEN
            CALL MPI_TYPE_INDEXED( Num, BlockSizesA(Off+1),        &
                 DisplacementsA(Off+1), DataType, Ptr, Ierror )
            CALL MPI_TYPE_COMMIT( Ptr, Ierror )
            Pattern%SendDesc(ipe)%type = Ptr
          ELSE
            Pattern%SendDesc(ipe)%type = MPI_DATATYPE_NULL
          ENDIF

          ALLOCATE( Pattern%SendDesc(ipe)%Displacements(Num) )
          ALLOCATE( Pattern%SendDesc(ipe)%BlockSizes(Num) )
          DO i=1, Num
            Pattern%SendDesc(ipe)%Displacements(i) = DisplacementsA(i+Off)
            Pattern%SendDesc(ipe)%BlockSizes(i)    = BlockSizesA(i+Off)
          ENDDO
          Pattern%SendDesc(ipe)%Nparcels  = &
            size (Pattern%SendDesc(ipe)%Displacements)
          Pattern%SendDesc(ipe)%Tot_Size = &
            sum  (Pattern%SendDesc(ipe)%Blocksizes)
          Max_Nparcels = max (Max_Nparcels, Pattern%SendDesc(ipe)%Nparcels)

        Off = Off + Num
      ENDDO

    endif !  (iam .lt. NpesA)

      CALL get_partneroffset( InComm, Pattern%SendDesc, Pattern%RecvDesc )


      DEALLOCATE( CountOut )
      DEALLOCATE( Count )

      DEALLOCATE( PeA )
      DEALLOCATE( GlobalA )
      DEALLOCATE( BlockSizesA )
      DEALLOCATE( DisplacementsA )

      DEALLOCATE( GlobalB )
      DEALLOCATE( BlockSizesB )
      DEALLOCATE( DisplacementsB )

      CPP_LEAVE_PROCEDURE( "PARPATTERNGHOSTTOGHOST" )
      RETURN
!EOC
      END SUBROUTINE ParPatternGhostToGhost
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE:   ParPatternFree --- Free the communication pattern
!
! !INTERFACE:
      SUBROUTINE ParPatternFree( InComm, Pattern )
!
! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER,  INTENT( IN )                 :: InComm  ! # of PEs
! !INPUT/OUTPUT PARAMETERS:
      TYPE(ParPatternType), INTENT( INOUT )  :: Pattern ! Comm Pattern
!
! !DESCRIPTION:
!     This routine frees a communication pattern.  
!
! !SYSTEM ROUTINES:
!     MPI_TYPE_FREE
!
! !BUGS:
!     The MPI_TYPE_FREE statement does not seem to work with FFC
!
! !REVISION HISTORY:
!   01.02.10   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER  ipe, GroupSize, Pointer, Ierror, method

      CPP_ENTER_PROCEDURE( "PARPATTERNFREE" )

      method = Pattern%RecvDesc(1)%method

!
! First request the needed ghost values from other processors.
!
! Free all the MPI derived types
!
      DO ipe=1, Pattern%Size
        Pointer = Pattern%SendDesc(ipe)%type
        IF ( Pointer /= MPI_DATATYPE_NULL ) THEN
          CALL MPI_TYPE_FREE( Pointer, Ierror )
        ENDIF
        Pointer = Pattern%RecvDesc(ipe)%type
        IF ( Pointer /= MPI_DATATYPE_NULL ) THEN
          CALL MPI_TYPE_FREE( Pointer, Ierror )
        ENDIF
      ENDDO

      DO ipe=1, size(Pattern%RecvDesc)
        DEALLOCATE( Pattern%RecvDesc(ipe)%Displacements )
        DEALLOCATE( Pattern%RecvDesc(ipe)%BlockSizes )
      ENDDO
      DO ipe=1, size(Pattern%SendDesc)
        DEALLOCATE( Pattern%SendDesc(ipe)%Displacements )
        DEALLOCATE( Pattern%SendDesc(ipe)%BlockSizes )
      ENDDO

      DEALLOCATE( Pattern%SendDesc )
      DEALLOCATE( Pattern%RecvDesc )

      CPP_LEAVE_PROCEDURE( "PARPATTERNFREE" )
      RETURN
!EOC
      END SUBROUTINE ParPatternFree
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParScatterReal --- Scatter slice to all PEs
!
! !INTERFACE:
      SUBROUTINE ParScatterReal ( InComm, Root, Slice, Decomp, Local )

! !USES:
      USE decompmodule, ONLY:  DecompType, Lists
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )          :: InComm       ! Communicator
      INTEGER, INTENT( IN )          :: Root         ! Root PE
      REAL(CPP_REAL8), INTENT( IN )  :: Slice(*)     ! Global Slice
      TYPE(DecompType), INTENT( IN ) :: Decomp       ! Decomp information

! !OUTPUT PARAMETERS:
      REAL(CPP_REAL8), INTENT( OUT ) :: Local(*)     ! Local Slice

! !DESCRIPTION:
!     Given a decomposition of the domain, dole out a slice 
!     (one-dimensional array) to all the constituent PEs as described
!     by the decomposition Decomp.
!
!
! !SYSTEM ROUTINES:
!     MPI_ISEND, MPI_RECV, MPI_COMM_RANK
!
! !REVISION HISTORY:
!   97.04.14   Sawyer     Creation
!   97.04.16   Sawyer     Cleaned up for walk-through
!   97.05.01   Sawyer     Use Decomp%Comm for all local info
!   97.05.18   Sawyer     DecompType has moved to ParUtilitiesTypes
!   97.05.29   Sawyer     Changed 2-D arrays to 1-D
!   97.07.03   Sawyer     Reformulated documentation
!   97.07.22   Sawyer     DecompType has moved to DecompModule
!   97.12.01   Sawyer     Changed MPI_SSEND to MPI_ISEND
!   97.12.05   Sawyer     Added InComm and Root as arguments
!   97.12.05   Sawyer     Added logic to support intercommunicators
!   98.01.24   Sawyer     Removed dependence on MPI derived types TESTED
!   98.02.05   Sawyer     Removed the use of intercommunicators
!   98.03.30   Sawyer     Stats dimension corrected: Gsize*MPI_STATUS_SIZE
!   99.01.19   Sawyer     Dropped assumed-size arrays
!   00.07.07   Sawyer     Removed "1D" references
!   00.07.23   Sawyer     Implementation with shared memory arenas
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:

      INTEGER Ierror, I, J, K, L, Iam, GroupSize
      INTEGER Status( MPI_STATUS_SIZE )
      Integer, allocatable :: Reqs(:), Stats(:)
      REAL(CPP_REAL8), ALLOCATABLE    :: SendBuf(:)
!
      CPP_ENTER_PROCEDURE( "PARSCATTERREAL" )
!
      CALL MPI_COMM_RANK( InComm, Iam, Ierror )
      CALL MPI_COMM_SIZE( InComm, GroupSize, Ierror )

      allocate (Reqs(GroupSize))
      allocate (Stats(GroupSize*MPI_STATUS_SIZE))

      IF ( Iam .EQ. Root ) THEN
        ALLOCATE( SendBuf( SUM( Decomp%NumEntries ) ) )
        L = 0
        DO I = 1, GroupSize
!
! Pick out the array sections to be sent.
! This is the inverse of the operation in ParGather
!
          DO J = 1, SIZE( Decomp%HEAD(I)%StartTags )
            DO K = Decomp%HEAD(I)%StartTags(J),Decomp%HEAD(I)%EndTags(J)
              L = L+1
              SendBuf(L) = Slice(K)
            ENDDO
          ENDDO
!
! This is a non-blocking send. SendBuf cannot be immediately deallocated
!
! WARNING: F90-MPI inconsistency: make sure the indexing below always works
!
          CALL MPI_ISEND( SendBuf(L-Decomp%NumEntries(I)+1),             &
                          Decomp%NumEntries(I), CPP_MPI_REAL8,           &
                          I-1, 0, InComm, Reqs(I), Ierror )

        ENDDO
      ENDIF

!
! All receive from the root.  
!
! The local array may be larger than that specified in the decomposition
!
      CALL MPI_RECV( Local, Decomp%NumEntries(Iam+1),                    &
                     CPP_MPI_REAL8,                                      &
                     Root, 0, InComm, Status, Ierror )
!
! Experience shows that we should wait for all the non-blocking
! PEs to check in, EVEN THOUGH THE MPI_RECV HAS COMPLETED !!
!
      IF ( Iam .EQ. Root ) THEN
        CALL MPI_WAITALL( GroupSize, Reqs, Stats, Ierror )
        DEALLOCATE( SendBuf )
      ENDIF

!
! The following may be needed on some platforms to avoid an MPI bug.
!
      CALL MPI_BARRIER( InComm, Ierror )

      deallocate (Reqs)
      deallocate (Stats)

      CPP_LEAVE_PROCEDURE( "PARSCATTERREAL" )
      RETURN
!EOC
      END SUBROUTINE ParScatterReal
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParScatterReal4 --- Scatter slice to all PEs
!
! !INTERFACE:
      SUBROUTINE ParScatterReal4 ( InComm, Root, Slice, Decomp, Local )

! !USES:
      USE decompmodule, ONLY:  DecompType, Lists
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )          :: InComm       ! Communicator
      INTEGER, INTENT( IN )          :: Root         ! Root PE
      REAL(CPP_REAL4), INTENT( IN )  :: Slice(*)     ! Global Slice
      TYPE(DecompType), INTENT( IN ) :: Decomp       ! Decomp information

! !OUTPUT PARAMETERS:
      REAL(CPP_REAL4), INTENT( OUT ) :: Local(*)     ! Local Slice

! !DESCRIPTION:
!     Given a decomposition of the domain, dole out a slice 
!     (one-dimensional array) to all the constituent PEs as described
!     by the decomposition Decomp.
!
!
! !SYSTEM ROUTINES:
!     MPI_ISEND, MPI_RECV, MPI_COMM_RANK
!
! !REVISION HISTORY:
!   97.04.14   Sawyer     Creation
!   97.04.16   Sawyer     Cleaned up for walk-through
!   97.05.01   Sawyer     Use Decomp%Comm for all local info
!   97.05.18   Sawyer     DecompType has moved to ParUtilitiesTypes
!   97.05.29   Sawyer     Changed 2-D arrays to 1-D
!   97.07.03   Sawyer     Reformulated documentation
!   97.07.22   Sawyer     DecompType has moved to DecompModule
!   97.12.01   Sawyer     Changed MPI_SSEND to MPI_ISEND
!   97.12.05   Sawyer     Added InComm and Root as arguments
!   97.12.05   Sawyer     Added logic to support intercommunicators
!   98.01.24   Sawyer     Removed dependence on MPI derived types TESTED
!   98.02.05   Sawyer     Removed the use of intercommunicators
!   98.03.30   Sawyer     Stats dimension corrected: Gsize*MPI_STATUS_SIZE
!   99.01.19   Sawyer     Dropped assumed-size arrays
!   00.07.07   Sawyer     Removed "1D" references
!   00.07.23   Sawyer     Implementation with shared memory arenas
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:

      INTEGER Ierror, I, J, K, L, Iam, GroupSize
      INTEGER Status( MPI_STATUS_SIZE )
      Integer, allocatable :: Reqs(:), Stats(:)
      REAL(CPP_REAL4), ALLOCATABLE    :: SendBuf(:)
!
      CPP_ENTER_PROCEDURE( "PARSCATTERREAL4" )
!
      CALL MPI_COMM_RANK( InComm, Iam, Ierror )
      CALL MPI_COMM_SIZE( InComm, GroupSize, Ierror )

      allocate (Reqs(GroupSize))
      allocate (Stats(GroupSize*MPI_STATUS_SIZE))

      IF ( Iam .EQ. Root ) THEN
        ALLOCATE( SendBuf( SUM( Decomp%NumEntries ) ) )
        L = 0
        DO I = 1, GroupSize
!
! Pick out the array sections to be sent.
! This is the inverse of the operation in ParGather
!
          DO J = 1, SIZE( Decomp%HEAD(I)%StartTags )
            DO K = Decomp%HEAD(I)%StartTags(J),Decomp%HEAD(I)%EndTags(J)
              L = L+1
              SendBuf(L) = Slice(K)
            ENDDO
          ENDDO
!
! This is a non-blocking send. SendBuf cannot be immediately deallocated
!
! WARNING: F90-MPI inconsistency: make sure the indexing below always works
!
          CALL MPI_ISEND( SendBuf(L-Decomp%NumEntries(I)+1),             &
                          Decomp%NumEntries(I), CPP_MPI_REAL4,           &
                          I-1, 0, InComm, Reqs(I), Ierror )

        ENDDO
      ENDIF

!
! All receive from the root.  
!
! The local array may be larger than that specified in the decomposition
!
      CALL MPI_RECV( Local, Decomp%NumEntries(Iam+1),                    &
                     CPP_MPI_REAL4,                                      &
                     Root, 0, InComm, Status, Ierror )
!
! Experience shows that we should wait for all the non-blocking
! PEs to check in, EVEN THOUGH THE MPI_RECV HAS COMPLETED !!
!
      IF ( Iam .EQ. Root ) THEN
        CALL MPI_WAITALL( GroupSize, Reqs, Stats, Ierror )
        DEALLOCATE( SendBuf )
      ENDIF

!
! The following may be needed on some platforms to avoid an MPI bug.
!
      CALL MPI_BARRIER( InComm, Ierror )

      deallocate (Reqs)
      deallocate (Stats)

      CPP_LEAVE_PROCEDURE( "PARSCATTERREAL4" )
      RETURN
!EOC
      END SUBROUTINE ParScatterReal4
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParScatterInt --- Scatter slice to all PEs
!
! !INTERFACE:
      SUBROUTINE ParScatterInt ( InComm, Root, Slice, Decomp, Local )

! !USES:
      USE decompmodule, ONLY:  DecompType, Lists
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )          :: InComm       ! Communicator
      INTEGER, INTENT( IN )          :: Root         ! Root PE
      INTEGER, INTENT( IN )          :: Slice(*)     ! Global Slice
      TYPE(DecompType), INTENT( IN ) :: Decomp       ! Decomp information

! !OUTPUT PARAMETERS:
      INTEGER, INTENT( OUT )         :: Local(*)     ! Local Slice

! !DESCRIPTION:
!     Given a decomposition of the domain, dole out a slice 
!     (one-dimensional array) to all the constituent PEs as described
!     by the decomposition Decomp.
!
!
! !SYSTEM ROUTINES:
!     MPI_ISEND, MPI_RECV, MPI_COMM_RANK
!
! !REVISION HISTORY:
!   97.04.14   Sawyer     Creation
!   97.04.16   Sawyer     Cleaned up for walk-through
!   97.05.01   Sawyer     Use Decomp%Comm for all local info
!   97.05.18   Sawyer     DecompType has moved to ParUtilitiesTypes
!   97.05.29   Sawyer     Changed 2-D arrays to 1-D
!   97.07.03   Sawyer     Reformulated documentation
!   97.07.22   Sawyer     DecompType has moved to DecompModule
!   97.12.01   Sawyer     Changed MPI_SSEND to MPI_ISEND
!   97.12.05   Sawyer     Added InComm and Root as arguments
!   97.12.05   Sawyer     Added logic to support intercommunicators
!   98.01.24   Sawyer     Removed dependence on MPI derived types TESTED
!   98.02.05   Sawyer     Removed the use of intercommunicators
!   98.03.30   Sawyer     Stats dimension corrected: Gsize*MPI_STATUS_SIZE
!   99.01.19   Sawyer     Dropped assumed-size arrays
!   00.07.07   Sawyer     Removed "1D" references
!   00.07.23   Sawyer     Implementation with shared memory arenas
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:

      INTEGER Ierror, I, J, K, L, Iam, GroupSize
      INTEGER Status( MPI_STATUS_SIZE )
      Integer, allocatable :: Reqs(:), Stats(:)
      INTEGER, ALLOCATABLE    :: SendBuf(:)
!
      CPP_ENTER_PROCEDURE( "PARSCATTERINT" )
!
      CALL MPI_COMM_RANK( InComm, Iam, Ierror )
      CALL MPI_COMM_SIZE( InComm, GroupSize, Ierror )

      allocate (Reqs(GroupSize))
      allocate (Stats(GroupSize*MPI_STATUS_SIZE))

      IF ( Iam .EQ. Root ) THEN
        ALLOCATE( SendBuf( SUM( Decomp%NumEntries ) ) )
        L = 0
        DO I = 1, GroupSize
!
! Pick out the array sections to be sent.
! This is the inverse of the operation in ParGather
!
          DO J = 1, SIZE( Decomp%HEAD(I)%StartTags )
            DO K = Decomp%HEAD(I)%StartTags(J),Decomp%HEAD(I)%EndTags(J)
              L = L+1
              SendBuf(L) = Slice(K)
            ENDDO
          ENDDO
!
! This is a non-blocking send. SendBuf cannot be immediately deallocated
!
! WARNING: F90-MPI inconsistency: make sure the indexing below always works
!
          CALL MPI_ISEND( SendBuf(L-Decomp%NumEntries(I)+1),              &
                          Decomp%NumEntries(I), CPP_MPI_INTEGER,          &
                          I-1, 0, InComm, Reqs(I), Ierror )

        ENDDO
      ENDIF

!
! All receive from the root.  
!
! The local array may be larger than that specified in the decomposition
!
      CALL MPI_RECV( Local, Decomp%NumEntries(Iam+1),                     &
                     CPP_MPI_INTEGER,                                     &
                     Root, 0, InComm, Status, Ierror )
!
! Experience shows that we should wait for all the non-blocking
! PEs to check in, EVEN THOUGH THE MPI_RECV HAS COMPLETED !!
!
      IF ( Iam .EQ. Root ) THEN
        CALL MPI_WAITALL( GroupSize, Reqs, Stats, Ierror )
        DEALLOCATE( SendBuf )
      ENDIF

!
! The following may be needed on some platforms to avoid an MPI bug.
!
      CALL MPI_BARRIER( InComm, Ierror )

      deallocate (Reqs)
      deallocate (Stats)

      CPP_LEAVE_PROCEDURE( "PARSCATTERINT" )
      RETURN
!EOC
      END SUBROUTINE ParScatterInt
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParGatherReal --- Gather Slice from all PEs
!
! !INTERFACE:  
      SUBROUTINE ParGatherReal ( InComm, Root, Local, Decomp, Slice )

! !USES:
      USE decompmodule, ONLY:  DecompType, Lists
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )          :: InComm       ! Communicator
      INTEGER, INTENT( IN )          :: Root         ! Root PE
      REAL(CPP_REAL8), INTENT( IN )  :: Local(*)     ! Local Slice
      TYPE(DecompType), INTENT( IN ) :: Decomp       ! Decomp information

! !OUTPUT PARAMETERS:
      REAL(CPP_REAL8), INTENT( OUT ) :: Slice(*)     ! Global Slice

! !DESCRIPTION:
!     Given a decomposition of the domain and a local portion of the
!     total slice on each PE, gather together the portions into a
!     global slice on the root PE
!
! !SYSTEM ROUTINES:
!     MPI_ISEND, MPI_RECV, MPI_COMM_RANK
!
! !REVISION HISTORY:
!   97.04.14   Sawyer     Creation
!   97.04.16   Sawyer     Cleaned up for walk-through
!   97.05.01   Sawyer     Use Decomp%Comm for all local info
!   97.05.18   Sawyer     DecompType has moved to ParUtilitiesTypes
!   97.05.29   Sawyer     Changed 2-D arrays to 1-D
!   97.07.03   Sawyer     Reformulated documentation
!   97.07.22   Sawyer     DecompType has moved to DecompModule
!   97.12.01   Sawyer     Changed MPI_SSEND to MPI_ISEND
!   97.12.05   Sawyer     Added InComm and Root as arguments
!   97.12.05   Sawyer     Added logic to support intercommunicators
!   98.01.24   Sawyer     Removed dependence on MPI derived types TESTED
!   98.01.29   Sawyer     Corrected assertions
!   98.02.05   Sawyer     Removed the use of intercommunicators
!   98.03.31   Sawyer     Stat dimension corrected: MPI_STATUS_SIZE
!   98.04.22   Sawyer     Local no longer assumed shape: Local(*)
!   99.01.19   Sawyer     Dropped assumed-size arrays
!   00.07.07   Sawyer     Removed "1D" references
!   00.07.23   Sawyer     Implementation with shared memory arenas
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER Ierror, I, J, K, L, Iam, GroupSize, Req
      INTEGER Status( MPI_STATUS_SIZE ), Stat( MPI_STATUS_SIZE )
      REAL(CPP_REAL8), ALLOCATABLE    :: RecvBuf(:)
!
      CPP_ENTER_PROCEDURE( "PARGATHERREAL" )
!
      CALL MPI_COMM_RANK( InComm, Iam, Ierror )
      CALL MPI_COMM_SIZE( InComm, GroupSize, Ierror )
!
! All PEs send their contribution to the root
!
      CALL MPI_ISEND( Local, Decomp%NumEntries(Iam+1),                   &
                      CPP_MPI_REAL8,                                     &
                      Root, Iam+3001, InComm, Req, Ierror )

      IF ( Iam .EQ. Root ) THEN
        ALLOCATE( RecvBuf( SUM( Decomp%NumEntries ) ) )
!
! On the Root PE receive from every other PE
!
        L = 0
        DO I = 1, GroupSize
!
! This is a blocking, synchronous recv.  All the
! sends should have been posted so it should not deadlock
!
! WARNING: F90-MPI inconsistency: make sure the indexing below always works
!
          CPP_ASSERT_F90( L .LT. SIZE( RecvBuf ) )
          CALL MPI_RECV( RecvBuf(L+1), Decomp%NumEntries(I),             &
                         CPP_MPI_REAL8, I-1, I+3000, InComm,             &
                         Status, Ierror )
!
! This is the simple reverse mapping of that in ParScatter
!
          DO J = 1, SIZE( Decomp%HEAD(I)%StartTags )
            DO K = Decomp%HEAD(I)%StartTags(J),Decomp%HEAD(I)%EndTags(J)
              L = L + 1
              Slice(K) = RecvBuf(L)
#if defined(DEBUG_PARGATHERREAL)
                PRINT *, " Entry ", L, RecvBuf(L), K, SIZE(Slice)
#endif
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE( RecvBuf )
      ENDIF
      CALL MPI_WAIT( Req, Stat, Ierror )
!
! The following may be needed on some platforms to avoid an MPI bug.
!
      CALL MPI_BARRIER( InComm, Ierror )

      CPP_LEAVE_PROCEDURE( "PARGATHERREAL" )
      RETURN
!EOC
      END SUBROUTINE ParGatherReal
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParGatherReal4 --- Gather Slice from all PEs
!
! !INTERFACE:  
      SUBROUTINE ParGatherReal4 ( InComm, Root, Local, Decomp, Slice )

! !USES:
      USE decompmodule, ONLY:  DecompType, Lists
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )          :: InComm       ! Communicator
      INTEGER, INTENT( IN )          :: Root         ! Root PE
      REAL(CPP_REAL4), INTENT( IN )  :: Local(*)     ! Local Slice
      TYPE(DecompType), INTENT( IN ) :: Decomp       ! Decomp information

! !OUTPUT PARAMETERS:
      REAL(CPP_REAL4), INTENT( OUT ) :: Slice(*)     ! Global Slice

! !DESCRIPTION:
!     Given a decomposition of the domain and a local portion of the
!     total slice on each PE, gather together the portions into a
!     global slice on the root PE
!
! !SYSTEM ROUTINES:
!     MPI_ISEND, MPI_RECV, MPI_COMM_RANK
!
! !REVISION HISTORY:
!   97.04.14   Sawyer     Creation
!   97.04.16   Sawyer     Cleaned up for walk-through
!   97.05.01   Sawyer     Use Decomp%Comm for all local info
!   97.05.18   Sawyer     DecompType has moved to ParUtilitiesTypes
!   97.05.29   Sawyer     Changed 2-D arrays to 1-D
!   97.07.03   Sawyer     Reformulated documentation
!   97.07.22   Sawyer     DecompType has moved to DecompModule
!   97.12.01   Sawyer     Changed MPI_SSEND to MPI_ISEND
!   97.12.05   Sawyer     Added InComm and Root as arguments
!   97.12.05   Sawyer     Added logic to support intercommunicators
!   98.01.24   Sawyer     Removed dependence on MPI derived types TESTED
!   98.01.29   Sawyer     Corrected assertions
!   98.02.05   Sawyer     Removed the use of intercommunicators
!   98.03.31   Sawyer     Stat dimension corrected: MPI_STATUS_SIZE
!   98.04.22   Sawyer     Local no longer assumed shape: Local(*)
!   99.01.19   Sawyer     Dropped assumed-size arrays
!   00.07.07   Sawyer     Removed "1D" references
!   00.07.23   Sawyer     Implementation with shared memory arenas
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER Ierror, I, J, K, L, Iam, GroupSize, Req
      INTEGER Status( MPI_STATUS_SIZE ), Stat( MPI_STATUS_SIZE )
      REAL(CPP_REAL4), ALLOCATABLE    :: RecvBuf(:)
!
      CPP_ENTER_PROCEDURE( "PARGATHERREAL4" )
!
      CALL MPI_COMM_RANK( InComm, Iam, Ierror )
      CALL MPI_COMM_SIZE( InComm, GroupSize, Ierror )
!
! All PEs send their contribution to the root
!
      CALL MPI_ISEND( Local, Decomp%NumEntries(Iam+1),                   &
                      CPP_MPI_REAL4,                                     &
                      Root, Iam+3001, InComm, Req, Ierror )

      IF ( Iam .EQ. Root ) THEN
        ALLOCATE( RecvBuf( SUM( Decomp%NumEntries ) ) )
!
! On the Root PE receive from every other PE
!
        L = 0
        DO I = 1, GroupSize
!
! This is a blocking, synchronous recv.  All the
! sends should have been posted so it should not deadlock
!
! WARNING: F90-MPI inconsistency: make sure the indexing below always works
!
          CPP_ASSERT_F90( L .LT. SIZE( RecvBuf ) )
          CALL MPI_RECV( RecvBuf(L+1), Decomp%NumEntries(I),             &
                         CPP_MPI_REAL4, I-1, I+3000, InComm,             &
                         Status, Ierror )
!
! This is the simple reverse mapping of that in ParScatter
!
          DO J = 1, SIZE( Decomp%HEAD(I)%StartTags )
            DO K = Decomp%HEAD(I)%StartTags(J),Decomp%HEAD(I)%EndTags(J)
              L = L + 1
              Slice(K) = RecvBuf(L)
#if defined(DEBUG_PARGATHERREAL4)
                PRINT *, " Entry ", L, RecvBuf(L), K, SIZE(Slice)
#endif
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE( RecvBuf )
      ENDIF
      CALL MPI_WAIT( Req, Stat, Ierror )
!
! The following may be needed on some platforms to avoid an MPI bug.
!
      CALL MPI_BARRIER( InComm, Ierror )
      CPP_LEAVE_PROCEDURE( "PARGATHERREAL4" )
      RETURN
!EOC
      END SUBROUTINE ParGatherReal4
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParGatherInt --- Gather Slice from all PEs
!
! !INTERFACE:  
      SUBROUTINE ParGatherInt ( InComm, Root, Local, Decomp, Slice )

! !USES:
      USE decompmodule, ONLY:  DecompType, Lists
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )          :: InComm       ! Communicator
      INTEGER, INTENT( IN )          :: Root         ! Root PE
      INTEGER, INTENT( IN )          :: Local(*)     ! Local Slice
      TYPE(DecompType), INTENT( IN ) :: Decomp       ! Decomp information

! !OUTPUT PARAMETERS:
      INTEGER, INTENT( OUT )         :: Slice(*)     ! Global Slice

! !DESCRIPTION:
!     Given a decomposition of the domain and a local portion of the
!     total slice on each PE, gather together the portions into a
!     global slice on the root PE
!
! !SYSTEM ROUTINES:
!     MPI_ISEND, MPI_RECV, MPI_COMM_RANK
!
! !REVISION HISTORY:
!   97.04.14   Sawyer     Creation
!   97.04.16   Sawyer     Cleaned up for walk-through
!   97.05.01   Sawyer     Use Decomp%Comm for all local info
!   97.05.18   Sawyer     DecompType has moved to ParUtilitiesTypes
!   97.05.29   Sawyer     Changed 2-D arrays to 1-D
!   97.07.03   Sawyer     Reformulated documentation
!   97.07.22   Sawyer     DecompType has moved to DecompModule
!   97.12.01   Sawyer     Changed MPI_SSEND to MPI_ISEND
!   97.12.05   Sawyer     Added InComm and Root as arguments
!   97.12.05   Sawyer     Added logic to support intercommunicators
!   98.01.24   Sawyer     Removed dependence on MPI derived types TESTED
!   98.01.29   Sawyer     Corrected assertions
!   98.02.05   Sawyer     Removed the use of intercommunicators
!   98.03.31   Sawyer     Stat dimension corrected: MPI_STATUS_SIZE
!   98.04.22   Sawyer     Local no longer assumed shape: Local(*)
!   99.01.19   Sawyer     Dropped assumed-size arrays
!   00.07.07   Sawyer     Removed "1D" references
!   00.07.23   Sawyer     Implementation with shared memory arenas
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER Ierror, I, J, K, L, Iam, GroupSize, Req
      INTEGER Status( MPI_STATUS_SIZE ), Stat( MPI_STATUS_SIZE )
      INTEGER, ALLOCATABLE    :: RecvBuf(:)
!
      CPP_ENTER_PROCEDURE( "PARGATHERINT" )
!
      CALL MPI_COMM_RANK( InComm, Iam, Ierror )
      CALL MPI_COMM_SIZE( InComm, GroupSize, Ierror )
!
! All PEs send their contribution to the root
!
      CALL MPI_ISEND( Local, Decomp%NumEntries(Iam+1), CPP_MPI_INTEGER,       &
                      Root, Iam+3001, InComm, Req, Ierror )

      IF ( Iam .EQ. Root ) THEN
        ALLOCATE( RecvBuf( SUM( Decomp%NumEntries ) ) )
!
! On the Root PE receive from every other PE
!
        L = 0
        DO I = 1, GroupSize
!
! This is a blocking, synchronous recv.  All the
! sends should have been posted so it should not deadlock
!
! WARNING: F90-MPI inconsistency: make sure the indexing below always works
!
          CPP_ASSERT_F90( L .LT. SIZE( RecvBuf ) )
          CALL MPI_RECV( RecvBuf(L+1), Decomp%NumEntries(I),                  &
                         CPP_MPI_INTEGER, I-1, I+3000, InComm,                &
                         Status, Ierror )
!
! This is the simple reverse mapping of that in ParScatter
!
          DO J = 1, SIZE( Decomp%HEAD(I)%StartTags )
            DO K = Decomp%HEAD(I)%StartTags(J),Decomp%HEAD(I)%EndTags(J)
              L = L + 1
              Slice(K) = RecvBuf(L)
#if defined(DEBUG_PARGATHERINT)
                PRINT *, " Entry ", L, RecvBuf(L), K, SIZE(Slice)
#endif
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE( RecvBuf )
      ENDIF
      CALL MPI_WAIT( Req, Stat, Ierror )
!
! The following may be needed on some platforms to avoid an MPI bug.
!
      CALL MPI_BARRIER( InComm, Ierror )

      CPP_LEAVE_PROCEDURE( "PARGATHERINT" )
      RETURN
!EOC
      END SUBROUTINE ParGatherInt
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParBeginTransferReal --- Start an ASYNC Real Transfer
!
! !INTERFACE:
      SUBROUTINE ParBeginTransferReal(InComm, NrInPackets, NrOutPackets, &
                                      Dest, Src, InBuf, InIA,            &
                                      OutBuf, OutIA )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )       :: InComm       ! Communicator
      INTEGER, INTENT( IN )       :: NrInPackets  ! Number of in packets
      INTEGER, INTENT( IN )       :: NrOutPackets ! Number of out packets
      INTEGER, INTENT( IN )       :: Dest(:)      ! PE destinations
      INTEGER, INTENT( IN )       :: Src(:)       ! PE sources
      REAL(CPP_REAL8), INTENT(IN) :: InBuf(:)     ! Input buffer
      INTEGER, INTENT( IN )       :: InIA(:)      ! In packet counter
      INTEGER, INTENT( IN )       :: OutIA(:)     ! Out packet counter

! !OUTPUT PARAMETERS:
      REAL(CPP_REAL8), INTENT( OUT ) :: OutBuf(:)  ! Output buffer

! !DESCRIPTION: 
!
!     This routine initiates an async. transfer of an array InBuf
!     partitioned into parcels defined by the arrays InIA and Dest
!     to an output array OutBuf on another PE. InIA(1) contains 
!     the number of reals to be sent to Dest(1), InIA(2) the number 
!     of reals to be sent to Dest(2), etc.  Similarly, the array
!     OutBuf on the calling PE is partitioned into parcels by OutIA
!     and Src, with OutIA(1) the number of reals anticipated from
!     Src(1), etc.  
!
!     The default implementation reads through the contiguous array 
!     InBuf and sends the parcels to the PEs designated with an 
!     asyncronous MPI\_ISEND.  Correspondingly it posts the receives 
!     with an asynchronous MPI\_IRECV.
!
!     Wait handles InHandle(:) and OutHandle(:) are in common block.
!
! !BUGS:
!
!     It is assumed that the buffers are passed to this routine by
!     reference!!!!!!!!!!
!
!     The buffers may not be accessed until after the call to 
!     ParEndTransferReal.
!
!
! !SYSTEM ROUTINES:
!     MPI_COMM_RANK, MPI_ISEND, MPI_IRECV
!
! !REVISION HISTORY:
!   97.09.26   Sawyer     Creation
!   97.12.05   Sawyer     Renamed Comm to InComm to avoid collisions
!   98.02.26   Sawyer     Added Dest, Src and Remote to clean up code
!   98.04.16   Sawyer     Number of packets become input arguments
!   98.09.04   Sawyer     Cleaned interface: handles in common, no Remote
!   99.03.04   Sawyer     Inlined ParCalculateRemote
!   99.06.01   Sawyer     Changed pointer arrays to INTEGER*8 for SGI
!   00.08.07   Sawyer     Implementation with shared memory arenas
!   01.09.27   Sawyer     Added multiple shared buffers for USE_MLP
!
!EOP
!-----------------------------------------------------------------------
!BOC

! !LOCAL VARIABLES:
      INTEGER Iam, GroupSize, Nr, Icnt, Packet, I, Ierr

      CPP_ENTER_PROCEDURE( "PARBEGINTRANSFERREAL" )
      CPP_ASSERT_F90( NrInPackets .LE. SIZE( Dest ) )
      CPP_ASSERT_F90( NrInPackets .LE. SIZE( InIA ) )
      CPP_ASSERT_F90( NrOutPackets .LE. SIZE( Src ) )
      CPP_ASSERT_F90( NrOutPackets .LE. SIZE( OutIA ) )

!
! Increment the ongoing transfer number
      BegTrf = MOD(BegTrf,MAX_TRF) + 1

      CALL MPI_COMM_RANK( InComm, Iam, Ierr )
      CALL MPI_COMM_SIZE( InComm, GroupSize, Ierr )

!
!     MPI: Irecv over all processes
!
      Icnt = 1
      DO Packet = 1, NrOutPackets
        Nr = OutIA( Packet )
        IF ( Nr .GT. 0 ) THEN
#if defined( DEBUG_PARBEGINTRANSFERREAL )
          PRINT *, "Iam ",Iam," posts recv ",Nr," from ", Src( Packet )
#endif
!
! Receive the buffers with MPI_Irecv. Non-blocking
!
          CPP_ASSERT_F90( Icnt+Nr-1 .LE. SIZE( OutBuf ) )
          CALL MPI_IRECV( OutBuf( Icnt ), Nr,                            &
                CPP_MPI_REAL8, Src( Packet ), Src( Packet ),             &
                InComm, OutHandle(Packet,1,BegTrf), Ierr )
        ELSE
          OutHandle(Packet,1,BegTrf) = MPI_REQUEST_NULL
        END IF
        Icnt = Icnt + Nr
      END DO
!
!     MPI: Isend over all processes
!
      Icnt = 1
      CPP_ASSERT_F90( NrInPackets .LE. SIZE( Dest ) )
      CPP_ASSERT_F90( NrInPackets .LE. SIZE( InIA ) )
      DO Packet = 1, NrInPackets
        Nr = InIA( Packet )
        IF ( Nr .GT. 0 ) THEN
#if defined( DEBUG_PARBEGINTRANSFERREAL )
          PRINT *,"Iam ",Iam," posts send ",Nr," to ",Dest( Packet )
#endif
!
!     Send the individual buffers with non-blocking sends
!
          CPP_ASSERT_F90( Icnt+Nr-1 .LE. SIZE( InBuf ) )
          CALL MPI_ISEND ( InBuf( Icnt ), Nr,                            &
                CPP_MPI_REAL8, Dest( Packet ), Iam,                      &
                InComm, InHandle(Packet,1,BegTrf), Ierr )
        ELSE
          InHandle(Packet,1,BegTrf) = MPI_REQUEST_NULL
        END IF
        Icnt = Icnt + Nr
      END DO
!
!
      CPP_LEAVE_PROCEDURE( "PARBEGINTRANSFERREAL" )
      RETURN
!EOC
      END SUBROUTINE ParBeginTransferReal
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParBeginTransferPattern1D --- Start ASYNC Pattern Transfer
!
! !INTERFACE:
      SUBROUTINE ParBeginTransferPattern1D( InComm, Pattern, InBuf, OutBuf )

! !USES:
      USE mod_comm, ONLY : mp_sendirr
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )       :: InComm       ! Communicator
      TYPE (ParPatternType), INTENT( IN )  :: Pattern   ! Comm Pattern
      REAL(CPP_REAL8), INTENT( IN )        :: InBuf(*)  ! Input buffer

! !OUTPUT PARAMETERS:
      REAL(CPP_REAL8), INTENT( OUT )       :: OutBuf(*) ! Output buffer

! !DESCRIPTION: 
!
!     This routine initiates an async. transfer of an array InBuf.
!     The communication pattern indicates the indices outgoing 
!     values of InBuf and  incoming values for OutBuf.  This routine
!     is fundamentally equivalent to ParBeginTransferReal; the use 
!     of a communication pattern is largely a performance enhancement, 
!     since it eliminates the need for intermediate buffering.
!     
!     Wait handles InHandle and OutHandle are module variables
!     The buffers may not be accessed until after the call to 
!     ParEndTransferReal.  
!
! !BUGS:
!
!     It is assumed that the buffers are passed to this routine by
!     reference.
!
! !REVISION HISTORY:
!   01.02.14   Sawyer     Creation from ParBeginTransferReal
!   01.09.27   Sawyer     Added multiple shared buffers for USE_MLP
!   02.08.13   Sawyer     Now uses mod_comm unless Use_Mpi_Types
!   03.06.24   Sawyer     All complexity now in mp_sendirr
!
!EOP
!-----------------------------------------------------------------------
!BOC

! !LOCAL VARIABLES:
      CPP_ENTER_PROCEDURE( "PARBEGINTRANSFERPATTERN1D" )

      CALL mp_sendirr( InComm,Pattern%SendDesc,Pattern%RecvDesc,InBuf,OutBuf )
!
      CPP_LEAVE_PROCEDURE( "PARBEGINTRANSFERPATTERN1D" )
      RETURN
!EOC
      END SUBROUTINE ParBeginTransferPattern1D
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParBeginTransferPattern1Dint --- Start ASYNC Pattern Transfer
!
! !INTERFACE:
      SUBROUTINE ParBeginTransferPattern1Dint( InComm, Pattern, InBuf, OutBuf )

! !USES:
      USE mod_comm, ONLY : mp_sendirr_i4
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )       :: InComm       ! Communicator
      TYPE (ParPatternType), INTENT( IN )  :: Pattern   ! Comm Pattern
      INTEGER, INTENT( IN )                :: InBuf(*)  ! Input buffer

! !OUTPUT PARAMETERS:
      INTEGER, INTENT( OUT )               :: OutBuf(*) ! Output buffer

! !DESCRIPTION: 
!
!     This routine initiates an async. transfer of an array InBuf.
!     The communication pattern indicates the indices outgoing 
!     values of InBuf and  incoming values for OutBuf.  This routine
!     is fundamentally equivalent to ParBeginTransferReal; the use 
!     of a communication pattern is largely a performance enhancement, 
!     since it eliminates the need for intermediate buffering.
!     
!     Wait handles InHandle and OutHandle are module variables
!     The buffers may not be accessed until after the call to 
!     ParEndTransferReal.  
!
! !BUGS:
!
!     It is assumed that the buffers are passed to this routine by
!     reference.
!
! !REVISION HISTORY:
!   01.02.14   Sawyer     Creation from ParBeginTransferReal
!   01.09.27   Sawyer     Added multiple shared buffers for USE_MLP
!   02.08.13   Sawyer     Now uses mod_comm unless Use_Mpi_Types
!   03.06.24   Sawyer     All complexity now in mp_sendirr_i4
! 
!EOP
!-----------------------------------------------------------------------
!BOC

      CPP_ENTER_PROCEDURE( "PARBEGINTRANSFERPATTERN1DINT" )

      CALL mp_sendirr_i4( InComm,Pattern%SendDesc,Pattern%RecvDesc,InBuf,OutBuf )

      CPP_LEAVE_PROCEDURE( "PARBEGINTRANSFERPATTERN1DINT" )
      RETURN
!EOC
      END SUBROUTINE ParBeginTransferPattern1Dint
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParBeginTransferPattern2D --- Start an ASYNC Pattern Transfer
!
! !INTERFACE:
      SUBROUTINE ParBeginTransferPattern2D( InComm, Pattern, InBuf, OutBuf )

! !USES:
      USE mod_comm, ONLY : mp_sendirr
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )       :: InComm       ! Communicator
      TYPE (ParPatternType), INTENT(IN)  :: Pattern      ! Comm Pattern
      REAL(CPP_REAL8), INTENT(IN)        :: InBuf(:,:)   ! Input buffer

! !OUTPUT PARAMETERS:
      REAL(CPP_REAL8), INTENT(OUT)       :: OutBuf(:,:)  ! Output buffer

! !DESCRIPTION: 
!
!     This routine initiates an async. transfer of an array InBuf.
!     The communication pattern indicates the indices outgoing 
!     values of InBuf and  incoming values for OutBuf.  This routine
!     is fundamentally equivalent to ParBeginTransferReal; the use 
!     of a communication pattern is largely a performance enhancement, 
!     since it eliminates the need for intermediate buffering.
!
!     Wait handles InHandle and OutHandle are module variables
!     The buffers may not be accessed until after the call to 
!     ParEndTransferReal.  
!
! !REVISION HISTORY:
!   01.10.01   Sawyer     Creation from ParBeginTransferPattern
!   02.08.13   Sawyer     Now uses mod_comm unless Use_Mpi_Types
!   03.06.24   Sawyer     All complexity now in mp_sendirr
! 
!EOP
!-----------------------------------------------------------------------
!BOC

      CPP_ENTER_PROCEDURE( "PARBEGINTRANSFERPATTERN2D" )

      CALL mp_sendirr( InComm,Pattern%SendDesc,Pattern%RecvDesc,InBuf,OutBuf )

      CPP_LEAVE_PROCEDURE( "PARBEGINTRANSFERPATTERN2D" )
      RETURN
!EOC
      END SUBROUTINE ParBeginTransferPattern2D
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParBeginTransferPattern3D --- Start an ASYNC Pattern Transfer
!
! !INTERFACE:
      SUBROUTINE ParBeginTransferPattern3D( InComm, Pattern, InBuf, OutBuf )

! !USES:
      USE mod_comm, ONLY : mp_sendirr
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )       :: InComm       ! Communicator
      TYPE (ParPatternType), INTENT(IN)  :: Pattern      ! Comm Pattern
      REAL(CPP_REAL8), INTENT(IN)        :: InBuf(:,:,:) ! Input buffer

! !OUTPUT PARAMETERS:
      REAL(CPP_REAL8), INTENT(OUT)       :: OutBuf(:,:,:)! Output buffer

! !DESCRIPTION: 
!
!     This routine initiates an async. transfer of an array InBuf.
!     The communication pattern indicates the indices outgoing 
!     values of InBuf and  incoming values for OutBuf.  This routine
!     is fundamentally equivalent to ParBeginTransferReal; the use 
!     of a communication pattern is largely a performance enhancement, 
!     since it eliminates the need for intermediate buffering.
!
!     Wait handles InHandle and OutHandle are module variables
!     The buffers may not be accessed until after the call to 
!     ParEndTransferReal.  
!
! !REVISION HISTORY:
!   01.10.01   Sawyer     Creation from ParBeginTransferPattern
!   02.08.13   Sawyer     Now uses mod_comm unless Use_Mpi_Types
!   03.06.24   Sawyer     All complexity now in mp_sendirr
! 
!EOP
!-----------------------------------------------------------------------
!BOC

      CPP_ENTER_PROCEDURE( "PARBEGINTRANSFERPATTERN3D" )

      CALL mp_sendirr( InComm,Pattern%SendDesc,Pattern%RecvDesc,InBuf,OutBuf )

      CPP_LEAVE_PROCEDURE( "PARBEGINTRANSFERPATTERN3D" )
      RETURN
!EOC
      END SUBROUTINE ParBeginTransferPattern3D
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParBeginTransferPattern4D --- Start an ASYNC Pattern Transfer
!
! !INTERFACE:
      SUBROUTINE ParBeginTransferPattern4D( InComm, Pattern, InBuf, OutBuf )

! !USES:
      USE mod_comm, ONLY : mp_sendirr
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )       :: InComm       ! Communicator
      TYPE (ParPatternType), INTENT(IN)  :: Pattern        ! Comm Pattern
      REAL(CPP_REAL8), INTENT(IN)        :: InBuf(:,:,:,:) ! Input buffer

! !OUTPUT PARAMETERS:
      REAL(CPP_REAL8), INTENT(OUT)       :: OutBuf(:,:,:,:)! Output buffer

! !DESCRIPTION: 
!
!     This routine initiates an async. transfer of an array InBuf.
!     The communication pattern indicates the indices outgoing 
!     values of InBuf and  incoming values for OutBuf.  This routine
!     is fundamentally equivalent to ParBeginTransferReal; the use 
!     of a communication pattern is largely a performance enhancement, 
!     since it eliminates the need for intermediate buffering.
!
!     Wait handles InHandle and OutHandle are module variables
!     The buffers may not be accessed until after the call to 
!     ParEndTransferReal.  
!
! !REVISION HISTORY:
!   02.12.19   Sawyer     Creation from ParBeginTransferPattern
!   03.06.24   Sawyer     All complexity now in mp_sendirr
! 
!EOP
!-----------------------------------------------------------------------
!BOC

      CPP_ENTER_PROCEDURE( "PARBEGINTRANSFERPATTERN4D" )

      CALL mp_sendirr( InComm,Pattern%SendDesc,Pattern%RecvDesc,InBuf,OutBuf )

      CPP_LEAVE_PROCEDURE( "PARBEGINTRANSFERPATTERN4D" )
      RETURN
!EOC
      END SUBROUTINE ParBeginTransferPattern4D
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParEndTransferReal --- Complete an ASYNC Real Transfer
!
! !INTERFACE:
      SUBROUTINE ParEndTransferReal( InComm, NrInPackets, NrOutPackets,  &
                                     Dest, Src, InBuf, InIA,             &
                                     OutBuf, OutIA )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )       :: InComm       ! Communicator
      INTEGER, INTENT( IN )       :: NrInPackets  ! Number of in packets
      INTEGER, INTENT( IN )       :: NrOutPackets ! Number of out packets
      INTEGER, INTENT( IN )       :: Dest(:)      ! PE destinations
      INTEGER, INTENT( IN )       :: Src(:)       ! PE sources
      REAL(CPP_REAL8), INTENT(IN) :: InBuf(:)     ! Input buffer
      INTEGER, INTENT( IN )       :: InIA(:)      ! Pointer array
      INTEGER, INTENT( IN )       :: OutIA(:)     ! Pointer array

! !INPUT/OUTPUT PARAMETERS:
      REAL(CPP_REAL8), INTENT( INOUT ) :: OutBuf(:)! Output buffer

! !DESCRIPTION: 
!
!     This routine completes an async. transfer of an array
!     partitioned into parcels defined by the array InIA.  In the 
!     MPI version, neither InBuf nor OutBuf is not used since
!     that information was utilized in ParBeginTransferReal.
!
!     The link between StartTransfer and EndTransfer is made possible
!     by the InHandle and OutHandle: they reflect the status of
!     the ongoing transfer.  When this routine completes, a valid
!     and accessible copy of the OutBuf is ready for use.
!
! !BUGS:
!
!     It is assumed that the buffers are passed to this routine by
!     reference! The buffers may not be accessed until after the 
!     completion of ParEndTransferReal.  
!
!
! !SYSTEM ROUTINES:
!     MPI_COMM_RANK, MPI_ISEND, MPI_IRECV
!
! !REVISION HISTORY:
!   97.09.26   Sawyer     Creation
!   97.12.05   Sawyer     Renamed Comm to InComm to avoid collisions
!   98.02.26   Sawyer     Count through packets, not PEs
!   98.04.16   Sawyer     Number of packets become input arguments
!   98.09.04   Sawyer     Cleaned interface: handles in common
!   99.03.05   Sawyer     Support for contiguous communicators in SHMEM
!   99.04.22   Sawyer     Bug fix: replaced MPI_WAIT with MPI_WAITALL
!   99.06.03   Sawyer     Bug fix: GroupSize in SHMEM_BARRIER
!   00.07.28   Sawyer     Implemented with shared memory arenas
!   01.09.27   Sawyer     Added multiple shared buffers for USE_MLP
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      INTEGER Iam, GroupSize, J, Offset, Packet, Ierr
      INTEGER InStats(NrInPackets*MPI_STATUS_SIZE)
      INTEGER OutStats(NrOutPackets*MPI_STATUS_SIZE)

      CPP_ENTER_PROCEDURE( "PARENDTRANSFERREAL" )

!
! Increment the receiver 
      EndTrf = MOD(EndTrf,MAX_TRF)+1

      CPP_ASSERT_F90( NrInPackets .LE. MAX_PAX )
      CALL MPI_WAITALL( NrInPackets, InHandle(:,1,EndTrf), InStats, Ierr )
 
      CPP_ASSERT_F90( NrOutPackets .LE. MAX_PAX )
      CALL MPI_WAITALL( NrOutPackets, OutHandle(:,1,EndTrf), OutStats, Ierr )
!
! WS 98.09.22 : This barrier needed to synchronize.
!
      CALL MPI_BARRIER( InComm, Ierr )

      CPP_LEAVE_PROCEDURE( "PARENDTRANSFERREAL" )
      RETURN
!EOC
      END SUBROUTINE ParEndTransferReal
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParEndTransferPattern1D --- Complete ASYNC Pattern Transfer
!
! !INTERFACE:
      SUBROUTINE ParEndTransferPattern1D( InComm, Pattern, InBuf, OutBuf )

! !USES:
      USE mod_comm, ONLY : mp_recvirr
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )       :: InComm       ! Communicator
      TYPE (ParPatternType), INTENT( IN )  :: Pattern   ! Comm Pattern
      REAL(CPP_REAL8), INTENT( IN )        :: InBuf(*)  ! Input buffer

! !INPUT/OUTPUT PARAMETERS:
      REAL(CPP_REAL8), INTENT( INOUT )     :: OutBuf(*) ! Output buffer

! !DESCRIPTION: 
!
!     This routine completes an async. transfer of an array communicated
!     with a communication pattern.  
!
!     The link between StartTransfer and EndTransfer is made possible
!     by the InHandle and OutHandle: they reflect the status of
!     the ongoing transfer.  When this routine completes, a valid
!     and accessible copy of the OutBuf is ready for use.
!     The buffers may not be accessed until after the 
!     completion of ParEndTransfer.  
!
! !BUGS:
!
!     It is assumed that the buffers are passed to this routine by
!     reference.
!
! !REVISION HISTORY:
!   01.02.14   Sawyer     Creation from ParEndTransferReal
!   02.08.13   Sawyer     Now uses mod_comm unless Use_Mpi_Types
!   03.06.24   Sawyer     All complexity now in mp_recvirr
!
!EOP
!-----------------------------------------------------------------------
!BOC

      CPP_ENTER_PROCEDURE( "PARENDTRANSFERPATTERN1D" )

      CALL mp_recvirr( InComm,Pattern%SendDesc,Pattern%RecvDesc,InBuf,OutBuf )

      CPP_LEAVE_PROCEDURE( "PARENDTRANSFERPATTERN1D" )
      RETURN
!EOC
      END SUBROUTINE ParEndTransferPattern1D
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParEndTransferPattern1Dint --- Complete ASYNC Pattern Transfer
!
! !INTERFACE:
      SUBROUTINE ParEndTransferPattern1Dint( InComm, Pattern, InBuf, OutBuf )

! !USES:
      USE mod_comm, ONLY : mp_recvirr_i4
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )       :: InComm       ! Communicator
      TYPE (ParPatternType), INTENT( IN )  :: Pattern   ! Comm Pattern
      INTEGER, INTENT( IN )                :: InBuf(*)  ! Input buffer

! !INPUT/OUTPUT PARAMETERS:
      INTEGER, INTENT( INOUT )             :: OutBuf(*) ! Output buffer

! !DESCRIPTION: 
!
!     This routine completes an async. transfer of an array communicated
!     with a communication pattern.  
!
!     The link between StartTransfer and EndTransfer is made possible
!     by the InHandle and OutHandle: they reflect the status of
!     the ongoing transfer.  When this routine completes, a valid
!     and accessible copy of the OutBuf is ready for use.
!     The buffers may not be accessed until after the 
!     completion of ParEndTransfer.  
!
! !BUGS:
!
!     It is assumed that the buffers are passed to this routine by
!     reference.
!
! !REVISION HISTORY:
!   01.02.14   Sawyer     Creation from ParEndTransferReal
!   02.08.13   Sawyer     Now uses mod_comm unless Use_Mpi_Types
!   03.06.24   Sawyer     All complexity now in mp_recvirr_i4
!
!EOP
!-----------------------------------------------------------------------
!BOC

      CPP_ENTER_PROCEDURE( "PARENDTRANSFERPATTERN1DINT" )

      CALL mp_recvirr_i4( InComm,Pattern%SendDesc,Pattern%RecvDesc,InBuf,OutBuf )

      CPP_LEAVE_PROCEDURE( "PARENDTRANSFERPATTERN1DINT" )
      RETURN
!EOC
      END SUBROUTINE ParEndTransferPattern1Dint
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParEndTransferPattern2D --- Complete an ASYNC Pattern Transfer
!
! !INTERFACE:
      SUBROUTINE ParEndTransferPattern2D( InComm, Pattern, InBuf, OutBuf )

! !USES:
      USE mod_comm, ONLY : mp_recvirr
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )       :: InComm       ! Communicator
      TYPE (ParPatternType), INTENT( IN )  :: Pattern       ! Comm Pattern
      REAL(CPP_REAL8), INTENT( IN )        :: InBuf(:,:)    ! Input buffer

! !INPUT/OUTPUT PARAMETERS:
      REAL(CPP_REAL8), INTENT( INOUT )     :: OutBuf(:,:)   ! Output buffer

! !DESCRIPTION: 
!
!     This routine completes an async. transfer of an array communicated
!     with a communication pattern.  
!
!     The link between StartTransfer and EndTransfer is made possible
!     by the InHandle and OutHandle: they reflect the status of
!     the ongoing transfer.  When this routine completes, a valid
!     and accessible copy of the OutBuf is ready for use.
!     The buffers may not be accessed until after the 
!     completion of ParEndTransfer.  
!
! !BUGS:
!
!     It is assumed that the buffers are passed to this routine by
!     reference.
!
! !REVISION HISTORY:
!   01.10.01   Sawyer     Creation from ParEndTransferPattern
!   02.08.13   Sawyer     Now uses mod_comm unless Use_Mpi_Types
!   03.06.24   Sawyer     All complexity now in mp_recvirr
!
!EOP
!-----------------------------------------------------------------------
!BOC

      CPP_ENTER_PROCEDURE( "PARENDTRANSFERPATTERN2D" )

      CALL mp_recvirr( InComm,Pattern%SendDesc,Pattern%RecvDesc,InBuf(:,:),OutBuf(:,:) )

      CPP_LEAVE_PROCEDURE( "PARENDTRANSFERPATTERN2D" )
      RETURN
!EOC
      END SUBROUTINE ParEndTransferPattern2D
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParEndTransferPattern3D --- Complete an ASYNC Pattern Transfer
!
! !INTERFACE:
      SUBROUTINE ParEndTransferPattern3D( InComm, Pattern, InBuf, OutBuf )

! !USES:
      USE mod_comm, ONLY : mp_recvirr
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )       :: InComm       ! Communicator
      TYPE (ParPatternType), INTENT( IN )  :: Pattern       ! Comm Pattern
      REAL(CPP_REAL8), INTENT( IN )        :: InBuf(:,:,:)  ! Input buffer

! !INPUT/OUTPUT PARAMETERS:
      REAL(CPP_REAL8), INTENT( INOUT )     :: OutBuf(:,:,:) ! Output buffer

! !DESCRIPTION: 
!
!     This routine completes an async. transfer of an array communicated
!     with a communication pattern.  
!
!     The link between StartTransfer and EndTransfer is made possible
!     by the InHandle and OutHandle: they reflect the status of
!     the ongoing transfer.  When this routine completes, a valid
!     and accessible copy of the OutBuf is ready for use.
!     The buffers may not be accessed until after the 
!     completion of ParEndTransfer.  
!
! !BUGS:
!
!     It is assumed that the buffers are passed to this routine by
!     reference.
!
! !REVISION HISTORY:
!   01.10.01   Sawyer     Creation from ParEndTransferPattern
!   02.08.13   Sawyer     Now uses mod_comm unless Use_Mpi_Types
!   03.06.24   Sawyer     All complexity now in mp_recvirr
!
!EOP
!-----------------------------------------------------------------------
!BOC

      CPP_ENTER_PROCEDURE( "PARENDTRANSFERPATTERN3D" )

      CALL mp_recvirr( InComm,Pattern%SendDesc,Pattern%RecvDesc,InBuf(:,:,:),OutBuf(:,:,:) )

      CPP_LEAVE_PROCEDURE( "PARENDTRANSFERPATTERN3D" )
      RETURN
!EOC
      END SUBROUTINE ParEndTransferPattern3D
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParEndTransferPattern4D --- Complete an ASYNC Pattern Transfer
!
! !INTERFACE:
      SUBROUTINE ParEndTransferPattern4D( InComm, Pattern, InBuf, OutBuf )

! !USES:
      USE mod_comm, ONLY : mp_recvirr
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )       :: InComm       ! Communicator
      TYPE (ParPatternType), INTENT( IN )  :: Pattern        ! Comm Pattern
      REAL(CPP_REAL8), INTENT( IN )        :: InBuf(:,:,:,:) ! Input buffer

! !INPUT/OUTPUT PARAMETERS:
      REAL(CPP_REAL8), INTENT( INOUT )     :: OutBuf(:,:,:,:)! Output buffer

! !DESCRIPTION: 
!
!     This routine completes an async. transfer of an array communicated
!     with a communication pattern.  
!
!     The link between StartTransfer and EndTransfer is made possible
!     by the InHandle and OutHandle: they reflect the status of
!     the ongoing transfer.  When this routine completes, a valid
!     and accessible copy of the OutBuf is ready for use.
!     The buffers may not be accessed until after the 
!     completion of ParEndTransfer.  
!
! !BUGS:
!
!     It is assumed that the buffers are passed to this routine by
!     reference.
!
! !REVISION HISTORY:
!   02.12.19   Sawyer     Creation from ParEndTransferPattern
!   03.06.24   Sawyer     All complexity now in mp_recvirr
!
!EOP
!-----------------------------------------------------------------------
!BOC

      CPP_ENTER_PROCEDURE( "PARENDTRANSFERPATTERN4D" )

      CALL mp_recvirr( InComm,Pattern%SendDesc,Pattern%RecvDesc,InBuf,OutBuf )

      CPP_LEAVE_PROCEDURE( "PARENDTRANSFERPATTERN4D" )
      RETURN
!EOC
      END SUBROUTINE ParEndTransferPattern4D
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParExchangeVectorReal --- Exchange a sparse packed vector
!
! !INTERFACE:  
      SUBROUTINE ParExchangeVectorReal ( InComm, LenInVector, InVector,  &
                                         LenOutVector, OutVector )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )   :: InComm            ! Communicator
      INTEGER, INTENT( IN )   :: LenInVector( * )  ! Length on each PE
      REAL(CPP_REAL8), INTENT( IN ):: InVector( * ) ! The input buffer

! !OUTPUT PARAMETERS:
      INTEGER, INTENT( OUT )  :: LenOutVector( * ) ! Length on each PE
      REAL(CPP_REAL8), INTENT( OUT ) :: OutVector( * ) ! The output buffer

! !DESCRIPTION:
!
!     This routine exchanges vectors stored in compressed format, i.e.,
!     in so-called compressed sparse row (CSR) format, with other
!     PEs.  In essence it first exchanges the lengths with
!     MPI\_Alltoall, then the exchange of the actual vectors (can be
!     different in size) using MPI\_AlltoallV.  Since the latter is
!     inefficient, it is simulated using MPI\_Isend and MPI\_Recv.
!
! !SYSTEM ROUTINES:
!     MPI_ISEND, MPI_RECV, MPI_WAITALL, MPI_ALLTOALL
!
! !REVISION HISTORY:
!   98.03.17   Sawyer     Creation from F77 version
!   98.03.30   Sawyer     Removed assumed shape arrays due to problems
!   99.01.18   Sawyer     Added barrier for safety
!   99.03.08   Sawyer     USE_SHMEM version for CRAY only; untested
!   99.06.01   Sawyer     USE_SHMEM version revised per comments from Tom
!   00.07.28   Sawyer     Implemented with shared memory arenas
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER :: i, iscnt, ircnt, nr, pe, icnt, Nsize, Iam, Ierr
      INTEGER :: Status(MPI_STATUS_SIZE)
      Integer, allocatable :: Reqs(:), Stats(:)

      CPP_ENTER_PROCEDURE( "PAREXCHANGEVECTORREAL" )

      CALL MPI_COMM_SIZE( InComm, Nsize, Ierr )
      CALL MPI_COMM_RANK( InComm, Iam, Ierr )

      allocate (Reqs(Nsize))
      allocate (Stats(Nsize*MPI_STATUS_SIZE))

#if defined( MY_ALLTOALL )
      DO pe = 0, Nsize-1
!
! Send the individual buffers with non-blocking sends
!
        nr = LenInVector( pe + 1 )
        CALL MPI_ISEND( nr, 1, CPP_MPI_INTEGER, pe, Iam+3000,             &
                        InComm, Reqs( pe+1 ), Ierr )
      ENDDO
      DO pe = 0, Nsize - 1
!
! Receive the buffers with MPI_Recv. Now we are blocking.
!
        CALL MPI_RECV( nr, 1, CPP_MPI_INTEGER, pe, pe+3000,               &
                       InComm, Status, Ierr )
        LenOutVector(pe + 1) = nr
      ENDDO
      CALL MPI_WAITALL( Nsize, Reqs, Stats, Ierr )
#else
      CALL MPI_ALLTOALL( LenInVector, 1, CPP_MPI_INTEGER,                 &
                         LenOutVector, 1, CPP_MPI_INTEGER,                &
                         InComm, Ierr )
#endif
!
! Over all processes
!
      icnt = 1
      DO pe = 0, Nsize-1
!
! Send the individual buffers with non-blocking sends
!
        nr = LenInVector( pe + 1 )
        IF ( nr .gt. 0 ) THEN
          CALL MPI_ISEND( InVector( icnt ), nr,                           &
                          CPP_MPI_REAL8, pe, Iam+2000,                    &
                          InComm, Reqs( pe+1 ), Ierr )
        ELSE
          Reqs( pe+1 ) = MPI_REQUEST_NULL
        ENDIF
        icnt = icnt + nr
      ENDDO

!
! Over all processes
!
      icnt = 1
      DO pe = 0, Nsize - 1
!
! Receive the buffers with MPI_Recv. Now we are blocking. 
!
        nr = LenOutVector(pe + 1)
        IF ( nr .gt. 0 ) THEN
          CALL MPI_RECV( OutVector( icnt ), nr,                          &
                         CPP_MPI_REAL8, pe, pe+2000,                     &
                         InComm, Status, Ierr )
        ENDIF
        icnt = icnt + nr
      ENDDO
      CALL MPI_WAITALL( Nsize, Reqs, Stats, Ierr )

      deallocate (Reqs)
      deallocate (Stats)

      CPP_LEAVE_PROCEDURE( "PAREXCHANGEVECTORREAL" )

      RETURN
!EOC
      END SUBROUTINE ParExchangeVectorReal
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParExchangeVectorReal4 --- Exchange a sparse packed vector
!
! !INTERFACE:  
      SUBROUTINE ParExchangeVectorReal4 ( InComm, LenInVector, InVector,&
                                          LenOutVector, OutVector )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )   :: InComm            ! Communicator
      INTEGER, INTENT( IN )   :: LenInVector( * )  ! Length on each PE
      REAL(CPP_REAL4), INTENT( IN ):: InVector( * ) ! The input buffer

! !OUTPUT PARAMETERS:
      INTEGER, INTENT( OUT )  :: LenOutVector( * ) ! Length on each PE
      REAL(CPP_REAL4), INTENT( OUT ) :: OutVector( * ) ! The output buffer

! !DESCRIPTION:
!
!     This routine exchanges vectors stored in compressed format, i.e.,
!     in so-called compressed sparse row (CSR) format, with other
!     PEs.  In essence it first exchanges the lengths with
!     MPI\_Alltoall, then the exchange of the actual vectors (can be
!     different in size) using MPI\_AlltoallV.  Since the latter is
!     inefficient, it is simulated using MPI\_Isend and MPI\_Recv.
!
! !SYSTEM ROUTINES:
!     MPI_ISEND, MPI_RECV, MPI_WAITALL, MPI_ALLTOALL
!
! !REVISION HISTORY:
!   98.03.17   Sawyer     Creation from F77 version
!   98.03.30   Sawyer     Removed assumed shape arrays due to problems
!   99.01.18   Sawyer     Added barrier for safety
!   99.03.08   Sawyer     USE_SHMEM version for CRAY only; untested
!   99.06.01   Sawyer     USE_SHMEM version revised per comments from Tom
!   00.07.28   Sawyer     Implemented with shared memory arenas
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER :: i, iscnt, ircnt, nr, pe, icnt, Nsize, Iam, Ierr
      INTEGER :: Status(MPI_STATUS_SIZE)
      Integer, allocatable :: Reqs(:), Stats(:)

      CPP_ENTER_PROCEDURE( "PAREXCHANGEVECTORREAL4" )

      CALL MPI_COMM_SIZE( InComm, Nsize, Ierr )
      CALL MPI_COMM_RANK( InComm, Iam, Ierr )

      allocate (Reqs(Nsize))
      allocate (Stats(Nsize*MPI_STATUS_SIZE))

#if defined( MY_ALLTOALL )
      DO pe = 0, Nsize-1
!
! Send the individual buffers with non-blocking sends
!
        nr = LenInVector( pe + 1 )
        CALL MPI_ISEND( nr, 1, CPP_MPI_INTEGER, pe, Iam+3000,             &
                        InComm, Reqs( pe+1 ), Ierr )
      ENDDO
      DO pe = 0, Nsize - 1
!
! Receive the buffers with MPI_Recv. Now we are blocking.
!
        CALL MPI_RECV( nr, 1, CPP_MPI_INTEGER, pe, pe+3000,               &
                       InComm, Status, Ierr )
        LenOutVector(pe + 1) = nr
      ENDDO
      CALL MPI_WAITALL( Nsize, Reqs, Stats, Ierr )
#else
      CALL MPI_ALLTOALL( LenInVector, 1, CPP_MPI_INTEGER,                 &
                         LenOutVector, 1, CPP_MPI_INTEGER,                &
                         InComm, Ierr )
#endif
!
! Over all processes
!
      icnt = 1
      DO pe = 0, Nsize-1
!
! Send the individual buffers with non-blocking sends
!
        nr = LenInVector( pe + 1 )
        IF ( nr .gt. 0 ) THEN
          CALL MPI_ISEND( InVector( icnt ), nr,                           &
                          CPP_MPI_REAL4, pe, Iam+2000,                    &
                          InComm, Reqs( pe+1 ), Ierr )
        ELSE
          Reqs( pe+1 ) = MPI_REQUEST_NULL
        ENDIF
        icnt = icnt + nr
      ENDDO

!
! Over all processes
!
      icnt = 1
      DO pe = 0, Nsize - 1
!
! Receive the buffers with MPI_Recv. Now we are blocking. 
!
        nr = LenOutVector(pe + 1)
        IF ( nr .gt. 0 ) THEN
          CALL MPI_RECV( OutVector( icnt ), nr,                          &
                         CPP_MPI_REAL4, pe, pe+2000,                     &
                         InComm, Status, Ierr )
        ENDIF
        icnt = icnt + nr
      ENDDO
      CALL MPI_WAITALL( Nsize, Reqs, Stats, Ierr )

      deallocate (Reqs)
      deallocate (Stats)

      CPP_LEAVE_PROCEDURE( "PAREXCHANGEVECTORREAL4" )

      RETURN
!EOC
      END SUBROUTINE ParExchangeVectorReal4
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParExchangeVectorInt --- Exchange a sparse packed vector
!
! !INTERFACE:  
      SUBROUTINE ParExchangeVectorInt ( InComm, LenInVector, InVector,   &
                                         LenOutVector, OutVector )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )   :: InComm            ! Communicator
      INTEGER, INTENT( IN )   :: LenInVector( * )  ! Length on each PE
      INTEGER, INTENT( IN )   :: InVector( * )     ! The input buffer

! !OUTPUT PARAMETERS:
      INTEGER, INTENT( OUT )  :: LenOutVector( * ) ! Length on each PE
      INTEGER, INTENT( OUT )  :: OutVector( * )    ! The output buffer

! !DESCRIPTION:
!
!     This routine exchanges vectors stored in compressed format, i.e.,
!     in so-called compressed sparse row (CSR) format, with other
!     PEs.  In essence it first exchanges the lengths with
!     MPI\_Alltoall, then the exchange of the actual vectors (can be
!     different in size) using MPI\_AlltoallV.  Since the latter is
!     inefficient, it is simulated using MPI\_Isend and MPI\_Recv.
!
! !SYSTEM ROUTINES:
!     MPI_ISEND, MPI_RECV, MPI_WAITALL, MPI_ALLTOALL
!
! !REVISION HISTORY:
!   98.03.17   Sawyer     Creation from F77 version
!   98.03.30   Sawyer     Removed assumed shape arrays due to problems
!   99.01.18   Sawyer     Added barrier for safety
!   99.03.08   Sawyer     USE_SHMEM version for CRAY only; untested
!   99.06.01   Sawyer     USE_SHMEM version revised per comments from Tom
!   00.07.28   Sawyer     Implemented with shared memory arenas
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER :: i, iscnt, ircnt, nr, pe, icnt, Nsize, Iam, Ierr
      INTEGER :: Status(MPI_STATUS_SIZE)
      Integer, allocatable :: Reqs(:), Stats(:)

      CPP_ENTER_PROCEDURE( "PAREXCHANGEVECTORINT" )

      CALL MPI_COMM_SIZE( InComm, Nsize, Ierr )
      CALL MPI_COMM_RANK( InComm, Iam, Ierr )

      allocate (Reqs(Nsize))
      allocate (Stats(Nsize*MPI_STATUS_SIZE))

#if defined( MY_ALLTOALL )
      DO pe = 0, Nsize-1
!
! Send the individual buffers with non-blocking sends
!
        nr = LenInVector( pe + 1 )
        CALL MPI_ISEND( nr, 1,                                           &
                        MPI_INTEGER, pe, Iam+3000,                       &
                        InComm, Reqs( pe+1 ), Ierr )
      ENDDO
      DO pe = 0, Nsize - 1
!
! Receive the buffers with MPI_Recv. Now we are blocking.
!
        CALL MPI_RECV( nr, 1,                                                 &
                       MPI_INTEGER, pe, pe+3000,                              &
                       InComm, Status, Ierr )
        LenOutVector(pe + 1) = nr
      ENDDO
      CALL MPI_WAITALL( Nsize, Reqs, Stats, Ierr )
#else
      CALL MPI_ALLTOALL( LenInVector, 1, CPP_MPI_INTEGER,                     &
                         LenOutVector, 1, CPP_MPI_INTEGER,                    &
                         InComm, Ierr )
#endif
!
! Over all processes
!
      icnt = 1
      DO pe = 0, Nsize-1
!
! Send the individual buffers with non-blocking sends
!
        nr = LenInVector( pe + 1 )
        IF ( nr .gt. 0 ) THEN
          CALL MPI_ISEND( InVector( icnt ), nr,                               &
                          CPP_MPI_INTEGER, pe, Iam+2000,                      &
                          InComm, Reqs( pe+1 ), Ierr )
        ELSE
          Reqs( pe+1 ) = MPI_REQUEST_NULL
        ENDIF
        icnt = icnt + nr
      ENDDO

!
! Over all processes
!
      icnt = 1
      DO pe = 0, Nsize - 1
!
! Receive the buffers with MPI_Recv. Now we are blocking. 
!
        nr = LenOutVector(pe + 1)
        IF ( nr .gt. 0 ) THEN
          CALL MPI_RECV( OutVector( icnt ), nr,                               &
                         CPP_MPI_INTEGER, pe, pe+2000,                        &
                         InComm, Status, Ierr )
        ENDIF
        icnt = icnt + nr
      ENDDO
      CALL MPI_WAITALL( Nsize, Reqs, Stats, Ierr )
!
! WS 98.09.22 : This barrier needed to synchronize.  Why?
!
      CALL MPI_BARRIER( InComm, Ierr )

      deallocate (Reqs)
      deallocate (Stats)

      CPP_LEAVE_PROCEDURE( "PAREXCHANGEVECTORINT" )

      RETURN
!EOC
      END SUBROUTINE ParExchangeVectorInt
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !ROUTINE: ParCollectiveBarrier --- Barrier: Simplest collective op.
!
! !INTERFACE:
      SUBROUTINE ParCollectiveBarrier( InComm )

! !USES:
      IMPLICIT NONE
! !INPUT PARAMETERS:
      INTEGER, INTENT( IN ) :: InComm   ! Communicator

! !DESCRIPTION:
!
!     This routine performs a barrier only within the communicator InComm
!     
! !REVISION HISTORY:
!   00.09.10   Sawyer     Creation
!
!EOP
!---------------------------------------------------------------------
!BOC
      INTEGER Ierror

      CALL MPI_Barrier(InComm, Ierror )

      RETURN
!EOC
      END SUBROUTINE ParCollectiveBarrier
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: ParCollective0D --- Perform global Collective of a scalar
!
! !INTERFACE:
      SUBROUTINE ParCollective0D( InComm, Op, Var )

! !USES:
      IMPLICIT NONE
! !INPUT PARAMETERS:
      INTEGER, INTENT( IN ) :: InComm   ! Communicator
      INTEGER, INTENT( IN ) :: Op       ! Operation (see header)

! !INPUT/OUTPUT PARAMETERS:
      REAL(CPP_REAL8), INTENT( INOUT ) :: Var  ! partial Var in, Var out

! !DESCRIPTION:
!
!     This utility makes a collective operation over all processes in 
!     communicator InComm.  
!     
! !REVISION HISTORY:
!   00.08.07   Sawyer     Creation
!
!EOP
!---------------------------------------------------------------------
!BOC
      INTEGER Ierror
      REAL(CPP_REAL8)    Tmp

      IF ( Op .EQ. BCSTOP ) THEN
        CALL MPI_BCAST( Var, 1, CPP_MPI_REAL8, 0, InComm, Ierror )
      ELSE
        CALL MPI_ALLREDUCE( Var, Tmp, 1, CPP_MPI_REAL8,                  &
                            Op, InComm, Ierror )
        Var = Tmp
      ENDIF

      RETURN
!EOC
      END SUBROUTINE ParCollective0D
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: ParCollective1D --- Perform component-wise global Collective of a vector
!
! !INTERFACE:
      SUBROUTINE ParCollective1D( InComm, Op, Im, Var )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN ) :: InComm   ! Communicator
      INTEGER, INTENT( IN ) :: Op       ! Operation (see header)
      INTEGER, INTENT( IN ) :: Im       ! Size of 1-D array

! !INPUT/OUTPUT PARAMETERS:
      REAL(CPP_REAL8), INTENT( INOUT ) :: Var(Im) ! partial Var in, Var out

! !DESCRIPTION:
!
!     This utility makes a collective operation over all processes in 
!     communicator InComm.  
!     
! !REVISION HISTORY:
!   00.08.07   Sawyer     Creation
!
!EOP
!---------------------------------------------------------------------
!BOC
      INTEGER Ierror
      REAL(CPP_REAL8)    Tmp(Im)

      IF ( Op .EQ. BCSTOP ) THEN
        CALL MPI_BCAST( Var, Im, CPP_MPI_REAL8, 0, InComm, Ierror )
      ELSE
        CALL MPI_ALLREDUCE( Var, Tmp, Im, CPP_MPI_REAL8,                 &
                            Op, InComm, Ierror )
        Var = Tmp
      ENDIF

      RETURN
!EOC
      END SUBROUTINE ParCollective1D
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: ParCollective1DReal4 --- Perform component-wise global Collective of a vector
!
! !INTERFACE:
      SUBROUTINE ParCollective1DReal4( InComm, Op, Im, Var )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN ) :: InComm   ! Communicator
      INTEGER, INTENT( IN ) :: Op       ! Operation (see header)
      INTEGER, INTENT( IN ) :: Im       ! Size of 1-D array

! !INPUT/OUTPUT PARAMETERS:
      REAL(CPP_REAL4), INTENT( INOUT ) :: Var(Im) ! partial Var in, Var out

! !DESCRIPTION:
!
!     This utility makes a collective operation over all processes in 
!     communicator InComm.  
!     
! !REVISION HISTORY:
!   00.08.07   Sawyer     Creation
!
!EOP
!---------------------------------------------------------------------
!BOC
      INTEGER Ierror
      REAL(CPP_REAL4)    Tmp(Im)

      IF ( Op .EQ. BCSTOP ) THEN
        CALL MPI_BCAST( Var, Im, CPP_MPI_REAL4, 0, InComm, Ierror )
      ELSE
        CALL MPI_ALLREDUCE( Var, Tmp, Im, CPP_MPI_REAL4,                 &
                            Op, InComm, Ierror )
        Var = Tmp
      ENDIF

      RETURN
!EOC
      END SUBROUTINE ParCollective1DReal4
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: ParCollective2D --- Perform component-wise collective operation
!
! !INTERFACE:
      SUBROUTINE ParCollective2D( InComm, Op, Im, Jm, Var )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN ) :: InComm     ! Communicator
      INTEGER, INTENT( IN ) :: Op         ! Operation (see header)
      INTEGER, INTENT( IN ) :: Im         ! First dimension of 2-D array
      INTEGER, INTENT( IN ) :: Jm         ! Second dimension of 2-D array

! !INPUT/OUTPUT PARAMETERS:
      REAL(CPP_REAL8), INTENT( INOUT ) :: Var(Im,Jm) ! partial Var in, Var out

! !DESCRIPTION:
!
!     This utility makes a collective operation over all processes in 
!     communicator InComm.  
!     
! !REVISION HISTORY:
!   00.08.07   Sawyer     Creation
!
!EOP
!---------------------------------------------------------------------
!BOC
      INTEGER Ierror
      REAL(CPP_REAL8)    Tmp(Im,Jm)

      IF ( Op .EQ. BCSTOP ) THEN
        CALL MPI_BCAST( Var, Im*Jm, CPP_MPI_REAL8, 0, InComm, Ierror )
      ELSE
        CALL MPI_ALLREDUCE( Var, Tmp, Im*Jm, CPP_MPI_REAL8,              &
                            Op, InComm, Ierror )
        Var = Tmp
      ENDIF

      RETURN
!EOC
      END SUBROUTINE ParCollective2D
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: ParCollective2DReal4 --- Perform component-wise collective operation
!
! !INTERFACE:
      SUBROUTINE ParCollective2DReal4( InComm, Op, Im, Jm, Var )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN ) :: InComm     ! Communicator
      INTEGER, INTENT( IN ) :: Op         ! Operation (see header)
      INTEGER, INTENT( IN ) :: Im         ! First dimension of 2-D array
      INTEGER, INTENT( IN ) :: Jm         ! Second dimension of 2-D array

! !INPUT/OUTPUT PARAMETERS:
      REAL(CPP_REAL4), INTENT( INOUT ) :: Var(Im,Jm) ! partial Var in, Var out

! !DESCRIPTION:
!
!     This utility makes a collective operation over all processes in 
!     communicator InComm.  
!     
! !REVISION HISTORY:
!   00.08.07   Sawyer     Creation
!
!EOP
!---------------------------------------------------------------------
!BOC
      INTEGER Ierror
      REAL(CPP_REAL4)    Tmp(Im,Jm)

      IF ( Op .EQ. BCSTOP ) THEN
        CALL MPI_BCAST( Var, Im*Jm, CPP_MPI_REAL4, 0, InComm, Ierror )
      ELSE
        CALL MPI_ALLREDUCE( Var, Tmp, Im*Jm, CPP_MPI_REAL4,              &
                            Op, InComm, Ierror )
        Var = Tmp
      ENDIF

      RETURN
!EOC
      END SUBROUTINE ParCollective2DReal4
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: ParCollective3D --- Perform component-wise global Collective of a vector
!
! !INTERFACE:
      SUBROUTINE ParCollective3D( InComm, Op, Im, Jm, Lm, Var )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN ) :: InComm     ! Communicator
      INTEGER, INTENT( IN ) :: Op         ! Operation (see header)
      INTEGER, INTENT( IN ) :: Im         ! First dimension of 3-D array
      INTEGER, INTENT( IN ) :: Jm         ! Second dimension of 3-D array
      INTEGER, INTENT( IN ) :: Lm         ! Third dimension of 3-D array

! !INPUT/OUTPUT PARAMETERS:
      REAL(CPP_REAL8), INTENT( INOUT ):: Var(Im,Jm,LM) ! partial Var in, Var out

! !DESCRIPTION:
!
!     This utility makes a collective operation over all processes in 
!     communicator InComm.  
!     
! !REVISION HISTORY:
!   00.08.07   Sawyer     Creation
!
!EOP
!---------------------------------------------------------------------
!BOC
      INTEGER Ierror
      REAL(CPP_REAL8) Tmp(Im,Jm,Lm)

      IF ( Op .EQ. BCSTOP ) THEN
        CALL MPI_BCAST( Var, Im*Jm*Lm, CPP_MPI_REAL8, 0, InComm, Ierror )
      ELSE
        CALL MPI_ALLREDUCE( Var, Tmp, Im*Jm*Lm, CPP_MPI_REAL8,           &
                            Op, InComm, Ierror )
        Var = Tmp
      ENDIF

      RETURN
!EOC
      END SUBROUTINE ParCollective3D
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: ParCollective0DInt --- Perform global Collective of a scalar
!
! !INTERFACE:
      SUBROUTINE ParCollective0DInt( InComm, Op, Var )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN ) :: InComm   ! Communicator
      INTEGER, INTENT( IN ) :: Op       ! Operation (see header)

! !INPUT/OUTPUT PARAMETERS:
      INTEGER, INTENT( INOUT ) :: Var   ! partial Var in, Var out

! !DESCRIPTION:
!
!     This utility makes a collective operation over all processes in 
!     communicator InComm.  
!     
! !REVISION HISTORY:
!   00.08.07   Sawyer     Creation
!
!EOP
!---------------------------------------------------------------------
!BOC
      INTEGER Ierror
      INTEGER    Tmp

      IF ( Op .EQ. BCSTOP ) THEN
        CALL MPI_BCAST( Var, 1, CPP_MPI_INTEGER, 0, InComm, Ierror )
      ELSE
        CALL MPI_ALLREDUCE( Var,Tmp,1,CPP_MPI_INTEGER,Op,InComm,Ierror )
        Var = Tmp
      ENDIF

      RETURN
!EOC
      END SUBROUTINE ParCollective0DInt
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: ParCollective0DStr --- Perform global Collective of a string
!
! !INTERFACE:
      SUBROUTINE ParCollective0DStr( InComm, Op, Var )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN ) :: InComm   ! Communicator
      INTEGER, INTENT( IN ) :: Op       ! Operation (see header)

! !INPUT/OUTPUT PARAMETERS:
      CHARACTER (LEN=*), INTENT( INOUT ) :: Var ! partial Var in, Var out

! !DESCRIPTION:
!
!     This utility makes a collective operation over all processes in
!     communicator InComm.
!
! !REVISION HISTORY:
!   00.08.07   Sawyer     Creation
!
!EOP
!---------------------------------------------------------------------
!BOC
      INTEGER Ierror, StrLen

      StrLen = LEN(Var)
      IF ( Op .EQ. BCSTOP ) THEN
        CALL MPI_BCAST( Var, StrLen, MPI_CHARACTER, 0, InComm, Ierror )
      ELSE
        write(iulog,*) "global reduction of string not supported"
      ENDIF

      RETURN
!EOC
      END SUBROUTINE ParCollective0DStr
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: ParCollective1DStr --- Perform global Collective of a string
!
! !INTERFACE:
      SUBROUTINE ParCollective1DStr( InComm, Op, Im, Var )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN ) :: InComm   ! Communicator
      INTEGER, INTENT( IN ) :: Op       ! Operation (see header)
      INTEGER, INTENT( IN ) :: Im       ! Size of 1-D array

! !INPUT/OUTPUT PARAMETERS:
      CHARACTER (LEN=*), INTENT( INOUT ) :: Var(:)   ! partial Var in, Var out

! !DESCRIPTION:
!
!     This utility makes a collective operation over all processes in 
!     communicator InComm.  
!     
! !REVISION HISTORY:
!   00.08.07   Sawyer     Creation
!
!EOP
!---------------------------------------------------------------------
!BOC
      INTEGER Ierror, StrLen

      StrLen = LEN(Var(1))
      IF ( Op .EQ. BCSTOP ) THEN
        CALL MPI_BCAST( Var, Im*StrLen, MPI_CHARACTER, 0, InComm, Ierror )
      ELSE
        write(iulog,*) "global reduction of string not supported"
      ENDIF

      RETURN
!EOC
      END SUBROUTINE ParCollective1DStr
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: ParCollective1DInt --- Perform component-wise global 
!                                  collective operations of int vector
!
! !INTERFACE:
      SUBROUTINE ParCollective1DInt( InComm, Op, Im, Var )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN ) :: InComm   ! Communicator
      INTEGER, INTENT( IN ) :: Op       ! Operation (see header)
      INTEGER, INTENT( IN ) :: Im       ! Size of 1-D array

! !INPUT/OUTPUT PARAMETERS:
      INTEGER, INTENT( INOUT ) :: Var(Im) ! partial Var in, Var out

! !DESCRIPTION:
!
!     This utility makes a collective operation over all processes in 
!     communicator InComm.  
!     
! !REVISION HISTORY:
!   00.08.07   Sawyer     Creation
!
!EOP
!---------------------------------------------------------------------
!BOC
      INTEGER Ierror
      INTEGER Tmp(Im)

      IF ( Op .EQ. BCSTOP ) THEN
        CALL MPI_BCAST( Var, Im, CPP_MPI_INTEGER, 0, InComm, Ierror )
      ELSE
        CALL MPI_ALLREDUCE( Var,Tmp,Im,CPP_MPI_INTEGER,Op,InComm,Ierror )
        Var = Tmp
      ENDIF

      RETURN
!EOC
      END SUBROUTINE ParCollective1DInt
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: ParCollective2DInt --- Perform component-wise collective op.
!
! !INTERFACE:
      SUBROUTINE ParCollective2DInt( InComm, Op, Im, Jm, Var )

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN ) :: InComm     ! Communicator
      INTEGER, INTENT( IN ) :: Op         ! Operation (see header)
      INTEGER, INTENT( IN ) :: Im         ! First dimension of 2-D array
      INTEGER, INTENT( IN ) :: Jm         ! Second dimension of 2-D array

! !INPUT/OUTPUT PARAMETERS:
      INTEGER, INTENT( INOUT ):: Var(Im,Jm) ! partial Var in, Var out

! !DESCRIPTION:
!
!     This utility makes a collective operation over all processes in 
!     communicator InComm.  
!     
! !REVISION HISTORY:
!   00.08.07   Sawyer     Creation
!
!EOP
!---------------------------------------------------------------------
!BOC
      INTEGER Ierror
      INTEGER Tmp(Im,Jm)

      IF ( Op .EQ. BCSTOP ) THEN
        CALL MPI_BCAST( Var, Im*Jm, CPP_MPI_INTEGER, 0, InComm, Ierror )
      ELSE
        CALL MPI_ALLREDUCE( Var, Tmp, Im*Jm, CPP_MPI_INTEGER,           &
                            Op, InComm, Ierror )
        Var = Tmp
      ENDIF

      RETURN
!EOC
      END SUBROUTINE ParCollective2DInt
!-----------------------------------------------------------------------
# ifdef _SMEMORY
!-----------------------------------------------------------------------
!BOP
! !IROUTINE: ParExchangeLength --- Exchange a sparse packed vector
!
! !INTERFACE:  
      SUBROUTINE ParExchangeLength ( InComm, LenInVector, LenOutVector)

! !USES:
      IMPLICIT NONE

! !INPUT PARAMETERS:
      INTEGER, INTENT( IN )   :: InComm            ! Communicator
      INTEGER, INTENT( IN )   :: LenInVector( * )  ! Length on each PE

! !OUTPUT PARAMETERS:
      INTEGER, INTENT( OUT )  :: LenOutVector( * ) ! Length on each PE

! !DESCRIPTION:
!
!     This routine exchanges vectors stored in compressed format, i.e.,
!     in so-called compressed sparse row (CSR) format, with other
!     PEs.  In essence it first exchanges the lengths with
!     MPI\_Alltoall, then the exchange of the actual vectors (can be
!     different in size) using MPI\_AlltoallV.  Since the latter is
!     inefficient, it is simulated using MPI\_Isend and MPI\_Recv.
!
! !SYSTEM ROUTINES:
!     MPI_ISEND, MPI_RECV, MPI_WAITALL, MPI_ALLTOALL
!
! !REVISION HISTORY:
!   98.03.17   Sawyer     Creation from F77 version
!   98.03.30   Sawyer     Removed assumed shape arrays due to problems
!   99.01.18   Sawyer     Added barrier for safety
!   99.03.08   Sawyer     USE_SHMEM version for CRAY only; untested
!   99.06.01   Sawyer     USE_SHMEM version revised per comments from Tom
!   00.07.28   Sawyer     Implemented with shared memory arenas
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER :: i, iscnt, ircnt, nr, pe, icnt, Nsize, Iam, Ierr
      INTEGER :: Status(MPI_STATUS_SIZE)
      Integer, allocatable :: Reqs(:), Stats(:)

      CPP_ENTER_PROCEDURE( "PAREXCHANGELENGTH" )

      CALL MPI_COMM_SIZE( InComm, Nsize, Ierr )
      CALL MPI_COMM_RANK( InComm, Iam, Ierr )

      allocate (Reqs(Nsize))
      allocate (Stats(Nsize*MPI_STATUS_SIZE))

#if defined( MY_ALLTOALL )
      DO pe = 0, Nsize-1
!
! Send the individual buffers with non-blocking sends
!
        nr = LenInVector( pe + 1 )
        CALL MPI_ISEND( nr, 1,                                           &
                        MPI_INTEGER, pe, Iam+3000,                       &
                        InComm, Reqs( pe+1 ), Ierr )
      ENDDO
      DO pe = 0, Nsize - 1
!
! Receive the buffers with MPI_Recv. Now we are blocking.
!
        CALL MPI_RECV( nr, 1,                                                 &
                       MPI_INTEGER, pe, pe+3000,                              &
                       InComm, Status, Ierr )
        LenOutVector(pe + 1) = nr
      ENDDO
      CALL MPI_WAITALL( Nsize, Reqs, Stats, Ierr )

      deallocate (Reqs)
      deallocate (Stats)

#else
      CALL MPI_ALLTOALL( LenInVector, 1, CPP_MPI_INTEGER,                     &
                         LenOutVector, 1, CPP_MPI_INTEGER,                    &
                         InComm, Ierr )
#endif
      CALL MPI_BARRIER( InComm, Ierr )


      CPP_LEAVE_PROCEDURE( "PAREXCHANGELENGTH" )

      RETURN
!EOC
      END SUBROUTINE ParExchangeLength
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:   ParCalcInfoDecompToGhost --- calculates info about the pattern
!
! !INTERFACE:
      subroutine ParCalcInfoDecompToGhost(InComm, DA,GB,Info)
!
! !USES:
      USE decompmodule, ONLY : DecompType,DecompInfo,DecompGlobalToLocal
      USE ghostmodule, ONLY : GhostType,GhostInfo
      IMPLICIT NONE

! !INPUT PARAMETERS:
      integer, intent(in)           :: InComm ! communicator
      type(DecompType), intent(in)  :: DA   ! Source Decomp Desc
      type(GhostType) , intent(in)  :: GB   ! Destination Ghost Desc

! !OUTPUT PARAMETERS:
      type (ParInfoType), intent(out) :: Info  ! Info structure
!
! !DESCRIPTION:
!     This routine calulcates the information about a communication 
!     pattern that transforms from one decomposition to another, 
!     i.e., a so-called "transpose".  This is a copy of an algorithm 
!     from the ParPatternDecompToGhost subroutine.
!
! !SYSTEM ROUTINES:
!    MPI_COMM_SIZE, MPI_COMM_RANK, MPI_ALLREDUCE
!
! !REVISION HISTORY:
!   07.09.04   Dennis     Creation based on algorithm in ParPatternDecompToGhost
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer :: nTags,oldpe,oldlocal,sTag,eTag,nCount
      integer :: j,pe,local,tag,ierr,iam,npes
      integer :: npesA,npesB,tmpA,tmp1B,tmp2B,tmp3B
      integer, allocatable :: sCount(:),rCount(:)
  
      call DecompInfo(DA,npesA,tmpA)
      call GhostInfo(GB,npesB,tmp1B,tmp2B,tmp3B)

      call MPI_COMM_SIZE(InComm,npes,ierr)
      call MPI_COMM_RANK(InComm,iam,ierr)

      allocate(sCount(npes),rCount(npes))
      sCount=0
      rCount=0
      if(iam .lt. npesB)  then 
! Parse through all the tags in the local segment
        nTags = SIZE(GB%Local%Head(iam+1)%StartTags)
        do j=1,nTags
          oldpe = -1
          oldlocal = 0 
          sTag = GB%Local%Head(iam+1)%StartTags(j)
          eTag = GB%Local%Head(iam+1)%EndTags(j)
          do tag = sTag,eTag
            if(tag > 0) then 
!
! Determine the index and PE of this entry on A. This might be inlined later
!
              call DecompGlobalToLocal(DA,tag,Local,Pe)
!
! If ipe-1 is my id, then this is an entry ipe will receive from Pe
!
              if( pe /= oldpe .or. local /= oldlocal+1) then 
                sCount(pe+1) = sCount(pe+1) + 1
              endif
              oldpe = pe  ! Update PE
              oldlocal = local  ! Update local index
            endif
          enddo
        enddo
      endif
      
! Calculate the length of receive segments
      call ParExchangeLength(InComm,sCount,rCount)
!  Record some information 
      Info%numSendSeg   = SUM(sCount)
      InFo%numSendNeigh = COUNT(sCount > 0) 

      Info%numRecvSeg   = SUM(rCount)
      InFo%numRecvNeigh = COUNT(rCount > 0) 
      nCount=MAX(Info%numSendSeg,Info%numRecvSeg)
      call MPI_ALLREDUCE(nCount,Info%maxNumSeg,1,INT4,MPI_MAX,InComm,ierr)

      deallocate(sCount,rCount)

      CPP_LEAVE_PROCEDURE( "PARCALCLENGTHDECOMPTOGHOST" )
      RETURN
!EOC
      end subroutine ParCalcInfoDecompToGhost
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP  
! !IROUTINE:   ParCalcInfoDecompToDecomp --- calculates info about the pattern
!
! !INTERFACE:
      subroutine ParCalcInfoDecompToDecomp(InComm, DA,DB,Info)
!     
! !USES:
      USE decompmodule, ONLY : DecompType,DecompInfo,DecompGlobalToLocal
      IMPLICIT NONE

! !INPUT PARAMETERS:
      integer, intent(in)           :: InComm  ! communicator
      type(DecompType), intent(in)  :: DA      ! Source Decomp Desc
      type(DecompType), intent(in)  :: DB      ! Destination Decomp Desc

! !OUTPUT PARAMETERS:
      type (ParInfoType), intent(out) :: Info  ! Info structure
!           
! !DESCRIPTION:
!     This routine calulcates the information about a communication
!     pattern that transforms from one decomposition to another,
!     i.e., a so-called "transpose".  This is a copy of an algorithm
!     from the ParPatternDecompToDecomp subroutine.
!
! !SYSTEM ROUTINES:
!    MPI_COMM_SIZE, MPI_COMM_RANK, MPI_ALLREDUCE
!
! !REVISION HISTORY:
!   07.09.04   Dennis     Creation based on algorithm in ParPatternDecompToDecomp
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer :: nCount,npes,iam,ierr
      integer :: nTags,oldpe,oldlocal,sTag,eTag
      integer :: j,pe,local,tag,tmpA,tmpB,npesA,npesB
      integer, allocatable :: sCount(:),rCount(:)

      call DecompInfo(DA,npesA,tmpA)
      call DecompInfo(DB,npesB,tmpB)

      call MPI_COMM_SIZE(InComm,npes,ierr)
      call MPI_COMM_RANK(InComm,iam,ierr)
  
      allocate(sCount(npes),rCount(npes))
      sCount=0
      rCount=0
      if(iam .lt. npesB)  then
! Parse through all the tags in the local segment
        nTags = SIZE(DB%Head(iam+1)%StartTags)
        do j=1,nTags
          oldpe = -1
          sTag = DB%Head(iam+1)%StartTags(j)
          eTag = DB%Head(iam+1)%EndTags(j)
          do tag = sTag,eTag
!
! Determine the index and PE of this entry on A. This might be inlined later
!
            call DecompGlobalToLocal(DA,tag,Local,Pe)
!
! If ipe-1 is my id, then this is an entry ipe will receive from Pe
!
            if( pe /= oldpe ) then 
              oldpe = pe
              sCount(pe+1) = sCount(pe+1) + 1
            endif
          enddo
        enddo
      endif
! Calculate the length of recieve segments    
      call ParExchangeLength(InComm,sCount,rCount)      
!  Record some information
      Info%numSendSeg   = SUM(sCount)
      InFo%numSendNeigh = COUNT(sCount > 0)

      Info%numRecvSeg   = SUM(rCount)
      InFo%numRecvNeigh = COUNT(rCount > 0)
      nCount=MAX(Info%numSendSeg,Info%numRecvSeg)
      call MPI_ALLREDUCE(nCount,Info%maxNumSeg,1,INT4,MPI_MAX,InComm,ierr)
  
      deallocate(sCount,rCount)

      CPP_LEAVE_PROCEDURE( "PARCALCINFODECOMPTODECOMP" )
      RETURN
!EOC
      end subroutine ParCalcInfoDecompToDecomp
!--------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:   ParCalcInfoGhostToDecomp --- calculates info about the pattern
!
! !INTERFACE:
      subroutine ParCalcInfoGhostToDecomp(InComm, GA,DB,Info)
!
! !USES:
      USE decompmodule, ONLY : DecompType,DecompInfo,DecompGlobalToLocal
      USE ghostmodule, ONLY : GhostType,GhostInfo
      IMPLICIT NONE

! !INPUT PARAMETERS:
      integer, intent(in)           :: InComm  ! communicator
      type(GhostType), intent(in)   :: GA      ! Source Ghost Desc
      type(DecompType), intent(in)  :: DB      ! Destination Decomp Desc

! !OUTPUT PARAMETERS:
      type (ParInfoType), intent(out) :: Info  ! Info structure
!
! !DESCRIPTION:
!     This routine calulcates the information about a communication
!     pattern that transforms from one decomposition to another,
!     i.e., a so-called "transpose".  This is a copy of an algorithm
!     from the ParPatternGhostToDecomp subroutine.
!
! !SYSTEM ROUTINES:
!    MPI_COMM_SIZE, MPI_COMM_RANK, MPI_ALLREDUCE
!
! !REVISION HISTORY:
!   07.09.04   Dennis     Creation based on algorithm in ParPatternGhostToDecomp
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer :: nTags,oldpe,oldlocal,sTag,eTag
      integer :: npes, nCount,iam,ierr
      integer :: j,pe,local,tag,npesA,npesB,tmpB,tmp1A,tmp2A,tmp3A
      integer, allocatable :: sCount(:),rCount(:)

      call GhostInfo(GA,npesA,tmp1A,tmp2A,tmp3A)
      call DecompInfo(DB,npesB,tmpB)

      call MPI_COMM_SIZE(InComm,npes,ierr)
      call MPI_COMM_RANK(InComm,iam,ierr)
  
      allocate(sCount(npes),rCount(npes))
      sCount=0
      rCount=0
      if(iam .lt. npesB) then 
! Parse through all the tags in the local segment
        nTags = SIZE(DB%Head(iam+1)%StartTags)
        do j=1,nTags
          oldpe = -1
          oldlocal = 0
          sTag = DB%Head(iam+1)%StartTags(j)
          eTag = DB%Head(iam+1)%EndTags(j)
          do tag = sTag,eTag
!
! Determine the index and PE of this entry on A. This might be inlined later
!
            call DecompGlobalToLocal(GA%Decomp,tag,Local,Pe)
!
! If ipe-1 is my id, then this is an entry ipe will receive from Pe
!
            if ((pe /= -1) .and. ( pe /= oldpe .or. local /= OldLocal+1 )) then 
              sCount(pe+1) = sCount(pe+1) + 1
            endif
            oldpe = pe
            oldlocal = local
          enddo
        enddo
      endif
! Calculate the lenght of recieve segments
      call ParExchangeLength(InComm,sCount,rCount)
!  Record some information
      Info%numSendSeg   = SUM(sCount)
      InFo%numSendNeigh = COUNT(sCount > 0)

      Info%numRecvSeg   = SUM(rCount)
      InFo%numRecvNeigh = COUNT(rCount > 0)
      nCount=MAX(Info%numSendSeg,Info%numRecvSeg)
      call MPI_ALLREDUCE(nCount,Info%maxNumSeg,1,INT4,MPI_MAX,InComm,ierr)

      deallocate(sCount,rCount)

      CPP_LEAVE_PROCEDURE( "PARCALCLENGTHGHOSTTODECOMP" )
      RETURN
!EOC
      end subroutine ParCalcInfoGhostToDecomp
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:   ParCalcInfoGhostToGhost --- calculates info about the pattern
!
! !INTERFACE:
      subroutine ParCalcInfoGhostToGhost(InComm, GA,GB,Info)
!
! !USES:
      USE decompmodule, ONLY : DecompGlobalToLocal
      USE ghostmodule, ONLY  : GhostType,GhostInfo
      IMPLICIT NONE

! !INPUT PARAMETERS:
      integer, intent(in)           :: InComm ! communicator
      type(GhostType), intent(in)   :: GA     ! Source Ghost Desc
      type(GhostType), intent(in)   :: GB     ! Destination Ghost Desc

! !OUTPUT PARAMETERS:
      type (ParInfoType), intent(out) :: Info  ! Info structure
!
! !DESCRIPTION:
!     This routine calulcates the information about a communication
!     pattern that transforms from one decomposition to another,
!     i.e., a so-called "transpose".  This is a copy of an algorithm
!     from the ParPatternGhostToGhost subroutine.
!
! !SYSTEM ROUTINES:
!    MPI_COMM_SIZE, MPI_COMM_RANK, MPI_ALLREDUCE
!
! !REVISION HISTORY:
!   07.09.04   Dennis     Creation based on algorithm in ParPatternGhostToGhost
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer :: nTags,oldpe,oldlocal,sTag,eTag,ierr,nCount
      integer :: j,pe,local,tag,npes,iam,npesA,npesB
      integer :: tmp1A,tmp2A,tmp3A,tmp1B,tmp2B,tmp3B
      integer, allocatable :: sCount(:),rCount(:)

      call GhostInfo(GA,npesA,tmp1A,tmp2A,tmp3A)
      call GhostInfo(GB,npesB,tmp1B,tmp2B,tmp3B)

      call MPI_COMM_SIZE(InComm,npes,ierr)
      call MPI_COMM_RANK(InComm,iam,ierr)
  
      allocate(sCount(npes),rCount(npes))
      sCount=0
      rCount=0 
      if(iam .lt. npesB) then 
! Parse through all the tags in the local segment
        nTags = SIZE(GB%Local%Head(iam+1)%StartTags)
        do j=1,nTags
          oldpe = -1
          oldlocal = 0
          sTag = GB%Local%Head(iam+1)%StartTags(j)
          eTag = GB%Local%Head(iam+1)%EndTags(j)
          do tag = sTag,eTag
            if (Tag > 0 ) THEN 
!
! Determine the index and PE of this entry on A. This might be inlined later
!
              call DecompGlobalToLocal(GA%Decomp,tag,Local,Pe)
!
! If ipe-1 is my id, then this is an entry ipe will receive from Pe
!
              if( pe /= oldpe .or. local /= OldLocal+1 ) then 
                sCount(pe+1)=sCount(pe+1)+1
              endif
              oldpe = pe
              oldlocal = local
            endif
          enddo
        enddo
      endif

! Calculate the length of receive segments
      call ParExchangeLength(InComm,sCount,rCount)
!  Record some information
      Info%numSendSeg   = SUM(sCount)
      InFo%numSendNeigh = COUNT(sCount > 0)

      Info%numRecvSeg   = SUM(rCount)
      InFo%numRecvNeigh = COUNT(rCount > 0)
      nCount=MAX(Info%numSendSeg,Info%numRecvSeg)
      call MPI_ALLREDUCE(nCount,Info%maxNumSeg,1,INT4,MPI_MAX,InComm,ierr)

      deallocate(sCount,rCount)

      CPP_LEAVE_PROCEDURE( "PARCALCINFOGHOSTTOGHOST" )
      RETURN
!EOC
      end subroutine ParCalcInfoGhostToGhost
# endif
#endif
      END MODULE parutilitiesmodule
