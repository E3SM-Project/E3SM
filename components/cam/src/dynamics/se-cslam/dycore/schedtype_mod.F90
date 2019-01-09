module schedtype_mod

  use metagraph_mod, only : MetaEdge_t

  implicit none
  private
  type, public :: Cycle_t
     integer                    :: tag
     integer                    :: dest
     integer                    :: source
     integer                    :: lengthP
     integer                    :: lengthP_ghost
     integer                    :: lengthS
     integer                    :: type
     integer                    :: ptrP
     integer                    :: ptrP_ghost
     integer                    :: ptrS
     logical                    :: onNode
     type (MetaEdge_t),pointer  :: edge
  end type Cycle_t

  type, public :: pgindex_t
      integer :: elemid
      integer :: edgeid
      integer :: mesgid
      integer :: lenP,lenS
      integer :: edgeType
  end type pgindex_t

  type, public :: Schedule_t
     integer                     ::  ncycles
     integer                     ::  nelemd
     integer                     :: placeholder  ! total integer count should be even
     integer                     :: nSendCycles
     integer                     :: nRecvCycles
     integer                     :: nInter       ! number of off-node or inter node communication cycles
     integer                     :: nIntra       ! number of on-node or intra node communication cycles
     integer                     :: padding
     integer,pointer             :: Local2Global(:)
     integer,pointer             :: destFull(:)
     integer,pointer             :: srcFull(:)
     type (Cycle_t), pointer     :: Cycle(:)          
     type (Cycle_t), pointer :: SendCycle(:)
     type (Cycle_t), pointer :: RecvCycle(:)
     type (Cycle_t), pointer :: MoveCycle(:)
     type (pgindex_t), pointer :: pIndx(:)
     type (pgindex_t), pointer :: gIndx(:)
     integer :: pPtr,gPtr
  end type Schedule_t

  type (Schedule_t), public, allocatable, target  :: Schedule(:)
  type (Schedule_t), public, allocatable, target  :: gSchedule(:)
  type (Schedule_t), public, allocatable, target  :: sSchedule(:)

  integer,public,parameter :: HME_Cardinal = 101
  integer,public,parameter :: HME_Ordinal  = 102


end module schedtype_mod
