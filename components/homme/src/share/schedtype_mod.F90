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
     type (MetaEdge_t),pointer  :: edge
  end type Cycle_t

  type, public :: pgindex_t
      integer :: elemid
      integer :: edgeid
  end type pgindex_t

  type, public :: Schedule_t
     integer                     ::  ncycles
     integer                     ::  nelemd
     integer                     :: placeholder  ! total integer count should be even
     integer                     :: nSendCycles
     integer                     :: nRecvCycles
     integer                     :: padding
     integer,pointer             :: Local2Global(:)
     type (Cycle_t), pointer     ::   Cycle(:)          
     type (Cycle_t), pointer :: SendCycle(:)
     type (Cycle_t), pointer :: RecvCycle(:)
     type (Cycle_t), pointer :: MoveCycle(:)
     type (pgindex_t), pointer :: pIndx(:)
     type (pgindex_t), pointer :: gIndx(:)
     integer :: pPtr,gPtr
  end type Schedule_t

  type (Schedule_t), public, allocatable, target  :: Schedule(:)


end module schedtype_mod
