#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#undef DEBUGPART
module schedule_mod 
  use metagraph_mod, only : MetaEdge_t
  use kinds, only : int_kind, iulog
  use schedtype_mod, only : Cycle_t, Schedule_t, schedule 
  use parallel_mod, only: iam

  implicit none 
  private 

  type, public :: GraphStats_t
     integer(kind=int_kind) :: offnode
     integer(kind=int_kind) :: onnode
     integer(kind=int_kind) :: LB
     integer(kind=int_kind) :: padding 
  end type GraphStats_t

  integer,public,parameter :: BNDRY_EXCHANGE_MESSAGE=10
  integer,private,allocatable,target  :: Global2Local(:)

  integer :: MinNelemd,MaxNelemd

  public :: genEdgeSched              ! Setup the communication schedule for the edge based boundary exchange
  public :: PrintSchedule, PrintCycle
  public :: CheckSchedule
  public :: FindBufferSlot
!  public :: MessageStats
!  public :: PrimMessageStats

contains

  !********************** GENCOMSCHED.F ******************************************
  subroutine genEdgeSched(elem, PartNumber,LSchedule,MetaVertex)
    use element_mod, only : element_t
    use metagraph_mod, only : metavertex_t
    use dimensions_mod, only : nelem, max_neigh_edges
    use gridgraph_mod, only : gridvertex_t, gridedge_t, assignment ( = )
#ifdef _MPI
    use parallel_mod, only : nComPoints, iam, mpi_status_size, rrequest, srequest, &
	status, npackpoints
#else
    use parallel_mod, only : nComPoints, iam
#endif
     use parallel_mod, only : haltmp, iam
    implicit none
    type(element_t), intent(inout)        :: elem(:)
    integer, intent(in)                :: PartNumber
    type (schedule_t), intent(inout)   :: LSchedule
    type (MetaVertex_t),intent(inout)  :: MetaVertex

    integer                       :: lengthP,lengthS,total_length,lengthp_ghost
    integer                       :: i,j,is,ir,ncycle
    integer                       :: il,ie,ig
    integer                       :: nelemd0
    integer                       :: jmd
    integer                       :: inbr
    integer                       :: nSched
    integer,allocatable           :: tmpP(:,:)
    integer,allocatable           :: tmpS(:,:)
    integer,allocatable           :: tmpP_ghost(:,:)
    integer                       :: nSend,nRecv,nedges 
    integer                       :: icycle
    integer			  :: iSched
    logical, parameter            :: VerbosePrint=.FALSE.
    logical, parameter            :: Debug=.FALSE.
    integer :: ierr
    integer :: l1,l2,l1id,l2id


    nSched=SIZE(schedule)
    ! ================================================
    ! allocate some arrays for the call to MPI_gatherv
    ! ================================================

    MinNelemd = nelem
    MaxNelemd = 0
    ! =====================================================
    ! It looks like this is only used in this routine...
    ! so no need to put it in the schedule data-structure
    ! =====================================================
    allocate(Global2Local(nelem))
    if(Debug) write(iulog,*)'genEdgeSched: point #1'
    iSched = PartNumber

    nelemd0 = MetaVertex%nmembers
    !   if(iam .eq. 1)  write(iulog,*)'genEdgeSched: Part # ',i,' has ',nelemd0, ' elements '
    !if(VerbosePrint) then
    !endif
    MaxNelemd = AMAX0(MaxNelemd,nelemd0)
    MinNelemd = AMIN0(MinNelemd,nelemd0)
    if(Debug) write(iulog,*)'genEdgeSched: point #2'

    if(Debug) write(iulog,*)'genEdgeSched: point #3'
    LSchedule%ncycles = MetaVertex%nedges
    LSchedule%nelemd  = nelemd0
    if(Debug) write(iulog,*)'genEdgeSched: point #4'

    !  Note the minus one is for the internal node
    nedges = MetaVertex%nedges
    if(2*(nedges/2) .eq. nedges) then
       nedges = nedges/2
    else
       nedges = (nedges-1)/2
    endif
    LSchedule%nSendCycles = nedges
    LSchedule%nRecvCycles = nedges
    if(Debug) write(iulog,*)'genEdgeSched: point #5'

    ! Temporary array to calculate the Buffer Slot
    allocate(tmpP(2,nedges+1))
    allocate(tmpS(2,nedges+1))
    allocate(tmpP_ghost(2,nedges+1))

    tmpP(1,:) = -1
    tmpP(2,:) = 0

    tmpS(1,:) = -1
    tmpS(2,:) = 0

    tmpP_ghost(1,:) = -1
    tmpP_ghost(2,:) = 0

    !  Allocate all the cycle structures
    allocate(LSchedule%SendCycle(nedges))
    allocate(LSchedule%RecvCycle(nedges))
    allocate(LSchedule%MoveCycle(1))

    ! Initialize the schedules...
    LSchedule%MoveCycle(1)%ptrP = 0
    LSchedule%MoveCycle(1)%ptrS = 0
    LSchedule%MoveCycle(1)%lengthP = 0
    LSchedule%MoveCycle(1)%lengthS = 0
    if(Debug) write(iulog,*)'genEdgeSched: point #6'

    !==================================================================
    !  Allocate and initalized the index translation arrays
    Global2Local = -1
    allocate(LSchedule%Local2Global(nelemd0))
    allocate(LSchedule%pIndx(max_neigh_edges*nelemd0))
    allocate(LSchedule%gIndx(max_neigh_edges*nelemd0))
    LSchedule%pIndx(:)%elemId=-1
    LSchedule%pIndx(:)%edgeId=-1
    LSchedule%gIndx(:)%elemId=-1
    LSchedule%gIndx(:)%edgeId=-1
    LSchedule%pPtr=1
    LSchedule%gPtr=1

    if(Debug) write(iulog,*)'genEdgeSched: point #7'

    do il=1,nelemd0
       ig     = MetaVertex%members(il)%number
       Global2Local(ig)=il
       LSchedule%Local2Global(il)=ig
       elem(il)%desc%putmapP=-1
       elem(il)%desc%getmapP=-1
       elem(il)%desc%putmapS=-1
       elem(il)%desc%getmapS=-1
       elem(il)%desc%putmapP_ghost=-1
       elem(il)%desc%getmapP_ghost=-1
       elem(il)%desc%reverse = .FALSE.
    enddo
    !==================================================================
    if(Debug) write(iulog,*)'genEdgeSched: point #8'



    total_length = 0
    ncycle = LSchedule%ncycles
    is=1
    ir=1
    do j=1,ncycle
       lengthP     =  MetaVertex%edges(j)%wgtP
       lengthS     =  MetaVertex%edges(j)%wgtS
       lengthP_ghost     =  MetaVertex%edges(j)%wgtP_ghost

       if((MetaVertex%edges(j)%HeadVertex == PartNumber) .AND. &
            (MetaVertex%edges(j)%TailVertex == PartNumber)) then
          inbr                            = PartNumber
          if(Debug) write(iulog,*)'genEdgeSched: point #9', iam
          LSchedule%MoveCycle%ptrP         = FindBufferSlot(inbr,lengthP,tmpP)
          LSchedule%MoveCycle%ptrS         = FindBufferSlot(inbr,lengthS,tmpS)
          LSchedule%MoveCycle%ptrP_ghost   = FindBufferSlot(inbr,lengthP_ghost,tmpP_ghost)
          call SetCycle(elem, LSchedule,LSchedule%MoveCycle(1),MetaVertex%edges(j))
          if(Debug) write(iulog,*)'genEdgeSched: point #10',iam
       else if (MetaVertex%edges(j)%TailVertex == PartNumber) then
          inbr                            = MetaVertex%edges(j)%HeadVertex
          if(Debug) write(iulog,*)'genEdgeSched: point #11', iam
          LSchedule%SendCycle(is)%ptrP     = FindBufferSlot(inbr,lengthP,tmpP)
          LSchedule%SendCycle(is)%ptrS     = FindBufferSlot(inbr,lengthS,tmpS)
          LSchedule%SendCycle(is)%ptrP_ghost= FindBufferSlot(inbr,lengthP_ghost,tmpP_ghost)
          call SetCycle(elem, LSchedule,LSchedule%SendCycle(is),MetaVertex%edges(j))
          if(Debug) write(iulog,*)'genEdgeSched: point #12',iam
          is = is+1
       else if (MetaVertex%edges(j)%HeadVertex == PartNumber) then
          inbr                            = MetaVertex%edges(j)%TailVertex
          if(Debug) write(iulog,*)'genEdgeSched: point #13',iam
          LSchedule%RecvCycle(ir)%ptrP     = FindBufferSlot(inbr,lengthP,tmpP)
          LSchedule%RecvCycle(ir)%ptrS     = FindBufferSlot(inbr,lengthS,tmpS)
          LSchedule%RecvCycle(ir)%ptrP_ghost= FindBufferSlot(inbr,lengthP_ghost,tmpP_ghost)
          call SetCycle(elem, LSchedule,LSchedule%RecvCycle(ir),MetaVertex%edges(j))
          if(Debug) write(iulog,*)'genEdgeSched: point #14',iam
          ir = ir+1
       endif
    enddo
    if(iam == 1) then 
       print *,'getMetaSchedule: tmpP: ',tmpP
    endif

    deallocate(tmpP)
    deallocate(tmpS)
    deallocate(tmpP_ghost)



    do ie=1,nelemd0
       ! compute number of neighbers for each element
       elem(ie)%desc%actual_neigh_edges=0
       do i=1,max_neigh_edges
          if (elem(ie)%desc%globalID(i)>0) then
             elem(ie)%desc%actual_neigh_edges=elem(ie)%desc%actual_neigh_edges+1
          endif
       enddo

       ! normally, we loop over max_neigh_edges, checking if there is an edge
       ! let's create a mapping so that we can loop over actual_neigh_edges
       ! sort in REVERSE global id order (so the ones with globalID=0 are last)
       do l1 = 1,max_neigh_edges-1
          do l2=l1+1,max_neigh_edges
             l1id=elem(ie)%desc%loc2buf(l1)
             l2id=elem(ie)%desc%loc2buf(l2)
             if (elem(ie)%desc%globalID(l2id) > elem(ie)%desc%globalID(l1id)) then
                ! swap index:
                l1id=elem(ie)%desc%loc2buf(l2)
                elem(ie)%desc%loc2buf(l2)=elem(ie)%desc%loc2buf(l1)
                elem(ie)%desc%loc2buf(l1)=l1id
             endif
          enddo
       enddo




       elem(ie)%vertex     = MetaVertex%members(ie)
       ig                  = MetaVertex%members(ie)%number
       elem(ie)%GlobalId   = ig
       elem(ie)%LocalId    = ie  
#if 0
       call LLInsertEdge(eroot,ig,jmd)
       !DBG write(iulog,*)'After call to LLInsertEdge in schedule: ie,ig ',ie,ig,jmd
#endif
    enddo

    deallocate(Global2Local)
    !S-JMD call CheckSchedule()
#ifdef _MPI
    !================================================================
    !     Allocate a couple of structures for bndry_exchange
    !        done here to remove it from the critical path
    !================================================================

    nComPoints=0

    nSend = nedges
    nRecv = nedges
    allocate(Rrequest(nRecv))
    allocate(Srequest(nSend))
    allocate(status(MPI_STATUS_SIZE,nRecv))

    !===============================================================
    !   Number of communication points ... to be used later to
    !    setup the size of the communication buffer for MPI_Ibsend
    !===============================================================
    do icycle=1,nSend
       nComPoints = nComPoints + LSchedule%SendCycle(icycle)%lengthP
    enddo
    nPackPoints = nComPoints + LSchedule%MoveCycle(1)%lengthP

#endif
#ifdef DEBUGPART
    call haltmp("genEdgeSched: Just testing the partitioning algorithms")
#endif
  end subroutine genEdgeSched
  subroutine CheckSchedule()
    implicit none 

    integer                             :: i,nSched,nbufferwords_1,nbufferwords_2
    type (Schedule_t), pointer          :: pSchedule

    nSched = SIZE(Schedule)

    do i=1,nSched
       pSchedule => Schedule(i)
       nbufferwords_1 = SUM(pSchedule%SendCycle(:)%lengthP)
       nbufferwords_2 = SUM(pSchedule%RecvCycle(:)%lengthP)
       if(nbufferwords_1 .ne. nbufferwords_2) then 
          write (*,100) i,nbufferwords_1, nbufferwords_2
       endif
    enddo
100 format('CheckSchedule: ERR IAM:',I3,' SIZEOF(SendBuffer):',I10,' != SIZEOF(RecvBuffer) :',I10)
110 format('CheckSchedule: ERR IAM:',I3,' LENGTH(SendBuffer):',I10,' != LENGTH(RecvBuffer) :',I10)

  end subroutine CheckSchedule
  subroutine PrintSchedule(Schedule)
    use gridgraph_mod, only : printgridedge

    implicit none
    type (Schedule_t),intent(in),target   :: Schedule(:)
    type (Schedule_t), pointer            :: pSchedule
    type (Cycle_t),pointer                :: pCycle

    integer               :: i,j,nSched

    nSched = SIZE(Schedule)

    write(*,*) '------NEW SCHEDULE FORMAT---------------------'
    do i=1,nSched
       pSchedule => Schedule(i)
       write(*,*)
       write(*,*) '----------------------------------------------'
       write(*,90) i,pSchedule%ncycles
       write(*,*) '----------------------------------------------'
       write(*,*) '-----------SEND-------------------------------'
       do j=1,pSchedule%nSendCycles
          pCycle => pSchedule%SendCycle(j)
          call PrintCycle(pCycle)
          call PrintGridEdge(pCycle%edge%members)
       enddo
       write(*,*) '-----------RECV-------------------------------'
       do j=1,pSchedule%nRecvCycles
          pCycle => pSchedule%RecvCycle(j)
          call PrintCycle(pCycle)
          call PrintGridEdge(pCycle%edge%members)
       enddo
       write(*,*) '-----------MOVE-------------------------------'
       pCycle => pSchedule%MoveCycle(1)
       call PrintCycle(pCycle)
       call PrintGridEdge(pCycle%edge%members)
    enddo
90  format('NODE # ',I2,2x,'NCYCLES ',I2)
97  format(10x,'EDGE #',I2,2x,'TYPE ',I1,2x,'G.EDGES',I4,2x,'WORDS ',I5,2x, &
         'SRC ',I3,2x,'DEST ',I3,2x,'PTR ',I4)
100 format(15x,I4,5x,I3,1x,'(',I1,') --',I1,'--> ',I3,1x,'(',I1,')')


  end subroutine PrintSchedule

  subroutine PrintCycle(Cycle)

    implicit none 
    type (Cycle_t),intent(in),target  ::  Cycle

    write(*,97) Cycle%edge%number,Cycle%type,Cycle%edge%nmembers, &
         Cycle%lengthP,Cycle%source, Cycle%dest,Cycle%ptrP

97  format(5x,'METAEDGE #',I2,2x,'TYPE ',I1,2x,'G.EDGES',I4,2x,'WORDS ',I5,2x, &
         'SRC ',I3,2x,'DEST ',I3,2x,'PTR ',I5)

  end subroutine PrintCycle

  subroutine SetCycle(elem, schedule,Cycle,Edge)
    use element_mod, only : element_t
    use dimensions_mod, only : max_corner_elem, max_neigh_edges
    use parallel_mod, only : abortmp, iam   
    implicit none 

    type(element_t), intent(inout)        :: elem(:)
    type (Schedule_t),intent(inout)          :: Schedule
    type (Cycle_t),intent(inout)             :: Cycle
    type (MetaEdge_t),intent(in),target      :: Edge
    integer                                  :: i,il,face, loc, dir


    do i=1,Edge%nmembers
       !   Setup send index
       il                     = Global2Local(Edge%members(i)%tail%number)
       face                   = Edge%members(i)%tail_face
       !need to convert the location of corner elements for getmap and putmap
       if (face.ge.5) then ! if a corner element
          dir = Edge%members(i)%tail_dir
          loc = MOD(dir,max_corner_elem) !this is the location within that direction
          dir = (dir - loc)/max_corner_elem !this is the direction (1-8)
          loc = dir + (dir-5)*(max_corner_elem-1)+loc
       else
          loc = face
       end if

       if(il .gt. 0) then 
          elem(il)%desc%putmapP(loc) = Edge%edgeptrP(i) + Cycle%ptrP - 1  ! offset, so start at 0
          elem(il)%desc%putmapS(loc) = Edge%edgeptrS(i) + Cycle%ptrS - 1
          elem(il)%desc%putmapP_ghost(loc) = Edge%edgeptrP_ghost(i) + Cycle%ptrP_ghost  ! index, start at 1
          elem(il)%desc%reverse(loc) = Edge%members(i)%reverse
          schedule%pIndx(schedule%pPtr)%elemid=il
          schedule%pIndx(schedule%pPtr)%edgeid=loc
          schedule%pPtr=schedule%pPtr+1
!NEWEDGEBUFF          if(iam == 1) then 
!NEWEDGEBUFF            print *,'putmap: elemid,loc: ',il,loc
!NEWEDGEBUFF          endif
       endif



       !   Setup receive index
       il                     = Global2Local(Edge%members(i)%head%number)
       face                   = Edge%members(i)%head_face
       !need to convert the location of corner elements for getmap and putmap
       if (face.ge.5) then !its a corner
          dir = Edge%members(i)%head_dir
          loc = MOD(dir,max_corner_elem) !this is the location within that direction
          dir = (dir - loc)/max_corner_elem !this is the direction (1-8)
          loc = dir + (dir-5)*(max_corner_elem-1)+loc
          if(loc>max_neigh_edges) then
             print *,__FILE__,__LINE__,iam,face,i,max_corner_elem,max_neigh_edges,edge%members(i)%head_face
             call abortmp('max_neigh_edges set too low.')
          end if
       else
          loc = face
       end if

       if(il .gt. 0) then 
          elem(il)%desc%getmapP(loc) = Edge%edgeptrP(i) + Cycle%ptrP - 1
          elem(il)%desc%getmapS(loc) = Edge%edgeptrS(i) + Cycle%ptrS - 1
          elem(il)%desc%getmapP_ghost(loc) = Edge%edgeptrP_ghost(i) + Cycle%ptrP_ghost 
          elem(il)%desc%globalID(loc) = Edge%members(i)%tail%number
          schedule%gIndx(schedule%gPtr)%elemid=il
          schedule%gIndx(schedule%gPtr)%edgeid=loc
          schedule%gPtr=schedule%gPtr+1
!NEWEDGEBUFF          if(iam == 1) then 
!NEWEDGEBUFF            print *,'getmap: elemid,loc: ',il,loc
!NEWEDGEBUFF          endif
       endif
    enddo
    Cycle%edge   => Edge
    Cycle%type   = Edge%type
    Cycle%dest   = Edge%HeadVertex
    Cycle%source = Edge%TailVertex
    Cycle%tag    = BNDRY_EXCHANGE_MESSAGE
    Cycle%lengthP = Edge%wgtP
    Cycle%lengthS = Edge%wgtS
    Cycle%lengthP_ghost = Edge%wgtP_ghost

  end subroutine SetCycle

  function FindBufferSlot(inbr,length,tmp) result(ptr)

    integer                          :: ptr
    integer,intent(in)               :: inbr,length
    integer,intent(inout)    :: tmp(:,:)

    integer                          :: i,n

    n = SIZE(tmp,2)

    ptr = 0
    do i=1,n
       if( tmp(1,i) == inbr) then 
          ptr = tmp(2,i)
          return	
       endif
       if( tmp(1,i) == -1 ) then  
          tmp(1,i) = inbr
          if(i .eq. 1) tmp(2,i) = 1
          ptr = tmp(2,i)
          if(i .ne. n) tmp(2,i+1) = ptr +length
          return
       endif
    enddo

  end function FindBufferSlot

end module schedule_mod
