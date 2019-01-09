module schedule_mod
  use metagraph_mod, only: MetaEdge_t
  use schedtype_mod, only: Cycle_t, Schedule_t, schedule, pgindex_t, HME_Ordinal,HME_Cardinal
  use parallel_mod,  only: parallel_t
  use cam_logfile,   only: iulog

  implicit none
  private

  type, public :: GraphStats_t
     integer :: offnode
     integer :: onnode
     integer :: LB
     integer :: padding
  end type GraphStats_t

  integer,public,parameter :: HME_CYCLE_SEND=1
  integer,public,parameter :: HME_CYCLE_RECV=2
  integer,public,parameter :: HME_CYCLE_MOVE=3
  integer,public,parameter :: HME_CYCLE_ANY =4


  integer,public,parameter :: BNDRY_EXCHANGE_MESSAGE=10
  integer,private,allocatable,target  :: Global2Local(:)

  integer :: MinNelemd,MaxNelemd

  public :: genEdgeSched              ! Setup the communication schedule for the edge based boundary exchange
  public :: PrintSchedule, PrintCycle
  public :: PrintIndex
  public :: CheckSchedule
  public :: FindBufferSlot

contains

  subroutine genEdgeSched(par,elem, PartNumber,LSchedule,MetaVertex)
    use element_mod,    only: element_t
    use metagraph_mod,  only: metavertex_t
    use dimensions_mod, only: nelem, max_neigh_edges
    use gridgraph_mod,  only: gridvertex_t, gridedge_t, assignment ( = )
    use cam_abortutils, only: endrun
    use spmd_utils,     only: mpi_status_size, mpi_info_null, mpi_success
    use parallel_mod,   only: nComPoints, rrequest, srequest, status, npackpoints

    type(parallel_t),    intent(inout) :: par
    type(element_t),     intent(inout) :: elem(:)
    integer,             intent(in)    :: PartNumber
    type (schedule_t),   intent(inout) :: LSchedule
    type (MetaVertex_t), intent(inout) :: MetaVertex

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
    integer           :: iSched
    logical, parameter            :: VerbosePrint=.FALSE.
    logical, parameter            :: Debug=.FALSE.
    character(len=*),       parameter :: subname = 'genEdgeSched'
    integer                           :: errorcode,errorlen
    character*(80)                    :: errorstring
    integer, allocatable :: intracommranks(:)
    integer :: numIntra, numInter, rank
    logical :: OnNode


    integer :: ierr
    integer :: l1,l2,l1id,l2id
    integer :: src,dest,wgt
    integer :: icIntra, icInter

    integer, allocatable :: srcFull(:), destFull(:),  srcweightFull(:), destweightFull(:)
    integer, allocatable :: srcInter(:),destInter(:), srcweightInter(:),destweightInter(:) 
    integer, allocatable :: srcIntra(:),destIntra(:), srcweightIntra(:),destweightIntra(:) 

    logical :: reorder
    integer :: sizeGroup, groupFull

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


    !  Allocate all the cycle structures
    allocate(LSchedule%SendCycle(nedges))
    allocate(LSchedule%RecvCycle(nedges))
    allocate(LSchedule%MoveCycle(1))

    ! Initialize the schedules...
    LSchedule%MoveCycle(1)%ptrP = 0
    LSchedule%MoveCycle(1)%ptrS = 0
    LSchedule%MoveCycle(1)%lengthP = 0
    if(Debug) write(iulog,*)'genEdgeSched: point #6'

    !==================================================================
    !  Allocate and initalized the index translation arrays
    Global2Local = -1
    allocate(LSchedule%Local2Global(nelemd0))
    allocate(LSchedule%pIndx(max_neigh_edges*nelemd0))
    allocate(LSchedule%gIndx(max_neigh_edges*nelemd0))

    LSchedule%pIndx(:)%elemId   = -1
    LSchedule%pIndx(:)%edgeId   = -1
    LSchedule%pIndx(:)%lenP     = -1
    LSchedule%pIndx(:)%lenS     = -1
    LSchedule%pIndx(:)%mesgid   = -1
    LSchedule%pIndx(:)%edgeType = -1

    LSchedule%gIndx(:)%elemId   = -1
    LSchedule%gIndx(:)%edgeId   = -1
    LSchedule%gIndx(:)%lenP     = -1
    LSchedule%gIndx(:)%lenS     = -1
    LSchedule%gIndx(:)%mesgid   = -1
    LSchedule%gIndx(:)%edgeType = -1

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
    !
    ! Send Cycle
    !
    is=1
    tmpP(1,:) = -1
    tmpP(2,:) = 0
    tmpS(1,:) = -1
    tmpS(2,:) = 0
    tmpP_ghost(1,:) = -1
    tmpP_ghost(2,:) = 0

    do j=1,ncycle
       lengthP     =  MetaVertex%edges(j)%wgtP
       lengthS     =  MetaVertex%edges(j)%wgtS
       lengthP_ghost     =  MetaVertex%edges(j)%wgtP_ghost

       if ((MetaVertex%edges(j)%TailVertex == PartNumber) .AND. &
                (MetaVertex%edges(j)%HeadVertex .ne. PartNumber) ) then
          inbr                            = MetaVertex%edges(j)%HeadVertex
          if(Debug) write(iulog,*)'genEdgeSched: point #11', par%rank
          LSchedule%SendCycle(is)%ptrP     = FindBufferSlot(inbr,lengthP,tmpP)
          LSchedule%SendCycle(is)%ptrS     = FindBufferSlot(inbr,lengthS,tmpS)
          LSchedule%SendCycle(is)%ptrP_ghost= FindBufferSlot(inbr,lengthP_ghost,tmpP_ghost)
          call SetCycle(par, elem, LSchedule,LSchedule%SendCycle(is),MetaVertex%edges(j), HME_CYCLE_SEND)
          if(Debug) write(iulog,*)'genEdgeSched: point #12',par%rank
          is = is+1
       endif
    enddo

    !
    ! Recv Cycle:  Note that by reinitializing the tmpP array we change the structure of the receive buffer
    !
    ir=1
    tmpP(1,:) = -1
    tmpP(2,:) = 0
    tmpS(1,:) = -1
    tmpS(2,:) = 0
    tmpP_ghost(1,:) = -1
    tmpP_ghost(2,:) = 0

    do j=1,ncycle
       lengthP     =  MetaVertex%edges(j)%wgtP
       lengthS     =  MetaVertex%edges(j)%wgtS
       lengthP_ghost     =  MetaVertex%edges(j)%wgtP_ghost

       if ( (MetaVertex%edges(j)%HeadVertex == PartNumber) .AND. &
               (MetaVertex%edges(j)%TailVertex .ne. PartNumber) ) then
          inbr                            = MetaVertex%edges(j)%TailVertex
          if(Debug) write(iulog,*)'genEdgeSched: point #13',par%rank
          LSchedule%RecvCycle(ir)%ptrP     = FindBufferSlot(inbr,lengthP,tmpP)
          LSchedule%RecvCycle(ir)%ptrS     = FindBufferSlot(inbr,lengthS,tmpS)
          LSchedule%RecvCycle(ir)%ptrP_ghost= FindBufferSlot(inbr,lengthP_ghost,tmpP_ghost)
          call SetCycle(par, elem, LSchedule,LSchedule%RecvCycle(ir),MetaVertex%edges(j),HME_CYCLE_RECV)
          if(Debug) write(iulog,*)'genEdgeSched: point #14',par%rank
          ir = ir+1
       endif
    enddo

    ! Put the move cycle at the end of the buffer.
    do j=1,ncycle
       lengthP     =  MetaVertex%edges(j)%wgtP
       lengthS     =  MetaVertex%edges(j)%wgtS
       lengthP_ghost     =  MetaVertex%edges(j)%wgtP_ghost

       if((MetaVertex%edges(j)%HeadVertex == PartNumber) .AND. &
            (MetaVertex%edges(j)%TailVertex == PartNumber)) then
          inbr                            = PartNumber
          if(Debug) write(iulog,*)'genEdgeSched: point #9', par%rank
          LSchedule%MoveCycle%ptrP         = FindBufferSlot(inbr,lengthP,tmpP)
          LSchedule%MoveCycle%ptrS         = FindBufferSlot(inbr,lengthS,tmpS)
          LSchedule%MoveCycle%ptrP_ghost   = FindBufferSlot(inbr,lengthP_ghost,tmpP_ghost)
          call SetCycle(par, elem, LSchedule,LSchedule%MoveCycle(1),MetaVertex%edges(j),HME_CYCLE_MOVE)
          if(Debug) write(iulog,*)'genEdgeSched: point #10',par%rank
       endif
    enddo

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
    enddo

    deallocate(Global2Local)

#ifdef SPMD
    !================================================================
    !     Allocate a couple of structures for bndry_exchange
    !        done here to remove it from the critical path
    !================================================================
    nComPoints = 0

    nSend = nedges
    nRecv = nedges
    allocate(Rrequest(nRecv))
    allocate(Srequest(nSend))
    allocate(status(MPI_STATUS_SIZE,nRecv))

    !===============================================================
    !   Number of communication points ... to be used later to
    !    setup the size of the communication buffer for MPI_Ibsend
    !===============================================================
    do icycle = 1, nSend
      nComPoints = nComPoints + LSchedule%SendCycle(icycle)%lengthP
    end do
    nPackPoints = nComPoints + LSchedule%MoveCycle(1)%lengthP
#if MPI_VERSION >= 3
   ! Create a communicator that only contains the on-node MPI ranks
   call MPI_Comm_split_type(par%comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, par%intracomm, ierr)
   
   call MPI_Comm_size(par%intracomm, par%intracommsize, ierr)
   call MPI_Comm_rank(par%intracomm, par%intracommrank, ierr)

   allocate(intracommranks(par%intracommsize))
   call MPI_Allgather(par%rank,1,MPIinteger_t,intracommranks,1,MPIinteger_t,par%intracomm,ierr)

   numIntra=0
   do icycle=1,nSend
      rank = LSchedule%SendCycle(icycle)%dest - 1
      onNode = isIntraComm(intracommranks,rank)
      LSchedule%SendCycle(icycle)%onNode = onNode
      if(onNode) then 
         numIntra=numIntra+1
      endif
   enddo
   do icycle=1,nRecv
      rank = LSchedule%RecvCycle(icycle)%source - 1
      onNode = isIntraComm(intracommranks,rank)
      LSchedule%RecvCycle(icycle)%onNode = onNode
   enddo
   numInter = nsend-numIntra 

   
   deallocate(intracommranks)
#else
   numIntra = 0
   numInter = nSend
   ! Mark all communications as off-node by default
   do icycle=1,nSend
      LSchedule%SendCycle(icycle)%onNode = .False.
   enddo
   do icycle=1,nRecv
      LSchedule%RecvCycle(icycle)%onNode = .False.
   enddo
#endif
    LSchedule%nInter = numInter
    LSchedule%nIntra = numIntra

    allocate(srcFull(nRecv), srcWeightFull(nRecv),destFull(nSend),destWeightFull(nSend))
    if(numInter>0) then 
      allocate(srcInter(numInter),srcWeightInter(numInter),destInter(numInter), destWeightInter(numInter))
    endif
    if(numIntra>0) then 
      allocate(srcIntra(numIntra),srcWeightIntra(numIntra),destIntra(numIntra), destWeightIntra(numIntra))
    endif

    icIntra=0
    icInter=0
    do icycle=1,nSend
       dest = LSchedule%SendCycle(icycle)%dest - 1
       wgt  = LSchedule%SendCycle(icycle)%lengthP
       destFull(icycle) = dest
       destWeightFull(icycle) = wgt
       if(LSchedule%SendCycle(icycle)%onNode) then 
          icIntra=icIntra+1
          destIntra(icIntra) = dest          
          destWeightIntra(icIntra) = wgt
       else
          icInter=icInter+1
          destInter(icInter) = dest          
          destWeightInter(icInter) = wgt
       endif
    enddo

    icIntra=0
    icInter=0
    do icycle=1,nRecv
       src = LSchedule%RecvCycle(icycle)%source - 1
       wgt = LSchedule%RecvCycle(icycle)%lengthP 
       srcFull(icycle) = src
       srcWeightFUll(icycle) = wgt
       if(LSchedule%RecvCycle(icycle)%onNode) then
          icIntra=icIntra+1
          srcIntra(icIntra) = src
          srcWeightIntra(icIntra) = wgt
       else
          icInter=icInter+1
          srcInter(icInter) = src
          srcWeightInter(icInter) = wgt
       endif
    enddo

    ! construct the FULL communication graph 
    reorder=.FALSE.
    call MPI_Dist_graph_create_adjacent(par%comm, nRecv,srcFull,srcWeightFull, &
         nSend,destFull,destWeightFull,MPI_INFO_NULL,reorder,par%commGraphFull,ierr)
    if(ierr .ne. MPI_SUCCESS) then
       errorcode=ierr
       call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
       print *,subname,': Error after call to MPI_dist_graph_create_adjacent(FULL) ',errorstring
    endif
    allocate(LSchedule%destFull(nSend),LSchedule%srcFull(nRecv))
    LSchedule%destFull(:) = destFull(:)
    LSchedule%srcFull(:)  = srcFull(:)
    ! construct the FULL communication -group- (for one-sided operations):
    call MPI_Comm_group(par%comm, groupFull, ierr)
    call MPI_group_incl(groupFull,nRecv,srcFull,par%groupGraphFull,ierr)
    if (ierr .ne. MPI_SUCCESS) then
       errorcode=ierr
       call MPI_Error_String(errorcode, errorstring, errorlen, ierr)
       print *,subname, ': Error after call to MPI_Comm_group (groupGraphFull) ', errorstring
    endif
    call MPi_Group_size(par%groupGraphFull,sizeGroup,ierr)
    if(Debug) write (*,199) par%rank,sizeGroup,nSend,nRecv

199 format ('RANK: ',i4,' genEdgeSched: size of groupGraphFUll is: ',i8,' nSend, nRecv: ',2(i4))
    deallocate(srcFull,srcWeightFull,destFull,destWeightFull)

    ! construct the INTER communication graph 
    reorder=.FALSE.
    if(numInter>0) then 
       call MPI_Dist_graph_create_adjacent(par%comm, numInter,srcInter,srcWeightInter, &
         numInter,destInter,destWeightInter,MPI_INFO_NULL,reorder,par%commGraphInter,ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_dist_graph_create_adjacent(INTER) ',errorstring
       endif
       deallocate(srcInter,srcWeightInter,destInter,destWeightInter)
    endif

    ! construct the INTRA communication graph 
    reorder=.FALSE.
    if(numIntra>0) then 
       call MPI_Dist_graph_create_adjacent(par%comm, numIntra,srcIntra,srcWeightIntra, &
         numIntra,destIntra,destWeightIntra,MPI_INFO_NULL,reorder,par%commGraphIntra,ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,subname,': Error after call to MPI_dist_graph_create_adjacent(INTRA) ',errorstring
       endif
       deallocate(srcIntra,srcWeightIntra,destIntra,destWeightIntra)
    endif

  200 format ('IAM: ',i4,': ', i2,' of',i2,' comms are interNode')
  201 format ('IAM: ',i4,': ', i2,' of',i2,' comms are intraNode')
#endif


  end subroutine genEdgeSched

  logical function isIntraComm(commranks,rank)

   
    integer, intent(in) :: commranks(:)
    integer, intent(in) :: rank

    integer :: i,nranks

    nranks = SIZE(commranks)
    isIntraComm = .FALSE.
    do i=1,nranks
        if(commranks(i) .eq. rank) then 
           isIntraComm=.TRUE.
        endif
    enddo

  end function isIntraComm

  subroutine CheckSchedule()

    integer                    :: i, nSched, nbufferwords_1, nbufferwords_2
    type (Schedule_t), pointer :: pSchedule

    nSched = SIZE(Schedule)

    do i = 1, nSched
      pSchedule => Schedule(i)
      nbufferwords_1 = SUM(pSchedule%SendCycle(:)%lengthP)
      nbufferwords_2 = SUM(pSchedule%RecvCycle(:)%lengthP)
      if(nbufferwords_1 .ne. nbufferwords_2) then
        write (iulog,100) i,nbufferwords_1, nbufferwords_2
      end if
    end do
100 format('CheckSchedule: ERR IAM:',I3,' SIZEOF(SendBuffer):',I10,' != SIZEOF(RecvBuffer) :',I10)

  end subroutine CheckSchedule

  subroutine PrintSchedule(Schedule)
    ! Debug subroutine for the schedule_t data-structure
    use gridgraph_mod, only : printgridedge

    type (Schedule_t),intent(in),target   :: Schedule(:)
    type (Schedule_t), pointer            :: pSchedule
    type (Cycle_t),pointer                :: pCycle

    integer               :: i,j,nSched

    nSched = SIZE(Schedule)

    write(6,*) '------NEW SCHEDULE FORMAT---------------------'
    do i=1,nSched
       pSchedule => Schedule(i)
       write(6,*)
       write(6,*) '----------------------------------------------'
       write(6,90) i,pSchedule%ncycles
       write(6,*) '----------------------------------------------'
       write(6,*) '-----------SEND-------------------------------'
       do j=1,pSchedule%nSendCycles
          pCycle => pSchedule%SendCycle(j)
          call PrintCycle(pCycle)
          call PrintGridEdge(pCycle%edge%members)
       enddo
       write(6,*) '-----------RECV-------------------------------'
       do j=1,pSchedule%nRecvCycles
          pCycle => pSchedule%RecvCycle(j)
          call PrintCycle(pCycle)
          call PrintGridEdge(pCycle%edge%members)
       enddo
       write(6,*) '-----------MOVE-------------------------------'
       pCycle => pSchedule%MoveCycle(1)
       call PrintCycle(pCycle)
       call PrintGridEdge(pCycle%edge%members)
    enddo
    write(6,*) '-----------Put Index--------------------'
    call PrintIndex(Schedule(1)%pIndx)
    write(6,*) '-----------Get Index--------------------'
    call PrintIndex(Schedule(1)%gIndx)

90  format('NODE # ',I2,2x,'NCYCLES ',I2)
97  format(10x,'EDGE #',I2,2x,'TYPE ',I1,2x,'G.EDGES',I4,2x,'WORDS ',I5,2x, &
         'SRC ',I3,2x,'DEST ',I3,2x,'PTR ',I4)
100 format(15x,I4,5x,I3,1x,'(',I1,') --',I1,'--> ',I3,1x,'(',I1,')')

  end subroutine PrintSchedule

  subroutine PrintIndex(Indx)
  ! Debugging subroutine for the pgindex_t data-structure
    
    !  type, public :: pgindex_t
    !     integer :: elemid
    !     integer :: edgeid
    !     integer :: mesgid
    !     integer :: lenP,lenS
    !  end type pgindex_t

    type (pgindex_t) :: Indx(:)

    integer :: i, len

    len = SIZE(Indx)

    write(6,*) ' elemID,  edgeID,  mesgID, lenP, lenS '
    do i=1,len
       write(6,1099) Indx(i)%elemid,Indx(i)%edgeid,Indx(i)%mesgid,Indx(i)%lenP,Indx(i)%lenS 
    enddo

1099 format(I4,5X,I4,5X,I4,5X,I2,4X,I2)

  end subroutine PrintIndex

  subroutine PrintCycle(Cycle)
  ! debug subroutine for the cycle_t data-structure
    type (Cycle_t),intent(in),target  ::  Cycle

    write(6,97) Cycle%edge%number,Cycle%type,Cycle%edge%nmembers, &
         Cycle%lengthP,Cycle%source, Cycle%dest,Cycle%ptrP

97  format(5x,'METAEDGE #',I2,2x,'TYPE ',I1,2x,'G.EDGES',I4,2x,'WORDS ',I5,2x, &
         'SRC ',I3,2x,'DEST ',I3,2x,'PTR ',I5)

  end subroutine PrintCycle

  subroutine SetCycle(par, elem, schedule,Cycle,Edge,ctype)
    use element_mod,    only: element_t
    use dimensions_mod, only: max_corner_elem, max_neigh_edges
    use cam_abortutils, only: endrun

    type(parallel_t),  intent(in)         :: par
    type(element_t),   intent(inout)      :: elem(:)
    type (Schedule_t), intent(inout)      :: Schedule
    type (Cycle_t),    intent(inout)      :: Cycle
    type (MetaEdge_t), intent(in), target :: Edge
    integer,           intent(in)         :: ctype
    integer                               :: i,il,face, loc, dir

    do i = 1, Edge%nmembers
      if((ctype == HME_CYCLE_SEND) .or. &
         (ctype == HME_CYCLE_MOVE) .or. &
         (ctype == HME_CYCLE_ANY)) then
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
          schedule%pIndx(schedule%pPtr)%mesgid=Edge%HeadVertex-1  ! convert this to 0-based
          schedule%pIndx(schedule%pPtr)%lenP  =Edge%members(i)%wgtP
          schedule%pIndx(schedule%pPtr)%lenS  =Edge%members(i)%wgtS
          if (face.ge.5) then 
             schedule%pIndx(schedule%pPtr)%edgeType = HME_Ordinal
          else
             schedule%pIndx(schedule%pPtr)%edgeType = HME_Cardinal
          endif
          schedule%pPtr=schedule%pPtr+1
        end if
      end if

      if((ctype == HME_CYCLE_RECV) .or. &
         (ctype == HME_CYCLE_MOVE) .or. &
         (ctype == HME_CYCLE_ANY)) then
        !   Setup receive index
        il                     = Global2Local(Edge%members(i)%head%number)
        face                   = Edge%members(i)%head_face
        !need to convert the location of corner elements for getmap and putmap
        if (face.ge.5) then !its a corner
          dir = Edge%members(i)%head_dir
          loc = MOD(dir,max_corner_elem) !this is the location within that direction
          dir = (dir - loc)/max_corner_elem !this is the direction (1-8)
          loc = dir + (dir-5)*(max_corner_elem-1)+loc
          if(loc > max_neigh_edges) then
            write(iulog, *) __FILE__,__LINE__,par%rank,face,i,max_corner_elem,max_neigh_edges,edge%members(i)%head_face
            call endrun('max_neigh_edges set too low.')
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
          schedule%gIndx(schedule%gPtr)%mesgid=Edge%TailVertex-1  ! convert this to 0-based
          schedule%gIndx(schedule%gPtr)%lenP  =Edge%members(i)%wgtP
          schedule%gIndx(schedule%gPtr)%lenS  =Edge%members(i)%wgtS
          if (face.ge.5) then 
             schedule%gIndx(schedule%gPtr)%edgeType = HME_Ordinal
          else
             schedule%gIndx(schedule%gPtr)%edgeType = HME_Cardinal
          endif
          schedule%gPtr=schedule%gPtr+1
        end if
      end if
    end do
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

    integer                :: ptr
    integer, intent(in)    :: inbr,length
    integer, intent(inout) :: tmp(:,:)

    integer                :: i,n

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
