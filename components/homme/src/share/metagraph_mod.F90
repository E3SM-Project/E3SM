#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!************************metagraph_mod.F****************************************

module metagraph_mod
!
!  Revisions:
!  Mark Taylor 2018/3  fix memory leak
!
  use kinds, only : int_kind, iulog
  use gridgraph_mod, only : gridvertex_t, gridedge_t, &
       deallocate_gridvertex_nbrs, allocate_gridvertex_nbrs, assignment ( = )

  use pio_types ! _EXTERNAL

  implicit none 
  private 

  type, public :: MetaEdge_t
     type (GridEdge_t),pointer     :: members(:)
     integer          ,pointer     :: edgeptrP(:)
     integer          ,pointer     :: edgeptrP_ghost(:)
     integer          ,pointer     :: edgeptrS(:)
     integer                       :: number
     integer                       :: type 
     integer                       :: wgtP       ! sum of lengths of all messages to pack for edges
     integer                       :: wgtP_ghost ! sum of lengths of all messages to pack for ghost cells
     integer                       :: wgtS
     integer                       :: HeadVertex ! processor number to send to
     integer                       :: TailVertex ! processor number to send from
     integer                       :: nmembers   ! number of messages to (un)pack (out)into this buffer
     integer                       :: padding    ! just to quite compiler
  end type MetaEdge_t

  type, public :: MetaVertex_t             ! one for each processor
     integer                       :: number     ! USELESS just the local processor number
     integer                       :: nmembers   ! number of elements on this processor 
     type (GridVertex_t),pointer   :: members(:) ! array of elements on this processor
     type (MetaEdge_t),pointer     :: edges(:)   ! description of messages to send/receive
     integer                       :: nedges     ! number of processors to communicate with (length of edges)
     integer                       :: padding    ! just to quite compiler
  end type MetaVertex_t


  ! public :: findedge
  public :: edge_uses_vertex
  public :: PrintMetaEdge, PrintMetaVertex
  public :: LocalElemCount
  !public :: MetaEdgeCount
  public :: initMetaGraph, destroyMetaGraph

  interface assignment ( = )
     module procedure copy_metaedge
  end interface


contains 
  ! =====================================
  ! copy vertex:
  ! copy device for overloading = sign.
  ! =====================================

  recursive subroutine copy_metaedge(edge2,edge1)

    type (MetaEdge_t), intent(out) :: edge2
    type (MetaEdge_t), intent(in)  :: edge1

    integer i

    edge2%number   = edge1%number
    edge2%type     = edge1%type
    edge2%wgtP      = edge1%wgtP
    edge2%wgtP_ghost      = edge1%wgtP_ghost
    edge2%nmembers = edge1%nmembers

    if (associated(edge1%members)) then
       allocate(edge2%members(edge2%nmembers))
       do i=1,edge2%nmembers
          edge2%members(i)=edge1%members(i)
       end do
    end if

    if (associated(edge1%edgeptrP)) then
       allocate(edge2%edgeptrP(edge2%nmembers))
       allocate(edge2%edgeptrS(edge2%nmembers))
       allocate(edge2%edgeptrP_ghost(edge2%nmembers))
       do i=1,edge2%nmembers
          edge2%edgeptrP(i)=edge1%edgeptrP(i)
          edge2%edgeptrS(i)=edge1%edgeptrS(i)
          edge2%edgeptrP_ghost(i)=edge1%edgeptrP_ghost(i)
       end do
    end if

    edge2%HeadVertex = edge1%HeadVertex
    edge2%TailVertex = edge1%TailVertex

  end subroutine copy_metaedge


! function findedge(mEdge,Edge) result(number)

!   type(MetaEdge_t), intent(inout) :: mEdge(:)
!   type(GridEdge_t), intent(in)    :: Edge
!   integer :: number

!   integer :: head,tail,exist

!   integer :: nedge
!   integer :: i

!   nedge=SIZE(mEdge)
!   number = 0
!   tail=Edge%tail%processor_number
!   head=Edge%head%processor_number

!   exist=0
!   do i=1,nedge
!      !       write(iulog,*)'mEdge(i)%number: ',mEdge(i)%number
!      if(mEdge(i)%number .ne. 0) then 
!         if  ((mEdge(i)%TailVertex==tail .and. mEdge(i)%HeadVertex==head) ) then
!            number=i
!         end if
!         exist=exist+1
!      endif
!   end do
!   if(number == 0) number = exist + 1

! end function findedge

! function MetaEdgeCount(Edge) result(nMedge)
!   implicit none

!   type (GridEdge_t),intent(in)  :: Edge(:)
!   integer                       :: nMedge

!   integer                       :: nedge,i,j,maxedges
!   integer                       :: head_processor_number,tail_processor_number
!   integer, allocatable          :: tmp(:,:)
!   logical                       :: found


!   nedge = SIZE(Edge)
!   maxedges=nedge

!   allocate(tmp(2,maxedges))
!   tmp=0
!   nMedge=0
!   do i=1,nedge
!      head_processor_number = Edge(i)%head%processor_number
!      tail_processor_number = Edge(i)%tail%processor_number
!      found = .FALSE.
!      do j=1,nMedge
!         if((tmp(1,j) .eq. head_processor_number).and. &
!              (tmp(2,j) .eq. tail_processor_number)) found=.TRUE.
!      enddo
!      if(.NOT. found) then 
!         nMedge=nMedge+1
!         tmp(1,nMedge) = head_processor_number
!         tmp(2,nMedge) = tail_processor_number
!      endif
!   enddo
!   !mem    write(iulog,*)'MetaEdgeCount: before call to deallocate(tmp)'
!   deallocate(tmp)

! end function MetaEdgeCount

  function LocalElemCount(Vertex) result(nelemd)
    implicit none

    type (MetaVertex_t),intent(in)  :: Vertex
    integer                         :: nelemd

    nelemd=Vertex%nmembers


  end function LocalElemCount


  function edge_uses_vertex(Vertex,Edge) result(log)

    implicit none
    type(MetaVertex_t), intent(in) :: Vertex
    type(MetaEdge_t),   intent(in) :: Edge
    logical :: log
    integer  :: number

    number = Vertex%number
    if(number == Edge%HeadVertex .or. number == Edge%TailVertex) then 
       log = .TRUE.
    else
       log = .FALSE. 
    endif

  end function edge_uses_vertex

  subroutine PrintMetaEdge(Edge)
    use gridgraph_mod, only : PrintGridEdge
    implicit none
    type (MetaEdge_t), intent(in) :: Edge(:)
    integer          :: i,nedge

    nedge = SIZE(Edge)
    do i=1,nedge
       print *
       write(iulog,90)  Edge(i)%number,Edge(i)%type,Edge(i)%wgtP,Edge(i)%nmembers, &
            Edge(i)%TailVertex, Edge(i)%HeadVertex
       if(associated(Edge(i)%members)) then
          call PrintGridEdge(Edge(i)%members)
       endif
    enddo
90  format('METAEDGE #',I4,2x,'TYPE ',I1,2x,'WGT ',I4,2x,'NUM ',I6,2x,'Processors ',I4,' ---> ',I4)

  end subroutine PrintMetaEdge

  subroutine PrintMetaVertex(Vertex)
    use gridgraph_mod, only : PrintGridVertex
    implicit none 
    type (MetaVertex_t), intent(in),target  :: Vertex

    integer               :: j


    write(iulog,*)
    write(iulog,95) Vertex%nmembers
    call PrintGridVertex(Vertex%members)
    write(iulog,96) Vertex%nedges
    if(associated(Vertex%edges)) then 
       do j=1,Vertex%nedges
          write(iulog,97) Vertex%edges(j)%number,     Vertex%edges(j)%type,              &
               Vertex%edges(j)%wgtP,        Vertex%edges(j)%HeadVertex,        &
               Vertex%edges(j)%TailVertex
       enddo
    endif

95  format(5x,I2,' Member Grid Vertices')
96  format(5x,I2,' Incident Meta Edges ')
97  format(10x,'METAEDGE #',I2,2x,'TYPE ',I1,2x,'WGT ',I4,2x,'Processors ',I2,' ---> ',I2)

  end subroutine PrintMetaVertex

  subroutine  initMetaGraph(ThisProcessorNumber,MetaVertex,GridVertex,GridEdge)
    use ll_mod, only : root_t, LLSetEdgeCount, LLFree, LLInsertEdge, LLGetEdgeCount, LLFindEdge
    use gridgraph_mod, only : GridEdge_type, printGridVertex
    !------------------
    !------------------
    implicit none

    integer(kind=int_kind), intent(in)      :: ThisProcessorNumber
    type (MetaVertex_t), intent(out)        :: MetaVertex
    type (GridVertex_t), intent(in),target  :: GridVertex(:)
    type (GridEdge_t),   intent(in),target  :: GridEdge(:)

    !type (MetaEdge_t), allocatable :: MetaEdge(:)
    integer                          :: nelem,nelem_edge, nedges  
    integer,allocatable              :: icount(:)
    integer                          :: ic,i,j,ii
    integer                          :: npart
    integer                          :: head_processor_number
    integer                          :: tail_processor_number
    integer :: nedge_active,enum
    logical :: found
    integer iTail, iHead, wgtP,wgtS

    type (root_t) :: mEdgeList ! root_t = C++ std::set<std::pair<int,int> >

    logical  :: Verbose = .FALSE.
    logical  :: Debug = .FALSE.


    if(Debug) write(iulog,*)'initMetagraph: point #1'
    !  Number of grid vertices
    nelem      = SIZE(GridVertex)
    !  Number of grid edges
    nelem_edge = SIZE(GridEdge)

    mEdgeList%number = ThisProcessorNumber
    NULLIFY(mEdgeList%first)
    call LLSetEdgeCount(0)

#if 0
    ! look for internal (move) edges and add them first 
    do i=1,nelem_edge
       tail_processor_number = GridEdge(i)%tail%processor_number
       head_processor_number = GridEdge(i)%head%processor_number
       if(tail_processor_number  .eq. ThisProcessorNumber .and.  &
          head_processor_number  .eq. ThisProcessorNumber ) then
          call LLInsertEdge(mEdgeList,tail_processor_number,head_processor_number,eNum)
       endif
    enddo
    ! next add the off processor (send/receive) edges 
    do i=1,nelem_edge
       tail_processor_number = GridEdge(i)%tail%processor_number
       head_processor_number = GridEdge(i)%head%processor_number
       if((tail_processor_number  .eq. ThisProcessorNumber .or.  &
          head_processor_number  .eq. ThisProcessorNumber) .and. &
          (tail_processor_number .ne. head_processor_number) ) then
          call LLInsertEdge(mEdgeList,tail_processor_number,head_processor_number,eNum)
       endif
    enddo

#else
    do i=1,nelem_edge
       tail_processor_number = GridEdge(i)%tail%processor_number
       head_processor_number = GridEdge(i)%head%processor_number
       if(tail_processor_number  .eq. ThisProcessorNumber .or.  &
          head_processor_number  .eq. ThisProcessorNumber ) then 
          call LLInsertEdge(mEdgeList,tail_processor_number,head_processor_number,eNum)
       endif
    enddo
#endif
    call LLGetEdgeCount(nedges)

    NULLIFY(MetaVertex%edges)
        
    allocate(MetaVertex%edges(nedges))

    ! Initalize the Meta Vertices to zero... probably should be done
    ! in a separate routine
    MetaVertex%nmembers=0  
    MetaVertex%number=0
    MetaVertex%nedges=0
    if(Debug) write(iulog,*)'initMetagraph: point #2'


    !  Give some identity to the Meta_vertex
    MetaVertex%number   = ThisProcessorNumber
    if(Debug) write(iulog,*)'initMetagraph: point #3'

    !  Look through all the small_vertices and determine the number of
    !  member vertices
    if(Debug) call PrintGridVertex(GridVertex)
    if(Debug) write(iulog,*)'initMetagraph: After call to PrintGridVertex point #3.1'
    if(Debug) write(iulog,*)'initMetaGraph: ThisProcessorNumber is ',ThisProcessorNumber

    do j=1,nelem  ! count number of elements on this processor
       if(GridVertex(j)%processor_number .eq. ThisProcessorNumber) then
         MetaVertex%nmembers = MetaVertex%nmembers + 1
       endif
    enddo

    if(Debug) write(iulog,*)'initMetagraph: point #4 '
    !  Allocate space for the members of the MetaVertices
    if(Debug) write(iulog,*)'initMetagraph: point #4.1 i,MetaVertex%nmembers',i,MetaVertex%nmembers
    allocate(MetaVertex%members(MetaVertex%nmembers))

    ! dont allocate %members - struct will be allocated in the overloaded assignement
    ! operation below
!    do j=1, MetaVertex%nmembers
!       call allocate_gridvertex_nbrs(MetaVertex%members(j))
!    end do
    
    if(Debug) write(iulog,*)'initMetagraph: point #5'

    !  Set the identity of the members of the MetaVertices
    ic=1
    do j=1,nelem
       if( GridVertex(j)%processor_number .eq. ThisProcessorNumber) then
          MetaVertex%members(ic) = GridVertex(j)
          ic=ic+1
       endif
    enddo

    nedges = SIZE(MetaVertex%edges)
    if(Debug) write(iulog,*)'initMetagraph: point #6 nedges',nedges
    !  Zero out all the edge numbers ... this should probably be
    !  move to some initalization routine
    MetaVertex%edges%number   = 0
    MetaVertex%edges%nmembers = 0
    MetaVertex%edges%wgtP     = 0
    MetaVertex%edges%wgtS     = 0
    MetaVertex%edges%wgtP_ghost     = 0
    do i=1,nedges
       NULLIFY(MetaVertex%edges(i)%members)
    enddo

    if(Debug) write(iulog,*)'initMetagraph: point #7'

    ! Insert all the grid edges into the Meta Edges
    do i=1, nelem_edge
       !  Which Meta Edge does this grid edge belong
       head_processor_number = GridEdge(i)%head%processor_number
       tail_processor_number = GridEdge(i)%tail%processor_number
       call LLFindEdge(mEdgeList,tail_processor_number,head_processor_number,j,found)
       if(found) then 

          !  Increment the number of grid edges contained in the grid edge
          !  and setup the pointers
          if(Debug) write(iulog,*)'initMetagraph: point #8'
          ii=GridEdge(i)%tail_face

          wgtP=Gridedge(i)%tail%nbrs_wgt(ii)
          wgtS=1

          MetaVertex%edges(j)%nmembers                = MetaVertex%edges(j)%nmembers+1
          MetaVertex%edges(j)%wgtP                    = MetaVertex%edges(j)%wgtP + wgtP
          MetaVertex%edges(j)%wgtS                    = MetaVertex%edges(j)%wgtS + wgtS

          MetaVertex%edges(j)%wgtP_ghost              = MetaVertex%edges(j)%wgtP_ghost + Gridedge(i)%tail%nbrs_wgt_ghost(ii)

          if(Debug) write(iulog,*)'initMetagraph: point #9'

          !  If this the first grid edge to be inserted into the Meta Edge
          !  do some more stuff

          if(MetaVertex%edges(j)%nmembers .eq. 1) then

             if(Debug) write(iulog,*)'initMetagraph: point #10'
             MetaVertex%edges(j)%number   = j                      ! its identity
             MetaVertex%edges(j)%type  = gridedge_type(GridEdge(i))  ! Type of grid edge

             if(Debug) write(iulog,*)'initMetagraph: point #11'

             !  Setup the pointer to the head and tail of the Vertex
             MetaVertex%edges(j)%HeadVertex          = head_processor_number
             MetaVertex%edges(j)%TailVertex          = tail_processor_number
             if(Debug) write(iulog,*)'initMetagraph: point #12'

             !  Determine the number of edges for the Meta_Vertex
             !  This is the number of processors to communicate with
             MetaVertex%nedges =  MetaVertex%nedges + 1
             if(Debug) write(iulog,*)'initMetagraph: point #13'
          endif
       endif
    enddo

    do i=1,nedges
       !  Allocate space for the member edges and edge index
       allocate(MetaVertex%edges(i)%members (MetaVertex%edges(i)%nmembers))
       allocate(MetaVertex%edges(i)%edgeptrP(MetaVertex%edges(i)%nmembers))
       allocate(MetaVertex%edges(i)%edgeptrS(MetaVertex%edges(i)%nmembers))
       allocate(MetaVertex%edges(i)%edgeptrP_ghost(MetaVertex%edges(i)%nmembers))
       MetaVertex%edges(i)%edgeptrP(:)=0
       MetaVertex%edges(i)%edgeptrS(:)=0
       MetaVertex%edges(i)%edgeptrP_ghost(:)=0
    enddo
    if(Debug) write(iulog,*)'initMetagraph: point #14'

    !  Insert the edges into the proper meta edges
    allocate(icount(nelem_edge))
    icount=1
    do i=1,nelem_edge
       head_processor_number = GridEdge(i)%head%processor_number
       tail_processor_number = GridEdge(i)%tail%processor_number
       call LLFindEdge(mEdgeList,tail_processor_number,head_processor_number,j,found)
       if(found) then 
          MetaVertex%edges(j)%members(icount(j)) = GridEdge(i)
          if(icount(j)+1 .le. MetaVertex%edges(j)%nmembers) then

             ii=GridEdge(i)%tail_face

             wgtP=Gridedge(i)%tail%nbrs_wgt(ii)
             MetaVertex%edges(j)%edgeptrP(icount(j)+1) = MetaVertex%edges(j)%edgeptrP(icount(j)) + wgtP
            
             wgtS = 1
             MetaVertex%edges(j)%edgeptrS(icount(j)+1) = MetaVertex%edges(j)%edgeptrS(icount(j)) + wgtS

             wgtP=Gridedge(i)%tail%nbrs_wgt_ghost(ii)
             MetaVertex%edges(j)%edgeptrP_ghost(icount(j)+1) = MetaVertex%edges(j)%edgeptrP_ghost(icount(j)) + wgtP
          endif
          if(Debug) write(iulog,*)'initMetagraph: point #15'
          icount(j)=icount(j)+1
       endif
    enddo
    deallocate(icount)
    if(Debug) write(iulog,*)'initMetagraph: point #16'

    if(Verbose) then
       print *
       write(iulog,*)"edge bundle list:(INITMETAGRAPH)"
       call PrintMetaEdge( MetaVertex%edges)
       write(iulog,*)'initmetagrap: Before last call to PrintMetaVertex'
       call PrintMetaVertex(MetaVertex)
    endif

    call LLFree(mEdgeList)

90  format('EDGE #',I2,2x,'TYPE ',I1,2x,'Processor Numbers ',I2,' ---> ',I2)
100 format(10x,I2,1x,'(',I1,') ---> ',I2,1x,'(',I1,')')

  end subroutine initMetaGraph

  subroutine destroyMetaGraph(MetaVertex)
    use gridgraph_mod, only: deallocate_gridvertex_nbrs

    type (MetaVertex_t), intent(inout) :: MetaVertex
    integer :: j

    do j = 1, MetaVertex%nmembers
       call deallocate_gridvertex_nbrs(MetaVertex%members(j))
    end do
    do j = 1, MetaVertex%nedges
       deallocate(MetaVertex%edges(j)%members)
       deallocate(MetaVertex%edges(j)%edgeptrP)
       deallocate(MetaVertex%edges(j)%edgeptrS)
       deallocate(MetaVertex%edges(j)%edgeptrP_ghost)
    end do
    deallocate(MetaVertex%edges)
    deallocate(MetaVertex%members)
  end subroutine destroyMetaGraph

end module metagraph_mod
