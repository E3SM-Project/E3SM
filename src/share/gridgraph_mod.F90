#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module GridGraph_mod
  !-------------------------
  use kinds, only : real_kind, iulog
  !-------------------------
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  !-----
  implicit none


  private

  integer, public, parameter :: num_neighbors=8 

  type, public :: ptr
      sequence
      integer, dimension(:), pointer :: n => NULL()  ! number of neighbor element
      integer, dimension(:), pointer :: f => NULL()  ! cube face number of neighbor element
  end type ptr

  type, public :: GridVertex_t
      sequence
      type (ptr)                :: nbrs(num_neighbors)
      integer                   :: wgtP(num_neighbors)   ! The weights for edges defined by neighbors array
      integer                   :: wgtP_ghost(num_neighbors)   ! The weights for edges defined by neighbors array
      integer                   :: face_number           ! which face of the cube this vertex is on
      integer                   :: number                ! element number
      integer                   :: processor_number      ! processor number
      integer                   :: SpaceCurve  ! index in Space-Filling curve
  end type GridVertex_t

  type, public :: EdgeIndex_t
      sequence
      integer, pointer            :: ixP(:)
      integer, pointer            :: iyP(:)
  end type EdgeIndex_t

  type, public :: GridEdge_t
      sequence
      integer                      :: head_face  ! needed if head vertex has shape (i.e. square)
      integer                      :: tail_face  ! needed if tail vertex has shape (i.e. square)
      integer                      :: head_ind   ! 
      integer                      :: tail_ind   !
      type (GridVertex_t),pointer  :: head       ! edge head vertex
      type (GridVertex_t),pointer  :: tail       ! edge tail vertex
      logical                      :: reverse

#ifdef TESTGRID
      integer                      :: wgtP        ! amount of information which must be transfered
      type (EdgeIndex_t)           :: HeadIndex
      type (EdgeIndex_t)           :: TailIndex
#endif
  end type GridEdge_t
  
! ==========================================
! Public Interfaces
! ==========================================

  public :: set_GridVertex_number
  public :: PrintGridVertex
 
  public :: initgridedge
  public :: gridedge_search
  public :: gridedge_type
  public :: grid_edge_uses_vertex
  public :: PrintGridEdge
  public :: CheckGridNeighbors
  public :: PrintChecksum

  public :: CreateSubGridGraph
  public :: FreeGraph

  interface assignment ( = )
      module procedure copy_gridedge
      module procedure copy_edgeindex
      module procedure copy_gridvertex
  end interface

contains

! =====================================
! copy edge:
! copy device for overloading = sign.
! =====================================

  recursive subroutine copy_gridedge(edge2,edge1)

    type (GridEdge_t), intent(out) :: edge2
    type (GridEdge_t), intent(in)  :: edge1


    edge2%tail_face = edge1%tail_face
    edge2%head_face = edge1%head_face
    edge2%reverse   = edge1%reverse

    if (associated(edge1%tail)) then
       edge2%tail=>edge1%tail
    end if
    if (associated(edge1%head)) then
       edge2%head=>edge1%head
    end if

#ifdef TESTGRID
    edge2%wgtV       = edge1%wgtV
    edge2%wgtP       = edge1%wgtP

    edge2%TailIndex = edge1%TailIndex
    edge2%HeadIndex = edge1%HeadIndex
#endif

  end subroutine copy_gridedge

  recursive subroutine copy_gridvertex(vertex2,vertex1)
        
    implicit none 

    type (GridVertex_t), intent(out)   :: vertex2
    type (GridVertex_t), intent(in)    :: vertex1

    integer                            :: i,j,n
   
    
!JMD     vertex2%size      = vertex1%size
 
     n = SIZE(vertex1%nbrs)
     do i=1,n
        vertex2%wgtP(i)  = vertex1%wgtP(i)
        vertex2%wgtP_ghost(i)  = vertex1%wgtP_ghost(i)
	if(associated(vertex2%nbrs(i)%n)) &
             deallocate(vertex2%nbrs(i)%n)
	if(associated(vertex2%nbrs(i)%f)) &
             deallocate(vertex2%nbrs(i)%f)
        allocate(vertex2%nbrs(i)%n(SIZE(vertex1%nbrs(i)%n)))
        allocate(vertex2%nbrs(i)%f(SIZE(vertex1%nbrs(i)%n)))
        do j=1,SIZE(vertex1%nbrs(i)%n)
           vertex2%nbrs(i)%n(j) = vertex1%nbrs(i)%n(j)
           vertex2%nbrs(i)%f(j) = vertex1%nbrs(i)%f(j)
        end do
     enddo
     vertex2%number     = vertex1%number
     vertex2%processor_number  = vertex1%processor_number 
     vertex2%SpaceCurve = vertex1%SpaceCurve 

  end subroutine copy_gridvertex

  recursive subroutine copy_edgeindex(index2,index1)
  
  type (EdgeIndex_t), intent(out) :: index2
  type (EdgeIndex_t), intent(in)  :: index1

  if(associated(index1%iyP)) then 
    index2%iyP => index1%iyP
  endif


  if(associated(index1%ixP)) then 
     index2%ixP => index1%ixP
  endif
 
  end subroutine copy_edgeindex

  subroutine FreeGraph(Vertex)

     implicit none
     type (GridVertex_t)           :: Vertex(:)
     integer                       :: i,nelem

     nelem = SIZE(Vertex)

!JMD     do i=1,nelem
!JMD        deallocate(Vertex(i)%wgtV)
!JMD        deallocate(Vertex(i)%wgtG)
!JMD        deallocate(Vertex(i)%nbrs)
!JMD     enddo

  end subroutine FreeGraph

!===========================
! search edge list for match
!===========================

  function gridedge_search(nvert1,nvert2,edge) result(number)

    integer, intent(in) :: nvert1
    integer, intent(in) :: nvert2
    type(GridEdge_t), intent(in) :: edge(:)
    integer :: number

    integer :: tmp
    integer :: head
    integer :: tail

    integer :: nedge
    integer :: i

    nedge=SIZE(edge)

    tail=nvert1
    head=nvert2

    if (tail > head) then
       tmp  = tail
       tail = head
       head = tmp
    end if

    do i=1,nedge
       if (edge(i)%tail%number==tail .and. edge(i)%head%number==head)then
          number=i
       end if
    end do

  end function gridedge_search


  function gridedge_type(edge) result(type)

    use params_mod, only : INTERNAL_EDGE, EXTERNAL_EDGE
    type (GridEdge_t), intent(in)  :: edge
    integer                        :: type

    if (edge%head%processor_number==edge%tail%processor_number) then
        type=INTERNAL_EDGE
    else
        type=EXTERNAL_EDGE
    endif

  end function gridedge_type







  function grid_edge_uses_vertex(Vertex,Edge) result(log)

    type(GridVertex_t), intent(in) :: Vertex
    type(GridEdge_t),   intent(in) :: Edge
    logical :: log
    integer  :: number

    number = Vertex%number
    if(number == Edge%head%number .or. number == Edge%tail%number) then
        log = .TRUE.
    else
        log = .FALSE.
    endif

  end function grid_edge_uses_vertex

  subroutine PrintChecksum(TestPattern,Checksum)

   use dimensions_mod, only : nlev, nelemd, np

   implicit none

   real(kind=real_kind), target,intent(in)   :: TestPattern(:,:,:,:)
   real(kind=real_kind), target,intent(in)   :: Checksum(:,:,:,:)

   integer                                  :: i,k,ix,iy

   print *
   write (iulog,*) 'checksums:'
   do i=1,nelemd
     !  Lets start out only looking at the first element
        write(iulog,*)
        do k=1,nlev
        do iy=1,np
        do ix=1,np
           write(iulog,*)INT(TestPattern(ix,iy,k,i))," checksum = ",INT(Checksum(ix,iy,k,i))
        enddo
        enddo
        enddo
   enddo


  end subroutine PrintChecksum

  subroutine CreateSubGridGraph(Vertex,SVertex,local2global)

    implicit none
        
    type (GridVertex_t),intent(in)         :: Vertex(:)
    type (GridVertex_t),intent(inout)      :: SVertex(:)
    integer,intent(in)                     :: local2global(:)

    integer                                :: nelem,nelem_s,n,ncount
    integer                                :: inbr,i,ig,j,k
    
    integer,allocatable                    :: global2local(:)
    logical, parameter    :: Debug = .FALSE.


    nelem   = SIZE(Vertex)
    nelem_s = SiZE(SVertex) 

    if(Debug) write(iulog,*)'CreateSubGridGraph: point #1'
    allocate(global2local(nelem))
    if(Debug) write(iulog,*)'CreateSubGridGraph: point #2'

    global2local(:) = 0
    do i=1,nelem_s
        ig = local2global(i)
        global2local(ig) = i
    enddo
    if(Debug) write(iulog,*)'CreateSubGridGraph: point #3'
    do i=1,nelem_s
        ig = local2global(i)

        if(Debug) write(iulog,*)'CreateSubGridGraph: point #4'
        call copy_gridvertex(SVertex(i),Vertex(ig))
        n = SIZE(SVertex(i)%nbrs(:))
        ! ==============================================
        ! Apply the correction to the neighbors list to 
        ! reflect new subgraph numbers 
        ! ==============================================
        ncount=0
        if(Debug) write(iulog,*)'CreateSubGridGraph: point #5'
        do j=1,n
           if(ASSOCIATED(Svertex(i)%nbrs(j)%n)) then 
              do k = 1,SIZE(Svertex(i)%nbrs(j)%n)
                 if(Debug) write(iulog,*)'CreateSubGridGraph: point #5.1 size(global2local) Svertex(i)%nbrs(j) ', &
                              size(global2local), Svertex(i)%nbrs(j)%n(k)
                 inbr = global2local(Svertex(i)%nbrs(j)%n(k))
                 if(Debug) write(iulog,*)'CreateSubGridGraph: point #5.2'
                 if(inbr .gt. 0) then 
                    Svertex(i)%nbrs(j)%n(k) = inbr
                    ncount = ncount+1
                 else 
                    Svertex(i)%wgtP(j)  = 0        
                    Svertex(i)%wgtP_ghost(j)  = 0        
                    DEALLOCATE(Svertex(i)%nbrs(j)%n)
                    DEALLOCATE(Svertex(i)%nbrs(j)%f)
                 endif
              enddo
           endif
           if(Debug) write(iulog,*)'CreateSubGridGraph: point #5.3'
        enddo
    if(Debug) write(iulog,*)'CreateSubGridGraph: point #6'
        Svertex(i)%number = i
     enddo
    if(Debug) write(iulog,*)'CreateSubGridGraph: point #7'
     deallocate(global2local)
    if(Debug) write(iulog,*)'CreateSubGridGraph: point #8'

  end subroutine CreateSubGridGraph

  subroutine PrintGridEdge(Edge)

    implicit none
    type (GridEdge_t), intent(in) :: Edge(:)

    integer           :: i,nedge,ii,wgtP

    nedge = SIZE(Edge)

    write(iulog,95)
    do i=1,nedge
          ii=Edge(i)%tail_face
          wgtP=Edge(i)%tail%wgtP(ii)       
          write(iulog,100) i, &
               Edge(i)%tail%number,Edge(i)%tail_face, wgtP, &
               Edge(i)%head%number,Edge(i)%head_face, gridedge_type(Edge(i))
    enddo
  95 format(5x,'GRIDEDGE #',3x,'Tail (face)',5x,'Head (face)',3x,'Type')
 100 format(10x,I6,8x,I4,1x,'(',I1,')  --',I2,'--> ',I6,1x,'(',I1,')',5x,'[',I1,']')

  end subroutine PrintGridEdge

! ==========================================
! set_GridVertex_neighbors:
!
! Set global element number for element elem
! ==========================================

  subroutine set_GridVertex_number(elem,number)

    type(GridVertex_t)         :: elem
    integer                 :: number

    elem%number=number

  end subroutine set_GridVertex_number

  subroutine PrintGridVertex(Vertex)

    implicit none 
    type (GridVertex_t), intent(in),target :: Vertex(:)
  
    integer        :: i,nvert

    nvert = SIZE(Vertex)
        
    write(iulog,98)
    do i=1,nvert
       write(iulog,99) Vertex(i)%number, Vertex(i)%processor_number, &
            Vertex(i)%nbrs(west)%n(:),Vertex(i)%wgtP(west), &
            Vertex(i)%nbrs(east)%n(:),Vertex(i)%wgtP(east), &
            Vertex(i)%nbrs(south)%n(:),Vertex(i)%wgtP(south), &
            Vertex(i)%nbrs(north)%n(:),Vertex(i)%wgtP(north), &
            Vertex(i)%nbrs(swest)%n(:),Vertex(i)%wgtP(swest), &
            Vertex(i)%nbrs(seast)%n(:),Vertex(i)%wgtP(seast), &
            Vertex(i)%nbrs(nwest)%n(:),Vertex(i)%wgtP(nwest), &
            Vertex(i)%nbrs(neast)%n(:),Vertex(i)%wgtP(neast)
    enddo
  98 format(5x,'GRIDVERTEX #',2x,'PART',2x,'DEG',4x,'W',8x,'E',8x, &
                'S',8x,'N',7x,'SW',7x,'SE',7x,'NW',7x,'NE')
  99 format(10x,I3,8x,I1,2x,I1,2x,8(2x,I3,1x,'(',I2,')'))


  end subroutine PrintGridVertex
  subroutine CheckGridNeighbors(Vertex)
  
  implicit none
  type (GridVertex_t), intent(in) :: Vertex(:)

  integer :: i,j,k,l,m,nnbrs,inbrs,nvert
  nvert = SIZE(Vertex)

  do i=1,nvert
        nnbrs = SIZE(Vertex(i)%nbrs)
        do j=1,nnbrs
           do l=1,SIZE(Vertex(i)%nbrs(j)%n)
              inbrs = Vertex(i)%nbrs(j)%n(l)
              if(inbrs > 0) then
                 do k=1,nnbrs
                    do m=1,SIZE(Vertex(i)%nbrs(k)%n)
                        if( inbrs .eq. Vertex(i)%nbrs(k)%n(m) .and. (j/=k .or. l/=m)) &
                           write(iulog,*)'CheckGridNeighbors: ERROR identical neighbors detected  for Vertex ',i
             
                    enddo
                 enddo
              endif
           enddo
        enddo
  enddo

  end subroutine CheckGridNeighbors


  subroutine initgridedge(GridEdge,GridVertex)
  use parallel_mod, only : abortmp

  implicit none

  type (GridEdge_t), intent(inout)       :: GridEdge(:)
  type (GridVertex_t), intent(in),target :: GridVertex(:)

  integer                                :: i,j,k,m,l,iptr,wgtV,wgtP
  integer                                :: nelem,nelem_edge,inbr
  logical                                :: Verbose=.FALSE.

  nelem      = SIZE(GridVertex)
  nelem_edge = SIZE(GridEdge)

  GridEdge(:)%reverse=.FALSE.

  iptr=1
  do j=1,nelem
     do i=1,num_neighbors           
        if((GridVertex(j)%wgtP(i) .gt. 0).and.(associated(GridVertex(j)%nbrs(i)%n))) then    ! Do this only if has a non-zero weight
           do l=1,SIZE(GridVertex(j)%nbrs(i)%n)
             if (nelem_edge<iptr) call abortmp('Error in initgridedge: Number of edges greater than expected.')
             GridEdge(iptr)%tail      => GridVertex(j)
             GridEdge(iptr)%tail_face =  i
             GridEdge(iptr)%tail_ind  =  l
             inbr                     =  GridVertex(j)%nbrs(i)%n(l)
             GridEdge(iptr)%head      => GridVertex(inbr)
#ifdef TESTGRID
             wgtV                     =  GridVertex(j)%wgtV(i)
             GridEdge(iptr)%wgtV      =  wgtV
             wgtP                     =  GridVertex(j)%wgtP(i)
             GridEdge(iptr)%wgtP      =  wgtP
             ! allocate the indirect addressing arrays
             allocate(GridEdge(iptr)%TailIndex%ixV(wgtV))
             allocate(GridEdge(iptr)%TailIndex%ixP(wgtP))

             allocate(GridEdge(iptr)%TailIndex%iyV(wgtV))
             allocate(GridEdge(iptr)%TailIndex%iyP(wgtP))

             allocate(GridEdge(iptr)%HeadIndex%ixV(wgtV))
             allocate(GridEdge(iptr)%HeadIndex%ixP(wgtP))

             allocate(GridEdge(iptr)%HeadIndex%iyV(wgtV))
             allocate(GridEdge(iptr)%HeadIndex%iyP(wgtP))
#endif
              ! ===========================================
              ! Need this aweful piece of code to determine
              ! which "face" of the neighbor element the
              ! edge links (i.e. the "head_face")
              ! ===========================================
 
              do k=1,num_neighbors
                 if (associated(GridVertex(inbr)%nbrs(k)%n)) then
                    do m=1,SIZE(GridVertex(inbr)%nbrs(k)%n)
                       if(GridVertex(inbr)%nbrs(k)%n(m) == GridVertex(j)%number) then
                          GridEdge(iptr)%head_face=k
                          GridEdge(iptr)%head_ind =m
                       endif
                    enddo
                 end if
              enddo
              iptr=iptr+1
           enddo
        endif
     end do
  end do
  if (nelem_edge+1 /= iptr) call abortmp('Error in initgridedge: Number of edges less than expected.')
  if (Verbose) then

     print *
     write(iulog,*)"element edge tail,head list: (TEST)"
     do i=1,nelem_edge
        write(iulog,*)GridEdge(i)%tail%number,GridEdge(i)%head%number
     end do

     print *
     write(iulog,*)"element edge tail_face, head_face list: (TEST)"
     do i=1,nelem_edge
        write(iulog,*)GridEdge(i)%tail_face,GridEdge(i)%head_face
     end do
  end if

  end subroutine initgridedge

end module GridGraph_mod
