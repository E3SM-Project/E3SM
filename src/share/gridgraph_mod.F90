#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module GridGraph_mod
  !-------------------------
  use kinds, only : real_kind, iulog
  !-------------------------------
  use dimensions_mod, only : max_neigh_edges
  !-------------------------
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  !-----
  implicit none


  private

  integer, public, parameter :: num_neighbors=8 ! for north, south, east, west, neast, nwest, seast, swest


  type, public :: GridVertex_t

      integer, pointer          :: nbrs(:)                     ! The numbers of the neighbor elements
      integer, pointer          :: nbrs_face(:)                ! The cube face number of the neighbor element (nbrs array)
      integer, pointer          :: nbrs_wgt(:)                 ! The weights for edges defined by nbrs array
      integer, pointer          :: nbrs_wgt_ghost(:)           ! The weights for edges defined by nbrs array
      integer                   :: nbrs_ptr(num_neighbors + 1) !index into the nbrs array for each neighbor direction

      integer                   :: face_number           ! which face of the cube this vertex is on
      integer                   :: number                ! element number
      integer                   :: processor_number      ! processor number
      integer                   :: SpaceCurve  ! index in Space-Filling curve
  end type GridVertex_t

  type, public :: EdgeIndex_t
      integer, pointer            :: ixP(:)
      integer, pointer            :: iyP(:)
  end type EdgeIndex_t

  type, public :: GridEdge_t
      integer                      :: head_face  ! needed if head vertex has shape (i.e. square)
      integer                      :: tail_face  ! needed if tail vertex has shape (i.e. square)
      integer                      :: head_dir   !which of 8 neighbor directions is the head
      integer                      :: tail_dir   !which of 8 neighbor directions is the tail
      type (GridVertex_t),pointer  :: head       ! edge head vertex
      type (GridVertex_t),pointer  :: tail       ! edge tail vertex
      logical                      :: reverse

  end type GridEdge_t
  
! ==========================================
! Public Interfaces
! ==========================================

  public :: set_GridVertex_number
  public :: PrintGridVertex
 
  public :: allocate_gridvertex_nbrs
  public :: deallocate_gridvertex_nbrs
  public :: initgridedge
  public :: gridedge_search
  public :: gridedge_type
  public :: grid_edge_uses_vertex
  public :: PrintGridEdge
  public :: CheckGridNeighbors
  public :: PrintChecksum

  public :: CreateSubGridGraph
  public :: FreeGraph

  public :: assignment ( = ) 

  interface assignment ( = )
      module procedure copy_gridedge
      module procedure copy_edgeindex
      module procedure copy_gridvertex
  end interface

contains

!======================================================================

  subroutine allocate_gridvertex_nbrs(vertex, dim)

    type (GridVertex_t), intent(inout)   :: vertex
    integer, optional, intent(in)        :: dim
    integer                              :: num

    if (present(dim)) then
       num = dim
    else
       num = max_neigh_edges
    end if

    allocate(vertex%nbrs(num))
    allocate(vertex%nbrs_face(num))
    allocate(vertex%nbrs_wgt(num))
    allocate(vertex%nbrs_wgt_ghost(num))
 

  end subroutine allocate_gridvertex_nbrs
!======================================================================

  subroutine deallocate_gridvertex_nbrs(vertex)

    type (GridVertex_t), intent(inout)   :: vertex

    deallocate(vertex%nbrs)
    deallocate(vertex%nbrs_face)
    deallocate(vertex%nbrs_wgt)
    deallocate(vertex%nbrs_wgt_ghost)
 
  end subroutine deallocate_gridvertex_nbrs

!======================================================================

! =====================================
! copy edge:
! copy device for overloading = sign.
! =====================================


  recursive subroutine copy_gridedge(edge2, edge1)

    type (GridEdge_t), intent(out) :: edge2
    type (GridEdge_t), intent(in)  :: edge1

    edge2%tail_face = edge1%tail_face
    edge2%head_face = edge1%head_face
    edge2%tail_dir = edge1%tail_dir
    edge2%head_dir = edge1%head_dir
    edge2%reverse   = edge1%reverse

    if (associated(edge1%tail)) then
       edge2%tail=>edge1%tail
    end if
    if (associated(edge1%head)) then
       edge2%head=>edge1%head
    end if

  end subroutine copy_gridedge

!======================================================================

  recursive subroutine copy_gridvertex(vertex2, vertex1)
        
    implicit none 

    type (GridVertex_t), intent(out)   :: vertex2
    type (GridVertex_t), intent(in)    :: vertex1

    integer                            :: i,j,n
   
     n = SIZE(vertex1%nbrs)

     if (associated(vertex2%nbrs)) then
        nullify(vertex2%nbrs)
     end if
     if (associated(vertex2%nbrs_face)) then
        nullify(vertex2%nbrs_face)
     end if
     if (associated(vertex2%nbrs_wgt)) then
        nullify(vertex2%nbrs_wgt)
     end if
     if (associated(vertex2%nbrs_wgt_ghost)) then
        nullify(vertex2%nbrs_wgt_ghost)
     end if

     call allocate_gridvertex_nbrs(vertex2)

     do i=1,n
        vertex2%nbrs(i) = vertex1%nbrs(i)
        vertex2%nbrs_face(i) = vertex1%nbrs_face(i)
        vertex2%nbrs_wgt(i)  = vertex1%nbrs_wgt(i)
        vertex2%nbrs_wgt_ghost(i)  = vertex1%nbrs_wgt_ghost(i)
     enddo

     do i=1, num_neighbors+1
        vertex2%nbrs_ptr(i) = vertex1%nbrs_ptr(i)
     enddo

     vertex2%face_number     = vertex1%face_number
     vertex2%number     = vertex1%number
     vertex2%processor_number  = vertex1%processor_number 
     vertex2%SpaceCurve = vertex1%SpaceCurve 

  end subroutine copy_gridvertex

!======================================================================
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

!======================================================================
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

!======================================================================

!===========================
! search edge list for match
!===========================

  function gridedge_search(nvert1, nvert2, edge) result(number)

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

!======================================================================

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

!======================================================================



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

!======================================================================

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

!======================================================================

  subroutine CreateSubGridGraph(Vertex, SVertex, local2global)

    implicit none
        
    type (GridVertex_t),intent(in)         :: Vertex(:)
    type (GridVertex_t),intent(inout)      :: SVertex(:)
    integer,intent(in)                     :: local2global(:)

    integer                                :: nelem,nelem_s,n,ncount,cnt,pos, orig_start
    integer                                :: inbr,i,ig,j,k, new_pos
    
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
        call copy_gridvertex(SVertex(i),Vertex(ig))  !svertex(i) = vertex(ig)

        n = SIZE(SVertex(i)%nbrs(:))
        ! ==============================================
        ! Apply the correction to the neighbors list to 
        ! reflect new subgraph numbers 
        ! ==============================================


        if(Debug) write(iulog,*)'CreateSubGridGraph: point #5'
        
        orig_start = 1

        do j=1,num_neighbors

           cnt = Svertex(i)%nbrs_ptr(j+1) - orig_start  !number of neighbors for this direction
           ncount = 0
           do k = 1, cnt
              pos = orig_start + k-1
              if(Debug) write(iulog,*)'CreateSubGridGraph: point #5.1 size(global2local) Svertex(i)%nbrs(j) ', &
                   size(global2local), Svertex(i)%nbrs(pos)

                 inbr = global2local(Svertex(i)%nbrs(pos))

                 if(Debug) write(iulog,*)'CreateSubGridGraph: point #5.2'
                 if(inbr .gt. 0) then 
                    new_pos = Svertex(i)%nbrs_ptr(j) + ncount

                    Svertex(i)%nbrs(new_pos) = inbr
                    Svertex(i)%nbrs_face(new_pos) = Svertex(i)%nbrs_face(pos)
                    Svertex(i)%nbrs_wgt(new_pos) = Svertex(i)%nbrs_wgt(pos)
                    Svertex(i)%nbrs_wgt_ghost(new_pos) = Svertex(i)%nbrs_wgt_ghost(pos)
                    ncount = ncount+1
                 endif
           enddo
           !set neighbors ptr
           orig_start =  Svertex(i)%nbrs_ptr(j+1);
           Svertex(i)%nbrs_ptr(j+1) =  Svertex(i)%nbrs_ptr(j) + ncount
           

           if(Debug) write(iulog,*)'CreateSubGridGraph: point #5.3'
        enddo !num_neighbors loop


        if(Debug) write(iulog,*)'CreateSubGridGraph: point #6'
        Svertex(i)%number = i
     enddo !nelem_s loop
     if(Debug) write(iulog,*)'CreateSubGridGraph: point #7'
     deallocate(global2local)
     if(Debug) write(iulog,*)'CreateSubGridGraph: point #8'

  end subroutine CreateSubGridGraph

!======================================================================

  subroutine PrintGridEdge(Edge)

    implicit none
    type (GridEdge_t), intent(in) :: Edge(:)

    integer           :: i,nedge,ii,wgtP

    nedge = SIZE(Edge)

    write(iulog,95)
    do i=1,nedge
          ii=Edge(i)%tail_face

          !map to correct location - for now all on same nbr side have same wgt, so take the first one
          ii = Edge(i)%tail%nbrs_ptr(ii)

          wgtP=Edge(i)%tail%nbrs_wgt(ii)       
          write(iulog,100) i, &
               Edge(i)%tail%number,Edge(i)%tail_face, wgtP, &
               Edge(i)%head%number,Edge(i)%head_face, gridedge_type(Edge(i))
    enddo
  95 format(5x,'GRIDEDGE #',3x,'Tail (face)',5x,'Head (face)',3x,'Type')
 100 format(10x,I6,8x,I4,1x,'(',I1,')  --',I2,'--> ',I6,1x,'(',I1,')',5x,'[',I1,']')

  end subroutine PrintGridEdge

!======================================================================
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

!======================================================================
  subroutine PrintGridVertex(Vertex)

    implicit none 
    type (GridVertex_t), intent(in),target :: Vertex(:)
  
    integer        :: i,nvert
    integer ::n_west, n_east, n_south, n_north, n_swest, n_seast, n_nwest, n_neast
    integer ::w_west, w_east, w_south, w_north, w_swest, w_seast, w_nwest, w_neast
    integer ::n, print_buf(90), nbr(8), j, k, start, cnt, nbrs_cnt(8)

    integer, pointer :: np(:)

    nbr = (/ west, east, south, north, swest, seast, nwest, neast/)

    nvert = SIZE(Vertex)
        
    write(iulog,98)
    do i=1,nvert

       print_buf(:) = 0
       nbrs_cnt(:) = 0
       np =>  Vertex(i)%nbrs_ptr !alias
       cnt = 1
       do j = 1,num_neighbors
          n = np(nbr(j)+1) - np(nbr(j)) !num neigbors in that directions
          start =  np(nbr(j)) !start in array
          nbrs_cnt(j) = n
          do k = 1, n
             print_buf(cnt) = Vertex(i)%nbrs(start+k-1)
             print_buf(cnt+1) = Vertex(i)%nbrs_wgt(start+k-1)
             print_buf(cnt+2) = Vertex(i)%nbrs_face(start+k-1)
             cnt = cnt + 3
          end do
       enddo

       write(iulog,991) Vertex(i)%number, Vertex(i)%processor_number, &
            Vertex(i)%face_number, &
            print_buf(1:cnt-1)

       write(iulog,992) nbrs_cnt(1:8)


    enddo
  98  format(5x,'GRIDVERTEX #',2x,'PART',2x,'DEG',4x,'W',9x,'E',9x, &
                'S',9x,'N',9x,'SW',9x,'SE',9x,'NW',9x,'NE')

  991  format(10x,I3,8x,I4,8x,I4,2x,30(1x,I4,1x,'(',I2,I2')'))
  992 format(30x,'nbrs_cnt:', 2x,8(1x,I4))

  end subroutine PrintGridVertex


!======================================================================

  subroutine CheckGridNeighbors(Vertex)
  
  implicit none
  type (GridVertex_t), intent(in) :: Vertex(:)

  integer :: i,j,k,l,m,nnbrs,inbrs,nvert
  nvert = SIZE(Vertex)

  do i=1,nvert
        nnbrs = SIZE(Vertex(i)%nbrs)
        do j=1,nnbrs
           inbrs = Vertex(i)%nbrs(j)
           if(inbrs > 0) then
              do k=1,nnbrs
                 if( inbrs .eq. Vertex(i)%nbrs(k) .and. (j/=k) ) &
                      write(iulog,*)'CheckGridNeighbors: ERROR identical neighbors detected  for Vertex ',i
                    
              enddo
           endif
        enddo
     enddo

  end subroutine CheckGridNeighbors

!======================================================================
  subroutine initgridedge(GridEdge,GridVertex)
  use parallel_mod, only : abortmp
  use dimensions_mod, only : max_corner_elem
  implicit none

  type (GridEdge_t), intent(inout)       :: GridEdge(:)
  type (GridVertex_t), intent(in),target :: GridVertex(:)

  integer                                :: i,j,k,iptr,m,n,wgtV,wgtP
  integer                                :: nelem,nelem_edge,inbr
  logical                                :: Verbose=.FALSE.
  integer                                :: mynbr_cnt, cnt, mystart, start

  nelem      = SIZE(GridVertex)
  nelem_edge = SIZE(GridEdge)

  GridEdge(:)%reverse=.FALSE.

  iptr=1
  do j=1,nelem
     do i=1,num_neighbors    
        mynbr_cnt = GridVertex(j)%nbrs_ptr(i+1) - GridVertex(j)%nbrs_ptr(i) !length of neighbor location  
        mystart = GridVertex(j)%nbrs_ptr(i) 
        do m=0,mynbr_cnt-1
           if((GridVertex(j)%nbrs_wgt(mystart + m) .gt. 0)) then    ! Do this only if has a non-zero weight
              if (nelem_edge<iptr) call abortmp('Error in initgridedge: Number of edges greater than expected.')
              GridEdge(iptr)%tail      => GridVertex(j)
              GridEdge(iptr)%tail_face =  mystart + m ! needs to be mystart + m (location in array)
              GridEdge(iptr)%tail_dir = i*max_corner_elem + m !conversion needed for setcycle
              inbr                     =  GridVertex(j)%nbrs(mystart+m)
              GridEdge(iptr)%head      => GridVertex(inbr)

              ! ===========================================
              ! Need this awful piece of code to determine
              ! which "face" of the neighbor element the
              ! edge links (i.e. the "head_face")
              ! ===========================================
              do k=1,num_neighbors
                 cnt = GridVertex(inbr)%nbrs_ptr(k+1) -GridVertex(inbr)%nbrs_ptr(k)                     
                 start = GridVertex(inbr)%nbrs_ptr(k)
                 do  n = 0, cnt-1
                    if(GridVertex(inbr)%nbrs(start+n) == GridVertex(j)%number) then
                       GridEdge(iptr)%head_face=start+n !needs to be start + n (location in array)
                       GridEdge(iptr)%head_dir=k*max_corner_elem+n !conversion (un-done in setcycle)
                    endif
                 enddo
              enddo
              iptr=iptr+1
           end if
        end do ! m loop
     end do !end i loop
  end do !end j loop
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
!======================================================================

end module GridGraph_mod
