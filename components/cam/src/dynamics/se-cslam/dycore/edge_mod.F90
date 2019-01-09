module edge_mod

  use shr_kind_mod,           only: r8=>shr_kind_r8, i8=>shr_kind_i8
  use dimensions_mod,         only: max_neigh_edges, nelemd
  use perf_mod,               only: t_startf, t_stopf, t_adj_detailf ! _EXTERNAL
  use thread_mod,             only: max_num_threads, omp_get_num_threads, omp_get_thread_num
  use coordinate_systems_mod, only: cartesian3D_t
  use schedtype_mod,          only: cycle_t, schedule_t, pgindex_t, schedule, HME_Ordinal,HME_Cardinal
  use cam_abortutils,         only: endrun
  use cam_logfile,            only: iulog
  use parallel_mod,           only: parallel_t, &
      MAX_ACTIVE_MSG, HME_status_size, BNDRY_TAG_BASE, HME_BNDRY_A2A, HME_BNDRY_P2P, &
      HME_BNDRY_A2AO
  use edgetype_mod,           only: edgedescriptor_t, edgebuffer_t,           &
      Longedgebuffer_t, initedgebuffer_callid, Ghostbuffer3D_t
  use element_mod,            only: element_t
  use gbarrier_mod,           only: gbarrier_init, gbarrier_delete
  use spmd_utils,             only: mpi_real8, mpi_integer, mpi_info_null, mpi_success

  implicit none
  private
  save

  ! 8-byte Integer routines
  public :: LongEdgeVpack, LongEdgeVunpackMIN

  ! 8-byte Real routines
  public :: zeroEdgeBuffer

  interface initEdgeBuffer
     module procedure initEdgeBuffer_r8
     module procedure initEdgeBuffer_i8
  end interface
  interface initEdgeSBuffer
     module procedure initEdgeSbuffer_r8
  end interface
  interface freeEdgeBuffer
     module procedure freeEdgeBuffer_r8
     module procedure freeEdgeBuffer_i8
  end interface
  interface freeGhostBuffer
     module procedure freeGhostBuffer_r8
  end interface

  public :: initEdgeBuffer
  public :: initEdgeSBuffer
  public :: freeEdgeBuffer

  public :: initGhostBuffer
  public :: ghostpack, ghostunpack
  public :: freeGhostBuffer

  !---------------------------------------------------------
  ! Pack/unpack routines that use the New format Edge buffer
  !---------------------------------------------------------

  public :: edgeVpack, edgeVunpack
  public :: edgeVunpackMIN, edgeVunpackMAX
  public :: edgeDGVpack, edgeDGVunpack
  public :: edgeVunpackVert


  public :: initGhostBuffer3D
  public :: FreeGhostBuffer3D
  public :: ghostVpack3D, ghostVunpack3D

  !----------------------------------------------------------------
  ! Pack/unpack routines that communicate a fixed number values
  ! per element.  This is used to communicate MIN/MAX values from
  ! neighboring elemeents
  !----------------------------------------------------------------
  interface edgeSpack
     module procedure edgeSpack_r8
  end interface
  public :: edgeSpack
  public :: edgeSunpackMIN, edgeSunpackMAX

  logical, private :: threadsafe=.true.

  real(kind=r8), parameter, public :: edgeDefaultVal = 1.11e+100_r8

! NOTE ON ELEMENT ORIENTATION
!
! Element orientation:  index V(i,j)
!
!           (1,np) NWEST      (np,np) NEAST
!
!           (1,1) SWEST       (np,1) SEAST
!
!
! for the edge neighbors:
!    we set the "reverse" flag if two elements who share an edge use a
!    reverse orientation.  The data is reversed during the *pack* stage
! For corner neighbors:
!    for edge buffers, there is no orientation because two corner neighbors
!    only share a single point.
!    For ghost cell data, there is a again two posible orientations. For
!    this case, we set the "reverse" flag if the corner element is using
!    the reverse orientation.  In this case, the data is reversed during the
!    *unpack* stage (not sure why)
!
! The edge orientation is set at startup.  The corner orientation is computed
! at run time, via the call to compute_ghost_corner_orientation()
! This routine only works for meshes with at most 1 corner element.  It's
! not called and the corner orientation flag is not set for unstructured meshes

!
! Christoph Erath
! pack/unpack partial element of data of size (nx,nx) with user specifed halo size nh
! user specifies the sizes when creating the buffer
! buffer has 1 extra dimension (as compared to subroutines above) for multiple tracers
! input/output arrays are cartesian, and thus assume at most 1 element at each corner
! hence currently only supports cube-sphere grids.
!
!
! routines which including element edge data
! (used for FVM arrays where edge data is not shared by neighboring elements)
! these routines pack/unpack element data with user specified halo size
                                  
  ! Wrap pointer so we can make an array of them.
  type :: wrap_ptr
    real (kind=r8), dimension(:,:), pointer :: ptr => null()
  end type wrap_ptr

  type(wrap_ptr) :: edgebuff_ptrs(0:1)

contains

  subroutine initEdgeSBuffer_r8(par,edge,elem,nlyr,bndry_type, nthreads)
    type (parallel_t),             intent(in)  :: par
    type (EdgeBuffer_t), target,   intent(out) :: edge
    type (element_t),              intent(in)  :: elem(:)
    integer,                       intent(in)  :: nlyr
    integer ,            optional, intent(in)  :: bndry_type
    integer,             optional, intent(in)  :: nthreads


    call initEdgeBuffer(par,edge,elem,nlyr,bndry_type=bndry_type, &
                         nthreads=nthreads,CardinalLength=1,OrdinalLength=1)

  end subroutine initEdgeSBuffer_r8

  subroutine initGhostBuffer(par,edge,elem,nlyr,ndepth, npoints,bndry_type,nthreads)

    type (parallel_t),             intent(in)  :: par
    type (Edgebuffer_t), target,   intent(out) :: edge
    type (element_t),              intent(in)  :: elem(:)
    integer,intent(in)                         :: nlyr,ndepth, npoints
    integer ,            optional, intent(in)  :: bndry_type
    integer,             optional, intent(in)  :: nthreads

    call initEdgeBuffer(par,edge,elem,nlyr,bndry_type=bndry_type, &
                         nthreads=nthreads,CardinalLength=ndepth*npoints,OrdinalLength=ndepth*ndepth)
    ! set some parameters need to support deep halos
    edge%ndepth  = ndepth 
    edge%npoints = npoints
    edge%lb      = 1 - edge%ndepth
    edge%ub      = edge%npoints + edge%ndepth

  end subroutine initGhostBuffer
   


  subroutine zeroEdgeBuffer(edge)

    type (EdgeBuffer_t), intent(inout) :: edge
    integer :: i

    do i=1,edge%nbuf
       edge%buf(i)     = 0.0d0
       edge%receive(i) = 0.0d0
    enddo

  end subroutine zeroEdgeBuffer

  ! =========================================
  ! initEdgeBuffer:
  !
  ! create an Real based communication buffer
  ! =========================================
  subroutine initEdgeBuffer_r8(par,edge,elem,nlyr, bndry_type,nthreads,CardinalLength,OrdinalLength)
    use dimensions_mod, only: np, nelemd, max_corner_elem
    use schedtype_mod,  only: cycle_t, schedule_t, schedule
    use mpi,            only: MPI_VERSION

    type (parallel_t),             intent(in)  :: par
    type (EdgeBuffer_t), target,   intent(out) :: edge
    type (element_t),              intent(in)  :: elem(:)
    integer,                       intent(in)  :: nlyr
    integer,             optional, intent(in)  :: bndry_type 
    integer,             optional, intent(in)  :: nthreads
    integer,             optional, intent(in)  :: CardinalLength
    integer,             optional, intent(in)  :: OrdinalLength

    ! Notes about the buf_ptr/receive_ptr options:
    !
    ! If an EdgeBuffer_t object is initialized from pre-existing storage
    ! (i.e. buf_ptr is provided and not null), it must *not* be freed,
    ! and must not be used if the underlying storage has been deallocated.
    !
    ! All these restrictions also applied to the old newbuf and newreceive
    ! options.

    ! Workaround for NAG bug.
    ! NAG 5.3.1 dies if you use pointer bounds remapping to set
    ! a pointer that is also a component. So remap to temporary,
    ! then use that to set component pointer.

    ! Local variables
    integer                    :: nbuf,ith
    integer                    :: nSendCycles, nRecvCycles
    integer                    :: icycle, ierr
    integer                    :: ie, i
    integer                    :: edgeid,elemid
    integer                    :: ptr,llen,moveLength, mLen, tlen
    type (Cycle_t),    pointer :: pCycle
    type (Schedule_t), pointer :: pSchedule
    integer                    :: dest, source, length, tag, iptr
    integer                    :: nlen, ithr

    integer :: len, lenP,lenS
    integer :: j,jj,il,mesgid, dst0,src0
    integer :: moveptr
    integer :: nbuf2,ilm1,iem1,lenm1
    integer,allocatable :: putmap2(:,:),getmap2(:,:)
    integer,allocatable :: scounts(:), rcounts(:)
    integer,allocatable :: sdispls(:), rdispls(:)
    integer :: nInter, nIntra
    integer :: icInter, icIntra

    integer :: maxnsend
    integer :: tmpnMesg
    integer :: wintmpnMesg, wintmpDest, wintmpDisp
    integer(kind=i8) :: winsize
    integer :: win
    integer :: sizeofreal
    integer, allocatable :: tmpDest(:),tmpDisp(:)
    integer :: nFull
    integer :: disp, one
    integer                           :: errorcode,errorlen
    integer :: CardinalLen, OrdinalLen
    character(len=80)                 :: errorstring
    character(len=80), parameter      :: subname='initedgeBuffer'

    if(present(bndry_type)) then 
      if ( MPI_VERSION >= 3 ) then
        edge%bndry_type = bndry_type
      else
        edge%bndry_type = HME_BNDRY_P2P
      endif
    else
       edge%bndry_type = HME_BNDRY_P2P
    endif

    ! Set the length of the cardinal and ordinal message lengths
    if(present(CardinalLength)) then 
       CardinalLen = CardinalLength
    else
       CardinalLen = np
    endif
    if(present(OrdinalLength)) then 
       OrdinalLen = OrdinalLength
    else
       OrdinalLen = 1
    endif

! DO NOT REMOVE THIS NEXT BARRIER
! MT: This initial barrier fixes a long standing issue with Intel compilers on
! two different platforms.  Without this barrier, edge buffers initialized from
! within the threaded region would not work in a reproducable way with certain
! thread combinations.  I cant explain why, but this fixes that issue on Edison
!$OMP BARRIER

    if (nlyr==0) return  ! tracer code might call initedgebuffer() with zero tracers


!$OMP MASTER
    !
    ! Keep a counter of how many times initedgebuffer is called.
    ! This is used to assign a unique message ID for the boundary exchange
    !
    initedgebuffer_callid=initedgebuffer_callid+1
    edge%id  = initedgebuffer_callid
    edge%tag = BNDRY_TAG_BASE + MODULO(edge%id, MAX_ACTIVE_MSG)

    allocate(edge%putmap(max_neigh_edges,nelemd))
    allocate(edge%getmap(max_neigh_edges,nelemd))
    allocate(edge%reverse(max_neigh_edges,nelemd))

    edge%putmap(:,:)=-1
    edge%getmap(:,:)=-1

    allocate(putmap2(max_neigh_edges,nelemd))
    allocate(getmap2(max_neigh_edges,nelemd))
    putmap2(:,:)=-1
    getmap2(:,:)=-1
    do ie=1,nelemd
       do i=1,max_neigh_edges
          edge%reverse(i,ie) = elem(ie)%desc%reverse(i)
       enddo
    enddo

    pSchedule  => Schedule(1)
    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles
    nInter      = pSchedule%nInter
    nIntra      = pSchedule%nIntra
    nFull       = nInter+nIntra

    edge%nInter=nInter
    edge%nIntra=nIntra

    if(nInter>0) then
       allocate(edge%rcountsInter(nInter),edge%rdisplsInter(nInter))
       allocate(edge%scountsInter(nInter),edge%sdisplsInter(nInter))
    endif
    if(nIntra>0) then
       allocate(edge%rcountsIntra(nIntra),edge%rdisplsIntra(nIntra))
       allocate(edge%scountsIntra(nIntra),edge%sdisplsIntra(nIntra))
    endif

    if (nSendCycles>0) then
       allocate(edge%scountsFull(nSendCycles),edge%sdisplsFull(nSendCycles))
       allocate(edge%Srequest(nSendCycles))
       edge%scountsFull(:) = 0
    endif
    !
    ! Setup the data-structures for the sends
    !
    j         = 1
    icycle    = 1
    dst0      = pSchedule%pIndx(j)%mesgid
    il        = pSchedule%pIndx(j)%edgeid
    ie        = pSchedule%pIndx(j)%elemid
    len       = CalcSegmentLength(pSchedule%pIndx(j),CardinalLen,OrdinalLen,nlyr)
    edge%putmap(il,ie) = 0
    if(nSendCycles>0) then 
        edge%sdisplsFull(icycle) = edge%putmap(il,ie)
        edge%scountsFull(icycle) = len
    endif
    ilm1 = il
    iem1 = ie
    lenm1 = len

    do j=2,SIZE(pSchedule%pIndx)
       il = pSchedule%pIndx(j)%edgeid
       ie = pSchedule%pIndx(j)%elemid
       mesgid = pSchedule%pIndx(j)%mesgid
       if(il>0 .and. ie >0) then
          len     = CalcSegmentLength(pSchedule%pIndx(j),CardinalLen,OrdinalLen,nlyr)
          edge%putmap(il,ie) = edge%putmap(ilm1,iem1)+lenm1
          if(mesgid .ne. par%rank) then  ! don't enter if this is a move cycle where (mesgid == par%rank)
             if(mesgid .ne. dst0) then
                icycle=icycle+1
                if (nSendCycles>0) edge%sdisplsFull(icycle) = edge%putmap(il,ie)
                dst0=mesgid
             endif
             if (nSendCycles>0) edge%scountsFull(icycle) = edge%scountsFull(icycle)+len
          endif
          ilm1=il
          iem1=ie
          lenm1=len
       endif
    enddo

    icInter=0
    icIntra=0
    do icycle=1,nSendCycles
       if(pSchedule%SendCycle(icycle)%onNode .eqv. .FALSE.) then
          icInter=icInter+1
          edge%sdisplsInter(icInter)=edge%sdisplsFull(icycle)
          edge%scountsInter(icInter)=edge%scountsFull(icycle)
       else
          icIntra=icIntra+1
          edge%sdisplsIntra(icIntra)=edge%sdisplsFull(icycle)
          edge%scountsIntra(icIntra)=edge%scountsFull(icycle)
       endif
    enddo

    if (nRecvCycles>0) then
       allocate(edge%rcountsFull(nRecvCycles),edge%rdisplsFull(nRecvCycles))
       allocate(edge%getDisplsFull(nRecvCycles),edge%putDisplsFull(nRecvCycles))
       edge%rcountsFull(:) = 0
       ! allocate the MPI Send/Recv request handles
       allocate(edge%Rrequest(nRecvCycles))
       allocate(edge%status(HME_status_size,nRecvCycles))
    endif

    !
    ! Setup the data-structures for the receives
    !
    j      = 1
    icycle = 1
    src0    = pSchedule%gIndx(j)%mesgid
    il      = pSchedule%gIndx(j)%edgeid
    ie      = pSchedule%gIndx(j)%elemid
    len     = CalcSegmentLength(pSchedule%gIndx(j),CardinalLen,OrdinalLen,nlyr)
    edge%getmap(il,ie) = 0
    if (nRecvCycles>0) then
       edge%rdisplsFull(icycle) = edge%getmap(il,ie)
       edge%rcountsFull(icycle) = len
    endif
    ilm1=il
    iem1=ie
    lenm1=len

    do j=2,SIZE(pSchedule%gIndx)
       il = pSchedule%gIndx(j)%edgeid
       ie = pSchedule%gIndx(j)%elemid
       mesgid = pSchedule%gIndx(j)%mesgid
       if(il>0 .and. ie >0) then
          len     = CalcSegmentLength(pSchedule%gIndx(j),CardinalLen,OrdinalLen,nlyr)
          edge%getmap(il,ie) = edge%getmap(ilm1,iem1)+lenm1
          if(mesgid .ne. par%rank) then ! don't enter if this is a move cycle where (mesgid == par%rank)
             if(mesgid .ne. src0) then
                if (nRecvCycles>0) edge%rdisplsFull(icycle+1) = edge%getmap(il,ie)
                icycle=icycle+1
                src0=mesgid
             endif
             if (nRecvCycles>0) edge%rcountsFull(icycle) = edge%rcountsFull(icycle)+len
          endif
          ilm1=il
          iem1=ie
          lenm1=len
       endif
    enddo


    !
    ! populate the Inter and Intra node communication data-structures
    !
    icInter=0
    icIntra=0
    do icycle=1,nRecvCycles
       if(pSchedule%RecvCycle(icycle)%onNode .eqv. .FALSE.) then
          icInter=icInter+1
          edge%rdisplsInter(icInter)=edge%rdisplsFull(icycle)
          edge%rcountsInter(icInter)=edge%rcountsFull(icycle)
       else
          icIntra=icIntra+1
          edge%rdisplsIntra(icIntra)=edge%rdisplsFull(icycle)
          edge%rcountsIntra(icIntra)=edge%rcountsFull(icycle)
       endif
    enddo


    ! Setup the data-structures for the on process moves
    ! Note that this assumes that the data to move is at
    ! the end of the  message buffer.
    if(nRecvCycles>0) then
       moveptr = edge%rdisplsFull(nRecvCycles)+edge%rcountsFull(nRecvCycles)+1
    else
       moveptr = 1
    endif
    moveLength = 0
    do j=1,SIZE(pSchedule%gIndx)
       il     = pSchedule%gIndx(j)%edgeid
       ie     = pSchedule%gIndx(j)%elemid
       mesgid = pSchedule%gIndx(j)%mesgid
       if(mesgid == par%rank) then
          len     = CalcSegmentLength(pSchedule%gIndx(j),CardinalLen,OrdinalLen,nlyr)
          moveLength = moveLength + len
       endif
    enddo

    ! decompose the move data between the available threads
    if(max_num_threads<=0) then
       nlen = 1
    else
       if(present(nthreads)) then
          if (nthreads > 0) then
             nlen = nthreads
          else
             nlen = max_num_threads
          end if
       else
          nlen = max_num_threads
       end if
    end if
    call gbarrier_init(edge%gbarrier, nlen)

    allocate(edge%moveLength(nlen))
    allocate(edge%movePtr(nlen))

    if (nlen > 1) then
       ! the master thread performs no data movement because it is busy with the
       ! MPI messaging
       edge%moveLength(1) = -1
       edge%movePtr(1) = 0

       ! Calculate the length of the local copy in bndy_exchange
       llen = ceiling(real(moveLength,kind=r8)/real(nlen-1,kind=r8))
       iptr = moveptr
       mLen = 0
       do i=2,nlen
         if( (mLen+llen) <= moveLength)  then
            tlen = llen
         else
            tlen = moveLength - mLen
         endif
         edge%moveLength(i) = tlen
         edge%movePtr(i)    = iptr
         iptr = iptr + tlen
         mLen = mLen + tLen
       enddo
    else
       edge%moveLength(1) = moveLength
       edge%movePtr(1) = moveptr
    endif

    ! Set the maximum length of the message buffer
    nbuf = movePtr+moveLength

    edge%nlyr=nlyr
    edge%nbuf=nbuf

    allocate(edge%receive(nbuf))
    allocate(edge%buf(nbuf))

21  format('RANK: ',i2, A,8(i6))

!$OMP END MASTER
! MT: This next barrier is also needed - threads cannot start using edge()
! until MASTER is done initializing it
!$OMP BARRIER

  end subroutine initEdgeBuffer_r8

  integer function CalcSegmentLength(pgIndx,CardinalLength,OrdinalLength,nlyr) result(len)

     type(pgindex_t) ::  pgIndx
     integer, intent(in) :: CardinalLength,OrdinalLength
     integer, intent(in) :: nlyr

     integer :: rem
     integer, parameter :: alignment=1  ! align on word boundaries
!     integer, parameter :: alignment=2  ! align on 2 word boundaries
!     integer, parameter :: alignment=8  ! align on 8 word boundaries

     select case(pgIndx%edgeType)
        CASE(HME_Cardinal)
          len = nlyr*CardinalLength
        CASE(HME_Ordinal)
          len = nlyr*OrdinalLength
     end select

     rem = MODULO(len,alignment)
     if(rem .ne. 0) then
        len = len + (alignment-rem)
     endif

  end function calcSegmentLength

  ! =========================================
  ! initEdgeBuffer:
  !
  ! create an Integer based communication buffer
  ! =========================================
  subroutine initEdgeBuffer_i8(edge,nlyr)
    use dimensions_mod, only : np, nelemd, max_corner_elem

    integer,                 intent(in)  :: nlyr
    type (LongEdgeBuffer_t), intent(out) :: edge

    ! Local variables
    integer :: nbuf

    ! sanity check for threading
    if (omp_get_num_threads()>1) then
       call endrun('ERROR: initEdgeBuffer must be called before threaded reagion')
    endif

    nbuf=4*(np+max_corner_elem)*nelemd
    edge%nlyr=nlyr
    edge%nbuf=nbuf
    allocate(edge%buf(nlyr,nbuf))
    edge%buf(:,:)=0

    allocate(edge%receive(nlyr,nbuf))
    edge%receive(:,:)=0

  end subroutine initEdgeBuffer_i8
  ! =========================================
  ! edgeDGVpack:
  !
  ! Pack edges of v into buf for DG stencil
  ! =========================================
  subroutine edgeDGVpack(edge,v,vlyr,kptr,ielem)
    use dimensions_mod, only: np

    type (EdgeBuffer_t)        :: edge
    integer,        intent(in) :: vlyr
    real (kind=r8), intent(in) :: v(np,np,vlyr)
    integer,        intent(in) :: kptr
    integer,        intent(in) :: ielem

    ! =========================================
    ! This code is just a wrapper call the
    !   normal oldedgeVpack
    ! =========================================
    call edgeVpack(edge,v,vlyr,kptr,ielem)

  end subroutine edgeDGVpack

  subroutine FreeGhostBuffer_r8(edge)
     type (EdgeBuffer_t), intent(inout) :: edge
     call FreeEdgeBuffer_r8(edge)
  end subroutine FreeGhostBuffer_r8
  ! ===========================================
  !  FreeEdgeBuffer:
  !
  !  Freed an edge communication buffer
  ! =========================================
  subroutine FreeEdgeBuffer_r8(edge)

    type (EdgeBuffer_t),intent(inout) :: edge

!$OMP BARRIER
!$OMP MASTER
    deallocate(edge%buf)
    deallocate(edge%receive)
    if(associated(edge%putmap))     deallocate(edge%putmap)
    if(associated(edge%getmap))     deallocate(edge%getmap)
    if(associated(edge%reverse))    deallocate(edge%reverse)
    if(associated(edge%moveLength)) deallocate(edge%moveLength)
    if(associated(edge%movePtr))    deallocate(edge%movePtr)

    ! All MPI communications
    if(associated(edge%rcountsFull)) deallocate(edge%rcountsFull)
    if(associated(edge%scountsFull)) deallocate(edge%scountsFull)
    if(associated(edge%sdisplsFull)) deallocate(edge%sdisplsFull)
    if(associated(edge%rdisplsFull)) deallocate(edge%rdisplsFull)

    ! Intra-node MPI Communication
    if(edge%nIntra>0) then
      if(associated(edge%rcountsIntra)) deallocate(edge%rcountsIntra)
      if(associated(edge%scountsIntra)) deallocate(edge%scountsIntra)
      if(associated(edge%sdisplsIntra)) deallocate(edge%sdisplsIntra)
      if(associated(edge%rdisplsIntra)) deallocate(edge%rdisplsIntra)
    endif

    ! Inter-node MPI Communication
    if(edge%nInter>0) then
      if(associated(edge%rcountsInter)) deallocate(edge%rcountsInter)
      if(associated(edge%scountsInter)) deallocate(edge%scountsInter)
      if(associated(edge%sdisplsInter)) deallocate(edge%sdisplsInter)
      if(associated(edge%rdisplsInter)) deallocate(edge%rdisplsInter)
    endif
    if(allocated(edge%rRequest)) deallocate(edge%rRequest)
    if(allocated(edge%sRequest)) deallocate(edge%sRequest)
    if(allocated(edge%status)) deallocate(edge%status)
    call gbarrier_delete(edge%gbarrier)

!$OMP END MASTER

  end subroutine FreeEdgeBuffer_r8

  ! ===========================================
  !  FreeEdgeBuffer:
  !
  !  Freed an edge communication buffer
  ! =========================================
  subroutine FreeEdgeBuffer_i8(edge)

    type (LongEdgeBuffer_t),intent(inout) :: edge

    edge%nbuf=0
    edge%nlyr=0
    deallocate(edge%buf)
    deallocate(edge%receive)

  end subroutine FreeEdgeBuffer_i8

  ! =========================================
  !
  !> @brief Pack edges of v into an edge buffer for boundary exchange.
  !
  !> This subroutine packs for one or more vertical layers into an edge
  !! buffer. If the buffer associated with edge is not large enough to
  !! hold all vertical layers you intent to pack, the method will
  !! halt the program with a call to endrum().
  !! @param[in] edge Edge Buffer into which the data will be packed.
  !! This buffer must be previously allocated with initEdgeBuffer().
  !! @param[in] v The data to be packed.
  !! @param[in] vlyr Number of vertical level coming into the subroutine
  !! for packing for input v.
  !! @param[in] kptr Vertical pointer to the place in the edge buffer where
  !! data will be located.
  ! =========================================
  subroutine edgeVpack(edge,v,vlyr,kptr,ielem)
    use dimensions_mod, only: np, max_corner_elem
    use control_mod,    only: north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t)        :: edge
    integer,        intent(in) :: vlyr
    real (kind=r8), intent(in) :: v(np,np,vlyr)
    integer,        intent(in) :: kptr
    integer,        intent(in) :: ielem

    ! Local variables
    integer :: i,k,ir,ll,iptr
    integer :: is,ie,in,iw,edgeptr

    is = edge%putmap(south,ielem)
    ie = edge%putmap(east,ielem)
    in = edge%putmap(north,ielem)
    iw = edge%putmap(west,ielem)
    if (edge%nlyr < (kptr+vlyr) ) then
       print *,'edge%nlyr = ',edge%nlyr
       print *,'kptr+vlyr = ',kptr+vlyr
       call endrun('edgeVpack: Buffer overflow: size of the vertical dimension must be increased!')
    endif

!dir$ ivdep
    do k=1,vlyr
       iptr = np*(kptr+k-1)
       do i=1,np
          edge%buf(iptr+ie+i)   = v(np ,i ,k) ! East
          edge%buf(iptr+is+i)   = v(i  ,1 ,k) ! South
          edge%buf(iptr+in+i)   = v(i  ,np,k) ! North
          edge%buf(iptr+iw+i)   = v(1  ,i ,k) ! West
       enddo
    enddo

    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing

    if(edge%reverse(south,ielem)) then
       do k=1,vlyr
          iptr = np*(kptr+k-1)+is
          do i=1,np
             ir = np-i+1
             edge%buf(iptr+ir)=v(i,1,k)
          enddo
       enddo
    endif

    if(edge%reverse(east,ielem)) then
       do k=1,vlyr
          iptr=np*(kptr+k-1)+ie
          do i=1,np
             ir = np-i+1
             edge%buf(iptr+ir)=v(np,i,k)
          enddo
       enddo
    endif

    if(edge%reverse(north,ielem)) then
       do k=1,vlyr
          iptr=np*(kptr+k-1)+in
          do i=1,np
             ir = np-i+1
             edge%buf(iptr+ir)=v(i,np,k)
          enddo
       enddo
    endif

    if(edge%reverse(west,ielem)) then
       do k=1,vlyr
          iptr=np*(kptr+k-1)+iw
          do i=1,np
             ir = np-i+1
             edge%buf(iptr+ir)=v(1,i,k)
          enddo
       enddo
    endif

! SWEST
    do ll=swest,swest+max_corner_elem-1
        if (edge%putmap(ll,ielem) /= -1) then
            edgeptr = edge%putmap(ll,ielem)+1
            do k=1,vlyr
                iptr = (kptr+k-1)+edgeptr
                if (iptr > size(edge%buf)) then
                   write(6, *) 'ERROR SW: ',size(edge%buf),iptr,edge%putmap(ll,ielem)
                   call endrun('pointer bounds ERROR SW')
                end if
                edge%buf(iptr) = v(1, 1, k)
            end do
        end if
    end do

! SEAST
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
        if (edge%putmap(ll,ielem) /= -1) then
            edgeptr = edge%putmap(ll,ielem)+1
            do k=1,vlyr
               iptr = (kptr+k-1)+edgeptr
               if (iptr > size(edge%buf)) then
                  write(6, *) 'ERROR SE: ',size(edge%buf),iptr,edge%putmap(ll,ielem)
                  call endrun('pointer bounds ERROR SE')
               end if
               edge%buf(iptr)=v(np, 1, k)
            end do
        end if
    end do

! NEAST
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        if (edge%putmap(ll,ielem) /= -1) then
            edgeptr = edge%putmap(ll,ielem)+1
            do k=1,vlyr
               iptr = (kptr+k-1)+edgeptr
               if (iptr > size(edge%buf)) then
                  write(6, *) 'ERROR NE: ',size(edge%buf),iptr,edge%putmap(ll,ielem)
                  call endrun('pointer bounds ERROR NE')
               end if
               edge%buf(iptr) = v(np, np, k)
            end do
        end if
    end do

! NWEST
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        if (edge%putmap(ll,ielem) /= -1) then
            edgeptr = edge%putmap(ll,ielem)+1
            do k=1,vlyr
               iptr = (kptr+k-1)+edgeptr
               if (iptr > size(edge%buf)) then
                  write(6, *) 'ERROR NW: ',size(edge%buf),iptr,edge%putmap(ll,ielem)
                  call endrun('pointer bounds ERROR NW')
               end if
               edge%buf(iptr) = v(1, np, k)
            end do
        end if
    end do

  end subroutine edgeVpack

  subroutine edgeSpack_r8(edge,v,vlyr,kptr,ielem)
    use dimensions_mod, only: np, max_corner_elem
    use control_mod,    only: north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t)        :: edge
    integer,        intent(in) :: vlyr
    real (kind=r8), intent(in) :: v(vlyr)
    integer,        intent(in) :: kptr
    integer,        intent(in) :: ielem

    ! Local variables
    integer        :: i,k,ir,ll,iptr
    integer        :: is,ie,in,iw,edgeptr
    real (kind=r8) :: tmp

    is = edge%putmap(south,ielem)
    ie = edge%putmap(east,ielem)
    in = edge%putmap(north,ielem)
    iw = edge%putmap(west,ielem)
    if (edge%nlyr < (kptr+vlyr) ) then
       call endrun('edgeSpack: Buffer overflow: size of the vertical dimension must be increased!')
    endif

    do k=1,vlyr
       iptr = kptr+k-1
       edge%buf(iptr+ie+1)   = v(k) ! East
       edge%buf(iptr+is+1)   = v(k) ! South
       edge%buf(iptr+in+1)   = v(k) ! North
       edge%buf(iptr+iw+1)   = v(k) ! West
    enddo

! SWEST
    do ll=swest,swest+max_corner_elem-1
        if (edge%putmap(ll,ielem) /= -1) then
            edgeptr=edge%putmap(ll,ielem)+1
            do k=1,vlyr
                iptr = (kptr+k-1)+edgeptr
                edge%buf(iptr)=v(k)
            end do
        end if
    end do

! SEAST
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
        if (edge%putmap(ll,ielem) /= -1) then
            edgeptr=edge%putmap(ll,ielem)+1
            do k=1,vlyr
                iptr = (kptr+k-1)+edgeptr
                edge%buf(iptr)=v(k)
            end do
        end if
    end do

! NEAST
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        if (edge%putmap(ll,ielem) /= -1) then
            edgeptr=edge%putmap(ll,ielem)+1
            do k=1,vlyr
                iptr = (kptr+k-1)+edgeptr
                edge%buf(iptr)=v(k)
            end do
        end if
    end do

! NWEST
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        if (edge%putmap(ll,ielem) /= -1) then
            edgeptr=edge%putmap(ll,ielem)+1
            do k=1,vlyr
                iptr = (kptr+k-1)+edgeptr
                edge%buf(iptr)=v(k)
            end do
        end if
    end do

  end subroutine edgeSpack_r8

  ! =========================================
  ! LongEdgeVpack:
  !
  ! Pack edges of v into buf...
  ! =========================================
  subroutine LongEdgeVpack(edge,v,vlyr,kptr,desc)
    use control_mod,    only: north, south, east, west, neast, nwest, seast, swest
    use dimensions_mod, only: np, max_corner_elem

    type (LongEdgeBuffer_t)             :: edge
    integer,                 intent(in) :: vlyr
    integer ,                intent(in) :: v(np,np,vlyr)
    integer,                 intent(in) :: kptr
    type (EdgeDescriptor_t), intent(in) :: desc

    ! Local variables
    logical, parameter :: UseUnroll = .TRUE.
    integer            :: i,k,ir,l
    integer            :: is,ie,in,iw

    if(.not. threadsafe) then
!$OMP BARRIER
       threadsafe=.true.
    end if

    is = desc%putmapP(south)
    ie = desc%putmapP(east)
    in = desc%putmapP(north)
    iw = desc%putmapP(west)

    if(MODULO(np,2) == 0 .and. UseUnroll) then
       do k=1,vlyr
          do i=1,np,2
             edge%buf(kptr+k,is+i)   = v(i  ,1 ,k)
             edge%buf(kptr+k,is+i+1) = v(i+1,1 ,k)
             edge%buf(kptr+k,ie+i)   = v(np ,i ,k)
             edge%buf(kptr+k,ie+i+1) = v(np ,i+1 ,k)
             edge%buf(kptr+k,in+i)   = v(i  ,np,k)
             edge%buf(kptr+k,in+i+1) = v(i+1  ,np,k)
             edge%buf(kptr+k,iw+i)   = v(1  ,i ,k)
             edge%buf(kptr+k,iw+i+1) = v(1  ,i+1 ,k)

          enddo
       end do
    else
       do k=1,vlyr
          do i=1,np
             edge%buf(kptr+k,is+i)   = v(i  ,1 ,k)
             edge%buf(kptr+k,ie+i)   = v(np ,i ,k)
             edge%buf(kptr+k,in+i)   = v(i  ,np,k)
             edge%buf(kptr+k,iw+i)   = v(1  ,i ,k)
          enddo
       end do

    endif


    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing

    if(desc%reverse(south)) then
       is = desc%putmapP(south)
       do k=1,vlyr
          do i=1,np
             ir = np-i+1
             edge%buf(kptr+k,is+ir)=v(i,1,k)
          enddo
       enddo
    endif

    if(desc%reverse(east)) then
       ie = desc%putmapP(east)
       do k=1,vlyr
          do i=1,np
             ir = np-i+1
             edge%buf(kptr+k,ie+ir)=v(np,i,k)
          enddo
       enddo
    endif

    if(desc%reverse(north)) then
       in = desc%putmapP(north)
       do k=1,vlyr
          do i=1,np
             ir = np-i+1
             edge%buf(kptr+k,in+ir)=v(i,np,k)
          enddo
       enddo
    endif

    if(desc%reverse(west)) then
       iw = desc%putmapP(west)
       do k=1,vlyr
          do i=1,np
             ir = np-i+1
             edge%buf(kptr+k,iw+ir)=v(1,i,k)
          enddo
       enddo
    endif

! SWEST
    do l=swest,swest+max_corner_elem-1
        if (desc%putmapP(l) /= -1) then
            do k=1,vlyr
                edge%buf(kptr+k,desc%putmapP(l)+1)=v(1  ,1 ,k)
            end do
        end if
    end do

! SEAST
    do l=swest+max_corner_elem,swest+2*max_corner_elem-1
        if (desc%putmapP(l) /= -1) then
            do k=1,vlyr
                edge%buf(kptr+k,desc%putmapP(l)+1)=v(np ,1 ,k)
            end do
        end if
    end do

! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        if (desc%putmapP(l) /= -1) then
            do k=1,vlyr
                edge%buf(kptr+k,desc%putmapP(l)+1)=v(np ,np,k)
            end do
        end if
    end do

! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        if (desc%putmapP(l) /= -1) then
            do k=1,vlyr
                edge%buf(kptr+k,desc%putmapP(l)+1)=v(1  ,np,k)
            end do
        end if
    end do

  end subroutine LongEdgeVpack

  subroutine edgeVunpack(edge,v,vlyr,kptr,ielem,rank)
    use dimensions_mod, only: np, max_corner_elem
    use control_mod,    only: north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t), intent(in)    :: edge
    integer,             intent(in)    :: vlyr
    real (kind=r8),      intent(inout) :: v(np,np,vlyr)
    integer,             intent(in)    :: kptr
    integer,             intent(in)    :: ielem
    integer, optional,   intent(in)    :: rank

    ! Local
    integer :: i,k,ll,iptr
    integer :: is,ie,in,iw,edgeptr
    integer :: ise,isw,ine,inw
    integer :: ks,ke,kblock
    logical :: done

    is=edge%getmap(south,ielem)
    ie=edge%getmap(east,ielem)
    in=edge%getmap(north,ielem)
    iw=edge%getmap(west,ielem)
    isw=edge%getmap(swest,ielem)
    ise=edge%getmap(seast,ielem)
    inw=edge%getmap(nwest,ielem)
    ine=edge%getmap(neast,ielem)

    !DIR$ IVDEP
    do k=1,vlyr
       iptr=np*(kptr+k-1)
       do i=1,np
          v(np ,i  ,k) = v(np ,i  ,k)+edge%receive(iptr+i+ie) ! East
          v(i  ,1  ,k) = v(i  ,1  ,k)+edge%receive(iptr+i+is) ! South
          v(i  ,np ,k) = v(i  ,np ,k)+edge%receive(iptr+i+in) ! North
          v(1  ,i  ,k) = v(1  ,i  ,k)+edge%receive(iptr+i+iw) ! West
       enddo
    enddo

! SWEST
    do ll=swest,swest+max_corner_elem-1
        if(edge%getmap(ll,ielem) /= -1) then
            edgeptr=edge%getmap(ll,ielem)+1
            do k=1,vlyr
                iptr = (kptr+k-1)+edgeptr
                v(1  ,1 ,k)=v(1 ,1 ,k)+edge%receive(iptr)
            enddo
        endif
    end do

! SEAST
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
        if(edge%getmap(ll,ielem) /= -1) then
            edgeptr=edge%getmap(ll,ielem)+1
            do k=1,vlyr
                iptr = (kptr+k-1)+edgeptr
                v(np ,1 ,k)=v(np,1 ,k)+edge%receive(iptr)
            enddo
        endif
    end do

! NEAST
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        if(edge%getmap(ll,ielem) /= -1) then
            edgeptr=edge%getmap(ll,ielem)+1
            do k=1,vlyr
                iptr = (kptr+k-1)+edgeptr
                v(np ,np,k)=v(np,np,k)+edge%receive(iptr)
            enddo
        endif
    end do

! NWEST
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        if(edge%getmap(ll,ielem) /= -1) then
            edgeptr=edge%getmap(ll,ielem)+1
            do k=1,vlyr
                iptr = (kptr+k-1)+edgeptr
                v(1  ,np,k)=v(1 ,np,k)+edge%receive(iptr)
            enddo
        endif
    end do


  end subroutine edgeVunpack
!
  subroutine edgeVunpackVert(edge,v,ielem)
    use control_mod,            only: north, south, east, west, neast, nwest, seast, swest
    use dimensions_mod,         only: np, max_corner_elem, ne
    use coordinate_systems_mod, only: cartesian3D_t

    type (EdgeBuffer_t),  intent(inout) :: edge
    type (cartesian3D_t), intent(inout)   :: v(:,:,:)
    integer,              intent(in)    :: ielem

    ! Local
    logical, parameter :: UseUnroll = .TRUE.
    integer            :: i,k,l, nce
    integer            :: is,ie,in,iw,ine,inw,isw,ise

    threadsafe=.false.

    if (max_corner_elem.ne.1 .and. ne==0) then
        ! MNL: this is used to construct the dual grid on the cube,
        !      currently only supported for the uniform grid. If
        !      this is desired on a refined grid, a little bit of
        !      work will be required.
        call endrun("edgeVunpackVert should not be called with unstructured meshes")
    end if

    is=edge%getmap(south,ielem)
    ie=edge%getmap(east,ielem)
    in=edge%getmap(north,ielem)
    iw=edge%getmap(west,ielem)


    ! N+S
    do i=1,np/2
       ! North
       v(3,i ,np)%x = edge%receive(in+i)
       v(3,i ,np)%y = edge%receive(np+in+i)
       v(3,i ,np)%z = edge%receive(2*np+in+i)

       ! South
       v(2,i ,1)%x  = edge%receive(is+i)
       v(2,i ,1)%y  = edge%receive(np+is+i)
       v(2,i ,1)%z  = edge%receive(2*np+is+i)
    enddo

    do i=np/2+1,np
       ! North
       v(4,i ,np)%x = edge%receive(in+i)
       v(4,i ,np)%y = edge%receive(np+in+i)
       v(4,i ,np)%z = edge%receive(2*np+in+i)
       ! South
       v(1,i ,1)%x  = edge%receive(is+i)
       v(1,i ,1)%y  = edge%receive(np+is+i)
       v(1,i ,1)%z  = edge%receive(2*np+is+i)
    enddo

    do i=1,np/2
       ! East
       v(3,np,i)%x = edge%receive(ie+i)
       v(3,np,i)%y = edge%receive(np+ie+i)
       v(3,np,i)%z = edge%receive(2*np+ie+i)
       ! West
       v(4,1,i)%x  = edge%receive(iw+i)
       v(4,1,i)%y  = edge%receive(np+iw+i)
       v(4,1,i)%z  = edge%receive(2*np+iw+i)
    end do

    do i=np/2+1,np
       ! East
       v(2,np,i)%x = edge%receive(ie+i)
       v(2,np,i)%y = edge%receive(np+ie+i)
       v(2,np,i)%z = edge%receive(2*np+ie+i)
       ! West
       v(1,1,i)%x  = edge%receive(iw+i)
       v(1,1,i)%y  = edge%receive(np+iw+i)
       v(1,1,i)%z  = edge%receive(2*np+iw+i)
    end do

! SWEST
    nce = max_corner_elem
    do l=swest,swest+max_corner_elem-1
       ! find the one active corner, then exist
        isw=edge%getmap(l,ielem)
        if(isw /= -1) then
            v(1,1,1)%x=edge%receive(isw+1)
            v(1,1,1)%y=edge%receive(nce+isw+1)
            v(1,1,1)%z=edge%receive(2*nce+isw+1)
            exit
        else
            v(1,1,1)%x=0_r8
            v(1,1,1)%y=0_r8
            v(1,1,1)%z=0_r8
        endif
    end do

! SEAST
    do l=swest+max_corner_elem,swest+2*max_corner_elem-1
       ! find the one active corner, then exist
        ise=edge%getmap(l,ielem)
        if(ise /= -1) then
            v(2,np,1)%x=edge%receive(ise+1)
            v(2,np,1)%y=edge%receive(nce+ise+1)
            v(2,np,1)%z=edge%receive(2*nce+ise+1)
            exit
        else
            v(2,np,1)%x=0_r8
            v(2,np,1)%y=0_r8
            v(2,np,1)%z=0_r8
        endif
    end do

! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
       ! find the one active corner, then exist
        ine=edge%getmap(l,ielem)
        if(ine /= -1) then
            v(3,np,np)%x=edge%receive(ine+1)
            v(3,np,np)%y=edge%receive(nce+ine+1)
            v(3,np,np)%z=edge%receive(2*nce+ine+1)
            exit
        else
            v(3,np,np)%x=0_r8
            v(3,np,np)%y=0_r8
            v(3,np,np)%z=0_r8
        endif
    end do

! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
       ! find the one active corner, then exist
        inw = edge%getmap(l,ielem)
        if(inw/= -1) then
            v(4,1,np)%x=edge%receive(inw+1)
            v(4,1,np)%y=edge%receive(nce+inw+1)
            v(4,1,np)%z=edge%receive(2*nce+inw+1)
            exit
        else
            v(4,1,np)%x=0_r8
            v(4,1,np)%y=0_r8
            v(4,1,np)%z=0_r8
        endif
    end do

    ! Fill the missing vertex info

    do i=2,np/2
       ! North
       v(4,i ,np)%x = v(3,i-1 ,np)%x
       v(4,i ,np)%y = v(3,i-1 ,np)%y
       v(4,i ,np)%z = v(3,i-1 ,np)%z
       ! South
       v(1,i ,1)%x  = v(2,i-1 ,1)%x
       v(1,i ,1)%y  = v(2,i-1 ,1)%y
       v(1,i ,1)%z  = v(2,i-1 ,1)%z
    enddo

    do i=np/2+1,np-1
       ! North
       v(3,i ,np)%x = v(4,i+1 ,np)%x
       v(3,i ,np)%y = v(4,i+1 ,np)%y
       v(3,i ,np)%z = v(4,i+1 ,np)%z
       ! South
       v(2,i ,1)%x  = v(1,i+1 ,1)%x
       v(2,i ,1)%y  = v(1,i+1 ,1)%y
       v(2,i ,1)%z  = v(1,i+1 ,1)%z
    enddo

    do i=2,np/2
       ! East
       v(2,np,i)%x = v(3,np,i-1)%x
       v(2,np,i)%y = v(3,np,i-1)%y
       v(2,np,i)%z = v(3,np,i-1)%z
       ! West
       v(1,1,i)%x  = v(4,1,i-1)%x
       v(1,1,i)%y  = v(4,1,i-1)%y
       v(1,1,i)%z  = v(4,1,i-1)%z
    end do

    do i=np/2+1,np-1
       ! East
       v(3,np,i)%x = v(2,np,i+1)%x
       v(3,np,i)%y = v(2,np,i+1)%y
       v(3,np,i)%z = v(2,np,i+1)%z
       ! West
       v(4,1,i)%x  = v(1,1,i+1)%x
       v(4,1,i)%y  = v(1,1,i+1)%y
       v(4,1,i)%z  = v(1,1,i+1)%z
    end do

  end subroutine edgeVunpackVert
  ! ========================================
  ! edgeDGVunpack:
  !
  ! Unpack edges from edge buffer into v...
  ! ========================================

  subroutine edgeDGVunpack(edge,v,vlyr,kptr,ielem)
    use dimensions_mod, only: np, max_corner_elem
    use control_mod,    only: north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t), intent(in)    :: edge
    integer,             intent(in)    :: vlyr
    real (kind=r8),      intent(inout) :: v(0:np+1,0:np+1,vlyr)
    integer,             intent(in)    :: kptr
    integer,             intent(in)    :: ielem

    ! Local
    integer :: i,k,iptr
    integer :: is,ie,in,iw

    threadsafe=.false.

    is=edge%getmap(south,ielem)
    ie=edge%getmap(east,ielem)
    in=edge%getmap(north,ielem)
    iw=edge%getmap(west,ielem)
    do k=1,vlyr
       iptr=np*(kptr+k-1)
       do i=1,np
          v(i   ,0   ,k)=edge%receive(iptr+is+i)
          v(np+1,i   ,k)=edge%receive(iptr+ie+i)
          v(i   ,np+1,k)=edge%receive(iptr+in+i)
          v(0   ,i   ,k)=edge%receive(iptr+iw+i)
       end do
    end do

    i = swest
    if(edge%getmap(i,ielem) /= -1) then
      do k=1,vlyr
        iptr=(kptr+k-1)
        v(0,0,k) = edge%receive(iptr+edge%getmap(i,ielem)+1)
      end do
    end if
    i = swest+max_corner_elem
    if(edge%getmap(i,ielem) /= -1) then
      do k=1,vlyr
        iptr=(kptr+k-1)
        v(np+1,0,k) = edge%receive(iptr+edge%getmap(i,ielem)+1)
      end do
    end if
    i = swest+3*max_corner_elem
    if(edge%getmap(i,ielem) /= -1) then
      do k=1,vlyr
        iptr=(kptr+k-1)
        v(np+1,np+1,k) = edge%receive(iptr+edge%getmap(i,ielem)+1)
      end do
    end if
    i = swest+2*max_corner_elem
    if(edge%getmap(i,ielem) /= -1) then
      do k=1,vlyr
        iptr=(kptr+k-1)
        v(0,np+1,k) = edge%receive(iptr+edge%getmap(i,ielem)+1)
      end do
    end if

  end subroutine edgeDGVunpack

  ! ========================================
  ! edgeVunpackMIN/MAX:
  !
  ! Finds the Min/Max edges from edge buffer into v...
  ! ========================================
  subroutine edgeVunpackMAX(edge,v,vlyr,kptr,ielem)
    use dimensions_mod, only: np, max_corner_elem
    use control_mod,    only: north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t), intent(in)    :: edge
    integer,             intent(in)    :: vlyr
    real (kind=r8),      intent(inout) :: v(np,np,vlyr)
    integer,             intent(in)    :: kptr
    integer,             intent(in)    :: ielem

    ! Local
    integer :: i,k,l,iptr
    integer :: is,ie,in,iw

    threadsafe=.false.

    is=edge%getmap(south,ielem)
    ie=edge%getmap(east,ielem)
    in=edge%getmap(north,ielem)
    iw=edge%getmap(west,ielem)
    do k=1,vlyr
       iptr=np*(kptr+k-1)
       do i=1,np
          v(np ,i  ,k) = MAX(v(np ,i  ,k),edge%receive(iptr+ie+i  ))
          v(i  ,1  ,k) = MAX(v(i  ,1  ,k),edge%receive(iptr+is+i  ))
          v(i  ,np ,k) = MAX(v(i  ,np ,k),edge%receive(iptr+in+i  ))
          v(1  ,i  ,k) = MAX(v(1  ,i  ,k),edge%receive(iptr+iw+i  ))
       end do
    end do

! SWEST
    do l=swest,swest+max_corner_elem-1
        if(edge%getmap(l,ielem) /= -1) then
            do k=1,vlyr
                v(1  ,1 ,k)=MAX(v(1 ,1 ,k),edge%receive((kptr+k-1)+edge%getmap(l,ielem)+1))
            enddo
        endif
    end do

! SEAST
    do l=swest+max_corner_elem,swest+2*max_corner_elem-1
        if(edge%getmap(l,ielem) /= -1) then
            do k=1,vlyr
                v(np ,1 ,k)=MAX(v(np,1 ,k),edge%receive((kptr+k-1)+edge%getmap(l,ielem)+1))
            enddo
        endif
    end do

! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        if(edge%getmap(l,ielem) /= -1) then
            do k=1,vlyr
                v(np ,np,k)=MAX(v(np,np,k),edge%receive((kptr+k-1)+edge%getmap(l,ielem)+1))
            enddo
        endif
    end do

! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        if(edge%getmap(l,ielem) /= -1) then
            do k=1,vlyr
                v(1  ,np,k)=MAX(v(1 ,np,k),edge%receive((kptr+k-1)+edge%getmap(l,ielem)+1))
            enddo
        endif
    end do

  end subroutine edgeVunpackMAX

  subroutine edgeSunpackMAX(edge,v,vlyr,kptr,ielem)
    use dimensions_mod, only: np, max_corner_elem
    use control_mod,    only: north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t), intent(in)    :: edge
    integer,             intent(in)    :: vlyr
    real (kind=r8),      intent(inout) :: v(vlyr)
    integer,             intent(in)    :: kptr
    integer,             intent(in)    :: ielem

    ! Local
    integer :: i,k,l,iptr
    integer :: is,ie,in,iw,edgeptr

    threadsafe=.false.

    is=edge%getmap(south,ielem)
    ie=edge%getmap(east,ielem)
    in=edge%getmap(north,ielem)
    iw=edge%getmap(west,ielem)
    do k=1,vlyr
       iptr=(kptr+k-1)
       v(k) = MAX(v(k),edge%receive(iptr+is+1),edge%receive(iptr+ie+1),edge%receive(iptr+in+1),edge%receive(iptr+iw+1))
    end do

! SWEST
    do l=swest,swest+max_corner_elem-1
        if(edge%getmap(l,ielem) /= -1) then
            edgeptr = edge%getmap(l,ielem)+1
            do k=1,vlyr
                iptr = (kptr+k-1)+edgeptr
                v(k)=MAX(v(k),edge%receive(iptr))
            enddo
        endif
    end do

! SEAST
    do l=swest+max_corner_elem,swest+2*max_corner_elem-1
        if(edge%getmap(l,ielem) /= -1) then
            edgeptr = edge%getmap(l,ielem)+1
            do k=1,vlyr
                iptr = (kptr+k-1)+edgeptr
                v(k)=MAX(v(k),edge%receive(iptr))
            enddo
        endif
    end do

! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        if(edge%getmap(l,ielem) /= -1) then
            edgeptr = edge%getmap(l,ielem)+1
            do k=1,vlyr
                iptr = (kptr+k-1)+edgeptr
                v(k)=MAX(v(k),edge%receive(iptr))
            enddo
        endif
    end do

! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        if(edge%getmap(l,ielem) /= -1) then
            edgeptr = edge%getmap(l,ielem)+1
            do k=1,vlyr
                iptr = (kptr+k-1)+edgeptr
                v(k)=MAX(v(k),edge%receive(iptr))
            enddo
        endif
    end do

  end subroutine edgeSunpackMAX

  subroutine edgeSunpackMIN(edge,v,vlyr,kptr,ielem)
    use dimensions_mod, only: np, max_corner_elem
    use control_mod,    only: north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t), intent(in)    :: edge
    integer,             intent(in)    :: vlyr
    real (kind=r8),      intent(inout) :: v(vlyr)
    integer,             intent(in)    :: kptr
    integer,             intent(in)    :: ielem

    ! Local
    integer :: i,k,l,iptr
    integer :: is,ie,in,iw,edgeptr

    threadsafe=.false.

    is=edge%getmap(south,ielem)
    ie=edge%getmap(east,ielem)
    in=edge%getmap(north,ielem)
    iw=edge%getmap(west,ielem)
    do k=1,vlyr
       iptr=(kptr+k-1)
       v(k) = MIN(v(k),edge%receive(iptr+is+1),edge%receive(iptr+ie+1),edge%receive(iptr+in+1),edge%receive(iptr+iw+1))
    end do

! SWEST
    do l=swest,swest+max_corner_elem-1
        if(edge%getmap(l,ielem) /= -1) then
            edgeptr = edge%getmap(l,ielem)+1
            do k=1,vlyr
                iptr = (kptr+k-1)+edgeptr
                v(k)=MiN(v(k),edge%receive(iptr))
            enddo
        endif
    end do

! SEAST
    do l=swest+max_corner_elem,swest+2*max_corner_elem-1
        if(edge%getmap(l,ielem) /= -1) then
            edgeptr = edge%getmap(l,ielem)+1
            do k=1,vlyr
                iptr = (kptr+k-1)+edgeptr
                v(k)=MIN(v(k),edge%receive(iptr))
            enddo
        endif
    end do

! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        if(edge%getmap(l,ielem) /= -1) then
            edgeptr = edge%getmap(l,ielem)+1
            do k=1,vlyr
                iptr = (kptr+k-1)+edgeptr
                v(k)=MIN(v(k),edge%receive(iptr))
            enddo
        endif
    end do

! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        if(edge%getmap(l,ielem) /= -1) then
            edgeptr = edge%getmap(l,ielem)+1
            do k=1,vlyr
                iptr = (kptr+k-1)+edgeptr
                v(k)=MIN(v(k),edge%receive(iptr))
            enddo
        endif
    end do

  end subroutine edgeSunpackMIN

  subroutine edgeVunpackMIN(edge,v,vlyr,kptr,ielem)
    use dimensions_mod, only: np, max_corner_elem
    use control_mod,    only: north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t), intent(in)    :: edge
    integer,             intent(in)    :: vlyr
    real (kind=r8),      intent(inout) :: v(np,np,vlyr)
    integer,             intent(in)    :: kptr
    integer,             intent(in)    :: ielem

    ! Local
    integer :: i,k,l,iptr
    integer :: is,ie,in,iw,edgeptr

    threadsafe=.false.

    is=edge%getmap(south,ielem)
    ie=edge%getmap(east,ielem)
    in=edge%getmap(north,ielem)
    iw=edge%getmap(west,ielem)
    do k=1,vlyr
       iptr = np*(kptr+k-1)
       do i=1,np
          v(np ,i  ,k) = MIN(v(np ,i  ,k),edge%receive(iptr+ie+i  ))
          v(i  ,1  ,k) = MIN(v(i  ,1  ,k),edge%receive(iptr+is+i  ))
          v(i  ,np ,k) = MIN(v(i  ,np ,k),edge%receive(iptr+in+i  ))
          v(1  ,i  ,k) = MIN(v(1  ,i  ,k),edge%receive(iptr+iw+i  ))
       end do
    end do

! SWEST
    do l=swest,swest+max_corner_elem-1
        if(edge%getmap(l,ielem) /= -1) then
            edgeptr=edge%getmap(l,ielem)+1
            do k=1,vlyr
                iptr=(kptr+k-1)+edgeptr
                v(1  ,1 ,k)=MIN(v(1 ,1 ,k),edge%receive(iptr))
            enddo
        endif
    end do

! SEAST
    do l=swest+max_corner_elem,swest+2*max_corner_elem-1
        if(edge%getmap(l,ielem) /= -1) then
            edgeptr=edge%getmap(l,ielem)+1
            do k=1,vlyr
                iptr=(kptr+k-1)+edgeptr
                v(np ,1 ,k)=MIN(v(np,1 ,k),edge%receive(iptr))
            enddo
        endif
    end do

! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        if(edge%getmap(l,ielem) /= -1) then
            edgeptr=edge%getmap(l,ielem)+1
            do k=1,vlyr
                iptr=(kptr+k-1)+edgeptr
                v(np ,np,k)=MIN(v(np,np,k),edge%receive(iptr))
            enddo
        endif
    end do

! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        if(edge%getmap(l,ielem) /= -1) then
            edgeptr=edge%getmap(l,ielem)+1
            do k=1,vlyr
                iptr=(kptr+k-1)+edgeptr
                v(1  ,np,k)=MIN(v(1 ,np,k),edge%receive(iptr))
            enddo
        endif
    end do

  end subroutine edgeVunpackMIN

  ! ========================================
  ! LongEdgeVunpackMIN:
  !
  ! Finds the Min edges from edge buffer into v...
  ! ========================================
  subroutine LongEdgeVunpackMIN(edge,v,vlyr,kptr,desc)
    use control_mod,    only: north, south, east, west, neast, nwest, seast, swest
    use dimensions_mod, only: np, max_corner_elem

    type (LongEdgeBuffer_t), intent(in)    :: edge
    integer,                 intent(in)    :: vlyr
    integer ,                intent(inout) :: v(np,np,vlyr)
    integer,                 intent(in)    :: kptr
    type (EdgeDescriptor_t), intent(in)    :: desc

    ! Local

    integer :: i,k,l
    integer :: is,ie,in,iw

    threadsafe=.false.

    is=desc%getmapP(south)
    ie=desc%getmapP(east)
    in=desc%getmapP(north)
    iw=desc%getmapP(west)
    do k=1,vlyr
       do i=1,np
          v(i  ,1  ,k) = MIN(v(i  ,1  ,k),edge%buf(kptr+k,is+i  ))
          v(np ,i  ,k) = MIN(v(np ,i  ,k),edge%buf(kptr+k,ie+i  ))
          v(i  ,np ,k) = MIN(v(i  ,np ,k),edge%buf(kptr+k,in+i  ))
          v(1  ,i  ,k) = MIN(v(1  ,i  ,k),edge%buf(kptr+k,iw+i  ))
       end do
    end do

! SWEST
    do l=swest,swest+max_corner_elem-1
        if(desc%getmapP(l) /= -1) then
            do k=1,vlyr
                v(1  ,1 ,k)=MIN(v(1 ,1 ,k),edge%buf(kptr+k,desc%getmapP(l)+1))
            enddo
        endif
    end do

! SEAST
    do l=swest+max_corner_elem,swest+2*max_corner_elem-1
        if(desc%getmapP(l) /= -1) then
            do k=1,vlyr
                v(np ,1 ,k)=MIN(v(np,1 ,k),edge%buf(kptr+k,desc%getmapP(l)+1))
            enddo
        endif
    end do

! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        if(desc%getmapP(l) /= -1) then
            do k=1,vlyr
                v(np ,np,k)=MIN(v(np,np,k),edge%buf(kptr+k,desc%getmapP(l)+1))
            enddo
        endif
    end do

! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        if(desc%getmapP(l) /= -1) then
            do k=1,vlyr
                v(1  ,np,k)=MIN(v(1 ,np,k),edge%buf(kptr+k,desc%getmapP(l)+1))
            enddo
        endif
    end do

  end subroutine LongEdgeVunpackMIN


subroutine ghostpack(edge,v,vlyr,kptr,ielem)
  
  use dimensions_mod, only : max_corner_elem
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  use edgetype_mod, only : EdgeDescriptor_t

  implicit none

  type (Edgebuffer_t)                :: edge
  integer,              intent(in)   :: vlyr
  integer,              intent(in)   :: kptr

  real (kind=r8),intent(in)   :: v(edge%lb:edge%ub,edge%lb:edge%ub,vlyr)
  integer,              intent(in)   :: ielem

  ! Local variables
  integer :: i,j,k,ir,l,itr,ktmp

  integer :: is,ie,in,iw,isw,ise,inw,ine
  integer :: nhc, npoints
  integer :: edgeptr,iptr

  is = edge%putmap(south,ielem)
  ie = edge%putmap(east,ielem)
  in = edge%putmap(north,ielem)
  iw = edge%putmap(west,ielem)
  if (edge%nlyr < (kptr+vlyr) ) then
       print *,'edge%nlyr = ',edge%nlyr
       print *,'kptr+vlyr = ',kptr+vlyr
       call endrun('ghostpack: Buffer overflow: size of the vertical dimension must be increased!')
  endif


  nhc     = edge%ndepth
  npoints = edge%npoints

    !DIR$ IVDEP
    do k=1,vlyr
      ktmp = nhc*(kptr+k-1)
      do j=1,nhc
        iptr = npoints*(ktmp + j - 1)
        do i=1,npoints
          edge%buf(iptr+is+i)   = v(i  ,j ,k)
          edge%buf(iptr+ie+i)   = v(npoints-j+1 ,i ,k)
          edge%buf(iptr+in+i)   = v(i  ,npoints-j+1,k)
          edge%buf(iptr+iw+i)   = v(j  ,i ,k)
        enddo
      end do
    end do


  !  This is really kludgy way to setup the index reversals
  !  But since it is so a rare event not real need to spend time optimizing
  !  Check if the edge orientation of the recieving element is different
  !  if it is, swap the order of data in the edge
  if(edge%reverse(south,ielem)) then
     !DIR$ IVDEP
     do k=1,vlyr
       ktmp = nhc*(kptr+k-1)
       do j=1,nhc
         iptr = npoints*(ktmp + j - 1)
         do i=1,npoints
           ir = npoints-i+1
           edge%buf(iptr+is+i)=v(ir,j,k)
         enddo
       enddo
     enddo
  endif

  if(edge%reverse(east,ielem)) then
     !DIR$ IVDEP
     do k=1,vlyr
       ktmp = nhc*(kptr+k-1)
       do j=1,nhc
         iptr = npoints*(ktmp + j - 1)
         do i=1,npoints
           ir = npoints-i+1
           edge%buf(iptr+ie+i)=v(npoints-j+1,ir,k)
          enddo
        enddo
      enddo
  endif

  if(edge%reverse(north,ielem)) then
     !DIR$ IVDEP
      do k=1,vlyr
        ktmp = nhc*(kptr+k-1)
        do j=1,nhc
         iptr = npoints*(ktmp + j - 1)
          do i=1,npoints
            ir = npoints-i+1
            edge%buf(iptr+in+i)=v(ir,npoints-j+1,k)
          enddo
        enddo
      enddo
  endif

  if(edge%reverse(west,ielem)) then
     !DIR$ IVDEP
     do k=1,vlyr
       ktmp = nhc*(kptr+k-1)
       do j=1,nhc
         iptr = npoints*(ktmp + j - 1)
         do i=1,npoints
            ir = npoints-i+1
            edge%buf(iptr+iw+i)=v(j,ir,k)
          enddo
        enddo
      enddo
  endif


  ! corners.  this is difficult because we dont know the orientaton
  ! of the corners, and this which (i,j) dimension maps to which dimension
! SWEST
  do l=swest,swest+max_corner_elem-1
     if (edge%putmap(l,ielem) /= -1) then
         isw = edge%putmap(l,ielem)
         !DIR$ IVDEP
         do k=1,vlyr
          ktmp = nhc*(kptr+k-1)
          do j=1,nhc
            iptr = nhc*(ktmp + j - 1)
            do i=1,nhc
              edge%buf(iptr+isw+i)=v(i  ,j ,k)
            enddo
          end do
        end do
     end if
  end do

! SEAST
  do l=swest+max_corner_elem,swest+2*max_corner_elem-1
     if (edge%putmap(l,ielem) /= -1) then
         ise = edge%putmap(l,ielem)
         !DIR$ IVDEP
         do k=1,vlyr
          ktmp = nhc*(kptr+k-1)
          do j=1,nhc
            iptr = nhc*(ktmp + j - 1)
            do i=1,nhc
              edge%buf(iptr+ise+i)=v(npoints-i+1 ,j ,k)
            enddo
          end do
        end do
     end if
  end do

! NEAST
  do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
     if (edge%putmap(l,ielem) /= -1) then
        ine = edge%putmap(l,ielem)
         !DIR$ IVDEP
        do k=1,vlyr
           ktmp = nhc*(kptr+k-1)
           do j=1,nhc
              iptr = nhc*(ktmp + j - 1)
              do i=1,nhc
                 edge%buf(iptr+ine+i)=v(npoints-i+1,npoints-j+1,k)
              enddo
            enddo
          end do
      end if
  end do

! NWEST
  do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
     if (edge%putmap(l,ielem) /= -1) then
        inw = edge%putmap(l,ielem)
        !DIR$ IVDEP
        do k=1,vlyr
          ktmp = nhc*(kptr+k-1)
          do j=1,nhc
            iptr = nhc*(ktmp + j - 1)
            do i=1,nhc
              edge%buf(iptr+inw+i)=v(i  ,npoints-j+1,k)
            enddo
          end do
        end do
     end if
  end do

end subroutine ghostpack

subroutine ghostunpack(edge,v,vlyr,kptr,ielem)
  use dimensions_mod, only : max_corner_elem
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  type (Edgebuffer_t),         intent(in)  :: edge

  integer,               intent(in)  :: vlyr
  integer,               intent(in)  :: kptr
  integer,               intent(in)  :: ielem

  real (kind=r8), intent(inout) :: v(edge%lb:edge%ub,edge%lb:edge%ub,vlyr)


  ! Local
  logical, parameter :: UseUnroll = .TRUE.
  integer :: i,j,k,l,itr, ktmp
  integer :: is,ie,in,iw,isw,ise,inw,ine
  integer :: nhc,npoints,iptr
  logical :: reverse

  threadsafe=.false.

  is=edge%getmap(south,ielem)
  ie=edge%getmap(east,ielem)
  in=edge%getmap(north,ielem)
  iw=edge%getmap(west,ielem)

  nhc     = edge%ndepth
  npoints = edge%npoints

  ! example for north buffer
  ! first row ('edge') goes in v(:,np+1,k)
  ! 2nd   row ('edge') goes in v(:,np+2,k)
  ! etc...
    !DIR$ IVDEP
    do k=1,vlyr
      ktmp = nhc*(kptr+k-1)
      do j=1,nhc
        iptr = npoints*(ktmp + j - 1)
        do i=1,npoints
          v(i  ,1-j  ,k)      = edge%receive(iptr+is+i)  ! South
          v(npoints+j ,i  ,k) = edge%receive(iptr+ie+i)  ! East
          v(i  ,npoints+j ,k) = edge%receive(iptr+in+i)  ! North
          v(1-j  ,i  ,k)      = edge%receive(iptr+iw+i)  ! West
        end do
      end do
    end do


! SWEST
  do l=swest,swest+max_corner_elem-1
     isw = edge%getmap(l,ielem)
     if(isw /= -1) then
        ! note the following is the the correct meaning of reverse in this code.  
        ! It is  best described as a transponse operation 
        if (edge%reverse(l,ielem)) then
           do k=1,vlyr
              ktmp = nhc*(kptr+k-1)
              do j=1,nhc
                 iptr = nhc*(ktmp + j - 1)
                 do i=1,nhc
                    v(1-j,1-i,k)=edge%receive(iptr+isw+i)
                 enddo
              enddo
             enddo
        else
           do k=1,vlyr
              ktmp = nhc*(kptr+k-1)
              do i=1,nhc
                 iptr = nhc*(ktmp + i - 1)
                 do j=1,nhc
                    v(1-j,1-i,k)=edge%receive(iptr+isw+j)
                 enddo
              enddo
             enddo
        endif
     else
         do k=1,vlyr
           do j=1,nhc
             do i=1,nhc
               v(1-i,1-j,k)=edgeDefaultVal
             enddo
           enddo
         enddo
     endif
  end do

! SEAST
  do l=swest+max_corner_elem,swest+2*max_corner_elem-1
     ise = edge%getmap(l,ielem)
     if(ise /= -1) then 
        if (edge%reverse(l,ielem)) then
           do k=1,vlyr
              ktmp = nhc*(kptr+k-1)
              do i=1,nhc
                 iptr = nhc*(ktmp + i - 1)
                 do j=1,nhc
                    v(npoints+i,1-j,k)=edge%receive(iptr+ise+j)
                 enddo
              enddo
            enddo
        else
           do k=1,vlyr
              ktmp = nhc*(kptr+k-1)
              do j=1,nhc
                 iptr = nhc*(ktmp + j - 1)
                 do i=1,nhc
                    v(npoints+i ,1-j ,k)=edge%receive(iptr+ise+i)
                 enddo
              enddo
            enddo
        endif
      else
         do k=1,vlyr
          do j=1,nhc
            do i=1,nhc
              v(npoints+i,1-j,k)=edgeDefaultVal
            enddo
          enddo
         enddo
     endif
  end do

! NEAST
  do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
     ine = edge%getmap(l,ielem)
     if(ine /= -1) then 
        if (edge%reverse(l,ielem)) then
           do k=1,vlyr
              ktmp = nhc*(kptr+k-1)
              do j=1,nhc
                 do i=1,nhc
                    iptr = nhc*(ktmp + i - 1)
                    v(npoints+i ,npoints+j,k)=edge%receive(iptr+ine+j)
                 enddo
              enddo
            enddo
        else
           do k=1,vlyr
              ktmp = nhc*(kptr+k-1)
              do j=1,nhc
                 iptr = nhc*(ktmp + j - 1)
                 do i=1,nhc
                    v(npoints+i ,npoints+j,k)=edge%receive(iptr+ine+i)
                 enddo
              enddo
            enddo
        endif
      else
         do k=1,vlyr
          do j=1,nhc
            do i=1,nhc
              v(npoints+i,npoints+j,k)=edgeDefaultVal
            enddo
          enddo
         enddo
     endif
  end do

! NWEST
  do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
     inw = edge%getmap(l,ielem)
     if(inw /= -1) then 
        if (edge%reverse(l,ielem)) then
           do k=1,vlyr
              ktmp = nhc*(kptr+k-1)
              do i=1,nhc
                 iptr = nhc*(ktmp + i - 1)
                 do j=1,nhc
                    v(1-i ,npoints+j,k)=edge%receive(iptr+inw+j)
                 enddo
              enddo
            enddo
        else
           do k=1,vlyr
              ktmp = nhc*(kptr+k-1)
              do j=1,nhc
                 iptr = nhc*(ktmp + j - 1)
                 do i=1,nhc
                    v(1-i ,npoints+j,k)=edge%receive(iptr+inw+i)
                 enddo
              enddo
            enddo
        endif
      else
         do k=1,vlyr
          do j=1,nhc
            do i=1,nhc
              v(1-i,npoints+j,k)=edgeDefaultVal
            enddo
          enddo
         enddo
     endif
  end do

end subroutine ghostunpack

  ! =========================================
  ! initGhostBuffer3d:
  ! Author: James Overfelt
  ! create an Real based communication buffer
  ! npoints is the number of points on one side
  ! nhc is the deep of the ghost/halo zone
  ! =========================================
  subroutine initGhostBuffer3d(ghost,nlyr,np,nhc_in)

    implicit none
    integer,intent(in)                  :: nlyr, np
    integer,intent(in),optional         :: nhc_in
    type (Ghostbuffer3d_t),intent(out)  :: ghost

    ! Local variables

    integer :: nbuf,nhc,i

    ! sanity check for threading
    if (omp_get_num_threads()>1) then
       call endrun('ERROR: initGhostBuffer must be called before threaded region')
    endif

    if (present(nhc_in)) then
       nhc=nhc_in
    else
       nhc = np-1
    endif

    nbuf=max_neigh_edges*nelemd

    ghost%nlyr    = nlyr
    ghost%nhc     = nhc
    ghost%np      = np
    ghost%nbuf    = nbuf
    ghost%elem_size = np*(nhc+1)
    allocate(ghost%buf    (np,(nhc+1),nlyr,nbuf))
    allocate(ghost%receive(np,(nhc+1),nlyr,nbuf))
    ghost%buf=0
    ghost%receive=0

  end subroutine initGhostBuffer3d

  ! =================================================================================
  ! GHOSTVPACK3D
  ! AUTHOR: James Overfelt (from a subroutine of Christoph Erath, ghostvpack2D)
  ! Pack edges of v into an ghost buffer for boundary exchange.
  !
  ! This subroutine packs for many vertical layers into an ghost
  ! buffer.
  ! If the buffer associated with edge is not large enough to
  ! hold all vertical layers you intent to pack, the method will
  ! halt the program with a call to endrun().
  ! INPUT:
  ! - ghost Buffer into which the data will be packed.
  !   This buffer must be previously allocated with initGhostBuffer().
  ! - v The data to be packed.
  ! - nhc deep of ghost/halo zone
  ! - npoints number of points on on side
  ! - kptr Vertical pointer to the place in the edge buffer where
  ! data will be located.
  ! =================================================================================
  subroutine ghostVpack3d(ghost, v, vlyr, kptr, desc)
    use dimensions_mod, only : max_corner_elem
    use control_mod,    only : north, south, east, west, neast, nwest, seast, swest
    use edgetype_mod, only : edgedescriptor_t, ghostbuffer3d_t
    implicit none

    type (Ghostbuffer3d_t)                :: ghost
    integer,              intent(in)      :: kptr,vlyr
    real (kind=r8),intent(in)      :: v(ghost%np, ghost%np, vlyr)
    type (EdgeDescriptor_t),intent(in)    :: desc

    integer                               :: nhc, np

    ! Local variables
    integer :: i,j,k,ir,l,e

    integer :: is,ie,in,iw

    if(.not. threadsafe) then
!$OMP BARRIER
       threadsafe=.true.
    end if
    ! Example convenction for buffer to the north:
    !  buf(:,,:,i,e)
    !   each "edge" is a row of data (i=1,np) in the element
    !     north most row of data goes into e=1
    !     next row of data goes into e=2
    !     ....
    !     south most row of data goes into e=np
    !   We need to pack this way to preserve the orientation
    !   so the data can be unpacked correctly

    ! note: we think of buf as dimensioned buf(k,is,i,e)
    ! but this array is flatted to:   buf(k,is+(i-1)+(e-1)*np)
    !
    nhc = ghost%nhc
    np = ghost%np
    is = desc%putmapP_ghost(south)
    ie = desc%putmapP_ghost(east)
    in = desc%putmapP_ghost(north)
    iw = desc%putmapP_ghost(west)

    do k=1,vlyr
       do j=1,nhc
          do i=1,np
             ghost%buf(i,j,kptr+k,is)   = v(i,    j+1 , k)
             ghost%buf(i,j,kptr+k,ie)   = v(np-j,   i , k)
             ghost%buf(i,j,kptr+k,in)   = v(i,   np-j , k)
             ghost%buf(i,j,kptr+k,iw)   = v(j+1,    i , k)
          enddo
       end do
    end do
    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing
    !  Check if the edge orientation of the recieving element is different
    !  if it is, swap the order of data in the edge
    if(desc%reverse(south)) then
       do k=1,vlyr
          do j=1,nhc
             do i=1,np
                ir = np-i+1
                ghost%buf(ir, j, kptr+k, is)=v(i, j+1, k)
             enddo
          enddo
       enddo
    endif

    if(desc%reverse(east)) then
       do k=1,vlyr
          do j=1,nhc
             do i=1,np
                ir = np-i+1
                ghost%buf(ir, j, kptr+k, ie)=v(np-j, i, k)
             enddo
          enddo
       enddo
    endif

    if(desc%reverse(north)) then
       do k=1,vlyr
          do j=1,nhc
              do i=1,np
                ir = np-i+1
                ghost%buf(ir, j, kptr+k, in)=v(i, np-j, k)
             enddo
          enddo
       enddo
    endif

    if(desc%reverse(west)) then
       do k=1,vlyr
          do j=1,nhc
             do i=1,np
                ir = np-i+1
                ghost%buf(ir, j, kptr+k, iw)=v(j+1, i, k)
             enddo
          enddo
       enddo
    endif

    ! corners.  this is difficult because we dont know the orientaton
    ! of the corners, and this which (i,j) dimension maps to which dimension
! SWEST
    do l=swest, swest+max_corner_elem-1
       if (desc%putmapP_ghost(l) /= -1) then
          do k=1,vlyr
             do j=1,nhc+1
                do i=1,nhc+1
                   ghost%buf(i, j, kptr+k, desc%putmapP_ghost(l))=v(i, j, k)
                enddo
             enddo
          enddo
       end if
    end do

! SEAST
    do l=swest+max_corner_elem,swest+2*max_corner_elem-1
       if (desc%putmapP_ghost(l) /= -1) then
          do k=1,vlyr
             do j=1,nhc+1
                do i=1,nhc+1
                   ghost%buf(i, j, kptr+k, desc%putmapP_ghost(l))=v(np-i+1, j, k)
                enddo
             enddo
          enddo
       end if
    end do

! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
       if (desc%putmapP_ghost(l) /= -1) then
          do k=1,vlyr
             do j=1,nhc+1
                do i=1,nhc+1
                   ghost%buf(i, j, kptr+k,desc%putmapP_ghost(l))=v(np-i+1, np-j+1, k)
                enddo
             enddo
          enddo
       end if
    end do

! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
       if (desc%putmapP_ghost(l) /= -1) then
          do k=1,vlyr
             do j=1,nhc+1
                do i=1,nhc+1
                   ghost%buf(i, j, kptr+k,desc%putmapP_ghost(l))=v(i, np-j+1, k)
                enddo
             enddo
          enddo
       end if
    end do
  end subroutine ghostVpack3d

  ! =================================================================================
  ! GHOSTVUNPACK3D
  ! AUTHOR: James Overfelt (from a subroutine of Christoph Erath,
  ! ghostVunpack2d)
  ! Unpack ghost points from ghost buffer into v...
  ! It is for cartesian points (v is only two dimensional).
  ! INPUT SAME arguments as for GHOSTVPACK
  ! =================================================================================

  subroutine ghostVunpack3d(g, v, vlyr, kptr, desc, sw, se, nw, ne, mult)
    use dimensions_mod, only : max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use edgetype_mod, only : edgedescriptor_t, ghostbuffer3d_t
    implicit none
    type (Ghostbuffer3d_t),         intent(in)  :: g

    integer,               intent(in)     :: kptr,vlyr
    real (kind=r8), intent(inout)  :: v (1-g%nhc : g%np+g%nhc, 1-g%nhc : g%np+g%nhc, vlyr)
    integer,               intent(out)    :: mult(5:8)
    real (kind=r8), intent(out)    :: sw(1-g%nhc : 1,          1-g%nhc : 1,          vlyr, max_corner_elem-1)
    real (kind=r8), intent(out)    :: se(   g%np : g%np+g%nhc, 1-g%nhc : 1,          vlyr, max_corner_elem-1)
    real (kind=r8), intent(out)    :: ne(   g%np : g%np+g%nhc,    g%np : g%np+g%nhc, vlyr, max_corner_elem-1)
    real (kind=r8), intent(out)    :: nw(1-g%nhc : 1,             g%np : g%np+g%nhc, vlyr, max_corner_elem-1)
    type (EdgeDescriptor_t)               :: desc

    integer                               :: nhc, np

    ! Local
    logical, parameter :: UseUnroll = .TRUE.
    integer :: i,j,k,l
    integer :: is,ie,in,iw,ic
    logical :: reverse

    threadsafe=.false.

    nhc     = g%nhc
    np      = g%np

    is=desc%getmapP_ghost(south)
    ie=desc%getmapP_ghost(east)
    in=desc%getmapP_ghost(north)
    iw=desc%getmapP_ghost(west)

! fill in optional values with edgeDefaultVal
    do k=1,vlyr
       do j=1,nhc
          do i=1,nhc
             v(1-i,   1-j, k)=edgeDefaultVal
             v(np+i , 1-j, k)=edgeDefaultVal
             v(np+i, np+j, k)=edgeDefaultVal
             v(1-i , np+j, k)=edgeDefaultVal
          enddo
       enddo
    enddo

    ! example for north buffer
    ! first row ('edge') goes in v(:,np+1)
    ! 2nd   row ('edge') goes in v(:,np+2)
    ! etc...

    do k=1,vlyr
       do j=1,nhc
          do i=1,np
             v(i  ,  1-j , k) = g%buf(i,j,kptr+k,is  )
             v(np+j ,  i , k) = g%buf(i,j,kptr+k,ie  )
             v(i  , np+j,  k) = g%buf(i,j,kptr+k,in  )
             v(1-j  ,  i , k) = g%buf(i,j,kptr+k,iw  )
          end do
       end do
    end do

    ! four sides are always just one
    mult(swest) = 0
    mult(seast) = 0
    mult(neast) = 0
    mult(nwest) = 0



! SWEST
    do l=swest, swest+max_corner_elem-1
       ic = desc%getmapP_ghost(l)
       if(ic /= -1) then
          reverse=desc%reverse(l)
          if (mult(swest) .eq. 0) then
            if (reverse) then
               do k=1,vlyr
                  do j=1,nhc
                     do i=1,nhc
                        v(1-i, 1-j, k)=g%buf(j+1, i+1, kptr+k, ic)
                     enddo
                  enddo
               enddo
             else
               do k=1,vlyr
                  do j=1,nhc
                     do i=1,nhc
                        v(1-i,1-j,k)=g%buf(i+1,j+1,kptr+k,ic)
                     enddo
                  enddo
               enddo
            endif
          else
            if (reverse) then
               do k=1,vlyr
                  do j=0,nhc
                     do i=0,nhc
                        sw(1-i,1-j,k,mult(swest))=g%buf(j+1,i+1,kptr+k,ic)
                     enddo
                  enddo
               enddo
             else
                do k=1,vlyr
                   do j=0,nhc
                      do i=0,nhc
                         sw(1-i,1-j,k,mult(swest))=g%buf(i+1,j+1,kptr+k,ic)
                     enddo
                  enddo
               enddo
            endif
          endif
          mult(swest) = mult(swest) + 1
       endif
    end do

! SEAST
    do l=swest+max_corner_elem,swest+2*max_corner_elem-1
       ic = desc%getmapP_ghost(l)
       if(ic /= -1) then
          reverse=desc%reverse(l)
          if (mult(seast) .eq. 0) then
            if (reverse) then
               do k=1,vlyr
                  do j=1,nhc
                     do i=1,nhc
                        v(np+i,1-j,k)=g%buf(j+1,i+1,kptr+k,ic)
                     enddo
                  enddo
               enddo
            else
               do k=1,vlyr
                  do j=1,nhc
                     do i=1,nhc
                        v(np+i ,1-j,k)=g%buf(i+1,j+1,kptr+k,ic)
                     enddo
                  enddo
               enddo
            endif
          else
            if (reverse) then
               do k=1,vlyr
                  do j=0,nhc
                     do i=0,nhc
                        se(np+i,1-j,k,mult(seast))=g%buf(j+1,i+1,kptr+k,ic)
                     enddo
                  enddo
               enddo
            else
               do k=1,vlyr
                  do j=0,nhc
                     do i=0,nhc
                        se(np+i ,1-j,k,mult(seast))=g%buf(i+1,j+1,kptr+k,ic)
                     enddo
                  enddo
               enddo
            endif
          endif
          mult(seast) = mult(seast) + 1
       endif
    end do


! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
       ic = desc%getmapP_ghost(l)
       if(ic /= -1) then
          reverse=desc%reverse(l)
          if (mult(neast) .eq. 0) then
            if (reverse) then
               do k=1,vlyr
                  do j=1,nhc
                     do i=1,nhc
                        v(np+i ,np+j,k)=g%buf(j+1,i+1,kptr+k,ic)
                     enddo
                  enddo
               enddo
            else
               do k=1,vlyr
                  do j=1,nhc
                     do i=1,nhc
                        v(np+i ,np+j,k)=g%buf(i+1,j+1,kptr+k,ic)
                     enddo
                  enddo
               enddo
            endif
          else
            if (reverse) then
               do k=1,vlyr
                  do j=0,nhc
                     do i=0,nhc
                        ne(np+i ,np+j,k,mult(neast))=g%buf(j+1,i+1,kptr+k,ic)
                     enddo
                  enddo
               enddo
            else
               do k=1,vlyr
                  do j=0,nhc
                     do i=0,nhc
                        ne(np+i ,np+j,k,mult(neast))=g%buf(i+1,j+1,kptr+k,ic)
                     enddo
                  enddo
               enddo
            endif
          endif
          mult(neast) = mult(neast) + 1
       endif
    end do

! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
       ic = desc%getmapP_ghost(l)
       if(ic /= -1) then
          reverse=desc%reverse(l)
          if (mult(nwest) .eq. 0) then
             if (reverse) then
                do k=1,vlyr
                   do j=1,nhc
                      do i=1,nhc
                         v(1-i ,np+j,k)=g%buf(j+1,i+1,kptr+k,ic)
                      enddo
                   enddo
                enddo
             else
                do k=1,vlyr
                   do j=1,nhc
                      do i=1,nhc
                         v(1-i ,np+j,k)=g%buf(i+1,j+1,kptr+k,ic)
                      enddo
                   enddo
                enddo
             endif
          else
             if (reverse) then
                do k=1,vlyr
                   do j=0,nhc
                      do i=0,nhc
                         nw(1-i ,np+j,k,mult(nwest))=g%buf(j+1,i+1,kptr+k,ic)
                      enddo
                   enddo
                enddo
             else
                do k=1,vlyr
                   do j=0,nhc
                      do i=0,nhc
                         nw(1-i ,np+j,k,mult(nwest))=g%buf(i+1,j+1,kptr+k,ic)
                      enddo
                   enddo
                enddo
             endif
          endif
          mult(nwest) = mult(nwest) + 1
       endif
    end do

  end subroutine ghostVunpack3d

  subroutine FreeGhostBuffer3D(buffer)
    use edgetype_mod, only : ghostbuffer3d_t
    implicit none
    type (Ghostbuffer3d_t),intent(inout) :: buffer

!$OMP BARRIER
!$OMP MASTER
    buffer%nbuf=0
    buffer%nlyr=0
    deallocate(buffer%buf)
    deallocate(buffer%receive)
!$OMP END MASTER

  end subroutine FreeGhostBuffer3D


End module edge_mod
