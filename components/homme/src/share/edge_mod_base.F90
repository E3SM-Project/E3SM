#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module edge_mod_base

  use kinds, only : int_kind, log_kind, real_kind
  use dimensions_mod, only : max_neigh_edges, nelemd
  use perf_mod, only: t_startf, t_stopf, t_adj_detailf ! _EXTERNAL
  use thread_mod, only: nthreadshoriz, omp_get_num_threads, omp_get_thread_num
  use coordinate_systems_mod, only : cartesian3D_t
  use schedtype_mod, only : cycle_t, schedule_t, schedule
  use parallel_mod, only : abortmp, haltmp, MPIreal_t, iam,parallel_t, &
      MAX_ACTIVE_MSG, HME_status_size, BNDRY_TAG_BASE
  use edgetype_mod, only : edgedescriptor_t, edgebuffer_t, &
      Longedgebuffer_t, Ghostbuffertr_t, Ghostbuffer3d_t, initedgebuffer_callid
  use element_mod, only : element_t


  implicit none
  private
  save

  ! 8-byte Integer routines 
  public :: initLongEdgeBuffer, FreeLongEdgeBuffer
  public :: LongEdgeVpack, LongEdgeVunpackMIN


  public :: initEdgeBuffer, initEdgeSBuffer, FreeEdgeBuffer

  !--------------------------------------------------------- 
  ! Pack/unpack routines that use the New format Edge buffer
  !--------------------------------------------------------- 
  public :: edgeVpack, edgeVunpack
  public :: edgeVunpackMIN, edgeVunpackMAX
  public :: edgeDGVpack, edgeDGVunpack
  public :: edgeVunpackVert
  
  !----------------------------------------------------------------
  ! Pack/unpack routines that communicate a fixed number values 
  ! per element.  This is used to communicate MIN/MAX values from
  ! neighboring elemeents
  !----------------------------------------------------------------
  public :: edgeSpack
  public :: edgeSunpackMIN, edgeSunpackMAX
  public :: edgerotate

  public :: buffermap

  logical, private :: threadsafe=.true.

  real(kind=real_kind), parameter, public :: edgeDefaultVal = 1.11e+100_real_kind

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
!
! Mark Taylor
! pack/unpack full element of data of size (nx,nx)  
! user specifies the size when creating the buffer 
! input/output arrays are cartesian, and will only unpack 1 corner element 
! (even if there are more when running with an unstructured grid) 
! This routine is used mostly for testing and to compute the orientation of
! an elements corner neighbors
!
  public :: initGhostBuffer3D      ! init/free buffers used by pack/unpack full and 3D
  public :: FreeGhostBuffer3D
  public :: ghostVpackfull       
  public :: ghostVunpackfull     
  ! same as above, except orientation of element data is preserved
  ! (so boundary data for two adjacent element may not match up)
  public :: ghostVpack_unoriented   
  public :: ghostVunpack_unoriented     


!
! James Overfelt
! pack/unpack user specifed halo region "nhc".  
! Does not include element edge data (assumes element edge data is C0)
! (appropriate for continuous GLL data where the edge data does not need to be sent)
! support for unstructed meshes via extra output arrays: sw,se,ne,nw
! This routine is currently used by surfaces_mod.F90 to construct the GLL dual grid
!
  public :: ghostVpack3d          ! pack/unpack specifed halo size (up to 1 element)
                                  ! should be identical to ghostVpack2d except for
                                  ! shape of input array
  public :: ghostVunpack3d        ! returns v including populating halo region of v
                                  ! "extra" corner elements are returned in arrays
                                  ! sw,se,ne,nw
! MT TODO: this routine works for unstructed data (where the corner orientation flag is
! not set).  So why dont we remove all the "reverse" checks in unpack?




!
! Christoph Erath
! pack/unpack partial element of data of size (nx,nx) with user specifed halo size nh
! user specifies the sizes when creating the buffer 
! buffer has 1 extra dimension (as compared to subroutines above) for multiple tracers
! input/output arrays are cartesian, and thus assume at most 1 element at each corner
! hence currently only supports cube-sphere grids.
!
! TODO: GhostBufferTR (init and type) should be removed - we only need GhostBuffer3D, 
! if we can fix
! ghostVpack2d below to pass vlyr*ntrac_d instead of two seperate arguments
!
  public :: initGhostBufferTR     ! ghostbufferTR_t
  public :: FreeGhostBufferTR     ! ghostbufferTR_t


! routines which including element edge data  
! (used for FVM arrays where edge data is not shared by neighboring elements)
! these routines pack/unpack element data with user specified halo size
!
! THESE ROUTINES SHOULD BE MERGED 
!
  public :: ghostVpack            ! input/output: 
  public :: ghostVunpack          ! v(1-nhc:npoints+nhc,1-nhc:npoints+nhc,vlyr,ntrac_d,timelevels)
                                  
  public :: ghostVpackR           ! used to pack/unpack SPELT "Rp".  What's this?
  public :: ghostVunpackR         ! v(1-nhc:npoints+nhc,1-nhc:npoints+nhc,vlyr,ntrac_d)


! routines which do NOT include element edge data
! (used for SPELT arrays and GLL point arrays, where edge data is shared and does not need
! to be sent/received.
! these routines pack/unpack element data with user specifed halo size
!
! THESE ROUTINES CAN ALL BE REPLACED BY ghostVpack3D (if we make extra corner data arrays
! an optional argument).  Or at least these should be merged to 1 routine
  public :: ghostVpack2d          ! input/output:  
  public :: ghostVunpack2d        ! v(1-nhc:npoints+nhc,1-nhc:npoints+nhc, vlyr, ntrac_d,timelevels)
                                   
                                  ! used to pack/unpack SPELT%sga.  what's this?
  public :: ghostVpack2d_single   ! input/output
  public :: ghostVunpack2d_single !   v(1-nhc:npoints+nhc,1-nhc:npoints+nhc)
                                  
                                  ! used to pack/unpack FV vertex data (velocity/grid)
  public :: ghostVpack2d_level    ! input/output
  public :: ghostVunpack2d_level  ! v(1-nhc:npoints+nhc,1-nhc:npoints+nhc, vlyr) 



  ! Wrap pointer so we can make an array of them.
  type :: wrap_ptr
     real (kind=real_kind), dimension(:,:), pointer :: ptr => null()
  end type wrap_ptr

  type(wrap_ptr) :: edgebuff_ptrs(0:1)

contains

  subroutine initEdgeSBuffer(par,edge,elem,nlyr,bptr,rptr)
    implicit none 
    type (parallel_t), intent(in) :: par
    type (EdgeBuffer_t), intent(out), target :: edge 
    type (element_t), intent(in) :: elem(:)
    integer, intent(in)          :: nlyr
    real(kind=real_kind), optional, pointer :: bptr(:), rptr(:)
    
    !local 
    logical (kind=log_kind), parameter :: nMethod =.TRUE.

!    call initEdgeBuffer(par,edge,elem,nlyr,NewMethod=nMethod,buf_ptr=bptr,receive_ptr=rptr)
    call initEdgeBuffer(par,edge,elem,nlyr,NewMethod=nMethod)

  end subroutine initEdgeSBuffer
  ! =========================================
  ! initEdgeBuffer:
  !
  ! create an Real based communication buffer
  ! =========================================
!IDEA   subroutine initEdgeBuffer(par,edge,elem,nlyr,buf_ptr, receive_ptr, NewMethod)
  subroutine initEdgeBuffer(par,edge,elem,nlyr, NewMethod, numthreads_in )
    use dimensions_mod, only : np, nelemd, max_corner_elem
    use schedtype_mod, only : cycle_t, schedule_t, schedule
    implicit none
    type (parallel_t), intent(in) :: par
    type (EdgeBuffer_t),intent(out), target :: edge
    type (element_t),intent(in)  :: elem(:)
    integer,intent(in)                :: nlyr
    logical (kind=log_kind), intent(in), optional :: NewMethod
    integer,intent(in), optional :: numthreads_in
    !
    ! Note: this routine is now thread safe.  but 'edge' should be shared by 
    ! all threads and should be instantiated outside the threaded region
    !
    ! edge buffer must be intialized for the number of threads that will be
    ! active when calling bndry_exchange.  
    ! default is 'nthreadshoriz', which can be overriden with the optional
    ! numthreads argument.  
    !
    ! 
    ! Local variables
    integer :: nbuf,ith
    integer :: nSendCycles, nRecvCycles
    integer :: icycle, ierr
    integer :: iam,ie, i 
    integer :: edgeid,elemid
    integer :: ptr,llen,moveLength, mLen, tlen 
    integer :: numthreads
    type (Cycle_t), pointer :: pCycle
    type (Schedule_t), pointer :: pSchedule
    integer :: dest, source, length, tag, iptr
    integer :: nlen, ithr

!    call t_adj_detailf(+3)
!    call t_startf('initedgebuffer')

    if(present(NewMethod)) then 
        nbuf=nlyr*4*(1+max_corner_elem)*nelemd
    else 
        nbuf=nlyr*4*(np+max_corner_elem)*nelemd
    endif

    if (present(numthreads_in)) then
       numthreads = numthreads_in
    else
       numthreads = nthreadshoriz
    end if

! DO NOT REMOVE THIS NEXT BARRIER
! MT: This initial barrier fixes a long standing issue with Intel compilers on
! two different platforms.  Without this barrier, edge buffers initialized from
! within the threaded region would not work in a reproducable way with certain
! thread combinations.  I cant explain why, but this fixes that issue on Edison
!$OMP BARRIER


!$OMP MASTER
    edge%nlyr=nlyr
    edge%nbuf=nbuf
!$OMP END MASTER
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

    iam = par%rank
    allocate(edge%putmap(max_neigh_edges,nelemd))
    allocate(edge%getmap(max_neigh_edges,nelemd))
    allocate(edge%reverse(max_neigh_edges,nelemd))

#if 0
if(present(NewMethod)) then 
    do ie=1,nelemd
       if(iam==1) then 
          print *,'IAM: ',iam,' ie: ',ie, ' generic ORDERED putmap: ', elem(ie)%desc%putmapS
       endif
    enddo
    do ie=1,nelemd
       if(iam==1) then 
          print *,'IAM: ',iam,' ie: ',ie, ' generic ORDERED getmap: ', elem(ie)%desc%getmapS
       endif
    enddo
endif
#endif
    do ie=1,nelemd
       do i=1,max_neigh_edges
          if(elem(ie)%desc%putmapP(i) == -1) then 
              edge%putmap(i,ie) = -1
          else
              if(present(NewMethod)) then 
                  edge%putmap(i,ie) = nlyr*elem(ie)%desc%putmapS(i)
              else
                  edge%putmap(i,ie) = nlyr*elem(ie)%desc%putmapP(i)
              endif
          endif
          if(elem(ie)%desc%getmapP(i) == -1) then 
              edge%getmap(i,ie) = -1
          else
              if(present(NewMethod)) then 
                  edge%getmap(i,ie) = nlyr*elem(ie)%desc%getmapS(i)
              else
                  edge%getmap(i,ie) = nlyr*elem(ie)%desc%getmapP(i)
              endif
          endif
          edge%reverse(i,ie) = elem(ie)%desc%reverse(i) 
       enddo
    enddo

    ! Determine the most optimal way to move data in the bndry_exchange call 
    pSchedule  => Schedule(1)
    if(present(NewMethod)) then 
        moveLength = nlyr*pSchedule%MoveCycle(1)%lengthS
        ptr       = nlyr*(pSchedule%MoveCycle(1)%ptrS -1) + 1 
    else
        moveLength = nlyr*pSchedule%MoveCycle(1)%lengthP
        ptr       = nlyr*(pSchedule%MoveCycle(1)%ptrP -1) + 1 
    endif

#if 0
if(present(NewMethod)) then 
    do ie=1,nelemd
       if(iam==1) then 
          print *,'IAM: ',iam,' ie: ',ie, ' ORDERED putmap: ', edge%putmap(:,ie)
       endif
    enddo
    do ie=1,nelemd
       if(iam==1) then 
          print *,'IAM: ',iam,' ie: ',ie,' ORDERED getmap: ', edge%getmap(:,ie)
       endif
    enddo
endif
#endif
    if(numthreads<=0) then 
       nlen=1
    else
       nlen=numthreads
    endif
     
    allocate(edge%moveLength(nlen))
    allocate(edge%movePtr(nlen))

    if (numthreads > 1) then 
       ! the master thread performs no data movement because it is busy with the
       ! MPI messaging 
       edge%moveLength(1) = -1
       edge%movePtr(1) = 0
       
       ! Calculate the length of the local copy in bndy_exchange
       llen = ceiling(real(moveLength,kind=real_kind)/real(numthreads-1,kind=real_kind))
       iptr = ptr
       mLen = 0
       do i=2,numthreads
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
       edge%movePtr(1) = ptr
    endif

    ! allocate the MPI Send/Recv request handles
    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles
    allocate(edge%Srequest(nSendCycles))
    allocate(edge%Rrequest(nRecvCycles))
    allocate(edge%status(HME_status_size,nRecvCycles))

    allocate(edge%receive(nbuf))   
    allocate(edge%buf(nbuf))

!    ithr = omp_get_thread_num()+1
!    ! first touch for message buffers
!    iptr   = edge%moveptr(ithr)
!    length = edge%moveLength(ithr)
!    if(length>0) then
!        edge%buf(iptr:iptr+length-1) = 0.0D0
!        edge%receive(iptr:iptr+length-1) = 0.0D0
!    endif

!    dont do this, to improve first touch data placement
!    edge%buf    (:)=0.0D0
!    edge%receive(:)=0.0D0

!$OMP END MASTER
! MT: This next barrier is also needed - threads cannot start using edge()
! until MASTER is done initializing it
!$OMP BARRIER

!JMD DEBUGGING print statements 
!if(present(NewMethod)) then 
!    nSendCycles = pSchedule%nSendCycles
!    nRecvCycles = pSchedule%nRecvCycles
!    do icycle=1,nRecvCycles
!       pCycle => pSchedule%RecvCycle(icycle)
!       length = nlyr * pCycle%lengthS
!       iptr   = nlyr * (pCycle%ptrS - 1) + 1
!       print *,'IAM: ', iam, 'RecvCycle: Pointer: ',iptr,' LENGTH: ',length
!    enddo
!    do icycle=1,nSendCycles
!       pCycle => pSchedule%SendCycle(icycle)
!       length = nlyr * pCycle%lengthS
!!!       iptr   = nlyr * (pCycle%ptrS - 1) + 1
!       print *,'IAM: ', iam, 'SendCycle: Pointer: ',iptr,' LENGTH: ',length
!    enddo
!    iptr   = nlyr*(pSchedule%MoveCycle(1)%ptrS - 1) + 1
!    length = nlyr*pSchedule%MoveCycle(1)%lengthS
!    print *,'IAM: ',iam,'MoveCycle: Pointers: ', iptr,' LENGTH: ',length 
!endif
!    print *,'IAM: ',iam,'MoveCycle: Pointers: ',edge%movePtr  
!    print *,'IAM: ',iam,'MoveCycle: LENGTH: ',edge%moveLength  

!    stop 'initNewEdgeBuffer'
!    allocate(edge%buf(nbuf))
!    allocate(edge%receive(nbuf))
!    edge%buf    (:)=0.0D0
!    edge%receive(:)=0.0D0

#ifdef MPI_PERSISTENT
!JMD
!JMD  This is old 2D message buffer stuff
!JMD
!    pSchedule => Schedule(1)
!    nSendCycles = pSchedule%nSendCycles
!    nRecvCycles = pSchedule%nRecvCycles
!!    print *,'iam: ',iam, ' nSendCycles: ',nSendCycles, ' nRecvCycles: ',
!!    nRecvCycles
!    allocate(edge%Rrequest(nRecvCycles))
!    allocate(edge%Srequest(nSendCycles))
!    do icycle=1,nSendCycles
!       pCycle => pSchedule%SendCycle(icycle)
!       dest   = pCycle%dest -1
!       length = nlyr * pCycle%lengthP
!       tag    = pCycle%tag
!       iptr   = pCycle%ptrP
!!       print *,'IAM: ',iam, ' length: ',length,' dest: ',dest,' tag: ',tag
!       call MPI_Send_init(edge%buf(1,iptr),length,MPIreal_t,dest,tag,par%comm, edge%Srequest(icycle),ierr)
!    enddo
!    do icycle=1,nRecvCycles
!       pCycle => pSchedule%RecvCycle(icycle)
!       source   = pCycle%source -1
!       length = nlyr * pCycle%lengthP
!       tag    = pCycle%tag
!       iptr   = pCycle%ptrP
!!       print *,'IAM: ',iam, 'length: ',length,' dest: ',source,' tag: ',tag
!       call MPI_Recv_init(edge%receive(1,iptr),length,MPIreal_t,source,tag,par%comm, edge%Rrequest(icycle),ierr)
!    enddo
#endif

!    call t_stopf('initedgebuffer')
!    call t_adj_detailf(-3)

  end subroutine initEdgeBuffer
  ! =========================================
  ! initLongEdgeBuffer:
  !
  ! create an Integer based communication buffer
  ! =========================================
  subroutine initLongEdgeBuffer(edge,nlyr)
    use dimensions_mod, only : np, nelemd, max_corner_elem
    implicit none
    integer,intent(in)                :: nlyr
    type (LongEdgeBuffer_t),intent(out) :: edge

    ! Local variables

    integer :: nbuf

    ! sanity check for threading
    if (omp_get_num_threads()>1) then
       call haltmp('ERROR: initLongEdgeBuffer must be called before threaded reagion')
    endif

    nbuf=4*(np+max_corner_elem)*nelemd
    edge%nlyr=nlyr
    edge%nbuf=nbuf
    allocate(edge%buf(nlyr,nbuf))
    edge%buf(:,:)=0

    allocate(edge%receive(nlyr,nbuf))
    edge%receive(:,:)=0

  end subroutine initLongEdgeBuffer
  ! =========================================
  ! edgeDGVpack:
  !
  ! Pack edges of v into buf for DG stencil
  ! =========================================
  subroutine edgeDGVpack(edge,v,vlyr,kptr,ielem)
    use dimensions_mod, only : np
    type (EdgeBuffer_t)                      :: edge
    integer,              intent(in)   :: vlyr
    real (kind=real_kind),intent(in)   :: v(np,np,vlyr)
    integer,              intent(in)   :: kptr
    integer,              intent(in)   :: ielem

    ! =========================================
    ! This code is just a wrapper call the 
    !   normal oldedgeVpack
    ! =========================================
    call edgeVpack(edge,v,vlyr,kptr,ielem)

  end subroutine edgeDGVpack

  ! ===========================================
  !  FreeEdgeBuffer:
  !
  !  Freed an edge communication buffer
  ! =========================================
  subroutine FreeEdgeBuffer(edge)
    implicit none
    type (EdgeBuffer_t),intent(inout) :: edge

#if (defined HORIZ_OPENMP)
!$OMP BARRIER
!$OMP MASTER
#endif
    deallocate(edge%buf)
    deallocate(edge%receive)
    deallocate(edge%putmap)
    deallocate(edge%getmap)
    deallocate(edge%reverse)
#if (defined HORIZ_OPENMP)
!$OMP END MASTER
#endif

  end subroutine FreeEdgeBuffer

  subroutine FreeGhostBuffer3D(buffer) 
    use edgetype_mod, only : ghostbuffer3d_t 
    implicit none
    type (Ghostbuffer3d_t),intent(inout) :: buffer

#if (defined HORIZ_OPENMP)
!$OMP BARRIER
!$OMP MASTER
#endif
    buffer%nbuf=0
    buffer%nlyr=0
    deallocate(buffer%buf)
    deallocate(buffer%receive)
#if (defined HORIZ_OPENMP)
!$OMP END MASTER
#endif

  end subroutine FreeGhostBuffer3D
  ! ===========================================
  !  FreeLongEdgeBuffer:
  !
  !  Freed an edge communication buffer
  ! =========================================
  subroutine FreeLongEdgeBuffer(edge) 
    implicit none
    type (LongEdgeBuffer_t),intent(inout) :: edge

    edge%nbuf=0
    edge%nlyr=0
    deallocate(edge%buf)
    deallocate(edge%receive)

  end subroutine FreeLongEdgeBuffer

  ! =========================================
  !
  !> @brief Pack edges of v into an edge buffer for boundary exchange.
  !
  !> This subroutine packs for one or more vertical layers into an edge 
  !! buffer. If the buffer associated with edge is not large enough to 
  !! hold all vertical layers you intent to pack, the method will 
  !! halt the program with a call to parallel_mod::haltmp().
  !! @param[in] edge Edge Buffer into which the data will be packed.
  !! This buffer must be previously allocated with initEdgeBuffer().
  !! @param[in] v The data to be packed.
  !! @param[in] vlyr Number of vertical level coming into the subroutine
  !! for packing for input v.
  !! @param[in] kptr Vertical pointer to the place in the edge buffer where 
  !! data will be located.
  ! =========================================
  subroutine edgeVpack(edge,v,vlyr,kptr,ielem)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t)                :: edge
    integer,              intent(in)   :: vlyr
    real (kind=real_kind),intent(in)   :: v(np,np,vlyr)
    integer,              intent(in)   :: kptr
    integer,              intent(in)   :: ielem
!    type (EdgeDescriptor_t),intent(in) :: desc

    ! Local variables
    integer :: i,k,ir,ll,llval,iptr

    integer :: is,ie,in,iw

    !call t_adj_detailf(+2)
    !call t_startf('edgeVpack')

    is = edge%putmap(south,ielem)
    ie = edge%putmap(east,ielem)
    in = edge%putmap(north,ielem)
    iw = edge%putmap(west,ielem)
    if (edge%nlyr < (kptr+vlyr) ) then
       print *,'edge%nlyr = ',edge%nlyr
       print *,'kptr+vlyr = ',kptr+vlyr
       call haltmp('edgeVpack: Buffer overflow: size of the vertical dimension must be increased!')
    endif

!dir$ ivdep
    do k=1,vlyr
       iptr = np*(kptr+k-1)
       do i=1,np
          edge%buf(iptr+is+i)   = v(i  ,1 ,k) ! South
          edge%buf(iptr+in+i)   = v(i  ,np,k) ! North
          edge%buf(iptr+iw+i)   = v(1  ,i ,k) ! West
          edge%buf(iptr+ie+i)   = v(np ,i ,k) ! East
       enddo
    enddo

    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing
    if(edge%reverse(south,ielem)) then
!dir$ ivdep
       do k=1,vlyr
          iptr = np*(kptr+k-1)+is
          do i=1,np
             edge%buf(iptr+np-i+1)=v(i,1,k)
          enddo
       enddo
    endif

    if(edge%reverse(east,ielem)) then
!dir$ ivdep
       do k=1,vlyr
          iptr=np*(kptr+k-1)+ie
          do i=1,np
             edge%buf(iptr+np-i+1)=v(np,i,k)
          enddo
       enddo
    endif

    if(edge%reverse(north,ielem)) then
!dir$ ivdep
       do k=1,vlyr
          iptr=np*(kptr+k-1)+in
          do i=1,np
             edge%buf(iptr+np-i+1)=v(i,np,k)
          enddo
       enddo
    endif

    if(edge%reverse(west,ielem)) then
!dir$ ivdep
       do k=1,vlyr
          iptr=np*(kptr+k-1)+iw
          do i=1,np
             edge%buf(iptr+np-i+1)=v(1,i,k)
          enddo
       enddo
    endif

! SWEST
    do ll=swest,swest+max_corner_elem-1
        llval=edge%putmap(ll,ielem)
        if (llval /= -1) then
!dir$ ivdep
            do k=1,vlyr
                edge%buf(kptr+k+llval)=v(1  ,1 ,k)
            end do
        end if
    end do

! SEAST
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
        llval=edge%putmap(ll,ielem)
        if (llval /= -1) then
!dir$ ivdep
            do k=1,vlyr
                edge%buf(kptr+k+llval)=v(np ,1 ,k)
            end do
        end if
    end do

! NEAST
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        llval=edge%putmap(ll,ielem)
        if (llval /= -1) then
!dir$ ivdep
            do k=1,vlyr
                edge%buf(kptr+k+llval)=v(np ,np,k)
            end do
        end if
    end do

! NWEST
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        llval=edge%putmap(ll,ielem)
        if (llval /= -1) then
!dir$ ivdep
            do k=1,vlyr
                edge%buf(kptr+k+llval)=v(1  ,np,k)
            end do
        end if
    end do

    !call t_stopf('edgeVpack')
    !call t_adj_detailf(-2)

  end subroutine edgeVpack

  subroutine edgeSpack(edge,v,vlyr,kptr,ielem)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t)                :: edge
    integer,              intent(in)   :: vlyr
    real (kind=real_kind),intent(in)   :: v(vlyr)
    integer,              intent(in)   :: kptr
    integer,              intent(in)   :: ielem
!    type (EdgeDescriptor_t),intent(in) :: desc

    ! Local variables
    integer :: i,k,ir,ll,llval,iptr

    integer :: is,ie,in,iw
    real (kind=real_kind) :: tmp

!pw call t_adj_detailf(+2)
!pw call t_startf('edgeSpack')

    is = edge%putmap(south,ielem)
    ie = edge%putmap(east,ielem)
    in = edge%putmap(north,ielem)
    iw = edge%putmap(west,ielem)
    if (edge%nlyr < (kptr+vlyr) ) then
       call haltmp('edgeSpack: Buffer overflow: size of the vertical dimension must be increased!')
    endif

!dir$ ivdep
    do k=1,vlyr
       edge%buf(kptr+k+ie) = v(k) ! East
       edge%buf(kptr+k+is) = v(k) ! South
       edge%buf(kptr+k+in) = v(k) ! North
       edge%buf(kptr+k+iw) = v(k) ! West
    enddo

! SWEST
    do ll=swest,swest+max_corner_elem-1
        llval=edge%putmap(ll,ielem)
        if (llval /= -1) then
!dir$ ivdep
            do k=1,vlyr
                edge%buf(kptr+k+llval)=v(k)
            end do
        end if
    end do

! SEAST
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
        llval=edge%putmap(ll,ielem)
        if (llval /= -1) then
!dir$ ivdep
            do k=1,vlyr
                edge%buf(kptr+k+llval)=v(k)
            end do
        end if
    end do

! NEAST
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        llval=edge%putmap(ll,ielem)
        if (llval /= -1) then
!dir$ ivdep
            do k=1,vlyr
                edge%buf(kptr+k+llval)=v(k)
            end do
        end if
    end do

! NWEST
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        llval=edge%putmap(ll,ielem)
        if (llval /= -1) then
!dir$ ivdep
            do k=1,vlyr
                edge%buf(kptr+k+llval)=v(k)
            end do
        end if
    end do

!pw call t_stopf('edgeSpack')
!pw call t_adj_detailf(-2)

  end subroutine edgeSpack

  ! =========================================
  ! LongEdgeVpack:
  !
  ! Pack edges of v into buf...
  ! =========================================
  subroutine LongEdgeVpack(edge,v,vlyr,kptr,desc)
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use dimensions_mod, only : np, max_corner_elem

    type (LongEdgeBuffer_t)            :: edge
    integer,              intent(in)   :: vlyr
    integer (kind=int_kind),intent(in)   :: v(np,np,vlyr)
    integer,              intent(in)   :: kptr
    type (EdgeDescriptor_t),intent(in) :: desc

    ! Local variables
    integer :: i,k,ir,l
    integer :: is,ie,in,iw

    if(.not. threadsafe) then
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
       threadsafe=.true.
    end if

    is = desc%putmapP(south)
    ie = desc%putmapP(east)
    in = desc%putmapP(north)
    iw = desc%putmapP(west)

!dir$ ivdep
    do k=1,vlyr
!dir$ ivdep
       do i=1,np
          edge%buf(kptr+k,is+i) = v(i  ,1 ,k)
          edge%buf(kptr+k,ie+i) = v(np ,i ,k)
          edge%buf(kptr+k,in+i) = v(i  ,np,k)
          edge%buf(kptr+k,iw+i) = v(1  ,i ,k)
       enddo
    end do

    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing
    if(desc%reverse(south)) then
!dir$ ivdep
       do k=1,vlyr
!dir$ ivdep
          do i=1,np
             edge%buf(kptr+k,is+np-i+1)=v(i,1,k)
          enddo
       enddo
    endif

    if(desc%reverse(east)) then
!dir$ ivdep
       do k=1,vlyr
!dir$ ivdep
          do i=1,np
             edge%buf(kptr+k,ie+np-i+1)=v(np,i,k)
          enddo
       enddo
    endif

    if(desc%reverse(north)) then
!dir$ ivdep
       do k=1,vlyr
!dir$ ivdep
          do i=1,np
             edge%buf(kptr+k,in+np-i+1)=v(i,np,k)
          enddo
       enddo
    endif

    if(desc%reverse(west)) then
!dir$ ivdep
       do k=1,vlyr
!dir$ ivdep
          do i=1,np
             edge%buf(kptr+k,iw+np-i+1)=v(1,i,k)
          enddo
       enddo
    endif

! SWEST
    do l=swest,swest+max_corner_elem-1
        if (desc%putmapP(l) /= -1) then
!dir$ ivdep
            do k=1,vlyr
                edge%buf(kptr+k,desc%putmapP(l)+1)=v(1  ,1 ,k)
            end do
        end if
    end do

! SEAST
    do l=swest+max_corner_elem,swest+2*max_corner_elem-1
        if (desc%putmapP(l) /= -1) then
!dir$ ivdep
            do k=1,vlyr
                edge%buf(kptr+k,desc%putmapP(l)+1)=v(np ,1 ,k)
            end do
        end if
    end do

! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        if (desc%putmapP(l) /= -1) then
!dir$ ivdep
            do k=1,vlyr
                edge%buf(kptr+k,desc%putmapP(l)+1)=v(np ,np,k)
            end do
        end if
    end do

! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        if (desc%putmapP(l) /= -1) then
!dir$ ivdep
            do k=1,vlyr
                edge%buf(kptr+k,desc%putmapP(l)+1)=v(1  ,np,k)
            end do
        end if
    end do

  end subroutine LongEdgeVpack

  subroutine edgeVunpack(edge,v,vlyr,kptr,ielem)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    type (EdgeBuffer_t),         intent(in)  :: edge

    integer,               intent(in)  :: vlyr
    real (kind=real_kind), intent(inout) :: v(np,np,vlyr)
!$dir assume_aligned v:64
    integer,               intent(in)  :: kptr
    integer,               intent(in)  :: ielem
    !type (EdgeDescriptor_t)            :: desc

    ! Local
    integer :: i,k,ll,iptr
    integer :: is,ie,in,iw
    integer :: ks,ke,kblock
    logical :: done
    integer :: getmapL

    !call t_adj_detailf(+2)
    !call t_startf('edgeVunpack')

    is=edge%getmap(south,ielem)
    ie=edge%getmap(east,ielem)
    in=edge%getmap(north,ielem)
    iw=edge%getmap(west,ielem)

!dir$ ivdep
    do k=1,vlyr
       iptr=np*(kptr+k-1)
       do i=1,np
          v(i  ,1  ,k) = v(i  ,1  ,k)+edge%receive(iptr+is+i) ! South
          v(i  ,np ,k) = v(i  ,np ,k)+edge%receive(iptr+in+i) ! North
          v(1  ,i  ,k) = v(1  ,i  ,k)+edge%receive(iptr+iw+i) ! West
          v(np ,i  ,k) = v(np ,i  ,k)+edge%receive(iptr+ie+i) ! East
       enddo
    enddo

! SWEST
    do ll=swest,swest+max_corner_elem-1
        getmapL = edge%getmap(ll,ielem)
        if(getmapL /= -1) then 
!dir$ ivdep
            do k=1,vlyr
                v(1  ,1 ,k)=v(1 ,1 ,k)+edge%receive((kptr+k-1)+getmapL+1)
            enddo
        endif
    end do

! SEAST
    do ll=swest+max_corner_elem,swest+2*max_corner_elem-1
        getmapL = edge%getmap(ll,ielem)
        if(getmapL /= -1) then 
!dir$ ivdep
            do k=1,vlyr
                v(np ,1 ,k)=v(np,1 ,k)+edge%receive((kptr+k-1)+getmapL+1)
            enddo
        endif
    end do

! NEAST
    do ll=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        getmapL = edge%getmap(ll,ielem)
        if(getmapL /= -1) then 
!dir$ ivdep
            do k=1,vlyr
                v(np ,np,k)=v(np,np,k)+edge%receive((kptr+k-1)+getmapL+1)
            enddo
        endif
    end do

! NWEST
    do ll=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        getmapL = edge%getmap(ll,ielem)
        if(getmapL /= -1) then 
!dir$ ivdep
            do k=1,vlyr
                v(1  ,np,k)=v(1 ,np,k)+edge%receive((kptr+k-1)+getmapL+1)
            enddo
        endif
    end do

    !call t_stopf('edgeVunpack')
    !call t_adj_detailf(-2)

  end subroutine edgeVunpack


  subroutine edgeVunpackVert(edge,v,ielem)
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use dimensions_mod, only : np, max_corner_elem, ne
    use coordinate_systems_mod, only: cartesian3D_t

    type (EdgeBuffer_t),   intent(inout)  :: edge
    type (cartesian3D_t), intent(out) :: v(:,:,:)
    integer, intent(in) :: ielem

    ! Local
    logical, parameter :: UseUnroll = .TRUE.
    integer :: i,k,l, nce
    integer :: is,ie,in,iw,ine,inw,isw,ise

    threadsafe=.false.

    if (max_corner_elem.ne.1 .and. ne==0) then
        ! MNL: this is used to construct the dual grid on the cube,
        !      currently only supported for the uniform grid. If
        !      this is desired on a refined grid, a little bit of
        !      work will be required.
        call haltmp("edgeVunpackVert should not be called with unstructured meshes")
    end if

    is=edge%getmap(south,ielem)
    ie=edge%getmap(east,ielem)
    in=edge%getmap(north,ielem)
    iw=edge%getmap(west,ielem)


    ! N+S
    do i=1,np/2
       ! North
       ! kptr = 0
       ! iptr=np*(kptr+k-1)+in
  
!      iptr = np*(1-1)+in 
!      v(3,i ,np)%x = edge%receive(1,in+i) 
       v(3,i ,np)%x = edge%receive(in+i) 

!      iptr = np*(2-1)+in 
!      v(3,i ,np)%y = edge%receive(2,in+i) 
       v(3,i ,np)%y = edge%receive(np+in+i) 

!      iptr = np*(3-1)+in 
!      v(3,i ,np)%z = edge%receive(3,in+i) 
       v(3,i ,np)%z = edge%receive(2*np+in+i) 

       ! South
       ! v(2,i ,1)%x  = edge%receive(1,is+i) 
       v(2,i ,1)%x  = edge%receive(is+i) 

       ! v(2,i ,1)%y  = edge%receive(2,is+i) 
       v(2,i ,1)%y  = edge%receive(np+is+i) 

       ! v(2,i ,1)%z  = edge%receive(3,is+i) 
       v(2,i ,1)%z  = edge%receive(2*np+is+i) 
    enddo

    do i=np/2+1,np
       ! North
       ! v(4,i ,np)%x = edge%receive(1,in+i) 
       v(4,i ,np)%x = edge%receive(in+i) 

       ! v(4,i ,np)%y = edge%receive(2,in+i) 
       v(4,i ,np)%y = edge%receive(np+in+i) 

       ! v(4,i ,np)%z = edge%receive(3,in+i) 
       v(4,i ,np)%z = edge%receive(2*np+in+i) 
       ! South
       ! v(1,i ,1)%x  = edge%receive(1,is+i) 
       v(1,i ,1)%x  = edge%receive(is+i) 
       ! v(1,i ,1)%y  = edge%receive(2,is+i) 
       v(1,i ,1)%y  = edge%receive(np+is+i) 
       ! v(1,i ,1)%z  = edge%receive(3,is+i)        
       v(1,i ,1)%z  = edge%receive(2*np+is+i)        
    enddo

    do i=1,np/2
       ! East
       ! v(3,np,i)%x = edge%receive(1,ie+i)
       v(3,np,i)%x = edge%receive(ie+i)
       ! v(3,np,i)%y = edge%receive(2,ie+i)
       v(3,np,i)%y = edge%receive(np+ie+i)
       ! v(3,np,i)%z = edge%receive(3,ie+i)       
       v(3,np,i)%z = edge%receive(2*np+ie+i)       
       ! West
       ! v(4,1,i)%x  = edge%receive(1,iw+i)
       v(4,1,i)%x  = edge%receive(iw+i)
       ! v(4,1,i)%y  = edge%receive(2,iw+i)
       v(4,1,i)%y  = edge%receive(np+iw+i)
       ! v(4,1,i)%z  = edge%receive(3,iw+i)
       v(4,1,i)%z  = edge%receive(2*np+iw+i)
    end do

    do i=np/2+1,np
       ! East
       ! v(2,np,i)%x = edge%receive(1,ie+i)
       v(2,np,i)%x = edge%receive(ie+i)

       ! v(2,np,i)%y = edge%receive(2,ie+i)
       v(2,np,i)%y = edge%receive(np+ie+i)

       ! v(2,np,i)%z = edge%receive(3,ie+i)       
       v(2,np,i)%z = edge%receive(2*np+ie+i)       
       ! West

       ! v(1,1,i)%x  = edge%receive(1,iw+i)
       v(1,1,i)%x  = edge%receive(iw+i)

       ! v(1,1,i)%y  = edge%receive(2,iw+i)
       v(1,1,i)%y  = edge%receive(np+iw+i)

       ! v(1,1,i)%z  = edge%receive(3,iw+i)
       v(1,1,i)%z  = edge%receive(2*np+iw+i)
    end do

! SWEST
    nce = max_corner_elem
    do l=swest,swest+max_corner_elem-1
       ! find the one active corner, then exist
        isw=edge%getmap(l,ielem)
        if(isw /= -1) then 
            ! v(1,1,1)%x=edge%receive(1,desc%getmapP(l)+1)
            v(1,1,1)%x=edge%receive(isw+1)
            ! v(1,1,1)%y=edge%receive(2,desc%getmapP(l)+1)
            v(1,1,1)%y=edge%receive(nce+isw+1)
            ! v(1,1,1)%z=edge%receive(3,desc%getmapP(l)+1)
            v(1,1,1)%z=edge%receive(2*nce+isw+1)
            exit 
        else
            v(1,1,1)%x=0_real_kind
            v(1,1,1)%y=0_real_kind
            v(1,1,1)%z=0_real_kind
        endif
    end do

! SEAST
    do l=swest+max_corner_elem,swest+2*max_corner_elem-1
       ! find the one active corner, then exist
        ise=edge%getmap(l,ielem)
        if(ise /= -1) then 
            ! v(2,np,1)%x=edge%receive(1,desc%getmapP(l)+1)
            v(2,np,1)%x=edge%receive(ise+1)
            ! v(2,np,1)%y=edge%receive(2,desc%getmapP(l)+1)
            v(2,np,1)%y=edge%receive(nce+ise+1)
            ! v(2,np,1)%z=edge%receive(3,desc%getmapP(l)+1)
            v(2,np,1)%z=edge%receive(2*nce+ise+1)
            exit
        else
            v(2,np,1)%x=0_real_kind
            v(2,np,1)%y=0_real_kind
            v(2,np,1)%z=0_real_kind
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
            v(3,np,np)%x=0_real_kind
            v(3,np,np)%y=0_real_kind
            v(3,np,np)%z=0_real_kind
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
            v(4,1,np)%x=0_real_kind
            v(4,1,np)%y=0_real_kind
            v(4,1,np)%z=0_real_kind
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
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t),         intent(in)  :: edge
    integer,               intent(in)  :: vlyr
    real (kind=real_kind), intent(out) :: v(0:np+1,0:np+1,vlyr)
    integer,               intent(in)  :: kptr
    integer,               intent(in)  :: ielem

    ! Local
    integer :: i,k,iptr,nce
    integer :: is,ie,in,iw

    threadsafe=.false.

    is=edge%getmap(south,ielem)
    ie=edge%getmap(east,ielem)
    in=edge%getmap(north,ielem)
    iw=edge%getmap(west,ielem)
!dir$ ivdep
    do k=1,vlyr
       iptr=np*(kptr+k-1)
       do i=1,np
          v(i   ,0   ,k)=edge%receive(iptr+is+i)
          v(np+1,i   ,k)=edge%receive(iptr+ie+i)
          v(i   ,np+1,k)=edge%receive(iptr+in+i)
          v(0   ,i   ,k)=edge%receive(iptr+iw+i)
       end do
    end do

    nce = max_corner_elem
!   this is probably broken.  nce should be 1?  MT 2016/2/9
    i = swest
    if(edge%getmap(i,ielem) /= -1) then
!dir$ ivdep
      do k=1,vlyr
        v(0,0,k) = edge%receive(nce*(kptr+k-1)+edge%getmap(i,ielem)+1)
      end do
    end if
    i = swest+max_corner_elem
    if(edge%getmap(i,ielem) /= -1) then
!dir$ ivdep
      do k=1,vlyr
        v(np+1,0,k) = edge%receive(nce*(kptr+k-1)+edge%getmap(i,ielem)+1)
      end do
    end if
    i = swest+3*max_corner_elem
    if(edge%getmap(i,ielem) /= -1) then
!dir$ ivdep
      do k=1,vlyr
        v(np+1,np+1,k) = edge%receive(nce*(kptr+k-1)+edge%getmap(i,ielem)+1)
      end do
    end if
    i = swest+2*max_corner_elem
    if(edge%getmap(i,ielem) /= -1) then
!dir$ ivdep
      do k=1,vlyr
        v(0,np+1,k) = edge%receive(nce*(kptr+k-1)+edge%getmap(i,ielem)+1)
      end do
    end if

  end subroutine edgeDGVunpack

  ! ========================================
  ! edgeVunpackMIN/MAX:
  !
  ! Finds the Min/Max edges from edge buffer into v...
  ! ========================================
  subroutine edgeVunpackMAX(edge,v,vlyr,kptr,ielem)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t),         intent(in)  :: edge
    integer,               intent(in)  :: vlyr
    real (kind=real_kind), intent(inout) :: v(np,np,vlyr)
    integer,               intent(in)  :: kptr
    integer,               intent(in)  :: ielem

    ! Local
    integer :: i,k,l,iptr
    integer :: is,ie,in,iw
    integer :: getmapL

    threadsafe=.false.

    is=edge%getmap(south,ielem)
    ie=edge%getmap(east,ielem)
    in=edge%getmap(north,ielem)
    iw=edge%getmap(west,ielem)
!dir$ ivdep
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
        getmapL = edge%getmap(l,ielem)
        if(getmapL /= -1) then 
!dir$ ivdep
            do k=1,vlyr
                v(1  ,1 ,k)=MAX(v(1 ,1 ,k),edge%receive(kptr+k+getmapL))
            enddo
        endif
    end do

! SEAST
    do l=swest+max_corner_elem,swest+2*max_corner_elem-1
        getmapL = edge%getmap(l,ielem)
        if(getmapL /= -1) then 
!dir$ ivdep
            do k=1,vlyr
                v(np ,1 ,k)=MAX(v(np,1 ,k),edge%receive(kptr+k+getmapL))
            enddo
        endif
    end do

! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        getmapL = edge%getmap(l,ielem)
        if(getmapL /= -1) then
!dir$ ivdep
            do k=1,vlyr
                v(np ,np,k)=MAX(v(np,np,k),edge%receive(kptr+k+getmapL))
            enddo
        endif
    end do

! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        getmapL = edge%getmap(l,ielem)
        if(getmapL /= -1) then 
!dir$ ivdep
            do k=1,vlyr
                v(1  ,np,k)=MAX(v(1 ,np,k),edge%receive(kptr+k+getmapL))
            enddo
        endif
    end do
    
  end subroutine edgeVunpackMAX

  subroutine edgeSunpackMAX(edge,v,vlyr,kptr,ielem)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t),         intent(in)  :: edge
    integer,               intent(in)  :: vlyr
    real (kind=real_kind), intent(inout) :: v(vlyr)
    integer,               intent(in)  :: kptr
    integer,               intent(in)  :: ielem

    ! Local
    integer :: i,k,l,iptr
    integer :: is,ie,in,iw
    integer :: getmapL

!pw call t_startf('edgeSunpack')
    threadsafe=.false.

    is=edge%getmap(south,ielem)
    ie=edge%getmap(east,ielem)
    in=edge%getmap(north,ielem)
    iw=edge%getmap(west,ielem)
!dir$ ivdep
    do k=1,vlyr
       iptr=(kptr+k-1)
       v(k) = MAX(v(k),edge%receive(iptr+is+1),edge%receive(iptr+ie+1),edge%receive(iptr+in+1),edge%receive(iptr+iw+1))
    end do

! SWEST
    do l=swest,swest+max_corner_elem-1
        getmapL = edge%getmap(l,ielem)
        if(getmapL /= -1) then 
!dir$ ivdep
            do k=1,vlyr
                iptr = (kptr+k-1)
                v(k)=MAX(v(k),edge%receive(kptr+k+getmapL))
            enddo
        endif
    end do

! SEAST
    do l=swest+max_corner_elem,swest+2*max_corner_elem-1
        getmapL = edge%getmap(l,ielem)
        if(getmapL /= -1) then 
!dir$ ivdep
            do k=1,vlyr
                iptr = (kptr+k-1)
                v(k)=MAX(v(k),edge%receive(kptr+k+getmapL))
            enddo
        endif
    end do

! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        getmapL = edge%getmap(l,ielem)
        if(getmapL /= -1) then
!dir$ ivdep
            do k=1,vlyr
                iptr = (kptr+k-1)
                v(k)=MAX(v(k),edge%receive(kptr+k+getmapL))
            enddo
        endif
    end do

! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        getmapL = edge%getmap(l,ielem)
        if(getmapL /= -1) then 
!dir$ ivdep
            do k=1,vlyr
                iptr = (kptr+k-1)
                v(k)=MAX(v(k),edge%receive(kptr+k+getmapL))
            enddo
        endif
    end do
!pw call t_stopf('edgeSunpack')
    
  end subroutine edgeSunpackMAX

  subroutine edgeSunpackMIN(edge,v,vlyr,kptr,ielem)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t),         intent(in)  :: edge
    integer,               intent(in)  :: vlyr
    real (kind=real_kind), intent(inout) :: v(vlyr)
    integer,               intent(in)  :: kptr
    integer,               intent(in)  :: ielem


    ! Local

    integer :: i,k,l,iptr
    integer :: is,ie,in,iw
    integer :: getmapL

!pw call t_startf('edgeSunpack')
    threadsafe=.false.

    is=edge%getmap(south,ielem)
    ie=edge%getmap(east,ielem)
    in=edge%getmap(north,ielem)
    iw=edge%getmap(west,ielem)
!dir$ ivdep
    do k=1,vlyr
       iptr=(kptr+k-1)
       v(k) = MIN(v(k),edge%receive(iptr+is+1),edge%receive(iptr+ie+1),edge%receive(iptr+in+1),edge%receive(iptr+iw+1))
    end do

! SWEST
    do l=swest,swest+max_corner_elem-1
        getmapL = edge%getmap(l,ielem)
        if(getmapL /= -1) then 
!dir$ ivdep
            do k=1,vlyr
                v(k)=MiN(v(k),edge%receive(kptr+k+getmapL))
            enddo
        endif
    end do

! SEAST
    do l=swest+max_corner_elem,swest+2*max_corner_elem-1
        getmapL = edge%getmap(l,ielem)
        if(getmapL /= -1) then 
!dir$ ivdep
            do k=1,vlyr
                v(k)=MIN(v(k),edge%receive(kptr+k+getmapL))
            enddo
        endif
    end do

! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        getmapL = edge%getmap(l,ielem)
        if(getmapL /= -1) then
!dir$ ivdep
            do k=1,vlyr
                v(k)=MIN(v(k),edge%receive(kptr+k+getmapL))
            enddo
        endif
    end do

! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        getmapL = edge%getmap(l,ielem)
        if(getmapL /= -1) then 
!dir$ ivdep
            do k=1,vlyr
                v(k)=MIN(v(k),edge%receive(kptr+k+getmapL))
            enddo
        endif
    end do
!pw call t_stopf('edgeSunpack')
    
  end subroutine edgeSunpackMIN

  subroutine edgeVunpackMIN(edge,v,vlyr,kptr,ielem)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest

    type (EdgeBuffer_t),         intent(in)  :: edge
    integer,               intent(in)  :: vlyr
    real (kind=real_kind), intent(inout) :: v(np,np,vlyr)
    integer,               intent(in)  :: kptr
    integer,               intent(in)  :: ielem

    ! Local
    integer :: i,k,l,iptr
    integer :: is,ie,in,iw
    integer :: getmapL

    threadsafe=.false.

    is=edge%getmap(south,ielem)
    ie=edge%getmap(east,ielem)
    in=edge%getmap(north,ielem)
    iw=edge%getmap(west,ielem)
!dir$ ivdep
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
        getmapL = edge%getmap(l,ielem)
        if(getmapL /= -1) then 
!dir$ ivdep
            do k=1,vlyr
                v(1  ,1 ,k)=MIN(v(1 ,1 ,k),edge%receive(kptr+k+getmapL))
            enddo
        endif
    end do

! SEAST
    do l=swest+max_corner_elem,swest+2*max_corner_elem-1
        getmapL = edge%getmap(l,ielem)
        if(getmapL /= -1) then 
!dir$ ivdep
            do k=1,vlyr
                v(np ,1 ,k)=MIN(v(np,1 ,k),edge%receive(kptr+k+getmapL))
            enddo
        endif
    end do

! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        getmapL = edge%getmap(l,ielem)
        if(getmapL /= -1) then 
!dir$ ivdep
            do k=1,vlyr
                v(np ,np,k)=MIN(v(np,np,k),edge%receive(kptr+k+getmapL))
            enddo
        endif
    end do

! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        getmapL = edge%getmap(l,ielem)
        if(getmapL /= -1) then 
!dir$ ivdep
            do k=1,vlyr
                v(1  ,np,k)=MIN(v(1 ,np,k),edge%receive(kptr+k+getmapL))
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
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use dimensions_mod, only : np, max_corner_elem

    type (LongEdgeBuffer_t),         intent(in)  :: edge
    integer,               intent(in)  :: vlyr
    integer (kind=int_kind), intent(inout) :: v(np,np,vlyr)
    integer,               intent(in)  :: kptr
    type (EdgeDescriptor_t),intent(in) :: desc

    ! Local
    integer :: i,k,l
    integer :: is,ie,in,iw
    integer :: getmapL

    threadsafe=.false.

    is=desc%getmapP(south)
    ie=desc%getmapP(east)
    in=desc%getmapP(north)
    iw=desc%getmapP(west)
!dir$ ivdep
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
        getmapL = desc%getmapP(l)
        if(getmapL /= -1) then 
!dir$ ivdep
            do k=1,vlyr
                v(1  ,1 ,k)=MIN(v(1 ,1 ,k),edge%buf(kptr+k,getmapL+1))
            enddo
        endif
    end do

! SEAST
    do l=swest+max_corner_elem,swest+2*max_corner_elem-1
        getmapL = desc%getmapP(l)
        if(getmapL /= -1) then 
!dir$ ivdep
            do k=1,vlyr
                v(np ,1 ,k)=MIN(v(np,1 ,k),edge%buf(kptr+k,getmapL+1))
            enddo
        endif
    end do

! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
        getmapL = desc%getmapP(l)
        if(getmapL /= -1) then 
!dir$ ivdep
            do k=1,vlyr
                v(np ,np,k)=MIN(v(np,np,k),edge%buf(kptr+k,getmapL+1))
            enddo
        endif
    end do

! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
        getmapL = desc%getmapP(l)
        if(getmapL /= -1) then 
!dir$ ivdep
            do k=1,vlyr
                v(1  ,np,k)=MIN(v(1 ,np,k),edge%buf(kptr+k,getmapL+1))
            enddo
        endif
    end do

  end subroutine LongEdgeVunpackMIN

  ! =============================
  ! edgerotate:
  !
  ! Rotate edges in buffer...
  ! =============================
  subroutine edgerotate(edge,vlyr,kptr,desc)
    use dimensions_mod, only : np
    type (EdgeBuffer_t)           :: edge         ! edge struct
    integer, intent(in)           :: vlyr         ! number of 2d vector fields to rotate
    integer, intent(in)           :: kptr         ! layer pointer into edge buffer
    type (EdgeDescriptor_t), intent(in) :: desc

    ! Local variables

    integer :: i,k,k1,k2
    integer :: irot,ia,nbr

    real(kind=real_kind), dimension(:,:,:), pointer :: R
    real(kind=real_kind)  :: tmp1,tmp2

    print *,'entered into edgerotate... Note that this code is currently not functional' 
    stop 'ERROR: edgerotate'
#ifdef _USEASSOCIATED
    if (associated(rot)) then
#else
       if (desc%use_rotation == 1) then
#endif

          do irot=1,SIZE(desc%rot)

             nbr  =  desc%rot(irot)%nbr
             R    => desc%rot(irot)%R

             ia=desc%putmapP(nbr)

             ! ========================================
             ! If nbr direction is (1-4) => is an edge
             ! ========================================

             if (nbr <= 4) then

                ! ========================================================
                !  Is an edge. Rotate it in place
                ! ========================================================

!                do i=1,np
!                   do k=1,vlyr,2
!                      k1 = kptr + k
!                      k2 = k1 + 1
!                      tmp1=R(1,1,i)*edge%buf(k1,ia+i) + R(1,2,i)*edge%buf(k2,ia+i)
!                      tmp2=R(2,1,i)*edge%buf(k1,ia+i) + R(2,2,i)*edge%buf(k2,ia+i)
!                      edge%buf(k1,ia+i)=tmp1
!                      edge%buf(k2,ia+i)=tmp2
!                   end do
!                end do

             else

                ! ===================================================
                ! Is an element corner point, but not a cube corner
                ! point, just rotate it in place.
                ! ===================================================

!                if (ia /= -1) then
!                   do k=1,vlyr,2
!                      k1 = kptr + k
!                      k2 = k1+1
!                      tmp1=R(1,1,1)*edge%buf(k1,ia+1) + R(1,2,1)*edge%buf(k2,ia+1)
!                      tmp2=R(2,1,1)*edge%buf(k1,ia+1) + R(2,2,1)*edge%buf(k2,ia+1)
!                      edge%buf(k1,ia+1)=tmp1
!                      edge%buf(k2,ia+1)=tmp2
!                   end do
!                end if

             end if

          end do

       endif

     end subroutine edgerotate

     ! =============================================
     ! buffermap:
     !
     ! buffermap translates element number, inum and
     ! element edge/corner, facet, into an edge buffer 
     ! memory location, loc.
     ! =============================================

     function buffermap(inum,facet) result(loc)
       use dimensions_mod, only : np
       integer, intent(in) :: inum   
       integer, intent(in) :: facet
       integer :: loc

       if (facet>4) then
          if (inum == -1) then
             loc = inum
          else
             loc=(inum-1)*(4*np+4)+4*np+(facet-5)
          end if
       else
          loc=(inum-1)*(4*np+4)+np*(facet-1)
       end if

     end function buffermap

  ! ===========================================
  !  FreeGhostBuffer:
  !  Author: Christoph Erath, Mark Taylor
  !  Freed an ghostpoints communication buffer
  ! =========================================
  subroutine FreeGhostBufferTR(ghost) 
    use edgetype_mod, only : GhostBuffertr_t 
    implicit none
    type (GhostBuffertr_t),intent(inout) :: ghost

#if (defined HORIZ_OPENMP)
!$OMP BARRIER
!$OMP MASTER
#endif
    ghost%nbuf=0
    ghost%nlyr=0
    deallocate(ghost%buf)
    deallocate(ghost%receive)
#if (defined HORIZ_OPENMP)
!$OMP END MASTER
#endif

  end subroutine FreeGhostBufferTR 
  ! =========================================

  ! =========================================
  !
  !> @brief Pack edges of v into an edge buffer for boundary exchange.
  !
  !> This subroutine packs for one or more vertical layers into an edge 
  !! buffer. If the buffer associated with edge is not large enough to 
  !! hold all vertical layers you intent to pack, the method will 
  !! halt the program with a call to parallel_mod::haltmp().
  !! @param[in] edge Ghost Buffer into which the data will be packed.
  !! This buffer must be previously allocated with initghostbufferfull().
  !! @param[in] v The data to be packed.
  !! @param[in] vlyr Number of vertical level coming into the subroutine
  !! for packing for input v.
  !! @param[in] kptr Vertical pointer to the place in the edge buffer where 
  !! data will be located.
  ! =========================================
  subroutine GhostVpackfull(edge,v,nc1,nc2,nc,vlyr,kptr,desc)

    use dimensions_mod, only : max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use edgetype_mod, only : ghostbuffer3d_t

    implicit none 
    type (Ghostbuffer3D_t)                      :: edge
    integer,              intent(in)   :: vlyr
    integer,              intent(in)   :: nc1,nc2,nc,kptr
    real (kind=real_kind),intent(in)   :: v(nc1:nc2,nc1:nc2,vlyr)
    type (EdgeDescriptor_t),intent(in) :: desc

    ! Local variables
    integer :: i,k,ir,l,e

    integer :: is,ie,in,iw

    if(.not. threadsafe) then
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
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
    is = desc%putmapP_ghost(south) 
    ie = desc%putmapP_ghost(east)  
    in = desc%putmapP_ghost(north) 
    iw = desc%putmapP_ghost(west)  
#if 1
    if (is>edge%nbuf) call abortmp('error is=')
    if (ie>edge%nbuf) call abortmp('error ie=')
    if (in>edge%nbuf) call abortmp('error in=')
    if (iw>edge%nbuf) call abortmp('error iw=')
    if (is<1) call abortmp('error is=0')
    if (ie<1) call abortmp('error ie=0')
    if (in<1) call abortmp('error in=0')
    if (iw<1) call abortmp('error iw=0')
#endif

!    print *,nc,is,ie,in,iw
    do k=1,vlyr
       do e=1,nc
       do i=1,nc
          edge%buf(i,e,kptr+k,is)   = v(i  ,e ,k)
          edge%buf(i,e,kptr+k,ie)   = v(nc-e+1 ,i ,k)
          edge%buf(i,e,kptr+k,in)   = v(i  ,nc-e+1,k)
          edge%buf(i,e,kptr+k,iw)   = v(e  ,i ,k)
       enddo
       end do
    end do

    !  This is really kludgy way to setup the index reversals
    !  But since it is so a rare event not real need to spend time optimizing
    !  Check if the edge orientation of the recieving element is different
    !  if it is, swap the order of data in the edge
    if(desc%reverse(south)) then
       is = desc%putmapP_ghost(south)
       do e=1,nc
       do k=1,vlyr
          do i=1,nc
             ir = nc-i+1
             edge%buf(ir,e,kptr+k,is)=v(i,e,k)
          enddo
       enddo
       enddo
    endif

    if(desc%reverse(east)) then
       ie = desc%putmapP_ghost(east)
       do e=1,nc
       do k=1,vlyr
          do i=1,nc
             ir = nc-i+1
             edge%buf(ir,e,kptr+k,ie)=v(nc-e+1,i,k)
          enddo
       enddo
       enddo
    endif

    if(desc%reverse(north)) then
       in = desc%putmapP_ghost(north)
       do e=1,nc
       do k=1,vlyr
          do i=1,nc
             ir = nc-i+1
             edge%buf(ir,e,kptr+k,in)=v(i,nc-e+1,k)
          enddo
       enddo
       enddo
    endif

    if(desc%reverse(west)) then
       iw = desc%putmapP_ghost(west)
       do e=1,nc
       do k=1,vlyr
          do i=1,nc
             ir = nc-i+1
             edge%buf(ir,e,kptr+k,iw)=v(e,i,k)
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
             do e=1,nc
                edge%buf(:,e,kptr+k,desc%putmapP_ghost(l))=v(1:nc  ,e ,k)
             enddo
          end do
       end if
    end do
    
! SEAST
    do l=swest+max_corner_elem, swest+2*max_corner_elem-1
       if (desc%putmapP_ghost(l) /= -1) then
          do k=1,vlyr
             do e=1,nc
                edge%buf(e,:,kptr+k,desc%putmapP_ghost(l))=v(nc-e+1 ,1:nc ,k)
             enddo
          end do
       end if
    end do
    
! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
       if (desc%putmapP_ghost(l) /= -1) then
          do k=1,vlyr
             do e=1,nc
                do i=1,nc
                   edge%buf(i,e,kptr+k,desc%putmapP_ghost(l))=v(nc-i+1,nc-e+1,k)
                enddo
             enddo
             end do
          end if
    end do
    
! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
       if (desc%putmapP_ghost(l) /= -1) then
          do k=1,vlyr
             do e=1,nc
                edge%buf(:,e,kptr+k,desc%putmapP_ghost(l))=v(1:nc  ,nc-e+1,k)
             enddo
          end do
       end if
    end do
  end subroutine GhostVpackfull

  ! ========================================
  ! edgeVunpack:
  !
  ! Unpack edges from edge buffer into v...
  ! ========================================

  subroutine GhostVunpackfull(edge,v,nc1,nc2,nc,vlyr,kptr,desc)
    use dimensions_mod, only : max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use edgetype_mod, only : Ghostbuffer3d_t
    implicit none 
    type (Ghostbuffer3D_t),         intent(in)  :: edge

    integer,               intent(in)  :: vlyr
    integer,               intent(in)  :: kptr,nc1,nc2,nc
    real (kind=real_kind), intent(inout) :: v(nc1:nc2,nc1:nc2,vlyr)
    type (EdgeDescriptor_t)            :: desc

    ! Local
    logical, parameter :: UseUnroll = .TRUE.
    integer :: i,k,l,e
    integer :: is,ie,in,iw,ic
    logical :: reverse

    ! make sure buffer is big enough:
    if ( (nc2-nc1+1) <  3*nc ) then
       call haltmp("GhostVunpack:  insufficient ghost cell region")
    endif


    threadsafe=.false.

    is=desc%getmapP_ghost(south) 
    ie=desc%getmapP_ghost(east)  
    in=desc%getmapP_ghost(north) 
    iw=desc%getmapP_ghost(west)  

    ! example for north buffer
    ! first row ('edge') goes in v(:,np+1,k)
    ! 2nd   row ('edge') goes in v(:,np+2,k)
    ! etc...
    do e=1,nc
    do k=1,vlyr
       do i=1,nc
          v(i  ,1-e  ,k) = edge%buf(i,e,kptr+k,is  )
          v(nc+e ,i  ,k) = edge%buf(i,e,kptr+k,ie  )
          v(i  ,nc+e ,k) = edge%buf(i,e,kptr+k,in  )
          v(1-e  ,i  ,k) = edge%buf(i,e,kptr+k,iw  )
       end do
    end do
    end do

! SWEST
    do l=swest,swest+max_corner_elem-1
       ic = desc%getmapP_ghost(l)
       if (ic /= -1 .and. l.eq.swest) then   ! unpack swest corner, IGNORE all other corners
       !if(ic /= -1) then 
          reverse=desc%reverse(l)
          if (reverse) then
             do k=1,vlyr
                !v(0  ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
                do e=1,nc
                   do i=1,nc
                      v(1-e,1-i,k)=edge%buf(i,e,kptr+k,ic)
                   enddo
                enddo
             enddo
          else
             do k=1,vlyr
                !v(0  ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
                do e=1,nc
                   do i=1,nc
                      v(1-e,1-i,k)=edge%buf(e,i,kptr+k,ic)
                   enddo
                enddo
             enddo
          endif
       endif
    end do
    
! SEAST
    do l=swest+max_corner_elem,swest+2*max_corner_elem-1
       ic = desc%getmapP_ghost(l)
       if (ic /= -1 .and. l.eq.seast) then   ! unpack seast corner, IGNORE all other corners
       !if(ic /= -1) then 
          reverse=desc%reverse(l)
          if (reverse) then
             do k=1,vlyr
                !v(nc+1 ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
                do e=1,nc
                   do i=1,nc
                      v(nc+i,1-e,k)=edge%buf(e,i,kptr+k,ic)
                   enddo
                enddo
             enddo
          else
             do k=1,vlyr
                !v(nc+1 ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
                do e=1,nc
                   do i=1,nc
                      v(nc+i ,1-e ,k)=edge%buf(i,e,kptr+k,ic)
                   enddo
                enddo
             enddo
          endif
       endif
    end do

! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
       ic = desc%getmapP_ghost(l)
       if (ic /= -1 .and. l.eq.neast) then   ! unpack neast corner, IGNORE all other corners
       !if(ic /= -1) then 
          reverse=desc%reverse(l)
          if (reverse) then
             do k=1,vlyr
                !v(nc+1 ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
                do e=1,nc
                   do i=1,nc
                      v(nc+i ,nc+e,k)=edge%buf(e,i,kptr+k,ic)
                   enddo
                enddo
             enddo
          else
             do k=1,vlyr
                !v(nc+1 ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
                do e=1,nc
                   do i=1,nc
                      v(nc+i ,nc+e,k)=edge%buf(i,e,kptr+k,ic)
                   enddo
                enddo
             enddo
          endif
       endif
    end do

! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
       ic = desc%getmapP_ghost(l)
       if (ic /= -1 .and. l.eq.nwest) then
          ! unpack nwest corner, IGNORE all other corners
          reverse=desc%reverse(l)
          if (reverse) then
             do k=1,vlyr
                !v(0  ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
                do e=1,nc
                   do i=1,nc
                      v(1-i ,nc+e,k)=edge%buf(e,i,kptr+k,ic)
                   enddo
                enddo
             enddo
          else
             do k=1,vlyr
                !v(0  ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
                do e=1,nc
                   do i=1,nc
                      v(1-i ,nc+e,k)=edge%buf(i,e,kptr+k,ic)
                   enddo
                enddo
             enddo
          endif
       endif
    end do

  end subroutine GhostVunpackfull






  ! =========================================
  !
  !> @brief Pack edges of v into an edge buffer for boundary exchange.
  !
  !> This subroutine packs for one or more vertical layers into an edge 
  !! buffer. If the buffer associated with edge is not large enough to 
  !! hold all vertical layers you intent to pack, the method will 
  !! halt the program with a call to parallel_mod::haltmp().
  !! @param[in] edge Ghost Buffer into which the data will be packed.
  !! This buffer must be previously allocated with initghostbuffer().
  !! @param[in] v The data to be packed.
  !! @param[in] vlyr Number of vertical level coming into the subroutine
  !! for packing for input v.
  !! @param[in] kptr Vertical pointer to the place in the edge buffer where 
  !! data will be located.
  ! =========================================
  subroutine GhostVpack_unoriented(edge,v,nc,vlyr,kptr,desc)
    use edgetype_mod, only : edgedescriptor_t, ghostbuffer3d_t 
    implicit none
    type (Ghostbuffer3D_t),intent(inout) :: edge
    integer,              intent(in)   :: vlyr
    integer,              intent(in)   :: nc
    real (kind=real_kind),intent(in)   :: v(nc,nc,vlyr)
    integer,              intent(in)   :: kptr
    type (EdgeDescriptor_t),intent(in) :: desc

    integer :: k,l,l_local,is

    if(.not. threadsafe) then
#if ( defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
       threadsafe=.true.
    end if

    do l_local=1,desc%actual_neigh_edges
       l=desc%loc2buf(l_local)
       is = desc%putmapP_ghost(l)
       do k=1,vlyr
          edge%buf(:,:,kptr+k,is) = v(:,:,k)  
       enddo
    end do

  end subroutine GhostVpack_unoriented


  ! ========================================
  ! edgeVunpack:
  !
  ! Unpack edges from edge buffer into v...
  ! ========================================
  subroutine GhostVunpack_unoriented(edge,v,nc,vlyr,kptr,desc,GlobalId,u)

    use edgetype_mod, only : Ghostbuffer3d_t, EdgeDescriptor_t
    implicit none

    type (Ghostbuffer3D_t),intent(inout)  :: edge
    integer,               intent(in)     :: vlyr
    integer,               intent(in)     :: nc
    real (kind=real_kind), intent(out)    :: v(nc,nc,vlyr,*)
    integer,               intent(in)     :: kptr
    type (EdgeDescriptor_t),intent(in)    :: desc
    integer(kind=int_kind),intent(in)     :: GlobalId
    real (kind=real_kind), intent(in)     :: u(nc,nc,vlyr)

    integer :: k,l,n,is,m,pid,gid

    threadsafe=.false.

    m=0
    gid = GlobalID
    do n=1,desc%actual_neigh_edges+1
       pid = desc%globalID(desc%loc2buf((m+1)))
       if (m==desc%actual_neigh_edges .OR. pid < gid) then
          gid = -1
          v(:,:,:,n) = u(:,:,:)
       else
          m = m+1
          l = desc%loc2buf(m)
          is = desc%getmapP_ghost(l)
          do k=1,vlyr
             v(:,:,k,n) = edge%buf(:,:,kptr+k,is) 
          enddo
       end if
    end do


  end subroutine GhostVunpack_unoriented











  ! =========================================
  ! initGhostBuffer:
  ! Author: Christoph Erath
  ! create an Real based communication buffer
  ! npoints is the number of points on one side
  ! nhc is the deep of the ghost/halo zone
  ! =========================================
  subroutine initGhostBufferTR(ghost,nlyr,ntrac,nhc, npoints)

    implicit none
    integer,intent(in)                :: nlyr,ntrac,nhc, npoints
    type (Ghostbuffertr_t),intent(out)  :: ghost

    ! Local variables

    integer :: nbuf

    ! make sure buffer is big enough:
    if (nhc>npoints ) then
       call haltmp("intGhostBuffer:  halo region can not be larger then element size")
    endif
    if (ntrac<1 ) then
       call haltmp("intGhostBuffer:  you have to consider at least one tracer")
    endif

    ! sanity check for threading
    if (omp_get_num_threads()>1) then
       call haltmp('ERROR: initGhostBuffer must be called before threaded region')
    endif

    nbuf=max_neigh_edges*nelemd
#if (defined HORIZ_OPENMP)
!$OMP MASTER
#endif
    ghost%nlyr=nlyr
    ghost%nbuf=nbuf
    allocate(ghost%buf(npoints,nhc,nlyr,ntrac,nbuf))
    ghost%buf=0

    allocate(ghost%receive(npoints,nhc,nlyr,ntrac,nbuf))
    ghost%receive=0
#if (defined HORIZ_OPENMP)
!$OMP END MASTER
!$OMP BARRIER
#endif

  end subroutine initGhostBufferTR



! =========================================
! Christoph Erath
!> Packs the halo zone from v
! =========================================

subroutine ghostVpack(edge,v,nhc,npoints,vlyr,ntrac,kptr,desc)
  
  use dimensions_mod, only : max_corner_elem
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  use edgetype_mod, only : EdgeDescriptor_t, Ghostbuffertr_t

  implicit none

  type (Ghostbuffertr_t)                      :: edge
  integer,              intent(in)   :: vlyr
  integer,              intent(in)   :: ntrac
  integer,              intent(in)   :: npoints,nhc, kptr
  
  real (kind=real_kind),intent(in)   :: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc,vlyr,ntrac)
  type (EdgeDescriptor_t),intent(in) :: desc

  ! Local variables
  integer :: i,j,k,ir,l,itr

  integer :: is,ie,in,iw

  if(.not. threadsafe) then
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
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
  is = desc%putmapP_ghost(south) 
  ie = desc%putmapP_ghost(east)  
  in = desc%putmapP_ghost(north) 
  iw = desc%putmapP_ghost(west)  

!    print *,nc,is,ie,in,iw
  do itr=1,ntrac
    do k=1,vlyr
      do i=1,npoints
        do j=1,nhc
          edge%buf(i,j,k,itr+kptr,is)   = v(i  ,j ,k,itr)
          edge%buf(i,j,k,itr+kptr,ie)   = v(npoints-j+1 ,i ,k,itr)
          edge%buf(i,j,k,itr+kptr,in)   = v(i  ,npoints-j+1,k,itr)
          edge%buf(i,j,k,itr+kptr,iw)   = v(j  ,i ,k,itr)
        enddo
      end do
    end do
  end do

  !  This is really kludgy way to setup the index reversals
  !  But since it is so a rare event not real need to spend time optimizing
  !  Check if the edge orientation of the recieving element is different
  !  if it is, swap the order of data in the edge
  if(desc%reverse(south)) then
!      is = desc%putmapP_ghost(south)
    do itr=1,ntrac
     do k=1,vlyr
       do i=1,npoints
         do j=1,nhc
           ir = npoints-i+1
           edge%buf(ir,j,k,itr+kptr,is)=v(i,j,k,itr)
         enddo
       enddo
     enddo
   enddo
  endif

  if(desc%reverse(east)) then
!      ie = desc%putmapP_ghost(east)
   do itr=1,ntrac
     do k=1,vlyr
       do i=1,npoints
         do j=1,nhc
           ir = npoints-i+1
           edge%buf(ir,j,k,itr+kptr,ie)=v(npoints-j+1,i,k,itr)
          enddo
        enddo
      enddo
    end do
  endif

  if(desc%reverse(north)) then
!      in = desc%putmapP_ghost(north)
    do itr=1,ntrac
      do k=1,vlyr
        do i=1,npoints
          do j=1,nhc
            ir = npoints-i+1
            edge%buf(ir,j,k,itr+kptr,in)=v(i,npoints-j+1,k,itr)
          enddo
        enddo
      enddo
    enddo
  endif

  if(desc%reverse(west)) then
!      iw = desc%putmapP_ghost(west)
   do itr=1,ntrac
     do k=1,vlyr
       do i=1,npoints
         do j=1,nhc
            ir = npoints-i+1
            edge%buf(ir,j,k,itr+kptr,iw)=v(j,i,k,itr)
          enddo
        enddo
      enddo
    enddo
  endif



  ! corners.  this is difficult because we dont know the orientaton
  ! of the corners, and this which (i,j) dimension maps to which dimension
! SWEST
  do l=swest,swest+max_corner_elem-1
     if (l.ne.swest) call abortmp('ERROR2: swest ghost cell update requires ne>0 cubed-sphere mesh')
     if (desc%putmapP_ghost(l) /= -1) then
       do itr=1,ntrac
         do k=1,vlyr
           ! edge%buf(1,1,kptr+k,desc%putmapP_ghost(l))=v(1,1 ,k)
          do i=1,nhc
            do j=1,nhc
              edge%buf(i,j,k,itr+kptr,desc%putmapP_ghost(l))=v(i  ,j ,k,itr)
            enddo
          end do
        end do
      enddo
     end if
  end do
  
! SEAST
  do l=swest+max_corner_elem,swest+2*max_corner_elem-1
     if (l.ne.seast) call abortmp('ERROR2: seast ghost cell update requires ne>0 cubed-sphere mesh')
     if (desc%putmapP_ghost(l) /= -1) then
       do itr=1,ntrac
         do k=1,vlyr
           ! edge%buf(1,1,kptr+k,desc%putmapP_ghost(l))=v(nc ,1 ,k)
          do i=1,nhc
            do j=1,nhc
              edge%buf(i,j,k,itr+kptr,desc%putmapP_ghost(l))=v(npoints-i+1 ,j ,k,itr)
            enddo
          end do
        end do
      enddo
     end if
  end do
  
! NEAST
  do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
     if (l.ne.neast) call abortmp('ERROR2: neast ghost cell update requires ne>0 cubed-sphere mesh')
     if (desc%putmapP_ghost(l) /= -1) then
       do itr=1,ntrac
        do k=1,vlyr
           !edge%buf(1,1,kptr+k,desc%putmapP_ghost(l))=v(nc ,nc,k)
           do i=1,nhc
              do j=1,nhc
                 edge%buf(i,j,k,itr+kptr,desc%putmapP_ghost(l))=v(npoints-i+1,npoints-j+1,k,itr)
              enddo
            enddo
          end do
        enddo
      end if
  end do
  
! NWEST
  do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
     if (l.ne.nwest) call abortmp('ERROR2: nwest ghost cell update requires ne>0 cubed-sphere mesh')
     if (desc%putmapP_ghost(l) /= -1) then
      do itr=1,ntrac
        do k=1,vlyr
           !edge%buf(1,1,kptr+k,desc%putmapP_ghost(l))=v(1  ,nc,k)
          do i=1,nhc
            do j=1,nhc
              edge%buf(i,j,k,itr+kptr,desc%putmapP_ghost(l))=v(i  ,npoints-j+1,k,itr)
            enddo
          end do
        end do
       enddo
     end if
  end do
end subroutine GhostVpack


! =========================================
! Christoph Erath
!> Packs the halo zone from v
! =========================================
! NOTE: I have to give timelevels as argument, because element_mod is not compiled first
! and the array call has to be done in this way because of performance reasons!!!
subroutine ghostVpackR(edge,v,nhc,npoints,vlyr,ntrac,kptr,desc)
  use dimensions_mod, only : max_corner_elem!, ntrac_d
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  use edgetype_mod, only : Ghostbuffertr_t
  implicit none 
  type (Ghostbuffertr_t)                      :: edge
  integer,              intent(in)   :: vlyr
  integer,              intent(in)   :: ntrac
  integer,              intent(in)   :: npoints,nhc,kptr
  
!  real (kind=real_kind),intent(in)   :: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc,vlyr,ntrac_d)
  real (kind=real_kind),intent(in)   :: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc,vlyr,ntrac)!phl
  type (EdgeDescriptor_t),intent(in) :: desc

  ! Local variables
  integer :: i,j,k,ir,l,itr

  integer :: is,ie,in,iw

  if(.not. threadsafe) then
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
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
  is = desc%putmapP_ghost(south) 
  ie = desc%putmapP_ghost(east)  
  in = desc%putmapP_ghost(north) 
  iw = desc%putmapP_ghost(west)  

!    print *,nc,is,ie,in,iw
  do itr=1,ntrac
    do k=1,vlyr
      do i=1,npoints
        do j=1,nhc
          edge%buf(i,j,kptr+k,itr,is)   = v(i  ,j ,k,itr)
          edge%buf(i,j,kptr+k,itr,ie)   = v(npoints-j+1 ,i ,k,itr)
          edge%buf(i,j,kptr+k,itr,in)   = v(i  ,npoints-j+1,k,itr)
          edge%buf(i,j,kptr+k,itr,iw)   = v(j  ,i ,k,itr)
        enddo
      end do
    end do
  end do

  !  This is really kludgy way to setup the index reversals
  !  But since it is so a rare event not real need to spend time optimizing
  !  Check if the edge orientation of the recieving element is different
  !  if it is, swap the order of data in the edge
  if(desc%reverse(south)) then
!      is = desc%putmapP_ghost(south)
    do itr=1,ntrac
     do k=1,vlyr
       do i=1,npoints
         do j=1,nhc
           ir = npoints-i+1
           edge%buf(ir,j,kptr+k,itr,is)=v(i,j,k,itr)
         enddo
       enddo
     enddo
   enddo
  endif

  if(desc%reverse(east)) then
!      ie = desc%putmapP_ghost(east)
   do itr=1,ntrac
     do k=1,vlyr
       do i=1,npoints
         do j=1,nhc
           ir = npoints-i+1
           edge%buf(ir,j,kptr+k,itr,ie)=v(npoints-j+1,i,k,itr)
          enddo
        enddo
      enddo
    end do
  endif

  if(desc%reverse(north)) then
!      in = desc%putmapP_ghost(north)
    do itr=1,ntrac
      do k=1,vlyr
        do i=1,npoints
          do j=1,nhc
            ir = npoints-i+1
            edge%buf(ir,j,kptr+k,itr,in)=v(i,npoints-j+1,k,itr)
          enddo
        enddo
      enddo
    enddo
  endif

  if(desc%reverse(west)) then
!      iw = desc%putmapP_ghost(west)
   do itr=1,ntrac
     do k=1,vlyr
       do i=1,npoints
         do j=1,nhc
            ir = npoints-i+1
            edge%buf(ir,j,kptr+k,itr,iw)=v(j,i,k,itr)
          enddo
        enddo
      enddo
    enddo
  endif



  ! corners.  this is difficult because we dont know the orientaton
  ! of the corners, and this which (i,j) dimension maps to which dimension
! SWEST
  do l=swest,swest+max_corner_elem-1
     if (l.ne.swest) call abortmp('ERROR2: swest ghost cell update requires ne>0 cubed-sphere mesh')
     if (desc%putmapP_ghost(l) /= -1) then
       do itr=1,ntrac
         do k=1,vlyr
           ! edge%buf(1,1,kptr+k,desc%putmapP_ghost(l))=v(1,1 ,k)
          do i=1,nhc
            do j=1,nhc
              edge%buf(i,j,kptr+k,itr,desc%putmapP_ghost(l))=v(i  ,j ,k,itr)
            enddo
          end do
        end do
      enddo
     end if
  end do
  
! SEAST
  do l=swest+max_corner_elem,swest+2*max_corner_elem-1
     if (l.ne.seast) call abortmp('ERROR2: seast ghost cell update requires ne>0 cubed-sphere mesh')
     if (desc%putmapP_ghost(l) /= -1) then
       do itr=1,ntrac
         do k=1,vlyr
           ! edge%buf(1,1,kptr+k,desc%putmapP_ghost(l))=v(nc ,1 ,k)
          do i=1,nhc
            do j=1,nhc
              edge%buf(i,j,kptr+k,itr,desc%putmapP_ghost(l))=v(npoints-i+1 ,j ,k,itr)
            enddo
          end do
        end do
      enddo
     end if
  end do
  
! NEAST
  do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
     if (l.ne.neast) call abortmp('ERROR2: neast ghost cell update requires ne>0 cubed-sphere mesh')
     if (desc%putmapP_ghost(l) /= -1) then
       do itr=1,ntrac
        do k=1,vlyr
           !edge%buf(1,1,kptr+k,desc%putmapP_ghost(l))=v(nc ,nc,k)
           do i=1,nhc
              do j=1,nhc
                 edge%buf(i,j,kptr+k,itr,desc%putmapP_ghost(l))=v(npoints-i+1,npoints-j+1,k,itr)
              enddo
            enddo
          end do
        enddo
      end if
  end do
  
! NWEST
  do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
     if (l.ne.nwest) call abortmp('ERROR2: nwest ghost cell update requires ne>0 cubed-sphere mesh')
     if (desc%putmapP_ghost(l) /= -1) then
      do itr=1,ntrac
        do k=1,vlyr
           !edge%buf(1,1,kptr+k,desc%putmapP_ghost(l))=v(1  ,nc,k)
          do i=1,nhc
            do j=1,nhc
              edge%buf(i,j,kptr+k,itr,desc%putmapP_ghost(l))=v(i  ,npoints-j+1,k,itr)
            enddo
          end do
        end do
       enddo
     end if
  end do
end subroutine GhostVpackR
! ========================================
! Christoph Erath
!
! Unpack the halo zone into v
! ========================================
subroutine ghostVunpack(edge,v,nhc,npoints,vlyr,ntrac,kptr,desc)
  use dimensions_mod, only : max_corner_elem
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  type (Ghostbuffertr_t),         intent(in)  :: edge

  integer,               intent(in)  :: vlyr
  integer,              intent(in)   :: ntrac
  integer,               intent(in)  :: kptr,nhc,npoints
  
  real (kind=real_kind), intent(inout) :: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc,vlyr,ntrac)
  
  type (EdgeDescriptor_t)            :: desc

  ! Local
  logical, parameter :: UseUnroll = .TRUE.
  integer :: i,j,k,l,itr
  integer :: is,ie,in,iw,ic
  logical :: reverse
  
  threadsafe=.false.

  is=desc%getmapP_ghost(south) 
  ie=desc%getmapP_ghost(east)  
  in=desc%getmapP_ghost(north) 
  iw=desc%getmapP_ghost(west)  

  ! example for north buffer
  ! first row ('edge') goes in v(:,np+1,k)
  ! 2nd   row ('edge') goes in v(:,np+2,k)
  ! etc...
  do itr=1,ntrac
    do k=1,vlyr
      do i=1,npoints
       do j=1,nhc
          v(i  ,1-j  ,k,itr)      = edge%buf(i,j,k,itr+kptr,is  )
          v(npoints+j ,i  ,k,itr) = edge%buf(i,j,k,itr+kptr,ie  )
          v(i  ,npoints+j ,k,itr) = edge%buf(i,j,k,itr+kptr,in  )
          v(1-j  ,i  ,k,itr)      = edge%buf(i,j,k,itr+kptr,iw  )
       end do
      end do
    end do
  enddo

! SWEST
  do l=swest,swest+max_corner_elem-1
     ic = desc%getmapP_ghost(l)
     if(ic /= -1) then 
        reverse=desc%reverse(l)
        if (reverse) then
          do itr=1,ntrac
           do k=1,vlyr
              !v(0  ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
              do j=1,nhc
                 do i=1,nhc
                    v(1-j,1-i,k,itr)=edge%buf(i,j,k,itr+kptr,ic)
                 enddo
              enddo
             enddo
           enddo
        else
          do itr=1,ntrac
           do k=1,vlyr
              !v(0  ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
              do j=1,nhc
                 do i=1,nhc
                    v(1-j,1-i,k,itr)=edge%buf(j,i,k,itr+kptr,ic)
                 enddo
              enddo
             enddo
           enddo
        endif
     else
       do itr=1,ntrac
         do k=1,vlyr
           do j=1,nhc
             do i=1,nhc
               v(1-i,1-j,k,itr)=edgeDefaultVal
             enddo
           enddo
         enddo
       end do    
     endif
  end do
  
! SEAST
  do l=swest+max_corner_elem,swest+2*max_corner_elem-1
     ic = desc%getmapP_ghost(l)
     if(ic /= -1) then 
        reverse=desc%reverse(l)
        if (reverse) then
          do itr=1,ntrac
           do k=1,vlyr
              !v(nc+1 ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
              do j=1,nhc
                 do i=1,nhc
                    v(npoints+i,1-j,k,itr)=edge%buf(j,i,k,itr+kptr,ic)
                 enddo
              enddo
            enddo
           enddo
        else
          do itr=1,ntrac
           do k=1,vlyr
              !v(nc+1 ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
              do j=1,nhc
                 do i=1,nhc
                    v(npoints+i ,1-j ,k,itr)=edge%buf(i,j,k,itr+kptr,ic)
                 enddo
              enddo
            enddo
           enddo
        endif
      else
        do itr=1,ntrac
         do k=1,vlyr
          do j=1,nhc
            do i=1,nhc
              v(npoints+i,1-j,k,itr)=edgeDefaultVal
            enddo
          enddo
         enddo
        end do  
     endif
  end do

! NEAST
  do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
     ic = desc%getmapP_ghost(l)
     if(ic /= -1) then 
        reverse=desc%reverse(l)
        if (reverse) then
          do itr=1,ntrac
           do k=1,vlyr
              !v(nc+1 ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
              do j=1,nhc
                 do i=1,nhc
                    v(npoints+i ,npoints+j,k,itr)=edge%buf(j,i,k,itr+kptr,ic)
                 enddo
              enddo
            enddo
           enddo
        else
          do itr=1,ntrac
           do k=1,vlyr
              !v(nc+1 ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
              do j=1,nhc
                 do i=1,nhc
                    v(npoints+i ,npoints+j,k,itr)=edge%buf(i,j,k,itr+kptr,ic)
                 enddo
              enddo
            enddo
           enddo
        endif
      else
        do itr=1,ntrac
         do k=1,vlyr        
          do j=1,nhc
            do i=1,nhc
              v(npoints+i,npoints+j,k,itr)=edgeDefaultVal
            enddo
          enddo
         enddo
        end do    
     endif
  end do

! NWEST
  do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
     ic = desc%getmapP_ghost(l)
     if(ic /= -1) then 
        reverse=desc%reverse(l)
        if (reverse) then
          do itr=1,ntrac
           do k=1,vlyr
              !v(0  ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
              do j=1,nhc
                 do i=1,nhc
                    v(1-i ,npoints+j,k,itr)=edge%buf(j,i,k,itr+kptr,ic)
                 enddo
              enddo
            enddo
           enddo
        else
          do itr=1,ntrac
           do k=1,vlyr
              !v(0  ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
              do j=1,nhc
                 do i=1,nhc
                    v(1-i ,npoints+j,k,itr)=edge%buf(i,j,k,itr+kptr,ic)
                 enddo
              enddo
            enddo
           enddo
        endif
      else
        do itr=1,ntrac
         do k=1,vlyr
          do j=1,nhc
            do i=1,nhc
              v(1-i,npoints+j,k,itr)=edgeDefaultVal
            enddo
          enddo
         enddo
        end do    
     endif
  end do

end subroutine ghostVunpack


! ========================================
! Christoph Erath
!
! Unpack the halo zone into v
! ========================================
! NOTE: I have to give timelevels as argument, because element_mod is not compiled first
! and the array call has to be done in this way because of performance reasons!!!
subroutine ghostVunpackR(edge,v,nhc,npoints,vlyr,ntrac,kptr,desc)
  use dimensions_mod, only : max_corner_elem!, ntrac_d
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  type (Ghostbuffertr_t),         intent(in)  :: edge

  integer,               intent(in)  :: vlyr
  integer,              intent(in)   :: ntrac
  integer,               intent(in)  :: kptr,nhc,npoints
  
!  real (kind=real_kind), intent(inout) :: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc,vlyr,ntrac_d)
  real (kind=real_kind), intent(inout) :: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc,vlyr,ntrac)!phl
  
  type (EdgeDescriptor_t)            :: desc

  ! Local
  logical, parameter :: UseUnroll = .TRUE.
  integer :: i,j,k,l,itr
  integer :: is,ie,in,iw,ic
  logical :: reverse
  
  threadsafe=.false.

  is=desc%getmapP_ghost(south) 
  ie=desc%getmapP_ghost(east)  
  in=desc%getmapP_ghost(north) 
  iw=desc%getmapP_ghost(west)  

  ! example for north buffer
  ! first row ('edge') goes in v(:,np+1,k)
  ! 2nd   row ('edge') goes in v(:,np+2,k)
  ! etc...
  do itr=1,ntrac
    do k=1,vlyr
      do i=1,npoints
       do j=1,nhc
          v(i  ,1-j  ,k,itr)      = edge%buf(i,j,kptr+k,itr,is  )
          v(npoints+j ,i  ,k,itr) = edge%buf(i,j,kptr+k,itr,ie  )
          v(i  ,npoints+j ,k,itr) = edge%buf(i,j,kptr+k,itr,in  )
          v(1-j  ,i  ,k,itr)      = edge%buf(i,j,kptr+k,itr,iw  )
       end do
      end do
    end do
  enddo

! SWEST
  do l=swest,swest+max_corner_elem-1
     ic = desc%getmapP_ghost(l)
     if(ic /= -1) then 
        reverse=desc%reverse(l)
        if (reverse) then
          do itr=1,ntrac
           do k=1,vlyr
              !v(0  ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
              do j=1,nhc
                 do i=1,nhc
                    v(1-j,1-i,k,itr)=edge%buf(i,j,kptr+k,itr,ic)
                 enddo
              enddo
             enddo
           enddo
        else
          do itr=1,ntrac
           do k=1,vlyr
              !v(0  ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
              do j=1,nhc
                 do i=1,nhc
                    v(1-j,1-i,k,itr)=edge%buf(j,i,kptr+k,itr,ic)
                 enddo
              enddo
             enddo
           enddo
        endif
     else
       do itr=1,ntrac
         do k=1,vlyr
           do j=1,nhc
             do i=1,nhc
               v(1-i,1-j,k,itr)=edgeDefaultVal
             enddo
           enddo
         enddo
       end do    
     endif
  end do
  
! SEAST
  do l=swest+max_corner_elem,swest+2*max_corner_elem-1
     ic = desc%getmapP_ghost(l)
     if(ic /= -1) then 
        reverse=desc%reverse(l)
        if (reverse) then
          do itr=1,ntrac
           do k=1,vlyr
              !v(nc+1 ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
              do j=1,nhc
                 do i=1,nhc
                    v(npoints+i,1-j,k,itr)=edge%buf(j,i,kptr+k,itr,ic)
                 enddo
              enddo
            enddo
           enddo
        else
          do itr=1,ntrac
           do k=1,vlyr
              !v(nc+1 ,0 ,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
              do j=1,nhc
                 do i=1,nhc
                    v(npoints+i ,1-j ,k,itr)=edge%buf(i,j,kptr+k,itr,ic)
                 enddo
              enddo
            enddo
           enddo
        endif
      else
        do itr=1,ntrac
         do k=1,vlyr
          do j=1,nhc
            do i=1,nhc
              v(npoints+i,1-j,k,itr)=edgeDefaultVal
            enddo
          enddo
         enddo
        end do  
     endif
  end do

! NEAST
  do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
     ic = desc%getmapP_ghost(l)
     if(ic /= -1) then 
        reverse=desc%reverse(l)
        if (reverse) then
          do itr=1,ntrac
           do k=1,vlyr
              !v(nc+1 ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
              do j=1,nhc
                 do i=1,nhc
                    v(npoints+i ,npoints+j,k,itr)=edge%buf(j,i,kptr+k,itr,ic)
                 enddo
              enddo
            enddo
           enddo
        else
          do itr=1,ntrac
           do k=1,vlyr
              !v(nc+1 ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
              do j=1,nhc
                 do i=1,nhc
                    v(npoints+i ,npoints+j,k,itr)=edge%buf(i,j,kptr+k,itr,ic)
                 enddo
              enddo
            enddo
           enddo
        endif
      else
        do itr=1,ntrac
         do k=1,vlyr        
          do j=1,nhc
            do i=1,nhc
              v(npoints+i,npoints+j,k,itr)=edgeDefaultVal
            enddo
          enddo
         enddo
        end do    
     endif
  end do

! NWEST
  do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
     ic = desc%getmapP_ghost(l)
     if(ic /= -1) then 
        reverse=desc%reverse(l)
        if (reverse) then
          do itr=1,ntrac
           do k=1,vlyr
              !v(0  ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
              do j=1,nhc
                 do i=1,nhc
                    v(1-i ,npoints+j,k,itr)=edge%buf(j,i,kptr+k,itr,ic)
                 enddo
              enddo
            enddo
           enddo
        else
          do itr=1,ntrac
           do k=1,vlyr
              !v(0  ,nc+1,k)=edge%buf(1,1,kptr+k,desc%getmapP_ghost(l))
              do j=1,nhc
                 do i=1,nhc
                    v(1-i ,npoints+j,k,itr)=edge%buf(i,j,kptr+k,itr,ic)
                 enddo
              enddo
            enddo
           enddo
        endif
      else
        do itr=1,ntrac
         do k=1,vlyr
          do j=1,nhc
            do i=1,nhc
              v(1-i,npoints+j,k,itr)=edgeDefaultVal
            enddo
          enddo
         enddo
        end do    
     endif
  end do

end subroutine ghostVunpackR



  ! =================================================================================
  ! GHOSTVPACK2D
  ! AUTHOR: Christoph Erath 
  ! Pack edges of v into an ghost buffer for boundary exchange.
  !
  ! This subroutine packs for one vertical layers into an ghost
  ! buffer. It is for cartesian points (v is only two dimensional). 
  ! If the buffer associated with edge is not large enough to 
  ! hold all vertical layers you intent to pack, the method will 
  ! halt the program with a call to parallel_mod::haltmp().
  ! INPUT: 
  ! - ghost Buffer into which the data will be packed.
  !   This buffer must be previously allocated with initGhostBuffer().
  ! - v The data to be packed.
  ! - nhc deep of ghost/halo zone
  ! - npoints number of points on on side
  ! - kptr Vertical pointer to the place in the edge buffer where 
  ! data will be located.
  ! =================================================================================
    subroutine ghostVpack2d(ghost,v,nhc, npoints,vlyr,ntrac,kptr,tn0,timelevels,desc)
      use dimensions_mod, only : max_corner_elem!, ntrac_d
      use control_mod, only : north, south, east, west, neast, nwest, seast, swest

      type (Ghostbuffertr_t)                  :: ghost
      integer,               intent(in)  :: vlyr
      integer,              intent(in)   :: ntrac
      integer,              intent(in)      :: nhc,npoints,kptr, tn0, timelevels
      
      real (kind=real_kind),intent(in)      :: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc, vlyr, ntrac,timelevels)!phl
!      real (kind=real_kind),intent(in)      :: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc, vlyr, ntrac_d,timelevels)
      type (EdgeDescriptor_t),intent(in)    :: desc

      ! Local variables
      integer :: i,j,ir,l,e
      integer :: itr,k
      integer :: is,ie,in,iw


      if(.not. threadsafe) then
#if (defined HORIZ_OPENMP)
  !$OMP BARRIER
#endif
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
      is = desc%putmapP_ghost(south) 
      ie = desc%putmapP_ghost(east)  
      in = desc%putmapP_ghost(north) 
      iw = desc%putmapP_ghost(west)  
      do itr=1,ntrac
       do k=1,vlyr
          do i=1,npoints
             do j=1,nhc
                ghost%buf(i,j,kptr+k,itr,is)   = v(i  ,j+1, k,itr,tn0)
                ghost%buf(i,j,kptr+k,itr,ie)   = v(npoints-j,i, k,itr,tn0 )
                ghost%buf(i,j,kptr+k,itr,in)   = v(i  ,npoints-j, k,itr,tn0)
                ghost%buf(i,j,kptr+k,itr,iw)   = v(j+1,i, k,itr,tn0 )
             enddo
          end do
        end do
      end do

      !  This is really kludgy way to setup the index reversals
      !  But since it is so a rare event not real need to spend time optimizing
      !  Check if the edge orientation of the recieving element is different
      !  if it is, swap the order of data in the edge
      if(desc%reverse(south)) then
  !        is = desc%putmapP_ghost(south)
        do itr=1,ntrac
          do k=1,vlyr
            do i=1,npoints
              do j=1,nhc
                ir = npoints-i+1
                ghost%buf(ir,j,kptr+k,itr,is)=v(i,j+1, k,itr,tn0)
               enddo
            enddo
          end do
        end do
      endif

      if(desc%reverse(east)) then
  !        ie = desc%putmapP_ghost(east)
        do itr=1,ntrac
          do k=1,vlyr  
            do i=1,npoints
              do j=1,nhc
                ir = npoints-i+1
                ghost%buf(ir,j,kptr+k,itr,ie)=v(npoints-j,i, k,itr,tn0)
              enddo
            enddo
          end do
        end do
      endif

      if(desc%reverse(north)) then
  !        in = desc%putmapP_ghost(north)
        do itr=1,ntrac
          do k=1,vlyr
            do i=1,npoints
              do j=1,nhc
                ir = npoints-i+1
                ghost%buf(ir,j,kptr+k,itr,in)=v(i,npoints-j, k,itr,tn0)
              enddo
            enddo
          end do
        end do
      endif

      if(desc%reverse(west)) then
  !        iw = desc%putmapP_ghost(west)
        do itr=1,ntrac
          do k=1,vlyr
            do i=1,npoints
              do j=1,nhc
                ir = npoints-i+1
                ghost%buf(ir,j,kptr+k,itr,iw)=v(j+1,i, k,itr,tn0)
              enddo
            enddo
          enddo
        enddo
      endif

      ! corners.  this is difficult because we dont know the orientaton
      ! of the corners, and this which (i,j) dimension maps to which dimension
  ! SWEST
      do l=swest,swest+max_corner_elem-1
         if (l.ne.swest) call abortmp('ERROR3: swest ghost cell update requires ne>0 cubed-sphere mesh')
           if (desc%putmapP_ghost(l) /= -1) then
             do itr=1,ntrac
               do k=1,vlyr
                 do i=1,nhc
                   do j=1,nhc
                     ghost%buf(i,j,kptr+k,itr,desc%putmapP_ghost(l))=v(i+1,j+1,k,itr,tn0 )
                   enddo
                 enddo     
               enddo
             enddo        
           end if
      end do

  ! SEAST
      do l=swest+max_corner_elem,swest+2*max_corner_elem-1
         if (l.ne.seast) call abortmp('ERROR3: seast ghost cell update requires ne>0 cubed-sphere mesh')
           if (desc%putmapP_ghost(l) /= -1) then
             do itr=1,ntrac
               do k=1,vlyr
                 do i=1,nhc             
                   do j=1,nhc
                     ghost%buf(i,j,kptr+k,itr,desc%putmapP_ghost(l))=v(npoints-i ,j+1,k,itr,tn0)
                   enddo
                 enddo
               enddo
             enddo             
           end if
      end do

  ! NEAST
      do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
         if (l.ne.neast) call abortmp('ERROR3: neast ghost cell update requires ne>0 cubed-sphere mesh')
           if (desc%putmapP_ghost(l) /= -1) then
             do itr=1,ntrac
               do k=1,vlyr
                 do i=1,nhc
                   do j=1,nhc
                     ghost%buf(i,j,kptr+k,itr,desc%putmapP_ghost(l))=v(npoints-i,npoints-j,k,itr,tn0)           
                   enddo
                 enddo
               enddo
             enddo
            end if
      end do

  ! NWEST
      do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
         if (l.ne.nwest) call abortmp('ERROR3: nwest ghost cell update requires ne>0 cubed-sphere mesh')
           if (desc%putmapP_ghost(l) /= -1) then
             do itr=1,ntrac
               do k=1,vlyr
                 do i=1,nhc
                   do j=1,nhc
                     ghost%buf(i,j,kptr+k,itr,desc%putmapP_ghost(l))=v(i+1,npoints-j,k,itr,tn0)
                   enddo       
                 enddo
               enddo
             enddo      
           end if
      end do   
    end subroutine ghostVpack2d

  ! =================================================================================
  ! GHOSTVUNPACK2D
  ! AUTHOR: Christoph Erath 
  ! Unpack ghost points from ghost buffer into v...
  ! It is for cartesian points (v is only two dimensional).
  ! INPUT SAME arguments as for GHOSTVPACK2d
  ! =================================================================================

  subroutine ghostVunpack2d(ghost,v,nhc,npoints,vlyr,ntrac,kptr,tn0,timelevels,desc)
    use dimensions_mod, only : max_corner_elem!, ntrac_d
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    type (Ghostbuffertr_t),         intent(in)  :: ghost

    integer,              intent(in)      :: nhc,npoints,kptr, vlyr,ntrac,tn0,timelevels

!    real (kind=real_kind), intent(inout)  :: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc, vlyr, ntrac_d,timelevels)
    real (kind=real_kind), intent(inout)  :: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc, vlyr, ntrac,timelevels)!phl
    type (EdgeDescriptor_t)               :: desc


    ! Local
    logical, parameter :: UseUnroll = .TRUE.
    integer :: i,j,l, itr, k
    integer :: is,ie,in,iw,ic
    logical :: reverse
    
    threadsafe=.false.

    is=desc%getmapP_ghost(south) 
    ie=desc%getmapP_ghost(east)  
    in=desc%getmapP_ghost(north) 
    iw=desc%getmapP_ghost(west)  

    ! example for north buffer
    ! first row ('edge') goes in v(:,np+1)
    ! 2nd   row ('edge') goes in v(:,np+2)
    ! etc...
    do itr=1,ntrac
      do k=1,vlyr
        do i=1,npoints
          do j=1,nhc
            v(i  ,1-j, k,itr,tn0  ) = ghost%buf(i,j,kptr+k,itr,is  )
            v(npoints+j ,i,k,itr,tn0 ) = ghost%buf(i,j,kptr+k,itr,ie  )
            v(i  ,npoints+j,k,itr,tn0 ) = ghost%buf(i,j,kptr+k,itr,in  )
            v(1-j  ,i,k,itr,tn0  ) = ghost%buf(i,j,kptr+k,itr,iw  )
          end do
        end do
      end do
    enddo

! SWEST
    do l=swest,swest+max_corner_elem-1
       ic = desc%getmapP_ghost(l)       
       if(ic /= -1) then 
          reverse=desc%reverse(l)
          if (reverse) then
            do itr=1,ntrac
              do k=1,vlyr
                do j=1,nhc
                   do i=1,nhc
                      v(1-i,1-j,k,itr,tn0)=ghost%buf(j,i,kptr+k,itr,ic)
                   enddo
                enddo
              enddo
            enddo
          else
            do itr=1,ntrac
              do k=1,vlyr
                do j=1,nhc
                   do i=1,nhc
                      v(1-i,1-j,k,itr,tn0)=ghost%buf(i,j,kptr+k,itr,ic)
                   enddo
                enddo
              enddo
            enddo
          endif
       else
         do itr=1,ntrac
           do k=1,vlyr
             do j=1,nhc
               do i=1,nhc
                 v(1-i,1-j,k,itr,tn0)=edgeDefaultVal
               enddo
             enddo
           enddo
         enddo      
       endif
    end do
    
! SEAST
    do l=swest+max_corner_elem,swest+2*max_corner_elem-1
       ic = desc%getmapP_ghost(l)
       if(ic /= -1) then 
          reverse=desc%reverse(l)
          if (reverse) then
            do itr=1,ntrac
              do k=1,vlyr
                do j=1,nhc
                   do i=1,nhc
                      v(npoints+i,1-j,k,itr,tn0)=ghost%buf(j,i,kptr+k,itr,ic)
                   enddo
                enddo
              enddo
            enddo
          else
            do itr=1,ntrac
              do k=1,vlyr
                do j=1,nhc
                   do i=1,nhc
                      v(npoints+i ,1-j,k,itr,tn0)=ghost%buf(i,j,kptr+k,itr,ic)
                   enddo
                enddo
              enddo
            enddo
          endif
       else
         do itr=1,ntrac
           do k=1,vlyr
              do j=1,nhc
                do i=1,nhc
                  v(npoints+i ,1-j,k,itr,tn0)=edgeDefaultVal
                enddo
              enddo
            enddo
         enddo
       endif
    end do


! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
       ic = desc%getmapP_ghost(l)
       
       if(ic /= -1) then 
          reverse=desc%reverse(l)
          if (reverse) then
            do itr=1,ntrac
              do k=1,vlyr
                do j=1,nhc
                   do i=1,nhc
                      v(npoints+i ,npoints+j,k,itr,tn0)=ghost%buf(j,i,kptr+k,itr,ic)
                   enddo
                enddo
              enddo
            enddo
          else
            do itr=1,ntrac
              do k=1,vlyr
                do j=1,nhc
                   do i=1,nhc
                      v(npoints+i ,npoints+j,k,itr,tn0)=ghost%buf(i,j,kptr+k,itr,ic)
                   enddo
                enddo
              enddo
            enddo
          endif
        else
          do itr=1,ntrac
            do k=1,vlyr
              do j=1,nhc
                do i=1,nhc
                  v(npoints+i ,npoints+j,k,itr,tn0)=edgeDefaultVal
                enddo
              enddo 
            enddo
          enddo   
       endif
    end do

! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
       ic = desc%getmapP_ghost(l)
       
       if(ic /= -1) then 
          reverse=desc%reverse(l)
          if (reverse) then
            do itr=1,ntrac
              do k=1,vlyr
                do j=1,nhc
                   do i=1,nhc
                      v(1-i ,npoints+j,k,itr,tn0)=ghost%buf(j,i,kptr+k,itr,ic)
                   enddo
                enddo
              enddo
            enddo
          else
            do itr=1,ntrac
              do k=1,vlyr
                do j=1,nhc
                   do i=1,nhc
                      v(1-i ,npoints+j,k,itr,tn0)=ghost%buf(i,j,kptr+k,itr,ic)
                   enddo
                enddo
              enddo
            enddo
          endif
        else
          do itr=1,ntrac
            do k=1,vlyr
              do j=1,nhc
                do i=1,nhc
                  v(1-i ,npoints+j,k,itr,tn0)=edgeDefaultVal
                enddo
              enddo
            enddo
          enddo
       endif
    end do

  end subroutine ghostVunpack2d



  ! =================================================================================
  ! GHOSTVPACK2D
  ! AUTHOR: Christoph Erath 
  ! Pack edges of v into an ghost buffer for boundary exchange.
  !
  ! This subroutine packs for one vertical layers into an ghost
  ! buffer. It is for cartesian points (v is only two dimensional). 
  ! If the buffer associated with edge is not large enough to 
  ! hold all vertical layers you intent to pack, the method will 
  ! halt the program with a call to parallel_mod::haltmp().
  ! INPUT: 
  ! - ghost Buffer into which the data will be packed.
  !   This buffer must be previously allocated with initGhostBuffer().
  ! - v The data to be packed.
  ! - nhc deep of ghost/halo zone
  ! - npoints number of points on on side
  ! - kptr Vertical pointer to the place in the edge buffer where 
  ! data will be located.
  ! =================================================================================
    subroutine ghostVpack2d_level(ghost,v,kptr,nhc, npoints,vlyr,desc)
      use dimensions_mod, only : max_corner_elem!, ntrac_d
      use control_mod, only : north, south, east, west, neast, nwest, seast, swest

      type (Ghostbuffertr_t)                  :: ghost
      integer,               intent(in)     :: vlyr,kptr
      integer,              intent(in)      :: nhc,npoints
      
      real (kind=real_kind),intent(in)      :: v(1-nhc:,1-nhc:,:)
      type (EdgeDescriptor_t),intent(in)    :: desc

      ! Local variables
      integer :: i,j,ir,l,e
      integer :: itr,k
      integer :: is,ie,in,iw


      if(.not. threadsafe) then
#if (defined HORIZ_OPENMP)
  !$OMP BARRIER
#endif
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
      is = desc%putmapP_ghost(south) 
      ie = desc%putmapP_ghost(east)  
      in = desc%putmapP_ghost(north) 
      iw = desc%putmapP_ghost(west)  
       do k=1,vlyr
          do i=1,npoints
             do j=1,nhc
                ghost%buf(i,j,k,kptr,is)   = v(i  ,j+1, k)
                ghost%buf(i,j,k,kptr,ie)   = v(npoints-j,i, k)
                ghost%buf(i,j,k,kptr,in)   = v(i  ,npoints-j, k)
                ghost%buf(i,j,k,kptr,iw)   = v(j+1,i, k)
             enddo
          end do
        end do

      !  This is really kludgy way to setup the index reversals
      !  But since it is so a rare event not real need to spend time optimizing
      !  Check if the edge orientation of the recieving element is different
      !  if it is, swap the order of data in the edge
      if(desc%reverse(south)) then
  !        is = desc%putmapP_ghost(south)
          do k=1,vlyr
            do i=1,npoints
              do j=1,nhc
                ir = npoints-i+1
                ghost%buf(ir,j,k,kptr,is)=v(i,j+1, k)
               enddo
            enddo
          end do
      endif

      if(desc%reverse(east)) then
  !        ie = desc%putmapP_ghost(east)
          do k=1,vlyr  
            do i=1,npoints
              do j=1,nhc
                ir = npoints-i+1
                ghost%buf(ir,j,k,kptr,ie)=v(npoints-j,i, k)
              enddo
            enddo
          end do
      endif

      if(desc%reverse(north)) then
  !        in = desc%putmapP_ghost(north)
          do k=1,vlyr
            do i=1,npoints
              do j=1,nhc
                ir = npoints-i+1
                ghost%buf(ir,j,k,kptr,in)=v(i,npoints-j, k)
              enddo
            enddo
          end do
      endif

      if(desc%reverse(west)) then
  !        iw = desc%putmapP_ghost(west)
          do k=1,vlyr
            do i=1,npoints
              do j=1,nhc
                ir = npoints-i+1
                ghost%buf(ir,j,k,kptr,iw)=v(j+1,i, k)
              enddo
            enddo
          enddo
      endif

      ! corners.  this is difficult because we dont know the orientaton
      ! of the corners, and this which (i,j) dimension maps to which dimension
  ! SWEST
      do l=swest,swest+max_corner_elem-1
         if (l.ne.swest) call abortmp('ERROR3: swest ghost cell update requires ne>0 cubed-sphere mesh')
           if (desc%putmapP_ghost(l) /= -1) then
               do k=1,vlyr
                 do i=1,nhc
                   do j=1,nhc
                     ghost%buf(i,j,k,kptr,desc%putmapP_ghost(l))=v(i+1,j+1,k)
                   enddo
                 enddo     
               enddo
           end if
      end do

  ! SEAST
      do l=swest+max_corner_elem,swest+2*max_corner_elem-1
         if (l.ne.seast) call abortmp('ERROR3: seast ghost cell update requires ne>0 cubed-sphere mesh')
           if (desc%putmapP_ghost(l) /= -1) then
               do k=1,vlyr
                 do i=1,nhc             
                   do j=1,nhc
                     ghost%buf(i,j,k,kptr,desc%putmapP_ghost(l))=v(npoints-i ,j+1,k)
                   enddo
                 enddo
               enddo
           end if
      end do

  ! NEAST
      do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
         if (l.ne.neast) call abortmp('ERROR3: neast ghost cell update requires ne>0 cubed-sphere mesh')
           if (desc%putmapP_ghost(l) /= -1) then
               do k=1,vlyr
                 do i=1,nhc
                   do j=1,nhc
                     ghost%buf(i,j,k,kptr,desc%putmapP_ghost(l))=v(npoints-i,npoints-j,k)           
                   enddo
                 enddo
               enddo
            end if
      end do

  ! NWEST
      do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
         if (l.ne.nwest) call abortmp('ERROR3: nwest ghost cell update requires ne>0 cubed-sphere mesh')
           if (desc%putmapP_ghost(l) /= -1) then
               do k=1,vlyr
                 do i=1,nhc
                   do j=1,nhc
                     ghost%buf(i,j,k,kptr,desc%putmapP_ghost(l))=v(i+1,npoints-j,k)
                   enddo       
                 enddo
               enddo
           end if
      end do   
    end subroutine ghostVpack2d_level

  ! =================================================================================
  ! GHOSTVUNPACK2D
  ! AUTHOR: Christoph Erath 
  ! Unpack ghost points from ghost buffer into v...
  ! It is for cartesian points (v is only two dimensional).
  ! INPUT SAME arguments as for GHOSTVPACK2d
  ! =================================================================================

  subroutine ghostVunpack2d_level(ghost,v,kptr,nhc,npoints,vlyr,desc)
    use dimensions_mod, only : max_corner_elem!, ntrac_d
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    type (Ghostbuffertr_t),         intent(in)  :: ghost

    integer,              intent(in)      :: nhc,npoints, vlyr,kptr

    real (kind=real_kind), intent(inout)  :: v(1-nhc:,1-nhc:,:)
    type (EdgeDescriptor_t)               :: desc


    ! Local
    logical, parameter :: UseUnroll = .TRUE.
    integer :: i,j,l, itr, k
    integer :: is,ie,in,iw,ic
    logical :: reverse
    
    threadsafe=.false.

    is=desc%getmapP_ghost(south) 
    ie=desc%getmapP_ghost(east)  
    in=desc%getmapP_ghost(north) 
    iw=desc%getmapP_ghost(west)  

    ! example for north buffer
    ! first row ('edge') goes in v(:,np+1)
    ! 2nd   row ('edge') goes in v(:,np+2)
    ! etc...
      do k=1,vlyr
        do i=1,npoints
          do j=1,nhc
            v(i  ,1-j, k) = ghost%buf(i,j,k,kptr,is  )
            v(npoints+j ,i,k ) = ghost%buf(i,j,k,kptr,ie  )
            v(i  ,npoints+j,k) = ghost%buf(i,j,k,kptr,in  )
            v(1-j  ,i,k) = ghost%buf(i,j,k,kptr,iw  )
          end do
        end do
      end do

! SWEST
    do l=swest,swest+max_corner_elem-1
       ic = desc%getmapP_ghost(l)       
       if(ic /= -1) then 
          reverse=desc%reverse(l)
          if (reverse) then
              do k=1,vlyr
                do j=1,nhc
                   do i=1,nhc
                      v(1-i,1-j,k)=ghost%buf(j,i,k,kptr,ic)
                   enddo
                enddo
              enddo
          else
              do k=1,vlyr
                do j=1,nhc
                   do i=1,nhc
                      v(1-i,1-j,k)=ghost%buf(i,j,k,kptr,ic)
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
       ic = desc%getmapP_ghost(l)
       if(ic /= -1) then 
          reverse=desc%reverse(l)
          if (reverse) then
              do k=1,vlyr
                do j=1,nhc
                   do i=1,nhc
                      v(npoints+i,1-j,k)=ghost%buf(j,i,k,kptr,ic)
                   enddo
                enddo
              enddo
          else
              do k=1,vlyr
                do j=1,nhc
                   do i=1,nhc
                      v(npoints+i ,1-j,k)=ghost%buf(i,j,k,kptr,ic)
                   enddo
                enddo
              enddo
          endif
       else
           do k=1,vlyr
              do j=1,nhc
                do i=1,nhc
                  v(npoints+i ,1-j,k)=edgeDefaultVal
                enddo
              enddo
            enddo
       endif
    end do


! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
       ic = desc%getmapP_ghost(l)
       
       if(ic /= -1) then 
          reverse=desc%reverse(l)
          if (reverse) then
              do k=1,vlyr
                do j=1,nhc
                   do i=1,nhc
                      v(npoints+i ,npoints+j,k)=ghost%buf(j,i,k,kptr,ic)
                   enddo
                enddo
              enddo
          else
              do k=1,vlyr
                do j=1,nhc
                   do i=1,nhc
                      v(npoints+i ,npoints+j,k)=ghost%buf(i,j,k,kptr,ic)
                   enddo
                enddo
              enddo
          endif
        else
            do k=1,vlyr
              do j=1,nhc
                do i=1,nhc
                  v(npoints+i ,npoints+j,k)=edgeDefaultVal
                enddo
              enddo 
            enddo
       endif
    end do

! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
       ic = desc%getmapP_ghost(l)
       
       if(ic /= -1) then 
          reverse=desc%reverse(l)
          if (reverse) then
              do k=1,vlyr
                do j=1,nhc
                   do i=1,nhc
                      v(1-i ,npoints+j,k)=ghost%buf(j,i,k,kptr,ic)
                   enddo
                enddo
              enddo
          else
              do k=1,vlyr
                do j=1,nhc
                   do i=1,nhc
                      v(1-i ,npoints+j,k)=ghost%buf(i,j,k,kptr,ic)
                   enddo
                enddo
              enddo
          endif
        else
            do k=1,vlyr
              do j=1,nhc
                do i=1,nhc
                  v(1-i ,npoints+j,k)=edgeDefaultVal
                enddo
              enddo
            enddo
       endif
    end do

  end subroutine ghostVunpack2d_level



  ! =================================================================================
  ! GHOSTVPACK2D
  ! AUTHOR: Christoph Erath 
  ! Pack edges of v into an ghost buffer for boundary exchange.
  !
  ! This subroutine packs for one vertical layers into an ghost
  ! buffer. It is for cartesian points (v is only two dimensional). 
  ! If the buffer associated with edge is not large enough to 
  ! hold all vertical layers you intent to pack, the method will 
  ! halt the program with a call to parallel_mod::haltmp().
  ! INPUT: 
  ! - ghost Buffer into which the data will be packed.
  !   This buffer must be previously allocated with initGhostBuffer().
  ! - v The data to be packed.
  ! - nhc deep of ghost/halo zone
  ! - npoints number of points on on side
  ! - kptr Vertical pointer to the place in the edge buffer where 
  ! data will be located.
  ! =================================================================================
    subroutine ghostVpack2d_single(ghost,v,nhc, npoints,desc)
      use dimensions_mod, only : max_corner_elem!, ntrac_d
      use control_mod, only : north, south, east, west, neast, nwest, seast, swest

      type (Ghostbuffertr_t)                  :: ghost
      integer,              intent(in)      :: nhc,npoints
      
      real (kind=real_kind),intent(in)      :: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc)
      type (EdgeDescriptor_t),intent(in)    :: desc

      ! Local variables
      integer :: i,j,ir,l,e
      integer :: itr,k
      integer :: is,ie,in,iw


      if(.not. threadsafe) then
#if (defined HORIZ_OPENMP)
  !$OMP BARRIER
#endif
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
      is = desc%putmapP_ghost(south) 
      ie = desc%putmapP_ghost(east)  
      in = desc%putmapP_ghost(north) 
      iw = desc%putmapP_ghost(west)  

      do i=1,npoints
         do j=1,nhc
            ghost%buf(i,j,1,1,is)   = v(i  ,j+1)
            ghost%buf(i,j,1,1,ie)   = v(npoints-j,i)
            ghost%buf(i,j,1,1,in)   = v(i  ,npoints-j)
            ghost%buf(i,j,1,1,iw)   = v(j+1,i)
         enddo
      end do

      !  This is really kludgy way to setup the index reversals
      !  But since it is so a rare event not real need to spend time optimizing
      !  Check if the edge orientation of the recieving element is different
      !  if it is, swap the order of data in the edge
      if(desc%reverse(south)) then
  !        is = desc%putmapP_ghost(south)
        do i=1,npoints
          do j=1,nhc
            ir = npoints-i+1
            ghost%buf(ir,j,1,1,is)=v(i,j+1)
           enddo
        enddo
      endif

      if(desc%reverse(east)) then
  !        ie = desc%putmapP_ghost(east) 
        do i=1,npoints
          do j=1,nhc
            ir = npoints-i+1
            ghost%buf(ir,j,1,1,ie)=v(npoints-j,i)
          enddo
        enddo
      endif

      if(desc%reverse(north)) then
  !        in = desc%putmapP_ghost(north)
        do i=1,npoints
          do j=1,nhc
            ir = npoints-i+1
            ghost%buf(ir,j,1,1,in)=v(i,npoints-j)
          enddo
        enddo
      endif

      if(desc%reverse(west)) then
  !        iw = desc%putmapP_ghost(west)
        do i=1,npoints
          do j=1,nhc
            ir = npoints-i+1
            ghost%buf(ir,j,1,1,iw)=v(j+1,i)
          enddo
        enddo
      endif

      ! corners.  this is difficult because we dont know the orientaton
      ! of the corners, and this which (i,j) dimension maps to which dimension
  ! SWEST
      do l=swest,swest+max_corner_elem-1
         if (l.ne.swest) call abortmp('ERROR3: swest ghost cell update requires ne>0 cubed-sphere mesh')
           if (desc%putmapP_ghost(l) /= -1) then
             do i=1,nhc
               do j=1,nhc
                 ghost%buf(i,j,1,1,desc%putmapP_ghost(l))=v(i+1,j+1)
               enddo
             enddo           
           end if
      end do

  ! SEAST
      do l=swest+max_corner_elem,swest+2*max_corner_elem-1
         if (l.ne.seast) call abortmp('ERROR3: seast ghost cell update requires ne>0 cubed-sphere mesh')
           if (desc%putmapP_ghost(l) /= -1) then
             do i=1,nhc             
               do j=1,nhc
                 ghost%buf(i,j,1,1,desc%putmapP_ghost(l))=v(npoints-i ,j+1)
               enddo
             enddo         
           end if
      end do

  ! NEAST
      do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
         if (l.ne.neast) call abortmp('ERROR3: neast ghost cell update requires ne>0 cubed-sphere mesh')
           if (desc%putmapP_ghost(l) /= -1) then
             do i=1,nhc
               do j=1,nhc
                 ghost%buf(i,j,1,1,desc%putmapP_ghost(l))=v(npoints-i,npoints-j)           
               enddo
             enddo
            end if
      end do

  ! NWEST
      do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
         if (l.ne.nwest) call abortmp('ERROR3: nwest ghost cell update requires ne>0 cubed-sphere mesh')
           if (desc%putmapP_ghost(l) /= -1) then
             do i=1,nhc
               do j=1,nhc
                 ghost%buf(i,j,1,1,desc%putmapP_ghost(l))=v(i+1,npoints-j)
               enddo       
             enddo     
           end if
      end do   
    end subroutine ghostVpack2d_single

  ! =================================================================================
  ! GHOSTVUNPACK2D
  ! AUTHOR: Christoph Erath 
  ! Unpack ghost points from ghost buffer into v...
  ! It is for cartesian points (v is only two dimensional).
  ! INPUT SAME arguments as for GHOSTVPACK2d
  ! =================================================================================

  subroutine ghostVunpack2d_single(ghost,v,nhc,npoints,desc)
    use dimensions_mod, only : max_corner_elem!, ntrac_d
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    type (Ghostbuffertr_t),         intent(in)  :: ghost

    integer,              intent(in)      :: nhc,npoints

    real (kind=real_kind), intent(inout)  :: v(1-nhc:npoints+nhc,1-nhc:npoints+nhc)
    type (EdgeDescriptor_t)               :: desc


    ! Local
    logical, parameter :: UseUnroll = .TRUE.
    integer :: i,j,l, itr, k
    integer :: is,ie,in,iw,ic
    logical :: reverse
    
    threadsafe=.false.

    is=desc%getmapP_ghost(south) 
    ie=desc%getmapP_ghost(east)  
    in=desc%getmapP_ghost(north) 
    iw=desc%getmapP_ghost(west)  

    ! example for north buffer
    ! first row ('edge') goes in v(:,np+1)
    ! 2nd   row ('edge') goes in v(:,np+2)
    ! etc...
    do i=1,npoints
      do j=1,nhc
        v(i  ,1-j) = ghost%buf(i,j,1,1,is  )
        v(npoints+j ,i) = ghost%buf(i,j,1,1,ie  )
        v(i  ,npoints+j) = ghost%buf(i,j,1,1,in  )
        v(1-j  ,i) = ghost%buf(i,j,1,1,iw  )
      end do
    end do

! SWEST
    do l=swest,swest+max_corner_elem-1
       ic = desc%getmapP_ghost(l)       
       if(ic /= -1) then 
          reverse=desc%reverse(l)
          if (reverse) then
            do j=1,nhc
               do i=1,nhc
                  v(1-i,1-j)=ghost%buf(j,i,1,1,ic)
               enddo
            enddo
          else
            do j=1,nhc
               do i=1,nhc
                  v(1-i,1-j)=ghost%buf(i,j,1,1,ic)
               enddo
            enddo
          endif
       else
         do j=1,nhc
           do i=1,nhc
             v(1-i,1-j)=edgeDefaultVal
           enddo
         enddo     
       endif
    end do
    
! SEAST
    do l=swest+max_corner_elem,swest+2*max_corner_elem-1
       ic = desc%getmapP_ghost(l)
       if(ic /= -1) then 
          reverse=desc%reverse(l)
          if (reverse) then
            do j=1,nhc
               do i=1,nhc
                  v(npoints+i,1-j)=ghost%buf(j,i,1,1,ic)
               enddo
            enddo
          else
            do j=1,nhc
               do i=1,nhc
                  v(npoints+i ,1-j)=ghost%buf(i,j,1,1,ic)
               enddo
            enddo
          endif
       else
          do j=1,nhc
            do i=1,nhc
              v(npoints+i ,1-j)=edgeDefaultVal
            enddo
          enddo
       endif
    end do


! NEAST
    do l=swest+3*max_corner_elem,swest+4*max_corner_elem-1
       ic = desc%getmapP_ghost(l)
       
       if(ic /= -1) then 
          reverse=desc%reverse(l)
          if (reverse) then
            do j=1,nhc
               do i=1,nhc
                  v(npoints+i ,npoints+j)=ghost%buf(j,i,1,1,ic)
               enddo
            enddo
          else
            do j=1,nhc
               do i=1,nhc
                  v(npoints+i ,npoints+j)=ghost%buf(i,j,1,1,ic)
               enddo
            enddo
          endif
        else
          do j=1,nhc
            do i=1,nhc
              v(npoints+i ,npoints+j)=edgeDefaultVal
            enddo
          enddo  
       endif
    end do

! NWEST
    do l=swest+2*max_corner_elem,swest+3*max_corner_elem-1
       ic = desc%getmapP_ghost(l)
       
       if(ic /= -1) then 
          reverse=desc%reverse(l)
          if (reverse) then
            do j=1,nhc
               do i=1,nhc
                  v(1-i ,npoints+j)=ghost%buf(j,i,1,1,ic)
               enddo
            enddo
          else
            do j=1,nhc
               do i=1,nhc
                  v(1-i ,npoints+j)=ghost%buf(i,j,1,1,ic)
               enddo
            enddo
          endif
        else
          do j=1,nhc
            do i=1,nhc
              v(1-i ,npoints+j)=edgeDefaultVal
            enddo
          enddo
       endif
    end do

  end subroutine ghostVunpack2d_single






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
       call haltmp('ERROR: initGhostBuffer must be called before threaded region')
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
  ! halt the program with a call to parallel_mod::haltmp().
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
    real (kind=real_kind),intent(in)      :: v(ghost%np, ghost%np, vlyr)
    type (EdgeDescriptor_t),intent(in)    :: desc
    
    integer                               :: nhc, np

    ! Local variables
    integer :: i,j,k,ir,l,e

    integer :: is,ie,in,iw

    if(.not. threadsafe) then
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
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
  ! AUTHOR: James Overfelt (from a subroutine of Christoph Erath, ghostVunpack2d)
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
    real (kind=real_kind), intent(inout)  :: v (1-g%nhc : g%np+g%nhc, 1-g%nhc : g%np+g%nhc, vlyr)
    integer,               intent(out)    :: mult(5:8)
    real (kind=real_kind), intent(out)    :: sw(1-g%nhc : 1,          1-g%nhc : 1,          vlyr, max_corner_elem-1)
    real (kind=real_kind), intent(out)    :: se(   g%np : g%np+g%nhc, 1-g%nhc : 1,          vlyr, max_corner_elem-1)
    real (kind=real_kind), intent(out)    :: ne(   g%np : g%np+g%nhc,    g%np : g%np+g%nhc, vlyr, max_corner_elem-1)
    real (kind=real_kind), intent(out)    :: nw(1-g%nhc : 1,             g%np : g%np+g%nhc, vlyr, max_corner_elem-1)
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


End module edge_mod_base

#if 0
#ifndef HAVE_F2003_PTR_BND_REMAP
!
! subroutine to allow sharing edge buffers
! this has to be outside a module to allow us to (F77 style) access the same chunk 
! of memory with a different shape
!
! some compilers dont allow the 'target' attribute to be used in a F77 subroutine
! such as cray.  if that is the case, try compiling with -DHAVE_F2003_PTR_BND_REMAP
!
subroutine remap_1D_ptr_buf(edge,nbuf,src_array)
use edgetype_mod, only       : EdgeBuffer_t ! _EXTERNAL                                                   
use kinds, only          : real_kind
! input                                                                                               
type (EdgeBuffer_t) :: edge
integer :: nbuf
real(kind=real_kind) , target :: src_array(nbuf)

edge%buf  => src_array

end subroutine remap_1D_ptr_buf

subroutine remap_1D_ptr_receive(edge,nbuf,src_array)
use edgetype_mod, only       : EdgeBuffer_t ! _EXTERNAL                                                   
use kinds, only          : real_kind
implicit none 
! input                                                                                               
type (EdgeBuffer_t) :: edge
integer :: nbuf
real(kind=real_kind) , target :: src_array(nbuf)

edge%receive  => src_array

end subroutine remap_1D_ptr_receive

#endif
#endif
