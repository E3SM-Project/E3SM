!This is where all of the PGI CUDA FORTRAN code will go, and these routines will be called from prim_advection_mod.
!This is compiled regardless, but PGI-specific calls are always wrapped in the USE_CUDA_FORTRAN ifdefs that are automagically
!activated when -Mcuda is specified during compilation with a PGI compiler. Thus, it will be ignored unless explicitly
!activated by the user
!
!As a general rule, all of the routines in here will be called within a threaded context (assuming ELEMENT_OPENMP is not
!deifned), and therefore, we enforce BARRIERS, MASTERS, and SINGLES from within these routines rather than outside them.
!This is to minimize the visible code impacts on the existing CPU code.

! Please pay attention to this all caps passive aggresive banner.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!                     !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!  STATUS INCOMPLETE  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!  DO NOT USE YET     !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!  UNTIL THIS BANNER  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!  IS REMOVED         !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!                     !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module cuda_mod
#if USE_CUDA_FORTRAN

! NP > 4 is not supported due to shared memory constraints
#if NP > 4
#error CUDA Fortran build only supported with NP <= 4
#endif

#define PAD 1

!Put everything CUDA-specific in here so it doesn't get compiled without -Mcuda enabled on a PGI compiler
  use cudafor
  use kinds          , only: real_kind
  use dimensions_mod , only: np,nlevp,nlev,qsize,qsize_d,max_corner_elem,max_neigh_edges,nelemd
  use element_mod    , only: timelevels
  use edge_mod       , only: EdgeBuffer_t
  implicit none
  private
  save

  !First listed are all externally accescible routines
  public :: cuda_mod_init
  public :: euler_step_cuda 
  public :: qdp_time_avg_cuda
  public :: advance_hypervis_scalar_cuda
  public :: copy_qdp_d2h
  public :: copy_qdp_h2d

  !This is from prim_advection_mod.F90
  type(EdgeBuffer_t) :: edgeAdv, edgeAdvQ3, edgeAdvQ2, edgeAdvDSS
  integer,parameter :: DSSeta = 1
  integer,parameter :: DSSomega = 2
  integer,parameter :: DSSdiv_vdp_ave = 3
  integer,parameter :: DSSno_var = -1

  !Device arrays
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:,:,:) :: qdp_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:,:)   :: qtens_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:,:)   :: qtens_biharmonic_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: spheremp_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: rspheremp_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:,:)   :: dinv_d
  real (kind=real_kind),device,allocatable,dimension(:,:)         :: deriv_dvv_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: variable_hyperviscosity_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: metdet_d
!  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: psdiss_biharmonic_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:)     :: dpdiss_biharmonic_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: rmetdet_d
  real (kind=real_kind),device,allocatable,dimension(:,:)         :: edgebuf_d
  real (kind=real_kind),device,allocatable,dimension(:)           :: hyai_d
  real (kind=real_kind),device,allocatable,dimension(:)           :: hybi_d
  logical              ,device,allocatable,dimension(:,:)         :: reverse_d
  integer              ,device,allocatable,dimension(:,:)         :: putmapP_d
  integer              ,device,allocatable,dimension(:,:)         :: getmapP_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:,:)   :: vstar_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:)     :: divdp_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:)     :: dp_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:,:)   :: dp_star_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: qmin_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: qmax_d
  integer              ,device,allocatable,dimension(:)           :: send_nelem_d
  integer              ,device,allocatable,dimension(:)           :: recv_nelem_d
  integer              ,device,allocatable,dimension(:,:)         :: send_indices_d
  integer              ,device,allocatable,dimension(:,:)         :: recv_indices_d
  integer              ,device,allocatable,dimension(:)           :: recv_internal_indices_d
  integer              ,device,allocatable,dimension(:)           :: recv_external_indices_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: recvbuf_d

  !PINNED Host arrays
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:,:,:) :: qdp_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:,:) :: vstar_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:,:) :: qtens_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:)   :: dp_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:)   :: divdp_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:,:) :: dp_star_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:) :: dp_np1_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:) :: sendbuf_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:) :: recvbuf_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:,:) :: qtens_biharmonic
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:) :: qmin
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:) :: qmax
!  real(kind=real_kind),pinned,allocatable,dimension(:,:,:) :: psdiss_biharmonic_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:) :: dpdiss_biharmonic_h

  !Normal Host arrays
  integer,allocatable,dimension(:)   :: send_nelem
  integer,allocatable,dimension(:)   :: recv_nelem
  integer,allocatable,dimension(:,:) :: send_indices
  integer,allocatable,dimension(:,:) :: recv_indices
  integer,allocatable,dimension(:)   :: recv_internal_indices
  integer,allocatable,dimension(:)   :: recv_external_indices
  integer :: recv_external_nelem
  integer :: recv_internal_nelem
  logical :: old_peu
  logical,allocatable,dimension(:)   :: d2h_done
  logical,allocatable,dimension(:)   :: msg_sent
  logical,allocatable,dimension(:)   :: msg_rcvd
  logical,allocatable,dimension(:)   :: h2d_done
  integer, parameter :: south_px = 1
  integer, parameter :: east_px  = 3
  integer, parameter :: north_px = 4
  integer, parameter :: west_px  = 2
  integer, parameter :: cuda_streams = 16
  integer(kind=cuda_stream_kind) :: streams(0:cuda_streams)
  integer(kind=cuda_stream_kind) :: streams2(0:cuda_streams)
  integer            :: nbuf
  integer            :: nmsg_rcvd
  integer            :: nmsg_sent
  integer, device :: max_neigh_edges_d, max_corner_elem_d

  type(cudaEvent) :: timer1, timer2


contains



  !The point of this is to initialize any data required in other routines of this module as well
  !as to run one initial CUDA kernel just to get those overheads out of the way so that subsequent
  !timing routines are accurage.
  subroutine cuda_mod_init(elem,deriv,hvcoord)
    use edge_mod      , only: initEdgeBuffer
    use schedule_mod  , only: schedule_t, cycle_t, schedule
    use edge_mod      , only: Edgebuffer_t
    use element_mod   , only: element_t
    use derivative_mod, only: derivative_t
    use hybvcoord_mod, only: hvcoord_t
    implicit none
    type(element_t)   , intent(in) :: elem(:)
    type(derivative_t), intent(in) :: deriv
    type(hvcoord_t)   , intent(in) :: hvcoord

    type (Cycle_t),pointer    :: pCycle
    type (Schedule_t),pointer :: pSchedule
    integer                   :: ie , ierr , icycle , iPtr , rank , nSendCycles , nRecvCycles , nlyr , mx_send_len , mx_recv_len , n
    real(kind=real_kind)      :: dinv_t(np , np , 2 , 2)
    type (dim3)               :: griddim , blockdim
    logical,allocatable,dimension(:,:) :: send_elem_mask
    logical,allocatable,dimension(:,:) :: recv_elem_mask
    logical,allocatable,dimension(:)   :: elem_computed
    integer :: total_work

    max_neigh_edges_d = max_neigh_edges
    max_corner_elem_d = max_corner_elem

#if (defined ELEMENT_OPENMP)
    write(*,*) 'ERROR: Do not use ELEMENT_OPENMP and CUDA FORTRAN'
    stop
#endif
!$OMP BARRIER
!$OMP MASTER
    write(*,*) "cuda_mod_init"

    write(*,*) "allocate arrays on device & host"
#ifdef _PREDICT
    pSchedule => Schedule(iam)
#else
    pSchedule => Schedule(1)
#endif
    nlyr = edgeAdv%nlyr
    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles
    mx_send_len = 0
    mx_recv_len = 0
    do icycle=1,nSendCycles
      if (pSchedule%SendCycle(icycle)%lengthP > mx_send_len) mx_send_len = pSchedule%SendCycle(icycle)%lengthP
    enddo
    do icycle=1,nRecvCycles
      if (pSchedule%RecvCycle(icycle)%lengthP > mx_recv_len) mx_recv_len = pSchedule%RecvCycle(icycle)%lengthP 
    enddo
    nbuf=4*(np+max_corner_elem)*nelemd  !inlined from edge_mod.F90, initEdgeBuffer()

    !Allocate the host and device arrays
    allocate( qmin_d                   (nlev,qsize_d                 ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( qmax_d                   (nlev,qsize_d                 ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( qdp_d                    (np,np,nlev,qsize_d,timelevels,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( qtens_d                  (np,np,nlev,qsize_d           ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( qtens_biharmonic_d       (np,np,nlev,qsize_d           ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( spheremp_d               (np,np                        ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( rspheremp_d              (np,np                        ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( dinv_d                   (np,np,2,2                    ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( deriv_dvv_d              (np,np                               ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( variable_hyperviscosity_d(np,np                        ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( metdet_d                 (np,np                        ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( rmetdet_d                (np,np                        ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( dpdiss_biharmonic_d      (np,np,nlev                   ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( vstar_d                  (np,np,nlev,2                 ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( divdp_d                  (np,np,nlev                   ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( dp_d                     (np,np,nlev                   ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( hyai_d                   (      nlev+1                        ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( hybi_d                   (      nlev+1                        ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( dp_star_d                (np,np,nlev,qsize_d           ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( reverse_d                (max_neigh_edges              ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( putmapP_d                (max_neigh_edges              ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( getmapP_d                (max_neigh_edges              ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( edgebuf_d                (nlev*qsize_d,nbuf                   ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( recv_internal_indices_d  (nelemd                              ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( recv_external_indices_d  (nelemd                              ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__

    allocate( qdp_h                    (np,np,nlev,qsize_d,timelevels,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( qmin                     (nlev,qsize_d                 ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( qmax                     (nlev,qsize_d                 ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( vstar_h                  (np,np,nlev,2                 ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( qtens_h                  (np,np,nlev,qsize_d           ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( qtens_biharmonic         (np,np,nlev,qsize_d           ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( dp_h                     (np,np,nlev                   ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( divdp_h                  (np,np,nlev                   ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( dp_star_h                (np,np,nlev,qsize_d           ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( dpdiss_biharmonic_h      (np,np,nlev                   ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( dp_np1_h                 (np,np,nlev                   ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__

    ! The PGI compiler with cuda enabled errors when allocating arrays of zero
    !   size - here when using only one MPI task
    if (nSendCycles > 0) then
      allocate( send_nelem_d             (       nSendCycles                  ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
      allocate( send_indices_d           (nelemd,nSendCycles                  ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
      allocate( sendbuf_h                (nlev*qsize_d,mx_send_len,nSendCycles) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
      allocate( send_elem_mask           (nelemd,nSendCycles                  ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
      allocate( send_nelem               (       nSendCycles                  ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
      allocate( send_indices             (nelemd,nSendCycles                  ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
      allocate( send_elem_mask           (nelemd,nSendCycles                  ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
      allocate( d2h_done                 (nSendCycles                         ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
      allocate( msg_sent                 (nSendCycles                         ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    endif

    if (nRecvCycles > 0) then
      allocate( recv_nelem_d             (       nRecvCycles                  ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
      allocate( recv_indices_d           (nelemd,nRecvCycles                  ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
      allocate( recvbuf_h                (nlev*qsize_d,mx_recv_len,nRecvCycles) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
      allocate( recvbuf_d                (nlev*qsize_d,mx_recv_len,nRecvCycles) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
      allocate( recv_elem_mask           (nelemd,nRecvCycles                  ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
      allocate( recv_nelem               (       nRecvCycles                  ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
      allocate( recv_indices             (nelemd,nRecvCycles                  ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
      allocate( recv_elem_mask           (nelemd,nRecvCycles                  ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
      allocate( h2d_done                 (nRecvCycles                         ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
      allocate( msg_rcvd                 (nRecvCycles                         ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    endif 

    allocate( recv_internal_indices    (nelemd                              ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( recv_external_indices    (nelemd                              ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( elem_computed            (nelemd                              ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    write(*,*) "send data from host to device"
    !Copy over data to the device
    ierr = cudaMemcpy( deriv_dvv_d , deriv%dvv    , size( deriv%dvv    ) , cudaMemcpyHostToDevice )
    ierr = cudaMemcpy( hyai_d      , hvcoord%hyai , size( hvcoord%hyai ) , cudaMemcpyHostToDevice )
    ierr = cudaMemcpy( hybi_d      , hvcoord%hybi , size( hvcoord%hybi ) , cudaMemcpyHostToDevice )
    do ie = 1,nelemd
      dinv_t(:,:,1,1) = elem(ie)%dinv(1,1,:,:)
      dinv_t(:,:,1,2) = elem(ie)%dinv(1,2,:,:)
      dinv_t(:,:,2,1) = elem(ie)%dinv(2,1,:,:)
      dinv_t(:,:,2,2) = elem(ie)%dinv(2,2,:,:)
      ierr = cudaMemcpy( qdp_d                    (1,1,1,1,1,ie) , elem(ie)%state%Qdp               , size(elem(ie)%state%Qdp              ) , cudaMemcpyHostToDevice ) ; if ( ierr .ne. 0 ) stop __LINE__
      ierr = cudaMemcpy( spheremp_d               (1,1      ,ie) , elem(ie)%spheremp                , size(elem(ie)%spheremp               ) , cudaMemcpyHostToDevice ) ; if ( ierr .ne. 0 ) stop __LINE__
      ierr = cudaMemcpy( rspheremp_d              (1,1      ,ie) , elem(ie)%rspheremp               , size(elem(ie)%rspheremp              ) , cudaMemcpyHostToDevice ) ; if ( ierr .ne. 0 ) stop __LINE__
      ierr = cudaMemcpy( dinv_d                   (1,1,1,1  ,ie) , dinv_t                           , size(dinv_t                          ) , cudaMemcpyHostToDevice ) ; if ( ierr .ne. 0 ) stop __LINE__
      ierr = cudaMemcpy( metdet_d                 (1,1      ,ie) , elem(ie)%metdet                  , size(elem(ie)%metdet                 ) , cudaMemcpyHostToDevice ) ; if ( ierr .ne. 0 ) stop __LINE__
      ierr = cudaMemcpy( rmetdet_d                (1,1      ,ie) , elem(ie)%rmetdet                 , size(elem(ie)%rmetdet                ) , cudaMemcpyHostToDevice ) ; if ( ierr .ne. 0 ) stop __LINE__
      ierr = cudaMemcpy( putmapP_d                (1        ,ie) , elem(ie)%desc%putmapP            , size(elem(ie)%desc%putmapP           ) , cudaMemcpyHostToDevice ) ; if ( ierr .ne. 0 ) stop __LINE__
      ierr = cudaMemcpy( getmapP_d                (1        ,ie) , elem(ie)%desc%getmapP            , size(elem(ie)%desc%getmapP           ) , cudaMemcpyHostToDevice ) ; if ( ierr .ne. 0 ) stop __LINE__
      ierr = cudaMemcpy( reverse_d                (1        ,ie) , elem(ie)%desc%reverse            , size(elem(ie)%desc%reverse           ) , cudaMemcpyHostToDevice ) ; if ( ierr .ne. 0 ) stop __LINE__
      ierr = cudaMemcpy( variable_hyperviscosity_d(1,1      ,ie) , elem(ie)%variable_hyperviscosity , size(elem(ie)%variable_hyperviscosity) , cudaMemcpyHostToDevice ) ; if ( ierr .ne. 0 ) stop __LINE__
    enddo

    write(*,*) "edgebuffers"
    !These have to be in a threaded region or they complain and die
!$OMP END MASTER
    call initEdgeBuffer(edgeAdv   ,qsize_d*nlev  )
    call initEdgeBuffer(edgeAdvDSS,      nlev  )
    call initEdgeBuffer(edgeAdvQ2 ,qsize_d*nlev*2)
    call initEdgeBuffer(edgeAdvQ3 ,qsize_d*nlev*3)
!$OMP MASTER

    write(*,*) "initial kernel"
    !This needs to run because we need accurate timing during out cuda profiler runs. Initial kernel runs incur overheads, so we do this here
    blockdim = dim3(1,1,1)
    griddim  = dim3(1,1,1)
    call warmup <<< griddim , blockdim >>> ( ie )
    ierr = cudaThreadSynchronize()

    do n = 0 , cuda_streams
      ierr = cudaStreamCreate(streams(n))
      ierr = cudaStreamCreate(streams2(n))
    enddo
    ierr = cudaDeviceSetCacheConfig(cudaFuncCachePreferShared)
    ierr = cudaEventCreate(timer1)
    ierr = cudaEventCreate(timer2)

    write(*,*) "Dividing elements among cycles in which they participate"
    !For efficient MPI, PCI-e, packing, and unpacking, we need to separate out the cycles by dependence. Once on cycle has packed, then stage the PCI-e D2H, MPI, PCI-e H2D, & internal unpack
    !We begin by testing what elements contribute to packing in what cycle's MPI data.
    do ie = 1,nelemd
      send_elem_mask(ie,:) = .false.
      do icycle = 1 , nSendCycles
        do n = 1 , max_neigh_edges
          if ( elem(ie)%desc%putmapP(n) >= pSchedule%SendCycle(icycle)%ptrP .and. &
               elem(ie)%desc%putmapP(n) <= pSchedule%SendCycle(icycle)%ptrP + pSchedule%SendCycle(icycle)%lengthP-1 ) then
            send_elem_mask(ie,icycle) = .true.
          endif
        enddo
      enddo
      recv_elem_mask(ie,:) = .false.
      do icycle = 1 , nRecvCycles
        do n = 1 , max_neigh_edges
          if ( elem(ie)%desc%getmapP(n) >= pSchedule%RecvCycle(icycle)%ptrP .and. &
               elem(ie)%desc%getmapP(n) <= pSchedule%RecvCycle(icycle)%ptrP + pSchedule%RecvCycle(icycle)%lengthP-1 ) then
            recv_elem_mask(ie,icycle) = .true.
          endif
        enddo
      enddo
    enddo
    edgebuf_d = 0.

    elem_computed = .false.   !elem_computed tells us whether an element has been touched by a cycle
    !This pass accumulates for each cycle incides participating in the MPI_Isend
    do icycle = 1 , nSendCycles
      send_nelem(icycle) = 0
      do ie = 1 , nelemd
        if ( send_elem_mask(ie,icycle) ) then
          send_nelem(icycle) = send_nelem(icycle) + 1
          send_indices(send_nelem(icycle),icycle) = ie
          elem_computed(ie) = .true.
        endif
      enddo
    enddo
    total_work = sum(send_nelem)
    do ie = 1 , nelemd
      if (.not. elem_computed(ie)) total_work = total_work + 1
    enddo
    !This pass adds to each cycle the internal elements not participating in MPI_Isend, so as to even distribute them across cycles.
    do icycle = 1 , nSendCycles
      do ie = 1 , nelemd
        if ( .not. elem_computed(ie) .and. send_nelem(icycle) < int(ceiling(total_work/dble(nSendCycles))) ) then
          send_nelem(icycle) = send_nelem(icycle) + 1
          send_indices(send_nelem(icycle),icycle) = ie
          elem_computed(ie) = .true.
        endif
      enddo
    enddo

    elem_computed = .false.
    !This pass accumulates for each cycle incides participating in the MPI_Irecv
    do icycle = 1 , nRecvCycles
      recv_nelem(icycle) = 0
      do ie = 1 , nelemd
        if ( recv_elem_mask(ie,icycle) .and. ( .not. elem_computed(ie) ) ) then
          recv_nelem(icycle) = recv_nelem(icycle) + 1
          recv_indices(recv_nelem(icycle),icycle) = ie
          elem_computed(ie) = .true.
        endif
      enddo
    enddo
    !This pass accumulates all elements from all cycles participating in MPI_Irecv into the recv_external_indices array
    recv_external_nelem = 0
    do icycle = 1 , nRecvCycles
      do ie = 1 , recv_nelem(icycle)
        recv_external_nelem = recv_external_nelem + 1
        recv_external_indices(recv_external_nelem) = recv_indices(ie,icycle)
      enddo
    enddo
    !This pass goes through all elements, and distributes evenly the elements not participating in MPI_Irecv 
    recv_internal_nelem = 0
    do ie = 1 , nelemd
      if ( .not. elem_computed(ie) ) then
        recv_internal_nelem = recv_internal_nelem + 1
        recv_internal_indices(recv_internal_nelem) = ie
      endif
    enddo
    !This pass adds to each cycle the internal elements not participating in MPI_Irecv, so as to even distribute them across cycles.
    do icycle = 1 , nRecvCycles
      do ie = 1 , nelemd
        if ( .not. elem_computed(ie) .and. recv_nelem(icycle) < int(ceiling(nelemd/dble(nRecvCycles))) ) then
          recv_nelem(icycle) = recv_nelem(icycle) + 1
          recv_indices(recv_nelem(icycle),icycle) = ie
          elem_computed(ie) = .true.
        endif
      enddo
    enddo

    old_peu = .false.
    do icycle = 1 , nSendCycles
      if (send_nelem(icycle) == 0) then
        write(*,*) 'WARNING: ZERO ELEMENT CYCLES EXIST. A BETTER DECOMPOSITION WILL RUN FASTER IN THE PACK-EXCHANGE-UNPACK.'
        old_peu = .true.
      endif
    enddo

    write(*,*) "Sending element & cycle informationt to device"
    ierr = cudaMemcpy(send_nelem_d           ,send_nelem           ,size(send_nelem           ),cudaMemcpyHostToDevice); if(ierr.ne.0) stop __LINE__
    ierr = cudaMemcpy(recv_nelem_d           ,recv_nelem           ,size(recv_nelem           ),cudaMemcpyHostToDevice); if(ierr.ne.0) stop __LINE__
    ierr = cudaMemcpy(send_indices_d         ,send_indices         ,size(send_indices         ),cudaMemcpyHostToDevice); if(ierr.ne.0) stop __LINE__
    ierr = cudaMemcpy(recv_indices_d         ,recv_indices         ,size(recv_indices         ),cudaMemcpyHostToDevice); if(ierr.ne.0) stop __LINE__
    ierr = cudaMemcpy(recv_internal_indices_d,recv_internal_indices,size(recv_internal_indices),cudaMemcpyHostToDevice); if(ierr.ne.0) stop __LINE__
    ierr = cudaMemcpy(recv_external_indices_d,recv_external_indices,size(recv_external_indices),cudaMemcpyHostToDevice); if(ierr.ne.0) stop __LINE__

    write(*,*)"done cuda_mod_init"
!$OMP END MASTER
!$OMP BARRIER
  end subroutine cuda_mod_init



  !Meaningless kernel just to get initial kernel overheads out of the way.
  attributes(global) subroutine warmup(a)
    integer,value :: a
    a = 2.0 * a
  end subroutine warmup



  subroutine copy_qdp_h2d( elem , nt )
    use element_mod, only: element_t
    implicit none
    type(element_t), intent(in   ) :: elem(:)
    integer        , intent(in   ) :: nt
    integer :: ierr , ie
!$OMP MASTER
    do ie = 1,nelemd
      ierr = cudaMemcpyAsync( qdp_d(1,1,1,1,nt,ie) , elem(ie)%state%qdp(1,1,1,1,nt) , size(elem(ie)%state%qdp(:,:,:,:,nt)) , cudaMemcpyHostToDevice , streams(1) )
    enddo
!$OMP END MASTER
  end subroutine copy_qdp_h2d



  subroutine copy_qdp_d2h( elem , nt )
    use element_mod, only: element_t
    implicit none
    type(element_t), intent(in   ) :: elem(:)
    integer        , intent(in   ) :: nt
    integer :: ierr , ie
!$OMP BARRIER
!$OMP MASTER
    ierr = cudaThreadSynchronize()
    do ie = 1,nelemd
      ierr = cudaMemcpy( elem(ie)%state%qdp(1,1,1,1,nt) , qdp_d(1,1,1,1,nt,ie) , size(elem(ie)%state%qdp(:,:,:,:,nt)) , cudaMemcpyDeviceToHost )
    enddo
    ierr = cudaThreadSynchronize()
!$OMP END MASTER
!$OMP BARRIER
  end subroutine copy_qdp_d2h



  subroutine euler_step_cuda( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
  use kinds             , only: real_kind
  use dimensions_mod    , only: np, npdg, nlev, qsize
  use hybrid_mod        , only: hybrid_t
  use element_mod       , only: element_t
  use derivative_mod    , only: derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
  use edge_mod          , only: edgevpack, edgevunpack
  use bndry_mod         , only: bndry_exchangev
  use hybvcoord_mod     , only: hvcoord_t
  use control_mod       , only: nu_q, nu_p, limiter_option
  use perf_mod          , only: t_startf, t_stopf  ! _EXTERNAL
  use viscosity_mod     , only: biharmonic_wk_scalar, biharmonic_wk_scalar_minmax, neighbor_minmax
  implicit none
  integer              , intent(in   )         :: np1_qdp, n0_qdp
  real (kind=real_kind), intent(in   )         :: dt
  type (element_t)     , intent(inout), target :: elem(:)
  type (hvcoord_t)     , intent(in   )         :: hvcoord
  type (hybrid_t)      , intent(in   )         :: hybrid
  type (derivative_t)  , intent(in   )         :: deriv
  integer              , intent(in   )         :: nets
  integer              , intent(in   )         :: nete
  integer              , intent(in   )         :: DSSopt
  integer              , intent(in   )         :: rhs_multiplier

  ! local
  real(kind=real_kind), dimension(np,np                       ) :: divdp, dpdiss
  real(kind=real_kind), dimension(np,np,2                     ) :: gradQ
  real(kind=real_kind), dimension(np,np,2,nlev                ) :: Vstar
  real(kind=real_kind), dimension(np,np  ,nlev                ) :: dp,dp_star
  real(kind=real_kind), pointer, dimension(:,:,:)               :: DSSvar
  real(kind=real_kind) :: dp0,qmintmp,qmaxtmp
  integer :: ie,q,i,j,k
  integer :: rhs_viss = 0

  integer :: ierr
  type(dim3) :: blockdim , griddim

  call t_startf('euler_step')

  if (rhs_multiplier == 0) then
    call copy_qdp_h2d( elem , n0_qdp)
  endif

  if (limiter_option == 8) then
    write(*,*) 'CUDA_MOD IS NOT INTENDED FOR USE WITH LIMITER_OPTION == 8 AT THIS TIME!'
    write(*,*) 'PLEASE USE LIMITER_OPTION == 0 WHEN THE GPU OPTION IS ENABLED!'
    stop
  endif
  if (nu_p > 0) then
    write(*,*) 'CUDA_MOD IS NOT INTENDED FOR USE WITH NU_P > 0 AT THIS TIME!'
    write(*,*) 'PLEASE USE NU_P == 0 WHEN THE GPU OPTION IS ENABLED!'
    stop
  endif

  rhs_viss = 0
  !   2D Advection step
  do ie = nets , nete
    ! Compute velocity used to advance Qdp 
    do k = 1 , nlev    !  Loop index added (AAM)
      ! derived variable divdp_proj() (DSS'd version of divdp) will only be correct on 2nd and 3rd stage
      ! but that's ok because rhs_multiplier=0 on the first stage:
      dp_h(:,:,k,ie) = elem(ie)%derived%dp(:,:,k) - rhs_multiplier * dt * elem(ie)%derived%divdp_proj(:,:,k) 
      vstar_h(:,:,k,1,ie) = elem(ie)%derived%vn0(:,:,1,k) / dp_h(:,:,k,ie)
      vstar_h(:,:,k,2,ie) = elem(ie)%derived%vn0(:,:,2,k) / dp_h(:,:,k,ie)
    enddo
  enddo
!$OMP BARRIER
!$OMP MASTER
  ierr = cudaMemcpyAsync( dp_d    , dp_h    , size( dp_h    ) , cudaMemcpyHostToDevice , streams(1) )
  ierr = cudaMemcpyAsync( vstar_d , vstar_h , size( vstar_h ) , cudaMemcpyHostToDevice , streams(1) )
  blockdim = dim3( np      , np     , nlev )
  griddim  = dim3( qsize_d , nelemd , 1    )
  call euler_step_kernel1<<<griddim,blockdim,0,streams(1)>>>( qdp_d , spheremp_d , qmin_d , qmax_d , dp_d , vstar_d , divdp_d , hybi_d ,               &
                                                              dpdiss_biharmonic_d , qtens_biharmonic_d , metdet_d , rmetdet_d , dinv_d , deriv_dvv_d , &
                                                              n0_qdp , np1_qdp , rhs_viss , dt , nu_p , nu_q , limiter_option , 1 , nelemd )
  if ( limiter_option == 4 ) then
    blockdim = dim3( np      , np     , nlev )
    griddim  = dim3( qsize_d , nelemd , 1    )
    call limiter2d_zero_kernel<<<griddim,blockdim,0,streams(1)>>>( qdp_d , 1 , nelemd , np1_qdp )
  endif
!$OMP END MASTER
  !This is a departure from the original order, adding an extra MPI communication. It's advantageous because it simplifies
  !the Pack-Exchange-Unpack procedure for us, since we're adding complexity to overlap MPI and packing
  if ( DSSopt /= DSSno_var ) then
    do ie = nets , nete
      if ( DSSopt == DSSeta         ) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
      if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega_p(:,:,:)
      if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)
      do k = 1 , nlev
        DSSvar(:,:,k) = elem(ie)%spheremp(:,:)*DSSvar(:,:,k) 
      enddo
      call edgeVpack(edgeAdvDSS,DSSvar(:,:,1:nlev),nlev,0,elem(ie)%desc)
    enddo
    call bndry_exchangeV(hybrid,edgeAdvDSS)
    do ie = nets , nete
      if ( DSSopt == DSSeta         ) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
      if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega_p(:,:,:)
      if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)
      call edgeVunpack(edgeAdvDSS,DSSvar(:,:,1:nlev),nlev,0,elem(ie)%desc)
      do k = 1 , nlev
        DSSvar(:,:,k)=DSSvar(:,:,k)*elem(ie)%rspheremp(:,:)
      enddo
    enddo
  endif
!$OMP MASTER
  call pack_exchange_unpack_stage(np1_qdp,hybrid,qdp_d,timelevels)
  blockdim = dim3( np*np   , nlev   , 1 )
  griddim  = dim3( qsize_d , nelemd , 1 )
  call euler_hypervis_kernel_last<<<griddim,blockdim>>>( qdp_d , rspheremp_d , 1 , nelemd , np1_qdp )
  ierr = cudaThreadSynchronize()
!$OMP END MASTER
!$OMP BARRIER

  call t_stopf('euler_step')
end subroutine euler_step_cuda



subroutine qdp_time_avg_cuda( elem , rkstage , n0_qdp , np1_qdp , limiter_option , nu_p , nets , nete )
  use element_mod, only: element_t
  implicit none
  type(element_t)     , intent(inout) :: elem(:)
  real(kind=real_kind), intent(in   ) :: nu_p
  integer             , intent(in   ) :: rkstage , n0_qdp , np1_qdp , nets , nete, limiter_option
  type(dim3) :: griddim , blockdim
  integer :: ierr, ie
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP BARRIER
!$OMP MASTER
  ierr = cudaThreadSynchronize()
  blockdim = dim3( np      , np     , nlev )
  griddim  = dim3( qsize_d , nelemd , 1    )
  call qdp_time_avg_kernel<<<griddim,blockdim>>>( qdp_d , rkstage , n0_qdp , np1_qdp , 1 , nelemd )
  ierr = cudaThreadSynchronize()
!$OMP END MASTER
!$OMP BARRIER
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine qdp_time_avg_cuda



subroutine advance_hypervis_scalar_cuda( edgeAdv , elem , hvcoord , hybrid , deriv , nt , nt_qdp , nets , nete , dt2 )
  !  hyperviscsoity operator for foward-in-time scheme
  !  take one timestep of:  
  !          Q(:,:,:,np) = Q(:,:,:,np) +  dt2*nu*laplacian**order ( Q )
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  use kinds          , only : real_kind
  use dimensions_mod , only : np, nlev
  use hybrid_mod     , only : hybrid_t
  use hybvcoord_mod  , only : hvcoord_t
  use element_mod    , only : element_t
  use derivative_mod , only : derivative_t
  use edge_mod       , only : EdgeBuffer_t, edgevpack, edgevunpack
  use bndry_mod      , only : bndry_exchangev
  use perf_mod       , only : t_startf, t_stopf, t_barrierf
  use control_mod    , only : rsplit, nu_q, hypervis_order, hypervis_subcycle_q
  implicit none
  type (EdgeBuffer_t)  , intent(inout)         :: edgeAdv
  type (element_t)     , intent(inout), target :: elem(:)
  type (hvcoord_t)     , intent(in   )         :: hvcoord
  type (hybrid_t)      , intent(in   )         :: hybrid
  type (derivative_t)  , intent(in   )         :: deriv
  integer              , intent(in   )         :: nt
  integer              , intent(in   )         :: nt_qdp
  integer              , intent(in   )         :: nets
  integer              , intent(in   )         :: nete
  real (kind=real_kind), intent(in   )         :: dt2
  ! local
  integer :: k , kptr , i , j , ie , ic , q , ierr
  real (kind=real_kind) :: dt
  type(dim3) :: griddim , blockdim
  if ( nu_q           == 0 ) return
  if ( hypervis_order /= 2 ) return
  call t_barrierf('sync_advance_hypervis_scalar', hybrid%par%comm)
  call t_startf('advance_hypervis_scalar_cuda')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dt = dt2 / hypervis_subcycle_q
  do ic = 1 , hypervis_subcycle_q
    ! Qtens = Q/dp   (apply hyperviscsoity to dp0 * Q, not Qdp)
    do ie = nets , nete
      do k = 1 , nlev
        dp_h(:,:,k,ie) = elem(ie)%derived%dp(:,:,k) - dt2*elem(ie)%derived%divdp_proj(:,:,k)
      enddo
    enddo
!$OMP BARRIER
!$OMP MASTER
    ierr = cudaMemcpyAsync( dp_d , dp_h , size( dp_h ) , cudaMemcpyHostToDevice , streams(1) )
    !KERNEL 1
    blockdim = dim3( np     , np , nlev )
    griddim  = dim3( nelemd , 1  , 1    )
    call hypervis_kernel1<<<griddim,blockdim,0,streams(1)>>>( qdp_d , qtens_d , dp_d , dinv_d , variable_hyperviscosity_d , spheremp_d , deriv_dvv_d , nets , nete , dt , nt_qdp )
    ierr = cudaThreadSynchronize()

    call pack_exchange_unpack_stage(1,hybrid,qtens_d,1)
    ierr = cudaThreadSynchronize()

    !KERNEL 2
    blockdim = dim3( np     , np , nlev )
    griddim  = dim3( nelemd , 1  , 1    )
    call hypervis_kernel2<<<griddim,blockdim,0,streams(1)>>>( qdp_d , qtens_d , dp_d , dinv_d , variable_hyperviscosity_d , spheremp_d , rspheremp_d , deriv_dvv_d , hyai_d , hybi_d , hvcoord%ps0 , nu_q , nets , nete , dt , nt_qdp )
    blockdim = dim3( np, np, nlev )
    griddim  = dim3( qsize_d, nelemd , 1 )
    call limiter2d_zero_kernel<<<griddim,blockdim,0,streams(1)>>>(qdp_d,nets,nete,nt_qdp)
    ierr = cudaThreadSynchronize()

    call pack_exchange_unpack_stage(nt_qdp,hybrid,qdp_d,timelevels)
    ierr = cudaThreadSynchronize()

    !KERNEL 3
    blockdim = dim3( np * np , nlev   , 1 )
    griddim  = dim3( qsize_d , nelemd , 1 )
    call euler_hypervis_kernel_last<<<griddim,blockdim>>>( qdp_d , rspheremp_d , nets , nete , nt_qdp )
    ierr = cudaThreadSynchronize()
!$OMP END MASTER
!$OMP BARRIER
  enddo

  call t_stopf('advance_hypervis_scalar_cuda')
  call copy_qdp_d2h( elem , nt_qdp )

end subroutine advance_hypervis_scalar_cuda



attributes(global) subroutine qdp_time_avg_kernel( qdp , rkstage , n0_qdp , np1_qdp , nets , nete )
  implicit none
  real(kind=real_kind), intent(inout) :: qdp( np , np , nlev , qsize_d , timelevels , nets:nete )
  integer, value      , intent(in   ) :: rkstage , n0_qdp , np1_qdp , nets , nete
  integer :: i, j, k, q, ie
  i  = threadidx%x
  j  = threadidx%y
  k  = threadidx%z
  q  = blockidx%x
  ie = blockidx%y
  qdp(i,j,k,q,np1_qdp,ie) = ( qdp(i,j,k,q,n0_qdp,ie) + (rkstage-1) * qdp(i,j,k,q,np1_qdp,ie) ) / rkstage
end subroutine qdp_time_avg_kernel



attributes(global) subroutine euler_step_kernel1( Qdp , spheremp , qmin , qmax , dp , vstar , divdp , hybi ,                   &
                                                  dpdiss_biharmonic , qtens_biharmonic , metdet , rmetdet , dinv , deriv_dvv , &
                                                  n0_qdp , np1_qdp , rhs_viss , dt , nu_p , nu_q , limiter_option , nets , nete )
  implicit none
  real(kind=real_kind), dimension(np,np,nlev,qsize_d,timelevels,nets:nete), intent(inout) :: Qdp
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: spheremp
  real(kind=real_kind), dimension(      nlev,qsize_d           ,nets:nete), intent(in   ) :: qmin
  real(kind=real_kind), dimension(      nlev,qsize_d           ,nets:nete), intent(in   ) :: qmax
  real(kind=real_kind), dimension(np,np,nlev                   ,nets:nete), intent(in   ) :: dp
  real(kind=real_kind), dimension(np,np,nlev,2                 ,nets:nete), intent(in   ) :: vstar
  real(kind=real_kind), dimension(np,np,nlev                   ,nets:nete), intent(in   ) :: divdp
  real(kind=real_kind), dimension(      nlev+1                           ), intent(in   ) :: hybi
  real(kind=real_kind), dimension(np,np,nlev                   ,nets:nete), intent(in   ) :: dpdiss_biharmonic
  real(kind=real_kind), dimension(np,np,nlev,qsize_d           ,nets:nete), intent(in   ) :: qtens_biharmonic
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: metdet
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: rmetdet
  real(kind=real_kind), dimension(np,np,2,2                    ,nets:nete), intent(in   ) :: dinv
  real(kind=real_kind), dimension(np,np                                  ), intent(in   ) :: deriv_dvv
  real(kind=real_kind), value                                             , intent(in   ) :: dt, nu_q, nu_p
  integer, value                                                          , intent(in   ) :: n0_qdp, np1_qdp, rhs_viss, limiter_option, nets, nete
  integer :: i , j , k , q , ie , ij , ijk
  real(kind=real_kind), dimension(np*np+1,nlev,2), shared :: vstar_s
  real(kind=real_kind), dimension(np*np+1,nlev  ), shared :: qtens_s
  real(kind=real_kind), dimension(np*np+1       ), shared :: spheremp_s
  real(kind=real_kind), dimension(np*np+1       ), shared :: deriv_dvv_s
  real(kind=real_kind), dimension(np*np+1       ), shared :: metdet_s
  real(kind=real_kind), dimension(np*np+1       ), shared :: rmetdet_s
  real(kind=real_kind), dimension(        nlev  ), shared :: qmin_s
  real(kind=real_kind), dimension(        nlev  ), shared :: qmax_s
  real(kind=real_kind), dimension(np*np+1,nlev  ), shared :: dp_star_s
  real(kind=real_kind), dimension(        nlev+1), shared :: hybi_s
  real(kind=real_kind) :: qtmp

  !Define the indices
  i  = threadidx%x
  j  = threadidx%y
  k  = threadidx%z
  q  = blockidx%x
  ie = blockidx%y
  ij  =             (j-1)*np+i
  ijk = (k-1)*np*np+(j-1)*np+i

  !Pre-load shared variables
  vstar_s(ij,k,:) = vstar(i,j,k,:,ie)
  if (ijk <= nlev) then
    qmin_s(ijk) = qmin(ijk,q,ie)
    qmax_s(ijk) = qmax(ijk,q,ie)
  endif
  if (ijk <= nlev+1) then
    hybi_s(ijk) = hybi(ijk)
  endif
  if (k == 1) then
    spheremp_s(ij) = spheremp(i,j,ie)
    metdet_s(ij) = metdet(i,j,ie)
    rmetdet_s(ij) = rmetdet(i,j,ie)
    deriv_dvv_s(ij) = deriv_dvv(i,j)
  endif
  call syncthreads()

  !Begin the kernel
  qtmp = Qdp(i,j,k,q,n0_qdp,ie)
  qtens_s(ij,k) = qtmp - dt * divergence_sphere( i , j , ie , k , ij , vstar_s , qtmp , metdet_s , rmetdet_s , dinv , deriv_dvv_s , nets , nete )
  if ( rhs_viss /= 0 ) qtens_s(ij,k) = qtens_s(ij,k) + qtens_biharmonic(i,j,k,q,ie)
  call syncthreads()
  Qdp(i,j,k,q,np1_qdp,ie) = spheremp_s(ij) * Qtens_s(ij,k)
end subroutine euler_step_kernel1



attributes(device) function divergence_sphere(i,j,ie,k,ij,Vstar,qtmp,metdet,rmetdet,dinv,deriv_dvv,nets,nete) result(dp_star)
  use physical_constants, only: rrearth
  implicit none
  integer,              intent(in) :: i, j, ie, k, ij, nets, nete
  real(kind=real_kind), intent(in) :: Dinv     (np*np,2,2,nets:nete)
  real(kind=real_kind), intent(in) :: metdet   (np*np+1)
  real(kind=real_kind), intent(in) :: rmetdet  (np*np+1)
  real(kind=real_kind), intent(in) :: deriv_dvv(np*np+1)
  real(kind=real_kind), intent(in) :: Vstar    (np*np+1,nlev,2)
  real(kind=real_kind), intent(in), value :: qtmp
  real(kind=real_kind)             :: dp_star
  real(kind=real_kind), shared :: gv(np*np,nlev,2)
  integer :: s
  real(kind=real_kind) :: vvtemp, divtemp
  gv(ij,k,1) = metdet(ij) * ( dinv(ij,1,1,ie) * Vstar(ij,k,1) * qtmp + dinv(ij,1,2,ie) * Vstar(ij,k,2) * qtmp )
  gv(ij,k,2) = metdet(ij) * ( dinv(ij,2,1,ie) * Vstar(ij,k,1) * qtmp + dinv(ij,2,2,ie) * Vstar(ij,k,2) * qtmp )
  call syncthreads()
  divtemp = 0.0d0
  vvtemp   = 0.0d0
  do s = 1 , np
    divtemp = divtemp + deriv_dvv((i-1)*np+s) * gv((j-1)*np+s,k,1)
    vvtemp  = vvtemp  + deriv_dvv((j-1)*np+s) * gv((s-1)*np+i,k,2)
  end do
  dp_star = ( divtemp + vvtemp ) * ( rmetdet(ij) * rrearth )
end function divergence_sphere



attributes(global) subroutine limiter2d_zero_kernel(Qdp,nets,nete,np1)
  use kinds, only : real_kind
  use dimensions_mod, only : np, nlev
  implicit none
  real (kind=real_kind), intent(inout) :: Qdp(np,np,nlev,qsize_d,timelevels,nets:nete)
  integer, value       , intent(in   ) :: nets,nete,np1
  integer :: i, j, k, q, ie, jj, tid, ind
  real (kind=real_kind) :: mass,mass_new
  real (kind=real_kind), shared :: Qdp_shared((np*np+PAD)*nlev)
  real (kind=real_kind), shared :: mass_shared(nlev)
  real (kind=real_kind), shared :: mass_new_shared(nlev)

  i  = threadidx%x
  j  = threadidx%y
  k  = threadidx%z
  q  = blockidx%x
  ie = blockidx%y

  tid = (threadidx%z-1)*(np*np    ) + (threadidx%y-1)*(np) + threadidx%x
  ind = (threadidx%z-1)*(np*np+PAD) + (threadidx%y-1)*(np) + threadidx%x

  Qdp_shared(ind) = Qdp(i,j,k,q,np1,ie)
  call syncthreads()

  if ( tid <= nlev ) then
    mass = 0.
    do jj = 1 , np*np
      mass = mass + Qdp_shared((tid-1)*(np*np+PAD)+jj)
    enddo
    mass_shared(tid) = mass
  endif
  call syncthreads()

  if ( mass_shared(k)  < 0 ) Qdp_shared(ind) = -Qdp_shared(ind)
  if ( Qdp_shared(ind) < 0 ) Qdp_shared(ind) = 0
  call syncthreads()

  if ( tid <= nlev ) then
    mass = 0.
    do jj = 1 , np*np
      mass = mass + Qdp_shared((tid-1)*(np*np+PAD)+jj)
    enddo
    mass_new_shared(tid) = mass
  endif
  call syncthreads()

  ! now scale the all positive values to restore mass
  if ( mass_new_shared(k) > 0 ) Qdp_shared(ind) =  Qdp_shared(ind) * abs(mass_shared(k)) / mass_new_shared(k)
  if ( mass_shared    (k) < 0 ) Qdp_shared(ind) = -Qdp_shared(ind)
  Qdp(i,j,k,q,np1,ie) = Qdp_shared(ind)
end subroutine limiter2d_zero_kernel



attributes(global) subroutine euler_hypervis_kernel_last( Qdp , rspheremp , nets , nete , np1 )
  implicit none
  real(kind=real_kind), dimension(np*np,nlev,qsize_d,timelevels,nets:nete), intent(inout) :: Qdp
  real(kind=real_kind), dimension(np*np                        ,nets:nete), intent(in   ) :: rspheremp
  integer, value                                                          , intent(in   ) :: nets , nete , np1
  integer :: i, k, q, ie
  i  = threadidx%x
  k  = threadidx%y
  q  = blockidx%x
  ie = blockidx%y
  Qdp(i,k,q,np1,ie) = rspheremp(i,ie) * Qdp(i,k,q,np1,ie)
end subroutine euler_hypervis_kernel_last



subroutine pack_exchange_unpack_stage(np1,hybrid,array_in,tl_in)
  use hybrid_mod, only : hybrid_t
  use schedule_mod, only : schedule_t, schedule, cycle_t
  use parallel_mod, only : abortmp, status, srequest, rrequest, mpireal_t, mpiinteger_t, iam
  implicit none
  include 'mpif.h'
  type(hybrid_t)              , intent(in   ) :: hybrid
  real(kind=real_kind), device, intent(inout) :: array_in(np,np,nlev,qsize_d,tl_in,nelemd)
  integer, value              , intent(in   ) :: np1 , tl_in
  ! local
  type(dim3)                :: griddim6,blockdim6
  integer                   :: icycle,nSendCycles,nRecvCycles,n, ierr
  type (Schedule_t),pointer :: pSchedule
  type (Cycle_t),pointer    :: pCycle
  integer                   :: dest,length,tag,iptr,source,nlyr,query_sum, npacked
  logical :: recvflag, internal_unpacked
  real :: time_milli
#ifdef _PREDICT
   pSchedule => Schedule(iam)
#else
   pSchedule => Schedule(1)
#endif
   nlyr = edgeAdv%nlyr
   nSendCycles = pSchedule%nSendCycles
   nRecvCycles = pSchedule%nRecvCycles
   d2h_done = .false.
   msg_sent = .false.
   h2d_done = .false.
   internal_unpacked = .false.
   nmsg_rcvd = 0
   nmsg_sent = 0


  ierr = cudaThreadSynchronize()
  do icycle=1,nRecvCycles
    pCycle => pSchedule%RecvCycle(icycle)
    source  = pCycle%source - 1
    length  = nlyr * pCycle%lengthP
    tag     = pCycle%tag
    call MPI_Irecv(recvbuf_h(1,1,icycle),length,MPIreal_t,source,tag,hybrid%par%comm,Rrequest(icycle),ierr)
  enddo
  if (old_peu) then
    blockdim6 = dim3( np      , np     , nlev )
    griddim6  = dim3( qsize_d , nelemd , 1    )
    call edgeVpack_kernel<<<griddim6,blockdim6>>>(edgebuf_d,array_in,putmapP_d,reverse_d,nbuf,0,1,nelemd,np1,tl_in)
    ierr = cudaThreadSynchronize()
  else
    do icycle = 1 , nSendCycles
      blockdim6 = dim3( np      , np                 , nlev )
      griddim6  = dim3( qsize_d , send_nelem(icycle) , 1    )
      call edgeVpack_kernel_stage<<<griddim6,blockdim6,0,streams2(icycle)>>>(edgebuf_d,array_in,putmapP_d,reverse_d,nbuf,0,1,nelemd,np1,send_indices_d,nSendCycles,icycle,tl_in)
    enddo
  endif
  do while ( nmsg_rcvd < nRecvCycles .or. nmsg_sent < nSendCycles .or. .not. internal_unpacked )
    !When this cycle's D2H memcpy is finished, call the MPI_Isend to shoot it over to the destination process
    do icycle = 1 , nSendCycles
      if ( .not. d2h_done(icycle) ) then
        if ( cudaStreamQuery(streams2(icycle)) == 0 ) then
          pCycle => pSchedule%SendCycle(icycle)
          iptr   =  pCycle%ptrP
          ierr = cudaMemcpyAsync(sendbuf_h(1,1,icycle),edgebuf_d(1,iptr),size(sendbuf_h(1:nlyr,1:pCycle%lengthP,icycle)),cudaMemcpyDeviceToHost,streams(icycle))
          d2h_done(icycle) = .true.
        endif
      endif
      if ( .not. msg_sent(icycle) ) then  !Only send once per cycle
        if ( d2h_done(icycle) ) then
          if ( cudaStreamQuery(streams(icycle)) == 0 ) then
            pCycle => pSchedule%SendCycle(icycle)
            dest   =  pCycle%dest - 1
            iptr   =  pCycle%ptrP
            length =  nlyr * pCycle%lengthP
            tag    =  pCycle%tag
            call MPI_Isend(sendbuf_h(1,1,icycle),length,MPIreal_t,dest,tag,hybrid%par%comm,Srequest(icycle),ierr)
            msg_sent(icycle) = .true.
            nmsg_sent = nmsg_sent + 1
          endif
        endif
      endif
    enddo
    if (.not. internal_unpacked) then
      if (nmsg_sent == nSendCycles) then
        blockdim6 = dim3( np      , np                  , nlev )
        griddim6  = dim3( qsize_d , recv_internal_nelem , 1    )
        call edgeVunpack_kernel_stage<<<griddim6,blockdim6>>>(edgebuf_d,array_in,getmapP_d,nbuf,0,1,nelemd,np1,recv_internal_indices_d,tl_in)
        internal_unpacked = .true.
      endif
    endif
    !When this cycle's MPI transfer is compliete, then call the D2H memcopy asynchronously
    do icycle = 1 , nRecvCycles
      if ( .not. h2d_done(icycle) ) then  !Only host to device once per cycle
        call MPI_Test(Rrequest(icycle),recvflag,status,ierr)
        if ( (ierr==MPI_SUCCESS) .and. recvflag ) then
          pCycle => pSchedule%RecvCycle(icycle)
          iptr   =  pCycle%ptrP
          ierr = cudaMemcpyAsync(recvbuf_d(1,1,icycle),recvbuf_h(1,1,icycle),size(recvbuf_h(1:nlyr,1:pCycle%lengthP,icycle)),cudaMemcpyHostToDevice,streams(icycle))
          h2d_done(icycle) = .true.
          nmsg_rcvd = nmsg_rcvd + 1 !This is how we close the polling loop, once every message has been received
        endif
      endif
    enddo
  enddo
  call MPI_WaitAll(nSendCycles,Srequest,status,ierr)
  do icycle = 1 , nRecvCycles
    pCycle => pSchedule%RecvCycle(icycle)
    iptr   =  pCycle%ptrP
    ierr = cudaMemcpyAsync(edgebuf_d(1,iptr),recvbuf_d(1,1,icycle),size(recvbuf_h(1:nlyr,1:pCycle%lengthP,icycle)),cudaMemcpyDeviceToDevice,streams(icycle))
  enddo
  ierr = cudaThreadSynchronize()

  ! Do not execute the following if recv_external_nelem = 0
  if (recv_external_nelem > 0) then
    blockdim6 = dim3( np      , np                  , nlev )
    griddim6  = dim3( qsize_d , recv_external_nelem , 1    )
    call edgeVunpack_kernel_stage<<<griddim6,blockdim6>>>(edgebuf_d,array_in,getmapP_d,nbuf,0,1,nelemd,np1,recv_external_indices_d,tl_in)
    ierr = cudaThreadSynchronize()
  endif

! ierr = cudaThreadSynchronize()
! do icycle=1,nRecvCycles
!   pCycle => pSchedule%RecvCycle(icycle)
!   source  = pCycle%source - 1
!   length  = nlyr * pCycle%lengthP
!   tag     = pCycle%tag
!   call MPI_Irecv(recvbuf_h(1,1,icycle),length,MPIreal_t,source,tag,hybrid%par%comm,Rrequest(icycle),ierr)
! enddo
! blockdim6 = dim3( np      , np     , nlev )
! griddim6  = dim3( qsize_d , nelemd , 1    )
! call edgeVpack_kernel<<<griddim6,blockdim6>>>(edgebuf_d,array_in,putmapP_d,reverse_d,nbuf,0,1,nelemd,np1,tl_in)
! ierr = cudaThreadSynchronize()
! do icycle = 1 , nSendCycles
!   pCycle => pSchedule%SendCycle(icycle)
!   iptr   =  pCycle%ptrP
!   ierr = cudaMemcpyAsync(sendbuf_h(1,1,icycle),edgebuf_d(1,iptr),size(sendbuf_h(1:nlyr,1:pCycle%lengthP,icycle)),cudaMemcpyDeviceToHost,streams(icycle))
! enddo
! ierr = cudaThreadSynchronize()
! do icycle = 1 , nSendCycles
!   pCycle => pSchedule%SendCycle(icycle)
!   dest   =  pCycle%dest - 1
!   iptr   =  pCycle%ptrP
!   length =  nlyr * pCycle%lengthP
!   tag    =  pCycle%tag
!   call MPI_Isend(sendbuf_h(1,1,icycle),length,MPIreal_t,dest,tag,hybrid%par%comm,Srequest(icycle),ierr)
! enddo
! call MPI_WaitAll(nRecvCycles,Rrequest,status,ierr)
! call MPI_WaitAll(nSendCycles,Srequest,status,ierr)
! !When this cycle's MPI transfer is compliete, then call the D2H memcopy asynchronously
! do icycle = 1 , nRecvCycles
!   pCycle => pSchedule%RecvCycle(icycle)
!   iptr   =  pCycle%ptrP
!   ierr = cudaMemcpyAsync(edgebuf_d(1,iptr),recvbuf_h(1,1,icycle),size(recvbuf_h(1:nlyr,1:pCycle%lengthP,icycle)),cudaMemcpyHostToDevice,streams(icycle))
! enddo
! blockdim6 = dim3( np      , np     , nlev )
! griddim6  = dim3( qsize_d , nelemd , 1    )
! ierr = cudaThreadSynchronize()
! call edgeVunpack_kernel<<<griddim6,blockdim6>>>(edgebuf_d,array_in,getmapP_d,nbuf,0,1,nelemd,np1,tl_in)
! ierr = cudaThreadSynchronize()

end subroutine pack_exchange_unpack_stage



attributes(global) subroutine edgeVpack_kernel_stage(edgebuf,v,putmapP,reverse,nbuf,kptr,nets,nete,nt,send_indices,nSendCycles,icycle,tl_in)
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  implicit none
  real (kind=real_kind), intent(  out) :: edgebuf(nlev*qsize_d,nbuf)
  integer              , intent(in   ) :: putmapP(max_neigh_edges_d,nets:nete)
  logical              , intent(in   ) :: reverse(max_neigh_edges_d,nets:nete)
  real (kind=real_kind), intent(in   ) :: v(np*np,nlev,qsize_d,tl_in,nets:nete)
  integer              , intent(in   ) :: send_indices(nets:nete,nSendCycles)
  integer, value       , intent(in   ) :: kptr,nets,nete,nt,nbuf,nSendCycles,icycle,tl_in
  integer :: i,j,k,q,l,offset,ij,ijk,ti,tj,tk,x,y,ir,  reverse_south, reverse_north, reverse_west, reverse_east, el
  integer, shared :: ic(max_corner_elem_d,4), direction(4), reverse_direction(4)
  real (kind=real_kind), shared :: vshrd(nlev+PAD,np,np)
  i  = threadidx%x
  j  = threadidx%y
  k  = threadidx%z
  q  = blockidx%x
  el = send_indices(blockidx%y,icycle)
  
  ij = (j-1)*np+i
  ijk = (k-1)*np*np + (j-1)*np + i -1
  tk = mod( ijk, nlev ) + 1
  ti = mod( ijk/nlev, np ) + 1
  tj = ijk/(nlev*np) + 1

  if( i+j+k == blockdim%x + blockdim%y + blockdim%z ) then
    direction(west_px)  = putmapP(west ,el)
    direction(east_px)  = putmapP(east ,el)
    direction(south_px) = putmapP(south,el)
    direction(north_px) = putmapP(north,el)
    reverse_direction(south_px) = reverse(south,el)
    reverse_direction(north_px) = reverse(north,el)
    reverse_direction(west_px)  = reverse(west,el)
    reverse_direction(east_px)  = reverse(east,el)
  endif
  if( ijk < max_corner_elem_d ) then
    ic(ijk+1,1) = putmapP(swest+ijk,el)+1
    ic(ijk+1,2) = putmapP(seast+ijk,el)+1
    ic(ijk+1,3) = putmapP(nwest+ijk,el)+1
    ic(ijk+1,4) = putmapP(neast+ijk,el)+1
  endif
  vshrd(k,i,j) = v(ij,k,q,nt,el)
  
  call syncthreads()
   
  offset = (q-1)*nlev + tk + kptr
  ir = np-ti+1
  if( 1==tj .or. 4==tj ) then
    if( reverse_direction(tj) ) then; edgebuf(offset,direction(tj)+ti) = vshrd(tk,ir,tj)
    else                            ; edgebuf(offset,direction(tj)+ti) = vshrd(tk,ti,tj); endif
  endif
  if( 2==tj ) then
    if( reverse_direction(2) ) then; edgebuf(offset,direction(tj)+ti) = vshrd(tk,1,ir)
    else                           ; edgebuf(offset,direction(tj)+ti) = vshrd(tk,1,ti); endif
  endif
  if( 3==tj ) then
    if( reverse_direction(3) ) then; edgebuf(offset,direction(tj)+ti) = vshrd(tk,4,ir)
    else                           ; edgebuf(offset,direction(tj)+ti) = vshrd(tk,4,ti); endif
  endif
  if( tj==1 ) then
    do l=1, max_corner_elem_d       
      x = mod(ti-1,2)*(np-1) + 1  ! we need to convert ti index from {1,2,3,4} to {(1,1),(4,1),(1,4),(4,4)}
      y = ((ti-1)/2)*(np-1) + 1   !   so, ti->(x,y)
      if( ic(l,ti) /= 0 ) edgebuf( offset, ic(l,ti) ) = vshrd(tk,x,y)
    enddo
  endif
end subroutine edgeVpack_kernel_stage



attributes(global) subroutine edgeVunpack_kernel_stage(edgebuf,v,getmapP,nbuf,kptr,nets,nete,nt,recv_indices,tl_in)
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  implicit none
  real (kind=real_kind), intent(in   ) :: edgebuf(nlev*qsize_d,nbuf)
  integer              , intent(in   ) :: getmapP(max_neigh_edges_d,nets:nete)
  real (kind=real_kind), intent(inout) :: v(np*np,nlev,qsize_d,tl_in,nets:nete)
  integer              , intent(in   ) :: recv_indices(nets:nete)
  integer, value       , intent(in   ) :: kptr,nets,nete,nt,nbuf,tl_in
  integer :: i,j,k,l,q,el,offset,ij,ti,tj,tk,ijk,x,y,tij
  integer, shared :: direction(4),is,ie,in,iw, ic(max_corner_elem_d,4)
  real (kind=real_kind), shared :: vshrd(nlev+PAD,np,np)
  real (kind=real_kind) :: v_before_update, neighbor_value
  
  i  = threadidx%x
  j  = threadidx%y
  k  = threadidx%z
  q  = blockidx%x
  el = recv_indices(blockidx%y)
  
  ij = (j-1)*np+i
  ijk = (k-1)*np*np + ij -1
  tk = mod( ijk, nlev ) + 1
  ti = mod( ijk/nlev, np ) + 1
  tj = ijk/(nlev*np) + 1
  
  if( i + j + k == np+np+nlev ) then
    direction(west_px)  = getmapP(west ,el)
    direction(east_px)  = getmapP(east ,el)
    direction(south_px) = getmapP(south,el)
    direction(north_px) = getmapP(north,el)
  endif
  if( ijk < max_corner_elem_d) then
    ic(ijk+1,1) = getmapP(swest+ijk,el)+1
    ic(ijk+1,2) = getmapP(seast+ijk,el)+1
    ic(ijk+1,3) = getmapP(nwest+ijk,el)+1
    ic(ijk+1,4) = getmapP(neast+ijk,el)+1
  endif
  vshrd(k,i,j) = 0.D0
  call syncthreads()
    
  offset = (q-1)*nlev + tk + kptr
  neighbor_value = edgebuf( offset, direction(tj)+ti ) ! load neighbor values into registers 
                                                       !  nlev x np consecutive threads contain all the face values
                                                       !   tj = 1:  south   
                                                       !   tj = 2:  west
                                                       !   tj = 3:  east
                                                       !   tj = 4:  north
                                                       
  ! combine the neighbor values in smem
  if( 1==tj .or. 4==tj ) vshrd(tk, ti, tj) = neighbor_value  ! add the south and north values to smem
  call syncthreads() ! this sync is needed to avoid race conditions (east/west share corners with sourth/north)
  if( 2==tj ) vshrd(tk,1,ti) = vshrd(tk,1,ti) + neighbor_value  ! update west
  if( 3==tj ) vshrd(tk,4,ti) = vshrd(tk,4,ti) + neighbor_value  ! update east
  call syncthreads()

  v_before_update = v(ij,k,q,nt,el) ! start loading the local value to be updated with neibhbor values
  
  ! update the "corner" columns
  if( tj==1 ) then
    do l=1, max_corner_elem_d       
      x = mod(ti-1,2)*(np-1) + 1  ! we need to convert ti index from {1,2,3,4} to {(1,1),(4,1),(1,4),(4,4)}
      y = ((ti-1)/2)*(np-1) + 1   !   so, ti->(x,y)
      if( ic(l,ti) /= 0 ) vshrd(tk,x,y) = vshrd(tk,x,y) + edgebuf( offset, ic(l,ti) )
    enddo
  endif
  call syncthreads()
    
  v(ij,k,q,nt,el) = v_before_update + vshrd(k,i,j)
end subroutine edgeVunpack_kernel_stage



attributes(global) subroutine edgeVpack_kernel(edgebuf,v,putmapP,reverse,nbuf,kptr,nets,nete,nt,tl_in)
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  implicit none
  real (kind=real_kind), intent(  out) :: edgebuf(nlev*qsize_d,nbuf)
  integer              , intent(in   ) :: putmapP(max_neigh_edges_d,nets:nete)
  logical              , intent(in   ) :: reverse(max_neigh_edges_d,nets:nete)
  real (kind=real_kind), intent(in   ) :: v(np*np,nlev,qsize_d,tl_in,nets:nete)
  integer, value       , intent(in   ) :: kptr,nets,nete,nt,nbuf,tl_in
  integer :: i,j,k,q,l,offset,ij,ijk,ti,tj,tk,x,y,ir,  reverse_south, reverse_north, reverse_west, reverse_east, el
  integer, shared :: ic(max_corner_elem_d,4), direction(4), reverse_direction(4)
  real (kind=real_kind), shared :: vshrd(nlev+PAD,np,np)
  i  = threadidx%x
  j  = threadidx%y
  k  = threadidx%z
  q  = blockidx%x
  el = blockidx%y
  
  ij = (j-1)*np+i
  ijk = (k-1)*np*np + (j-1)*np + i -1
  tk = mod( ijk, nlev ) + 1
  ti = mod( ijk/nlev, np ) + 1
  tj = ijk/(nlev*np) + 1

  if( i+j+k == blockdim%x + blockdim%y + blockdim%z ) then
    direction(west_px)  = putmapP(west ,el)
    direction(east_px)  = putmapP(east ,el)
    direction(south_px) = putmapP(south,el)
    direction(north_px) = putmapP(north,el)
    reverse_direction(south_px) = reverse(south,el)
    reverse_direction(north_px) = reverse(north,el)
    reverse_direction(west_px)  = reverse(west,el)
    reverse_direction(east_px)  = reverse(east,el)
  endif
  if( ijk < max_corner_elem_d ) then
    ic(ijk+1,1) = putmapP(swest+ijk,el)+1
    ic(ijk+1,2) = putmapP(seast+ijk,el)+1
    ic(ijk+1,3) = putmapP(nwest+ijk,el)+1
    ic(ijk+1,4) = putmapP(neast+ijk,el)+1
  endif
  vshrd(k,i,j) = v(ij,k,q,nt,el)
  
  call syncthreads()
   
  offset = (q-1)*nlev + tk + kptr
  ir = np-ti+1
  if( 1==tj .or. 4==tj ) then
    if( reverse_direction(tj) ) then; edgebuf(offset,direction(tj)+ti) = vshrd(tk,ir,tj)
    else                            ; edgebuf(offset,direction(tj)+ti) = vshrd(tk,ti,tj); endif
  endif
  if( 2==tj ) then
    if( reverse_direction(2) ) then; edgebuf(offset,direction(tj)+ti) = vshrd(tk,1,ir)
    else                           ; edgebuf(offset,direction(tj)+ti) = vshrd(tk,1,ti); endif
  endif
  if( 3==tj ) then
    if( reverse_direction(3) ) then; edgebuf(offset,direction(tj)+ti) = vshrd(tk,4,ir)
    else                           ; edgebuf(offset,direction(tj)+ti) = vshrd(tk,4,ti); endif
  endif
  if( tj==1 ) then
    do l=1, max_corner_elem_d       
      x = mod(ti-1,2)*(np-1) + 1  ! we need to convert ti index from {1,2,3,4} to {(1,1),(4,1),(1,4),(4,4)}
      y = ((ti-1)/2)*(np-1) + 1   !   so, ti->(x,y)
      if( ic(l,ti) /= 0 ) edgebuf( offset, ic(l,ti) ) = vshrd(tk,x,y)
    enddo
  endif
end subroutine edgeVpack_kernel



attributes(global) subroutine edgeVunpack_kernel(edgebuf,v,getmapP,nbuf,kptr,nets,nete,nt,tl_in)
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  real (kind=real_kind), intent(in   ) :: edgebuf(nlev*qsize_d,nbuf)
  integer              , intent(in   ) :: getmapP(max_neigh_edges_d,nets:nete)
  real (kind=real_kind), intent(inout) :: v(np*np,nlev,qsize_d,tl_in,nets:nete)
  integer, value       , intent(in   ) :: kptr,nets,nete,nt,nbuf,tl_in
  integer :: i,j,k,l,q,el,offset,ij,ti,tj,tk,ijk,x,y,tij
  integer, shared :: direction(4),is,ie,in,iw, ic(max_corner_elem_d,4)
  real (kind=real_kind), shared :: vshrd(nlev+PAD,np,np)
  real (kind=real_kind) :: v_before_update, neighbor_value
  
  i  = threadidx%x
  j  = threadidx%y
  k  = threadidx%z
  q  = blockidx%x
  el = blockidx%y
  
  ij = (j-1)*np+i
  ijk = (k-1)*np*np + ij -1
  tk = mod( ijk, nlev ) + 1
  ti = mod( ijk/nlev, np ) + 1
  tj = ijk/(nlev*np) + 1
  
  if( i + j + k == np+np+nlev ) then
    direction(west_px)  = getmapP(west ,el)
    direction(east_px)  = getmapP(east ,el)
    direction(south_px) = getmapP(south,el)
    direction(north_px) = getmapP(north,el)
  endif
  if( ijk < max_corner_elem_d) then
    ic(ijk+1,1) = getmapP(swest+ijk,el)+1
    ic(ijk+1,2) = getmapP(seast+ijk,el)+1
    ic(ijk+1,3) = getmapP(nwest+ijk,el)+1
    ic(ijk+1,4) = getmapP(neast+ijk,el)+1
  endif
  vshrd(k,i,j) = 0.D0
  call syncthreads()
    
  offset = (q-1)*nlev + tk + kptr
  neighbor_value = edgebuf( offset, direction(tj)+ti ) ! load neighbor values into registers 
                                                       !  nlev x np consecutive threads contain all the face values
                                                       !   tj = 1:  south   
                                                       !   tj = 2:  west
                                                       !   tj = 3:  east
                                                       !   tj = 4:  north
                                                       
  ! combine the neighbor values in smem
  if( 1==tj .or. 4==tj ) vshrd(tk, ti, tj) = neighbor_value  ! add the south and north values to smem
  call syncthreads() ! this sync is needed to avoid race conditions (east/west share corners with sourth/north)
  if( 2==tj ) vshrd(tk,1,ti) = vshrd(tk,1,ti) + neighbor_value  ! update west
  if( 3==tj ) vshrd(tk,4,ti) = vshrd(tk,4,ti) + neighbor_value  ! update east
  call syncthreads()

  v_before_update = v(ij,k,q,nt,el) ! start loading the local value to be updated with neibhbor values
  
  ! update the "corner" columns
  if( tj==1 ) then
    do l=1, max_corner_elem_d       
      x = mod(ti-1,2)*(np-1) + 1  ! we need to convert ti index from {1,2,3,4} to {(1,1),(4,1),(1,4),(4,4)}
      y = ((ti-1)/2)*(np-1) + 1   !   so, ti->(x,y)
      if( ic(l,ti) /= 0 ) vshrd(tk,x,y) = vshrd(tk,x,y) + edgebuf( offset, ic(l,ti) )
    enddo
  endif
  call syncthreads()
    
  v(ij,k,q,nt,el) = v_before_update + vshrd(k,i,j)
end subroutine edgeVunpack_kernel



attributes(global) subroutine hypervis_kernel1( Qdp , qtens , dp , dinv , variable_hyperviscosity , spheremp , deriv_dvv , nets , nete , dt , nt )
  implicit none
  real(kind=real_kind), dimension(np,np,nlev,qsize_d,timelevels,nets:nete), intent(in   ) :: Qdp
  real(kind=real_kind), dimension(np,np,nlev,qsize_d           ,nets:nete), intent(  out) :: Qtens
  real(kind=real_kind), dimension(np,np,nlev                   ,nets:nete), intent(in   ) :: dp
  real(kind=real_kind), dimension(np,np,4                      ,nets:nete), intent(in   ) :: dinv
  real(kind=real_kind), dimension(np,np                                  ), intent(in   ) :: deriv_dvv
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: variable_hyperviscosity
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: spheremp
  real(kind=real_kind), value                                             , intent(in   ) :: dt
  integer, value                                                          , intent(in   ) :: nets , nete , nt
  real(kind=real_kind), dimension(np,np,2,nlev), shared :: s
  real(kind=real_kind), dimension(np,np,4  ), shared :: dinv_s
  real(kind=real_kind), dimension(np,np), shared :: variable_hyperviscosity_s
  real(kind=real_kind), dimension(np,np), shared :: spheremp_s
  integer :: i, j, k, q, ie, iz
  i  = threadidx%x
  j  = threadidx%y
  k  = threadidx%z
  ie = blockidx%x
  iz = threadidx%z
  if (k <= 4) dinv_s(i,j,k) = dinv(i,j,k,ie)
  if (k == 1) then
    variable_hyperviscosity_s(i,j) = variable_hyperviscosity(i,j,ie)
    spheremp_s(i,j) = spheremp(i,j,ie)
  endif
  do q = 1 , qsize_d
    s(i,j,1,iz) = Qdp(i,j,k,q,nt,ie) / dp(i,j,k,ie)
    call syncthreads()
    qtens(i,j,k,q,ie) = laplace_sphere_wk(i,j,ie,iz,s,dinv_s,spheremp_s,variable_hyperviscosity_s,deriv_dvv,nets,nete)
  enddo
end subroutine hypervis_kernel1



attributes(global) subroutine hypervis_kernel2( Qdp , qtens , dp , dinv , variable_hyperviscosity , spheremp , rspheremp , deriv_dvv , hyai , hybi , ps0 , nu_q , nets , nete , dt , nt )
  implicit none
  real(kind=real_kind), dimension(np,np,nlev,qsize_d,timelevels,nets:nete), intent(inout) :: Qdp
  real(kind=real_kind), dimension(np,np,nlev,qsize_d           ,nets:nete), intent(in   ) :: Qtens
  real(kind=real_kind), dimension(np,np,nlev                   ,nets:nete), intent(in   ) :: dp
  real(kind=real_kind), dimension(np,np                                  ), intent(in   ) :: deriv_dvv
  real(kind=real_kind), dimension(      nlev+1                           ), intent(in   ) :: hyai
  real(kind=real_kind), dimension(      nlev+1                           ), intent(in   ) :: hybi
  real(kind=real_kind), dimension(np,np,4                      ,nets:nete), intent(in   ) :: dinv
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: variable_hyperviscosity
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: spheremp
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: rspheremp
  real(kind=real_kind), value                                             , intent(in   ) :: dt, ps0, nu_q
  integer, value                                                          , intent(in   ) :: nets , nete , nt
  real(kind=real_kind), dimension(np,np,2,nlev), shared :: s
  real(kind=real_kind), dimension(np,np,4  ), shared :: dinv_s
  real(kind=real_kind), dimension(np,np), shared :: variable_hyperviscosity_s
  real(kind=real_kind), dimension(np,np), shared :: spheremp_s
  real(kind=real_kind), dimension(np,np), shared :: rspheremp_s
  real(kind=real_kind) :: dp0
  integer :: i, j, k, q, ie, iz
  i  = threadidx%x
  j  = threadidx%y
  k  = threadidx%z
  ie = blockidx%x
  iz = threadidx%z
  if (k == 1) then
    dinv_s(i,j,:) = dinv(i,j,:,ie)
    variable_hyperviscosity_s(i,j) = variable_hyperviscosity(i,j,ie)
    spheremp_s(i,j) = spheremp(i,j,ie)
    rspheremp_s(i,j) = rspheremp(i,j,ie)
  endif
  call syncthreads()
  dp0 = dt*nu_q*( ( hyai(k+1) - hyai(k) )*ps0 + ( hybi(k+1) - hybi(k) )*ps0 )
  do q = 1 , qsize_d
    s(i,j,1,iz) = rspheremp_s(i,j)*qtens(i,j,k,q,ie)
    call syncthreads()
    Qdp(i,j,k,q,nt,ie) = Qdp(i,j,k,q,nt,ie)*spheremp_s(i,j)-dp0*laplace_sphere_wk(i,j,ie,iz,s,dinv_s,spheremp_s,variable_hyperviscosity_s,deriv_dvv,nets,nete)
  enddo
end subroutine hypervis_kernel2



attributes(device) function laplace_sphere_wk(i,j,ie,iz,s,dinv,spheremp,variable_hyperviscosity,deriv_dvv,nets,nete) result(lapl)
  implicit none
  integer,                                              intent(in) :: nets, nete, i, j, ie, iz
  real(kind=real_kind), dimension(np,np,2,nlev)       , intent(inout) :: s
  real(kind=real_kind), dimension(np,np,2,2          ), intent(in) :: dinv
  real(kind=real_kind), dimension(np,np              ), intent(in) :: deriv_dvv
  real(kind=real_kind), dimension(np,np              ), intent(in) :: variable_hyperviscosity
  real(kind=real_kind), dimension(np,np              ), intent(in) :: spheremp
  real(kind=real_kind)                                             :: lapl
  integer :: l
  real(kind=real_kind) :: dsdx00 , dsdy00, tmp1, tmp2
  real(kind=real_kind), dimension(2) :: ds
  ds = gradient_sphere(i,j,ie,iz,s,dinv,variable_hyperviscosity,deriv_dvv,nets,nete)
  lapl = divergence_sphere_wk(i,j,ie,iz,ds,s,dinv,spheremp,deriv_dvv,nets,nete)
end function laplace_sphere_wk



attributes(device) function gradient_sphere(i,j,ie,iz,s,dinv,variable_hyperviscosity,deriv_dvv,nets,nete) result(ds)
  use physical_constants, only: rrearth
  implicit none
  integer,                                              intent(in) :: nets, nete, i, j, ie, iz
  real(kind=real_kind), dimension(np,np,2,nlev)       , intent(in) :: s
  real(kind=real_kind), dimension(np,np,2,2          ), intent(in) :: dinv
  real(kind=real_kind), dimension(np,np              ), intent(in) :: deriv_dvv
  real(kind=real_kind), dimension(np,np              ), intent(in) :: variable_hyperviscosity
  real(kind=real_kind), dimension(2)                               :: ds
  integer :: l
  real(kind=real_kind) :: dsdx00 , dsdy00, tmp1, tmp2
  dsdx00 = 0.0d0
  dsdy00 = 0.0d0
  do l = 1 , np
    dsdx00 = dsdx00 + deriv_dvv(l,i)*s(l,j,1,iz)
    dsdy00 = dsdy00 + deriv_dvv(l,j)*s(i,l,1,iz)
  enddo
  ds(1) = ( dinv(i,j,1,1)*dsdx00 + dinv(i,j,2,1)*dsdy00 ) * rrearth * variable_hyperviscosity(i,j)
  ds(2) = ( dinv(i,j,1,2)*dsdx00 + dinv(i,j,2,2)*dsdy00 ) * rrearth * variable_hyperviscosity(i,j)
end function gradient_sphere



attributes(device) function divergence_sphere_wk(i,j,ie,iz,tmp,s,dinv,spheremp,deriv_dvv,nets,nete) result(lapl)
  use physical_constants, only: rrearth
  implicit none
  integer,                                              intent(in   ) :: nets, nete, i, j, ie, iz
  real(kind=real_kind), dimension(2),                   intent(in   ) :: tmp
  real(kind=real_kind), dimension(np,np,2,nlev)       , intent(inout) :: s
  real(kind=real_kind), dimension(np,np,2,2          ), intent(in   ) :: dinv
  real(kind=real_kind), dimension(np,np              ), intent(in   ) :: deriv_dvv
  real(kind=real_kind), dimension(np,np              ), intent(in   ) :: spheremp
  real(kind=real_kind)                                                :: lapl
  integer :: l
  s(i,j,1,iz) = ( dinv(i,j,1,1)*tmp(1) + dinv(i,j,1,2)*tmp(2) )
  s(i,j,2,iz) = ( dinv(i,j,2,1)*tmp(1) + dinv(i,j,2,2)*tmp(2) )
  call syncthreads()
  lapl = 0.0d0
  do l = 1 , np
    lapl = lapl - (spheremp(l,j)*s(l,j,1,iz)*deriv_dvv(i,l)+spheremp(i,l)*s(i,l,2,iz)*deriv_dvv(j,l)) * rrearth
  enddo
end function divergence_sphere_wk



#endif
end module cuda_mod


