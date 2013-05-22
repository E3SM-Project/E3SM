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
  public :: vertical_remap_cuda
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
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:)     :: dp_star_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: qmin_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: qmax_d
  integer              ,device,allocatable,dimension(:)           :: send_nelem_d
  integer              ,device,allocatable,dimension(:)           :: recv_nelem_d
  integer              ,device,allocatable,dimension(:,:)         :: send_indices_d
  integer              ,device,allocatable,dimension(:,:)         :: recv_indices_d
  integer              ,device,allocatable,dimension(:)           :: recv_internal_indices_d
  integer              ,device,allocatable,dimension(:)           :: recv_external_indices_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: recvbuf_d
  real(kind=real_kind) ,device,allocatable,dimension(:,:,:,:)     :: dpo_d
  real(kind=real_kind) ,device,allocatable,dimension(:,:,:,:,:)   :: ppmdx_d
  real(kind=real_kind) ,device,allocatable,dimension(:,:,:,:)     :: z1_d
  real(kind=real_kind) ,device,allocatable,dimension(:,:,:,:)     :: z2_d
  integer              ,device,allocatable,dimension(:,:,:,:)     :: kid_d

  !PINNED Host arrays
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:,:,:) :: qdp_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:,:)   :: vstar_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:,:)   :: qtens_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:)     :: dp_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:)     :: divdp_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:)     :: dp_star_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:)     :: dp_np1_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:)       :: sendbuf_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:)       :: recvbuf_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:,:)   :: qtens_biharmonic
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:)       :: qmin
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:)       :: qmax
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:)     :: dpdiss_biharmonic_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:)     :: dpo_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:,:)   :: ppmdx_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:)     :: z1_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:)     :: z2_h
  integer             ,pinned,allocatable,dimension(:,:,:,:)     :: kid_h

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
  integer :: cuda_streams
  integer(kind=cuda_stream_kind), allocatable :: streams(:)
  integer(kind=cuda_stream_kind), allocatable :: streams2(:)
  integer            :: nbuf
  integer            :: nmsg_rcvd
  integer            :: nmsg_sent
  integer, device :: max_neigh_edges_d, max_corner_elem_d

  type(cudaEvent) :: timer1, timer2
  integer, parameter :: gs = 2


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

    cuda_streams = nelemd
    allocate(streams (0:cuda_streams))
    allocate(streams2(0:cuda_streams))

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
    allocate( dp_star_d                (np,np,nlev                   ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( reverse_d                (max_neigh_edges              ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( putmapP_d                (max_neigh_edges              ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( getmapP_d                (max_neigh_edges              ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( edgebuf_d                (nlev*qsize_d,nbuf                   ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( recv_internal_indices_d  (nelemd                              ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( recv_external_indices_d  (nelemd                              ) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( dpo_d                    (np,np,   1-gs:nlev+gs        ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( ppmdx_d                  (np,np,10,   0:nlev+1         ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( z1_d                     (np,np,        nlev           ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( z2_d                     (np,np,        nlev           ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( kid_d                    (np,np,        nlev           ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__

    allocate( qdp_h                    (np,np,nlev,qsize_d,timelevels,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( qmin                     (nlev,qsize_d                 ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( qmax                     (nlev,qsize_d                 ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( vstar_h                  (np,np,nlev,2                 ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( qtens_h                  (np,np,nlev,qsize_d           ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( qtens_biharmonic         (np,np,nlev,qsize_d           ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( dp_h                     (np,np,nlev                   ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( divdp_h                  (np,np,nlev                   ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( dp_star_h                (np,np,nlev                   ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( dpdiss_biharmonic_h      (np,np,nlev                   ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( dp_np1_h                 (np,np,nlev                   ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( dpo_h                    (np,np,   1-gs:nlev+gs        ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( ppmdx_h                  (np,np,10,   0:nlev+1         ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( z1_h                     (np,np,        nlev           ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( z2_h                     (np,np,        nlev           ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__
    allocate( kid_h                    (np,np,        nlev           ,nelemd) , stat = ierr ) ; if ( ierr .ne. 0 ) stop __LINE__

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
!$OMP BARRIER
!$OMP MASTER
    ierr = cudaThreadSynchronize()
    do ie = 1,nelemd
      ierr = cudaMemcpyAsync( qdp_d(1,1,1,1,nt,ie) , elem(ie)%state%qdp(1,1,1,1,nt) , size(elem(ie)%state%qdp(:,:,:,:,nt)) , cudaMemcpyHostToDevice , streams(1) )
    enddo
    ierr = cudaThreadSynchronize()
!$OMP END MASTER
!$OMP BARRIER
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
  use parallel_mod      , only: iam
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
!   call copy_qdp_h2d( elem , n0_qdp)
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
  use parallel_mod   , only : iam
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
      x = mod(ti-1,2)*(np-1) + 1  ! we need to convert ti index from {1,2,3,4} to {(1,1),(4,1),(1,4),(4,4)}€;€;€;€;
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



subroutine vertical_remap_cuda(elem,fvm,hvcoord,dt,np1,np1_qdp,nets,nete)
  use kinds, only : real_kind
  use hybvcoord_mod, only : hvcoord_t
  use control_mod, only :  rsplit
  use parallel_mod, only : abortmp, iam
  use element_mod, only: element_t
  use dimensions_mod, only: nc, ntrac
  use perf_mod, only: t_startf, t_stopf
#if defined(_SPELT)
  use spelt_mod, only: spelt_struct
#else
  use fvm_control_volume_mod, only : fvm_struct
#endif    
#if defined(_SPELT)
  type(spelt_struct), intent(inout) :: fvm(:)
  real (kind=real_kind) :: cdp(1:nep,1:nep,nlev,ntrac-1) 
  real (kind=real_kind)  :: psc(nep,nep), dpc(nep,nep,nlev),dpc_star(nep,nep,nlev)
#else
  type(fvm_struct), intent(inout) :: fvm(:)
  real (kind=real_kind) :: cdp(1:nc,1:nc,nlev,ntrac-1) 
  real (kind=real_kind)  :: psc(nc,nc), dpc(nc,nc,nlev),dpc_star(nc,nc,nlev)
#endif
  type (element_t), intent(inout)   :: elem(:)
  type (hvcoord_t)                  :: hvcoord
  real (kind=real_kind)             :: dt
  integer :: ie,i,j,k,np1,nets,nete,np1_qdp,ierr
  real (kind=real_kind), dimension(np,np,nlev,2)  :: ttmp
  call t_startf('vertical_remap_cuda')
  do ie=nets,nete
!    ! SET VERTICAL VELOCITY TO ZERO FOR DEBUGGING
!    elem(ie)%derived%eta_dot_dpdn(:,:,:)=0
    if (rsplit==0) then
      ! compute dp_star from eta_dot_dpdn():
      do k=1,nlev
        dp_h(:,:,k,ie) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
        dp_star_h(:,:,k,ie) = dp_h(:,:,k,ie) + dt*(elem(ie)%derived%eta_dot_dpdn(:,:,k+1) - elem(ie)%derived%eta_dot_dpdn(:,:,k)) 
      enddo
      if (minval(dp_star_h(:,:,:,ie))<0) call abortmp('negative layer thickness.  timestep or remap time too large')
    else
      !  REMAP u,v,T from levels in dp3d() to REF levels
      ! update final ps_v 
      elem(ie)%state%ps_v(:,:,np1) = hvcoord%hyai(1)*hvcoord%ps0 + sum(elem(ie)%state%dp3d(:,:,:,np1),3)
      do k=1,nlev
        dp_h(:,:,k,ie) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
        dp_star_h(:,:,k,ie) = elem(ie)%state%dp3d(:,:,k,np1)
      enddo
      if (minval(dp_star_h(:,:,:,ie))<0) call abortmp('negative layer thickness.  timestep or remap time too large')
      ! remap the dynamics:  
#undef REMAP_TE
#ifdef REMAP_TE
      ! remap u,v and cp*T + .5 u^2 
      ttmp(:,:,:,1)=(elem(ie)%state%v(:,:,1,:,np1)**2 + elem(ie)%state%v(:,:,2,:,np1)**2)/2 + elem(ie)%state%t(:,:,:,np1)*cp
#else
      ttmp(:,:,:,1)=elem(ie)%state%t(:,:,:,np1)
#endif
      ttmp(:,:,:,1)=ttmp(:,:,:,1)*dp_star_h(:,:,:,ie)
      call remap_Q_ppm(ttmp,np,1,dp_star_h(:,:,:,ie),dp_h(:,:,:,ie))
      elem(ie)%state%t(:,:,:,np1)=ttmp(:,:,:,1)/dp_h(:,:,:,ie)
      ttmp(:,:,:,1)=elem(ie)%state%v(:,:,1,:,np1)*dp_star_h(:,:,:,ie)
      ttmp(:,:,:,2)=elem(ie)%state%v(:,:,2,:,np1)*dp_star_h(:,:,:,ie)
      call remap_Q_ppm(ttmp,np,2,dp_star_h(:,:,:,ie),dp_h(:,:,:,ie)) 
      elem(ie)%state%v(:,:,1,:,np1)=ttmp(:,:,:,1)/dp_h(:,:,:,ie)
      elem(ie)%state%v(:,:,2,:,np1)=ttmp(:,:,:,2)/dp_h(:,:,:,ie)
#ifdef REMAP_TE
      ! back out T from TE
      elem(ie)%state%t(:,:,:,np1) = ( elem(ie)%state%t(:,:,:,np1) - ( (elem(ie)%state%v(:,:,1,:,np1)**2 + elem(ie)%state%v(:,:,2,:,np1)**2)/2))/cp
#endif
    endif
  enddo
!$OMP BARRIER
!$OMP MASTER
  do ie=1,nelemd
    ! remap the tracers from lagrangian levels (dp_star)  to REF levels dp
    if (qsize>0) call remap_Q_ppm_cuda(elem,np,qsize,dp_star_h(:,:,:,ie),dp_h(:,:,:,ie),np1_qdp,ie)
  enddo
  ierr = cudaThreadSynchronize()
!$OMP END MASTER
!$OMP BARRIER
  do ie = nets , nete
    elem(ie)%state%qdp(:,:,:,:,np1_qdp) = qdp_h(:,:,:,:,np1_qdp,ie)
  enddo
! call copy_qdp_d2h( elem , np1_qdp )
  call t_stopf('vertical_remap_cuda')  
end subroutine vertical_remap_cuda



!This uses the exact same model and reference grids and data as remap_Q, but it interpolates
!using PPM instead of splines.
subroutine remap_Q_ppm_cuda(elem,nx,qsize,dp1,dp2,np1_qdp,ie)
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass
  !
  use perf_mod, only    : t_startf, t_stopf
  use control_mod, only : prescribed_wind, vert_remap_q_alg
  use element_mod, only : element_t
  implicit none
  type(element_t)      , intent(inout) :: elem(:)
  integer,intent(in) :: nx,qsize
  real (kind=real_kind), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
  integer, value       , intent(in) :: np1_qdp,ie
  ! Local Variables
  real(kind=real_kind), dimension(       nlev+2 ) :: pio    !Pressure at interfaces for old grid
  real(kind=real_kind), dimension(       nlev+1 ) :: pin    !Pressure at interfaces for new grid
  real(kind=real_kind), dimension(  1-gs:nlev+gs) :: dpn    !change in pressure over a cell for old grid
  real(kind=real_kind), dimension(nx,nx, nlev   ) :: z1, z2
  real(kind=real_kind) :: ppmdx(nx,nx,10,0:nlev+1)  !grid spacings
  integer :: i, j, k, q, kk, kid(nx,nx,nlev), ierr
  type(dim3) :: griddim, blockdim

  call t_startf('remap_Q_ppm_cuda')

  if ( vert_remap_q_alg /= 1 .and. vert_remap_q_alg /= 2 ) then
    write(*,*) 'MUST USE  VERT_REMAP_Q_ALG == 1 OR 2 WITH THE CUDA PORT.'
    stop
  endif

  do j = 1 , nx
    do i = 1 , nx
      pin(1)=0
      pio(1)=0
      do k=1,nlev
         dpn(k)=dp2(i,j,k)
         dpo_h(i,j,k,ie)=dp1(i,j,k)
         pin(k+1)=pin(k)+dpn(k)
         pio(k+1)=pio(k)+dpo_h(i,j,k,ie)
      enddo
      pio(nlev+2) = pio(nlev+1) + 1.  !This is here to allow an entire block of k threads to run in the remapping phase.
                                      !It makes sure there's an old interface value below the domain that is larger.
      pin(nlev+1) = pio(nlev+1)       !The total mass in a column does not change.
                                      !Therefore, the pressure of that mass cannot either.
      !Fill in the ghost regions with mirrored values. if vert_remap_q_alg is defined, this is of no consequence.
      do k = 1 , gs
        dpo_h(i,j,1   -k,ie) = dpo_h(i,j,       k,ie)
        dpo_h(i,j,nlev+k,ie) = dpo_h(i,j,nlev+1-k,ie)
      enddo
      !Compute remapping intervals once for all tracers. Find the old grid cell index in which the
      !k-th new cell interface resides. Then integrate from the bottom of that old cell to the new
      !interface location. In practice, the grid never deforms past one cell, so the search can be
      !simplified by this. Also, the interval of integration is usually of magnitude close to zero
      !or close to dpo because of minimial deformation.
      !Numerous tests confirmed that the bottom and top of the grids match to machine precision, so
      !I set them equal to each other.
      do k = 1 , nlev
        kk = k  !Keep from an order n^2 search operation by assuming the old cell index is close.
        !Find the index of the old grid cell in which this new cell's bottom interface resides.
        do while ( pio(kk) <= pin(k+1) )
          kk = kk + 1
        enddo
        kk = kk - 1                   !kk is now the cell index we're integrating over.
        if (kk == nlev+1) kk = nlev   !This is to keep the indices in bounds.
                                      !Top bounds match anyway, so doesn't matter what coefficients are used
        kid_h(i,j,k,ie) = kk                   !Save for reuse
        z1_h(i,j,k,ie) = -0.5D0                !This remapping assumes we're starting from the left interface of an old grid cell
                                      !In fact, we're usually integrating very little or almost all of the cell in question
        z2_h(i,j,k,ie) = ( pin(k+1) - ( pio(kk) + pio(kk+1) ) * 0.5 ) / dpo_h(i,j,kk,ie)  !PPM interpolants are normalized to an independent
                                                                        !coordinate domain [-0.5,0.5].
      enddo
      !This turned out a big optimization, remembering that only parts of the PPM algorithm depends on the data, namely the
      !limiting. So anything that depends only on the grid is pre-computed outside the tracer loop.
      ppmdx_h(i,j,:,:,ie) = compute_ppm_grids_d( dpo_h(i,j,:,ie) )
      !From here, we loop over tracers for only those portions which depend on tracer data, which includes PPM limiting and
      !mass accumulation
    enddo
  enddo
  griddim  = dim3( qsize , 1  , 1      )
  blockdim = dim3( np    , np , nlev+2 )
  ierr = cudaMemcpyAsync( dpo_d(:,:,:,ie)     , dpo_h(:,:,:,ie)     , size(dpo_h(:,:,:,ie)    ) , cudaMemcpyHostToDevice , streams(ie) )
  ierr = cudaMemcpyAsync( ppmdx_d(:,:,:,:,ie) , ppmdx_h(:,:,:,:,ie) , size(ppmdx_h(:,:,:,:,ie)) , cudaMemcpyHostToDevice , streams(ie) )
  ierr = cudaMemcpyAsync( z1_d(:,:,:,ie)      , z1_h(:,:,:,ie)      , size(z1_h(:,:,:,ie)     ) , cudaMemcpyHostToDevice , streams(ie) )
  ierr = cudaMemcpyAsync( z2_d(:,:,:,ie)      , z2_h(:,:,:,ie)      , size(z2_h(:,:,:,ie)     ) , cudaMemcpyHostToDevice , streams(ie) )
  ierr = cudaMemcpyAsync( kid_d(:,:,:,ie)     , kid_h(:,:,:,ie)     , size(kid_h(:,:,:,ie)    ) , cudaMemcpyHostToDevice , streams(ie) )
  call remap_Q_ppm_cuda_kernel<<<griddim,blockdim,0,streams(ie)>>>( Qdp_d(:,:,:,:,np1_qdp,ie) , dpo_d(:,:,:,ie) , ppmdx_d(:,:,:,:,ie) , z1_d(:,:,:,ie) , z2_d(:,:,:,ie) , kid_d(:,:,:,ie) , vert_remap_q_alg )
  ierr = cudaMemcpyAsync( qdp_h(:,:,:,:,np1_qdp,ie) , qdp_d(:,:,:,:,np1_qdp,ie) , size(elem(ie)%state%qdp(:,:,:,:,np1_qdp) ) , cudaMemcpyDeviceToHost , streams(ie) )
  call t_stopf('remap_Q_ppm_cuda')
end subroutine remap_Q_ppm_cuda



attributes(global) subroutine remap_Q_ppm_cuda_kernel( Qdp , dpo , ppmdx , z1 , z2 , kid , vert_remap_q_alg )
  implicit none
  real (kind=real_kind), intent(inout) :: Qdp   (np,np,        nlev   ,qsize_d)
  real(kind=real_kind) , intent(in   ) :: dpo   (np,np,   1-gs:nlev+gs        )
  real(kind=real_kind) , intent(in   ) :: ppmdx (np,np,10,   0:nlev+1         )
  real(kind=real_kind) , intent(in   ) :: z1    (np,np,        nlev           )
  real(kind=real_kind) , intent(in   ) :: z2    (np,np,        nlev           )
  integer              , intent(in   ) :: kid   (np,np,        nlev           )
  integer , value      , intent(in   ) :: vert_remap_q_alg
  integer ::  i , j ,k, q, kk, nx
  real(kind=real_kind), shared :: masso (np,np,       nlev+1 ) !Accumulate mass up to each interface
  real(kind=real_kind), shared :: ao    (np,np,  1-gs:nlev+gs) !Tracer value on old grid
  real(kind=real_kind), shared :: coefs (np,np,3,     nlev   ) !PPM coefficients within each cell
  real(kind=real_kind), shared :: massn1(np,np               )
  real(kind=real_kind), shared :: massn2(np,np               )
  real(kind=real_kind), shared :: ai    (np,np,0:nlev  )                     !fourth-order accurate, then limited interface values
  real(kind=real_kind), shared :: dma   (np,np,0:nlev+1)                     !An expression from Collela's '84 publication
  nx = np
  q = blockidx%x
  i = threadidx%x
  j = threadidx%y
  !Accumulate the old mass up to old grid cell interface locations to simplify integration
  !during remapping. Also, divide out the grid spacing so we're working with actual tracer
  !values and can conserve mass. The option for ifndef ZEROHORZ I believe is there to ensure
  !tracer consistency for an initially uniform field. I copied it from the old remap routine.
  masso(i,j,1) = 0.
  call syncthreads()
  if (threadidx%z == 1) then
    do k = 1 , nlev
      ao(i,j,k) = Qdp(i,j,k,q)
      masso(i,j,k+1) = masso(i,j,k) + ao(i,j,k) !Accumulate the old mass. This will simplify the remapping
      ao(i,j,k) = ao(i,j,k) / dpo(i,j,k)        !Divide out the old grid spacing because we want the tracer mixing ratio, not mass.
    enddo
  endif
  call syncthreads()
  !Fill in ghost values. Ignored if vert_remap_q_alg == 2
  if (threadidx%z == 1) then
    do k = 1 , gs
      ao(i,j,1   -k) = ao(i,j,       k)
      ao(i,j,nlev+k) = ao(i,j,nlev+1-k)
    enddo
  endif
  call syncthreads()
  !Compute monotonic and conservative PPM reconstruction over every cell
  coefs(:,:,:,:) = compute_ppm_d( ao , ppmdx , ai , dma , vert_remap_q_alg , i , j )
  call syncthreads()
  !Compute tracer values on the new grid by integrating from the old cell bottom to the new
  !cell interface to form a new grid mass accumulation. Taking the difference between
  !accumulation at successive interfaces gives the mass inside each cell. Since Qdp is
  !supposed to hold the full mass this needs no normalization.
  if (threadidx%z == 1) then
    massn1(:,:) = 0.
    do k = 1 , nlev
      kk = kid(i,j,k)
      massn2(i,j) = masso(i,j,kk) + ( coefs(i,j,1,kk) * (z2(i,j,k)      - z1(i,j,k)     )         + &
                                      coefs(i,j,2,kk) * (z2(i,j,k) ** 2 - z1(i,j,k) ** 2) / 0.2D1 + &
                                      coefs(i,j,3,kk) * (z2(i,j,k) ** 3 - z1(i,j,k) ** 3) / 0.3D1 ) * dpo(i,j,kk)
      Qdp(i,j,k,q) = massn2(i,j) - massn1(i,j)
      massn1(i,j) = massn2(i,j)
    enddo
  endif
end subroutine remap_Q_ppm_cuda_kernel



!This computes a limited parabolic interpolant using a net 5-cell stencil, but the stages of computation are broken up into 3 stages
attributes(device) function compute_ppm_d( a , dx , ai , dma , vert_remap_q_alg , ii , jj )    result(coefs)
  implicit none
  real(kind=real_kind), intent(in   ) :: a    (np,np,    -1:nlev+2)  !Cell-mean values
  real(kind=real_kind), intent(in   ) :: dx   (np,np,10,  0:nlev+1)  !grid spacings
  real(kind=real_kind), intent(inout) :: ai   (np,np,0:nlev  )                     !fourth-order accurate, then limited interface values
  real(kind=real_kind), intent(inout) :: dma  (np,np,0:nlev+1)                     !An expression from Collela's '84 publication
  integer , value     , intent(in   ) :: vert_remap_q_alg , ii , jj
  real(kind=real_kind) ::                coefs(np,np,0:2,   nlev  )  !PPM coefficients (for parabola)
  real(kind=real_kind) :: da                                !Ditto
  ! Hold expressions based on the grid (which are cumbersome).
  real(kind=real_kind) :: al, ar                            !Left and right interface values for cell-local limiting
  integer :: j, nx
  integer :: indB, indE
  nx = np
  j = threadidx%z-1

  call syncthreads()
  ! Stage 1: Compute dma for each cell, allowing a 1-cell ghost stencil below and above the domain
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-1
  else
    indB = 0
    indE = nlev+1
  endif
  if ( j >= indB .and. j <= indE ) then
    da = dx(ii,jj,1,j) * ( dx(ii,jj,2,j) * ( a(ii,jj,j+1) - a(ii,jj,j) ) + dx(ii,jj,3,j) * ( a(ii,jj,j) - a(ii,jj,j-1) ) )
    dma(ii,jj,j) = minval( (/ abs(da) , 2. * abs( a(ii,jj,j) - a(ii,jj,j-1) ) , 2. * abs( a(ii,jj,j+1) - a(ii,jj,j) ) /) ) * sign(1.D0,da)
    if ( ( a(ii,jj,j+1) - a(ii,jj,j) ) * ( a(ii,jj,j) - a(ii,jj,j-1) ) <= 0. ) dma(ii,jj,j) = 0.
  endif
  call syncthreads()

  ! Stage 2: Compute ai for each cell interface in the physical domain (dimension nlev+1)
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-2
  else
    indB = 0
    indE = nlev
  endif
  if ( j >= indB .and. j <= indE ) then
    ai(ii,jj,j) = a(ii,jj,j) + dx(ii,jj,4,j) * ( a(ii,jj,j+1) - a(ii,jj,j) ) + dx(ii,jj,5,j) * ( dx(ii,jj,6,j) * ( dx(ii,jj,7,j) - dx(ii,jj,8,j) ) &
                  * ( a(ii,jj,j+1) - a(ii,jj,j) ) - dx(ii,jj,9,j) * dma(ii,jj,j+1) + dx(ii,jj,10,j) * dma(ii,jj,j) )
  endif
  call syncthreads()

  ! Stage 3: Compute limited PPM interpolant over each cell in the physical domain
  ! (dimension nlev) using ai on either side and ao within the cell.
  if (vert_remap_q_alg == 2) then
    indB = 3
    indE = nlev-2
  else
    indB = 1
    indE = nlev
  endif
  if ( j >= indB .and. j <= indE ) then
    al = ai(ii,jj,j-1)
    ar = ai(ii,jj,j  )
    if ( (ar - a(ii,jj,j)) * (a(ii,jj,j) - al) <= 0. ) then
      al = a(ii,jj,j)
      ar = a(ii,jj,j)
    endif
    if ( (ar - al) * (a(ii,jj,j) - (al + ar)/2.) >  (ar - al)**2/6. ) al = 3.*a(ii,jj,j) - 2. * ar
    if ( (ar - al) * (a(ii,jj,j) - (al + ar)/2.) < -(ar - al)**2/6. ) ar = 3.*a(ii,jj,j) - 2. * al
    !Computed these coefficients from the edge values and cell mean in Maple. Assumes normalized coordinates: xi=(x-x0)/dx
    coefs(ii,jj,0,j) = 1.5 * a(ii,jj,j) - ( al + ar ) / 4.
    coefs(ii,jj,1,j) = ar - al
    coefs(ii,jj,2,j) = -6. * a(ii,jj,j) + 3. * ( al + ar )
  endif
  call syncthreads()

  !If we're not using a mirrored boundary condition, then make the two cells bordering the top and bottom
  !material boundaries piecewise constant. Zeroing out the first and second moments, and setting the zeroth
  !moment to the cell mean is sufficient to maintain conservation.
  if (vert_remap_q_alg == 2) then
    if ( j == 0 ) then
      coefs(ii,jj,0,1:2) = a(ii,jj,1:2)
      coefs(ii,jj,1:2,1:2) = 0.
      coefs(ii,jj,0,nlev-1:nlev) = a(ii,jj,nlev-1:nlev)
      coefs(ii,jj,1:2,nlev-1:nlev) = 0.D0
    endif
  endif
  call syncthreads()
end function compute_ppm_d



!THis compute grid-based coefficients from Collela & Woodward 1984.
function compute_ppm_grids_d( dx )   result(rslt)
  use control_mod, only: vert_remap_q_alg
  implicit none
  real(kind=real_kind), intent(in) :: dx(-1:nlev+2)  !grid spacings
  real(kind=real_kind)             :: rslt(10,0:nlev+1)  !grid spacings
  integer :: j
  integer :: indB, indE

  !Calculate grid-based coefficients for stage 1 of compute_ppm
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-1
  else
    indB = 0
    indE = nlev+1
  endif
  do j = indB , indE
    rslt( 1,j) = dx(j) / ( dx(j-1) + dx(j) + dx(j+1) )
    rslt( 2,j) = ( 2.*dx(j-1) + dx(j) ) / ( dx(j+1) + dx(j) )
    rslt( 3,j) = ( dx(j) + 2.*dx(j+1) ) / ( dx(j-1) + dx(j) )
  enddo

  !Caculate grid-based coefficients for stage 2 of compute_ppm
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-2
  else
    indB = 0
    indE = nlev
  endif
  do j = indB , indE
    rslt( 4,j) = dx(j) / ( dx(j) + dx(j+1) )
    rslt( 5,j) = 1. / sum( dx(j-1:j+2) )
    rslt( 6,j) = ( 2. * dx(j+1) * dx(j) ) / ( dx(j) + dx(j+1 ) )
    rslt( 7,j) = ( dx(j-1) + dx(j  ) ) / ( 2. * dx(j  ) + dx(j+1) )
    rslt( 8,j) = ( dx(j+2) + dx(j+1) ) / ( 2. * dx(j+1) + dx(j  ) )
    rslt( 9,j) = dx(j  ) * ( dx(j-1) + dx(j  ) ) / ( 2.*dx(j  ) +    dx(j+1) )
    rslt(10,j) = dx(j+1) * ( dx(j+1) + dx(j+2) ) / (    dx(j  ) + 2.*dx(j+1) )
  enddo
end function compute_ppm_grids_d


























!This uses the exact same model and reference grids and data as remap_Q, but it interpolates
!using PPM instead of splines.
subroutine remap_Q_ppm(Qdp,nx,qsize,dp1,dp2)
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass
  !
  use control_mod, only        : prescribed_wind, vert_remap_q_alg
  use perf_mod, only: t_startf, t_stopf
  implicit none
  integer,intent(in) :: nx,qsize
  real (kind=real_kind), intent(inout) :: Qdp(nx,nx,nlev,qsize)
  real (kind=real_kind), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
  ! Local Variables
  integer, parameter :: gs = 2                              !Number of cells to place in the ghost region
  real(kind=real_kind), dimension(       nlev+2 ) :: pio    !Pressure at interfaces for old grid
  real(kind=real_kind), dimension(       nlev+1 ) :: pin    !Pressure at interfaces for new grid
  real(kind=real_kind), dimension(       nlev+1 ) :: masso  !Accumulate mass up to each interface
  real(kind=real_kind), dimension(  1-gs:nlev+gs) :: ao     !Tracer value on old grid
  real(kind=real_kind), dimension(  1-gs:nlev+gs) :: dpo    !change in pressure over a cell for old grid
  real(kind=real_kind), dimension(  1-gs:nlev+gs) :: dpn    !change in pressure over a cell for old grid
  real(kind=real_kind), dimension(3,     nlev   ) :: coefs  !PPM coefficients within each cell
  real(kind=real_kind), dimension(       nlev   ) :: z1, z2
  real(kind=real_kind) :: ppmdx(10,0:nlev+1)  !grid spacings
  real(kind=real_kind) :: mymass, massn1, massn2
  integer :: i, j, k, q, kk, kid(nlev)

  call t_startf('remap_Q_ppm')
  do j = 1 , nx
    do i = 1 , nx
      
      pin(1)=0
      pio(1)=0
      do k=1,nlev
         dpn(k)=dp2(i,j,k)
         dpo(k)=dp1(i,j,k)
         pin(k+1)=pin(k)+dpn(k)
         pio(k+1)=pio(k)+dpo(k)
      enddo



      pio(nlev+2) = pio(nlev+1) + 1.  !This is here to allow an entire block of k threads to run in the remapping phase.
                                      !It makes sure there's an old interface value below the domain that is larger.
      pin(nlev+1) = pio(nlev+1)       !The total mass in a column does not change.
                                      !Therefore, the pressure of that mass cannot either.
      !Fill in the ghost regions with mirrored values. if vert_remap_q_alg is defined, this is of no consequence.
      do k = 1 , gs
        dpo(1   -k) = dpo(       k)
        dpo(nlev+k) = dpo(nlev+1-k)
      enddo

      !Compute remapping intervals once for all tracers. Find the old grid cell index in which the
      !k-th new cell interface resides. Then integrate from the bottom of that old cell to the new
      !interface location. In practice, the grid never deforms past one cell, so the search can be
      !simplified by this. Also, the interval of integration is usually of magnitude close to zero
      !or close to dpo because of minimial deformation.
      !Numerous tests confirmed that the bottom and top of the grids match to machine precision, so
      !I set them equal to each other.
      do k = 1 , nlev
        kk = k  !Keep from an order n^2 search operation by assuming the old cell index is close.
        !Find the index of the old grid cell in which this new cell's bottom interface resides.
        do while ( pio(kk) <= pin(k+1) )
          kk = kk + 1
        enddo
        kk = kk - 1                   !kk is now the cell index we're integrating over.
        if (kk == nlev+1) kk = nlev   !This is to keep the indices in bounds.
                                      !Top bounds match anyway, so doesn't matter what coefficients are used
        kid(k) = kk                   !Save for reuse
        z1(k) = -0.5D0                !This remapping assumes we're starting from the left interface of an old grid cell
                                      !In fact, we're usually integrating very little or almost all of the cell in question
        z2(k) = ( pin(k+1) - ( pio(kk) + pio(kk+1) ) * 0.5 ) / dpo(kk)  !PPM interpolants are normalized to an independent
                                                                        !coordinate domain [-0.5,0.5].
      enddo

      !This turned out a big optimization, remembering that only parts of the PPM algorithm depends on the data, namely the
      !limiting. So anything that depends only on the grid is pre-computed outside the tracer loop.
      ppmdx(:,:) = compute_ppm_grids( dpo )

      !From here, we loop over tracers for only those portions which depend on tracer data, which includes PPM limiting and
      !mass accumulation
      do q = 1 , qsize
        !Accumulate the old mass up to old grid cell interface locations to simplify integration
        !during remapping. Also, divide out the grid spacing so we're working with actual tracer
        !values and can conserve mass. The option for ifndef ZEROHORZ I believe is there to ensure
        !tracer consistency for an initially uniform field. I copied it from the old remap routine.
        masso(1) = 0.
        do k = 1 , nlev
          ao(k) = Qdp(i,j,k,q)
          masso(k+1) = masso(k) + ao(k) !Accumulate the old mass. This will simplify the remapping
          ao(k) = ao(k) / dpo(k)        !Divide out the old grid spacing because we want the tracer mixing ratio, not mass.
        enddo
        !Fill in ghost values. Ignored if vert_remap_q_alg == 2
        do k = 1 , gs
          ao(1   -k) = ao(       k)
          ao(nlev+k) = ao(nlev+1-k)
        enddo
        !Compute monotonic and conservative PPM reconstruction over every cell
        coefs(:,:) = compute_ppm( ao , ppmdx )
        !Compute tracer values on the new grid by integrating from the old cell bottom to the new
        !cell interface to form a new grid mass accumulation. Taking the difference between
        !accumulation at successive interfaces gives the mass inside each cell. Since Qdp is
        !supposed to hold the full mass this needs no normalization.
        massn1 = 0.
        do k = 1 , nlev
          kk = kid(k)
          massn2 = masso(kk) + integrate_parabola( coefs(:,kk) , z1(k) , z2(k) ) * dpo(kk)
          Qdp(i,j,k,q) = massn2 - massn1
          massn1 = massn2
        enddo
      enddo
    enddo
  enddo
  call t_stopf('remap_Q_ppm')
end subroutine remap_Q_ppm


!=======================================================================================================! 


!THis compute grid-based coefficients from Collela & Woodward 1984.
function compute_ppm_grids( dx )   result(rslt)
  use control_mod, only: vert_remap_q_alg
  implicit none
  real(kind=real_kind), intent(in) :: dx(-1:nlev+2)  !grid spacings
  real(kind=real_kind)             :: rslt(10,0:nlev+1)  !grid spacings
  integer :: j
  integer :: indB, indE

  !Calculate grid-based coefficients for stage 1 of compute_ppm
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-1
  else
    indB = 0
    indE = nlev+1
  endif
  do j = indB , indE
    rslt( 1,j) = dx(j) / ( dx(j-1) + dx(j) + dx(j+1) )
    rslt( 2,j) = ( 2.*dx(j-1) + dx(j) ) / ( dx(j+1) + dx(j) )
    rslt( 3,j) = ( dx(j) + 2.*dx(j+1) ) / ( dx(j-1) + dx(j) )
  enddo

  !Caculate grid-based coefficients for stage 2 of compute_ppm
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-2
  else
    indB = 0
    indE = nlev
  endif
  do j = indB , indE
    rslt( 4,j) = dx(j) / ( dx(j) + dx(j+1) )
    rslt( 5,j) = 1. / sum( dx(j-1:j+2) )
    rslt( 6,j) = ( 2. * dx(j+1) * dx(j) ) / ( dx(j) + dx(j+1 ) )
    rslt( 7,j) = ( dx(j-1) + dx(j  ) ) / ( 2. * dx(j  ) + dx(j+1) )
    rslt( 8,j) = ( dx(j+2) + dx(j+1) ) / ( 2. * dx(j+1) + dx(j  ) )
    rslt( 9,j) = dx(j  ) * ( dx(j-1) + dx(j  ) ) / ( 2.*dx(j  ) +    dx(j+1) )
    rslt(10,j) = dx(j+1) * ( dx(j+1) + dx(j+2) ) / (    dx(j  ) + 2.*dx(j+1) )
  enddo
end function compute_ppm_grids

!=======================================================================================================! 



!This computes a limited parabolic interpolant using a net 5-cell stencil, but the stages of computation are broken up into 3 stages
function compute_ppm( a , dx )    result(coefs)
  use control_mod, only: vert_remap_q_alg
  implicit none
  real(kind=real_kind), intent(in) :: a    (    -1:nlev+2)  !Cell-mean values
  real(kind=real_kind), intent(in) :: dx   (10,  0:nlev+1)  !grid spacings
  real(kind=real_kind) ::             coefs(0:2,   nlev  )  !PPM coefficients (for parabola)
  real(kind=real_kind) :: ai (0:nlev  )                     !fourth-order accurate, then limited interface values
  real(kind=real_kind) :: dma(0:nlev+1)                     !An expression from Collela's '84 publication
  real(kind=real_kind) :: da                                !Ditto
  ! Hold expressions based on the grid (which are cumbersome).
  real(kind=real_kind) :: dx1, dx2, dx3, dx4, dx5, dx6, dx7, dx8, dx9, dx10
  real(kind=real_kind) :: al, ar                            !Left and right interface values for cell-local limiting
  integer :: j
  integer :: indB, indE

  ! Stage 1: Compute dma for each cell, allowing a 1-cell ghost stencil below and above the domain
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-1
  else
    indB = 0
    indE = nlev+1
  endif
  do j = indB , indE
    da = dx(1,j) * ( dx(2,j) * ( a(j+1) - a(j) ) + dx(3,j) * ( a(j) - a(j-1) ) )
    dma(j) = minval( (/ abs(da) , 2. * abs( a(j) - a(j-1) ) , 2. * abs( a(j+1) - a(j) ) /) ) * sign(1.D0,da)
    if ( ( a(j+1) - a(j) ) * ( a(j) - a(j-1) ) <= 0. ) dma(j) = 0.
  enddo

  ! Stage 2: Compute ai for each cell interface in the physical domain (dimension nlev+1)
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-2
  else
    indB = 0
    indE = nlev
  endif
  do j = indB , indE
    ai(j) = a(j) + dx(4,j) * ( a(j+1) - a(j) ) + dx(5,j) * ( dx(6,j) * ( dx(7,j) - dx(8,j) ) &
         * ( a(j+1) - a(j) ) - dx(9,j) * dma(j+1) + dx(10,j) * dma(j) )
  enddo

  ! Stage 3: Compute limited PPM interpolant over each cell in the physical domain
  ! (dimension nlev) using ai on either side and ao within the cell.
  if (vert_remap_q_alg == 2) then
    indB = 3
    indE = nlev-2
  else
    indB = 1
    indE = nlev
  endif
  do j = indB , indE
    al = ai(j-1)
    ar = ai(j  )
    if ( (ar - a(j)) * (a(j) - al) <= 0. ) then
      al = a(j)
      ar = a(j)
    endif
    if ( (ar - al) * (a(j) - (al + ar)/2.) >  (ar - al)**2/6. ) al = 3.*a(j) - 2. * ar
    if ( (ar - al) * (a(j) - (al + ar)/2.) < -(ar - al)**2/6. ) ar = 3.*a(j) - 2. * al
    !Computed these coefficients from the edge values and cell mean in Maple. Assumes normalized coordinates: xi=(x-x0)/dx
    coefs(0,j) = 1.5 * a(j) - ( al + ar ) / 4.
    coefs(1,j) = ar - al
    coefs(2,j) = -6. * a(j) + 3. * ( al + ar )
  enddo

  !If we're not using a mirrored boundary condition, then make the two cells bordering the top and bottom
  !material boundaries piecewise constant. Zeroing out the first and second moments, and setting the zeroth
  !moment to the cell mean is sufficient to maintain conservation.
  if (vert_remap_q_alg == 2) then
    coefs(0,1:2) = a(1:2)
    coefs(1:2,1:2) = 0.
    coefs(0,nlev-1:nlev) = a(nlev-1:nlev)
    coefs(1:2,nlev-1:nlev) = 0.D0
  endif
end function compute_ppm

!=======================================================================================================! 


!Simple function computes the definite integral of a parabola in normalized coordinates, xi=(x-x0)/dx,
!given two bounds. Make sure this gets inlined during compilation.
function integrate_parabola( a , x1 , x2 )    result(mass)
  implicit none
  real(kind=real_kind), intent(in) :: a(0:2)  !Coefficients of the parabola
  real(kind=real_kind), intent(in) :: x1      !lower domain bound for integration
  real(kind=real_kind), intent(in) :: x2      !upper domain bound for integration
  real(kind=real_kind)             :: mass
  mass = a(0) * (x2 - x1) + a(1) * (x2 ** 2 - x1 ** 2) / 0.2D1 + a(2) * (x2 ** 3 - x1 ** 3) / 0.3D1
end function integrate_parabola






#endif
end module cuda_mod


