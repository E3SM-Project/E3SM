!This is where all of the PGI CUDA FORTRAN code will go, and these routines will be called from prim_advection_mod.
!This is compiled regardless, but PGI-specific calls are always wrapped in the USE_CUDA_FORTRAN ifdefs that are automagically
!activated when -Mcuda is specified during compilation with a PGI compiler. Thus, it will be ignored unless explicitly
!activated by the user
!
!As a general rule, all of the routines in here will be called within a threaded context (assuming ELEMENT_OPENMP is not
!deifned), and therefore, we enforce BARRIERS, MASTERS, and SINGLES from within these routines rather than outside them.
!This is to minimize the visible code impacts on the existing CPU code.

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!#define CUDA_DEBUG

#ifdef CUDA_DEBUG
#define _CHECK(line) ierr = cudaThreadSynchronize(); if (ierr .ne. 0) stop line
#else
#define _CHECK(line)
#endif


module cuda_mod
#if USE_CUDA_FORTRAN

! NP > 4 is not supported due to shared memory constraints
!BSINGH -  NAG compiler doesn't support "#error" CPP directive.
!Therefore, getting rid of "#" preceding "error" so that the code still fails at compile time
#if NP > 4
error CUDA Fortran build only supported with NP <= 4
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
  public :: cudaThreadSynchronize_wrap
  public :: precompute_divdp_cuda

  integer,parameter :: numk_eul = 6
  integer,parameter :: numk_hyp = 6
  integer,parameter :: numk_lim2d = 15

  !This is from prim_advection_mod.F90
  type(EdgeBuffer_t) :: edgeAdv, edgeAdv3, edgeAdvQ3, edgeAdvQ2, edgeAdvDSS
  integer,parameter :: DSSeta = 1
  integer,parameter :: DSSomega = 2
  integer,parameter :: DSSdiv_vdp_ave = 3
  integer,parameter :: DSSno_var = -1

  !Device arrays
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:,:,:) :: qdp_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:)     :: qdp1_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:,:)   :: qtens_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:,:)   :: qtens_biharmonic_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: spheremp_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:)     :: dpdiss_ave_d
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
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:)     :: divdp_proj_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:)     :: dp_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:,:)     :: dp_star_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: qmin_d
  real (kind=real_kind),device,allocatable,dimension(:,:,:)       :: qmax_d
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
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:)     :: qdp1_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:,:)   :: vstar_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:,:)   :: qtens_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:)     :: dp_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:)     :: divdp_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:)     :: divdp_proj_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:)     :: dp_star_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:)     :: dp_np1_h
  real(kind=real_kind),pinned,allocatable,dimension(:,:,:,:)     :: dpdiss_ave_h
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
  real(kind=real_kind), device :: rrearth_d


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
    use physical_constants, only: rrearth
    implicit none
    type(element_t)   , intent(in) :: elem(:)
    type(derivative_t), intent(in) :: deriv
    type(hvcoord_t)   , intent(in) :: hvcoord

    type (Cycle_t),pointer    :: pCycle
    type (Schedule_t),pointer :: pSchedule
    integer                   :: ie , ierr , icycle , iPtr , rank , nSendCycles , nRecvCycles , nlyr , mx_send_len , mx_recv_len , n
!    real(kind=real_kind)      :: dinv_t(np , np , 2 , 2)
    type (dim3)               :: griddim , blockdim
    logical,allocatable,dimension(:,:) :: send_elem_mask
    logical,allocatable,dimension(:,:) :: recv_elem_mask
    logical,allocatable,dimension(:)   :: elem_computed
    integer :: total_work

    rrearth_d = rrearth

    max_neigh_edges_d = max_neigh_edges
    max_corner_elem_d = max_corner_elem

    cuda_streams = nelemd*8
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
    allocate( qmin_d                   (nlev,qsize_d                 ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qmax_d                   (nlev,qsize_d                 ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qdp_d                    (np,np,nlev,qsize_d,timelevels,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qdp1_d                   (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dpdiss_ave_d             (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qtens_d                  (np,np,nlev,qsize_d           ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qtens_biharmonic_d       (np,np,nlev,qsize_d           ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( spheremp_d               (np,np                        ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( rspheremp_d              (np,np                        ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dinv_d                   (np,np,2,2                    ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( deriv_dvv_d              (np,np                               ) , stat = ierr ); _CHECK(__LINE__)
    allocate( variable_hyperviscosity_d(np,np                        ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( metdet_d                 (np,np                        ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( rmetdet_d                (np,np                        ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dpdiss_biharmonic_d      (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( vstar_d                  (np,np,nlev,2                 ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( divdp_d                  (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( divdp_proj_d             (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dp_d                     (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( hyai_d                   (      nlev+1                        ) , stat = ierr ); _CHECK(__LINE__)
    allocate( hybi_d                   (      nlev+1                        ) , stat = ierr ); _CHECK(__LINE__)
    allocate( dp_star_d                (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( reverse_d                (max_neigh_edges              ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( putmapP_d                (max_neigh_edges              ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( getmapP_d                (max_neigh_edges              ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( edgebuf_d                (nlev*qsize_d,nbuf                   ) , stat = ierr ); _CHECK(__LINE__)
    allocate( recv_internal_indices_d  (nelemd                              ) , stat = ierr ); _CHECK(__LINE__)
    allocate( recv_external_indices_d  (nelemd                              ) , stat = ierr ); _CHECK(__LINE__)
    allocate( dpo_d                    (np,np,   1-gs:nlev+gs        ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( ppmdx_d                  (np,np,10,   0:nlev+1         ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( z1_d                     (np,np,        nlev           ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( z2_d                     (np,np,        nlev           ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( kid_d                    (np,np,        nlev           ,nelemd) , stat = ierr ); _CHECK(__LINE__)

    allocate( qdp_h                    (np,np,nlev,qsize_d,timelevels,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qdp1_h                   (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dpdiss_ave_h             (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qmin                     (nlev,qsize_d                 ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qmax                     (nlev,qsize_d                 ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( vstar_h                  (np,np,nlev,2                 ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qtens_h                  (np,np,nlev,qsize_d           ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( qtens_biharmonic         (np,np,nlev,qsize_d           ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dp_h                     (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( divdp_h                  (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( divdp_proj_h             (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dp_star_h                (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dpdiss_biharmonic_h      (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dp_np1_h                 (np,np,nlev                   ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( dpo_h                    (np,np,   1-gs:nlev+gs        ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( ppmdx_h                  (np,np,10,   0:nlev+1         ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( z1_h                     (np,np,        nlev           ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( z2_h                     (np,np,        nlev           ,nelemd) , stat = ierr ); _CHECK(__LINE__)
    allocate( kid_h                    (np,np,        nlev           ,nelemd) , stat = ierr ); _CHECK(__LINE__)

    ! The PGI compiler with cuda enabled errors when allocating arrays of zero
    !   size - here when using only one MPI task
    if (nSendCycles > 0) then
      allocate( sendbuf_h                (nlev*qsize_d,mx_send_len,nSendCycles) , stat = ierr ); _CHECK(__LINE__)
      allocate( send_elem_mask           (nelemd,nSendCycles                  ) , stat = ierr ); _CHECK(__LINE__)
      allocate( send_nelem               (       nSendCycles                  ) , stat = ierr ); _CHECK(__LINE__)
      allocate( send_indices             (nelemd,nSendCycles                  ) , stat = ierr ); _CHECK(__LINE__)
      allocate( d2h_done                 (nSendCycles                         ) , stat = ierr ); _CHECK(__LINE__)
      allocate( msg_sent                 (nSendCycles                         ) , stat = ierr ); _CHECK(__LINE__)
    endif

    if (nRecvCycles > 0) then
      allocate( recvbuf_h                (nlev*qsize_d,mx_recv_len,nRecvCycles) , stat = ierr ); _CHECK(__LINE__)
      allocate( recvbuf_d                (nlev*qsize_d,mx_recv_len,nRecvCycles) , stat = ierr ); _CHECK(__LINE__)
      allocate( recv_elem_mask           (nelemd,nRecvCycles                  ) , stat = ierr ); _CHECK(__LINE__)
      allocate( recv_nelem               (       nRecvCycles                  ) , stat = ierr ); _CHECK(__LINE__)
      allocate( recv_indices             (nelemd,nRecvCycles                  ) , stat = ierr ); _CHECK(__LINE__)
      allocate( h2d_done                 (nRecvCycles                         ) , stat = ierr ); _CHECK(__LINE__)
      allocate( msg_rcvd                 (nRecvCycles                         ) , stat = ierr ); _CHECK(__LINE__)
    endif 

    allocate( recv_internal_indices    (nelemd                                ) , stat = ierr ); _CHECK(__LINE__)
    allocate( recv_external_indices    (nelemd                                ) , stat = ierr ); _CHECK(__LINE__)
    allocate( elem_computed            (nelemd                                ) , stat = ierr ); _CHECK(__LINE__)
    write(*,*) "send data from host to device"
    !Copy over data to the device
    ierr = cudaMemcpy( deriv_dvv_d , deriv%dvv    , size( deriv%dvv    ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
    ierr = cudaMemcpy( hyai_d      , hvcoord%hyai , size( hvcoord%hyai ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
    ierr = cudaMemcpy( hybi_d      , hvcoord%hybi , size( hvcoord%hybi ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
    do ie = 1,nelemd
!      dinv_t(:,:,1,1) = elem(ie)%dinv(1,1,:,:)
!      dinv_t(:,:,1,2) = elem(ie)%dinv(1,2,:,:)
!      dinv_t(:,:,2,1) = elem(ie)%dinv(2,1,:,:)
!      dinv_t(:,:,2,2) = elem(ie)%dinv(2,2,:,:)
      ierr = cudaMemcpy( dinv_d                   (1,1,1,1  ,ie) , elem(ie)%state%dinv              , size(elem(ie)%state%Qdp              ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)

      ierr = cudaMemcpy( qdp_d                    (1,1,1,1,1,ie) , elem(ie)%state%Qdp               , size(elem(ie)%state%Qdp              ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
      ierr = cudaMemcpy( spheremp_d               (1,1      ,ie) , elem(ie)%spheremp                , size(elem(ie)%spheremp               ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
      ierr = cudaMemcpy( rspheremp_d              (1,1      ,ie) , elem(ie)%rspheremp               , size(elem(ie)%rspheremp              ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
!      ierr = cudaMemcpy( dinv_d                   (1,1,1,1  ,ie) , dinv_t                           , size(dinv_t                          ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
      ierr = cudaMemcpy( metdet_d                 (1,1      ,ie) , elem(ie)%metdet                  , size(elem(ie)%metdet                 ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
      ierr = cudaMemcpy( rmetdet_d                (1,1      ,ie) , elem(ie)%rmetdet                 , size(elem(ie)%rmetdet                ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
      ierr = cudaMemcpy( putmapP_d                (1        ,ie) , elem(ie)%desc%putmapP            , size(elem(ie)%desc%putmapP           ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
      ierr = cudaMemcpy( getmapP_d                (1        ,ie) , elem(ie)%desc%getmapP            , size(elem(ie)%desc%getmapP           ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
      ierr = cudaMemcpy( reverse_d                (1        ,ie) , elem(ie)%desc%reverse            , size(elem(ie)%desc%reverse           ) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
      ierr = cudaMemcpy( variable_hyperviscosity_d(1,1      ,ie) , elem(ie)%variable_hyperviscosity , size(elem(ie)%variable_hyperviscosity) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
    enddo

    write(*,*) "edgebuffers"
    !These have to be in a threaded region or they complain and die
!$OMP END MASTER
    call initEdgeBuffer(edgeAdv   ,qsize_d*nlev  )
    call initEdgeBuffer(edgeAdvDSS,      3*nlev  )
    call initEdgeBuffer(edgeAdvQ2 ,qsize_d*nlev*2)
    call initEdgeBuffer(edgeAdvQ3 ,qsize_d*nlev*3)
!$OMP MASTER

    write(*,*) "initial kernel"
    !This needs to run because we need accurate timing during out cuda profiler runs. Initial kernel runs incur overheads, so we do this here
    blockdim = dim3(1,1,1)
    griddim  = dim3(1,1,1)
    call warmup <<< griddim , blockdim >>> ( ie ); _CHECK(__LINE__)
    ierr = cudaThreadSynchronize()

    do n = 0 , cuda_streams
      ierr = cudaStreamCreate(streams(n)); _CHECK(__LINE__)
      ierr = cudaStreamCreate(streams2(n)); _CHECK(__LINE__)
    enddo
    ierr = cudaDeviceSetCacheConfig(cudaFuncCachePreferShared); _CHECK(__LINE__)
    ierr = cudaEventCreate(timer1); _CHECK(__LINE__)
    ierr = cudaEventCreate(timer2); _CHECK(__LINE__)

    write(*,*) "Dividing elements among cycles in which they participate"
    !For efficient MPI, PCI-e, packing, and unpacking, we need to separate out the cycles by dependence. Once on cycle has packed, then stage the PCI-e D2H, MPI, PCI-e H2D, & internal unpack
    !We begin by testing what elements contribute to packing in what cycle's MPI data.
    do ie = 1,nelemd
      if (nRecvCycles > 0) recv_elem_mask(ie,:) = .false.
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

    write(*,*) "Sending element & cycle informationt to device"
    ierr = cudaMemcpy(recv_internal_indices_d,recv_internal_indices,size(recv_internal_indices),cudaMemcpyHostToDevice); _CHECK(__LINE__)
    ierr = cudaMemcpy(recv_external_indices_d,recv_external_indices,size(recv_external_indices),cudaMemcpyHostToDevice); _CHECK(__LINE__)

    write(*,*)"done cuda_mod_init"
!$OMP END MASTER
!$OMP BARRIER
  end subroutine cuda_mod_init



  !Meaningless kernel just to get initial kernel overheads out of the way.
  attributes(global) subroutine warmup(a)
    integer,value :: a
    a = 2.0 * a
  end subroutine warmup



  function cudaThreadSynchronize_wrap() result(ierr)
    integer :: ierr
    ierr = cudaThreadSynchronize()
  end function cudaThreadSynchronize_wrap



  subroutine copy_qdp_h2d( elem , nt )
    use element_mod, only: element_t
    use perf_mod, only: t_startf, t_stopf
    implicit none
    type(element_t), intent(in   ) :: elem(:)
    integer        , intent(in   ) :: nt
    integer :: ierr , ie
    !$OMP BARRIER
    !$OMP MASTER
    ierr = cudaThreadSynchronize()
    call t_startf('CUDA QDP H2D')
    do ie = 1,nelemd
      ierr = cudaMemcpy( qdp_d(1,1,1,1,nt,ie) , elem(ie)%state%qdp(1,1,1,1,nt) , size(elem(ie)%state%qdp(:,:,:,:,nt)) , cudaMemcpyHostToDevice ); _CHECK(__LINE__)
    enddo
    ierr = cudaThreadSynchronize()
    call t_stopf('CUDA QDP H2D')
    !$OMP END MASTER
    !$OMP BARRIER
  end subroutine copy_qdp_h2d



  subroutine copy_qdp_d2h( elem , nt )
    use element_mod, only: element_t
    use perf_mod, only: t_startf, t_stopf
    implicit none
    type(element_t), intent(in   ) :: elem(:)
    integer        , intent(in   ) :: nt
    integer :: ierr , ie
    !$OMP BARRIER
    !$OMP MASTER
    ierr = cudaThreadSynchronize()
    call t_startf('CUDA QDP D2H')
    do ie = 1,nelemd
      ierr = cudaMemcpy( elem(ie)%state%qdp(1,1,1,1,nt) , qdp_d(1,1,1,1,nt,ie) , size(elem(ie)%state%qdp(:,:,:,:,nt)) , cudaMemcpyDeviceToHost ); _CHECK(__LINE__)
    enddo
    ierr = cudaThreadSynchronize()
    call t_stopf('CUDA QDP D2H')
    !$OMP END MASTER
    !$OMP BARRIER
  end subroutine copy_qdp_d2h



  subroutine precompute_divdp_cuda( elem , hybrid , deriv , dt , nets , nete , n0_qdp )
    use element_mod, only: element_t
    use derivative_mod, only: derivative_t, divergence_sphere
    use hybrid_mod, only: hybrid_t
    use edge_mod, only: edgeVpack,edgeVunpack
    use bndry_mod, only: bndry_exchangeV
    use perf_mod, only: t_startf, t_stopf
    implicit none
    type(element_t)      , intent(inout) :: elem(:)
    type (hybrid_t)      , intent(in   ) :: hybrid
    type (derivative_t)  , intent(in   ) :: deriv
    real(kind=real_kind) , intent(in   ) :: dt
    integer              , intent(in   ) :: nets , nete , n0_qdp
    integer :: ie , k
    integer :: ierr
    type(dim3) :: blockdim , griddim

    call t_startf('precompute_divdp_cuda')

    do ie = nets , nete
      qdp1_h(:,:,:,ie) = elem(ie)%state%Qdp(:,:,:,1,n0_qdp)
      dpdiss_ave_h(:,:,:,ie) = elem(ie)%derived%dpdiss_ave(:,:,:)
    enddo
    !$OMP BARRIER
    !$OMP MASTER
    ierr = cudaMemcpyAsync( dpdiss_ave_d , dpdiss_ave_h , size( dpdiss_ave_h ) , cudaMemcpyHostToDevice , streams(1) ); _CHECK(__LINE__)
    ierr = cudaMemcpyAsync( qdp1_d , qdp1_h , size( qdp1_h ) , cudaMemcpyHostToDevice , streams(1) ); _CHECK(__LINE__)
    blockdim = dim3( np     , np , nlev )
    griddim  = dim3( nelemd , 1  , 1    )
    call unpack_qdp1<<<griddim,blockdim,0,streams(1)>>>( qdp1_d , qdp_d , n0_qdp , nets , nete ); _CHECK(__LINE__)
    !$OMP END MASTER

    do ie = nets , nete
      do k = 1 , nlev   ! div( U dp Q), 
        elem(ie)%derived%divdp(:,:,k) = divergence_sphere(elem(ie)%derived%vn0(:,:,:,k),deriv,elem(ie))
        elem(ie)%derived%eta_dot_dpdn(:,:,k) = elem(ie)%spheremp(:,:) * elem(ie)%derived%eta_dot_dpdn(:,:,k)
        elem(ie)%derived%omega_p     (:,:,k) = elem(ie)%spheremp(:,:) * elem(ie)%derived%omega_p     (:,:,k)
        elem(ie)%derived%divdp_proj  (:,:,k) = elem(ie)%spheremp(:,:) * elem(ie)%derived%divdp       (:,:,k)
      enddo
      call edgeVpack(edgeAdvDSS,elem(ie)%derived%eta_dot_dpdn(:,:,1:nlev),nlev,0*nlev,elem(ie)%desc)
      call edgeVpack(edgeAdvDSS,elem(ie)%derived%omega_p     (:,:,1:nlev),nlev,1*nlev,elem(ie)%desc)
      call edgeVpack(edgeAdvDSS,elem(ie)%derived%divdp_proj  (:,:,1:nlev),nlev,2*nlev,elem(ie)%desc)
    enddo
    call t_startf('eus_cuda_bexchV')
    call bndry_exchangeV(hybrid,edgeAdvDSS)
    call t_stopf('eus_cuda_bexchV')
    do ie = nets , nete
      call edgeVunpack(edgeAdvDSS,elem(ie)%derived%eta_dot_dpdn(:,:,1:nlev),nlev,0*nlev,elem(ie)%desc)
      call edgeVunpack(edgeAdvDSS,elem(ie)%derived%omega_p     (:,:,1:nlev),nlev,1*nlev,elem(ie)%desc)
      call edgeVunpack(edgeAdvDSS,elem(ie)%derived%divdp_proj  (:,:,1:nlev),nlev,2*nlev,elem(ie)%desc)
      do k = 1 , nlev
        elem(ie)%derived%eta_dot_dpdn(:,:,k) = elem(ie)%derived%eta_dot_dpdn(:,:,k) * elem(ie)%rspheremp(:,:)
        elem(ie)%derived%omega_p     (:,:,k) = elem(ie)%derived%omega_p     (:,:,k) * elem(ie)%rspheremp(:,:)
        elem(ie)%derived%divdp_proj  (:,:,k) = elem(ie)%derived%divdp_proj  (:,:,k) * elem(ie)%rspheremp(:,:)
      enddo
    enddo

    call t_stopf('precompute_divdp_cuda')


  end subroutine precompute_divdp_cuda



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
  real(kind=real_kind) :: dp0,qmintmp,qmaxtmp
  integer :: ie,q,i,j,k
  integer :: rhs_viss = 0

  integer :: ierr
  type(dim3) :: blockdim , griddim

  !call t_startf('euler_step')
  call t_startf('euler_step_cuda')

  if (limiter_option == 8) then
    write(*,*) 'CUDA_MOD IS NOT INTENDED FOR USE WITH LIMITER_OPTION == 8 AT THIS TIME!'
    write(*,*) 'PLEASE USE LIMITER_OPTION == 0 WHEN THE GPU OPTION IS ENABLED!'
    stop
  endif

  rhs_viss = 0
  !   2D Advection step
  do ie = nets , nete
    ! Compute velocity used to advance Qdp 
    ! derived variable divdp_proj() (DSS'd version of divdp) will only be correct on 2nd and 3rd stage
    ! but that's ok because rhs_multiplier=0 on the first stage:
    dp_h(:,:,:,ie) = elem(ie)%derived%dp(:,:,:) - rhs_multiplier * dt * elem(ie)%derived%divdp_proj(:,:,:) 
    vstar_h(:,:,:,1,ie) = elem(ie)%derived%vn0(:,:,1,:) / dp_h(:,:,:,ie)
    vstar_h(:,:,:,2,ie) = elem(ie)%derived%vn0(:,:,2,:) / dp_h(:,:,:,ie)
  enddo
!$OMP BARRIER
!$OMP MASTER
  ierr = cudaMemcpyAsync( dp_d    , dp_h    , size( dp_h    ) , cudaMemcpyHostToDevice , streams(1) ); _CHECK(__LINE__)
  ierr = cudaMemcpyAsync( vstar_d , vstar_h , size( vstar_h ) , cudaMemcpyHostToDevice , streams(1) ); _CHECK(__LINE__)
  blockdim = dim3( np*np*numk_eul , 1 , 1 )
  griddim  = dim3( int(ceiling(dble(nlev)/numk_eul))*qsize_d*nelemd , 1 , 1 )
  call euler_step_kernel1<<<griddim,blockdim,0,streams(1)>>>( qdp_d , spheremp_d , vstar_d , qtens_biharmonic_d , metdet_d , rmetdet_d , dinv_d , deriv_dvv_d , &
                                                              n0_qdp , np1_qdp , rhs_viss , dt , 1 , nelemd ); _CHECK(__LINE__)
  if ( limiter_option == 4 ) then
  blockdim = dim3( np*np*numk_lim2d , 1 , 1 )
  griddim  = dim3( int(ceiling(dble(nlev)/numk_lim2d))*qsize_d*nelemd , 1 , 1 )
    call limiter2d_zero_kernel<<<griddim,blockdim,0,streams(1)>>>( qdp_d , 1 , nelemd , np1_qdp ); _CHECK(__LINE__)
  endif
  call t_startf('eus_cuda_peus')
  call pack_exchange_unpack_stage(np1_qdp,hybrid,qdp_d,timelevels)
  call t_stopf('eus_cuda_peus')
  blockdim = dim3( np*np   * nlev  , 1 , 1 )
  griddim  = dim3( qsize_d , nelemd , 1 )
  call euler_hypervis_kernel_last<<<griddim,blockdim,0,streams(1)>>>( qdp_d , rspheremp_d , 1 , nelemd , np1_qdp ); _CHECK(__LINE__)
!$OMP END MASTER
!$OMP BARRIER

  call t_stopf('euler_step_cuda')
  !call t_stopf('euler_step')
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
  blockdim = dim3( np      , np     , nlev )
  griddim  = dim3( qsize_d , nelemd , 1    )
  call qdp_time_avg_kernel<<<griddim,blockdim,0,streams(1)>>>( qdp_d , rkstage , n0_qdp , np1_qdp , 1 , nelemd ); _CHECK(__LINE__)
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
  use control_mod    , only : rsplit, nu_q, hypervis_order, hypervis_subcycle_q , nu_p
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
    ierr = cudaMemcpyAsync( dp_d         , dp_h         , size( dp_h )         , cudaMemcpyHostToDevice , streams(1) ); _CHECK(__LINE__)

    !KERNEL 1
    blockdim = dim3( np*np*numk_hyp , 1 , 1 )
    griddim  = dim3( int(ceiling(dble(nlev)/numk_hyp))*qsize_d*nelemd , 1 , 1 )
    call hypervis_kernel1<<<griddim,blockdim,0,streams(1)>>>( qdp_d , qtens_d , dp_d , dinv_d , variable_hyperviscosity_d , dpdiss_ave_d , spheremp_d , &
                                                              deriv_dvv_d , hyai_d , hybi_d , hvcoord%ps0 , 1 , nelemd , dt , nt_qdp , nu_p ); _CHECK(__LINE__)
    call t_startf('ahs_cuda_peus1')
    call pack_exchange_unpack_stage(1,hybrid,qtens_d,1)
    call t_stopf('ahs_cuda_peus1')

    !KERNEL 2
    blockdim = dim3( np*np*numk_hyp , 1 , 1 )
    griddim  = dim3( int(ceiling(dble(nlev)/numk_hyp))*qsize_d*nelemd , 1 , 1 )
    call hypervis_kernel2<<<griddim,blockdim,0,streams(1)>>>( qdp_d , qtens_d , dp_d , dinv_d , variable_hyperviscosity_d , spheremp_d , rspheremp_d , deriv_dvv_d , &
                                                              nu_q , 1 , nelemd , dt , nt_qdp ); _CHECK(__LINE__)
    blockdim = dim3( np*np*numk_lim2d , 1 , 1 )
    griddim  = dim3( int(ceiling(dble(nlev)/numk_lim2d))*qsize_d*nelemd , 1 , 1 )
    call limiter2d_zero_kernel<<<griddim,blockdim,0,streams(1)>>>(qdp_d,nets,nete,nt_qdp); _CHECK(__LINE__)

    call t_startf('ahs_cuda_peus2')
    call pack_exchange_unpack_stage(nt_qdp,hybrid,qdp_d,timelevels)
    call t_stopf('ahs_cuda_peus2')

    !KERNEL 3
    blockdim = dim3( np * np * nlev , 1 , 1 )
    griddim  = dim3( qsize_d , nelemd , 1 )
    call euler_hypervis_kernel_last<<<griddim,blockdim,0,streams(1)>>>( qdp_d , rspheremp_d , nets , nete , nt_qdp ); _CHECK(__LINE__)
    blockdim = dim3( np , np , nlev )
    griddim  = dim3( nelemd , 1 , 1 )
    call pack_qdp1<<<griddim,blockdim,0,streams(1)>>>( qdp1_d , qdp_d , nt_qdp , 1 , nelemd ); _CHECK(__LINE__)
    ierr = cudaMemcpyAsync( qdp1_h , qdp1_d , size( qdp1_h ) , cudaMemcpyDeviceToHost , streams(1) ); _CHECK(__LINE__)
    ierr = cudaStreamSynchronize(streams(1))
!$OMP END MASTER
!$OMP BARRIER
    do ie = nets , nete
      elem(ie)%state%Qdp(:,:,:,1,nt_qdp) = qdp1_h(:,:,:,ie)
    enddo
  enddo

  call t_stopf('advance_hypervis_scalar_cuda')

end subroutine advance_hypervis_scalar_cuda



attributes(global) subroutine pack_qdp1( qdp1 , qdp , nt , nets , nete )
  implicit none
  real(kind=real_kind), intent(  out) :: qdp1( np , np , nlev                        , nets:nete )
  real(kind=real_kind), intent(in   ) :: qdp ( np , np , nlev , qsize_d , timelevels , nets:nete )
  integer, value      , intent(in   ) :: nt , nets , nete
  integer :: i,j,k,ie
  i  = threadidx%x
  j  = threadidx%y
  k  = threadidx%z
  ie = blockidx%x
  qdp1(i,j,k,ie) = qdp(i,j,k,1,nt,ie)
end subroutine pack_qdp1



attributes(global) subroutine unpack_qdp1( qdp1 , qdp , nt , nets , nete )
  implicit none
  real(kind=real_kind), intent(in   ) :: qdp1( np , np , nlev                        , nets:nete )
  real(kind=real_kind), intent(  out) :: qdp ( np , np , nlev , qsize_d , timelevels , nets:nete )
  integer, value      , intent(in   ) :: nt , nets , nete
  integer :: i,j,k,ie
  i  = threadidx%x
  j  = threadidx%y
  k  = threadidx%z
  ie = blockidx%x
  qdp(i,j,k,1,nt,ie) = qdp1(i,j,k,ie)
end subroutine unpack_qdp1



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



attributes(global) subroutine euler_step_kernel1( Qdp , spheremp , vstar , qtens_biharmonic , metdet , rmetdet , dinv , deriv_dvv , n0_qdp , np1_qdp , rhs_viss , dt , nets , nete )
  implicit none
  real(kind=real_kind), dimension(np,np,nlev,qsize_d,timelevels,nets:nete), intent(inout) :: Qdp
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: spheremp
  real(kind=real_kind), dimension(np,np,nlev,2                 ,nets:nete), intent(in   ) :: vstar
  real(kind=real_kind), dimension(np,np,nlev,qsize_d           ,nets:nete), intent(in   ) :: qtens_biharmonic
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: metdet
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: rmetdet
  real(kind=real_kind), dimension(np,np,2,2                    ,nets:nete), intent(in   ) :: dinv
  real(kind=real_kind), dimension(np,np                                  ), intent(in   ) :: deriv_dvv
  real(kind=real_kind), value                                             , intent(in   ) :: dt
  integer, value                                                          , intent(in   ) :: n0_qdp, np1_qdp, rhs_viss, nets, nete
  integer :: ks
  integer :: i , j , k , kk , q , ie , ij
  real(kind=real_kind), shared :: deriv_dvv_s(np*np+1)
  real(kind=real_kind) :: qtmp, qtens

  ks = int(ceiling(dble(nlev)/numk_eul))

  !Define the indices
  i  = modulo( threadidx%x-1    ,np)+1
  j  = modulo((threadidx%x-1)/np,np)+1
  kk = (threadidx%x-1)/(np*np)+1
  k  = modulo(blockidx%x-1,ks)*numk_eul + kk
  q  = modulo((blockidx%x-1)/ks,qsize_d)+1
  ie =       ((blockidx%x-1)/ks)/qsize_d+1
  ij   =              (j-1)*np+i
  if (k  <= nlev .and. q  <= qsize_d .and. ie <= nete ) then
    if (kk == 1) deriv_dvv_s(ij) = deriv_dvv(i,j)
    qtmp = Qdp(i,j,k,q,n0_qdp,ie)
    qtens = qtmp - dt * divergence_sphere( i , j , k , q , ie , kk , ij , vstar , qtmp , metdet , rmetdet , dinv , deriv_dvv_s , nets , nete )
    if ( rhs_viss /= 0 ) qtens = qtens + qtens_biharmonic(i,j,k,q,ie)
    Qdp(i,j,k,q,np1_qdp,ie) = spheremp(i,j,ie) * qtens
  endif
end subroutine euler_step_kernel1



attributes(device) function divergence_sphere(i,j,k,q,ie,kk,ij,Vstar,qtmp,metdet,rmetdet,dinv,deriv_dvv,nets,nete) result(dp_star)
  implicit none
  integer,              intent(in) :: i, j, k, q, ie, kk, ij, nets, nete
  real(kind=real_kind), intent(in) :: Dinv     (np,np,2,2,nets:nete)
  real(kind=real_kind), intent(in) :: metdet   (np,np,nets:nete)
  real(kind=real_kind), intent(in) :: rmetdet  (np,np,nets:nete)
  real(kind=real_kind), intent(in) :: deriv_dvv(np*np+1)
  real(kind=real_kind), intent(in) :: Vstar    (np,np,nlev,2,nets:nete)
  real(kind=real_kind), intent(in), value :: qtmp
  real(kind=real_kind)             :: dp_star
  real(kind=real_kind), shared :: gv(np*np+1,numk_eul,2)
  integer :: s
  real(kind=real_kind) :: vvtemp, divtemp, vs1tmp, vs2tmp
  vs1tmp = vstar(i,j,k,1,ie) * metdet(i,j,ie) * qtmp
  vs2tmp = vstar(i,j,k,2,ie) * metdet(i,j,ie) * qtmp
  gv(ij,kk,1) = dinv(i,j,1,1,ie) * vs1tmp + dinv(i,j,1,2,ie) * vs2tmp
  gv(ij,kk,2) = dinv(i,j,2,1,ie) * vs1tmp + dinv(i,j,2,2,ie) * vs2tmp
  call syncthreads()
  divtemp = 0.0d0
  vvtemp   = 0.0d0
  do s = 1 , np
    divtemp = divtemp + deriv_dvv((i-1)*np+s) * gv((j-1)*np+s,kk,1)
    vvtemp  = vvtemp  + deriv_dvv((j-1)*np+s) * gv((s-1)*np+i,kk,2)
  end do
  dp_star = ( divtemp + vvtemp ) * ( rmetdet(i,j,ie) * rrearth_d )
end function divergence_sphere



attributes(global) subroutine limiter2d_zero_kernel(Qdp,nets,nete,np1)
  use kinds, only : real_kind
  use dimensions_mod, only : np, nlev
  implicit none
  real (kind=real_kind), intent(inout) :: Qdp(np,np,nlev,qsize_d,timelevels,nets:nete)
  integer, value       , intent(in   ) :: nets,nete,np1
  integer :: i, j, k, kk, q, ie, jj, ij, ijk, ks
  real (kind=real_kind) :: mass,mass_new
  real (kind=real_kind), shared :: Qdp_shared(np*np+PAD,nlev)
  real (kind=real_kind), shared :: mass_shared(nlev)
  real (kind=real_kind), shared :: mass_new_shared(nlev)

  ks = int(ceiling(dble(nlev)/numk_lim2d))

  i  = modulo( threadidx%x-1    ,np)+1
  j  = modulo((threadidx%x-1)/np,np)+1
  kk = (threadidx%x-1)/(np*np)+1
  k  = modulo(blockidx%x-1,ks)*numk_lim2d + kk
  q  = modulo((blockidx%x-1)/ks,qsize_d)+1
  ie =       ((blockidx%x-1)/ks)/qsize_d+1
  ij  = (j-1)*np+i
  ijk = (kk-1)*np*np+ij

  Qdp_shared(ij,kk) = Qdp(i,j,k,q,np1,ie)
  call syncthreads()

  if ( ijk <= numk_lim2d ) then
    mass = 0.
    do jj = 1 , np*np
      mass = mass + Qdp_shared(jj,ijk)
    enddo
    mass_shared(ijk) = mass
  endif
  call syncthreads()

  if ( mass_shared(kk)   < 0 ) Qdp_shared(ij,kk) = -Qdp_shared(ij,kk)
  if ( Qdp_shared(ij,kk) < 0 ) Qdp_shared(ij,kk) = 0
  call syncthreads()

  if ( ijk <= numk_lim2d ) then
    mass = 0.
    do jj = 1 , np*np
      mass = mass + Qdp_shared(jj,ijk)
    enddo
    mass_new_shared(ijk) = mass
  endif
  call syncthreads()

  ! now scale the all positive values to restore mass
  if ( mass_new_shared(kk) > 0 ) Qdp_shared(ij,kk) =  Qdp_shared(ij,kk) * abs(mass_shared(kk)) / mass_new_shared(kk)
  if ( mass_shared    (kk) < 0 ) Qdp_shared(ij,kk) = -Qdp_shared(ij,kk)
  Qdp(i,j,k,q,np1,ie) = Qdp_shared(ij,kk)
end subroutine limiter2d_zero_kernel



attributes(global) subroutine euler_hypervis_kernel_last( Qdp , rspheremp , nets , nete , np1 )
  implicit none
  real(kind=real_kind), dimension(np*np*nlev,qsize_d,timelevels,nets:nete), intent(inout) :: Qdp
  real(kind=real_kind), dimension(np*np                        ,nets:nete), intent(in   ) :: rspheremp
  integer, value                                                          , intent(in   ) :: nets , nete , np1
  integer :: ijk, ij, q, ie
  ijk  = threadidx%x
  ij   = modulo(threadidx%x-1,np*np)+1
  q  = blockidx%x
  ie = blockidx%y
  Qdp(ijk,q,np1,ie) = rspheremp(ij,ie) * Qdp(ijk,q,np1,ie)
end subroutine euler_hypervis_kernel_last



subroutine pack_exchange_unpack_stage(np1,hybrid,array_in,tl_in)
  use hybrid_mod, only : hybrid_t
  use schedule_mod, only : schedule_t, schedule, cycle_t
  use parallel_mod, only : abortmp, status, srequest, rrequest, mpireal_t, mpiinteger_t, iam
  use perf_mod, only: t_startf, t_stopf
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

  ierr = cudaStreamSynchronize(streams(1))
  call t_startf('CUDA PEU NEWER')
  do icycle=1,nRecvCycles
    pCycle => pSchedule%RecvCycle(icycle)
    source  = pCycle%source - 1
    length  = nlyr * pCycle%lengthP
    tag     = pCycle%tag
    call MPI_Irecv(recvbuf_h(1,1,icycle),length,MPIreal_t,source,tag,hybrid%par%comm,Rrequest(icycle),ierr)
  enddo
  if ( recv_external_nelem > 0 ) then
    blockdim6 = dim3( max(np,max_corner_elem_d) , nlev , 1 )
    griddim6  = dim3( qsize_d , recv_external_nelem , 1    )
    call edgeVpack_kernel_stage<<<griddim6,blockdim6,0,streams(1)>>>(edgebuf_d,array_in,putmapP_d,reverse_d,nbuf,0,1,nelemd,np1,recv_external_indices_d,tl_in); _CHECK(__LINE__)
  endif
  do icycle = 1 , nSendCycles
    pCycle => pSchedule%SendCycle(icycle)
    iptr   =  pCycle%ptrP
    ierr = cudaMemcpyAsync(sendbuf_h(1,1,icycle),edgebuf_d(1,iptr),size(sendbuf_h(1:nlyr,1:pCycle%lengthP,icycle)),cudaMemcpyDeviceToHost,streams(1)); _CHECK(__LINE__)
  enddo
  if ( recv_internal_nelem > 0 ) then
    blockdim6 = dim3( max(np,max_corner_elem_d) , nlev , 1 )
    griddim6  = dim3( qsize_d , recv_internal_nelem , 1    )
    call edgeVpack_kernel_stage<<<griddim6,blockdim6,0,streams(2)>>>(edgebuf_d,array_in,putmapP_d,reverse_d,nbuf,0,1,nelemd,np1,recv_internal_indices_d,tl_in); _CHECK(__LINE__)
    blockdim6 = dim3( max(np,max_corner_elem_d) , nlev , 1 )
    griddim6  = dim3( qsize_d , recv_internal_nelem , 1    )
    call edgeVunpack_kernel_stage<<<griddim6,blockdim6,0,streams(2)>>>(edgebuf_d,array_in,getmapP_d,nbuf,0,1,nelemd,np1,recv_internal_indices_d,tl_in); _CHECK(__LINE__)
  endif
  do icycle = 1 , nSendCycles
    pCycle => pSchedule%SendCycle(icycle)
    iptr   =  pCycle%ptrP
    ierr = cudaMemcpyAsync(sendbuf_h(1,1,icycle),edgebuf_d(1,iptr),size(sendbuf_h(1:nlyr,1:pCycle%lengthP,icycle)),cudaMemcpyDeviceToHost,streams(1)); _CHECK(__LINE__)
  enddo
  ierr = cudaStreamSynchronize(streams(1))
  do icycle = 1 , nSendCycles
    pCycle => pSchedule%SendCycle(icycle)
    dest   =  pCycle%dest - 1
    iptr   =  pCycle%ptrP
    length =  nlyr * pCycle%lengthP
    tag    =  pCycle%tag
    call MPI_Isend(sendbuf_h(1,1,icycle),length,MPIreal_t,dest,tag,hybrid%par%comm,Srequest(icycle),ierr)
  enddo
  call MPI_WaitAll(nRecvCycles,Rrequest,status,ierr)
  !When this cycle's MPI transfer is compliete, then call the D2H memcopy asynchronously
  do icycle = 1 , nRecvCycles
    pCycle => pSchedule%RecvCycle(icycle)
    iptr   =  pCycle%ptrP
    ierr = cudaMemcpyAsync(edgebuf_d(1,iptr),recvbuf_h(1,1,icycle),size(recvbuf_h(1:nlyr,1:pCycle%lengthP,icycle)),cudaMemcpyHostToDevice,streams(1)); _CHECK(__LINE__)
  enddo
  call MPI_WaitAll(nSendCycles,Srequest,status,ierr)
  if ( recv_external_nelem > 0 ) then
    blockdim6 = dim3( max(np,max_corner_elem_d) , nlev , 1 )
    griddim6  = dim3( qsize_d , recv_external_nelem , 1    )
    call edgeVunpack_kernel_stage<<<griddim6,blockdim6,0,streams(1)>>>(edgebuf_d,array_in,getmapP_d,nbuf,0,1,nelemd,np1,recv_external_indices_d,tl_in); _CHECK(__LINE__)
  endif
  ierr = cudaStreamSynchronize(streams(2))
  call t_stopf('CUDA PEU NEWER')

end subroutine pack_exchange_unpack_stage



attributes(global) subroutine edgeVpack_kernel_stage(edgebuf,v,putmapP,reverse,nbuf,kptr,nets,nete,nt,indices,tl_in)
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  implicit none
  real(kind=real_kind)   ,intent(in   ) :: v(np*np*nlev,qsize_d,tl_in,nets:nete)
  real(kind=real_kind)   ,intent(  out) :: edgebuf(nlev*qsize_d,nbuf)
  integer, value         ,intent(in   ) :: nbuf
  integer, value         ,intent(in   ) :: kptr
  integer, value         ,intent(in   ) :: nt
  integer                ,intent(in   ) :: putmapP(max_neigh_edges_d,nets:nete)
  logical                ,intent(in   ) :: reverse(max_neigh_edges_d,nets:nete)
  integer                ,intent(in   ) :: indices(nets:nete)
  integer, value         ,intent(in   ) :: tl_in
  integer, value         ,intent(in   ) :: nets
  integer, value         ,intent(in   ) :: nete
  integer :: i,k,ir,ll,kq,ie,q,ii,ik,koff,ind,ind_k,ind_ij
  real(kind=real_kind), shared :: v_s((np*np+1)*nlev)
  integer             , shared :: put_s(20)
  logical             , shared :: rev_s(4)

  i  = threadidx%x
  k  = threadidx%y
  q  = blockidx%x
  ie = indices(blockidx%y)
  kq = (q-1)*nlev+k + kptr
  ik = (k-1)*np+i
  koff = (k-1)*(np*np+1)

  !Efficiently load v, reverse, and putmapP into shared memory to re-use and to efficiently retreive DRAM memory.
  if (i <= np) then
    do ii = 1 , np  !Though this looping structure is strange, it allows the kernel to access contiguous np*nlev chunks of v.
      ind = (ii-1)*np*nlev+ik
      ind_k = (ind-1)/(np*np) + 1
      ind_ij = modulo(ind-1,np*np) + 1
      v_s((ind_k-1)*(np*np+1) + ind_ij) = v(ind,q,nt,ie)
    enddo
  endif
  if (ik <= 4                ) then
    rev_s(ik) = reverse(ik,ie)
  elseif (ik <= max_neigh_edges_d+4) then
    put_s(ik-4) = putmapP(ik-4,ie)
  endif
  call syncthreads()

  !Begin the packing process
  if (i <= np) then
    edgebuf(kq,put_s(west )+i) = v_s(koff+(i -1)*np+1 )
    edgebuf(kq,put_s(east )+i) = v_s(koff+(i -1)*np+np)
    edgebuf(kq,put_s(south)+i) = v_s(koff          +i )
    edgebuf(kq,put_s(north)+i) = v_s(koff+(np-1)*np+i )
  endif
  call syncthreads()
  if (i <= np) then
    ir = np-i+1
    if (rev_s(south)) edgebuf(kq,put_s(south)+ir) = v_s(koff          +i )
    if (rev_s(east )) edgebuf(kq,put_s(east )+ir) = v_s(koff+(i -1)*np+np)
    if (rev_s(north)) edgebuf(kq,put_s(north)+ir) = v_s(koff+(np-1)*np+i )
    if (rev_s(west )) edgebuf(kq,put_s(west )+ir) = v_s(koff+(i -1)*np+1 )
  endif
  if (i <= max_corner_elem_d) then
    ll = put_s(swest+0*max_corner_elem_d+i-1); if (ll /= -1) edgebuf(kq,ll+1) = v_s(koff          +1 )
    ll = put_s(swest+1*max_corner_elem_d+i-1); if (ll /= -1) edgebuf(kq,ll+1) = v_s(koff          +np)
    ll = put_s(swest+2*max_corner_elem_d+i-1); if (ll /= -1) edgebuf(kq,ll+1) = v_s(koff+(np-1)*np+1 )
    ll = put_s(swest+3*max_corner_elem_d+i-1); if (ll /= -1) edgebuf(kq,ll+1) = v_s(koff+(np-1)*np+np)
  endif
end subroutine edgeVpack_kernel_stage



attributes(global) subroutine edgeVunpack_kernel_stage(edgebuf,v,getmapP,nbuf,kptr,nets,nete,nt,recv_indices,tl_in)
  use edge_mod, only: EdgeDescriptor_t, EdgeBuffer_t
  use dimensions_mod, only : np, max_corner_elem
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  implicit none
  real(kind=real_kind)   ,intent(in   ) :: edgebuf(nlev*qsize_d,nbuf)
  real(kind=real_kind)   ,intent(inout) :: v(np*np*nlev,qsize_d,tl_in,nets:nete)
  integer                ,intent(in   ) :: getmapP(max_neigh_edges_d,nets:nete)
  integer, value         ,intent(in   ) :: nbuf
  integer, value         ,intent(in   ) :: kptr
  integer, value         ,intent(in   ) :: nets
  integer, value         ,intent(in   ) :: nete
  integer, value         ,intent(in   ) :: nt
  integer                ,intent(in   ) :: recv_indices(nets:nete)
  integer, value         ,intent(in   ) :: tl_in
  real(kind=real_kind), shared :: v_s((np*np+1)*nlev)
  integer             , shared :: get_s(20)
  integer :: i,k,ll,q,ie,kq,ik,koff,ind,ind_k,ind_ij,ii
  i  = threadidx%x
  k  = threadidx%y
  q  = blockidx%x
  ie = recv_indices(blockidx%y)
  kq = (q-1)*nlev+k + kptr
  ik = (k-1)*np+i
  koff = (k-1)*(np*np+1)

  !Efficiently load v, and getmapP into shared memory to re-use and to efficiently retreive DRAM memory.
  if (i <= np) then
    do ii = 1 , np  !Though this looping structure is strange, it allows the kernel to access contiguous np*nlev chunks of v.
      ind = (ii-1)*np*nlev+ik
      ind_k = (ind-1)/(np*np) + 1
      ind_ij = modulo(ind-1,np*np) + 1
      v_s((ind_k-1)*(np*np+1) + ind_ij) = v(ind,q,nt,ie)
    enddo
  endif
  if (ik <= max_neigh_edges_d) get_s(ik) = getmapP(ik,ie)
  call syncthreads()

  if (i <= np) then
    v_s(koff+(1 -1)*np+i ) = v_s(koff+(1 -1)*np+i ) + edgebuf(kptr+kq,get_s(south)+i)
    v_s(koff+(np-1)*np+i ) = v_s(koff+(np-1)*np+i ) + edgebuf(kptr+kq,get_s(north)+i)
  endif
  call syncthreads()
  if (i <= np) then
    v_s(koff+(i -1)*np+1 ) = v_s(koff+(i -1)*np+1 ) + edgebuf(kptr+kq,get_s(west )+i)
    v_s(koff+(i -1)*np+np) = v_s(koff+(i -1)*np+np) + edgebuf(kptr+kq,get_s(east )+i)
  endif
  call syncthreads()
  if (i <= max_corner_elem_d) then
    ll = get_s(swest+0*max_corner_elem_d+i-1); if(ll /= -1) v_s(koff+(1 -1)*np+1 ) = v_s(koff+(1 -1)*np+1 ) + edgebuf(kptr+kq,ll+1)
    ll = get_s(swest+1*max_corner_elem_d+i-1); if(ll /= -1) v_s(koff+(1 -1)*np+np) = v_s(koff+(1 -1)*np+np) + edgebuf(kptr+kq,ll+1)
    ll = get_s(swest+2*max_corner_elem_d+i-1); if(ll /= -1) v_s(koff+(np-1)*np+1 ) = v_s(koff+(np-1)*np+1 ) + edgebuf(kptr+kq,ll+1)
    ll = get_s(swest+3*max_corner_elem_d+i-1); if(ll /= -1) v_s(koff+(np-1)*np+np) = v_s(koff+(np-1)*np+np) + edgebuf(kptr+kq,ll+1)
  endif
  call syncthreads()

  !Efficiently store v_s into DRAM memory.
  if (i <= np) then
    do ii = 1 , np  !Though this looping structure is strange, it allows the kernel to access contiguous np*nlev chunks of v.
      ind = (ii-1)*np*nlev+ik
      ind_k = (ind-1)/(np*np) + 1
      ind_ij = modulo(ind-1,np*np) + 1
      v(ind,q,nt,ie) = v_s((ind_k-1)*(np*np+1) + ind_ij)
    enddo
  endif
end subroutine edgeVunpack_kernel_stage



attributes(global) subroutine hypervis_kernel1( Qdp , qtens , dp , dinv , variable_hyperviscosity , dpdiss_ave , spheremp , deriv_dvv , hyai , hybi , ps0 , nets , nete , dt , nt , nu_p )
  implicit none
  real(kind=real_kind), dimension(np,np,nlev,qsize_d,timelevels,nets:nete), intent(in   ) :: Qdp
  real(kind=real_kind), dimension(np,np,nlev,qsize_d           ,nets:nete), intent(  out) :: Qtens
  real(kind=real_kind), dimension(np,np,nlev                   ,nets:nete), intent(in   ) :: dp
  real(kind=real_kind), dimension(np,np,4                      ,nets:nete), intent(in   ) :: dinv
  real(kind=real_kind), dimension(np,np                                  ), intent(in   ) :: deriv_dvv
  real(kind=real_kind), dimension(      nlev+1                           ), intent(in   ) :: hyai
  real(kind=real_kind), dimension(      nlev+1                           ), intent(in   ) :: hybi
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: variable_hyperviscosity
  real(kind=real_kind), dimension(np,np,nlev                   ,nets:nete), intent(in   ) :: dpdiss_ave
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: spheremp
  real(kind=real_kind), value                                             , intent(in   ) :: dt , ps0 , nu_p
  integer, value                                                          , intent(in   ) :: nets , nete , nt
  real(kind=real_kind), dimension(np*np+1,numk_hyp,2), shared :: s
  real(kind=real_kind), dimension(np*np+1,4  ), shared :: dinv_s
  real(kind=real_kind), dimension(np,np), shared :: spheremp_s
  real(kind=real_kind), dimension(np,np), shared :: deriv_dvv_s
  real(kind=real_kind) :: dp0
  integer :: i, j, k, q, ie, kk, ks, ij

  ks = int(ceiling(dble(nlev)/numk_hyp))

  i  = modulo( threadidx%x-1    ,np)+1
  j  = modulo((threadidx%x-1)/np,np)+1
  kk = (threadidx%x-1)/(np*np)+1
  k  = modulo(blockidx%x-1,ks)*numk_hyp + kk
  q  = modulo((blockidx%x-1)/ks,qsize_d)+1
  ie =       ((blockidx%x-1)/ks)/qsize_d+1
  ij = (j-1)*np+i

  if (k  > nlev   ) return
  if (q  > qsize_d) return
  if (ie > nete   ) return

  if (kk == 1) then
    dinv_s(ij,:) = dinv(i,j,:,ie)
    spheremp_s(i,j) = spheremp(i,j,ie)
    deriv_dvv_s(i,j) = deriv_dvv(i,j)
  endif
  if (nu_p>0) then
    s(ij,kk,1) = dpdiss_ave(i,j,k,ie) * Qdp(i,j,k,q,nt,ie) / dp(i,j,k,ie)
  else
    dp0 = ( ( hyai(k+1) - hyai(k) )*ps0 + ( hybi(k+1) - hybi(k) )*ps0 )
    s(ij,kk,1) = dp0 * Qdp(i,j,k,q,nt,ie) / dp(i,j,k,ie)
  endif
  call syncthreads()
  qtens(i,j,k,q,ie) = laplace_sphere_wk(i,j,ie,kk,s,dinv_s,spheremp_s,variable_hyperviscosity,deriv_dvv_s,nets,nete)
end subroutine hypervis_kernel1



attributes(global) subroutine hypervis_kernel2( Qdp , qtens , dp , dinv , variable_hyperviscosity , spheremp , rspheremp , deriv_dvv , nu_q , nets , nete , dt , nt )
  implicit none
  real(kind=real_kind), dimension(np,np,nlev,qsize_d,timelevels,nets:nete), intent(inout) :: Qdp
  real(kind=real_kind), dimension(np,np,nlev,qsize_d           ,nets:nete), intent(in   ) :: Qtens
  real(kind=real_kind), dimension(np,np,nlev                   ,nets:nete), intent(in   ) :: dp
  real(kind=real_kind), dimension(np,np                                  ), intent(in   ) :: deriv_dvv
  real(kind=real_kind), dimension(np,np,4                      ,nets:nete), intent(in   ) :: dinv
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: variable_hyperviscosity
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: spheremp
  real(kind=real_kind), dimension(np,np                        ,nets:nete), intent(in   ) :: rspheremp
  real(kind=real_kind), value                                             , intent(in   ) :: dt, nu_q
  integer, value                                                          , intent(in   ) :: nets , nete , nt
  real(kind=real_kind), dimension(np*np+1,numk_hyp,2), shared :: s
  real(kind=real_kind), dimension(np*np+1,4  ), shared :: dinv_s
  real(kind=real_kind), dimension(np,np), shared :: spheremp_s
  real(kind=real_kind), dimension(np,np), shared :: deriv_dvv_s
  integer :: i, j, k, q, ie, kk, ks, ij

  ks = int(ceiling(dble(nlev)/numk_hyp))

  i  = modulo( threadidx%x-1    ,np)+1
  j  = modulo((threadidx%x-1)/np,np)+1
  kk = (threadidx%x-1)/(np*np)+1
  k  = modulo(blockidx%x-1,ks)*numk_hyp + kk
  q  = modulo((blockidx%x-1)/ks,qsize_d)+1
  ie =       ((blockidx%x-1)/ks)/qsize_d+1
  ij = (j-1)*np+i

  if (k  > nlev   ) return
  if (q  > qsize_d) return
  if (ie > nete   ) return

  if (kk == 1) then
    dinv_s(ij,:) = dinv(i,j,:,ie)
    spheremp_s(i,j) = spheremp(i,j,ie)
    deriv_dvv_s(i,j) = deriv_dvv(i,j)
  endif
  s(ij,kk,1) = rspheremp(i,j,ie)*qtens(i,j,k,q,ie)
  call syncthreads()
  Qdp(i,j,k,q,nt,ie) = Qdp(i,j,k,q,nt,ie)*spheremp_s(i,j)-dt*nu_q*laplace_sphere_wk(i,j,ie,kk,s,dinv_s,spheremp_s,variable_hyperviscosity,deriv_dvv_s,nets,nete)
end subroutine hypervis_kernel2



attributes(device) function laplace_sphere_wk(i,j,ie,k,s,dinv,spheremp,variable_hyperviscosity,deriv_dvv,nets,nete) result(lapl)
  implicit none
  integer,                                              intent(in) :: nets, nete, i, j, ie, k
  real(kind=real_kind), dimension(np*np+1,numk_hyp,2) , intent(inout) :: s
  real(kind=real_kind), dimension(np*np+1,2,2        ), intent(in) :: dinv
  real(kind=real_kind), dimension(np,np              ), intent(in) :: deriv_dvv
  real(kind=real_kind), dimension(np,np,nets:nete    ), intent(in) :: variable_hyperviscosity
  real(kind=real_kind), dimension(np,np              ), intent(in) :: spheremp
  real(kind=real_kind)                                             :: lapl
  integer :: l
  real(kind=real_kind) :: dsdx00 , dsdy00, tmp1, tmp2
  real(kind=real_kind), dimension(2) :: ds
  ds = gradient_sphere(i,j,ie,k,s,dinv,variable_hyperviscosity,deriv_dvv,nets,nete)
  lapl = divergence_sphere_wk(i,j,ie,k,ds,s,dinv,spheremp,deriv_dvv,nets,nete)
end function laplace_sphere_wk



attributes(device) function gradient_sphere(i,j,ie,k,s,dinv,variable_hyperviscosity,deriv_dvv,nets,nete) result(ds)
  implicit none
  integer,                                              intent(in) :: nets, nete, i, j, ie, k
  real(kind=real_kind), dimension(np*np+1,numk_hyp,2) , intent(in) :: s
  real(kind=real_kind), dimension(np*np+1,2,2        ), intent(in) :: dinv
  real(kind=real_kind), dimension(np,np              ), intent(in) :: deriv_dvv
  real(kind=real_kind), dimension(np,np,nets:nete    ), intent(in) :: variable_hyperviscosity
  real(kind=real_kind), dimension(2)                               :: ds
  integer :: l
  real(kind=real_kind) :: dsdx00 , dsdy00, tmp1, tmp2
  dsdx00 = 0.0d0
  dsdy00 = 0.0d0
  do l = 1 , np
    dsdx00 = dsdx00 + deriv_dvv(l,i)*s((j-1)*np+l,k,1)
    dsdy00 = dsdy00 + deriv_dvv(l,j)*s((l-1)*np+i,k,1)
  enddo
  ds(1) = ( dinv((j-1)*np+i,1,1)*dsdx00 + dinv((j-1)*np+i,2,1)*dsdy00 ) * rrearth_d * variable_hyperviscosity(i,j,ie)
  ds(2) = ( dinv((j-1)*np+i,1,2)*dsdx00 + dinv((j-1)*np+i,2,2)*dsdy00 ) * rrearth_d * variable_hyperviscosity(i,j,ie)
end function gradient_sphere



attributes(device) function divergence_sphere_wk(i,j,ie,k,tmp,s,dinv,spheremp,deriv_dvv,nets,nete) result(lapl)
  implicit none
  integer,                                              intent(in   ) :: nets, nete, i, j, ie, k
  real(kind=real_kind), dimension(2),                   intent(in   ) :: tmp
  real(kind=real_kind), dimension(np*np+1,numk_hyp,2) , intent(inout) :: s
  real(kind=real_kind), dimension(np*np+1,2,2        ), intent(in   ) :: dinv
  real(kind=real_kind), dimension(np,np              ), intent(in   ) :: deriv_dvv
  real(kind=real_kind), dimension(np,np              ), intent(in   ) :: spheremp
  real(kind=real_kind)                                                :: lapl
  integer :: l, ij
  ij = (j-1)*np+i
  s(ij,k,1) = ( dinv((j-1)*np+i,1,1)*tmp(1) + dinv((j-1)*np+i,1,2)*tmp(2) ) * spheremp(i,j)
  s(ij,k,2) = ( dinv((j-1)*np+i,2,1)*tmp(1) + dinv((j-1)*np+i,2,2)*tmp(2) ) * spheremp(i,j)
  call syncthreads()
  lapl = 0.0d0
  do l = 1 , np
    lapl = lapl - (s((j-1)*np+l,k,1)*deriv_dvv(i,l)+&
                   s((l-1)*np+i,k,2)*deriv_dvv(j,l)) * rrearth_d
  enddo
end function divergence_sphere_wk



#endif
end module cuda_mod


