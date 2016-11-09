
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module bndry_mod
  use bndry_mod_base, only: bndry_exchangeV, ghost_exchangeVfull, compute_ghost_corner_orientation, ghost_exchangeV, bndry_exchangeS, bndry_exchangeS_start, bndry_exchangeS_finish, sort_neighbor_buffer_mapping
  use parallel_mod, only : syncmp,parallel_t,abortmp,iam
  use edgetype_mod, only : Ghostbuffertr_t, Ghostbuffer3D_t,Edgebuffer_t,LongEdgebuffer_t
  use thread_mod, only : omp_in_parallel, omp_get_thread_num, omp_get_num_threads
  use kinds, only: real_kind
  implicit none
  private
  integer, parameter, private :: maxCycles = 20
  integer, parameter, private :: maxChunks = 64
  real(kind=real_kind), parameter, private :: chunk_denom = 1.e5

  type send_stager_t
    integer :: nUpdateHost, nSendComp
    logical :: updateHost(maxchunks), sendComp(maxchunks)
    integer :: beg(maxchunks), end(maxchunks), len(maxchunks), asyncid(maxchunks), tag(maxchunks), req(maxchunks)
  end type send_stager_t

  type recv_stager_t
    integer :: nRecvComp, nUpdateDev
    logical :: recvComp(maxchunks), updateDev(maxchunks)
    integer :: beg(maxchunks), end(maxchunks), len(maxchunks), asyncid(maxchunks), tag(maxchunks), req(maxchunks)
  end type recv_stager_t

  type(send_stager_t), private :: stg_send(maxCycles)
  type(recv_stager_t), private :: stg_recv(maxCycles)

  public :: bndry_exchangeV, ghost_exchangeVfull, compute_ghost_corner_orientation, ghost_exchangeV, bndry_exchangeS, bndry_exchangeS_start, bndry_exchangeS_finish, sort_neighbor_buffer_mapping
  public :: bndry_exchangeS_simple_overlap
  public :: bndry_exchangeV_timing
  public :: bndry_exchangeV_simple_overlap
  public :: bndry_exchangeV_finer_overlap

contains

  subroutine bndry_exchangeS_simple_overlap(hybrid,buffer)
    use hybrid_mod       , only : hybrid_t
    use kinds            , only : log_kind
    use edgetype_mod     , only : Edgebuffer_t
    use schedtype_mod    , only : schedule_t, cycle_t, schedule
    use parallel_mod     , only : abortmp, status, mpireal_t, mpiinteger_t, mpi_success
    use mpi              , only : MPI_REQUEST_NULL
    use openacc_utils_mod, only : update_host_async, update_device_async, copy_ondev_async, acc_async_test_wrap
    implicit none
    type (hybrid_t)           :: hybrid
    type (EdgeBuffer_t)       :: buffer
    type (Schedule_t),pointer :: pSchedule
    type (Cycle_t),pointer    :: pCycle
    integer                   :: dest,length,tag
    integer                   :: icycle,ierr
    integer                   :: iptr,source,nlyr
    integer                   :: nSendCycles,nRecvCycles
    integer                   :: errorcode,errorlen
    character*(80) errorstring
    integer        :: i, ithr
    integer :: nUpdateHost, nSendComp, nRecvComp
    logical :: updateHost(maxCycles), sendComp(maxCycles), recvComp(maxCycles)
    logical :: mpiflag
    !$OMP BARRIER
    if(hybrid%ithr == 0) then 
      !$acc wait
      pSchedule => Schedule(1)
      nlyr = buffer%nlyr
      nSendCycles = pSchedule%nSendCycles
      nRecvCycles = pSchedule%nRecvCycles
      if (max(nRecvCycles,nSendCycles) > maxCycles) then
        write(*,*) 'ERROR: Must increase maxCycles'
        stop
      endif
      nUpdateHost = 0
      nSendComp   = 0
      nRecvComp   = 0
      updateHost(1:nSendCycles) = .false.
      sendComp  (1:nSendCycles) = .false.
      recvComp  (1:nRecvCycles) = .false.
      buffer%Srequest(:) = MPI_REQUEST_NULL
      buffer%Rrequest(:) = MPI_REQUEST_NULL
      !==================================================
      !  Post the Receives 
      !==================================================
      do icycle=1,nRecvCycles
        pCycle => pSchedule%RecvCycle(icycle)
        source =  pCycle%source - 1
        length =  nlyr * pCycle%lengthS
!       tag    =  pCycle%tag
        tag    =  buffer%tag
        iptr   =  pCycle%ptrS
        call MPI_Irecv(buffer%receive(1+nlyr*(iptr-1)),length,MPIreal_t,source,tag,hybrid%par%comm,buffer%Rrequest(icycle),ierr)
        if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,'bndry_exchangeV: Error after call to MPI_Irecv: ',errorstring
        endif
      enddo    ! icycle

      !Copy internal data to receive buffer
      do ithr = 1 , hybrid%NThreads
        iptr   = buffer%moveptr(ithr)
        length = buffer%moveLength(ithr)
        if(length>0) call copy_ondev_async(buffer%receive(iptr),buffer%buf(iptr),length,maxCycles*2+ithr)
      enddo

      !Launch PCI-e copies
      do icycle = 1 , nSendCycles
        pCycle => pSchedule%SendCycle(icycle)
        iptr   =  pCycle%ptrS
        if (pCycle%lengthS > 0) call update_host_async(buffer%buf(1+nlyr*(iptr-1)),nlyr*pCycle%lengthS,icycle)
      enddo
      !Initiate polling loop for MPI_Isend and PCI-e returns after data is received
      do while (nRecvComp < nRecvCycles .or. nSendComp < nSendCycles .or. nUpdateHost < nSendCycles)
        !If there are host updates yet pending, test to see if each cycle is done. If a cycle is done, so the mpi_isend
        if (nUpdateHost < nSendCycles) then
          do icycle = 1 , nSendCycles
            if (.not. updateHost(icycle)) then
              if (acc_async_test_wrap(icycle)) then
                pCycle => pSchedule%SendCycle(icycle)
                dest   =  pCycle%dest - 1
                length =  nlyr * pCycle%lengthS
!               tag    =  pCycle%tag
                tag    =  buffer%tag
                iptr   =  pCycle%ptrS
                call MPI_Isend(buffer%buf(1+nlyr*(iptr-1)),length,MPIreal_t,dest,tag,hybrid%par%comm,buffer%Srequest(icycle),ierr)
                updateHost(icycle) = .true.
                nUpdateHost = nUpdateHost + 1
              endif
            endif
          enddo
        endif
        !If there are mpi_isend's still pending, test to see if each cycle is completed. This is for bookkeeping. I cannot copy from receive to buf until all sends are completed
        if (nSendComp < nSendCycles) then
          do icycle = 1 , nSendCycles
            if (updateHost(icycle) .and. (.not. sendComp(icycle))) then
              call MPI_Test(buffer%Srequest(icycle),mpiflag,status(:,icycle),ierr)
              if (mpiflag) then
                sendComp(icycle) = .true.
                nSendComp = nSendComp + 1
              endif
            endif
          enddo
        endif
        !if there are mpi_irecv's still pending, test to see if each cycle is completed. If it is, the send receive buffer to device
        if (nRecvComp < nRecvCycles) then
          do icycle = 1 , nRecvCycles
            if (.not. recvComp(icycle)) then
              call MPI_Test(buffer%Rrequest(icycle),mpiflag,status(:,icycle),ierr)
              if (mpiflag) then
                pCycle => pSchedule%RecvCycle(icycle)
                iptr   =  pCycle%ptrS
                call update_device_async(buffer%receive(1+nlyr*(iptr-1)),nlyr*pCycle%lengthS,maxCycles+icycle)
                recvComp(icycle) = .true.
                nRecvComp = nRecvComp + 1
              endif
            endif
          enddo
        endif
      enddo
      !$acc wait
    endif  ! if (hybrid%ithr == 0)
    !$OMP BARRIER
  end subroutine bndry_exchangeS_simple_overlap

  subroutine bndry_exchangeV_timing(hybrid,buffer)
    use hybrid_mod       , only : hybrid_t
    use kinds            , only : log_kind
    use edgetype_mod     , only : Edgebuffer_t
    use schedtype_mod    , only : schedule_t, cycle_t, schedule
    use parallel_mod     , only : abortmp, status, srequest, rrequest, mpireal_t, mpiinteger_t, mpi_success
    use mpi              , only : MPI_REQUEST_NULL
    use openacc_utils_mod, only : update_host_async, update_device_async, copy_ondev_async, acc_async_test_wrap
    use perf_mod         , only : t_startf, t_stopf
    implicit none
    type (hybrid_t)           :: hybrid
    type (EdgeBuffer_t)       :: buffer
    type (Schedule_t),pointer :: pSchedule
    type (Cycle_t),pointer    :: pCycle
    integer                   :: dest,length,tag
    integer                   :: icycle,ierr
    integer                   :: iptr,source,nlyr
    integer                   :: nSendCycles,nRecvCycles
    integer                   :: errorcode,errorlen
    character*(80) errorstring
    integer        :: i
    integer :: nUpdateHost, nSendComp, nRecvComp, nUpdateDev
    logical :: updateHost(maxCycles), sendComp(maxCycles), recvComp(maxCycles), updateDev(maxCycles)
    logical :: mpiflag
    integer :: ithr
    !$OMP BARRIER
    if(hybrid%ithr == 0) then 
      !$acc wait
      call t_startf('bndry_timing')
      pSchedule => Schedule(1)
      nlyr = buffer%nlyr
      nSendCycles = pSchedule%nSendCycles
      nRecvCycles = pSchedule%nRecvCycles
      if (max(nRecvCycles,nSendCycles) > maxCycles) then
        write(*,*) 'ERROR: Must increase maxCycles'
        stop
      endif
      nUpdateHost = 0
      nSendComp   = 0
      nRecvComp   = 0
      nUpdateDev  = 0
      updateHost(1:nSendCycles) = .false.
      sendComp  (1:nSendCycles) = .false.
      recvComp  (1:nRecvCycles) = .false.
      updateDev (1:nRecvCycles) = .false.
      Srequest(:) = MPI_REQUEST_NULL
      Rrequest(:) = MPI_REQUEST_NULL
      !==================================================
      !  Post the Receives 
      !==================================================
      do icycle=1,nRecvCycles
        pCycle => pSchedule%RecvCycle(icycle)
        source =  pCycle%source - 1
        length =  nlyr * pCycle%lengthP
        tag    =  pCycle%tag
        iptr   =  pCycle%ptrP
        call MPI_Irecv(buffer%receive(1+nlyr*(iptr-1)),length,MPIreal_t,source,tag,hybrid%par%comm,Rrequest(icycle),ierr)
        if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,'bndry_exchangeV: Error after call to MPI_Irecv: ',errorstring
        endif
      enddo    ! icycle

      !Copy internal data to receive buffer
      call t_startf('bndry_timing_internal_on_dev_copy')
      do ithr = 1 , hybrid%NThreads
        iptr   = buffer%moveptr(ithr)
        length = buffer%moveLength(ithr)
        if(length>0) call copy_ondev_async(buffer%receive(iptr),buffer%buf(iptr),length,1)
      enddo
      !$acc wait(1)
      call t_stopf('bndry_timing_internal_on_dev_copy')

      !Launch PCI-e copies
      call t_startf('bndry_timing_pcie_d2h')
      do icycle = 1 , nSendCycles
        pCycle => pSchedule%SendCycle(icycle)
        iptr   =  pCycle%ptrP
        if (pCycle%lengthP > 0) call update_host_async(buffer%buf(1+nlyr*(iptr-1)),nlyr*pCycle%lengthP,1)
      enddo
      !$acc wait(1)
      call t_stopf('bndry_timing_pcie_d2h')

      do icycle = 1 , nSendCycles
        pCycle => pSchedule%SendCycle(icycle)
        dest   =  pCycle%dest - 1
        length =  nlyr * pCycle%lengthP
        tag    =  pCycle%tag
        iptr   =  pCycle%ptrP
        call MPI_Isend(buffer%buf(1+nlyr*(iptr-1)),length,MPIreal_t,dest,tag,hybrid%par%comm,Srequest(icycle),ierr)
      enddo

      call MPI_Waitall(nSendCycles,Srequest,status,ierr)
      call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)

      call t_startf('bndry_timing_pcie_h2d')
      do icycle = 1 , nRecvCycles
        pCycle => pSchedule%RecvCycle(icycle)
        iptr   =  pCycle%ptrP
        call update_device_async(buffer%receive(1+nlyr*(iptr-1)),nlyr*pCycle%lengthP,1)
        recvComp(icycle) = .true.
        nRecvComp = nRecvComp + 1
      enddo
      !$acc wait(1)
      call t_stopf('bndry_timing_pcie_h2d')
      call t_stopf('bndry_timing')
    endif  ! if (hybrid%ithr == 0)
    !$OMP BARRIER
  end subroutine bndry_exchangeV_timing

  subroutine bndry_exchangeV_simple_overlap(hybrid,buffer)
    use hybrid_mod       , only : hybrid_t
    use kinds            , only : log_kind
    use edgetype_mod     , only : Edgebuffer_t
    use schedtype_mod    , only : schedule_t, cycle_t, schedule
    use parallel_mod     , only : abortmp, status, srequest, rrequest, mpireal_t, mpiinteger_t, mpi_success
    use mpi              , only : MPI_REQUEST_NULL
    use openacc_utils_mod, only : update_host_async, update_device_async, copy_ondev_async, acc_async_test_wrap
    implicit none
    type (hybrid_t)           :: hybrid
    type (EdgeBuffer_t)       :: buffer
    type (Schedule_t),pointer :: pSchedule
    type (Cycle_t),pointer    :: pCycle
    integer                   :: dest,length,tag
    integer                   :: icycle,ierr
    integer                   :: iptr,source,nlyr
    integer                   :: nSendCycles,nRecvCycles
    integer                   :: errorcode,errorlen
    character*(80) errorstring
    integer        :: i
    integer :: nUpdateHost, nSendComp, nRecvComp, nUpdateDev
    logical :: updateHost(maxCycles), sendComp(maxCycles), recvComp(maxCycles), updateDev(maxCycles)
    logical :: mpiflag
    integer :: ithr
    !$OMP BARRIER
    if(hybrid%ithr == 0) then 
      !$acc wait
      pSchedule => Schedule(1)
      nlyr = buffer%nlyr
      nSendCycles = pSchedule%nSendCycles
      nRecvCycles = pSchedule%nRecvCycles
      if (max(nRecvCycles,nSendCycles) > maxCycles) then
        write(*,*) 'ERROR: Must increase maxCycles'
        stop
      endif
      nUpdateHost = 0
      nSendComp   = 0
      nRecvComp   = 0
      nUpdateDev  = 0
      updateHost(1:nSendCycles) = .false.
      sendComp  (1:nSendCycles) = .false.
      recvComp  (1:nRecvCycles) = .false.
      updateDev (1:nRecvCycles) = .false.
      Srequest(:) = MPI_REQUEST_NULL
      Rrequest(:) = MPI_REQUEST_NULL
      !==================================================
      !  Post the Receives 
      !==================================================
      do icycle=1,nRecvCycles
        pCycle => pSchedule%RecvCycle(icycle)
        source =  pCycle%source - 1
        length =  nlyr * pCycle%lengthP
        tag    =  pCycle%tag
        iptr   =  pCycle%ptrP
        call MPI_Irecv(buffer%receive(1+nlyr*(iptr-1)),length,MPIreal_t,source,tag,hybrid%par%comm,Rrequest(icycle),ierr)
        if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,'bndry_exchangeV: Error after call to MPI_Irecv: ',errorstring
        endif
      enddo    ! icycle

      !Copy internal data to receive buffer
      do ithr = 1 , hybrid%NThreads
        iptr   = buffer%moveptr(ithr)
        length = buffer%moveLength(ithr)
        if(length>0) call copy_ondev_async(buffer%receive(iptr),buffer%buf(iptr),length,maxCycles*2+ithr)
      enddo

      !Launch PCI-e copies
      do icycle = 1 , nSendCycles
        pCycle => pSchedule%SendCycle(icycle)
        iptr   =  pCycle%ptrP
        if (pCycle%lengthP > 0) call update_host_async(buffer%buf(1+nlyr*(iptr-1)),nlyr*pCycle%lengthP,icycle)
      enddo
      !Initiate polling loop for MPI_Isend and PCI-e returns after data is received
      do while (nRecvComp < nRecvCycles .or. nSendComp < nSendCycles .or. nUpdateHost < nSendCycles)
        !If there are host updates yet pending, test to see if each cycle is done. If a cycle is done, so the mpi_isend
        if (nUpdateHost < nSendCycles) then
          do icycle = 1 , nSendCycles
            if (.not. updateHost(icycle)) then
              if (acc_async_test_wrap(icycle)) then
                pCycle => pSchedule%SendCycle(icycle)
                dest   =  pCycle%dest - 1
                length =  nlyr * pCycle%lengthP
                tag    =  pCycle%tag
                iptr   =  pCycle%ptrP
                call MPI_Isend(buffer%buf(1+nlyr*(iptr-1)),length,MPIreal_t,dest,tag,hybrid%par%comm,Srequest(icycle),ierr)
                updateHost(icycle) = .true.
                nUpdateHost = nUpdateHost + 1
              endif
            endif
          enddo
        endif
        !If there are mpi_isend's still pending, test to see if each cycle is completed. This is for bookkeeping. I cannot copy from receive to buf until all sends are completed
        if (nSendComp < nSendCycles) then
          do icycle = 1 , nSendCycles
            if (updateHost(icycle) .and. (.not. sendComp(icycle))) then
              call MPI_Test(Srequest(icycle),mpiflag,status(:,icycle),ierr)
              if (mpiflag) then
                sendComp(icycle) = .true.
                nSendComp = nSendComp + 1
              endif
            endif
          enddo
        endif
        !if there are mpi_irecv's still pending, test to see if each cycle is completed. If it is, the send receive buffer to device
        if (nRecvComp < nRecvCycles) then
          do icycle = 1 , nRecvCycles
            if (.not. recvComp(icycle)) then
              call MPI_Test(Rrequest(icycle),mpiflag,status(:,icycle),ierr)
              if (mpiflag) then
                pCycle => pSchedule%RecvCycle(icycle)
                iptr   =  pCycle%ptrP
                call update_device_async(buffer%receive(1+nlyr*(iptr-1)),nlyr*pCycle%lengthP,maxCycles+icycle)
                recvComp(icycle) = .true.
                nRecvComp = nRecvComp + 1
              endif
            endif
          enddo
        endif
      enddo
      !$acc wait
    endif  ! if (hybrid%ithr == 0)
    !$OMP BARRIER
  end subroutine bndry_exchangeV_simple_overlap

  subroutine bndry_exchangeV_finer_overlap(hybrid,buffer)
    use hybrid_mod       , only : hybrid_t
    use kinds            , only : log_kind
    use edgetype_mod     , only : Edgebuffer_t
    use schedtype_mod    , only : schedule_t, cycle_t, schedule
    use parallel_mod     , only : abortmp, status, srequest, rrequest, mpireal_t, mpiinteger_t, mpi_success
    use mpi              , only : MPI_REQUEST_NULL
    implicit none
    type (hybrid_t)           :: hybrid
    type (EdgeBuffer_t)       :: buffer
    type (Schedule_t),pointer :: pSchedule
    type (Cycle_t),pointer    :: pCycle
    integer                   :: dest,length,tag
    integer                   :: icycle,ierr
    integer                   :: iptr,source,nlyr
    integer                   :: nSendCycles,nRecvCycles
    integer                   :: errorcode,errorlen
    character*(80) errorstring
    integer        :: i
    integer :: nSendComp, nRecvComp
    logical :: sendComp(maxCycles), recvComp(maxCycles)
    logical :: dummy
    integer :: nchunks
    
    !$OMP BARRIER
    if(hybrid%ithr == 0) then 
      !$acc wait
      pSchedule => Schedule(1)
      nlyr = buffer%nlyr
      nSendCycles = pSchedule%nSendCycles
      nRecvCycles = pSchedule%nRecvCycles
      if (max(nRecvCycles,nSendCycles) > maxCycles) then
        write(*,*) 'ERROR: Must increase maxCycles'
        stop
      endif
      nSendComp   = 0
      nRecvComp   = 0
      sendComp  (1:nSendCycles) = .false.
      recvComp  (1:nRecvCycles) = .false.

      !  Post the Receives 
      do icycle=1,nRecvCycles
        pCycle => pSchedule%RecvCycle(icycle)
        source =  pCycle%source - 1
        length =  nlyr * pCycle%lengthP
        tag    =  pCycle%tag
        iptr   =  pCycle%ptrP
        nchunks = int(ceiling(length/chunk_denom))
        dummy = mpi_irecv_openacc_stage(buffer%receive(1+nlyr*(iptr-1)), length, source, tag, hybrid%par%comm, ierr, nchunks, maxCycles+icycle, .true. , icycle , .false. , buffer%buf(1+nlyr*(iptr-1)) )
      enddo    ! icycle

      !Launch PCI-e copies
      do icycle = 1 , nSendCycles
        pCycle => pSchedule%SendCycle(icycle)
        iptr   =  pCycle%ptrP
        length =  nlyr * pCycle%lengthP
        dest   =  pCycle%dest - 1
        tag    =  pCycle%tag
        nchunks = int(ceiling(length/chunk_denom))
        dummy = mpi_isend_openacc_stage(buffer%buf(1+nlyr*(iptr-1)), length, dest, tag, hybrid%par%comm, ierr, nchunks, icycle, .true. , icycle)
      enddo

      !Initiate polling loop for MPI_Isend and PCI-e returns after data is received
      do while (nRecvComp < nRecvCycles .or. nSendComp < nSendCycles)
        !If there are mpi_isend's still pending, test to see if each cycle is completed. This is for bookkeeping. I cannot copy from receive to buf until all sends are completed
        if (nSendComp < nSendCycles) then
          do icycle = 1 , nSendCycles
            if (.not. sendComp(icycle)) then
              pCycle => pSchedule%SendCycle(icycle)
              iptr   =  pCycle%ptrP
              length =  nlyr * pCycle%lengthP
              dest   =  pCycle%dest - 1
              tag    =  pCycle%tag
              nchunks = int(ceiling(length/chunk_denom))
              if (mpi_isend_openacc_stage(buffer%buf(1+nlyr*(iptr-1)), length, dest, tag, hybrid%par%comm, ierr, nchunks, icycle, .false. , icycle )) then
                sendComp(icycle) = .true.
                nSendComp = nSendComp + 1
              endif
            endif
          enddo
        endif
        !if there are mpi_irecv's still pending, test to see if each cycle is completed. If it is, the send receive buffer to device
        if (nRecvComp < nRecvCycles) then
          do icycle = 1 , nRecvCycles
            if (.not. recvComp(icycle)) then
              pCycle => pSchedule%RecvCycle(icycle)
              source =  pCycle%source - 1
              length =  nlyr * pCycle%lengthP
              tag    =  pCycle%tag
              iptr   =  pCycle%ptrP
              nchunks = int(ceiling(length/chunk_denom))
              if (mpi_irecv_openacc_stage(buffer%receive(1+nlyr*(iptr-1)), length, source, tag, hybrid%par%comm, ierr, nchunks, maxCycles+icycle, .false. , icycle , &
                                        & nSendComp == nSendCycles , buffer%buf(1+nlyr*(iptr-1)) )) then
                recvComp(icycle) = .true.
                nRecvComp = nRecvComp + 1
              endif
            endif
          enddo
        endif
      enddo
      !$acc wait
    endif  ! if (hybrid%ithr == 0)
    !$OMP BARRIER
  end subroutine bndry_exchangeV_finer_overlap

  !Consider this to be the inside of a polling loop for sending something over MPI that currently resides on-device in OpenACC
  !You'll call this function over and over and over again until it returns .true.
  !On the first call, set reset=.true., and it will treat it as the first call and perform the
  !PCI-e update host calls. After that, use reset=.false. and make sure all input variables are the same
  !If an error occurs, I'll try to pass it along with ierror.
  !When the isends are completed, the function will return .true. Until then, it will return .false.
  function mpi_isend_openacc_stage(buf, count, dest, tag_root, comm, ierror, nchunks, async_root, reset, myid)   result(finished)
    use mpi              , only: MPI_REQUEST_NULL, MPI_STATUS_SIZE, mpi_test
    use parallel_mod     , only: mpireal_t
    use openacc_utils_mod, only: update_host_async, acc_async_test_wrap
    implicit none
    integer             , intent(in   ) :: count        !number of elements in buffer
    real(kind=real_kind), intent(in   ) :: buf(count)   !buffer from which to send data
    integer             , intent(in   ) :: dest         !the MPI rank I'm sending data to
    integer             , intent(in   ) :: tag_root     !tag of the original send (I'll alter this for internal mpi_isend calls)
    integer             , intent(in   ) :: comm         !Communicator to use
    integer             , intent(  out) :: ierror       !Try to return errors to the user
    integer             , intent(in   ) :: nchunks      !Number of chunks to break data into (number of mpi_isend calls)
    integer             , intent(in   ) :: async_root   !Root asyncid to use
    logical             , intent(in   ) :: reset        !First call in a given mpi_isend should always be a reset
    integer             , intent(in   ) :: myid         !the unique id for this series of mpi_isend's
    logical                             :: finished     !Return: Let user know if the mpi_isend's are all completed.
    integer :: i
    logical :: mpiflag
    integer :: status(MPI_STATUS_SIZE)
    if (reset) then
      stg_send(myid)%nUpdateHost = 0
      stg_send(myid)%nSendComp = 0
      do i = 1 , nchunks
        stg_send(myid)%updateHost(i) = .false.
        stg_send(myid)%sendComp(i) = .false.
        stg_send(myid)%req(i) = MPI_REQUEST_NULL
        stg_send(myid)%beg(i) = count*(i-1)/nchunks+1
        stg_send(myid)%end(i) = count*(i  )/nchunks
        stg_send(myid)%len(i) = stg_send(myid)%end(i) - stg_send(myid)%beg(i) + 1
        stg_send(myid)%asyncid(i) = async_root*maxchunks+i
        stg_send(myid)%tag(i) = tag_root*maxchunks+i
        call update_host_async( buf( stg_send(myid)%beg(i) ) , stg_send(myid)%len(i) , stg_send(myid)%asyncid(i) )
      enddo
    endif
    !If there are host updates yet pending, test to see if each chunk is done. If achunke is done, do the mpi_isend
    if (stg_send(myid)%nUpdateHost < nchunks) then
      do i = 1 , nchunks
        if (.not. stg_send(myid)%updateHost(i)) then
          if (acc_async_test_wrap(stg_send(myid)%asyncid(i))) then
            call MPI_Isend(buf(stg_send(myid)%beg(i)),stg_send(myid)%len(i),MPIreal_t,dest,stg_send(myid)%tag(i),comm,stg_send(myid)%req(i),ierror)
            stg_send(myid)%updateHost(i) = .true.
            stg_send(myid)%nUpdateHost = stg_send(myid)%nUpdateHost + 1
          endif
        endif
      enddo
    endif
    !If there are mpi_isend's still pending, test to see if each chunk is completed.
    if (stg_send(myid)%nSendComp < nchunks) then
      do i = 1 , nchunks
        if ( stg_send(myid)%updateHost(i) .and. (.not. stg_send(myid)%sendComp(i)) ) then
          call MPI_Test(stg_send(myid)%req(i),mpiflag,status,ierror)
          if (mpiflag) then
            stg_send(myid)%sendComp(i) = .true.
            stg_send(myid)%nSendComp = stg_send(myid)%nSendComp + 1
          endif
        endif
      enddo
    endif
    finished = .false.
    if (stg_send(myid)%nSendComp == nchunks) finished = .true.
  end function mpi_isend_openacc_stage

  function mpi_irecv_openacc_stage(buf, count, source, tag_root, comm, ierror, nchunks, async_root, reset, myid, sends_done, sendbuf)   result(finished)
    use mpi              , only: MPI_REQUEST_NULL, MPI_STATUS_SIZE, mpi_test
    use parallel_mod     , only: mpireal_t
    use openacc_utils_mod, only: update_device_async, copy_ondev_async
    implicit none
    integer             , intent(in   ) :: count        !number of elements in buffer
    real(kind=real_kind), intent(in   ) :: buf(count)   !buffer in which to receive data
    integer             , intent(in   ) :: source       !the MPI rank I'm receiving data from
    integer             , intent(in   ) :: tag_root     !tag of the original send (I'll alter this for internal mpi_isend calls)
    integer             , intent(in   ) :: comm         !Communicator to use
    integer             , intent(  out) :: ierror       !Try to return errors to the user
    integer             , intent(in   ) :: nchunks      !Number of chunks to break data into (number of mpi_isend calls)
    integer             , intent(in   ) :: async_root   !Root asyncid to use
    logical             , intent(in   ) :: reset        !First call in a given mpi_isend should always be a reset
    integer             , intent(in   ) :: myid         !the unique id for this series of mpi_irecv's
    logical             , intent(in   ) :: sends_done   !the unique id for this series of mpi_irecv's
    real(kind=real_kind), intent(inout) :: sendbuf(count)   !buffer from which to send data
    logical                             :: finished     !Return: Let user know if the mpi_isend's are all completed.
    integer :: i
    logical :: mpiflag
    integer :: status(MPI_STATUS_SIZE)
    if (reset) then
      stg_recv(myid)%nRecvComp = 0
      stg_recv(myid)%nUpdateDev = 0
      do i = 1 , nchunks
        stg_recv(myid)%recvComp(i) = .false.
        stg_recv(myid)%updateDev(i) = .false.
        stg_recv(myid)%req(i) = MPI_REQUEST_NULL
        stg_recv(myid)%beg(i) = count*(i-1)/nchunks+1
        stg_recv(myid)%end(i) = count*(i  )/nchunks
        stg_recv(myid)%len(i) = stg_recv(myid)%end(i) - stg_recv(myid)%beg(i) + 1
        stg_recv(myid)%asyncid(i) = async_root*maxchunks+i
        stg_recv(myid)%tag(i) = tag_root*maxchunks+i
        call MPI_Irecv(buf(stg_recv(myid)%beg(i)),stg_recv(myid)%len(i),MPIreal_t,source,stg_recv(myid)%tag(i),comm,stg_recv(myid)%req(i),ierror)
      enddo
    endif
    !if there are mpi_irecv's still pending, test to see if each cycle is completed. If it is, the send receive buffer to device
    if (stg_recv(myid)%nRecvComp < nchunks) then
      do i = 1 , nchunks
        if (.not. stg_recv(myid)%recvComp(i)) then
          call MPI_Test(stg_recv(myid)%req(i),mpiflag,status,ierror)
          if (mpiflag) then
            call update_device_async(buf(stg_recv(myid)%beg(i)),stg_recv(myid)%len(i),stg_recv(myid)%asyncid(i))
            stg_recv(myid)%recvComp(i) = .true.
            stg_recv(myid)%nRecvComp = stg_recv(myid)%nRecvComp + 1
          endif
        endif
      enddo
    endif
    !if all sends are completed,
    if (sends_done) then
      do i = 1 , nchunks
        if ( stg_recv(myid)%recvComp(i) .and. ( .not. stg_recv(myid)%updateDev(i) ) ) then
          call copy_ondev_async(sendbuf(stg_recv(myid)%beg(i)),buf(stg_recv(myid)%beg(i)),stg_recv(myid)%len(i),stg_recv(myid)%asyncid(i))
          stg_recv(myid)%updateDev(i) = .true.
          stg_recv(myid)%nUpdateDev = stg_recv(myid)%nUpdateDev + 1
        endif
      enddo
    endif
    finished = .false.
    if (stg_recv(myid)%nUpdateDev == nchunks) finished = .true.
  end function mpi_irecv_openacc_stage

end module bndry_mod

