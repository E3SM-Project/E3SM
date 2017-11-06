
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module bndry_mod
  use bndry_mod_base, only: bndry_exchangeV, ghost_exchangeVfull, compute_ghost_corner_orientation, bndry_exchangeS, bndry_exchangeS_start, bndry_exchangeS_finish, sort_neighbor_buffer_mapping
  use parallel_mod, only : syncmp,parallel_t,abortmp,iam
  use edgetype_mod, only : Ghostbuffer3D_t,Edgebuffer_t,LongEdgebuffer_t
  use kinds, only: real_kind
  implicit none
  private
  integer, parameter, private :: maxCycles = 20

  public :: bndry_exchangeV, ghost_exchangeVfull, compute_ghost_corner_orientation, bndry_exchangeS, bndry_exchangeS_start, bndry_exchangeS_finish, sort_neighbor_buffer_mapping
  public :: bndry_exchangeS_simple_overlap
  public :: bndry_exchangeV_timing
  public :: bndry_exchangeV_simple_overlap

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

end module bndry_mod
