
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
  public :: bndry_exchangeV_simple_overlap
  public :: bndry_exchangeV_minimize_latency

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
      do ithr = 1 , hybrid%hthreads
        iptr   = nlyr*buffer%moveptr0(ithr) + 1
        length = nlyr*buffer%moveLength(ithr)
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
      do ithr = 1 , hybrid%hthreads
        iptr   = nlyr*buffer%moveptr0(ithr) + 1
        length = nlyr*buffer%moveLength(ithr)
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

  subroutine bndry_exchangeV_minimize_latency(hybrid,buffer,asyncid)
    use hybrid_mod       , only : hybrid_t
    use kinds            , only : log_kind
    use edgetype_mod     , only : Edgebuffer_t
    use schedtype_mod    , only : schedule_t, cycle_t, schedule
    use parallel_mod     , only : abortmp, status, srequest, rrequest, mpireal_t, mpiinteger_t, mpi_success
    use mpi              , only : MPI_REQUEST_NULL
    use openacc_utils_mod, only : update_host_async, update_device_async, copy_ondev_async, acc_async_test_wrap
    implicit none
    type (hybrid_t)     :: hybrid
    type (EdgeBuffer_t) :: buffer
    integer, intent(in) :: asyncid
    type (Schedule_t),pointer :: pSchedule
    type (Cycle_t),pointer    :: pCycle
    integer                   :: dest,length,tag,icycle,ierr,iptr,source,nlyr,nSendCycles,nRecvCycles,errorcode,errorlen,nbuf
    character*(80) errorstring
    !$OMP BARRIER
    !$omp master
    nbuf = buffer%nbuf
    !Assume latency is dominating, so reduce transactions as much as possible by transferring more data
    !Also keep asynchronous packed together to avoid launch overheads
    call copy_ondev_async(buffer%receive,buffer%buf,product(shape(buffer%buf)),asyncid)
    !$acc update host(buffer%receive) async(asyncid)
    !$acc update host(buffer%buf) async(asyncid)

    pSchedule => Schedule(1)
    nlyr = buffer%nlyr
    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles
    Srequest(:) = MPI_REQUEST_NULL
    Rrequest(:) = MPI_REQUEST_NULL

    !Must have the receive buffer transfer completed before the MPI receives start. Otherwise MPI data could be overwritten
    !$acc wait(asyncid)
    !==================================================
    !  Post the Receives
    !==================================================
    do icycle=1,nRecvCycles
       pCycle => pSchedule%RecvCycle(icycle)
       source =  pCycle%source - 1
       length =  nlyr * pCycle%lengthP
       tag    =  buffer%tag
       iptr   =  nlyr * (pCycle%ptrP -1) + 1
       call MPI_Irecv(buffer%receive(iptr),length,MPIreal_t,source,tag,hybrid%par%comm,buffer%Rrequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,'bndry_exchangeV: Error after call to MPI_Irecv: ',errorstring
       endif
    end do    ! icycle
    do icycle=1,nSendCycles
       pCycle => pSchedule%SendCycle(icycle)
       dest   =  pCycle%dest - 1
       length =  nlyr * pCycle%lengthP
       tag    = buffer%tag
       iptr   = nlyr * (pCycle%ptrP - 1) + 1
       call MPI_Isend(buffer%buf(iptr),length,MPIreal_t,dest,tag,hybrid%par%comm,buffer%Srequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,'bndry_exchangeV: Error after call to MPI_Isend: ',errorstring
       endif
    end do    ! icycle
    call MPI_Waitall(nSendCycles,buffer%Srequest,buffer%status,ierr)
    !$acc update device(buffer%receive) async(asyncid)
    call MPI_Waitall(nRecvCycles,buffer%Rrequest,buffer%status,ierr)
    !$omp end master
    !$OMP BARRIER
  end subroutine bndry_exchangeV_minimize_latency

end module bndry_mod
