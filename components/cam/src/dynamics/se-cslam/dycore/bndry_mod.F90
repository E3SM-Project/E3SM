module bndry_mod
  use shr_kind_mod,   only: r8=>shr_kind_r8, i8=>shr_kind_i8
  use parallel_mod,   only: HME_BNDRY_A2A, HME_BNDRY_A2AO
  use thread_mod,     only: omp_in_parallel, omp_get_thread_num
  use gbarrier_mod,   only: gbarrier
  use cam_abortutils, only: endrun
  use cam_logfile,     only: iulog


  implicit none
  private

  interface bndry_exchange
     module procedure bndry_exchange_threaded
     module procedure bndry_exchange_nonthreaded
     module procedure long_bndry_exchange_nonth
  end interface
  public :: bndry_exchange

  interface ghost_exchange
     module procedure ghost_exchange_threaded
     module procedure ghost_exchange_nonthreaded
  end interface
  public :: ghost_exchange

  interface bndry_exchange_start
     module procedure bndry_exchange_threaded_start
     module procedure bndry_exchange_nonthreaded_start
  end interface
  public :: bndry_exchange_start

  interface bndry_exchange_finish
     module procedure bndry_exchange_threaded_finish
     module procedure bndry_exchange_nonthreaded_finish
  end interface
  public :: bndry_exchange_finish


  public :: compute_ghost_corner_orientation
  public :: ghost_exchangeVfull
  public :: copyBuffer

contains

  subroutine bndry_exchange_a2a(par,nthreads,ithr,buffer,location)
    use edgetype_mod,  only: Edgebuffer_t
    use schedtype_mod, only: schedule_t, cycle_t, schedule
    use thread_mod,    only: omp_in_parallel, omp_get_thread_num
    use perf_mod,      only: t_startf, t_stopf
    use spmd_utils,    only: mpi_real8, mpi_success
    use parallel_mod,  only: parallel_t
    use perf_mod,      only: t_startf, t_stopf

    type (parallel_t)            :: par
    integer, intent(in)          :: nthreads
    integer                      :: ithr ! The OpenMP thread ID
    type (EdgeBuffer_t)          :: buffer
    character(len=*),  optional  :: location

    type (Schedule_t), pointer   :: pSchedule
    type (Cycle_t),    pointer   :: pCycle
    integer                      :: icycle,ierr
    integer                      :: length
    integer                      :: iptr,source,nlyr
    integer                      :: nSendCycles,nRecvCycles
    integer                      :: errorcode,errorlen
    character*(80)               :: errorstring
    character(len=*),  parameter :: subname = 'bndry_exchange_a2a'
    character(len=80)            :: locstring
    logical                      :: ompthreadMissmatch

    integer                      :: i,j
    integer                      :: request

! Neighborhood collectives are only in MPI3 and up
#ifdef SPMD
#if MPI_VERSION >= 3

   if(ithr == 0) then

      call MPI_Ineighbor_Alltoallv(buffer%buf,buffer%scountsFull,buffer%sdisplsFull,Mpi_real8, &
                     buffer%receive,buffer%rcountsFull,buffer%rdisplsFull,Mpi_real8,par%commGraphFull,request,ierr)
      if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        write(iulog,*) subname,': Error after call to MPI_Ineighbor_alltoallv: ',errorstring
      endif

      if(present(location)) then
        locstring = TRIM(subname) // ': ' // TRIM(location)
      else
        locstring = TRIM(subname)
      endif
      ! location 1 for copyBuffer
      call t_startf('bndry_copy')
      call copyBuffer(nthreads,ithr,buffer,locstring)
      call t_stopf('bndry_copy')

      call MPI_wait(request,lstatus,ierr)
      call t_stopf('bndry_a2a')
   else

      if(present(location)) then
        locstring = TRIM(subname) // ': ' // TRIM(location)
      else
        locstring = TRIM(subname)
      endif
      call t_startf('bndry_copy')
      call copyBuffer(nthreads,ithr,buffer,locstring)
      call t_stopf('bndry_copy')

   endif
#else
    call endrun('bndry_exchange_a2a requires MPI-3 feature support')
#endif
#endif

  end subroutine bndry_exchange_a2a

  subroutine copyBuffer(nthreads,ithr,buffer,location)
    use edgetype_mod, only : Edgebuffer_t
    integer :: nthreads
    integer :: ithr
    type (EdgeBuffer_t)          :: buffer
    character(len=80)            :: location
    logical ::  ompThreadMissmatch
    integer lenMovePtr, iptr,length,i,j

    ompThreadMissmatch = .false.
    lenMovePtr = size(buffer%moveptr)
    if ( lenMOveptr .ne. nthreads) then
      ompthreadMissmatch = .true.
      write(*,30) TRIM(location), lenMoveptr, nthreads
    endif

    if (.not. ompthreadMissmatch) then
      iptr   = buffer%moveptr(ithr+1)
      length = buffer%moveLength(ithr+1)
      if(length>0) then
        do i=0,length-1
           buffer%receive(iptr+i) = buffer%buf(iptr+i)
        enddo
      endif
    else if(ompthreadMissmatch .and. ithr == 0) then
       do j=1,lenMovePtr
          iptr   = buffer%moveptr(j)
          length = buffer%moveLength(j)
          if(length>0) then
             do i=0,length-1
                buffer%receive(iptr+i) = buffer%buf(iptr+i)
             enddo
          endif
       enddo
    endif
30  format(a,'Potential performance issue: ',a,'LenMoveptr,nthreads: ',2(i3))
  end subroutine copyBuffer

  subroutine bndry_exchange_a2ao(par,nthreads,ithr,buffer,location)
    use edgetype_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
    use perf_mod, only : t_startf, t_stopf
    use spmd_utils,   only: mpi_real8, mpi_success, mpi_status_size
    use parallel_mod, only: parallel_t
    use perf_mod, only : t_startf, t_stopf

    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr  ! The OpenMP thread ID
    type (EdgeBuffer_t)               :: buffer
    character(len=*), optional        :: location

    integer                           :: ierr
    integer                           :: errorcode,errorlen
    character(len=80)                 :: errorstring
    character(len=*),       parameter :: subname = 'bndry_exchange_a2ao'
    character(len=80)                 :: locstring

    integer :: requestIntra,requestInter
    integer :: lstatus(MPI_status_size)

! Neighborhood collectives are only in MPI3 and up
#ifdef SPMD
#if MPI_VERSION >= 3

   if(ithr == 0) then

      call t_startf('bndry_a2ao')
      ! Start Inter-node communication
      call MPI_Ineighbor_Alltoallv(buffer%buf,buffer%scountsInter,buffer%sdisplsInter,MPI_real8, &
           buffer%receive,buffer%rcountsInter,buffer%rdisplsInter,MPI_real8,par%commGraphInter,requestInter,ierr)
      if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        write(iulog,*) subname,': Error after call to MPI_Ineighbor_alltoallv: ',errorstring
      endif
      ! Start Intra-node communication
      call MPI_Ineighbor_Alltoallv(buffer%buf,buffer%scountsIntra,buffer%sdisplsIntra,MPI_real8, &
           buffer%receive,buffer%rcountsIntra,buffer%rdisplsIntra,MPI_real8,par%commGraphIntra,requestIntra,ierr)
      if(ierr .ne. MPI_SUCCESS) then
        errorcode=ierr
        call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
        write(iulog,*) subname,': Error after call to MPI_Ineighbor_alltoallv: ',errorstring
      endif

      if(present(location)) then
        locstring = TRIM(subname) // ': ' // TRIM(location)
      else
        locstring = TRIM(subname)
      endif
      ! Finish the Intra-node communication
      call MPI_wait(requestIntra,lstatus,ierr)

      ! location 3 for copyBuffer
      call t_startf('bndry_copy')
      call copyBuffer(nthreads,ithr,buffer,locstring)
      call t_stopf('bndry_copy')

      ! Finish the Inter-node communication
      call MPI_wait(requestInter,lstatus,ierr)
      call t_stopf('bndry_a2ao')

   else

      if(present(location)) then
        locstring = TRIM(subname) // ': ' // TRIM(location)
      else
        locstring = TRIM(subname)
      endif
      !Copy buffer for ithr!=0
      call t_startf('bndry_copy')
      call copyBuffer(nthreads,ithr,buffer,locstring)
      call t_stopf('bndry_copy')

   endif
#else
    call endrun('bndry_exchange_a2ao requires MPI-3 feature support')
#endif
#endif

  end subroutine bndry_exchange_a2ao

  subroutine bndry_exchange_p2p(par,nthreads,ithr,buffer,location)
    use edgetype_mod,  only: Edgebuffer_t
    use schedtype_mod, only: schedule_t, cycle_t, schedule
    use thread_mod,    only: omp_in_parallel, omp_get_thread_num
    use spmd_utils,    only: mpi_real8, mpi_success
    use parallel_mod,  only: parallel_t
    use perf_mod,      only: t_startf, t_stopf

    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr
    type (EdgeBuffer_t)               :: buffer
    character(len=*), optional        :: location

    type (Schedule_t),pointer         :: pSchedule
    type (Cycle_t),pointer            :: pCycle
    integer                           :: dest,length,tag
    integer                           :: icycle,ierr
    integer                           :: iptr,source,nlyr
    integer                           :: nSendCycles,nRecvCycles
    integer                           :: errorcode,errorlen
    character*(80)                    :: errorstring
    character(len=*), parameter       :: subname = 'bndry_exchange_p2p'
    character(len=80)                 :: locstring
    logical, parameter :: Debug=.FALSE.

    integer                           :: i,j
    logical :: ompthreadMissmatch
    integer :: lenMovePtr

    pSchedule => Schedule(1)
    nlyr = buffer%nlyr
    ompthreadMissmatch = .FALSE.

    lenMovePtr = size(buffer%moveptr)

  if(ithr == 0) then
    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles


    !==================================================
    !  Fire off the sends
    !==================================================

    do icycle=1,nSendCycles
       pCycle      => pSchedule%SendCycle(icycle)
       dest            = pCycle%dest - 1
       length          = buffer%scountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%sdisplsFull(icycle) + 1
       if(Debug) write(iulog,*) subname,': MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
       call MPI_Isend(buffer%buf(iptr),length,Mpi_real8,dest,tag,par%comm,buffer%Srequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          write(iulog,*) subname,': Error after call to MPI_Isend: ',errorstring
       endif
    end do    ! icycle

    !==================================================
    !  Post the Receives
    !==================================================
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       source          = pCycle%source - 1
       length          = buffer%rcountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%rdisplsFull(icycle) + 1
       if(Debug) write(iulog,*) subname,': MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
       call MPI_Irecv(buffer%receive(iptr),length,Mpi_real8, &
            source,tag,par%comm,buffer%Rrequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          write(iulog,*) subname,': Error after call to MPI_Irecv: ',errorstring
       endif
    end do    ! icycle
    if(present(location)) then
      locstring = TRIM(subname) // ': ' // TRIM(location)
    else
      locstring = TRIM(subname)
    endif
    call t_startf('bndry_copy')
    call copyBuffer(nthreads,ithr,buffer,locstring)
    call t_stopf('bndry_copy')
    if (nSendCycles>0) call MPI_Waitall(nSendCycles,buffer%Srequest,buffer%status,ierr)
    if (nRecvCycles>0) call MPI_Waitall(nRecvCycles,buffer%Rrequest,buffer%status,ierr)
  else
    if(present(location)) then
      locstring = TRIM(subname) // ': ' // TRIM(location)
    else
      locstring = TRIM(subname)
    endif
    call t_startf('bndry_copy')
    call copyBuffer(nthreads,ithr,buffer,locstring)
    call t_stopf('bndry_copy')
  endif

  end subroutine bndry_exchange_p2p

  subroutine bndry_exchange_p2p_start(par,nthreads,ithr,buffer,location)

    use edgetype_mod,  only: Edgebuffer_t
    use schedtype_mod, only: schedule_t, cycle_t, schedule
    use thread_mod,    only: omp_in_parallel, omp_get_thread_num
    use spmd_utils,    only: mpi_real8, mpi_success
    use parallel_mod,  only: parallel_t

    type (parallel_t)                 :: par
    integer, intent(in)               :: nthreads
    integer                           :: ithr
    type (EdgeBuffer_t)               :: buffer
    character (len=*), optional :: location

    type (Schedule_t),pointer         :: pSchedule
    type (Cycle_t),pointer            :: pCycle
    integer                           :: dest,length,tag
    integer                           :: icycle,ierr
    integer                           :: iptr,source,nlyr
    integer                           :: nSendCycles,nRecvCycles
    integer                           :: errorcode,errorlen
    character*(80)                    :: errorstring
    character(len=*),       parameter :: subname = 'bndry_exchange_p2p_start'
    logical,                parameter :: Debug=.FALSE.

    integer                           :: i,j, lenMovePtr
    logical :: ompthreadMissmatch

    pSchedule => Schedule(1)
    nlyr = buffer%nlyr
    ompthreadMissmatch = .FALSE.

    lenMovePtr = size(buffer%moveptr)

  if(ithr == 0) then
    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles

    !==================================================
    !  Fire off the sends
    !==================================================

    do icycle=1,nSendCycles
       pCycle      => pSchedule%SendCycle(icycle)
       dest            = pCycle%dest - 1
       length          = buffer%scountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%sdisplsFull(icycle) + 1
       if(Debug) write(iulog,*) subname,': MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
       call MPI_Isend(buffer%buf(iptr),length,Mpi_real8,dest,tag,par%comm,buffer%Srequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          write(iulog,*) subname,': Error after call to MPI_Isend: ',errorstring
       endif
    end do    ! icycle

    !==================================================
    !  Post the Receives
    !==================================================
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       source          = pCycle%source - 1
       length          = buffer%rcountsFull(icycle)
       tag             = buffer%tag
       iptr            = buffer%rdisplsFull(icycle) + 1
       if(Debug) write(iulog,*) subname,': MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
       call MPI_Irecv(buffer%receive(iptr),length,Mpi_real8, &
            source,tag,par%comm,buffer%Rrequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          write(iulog,*) subname,': Error after call to MPI_Irecv: ',errorstring
       endif
    end do    ! icycle
  endif

  end subroutine bndry_exchange_p2p_start

  subroutine bndry_exchange_p2p_finish(par,nthreads,ithr,buffer,location)
    use edgetype_mod,  only: Edgebuffer_t
    use schedtype_mod, only: schedule_t, cycle_t, schedule
    use thread_mod,    only: omp_in_parallel, omp_get_thread_num
    use parallel_mod,  only: parallel_t
    use perf_mod,      only: t_startf, t_stopf


    type (parallel_t)            :: par
    integer, intent(in)          :: nthreads
    integer                      :: ithr
    type (EdgeBuffer_t)          :: buffer
    character(len=*),  optional  :: location

    type (Schedule_t), pointer   :: pSchedule
    type (Cycle_t),    pointer   :: pCycle
    integer                      :: dest,length,tag
    integer                      :: icycle,ierr
    integer                      :: iptr,source,nlyr
    integer                      :: nSendCycles,nRecvCycles
    integer                      :: errorcode,errorlen
    character*(80)               :: errorstring
    character(len=*),  parameter :: subname = 'bndry_exchange_p2p_finish'
    character(len=80)            :: locstring

    integer                      :: i,j
    logical                      :: ompthreadMissmatch
    integer                      :: lenMovePtr


  pSchedule => Schedule(1)
  if(present(location)) then
    locstring = TRIM(subname) // ': ' // TRIM(location)
  else
    locstring = TRIM(subname)
  endif
  call t_startf('bndry_copy')
  call copyBuffer(nthreads,ithr,buffer,locstring)
  call t_stopf('bndry_copy')

  if(ithr == 0) then

    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles

    if (nSendCycles>0) call MPI_Waitall(nSendCycles,buffer%Srequest,buffer%status,ierr)
    if (nRecvCycles>0) call MPI_Waitall(nRecvCycles,buffer%Rrequest,buffer%status,ierr)

  endif

  end subroutine bndry_exchange_p2p_finish

  subroutine long_bndry_exchange_nonth(par,buffer)
    use edgetype_mod,  only: LongEdgebuffer_t
    use schedtype_mod, only: schedule_t, cycle_t, schedule
    use thread_mod,    only: omp_in_parallel
    use parallel_mod,  only: parallel_t, status, srequest, rrequest
    use spmd_utils,    only: mpi_integer, mpi_success

    type (parallel_t)            :: par
    type (LongEdgeBuffer_t)      :: buffer

    type (Schedule_t), pointer   :: pSchedule
    type (Cycle_t),    pointer   :: pCycle
    integer                      :: dest,length,tag
    integer                      :: icycle,ierr
    integer                      :: iptr,source,nlyr
    integer                      :: nSendCycles,nRecvCycles
    integer                      :: errorcode,errorlen
    character*(80)               :: errorstring
    character(len=*),  parameter :: subname = 'long_bndry_exchange_nonth'

    integer                      :: i

#ifdef SPMD
    if(omp_in_parallel()) then
       print *,subname,': Warning you are calling a non-thread safe'
       print *,'         routine inside a threaded region....     '
       print *,'                Results are not predictable!!            '
    endif


    ! Setup the pointer to proper Schedule
    pSchedule => Schedule(1)
    nlyr = buffer%nlyr

    nSendCycles = pSchedule%nSendCycles
    nRecvCycles = pSchedule%nRecvCycles


    !==================================================
    !  Fire off the sends
    !==================================================

    do icycle=1,nSendCycles
       pCycle      => pSchedule%SendCycle(icycle)
       dest            = pCycle%dest - 1
       length      = nlyr * pCycle%lengthP
       tag             = pCycle%tag
       iptr            = pCycle%ptrP

       call MPI_Isend(buffer%buf(1,iptr),length,Mpi_integer,dest,tag,par%comm,Srequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          write(iulog,*) subname,': Error after call to MPI_Isend: ',errorstring
       endif
    end do    ! icycle

    !==================================================
    !  Post the Receives
    !==================================================
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       source          = pCycle%source - 1
       length      = nlyr * pCycle%lengthP
       tag             = pCycle%tag
       iptr            = pCycle%ptrP

       call MPI_Irecv(buffer%receive(1,iptr),length,Mpi_integer, &
            source,tag,par%comm,Rrequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          write(iulog,*) subname,': Error after call to MPI_Irecv: ',errorstring
       endif
    end do    ! icycle


    !==================================================
    !  Wait for all the receives to complete
    !==================================================

    if (nSendCycles>0) call MPI_Waitall(nSendCycles,Srequest,status,ierr)
    if (nRecvCycles>0) call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       length             = pCycle%lengthP
       iptr            = pCycle%ptrP
       do i=0,length-1
          buffer%buf(1:nlyr,iptr+i) = buffer%receive(1:nlyr,iptr+i)
       enddo
    end do   ! icycle

#endif

  end subroutine long_bndry_exchange_nonth
  !********************************************************************************
  !
  !********************************************************************************


 subroutine ghost_exchange_threaded(hybrid,buffer,location)
    use hybrid_mod,   only : hybrid_t
    use edgetype_mod, only : Edgebuffer_t

    implicit none

    type (hybrid_t)             :: hybrid
    type (EdgeBuffer_t)         :: buffer
    character(len=*), optional :: location 

    call bndry_exchange_threaded(hybrid,buffer,location)
 end subroutine ghost_exchange_threaded

 subroutine bndry_exchange_threaded(hybrid,buffer,location)
    use hybrid_mod, only : hybrid_t
    use edgetype_mod, only : Edgebuffer_t
    use perf_mod, only: t_startf, t_stopf, t_adj_detailf
    implicit none

    type (hybrid_t)             :: hybrid
    type (EdgeBuffer_t)         :: buffer
    character(len=*), optional  :: location

    character(len=*), parameter :: subname = 'bndry_exchange_threaded'

    call gbarrier(buffer%gbarrier, hybrid%ithr)
    if(buffer%bndry_type == HME_BNDRY_A2A) then
       call bndry_exchange_a2a(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_A2AO) then
       call bndry_exchange_a2ao(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    else
       call bndry_exchange_p2p(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    endif
    call gbarrier(buffer%gbarrier, hybrid%ithr)

 end subroutine bndry_exchange_threaded

 subroutine bndry_exchange_threaded_start(hybrid,buffer,location)
    use hybrid_mod, only : hybrid_t
    use edgetype_mod, only : Edgebuffer_t
    use perf_mod, only: t_startf, t_stopf, t_adj_detailf
    implicit none

    type (hybrid_t)             :: hybrid
    type (EdgeBuffer_t)         :: buffer
    character(len=*), optional  :: location

    character(len=*), parameter :: subname = 'bndry_exchange_threaded_start'

    call gbarrier(buffer%gbarrier, hybrid%ithr)
    call bndry_exchange_p2p_start(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)

 end subroutine bndry_exchange_threaded_start

 subroutine bndry_exchange_threaded_finish(hybrid,buffer,location)
    use hybrid_mod, only : hybrid_t
    use edgetype_mod, only : Edgebuffer_t
    use perf_mod, only: t_startf, t_stopf, t_adj_detailf
    implicit none

    type (hybrid_t)             :: hybrid
    type (EdgeBuffer_t)         :: buffer
    character(len=*), optional  :: location

    character(len=*), parameter :: subname = 'bndry_exchange_threaded_finish'

    call bndry_exchange_p2p_finish(hybrid%par,hybrid%nthreads,hybrid%ithr,buffer,location)
    call gbarrier(buffer%gbarrier, hybrid%ithr)

 end subroutine bndry_exchange_threaded_finish

 subroutine ghost_exchange_nonthreaded(par,buffer,location)
    use parallel_mod,   only : parallel_t
    use edgetype_mod, only : Edgebuffer_t
    type (parallel_t)          :: par
    type (EdgeBUffer_t)        :: buffer
    character(len=*), optional :: location 
    call bndry_exchange_nonthreaded(par,buffer,location)
 end subroutine ghost_exchange_nonthreaded

 subroutine bndry_exchange_nonthreaded(par,buffer,location)
    use parallel_mod, only : parallel_t
    use edgetype_mod, only : Edgebuffer_t
    use perf_mod, only: t_startf, t_stopf, t_adj_detailf
    implicit none

    type (parallel_t)           :: par
    type (EdgeBuffer_t)         :: buffer
    character(len=*), optional  :: location

    integer                     :: ithr
    integer                     :: nthreads
    character(len=*), parameter :: subname = 'bndry_exchange_nonthreaded'

    !$OMP BARRIER
    ithr=0
    nthreads = 1
    if(buffer%bndry_type == HME_BNDRY_A2A) then
       call bndry_exchange_a2a(par,nthreads,ithr,buffer,location)
    else if (buffer%bndry_type == HME_BNDRY_A2AO) then
       call bndry_exchange_a2ao(par,nthreads,ithr,buffer,location)
    else
       call bndry_exchange_p2p(par,nthreads,ithr,buffer,location)
    endif
    !$OMP BARRIER

  end subroutine bndry_exchange_nonthreaded

 subroutine bndry_exchange_nonthreaded_start(par,buffer,location)
    use parallel_mod, only : parallel_t
    use edgetype_mod, only : Edgebuffer_t
    use perf_mod, only: t_startf, t_stopf, t_adj_detailf
    implicit none

    type (parallel_t)           :: par
    type (EdgeBuffer_t)         :: buffer
    character (len=*), optional :: location

    integer                     :: ithr
    integer                     :: nthreads
    character(len=*), parameter :: subname = 'bndry_exchange_nonthreaded_start'

    !$OMP BARRIER
    ithr=0
    nthreads=1
    call bndry_exchange_p2p_start(par,nthreads,ithr,buffer,location)

 end subroutine bndry_exchange_nonthreaded_start

 subroutine bndry_exchange_nonthreaded_finish(par,buffer,location)
    use parallel_mod, only : parallel_t
    use edgetype_mod, only : Edgebuffer_t
    use perf_mod, only: t_startf, t_stopf, t_adj_detailf
    implicit none

    type (parallel_t)                 :: par
    integer                           :: ithr
    type (EdgeBuffer_t)               :: buffer
    character (len=*), optional :: location
    integer :: nthreads

    character(len=*), parameter :: subname = 'bndry_exchange_nonthreaded_finish'

    ithr=0
    nthreads=1
    call bndry_exchange_p2p_finish(par,nthreads,ithr,buffer,location)
    !$OMP BARRIER

  end subroutine bndry_exchange_nonthreaded_finish

  subroutine compute_ghost_corner_orientation(hybrid,elem,nets,nete)
!
!  this routine can NOT be called in a threaded region because then each thread
!  will have its on ghostbuffer.   initghostbufer3D() should detect this and abort.
!
  use dimensions_mod, only: nelemd, np
  use parallel_mod, only : syncmp
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use edgetype_mod, only : edgebuffer_t
  use edge_mod, only : ghostpack, ghostunpack, &
       initghostbuffer,freeghostbuffer

  use control_mod, only : north,south,east,west,neast, nwest, seast, swest

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  integer :: nets,nete
  type (edgeBuffer_t)    :: ghostbuf_cv

  real (kind=r8) :: cin(-1:4,-1:4,1,nets:nete)  !CE: fvm tracer
  real (kind=r8) :: cout(-1:4,-1:4,1,nets:nete)  !CE: fvm tracer
  integer :: i,j,ie,kptr,np1,np2,nc,nc1,nc2,k,nlev
  logical :: fail,fail1,fail2
  real (kind=r8) :: tol = 0.1_r8
  call syncmp(hybrid%par)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! first test on the Gauss Grid with same number of ghost cells:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nc=2   ! test using GLL interior points
  nc1=-1
  nc2=4

  nlev=1

  if (hybrid%nthreads > 1) then
     call endrun('ERROR: compute_ghost_corner_orientation must be called before threaded region')
  endif
  call initghostbuffer(hybrid%par,ghostbuf_cv,elem,nlev,nc,nc)


  cin = 0._r8
  do ie=nets,nete
     cin(1,1,1,ie)=  elem(ie)%gdofp(1,1)
     cin(nc,nc,1,ie)=  elem(ie)%gdofp(np,np)
     cin(1,nc,1,ie)=   elem(ie)%gdofp(1,np)
     cin(nc,1,1,ie)=  elem(ie)%gdofp(np,1)
  enddo
  cout=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  run ghost exchange on c array to get corner orientation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie=nets,nete
     kptr=0
     call ghostpack(ghostbuf_cv, cin(:,:,:,ie),nlev,kptr,ie)
  end do
  call ghost_exchange(hybrid,ghostbuf_cv)
  do ie=nets,nete
     kptr=0
     call ghostunpack(ghostbuf_cv, cout(:,:,:,ie),nlev,kptr,ie)
  enddo

!       nc +--------+
!        ^ | nw  ne |
!     j  | |        |
!        1 | sw  se |
!          +--------+
!           1 --> nc
!              i

! check SW corner
  do ie=nets,nete
     fail1=.false.
     fail2=.false.
     if ( elem(ie)%desc%putmapP_ghost(swest) /= -1) then
        if (abs(cout(nc1,1,1,ie)-cout(nc1,0,1,ie)) .gt. tol )  fail1=.true.
        if (abs(cout(1,nc1,1,ie)-cout(0,nc1,1,ie)).gt.tol) fail2=.true.
     endif
     if (fail1 .neqv. fail2 ) call endrun( 'ghost exchange SW orientation failure')
     if (fail1) then
        elem(ie)%desc%reverse(swest)=.true.
     endif
  enddo
! check SE corner
  do ie=nets,nete
     fail1=.false.
     fail2=.false.
     if ( elem(ie)%desc%putmapP_ghost(seast) /= -1) then
        if (abs(cout(nc2,1,1,ie)-cout(nc2,0,1,ie)) .gt. tol )  fail1=.true.
        if (abs(cout(nc+1,nc1,1,ie)-cout(nc,nc1,1,ie)).gt.tol) fail2=.true.
     endif
     if (fail1 .neqv. fail2 ) call endrun('ghost exchange SE orientation failure')
     if (fail1) then
        elem(ie)%desc%reverse(seast)=.true.
     endif
  enddo
! check NW corner
  do ie=nets,nete
     fail1=.false.
     fail2=.false.
     if ( elem(ie)%desc%putmapP_ghost(nwest) /= -1) then
        if (abs(cout(nc1,nc+1,1,ie)-cout(nc1,nc,1,ie)) .gt. tol )  fail1=.true.
        if (abs(cout(1,nc2,1,ie)-cout(0,nc2,1,ie)).gt.tol) fail2=.true.
     endif
     if (fail1 .neqv. fail2 ) call endrun( 'ghost exchange NW orientation failure')
     if (fail1) then
        elem(ie)%desc%reverse(nwest)=.true.
     endif
  enddo
! check NE corner
  do ie=nets,nete
     fail1=.false.
     fail2=.false.
     if ( elem(ie)%desc%putmapP_ghost(neast) /= -1) then
        if (abs(cout(nc2,nc+1,1,ie)-cout(nc2,nc,1,ie)) .gt. tol )  fail1=.true.
        if (abs(cout(nc+1,nc2,1,ie)-cout(nc,nc2,1,ie)).gt.tol) fail2=.true.
     endif
     if (fail1 .neqv. fail2 ) call endrun( 'ghost exchange NE orientation failure')
     if (fail1) then
        elem(ie)%desc%reverse(neast)=.true.
     endif
  enddo
  call freeghostbuffer(ghostbuf_cv)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  end ghost exchange corner orientation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end subroutine
  subroutine ghost_exchangeVfull(par,ithr,buffer)
!
!   MT 2011:  derived from bndry_exchange, but copies an entire
!             element of ghost cell information, including corner
!             elements.  Requres cubed-sphere grid
!
    use hybrid_mod, only : hybrid_t
    use edgetype_mod,   only: Ghostbuffer3D_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use dimensions_mod, only: nelemd
    use parallel_mod, only : status, srequest, rrequest, parallel_t
    use spmd_utils,  only: mpi_integer, mpi_success,mpi_real8

    implicit none
    type (parallel_t)                :: par
    integer                          :: ithr     ! hybrid%ithr 0 if called outside threaded region

    type (GhostBuffer3D_t)           :: buffer

    type (Schedule_t),pointer        :: pSchedule
    type (Cycle_t),pointer           :: pCycle
    integer                          :: dest,length,tag
    integer                          :: icycle,ierr
    integer                          :: iptr,source,nlyr
    integer                          :: nSendCycles,nRecvCycles
    integer                          :: errorcode,errorlen
    character(len=*), parameter      :: subname = 'ghost_exchangeVfull'
    character*(80) errorstring

    integer                          :: i,i1,i2

    !$OMP BARRIER
    if(ithr == 0) then


#ifdef SPMD
       ! Setup the pointer to proper Schedule
       pSchedule => Schedule(1)
       nlyr = buffer%nlyr

       nSendCycles = pSchedule%nSendCycles
       nRecvCycles = pSchedule%nRecvCycles

       !==================================================
       !  Fire off the sends
       !==================================================
       do icycle=1,nSendCycles
          pCycle      => pSchedule%SendCycle(icycle)
          dest            = pCycle%dest - 1
          length      = nlyr * pCycle%lengthP_ghost * buffer%elem_size
          tag             = pCycle%tag
          iptr            = pCycle%ptrP_ghost

          call MPI_Isend(buffer%buf(1,1,1,iptr),length,MPI_real8,dest,tag,par%comm,Srequest(icycle),ierr)
          if(ierr .ne. MPI_SUCCESS) then
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,subname,': Error after call to MPI_Isend: ',errorstring
          endif
       end do    ! icycle

       !==================================================
       !  Post the Receives
       !==================================================
       do icycle=1,nRecvCycles
          pCycle         => pSchedule%RecvCycle(icycle)
          source          = pCycle%source - 1
          length      = nlyr * pCycle%lengthP_ghost * buffer%elem_size
          tag             = pCycle%tag
          iptr            = pCycle%ptrP_ghost

          call MPI_Irecv(buffer%receive(1,1,1,iptr),length,MPI_real8, &
               source,tag,par%comm,Rrequest(icycle),ierr)
          if(ierr .ne. MPI_SUCCESS) then
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,subname,': Error after call to MPI_Irecv: ',errorstring
          endif
       end do    ! icycle


       !==================================================
       !  Wait for all the receives to complete
       !==================================================

       call MPI_Waitall(nSendCycles,Srequest,status,ierr)
       call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)

       do icycle=1,nRecvCycles
          pCycle         => pSchedule%RecvCycle(icycle)
          length             = pCycle%lengthP_ghost
          iptr            = pCycle%ptrP_ghost
          do i=0,length-1
             buffer%buf(:,:,1:nlyr,iptr+i) = buffer%receive(:,:,1:nlyr,iptr+i)
          enddo
       end do   ! icycle


#endif
    endif  ! if (hybrid%ithr == 0)
    !$OMP BARRIER

  end subroutine ghost_exchangeVfull


end module bndry_mod
