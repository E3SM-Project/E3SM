#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module bndry_mod
  use parallel_mod, only : abortmp
  implicit none
  private
  public :: bndry_exchangeV, ghost_exchangeVfull, compute_ghost_corner_orientation
  public :: ghost_exchangeV, ghost_exchangev3d

  interface bndry_exchangeV
     module procedure bndry_exchangeV_nonth
     module procedure long_bndry_exchangeV_nonth
     module procedure bndry_exchangeV_thsave 
  end interface

contains 

  subroutine bndry_exchangeV_nonth(par,buffer)
    use kinds, only : log_kind
    use edge_mod, only : Edgebuffer_t
    use schedule_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
    type (parallel_t)              :: par
    type (EdgeBuffer_t)            :: buffer

    type (Schedule_t),pointer                     :: pSchedule
    type (Cycle_t),pointer                        :: pCycle
    integer                                       :: dest,length,tag
    integer                                       :: icycle,ierr
    integer                                       :: iptr,source,nlyr
    integer                                       :: nSendCycles,nRecvCycles
    integer                                       :: errorcode,errorlen
    character*(80) errorstring

    logical(kind=log_kind),parameter              :: Debug=.FALSE.

    integer        :: i

#ifdef _MPI
    if(omp_in_parallel()) then 
       print *,'bndry_exchangeV: Warning you are calling a non-thread safe'
       print *,'		 routine inside a threaded region....     '
       print *,'                Results are not predictable!!            '
    endif


    ! Setup the pointer to proper Schedule
#ifdef _PREDICT
    pSchedule => Schedule(iam)
#else
    pSchedule => Schedule(1)
#endif
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
       !DBG if(Debug) print *,'bndry_exchangeV: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
       call MPI_Isend(buffer%buf(1,iptr),length,MPIreal_t,dest,tag,par%comm,Srequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,'bndry_exchangeV: Error after call to MPI_Isend: ',errorstring
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
       !DBG if(Debug) print *,'bndry_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
       call MPI_Irecv(buffer%receive(1,iptr),length,MPIreal_t, &
            source,tag,par%comm,Rrequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,'bndry_exchangeV: Error after call to MPI_Irecv: ',errorstring
       endif
    end do    ! icycle


    !==================================================
    !  Wait for all the receives to complete
    !==================================================

    call MPI_Waitall(nSendCycles,Srequest,status,ierr)
    call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       length             = pCycle%lengthP
       iptr            = pCycle%ptrP
       do i=0,length-1
          buffer%buf(1:nlyr,iptr+i) = buffer%receive(1:nlyr,iptr+i)
       enddo
    end do   ! icycle


#endif

  end subroutine bndry_exchangeV_nonth

  subroutine long_bndry_exchangeV_nonth(par,buffer)
    use kinds, only : log_kind
    use edge_mod, only : LongEdgebuffer_t
    use schedule_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel
#ifdef _MPI
    use parallel_mod, only : parallel_t, abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success
#else
    use parallel_mod, only : parallel_t, abortmp
#endif
    type (parallel_t)              :: par
    type (LongEdgeBuffer_t)            :: buffer

    type (Schedule_t),pointer                     :: pSchedule
    type (Cycle_t),pointer                        :: pCycle
    integer                                       :: dest,length,tag
    integer                                       :: icycle,ierr
    integer                                       :: iptr,source,nlyr
    integer                                       :: nSendCycles,nRecvCycles
    integer                                       :: errorcode,errorlen
    character*(80) errorstring

    logical(kind=log_kind),parameter              :: Debug=.FALSE.

    integer        :: i

#ifdef _MPI
    if(omp_in_parallel()) then 
       print *,'bndry_exchangeV: Warning you are calling a non-thread safe'
       print *,'		 routine inside a threaded region....     '
       print *,'                Results are not predictable!!            '
    endif


    ! Setup the pointer to proper Schedule
#ifdef _PREDICT
    pSchedule => Schedule(iam)
#else
    pSchedule => Schedule(1)
#endif
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
       !DBG if(Debug) print *,'bndry_exchangeV: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
       call MPI_Isend(buffer%buf(1,iptr),length,MPIinteger_t,dest,tag,par%comm,Srequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,'bndry_exchangeV: Error after call to MPI_Isend: ',errorstring
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
       !DBG if(Debug) print *,'bndry_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
       call MPI_Irecv(buffer%receive(1,iptr),length,MPIinteger_t, &
            source,tag,par%comm,Rrequest(icycle),ierr)
       if(ierr .ne. MPI_SUCCESS) then
          errorcode=ierr
          call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
          print *,'bndry_exchangeV: Error after call to MPI_Irecv: ',errorstring
       endif
    end do    ! icycle


    !==================================================
    !  Wait for all the receives to complete
    !==================================================

    call MPI_Waitall(nSendCycles,Srequest,status,ierr)
    call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)
    do icycle=1,nRecvCycles
       pCycle         => pSchedule%RecvCycle(icycle)
       length             = pCycle%lengthP
       iptr            = pCycle%ptrP
       do i=0,length-1
          buffer%buf(1:nlyr,iptr+i) = buffer%receive(1:nlyr,iptr+i)
       enddo
    end do   ! icycle


#endif

  end subroutine long_bndry_exchangeV_nonth
  !********************************************************************************
  !
  !********************************************************************************
  subroutine bndry_exchangeV_thsave(hybrid,buffer)
    use hybrid_mod, only : hybrid_t
    use kinds, only : log_kind
    use edge_mod, only : Edgebuffer_t
    use schedule_mod, only : schedule_t, cycle_t, schedule
    use dimensions_mod, only: nelemd, np
    use perf_mod, only: t_startf, t_stopf ! _EXTERNAL
#ifdef _MPI
    use parallel_mod, only : abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success
#else
    use parallel_mod, only : abortmp
#endif
    implicit none

    type (hybrid_t)                   :: hybrid
    type (EdgeBuffer_t)               :: buffer

    type (Schedule_t),pointer                     :: pSchedule
    type (Cycle_t),pointer                        :: pCycle
    integer                                       :: dest,length,tag
    integer                                       :: icycle,ierr
    integer                                       :: iptr,source,nlyr
    integer                                       :: nSendCycles,nRecvCycles
    integer                                       :: errorcode,errorlen
    character*(80) errorstring

    integer        :: i
    logical(kind=log_kind),parameter      :: Debug = .FALSE.


    call t_startf('bndry_exchange')
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif
    if(hybrid%ithr == 0) then 

#ifdef _MPI
       ! Setup the pointer to proper Schedule
#ifdef _PREDICT
       pSchedule => Schedule(iam)
#else
       pSchedule => Schedule(1)
#endif
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
          !DBG if(Debug) print *,'bndry_exchangeV: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
          call MPI_Isend(buffer%buf(1,iptr),length,MPIreal_t,dest,tag,hybrid%par%comm,Srequest(icycle),ierr)
          if(ierr .ne. MPI_SUCCESS) then
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,'bndry_exchangeV: Error after call to MPI_Isend: ',errorstring
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
          !DBG if(Debug) print *,'bndry_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
          call MPI_Irecv(buffer%receive(1,iptr),length,MPIreal_t, &
               source,tag,hybrid%par%comm,Rrequest(icycle),ierr)
          if(ierr .ne. MPI_SUCCESS) then
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,'bndry_exchangeV: Error after call to MPI_Irecv: ',errorstring
          endif
       end do    ! icycle


       !==================================================
       !  Wait for all the receives to complete
       !==================================================

       call MPI_Waitall(nSendCycles,Srequest,status,ierr)
       call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)

       do icycle=1,nRecvCycles
          pCycle         => pSchedule%RecvCycle(icycle)
          length             = pCycle%lengthP
          iptr            = pCycle%ptrP
          do i=0,length-1
             buffer%buf(1:nlyr,iptr+i) = buffer%receive(1:nlyr,iptr+i)
          enddo
       end do   ! icycle


#endif
    endif  ! if (hybrid%ithr == 0)
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif
    call t_stopf('bndry_exchange')

  end subroutine bndry_exchangeV_thsave





  subroutine ghost_exchangeVfull(hybrid,buffer,nc)
!
!   MT 2011:  derived from bndry_exchange, but copies an entire
!             element of ghost cell information, including corner
!             elements.  Requres cubed-sphere grid
!
    use hybrid_mod, only : hybrid_t
    use kinds, only : log_kind
    use edge_mod, only : Ghostbuffer_t
    use schedule_mod, only : schedule_t, cycle_t, schedule
    use dimensions_mod, only: nelemd
#ifdef _MPI
    use parallel_mod, only : abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success
#else
    use parallel_mod, only : abortmp
#endif
    implicit none

    type (hybrid_t)                   :: hybrid
    type (GhostBuffer_t)               :: buffer
    integer :: nc

    type (Schedule_t),pointer                     :: pSchedule
    type (Cycle_t),pointer                        :: pCycle
    integer                                       :: dest,length,tag
    integer                                       :: icycle,ierr
    integer                                       :: iptr,source,nlyr
    integer                                       :: nSendCycles,nRecvCycles
    integer                                       :: errorcode,errorlen
    character*(80) errorstring

    integer        :: i,i1,i2
    logical(kind=log_kind),parameter      :: Debug = .FALSE.


#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif
    if(hybrid%ithr == 0) then 

#ifdef _MPI
       ! Setup the pointer to proper Schedule
#ifdef _PREDICT
       pSchedule => Schedule(iam)
#else
       pSchedule => Schedule(1)
#endif
       nlyr = buffer%nlyr

       nSendCycles = pSchedule%nSendCycles
       nRecvCycles = pSchedule%nRecvCycles

       !==================================================
       !  Fire off the sends
       !==================================================
       do icycle=1,nSendCycles
          pCycle      => pSchedule%SendCycle(icycle)
          dest            = pCycle%dest - 1
          length      = nlyr * pCycle%lengthP_ghost*nc*nc
          tag             = pCycle%tag
          iptr            = pCycle%ptrP_ghost
          !print *,'ghost_exchangeV: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
          call MPI_Isend(buffer%buf(1,1,1,iptr),length,MPIreal_t,dest,tag,hybrid%par%comm,Srequest(icycle),ierr)
          if(ierr .ne. MPI_SUCCESS) then
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,'bndry_exchangeV: Error after call to MPI_Isend: ',errorstring
          endif
       end do    ! icycle

       !==================================================
       !  Post the Receives 
       !==================================================
       do icycle=1,nRecvCycles
          pCycle         => pSchedule%RecvCycle(icycle)
          source          = pCycle%source - 1
          length      = nlyr * pCycle%lengthP_ghost*nc*nc
          tag             = pCycle%tag
          iptr            = pCycle%ptrP_ghost
          !print *,'ghost_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
          call MPI_Irecv(buffer%receive(1,1,1,iptr),length,MPIreal_t, &
               source,tag,hybrid%par%comm,Rrequest(icycle),ierr)
          if(ierr .ne. MPI_SUCCESS) then
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,'bndry_exchangeV: Error after call to MPI_Irecv: ',errorstring
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
             do i2=1,nc
                do i1=1,nc
                   buffer%buf(i1,i2,1:nlyr,iptr+i) = buffer%receive(i1,i2,1:nlyr,iptr+i)
                enddo
             enddo
          enddo
       end do   ! icycle


#endif
    endif  ! if (hybrid%ithr == 0)
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif

  end subroutine ghost_exchangeVfull

  ! ===========================================
  !  GHOST_EXCHANGEV:
  !  Author: Christoph Erath
  !  derived from bndry_exchange, but copies an entire
  !             element of ghost cell information, including corner
  !             elements.  Requres cubed-sphere grid
  ! =========================================
 subroutine ghost_exchangeV(hybrid,buffer,nhc,npoints,ntrac)
!
!   2011:  derived from bndry_exchange, but copies an entire
!             element of ghost cell information, including corner
!             elements.  Requres cubed-sphere grid
!
    use hybrid_mod, only : hybrid_t
    use kinds, only : log_kind
    use edge_mod, only : Ghostbuffertr_t
    use schedule_mod, only : schedule_t, cycle_t, schedule
    use dimensions_mod, only: nelemd
#ifdef _MPI
    use parallel_mod, only : abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success
#else
    use parallel_mod, only : abortmp
#endif
    implicit none

    type (hybrid_t)                   :: hybrid
    type (GhostBuffertr_t)               :: buffer
    integer :: nhc,npoints,ntrac

    type (Schedule_t),pointer                     :: pSchedule
    type (Cycle_t),pointer                        :: pCycle
    integer                                       :: dest,length,tag
    integer                                       :: icycle,ierr
    integer                                       :: iptr,source,nlyr
    integer                                       :: nSendCycles,nRecvCycles
    integer                                       :: errorcode,errorlen
    character*(80) errorstring

    integer        :: i,i1,i2
    logical(kind=log_kind),parameter      :: Debug = .FALSE.


#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif
    if(hybrid%ithr == 0) then 

#ifdef _MPI
       ! Setup the pointer to proper Schedule
#ifdef _PREDICT
       pSchedule => Schedule(iam)
#else
       pSchedule => Schedule(1)
#endif
       nlyr = buffer%nlyr
              
       nSendCycles = pSchedule%nSendCycles
       nRecvCycles = pSchedule%nRecvCycles

       !==================================================
       !  Fire off the sends
       !==================================================
       do icycle=1,nSendCycles
          pCycle      => pSchedule%SendCycle(icycle)
          dest            = pCycle%dest - 1
          length      = nlyr * ntrac * pCycle%lengthP_ghost*nhc*npoints
          tag             = pCycle%tag
          iptr            = pCycle%ptrP_ghost
          !print *,'ghost_exchangeV: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
          call MPI_Isend(buffer%buf(1,1,1,1,iptr),length,MPIreal_t,dest,tag,hybrid%par%comm,Srequest(icycle),ierr)
          if(ierr .ne. MPI_SUCCESS) then
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,'bndry_exchangeV: Error after call to MPI_Isend: ',errorstring
          endif
       end do    ! icycle

       !==================================================
       !  Post the Receives 
       !==================================================
       do icycle=1,nRecvCycles
          pCycle         => pSchedule%RecvCycle(icycle)
          source          = pCycle%source - 1
          length      = nlyr * ntrac * pCycle%lengthP_ghost*nhc*npoints
          tag             = pCycle%tag
          iptr            = pCycle%ptrP_ghost
          !print *,'ghost_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
          call MPI_Irecv(buffer%receive(1,1,1,1,iptr),length,MPIreal_t, &
               source,tag,hybrid%par%comm,Rrequest(icycle),ierr)
          if(ierr .ne. MPI_SUCCESS) then
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,'bndry_exchangeV: Error after call to MPI_Irecv: ',errorstring
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
             do i2=1,nhc
                do i1=1,npoints
                   buffer%buf(i1,i2,1:nlyr,1:ntrac,iptr+i) = buffer%receive(i1,i2,1:nlyr,1:ntrac,iptr+i)
                enddo
             enddo
          enddo
       end do   ! icycle


#endif
    endif  ! if (hybrid%ithr == 0)
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif

  end subroutine ghost_exchangeV



  ! ===========================================
  !  GHOST_EXCHANGEV3d:
  !  Author: James overflet
  !  derived from bndry_exchange, but copies an entire
  !             element of ghost cell information, including corner
  !             elements.  Requres cubed-sphere grid
  ! =========================================
 subroutine ghost_exchangeV3d(hybrid, buffer)
!
!   2011:  derived from bndry_exchange, but copies an entire
!             element of ghost cell information, including corner
!             elements.  Requres cubed-sphere grid
!
    use hybrid_mod,    only : hybrid_t
    use kinds,         only : log_kind
    use edge_mod,      only : Ghostbuffer3d_t
    use schedule_mod,  only : schedule_t, cycle_t, schedule
    use dimensions_mod, only: nelemd
#ifdef _MPI
    use parallel_mod, only : abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success
#else
    use parallel_mod, only : abortmp
#endif
    implicit none

    type (hybrid_t)                   :: hybrid
    type (GhostBuffer3d_t)            :: buffer
    integer :: nhc,np

    type (Schedule_t),pointer                     :: pSchedule
    type (Cycle_t),pointer                        :: pCycle
    integer                                       :: dest,length,tag
    integer                                       :: icycle,ierr
    integer                                       :: iptr,source,nlyr
    integer                                       :: nSendCycles,nRecvCycles
    integer                                       :: errorcode,errorlen
    character*(80) errorstring

    integer        :: i,i1,i2
    logical(kind=log_kind),parameter      :: Debug = .FALSE.


#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif
    if(hybrid%ithr == 0) then 

#ifdef _MPI
       ! Setup the pointer to proper Schedule
#ifdef _PREDICT
       pSchedule => Schedule(iam)
#else
       pSchedule => Schedule(1)
#endif
       nlyr    = buffer%nlyr
       nhc     = buffer%nhc
       np      = buffer%np


       nSendCycles = pSchedule%nSendCycles
       nRecvCycles = pSchedule%nRecvCycles

       !==================================================
       !  Fire off the sends
       !==================================================
       do icycle=1,nSendCycles
          pCycle      => pSchedule%SendCycle(icycle)
          dest            = pCycle%dest - 1
          length          = nlyr * pCycle%lengthP_ghost*(nhc+1)*np
          tag             = pCycle%tag
          iptr            = pCycle%ptrP_ghost
          !print *,'ghost_exchangeV3d: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
          call MPI_Isend(buffer%buf(1,1,1,iptr),length,MPIreal_t,dest,tag,hybrid%par%comm,Srequest(icycle),ierr)
          if(ierr .ne. MPI_SUCCESS) then
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,'ghost_exchangeV3d: Error after call to MPI_Isend: ',errorstring
          endif
       end do    ! icycle

       !==================================================
       !  Post the Receives 
       !==================================================
       do icycle=1,nRecvCycles
          pCycle         => pSchedule%RecvCycle(icycle)
          source          = pCycle%source - 1
          length          = nlyr * pCycle%lengthP_ghost*(nhc+1)*np
          tag             = pCycle%tag
          iptr            = pCycle%ptrP_ghost
          !print *,'ghost_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
          call MPI_Irecv(buffer%receive(1,1,1,iptr),length,MPIreal_t, &
               source,tag,hybrid%par%comm,Rrequest(icycle),ierr)
          if(ierr .ne. MPI_SUCCESS) then
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,'ghost_exchangeV3d: Error after call to MPI_Irecv: ',errorstring
          endif
       end do    ! icycle


       !==================================================
       !  Wait for all the receives to complete
       !==================================================

       call MPI_Waitall(nSendCycles,Srequest,status,ierr)
       call MPI_Waitall(nRecvCycles,Rrequest,status,ierr)

       do icycle=1,nRecvCycles
          pCycle         => pSchedule%RecvCycle(icycle)
          length          = pCycle%lengthP_ghost
          iptr            = pCycle%ptrP_ghost
          do i=0,length-1
             buffer%buf(:, :, :, iptr+i) = buffer%receive(:, :, :, iptr+i)
!            do i2=1,nhc
!               do i1=1,np
!                  buffer%buf(i1, i2, 1:nlyr, iptr+i) = buffer%receive(i1, i2, 1:nlyr, iptr+i)
!               enddo
!            enddo
          enddo
       end do   ! icycle


#endif
    endif  ! if (hybrid%ithr == 0)
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif

  end subroutine ghost_exchangeV3d





  subroutine compute_ghost_corner_orientation(hybrid,elem,nets,nete)
  use kinds, only : real_kind
  use dimensions_mod, only: nelemd, np
  use parallel_mod, only : syncmp
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use edge_mod, only : ghostbuffer_t, ghostvpackfull, ghostvunpackfull, initghostbufferfull,&
       freeghostbuffer
  use control_mod, only : north,south,east,west,neast, nwest, seast, swest

  implicit none

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  integer :: nets,nete
  type (ghostBuffer_t)   :: ghostbuf,ghostbuf_cv

  real (kind=real_kind) :: cin(2,2,1,nets:nete)  !CE: fvm tracer
  real (kind=real_kind) :: cout(-1:4,-1:4,1,nets:nete)  !CE: fvm tracer
  integer :: i,j,ie,kptr,np1,np2,nc,nc1,nc2,k,nlev
  logical :: fail,fail1,fail2
  real (kind=real_kind) :: tol=.1
  call syncmp(hybrid%par)
!   if (hybrid%par%masterproc) print *,'computing ghost cell corner orientations'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! first test on the Gauss Grid with same number of ghost cells:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nc=2   ! test using GLL interior points
  nc1=-1
  nc2=4

  nlev=1
  call initghostbufferfull(ghostbuf_cv,nlev,nc)

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
     call ghostVpackfull(ghostbuf_cv, cin(:,:,:,ie),1,nc,nc,nlev,kptr,elem(ie)%desc)
  end do
  call ghost_exchangeVfull(hybrid,ghostbuf_cv,nc)
  do ie=nets,nete
     kptr=0
     call ghostVunpackfull(ghostbuf_cv, cout(:,:,:,ie), nc1,nc2,nc,nlev, kptr, elem(ie)%desc)
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
     if (fail1 .neqv. fail2 ) call abortmp( 'ghost exchange SW orientation failure')
     if (fail1) then
        elem(ie)%desc%reverse(swest)=.true.
        !print *,'reversion sw orientation ie',ie
        !print *,elem(ie)%desc%reverse(nwest),elem(ie)%desc%reverse(north),elem(ie)%desc%reverse(neast)
        !print *,elem(ie)%desc%reverse(west),' ',elem(ie)%desc%reverse(east)
        !print *,elem(ie)%desc%reverse(swest),elem(ie)%desc%reverse(south),elem(ie)%desc%reverse(seast)
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
     if (fail1 .neqv. fail2 ) call abortmp('ghost exchange SE orientation failure')
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
     if (fail1 .neqv. fail2 ) call abortmp( 'ghost exchange NW orientation failure')
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
     if (fail1 .neqv. fail2 ) call abortmp( 'ghost exchange NE orientation failure')
     if (fail1) then
        elem(ie)%desc%reverse(neast)=.true.
     endif
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  end ghost exchange corner orientation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end subroutine




end module bndry_mod
