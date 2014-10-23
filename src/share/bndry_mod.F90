#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module bndry_mod
  use parallel_mod, only : abortmp
  use edge_mod, only : Ghostbuffer3D_t
  implicit none
  private
  public :: bndry_exchangeV, ghost_exchangeVfull, compute_ghost_corner_orientation
  public :: ghost_exchangeV
!  public :: ghost_exchangev3d
  public :: sort_neighbor_buffer_mapping

  interface bndry_exchangeV
     module procedure bndry_exchangeV_nonth
     module procedure long_bndry_exchangeV_nonth
     module procedure bndry_exchangeV_thsave 
  end interface

contains 

  subroutine bndry_exchangeV_nonth(par,buffer)
    use kinds, only : log_kind
    use edge_mod, only : Edgebuffer_t
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use thread_mod, only : omp_in_parallel, omp_get_thread_num
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
    if(omp_get_thread_num() > 0) then
       print *,'bndry_exchangeV: Warning you are calling a non-thread safe'
       print *,'		 routine inside a threaded region....     '
       print *,'                Results are not predictable!!            '
    endif

#ifdef MPI_PERSISTENT
      pSchedule => Schedule(1)
      nlyr = buffer%nlyr
      nSendCycles = SIZE(buffer%Srequest)
      nRecvCycles = SIZE(buffer%Rrequest)
!      print *,'bndry_exchange: nSendCycles: ',nSendCycles, ' nRecvCycles:
!      ',nRecvCycles
      call MPI_startall(nRecvCycles,buffer%Rrequest,ierr)
      call MPI_startall(nSendCycles,buffer%Srequest,ierr)

      call MPI_Waitall(nSendCycles,buffer%Srequest,status,ierr)
      call MPI_Waitall(nRecvCycles,buffer%Rrequest,status,ierr)

#else

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

#endif

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
    use schedtype_mod, only : schedule_t, cycle_t, schedule
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
    use edge_mod, only : Edgebuffer_t
    use perf_mod, only: t_startf, t_stopf, t_adj_detailf
    implicit none

    type (hybrid_t)                   :: hybrid
    type (EdgeBuffer_t)               :: buffer

    call t_adj_detailf(+2)
    call t_startf('bndry_exchange')
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
    !$OMP MASTER
    ! Only the root thread calls the bndry_exchangeV_nonth
#endif
    call bndry_exchangeV_nonth(hybrid%par,buffer)
#if (defined HORIZ_OPENMP)
    !$OMP END MASTER
    !$OMP BARRIER
#endif
    call t_stopf('bndry_exchange')
    call t_adj_detailf(-2)

  end subroutine bndry_exchangeV_thsave



  subroutine ghost_exchangeVfull(par,ithr,buffer)
!
!   MT 2011:  derived from bndry_exchange, but copies an entire
!             element of ghost cell information, including corner
!             elements.  Requres cubed-sphere grid
!
    use hybrid_mod, only : hybrid_t
    use kinds, only : log_kind
    use schedtype_mod, only : schedule_t, cycle_t, schedule
    use dimensions_mod, only: nelemd
#ifdef _MPI
    use parallel_mod, only : abortmp, status, srequest, rrequest, &
         mpireal_t, mpiinteger_t, mpi_success, parallel_t
#else
    use parallel_mod, only : abortmp, parallel_t
#endif
    implicit none
    type (parallel_t)              :: par
    integer                        :: ithr     ! hybrid%ithr 0 if called outside threaded region

!    type (hybrid_t)                   :: hybrid
    type (GhostBuffer3D_t)               :: buffer

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
    if(ithr == 0) then 


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
          length      = nlyr * pCycle%lengthP_ghost * buffer%elem_size
          tag             = pCycle%tag
          iptr            = pCycle%ptrP_ghost
          !print *,'ghost_exchangeV: MPI_Isend: DEST:',dest,'LENGTH:',length,'TAG: ',tag
          call MPI_Isend(buffer%buf(1,1,1,iptr),length,MPIreal_t,dest,tag,par%comm,Srequest(icycle),ierr)
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
          length      = nlyr * pCycle%lengthP_ghost * buffer%elem_size
          tag             = pCycle%tag
          iptr            = pCycle%ptrP_ghost
          !print *,'ghost_exchangeV: MPI_Irecv: SRC:',source,'LENGTH:',length,'TAG: ',tag
          call MPI_Irecv(buffer%receive(1,1,1,iptr),length,MPIreal_t, &
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
          length             = pCycle%lengthP_ghost
          iptr            = pCycle%ptrP_ghost
          do i=0,length-1
             buffer%buf(:,:,1:nlyr,iptr+i) = buffer%receive(:,:,1:nlyr,iptr+i)
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
    use schedtype_mod, only : schedule_t, cycle_t, schedule
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

  subroutine compute_ghost_corner_orientation(hybrid,elem,nets,nete)
!
!  this routine can NOT be called in a threaded region because then each thread
!  will have its on ghostbuffer.   initghostbufer3D() should detect this and abort.
!
  use kinds, only : real_kind
  use dimensions_mod, only: nelemd, np
  use parallel_mod, only : syncmp
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use edge_mod, only : ghostbuffer3D_t, ghostvpackfull, ghostvunpackfull, &
       initghostbuffer3D,freeghostbuffer3D
  use control_mod, only : north,south,east,west,neast, nwest, seast, swest

  implicit none

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  integer :: nets,nete
  type (ghostBuffer3D_t)   :: ghostbuf_cv

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
  call initghostbuffer3D(ghostbuf_cv,nlev,nc)


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
  call ghost_exchangeVfull(hybrid%par,hybrid%ithr,ghostbuf_cv)
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
  call freeghostbuffer3D(ghostbuf_cv)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  end ghost exchange corner orientation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end subroutine





  subroutine sort_neighbor_buffer_mapping(par,elem,nets,nete)
!
!  gather global ID's of all neighbor elements.  Then create sorted (in global ID numbering)
!  mapping between edge buffer for each neighbor and a local map.  
!
!  this routine can NOT be called in a threaded region because then each thread
!  will have its on ghostbuffer.   initghostbufer3D() should detect this and abort.
!
!  also return num_neigh(ie) = number of neighbors (including onself) for element ie
!  
!
  use kinds, only : real_kind
  use dimensions_mod, only: nelemd, np, max_neigh_edges
  use parallel_mod, only : syncmp, parallel_t
  use element_mod, only : element_t
  use edge_mod, only : ghostbuffer3D_t, ghostvpack_unoriented, ghostvunpack_unoriented, &
       initghostbuffer3D,freeghostbuffer3D
  use control_mod, only : north,south,east,west,neast, nwest, seast, swest
  use coordinate_systems_mod, only: cartesian3D_t
  implicit none

  type (parallel_t)      , intent(in) :: par
  type (element_t)     , intent(inout), target :: elem(:)
  integer :: nets,nete
  type (ghostBuffer3D_t)   :: ghostbuf_cv

  real (kind=real_kind) :: cin(2,2,4,nets:nete)                    ! 1x1 element input data
  real (kind=real_kind) :: cout(2,2,4,max_neigh_edges+1,nets:nete)   ! 1x1 element output data
  real (kind=real_kind) :: u   (2,2,4)   
  integer :: i,j,ie,kptr,np1,np2,nc,k,nlev,patch_size,l,l2,sum1,sum2,m
  logical :: fail,fail1,fail2
  real (kind=real_kind) :: tol=.1


  if (par%masterproc) print *,'creating sorted ghost cell neigbor map...' 
  if (par%masterproc) print *,'checking ghost cell neighbor buffer sorting...' 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! first test on the Gauss Grid with same number of ghost cells:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nc=2
  nlev=4
  call initghostbuffer3D(ghostbuf_cv,nlev,nc)

  call syncmp(par)

  do ie=nets,nete
     cin(:,:,nlev,ie)=  elem(ie)%GlobalID
     k=0
     do i=1,nc
     do j=1,nc
        k=k+1
        cin(i,j,1,ie) = elem(ie)%corners3D(k)%x
        cin(i,j,2,ie) = elem(ie)%corners3D(k)%y
        cin(i,j,3,ie) = elem(ie)%corners3D(k)%z
     enddo
     enddo
  enddo
  cout=-1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  run ghost exchange to get global ID of all neighbors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie=nets,nete
     kptr=0
     call ghostVpack_unoriented(ghostbuf_cv, cin(:,:,:,ie),nc,nlev,kptr,elem(ie)%desc)
  end do

  ! check for array out of bouds overwriting 
  if (int(maxval(  cout(:,:,:,:,:))) /= -1 ) then
     call abortmp('ghost excchange unoriented failure ob1')
  endif
  call ghost_exchangeVfull(par,0,ghostbuf_cv)
  if (int(maxval(  cout(:,:,:,:,:))) /= -1 ) then
     call abortmp('ghost excchange unoriented failure ob2')
  endif

  do ie=nets,nete
     kptr=0

     call ghostVunpack_unoriented(ghostbuf_cv, cout(:,:,:,:,ie),nc,nlev, kptr, elem(ie)%desc, elem(ie)%GlobalId,cin(:,:,:,ie))

     ! check that we get the count of real neighbors correct
     patch_size=0

     do l=1,max_neigh_edges+1
        if (int(cout(1,1,nlev,l,ie)) /= -1 ) then
           patch_size = patch_size + 1
        endif
     enddo

     if (elem(ie)%desc%actual_neigh_edges+1 /= patch_size) then
        print *,'desc  actual_neigh_edges: ',elem(ie)%desc%actual_neigh_edges
        print *,'check patch_size: ',patch_size
        call abortmp( 'ghost exchange unoriented failure 1')
     endif

     ! check that all non-neighbors stayed -1
     do l=patch_size+1,max_neigh_edges+1
     if (int(cout(1,1,nlev,l,ie)) /= -1 ) then
        call abortmp( 'ghost exchange unoriented failure 2')
     endif
     enddo

     ! i am too lazy to check if all id's are identical since they are in
     ! different order.  check if there sum is identical
     sum1 = sum(int(cout(1,1,nlev,1:patch_size,ie)))
     sum2 = elem(ie)%globalID
     do l=1,max_neigh_edges
        if (elem(ie)%desc%globalID(l)>0) sum2 = sum2 + elem(ie)%desc%globalID(l)
     enddo
     if (sum1 /= sum2 ) then
        print *,int(cin(1,1,nlev,ie)),elem(ie)%desc%actual_neigh_edges,patch_size
        write(*,'(a,99i5)') 'ghost=',int(cout(1,1,nlev,1:patch_size,ie))
        write(*,'(a,99i5)') 'desc =',elem(ie)%desc%globalID(:)

        print *,'cout sum of all neighbor global ids:',sum1 
        print *,'desc sum of all neighbor global ids:',sum2  
        call abortmp( 'ghost exchange unoriented failure 3')        
     endif

     ALLOCATE(elem(ie)%desc%neigh_corners(4,patch_size))
     ! unpack corner data into array
     do l=1,patch_size
        k=0
        do i=1,nc
        do j=1,nc
           k=k+1
           elem(ie)%desc%neigh_corners(k,l)%x = cout(i,j,1,l,ie) 
           elem(ie)%desc%neigh_corners(k,l)%y = cout(i,j,2,l,ie) 
           elem(ie)%desc%neigh_corners(k,l)%z = cout(i,j,3,l,ie) 
        enddo
        enddo
     enddo
  enddo

  call freeghostbuffer3D(ghostbuf_cv)
  if (par%masterproc) print *,'passed.'
  end subroutine




end module bndry_mod
