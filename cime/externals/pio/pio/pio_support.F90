#define __PIO_FILE__ "pio_support.F90"
!>
!! @file pio_support.F90
!! @brief internal code for compiler workarounds, aborts and debug functions
!!
!! $Revision$
!! $LastChangedDate$
!<
!>
!! \def _NO_MPI_RSEND
!! Code added as a work around for poor rsend performance on cray systems with
!! Gemini interconnect
!<
#ifdef BGP
#define BGx
#endif
#ifdef BGL
#define BGx
#endif
#ifdef BGQ
#define BGx
#endif
#ifdef _NO_MPI_RSEND
#define MPI_RSEND MPI_SEND
#define mpi_rsend mpi_send
#define MPI_IRSEND MPI_ISEND
#define mpi_irsend mpi_isend
#endif

module pio_support
  use pio_kinds
#ifndef NO_MPIMOD
  use mpi !_EXTERNAL
#endif
  implicit none
  private
#ifdef NO_MPIMOD
  include 'mpif.h'    ! _EXTERNAL
#endif
  public :: piodie
  public :: CheckMPIreturn
  public :: pio_readdof
  public :: pio_writedof
  public :: pio_fc_gather_offset
#ifdef NO_MPI2
  public :: MPI_TYPE_CREATE_INDEXED_BLOCK
#endif


  logical, public :: Debug=.FALSE.
  logical, public :: DebugIO=.FALSE.
  logical, public :: DebugAsync=.FALSE.
  integer,private,parameter :: versno = 1001

  character(len=*), parameter :: modName='pio_support'

contains

  subroutine piodie (file,line, msg, ival1, msg2, ival2, msg3, ival3, mpirank)
#ifdef CPRINTEL
  ! tracebackqq uses optional arguments, so *must* have an explicit
  ! interface.
  use ifcore, only: tracebackqq
#endif
    !-----------------------------------------------------------------------
    ! Purpose:
    !
    ! Abort the model for abnormal termination
    !
    ! Author: Jim Edwards
    !
    ! Change History
    ! 20070608 R. Loy  added optional args
    !-----------------------------------------------------------------------
    ! $Id$
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    character(len=*), intent(in) :: file
    integer,intent(in) :: line
    character(len=*), intent(in), optional :: msg,msg2,msg3
    integer,intent(in),optional :: ival1,ival2,ival3, mpirank

    character(len=*), parameter :: subName=modName//'::pio_die'
    integer :: ierr, myrank=-1
    
    if(present(mpirank)) myrank=mpirank

    if (present(ival3)) then
       write(6,*) subName,':: myrank=',myrank,': ERROR: ',file,':',line,': ', &
            msg,ival1,msg2,ival2,msg3,ival3
    else if (present(msg3)) then
       write(6,*) subName,':: myrank=',myrank,': ERROR: ',file,':',line,': ', &
            msg,ival1,msg2,ival2,msg3
    else if (present(ival2)) then
       write(6,*) subName,':: myrank=',myrank,': ERROR: ',file,':',line,': ',msg,ival1,msg2,ival2
    else if (present(msg2)) then
       write(6,*) subName,':: myrank=',myrank,': ERROR: ',file,':',line,': ',msg,ival1,msg2
    else if (present(ival1)) then
       write(6,*) subName,':: myrank=',myrank,': ERROR: ',file,':',line,': ',msg,ival1
    else if (present(msg)) then
       write(6,*) subName,':: myrank=',myrank,': ERROR: ',file,':',line,': ',msg
    else
       write(6,*) subName,':: myrank=',myrank,': ERROR: ',file,':',line,': (no message)'
    endif


#if defined(CPRXLF) && !defined(BGx)
  close(5)    ! needed to prevent batch jobs from hanging in xl__trbk
  call xl__trbk()
#elif defined(CPRINTEL)


  ! An exit code of -1 is a special value that prevents this subroutine
  ! from aborting the run.
  call tracebackqq(user_exit_code=-1)

#endif

    ! passing an argument of 1 to mpi_abort will lead to a STOPALL output 
    ! error code of 257
    call mpi_abort (MPI_COMM_WORLD, 1, ierr)  

#ifdef CPRNAG
    stop
#else
    call abort
#endif


  end subroutine piodie

!=============================================
!  CheckMPIreturn:
!
!      Check and prints an error message
!  if an error occured in a MPI subroutine.
!=============================================
  subroutine CheckMPIreturn(locmesg, errcode, file, line)

     character(len=*), intent(in) :: locmesg
     integer(i4), intent(in) :: errcode
     character(len=*),optional :: file
     integer, intent(in),optional :: line
     character(len=MPI_MAX_ERROR_STRING) :: errorstring

     integer(i4) :: errorlen

     integer(i4) :: ierr
     if (errcode .ne. MPI_SUCCESS) then
        call MPI_Error_String(errcode,errorstring,errorlen,ierr)
        write(*,*) TRIM(ADJUSTL(locmesg))//errorstring(1:errorlen)
        if(present(file).and.present(line)) then
           call piodie(file,line)
        endif
     end if
  end subroutine CheckMPIreturn

  subroutine pio_writedof (file, DOF, comm, punit)
    !-----------------------------------------------------------------------
    ! Purpose:
    !
    ! Write a DOF to standard format
    !
    ! Author: T Craig
    !
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    character(len=*),intent(in) :: file
    integer(kind=pio_offset)  ,intent(in) :: dof(:)
    integer         ,intent(in) :: comm
    integer,optional,intent(in) :: punit

    character(len=*), parameter :: subName=modName//'::pio_writedof'
    integer ierr, myrank, npes, m, n, unit
    integer(kind=pio_offset), pointer :: wdof(:)
    integer(kind=pio_offset), pointer :: sdof1d(:)
    integer(kind=pio_offset) :: sdof, sdof_tmp(1)
    integer          :: status(MPI_STATUS_SIZE)
    integer, parameter :: masterproc = 0

    integer :: pio_offset_kind                     ! kind of pio_offset
#ifndef _NO_FLOW_CONTROL
    integer :: &
      rcv_request    ,&! request id
      hs = 1           ! MPI handshaking variable
#endif

    unit = 81
    if (present(punit)) then
       unit = punit
    endif

    call MPI_COMM_SIZE(comm,npes,ierr)
    call CheckMPIReturn(subName,ierr)
    call MPI_COMM_RANK(comm,myrank,ierr)
    call CheckMPIReturn(subName,ierr)
    sdof = size(dof)

    allocate(sdof1d(0:npes-1))
    sdof1d = -1
    sdof_tmp(1) = sdof

    if(kind(sdof_tmp) == kind(comm)) then
       pio_offset_kind = MPI_INTEGER
    else
       pio_offset_kind = MPI_INTEGER8
    end if



    call pio_fc_gather_offset(sdof_tmp, 1, PIO_OFFSET_KIND, &
       sdof1d, 1, PIO_OFFSET_KIND,masterproc,comm)

    if (myrank == masterproc) then
       write(6,*) subName,': writing file ',trim(file),' unit=',unit
       open(unit,file=file)
       write(unit,*) versno,npes
    endif

    do n = 0,npes-1
       if (myrank == masterproc) then
          allocate(wdof(sdof1d(n)))
       endif
       if (myrank == masterproc .and. n == masterproc) then
          wdof = dof
       else
          if (myrank == n .and. sdof > 0) then
#ifndef _NO_FLOW_CONTROL
             call MPI_RECV(hs,1,MPI_INTEGER,masterproc,n,comm,status,ierr)
             if (ierr /= MPI_SUCCESS) call piodie(__PIO_FILE__,__LINE__,' pio_writedof mpi_recv')
#endif
             call MPI_SEND(dof,int(sdof),PIO_OFFSET_KIND,masterproc,n,comm,ierr)
             if (ierr /= MPI_SUCCESS) call piodie(__PIO_FILE__,__LINE__,' pio_writedof mpi_send')
          endif
          if (myrank == masterproc .and. sdof1d(n) > 0) then
#ifndef _NO_FLOW_CONTROL
             call MPI_IRECV(wdof,int(sdof1d(n)),PIO_OFFSET_KIND,n,n,comm,rcv_request,ierr)
             if (ierr /= MPI_SUCCESS) call piodie(__PIO_FILE__,__LINE__,' pio_writedof mpi_irecv')
             call MPI_SEND(hs,1,MPI_INTEGER,n,n,comm,ierr)
             if (ierr /= MPI_SUCCESS) call piodie(__PIO_FILE__,__LINE__,' pio_writedof mpi_send')
             call MPI_WAIT(rcv_request,status,ierr)
             if (ierr /= MPI_SUCCESS) call piodie(__PIO_FILE__,__LINE__,' pio_writedof mpi_wait')
#else
             call MPI_RECV(wdof,int(sdof1d(n)),PIO_OFFSET_KIND,n,n,comm,status,ierr)
             if (ierr /= MPI_SUCCESS) call piodie(__PIO_FILE__,__LINE__,' pio_writedof mpi_recv')
#endif
          endif
       endif
       if (myrank == masterproc) then
          write(unit,*) n,sdof1d(n)
          do m = 1,sdof1d(n)
             write(unit,*) wdof(m)
          enddo
          deallocate(wdof)
       endif
    enddo

    if (myrank == masterproc) then
       close(unit)
    endif

    deallocate(sdof1d)

  end subroutine pio_writedof

  subroutine pio_readdof (file, DOF, comm, punit)
    !-----------------------------------------------------------------------
    ! Purpose:
    !
    ! Read a DOF to standard format
    !
    ! Author: T Craig
    !
    ! Change History
    ! 
    !-----------------------------------------------------------------------
    ! $Id$
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    character(len=*),intent(in) :: file
    integer(kind=pio_offset),pointer:: dof(:)
    integer         ,intent(in) :: comm
    integer,optional,intent(in) :: punit

    character(len=*), parameter :: subName=modName//'::pio_readdof'
    integer :: ierr, myrank, npes, m, n, unit, rn
    integer(kind=pio_offset) :: sdof
    integer :: rversno, rnpes
    integer(kind=pio_offset), pointer :: wdof(:)
    integer, parameter :: masterproc = 0
    integer :: status(MPI_STATUS_SIZE)
    integer :: pio_offset_kind                        ! kind of pio_offset

    unit = 81
    if (present(punit)) then
       unit = punit
    endif

    call MPI_COMM_SIZE(comm,npes,ierr)
    call CheckMPIReturn(subName,ierr)
    call MPI_COMM_RANK(comm,myrank,ierr)
    call CheckMPIReturn(subName,ierr)

    if(kind(sdof) == kind(comm)) then
       pio_offset_kind = MPI_INTEGER
    else
       pio_offset_kind = MPI_INTEGER8
    end if


    allocate(dof(0))   ! default for pes with no dof

    if (myrank == masterproc) then
       write(6,*) subName,': reading file ',trim(file),' unit=',unit
       open(unit,file=file,status='old')
       read(unit,*) rversno,rnpes
       write(6,*) subName,': reading file ',trim(file),' versno=',rversno
       if (rnpes /= npes) then
          call piodie(__PIO_FILE__,__LINE__,'pio_readdof npes incorrect')
       endif

       do n = 0,npes-1
          read(unit,*) rn,sdof
          if (rn /= n) then
             call piodie(__PIO_FILE__,__LINE__,'pio_readdof rn out of sync')
          endif
          allocate(wdof(sdof))
          do m = 1,sdof
             read(unit,*) wdof(m)
          enddo
          if (n == masterproc) then
             deallocate(dof)
             allocate(dof(sdof))
             dof = wdof
          else
             call MPI_SEND(sdof,1,PIO_OFFSET_KIND,n,n,comm,ierr)
             if (ierr /= MPI_SUCCESS) call piodie(__PIO_FILE__,__LINE__,' pio_readdof mpi_send1')
             if (sdof > 0) then
                call MPI_SEND(wdof,int(sdof),PIO_OFFSET_KIND,n,npes+n,comm,ierr)
                if (ierr /= MPI_SUCCESS) call piodie(__PIO_FILE__,__LINE__,' pio_readdof mpi_send2')
             endif
          endif
          deallocate(wdof)
       enddo
       close(unit)
    else
       call MPI_RECV(sdof,1,PIO_OFFSET_KIND,masterproc,myrank,comm,status,ierr)
       if (ierr /= MPI_SUCCESS) call piodie(__PIO_FILE__,__LINE__,' pio_readdof mpi_recv1')
       if (sdof > 0) then
          deallocate(dof)
          allocate(dof(sdof))
          call MPI_RECV(dof,int(sdof),PIO_OFFSET_KIND,masterproc,npes+myrank,comm,status,ierr)
          if (ierr /= MPI_SUCCESS) call piodie(__PIO_FILE__,__LINE__,' pio_readdof mpi_recv2')
       endif
    endif

  end subroutine pio_readdof

#ifdef NO_MPI2

  subroutine MPI_TYPE_CREATE_INDEXED_BLOCK(count, blen, disp, oldtype, newtype, ierr)
    integer, intent(in)  :: count
    integer, intent(in)  :: blen
    integer, intent(in)  :: disp(:)
    integer, intent(in)  :: oldtype
    integer, intent(out) :: newtype
    integer :: ierr

    integer :: dblens(count)

    dblens = blen
#ifndef _MPISERIAL
    call mpi_type_indexed(count, dblens, disp, oldtype, newtype, ierr)
#endif
  end subroutine MPI_TYPE_CREATE_INDEXED_BLOCK
#endif

!
!========================================================================
!

   subroutine pio_fc_gather_offset ( sendbuf, sendcnt, sendtype, &
                                  recvbuf, recvcnt, recvtype, &
                                  root, comm, flow_cntl )

!----------------------------------------------------------------------- 
! 
!> Purpose: 
!!   Gather collective with additional flow control, so as to 
!!   be more robust when used with high process counts. 
!!
!! Method: 
!!   If flow_cntl optional parameter 
!!     < 0: use MPI_Gather
!!     >= 0: use point-to-point with handshaking messages and 
!!           preposting receive requests up to 
!!           max(min(1,flow_cntl),max_gather_block_size) 
!!           ahead if optional flow_cntl parameter is present.
!!           Otherwise, fc_gather_flow_cntl is used in its place.
!!     Default value is 64.
!! 
!! Author of original version:  P. Worley
!! Ported from CAM: P. Worley, Jan 2010
!< 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
   implicit none

!---------------------------Parameters ---------------------------------
!
   integer, parameter :: max_gather_block_size = 64

!---------------------------Input arguments--------------------------
!
   integer(kind=pio_offset), intent(in)  :: sendbuf(:)       ! outgoing message buffer
   integer, intent(in)  :: sendcnt          ! size of send buffer
   integer, intent(in)  :: sendtype         ! MPI type of send buffer
   integer, intent(in)  :: recvcnt          ! size of receive buffer
   integer, intent(in)  :: recvtype         ! MPI type of receive buffer
   integer, intent(in)  :: root             ! gather destination
   integer, intent(in)  :: comm             ! MPI communicator
   integer,optional, intent(in):: flow_cntl ! flow control variable

!---------------------------Output arguments--------------------------
!
   integer(kind=pio_offset), intent(out) :: recvbuf(*)       ! incoming message buffer
!
!---------------------------Local workspace---------------------------------
!
   character(len=*), parameter :: subName=modName//'::pio_fc_gather_int'

   logical :: fc_gather                     ! use explicit flow control?
   integer :: hs                            ! handshake variable
   integer :: gather_block_size             ! number of preposted receive requests

   integer :: nprocs                        ! size of communicator
   integer :: mytask                        ! MPI task id with communicator
   integer :: mtag                          ! MPI message tag
   integer :: p, i                          ! loop indices
   integer :: displs                        ! offset into receive buffer
   integer :: count, preposts, head, tail   ! variables controlling recv-ahead logic

   integer :: rcvid(max_gather_block_size)  ! receive request ids

   integer :: ier                           ! return error status    
   integer :: status(MPI_STATUS_SIZE)       ! MPI status 

!
!-------------------------------------------------------------------------------------
!
   if ( present(flow_cntl) ) then
      if (flow_cntl >= 0) then
         gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
         fc_gather = .true.
      else
         fc_gather = .false.
      endif
   else
#ifndef _NO_FLOW_CONTROL
      gather_block_size = max(1,max_gather_block_size)
      fc_gather = .true.
#else
      fc_gather = .false.
#endif
   endif

   if (fc_gather) then
 
      ! Determine task id and size of communicator
      call mpi_comm_rank (comm, mytask, ier)
      call mpi_comm_size (comm, nprocs, ier)

      ! Initialize tag and hs variable
#ifdef _NO_PIO_SWAPM_TAG_OFFSET
      mtag = 0
#else
      mtag = 2*nprocs
#endif
      hs = 1

      if (root .eq. mytask) then

! prepost gather_block_size irecvs, and start receiving data
         preposts = min(nprocs-1, gather_block_size)
         head = 0
         count = 0
         do p=0, nprocs-1
            if (p .ne. root) then
               if (recvcnt > 0) then
                  count = count + 1
                  if (count > preposts) then
                     tail = mod(head,preposts) + 1
                     call mpi_wait (rcvid(tail), status, ier)
                  end if
                  head = mod(head,preposts) + 1
                  displs = p*recvcnt
                  call mpi_irecv ( recvbuf(displs+1), recvcnt, &
                                   recvtype, p, mtag, comm, rcvid(head), &
                                   ier )
                  call mpi_send ( hs, 1, MPI_INTEGER, p, mtag, comm, ier )
               end if
            end if
         end do

! copy local data
         displs = mytask*recvcnt
         do i=1,sendcnt
            recvbuf(displs+i) = sendbuf(i)
         enddo

! wait for final data
         do i=1,min(count,preposts)
            call mpi_wait (rcvid(i), status, ier)
         enddo

      else

         if (sendcnt > 0) then
            call mpi_recv  ( hs, 1, MPI_INTEGER, root, mtag, comm, &
                             status, ier )
            call mpi_rsend ( sendbuf, sendcnt, sendtype, root, mtag, &
                             comm, ier )
         end if

      endif
      call CheckMPIReturn(subName,ier)

   else
 
      call mpi_gather (sendbuf, sendcnt, sendtype, &
                       recvbuf, recvcnt, recvtype, &
                       root, comm, ier)
      call CheckMPIReturn(subName,ier)

   endif

   return

   end subroutine pio_fc_gather_offset

end module pio_support
