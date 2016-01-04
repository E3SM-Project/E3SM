module spmd_utils

!----------------------------------------------------------------------- 
! 
! Purpose: This module is responsible for miscellaneous SPMD utilities
!          and information that are shared between dynamics and 
!          physics packages.  
! 
! Author:
!   Original routines:  CMS
!   Module:             T. Henderson, December 2003
!   swap routines:      P. Worley
!   fc routines:        P. Worley
!   SMP node id logic:  P. Worley
!
! $Id$
! 
!-----------------------------------------------------------------------

!
! Performance bug work around for Gemini interconnect
!
#ifdef _NO_MPI_RSEND
#define mpi_rsend mpi_send
#define mpi_irsend mpi_isend
#endif

!-----------------------------------------------------------------------
!- use statements ------------------------------------------------------
!-----------------------------------------------------------------------
   use shr_kind_mod,     only: r8 => shr_kind_r8
   use cam_abortutils,   only: endrun

#if ( defined SPMD )
   use mpishorthand, only: mpiint, mpii8, mpichar, mpilog, mpipk,      &
                           mpic16, mpir8, mpir4, mpicom, mpimax
#endif
   use cam_logfile,  only: iulog

!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
   include 'mpif.h'          
   private                   ! Make the default access private
   save
!
! Forward from mpishorthand.F with the idea of phasing out use of and removing that file
!
#ifndef SPMD
   integer :: mpi_status_ignore     ! Needs to be defined in mpi-serial
   integer :: mpir8
#endif
!
!  Forward these from mpif.h (or mpi.mod), the idea being that this should
!  be the only module that uses mpi directly, the rest of cam should use spmd_utils
!
   public :: mpi_max_processor_name,                     &
             mpi_integer, mpi_integer8, mpi_character,   &
             mpi_logical, mpi_real8, mpi_real4,          &
             mpi_complex16,                              &
             mpi_packed, mpi_max, mpi_min,               &
             mpi_comm_null, mpi_group_null,              &
             mpi_undefined, mpi_status_size, mpi_success,&
             mpi_status_ignore, mpi_double_precision, mpi_sum, mpir8





!-----------------------------------------------------------------------
! Public interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
   public pair      ! $$$here...  originally from eul|sld/spmd_dyn
   public ceil2     ! $$$here...  originally from eul|sld/spmd_dyn
   public spmdinit
   public spmd_utils_readnl
#if ( defined SPMD )
   public swapm
   public fc_gatherv
   public fc_gathervr4
   public fc_gathervint
   public fc_gathervc
   public altalltoallv
#endif

!-----------------------------------------------------------------------
! Public data ----------------------------------------------------------
!-----------------------------------------------------------------------
! physics-motivated dynamics decomposition request
   logical, parameter :: def_mirror = .false.                 ! default
   logical, public    :: phys_mirror_decomp_req = def_mirror 
                         ! flag indicating whether latitudes and their
                         ! reflections across the equator should be 
                         ! assigned to consecutive processes

#if (defined SPMD)
   public :: mpicom
   public :: mpichar
#else
   integer, public              :: mpicom
   integer, public              :: mpichar
#endif
   logical, public              :: masterproc
   integer, public              :: masterprocid
   integer, public              :: iam
   integer, public              :: npes
   integer, public              :: nsmps
   integer, allocatable, public :: proc_smp_map(:)
   integer, parameter           :: DEFAULT_MASTERPROC=0 
                                      ! the value of iam which is assigned 
                                      ! the masterproc duties

!-----------------------------------------------------------------------
! Private data ---------------------------------------------------------
!-----------------------------------------------------------------------
! Swap communication protocol options (reduced set):
!  3, 5:                  nonblocking send
!  2, 3, 4, 5:            nonblocking receive
!  4, 5:                  ready send
   integer, private, parameter :: min_comm_protocol =  2
   integer, private, parameter :: max_comm_protocol =  5
   integer, private, parameter :: def_comm_protocol =  4        ! default
   integer, public :: swap_comm_protocol = def_comm_protocol

! Swap communication maximum request count:
! = -1,0: do not limit number of outstanding send/receive requests
!    > 0: do not allow more than swap_comm_maxreq outstanding
!         nonblocking send requests or nonblocking receive requests
   integer, private, parameter :: def_comm_maxreq = 128        ! default
   integer, public :: swap_comm_maxreq = def_comm_maxreq

! Flow-controlled gather option:
!   < 0: use MPI_Gather
!  >= 0: use point-to-point with handshaking messages and 
!        preposting receive requests up to 
!         min(max(1,fc_gather_flow_cntl),max_gather_block_size) 
!        ahead
   integer, private, parameter :: max_gather_block_size = 64 ! max and default
   integer, public :: fc_gather_flow_cntl = max_gather_block_size

!-----------------------------------------------------------------------
! Subroutines and functions --------------------------------------------
!-----------------------------------------------------------------------
contains

!========================================================================

   integer function pair(np,p,k)

      integer np,p,k,q
      q = ieor(p,k)
      if(q.gt.np-1) then
         pair = -1
      else
         pair = q
      endif
      return

   end function pair

!========================================================================

  integer function ceil2(n)
     integer n,p
     p=1
     do while(p.lt.n)
        p=p*2
     enddo
     ceil2=p
     return
  end function ceil2

!========================================================================
  
  subroutine spmdinit( mpicom_atm )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: MPI initialization routine:  
    ! 
    ! Method: get number of cpus, processes, tids, etc
    !         dynamics and physics decompositions are set up later
    ! 
    ! Author: CCM Core Group
    ! 
    !-----------------------------------------------------------------------

    implicit none
    integer, intent(in) :: mpicom_atm

#if ( defined SPMD )
    !
    ! Local workspace
    !
    integer i,j,c             ! indices
    integer npthreads         ! thread status
    integer ier               ! return error status    
    integer length            ! length of name
    integer max_len           ! maximum name length
    integer, allocatable :: lengths(:)! max lengths of names for use in gatherv
    integer, allocatable :: displs(:) ! offsets for use in gatherv
    logical done
    character, allocatable                             :: proc_name(:)  ! processor name, this task
    character, allocatable                             :: proc_names(:) ! processor names, all tasks
    character(len=mpi_max_processor_name)              :: tmp_name      ! temporary storage
    character(len=mpi_max_processor_name), allocatable :: smp_names(:)  ! SMP name
    logical mpi_running       ! returned value indicates if MPI_INIT has been called

    !---------------------------------------------------------------------------
    !
    ! Determine CAM MPI communicator group
    !
    mpicom  = mpicom_atm
    !
    ! Set mpishorthand variables.  Need to set as variables rather than parameters since
    ! some MPI implementations set values for MPI tags at run time
    !
    mpiint  = mpi_integer
    mpii8   = mpi_integer8
    mpichar = mpi_character
    mpilog  = mpi_logical
    mpir4   = mpi_real4
    mpir8   = mpi_real8
    mpic16  = mpi_complex16
    mpipk   = mpi_packed
    mpimax  = mpi_max
    !
    ! Get my id  
    !
    call mpi_comm_rank (mpicom, iam, ier) 
    masterprocid = DEFAULT_MASTERPROC
    if (iam == DEFAULT_MASTERPROC) then 
       masterproc = .true.
    else
       masterproc = .false.
    end if
    !
    ! Get number of processors
    !
    max_len = mpi_max_processor_name
    call mpi_comm_size (mpicom, npes, ier)
    allocate ( displs(npes) )
    allocate ( lengths(npes) )
    allocate ( proc_name(max_len) )
    allocate ( proc_names(max_len*npes) )
 
    !
    ! Get processor names and send to root. 
    !
    call mpi_get_processor_name (tmp_name, length, ier)
    proc_name(:) = ' '
    do i = 1, length
       proc_name(i) = tmp_name(i:i)
    end do

    proc_names(:) = ' '
    lengths(:) = max_len
    do i=1,npes
       displs(i) = (i-1)*max_len
    enddo
    call fc_gathervc (proc_name,  max_len, mpichar, &
                      proc_names, lengths, displs, mpichar, &
                      0, mpicom, flow_cntl=-1)
    if (masterproc) then
       write(iulog,*) npes, 'pes participating in computation'
       write(iulog,*) '-----------------------------------'
       write(iulog,*) 'TASK#  NAME'
       do i=0,min(npes-1,256)  ! dont print too many of these
          do c=1,max_len
             tmp_name(c:c) = proc_names(i*max_len+c)
          enddo
          write(iulog,'(i3,2x,a)') i,trim(tmp_name)
       end do
       if(npes-1>256) then
          write(iulog,*) '... list truncated at 256'
       end if
    end if
    !
    ! Identify SMP nodes and process/SMP mapping.
    ! (Assume that processor names are SMP node names on SMP clusters.)
    !
    allocate ( proc_smp_map(0:npes-1) )
    if (masterproc) then
       allocate ( smp_names(0:npes-1) )
       smp_names(:) = ' '
       proc_smp_map(:) = -1
       !
       nsmps = 1
       do c=1,max_len
          tmp_name(c:c) = proc_names(c)
       enddo
       smp_names(0) = trim(tmp_name)
       proc_smp_map(0) = 0
       !
       do i=1,npes-1
          do c=1,max_len
             tmp_name(c:c) = proc_names(i*max_len+c)
          enddo

          j = 0
          done = .false.
          do while ((.not. done) .and. (j < nsmps))
             if (smp_names(j) .eq. trim(tmp_name)) then
                proc_smp_map(i) = j
                done = .true.
             endif
             j = j + 1
          enddo

          if (.not. done) then
             smp_names(nsmps) = trim(tmp_name)
             proc_smp_map(i) = nsmps
             nsmps = nsmps + 1
          endif

       enddo
       deallocate(smp_names)
    endif
    call mpibcast(nsmps, 1, mpiint, 0, mpicom)
    call mpibcast(proc_smp_map, npes, mpiint, 0, mpicom)
    !
    deallocate(displs)
    deallocate(lengths)
    deallocate(proc_name)
    deallocate(proc_names)

#else
    !	
    ! spmd is not defined
    !
    mpicom = mpicom_atm
    iam = 0
    masterprocid = 0
    masterproc = .true.
    npes = 1
    nsmps = 1
    allocate ( proc_smp_map(0:0) )
    proc_smp_map(:) = -1

#endif   

  end subroutine spmdinit

#if (defined SPMD)
!
!========================================================================
!
   subroutine swapm (steps, nprocs, swapids,               &
                     sndbuf, sbuf_siz, sndlths, sdispls,   &
                     rcvbuf, rbuf_siz, rcvlths, rdispls,   &
                     comm, comm_protocol, comm_maxreq      )

!----------------------------------------------------------------------- 
! 
! Purpose: 
!   Reduced version of original swapm (for swap of multiple messages 
!   using MPI point-to-point routines), more efficiently implementing a 
!   subset of the swap protocols.
! 
! Method: 
! comm_protocol:
!  = 3 or 5: use nonblocking send
!  = 2 or 4: use blocking send
!  = 4 or 5: use handshaking protocol
! comm_maxreq:
!  =-1,0: do not limit number of outstanding send/receive requests
!     >0: do not allow more than min(comm_maxreq, steps) outstanding
!         nonblocking send requests or nonblocking receive requests
!
! Author of original version:  P. Worley
! Ported to CAM: P. Worley, December 2003
! Simplified version: P. Worley, October, 2008
! 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
   implicit none
!---------------------------Input arguments--------------------------
!
   integer, intent(in)   :: steps              ! number of swaps to initiate
   integer, intent(in)   :: nprocs             ! size of communicator
   integer, intent(in)   :: sbuf_siz           ! size of send buffer
   integer, intent(in)   :: rbuf_siz           ! size of receive buffer
   integer, intent(in)   :: swapids(steps)     ! MPI process id of swap partners

   integer, intent(in)   :: sndlths(0:nprocs-1)! length of outgoing message
   integer, intent(in)   :: sdispls(0:nprocs-1)! offset from beginning of send
                                               !  buffer where outgoing messages
                                               !  should be sent from
   integer, intent(in)   :: rcvlths(0:nprocs-1)! length of incoming messages
   integer, intent(in)   :: rdispls(0:nprocs-1)! offset from beginning of receive 
                                               !  buffer where incoming messages
                                               !  should be placed
   real(r8), intent(in)  :: sndbuf(sbuf_siz)   ! outgoing message buffer
   real(r8), intent(out) :: rcvbuf(rbuf_siz)   ! incoming message buffer

   integer, intent(in)   :: comm               ! MPI communicator
   integer, intent(in)   :: comm_protocol      ! swap_comm protocol
   integer, intent(in)   :: comm_maxreq        ! maximum number of outstanding 
                                               !  nonblocking requests

!
!---------------------------Local workspace-----------------------------
!
   integer :: p                                ! process index
   integer :: istep                            ! loop index
   integer :: offset_s                         ! index of message beginning in 
                                               !  send buffer
   integer :: offset_r                         ! index of message beginning in 
                                               !  receive buffer
   integer :: sndids(steps)                    ! send request ids
   integer :: rcvids(steps)                    ! receive request ids
   integer :: hs_rcvids(steps)                 ! handshake receive request ids

   integer :: maxreq, maxreqh                  ! maximum number of outstanding 
                                               !  nonblocking requests (and half)
   integer :: hs_s, hs_r(steps)                ! handshake variables (send/receive)
   integer :: rstep                            ! "receive" step index

   logical :: handshake, sendd                 ! protocol option flags

   integer :: ier                              ! return error status    
   integer :: status(MPI_STATUS_SIZE)          ! MPI status 
!
!-------------------------------------------------------------------------------------
!
   if (steps .eq. 0) return

   ! identify communication protocol
   if ((comm_protocol < 2) .or. (comm_protocol > 5)) then
      sendd = .true.
      handshake = .true.
   else
      if ((comm_protocol .eq. 4) .or. (comm_protocol .eq. 5)) then
         handshake = .true.
      else
         handshake = .false.
      endif

      if ((comm_protocol .eq. 2) .or. (comm_protocol .eq. 4)) then
         sendd = .true.
      else
         sendd = .false.
      endif
   endif

   ! identify maximum number of outstanding nonblocking requests to permit
   if (steps .eq. 1) then
      maxreq  = 1
      maxreqh = 1
   else
      if (comm_maxreq >= -1) then
         maxreq = comm_maxreq
      else
         maxreq = steps
      endif

      if ((maxreq .le. steps) .and. (maxreq > 0)) then
         if (maxreq > 1) then
            maxreqh = maxreq/2
         else
            maxreq  = 2
            maxreqh = 1
         endif
      else
         maxreq  = steps
         maxreqh = steps
      endif
   endif

! Four protocol options:
!  (1) handshaking + blocking sends
   if ((handshake) .and. (sendd)) then

      ! Initialize handshake variable
      hs_s = 1

      ! Post initial handshake receive requests
      do istep=1,maxreq
         p = swapids(istep)
         if (sndlths(p) > 0) then
            call mpi_irecv( hs_r(istep), 1, mpiint, p, iam, comm, &
                            hs_rcvids(istep), ier )
         endif
      enddo

      ! Post initial receive requests
      do istep=1,maxreq
         p = swapids(istep)
         if (rcvlths(p) > 0) then
            offset_r = rdispls(p)+1
            call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, p, &
                            comm, rcvids(istep), ier )
            call mpi_send ( hs_s, 1, mpiint, p, p, comm, ier )
         endif
      enddo
      rstep = maxreq

      ! Send (and start receiving) data 
      do istep=1,steps
         p = swapids(istep)

         ! Submit new rsend request
         if (sndlths(p) > 0) then
            offset_s = sdispls(p)+1
            call mpi_wait  ( hs_rcvids(istep), MPI_STATUS_IGNORE, ier )
            call mpi_rsend ( sndbuf(offset_s), sndlths(p), mpir8, p, iam, &
                             comm, ier )
         endif

         if (istep > maxreqh) then

            ! Wait for oldest irecv request to complete
            p = swapids(istep-maxreqh)
            if (rcvlths(p) > 0) then
               call mpi_wait( rcvids(istep-maxreqh), status, ier )
            endif

            if (rstep < steps) then
               rstep = rstep + 1
               p = swapids(rstep)

               ! Submit a new handshake irecv request
               if (sndlths(p) > 0) then
                  call mpi_irecv( hs_r(rstep), 1, mpiint, p, iam, comm, &
                                  hs_rcvids(rstep), ier )
               endif

               ! Submit a new irecv request
               if (rcvlths(p) > 0) then
                  offset_r = rdispls(p)+1
                  call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, p, &
                                  comm, rcvids(rstep), ier )
                  call mpi_send ( hs_s, 1, mpiint, p, p, comm, ier )
               endif
            endif

         endif
!
      enddo

      ! wait for rest of receive requests to complete
      do istep=steps-maxreqh+1,steps
         p = swapids(istep)
         if (rcvlths(p) > 0) then
            call mpi_wait( rcvids(istep), status, ier )
         endif
      enddo

!  (2) handshaking + nonblocking sends
   elseif ((handshake) .and. (.not. sendd)) then

      ! Initialize handshake variable
      hs_s = 1

      ! Post initial handshake receive requests
      do istep=1,maxreq
         p = swapids(istep)
         if (sndlths(p) > 0) then
            call mpi_irecv( hs_r(istep), 1, mpiint, p, iam, comm, &
                            hs_rcvids(istep), ier )
         endif
      enddo

      ! Post initial receive requests
      do istep=1,maxreq
         p = swapids(istep)
         if (rcvlths(p) > 0) then
            offset_r = rdispls(p)+1
            call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, p, &
                            comm, rcvids(istep), ier )
            call mpi_send ( hs_s, 1, mpiint, p, p, comm, ier )
         endif
      enddo
      rstep = maxreq

      ! Send (and start receiving) data 
      do istep=1,steps
         p = swapids(istep)

         ! Submit new irsend request
         if (sndlths(p) > 0) then
            offset_s = sdispls(p)+1
            call mpi_wait  ( hs_rcvids(istep), MPI_STATUS_IGNORE, ier )
            call mpi_irsend( sndbuf(offset_s), sndlths(p), mpir8, p, iam, &
                             comm, sndids(istep), ier )
         endif

         if (istep > maxreqh) then

            ! Wait for oldest irecv request to complete
            p = swapids(istep-maxreqh)
            if (rcvlths(p) > 0) then
               call mpi_wait( rcvids(istep-maxreqh), status, ier )
            endif

            if (rstep < steps) then
               rstep = rstep + 1
               p = swapids(rstep)

               ! Submit a new handshake irecv request
               if (sndlths(p) > 0) then
                  call mpi_irecv( hs_r(rstep), 1, mpiint, p, iam, comm, &
                                  hs_rcvids(rstep), ier )
               endif

               ! Submit a new irecv request
               if (rcvlths(p) > 0) then
                  offset_r = rdispls(p)+1
                  call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, p, &
                                  comm, rcvids(rstep), ier )
                  call mpi_send ( hs_s, 1, mpiint, p, p, comm, ier )
               endif
            endif

            ! Wait for outstanding i(r)send request to complete
            p = swapids(istep-maxreqh)
            if (sndlths(p) > 0) then
               call mpi_wait( sndids(istep-maxreqh), status, ier )
            endif

         endif

      enddo

      ! wait for rest of send and receive requests to complete
      do istep=steps-maxreqh+1,steps
         p = swapids(istep)
         if (rcvlths(p) > 0) then
            call mpi_wait( rcvids(istep), status, ier )
         endif
         if (sndlths(p) > 0) then
            call mpi_wait( sndids(istep), status, ier )
         endif
      enddo

!  (3) no handshaking + blocking sends
   elseif ((.not. handshake) .and. (sendd)) then

      ! Post receive requests
      do istep=1,maxreq
         p = swapids(istep)
         if (rcvlths(p) > 0) then
            offset_r = rdispls(p)+1
            call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, p, &
                            comm, rcvids(istep), ier )
         endif
      enddo
      rstep = maxreq

      ! Send (and start receiving) data 
      do istep=1,steps
         p = swapids(istep)

         ! Submit new send request
         if (sndlths(p) > 0) then
            offset_s = sdispls(p)+1
            call mpi_send( sndbuf(offset_s), sndlths(p), mpir8, p, iam, &
                           comm, ier )
         endif

         if (istep > maxreqh) then

            ! Wait for oldest irecv request to complete
            p = swapids(istep-maxreqh)
            if (rcvlths(p) > 0) then
               call mpi_wait( rcvids(istep-maxreqh), status, ier )
            endif

            ! Submit a new irecv request
            if (rstep < steps) then
               rstep = rstep + 1
               p = swapids(rstep)
               if (rcvlths(p) > 0) then
                  offset_r = rdispls(p)+1
                  call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, p, &
                                  comm, rcvids(rstep), ier )
               endif
            endif

         endif

      enddo

      ! wait for rest of send and receive requests to complete
      do istep=steps-maxreqh+1,steps
         p = swapids(istep)
         if (rcvlths(p) > 0) then
            call mpi_wait( rcvids(istep), status, ier )
         endif
      enddo

!  (4) no handshaking + nonblocking sends
   elseif ((.not. handshake) .and. (.not. sendd)) then

      ! Post receive requests
      do istep=1,maxreq
         p = swapids(istep)
         if (rcvlths(p) > 0) then
            offset_r = rdispls(p)+1
            call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, p, &
                            comm, rcvids(istep), ier )
         endif
      enddo
      rstep = maxreq

      ! Send (and start receiving) data 
      do istep=1,steps
         p = swapids(istep)

         ! Submit new isend request
         if (sndlths(p) > 0) then
            offset_s = sdispls(p)+1
            call mpi_isend( sndbuf(offset_s), sndlths(p), mpir8, p, iam, &
                            comm, sndids(istep), ier )
         endif

         if (istep > maxreqh) then

            ! Wait for oldest irecv request to complete
            p = swapids(istep-maxreqh)
            if (rcvlths(p) > 0) then
               call mpi_wait( rcvids(istep-maxreqh), status, ier )
            endif

            ! Submit a new irecv request
            if (rstep < steps) then
               rstep = rstep + 1
               p = swapids(rstep)
               if (rcvlths(p) > 0) then
                  offset_r = rdispls(p)+1
                  call mpi_irecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, p, &
                                  comm, rcvids(rstep), ier )
               endif
            endif

            ! Wait for outstanding i(r)send request to complete
            p = swapids(istep-maxreqh)
            if (sndlths(p) > 0) then
               call mpi_wait( sndids(istep-maxreqh), status, ier )
            endif

         endif

      enddo

      ! wait for rest of send and receive requests to complete
      do istep=steps-maxreqh+1,steps
         p = swapids(istep)
         if (rcvlths(p) > 0) then
            call mpi_wait( rcvids(istep), status, ier )
         endif
         if (sndlths(p) > 0) then
            call mpi_wait( sndids(istep), status, ier )
         endif
      enddo

   endif

   return

   end subroutine swapm
!
!========================================================================

!----------------------------------------------------------------------- 
! 
! Purpose: gather collective with additional flow control, so as to 
!          be more robust when used with high process counts. 
!          If flow_cntl optional parameter 
!           < 0: use MPI_Gather
!          >= 0: use point-to-point with handshaking messages and 
!                preposting receive requests up to 
!                 min(max(1,flow_cntl),max_gather_block_size) 
!                ahead if optional flow_cntl parameter is present.
!                Otherwise, fc_gather_flow_cntl is used in its place.
!          Default value is 64.
! 
! Entry points:
!      fc_gatherv       functionally equivalent to mpi_gatherv
!      fc_gathervr4     functionally equivalent to mpi_gatherv for real*4 data
!      fc_gathervint    functionally equivalent to mpi_gatherv for integer data
!      fc_gathervc      functionally equivalent to mpi_gatherv for character data
!
! Author: P. Worley
!-----------------------------------------------------------------------

!
!========================================================================
!
   subroutine fc_gatherv (sendbuf, sendcnt, sendtype, &
                         recvbuf, recvcnts, displs, recvtype, &
                         root, comm, flow_cntl )
!
! Collects different messages from each process on masterproc
!
   use shr_kind_mod,   only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils, only: endrun
   use cam_logfile,    only: iulog

#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   real (r8), intent(in)  :: sendbuf(*)
   real (r8), intent(out) :: recvbuf(*)
   integer, intent(in) :: displs(*)
   integer, intent(in) :: sendcnt
   integer, intent(in) :: sendtype
   integer, intent(in) :: recvcnts(*)
   integer, intent(in) :: recvtype
   integer, intent(in) :: root
   integer, intent(in) :: comm
   integer, optional, intent(in) :: flow_cntl

   real (r8) :: signal
   logical fc_gather         ! use explicit flow control?
   integer gather_block_size ! number of preposted receive requests

   integer :: mytid, mysize, mtag, p, q, i, count
   integer :: preposts, head, tail
   integer :: rcvid(max_gather_block_size)
   integer :: status(MPI_STATUS_SIZE)
   integer ier               ! MPI error code

   if ( present(flow_cntl) ) then
      if (flow_cntl >= 0) then
         gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
         fc_gather = .true.
      else
         fc_gather = .false.
      endif
   else
      if (fc_gather_flow_cntl >= 0) then
         gather_block_size = min(max(1,fc_gather_flow_cntl),max_gather_block_size)
         fc_gather = .true.
      else
         fc_gather = .false.
      endif
   endif

   if (fc_gather) then
 
#if defined( WRAP_MPI_TIMING )
      call t_startf ('fc_gatherv_r8')
#endif
      call mpi_comm_rank (comm, mytid, ier)
      call mpi_comm_size (comm, mysize, ier)
      mtag = 0
      if (root .eq. mytid) then

! prepost gather_block_size irecvs, and start receiving data
         preposts = min(mysize-1, gather_block_size)
         head = 0
         count = 0
         do p=0, mysize-1
            if (p .ne. root) then
               q = p+1
               if (recvcnts(q) > 0) then
                  count = count + 1
                  if (count > preposts) then
                     tail = mod(head,preposts) + 1
                     call mpi_wait (rcvid(tail), status, ier)
                  end if
                  head = mod(head,preposts) + 1
                  call mpi_irecv ( recvbuf(displs(q)+1), recvcnts(q), &
                                   recvtype, p, mtag, comm, rcvid(head), &
                                   ier )
                  call mpi_send ( signal, 1, mpir8, p, mtag, comm, ier )
               end if
            end if
         end do

! copy local data
         q = mytid+1
         do i=1,sendcnt
            recvbuf(displs(q)+i) = sendbuf(i)
         enddo

! wait for final data
         do i=1,min(count,preposts)
            call mpi_wait (rcvid(i), status, ier)
         enddo

      else

         if (sendcnt > 0) then
            call mpi_recv ( signal, 1, mpir8, root, mtag, comm, &
                            status, ier )
            call mpi_rsend ( sendbuf, sendcnt, sendtype, root, mtag, &
                             comm, ier )
         end if

      endif
      if (ier /= mpi_success) then
         write(iulog,*)'fc_gatherv_r8 failed ier=',ier
         call endrun
      end if
#if defined( WRAP_MPI_TIMING )
      call t_stopf ('fc_gatherv_r8')
#endif

   else
 
#if defined( WRAP_MPI_TIMING )
      call t_startf ('mpi_gatherv')
#endif
      call mpi_gatherv (sendbuf, sendcnt, sendtype, &
                        recvbuf, recvcnts, displs, recvtype, &
                        root, comm, ier)
      if (ier /= mpi_success) then
         write(iulog,*)'mpi_gatherv failed ier=',ier
         call endrun
      end if
#if defined( WRAP_MPI_TIMING )
      call t_stopf ('mpi_gatherv')
#endif

   endif

   return
   end subroutine fc_gatherv
!
!========================================================================
!
   subroutine fc_gathervr4 (sendbuf, sendcnt, sendtype, &
                           recvbuf, recvcnts, displs, recvtype, &
                           root, comm, flow_cntl )
!
! Collects different messages from each process on masterproc
!
   use shr_kind_mod,   only: r4 => shr_kind_r4, r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils, only: endrun
   use cam_logfile,    only: iulog

#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   real (r4), intent(in)  :: sendbuf(*)
   real (r4), intent(out) :: recvbuf(*)
   integer, intent(in) :: displs(*)
   integer, intent(in) :: sendcnt
   integer, intent(in) :: sendtype
   integer, intent(in) :: recvcnts(*)
   integer, intent(in) :: recvtype
   integer, intent(in) :: root
   integer, intent(in) :: comm
   integer, optional, intent(in) :: flow_cntl

   real (r8) :: signal
   logical fc_gather         ! use explicit flow control?
   integer gather_block_size ! number of preposted receive requests

   integer :: mytid, mysize, mtag, p, q, i, count
   integer :: preposts, head, tail
   integer :: rcvid(max_gather_block_size)
   integer :: status(MPI_STATUS_SIZE)
   integer ier               ! MPI error code

   if ( present(flow_cntl) ) then
      if (flow_cntl >= 0) then
         gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
         fc_gather = .true.
      else
         fc_gather = .false.
      endif
   else
      if (fc_gather_flow_cntl >= 0) then
         gather_block_size = min(max(1,fc_gather_flow_cntl),max_gather_block_size)
         fc_gather = .true.
      else
         fc_gather = .false.
      endif
   endif

   if (fc_gather) then
 
#if defined( WRAP_MPI_TIMING )
      call t_startf ('fc_gatherv_r4')
#endif
      call mpi_comm_rank (comm, mytid, ier)
      call mpi_comm_size (comm, mysize, ier)
      mtag = 0
      if (root .eq. mytid) then

! prepost gather_block_size irecvs, and start receiving data
         preposts = min(mysize-1, gather_block_size)
         head = 0
         count = 0
         do p=0, mysize-1
            if (p .ne. root) then
               q = p+1
               if (recvcnts(q) > 0) then
                  count = count + 1
                  if (count > preposts) then
                     tail = mod(head,preposts) + 1
                     call mpi_wait (rcvid(tail), status, ier)
                  end if
                  head = mod(head,preposts) + 1
                  call mpi_irecv ( recvbuf(displs(q)+1), recvcnts(q), &
                                   recvtype, p, mtag, comm, rcvid(head), &
                                   ier )
                  call mpi_send ( signal, 1, mpir8, p, mtag, comm, ier )
               end if
            end if
         end do

! copy local data
         q = mytid+1
         do i=1,sendcnt
            recvbuf(displs(q)+i) = sendbuf(i)
         enddo

! wait for final data
         do i=1,min(count,preposts)
            call mpi_wait (rcvid(i), status, ier)
         enddo

      else

         if (sendcnt > 0) then
            call mpi_recv ( signal, 1, mpir8, root, mtag, comm, &
                            status, ier )
            call mpi_rsend ( sendbuf, sendcnt, sendtype, root, mtag, &
                             comm, ier )
          end if

      endif
      if (ier /= mpi_success) then
         write(iulog,*)'fc_gatherv_r4 failed ier=',ier
         call endrun
      end if
#if defined( WRAP_MPI_TIMING )
      call t_stopf ('fc_gatherv_r4')
#endif

   else
 
#if defined( WRAP_MPI_TIMING )
      call t_startf ('mpi_gatherv')
#endif
      call mpi_gatherv (sendbuf, sendcnt, sendtype, &
                        recvbuf, recvcnts, displs, recvtype, &
                        root, comm, ier)
      if (ier /= mpi_success) then
         write(iulog,*)'mpi_gatherv failed ier=',ier
         call endrun
      end if
#if defined( WRAP_MPI_TIMING )
      call t_stopf ('mpi_gatherv')
#endif

   endif

   return
   end subroutine fc_gathervr4
!
!========================================================================
!
   subroutine fc_gathervint (sendbuf, sendcnt, sendtype, &
                            recvbuf, recvcnts, displs, recvtype, &
                            root, comm, flow_cntl )
!
! Collects different messages from each process on masterproc
!
   use shr_kind_mod,   only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils, only: endrun
   use cam_logfile,    only: iulog

#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   integer, intent(in)  :: sendbuf(*)
   integer, intent(out) :: recvbuf(*)
   integer, intent(in) :: displs(*)
   integer, intent(in) :: sendcnt
   integer, intent(in) :: sendtype
   integer, intent(in) :: recvcnts(*)
   integer, intent(in) :: recvtype
   integer, intent(in) :: root
   integer, intent(in) :: comm
   integer, optional, intent(in) :: flow_cntl

   real (r8) :: signal
   logical fc_gather         ! use explicit flow control?
   integer gather_block_size ! number of preposted receive requests

   integer :: mytid, mysize, mtag, p, q, i, count
   integer :: preposts, head, tail
   integer :: rcvid(max_gather_block_size)
   integer :: status(MPI_STATUS_SIZE)
   integer ier               ! MPI error code

   if ( present(flow_cntl) ) then
      if (flow_cntl >= 0) then
         gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
         fc_gather = .true.
      else
         fc_gather = .false.
      endif
   else
      if (fc_gather_flow_cntl >= 0) then
         gather_block_size = min(max(1,fc_gather_flow_cntl),max_gather_block_size)
         fc_gather = .true.
      else
         fc_gather = .false.
      endif
   endif

   if (fc_gather) then
 
#if defined( WRAP_MPI_TIMING )
      call t_startf ('fc_gatherv_int')
#endif
      call mpi_comm_rank (comm, mytid, ier)
      call mpi_comm_size (comm, mysize, ier)
      mtag = 0
      if (root .eq. mytid) then

! prepost gather_block_size irecvs, and start receiving data
         preposts = min(mysize-1, gather_block_size)
         head = 0
         count = 0
         do p=0, mysize-1
            if (p .ne. root) then
               q = p+1
               if (recvcnts(q) > 0) then
                  count = count + 1
                  if (count > preposts) then
                     tail = mod(head,preposts) + 1
                     call mpi_wait (rcvid(tail), status, ier)
                  end if
                  head = mod(head,preposts) + 1
                  call mpi_irecv ( recvbuf(displs(q)+1), recvcnts(q), &
                                   recvtype, p, mtag, comm, rcvid(head), &
                                   ier )
                  call mpi_send ( signal, 1, mpir8, p, mtag, comm, ier )
               end if
            end if
         end do

! copy local data
         q = mytid+1
         do i=1,sendcnt
            recvbuf(displs(q)+i) = sendbuf(i)
         enddo

! wait for final data
         do i=1,min(count,preposts)
            call mpi_wait (rcvid(i), status, ier)
         enddo

      else

         if (sendcnt > 0) then
            call mpi_recv ( signal, 1, mpir8, root, mtag, comm, &
                            status, ier )
            call mpi_rsend ( sendbuf, sendcnt, sendtype, root, mtag, &
                             comm, ier )
          end if

      endif
      if (ier /= mpi_success) then
         write(iulog,*)'fc_gatherv_int failed ier=',ier
         call endrun
      end if
#if defined( WRAP_MPI_TIMING )
      call t_stopf ('fc_gatherv_int')
#endif

   else
 
#if defined( WRAP_MPI_TIMING )
      call t_startf ('mpi_gatherv')
#endif
      call mpi_gatherv (sendbuf, sendcnt, sendtype, &
                        recvbuf, recvcnts, displs, recvtype, &
                        root, comm, ier)
      if (ier /= mpi_success) then
         write(iulog,*)'mpi_gatherv failed ier=',ier
         call endrun
      end if
#if defined( WRAP_MPI_TIMING )
      call t_stopf ('mpi_gatherv')
#endif

   endif

   return
   end subroutine fc_gathervint
!
!========================================================================
!
   subroutine fc_gathervc (sendbuf, sendcnt, sendtype, &
                           recvbuf, recvcnts, displs, recvtype, &
                           root, comm, flow_cntl )
!
! Collects different messages from each process on masterproc
!
   use shr_kind_mod,   only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils, only: endrun
   use cam_logfile,    only: iulog

#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   character, intent(in)  :: sendbuf(*)
   character, intent(out) :: recvbuf(*)
   integer, intent(in) :: displs(*)
   integer, intent(in) :: sendcnt
   integer, intent(in) :: sendtype
   integer, intent(in) :: recvcnts(*)
   integer, intent(in) :: recvtype
   integer, intent(in) :: root
   integer, intent(in) :: comm
   integer, optional, intent(in) :: flow_cntl

   real (r8) :: signal
   logical fc_gather         ! use explicit flow control?
   integer gather_block_size ! number of preposted receive requests

   integer :: mytid, mysize, mtag, p, q, i, count
   integer :: preposts, head, tail
   integer :: rcvid(max_gather_block_size)
   integer :: status(MPI_STATUS_SIZE)
   integer ier               ! MPI error code

   if ( present(flow_cntl) ) then
      if (flow_cntl >= 0) then
         gather_block_size = min(max(1,flow_cntl),max_gather_block_size)
         fc_gather = .true.
      else
         fc_gather = .false.
      endif
   else
      if (fc_gather_flow_cntl >= 0) then
         gather_block_size = min(max(1,fc_gather_flow_cntl),max_gather_block_size)
         fc_gather = .true.
      else
         fc_gather = .false.
      endif
   endif

   if (fc_gather) then
 
#if defined( WRAP_MPI_TIMING )
      call t_startf ('fc_gatherv_char')
#endif
      call mpi_comm_rank (comm, mytid, ier)
      call mpi_comm_size (comm, mysize, ier)
      mtag = 0
      if (root .eq. mytid) then

! prepost gather_block_size irecvs, and start receiving data
         preposts = min(mysize-1, gather_block_size)
         head = 0
         count = 0
         do p=0, mysize-1
            if (p .ne. root) then
               q = p+1
               if (recvcnts(q) > 0) then
                  count = count + 1
                  if (count > preposts) then
                     tail = mod(head,preposts) + 1
                     call mpi_wait (rcvid(tail), status, ier)
                  end if
                  head = mod(head,preposts) + 1
                  call mpi_irecv ( recvbuf(displs(q)+1), recvcnts(q), &
                                   recvtype, p, mtag, comm, rcvid(head), &
                                   ier )
                  call mpi_send ( signal, 1, mpir8, p, mtag, comm, ier )
               end if
            end if
         end do

! copy local data
         q = mytid+1
         do i=1,sendcnt
            recvbuf(displs(q)+i) = sendbuf(i)
         enddo

! wait for final data
         do i=1,min(count,preposts)
            call mpi_wait (rcvid(i), status, ier)
         enddo

      else

         if (sendcnt > 0) then
            call mpi_recv ( signal, 1, mpir8, root, mtag, comm, &
                            status, ier )
            call mpi_rsend ( sendbuf, sendcnt, sendtype, root, mtag, &
                             comm, ier )
          end if

      endif
      if (ier /= mpi_success) then
         write(iulog,*)'fc_gatherv_char failed ier=',ier
         call endrun
      end if
#if defined( WRAP_MPI_TIMING )
      call t_stopf ('fc_gatherv_char')
#endif

   else
 
#if defined( WRAP_MPI_TIMING )
      call t_startf ('mpi_gatherv')
#endif
      call mpi_gatherv (sendbuf, sendcnt, sendtype, &
                        recvbuf, recvcnts, displs, recvtype, &
                        root, comm, ier)
      if (ier /= mpi_success) then
         write(iulog,*)'mpi_gatherv failed ier=',ier
         call endrun
      end if
#if defined( WRAP_MPI_TIMING )
      call t_stopf ('mpi_gatherv')
#endif

   endif

   return
   end subroutine fc_gathervc
!
!========================================================================
#endif

!----------------------------------------------------------------------- 
! 
! Purpose: implementations of MPI_Alltoall using different messaging
!          layers and different communication protocols, controlled
!          by option argument:
!  0: use mpi_alltoallv
!  1: use point-to-point MPI-1 two-sided implementation
!  2: use point-to-point MPI-2 one-sided implementation if supported, 
!       otherwise use MPI-1 implementation
!  3: use Co-Array Fortran implementation if supported, 
!       otherwise use MPI-1 implementation
!  otherwise use mpi_sendrecv implementation
! 
! Entry points:
!      altalltoallv
!
! Author: P. Worley
!-----------------------------------------------------------------------

#if (defined SPMD)
!****************************************************************
   subroutine altalltoallv (option, mytid, nprocs, steps, dests, &
                 sendbuf, sbuf_siz, sendcnts, sdispls, sendtype, &
                 recvbuf, rbuf_siz, recvcnts, rdispls, recvtype, &
                 msgtag, pdispls, desttype, recvwin, comm)
!
! All-to-all scatter/gather implemented using Co-Array
! Fortran one-sided commands, MPI-2 one sided commands,
! SWAP module MPI-1 commands, MPI_ALLTOALLV or MPI_SENDRECV.
!
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   integer, intent(in) :: option               ! 0: mpi_alltoallv
                                               ! 1: swap package
                                               ! 2: mpi2 
                                               ! 3: co-array fortran
                                       ! otherwise: sendrecv
   integer, intent(in) :: mytid
   integer, intent(in) :: nprocs
   integer, intent(in) :: steps
   integer, intent(in) :: dests(steps)
   integer, intent(in) :: sbuf_siz
   integer, intent(in) :: sendcnts(0:nprocs-1)
   integer, intent(in) :: sdispls(0:nprocs-1)
   integer, intent(in) :: sendtype
   integer, intent(in) :: rbuf_siz
   integer, intent(in) :: recvcnts(0:nprocs-1)
   integer, intent(in) :: rdispls(0:nprocs-1)
   integer, intent(in) :: recvtype
   integer, intent(in) :: msgtag
   integer, intent(in) :: pdispls(0:nprocs-1)   ! displacement at 
                                                !  destination
   integer, intent(in) :: desttype
   integer, intent(in) :: recvwin
   integer, intent(in) :: comm

#if (defined CAF)
   real (r8), intent(in)  :: sendbuf(sbuf_siz)[*]
   real (r8), intent(out) :: recvbuf(rbuf_siz)[*]

   integer :: istart, iend, jstart, jend
#else
   real (r8), intent(in)  :: sendbuf(sbuf_siz)
   real (r8), intent(out) :: recvbuf(rbuf_siz)
#endif

   integer :: loption          ! local copy of option
   integer :: dest             ! MPI remote process id
   integer :: ier              ! MPI error code
   integer :: i                ! loop index
   integer :: sndids(steps)    ! nonblocking MPI send request ids
   integer :: rcvids(steps)    ! nonblocking MPI recv request ids
   integer :: status(MPI_STATUS_SIZE)
#if ( defined MPI2)
   integer(kind=MPI_ADDRESS_KIND) :: ddispls
#endif

!-----------------------------------------------------------------------
   loption = option

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  using MPI library collective MPI_ALLTOALLV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (loption .eq. 0) then

#if defined( WRAP_MPI_TIMING )
      call t_startf ('mpi_alltoallv')
#endif
      call mpi_alltoallv (sendbuf, sendcnts, sdispls, sendtype, &
                          recvbuf, recvcnts, rdispls, recvtype, &
                          comm, ier)
!
! test for error
      if (ier/=mpi_success) then
         write(iulog,*)'altalltoallv (mpi_alltoallv) failed ier=',ier
         call endrun
      end if
#if defined( WRAP_MPI_TIMING )
      call t_stopf ('mpi_alltoallv')
#endif

   else

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Co-Array Fortran implementation of alltoallv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (loption .eq. 3) then

#if ( defined CAF )
#if defined( WRAP_MPI_TIMING )
         call t_startf ('caf_alltoallv')
#endif
         if (this_image() .ne. (mytid+1)) then
            call endrun('altalltoallv (caf_alltoallv) failed: MPI id .ne. CAF id')
         endif

         call sync_images()

!DIR$ CONCURRENT
         do i = 1, steps
            dest = dests(i)
            if (sendcnts(dest) > 0) then
               istart = sdispls(dest)+1
               iend   = istart+sendcnts(dest)-1
               jstart = pdispls(dest)+1
               jend   = jstart+sendcnts(dest)-1
               recvbuf(jstart:jend)[dest+1] = sendbuf(istart:iend)
            end if
         end do

         call sync_images()
#if defined( WRAP_MPI_TIMING )
         call t_stopf ('caf_alltoallv')
#endif
#else
         loption = -1
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  MPI-2 one-sided implementation of alltoallv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif (loption .eq. 2) then
#ifdef MPI2
#if defined( WRAP_MPI_TIMING )
         call t_startf ('mpi2_alltoallv')
#endif
         call mpi_win_fence(0,recvwin,ier)
         do i=1, steps
            dest = dests(i)
            if (sendcnts(dest) > 0) then
               ddispls = pdispls(dest)
               call mpi_put(sendbuf(sdispls(dest)+1), sendcnts(dest), sendtype, &
                            dest, ddispls, sendcnts(dest), desttype, &
                            recvwin, ier)
            endif
         end do
!
! wait for completion
         call mpi_win_fence(0,recvwin,ier)
         if (ier/=mpi_success) then
            write(iulog,*)'altalltoallv (mpi2_alltoallv) failed ier=',ier
            call endrun
         end if
#if defined( WRAP_MPI_TIMING )
         call t_stopf ('mpi2_alltoallv')
#endif
#else
         loption = -1
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  MPI-1 two-sided implementation of alltoallv
!  using SWAP routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif (loption .eq. 1) then
#if defined( WRAP_MPI_TIMING )
         call t_startf ('swap_alltoallv')
#endif

         call swapm(steps, nprocs, dests,                      &
                    sendbuf, sbuf_siz, sendcnts, sdispls,      &
                    recvbuf, rbuf_siz, recvcnts, rdispls,      &
                    comm, swap_comm_protocol, swap_comm_maxreq )
!
#if defined( WRAP_MPI_TIMING )
         call t_stopf ('swap_alltoallv')
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Anything else defined to be MPI_SENDRECV implementation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else
!
         loption = -1
!
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  MPI_SENDRECV implementation of alltoallv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (loption .eq. -1) then
#if defined( WRAP_MPI_TIMING )
         call t_startf ('mpi1_alltoallv')
#endif
         do i=1, steps
            dest = dests(i)
            call mpi_sendrecv (sendbuf(sdispls(dest)+1), sendcnts(dest), &
                               sendtype, dest, msgtag,                   &
                               recvbuf(rdispls(dest)+1), recvcnts(dest), &
                               recvtype, dest, msgtag,                   &
                               comm, status, ier)
         end do
!
! test for error
         if (ier/=mpi_success) then
            write(iulog,*)'altalltoallv (mpi1_alltoallv) failed ier=',ier
            call endrun
         end if

#if defined( WRAP_MPI_TIMING )
         call t_stopf ('mpi1_alltoallv')
#endif
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Local copy (if necessary)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (sendcnts(mytid) > 0) then
         do i=1,sendcnts(iam)
            recvbuf(rdispls(mytid)+i) = sendbuf(sdispls(mytid)+i)
         enddo
      endif
!
   endif
!
   return
   end subroutine altalltoallv

#endif
   
   subroutine spmd_utils_readnl(nlfile)
!----------------------------------------------------------------------- 
! 
! Purpose: 
!   Read spmd utils namelist to set swap communication protocol options as
!   well as the flow control gather options
! 
! Method: 
! spmd_utils_readnl:
!
! Author of original version:  J. Truesdale
! 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
     use namelist_utils,  only: find_group_name
     use units,           only: getunit, freeunit
     use mpishorthand
     
     implicit none
!---------------------------Input arguments--------------------------
!
     character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

#if ( defined SPMD )
!---------------------------Local variables--------------------------
!
     integer :: unitn, ierr
     character(len=*), parameter :: subname = 'spmd_utils_readnl'
     
     namelist /spmd_utils_nl/ swap_comm_protocol,swap_comm_maxreq,fc_gather_flow_cntl

!-----------------------------------------------------------------------------

     if (masterproc) then
        unitn = getunit()
        open( unitn, file=trim(nlfile), status='old' )
        call find_group_name(unitn, 'spmd_utils_nl', status=ierr)
        if (ierr == 0) then
           read(unitn, spmd_utils_nl, iostat=ierr)
           if (ierr /= 0) then
              call endrun(subname // ':: ERROR reading namelist')
           end if
           write(iulog,*) 'Read in spmd_utils_nl namelist from: ', trim(nlfile)
        end if
        close(unitn)
        call freeunit(unitn)
        
           
        if ((swap_comm_protocol < min_comm_protocol) .or. &
             (swap_comm_protocol > max_comm_protocol)) then
           write(iulog,*)                                        &
                'SPMD_UTILS_READNL:  ERROR:  swap_comm_protocol=', &
                swap_comm_protocol, ' is out of range.'
           write(iulog,*)                                        &
                '  It must be between ', min_comm_protocol,' and ',&
                max_comm_protocol
           write(iulog,*)                                        &
                '  Using default value.'
           swap_comm_protocol = def_comm_protocol
        endif
        
        write(iulog,*) 'SPMD SWAP_COMM OPTIONS: '
        write(iulog,*) '  swap_comm_protocol = ', swap_comm_protocol
        write(iulog,*) '  swap_comm_maxreq   = ', swap_comm_maxreq         
        write(iulog,*) 'SPMD FLOW CONTROL GATHER OPTION: '
        write(iulog,*) '  fc_gather_flow_cntl = ', fc_gather_flow_cntl
     endif
        
     ! Broadcast namelist variables
     call mpibcast (swap_comm_protocol , 1,   mpiint ,  0, mpicom)
     call mpibcast (swap_comm_maxreq   , 1,   mpiint ,  0, mpicom)
     call mpibcast (fc_gather_flow_cntl, 1,   mpiint ,  0, mpicom)
#endif
      
   end subroutine spmd_utils_readnl
 
 end module spmd_utils

