!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$ 
!BOP -------------------------------------------------------------------
!
! !MODULE: m_GlobalSegMapComms - GlobalSegMap Communications Support
!
! !DESCRIPTION:
!
! This module provides communications support for the {\tt GlobalSegMap} 
! datatype.  Both blocking and non-blocking point-to-point communications 
! are supported (i.e., analogues to {\tt MPI\_SEND()/MPI\_RECV()} and 
! {\tt MPI\_ISEND()/MPI\_IRECV()}).  A broadcast method is also supplied.
!
! !INTERFACE:

 module m_GlobalSegMapComms

      implicit none

      private   ! except

! !PUBLIC MEMBER FUNCTIONS:

      public :: send
      public :: recv
      public :: isend
      public :: irecv
      public :: bcast

      interface bcast ; module procedure bcast_ ; end interface
      interface send ; module procedure send_ ; end interface
      interface recv ; module procedure recv_ ; end interface
      interface isend ; module procedure isend_ ; end interface
      interface irecv ; module procedure irecv_ ; end interface

! !REVISION HISTORY:
! 11Aug03 - J.W. Larson <larson@mcs.anl.gov> - initial version
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_GlobalSegMapComms'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: send_ - Point-to-point Blocking Send of a GlobalSegMap
!
! !DESCRIPTION:
! This routine performs a blocking send of a {\tt GlobalSegMap} (the 
! input argument {\tt outgoingGSMap}) to the processor whose rank is 
! {\tt destID} on the MPI communicator associated with the integer 
! handle {\tt comm}.  The input {\tt INTEGER} argument {\tt TagBase} 
! is used to generate tags for the messages associated with this operation; 
! there are two messages involved, so the user should avoid using tag 
! values {\tt TagBase} and {\tt TagBase + 1}.  The success (failure) of 
! this operation is reported in the zero (non-zero) value of the optional 
! {\tt INTEGER} output variable {\tt status}.
!
! !INTERFACE:

 subroutine send_(outgoingGSMap, destID, TagBase, comm, status)

!
! !USES:
!
      use m_mpif90
      use m_die, only : MP_perr_die, die
      use m_stdio

      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_ngseg => ngseg
      use m_GlobalSegMap, only : GlobalSegMap_comp_id => comp_ID
      use m_GlobalSegMap, only : GlobalSegMap_gsize => gsize

      implicit none

! !INPUT PARAMETERS:

  type(GlobalSegMap), intent(IN)  :: outgoingGSMap
  integer,            intent(IN)  :: destID
  integer,            intent(IN)  :: TagBase
  integer,            intent(IN)  :: comm

! !OUTPUT PARAMETERS: 

  integer, optional,  intent(OUT) :: status

! !REVISION HISTORY:
! 13Aug03 - J.W. Larson <larson@mcs.anl.gov> - API and initial version.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::send_'

  integer :: BufferSize, nSegs, gSize, compID, ierr
  integer :: i, Offset
  integer, dimension(:), pointer :: SendBuffer

  if(present(status)) status = 0 ! the success value

       ! First order of business--determine buffer size

  nSegs = GlobalSegMap_ngseg(outgoingGSMap)

  BufferSize = 3 * nSegs + 3 ! add in %comp_id, %ngseg, and %gsize

       ! Next, send the buffer size to destID so it can prepare a
       ! receive buffer of the correct size

  call MPI_SEND(BufferSize, 1, MP_Type(BufferSize), destID, TagBase, &
                comm, ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_, 'Send buffer size failed',ierr)
  endif

       ! Query outgoingGSMap for component ID number and gsize

  compID = GlobalSegMap_comp_id(outgoingGSMap)
  gSize = GlobalSegMap_gsize(outgoingGSMap)

       ! Allocate send buffer

  allocate(SendBuffer(BufferSize), stat=ierr)
  if(ierr /= 0) then
     write(stderr,'(2a,i8)') myname_, &
	  ':: FATAL--allocation of send buffer failed with ierr=',ierr
     call die(myname_)
  endif

       ! Pack send buffer

  SendBuffer(1) = compID
  SendBuffer(2) = nSegs
  SendBuffer(3) = gsize

  Offset = 3

  do i=1,nSegs
     SendBuffer(i+Offset) = outGoingGSMap%start(i)
  end do

  Offset = nSegs + 3

  do i=1,nSegs
     SendBuffer(i+Offset) = outGoingGSMap%length(i)
  end do

  Offset = 2*nSegs + 3

  do i=1,nSegs
     SendBuffer(i+Offset) = outGoingGSMap%pe_loc(i)
  end do

       ! Send the buffer.  And, we should be finished!

  call MPI_SEND(SendBuffer, BufferSize, MP_Type(SendBuffer(1)), &
                destID, TagBase+1, comm, ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_, 'Send segment information buffer failed',ierr)
  endif

       ! Clean up.

  deallocate(SendBuffer, stat=ierr)
  if(ierr /= 0) then
     write(stderr,'(2a,i8)') myname_, &
	  ':: FATAL--deallocation of send buffer failed with ierr=',ierr
     if(present(status)) then
	status = ierr
	return
     else
	call die(myname_)
     endif
  endif

 end subroutine send_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: recv_ - Point-to-point Blocking Receive of a GlobalSegMap
!
! !DESCRIPTION:
! This routine performs a blocking receive of a {\tt GlobalSegMap} (the 
! input argument {\tt incomingGSMap}) from the processor whose rank is 
! {\tt sourceID} on the MPI communicator associated with the integer 
! handle {\tt comm}.  The input {\tt INTEGER} argument {\tt TagBase} 
! is used to generate tags for the messages associated with this operation; 
! there are two messages involved, so the user should avoid using tag 
! values {\tt TagBase} and {\tt TagBase + 1}.  The success (failure) of 
! this operation is reported in the zero (non-zero) value of the optional 
! {\tt INTEGER} output variable {\tt status}.
!
! !INTERFACE:

 subroutine recv_(incomingGSMap, sourceID, TagBase, comm, status)

!
! !USES:
!
      use m_mpif90
      use m_die, only : MP_perr_die, die
      use m_stdio

      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_init => init

      implicit none

! !INPUT PARAMETERS:

  integer,            intent(IN)  :: sourceID
  integer,            intent(IN)  :: TagBase
  integer,            intent(IN)  :: comm

! !OUTPUT PARAMETERS: 

  type(GlobalSegMap), intent(OUT) :: incomingGSMap
  integer, optional,  intent(OUT) :: status

! !REVISION HISTORY:
! 13Aug03 - J.W. Larson <larson@mcs.anl.gov> - API and initial version.

!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::recv_'

  integer :: BufferSize, nSegs, gSize, compID, ierr
  integer :: MPstatus(MP_STATUS_SIZE)
  integer, dimension(:), pointer :: RecvBuffer

  if(present(status)) status = 0 ! the success value

       ! Receive buffer size

  call MPI_RECV(BufferSize, 1, MP_Type(BufferSize), sourceID, TagBase, &
                comm, MPstatus, ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_, 'Receive buffer size failed',ierr)
  endif

       ! Allocate receive buffer

  allocate(RecvBuffer(BufferSize), stat=ierr)
  if(ierr /= 0) then
     write(stderr,'(2a,i8)') myname_, &
	  ':: FATAL--allocation of receive buffer failed with ierr=',ierr
     call die(myname_)
  endif

       ! Receive the buffer...

  call MPI_RECV(RecvBuffer, BufferSize, MP_Type(RecvBuffer(1)), sourceID, &
                TagBase+1, comm, MPstatus, ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_, 'Receive of buffer failed',ierr)
  endif

       ! Check buffer size for compatibility with a GlobalSegMap; i.e., 
       ! is BufferSize-3 a whole number multiple of 3?  If not, something
       ! is wrong and the process should be aborted.

  if(mod(BufferSize-3,3) /= 0) then
     write(stderr,*) myname_, &
	  ':: FATAL BufferSize-3 must be a whole number multiple of 3.',&
	  '  BufferSize=',BufferSize
     if(present(status)) then
	status = -2
     else
	call die(myname_)
     endif
  endif

       ! Create the outgoing GlobalSegMap

  call GlobalSegMap_init(incomingGSMap, RecvBuffer(1), &
                         RecvBuffer(2), &
                         RecvBuffer(3), &
			 RecvBuffer(4:BufferSize))

       ! Clean Up.

  deallocate(RecvBuffer, stat=ierr)
  if(ierr /= 0) then
     write(stderr,'(2a,i8)') myname_, &
	  ':: FATAL--allocation of receive buffer failed with ierr=',ierr
     call die(myname_)
  endif

 end subroutine recv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: isend_ - Point-to-point Non-blocking Send of a GlobalSegMap
!
! !DESCRIPTION:
! This routine performs a non-blocking send of a {\tt GlobalSegMap} (the 
! input argument {\tt outgoingGSMap}) to the processor whose rank is 
! {\tt destID} on the MPI communicator associated with the integer 
! handle {\tt comm}.  The input {\tt INTEGER} argument {\tt TagBase} 
! is used to generate tags for the messages associated with this operation; 
! there are four messages involved, so the user should avoid using tag 
! values {\tt TagBase} and {\tt TagBase + 3}.  All four messages are non-
! blocking, and the request handles for them are returned in the output
! {\tt INTEGER} array {\tt reqHandle}, which can be checked for completion 
! using any of MPI's wait functions.  The success (failure) of 
! this operation is reported in the zero (non-zero) value of the optional 
! {\tt INTEGER} output variable {\tt status}.
!
! {\bf N.B.}:  The array {\tt reqHandle} represents allocated memory that 
! must be deallocated when it is no longer needed.  Failure to do so will 
! create a memory leak.
!
! !INTERFACE:

 subroutine isend_(outgoingGSMap, destID, TagBase, comm, reqHandle, status)

!
! !USES:
!
      use m_mpif90
      use m_die, only : MP_perr_die,die
      use m_stdio

      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_ngseg => ngseg
      use m_GlobalSegMap, only : GlobalSegMap_comp_id => comp_ID
      use m_GlobalSegMap, only : GlobalSegMap_gsize => gsize

      implicit none

! !INPUT PARAMETERS:

  type(GlobalSegMap),    intent(IN)  :: outgoingGSMap
  integer,               intent(IN)  :: destID
  integer,               intent(IN)  :: TagBase
  integer,               intent(IN)  :: comm

! !OUTPUT PARAMETERS: 

  integer, dimension(:), pointer     :: reqHandle
  integer, optional,     intent(OUT) :: status

! !REVISION HISTORY:
! 13Aug03 - J.W. Larson <larson@mcs.anl.gov> - API and initial version.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::isend_'

  integer :: ierr
  integer :: SendBuffer(3)

  if(present(status)) status = 0 ! the success value

  allocate(reqHandle(4), stat=ierr)
  if(ierr /= 0) then
     write(stderr,'(2a,i8)') myname_, &
	  'FATAL--allocation of send buffer failed with ierr=',ierr
     call die(myname_)
  endif

       ! First order of business--query outgoingGSMap for
       ! number of segments, component ID number, and grid
       ! size, packing them in SendBuffer(:).

  SendBuffer(1) = GlobalSegMap_comp_id(outgoingGSMap)
  SendBuffer(2) = GlobalSegMap_ngseg(outgoingGSMap)
  SendBuffer(3) = GlobalSegMap_gsize(outgoingGSMap)

       ! Next, send the buffer size to destID so it can prepare a
       ! receive buffer of the correct size.

  call MPI_ISEND(SendBuffer, 3, MP_Type(SendBuffer(1)), destID, TagBase, &
                comm, reqHandle(1), ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_, 'Send buffer size failed',ierr)
  endif

       ! Send segment information data (3 messages)

  call MPI_ISEND(outgoingGSMap%start, SendBuffer(2), &
                MP_Type(outgoingGSMap%start(1)), &
                destID, TagBase+1, comm, reqHandle(2), ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_, 'Send outgoingGSMap%start failed',ierr)
  endif

  call MPI_ISEND(outgoingGSMap%length, SendBuffer(2), &
                MP_Type(outgoingGSMap%length(1)), &
                destID, TagBase+2, comm, reqHandle(3), ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_, 'Send outgoingGSMap%length failed',ierr)
  endif

  call MPI_ISEND(outgoingGSMap%pe_loc, SendBuffer(2), &
                MP_Type(outgoingGSMap%pe_loc(1)), &
                destID, TagBase+3, comm, reqHandle(4), ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_, 'Send outgoingGSMap%pe_loc failed',ierr)
  endif

 end subroutine isend_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: irecv_ - Point-to-point Non-blocking Receive of a GlobalSegMap
!
! !DESCRIPTION:
! This routine performs a non-blocking receive of a {\tt GlobalSegMap} (the 
! input argument {\tt outgoingGSMap}) from the processor whose rank is 
! {\tt sourceID} on the MPI communicator associated with the integer 
! handle {\tt comm}.  The input {\tt INTEGER} argument {\tt TagBase} 
! is used to generate tags for the messages associated with this operation; 
! there are four messages involved, so the user should avoid using tag 
! values {\tt TagBase} and {\tt TagBase + 3}.  Three of the four messages 
! are non-blocking (the first message defines receive buffer sizes and thus
! must block), and the request handles for them are returned in the output
! {\tt INTEGER} array {\tt reqHandle}, which can be checked for completion 
! using any of MPI's wait functions.  The success (failure) of this operation 
! is reported in the zero (non-zero) value of the optional {\tt INTEGER} 
! output variable {\tt status}.
!
! {\bf N.B.}:  The array {\tt reqHandle} represents allocated memory that 
! must be deallocated when it is no longer needed.  Failure to do so will 
! create a memory leak.
!
! !INTERFACE:

 subroutine irecv_(incomingGSMap, sourceID, TagBase, comm, reqHandle, status)

!
! !USES:
!
      use m_mpif90
      use m_die, only : MP_perr_die, die
      use m_stdio

      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_init => init

      implicit none

! !INPUT PARAMETERS:

  integer,               intent(IN)  :: sourceID
  integer,               intent(IN)  :: TagBase
  integer,               intent(IN)  :: comm

! !OUTPUT PARAMETERS: 

  type(GlobalSegMap),    intent(OUT) :: incomingGSMap
  integer, dimension(:), pointer     :: reqHandle
  integer, optional,     intent(OUT) :: status

! !REVISION HISTORY:
! 13Aug03 - J.W. Larson <larson@mcs.anl.gov> - API and initial version.

!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::irecv_'

  integer :: ierr
  integer :: MPstatus(MP_STATUS_SIZE)
  integer :: RecvBuffer(3)

  if(present(status)) status = 0 ! the success value

       ! Receive the GlobalSegMap's basic constants:  component id,
       ! grid size, and number of segments.  The number of segments
       ! is needed to construct the arrays into which segment 
       ! information will be received.  Thus, this receive blocks.

  call MPI_RECV(RecvBuffer, 3, MP_Type(RecvBuffer(1)), sourceID, &
                TagBase, comm, MPstatus, ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_, 'Receive of RecvBuffer failed',ierr)
  endif

       ! Create Empty GlobaSegMap into which segment information
       ! will be received

  call GlobalSegMap_init(incomingGSMap, RecvBuffer(1), RecvBuffer(2), &
                         RecvBuffer(3))

       ! Set up request handle array

  allocate(reqHandle(3), stat=ierr)
  if(ierr /= 0) then
     write(stderr,'(2a,i8)') myname_, &
	  'FATAL--allocation of reqHandle(:) failed with ierr=',ierr
     call die(myname_)
  endif

       ! Receive segment information data (3 messages)

  call MPI_IRECV(incomingGSMap%start, RecvBuffer(2), &
                MP_Type(incomingGSMap%start(1)), &
                sourceID, TagBase+1, comm, reqHandle(1), ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_, 'Send incomingGSMap%start failed',ierr)
  endif

  call MPI_IRECV(incomingGSMap%length, RecvBuffer(2), &
                MP_Type(incomingGSMap%length(1)), &
                sourceID, TagBase+2, comm, reqHandle(2), ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_, 'Send incomingGSMap%length failed',ierr)
  endif

  call MPI_IRECV(incomingGSMap%pe_loc, RecvBuffer(2), &
                MP_Type(incomingGSMap%pe_loc(1)), &
                sourceID, TagBase+3, comm, reqHandle(3), ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_, 'Send incomingGSMap%pe_loc failed',ierr)
  endif

 end subroutine irecv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: bcast_ - broadcast a GlobalSegMap object
!
! !DESCRIPTION:
!
! The routine {\tt bcast\_()} takes the input/output {\em GlobalSegMap} 
! argument {\tt GSMap} (on input valid only on the {\tt root} process,
! on output valid on all processes) and broadcasts it to all processes
! on the communicator associated with the F90 handle {\tt comm}.  The
! success (failure) of this operation is returned as a zero (non-zero) 
! value of the optional output {\tt INTEGER} argument {\tt status}.
!
! !INTERFACE:

 subroutine bcast_(GSMap, root, comm, status)

!
! !USES:
!
      use m_mpif90
      use m_die, only : MP_perr_die,die
      use m_stdio

      use m_GlobalSegMap, only : GlobalSegMap

      implicit none

! !INPUT PARAMETERS:

      integer,            intent(in)     :: root
      integer,            intent(in)     :: comm

! !INPUT/OUTPUT PARAMETERS: 

      type(GlobalSegMap), intent(inout)  :: GSMap  ! Output GlobalSegMap

! !OUTPUT PARAMETERS: 

      integer, optional,  intent(out)    :: status ! global vector size

! !REVISION HISTORY:
! 17Oct01 - J.W. Larson <larson@mcs.anl.gov> - Initial version.
! 11Aug03 - J.W. Larson <larson@mcs.anl.gov> - Relocated from original
!           location in m_GlobalSegMap.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::bcast_'

  integer :: myID, ierr, n
  integer, dimension(:), allocatable :: IntBuffer

       ! Step One:  which process am I?

  call MP_COMM_RANK(comm, myID, ierr)
  if(ierr /= 0) call MP_perr_die(myname_,'MP_comm_rank()',ierr)

       ! Step Two:  Broadcast the scalar bits of the GlobalSegMap from
       ! the root.

  allocate(IntBuffer(3), stat=ierr) ! allocate buffer space (all PEs)
  if(ierr /= 0) then
    if(.not. present(status)) then
       call die(myname_,'allocate(IntBuffer)',ierr)
    else
       write(stderr,*) myname_,':: error during allocate(IntBuffer)'
       status = 2
       return
    endif
  endif

  if(myID == root) then ! pack the buffer
     IntBuffer(1) = GSMap%comp_id
     IntBuffer(2) = GSMap%ngseg
     IntBuffer(3) = GSMap%gsize
  endif

  call MPI_BCAST(IntBuffer, 3, MP_type(IntBuffer(1)), root, comm, ierr)
  if(ierr /= 0) call MP_perr_die(myname_,'MPI_BCAST(IntBuffer)',ierr)

  if(myID /= root) then ! unpack from buffer to GSMap
     GSMap%comp_id = IntBuffer(1)
     GSMap%ngseg = IntBuffer(2)
     GSMap%gsize = IntBuffer(3)
  endif

  deallocate(IntBuffer, stat=ierr) ! deallocate buffer space
  if(ierr /= 0) then
    if(.not. present(status)) then
       call die(myname_,'deallocate(IntBuffer)',ierr)
    else
       write(stderr,*) myname_,':: error during deallocate(IntBuffer)'
       status = 4
       return
    endif
  endif

       ! Step Three:  Broadcast the vector bits of GSMap from the root.
       ! Pack them into one big array to save latency costs associated
       ! with multiple broadcasts.

  allocate(IntBuffer(3*GSMap%ngseg), stat=ierr) ! allocate buffer space (all PEs)
  if(ierr /= 0) then
    if(.not. present(status)) then
       call die(myname_,'second allocate(IntBuffer)',ierr)
    else
       write(stderr,*) myname_,':: error during second allocate(IntBuffer)'
       status = 5
       return
    endif
  endif 

  if(myID == root) then ! pack outgoing broadcast buffer
     do n=1,GSMap%ngseg
	IntBuffer(n) = GSMap%start(n)
	IntBuffer(GSMap%ngseg+n) = GSMap%length(n)
	IntBuffer(2*GSMap%ngseg+n) = GSMap%pe_loc(n)
     end do
  endif

  call MPI_BCAST(IntBuffer, 3*GSMap%ngseg, MP_Type(IntBuffer(1)), root, comm, ierr)
  if(ierr /= 0) call MP_perr_die(myname_,'Error in second MPI_BCAST(IntBuffer)',ierr)

  if(myID /= root) then ! Allocate GSMap%start, GSMap%length,...and fill them

     allocate(GSMap%start(GSMap%ngseg), GSMap%length(GSMap%ngseg), &
	      GSMap%pe_loc(GSMap%ngseg), stat=ierr)
     if(ierr /= 0) then
	if(.not. present(status)) then
	   call die(myname_,'off-root allocate(GSMap%start...)',ierr)
	else
	   write(stderr,*) myname_,':: error during off-root allocate(GSMap%start...)'
	   status = 7
	   return
	endif
     endif

     do n=1,GSMap%ngseg ! unpack the buffer into the GlobalSegMap
	GSMap%start(n) = IntBuffer(n)
	GSMap%length(n) = IntBuffer(GSMap%ngseg+n)
	GSMap%pe_loc(n) = IntBuffer(2*GSMap%ngseg+n)
     end do

  endif

       ! Clean up buffer space:

  deallocate(IntBuffer, stat=ierr) 
  if(ierr /= 0) then
    if(.not. present(status)) then
       call die(myname_,'second deallocate(IntBuffer)',ierr)
    else
       write(stderr,*) myname_,':: error during second deallocate(IntBuffer)'
       status = 8
       return
    endif
  endif 

 end subroutine bcast_

 end module m_GlobalSegMapComms
