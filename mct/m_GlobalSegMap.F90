!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_GlobalSegMap - a nontrivial 1-D decomposition of an array.
!
! !DESCRIPTION:
! Consider the problem of the 1-dimensional decomposition of an array 
! across multiple processes.  If each process owns only one contiguous 
! segment, then the {\tt GlobalMap} (see {\tt m\_GlobalMap} or details) 
! is sufficient to describe the decomposition.  If, however, each  
! process owns multiple, non-adjacent segments of the array, a more 
! sophisticated approach is needed.   The {\tt GlobalSegMap} data type 
! allows one to describe a one-dimensional decomposition of an array
! with each process owning multiple, non-adjacent segments of the array.
!
! In the current implementation of the {\tt GlobalSegMap}, there is no
! santity check to guarantee that 
!$${\tt GlobalSegMap\%gsize} = \sum_{{\tt i}=1}^{\tt ngseg} 
! {\tt GlobalSegMap\%length(i)} . $$
! The reason we have not implemented such a check is to allow the user
! to use the {\tt GlobalSegMap} type to support decompositions of both 
! {\em haloed} and {\em masked} data.
!
! !INTERFACE:

 module m_GlobalSegMap

      implicit none

      private	! except

      public :: GlobalSegMap	! The class data structure
      public :: init            ! Create
      public :: clean           ! Destroy
      public :: comp_id         ! Return component ID number
      public :: gsize           ! Return global vector size (excl. halos)
      public :: lsize           ! Return local storage size (incl. halos)
      public :: ngseg           ! Return global number of segments
      public :: nlseg           ! Return local number of segments
      public :: Sort            ! compute index permutation to re-order
                                ! GlobalSegMap%start, GlobalSegMap%length,
                                ! and GlobalSegMap%pe_loc
      public :: Permute         ! apply index permutation to re-order 
                                ! GlobalSegMap%start, GlobalSegMap%length,
                                ! and GlobalSegMap%pe_loc
      public :: SortPermute     ! compute index permutation and apply it to
                                ! re-order the GlobalSegMap components
                                ! GlobalSegMap%start, GlobalSegMap%length,
                                ! and GlobalSegMap%pe_loc

    type GlobalSegMap
      integer :: comp_id			! Component ID number
      integer :: ngseg				! No. of Global segments
      integer :: gsize				! No. of Global elements
      integer :: lsize				! No. of Local elements
      integer,dimension(:),pointer :: start	! global seg. start index
      integer,dimension(:),pointer :: length	! segment lengths
      integer,dimension(:),pointer :: pe_loc	! PE locations
    end type GlobalSegMap

    interface init ; module procedure	&
	initd_,	&	! initialize from all PEs
	initr_		! initialize from the root
    end interface
    interface clean ; module procedure clean_ ; end interface
    interface comp_id  ; module procedure comp_id_  ; end interface
    interface gsize ; module procedure gsize_ ; end interface
    interface lsize ; module procedure lsize_ ; end interface
    interface ngseg ; module procedure ngseg_ ; end interface
    interface nlseg ; module procedure nlseg_ ; end interface
    interface rank  ; module procedure &
	rank1_ , &	! single rank case
	rankm_	        ! degenerate (multiple) ranks for halo case
    end interface
    interface Sort ; module procedure Sort_ ; end interface
    interface Permute ; module procedure &
	PermuteInPlace_ 
    end interface
    interface SortPermute ; module procedure &
	PermuteInPlace_ 
    end interface


! !REVISION HISTORY:
! 	28Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
! 	26Jan01 - J.W. Larson <larson@mcs.anl.gov> - replaced the component
!                 GlobalSegMap%comm with GlobalSegMap%comp_id.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_GlobalSegMap'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initd_ - define the map from distributed data
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine initd_(GSMap, start, length, root, my_comm, &
                   comp_id, pe_loc, gsize)

!
! !USES:
!
      use m_mpif90
      use m_die
      use m_stdio

      implicit none

      type(GlobalSegMap),intent(out)  :: GSMap   ! Output GlobalSegMap

      integer,dimension(:),intent(in) :: start   ! segment local start index 
      integer,dimension(:),intent(in) :: length  ! the distributed sizes
      integer,intent(in)              :: root    ! root on my_com
      integer,intent(in)              :: my_comm ! local communicatior
      integer,intent(in)              :: comp_id ! component model ID
      integer,dimension(:), pointer, optional :: pe_loc ! process location
      integer,intent(in), optional    :: gsize   ! global vector size
                                                 ! (optional).  It can
                                                 ! be computed by this 
                                                 ! routine if no haloing
                                                 ! is assumed.

! !REVISION HISTORY:
! 	29Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
! 	14Nov00 - J.W. Larson <larson@mcs.anl.gov> - final working version
! 	09Jan01 - J.W. Larson <larson@mcs.anl.gov> - repaired:  a subtle 
!                 bug concerning the usage of the argument pe_loc (result 
!                 was the new pointer variable my_pe_loc); a mistake in 
!                 the tag arguments to MPI_IRECV; a bug in the declaration
!                 of the array status used by MPI_WAITALL.
! 	26Jan01 - J.W. Larson <larson@mcs.anl.gov> - replaced optional 
!                 argument gsm_comm with required argument comp_id.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initd_'
  integer :: nPEs, myID, ier, l, i, ngseg

        ! arrays allocated on the root to which data are gathered
  integer, dimension(:), allocatable :: root_start, root_length, root_pe_loc
        ! arrays allocated on the root to coordinate gathering of 
        ! data and non-blocking receives by the root
  integer, dimension(:), allocatable :: counts, displs, reqs
  integer, dimension(:,:), allocatable :: status
        ! data and non-blocking receives by the root
  integer, dimension(:), pointer :: my_pe_loc

        ! Determine local process ID:

  call MP_COMM_RANK(my_comm, myID, ier)

  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_rank()',ier)


        ! Check consistency of sizes of input arrays:

  if(size(length) /= size(start)) then
     ier = -1
     call MP_perr_die(myname_,'length/start array size mismatch',ier)
  endif

  if(present(pe_loc)) then
     if(size(pe_loc) /= size(start)) then
	ier = -1
	call MP_perr_die(myname_,'pe_loc/start array size mismatch',ier)
     endif
  endif

        ! Store in the variable ngseg the local size 
        ! array start(:)

  ngseg = size(start)

        ! If the argument pe_loc is not present, then we are
        ! initializing the GlobalSegMap on the communicator 
        ! my_comm.  We will need pe_loc to be allocated and
        ! with local size given by the input value of ngseg,
        ! and then initialize it with the local process id myID.

  if(present(pe_loc)) then
     my_pe_loc => pe_loc
  else
     allocate(my_pe_loc(ngseg), stat=ier)
     if(ier /= 0) call MP_perr_die(myname_,'allocate(my_pe_loc)',ier)
     my_pe_loc = myID
  endif

  call MPI_COMM_SIZE(my_comm, npes, ier)
  if(ier /= 0) call MP_perr_die(myname_,'MPI_COMM_SIZE()',ier)

        ! Allocate an array of displacements (displs) and counts
        ! to hold the local values of ngseg on the root

  if(myID == root) then
     allocate(counts(0:npes-1), displs(0:npes-1), reqs(0:npes-1), &
	      status(MP_STATUS_SIZE,0:npes-1), stat=ier)
     if (ier /= 0) then  
	call MP_perr_die(myname_, 'allocate(counts,...',ier)
     endif
  endif

        ! Send local number of segments to the root.
        ! Here, the root posts non-blocking receives,
        ! and a barrier brings all the processes together
        ! once the root has good data.

  if(myID == root) then
     do i=0,npes-1
	call MPI_IRECV(counts(i), 1, MP_INTEGER, i, i, &
	     my_comm, reqs(i), ier)
     end do
  endif

  call MPI_SEND(ngseg, 1, MP_INTEGER, root, myID, my_comm, ier)
  if(ier /= 0) call MP_perr_die(myname_,'MPI_COMM_SIZE()',ier)

  if(myID == root) then
     call MPI_WAITALL(size(reqs), reqs, status, ier)
     if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL()',ier)
  endif

  call MPI_BARRIER(my_comm, ier)
  if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL()',ier)

        ! On the root compute the value of ngseg, along with
        ! the entries of counts and displs.

  if(myID == root) then
     ngseg = 0
     do i=0,npes-1
	ngseg = ngseg + counts(i)
        if(i == 0) then
	   displs(i) = 1
	else
	   displs(i) = displs(i-1) + counts(i-1)
	endif
     end do
  endif

        ! Now only the root has the correct value of ngseg.

        ! On the root, allocate memory for the arrays root_start,
        ! and root_length.  If the argument pe_loc is present,
        ! allocate root_pe_loc, too.
  
  if(myID == root) then

     allocate(root_start(ngseg), root_length(ngseg), &
	      root_pe_loc(ngseg), stat=ier)
     if (ier /= 0) then
	call MP_perr_die(myname_, 'allocate(root_start...',ier)
     endif

  endif

        ! Now, each process sends its values of start(:) to fill in 
        ! the appropriate portion of root_start(:y) on the root--post
        ! non-blocking receives on the root first, then the individual
        ! sends, followed by a call to MPI_WAITALL().
   

  if(myID == root) then
     do i=0,npes-1
	call MPI_IRECV(root_start(displs(i)), counts(i), MP_INTEGER, &
	     i, i, my_comm, reqs(i), ier)
     end do
  endif

  call MPI_SEND(start, size(start), MP_INTEGER, root, myID, my_comm, ier)
  if(ier /= 0) call MP_perr_die(myname_,'MPI_COMM_SIZE()',ier)

  if(myID == root) then
     call MPI_WAITALL(size(reqs), reqs, status, ier)
     if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL()',ier)
  endif

  call MPI_BARRIER(my_comm, ier)
  if(ier /= 0) call MP_perr_die(myname_,'MPI_BARRIER()',ier)

        ! Next, each process sends its values of length(:) to fill in 
        ! the appropriate portion of root_length(:) on the root--post
        ! non-blocking receives on the root first, then the individual
        ! sends, followed by a call to MPI_WAITALL().
   

  if(myID == root) then
     do i=0,npes-1
	call MPI_IRECV(root_length(displs(i)), counts(i), MP_INTEGER, &
	     i, i, my_comm, reqs(i), ier)
     end do
  endif

  call MPI_SEND(length, size(length), MP_INTEGER, root, myID, my_comm, ier)
  if(ier /= 0) call MP_perr_die(myname_,'MPI_COMM_SIZE()',ier)

  if(myID == root) then
     call MPI_WAITALL(size(reqs), reqs, status, ier)
     if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL()',ier)
  endif

  call MPI_BARRIER(my_comm, ier)
  if(ier /= 0) call MP_perr_die(myname_,'MPI_BARRIER()',ier)

        ! Finally, if the argument pe_loc is present, each process sends 
        ! its values of pe_loc(:) to fill in the appropriate portion of 
        ! root_pe_loc(:) on the root--post non-blocking receives on the 
        ! root first, then the individual sends, followed by a call to 
        ! MPI_WAITALL().  
   
  if(myID == root) then
     do i=0,npes-1
	call MPI_IRECV(root_pe_loc(displs(i)), counts(i), MP_INTEGER, &
	     i, i, my_comm, reqs(i), ier)
     end do
  endif

  call MPI_SEND(my_pe_loc, size(my_pe_loc), MP_INTEGER, root, myID, &
       my_comm, ier)
  if(ier /= 0) call MP_perr_die(myname_,'MPI_COMM_SIZE()',ier)

  if(myID == root) then
     call MPI_WAITALL(size(reqs), reqs, status, ier)
     if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL()',ier)
  endif

  call MPI_BARRIER(my_comm, ier)
  if(ier /= 0) call MP_perr_die(myname_,'MPI_BARRIER()',ier)

        ! Now, we have everything on the root needed to call initr_().

  if(present(gsize)) then
     call initr_(GSMap, ngseg, root_start, root_length, &
                 root_pe_loc, root, my_comm, comp_id, gsize)
  else
     call initr_(GSMap, ngseg, root_start, root_length, &
                 root_pe_loc, root, my_comm, comp_id)
  endif


        ! Clean up the array pe_loc(:) if it was allocated

  if(.not. present(pe_loc)) then
     deallocate(my_pe_loc, stat=ier)
     if(ier /= 0) call MP_perr_die(myname_, 'deallocate(my_pe_loc)', ier)
  endif

        ! Clean up arrays allocated on the root process:

  if(myID == root) then

        ! Clean up the arrays root_start(:), et cetera...

     deallocate(root_start, root_length, root_pe_loc, stat=ier)
     if(ier /= 0) then
	call MP_perr_die(myname_, 'deallocate(root_start,...)', ier)
     endif

        ! Clean up the arrays counts(:) and displs(:)

     deallocate(counts, displs, reqs, status, stat=ier)
     if(ier /= 0) then
	call MP_perr_die(myname_, 'deallocate(counts,...)', ier)
     endif

  endif

 end subroutine initd_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initr_ initialize the map from the root
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine initr_(GSMap, ngseg, start, length, pe_loc, root,  &
                   my_comm, comp_id, gsize)
!
! !Uses:
!
      use m_mpif90
      use m_die
      use m_stdio
 
     implicit none

      type(GlobalSegMap),intent(out)  :: GSMap   ! Output GlobalSegMap
      integer, intent(in)             :: ngseg   ! no. of global segments
      integer,dimension(:),intent(in) :: start   ! segment local start index 
      integer,dimension(:),intent(in) :: length  ! the distributed sizes
      integer,dimension(:),intent(in) :: pe_loc  ! process location
      integer,intent(in)              :: root    ! root on my_com
      integer,intent(in)              :: my_comm ! local communicatior
      integer,intent(in)              :: comp_id ! component id number
      integer,intent(in), optional    :: gsize   ! global vector size
                                                 ! (optional).  It can
                                                 ! be computed by this 
                                                 ! routine if no haloing
                                                 ! is assumed.
! !REVISION HISTORY:
! 	29Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
! 	09Nov00 - J.W. Larson <larson@mcs.anl.gov> - final working version
! 	10Jan01 - J.W. Larson <larson@mcs.anl.gov> - minor bug fix
! 	12Jan01 - J.W. Larson <larson@mcs.anl.gov> - minor bug fix regarding
!                                                    disparities in ngseg on
!                                                    the root and other 
!                                                    processes
! 	26Jan01 - J.W. Larson <larson@mcs.anl.gov> - replaced optional 
!                 argument gsm_comm with required argument comp_id.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initr_'
  integer :: myID,ier,l,i

        ! Determine the local process ID myID:

  call MPI_COMM_RANK(my_comm, myID, ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_COMM_RANK()',ier)

        ! Argument checking:  check to make sure the arrays
        ! start, length, and pe_loc each have ngseg elements.
        ! If not, stop with an error.  This is done on the 
        ! root process since it owns the initialization data.

  if(myID == root) then
    if( size(start(:)) /= ngseg ) then
      write(stderr,'(2a,2(a,i4))') myname_,	&
	': _root_ argument error',		&
	', size(start) =',size(start),	&
	', ngseg =',ngseg
      call die(myname_)
    endif
    if( size(length(:)) /= ngseg ) then
      write(stderr,'(2a,2(a,i4))') myname_,	&
	': _root_ argument error',		&
	', size(length) =',size(length),	&
	', ngseg =',ngseg
      call die(myname_)
    endif
    if( size(pe_loc(:)) /= ngseg ) then
      write(stderr,'(2a,2(a,i4))') myname_,	&
	': _root_ argument error',		&
	', size(pe_loc) =',size(pe_loc),	&
	', ngseg =',ngseg
      call die(myname_)
    endif
  endif

        ! Initialize GSMap%ngseg and GSMap%comp_id on the root:

  if(myID == root) then
     GSMap%ngseg = ngseg
     GSMap%comp_id = comp_id
  endif

        ! Broadcast the value of GSMap%ngseg

  call MPI_BCAST(GSMap%ngseg, 1, MP_INTEGER, root, my_comm, ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_BCAST(GSmap%ngseg)',ier)

        ! Broadcast the value of GSMap%comp_id

  call MPI_BCAST(GSMap%comp_id, 1, MP_INTEGER, root, my_comm, ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_BCAST(GSmap%comp_id)',ier)

        ! Allocate the components GSMap%start(:), GSMap%length(:),
        ! and GSMap%pe_loc(:)

  allocate(GSMap%start(GSMap%ngseg), GSMap%length(GSMap%ngseg), &
           GSMap%pe_loc(GSMap%ngseg), stat = ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(GSmap%start(:),...',ier)

#ifdef MALL_ON
	call mall_ci(size(transfer(GSMap%start,(/1/))),myname_)
	call mall_ci(size(transfer(GSMap%length,(/1/))),myname_)
	call mall_ci(size(transfer(GSMap%pe_loc,(/1/))),myname_)
#endif

        ! On the root process, initialize GSMap%start(:), GSMap%length(:),
        ! and GSMap%pe_loc(:) with the data contained in start(:), 
        ! length(:) and pe_loc(:), respectively

  if(myID == root) then
     GSMap%start(1:GSMap%ngseg) = start(1:GSMap%ngseg)
     GSMap%length(1:GSMap%ngseg) = length(1:GSMap%ngseg)
     GSMap%pe_loc(1:GSMap%ngseg) = pe_loc(1:GSMap%ngseg)
  endif

        ! Broadcast the root values of GSMap%start(:), GSMap%length(:),
        ! and GSMap%pe_loc(:)

  call MPI_BCAST(GSMap%start, GSMap%ngseg, MP_INTEGER, root, my_comm, ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_BCAST(GSMap%start)',ier)

  call MPI_BCAST(GSMap%length, GSMap%ngseg, MP_INTEGER, root, my_comm, ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_BCAST(GSMap%length)',ier)

  call MPI_BCAST(GSMap%pe_loc, GSMap%ngseg, MP_INTEGER, root, my_comm, ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_BCAST(GSMap%pe_loc)',ier)

        ! If the argument gsize is present, use the root value to
        ! set GSMap%gsize and broadcast it.  If it is not present,
        ! this will be computed by summing the entries of GSM%length(:).
        ! Again, note that if one is storing halo points, the sum will
        ! produce a result larger than the actual global vector.  If
        ! halo points are to be used in the mapping we advise strongly
        ! that the user specify the value gsize as an argument.

  if(present(gsize)) then
     if(myID == root) then
	GSMap%gsize = gsize
     endif
     call MPI_BCAST(GSMap%gsize, 1, MP_INTEGER, root, my_comm, ier)
     if(ier/=0) call MP_perr_die(myname_, 'MPI_BCAST(GSMap%gsize)', ier)
  else
     GSMap%gsize = 0
     do i=1,GSMap%ngseg
	GSMap%gsize = GSMap%gsize + GSMap%length(i)
     end do
  endif

        ! Compute the local size of the distributed vector by summing
        ! the entries of GSMap%length(:) whose corresponding values in
        ! GSMap%pe_loc(:) equal the local process ID.  This automatically
        ! takes into account haloing (if present).

  GSMap%lsize = 0

  do i=1,GSMap%ngseg
     if(GSMap%pe_loc(i) == myID) then
	GSMap%lsize = GSMap%lsize + GSMap%length(i)
     endif
  end do

 end subroutine initr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean the map
!
! !DESCRIPTION:
! This routine deallocates the array components of the {\tt GlobalSegMap}
! argument {\tt GSMap}: {\tt GSMap\%start}, {\tt GSMap\%length}, and
! {\tt GSMap\%pe\_loc}.  It also zeroes out the values of the integer
! components {\tt GSMap\%ngseg}, {\tt GSMap\%comp_id}, {\tt GSMap\%gsize},
! and {\tt GSMap\%lsize}.
!
! !INTERFACE:

    subroutine clean_(GSMap)
!
! !USES:
!
      use m_die

      implicit none
 
      type(GlobalSegMap),intent(inout) :: GSMap

! !REVISION HISTORY:
! 	29Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

#ifdef MALL_ON
	if( .not.associated(GSMap%start)  .or.		&
	    .not.associated(GSMap%length) .or.		&
	    .not.associated(GSMap%pe_loc) )		&
		write(stderr,'(2a)') myname_,	&
		  ': trying to clean uninitialized variable'
	call mall_co(size(transfer(GSMap%start,(/1/))),myname_)
	call mall_co(size(transfer(GSMap%length,(/1/))),myname_)
	call mall_co(size(transfer(GSMap%pe_loc,(/1/))),myname_)
#endif

  deallocate(GSMap%start, GSMap%length, GSMap%pe_loc, stat=ier)
  if(ier /= 0) call perr_die(myname_,'deallocate(GSMap%start,...)',ier)

  GSMap%ngseg = 0
  GSMap%comp_id  = 0
  GSMap%gsize = 0
  GSMap%lsize = 0

 end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ngseg_ - Return the global number of segments from the map
!
! !DESCRIPTION:
! The function {\tt ngseg\_()} returns the global number of vector
! segments in the {\tt GlobalSegMap} argument {\tt GSMap}.  This is
! merely the value of {\tt GSMap\%ngseg}.
!
! !INTERFACE:

 function ngseg_(GSMap)

      implicit none

      type(GlobalSegMap),intent(in) :: GSMap
      integer :: ngseg_

! !REVISION HISTORY:
! 	29Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ngseg_'

  ngseg_=GSMap%ngseg

 end function ngseg_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: nlseg_ - Return the global number of segments from the map
!
! !DESCRIPTION:
! The function {\tt nlseg\_()} returns the number of vector segments 
! in the {\tt GlobalSegMap} argument {\tt GSMap} that reside on the 
! process specified by the input argument {\tt pID}.  This is the 
! number of entries {\tt GSMap\%pe\_loc} whose value equals {\tt pID}.
!
! !INTERFACE:

 function nlseg_(GSMap, pID)

      implicit none

      type(GlobalSegMap),intent(in) :: GSMap
      integer,           intent(in) :: pID

      integer :: nlseg_

! !REVISION HISTORY:
! 	29Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::nlseg_'
  integer :: i, nlocseg

        ! Initialize the number of segments residing on pID, nlocseg

  nlocseg = 0

        ! Compute the number of segments residing on pID, nlocseg

  do i=0,GSMap%ngseg
     if(GSMap%pe_loc(i) == pID) then
	nlocseg = nlocseg + 1
     endif
  end do

        ! Return the total

  nlseg_ = nlocseg

 end function nlseg_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: comp_id_ - Return the commponent ID from the GlobalSegMap.
!
! !DESCRIPTION:
! The function {\tt comp\_id\_()} returns component ID number stored in
! {\tt GSMap\%comp\_id}.
!
! !INTERFACE:

 function comp_id_(GSMap)

      implicit none

      type(GlobalSegMap),intent(in) :: GSMap
      integer :: comp_id_

! !REVISION HISTORY:
! 	29Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
! 	26Jan01 - J.W. Larson <larson@mcs.anl.gov> - renamed comp_id_
!                 to fit within MCT_World component ID context.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::comp_id_'

  comp_id_ = GSMap%comp_id

 end function comp_id_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gsize_ - Return the global vector size from the GlobalSegMap.
!
! !DESCRIPTION:
! The function {\tt gsize\_()} takes the input {\tt GlobalSegMap} 
! arguement {\tt GSMap} and returns the global vector length stored
! in {\tt GlobalSegMap\%gsize}.
!
! !INTERFACE:

 function gsize_(GSMap)

      implicit none

      type(GlobalSegMap),intent(in) :: GSMap
      integer :: gsize_

! !REVISION HISTORY:
! 	29Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gsize_'

  gsize_=GSMap%gsize

 end function gsize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lsize_ - find the local storage size from the map
!
! !DESCRIPTION:
!
! !INTERFACE:

 function lsize_(GSMap)

      implicit none

      type(GlobalSegMap),intent(in) :: GSMap
      integer :: lsize_

! !REVISION HISTORY:
! 	29Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lsize_'

  lsize_=GSMap%lsize

 end function lsize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rank1_ - rank which process owns a datum with given global 
! index.
!
! !DESCRIPTION:
! This routine assumes that there is one process that owns the datum with
! a given global index.  It should not be used when the input 
! {\tt GlobalSegMap} argument {\tt GSMap} has been built to incorporate
! halo points.
!
! !INTERFACE:

    subroutine rank1_(GSMap, i_g, rank)

      implicit none

      type(GlobalSegMap), intent(in) :: GSMap   ! input GlobalSegMap
      integer,            intent(in) :: i_g	! a global index
      integer,           intent(out) :: rank    ! the pe on which this
                                                ! element resides
! !REVISION HISTORY:
! 	29Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rank1_'
  integer :: i,ilc,ile

        ! Initially, set the rank to -1 (invalid).
  rank=-1	

  do i=1,size(GSMap%start)
    ilc = GSMap%start(i)
    ile = ilc + GSMap%length(i)

		! If i_g in (ilc,ile].  Note that i_g := [1:..]

    if(ilc < i_g .and. i_g <= ile) then
      rank = GSMap%pe_loc(i)
      return
    endif
  end do

 end subroutine rank1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rankm_ - rank which processes own a datum with given global 
! index.
!
! !DESCRIPTION:
! This routine assumes that there may be more than one process that owns 
! the datum with a given global index.  This routine should  be used when 
! the input {\tt GlobalSegMap} argument {\tt GSMap} has been built to 
! incorporate ! halo points.  {\em Nota Bene}:  The output array {\tt rank} 
! is allocated in this routine and must be deallocated by the routine calling 
! {\tt rankm\_()}.  Failure to do so could result in a memory leak.
!
! !INTERFACE:

    subroutine rankm_(GSMap, i_g, num_loc, rank)

      implicit none

      type(GlobalSegMap), intent(in) :: GSMap   ! input GlobalSegMap
      integer,            intent(in) :: i_g	! a global index
      integer,           intent(out) :: num_loc ! the number of processes
                                                ! which own element i_g
      integer, dimension(:), pointer :: rank    ! the process(es) on which 
                                                ! element i_g resides
! !REVISION HISTORY:
! 	29Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rankm_'
  integer :: i, ilc, ile, ier, n

        ! First sweep:  determine the number of processes num_loc 
        ! that own the given datum:

  num_loc = 0

  do i=1,size(GSMap%start)

    ilc = GSMap%start(i)
    ile = ilc + GSMap%length(i) - 1

		! If i_g in [ilc,ile].  Note that i_g := [1:..]

    if(ilc <= i_g .and. i_g <= ile) then
      num_loc = num_loc + 1
    endif

  end do

  if(num_loc == 0) then

        ! If i_g is nowhere to be found in GSMap, set num_loc to
        ! unity and return a null value for rank

     num_loc = 1
     allocate(rank(num_loc), stat=ier)
     rank = -1	! null value
     return

  else
        ! Allocate output array rank(1:num_loc)

     allocate(rank(num_loc), stat=ier)

        ! Second sweep:  fill in the entries to rank(:)

     n = 0 ! counter

     do i=1,size(GSMap%start)

	ilc = GSMap%start(i)
	ile = ilc + GSMap%length(i) - 1

	! If i_g in [ilc,ile].  Note that i_g := [1:..]

	if(ilc <= i_g .and. i_g <= ile) then
	   n = n + 1
	   rank(n) = GSMap%pe_loc(i)
	endif

     end do

  endif

 end subroutine rankm_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: Sort_ - generate index permutation for GlobalSegMap.
!
! !DESCRIPTION:
! {\tt Sort\_()} uses the supplied keys {\tt key1} and {\tt key2} to 
! generate a permutation {\tt perm} that will put the entries of the 
! components {\tt GlobalSegMap\%start}, {\tt GlobalSegMap\%length} and
! {\tt GlobalSegMap\%pe\_loc} in {\em ascending} lexicographic order.
!
! {\bf N.B.:} {\tt Sort\_()} returns an allocated array {\tt perm(:)}.  It
! the user must deallocate this array once it is no longer needed.  Failure
! to do so could create a memory leak.
!
! !INTERFACE:

    subroutine Sort_(GSMap, key1, key2, perm)
!
! !USES:
!
      use m_die ,          only : die
      use m_SortingTools , only : IndexSet
      use m_SortingTools , only : IndexSort

      implicit none

      type(GlobalSegMap), intent(in) :: GSMap   ! input GlobalSegMap
      integer, dimension(:), intent(in)           :: key1 ! first sort key
      integer, dimension(:), intent(in), optional :: key2 ! second sort key
      integer, dimension(:), pointer :: perm    ! output index permutation
! !REVISION HISTORY:
! 	02Feb01 - J.W. Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::Sort_'

  integer :: ierr, length

  length = ngseg_(GSMap)

        ! Argument checking.  are key1 and key2 (if supplied) the
        ! same length as the components of GSMap?  If not, stop with
        ! an error.

  ierr = 0

  if(size(key1) /= length) then
     ierr = 1
     call MP_perr_die(myname_,'key1 GSMap size mismatch',ierr)
  endif

  if(present(key2)) then
     if(size(key2) /= length) then
        ierr = 2
	call MP_perr_die(myname_,'key2 GSMap size mismatch',ierr)
     endif
     if(size(key1) /= size(key2)) then
        ierr = 3
	call MP_perr_die(myname_,'key1 key2 size mismatch',ierr)
     endif
  endif

        ! allocate space for permutation array perm(:)

  allocate(perm(length), stat=ierr)
  if(ierr /= 0) call MP_perr_die(myname_,'allocate(perm)',ierr)

        ! Initialize perm(i)=i, for i=1,length

  call IndexSet(perm)
 
        ! Index permutation is achieved by successive calls to IndexSort(),
        ! with the keys supplied one at a time in the order reversed from
        ! the desired sort order.

  if(present(key2)) then
     call IndexSort(length, perm, key2, descend=.false.)
  endif

  call IndexSort(length, perm, key1, descend=.false.)

        ! Yes, it is that simple.  The desired index permutation is now
        ! stored in perm(:)

 end subroutine Sort_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: PermuteInPlace_ - apply index permutation to GlobalSegMap.
!
! !DESCRIPTION:
! {\tt PermuteInPlace\_()} uses a supplied index permutation {\tt perm} 
! to re-order {\tt GlobalSegMap\%start}, {\tt GlobalSegMap\%length} and
! {\tt GlobalSegMap\%pe\_loc}.
!
! !INTERFACE:

    subroutine PermuteInPlace_(GSMap, perm)
!
! !USES:
!
      use m_die ,          only : die
      use m_SortingTools , only : IndexSet
      use m_SortingTools , only : Permute

      implicit none

      type(GlobalSegMap), intent(inout) :: GSMap
      integer, dimension(:), intent(in) :: perm

! !REVISION HISTORY:
!       02Feb01 - J.W. Larson <larson@mcs.anl.gov> - initial version.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::PermuteInPlace_'

  integer :: length, ierr

  length = ngseg_(GSMap)

        ! Argument checking.  Do the components of GSMap
        ! (e.g. GSMap%start) have the same length as the
        ! permutation array perm?  If not, stop with an error.

  ierr = 0

  if(size(perm) /= length) then
     ierr = 1
     call MP_perr_die(myname_,'perm GSMap size mismatch',ierr)
  endif

        ! In-place index permutation using perm(:) :

  call Permute(GSMap%start,perm,length)
  call Permute(GSMap%length,perm,length)
  call Permute(GSMap%pe_loc,perm,length)

        ! Now, the components of GSMap are ordered according to
        ! perm(:).

 end subroutine PermuteInPlace_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: SortPermuteInPlace_ - Sort in-place GlobalSegMap components.
!
! !DESCRIPTION:
! {\tt SortPermuteInPlace\_()} uses a the supplied key(s) to generate 
! and apply an index permutation that will place the {\tt GlobalSegMap}
! components {\tt GlobalSegMap\%start}, {\tt GlobalSegMap\%length} and
! {\tt GlobalSegMap\%pe\_loc} in lexicographic order.
!
! !INTERFACE:

    subroutine SortPermuteInPlace_(GSMap, key1, key2)
!
! !USES:
!
      use m_die ,          only : die

      implicit none

      type(GlobalSegMap), intent(inout) :: GSMap
      integer, dimension(:), intent(in)           :: key1
      integer, dimension(:), intent(in), optional :: key2
! !REVISION HISTORY:
!       02Feb01 - J.W. Larson <larson@mcs.anl.gov> - initial version.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::SortPermuteInPlace_'

  integer :: length, ierr
  integer, dimension(:), pointer :: perm

   length = ngseg_(GSMap)

        ! Argument checking.  are key1 and key2 (if supplied) the
        ! same length as the components of GSMap?  If not, stop with
        ! an error.
  ierr = 0
  if(size(key1) /= length) then
     ierr = 1
     call MP_perr_die(myname_,'key1 GSMap size mismatch',ierr)
  endif

  if(present(key2)) then
     if(size(key2) /= length) then
        ierr = 2
	call MP_perr_die(myname_,'key2 GSMap size mismatch',ierr)
     endif
     if(size(key1) /= size(key2)) then
        ierr = 3
	call MP_perr_die(myname_,'key1 key2 size mismatch',ierr)
     endif
  endif

        ! Generate desired index permutation:      

  if(present(key2)) then
     call Sort_(GSMap, key1, key2, perm)
  else
     call Sort_(GSMap, key1=key1, perm=perm)
  endif

        ! Apply index permutation:      

  call PermuteInPlace_(GSMap, perm)

        ! Now the components of GSMap have been re-ordered.
        ! Deallocate the index permutation array perm(:)

  deallocate(perm, stat=ierr)
  if(ierr /= 0) call MP_perr_die(myname_,'deallocate(perm...)',ierr)

 end subroutine SortPermuteInPlace_

 end module m_GlobalSegMap

