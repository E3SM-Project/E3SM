!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_LocalSegMap - a nontrivial 1-D decomposition of an array.
!
! !DESCRIPTION:
! Consider the problem of the 1-dimensional decomposition of an array 
! across multiple processes.  If each process owns only one contiguous 
! segment, then the {\tt GlobalMap} (see {\tt m\_GlobalMap} for details) 
! is sufficient to describe the decomposition.  If, however, each 
! process owns multiple, non-adjacent segments of the array, a more 
! sophisticated approach is needed.   The {\tt LocalSegMap} data type 
! allows one to describe a one-dimensional decomposition of an array
! with each process owning multiple, non-adjacent segments of the array.
!
! !INTERFACE:

 module m_LocalSegMap

      implicit none

      private	! except

      public :: LocalSegMap		! The class data structure
      public :: lsize
      public :: init
      public :: clean
      public :: local_segs

    type LocalSegMap
      integer :: nseg				! No. of Segments
      integer,dimension(:),pointer :: start	! Segment start indices
      integer,dimension(:),pointer :: length	! Segment lengths
      integer,dimension(:),pointer :: gstart	! Global start indices
    end type LocalSegMap

    interface gsize; module procedure gsize_; end interface
    interface init ; module procedure	&
	initd_,	&	! initialize from all PEs
	initr_		! initialize from the root
    end interface
    interface clean; module procedure clean_; end interface
    interface rank ; module procedure rank_ ; end interface

! !REVISION HISTORY:
! 	28Sep98 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_LocalSegMap'

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

    subroutine initd_(LSMap,ln,comm)
!
! !USES:
!
      use m_mpif90
      use m_die

      implicit none

      type(LocalSegMap),intent(out) :: LSMap
      integer,intent(in) :: ln	! the local size
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	29Sep98 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initd_'
  integer :: nPEs,myID,ier,l,i

  call MP_comm_size(comm,nPEs,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_size()',ier)

  call MP_comm_rank(comm,myID,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_rank()',ier)

  allocate(LSMap%counts(0:nPEs-1),LSMap%displs(0:nPEs-1),stat=ier)
  if(ier /= 0) call perr_die(myname_,'allocate()',ier)

#ifdef MALL_ON
	call mall_ci(size(transfer(LSMap%counts,(/1/))),myname_)
	call mall_ci(size(transfer(LSMap%displs,(/1/))),myname_)
#endif

  call MPI_allgather(ln,1,MP_INTEGER,LSMap%counts,1,MP_INTEGER,comm,ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_allgather()',ier)

  l=0
  do i=0,nPEs-1
    LSMap%displs(i)=l
    l=l+LSMap%counts(i)
  end do

  LSMap%lsize=LSMap%counts(myID)	! the local size
  LSMap%lleft=LSMap%displs(myID)
  LSMap%gsize=l	! the Local size

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

 subroutine initr_(LSMap,lns,root,comm)
!
! !USES:
!
      use m_mpif90
      use m_die
      use m_stdio
 
     implicit none

      type(LocalSegMap),intent(out) :: LSMap
      integer,dimension(:),intent(in) :: lns	! the distributed sizes
      integer,intent(in) :: root
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	29Sep98 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initr_'
  integer :: nPEs,myID,ier,l,i

  call MP_comm_size(comm,nPEs,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_size()',ier)

  call MP_comm_rank(comm,myID,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_rank()',ier)

  allocate(LSMap%counts(0:nPEs-1),LSMap%displs(0:nPEs-1),stat=ier)
  if(ier /= 0) call perr_die(myname_,'allocate()',ier)

#ifdef MALL_ON
	call mall_ci(size(transfer(LSMap%counts,(/1/))),myname_)
	call mall_ci(size(transfer(LSMap%displs,(/1/))),myname_)
#endif

  if(myID == root) then
    if(size(lns(:)) /= nPEs) then
      write(stderr,'(2a,2(a,i4))') myname_,	&
	': _root_ argument error',		&
	', size(lns) =',size(lns),		&
	', nPEs =',nPEs
      call die(myname_)
    endif

    LSMap%counts(:)=lns(:)
  endif

  call MPI_bcast(LSMap%counts,nPEs,MP_INTEGER,root,comm,ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_bcast()',ier)

  l=0
  do i=0,nPEs-1
    LSMap%displs(i)=l
    l=l+LSMap%counts(i)
  end do

  LSMap%lsize=LSMap%counts(myID)	! the local size
  LSMap%lleft=LSMap%displs(myID)
  LSMap%gsize=l	! the Local size

 end subroutine initr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean the map
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(LSMap)
!
! !USES:
!
      use m_die

      implicit none
 
      type(LocalSegMap),intent(inout) :: LSMap

! !REVISION HISTORY:
! 	29Sep98 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

#ifdef MALL_ON
	if( .not.associated(LSMap%counts) .or.		&
	    .not.associated(LSMap%displs) )		&
		write(stderr,'(2a)') myname_,	&
		  ': trying to clean uninitialized variable'
	call mall_co(size(transfer(LSMap%counts,(/1/))),myname_)
	call mall_co(size(transfer(LSMap%displs,(/1/))),myname_)
#endif

  deallocate(LSMap%counts,LSMap%displs,stat=ier)
  if(ier /= 0) call perr_die(myname_,'deallocate()',ier)

  LSMap%lsize=0
  LSMap%lleft=0
  LSMap%gsize=0

 end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lsize_ - find the Local size from the map
!
! !DESCRIPTION:
!
! !INTERFACE:

 function gsize_(LSMap)

      implicit none

      type(LocalSegMap),intent(in) :: LSMap
      integer :: lsize_

! !REVISION HISTORY:
! 	29Sep98 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lsize_'

  lsize_=LSMap%gsize???

 end function lsize_

 end module m_LocalSegMap

