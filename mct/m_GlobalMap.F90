!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_GlobalMap - a size description of a distributed array
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_GlobalMap
      implicit none
      private	! except

      public :: GlobalMap		! The class data structure
      public :: gsize
      public :: lsize
      public :: init
      public :: init_remote
      public :: clean
      public :: rank
      public :: local_range

    type GlobalMap
      integer :: gsize				! the Global size
      integer :: lsize				! my local size
      integer :: lleft				! my left index
      integer,dimension(:),pointer :: counts	! all local sizes
      integer,dimension(:),pointer :: displs	! PE ordered locations
    end type GlobalMap

    interface gsize; module procedure gsize_; end interface
    interface lsize; module procedure lsize_; end interface
    interface init ; module procedure	&
	initd_,	&	! initialize from all PEs
	initr_		! initialize from the root
    end interface
    interface init_remote; module procedure init_remote_; end interface
    interface clean; module procedure clean_; end interface
    interface rank ; module procedure rank_ ; end interface
    interface local_range; module procedure range_; end interface

! !REVISION HISTORY:
! 	21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	09Nov00 - J.W. Larson <larson@mcs.anl.gov> - added init_remote
!                 interface.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_GlobalMap'


contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initd_ - define the map from distributed data
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initd_(GMap,ln,comm)
      use m_mpif90
      use m_die
      implicit none
      type(GlobalMap),intent(out) :: GMap
      integer,intent(in) :: ln	! the local size
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initd_'
  integer :: nPEs,myID,ier,l,i

  call MP_comm_size(comm,nPEs,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_size()',ier)

  call MP_comm_rank(comm,myID,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_rank()',ier)

  allocate(GMap%counts(0:nPEs-1),GMap%displs(0:nPEs-1),stat=ier)
  if(ier /= 0) call perr_die(myname_,'allocate()',ier)

#ifdef MALL_ON
	call mall_ci(size(transfer(GMap%counts,(/1/))),myname_)
	call mall_ci(size(transfer(GMap%displs,(/1/))),myname_)
#endif

  call MPI_allgather(ln,1,MP_INTEGER,GMap%counts,1,MP_INTEGER,comm,ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_allgather()',ier)

  l=0
  do i=0,nPEs-1
    GMap%displs(i)=l
    l=l+GMap%counts(i)
  end do

  GMap%lsize=GMap%counts(myID)	! the local size
  GMap%lleft=GMap%displs(myID)
  GMap%gsize=l	! the global size

 end subroutine initd_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initr_ initialize the map from the root
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initr_(GMap,lns,root,comm)
      use m_mpif90
      use m_die
      use m_stdio
      implicit none
      type(GlobalMap),intent(out) :: GMap
      integer,dimension(:),intent(in) :: lns	! the distributed sizes
      integer,intent(in) :: root
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	29May98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initr_'
  integer :: nPEs,myID,ier,l,i

  call MP_comm_size(comm,nPEs,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_size()',ier)

  call MP_comm_rank(comm,myID,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_rank()',ier)

  allocate(GMap%counts(0:nPEs-1),GMap%displs(0:nPEs-1),stat=ier)
  if(ier /= 0) call perr_die(myname_,'allocate()',ier)

#ifdef MALL_ON
	call mall_ci(size(transfer(GMap%counts,(/1/))),myname_)
	call mall_ci(size(transfer(GMap%displs,(/1/))),myname_)
#endif

  if(myID == root) then
    if(size(lns(:)) /= nPEs) then
      write(stderr,'(2a,2(a,i4))') myname_,	&
	': _root_ argument error',		&
	', size(lns) =',size(lns),		&
	', nPEs =',nPEs
      call die(myname_)
    endif

    GMap%counts(:)=lns(:)
  endif

  call MPI_bcast(GMap%counts,nPEs,MP_INTEGER,root,comm,ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_bcast()',ier)

  l=0
  do i=0,nPEs-1
    GMap%displs(i)=l
    l=l+GMap%counts(i)
  end do

  GMap%lsize=GMap%counts(myID)	! the local size
  GMap%lleft=GMap%displs(myID)
  GMap%gsize=l	! the global size

 end subroutine initr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_remote_ initialize remote GlobalMap from the root
!
! !DESCRIPTION:
! This method takes data describing a GlobalMap for another communicator
!
! !INTERFACE:

 subroutine init_remote_(GMap, remote_lns, remote_npes, my_root, &
	                    my_comm, remote_comm)
      use m_mpif90
      use m_die
      use m_stdio
      implicit none
      type(GlobalMap),intent(out) :: GMap
      integer, dimension(:) :: remote_lns     ! distributed sizes on 
                                              ! the remote communicator
      integer            :: remote_npes       ! number of processes
                                              ! on remote communicator
      integer,intent(in) :: my_root           ! my root
      integer,intent(in) :: my_comm           ! my communicator
      integer,intent(in) :: remote_comm       ! remote communicator

! !REVISION HISTORY:
! 	08Nov00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_remote_'
  integer :: nPEs,myID,ier,l,i


        ! Which processor am I on communicator my_comm?  Store
        ! the answer in myID:

  call MP_comm_rank(my_comm, myID, ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_rank()',ier)

        ! allocate counts and displacements component arrays
        ! for the sake of compactness, store the value of remote_npes
        ! in the more tersely named variable nPEs.

  nPEs = remote_npes

  allocate(GMap%counts(0:nPEs-1),GMap%displs(0:nPEs-1),stat=ier)
  if(ier /= 0) call perr_die(myname_,'allocate()',ier)

#ifdef MALL_ON
	call mall_ci(size(transfer(GMap%counts,(/1/))),myname_)
	call mall_ci(size(transfer(GMap%displs,(/1/))),myname_)
#endif

        ! On the Root processor, check the size of remote_lns(:)
        ! to see it is equal to nPEs, the number of remote processes,
        ! then store it as GMap%counts and broadcast it.

  if(myID == my_root) then
    if(size(remote_lns(:)) /= nPEs) then
      write(stderr,'(2a,2(a,i4))') myname_,	 &
	': _root_ argument error',		 &
	', size(remote_lns) =',size(remote_lns), &
	', nPEs =',nPEs
      call die(myname_)
    endif

    GMap%counts(:)=remote_lns(:)
  endif

  call MPI_bcast(GMap%counts, nPEs, MP_INTEGER, my_root, my_comm, ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_bcast()',ier)

        ! Now, on each processor of my_comm, compute from 
        ! GMap%counts(:) the entries of GMap%displs(:)

  l=0
  do i=0,nPEs-1
    GMap%displs(i)=l
    l=l+GMap%counts(i)
  end do

  GMap%lsize=GMap%counts(myID)	! the local size
  GMap%lleft=GMap%displs(myID)
  GMap%gsize=l	! the global size

 end subroutine init_remote_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean the map
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(GMap)
      use m_die
      implicit none
      type(GlobalMap),intent(inout) :: GMap

! !REVISION HISTORY:
! 	21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

#ifdef MALL_ON
	if( .not.associated(GMap%counts) .or.		&
	    .not.associated(GMap%displs) )		&
		write(stderr,'(2a)') myname_,	&
		  ': trying to clean uninitialized variable'
	call mall_co(size(transfer(GMap%counts,(/1/))),myname_)
	call mall_co(size(transfer(GMap%displs,(/1/))),myname_)
#endif

  deallocate(GMap%counts,GMap%displs,stat=ier)
  if(ier /= 0) call perr_die(myname_,'deallocate()',ier)

  GMap%lsize=0
  GMap%lleft=0
  GMap%gsize=0

end subroutine clean_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lsize_ - find the local size from the map
!
! !DESCRIPTION:
!
! !INTERFACE:

    function lsize_(GMap)
      implicit none
      type(GlobalMap),intent(in) :: GMap
      integer :: lsize_

! !REVISION HISTORY:
! 	21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lsize_'

  lsize_=GMap%lsize

end function lsize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gsize_ - find the global size from the map
!
! !DESCRIPTION:
!
! !INTERFACE:

    function gsize_(GMap)
      implicit none
      type(GlobalMap),intent(in) :: GMap
      integer :: gsize_

! !REVISION HISTORY:
! 	21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gsize_'

  gsize_=GMap%gsize

end function gsize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rank_ - rank (which PE) of a given global index
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine rank_(GMap,i_g,rank)
      implicit none
      type(GlobalMap),intent(in) :: GMap
      integer, intent(in)  :: i_g	! a global index
      integer, intent(out) :: rank

! !REVISION HISTORY:
! 	05May98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rank_'
  integer :: i,ilc,ile

  rank=-1	! if nowhere fits
  do i=0,size(GMap%displs)-1
    ilc=GMap%displs(i)
    ile=ilc+GMap%counts(i)

		! If i_g in (ilc,ile].  Note that i_g := [1:..]

    if(ilc < i_g .and. i_g <= ile) then
      rank=i
      return
    endif
  end do

end subroutine rank_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: range_ - the range of the local indices
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine range_(GMap,lbnd,ubnd)
      implicit none
      type(GlobalMap),intent(in) :: GMap
      integer,intent(out) :: lbnd
      integer,intent(out) :: ubnd

! !REVISION HISTORY:
! 	05May98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::range_'

  lbnd=GMap%lleft+1
  ubnd=GMap%lleft+GMap%lsize

end subroutine range_

end module m_GlobalMap
