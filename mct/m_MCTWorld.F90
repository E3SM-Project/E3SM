!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_MCTWorld -- MCTWorld Class
!
! !DESCRIPTION:
! Summarize the number of things MCT is hooked up to and how they are
! distributed on processors.  Every component MCT is communicating
! with must call this routine once.
!
! !INTERFACE:

 module m_MCTWorld
!
! !USES:
      use m_List, only : List   ! Support for List components.

      implicit none

      private   ! except

      public :: MCTWorld        ! The class data structure

      public :: init            ! Create a MCTWorld
      public :: clean           ! Destroy a MCTWorld

      public :: NumComponents      ! Number of Components in the MCTWorld
      public :: ComponentNumProcs  ! Number of processes owned by a given
                                   ! component

      public :: ComponentToWorldRank ! Given the rank of a process on a 
                                     ! component, return its rank on the 
                                     ! world communicator
      public :: ComponentRootRank    ! Return the rank on the world 
                                     ! communicator of the root process of 
                                     ! a component

      public :: ThisMCTWorld   ! Instantiation of the MCTWorld

!  
    type MCTWorld
      integer :: myid	       ! my unique id 
      integer :: ncomps	       !  number of components
      integer :: mynprocs      !  number of procs
      integer,dimension(:),pointer :: allids	   ! unique id for each 
                                                   ! component
      integer,dimension(:),pointer :: nprocspid	   ! number of processes 
                                                   ! each component is on
      integer,dimension(:,:),pointer :: idGprocid  ! translate between local 
                                                   ! and global ranks for each 
                                                   ! component
    end type MCTWorld

    type(MCTWorld) :: ThisMCTWorld	! declare an MCTWorld

    interface init  ; module procedure init_  ; end interface
    interface clean ; module procedure clean_ ; end interface
    interface NumComponents ; module procedure &
	 NumComponents_ 
    end interface
    interface ComponentNumProcs ; module procedure &
       ComponentNumProcs_ 
    end interface
    interface ComponentToWorldRank ; module procedure &
       ComponentToWorldRank_ 
    end interface
    interface ComponentRootRank ; module procedure &
       ComponentRootRank_ 
    end interface

! !REVISION HISTORY:
!      19Jan01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
!       5Feb01 - J. Larson <larson@mcs.anl.gov> - added query and
!                local-to-global mapping services NumComponents, 
!                ComponentNumProcs, ComponentToWorldRank, and 
!                ComponentRootRank
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_MCTWorld'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize MCTWorld
!
! !DESCRIPTION:
! Initialize MCTWorld using the total number of components and
! a unique component id (provided by MPH).  Also need the local 
! communicator.
!
! !INTERFACE:

 subroutine init_(ncomps,myid,mycomm)
!
! !USES:
!
      use m_mpif90
      use m_die
      use m_stdio

      implicit none

      integer, intent(in)	       :: ncomps  ! number of components
      integer, intent(in)	       :: myid    ! my component id
      integer, intent(in)	       :: mycomm  ! my communicator

! !REVISION HISTORY:
!       19Jan01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: ier,myGid,myLid,i,mysize,Gsize,j

! arrays allocated on the root to coordinate gathring of data
! and non-blocking receives by the root
  integer, dimension(:), allocatable :: compids,reqs,nprocs,Gprocids
  integer, dimension(:,:),allocatable :: status
  integer, dimension(:,:),pointer :: tmparray
  integer,dimension(:),pointer :: apoint
! ------------------------------------------------------------------

! make sure this has not been called already
  if(associated(ThisMCTWorld%allids) ) then
     write(stderr,'(2a)') myname_, &
      'Trying to initialize MCTWorld twice'
      call die(myname_)
  endif

! determine size on local communicator
  call MP_comm_size(mycomm,mysize,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_size()',ier)

! determine overall size
  call MP_comm_size(MP_COMM_WORLD,Gsize,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_size()',ier)

!  set some parts of the MCTWorld
  ThisMCTWorld%ncomps = ncomps
  ThisMCTWorld%myid = myid
  ThisMCTWorld%mynprocs = mysize

! determine my rank in comm_world
  call MP_comm_rank(MP_COMM_WORLD,myGid,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_rank()',ier)

! determine my rank in local comm
  call MP_comm_rank(mycomm,myLid,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_rank()',ier)


! allocate space on global root
  if(myGid == 0) then
     allocate(nprocs(ncomps),compids(ncomps),&
     reqs(ncomps),status(MP_STATUS_SIZE,ncomps),stat=ier)
     if (ier /= 0) then
        call MP_perr_die(myname_, 'allocate(compids,...)',ier)
     endif
  endif

!!!!!!!!!!!!!!!!!!
!  Gather the component id from the root of each component
!!!!!!!!!!!!!!!!!!
!
!  First on the global root, post a receive for each component
  if(myGid == 0) then
    do i=1,ncomps
       call MPI_IRECV(compids(i), 1, MP_INTEGER, MP_ANY_SOURCE,i, &
	 MP_COMM_WORLD, reqs(i), ier)
       if(ier /= 0) call MP_perr_die(myname_,'MPI_IRECV()',ier)
    enddo
  endif

!  The root on each component sends
  if(myLid == 0) then
    call MPI_SEND(myid,1,MP_INTEGER,0,myid,MP_COMM_WORLD,ier)
    if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL()',ier)
  endif

!  Global root waits for all sends
  if(myGid == 0) then
    call MPI_WAITALL(size(reqs), reqs, status, ier)
    if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL()',ier)
  endif

!  declare space for the component ids
  allocate(ThisMCTWorld%allids(ncomps),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(MCTWorld%allids(:),...',ier)

!  Now store the component ids in the World description and Broadcast
  if(myGid == 0) then
    ThisMCTWorld%allids(1:ncomps) = compids(1:ncomps)
  endif
  call MPI_BCAST(ThisMCTWorld%allids, ncomps, MP_INTEGER, 0, MP_COMM_WORLD, ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_BCast allids',ier)
!!!!!!!!!!!!!!!!!!
! end of component id 
!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!
!  Gather the number of procs from the root of each component
!!!!!!!!!!!!!!!!!!
!
!  First on the global root, post a receive for each component
  if(myGid == 0) then
    do i=1,ncomps
       call MPI_IRECV(nprocs(i), 1, MP_INTEGER, MP_ANY_SOURCE,i, &
	 MP_COMM_WORLD, reqs(i), ier)
       if(ier /= 0) call MP_perr_die(myname_,'MPI_IRECV()',ier)
    enddo
  endif

!  The local root on each component sends
  if(myLid == 0) then
    call MPI_SEND(mysize,1,MP_INTEGER,0,myid,MP_COMM_WORLD,ier)
    if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL()',ier)
  endif

!  Global root waits for all sends
  if(myGid == 0) then
    call MPI_WAITALL(size(reqs), reqs, status, ier)
    if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL()',ier)
  endif

!  Allocate space for the nprocs info
  allocate(ThisMCTWorld%nprocspid(ncomps),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(MCTWorld%nprocspid(:),...',ier)
!  Now store the component ids in the World description and Broadcast
  if(myGid == 0) then
    ThisMCTWorld%nprocspid(1:ncomps) = nprocs(1:ncomps)
  endif
  call MPI_BCAST(ThisMCTWorld%nprocspid, ncomps, MP_INTEGER, 0, MP_COMM_WORLD, ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_BCast nprocspid',ier)

!!!!!!!!!!!!!!!!!!
! end of nprocs
!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! make the master list of global proc ids
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! allocate space on local root to hold global ids
  if(myLid == 0) then
    allocate(Gprocids(mysize),stat=ier)
    if(ier/=0) call MP_perr_die(myname_,'allocate(Gprocids)',ier)
  endif

! gather over the LOCAL comm
  call MPI_GATHER(myGid,1,MP_INTEGER,Gprocids,1,MP_INTEGER,0,mycomm,ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_GATHER Gprocids',ier)

! Now get ready to gather from each local root
! allocate space in MCTWorld
  allocate(ThisMCTWorld%idGprocid(ncomps,0:Gsize-1),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(ThisMCTWorld%idGprocid)',ier)

! allocate a tmp array for the receive
  allocate(tmparray(0:Gsize-1,ncomps),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(tmparray)',ier)

! fill tmparray with a bad rank value for later error checking
  tmparray = -1

!!!!!!!!!!!!!!!!!!
!  Gather the Gprocids from each local root
!!!!!!!!!!!!!!!!!!
!
!  First on the global root, post a receive for each component
  if(myGid == 0) then
    do i=1,ncomps
       apoint => tmparray(0:Gsize-1,i)
       call MPI_IRECV(apoint, ThisMCTWorld%nprocspid(i),MP_INTEGER, &
       MP_ANY_SOURCE,i,MP_COMM_WORLD, reqs(i), ier)
       if(ier /= 0) call MP_perr_die(myname_,'MPI_IRECV()',ier)
    enddo
  endif

!  The root on each component sends
  if(myLid == 0) then
    call MPI_SEND(Gprocids,mysize,MP_INTEGER,0,myid,MP_COMM_WORLD,ier)
    if(ier /= 0) call MP_perr_die(myname_,'MPI_SEND(Gprocids)',ier)
  endif

!  Global root waits for all sends
  if(myGid == 0) then
    call MPI_WAITALL(size(reqs), reqs, status, ier)
    if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL(Gprocids)',ier)
  endif

!  Now store the Gprocids in the World description and Broadcast

  if(myGid == 0) then
    ThisMCTWorld%idGprocid = transpose(tmparray)
  endif

  call MPI_BCAST(ThisMCTWorld%idGprocid, ncomps*Gsize,MP_INTEGER, 0, MP_COMM_WORLD, ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_BCast Gprocids',ier)
!!!!!!!!!!!!!!!!!!
! end of Gprocids
!!!!!!!!!!!!!!!!!!

! if(myGid==17) then
!      do i=1,ThisMCTWorld%ncomps
!       do j=1,ThisMCTWorld%nprocspid(i)
!     write(*,*)'MCTK',myGid,i,j-1,ThisMCTWorld%idGprocid(i,j-1)
!    enddo
!   enddo
! endif

 if(myGid == 0) then
  deallocate(compids,reqs,status,nprocs,tmparray,stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'deallocate(compids,..)',ier)
 endif
 if(myLid == 0) then
  deallocate(Gprocids,stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'deallocate(Gprocids)',ier)
 endif

 end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - Destroy a MCTWorld
!
! !DESCRIPTION:
! This routine deallocates the array componets of the {\tt ThisMCTWorld}
! It also zeros out the integer components.
!
! !INTERFACE:

    subroutine clean_()
!
! !USES:
!
      use m_die

      implicit none

! !REVISION HISTORY:
!       19Jan01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  deallocate(ThisMCTWorld%allids,ThisMCTWorld%nprocspid,ThisMCTWorld%idGprocid,stat=ier)
  if(ier /= 0) call perr_die(myname_,'deallocate(MCTW,...)',ier)

  ThisMCTWorld%myid = 0
  ThisMCTWorld%ncomps = 0
  ThisMCTWorld%mynprocs = 0

 end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: NumComponents_ - Determine number of components in World.
!
! !DESCRIPTION:
! The function {\tt NumComponents\_} takes an input {\tt MCTWorld} 
! argument {\tt World}, and returns the number of component models 
! present.
!
! !INTERFACE:

 integer function NumComponents_(World)
!
! !USES:
!
      use m_die
      use m_stdio

      implicit none

      type(MCTWorld), intent(in)      :: World

! !REVISION HISTORY:
!       05Feb01 - J. Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::NumComponents_'

  integer :: ncomps

  ncomps = World%ncomps

  if(ncomps <= 0) then
     write(stderr,'(2a,1i)') myname,":: invalid no. of components = ",ncomps
     call MP_perr_die(myname_,'ncomps = ',ncomps)
  endif

  NumComponents_ = ncomps

 end function NumComponents_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ComponentNumProcs_ - Number of processes a component owns.
!
! !DESCRIPTION:
! The function {\tt ComponentNumProcs\_} takes an input {\tt MCTWorld} 
! argument {\tt World}, and a component ID {\tt comp\_id}, and returns 
! the number of processes owned by that component.
!
! !INTERFACE:

 integer function ComponentNumProcs_(World, comp_id)
!
! !USES:
!
      use m_die
      use m_stdio

      implicit none

      type(MCTWorld), intent(in)      :: World
      integer,        intent(in)      :: comp_id

! !REVISION HISTORY:
!       05Feb01 - J. Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::ComponentNumPros_'

  integer :: mynprocs

  mynprocs = World%mynprocs

  if(mynprocs <= 0) then
     write(stderr,'(2a,1i)') myname,":: invalid no. of processes = ",mynprocs
     call MP_perr_die(myname_,'Number of processes = ',mynprocs)
  endif

  ComponentNumProcs_ = mynprocs

 end function ComponentNumProcs_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ComponentToWorldRank_ - Determine rank on COMM_WORLD.
!
! !DESCRIPTION:
! The function {\tt ComponentToWorldRank\_} takes an input component ID 
! {\tt comp\_id} and input rank on that component communicator 
! {\tt comp\_rank}, and returns the rank of that process on the world 
! communicator of {\tt MCTWorld}.  The optional {\tt LOGICAL} argument 
! {\tt check} is used to control exhaustive argument validity checking 
! against {\tt World}.  If the argument {\tt check} is not provided, 
! no checking, which can hamper performance will be done.
!
! !INTERFACE:

 integer function ComponentToWorldRank_(comp_rank, comp_id, World, check)
!
! !USES:
!
      use m_die
      use m_stdio

      implicit none

      integer, intent(in)	     :: comp_rank ! process rank on
                                                  ! the communicator
                                                  ! associated with
                                                  ! comp_id
      integer, intent(in)	     :: comp_id   ! component id
      type(MCTWorld), intent(in)     :: World     ! World
      logical, intent(in), optional  :: check     ! exhaustive checking
                                                  ! flag

! !REVISION HISTORY:
!       05Feb01 - J. Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::ComponentToWorldRank_'

  logical :: valid
  integer :: n, world_rank


      ! Do we want the potentially time-consuming argument checks?
      ! The first time we use this function during execution on a
      ! given set of components and component ranks, we will.  In
      ! later invocations, these argument checks are probably not
      ! necessary (unless one alters MCTWorld), and impose a cost
      ! one may wish to avoid.

  if(present(check)) then
     if(check) then

      ! Check argument comp_id for validity--assume initially it is not...

	valid = .false.
	n = 0

	CHECK_COMP_ID_LOOP: do 

	   n = n + 1
	   if(n > NumComponents_(World)) EXIT
	   if((comp_id) == World%allids(n)) then
	      valid = .true.
	      EXIT
	   endif

	end do CHECK_COMP_ID_LOOP

	if(.not. valid) then
	   write(stderr,'(2a,1i)') myname,":: invalid component id no. = ",&
		comp_id
	   call MP_perr_die(myname_,'invalid comp_id = ',comp_id)
	endif

      ! Check argument comp_rank for validity on the communicator associated
      ! with comp_id.  Assume initialy it is invalid.

	valid = .false.

	if((0 <= comp_rank) .or. &
	     (comp_rank < ComponentNumProcs_(World, comp_id))) then
	   valid = .true.
	endif

	if(.not. valid) then
	   write(stderr,'(2a,1i,1a,1i)') myname,":: invalid process ID. = ", &
		comp_rank, "on component ",comp_id
	   call MP_perr_die(myname_,'invalid comp_rank = ',comp_rank)
	endif

     endif ! if(check)

  endif ! if((present(check)) .and. check)

      ! If we have reached this point, the input data are valid.
      ! Return the global rank for comp_rank on component comp_id

  world_rank = World%idGprocid(comp_id, comp_rank)

  if(world_rank < 0) then
     write(stderr,'(2a,1i)') myname,":: negative world rank = ",world_rank
     call MP_perr_die(myname_,'negative world rank = ',world_rank)
  endif    

  ComponentToWorldRank_ = world_rank

 end function ComponentToWorldRank_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ComponentRootRank_ - Rank of component root on COMM_WORLD.
!
! !DESCRIPTION:
! The function {\tt ComponentRootRank\_} takes an input component ID 
! {\tt comp\_id} and input {\tt MCTWorld} variable {\tt World}, and
! returns the global rank of the root of this component.  The optional 
! {\tt LOGICAL} argument {\tt check} is used to control exhaustive 
! argument validity checking against {\tt World}.  If the argument 
! {\tt check} is not provided, no checking, which can hamper performance 
! is done.
!
! !INTERFACE:

 integer function ComponentRootRank_(comp_id, World, check)
!
! !USES:
!
      use m_die
      use m_stdio

      implicit none

      integer, intent(in)	     :: comp_id   ! component id
      type(MCTWorld), intent(in)     :: World     ! World
      logical, intent(in), optional  :: check     ! exhaustive checking
                                                  ! flag

! !REVISION HISTORY:
!       05Feb01 - J. Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::ComponentRootRank_'

  integer :: world_comp_root

      ! Call ComponentToWorldRank_ assuming the root on a remote component
      ! has rank zero on the communicator associated with that component.

  if(present(check)) then
     world_comp_root = ComponentToWorldRank_(0, comp_id, World, check)
  else
     world_comp_root = ComponentToWorldRank_(0, comp_id, World)
  endif

  if(world_comp_root < 0) then
     write(stderr,'(2a,1i)') myname,":: negative world rank = ",& 
	  world_comp_root
     call MP_perr_die(myname_,'invalid root id = ',world_comp_root)
  endif    

  ComponentRootRank_ = world_comp_root

 end function ComponentRootRank_

 end module m_MCTWorld

