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
      public :: ThisMCTWorld 

!  
    type MCTWorld
      integer :: myid	       ! my unique id 
      integer :: ncomps	       !  number of components
      integer :: mynprocs      !  number of procs
      integer,dimension(:),pointer :: allids	       ! unique id for each component
      integer,dimension(:),pointer :: nprocspid	       ! number of processes each component is on
      integer,dimension(:,:),pointer :: idGprocid 	       ! translate between local and global ranks for each component
    end type MCTWorld

    type(MCTWorld) :: ThisMCTWorld	! declare an MCTWorld

    interface init  ; module procedure init_  ; end interface
    interface clean ; module procedure clean_ ; end interface

! !REVISION HISTORY:
!      19Jan01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
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

 end module m_MCTWorld

