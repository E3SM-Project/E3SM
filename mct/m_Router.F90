!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Router -- Router class 
!
! !DESCRIPTION:
! The Router data type contains all the information needed
! to send an AttrVect from a component on M MPI-processes to a component
! on N MPI-processes.  
!
! !INTERFACE:

 module m_Router
!
! !USES:
!
      use m_Navigator, only : Navigator
      use m_Navigator, only : Navigator_init => init
      use m_Navigator, only : Navigator_clean => clean
!     use MPH_module

      implicit none

      private   ! except

      public :: Router	        ! The class data structure
      public :: init            ! Create a Router
      public :: clean           ! Destroy a Router

!% if I m the sending model, then on return pe_list is the
!% processor ranks in the remote model to send to, num_segs
!% is how many segments of my GlobalSegMap to send to each, seg_starts
!% is the starting place for each segment and each destination processor
!% and seg_lengths is the total length of each segment, for each processor
!
!% if I m the receiving model, the on return pe_list is the
!% processor ranks in the remote model to receive from, num_segs
!% is how many segments to receive, seg_starts is the starting place
!% for each segment and each sending processor
!  
    type Router
      integer :: sendid	       ! id of sending component
      integer :: recvid	       ! id of receiving component
      character*4 :: type	       ! 'send' or 'recv'
      integer,dimension(:),pointer :: pe_list ! processor ranks of send/receive
      integer,dimension(:),pointer :: num_segs ! number of segments to send/receive
      integer,dimension(:,:),pointer :: seg_starts ! starting index
      integer,dimension(:,:),pointer :: seg_lengths ! total length

      real,dimension(:,:),pointer :: Rbuffer
      integer,dimension(:,:),pointer :: Ibuffer
      Type(Navigator) ::  buffer_nav
    end type Router


    interface init  ; module procedure init_  ; end interface
    interface clean ; module procedure clean_ ; end interface

! !REVISION HISTORY:
!      15Jan01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
!      22Jan01 - J. Larson <larson@mcs.anl.gov> - minor modification
!                for port to SunOS platform:  made more explicit the
!                use blocks for m_Navigator to alleviate confusion in
!                interface declarations.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_Router'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize a Router
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine init_(comp1,comp2,GSMap,mycomm,Rout )
!
! !USES:
!
      use m_GlobalSegMap, only :GlobalSegMap
      use m_GlobalSegMap, only :nlseg
      use m_GlobalSegMap, only :lsize
      use m_MCTWorld,only :MCTWorld
      use m_MCTWorld,only :ThisMCTWorld
      use m_mpif90
      use m_die

      implicit none

      type(Router), intent(out)        :: Rout
      integer, intent(in)	       :: comp1
      integer, intent(in)	       :: comp2
      integer, intent(in)	       :: mycomm
      type(GlobalSegMap),intent(in)    :: GSMap  ! of the calling comp

! !REVISION HISTORY:
!       15Jan01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
!       01Feb01 - R. Jacob <jacob@mcs.anl.gov> - initialize some parts
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::init_'
  integer			:: myPid,ier,totalsegs,i,j,k
  integer			:: mysize,recvsize,rngseg,tag
  integer			:: sendsize
  integer 	:: status(MP_STATUS_SIZE)
  integer,dimension(:),allocatable :: reqs,rstarts,rlengths
  integer,dimension(:),allocatable :: rpe_locs
  integer,dimension(:),allocatable :: mystarts,mylengths
  integer,dimension(:),allocatable :: myvector

  integer,dimension(:,:),pointer :: recvstarts
  integer,dimension(:),pointer :: precvstarts

  call MP_comm_rank(mycomm,myPid,ier)
  call MP_comm_size(mycomm,mysize,ier)

  totalsegs = nlseg(GSMap,myPid)

! build a mystart and mylength array on all processors
  allocate(mystarts(totalsegs),mylengths(totalsegs), &
   myvector(lsize(GSMap)), stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(mystarts,..)',ier)


  j=1
  do i=1,GSMap%ngseg
    if(GSMap%pe_loc(i)==myPid) then
      mystarts(j)=GSMap%start(i)
      mylengths(j)=GSMap%length(i)
      j=j+1
    endif
  enddo

! form one long vector which is all local GSMap points
  do i=1,lsize(GSMap)
   do j=1,totalsegs
    do k=1,mylengths(j)
    myvector(i)=mystarts(j)+k-1
    enddo
   enddo
  enddo


!NOTE  at this point the init routine follows two different
! paths.  The sending component will first learn about the
! receiving components global seg map and figure out what
! portion of its own map it must send.  It will then tell
! the receiving component what its going to send.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! comp1 is the sender
  if(ThisMCTWorld%myid == comp1) then

    Rout%sendid = comp1
    Rout%recvid = comp2
    Rout%type = "send"
    recvsize=ThisMCTWorld%nprocspid(comp2)

! local roots on sender and reciver talk to each other
    if(myPid ==0) then

! learn number of global segments on receiver
       call MPI_RECV(rngseg,1,MP_INTEGER,ThisMCTWorld%idGprocid(comp2,0), &
	3010,MP_COMM_WORLD,status,ier)
       if(ier /= 0) call MP_perr_die(myname_,'MPI_RECV(rngseg)',ier)
    endif
    call MPI_BCAST(rngseg,1,MP_INTEGER,0,mycomm,ier)
    if(ier /= 0) call MP_perr_die(myname_,'MPI_BCAST(rngseg)',ier)

    allocate(rstarts(rngseg),rlengths(rngseg),rpe_locs(rngseg),stat=ier)
    if(ier/=0) call MP_perr_die(myname_,'allocate(rngseg)',ier)

! get starts of receivers GSMap
    if(myPid ==0) then

       call MPI_RECV(rstarts,rngseg,MP_INTEGER,ThisMCTWorld%idGprocid(comp2,0), &
	3020,MP_COMM_WORLD,status,ier)
       if(ier /= 0) call MP_perr_die(myname_,'MPI_RECV(rstarts)',ier)
    endif
    call MPI_BCAST(rstarts,rngseg,MP_INTEGER,0,mycomm,ier)
    if(ier /= 0) call MP_perr_die(myname_,'MPI_BCAST(rstarts)',ier)

! get lengths of receivers GSMap
    if(myPid ==0) then

       call MPI_RECV(rlengths,rngseg,MP_INTEGER,ThisMCTWorld%idGprocid(comp2,0), &
	3030,MP_COMM_WORLD,status,ier)
       if(ier /= 0) call MP_perr_die(myname_,'MPI_RECV(rlengths)',ier)
    endif
    call MPI_BCAST(rlengths,rngseg,MP_INTEGER,0,mycomm,ier)
    if(ier /= 0) call MP_perr_die(myname_,'MPI_BCAST(rlengths)',ier)

! get pe_locs of receivers GSMap
    if(myPid ==0) then

       call MPI_RECV(rpe_locs,rngseg,MP_INTEGER,ThisMCTWorld%idGprocid(comp2,0), &
	3040,MP_COMM_WORLD,status,ier)
       if(ier /= 0) call MP_perr_die(myname_,'MPI_RECV(rlengths)',ier)
    endif
    call MPI_BCAST(rpe_locs,rngseg,MP_INTEGER,0,mycomm,ier)
    if(ier /= 0) call MP_perr_die(myname_,'MPI_BCAST(rlengths)',ier)

! begin looping through other components processors
! looking for points in common




  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! comp2 is the receiver
  else if(ThisMCTWorld%myid == comp2) then
  Rout%sendid = comp1
  Rout%recvid = comp2
  Rout%type = "recv"


! local roots on sender and receiver talk to each other
    if(myPid ==0) then

! send ngseg to sending component
     call MPI_SEND(GSMap%ngseg,1,MP_INTEGER, &
     ThisMCTWorld%idGprocid(comp1,0),3010,MP_COMM_WORLD,status,ier)
     if(ier /= 0) call MP_perr_die(myname_,'MPI_SEND(ngseg)',ier)

! send starts
     call MPI_SEND(GSMap%start,GSMap%ngseg,MP_INTEGER, &
     ThisMCTWorld%idGprocid(comp1,0),3020,MP_COMM_WORLD,status,ier)
     if(ier /= 0) call MP_perr_die(myname_,'MPI_SEND(start)',ier)

! send lengths
     call MPI_SEND(GSMap%length,GSMap%ngseg,MP_INTEGER, &
     ThisMCTWorld%idGprocid(comp1,0),3030,MP_COMM_WORLD,status,ier)
     if(ier /= 0) call MP_perr_die(myname_,'MPI_SEND(length)',ier)

! send pe_locs
     call MPI_SEND(GSMap%pe_loc,GSMap%ngseg,MP_INTEGER, &
     ThisMCTWorld%idGprocid(comp1,0),3040,MP_COMM_WORLD,status,ier)
     if(ier /= 0) call MP_perr_die(myname_,'MPI_SEND(pe_loc)',ier)

    endif


  endif


 end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - Destroy a Router
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(Rout)
!
! !USES:
!
      use m_die

      implicit none

      type(Router), intent(inout) :: Rout

! !REVISION HISTORY:
!       15Jan01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
!       31Jan01 - R. Jacob <jacob@mcs.anl.gov> - actual code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  deallocate(Rout%pe_list,Rout%num_segs,Rout%seg_starts, stat=ier)
  if(ier /= 0) call perr_die(myname_,'deallocate(Rout,...)',ier)

  Rout%sendid = 0
  Rout%recvid = 0
  Rout%type = ""


 end subroutine clean_

 end module m_Router

