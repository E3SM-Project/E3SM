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
      integer :: nprocs	       ! number of procs to talk to
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
!       02Feb01 - R. Jacob <jacob@mcs.anl.gov> - initialize the send
!		portion
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::init_'
  integer			:: myPid,ier,nlsegs,i,j,k,m
  integer			:: mysize,recvsize,rngseg,tag,pr
  integer			:: sendsize,mysize,count
  integer 	:: status(MP_STATUS_SIZE)
  integer,dimension(:),allocatable :: reqs,rstarts,rlengths
  integer,dimension(:),allocatable :: rpe_locs,rlsizes
  integer,dimension(:),allocatable :: mystarts,mylengths
  integer,dimension(:),allocatable :: myvector,lsizes,rvector
  logical,dimension(:),allocatable :: hitparade,tmppe_list

  integer,dimension(:,:),pointer :: recvstarts
  integer,dimension(:,:),pointer :: tmpsegcount,tmpsegstart
  integer,dimension(:),pointer :: precvstarts
  integer,dimension(:),pointer :: beginseg,targproc,numpoints
  logical :: firstpoint
  integer :: maxsegcount,mmax

  call MP_comm_rank(mycomm,myPid,ier)
  call MP_comm_size(mycomm,mysize,ier)
  firstpoint = .TRUE.

  nlsegs = nlseg(GSMap,myPid)
  mysize = lsize(GSMap)

! build a mystart and mylength array on all processors
  allocate(mystarts(nlsegs),mylengths(nlsegs), &
   myvector(mysize), stat=ier)
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
  i=1
   do j=1,nlsegs
    do k=1,mylengths(j)
    myvector(i)=mystarts(j)+k-1
    i=i+1
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

    allocate(rstarts(rngseg),rlengths(rngseg), &
    rlsizes(ThisMCTWorld%nprocspid(comp2)),rpe_locs(rngseg),stat=ier)
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
       if(ier /= 0) call MP_perr_die(myname_,'MPI_RECV(rpe_locs)',ier)
    endif
    call MPI_BCAST(rpe_locs,rngseg,MP_INTEGER,0,mycomm,ier)
    if(ier /= 0) call MP_perr_die(myname_,'MPI_BCAST(rpe_locs)',ier)

! get sizes of receivers GSMap
    if(myPid ==0) then

       call MPI_RECV(rlsizes,ThisMCTWorld%nprocspid(comp2),MP_INTEGER,&
       ThisMCTWorld%idGprocid(comp2,0),3050,MP_COMM_WORLD,status,ier)
       if(ier /= 0) call MP_perr_die(myname_,'MPI_RECV(sizes)',ier)
    endif
    call MPI_BCAST(rlsizes,ThisMCTWorld%nprocspid(comp2),MP_INTEGER,0,mycomm,ier)
    if(ier /= 0) call MP_perr_die(myname_,'MPI_BCAST(sizes)',ier)
!    write(*,*)"ROUTER",myPid,rlsizes

    count =0
    mmax=0
    allocate(hitparade(mysize),stat=ier)
    if(ier/=0) call MP_perr_die(myname_,'allocate(hitparade,..)',ier)
    allocate(tmpsegcount(ThisMCTWorld%nprocspid(comp2),mysize),&
             tmpsegstart(ThisMCTWorld%nprocspid(comp2),mysize),&
    		tmppe_list(ThisMCTWorld%nprocspid(comp2)),stat=ier)
    if(ier/=0) call MP_perr_die(myname_,'allocate(tmpsegcount,..)',ier)
    tmpsegcount=0
    tmpsegstart=0

!  check each processor to see what I need to send
    do i=1,ThisMCTWorld%nprocspid(comp2)

      hitparade=.FALSE.

! build the vector of GSMap points handled on remote proc i-1
      allocate(rvector(rlsizes(i)),stat=ier)
      if(ier/=0) call MP_perr_die(myname_,'allocate(mysize,..)',ier)

      m=1

      do j=1,rngseg
        if(rpe_locs(j)== (i-1)) then
          do k=1,rlengths(j)
            rvector(m)=rstarts(j)+k-1
           m=m+1
          enddo
        endif
      enddo

! check every point in my local vector against every point of
! the remote procs vector.  Record hits.
      do j=1,mysize
	do k=1,rlsizes(i)
	  if(myvector(j)==rvector(k)) hitparade(j)=.TRUE.
        enddo
      enddo

      if(myPid==0)then
! 	write(*,*)"ROUTER",i,hitparade,rlsizes(i),rvector(1),rvector(rlsizes(i))
      endif

! if no hits found, go look at the next processor.
      if(.NOT.ANY(hitparade)) then
        deallocate(rvector,stat=ier)
        tmppe_list(i)=.FALSE.
	CYCLE
      endif
      
! if reached this point, then this processor has some points
      count=count+1
      tmppe_list(i)=.TRUE.

!  now find and count them.
      firstpoint=.TRUE.
      m=0

! search my vector again
      do j=1,mysize
	if(hitparade(j)) then
!       found a segment.  note firstpoint and start counting
	  if(firstpoint) then 
	    m=m+1
	    tmpsegstart(count,m)=myvector(j)
	    tmpsegcount(count,m)=1
	    firstpoint=.FALSE.
          else
	    tmpsegcount(count,m)=tmpsegcount(count,m)+1
	  endif
        else
! if firstpoint is false, then I just finished counting a
! segment and need to reset it in case I find another one.
	  if(.NOT.firstpoint) then
	    firstpoint=.TRUE.
          endif
        endif
      enddo

      if(myPid==0)then
	 do pr=1,m
! 	write(*,*)"ROUTERD",i,pr,tmpsegcount(i,pr),tmpsegstart(i,pr)
	enddo
      endif


      mmax=MAX(m,mmax)
      deallocate(rvector,stat=ier)
    enddo

! start loading up the Router
    Rout%nprocs = count

    maxsegcount=mmax
!   if(myPid==0) write(*,*)"COUNT",count,maxsegcount

    allocate(Rout%pe_list(count),Rout%num_segs(count), &
    Rout%seg_starts(count,maxsegcount), &
    Rout%seg_lengths(count,maxsegcount),stat=ier)
    if(ier/=0) call MP_perr_die(myname_,'allocate(Rout..)',ier)
    
    m=0
    do i=1,ThisMCTWorld%nprocspid(comp2)
      if(tmppe_list(i))then 
      m=m+1
! load processor rank in MPI_COMM_WORLD
      Rout%pe_list(m)=ThisMCTWorld%idGprocid(comp2,i-1)
      endif
    enddo

    do i=1,count
      k=0
      do j=1,mysize
!	if(myPid==0)write(*,*)"RRR",i,j,tmpsegcount(i,j)
	if(tmpsegcount(i,j) /= 0) then
	 k=k+1
	 Rout%seg_starts(i,k)=tmpsegstart(i,j)
	 Rout%seg_lengths(i,k)=tmpsegcount(i,j)
	endif
      enddo
      Rout%num_segs(i)=k
    enddo

    if(myPid==0) then
     do i=1,Rout%nprocs
!      write(*,*)"ROUTERE",i,Rout%pe_list(i),Rout%num_segs(i)
      do j=1,Rout%num_segs(i)
!       write(*,*)"ROUTEREE",i,j,Rout%seg_starts(i,j),Rout%seg_lengths(i,j)
      enddo
     enddo
    endif
       
      
    deallocate(tmpsegstart,tmpsegcount,tmppe_list,stat=ier)




  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! comp2 is the receiver
  else if(ThisMCTWorld%myid == comp2) then
  Rout%sendid = comp1
  Rout%recvid = comp2
  Rout%type = "recv"

    allocate(lsizes(ThisMCTWorld%nprocspid(comp2)),stat=ier)
    if(ier/=0) call MP_perr_die(myname_,'allocate(lsizes,..)',ier)

    
    mysize=lsize(GSMap)
    call MPI_GATHER(mysize,1,MP_INTEGER,lsizes,1,MP_INTEGER,0,mycomm,ier)
    if(ier /= 0) call MP_perr_die(myname_,'MPI_GATHER(lsize)',ier)

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

! send sizes
     call MPI_SEND(lsizes,ThisMCTWorld%nprocspid(comp2),MP_INTEGER, &
     ThisMCTWorld%idGprocid(comp1,0),3050,MP_COMM_WORLD,status,ier)
     if(ier /= 0) call MP_perr_die(myname_,'MPI_SEND(sizes)',ier)


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

