!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Router -- Router class 
!
! !DESCRIPTION:
! The Router data type contains all the information needed
! to send an AttrVect between a component on M MPI-processes and a component
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

!\end{verbatim}
!% On return, pe_list is the processor ranks of the other
!% component to receive from/send to.  num_segs is the
!% number of segments out of my local AttrVect which must
!% be sent/received.  (In general, these wont coincide exactly
!% with the segments used to define the GlobalMap)
!% seg_start is the start *in the local AttrVect* of each segment
!%  (start goes from 1 to lsize(GSMap))
!% and seg_lengths is the length.
!\begin{verbatim}
!  
    type Router
      integer :: comp1id	       ! myid
      integer :: comp2id	       ! id of second component
      character*4 :: type	       ! '1way' or '2way'
      integer :: nprocs	       ! number of procs to talk to
      integer :: maxsize		! maximum amount of data going to a processor
      integer,dimension(:),pointer :: pe_list ! processor ranks of send/receive in MPI_COMM_WORLD
      integer,dimension(:),pointer :: num_segs ! number of segments to send/receive
      integer,dimension(:),pointer :: locsize ! total of seg_lengths for a proc
      integer,dimension(:,:),pointer :: seg_starts ! starting index
      integer,dimension(:,:),pointer :: seg_lengths ! total length

!     real,dimension(:,:),pointer :: Rbuffer
!     integer,dimension(:,:),pointer :: Ibuffer
!     Type(Navigator) ::  buffer_nav
    end type Router


    interface init  ; module procedure init_  ; end interface
    interface clean ; module procedure clean_ ; end interface

! !REVISION HISTORY:
!      15Jan01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
!      22Jan01 - J. Larson <larson@mcs.anl.gov> - minor modification
!                for port to SunOS platform:  made more explicit the
!                use blocks for m_Navigator to alleviate confusion in
!                interface declarations.
!      08Feb01 - R. Jacob <jacob@mcs.anl.gov> add locsize and maxsize 
!                to Router type
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

 subroutine init_(othercomp,GSMap,mycomm,Rout )
!
! !USES:
!
      use m_GlobalSegMap, only :GlobalSegMap
      use m_GlobalSegMap, only :nlseg
      use m_GlobalSegMap, only :lsize
      use m_GlobalToLocal, only :GlobalToLocalIndex
      use m_MCTWorld,only :MCTWorld
      use m_MCTWorld,only :ThisMCTWorld
      use m_mpif90
      use m_die

      implicit none

      type(Router), intent(out)        :: Rout
      integer, intent(in)	       :: othercomp
      integer, intent(in)	       :: mycomm
      type(GlobalSegMap),intent(in)    :: GSMap  ! of the calling comp

! !REVISION HISTORY:
!       15Jan01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
!       01Feb01 - R. Jacob <jacob@mcs.anl.gov> - initialize some parts
!       02Feb01 - R. Jacob <jacob@mcs.anl.gov> - initialize the send
!       06Feb01 - R. Jacob <jacob@mcs.anl.gov> - Finish initialization
!                 of the Router.  Router now works both ways.
!       08Feb01 - R. Jacob <jacob@mcs.anl.gov> - use GlobaltoLocalIndex
!                 to load local index values into Router. Init locsize
!                 and maxsize.  add deallocate statements.
!	22Mar01 - R. Jacob <jacob@mcs.anl.gov> - only use other components
!                 id when initializing
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::init_'
  integer			:: myPid,ier,nlsegs,i,j,k,m
  integer			:: mysize,recvsize,rngseg,tag,pr
  integer			:: sendsize,count,prsize
  integer			:: lmaxsize,totallength
  integer 	:: status(MP_STATUS_SIZE)
  integer,dimension(:),allocatable :: rstarts,rlengths
  integer,dimension(:),allocatable :: rpe_locs,rlsizes,rvector
  integer,dimension(:),allocatable :: mystarts,mylengths
  integer,dimension(:),allocatable :: myvector,lsizes
  logical,dimension(:),allocatable :: hitparade,tmppe_list

  integer,dimension(:,:),pointer :: tmpsegcount,tmpsegstart
  logical :: firstpoint
  integer :: maxsegcount,mmax,othercomp,mycomp

  call MP_comm_rank(mycomm,myPid,ier)
  call MP_comm_size(mycomm,prsize,ier)
  firstpoint = .TRUE.

  nlsegs = nlseg(GSMap,myPid)
  mysize = lsize(GSMap,mycomm)

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

  deallocate(mystarts,mylengths, stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'deallocate(mystarts,..)',ier)

  

  mycomp = ThisMCTWorld%myid
  Rout%comp1id = ThisMCTWorld%myid
  Rout%comp2id = othercomp
  Rout%type = "2way"
  recvsize=ThisMCTWorld%nprocspid(othercomp)

!!!!!!!!!!!!!!!!!Begin exchange of global map data 

! learn number of global segments on othercomp
  if(myPid ==0) then
     call MPI_SENDRECV(GSMap%ngseg,1,MP_INTEGER, &
      ThisMCTWorld%idGprocid(othercomp,0),3010, &
      rngseg,1,MP_INTEGER,ThisMCTWorld%idGprocid(othercomp,0), &
      3010,MP_COMM_WORLD,status,ier)
     if(ier /= 0) call MP_perr_die(myname_,'MPI_SENDRECV(ngseg)',ier)
  endif

  call MPI_BCAST(rngseg,1,MP_INTEGER,0,mycomm,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MPI_BCAST(rngseg)',ier)

! use rngseg to allocate memory for other arrays.
  allocate(rstarts(rngseg),rlengths(rngseg), &
  rlsizes(ThisMCTWorld%nprocspid(othercomp)),rpe_locs(rngseg),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(rngseg)',ier)

! get starts of receivers GSMap
  if(myPid ==0) then
     call MPI_SENDRECV(GSMap%start,GSMap%ngseg,MP_INTEGER, &
      ThisMCTWorld%idGprocid(othercomp,0),3020, &
      rstarts,rngseg,MP_INTEGER,ThisMCTWorld%idGprocid(othercomp,0), &
      3020,MP_COMM_WORLD,status,ier)
     if(ier /= 0) call MP_perr_die(myname_,'MPI_SENDRECV(starts)',ier)
  endif

  call MPI_BCAST(rstarts,rngseg,MP_INTEGER,0,mycomm,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MPI_BCAST(rstarts)',ier)
  
! get lengths of receivers GSMap
  if(myPid ==0) then
     call MPI_SENDRECV(GSMap%length,GSMap%ngseg,MP_INTEGER, &
      ThisMCTWorld%idGprocid(othercomp,0),3030, &
      rlengths,rngseg,MP_INTEGER,ThisMCTWorld%idGprocid(othercomp,0), &
      3030,MP_COMM_WORLD,status,ier)
     if(ier /= 0) call MP_perr_die(myname_,'MPI_SENDRECV(lengths)',ier)
  endif

  call MPI_BCAST(rlengths,rngseg,MP_INTEGER,0,mycomm,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MPI_BCAST(rlengths)',ier)


! get pe_locs of receivers GSMap
  if(myPid ==0) then
     call MPI_SENDRECV(GSMap%pe_loc,GSMap%ngseg,MP_INTEGER, &
      ThisMCTWorld%idGprocid(othercomp,0),3040, &
      rpe_locs,rngseg,MP_INTEGER,ThisMCTWorld%idGprocid(othercomp,0), &
      3040,MP_COMM_WORLD,status,ier)
     if(ier /= 0) call MP_perr_die(myname_,'MPI_SENDRECV(pe_locs)',ier)
  endif

  call MPI_BCAST(rpe_locs,rngseg,MP_INTEGER,0,mycomm,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MPI_BCAST(rpe_locs)',ier)


! gather the local GSMap size on each processor
  allocate(lsizes(ThisMCTWorld%nprocspid(mycomp)),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(lsizes,..)',ier)
    
  call MPI_GATHER(mysize,1,MP_INTEGER,lsizes,1,MP_INTEGER,0,mycomm,ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_GATHER(lsizes)',ier)

! get sizes of receivers GSMap
  if(myPid ==0) then
     call MPI_SENDRECV(lsizes,ThisMCTWorld%nprocspid(mycomp),MP_INTEGER, &
      ThisMCTWorld%idGprocid(othercomp,0),3050,rlsizes, &
      ThisMCTWorld%nprocspid(othercomp),MP_INTEGER, &
      ThisMCTWorld%idGprocid(othercomp,0),3050,MP_COMM_WORLD,status,ier)
     if(ier /= 0) call MP_perr_die(myname_,'MPI_SENDRECV(sizes)',ier)
  endif

  call MPI_BCAST(rlsizes,ThisMCTWorld%nprocspid(othercomp),MP_INTEGER,0,mycomm,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MPI_BCAST(sizes)',ier)


  if(ThisMCTWorld%myid == mycomp) then
!   write(*,*)"mycomp says",myPid,rngseg
  else
!   write(*,*)"othercomp says",myPid,rngseg
  endif


! write(*,*)"ROUTER",myPid,rlsizes

  count =0
  mmax=0
  allocate(hitparade(mysize),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(hitparade,..)',ier)

  allocate(tmpsegcount(ThisMCTWorld%nprocspid(othercomp),mysize),&
           tmpsegstart(ThisMCTWorld%nprocspid(othercomp),mysize),&
 	   tmppe_list(ThisMCTWorld%nprocspid(othercomp)),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(tmpsegcount,..)',ier)

  tmpsegcount=0
  tmpsegstart=0

!  check each processor to see what I have in common
    do i=1,ThisMCTWorld%nprocspid(othercomp)

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

    if((myPid==0) .and. (ThisMCTWorld%myid == mycomp)) then
! 	write(*,*)"ROUTER",i,hitparade,rlsizes(i),rvector(1),rvector(rlsizes(i))
      endif

! if no hits found, go look at the next processor.
      if(.NOT.ANY(hitparade)) then
        deallocate(rvector,stat=ier)
        if(ier/=0) call MP_perr_die(myname_,'deallocate()',ier)
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
 
      if((myPid==0) .and. (ThisMCTWorld%myid == mycomp)) then
        do pr=1,m
! 	 write(*,*)"ROUTERD",i,pr,tmpsegcount(i,pr),tmpsegstart(i,pr)
        enddo
      endif

      mmax=MAX(m,mmax)
      deallocate(rvector,stat=ier)
      if(ier/=0) call MP_perr_die(myname_,'deallocate()',ier)
    enddo

! start loading up the Router with data
    Rout%nprocs = count

    maxsegcount=mmax
!   if((myPid==0) .and. (ThisMCTWorld%myid == mycomp))  write(*,*)"COUNT",count,maxsegcount

    allocate(Rout%pe_list(count),Rout%num_segs(count), &
    Rout%seg_starts(count,maxsegcount), &
    Rout%seg_lengths(count,maxsegcount), &
    Rout%locsize(count),stat=ier)
    if(ier/=0) call MP_perr_die(myname_,'allocate(Rout..)',ier)
    
    m=0
    do i=1,ThisMCTWorld%nprocspid(othercomp)
      if(tmppe_list(i))then 
      m=m+1
! load processor rank in MPI_COMM_WORLD
      Rout%pe_list(m)=ThisMCTWorld%idGprocid(othercomp,i-1)
      endif
    enddo

    lmaxsize=0
    do i=1,count
      k=0
      totallength=0
      do j=1,mysize
!	if(myPid==0)write(*,*)"RRR",i,j,tmpsegcount(i,j)
	if(tmpsegcount(i,j) /= 0) then
	 k=k+1
 	 Rout%seg_starts(i,k)=GlobalToLocalIndex(GSMap,tmpsegstart(i,j),mycomm)
! 	 Rout%seg_starts(i,k)=tmpsegstart(i,j)
	 Rout%seg_lengths(i,k)=tmpsegcount(i,j)
	 totallength=totallength+Rout%seg_lengths(i,k)
	endif
      enddo
      Rout%num_segs(i)=k

      Rout%locsize(i)=totallength
      lmaxsize=MAX(lmaxsize,totallength)
    enddo

    Rout%maxsize=lmaxsize

    if((myPid==0) .and. (ThisMCTWorld%myid == mycomp)) then
     do i=1,Rout%nprocs
!      write(*,*)"ROUTERE",i,Rout%pe_list(i),Rout%num_segs(i),Rout%locsize(i),Rout%maxsize
      do j=1,Rout%num_segs(i)
!       write(*,*)"ROUTEREE",i,j,Rout%seg_starts(i,j),Rout%seg_lengths(i,j)
      enddo
     enddo
    endif
       
      
  deallocate(tmpsegstart,tmpsegcount,tmppe_list, &
    myvector,rstarts,rlengths,rlsizes,lsizes,hitparade,&
    stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'deallocate()',ier)
!endif

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
!       08Feb01 - R. Jacob <jacob@mcs.anl.gov> - add code to clean
!                 the maxsize and locsize
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  deallocate(Rout%pe_list,Rout%num_segs,Rout%seg_starts, &
  Rout%locsize,Rout%seg_lengths,stat=ier)
  if(ier /= 0) call perr_die(myname_,'deallocate(Rout,...)',ier)

  Rout%comp1id = 0
  Rout%comp2id = 0
  Rout%type = ""
  Rout%nprocs = 0
  Rout%maxsize = 0


 end subroutine clean_

 end module m_Router

