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


    interface init  ; module procedure  &
	initd_, &	! initialize a Router between two seperate components
	initp_ 		! initialize a Router locally with two GSMaps
    end interface

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
! !IROUTINE: initd_ - initialize a Router between two seperate components
!
! !DESCRIPTION:
! The routine {\tt initd\_()} exchanges the {\tt GSMap} with the
! component identified by {\tt othercomp} and then calls {\tt initp\_()}
! to build a Router between them.
!
! !INTERFACE:

 subroutine initd_(othercomp,GSMap,mycomm,Rout )
!
! !USES:
!
      use m_GlobalSegMap, only :GlobalSegMap
      use m_ExchangeMaps,only: MCT_ExGSMap => ExchangeMap
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
!	25Apr01 - R. Jacob <jacob@mcs.anl.gov> - Eliminate early 
!                 custom code to exchange GSMap components and instead
!                 the more general purpose routine in m_ExchangeMaps.
!                 Use new subroutine OrderedPoints in m_GlobalSegMap
!                 to construct the vector of local and remote GSMaps.
!                 Clean-up code a little.
!       03May01 - R. Jacob <jacob@mcs.anl.gov> - rename to initd and
!                 move most of code to new initp routine
!
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::initd_'

  type(GlobalSegMap)    :: RGSMap  !  the other GSMap
  integer ::		   ier

!--------------------------begin code-----------------------

!!!!!!!!!!!!!!!!!Exchange of global map data 

  call MCT_ExGSMap(GSMap,mycomm,RGSMap,othercomp,ier)
  if(ier /= 0) call MP_perr_die(myname_,'ExGSMap',ier)

!!!!!!!!!!!!!!!!!Begin comparison of globalsegmaps

  call initp_(GSMap,RGSMap, mycomm, Rout)

 end subroutine initd_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initp_ - initialize a Router from two GlobalSegMaps
!
! !DESCRIPTION:
!
! Given two GlobalSegmentMaps {\tt GSMap} and {\tt RGSMap}, intialize a
! Router between them.
!
! !INTERFACE:

 subroutine initp_(GSMap,RGSMap,mycomm,Rout )
!
! !USES:
!
      use m_GlobalSegMap, only :GlobalSegMap
      use m_GlobalSegMap, only :OrderedPoints
      use m_GlobalSegMap, only :ProcessStorage
      use m_GlobalSegMap, only : GSMap_comp_id => comp_id
      use m_GlobalToLocal, only :GlobalToLocalIndex
      use m_MCTWorld,only :MCTWorld
      use m_MCTWorld,only :ThisMCTWorld
      use m_mpif90
      use m_die

      implicit none
      type(GlobalSegMap), intent(in)	:: GSMap
      type(GlobalSegMap), intent(in)	:: RGSMap
      integer	     , intent(in)	:: mycomm
      type(Router), intent(out)		:: Rout

! !REVISION HISTORY:
!       03May01 - R.L. Jacob <jacob@mcs.anl.gov> - Initial code brought
!                 in from old init routine.
!
!EOP -------------------------------------------------------------------

  character(len=*),parameter :: myname_=myname//'::ExGMapGMap_'
  integer			:: ier,i,j,k,m
  integer			:: mysize,myPid
  integer			:: count
  integer			:: lmaxsize,totallength,rlsize
  integer,dimension(:),pointer	   :: myvector,rvector
  logical,dimension(:),allocatable :: hitparade,tmppe_list
  integer,dimension(:,:),pointer :: tmpsegcount,tmpsegstart
  logical :: firstpoint
  integer :: maxsegcount,mmax,othercomp

  call MP_comm_rank(mycomm,myPid,ier)
  if(ier/=0) call MP_perr_die(myname_,'MP_comm_rank',ier)

  count =0
  mmax=0
  firstpoint = .TRUE.

  mysize = ProcessStorage(GSMap,myPid)
  othercomp = GSMap_comp_id(RGSMap)

! allocate space for searching
  allocate(hitparade(mysize),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(hitparade,..)',ier)

  allocate(tmpsegcount(ThisMCTWorld%nprocspid(othercomp),mysize),&
           tmpsegstart(ThisMCTWorld%nprocspid(othercomp),mysize),&
 	   tmppe_list(ThisMCTWorld%nprocspid(othercomp)),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(tmpsegcount,..)',ier)

  tmpsegcount=0
  tmpsegstart=0

! form the vector of points this processor owns.
  call OrderedPoints(GSMap,myPid,myvector) 

!  check each processor of remote map to see what I have in common
  do i=1,ThisMCTWorld%nprocspid(othercomp)

      hitparade=.FALSE.

! build the vector of GSMap points handled on remote proc i-1

      call OrderedPoints(RGSMap,i-1,rvector) 
      rlsize=ProcessStorage(RGSMap,i-1)


! check every point in my local vector against every point of
! the remote procs vector.  Record hits.
      do j=1,mysize
	do k=1,rlsize
	  if(myvector(j)==rvector(k)) hitparade(j)=.TRUE.
        enddo
      enddo

!     if(myPid==0) then
! 	write(*,*)"ROUTER",i,hitparade,rvector(1),rvector(rlsize)
!     endif

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
!----end of search through my vector again
 
!     if(myPid==0) then
!       do k=1,m
! 	 write(*,*)"ROUTERD",i,k,tmpsegcount(i,k),tmpsegstart(i,k)
!       enddo
!     endif

      mmax=MAX(m,mmax)
      deallocate(rvector,stat=ier)
      if(ier/=0) call MP_perr_die(myname_,'deallocate()',ier)
  enddo

!!!!!!!!!!!!!!!!!!!!end of search through remote GSMap

! start loading up the Router with data

  Rout%comp1id = GSMap_comp_id(GSMap)
  Rout%comp2id = othercomp
  Rout%type = "2way"
  Rout%nprocs = count

    maxsegcount=mmax
!   if(myPid==0)  write(*,*)"COUNT",count,maxsegcount

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

    if(myPid==0) then
     do i=1,Rout%nprocs
!      write(*,*)"ROUTERE",i,Rout%pe_list(i),Rout%num_segs(i),Rout%locsize(i),Rout%maxsize
      do j=1,Rout%num_segs(i)
!       write(*,*)"ROUTEREE",i,j,Rout%seg_starts(i,j),Rout%seg_lengths(i,j)
      enddo
     enddo
    endif
       
      
  deallocate(tmpsegstart,tmpsegcount,tmppe_list, &
    myvector,hitparade,stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'deallocate()',ier)

 end subroutine initp_

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

