!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Router -- Router class 
!
! !DESCRIPTION:
! The Router data type contains all the information needed
! to send an AttrVect between a component on M MPI-processes and a component
! on N MPI-processes.   This module defines the Router datatype and provides
! methods to create and destroy one.
!
! !INTERFACE:

 module m_Router

      use m_realkinds, only : FP

      implicit none

      private   ! except

! !declare a private pointer structure for the real data
      type :: rptr
        sequence
        real(FP),dimension(:),pointer :: pr
      end type

! !declare a private pointer structure for the integer data
      type :: iptr
        sequence
        integer,dimension(:),pointer :: pi
      end type

! !PUBLIC TYPES:
      public :: Router	        ! The class data structure

      public :: rptr,iptr       ! pointer types used in Router
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

    type Router
      sequence
      integer :: comp1id                           ! myid
      integer :: comp2id                           ! id of second component
      integer :: nprocs	                           ! number of procs to talk to
      integer :: maxsize                           ! maximum amount of data going to a processor
      integer :: lAvsize                           ! The local size of AttrVect which can be 
                                                   ! used with this Router in MCT_Send/MCT_Recv
      integer :: numiatt                           ! Number of integer attributes currently in use
      integer :: numratt                           ! Number of real attributes currently in use
      integer,dimension(:),pointer   :: pe_list    ! processor ranks of send/receive in MCT_comm
      integer,dimension(:),pointer   :: num_segs   ! number of segments to send/receive
      integer,dimension(:),pointer   :: locsize    ! total of seg_lengths for a proc
      integer,dimension(:,:),pointer :: seg_starts ! starting index
      integer,dimension(:,:),pointer :: seg_lengths! total length
      type(rptr),dimension(:),pointer :: rp1       ! buffer to hold real data
      type(iptr),dimension(:),pointer :: ip1       ! buffer to hold integer data
      integer,dimension(:),pointer   :: ireqs,rreqs  ! buffer for MPI_Requests
      integer,dimension(:,:),pointer :: istatus,rstatus  ! buffer for MPI_Status
    end type Router

! !PUBLIC MEMBER FUNCTIONS:
      public :: init            ! Create a Router
      public :: clean           ! Destroy a Router


    interface init  ; module procedure  &
        initd_, &       ! initialize a Router between two seperate components
        initp_ 	        ! initialize a Router locally with two GSMaps
    end interface
    interface clean ; module procedure clean_ ; end interface

! !REVISION HISTORY:
! 15Jan01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
! 08Feb01 - R. Jacob <jacob@mcs.anl.gov> add locsize and maxsize 
!           to Router type
! 25Sep02 - R. Jacob <jacob@mcs.anl.gov> Remove type string.  Add lAvsize
! 23Jul03 - R. Jacob <jacob@mcs.anl.gov> Add status and reqs arrays used
!           in send/recv to the Router datatype.
! 24Jul03 - R. Jacob <jacob@mcs.anl.gov> Add real and integer buffers
!           for send/recv to the Router datatype.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='MCT::m_Router'

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
! to build a Router {\tt Rout} between them.
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

! !INPUT PARAMETERS:
!
      integer, intent(in)	       :: othercomp
      integer, intent(in)	       :: mycomm
      type(GlobalSegMap),intent(in)    :: GSMap     ! of the calling comp

! !OUTPUT PARAMETERS:
!
      type(Router), intent(out)        :: Rout

! !REVISION HISTORY:
! 15Jan01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
! 06Feb01 - R. Jacob <jacob@mcs.anl.gov> - Finish initialization
!           of the Router.  Router now works both ways.
! 25Apr01 - R. Jacob <jacob@mcs.anl.gov> - Eliminate early 
!           custom code to exchange GSMap components and instead
!           the more general purpose routine in m_ExchangeMaps.
!           Use new subroutine OrderedPoints in m_GlobalSegMap
!           to construct the vector of local and remote GSMaps.
!           Clean-up code a little.
! 03May01 - R. Jacob <jacob@mcs.anl.gov> - rename to initd and
!           move most of code to new initp routine
!
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::initd_'

  type(GlobalSegMap)    :: RGSMap  !  the other GSMap
  integer ::		   ier

!--------------------------begin code-----------------------

!!!!!!!!!!!!!!!!!Exchange of global map data 

  call MCT_ExGSMap(GSMap,mycomm,RGSMap,othercomp,ier)
  if(ier /= 0) call die(myname_,'ExGSMap',ier)

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
! Router {\tt Rout} between them.  Use local communicator {\tt mycomm}.
!
! !INTERFACE:

 subroutine initp_(GSMap,RGSMap,mycomm,Rout )
!
! !USES:
!
      use m_GlobalSegMap,  only :GlobalSegMap
      use m_GlobalSegMap,  only :OrderedPoints
      use m_GlobalSegMap,  only :ProcessStorage
      use m_GlobalSegMap,  only : GSMap_comp_id => comp_id
      use m_GlobalToLocal, only :GlobalToLocalIndex
      use m_MCTWorld,      only :MCTWorld
      use m_MCTWorld,      only :ThisMCTWorld
      use m_mpif90
      use m_die

      implicit none

! !INPUT PARAMETERS:
!
      type(GlobalSegMap), intent(in)	:: GSMap
      type(GlobalSegMap), intent(in)	:: RGSMap
      integer	     ,    intent(in)	:: mycomm

! !OUTPUT PARAMETERS:
!
      type(Router),      intent(out)	:: Rout

! !REVISION HISTORY:
! 03May01 - R.L. Jacob <jacob@mcs.anl.gov> - Initial code brought
!           in from old init routine.
! 31Jul01 - Jace A Mogill <mogill@cray.com>
!           Rewrote to reduce number of loops and temp storage
!EOP -------------------------------------------------------------------

  character(len=*),parameter :: myname_=myname//'::initp_'
  integer			     :: ier,i,j,k,m
  integer			     :: mysize,myPid,othercomp
  integer			     :: lmaxsize,totallength
  integer                            :: maxsegcount,count
  logical, dimension(:), allocatable :: tmppe_list
  integer, dimension(:,:), pointer   :: tmpsegcount,tmpsegstart

  integer :: procn          ! Index over processors
  integer :: my_segn        ! Index into local segments
  integer :: r_segn         ! Index into remote segments
  integer :: my_left        ! Left point in local segment (global memory)
  integer :: my_right       ! Right point in local segment (global memory)
  integer :: r_left         ! Left point in remote segment (global memory)
  integer :: r_right        ! Right point in remote segment (global memory)
  integer :: v_left         ! Leftmost point in overlap (local memory)
  integer :: v_right        ! Rightmost point in overlap (local memory)
  integer :: previous_right ! Rightmost point in previous overlapped segment (local memory)
  integer :: nsegs_overlap  ! Number of segments that overlap between two procs


  call MP_comm_rank(mycomm,myPid,ier)
  if(ier/=0) call MP_perr_die(myname_,'MP_comm_rank',ier)

  mysize = ProcessStorage(GSMap,myPid)
  othercomp = GSMap_comp_id(RGSMap)

! allocate space for searching
  allocate(tmpsegcount(ThisMCTWorld%nprocspid(othercomp),mysize),&
           tmpsegstart(ThisMCTWorld%nprocspid(othercomp),mysize),&
 	   tmppe_list(ThisMCTWorld%nprocspid(othercomp)),stat=ier)
  if(ier/=0) call die(myname_,'allocate(tmpsegcount,..)',ier)

  tmpsegcount=0
  tmpsegstart=0
  count =0
  maxsegcount=0
!.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  

  do procn = 1, ThisMCTWorld%nprocspid(othercomp)
     nsegs_overlap = 0
     tmppe_list(procn) = .FALSE.
     do my_segn = 1, GSMap%ngseg
        if(GSMap%pe_loc(my_segn) == myPid) then
           do r_segn = 1, RGSMap%ngseg
              if(RGSMap%pe_loc(r_segn) == procn - 1) then
                 my_left = GSMap%start(my_segn)
                 my_right= GSMap%length(my_segn) + my_left - 1
                 r_left  = RGSMap%start(r_segn)
                 r_right = RGSMap%length(r_segn) + r_left - 1
                 if( .not. (my_right < r_left   .or.  &
                            my_left  > r_right  .or.  &
                            my_left  > r_right  .or.  &
                            my_right < r_left) ) then
                    if(nsegs_overlap == 0) then
                       count = count + 1
		       v_right = -9999
                       tmppe_list(procn) = .TRUE.
                    endif
		    v_left = GlobalToLocalIndex(GSMap,max(my_left, r_left),mycomm)
		    previous_right = v_right
		    v_right = v_left + min(my_right, r_right) - max(my_left, r_left)
		    if(v_left == previous_right+1) then
		       tmpsegcount(count, nsegs_overlap) = &
		       tmpsegcount(count, nsegs_overlap) + v_right - v_left + 1
		    else
		       nsegs_overlap = nsegs_overlap + 1
		       tmpsegstart(count, nsegs_overlap) = v_left
		       tmpsegcount(count, nsegs_overlap) = v_right - v_left + 1
		    endif
                 endif
              endif
           enddo
        endif
     enddo
     maxsegcount = max(nsegs_overlap, maxsegcount)
  enddo

!.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  


!!!!!!!!!!!!!!!!!!!!end of search through remote GSMap

! start loading up the Router with data

  Rout%comp1id = GSMap_comp_id(GSMap)
  Rout%comp2id = othercomp
  Rout%nprocs = count
  Rout%numiatt = 0
  Rout%numratt = 0

  allocate(Rout%pe_list(count),Rout%num_segs(count), &
    Rout%seg_starts(count,maxsegcount), &
    Rout%seg_lengths(count,maxsegcount), &
    Rout%locsize(count),stat=ier)
  if(ier/=0) call die(myname_,'allocate(Rout..)',ier)

  allocate(Rout%istatus(MP_STATUS_SIZE,count), &
             Rout%rstatus(MP_STATUS_SIZE,count), &
	     Rout%rreqs(count),Rout%ireqs(count),stat=ier)
  if(ier/=0) call die(myname_,'allocate(status,reqs,...)',ier)

! allocate the number of pointers needed
  allocate(Rout%ip1(count),stat=ier)
  if(ier/=0) call die(myname_,'allocate(ip1)',ier)

! allocate the number of pointers needed
  allocate(Rout%rp1(count),stat=ier)
  if(ier/=0) call die(myname_,'allocate(rp1)',ier)


    
  m=0
  do i=1,ThisMCTWorld%nprocspid(othercomp)
      if(tmppe_list(i))then 
      m=m+1
      ! load processor rank in MCT_comm
      Rout%pe_list(m)=ThisMCTWorld%idGprocid(othercomp,i-1)
      endif
    enddo

    lmaxsize=0
    do i=1,count
      totallength=0
      do j=1,maxsegcount
	if(tmpsegcount(i,j) /= 0) then
	 Rout%num_segs(i)=j
 	 Rout%seg_starts(i,j)=tmpsegstart(i,j)
	 Rout%seg_lengths(i,j)=tmpsegcount(i,j)
	 totallength=totallength+Rout%seg_lengths(i,j)
	endif
      enddo
      Rout%locsize(i)=totallength
      lmaxsize=MAX(lmaxsize,totallength)
    enddo

    Rout%maxsize=lmaxsize
    Rout%lAvsize=mysize

      
  deallocate(tmpsegstart,tmpsegcount,tmppe_list,stat=ier)
  if(ier/=0) call die(myname_,'deallocate()',ier)

 end subroutine initp_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - Destroy a Router
!
! !DESCRIPTION:
! Deallocate Router internal data structures and set integer parts to zero.
!
! !INTERFACE:

    subroutine clean_(Rout,stat)
!
! !USES:
!
      use m_die

      implicit none

!INPUT/OUTPUT PARAMETERS:
      type(Router),      intent(inout) :: Rout

!OUTPUT PARAMETERS:
      integer, optional, intent(out)   :: stat

! !REVISION HISTORY:
! 15Jan01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
! 08Feb01 - R. Jacob <jacob@mcs.anl.gov> - add code to clean
!           the maxsize and locsize
! 01Mar02 - E.T. Ong <eong@mcs.anl.gov> removed the die to prevent
!           crashes and added stat argument.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  deallocate(Rout%pe_list,Rout%num_segs,Rout%seg_starts, &
  Rout%locsize,Rout%seg_lengths,stat=ier)

  deallocate(Rout%rreqs,Rout%ireqs,Rout%rstatus,&
   Rout%istatus,stat=ier)

  deallocate(Rout%ip1,Rout%rp1,stat=ier)

  if(present(stat)) then
     stat=ier
  else
     if(ier /= 0) call warn(myname_,'deallocate(Rout%...)',ier)
  endif

  Rout%comp1id = 0
  Rout%comp2id = 0
  Rout%nprocs = 0
  Rout%maxsize = 0
  Rout%lAvsize = 0
  Rout%numiatt = 0
  Rout%numratt = 0


 end subroutine clean_

 end module m_Router

