!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: MCT_Send
!
! !DESCRIPTION:
! Send the local AttrVect to another component using a Router.
! Will send entire AttrVect.
! Requires a corresponding MCT\_Recv to be called on the other component.
!
! !INTERFACE:

 subroutine MCT_Send(aV, Rout,iList,rList)

!
! !USES:
!
      use m_Router,only  : Router
      use m_AttrVect,only : AttrVect
      use m_AttrVect,only : nIAttr,nRAttr
      use m_MCTWorld,only :MCTWorld
      use m_MCTWorld,only :ThisMCTWorld
      use m_list,only:	List
      use m_list,only:  nitem
      use m_mpif90
      use m_die
      use m_stdio

      implicit none
      Type(AttrVect),intent(in) :: 	aV	
      Type(Router),intent(in) ::	Rout
      Type(List),optional,intent(in) ::	iList
      Type(List),optional,intent(in) ::	rList

! !REVISION HISTORY:
!      07Feb01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
!      08Feb01 - R. Jacob <jacob@mcs.anl.gov> - First working code
!      18May01 - R. Jacob <jacob@mcs.anl.gov> - use MP_Type to
!                determine type in mpi_send
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_='MCT_Send'
  integer ::	numi,numr,i,proc,buffsize,j,k
  integer ::	ier,nseg,mycomp,othercomp,tag
  integer ::    mp_Type_rp1
  integer,dimension(:),pointer	:: ireqs,rreqs
  integer,dimension(:,:),allocatable	:: istatus,rstatus

! declare a pointer structure for the real data
  type :: rptr
    real,dimension(:),pointer :: pr
  end type

! declare a pointer structure for the integer data
  type :: iptr
    integer,dimension(:),pointer :: pi
  end type

! declare arrays of pointers to hold data to be sent
  type(rptr),dimension(:),allocatable :: rp1
  type(iptr),dimension(:),allocatable :: ip1

!--------------------------------------------------------

!  I am a sender.  See if router was inialized so that I am
! the sender
  if(ThisMCTWorld%myid == Rout%comp1id) then
    mycomp=Rout%comp1id
    othercomp=Rout%comp2id

! if not then it must be a 2way router
  elseif(ThisMCTWorld%myid == Rout%comp2id) then
     if(Rout%type /= '2way') then
      write(stderr,'(2a)') myname_, &
       'MCTERROR:  Wrong direction with one-way Router',  &
       '  Router is one way from ',Rout%comp1id,' to ',Rout%comp2id
      call MP_perr_die(myname_,'Wrong Way Router',ier)
    endif
    mycomp=Rout%comp2id
    othercomp=Rout%comp1id

! if neither than bad Router
  else
     write(stderr,'(2a)') myname_, &
     'MCTERROR:  Wrong Router passed to MCT_Send.',  &
     '  My compid = ',ThisMCTWorld%myid,             &
     '  Router is between ',Rout%comp1id,' and ',Rout%comp2id
     call MP_perr_die(myname_,'Bad Router',ier)
  endif


!  find total number of real and integer vectors
! for now, assume we are sending all of them
  numi = nIAttr(aV)
  numr = nRAttr(aV)
! if(ThisMCTWorld%mylrank==0)write(*,*)'Send',numi,numr


!!!!!!!!!!!!!! IF SENDING INTEGER DATA
  if(numi .ge. 1) then

! allocate the number of pointers needed
  allocate(ip1(Rout%nprocs),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(ip1)',ier)

! allocate MPI request array
  allocate(ireqs(Rout%nprocs),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(ireqs)',ier)

! allocate MPI status array
  allocate(istatus(MP_STATUS_SIZE,Rout%nprocs),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(istatus)',ier)

! allocate buffers to hold all outgoing data
   do proc=1,Rout%nprocs
    allocate(ip1(proc)%pi(Rout%locsize(proc)*numi),stat=ier)
    if(ier/=0) call MP_perr_die(myname_,'allocate(ip1%pi)',ier)
   enddo
  endif

!!!!!!!!!!!!!! IF SENDING REAL DATA
  if(numr .ge. 1) then

! allocate the number of pointers needed
  allocate(rp1(Rout%nprocs),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(rp1)',ier)

! allocate MPI request array
  allocate(rreqs(Rout%nprocs),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(rreqs)',ier)

! allocate MPI status array
  allocate(rstatus(MP_STATUS_SIZE,Rout%nprocs),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(rstatus)',ier)

! allocate buffers to hold all outgoing data
   do proc=1,Rout%nprocs
    allocate(rp1(proc)%pr(Rout%locsize(proc) *numr),stat=ier)
    if(ier/=0) call MP_perr_die(myname_,'allocate(rp1%pr)',ier)
   enddo

  endif
  if(ThisMCTWorld%mylrank==0) then
  do proc=1,Rout%nprocs
!  write(*,*)'Send',size(rp1(proc)%pr),Rout%pe_list(proc),ThisMCTWorld%mygrank
  enddo
  endif



! Load data going to each processor
  do proc=1,Rout%nprocs


! Load the integer data
    if(numi .ge. 1) then

      k=1
! load all integer vectors to be sent into buffer
      do j=1,numi

! load the correct pieces of this particular integer vector
       do nseg=1,Rout%num_segs(proc)
	 do i=0,Rout%seg_lengths(proc,nseg)-1
           ip1(proc)%pi(k)=aV%iAttr(j,Rout%seg_starts(proc,nseg)+i)
	   k=k+1
	 enddo
       enddo

      enddo
    endif


! Load the real data
    if(numr .ge. 1) then

      k=1
! load all real vectors to be sent into buffer
      do j=1,numr

! load the correct pieces of this particular real vector
       do nseg=1,Rout%num_segs(proc)
	 do i=0,Rout%seg_lengths(proc,nseg)-1
           rp1(proc)%pr(k)=aV%rAttr(j,Rout%seg_starts(proc,nseg)+i)
	   k=k+1
	 enddo
       enddo

      enddo
    endif
  enddo

! send the integer data
  if(numi .ge. 1) then
    do proc=1,Rout%nprocs

! corresponding tag logic must be in MCT_Recv
      tag = 100000*mycomp + 1000*ThisMCTWorld%mygrank + &
	    500 + Rout%pe_list(proc)
      call MPI_ISEND(ip1(proc)%pi,Rout%locsize(proc)*numi,MP_INTEGER,&
	 Rout%pe_list(proc),tag,MP_COMM_WORLD,ireqs(proc),ier)
      if(ier /= 0) call MP_perr_die(myname_,'MPI_ISEND(ints)',ier)
    enddo
  endif

! send the real data
  if(numr .ge. 1) then
    mp_Type_rp1=MP_Type(rp1(1)%pr)
    do proc=1,Rout%nprocs
! corresponding tag logic must be in MCT_Recv
      tag = 100000*mycomp + 1000*ThisMCTWorld%mygrank + &
	    700 + Rout%pe_list(proc)
!     write(*,*)"SENDR",ThisMCTWorld%mygrank,Rout%pe_list(proc),tag
      call MPI_ISEND(rp1(proc)%pr,Rout%locsize(proc)*numr,mp_Type_rp1, &
	 Rout%pe_list(proc),tag,MP_COMM_WORLD,rreqs(proc),ier)
      if(ier /= 0) call MP_perr_die(myname_,'MPI_ISEND(reals)',ier)
    enddo
  endif

! wait for all sends to complete
  if(numi .ge. 1) then
   call MPI_WAITALL(Rout%nprocs,ireqs,istatus,ier)
   if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL(ints)',ier)
  endif

  if(numr .ge. 1)   then
   call MPI_WAITALL(Rout%nprocs,rreqs,rstatus,ier)
   if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL(reals)',ier)
  endif

end subroutine
