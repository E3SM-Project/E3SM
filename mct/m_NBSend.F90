!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_NBSend -- Non blocking send module
!
! !DESCRIPTION:
! Provide support for a nonblocking version of MCT\_Send.  
!
! !INTERFACE:
 module m_NBSend
!
! !USES:

      implicit none

      private	! except

      public :: MCTReqs     ! the reqs datatype

      public :: MCT_ISend   ! the non blocking MCT_Send

      public :: MCT_Wait    ! Wait for the nonblocking send to finish

      type MCTReqs
	integer,dimension(:),pointer :: ireqs   ! the integer sends
	integer,dimension(:),pointer :: rreqs   ! the real sends
      end type MCTReqs

      interface MCT_ISend ; module procedure MCT_ISend_ ; end interface
      interface MCT_Wait ; module procedure MCT_Wait_ ; end interface

! !REVISION HISTORY:
!      07Jun01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_NBSend'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: MCT_ISend_
!
! !DESCRIPTION:
! Send the local AttrVect to another component using a Router.
! Will send entire AttrVect. Requires a corresponding MCT\_Recv 
! or MCT\_IRecv to be called on the other component.
! Returns immediately after posting the send with a MCT\_Req data type.
!
! !INTERFACE:

 subroutine MCT_ISend_(aV, Rout,iList,rList,Reqs)

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
      Type(AttrVect),intent(in) :: 	aV     ! the Attribute vector to send	
      Type(Router),intent(in) ::	Rout   ! the router to use
      Type(MCTReqs),intent(inout) ::	Reqs   ! the returned list of MPI requests
      Type(List),optional,intent(in) ::	iList  ! optional list of integer attributes to send
      Type(List),optional,intent(in) ::	rList  ! optional list of real attributes to send

! !REVISION HISTORY:
!      07Jun01 - R. Jacob <jacob@mcs.anl.gov> - first prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'MCT_ISend_'
  integer ::	numi,numr,i,proc,buffsize,j,k
  integer ::	ier,nseg,mycomp,othercomp,tag,myproc
  integer ::    mp_Type_rp1
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

! call MP_Comm_rank(MP_COMM_WORLD,myproc,ier)
  mycomp=Rout%comp1id
  othercomp=Rout%comp2id

!  find total number of real and integer vectors
! for now, assume we are sending all of them
  numi = nIAttr(aV)
  numr = nRAttr(aV)


!!!!!!!!!!!!!! IF SENDING INTEGER DATA
  if(numi .ge. 1) then

! allocate the number of pointers needed
  allocate(ip1(Rout%nprocs),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(ip1)',ier)

! allocate MPI request array
  allocate(Reqs%ireqs(Rout%nprocs),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(Reqs%ireqs)',ier)
! write(*,*)"Isend",Rout%nprocs,size(Reqs%ireqs)

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
  allocate(Reqs%rreqs(Rout%nprocs),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(Reqs%rreqs)',ier)

! allocate MPI status array
  allocate(rstatus(MP_STATUS_SIZE,Rout%nprocs),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(rstatus)',ier)

! allocate buffers to hold all outgoing data
   do proc=1,Rout%nprocs
    allocate(rp1(proc)%pr(Rout%locsize(proc) *numr),stat=ier)
    if(ier/=0) call MP_perr_die(myname_,'allocate(rp1%pr)',ier)
   enddo

  endif
! if(myproc==0) then
! do proc=1,Rout%nprocs
!  write(*,*)'Send',size(rp1(proc)%pr),Rout%pe_list(proc),ThisMCTWorld%mygrank
! enddo
! endif



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
      call MPI_ISEND(ip1(proc)%pi(1),Rout%locsize(proc)*numi,MP_INTEGER,&
	 Rout%pe_list(proc),tag,MP_COMM_WORLD,Reqs%ireqs(proc),ier)
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
      call MPI_ISEND(rp1(proc)%pr(1),Rout%locsize(proc)*numr,mp_Type_rp1, &
	 Rout%pe_list(proc),tag,MP_COMM_WORLD,Reqs%rreqs(proc),ier)
      if(ier /= 0) call MP_perr_die(myname_,'MPI_ISEND(reals)',ier)
    enddo
  endif

end subroutine MCT_ISend_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: MCT_Wait_
!
! !DESCRIPTION:
! Wait for an MCT\_ISend to complete
!
! !INTERFACE:

 subroutine MCT_Wait_(Rout,Reqs)

!
! !USES:
!
      use m_Router,only  : Router
      use m_mpif90
      use m_die
      use m_stdio

      implicit none
      Type(Router),intent(in) ::        Rout
      Type(MCTReqs),intent(in) ::	Reqs
!      integer,dimension(:),pointer   :: rreqs

! !REVISION HISTORY:
!      07Jun01 - R. Jacob <jacob@mcs.anl.gov> - first prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'MCT_Wait_'
  integer ::	ier
  integer,dimension(:,:),allocatable	:: istatus,rstatus

!--------------------------------------------------------

  if(associated(Reqs%ireqs)) then
! allocate MPI status array
  allocate(istatus(MP_STATUS_SIZE,Rout%nprocs),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(istatus)',ier)
!  write(*,*)"wait",Rout%nprocs,size(Reqs%ireqs)
   call MPI_WAITALL(Rout%nprocs,Reqs%ireqs,istatus,ier)
   if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL(ints)',ier)

! done waiting, free up ireqs
! deallocate(Reqs%ireqs,stat=ier)
! if(ier /= 0) call MP_perr_die(myname_,'deallocate(Reqs%ireqs)',ier)
  
  endif

  if(associated(Reqs%rreqs))   then
! allocate MPI status array
  allocate(rstatus(MP_STATUS_SIZE,Rout%nprocs),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(rstatus)',ier)

   call MPI_WAITALL(Rout%nprocs,Reqs%rreqs,rstatus,ier)
   if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL(reals)',ier)

! done waiting, free up Reqs%rreqs
! call deallocate(Reqs%rreqs,stat=ier)
! if(ier /= 0) call MP_perr_die(myname_,'deallocate(Reqs%rreqs)',ier)

  endif


end subroutine MCT_Wait_

end module m_NBSend
