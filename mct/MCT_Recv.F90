!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: MCT_Recv
!
! !DESCRIPTION:
! Recieve into the AttrVect the data coming from the component
! specified in the Router.  An Error will result if the
! attribute list of the incoming data doesnt match any of
! the attributes in the argument AttrVect.
! Requires a corresponding MCT\_Send to be called on the other component.
!
! !INTERFACE:

 subroutine MCT_Recv(aV, Rout)
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
      Type(AttrVect),intent(inout) :: 	aV
      Type(Router),intent(in) ::	Rout

! !REVISION HISTORY:
!      07Feb01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
!      08Feb01 - R. Jacob <jacob@mcs.anl.gov> - first working code
!      18May01 - R. Jacob <jacob@mcs.anl.gov> - use MP_Type to
!                determine type in mpi_recv
!      07Jun01 - R. Jacob <jacob@mcs.anl.gov> - remove logic to
!                check "direction" of Router.  remove references
!                to ThisMCTWorld%mylrank
!      03Aug01 - E.T. Ong <eong@mcs.anl.gov> - explicity specify starting
!                address in MPI_RECV
!      27Nov01 - E.T. Ong <eong@mcs.anl.gov> - deallocated to prevent
!                memory leaks
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_='MCT_Recv'
  integer ::    numi,numr,i,proc,buffsize,j,k
  integer ::    ier,nseg,mycomp,othercomp,tag
  integer :: mp_Type_rp1
  integer,dimension(:),pointer  :: ireqs,rreqs
  integer,dimension(:,:),allocatable    :: istatus,rstatus

! declare a pointer structure for the real data
  type :: rptr
    real,dimension(:),pointer :: pr
  end type

! declare a pointer structure for the integer data
  type :: iptr
    integer,dimension(:),pointer :: pi
  end type

! declare arrays of pointers to hold data to be received
  type(rptr),dimension(:),allocatable :: rp1
  type(iptr),dimension(:),allocatable :: ip1

!--------------------------------------------------------

  mycomp=Rout%comp1id
  othercomp=Rout%comp2id


!  find total number of real and integer vectors
! for now, assume we are receiving all of them
  numi = nIAttr(aV)
  numr = nRAttr(aV)

!!!!!!!!!!!!!! IF RECEIVING INTEGER DATA
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

! allocate buffers to hold all incoming data
   do proc=1,Rout%nprocs
    allocate(ip1(proc)%pi(Rout%locsize(proc)*numi),stat=ier)
    if(ier/=0) call MP_perr_die(myname_,'allocate(ip1%pi)',ier)
   enddo
  endif

!!!!!!!!!!!!!! IF RECEIVING REAL DATA
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

! allocate buffers to hold all incoming data
   do proc=1,Rout%nprocs
    allocate(rp1(proc)%pr(Rout%locsize(proc)*numr),stat=ier)
    if(ier/=0) call MP_perr_die(myname_,'allocate(rp1%pr)',ier)
   enddo

  endif
! if(myproc==0) then
! do proc=1,Rout%nprocs
! write(*,*)'Recv',size(rp1(proc)%pr),Rout%pe_list(proc)
! enddo
! endif


! receive the integer data
  if(numi .ge. 1) then
    do proc=1,Rout%nprocs

! corresponding tag logic must be in MCT_Send
      tag = 100000*othercomp + 1000*Rout%pe_list(proc) + &
            500 + ThisMCTWorld%mygrank
      call MPI_IRECV(ip1(proc)%pi(1),Rout%locsize(proc)*numi,MP_INTEGER,&
         Rout%pe_list(proc),tag,MP_COMM_WORLD,ireqs(proc),ier)
      if(ier /= 0) call MP_perr_die(myname_,'MPI_IRECV(ints)',ier)
    enddo
  endif

! receive the real data
  if(numr .ge. 1) then
    mp_Type_rp1=MP_Type(rp1(1)%pr(1))
    do proc=1,Rout%nprocs

! corresponding tag logic must be in MCT_Send
      tag = 100000*othercomp + 1000*Rout%pe_list(proc) + &
            700 + ThisMCTWorld%mygrank
!     write(*,*)"RECVR",ThisMCTWorld%mygrank,Rout%pe_list(proc),tag
      call MPI_IRECV(rp1(proc)%pr(1),Rout%locsize(proc)*numr,mp_Type_rp1,&
         Rout%pe_list(proc),tag,MP_COMM_WORLD,rreqs(proc),ier)
      if(ier /= 0) call MP_perr_die(myname_,'MPI_IRECV(reals)',ier)
    enddo
  endif


! wait for all recieves to complete
!  should redo this to load data as its received -RLJ
  if(numi .ge. 1)  then
    call MPI_WAITALL(Rout%nprocs,ireqs,istatus,ier)
    if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL(ints)',ier)
   endif

  if(numr .ge. 1) then
    call MPI_WAITALL(Rout%nprocs,rreqs,rstatus,ier)
    if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL(reals)',ier)
  endif

! Load data which came from each processor
  do proc=1,Rout%nprocs


! Load the integer data
    if(numi .ge. 1) then

      k=1
! load all integer vectors from the receive buffer
      do j=1,numi

! load the correct pieces of this particular integer vector
       do nseg=1,Rout%num_segs(proc)
         do i=0,Rout%seg_lengths(proc,nseg)-1
           aV%iAttr(j,Rout%seg_starts(proc,nseg)+i)=ip1(proc)%pi(k)
           k=k+1
         enddo
       enddo

      enddo
    endif

! Load the real data
    if(numr .ge. 1) then

      k=1
! load all real vectors from the receive buffer
      do j=1,numr

! load the correct pieces of this particular real vector
       do nseg=1,Rout%num_segs(proc)
         do i=0,Rout%seg_lengths(proc,nseg)-1
           aV%rAttr(j,Rout%seg_starts(proc,nseg)+i)=rp1(proc)%pr(k)
           k=k+1
         enddo
       enddo

      enddo
    endif
  enddo

    ! Deallocate memory to prevent leaks!

    if(numi .ge. 1) then
       do proc=1,Rout%nprocs
	  deallocate(ip1(proc)%pi,stat=ier)
	  if(ier/=0) call MP_perr_die(myname_,'deallocate(ip1%pi)',0)
       enddo

       deallocate(ip1,stat=ier)
       if(ier/=0) call MP_perr_die(myname_,'deallocate(ip1)',0)
       
       deallocate(ireqs,stat=ier)
       if(ier/=0) call MP_perr_die(myname_,'deallocate(ireqs)',0)

       deallocate(istatus,stat=ier)
       if(ier/=0) call MP_perr_die(myname_,'deallocate(istatus)',0)
    endif

    if(numr .ge. 1) then
       do proc=1,Rout%nprocs
	  deallocate(rp1(proc)%pr,stat=ier)
	  if(ier/=0) call MP_perr_die(myname_,'deallocate(rp1%pi)',0)
       enddo

       deallocate(rp1,stat=ier)
       if(ier/=0) call MP_perr_die(myname_,'deallocate(rp1)',0)

       deallocate(rreqs,stat=ier)
       if(ier/=0) call MP_perr_die(myname_,'deallocate(rreqs)',0)

       deallocate(rstatus,stat=ier)
       if(ier/=0) call MP_perr_die(myname_,'deallocate(rstatus)',0)
    endif

end subroutine
