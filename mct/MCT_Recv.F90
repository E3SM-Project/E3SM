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

 subroutine MCT_Recv(aV, Rout, iList, rList)
!
! !USES:
!
      use m_MCTWorld, only : MCTWorld
      use m_MCTWorld, only : ThisMCTWorld
      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : nIAttr,nRAttr
      use m_Router,   only : Router
      use m_List,     only : List

      use m_mpif90
      use m_die
      use m_stdio

      implicit none

      Type(AttrVect),       intent(inout) :: aV
      Type(Router),         intent(in)    :: Rout
      Type(List), optional, intent(in)    :: iList
      Type(List), optional, intent(in)    :: rList

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
!      15Feb02 - R. Jacob <jacob@mcs.anl.gov> - Use MCT_comm
!      26Mar02 - E. Ong <eong@mcs.anl.gov> - Apply faster copy order.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_='MCT_Recv'
  integer ::	numi,numr,i,j,k,ier
  integer ::    mycomp,othercomp
  integer ::    AttrIndex,VectIndex,seg_start,seg_end
  integer ::    proc,numprocs,nseg,tag
  integer ::    mp_Type_rp2
  integer, dimension(:), pointer        :: ireqs,rreqs
  integer, dimension(:,:), allocatable  :: istatus,rstatus

! declare a pointer structure for the real data
  type :: rptr
    real,dimension(:),pointer :: pr
  end type

! declare a pointer structure for the integer data
  type :: iptr
    integer,dimension(:),pointer :: pi
  end type

! declare arrays of pointers to hold data to be received
  type(rptr),dimension(:),allocatable :: rp2
  type(iptr),dimension(:),allocatable :: ip2

!--------------------------------------------------------

  mycomp=Rout%comp1id
  othercomp=Rout%comp2id

!  find total number of real and integer vectors
! for now, assume we are receiving all of them
  numi = nIAttr(aV)
  numr = nRAttr(aV)

!!!!!!!!!!!!!! IF SENDING INTEGER DATA
  if(numi .ge. 1) then

! allocate the number of pointers needed
  allocate(ip2(Rout%nprocs),stat=ier)
  if(ier/=0) call die(myname_,'allocate(ip2)',ier)

! allocate buffers to hold all outgoing data
  do proc=1,Rout%nprocs
     allocate(ip2(proc)%pi(Rout%locsize(proc)*numi),stat=ier)
     if(ier/=0) call die(myname_,'allocate(ip2%pi)',ier)
  enddo

! allocate MPI request array
  allocate(ireqs(Rout%nprocs),stat=ier)
  if(ier/=0) call die(myname_,'allocate(ireqs)',ier)

! allocate MPI status array
  allocate(istatus(MP_STATUS_SIZE,Rout%nprocs),stat=ier)
  if(ier/=0) call die(myname_,'allocate(istatus)',ier)

  endif

!!!!!!!!!!!!!! IF RECEIVING REAL DATA
  if(numr .ge. 1) then

! allocate the number of pointers needed
  allocate(rp2(Rout%nprocs),stat=ier)
  if(ier/=0) call die(myname_,'allocate(rp2)',ier)

! allocate buffers to hold all incoming data
  do proc=1,Rout%nprocs
     allocate(rp2(proc)%pr(Rout%locsize(proc)*numr),stat=ier)
     if(ier/=0) call die(myname_,'allocate(rp2%pr)',ier)
  enddo

! allocate MPI request array
  allocate(rreqs(Rout%nprocs),stat=ier)
  if(ier/=0) call die(myname_,'allocate(rreqs)',ier)

! allocate MPI status array
  allocate(rstatus(MP_STATUS_SIZE,Rout%nprocs),stat=ier)
  if(ier/=0) call die(myname_,'allocate(rstatus)',ier)

  mp_Type_rp2=MP_Type(rp2(1)%pr(1))

  endif

  ! Post all MPI_IRECV
  do proc=1,Rout%nprocs

     ! receive the integer data
     if(numi .ge. 1) then

	! corresponding tag logic must be in MCT_Send
	tag = 100000*othercomp + 1000*Rout%pe_list(proc) + &
              500 + ThisMCTWorld%mygrank

	if( Rout%num_segs(proc) > 1 ) then

	   call MPI_IRECV(ip2(proc)%pi(1), &
		Rout%locsize(proc)*numi,MP_INTEGER,Rout%pe_list(proc), &
		tag,ThisMCTWorld%MCT_comm,ireqs(proc),ier)

	else

	   call MPI_IRECV(aV%iAttr(1,Rout%seg_starts(proc,1)), &
		Rout%locsize(proc)*numi,MP_INTEGER,Rout%pe_list(proc), &
		tag,ThisMCTWorld%MCT_comm,ireqs(proc),ier)

	endif

	if(ier /= 0) call MP_perr_die(myname_,'MPI_IRECV(ints)',ier)

     endif

     ! receive the real data
     if(numr .ge. 1) then

	! corresponding tag logic must be in MCT_Send
	tag = 100000*othercomp + 1000*Rout%pe_list(proc) + &
              700 + ThisMCTWorld%mygrank

	if( Rout%num_segs(proc) > 1 ) then

	   call MPI_IRECV(rp2(proc)%pr(1), &
		Rout%locsize(proc)*numr,mp_Type_rp2,Rout%pe_list(proc), &
		tag,ThisMCTWorld%MCT_comm,rreqs(proc),ier)

	else

	   call MPI_IRECV(aV%rAttr(1,Rout%seg_starts(proc,1)), &
		Rout%locsize(proc)*numr,mp_Type_rp2,Rout%pe_list(proc), &
		tag,ThisMCTWorld%MCT_comm,rreqs(proc),ier)

	endif

	if(ier /= 0) call MP_perr_die(myname_,'MPI_IRECV(reals)',ier)

     endif

  enddo

  ! wait for all recieves to complete
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

     if( Rout%num_segs(proc) > 1 ) then

	j=1
	k=1

	! load the correct pieces of the integer and real vectors
	do nseg = 1,Rout%num_segs(proc)
	   seg_start = Rout%seg_starts(proc,nseg)
	   seg_end = seg_start + Rout%seg_lengths(proc,nseg)-1
	   do VectIndex = seg_start,seg_end
	      do AttrIndex = 1,numi
		 aV%iAttr(AttrIndex,VectIndex)=ip2(proc)%pi(j)
		 j=j+1
	      enddo
	      do AttrIndex = 1,numr
		 aV%rAttr(AttrIndex,VectIndex)=rp2(proc)%pr(k)
		 k=k+1
	      enddo
	   enddo
	enddo
	
     endif

  enddo

!........................WAITANY METHOD................................
!
!....NOTE: Make status argument a 1-dimensional array
!  ! Load data which came from each processor
!  do numprocs = 1,Rout%nprocs
!     ! Load the integer data
!     if(numi .ge. 1) then
!	call MPI_WAITANY(Rout%nprocs,ireqs,proc,istatus,ier)
!	if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITANY(ints)',ier)
!	j=1
!	! load the correct pieces of the integer vectors
!	do nseg = 1,Rout%num_segs(proc)
!	   seg_start = Rout%seg_starts(proc,nseg)
!	   seg_end = seg_start + Rout%seg_lengths(proc,nseg)-1
!	   do VectIndex = seg_start,seg_end
!	      do AttrIndex = 1,numi
!		 aV%iAttr(AttrIndex,VectIndex)=ip2(proc)%pi(j)
!		 j=j+1
!	      enddo
!	   enddo
!	enddo
!     endif
!     ! Load the real data
!     if(numr .ge. 1) then
!	call MPI_WAITANY(Rout%nprocs,rreqs,proc,rstatus,ier)
!	if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITANY(reals)',ier)
!	k=1
!	! load the correct pieces of the real vectors
!	do nseg = 1,Rout%num_segs(proc)
!	   seg_start = Rout%seg_starts(proc,nseg)
!	   seg_end = seg_start + Rout%seg_lengths(proc,nseg)-1
!	   do VectIndex = seg_start,seg_end
!	      do AttrIndex = 1,numr
!		 aV%rAttr(AttrIndex,VectIndex)=rp2(proc)%pr(k)
!		 k=k+1
!	      enddo
!	   enddo
!	enddo
!     endif
!  enddo
!........................................................................

  ! Deallocate all structures
  if(numi .ge. 1) then

     ! Deallocate the receive buffers
     do proc=1,Rout%nprocs
	deallocate(ip2(proc)%pi,stat=ier)
	if(ier/=0) call die(myname_,'deallocate(ip2%pi)',ier)
     enddo

     deallocate(ip2,stat=ier)
     if(ier/=0) call die(myname_,'deallocate(ip2)',ier)

     ! deallocate MPI request array
     deallocate(ireqs,stat=ier)
     if(ier/=0) call die(myname_,'deallocate(ireqs)',ier)

     ! deallocatAe MPI status array
     deallocate(istatus,stat=ier)
     if(ier/=0) call die(myname_,'deallocate(istatus)',ier)

  endif

  if(numr .ge. 1) then

     ! Deallocate the receive buffers
     do proc=1,Rout%nprocs
	deallocate(rp2(proc)%pr,stat=ier)
	if(ier/=0) call die(myname_,'deallocate(rp2%pr)',ier)
     enddo

     deallocate(rp2,stat=ier)
     if(ier/=0) call die(myname_,'deallocate(rp2)',ier)

     ! deallocate MPI request array
     deallocate(rreqs,stat=ier)
     if(ier/=0) call die(myname_,'deallocate(rreqs)',ier)

     ! deallocatAe MPI status array
     deallocate(rstatus,stat=ier)
     if(ier/=0) call die(myname_,'deallocate(rstatus)',ier)

  endif


end subroutine MCT_Recv




