!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$ 
!BOP -------------------------------------------------------------------
!
! !IROUTINE: MCT_Send - Distributed send of an Attribute Vector
!
! !DESCRIPTION:
! Send the the data in the {\tt AttrVect} {\tt aV} to the 
! component specified in the {\tt Router} {\tt Rout}.  An error will 
! result if the size of the attribute vector does not match the size
! parameter stored in the {\tt Router}.
!
! Requires a corresponding {\tt MCT\_Recv} to be called on the other component.
!
! {\bf N.B.:} The {\tt AttrVect} argument in the corresponding
! {\tt MCT\_Recv} call is assumed to have exactly the same attributes
! in exactly the same order as {\tt aV}.
!
! !INTERFACE:

 subroutine MCT_Send(aV, Rout)

!
! !USES:
!
      use m_MCTWorld, only : MCTWorld
      use m_MCTWorld, only : ThisMCTWorld
      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : nIAttr,nRAttr
      use m_AttrVect, only : lsize
      use m_Router,   only : Router
      use m_List,     only : List

      use m_mpif90
      use m_die
      use m_stdio

      implicit none

! !INPUT PARAMETERS:
!

      Type(AttrVect),       intent(in) :: aV	
      Type(Router),         intent(in) :: Rout

! !REVISION HISTORY:
! 07Feb01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
! 08Feb01 - R. Jacob <jacob@mcs.anl.gov> - First working code
! 18May01 - R. Jacob <jacob@mcs.anl.gov> - use MP_Type to
!           determine type in mpi_send
! 07Jun01 - R. Jacob <jacob@mcs.anl.gov> - remove logic to
!           check "direction" of Router.  remove references
!           to ThisMCTWorld%mylrank
! 03Aug01 - E. Ong <eong@mcs.anl.gov> - Explicitly specify the starting
!           address in mpi_send.  
! 15Feb02 - R. Jacob <jacob@mcs.anl.gov> - Use MCT_comm
! 26Mar02 - E. Ong <eong@mcs.anl.gov> - Apply faster copy order
! 26Sep02 - R. Jacob <jacob@mcs.anl.gov> - Check Av against Router lAvsize
! 05Nov02 - R. Jacob <jacob@mcs.anl.gov> - Remove iList, rList arguments.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_='MCT_Send'
  integer ::	numi,numr,i,j,k,ier
  integer ::    mycomp,othercomp
  integer ::    AttrIndex,VectIndex,seg_start,seg_end
  integer ::    proc,nseg,tag
  integer ::    mp_Type_rp1
  integer, dimension(:), pointer	     :: ireqs,rreqs
  integer, dimension(:,:), allocatable       :: istatus,rstatus

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

!check Av size against Router
!
  if(lsize(aV) /= Rout%lAvsize) then
    write(stderr,'(2a)') myname_, &
    ' MCTERROR:  AV size not appropriate for this Router...exiting'
    call die(myname_)
  endif

! get ids of components involved in this communication
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
  if(ier/=0) call die(myname_,'allocate(ip1)',ier)

! allocate buffers to hold all outgoing data
  do proc=1,Rout%nprocs
     allocate(ip1(proc)%pi(Rout%locsize(proc)*numi),stat=ier)
     if(ier/=0) call die(myname_,'allocate(ip1%pi)',ier)
  enddo

! allocate MPI request array
  allocate(ireqs(Rout%nprocs),stat=ier)
  if(ier/=0) call die(myname_,'allocate(ireqs)',ier)

! allocate MPI status array
  allocate(istatus(MP_STATUS_SIZE,Rout%nprocs),stat=ier)
  if(ier/=0) call die(myname_,'allocate(istatus)',ier)

  endif

!!!!!!!!!!!!!! IF SENDING REAL DATA
  if(numr .ge. 1) then

! allocate the number of pointers needed
  allocate(rp1(Rout%nprocs),stat=ier)
  if(ier/=0) call die(myname_,'allocate(rp1)',ier)

! allocate buffers to hold all outgoing data
  do proc=1,Rout%nprocs
     allocate(rp1(proc)%pr(Rout%locsize(proc)*numr),stat=ier)
     if(ier/=0) call die(myname_,'allocate(rp1%pr)',ier)
  enddo

! allocate MPI request array
  allocate(rreqs(Rout%nprocs),stat=ier)
  if(ier/=0) call die(myname_,'allocate(rreqs)',ier)

! allocate MPI status array
  allocate(rstatus(MP_STATUS_SIZE,Rout%nprocs),stat=ier)
  if(ier/=0) call die(myname_,'allocate(rstatus)',ier)

  mp_Type_rp1=MP_Type(rp1(1)%pr(1))

  endif


  ! Load data going to each processor
  do proc = 1,Rout%nprocs
    
     if( Rout%num_segs(proc) > 1 ) then

	j=1
	k=1

	! load the correct pieces of the integer and real vectors
	do nseg = 1,Rout%num_segs(proc)
	   seg_start = Rout%seg_starts(proc,nseg)
	   seg_end = seg_start + Rout%seg_lengths(proc,nseg)-1
	   do VectIndex = seg_start,seg_end
	      do AttrIndex = 1,numi
		 ip1(proc)%pi(j) = aV%iAttr(AttrIndex,VectIndex)
		 j=j+1
	      enddo
	      do AttrIndex = 1,numr
		 rp1(proc)%pr(k) = aV%rAttr(AttrIndex,VectIndex)
		 k=k+1
	      enddo
	   enddo
	enddo

     endif

     ! Send the integer data
     if(numi .ge. 1) then

	! corresponding tag logic must be in MCT_Recv
	tag = 100000*mycomp + 1000*ThisMCTWorld%mygrank + &
	      500 + Rout%pe_list(proc)

	if( Rout%num_segs(proc) > 1 ) then

	   call MPI_ISEND(ip1(proc)%pi(1), &
		Rout%locsize(proc)*numi,MP_INTEGER,Rout%pe_list(proc), &
		tag,ThisMCTWorld%MCT_comm,ireqs(proc),ier)
	   
	else

	   call MPI_ISEND(aV%iAttr(1,Rout%seg_starts(proc,1)), & 
		Rout%locsize(proc)*numi,MP_INTEGER,Rout%pe_list(proc), &
		tag,ThisMCTWorld%MCT_comm,ireqs(proc),ier)

	endif

	if(ier /= 0) call MP_perr_die(myname_,'MPI_ISEND(ints)',ier)

     endif

     ! Send the real data
     if(numr .ge. 1) then

       ! corresponding tag logic must be in MCT_Recv
       tag = 100000*mycomp + 1000*ThisMCTWorld%mygrank + &
	     700 + Rout%pe_list(proc)

       if( Rout%num_segs(proc) > 1 ) then

	  call MPI_ISEND(rp1(proc)%pr(1), &
	       Rout%locsize(proc)*numr,mp_Type_rp1,Rout%pe_list(proc), &
	       tag,ThisMCTWorld%MCT_comm,rreqs(proc),ier)

       else

	  call MPI_ISEND(aV%rAttr(1,Rout%seg_starts(proc,1)), &
	       Rout%locsize(proc)*numr,mp_Type_rp1,Rout%pe_list(proc), &
	       tag,ThisMCTWorld%MCT_comm,rreqs(proc),ier)

       endif

       if(ier /= 0) call MP_perr_die(myname_,'MPI_ISEND(reals)',ier)

    endif

  enddo


  ! wait for all sends to complete
  if(numi .ge. 1) then

     call MPI_WAITALL(Rout%nprocs,ireqs,istatus,ier)
     if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL(ints)',ier)

     do proc=1,Rout%nprocs
	deallocate(ip1(proc)%pi,stat=ier)
	if(ier/=0) call die(myname_,'deallocate(ip1%pi)',ier)
     enddo

     deallocate(ip1,stat=ier)
     if(ier/=0) call die(myname_,'deallocate(ip1)',ier)
       
     deallocate(ireqs,stat=ier)
     if(ier/=0) call die(myname_,'deallocate(ireqs)',ier)

     deallocate(istatus,stat=ier)
     if(ier/=0) call die(myname_,'deallocate(istatus)',ier)

  endif

  if(numr .ge. 1) then

     call MPI_WAITALL(Rout%nprocs,rreqs,rstatus,ier)
     if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL(reals)',ier)

     do proc=1,Rout%nprocs
	deallocate(rp1(proc)%pr,stat=ier)
	if(ier/=0) call die(myname_,'deallocate(rp1%pi)',ier)
     enddo

     deallocate(rp1,stat=ier)
     if(ier/=0) call die(myname_,'deallocate(rp1)',ier)

     deallocate(rreqs,stat=ier)
     if(ier/=0) call die(myname_,'deallocate(rreqs)',ier)

     deallocate(rstatus,stat=ier)
     if(ier/=0) call die(myname_,'deallocate(rstatus)',ier)

  endif


end subroutine MCT_Send


