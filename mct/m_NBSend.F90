!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$ 
!BOP -------------------------------------------------------------------
!
! !MODULE: m_NBSend -- Module for non-blocking version of send_ in m_Transfer
!
! !DESCRIPTION:
! This module provides functions and data types to support a non-blocking
! version of {\tt send\_}.
!
! {\bf N.B.:} This module may be deleted in future versions.
!
! !SEE ALSO:
! m_Transfer
!
! !INTERFACE:
 module m_NBSend
!
! !USES:

      implicit none

      private	! except

! !PUBLIC TYPES:

      public :: MCTReqs     ! the reqs datatype

      type MCTReqs
        integer,dimension(:),pointer :: ireqs   ! the integer sends
        integer,dimension(:),pointer :: rreqs   ! the real sends
      end type MCTReqs

! !PUBLIC MEMBER FUNCTIONS:

      public :: isend   ! the non blocking send_

      public :: wait    ! Wait for the nonblocking send to finish

      interface isend ; module procedure isend_ ; end interface
      interface wait ; module procedure wait_ ; end interface

! !DEFINED PARAMETERS:
!
  integer,parameter                    :: DefaultTag = 600


! !REVISION HISTORY:
! 07Jun01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
! 14Feb02 - R. Jacob <jacob@mcs.anl.gov> - Use MCT_comm instead
!           of MP_COMM_WORLD
!EOP ___________________________________________________________________

      ! declare a pointer structure for the real data
      type :: rptr
	 real,dimension(:),pointer :: pr
      end type rptr

      ! declare a pointer structure for the integer data
      type :: iptr
	 integer,dimension(:),pointer :: pi
      end type iptr

      ! declare arrays of pointers to hold data to be sent
      type(rptr),dimension(:),allocatable :: rp1
      type(iptr),dimension(:),allocatable :: ip1

  character(len=*),parameter :: myname='m_NBSend'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: isend_
!
! !DESCRIPTION:
! Send the the data in the {\tt AttrVect} {\tt aV} to the 
! component specified in the {\tt Router} {\tt Rout}.  An error will 
! result if the size of the attribute vector does not match the size
! parameter stored in the {\tt Router}.
!
! Returns immediately after posting the send with a {\tt MCT\_Req} data type
! {\tt Reqs}
!
! Requires a corresponding {\tt recv\_} to be called on the other component.
!
! The optional argument {\tt Tag} can be used to set the tag value used in
! the data transfer.  DefaultTag will be used otherwise. {\tt Tag} must be
! the same in the matching {\tt recv\_}
!
!
! {\bf N.B.:} The {\tt AttrVect} argument in the corresponding
! {\tt recv\_} call is assumed to have exactly the same attributes
! in exactly the same order as {\tt aV}.
!
! {\bf N.B.:} Currently, only one instance of ISend\_ can be outstanding at a time
! within an application.  This is because the buffers holding the data are currently a
! private data member in this module.
!
! !INTERFACE:

 subroutine isend_(aV, Rout, Reqs, Tag)

!
! !USES:
!
      use m_Router,only   : Router
      use m_AttrVect,only : AttrVect
      use m_AttrVect,only : nIAttr,nRAttr
      use m_AttrVect,only : lsize
      use m_MCTWorld,only : MCTWorld
      use m_MCTWorld,only : ThisMCTWorld

      use m_mpif90
      use m_die
      use m_stdio

      implicit none

! !INPUT PARAMETERS:

      Type(AttrVect),    intent(in) :: 	aV     ! the Attribute vector to send	
      Type(Router),      intent(in) ::	Rout   ! the router to use
      integer,optional,  intent(in) ::	Tag    ! message tag

! !OUTPUT PARAMETERS: 

      Type(MCTReqs),intent(inout) ::	Reqs   ! the returned list of MPI requests

! !REVISION HISTORY:
! 07Jun01 - R. Jacob <jacob@mcs.anl.gov> - first prototype
! 28Mar01 - E. Ong <eong@mcs.anl.gov> - changed copy order to correspond
!           to MCT_Recv
! 06Nov02 - R. Jacob <jacob@mcs.anl.gov> - Remove iList and rList arguments.
!           Add check with Router lsize.
! 11Nov02 - R. Jacob <jacob@mcs.anl.gov> - Add Tag optional argument
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'isend_'
  integer ::	numi,numr,i,j,k,ier
  integer ::    mycomp,othercomp
  integer ::    AttrIndex,VectIndex,seg_start,seg_end
  integer ::    proc,nseg,mytag
  integer ::    mp_Type_rp1

!--------------------------------------------------------

!check Av size against Router
!
  if(lsize(aV) /= Rout%lAvsize) then
    write(stderr,'(2a)') myname_, &
    ' MCTERROR:  AV size not appropriate for this Router...exiting'
    call die(myname_)
  endif


  mycomp=Rout%comp1id
  othercomp=Rout%comp2id

!  find total number of real and integer vectors
! for now, assume we are sending all of them
  numi = nIAttr(aV)
  numr = nRAttr(aV)

! Nullify the pointers
  nullify(Reqs%ireqs)
  nullify(Reqs%rreqs)

!!!!!!!!!!!!!! IF SENDING INTEGER DATA
  if(numi .ge. 1) then

! allocate the number of pointers needed
  allocate(ip1(Rout%nprocs),stat=ier)
  if(ier/=0) call die(myname_,'allocate(ip1)',ier)

! allocate MPI request array
  allocate(Reqs%ireqs(Rout%nprocs),stat=ier)
  if(ier/=0) call die(myname_,'allocate(Reqs%ireqs)',ier)
! write(*,*)"Isend",Rout%nprocs,size(Reqs%ireqs)

! allocate buffers to hold all outgoing data
   do proc=1,Rout%nprocs
    allocate(ip1(proc)%pi(Rout%locsize(proc)*numi),stat=ier)
    if(ier/=0) call die(myname_,'allocate(ip1%pi)',ier)
   enddo
  endif

!!!!!!!!!!!!!! IF SENDING REAL DATA
  if(numr .ge. 1) then

! allocate the number of pointers needed
  allocate(rp1(Rout%nprocs),stat=ier)
  if(ier/=0) call die(myname_,'allocate(rp1)',ier)

! allocate MPI request array
  allocate(Reqs%rreqs(Rout%nprocs),stat=ier)
  if(ier/=0) call die(myname_,'allocate(Reqs%rreqs)',ier)

! allocate buffers to hold all outgoing data
   do proc=1,Rout%nprocs
    allocate(rp1(proc)%pr(Rout%locsize(proc) *numr),stat=ier)
    if(ier/=0) call die(myname_,'allocate(rp1%pr)',ier)
   enddo

! define the real type
   mp_Type_rp1=MP_Type(rp1(1)%pr(1))

  endif

  ! Load data going to each processor
  do proc = 1,Rout%nprocs
    
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

     ! Send the integer data
     if(numi .ge. 1) then

        ! set tag
        mytag = DefaultTag
        if(present(Tag)) mytag=Tag

	call MPI_ISEND(ip1(proc)%pi(1),Rout%locsize(proc)*numi,MP_INTEGER,&
	     Rout%pe_list(proc),mytag,ThisMCTWorld%MCT_comm,Reqs%ireqs(proc),ier)

	if(ier /= 0) call MP_perr_die(myname_,'MPI_ISEND(ints)',ier)

     endif

     ! Send the real data
     if(numr .ge. 1) then

       ! set tag
       mytag = DefaultTag + 1
       if(present(Tag)) mytag=Tag +1

       call MPI_ISEND(rp1(proc)%pr(1),Rout%locsize(proc)*numr,mp_Type_rp1,&
	    Rout%pe_list(proc),mytag,ThisMCTWorld%MCT_comm,Reqs%rreqs(proc),ier)

       if(ier /= 0) call MP_perr_die(myname_,'MPI_ISEND(reals)',ier)

    endif

  enddo


end subroutine isend_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: wait_  Wait for an MCT_ISend to complete.
!
! !DESCRIPTION:
! Wait for the {\tt isend\_} on the {\tt Router} {\tt Rout} and handled by
! {\tt Reqs} to complete.  Deallocate internal memory buffers and memory in {\tt Reqs}
! when message has been received.
!
! !INTERFACE:

 subroutine wait_(Rout,Reqs)

!
! !USES:
!
      use m_Router,only  : Router
      use m_mpif90
      use m_die
      use m_stdio


      implicit none

! !INPUT PARAMETERS:

      Type(Router),intent(in) ::        Rout

! !INPUT/OUTPUT PARAMETERS:

      Type(MCTReqs),intent(inout) ::	Reqs

! !REVISION HISTORY:
! 07Jun01 - R. Jacob <jacob@mcs.anl.gov> - first prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'wait_'
  integer ::	ier
  integer :: proc
  integer,dimension(:,:),allocatable	:: istatus,rstatus

!--------------------------------------------------------

  if(associated(Reqs%ireqs)) then

     ! allocatAe MPI status array
     allocate(istatus(MP_STATUS_SIZE,Rout%nprocs),stat=ier)
     if(ier/=0) call die(myname_,'allocate(istatus)',ier)

     !  write(*,*)"wait",Rout%nprocs,size(Reqs%ireqs)
     call MPI_WAITALL(Rout%nprocs,Reqs%ireqs,istatus,ier)
     if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL(ints)',ier)

     ! deallocatAe MPI status array
     deallocate(istatus,stat=ier)
     if(ier/=0) call die(myname_,'deallocate(istatus)',ier)

     ! done waiting, free up ireqs
     deallocate(Reqs%ireqs,stat=ier)
     if(ier /= 0) call die(myname_,'deallocate(Reqs%ireqs)',ier)

    ! Deallocate the send buffers
     do proc=1,Rout%nprocs
	deallocate(ip1(proc)%pi,stat=ier)
	if(ier/=0) call die(myname_,'deallocate(ip1%pr)',ier)
     enddo

     deallocate(ip1,stat=ier)
     if(ier/=0) call die(myname_,'deallocate(ip1)',ier)

  endif

  if(associated(Reqs%rreqs)) then

     ! allocate MPI status array
     allocate(rstatus(MP_STATUS_SIZE,Rout%nprocs),stat=ier)
     if(ier/=0) call die(myname_,'allocate(rstatus)',ier)

     call MPI_WAITALL(Rout%nprocs,Reqs%rreqs,rstatus,ier)
     if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL(reals)',ier)

     ! deallocate MPI status array
     deallocate(rstatus,stat=ier)
     if(ier/=0) call die(myname_,'deallocate(rstatus)',ier)

     ! done waiting, free up Reqs%rreqs
     deallocate(Reqs%rreqs,stat=ier)
     if(ier /= 0) call die(myname_,'deallocate(Reqs%rreqs)',ier)

    ! Deallocate the send buffers
     do proc=1,Rout%nprocs
	deallocate(rp1(proc)%pr,stat=ier)
	if(ier/=0) call die(myname_,'deallocate(rp1%pr)',ier)
     enddo

     deallocate(rp1,stat=ier)
     if(ier/=0) call die(myname_,'deallocate(rp1)',ier)

  endif


end subroutine wait_

end module m_NBSend
