!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: MCT_Send
!
! !DESCRIPTION:
! Send the local AttrVect to another component using a Router.
! Will send entire AttrVect or only the integer and real attributes
! specified in iList and rList. (NOTE:  the latter is not yet supported)
! Requires a corresponding MCT_Recv to be called on the other component.
!
! !INTERFACE:

 subroutine MCT_Send(aV, Rout,iList,rList)
!
! !USES:
!
      use m_Router,only  : Router
      use m_AttrVect,only : AttrVect
      use m_AttrVect,only : nIAttr,nRAttr
      use m_list,only:	List
      use m_list,only:  nitem
      use m_mpif90
      use m_die

      implicit none
      Type(AttrVect),intent(in) :: 	aV	
      Type(Router),intent(in) ::	Rout
      Type(List),optional,intent(in) ::	iList
      Type(List),optional,intent(in) ::	rList

! !REVISION HISTORY:
!      07Feb01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_='MCT_Send'
  integer ::	numi,numr,i,proc,buffsize,j,k
  integer ::	ier,nseg
  integer,dimension(:),pointer	:: ireqs,rreqs
  integer,dimension(:),pointer	:: Ibuffer
  real,dimension(:),pointer	:: Rbuffer

!  find total number of real and integer vectors
  numi = nIAttr(aV)
  numr = nRAttr(aV)

! Send data to each processor
  allocate(ireqs(Rout%nprocs),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(ireqs)',ier)
  allocate(rreqs(Rout%nprocs),stat=ier)
  if(ier/=0) call MP_perr_die(myname_,'allocate(rreqs)',ier)

  do proc=1,Rout%nprocs

!  find out how much data to send to this processor
    buffsize=0
    do i=1,Rout%num_segs(proc)
      buffsize=buffsize+Rout%seg_lengths(proc,i)
    enddo

! Send the integer data
    if(numi .ge. 1) then
      allocate(Ibuffer(buffsize*numi),stat=ier)
      if(ier/=0) call MP_perr_die(myname_,'allocate(mystarts,..)',ier)

      k=1
! load all integer vectors to be sent into buffer
      do j=1,numi

! load the correct pieces of this particular integer vector
       do nseg=1,Rout%num_segs(proc)
	 do i=0,Rout%seg_lengths(proc,nseg)-1
           Ibuffer(k)=aV%iAttr(j,Rout%seg_starts(proc,nseg)+i)
	   k=k+1
	 enddo
       enddo

      enddo

! send it
      call MPI_ISEND(Ibuffer,buffsize,MP_INTEGER,Rout%pe_list(proc),&
	   1000,MP_COMM_WORLD,ireqs(proc))
    endif


! Send the real data
    if(numr .ge. 1) then
      allocate(Rbuffer(buffsize*numr),stat=ier)
      if(ier/=0) call MP_perr_die(myname_,'allocate(mystarts,..)',ier)

      k=1
! load all integer vectors to be sent into buffer
      do j=1,numr

! load the correct pieces of this particular integer vector
       do nseg=1,Rout%num_segs(proc)
	 do i=0,Rout%seg_lengths(proc,nseg)-1
           Rbuffer(k)=aV%rAttr(j,Rout%seg_starts(proc,nseg)+i)
	   k=k+1
	 enddo
       enddo

      enddo

! send it
      call MPI_ISEND(Rbuffer,buffsize,MP_REAL,Rout%pe_list(proc),&
	   1000,MP_COMM_WORLD,rreqs(proc))
   endif

  enddo

end subroutine
