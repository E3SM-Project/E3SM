!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$ 
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Rearranger -- Remaps an AttrVect within a group of processes
!
! !DESCRIPTION:
! This module provides routines and datatypes for rearranging data
! between two {\tt Attribute Vectors} defined on the same grid but
! with two different {\tt GlobalSegMaps}.  ''Rearrange'' is a
! generalized form of a parallel matrix transpose.
! A parallel matrix transpose can take advantage of symmetry in the
! data movement algorithm.  An MCT Rearranger makes no assumptions
! about symmetry.
!
! When data needs to move between two components and the components
! share any processors, use m\_Rearranger.  If the components are on
! distinct sets of processors, use m\_Transfer.
!
! !SEE ALSO:
!  m_Transfer
! 
!
! !INTERFACE:

 module m_Rearranger

!
! !USES:

      use m_Router, only : Router

      implicit none

      private	! except

! !PUBLIC DATA MEMBERS:

      public :: Rearranger  ! The class data structure

      type :: Rearranger
         private
         type(Router) :: SendRouter
         type(Router) :: RecvRouter
         integer,dimension(:,:),pointer :: LocalPack
         integer :: LocalSize
      end type Rearranger

! !PUBLIC MEMBER FUNCTIONS:

      public :: init         ! creation method

      public :: rearrange    ! the rearrange routine

      public :: clean        ! destruction method

      interface init      ; module procedure init_      ; end interface
      interface Rearrange ; module procedure Rearrange_ ; end interface
      interface clean     ; module procedure clean_     ; end interface

! !DEFINED PARAMETERS:

  integer,parameter                    :: DefaultTag = 500


! !REVISION HISTORY:
! 31Jan02 - E.T. Ong <eong@mcs.anl.gov> - initial prototype
! 04Jun02 - E.T. Ong <eong@mcs.anl.gov> - changed local copy structure to
!           LocalSize. Made myPid a global process in MCTWorld.
! 27Sep02 - R. Jacob <jacob@mcs.anl.gov> - Remove SrcAVsize and TrgAVsize
!           and use Router%lAvsize instead for sanity check.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='MCT::m_Rearranger'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: Init_ - Initialize a Rearranger
!
! !DESCRIPTION:
! This routine takes two {\tt GlobalSegMap} inputs, {\tt SourceGSMap}
! and {\tt TargetGSMap} and build a Rearranger {\tt OutRearranger}
! between them. {\tt myComm} is used for the internal communication.
!
! !INTERFACE:

 subroutine init_(SourceGSMap,TargetGSMap,myComm,OutRearranger)

!
! !USES:
!

   use m_MCTWorld,     only : ThisMCTWorld
   use m_GlobalSegMap, only : GlobalSegMap
   use m_GlobalSegMap, only : GSMap_lsize => lsize
   use m_Router,       only : Router     
   use m_Router,       only : Router_init => init
   use m_mpif90
   use m_die
   use m_stdio

   implicit none
  
! !INPUT PARAMETERS:
!
   type(GlobalSegMap), intent(in)            :: SourceGSMap, TargetGSMap
   integer,            intent(in)            :: myComm

! !OUTPUT PARAMETERS:
!
   type(Rearranger),   intent(out)           :: OutRearranger

! !REVISION HISTORY:
! 31Jan02 - E.T. Ong <eong@mcs.anl.gov> - initial prototype
! 20Mar02 - E.T. Ong <eong@mcs.anl.gov> - working code
! 05Jun02 - E.T. Ong <eong@mcs.anl.gov> - Use LocalPack
!EOP ___________________________________________________________________
   character(len=*),parameter :: myname_=myname//'::init_'
   integer,dimension(:,:),pointer :: temp_seg_starts,temp_seg_lengths
   integer,dimension(:),pointer :: temp_pe_list,temp_numsegs,temp_locsize
   integer :: temp_maxsize,temp_nprocs,maxsegcount
   integer :: procindex,nprocs,nseg,len,myPid
   integer :: src_seg_start,src_seg_length,trg_seg_start,trg_seg_length
   integer :: i,j,k,l,m,n,ier
   logical :: SendingToMyself,ReceivingFromMyself

   ! Initialize Router component of Rearranger
   call Router_init(SourceGSMap,TargetGSMap,myComm,OutRearranger%SendRouter)
   call Router_init(TargetGSMap,SourceGSMap,myComm,OutRearranger%RecvRouter)

   ! SANITY CHECK: Make sure that if SendRouter is sending to self, then,
   ! by definition, RecvRouter is also receiving from self. If this is not 
   ! true, then write to stderr and die. 

   call MP_comm_rank(ThisMCTWorld%MCT_comm,myPid,ier)
   if(ier/=0) call MP_perr_die(myname_,'MP_comm_rank',ier)

   SendingToMyself = .false.
   ReceivingFromMyself = .false.

   do i=1,OutRearranger%SendRouter%nprocs
      if(OutRearranger%SendRouter%pe_list(i) == myPid) then
	 SendingToMyself = .true.
      endif
   enddo

   do i=1,OutRearranger%RecvRouter%nprocs
      if(OutRearranger%RecvRouter%pe_list(i) == myPid) then
	 ReceivingFromMyself = .true.
      endif
   enddo

   if( SendingToMyself.or.ReceivingFromMyself ) then
      if( .not. (SendingToMyself.and.ReceivingFromMyself) ) then
	 call die(myname_,"SendRouter is not compatible with RecvRouter")
      endif
   endif


   ! If not sending to nor receiving from own processor then initialize 
   ! the rearranger so that no local copy can be made. Then end the routine. 

   if( .not. (SendingToMyself.or.ReceivingFromMyself) ) then
      nullify(OutRearranger%LocalPack)
      allocate(OutRearranger%LocalPack(0,0),stat=ier)
      if(ier/=0) call die(myname_,'allocate(OutRearranger%LocalPack(0,0))',ier)
      OutRearranger%LocalSize=0
   endif


   ! Start the process of Router modification: Router information for 
   ! the local processor is extracted out and put into the local copy 
   ! structure- Rearranger%LocalPack. Router structures are then reassigned 
   ! to exclude the local copy information.


   ! Operate on SendRouter and create local copy structures.

   if( SendingToMyself.and.ReceivingFromMyself ) then

   temp_nprocs = OutRearranger%SendRouter%nprocs-1
   maxsegcount = SIZE(OutRearranger%SendRouter%seg_starts,2)

   ! Allocate temporary Router structures to be used for modifying SendRouter 
   nullify(temp_seg_starts,temp_seg_lengths,temp_pe_list, &
           temp_numsegs,temp_locsize)
   allocate(temp_seg_starts(temp_nprocs,maxsegcount), &
            temp_seg_lengths(temp_nprocs,maxsegcount), &
	    temp_pe_list(temp_nprocs), &
            temp_numsegs(temp_nprocs), &
            temp_locsize(temp_nprocs), stat=ier)
   if(ier/=0) call die(myname_,'allocate(temp_seg_starts...)',ier)

   temp_maxsize=0
   procindex=0
   nullify(OutRearranger%LocalPack)

   ! Start assigning Rearranger copy structures and  
   ! non-local Router components
   do i=1,OutRearranger%SendRouter%nprocs

      ! Gather local copy information 
      if(OutRearranger%SendRouter%pe_list(i) == myPid) then

	 ! Allocate Rearranger copy structure
	 allocate(OutRearranger%LocalPack(2, &
                  OutRearranger%SendRouter%locsize(i)),stat=ier)
	 if(ier/=0) call die(myname_,'allocate(OutRearranger%LocalPack)',ier)
	 OutRearranger%LocalPack = 0

	 m=0
	 do nseg = 1,OutRearranger%SendRouter%num_segs(i)
	    src_seg_start = OutRearranger%SendRouter%seg_starts(i,nseg)
	    src_seg_length = OutRearranger%SendRouter%seg_lengths(i,nseg)-1
	    do len=0,src_seg_length
	       m=m+1
	       OutRearranger%LocalPack(2,m) = src_seg_start+len
	    enddo
	 enddo

      else

	 ! Gather non-local Router information
	 procindex = procindex+1
	 temp_seg_starts(procindex,1:maxsegcount) = &
            OutRearranger%SendRouter%seg_starts(i,1:maxsegcount)
	 temp_seg_lengths(procindex,1:maxsegcount) = &
            OutRearranger%SendRouter%seg_lengths(i,1:maxsegcount)
	 temp_pe_list(procindex) = OutRearranger%SendRouter%pe_list(i)
	 temp_numsegs(procindex) = OutRearranger%SendRouter%num_segs(i)
	 temp_locsize(procindex) = OutRearranger%SendRouter%locsize(i)
	 temp_maxsize = max(temp_locsize(procindex),temp_maxsize)

      endif

   enddo

   ! Copy SendRouter components back in

   ! Deallocate existing SendRouter components
   deallocate(OutRearranger%SendRouter%seg_starts,&
              OutRearranger%SendRouter%seg_lengths, &
	      OutRearranger%SendRouter%pe_list, &
              OutRearranger%SendRouter%num_segs, &
              OutRearranger%SendRouter%locsize,stat=ier)
   if(ier/=0) call die(myname_, &
                  'deallocate(OutRearranger%SendRouter%seg_starts...)',ier)

   ! Re-allocate SendRouter components
   allocate(OutRearranger%SendRouter%seg_starts(temp_nprocs,maxsegcount), &
            OutRearranger%SendRouter%seg_lengths(temp_nprocs,maxsegcount), &
	    OutRearranger%SendRouter%pe_list(temp_nprocs), &
	    OutRearranger%SendRouter%num_segs(temp_nprocs), &
            OutRearranger%SendRouter%locsize(temp_nprocs),stat=ier)
   if(ier/=0) call die(myname_, &
                   'allocate(OutRearranger%SendRouter%seg_starts...)',ier)      

   ! Copy back in the spliced router information
   OutRearranger%SendRouter%nprocs = temp_nprocs
   OutRearranger%SendRouter%seg_starts(1:temp_nprocs,1:maxsegcount) = &
      temp_seg_starts(1:temp_nprocs,1:maxsegcount)
   OutRearranger%SendRouter%seg_lengths(1:temp_nprocs,1:maxsegcount) = &
      temp_seg_lengths(1:temp_nprocs,1:maxsegcount)
   OutRearranger%SendRouter%pe_list(1:temp_nprocs) = &
      temp_pe_list(1:temp_nprocs)
   OutRearranger%SendRouter%num_segs(1:temp_nprocs) = &
      temp_numsegs(1:temp_nprocs)
   OutRearranger%SendRouter%locsize(1:temp_nprocs) = &
      temp_locsize(1:temp_nprocs)
   OutRearranger%SendRouter%maxsize = temp_maxsize

   deallocate(temp_seg_starts,temp_seg_lengths,temp_pe_list, &
              temp_numsegs,temp_locsize,stat=ier)
   if(ier/=0) call die(myname_,'deallocate(temp_seg_starts...)',ier)      


   ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::


   ! Operate on RecvRouter and create local copy structures.

   temp_nprocs = OutRearranger%RecvRouter%nprocs-1
   maxsegcount = SIZE(OutRearranger%RecvRouter%seg_starts,2)

  ! Allocate temporary Router structures to be used for modifying RecvRouter
   nullify(temp_seg_starts,temp_seg_lengths,temp_pe_list, &
           temp_numsegs,temp_locsize)
   allocate(temp_seg_starts(temp_nprocs,maxsegcount), &
            temp_seg_lengths(temp_nprocs,maxsegcount), &
	    temp_pe_list(temp_nprocs),temp_numsegs(temp_nprocs), &
            temp_locsize(temp_nprocs),stat=ier)
   if(ier/=0) call die(myname_,'allocate(temp_seg_starts...)',ier)

   temp_maxsize=0
   procindex = 0

   ! Start assigning Rearranger copy structures and  
   ! non-local Router components
   do i=1,OutRearranger%RecvRouter%nprocs

      ! Gather local copy information 
      if(OutRearranger%RecvRouter%pe_list(i) == myPid) then

	 ! Senity Check for Router%locsize
	 if( (SIZE(OutRearranger%LocalPack,2) /= &
              OutRearranger%RecvRouter%locsize(i)) ) then
	    call die(myname_, &
                     'Router Error: Local RecvRouter%locsize(myPid) /=  &
                     & Local SendRouter%locsize(myPid)')
	 endif

	 OutRearranger%LocalSize = OutRearranger%RecvRouter%locsize(i)

	 m=0
	 do nseg = 1,OutRearranger%RecvRouter%num_segs(i)
	    trg_seg_start = OutRearranger%RecvRouter%seg_starts(i,nseg)
	    trg_seg_length = OutRearranger%RecvRouter%seg_lengths(i,nseg)-1
	    do len=0,trg_seg_length
	       m=m+1
	       OutRearranger%LocalPack(1,m) = trg_seg_start+len
	    enddo
	 enddo
	 
      else

	 ! Gather non-local Router information
	 procindex = procindex+1
	 temp_seg_starts(procindex,1:maxsegcount) = &
	    OutRearranger%RecvRouter%seg_starts(i,1:maxsegcount)
	 temp_seg_lengths(procindex,1:maxsegcount) = &
            OutRearranger%RecvRouter%seg_lengths(i,1:maxsegcount)
	 temp_pe_list(procindex) = OutRearranger%RecvRouter%pe_list(i)
	 temp_numsegs(procindex) = OutRearranger%RecvRouter%num_segs(i)
	 temp_locsize(procindex) = OutRearranger%RecvRouter%locsize(i)
	 temp_maxsize = max(temp_locsize(procindex),temp_maxsize)

      endif

   enddo

   ! Copy RecvRouter components back in

   ! Deallocate existing SendRouter components
   deallocate(OutRearranger%RecvRouter%seg_starts, &
              OutRearranger%RecvRouter%seg_lengths, &
	      OutRearranger%RecvRouter%pe_list, &
              OutRearranger%RecvRouter%num_segs, &
              OutRearranger%RecvRouter%locsize,stat=ier)
   if(ier/=0) call die(myname_, &
                   'deallocate(OutRearranger%RecvRouter%seg_starts...)',ier)

   ! Re-allocate RecvRouter components
   allocate(OutRearranger%RecvRouter%seg_starts(temp_nprocs,maxsegcount), &
            OutRearranger%RecvRouter%seg_lengths(temp_nprocs,maxsegcount), &
	    OutRearranger%RecvRouter%pe_list(temp_nprocs), &
	    OutRearranger%RecvRouter%num_segs(temp_nprocs), &
            OutRearranger%RecvRouter%locsize(temp_nprocs),stat=ier)
   if(ier/=0) call die(myname_, &
                   'allocate(OutRearranger%RecvRouter%seg_starts...)',ier)      

    ! Copy back in the spliced router information
   OutRearranger%RecvRouter%nprocs = temp_nprocs
   OutRearranger%RecvRouter%seg_starts(1:temp_nprocs,1:maxsegcount) =  &
      temp_seg_starts(1:temp_nprocs,1:maxsegcount)
   OutRearranger%RecvRouter%seg_lengths(1:temp_nprocs,1:maxsegcount) = &
      temp_seg_lengths(1:temp_nprocs,1:maxsegcount)
   OutRearranger%RecvRouter%pe_list(1:temp_nprocs) = &
      temp_pe_list(1:temp_nprocs)
   OutRearranger%RecvRouter%num_segs(1:temp_nprocs) = &
      temp_numsegs(1:temp_nprocs)
   OutRearranger%RecvRouter%locsize(1:temp_nprocs) = &
      temp_locsize(1:temp_nprocs)
   OutRearranger%RecvRouter%maxsize = temp_maxsize

   deallocate(temp_seg_starts,temp_seg_lengths,temp_pe_list, &
              temp_numsegs,temp_locsize,stat=ier)
   if(ier/=0) call die(myname_,'deallocate(temp_seg_starts...)',ier)
   
   endif

 end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - Clean a Rearranger
!
! !DESCRIPTION:
! This routine deallocates allocated memory associated with the 
! input/output {\tt Rearranger} argument {\tt ReArr}.  The success 
! (failure) of this operation is reported in the zero (nonzero) value of
! the optional output {\tt INTEGER} argument {\tt status}.
!
! !INTERFACE:

 subroutine clean_(ReArr, status)

!
! !USES:
!
   use m_Router,only : Router     
   use m_Router,only : Router_clean => clean
   use m_mpif90
   use m_die
   use m_stdio

   implicit none
  
! !INPUT/OUTPUT PARAMETERS:
!
   type(Rearranger),    intent(inout)           :: ReArr

! !OUTPUT PARAMETERS:
!
   integer, optional,   intent(out)             :: status
   
! !REVISION HISTORY:
! 31Jan02 - E.T. Ong <eong@mcs.anl.gov> - initial prototype
! 20Mar02 - E.T. Ong <eong@mcs.anl.gov> - working code
!EOP ___________________________________________________________________
   character(len=*),parameter :: myname_=myname//'::clean_'
   integer :: ier

        ! Set output status flag (if present) to zero, which assumes
        ! success.

   if(present(status)) status = 0

        ! Clean up send and receive Routers:

   call Router_clean(ReArr%SendRouter,ier)
   if(ier /= 0) then
      write(stderr,'(2a,i8)') myname_, &
	   ':: ERROR--Router_clean(ReArr%SendRouter) failed with ier=',ier
      if(present(status)) then
	 status = ier
	 return
      else
	 call die(myname_)
      endif
   endif

   call Router_clean(ReArr%RecvRouter,ier)
   if(ier /= 0) then
      write(stderr,'(2a,i8)') myname_, &
	   ':: ERROR--Router_clean(ReArr%RecvRouter) failed with ier=',ier
      if(present(status)) then
	 status = ier
	 return
      else
	 call die(myname_)
      endif
   endif

       ! Clean up Local on-PE copy buffer:

   if(associated(ReArr%LocalPack)) then
      deallocate(ReArr%LocalPack, stat=ier)
      if(ier /= 0) then
	 write(stderr,'(2a,i8)') myname_, &
	      ':: ERROR--deallocate(ReArr%LocalPack) failed with stat=',ier
	 if(present(status)) then
	    status=ier
	 else
	    call die(myname_)
	 endif
      endif
   endif

 end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rearrange_ - Rearrange data between two Attribute Vectors
!
! !DESCRIPTION: 
! This subroutine will take data in the {\tt SourceAv} Attribute
! Vector and rearrange it to match the GlobalSegMap used to define
! the {\tt TargetAv} Attribute Vector using the Rearrnger
! {\tt InRearranger}.
!
! The optional argument {\tt Tag} can be used to set the tag value used in 
! the rearrangement.  DefaultTag will be used otherwise.
!
! If the optional argument {\tt Sum} is present, data for the same
! physical point coming from two or more processes will be summed.
! Otherwise, data is overwritten.
!
! If the optional argument {\tt Vector} is present and true,
! vector architecture-friendly parts of this routine will be invoked.
!
! The size of the {\tt SourceAv} and {\tt TargetAv}
! argument must match those stored in the {\tt InRearranger} or
! and error will result.
!
! {\bf N.B.:} {\tt SourceAv} and {\tt TargetAv} are
! assumed to have exactly the same attributes
! in exactly the same order.
!
! !INTERFACE:

 subroutine rearrange_(SourceAV,TargetAV,InRearranger,Tag,Sum,Vector)

!
! !USES:
!

   use m_MCTWorld,only :MCTWorld
   use m_MCTWorld,only :ThisMCTWorld
   use m_AttrVect,  only : AttrVect
   use m_AttrVect,  only : AttrVect_init => init
   use m_AttrVect,  only : AttrVect_lsize => lsize
   use m_AttrVect,  only : AttrVect_zero => zero
   use m_AttrVect,  only : nIAttr,nRAttr
   use m_Router,    only : Router     
   use m_realkinds, only : FP
   use m_mpif90
   use m_die
   use m_stdio

   implicit none
  
! !INPUT/OUTPUT PARAMETERS:
!
   type(AttrVect),             intent(inout)   :: TargetAV
   
! !INPUT PARAMETERS:
!
   type(AttrVect),             intent(in)      :: SourceAV
   type(Rearranger), target,   intent(in)      :: InRearranger
   integer,          optional, intent(in)      :: Tag
   logical,          optional, intent(in)      :: Sum
   logical,          optional, intent(in)      :: Vector

! !REVISION HISTORY:
! 31Jan02 - E.T. Ong <eong@mcs.anl.gov> - initial prototype
! 20Mar02 - E.T. Ong <eong@mcs.anl.gov> - working code
! 08Jul02 - E.T. Ong <eong@mcs.anl.gov> - change intent of Target,Source
! 29Oct03 - R. Jacob <jacob@mcs.anl.gov> - add optional argument vector
!           to control use of vector-friendly mods provided by Fujitsu.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'Rearrange_'
  integer ::	numi,numr,i,j,k,ier
  integer ::    VectIndex,AttrIndex,seg_start,seg_end
  integer ::    localindex,SrcVectIndex,TrgVectIndex,IAttrIndex,RAttrIndex
  integer ::    proc,numprocs,nseg
  integer ::    mp_Type_rp
  integer ::    mytag
  logical ::    usevector
!-----------------------------------------------------------------------

   ! DECLARE STRUCTURES FOR MPI ARGUMENTS.

   ! declare a pointer structure for the real data
   type :: rptr
      real(FP),dimension(:),pointer :: pr
   end type rptr

   ! declare a pointer structure for the integer data
   type :: iptr
      integer,dimension(:),pointer :: pi
   end type iptr

   ! declare arrays of pointers to hold data to be sent
   type(rptr),dimension(:),allocatable :: rp1
   type(iptr),dimension(:),allocatable :: ip1

   ! declare arrays of pointers to hold data to be received
   type(rptr),dimension(:),allocatable :: rp2
   type(iptr),dimension(:),allocatable :: ip2

   ! Structure to hold MPI request information for sends
   integer,dimension(:),allocatable :: send_ireqs
   integer,dimension(:),allocatable :: send_rreqs

   ! Structure to hold MPI request information for sends
   integer,dimension(:),allocatable :: recv_ireqs
   integer,dimension(:),allocatable :: recv_rreqs

   ! Structure to hold MPI status information for sends 
   integer,dimension(:,:),allocatable :: send_istatus,send_rstatus

   ! Structure to hold MPI status information for sends 
   integer,dimension(:),allocatable :: recv_istatus,recv_rstatus

   ! Pointer structure to make Router access simpler
   type(Router), pointer :: SendRout, RecvRout

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   ! CHECK ARGUMENTS 

   ! Check the size of the Source AttrVect
   if(InRearranger%SendRouter%lAvsize /= AttrVect_lsize(SourceAV)) then
      call warn(myname_,"SourceAV size is not appropriate for this Rearranger")
      call die(myname_,"InRearranger%SendRouter%lAvsize",InRearranger%SendRouter%lAvsize, &
	        "AttrVect_lsize(SourceAV)", AttrVect_lsize(SourceAV))
   endif

   ! Check the size of the Target AttrVect
   if(InRearranger%RecvRouter%lAvsize /= AttrVect_lsize(TargetAV)) then
      call warn(myname_,"TargetAV size is not appropriate for this Rearranger")
      call die(myname_,"InRearranger%RecvRouter%lAvsize",InRearranger%RecvRouter%lAvsize, &
	        "AttrVect_lsize(TargetAV)", AttrVect_lsize(TargetAV))
   endif

   ! Check the number of integer attributes 
   if(nIAttr(SourceAV) /= nIAttr(TargetAV)) then
      call warn(myname_, &
                "Number of attributes in SourceAV and TargetAV do not match")
      call die(myname_,"nIAttr(SourceAV)", nIAttr(SourceAV), &
                        "nIAttr(TargetAV)", nIAttr(TargetAV))
   endif

   ! Check the number of real attributes
   if(nRAttr(SourceAV) /= nRAttr(TargetAV)) then
      call warn(myname_, &
      "Number of attributes in SourceAV and TargetAV do not match")
      call die(myname_,"nRAttr(SourceAV)", nRAttr(SourceAV), &
                        "nRAttr(TargetAV)", nRAttr(TargetAV))
   endif

   usevector=.false.
   if(present(Vector)) then
    if(Vector) usevector=.true.
   endif

   ! ASSIGN VARIABLES

   ! Get the number of integer and real attributes
   numi = nIAttr(SourceAV)
   numr = nRAttr(SourceAV)

   ! Assign the pointers
   nullify(SendRout,RecvRout)
   SendRout => InRearranger%SendRouter
   RecvRout => InRearranger%RecvRouter

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  ! ALLOCATE DATA STRUCTURES !

  ! IF SENDING DATA
  if(SendRout%nprocs > 0) then

     ! IF SENDING INTEGER DATA
     if(numi .ge. 1) then

	! allocate the number of pointers needed to send
	allocate(ip1(SendRout%nprocs),stat=ier)
	if(ier/=0) call die(myname_,'allocate(ip1)',ier)

	! allocate buffers to hold all outgoing data
	do proc=1,SendRout%nprocs
	   allocate(ip1(proc)%pi(SendRout%locsize(proc)*numi),stat=ier)
	   if(ier/=0) call die(myname_,'allocate(ip1%pi)',ier)
	enddo

	! allocate MPI send request array
	allocate(send_ireqs(SendRout%nprocs),stat=ier)
	if(ier/=0) call die(myname_,'allocate(send_ireqs)',ier)

	! allocatAe MPI status array
	allocate(send_istatus(MP_STATUS_SIZE,SendRout%nprocs),stat=ier)
	if(ier/=0) call die(myname_,'allocate(istatus)',ier)

     endif

     ! IF SENDING REAL DATA
     if(numr .ge. 1) then

	! allocate the number of pointers needed to send
	allocate(rp1(SendRout%nprocs),stat=ier)
	if(ier/=0) call die(myname_,'allocate(rp1)',ier)

	! allocate buffers to hold all outgoing data
	do proc=1,SendRout%nprocs
	   allocate(rp1(proc)%pr(SendRout%locsize(proc)*numr),stat=ier)
	   if(ier/=0) call die(myname_,'allocate(rp1%pr)',ier)
	enddo
	
	! allocate MPI send request array
	allocate(send_rreqs(SendRout%nprocs),stat=ier)
	if(ier/=0) call die(myname_,'allocate(send_rreqs)',ier)

	! allocate MPI status array
	allocate(send_rstatus(MP_STATUS_SIZE,SendRout%nprocs),stat=ier)
	if(ier/=0) call die(myname_,'allocate(rstatus)',ier)

	mp_Type_rp=MP_Type(rp1(1)%pr(1))

     endif

  endif

  ! IF RECEVING DATA
  if(RecvRout%nprocs > 0) then

     ! IF RECEIVING INTEGER DATA
     if(numi .ge. 1) then

	! allocate the number of pointers needed to receive
	allocate(ip2(RecvRout%nprocs),stat=ier)
	if(ier/=0) call die(myname_,'allocate(ip2)',ier)

	! allocate buffers to hold all incoming data
	do proc=1,RecvRout%nprocs
	   allocate(ip2(proc)%pi(RecvRout%locsize(proc)*numi),stat=ier)
	   if(ier/=0) call die(myname_,'allocate(ip2%pi)',ier)
	enddo
     
	! allocate MPI send request array
	allocate(recv_ireqs(RecvRout%nprocs),stat=ier)
	if(ier/=0) call die(myname_,'allocate(recv_ireqs)',ier)

	! allocate MPI integer status array
	allocate(recv_istatus(MP_STATUS_SIZE),stat=ier)
	if(ier/=0) call die(myname_,'allocate(istatus)',ier)

     endif

     ! IF RECEIVING REAL DATA
     if(numr .ge. 1) then

	! allocate the number of pointers needed to receive
	allocate(rp2(RecvRout%nprocs),stat=ier)
	if(ier/=0) call die(myname_,'allocate(rp2)',ier)

	! allocate buffers to hold all incoming data
	do proc=1,RecvRout%nprocs
	   allocate(rp2(proc)%pr(RecvRout%locsize(proc)*numr),stat=ier)
	   if(ier/=0) call die(myname_,'allocate(rp2%pr)',ier)
	enddo
     
	! allocate MPI receive request array
	allocate(recv_rreqs(RecvRout%nprocs),stat=ier)
	if(ier/=0) call die(myname_,'allocate(recv_rreqs)',ier)

	! allocate MPI real status array
	allocate(recv_rstatus(MP_STATUS_SIZE),stat=ier)
	if(ier/=0) call die(myname_,'allocate(rstatus)',ier)

	mp_Type_rp=MP_Type(rp2(2)%pr(1))

     endif

  endif

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  ! POST MPI_IRECV

  ! Load data coming from each processor
  do proc = 1,RecvRout%nprocs
    
     ! receive the integer data
     if(numi .ge. 1) then

        ! set tag
        mytag = DefaultTag
        if(present(Tag)) mytag=Tag

	if( (RecvRout%num_segs(proc) > 1) .or. present(Sum) ) then

	   call MPI_IRECV(ip2(proc)%pi(1),                        &
		          RecvRout%locsize(proc)*numi,MP_INTEGER, &
			  RecvRout%pe_list(proc),mytag,             &
			  ThisMCTWorld%MCT_comm,recv_ireqs(proc),ier)

	else

	   call MPI_IRECV(TargetAV%iAttr(1,RecvRout%seg_starts(proc,1)), &
		          RecvRout%locsize(proc)*numi,MP_INTEGER,        &
                          RecvRout%pe_list(proc),mytag,                    &
                          ThisMCTWorld%MCT_comm,recv_ireqs(proc),ier)

	endif

	if(ier /= 0) call MP_perr_die(myname_,'MPI_IRECV(ints)',ier)

     endif

     ! receive the real data
     if(numr .ge. 1) then

        ! set tag
        mytag = DefaultTag + 1
        if(present(Tag)) mytag=Tag +1

	if( (RecvRout%num_segs(proc) > 1) .or. present(Sum) ) then

	   call MPI_IRECV(rp2(proc)%pr(1),                        &
		          RecvRout%locsize(proc)*numr,mp_Type_rp, &
			  RecvRout%pe_list(proc),mytag,             &
			  ThisMCTWorld%MCT_comm,recv_rreqs(proc),ier)

	else

	   call MPI_IRECV(TargetAV%rAttr(1,RecvRout%seg_starts(proc,1)), &
		          RecvRout%locsize(proc)*numr,mp_Type_rp,        &
			  RecvRout%pe_list(proc),mytag,                    &
			  ThisMCTWorld%MCT_comm,recv_rreqs(proc),ier)

	endif

	if(ier /= 0) call MP_perr_die(myname_,'MPI_IRECV(reals)',ier)

     endif

  enddo

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  ! POST MPI_ISEND

  ! Load data going to each processor
  do proc = 1,SendRout%nprocs
    
     if( SendRout%num_segs(proc) > 1 ) then

	j=1
	k=1

	! load the correct pieces of the integer and real vectors
	do nseg = 1,SendRout%num_segs(proc)
	   seg_start = SendRout%seg_starts(proc,nseg)
	   seg_end = seg_start + SendRout%seg_lengths(proc,nseg)-1
	   do VectIndex = seg_start,seg_end
	      do AttrIndex = 1,numi
		 ip1(proc)%pi(j) = SourceAV%iAttr(AttrIndex,VectIndex)
		 j=j+1
	      enddo
	      do AttrIndex = 1,numr
		 rp1(proc)%pr(k) = SourceAV%rAttr(AttrIndex,VectIndex)
		 k=k+1
	      enddo
	   enddo
	enddo

     endif

     ! send the integer data
     if(numi .ge. 1) then

        ! set tag
        mytag = DefaultTag
        if(present(Tag)) mytag=Tag

	if( SendRout%num_segs(proc) > 1 ) then

	   call MPI_ISEND(ip1(proc)%pi(1),                        &
		          SendRout%locsize(proc)*numi,MP_INTEGER, &
			  SendRout%pe_list(proc),mytag,             &
			  ThisMCTWorld%MCT_comm,send_ireqs(proc),ier)

	else

	   call MPI_ISEND(SourceAV%iAttr(1,SendRout%seg_starts(proc,1)), &
		          SendRout%locsize(proc)*numi,MP_INTEGER,        &
                          SendRout%pe_list(proc),mytag,                  &
                          ThisMCTWorld%MCT_comm,send_ireqs(proc),ier)

	endif

	if(ier /= 0) call MP_perr_die(myname_,'MPI_ISEND(ints)',ier)

     endif

     ! send the real data
     if(numr .ge. 1) then

        ! set tag
        mytag = DefaultTag +1
        if(present(Tag)) mytag=Tag +1

	if( SendRout%num_segs(proc) > 1 ) then

	   call MPI_ISEND(rp1(proc)%pr(1),                        &
		          SendRout%locsize(proc)*numr,mp_Type_rp, &
			  SendRout%pe_list(proc),mytag,             &
			  ThisMCTWorld%MCT_comm,send_rreqs(proc),ier)

	else

	   call MPI_ISEND(SourceAV%rAttr(1,SendRout%seg_starts(proc,1)), &
		          SendRout%locsize(proc)*numr,mp_Type_rp,        &
			  SendRout%pe_list(proc),mytag,                    &
			  ThisMCTWorld%MCT_comm,send_rreqs(proc),ier)

	endif

	if(ier /= 0) call MP_perr_die(myname_,'MPI_ISEND(reals)',ier)

     endif

  enddo

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  ! ZERO TARGETAV WHILE WAITING FOR MESSAGES TO COMPLETE

  if(present(Sum)) call AttrVect_zero(TargetAV)

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  ! LOAD THE LOCAL PIECES OF THE INTEGER AND REAL VECTOR

  if(usevector) then
    do IAttrIndex=1,numi
!CDIR SELECT(VECTOR)
     do localindex=1,InRearranger%LocalSize
        TrgVectIndex = InRearranger%LocalPack(1,localindex)
        SrcVectIndex = InRearranger%LocalPack(2,localindex)
        TargetAV%iAttr(IAttrIndex,TrgVectIndex) = &
             SourceAV%iAttr(IAttrIndex,SrcVectIndex)
      enddo
    enddo
    do RAttrIndex=1,numr
!CDIR SELECT(VECTOR)
     do localindex=1,InRearranger%LocalSize
        TrgVectIndex = InRearranger%LocalPack(1,localindex)
        SrcVectIndex = InRearranger%LocalPack(2,localindex)
        TargetAV%rAttr(RAttrIndex,TrgVectIndex) = &
             SourceAV%rAttr(RAttrIndex,SrcVectIndex)
     enddo
    enddo

  else
    do localindex=1,InRearranger%LocalSize
     TrgVectIndex = InRearranger%LocalPack(1,localindex)
     SrcVectIndex = InRearranger%LocalPack(2,localindex)
     do IAttrIndex=1,numi
	TargetAV%iAttr(IAttrIndex,TrgVectIndex) = &
	     SourceAV%iAttr(IAttrIndex,SrcVectIndex)
     enddo
     do RAttrIndex=1,numr
	TargetAV%rAttr(RAttrIndex,TrgVectIndex) = &
	     SourceAV%rAttr(RAttrIndex,SrcVectIndex)
     enddo
    enddo
  endif

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  ! WAIT FOR THE NONBLOCKING SENDS TO COMPLETE

  if(SendRout%nprocs > 0) then

     if(numi .ge. 1) then

	call MPI_WAITALL(SendRout%nprocs,send_ireqs,send_istatus,ier)
	if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL(ints)',ier)

	! deallocatAe MPI status array
	deallocate(send_istatus,stat=ier)
	if(ier/=0) call die(myname_,'deallocate(istatus)',ier)

	! done waiting, free up ireqs
	deallocate(send_ireqs,stat=ier)
	if(ier /= 0) call die(myname_,'deallocate(Reqs%ireqs)',ier)

	! Deallocate the send buffers
	do proc=1,SendRout%nprocs
	   deallocate(ip1(proc)%pi,stat=ier)
	   if(ier/=0) call die(myname_,'deallocate(ip1%pr)',ier)
	enddo

	deallocate(ip1,stat=ier)
	if(ier/=0) call die(myname_,'deallocate(ip1)',ier)

     endif

     if(numr .ge. 1) then

	call MPI_WAITALL(SendRout%nprocs,send_rreqs,send_rstatus,ier)
	if(ier /= 0) call MP_perr_die(myname_,'MPI_WAITALL(reals)',ier)

	! deallocate MPI status array
	deallocate(send_rstatus,stat=ier)
	if(ier/=0) call die(myname_,'deallocate(rstatus)',ier)

	! done waiting, free up Reqs%rreqs
	deallocate(send_rreqs,stat=ier)
	if(ier /= 0) call die(myname_,'deallocate(Reqs%rreqs)',ier)

	! Deallocate the send buffers
	do proc=1,SendRout%nprocs
	   deallocate(rp1(proc)%pr,stat=ier)
	   if(ier/=0) call die(myname_,'deallocate(rp1%pr)',ier)
	enddo

	deallocate(rp1,stat=ier)
	if(ier/=0) call die(myname_,'deallocate(rp1)',ier)

     endif

  endif

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  ! WAIT FOR THE NONBLOCKING RECEIVES TO COMPLETE AND UNPACK BUFFER

  do numprocs = 1,RecvRout%nprocs

     if(numi .ge. 1) then

	call MPI_WAITANY(RecvRout%nprocs,recv_ireqs,proc,recv_istatus,ier)

	j=1

	if(present(Sum)) then

	   ! load the correct pieces of the integer vectors
	   do nseg = 1,RecvRout%num_segs(proc)
	      seg_start = RecvRout%seg_starts(proc,nseg)
	      seg_end = seg_start + RecvRout%seg_lengths(proc,nseg)-1
	      do VectIndex = seg_start,seg_end
		 do AttrIndex = 1,numi
		    TargetAV%iAttr(AttrIndex,VectIndex)= &
                    TargetAV%iAttr(AttrIndex,VectIndex) + ip2(proc)%pi(j)
		    j=j+1
		 enddo
	      enddo
	   enddo
	   
	else

	   if( RecvRout%num_segs(proc) > 1 ) then

	      ! load the correct pieces of the integer vectors
	      do nseg = 1,RecvRout%num_segs(proc)
		 seg_start = RecvRout%seg_starts(proc,nseg)
		 seg_end = seg_start + RecvRout%seg_lengths(proc,nseg)-1
		 do VectIndex = seg_start,seg_end
		    do AttrIndex = 1,numi
		       TargetAV%iAttr(AttrIndex,VectIndex)=ip2(proc)%pi(j)
		       j=j+1
		    enddo
		 enddo
	      enddo
	
	   endif

	endif

     endif

     if(numr .ge. 1) then

	call MPI_WAITANY(RecvRout%nprocs,recv_rreqs,proc,recv_rstatus,ier)

	k=1

	if(present(Sum)) then

	   ! load the correct pieces of the integer vectors
	   do nseg = 1,RecvRout%num_segs(proc)
	      seg_start = RecvRout%seg_starts(proc,nseg)
	      seg_end = seg_start + RecvRout%seg_lengths(proc,nseg)-1
	      do VectIndex = seg_start,seg_end
		 do AttrIndex = 1,numr
		    TargetAV%rAttr(AttrIndex,VectIndex) = &
                    TargetAV%rAttr(AttrIndex,VectIndex) + rp2(proc)%pr(k)
		    k=k+1
		 enddo
	      enddo
	   enddo
	   
	else

	   if( RecvRout%num_segs(proc) > 1 ) then

	      ! load the correct pieces of the integer vectors
	      do nseg = 1,RecvRout%num_segs(proc)
		 seg_start = RecvRout%seg_starts(proc,nseg)
		 seg_end = seg_start + RecvRout%seg_lengths(proc,nseg)-1
		 do VectIndex = seg_start,seg_end
		    do AttrIndex = 1,numr
		       TargetAV%rAttr(AttrIndex,VectIndex)=rp2(proc)%pr(k)
		       k=k+1
		    enddo
		 enddo
	      enddo
	   
	   endif

	endif

     endif

  enddo

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  ! DEALLOCATE ALL STRUCTURES

  if(RecvRout%nprocs > 0) then

     if(numi .ge. 1) then

	! deallocatAe MPI status array
	deallocate(recv_istatus,stat=ier)
	if(ier/=0) call die(myname_,'deallocate(istatus)',ier)

	! done waiting, free up ireqs
	deallocate(recv_ireqs,stat=ier)
	if(ier /= 0) call die(myname_,'deallocate(Reqs%ireqs)',ier)

	! Deallocate the send buffers
	do proc=1,RecvRout%nprocs
	   deallocate(ip2(proc)%pi,stat=ier)
	   if(ier/=0) call die(myname_,'deallocate(ip1%pr)',ier)
	enddo

	deallocate(ip2,stat=ier)
	if(ier/=0) call die(myname_,'deallocate(ip1)',ier)

     endif

     if(numr .ge. 1) then

	! deallocate MPI status array
	deallocate(recv_rstatus,stat=ier)
	if(ier/=0) call die(myname_,'deallocate(rstatus)',ier)

	! done waiting, free up Reqs%rreqs
	deallocate(recv_rreqs,stat=ier)
	if(ier /= 0) call die(myname_,'deallocate(Reqs%rreqs)',ier)

	! Deallocate the send buffers
	do proc=1,RecvRout%nprocs
	   deallocate(rp2(proc)%pr,stat=ier)
	   if(ier/=0) call die(myname_,'deallocate(rp1%pr)',ier)
	enddo

	deallocate(rp2,stat=ier)
	if(ier/=0) call die(myname_,'deallocate(rp1)',ier)

     endif

  endif

  nullify(SendRout,RecvRout)


 end subroutine rearrange_




end module m_Rearranger





