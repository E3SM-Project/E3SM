!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_AttrVectComms - Communications methods for the AttrVect
!
! !DESCRIPTION:
!
! In this module, we define communications methods specific to the 
! {\tt AttrVect} class (see the module {\tt m\_AttrVect} for more 
! information about this class and its methods).
!
! !INTERFACE:
 module m_AttrVectComms
!
! !USES:
!
      use m_AttrVect ! AttrVect class and its methods


      implicit none

      private	! except

      public :: gather		! gather all local vectors to the root
      public :: scatter		! scatter from the root to all PEs
      public :: bcast		! bcast from root to all PEs

    interface gather ; module procedure &
	      GM_gather_, &
	      GSM_gather_ 
    end interface
    interface scatter ; module procedure &
	      GM_scatter_, &
	      GSM_scatter_ 
    end interface
    interface bcast  ; module procedure bcast_  ; end interface

! !REVISION HISTORY:
! 	27Oct00 - J.W. Larson <larson@mcs.anl.gov> - relocated routines
!                 from m_AttrVect to create this module.
!       15Jan01 - J.W. Larson <larson@mcs.anl.gov> - Added APIs for 
!                 GSM_gather_() and GSM_scatter_().
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_AttrVectComms'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GM_gather_ - gather a vector using input GlobalMap.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine GM_gather_(iV,oV,GMap,root,comm,stat)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90
      use m_GlobalMap, only : GlobalMap
      use m_GlobalMap, only : GlobalMap_lsize => lsize
      use m_GlobalMap, only : GlobalMap_gsize => gsize
      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_nIAttr => nIAttr
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr

      implicit none

      type(AttrVect),intent(in)  :: iV
      type(AttrVect),intent(out) :: oV
      type(GlobalMap) ,intent(in)  :: GMap
      integer, intent(in) :: root
      integer, intent(in) :: comm
      integer, optional,intent(out) :: stat

! !REVISION HISTORY:
! 	15Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	27Oct00 - J.W. Larson <larson@mcs.anl.gov> - relocated from
!                 m_AttrVect
!       15Jan01 - J.W. Larson <larson@mcs.anl.gov> - renamed GM_gather_
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GM_gather_'
  integer :: nIA,nRA,niV,noV,ier
  integer :: myID

  if(present(stat)) stat=0

  call MP_comm_rank(comm,myID,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MP_comm_rank()',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

	! Verify the input: a _scatterd_ vector

  niV=GlobalMap_lsize(GMap)
  noV=AttrVect_lsize(iV)

  if(niV /= noV) then
    write(stderr,'(2a,i4,a,i4)') myname_,	&
	': invalid input, lsize(GMap) =',niV,	&
	', lsize(iV) =',noV
    if(.not.present(stat)) call die(myname_)
    stat=-1
    return
  endif

  noV=GlobalMap_gsize(GMap) ! the gathered local size, as for the output
  if(myID /= root) nov=0
  call AttrVect_init(oV,iV,noV)

  niV=GlobalMap_lsize(GMap) ! the scattered local size, as for the input

  nIA=AttrVect_nIAttr(oV)	! number of INTEGER attributes
  nRA=AttrVect_nRAttr(oV)	! number of REAL attributes

  call MPI_gatherv(iV%iAttr,niV*nIA,MP_INTEGER,			&
	oV%iAttr,GMap%counts*nIA,GMap%displs*nIA,MP_INTEGER,	&
	root,comm,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MPI_gatherv(iAttr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

  call MPI_gatherv(iV%rAttr,niV*nRA,MP_REAL,		&
	oV%rAttr,GMap%counts*nRA,GMap%displs*nRA,MP_REAL,	&
	root,comm,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MPI_gatherv(rAttr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

 end subroutine GM_gather_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GSM_gather_ - gather a vector using input GlobalSegMap.
!
! !DESCRIPTION:
! The routine {\tt GSM\_gather\_()} takes a distributed input 
! {\tt AttrVect} argument {\tt iV}, whose decomposition is described 
! by the input {\tt GlobalSegMap} argument {\tt GSMap}, and gathers 
! it to the output {\tt AttrVect} argument {\tt oV}.  The gathered 
! {\tt AttrVect} {\tt oV} is valid only on the root process specified 
! by the input argument {\tt root}.  The communicator used to gather
! the data is specified by the argument {\tt comm}.  The success (failure)
! is reported in the zero (non-zero) value of the output argument 
! {\tt stat}.
!
! {\tt GSM\_gather\_()} converts the problem of gathering data 
! according to a {\tt GlobalSegMap} into the simpler problem of 
! gathering data as specified by a {\tt GlobalMap}.  The {\tt GlobalMap}
! variable {\tt GMap} is created based on the local storage requirements 
! for each distributed piece of {\tt iV}.  On the root, a complete 
! (including halo points) gathered copy of {\tt iV} is collected into 
! the temporary {\tt AttrVect} variable {\tt workV} (the length of
! {\tt workV} is {\tt GlobalSegMap\_GlobalStorage(GSMap)}).  The 
! variable {\tt workV} is segmented by process, and segments are 
! copied into it by process, but ordered in the same order the segments
! appear in {\tt GSMap}.  Once {\tt workV} is loaded, the data are 
! copied segment-by-segment to their appropriate locations in the output 
! {\tt AttrVect} {\tt oV}.
!
! !INTERFACE:

 subroutine GSM_gather_(iV, oV, GSMap, root, comm, stat)
!
! !USES:
!
! Message-passing environment utilities (mpeu) modules:
      use m_stdio
      use m_die,          only : MP_perr_die
      use m_mpif90
! GlobalSegMap and associated services:
      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_comp_id => comp_id
      use m_GlobalSegMap, only : GlobalSegMap_ngseg => ngseg
      use m_GlobalSegMap, only : GlobalSegMap_lsize => lsize
      use m_GlobalSegMap, only : GlobalSegMap_gsize => gsize
      use m_GlobalSegMap, only : GlobalSegMap_haloed => haloed
      use m_GlobalSegMap, only : GlobalSegMap_ProcessStorage => ProcessStorage
      use m_GlobalSegMap, only : GlobalSegMap_GlobalStorage => GlobalStorage
! AttrVect and associated services:
      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_nIAttr => nIAttr
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr
! GlobalMap and associated services:
      use m_GlobalMap, only : GlobalMap
      use m_GlobalMap, only : GlobalMap_init => init
      use m_GlobalMap, only : GlobalMap_clean => clean

      implicit none

      type(AttrVect),     intent(in)  :: iV
      type(AttrVect),     intent(out) :: oV
      type(GlobalSegMap) ,intent(in)  :: GSMap
      integer,            intent(in)  :: root
      integer,            intent(in)  :: comm
      integer, optional,  intent(out) :: stat

! !REVISION HISTORY:
! 	15Jan01 - J.W. Larson <larson@mcs.anl.gov> - API specification.
! 	25Feb01 - J.W. Larson <larson@mcs.anl.gov> - Prototype code.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GSM_gather_'

! Temporary workspace AttrVect:
  type(AttrVect) :: workV
! Component ID and number of segments for GSMap:
  integer :: comp_id, ngseg
! Total length of GSMap segments laid end-to-end:
  integer :: global_storage
! Error Flag
  integer :: ierr
! Number of processes on communicator, and local rank:
  integer :: NumProcs, myID
! Total local storage on each pe according to GSMap:
  integer, dimension(:), allocatable :: lns
! Temporary GlobalMap used to scatter the segmented (by pe) data
  type(GlobalMap) :: workGMap
! Loop counters and temporary indices:
  integer :: m, n, ilb, iub, olb, oub, pe
! workV segment tracking index array:
  integer, dimension(:), allocatable :: current_pos

       ! Initial Check:  If GSMap contains halo points, die
       
  if(GlobalSegMap_haloed(GSMap)) then
     ierr = 1
     call MP_perr_die(myname_,"Input GlobalSegMap haloed--not allowed", &
	              ierr)
  endif

       ! Which process am I?

  call MPI_COMM_RANK(comm, myID, ierr)

  if(myID == root) then

       ! How many processes are there on this communicator?

     call MPI_COMM_SIZE(comm, NumProcs, ierr)

       ! Allocate a precursor to a GlobalMap accordingly...

     allocate(lns(0:NumProcs-1), stat=ierr)

       ! And Load it...

     do n=0,NumProcs-1
        lns(n) = GlobalSegMap_ProcessStorage(GSMap, comm)
     end do

  endif ! if(myID == root)

       ! Determine the component id of GSMap:

  comp_id = GlobalSegMap_comp_id(GSMap)

       ! Create working GlobalMap workGMap (used for the gather):

  call GlobalMap_init(workGMap, comp_id, lns, root, comm)

       ! Gather the Data process-by-process to workV...

  call GM_gather_(iV, workV, workGMap, root, comm, stat)

       ! On the root, initialize oV, and load the contents of
       !workV into it...

  if(myID == root) then

       ! On the root, allocate current position index for
       ! each process chunk:

     allocate(current_pos(0:NumProcs-1), stat=ierr)

       ! Initialize current_pos(:) using GMap%displs(:)

     do n=0,NumProcs-1
	current_pos(n) = workGMap%displs(n)
     end do

       ! Load each segment of iV into its appropriate segment
       ! of workV:

     ngseg = GlobalSegMap_ngseg(GSMap)

     do n=1,ngseg

       ! Determine which process owns segment n:

	pe = GSMap%pe_loc(n)

       ! Input map (lower/upper indicess) of segment of iV:

	ilb = current_pos(pe)
	iub = current_pos(pe) + GSMap%length(n) - 1

       ! Output map of (lower/upper indicess) segment of workV:

	olb = GSMap%start(n)
	oub = GSMap%start(n) + GSMap%length(n) - 1

       ! Increment current_pos(n) for next time:

	current_pos(pe) = current_pos(pe) + GSMap%length(n)

       ! Now we are equipped to do the copy:

	do m=1,AttrVect_nIAttr(oV)
   	   oV%iAttr(m,olb:oub) = workV%iAttr(m,ilb:iub)
	end do

	do m=1,AttrVect_nRAttr(iV)
   	   oV%rAttr(m,olb:oub) = workV%rAttr(m,ilb:iub)
	end do

     end do ! do n=1,ngseg

       ! Clean up current_pos, which was only allocated on the root

     deallocate(current_pos, stat=ierr)

  endif ! if(myID == root)

       ! At this point, we are finished.  The data have been gathered
       ! to oV

       ! Finally, clean up allocated structures:

  call AttrVect_clean(workV)
  call GlobalMap_clean(workGMap)
  deallocate(lns, stat=ierr)

 end subroutine GSM_gather_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GM_scatter_ - scatter a vecter using input GlobalMap.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine GM_scatter_(iV, oV, GMap, root, comm, stat)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90
      use m_GlobalMap, only : GlobalMap
      use m_GlobalMap, only : GlobalMap_lsize => lsize
      use m_GlobalMap, only : GlobalMap_gsize => gsize
      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect,  only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_nIAttr => nIAttr
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr

      implicit none

      type(AttrVect),intent(in)  :: iV
      type(AttrVect),intent(out) :: oV
      type(GlobalMap) ,intent(in)  :: GMap
      integer, intent(in) :: root
      integer, intent(in) :: comm
      integer, optional,intent(out) :: stat

! !REVISION HISTORY:
! 	21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	27Oct00 - J.W. Larson <larson@mcs.anl.gov> - relocated from
!                 m_AttrVect
! 	15Jan01 - J.W. Larson <larson@mcs.anl.gov> - renamed  GM_scatter_
! 	08Feb01 - J.W. Larson <larson@mcs.anl.gov> - add logic to prevent
!                 empty calls (i.e. no data in buffer) to MPI_SCATTERV()
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GM_scatter_'
  integer :: nIA,nRA,niV,noV,ier
  integer :: myID

  if(present(stat)) stat=0

  call MP_comm_rank(comm,myID,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MP_comm_rank()',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

	! Verify the input: a _gathered_ vector

  niV = GlobalMap_gsize(GMap)  ! the _gathered_ local size
  if(myID /= root) niV=0

  noV = AttrVect_lsize(iV)
  if(niV /= noV) then
    write(stderr,'(2a,i4,a,i4)') myname_,	&
	': invalid input, rsize(GMap) =',niV,	&
	', lsize(iV) =',noV
    if(.not.present(stat)) call die(myname_)
    stat=-1
    return
  endif

  noV = GlobalMap_lsize(GMap) ! the _scatterd_ local size
  call AttrVect_initv(oV,iV,noV)

  nIA = AttrVect_nIAttr(iV)	! number of INTEGER attributes
  nRA = AttrVect_nRAttr(iV)	! number of REAL attributes

  if(nIA > 0) then
     call MPI_scatterv(iV%iAttr(1,1),GMap%counts*nIA,	&
	  GMap%displs*nIA,MP_INTEGER,			&
	  oV%iAttr(1,1),noV*nRA,MP_INTEGER,root,comm,ier )
     if(ier /= 0) then
	call MP_perr(myname_,'MPI_scatterv(iAttr)',ier)
	if(.not.present(stat)) call die(myname_)
	stat=ier
	return
     endif
  endif

  if(nRA > 0) then
     call MPI_scatterv(iV%rAttr(1,1),GMap%counts*nRA,	&
	  GMap%displs*nRA,MP_REAL,				&
	  oV%rAttr(1,1),noV*nRA,MP_REAL,root,comm,ier )
     if(ier /= 0) then
	call MP_perr(myname_,'MPI_scatterv(rAttr)',ier)
	if(.not.present(stat)) call die(myname_)
	stat=ier
	return
     endif
  endif

 end subroutine GM_scatter_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GSM_scatter_ - scatter a vecter using input GlobalSegMap.
!
! !DESCRIPTION:
! The routine {\tt GSM\_scatter\_} takes an input {\tt AttrVect} type
! {\tt iV} located on the root, and scatters it to a distributed 
! {\tt AttrVect} {\tt oV}.  The input {\tt GlobalSegMap} argument 
! {\tt GSMap} dictates how {\tt iV} is scattered to {\tt oV}.  The 
! success (failure) of this routine is reported in the zero (non-zero)
! value of the output argument {\tt stat}.
!
! {\tt GSM\_scatter\_()} converts the problem of scattering data 
! according to a {\tt GlobalSegMap} into the simpler problem of 
! scattering data as specified by a {\tt GlobalMap}.  The {\tt GlobalMap}
! variable {\tt GMap} is created based on the local storage requirements 
! for each distributed piece of {\tt iV}.  On the root, a complete 
! (including halo points) copy of {\tt iV} is stored in 
! the temporary {\tt AttrVect} variable {\tt workV} (the length of
! {\tt workV} is {\tt GlobalSegMap\_GlobalStorage(GSMap)}).  The 
! variable {\tt workV} is segmented by process, and segments are 
! copied into it by process, but ordered in the same order the segments
! appear in {\tt GSMap}.  Once {\tt workV} is loaded, the data are 
! scattered to the output {\tt AttrVect} {\tt oV} by a call to the
! routine {\tt GM\_scatter\_()} defined in this module, with {\tt workV}
! and {\tt GMap} as the input arguments.
!
! {\bf N.B.:}  This algorithm  assumes that memory access times are much 
! shorter than message-passing transmission times.
!
! !INTERFACE:

 subroutine GSM_scatter_(iV, oV, GSMap, root, comm, stat)
!
! !USES:
!
! Environment utilities from mpeu:
      use m_stdio
      use m_die
      use m_mpif90
! GlobalSegMap and associated services:
      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_comp_id => comp_id
      use m_GlobalSegMap, only : GlobalSegMap_ngseg => ngseg
      use m_GlobalSegMap, only : GlobalSegMap_lsize => lsize
      use m_GlobalSegMap, only : GlobalSegMap_gsize => gsize
      use m_GlobalSegMap, only : GlobalSegMap_GlobalStorage => GlobalStorage
      use m_GlobalSegMap, only : GlobalSegMap_ProcessStorage => ProcessStorage
! AttrVect and associated services:
      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect,  only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_nIAttr => nIAttr
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr
! GlobalMap and associated services:
      use m_GlobalMap, only : GlobalMap
      use m_GlobalMap, only : GlobalMap_init => init
      use m_GlobalMap, only : GlobalMap_clean => clean

      implicit none

      type(AttrVect),intent(in)  :: iV
      type(AttrVect),intent(out) :: oV
      type(GlobalSegMap) ,intent(in)  :: GSMap
      integer, intent(in) :: root
      integer, intent(in) :: comm
      integer, optional,intent(out) :: stat

! !REVISION HISTORY:
! 	15Jan01 - J.W. Larson <larson@mcs.anl.gov> - API specification.
! 	08Feb01 - J.W. Larson <larson@mcs.anl.gov> - Initial code.
! 	25Feb01 - J.W. Larson <larson@mcs.anl.gov> - Bug fix--replaced
!                 call to GlobalSegMap_lsize with call to the new fcn.
!                 GlobalSegMap_ProcessStorage().
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GSM_scatter_'

! Temporary workspace AttrVect:
  type(AttrVect) :: workV
! Component ID and number of segments for GSMap:
  integer :: comp_id, ngseg
! Total length of GSMap segments laid end-to-end:
  integer :: global_storage
! Error Flag
  integer :: ierr
! Number of processes on communicator, and local rank:
  integer :: NumProcs, myID
! Total local storage on each pe according to GSMap:
  integer, dimension(:), allocatable :: lns
! Temporary GlobalMap used to scatter the segmented (by pe) data
  type(GlobalMap) :: GMap
! Loop counters and temporary indices:
  integer :: m, n, ilb, iub, olb, oub, pe
! workV segment tracking index array:
  integer, dimension(:), allocatable :: current_pos


       ! Which process am I?

  call MPI_COMM_RANK(comm, myID, ierr)

  if(myID == root) then

       ! Determine the _total_ storage space for all the 
       ! segments of GSMap laid end-to-end:

     global_storage = GlobalSegMap_GlobalStorage(GSMap)

       ! Initialize a work AttrVect type of the above 
       ! length, and with the same attribute lists as iV:

     call AttrVect_init(workV, iV, global_storage)

       ! How many processes are there on this communicator?

     call MPI_COMM_SIZE(comm, NumProcs, ierr)

       ! Allocate a precursor to a GlobalMap accordingly...

     allocate(lns(0:NumProcs-1), stat=ierr)

       ! And Load it...

     do n=0,NumProcs-1
	lns(n) = GlobalSegMap_ProcessStorage(GSMap, n)
     end do

  endif ! if(myID == root)

       ! Create a GlobalMap describing the 1-D decomposition 
       ! of workV:

  comp_id = GlobalSegMap_comp_id(GSMap)
  call GlobalMap_init(GMap, comp_id, lns, root, comm)

       ! On the root, load workV:

  if(myID == root) then

       ! On the root, allocate current position index for
       ! each process chunk:

     allocate(current_pos(0:NumProcs-1), stat=ierr)

       ! Initialize current_pos(:) using GMap%displs(:)

     do n=0,NumProcs-1
	current_pos(n) = GMap%displs(n)
     end do

       ! Load each segment of iV into its appropriate segment
       ! of workV:

     ngseg = GlobalSegMap_ngseg(GSMap)

     do n=1,ngseg

       ! Determine which process owns segment n:

	pe = GSMap%pe_loc(n)

       ! Input map (lower/upper indicess) of segment of iV:

	ilb = GSMap%start(n)
	iub = GSMap%start(n) + GSMap%length(n) - 1

       ! Output map of (lower/upper indicess) segment of workV:

	olb = current_pos(pe)
	oub = current_pos(pe) + GSMap%length(n) - 1

       ! Increment current_pos(n) for next time:

	current_pos(pe) = current_pos(pe) + GSMap%length(n)

       ! Now we are equipped to do the copy:

	do m=1,AttrVect_nIAttr(iV)
   	   workV%iAttr(m,olb:oub) = iV%iAttr(m,ilb:iub)
	end do

	do m=1,AttrVect_nRAttr(iV)
   	   workV%rAttr(m,olb:oub) = iV%rAttr(m,ilb:iub)
	end do

     end do ! do n=1,ngseg

       ! Clean up current_pos, which was only allocated on the root

     deallocate(current_pos, stat=ierr)


  endif ! if(myID == root)

       ! Now we are in business...we have:  1) an AttrVect laid out
       ! in contiguous segments, each segment corresponding to a 
       ! process, and in the same order dictated by GSMap; 
       ! 2) a GlobalMap telling us which segment of workV goes to 
       ! which process.  Thus, we can us GM_scatter_() to achieve 
       ! our goal.

  call GM_scatter_(workV, oV, GMap, root, comm, stat)

       ! Finally, clean up allocated structures:

  call AttrVect_clean(workV)
  call GlobalMap_clean(GMap)
  deallocate(lns, stat=ierr)

 end subroutine GSM_scatter_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: bcast_ - broadcast from the root to all PEs
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine bcast_(aV,root,comm,stat)
!
! !USES:
!
      use m_die, only : die, perr
      use m_mpif90
      use m_String, only : String,bcast,char
      use m_String, only : String_bcast => bcast
      use m_List, only : List_get => get
      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_nIAttr => nIAttr
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr

      implicit none

      type(AttrVect) :: aV	! (IN) on the root, (OUT) elsewhere
      integer,intent(in) :: root
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	27Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	27Oct00 - J.W. Larson <larson@mcs.anl.gov> - relocated from
!                 m_AttrVect
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::bcast_'
  type(String) :: iLStr,rLStr
  integer :: nIA, nRA, lsize
  integer :: myID
  integer :: ier

  if(present(stat)) stat=0

  call MP_comm_rank(comm,myID,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MP_comm_rank()',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

	! Convert the two Lists to two Strings

  if(myID == root)	&
    call List_get(iLStr,aV%iList)

  call String_bcast(iLStr,root,comm,stat=ier)	! bcast.String()
  if(ier /= 0) then
    call perr(myname_,'bcast.String(iLstr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

  if(myID == root)	&
    call List_get(rLStr,aV%rList)

  call String_bcast(rLStr,root,comm,stat=ier)	! bcast.String()
  if(ier /= 0) then
    call perr(myname_,'bcast.String(rLstr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

  if(myID == root) lsize = AttrVect_lsize(aV)

	! Set lsize for all PEs

  call MPI_bcast(lsize,1,MP_INTEGER,root,comm,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MPI_bcast(lsize)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

  if(myID /= root) then
     call AttrVect_init(aV,iList=char(iLStr),rList=char(rLStr), &
	                lsize=lsize)
  endif

  nIA = AttrVect_nIAttr(aV)
  nRA = AttrVect_nRAttr(aV)

  call MPI_bcast(aV%iAttr,nIA*lsize,MP_INTEGER,root,comm,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MPI_bcast(iAttr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

  call MPI_bcast(aV%rAttr,nRA*lsize,MP_REAL,   root,comm,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MPI_bcast(rAttr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

 end subroutine bcast_

 end module m_AttrVectComms
