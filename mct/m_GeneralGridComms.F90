!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_GeneralGridComms - Communications methods for the GeneralGrid
!
! !DESCRIPTION:
!
! In this module, we define communications methods specific to the 
! {\tt GeneralGrid} class (see the module {\tt m\_GeneralGrid} for more 
! information about this class and its methods).
!
! !INTERFACE:
 module m_GeneralGridComms
!
! !USES:
!
      use m_GeneralGrid ! GeneralGrid class and its methods


      implicit none

      private   ! except

      public :: gather          ! gather all local vectors to the root
      public :: scatter         ! scatter from the root to all PEs
      public :: bcast           ! bcast from root to all PEs

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
!       27Apr01 - J.W. Larson <larson@mcs.anl.gov> - Initial module/APIs
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_GeneralGridComms'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: send_ - Point-to-point send for the GeneralGrid.
!
! !DESCRIPTION:  The point-to-point send routine {\tt send\_()} sends 
! the input {\tt GeneralGrid} argument {\tt iGGrid} to process {\tt dest} 
! on the communicator associated with the F90 integer handle {\tt comm}.
! The message is identified by the tag defined by the {\tt INTEGER} 
! argument {\tt TagBase}.  The success (failure) of this operation 
! corresponds to a zero (nonzero) value for the output {\tt INTEGER} 
! flag {\tt status}. 
!
! !INTERFACE:

 subroutine send_(iGGrid, dest, TagBase, comm, status)

!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90

      use m_GeneralGrid, only : GeneralGrid
      use m_GeneralGrid, only : GeneralGrid_init => init
      use m_GeneralGrid, only : GeneralGrid_lsize => lsize

      implicit none

! !INPUT PARAMETERS: 
!
      type(GeneralGrid), intent(in) :: iGGrid
      integer,           intent(in) :: dest
      integer,           intent(in) :: TagBase
      integer,           intent(in) :: comm

! !OUTPUT PARAMETERS: 
!
      integer, optional, intent(out) :: status

! !REVISION HISTORY:
!       04Jun01 - J.W. Larson <larson@mcs.anl.gov> - API Specification.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::send_'

  integer :: ierr
  logical :: HeaderAssoc(6)

      ! Step 1. Check elements of the GeneralGrid header to see 
      ! which components of it are allocated.  Load the results
      ! into HeaderAssoc(:).

  HeaderAssoc(1) = associated(iGGrid%coordinate_list%bf)
  HeaderAssoc(2) = associated(iGGrid%coordinate_sort_order%bf)
  HeaderAssoc(3) = associated(iGGrid%descend)
  HeaderAssoc(4) = associated(iGGrid%weight_list%bf)
  HeaderAssoc(5) = associated(iGGrid%other_list%bf)
  HeaderAssoc(6) = associated(iGGrid%index_list%bf)

  call MPI_SEND(HeaderAssoc, 6, MP_LOGICAL, dest, TagBase, comm, ierr)
  if(ierr /= 0) then
     if(present(status)) then
        write(stderr,*) myname_,':: MPI_SEND(HeaderAssoc...'
        status = ierr
        return
     else
        call MP_perr_die(myname_,':: MPI_SEND(HeaderAssoc...',ierr)
     endif
  endif

  if(HeaderAssoc(1)) then
    call List_send(iGGrid%coordinate_list, dest, TagBase+1, comm, ierr)
    if(ierr /= 0) then
       if(present(status)) then
          write(stderr,*) myname_,':: call List_send(iGGrid%coordinate_list...'
          status = ierr
          return
       else
          call MP_perr_die(myname_,':: call List_send(iGGrid%coordinate_list...',ierr)
       endif
    endif
  endif

  if(HeaderAssoc(2)) then
    call List_send(iGGrid%coordinate_sort_order, dest, TagBase+3, comm, ierr)
    if(ierr /= 0) then
       if(present(status)) then
          write(stderr,*) myname_,':: call List_send(iGGrid%coordinate_sort_order...'
          status = ierr
          return
       else
          call MP_perr_die(myname_,':: call List_send(iGGrid%coordinate_sort_order...',ierr)
       endif
    endif
  endif
     
  if(HeaderAssoc(3)) then

       ! 

    call List_send(iGGrid%coordinate_sort_order, dest, TagBase+3, comm, ierr)
    if(ierr /= 0) then
       if(present(status)) then
          write(stderr,*) myname_,':: call List_send(iGGrid%coordinate_sort_order...'
          status = ierr
          return
       else
          call MP_perr_die(myname_,':: call List_send(iGGrid%coordinate_sort_order...',ierr)
       endif
    endif
  endif

 end subroutine send_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: recv_ - Point-to-point recv for the GeneralGrid.
!
! !DESCRIPTION:  The point-to-point receive routine {\tt recv\_()} 
! receives the output {\tt GeneralGrid} argument {\tt oGGrid} from process 
! {\tt source} on the communicator associated with the F90 integer handle 
! {\tt comm}.  The message is identified by the tag defined by the 
! {\tt INTEGER} argument {\tt tag}.  The success (failure) of this operation 
! corresponds to a zero (nonzero) value for the output {\tt INTEGER} flag 
! {\tt status}. 
!
! {\bf N.B.}:  This routine assumes that the {\tt GeneralGrid} argument
! {\tt oGGrid} is uninitialized on input; that is, all the {\tt List} 
! components are blank, the {\tt LOGICAL} array {\tt oGGrid\%descend} is
! unallocated, and the {\tt AttrVect} component {\tt oGGrid\%data} is
! uninitialized.
!
! !INTERFACE:

 subroutine recv_(oGGrid, source, tag_base, comm, status)

!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90

      use m_GeneralGrid, only : GeneralGrid
      use m_GeneralGrid, only : GeneralGrid_init => init
      use m_GeneralGrid, only : GeneralGrid_lsize => lsize

      implicit none

! !INPUT PARAMETERS: 
!
      integer,           intent(in) :: source
      integer,           intent(in) :: tag_base
      integer,           intent(in) :: comm

! !OUTPUT PARAMETERS: 
!
      type(GeneralGrid), intent(out) :: oGGrid
      integer, optional, intent(out) :: status

! !REVISION HISTORY:
!       04Jun01 - J.W. Larson <larson@mcs.anl.gov> - API Specification.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::recv_'

 end subroutine recv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GM_gather_ - gather a GeneralGrid using input GlobalMap.
!
! !DESCRIPTION:  {\tt GM\_gather\_()} takes an input {\tt GeneralGrid} 
! argument {\tt iG} whose decomposition on the communicator associated 
! with the F90 handle {\tt comm} is described by the {\tt GlobalMap} 
! argument {\tt GMap}, and gathers it to the {\tt GeneralGrid} output
! argument {\tt oG} on the {\tt root}.  The success (failure) of this 
! operation is reported as a zero (nonzero) value in the optional 
! {\tt INTEGER} output argument {\tt stat}.

! {\bf N.B.}:  An important assumption made here is that the distributed
! {\tt GeneralGrid} {\tt iG} has been initialized with the same 
! coordinate system, sort order, other real attributes, and the same 
! indexing attributes for all processes on {\tt comm}.
!
! {\bf N.B.}:  Once the gridpoint data of the {\tt GeneralGrid} are assembled 
! on the {\tt root}, they are stored in the order determined by the input 
! {\tt GlobalMap} {\tt GMap}.  The user may need to sorted these gathered
! data to order them in  accordance with the {\tt coordinate\_sort\_order} 
! attribute of {\tt iG}.
!
! {\bf N.B.}:  The output {\tt GeneralGrid} {\tt oG} represents allocated
! memory on the {\tt root}.  When the user no longer needs {\tt oG} it
! should be deallocated using {\tt GeneralGrid\_clean()} to avoid a memory
! leak
!
! !INTERFACE:
!
 subroutine GM_gather_(iG, oG, GMap, root, comm, stat)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90

      use m_GlobalMap, only : GlobalMap
      use m_GlobalMap, only : GlobalMap_gsize => gsize

      use m_GeneralGrid, only : GeneralGrid
      use m_GeneralGrid, only : GeneralGrid_init => init

      implicit none

! !INPUT PARAMETERS: 
!
      type(GeneralGrid), intent(in)  :: iG
      type(GlobalMap),   intent(in)  :: GMap
      integer,           intent(in)  :: root
      integer,           intent(in)  :: comm

! !OUTPUT PARAMETERS: 
!
      type(GeneralGrid), intent(out) :: oG
      integer, optional, intent(out) :: stat

! !REVISION HISTORY:
!       27Apr01 - J.W. Larson <larson@mcs.anl.gov> - API Specification.
!       02May01 - J.W. Larson <larson@mcs.anl.gov> - Initial code.
!EOP ___________________________________________________________________

 character(len=*),parameter :: myname_=myname//'::GM_gather_'
!Process ID
 integer :: myID
!Error flag
 integer :: ierr
!Number of points on the _Gathered_ grid:
 integer :: length

       ! Which process am I?

  call MPI_COMM_RANK(comm, myID, ierr)

  if(ierr /= 0) then
    call MP_perr(myname_,'MPI_COMM_RANK()',ierr)
    if(.not.present(stat)) call die(myname_)
    stat=ierr
    return
  endif

  if(myID == root) then ! prepare oG:

       ! The length of the _gathered_ GeneralGrid oG is determined by 
       ! the GlobalMap function GlobalMap_gsize()

     length = GlobalMap_gsize(GMap)

       ! Initialize attributes of oG from iG, and length

     call GeneralGrid_init(oG, iG, length)

  endif

       ! Gather gridpoint data in iG%data to oG%data

  call AttrVect_Gather(oG%data, iG%data, GMap, root, comm, ierr)

  if(ierr /= 0) then
    call MP_perr(myname_,'MPI_COMM_RANK()',ierr)
    if(.not.present(stat)) call die(myname_)
    stat=ierr
    return
  else
    if(present(stat)) stat=ierr
  endif

 end subroutine GM_gather_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GSM_gather_ - gather a GeneralGrid using input GlobalSegMap.
!
! !DESCRIPTION:  {\tt GMS\_gather\_()} takes an input {\tt GeneralGrid} 
! argument {\tt iG} whose decomposition on the communicator associated 
! with the F90 handle {\tt comm} is described by the {\tt GlobalSegMap} 
! argument {\tt GSMap}, and gathers it to the {\tt GeneralGrid} output
! argument {\tt oG} on the {\tt root}.  The success (failure) of this 
! operation is reported as a zero (nonzero) value in the optional 
! {\tt INTEGER} output argument {\tt stat}.
!
! {\bf N.B.}:  An important assumption made here is that the distributed
! {\tt GeneralGrid} {\tt iG} has been initialized with the same 
! coordinate system, sort order, other real attributes, and the same 
! indexing attributes for all processes on {\tt comm}.
!
! {\bf N.B.}:  Once the gridpoint data of the {\tt GeneralGrid} are assembled 
! on the {\tt root}, they are stored in the order determined by the input 
! {\tt GlobalSegMap} {\tt GSMap}.  The user may need to sorted these gathered
! data to order them in  accordance with the {\tt coordinate\_sort\_order} 
! attribute of {\tt iG}.
!
! {\bf N.B.}:  The output {\tt GeneralGrid} {\tt oG} represents allocated
! memory on the {\tt root}.  When the user no longer needs {\tt oG} it
! should be deallocated using {\tt GeneralGrid\_clean()} to avoid a memory
! leak
!
! !INTERFACE:

 subroutine GSM_gather_(iG, oG, GSMap, root, comm, stat)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90

      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_lsize => lsize
      use m_GlobalSegMap, only : GlobalSegMap_gsize => gsize

      use m_GeneralGrid, only : GeneralGrid
      use m_GeneralGrid, only : GeneralGrid_init => init
      use m_GeneralGrid, only : GeneralGrid_lsize => lsize

      implicit none

! !INPUT PARAMETERS: 
!
      type(GeneralGrid),  intent(in)  :: iG
      type(GlobalSegMap), intent(in)  :: GSMap
      integer,            intent(in)  :: root
      integer,            intent(in)  :: comm

! !OUTPUT PARAMETERS: 
!
      type(GeneralGrid),  intent(out) :: oG
      integer, optional,  intent(out) :: stat

! !REVISION HISTORY:
!       27Apr01 - J.W. Larson <larson@mcs.anl.gov> - API Specification.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GSM_gather_'

!Process ID
 integer :: myID
!Error flag
 integer :: ierr
!Number of points on the _Gathered_ grid:
 integer :: length

       ! Which process am I?

  call MPI_COMM_RANK(comm, myID, ierr)

  if(ierr /= 0) then
    call MP_perr(myname_,'MPI_COMM_RANK()',ierr)
    if(.not.present(stat)) call die(myname_)
    stat=ierr
    return
  endif

  if(myID == root) then ! prepare oG:

       ! The length of the _gathered_ GeneralGrid oG is determined by 
       ! the GlobalMap function GlobalSegMap_gsize()

     length = GlobalSegMap_gsize(GSMap)

       ! Initialize attributes of oG from iG, and length

     call GeneralGrid_init(oG, iG, length)

  endif

       ! Gather gridpoint data in iG%data to oG%data

  call AttrVect_Gather(oG%data, iG%data, GSMap, root, comm, ierr)

  if(ierr /= 0) then
    call MP_perr(myname_,'MPI_COMM_RANK()',ierr)
    if(.not.present(stat)) call die(myname_)
    stat=ierr
    return
  else
    if(present(stat)) stat=ierr
  endif

 end subroutine GSM_gather_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GM_scatter_ - scatter a GeneralGrid using input GlobalMap.
!
! !DESCRIPTION:  {\tt GM\_scatter\_()} takes an input {\tt GeneralGrid} 
! argument {\tt iG} (valid only on the {\tt root} process), and scatters 
! it to the distributed {\tt GeneralGrid} variable {\tt oG}.  The 
! {\tt GeneralGrid} {\tt oG} is distributed on the communicator 
! associated with the F90  handle {\tt comm} using the domain 
! decomposition described by the {\tt GlobalMap} argument {\tt GMap}.
! The success (failure) of this operation is reported as a zero (nonzero) 
! value in the optional {\tt INTEGER} output argument {\tt stat}.
!
! {\bf N.B.}:  The output {\tt GeneralGrid} {\tt oG} represents allocated
! memory on the {\tt root}.  When the user no longer needs {\tt oG} it
! should be deallocated using {\tt GeneralGrid\_clean()} to avoid a memory
! leak.
!
! !INTERFACE:

 subroutine GM_scatter_(iG, oG, GMap, root, comm, stat)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90

      use m_List, only : List_bcast => bcast
      use m_List, only : assignment(=)

      use m_GlobalMap, only : GlobalMap
      use m_GlobalMap, only : GlobalMap_lsize => lsize
      use m_GlobalMap, only : GlobalMap_gsize => gsize

      use m_GeneralGrid, only : GeneralGrid
      use m_GeneralGrid, only : GeneralGrid_init => init
      use m_GeneralGrid, only : GeneralGrid_lsize => lsize

      implicit none

! !INPUT PARAMETERS: 
!
      type(GeneralGrid), intent(in)  :: iG
      type(GlobalMap),   intent(in)  :: GMap
      integer,           intent(in)  :: root
      integer,           intent(in)  :: comm

! !OUTPUT PARAMETERS: 
!
      type(GeneralGrid), intent(out) :: oG
      integer, optional, intent(out) :: stat

! !REVISION HISTORY:
!       27Apr01 - J.W. Larson <larson@mcs.anl.gov> - API Specification.
!       04Jun01 - J.W. Larson <larson@mcs.anl.gov> - Changed comms model
!                 to MPI-style (i.e. iG valid on root only).
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GM_scatter_'

  logical :: DescendAssoc
  integer :: DescendSize
  integer :: ierr, myID

       ! Step 1.  Determine process ID number myID

  call MPI_COMM_RANK(comm, myID, ierr)
  if(ierr /= 0) then
     if(present(stat)) then
	write(stderr,*) myname_,':: Error calling MPI_COMM_RANK(comm...'
	stat = ierr
	return
     else
	call MP_perr_die(myname_,'MPI_COMM_RANK(comm...',ierr)
     endif
  endif

       ! Step 2.  On the root, initialize the List and LOGICAL 
       ! attributes of the GeneralGrid variable iG to oG.

  if(myID == root) then
     call copyGeneralGridHeader_(iG, oG)
     if(ierr /= 0) then
	if(present(stat)) then
	   write(stderr,*) myname_,':: Error calling copyGeneralGridHeader_(...'
	   stat = ierr
	   return
	else
	   call MP_perr_die(myname_,'copyGeneralGridHeader_(comm...',ierr)
	endif
     endif
  endif

       ! Step 3.  Broadcast from the root the List and LOGICAL 
       ! attributes of the GeneralGrid variable oG.

  call bcastGeneralGridHeader_(oG, root, comm, ierr)
  if(ierr /= 0) then
     if(present(stat)) then
	write(stderr,*) myname_,':: Error calling bcastGeneralGridHeader_(...'
	stat = ierr
	return
     else
	call MP_perr_die(myname_,'bcastGeneralGridHeader_(oG...',ierr)
     endif
  endif


       ! Step 4.  Using the GeneralMap GMap, scatter the AttrVect 
       ! portion of the input GeneralGrid iG to the GeneralGrid oG.

  call AttrVect_scatter(iG%data, oG%data, GMap, root, comm, ierr)
  if(ierr /= 0) then
     if(present(stat)) then
	write(stderr,*) myname_,':: Error calling AttrVect_scatter(iG%data...'
	stat = ierr
	return
     else
	call MP_perr_die(myname_,'call AttrVect_scatter(iG%data...',ierr)
     endif
  endif

       ! The GeneralGrid scatter is now complete.

 end subroutine GM_scatter_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GSM_scatter_ - scatter a GeneralGrid using input GlobalSegMap.
!
! !DESCRIPTION:  {\tt GM\_scatter\_()} takes an input {\tt GeneralGrid} 
! argument {\tt iG} (valid only on the {\tt root} process), and scatters 
! it to the distributed {\tt GeneralGrid} variable {\tt oG}.  The 
! {\tt GeneralGrid} {\tt oG} is distributed on the communicator 
! associated with the F90  handle {\tt comm} using the domain 
! decomposition described by the {\tt GlobalSegMap} argument {\tt GSMap}.
! The success (failure) of this operation is reported as a zero (nonzero) 
! value in the optional {\tt INTEGER} output argument {\tt stat}.
!
! {\bf N.B.}:  The output {\tt GeneralGrid} {\tt oG} represents allocated
! memory on the {\tt root}.  When the user no longer needs {\tt oG} it
! should be deallocated using {\tt GeneralGrid\_clean()} to avoid a memory
! leak.
!
! !INTERFACE:

 subroutine GSM_scatter_(iG, oG, GSMap, root, comm, stat)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90

      use m_List, only : List_bcast => bcast
      use m_List, only : assignment(=)

      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_lsize => lsize
      use m_GlobalSegMap, only : GlobalSegMap_gsize => gsize

      use m_GeneralGrid, only : GeneralGrid
      use m_GeneralGrid, only : GeneralGrid_init => init
      use m_GeneralGrid, only : GeneralGrid_lsize => lsize

      implicit none

! !INPUT PARAMETERS: 
!
      type(GeneralGrid),  intent(in)  :: iG
      type(GlobalSegMap), intent(in)  :: GSMap
      integer,            intent(in)  :: root
      integer,            intent(in)  :: comm

! !OUTPUT PARAMETERS: 
!
      type(GeneralGrid),  intent(out) :: oG
      integer, optional,  intent(out) :: stat

! !REVISION HISTORY:
!       27Apr01 - J.W. Larson <larson@mcs.anl.gov> - API Specification.
!       04Jun01 - J.W. Larson <larson@mcs.anl.gov> - Initial code.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GSM_scatter_'

  integer :: ierr, myID

       ! Step 1.  Determine process ID number myID

  call MPI_COMM_RANK(comm, myID, ierr)
  if(ierr /= 0) then
     if(present(stat)) then
	write(stderr,*) myname_,':: Error calling MPI_COMM_RANK(comm...'
	stat = ierr
	return
     else
	call MP_perr_die(myname_,'MPI_COMM_RANK(comm...',ierr)
     endif
  endif

       ! Step 2.  On the root, initialize the List and LOGICAL 
       ! attributes of the GeneralGrid variable iG to oG.

  if(myID == root) then
     call copyGeneralGridHeader_(iG, oG)
     if(ierr /= 0) then
	if(present(stat)) then
	   write(stderr,*) myname_,':: Error calling copyGeneralGridHeader_(...'
	   stat = ierr
	   return
	else
	   call MP_perr_die(myname_,'copyGeneralGridHeader_(comm...',ierr)
	endif
     endif
  endif

       ! Step 3.  Broadcast from the root the List and LOGICAL 
       ! attributes of the GeneralGrid variable oG.

  call bcastGeneralGridHeader_(oG, root, comm, ierr)
  if(ierr /= 0) then
     if(present(stat)) then
	write(stderr,*) myname_,':: Error calling bcastGeneralGridHeader_(...'
	stat = ierr
	return
     else
	call MP_perr_die(myname_,'bcastGeneralGridHeader_(oG...',ierr)
     endif
  endif

       ! Step 4.  Using the GeneralSegMap GSMap, scatter the AttrVect 
       ! portion of the input GeneralGrid iG to the GeneralGrid oG.

  call AttrVect_scatter(iG%data, oG%data, GSMap, root, comm, ierr)
  if(ierr /= 0) then
     if(present(stat)) then
	write(stderr,*) myname_,':: Error calling AttrVect_scatter(iG%data...'
	stat = ierr
	return
     else
	call MP_perr_die(myname_,'call AttrVect_scatter(iG%data...',ierr)
     endif
  endif

       ! The GeneralGrid scatter is now complete.

 end subroutine GSM_scatter_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: bcast_ - Broadcast a GeneralGrid.
!
! !DESCRIPTION:  {\tt bcast\_()} takes an input {\tt GeneralGrid} 
! argument {\tt ioG} (valid only on the {\tt root} process), and 
! broadcasts it to all processes on the communicator associated with the 
! F90  handle {\tt comm}.  The success (failure) of this operation is 
! reported as a zero (nonzero) value in the optional {\tt INTEGER} 
! output argument {\tt stat}.
!
! {\bf N.B.}:  On the non-root processes, the output {\tt GeneralGrid} 
! {\tt ioG} represents allocated memory.  When the user no longer needs 
! {\tt ioG} it should be deallocated by invoking {\tt GeneralGrid\_clean()}. 
! Failure to do so risks a memory leak.
!
! !INTERFACE:

 subroutine bcast_(ioG, root, comm, stat)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90

      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_lsize => lsize
      use m_GlobalSegMap, only : GlobalSegMap_gsize => gsize

      use m_GeneralGrid, only : GeneralGrid
      use m_GeneralGrid, only : GeneralGrid_init => init
      use m_GeneralGrid, only : GeneralGrid_lsize => lsize

      implicit none

! !INPUT PARAMETERS: 
!
      integer,           intent(in)    :: root
      integer,           intent(in)    :: comm

! !INPUT/OUTPUT PARAMETERS: 
!
      type(GeneralGrid), intent(inout) :: ioG

! !OUTPUT PARAMETERS: 
!
      integer, optional, intent(out)   :: stat

! !REVISION HISTORY:
!       27Apr01 - J.W. Larson <larson@mcs.anl.gov> - API Specification.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::bcast_'

  integer :: ierr, myID

       ! Step 1.  Determine process ID number myID

  call MPI_COMM_RANK(comm, myID, ierr)
  if(ierr /= 0) then
     if(present(stat)) then
	write(stderr,*) myname_,':: Error calling MPI_COMM_RANK(comm...'
	stat = ierr
	return
     else
	call MP_perr_die(myname_,'MPI_COMM_RANK(comm...',ierr)
     endif
  endif

       ! Step 2.  Broadcast from the root the List and LOGICAL 
       ! attributes of the GeneralGrid variable ioG.

  call bcastGeneralGridHeader_(ioG, root, comm, ierr)
  if(ierr /= 0) then
     if(present(stat)) then
	write(stderr,*) myname_,':: Error calling bcastGeneralGridHeader_(...'
	stat = ierr
	return
     else
	call MP_perr_die(myname_,'bcastGeneralGridHeader_(oG...',ierr)
     endif
  endif

       ! Step 3.  Broadcast ioG%data from the root.

  call AttrVect_bcast(ioG%data, root, comm, ierr)
  if(ierr /= 0) then
     if(present(stat)) then
	write(stderr,*) myname_,':: Error calling AttrVect_scatter(iG%data...'
	stat = ierr
	return
     else
	call MP_perr_die(myname_,'call AttrVect_scatter(iG%data...',ierr)
     endif
  endif

       ! The GeneralGrid broadcast is now complete.

 end subroutine bcast_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: bcastGeneralGridHeader_ - Broadcast the GeneralGrid Header.
!
! !DESCRIPTION:  This routine broadcasts the header information from 
! the input {\tt GeneralGrid} argument {\tt ioGGrid} (on input valid 
! on the {\tt root} only).  This broadcast is from the {\tt root} to 
! all processes on the communicator associated with the fortran 90 
! {\tt INTEGER} handle {\tt comm}.  The success (failure) of this operation 
! corresponds to a zero (nonzero) value for the output {\tt INTEGER} flag 
! {\tt stat}. 
!
! The {\em header information} in a {\tt GeneralGrid} variable comprises 
! all the non-{\tt AttrVect} components of the {\tt GeneralGrid}; that 
! is, everything except the gridpoint coordinate, geometry, and index 
! data stored in {\tt iGGrid\%data}.  This information includes:
! \begin{enumerate}
! \item The coordinates in {\tt iGGrid\%coordinate\_list}
! \item The coordinate sort order in {\tt iGGrid\%coordinate\_sort\_order}
! \item The area/volume weights in {\tt iGGrid\%weight\_list}
! \item Other {\tt REAL} geometric information in {\tt iGGrid\%other\_list}
! \item Indexing information in {\tt iGGrid\%index\_list}
! \item The {\tt LOGICAL} descending/ascending order sort flags in 
! {\tt iGGrid\%descend(:)}.
! \end{enumerate}
!
! !INTERFACE:

 subroutine bcastGeneralGridHeader_(ioGGrid, root, comm, stat)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90

      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_lsize => lsize
      use m_GlobalSegMap, only : GlobalSegMap_gsize => gsize

      use m_GeneralGrid, only : GeneralGrid
      use m_GeneralGrid, only : GeneralGrid_init => init
      use m_GeneralGrid, only : GeneralGrid_lsize => lsize

      implicit none

! !INPUT PARAMETERS: 
!
      integer,           intent(in)    :: root
      integer,           intent(in)    :: comm

! !INPUT/OUTPUT PARAMETERS: 
!
      type(GeneralGrid), intent(inout) :: ioGGrid

! !OUTPUT PARAMETERS: 
!
      integer, optional, intent(out)   :: stat

! !REVISION HISTORY:
!       05Jun01 - J.W. Larson <larson@mcs.anl.gov> - Initial code.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::bcastGeneralGridHeader_'
! Process ID
  integer :: myID
! Error flag
  integer :: ierr
! Logical flag--is ioGGrid%descend associated?
  logical :: DescendAssoc
! Size of array ioGGrid%descend(:)
  integer :: DescendSize

       ! Determine process ID number myID

  call MPI_COMM_RANK(comm, myID, ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,'MPI_COMM_RANK(comm...',ierr)
  endif

       ! Step 1. Broadcast List attributes of the GeneralGrid.

  call List_bcast(ioGGrid%coordinate_list, root, comm, ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,'List_bcast(ioGGrid%coordinate_list...',ierr)
  endif

  call List_bcast(ioGGrid%coordinate_sort_order, root, comm, ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,'List_bcast(ioGGrid%coordinate_sort_order...', &
	              ierr)
  endif

  call List_bcast(ioGGrid%weight_list, root, comm, ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,'List_bcast(ioGGrid%weight_list...',ierr)
  endif

  call List_bcast(ioGGrid%other_list, root, comm, ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,'List_bcast(ioGGrid%other_list...',ierr)
  endif

  call List_bcast(ioGGrid%index_list, root, comm, ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,'List_bcast(ioGGrid%index_list...',ierr)
  endif

       ! Step 2. If applicable Broadcast ioGGrid%descend from the root.

  if(myID == root) then
     DescendAssoc = associated(ioGGrid%descend)
  endif

  call MPI_BCAST(DescendAssoc, 1, MP_LOGICAL, root, comm, ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,'MPI_BCAST(DescendAssoc...',ierr)
  endif

       ! If ioGGrid%descend is associated on the root, prepare and
       ! execute its broadcast

  if(DescendAssoc) then

       ! On the root, get the size of ioGGrid%descend(:)

     if(myID == root) then
	DescendSize = size(ioGGrid%descend)
     endif

       ! Broadcast the size of ioGGrid%descend(:) from the root.

     call MPI_BCAST(DescendSize, 1, MP_INTEGER, root, comm, ierr)
     if(ierr /= 0) then
	call MP_perr_die(myname_,'MPI_BCAST(DescendSize...',ierr)
     endif

       ! Off the root, allocate ioGGrid%descend(:)

     if(myID /= root) then
	allocate(ioGGrid%descend(DescendSize), stat=ierr)
	if(ierr /= 0) then
	   call MP_perr_die(myname_,'allocate(ioGGrid%descend(...',ierr)
	endif
     endif

       ! Finally, broadcast ioGGrid%descend(:) from the root

     call MPI_BCAST(ioGGrid%descend, DescendSize, MP_LOGICAL, root, &
                    comm, ierr)
     if(ierr /= 0) then
	call MP_perr_die(myname_,'MPI_BCAST(ioGGrid%descend...',ierr)
     endif

  endif

       ! The broadcast of the GeneralGrid Header from the &
       ! root is complete.

 end subroutine bcastGeneralGridHeader_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: copyGeneralGridHeader_ - Copy the GeneralGrid Header.
!
! !DESCRIPTION:  This routine copies the header information from the 
! input {\tt GeneralGrid} argument {\tt iGGrid} to the output 
! {\tt GeneralGrid} argument {\tt oGGrid}.  The {\em header information} 
! in a {\tt GeneralGrid} variable comprises all the non-{\tt AttrVect} 
! components of the {\tt GeneralGrid}; that is, everything except the 
! gridpoint coordinate, geometry, and index data stored in 
! {\tt iGGrid\%data}.  This information includes:
! \begin{enumerate}
! \item The coordinates in {\tt iGGrid\%coordinate\_list}
! \item The coordinate sort order in {\tt iGGrid\%coordinate\_sort\_order}
! \item The area/volume weights in {\tt iGGrid\%weight\_list}
! \item Other {\tt REAL} geometric information in {\tt iGGrid\%other\_list}
! \item Indexing information in {\tt iGGrid\%index\_list}
! \item The {\tt LOGICAL} descending/ascending order sort flags in 
! {\tt iGGrid\%descend(:)}.
! \end{enumerate}
!
! !INTERFACE:

 subroutine copyGeneralGridHeader_(iGGrid, oGGrid)
!
! !USES:
!
      use m_stdio
      use m_die

      use m_List, only : List_nitem => nitem

      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_lsize => lsize
      use m_GlobalSegMap, only : GlobalSegMap_gsize => gsize

      use m_GeneralGrid, only : GeneralGrid
      use m_GeneralGrid, only : GeneralGrid_init => init
      use m_GeneralGrid, only : GeneralGrid_lsize => lsize

      implicit none

! !INPUT PARAMETERS: 
!
      type(GeneralGrid), intent(in)  :: iGGrid

! !OUTPUT PARAMETERS: 
!
      type(GeneralGrid), intent(out) :: oGGrid

! !REVISION HISTORY:
!       05Jun01 - J.W. Larson <larson@mcs.anl.gov> - Initial code.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::copyGeneralGridHeader_'

  logical :: DescendAssoc
  integer :: DescendSize, i, ierr

       ! Step 1. Copy GeneralGrid List attributes from iGGrid 
       ! to oGGrid.

  oGGrid%coordinate_list = iGGrid%coordinate_list
  oGGrid%coordinate_sort_order = iGGrid%coordinate_sort_order
  oGGrid%weight_list = iGGrid%weight_list
  oGGrid%other_list = iGGrid%other_list
  oGGrid%index_list = iGGrid%index_list

       ! Step 2. If necessary, copy over the LOGICAL flags in
       ! iGGrid%descend(:).

  DescendAssoc = .FALSE.

       ! Is iGGrid%descend associated?  If so, does the size of 
       ! iGGrid%descend match the number of coordinates?

  DescendAssoc = associated(iGGrid%descend)
  if(DescendAssoc) then
     DescendSize = size(iGGrid%descend)
     if((DescendSize) /= List_nitem(iGGrid%coordinate_list)) then
	write(stderr,*) myname,':: size(iGGrid%descend)/coordinate count mismatch'
	write(stderr,*) myname_,':: size(iGGrid%descend) = ',size(iGGrid%descend)
	write(stderr,*) myname_,':: List_nitem(iGGrid%coordinate_list = ', &
	     List_nitem(iGGrid%coordinate_list)
	call MP_perr_die(myname_,'size(iGGrid%descend)/coordinate count mismatch', &
	     size(iGGrid%descend)-List_nitem(iGGrid%coordinate_list))
     endif
  endif

       ! If iGGrid%descend is associated, copy its entries to 
       ! oGGrid%descend.

  if(DescendAssoc) then
     allocate(oGGrid%descend(DescendSize), stat=ierr)
     if(ierr /= 0) then 
        call MP_perr_die(myname_,'allocate(iGGrid%descend(...',ierr)
     endif
     do i=1,DescendSize
	oGGrid%descend(i) = iGGrid%descend(i)
     end do
  endif

       ! The GeneralGrid header copy is now complete.

 end subroutine copyGeneralGridHeader_

 end module m_GeneralGridComms








