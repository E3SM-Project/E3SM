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
! !IROUTINE: GM_gather_ - gather a GeneralGrid using input GlobalMap.
!
! !DESCRIPTION:  {\tt GM\_gather\_()} takes an input {\tt GeneralGrid} 
! argument {\tt iG} whose decomposition on the communicator associated 
! with the F90 handle {\tt comm} is described by the {\tt GlobalMap} 
! argument {\tt GMap}, and gathers it to the {\tt GeneralGrid} output
! argument {\tt oG} on the {\tt root}.  
!
! {\bf N.B.}:  An important assumption made here is that the distribute 
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
! !DESCRIPTION:
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
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine GM_scatter_(iG, oG, GMap, root, comm, stat)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90

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
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GM_scatter_'

  character(len=1), dimension(:), pointer :: coordinate_list
  character(len=1), dimension(:), pointer :: coordinate_sort_order
  character(len=1), dimension(:), pointer :: weight_list
  character(len=1), dimension(:), pointer :: other_list
  character(len=1), dimension(:), pointer :: index_list
  integer :: nchars

  integer :: ierr, myID

       ! Determine process ID number myID

  call MPI_COMM_RANK(comm, myID, ierr)

  if(myID == root) nchars = size(iG%coordinate_list%bf)

  call MPI_BCAST(nchars, 1, MP_INTEGER, root, comm, ierr)

  if(myID == root) nchars = size(iG%coordinate_sort_order%bf)

  call AttrVect_scatter(iG, oG, Gmap, root, comm, ierr)

 end subroutine GM_scatter_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GSM_scatter_ - scatter a GeneralGrid using input GlobalSegMap.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine GSM_scatter_(iG, oG, GSMap, root, comm, stat)
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

  character(len=*),parameter :: myname_=myname//'::GSM_scatter_'

 end subroutine GSM_scatter_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: bcast_ - Broadcast a GeneralGrid.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine bcast_(iG, root, comm, stat)
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
      type(GeneralGrid), intent(inout) :: iG

! !OUTPUT PARAMETERS: 
!
      integer, optional, intent(out)   :: stat

! !REVISION HISTORY:
!       27Apr01 - J.W. Larson <larson@mcs.anl.gov> - API Specification.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::bcast_'

 end subroutine bcast_

 end module m_GeneralGridComms








