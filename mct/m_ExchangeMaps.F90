!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$ 
!BOP -------------------------------------------------------------------
!
! !MODULE: m_ExchangeMaps - Exchange of Global Mapping Objects.
!
! !DESCRIPTION:
! This module contains tools for exchanging between two communicators 
! comm_A and comm_B the global mapping index objects on the two 
! communicators.
!
! !INTERFACE:

 module m_ExchangeMaps

      implicit none

      private   ! except
!
! !PUBLIC MEMBER FUNCTIONS:
!
      public :: MapHandshake
      public :: LoadMapPars
      public :: ExchangeMap

    interface LoadMapPars ; module procedure   &
	 LoadGlobalMapPars_, &
	 LoadGlobalSegMapPars_
    end interface

    interface MapHandshake ; module procedure   &
	 MapHandshake_
    end interface

    interface ExchangeMap ; module procedure   &
        ExGMapGMap_,   &  ! GlobalMap for GlobalMap
        ExGSMapGSMap_, &  ! GlobalSegMap for GlobalSegMap
        ExGMapGSMap_,  &  ! GlobalMap for GlobalSegMap
        ExGSMapGMap_      ! GlobalSegMap for GlobalMap
    end interface

! !PUBLIC DATA MEMBERS:

                                ! Map handshaking is implemented as
                                ! the exchange of an INTEGER array of
                                ! flags / values.  These are defined
                                ! below:

    public :: NumHandshakePars  ! Number of handshake parameters; i.e.
                                ! size of exhcanged parameters array
    public :: ComponentIDIndex  ! location of component ID number
    public :: MapTypeIndex      ! storage location of map type flag
    public :: NumMapTypes       ! number of defined map types
    public :: GlobalMapFlag     ! value of map type flag signifying
                                ! GlobalMap
    public :: GlobalSegMapFlag  ! value of map type flag signifying
                                ! GlobalSegMap
    public :: GsizeIndex        ! storage location of gsize for map
    public :: NumSegIndex       ! storage location of the number of
                                ! segments in the map:  Number of 
                                ! processes for a GlobalMap; Number of
                                ! segments for GlobalSegMap
! !REVISION HISTORY:
!       03Feb01 - J.W. Larson <larson@mcs.anl.gov> - initial module
!       03Aug01 - E.T. Ong <eong@mcs.anl.gov> - in ExGSMapGSMap,
!                 call GlobalSegMap_init with actual shaped arrays
!                 for non-root processes to satisfy Fortran 90 standard.
!                 See comments in subroutine.
!       15Feb02 - R. Jacob <jacob@mcs.anl.gov> - use MCT_comm instead of
!		  MP_COMM_WORLD
!EOP ___________________________________________________________________


  character(len=*),parameter :: myname='m_ExchangeMaps'

! Map Handshaking Parameters:  Map handshaking occurs via 
! exchange of an array of INTEGER flags.

  ! Number of Handshaking Parameters

  integer, parameter :: NumHandshakePars = 4

  ! ComponentIDIndex defines the storage location of the  flag 
  ! signifying the component number in MCTWorld

  integer, parameter :: ComponentIDIndex = 1

  ! MapTypeIndex defines the storage location in the handshake array 
  ! of the type of map offered for exchange

  integer, parameter :: MapTypeIndex = 2

          ! NumMapTypes is the number of legitimate MapTypeIndex Values:

  integer, parameter :: NumMapTypes = 2

          ! Recognized MapTypeIndex Values:

  integer, parameter :: GlobalMapFlag = 1
  integer, parameter :: GlobalSegMapFlag = 2

  ! GsizeIndex defines the location of the grid size (number of points)
  ! for the map.  This size is 

  integer, parameter :: GsizeIndex = 3

  ! NumSegIndex defines the location of the number of segments in the
  ! map.  For a GlobalMap, this is the number of processes in the map.
  ! For a GlobalSegMap, this is the number of global segments (ngseg).

  integer, parameter :: NumSegIndex = 4

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: MapHandshake_ - Exchange Map descriptors.
!
! !DESCRIPTION:
! This routine takes input Map descriptors stored in the {\tt INTEGER} 
! array {\tt LocalMapPars}, the local communicator on which this map is
! defined ({\tt LocalComm}), and the remote component ID 
! {\tt RemoteCompID}, and effects an exchange of map descriptors with
! the remote component, which are returned in the {\tt INTEGER} array
! {\tt RemoteMapPars}.
!
! {\bf N.B.: } The values present in {\tt LocalMapPars} need to be valid
! only on the root of {\tt LocalComm}.  Likewise, the returned values in
! {\tt RemoteMapPars} will be valid on the root of {\tt LocalComm}.
!
! !INTERFACE:

 subroutine MapHandshake_(LocalMapPars, LocalComm, RemoteCompID, &
                          RemoteMapPars)

!
! !USES:
!
      use m_mpif90
      use m_die,      only : MP_perr_die
      use m_stdio
      use m_MCTWorld, only : ThisMCTWorld
      use m_MCTWorld, only : ComponentRootRank

      implicit none
!
! !INPUT PARAMETERS: 
!
      integer, intent(in)  :: LocalMapPars(NumHandshakePars)
      integer, intent(in)  :: LocalComm
      integer, intent(in)  :: RemoteCompID
!
! !OUTPUT PARAMETERS:
!
      integer, intent(out) :: RemoteMapPars(NumHandshakePars)

! !REVISION HISTORY:
!       06Feb01 - J.W. Larson <larson@mcs.anl.gov> - API specification.
!       20Apr01 - R.L. Jacob  <jacob@mcs.anl.gov> - add status argument
!                 to MPI_RECV
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::MapHandshake_'

  integer :: ierr, myID, RemoteRootID, SendTag, RecvTag
  integer,dimension(MP_STATUS_SIZE) :: status

  call MP_COMM_RANK(LocalComm, myID, ierr)
  if(ierr /= 0) call MP_perr_die(myname_,'call MP_COMM_RANK()',ierr)

  RemoteRootID = ComponentRootRank(RemoteCompID, ThisMCTWorld, check=.true.)

  if(myID == 0) then ! I am the root on LocalComm

      ! Compute send/receive tags:

     SendTag = 10 * LocalMapPars(ComponentIDIndex) + RemoteCompID
     RecvTag = LocalMapPars(ComponentIDIndex) + 10 * RemoteCompID

      ! Post send to RemoteRootID:

     call MPI_SEND(LocalMapPars, NumHandshakePars, MP_INTEGER, &
	           RemoteRootID, SendTag, ThisMCTWorld%MCT_comm, ierr)
     if(ierr /= 0) call MP_perr_die(myname_,'call MPI_SEND()',ierr)

      ! Post receive from RemoteRootID:

     call MPI_RECV(RemoteMapPars, NumHandshakePars, MP_INTEGER, &
	           RemoteRootID, RecvTag, ThisMCTWorld%MCT_comm, status, ierr)
     if(ierr /= 0) call MP_perr_die(myname_,'call MPI_RECV()',ierr)

  endif ! if(myID == 0)

 end subroutine MapHandshake_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: LoadGlobalMapPars_ - Load GlobalMap descriptors.
!
! !DESCRIPTION:
! This routine takes an input {\tt GlobalMap} variable {\tt Gmap}, and
! loads its descriptors the output {\tt INTEGER} array {\tt MapPars}.
! The dimensions of this array, and loading order are all defined in
! the declaration section of this module.
!
! !INTERFACE:

 subroutine LoadGlobalMapPars_(GMap, MapPars)

!
! !USES:
!
      use m_mpif90
      use m_die
      use m_stdio
      use m_GlobalMap, only : GlobalMap
      use m_GlobalMap, only : GlobalMap_comp_id => comp_id
      use m_GlobalMap, only : GlobalMap_gsize => gsize
!      use m_GlobalMap, only : GlobalMap_nprocs => nprocs

      implicit none
!
! !INPUT PARAMETERS: 
!
      type(GlobalMap), intent(in)  :: GMap
!
! !OUTPUT PARAMETERS: 
!
      integer,         intent(out) :: MapPars(NumHandshakePars)

! !REVISION HISTORY:
!       06Feb01 - J.W. Larson <larson@mcs.anl.gov> - Initial version.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::LoadGlobalMapPars_'

  MapPars(ComponentIDIndex) = GlobalMap_comp_id(GMap)
  MapPars(MapTypeIndex) = GlobalMapFlag
  MapPars(GsizeIndex) = GlobalMap_gsize(GMap)
!  MapPars(NumSegIndex) = GlobalMap_nprocs(GSMap)

 end subroutine LoadGlobalMapPars_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: LoadGlobalSegMapPars_ - Load GlobalSegMap descriptors.
!
! !DESCRIPTION:
! This routine takes an input {\tt GlobalSegMap} variable {\tt Gmap}, and
! loads its descriptors the output {\tt INTEGER} array {\tt MapPars}.
! The dimensions of this array, and loading order are all defined in
! the declaration section of this module.
!
! !INTERFACE:

 subroutine LoadGlobalSegMapPars_(GSMap, MapPars)

!
! !USES:
!
      use m_mpif90
      use m_die
      use m_stdio
      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_comp_id => comp_id
      use m_GlobalSegMap, only : GlobalSegMap_gsize => gsize
      use m_GlobalSegMap, only : GlobalSegMap_ngseg => ngseg


      implicit none
!
! !INPUT PARAMETERS: 
!
      type(GlobalSegMap), intent(in)  :: GSMap
!
! !OUTPUT PARAMETERS: 
!
      integer,            intent(out) :: MapPars(NumHandshakePars)

! !REVISION HISTORY:
!       06Feb01 - J.W. Larson <larson@mcs.anl.gov> - Initial version.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::LoadGlobalSegMapPars_'

  MapPars(ComponentIDIndex) = GlobalSegMap_comp_id(GSMap)
  MapPars(MapTypeIndex) = GlobalSegMapFlag
  MapPars(GsizeIndex) = GlobalSegMap_gsize(GSMap)
  MapPars(NumSegIndex) = GlobalSegMap_ngseg(GSMap)

 end subroutine LoadGlobalSegMapPars_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ExGMapGMap_ - Trade of GlobalMap structures
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine ExGMapGMap_(LocalGMap, LocalComm, RemoteGMap, &
                        Remote_comp_id, ierr) 

!
! !USES:
!
      use m_mpif90
      use m_die
      use m_stdio
      use m_GlobalMap, only : GlobalMap
      use m_GlobalMap, only : GlobalMap_init => init

      implicit none

      type(GlobalMap), intent(in)  :: LocalGMap      ! Local GlobalMap
      integer,         intent(in)  :: LocalComm      ! Local Communicator
      type(GlobalMap), intent(out) :: RemoteGMap     ! Remote GlobalMap
      integer        , intent(in)  :: Remote_comp_id ! Remote component id
      integer        , intent(out) :: ierr           ! Error Flag

! !REVISION HISTORY:
!       03Feb01 - J.W. Larson <larson@mcs.anl.gov> - API specification.
!       20Apr01 - R.L. Jacob  <larson@mcs.anl.gov> - bug fix.  Use
!                 NumSegIndex instead of GsizeIndex for start,length,
!                 pe_loc arrays.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ExGMapGMap_'

 end subroutine ExGMapGMap_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ExGSMapGSMap_ - Trade of GlobalSegMap structures.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine ExGSMapGSMap_(LocalGSMap, LocalComm, RemoteGSMap, &
                          RemoteCompID, ierr) 

!
! !USES:
!
      use m_mpif90
      use m_die
      use m_stdio
      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_init => init

      use m_MCTWorld, only : ThisMCTWorld
      use m_MCTWorld, only : ComponentRootRank

      implicit none

      type(GlobalSegMap), intent(in)  :: LocalGSMap  ! Local GlobalSegMap
      integer,            intent(in)  :: LocalComm   ! Local Communicator
      type(GlobalSegMap), intent(out) :: RemoteGSMap ! Remote GlobalSegMap
      integer        , intent(in)  :: RemoteCompID   ! Remote component id
      integer        , intent(out) :: ierr           ! Error Flag

! !REVISION HISTORY:
!       03Feb01 - J.W. Larson <larson@mcs.anl.gov> - API specification.
!       07Feb01 - J.W. Larson <larson@mcs.anl.gov> - First full version.
!       20Apr01 - R.L. Jacob  <jacob@mcs.anl.gov> - add status argument
!                 to MPI_RECV
!       25Apr01 - R.L. Jacob  <jacob@mcs.anl.gov> - set SendTag and
!                 RecvTag values
!       03May01 - R.L. Jacob <jacob@mcs.anl.gov> - change MPI_SEND to
!                 MPI_ISEND to avoid possible buffering problems seen
!                 on IBM SP.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ExGSMapGSMap_'

! root ID on local communicator:
  integer, parameter :: root = 0
! Storage for local and remote map descriptors:
  integer :: LocalMapPars(NumHandshakePars)
  integer :: RemoteMapPars(NumHandshakePars)
! Send and Receive Buffers
  integer, dimension(:), allocatable :: SendBuf
  integer, dimension(:), allocatable :: RecvBuf
! Send and Receive Tags
  integer :: SendTag, RecvTag
! Storage arrays for Remote GlobalSegMap data:
  integer, dimension(:), allocatable :: start, length, pe_loc

  integer :: myID, ngseg, remote_root,req
  integer :: local_ngseg, remote_ngseg
  integer,dimension(MP_STATUS_SIZE) :: status,wstatus

      ! Determine rank on local communicator:

  call MP_COMM_RANK(LocalComm, myID, ierr)
  if(ierr /= 0) call MP_perr_die(myname_,'call MP_COMM_RANK()',ierr)

      ! If the root, exchange map handshake descriptors,
      ! and information needed to initialize the remote map 
      ! on the local communicator.

  if(myID == root) then

     call LoadGlobalSegMapPars_(LocalGSMap, LocalMapPars)

     call MapHandshake_(LocalMapPars, LocalComm, RemoteCompID, &
                        RemoteMapPars)

      ! Consistency Checks between LocalMapPars and RemoteMapPars:

     if(LocalMapPars(MapTypeIndex) /= RemoteMapPars(MapTypeIndex)) then
	ierr = 2
	call die(myname_,'Map Type mismatch',ierr)
     endif

     if(LocalMapPars(GsizeIndex) /= RemoteMapPars(GsizeIndex)) then
	ierr = 3
	call die(myname_,'Map Grid Size mismatch',ierr)
     endif

     if(RemoteCompID /= RemoteMapPars(ComponentIDIndex)) then
	ierr = 4
	call die(myname_,'Component ID mismatch',ierr)
     endif

      ! SendBuf will hold the arrays LocalGSMap%start, LocalGSMap%length,
      ! and LocalGSMap%pe_loc in that order.

     allocate(SendBuf(3*LocalMapPars(NumSegIndex)), stat=ierr)
     if(ierr /= 0) call die(myname_,'allocate(SendBuf...)',ierr)

      ! RecvBuf will hold the arrays RemoteGSMap%start, RemoteGSMap%length,
      ! and RemoteGSMap%pe_loc in that order.

     allocate(RecvBuf(3*RemoteMapPars(NumSegIndex)), stat=ierr)
     if(ierr /= 0) call die(myname_,'allocate(RecvBuf...)',ierr)

      ! Load SendBuf in the order described above:
     local_ngseg = LocalMapPars(NumSegIndex)
     SendBuf(1:local_ngseg) = &
	                  LocalGSMap%start(1:local_ngseg)
     SendBuf(local_ngseg+1:2*local_ngseg) = &
	                  LocalGSMap%length(1:local_ngseg)
     SendBuf(2*local_ngseg+1:3*local_ngseg) = &
	                  LocalGSMap%pe_loc(1:local_ngseg)

      ! Determine the remote component root:

     remote_root = ComponentRootRank(RemoteMapPars(ComponentIDIndex), &
	                             ThisMCTWorld, check = .true.)

     SendTag = 10 * LocalMapPars(ComponentIDIndex) + RemoteCompID
     RecvTag = LocalMapPars(ComponentIDIndex) + 10 * RemoteCompID

      ! Send off SendBuf to the remote component root:

     call MPI_ISEND(SendBuf(1), 3*LocalMapPars(NumSegIndex), MP_INTEGER, &
	           remote_root, SendTag, ThisMCTWorld%MCT_comm, req, ierr)
     if(ierr /= 0) call MP_perr_die(myname_,'MPI_SEND(SendBuf...',ierr)

      ! Receive RecvBuf from the remote component root:

     call MPI_RECV(RecvBuf, 3*RemoteMapPars(NumSegIndex), MP_INTEGER, &
	           remote_root, RecvTag, ThisMCTWorld%MCT_comm, status, ierr)
     if(ierr /= 0) call MP_perr_die(myname_,'MPI_Recv(RecvBuf...',ierr)

     call MPI_WAIT(req,wstatus,ierr)
     if(ierr /= 0) call MP_perr_die(myname_,'MPI_WAIT(SendBuf..',ierr)

      ! Allocate arrays start(:), length(:), and pe_loc(:)

     allocate(start(RemoteMapPars(NumSegIndex)),  &
	      length(RemoteMapPars(NumSegIndex)), &
	      pe_loc(RemoteMapPars(NumSegIndex)), stat=ierr)
     if(ierr /= 0) call die(myname_,'allocate(start...',ierr)

      ! Unpack RecvBuf into arrays start(:), length(:), and pe_loc(:)
     remote_ngseg = RemoteMapPars(NumSegIndex)
     start(1:remote_ngseg) = RecvBuf(1:remote_ngseg)
     length(1:remote_ngseg) = &
                        RecvBuf(remote_ngseg+1:2*remote_ngseg)
     pe_loc(1:remote_ngseg) = &
                        RecvBuf(2*remote_ngseg+1:3*remote_ngseg)

  endif ! if(myID == root)

        ! Non-root processes call GlobalSegMap_init with start, 
        ! length, and pe_loc, although these arguments are 
        ! not used in the subroutine. Since these correspond to dummy 
        ! shaped array arguments in GlobalSegMap_init, the Fortran 90 
        ! standard dictates that the actual arguments must contain 
        ! complete shape information. Therefore, these array arguments 
        ! must be allocated on all processes.

  if(myID /= root) then
     
      allocate(start(1), length(1), pe_loc(1), stat=ierr)
      if(ierr /= 0) call die(myname_,'non-root allocate(start...',ierr)

  endif
  

      ! Initialize the Remote GlobalSegMap RemoteGSMap

  call GlobalSegMap_init(RemoteGSMap, RemoteMapPars(NumSegIndex), &
                         start, length, pe_loc, root, LocalComm,  &
			 RemoteCompID, RemoteMapPars(GsizeIndex))


      ! Deallocate allocated arrays

  deallocate(start, length, pe_loc, stat=ierr)
  if(ierr /= 0) then
     call die(myname_,'deallocate(start...',ierr)
  endif

      ! Deallocate allocated arrays on the root:

  if(myID == root) then
     
     deallocate(SendBuf, RecvBuf, stat=ierr)
     if(ierr /= 0) then
	call die(myname_,'deallocate(SendBuf...',ierr)
     endif

  endif ! if(myID == root)

 end subroutine ExGSMapGSMap_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ExGMapGSMap_ - Trade of GlobalMap for GlobalSegMap.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine ExGMapGSMap_(LocalGMap, LocalComm, RemoteGSMap, &
                         Remote_comp_id, ierr) 

!
! !USES:
!
      use m_mpif90
      use m_die
      use m_stdio

      use m_GlobalMap, only : GlobalMap
      use m_GlobalMap, only : GlobalMap_init => init
      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_init => init

      implicit none

      type(GlobalMap),    intent(in)  :: LocalGMap   ! Local GlobalMap
      integer,            intent(in)  :: LocalComm   ! Local Communicator
      type(GlobalSegMap), intent(out) :: RemoteGSMap ! Remote GlobalSegMap
      integer        , intent(in)  :: Remote_comp_id ! Remote component id
      integer        , intent(out) :: ierr           ! Error Flag

! !REVISION HISTORY:
!       03Feb01 - J.W. Larson <larson@mcs.anl.gov> - API specification.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ExGMapGSMap_'

 end subroutine ExGMapGSMap_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ExGSMapGMap_ - Trade of GlobalSegMap structures.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine ExGSMapGMap_(LocalGSMap, LocalComm, RemoteGMap, &
                         Remote_comp_id, ierr) 

!
! !USES:
!
      use m_mpif90
      use m_die
      use m_stdio

      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_init => init
      use m_GlobalMap, only : GlobalMap
      use m_GlobalMap, only : GlobalMap_init => init

      implicit none

      type(GlobalSegMap), intent(in)  :: LocalGSMap  ! Local GlobalSegMap
      integer,            intent(in)  :: LocalComm   ! Local Communicator
      type(GlobalMap),    intent(out) :: RemoteGMap  ! Remote GlobalMap
      integer        , intent(in)  :: Remote_comp_id ! Remote component id
      integer        , intent(out) :: ierr           ! Error Flag

! !REVISION HISTORY:
!       03Feb01 - J.W. Larson <larson@mcs.anl.gov> - API specification.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ExGSMapGMap_'

 end subroutine ExGSMapGMap_

 end module m_ExchangeMaps







