!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_GlobalToLocal - a set of Global to Local index Mappings.
!
! !DESCRIPTION:
!
! !INTERFACE:

 module m_GlobalToLocal

      implicit none

      private   ! except

      public :: GlobalToLocalIndices ! Translate Global to Local indices
                                     ! (i.e. recover local starts/lengths 
                                     ! of distributed data segments).
                                     
      public :: GlobalToLocalIndex   ! Translate Global to Local index
                                     ! (i.e. recover local index for a
                                     ! point from its global index). 

    interface GlobalToLocalIndices ; module procedure   &
        GlobalSegMapToIndices_,  &   ! local arrays of starts/lengths
        GlobalSegMapToNavigator_     ! return local indices as Navigator
    end interface

    interface GlobalToLocalIndex ; module procedure &
	GlobalSegMapToIndex_, &
	GlobalMapToIndex_
    end interface

! !REVISION HISTORY:
!       02Feb01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=len('m_GlobalToLocal')),parameter :: myname='m_GlobalToLocal'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GlobalSegMapToIndices_ - Return _local_ indices in arrays.
!
! !DESCRIPTION:  {\tt GlobalSegMapToIndices\_()} takes a user-supplied
! {\tt GlobalSegMap} data type {\tt GSMap}, and input communicator 
! {\tt comm} to translate the global directory of segment locations into
! local indices for referencing the on-pe storage of the mapped distributed
! data.
!
! {\bf N.B.:}  This routine returns two allocated arrays---{\tt start(:)} 
! and {\tt length(:)}---which must be deallocated once the user no longer
! needs them.  Failure to do this will create a memory leak.
!
! !INTERFACE:

 subroutine GlobalSegMapToIndices_(GSMap, comm, start, length)

!
! !USES:
!
      use m_mpif90
      use m_die,          only : MP_perr_die
      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_ngseg => ngseg
      use m_GlobalSegMap, only : GlobalSegMap_nlseg => nlseg

      implicit none

      type(GlobalSegMap),intent(in) :: GSMap   ! Output GlobalSegMap
      integer,           intent(in) :: comm    ! communicator handle

      integer,dimension(:), pointer :: start  ! local segment start indices
      integer,dimension(:), pointer :: length ! local segment sizes

! !REVISION HISTORY:
!       02Feb01 - J.W. Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GlobalSegMapToIndices_'

  integer :: myID, ierr, ngseg, nlseg, n, count
 
          ! determine local process id myID

  call MP_COMM_RANK(comm, myID, ierr)
  if(ierr /= 0) call MP_perr_die(myname_,'MP_COMM_RANK',ierr)

          ! determine number of global segments ngseg:

  ngseg = GlobalSegMap_ngseg(GSMap)

          ! determine number of local segments on process myID nlseg:

  nlseg = GlobalSegMap_nlseg(GSMap, myID)

          ! allocate arrays start(:) and length(:) to store local
          ! segment information.

  allocate(start(nlseg), length(nlseg), stat=ierr)
  if(ierr /= 0) call MP_perr_die(myname_,'allocate(start...',ierr)

          ! Loop over GlobalSegMap%pe_loc(:) values to isolate
          ! global index values of local data.  Record number of
          ! matches in the INTEGER count.

  count = 0
  do n=1, ngseg
     if(GSMap%pe_loc(n) == myID) then
        count = count + 1
        if(count > nlseg) then
           ierr = 2
           call MP_perr_die(myname_,'too many pe matches',ierr)
	endif
	start(count) = GSMap%start(count)
	length(count) = GSMap%length(count)
     endif
  end do

  if(count < nlseg) then
     ierr = 3
     call MP_perr_die(myname_,'too few pe matches',ierr)
  endif

          ! translate global start indices to their local 
          ! values, based on their storage order and number
          ! of elements in each segment

  do n=1, count
     if(n == 1) then
	start(n) = 1
     else
	start(n) = start(n-1) + length(n-1)
     endif
  end do

 end subroutine GlobalSegMapToIndices_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GlobalSegMapToIndex_ - translate global to local index.
!
! !DESCRIPTION:  {\tt GlobalSegMapToIndex\_()} takes a user-supplied
! {\tt GlobalSegMap} data type {\tt GSMap}, input global index value 
! {\tt i\_g}, and input communicator {\tt comm} and returns a positive 
! local index value if the datum {\tt i\_g} is on the local communicator.
! If the datum {\tt i\_g} is not local, a value of {\tt -1} is returned.
!
! !INTERFACE:

 integer function GlobalSegMapToIndex_(GSMap, i_g, comm)

!
! !USES:
!
      use m_mpif90
      use m_die,          only : MP_perr_die
      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_ngseg => ngseg
      use m_GlobalSegMap, only : GlobalSegMap_nlseg => nlseg

      implicit none

      type(GlobalSegMap),intent(in)  :: GSMap ! Output GlobalSegMap
      integer,           intent(in)  :: i_g   ! global index
      integer,           intent(in)  :: comm  ! communicator handle

! !REVISION HISTORY:
!       02Feb01 - J.W. Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GlobalSegMapToIndex_'

  integer :: myID
  integer :: count, ierr, ngseg, nlseg, n
  integer :: lower_bound, upper_bound
  integer :: local_start, local_index
  logical :: found

  ! Determine local process id myID:

  call MP_COMM_RANK(comm, myID, ierr)
  if(ierr /= 0) call MP_perr_die(myname_,'MP_COMM_RANK()',ierr)

  ! Extract the global number of segments in GSMap

  ngseg = GlobalSegMap_ngseg(GSMap)

  ! Extract the global number of segments in GSMap for myID

  nlseg = GlobalSegMap_nlseg(GSMap, myID)

  ! set the counter count, which records the number of times myID
  ! matches entries in GSMap%pe_loc(:)

  count = 0

  ! set local_start, which is the current local storage segment
  ! starting position

  local_start = 1

  ! set logical flag found to signify we havent found i_g:

  found = .false.

  n = 0

  SEARCH_LOOP: do 
     
     n = n+1
     if (n > ngseg) EXIT

     if(GSMap%pe_loc(n) == myID) then

  ! increment / check the pe_loc match counter

        count = count + 1
        if(count > nlseg) then
           ierr = 2
           call MP_perr_die(myname_,'too many pe matches',ierr)
	endif

  ! is i_g in this segment?

        lower_bound = GSMap%start(n)
        upper_bound = GSMap%start(n) + GSMap%length(n) - 1

        if((lower_bound <= i_g) .and. (i_g <= upper_bound)) then
	   local_index = local_start + (i_g - GSMap%start(n))
	   found = .true.
	   EXIT
	else
	   local_start = local_start + GSMap%length(n)
        endif

     endif
  end do SEARCH_LOOP

  ! We either found the local index, or have exhausted our options.

  if(found) then
     GlobalSegMapToIndex_ = local_index
  else
     GlobalSegMapToIndex_ = -1
  endif

 end function GlobalSegMapToIndex_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GlobalMapToIndex_ - translate global to local index.
!
! !DESCRIPTION:  {\tt GlobalMapToIndex\_()} takes a user-supplied
! {\tt GlobalMap} data type {\tt GMap}, input global index value 
! {\tt i\_g}, and input communicator {\tt comm} and returns a positive 
! local index value if the datum {\tt i\_g} is on the local communicator.
! If the datum {\tt i\_g} is not local, a value of {\tt -1} is returned.
!
! !INTERFACE:

 integer function GlobalMapToIndex_(GMap, i_g, comm)

!
! !USES:
!
      use m_mpif90
      use m_die,          only : MP_perr_die
      use m_GlobalMap, only : GlobalMap

      implicit none

      type(GlobalMap),intent(in)  :: GMap     ! Input GlobalMap
      integer,           intent(in)  :: i_g   ! global index
      integer,           intent(in)  :: comm  ! communicator handle

! !REVISION HISTORY:
!       02Feb01 - J.W. Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GlobalMapToIndex_'

  integer :: myID
  integer :: count, ierr, ngseg, nlseg, n
  integer :: lower_bound, upper_bound
  integer :: local_start, local_index
  logical :: found

  ! Determine local process id myID:

  call MP_COMM_RANK(comm, myID, ierr)
  if(ierr /= 0) call MP_perr_die(myname_,'MP_COMM_RANK()',ierr)

  ! Initialize logical "point located" flag found as false

  found = .false.

  lower_bound = GMap%displs(myID)
  upper_bound = GMap%displs(myID) + GMap%counts(myID) - 1

  if((lower_bound <= i_g) .and. (i_g <= upper_bound)) then
     found = .true.
     local_index = i_g - lower_bound + 1
  endif

  if(found) then
     GlobalMapToIndex_ = local_index
  else
     GlobalMapToIndex_ = -1
  endif

 end function GlobalMapToIndex_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GlobalSegMapToNavigator_ - Return _local_ Navigator.
!
! !DESCRIPTION:  {\tt GlobalSegMapToNavigator\_()} takes a user-supplied
! {\tt GlobalSegMap} data type {\tt GSMap}, and input communicator 
! {\tt comm} to translate the global directory of segment locations into
! local indices for referencing the on-pe storage of the mapped distributed
! data.  These data are returned in the form of a {\tt Navigator} data 
! type {Nav}.
!
! {\bf N.B.:}  This routine returns a {\tt Navigator} variable {\tt Nav},
! which must be deallocated once the user no longer needs it.  Failure to 
! do this will create a memory leak.
!
! !INTERFACE:

 subroutine GlobalSegMapToNavigator_(GSMap, comm, oNav)

!
! !USES:
!
      use m_mpif90
      use m_die,          only : MP_perr_die
      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_ngseg => ngseg
      use m_GlobalSegMap, only : GlobalSegMap_nlseg => nlseg
      use m_Navigator, only    : Navigator
      use m_Navigator, only    : Navigator_init => init

      implicit none

      type(GlobalSegMap),intent(in) :: GSMap   ! Input GlobalSegMap
      integer,           intent(in) :: comm    ! communicator handle

      type(Navigator),  intent(out) :: oNav    ! Output Navigator

! !REVISION HISTORY:
!       02Feb01 - J.W. Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GlobalSegMapToNavigator_'

  integer :: myID, ierr, ngseg, nlseg, n, count
 
          ! determine local process id myID

  call MP_COMM_RANK(comm, myID, ierr)
  if(ierr /= 0) call MP_perr_die(myname_,'MP_COMM_RANK',ierr)

          ! determine number of global segments ngseg:

  ngseg = GlobalSegMap_ngseg(GSMap)

          ! determine number of local segments on process myID nlseg:

  nlseg = GlobalSegMap_nlseg(GSMap, myID)

          ! Allocate space for the Navigator oNav:

  call Navigator_init(oNav, nlseg, ierr)
  if(ierr /= 0) call MP_perr_die(myname_,'Navigator_init',ierr)

  call GlobalSegMapToIndices_(GSMap, comm, oNav%displs, oNav%counts)

 end subroutine GlobalSegMapToNavigator_

 end module m_GlobalToLocal
