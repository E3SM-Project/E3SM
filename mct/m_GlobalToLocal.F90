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

      public :: GlobalSegMapToLocal ! Translate Global indices from 
                                    ! the GlobalSegMap to local indices.

    interface GlobalSegMapToLocal ; module procedure   &
        GlobalSegMapToIndices_      ! local arrays of starts/lengths
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

 end module m_GlobalToLocal
