!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_ConvertMaps - Conversion between GlobalMap and GlobalSegMap.
!
! !DESCRIPTION:
!
! This module contains routines to convert between the {\tt GlobalMap} 
! and {\tt GlobalSegMap} types.  Since the {\tt GlobalMap} is a 1-D 
! decomposition with one contiguous segment per process, it is always 
! possible to create a {\tt GlobalSegMap} containing the same decomposition
! information.  In the unusual case that a {\tt GlobalSegMap} contains 
! {\em at most} one segment per process, it is possible to create a 
! {\tt GlobalMap} describing the same decomposition.
!
! !INTERFACE:

 module m_ConvertMaps
!
! !USES:
!
      use m_GlobalMap,    only : GlobalMap
      use m_GlobalSegMap, only : GlobalSegMap

      implicit none

      private   ! except

      public :: GlobalMapToGlobalSegMap
      public :: GlobalSegMapToGlobalMap


    interface GlobalMapToGlobalSegMap ; module procedure &
        GlobalMapToGlobalSegMap_
    end interface
    interface GlobalSegMapToGlobalMap ; module procedure &
        GlobalSegMapToGlobalMap_
    end interface

! !REVISION HISTORY:
!       12Feb01 - J.W. Larson <larson@mcs.anl.gov> - initial module
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_ConvertMap'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GlobalMapToGlobalSegMap_() -- GlobalMap -> GlobalSegMap.
!
! !DESCRIPTION:
! This routine takes an input {\tt GlobalMap} argument {\tt GMap}, and
! converts its decomposition information into the output {\tt GlobalSegMap}
! argument {\tt GSMap}.
!
! {\bf N.B.:}  This routine creates an allocated structure {\tt GSMap}.
! The user is responsible for deleting this structure using the {\tt clean()}
! method for the {\tt GlobalSegMap} when {\tt GSMap} is no longer needed.
! Failure to do so will create a memory leak.
!
! !INTERFACE:

 subroutine GlobalMapToGlobalSegMap_(GMap, GSMap)
!
! !USES:
!
      use m_stdio, only : stderr
      use m_die,   only : MP_perr_die

      use m_GlobalMap,    only : GlobalMap

      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_init => init

      use m_MCTWorld, only : ThisMCTWorld
      use m_MCTWorld, only : MCTWorld_ComponentNumProcs => ComponentNumProcs

      implicit none

      type(GlobalMap),    intent(in)  :: GMap
      type(GlobalSegMap), intent(out) :: GSMap

! !REVISION HISTORY:
!       12Feb01 - J.W. Larson <larson@mcs.anl.gov> - Prototype code.
!       24Feb01 - J.W. Larson <larson@mcs.anl.gov> - Finished code.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GlobalMapToGlobalSegMap_'

  integer :: ierr, n, NumProcs
  integer, dimension(:), allocatable :: start, length, pe_loc

       ! Sanity Check -- is GMap the right size?

  NumProcs = MCTWorld_ComponentNumProcs(ThisMCTWorld, GMap%comp_id)
  if(NumProcs /= size(GMap%displs)) then
     write(stderr,'(a,)') myname_,":: Size mismatch-NumProcs = ", &
	  NumProcs,"size(GMap%displs) = ",size(GMap%displs)
     call MP_perr_die(myname_,"component/GlobalMap size mismatch",NumProcs)
  endif

       ! Allocate space for process location

  allocate(start(NumProcs), length(NumProcs), pe_loc(NumProcs), stat=ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,"allocate(start(NumProcs...",ierr)
  endif

       ! Load the arrays:

  do n=1,NumProcs
     start(n) = GMap%displs(n-1)
     length(n) = GMap%counts(n-1)
     pe_loc(n) = n-1
  end do

  call GlobalSegMap_init(GSMap, GMap%comp_id, NumProcs, GMap%gsize, &
                         start, length, pe_loc)

       ! Clean up...

  deallocate(start, length, pe_loc, stat=ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,"deallocate(start,...",ierr)
  endif

 end subroutine GlobalMapToGlobalSegMap_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GlobalSegMapToGlobalMap_() -- GlobalSegMap -> GlobalMap.
!
! !DESCRIPTION:
! This routine takes an input {\tt GlobalSegMap} argument {\tt GSMap}, 
! and examines it to determine whether or not it may be expressed in 
! {\tt GlobalMap} form.  This condition is equivalent to each process in
! the {\tt GSMap} owning {\em at most} one segment.  If this condition is
! satisfied, {\tt GlobalSegMapToGlobalMap\_()} creates an output 
! {\tt GlobalMap} argument {\tt GMap} describing the sam decomposition 
! as {\tt GSMap}.
!
! {\bf N.B.:}  This routine creates an allocated structure {\tt GMap}.
! The user is responsible for deleting this structure using the {\tt clean()}
! method for the {\tt GlobalMap} when {\tt GMap} is no longer needed.
! Failure to do so will create a memory leak.
!
! !INTERFACE:

 subroutine GlobalSegMapToGlobalMap_(GSMap, GMap, NumPes)
!
! !USES:
!
      use m_stdio, only : stderr
      use m_die,   only : MP_perr_die

      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_nlseg => nlseg
      use m_GlobalSegMap, only : GlobalSegMap_active_pes => active_pes

      use m_GlobalMap,    only : GlobalMap
      use m_GlobalMap,    only : GlobalMap_init => init

      implicit none

      type(GlobalSegMap), intent(in)  :: GSMap
      type(GlobalMap),    intent(out) :: GMap
      integer,   intent(in), optional :: NumPes

! !REVISION HISTORY:
!       12Feb01 - J.W. Larson <larson@mcs.anl.gov> - Prototype code.
!       24Feb01 - J.W. Larson <larson@mcs.anl.gov> - Further changes,
!                 but still incomplete--do not call this routine.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GlobalSegMapToGlobalMap_'

  integer :: ierr, n
  integer :: nlseg, NumActive, NumProcs
  integer, dimension(:), allocatable :: NumSegs
  integer, dimension(:), pointer :: pe_list

  if(present(NumPes)) then

     allocate(NumSegs(0:NumPes-1), stat=ierr)
     if(ierr /= 0) then
	call MP_perr_die(myname_,'allocate(NumSegs(1:NumPes-1))=',ierr)
     endif

     do n=0,NumPes-1

       ! Is there at most one segment per process?  If not, die.

	NumSegs(n) = GlobalSegMap_nlseg(GSMap, n)
	if(NumSegs(n) > 1) then
	   call MP_perr_die(myname_,'To many segments, nlseg=',nlseg)
        endif

       ! Are there any segment start indices out of order?
       ! If so, die.

	if(GSMap%start(n+1) <= GSMap%start(n)) then
	   call MP_perr_die(myname_,'segment start indices out of order',n)
	endif

     end do

  else

     call GlobalSegMap_active_pes(GSMap, NumActive, pe_list)

       ! Find the number of processes...

     do n=1,NumActive
	if(n == 1) then
	   NumProcs = pe_list(1)
	else
	   if( pe_list(n) > NumProcs) NumProcs = pe_list(n)
	endif
     end do

     allocate(NumSegs(0:NumActive-1), stat=ierr)
     if(ierr /= 0) then
	call MP_perr_die(myname_,'allocate(NumSegs(1:NumActive-1))=',ierr)
     endif

     do n=0,NumActive-1

       ! Is there at most one segment per process?  If not, die.

	NumSegs(n) = GlobalSegMap_nlseg(GSMap, n)
	if(NumSegs(n) > 1) then
	   call MP_perr_die(myname_,'To many segments, nlseg=',nlseg)
        endif

       ! Are there any segment start indices out of order?
       ! If so, die.

	if(GSMap%start(n+1) <= GSMap%start(n)) then
	   call MP_perr_die(myname_,'segment start indices out of order',n)
	endif

     end do
 
 endif ! if(present(NumPes))

 end subroutine GlobalSegMapToGlobalMap_

 end module m_ConvertMaps
