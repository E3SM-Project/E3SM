!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
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

      public :: ExchangeMap

    interface ExchangeMap ; module procedure   &
        ExGMapGMap_,   &  ! GlobalMap for GlobalMap
        ExGSMapGSMap_, &  ! GlobalSegMap for GlobalSegMap
        ExGMapGSMap_,  &  ! GlobalMap for GlobalSegMap
        ExGSMapGMap_      ! GlobalSegMap for GlobalMap
    end interface

!       03Feb01 - J.W. Larson <larson@mcs.anl.gov> - initial module
!EOP ___________________________________________________________________


  character(len=*),parameter :: myname='m_ExchangeMaps'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ExGMapGMap_ - Trade of GlobalMap structures
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine ExGMapGMap_(LocalGMap, RemoteGMap, Remote_comp_id, ierr) 

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
      type(GlobalMap), intent(out) :: RemoteGMap     ! Remote GlobalMap
      integer        , intent(in)  :: Remote_comp_id ! Remote component id
      integer        , intent(out) :: ierr           ! Error Flag

! !REVISION HISTORY:
!       03Feb01 - J.W. Larson <larson@mcs.anl.gov> - API specification.
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

 subroutine ExGSMapGSMap_(LocalGSMap, RemoteGSMap, Remote_comp_id, ierr) 

!
! !USES:
!
      use m_mpif90
      use m_die
      use m_stdio
      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_init => init

      implicit none

      type(GlobalSegMap), intent(in)  :: LocalGSMap   ! Local GlobalSegMap
      type(GlobalSegMap), intent(out) :: RemoteGSMap  ! Remote GlobalSegMap
      integer        , intent(in)  :: Remote_comp_id  ! Remote component id
      integer        , intent(out) :: ierr            ! Error Flag

! !REVISION HISTORY:
!       03Feb01 - J.W. Larson <larson@mcs.anl.gov> - API specification.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ExGSMapGSMap_'

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

 subroutine ExGMapGSMap_(LocalGMap, RemoteGSMap, Remote_comp_id, ierr) 

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
      type(GlobalSegMap), intent(out) :: RemoteGSMap  ! Remote GlobalSegMap
      integer        , intent(in)  :: Remote_comp_id  ! Remote component id
      integer        , intent(out) :: ierr            ! Error Flag

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

 subroutine ExGSMapGMap_(LocalGSMap, RemoteGMap, Remote_comp_id, ierr) 

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
      type(GlobalMap),    intent(out) :: RemoteGMap  ! Remote GlobalMap
      integer        , intent(in)  :: Remote_comp_id ! Remote component id
      integer        , intent(out) :: ierr           ! Error Flag

! !REVISION HISTORY:
!       03Feb01 - J.W. Larson <larson@mcs.anl.gov> - API specification.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ExGSMapGMap_'

 end subroutine ExGSMapGMap_

 end module m_ExchangeMaps
