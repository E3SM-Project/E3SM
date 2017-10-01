      module ESMF_ShrTimeMod
!
!==============================================================================
!
! This file contains types and methods that are shared in the hierarchy
!
!------------------------------------------------------------------------------
! INCLUDES

!==============================================================================
!BOPI
! !MODULE: ESMF_ShrTimeMod
!
! !DESCRIPTION:
!
!------------------------------------------------------------------------------
! !USES:
      ! inherit from ESMF base class
      use ESMF_BaseMod

      ! inherit from base time class
      use ESMF_BaseTimeMod
      use ESMF_CalendarMod

      implicit none
!
!------------------------------------------------------------------------------
! !PRIVATE TYPES:
      private
!------------------------------------------------------------------------------
!     ! ESMF_Time
!
!     ! F90 class type to match C++ Time class in size only;
!     !  all dereferencing within class is performed by C++ implementation

     type ESMF_Time
       type(ESMF_BaseTime) :: basetime           ! inherit base class
       ! time instant is expressed as year + basetime
       integer :: YR
       type(ESMF_Calendar), pointer :: calendar => null() ! associated calendar
     end type

     public ESMF_Time
!==============================================================================
end module ESMF_ShrTimeMod
