module dycore
!
! Data and utility routines related to the dycore
!
   implicit none

PRIVATE

   public :: dycore_is, get_resolution

CONTAINS

   logical function dycore_is (name)
!
! Input arguments
!
      character(len=*) :: name
      
      dycore_is = .false.
      if (name == 'unstructured' .or. name == 'UNSTRUCTURED' .or. &
           name == 'se' .or. name == 'SE') then
         dycore_is = .true.
      end if
      
      return
   end function dycore_is

   character(len=7) function get_resolution()

!     use pmgrid, only: plat

!     select case ( plat )
!     case ( 8 )
!        get_resolution = 'T5'
!     case ( 32 )
!        get_resolution = 'T21'
!     case ( 48 )
!        get_resolution = 'T31'
!     case ( 64 )
!        get_resolution = 'T42'
!     case ( 128 )
!        get_resolution = 'T85'
!     case ( 256 )
!        get_resolution = 'T170'
!     case default
!        get_resolution = 'UNKNOWN'
! This forces the physics settings to be the same as those used for T85 eul
        get_resolution = 'T85'
!     end select

     return
   end function get_resolution

end module dycore


