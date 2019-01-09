module dycore

implicit none
private

public :: dycore_is

!=========================================================================================
CONTAINS
!=========================================================================================

logical function dycore_is (name)

   character(len=*) :: name

   dycore_is = .false.
   if (name == 'unstructured' .or. name == 'UNSTRUCTURED' .or. &
      name == 'se' .or. name == 'SE') then
      dycore_is = .true.
   end if

end function dycore_is

!=========================================================================================

end module dycore
