module add_mod
   implicit none

contains

   real function add(x, y) result (sum)
      implicit none
      real, intent(in) :: x
      real, intent(in) :: y

      sum = x + y

   end function add
end module add_mod
