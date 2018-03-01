module LinearInterpolator_mod
   implicit none
   private

   public :: LinearInterpolator

   type Node
      real :: location
      real :: value
   end type Node

   type LinearInterpolator
      private
      integer :: numNodes
      type(Node), allocatable :: nodes(:)
   contains
      procedure :: getBracket
      procedure :: interpolate
   end type LinearInterpolator

   interface LinearInterpolator
      module procedure newLinearInterpolator
   end interface LinearInterpolator

contains

   function newLinearInterpolator(x, y) result(interpolator)
      type (LinearInterpolator) :: interpolator
      real, intent(in) :: x(:)
      real, intent(in) :: y(:)

      integer :: numNodes

      numNodes = size(x) ! must be same as y 

      interpolator%numNodes = numNodes
      allocate(interpolator%nodes(numNodes))

      interpolator%nodes%location = x
      interpolator%nodes%value = y
      
   end function newLinearInterpolator

   pure function getBracket(this, at) result(upperLower)
      integer :: upperLower(2)
      class (LinearInterpolator), intent(in) :: this
      real, intent(in) :: at

      integer :: i

      upperLower(1) = 1
      upperLower(2) = this%numNodes

      do i = 1, this%numNodes
         if (this%nodes(i)%location > at) then
            upperLower(1) = i - 1
            exit
         else if (this%nodes(i)%location == at) then
            upperLower(1) = i
            exit
         end if
      end do

      do i = this%numNodes, 1, -1
         if (this%nodes(i)%location < at) then
            upperLower(2) = i + 1
            exit
         else if (this%nodes(i)%location == at) then
            upperLower(2) = i
            exit
         end if
      end do

   end function getBracket

   pure function interpolate(this, at) result(value)
      real :: value
      class (LinearInterpolator), intent(in) :: this
      real, intent(in) :: at

      integer :: upperLower(2)
      real :: x1, x2, dx
      real :: y1, y2
      real :: f1, f2

      upperLower = this%getBracket(at)
      
      if (upperLower(1) == upperLower(2)) then ! node point
         value = this%nodes(upperLower(1))%value
      else
         x1 = this%nodes(upperLower(1))%location
         x2 = this%nodes(upperLower(2))%location
         y1 = this%nodes(upperLower(1))%value
         y2 = this%nodes(upperLower(2))%value

         dx = x2 - x1
         f1 = (x2 - at) / dx
         f2 = 1 - f1

         value = f1 * y1 + f2 * y2
      end if
      
   end function interpolate
   
end module LinearInterpolator_mod
