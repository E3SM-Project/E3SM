module MpiTestParameter_mod
   use AbstractTestParameter_mod
   implicit none
   private

   public :: MpiTestParameter
   public :: newMpiTestParameter

   ! We can request as many processes as we like.
   ! If available, actual = request.  Otherwise throw an exception at run()
   type, extends(AbstractTestParameter) :: MpiTestParameter
      integer :: numProcessesRequested
   contains
      procedure :: setNumProcessesRequested
      procedure :: getNumProcessesRequested
      procedure :: toString
      procedure :: toStringActual
   end type MpiTestParameter

!!$   interface MpiTestParameter
!!$      module procedure :: newMpiTestParameter
!!$   end interface MpiTestParameter

contains

   ! Note that npes requested may not be available. 
   function newMpiTestParameter(numProcessesRequested) result(testParameter)
      type (MpiTestParameter) :: testParameter
      integer, intent(in) :: numProcessesRequested
      
      call testParameter%setNumProcessesRequested(numProcessesRequested)
      
   end function newMpiTestParameter

   pure subroutine setNumProcessesRequested(this, numProcessesRequested)
      class (MpiTestParameter), intent(inout) :: this
      integer, intent(in) :: numProcessesRequested
      this%numProcessesRequested = numProcessesRequested
   end subroutine setNumProcessesRequested


   ! This function ensures that "npes = #" is included in the message string 
   ! for each exception.   It should rarely be overridden.
   function toStringActual(this) result(string)
      class (MpiTestParameter), intent(in) :: this
      character(:), allocatable :: string

      character(len=8) :: npesString
      character(:), allocatable :: tmp

      write(npesString,'(i0)') this%numProcessesRequested

      string = 'npes=' // trim(npesString) 
      tmp = this%toString()

      if (len_trim(tmp) > 0) then
         string = string // ' :: ' // trim(tmp)
      end if

   end function toStringActual

   
   ! Provide a default empty string.  It is expected that this function
   ! will be overridden for user defined test cases.
   function toString(this) result(string)
      class (MpiTestParameter), intent(in) :: this
      character(:), allocatable :: string

      string = ''

   end function toString

   pure integer function getNumProcessesRequested(this) result(numProcessesRequested)
      class (MpiTestParameter), intent(in) :: this
      numProcessesRequested = this%numProcessesRequested
   end function getNumProcessesRequested


end module MpiTestParameter_mod
