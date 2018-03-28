module ice_warnings

  use ice_kinds_mod

  implicit none

  private
  save

  ! private warning messages
  character(len=char_len_long), dimension(:), allocatable :: &
       warnings

  integer :: &
       nWarnings

  public :: &
       add_warning, &
       reset_warnings, &
       get_number_warnings, &
       get_warning

!=======================================================================

contains

!=======================================================================

  subroutine add_warning(warning)

    character(len=*), intent(in) :: &
         warning ! warning to add to array of warnings

    ! number of array elements to increase size of warnings array if that array has run out of space
    integer, parameter :: &
         nWarningsBuffer = 10

    ! temporary array to store previous warnings while warning array is increased in size
    character(len=char_len_long), dimension(:), allocatable :: &
         warningsTmp

    integer :: &
         nWarningsArray, & ! size of warnings array at start
         iWarning ! warning index

    ! check if warnings array is not allocated
    if (.not. allocated(warnings)) then

       ! allocate warning array with number of buffer elements
       allocate(warnings(nWarningsBuffer))

       ! set initial number of nWarnings
       nWarnings = 0

    ! already allocated
    else

       ! find the size of the warnings array at the start
       nWarningsArray = size(warnings)
       
       ! check to see if need more space in warnings array
       if (nWarnings + 1 > nWarningsArray) then
       
          ! allocate the temporary warning storage
          allocate(warningsTmp(nWarningsArray))

          ! copy the warnings to temporary storage
          do iWarning = 1, nWarningsArray
             warningsTmp(iWarning) = trim(warnings(iWarning))
          enddo ! iWarning

          ! increase the size of the warning array by the buffer size
          deallocate(warnings)
          allocate(warnings(nWarningsArray + nWarningsBuffer))

          ! copy back the temporary stored warnings
          do iWarning = 1, nWarningsArray
             warnings(iWarning) = trim(warningsTmp(iWarning))
          enddo ! iWarning

          ! deallocate the temporary storage
          deallocate(warningsTmp)

       endif
          
    endif

    ! increase warning number
    nWarnings = nWarnings + 1

    ! add the new warning
    warnings(nWarnings) = trim(warning)

  end subroutine add_warning

!=======================================================================

  subroutine reset_warnings()

    nWarnings = 0

  end subroutine reset_warnings

!=======================================================================

  function get_number_warnings() result(nWarningsOut)

    integer :: nWarningsOut

    nWarningsOut = nWarnings

  end function get_number_warnings

!=======================================================================

  function get_warning(iWarning) result(warning)

    integer, intent(in) :: iWarning

    character(len=char_len_long) :: warning

    if (iWarning <= nWarnings) then
       warning = warnings(iWarning)
    else
       warning = ""
    endif

  end function get_warning

!=======================================================================

end module ice_warnings
