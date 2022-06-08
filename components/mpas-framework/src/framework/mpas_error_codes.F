! Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_error_codes
!
!> \MPAS handle combined mpas error codes
!> \author Adrian K. Turner
!> \date   8th June 2022
!> \details
!>  This module defines functions to allow mpas error codes to include
!>  two separate codes encoded within the single integer value. The
!>  first half of the bits within the integer error code encode a
!>  primary error code and the second half of the integer bits define
!>  a secondary error code. This functionality allows the error codes
!>  of included libraries (such as PIO) to be passed through to the
!>  calling MPAS cores along with MPAS error codes.
!
!-----------------------------------------------------------------------

module mpas_error_codes

  implicit none
  private

  !--------------------------------------------------------------------
  ! Public member functions
  !--------------------------------------------------------------------
  public :: &
       mpas_combine_error_codes, &
       mpas_primary_error_code, &
       mpas_secondary_error_code, &
       mpas_change_primary_error_code


!***********************************************************************


contains


!***********************************************************************
!
!  routine mpas_combine_error_codes
!
!> \brief   Combine a primary and secondart error code into a single code
!> \author  A. K. Turner
!> \date    8th June 2022
!> \details
!>  A primary and secondart integer error codes are combined into a
!>  single integer error code, with the first half of bits in the
!>  combined error code encoding the primary error code and the second
!>  half of bits in the combined error code encoding the secondary
!>  error code.
!
!-----------------------------------------------------------------------

  function mpas_combine_error_codes(primaryErrorCode, secondaryErrorCode) result(combinedErrorCode)

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------

    integer, intent(in) :: &
         primaryErrorCode, & !< Input: Primary error code
         secondaryErrorCode  !< Input: Secondary error code

    integer :: &
         combinedErrorCode

    integer :: &
         bitSize

    bitSize = bit_size(primaryErrorCode)

    call mvbits(primaryErrorCode,   0, bitSize/2, combinedErrorCode, 0)
    call mvbits(secondaryErrorCode, 0, bitSize/2, combinedErrorCode, bitSize/2)

  end function mpas_combine_error_codes

!***********************************************************************
!
!  routine mpas_primary_error_code
!
!> \brief   Extract the primary error code from a combined error code
!> \author  A. K. Turner
!> \date    8th June 2022
!> \details
!>  Extract the first half of the bits in the combined error code and
!>  extract the primary error code from it.
!
!-----------------------------------------------------------------------

  function mpas_primary_error_code(combinedErrorCode) result(primaryErrorCode)

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------

    integer, intent(in) :: &
         combinedErrorCode !< Input: Combined error code

    integer :: &
         primaryErrorCode

    integer :: &
         bitSize

    bitSize = bit_size(combinedErrorCode)

    primaryErrorCode = 0
    call mvbits(combinedErrorCode, 0, bitSize/2, primaryErrorCode, 0)

    if (btest(combinedErrorCode, bitSize/2-1)) then
       call mvbits(not(0), 0, bitSize/2, primaryErrorCode, bitSize/2)
    endif

  end function mpas_primary_error_code

!***********************************************************************
!
!  routine mpas_secondary_error_code
!
!> \brief   Extract the secondary error code from a combined error code
!> \author  A. K. Turner
!> \date    8th June 2022
!> \details
!>  Extract the first half of the bits in the combined error code and
!>  extract the secondary error code from it.
!
!-----------------------------------------------------------------------

  function mpas_secondary_error_code(combinedErrorCode) result(secondaryErrorCode)

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------

    integer, intent(in) :: &
         combinedErrorCode !< Input: Combined error code

    integer :: &
         secondaryErrorCode

    integer :: &
         bitSize

    bitSize = bit_size(combinedErrorCode)

    secondaryErrorCode = 0
    call mvbits(combinedErrorCode, bitSize/2, bitSize/2, secondaryErrorCode, 0)

    if (btest(combinedErrorCode, bitSize-1)) then
       call mvbits(not(0), 0, bitSize/2, secondaryErrorCode, bitSize/2)
    endif

  end function mpas_secondary_error_code

!***********************************************************************
!
!  routine mpas_change_primary_error_code
!
!> \brief   Change the primary error code of a combined error code
!> \author  A. K. Turner
!> \date    8th June 2022
!> \details
!>  Change the primary error encoded in a combined error code and return
!>  the changed combined error code.
!
!-----------------------------------------------------------------------

  function mpas_change_primary_error_code(combinedErrorCode, newPrimaryErrorCode) result(newCombinedErrorCode)

    !-----------------------------------------------------------------
    ! input variables
    !-----------------------------------------------------------------

    integer, intent(in) :: &
         combinedErrorCode, & !< Input: Combined error code to modify
         newPrimaryErrorCode  !< Input: Primary error code to update with

    integer :: &
         newCombinedErrorCode

    newCombinedErrorCode = mpas_combine_error_codes(newPrimaryErrorCode, mpas_secondary_error_code(combinedErrorCode))

  end function mpas_change_primary_error_code

  !-----------------------------------------------------------------------

end module mpas_error_codes
