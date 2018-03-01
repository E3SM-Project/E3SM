!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: ThrowFundamentalTypes
!
!> @brief
!! <BriefDescription>
!!
!! @author
!! Tom Clune, NASA/GSFC 
!!
!! @date
!! 07 Nov 2013
!! 
!! @note <A note here.>
!! <Or starting here...>
!
! REVISION HISTORY:
!
! 07 Nov 2013 - Added the prologue for the compliance with Doxygen. 
!
!-------------------------------------------------------------------------------

module ThrowFundamentalTypes_mod

  use Params_mod
  use StringConversionUtilities_mod
  use Exception_mod
  use SourceLocation_mod

  implicit none
  private

  public :: locationFormat
  public :: throwNonConformable
  public :: throwDifferentValues
  public :: throwDifferentValuesWithLocation

!mlr-!!! Can we put these in one place?
!mlr-   integer, parameter :: MAX_LEN_MSG   = 1000
!mlr-   integer, parameter :: MAX_LEN_FLOAT = 25
!mlr-   integer, parameter :: MAX_LEN_INT   = 15
!mlr-  
!mlr-   integer, parameter :: L_INFINITY_NORM = 0
!mlr-   integer, parameter :: L1_NORM         = 1
!mlr-   integer, parameter :: L2_NORM         = 2

!  interface locationFormat
!     module procedure locationFormat
!  end interface locationFormat

  interface throwDifferentValues
     module procedure throwDifferentValues_ii
     module procedure throwDifferentValues_ir
     module procedure throwDifferentValues_rr
  end interface throwDifferentValues

  interface throwDifferentValuesWithLocation
     module procedure throwDifferentValuesWithLocation_ii
     module procedure throwDifferentValuesWithLocation_ir
     module procedure throwDifferentValuesWithLocation_rr
  end interface throwDifferentValuesWithLocation

contains

  ! Consider promoting to module level scope.
  subroutine throwNonConformable(shapeExpected, shapeFound, location)
    integer, intent(in) :: shapeExpected(:)
    integer, intent(in) :: shapeFound(:)
    type (SourceLocation), optional, intent(in) :: location

    call throw( &
         & 'Assertion failed: non-conformable real arrays.' // new_line('$') //&
         & '    expected shape: <['//trim(toString(shapeExpected))//']>' // new_line('$') //&
         & '   but found shape: <['//trim(toString(shapeFound))//']>', &
         & location=location &
         & )
  end subroutine throwNonConformable

  subroutine compareElements(expected, found, i1, i2, location)
    real, intent(in) :: expected, found
    integer, intent(in) :: i1, i2
    type (SourceLocation), optional, intent(in) :: location

    ! the test
    if (expected /= found) then
       call throwDifferentValues(expected, found, i1, i2, 0.0, location=location)
    end if
  end subroutine compareElements

  subroutine throwDifferentValues_ii(iExpected, iFound, i1, i2, tolerance, location)
    integer, intent(in) :: iExpected, iFound
    integer, intent(in) :: i1, i2
    real, intent(in) :: tolerance
    type (SourceLocation), optional, intent(in) :: location

    ! Check with team to see if this is okay.
    call throwDifferentValues_rr(real(iExpected), real(iFound), i1, i2, tolerance, &
         & location=location)

  end subroutine throwDifferentValues_ii

  subroutine throwDifferentValues_ir(iExpected, found, i1, i2, tolerance, location)
    integer, intent(in) :: iExpected
    real, intent(in) :: found
    integer, intent(in) :: i1, i2
    real, intent(in) :: tolerance
    type (SourceLocation), optional, intent(in) :: location

    ! Check with team to see if this is okay.
    call throwDifferentValues_rr(real(iExpected), found, i1, i2, tolerance, &
         & location=location)

  end subroutine throwDifferentValues_ir

  subroutine throwDifferentValues_rr(expected, found, i1, i2, tolerance, &
       & location )
    real, intent(in) :: expected, found
    integer, intent(in) :: i1, i2
    real, intent(in) :: tolerance
    type (SourceLocation), optional, intent(in) :: location

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80

    ! "locationInArray" is not used in the original AssertEqual code.
    character(len=MAXLEN_SHAPE) :: locationInArray
    write(locationInArray,'("[",i0,", ",i0," ]")') i1, i2

    call throw( &
         & 'Assertion failed:  unequal real 2D arrays.' // new_line('$') // &
         & '  First difference at element <' // locationInArray // '>' // &
         & trim(valuesReport(expected, found)) // &
         & trim(differenceReport(found - expected, tolerance)), &
         & location=location &
!         & trim(differenceReport(found - expected, 0.)) &
         & )

  end subroutine throwDifferentValues_rr

  subroutine throwDifferentValuesWithLocation_ii( &
       & iExpected, iFound, iLocation, tolerance, location )
    integer, intent(in) :: iExpected, iFound
    integer, intent(in) :: iLocation(:)
    real, intent(in) :: tolerance
    type (SourceLocation), optional, intent(in) :: location

    ! Check with team to see if this is okay.
    call throwDifferentValuesWithLocation_rr( &
         & real(iExpected), real(iFound), iLocation, tolerance, location=location )

  end subroutine throwDifferentValuesWithLocation_ii

  subroutine throwDifferentValuesWithLocation_ir( &
       & iExpected, found, iLocation, tolerance, location)
    integer, intent(in) :: iExpected
    real, intent(in) :: found
    integer, intent(in) :: iLocation(:)
    real, intent(in) :: tolerance
    type (SourceLocation), optional, intent(in) :: location

    ! Check with team to see if this is okay. ! Answer: meh...
    call throwDifferentValuesWithLocation_rr( &
         & real(iExpected), found, iLocation, tolerance, location=location)

  end subroutine throwDifferentValuesWithLocation_ir

  function locationFormat(iLocation) result (fmt)
    integer, intent(in) :: iLocation(:)

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80*2

    character(len=MAXLEN_SHAPE) :: fmt
    integer :: iLocationSize

    iLocationSize = size(iLocation)

    if (iLocationSize .eq. 0) then
       fmt = '("[", i0, "]")'
    else if (iLocationSize .eq. 1) then
       fmt = '("[", i0, "]")'
    else
       write(fmt,*) '("[",',iLocationSize-1,'(i0,", "), i0, "]")'
    end if

  end function locationFormat

  subroutine throwDifferentValuesWithLocation_rr( &
       & expected, found, iLocation, tolerance, location)
    real, intent(in) :: expected, found
    integer, intent(in) :: iLocation(:)
!    integer :: iLocationSize
    real, intent(in) :: tolerance
    type (SourceLocation), optional, intent(in) :: location

    !mlr maybe move this to a larger scope...
    integer, parameter :: MAXLEN_SHAPE = 80*2

    ! "location" is not used in the original AssertEqual code.
    character(len=MAXLEN_SHAPE) :: locationInArray

    write(locationInArray, locationFormat(iLocation)) iLocation

    call throw( &
         & 'Assertion failed: unequal arrays.' // new_line('$') // &
         & '  First difference at element <' // trim(locationInArray) // '>' // &
         & trim(valuesReport(expected, found)) // &
         & trim(differenceReport(found - expected, tolerance)), &
         & location=location &
         & )

  end subroutine throwDifferentValuesWithLocation_rr

   character(len=MAXLEN_MESSAGE) function valuesReport(expected, found)
      real, intent(in) :: expected
      real, intent(in) :: found

      valuesReport = 'expected: <' // trim(toString(expected)) // '> but found: <' // trim(toString(found)) // '>'
   end function valuesReport

   character(len=MAXLEN_MESSAGE) function differenceReport(difference, tolerance)
      real, intent(in) :: difference
      real, intent(in) :: tolerance
      differenceReport = '    difference: |' // trim(toString(difference)) // '| > tolerance:' // trim(toString(tolerance))
   end function differenceReport

end module ThrowFundamentalTypes_mod
