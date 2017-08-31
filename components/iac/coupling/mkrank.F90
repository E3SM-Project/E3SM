!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: mkrank
!
! !INTERFACE:
subroutine mkrank (n, a, miss, iv, num)
!
! !DESCRIPTION:
! Return indices of largest [num] values in array [a]
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
!
! !ARGUMENTS:
  implicit none
  integer , intent(in) :: n        !array length
  real(r8), intent(in) :: a(0:n)   !array to be ranked
  integer , intent(in) :: miss     !missing data value
  integer , intent(in) :: num      !number of largest values requested
  integer , intent(out):: iv(num)  !index to [num] largest values in array [a]
!
! !CALLED FROM:
! subroutine mkpft
! subroutine mksoicol
! subroutine mksoitex
!
! !REVISION HISTORY:
! Author: Gordon Bonan
!
!
! !LOCAL VARIABLES:
!EOP
  real(r8) a_max       !maximum value in array
  integer i            !array index
  real(r8) delmax      !tolerance for finding if larger value
  integer m            !do loop index
  integer k            !do loop index
  logical exclude      !true if data value has already been chosen
!-----------------------------------------------------------------------

  delmax = 1.e-06

  ! Find index of largest non-zero number

  iv(1) = miss
  a_max = -9999.

  do i = 0, n
     if (a(i)>0. .and. (a(i)-a_max)>delmax) then
        a_max = a(i)
        iv(1)  = i
     end if
  end do

  ! iv(1) = miss indicates no values > 0. this is an error

  if (iv(1) == miss) then
     write (6,*) 'MKRANK error: iv(1) = missing'
     call abort()
  end if

  ! Find indices of the next [num]-1 largest non-zero number.
  ! iv(m) = miss if there are no more values > 0

  do m = 2, num
     iv(m) = miss
     a_max = -9999.
     do i = 0, n

        ! exclude if data value has already been chosen

        exclude = .false.
        do k = 1, m-1
           if (i == iv(k)) exclude = .true.
        end do

        ! if not already chosen, see if it is the largest of
        ! the remaining values

        if (.not. exclude) then
           if (a(i)>0. .and. (a(i)-a_max)>delmax) then
              a_max = a(i)
              iv(m)  = i
           end if
        end if
     end do
  end do

  return
end subroutine mkrank
